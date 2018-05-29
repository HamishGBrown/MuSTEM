!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D’Alfonso, S.D. Findlay, B. D. Forbes
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!   
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!                       
!--------------------------------------------------------------------------------

function qep_pacbed_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly
    use m_lens, only: imaging
    use m_slicing, only: n_slices
    use m_qep, only: n_qep_grates, quick_shift, phase_ramp_shift
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(3 + n_slices + n_slices*n_qep_grates) + 2
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_pacbed
    
    use global_variables
    use m_lens
    use m_user_input
    use m_qep
    use m_precision
    use output
    use cuda_array_library
    use cudafor
    use cuda_ms
    use CUFFT
	use cufft_wrapper
    use cuda_potential
    use m_slicing
    use cuda_setup
	use m_string
    use m_probe_scan, only: nysample, nxsample, probe_positions
    use m_tilt, only: tilt_wave_function
    use m_multislice, only: make_qep_grates, setup_propagators
    use m_potential
    
    implicit none
    
    !dummy variables
    integer(4) ::  i, j, i_qep_pass
    integer(4) ::  shifty,shiftx
    integer(4) ::  ny, nx,z_indx(1),length
    
    !random variables
    integer(4) :: idum
    !real(fp_kind) :: ran1
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,psi_initial,psi_elastic
         
    real(fp_kind),dimension(nopiy, nopix) :: fourDSTEM_pattern,cbed,fourDSTEM_pattern_el
	real(fp_kind),dimension(nopiy,nopix,nz)::pacbed_pattern
    
    !diagnostic variables
    real(fp_kind) :: intensity,t1, delta
    
    !output variables
    character(120) :: filename

    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_initial_d, psi_d, psi_out_d,fourDSTEM_el_pattern_d
    real(fp_kind),device,dimension(nopiy,nopix) :: temp_cbed_d, pacbed_pattern_d,fourDSTEM_pattern_d,temp_d
	real(fp_kind),device,dimension(nopiy,nopix),allocatable :: cbed_d(:,:,:)
	complex(fp_kind),device,allocatable :: prop_d(:,:,:), transf_d(:,:,:,:),psi_elastic_d(:,:,:)
	complex(fp_kind),device,allocatable :: trans_d(:,:), shift_arrayx_d(:,:), shift_arrayy_d(:,:), shift_array_d(:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d, inverse_sinc_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d

	logical:: fourDSTEM,elfourd

    real(fp_kind) :: qep_pacbed_GPU_memory
    
    
    
    call GPU_memory_message(qep_pacbed_GPU_memory(), on_the_fly)
    
    
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)
    
	write(*,*) 'Output diffraction patterns for each scan position?'
    write(*,*) '<1> Output the total diffraction pattern (including elastically and '
	write(*,*) '    inelastically scattered contributions) for each scan position' 
    write(*,*) '<2> Output both the total diffraction pattern and also the elastically'
	write(*,*) '    scattered contribution for each scan position'
    write(*,*) '<3> Do not output diffraction patterns for each scan position'
	call get_input('Diffraction pattern for each probe position choice',idum)
	write(*,*)
	
	fourDSTEM = (idum == 1.or.idum ==2)
	elfourd = idum == 2
	
	if(elfourd) allocate(psi_elastic_d(nopiy,nopix,nz))
	allocate(cbed_d(nopiy,nopix,nz))
	
    call precalculate_scattering_factors
    idum = seed_rng()
    if (on_the_fly) then
        call cuda_setup_many_phasegrate
        
    else
        
        call make_qep_grates(idum)
    
    endif
    
 	call setup_propagators
    
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) then
        allocate(trans_d(nopiy,nopix))
    endif
    
	t1 = secnds(0.0)
    
    
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)
    
	! Plan the fourier transforms
    if(fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
    ! Copy host arrays to the device
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    
    if(on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix))
        allocate(inverse_sinc_d(nopiy,nopix))
        allocate(fz_d(nopiy,nopix,nt))
        fz_d = fz 
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
  	    transf_d = qep_grates
        
        if(phase_ramp_shift) then 
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d = shift_arrayx
            shift_arrayy_d = shift_arrayy
        endif
    endif
	
    intensity = 1.0_fp_kind
    pacbed_pattern_d = 0.0_fp_kind
    
	do ny = 1, nysample
    do nx = 1, nxsample
        write(6,903,ADVANCE='NO') achar(13), ny, nysample, nx, nxsample, intensity
903     format(a1, ' y:', i3, '/', i3, ' x:', i3, '/', i3, '  Intensity:', f6.3, ' (to monitor BWL)')	
    
        call make_stem_wfn(psi_initial, probe_df(1), probe_positions(:,ny,nx),probe_aberrations)
        
        call tilt_wave_function(psi_initial)
        
        psi_initial_d = psi_initial
        
		!Reset 4DSTEM pattern to zero
		fourDSTEM_pattern_d = 0_fp_kind
		fourDSTEM_el_pattern_d = 0_fp_kind
		psi_elastic_d = 0
		cbed_d = 0

        do i_qep_pass = 1, n_qep_passes 
        
            ! Reset wavefunction
            psi_d = psi_initial_d
            
			

            do i = 1,maxval(ncells)
	            do j = 1, n_slices
                    
                    ! Phase grate
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    if(on_the_fly) then
                        call cuda_fph_make_potential(trans_d,ccd_slice_array(j),tau_slice,nat_slice(:,j),j,prop_distance(j),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(quick_shift) then
                        shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                        shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
                        call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,j),trans_d,nopiy,nopix,shifty,shiftx)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(phase_ramp_shift) then                       !randomly shift phase grate
                        shiftx = floor(ifactorx*ran1(idum)) + 1
                        shifty = floor(ifactory*ran1(idum)) + 1
                        call cuda_make_shift_array<<<blocks,threads>>>(shift_array_d,shift_arrayy_d(:,shifty),shift_arrayx_d(:,shiftx),nopiy,nopix)     !make the qspace shift array
                        call cuda_multiplication<<<blocks,threads>>>(transf_d(:,:,nran,j),shift_array_d, trans_d,1.0_fp_kind,nopiy,nopix) !multiply by the qspace shift array
                        call cufftExec(plan,trans_d,trans_d,CUFFT_INVERSE)                                                                    !inverse fourier transform
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,sqrt(normalisation),nopiy,nopix)              !do the phase grate multiplication
                    else
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,nran,j), psi_out_d,1.0_fp_kind,nopiy,nopix)
                    endif
                    
                    ! Propagate
			    	call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,j), psi_out_d,normalisation,nopiy,nopix)
                    call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)
                    
		        enddo ! End loop over slices
				
				!If this thickness is an output thickness then accumulate pacbed and output relevant diffraction patterns
				if (any(i==ncells)) then
					!Transform into diffraction space
					call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
					z_indx = minloc(abs(ncells-i))

					! Accumulate elastic wave function
					if(elfourd) call cuda_addition<<<blocks,threads>>>(psi_elastic_d(:,:,z_indx(1)),psi_out_d,psi_elastic_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
            

					! Accumulate diffaction pattern
					call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
					call cuda_addition<<<blocks,threads>>>(cbed_d(:,:,z_indx(1)),temp_d,cbed_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)

				endif
				
            enddo ! End loop over cells
          
		
            
	    enddo ! End loop over QEP passes
        
        intensity = get_sum(psi_d)
		length = ceiling(log10(maxval(zarray)))
		do i=1,nz
			fourDSTEM_pattern = cbed_d(:,:,i)
			!Output 4D STEM diffraction pattern
			if(fourDSTEM) then
					!Output total (elastic and inelastic) diffraction pattern
					filename = trim(adjustl(output_prefix))
					if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(i)))//'_A'
					call binary_out_unwrap(nopiy, nopix, fourDSTEM_pattern/n_qep_passes, trim(adjustl(filename)) //'_pp_'//&
										  &to_string(nx)//'_'//to_string(ny)//'_Diffraction_pattern',write_to_screen=.false.)
			endif

			pacbed_pattern(:,:,i) = pacbed_pattern(:,:,i) + fourDSTEM_pattern

			if (output_thermal.and.fourDSTEM) then 
					!Output elastic only diffraction pattern
					filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_Elastic_Diffraction_pattern'
					call cuda_mod<<<blocks,threads>>>(psi_elastic_d(:,:,i), fourDSTEM_pattern_d, normalisation, nopiy, nopix)
					fourDSTEM_pattern_el = fourDSTEM_pattern_d
					call binary_out_unwrap(nopiy, nopix, fourDSTEM_pattern_el/n_qep_passes**2, filename,write_to_screen=.false.)
			endif
		
		
		enddo

	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	
    ! Copy PACBED pattern to CPU
    !pacbed_pattern = pacbed_pattern_d
    
    ! QEP normalisation
    pacbed_pattern = pacbed_pattern / n_qep_passes
    
    delta = secnds(t1)
    
    write(*,*) 
    write(*,*) 
    
    write(*,*) 'Calculation is finished.'
    write(*,*) 
    write(*,*) 'Time elapsed ', delta, ' seconds.'
    write(*,*)    
    
	if(timing) then
		open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
		write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
		close(9834)
	endif
    
    if (fp_kind.eq.8) then
        write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	else
        write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	endif
    
    write(*,*)

    length = ceiling(log10(maxval(zarray)))
	do i=1,nz
		filename = trim(adjustl(output_prefix))
		if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A_'
		filename = trim(adjustl(filename))//'_PACBED_Pattern'
		call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i),filename)
	enddo
    
end subroutine qep_pacbed
