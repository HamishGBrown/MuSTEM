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
	use cufft_wrapper
    use m_slicing
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
    complex(fp_kind),dimension(nopiy,nopix) :: psi,psi_initial,trans,psi_out
    complex(fp_kind),dimension(nopiy,nopix,nz) :: psi_elastic
         
    real(fp_kind),dimension(nopiy, nopix) :: fourDSTEM_pattern,fourDSTEM_el_pattern,temp
	real(fp_kind),dimension(nopiy,nopix,nz) :: pacbed_pattern,cbed
    real(fp_kind)::ccd_slice(n_slices)
    
    !diagnostic variables
    real(fp_kind) :: intensity,t1, delta
    
    !output variables
    character(120) :: filename

	logical:: fourDSTEM,elfourd
    
    ccd_slice = relm / (tp * ak * ss_slice(7,:))
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)
    
	write(*,*) 'Output diffraction patterns for each scan position?'
    write(*,*) '<1> Output the diffraction pattern for each scan position'
    write(*,*) '<2> Do not output diffraction patterns for each scan position'
	call get_input('Diffraction pattern for each probe position? <1> y <2> n',idum)
	write(*,*)
	
	fourDSTEM = idum == 1
	
    call precalculate_scattering_factors
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    if (on_the_fly) then
        ! call cuda_setup_many_phasegrate
        
    else
        idum = seed_rng()
        call make_qep_grates(idum)
    
    endif
    
 	call setup_propagators
       
	t1 = secnds(0.0)
    
    
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)
        	
    intensity = 1.0_fp_kind
    pacbed_pattern = 0.0_fp_kind
    
	do ny = 1, nysample
    do nx = 1, nxsample
#ifdef gpu       
        write(6,903,ADVANCE='NO') achar(13), ny, nysample, nx, nxsample, intensity
903     format(a1, ' y:', i3, '/', i3, ' x:', i3, '/', i3, '  Intensity:', f6.3, ' (to monitor BWL)')	
#else
        write(6,900) ny, nysample, nx, nxsample, intensity
900     format(1h+,  ' y:', i3, '/', i3, ' x:', i3, '/', i3, '  Intensity:', f6.3, ' (to monitor BWL)')
#endif
    
        call make_stem_wfn(psi_initial, probe_df(1), probe_positions(:,ny,nx),probe_aberrations)
        
        call tilt_wave_function(psi_initial)
        
        psi_initial = psi_initial
        
		!Reset 4DSTEM pattern to zero
		fourDSTEM_pattern = 0_fp_kind
		fourDSTEM_el_pattern = 0_fp_kind
		psi_elastic = 0
		cbed = 0

        do i_qep_pass = 1, n_qep_passes 
        
            ! Reset wavefunction
            psi = psi_initial
            
			

            do i = 1,maxval(ncells)
	            do j = 1, n_slices
                    
                    ! Phase grate
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    if(on_the_fly) then
						!call make_qep_potential(trans, tau_slice, nat_slice, ss_slice(7,j))
						!psi_out = psi*trans
                    elseif(quick_shift) then
                        shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                        shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
						trans = cshift(cshift(qep_grates(:,:,nran,j),shifty,dim=1),shiftx,dim=2)
						psi_out = psi*trans
                    elseif(phase_ramp_shift) then                       !randomly shift phase grate
                        shiftx = floor(ifactorx*ran1(idum)) + 1
                        shifty = floor(ifactory*ran1(idum)) + 1
						call phase_shift_array(qep_grates(:,:,nran,j),trans,shift_arrayy,shift_arrayx)
                        psi_out = psi*trans
                    else
						psi_out = psi*qep_grates(:,:,nran,j)
                    endif
						
					call fft2(nopiy,nopix,psi_out,nopiy,psi,nopiy)
					psi_out = prop(:,:,j)*psi
					call ifft2(nopiy,nopix,psi_out,nopiy,psi,nopiy)
                                        
		        enddo ! End loop over slices
				
				!If this thickness is an output thickness then accumulate pacbed and output relevant diffraction patterns
				if (any(i==ncells)) then
					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
					z_indx = minloc(abs(ncells-i))

					! Accumulate elastic wave function
					if(elfourd)  psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi_out

					! Accumulate diffaction pattern
					temp = abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp
				endif
				
            enddo ! End loop over cells
            
	    enddo ! End loop over QEP passes
        
        intensity = sum(abs(psi)**2)
		length = ceiling(log10(maxval(zarray)))
		do i=1,nz
			fourDSTEM_pattern = cbed(:,:,i)
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
					fourDSTEM_el_pattern = abs(psi_elastic(:,:,i))**2
					call binary_out_unwrap(nopiy, nopix, fourDSTEM_el_pattern/n_qep_passes**2, filename,write_to_screen=.false.)
			endif
		
		
		enddo

	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	
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
