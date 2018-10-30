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

function qep_tem_GPU_memory(n_qep_grates, quick_shift, phase_ramp_shift) result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly
    use m_lens, only: imaging,imaging_ndf
    use m_multislice
    
    implicit none
    
	integer*4,intent(in)::n_qep_grates
	logical,intent(in)::quick_shift,phase_ramp_shift
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(3 + n_slices + n_slices*n_qep_grates) + 3
    
    array_count = array_count + 2 + 1 +imaging_ndf
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function
    


subroutine qep_tem
    
    use m_precision;use m_numerical_tools;use global_variables;use m_lens;use cufft_wrapper
#ifdef GPU	
    use cudafor, only: dim3;use cuda_ms;use cuda_potential
    use cuda_setup, only: GPU_memory_message;use cuda_array_library
#endif
    use output;use m_tilt;use m_multislice;use m_potential;use m_string;use m_hn0
    
    implicit none
    
    !dummy variables
    integer(4) :: i_cell, i_slice, i_qep_pass,i,j,ntilt
    integer(4) :: shifty, shiftx,z_indx(1),starting_slice
    
    !random variables
    integer(4) :: idum
    !real(fp_kind) :: ran1
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix,n_slices)::prop
    complex(fp_kind),dimension(nopiy,nopix,n_qep_grates,n_slices)::projected_potential,qep_grates
    complex(fp_kind),dimension(nopiy,nopix) :: psi,psi_initial,psi_out,psi_temp,trans
	complex(fp_kind),dimension(nopiy,nopix,imaging_ndf)::ctf
    complex(fp_kind)::psi_elastic(nopiy,nopix,nz)
	
	
    !output
    real(fp_kind),dimension(nopiy,nopix,nz) :: cbed,total_intensity
	real(fp_kind)::tem_image(nopiy,nopix,nz,imaging_ndf)
    real(fp_kind),dimension(nopiy,nopix) :: image,temp
	character*120 ::fnam_df
	integer :: lengthdf,lengthz
	logical::manyz,manytilt,factorized_propagator
    
    !diagnostic variables
    real(fp_kind) :: intensity, t1, delta
#ifdef GPU    
    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d, psi_initial_d,psi_out_d,psi_temp_d
    complex(fp_kind),device :: psi_elastic_d(nopiy,nopix,nz)
	complex(fp_kind),device,allocatable :: prop_d(:,:,:), transf_d(:,:,:,:), ctf_d(:,:,:), shift_arrayx_d(:,:), shift_arrayy_d(:,:), shift_array_d(:,:), trans_d(:,:),lens_ctf(:,:,:)
    real(fp_kind),device,allocatable :: tem_image_d(:,:,:,:)
    real(fp_kind),device,dimension(nopiy,nopix) :: temp_d,temp_image_d
    real(fp_kind),device,dimension(nopiy,nopix,nz) :: cbed_d, total_intensity_d
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d, inverse_sinc_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d
	
	!Double channeling variables
	complex(fp_kind),device,allocatable,dimension(:,:) ::psi_inel_d,shiftarray,tmatrix_d,q_tmatrix_d
	complex(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states_d,Hn0_shifty_coord_d,Hn0_shiftx_coord_d
	real(fp_kind),device,allocatable,dimension(:,:)::cbed_inel_dc_d
	real(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states
	real(fp_kind),device,allocatable,dimension(:,:,:,:)::eftem_image_d
	integer::jj,l,i_target,k,ii,i_df
#endif
    real(fp_kind) :: qep_tem_GPU_memory
    
    character*100::filename
#ifdef GPU    
    call GPU_memory_message(qep_tem_GPU_memory(n_qep_grates, quick_shift, phase_ramp_shift), on_the_fly)
#endif    

    manyz = nz>1
	manytilt = n_tilts_total>1
    lengthz = calculate_padded_string_length(zarray,nz)
    lengthdf = calculate_padded_string_length(imaging_df,imaging_ndf)
	
    call command_line_title_box('Pre-calculation setup')
    call precalculate_scattering_factors
	idum = seed_rng()
#ifdef GPU						
    if (on_the_fly) then
        call cuda_setup_many_phasegrate
    else
	
	 if(double_channeling) then
        allocate(tmatrix_states_d(nopiy,nopix,nstates),psi_inel_d(nopiy,nopix),cbed_inel_dc_d(nopiy,nopix),tmatrix_states(nopiy,nopix,nstates))
		allocate(shiftarray(nopiy,nopix),tmatrix_d(nopiy,nopix),q_tmatrix_d(nopiy,nopix))
        tmatrix_states_d = setup_ms_hn0_tmatrices(nopiy,nopix,nstates)*alpha_n
        allocate(Hn0_shifty_coord_d(nopiy,maxval(natoms_slice_total),n_slices))
        allocate(Hn0_shiftx_coord_d(nopix,maxval(natoms_slice_total),n_slices))
        Hn0_shiftx_coord_d = Hn0_shiftx_coord
        Hn0_shifty_coord_d = Hn0_shifty_coord
		allocate(eftem_image_d(nopiy,nopix,imaging_ndf,nz));eftem_image_d=0
    endif 

#endif        
    if(.not. load_grates) projected_potential = make_qep_grates(idum)
    call load_save_add_grates(idum,projected_potential,nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice)
        
#ifdef GPU        
    endif
#else
	on_the_fly=.false.
#endif

	if(pw_illum) then
		do i=1,imaging_ndf
			ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,imaging_aberrations,imaging_apodisation)
		enddo
#ifdef GPU		
		allocate(ctf_d(nopiy,nopix,imaging_ndf))
		ctf_d = ctf
	    allocate(tem_image_d(nopiy,nopix,nz,imaging_ndf))
		tem_image_d = 0.0_fp_kind
#endif
	endif

    call command_line_title_box('Calculation running')
#ifdef GPU
    if (fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
	allocate(prop_d(nopiy,nopix,n_slices))
	prop_d = prop

    if (on_the_fly) then
        if(allocated(bwl_mat_d)) deallocate(bwl_mat_d)
        if(allocated(inverse_sinc_d)) deallocate(inverse_sinc_d)
        if(allocated(fz_d)) deallocate(fz_d)
        allocate(bwl_mat_d(nopiy,nopix))
        allocate(inverse_sinc_d(nopiy,nopix))
        allocate(fz_d(nopiy,nopix,nt))
        fz_d = fz
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
        transf_d = qep_grates
        
        if (qep_mode == 3) then
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d=shift_arrayx
            shift_arrayy_d=shift_arrayy
        endif
    endif
    
    if (any([on_the_fly,quick_shift,phase_ramp_shift])) allocate(trans_d(nopiy,nopix))
#endif    
	t1 = secnds(0.0)

	if (pw_illum) then
	      psi_initial = 1.0_fp_kind / sqrt(float(npixels))
	else
        psi_initial = make_ctf(probe_initial_position,probe_df(1),probe_cutoff,probe_aberrations,probe_apodisation)
        call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
        psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
	endif
    
    call tilt_wave_function(psi_initial)
    do ntilt = 1,n_tilts_total
        
	do i = 1, n_slices
	    call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory)
	    prop(:,:,i) = prop(:,:,i) * bwl_mat
        if(.not.on_the_fly) then
		qep_grates(:,:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,:,i))
		do j=1,n_qep_grates
		call fft2(nopiy,nopix,qep_grates(:,:,j,i),nopiy,psi,nopiy)
		if(qep_mode == 3) qep_grates(:,:,j,i)= psi*bwl_mat
		if(qep_mode .ne. 3) call ifft2(nopiy,nopix,psi*bwl_mat,nopiy,qep_grates(:,:,j,i),nopiy)
		enddo
		endif
    enddo
    
#ifdef GPU
prop_d=prop
if(.not.on_the_fly) transf_d = qep_grates
psi_initial_d = psi_initial
        
	! Set accumulators to zero
    cbed_d = 0.0_fp_kind
    psi_elastic_d = 0.0_fp_kind
    total_intensity_d = 0.0_fp_kind

    
    do i_qep_pass = 1, n_qep_passes 
    
        ! Reset wavefunction
        psi_d = psi_initial_d
        
        do i_cell = 1,maxval(ncells)
	        do i_slice = 1, n_slices
				
				if(double_channeling) then
                do i_target = 1, natoms_slice_total(i_slice) ! Loop over targets
                ! Calculate inelastic transmission matrix
					call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,i_slice),Hn0_shiftx_coord_d(:,i_target,i_slice),nopiy,nopix)
					
					do k = 1, nstates
							
					call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
                    call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
					call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,sqrt(normalisation),nopiy,nopix)
					starting_slice = i_slice
					do ii = i_cell, n_cells
                        do jj = starting_slice, n_slices;nran = floor(n_qep_grates*ran1(idum)) + 1
							shiftx = floor(ifactorx*ran1(idum));shifty = floor(ifactory*ran1(idum))
							if(on_the_fly) then
								call cuda_fph_make_potential(trans_d,ccd_slice_array(jj),tau_slice,nat_slice(:,jj),jj,prop_distance(jj),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
								call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							elseif(qep_mode == 2) then
								call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,jj),trans_d,nopiy,nopix,shifty* nopiy_ucell,shiftx* nopix_ucell)
								call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							elseif(qep_mode == 3) then                       !randomly shift phase grate
								call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,jj),trans_d,shift_arrayy_d(:,shifty+1),shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
								call cuda_multislice_iteration(psi_inel_d, trans_d, prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							else
								call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,nran,jj), prop_d(:,:,jj), normalisation, nopiy, nopix,plan)
							endif;enddo
							
						if (any(ii==ncells)) then
							z_indx = minloc(abs(ncells-ii))
							call cufftExec(plan, psi_inel_d, tmatrix_d, CUFFT_FORWARD)
								
							! Accumulate EFTEM images
							do l = 1, imaging_ndf
								call cuda_image(tmatrix_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
									 
								call cuda_addition<<<blocks,threads>>>(eftem_image_d(:,:,l,z_indx(1)), temp_d, eftem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
							enddo
						endif
						starting_slice =1
                    enddo
					enddo
				enddo
			endif  ! End loop over cells,targets and states and end double_channeling section

				
                ! QEP multislice
				nran = floor(n_qep_grates*ran1(idum)) + 1
				shiftx = floor(ifactorx*ran1(idum));shifty = floor(ifactory*ran1(idum))
				if(on_the_fly) then
					call cuda_fph_make_potential(trans_d,ccd_slice_array(i_slice),tau_slice,nat_slice(:,i_slice),i_slice,prop_distance(i_slice),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
					call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,i_slice), normalisation, nopiy, nopix,plan)
				elseif(qep_mode == 2) then
					call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,i_slice),trans_d,nopiy,nopix,shifty* nopiy_ucell,shiftx* nopix_ucell)
					call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,i_slice), normalisation, nopiy, nopix,plan)
				elseif(qep_mode == 3) then                       !randomly shift phase grate
					call cuda_phase_shift_from_1d_factor_arrays(transf_d(:,:,nran,i_slice),trans_d,shift_arrayy_d(:,shifty+1),shift_arrayx_d(:,shiftx+1),nopiy,nopix,plan)
					call cuda_multislice_iteration(psi_d, trans_d, prop_d(:,:,i_slice), normalisation, nopiy, nopix,plan)
				else
					call cuda_multislice_iteration(psi_d, transf_d(:,:,nran,i_slice), prop_d(:,:,i_slice), normalisation, nopiy, nopix,plan)
				endif


            enddo ! End loop over slices
			
				!If this thickness is an output thickness then accumulate relevent TEM images and diffraction patterns
				if (any(i_cell==ncells)) then
					
					!Transform into diffraction space
					call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
					call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
					z_indx = minloc(abs(ncells-i_cell))

					! Accumulate elastic wave function
					call cuda_addition<<<blocks,threads>>>(psi_elastic_d(:,:,z_indx(1)),psi_d,psi_elastic_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)

					! Accumulate exit surface intensity
					call cuda_mod<<<blocks,threads>>>(psi_d, temp_image_d, 1.0_fp_kind, nopiy, nopix)
					call cuda_addition<<<blocks,threads>>>(total_intensity_d(:,:,z_indx(1)), temp_image_d, total_intensity_d(:,:,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
					
					! Accumulate diffaction pattern
					call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
					call cuda_addition<<<blocks,threads>>>(cbed_d(:,:,z_indx(1)),temp_d,cbed_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
					
					if(pw_illum) then
					! Accumulate image
					do i=1,imaging_ndf
						psi_temp_d = psi_out_d
						call cuda_multiplication<<<blocks,threads>>>(psi_temp_d, ctf_d(:,:,i), psi_temp_d, sqrt(normalisation), nopiy, nopix)
						call cufftExec(plan, psi_temp_d, psi_temp_d, CUFFT_INVERSE)
						call cuda_mod<<<blocks,threads>>>(psi_temp_d, temp_image_d, normalisation, nopiy, nopix)
						call cuda_addition<<<blocks,threads>>>(tem_image_d(:,:,z_indx(1),i), temp_image_d, tem_image_d(:,:,z_indx(1),i), 1.0_fp_kind, nopiy, nopix)
					enddo
					endif

				endif
        enddo ! End loop over cells
        
        intensity = get_sum(psi_d)
		write(6,900,advance='no') achar(13), i_qep_pass, n_qep_passes, intensity
    900 format(a1, 1x, 'QEP pass:', i4, '/', i4, ' Intensity: ', f8.3)	
#else
    cbed = 0.0_fp_kind
    psi_elastic = 0.0_fp_kind
    total_intensity = 0.0_fp_kind
    
    do i_qep_pass = 1, n_qep_passes 
    
        ! Reset wavefunction
		psi = psi_initial
        
        do i_cell = 1,maxval(ncells)
	        do i_slice = 1, n_slices
            
                ! Phase grate
				nran = floor(n_qep_grates*ran1(idum)) + 1
				
                  call qep_multislice_iteration(psi,prop(:,:,i_slice),qep_grates(:,:,:,i_slice),nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,qep_mode,shift_arrayy,shift_arrayx)


            enddo ! End loop over slices
			
				!If this thickness is an output thickness then accumulate relevent TEM images and diffraction patterns
				if (any(i_cell==ncells)) then
					
					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
					temp = abs(psi_out)**2
					z_indx = minloc(abs(ncells-i_cell))

					! Accumulate elastic wave function
					psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi

					! Accumulate exit surface intensity
					total_intensity(:,:,z_indx(1)) = total_intensity(:,:,z_indx(1)) + abs(psi)**2
					
					! Accumulate diffaction pattern
					temp=abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp
					
					! Accumulate image
                    if(pw_illum) then
					do i=1,imaging_ndf
						psi_temp = ctf(:,:,i)*psi_out
						call ifft2(nopiy,nopix,psi_temp,nopiy,psi_temp,nopiy)
						tem_image(:,:,z_indx(1),i) = tem_image(:,:,z_indx(1),i)+abs(psi_temp)**2
                    enddo
                    endif

				endif
        enddo ! End loop over cells
        
        intensity = sum(abs(psi)**2)
		
        write(6,900) i_qep_pass, n_qep_passes, intensity
900     format(1h+,  'QEP pass:', i4, '/', i4, ' Intensity: ', f8.3)
#endif
	enddo ! End loop over QEP passes
#ifdef GPU    
    cbed = cbed_d
    psi_elastic = psi_elastic_d
    total_intensity = total_intensity_d
	if(pw_illum) tem_image = tem_image_d
#endif    
    delta = secnds(t1)
    
    write(*,*)
    write(*,*)
    write(*,*) 'Calculation is finished.'
    write(*,*)
    write(*,*) 'Time elapsed: ', delta, ' seconds.'
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
    
    cbed = cbed / n_qep_passes
    psi_elastic = psi_elastic / n_qep_passes
    
	if(pw_illum) then
		total_intensity = total_intensity / n_qep_passes
		tem_image = tem_image / n_qep_passes
	endif
	
	do i=1,nz
		filename = trim(adjustl(output_prefix))
		if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),lengthz)//'_A'
		
		if(.not.output_thermal) then
			call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffractionPattern')
			if(pw_illum) call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity')
		else
			call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffractionPattern_Total')
			call fft2(nopiy, nopix, psi_elastic(:,:,i), nopiy, psi, nopiy)
			image = abs(psi)**2
			call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffractionPattern_Elastic')
			
			image = cbed(:,:,i) - image
			call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffractionPattern_TDS')
			
			if(pw_illum) then
				call binary_out(nopiy, nopix, abs(psi_elastic(:,:,i))**2, trim(filename)//'_ExitSurface_Intensity_Elastic')
				call binary_out(nopiy, nopix, atan2(imag(psi_elastic(:,:,i)), real(psi_elastic(:,:,i))), trim(filename)//'_ExitSurface_Phase_Elastic')
				call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity_Total')
			
				total_intensity(:,:,i) = total_intensity(:,:,i) - abs(psi_elastic(:,:,i))**2
				call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity_TDS')
			endif
		endif
		if(pw_illum) then
		do j=1,imaging_ndf
			! Elastic image
			call fft2 (nopiy, nopix, psi_elastic(:,:,i), nopiy, psi, nopiy)
			psi = psi * ctf(:,:,j)
			call ifft2 (nopiy, nopix, psi, nopiy, psi, nopiy)
			image = abs(psi)**2
				
			if(imaging_ndf>1) then
                fnam_df = trim(adjustl(filename))//'_Defocus_'//zero_padded_int(int(imaging_df(j)),lengthdf)//'_Ang'
            else
                fnam_df = trim(adjustl(filename))
            endif
			if(output_thermal) then
			call binary_out(nopiy, nopix, image, trim(fnam_df)//'_Image_Elastic')
				  
			! Inelastic image
			image = tem_image(:,:,i,j) - image
			call binary_out(nopiy, nopix, image, trim(fnam_df)//'_Image_TDS')
				  
			! Total image
			call binary_out(nopiy, nopix, tem_image(:,:,i,j), trim(fnam_df)//'_Image_Total')
			else
            call binary_out(nopiy, nopix, tem_image(:,:,i,j), trim(fnam_df)//'_Image')
            endif
#ifdef GPU            
			if(double_channeling) then
				image = eftem_image_d(:,:,j,i)
				call output_TEM_result(output_prefix,tile_out_image(image,ifactory,ifactorx),'energy_filtered_image',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i)&
								&,lengthz=lengthz,lengthdf=lengthdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(j))
            endif
#endif
		enddo
		endif
	
    enddo
    
    enddo !end loop over tilts
end subroutine qep_tem