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

function absorptive_tem_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix
    use m_multislice
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 4 + 4*n_slices
    
    required_memory = array_count * array_size
    
end function
    


subroutine absorptive_tem
    
    use global_variables
    use m_lens
    use m_user_input
    use m_absorption
    use m_precision
    use output
	use cufft_wrapper
#ifdef GPU
    use CUFFT
    use cuda_array_library
    use cudafor
    use cuda_ms
    use cuda_potential
    use cuda_setup
#endif
    use m_tilt
    use m_multislice
    use m_potential
    use m_string
    implicit none
    
    !dummy variables
    integer(4) :: i_cell,i_slice,z_indx(1),lengthz,lengthdf,i,idum,ntilt
    
    !random variables
    integer(4) :: count
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix)  :: psi,psi_initial,psi_out
    complex(fp_kind) :: lens_ctf(nopiy,nopix,imaging_ndf)
    complex(fp_kind),dimension(nopiy,nopix,n_slices)::projected_potential,prop,transf_absorptive
    
    !output
    real(fp_kind),dimension(nopiy,nopix) :: cbed,image,tem_image,temp_image
    
    !diagnostic variables
    real(fp_kind) :: intensity, t1, delta
    
    !output variables
    character(120) :: filename,fnam_df
    
#ifdef GPU
    !device variables
	integer :: plan
	complex(fp_kind),device,allocatable :: prop_d(:,:,:),transf_d(:,:,:)
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d, psi_out_d
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,trans_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_dwf_d,fz_abs_d
        
    real(fp_kind) :: absorptive_tem_GPU_memory

    call GPU_memory_message(absorptive_tem_GPU_memory(), on_the_fly)
#else
	on_the_fly=.false.
#endif    
    
    call command_line_title_box('Pre-calculation setup')
	do i=1,imaging_ndf
        if(pw_illum) lens_ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,imaging_aberrations,imaging_apodisation)
	enddo
    ! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
	
    if (pw_illum) then
        psi_initial = 1.0_fp_kind/sqrt(float(nopiy*nopix))
	else
        psi_initial = make_ctf(probe_initial_position,probe_df(1),probe_cutoff,probe_aberrations,probe_apodisation)
        call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
        psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
    endif
    call tilt_wave_function(psi_initial)
    if(.not.load_grates) then
    	projected_potential = make_absorptive_grates(nopiy,nopix,n_slices)
    endif
    call load_save_add_grates(projected_potential,nopiy,nopix,n_slices)
    

	t1 = secnds(0.0)
#ifdef GPU
	! Plan the fourier transforms
    if (fp_kind.eq.8)then
          call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
          call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
    !Copy host arrays to the device
    if (on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix),inverse_sinc_d(nopiy,nopix))
        allocate(trans_d(nopiy,nopix),fz_d(nopiy,nopix,nt))
        allocate(fz_dwf_d(nopiy,nopix,nt),fz_abs_d(nopiy,nopix,nt))
        fz_d=fz 
        fz_dwf_d=fz_dwf
        fz_abs_d = ci*absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak  !make the potential absorptive
        inverse_sinc_d=inverse_sinc
        bwl_mat_d = bwl_mat
    else
	    allocate(transf_d(nopiy,nopix,n_slices))
    endif
    allocate(prop_d(nopiy,nopix,n_slices))
#endif
    call command_line_title_box('Calculation running')
    
    do ntilt=1,n_tilts_total
	!For each specimen tilt we have to redo the transmission function
	!and propagators
#ifdef GPU								
    if (on_the_fly) call cuda_setup_many_phasegrate
#endif
    
    do i = 1, n_slices
	    call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory)
        
	    prop(:,:,i) = prop(:,:,i) * bwl_mat
        if(.not.on_the_fly) then
		transf_absorptive(:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,i))
	endif
        ! Bandwith limit the phase grate, psi is used for temporary storage
        call fft2(nopiy, nopix, transf_absorptive(:,:,i), nopiy, psi, nopiy)
        psi = psi * bwl_mat
        call ifft2(nopiy, nopix, psi, nopiy, transf_absorptive(:,:,i), nopiy)
    enddo
   lengthz = ceiling(log10(maxval(abs(zarray))))
   lengthdf = ceiling(log10(maxval(abs(imaging_df))))
   if(any(imaging_df<0)) lengthdf = lengthdf+1
#ifdef GPU
    if(.not.on_the_fly) transf_d=transf_absorptive
    psi_d = psi_initial
    prop_d = prop
    do i_cell = 1, maxval(ncells)
        do i_slice = 1, n_slices
			
            if(on_the_fly) then
                call cuda_make_abs_potential(trans_d,ccd_slice_array(i_slice),tau_slice(:,:,:,i_slice),nat_slice(:,i_slice),prop_distance(i_slice),plan,fz_d,fz_dwf_d,fz_abs_d,inverse_sinc_d,bwl_mat_d,Volume_array(i_slice))
                call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
            else
                call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,i_slice), psi_out_d,1.0_fp_kind,nopiy,nopix)
            endif
            
            call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
            call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,i_slice), psi_out_d,normalisation,nopiy,nopix)
            call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)
        enddo ! End loop over slices
        
		!If this thickness corresponds to any of the output values then output images
		if (any(i_cell==ncells)) then
			psi = psi_d
			call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
			cbed = abs(psi)**2
			z_indx = minloc(abs(ncells-i_cell))
			filename = trim(adjustl(output_prefix))
			if (nz>1) filename = trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(z_indx(1))),lengthz)//'_A'
            if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
			call binary_out_unwrap(nopiy, nopix, cbed, trim(adjustl(filename))//'_DiffractionPattern')
				
			if(pw_illum) then
			do i=1,imaging_ndf
				psi = psi_d
				call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
				psi = psi*lens_ctf(:,:,i)
				call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
				tem_image = abs(psi)**2
				fnam_df = trim(adjustl(filename))// '_Image'
				if(imaging_ndf>1) fnam_df = trim(adjustl(fnam_df))//'_Defocus_'//zero_padded_int(int(imaging_df(i)),lengthdf)//'_Ang'
				call binary_out(nopiy, nopix, tem_image, fnam_df)
			enddo
			endif
		endif
        intensity = get_sum(psi_d)
	    write(6,900,advance='no') achar(13), i_cell, intensity
900     format(a1, 1x, 'Cell: ', i5, ' Intensity: ', f12.6)
#else
    psi = psi_initial
    do i_cell = 1, maxval(ncells)
        do i_slice = 1, n_slices	
            call multislice_iteration(psi,prop(:,:,i_slice),transf_absorptive(:,:,i_slice),nopiy,nopix);												  
        enddo ! End loop over slices
        
		!If this thickness corresponds to any of the output values then output images
		if (any(i_cell==ncells)) then
			  
			call fft2(nopiy, nopix, psi, nopiy, psi_out, nopiy)
			cbed = abs(psi_out)**2
			z_indx = minloc(abs(ncells-i_cell))
			filename = trim(adjustl(output_prefix))
			if (nz>1) filename = trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(z_indx(1))),lengthz)//'_A'
            if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
			call binary_out_unwrap(nopiy, nopix, cbed, trim(adjustl(filename))//'_DiffractionPattern')
            
            
			if(pw_illum) then
            call binary_out(nopiy, nopix, abs(psi)**2, trim(adjustl(filename))//'_exit_surface_intensity')
            call binary_out(nopiy, nopix, atan2(imag(psi),real(psi)), trim(adjustl(filename))//'_exit_surface_phase')
			do i=1,imaging_ndf
			   
				call fft2(nopiy, nopix, psi, nopiy, psi_out, nopiy)
																		  
				psi_out = psi_out*lens_ctf(:,:,i)
				call ifft2(nopiy, nopix, psi_out, nopiy, psi_out, nopiy)
				tem_image = abs(psi_out)**2
				fnam_df = trim(adjustl(filename))// '_Image'
				if(imaging_ndf>1) fnam_df = trim(adjustl(fnam_df))//'_Defocus_'//zero_padded_int(int(imaging_df(i)),lengthdf)//'_Ang'
				call binary_out(nopiy, nopix, tem_image, fnam_df)
            enddo
            endif
		endif
        intensity = sum(abs(psi)**2)
        
        if(n_tilts_total<2) write(6,900) i_cell, intensity
        if(n_tilts_total>1) write(6,901) i_cell, ntilt,n_tilts_total, intensity
900     format(1h+,1x, 'Cell: ', i5, ' Intensity: ', f12.6)	        
901     format(1h+,1x, 'Cell: ', i5, ' tilt:',i5,'/',i5,' Intensity: ', f12.6)	
#endif	
        
	enddo ! End loop over cells
    enddo ! End loop over tilts
    
    
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
    

end subroutine absorptive_tem
