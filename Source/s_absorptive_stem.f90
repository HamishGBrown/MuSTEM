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

function absorptive_stem_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, adf, eels,zarray,nz,ncells
    use m_multislice
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 6 + 4*n_slices
    if (adf) array_count = array_count + 1
    if (adf) array_count = array_count + 1 + n_slices
    !if (ionization) array_count = array_count + 1 + n_slices
    if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
end function
    


subroutine absorptive_stem(STEM,ionization,PACBED)
    
    use m_precision, only: fp_kind
    use global_variables!, only: nopiy, nopix, normalisation, nt, n_cells, fz, fz_dwf, inverse_sinc, bwl_mat, prop, ndet, ionization, eels, adf, on_the_fly, inner, outer, ci,zarray,nz,ncells
    use m_lens!, only: probe_df, probe_ndf, make_stem_wfn
    use m_absorption!, only: fz_abs, transf_absorptive, calculate_absorption_mu, calculate_local_adf_mu
#ifdef GPU
    use cuda_array_library, only: cuda_multiplication, cuda_addition, cuda_mod, blocks, threads
    use cudafor, only: dim3
    use cuda_ms, only: cuda_stem_detector, get_sum
    use cuda_potential, only: volume_array, ccd_slice_array, cuda_setup_many_phasegrate, cuda_make_adf_potential, cuda_make_ion_potential, cuda_make_abs_potential
    use cufft, only: cufft_z2z, cufft_c2c, cufft_forward, cufft_inverse, cufftplan, cufftexec
    use cuda_setup, only: GPU_memory_message
#endif
    use m_crystallography    
    use m_multislice
    use m_tilt
    use m_potential
    use m_string
	use output
	use cufft_wrapper
    
    implicit none
    
    
    logical,intent(in)::pacbed,stem,ionization
    
    !dummy variables
    integer(4) :: i, j, k, i_df, ny, nx,z_indx(1),ii,length,idum,ntilt,lengthdf,idet

    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,qpsi,psi_out
    complex(fp_kind),dimension(nopiy,nopix,n_slices)::projected_potential,prop,transf_absorptive

    !output/detectors
    real(fp_kind),dimension(nopiy,nopix) :: image,temp_image,psi_intensity,temp
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz,ndet) :: stem_image,stem_elastic_image&
															 &,stem_inelastic_image

    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz)::eels_correction_image
	real(fp_kind),allocatable :: probe_intensity(:,:,:),stem_ion_image(:,:,:,:,:),pacbed_pattern(:,:,:)&
                               &,adf_image(:,:,:),ion_image(:,:,:),masks(:,:,:)
	real(8)::thmin,thmax

    !diagnostic variables
    real(fp_kind) :: intens, t1, delta
    
    !output variables
    character(120) :: filename,fnam,fnam_det
    logical:: many_df
    
#ifdef GPU
    !device variables
    integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d,psi_out_d
    real(fp_kind),device,dimension(nopiy,nopix) :: temp_d
    real(fp_kind),device, allocatable,dimension(:,:) :: psi_intensity_d,eels_correction_detector_d,inelastic_potential_d
    real(fp_kind),device,allocatable,dimension(:,:,:) :: adf_image_d,ion_image_d,pacbed_pattern_d,masks_d
	real(fp_kind),device,allocatable,dimension(:,:,:,:) :: adf_potential_d,ion_potential_d
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,trans_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_dwf_d,fz_abs_d,fz_mu_d, prop_d,transf_d
	complex(fp_kind),device,allocatable,dimension(:,:,:,:) :: fz_adf_d
#endif    
    integer,parameter :: iu = 8945
    
    real(fp_kind) :: absorptive_stem_GPU_memory
    
#ifdef GPU    
    call GPU_memory_message(absorptive_stem_GPU_memory(), on_the_fly)
#endif
	
    call command_line_title_box('Pre-calculation setup')
    
    ! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
#ifdef GPU
#else
	on_the_fly = .false.
#endif
    if (on_the_fly) then
#ifdef GPU
        ! Setup the atom co-ordinate for on the fly potentials
        call cuda_setup_many_phasegrate()                   
#endif        
    else
        if(.not.load_grates) then
		projected_potential = make_absorptive_grates(nopiy,nopix,n_slices)
	endif
        call load_save_add_grates(projected_potential,nopiy,nopix,n_slices)
        if(ionization.or.stem) call make_local_inelastic_potentials(ionization)
    endif
    
    
    if(pacbed) then
        allocate(pacbed_pattern(nopiy,nopix,nz)); pacbed_pattern= 0
        
    endif
    many_df = probe_ndf .gt. 1
    length = ceiling(log10(maxval(zarray)))
    
		
    lengthdf = ceiling(log10(maxval(abs(probe_df))))
    if(any(probe_df<0)) lengthdf = lengthdf+1
    
   
    ! Make detector mask arrays (adds the elastic contribution)        
	if(stem) then
        allocate(masks(nopiy,nopix,ndet))
		do i=1,ndet/nseg
		do j=1,nseg
			if(nseg>1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),2*pi*j/nseg-seg_det_offset,2*pi/nseg)
			if(nseg==1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))
		enddo
		enddo
		if (adf) then
		    allocate(adf_image(nopiy,nopix,ndet))
		    stem_inelastic_image = 0.0_fp_kind
        endif
    endif
		if (ionization) then
			allocate(ion_image(nopiy,nopix,num_ionizations))
			allocate(stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
			stem_ion_image = 0.0_fp_kind
			if (.not.EDX) eels_correction_image = 0.0_fp_kind
		endif
#ifdef GPU
    !plan the fourier transforms
    if (fp_kind.eq.8) then
        call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
    else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
    ! Copy host arrays to the device
    if (on_the_fly) then
        ! Set up device variables for on-the-fly potentials
        
        allocate(bwl_mat_d(nopiy,nopix),inverse_sinc_d(nopiy,nopix)&
				&,trans_d(nopiy,nopix),fz_d(nopiy,nopix,nt)&
				&,fz_dwf_d(nopiy,nopix,nt),fz_abs_d(nopiy,nopix,nt))
		
        fz_d = fz 
        fz_dwf_d = fz_dwf
		if(stem.and.adf)  then
			allocate(fz_adf_d(nopiy,nopix,nt,ndet))
			do k=1,ndet/nseg
			  thmin =  atan(inner((k-1)*nseg+1)/ak)
			  thmax =  atan(outer((k-1)*nseg+1)/ak)
		  
			  !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
			  !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
			  fz_adf_d(:,:,:,(k-1)*nseg+1:k*nseg) = spread(cmplx(absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog,thmin,thmax)/nseg,0.0_fp_kind)*4*pi,dim=4,ncopies=nseg)/nseg
			enddo
		endif
		
		fz_abs_d=0
        if(include_absorption) fz_abs_d = ci*absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak  !make the potential absorptive
		
        if (ionization) then
           allocate(fz_mu_d(nopiy,nopix,num_ionizations))
            fz_mu_d = ionization_mu
        endif
        
        if (ionization.or.adf) allocate(inelastic_potential_d(nopiy,nopix))
        
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
        
    else
        ! Set up device variables for precalculated potentials
        if (stem.and.adf) then
            allocate(adf_potential_d(nopiy,nopix,n_slices,ndet))
            adf_potential_d = adf_potential
        endif
        
        if(ionization) then
            allocate(ion_potential_d(nopiy,nopix,num_ionizations,n_slices))
            ion_potential_d = ionization_potential
        endif
        
    endif
    if (adf.or.ionization) allocate(psi_intensity_d(nopiy,nopix))
    if (stem.and.adf) allocate(adf_image_d(nopiy,nopix,ndet))
    if (ionization) allocate(ion_image_d(nopiy,nopix,num_ionizations))
	if (pacbed) allocate(pacbed_pattern_d(nopiy,nopix,nz)) 
    allocate(prop_d(nopiy,nopix,n_slices),transf_d(nopiy,nopix,n_slices))
    prop_d = prop
    if (stem) allocate(masks_d(nopiy,nopix,ndet))
    if (stem) masks_d = masks
    
    if(ionization.and.(.not.EDX)) then
        allocate(eels_correction_detector_d(nopiy,nopix))
        eels_correction_detector_d = eels_correction_detector
    endif
#endif

    call command_line_title_box('Calculation running')
    t1 = secnds(0.0)
    
    intens = 1.0_fp_kind

    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
    
    do ntilt=1,n_tilts_total
	!Have to redo transmission functions and propagator for each tilt
	
    
	do i = 1, n_slices
        transf_absorptive(:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,i))
	    call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory)
	    prop(:,:,i) = prop(:,:,i) * bwl_mat
        !! Bandwith limit the phase grate, psi is used for temporary storage
        call fft2(nopiy, nopix, transf_absorptive(:,:,i), nopiy, psi, nopiy)
        psi = psi * bwl_mat
        call ifft2(nopiy, nopix, psi, nopiy, transf_absorptive(:,:,i), nopiy)
    enddo
    
	if(pacbed) pacbed_pattern = 0
#ifdef GPU
    if(pacbed) pacbed_pattern_d = 0
	prop_d = prop
	transf_d = transf_absorptive
#endif
    do ny = 1, nysample
    do nx = 1, nxsample
    do i_df = 1, probe_ndf                                                                  !loop over probe points
#ifdef GPU
        
    if(n_tilts_total>1) then
    write(6, 902, advance='no') achar(13), i_df,probe_ndf,ny, nysample, nx, nxsample,ntilt,n_tilts_total, intens
    else
    write(6, 901, advance='no') achar(13), i_df,probe_ndf,ny, nysample, nx, nxsample, intens
    endif
901     format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3)	
902     format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,' tilt:',i3,'/',i3'  Intensity:', f6.3)	
#else
        if(n_tilts_total>1) then
        write(6, 902) i_df,probe_ndf,ny, nysample, nx, nxsample,ntilt,n_tilts_total, intens
        else
        write(6, 901) i_df,probe_ndf,ny, nysample, nx, nxsample, intens
        endif
        
901     format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3)	
902     format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,' tilt:',i3,'/',i3'  Intensity:', f6.3)	
#endif  

        if (output_probe_intensity) probe_intensity = 0_fp_kind
		
		!Make STEM probe
		psi = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
        psi = psi/sqrt(sum(abs(psi)**2))
		
        call tilt_wave_function(psi)

#ifdef GPU
        if (adf) adf_image_d = 0.0_fp_kind
        if (ionization) ion_image_d = 0.0_fp_kind
        psi_d = psi
        do i = 1,maxval(ncells);do j = 1, n_slices
            
			! Calculate inelastic cross sections
			if ((stem.and.adf).or.ionization) call cuda_mod<<<blocks,threads>>>(psi_d,psi_intensity_d,1.0_fp_kind,nopiy,nopix)
			if((stem.and.adf)) then
				do k=1,ndet
				if(on_the_fly) then
					call cuda_make_adf_potential(inelastic_potential_d,tau_slice(:,:,:,j),nat_slice(:,j),plan,fz_adf_d(:,:,:,k),inverse_sinc_d,Volume_array(j))
					call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
				else
					call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,adf_potential_d(:,:,j,k), temp_d,prop_distance(j),nopiy,nopix)     !overlap
				endif
				call cuda_addition<<<blocks,threads>>>(temp_d,adf_image_d(:,:,k),adf_image_d(:,:,k),1.0_fp_kind,nopiy,nopix)                          !depth sum
				enddo
			endif
			if(ionization) then
				do ii=1,num_ionizations
					if(on_the_fly) then
						call cuda_make_ion_potential(inelastic_potential_d,tau_slice(:,atm_indices(ii),:,j),nat_slice(atm_indices(ii),j),plan,&
														&fz_mu_d(:,:,ii),inverse_sinc_d,Volume_array(j))
						call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
					else
						call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,ion_potential_d(:,:,ii,j), temp_d,prop_distance(j),nopiy,nopix)     !overlap
					endif
					call cuda_addition<<<blocks,threads>>>(temp_d,ion_image_d(:,:,ii),ion_image_d(:,:,ii),1.0_fp_kind,nopiy,nopix)                          !depth sum
				enddo
			endif
			! Transmit through slice potential
			if(on_the_fly) then
				call cuda_make_abs_potential(trans_d,ccd_slice_array(j),tau_slice(:,:,:,j),nat_slice(:,j),prop_distance(j),plan,fz_d,fz_dwf_d,fz_abs_d,inverse_sinc_d,bwl_mat_d,Volume_array(j))
				call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_d,1.0_fp_kind,nopiy,nopix)
			else
				call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,j), psi_d,1.0_fp_kind,nopiy,nopix)
			endif
			! Propagate to next slice
			call cufftExec(plan,psi_d,psi_d,CUFFT_FORWARD)
			call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,j), psi_d,normalisation,nopiy,nopix)
			call cufftExec(plan,psi_d,psi_d,CUFFT_INVERSE)
			
			if (output_probe_intensity) then
				k = (i-1)*n_slices+j
				if (output_cell_list(k)) then
					psi = psi_d
					probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
				endif
			endif
		enddo;
		!If this thickness corresponds to any of the output values then output images
			if (any(i==ncells)) then
				z_indx = minloc(abs(ncells-i))
				call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
				call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
				
				if(pacbed.and.(i_df==1)) call cuda_addition<<<blocks,threads>>>(temp_d,pacbed_pattern_d(:,:,z_indx(1)),pacbed_pattern_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
				do ii=1,ndet
					if(stem) stem_elastic_image(ny,nx,i_df,z_indx(1),ii) = cuda_stem_detector(temp_d,masks_d(:,:,ii))
					if(stem.and.adf) stem_inelastic_image(ny,nx,i_df,z_indx(1),ii) = get_sum(adf_image_d(:,:,ii))
				enddo
				if(ionization) then
					do ii=1,num_ionizations
						 stem_ion_image(ny,nx,i_df,z_indx(1),ii) = get_sum(ion_image_d(:,:,ii))
					enddo
					if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = cuda_stem_detector(temp_d,eels_correction_detector_d)
				endif
				!Output 4D STEM diffraction pattern
				if(fourDSTEM) then
						
						temp = temp_d
						filename = trim(adjustl(output_prefix))
						if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(z_indx(1))))//'_A'
						if(many_df) filename = trim(adjustl(filename)) // '_Defocus_'//zero_padded_int(int(probe_df(i_df)),lengthdf)//'_Ang'
                        if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
						filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_abs_Diffraction_pattern'
						call binary_out_unwrap(nopiy, nopix, temp, filename,write_to_screen=.false.,nopiyout=nopiyout,nopixout=nopixout)
				endif
			endif
		
		enddo ! End loop over cells
        intens = get_sum(psi_d)
#else
        if (adf.and.stem) adf_image = 0.0_fp_kind
        if (ionization) ion_image = 0.0_fp_kind
        do i = 1,maxval(ncells);do j = 1, n_slices
			! Calculate inelastic cross sections
			if ((stem.and.adf).or.ionization) psi_intensity = abs(psi)**2
			
			if(stem.and.adf) then
                do k=1,ndet;adf_image(:,:,k) = psi_intensity*adf_potential(:,:,j,k)*prop_distance(j) + adf_image(:,:,k); enddo
            endif    
            
			if(ionization) ion_image = ion_image + spread(psi_intensity,ncopies=num_ionizations,dim=3) * ionization_potential(:,:,:,j) * prop_distance(j)
			call multislice_iteration(psi,prop(:,:,j),transf_absorptive(:,:,j),nopiy,nopix)
			
			if (output_probe_intensity) then
				k = (i-1)*n_slices+j
				if (output_cell_list(k)) then
					probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
				endif
			endif
        enddo
		
			!If this thickness corresponds to any of the output values then output images
			if (any(i==ncells)) then
				z_indx = minloc(abs(ncells-i))
				call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
				temp = abs(psi_out)**2
				if(stem) then
                    do k=1,ndet
                        stem_elastic_image(ny,nx,i_df,z_indx(1),k) = sum(masks(:,:,k)*temp)
                    enddo
                endif
				if(pacbed.and.(i_df==1)) pacbed_pattern(:,:,z_indx(1)) = pacbed_pattern(:,:,z_indx(1)) + temp
				if((stem.and.adf)) then
                    do k=1,ndet; stem_inelastic_image(ny,nx,i_df,z_indx(1),k) = sum(adf_image(:,:,k)); enddo
                endif
				if(ionization) then
					do ii=1,num_ionizations
						stem_ion_image(ny,nx,i_df,z_indx(1),ii) = sum(ion_image(:,:,ii))
                    enddo
                    !call binary_out(nopiy,nopix,eels_correction_detector,'eels_correction_detector')
                    !call binary_out(nopiy,nopix,temp,'temp')
                    !stop
					if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(temp*eels_correction_detector)
				endif
				 
				!Output 4D STEM diffraction pattern
				if(fourDSTEM) then
						filename = trim(adjustl(output_prefix))
						if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(z_indx(1))))//'_A'
						if(many_df) filename = trim(adjustl(filename)) // '_Defocus_'//zero_padded_int(int(probe_df(i_df)),lengthdf)//'_Ang'
                        if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
						filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_abs_Diffraction_pattern'
						call binary_out_unwrap(nopiy, nopix, temp, filename,write_to_screen=.false.,nopiyout=nopiyout,nopixout=nopixout)
				endif
				
			endif
		enddo ! End loops over cells and slices
        intens = sum(abs(psi)**2)
#endif		
		if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,1,probe_ndf,nysample,nxsample)        
    enddo                                                                         !end loop over defocus
    enddo                                                                         !end loop over probe position ny
    enddo                                                                         !end loop over probe position nx
    
    delta = secnds(t1)
    
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
    fnam = trim(adjustl(output_prefix))
    if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
    do idet = 1, ndet
    if(output_thermal) then    
        
        filename = trim(adjustl(fnam)) // '_DiffPlaneElastic_Detector'
        call add_zero_padded_int(filename, fnam_det, idet, 3)
		if (n_tilts_total>1) fnam_det = trim(adjustl(fnam_det))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
        if(stem) call output_stem_image(stem_elastic_image(:,:,:,:,idet),fnam_det,probe_df)
        
        if (stem.and.adf) then
            stem_image(:,:,:,:,idet) = stem_elastic_image(:,:,:,:,idet) + stem_inelastic_image(:,:,:,:,idet)
    
            filename = trim(adjustl(fnam)) // '_DiffPlaneTotal_Detector'
            call add_zero_padded_int(filename, fnam_det, idet, 3)
            if (n_tilts_total>1) fnam_det = trim(adjustl(fnam_det))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
            call output_stem_image(stem_image(:,:,:,:,idet),fnam_det,probe_df)
    
            filename = trim(adjustl(fnam)) // '_DiffPlaneTDS_Detector'
            call add_zero_padded_int(filename, fnam_det, idet, 3)
            if (n_tilts_total>1) fnam_det = trim(adjustl(fnam_det))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
            call output_stem_image(stem_inelastic_image(:,:,:,:,idet),fnam_det,probe_df)
        endif   
    elseif(stem) then
        if (adf) stem_image = stem_elastic_image + stem_inelastic_image
        if (.not.adf) stem_image = stem_elastic_image
        filename = trim(adjustl(fnam)) // '_DiffPlane'
        call add_zero_padded_int(filename, fnam_det, idet, 3)
        if (n_tilts_total>1) fnam_det = trim(adjustl(fnam_det))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
        call output_stem_image(stem_image(:,:,:,:,idet),fnam_det,probe_df)
    endif
    enddo
    if(ionization) then
        do ii=1,num_ionizations
		
        if(EDX) then
            filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))&
						   &//'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        else
            filename = trim(adjustl(fnam))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))&
						  &//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)

            filename =  trim(adjustl(fnam))//'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))&
						   &//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii)*eels_correction_image, filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(fnam)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
       
    endif
    
    if (pacbed) then
#ifdef GPU
	pacbed_pattern = pacbed_pattern_d
#endif	
		do i=1,nz
			filename = trim(adjustl(output_prefix))
			if(nz>1) filename = trim(adjustl(output_prefix))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
			if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
			filename = trim(adjustl(filename))//'_PACBED_Pattern'
			call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/nysample/nxsample,filename)
		enddo
    endif
    enddo !End loop over tilts
end subroutine absorptive_stem
