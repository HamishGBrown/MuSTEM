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
    use m_slicing, only: n_slices
    
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
    use global_variables
    use m_lens, only: probe_df, probe_ndf, make_stem_wfn,probe_aberrations
    use m_absorption, only: fz_abs, transf_absorptive, calculate_absorption_mu, calculate_local_adf_mu
    use m_slicing
    use m_probe_scan, only: nysample, nxsample, probe_positions, scan_quarter
    use m_tilt, only: tilt_wave_function
    use m_multislice
    use m_potential
    use m_string
	use output
	use cufft_wrapper
    
    implicit none
    
    logical,intent(in)::pacbed,stem,ionization
    
    !dummy variables
    integer(4) :: i,ii, j, k, i_df, ny, nx,z_indx(1),nopiyout,nopixout,length

    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,qpsi,psi_out

    !output/detectors
    real(fp_kind),dimension(nopiy,nopix) :: image,masks,temp_image,psi_intensity,&
                                            & inelastic_potential,temp
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz) :: stem_image,stem_elastic_image,&
                & stem_inelastic_image,eels_correction_image
    real(fp_kind),allocatable::adf_image(:,:),probe_intensity(:,:,:),ion_image(:,:,:),stem_ion_image(:,:,:,:,:)&
                              &,pacbed_pattern(:,:,:)
    logical::fourdstem,many_df
    
    !diagnostic variables
    real(fp_kind) :: intens, t1, delta
    
    !output variables
    character(120) :: filename

    integer,parameter :: iu = 8945
    
    real(fp_kind) :: absorptive_stem_GPU_memory
    
    write(*,*) '|----------------------------------|'
    write(*,*) '|      Pre-calculation setup       |'
    write(*,*) '|----------------------------------|'
    write(*,*)
    
    !call calculate_absorption_mu        
    call calculate_local_adf_mu
    
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    if (on_the_fly) then
        
        ! These make_fz_*() routines would normally be called in make_local_inelastic_potentials()
        ! But for on_the_fly we don't call that routine.
        if(adf.and.stem) call make_fz_adf
        !This step is now done automatically
        !if(ionization) call make_fz_EELS_EDX
        
    else
        call make_absorptive_grates
        call make_local_inelastic_potentials(ionization)
        
    endif
    
    call setup_propagators
    fourdSTEM= .false.
    if(pacbed) then
        allocate(pacbed_pattern(nopiy,nopix,nz))
        pacbed_pattern= 0
        call fourD_STEM_options(fourdSTEM,nopiyout,nopixout,nopiy,nopix)
    endif
    many_df = probe_ndf .gt. 1
    length = ceiling(log10(maxval(zarray)))
    
    ! Make detector mask arrays (adds the elastic contribution)        
    if(stem) then
        call make_detector_mask(inner(1),outer(1),masks)
        stem_elastic_image = 0.0_fp_kind
        if (adf) then
		    allocate(adf_image(nopiy,nopix))
		    stem_inelastic_image = 0.0_fp_kind
        endif
    endif
    
    if (ionization) then
        allocate(ion_image(nopiy,nopix,num_ionizations))
        allocate(stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
        stem_ion_image = 0.0_fp_kind
        if (.not.EDX) eels_correction_image = 0.0_fp_kind
    endif
    

    write(*,*) '|--------------------------------|'
    write(*,*) '|      Calculation running       |'
    write(*,*) '|--------------------------------|'
    write(*,*)

    t1 = secnds(0.0)
    
    intens = 1.0_fp_kind

    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
        
    
    do ny = 1, nysample
    do nx = 1, nxsample
    
#ifdef gpu       
    write(6, 901, advance='no') achar(13), ny, nysample, nx, nxsample, intens
901 format(a1, 1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)	
#else
        write(6,900) ny, nysample, nx, nxsample, intens
900     format(1h+,1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)	
#endif
    
    do i_df = 1, probe_ndf                                                                  !loop over probe points
        
        
        if (adf.and.stem) adf_image = 0.0_fp_kind
        if (ionization) ion_image = 0.0_fp_kind
        if (output_probe_intensity) probe_intensity = 0_fp_kind
        
        call make_stem_wfn(psi, probe_df(i_df), probe_positions(:,ny,nx),probe_aberrations)  
        call tilt_wave_function(psi)
                   
        do i = 1,maxval(ncells)
            
            do j = 1, n_slices
            
                ! Calculate inelastic cross sections
            
                if ((stem.and.adf).or.ionization) psi_intensity = abs(psi)**2
                
                if(stem.and.adf) adf_image = psi_intensity*adf_potential(:,:,j)*prop_distance(j) + adf_image
                if(ionization) ion_image = ion_image + spread(psi_intensity,ncopies=num_ionizations,dim=3) * ionization_potential(:,:,:,j) * prop_distance(j)
                call multislice_iteration(psi,prop(:,:,j),transf_absorptive(:,:,j),nopiy,nopix)
                
				!If this thickness corresponds to any of the output values then output images
				if (any(i==ncells)) then
                    z_indx = minloc(abs(ncells-i))
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
					temp = abs(psi_out)**2
					if(stem) stem_elastic_image(ny,nx,i_df,z_indx(1)) = sum(masks*temp)
                    if(pacbed.and.(i_df==1)) pacbed_pattern(:,:,z_indx(1)) = pacbed_pattern(:,:,z_indx(1)) + temp
					if((stem.and.adf)) stem_inelastic_image(ny,nx,i_df,z_indx(1)) = sum(adf_image)
					if(ionization) then
                        do ii=1,num_ionizations
                            stem_ion_image(ny,nx,i_df,z_indx(1),ii) = sum(ion_image(:,:,ii))
                        enddo
                        if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(temp*eels_correction_detector)
                    endif
					
                    !Output 4D STEM diffraction pattern
					if(fourDSTEM) then
							filename = trim(adjustl(output_prefix))
							if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(z_indx(1))))//'_A'
                            if(many_df) filename = trim(adjustl(filename)) // '_Defocus_'//zero_padded_int(int(probe_df(i_df)),length)//'_Ang'
							filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_abs_Diffraction_pattern'
							call binary_out_unwrap(nopiy, nopix, temp, filename,write_to_screen=.false.,nopiyout=nopiyout,nopixout=nopixout)
					endif
                    
                endif
    			if (output_probe_intensity) then
					k = (i-1)*n_slices+j
					if (output_cell_list(k)) then
						probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
					endif
				endif
            enddo ! End loop over slices
        enddo ! End loop over cells
                
        intens = sum(abs(psi)**2)
        
        
        if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,1)
    enddo                                                                         !end loop over defocus
    enddo                                                                         !end loop over probe position ny
    enddo                                                                         !end loop over probe position nx
    if (output_probe_intensity) close(iu)
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
   
    if(output_thermal) then    
        filename = trim(adjustl(output_prefix)) // '_DiffPlaneElastic'
        if(stem) call output_stem_image(stem_elastic_image,filename,probe_df)
    
        if (stem.and.adf) then
            stem_image = stem_elastic_image + stem_inelastic_image
    
            filename = trim(adjustl(output_prefix)) // '_DiffPlaneTotal'
            call output_stem_image(stem_image,filename,probe_df)
    
            filename = trim(adjustl(output_prefix)) // '_DiffPlaneTDS'
            call output_stem_image(stem_inelastic_image,filename,probe_df)
        endif   
    elseif(stem) then
        if (adf) stem_image = stem_elastic_image + stem_inelastic_image
        if (.not.adf) stem_image = stem_elastic_image
        filename = trim(adjustl(output_prefix)) // '_DiffPlane'
        call output_stem_image(stem_image,filename,probe_df)
    endif
    
    
    if(ionization) then
        do ii=1,num_ionizations
        if(EDX) then
            filename = trim(adjustl(output_prefix)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        else
            filename = trim(adjustl(output_prefix))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image

            filename =  trim(adjustl(output_prefix)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(output_prefix)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
    endif
    
    if (pacbed) then
		do i=1,nz
			filename = trim(adjustl(output_prefix))
			if(nz>1) filename = trim(adjustl(output_prefix))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
			filename = trim(adjustl(filename))//'_PACBED_Pattern'
			call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i),filename)
		enddo
    endif    

end subroutine absorptive_stem
