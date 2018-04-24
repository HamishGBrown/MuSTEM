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
    use global_variables, only: nopiy, nopix, adf, ionization, eels,zarray,nz,ncells
    use m_slicing, only: n_slices
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 6 + 4*n_slices
    if (adf.or.ionization) array_count = array_count + 1
    if (adf) array_count = array_count + 1 + n_slices
    if (ionization) array_count = array_count + 1 + n_slices
    if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
end function
    


subroutine absorptive_stem
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, normalisation, nt, n_cells, fz, fz_dwf, inverse_sinc, bwl_mat, prop, ndet, ionization, eels, adf, on_the_fly, inner, outer, ci,zarray,nz,ncells
    use m_lens!, only: probe_df, probe_ndf, make_stem_wfn
    use m_absorption, only: fz_abs, transf_absorptive, calculate_absorption_mu, calculate_local_adf_mu
    use output, only: output_prefix, output_stem_image,timing
    use cuda_array_library, only: cuda_multiplication, cuda_addition, cuda_mod, blocks, threads
    use cudafor, only: dim3
    use cuda_ms, only: cuda_stem_detector, get_sum
    use cuda_potential, only: volume_array, ccd_slice_array, cuda_setup_many_phasegrate, cuda_make_adf_potential, cuda_make_ion_potential, cuda_make_abs_potential
    use cufft, only: cufft_z2z, cufft_c2c, cufft_forward, cufft_inverse, cufftplan, cufftexec
    use local_ionization
    use m_slicing, only: n_slices, nat_slice, tau_slice, prop_distance
    use cuda_setup, only: GPU_memory_message
    use m_probe_scan, only: nysample, nxsample, probe_positions, scan_quarter
    use m_tilt, only: tilt_wave_function
    use m_multislice!, only: make_absorptive_grates, output_probe_intensity, output_cell_list, setup_propagators
    use m_potential, only: precalculate_scattering_factors
    use m_string, only: to_string
	use output
	use cufft_wrapper
    
    implicit none
    
    !dummy variables
    integer(4) :: i, j, k, i_df, ny, nx,z_indx(1),ii

    !probe variables
    complex(fp_kind) :: psi(nopiy,nopix),qpsi(nopiy,nopix)

    !output/detectors
    real(fp_kind),dimension(nopiy,nopix) :: image,ion_image,masks,temp_image
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz) :: stem_image,stem_elastic_image&
															 &,stem_inelastic_image,eels_correction_image
	real(fp_kind),allocatable :: probe_intensity(:,:,:),stem_ion_image(:,:,:,:,:)

    !diagnostic variables
    real(fp_kind) :: intens, t1, delta
    
    !output variables
    character(120) :: filename

    !device variables
    integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d,psi_out_d
    real(fp_kind),device,dimension(nopiy,nopix) :: temp_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: prop_d,transf_d
    real(fp_kind),device, allocatable,dimension(:,:) :: masks_d,adf_image_d,psi_intensity_d,eels_correction_detector_d,inelastic_potential_d
    real(fp_kind),device,allocatable :: adf_potential_d(:,:,:),ion_potential_d(:,:,:,:),ion_image_d(:,:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,trans_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_dwf_d,fz_abs_d,fz_adf_d,fz_mu_d
    
    integer,parameter :: iu = 8945
    
    real(fp_kind) :: absorptive_stem_GPU_memory
    
    
    
    call GPU_memory_message(absorptive_stem_GPU_memory(), on_the_fly)

    
    
    write(*,*) '|----------------------------------|'
    write(*,*) '|      Pre-calculation setup       |'
    write(*,*) '|----------------------------------|'
    write(*,*)
	
    call calculate_local_adf_mu
    
    ! Precalculate the scattering factors on a grid
    !call precalculate_scattering_factors()

    if (on_the_fly) then
        ! Setup the atom co-ordinate for on the fly potentials
        call cuda_setup_many_phasegrate()                   
        
        ! These make_fz_*() routines would normally be called in make_local_inelastic_potentials()
        ! But for on_the_fly we don't call that routine.
        if(adf) call make_fz_adf
        !if(ionization) call make_fz_EELS_EDX
        
    else
        call make_absorptive_grates
        call make_local_inelastic_potentials
    endif
    
    call setup_propagators
    
    ! Make detector mask arrays (adds the elastic contribution)        
    call make_detector_mask(inner(1),outer(1),masks)
    
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
        fz_abs_d = ci*fz_abs  !make the potential absorptive
        if (ionization) then
           allocate(fz_mu_d(nopiy,nopix,num_ionizations))
            fz_mu_d = ionization_mu
        endif
        
        if (adf) then
            allocate(fz_adf_d(nopiy,nopix,nt))
            fz_adf_d = fz_adf
        endif
        
        if (ionization.or.adf) allocate(inelastic_potential_d(nopiy,nopix))
        
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
        
    else
        ! Set up device variables for precalculated potentials
    
        allocate(transf_d(nopiy,nopix,n_slices))
        transf_d = transf_absorptive
        
        if (adf) then
            allocate(adf_potential_d(nopiy,nopix,n_slices))
            adf_potential_d = adf_potential
        endif
        
        if(ionization) then
            allocate(ion_potential_d(nopiy,nopix,num_ionizations,n_slices))
            ion_potential_d = ionization_potential
        endif
        
    endif
    
    if (adf.or.ionization) allocate(psi_intensity_d(nopiy,nopix))
    if (adf) allocate(adf_image_d(nopiy,nopix))
    if (ionization) allocate(ion_image_d(nopiy,nopix,num_ionizations),stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
    
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    
    allocate(masks_d(nopiy,nopix))
    masks_d = masks
    
    if(ionization.and.(.not.EDX)) then
        allocate(eels_correction_detector_d(nopiy,nopix))
        eels_correction_detector_d = eels_correction_detector
		eels_correction_image = 0.0_fp_kind
    endif

    stem_elastic_image = 0.0_fp_kind
    if (adf) stem_inelastic_image = 0.0_fp_kind
    if (ionization) stem_ion_image = 0.0_fp_kind

    write(*,*) '|--------------------------------|'
    write(*,*) '|      Calculation running       |'
    write(*,*) '|--------------------------------|'
    write(*,*)

    t1 = secnds(0.0)
    
    intens = 1.0_fp_kind

    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))

    do ny = 1, nysample
    do nx = 1, nxsample
    
    write(6, 901, advance='no') achar(13), ny, nysample, nx, nxsample, intens
901 format(a1, 1x, 'y:', i3, '/', i3, ' x:', i3, '/', i3, ' Intensity: ', f12.6)

    if (scan_quarter) then
        if (probe_positions(1,ny,nx).gt.0.501_fp_kind .or. probe_positions(2,ny,nx).gt.0.501_fp_kind) cycle
    endif
    
    do i_df = 1, probe_ndf                                                                  !loop over probe points
    
        if (adf) adf_image_d = 0.0_fp_kind
        if (ionization) ion_image_d = 0.0_fp_kind
        if (output_probe_intensity) probe_intensity = 0_fp_kind

        call make_stem_wfn(psi, probe_df(i_df), probe_positions(:,ny,nx),probe_aberrations) 
		
        call tilt_wave_function(psi)

        psi_d = psi
        
        do i = 1,maxval(ncells)
            
            do j = 1, n_slices
            
                ! Calculate inelastic cross sections
            
                if (adf.or.ionization) then
                    call cuda_mod<<<blocks,threads>>>(psi_d,psi_intensity_d,1.0_fp_kind,nopiy,nopix)
                endif
                
                if(adf) then
                    if(on_the_fly) then
                        call cuda_make_adf_potential(inelastic_potential_d,tau_slice(:,:,:,j),nat_slice(:,j),plan,fz_adf_d,inverse_sinc_d,Volume_array(j))
                        call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
                    else
                        call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,adf_potential_d(:,:,j), temp_d,prop_distance(j),nopiy,nopix)     !overlap
                    endif
                    call cuda_addition<<<blocks,threads>>>(temp_d,adf_image_d,adf_image_d,1.0_fp_kind,nopiy,nopix)                          !depth sum
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
                
				!If this thickness corresponds to any of the output values then output images
				if (any(i==ncells)) then
					call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
					call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
					z_indx = minloc(abs(ncells-i))
					stem_elastic_image(ny,nx,i_df,z_indx(1)) = cuda_stem_detector(temp_d,masks_d)
					if(adf) stem_inelastic_image(ny,nx,i_df,z_indx(1)) = get_sum(adf_image_d)
					if(ionization) then
						do ii=1,num_ionizations
							 stem_ion_image(ny,nx,i_df,z_indx(1),ii) = get_sum(ion_image_d(:,:,ii))
						enddo
						if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = cuda_stem_detector(temp_d,eels_correction_detector_d)
					endif
				endif
				if (output_probe_intensity) then
					k = (i-1)*n_slices+j
					if (output_cell_list(k)) then
						psi = psi_d
						probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
					endif
				endif
            enddo ! End loop over slices
            
        enddo ! End loop over cells
        
        intens = get_sum(psi_d)
		
		if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,1)        
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
   
    if(output_thermal) then
    filename = trim(adjustl(output_prefix)) // '_DiffPlaneElastic'
    call output_stem_image(stem_elastic_image,filename,probe_df)
    
    if (adf) then
        stem_image = stem_elastic_image + stem_inelastic_image
    
        filename = trim(adjustl(output_prefix)) // '_DiffPlaneTotal'
        call output_stem_image(stem_image,filename,probe_df)
    
        filename = trim(adjustl(output_prefix)) // '_DiffPlaneTDS'
        call output_stem_image(stem_inelastic_image,filename,probe_df)
    endif   
    else
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

            filename =  trim(adjustl(output_prefix)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii)*eels_correction_image, filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(output_prefix)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
       
    endif

end subroutine absorptive_stem
