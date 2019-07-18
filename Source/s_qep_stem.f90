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

function qep_stem_GPU_memory(n_qep_grates, quick_shift, phase_ramp_shift) result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly, ndet
    use m_lens, only: imaging
    use m_multislice
    
    implicit none
    
    logical,intent(in)::quick_shift,phase_ramp_shift
    integer*4,intent(in)::n_qep_grates
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(4 + n_slices + n_slices*n_qep_grates) + 2 + ndet
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
   ! if (ionization) array_count = array_count + 1 + n_slices
    !if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_stem(STEM,ionization,PACBED)
    
    use global_variables
    use m_lens
    use m_user_input
    use m_precision
    use output
#ifdef GPU	
    use cuda_array_library
    use cudafor
    use cuda_ms
    use CUFFT
    use cuda_potential
    use cuda_setup
#endif
	use cufft_wrapper
	use m_crystallography
    use m_tilt
    use m_string
    use m_multislice!, only: make_qep_grates, output_probe_intensity, output_cell_list, cell_map, output_thickness_list, setup_propagators
    use m_potential
	use m_numerical_tools
    
    implicit none
    
	logical,intent(in)::STEM,ionization,PACBED
    !dummy variables
    integer(4) ::  i,j,l,m,i_qep_pass,iz,k,ii,jj,ntilt
    integer(4) ::  count
    integer(4) ::  shifty,shiftx
    integer(4) ::  ny,nx,i_df,idet
    integer(4) ::  total_slices
    
    !random variables
    integer(4) :: idum,nsliceoutput,z_indx(1)
    !real(fp_kind) :: ran1
    complex(fp_kind) :: temp_transf(nopiy,nopix),prop(nopiy,nopix,n_slices)
       
    !probe variables
    complex(fp_kind),dimension(nopiy,nopix) :: psi,trans,psi_initial,psi_out,temp_image
     complex(fp_kind),dimension(nopiy,nopix,nz)::psi_elastic
     complex(fp_kind),dimension(nopiy,nopix,n_qep_grates,n_slices)::projected_potential,qep_grates
    !output
    real(fp_kind),dimension(nopiy,nopix) :: image,temp
    
    !STEM image variables
    real(fp_kind) :: masks(nopiy,nopix,ndet),cbed(nopiy,nopix,nz)               !detector masks
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,ndet,nz) :: stem_image,stem_elastic_image&
																  &,stem_inelastic_image
	real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz) :: eels_correction_image
    real(fp_kind),allocatable :: ion_image(:,:,:),stem_ion_image(:,:,:,:,:),pacbed_pattern(:,:,:),pacbed_elastic(:,:,:)
    
    !diagnostic variables
    real(fp_kind) :: intensity
    real(fp_kind) :: t1, delta
    
    !output variables
    character(120) :: fnam, fnam_temp, fnam_det
#ifdef GPU
    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d,psi_out_d,psi_initial_d,trans_d
	complex(fp_kind),device,dimension(nopiy,nopix,nz)::psi_elastic_d
    real(fp_kind),device :: cbed_d(nopiy,nopix,nz),temp_d(nopiy,nopix)
	complex(fp_kind),device,allocatable :: prop_d(:,:,:),transf_d(:,:,:,:)
	complex(fp_kind),device,allocatable,dimension(:,:) ::shift_arrayx_d,shift_arrayy_d,shift_array_d
    real(fp_kind),device,allocatable,dimension(:,:,:) :: masks_d,ion_image_d,pacbed_pattern_d
    real(fp_kind),device,allocatable :: eels_correction_detector_d(:,:),ion_potential_d(:,:,:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,fz_dwf_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_mu_d
    real(fp_kind),device,allocatable :: inelastic_potential_d(:,:)
#endif

    real(fp_kind),allocatable :: probe_intensity(:,:,:)
    character(1024) :: filename
    
    real(fp_kind)::qep_stem_GPU_memory,const
    logical::elfourd
    integer*4::length,lengthdf
    
#ifdef GPU    
    call GPU_memory_message(qep_stem_GPU_memory(), on_the_fly)
#endif
    
    
    call command_line_title_box('Pre-calculation setup')

    if(pacbed) then
        allocate(pacbed_pattern(nopiy,nopix,nz))
        pacbed_pattern=0
        if(output_thermal) then
            allocate(pacbed_elastic(nopiy,nopix,nz))
            pacbed_elastic = 0
        endif
	    elfourd = fourdSTEM.and.output_thermal
        length = ceiling(log10(maxval(zarray)))
    endif
    
! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
    idum = seed_rng()
#ifdef GPU	
    if (on_the_fly) then
        call cuda_setup_many_phasegrate()               !setup the atom co-ordinate for on the fly potentials (the arguments are simply so that
	else
#else
    !Generally not practical for CPU calculations
    on_the_fly = .false.
    if (.not.on_the_fly) then
#endif
        
        projected_potential = make_qep_grates(idum)
        if(.not. load_grates) then
		call load_save_add_grates(idum,projected_potential,nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice)
	endif
	
        call make_local_inelastic_potentials(ionization)          !setup the REAL space inelastic potential (ionization and adf) for QUEP ADF is disabled
    
    endif
    
    
	if(stem) then
		do i=1,ndet/nseg
		do j=1,nseg
			if(nseg>1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),2*pi*j/nseg-seg_det_offset,2*pi/nseg)
			if(nseg==1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))
		enddo
		enddo
	endif
	
    
	t1 = secnds(0.0)    
	
    call command_line_title_box('Calculation running')
	
    if(stem) stem_image = 0.0_fp_kind
    if(stem) stem_elastic_image = 0.0_fp_kind
#ifdef GPU
	! Plan the fourier transforms
    if(fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    ! Copy host arrays to the device
    if (stem.and.(ndet.gt.0)) then
        allocate(masks_d(nopiy,nopix,ndet))
        masks_d = masks
    endif

    if (ionization) then
		allocate(ion_image_d(nopiy,nopix,num_ionizations))
		if(.not.EDX) then
			allocate(eels_correction_detector_d(nopiy,nopix))
			eels_correction_detector_d=eels_correction_detector
		endif
    endif
    if(pacbed) then
		allocate(pacbed_pattern_d(nopiy,nopix,nz))
		pacbed_pattern_d=0
	endif
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    if (on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix),inverse_sinc_d(nopiy,nopix),fz_d(nopiy,nopix,nt))
        fz_d = fz 
        if(ionization) then
            allocate(inelastic_potential_d(nopiy,nopix))
            allocate(fz_mu_d(nopiy,nopix,num_ionizations))
            fz_mu_d = ionization_mu
        endif
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
  	    transf_d = qep_grates
        if (ionization) then
            allocate(ion_potential_d(nopiy,nopix,num_ionizations,n_slices))
            ion_potential_d = ionization_potential
        endif
        if (qep_mode==3) then
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d = shift_arrayx
            shift_arrayy_d = shift_arrayy
        endif
    endif
#else
	if (ionization) allocate(ion_image(nopiy,nopix,num_ionizations)) 
#endif
    if (ionization) then
        allocate(stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
		stem_ion_image = 0.0_fp_kind
		if (.not.EDX) eels_correction_image = 0.0_fp_kind
        !call binary_out_unwrap(nopiy,nopix,eels_correction_detector,'eels_correction_detector')
	endif
    if (output_probe_intensity) allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
    intensity = 1.0d0

    do ntilt=1,n_tilts_total
        
        
 	    do i = 1, n_slices
	        call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz(1),ss,ig1,ig2,claue(:,1),ifactorx,ifactory)
	        prop(:,:,i) = prop(:,:,i) * bwl_mat
            qep_grates(:,:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,:,i))
			do j=1,n_qep_grates
			call fft2(nopiy,nopix,qep_grates(:,:,j,i),nopiy,psi,nopiy)
			if(qep_mode.eq.3) qep_grates(:,:,j,i)= psi*bwl_mat
			if(qep_mode.ne.3) call ifft2(nopiy,nopix,psi*bwl_mat,nopiy,qep_grates(:,:,j,i),nopiy)
			enddo
        enddo
        
#ifdef GPU
        prop_d=prop
        if(.not.on_the_fly) then
		transf_d = qep_grates
	endif
#endif        
   lengthdf = ceiling(log10(maxval(abs(probe_df))))
   if(any(probe_df<0)) lengthdf = lengthdf+1
	do i_df = 1, probe_ndf
	do ny = 1, nysample
    do nx = 1, nxsample
#ifdef GPU
    write(6,903,ADVANCE='NO') achar(13), i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
903     format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#else
        write(6,900) i_df, probe_ndf, ny, nysample, nx, nxsample, intensity
900     format(1h+,1x,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
#endif        
    

		!Make STEM probe
		psi_initial = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
        call ifft2(nopiy, nopix, psi_initial, nopiy, psi_initial, nopiy)
        psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
        
        call tilt_wave_function(psi_initial)
        if (output_probe_intensity) probe_intensity = 0.0_fp_kind
        cbed=0_fp_kind
#ifdef GPU
        cbed_d=0.0_fp_kind
        psi_elastic_d=0.0_fp_kind
        psi_initial_d = psi_initial
        do i_qep_pass = 1, n_qep_passes 
			! Reset wavefunction
            psi_d = psi_initial_d
			if (ionization) ion_image_d=0.0_fp_kind
        
            
            do i = 1,maxval(ncells)
	            do j = 1, n_slices
                    ! Accumulate ionization cross section
                    if(ionization) then
                        
					  do ii=1,num_ionizations
							call cuda_mod<<<blocks,threads>>>(psi_d,temp_d,1.0_fp_kind,nopiy,nopix) 
							if(on_the_fly) then
								call cuda_make_ion_potential(inelastic_potential_d,tau_slice(:,atm_indices(ii),:,j),nat_slice(atm_indices(ii),j),plan,&
																&fz_mu_d(:,:,ii),inverse_sinc_d,Volume_array(j))
								call cuda_multiplication<<<blocks,threads>>>(temp_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
							else
								call cuda_multiplication<<<blocks,threads>>>(temp_d,ion_potential_d(:,:,ii,j), temp_d,prop_distance(j),nopiy,nopix)     !overlap
							endif
							call cuda_addition<<<blocks,threads>>>(temp_d,ion_image_d(:,:,ii),ion_image_d(:,:,ii),1.0_fp_kind,nopiy,nopix)                          !depth sum
						enddo
					endif
                    ! Phase grate
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    if(on_the_fly) then
                        call cuda_fph_make_potential(trans_d,ccd_slice_array(j),tau_slice,nat_slice(:,j),j,prop_distance(j),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(qep_mode == 2) then
                        shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                        shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
                        call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,j),trans_d,nopiy,nopix,shifty,shiftx)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(qep_mode == 3) then                       !randomly shift phase grate
                        shiftx = floor(ifactorx*ran1(idum)) + 1
                        shifty = floor(ifactory*ran1(idum)) + 1
                        call cuda_make_shift_array<<<blocks,threads>>>(shift_array_d,shift_arrayy_d(:,shifty),shift_arrayx_d(:,shiftx),nopiy,nopix)     !make the qspace shift array
                        call cuda_multiplication<<<blocks,threads>>>(transf_d(:,:,nran,j),shift_array_d, trans_d,1.0_fp_kind,nopiy,nopix) !multiply by the qspace shift array
                        call cufftExec(plan,trans_d,trans_d,CUFFT_INVERSE)!inverse fourier transform
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,sqrt(normalisation),nopiy,nopix)              !do the phase grate multiplication
                    else
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,nran,j), psi_out_d,1.0_fp_kind,nopiy,nopix)
                    endif
                    ! Propagate
			    	call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,j), psi_out_d,normalisation,nopiy,nopix)
                    call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)
					if (output_probe_intensity) then
						k = (i-1)*n_slices+j
						if (output_cell_list(k)) then
							psi = psi_d
							probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
						endif
					endif
		        enddo ! End loop over slices
                !If this thickness corresponds to any of the output values then accumulate diffraction pattern
				if (any(i==ncells)) then
					
					call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
					call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
					z_indx = minloc(abs(ncells-i))
					
					! Accumulate elastic wave function
					call cuda_addition<<<blocks,threads>>>(psi_elastic_d(:,:,z_indx(1)),psi_d,psi_elastic_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)

					! Accumulate diffaction pattern
					!call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
					call cuda_addition<<<blocks,threads>>>(cbed_d(:,:,z_indx(1)),temp_d,cbed_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
					
					if(ionization) then
					do ii=1,num_ionizations
						stem_ion_image(ny,nx,i_df,z_indx(1),ii) = stem_ion_image(ny,nx,i_df,z_indx(1),ii)+get_sum(ion_image_d(:,:,ii))
					enddo
					endif
				endif
            enddo ! End loop over cells
		enddo ! End loop over QEP passes
        
        intensity = get_sum(psi_d)
        cbed= cbed_d
		if (output_thermal.and.fourDSTEM) psi_elastic=psi_elastic_d
		do iz=1,nz

			! Integrate the diffraction pattern
			do idet = 1, ndet
				stem_image(ny,nx,i_df,idet,iz) = cuda_stem_detector(cbed_d(:,:,iz),masks_d(:,:,idet))
			enddo
			if (ionization.and.(.not.EDX)) eels_correction_image(ny,nx,i_df,iz) = cuda_stem_detector(cbed_d(:,:,iz),eels_correction_detector_d)
        
			! Integrate the elastic diffraction pattern
			call cufftExec(plan,psi_elastic_d(:,:,iz),psi_out_d,CUFFT_FORWARD)
			call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix) 
			do idet = 1, ndet
				stem_elastic_image(ny,nx,i_df,idet,iz)=cuda_stem_detector(temp_d,masks_d(:,:,idet))
			enddo
			psi_out = psi_out_d*sqrt(normalisation)
			
#else
		psi_elastic=0_fp_kind
			do i_qep_pass = 1, n_qep_passes 
			if (ionization) ion_image=0.0_fp_kind
            ! Reset wavefunction
            psi = psi_initial
            
            do i = 1,maxval(ncells)
	            do j = 1, n_slices
                    ! Accumulate ionization cross section
					if(ionization) then
                        do ii=1,num_ionizations
					        temp = abs(psi)**2 * ionization_potential(:,:,ii,j) * prop_distance(j)
					        ion_image(:,:,ii) = temp+ion_image(:,:,ii)
                        enddo
					endif
                    
                  call qep_multislice_iteration(psi,prop(:,:,j),qep_grates(:,:,:,j),nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,qep_mode,shift_arrayy,shift_arrayx)
                    
					if (output_probe_intensity) then
						k = (i-1)*n_slices+j
						if (output_cell_list(k)) then
							probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
						endif
					endif
		        enddo ! End loop over slices
				
                !If this thickness corresponds to any of the output values then accumulate diffration pattern
				if (any(i==ncells)) then
					z_indx = minloc(abs(ncells-i))
                    
                    ! Accumulate elastic wave function - this will be Fourier transformed later
					psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi
                    
					!Transform into diffraction space
					call fft2(nopiy,nopix,psi,nopiy,psi_out,nopiy)
                    
					! Accumulate diffaction pattern
					temp = abs(psi_out)**2
					cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp

					if(ionization) stem_ion_image(ny,nx,i_df,z_indx(1),:) = stem_ion_image(ny,nx,i_df,z_indx(1),:)+ sum(sum(ion_image,dim=2),dim=1)
                endif
                
            enddo ! End loop over cells
		enddo !End loop over QEP passes
		intensity = sum(abs(psi)**2)
        
		do iz=1,nz
			
        
			! Integrate the elastic diffraction pattern
			call fft2(nopiy,nopix,psi_elastic(:,:,iz),nopiy,psi_out,nopiy)
            if(stem) stem_image(ny,nx,i_df,1:ndet,iz) = sum(sum(spread(cbed(:,:,iz),dim=3,ncopies=ndet)*masks(:,:,:),dim=1),dim=1)
            if(stem) stem_elastic_image(ny,nx,i_df,1:ndet,iz)= sum(sum(spread(abs(psi_out)**2,dim=3,ncopies =ndet)*masks(:,:,1:ndet),dim=1),dim=1)
			if(ionization.and.(.not.EDX)) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(cbed(:,:,iz)*eels_correction_detector)
#endif            
            if(pacbed) then
			    !Output 4D STEM diffraction pattern
			    if(fourDSTEM) then
					    !Output total (elastic and inelastic) diffraction pattern
					    filename = trim(adjustl(output_prefix))
						if (probe_ndf>1) filename = trim(adjustl(filename))//defocus_string(probe_df(i_df),lengthdf)
					    if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(iz)))//'_A'
					    call binary_out_unwrap(nopiy, nopix, cbed(:,:,iz)/n_qep_passes, trim(adjustl(filename)) //'_pp_'//&
										      &to_string(nx)//'_'//to_string(ny)//'_Diffraction_pattern',write_to_screen=.false.&
                                              &,nopiyout=nopiyout,nopixout=nopixout)
			    endif

			    if(i_df==1) pacbed_pattern(:,:,iz) = pacbed_pattern(:,:,iz) + cbed(:,:,iz)/n_qep_passes
                if(output_thermal) pacbed_elastic(:,:,iz) =pacbed_elastic(:,:,iz) + abs(psi_out)**2/n_qep_passes**2

			    if (output_thermal.and.fourDSTEM) then 
					    !Output elastic only diffraction pattern
                        filename = trim(adjustl(output_prefix))
						if (probe_ndf>1) filename = trim(adjustl(filename))//defocus_string(probe_df(i_df),lengthdf)
					    if (nz>1) filename = trim(adjustl(filename))//'_z='//to_string(int(zarray(iz)))//'_A'
					    filename = trim(adjustl(filename))//'_pp_'//to_string(nx)//'_'//to_string(ny)//'_Elastic_Diffraction_pattern'
					    call binary_out_unwrap(nopiy, nopix, abs(psi_out)**2/n_qep_passes**2, filename,write_to_screen=.false.&
                                              &,nopiyout=nopiyout,nopixout=nopixout)
                endif
            endif
		enddo
			

                
		if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample)
	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	enddo ! End loop over defocus series
	
    ! QEP normalisation
    
    if(ionization) then
		stem_ion_image = stem_ion_image/float(n_qep_passes)
		if(.not.EDX) eels_correction_image = eels_correction_image/float(n_qep_passes)
	endif

    if (ndet.gt.0) then
	    stem_image = stem_image/float(n_qep_passes)
        stem_elastic_image = stem_elastic_image/float(n_qep_passes*n_qep_passes)
    endif

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
    fnam = trim(adjustl(output_prefix))
    if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
    stem_inelastic_image = stem_image - stem_elastic_image
    do idet = 1, ndet
        
        if(output_thermal) then
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneElastic_Detector'
			call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
			call output_stem_image(stem_elastic_image(:,:,:,idet,:),fnam_det,probe_df)
        
			
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneTDS_Detector'
			call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
			call output_stem_image(stem_inelastic_image(:,:,:,idet,:),fnam_det,probe_df)
        
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlaneTotal_Detector'
		else
			fnam_temp = trim(adjustl(fnam)) // '_DiffPlane_Detector'
		endif

        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_image(:,:,:,idet,:),fnam_det,probe_df)
    enddo

    !ionization
    if(ionization) then
        do ii=1,num_ionizations
        if(EDX) then
            filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_shell_EDX'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        else
            filename = trim(adjustl(fnam))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image

            filename =  trim(adjustl(fnam)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'//trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'            
            call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
        endif
        enddo
        if(.not.EDX) then
            filename = trim(adjustl(fnam)) // '_EELS_CorrectionMap' 
            call output_stem_image(eels_correction_image, filename,probe_df)
        endif
       
    endif
    
    if(pacbed) then
            const = float(nysample*nxsample)
        	do i=1,nz
            
		    filename = trim(adjustl(output_prefix))
            if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
		    if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A_'
		    call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const,trim(adjustl(filename))//'_PACBED_Pattern')
        
            if(output_thermal) then
                call binary_out_unwrap(nopiy,nopix,PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//'_elastic_PACBED_Pattern')    
                call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const-PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//'_thermal_PACBED_Pattern')
            endif
	    enddo
    endif
    enddo	
end subroutine qep_stem
