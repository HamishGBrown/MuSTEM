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
    use global_variables, only: nopiy, nopix, adf, eels
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

#ifdef GPU
    use cuda_array_library, only: cuda_multiplication, cuda_addition, cuda_mod, blocks, threads
    use cudafor, only: dim3
    use cuda_ms
    use cuda_potential
    use cufft, only: cufft_z2z, cufft_c2c, cufft_forward, cufft_inverse, cufftplan, cufftexec
    use cuda_setup, only: GPU_memory_message
#endif
    use m_crystallography; use m_multislice;use m_tilt;use m_potential;use m_string;use output;use FFTW3
	use m_Hn0; use m_absorption;use m_lens;use global_variables; use m_precision

    implicit none


    logical,intent(in)::pacbed,stem,ionization

    !dummy variables
    integer(4) :: i, j, k, l,i_df, ny, nx,z_indx(1),ii,length,idum,ntilt,lengthdf,idet,lengthimdf,starting_slice

    !probe variables
	complex(fp_kind),allocatable::propy(:),propx(:)
    complex(fp_kind),dimension(nopiy,nopix) :: psi,qpsi,psi_out,prop
	complex(fp_kind),dimension(:,:,:),allocatable::ctf,Vg,projected_potential,transf_absorptive

    !output/detectors
    real(fp_kind),dimension(nopiy,nopix) :: psi_intensity,temp
    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz,ndet) :: stem_image,stem_elastic_image&
                                                             &,stem_inelastic_image

    real(fp_kind),dimension(nysample,nxsample,probe_ndf,nz)::eels_correction_image
	real(fp_kind),allocatable :: probe_intensity(:,:,:),stem_ion_image(:,:,:,:,:),pacbed_pattern(:,:,:)&
                               &,adf_image(:,:,:),ion_image(:,:,:),masks(:,:,:),istem_image(:,:,:,:),adf_cont(:)
	real(8)::thmin,thmax


    !diagnostic variables
    real(fp_kind) :: intens, t1, delta,result

    !output variables
    character(120) :: filename,fnam,fnam_det
    logical:: many_df,manyz,dodf,manytilt,factorized_propagator

#ifdef GPU
  !device variables
  integer :: plan
  complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d!,psi_out_d

  real(fp_kind),device, allocatable,dimension(:,:) :: psi_intensity_d,eels_correction_detector_d,inelastic_potential_d,temp_d
  real(fp_kind),device,allocatable,dimension(:,:,:) :: adf_image_d,ion_image_d,pacbed_pattern_d,masks_d
  real(fp_kind),device,allocatable,dimension(:,:,:,:) :: adf_potential_d,ion_potential_d

  !device variables for on the fly potentials
  complex(fp_kind),device,allocatable,dimension(:) :: propy_d,propx_d
  complex(fp_kind),device,allocatable,dimension(:,:) :: inverse_sinc_d,trans_d,prop_d
  complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d,fz_dwf_d,fz_abs_d,fz_mu_d,transf_d
  complex(fp_kind),device,allocatable,dimension(:,:,:,:) :: fz_adf_d

  !Double channeling variables
  complex(fp_kind),device,allocatable,dimension(:,:) ::psi_inel_d,shiftarray,tmatrix_d,q_tmatrix_d
  complex(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states_d,Hn0_shifty_coord_d,Hn0_shiftx_coord_d,ctf_d
  real(fp_kind),device,allocatable,dimension(:,:)::cbed_inel_dc_d
  real(fp_kind),device,allocatable,dimension(:,:,:)::Hn0_eels_detector_d
  real(fp_kind),device,allocatable,dimension(:,:,:,:)::efistem_image_d,istem_image_d
  real(fp_kind),allocatable,dimension(:,:,:,:,:)::Hn0_eels_dc
  real(fp_kind),allocatable,dimension(:,:,:)::tmatrix_states
  integer::jj,i_target
#endif
    integer,parameter :: iu = 8945

#ifdef GPU
    real(fp_kind) :: absorptive_stem_GPU_memory

	manyz = nz>1
	manytilt = n_tilts_total>1

#ifdef GPU
    call GPU_memory_message(absorptive_stem_GPU_memory(), on_the_fly)
#else
	on_the_fly = .false.
#endif

    call command_line_title_box('Pre-calculation setup')

    ! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
    if (on_the_fly) then
#ifdef GPU
        ! Setup the atom co-ordinate for on the fly potentials
        call cuda_setup_many_phasegrate()
		if(stem) allocate(adf_cont(ndet))
#endif
    else
#ifdef GPU
        allocate(temp_d(nopiy,nopix))
#endif
        allocate(projected_potential(nopiy,nopix,n_slices),transf_absorptive(nopiy,nopix,n_slices))
		if(.not.load_grates) projected_potential = make_absorptive_grates(nopiy,nopix,n_slices)
        call load_save_add_grates(projected_potential,nopiy,nopix,n_slices)
        if(ionization.or.stem) call make_local_inelastic_potentials(ionization)
    endif


    if(pacbed) then;allocate(pacbed_pattern(nopiy,nopix,nz)); pacbed_pattern= 0;endif

#ifdef GPU
	if(double_channeling) then
        allocate(tmatrix_states_d(nopiy,nopix,nstates),psi_inel_d(nopiy,nopix),cbed_inel_dc_d(nopiy,nopix),tmatrix_states(nopiy,nopix,nstates))
		allocate(shiftarray(nopiy,nopix),tmatrix_d(nopiy,nopix),q_tmatrix_d(nopiy,nopix))
        tmatrix_states_d = setup_ms_hn0_tmatrices(nopiy,nopix,nstates)*alpha_n
        allocate(Hn0_shifty_coord_d(nopiy,maxval(natoms_slice_total),n_slices))
        allocate(Hn0_shiftx_coord_d(nopix,maxval(natoms_slice_total),n_slices))
        Hn0_shiftx_coord_d = Hn0_shiftx_coord
        Hn0_shifty_coord_d = Hn0_shifty_coord
		allocate(Hn0_eels_dc(nysample,nxsample,probe_ndf,nz,numeels))
		allocate(Hn0_eels_detector_d(nopiy,nopix,numeels),Hn0_eels_detector(nopiy,nopix,numeels))
		do l=1,numeels
			Hn0_eels_detector(:,:,l) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,0.0_fp_kind,outerrad(l))

		enddo
		Hn0_eels_detector_d = Hn0_eels_detector
		Hn0_eels_dc = 0.0_fp_kind
		if(istem) allocate(efistem_image_d(nopiy,nopix,imaging_ndf,nz));efistem_image_d=0
    endif
#endif

	if(istem) then; allocate(ctf(nopiy,nopix,imaging_ndf),istem_image(nopiy,nopix,imaging_ndf,nz));istem_image=0
			lengthimdf = calculate_padded_string_length(imaging_df,imaging_ndf)
			do i=1,imaging_ndf
				ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,imaging_aberrations,imaging_apodisation)
			enddo
#ifdef GPU
		allocate(ctf_d(nopiy,nopix,imaging_ndf));ctf_d=ctf
#endif
		endif

    many_df = probe_ndf .gt. 1
    length = calculate_padded_string_length(zarray,nz)
    lengthdf = calculate_padded_string_length(probe_df,probe_ndf)
    ! Make detector mask arrays (adds the elastic contribution)
	if(stem) then
        allocate(masks(nopiy,nopix,ndet))
        do i=1,ndet/nseg
        do j=1,nseg
            if(nseg>1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),&
                                                                            2*pi*j/nseg-seg_det_offset,2*pi/nseg)
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
	factorized_propagator = all(abs(claue)<1e-2).and.on_the_fly.and.even_slicing
#ifdef GPU
    !plan the fourier transforms
    if (fp_kind.eq.8) then
        call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
    else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif

    ! Copy host arrays to the device
    if (on_the_fly) then
        allocate(Vg(nopiy,nopix,nt))!,Vg_d(nopiy,nopix,nt)

	    Vg = fz*fz_dwf
		if(include_absorption) Vg = Vg +ci*absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm&
			&, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak*tp * ak*ss(7)/relm  !make the potential absorptive
		Vg = Vg*spread(inverse_sinc,ncopies=nt,dim=3)

		if(stem.and.adf)  then
			if(allocated(fz_adf))deallocate(fz_adf)
			allocate(fz_adf(nopiy,nopix,nt,ndet))
			do k=1,ndet/nseg
			  thmin =  atan(inner((k-1)*nseg+1)/ak)
			  thmax =  atan(outer((k-1)*nseg+1)/ak)

			  !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
			  !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
			  fz_adf(:,:,:,(k-1)*nseg+1:k*nseg) = spread(cmplx(absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss&
															&,atf,nat, ak, relm, orthog,thmin,thmax)/nseg,0.0_fp_kind)*4*pi,dim=4,ncopies=nseg)/nseg
			enddo
			fz_adf =fz_adf*spread(spread(inverse_sinc,ncopies=nt,dim=3),ncopies=ndet,dim=4)
		endif

		if(factorized_propagator) then
			allocate(propy_d(nopiy),propx_d(nopix),propy(nopiy),propx(nopix))
			call make_propagator_components(nopiy,nopix,propy,propx,prop_distance(1),ak1,ss,ig1,ig2,ifactorx,ifactory)
			propy_d = propy*bwl_mat(:,1)
			propx_d = propx*bwl_mat(1,:)
		endif

        if (ionization) then
           allocate(fz_mu_d(nopiy,nopix,num_ionizations))
            fz_mu_d = ionization_mu
        endif

        if (ionization) allocate(inelastic_potential_d(nopiy,nopix))
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
        allocate(transf_d(nopiy,nopix,n_slices))
    endif

    if(.not.factorized_propagator) allocate(prop_d(nopiy,nopix))
    if ((stem.and.adf.and.(.not.on_the_fly)).or.ionization) allocate(psi_intensity_d(nopiy,nopix))
    if (stem.and.adf.and.(.not.on_the_fly)) allocate(adf_image_d(nopiy,nopix,ndet))
    if (ionization) allocate(ion_image_d(nopiy,nopix,num_ionizations))
	  if (pacbed) allocate(pacbed_pattern_d(nopiy,nopix,nz))
    if (stem.and.(.not.on_the_fly)) allocate(masks_d(nopiy,nopix,ndet))
    if (stem.and.(.not.on_the_fly)) masks_d = masks

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

    if(even_slicing.and.(.not.factorized_propagator)) then
		call make_propagator(nopiy,nopix,prop,prop_distance(1),Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory,exponentiate = .true.)
		prop = prop*bwl_mat

	elseif(.not.factorized_propagator) then
		call make_propagator(nopiy,nopix,prop,1.0_fp_kind,Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory,exponentiate = .false.)
	endif

	if(.not.on_the_fly) then
	do i = 1, n_slices
        transf_absorptive(:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,i))
        !! Bandwith limit the phase grate, psi is used for temporary storage
        call inplace_fft(nopiy, nopix, transf_absorptive(:,:,i),norm=.true.)
        transf_absorptive(:,:,i) = transf_absorptive(:,:,i) * bwl_mat
        call inplace_ifft(nopiy, nopix, transf_absorptive(:,:,i),norm=.true.)
    enddo
    endif

	pacbed_pattern = 0
#ifdef GPU
    !if(pacbed) pacbed_pattern_d = 0
	if(.not.on_the_fly) transf_d = transf_absorptive
	if(.not.factorized_propagator) prop_d = prop
#endif
    do ny = 1, nysample
    do nx = 1, nxsample
    do i_df = 1, probe_ndf                                                                  !loop over probe points
#ifdef GCC

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
        call inplace_ifft(nopiy, nopix, psi)
        psi = psi/sqrt(sum(abs(psi)**2))

        call tilt_wave_function(psi)

#ifdef GPU
        if (adf.and.(.not.on_the_fly)) adf_image_d = 0.0_fp_kind
		if(stem.and.adf.and.on_the_fly) adf_cont=0
        if (ionization) ion_image_d = 0.0_fp_kind
        psi_d = psi
        do i = 1,maxval(ncells);do j = 1, n_slices

			! Calculate inelastic cross sections
			if ((stem.and.adf.and.(.not.on_the_fly)).or.ionization) call cuda_mod<<<blocks,threads>>>(psi_d,psi_intensity_d,1.0_fp_kind,nopiy,nopix)
			if((stem.and.adf)) then
				do k=1,ndet

				if(on_the_fly) then
					adf_cont(k) = adf_cont(k) + cuda_adf_crossection_on_the_fly(psi_d,tau_slice(:,:,:,j),nat_slice(:,j),plan,fz_adf(:,:,:,k),Volume_array(j),prop_distance(j))
				else
					call cuda_multiplication<<<blocks,threads>>>(psi_intensity_d,adf_potential_d(:,:,j,k), temp_d,prop_distance(j),nopiy,nopix)     !overlap
					call cuda_addition<<<blocks,threads>>>(temp_d,adf_image_d(:,:,k),adf_image_d(:,:,k),1.0_fp_kind,nopiy,nopix)                          !depth sum
				endif
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

			!Double channeling
			if(double_channeling) then

                do i_target = 1, natoms_slice_total(j) ! Loop over targets

                    ! Calculate inelastic transmission matrix
					call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,j),Hn0_shiftx_coord_d(:,i_target,j),nopiy,nopix)

					do k = 1, nstates

					call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
                    call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
					call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,sqrt(normalisation),nopiy,nopix)
					starting_slice=j
                        ! Scatter the inelastic wave through the remaining cells
                        do ii = i, n_cells
                            do jj = starting_slice, n_slices;call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,jj),prop_d,prop_distance(jj),even_slicing, normalisation, nopiy, nopix,plan);enddo
							starting_slice=1
							if (any(ii==ncells)) then
								z_indx = minloc(abs(ncells-ii))
								call cufftExec(plan, psi_inel_d, q_tmatrix_d, CUFFT_FORWARD)

								! Accumulate the EELS images
								do l=1,numeels
									Hn0_eels_dc(ny,nx,i_df,z_indx(1),l) = Hn0_eels_dc(ny,nx,i_df,z_indx(1),l)+cuda_stem_detector(q_tmatrix_d, Hn0_eels_detector_d(:,:,l))
								enddo
								! Accumulate EFISTEM images
								if (istem.and.i_df==1) then;do l = 1, imaging_ndf
									call cuda_image(q_tmatrix_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
									call cuda_addition<<<blocks,threads>>>(efistem_image_d(:,:,l,z_indx(1)), temp_d, efistem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
								enddo;endif;
								!stop
							endif
              enddo
						enddo
				enddo
			endif  ! End loop over cells,targets and states and end double_channeling section
			! Transmit through slice potential
			if(on_the_fly) then
				if(.not.factorized_propagator) call cuda_on_the_fly_abs_multislice(psi_d,ccd_slice_array(j),tau_slice(:,:,:,j),nat_slice(:,j),prop_distance(j),even_slicing,plan,Vg,Volume_array(j),prop_d = prop_d)
				if(factorized_propagator) call cuda_on_the_fly_abs_multislice(psi_d,ccd_slice_array(j),tau_slice(:,:,:,j),nat_slice(:,j),prop_distance(j),even_slicing,plan,Vg,Volume_array(j),propy_d, propx_d)
            else
				call cuda_multislice_iteration_new(psi_d, transf_d(:,:,j), prop_d,prop_distance(j),even_slicing, normalisation, nopiy, nopix,plan)
            endif
			if (output_probe_intensity) then
				k = (i-1)*n_slices+j
				if (output_cell_list(k)) then
					psi = psi_d
					probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
				endif
			endif
		enddo;
		!!If this thickness corresponds to any of the output values then output images
		!	if (any(i==ncells)) then
		!
		!		z_indx = minloc(abs(ncells-i))
		!
		!		call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
		!		call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
		!		!if(pacbed.and.(i_df==1)) call cuda_addition<<<blocks,threads>>>(temp_d,pacbed_pattern_d(:,:,z_indx(1)),pacbed_pattern_d(:,:,z_indx(1)),1.0_fp_kind,nopiy,nopix)
		!		do ii=1,ndet
		!			if(stem.and.nseg>1) stem_elastic_image(ny,nx,i_df,z_indx(1),ii) = cuda_stem_detector(temp_d,masks_d(:,:,ii))
		!			if(stem.and.nseg<2) then
		!			write(*,*) 'hi'
		!			psi=psi_d
		!			write(*,*) 'by'
		!			call cuda_stem_detector_wavefunction_on_the_fly<<<blocks,threads>>>(psi_out_d, result,inner(ii),outer(ii),ifactory*a0(1),ifactorx*a0(2),nopiy,nopix)
		!			write(*,*) 'hi'
		!			write(*,*) result,inner(ii),outer(ii),ifactory*a0(1),ifactorx*a0(2),nopiy,nopix
		!			psi=psi_out_d
		!			call binary_out(nopiy,nopix,psi,'psi')
		!			stop
		!			write(*,*) 'my name is'
		!			stem_elastic_image(ny,nx,i_df,z_indx(1),ii) = result
		!
		!			endif
		!			if(stem.and.adf.and.on_the_fly) stem_inelastic_image(ny,nx,i_df,z_indx(1),ii) = adf_cont(ii)
		!			if(stem.and.adf.and.(.not.on_the_fly)) stem_inelastic_image(ny,nx,i_df,z_indx(1),ii) = get_sum(adf_image_d(:,:,ii))
		!		enddo
		!		write(*,*) ionization
		!		if(ionization) then
		!			do ii=1,num_ionizations
		!				 stem_ion_image(ny,nx,i_df,z_indx(1),ii) = get_sum(ion_image_d(:,:,ii))
		!			enddo
		!			if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = cuda_stem_detector(temp_d,eels_correction_detector_d)
		!		endif
		!
		!		write(*,*) istem.and.i_df==1
		!		if (istem.and.i_df==1) then;do l = 1, imaging_ndf
		!			call cuda_image(psi_out_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)
		!			call cuda_addition<<<blocks,threads>>>(istem_image_d(:,:,l,z_indx(1)), temp_d, istem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
		!		enddo;endif;
		!		!Output 4D STEM diffraction pattern
		!		write(*,*) fourDSTEM
		!		if(fourDSTEM) then
		!				temp = temp_d
		!				call output_TEM_result(output_prefix,temp,'Diffraction_pattern',nopiy,nopix,manyz,many_df,manytilt,zarray(z_indx(1)),probe_df(i_df)&
		!				&,length,lengthdf,tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),nopiyout=nopiyout,nopixout=nopixout,pp=[ny,nx],write_to_screen=.false.)
		!		endif
	   !
		!		write(*,*) 'my name is'
		!		write(*,*) 'hi'
		!			psi=psi_d
		!			write(*,*) 'by'
		!	endif
		!enddo ! End loop over cells
		!write(*,*) 'hi'
       !intens = get_sum(psi_d)
		!write(*,*) 'my name is'
#else
        if (adf.and.stem) adf_image = 0.0_fp_kind
        if (ionization) ion_image = 0.0_fp_kind
        do i = 1,maxval(ncells);do j = 1, n_slices
			! Calculate inelastic cross sections
			if ((stem.and.adf).or.ionization) psi_intensity = abs(psi)**2


			if(stem.and.adf.and.on_the_fly) then

			elseif(stem.and.adf) then
        do k=1,ndet
          adf_image(:,:,k) = psi_intensity*adf_potential(:,:,j,k)*prop_distance(j)&
                             & + adf_image(:,:,k)
        enddo
      endif

			if(ionization) ion_image = ion_image + spread(psi_intensity,ncopies=num_ionizations,dim=3) &
                              &* ionization_potential(:,:,:,j) * prop_distance(j)
			call multislice_iteration(psi,prop,transf_absorptive(:,:,j),prop_distance(j),even_slicing,nopiy,nopix)

			if (output_probe_intensity) then
				k = (i-1)*n_slices+j
				if (output_cell_list(k)) then
					probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
				endif
			endif

      enddo !End loop over slices
#endif
			!If this thickness corresponds to any of the output values then output images
			if (any(i==ncells)) then
#ifdef GPU
				psi=psi_d
				if (adf.and.(.not.on_the_fly)) adf_image= adf_image_d
#endif
				z_indx = minloc(abs(ncells-i))

				temp = abs(fft(nopiy,nopix,psi,norm=.true.))**2
				if(stem) then
                    do k=1,ndet
                        stem_elastic_image(ny,nx,i_df,z_indx(1),k) = sum(masks(:,:,k)*temp)
                    enddo
                endif
				if(pacbed.and.(i_df==1)) pacbed_pattern(:,:,z_indx(1)) = pacbed_pattern(:,:,z_indx(1)) + temp
				if(stem.and.adf.and.on_the_fly) then
#ifdef GPU
					do k=1,ndet; stem_inelastic_image(ny,nx,i_df,z_indx(1),k) = adf_cont(k); enddo
#endif
				elseif(stem.and.adf) then
                    do k=1,ndet; stem_inelastic_image(ny,nx,i_df,z_indx(1),k) = sum(adf_image(:,:,k)); enddo
                endif
#ifdef GPU
				if (istem.and.i_df==1) then;do l = 1, imaging_ndf
					istem_image(:,:,l,z_indx(1)) = istem_image(:,:,l,z_indx(1)) + make_image(nopiy,nopix,psi_out,ctf(:,:,l),.false.)
				enddo;endif;
#endif
				if(ionization) then
					do ii=1,num_ionizations
						stem_ion_image(ny,nx,i_df,z_indx(1),ii) = sum(ion_image(:,:,ii))
          			enddo
					if(.not.EDX) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(temp*eels_correction_detector)
				endif

				!Output 4D STEM diffraction pattern
				if(fourDSTEM)  call output_TEM_result(output_prefix,temp,'Diffraction_pattern',nopiy,nopix,manyz,many_df,manytilt,&
													  &zarray(z_indx(1)),probe_df(i_df),length,lengthdf,tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),&
													  &nopiyout=nopiyout,nopixout=nopixout,pp=[ny,nx],write_to_screen=.false.)
			endif
		enddo ! End loops over cells
#ifdef GPU
		intens = get_sum(psi_d)
#else
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
        write(9834, '(a, g9.5, a, /)') 'The multislice calculation took ', delta, 'seconds.'
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
		do i=1,nz
			pacbed_pattern(:,:,i) = pacbed_pattern(:,:,i)/nysample/nxsample
			call output_TEM_result(output_prefix,pacbed_pattern(:,:,i),'PACBED',nopiy,nopix,manyz,.false.,manytilt,zarray(i)&
								&,lengthz=length,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),nopiyout=nopiyout,nopixout=nopixout)
		enddo
    endif

	if(istem) then
		do i=1,nz;do l=1,imaging_ndf
			temp = istem_image(:,:,l,i)
			call output_TEM_result(output_prefix,tile_out_image(temp/nysample/nxsample,ifactory,ifactorx),'ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i)&
								&,lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(l))
		enddo;enddo;endif
#ifdef GPU
			if(double_channeling.and.istem) then
			do i=1,nz;do l=1,imaging_ndf
			temp = efistem_image_d(:,:,l,i)
			call output_TEM_result(output_prefix,tile_out_image(temp/nysample/nxsample,ifactory,ifactorx),'energy_filtered_ISTEM',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(i)&
								&,lengthz=length,lengthdf=lengthimdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(l))
		enddo;enddo;endif
#endif
if(double_channeling) then
	do l=1,numeels
		filename =  trim(adjustl(fnam))//'_double_channeling_EELS_'//zero_padded_int(l,2)
		call output_stem_image(Hn0_eels_dc(:,:,:,:,l), filename,probe_df)
	enddo
	endif
#endif

    enddo !End loop over tilts
end subroutine absorptive_stem
