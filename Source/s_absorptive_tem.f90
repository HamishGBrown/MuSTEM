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
	use FFTW3
#ifdef GPU
    use CUFFT
    use cuda_array_library
    use cudafor
    use cuda_ms
    use cuda_potential
    use cuda_setup
#endif
    use m_tilt;use m_multislice;use m_potential;use m_string;use m_hn0
    implicit none

    !dummy variables
    integer(4) :: i_cell,i_slice,z_indx(1),lengthz,lengthdf,i,idum,ntilt,starting_slice

    !random variables
    integer(4) :: count

    !probe variables
	complex(fp_kind),allocatable::Vg(:,:,:),propy(:),propx(:)
    complex(fp_kind),dimension(nopiy,nopix)  :: psi,psi_initial,psi_out,prop
    complex(fp_kind) :: lens_ctf(nopiy,nopix,imaging_ndf)
    complex(fp_kind),dimension(nopiy,nopix,n_slices)::projected_potential,transf_absorptive

    !output
    real(fp_kind),dimension(nopiy,nopix) :: cbed,image,tem_image,temp_image

    !diagnostic variables
    real(fp_kind) :: intensity, t1, delta

    !output variables
    character(120) :: filename,fnam_df
	logical::manyz,manytilt,factorized_propagator

#ifdef GPU
    !device variables
	integer :: plan
	complex(fp_kind),device,allocatable :: prop_d(:,:),transf_d(:,:,:),propy_d(:),propx_d(:)
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d
	real(fp_kind),device,dimension(:,:),allocatable::temp_d

    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: trans_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: Vg_d

	!Double channeling variables
	complex(fp_kind),device,allocatable,dimension(:,:) ::psi_inel_d,shiftarray,tmatrix_d,q_tmatrix_d
	complex(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states_d,Hn0_shifty_coord_d,Hn0_shiftx_coord_d,ctf_d
	real(fp_kind),device,allocatable,dimension(:,:)::cbed_inel_dc_d
	real(fp_kind),device,allocatable,dimension(:,:,:)::tmatrix_states
	real(fp_kind),device,allocatable,dimension(:,:,:,:)::eftem_image_d
	integer::jj,l,i_target,j,k,ii,i_df

    real(fp_kind) :: absorptive_tem_GPU_memory


    call GPU_memory_message(absorptive_tem_GPU_memory(), on_the_fly)
#else
	on_the_fly=.false.
#endif
    manyz = nz>1
	manytilt = n_tilts_total>1
    lengthz = calculate_padded_string_length(zarray,nz)
    if (pw_illum) lengthdf = calculate_padded_string_length(imaging_df,imaging_ndf)

    call command_line_title_box('Pre-calculation setup')
	do i=1,imaging_ndf
    if(pw_illum) then
      lens_ctf(:,:,i) =  make_ctf(real([0,0,0],kind=fp_kind),imaging_df(i),&
                        &imaging_cutoff,imaging_aberrations,imaging_apodisation)
    endif
	enddo
    ! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()

    if (pw_illum) then
        psi_initial = 1.0_fp_kind/sqrt(float(nopiy*nopix))
	else
        psi_initial = make_ctf(probe_initial_position,probe_df(1),probe_cutoff,probe_aberrations,probe_apodisation)
        call inplace_ifft(nopiy, nopix, psi_initial)
        psi_initial = psi_initial/sqrt(sum(abs(psi_initial)**2))
    endif
    call tilt_wave_function(psi_initial)
    if(.not.load_grates.and.(.not.on_the_fly)) projected_potential = make_absorptive_grates(nopiy,nopix,n_slices)
    if(.not.on_the_fly) call load_save_add_grates(projected_potential,nopiy,nopix,n_slices)


	t1 = secnds(0.0)
#ifdef GPU
	! Plan the fourier transforms
    if (fp_kind.eq.8)then
          call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
          call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
	write(*,*) double_channeling
	if(double_channeling) then
        allocate(tmatrix_states_d(nopiy,nopix,nstates),psi_inel_d(nopiy,nopix),cbed_inel_dc_d(nopiy,nopix),tmatrix_states(nopiy,nopix,nstates))
		allocate(shiftarray(nopiy,nopix),tmatrix_d(nopiy,nopix),q_tmatrix_d(nopiy,nopix),temp_d(nopiy,nopix))
        tmatrix_states_d = setup_ms_hn0_tmatrices(nopiy,nopix,nstates)*alpha_n
        allocate(Hn0_shifty_coord_d(nopiy,maxval(natoms_slice_total),n_slices))
        allocate(Hn0_shiftx_coord_d(nopix,maxval(natoms_slice_total),n_slices))
        Hn0_shiftx_coord_d = Hn0_shiftx_coord
        Hn0_shifty_coord_d = Hn0_shifty_coord
		allocate(eftem_image_d(nopiy,nopix,imaging_ndf,nz),ctf_d(nopiy,nopix,imaging_ndf),temp_d(nopiy,nopix)); ctf_d = lens_ctf;
		eftem_image_d=0
    endif

    !Copy host arrays to the device
	factorized_propagator = all(abs(claue)<1e-2).and.on_the_fly.and.even_slicing
    if (on_the_fly) then
		allocate(Vg(nopiy,nopix,nt))!,Vg_d(nopiy,nopix,nt)trans_d(nopiy,nopix),

	    Vg = fz*fz_dwf
		if(include_absorption) Vg = Vg +ci*absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm&
			&, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak*tp * ak*ss(7)/relm  !make the potential absorptive
		Vg = Vg*spread(inverse_sinc,ncopies=nt,dim=3)
		if(factorized_propagator) then
			allocate(propy_d(nopiy),propx_d(nopix),propy(nopiy),propx(nopix))
			call make_propagator_components(nopiy,nopix,propy,propx,prop_distance(1),ak1,ss,ig1,ig2,ifactorx,ifactory)
			propy_d = propy*bwl_mat(:,1)
			propx_d = propx*bwl_mat(1,:)
		endif
		call cuda_setup_many_phasegrate
    else
	    allocate(transf_d(nopiy,nopix,n_slices))
    endif
    if(.not.factorized_propagator) allocate(prop_d(nopiy,nopix))
#endif
    call command_line_title_box('Calculation running')

    do ntilt=1,n_tilts_total
	!For each specimen tilt we have to redo the transmission function
	!and propagators

    if(even_slicing.and.(.not.factorized_propagator)) then
		call make_propagator(nopiy,nopix,prop,prop_distance(1),Kz(ntilt),ss,ig1,ig2&
                       &,claue(:,ntilt),ifactorx,ifactory,exponentiate = .true.)
		prop = prop*bwl_mat
	elseif(.not.factorized_propagator) then
		call make_propagator(nopiy,nopix,prop,1.0_fp_kind,Kz(ntilt),ss,ig1,ig2,claue(:,ntilt),ifactorx,ifactory,exponentiate = .false.)
	endif
	if(.not.on_the_fly) then
    do i = 1, n_slices
        transf_absorptive(:,:,i) = exp(ci*pi*a0_slice(3,i)/Kz(ntilt)*projected_potential(:,:,i))
        ! Bandwith limit the phase grate, psi is used for temporary storage
        write(*,*) sum(transf_absorptive(:,:,i))/nopiy/nopix
        call inplace_fft(nopiy, nopix, transf_absorptive(:,:,i))
        transf_absorptive(:,:,i) = bwl_mat * transf_absorptive(:,:,i)/nopiy/nopix
        call inplace_ifft(nopiy, nopix,transf_absorptive(:,:,i))
        write(*,*) sum(transf_absorptive(:,:,i))/nopiy/nopix
    enddo
	endif
#ifdef GPU
    if(.not.on_the_fly) transf_d=transf_absorptive
    psi_d = psi_initial
    if(.not.factorized_propagator) prop_d = prop
    do i_cell = 1, maxval(ncells)
		intensity = get_sum(psi_d)

900     format(a1, 1x, 'Cell: ', i5, ' Intensity: ', f12.6)
        do i_slice = 1, n_slices
			k =  i_slice+(i_cell-1)*n_slices
			if( modulo(k,10)==0) write(6,900,advance='no') achar(13), k, intensity
			if(double_channeling) then
                do i_target = 1, natoms_slice_total(i_slice) ! Loop over targets
                ! Calculate inelastic transmission matrix
					call cuda_make_shift_array<<<blocks,threads>>>(shiftarray,Hn0_shifty_coord_d(:,i_target,i_slice),Hn0_shiftx_coord_d(:,i_target,i_slice),nopiy,nopix)

					do k = 1, nstates

					call cuda_multiplication<<<blocks,threads>>>(tmatrix_states_d(:,:,k),shiftarray, q_tmatrix_d,1.0_fp_kind,nopiy,nopix)
                    call cufftExec(plan,q_tmatrix_d,tmatrix_d,CUFFT_INVERSE)
					call cuda_multiplication<<<blocks,threads>>>(psi_d,tmatrix_d,psi_inel_d,sqrt(normalisation),nopiy,nopix)
					starting_slice = i_slice
                    ! Scatter the inelastic wave through the remaining cells
                    do ii = i_cell, n_cells
                        do jj = starting_slice, n_slices;call cuda_multislice_iteration(psi_inel_d, transf_d(:,:,jj),prop_d,prop_distance(jj),even_slicing, normalisation, nopiy, nopix,plan);enddo

						if (any(ii==ncells)) then
							z_indx = minloc(abs(ncells-ii))
							call cufftExec(plan, psi_inel_d, tmatrix_d, CUFFT_FORWARD)

							! Accumulate EFTEM images
							do l = 1, imaging_ndf
								call cuda_image(tmatrix_d,ctf_d(:,:,l),temp_d,normalisation, nopiy, nopix,plan,.false.)

								call cuda_addition<<<blocks,threads>>>(eftem_image_d(:,:,l,z_indx(1)), temp_d, eftem_image_d(:,:,l,z_indx(1)), 1.0_fp_kind, nopiy, nopix)
							enddo
						endif
						starting_slice=1
                    enddo
					enddo
				enddo
			endif  ! End loop over cells,targets and states and end double_channeling section

            if(on_the_fly) then
				if(.not.factorized_propagator) call cuda_on_the_fly_abs_multislice(psi_d,ccd_slice_array(i_slice),tau_slice(:,:,:,i_slice),nat_slice(:,i_slice),prop_distance(i_slice),even_slicing,plan,Vg,Volume_array(i_slice),prop_d = prop_d)
				if(factorized_propagator) call cuda_on_the_fly_abs_multislice(psi_d,ccd_slice_array(i_slice),tau_slice(:,:,:,i_slice),nat_slice(:,i_slice),prop_distance(i_slice),even_slicing,plan,Vg,Volume_array(i_slice),propy_d, propx_d)

            else
				call cuda_multislice_iteration(psi_d, transf_d(:,:,i_slice), prop_d,prop_distance(i_slice),even_slicing, normalisation, nopiy, nopix,plan)
            endif
        enddo ! End loop over slices

		!If this thickness corresponds to any of the output values then output images
		if (any(i_cell==ncells)) then
			psi = psi_d
		!
		!	call fft2(nopiy, nopix, psi, nopiy, psi, nopiy)
		!	cbed = abs(psi)**2
		!	z_indx = minloc(abs(ncells-i_cell))
		!	call output_TEM_result(output_prefix,cbed,'Diffraction_pattern',nopiy,nopix,manyz,.false.,manytilt,z=zarray(z_indx(1))&
		!						&,lengthz=lengthz,lengthdf=lengthdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)&
		!						&,nopiyout=nopiy*2/3,nopixout=nopix*2/3)
		!	if(pw_illum) then
		!	do i=1,imaging_ndf
		!		tem_image = make_image(nopiy,nopix,psi,lens_ctf(:,:,i),.false.)
		!		call output_TEM_result(output_prefix,tem_image,'Image',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(z_indx(1))&
		!						&,lengthz=lengthz,lengthdf=lengthdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(i))
				if(double_channeling) then
				do l = 1, imaging_ndf
				tem_image = eftem_image_d(:,:,l,z_indx(1))
				call output_TEM_result(output_prefix,tile_out_image(tem_image,ifactory,ifactorx),'energy_filtered_image',nopiy,nopix,manyz,imaging_ndf>1,manytilt,z=zarray(z_indx(1))&
								&,lengthz=lengthz,lengthdf=lengthdf,tiltstring = tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2),df = imaging_df(l))
				enddo
				endif
		!	enddo
		!	endif
		endif

#else
	psi = psi_initial
    do i_cell = 1, maxval(ncells)
		intensity = sum(abs(psi)**2)
        if(n_tilts_total<2) write(6,900) i_cell, intensity
        if(n_tilts_total>1) write(6,901) i_cell, ntilt,n_tilts_total, intensity
900     format(1h+,1x, 'Cell: ', i5, ' Intensity: ', f12.6)
901     format(1h+,1x, 'Cell: ', i5, ' tilt:',i5,'/',i5,' Intensity: ', f12.6)

        do i_slice = 1, n_slices
            call multislice_iteration(psi,prop,transf_absorptive(:,:,i_slice),prop_distance(i_slice),even_slicing,nopiy,nopix);
        enddo ! End loop over slices

#endif
		!If this thickness cprop_dorresponds to any of the output values then output images
		if (any(i_cell==ncells)) then

			cbed = abs(fft(nopiy, nopix, psi,norm=.true.))**2
			z_indx = minloc(abs(ncells-i_cell))
			filename = trim(adjustl(output_prefix))
			if (nz>1) filename = trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(z_indx(1))),lengthz)//'_A'
            if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
			call binary_out_unwrap(nopiy, nopix, cbed, trim(adjustl(filename))//'_DiffractionPattern')


			if(pw_illum) then
            call binary_out(nopiy, nopix, abs(psi)**2, trim(adjustl(filename))//'_Exit_surface_intensity')
            call binary_out(nopiy, nopix, atan2(imag(psi),real(psi)), trim(adjustl(filename))//'_Exit_surface_phase')
			do i=1,imaging_ndf
				tem_image = make_image(nopiy,nopix,psi,lens_ctf(:,:,i))
				fnam_df = trim(adjustl(filename))// '_Image'
				if(imaging_ndf>1) fnam_df = trim(adjustl(fnam_df))//'_Defocus_'//zero_padded_int(int(imaging_df(i)),lengthdf)//'_Ang'
				call binary_out(nopiy, nopix, tem_image, fnam_df)
            enddo
            endif
		endif

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
		write(9834, '(a, g9.4, a, /)') 'The multislice calculation took ', delta, 'seconds.'
		close(9834)
	endif


end subroutine absorptive_tem
