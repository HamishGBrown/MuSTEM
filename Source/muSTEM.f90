!--------------------------------------------------------------------------------
!   Program: MU_STEM (GPU VERSION)
!
!   Description:    Calculate (S)TEM images and diffraction patterns using the
!                   multislice approach.
!                   Includes the contributions of TDS using either the Hall & Hirsch
!                   absorptive model or the Quantum Excitation of Phonons model.
!                   Both plane wave and convergent beam illumination may be used.
!                   STEM EDX/EELS images can be calculated within the local approximation.
!
!                   Outputs are big-endian floating point numbers in either
!                   32-bit or 64-bit precision, depending on which precision
!                   is chosen in mod_precision.f90.
!
!   Maintainer:     Hamish Brown
!   Email:          hamish.brown@monash.edu
!   Date:           August 2017
!   Requirements:   PGI Fortran
!
!   version:        5.1
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

    program MU_STEM

        use m_user_input
        use global_variables!, only: high_accuracy, nt, atf, nat, atomf, volts, ss, qep, adf, constants, nopiy, nopix,output_thermal,ionic
        use m_lens
#ifdef GPU
        use cuda_setup, only: setup_GPU
        use cuda_array_library, only: set_blocks
#endif
        use output
        use m_tilt, only: prompt_tilt
        use m_absorption
        use m_potential
        use m_multislice
        use m_string
        use m_tilt
        use m_string
        use m_electron
        use fftw3

        implicit none

        integer :: i_illum, i_tds_model, i_cb_calc_type,ifile,nfiles,i_arg,z_indx(1),nx,ny,iz,i_df,i_qep,i,idet,idum,ii,j,k,length&
                    &,length_df,nat_,ntilt,idf

        logical :: nopause = .false.,there,ionization,stem,pacbed,diffraction,elfourd,docbed
        character(512)::command_argument
        character(120)::fnam,filename,addendum,fnam_temp,fnam_df
        real(fp_kind)::const,intensity
        real(8)::thmin,thmax
        real(fp_kind),allocatable :: pacbed_elastic(:,:,:),pacbed_pattern(:,:,:),projected_potential(:,:,:,:),masks(:,:,:)&
                                    &,ion_image(:,:,:),probe_intensity(:,:,:),cbed(:,:,:),stem_ion_image(:,:,:,:,:),&
                                    &stem_image(:,:,:,:,:),stem_elastic_image(:,:,:,:,:),eels_correction_image(:,:,:,:),&
                                    &stem_inelastic_image(:,:,:,:,:),total_intensity(:,:,:),tem_image(:,:,:,:),image(:,:),&
                                    &temp(:,:),adf_potential(:,:,:,:),ionization_potential(:,:,:,:),adf_image(:,:,:)
        complex(fp_kind),allocatable ::lens_ctf(:,:,:),prop(:,:,:),transmission(:,:,:,:),psi(:,:),psi_elastic(:,:,:),psi_out(:,:)&
                                       ,fz_adf(:,:,:,:)
        type(t_slice)::slice
        type(C_PTR) :: forward_plan,inverse_plan


        write(6,109)
        109     format(&
       &1x,'|----------------------------------------------------------------------------|',/,&
       &1x,'|              Melbourne University (scanning) transmission electron         |',/,&
       &1x,'|                            microscopy computing suite                      |',/,&
       &1x,'|      __       __  __    __   ______  ________  ________  __       __       |',/,&
       &1x,'|     |  \     /  \|  \  |  \ /      \|        \|        \|  \     /  \      |',/,&
       &1x,'|     | $$\   /  $$| $$  | $$|  $$$$$$\\$$$$$$$$| $$$$$$$$| $$\   /  $$      |',/,&
       &1x,'|     | $$$\ /  $$$| $$  | $$| $$___\$$  | $$   | $$__    | $$$\ /  $$$      |',/,&
       &1x,'|     | $$$$\  $$$$| $$  | $$ \$$    \   | $$   | $$  \   | $$$$\  $$$$      |',/,&
       &1x,'|     | $$\$$ $$ $$| $$  | $$ _\$$$$$$\  | $$   | $$$$$   | $$\$$ $$ $$      |',/,&
       &1x,'|     | $$ \$$$| $$| $$__/ $$|  \__| $$  | $$   | $$_____ | $$ \$$$| $$      |',/,&
       &1x,'|     | $$  \$ | $$ \$$    $$ \$$    $$  | $$   | $$     \| $$  \$ | $$      |',/,&
       &1x,'|      \$$      \$$  \$$$$$$   \$$$$$$    \$$    \$$$$$$$$ \$$      \$$      |',/,&
       &1x,'|                                                                            |',/,&
       &1x,"|       Copyright (C) 2017 L.J. Allen, H.G. Brown, A.J. D'Alfonso,           |",/,&
       &1x,'|              S.D. Findlay, B.D. Forbes                                     |',/,&
       &1x,'|       email: hamish.brown@monash.edu                                       |',/,&
       &1x,'|       This program comes with ABSOLUTELY NO WARRANTY;                      |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|       This program is licensed to you under the terms of the GNU           |',/,&
       &1x,'|       General Public License Version 3 as published by the Free            |',/,&
       &1x,'|       Software Foundation.                                                 |',/,&
       &1x,'|                                                                            |',/,&
#ifdef GPU
       &1x,'|       GPU Version 5.3                                                      |',/,&
#else
       &1x,'|       CPU only Version 5.3                                                 |',/,&
#endif
       &1x,'|                                                                            |',/,&
       &1x,'|       Note: pass the argument "nopause" (without quotation marks)          |',/,&
       &1x,'|             e.g. muSTEM.exe nopause                                        |',/,&
       &1x,'|             to avoid pauses.                                               |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|----------------------------------------------------------------------------|',/)

        ! Process command line arguments
        ! Process command line arguments
        do i_arg = 1, command_argument_count()
            call get_command_argument(i_arg, command_argument)
            select case (trim(adjustl(command_argument)))
            case ('nopause')
                nopause = .true.
            case ('timing')
                timing = .true.
            case ('ionic')
                ionic = .true.
            end select
        enddo

        if (.not. nopause) then
            write(*,*) ' Press enter to continue.'
            read(*,*)
        endif



        ! Set up user input routines, nfiles is the number
        ! of user input files to play if "play all" is inputted
        nfiles = init_input()

        do ifile=1,nfiles
            !If play or play all open relevant user input file.
            if(input_file_number.ne.5) then
                fnam = get_driver_file(ifile)
                inquire(file=fnam,exist = there)
                if (there) then
                    open(unit=in_file_number, file=fnam, status='old')
                else
                    write(*,*) "Couldn't find user input file: ",trim(adjustl(fnam))
                    cycle
                endif
            endif
        ! Set up CPU multithreading
        call setup_threading
#ifdef GPU
        ! Set up GPU
        call setup_GPU
#else
        !open(6,carriagecontrol ='fortran')
#endif

        call command_line_title_box('Dataset output')
        ! Ask for the prefix for dataset output
        call setup_output_prefix

        call command_line_title_box('Specimen details')
        ! Read in the xtl file
        call set_xtl_global_params
        call validate_xtl(deg)

        ! Calculate the mean inner potential
        call set_volts(nt, atf, nat, atomf, volts, ss)

        ! Set electron quantities
        call constants

        ! Ask for slicing scheme
        call setup_slicing_depths(slice)


        ! Ask for thickness
        call setup_specimen_thickness

        ! Set up the potential calculation method
        call prompt_high_accuracy

        ! Set the unit cell tiling and grid size
        call set_tiling_grid

        i_illum = -1
        do while(i_illum<1.or.i_illum>2)
            call command_line_title_box('Illumination')
            write(*,*) '<1> Plane wave      (including HRTEM images and diffraction patterns)'
            write(*,*) '<2> Convergent-beam (including STEM images and CBED patterns)'

            call get_input('Calculation type', i_illum)
            write(*,*)
        enddo

        ! Set illumination
        pw_illum = (i_illum == 1)



        ! For convergent-beam, set up the probe forming lens
        if (.not.pw_illum) call setup_lens_parameters('Probe',probe_aberrations,probe_cutoff)


        i_tds_model=-1
        do while(i_tds_model<1.or.i_tds_model>2)
            call command_line_title_box('Thermal scattering model')
            write(*,*) '<1> Quantum Excitation of Phonons (QEP) model'
            write(*,*) '    (Accurately accounts for inelastic scattering'
            write(*,*) '     due to phonon excitation)'
            write(*,*) '<2> Absorptive model'
            write(*,*) '     (Calculates faster than QEP but only approximates'
            write(*,*) '      inelastic scattering due to phonon excitation)'

            call get_input('<1> QEP <2> ABS', i_tds_model)
            write(*,*)
        enddo

        qep = (i_tds_model == 1)

        ! Prompt for including absorptive potential
        include_absorption =.false.
        if (.not. qep) call prompt_include_absorption
        if (.not. qep) nran=1
        if(((.not.qep).and.include_absorption).or.qep) then
            write(*,*) ' Options for output of inelastically and elastically scattered components'
            write(*,*) ' <1> Only output total signal (ie that measured in experiment)'
            write(*,*) ' <2> Seperately output elastic, inelastic and total signal'
            if(.not.qep) write(*,*) 'Note: option <2> only applies to STEM imaging'
        write(*,*)
        call get_input('Elastic and inelastic scattering output choice', i_tds_model)
        write(*,*)
        output_thermal = i_tds_model==2
        endif


        ! Prompt user for a tilt for either kind of illumination
        call prompt_tilt
#ifdef GPU
        ! Setup the CUDA thread hierachy for nopiy, nopix arrays
        call set_blocks(nopiy, nopix)
#endif
        ! Calculate the slices of the supercell
        call initialize_slicing(slice,nat,ifactory,ifactorx,nt,a0,deg,tau,nm)

        ! Calculate the bandwidth limiting matrix
        call make_bwl_mat

        ! Ask for QEP parameters
        call setup_qep_parameters(qep,n_qep_grates,n_qep_passes,nran,quick_shift,ifactory,ifactorx)

        ! Save/load transmission functions
        call prompt_save_load_grates(slice)

        ionization= .false.;STEM = .false.;pacbed=.false.;docbed=.false.


        if (pw_illum) then
            ! Set up the imaging lens
            nysample=1;nxsample=1;imaging=.true.;probe_ndf =1
             call setup_lens_parameters('Image',imaging_aberrations,imaging_cutoff)
             pacbed=.false.;stem=.false.
        else
            imaging_ndf = 1
            call command_line_title_box('Calculation type')
            i_cb_calc_type = -1
            do while(i_cb_calc_type>2 .or. i_cb_calc_type<1)

            write(*,*) 'Choose a calculation type:'
115         write(*,*) '<1> CBED pattern'
            write(*,*) '<2> STEM (BF/ABF/ADF/EELS/EDX/PACBED/4D-STEM)'

            call get_input('<1> CBED <2> STEM/PACBED', i_cb_calc_type)
            write(*,*)

            if(i_cb_calc_type==1) then
                ! CBED pattern
                nysample=1
                nxsample=1
                allocate(probe_positions(3,1,1))
                probe_positions(:,1,1) = place_probe()
                docbed = .true.
            else if(i_cb_calc_type==2) then
                ! STEM images
                call STEM_options(STEM,ionization,PACBED)
                diffraction= .false.
                if(pacbed) call fourD_STEM_options(diffraction,nopiyout,nopixout,nopiy,nopix)
                call setup_probe_scan(PACBED.and.(.not.(ionization.or.STEM)))
                call prompt_output_probe_intensity(slice)
                if(STEM) call setup_integration_measurements
                adf = STEM.and.(.not.qep)

                ! Precalculate the scattering factors on a grid
                call precalculate_scattering_factors()
                if(ionization) call setup_ionization(EDX)

            endif
        enddo
        endif

        call command_line_title_box('Pre-calculation setup')
        !Allocate wave function and array in which to store temporary real space arrays
        allocate(psi(nopiy,nopix),temp(nopiy,nopix))

        !Seed random number generator
        idum = seed_rng()

        !Allocate PACBED variables
        if(pacbed) then
            allocate(pacbed_pattern(nopiy,nopix,nz))
            pacbed_pattern=0
            if(qep.and.output_thermal) then
                allocate(pacbed_elastic(nopiy,nopix,nz))
                pacbed_elastic = 0
            endif
            elfourd = fourdSTEM.and.output_thermal
        endif

        write(*,*) 'Planning Fast Fourier Transforms...'
        forward_plan = fft_plan(nopiy,nopix,psi,-1,0)
        inverse_plan = fft_plan(nopiy,nopix,psi,+1,0)


        ! Precalculate the scattering factors on a grid
        write(*,*) 'Precalculating scattering factors...'
        call precalculate_scattering_factors()
        write(*,*) 'Making transmission functions...'
        allocate(projected_potential(nopiy,nopix,n_qep_grates,slice%n_slices)&
                       ,transmission(nopiy,nopix,n_qep_grates,slice%n_slices)&
                       ,        prop(nopiy,nopix,             slice%n_slices)&
                       ,        cbed(nopiy,nopix,nz)&
                    ,total_intensity(nopiy,nopix,nz)&
                    ,    psi_elastic(nopiy,nopix,nz)&
                    )

        if(qep) then
            projected_potential(:,:,:,:) = make_qep_grates(slice,nopiy,nopix,n_qep_grates,idum)
        else

            projected_potential(:,:,1,:) = make_absorptive_grates(slice,nopiy,nopix)
        endif

        !Make detectors for STEM if required
        if(stem) then

            allocate(masks(nopiy,nopix,ndet),stem_image(nysample,nxsample,probe_ndf,ndet,nz))
            if(stem.and.((qep.and.output_thermal).or.((.not.qep).and.complex_absorption))) then
                allocate(stem_elastic_image(nysample,nxsample,probe_ndf,ndet,nz))
            endif
            write(*,*) size(stem_elastic_image)
            write(*,*) 'Making detectors...'
        do i=1,ndet/nseg
        do j=1,nseg
            if(nseg>1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i),&
                                                                             2*pi*j/nseg-seg_det_offset,2*pi/nseg)
            if(nseg==1) masks(:,:,(i-1)*nseg+j) = make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))
        enddo
        enddo
        endif

        !Make ionization potentials if required
        if(ionization) then
            write(*,*) 'Setting up ionization potentials...'
            allocate(ionization_potential(nopiy,nopix,num_ionizations,slice%n_slices)&
                   &,ion_image(nopiy,nopix,num_ionizations),stem_ion_image(nysample,nxsample,probe_ndf,nz,num_ionizations))
            stem_ion_image = 0

            do i=1,num_ionizations;do j=1,slice%n_slices
                nat_ = slice%nat(atm_indices(i),j)
                ionization_potential(:,:,i,j) = real(potential_from_scattering_factors(ionization_mu(:,:,i),&
                                                         slice%tau(:,atm_indices(i),1:nat_,j),nat_,nopiy,nopix,&
                                                         high_accuracy)/slice%ss(7,j))
            enddo;enddo
        endif

        !Make local phonon scattering potentials if required
        write(*,*) all([stem,.not.qep,include_absorption]),stem,.not.qep,include_absorption
        if(all([stem,.not.qep,include_absorption])) then
            allocate(adf_potential(nopiy,nopix,slice%n_slices,ndet),fz_adf(nopiy,nopix,nt,ndet),adf_image(nopiy,nopix,ndet))
            adf_potential = 0
            write(*,*) 'Calculating phonon inelastic scattering potentials...'

            do k=1,ndet/nseg
                    thmin =  atan(inner((k-1)/nseg+1)/ak)
                    thmax =  atan(outer((k-1)/nseg+1)/ak)
                    !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
                    !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
                    fz_adf(:,:,:,(k-1)*nseg+1:k*nseg) = spread(absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,&
                                                            nt,a0,ss,atf,nat, ak, relm, orthog,thmin,thmax),dim=4,ncopies=nseg)/nseg
            enddo
            do j=1,slice%n_slices;do i=1,nt
                nat_ = slice%nat(i,j)
                write(*,*) size(adf_potential),size(fz_adf)
                do k=1,ndet
                    adf_potential(:,:,j,k)= adf_potential(:,:,j,k) + real(potential_from_scattering_factors(fz_adf(:,:,i,k)&
                                           ,slice%tau(:,i,:nat_,j),nat_,nopiy,nopix,high_accuracy)/slice%ss(7,j)*ss(7)*4*pi)
                enddo
            enddo;enddo
            do k=1,nt

                call binary_out(nopiy,nopix,adf_potential(:,:,1,k),'scattering_factors_'//to_string(k))

            enddo
            write(*,*) 'made'
            deallocate(fz_adf)
        endif

        if(imaging) then
            allocate(lens_ctf(nopiy,nopix,imaging_ndf),tem_image(nopiy,nopix,nz,imaging_ndf))
            write(*,*) 'Setting up imaging lenses....'
            do i=1,imaging_ndf
                if(pw_illum) lens_ctf(:,:,i) =  make_ctf([0.0_fp_kind,0.0_fp_kind,0.0_fp_kind],imaging_df(i),imaging_cutoff,&
                                                                                imaging_aberrations,imaging_apodisation)
            enddo
        endif

        length = ceiling(log10(maxval(zarray)))
        call command_line_title_box('Calculation')



        !Loop over specimen tilts
        do ntilt=1,n_tilts_total

        !Make transmission function and propagator
        do i = 1, slice%n_slices
            call make_propagator(nopiy,nopix,prop(:,:,i),slice%a0(3,i),Kz(1),ss,claue(:,1),ifactorx,ifactory)
            prop(:,:,i) = prop(:,:,i) * bwl_mat
            !stop
            transmission(:,:,:,i) = exp(ci*pi*slice%a0(3,i)/Kz(ntilt)*projected_potential(:,:,:,i))/nopiy/nopix

            do j=1,n_qep_grates
            call inplace_fft(nopiy,nopix,transmission(:,:,j,i),norm=.true.)
            transmission(:,:,j,i)= transmission(:,:,j,i)*bwl_mat
            if(qep_mode.ne.3) call inplace_ifft(nopiy,nopix,transmission(:,:,j,i),norm=.true.)
            enddo
        enddo
        !Loop over...
        do i_df = 1, probe_ndf !Probe defocii
        do ny = 1, nysample !Probe scan positions in y and...
        do nx = 1, nxsample !... x directions
        psi_elastic=0;cbed=0;total_intensity=0;
        do i_qep=1, n_qep_passes !QEP multislice passes (1 in case of absorptive calculation)
        !Set up illumination
        if(pw_illum) then
            psi = 1.0_fp_kind/sqrt(float(nopiy*nopix))
        else
            psi = make_ctf(probe_positions(:,ny,nx),probe_df(i_df),probe_cutoff,probe_aberrations,probe_apodisation)
            call inplace_ifft(nopiy, nopix, psi)
            psi = psi/sqrt(sum(abs(psi)**2))
        endif
        call tilt_wave_function(psi)

        !Initialize variables
        if (ionization) ion_image = 0
        if(all([stem,complex_absorption,include_absorption])) adf_image =0

        do i = 1,maxval(ncells)
            call multislice_progress([i_df,ny,nx,i_qep,i],[probe_ndf,nysample,nxsample,n_qep_passes,maxval(ncells)],&
                                &['df   ','ny   ','nx   ','nqep ','slice'],sum(abs(psi)**2))
                do j = 1, slice%n_slices
                    ! Accumulate ionization cross section
                    if (ionization.or.all([stem,complex_absorption,include_absorption])) temp = abs(psi)**2
                    if(ionization) then
                        do ii=1,num_ionizations
                            ion_image(:,:,ii) = temp* ionization_potential(:,:,ii,j) * slice%a0(3,j)+ion_image(:,:,ii)

                        enddo
                    endif

                    if(all([stem,complex_absorption,include_absorption])) then
                        do k=1,ndet;adf_image(:,:,k) = abs(psi)**2*adf_potential(:,:,j,k)* slice%a0(3,j) + adf_image(:,:,k); enddo
                    endif
                !Apply multislice operation, the qep multislice function is called, if an absorptive calculation was
                !requested then n_qep_grates is set to zero and qep_mode is set to 4 (no shifting)
                  call qep_multislice_iteration(psi,prop(:,:,j),transmission(:,:,:,j),nopiy,nopix,ifactory,ifactorx,idum,&
                                              n_qep_grates,qep_mode,shift_arrayy,shift_arrayx,forward_plan,inverse_plan)

                    if (output_probe_intensity) then
                        k = (i-1)*slice%n_slices+j
                        if (output_cell_list(k)) then
                            probe_intensity(:,:,cell_map(k)) = probe_intensity(:,:,cell_map(k)) + abs(psi)**2
                        endif
                    endif
                enddo ! End loop over slices

                !If this thickness corresponds to any of the output values then accumulate diffraction pattern
                if (any(i==ncells)) then
                    z_indx = minloc(abs(ncells-i))

                    ! Accumulate elastic wave function - this will be Fourier transformed later
                    psi_elastic(:,:,z_indx(1)) = psi_elastic(:,:,z_indx(1)) + psi

                    !Accumulate exit surface intensity
                    total_intensity(:,:,z_indx(1)) = total_intensity(:,:,z_indx(1)) + abs(psi)**2

                    !Transform into diffraction space
                    psi_out = fft(nopiy,nopix,psi,norm=.true.)

                    ! Accumulate diffaction pattern
                    temp = abs(psi_out)**2
                    cbed(:,:,z_indx(1)) = cbed(:,:,z_indx(1)) + temp

                    if(imaging) then
                    do idf=1,imaging_ndf
                        temp = abs(ifft(nopiy,nopix,lens_ctf(:,:,idf)*psi_out,norm=.true.))
                        tem_image(:,:,z_indx(1),idf) = tem_image(:,:,z_indx(1),idf)+temp**2
                    enddo
                    endif

                    if(ionization) stem_ion_image(ny,nx,i_df,z_indx(1),:) = stem_ion_image(ny,nx,i_df,z_indx(1),:)&
                                                                          + sum(sum(ion_image,dim=2),dim=1)
                endif

            enddo ! End loop over cells
        enddo !End loop over QEP passes
        if(n_qep_passes>1) psi_elastic = psi_elastic/n_qep_passes
        if(n_qep_passes>1) cbed = cbed/n_qep_passes
        if(n_qep_passes>1) total_intensity = total_intensity/n_qep_passes
        if(n_qep_passes>1 .and. imaging) tem_image = tem_image/n_qep_passes

        do iz=1,nz


            ! Integrate the elastic diffraction pattern

            if(qep) psi_out = fft(nopiy,nopix,psi_elastic(:,:,iz))
            if(stem) then
                if(qep) then
                    stem_image(ny,nx,i_df,1:ndet,iz) = sum(sum(spread(cbed(:,:,iz),dim=3,ncopies=ndet)*masks(:,:,:),dim=1),dim=1)
                    stem_elastic_image(ny,nx,i_df,1:ndet,iz)= sum(sum(spread(abs(psi_out)**2,dim=3,ncopies =ndet)&
                                                                                  *masks(:,:,1:ndet),dim=1),dim=1)
                else
                    stem_image(ny,nx,i_df,1:ndet,iz) = sum(sum(spread(cbed(:,:,iz)&
                                                                 &,dim=3,ncopies=ndet)*masks(:,:,:),dim=1),dim=1)
                    if(include_absorption) then

                    stem_elastic_image(ny,nx,i_df,1:ndet,iz) = stem_image(ny,nx,i_df,1:ndet,iz)
                    stem_image(ny,nx,i_df,1:ndet,iz) = stem_image(ny,nx,i_df,1:ndet,iz) &
                                                                     &+ sum(sum(adf_image(:,:,1:ndet),1),1)
                    endif
                endif
            endif
            if(ionization.and.(.not.EDX)) eels_correction_image(ny,nx,i_df,z_indx(1)) = sum(cbed(:,:,iz)*eels_correction_detector)

            addendum = ''
            if (probe_ndf>1) addendum = triml(addendum)//defocus_string(probe_df(i_df),length_df)
            if (nz>1) addendum = triml(addendum)//'_z='//to_string(int(zarray(iz)))//'_A'
            if((nysample>1).or.(nxsample>1)) addendum = triml(addendum)//'_pp_'//&
                                          &to_string(nx)//'_'//to_string(ny)
            !Output 4D STEM diffraction pattern
            if(diffraction) call binary_out_unwrap(nopiy, nopix, cbed(:,:,iz), triml(output_prefix)//&
                                          &triml(addendum) //'_Diffraction_pattern',write_to_screen=.false. &
                                          &,nopiyout=nopiyout,nopixout=nopixout)

            if(pacbed) then

                if(i_df==1) pacbed_pattern(:,:,iz) = pacbed_pattern(:,:,iz) + cbed(:,:,iz)
                if(qep.and.output_thermal) pacbed_elastic(:,:,iz) =pacbed_elastic(:,:,iz) + abs(psi_out)**2

            endif
        enddo



        if (output_probe_intensity) call probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,&
                                                                 probe_ndf,nysample,nxsample)
    enddo ! End loop over x probe positions
    enddo ! End loop over y probe positions
    enddo ! End loop over defocus series

    if(stem) then
        fnam = trim(adjustl(output_prefix))
        if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
        if(output_thermal) then
            allocate(stem_inelastic_image(nysample,nxsample,probe_ndf,ndet,nz))
            stem_inelastic_image = stem_image - stem_elastic_image
        endif
        do idet = 1, ndet

            if(output_thermal) then
                fnam_temp = triml(fnam) // '_DiffPlaneElastic_Detector'//zero_padded_int(idet,2)

                call output_stem_image(stem_elastic_image(:,:,:,idet,:),fnam_temp,probe_df)


                fnam_temp = triml(fnam) // '_DiffPlaneTDS_Detector_'//zero_padded_int(idet, 2)
                call output_stem_image(stem_inelastic_image(:,:,:,idet,:),fnam_temp ,probe_df)

                fnam_temp = triml(fnam) // '_DiffPlaneTotal_Detector'
            else
                fnam_temp = triml(fnam) // '_DiffPlane_Detector'
            endif

            call output_stem_image(stem_image(:,:,:,idet,:),triml(fnam_temp)//zero_padded_int(idet,2),probe_df)
        enddo
        if(output_thermal) deallocate(stem_inelastic_image)
        !ionization
        if(ionization) then
            do ii=1,num_ionizations
            if(EDX) then
                filename = trim(adjustl(fnam)) // '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'&
                                                     //trim(adjustl(Ion_description(ii)))//'_shell_EDX'
                call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            else
                filename = trim(adjustl(fnam))// '_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'&
                                                   //trim(adjustl(Ion_description(ii)))//'_orbital_EELS'
                call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
                stem_ion_image(:,:,:,:,ii) = stem_ion_image(:,:,:,:,ii)*eels_correction_image

                filename =  trim(adjustl(fnam)) //'_'//trim(adjustl(substance_atom_types(atm_indices(ii))))//'_'&
                                                 //trim(adjustl(Ion_description(ii)))//'_orbital_EELS_Corrected'
                call output_stem_image(stem_ion_image(:,:,:,:,ii), filename,probe_df)
            endif
            enddo
            if(.not.EDX) then
                filename = trim(adjustl(fnam)) // '_EELS_CorrectionMap'
                call output_stem_image(eels_correction_image, filename,probe_df)
            endif

        endif
    endif

    if(pacbed) then
            const = float(nysample*nxsample)
            do i=1,nz

            filename = trim(adjustl(output_prefix))
            if (n_tilts_total>1) fnam = trim(adjustl(output_prefix))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)
            if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A_'
            call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const,trim(adjustl(filename))//'_PACBED_Pattern')

            if(qep.and.output_thermal) then
                call binary_out_unwrap(nopiy,nopix,PACBED_elastic(:,:,i)/const,trim(adjustl(filename))//'_elastic_PACBED_Pattern')
                call binary_out_unwrap(nopiy,nopix,pacbed_pattern(:,:,i)/const-PACBED_elastic(:,:,i)/const,trim(adjustl(filename))&
                                                                                                      //'_thermal_PACBED_Pattern')
            endif
        enddo
    endif
    if(pw_illum.or.docbed) then
        do i=1,nz
            filename = trim(adjustl(output_prefix))
            if(nz>1) filename=trim(adjustl(filename))//'_z='//zero_padded_int(int(zarray(i)),length)//'_A'
            if (n_tilts_total>1) filename = trim(adjustl(filename))//tilt_description(claue(:,ntilt),ak1,ss,ig1,ig2)

            if((.not.output_thermal).or.(.not.qep)) then
                call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffPlane')
                if(pw_illum) call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity')
            else
                call binary_out_unwrap(nopiy, nopix, cbed(:,:,i), trim(filename)//'_DiffPlaneTotal')
                psi = fft(nopiy, nopix, psi_elastic(:,:,i),norm=.true.)
                image = abs(psi)**2
                call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffPlaneElastic')

                image = cbed(:,:,i) - image
                call binary_out_unwrap(nopiy, nopix, image, trim(filename)//'_DiffPlaneTDS')

                if(pw_illum) then
                    call binary_out(nopiy, nopix, abs(psi_elastic(:,:,i))**2, &
                              &trim(filename)//'_ExitSurface_Intensity_Elastic')
                    call binary_out(nopiy, nopix, atan2(imag(psi_elastic(:,:,i)), real(psi_elastic(:,:,i))), trim(filename)&
                                                                                             //'_ExitSurface_Phase_Elastic')
                    call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity_Total')

                    total_intensity(:,:,i) = total_intensity(:,:,i) - abs(psi_elastic(:,:,i))**2
                    call binary_out(nopiy, nopix, total_intensity(:,:,i), trim(filename)//'_ExitSurface_Intensity_TDS')
                endif
            endif
            if(pw_illum) then
            do j=1,imaging_ndf
                ! Elastic image
                psi = fft(nopiy, nopix, psi_elastic(:,:,i),norm=.true.)* lens_ctf(:,:,j)
                call inplace_ifft (nopiy, nopix, psi,norm=.true.)
                image = abs(psi)**2

                if(imaging_ndf>1) then
                    fnam_df = trim(adjustl(filename))//'_Defocus_'//zero_padded_int(int(imaging_df(j)),length_df)//'_Ang'
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
            enddo
            endif
        enddo
    endif

    enddo !End loop over tilts

    deallocate(projected_potential,transmission,prop,cbed,total_intensity,psi_elastic)
    if(stem) deallocate(masks,stem_image)
    if(stem.and.all([stem,complex_absorption,include_absorption])) deallocate(adf_potential)
    if(stem.and.output_thermal) deallocate(stem_elastic_image)
    if(ionization) deallocate(stem_ion_image,ionization_mu,ionization_potential)
    if(ionization.and.(.not.EDX)) deallocate(eels_correction_image)
    if(pacbed) deallocate(pacbed_pattern)
    if(pacbed.and.(qep.and.output_thermal)) deallocate(pacbed_elastic)
    if(imaging) deallocate(tem_image,lens_ctf)
    enddo !End loop over instruction files
    end program Mu_STEM



    subroutine setup_threading()

        use m_string, only: to_string,command_line_title_box

        implicit none

        integer*4 :: num_threads
        integer*4 :: omp_get_max_threads, omp_get_num_procs

        !num_threads = omp_get_num_procs()

        !call omp_set_num_threads(num_threads)

        call command_line_title_box('CPU multithreading')
        write(*,*) 'The number of threads being used on the CPU is: ' // to_string(num_threads)
        write(*,*)

    end subroutine
      subroutine reset_allocatable_variables()

        !The downside of using global variables... :(
        use global_variables
        use m_lens
        use m_potential
        use m_multislice
#ifdef GPU
        use cuda_potential
#endif
        if(allocated(Kz)                      ) deallocate(Kz)
        if(allocated(claue)                   ) deallocate(claue)
        if(allocated(nat)                     ) deallocate(nat)
        if(allocated(tau)                     ) deallocate(tau)
        if(allocated(atf)                     ) deallocate(atf)
        if(allocated(atomf)                   ) deallocate(atomf)
        if(allocated(fx)                      ) deallocate(fx)
        if(allocated(dz)                      ) deallocate(dz)
        if(allocated(zarray)                  ) deallocate(zarray)
        if(allocated(ncells)                  ) deallocate(ncells)
        if(allocated(fz    )                  ) deallocate(fz)
        if(allocated(fz_DWF)                  ) deallocate(fz_DWF)
        if(allocated(sinc)                    ) deallocate(sinc)
        if(allocated(inverse_sinc)            ) deallocate(inverse_sinc)
        if(allocated(substance_atom_types    )) deallocate(substance_atom_types)
        if(allocated(outer                   )) deallocate(outer)
        if(allocated(inner                   )) deallocate(inner)
        if(allocated(probe_df                )) deallocate(probe_df          )
        if(allocated(imaging_df              )) deallocate(imaging_df          )
        if(allocated(atm_indices             )) deallocate(atm_indices          )
#ifdef GPU
        if(allocated(ccd_slice_array         )) deallocate(ccd_slice_array)
        if(allocated(Volume_array            )) deallocate(Volume_array)
#endif
        end subroutine
