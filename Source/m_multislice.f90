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

module m_multislice

    use m_precision, only: fp_kind

    implicit none

    logical :: save_grates = .false.
    logical :: load_grates = .false.
    character(1024) :: grates_filename

    integer*4::qep_mode
    logical :: output_probe_intensity = .false.,additional_transmission_function = .false.,pure_phase
    logical,allocatable :: output_cell_list(:)
    real(fp_kind),allocatable :: output_thickness_list(:)
    integer,allocatable :: cell_map(:)

    character*200,allocatable::amplitude_fnam(:),phase_fnam(:)


    complex(fp_kind),allocatable,dimension(:,:) :: shift_arrayx, shift_arrayy

    integer(4) :: n_slices                             !number of subslices in the unit cell
    integer(4) :: maxnat_slice                            !maximum number of atoms for a particular type in the supercell
    real(fp_kind), allocatable :: a0_slice(:,:)
    integer(4),    allocatable :: nat_slice(:,:)             !number of atoms for each type in the slice (supercell)
    integer(4),    allocatable :: nat_slice_unitcell(:,:)    !number of atoms for each type in the slice (unit cell)
    real(fp_kind), allocatable :: tau_slice(:,:,:,:)      !the atom positions for each potential subslice (supercell)
    real(fp_kind), allocatable :: tau_slice_unitcell(:,:,:,:) !the atom positions for each potential subslice (unit cell)
    real(fp_kind), allocatable :: prop_distance(:)        !Unit cell subslice propagation distance
    real(fp_kind), allocatable :: depths(:)               !Unit cell subslice depths (i.e. same as above)
    real(fp_kind), allocatable :: ss_slice(:,:)

    interface load_save_add_grates
        module procedure load_save_add_grates_qep,load_save_add_grates_abs
    end interface

    interface multislice_iteration
        module procedure multislice_iteration_new,multislice_iteration_old
    end interface

    contains

    !This subroutine samples from the available phase grates and then performs one iteration of the multislice algorithm (called in CPU versions only)

    subroutine qep_multislice_iteration(psi,propagator,transmission,nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,mode&
                                                                    ,shift_arrayy,shift_arrayx,forward_plan,inverse_plan)
        use m_numerical_tools,only:ran1
        use, intrinsic :: iso_c_binding
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix,n_qep_grates)
        complex(fp_kind),intent(in)::shift_arrayy(nopiy),shift_arrayx(nopix)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,mode,ifactory,ifactorx
        type(c_ptr),intent(in),optional::forward_plan,inverse_plan
        integer*4,intent(inout)::idum

        integer*4::shifty,shiftx,nran
        complex(fp_kind)::trans(nopiy,nopix)

        ! Phase grate
        nran = floor(n_qep_grates*ran1(idum)) + 1

        if(mode == 1) then !On the fly calculation
            !call make_qep_potential(trans, tau_slice, nat_slice, ss_slice(7,j))
            !psi_out = psi*trans
        elseif(mode == 2) then !Quick shift
            shiftx = floor(ifactorx*ran1(idum)) * nopix/ifactorx
            shifty = floor(ifactory*ran1(idum)) * nopiy/ifactory
            trans = cshift(cshift(transmission(:,:,nran),shifty,dim=1),shiftx,dim=2)
        elseif(mode == 3) then      !Phase ramp shift
            shiftx = floor(ifactorx*ran1(idum)) + 1
            shifty = floor(ifactory*ran1(idum)) + 1
            call phase_shift_array(transmission(:,:,nran),trans,shift_arrayy,shift_arrayx)
        else
            trans = transmission(:,:,nran)
        endif
        if(all([present(forward_plan),present(inverse_plan)])) then
            call multislice_iteration(psi,propagator,trans,nopiy,nopix,forward_plan,inverse_plan)
        else
            call multislice_iteration(psi,propagator,trans,nopiy,nopix)
        end if

    end subroutine

    !This subroutine performs one iteration of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration_old(psi,propagator,transmission,nopiy,nopix,forward_plan,inverse_plan)
	    use FFTW3
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix
        type(C_PTR),intent(in),optional :: forward_plan,inverse_plan


        ! Transmit through slice potential (assumes transmission function divided by 1/nopiy/nopix)
        psi = psi*transmission
        ! Propagate to next slice
        if(present(forward_plan)) then
            call inplace_fft(nopiy,nopix,psi,forward_plan)
        else
            call inplace_fft(nopiy,nopix,psi)
        endif
        psi = psi*propagator

        if(present(inverse_plan)) then
            call inplace_ifft(nopiy,nopix,psi,inverse_plan)
        else
            call inplace_ifft(nopiy,nopix,psi)
        endif
    end subroutine

    !This subroutine performs many iterations of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration_many_slices(psi,propagator,transmission,nopiy,nopix,ncells,nslices,forward_plan,inverse_plan)
        use FFTW3
        use output
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in),dimension(nopiy,nopix,nslices)::propagator,transmission
        integer*4,intent(in)::nopiy,nopix,ncells,nslices
        integer*4::icell,islice
        type(C_PTR),intent(in),optional :: forward_plan,inverse_plan

        do icell=1,ncells;do islice=1,nslices
        ! Transmit through slice potential (assumes transmission function divided by 1/nopiy/nopix)
        psi = psi*transmission(:,:,islice)
        ! Propagate to next slice
        if(present(forward_plan)) then
            call inplace_fft(nopiy,nopix,psi,forward_plan)
        else
            call inplace_fft(nopiy,nopix,psi)
        endif
        psi = psi*propagator(:,:,islice)

        if(present(inverse_plan)) then
            call inplace_ifft(nopiy,nopix,psi,inverse_plan)
        else
            call inplace_ifft(nopiy,nopix,psi)
        endif
        enddo;enddo
    end subroutine

    !This subroutine performs one iteration of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration_new(psi,propagator,transmission,prop_distance,even_slicing,nopiy,nopix)
	    use FFTW3
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix
        real(fp_kind),intent(in)::prop_distance
        logical,intent(in)::even_slicing

        ! Transmit through slice potential
		psi = psi*transmission

        ! Propagate to next slice
		call fft2(nopiy,nopix,psi,nopiy,psi,nopiy)
		if(even_slicing) then
            psi = psi*propagator
        else
            psi = psi*exp(prop_distance*propagator)
        endif
		call ifft2(nopiy,nopix,psi,nopiy,psi,nopiy)

    end subroutine

    subroutine prompt_output_probe_intensity

        use m_user_input, only: get_input
        use global_variables, only: thickness, nopiy, nopix
        use output, only: output_prefix,split_filepath
        use m_string, only: to_string,command_line_title_box

        implicit none

        integer :: i_output
        real(fp_kind) :: thickness_interval
        character(1024)::dir,fnam
        logical::invalid_thickness

        call command_line_title_box('Output probe intensity')
        write(*,*) 'The probe intensity as a function of thickness'
        write(*,*) 'can be outputted to file at each probe position.'
        write(*,*) 'The user is advised that the outputted dataset'
        write(*,*) 'may be very large.'
        write(*,*)
        write(*,*) '<0> Proceed'
        write(*,*) '<1> Output probe intensity'
        call get_input('<0> Proceed <1> Output probe intensity', i_output)
        write(*,*)

        output_probe_intensity = i_output ==1

        if(output_probe_intensity) then

            invalid_thickness = .true.
            do while(invalid_thickness)
                write(*,*) 'At what thickness interval (in Angstroms)',char(10),' should intensities be outputted?'
                call get_input('At what thickness interval should intensities be outputted?', thickness_interval)
                write(*,*)
                invalid_thickness = thickness_interval.le.0.0_fp_kind .or. thickness_interval.gt.thickness
                if(invalid_thickness) write(*,*) 'ERROR: invalid thickness.'
            enddo

                call split_filepath(output_prefix,dir,fnam)
                call system('mkdir '//trim(adjustl(dir))//'\Probe_intensity')
                call generate_cell_list(thickness_interval)
                call write_thicknesss_to_file

                write(*,*) 'The probe intensities will be written to the files'
                write(*,*)
                write(*,*) '  '//trim(adjustl(dir))//'\Probe_intensity' // trim(adjustl(fnam)) // '_ProbeIntensity*.bin'
                write(*,*)
                if (fp_kind.eq.4) write(*,*) 'as 32-bit big-endian floats.'
                if (fp_kind.eq.8) write(*,*) 'as 64-bit big-endian floats.'

                write(*,*) 'Each file contains a sequence of ' // to_string(nopiy) // 'x' // to_string(nopix) // ' arrays.'
                write(*,*)
        end if

    end subroutine



    subroutine generate_cell_list(thickness_interval)

        use global_variables, only: ncells, a0

        implicit none

        real(fp_kind) :: thickness_interval

        integer :: count,i,j
        real(fp_kind) :: t,tout

        if(allocated(output_cell_list)) deallocate(output_cell_list)
        allocate(output_cell_list(maxval(ncells)*n_slices))
        output_cell_list = .false.

        if(allocated(cell_map)) deallocate(cell_map)
        allocate(cell_map(maxval(ncells)*n_slices))
        cell_map = 0

        t = 0.0_fp_kind
        tout = thickness_interval

        count = 0

        do i= 1,maxval(ncells)
        do j=1,n_slices
            if((t+depths(j+1)*a0(3).gt.tout)) then
                output_cell_list((i-1)*n_slices+j) =.true.
                count = count + 1
                tout = (floor((t+depths(j+1)*a0(3))/thickness_interval)+1)*thickness_interval
            endif

        enddo
        t = t+a0(3)
        enddo
        if(allocated(output_thickness_list)) deallocate(output_thickness_list)
        allocate(output_thickness_list(count))

        t = 0.0_fp_kind
        tout = thickness_interval

        count = 0
        do i= 1,maxval(ncells)
        do j=1,n_slices
            if((t+depths(j+1)*a0(3).gt.tout)) then
                count = count + 1
                output_thickness_list(count) = t+depths(j+1)*a0(3)
                cell_map((i-1)*n_slices+j) = count
                tout = (floor((t+depths(j+1)*a0(3))/thickness_interval)+1)*thickness_interval
            endif
        enddo
        t = t+a0(3)
        enddo

    end subroutine



    subroutine write_thicknesss_to_file

        use output, only: output_prefix

        implicit none

        integer :: i,j
        character(1024) :: filename,dir,fnam

        j = index(output_prefix,'/',back=.true.)
        j = max(j,index(output_prefix,'\\',back=.true.))

        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'\Probe_intensity'//trim(adjustl(fnam))//'_probe_intensity_thicknesss.txt'
        else
            filename = 'Probe_intensity\\'//trim(adjustl(output_prefix))//'_probe_intensity_thicknesss.txt'
        endif

        write(*,*) 'The thicknesses at which the probe intensity'
        write(*,*) 'is being outputted have been written to'
        write(*,*)
        write(*,*) '  ' // trim(filename)
        write(*,*)

        open(unit=8734, file=filename)

        do i = 1, size(output_thickness_list)
            write(8734, *) output_thickness_list(i)
        enddo

        close(8734)



    end subroutine

    subroutine probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample)

    use output, only: output_prefix,quad_shift
    use global_variables, only: nopiy,nopix
    use m_string

    real(fp_kind),intent(in)::probe_intensity(nopiy,nopix,size(output_thickness_list))
    integer*4,intent(in)::i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample

    integer*4::j,z
    character*1024::filename,fnam,dir

        j = index(output_prefix,'/',back=.true.)
        j = max(j,index(output_prefix,'\\',back=.true.))
        z = size(output_thickness_list)

        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'\\Probe_intensity\\'//trim(adjustl(fnam))//'_ProbeIntensity'
        else
            filename = 'Probe_intensity\\'//trim(adjustl(output_prefix))//'_ProbeIntensity'
        endif

        if (probe_ndf.gt.1) filename = trim(filename) // '_df' // to_string(i_df)
        if (nysample.gt.1) filename = trim(filename) // '_ny' // to_string(ny)
        if (nxsample.gt.1) filename = trim(filename) // '_nx' // to_string(nx)
        filename = trim(filename) // '_'//to_string(nopiy)//'x'//to_string(nopix)//'x'//to_string(z)//'.bin'
        open(4985, file=filename, form='unformatted', convert='big_endian')
        do j=1,z
            write(4985) quad_shift(probe_intensity(:,:,j),nopiy,nopix)/ n_qep_passes
        enddo
        close(4985)

    end subroutine

    subroutine prompt_save_load_grates

        use m_user_input, only: get_input
        use output, only: output_prefix
        use m_string

        implicit none

        integer :: i_save_load, i_retry,i
        logical :: exists,retry

        i_save_load = -1
        do while(i_save_load<0.or.i_save_load>2)


        call command_line_title_box('Save/load transmission functions')
        write(*,*) 'Warning: the files outputted when saving may be very large.'
        write(*,*)
        write(*,*) '<0> Proceed without saving or loading'
        write(*,*) '<1> Save transmission functions'
        write(*,*) '<2> Load transmission functions'
        write(*,*) '<3> Add additional transmission function '
        write(*,*) '    (eg. from magnetic structure) from file'

		call get_input('<0> continue <1> save <2> load', i_save_load)

323     format( '  - The xtl file',/,'  - The slicing of the unit cell',/,&
               &'  - The choice of thermal scattering model (QEP vs. absorptive)',/,&
               &'  - The tiling of the unit cell',/,'  - The number of pixels',/,&
               &'  - (For absorptive model: whether absorption is included)',/,&
               &'  - (For QEP model: the number of distinct transmission functions)',/,&
               &'  - (For QEP model: phase ramp shift choice)',/)

        select case (i_save_load)
            case (0)
                return

            case (1)
                save_grates = .true.
                grates_filename = trim(adjustl(output_prefix)) // '_transmission_functions.bin'

                write(*,*) 'The transmission functions will be saved to the file'
                write(*,*)
                write(*,*) '    ' // trim(grates_filename)
                write(*,*)
                write(*,*) 'They can be loaded for later calculations provided'
                write(*,*) 'the following parameters are identical:'
                write(6,323)

            case (2)
                write(*,*) 'It is up to the user to ensure that the parameters used'
                write(*,*) 'to create the loaded transmission functions are consistent'
                write(*,*) 'with those of the current calculation:'
                write(6,323)
                retry=.true.
                do while(retry)
               write(*,*) 'Enter filename of transmission functions:'
                call get_input('filename of transmission functions', grates_filename)
                write(*,*)

                inquire(file=grates_filename, exist=exists)
                retry=.not.exists
                if (.not.exists) then
                    write(*,*) 'ERROR: cannot find this file.'
                    write(*,*) '<1> Enter again'
                    write(*,*) '<2> Proceed without loading'
                    call get_input('<1> Enter again <2> Proceed without loading', i_retry)
                    write(*,*)
                    retry = i_retry==1
                    if(.not.retry) return


                endif
                enddo
                load_grates = .true.
            case (3)
                additional_transmission_function=.true.
                pure_phase=.true.

                if(.not.allocated(amplitude_fnam)) allocate(amplitude_fnam(n_slices),phase_fnam(n_slices))
                do i=1,n_slices
                    if(.not.pure_phase) then
                        write(*,*) char(10),' Please input filename of amplitude of additional transmission function for slice ',i
                        if(i==1) then
                            write(*,*) 'If your intended transmission function is a pure phase object (ie. the '
                            write(*,*) 'amplitude of theadditional transmission function is everywhere 1), input'
                            write(*,*) '-1 and no amplitude file willbe loaded nor will you be prompted for the '
                            write(*,*) 'amplitude for the remaining slices'
                        endif
                        call get_input('Amplitude of additional transmission function',amplitude_fnam(i))
                        pure_phase = trim(adjustl(amplitude_fnam(i)))=='-1'
                    endif
                    write(*,*) char(10),' Please input filename of phase of additional transmission function for slice ',i
                    call get_input('Phase of additional transmission function',phase_fnam(i))
                enddo



        end select
        enddo
    end subroutine

 subroutine load_save_add_grates_qep(idum,qep_grates,nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice)
        use m_numerical_tools, only: gasdev
        use output
        use global_variables, only: ak,pi
        integer(4),intent(inout):: idum
        complex(fp_kind),intent(inout)::qep_grates(nopiy,nopix,n_qep_grates,n_slices)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice(nt,n_slices)

        integer*4::j,i,m,n
        real(fp_kind)::junk

        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)


        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            read(3984) qep_grates
            close(3984)

            ! Make sure the random number sequence is as if grates were calculated
            ! So call gasdev as many times as it would usually be called
            do j = 1, n_slices;do i = 1, n_qep_grates;do m = 1, nt;do n = 1, nat_slice(m,j)*2
                junk = gasdev(idum)
            enddo;enddo;enddo;enddo

            return
        endif

        if(additional_transmission_function) then
            write(*,*) 'Adding addition transmission function to file...'
            amplitude = 1
            do j=1,n_slices
                if(.not.pure_phase) call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                call binary_in(nopiy,nopix,phase,phase_fnam(j))
                qep_grates(:,:,:,j) = qep_grates(:,:,:,j)+spread(transpose(phase),dim=3,ncopies=n_qep_grates)/pi/a0_slice(3,j)*ak
                if(.not.pure_phase) then
                    call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                    qep_grates(:,:,:,j) = qep_grates(:,:,:,j)&
                                        +spread(cmplx(0_fp_kind,log(transpose(amplitude))),dim=3,ncopies=n_qep_grates)
                endif
            enddo
            qep_mode=4
        endif

        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            write(3984) qep_grates*pi*spread(spread(spread(a0_slice(3,1:n_slices),dim=1,ncopies=nopiy),dim=2,ncopies=nopix),&
                                                                                &dim=3,ncopies=n_qep_grates)/ak;close(3984)
        endif
    end subroutine

subroutine load_save_add_grates_abs(abs_grates,nopiy,nopix,n_slices)
        use global_variables, only: pi,ak
        use m_numerical_tools, only: gasdev
        use output
        complex(fp_kind),intent(inout)::abs_grates(nopiy,nopix,n_slices)
        integer*4,intent(in)::nopiy,nopix,n_slices

        integer*4::j

        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)


        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            read(3984) abs_grates
            close(3984)
            return
        endif

        if(additional_transmission_function) then
            write(*,*) 'Adding addition transmission function from file...'
            amplitude = 1
            do j=1,n_slices
                call binary_in(nopiy,nopix,phase,phase_fnam(j))
                abs_grates(:,:,j) = abs_grates(:,:,j)+transpose(phase)/pi/a0_slice(3,j)*ak
                if(.not.pure_phase) then
                    call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                    abs_grates(:,:,j) = abs_grates(:,:,j)+cmplx(0_fp_kind,log(transpose(amplitude)))
                endif
            enddo
        endif

        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            write(3984) abs_grates*pi*spread(spread(a0_slice(3,1:n_slices),dim=1,ncopies=nopiy),dim=2,ncopies=nopix)/ak;close(3984)
        endif
    end subroutine


    function make_detector(nopiy,nopix,ifactory,ifactorx,ss,kmin,kmax,phi,delphi) result(detector)
        use m_crystallography
        !makes a detector on an array nopiy x nopix with Fourier space pixel size deltaky x deltakx
        !which measures from kmin to kmax, supplying phi and delphi (detector orientation
        !and angular range in radians) will create a detector segment
        integer*4,intent(in)::nopiy,nopix,ifactory,ifactorx
        real(fp_kind),intent(in)::kmin,kmax,phi,delphi,ss(7)
        optional::phi,delphi

        real(fp_kind)::phi_(2)=0_fp_kind,g_vec_array(3,nopiy,nopix),ang,kabs,detector(nopiy,nopix),pi
        integer*4::y,x
        logical::segment,notwrapped

        pi = 4.0_fp_kind*atan(1.0_fp_kind)
        !If the phi and delphi variables are present, then make a segmented detector
        segment = present(phi).and.present(delphi)
        if(segment) then
            !Force angle to be between 0 and 2pi
            phi_= mod([phi,phi+delphi],2*pi)
            do y=1,2
                if (phi_(y)<0) phi_(y) = phi_(y) + 2*pi
            enddo
            notwrapped = phi_(2)>phi_(1)
        endif
        detector = 0

        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        detector = 0
        do y=1,nopiy;do x=1,nopix
            kabs = trimr(g_vec_array(:,y,x),ss)
            if((kabs.ge.kmin).and.(kabs.le.kmax)) then
                detector(y,x) =1
                if(segment) then
                    !Force angle to be between 0 and 2pi
                    ang = modulo(atan2(g_vec_array(1,y,x),g_vec_array(2,y,x)),2*pi)
                    if (notwrapped) then
                        if(.not.(ang>phi_(1).and.ang<=phi_(2))) detector(y,x)=0
                    else
                        if(.not.(ang>phi_(1).or.ang<=phi_(2))) detector(y,x)=0
                    endif
            endif;endif
            enddo;enddo;
    end function

    subroutine setup_qep_parameters(n_qep_grates,n_qep_passes,nran,quick_shift,ifactory,ifactorx)

        use m_user_input, only: get_input
        use m_string, only: to_string,command_line_title_box

        implicit none

        integer*4,intent(out)::n_qep_grates,n_qep_passes,nran
        integer*4,intent(in)::ifactory,ifactorx
        logical,intent(in)::quick_shift


        integer :: i
        call command_line_title_box('QEP model parameters')

        qep_mode=1
        do while(qep_mode<2.or.qep_mode>4)
        write(*,*) 'Enter the number of distinct transmission functions to calculate:'
        call get_input("Number of phase grates calculated", n_qep_grates)
        write(*,*)

            write(*,*) 'As implemented in muSTEM, the QEP multislice algorithm will randomly shift the'
            write(*,*) 'unit cells in your supercell around to effectively generate new random '
            write(*,*) 'transmission functions. This means the QEP calculation samples a larger '
            write(*,*) 'effective number of random transmission functions than would otherwise be the'
            write(*,*) 'case.'
            write(*,*)
            write(*,*) 'If your choice of unit cell tiling and pixel grid size is such that each unit '
            write(*,*) 'cell is an integer number of pixels, unit cell shifting can be achieved by '
            write(*,*) 'circular shift of the transmission function arrays in memory (quick shifting).'
            write(*,*) 'Otherwise the unit cells will have to be shifted with sub-pixel precision '
            write(*,*) 'using the Fourier shift algorithm which will impact calculation speed. '
            write(*,*)
            write(*,*) 'A third option is to not shift the unit cells around at all, but the user is'
            write(*,*) 'advised that this requires a much larger number of distinct transmission '
            write(*,*) 'functions to achieve convergence than would usually be the case.'
            write(*,*)
        if (quick_shift) then
            write(6,100) to_string(ifactorx*ifactory), to_string(n_qep_grates), to_string(ifactorx*ifactory*n_qep_grates)
100         format(' The choice of tiling and grid size permits quick shifting.',/,&
                  &' The effective number of transmission functions used in  ',/,&
                  &' calculations will be ', a, ' * ', a, ' = ', a, '.', /)
            qep_mode=2
        else
            write(6,101)
        101 format( ' Your choice of tiling and grid size does not permit quick shifting ', /, &
                    &' of the precalculated transmission functions. Shifting using the ', /, &
                    &' Fourier shift algorithm can be performed but is time consuming. ', /, &
                    &' You may wish to go back and calculate more distinct transmission ', /, &
                    &' functions, or simply proceed without using phase ramp shifting. ' /)

            i=-1
            do while(i<1.or.i>3)
                write(6,111)
            111 format(  ' <1> Go back and choose a larger number', /, &
                        &' <2> Proceed with phase ramp shifting', /, &
                        &' <3> Proceed without phase ramp shifting', / &
                        )
                call get_input("<1> choose more <2> phase ramp <3> no phase ramp", i)
                write(*,*)

                if (i.eq.2 .and. (ifactory.gt.1 .or. ifactorx.gt.1)) then
                    qep_mode=3
                    call setup_phase_ramp_shifts

                elseif (i.eq.3) then
                    qep_mode=4
                endif
            enddo
        endif
        enddo
        write(6,*) 'Enter the number of passes to perform for QEP calculation:'
        write(*,*) 'Warning: using only a single pass is usually NOT sufficient.'
        call get_input("Number of Monte Carlo calculated", n_qep_passes )
        write(*,*)

        write(6,*) 'A pseudo-random number generator is used to randomly select'
        write(6,*) 'atomic displacements in the QEP model. You are given the'
        write(6,*) 'option of choosing the seed for this random number generator.'
        write(6,*) 'Any integer number is acceptable and you will need to reuse'
        write(6,*) 'this seed if you want to reproduce calculations.'
        call get_input("Rrandom number generator seed", nran )
        write(*,*)

        end subroutine

        subroutine setup_phase_ramp_shifts
            use global_variables, only: nopiy, nopix, ifactory, ifactorx
            implicit none

            integer(4) :: i
            real(fp_kind) :: r_coord

            if(allocated(shift_arrayy)) deallocate(shift_arrayy)
            if(allocated(shift_arrayx)) deallocate(shift_arrayx)
            allocate(shift_arrayy(nopiy,ifactory))
            allocate(shift_arrayx(nopix,ifactorx))

            do i = 1,ifactory
                r_coord=float(i)/float(ifactory)
                call make_shift_oned(shift_arrayy(:,i),nopiy,r_coord)
            enddo

            do i = 1,ifactorx
                r_coord=float(i)/float(ifactorx)
                call make_shift_oned(shift_arrayx(:,i),nopix,r_coord)
            enddo

        end subroutine
    subroutine setup_slicing_depths()

        use m_user_input, only: get_input
        use m_precision, only: fp_kind
        use m_string,only:command_line_title_box

        implicit none

        integer :: i_slice
        call command_line_title_box('Unit cell slicing')

    22  write(6,23) char(143)
    23  format(' Do you wish to slice the unit cell in the beam direction?', /, &
              &' This may be a good idea if the unit cell is larger than 2 ', a1, /, &
              &' in the z-direction.', /, &
              &' <1> Yes <2> No ')
        call get_input('Slice unit cell <1> yes <2> no', i_slice)
        write(*,*)

        if (i_slice.eq.1) then
            call calculate_depths_slicing

        elseif (i_slice.eq.2) then
            call calculate_depths_no_slicing

        else
            goto 22

        endif

    end subroutine

    subroutine calculate_depths_no_slicing
        implicit none

        n_slices = 1

        allocate(depths(2))

        depths(1) = 0.0_fp_kind
        depths(2) = 1.0_fp_kind

    end subroutine

    subroutine calculate_depths_slicing
        use m_user_input

        implicit none

        integer(4) i, ichoice

        n_slices = -1
        do while(n_slices<1)
        write(*,*) 'Enter the number of slices per unit cell in the beam direction:'
        call get_input("Number of distinct potentials ", n_slices)
        write(*,*)

        if (n_slices.eq.1) then; call calculate_depths_no_slicing
        elseif (n_slices.gt.1) then
            if (allocated(depths)) deallocate(depths)
            allocate(depths(n_slices+1))

            depths(n_slices+1) = 1.0_fp_kind

            write(6,10)
         10 format( ' You will now be asked to enter the depths at which slicing',/,&
                    &' will take place. These should be entered as fractions of',/,&
                    &' the unit cell. It is the front of the slice which should',/,&
                    &' be entered. (e.g. three even slices would be entered as',/,&
                    &' 0.0, 0.3333, 0.6666 ).', /)

            ichoice = 0
            do while(ichoice<1.or.ichoice>2)
            write(6,16)
    16      format(' How would you like to specify the slice depths?', /, &
                  &' <1> Manually ', /, &
                  &' <2> Automatically (uniformly spaced)')
            call get_input("<1> manual or <2> auto slicing", ichoice)
            write(*,*)

            if(ichoice.eq.1) then
                ! Manual slicing

                do i = 1, n_slices
                    write(6,20,ADVANCE='NO') achar(13), i
        20          format( a1,' Enter the fractional depth of slice number ', i4, ':  ')
                    call get_input("depths", depths(i))
                enddo

            elseif (ichoice.eq.2) then
                ! Automatic slicing

    21          format( ' Fractional depth of slice ', i4, ':  ', f7.4)
                do i = 1, n_slices
                    depths(i) = float( (i-1)) / float(n_slices)
                    write(6,21) i, depths(i)
                enddo
            endif
            enddo
        endif
        enddo

        write(*,*)

    end subroutine

	function calculate_nat_slice( depths, tau, nm, nt, nat, n_slices, ifactorx, ifactory ) result(nat_slice)
   !Calculate the number of atoms within each slice
   !depths   - an array describing the depth of each unit cell
   !tau      - an array dimensions 3 x # atom types x max number of atoms of a given type
   !           containing the fractional coordinates of all the atoms in the unit cel
   !nm       - max number of atom of a given type
   !nt       - Number of different types of atom
   !nat      - An array dimensions #atom types which contains the number of atoms of each type
   !n_slices - The sub-slicing chosen by the user
   !ifactory - Unit cell tiling in y direction
   !ifactorx - Unit cell tiling in x direction
    implicit none

	real(fp_kind),intent(in):: depths(n_slices), tau(3,nt,nm)
	integer(4),intent(in):: nm, nt, ifactorx, ifactory, n_slices,nat(nt)

	integer(4) nat_slice(nt,n_slices),i,kk
	real(fp_kind) diff, tol,factorx, factory

	tol = 1.0d-4

	do i = 1, nt !Loop over atom types
		!Calculate atoms in all slices bar the final one
		do kk=1,n_slices-1;
			nat_slice(i,kk) = count((tau(3,i,1:nat(i)) .lt. (depths(kk+1)-tol)) &
							&.and.(tau(3,i,1:nat(i)) .ge. (depths(kk)-tol))  )&
							&*ifactory*ifactorx

		enddo
		!Calculate the number of atoms in the final slice
		nat_slice(i,n_slices) = count(tau(3,i,1:nat(i)) .ge. (depths(n_slices)-tol))*ifactory*ifactorx


	enddo
    end function

    subroutine calculate_slices
    !This subroutine is called after the slicing depths and supercell tiling has been chosen
    !It allocates the global variable arrays that contain the fractional coordinates of all the
    !atoms and the arrays which describe the dimensions of the atoms
        use global_variables, only: nat, ifactory, ifactorx, nt, a0, deg, tau, nm,even_slicing
        use m_crystallography, only: cryst

        implicit none

        integer :: j,jj,kk,maxnat_slice_uc

        !Allocate the arrays which contain the number of atoms in the supercell slices
        !and the single unit cell slices
        allocate(nat_slice(nt,n_slices),nat_slice_unitcell(nt,n_slices))
        nat_slice = calculate_nat_slice( depths,tau,nm,nt,nat,n_slices,ifactorx,ifactory )
		nat_slice_unitcell = nat_slice / dfloat(ifactorx*ifactory)
		maxnat_slice=maxval(nat_slice)
		maxnat_slice_uc=maxnat_slice / dfloat(ifactorx*ifactory)

        allocate(tau_slice(3,nt,maxnat_slice,n_slices),&
                &a0_slice(3,n_slices),ss_slice(7,n_slices),prop_distance(n_slices),&
                 tau_slice_unitcell(3,nt,maxnat_slice_uc,n_slices))

        do j = 1, n_slices
          a0_slice(1:2,j) = a0(1:2)*[ifactorx,ifactory]

          if (j .eq. n_slices) a0_slice(3,j) = (1.0_fp_kind-depths(j))*a0(3)
          if (j .lt. n_slices) a0_slice(3,j) = (depths(j+1)-depths(j))*a0(3)

          call cryst(a0_slice(:,j),deg,ss_slice(:,j))

          ! Make supercell cell subsliced
          call calculate_tau_nat_slice( depths(j),depths(j+1),tau,tau_slice(:,:,:,j),nm,nt,nat,nat_slice(:,j),ifactorx,ifactory,&
		  &                             maxnat_slice )

          ! Make unit cell subsliced
          call make_mod_tau_unitcell( depths(j),depths(j+1),tau,tau_slice_unitcell(:,:,:,j),nm,nt,nat,nat_slice_unitcell(:,j) )

        enddo

        prop_distance = a0_slice(3,:)
		 even_slicing = all(abs(prop_distance - sum(prop_distance)/n_slices)<1e-3)

    end subroutine

	subroutine calculate_tau_nat_slice( depth1, depth2, tau, tau_slice, nm, nt, nat, nat_slice, ifactorx, ifactory, maxnat_slice )
	!Calculate
    implicit none

	integer(4) nm, nt, ifactorx, ifactory, maxnat_slice
	integer(4) nat(nt)
	integer(4) nat_slice(nt)
	real(fp_kind) depth2, depth1, tau(3,nt,nm)
	real(fp_kind) tau_slice(3,nt,maxnat_slice)
	integer(4) i,j,jj, mm, nn
	real(fp_kind) diff, tol
	real(fp_kind) factorx, factory

	tol = 1.0d-4 !Numerical Tolerance
	diff = depth2 - depth1 !Slice thickness (difference)

	factorx = float( ifactorx ) !Make floating point versions of ifactory and ifactorx
	factory = float( ifactory )

	do i = 1, nt

	   jj = 1

	   do mm = 1, ifactory
	   do nn = 1, ifactorx
	      do j = 1, nat(i)

	         if( (tau(3,i,j) .lt. (depth2-tol)) .and.(tau(3,i,j) .ge. (depth1-tol)) ) then
	            tau_slice(1:2,i,jj) = tau(1:2,i,j)/[factorx,factory] + [float(nn-1)/factorx,float(mm-1)/factory]
			    tau_slice(3,i,jj) = (tau(3,i,j)-depth1)/diff
	            jj = jj + 1

	         elseif(diff.eq.1.0) then
	            tau_slice(1:2,i,jj) = tau(1:2,i,j)/[factorx,factory] + [float(nn-1)/factorx,float(mm-1)/factory]
			    tau_slice(3,i,jj) = (tau(3,i,j)-depth1)/diff

	         endif
	      enddo
		enddo
	    enddo
	    nat_slice(i) = jj-1

	enddo

    end subroutine

!--------------------------------------------------------------------------------------
    subroutine make_mod_tau_unitcell( depth1, depth2, tau, mod_tau, nm, nt, nat, nat2)
    !subroutine make_mod_tau_force() is called if natural depths
      !were NOT chosen.  In this case the current and next depth are
      !required to select which of the tau fit within this slice.
      !Note that a tolerance is employed -- this measures how far
      !beyond the limits an atom can be and still be counted.
      !CAREFUL to realise that false holz can be an issue using this forced
      !slicing.
    use m_precision
    implicit none

	integer(4) nm, nt
	integer(4) nat(nt),nat2(nt)
	real(fp_kind) depth2, depth1, tau(3,nt,nm),mod_tau(3,nt,nm)
	integer(4) i,j,jj
	real(fp_kind) diff, tol
	tol = 1.0e-4_fp_kind

	diff = depth2-depth1


	do i=1,nt
	    jj=1
	    do j=1,nat(i)
	        if( (tau(3,i,j) .lt. (depth2-tol)) .and.(tau(3,i,j) .ge. (depth1-tol)) ) then
	            mod_tau(1:2,i,jj) = tau(1:2,i,j)
			    mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff
			    jj = jj+1
	        elseif(diff.eq.1.0) then
	            mod_tau(1:2,i,jj) = tau(1:2,i,j)
			    mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff
	        endif
	    enddo
	    nat2(i)=jj-1

	enddo

	end subroutine

    subroutine get_cbed_detector_info(ndet,k,outer,inner)

        use m_precision, only: fp_kind
        use m_user_input, only: get_input
        use m_string, only: to_string

        implicit none

        integer(4)   ndet,ichoice,i,mrad
        real(fp_kind),dimension(ndet) :: outer,inner,inner_mrad,outer_mrad
		optional::inner
        real(fp_kind) :: k,dummy
		logical::getinner

		getinner = present(inner)

        write(*,*) 'Select a method for choosing detector angles:'
        write(*,*) '<1> Manual',char(10),' <2> Automatic'

        call get_input("manual detector <1> auto <2>",ichoice)
        write(*,*)
        mrad =-1
		do while(.not.((mrad==1).or.(mrad==2)))
			write(*,*) 'Select how angles will be specified:'
			write(*,*) '<1> mrad',char(10),' <2> inverse Angstroms'
			call get_input("<1> mrad <2> inv A", mrad)
			write(*,*)
		enddo

        if(ichoice.eq.1) then
              do i = 1, ndet
                    write(*,*) char(10),' Detector ', to_string(i)
					if(getinner) then
						write(*,*)' Inner angle:'
						call get_input("inner",dummy)
						if(mrad.eq.1) inner_mrad(i) = dummy
						if(mrad.eq.2) inner(i) =  dummy
					endif
					write(*,*) "Outer angle:"
                    call get_input("outer",dummy)
                    if(mrad.eq.1) outer_mrad(i) = dummy
					if(mrad.eq.2) outer(i) =  dummy
              enddo
        else
			  if(getinner) then
				  write(*,*) "Initial inner angle:"
				  call get_input("initial inner angle", dummy)
				  if(mrad.eq.1) inner_mrad(1) = dummy
				  if(mrad.eq.2) inner(1) = dummy
			  endif

              write(*,*) "Initial outer angle:"
              call get_input("initial outer angle", dummy)
			  if(mrad.eq.1) outer_mrad(1) = dummy
			  if(mrad.eq.2) outer(1) = dummy

              if(getinner) write(*,*) "Increment (both angles incremented by this amount):"
			  if(.not.getinner)  write(*,*) "Increment:"
              call get_input("increment", dummy)
              write(*,*)

			  if(mrad.eq.1) then
				if(getinner) inner_mrad(2:) = (/((i-1)*dummy+inner_mrad(1), i=2,ndet,1)/)
				outer_mrad(2:) = (/((i-1)*dummy+outer_mrad(1), i=2,ndet)/)
			  else
				if(getinner) inner(2:) = (/((i-1)*dummy+inner(1), i=2,ndet)/)
				outer(2:) = (/((i-1)*dummy+outer(1), i=2,ndet)/)
			  endif

        endif
		if(mrad.eq.2) then
			outer_mrad = 1000*atan(outer/k)
			if(getinner) inner_mrad = 1000*atan(inner/k)
		else
			if(getinner) inner = k*tan(inner_mrad/1000.0_fp_kind)
			outer = k*tan(outer_mrad/1000.0_fp_kind)
		endif

        write(*,*) 'Summary of diffraction plane detectors:'
        write(*,*)

        if(getinner) write(*,*) '         inner    outer'
		if(.not.getinner) write(*,*) ' outer angle'
        write(*,*) '  --------------------------------'
        if(getinner) then
		do i = 1, ndet
            write(*,50) i, inner_mrad(i), outer_mrad(i)
            write(*,55) inner(i), outer(i), char(143)
            write(*,60)
        enddo
		else
		do i = 1, ndet
            write(*,65) i, outer_mrad(i)
            write(*,70) outer(i), char(143)
            write(*,60)
        enddo
		endif

50          format(1x, i5, ' | ', f6.2, ' | ', f6.2, '   (mrad)')
55          format(1x, 5x, ' | ', f6.2, ' | ', f6.2, '   (', a1, '^-1)')
60          format(1x, 5x, ' | ', 6x, ' | ')
65          format(1x, i5, ' | ', f6.2,  '   (mrad)')
70          format(1x, 5x, ' | ', f6.2, '   (', a1, '^-1)')

        write(*,*)



    end subroutine


	function make_image(nopiy,nopix,psi,ctf,rspacein) result(image)
		use FFTW3
		implicit none
		integer*4,intent(in)::nopiy,nopix
		complex(fp_kind),dimension(nopiy,nopix),intent(in)::ctf,psi
		logical,intent(in),optional::rspacein
		real(fp_kind)::image(nopiy,nopix)

		complex(fp_kind),dimension(nopiy,nopix)::psi_temp2,psi_temp

		if(present(rspacein)) then
			if(rspacein)call fft2(nopiy,nopix,psi,nopiy,psi_temp,nopiy)
			if(.not.rspacein) psi_temp = psi
		else
			call fft2(nopiy,nopix,psi,nopiy,psi_temp,nopiy)
		endif

		psi_temp = psi_temp*ctf
		call ifft2(nopiy,nopix,psi_temp,nopiy,psi_temp2,nopiy)
		image = abs(psi_temp2)**2

	end function

end module


    do i=1,nt
        jj=1
        do j=1,nat(i)
            if( (tau(3,i,j) .lt. (depth2-tol)) .and.(tau(3,i,j) .ge. (depth1-tol)) ) then
                mod_tau(1,i,jj) = tau(1,i,j)
                mod_tau(2,i,jj) = tau(2,i,j)
                mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff
                jj = jj+1
            elseif(diff.eq.1.0) then
                mod_tau(1,i,jj) = tau(1,i,j)
                mod_tau(2,i,jj) = tau(2,i,j)
                mod_tau(3,i,jj) = (tau(3,i,j)-depth1)/diff
            endif
        enddo
        nat2(i)=jj-1
    enddo

    end subroutine



end module
