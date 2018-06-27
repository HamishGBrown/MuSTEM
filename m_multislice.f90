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
    
    logical :: output_probe_intensity = .false.,additional_transmission_function = .false.,pure_phase
    logical,allocatable :: output_cell_list(:)
    real(fp_kind),allocatable :: output_thickness_list(:)
    integer,allocatable :: cell_map(:)
     
    character*200,allocatable::amplitude_fnam(:),phase_fnam(:) 
    

    complex(fp_kind),allocatable,dimension(:,:) :: shift_arrayx, shift_arrayy        
    
    
    interface load_save_add_grates
        module procedure load_save_add_grates_qep,load_save_add_grates_abs
    end interface
    
    contains
    
    !This subroutine samples from the available phase grates and then performs one iteration of the multislice algorithm (called in CPU versions only)
    
    subroutine qep_multislice_iteration(psi,propagator,transmission,nopiy,nopix,ifactory,ifactorx,idum,n_qep_grates,mode,shift_arrayy,shift_arrayx)
        use m_numerical_tools,only:ran1
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix,n_qep_grates),shift_arrayy(nopiy),shift_arrayx(nopix)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,mode,ifactory,ifactorx
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
            call multislice_iteration(psi,propagator,trans,nopiy,nopix)
        elseif(mode == 3) then      !Phase ramp shift
            shiftx = floor(ifactorx*ran1(idum)) + 1
            shifty = floor(ifactory*ran1(idum)) + 1
			call phase_shift_array(transmission(:,:,nran),trans,shift_arrayy,shift_arrayx)
            call multislice_iteration(psi,propagator,trans,nopiy,nopix)
        else
            call multislice_iteration(psi,propagator,transmission(:,:,nran),nopiy,nopix)
        endif
    end subroutine
    
    !This subroutine performs one iteration of the multislice algorithm (called in CPU versions only)
    !Probe (psi) input and output is in real space
    subroutine multislice_iteration(psi,propagator,transmission,nopiy,nopix)
	    use cufft_wrapper
        complex(fp_kind),intent(inout)::psi(nopiy,nopix)
        complex(fp_kind),intent(in)::propagator(nopiy,nopix),transmission(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix
        
        ! Transmit through slice potential
		psi = psi*transmission
                
        ! Propagate to next slice
		call fft2(nopiy,nopix,psi,nopiy,psi,nopiy)
		psi = psi*propagator
		call ifft2(nopiy,nopix,psi,nopiy,psi,nopiy)
    
    end subroutine
    
    subroutine prompt_output_probe_intensity
    
        use m_user_input, only: get_input
        use global_variables, only: thickness, nopiy, nopix
        use output, only: output_prefix,split_filepath
        use m_string, only: to_string
        
        implicit none
        
        integer :: i_output,j
        real(fp_kind) :: thickness_interval
        character(1024)::dir,fnam
        
        write(*,*) '|-----------------------------------|'
        write(*,*) '|      Output probe intensity       |'
        write(*,*) '|-----------------------------------|'
        write(*,*)
        write(*,*) 'The probe intensity as a function of thickness'
        write(*,*) 'can be outputted to file at each probe position.'
        write(*,*) 'The user is advised that the outputted dataset'
        write(*,*) 'may be very large.'
        write(*,*) 
10      write(*,*) '<0> Proceed'
        write(*,*) '<1> Output probe intensity'
        call get_input('<0> Proceed <1> Output probe intensity', i_output)
        write(*,*)
        
        output_probe_intensity = i_output ==1
        
        if(output_probe_intensity) then
                
20              write(*,*) 'At what thickness interval (in Angstroms)'
                write(*,*) 'should intensities be outputted?'
                call get_input('At what thickness interval should intensities be outputted?', thickness_interval)
                write(*,*)
                
                if (thickness_interval.le.0.0_fp_kind .or. thickness_interval.gt.thickness) then
                    write(*,*) 'ERROR: invalid thickness.'
                    goto 20
                endif
                call split_filepath(output_prefix,dir,fnam)
                
                call system('mkdir '//trim(adjustl(dir))//'Probe_intensity')
                
                call generate_cell_list(thickness_interval)
                
                call write_thicknesss_to_file
                
                write(*,*) 'The probe intensities will be written to the files'
                write(*,*)
                write(*,*) '  '//trim(adjustl(dir))//'Probe_intensity' // trim(adjustl(fnam)) // '_ProbeIntensity*.bin'
                write(*,*)
                if (fp_kind.eq.4) then
                    write(*,*) 'as 32-bit big-endian floats.'
                elseif (fp_kind.eq.8) then
                    write(*,*) 'as 64-bit big-endian floats.'
                endif
                write(*,*) 'Each file contains a sequence of ' // to_string(nopiy) // 'x' // to_string(nopix) // ' arrays.'
                write(*,*)
        end if
        
    end subroutine
    
    
    
    subroutine generate_cell_list(thickness_interval)
    
        use global_variables, only: thickness, ncells, a0
		use m_slicing, only: n_slices,depths
        
        implicit none
        
        real(fp_kind) :: thickness_interval
        
        integer :: i_cell, count,i,j
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
        j = max(j,index(output_prefix,'\',back=.true.))
        
        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'Probe_intensity'//trim(adjustl(fnam))//'_probe_intensity_thicknesss.txt'
        else
            filename = 'Probe_intensity\'//trim(adjustl(output_prefix))//'_probe_intensity_thicknesss.txt'
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
    
    subroutine probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes)
    
    use output, only: output_prefix,quad_shift
    use global_variables, only: nopiy,nopix
    use m_probe_scan, only: nysample,nxsample
    use m_lens, only:probe_ndf
    use m_string
    
    real(fp_kind),intent(in)::probe_intensity(nopiy,nopix,size(output_thickness_list))
    integer*4,intent(in)::i_df,ny,nx,n_qep_passes
    
    integer*4::j
    character*1024::filename,fnam,dir
    
        j = index(output_prefix,'/',back=.true.)
        j = max(j,index(output_prefix,'\',back=.true.))
        
        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'Probe_intensity\'//trim(adjustl(fnam))//'_ProbeIntensity'
        else
            filename = 'Probe_intensity\'//trim(adjustl(output_prefix))//'_ProbeIntensity'
        endif
            
        if (probe_ndf.gt.1) filename = trim(filename) // '_df' // to_string(i_df)
        if (nysample.gt.1) filename = trim(filename) // '_ny' // to_string(ny)
        if (nxsample.gt.1) filename = trim(filename) // '_nx' // to_string(nx)
        filename = trim(filename) // '.bin'
        open(4985, file=filename, form='binary', convert='big_endian')
        do j=1,size(output_thickness_list)
            write(4985) quad_shift(probe_intensity(:,:,j),nopiy,nopix)/ n_qep_passes
        enddo
        close(4985)
    
    end subroutine
    
    subroutine prompt_save_load_grates
    
        use m_user_input, only: get_input
        use output, only: output_prefix
        use m_slicing, only: n_slices
        
        implicit none
        
        integer :: i_save_load, i_retry,i
        logical :: exists
        
        i_save_load = -1
        do while(i_save_load<0.or.i_save_load>2)
            
        write(*,*) '|---------------------------------------------|'
        write(*,*) '|      Save/load transmission functions       |'
        write(*,*) '|---------------------------------------------|'
        write(*,*)
        write(*,*) 'Warning: the files outputted when saving may be very large.'
        write(*,*)
        write(*,*) '<0> Proceed without saving or loading'
        write(*,*) '<1> Save transmission functions'
        write(*,*) '<2> Load transmission functions'
        write(*,*) '<3> Add additional transmission function '
        write(*,*) '    (eg. from magnetic structure) from file'
        call get_input('<0> continue <1> save <2> load', i_save_load)
        write(*,*)
        
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
                write(*,*) '  - The xtl file'
                write(*,*) '  - The slicing of the unit cell'
                write(*,*) '  - The choice of thermal scattering model (QEP vs. absorptive)'
                write(*,*) '  - The tiling of the unit cell'
                write(*,*) '  - The number of pixels'
                write(*,*) '  - (For absorptive model: whether absorption is included)'
                write(*,*) '  - (For QEP model: the number of distinct transmission functions)'
                write(*,*) '  - (For QEP model: phase ramp shift choice)'
                write(*,*)
                
            case (2)
                write(*,*) 'It is up to the user to ensure that the parameters used'
                write(*,*) 'to create the loaded transmission functions are consistent'
                write(*,*) 'with those of the current calculation:'
                write(*,*) '  - The xtl file'
                write(*,*) '  - The slicing of the unit cell'
                write(*,*) '  - The choice of thermal scattering model (QEP vs. absorptive)'
                write(*,*) '  - The tiling of the unit cell'
                write(*,*) '  - The number of pixels'
                write(*,*) '  - (For absorptive model: whether absorption is included)'
                write(*,*) '  - (For QEP model: the number of distinct transmission functions)'
                write(*,*) '  - (For QEP model: phase ramp shift choice)'
                write(*,*)
2               write(*,*) 'Enter filename of transmission functions:'
                call get_input('filename of transmission functions', grates_filename)
                write(*,*)
                
                inquire(file=grates_filename, exist=exists)
                
                if (.not.exists) then
                    write(*,*) 'ERROR: cannot find this file.'
3                   write(*,*) '<1> Enter again'
                    write(*,*) '<2> Proceed without loading'
                    call get_input('<1> Enter again <2> Proceed without loading', i_retry)
                    write(*,*)
                    
                    select case (i_retry)
                        case (1)
                            goto 2
                        case (2)
                            return
                        case default
                            goto 3
                    end select
                    
                endif
                
                load_grates = .true.
            case (3)
                additional_transmission_function=.true.
                pure_phase=.false.
                
                if(.not.allocated(amplitude_fnam)) allocate(amplitude_fnam(n_slices),phase_fnam(n_slices))
                do i=1,n_slices
                    if(.not.pure_phase) then
                        write(*,*) char(10),' Please input filename of amplitude of additional transmission function for slice ',i
                        if(i==1) then
                            write(*,*) 'If your intended transmission function is a pure phase object (ie. the amplitude of the'
                            write(*,*) 'additional transmission function is everywhere 1), input -1 and no amplitude file will'
                            write(*,*) 'be loaded nor will you be prompted for the amplitude for the remaining slices'
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
	    integer(4),intent(inout):: idum
        complex(fp_kind),intent(inout)::qep_grates(nopiy,nopix,n_qep_grates,n_slices)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,n_slices,nt,nat_slice(nt,n_slices)
        
        integer*4::j,i,m,n
        real(fp_kind)::junk
        
        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)
        
        
        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
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
                qep_grates(:,:,:,j) = qep_grates(:,:,:,j)*spread(transpose(amplitude)*exp( cmplx(0,1)*transpose(phase)),dim=3,ncopies=n_qep_grates)
            enddo
        endif
        
        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            write(3984) qep_grates;close(3984)
        endif    
    end subroutine

subroutine load_save_add_grates_abs(abs_grates,nopiy,nopix,n_slices)
        use m_numerical_tools, only: gasdev
        use output
        complex(fp_kind),intent(inout)::abs_grates(nopiy,nopix,n_slices)
        integer*4,intent(in)::nopiy,nopix,n_slices
        
        integer*4::j,i,m,n
        real(fp_kind)::junk
        
        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)
        
        
        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            read(3984) abs_grates
            close(3984)
            return
        endif
    
        if(additional_transmission_function) then
            write(*,*) 'Adding addition transmission function to file...'
            amplitude = 1
            do j=1,n_slices
                if(.not.pure_phase) call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                call binary_in(nopiy,nopix,phase,phase_fnam(j))
                abs_grates(:,:,j) = abs_grates(:,:,j)*transpose(amplitude)*exp( cmplx(0,1)*transpose(phase))
            enddo
        endif
        
        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            write(3984) abs_grates;close(3984)
        endif    
    end subroutine
    

	function make_detector(nopiy,nopix,kmin,kmax,deltaky,deltakx,phi,delphi) result(detector)
        use global_variables,only:pi
		!makes a detector on an array nopiy x nopix with Fourier space pixel size deltaky x deltakx
		!which measures from kmin to kmax, supplying phi and delphi (detector orientation 
		!and angular range in radians) will create a detector segment
		integer*4,intent(in)::nopiy,nopix
		real(fp_kind),intent(in)::kmin,kmax,deltaky,deltakx,phi,delphi
		optional::phi,delphi

		real(fp_kind)::detector(nopiy,nopix)

		real(fp_kind)::ky(nopiy),kx(nopix),vec(2),phi_(2)
		real(fp_kind),dimension(nopiy,nopix)::phigrid,kabs,phiarray
		integer*4::y,x,i,lower(nopiy,nopix),upper(nopiy,nopix)
		logical::segment

		!If the phi and delphi variables are present, then make a segmented detector
		segment = present(phi).and.present(delphi)
		detector = 0

		!Generate k space arrays
		ky = kspace_array(nopiy)*deltaky
		kx = kspace_array(nopix)*deltakx

		!Generate the radial part of the detector
		kabs = sqrt(real(spread(ky**2,dim=2,ncopies = nopix)&
		               &+spread(kx**2,dim=1,ncopies = nopiy),kind=fp_kind))
		detector = merge(1,0,(kabs.ge.kmin) .and. (kabs.le.kmax))
	

		!If not a segmented detector then return at this point
		if(.not.segment) return
        phiarray = atan2(real(spread(ky,dim = 2,ncopies = nopix)),real(spread(kx,dim = 1,ncopies = nopiy)))+pi
	    phi_= mod([phi,phi+delphi],2*pi)
        do i=1,2
            if (phi_(i)<0) phi_(i) = phi_(i) + 2*pi
        enddo

        upper = merge(1,0,phiarray>phi_(1))
	    lower = merge(1,0,phiarray<=phi_(2))
    
    
	    if (phi_(2)>phi_(1)) then
		    detector = detector*lower*upper
        
		    detector = detector*upper
	    else
		    detector = merge(detector,0.0_fp_kind,((lower==1).or.(upper==1)))
        endif

	end function   

	function kspace_array(nopiy) result(karray)

		integer*4,intent(in)::nopiy
		integer*4::karray(nopiy)
	
		integer*4::i
	
		karray = [((i - nopiy/2),i=0,nopiy-1)]
        karray = cshift(karray,-nopiy/2)
	end function     
    
    subroutine setup_qep_parameters(n_qep_grates,n_qep_passes,phase_ramp_shift,nran,quick_shift,ifactory,ifactorx)
    
        use m_user_input, only: get_input
        use m_string, only: to_string
        
        implicit none
        
        integer*4,intent(out)::n_qep_grates,n_qep_passes,nran
        logical,intent(out)::phase_ramp_shift
        integer*4,intent(in)::ifactory,ifactorx
        logical,intent(in)::quick_shift
        
    
        integer :: i
        
        write(*,*) '|--------------------------------|'
	    write(*,*) '|      QEP model parameters      |'
	    write(*,*) '|--------------------------------|'
        write(*,*)
             	
    10  write(*,*) 'Enter the number of distinct transmission functions to calculate:'
        call get_input("Number of phase grates calculated", n_qep_grates)
        write(*,*)
    
        if (quick_shift) then
            write(6,100) to_string(ifactorx*ifactory), to_string(n_qep_grates), to_string(ifactorx*ifactory*n_qep_grates)
    100     format(' The choice of tiling and grid size permits quick shifting.',/,&
                  &' The effective number of transmission functions used in  ',/,&
                  &' calculations will be ', a, ' * ', a, ' = ', a, '.', /)
        else
            write(6,101)
        101 format( ' The choice of tiling and grid size does not permit quick shifting ', /, &
                    &' of the precalculated transmission functions. Shifting using ', /, &
                    &' phase ramps can be performed but is time consuming. You may wish ', /, &
                    &' to go back and calculate more distinct transmission functions, or ', /, &
                    &' simply proceed without using phase ramp shifting. ' /)
                
        110 write(6,111)
        111 format(  ' <1> Go back and choose a larger number', /, &
                    &' <2> Proceed with phase ramp shifting', /, &
                    &' <3> Proceed without phase ramp shifting', / &                
                    )
            call get_input("<1> choose more <2> phase ramp <3> no phase ramp", i)
            write(*,*)
        
            if (i.eq.1) then
                goto 10
            
            elseif (i.eq.2 .and. (ifactory.gt.1 .or. ifactorx.gt.1)) then
                phase_ramp_shift = .true.
                
                call setup_phase_ramp_shifts
                        
            elseif (i.eq.3) then
                phase_ramp_shift = .false.
            
            else
                goto 110
            
            endif
              
        endif

        write(6,*) 'Enter the number of passes to perform for QEP calculation:'
        write(*,*) 'Warning: using only a single pass is usually NOT sufficient.'
        call get_input("Number of Monte Carlo calculated", n_qep_passes )
        write(*,*)
    
        write(6,*) 'Enter the starting position of the random number sequence:'
        call get_input("Number of ran1 discarded", nran )
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
        
        
 
	    
end module

