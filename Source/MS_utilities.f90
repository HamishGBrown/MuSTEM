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

    !Checks the list of atom labels and will append 2,3,4 etc. if there are duplicates
    !Note that this method is not robust to cases where there are three labels of the form
    !eg. Sr, Sr2, Sr
    subroutine correct_duplicate_atom_labels(substance_atom_types,nt)
    use m_string
    implicit none
    character*10::substance_atom_types(nt)
    integer*4,intent(in)::nt

    integer*4:: unique_list(nt),j,i,ii
    logical::unique


    j=1
    do i=1,nt
        unique = .true.

        do ii=1,i-1
            unique = .not.(trim(adjustl(substance_atom_types(ii)))==trim(adjustl(substance_atom_types(i))))
            if(.not.unique) then
                unique_list(i)=unique_list(ii)+1
                exit
            endif
        enddo

        if(unique) unique_list(i) = 1
    enddo

    do i=1,nt
        if(unique_list(i)>1) substance_atom_types(i)=trim(adjustl(substance_atom_types(i)))//to_string(unique_list(i))
    enddo
    end subroutine

      subroutine set_xtl_global_params()
      use m_precision
      use global_variables
      use m_electron
      use m_user_input
      use m_string, only:is_numeric
      use m_crystallography, only: cryst, zone, subuvw, subhkl, rshkl, angle, trimr, trimi, rsd

      implicit none

      integer(4) iunit,nm_curr

      character*120     xtl_fnam
      character*20      junk
      character*2       junk4
      logical:: contains_kev
      real(fp_kind) junk2(6), junk3, junk5(1:3)

      real(fp_kind) ag1,ag2,ct,proj,temp_vec_length
      real(fp_kind) uvwm1,uvwm2,gg(3),rzone(3)

      integer(4) izone1(3),i,j,jm,icount,itrue

      iunit = 15
100   write(6,101)
101   format(1x, 'Please enter the name of the input .xtl file.')
      call get_input("Input crystal file name", xtl_fnam)
      open(unit=iunit,file=trim(adjustl(xtl_fnam)),status='old',err=998)

      !First check if accelerating voltage is in xtl file
      do i=1,4
          read(iunit,*) junk
      enddo
      !If the fourth line is numeric (ie number of atoms) then
      !the accelerating voltage is in the xtl file
      contains_kev = is_numeric(junk)
      if(.not.contains_keV) then
        write(*,*) char(10),'Please input the probe accelerating voltage in kV'
        call get_input("Probe accelerating voltage (kV)",ekv)
      endif
      rewind(iunit)

      nm = 0
      ! Determine the maximum number of atoms an types to allocate the
      ! pertinent automatic objects stored in the module
      ! global_variables.
      read(iunit,*) junk
      read(iunit,*) junk2(1:6)
      if(contains_kev) read(iunit,*) junk3
      read(iunit,*) nt
      do i = 1, nt
            read(iunit,*) junk4
            read(iunit,*) nm_curr
            if (nm_curr.gt.nm) nm = nm_curr
            do j=1,nm_curr
                  read(iunit,*) junk5(1:3)
            enddo
      enddo
      rewind(iunit)

      allocate(atf(3,nt),tau(3,nt,nm),nat(nt),atomf(13,nt),substance_atom_types(nt),fx(nt),dz(nt))


102   format( a20 )
      read(iunit,102) substance
      write(6,111) substance
111   format(/,1x,a20)

      read(iunit,*) a0(1:3), deg(1:3)
      write(6,121) a0(1:3), char(143), deg(1:3)
121   format(4x,' a = ',f9.4,7x,'b = ',f9.4,6x,'c = ',f9.4,1x,a1,/,&
     &' alpha = ',f9.4,2x,'  beta = ',f9.4,2x,'gamma = ',f9.4,&
     &' degrees',/)

      if(contains_kev) then
      read(iunit,*) ekv
      write(6,131) ekv
131   format(' Incident beam energy = ',f12.3,' keV',/)
      endif

      read(iunit,*) nt
      write(6,141) nt
141   format(' Number of atom atom types = ',i2,/)

      do i = 1, nt
            read(iunit,*) substance_atom_types(i)
            if (ionic) then
                read(iunit,*) nat(i), atf(1:3,i),dz(i)
            else
                read(iunit,*) nat(i), atf(1:3,i)
            endif
            if(nat(i).gt.nm.or.nat(i).lt.0) then
                  write(6,171) i,nat(i)
171               format(' Type ',i3,2x,' nat = ',i5,' this is wrong - EXIT')
                  go to 999
            endif

            write(6,181) i,substance_atom_types(i),atf(1:3,i)
181         format(' Type ',i2,2x,a10,8X,' Z = ',F5.0,/,' Occupancy = ',f6.3,2x,' <us> ** 2 = ',g12.5)
            jm = nat(i)
            write(6,191)
191         format(/,'      x    ',5x,'y',9x,'z',/)

            icount = 0
            do j = 1, jm
                  icount = icount + 1
                  read(iunit,*) tau(1:3,i,j)

                  if (icount.le.100) then
                    write(6,201) icount,tau(1:3,i,j)
201                 format(i5, 3f10.6)
                  elseif (icount.eq.101) then
                    write(*,*) 'Number of atoms exceeds 100.'
                    write(*,*) 'See xtl file for full list.'
                    write(*,*)
                  endif

                  if (any(tau(:,i,j)<0) .or. any(tau(:,i,j)>1)) then
                    write(*,*) 'ERROR: fractional coordinates must be between 0 and 1'
                    write(*,*) 'Program will now halt.'
                    read(*,*)
                    stop
                  endif
                  !if(icount.gt.12) then
                        !icount = 0
                        !write(6,191)
                  !endif
            enddo
      enddo
      call correct_duplicate_atom_labels(substance_atom_types,nt)
      call cryst(a0,deg,ss)   !Establish the triclinic information
      icount = 0

      ! Force to 001 convention
      izone1 = [0, 0, 1]
      ig1 = [1, 0, 0]
      ig2 = [0, 1, 0]

      write(6,231) izone1(1:3)
231   format(' Zone Axis = [ ', 3i4, ' ]')

      write(6,241) ig1(1:3)
241   format(' X-scan reciprocal lattice vector = ( ', 3i4, ' )')

      write(6,251) ig2(1:3)
251   format(' Y-scan reciprocal lattice vector = ( ', 3i4, ' )')

      ag1 = trimi(ig1,ss)
      ag2 = trimi(ig2,ss)

      call angle(ig1,ig2,ss,thetad)

      if(abs(thetad).lt.0.1_fp_kind.or.abs(thetad-180.0_fp_kind).lt.0.1_fp_kind) then
            write(6,261)
261         format(' Scan vectors form a linear array',/,&
     &          ' You need two vectors that are NOT co-linear',/,&
     &          ' Please terminate the program and ammend your',&
     &          ' .xtl file',/)
      endif

      call zone(ig1, ig2, izone)

      itrue = 0
      do i = 1,3
            if(izone1(i).ne.izone(i)) itrue = 1
      enddo

      if(itrue.eq.1) then
            write(6,271) izone1(1:3),izone(1:3)
271         format(' Entered zone = [',3i5,'] is incorrect',/,&
     &          ' New zone defined by scan vectors is [',3i5,' ]')
      endif

      call subuvw(ig1,uvw1,a0,deg,ss)
      call subuvw(ig2,uvw2,a0,deg,ss)
      call subhkl(izone,gg,a0,deg,ss)
      uvwm1 = rsd(uvw1,a0,deg)
      uvwm2 = rsd(uvw2,a0,deg)
      ! find orthogonal reciprocal lattice vectors
      !
      ct = 180.0_fp_kind/(atan(1.0_fp_kind)*4.0_fp_kind)
      ct = cos(thetad/ct)
      if(ct.gt.0.9999_fp_kind) ct = 1.0_fp_kind
      proj = ag2 * ct / ag1
      do i = 1, 3
            orthog(i,1) = real(ig1(i),fp_kind)
            orthog(i,2) = real(ig2(i),fp_kind) - proj * real(ig1(i),fp_kind)
            orthog(i,3) = gg(i)
      enddo
      rzone(1:3) = real(izone(1:3),fp_kind)
      call rshkl(rzone,surfn,a0,deg,ss)
      temp_vec_length = trimr(surfn,ss)
      !
      ! make unit vector surfn for surface normal
      surfn(1:3) = surfn(1:3) / temp_vec_length

      do i = 1,nt
            atomf(1:13,i)=xrayFF(2:14,int(atf(1,i)))
      enddo

      write(6,321)
321   format(' Successfully read all atomic X-ray scattering factors',/,&
     &' from D. Waasmaier & A. Kirfel Acta. Cryst. A51, 416 (1995)',/,  &
     &' Mott formula used for conversion to electron form factors', /)

        i_xtl=1     !set the flag ixtl=1 as in the file has been read
      close(iunit)
      goto 999
998   write(6,*) 'ERROR: Crystal file ', trim(xtl_fnam), ' does not exist!'
      goto 100
999   continue

      return
      end

    subroutine validate_xtl(deg)
        use m_precision
        implicit none
        real(fp_kind),intent(in)::deg(3)

        if (any(abs(deg - 90) .gt. 1e-3)) then
            write(*,*) ' You have set one or more triclinic angles to something other than 90 degrees.'
            write(*,*) ' This program only works when the inputted crystal structure is given in terms'
            write(*,*) ' of an orthorhombic coordinate system. Please see the manual for more details.'
            write(*,*)
            write(*,*) ' The program will now halt.'
            write(*,*)
            read(*,*)
            stop
        endif

    end subroutine

    subroutine set_volts(nt, atf, nat, atomf, volts, ss)

        use m_precision, only: fp_kind
        use m_electron, only:elsa_ext

        implicit none

        integer(4) :: nt, nat(nt)
        real(fp_kind) :: atf(3,nt), atomf(13,nt), volts, ss(7)

        integer :: m

        real(fp_kind),parameter :: evconv = 47.87801_fp_kind

        if(abs(ss(7)).lt.1.0e-10_fp_kind) then
            volts = 0.0_fp_kind

        else
            volts = 0.0_fp_kind
            do m = 1, nt
                volts = volts + nat(m)*elsa_ext(nt, m, atomf, 0.0_fp_kind)*atf(2,m)
            enddo
            volts = volts * evconv / ss(7)
        endif

        if(abs(volts).gt.50.0_fp_kind) then
            write(6,102) volts
102         format(/,' Inner potential = ',g16.9,' volts. This is unrealistic.',/,&
                   & ' Resetting inner potential to 20 volts.')
            volts = 20.0_fp_kind
        endif

        write(6,101) volts
101     format(/,' Inner potential = ', g16.9,' volts')

    end

    subroutine set_tiling_grid()

        use m_precision, only: fp_kind
        use m_user_input, only: get_input
        use global_variables, only: ifactory, ifactorx, nopiy, nopix, nopiy_ucell, nopix_ucell, &
                                    ig1, ig2, ss, a0, deg, uvw1, uvw2, ak1, deltay, deltax, npixels, normalisation
        use m_lens, only: pw_illum
        use m_potential, only: quick_shift
        use m_crystallography, only: trimi, rsd
        use m_string

        implicit none

        integer(4) ich
        real(fp_kind) sitey(3),sitex(3)
        logical::chosen

        real(fp_kind) max_qy,max_qx,max_mrady,max_mradx

        call command_line_title_box('Unit cell tiling and grid size')

        chosen=.false.
        do while(.not.chosen)
    131 format( ' Enter the integer by which to tile the unit cell in the ', a1,' direction:')
        write(6,131) 'x'
        call get_input('Tile supercell x', ifactorx)
        write(*,131) 'y'
        call get_input('Tile supercell y', ifactory)
        write(*,*)

        if (pw_illum) then
            if (a0(1)*ifactorx .ne. a0(2)*ifactory) then
                write(*,10) a0(1)*ifactorx, char(143), a0(2)*ifactory, char(143)
10              format(1x, 'Warning: the supercell has the non-square dimensions ', f8.2, 1x, a1, ' x ', f8.2, 1x, a1)
                write(*,*) 'If a square grid of pixels is specified, the outputted'
                write(*,*) 'images and diffraction patterns will have the wrong'
                write(*,*) 'aspect ratio. The user may wish to choose the number of'
                write(*,*) 'pixels to avoid this, whilst still taking into consideration'
                write(*,*) 'the other guidelines mentioned in the manual regarding'
                write(*,*) 'choice of grid size.'
                write(*,*)
            endif
        endif

        write(6,111) 'x'
    111 format(' Enter number of pixels in the ', a1, ' direction:')
        call get_input('Number of pixels in x',nopix)

        write(6,111) 'y'
        call get_input('Number of pixels in y',nopiy)
        write(*,*)

        !----------------------------------------------------------------
        !Summarise the choices thus far
        !----------------------------------------------------------------

        nopiy_ucell = nopiy / ifactory;nopix_ucell = nopix / ifactorx

        max_qy = trimi(ig2,ss)*float(nopiy_ucell)/3.0_fp_kind
        max_qx = trimi(ig1,ss)*float(nopix_ucell)/3.0_fp_kind
        max_mrady = atan(max_qy/ak1)*1000.0_fp_kind
        max_mradx = atan(max_qx/ak1)*1000.0_fp_kind


        ! Test if quick shift is possible
        ! (and if unit cell is actually tiled in both directions)
        quick_shift = mod(nopiy,ifactory).eq.0 .and. (mod(nopix,ifactorx)).eq.0 .and. (ifactory.gt.1 .or. ifactorx.gt.1)

            write(6,161) nopix_ucell, nopiy_ucell
            if (.not.quick_shift) write(6,162)
            write(6,163) ifactorx, ifactory, nopix, nopiy, max_qx, char(143),&
                         &max_qy, char(143),max_mradx, max_mrady

161  format(  '                                x               y', /, &
            & ' ---------------------------------------------------------', /, &
            & ' Pixels per unit cell  | ',i13,' | ',i13)
162 format(   '    (Not integer! Quick-shifting not possible for QEP calculations.)')
163 format(   ' Tiling of unit cell   | ',i13,' | ',i13,/,&
            & ' Pixels in super cell  | ',i13,' | ',i13,/,&
            & ' Max scattering vector | ',f9.1,' ', a1, '-1 | ', f9.1, ' ', a1, '-1', /, &
            & '                       | (',f6.1,' mrad) | (', f6.1, ' mrad)', /, &
            & ' ---------------------------------------------------------', /, &
            & ' <1> Continue', /, &
            & ' <2> Change')


        call get_input('<1> Continue <2> Change', ich)
        write(*,*)
        chosen = ich==1
        enddo

        sitey = uvw2 / (float(nopiy)/float(ifactory))
        sitex = uvw1 / (float(nopix)/float(ifactorx))
        deltay = rsd(sitey,a0,deg)
        deltax = rsd(sitex,a0,deg)
        npixels = nopiy*nopix
        normalisation = 1.0_fp_kind/float(npixels)

    end subroutine



    subroutine setup_integration_measurements()

        use m_precision, only: fp_kind
        use m_user_input, only: get_input
        use global_variables
        use m_crystallography, only: trimi
        use m_string
        use output
        use m_multislice

        implicit none
        character*100::dstring
        integer*4::comaindex,i,j
        logical::outputdetectors

        call command_line_title_box(' Diffraction plane detectors')
        write(6,*) 'Enter the number of detectors in the diffraction plane:'
        write(6,*) '(e.g. 3 if you wish to simulate BF/ABF/ADF simultaneously.)'
        write(6,*) "To segment detectors input a comma ',' and then the number"
        write(6,*) "of angular segments (eg. for 4 rings and 4 quadrants input '4,4')."
        write(6,*) 'To output the detectors for inspection end the input with a '
        write(6,*) "question mark ('?'), ie '4,4?'."
        call get_input('Number of detectors', dstring)
        write(*,*)

        !Check if there is a ',' which indicates the user
        !wants segmented detectors
        comaindex = index(dstring,',')
        segments = (comaindex.ne.0)

        !Check if there is a ?, which indicates the user would
        !like to output the detectors
        outputdetectors = (index(dstring,'?').ne.0)
        !Remove the ? so that it doesn't cause problems later
        if (outputdetectors) dstring = dstring(:index(dstring,'?')-1)

        !Read detector string
        if(segments) then
            read(dstring(:comaindex),*) ndet
            read(dstring(comaindex+1:),*) nseg
            ndet = ndet*nseg
			if(nseg>1) then
				write(*,*) 'Please input orientation offset for segmented detectors in degrees'
				call get_input('Segment orientation offset (degrees)', seg_det_offset)
				seg_det_offset = seg_det_offset/180*pi !Convert from degrees to mrad
			endif
		else
			read(dstring,*) ndet
			nseg = 1
		endif

        if(allocated(outer)) deallocate(outer)
        if(allocated(inner)) deallocate(inner)
        allocate(inner(ndet/nseg),outer(ndet/nseg))

        if (ndet.eq.0) return

        call get_cbed_detector_info(ndet/nseg,ak1,outer,inner)

        if(outputdetectors) then
            do i=1,ndet/nseg
            do j=1,nseg
                if(nseg>1) call binary_out_unwrap(nopiy,nopix,make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i)&
                                                  ,2*pi*j/nseg-seg_det_offset,2*pi/nseg),'detector_'//to_string((i-1)*nseg+j))
                if(nseg==1) call binary_out_unwrap(nopiy,nopix,make_detector(nopiy,nopix,ifactory,ifactorx,ss,inner(i),outer(i))&
                                                                                        ,'detector_'//to_string((i-1)*nseg+j))
            enddo
            enddo
        endif

    end





    subroutine setup_specimen_thickness()

        use global_variables, only: thickness, n_cells, a0,nz,zarray,ncells
        use m_user_input, only: get_input
        use m_precision, only: fp_kind
		use m_string, only: read_sequence_string,command_line_title_box,series_prompt

        implicit none

        integer*4::i
        character*120::thickness_string

        call command_line_title_box('Specimen thickness')
        write(6,11)
    11  format( ' Enter the specimen thickness in Angstroms:')
		call Series_prompt('thickness')
		call get_input("Thickness", thickness_string)

		call read_sequence_string(thickness_string,120,nz,minstep=a0(3))
		allocate(zarray(nz),ncells(nz))
		call read_sequence_string(thickness_string,120,nz,zarray,a0(3))
		ncells = nint(zarray/a0(3))
		zarray = ncells*a0(3)

        call read_sequence_string(thickness_string,120,nz,minstep=a0(3))
        allocate(zarray(nz),ncells(nz))
        call read_sequence_string(thickness_string,120,nz,zarray,a0(3))
        ncells = nint(zarray/a0(3))
        zarray = ncells*a0(3)

        thickness = zarray(nz)

        n_cells = nint(thickness/a0(3))
        thickness = n_cells*a0(3)
        do i=1,nz
            write(6, 15) ncells(i), zarray(i), char(143)
        enddo
    15  format(' This corresponds to ', i5, ' unit cells with a total thickness of ', f6.1, ' ', a1, '.')

        write(*,*)

    end subroutine

    subroutine fourD_STEM_options(fourdSTEM,nopiyout,nopixout,nopiy,nopix)
        use m_user_input
        use m_string
        integer*4,intent(in)::nopiy,nopix
        integer*4,intent(out)::nopiyout,nopixout
        logical,intent(out)::fourdSTEM

        integer*4::idum


        call command_line_title_box('4D-STEM options')

        write(*,*) 'Output diffraction patterns for each scan position?'
        write(*,*) '<1> Yes',char(10),' <2> Output cropped diffraction patterns (saves memory)',char(10),' <3> No'&
					,'<4> No, but output cropped PACBED pattern',char(10)
		call get_input('<1> Diffraction pattern for each probe position',idum)
		fourDSTEM = (idum == 1).or.(idum ==2)
        if(idum==2.or.idum==4) then
            write(*,*) 'Please input number of y pixels in diffration pattern output'
            call get_input('diffraction pattern y pixels',nopiyout)

            write(*,*) 'Please input number of x pixels in diffration pattern output'
            call get_input('diffraction pattern x pixels',nopixout)
        else
            nopiyout = nopiy/3*2
            nopixout = nopix/3*2
        endif
    end subroutine


    !--------------------------------------------------------------------------------------
    subroutine make_bwl_mat()
    !   subroutine make_bwl_mat forms a real matrix with entries 1 for
    !   reciprocal space vectors with magnitude inside the cutoff for
    !   band-width limitting, and 0 outside.  Multiplying a reciprocal
    !   space array with the resultant bwl_mat *is* applying band-width
    !   limitting.

    use global_variables
    use m_precision
    use m_crystallography, only: rsd

    implicit none

    real(fp_kind) return_array(nopiy, nopix)

    integer(4) middlex, middley
    real(fp_kind) rad, xstep, ystep,a1,a2
    real(fp_kind) xstep2, ystep2, xpos, ypos
    integer(4) i, j

    if(allocated(bwl_mat)) deallocate(bwl_mat)
    allocate(bwl_mat(nopiy, nopix))
    middlex = (nopix+1)/2
    middley = (nopiy+1)/2
    a1 = rsd(uvw1,a0,deg)
    a2 = rsd(uvw2,a0,deg)
    ystep = 1.0_fp_kind/(float(ifactory)*a2)
    xstep = 1.0_fp_kind/(float(ifactorx)*a1)

    ystep2 = ystep*ystep
    xstep2 = xstep*xstep

    bwl_rad = min( ystep*nopiy, xstep*nopix )
    bwl_rad = 1.0_fp_kind*bwl_rad/3.0_fp_kind

    do i=1, nopiy
          do j=1, nopix

                ypos = float(middley - i)
                xpos = float(middlex - j)

                rad = sqrt(ypos*ypos*ystep2 + xpos*xpos*xstep2)

                if( rad .gt. bwl_rad ) then
                      bwl_mat(i,j) = 0.0_fp_kind
                else
                      bwl_mat(i,j) = 1.0_fp_kind
                endif

          enddo
    enddo

    return_array=cshift(bwl_mat,SHIFT=-middlex,DIM=2)
    bwl_mat=cshift(return_array,SHIFT=-middley,DIM=1)

    return
    end


    !-------------------------------------------------------------------------------------
    ! subroutine to make the 1D factorisation array to speed up the array shift
    !
    subroutine make_shift_oned(shift, dim_shift, coord)

        use m_precision, only: fp_kind
        use global_variables, only: pi
        use m_crystallography, only: make_g_vec_array

        implicit none

        integer(4) :: i,half_shift
        integer(4),intent(in) :: dim_shift
        real(fp_kind),intent(in) :: coord
        real(fp_kind) :: q
        complex(fp_kind),dimension(dim_shift),intent(out) :: shift

        half_shift = (dim_shift-1)/2-1
        do i = 1, dim_shift
            q = float(mod( i+half_shift, dim_shift) - half_shift -1)
            shift(i) = exp(cmplx(0.0_fp_kind, -2.0_fp_kind*pi*q*coord))
        enddo

    end

    !----------------------------------------------------------------
    !Subroutine to do the fastest shift array from the precomputed
    !1D factorised shift arrays
    subroutine phase_shift_array(input, output_, shifty, shiftx)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix
        use FFTW3
        use output

        implicit none

        complex(fp_kind),intent(in) :: shifty(nopiy), shiftx(nopix),input(nopiy,nopix)
        complex(fp_kind),dimension(nopiy,nopix),intent(out) :: output_

        output_ = input*spread(shiftx,dim=1,ncopies=nopiy)*spread(shifty,dim=2,ncopies = nopix)
        call inplace_ifft(nopiy, nopix, output_,norm=.true.)

    end

    subroutine STEM_options(STEM,ionization,PACBED,istem,double_channeling)
        use m_string
        use m_user_input
        logical,intent(out)::STEM,ionization,PACBED,istem
		logical,intent(inout)::double_channeling
		logical::dc_init

        integer*4::i

        STEM= .false.
        ionization = .false.
        PACBED = .false.
		istem = .false.
		dc_init = double_channeling
		double_channeling=.false.

        i=-1

        do while(i.ne.0)
            call command_line_title_box('STEM modes')
            write(*,*)'muSTEM offers a number of options for STEM modes,'
            write(*,*)' you can choose to do any number of them simultaneously'
            write(*,*)'-----------------------------------------------------------'
            write(*,*)'Option                                       | Included(y/n)'
            write(*,*)'-----------------------------------------------------------'
            write(*,*)'<1> Conventional STEM (ADF,ABF,BF etc.)      | ',logical_to_yn(STEM)
            write(*,*)'<2> Ionization based STEM (EELS and EDX)     | ',logical_to_yn(ionization)
            write(*,*)'<3> Diffraction (PACBED and 4D-STEM)         | ',logical_to_yn(PACBED)
			write(*,*)'<4> Imaging STEM (iSTEM)                     | ',logical_to_yn(istem)
			if(dc_init)&
		   &write(*,*)'<5> STEM EELS with double channeling         | ',logical_to_yn(double_channeling)
            write(*,*)'<0> Continue',char(10)
            write(*,*)'-----------------------------------------------------------'

            call get_input('STEM modes',i)
            write(*,*)

            if(i==1) STEM = .not.STEM
            if(i==2) ionization = .not.ionization
            if(i==3) PACBED = .not.PACBED
			if(i==4) istem = .not.istem
			if(dc_init.and.i==5) double_channeling = .not.double_channeling
            if(i==0.and.(.not.any([STEM,ionization,PACBED,double_channeling,istem]))) then
                write(*,*) "You must choose at least one imaging mode to proceed"
                i=-1
            endif

        enddo




    end subroutine
