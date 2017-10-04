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


      subroutine set_xtl_global_params()
      use m_precision
      use global_variables
      use m_xray_factors
      use m_user_input
      use m_crystallography, only: cryst, zone, subuvw, subhkl, rshkl, angle, trimr, trimi, rsd
      
      implicit none
      
      integer(4) iunit,nm_curr
      
      character*120     xtl_fnam
      character*20      junk
      character*2       junk4
      real(fp_kind) junk2(6), junk3, junk5(1:3)
      
      integer(4) nb
      real(fp_kind) ag1,ag2,ct,proj,temp_vec_length
      real(fp_kind) uvwm1,uvwm2,gg(3),rzone(3)
              
      integer(4) izone1(3),i,j,jm,icount,itrue,zindex,k

      iunit = 15
100   write(6,101)
101   format(1x, 'Please enter the name of the input .xtl file.')
      call get_input("Input crystal file name", xtl_fnam)
      open(unit=iunit,file=xtl_fnam,status='old',err=998)
      
      nm = 0
      ! Determine the maximum number of atoms an types to allocate the
      ! pertinent automatic objects stored in the module
      ! global_variables.	
      read(iunit,*) junk
      read(iunit,*) junk2(1:6)
      read(iunit,*) junk3
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
        
      if (i_xtl.eq.1) then
            deallocate( atf,tau,atomf,nat,substance_atom_types,fx)
      endif
      allocate(atf(3,nt),tau(3,nt,nm),nat(nt),atomf(13,nt),substance_atom_types(nt),fx(nt))


102   format( a20 )
      read(iunit,102) substance
      write(6,111) substance
111   format(/,1x,a20)

      read(iunit,*) a0(1:3), deg(1:3)
      write(6,121) a0(1:3), char(143), deg(1:3)
121   format(4x,' a = ',f9.4,7x,'b = ',f9.4,6x,'c = ',f9.4,1x,a1,/,&
     &' alpha = ',f9.4,2x,'  beta = ',f9.4,2x,'gamma = ',f9.4,&
     &' degrees',/)

      read(iunit,*) ekv
      write(6,131) ekv
131   format(' Incident beam energy = ',f12.3,' keV',/)

      read(iunit,*) nt
      write(6,141) nt
141   format(' Number of atom atom types = ',i2,/)

      do i = 1, nt
            read(iunit,*) substance_atom_types(i)
            read(iunit,*) nat(i), atf(1:3,i)
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
                    pause
                    stop
                  endif
                  !if(icount.gt.12) then
                        !icount = 0
                        !write(6,191)
                  !endif
            enddo
      enddo
      call cryst(a0,deg,ss)   !Establish the triclinic information
      read(iunit,*) junk
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
            orthog(i,1) = dble(ig1(i))
            orthog(i,2) = dble(ig2(i)) - proj * dble(ig1(i))
            orthog(i,3) = gg(i)
      enddo
      rzone(1:3) = dble(izone(1:3))
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
      
      goto 999
998   write(6,*) 'ERROR: Crystal file ', trim(xtl_fnam), ' does not exist!'
      goto 100     
999   continue
      
      return
      end
      
      
      
    subroutine validate_xtl()
    
        use global_variables, only: deg, ig1, ig2
        
        implicit none
        
        if (any(abs(deg - 90) .gt. 1e-3)) then
            write(*,*) ' You have set one or more triclinic angles to something other than 90 degrees.'
            write(*,*) ' This program only works when the inputted crystal structure is given in terms'
            write(*,*) ' of an orthorhombic coordinate system. Please see the manual for more details.'
            write(*,*)
            write(*,*) ' The program will now halt.'
            write(*,*)
            pause
            stop
            
        endif
                
    end subroutine
      
      
      
    subroutine set_volts(nt, atf, nat, atomf, volts, ss)

        use m_precision, only: fp_kind
        use m_elsa, only:elsa_ext

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
        use m_qep, only: quick_shift
        use m_crystallography, only: trimi, rsd
        
        implicit none
        
        integer(4) ich, i, j
        real(fp_kind) sitey(3),sitex(3),site(3)
        
        real(fp_kind) max_qy,max_qx,max_mrady,max_mradx
        real(fp_kind) k

        write(*,*) '|------------------------------------------|'
        write(*,*) '|      Unit cell tiling and grid size      |'
        write(*,*) '|------------------------------------------|'
        write(*,*)
        
    !    write(6,100)
    !100 format(/,&
    !           &' You will now be asked how the unit cell should be tiled in each direction, ', /, &
    !           &' as well as the number of pixels in each direction. The number of pixels ', /, &
    !           &' should preferably be a product of powers of small primes, in order to speed ', / &
    !           &' up the calculations which employ the Fast Fourier Transform. If possible also ', / &
    !           &' try to ensure that the number of pixels per unit cell is an integer, ', /, &
    !           &' as this may further speed up calculations.                        ', /, &
    !           )
                   
    131 format( ' Enter the integer by which to tile the unit cell in the ', a1,' direction:')
    110 write(6,131) 'x'
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
        
        nopiy_ucell = nopiy / ifactory
        nopix_ucell = nopix / ifactorx

        max_qy = trimi(ig2,ss)*float(nopiy_ucell)/3.0_fp_kind
        max_qx = trimi(ig1,ss)*float(nopix_ucell)/3.0_fp_kind
        max_mrady = atan(max_qy/ak1)*1000.0_fp_kind
        max_mradx = atan(max_qx/ak1)*1000.0_fp_kind
        
    165 continue
        if (mod(nopiy,ifactory).eq.0 .and. mod(nopix, ifactorx).eq.0) then
            write(6,161) nopix_ucell, nopiy_ucell, &
                         ifactorx, ifactory, &
                         nopix, nopiy, &
                         max_qx, char(143), max_qy, char(143), &
                         max_mradx, max_mrady
        161 format(  '                                x               y', /, &
            &        ' ---------------------------------------------------------', /, &
            &        ' Pixels per unit cell  | ',i13,' | ',i13,/,&
            &        ' Tiling of unit cell   | ',i13,' | ',i13,/,&
            &        ' Pixels in super cell  | ',i13,' | ',i13,/,&
            &        ' Max scattering vector | ',f9.1,' ', a1, '-1 | ', f9.1, ' ', a1, '-1', /, &
            &        '                       | (',f6.1,' mrad) | (', f6.1, ' mrad)', /, &
            &        ' ---------------------------------------------------------', /, &
            &        ' <1> Continue', /, &
            &        ' <2> Change')
            
        else
            write(6,162) float(nopiy)/ifactory, float(nopix)/ifactorx, &
                         ifactorx, ifactory, &
                         nopix, nopiy, &
                         max_qx, char(143), max_qy, char(143), &
                         max_mradx, max_mrady
        162 format(  '                                x               y', /, &
            &        ' ---------------------------------------------------------', /, &
            &        ' Pixels per unit cell  | ',f13.1,' | ',f13.1,/,&
            &        '    (Not integer! Quick-shifting not possible for QEP calculations.)', /, &
            &        ' Tiling of unit cell   | ',i13,' | ',i13,/,&
            &        ' Pixels in super cell  | ',i13,' | ',i13,/,&
            &        ' Max scattering vector | ',f9.1,' ', a1, '-1 | ', f9.1, ' ', a1, '-1', /, &
            &        '                       | (',f6.1,' mrad) | (', f6.1, ' mrad)', /, &
            &        ' ---------------------------------------------------------', /, &
            &        ' <1> Continue', /, &
            &        ' <2> Change')
            
        endif    
        
        call get_input('<1> Continue <2> Change', ich)
        write(*,*)
        
        if (ich.eq.2) then
            goto 110

        elseif (ich.ne.1) then
            goto 165

        endif

        ! Test if quick shift is possible
        ! (and if unit cell is actually tiled in both directions)
        if (mod(nopiy,ifactory).eq.0 .and. (mod(nopix,ifactorx)).eq.0 .and. ifactory.gt.1 .and. ifactorx.gt.1) then 
            quick_shift = .true.

        else
            quick_shift = .false.

        endif

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
        use global_variables, only: ak1, inner, outer, ndet
        
        implicit none
    
        write(*,*) '|-----------------------------------------|'
	    write(*,*) '|      Diffraction plane detectors        |'
	    write(*,*) '|-----------------------------------------|'
        write(*,*)
    
        write(6,*) 'Enter the number of detectors in the diffraction plane:'
        write(6,*) '(e.g. 3 if you wish to simulate BF/ABF/ADF simultaneously.)'
        call get_input('Number of detectors', ndet)
        write(*,*)
    
        if(allocated(outer)) deallocate(outer)
        if(allocated(inner)) deallocate(inner)
        allocate(inner(ndet),outer(ndet))
        
        if (ndet.eq.0) return
        
        call get_cbed_detector_info(outer,inner,ndet,ak1)
    
    end
    

    
    subroutine get_cbed_detector_info(outer,inner,ndet,k)
    
        use m_precision, only: fp_kind
        use m_user_input, only: get_input
        use m_string, only: to_string
        
        implicit none
    
        integer(4)   ndet,ichoice,i,mrad
        real(fp_kind)      outer(ndet),inner(ndet),k,dummy
    
        real(fp_kind) :: inner_mrad(ndet), outer_mrad(ndet)
        
        write(*,*) 'Select a method for choosing inner and outer angles:'
        write(*,*) '<1> Manual'
        write(*,*) '<2> Automatic'
    
        call get_input("manual detector <1> auto <2>",ichoice)
        write(*,*)
    
        write(*,*) 'Select how angles will be specified:'
        write(*,*) '<1> mrad'
        write(*,*) '<2> inverse Angstroms'
        call get_input("<1> mrad <2> inv A", mrad)
        write(*,*)
    
        if(ichoice.eq.1) then
              do i = 1, ndet
                    write(*,*) 'Detector ', to_string(i)
                    write(*,*) 'Inner angle:'
                    call get_input("inner",dummy)
                    if(mrad.eq.1) then 
                          inner(i) = k*tan(dummy/1000.0_fp_kind)
                          inner_mrad(i) = dummy
                    else
                          inner(i) = dummy
                          inner_mrad(i) = 1000*atan(dummy/k)
                    endif
                
                    write(*,*) "Outer angle:"
                    call get_input("outer",dummy)
                    if(mrad.eq.1) then 
                          outer(i) = k*tan(dummy/1000.0_fp_kind)
                          outer_mrad(i) = dummy
                    else
                          outer(i) = dummy
                          outer_mrad(i) = 1000*atan(dummy/k)
                    endif
                    
                    write(*,*)
                    
              enddo
        else
              write(*,*) "Initial inner angle:"
              call get_input("initial inner angle", dummy)
              if(mrad.eq.1) then
                    inner(1) = k*tan(dummy/1000.0_fp_kind)
                    inner_mrad(1) = dummy
              else
                    inner(1) = dummy
                    inner_mrad(1) = 1000*atan(dummy/k)
              endif
              
              write(*,*) "Initial outer angle:"
              call get_input("initial outer angle", dummy)
              if(mrad.eq.1) then
                    outer(1) = k*tan(dummy/1000.0_fp_kind)
                    outer_mrad(1) = dummy
              else
                    outer(1) = dummy
                    outer_mrad(1) = 1000*atan(dummy/k)
              endif
              
              write(*,*) "Increment (both angles incremented by this amount):"
              call get_input("increment", dummy)
              do i = 1, ndet-1
                    if(mrad.eq.1) then
                          inner(1+i) = inner(i) + k*tan(i * dummy/1000.0_fp_kind)
                          outer(1+i) = outer(i) + k*tan(i * dummy/1000.0_fp_kind)
                          
                          inner_mrad(1+i) = inner_mrad(i) + i * dummy
                          outer_mrad(1+i) = outer_mrad(i) + i * dummy
                    else
                          inner(1+i) = inner(i) + i * dummy
                          outer(1+i) = outer(i) + i * dummy
                          
                          inner_mrad(1+i) = inner_mrad(i) + i * 1000 * atan(dummy/k)
                          outer_mrad(1+i) = outer_mrad(i) + i * 1000 * atan(dummy/k)
                         
                    endif
              enddo
              
              write(*,*)
        endif
        
        write(*,*) 'Summary of diffraction plane detectors:'
        write(*,*)
        
        write(*,*) '         inner    outer'
        write(*,*) '  --------------------------------'
        do i = 1, ndet
            write(*,50) i, inner_mrad(i), outer_mrad(i)
50          format(1x, i5, ' | ', f6.2, ' | ', f6.2, '   (mrad)')
            write(*,55) inner(i), outer(i), char(143)
55          format(1x, 5x, ' | ', f6.2, ' | ', f6.2, '   (', a1, '^-1)')
            write(*,60) 
60          format(1x, 5x, ' | ', 6x, ' | ')
        enddo
        
        write(*,*)
        
    end subroutine


    
    subroutine setup_specimen_thickness()
    
        use global_variables, only: thickness, n_cells, a0
        use m_user_input, only: get_input
        use m_precision, only: fp_kind

        implicit none
      
        write(*,*) '|-------------------------------|'
	    write(*,*) '|      Specimen thickness       |'
	    write(*,*) '|-------------------------------|'
        write(*,*)
      
    10  write(6,11)
    11  format( ' Enter the specimen thickness in Angstroms:')
        call get_input("Thickness", thickness)
        write(*,*)
    
        n_cells = nint(thickness/a0(3))
        thickness = n_cells*a0(3)
    
    15  format(' This corresponds to ', i5, ' unit cells with a total thickness of ', f6.1, ' ', a1, '.')
        write(6, 15) n_cells, thickness, char(143)
        write(*,*)

    end subroutine
    

    
    !--------------------------------------------------------------------------------------
    subroutine make_bwl_mat()
    !	subroutine make_bwl_mat forms a real matrix with entries 1 for
    !	reciprocal space vectors with magnitude inside the cutoff for
    !	band-width limitting, and 0 outside.  Multiplying a reciprocal
    !	space array with the resultant bwl_mat *is* applying band-width
    !	limitting.

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
    

    !--------------------------------------------------------------------------------
    !     subroutine make mask
    !     makes a mask based on an input inne and outer q radius in inv A
    !--------------------------------------------------------------------------------

    subroutine make_detector_mask(inner_rad,outer_rad,mask)
    use global_variables
    use m_precision
    use m_crystallography, only: trimr
        
    implicit none
    integer(4)   m1,m2,ny,nx,shifty,shiftx
    
    real(fp_kind)      akr,inner_rad,outer_rad,kx(3),ky(3),kr(3)
    real(fp_kind)      ig1_temp(3),ig2_temp(3)
    real(fp_kind)      mask(nopiy,nopix)
       
    mask=0.0_fp_kind

    shifty = (nopiy-1)/2-1
    shiftx = (nopix-1)/2-1
    ig1_temp=float(ig1)/float(ifactorx)
    ig2_temp=float(ig2)/float(ifactory)

    do nx=1,nopix
          m1 = float(mod( nx+shiftx, nopix) - shiftx -1)
          kx = m1 * ig1_temp
          !$OMP PARALLEL PRIVATE(m2,ky,kr,akr), SHARED(mask) 
          do ny=1,nopiy
                m2 = float(mod( ny+shifty, nopiy) - shifty -1)
                ky = m2 * ig2_temp
                kr = kx + ky 
                akr = trimr(kr,ss)
                if ( (akr.le.outer_rad).and.(akr.ge.inner_rad)) mask(ny,nx)=1.0_fp_kind
           enddo
           !$OMP END PARALLEL
    enddo

    end


    
    !-------------------------------------------------------------------------------------
    ! subroutine to make the 1D factorisation array to speed up the array shift
    !
    subroutine make_shift_oned(shift, dim_shift, coord)
    
        use m_precision, only: fp_kind
        use global_variables, only: pi
    
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

    
    
    !make the shift array from the precomputed 1d factorisations
    subroutine make_shift_array(shift_array, shifty, shiftx)
    
        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix
        
        implicit none
    
        integer(4) :: i, j
        complex(fp_kind),intent(in) :: shifty(nopiy), shiftx(nopix)
        complex(fp_kind),dimension(nopiy,nopix),intent(out) :: shift_array
    
    
        do j = 1, nopix
            do i = 1, nopiy
                shift_array(i,j) = shiftx(j)*shifty(i)
            enddo
        enddo
    
    end
    
    
    
    !----------------------------------------------------------------
    !Subroutine to do the fastest shift array from the precomputed
    !1D factorised shift arrays
    subroutine phase_shift_array(input, output, shifty, shiftx)
    
        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix
        use cufft_wrapper, only: ifft2
        
        implicit none
    
        complex(fp_kind),intent(in) :: shifty(nopiy), shiftx(nopix)
        complex(fp_kind),dimension(nopiy,nopix),intent(in) :: input
        complex(fp_kind),dimension(nopiy,nopix),intent(out) :: output
        complex(fp_kind),dimension(nopiy,nopix) :: shift_array
    
    
        call make_shift_array(shift_array, shifty, shiftx)
        output = input*shift_array
        call ifft2(nopiy, nopix, output, nopiy, output, nopiy)
    
    end
    

    