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

module m_probe_scan
    
    use m_precision, only: fp_kind
    
    implicit none
    
    integer(4) :: nxsample, nysample          !number of probe positions
    real(fp_kind), allocatable :: probe_positions(:,:,:) !matrix containing the probe position
    real(fp_kind) :: delx, dely                !stepsize for the probe position in x and y
    real(fp_kind) :: probe_initial_position(3) = [0.0_fp_kind, 0.0_fp_kind, 0.0_fp_kind]
    
    logical :: scan_quarter

    
    
    contains
    
    
    
    subroutine setup_probe_scan(interp_setup)
    
        use global_variables, only: uvw1, uvw2
        
        implicit none

        real(fp_kind) :: fract, origin(3)
        logical,intent(in),optional::interp_setup
        logical::interp_
        
        
        interp_ = .true.
        if (present(interp_setup)) interp_ = interp_setup
        
        write(*,*) '|-------------------------------|'
        write(*,*) '|      Probe scan details       |'
        write(*,*) '|-------------------------------|'
        write(*,*)
                
        write(*,*) 'Warning: changing the following parameters may '
        write(*,*) 'result in incorrect interpolation of STEM images.'
        write(*,*)
        
        call setup_scan_geometry(fract, origin)
        
        call setup_sampling(fract)
    
        call calculate_probe_positions(uvw1, uvw2, origin)
        
        if(interp_) call setup_stem_image_interpolation

    end subroutine

    
    
    subroutine setup_scan_geometry(fract, origin)
    
        use m_user_input, only: get_input
        use global_variables, only: a0, deg, ss, uvw1, uvw2, thetad, izone
        use m_crystallography, only: zone, subuvw, angle, rsd
        
        implicit none

        real(fp_kind),intent(out) :: fract, origin(3)
        
        integer(4) :: ich
        integer(4)   ig1a(3), ig2a(3)
        real(fp_kind) a1, a2, r1(3), r2(3), thetad2
    
        integer :: i_scan_quarter
        
        origin = 0.0_fp_kind
    
        a1 = rsd(uvw1, a0, deg)
        a2 = rsd(uvw2, a0, deg)
        r1 = uvw1
        r2 = uvw2
        thetad2 = thetad 
                            
        fract = 1.0                        

    102 write(6,103) r1, a1, char(143), r2, a2, char(143), thetad2, origin(1), origin(2)
        103 format(' The probe scan vectors are: ', /,                  & 
        &       ' x = ', 3g12.5, ' mag = ', g12.5, 1x, a1, /,         &
        &       ' y = ', 3g12.5, ' mag = ', g12.5, 1x, a1, /,         &
        &       ' The angle between these scan vectors is ', f12.2, ' degrees', /, &
        &       ' The inital (fractional) position is ', g12.5, ', ', g12.5, /, &
        &       ' <1> Accept', /,                  &
        &       ' <2> Change size of x and y.', /,      &
        &       ' <3> Change orientation of x and y.', /,     & 
        &       ' <4> Change the initial position.' )

        call get_input("<1>default<2>size<3>dirn<4>origin", ich)
        write(*,*)
        
        if(ich.eq.2) then             !changing size of x and y vectors
              write(6,111)
    111       format(' Enter fractional increase in x and y')
              call get_input("Enter fractional increase in x and y", fract)

              r1 = fract * r1
              r2 = fract * r2
              
              a1 = fract * a1
              a2 = fract * a2
                        
              goto 102                !Summarise and repeat probe position choices
          
        elseif(ich.eq.3) then
    121       format( /, ' Please enter a new x-scan vector.' )
              write(6,121)
              call get_input("x-scan vector", ig1a, 3)
              call zone(izone, ig1a, ig2a)          !get an orthogonal vector to the zone axis and ig1
              call angle(ig1a, ig2a, ss, thetad2)    !calculate angle between scan vectors
              call subuvw(ig1a, r1, a0, deg, ss)      !calculate the real space scan vector from the 'ig1' given
              call subuvw(ig2a, r2, a0, deg, ss)      !calculate the real space scan length from the 'ig2' given

              goto 102
          
        elseif(ich.eq.4) then
            call place_probe(origin)
            
            goto 102
            
        elseif (ich.ne.1) then
            goto 102
            
        endif

        uvw1 = r1
        uvw2 = r2
        
    end subroutine

    
    
    subroutine unwrap_quarter_image(image, ny, nx, image_out)
    
        use m_precision, only: fp_kind
        
        implicit none
        
        integer :: ny, nx
        real(fp_kind) :: image(ny, nx), image_out(ny, nx)
        
        integer :: i_y, i_x
        
        do i_y = 1, ny
        do i_x = 1, nx
            image_out(i_y, i_x) = image(f(i_y, ny), f(i_x, nx))
        enddo
        enddo
        
        contains
        
        integer function f(i, n)
            implicit none
            
            integer :: i, n
            integer :: shift
            
            shift = (n-1)/2
            
            f = abs(mod(i-1+shift, n) - shift) + 1
            
        end function f
        
    end subroutine
    
    
    
    subroutine place_probe(xyposn)

        use m_precision, only: fp_kind
        use m_user_input, only: get_input
    
        implicit none
    
        real(fp_kind) :: xyposn(3)

        write(*,*) 'Enter the co-ordinates "u v" at which to place the probe:'
        write(*,*) '(as fractions of the unit cell side lengths)'
        call get_input('Initial probe position', xyposn(1:2), 2)
        write(*,*) 
   
        xyposn(3) = 0.0_fp_kind
    
    end subroutine
    
    
    
    subroutine setup_sampling(fract)
    
        use m_user_input, only: get_input
        use m_lens, only: probe_cutoff
        use global_variables, only: a0
        
        implicit none
        
        real(fp_kind) :: fract
        
        real(fp_kind) :: min_step
        integer :: ich
        
        min_step = nyquist_step(probe_cutoff)
        nysample = nyquist_sampling(probe_cutoff, fract*a0(2))
        nxsample = nyquist_sampling(probe_cutoff, fract*a0(1))
        
    97  write(6,98) probe_cutoff, min_step, nxsample, nysample
    98  format(   ' The maximum spatial frequency allowed by the probe is ', f5.2, ' A-1.', /, &
        &         ' The STEM image has a bandwidth limit of twice that frequency. ', /, &
        &         ' This corresponds to a minimum sampling of ', f6.2, ' positions per Angstrom,', /, &
        &         ' which is ', i4, ' x-positions and ', i4, ' y-positions per unit cell.', /, &
        &         ' <1> Accept', /, ' <2> Change')
        call get_input('<1> Accept probe sampling', ich)
        write(*,*)
        
        if(ich.eq.2) then
            write(*,*) 'Enter the number of probe positions in the x direction.'
            call get_input('nxsample', nxsample)
              
            write(*,*) 'Enter the number of probe positions in the y direction.' 
            call get_input('nysample', nysample)
            write(*,*)
              
        elseif (ich.ne.1) then
            goto 97
            
        endif
        
    end subroutine
    
    

    integer function nyquist_sampling(qmax, L) result(n)
    
        implicit none
        
        real(fp_kind) :: qmax, L
        
        n = ceiling(4 * qmax * L)
        
    end function
    
    real(fp_kind) function nyquist_step(qmax) result(step)
    
        implicit none
        
        real(fp_kind) :: qmax
        
        step = 4 * qmax

    end function
    
    
    
    subroutine setup_stem_image_interpolation()
    
        use m_user_input, only: get_input
        use global_variables, only: tiley, tilex, output_nopiy, output_nopix
    
        implicit none

        integer(4) :: out_max
        
        if((nysample.gt.1).and.(nxsample.gt.1)) then
            write(6,99) 
    99      format(' Enter the maximum number of pixels to interpolate the output image to.')
            call get_input('output interpolation max pixels', out_max)
            write(*,*) 
            
            
            write(*,*) 'Enter the tiling in x and y for interpolation output'
            call get_input('output interpolation tilex', tilex)
            call get_input('output interpolation tiley', tiley)
            write(*,*)
            
            if(tilex.lt.1) tilex = 1
            if(tiley.lt.1) tiley = 1

            if(mod(out_max, 2).ne.0) out_max = out_max + 1
            
            if(out_max.lt.max(nxsample*tilex, nysample*tiley)) then
				write(6,105) out_max,max(nxsample*tilex, nysample*tiley)
105			format(' The choice of ',i4,' output pixels means that sampling of the output STEM image',/,&
				& ' would fall below the Nyquist  criterion for your choice of probe parameters.',/,&
				&' The maximum number of output pixels has been increased to ',i4,' to avoid ',/,&
				&' aliasing of the output.',/)
				out_max = max(nxsample*tilex, nysample*tiley)
			endif

            if((nxsample*tilex).ge.(nysample*tiley)) then
                output_nopix = out_max
                output_nopiy = int( float(nysample*tiley)/float(nxsample*tilex)*output_nopix)
            else
                output_nopiy = out_max
                output_nopix = int( float(nxsample*tilex)/float(nysample*tiley)*output_nopiy)
            endif
          
        endif
    
    end subroutine

    
    
    subroutine setup_probe_scan_pacbed()
    
        use global_variables, only: a0, uvw1, uvw2
        use m_lens, only: probe_cutoff
    
        implicit none

        real(fp_kind) :: origin(3)
    
        origin = 0.0_fp_kind
    
        nysample = nyquist_sampling(probe_cutoff, a0(2))
        nxsample = nyquist_sampling(probe_cutoff, a0(1))
        
        call calculate_probe_positions(uvw1, uvw2, origin)
                
    end subroutine

    
    
    subroutine calculate_probe_positions(r1, r2, origin)
    
        use global_variables, only: a0, deg
        use m_crystallography, only: rsd
        
        implicit none
    
        real(fp_kind) :: r1(3), r2(3), origin(3)
        
        integer :: ny, nx
        real(fp_kind) :: sitey(3), sitex(3)
    
        if(allocated(probe_positions)) deallocate(probe_positions)
        allocate(probe_positions(3,nysample,nxsample))

        do ny = 1, nysample
        
            sitey = (ny-1) * r2 / nysample
            
            do nx = 1, nxsample
            
                sitex = (nx-1) * r1 / nxsample
                
                probe_positions(:,ny,nx) = sitey + sitex + origin
                
            enddo
            
        enddo

        delx = rsd(r1, a0, deg) / nxsample      !get the step size in the x scan direction
        dely = rsd(r2, a0, deg) / nysample      !get the step size in the y scan direction

    end subroutine

end module
