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

	module m_qep
    
    use m_precision, only: fp_kind
	use m_numerical_tools, only: ran1,gasdev
    
    implicit none
    
    integer(4) :: n_qep_grates
    integer(4) :: n_qep_passes
    integer(4) :: nran ! Start of random number sequence
    
    complex(fp_kind),allocatable :: qep_grates(:,:,:,:)
    
    logical :: phase_ramp_shift
    logical(4) :: quick_shift

    complex(fp_kind),allocatable,dimension(:,:) :: shift_arrayx, shift_arrayy                                                         
    

    
    contains
    
    
    
    subroutine setup_qep_parameters
    
        use global_variables, only: ifactory, ifactorx
        use m_user_input, only: get_input
        use m_string, only: to_string
        
        implicit none
    
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
        101 format( ,' The choice of tiling and grid size does not permit quick shifting ', /, &
                    &' of the precalculated transmission functions. Shifting using ', /, &
                    &' phase ramps can be performed but is time consuming. You may wish ', /, &
                    &' to go back and calculate more distinct transmission functions, or ', /, &
                    &' simply proceed without using phase ramp shifting. ', /)
                
        110 write(6,111)
        111 format(  ' <1> Go back and choose a larger number', /, &
                    &' <2> Proceed with phase ramp shifting', /, &
                    &' <3> Proceed without phase ramp shifting', /, &                
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
        
        
        
        integer function seed_rng() result(idum)
            
            use m_precision, only: fp_kind
            
            implicit none

            integer :: i
            real(fp_kind) :: random
            
	        idum = -1
            
	        do i = 1, nran
		        random = ran1(idum)
	        enddo
            
        end function seed_rng
        
        
        
    subroutine displace(tauin, tauout, urms, a0, idum)

        use m_precision
		
    
	    implicit none

	    real(fp_kind) tauin(3), tauout(3), a0(3)
	    real(fp_kind) urms
	    integer(4) idum
	    real(fp_kind) xran, yran, zran
	    !real(fp_kind) gasdev

	    xran = gasdev(idum)
	    yran = gasdev(idum)
	    !zran = gasdev(idum)

	    xran = xran * urms / a0(1)
	    yran = yran * urms / a0(2)
	    !zran = zran * urms / a0(3)

	    tauout(1) = tauin(1) + xran
	    tauout(2) = tauin(2) + yran
	    tauout(3) = tauin(3)
    !	tauout(3) = tauin(3) + zran
	
    end subroutine displace
	
  

    end module
    
    
    
    
    
    
    
    
