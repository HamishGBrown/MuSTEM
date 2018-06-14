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

module m_tilt

    use m_precision, only: fp_kind
    use global_variables, only: ifactory,ifactorx
    implicit none

    logical :: tilt_illumination = .false.
    real(fp_kind) :: alpha, beta
    
    contains
    
    
    
    subroutine prompt_tilt
        
        use m_user_input, only: get_input
        
        implicit none
        
        integer :: i_tilt, i_cshift
        
10      write(*,*) 'The illumination can be tilted off the beam axis.'
        write(*,*) '<1> No beam tilt',char(10),' <2> Beam tilt'
        call get_input('<1> No beam tilt <2> Beam tilt', i_tilt)
        write(*,*)
        
        if (i_tilt.eq.1) then
            tilt_illumination = .false.
            
            
        elseif (i_tilt.eq.2) then
            tilt_illumination = .true.
			
            call setup_tilt
        
        else
            goto 10
            
        endif

20      write(*,*) 'The specimen can be tilted off the specimen axis.'
        write(*,*) '<1> No specimen tilt',char(10),' <2> Specimen tilt'
        call get_input('<1> No specimen tilt <2> Specimen tilt', i_tilt)
        write(*,*)
        
        if (i_tilt.eq.2) call setup_specimen_tilt
        
    end subroutine
    
    
    
    subroutine setup_tilt
    
        use m_user_input, only: get_input
        use global_variables, only: bvec,ak1,ig1,ig2,ss
        use m_crystallography, only:trimi
        
        implicit none
        
        real(fp_kind)::tilt_theta,tilt_phi,tilt_phi_,tilt_theta_,bt_y,bt_x
        
        write(*,*) 'Please enter the beam tilt in mrad:',char(10)
        
        call get_input('Beam tilt in mrad', tilt_theta)
        
        write(*,*) char(10),'Please enter the azimuth of the beam tilt, measured from'
        write(*,*) '[100] ("East") clockwise to [010] ("South") in mrad:',char(10)
        call get_input('Beam tilt azimuth in mrad', tilt_phi)

		!Convert angles from mrad -> rad
		tilt_theta = tilt_theta*1e-3_fp_kind
		tilt_phi = tilt_phi*1e-3_fp_kind
		
		!calculate tilt vector components
		!The negative is to keep with the convention of paring x with the 
		!second array dimension and y with the first dimension
		!whilst remaining true to the description of the direction of the 
		!azimuth
		bt_x = ak1*sin(tilt_theta)*cos(tilt_phi)/trimi(ig1,ss)
		bt_y = ak1*sin(tilt_theta)*sin(tilt_phi)/trimi(ig2,ss)
		
		!Round tilt vector to an integer number of pixels
		!This avoids a boundary discontinuity in the tilted beam
		bt_x = float(nint(bt_x * ifactorx))/ifactorx
		bt_y = float(nint(bt_y * ifactory))/ifactory
		
		!Store tilt vector as a vector
		bvec = [bt_y,bt_x,0]
		
		!Recalculate theta and phi for output only
        tilt_phi_ = -atan2(bvec(1),bvec(2))*1e3
        tilt_theta_ = asin(sqrt(sum((bvec*[trimi(ig1,ss),trimi(ig2,ss),0_fp_kind])**2))/ak1)*1e3

        write(6,30) tilt_theta*1e3_fp_kind,tilt_theta_,tilt_phi*1e3_fp_kind,tilt_phi_
                
30     format(/,&
        &1x,'To ensure that the illumination wave function is continuous',/,&
        &1x,'at the boundary of the simulation grid, the beam tilt has',/,&
        &1x,'been rounded from ',f5.1,' to ',f5.1,' mrad and the azimuth has ',/,&
        &1x,'been rounded from ',f7.1,' to ',f7.1,' mrad. If you would prefer ',/,&
        &1x,'that the tilt used in simulation was closer to the value ',/,&
        &1x,'inputted, please consider increasing the dimensions of the',/,&
        &1x,'supercell used in the simulation by increasing the unit cell',/,&
        &1x,'tiling.',/)
        
        
    end subroutine
     
    subroutine setup_specimen_tilt
        
        use global_variables, only: claue,ak1,Kz,ss,ig1,ig2
        use m_crystallography
        use m_user_input, only: get_input
        
        implicit none
        
        real(fp_kind)::tilt_theta,tilt_phi
        
        
    
        write(*,*) 'Please enter the specimen tilt in mrad:',char(10)
        call get_input('Specimen tilt in mrad', tilt_theta)
        
        write(*,*) char(10),'Please enter the azimuth of the specimen tilt, measured from'
        write(*,*) '[100] ("East") clockwise to [010] ("South") in mrad:',char(10)
        call get_input('Specimen tilt azimuth in mrad', tilt_phi)
        write(*,*)

        Kz   = ak1 * cos( tilt_theta*1e-3_fp_kind )
        claue  = ak1 * [ sin(tilt_theta*1e-3_fp_kind)*cos(tilt_phi*1e-3_fp_kind)/trimi(ig1,ss), sin(tilt_theta*1e-3_fp_kind)*sin(tilt_phi*1e-3_fp_kind)/trimi(ig2,ss),0_fp_kind]
    end subroutine      
    
    
    subroutine tilt_wave_function(psi)
        
        use global_variables, only: ifactory, ifactorx, pi,bvec
        use output
        implicit none
        
        complex(fp_kind) :: psi(:,:)
        
        !complex(fp_kind),allocatable:: psi2(:,:)
        integer :: nopiy, nopix
        real(fp_kind) :: shift(3), shift_frac(3),bvec_(3)
        integer :: ny, nx
        
        nopiy = size(psi, 1)
        nopix = size(psi, 2)

        bvec_ = bvec/[nopiy,nopix,1]*[ifactory,ifactorx,0]
        
        !$OMP PARALLEL DO PRIVATE(nx, ny)
        do nx = 1, nopix;do ny = 1, nopiy
            psi(ny,nx) = psi(ny,nx) * exp(cmplx(0.0_fp_kind, 2*pi*dot_product([ny-1, nx-1, 0], bvec_), fp_kind))
        enddo;enddo
        !$OMP END PARALLEL DO
    end subroutine

    
    function fftfreq(n)
    integer*4,intent(in)::n
    integer*4::fftfreq(n),i
    
    if (mod(n,2)==0) then
        !Even case
        do i=1,n
            fftfreq(i) = -n/2+i
        enddo
    else
        !Odd case
        do i=1,n
            fftfreq(i) = -n/2+i-1
        enddo
    endif
    end function
    
end module
