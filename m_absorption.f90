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

module m_absorption
    
    use m_precision, only: fp_kind
	use quadpack
    
    implicit none
    
    ! Flags
    logical :: complex_absorption
    logical :: include_absorption
    
    ! Local effective TDS calculation
    real(fp_kind) :: theta_min, theta_max
    real(fp_kind),allocatable :: adfbr(:,:)
    
    ! Absorptive potential
    real(fp_kind),allocatable :: tdsbr(:,:)
    complex(fp_kind),allocatable :: fz_abs(:,:,:)
    complex(fp_kind),allocatable :: transf_absorptive(:,:,:)  
    
    ! Domain for get_inelffs() calculation
    real(fp_kind) :: kstep(3), max_dist, min_step, xstep, ystep
    integer :: max_int
    
    ! Variables used in get_inelffs() integration
    real(fp_kind) :: qxunit(3), qyunit(3), qzunit(3), qz(3)
    real(fp_kind) :: qkz(3)
    real(8) :: ct, st
    real(8) :: dwf
    real(fp_kind) :: tp2
    real(fp_kind) :: g(3)
    real(fp_kind) :: delta_kstep
    integer :: i_species

    real(8),parameter :: abs_tol = 1e-10,rel_tol = 5e-6
	
	!real(fp_kind),parameter :: abs_tol = 1e-10,rel_tol = 5e-6
    
    
    contains


    
    subroutine prompt_include_absorption
    
        use m_user_input, only: get_input
        
        implicit none
        
        integer :: choice
        
        write(*,*) 'You can choose to omit the absorptive potential'
        write(*,*) 'so that only a thermally smeared elastic potential'
        write(*,*) 'is used (typically just for checking numerics).'
206     write(*,*) '<1> Include absorption'
        write(*,*) '<2> Do not include absorption'
        call get_input("<1> Absorption <2> No absorption", choice) 
        write(*,*)
        
		include_absorption = .not.(choice==2)
                
    end subroutine

    
    
    subroutine calculate_absorption_mu
        
        use global_variables, only: nt, ss, atf, nat, ak, relm, orthog
        
        implicit none
        
        real(8),parameter :: pi = 4.0d0*atan(1.0d0)
		!real(fp_kind),parameter :: pi = 4.0d0*atan(1.0d0)
        
        call setup_absorptive_array
        
        if (allocated(tdsbr)) deallocate(tdsbr)
        allocate(tdsbr(max_int,nt))
          
        if(include_absorption) then
            call get_inelffs(tdsbr, max_int, kstep, ss, atf, nat, ak, relm, orthog, 0.0_8, pi)
          
        else
            tdsbr = 0.0_fp_kind
       
        endif
    
    end subroutine
    
    
    
    subroutine setup_absorptive_array
    
        use global_variables, only: ig1, ig2, ifactorx, ifactory, nopix, nopiy, a0
        implicit none
    
        xstep = sqrt( dfloat( dot_product(ig1,ig1) ) ) / float(ifactorx)
        ystep = sqrt( dfloat( dot_product(ig2,ig2) ) ) / float(ifactory)
        min_step = min( xstep, ystep )

        if (abs(min_step-xstep).lt.1.0e-10_fp_kind) then
           kstep(1:3) = float( ig1(1:3) ) / float(ifactorx)
           delta_kstep = kstep(1)/a0(1)
           
        else
           kstep(1:3) = float( ig2(1:3) ) / float(ifactory)
           delta_kstep = kstep(2)/a0(2)
           
        endif

        max_dist = sqrt( (xstep/2.0_fp_kind * nopix)**2 + (ystep/2.0_fp_kind * nopiy)**2 )
        max_int = nint( max_dist / min_step ) + 2
    
    end subroutine
    
    
    
    subroutine setup_local_diffraction_plane_geometry
    
        use m_user_input, only: get_input
        use global_variables, only: ak1, inner, outer
        
        implicit none
        
        write(*,*) '|---------------------------------------|'
        write(*,*) '|      Diffraction plane detector       |'
        write(*,*) '|---------------------------------------|'
        write(*,*)
        
        write(*,*) 'Inner angle (mrad):'
        call get_input('Diffraction plane detector inner angle (mrad)', theta_min)
        
        write(*,*) 'Outer angle (mrad):'
        call get_input('Diffraction plane detector outer angle (mrad)', theta_max)
    
        write(*,*)
        
        theta_min = theta_min / 1000
        theta_max = theta_max / 1000
        
        if(allocated(inner)) deallocate(inner)
        if(allocated(outer)) deallocate(outer)
        allocate(inner(1))
        allocate(outer(1))
        
        inner(1) = ak1*tan(theta_min) 
        outer(1) = ak1*tan(theta_max) 

    end subroutine
    


    subroutine calculate_local_adf_mu
    
        use global_variables, only: nt, atf, ss, atf, nat, ak, relm, orthog
        
        implicit none
    
        call setup_absorptive_array
        
        if (allocated(adfbr)) deallocate(adfbr)
        allocate( adfbr(max_int,nt) )
    
        if(sum(atf(3,:)).eq.0.0_fp_kind) then
            adfbr = 0.0_fp_kind
          
        else
            call get_inelffs(adfbr, max_int, kstep, ss, atf, nat, ak, relm, orthog, real(theta_min, kind=8), real(theta_max, kind=8))
            
        endif

    end subroutine
    
    
    
    subroutine get_inelffs(tdsbr_t, max_int, kstep, ss, atf, nat, ak, relm, orthog, thmin, thmax)
		!use quadpack
        ! Calculate TDS form factors
        ! These are actually M_g = mu_g/4pi
		
		
        use global_variables, only: nt
        use m_crystallography, only: trimi, trimr
		use m_string
		use output, only: output_prefix,timing
        
        implicit none

        integer*4 :: max_int
        integer*4 :: nat(nt)
        real(fp_kind) :: tdsbr_t(max_int,nt)
        real(fp_kind) :: atf(3,nt), orthog(3,3), kstep(3), ss(7), relm, ak
        real(8) :: thmin, thmax
    
        real(fp_kind) :: zx(3), zy(3), zz(3)
        real(fp_kind) :: g2, mu_0
		real(8)::abserr

        integer :: i, ipa,m,ier
    
        ! Double precision variables for accuracy
        real(8),parameter :: pi = 4.0d0*atan(1.0d0)
        real(8) :: sum1
		!real(8) :: sum2,sum3,diff,mu(1000)
		real(fp_kind) :: t1, delta
    

        ! Make sure orthogonal coordinate system has been previously set up
        ! (should have been done when reading in xtl file)
        do i = 1, 3
            if (all(orthog(:,i).eq.0.0_fp_kind)) then
                write(*,*) 'Orthogonal coordinate system not set up.'
                write(*,*) 'Calculation cannot proceed.'
                pause
                stop
            endif
        enddo

        zx = orthog(:,1)
        zy = orthog(:,2)
        zz = orthog(:,3)
    
        call validate_coordinate_system(zx, zy, zz)
    
        ! Calculate unit vectors
        qxunit = zx / trimr(zx,ss)
        qyunit = zy / trimr(zy,ss)
        qzunit = zz / trimr(zz,ss)
    
        tp2 = 2.0_fp_kind * pi * pi

        ! Set up q vector z component
        qkz = ak * qzunit
        ! Note that the wavenumber corrected for refraction is used.
        ! This is appropriate since we are considering inelastic
        ! scattering (i.e. phonon excitation) within the mean inner
        ! potential of the specimen.

        ! Accumulator for mu_0
        mu_0 = 0.0_fp_kind
		t1 = secnds(0.0)

        ! Loop over reciprocal lattice vectors g
        do ipa = 1, max_int

            ! Status
#ifdef GPU            
        131 format(a1,' Calculating absorptive scattering factor ', i4, ' of ', i4, '...' )
            write(6,131, advance='no') achar(13),ipa, max_int
#else
        131 format(1h+,' Calculating absorptive scattering factor ', i4, ' of ', i4, '...' )
            write(6,131) ipa, max_int
#endif
            !flush(unit=6)

            g = (ipa-1) * kstep

            g2 = trimr(g,ss)**2

            ! Loop over atomic species
            do i_species = 1, nt
        
                ! Set global dwf which is used in tds_calc_phi()
                dwf = exp( - tp2 * g2 * atf(3, i_species) )
												
                ! Integrate over theta and phi
				call qag ( tds_calc_theta, thmin, thmax, abs_tol, rel_tol, 1, sum1, abserr, m, ier )
				
                ! Fractional occupancy
                sum1 = sum1 * atf(2,i_species)

                ! Store result (not yet in correct units)
                tdsbr_t(ipa,i_species) = sum1
            
                ! Accumulate mu_0
                if (all(g.eq.0.0_fp_kind)) then
                    mu_0 = mu_0 + sum1*nat(i_species) 
                endif

            enddo

        enddo

        write(*,*)
        write(*,*)

		
		if(timing) then
			delta = secnds(t1)
			open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
			write(9834, '(a, g, a, /)') 'Calculation of absorptive form factors took ', delta, 'seconds.'
			close(9834)
		endif
		
        ! 1/Vc factor: 1/A for 2D Fourier series and 1/dz to get average potential
        ! Note that this is the unit cell volume and is corrected to be 
        ! the slice volume elsewhere in the program
        tdsbr_t = tdsbr_t / ss(7)
    
        ! Multiply by gamma^2, since there are two scattering factors
        ! to be corrected
        tdsbr_t = tdsbr_t * relm**2

        ! M'_g = mu'_g/(4*pi)
        tdsbr_t = tdsbr_t / (4 * pi)

        ! Display kinematic mean free path 
        if (mu_0.eq.0.0_fp_kind) then
            write(*,*) 'Kinematic mean free path for TDS is effectively infinite.'

        else
            write(*,211) 1/mu_0, char(143)
    211     format(' Kinematic mean free path for TDS = ', f10.2, 1x, a1)

        endif
    
        write(*,*)

    end subroutine

    
    
    subroutine validate_coordinate_system(zx, zy, zz)
    
        use m_crystallography, only: angler
        use global_variables, only: ss
        
        implicit none
        
        real(fp_kind) :: zx(3), zy(3), zz(3)
        
        real(fp_kind) :: degree
        
        call angler(zx,zy,ss,degree)

        if(abs(degree-90.0_fp_kind).gt.0.1_fp_kind) then
            write(6,821) degree
      821   format(' The angle between the x and y directions is ', f12.3, ' degrees.', /, &
                  &' These vectors are not mutually perpendicular.')
            write(*,*) 'Calculation cannot proceed.'
            pause
            stop
       
        endif

        call angler(zx,zz,ss,degree)

        if(abs(degree-90.0_fp_kind).gt.0.1_fp_kind) then
           write(6,831) degree
      831    format(' The angle between the x and z directions is ', f12.3, ' degrees.', /, &
                  &' These vectors are not mutually perpendicular.')
            write(*,*) 'Calculation cannot proceed.'
            pause
            stop
        endif

        call angler(zy,zz,ss,degree)
        if(abs(degree-90.0_fp_kind).gt.0.1_fp_kind) then
           write(6,841)
      841    format(' The angle between the y and z directions is ', f12.3, ' degrees.', /, &
                  &' These vectors are not mutually perpendicular.')
            write(*,*) 'Calculation cannot proceed.'
            pause
            stop
        endif

    end subroutine
    
    

    function tds_calc_theta(theta)
    
        use global_variables, only: ak
        
        implicit none

        real(8) :: tds_calc_theta
		!real(fp_kind) :: tds_calc_theta
    
        real(8) :: theta
    
        real(8) :: total_one, total_two,abserr
        real(8),parameter :: pi = 4.0d0*atan(1.0d0)
        real(8),parameter :: twopi = 8.0d0*atan(1.0d0)

		integer*4:: m,ier

        st = sin(theta)
        ct = cos(theta)
    
        qz = ak * qzunit * ct - qkz

		call qag ( tds_calc_phi, 0.0_8, pi, abs_tol, rel_tol, 1, total_one, abserr, m, ier )
		call qag ( tds_calc_phi, pi,      2*pi, abs_tol, rel_tol, 1, total_two, abserr, m, ier )
    
        tds_calc_theta = (total_one + total_two) * st
    
    end function tds_calc_theta


    
    function tds_calc_phi(phi)
  
        use global_variables, only: ss, nt, atomf, atf, ak,ionic,dz
        use m_crystallography, only: trimr        
        use m_elsa, only: elsa_ext,peng_ionic_ff
        
        implicit none
    
        real(8) :: tds_calc_phi
    
        real(8) :: phi
    
        real(fp_kind) :: qx(3), qy(3), q1(3), q2(3)
        real(8) :: q1_sq, q2_sq
        real(8) :: f1, f2
        real(8) :: dwf1

        qx = ak * qxunit * st * cos(phi)
        qy = ak * qyunit * st * sin(phi)
    
        q1 = qx + qy + qz
        q2 = qx + qy + qz - g 
    
        q1_sq = trimr(q1,ss)**2
    
        q2_sq = trimr(q2,ss)**2
		
		if(ionic) then
			f1 = Peng_ionic_FF(real(q1_sq/4,kind = fp_kind),nint(atf(1,i_species)),dZ(i_species))
			f2 = Peng_ionic_FF(real(q2_sq/4,kind = fp_kind),nint(atf(1,i_species)),dZ(i_species))
        else
			f1 = elsa_ext(nt,i_species,dble(atomf),q1_sq/4)
			f2 = elsa_ext(nt,i_species,dble(atomf),q2_sq/4)
		endif
        
    
        dwf1 = exp( - tp2 * atf(3,i_species) * (q1_sq + q2_sq) )
  
        tds_calc_phi = f1*f2*(dwf-dwf1)

    end function tds_calc_phi



    end module
