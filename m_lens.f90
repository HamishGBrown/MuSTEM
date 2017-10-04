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


module m_lens
    
    use m_precision, only: fp_kind

    implicit none

    logical :: pw_illum, cb_illum
    logical :: imaging = .false.
    
    ! Probe forming lens
    real(fp_kind) :: probe_Cs, probe_C5, probe_df, probe_cutoff
    
    ! Image forming lens
    real(fp_kind) :: imaging_Cs, imaging_C5, imaging_df, imaging_cutoff
    
    ! Defocus series
    integer :: n_df
    real(fp_kind) :: delta_df
    real(fp_kind),allocatable :: defoci(:)
    
    
    contains

    
    
    subroutine setup_probe_lens_parameters()
    
        use global_variables, only: ak1
        use m_user_input, only: get_input
        
        implicit none
    
        integer(4) nflag

        real(fp_kind) :: probe_Cs_mm, probe_C5_mm
        
        write(*,*) '|------------------------------------------|'
        write(*,*) '|      Probe forming lens parameters       |'
        write(*,*) '|------------------------------------------|'
        write(*,*)
        
        probe_Cs_mm = set_cs()
        probe_Cs = probe_Cs_mm*1.0e7_fp_kind
    
        probe_C5_mm = set_c5()
        probe_C5 = probe_C5_mm*1.0e7_fp_kind
    
        probe_cutoff = set_cutoff(probe_Cs)
        probe_df = set_defocus(probe_Cs)
    
        if (abs(probe_df).lt.1e-4) probe_df = 0.0_fp_kind
                
     10 write(6,11) atan(probe_cutoff/ak1)*1000.0_fp_kind, char(143), probe_cutoff, char(143), probe_df, probe_Cs_mm, probe_C5_mm

     11 format(   ' Current simulation parameters:',/             &
                 &' -----------------------------------------------', /, &
                 &' <1> Aperture cutoff (mrad)      ', t40, g11.4,/, &
                 &'     Aperture cutoff (', a1, '^-1)              ',t40, g11.4,/, &
                 &' <2> Defocus         (', a1, ')         ',t40, g11.4,/, &
                 &' <3> Cs              (mm)                ',t40, g11.4,/, &
                 &' <4> C5              (mm)                ',t40, g11.4,/, &
                 &' -----------------------------------------------', /, &
                 &' Select option to change any parameter or <0> to continue')

        call get_input('Change probe lens parameters', nflag)
        write(*,*)
        
        if (nflag.eq.1) then
            probe_cutoff = set_cutoff(probe_Cs)
            goto 10
            
        elseif (nflag.eq.2) then
            probe_df = set_defocus(probe_Cs)
            goto 10
            
        elseif (nflag.eq.3) then
            probe_Cs_mm = set_cs()
            probe_Cs = probe_Cs_mm*1.0e7_fp_kind
            goto 10
            
        elseif (nflag.eq.4) then
            probe_C5_mm = set_c5()
            probe_Cs = probe_C5_mm*1.0e7_fp_kind
            goto 10
            
        elseif (nflag.ne.0) then
            goto 10
            
        endif

    end  subroutine



    subroutine setup_imaging_lens_parameters()
    
        use global_variables, only: ak1
        use m_user_input, only: get_input
        
        implicit none

        integer(4) nflag
        
        real(fp_kind) :: imaging_Cs_mm, imaging_C5_mm

        write(*,*) '|------------------------------------------|'
        write(*,*) '|      Image forming lens parameters       |'
        write(*,*) '|------------------------------------------|'
        write(*,*)
        
        imaging=.true.

        imaging_Cs_mm = set_cs()
        imaging_Cs = imaging_Cs_mm*1.0e7_fp_kind

        imaging_C5_mm = set_c5()
        imaging_C5 = imaging_C5_mm*1.0e7_fp_kind
        
        imaging_cutoff = set_cutoff(imaging_Cs)
        imaging_df = set_defocus(imaging_Cs)

        if (abs(imaging_df).lt.1e-4) imaging_df = 0.0_fp_kind
        
     10 write(6,11) atan(imaging_cutoff/ak1)*1000.0_fp_kind, char(143), imaging_cutoff, char(143), imaging_df, imaging_Cs_mm, imaging_C5_mm

     11 format(   ' Current simulation parameters:',/             &
                 &' -----------------------------------------------', /, &
                 &' <1> Aperture cutoff (mrad)      ', t40, g11.4,/, &
                 &'     Aperture cutoff (', a1, '^-1)              ',t40, g11.4,/, &
                 &' <2> Defocus         (', a1, ')         ',t40, g11.4,/, &
                 &' <3> Cs              (mm)                ',t40, g11.4,/, &
                 &' <4> C5              (mm)                ',t40, g11.4,/, &
                 &' -----------------------------------------------', /, &
                 &' Select option to change any parameter or <0> to continue')

        call get_input('Change imaging lens parameters', nflag)
        write(*,*)
        
        if (nflag.eq.1) then
            imaging_cutoff = set_cutoff(imaging_Cs)
            goto 10

        elseif (nflag.eq.2) then
            imaging_df = set_defocus(imaging_Cs)
            if (abs(imaging_df).lt.1e-4) imaging_df = 0.0_fp_kind
            goto 10

        elseif (nflag.eq.3) then
            imaging_Cs_mm = set_cs()
            imaging_Cs = imaging_Cs_mm*1.0e7_fp_kind
            goto 10

        elseif (nflag.eq.4) then
            imaging_C5_mm = set_c5()
            imaging_C5 = imaging_C5_mm*1.0e7_fp_kind
            goto 10
 
        elseif (nflag.ne.0) then
            goto 10

        endif

    end subroutine



    function set_cutoff(cstemp)
    
        use global_variables, only: ak1
        use m_user_input, only: get_input
    
        implicit none

        real(fp_kind) set_cutoff, cstemp, mrad_opt, mrad_cutoff

        set_cutoff = 1.51_fp_kind*abs(cstemp)**(-0.25_fp_kind)*ak1**(0.75_fp_kind)
        mrad_opt = atan(set_cutoff/ak1)*1000.0_fp_kind
        if (mrad_opt .gt. 1500) mrad_opt = 1500
        
100     write(*,*) 'Enter the aperture cutoff in mrad:' 
        if (abs(cstemp).gt.1e-3) then
            write(*,'(1x, a, f7.2, a)') '(The optimal cutoff for the specified Cs is ', mrad_opt, ' mrad)'
        endif
        
        call get_input('aperture cutoff', mrad_cutoff)
        write(*,*)

        if(mrad_cutoff < 0.0_fp_kind) then
            goto 100
        endif

        set_cutoff = ak1*tan(mrad_cutoff/1000.0_fp_kind)
        
    end function
    
    
    
    function set_defocus(cstemp)
    
        use global_variables, only: ak1
        use m_user_input, only: get_input
        
        implicit none
    
        real(fp_kind) set_defocus
        real(fp_kind) scherzer_df
        real(fp_kind) cstemp
        integer(4)   def_choice

        scherzer_df = -1.0_fp_kind*sign(1.0_fp_kind,Cstemp)*sqrt(4.0_fp_kind*abs(Cstemp)/(3.0_fp_kind*ak1))
        
        write(*,*) 'Enter the defocus in Angstroms:'
        if (abs(cstemp).gt.1e-3) write(*,'(1x, a, f7.2, a)') '(The optimal Scherzer defocus for the specified Cs is ', scherzer_df, ' Angstroms)'
        call get_input('Defocus', set_defocus)
        write(*,*)
                
    end function
    
    
    
    function set_cs()
    
        use m_user_input, only: get_input
        
        implicit none
    
        real(fp_kind) set_cs

        write(6,*)' Enter the coefficient of spherical aberration (Cs) in mm:' 
        call get_input('Cs coefficient', set_cs)
        write(*,*)

    end function
    
    
    
    function set_c5()
    
        use m_user_input, only: get_input
        
        implicit none
    
        real(fp_kind) set_c5

        write(6,*)' Enter the coefficient of fifth order spherical aberration (C5) in mm:' 
        call get_input('C5 coefficient', set_c5)
        write(*,*)

    end function
    
    
    
    subroutine setup_defocus_series
    
        use m_user_input, only: get_input
        use global_variables, only: thickness
        
        implicit none
    
        integer(4) iok, i_prompt
        integer :: i_df
        
        write(*,*) '|---------------------------|'
        write(*,*) '|      Defocus series       |'
        write(*,*) '|---------------------------|'
        write(*,*)
        
     1  write(*,*) '<1> Setup a defocus series'
        write(*,*) '<2> Continue with single defocus'
        call get_input('<1> Setup defocus series', i_prompt)
        write(*,*)
        
        if (i_prompt.eq.1) then
        
         10 write(*,*) 'Enter the number of defocus steps you wish to use:' 
            call get_input("Number of defocus steps", n_df)
            write(*,*)
    
            if (n_df.eq.1) then
                delta_df = 0.0d0
        
            elseif (n_df.gt.1) then
                delta_df = thickness/n_df
        
             49 format(' The default defocus step size is ', f12.4,' Angstrom ', /, &
                      &' <1> Accept ', /, &
                      &' <2> Change ')
             50 write(6,49) delta_df
                call get_input("Accept defocus step size <1> yes <2> no", iok)
                write(*,*)
                
                if (iok.eq.1) then
                    continue
            
                elseif (iok.eq.2) then
                    write(*,*) 'Enter the MAGNITUDE of the defocus step size:'
                    call get_input("Magnitude of defocus step size", delta_df)
                    write(*,*)
                    
                else
                    goto 50
            
                endif
                
            else
                goto 10
        
            endif
            
        elseif (i_prompt.eq.2) then
            n_df = 1
            delta_df = 0.0d0
        
        else
            goto 1
        
        endif
        
        allocate(defoci(n_df))
        
        do i_df = 1, n_df
            defoci(i_df) = probe_df - (i_df-1)*delta_df
        enddo
          
    end subroutine
    


    subroutine make_ctf(ctf, temp_lens_defocus)

        use global_variables, only: nopiy, nopix, ig1, ig2, ifactory, ifactorx, pi, ss
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: ctf(nopiy,nopix)
        real(fp_kind) :: temp_lens_defocus, xyposn(3)

        integer(4) :: ny, nx, shifty, shiftx
        real(fp_kind) :: kx(3), ky(3), kr(3), akr, m1, m2
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
    
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1

        ig1_temp = float(ig1)/float(ifactorx)
        ig2_temp = float(ig2)/float(ifactory)

        !$OMP PARALLEL DO PRIVATE(nx, m1, kx, ny, m2, ky, kr, akr)
        do nx = 1, nopix
              m1 = float(mod( nx+shiftx, nopix) - shiftx -1)
              kx = m1 * ig1_temp

              do ny = 1, nopiy
                    m2 = float(mod( ny+shifty, nopiy) - shifty -1)
                    ky = m2 * ig2_temp
                    kr = kx + ky 
                    akr = trimr(kr,ss)

                    if (akr.le.imaging_cutoff) then
                          ctf(ny,nx) = exp(cmplx(0.0_fp_kind, -1*aberrate(akr,temp_lens_defocus,imaging_Cs,imaging_C5), fp_kind))

                    else
                          ctf(ny,nx) = 0.0_fp_kind

                    endif
               enddo
        enddo
        !$OMP END PARALLEL DO
    
    end  subroutine


    
    function aberrate(q, df, Cs, C5)
        ! Calculate the aberration function chi(q) for the transfer function \exp[-i chi(q)]
    
        use global_variables, only: pi, ak1
        
        implicit none
    
        real(fp_kind) :: aberrate
        
        real(fp_kind) :: q, df, Cs, C5
        
        real(fp_kind) :: q2
        
        q2 = q**2
        aberrate = pi*df*q2/ak1 + pi*0.5_fp_kind*Cs*q2*q2/(ak1**3) + pi*((1.0_fp_kind/3.0_fp_kind)*C5) * (q2*q2*q2)/(ak1**5)
    
    end function
    
    
    
    subroutine make_stem_wfn(psi, defocus, xyposn)
    
        use global_variables, only: nopiy, nopix, ifactory, ifactorx, ig1, ig2, ss, pi
        use CUFFT_wrapper, only: ifft2
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: psi(nopiy,nopix)
        real(fp_kind) :: defocus, xyposn(3)

        integer(4) :: m1,m2,ny,nx,shifty,shiftx
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
        real(fp_kind) :: kx(3), ky(3), kr(3), akr
        real(fp_kind) :: norm

        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1

        ig1_temp = float(ig1)/float(ifactorx)
        ig2_temp = float(ig2)/float(ifactory)

        !$OMP PARALLEL DO PRIVATE(nx, m1, kx, ny, m2, ky, kr, akr)
        do nx = 1, nopix
            m1 = float(mod( nx+shiftx, nopix) - shiftx -1)
            kx = m1 * ig1_temp

            do ny = 1, nopiy
                m2 = float(mod( ny+shifty, nopiy) - shifty -1)
                ky = m2 * ig2_temp
                kr = kx + ky 
                akr = trimr(kr,ss)

                if(akr.le.probe_cutoff) then
                    psi(ny,nx) = exp(cmplx(0.0_fp_kind, -1*aberrate(akr,defocus,probe_Cs,probe_C5), fp_kind))
                    psi(ny,nx) = psi(ny,nx) * exp(cmplx(0.0_fp_kind, -2*pi*dot_product(kr, xyposn), fp_kind))

                else
                    psi(ny,nx) = 0.0_fp_kind

                endif
            enddo
        enddo
        !$OMP END PARALLEL DO
    
        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
        
        norm = sqrt(sum(abs(psi)**2))
        psi = psi/norm

    end subroutine

    
    
end module
