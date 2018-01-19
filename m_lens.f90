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
    real(fp_kind) :: probe_Cs, probe_C5, probe_cutoff
	real(fp_kind) :: probe_fa2, probe_phia2, probe_fa3, probe_phia3, probe_fc3, probe_phic3
    real(fp_kind),allocatable :: probe_df(:)
    integer :: probe_ndf
    
    ! Image forming lens
    real(fp_kind) :: imaging_cs,imaging_c5,imaging_fa2,imaging_phia2,imaging_fa3
	real(fp_kind) :: imaging_phia3,imaging_fc3,imaging_phic3,imaging_cutoff
    real(fp_kind),allocatable :: imaging_df(:)
    integer :: imaging_ndf
    
    
    contains

    
    
    subroutine setup_lens_parameters(string,cs,c5,fa2,phia2,fa3,phia3,fc3,phic3,cutoff)
		
		use global_variables, only: ak1,nopiy,nopix
        use global_variables, only: ak1
        use m_user_input, only: get_input
		use output
		use CUFFT_wrapper
		use m_string
        
        implicit none
    
        integer(4) nflag
		
		character*(*),intent(in) :: string
        real(fp_kind),intent(out) :: cs,c5,fa2,phia2,fa3,phia3,fc3,phic3,cutoff
		
		real(fp_kind) :: stigmatism(2),cs_mm,c5_mm,xyposn(3),scherzer_df
		real(fp_kind),allocatable :: df(:)
		integer(4) ::input,ndf,i
		complex(fp_kind)::probe(nopiy,nopix)
		character(120)::set_defocus
        
        write(*,*) '|------------------------------------------|'
		write(6,5) trim(adjustl(string))
5       format(' |      ',a5,' forming lens parameters       |')
        write(*,*) '|------------------------------------------|'
        write(*,*)
        
		
		fa2=0
		phia2=0
		fa3=0
		phia3=0
		fc3=0
		phic3=0
		Cs=0
		C5=0

        cutoff = set_cutoff(Cs)
        ! df = set_defocus(Cs)
		scherzer_df = -1.0_fp_kind*sign(1.0_fp_kind,Cs)*sqrt(4.0_fp_kind*abs(Cs)/(3.0_fp_kind*ak1))
		write(*,*) 'Enter the defocus in Angstroms:'
        if (abs(Cs).gt.1e-3) write(*,'(1x, a, f7.2, a)') '(The optimal Scherzer defocus for the specified Cs is ', scherzer_df, ' Angstroms)'
        
		call get_input('Defocus', set_defocus)
		call read_sequence_string(set_defocus,120,ndf)
		
		allocate(df(ndf))
		call read_sequence_string(set_defocus,120,ndf,df)
		

		
    
        !if (abs(df).lt.1e-4) df = 0.0_fp_kind
                
     10 write(6,11) atan(cutoff/ak1)*1000.0_fp_kind, char(143), cutoff, char(143), df(1)   
		do i=1,ndf-1
			write(6,12) df(i+1)
		enddo
	   write(6,13) Cs_mm, C5_mm,char(143), fa2,char(232), phia2, char(143), &
	               & fa3, char(232), phia3, char(143),  fc3, char(232), phic3


     11 format(   ' Current simulation parameters:',/,             &
                 &' -----------------------------------------------', /, &
                 &' <1> Aperture cutoff (mrad)      ', t40, g11.4,/, &
                 &'     Aperture cutoff (', a1, '^-1)              ',t40, g11.4,/, &
                 &' <2> Defocus         (', a1, ')         ',t40, g11.4)
     12 format(   '                                        ',t40, g11.4)
	 13 format(   ' <3> Cs              (mm)                ',t40, g11.4,/, &
                 &' <4> C5              (mm)                ',t40, g11.4,/, &
                 &' <5> Two-fold astigmatism',/,&
				 &'     fa2             (', a1, ')          ',t40, g11.4,/, &
				 &'     ',a1,'a2             (rad)        't40, g11.4,/, &
                 &' <6> Three-fold astigmatism',/,&
				 &'     fa3             (', a1, ')          ',t40, g11.4,/, &
				 &'     ',a1,'a3             (rad)        ',t40, g11.4,/, &
                 &' <7> Coma',/,&
				 &'     fc3             (', a1, ')          ',t40, g11.4,/, &
				 &'     ',a1,'c3             (rad)        ',t40, g11.4,/,&
				 &' <8> Output lens contrast transfer function ',t40,/,&
                 &' -----------------------------------------------', /, &
                 &' Select option to change any parameter or <0> to continue')

        call get_input('Change '//to_lower(trim(adjustl(string)))//' forming lens parameters', nflag)
        write(*,*)
        
        if (nflag.eq.1) then
            cutoff = set_cutoff(Cs)
            goto 10
            
        elseif (nflag.eq.2) then
            ! df = set_defocus(Cs)
			scherzer_df = -1.0_fp_kind*sign(1.0_fp_kind,Cs)*sqrt(4.0_fp_kind*abs(Cs)/(3.0_fp_kind*ak1))
			write(*,*) 'Enter the defocus in Angstroms:'
			if (abs(Cs).gt.1e-3) write(*,'(1x, a, f7.2, a)') '(The optimal Scherzer defocus for the specified Cs is ', scherzer_df, ' Angstroms)'
			
			call get_input('Defocus', set_defocus)
			call read_sequence_string(set_defocus,120,ndf)
			
			allocate(df(ndf))
			call read_sequence_string(set_defocus,120,ndf,df)
            goto 10
            
        elseif (nflag.eq.3) then
            Cs_mm = set_cs()
            Cs = Cs_mm*1.0e7_fp_kind
            goto 10
            
        elseif (nflag.eq.4) then
            C5_mm = set_c5()
            Cs = C5_mm*1.0e7_fp_kind
            goto 10
  			
        elseif (nflag.eq.5) then
			stigmatism = set_astigmatism('two-fold astigmatism')
			fa2 = stigmatism(1)
			phia2 = stigmatism(2)
			goto 10
		
		elseif (nflag.eq.6) then
			stigmatism = set_astigmatism('three-fold astigmatism')
			fa3 = stigmatism(1)
			phia3 = stigmatism(2)
			goto 10
			
		elseif (nflag.eq.7) then
			stigmatism = set_astigmatism('coma')
			fc3 = stigmatism(1) 
			phic3 = stigmatism(2)    
			goto 10
		elseif(nflag.eq.8) then
			xyposn = [0,0,0]
			probe = make_ctf(xyposn,cs,c5,fa2,phia2,fa3,phia3,fc3,phic3,df(1),cutoff)
			call binary_out_unwrap(nopiy,nopix,atan2(imag(probe),real(probe))*abs(probe)**2,&
			                      &trim(adjustl(output_prefix)) //'_'//trim(adjustl(string))//'_forming_lens_ctf_phase')
			call ifft2(nopiy,nopix,probe,nopiy,probe,nopiy)
			call binary_out_unwrap(nopiy,nopix,abs(probe)**2,trim(adjustl(output_prefix))//'_'//trim(adjustl(string))&
			                                                      &//'_forming_lens_ctf_real_space_intensity')

		    goto 10	          
        elseif (nflag.ne.0) then
            goto 10
            
        endif
		
		if(trim(adjustl(string))=='Image') then
			imaging_ndf=ndf
			allocate(imaging_df(imaging_ndf))
			call read_sequence_string(set_defocus,120,imaging_ndf,imaging_df)
		else
			probe_ndf=ndf
			allocate(probe_df(probe_ndf))
			call read_sequence_string(set_defocus,120,probe_ndf,probe_df)
		endif
    end  subroutine    
    subroutine setup_probe_lens_parameters()
		
		use global_variables, only: ak1,nopiy,nopix
        use global_variables, only: ak1
        use m_user_input, only: get_input
		use output
		use CUFFT_wrapper
        
        implicit none
    
        integer(4) nflag

        real(fp_kind) :: probe_Cs_mm, probe_C5_mm,stigmatism(2)
		integer(4) ::input 
		complex(fp_kind)::probe(nopiy,nopix)
        call setup_lens_parameters('Probe',probe_cs,probe_c5,probe_fa2,probe_phia2,probe_fa3,&
								&probe_phia3,probe_fc3,probe_phic3,probe_cutoff)
		return

    end  subroutine



    subroutine setup_imaging_lens_parameters()
    
        use global_variables, only: ak1
        use m_user_input, only: get_input
        
        implicit none

        integer(4) nflag
        
        real(fp_kind) :: imaging_Cs_mm, imaging_C5_mm
		imaging=.true.
		call setup_lens_parameters('Image',imaging_cs,imaging_c5,imaging_fa2,imaging_phia2,imaging_fa3,&
								&imaging_phia3,imaging_fc3,imaging_phic3,imaging_cutoff)
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
    
	function set_astigmatism(tag)
		use m_user_input, only: get_input
        
        implicit none
    
		character(*),intent(in)::tag
        real(fp_kind) set_astigmatism(2)

        write(6,*)' Enter the coefficient of '//trim(adjustl(tag))//' in '//char(143)//':' 
        call get_input(trim(adjustl(tag))//' coefficient', set_astigmatism(1))
        write(*,*)
        write(6,*)' Enter the azimuthal orientation of '//trim(adjustl(tag))//' astigmatism ('//char(237)//') in radians:' 
        call get_input(trim(adjustl(tag))//' azimuth', set_astigmatism(2))
        write(*,*)
	end function			 

    subroutine make_lens_ctf(lens_ctf, temp_lens_defocus)

        use global_variables, only: nopiy, nopix, ig1, ig2, ifactory, ifactorx, pi, ss
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: lens_ctf(nopiy,nopix)
        real(fp_kind) :: temp_lens_defocus, xyposn(3)

        integer(4) :: ny, nx, shifty, shiftx
        real(fp_kind) :: kx(3), ky(3), kr(3), akr, m1, m2
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
    
		xyposn = [0,0,0]
		lens_ctf = make_ctf(xyposn,imaging_cs,imaging_c5,imaging_fa2,imaging_phia2,imaging_fa3,imaging_phia3,imaging_fc3,imaging_phic3,temp_lens_defocus,imaging_cutoff)
    end  subroutine

    function make_ctf(xyposn,cs,c5,fa2,phia2,fa3,phia3,fc3,phic3,df,cutoff) result(ctf)

        use global_variables, only: nopiy, nopix, ig1, ig2, ifactory, ifactorx, pi, ss
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: ctf(nopiy,nopix)
        real(fp_kind) :: xyposn(3),cs,c5,fa2,phia2,fa3,phia3,fc3,phic3,df,cutoff

        integer(4) :: ny, nx, shifty, shiftx
        real(fp_kind) :: kx(3), ky(3), kr(3), akr, m1, m2,phi
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
    
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1

        ig1_temp = float(ig1)/float(ifactorx)
        ig2_temp = float(ig2)/float(ifactory)

        !$OMP PARALLEL DO PRIVATE(nx, m1, kx, ny, m2, ky, kr, akr,phi)
        do nx = 1, nopix
              m1 = float(mod( nx+shiftx, nopix) - shiftx -1)
              kx = m1 * ig1_temp

              do ny = 1, nopiy
                    m2 = float(mod( ny+shifty, nopiy) - shifty -1)
                    ky = m2 * ig2_temp
                    kr = kx + ky 
                    akr = trimr(kr,ss)

                    if (akr.le.cutoff) then
						  phi = atan2(kr(2),kr(1))
                          ctf(ny,nx) = exp(cmplx(0.0_fp_kind, -1*aberrate(akr,df,Cs,C5,phi,fa2,phia2,&
																		& fa3, phia3, fc3, phic3), fp_kind))
						  ctf(ny,nx) = ctf(ny,nx) * exp(cmplx(0.0_fp_kind, -2*pi*dot_product(kr, xyposn), fp_kind))
                    else
                          ctf(ny,nx) = 0.0_fp_kind

                    endif
               enddo
        enddo
        !$OMP END PARALLEL DO
    
    end  function


    
    function aberrate(q, df, Cs, C5,phi, fa2, phia2, fa3, phia3, fc3, phic3)
        ! Calculate the aberration function chi(q) for the transfer function \exp[-i chi(q)]
    
        use global_variables, only: pi, ak1
        
        implicit none
    
        real(fp_kind) :: aberrate
        
									  
		
        real(fp_kind),intent(in) :: q,phi, df, Cs, C5,fa2,phia2,fa3,phia3,fc3,phic3
		optional:: phi,Cs, C5, fa2, phia2, fa3, phia3, fc3, phic3
		
        real(fp_kind) :: q2,phi_
		
		phi_=0
		if(present(phi)) phi_=phi
        
        q2 = q**2
        aberrate = pi*df*q2/ak1
		if(present(cs)) aberrate = aberrate + pi*0.5_fp_kind*Cs*q2*q2/(ak1**3) 
		if(present(c5)) aberrate = aberrate + pi*((1.0_fp_kind/3.0_fp_kind)*C5) * (q2*q2*q2)/(ak1**5)
		if(present(fa2).and.present(phia2)) aberrate = aberrate + pi*fa2*q2*sin(2*(phi_-phia2))/ak1
		if(present(fa3).and.present(phia3)) aberrate = aberrate + 2*pi/3*fa3*q2*abs(q)*sin(3*(phi_-phia3))/(ak1**2)
		if(present(fc3).and.present(phic3)) aberrate = aberrate + 2*pi/3*fc3*q2*abs(q)*sin((phi_-phic3))/(ak1**2)
    
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
        real(fp_kind) :: norm,phi

		psi = make_ctf(xyposn,probe_cs,probe_c5,probe_fa2,probe_phia2,probe_fa3,probe_phia3,probe_fc3,probe_phic3,defocus,probe_cutoff)
        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
        
        norm = sqrt(sum(abs(psi)**2))

        psi = psi/norm
    end subroutine

end module
