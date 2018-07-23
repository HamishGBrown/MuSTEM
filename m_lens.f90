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
    
    type aberration_coefficient
		character(len=2)::Haider
		character(len=10)::Krivanek
		character(len=17)::Description
		real(fp_kind)::amplitude,angle
		integer*4::m,n
    contains
		procedure :: initialise
		procedure :: set
    end type aberration_coefficient
    
    ! Probe forming lens
    type(aberration_coefficient)::probe_aberrations(14),imaging_aberrations(14)
    real(fp_kind),allocatable :: probe_df(:),imaging_df(:)
    real(fp_kind)::probe_cutoff,imaging_cutoff
    real(fp_kind)::probe_apodisation,imaging_apodisation
    
    integer :: probe_ndf,imaging_ndf

    
	
    contains
	subroutine initialise(ab,Krivanek,Haider,description,amplitude,angle,n,m)
		class(aberration_coefficient),intent(out)::ab
		character(len=*),intent(in)::Krivanek
		character(len=2),intent(in)::Haider
		character(len=*),intent(in)::Description
		real(fp_kind),intent(in)::amplitude,angle
		integer*4,intent(in)::m,n
		
		ab%Haider = Haider
		ab%Krivanek = trim(adjustl(Krivanek))
		ab%Description = trim(adjustl(Description))
		ab%amplitude = amplitude
		ab%angle = angle
		ab%m = m
		ab%n =n
	end subroutine
		
	subroutine set(ab)
        use m_user_input, only: get_input
		class(aberration_coefficient),intent(inout)::ab
	
		character(:),allocatable::tag
		
		
		tag = trim(adjustl(ab%description))
        write(6,*)' Enter the amplitude of '//tag//' in '//char(143)//':' 
        call get_input(tag//' coefficient', ab%amplitude)
        write(*,*)
		if(ab%m>0) then
			write(6,*)' Enter the azimuthal orientation of '//tag//' ('//char(237)//') in radians:' 
			call get_input(tag//' azimuth', ab%angle)
		endif
        write(*,*)
	end subroutine
    
    function chi(ab,q,phi,df)
        ! Calculate the aberration function chi(q) for the transfer function \exp[-i chi(q)]
    
        use global_variables, only: pi, ak1
        
        implicit none
    
        real(fp_kind) :: chi
		type(aberration_coefficient),intent(in)::ab(14)
        real(fp_kind),intent(in) :: q,phi, df
        
		integer*4:: naberrations,i
        real(fp_kind) :: q2
		        
        q2 = q**2
        chi = pi*df*q2/ak1
		
		do i=1,14
			chi = chi + 2*pi*ak1*(sqrt(q2)/ak1)**(ab(i)%n+1)*ab(i)%amplitude/(ab(i)%n+1)*cos(ab(i)%m*(phi-ab(i)%angle))
		enddo
        
    end function
	
    subroutine setup_lens_parameters(string,aberrations,cutoff)
		
		use global_variables, only: ak1,nopiy,nopix
        use m_user_input, only: get_input
		use output
		use CUFFT_wrapper
		use m_string
        
        implicit none
    
        integer(4) nflag
		
		character*(*),intent(in) :: string
        type(aberration_coefficient),intent(out)::aberrations(14)
        real(fp_kind),intent(out):: cutoff
		
		real(fp_kind) :: stigmatism(2),cs_mm,c5_mm,xyposn(3),scherzer_df,cs,apodisation
		real(fp_kind),allocatable :: df(:)
		integer(4) ::input,ndf,i
		complex(fp_kind)::probe(nopiy,nopix)
		character(120)::set_defocus
		
        write(*,*) '|------------------------------------------|'
		write(6,5) trim(adjustl(string))
5       format('  |      ',a5,' forming lens parameters       |')
        write(*,*) '|------------------------------------------|'
        write(*,*)
        
		call aberrations(1 )%initialise('C12'      ,'A1','2-Fold astig.    ',0.0_fp_kind,0.0_fp_kind,1,2)
		call aberrations(2 )%initialise('C21'      ,'B2','Axial coma       ',0.0_fp_kind,0.0_fp_kind,2,1)
		call aberrations(3 )%initialise('C23'      ,'A2','3-Fold astig.    ',0.0_fp_kind,0.0_fp_kind,2,3)
		call aberrations(4 )%initialise('C30 = CS' ,'C3','3rd order spher. ',0.0_fp_kind,0.0_fp_kind,3,0)
		call aberrations(5 )%initialise('C32'      ,'S3','Axial star aber. ',0.0_fp_kind,0.0_fp_kind,3,2)
		call aberrations(6 )%initialise('C34'      ,'A3','4-Fold astig.    ',0.0_fp_kind,0.0_fp_kind,3,4)
		call aberrations(7 )%initialise('C41'      ,'B4','4th order coma   ',0.0_fp_kind,0.0_fp_kind,4,1)
		call aberrations(8 )%initialise('C43'      ,'D4','3-Lobe aberr.    ',0.0_fp_kind,0.0_fp_kind,4,3)
		call aberrations(9 )%initialise('C45'      ,'A4','5-Fold astig     ',0.0_fp_kind,0.0_fp_kind,4,5)
		call aberrations(10)%initialise('C50 = CS5','C5','5th order spher. ',0.0_fp_kind,0.0_fp_kind,5,0)
		call aberrations(11)%initialise('C52'      ,'S5','5th order star   ',0.0_fp_kind,0.0_fp_kind,5,2)
		call aberrations(12)%initialise('C54'      ,'R5','5th order rosette',0.0_fp_kind,0.0_fp_kind,5,4)
		call aberrations(13)%initialise('C56'      ,'A5','6-Fold astig.    ',0.0_fp_kind,0.0_fp_kind,5,6)
                   
        cutoff = set_cutoff(0.0_fp_kind)
		
        write(*,*) 'Enter the defocus in Angstroms:'
		call get_input('Defocus', set_defocus)
		call read_sequence_string(set_defocus,120,ndf)
		
		allocate(df(ndf))
		call read_sequence_string(set_defocus,120,ndf,df)
        nflag = -1
        apodisation = -1
        do while(nflag.ne.0)
        write(6,11) atan(cutoff/ak1)*1000.0_fp_kind, char(143), cutoff, char(143),df(1)   
		do i=1,ndf-1
			write(6,12) df(i+1)
		enddo
		do i=1,13
			if(aberrations(i)%n==0) then 
				write(6,13) i+2,aberrations(i)%Krivanek,aberrations(i)%Haider,aberrations(i)%Description,aberrations(i)%amplitude
			else
				write(6,13) i+2,aberrations(i)%Krivanek,aberrations(i)%Haider,aberrations(i)%Description,aberrations(i)%amplitude,aberrations(i)%angle
			endif
		enddo
		write(6,14)
		
     11 format(   ' Current simulation parameters:',/,             &
                 &' -----------------------------------------------', /, &
                 &' < 1> Aperture cutoff (mrad)      ', t40, g11.4,/, &
                 &'      Aperture cutoff (', a1, '^-1)              ',t40, g11.4,/, &
                 &' -----------------------------------------------', /, &
				 &'      Symbol convention|                 |             |          ',/, &
				 &'      Krivanek  |Haider|Description      | Size (',a1,')    | Angle(rad)',/,&
				 &' < 2> C10       | C1   |Defocus          |' g13.4,'|')
	 12 format(   '                                         |' g13.4,'|')			 
	 13	format(   1x,'<',i2,'> ',a10,    '| ',a2,3x,'|',a17,          '|',g13.4,    '|',g13.4)
14      format(   ' <16> Output lens contrast transfer function ',t40,/,&
                 &' <17> Apodize probe (remove Airy tails)      ',t40,/,&
                 &' -----------------------------------------------', /, &
                 &' Select option to change any parameter or <0> to continue')
        
        call get_input('Change '//to_lower(trim(adjustl(string)))//' forming lens parameters', nflag)
        write(*,*)
        
        cs =aberrations(4)%amplitude
        if (nflag.eq.1) cutoff = set_cutoff(cs)    
        if (nflag.eq.2) then
			scherzer_df = -1.0_fp_kind*sign(1.0_fp_kind,Cs)*sqrt(4.0_fp_kind*abs(Cs)/(3.0_fp_kind*ak1))
			write(*,*) 'Enter the defocus in Angstroms:'
			if (abs(Cs).gt.1e-3) write(*,'(1x, a, f7.2, a)') '(The optimal Scherzer defocus for the specified Cs is ', scherzer_df, ' Angstroms)'
			
			call get_input('Defocus', set_defocus)
			call read_sequence_string(set_defocus,120,ndf)
			
            if(allocated(df)) deallocate(df)
			allocate(df(ndf))
			call read_sequence_string(set_defocus,120,ndf,df)
        endif
        
        if (nflag.gt.2.and.nflag.lt.16) call aberrations(nflag-2)%set()
		
		if(nflag.eq.16) then
			xyposn = [0,0,0]
            probe = make_ctf(xyposn,df(1),cutoff,aberrations,apodisation)
			call binary_out_unwrap(nopiy,nopix,atan2(imag(probe),real(probe))*abs(probe)**2,&
			                      &trim(adjustl(output_prefix)) //'_'//trim(adjustl(string))//'_forming_lens_ctf_phase')
			call ifft2(nopiy,nopix,probe,nopiy,probe,nopiy)
			call binary_out_unwrap(nopiy,nopix,abs(probe)**2,trim(adjustl(output_prefix))//'_'//trim(adjustl(string))&
			                                                      &//'_forming_lens_ctf_real_space_intensity')
            
        endif
        
        if(nflag.eq.17) apodisation = set_probe_apodisation(cutoff)
        enddo
		if(to_lower(trim(adjustl(string)))=='image') then
			imaging_ndf=ndf
            if(allocated(imaging_df)) deallocate(imaging_df)
			allocate(imaging_df(imaging_ndf))
			call read_sequence_string(set_defocus,120,imaging_ndf,imaging_df)
            imaging_apodisation = apodisation
		else
            if(allocated(probe_df)) deallocate(probe_df)
			probe_ndf=ndf
			allocate(probe_df(probe_ndf))
			call read_sequence_string(set_defocus,120,probe_ndf,probe_df)
            probe_apodisation = apodisation
		endif
    end  subroutine
    
    function set_probe_apodisation(cutoff)
        use global_variables,only:pi
        use m_user_input
        real(fp_kind),intent(in)::cutoff
        real(fp_kind)::set_probe_apodisation
        20      format( ' This options allows the use of a real space aperture to remove the airy tails', /, &
               &' of the probe.' /, &
               &' Please input the radius of the aperture in Angstrom, for an aberration free,', /, &
               &' a probe with cutoff of ',g11.4,1x, a1, '-1, the first, second and third radial minima', /, &
               &' are located at ',g11.4,', ',g11.4,' and ',g11.4,1x, a1, '.',/,&
               &' To disable probe apodisation input a negative number.')
        write(6,20) cutoff,char(143),3.8317/(cutoff*2*pi),7.0156/(cutoff*2*pi),10.1735/(cutoff*2*pi),char(143)
        call get_input('Probe real space cutoff',set_probe_apodisation)
    end function
    
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

    function make_ctf(xyposn,df,cutoff,aberrations,apodisation) result(ctf)

        use global_variables, only: nopiy, nopix, ifactory, ifactorx, pi, ss
        use m_multislice,only:make_detector
        use m_crystallography, only: trimr
        use m_potential, only: make_g_vec_array
        !use output
        use CUFFT_wrapper
        
        implicit none

        complex(fp_kind),dimension(nopiy,nopix) :: ctf,real_space_aperture
        real(fp_kind),intent(in) :: xyposn(3),df,cutoff,apodisation
        real(fp_kind)::deltay,deltax
		type(aberration_coefficient),intent(in)::aberrations(14)
        optional::apodisation

        integer(4) :: ny, nx
        real(fp_kind) :: kr(3), akr, m1, m2,phi, g_vec_array(3,nopiy,nopix)
        
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        ctf=0.0_fp_kind
        !$OMP PARALLEL DO PRIVATE(nx, m1, ny, m2, kr, akr,phi)
        do nx = 1, nopix;do ny = 1, nopiy
            kr = g_vec_array(:,ny,nx)
            akr = trimr(kr,ss)

            if (akr.le.cutoff) then
					phi = atan2(kr(2),kr(1))
                    ctf(ny,nx) = exp(cmplx(0.0_fp_kind, -1*chi(aberrations,akr,phi,df), fp_kind))
            endif
               
        enddo;enddo
        !$OMP END PARALLEL DO
        
        if(present(apodisation)) then
        if(apodisation>0) then
            deltay = 1/trimr(g_vec_array(:,2,1)-g_vec_array(:,1,1),ss)/nopiy
            deltax = 1/trimr(g_vec_array(:,1,2)-g_vec_array(:,1,1),ss)/nopix
            real_space_aperture =make_detector(nopiy,nopix,0.0_fp_kind,apodisation,deltay,deltax)
            call ifft2(nopiy,nopix,ctf,nopiy,ctf,nopiy)
            ctf = ctf*real_space_aperture
            call fft2(nopiy,nopix,ctf,nopiy,ctf,nopiy)
        endif
        endif
        !$OMP PARALLEL DO PRIVATE(nx, m1, ny, m2, kr, akr,phi)
        do nx = 1, nopix;do ny = 1, nopiy
            ctf(ny,nx) = ctf(ny,nx) * exp(cmplx(0.0_fp_kind, -2*pi*dot_product(g_vec_array(:,ny,nx), xyposn), fp_kind))
        enddo;enddo
        !$OMP END PARALLEL DO
    end  function
    subroutine make_lens_ctf(lens_ctf, temp_lens_defocus,aberrations)

        use global_variables, only: nopiy, nopix, ig1, ig2, ifactory, ifactorx, pi, ss
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: lens_ctf(nopiy,nopix)
        real(fp_kind) :: temp_lens_defocus, xyposn(3)
		type(aberration_coefficient),intent(in)::aberrations(14)
        

        integer(4) :: ny, nx, shifty, shiftx
        real(fp_kind) :: kx(3), ky(3), kr(3), akr, m1, m2
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
    
		xyposn = [0,0,0]
		lens_ctf = make_ctf(xyposn,temp_lens_defocus,imaging_cutoff,aberrations,imaging_apodisation)
    end  subroutine


    subroutine make_stem_wfn(psi, defocus, xyposn,aberrations)
    
        use global_variables, only: nopiy, nopix, ifactory, ifactorx, ig1, ig2, ss, pi
        use CUFFT_wrapper, only: ifft2
        use m_crystallography, only: trimr
        
        implicit none

        complex(fp_kind) :: psi(nopiy,nopix)
        real(fp_kind),intent(in) :: defocus, xyposn(3)
		type(aberration_coefficient),intent(in)::aberrations(14)

        integer(4) :: m1,m2,ny,nx,shifty,shiftx
        real(fp_kind) :: ig1_temp(3), ig2_temp(3)
        real(fp_kind) :: kx(3), ky(3), kr(3), akr
        real(fp_kind) :: norm,phi

		psi = make_ctf(xyposn,defocus,probe_cutoff,aberrations,probe_apodisation)
        call ifft2(nopiy, nopix, psi, nopiy, psi, nopiy)
        
        norm = sqrt(sum(abs(psi)**2))

        psi = psi/norm
    end subroutine

end module
