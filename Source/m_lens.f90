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
    
    !Future work: make a lens object that contains the aberrations, defocus, cutoff and apodisation
    !and has calculate ctf as a bound method
    ! Probe forming lens
    type(aberration_coefficient)::probe_aberrations(14),imaging_aberrations(14)
    real(fp_kind),allocatable :: probe_df(:),imaging_df(:)
    real(fp_kind)::probe_cutoff,imaging_cutoff
    real(fp_kind)::probe_apodisation,imaging_apodisation
    
    integer :: probe_ndf,imaging_ndf

    integer(4) :: nxsample, nysample          !number of probe positions
    real(fp_kind), allocatable :: probe_positions(:,:,:) !matrix containing the probe position
    real(fp_kind) :: delx, dely                !stepsize for the probe position in x and y
    real(fp_kind) :: probe_initial_position(3) = [0.0_fp_kind, 0.0_fp_kind, 0.0_fp_kind]    
	
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
		

        call command_line_title_box(string//' forming lens')
        
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
		call Series_prompt('defocus')
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

        use global_variables, only: nopiy, nopix, ifactory, ifactorx, pi, ss,a0
        use m_multislice,only:make_detector
        use m_crystallography, only: trimr,make_g_vec_array
        use output
        use CUFFT_wrapper
        
        implicit none

        complex(fp_kind),dimension(nopiy,nopix) :: ctf,real_space_aperture
        real(fp_kind),intent(in) :: xyposn(3),df,cutoff,apodisation
        real(fp_kind)::deltay,deltax
		type(aberration_coefficient),intent(in)::aberrations(14)
        optional::apodisation

        integer(4) :: ny, nx
        real(fp_kind) :: kr(3), akr, m1, m2,phi, g_vec_array(3,nopiy,nopix),x,y,r
        
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
            real_space_aperture = 0
            do nx = 1, nopix;x = (nx-nopix/2-2)*a0(2)*ifactorx/nopix; do ny = 1, nopiy
                y = (ny-nopiy/2-2)*a0(1)*ifactory/nopiy;r = sqrt(x**2+y**2)
                if (r.le.apodisation) real_space_aperture(ny,nx) =1
            enddo;enddo
            real_space_aperture = quad_shift(real_space_aperture,nopiy,nopix)
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
    
    subroutine setup_probe_scan(PACBED_only)
    
        use global_variables, only: uvw1, uvw2,tiley,tilex, output_nopiy, output_nopix,interpolation
        use m_string
        
        implicit none

        real(fp_kind) :: fract(2), origin(3)
        logical,intent(in),optional::PACBED_only
        logical::PACBED_only_
        
        
        PACBED_only_ = .false.
        if (present(PACBED_only)) PACBED_only_=PACBED_only
        interpolation = .not.PACBED_only_
        
        call command_line_title_box('Probe scan details')
        write(*,*) 'Warning: changing the following parameters will '
        write(*,*) 'disable interpolation of STEM images.',char(10)
        
        call setup_scan_geometry(fract, origin,interpolation,PACBED_only_)
            
        call calculate_probe_positions(uvw1, uvw2, origin)
        
        if(interpolation) then
            call setup_stem_image_interpolation
        else
            output_nopix = nxsample
            output_nopiy = nysample
            tiley = 1
            tilex = 1
        endif
    end subroutine

    subroutine reset_scan(origin,a1,a2,r1,r2,thetad2,fract,min_step,interpolation,PACBED_only)
        use global_variables, only: a0, deg, ss, uvw1, uvw2, thetad, izone,fourdstem
        use m_crystallography, only: zone, subuvw, angle, rsd
        
        real(fp_kind),intent(out) :: fract(2), origin(3),a1, a2, r1(3), r2(3), thetad2,min_step
        logical,intent(out)::interpolation
        logical,intent(in)::PACBED_only
        
        origin = 0.0_fp_kind
        a1 = rsd(uvw1, a0, deg)
        a2 = rsd(uvw2, a0, deg)
        r1 = uvw1
        r2 = uvw2
        thetad2 = thetad       
        fract = 1.0
        min_step = nyquist_step(probe_cutoff)
        nysample = nyquist_sampling(probe_cutoff, fract(2)*a0(2),PACBED_only.and.(.not.fourDSTEM))
        nxsample = nyquist_sampling(probe_cutoff, fract(1)*a0(1),PACBED_only.and.(.not.fourDSTEM))
        interpolation = .not.PACBED_only
    end subroutine
    
    subroutine setup_scan_geometry(fract, origin,interpolation,PACBED_only)
    
        use m_user_input, only: get_input
        use global_variables, only: a0, deg, ss, uvw1, uvw2, thetad, izone,fourDSTEM
        use m_crystallography, only: zone, subuvw, angle, rsd
        
        implicit none

        real(fp_kind),intent(out) :: fract(2), origin(3)
        logical,intent(in)::PACBED_only
        logical,intent(inout)::interpolation
        
        integer(4) :: ich
        integer(4)   ig1a(3), ig2a(3),nysample_,nxsample_
        real(fp_kind) a1, a2, r1(3), r2(3), thetad2,min_step
        character(8)::able_string
        character(7)::PACBED_or_STEM
    
        integer :: i_scan_quarter
        
        call reset_scan(origin,a1,a2,r1,r2,thetad2,fract,min_step,interpolation,PACBED_only)
        nysample_ = nyquist_sampling(probe_cutoff, a0(2),PACBED_only.and.(.not.fourDSTEM))
        nxsample_ = nyquist_sampling(probe_cutoff, a0(1),PACBED_only.and.(.not.fourDSTEM))
        ich=-1
        if(fourDSTEM.and.PACBED_only) then
            PACBED_or_STEM = '4D-STEM'
        elseif(PACBED_only) then
            PACBED_or_STEM = 'PACBED'
        else
            PACBED_or_STEM = 'STEM'
        endif
     do while(ich.ne.1)
         if(interpolation) able_string = 'enabled'
         if(.not.interpolation) able_string = 'disabled'
          write(6,103) r1, a1, char(143), r2, a2, char(143), thetad2, origin(1), origin(2), probe_cutoff
          if(.not.PACBED_only) write(6,104)
          write(6,105) PACBED_or_STEM,min_step,nxsample_, nysample_, ceiling(nxsample/fract(1)),ceiling(nysample/fract(2)),nxsample, nysample
          if(.not.PACBED_only) write(6,106) able_string
          write(6,107)
103         format(/,' The probe scan vectors are: ', /,                  & 
            &       ' x = ', 3g12.5, ' mag = ', g12.5, 1x, a1, /,         &
            &       ' y = ', 3g12.5, ' mag = ', g12.5, 1x, a1, /,         &
            &       ' The angle between these scan vectors is ', f12.2, ' degrees', /, &
            &       ' The inital (fractional) position is ', g12.5, ', ', g12.5, /,/, &
			&		' The maximum spatial frequency allowed by the probe is ', f5.2, ' A-1.')
104			format( ' The STEM image has a bandwidth limit of twice that frequency. ')
105			format( ' This corresponds to minimum ',a7,' sampling of ', f6.2, ' positions per Angstrom:', /, &
			&        i4, ' x-positions and ', i4, ' y-positions per unit cell. Currently, sampling is',/,&
            &        i4, ' x-positions and ', i4, ' y-positions per unit cell for a total of' /, &
            &        i4, ' x-positions and ', i4, ' y-positions.')
106         format( ' Interpolation of STEM images is currently ',a)
107         format( ' <1> Accept', /,                  &
            &       ' <2> Change size of x and y (also adjusts sampling).', /,      &
            &       ' <3> Change orientation of x and y.', /,     & 
            &       ' <4> Change the initial position.',/,&
			&       ' <5> Change the probe position sampling.',/,&
            &       ' <6> Reset.',/,&
            &       ' <7> Output probe positions.')

            call get_input("Probe scan menu choice", ich)
        
            write(*,*)
        
            if(ich.eq.2) then             !changing size of x and y vectors
                  write(6,111)
        111       format(' Enter fractional increase in both x and y.',/, 'For example, to double the size of the scan enter 2 2')
                  call get_input("Enter fractional increase in x and y", fract(1),fract(2))

                  r1 = fract(1) * r1
                  r2 = fract(2) * r2      
                  a1 = fract(1) * a1
                  a2 = fract(2) * a2
                  min_step = nyquist_step(probe_cutoff)
                  nxsample = nyquist_sampling(probe_cutoff, a1,PACBED_only.and.(.not.fourDSTEM))
                  nysample = nyquist_sampling(probe_cutoff, a2,PACBED_only.and.(.not.fourDSTEM))
                  interpolation=.false.
            elseif(ich.eq.3) then
        121       format( /, ' Please enter a new x-scan vector.' )
                  write(6,121)
                  call get_input("x-scan vector", ig1a, 3)
                  call zone(izone, ig1a, ig2a)          !get an orthogonal vector to the zone axis and ig1
                  call angle(ig1a, ig2a, ss, thetad2)    !calculate angle between scan vectors
                  call subuvw(ig1a, r1, a0, deg, ss)      !calculate the real space scan vector from the 'ig1' given
                  call subuvw(ig2a, r2, a0, deg, ss)      !calculate the real space scan length from the 'ig2' given
                  interpolation=.false.
            elseif(ich.eq.4) then
                call place_probe(origin)
                  interpolation=.false.
            elseif(ich.eq.5) then
				write(*,*) 'Enter the number of probe positions in the x direction.'
				call get_input('nxsample', nxsample)
				write(*,*) 'Enter the number of probe positions in the y direction.' 
				call get_input('nysample', nysample)
				write(*,*)
				interpolation = .false.
            elseif(ich.eq.6) then
                call reset_scan(origin,a1,a2,r1,r2,thetad2,fract,min_step,interpolation,PACBED_only)
			elseif(ich.eq.7) then 
                call calculate_probe_positions(r1, r2, origin)
                call plot_scan(probe_positions,nysample,nxsample)
            endif
            enddo
            uvw1 = r1
            uvw2 = r2
        
    end subroutine  
    
    subroutine plot_scan(probe_positions,nysample,nxsample)
        real(fp_kind),intent(in)::probe_positions(3,nysample,nxsample)
        integer*4,intent(in)::nysample,nxsample
        integer*4::ny,nx
        
        write(*,*) 'Writing probe positions to file "probe_positions.txt"'
       open(unit=52,file='probe_positions.txt')
467    format(3(f9.4))    
        do ny = 1, nysample;do nx = 1, nxsample
            write(52,467) probe_positions(:,ny,nx)
        enddo;enddo
        close(52)
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

    integer function nyquist_sampling(qmax, L,PACBED) result(n)
    
        implicit none
        
        real(fp_kind) :: qmax, L
        logical,intent(in),optional::PACBED
        logical::PACBED_
        
        PACBED_=.false.
        if(present(PACBED)) PACBED_=PACBED
        
        if(PACBED_) n = ceiling(2 * qmax * L)
        if(.not.PACBED_) n = ceiling(4 * qmax * L)
        
    end function
    
    real(fp_kind) function nyquist_step(qmax,PACBED) result(step)
    
        implicit none
        
        real(fp_kind) :: qmax
        logical,intent(in),optional::PACBED
        logical::PACBED_
        
        PACBED_=.false.
        if(present(PACBED)) PACBED_=PACBED
        
        if(PACBED_) step = 2 * qmax
        if(.not.PACBED_) step = 4 * qmax

    end function
    
    subroutine setup_stem_image_interpolation()
    
        use m_user_input, only: get_input
        use global_variables, only: tiley, tilex, output_nopiy, output_nopix,interpolation
    
        implicit none

        integer(4) :: out_max
        
        if((nysample.gt.1).and.(nxsample.gt.1)) then
            write(6,99) 
            99 format(' Enter the maximum number of pixels to interpolate the output image to.',/,&
                      &' Enter a negative number to disable interpolation.')
            call get_input('output interpolation max pixels', out_max)
            write(*,*)
            if (out_max<0) then
                write(*,*) 'Interpolation has been disabled.',char(10)
                interpolation = .false.
                output_nopix = nxsample
                output_nopiy = nysample
                tiley = 1
                tilex = 1
                return
            endif
            
            
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
				&' undersampling of the output.',/)
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
