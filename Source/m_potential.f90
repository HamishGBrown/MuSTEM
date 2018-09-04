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
 
module m_potential

    use m_precision, only: fp_kind
    implicit none
    
      integer(4) :: num_ionizations
      integer(4),allocatable::  atm_indices(:)
      character(2),allocatable::ion_description(:)
      logical:: EDX
      
      complex(fp_kind), allocatable :: ionization_mu(:,:,:)        !the ionization scattering factor array, calculated on the grid (supercell)
      complex(fp_kind), allocatable :: fz_adf(:,:,:,:)        !the adf scattering factor array, calculated on the grid (supercell)
      real(fp_kind),    allocatable :: adf_potential(:,:,:,:)
      real(fp_kind),    allocatable :: ionization_potential(:,:,:,:)
      real(fp_kind),    allocatable :: eels_correction_detector(:,:)
	 complex(fp_kind), allocatable :: inverse_sinc_new(:,:)
    
    integer(4) :: n_qep_grates,n_qep_passes,nran ! Start of random number sequence
    
    logical :: phase_ramp_shift
    logical(4) :: quick_shift

    interface
        subroutine make_site_factor_generic(site_factor, tau)
            use m_precision, only: fp_kind
            complex(fp_kind),intent(out) :: site_factor(:, :)
            real(fp_kind),intent(in) :: tau(:,:)
        end subroutine make_site_factor_generic
    end interface
    
    contains
    
    subroutine prompt_high_accuracy
        
        use m_user_input, only: get_input
        use global_variables, only: high_accuracy
        use m_string
        
        implicit none
        
        integer :: i
        
        call command_line_title_box('Potential calculation method')
        write(*,*) 'Two choices are available for the calculation of potentials.'
        write(*,*) 'The reciprocal space method is accurate but may be slower.'
        write(*,*) 'The hybrid method due to Van den Broek et al. is faster but'
        write(*,*) 'is an approximate approach. Note that if "on-the-fly"'
        write(*,*) 'scattering potentials are used the calculation defaults to the'
        write(*,*) 'hybrid approach. '
        write(*,*) '(Van den Broek et al., Ultramicroscopy 158 (2015) pp. 89-97)'
        write(*,*)
        write(*,*) 'Note: if there is insufficient GPU memory, this choice will be overridden.'
        write(*,*)
    10  write(*,*) 'Please choose a method:'
        write(*,*) '<1> Reciprocal space (accuracy)'
        write(*,*) '<2> Hybrid (speed)'
        call get_input('Scattering factor accuracy', i)   
        write(*,*) 
    
        if (i.eq.1) then
            high_accuracy = .true.
            
        elseif (i.eq.2) then
            high_accuracy = .false.
            
        else
            goto 10
            
        endif
        
    end subroutine prompt_high_accuracy
    

    
    subroutine precalculate_scattering_factors
        
        use m_crystallography
	    use m_precision, only: fp_kind
        use global_variables
        use m_electron, only: elsa_ext,peng_ionic_ff,element
        use m_absorption!, only: complex_absorption, setup_absorptive_array, max_int, delta_kstep, tdsbr, fz_abs,calculate_absorption_mu,include_absorption
		use m_numerical_tools, only: cubspl,ppvalu
        use output
        use m_string

	    implicit none
    
        integer(4) :: i, j, k,Z
        real(fp_kind) :: el_scat,ax,ay,g2,s2,sky,skx
        real(fp_kind) :: xkstep, temp
        real(fp_kind),allocatable :: xdata(:),tdsbrc(:,:,:) 
        real(fp_kind) :: factor, eps, g_vec_array(3,nopiy,nopix)

        if(allocated(sinc)) deallocate(sinc)
        if(allocated(inverse_sinc)) deallocate(inverse_sinc)
        if(allocated(inverse_sinc_new)) deallocate(inverse_sinc_new)
        if(allocated(fz)) deallocate(fz)
        if(allocated(fz_DWF)) deallocate(fz_DWF)
    
        allocate(fz(nopiy,nopix,nt))
        allocate(sinc(nopiy,nopix))
        allocate(inverse_sinc(nopiy,nopix))
        allocate(inverse_sinc_new(nopiy,nopix))
        allocate(fz_DWF(nopiy,nopix,nt))
    
        ax = (a0(1)*float(ifactorx))/(float(nopix)*2.0_fp_kind)
        ay = (a0(2)*float(ifactory))/(float(nopiy)*2.0_fp_kind)
        factor = 1.0_fp_kind
        eps = tiny(0.0_fp_kind)
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        
	    do i = 1, nopiy;do j = 1, nopix
            sky = trimr([0.0_fp_kind,g_vec_array(2,i,j),0.0_fp_kind],ss)
            skx = trimr([g_vec_array(1,i,j),0.0_fp_kind,0.0_fp_kind],ss)
            g2 =  trimr(g_vec_array(:,i,j),ss)**2
            s2 = g2 / 4.0_fp_kind
                
            do k = 1, nt
                ! Multiply by fractional occupancy
                if (.not. ionic) el_scat = elsa_ext(nt,k,atomf,s2) * atf(2,k)   
			    if(ionic) el_scat = Peng_ionic_FF(s2,nint(atf(1,k)),dZ(k)) * atf(2,k)        
                    
                
	            ! Fill the potential matrix. Note: these are U(g)/2K
                fz(i,j,k) = cmplx( el_scat, 0.0_fp_kind ,fp_kind)
                fz_DWF(i,j,k) = cmplx( exp( -tp**2.0_fp_kind*g2*atf(3,k) / 2.0_fp_kind ), 0.0_fp_kind,fp_kind ) 
            enddo
            
            !Sinc 
            sinc(i,j) = cmplx((sin(tp*skx*ax)+eps)/(tp*skx*ax+eps)*((sin(tp*sky*ay)+eps)/(tp*sky*ay+eps)),0.0_fp_kind,fp_kind)
            inverse_sinc(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)/(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)
            inverse_sinc_new(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)/(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)            
        enddo; enddo
        ! Currently have U(g)/2K, so multiply by 2K
        fz = 2*ak*fz
        
        ! Normalise the sinc functions
        inverse_sinc = inverse_sinc*float(nopiy)*float(nopix)
        sinc = sinc / (float(nopiy)*float(nopix))
    end subroutine precalculate_scattering_factors
    
    subroutine make_site_factor_matmul(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, tp
        use m_crystallography
    
	    implicit none
    
        !output
        complex(fp_kind),intent(out) :: site_factor(:, :)
        
        !input
        real(fp_kind),intent(in) :: tau(:,:)
    
        integer :: i, j
    
        integer :: g_vec_array(3,nopiy,nopix)
        
        call make_g_vec_array(g_vec_array)
        
        !$OMP PARALLEL PRIVATE(i, j)
        !$OMP DO
	    do i = 1, nopiy;do j = 1, nopix
                site_factor(i,j) = sum(exp(cmplx(0.0_fp_kind, -tp*matmul(g_vec_array(:,i,j), tau), fp_kind)))
	   	 enddo;enddo
	    !$OMP END DO
        !$OMP END PARALLEL
        
    end subroutine make_site_factor_matmul
    
    
#ifdef GPU    
    subroutine make_site_factor_cuda(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, tp
        use cuda_potential, only: cuda_site_factor
        use cuda_array_library, only: blocks, threads
		use m_crystallography,only:make_g_vec_array
        
	    implicit none
    
        !output
        complex(fp_kind),intent(out) :: site_factor(:, :)
        
        !input
        real(fp_kind),intent(in) :: tau(:,:)
    
        integer :: g_vec_array(3,nopiy,nopix)
    
        integer,device :: g_vec_array_d(3,nopiy,nopix)
        real(fp_kind),device :: tau_d(3, size(tau,2))
        complex(fp_kind),device :: site_factor_d(nopiy,nopix)
                          
        call make_g_vec_array(g_vec_array)
		g_vec_array_d = g_vec_array
        tau_d = tau
        
        call cuda_site_factor<<<blocks,threads>>>(site_factor_d, tau_d, g_vec_array_d, nopiy, nopix)        

        site_factor = site_factor_d
    end subroutine make_site_factor_cuda
#endif    
    
    subroutine make_site_factor_hybrid(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix
        use CUFFT_wrapper, only: fft2
        
        implicit none
    
        complex(fp_kind),intent(out) :: site_factor(:,:)
        
        real(fp_kind),intent(in) :: tau(:,:)
        
        integer :: j
        integer :: xpixel, ypixel
        real(fp_kind) :: xpos, ypos, fracx, fracy
    
        site_factor = 0.0_fp_kind
        
        !$OMP PARALLEL PRIVATE(xpos, ypos, j, xpixel,ypixel,fracx,fracy)
        !$OMP DO
        do j = 1, size(tau, 2)
            xpos = tau(1,j)*nopix
            ypos = tau(2,j)*nopiy
            
            ! Ensure that the pixel positions are in range
            
            if (ceiling(xpos).gt.nopix) then
                xpos = xpos - float(nopix)
            elseif (floor(xpos).lt.1) then
                xpos = xpos + float(nopix)
            endif
            
            if (ceiling(ypos).gt.nopiy) then
                ypos = ypos - float(nopiy)
            elseif (floor(ypos).lt.1) then
                ypos = ypos + float(nopiy)
            endif
            
            !fraction of the pixel top right
            xpixel = ceiling(xpos)
            ypixel = ceiling(ypos)
            fracx = mod(xpos, 1.0_fp_kind)
            fracy = mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel top left
            xpixel = floor(xpos)
            ypixel = ceiling(ypos) 
            fracx = 1.0_fp_kind - mod(xpos, 1.0_fp_kind)
            fracy = mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel bottom right
            xpixel = ceiling(xpos)
            ypixel = floor(ypos)
            fracx = mod(xpos, 1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos,1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
            
            !fraction of the pixel bottom left
            xpixel = floor(xpos)
            ypixel = floor(ypos)
            fracx = 1.0_fp_kind - mod(xpos, 1.0_fp_kind)
            fracy = 1.0_fp_kind - mod(ypos, 1.0_fp_kind)
            
            call pixel_check(xpixel, ypixel)
            site_factor(ypixel,xpixel) = site_factor(ypixel,xpixel) + fracx*fracy
        enddo
        !$OMP end do
        !$OMP end parallel
        !fix pixel offset
        site_factor = cshift(site_factor,SHIFT = -1,DIM=1)
        site_factor = cshift(site_factor,SHIFT = -1,DIM=2)

        call fft2(nopiy, nopix, site_factor, nopiy, site_factor, nopiy)
        site_factor = site_factor * inverse_sinc_new * sqrt(float(nopiy)*float(nopix))
        
        contains
        
        subroutine pixel_check(x, y)
            ! Wrap pixel coordinates around so that they remain in range.
    
            implicit none
    
            integer(4) :: x,y
    
            if(x.eq.0) x = nopix
            if(x.eq.nopix+1) x = 1
            if(y.eq.0) y = nopiy
            if(y.eq.nopiy+1) y = 1
    
        end subroutine pixel_check
    
    end subroutine make_site_factor_hybrid

    subroutine setup_inelastic_ionization_types()
        use global_variables    
        use m_user_input
        use m_string
        implicit none
		integer*4::i_eels
      
        call command_line_title_box('Ionization')
		i_eels = 0
		do while(i_eels<1.or.i_eels>2)
			write(*,*) char(10),' <1> EELS',char(10),' <2> EDX',char(10)
			call get_input('Ionization choice', i_eels)
        enddo
        
		EDX = i_eels.eq.2
		call local_potential(EDX)

      end subroutine
	
    function get_ionization_shell_line(shell,atno) result(lineno)
        !Get the line for the given ionization shell and atom number,returns 
        !-1 if that shell is not contained in the parameterization
        use m_string
        character(2),intent(in)::shell
		integer*4,intent(in)::atno
        
        character*31::filename,line
        character(len=:),allocatable ::string 
        
        integer*4::reason,lineno,l
    
        !open the pertinent data files
        
        filename = 'ionization_data\EELS_EDX_'//shell//'.dat'
		open(unit=35,file=filename,status='old',err=970)
        
        l = len('Z = '//to_string(int(atno))//' ')
        allocate(character(l)::string)
        string = 'Z = '//to_string(int(atno))//' '
        lineno = 1
        DO
           READ(35,960,IOSTAT=Reason) line
960        format(a31)
           IF (Reason > 0)  THEN
                write(*,*) 'Problem reading ',filename           
           !If end of file is reach return -1 (shell is not in parametrization)
           ELSE IF (Reason < 0) THEN
              lineno = -1
              close(35)
              return
           ELSE
              !If the substring Z = atno, then return this line 
              if(index(line,string)>0) then
                  close(35)
                  return
              endif
              lineno = lineno+1
           END IF
        END DO
        close(35)
        return
970 write(*,*) 'Problem reading ',filename           
    end function
        
	function get_ionization_parameters(shell,atno,DE,EDX) result (EELS_EDX_params)
        !Read ionisation form factor parameters from ionization_data files
        !shell should be a string describing the orbital ie, '1s', '2s', '2p' etc
        !atno is the atomic number
        !DE is the energy window for EELS and is ignored if the EDX parameterization is requested
        !EDX is a boolean variable, pass .TRUE. for EDX parameterization and .FALSE. for EELS
        use m_numerical_tools
		use m_string
        use global_variables,only: ekv,ak1,nt,atf,substance_atom_types
        
        character(2),intent(in)::shell
		integer*4,intent(in)::atno
		real(fp_kind),intent(in)::DE
		logical,intent(in):: EDX
		
		real(fp_kind)::params(29,8,5),EELS_PARAM_SET2(5,29),xdata(8),bscoef_(4,8)
		real(fp_kind) :: dedata(5),bscoef2_(4,5),p(29),EELS_EDX_params(29)
        
        integer*4::iatom,atno_check,i,ii,m,ishell,iz,lineno,n,mm
		character(10) junk
		character(1) cjunk1,cjunk2,cjunk3

		m = 5
		if(EDX) m=1
        
        n = str2int(shell(1:2))
        lineno = get_ionization_shell_line(shell,atno)
		
		!open the pertinent data files and read to relevant line        
		open(unit=16,file='ionization_data\EELS_EDX_'//shell//'.dat',status='old',err=970)
		do iz = 1,lineno
		   read(16,*) junk
        enddo
        p = 0
        !Later parametrizations only contain 28 datapoints
        mm=29;if (n>2) mm=28
        !Read parameters
		do i=1,8 !Loop over accelerating voltages
		  read(16,*) junk ! E=xx kev header
		  do ii=1,6 !Loop over energy loss above threshhold (EELS) and EDX 
			   read(16,*) p(1:mm)
			   if((.not.EDX).and.(ii<6)) params(:,i,ii) = p(:) !EELS is the first 5 lines
			   if(EDX.and.(ii==6)) params(:,i,1) = p(:) !EDX is the last line
		  enddo
        enddo
        !Can close parameters file now
        close(16)
	  !Interpolate to accelerating voltage used
      !data in files is in steps of 50 keV with 8 points
      !this is stored in xdata
      xdata =(/(i*50, i=1,8,1)/)

       do ii=1,m
          do i=1,mm
			 bscoef_(1,:) = params(i,1:8,ii)
			 call cubspl ( xdata, bscoef_(:,:), 8, 0, 0 )
			 EELS_param_set2(ii,i) = ppvalu(xdata,bscoef_(:,:),7,4,ekv,0)
          enddo
       enddo
	   
	   !If EDX then no energy window interpolation is needed
	  if (EDX) then
		EELS_EDX_params = EELS_param_set2(1,:)
		return
	  endif
	  
      ! contained within EELS_param_set2(i,ii) is the 29 data points (first index) intepolated
      ! to the correct incident energy there are 5 rows Interpolate to energy window desired
      dedata = real([1,10,25,50,100],kind=fp_kind)

      !f(s)/DE is mostly flat and interpolates more simply
      EELS_EDX_params=0
      do i=1,mm
			do ii=1,5
				bscoef2_(1,ii) = EELS_param_set2(ii,i) / dedata(ii)
			enddo
			call cubspl ( dedata, bscoef2_(:,:), 5, 0, 0 )
			EELS_EDX_params(i) = DE*ppvalu(dedata,bscoef2_(:,:),4,4,DE,0)
      enddo
		
		return
970 	write(*,*) ' Cannot access data file EELS_EDX_'//shell//'.dat'
		stop
	end function
	
      !********************************************************************************
      !     subroutine EELS_local_potential()
      !     reads in the scattering factors and performs the interpolation
      !     necessary to accommodate arbitrary geometries and energy windows
      !********************************************************************************
      subroutine local_potential(EDX)

      use m_string
      use m_numerical_tools, only: cubspl,ppvalu
        use m_multislice
	 use global_variables
        use m_user_input

      implicit none
	  logical,intent(in)::EDX

      integer(4) i,ii,iii,j,kval,m,nchoices,ZZ,nshells,norbitals,k
      integer(4),allocatable::available_shells(:),available_atoms(:)

      real(fp_kind),allocatable:: DE(:)
      real(fp_kind) eels_inner,eels_outer
      character(2) shell_name_EELS(9),orb
      character(3) shells(9)!(9)
      character(13):: contributions(4)
      logical,allocatable::choices(:)
      logical::k_shell,l_shell,EDXpresent(nt,4)

      shell_name_EELS = ['1s','2s','2p','3s','3p','3d','4s','4p','4d']
      shells = ['K','L1','L23','M1','M23','M45','N1','N23','N45']
      contributions = ['1s','2s and 2p','3s, 3p and 3d','4s, 4p and 4d']
      norbitals = size(shell_name_EELS)
      nshells = size(shells)
	  !If EELS, setup detectors
	  if(.not.EDX) then
			write(6,91)
	   91 format(1x,'The EELS calculations assume the local approximation, which', /, &
				&1x,'may be inappropriate when the the EELS detector does not', /, &
				&1x,'have a very large acceptance angle. To account for the', /, &
				&1x,'finite detector size, a correction is applied.', /, &
				&1x,'For more details see Y. Zhu et al. APL 103 (2013) 141908.', /)
		
		  eels_inner = 0.0_fp_kind
		  
		  write(*,*) 'EELS detector outer angle (mrad):'
		  call get_input('Outer EELS angle', eels_outer)
		  
		  write(*,*)

		  eels_inner = ak1*tan(eels_inner/1000.0_fp_kind)
		  eels_outer = ak1*tan(eels_outer/1000.0_fp_kind)
		  if(allocated(eels_correction_detector)) deallocate(eels_correction_detector)
		  allocate(eels_correction_detector(nopiy,nopix))
          eels_correction_detector = make_detector(nopiy,nopix,ifactory,ifactorx,ss,eels_inner,eels_outer)
      endif  
!Count available orbitals by checking what is available for the given atoms in
!the parametrization files
	ii=0
        
    do i = 1, nt; ZZ=nint(ATF(1,i)); do j=1,norbitals
        if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) ii=ii+1
    enddo; enddo
    !endif
    
	nchoices = ii
	allocate(choices(ii),DE(ii))
	DE = 0
    choices = .false.
	
	kval = -1
    write(*,*) char(10),' ',char(230),'STEM calculates EDX and EELS signals for ionization of electrons to the '
    write(*,*) 'continuum, at this point bound->bound (white line) transitions are not taken'
    write(*,*) 'into account.',char(10)
    write(*,*) 'K, L, M and N shell ionizations are available though users should be aware '
    write(*,*) 'that quantitative agreement between simulation and theory has only been '
    write(*,*) 'demonstrated for K and L shells (see Y. Zhu and C. Dwyer, Microsc. Microanal.'
    write(*,*) '20 (2014), 1070-1077)',char(10)
	do while ((kval.ne.0).or.all(.not.choices))
100   format(/,' Ionization choices',/,/,'Index  Atom| Z  |',a,'| Included(y/n)'/&
			&,'-----------------------------------------------------------')
	if(EDX) write(*,100) ' Orbital | Shell '
110 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,5x,'|',1x,a3,3x,'|',1x,a1,6x)
    
    if(.not.EDX) write(*,100) ' Orbital | Shell | Window (eV)'
111 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,5x,'|',1x,a3,3x,'|',1x,f5.1,6x,'|',1x,a1,6x)
120 format(' < 0> continue')
    
    !Display choices for EDX and EELS ionizations
	ii=1
      do i = 1, nt;ZZ=nint(ATF(1,i)); do j=1,norbitals
        
            if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) then
                if(EDX) write(*,110) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),shell_name_EELS(j),shells(j),logical_to_yn(choices(ii))
                if(.not.EDX) write(*,111) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),shell_name_EELS(j),shells(j),DE(ii),logical_to_yn(choices(ii))
                ii=ii+1
            endif
        enddo;enddo
      
      write(*,120)
	  
      call get_input('Shell choice <0> continue', kval)
	  !Update choice
	  if ((kval.gt.0).and.(kval.le.nchoices)) then
        choices(kval) = .not.choices(kval)
        
        !If EELS get energy window
        if (.not.EDX) then
            DE(kval) =-1
            do while ((DE(kval).lt.0).or.(DE(kval).gt.100)) 
                write(*,*) 'Enter EELS energy window above threshold in eV (between 1 and 100 ev):',char(10)
                call get_input('Energy window', DE(kval))
            enddo 
        end if
      end if
    enddo
    
	num_ionizations = count(choices)
    
	allocate(ionization_mu(nopiy,nopix,num_ionizations),atm_indices(num_ionizations),Ion_description(num_ionizations))
    ionization_mu = 0
	ii=1
    iii=1
	!Now read in EELS or EDX parameters
    do i = 1, nt;ZZ=nint(ATF(1,i))
    do j=1,norbitals
        if(get_ionization_shell_line(shell_name_EELS(j),ZZ)>-1) then; if(choices(ii)) then
            ionization_mu(:,:,iii) = make_fz_EELS_EDX(shell_name_EELS(j),zz,DE(ii),EDX)* atf(2,i)*fz_DWF(:,:,i)
            atm_indices(iii) = i
            Ion_description(iii) = shell_name_EELS(j)
            iii=iii+1;endif; ii=ii+1; endif
        enddo
    enddo;
 
    end subroutine
      
    !Subrotuine to make the Fz_mu needs to have prefactors accounted for (volume fo the unit cell etc.)
    !needs to be multiplied by the DWF for the pertinent atom type
    function make_fz_EELS_EDX(orbital,zz,DE,EDX) result(fz_mu)
	use m_precision
    use global_variables
	use m_numerical_tools, only: cubspl,ppvalu
    use m_crystallography,only:trimr,make_g_vec_array
	implicit none
    
    character(2),intent(in)::orbital
    integer*4,intent(in)::zz
    real(fp_kind),intent(in)::DE
    logical,intent(in)::EDX
    real(fp_kind):: g_vec_array(3,nopiy,nopix)
    
    complex(fp_kind):: fz_mu(nopiy,nopix)
    
    !dummy variables
    integer(4) i,j
    real(fp_kind) sval

    real(fp_kind) svals(29),EELS_EDX_bscoef(4,29)
    !DATA POINTS USED FOR THE INTERPOLATION S-VALUES (q/2)
    data svals / 0.0_fp_kind,0.025_fp_kind,0.05_fp_kind,0.1_fp_kind,0.2_fp_kind,0.3_fp_kind,0.4_fp_kind,0.5_fp_kind,0.625_fp_kind,&
               & 0.75_fp_kind,0.875_fp_kind,1.0_fp_kind,1.5_fp_kind,2.0_fp_kind,2.5_fp_kind,3.0_fp_kind,3.5_fp_kind,4.0_fp_kind,  &
               & 5.0_fp_kind,6.0_fp_kind,7.0_fp_kind,8.0_fp_kind,9.0_fp_kind,10.0_fp_kind,12.0_fp_kind,14.0_fp_kind,16.0_fp_kind, &
               & 18.0_fp_kind,20.0_fp_kind /
    
    write(*,*) 'Making the ionization inelastic scattering factor grid, please wait...',char(10)
    	
	!pppack interpolation
	EELS_EDX_bscoef(1,:)= get_ionization_parameters(orbital,zz,DE,EDX)
	call cubspl(svals,EELS_EDX_bscoef(:,:), 29, 0, 0 )
    
    fz_mu = 0.0_fp_kind
    call make_g_vec_array(g_vec_array,ifactory,ifactorx)
    !!$OMP PARALLEL PRIVATE(i, j, m2, m1, sky, skx, tempval, sval), SHARED(fz_mu) 
    !!$OMP DO
	do i=1, nopiy;do j=1, nopix
        sval = trimr(g_vec_array(:,i,j),ss) / 2.0_fp_kind
        if (sval.le.20.0_fp_kind) fz_mu(i,j) = cmplx(ppvalu(svals,EELS_EDX_bscoef(:,:),28,4,sval,0),0.0_fp_kind ,fp_kind) / (tp * ak1)
    enddo;enddo
    !!$OMP END DO
	!!$OMP END PARALLEL

    return
    end function

      !--------------------------------------------------------------------------------------
      !   make_mu_matrix() makes the mu matrices for each HOLZ slice
      !   subroutine to take the unit cell input, 
      !   and slice based on holz
      subroutine make_local_inelastic_potentials(ionization)
      
      use m_multislice
      use global_variables!, only:adf,nopiy,nopix,high_accuracy,nt,ss,ndet,ig1,ig2,ifactory
      use m_absorption
      !use m_string
      !use output
      
      implicit none
      
      logical,intent(in)::ionization
      integer(4)   i,j,k,nat_
      real(fp_kind) :: potential_matrix_complex(nopiy,nopix),vol
      real(8)::thmin,thmax,phmin,phmax
      complex(fp_kind)::fz_adf(nopiy,nopix,nt,ndet)

      write(6,134)
134   format(/,' Calculating effective inelastic potentials.',/)

      if(allocated(adf_potential)) deallocate(adf_potential)
      if(allocated(ionization_potential)) deallocate(ionization_potential)

      allocate(adf_potential(nopiy,nopix,n_slices,ndet))            !the adf potential
      adf_potential=0
      if(ionization) allocate(ionization_potential(nopiy,nopix,num_ionizations,n_slices))     !the ionization potential
      
      !if(adf) call make_fz_adf()
	  if(adf.and.complex_absorption) then  
      do k=1,ndet/nseg
          thmin =  atan(inner((k-1)/nseg+1)/ak)
          thmax =  atan(outer((k-1)/nseg+1)/ak)
          !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
          !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
          fz_adf(:,:,:,(k-1)*nseg+1:k*nseg) = spread(absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog,thmin,thmax),dim=4,ncopies=nseg)/nseg
          
      enddo
      endif

      do j = 1, n_slices
	      vol = ss_slice(7,j)
	      !calculate the ionization potential
	      if(ionization) then
              
				do i=1,num_ionizations
                    nat_ = nat_slice(atm_indices(i),j)
                    ionization_potential(:,:,i,j) = real(potential_from_scattering_factors(ionization_mu(:,:,i),tau_slice(:,atm_indices(i),:nat_,j),nat_,nopiy,nopix,high_accuracy)/vol)
				enddo
               
          endif  
	      !calculate the ADF potential
            if(adf.and.complex_absorption) then  
                
                  do i=1,nt
                      nat_ = nat_slice(i,j)
                      do k=1,ndet
                        adf_potential(:,:,j,k)= adf_potential(:,:,j,k) + real(potential_from_scattering_factors(fz_adf(:,:,i,k),tau_slice(:,i,:nat_,j),nat_,nopiy,nopix,high_accuracy)/vol*ss(7)*4*pi)
                      enddo
                  enddo
            endif
      enddo	!ends loop over the number of potential subslices
      
      end subroutine
      
    !--------------------------------------------------------------------------------------
    function potential_from_scattering_factors(scattering_factor,atom_posn,nat_layer,nopiy,nopix,high_accuracy) result(slice_potential)
    use m_precision
    use cufft_wrapper
    
    implicit none
    
    integer(4),intent(in) :: nat_layer,nopiy,nopix
    real(fp_kind),intent(in) :: atom_posn(3,nat_layer)
    complex(fp_kind),intent(in)::scattering_factor(nopiy,nopix)
	logical,intent(in),optional::high_accuracy
    
    complex(fp_kind),dimension(nopiy,nopix) :: potential, site_term,slice_potential
	logical::high_accuracy_

    procedure(make_site_factor_generic),pointer :: make_site_factor
    
	high_accuracy_ = .false.;if(present(high_accuracy)) high_accuracy_= high_accuracy
#ifdef GPU
    make_site_factor => make_site_factor_cuda
#else
    make_site_factor => make_site_factor_matmul
#endif        
    
    slice_potential = 0.0_fp_kind
    
    if (nat_layer.ne.0) then
        if (high_accuracy_) then
            call make_site_factor(site_term, atom_posn)
        else
            call make_site_factor_hybrid(site_term, atom_posn)        
        endif
        slice_potential = site_term*scattering_factor
        ! Get realspace potential
        call ifft2(nopiy,nopix,slice_potential,nopiy,slice_potential,nopiy)
        slice_potential = slice_potential*sqrt(float(nopiy*nopix))
    endif
    
    end function
    function make_absorptive_grates(nopiy,nopix,n_slices) result(projected_potential)
    
        use m_precision, only: fp_kind
	    use cufft_wrapper, only: fft2, ifft2
        use global_variables, only: ig1,ig2,ifactory,ifactorx,nt, relm, tp, ak, atf, high_accuracy, ci, pi, bwl_mat,fz,fz_DWF,ss,a0,nat,orthog
        use m_absorption!, only: transf_absorptive,fz_abs
        use m_multislice
        use m_string, only: to_string
        use output
        
        implicit none
        
        integer*4,intent(in)::nopiy,nopix,n_slices
        complex(fp_kind)::projected_potential(nopiy,nopix,n_slices)
        
        integer(4) :: j, m, n,nat_layer
        real(fp_kind) :: ccd_slice,V_corr
        complex(fp_kind),dimension(nopiy,nopix) :: scattering_pot,temp,effective_scat_fact
        complex(fp_kind)::fz_abs(nopiy,nopix,nt)
    
        real(fp_kind) :: t1, delta,amplitude(nopiy,nopix),phase(nopiy,nopix)
    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        projected_potential= 0 
        t1 = secnds(0.0_fp_kind)
        fz_abs=0
        if(include_absorption) fz_abs = absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,nt,a0,ss,atf,nat, ak, relm, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak
            
        do j = 1, n_slices
	        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
    198	    write(6,199) to_string(sum(nat_slice(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a, /) 

		    ccd_slice = relm / (tp * ak * ss_slice(7,j))
            V_corr = ss(7)/ss_slice(7,j)
            do m=1,nt
                nat_layer = nat_slice(m,j)
                effective_scat_fact = CCD_slice*fz(:,:,m)*fz_DWF(:,:,m)+cmplx(0,1)*fz_abs(:,:,m)*V_corr
                projected_potential(:,:,j) = projected_potential(:,:,j)+potential_from_scattering_factors(effective_scat_fact,tau_slice(:,m,:nat_layer,j),nat_layer,nopiy,nopix,high_accuracy)
            enddo
                
	    enddo ! End loop over slices
	
	    delta = secnds(t1)
        
		if(timing) then
            write(*,*) 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
            write(*,*)
			open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
			write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
			close(9834)
        endif    
        
    end function make_absorptive_grates
    
           
        integer function seed_rng() result(idum)
            
	        use m_numerical_tools, only: ran1
            use m_precision, only: fp_kind
            
            implicit none

            integer :: i
            real(fp_kind) :: random
            
	        idum = -1
            
	        do i = 1, nran
		        random = ran1(idum)
	        enddo
            
        end function seed_rng

      
    subroutine make_propagator(nopiy,nopix,prop,dz,ak1,ss,ig1,ig2,claue,ifactorx,ifactory)

        use m_precision, only: fp_kind
        use m_crystallography, only: trimr,make_g_vec_array
            
        implicit none

        integer(4) :: nopiy,nopix
        complex(fp_kind) :: prop(nopiy,nopix)        
        real(fp_kind) :: ak1, ss(7), claue(3), dz, g_vec_array(3,nopiy,nopix)
        integer(4) :: ifactorx, ifactory, ig1(3), ig2(3)

        real(fp_kind),parameter :: pi = atan(1.0d0)*4.0d0
        integer(4) :: ny, nx
        
        
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)

        do ny = 1, nopiy;do nx = 1, nopix
            prop(ny,nx) = exp(cmplx(0.0d0, -pi*dz*trimr(g_vec_array(:,ny,nx)-claue,ss)**2/ak1, fp_kind ))
        enddo;enddo

    end subroutine
       
    function make_qep_grates(idum) result(projected_potential)
    
        use m_precision, only: fp_kind
	    use cufft_wrapper, only: fft2, ifft2
        use global_variables, only: nopiy, nopix, nt, relm, tp, ak, ak1, atf, high_accuracy, ci, pi, bwl_mat,fz
        use m_multislice
        use m_string, only: to_string
        use output, only: output_prefix,timing,binary_in
        use m_numerical_tools, only: displace

        implicit none
        
        integer(4),intent(inout) :: idum
    
        complex(fp_kind) :: projected_potential(nopiy,nopix,n_qep_grates,n_slices),temp(nopiy,nopix),scattering_pot(nopiy,nopix,nt)
		integer(4), allocatable :: handled(:,:)
		integer(4):: save_list(2,nt),match_count, i, j, m, n,ii,jj,jjj,kk,iii
        real(fp_kind) :: tau_holder(3),tau_holder2(3),ccd_slice,ums,amplitude(nopiy,nopix),phase(nopiy,nopix)
        real(fp_kind) :: mod_tau(3,nt,maxnat_slice,n_slices,n_qep_grates),t1, delta
        logical::fracocc
    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        
	    
 	 	!	Search for fractional occupancy
         fracocc = any(atf(2,:).lt.0.99d0)
        
        t1 = secnds(0.0_fp_kind)

        do j = 1, n_slices
	        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
    198	    write(6,199) to_string(sum(nat_slice(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a) 

		    ccd_slice = relm / (tp * ak * ss_slice(7,j))

	        do i = 1, n_qep_grates
    200	        format(a1, 1x, i3, '/', i3)
	            write(6,200, advance='no') achar(13), i, n_qep_grates
          
                ! Randomly displace the atoms
				if (.not.fracocc) then
 	            do m = 1, nt
	                do n = 1, nat_slice(m,j)
			            call displace(tau_slice(1:3,m,n,j),mod_tau(1:3,m,n,j,i),sqrt(atf(3,m)),a0_slice,idum)
	                enddo
                enddo
				else
					allocate( handled(nt,maxnat_slice) )
					handled = 0
					do ii=1, nt
					 do jj = 1, nat_slice(ii,j)
						 if (handled(ii,jj).eq.1) cycle
						 tau_holder(1:3) = tau_slice(1:3,ii,jj,j)

						 save_list = 0
						 match_count = 0
						 ums = atf(3,ii)
						 do iii=ii+1,nt
						 do jjj=1,nat_slice(iii,j)
							if (same_site(tau_holder,tau_slice(1:3,iii,jjj,j))) then
							   match_count = match_count+1
							   save_list(1,match_count)=iii
							   save_list(2,match_count)=jjj
							   ums = ums + atf(3,iii)
							   cycle
							endif
						 enddo
						 enddo

						 ums = ums / dfloat(match_count+1)
					   call displace(tau_holder(1:3),tau_holder2(1:3),sqrt(ums),a0_slice,idum)
						 mod_tau(1:3,ii,jj,j,i) = tau_holder2(1:3)
						 handled(ii,jj) = 1
						 do kk=1,match_count
							 mod_tau(1:3,save_list(1,kk),save_list(2,kk),j,i)&
				                                              &= tau_holder2(1:3)
							 handled(save_list(1,kk),save_list(2,kk)) = 1
						 enddo
						   
					 enddo
					 enddo

					 deallocate( handled )
				endif
				
				projected_potential(:,:,i,j) = 0
				do m = 1, nt
					projected_potential(:,:,i,j) = projected_potential(:,:,i,j)+real(potential_from_scattering_factors(CCD_slice*fz(:,:,m)&
												&,mod_tau(:,m,1:nat_slice(m,j),j,i),nat_slice(m,j),nopiy,nopix,high_accuracy))
                enddo
	        enddo ! End loop over grates
        
            write(*,*)
            write(*,*)
        
	    enddo ! End loop over slices
	
	    delta = secnds(t1)
        
        write(*,*) 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
        write(*,*)

    	if(timing) then
			open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
			write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
			close(9834)
        endif    
    
    end function make_qep_grates

	logical(4) function same_site(site1,site2)
      
      implicit none
      
      real(fp_kind) site1(3),site2(3)
      real(fp_kind) tol
      
      tol = 1.0d-6
      same_site = all(abs(site1-site2).lt.tol)
      
      return
      end function
end module
