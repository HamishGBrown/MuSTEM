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
      complex(fp_kind), allocatable :: fz_adf(:,:,:)        !the adf scattering factor array, calculated on the grid (supercell)
      real(fp_kind),    allocatable :: adf_potential(:,:,:)
      real(fp_kind),    allocatable :: ionization_potential(:,:,:,:)
      real(fp_kind),    allocatable :: eels_correction_detector(:,:)
    

    interface
        subroutine make_site_factor_generic(site_factor, tau)
            use m_precision, only: fp_kind
            complex(fp_kind) :: site_factor(:, :)
            real(fp_kind) :: tau(:,:)
        end subroutine make_site_factor_generic
    end interface
    
    interface make_g_vec_array
        module procedure make_g_vec_array_real,make_g_vec_array_int
    end interface
    
    complex(fp_kind), allocatable :: inverse_sinc_new(:,:)
    
    contains
    
    
    subroutine prompt_high_accuracy
        
        use m_user_input, only: get_input
        use global_variables, only: high_accuracy
        
        implicit none
        
        integer :: i
        
        write(*,*) '|------------------------------------------|'
	    write(*,*) '|      Potential calculation method        |'
	    write(*,*) '|------------------------------------------|'
        write(*,*)
    
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
    
    subroutine make_g_vec_array_real(g_vec_array,ifactory,ifactorx)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, ig1, ig2
        
	    implicit none
    
        real(fp_kind) :: g_vec_array(3,nopiy,nopix)
        integer,intent(in),optional:: ifactory,ifactorx
        
        integer :: shiftx, shifty, m1, m2, i, j
        real(fp_kind)::ifactory_,ifactorx_
        
        ifactory_= 1.0_fp_kind; if(present(ifactory)) ifactory_=ifactory
        ifactorx_= 1.0_fp_kind; if(present(ifactorx)) ifactorx_=ifactorx
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1
        
        !$OMP PARALLEL PRIVATE(i, m2, j, m1)
        !$OMP DO
	    do i = 1, nopiy
	        m2 = mod( i+shifty, nopiy) - shifty -1
	        do j = 1, nopix
	            m1 = mod( j+shiftx, nopix) - shiftx -1
	            g_vec_array(:,i,j) = m1 * ig1/ifactorx_ + m2 * ig2/ifactory_                
	   	    enddo
        enddo
	    !$OMP END DO
        !$OMP END PARALLEL

    end subroutine make_g_vec_array_real
    
    
    
    subroutine make_g_vec_array_int(g_vec_array)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, ig1, ig2
        
	    implicit none
    
        integer :: g_vec_array(3,nopiy,nopix)
        
        integer :: shiftx, shifty, m1, m2, i, j
        
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1
        
        !$OMP PARALLEL PRIVATE(i, m2, j, m1)
        !$OMP DO
	    do i = 1, nopiy
	        m2 = mod( i+shifty, nopiy) - shifty -1
	        do j = 1, nopix
	            m1 = mod( j+shiftx, nopix) - shiftx -1
	            g_vec_array(:,i,j) = m1 * ig1 + m2 * ig2
	   	    enddo
        enddo
	    !$OMP END DO
        !$OMP END PARALLEL

    end subroutine make_g_vec_array_int
    
    
    
    subroutine precalculate_scattering_factors
        
        use m_crystallography
	    use m_precision, only: fp_kind
        use global_variables
        use m_electron, only: elsa_ext,peng_ionic_ff
        use m_absorption, only: complex_absorption, setup_absorptive_array, max_int, delta_kstep, tdsbr, fz_abs,calculate_absorption_mu
		use m_numerical_tools, only: cubspl,ppvalu

	    implicit none
    
        integer(4) :: i, j, k
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
    
        if(complex_absorption) then
            call setup_absorptive_array
            call calculate_absorption_mu
            ! Calculate TDS form factors
            if(allocated(xdata)) deallocate(xdata)
            if(allocated(fz_abs)) deallocate(fz_abs)
        
            allocate(xdata(max_int))
            allocate(fz_abs(nopiy,nopix,nt))


			if(allocated(tdsbrc)) deallocate(tdsbrc)
			allocate(tdsbrc(4,max_int,nt))

            ! Set up interpolation domain
            xkstep = delta_kstep
            do i = 1, max_int
                xdata(i) = (i-1)* xkstep
            enddo
            
            ! Set up spline coefficients for each atomic species
            do k = 1, nt
				tdsbrc(1,:,k) = tdsbr(:,k)
				call cubspl ( xdata, tdsbrc(:,:,k), max_int, 0, 0 )
            enddo
            
        endif

        write(*,*) 'Calculating scattering factors...'
        write(*,*) 
    
        ax = (a0(1)*float(ifactorx))/(float(nopix)*2.0_fp_kind)
        ay = (a0(2)*float(ifactory))/(float(nopiy)*2.0_fp_kind)
        factor = 1.0_fp_kind
        eps = tiny(0.0_fp_kind)
        
        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        
	    do i = 1, nopiy;do j = 1, nopix
            skx = trimr([0.0_fp_kind,g_vec_array(2,i,j),0.0_fp_kind],ss)
            sky = trimr([g_vec_array(1,i,j),0.0_fp_kind,0.0_fp_kind],ss)
            g2 =  trimr(g_vec_array(:,i,j),ss)**2
            s2 = g2 / 4.0_fp_kind
                
            do k = 1, nt
                ! Multiply by fractional occupancy
                if (.not. ionic) el_scat = elsa_ext(nt,k,atomf,s2) * atf(2,k)   
			    if(ionic) el_scat = Peng_ionic_FF(s2,nint(atf(1,k)),dZ(k)) * atf(2,k)        
                    
                
	            ! Fill the potential matrix. Note: these are U(g)/2K
                fz(i,j,k) = cmplx( el_scat, 0.0_fp_kind ,fp_kind)
                fz_DWF(i,j,k) = cmplx( exp( -tp**2.0_fp_kind*g2*atf(3,k) / 2.0_fp_kind ), 0.0_fp_kind,fp_kind ) 
                if(complex_absorption) then
						
				    temp = ppvalu(xdata,tdsbrc(:,:max_int-1,k),max_int-1,4,sqrt(g2),0)
						
                    fz_abs(i,j,k) = cmplx(temp,0.0_fp_kind,fp_kind)
                endif
            enddo
            
            !Sinc 
            sinc(i,j) = cmplx((sin(tp*skx*ax)+eps)/(tp*skx*ax+eps)*((sin(tp*sky*ay)+eps)/(tp*sky*ay+eps)),0.0_fp_kind,fp_kind)
            inverse_sinc(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)/(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)
            inverse_sinc_new(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)/(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)            
        enddo; enddo

        ! Currently have U(g)/2K, so multiply by 2K
        fz = 2*ak*fz
        fz_abs = 2*ak*fz_abs
        
        ! Normalise the sinc functions
        inverse_sinc = inverse_sinc*float(nopiy)*float(nopix)
        sinc = sinc / (float(nopiy)*float(nopix))
        
    end subroutine precalculate_scattering_factors
    
    subroutine make_site_factor_matmul(site_factor, tau)

        use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, tp
    
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
    

     
    subroutine make_qep_potential(potential, tau_slice, nat_slice, CCD_slice, make_site_factor)

	    use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, nt, fz
	    use cufft_wrapper, only: ifft2
        use m_slicing, only: maxnat_slice
    
	    implicit none
    
        !output
        complex(fp_kind) :: potential(nopiy, nopix)
    
        !input
        real(fp_kind) :: tau_slice(3,nt,maxnat_slice) !atom locations for the supercell slice
        integer(4) :: nat_slice(nt)       !number of atoms in this slice (of each type)
        real(fp_kind) :: CCD_slice        
        procedure(make_site_factor_generic),pointer :: make_site_fact
        procedure(make_site_factor_generic),pointer,optional:: make_site_factor
        
        integer :: m
        complex(fp_kind) :: site_factor(nopiy,nopix)
           
        potential = 0.0_fp_kind
        
        if(present(make_site_factor)) then
            make_site_fact => make_site_factor
        else
#ifdef GPU
            make_site_fact => make_site_factor_cuda
#else
            make_site_fact => make_site_factor_matmul
#endif
        endif
        
        do m = 1, nt
            if (nat_slice(m)==0) cycle
               
            call make_site_fact(site_factor, tau_slice(:,m,1:nat_slice(m)))
                       
            potential = potential + site_factor*CCD_slice*fz(:,:,m)
        
        enddo
    
        call ifft2(nopiy, nopix, potential, nopiy, potential, nopiy)
        potential = potential * sqrt(float(nopiy*nopix))
            
    end subroutine make_qep_potential
    

        
    subroutine make_absorptive_potential(potential, tau_slice, nat_slice, CCD_slice, volume_slice, make_site_factor)

	    use m_precision, only: fp_kind
        use global_variables, only: nopiy, nopix, nt, fz, fz_DWF, ss
        use m_absorption, only: fz_abs
	    use cufft_wrapper, only: ifft2
        use m_slicing, only: maxnat_slice
    
	    implicit none
    
        !output
        complex(fp_kind) :: potential(nopiy, nopix)
    
        !input
        real(fp_kind) :: tau_slice(3,nt,maxnat_slice) !atom locations for the supercell slice
        integer(4) :: nat_slice(nt)       !number of atoms in this slice (of each type)
        real(fp_kind) :: CCD_slice, volume_slice    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        
        integer :: m
        complex(fp_kind) :: site_factor(nopiy,nopix)
        real(fp_kind) :: V_corr
        
        V_corr = ss(7)/volume_slice
    
        potential = 0.0_fp_kind
    
        do m = 1, nt
            if (nat_slice(m)==0) cycle
               
            call make_site_factor(site_factor, tau_slice(:,m,1:nat_slice(m)))
                       
            potential = potential + site_factor*(CCD_slice*fz(:,:,m)*fz_DWF(:,:,m)+cmplx(0,1)*fz_abs(:,:,m)*V_corr)

        
        enddo 
    
        call ifft2(nopiy, nopix, potential, nopiy, potential, nopiy)
        potential = potential * sqrt(float(nopiy*nopix))


    end subroutine make_absorptive_potential
    
      
    subroutine setup_inelastic_ionization_types()
        use global_variables    
        use m_user_input
        implicit none
		integer*4::i_eels
      
        write(*,*) '|-----------------------|'
        write(*,*) '|      Ionization       |'
        write(*,*) '|-----------------------|'
        
		i_eels = 0
		do while(i_eels<1.or.i_eels>3)
			write(*,*) char(10),'<1> EELS',char(10),'<2> EDX',char(10),'<3> No ionization',char(10)
			call get_input('Ionization choice', i_eels)
		enddo
      
		ionization = i_eels.ne.3
		EDX = i_eels.eq.2
		if (ionization) call local_potential(EDX)

      end subroutine
      
	function logical_to_yn(input) result(str)
		logical,intent(in)::input
		character(1)::str
		if (input) str = 'y'
		if (.not.input) str = 'n'
	end function     

	
	function get_ionization_parameters(shell,atno,DE,EDX) result (EELS_EDX_params)
        !Read ionisation form factor parameters from ionization_data files
        use m_numerical_tools
		use m_string
        use global_variables,only: ekv,ak1,nt,atf,substance_atom_types
        
        character(2),intent(in)::shell
		integer*4,intent(in)::atno
		real(fp_kind),intent(in)::DE
		logical,intent(in):: EDX
		
		real(fp_kind)::params(29,8,5),EELS_PARAM_SET2(5,29),xdata(8),bscoef_(4,8)
		real(fp_kind) :: dedata(5),bscoef2_(4,5),p(29),EELS_EDX_params(29)
        
        integer*4::iatom,atno_check,i,ii,m,ishell,iz
		character(10) junk
		character(1) cjunk1,cjunk2,cjunk3

		m=5
		if(EDX) m=1
        
		!open the pertinent data files
		open(unit=16,file='ionization_data\EELS_EDX_'//shell//'.dat',status='old',err=970)
        
        ishell= str2int(shell(1:1))
		!The first atom in the parameters file depends on the shell
		if (ishell.eq.1) iatom = 6
		if (ishell.eq.2) iatom = 20

		!Skip elements before
		do iz = 1,58*(atno-iatom)
		   read(16,*) junk
        enddo
        
        !Read element header
		read(16,*) cjunk1,cjunk2,cjunk3,atno_check
		if (atno_check.ne.atno) then
            pause ' Something has gone wrong.  Debugging required.'
        endif
        
        !Read parameters
		do i=1,8 !Loop over accelerating voltages
		  read(16,*) junk ! E=xx kev header
		  do ii=1,6 !Loop over energy loss above threshhold (EELS) and EDX 
			   read(16,*) p(1:29)
			   if((.not.EDX).and.(ii<6)) params(:,i,ii) = p !EELS is the first 5 lines
			   if(EDX.and.(ii==6)) params(:,i,1) = p !EDX is the last line
		  enddo
        enddo
        !Can close parameters file now
        close(16)
	  !Interpolate to accelerating voltage used
      !data in files is in steps of 50 keV with 8 points
      !this is stored in xdata
      xdata =(/(i*50, i=1,8,1)/)

       do ii=1,m
          do i=1,29
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
      do i=1,29
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

      use m_string, only: to_string
      use m_numerical_tools, only: cubspl,ppvalu
	 use global_variables
        use m_user_input

      implicit none
	  logical,intent(in)::EDX

      integer(4) i,ii,iii,j,kval,m
      integer(4) nchoices,ZZ

      real(fp_kind),allocatable:: DE(:)
      real(fp_kind) eels_inner,eels_outer
      character(2) shell_name_EELS(3)
      logical,allocatable::choices(:)
      logical::k_shell,l_shell

      shell_name_EELS = ['1s','2s','2p']

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
		  call make_detector_mask(eels_inner,eels_outer,eels_correction_detector)     
	 endif  
!Count available orbitals
	ii=0
	do i = 1, nt
		if( 5<ATF(1,i) .and. ATF(1,i)<51) ii = ii +1 !K shell
		if(19<ATF(1,i) .and. ATF(1,i)<87) then
			if(EDX) ii = ii +1 !L shell
			if(.not.EDX) ii = ii+2 !2s and 2p orbitals
		endif
	enddo
	nchoices = ii
	allocate(choices(ii),DE(ii))
	DE = 0
    choices = .false.
	
	kval = -1
	do while ((kval.ne.0).or.all(.not.choices))
100   format(/,' Ionization choices',/,/,'Index  Atom| Z  |',a,'|Included(y/n)'/&
			&,'------------------------------------------------')
	if(EDX) write(*,100) ' shell '
    if(.not.EDX) write(*,100) 'orbital|Window (eV)|'
110 format(1x,'<',i2,'>',2x,a2,2x,'|',1x,i2,1x,'|',2x,a2,3x,'|',a3,6x,a1,6x)
111 format(1x,'<',i2,'>',6x,      '|'4x,        '|',2x,a2,3x,'|',a3,6x,a1,6x)
120 format(' < 0> continue')

	ii=1
      do i = 1, nt
		ZZ = ATF(1,i)
		K_shell = 5<ATF(1,i) .and. ATF(1,i)<51
		L_shell = 19<ATF(1,i) .and. ATF(1,i)<87
        if(K_shell) then
			if(EDX) write(*,110) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),'K',' ',logical_to_yn(choices(ii))
			if(.not.EDX) write(*,110) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),'1s',to_string(DE(ii))//'          ',logical_to_yn(choices(ii))
			ii = ii +1
		endif	
		!With EDX only the L shell as a whole is offered
		if(L_shell.and.EDX) then
			if(.not.K_shell) write(*,110) ii,trim(adjustl(substance_atom_types(i))),ZZ,'L',' ',logical_to_yn(choices(ii))
			if(K_shell) write(*,111) ii,'L',' ',logical_to_yn(choices(ii))
			ii=ii+1
		!with EELS individual orbitals can be ionized
		elseif(L_shell.and.(.not.EDX)) then
			if(.not.K_shell) write(*,110) ii,trim(adjustl(substance_atom_types(i))),ZZ,'2s',to_string(DE(ii))//'          ',logical_to_yn(choices(ii))
			if(K_shell) write(*,111) ii,'2s',to_string(DE(ii))//'          ',logical_to_yn(choices(ii))
			ii = ii+1
			write(*,111) ii,'2p',to_string(DE(ii))//'          ',logical_to_yn(choices(ii))
			ii = ii+1
		endif
		
      enddo
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
	if (EDX) m =1
	if (.not.EDX) m = 5
	allocate(ionization_mu(nopiy,nopix,num_ionizations),atm_indices(num_ionizations),Ion_description(num_ionizations))
	!Now read in EELS or EDX parameters
	ii=1
    iii=1
	do i = 1, nt
		ZZ = ATF(1,i)
		!K-shell
		if(5<ATF(1,i) .and. ATF(1,i)<51) then
			if(choices(ii)) then
                
                atm_indices(iii) = i
                if(EDX) Ion_description(iii) = 'K'
                if(.not.EDX) Ion_description(iii) = '1s'
                ionization_mu(:,:,iii) = make_fz_EELS_EDX('1s',zz,DE(ii),EDX)* atf(2,i)*fz_DWF(:,:,i)
                iii= iii+1
            endif
			ii=ii+1
		endif
		
		if(19<ATF(1,i) .and. ATF(1,i)<87) then
			if(EDX) then
				if(choices(ii))	then
                    atm_indices(iii) = i
                    Ion_description(iii) = 'L'
					ionization_mu(:,:,iii) = make_fz_EELS_EDX('2s',zz,DE(ii),EDX)* atf(2,i)*fz_DWF(:,:,i)&
										   &+ make_fz_EELS_EDX('2p',zz,DE(ii),EDX)* atf(2,i)*fz_DWF(:,:,i)
                    iii= iii+1
                endif
                ii =ii +1
            else
                !loop over 2s and 2p
                do j=2,3
				    if(choices(ii))	then
                        atm_indices(iii) = i
                        Ion_description(iii) = shell_name_EELS(j)
                        ionization_mu(:,:,iii) = make_fz_EELS_EDX(shell_name_EELS(j),zz,DE(ii),EDX)* atf(2,i)*fz_DWF(:,:,i)
                        iii = iii +1
                    endif
				    ii=ii+1
                enddo
			endif
		endif
		
	enddo



      end subroutine
      
    !Subrotuine to make the Fz_mu needs to have prefactors accounted for (volume fo the unit cell etc.)
    !needs to be multiplied by the DWF for the pertinent atom type
    function make_fz_EELS_EDX(orbital,zz,DE,EDX) result(fz_mu)
	use m_precision
    use global_variables
	use m_numerical_tools, only: cubspl,ppvalu
    use m_crystallography,only:trimr
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
        if (sval.le.20.0_fp_kind) fz_mu(i,j) = cmplx(ppvalu(svals,EELS_EDX_bscoef(:,:),28,4,sval,0),0.0_fp_kind ,fp_kind) / (tp * ak1) !multiply by fractional occupancy 
    enddo;enddo
    !!$OMP END DO
	!!$OMP END PARALLEL

    return
    end function
      
    subroutine make_fz_adf()
	use m_precision
    use global_variables
    use m_absorption
    !use m_potential, only: make_g_vec_array
    use m_crystallography, only: trimr
	use m_numerical_tools, only: cubspl,ppvalu
	implicit none
    
    !dummy variables
    integer(4) k,i,j
    real(fp_kind) temp,g_vec_array(3,nopiy,nopix)
    real(fp_kind) xkstep,xdata(max_int)
    real(fp_kind) adfbrcoeff_(4,max_int,nt)

    
    write(*,*) 'Making the ADF inelastic scattering factor grid, please wait...'
    write(*,*) 
    
    if(allocated(fz_adf)) deallocate(fz_adf)
    allocate(fz_adf(nopiy,nopix,nt))

    !Set up interpolation routines
    xkstep = delta_kstep
    do i=1,max_int
          xdata(i) = (i-1)* xkstep
    enddo
    do k = 1,nt
		  adfbrcoeff_(1,:,k)= adfbr(:,k)
		  call cubspl(xdata,adfbrcoeff_(:,:,k), max_int, 0, 0 )
    enddo
    
    call make_g_vec_array(g_vec_array,ifactory,ifactorx)
    !!$OMP PARALLEL PRIVATE(i, j, k,m2,m1,sky,skx,m,g2,s2,temp)
    !!$OMP DO
    !calculate adf form factor for each pixel and each atom type in the supercell
	do i=1, nopiy; do j=1, nopix; do k = 1, nt
		temp =  ppvalu(xdata,adfbrcoeff_(:,:,k),max_int-1,4,trimr(g_vec_array(:,i,j),ss) / 2.0_fp_kind,0)
        fz_adf(i,j,k) = cmplx(temp,0.0_fp_kind,fp_kind)* atf(2,k) !multiply by fractional occupancy   
    enddo;enddo;enddo
    !!$OMP END DO
	!!$OMP END PARALLEL
    
      fz_adf = 2.0_fp_kind*tp*fz_adf

    return
    end subroutine 

      !--------------------------------------------------------------------------------------
      !   make_mu_matrix() makes the mu matrices for each HOLZ slice
      !   subroutine to take the unit cell input, 
      !   and slice based on holz
      subroutine make_local_inelastic_potentials()
      
      use m_slicing
      use m_string, only: to_string
      use global_variables, only:adf,ionization,nopiy,nopix
      
      implicit none
      
      integer(4)   i,j
      real(fp_kind) :: potential_matrix_complex(nopiy,nopix)

       
      !write(*,*) 'This will setup that ADF mu`s please enter the same slicing as before'
      write(6,134)
134   format(/,' Calculating effective inelastic potentials.',/)

      if(allocated(adf_potential)) deallocate(adf_potential)
      if(allocated(ionization_potential)) deallocate(ionization_potential)

      allocate(adf_potential(nopiy,nopix,n_slices))            !the adf potential
      allocate(ionization_potential(nopiy,nopix,num_ionizations,n_slices))     !the ionization potential
      
      if(adf) call make_fz_adf()
      
      do j = 1, n_slices
	      
	      !calculate the ionization potential
	      if(ionization) then
				do i=1,num_ionizations
					ionization_potential(:,:,i,j)= make_ion_potential(ionization_mu(:,:,i),tau_slice(:,atm_indices(i),:,j),nat_slice(atm_indices(i),j),ss_slice(7,j))
				enddo
               
	      endif  
	      !calculate the ADF potential    
            if(adf) then
                  call make_adf_potential(potential_matrix_complex,fz_adf,tau_slice(:,:,:,j),nat_slice(:,j),ss_slice(7,j))
                  adf_potential(:,:,j)= potential_matrix_complex
            endif
      enddo	!ends loop over the number of potential subslices
      
      end subroutine
      
      
      
    !--------------------------------------------------------------------------------------
    subroutine make_adf_potential(slice_potential_out,scattering_factor,tau_ss_in,nat_layer,volume)
    use global_variables
    use m_precision
	use cufft_wrapper
    use m_slicing, only: maxnat_slice
!    use m_potential!, only: make_site_factor_cuda, make_site_factor_hybrid
    
    implicit none
    
    integer(4) :: nat_layer(nt)
    complex(fp_kind),dimension(nopiy,nopix) :: slice_potential,site_term 
    complex(fp_kind),dimension(nopiy,nopix,nt) ::  scattering_factor
    real(fp_kind),intent(out)::slice_potential_out(nopiy,nopix)
    real(fp_kind) :: tau_ss_in(3,nt,maxnat_slice)  !just changed this from nat_slice?
    real(fp_kind) :: volume,V_corr
    integer(4) :: i
    
    procedure(make_site_factor_generic),pointer :: make_site_factor

#ifdef GPU
    make_site_factor => make_site_factor_cuda
#else
    make_site_factor => make_site_factor_matmul
#endif      
    
    slice_potential = 0.0_fp_kind
    V_corr = ss(7)/volume
    do i = 1, nt
        if (nat_layer(i)==0) cycle
        
        if (high_accuracy) then
            call make_site_factor(site_term, tau_ss_in(:,i,1:nat_layer(i)))
            
        else
            call make_site_factor_hybrid(site_term,  tau_ss_in(:,i,1:nat_layer(i)))
            
        endif
        
        slice_potential = slice_potential + site_term*scattering_factor(:,:,i)
    enddo
    
    slice_potential = slice_potential*V_corr
        
    ! Get real space potential
    call ifft2(nopiy,nopix,slice_potential,nopiy,slice_potential,nopiy)
    slice_potential = slice_potential * sqrt(float(nopiy)*float(nopix))
    
    ! Force that imaginary term is zero
    slice_potential_out = real(slice_potential,fp_kind)
    
    end subroutine
    
    !--------------------------------------------------------------------------------------
    function make_ion_potential(scattering_factor,tau_ss_in,nat_layer,volume) result(slice_potential_out)
    use global_variables
    use m_precision
    !use CUFFT
    use cufft_wrapper
    use m_slicing, only: maxnat_slice
!    use m_potential!, only: make_site_factor_cuda, make_site_factor_hybrid
    
    implicit none
    
    integer(4),intent(in) :: nat_layer
    real(fp_kind),intent(in) :: tau_ss_in(3,nat_layer),volume
    complex(fp_kind),intent(in)::scattering_factor(nopiy,nopix)
    
    complex(fp_kind),dimension(nopiy,nopix) :: slice_potential, site_term
    real(fp_kind)::slice_potential_out(nopiy,nopix)

    procedure(make_site_factor_generic),pointer :: make_site_factor
    
#ifdef GPU
    make_site_factor => make_site_factor_cuda
#else
    make_site_factor => make_site_factor_matmul
#endif        
    
    slice_potential = 0.0_fp_kind
    
    if (nat_layer.ne.0) then
    
        if (high_accuracy) then
            call make_site_factor(site_term, tau_ss_in)
        else
            call make_site_factor_hybrid(site_term, tau_ss_in)        
        endif
    
        slice_potential = site_term*scattering_factor/Volume
    
        ! Get realspace potential
        call ifft2(nopiy,nopix,slice_potential,nopiy,slice_potential,nopiy)
        slice_potential = slice_potential*sqrt(float(nopiy*nopix))
    endif
    
    ! Force that imaginary term is zero
    slice_potential_out = real(slice_potential)
    
    end function
        
end module
