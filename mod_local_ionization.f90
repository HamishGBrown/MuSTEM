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


      module local_ionization
      use global_variables
      use m_absorption
      use m_precision
      use m_user_input
      !use CUFFT
	  use cufft_wrapper
	  !use numerical_tools!, only: cubspl,ppvalu
      implicit none
      save 

      integer(4) :: num_ionizations
      integer(4),allocatable::  atm_indices(:)
      character(2),allocatable::ion_description(:)
      logical:: EDX
      
      complex(fp_kind), allocatable :: ionization_mu(:,:,:)        !the ionization scattering factor array, calculated on the grid (supercell)
      complex(fp_kind), allocatable :: fz_adf(:,:,:)        !the adf scattering factor array, calculated on the grid (supercell)
      real(fp_kind),    allocatable :: adf_potential(:,:,:)
      real(fp_kind),    allocatable :: ionization_potential(:,:,:,:)
      real(fp_kind),    allocatable :: eels_correction_detector(:,:)
      
      contains
      
      
      
    subroutine setup_inelastic_ionization_types()
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
        
        character(2),intent(in)::shell
		integer*4,intent(in)::atno
		real(fp_kind),intent(in)::DE
		logical,intent(in):: EDX
		
		real(fp_kind)::params(29,8,5),EELS_PARAM_SET2(5,29),xdata(8),bscoef_(4,8),fdata(8),bscoef(8)
		real(fp_kind) :: dedata(5),bscoef2_(4,5)
        
		real(fp_kind),dimension(29)::p,EELS_EDX_params
        integer*4::ishell,iatom,iz,atno_check,i,ii,m
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
	end function
	
      !********************************************************************************
      !     subroutine EELS_local_potential()
      !     reads in the scattering factors and performs the interpolation
      !     necessary to accommodate arbitrary geometries and energy windows

      subroutine local_potential(EDX)

      use m_string, only: to_string
      use m_numerical_tools, only: cubspl,ppvalu
	 

      implicit none
	  logical,intent(in)::EDX

      integer(4) i,ii,iii,j,ishell,atno,kval,m
      integer(4) iz,nchoices,iok,ZZ

      real(fp_kind),allocatable:: DE(:)
      real(fp_kind) rtemp,eels_inner,eels_outer
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
		if( 5<ATF(1,i) .and. ATF(1,i)<51) ii = ii +1
		if(19<ATF(1,i) .and. ATF(1,i)<61) ii = ii +2
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
		L_shell = 19<ATF(1,i) .and. ATF(1,i)<61
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
		
		if(19<ATF(1,i) .and. ATF(1,i)<61) then
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
	implicit none
    
    character(2),intent(in)::orbital
    integer*4,intent(in)::zz
    real(fp_kind),intent(in)::DE
    logical,intent(in)::EDX
    
    complex(fp_kind):: fz_mu(nopiy,nopix)
    
    !dummy variables
    integer(4) shiftx, shifty,m1,m2,i,j,k,n,m,nmax
    real(fp_kind) s,s2, g2, sky, skx, ax, ay
    real(fp_kind) xlen,ylen,sval,tempval

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

    ax=a0(1)*float(ifactorx)/float(nopix)
    ay=a0(2)*float(ifactory)/float(nopiy)
    xlen = a0(1)*float(ifactorx)
    ylen = a0(2)*float(ifactory)
    
    shifty = (nopiy-1)/2-1
    shiftx = (nopix-1)/2-1
    !!$OMP PARALLEL PRIVATE(i, j, m2, m1, sky, skx, tempval, sval), SHARED(fz_mu) 
    !!$OMP DO
	do i=1, nopiy
	    m2 = mod( i+shifty, nopiy) - shifty -1
        sky = float(m2)/ylen 
	    do j=1, nopix
	        m1 = mod( j+shiftx, nopix) - shiftx -1
            skx = float(m1)/xlen
            sval = sqrt(sky**2.0_fp_kind+skx**2.0_fp_kind) / 2.0_fp_kind
            if (sval.le.20.0_fp_kind) then
				tempval =  ppvalu(svals,EELS_EDX_bscoef(:,:),28,4,sval,0)
            else
                tempval = 0.0_fp_kind
            endif
            fz_mu(i,j) = cmplx(tempval,0.0_fp_kind ,fp_kind) / (tp * ak1) !multiply by fractional occupancy 
           
        enddo
    enddo
    !!$OMP END DO
	!!$OMP END PARALLEL

    return
    end function
    
      
    subroutine make_fz_adf()
	use m_precision
    use global_variables
    use m_absorption
	use m_numerical_tools, only: cubspl,ppvalu
	implicit none
    
    !dummy variables
    integer(4) shiftx, shifty,m1,m2,i,j,k,n,m,nmax
    real(fp_kind) s,s2, g2, sky, skx, ax, ay
    real(fp_kind) xlen,ylen,temp
    real(fp_kind) xkstep,xdata(max_int)
    real(fp_kind) adfbrcoeff(max_int,nt),adfbrcoeff_(4,max_int,nt)

    
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
    
    ax=a0(1)*float(ifactorx)/float(nopix)
    ay=a0(2)*float(ifactory)/float(nopiy)
    xlen = a0(1)*float(ifactorx)
    ylen = a0(2)*float(ifactory)
    
    shifty = (nopiy-1)/2-1
    shiftx = (nopix-1)/2-1
    !!$OMP PARALLEL PRIVATE(i, j, k,m2,m1,sky,skx,m,g2,s2,temp)
    !!$OMP DO
	do i=1, nopiy
	    m2 = mod( i+shifty, nopiy) - shifty -1
        sky = float(m2)/ylen 
	    do j=1, nopix
	        m1 = mod( j+shiftx, nopix) - shiftx -1
            skx = float(m1)/xlen
            g2 =  sky**2.0_fp_kind+skx**2.0_fp_kind
            s2 = g2 / 4.0_fp_kind
            !calculate adf form factor each atom type in the supercell 
            do k = 1, nt
				temp =  ppvalu(xdata,adfbrcoeff_(:,:,k),max_int-1,4,sqrt(g2),0)
                fz_adf(i,j,k) = cmplx(temp,0.0_fp_kind,fp_kind)* atf(2,k) !multiply by fractional occupancy   
            enddo
        enddo
    enddo
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
      
      implicit none
      
      integer(4)   i,j,m,n,ny,nx        
      integer(4)   nmax,idum
      real(fp_kind)  ss_temp(7)
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
    use m_potential!, only: make_site_factor_cuda, make_site_factor_hybrid
    
    implicit none
    
    integer(4) :: nat_layer(nt)
    complex(fp_kind),dimension(nopiy,nopix) :: slice_potential,potential,site_term,atom_mask 
    complex(fp_kind),dimension(nopiy,nopix,nt) ::  scattering_factor
    real(fp_kind),intent(out)::slice_potential_out(nopiy,nopix)
    real(fp_kind) :: tau_ss_in(3,nt,maxnat_slice)  !just changed this from nat_slice?
    real(fp_kind) :: volume,V_corr
    integer(4) :: i,j
    
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
    use m_potential!, only: make_site_factor_cuda, make_site_factor_hybrid
    
    implicit none
    
    integer(4),intent(in) :: nat_layer
    real(fp_kind),intent(in) :: tau_ss_in(3,nat_layer),volume
    complex(fp_kind),intent(in)::scattering_factor(nopiy,nopix)
    
    complex(fp_kind),dimension(nopiy,nopix) :: slice_potential, site_term, potential
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
