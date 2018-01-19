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

      real(fp_kind) :: EELS_EDX_params(29)
      real(fp_kind),parameter :: yp1=1.0e31_fp_kind,ypn=1.0e31_fp_kind
      real(fp_kind) :: ZZ
      !interpolation data
      real(fp_kind) :: xdata(8),bscoef_(4,8),fdata(8),bscoef(8)
      real(fp_kind) :: dedata(5),bscoef2_(4,5)!,fedata(5),bscoef2(5)
      integer(4) :: kval
      integer(4) :: i_eels
      
      complex(fp_kind), allocatable :: fz_mu(:,:)        !the MU scattering factor array, calculated on the grid (supercell)
      complex(fp_kind), allocatable :: fz_adf(:,:,:)        !the adf scattering factor array, calculated on the grid (supercell)
      real(fp_kind), allocatable :: adf_potential(:,:,:)
      real(fp_kind), allocatable :: ionization_potential(:,:,:)
      real(fp_kind), allocatable :: eels_correction_detector(:,:)
      
      contains
      
      
      
    subroutine setup_inelastic_ionization_types()
    
        implicit none
      
        write(*,*) '|-----------------------|'
        write(*,*) '|      Ionization       |'
        write(*,*) '|-----------------------|'
        write(*,*)
        
10      write(*,*) '<1> EELS'
        write(*,*) '<2> EDX'
        write(*,*) '<3> No ionization'
        call get_input('Ionization choice', i_eels)
        write(*,*)
      
        if(i_eels.eq.1) then
            ionization = .true.
            EELS = .true.
            call EELS_local_potential()
            
        elseif(i_eels.eq.2) then
            ionization = .true.
            call EDX_local_potential()
            
        elseif(i_eels.eq.3) then
            ionization = .false.
            
        else
            goto 10
        
        endif

      end subroutine
      
      
      
      !********************************************************************************
      !     subroutine EELS_local_potential()
      !     reads in the scattering factors and performs the interpolation
      !     necessary to accomodate arbitrary geometries and energy windows

      subroutine EELS_local_potential()

      use m_string, only: to_string
      use m_numerical_tools, only: cubspl,ppvalu
	 

      implicit none

      integer(4) i,ii,ishell,atno,atno_check
      integer(4) iz
      integer(4) iok
      real(fp_kind) EELS_param_set(29,8,5)
      real(fp_kind) EELS_param_set2(29,5)
      real(fp_kind) rtemp,DE
      real(fp_kind) eels_inner,eels_outer
      character(10) junk
      character(1) cjunk1,cjunk2,cjunk3
      character(2) shell_name_EELS(3)

      shell_name_EELS(1) = '1s'
      shell_name_EELS(2) = '2s'
      shell_name_EELS(3) = '2p'

      shell_name_EELS(1) = '1s'
      shell_name_EELS(2) = '2s'
      shell_name_EELS(3) = '2p'

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

100   write(*,*) 'Select the atom to be ionized:'
      
      do i = 1, nt
        write(*,102) to_string(i), substance_atom_types(i)
102     format(1x, '<', a, '> ', a)
      enddo
      
      call get_input('Index for ionized atom type', kval)
      write(*,*)
      
      if (kval.lt.1 .or. kval.gt.nt) goto 100
      
	 ZZ = ATF(1,kval)
      atno = nint(zz)

      write(6,131)
  131 format(1x, 'Choose the target electron shell: ', /, &
            &1x, '<1> 1s', /, &
            &1x, '<2> 2s', /, &
            &1x, '<3> 2p')
      call get_input('<1> 1s, <2> 2s, <3> 2p',ishell)
      write(*,*)
      
      if(((ishell.eq.2).or.(ishell.eq.3)).and.((atno.lt.20).or.(atno.gt.86))) then
            write(6,141)
  141       format(/,' L-shell only supported for 20 <= Z <= 86 ',/,' Please select a different element.',/)
            goto 100
      elseif((ishell.eq.1).and.((atno.lt.6).or.(atno.gt.50))) then
            write(6,151)
  151       format(/,' K-shell only supported for 6 <= Z <= 50 ',/,' Please select a different element.',/)
            goto 100
      endif

      write(6,161) trim(substance_atom_types(kval)), to_string(atno), shell_name_EELS(ishell)
  161 format(1x, 'Ionizations occur on ', a, ' atoms with Z = ', a, /, &
            &1x, 'Shell for ionization: ', a2, /, &
            &1x, '<1> Continue', /, &
            &1x, '<2> Change')
      call get_input('<1> OK, <2> change',iok)
      write(*,*)
      
      if(iok.ne.1) goto 100

      !open the pertinent data files
      if (ishell.eq.1) then
         open(unit=16,file='ionization_data\EELS_EDX_1s.dat',status='old',err=970)
      elseif (ishell.eq.2) then
         open(unit=16,file='ionization_data\EELS_EDX_2s.dat',status='old',err=980)
      elseif (ishell.eq.3) then
         open(unit=16,file='ionization_data\EELS_EDX_2p.dat',status='old',err=990)
      else
         write(*,*) 'Somehow you are in the wrong place, putting you somewhere sensible.'
         goto 100
      endif
 
      !readin the pertinent data
      if ((ishell.eq.2).or.(ishell.eq.3)) then
          do iz = 1,58*(atno-20)
               read(16,*) junk
          enddo
      else
          do iz = 1,58*(atno-6)
               read(16,*) junk
          enddo
      endif

      read(16,*) cjunk1,cjunk2,cjunk3,atno_check
      if (atno_check.ne.atno) then
          write(6,201)
  201     format(/,' Something has gone wrong.  Debugging required.',/)
          stop
      endif

      do i=1,8
          read(16,*) junk ! E=xx kev header
          do ii=1,5
               read(16,*) EELS_param_set(1:29,i,ii)
          enddo
          read(16,*) junk ! Don't read in EDX params
      enddo

        write(6,205)
  205   format(' Successfully read EELS electron scattering factors.',/)

      !Interpolate to accelerating voltage used
      !data in files is in steps of 50 keV with 8 points
      !this is stored in xdata
       do i=1,8
          xdata(i) = 50.0_fp_kind*i
       enddo

       do ii=1,5
          do i=1,29
			 bscoef_(1,:) = EELS_param_set(i,1:8,ii)
			 call cubspl ( xdata, bscoef_(:,:), 8, 0, 0 )
			 EELS_param_set2(i,ii) = ppvalu(xdata,bscoef_(:,:),7,4,ekv,0)
          enddo
       enddo

      ! contained within EELS_param_set2(i,ii) is the 29 data points (first index) intepolated
      ! to the correct incident energy
      ! there are 5 rows

      ! Interpolate to energy window desired
      dedata(1) = 1.0_fp_kind
      dedata(2) = 10.0_fp_kind
      dedata(3) = 25.0_fp_kind
      dedata(4) = 50.0_fp_kind
      dedata(5) = 100.0_fp_kind

  210 write(6,211)
  211 format(1x, 'Enter EELS energy window above threshold in eV (between 1 and 100 ev):')
      call get_input('Energy window', DE)
      write(*,*)
      
      !check to make sure the energy window is between 1-100 ev
      if ((DE.gt.100.0_fp_kind).or.(DE.lt.1.0_fp_kind)) goto 210

      !f(s)/DE is mostly flat and interpolates more simply
      do i=1,29
			do ii=1,5
				bscoef2_(1,ii) = EELS_param_set2(i,ii) / dedata(ii)
			enddo
			call cubspl ( dedata, bscoef2_(:,:), 5, 0, 0 )
			EELS_EDX_params(i) = DE*ppvalu(dedata,bscoef2_(:,:),4,4,DE,0)

      enddo

      goto 999
      !file open error messages
  970 close(16)
      write(6,971)
      stop
  980 close(16)
      write(6,981)
      stop
  990 close(16)
      write(6,991)
      stop

  971 format(/,' Cannot access data file EELS_EDX_1s.dat.')
  981 format(/,' Cannot access data file EELS_EDX_2s.dat.')
  991 format(/,' Cannot access data file EELS_EDX_2p.dat.')
      
  999 continue
      return

      end subroutine
      
      
      
      !********************************************************************************
      !     subroutine EDX_local_potential()
      !     reads in the scattering factors and performs the interpolation
      !     necessary to accomodate arbitrary geometries

      subroutine EDX_local_potential()
      
      use m_string, only: to_string
	  use m_numerical_tools, only: cubspl,ppvalu
      
      implicit none
      
      integer(4) i,ii,ishell,atno,atno_check,iz
      integer(4) iok
      real(fp_kind) rtemp,DE
      real(fp_kind) EDX_param_set(29,8)
      real(fp_kind) EDX_params_temp1(29),EDX_params_temp2(29)
      character(1)  shell_name_EDX(2)
      character(10) junk
      character(1)  cjunk1,cjunk2,cjunk3

      shell_name_EDX(1) = 'K'
      shell_name_EDX(2) = 'L'
      
100   write(*,*) 'Select the atom to be ionized:'
      
      do i = 1, nt
        write(*,102) to_string(i), substance_atom_types(i)
102     format(1x, '<', a, '> ', a)
      enddo
      
      call get_input('Index for ionized atom type', kval)
      write(*,*)
      
      if (kval.lt.1 .or. kval.gt.nt) goto 100
      
	ZZ = ATF(1,kval)
      atno = nint(zz)
      
      write(6,221)
  221 format(1x, 'Choose the target electron shell: ', /, &
            &1x, '<1> K-shell ionization', /, &
            &1x, '<2> L-shell ionization')
      call get_input('<1> K-shell <2> L-shell', ishell)
      write(*,*)
      
      if((ishell.eq.2).and.((atno.lt.20).or.(atno.gt.86))) then
            write(6,141)
  141       format(/,' L-shell only supported for 20 <= Z <= 60 ',/,' Please select a different element.',/)
            goto 100
      elseif((ishell.eq.1).and.((atno.lt.6).or.(atno.gt.50))) then
            write(6,151)
  151       format(/,' K-shell only supported for 6 <= Z <= 50 ',/,' Please select a different element.',/)
            goto 100
      endif

      write(6,161) trim(substance_atom_types(kval)), to_string(atno), shell_name_EDX(ishell)
  161 format(1x, 'Ionizations occur on ', a, ' atoms with Z = ', a, /, &
            &1x, 'Shell for ionization: ', a2, /, &
            &1x, '<1> Continue', /, &
            &1x, '<2> Change')
      call get_input('<1> OK, <2> change',iok)
      write(*,*)
      
      if (iok.ne.1) go to 100

      !EDX we add the contribution from the L shell electrons as the 
      !EDX spectromoters can not discriminate between edges
      if (ishell.eq.2) then
         open(unit=16,file='ionization_data\EELS_EDX_2s.dat',status='old',err=240)
         open(unit=17,file='ionization_data\EELS_EDX_2p.dat',status='old',err=250)
      elseif(ishell.eq.1) then
         open(unit=16,file='ionization_data\EELS_EDX_1s.dat',status='old',err=230)
      endif

      !read data from file
      !cycle through the header information
      if (ishell.eq.2) then
          do iz = 1,58*(atno-20)
               read(16,*) junk     !2s file
               read(17,*) junk     !2p file
	     enddo
      elseif(ishell.eq.1) then
          do iz = 1,58*(atno-6)
               read(16,*) junk     !1s file
	     enddo
      endif

      read(16,*) cjunk1,cjunk2,cjunk3,atno_check
      if (atno_check.ne.atno) then
           write(6,201)
           stop
      endif

      if (ishell.eq.2) then
          read(17,*) cjunk1,cjunk2,cjunk3,atno_check
	     if (atno_check.ne.atno) then
               write(6,201)
  201          format(/,' Something has gone wrong.  Debugging required.',/)
          endif
      endif

      do i=1,8
          do ii=1,6 ! E=xx kev header and EELS lines (the 5 lines post the Energy label are the energy window info for EELS)
               read(16,*) junk 
               if (ishell.eq.2) then    !read in the 2p info if looking at L-shell
                    read(17,*) junk 
               endif
	    enddo
          !read in the EDX points 29 of them 
          read(16,*) EDX_params_temp1(1:29)
          if (ishell.eq.2) then
               read(17,*) EDX_params_temp2(1:29)  !remember to add the 2p contribution
          elseif(ishell.eq.1) then
               EDX_params_temp2(1:29) = 0.0_fp_kind    !if we only have the K shell set the contribution to temp2=0.0
          endif

          EDX_param_set(1:29,i) = EDX_params_temp1(1:29) + EDX_params_temp2(1:29)
      enddo

      write(6,261)
  261 format(' Successfully read EDX electron scattering factors.',/)

      !Interpolate to accelerating voltage used
      do i=1,8
          xdata(i) = 50.0_fp_kind*i
      enddo

      do i=1,29
		   bscoef_(1,:) = EDX_param_set(i,1:8)
		   call cubspl ( xdata, bscoef_(:,:), 8, 0, 0 )
	       EELS_EDX_params(i) = ppvalu(xdata,bscoef_(:,:),7,4,ekv,0)
      enddo
	  
      goto 260
      !file open error messages
  230 close(16)
      write(6,231)
      stop
  240 close(16)
      write(6,241)
      stop
  250 close(17)
      write(6,251)
      stop
        
  231 format(/,' Cannot access data file EELS_EDX_1s.dat.')
  241 format(/,' Cannot access data file EELS_EDX_2s.dat.')
  251 format(/,' Cannot access data file EELS_EDX_2p.dat.')
  
  
  260 continue
      return
      end subroutine


      
      
      
    !Subrotuinet to make the Fz_mu needs to have prefactors accounted for (volume fo the unit cell etc.)
    !needs to be multiplied by the DWF for the pertinent atom type
      
    subroutine make_fz_EELS_EDX()
	use m_precision
    use global_variables
	use m_numerical_tools, only: cubspl,ppvalu
	implicit none
    
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
    
    write(*,*) 'Making the ionization inelastic scattering factor grid, please wait...'
    write(*,*)
    
    if(allocated(fz_mu)) deallocate(fz_mu)
    allocate(fz_mu(nopiy,nopix))
    	
	!pppack interpolation
	EELS_EDX_bscoef(1,:)= EELS_EDX_params
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
            fz_mu(i,j) = cmplx(tempval,0.0_fp_kind ,fp_kind)* atf(2,kval) / (tp * ak1) !multiply by fractional occupancy 
           
        enddo
    enddo
    !!$OMP END DO
	!!$OMP END PARALLEL

    return
    end subroutine
    
      
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
      allocate(ionization_potential(nopiy,nopix,n_slices))     !the ionization potential
      
      if(adf) call make_fz_adf()
      if(ionization) call make_fz_EELS_EDX()
      
      do j = 1, n_slices
	      
	      !calculate the ionization potential
	      if(ionization) then
                call make_ion_potential(potential_matrix_complex,fz_mu,tau_slice(:,:,:,j),nat_slice(:,j),ss_slice(7,j))
	            ionization_potential(:,:,j)= potential_matrix_complex
               
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
    subroutine make_ion_potential(slice_potential_out,scattering_factor,tau_ss_in,nat_layer,volume)
    use global_variables
    use m_precision
    !use CUFFT
    use cufft_wrapper
    use m_slicing, only: maxnat_slice
    use m_potential!, only: make_site_factor_cuda, make_site_factor_hybrid
    
    implicit none
    
    integer(4) :: nat_layer(nt)
    complex(fp_kind),dimension(nopiy,nopix) :: slice_potential, site_term, potential, scattering_factor 
    real(fp_kind),intent(out)::slice_potential_out(nopiy,nopix)
    real(fp_kind) :: tau_ss_in(3,nt,maxnat_slice)
    real(fp_kind) :: volume

    procedure(make_site_factor_generic),pointer :: make_site_factor
    
#ifdef GPU
    make_site_factor => make_site_factor_cuda
#else
    make_site_factor => make_site_factor_matmul
#endif        
    
    slice_potential = 0.0_fp_kind
    
    if (nat_layer(kval)==0) then
        slice_potential = 0.0_fp_kind
        return
    endif
    
    if (high_accuracy) then
        call make_site_factor(site_term, tau_ss_in(:,kval,1:nat_layer(kval)))
    else
        call make_site_factor_hybrid(site_term, tau_ss_in(:,kval,1:nat_layer(kval)))        
    endif
    
    slice_potential = site_term*scattering_factor*fz_DWF(:,:,kval)/Volume
    
    ! Get realspace potential
    call ifft2(nopiy,nopix,slice_potential,nopiy,slice_potential,nopiy)
    slice_potential = slice_potential*sqrt(float(nopiy*nopix))
    
    ! Force that imaginary term is zero
    slice_potential_out = real(slice_potential)
    
    end subroutine
    
    
    
      end module
