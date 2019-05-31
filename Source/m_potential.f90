!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. D�Alfonso, S.D. Findlay, B. D. Forbes
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

     real(fp_kind),allocatable:: eels_correction_detector(:,:)
     complex(fp_kind), allocatable :: inverse_sinc_new(:,:),ionization_mu(:,:,:)

    integer(4) :: n_qep_grates,n_qep_passes,nran ! Start of random number sequence

    logical :: phase_ramp_shift
    logical(4) :: quick_shift

    type t_slice
        !number of subslices and maximum number of atoms for a particular type in the supercell
        integer(4) :: n_slices, maxnat
        !Dimensions, atom coordinates, depths and crystallographic information for each slice of the supercell
        real(fp_kind), allocatable :: a0(:,:), tau(:,:,:,:), depths(:) , ss(:,:)
        !number of atoms for each type in the slice (supercell)
        integer(4),    allocatable :: nat(:,:)
    end type

    interface load_save_add_grates
        module procedure load_save_add_grates_qep,load_save_add_grates_abs
    end interface
    interface
        subroutine make_site_factor_generic(site_factor, tau)
            use m_precision, only: fp_kind
            complex(fp_kind),intent(out) :: site_factor(:, :)
            real(fp_kind),intent(in) :: tau(:,:)
        end subroutine make_site_factor_generic
    end interface

    logical :: save_grates = .false.
    logical :: load_grates = .false.
    character(1024) :: grates_filename

    integer*4::qep_mode
    logical :: output_probe_intensity = .false.,additional_transmission_function = .false.,pure_phase
    logical,allocatable :: output_cell_list(:)
    real(fp_kind),allocatable :: output_thickness_list(:)
    integer,allocatable :: cell_map(:)

    character*200,allocatable::amplitude_fnam(:),phase_fnam(:)


    complex(fp_kind),allocatable,dimension(:,:) :: shift_arrayx, shift_arrayy
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
        use m_string

        implicit none

        integer(4) :: i, j, k
        real(fp_kind) :: el_scat,ax,ay,g2,s2,sky,skx
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
            sinc(i,j) = cmplx((sin(tp*skx*ax)+eps)/(tp*skx*ax+eps)*((sin(tp*sky*ay)+eps)&
                                                  /(tp*sky*ay+eps)),0.0_fp_kind,fp_kind)
            inverse_sinc(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)&
                                              /(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)
            inverse_sinc_new(i,j) = cmplx((tp*skx*ax+eps)/(sin(tp*skx*ax)+eps)*(tp*sky*ay+eps)&
                                              /(sin(tp*sky*ay)+eps),0.0_fp_kind,fp_kind)
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
        use FFTW3
        use output

        implicit none

        complex(fp_kind),intent(out) :: site_factor(:,:)

        real(fp_kind),intent(in) :: tau(:,:)

        integer :: j
        integer :: yceil,yfloor,xceil,xfloor
        real(fp_kind) :: xpos, ypos, modx,mody

        site_factor = 0.0_fp_kind

        !$OMP PARALLEL PRIVATE(xpos, ypos, j,modx,mody,xceil,yceil,xfloor,yfloor)
        !$OMP DO
        do j = 1, size(tau, 2)
            !Convert fractional coordinates to pixel coordinates on grid [nopiy nopix]
            xpos = tau(1,j)*nopix
            ypos = tau(2,j)*nopiy

            !Find the 2x2 set of pixels surrounding this position, using modulo to ensure
            !that these pixels remain within the simulation grid.
            xceil  = modulo(ceiling(xpos),nopix)+1
            xfloor = modulo(floor  (xpos),nopix)+1
            yceil  = modulo(ceiling(ypos),nopiy)+1
            yfloor = modulo(floor  (ypos),nopiy)+1
            modx = mod(xpos, 1.0_fp_kind)
            mody = mod(ypos, 1.0_fp_kind)
            !Fraction of the coordinate in top right pixel
            site_factor(yceil,xceil) = site_factor(yceil,xceil)+modx*mody
            !Fraction of the coordinate in top left pixel
            site_factor(yceil,xfloor) = site_factor(yceil,xfloor)+(1_fp_kind-modx)*mody
            !Fraction of the coordinate in bottom right pixel
            site_factor(yfloor,xceil) = site_factor(yfloor,xceil)+modx*(1_fp_kind-mody)
            !Fraction of the coordinate in bottom right pixel
            site_factor(yfloor,xfloor) = site_factor(yfloor,xfloor) +(1_fp_kind - modx)*(1_fp_kind-mody)
        enddo
        !$OMP end do
        !$OMP end parallel
        site_factor = cshift(cshift(site_factor,1,dim=1),1,dim=2)
        call inplace_fft(nopiy, nopix, site_factor,norm=.false.)
        site_factor = site_factor * inverse_sinc_new

    end subroutine make_site_factor_hybrid

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
        lineno=1
        filename = 'ionization_data\EELS_EDX_'//shell//'.dat'
        open(unit=35,file=filename,status='old',err=970)

        l = len('Z = '//to_string(int(atno))//' ')
        allocate(character(l)::string)
        string = 'Z = '//to_string(int(atno))//' '
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
        use global_variables,only: ekv

        character(2),intent(in)::shell
        integer*4,intent(in)::atno
        real(fp_kind),intent(in)::DE
        logical,intent(in):: EDX

        real(fp_kind)::params(29,8,5),EELS_PARAM_SET2(5,29),xdata(8),bscoef_(4,8)
        real(fp_kind) :: dedata(5),bscoef2_(4,5),p(29),EELS_EDX_params(29)

        integer*4::i,ii,m,iz,lineno,n,mm
        character(10) junk

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
970     write(*,*) ' Cannot access data file EELS_EDX_'//shell//'.dat'
        stop
    end function

      !********************************************************************************
      !     subroutine ionization_local_potential()
      !     reads in the scattering factors and performs the interpolation
      !     necessary to accommodate arbitrary geometries and energy windows
      !********************************************************************************
      subroutine setup_ionization(EDX)

      use m_string
      use m_numerical_tools, only: cubspl,ppvalu
      use m_multislice
      use global_variables
      use m_user_input

      implicit none
      logical,intent(out)::EDX

      integer(4) i,ii,iii,j,kval,nchoices,ZZ,nshells,norbitals,i_eels

      real(fp_kind),allocatable:: DE(:)
      real(fp_kind) eels_inner,eels_outer
      character(2) shell_name_EELS(9)
      character(3) shells(9)!(9)
      character(13):: contributions(4)
      logical,allocatable::choices(:)

      shell_name_EELS = ['1s','2s','2p','3s','3p','3d','4s','4p','4d']
      shells = ['K  ','L1 ','L23','M1 ','M23','M45','N1 ','N23','N45']
      contributions = ['1s           ','2s and 2p    ','3s, 3p and 3d','4s, 4p and 4d']
      norbitals = size(shell_name_EELS)
      nshells = size(shells)


          call command_line_title_box('Ionization')
        i_eels = 0
        do while(i_eels<1.or.i_eels>2)
            write(*,*) char(10),' <1> EELS',char(10),' <2> EDX',char(10)
            call get_input('Ionization choice', i_eels)
        enddo

        EDX = i_eels.eq.2

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
                if(EDX) write(*,110) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),shell_name_EELS(j)&
                                                                    ,shells(j),logical_to_yn(choices(ii))
                if(.not.EDX) write(*,111) ii,trim(adjustl(substance_atom_types(i))),int(ZZ),shell_name_EELS(j)&
                                                             ,shells(j),DE(ii),logical_to_yn(choices(ii))
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
        if (sval.le.20.0_fp_kind) fz_mu(i,j) = cmplx(ppvalu(svals,EELS_EDX_bscoef(:,:),28,4,sval,0)&
                                                                ,0.0_fp_kind ,fp_kind) / (tp * ak1)
    enddo;enddo
    !!$OMP END DO
    !!$OMP END PARALLEL

    return
    end function

      !--------------------------------------------------------------------------------------
      !   make_mu_matrix() makes the mu matrices for each HOLZ slice
      !   subroutine to take the unit cell input,
!      subroutine make_local_inelastic_potentials(slice,ionization,adf_potential,ionization_potential)
!      !   and slice based on holz
!
!      use global_variables
!      use m_absorption
!      use m_string
!      use output
!
!      implicit none
!      class(t_slice),intent(in)::slice
!      logical,intent(in)::ionization
!      integer(4)   i,j,k,nat_
!      real(fp_kind) :: vol,adf_potential(nopiy,nopix,n_slices,ndet),ionization_potential(nopiy,nopix,num_ionizations,n_slices)
!      real(8)::thmin,thmax
!      complex(fp_kind)::fz_adf(nopiy,nopix,nt,ndet)
!
!      write(6,134)
!134   format(/,' Calculating effective inelastic potentials.',/)
!
!      do k=1,ndet/nseg
!          thmin =  atan(inner((k-1)/nseg+1)/ak)
!          thmax =  atan(outer((k-1)/nseg+1)/ak)
!          !Note that the absorptive calculations do not take into account the directionality of inelastic scattering, the absorptive scattering
!          !factors are assumed isotropic and this is only an approximation for inelastic scattering to segmented detectors
!          fz_adf(:,:,:,(k-1)*nseg+1:k*nseg) = spread(absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,nopix,&
!                                                     nt,a0,ss,atf,nat, ak, relm, orthog,thmin,thmax),dim=4,ncopies=nseg)/nseg
!      enddo
!      endif
!
!      do j = 1, n_slices
!          vol = ss_slice(7,j)
!          !calculate the ionization potential
!          if(ionization) then
!
!                do i=1,num_ionizations
!                    nat_ = slice%nat(atm_indices(i),j)
!                    ionization_potential(:,:,i,j) = real(potential_from_scattering_factors(ionization_mu(:,:,i),&
!                                                              tau_slice(:,atm_indices(i),:nat_,j),nat_,nopiy,nopix,&
!                                                              high_accuracy)/vol)
!                enddo
!
!          endif
!          !calculate the ADF potential
!            if(adf.and.complex_absorption) then
!
!                  do i=1,nt
!                      nat_ = nat_slice(i,j)
!                      do k=1,ndet
!                        adf_potential(:,:,j,k)= adf_potential(:,:,j,k) + real(potential_from_scattering_factors(fz_adf(:,:,i,k)&
!                                                        ,tau_slice(:,i,:nat_,j),nat_,nopiy,nopix,high_accuracy)/vol*ss(7)*4*pi)
!                      enddo
!                  enddo
!            endif
!      enddo !ends loop over the number of potential subslices
!
!      end subroutine

    !--------------------------------------------------------------------------------------
    function potential_from_scattering_factors(scattering_factor,atom_posn,nat_layer,nopiy,nopix,high_accuracy)&
        result(slice_potential)
    use m_precision
    use FFTW3

    implicit none

    integer(4),intent(in) :: nat_layer,nopiy,nopix
    real(fp_kind),intent(in) :: atom_posn(3,nat_layer)
    complex(fp_kind),intent(in)::scattering_factor(nopiy,nopix)
    logical,intent(in),optional::high_accuracy

    complex(fp_kind),dimension(nopiy,nopix) ::  site_term,slice_potential
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
        call inplace_ifft(nopiy,nopix,slice_potential,norm=.false.)
    endif

    end function
    function make_absorptive_grates(slice,nopiy,nopix) result(projected_potential)

        use m_precision, only: fp_kind
        use global_variables, only: ig1,ig2,ifactory,ifactorx,nt, relm, tp, ak, atf,&
                           high_accuracy,  pi, fz,fz_DWF,ss,a0,nat,orthog
        use m_absorption!, only: transf_absorptive,fz_abs
        use m_string, only: to_string
        use output

        implicit none

        integer*4,intent(in)::nopiy,nopix
        type(t_slice),intent(in)::slice
        complex(fp_kind)::projected_potential(nopiy,nopix,slice%n_slices)

        integer(4) :: j, m, nat_layer
        real(fp_kind) :: ccd_slice,V_corr
        complex(fp_kind),dimension(nopiy,nopix) :: effective_scat_fact
        complex(fp_kind)::fz_abs(nopiy,nopix,nt)

        real(fp_kind) :: t1, delta

        projected_potential= 0
        t1 = secnds(0.0_fp_kind)
        fz_abs=0
        if(load_grates) then
            call load_save_add_grates(projected_potential,nopiy,nopix,slice)
        else
            if(include_absorption) fz_abs = absorptive_scattering_factors(ig1,ig2,ifactory,ifactorx,nopiy,&
                                            nopix,nt,a0,ss,atf,nat, ak, relm, orthog, 0.0_8, 4.0d0*atan(1.0d0))*2*ak

            do j = 1, slice%n_slices


                ccd_slice = relm / (tp * ak * slice%ss(7,j))
                V_corr = ss(7)/slice%ss(7,j)
                do m=1,nt
                    nat_layer = slice%nat(m,j)
                    effective_scat_fact = CCD_slice*fz(:,:,m)*fz_DWF(:,:,m)+cmplx(0,1)*fz_abs(:,:,m)*V_corr
                    projected_potential(:,:,m) = projected_potential(:,:,m)+potential_from_scattering_factors(effective_scat_fact,&
                                                 slice%tau(:,m,1:nat_layer,j),nat_layer,nopiy,nopix,high_accuracy)

                enddo
                write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', &
                            to_string(j), '/', to_string(slice%n_slices), '...'
                write(6,199) to_string(sum(slice%nat(:,j)))
        199     format(1x, 'Number of atoms in this slice: ', a, /)


            enddo ! End loop over slices
            call load_save_add_grates(projected_potential,nopiy,nopix,slice)
        endif
        delta = secnds(t1)

        if(timing) then
            write(*,*) 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
            write(*,*)
            open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
            write(9834, '(a, g9.5, a, /)') 'The calculation of transmission functions for the absorptive model took '&
                                                                                                  , delta, 'seconds.'
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
                random = real(ran1(idum),fp_kind)
            enddo

        end function seed_rng


    subroutine make_propagator(nopiy,nopix,prop,dz,ak1,ss,claue,ifactorx,ifactory)

        use m_precision, only: fp_kind
        use m_crystallography, only: trimr,make_g_vec_array

        implicit none

        integer(4) :: nopiy,nopix
        complex(fp_kind) :: prop(nopiy,nopix)
        real(fp_kind) :: ak1, ss(7), claue(3), dz, g_vec_array(3,nopiy,nopix)
        integer(4) :: ifactorx, ifactory

        real(fp_kind),parameter :: pi = atan(1.0_fp_kind)*4.0_fp_kind
        integer(4) :: ny, nx


        call make_g_vec_array(g_vec_array,ifactory,ifactorx)

        do ny = 1, nopiy;do nx = 1, nopix
            prop(ny,nx) = exp(cmplx(0.0d0, -pi*dz*trimr(g_vec_array(:,ny,nx)-claue,ss)**2/ak1, fp_kind ))
        enddo;enddo

    end subroutine

    function make_qep_grates(slice,nopiy,nopix,n_qep_grates,idum) result(projected_potential)

        use m_precision, only: fp_kind
        use global_variables, only: nt, relm, tp, ak, atf, high_accuracy, pi, fz
        use m_string, only: to_string
        use output, only: output_prefix,timing,binary_in
        use m_numerical_tools, only: displace

        implicit none

        class(t_slice),intent(in)::slice
        integer(4),intent(inout) :: idum
        integer*4,intent(in)::nopiy,nopix,n_qep_grates
        complex(fp_kind) :: projected_potential(nopiy,nopix,n_qep_grates,slice%n_slices)
        integer(4), allocatable :: handled(:,:)
        integer(4):: save_list(2,nt),match_count, i, j, m, n,ii,jj,jjj,kk,iii
        real(fp_kind) :: tau_holder(3),tau_holder2(3),ccd_slice,ums
        real(fp_kind) :: mod_tau(3,nt,slice%maxnat,slice%n_slices,n_qep_grates),t1, delta
        logical::fracocc

        procedure(make_site_factor_generic),pointer :: make_site_factor


        !Search for fractional occupancy, fractional occupancy would require that we shift
        !atoms sharing the same site (but with different fractional occupancy) by the same
        !amount
         fracocc = any(atf(2,:).lt.0.99d0)

        t1 = secnds(0.0_fp_kind)

        do j = 1, slice%n_slices
            write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', &
                                           to_string(j), '/', to_string(slice%n_slices), '...'

            write(6,199) to_string(sum(slice%nat(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a)

            ccd_slice = relm / (tp * ak * slice%ss(7,j))

            do i = 1, n_qep_grates
    200         format(a1, 1x, i3, '/', i3)
                write(6,200, advance='no') achar(13), i, n_qep_grates

                ! Randomly displace the atoms
                if (.not.fracocc) then
                do m = 1, nt
                    do n = 1,slice%nat(m,j)
                        call displace(slice%tau(1:3,m,n,j),mod_tau(1:3,m,n,j,i),sqrt(atf(3,m)),slice%a0,idum)
                    enddo
                enddo
                else
                    allocate( handled(nt,slice%maxnat) )
                    handled = 0
                    do ii=1, nt
                     do jj = 1, slice%nat(ii,j)
                         if (handled(ii,jj).eq.1) cycle
                         tau_holder(1:3) = slice%tau(1:3,ii,jj,j)

                         save_list = 0
                         match_count = 0
                         ums = atf(3,ii)
                         do iii=ii+1,nt
                         do jjj=1,slice%nat(iii,j)
                            if (same_site(tau_holder,slice%tau(1:3,iii,jjj,j))) then
                               match_count = match_count+1
                               save_list(:,match_count)=[iii,jjj]
                               ums = ums + atf(3,iii)
                               cycle
                            endif
                         enddo
                         enddo

                         ums = ums / real(match_count+1,kind=fp_kind)
                       call displace(tau_holder(1:3),tau_holder2(1:3),sqrt(ums),slice%a0,idum)
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
                    projected_potential(:,:,i,j) = projected_potential(:,:,i,j)&
                                                  +real(potential_from_scattering_factors(CCD_slice*fz(:,:,m)&
                                                 &,mod_tau(:,m,1:slice%nat(m,j),j,i),slice%nat(m,j),nopiy,nopix,high_accuracy))
                enddo
            enddo ! End loop over grates

            write(*,*)

        enddo ! End loop over slices

        delta = secnds(t1)

        write(*,*) 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
        write(*,*)

        if(timing) then
            open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
            write(9834, '(a, g9.5, a, /)') 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
            close(9834)
        endif

    end function make_qep_grates

    logical(4) function same_site(site1,site2)

      implicit none

      real(fp_kind) site1(3),site2(3)
      real(fp_kind) tol

      tol = 1.0e-6_fp_kind
      same_site = all(abs(site1-site2).lt.tol)

      return
      end function


    subroutine prompt_output_probe_intensity(slice)

        use m_user_input, only: get_input
        use global_variables, only: thickness, nopiy, nopix
        use output, only: output_prefix,split_filepath
        use m_string, only: to_string,command_line_title_box

        implicit none

        type(t_slice),intent(in)::slice
        integer :: i_output
        real(fp_kind) :: thickness_interval
        character(1024)::dir,fnam
        logical::invalid_thickness

        call command_line_title_box('Output probe intensity')
        write(*,*) 'The probe intensity as a function of thickness'
        write(*,*) 'can be outputted to file at each probe position.'
        write(*,*) 'The user is advised that the outputted dataset'
        write(*,*) 'may be very large.'
        write(*,*)
        write(*,*) '<0> Proceed'
        write(*,*) '<1> Output probe intensity'
        call get_input('<0> Proceed <1> Output probe intensity', i_output)
        write(*,*)

        output_probe_intensity = i_output ==1

        if(output_probe_intensity) then

            invalid_thickness = .true.
            do while(invalid_thickness)
                write(*,*) 'At what thickness interval (in Angstroms)',char(10),' should intensities be outputted?'
                call get_input('At what thickness interval should intensities be outputted?', thickness_interval)
                write(*,*)
                invalid_thickness = thickness_interval.le.0.0_fp_kind .or. thickness_interval.gt.thickness
                if(invalid_thickness) write(*,*) 'ERROR: invalid thickness.'
            enddo

                call split_filepath(output_prefix,dir,fnam)
                call system('mkdir '//trim(adjustl(dir))//'\Probe_intensity')
                call generate_cell_list(thickness_interval,slice)
                call write_thicknesss_to_file

                write(*,*) 'The probe intensities will be written to the files'
                write(*,*)
                write(*,*) '  '//trim(adjustl(dir))//'\Probe_intensity' // trim(adjustl(fnam)) // '_ProbeIntensity*.bin'
                write(*,*)
                if (fp_kind.eq.4) write(*,*) 'as 32-bit big-endian floats.'
                if (fp_kind.eq.8) write(*,*) 'as 64-bit big-endian floats.'

                write(*,*) 'Each file contains a sequence of ' // to_string(nopiy) // 'x' // to_string(nopix) // ' arrays.'
                write(*,*)
        end if

    end subroutine



    subroutine generate_cell_list(thickness_interval,slice)

        use global_variables, only: ncells, a0

        implicit none

        type(t_slice)::slice
        real(fp_kind) :: thickness_interval

        integer :: count,i,j
        real(fp_kind) :: t,tout

        if(allocated(output_cell_list)) deallocate(output_cell_list)
        allocate(output_cell_list(maxval(ncells)*slice%n_slices))
        output_cell_list = .false.

        if(allocated(cell_map)) deallocate(cell_map)
        allocate(cell_map(maxval(ncells)*slice%n_slices))
        cell_map = 0

        t = 0.0_fp_kind
        tout = thickness_interval

        count = 0

        do i= 1,maxval(ncells)
        do j=1,slice%n_slices
            if((t+slice%depths(j+1)*a0(3).gt.tout)) then
                output_cell_list((i-1)*slice%n_slices+j) =.true.
                count = count + 1
                tout = (floor((t+slice%depths(j+1)*a0(3))/thickness_interval)+1)*thickness_interval
            endif

        enddo
        t = t+a0(3)
        enddo
        if(allocated(output_thickness_list)) deallocate(output_thickness_list)
        allocate(output_thickness_list(count))

        t = 0.0_fp_kind
        tout = thickness_interval

        count = 0
        do i= 1,maxval(ncells)
        do j=1,slice%n_slices
            if((t+slice%depths(j+1)*a0(3).gt.tout)) then
                count = count + 1
                output_thickness_list(count) = t+slice%depths(j+1)*a0(3)
                cell_map((i-1)*slice%n_slices+j) = count
                tout = (floor((t+slice%depths(j+1)*a0(3))/thickness_interval)+1)*thickness_interval
            endif
        enddo
        t = t+a0(3)
        enddo

    end subroutine



    subroutine write_thicknesss_to_file

        use output, only: output_prefix

        implicit none

        integer :: i,j
        character(1024) :: filename,dir,fnam

        j = index(output_prefix,'/',back=.true.)
        j = max(j,index(output_prefix,'\\',back=.true.))

        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'\Probe_intensity'//trim(adjustl(fnam))//'_probe_intensity_thicknesss.txt'
        else
            filename = 'Probe_intensity\\'//trim(adjustl(output_prefix))//'_probe_intensity_thicknesss.txt'
        endif

        write(*,*) 'The thicknesses at which the probe intensity'
        write(*,*) 'is being outputted have been written to'
        write(*,*)
        write(*,*) '  ' // trim(filename)
        write(*,*)

        open(unit=8734, file=filename)

        do i = 1, size(output_thickness_list)
            write(8734, *) output_thickness_list(i)
        enddo

        close(8734)



    end subroutine

    subroutine probe_intensity_to_file(probe_intensity,i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample)

    use output, only: output_prefix,quad_shift
    use global_variables, only: nopiy,nopix
    use m_string

    real(fp_kind),intent(in)::probe_intensity(nopiy,nopix,size(output_thickness_list))
    integer*4,intent(in)::i_df,ny,nx,n_qep_passes,probe_ndf,nysample,nxsample

    integer*4::j,z
    character*1024::filename,fnam,dir

        j = index(output_prefix,'/',back=.true.)
        j = max(j,index(output_prefix,'\\',back=.true.))
        z = size(output_thickness_list)

        if(j>0) then
            dir = trim(adjustl(output_prefix(:j)))
            fnam = trim(adjustl(output_prefix(j:)))
            filename = trim(adjustl(dir))//'\\Probe_intensity\\'//trim(adjustl(fnam))//'_ProbeIntensity'
        else
            filename = 'Probe_intensity\\'//trim(adjustl(output_prefix))//'_ProbeIntensity'
        endif

        if (probe_ndf.gt.1) filename = trim(filename) // '_df' // to_string(i_df)
        if (nysample.gt.1) filename = trim(filename) // '_ny' // to_string(ny)
        if (nxsample.gt.1) filename = trim(filename) // '_nx' // to_string(nx)
        filename = trim(filename) // '_'//to_string(nopiy)//'x'//to_string(nopix)//'x'//to_string(z)//'.bin'
        open(4985, file=filename, form='unformatted', convert='big_endian')
        do j=1,z
            write(4985) quad_shift(probe_intensity(:,:,j),nopiy,nopix)/ n_qep_passes
        enddo
        close(4985)

    end subroutine

    subroutine prompt_save_load_grates(slice)

        use m_user_input, only: get_input
        use output, only: output_prefix
        use m_string

        implicit none
        type(t_slice),intent(in)::slice
        integer :: i_save_load, i_retry,i
        logical :: exists,retry

        i_save_load = -1
        do while(i_save_load<0.or.i_save_load>2)


        call command_line_title_box('Save/load transmission functions')
        write(*,*) 'Warning: the files outputted when saving may be very large.'
        write(*,*)
        write(*,*) '<0> Proceed without saving or loading'
        write(*,*) '<1> Save transmission functions'
        write(*,*) '<2> Load transmission functions'
        write(*,*) '<3> Add additional transmission function '
        write(*,*) '    (eg. from magnetic structure) from file'
        call get_input('<0> continue <1> save <2> load', i_save_load)
        write(*,*)
323     format( '  - The xtl file',/,'  - The slicing of the unit cell',/,&
               &'  - The choice of thermal scattering model (QEP vs. absorptive)',/,&
               &'  - The tiling of the unit cell',/,'  - The number of pixels',/,&
               &'  - (For absorptive model: whether absorption is included)',/,&
               &'  - (For QEP model: the number of distinct transmission functions)',/,&
               &'  - (For QEP model: phase ramp shift choice)',/)

        select case (i_save_load)
            case (0)
                return

            case (1)
                save_grates = .true.
                grates_filename = trim(adjustl(output_prefix)) // '_transmission_functions.bin'

                write(*,*) 'The transmission functions will be saved to the file'
                write(*,*)
                write(*,*) '    ' // trim(grates_filename)
                write(*,*)
                write(*,*) 'They can be loaded for later calculations provided'
                write(*,*) 'the following parameters are identical:'
                write(6,323)

            case (2)
                write(*,*) 'It is up to the user to ensure that the parameters used'
                write(*,*) 'to create the loaded transmission functions are consistent'
                write(*,*) 'with those of the current calculation:'
                write(6,323)
                retry=.true.
                do while(retry)
               write(*,*) 'Enter filename of transmission functions:'
                call get_input('filename of transmission functions', grates_filename)
                write(*,*)

                inquire(file=grates_filename, exist=exists)
                retry=.not.exists
                if (.not.exists) then
                    write(*,*) 'ERROR: cannot find this file.'
                    write(*,*) '<1> Enter again'
                    write(*,*) '<2> Proceed without loading'
                    call get_input('<1> Enter again <2> Proceed without loading', i_retry)
                    write(*,*)
                    retry = i_retry==1
                    if(.not.retry) return


                endif
                enddo
                load_grates = .true.
            case (3)
                additional_transmission_function=.true.
                pure_phase=.true.

                if(.not.allocated(amplitude_fnam)) allocate(amplitude_fnam(slice%n_slices),phase_fnam(slice%n_slices))
                do i=1,slice%n_slices
                    if(.not.pure_phase) then
                        write(*,*) char(10),' Please input filename of amplitude of additional transmission function for slice ',i
                        if(i==1) then
                            write(*,*) 'If your intended transmission function is a pure phase object (ie. the '
                            write(*,*) 'amplitude of theadditional transmission function is everywhere 1), input'
                            write(*,*) '-1 and no amplitude file willbe loaded nor will you be prompted for the '
                            write(*,*) 'amplitude for the remaining slices'
                        endif
                        call get_input('Amplitude of additional transmission function',amplitude_fnam(i))
                        pure_phase = trim(adjustl(amplitude_fnam(i)))=='-1'
                    endif
                    write(*,*) char(10),' Please input filename of phase of additional transmission function for slice ',i
                    call get_input('Phase of additional transmission function',phase_fnam(i))
                enddo



        end select
        enddo
    end subroutine

 subroutine load_save_add_grates_qep(idum,qep_grates,nopiy,nopix,n_qep_grates,slice,nt)
        use m_numerical_tools, only: gasdev
        use output
        use global_variables, only: ak,pi
        integer(4),intent(inout):: idum
        type(t_slice) :: slice
        complex(fp_kind),intent(inout)::qep_grates(nopiy,nopix,n_qep_grates,slice%n_slices)
        integer*4,intent(in)::nopiy,nopix,n_qep_grates,nt

        integer*4::j,i,m,n
        real(fp_kind)::junk

        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)


        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            read(3984) qep_grates
            close(3984)

            ! Make sure the random number sequence is as if grates were calculated
            ! So call gasdev as many times as it would usually be called
            do j = 1, slice%n_slices;do i = 1, n_qep_grates;do m = 1, nt;do n = 1, slice%nat(m,j)*2
                junk = gasdev(idum)
            enddo;enddo;enddo;enddo

            return
        endif

        if(additional_transmission_function) then
            write(*,*) 'Adding addition transmission function to file...'
            amplitude = 1
            do j=1,slice%n_slices
                if(.not.pure_phase) call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                call binary_in(nopiy,nopix,phase,phase_fnam(j))
                qep_grates(:,:,:,j) = qep_grates(:,:,:,j)+spread(transpose(phase),dim=3,ncopies=n_qep_grates)/pi/slice%a0(3,j)*ak
                if(.not.pure_phase) then
                    call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                    qep_grates(:,:,:,j) = qep_grates(:,:,:,j)&
                                        +spread(cmplx(0_fp_kind,log(transpose(amplitude))),dim=3,ncopies=n_qep_grates)
                endif
            enddo
            qep_mode=4
        endif

        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            write(3984) qep_grates*pi*spread(spread(spread(slice%a0(3,1:slice%n_slices),dim=1,ncopies=nopiy)&
                                                        &,dim=2,ncopies=nopix),dim=3,ncopies=n_qep_grates)/ak;close(3984)
        endif
    end subroutine

subroutine load_save_add_grates_abs(abs_grates,nopiy,nopix,slice)
        use global_variables, only: pi,ak
        use m_numerical_tools, only: gasdev
        use output
        type(t_slice)::slice
        complex(fp_kind),intent(inout)::abs_grates(nopiy,nopix,slice%n_slices)
        integer*4,intent(in)::nopiy,nopix


        integer*4::j

        real(fp_kind)::amplitude(nopiy,nopix),phase(nopiy,nopix)

        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            read(3984) abs_grates
            close(3984)
            return
        endif

        if(additional_transmission_function) then
            write(*,*) 'Adding addition transmission function from file...'
            amplitude = 1
            do j=1,slice%n_slices
                call binary_in(nopiy,nopix,phase,phase_fnam(j))
                abs_grates(:,:,j) = abs_grates(:,:,j)+transpose(phase)/pi/slice%a0(3,j)*ak
                if(.not.pure_phase) then
                    call binary_in(nopiy,nopix,amplitude,amplitude_fnam(j))
                    abs_grates(:,:,j) = abs_grates(:,:,j)+cmplx(0_fp_kind,log(transpose(amplitude)))
                endif
            enddo
        endif

        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...',char(10)
            open(unit=3984, file=grates_filename, form='unformatted', convert='big_endian')
            write(3984) abs_grates*pi*spread(spread(slice%a0(3,1:slice%n_slices),dim=1,ncopies=nopiy),dim=2,ncopies=nopix)/ak
            close(3984)
        endif

    end subroutine


    function make_detector(nopiy,nopix,ifactory,ifactorx,ss,kmin,kmax,phi,delphi) result(detector)
        use m_crystallography
        !makes a detector on an array nopiy x nopix with Fourier space pixel size deltaky x deltakx
        !which measures from kmin to kmax, supplying phi and delphi (detector orientation
        !and angular range in radians) will create a detector segment
        integer*4,intent(in)::nopiy,nopix,ifactory,ifactorx
        real(fp_kind),intent(in)::kmin,kmax,phi,delphi,ss(7)
        optional::phi,delphi

        real(fp_kind)::phi_(2)=0_fp_kind,g_vec_array(3,nopiy,nopix),ang,kabs,detector(nopiy,nopix),pi
        integer*4::y,x
        logical::segment,notwrapped

        pi = 4.0_fp_kind*atan(1.0_fp_kind)
        !If the phi and delphi variables are present, then make a segmented detector
        segment = present(phi).and.present(delphi)
        if(segment) then
            !Force angle to be between 0 and 2pi
            phi_= mod([phi,phi+delphi],2*pi)
            do y=1,2
                if (phi_(y)<0) phi_(y) = phi_(y) + 2*pi
            enddo
            notwrapped = phi_(2)>phi_(1)
        endif
        detector = 0

        call make_g_vec_array(g_vec_array,ifactory,ifactorx)
        detector = 0
        do y=1,nopiy;do x=1,nopix
            kabs = trimr(g_vec_array(:,y,x),ss)
            if((kabs.ge.kmin).and.(kabs.le.kmax)) then
                detector(y,x) =1
                if(segment) then
                    !Force angle to be between 0 and 2pi
                    ang = modulo(atan2(g_vec_array(1,y,x),g_vec_array(2,y,x)),2*pi)
                    if (notwrapped) then
                        if(.not.(ang>phi_(1).and.ang<=phi_(2))) detector(y,x)=0
                    else
                        if(.not.(ang>phi_(1).or.ang<=phi_(2))) detector(y,x)=0
                    endif
            endif;endif
            enddo;enddo;
    end function

    subroutine setup_qep_parameters(qep,n_qep_grates,n_qep_passes,nran,quick_shift,ifactory,ifactorx)

        use m_user_input, only: get_input
        use m_string, only: to_string,command_line_title_box

        implicit none

        integer*4,intent(out)::n_qep_grates,n_qep_passes,nran
        integer*4,intent(in)::ifactory,ifactorx
        logical,intent(in)::quick_shift,qep


        integer :: i

        if (.not.qep) then
            n_qep_grates = 1
            n_qep_passes = 1
            return
        endif

        call command_line_title_box('QEP model parameters')

        qep_mode=1
        do while(qep_mode<2.or.qep_mode>4)
        write(*,*) 'Enter the number of distinct transmission functions to calculate:'
        call get_input("Number of phase grates calculated", n_qep_grates)
        write(*,*)

            write(*,*) 'As implemented in muSTEM, the QEP multislice algorithm will randomly shift the'
            write(*,*) 'unit cells in your supercell around to effectively generate new random '
            write(*,*) 'transmission functions. This means the QEP calculation samples a larger '
            write(*,*) 'effective number of random transmission functions than would otherwise be the'
            write(*,*) 'case.'
            write(*,*)
            write(*,*) 'If your choice of unit cell tiling and pixel grid size is such that each unit '
            write(*,*) 'cell is an integer number of pixels, unit cell shifting can be achieved by '
            write(*,*) 'circular shift of the transmission function arrays in memory (quick shifting).'
            write(*,*) 'Otherwise the unit cells will have to be shifted with sub-pixel precision '
            write(*,*) 'using the Fourier shift algorithm which will impact calculation speed. '
            write(*,*)
            write(*,*) 'A third option is to not shift the unit cells around at all, but the user is'
            write(*,*) 'advised that this requires a much larger number of distinct transmission '
            write(*,*) 'functions to achieve convergence than would usually be the case.'
            write(*,*)
        if (quick_shift) then
            write(6,100) to_string(ifactorx*ifactory), to_string(n_qep_grates), to_string(ifactorx*ifactory*n_qep_grates)
100         format(' The choice of tiling and grid size permits quick shifting.',/,&
                  &' The effective number of transmission functions used in  ',/,&
                  &' calculations will be ', a, ' * ', a, ' = ', a, '.', /)
            qep_mode=2
        else
            write(6,101)
        101 format( ' Your choice of tiling and grid size does not permit quick shifting ', /, &
                    &' of the precalculated transmission functions. Shifting using the ', /, &
                    &' Fourier shift algorithm can be performed but is time consuming. ', /, &
                    &' You may wish to go back and calculate more distinct transmission ', /, &
                    &' functions, or simply proceed without using phase ramp shifting. ' /)

            i=-1
            do while(i<1.or.i>3)
                write(6,111)
            111 format(  ' <1> Go back and choose a larger number', /, &
                        &' <2> Proceed with phase ramp shifting', /, &
                        &' <3> Proceed without phase ramp shifting', / &
                        )
                call get_input("<1> choose more <2> phase ramp <3> no phase ramp", i)
                write(*,*)

                if (i.eq.2 .and. (ifactory.gt.1 .or. ifactorx.gt.1)) then
                    qep_mode=3
                    call setup_phase_ramp_shifts

                elseif (i.eq.3) then
                    qep_mode=4
                endif
            enddo
        endif
        enddo
        write(6,*) 'Enter the number of passes to perform for QEP calculation:'
        write(*,*) 'Warning: using only a single pass is usually NOT sufficient.'
        call get_input("Number of Monte Carlo calculated", n_qep_passes )
        write(*,*)

        write(6,*) 'Enter the starting position of the random number sequence:'
        call get_input("Number of ran1 discarded", nran )
        write(*,*)

        end subroutine

    subroutine setup_slicing_depths(slice)

        use m_user_input, only: get_input
        use m_precision, only: fp_kind
        use m_string,only:command_line_title_box

        implicit none

        integer :: i_slice
        type(t_slice)::slice
        call command_line_title_box('Unit cell slicing')

    22  write(6,23) char(143)
    23  format(' Do you wish to slice the unit cell in the beam direction?', /, &
              &' This may be a good idea if the unit cell is larger than 2 ', a1, /, &
              &' in the z-direction.', /, &
              &' <1> Yes <2> No ')
        call get_input('Slice unit cell <1> yes <2> no', i_slice)
        write(*,*)

        if (i_slice.eq.1) then
            call calculate_depths_slicing(slice)

        elseif (i_slice.eq.2) then
            call calculate_depths_no_slicing(slice)

        else
            goto 22

        endif

    end subroutine

    subroutine calculate_depths_no_slicing(slice)
        implicit none
        type(t_slice)::slice
        slice%n_slices = 1

        allocate(slice%depths(2))

        slice%depths(1) = 0.0_fp_kind
        slice%depths(2) = 1.0_fp_kind

    end subroutine

    subroutine calculate_depths_slicing(slice)
        use m_user_input

        implicit none
        type(t_slice)::slice
        integer(4) i, ichoice

        slice%n_slices = -1
        do while(slice%n_slices<1)
        write(*,*) 'Enter the number of slices per unit cell in the beam direction:'
        call get_input("Number of distinct potentials ", slice%n_slices)
        write(*,*)

        if (slice%n_slices.eq.1) then; call calculate_depths_no_slicing(slice)
        elseif (slice%n_slices.gt.1) then
            if (allocated(slice%depths)) deallocate(slice%depths)
            allocate(slice%depths(slice%n_slices+1))

            slice%depths(slice%n_slices+1) = 1.0_fp_kind

            write(6,10)
         10 format( ' You will now be asked to enter the depths at which slicing',/,&
                    &' will take place. These should be entered as fractions of',/,&
                    &' the unit cell. It is the front of the slice which should',/,&
                    &' be entered. (e.g. three even slices would be entered as',/,&
                    &' 0.0, 0.3333, 0.6666 ).', /)

            ichoice = 0
            do while(ichoice<1.or.ichoice>2)
            write(6,16)
    16      format(' How would you like to specify the slice depths?', /, &
                  &' <1> Manually ', /, &
                  &' <2> Automatically (uniformly spaced)')
            call get_input("<1> manual or <2> auto slicing", ichoice)
            write(*,*)

            if(ichoice.eq.1) then
                ! Manual slicing

                do i = 1, slice%n_slices
                    write(6,20,ADVANCE='NO') achar(13), i
        20          format( a1,' Enter the fractional depth of slice number ', i4, ':  ')
                    call get_input("depths", slice%depths(i))
                enddo

            elseif (ichoice.eq.2) then
                ! Automatic slicing

    21          format( ' Fractional depth of slice ', i4, ':  ', f7.4)
                do i = 1, slice%n_slices
                    slice%depths(i) = float( (i-1)) / float(slice%n_slices)
                    write(6,21) i, slice%depths(i)
                enddo
            endif
            enddo
        endif
        enddo

        write(*,*)

    end subroutine

    subroutine initialize_slicing(slice,nat,ifactory,ifactorx,nt,a0,deg,tau,nm)
        use m_crystallography, only: cryst

        implicit none
        type(t_slice)::slice
        integer*4,intent(in)::ifactory,ifactorx,nt,nat(nt),nm
        real(fp_kind),intent(in)::deg(3),a0(3),tau(3,nt,nm)
        integer :: i,j,k,mm,nn
        real(fp_kind)::diff,factorx,factory,tol

        !Make floating point copies of ifactory and ifactorx
        factorx = real(ifactorx,fp_kind)
        factory = real(ifactory,fp_kind)
        !Prompt user for slicing of unit cell

        tol=1e-4_fp_kind

        slice%maxnat = maxval(nat) * ifactory * ifactorx

        allocate(slice%nat(nt,slice%n_slices),slice%tau(3,nt,slice%maxnat,slice%n_slices),&
                &slice%a0(3,slice%n_slices),slice%ss(7,slice%n_slices))

        do j = 1, slice%n_slices

            !Calculate slice dimensions in Angstrom, this is given by the product of the unit cell dimensions
            !and the crystal tiling for x and y directions
            slice%a0(1,j) = a0(1)*ifactorx
            slice%a0(2,j) = a0(2)*ifactory

            !In the z-direction the dimensions of the slice are given by the user supplied slicing of the unit
            !cell
            if (j .eq. slice%n_slices) then
                slice%a0(3,j) = (1.0_fp_kind-slice%depths(j))*a0(3)
            else
                slice%a0(3,j) = (slice%depths(j+1)-slice%depths(j))*a0(3)
            endif

            !This routines calculates quantities used to calculate size of reciprocal lattice vectors and volume etc. for this slice
            call cryst(slice%a0(:,j),deg,slice%ss(:,j))

            !Store the slice thickness for convenience
            diff = slice%depths(j+1)-slice%depths(j)

            !Now tile out the atomic coordinates, loop over each atomic species
            do i = 1, nt
               !Zero the number of atoms in this slice
               slice%nat(i,j) = 0

               do mm = 1, ifactory
               do nn = 1, ifactorx
                  do k = 1, nat(i)

                    !Check if atomic coordinate sits within this slice, if so add to the list with appropriate adjustments
                    !to the fractional coordinate, since the effective unit cell is typically larger.
                     if( (tau(3,i,k) .lt. (slice%depths(j+1)-tol)) .and.(tau(3,i,k) .ge. (slice%depths(j)-tol)) ) then
                        slice%nat(i,j) = slice%nat(i,j) + 1
                        slice%tau(1:2,i,slice%nat(i,j),j) = (tau(1:2,i,k) + float([nn,mm]-1))/[factorx,factory]
                        slice%tau(3,i,slice%nat(i,j),j) = (tau(3,i,k)-slice%depths(j))/diff

                    !Possible bug, this bit of code will only be used if a single slice is requested and atoms have z fractional
                    !coordinate z>1-1e-4 the atom counter will not be increased by one and the atom will overwrite the previous
                    !entry.
                     elseif(diff.eq.1.0) then
                        slice%tau(1:2,i,slice%nat(i,j),j) = (tau(1:2,i,k) + float([nn,mm]-1))/[factorx,factory]
                        slice%tau(3,i,slice%nat(i,j),j) = (tau(3,i,k)-slice%depths(j))/diff
                     endif
                  enddo
                enddo
                enddo
            enddo
        enddo
    end subroutine

            subroutine setup_phase_ramp_shifts
            use global_variables, only: nopiy, nopix, ifactory, ifactorx
            implicit none

            integer(4) :: i
            real(fp_kind) :: r_coord

            if(allocated(shift_arrayy)) deallocate(shift_arrayy)
            if(allocated(shift_arrayx)) deallocate(shift_arrayx)
            allocate(shift_arrayy(nopiy,ifactory))
            allocate(shift_arrayx(nopix,ifactorx))

            do i = 1,ifactory
                r_coord=float(i)/float(ifactory)
                call make_shift_oned(shift_arrayy(:,i),nopiy,r_coord)
            enddo

            do i = 1,ifactorx
                r_coord=float(i)/float(ifactorx)
                call make_shift_oned(shift_arrayx(:,i),nopix,r_coord)
            enddo

        end subroutine
end module
