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
	!use output
    
    implicit none

    interface
        subroutine make_site_factor_generic(site_factor, tau)
            use m_precision, only: fp_kind
            complex(fp_kind) :: site_factor(:, :)
            real(fp_kind) :: tau(:,:)
        end subroutine make_site_factor_generic
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
    
    
    
    subroutine make_g_vec_array(g_vec_array)

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

    end subroutine make_g_vec_array
    
    
    
    subroutine precalculate_scattering_factors
    
	    use m_precision, only: fp_kind
        use global_variables, only: ifactory, ifactorx, nopiy, nopix, nt, a0, atf, atomf, ak, fz, fz_dwf, sinc, inverse_sinc, tp
        use m_elsa, only: elsa_ext
        use m_absorption, only: complex_absorption, setup_absorptive_array, max_int, delta_kstep, tdsbr, fz_abs,calculate_absorption_mu
		use m_numerical_tools, only: cubspl,ppvalu

	    implicit none
    
        integer(4) :: shiftx, shifty, m1, m2, i, j, k
        real(fp_kind) :: s, s2, g2, sky, skx, ax, ay
        real(fp_kind) :: el_scat
        real(fp_kind) :: xkstep, temp
        real(fp_kind),allocatable :: tdsbrcoeff(:,:), xdata(:),tdsbrc(:,:,:) 
        real(fp_kind) :: factor, eps
        real(fp_kind) :: xlen, ylen

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
            !if(allocated(tdsbrcoeff)) deallocate(tdsbrcoeff)
            if(allocated(xdata)) deallocate(xdata)
            if(allocated(fz_abs)) deallocate(fz_abs)
        
            allocate(xdata(max_int))
           !allocate(tdsbrcoeff(max_int,nt))
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
        xlen = a0(1)*float(ifactorx)
        ylen = a0(2)*float(ifactory)
        factor = 1.0_fp_kind
        eps = tiny(0.0_fp_kind)

        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1
        
        !!$OMP PARALLEL PRIVATE(i, m2, sky, j, m2, skx, g2, s2, k, el_scat, temp)
        !!$OMP DO
	    do i = 1, nopiy
			
	        m2 = mod( i+shifty, nopiy) - shifty -1
            sky = float(m2)/ylen 
	        do j = 1, nopix
	            m1 = mod( j+shiftx, nopix) - shiftx -1
                skx = float(m1)/xlen
                g2 =  sky**2.0_fp_kind+skx**2.0_fp_kind
                s2 = g2 / 4.0_fp_kind
                
                do k = 1, nt
                    ! Multiply by fractional occupancy
                    el_scat = elsa_ext(nt,k,atomf,s2) * atf(2,k)    
                    
                
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
            enddo
        enddo
        !!$OMP END DO
	    !!$OMP END PARALLEL

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
	    do i = 1, nopiy
	        do j = 1, nopix
                site_factor(i,j) = sum(exp(cmplx(0.0_fp_kind, -tp*matmul(g_vec_array(:,i,j), tau), fp_kind)))
	   	    enddo
        enddo
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
        
        integer :: i, j
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
    

        
end module
