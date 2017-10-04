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

module cuda_ms
    
    use global_variables
    use cuda_array_library
    use cudafor
    
    implicit none

    interface get_sum
        module procedure get_sum_complex
        module procedure get_sum_real
    end interface

    interface cuda_stem_detector
        module procedure cuda_stem_detector_wavefunction
        module procedure cuda_stem_detector_cbed
    end interface
    
    contains


    
    attributes(global) subroutine cuda_phase_shift(shift_array, ifactory, ifactorx, ig1, ig2, n, m, coord)
    
        implicit none
    
        integer(4), value :: n,m,ix,iy
        integer(4), value :: ifactory,ifactorx
        integer(4) :: m1,m2
        integer(4),device,dimension(3) :: ig1,ig2,kx,ky
        real(fp_kind),device,dimension(3) :: kr
        real(fp_kind),device,dimension(3) :: coord
        real(fp_kind), value :: tp

        complex(fp_kind),dimension(n,m) :: shift_array

        tp = atan(1.0_fp_kind)*8.0_fp_kind
        
        !calculate the number of threads and blocks
        ix = (blockIdx%x-1)*blockDim%x + threadIdx%x
        iy = (blockIdx%y-1)*blockDim%y + threadIdx%y

        if(ix <= n) then
            if (iy <= m) then
                m1 = mod(ix+(n-1)/2-1,n)-(n-1)/2
                ky = m1*ig1
                m2 = mod(iy+(m-1)/2-1,m)-(m-1)/2           
                kx = m2*ig2
                kr = kx/float(ifactorx)+ky/float(ifactory)
                shift_array(ix,iy) = exp(cmplx(0.0_fp_kind,-tp*dot_product( kr, coord) ))
            endif
        endif   

    end subroutine cuda_phase_shift

    
    
    attributes(host) function cuda_stem_detector_wavefunction(psi_d, mask_d)
    
        use cuda_array_library, only: blocks, threads
    
        implicit none
    
        real(fp_kind) :: cuda_stem_detector_wavefunction
        complex(fp_kind), device, dimension(nopiy,nopix) :: psi_d
        real(fp_kind), device, dimension(nopiy,nopix) :: mask_d, cbed_d

        call cuda_mod<<<blocks,threads>>>(psi_d,cbed_d,normalisation,nopiy,nopix)
        call cuda_multiplication<<<blocks,threads>>>(cbed_d,mask_d, cbed_d ,1.0_fp_kind,nopiy,nopix)

        cuda_stem_detector_wavefunction = get_sum(cbed_d)
        
    end function cuda_stem_detector_wavefunction
    
    
    
    attributes(host) function cuda_stem_detector_cbed(cbed_d, mask_d)
    
        use cuda_array_library, only: blocks, threads
        
        implicit none
    
        real(fp_kind) :: cuda_stem_detector_cbed
        real(fp_kind), device, dimension(nopiy,nopix) :: mask_d, image_d, cbed_d

        call cuda_multiplication<<<blocks,threads>>>(cbed_d,mask_d, image_d ,1.0_fp_kind,nopiy,nopix)
    
        cuda_stem_detector_cbed = get_sum(image_d)
    
    end function cuda_stem_detector_cbed
    
    

    attributes(host) function get_sum_complex(psi_d)
    
        implicit none
    
        integer(4) l,m
        real(fp_kind) :: get_sum_complex
        real(fp_kind),device :: get_sum_complex_d
        complex(fp_kind), device, dimension(nopiy,nopix) :: psi_d
    
        get_sum_complex_d = 0.0_fp_kind
    
        !$cuf kernel do (2) <<<*,*>>>
        do m = 1, nopix
            do l = 1, nopiy
                get_sum_complex_d = get_sum_complex_d + abs(psi_d(l,m))**2
            enddo
        enddo
        
        get_sum_complex = get_sum_complex_d
    
    end function get_sum_complex
    
    
    
    attributes(host) function get_sum_real(psi_d)
    
        implicit none
        
        integer(4) l, m
        real(fp_kind) :: get_sum_real
        real(fp_kind), device, dimension(nopiy,nopix) :: psi_d
        real(fp_kind),device :: get_sum_real_d
    
        get_sum_real_d = 0.0_fp_kind
    
        !$cuf kernel do (2) <<<*,*>>>
        do m = 1, nopix
          do l = 1, nopiy
              get_sum_real_d = get_sum_real_d + psi_d(l,m)
          enddo
        enddo
        
        get_sum_real = get_sum_real_d
    
    end function get_sum_real

    

end module
