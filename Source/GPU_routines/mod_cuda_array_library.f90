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

module cuda_array_library
    
    use m_precision
    use cudafor, only: dim3
    
    type(dim3),parameter :: threads = dim3(32, 32, 1)
    type(dim3) :: blocks
    
    interface cuda_scale
        module procedure cuda_complex_scale
        module procedure cuda_real_scale
    end interface cuda_scale
    
    interface cuda_multiplication
	    module procedure cuda_complex_multiplication
	    module procedure cuda_real_multiplication
		module procedure cuda_complex_factorized_multiplication
		module procedure cuda_cshifted_complex_multiplication
	end interface cuda_multiplication

    interface cuda_addition
	    module procedure cuda_complex_addition
	    module procedure cuda_real_addition
	end interface cuda_addition

    interface cuda_subtraction
	    module procedure cuda_complex_subtraction
	    module procedure cuda_real_subtraction
	end interface cuda_subtraction

    interface cuda_cshift
        module procedure cuda_cshift_complex
        module procedure cuda_cshift_real
    end interface cuda_cshift

    
    
    contains

    
    
    subroutine set_blocks(nopiy, nopix)

        implicit none
        
        integer :: nopiy, nopix
        
        blocks = dim3(ceiling(float(nopiy)/threads%x), ceiling(float(nopix)/threads%y), 1)
        
    end subroutine
    
    
    
    attributes(global) subroutine cuda_copy_real_to_cmplx(a, b, n, m)
    
        implicit none
     
        integer(4),value :: n, m, ix, iy
        real(fp_kind),dimension(n,m) :: a
        complex(fp_kind),dimension(n,m) :: b
        
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                b(iy,ix) = a(iy,ix)
            endif
        endif   
    
    end subroutine cuda_copy_real_to_cmplx

    
    
    attributes(global) subroutine cuda_complex_scale(dTemp, scale, fTemp, n, m)
    
        implicit none
     
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale
        complex(fp_kind),dimension(n,m) :: dTemp, fTemp
        
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                fTemp(iy,ix) = dTemp(iy,ix) * scale      
            endif
        endif   
    
    end subroutine cuda_complex_scale

    
    
    attributes(global) subroutine cuda_real_scale(dTemp, scale, fTemp, n, m)
    
        implicit none
    
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale
        real(fp_kind),dimension(n,m) :: dTemp, fTemp
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (iy <= n) then
            if (ix <= m) then
               fTemp(iy,ix) = dTemp(iy,ix) * scale   
            endif
        endif   
    
    end subroutine cuda_real_scale
    
    
    attributes(global) subroutine cuda_complex_factorized_multiplication(dTemp, e1Temp, e2Temp,fTemp, scale, n, m)
    
        implicit none
     
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale
        complex(fp_kind),dimension(n,m) :: dTemp, fTemp
		complex(fp_kind)::e1Temp(n),e2Temp(m)
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                fTemp(iy,ix) = dTemp(iy,ix) * e1Temp(iy)*e2Temp(ix) * scale      
            endif
        endif   
    
    end subroutine cuda_complex_factorized_multiplication
    
    
    attributes(global) subroutine cuda_complex_multiplication(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
     
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale
        complex(fp_kind),dimension(n,m) :: dTemp, eTemp, fTemp
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                fTemp(iy,ix) = dTemp(iy,ix) * eTemp(iy,ix) * scale      
            endif
        endif   
    
    end subroutine cuda_complex_multiplication

    attributes(global) subroutine cuda_cshifted_complex_multiplication(dTemp, eTemp, fTemp, scale, n, m,nshifty,nshiftx)
    
        implicit none
     
        integer(4),value :: n, m, ix, iy,nshifty,nshiftx
        real(fp_kind),value :: scale
        complex(fp_kind),dimension(n,m) :: dTemp, eTemp, fTemp
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (ix <= m) then
             if (iy <= n) then
                fTemp(iy,ix) = dTemp(iy,ix) * eTemp(modulo(iy-1+nshifty,n)+1,modulo(ix-1+nshiftx,m)+1) * scale      
            endif
        endif   
    
    end subroutine cuda_cshifted_complex_multiplication
    
    

    attributes(global) subroutine cuda_real_multiplication(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        integer(4),value :: n,m,ix,iy
        real(fp_kind),value :: scale
        real(fp_kind),dimension(n,m) :: dTemp,eTemp,fTemp

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (iy <= n) then
            if (ix <= m) then
               fTemp(iy,ix) = dTemp(iy,ix) * eTemp(iy,ix) * scale   
            endif
        endif   
    
    end subroutine cuda_real_multiplication

    
    
    attributes(global) subroutine cuda_multiplication_correlation(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: dTemp, eTemp, fTemp
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = scale * conjg(eTemp(iy,ix)) * dTemp(iy,ix)       
            endif
        endif
    
    end subroutine cuda_multiplication_correlation

    
    
    attributes(global) subroutine cuda_complex_addition(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: dTemp, eTemp, fTemp
        real(fp_kind),value :: scale
        integer(4),value :: n, m, ix, iy

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = dTemp(iy,ix) + scale * eTemp(iy,ix)     
            endif
        endif
    
    end subroutine cuda_complex_addition
    
    
    
    attributes(global) subroutine cuda_complex_subtraction(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: dTemp, eTemp, fTemp
        real(fp_kind),value :: scale
        integer(4),value :: n, m, ix, iy

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = dTemp(iy,ix) - scale * eTemp(iy,ix)     
            endif
        endif
    
    end subroutine cuda_complex_subtraction
           
    
    
    attributes(global) subroutine cuda_real_addition(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        real(fp_kind),dimension(n,m),device :: dTemp, eTemp, fTemp
        real(fp_kind),value :: scale
        integer(4),value :: n, m, ix, iy

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = dTemp(iy,ix) + scale * eTemp(iy,ix)    
            endif
        endif
    
    end subroutine cuda_real_addition

    
    
    attributes(global) subroutine cuda_real_subtraction(dTemp, eTemp, fTemp, scale, n, m)
    
        implicit none
    
        real(fp_kind),dimension(n,m),device :: dTemp, eTemp, fTemp
        real(fp_kind),value :: scale
        integer(4),value :: n, m, ix, iy

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = dTemp(iy,ix) - scale * eTemp(iy,ix) 
            endif
        endif
    
    end subroutine cuda_real_subtraction
    
    
    
    attributes(global) subroutine cuda_mod(dTemp, fTemp, scale, n, m)
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: dTemp
        real(fp_kind),dimension(n,m),device :: fTemp
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                fTemp(iy,ix) = scale * conjg(dTemp(iy,ix)) * dTemp(iy,ix)   
            endif
        endif
    
    end subroutine cuda_mod
    
    

    attributes(global) subroutine cuda_cshift_real(in_array, out_array, n, m, yshift, xshift)
    
        implicit none

        integer(4),value :: n, m, xshift, yshift
        integer(4) :: ix, iy
	    real(fp_kind),device,dimension(n,m) :: in_array, out_array

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (ix <= m) then
            if (iy <= n) then
                out_array(iy,ix) = in_array(modulo(iy+yshift-1,n)+1,modulo(ix+xshift-1,m)+1)
            endif
        endif
    
	end subroutine cuda_cshift_real
    
    

    attributes(global) subroutine cuda_cshift_complex(in_array, out_array, n, m, yshift, xshift)
    
        implicit none

        integer(4),value :: n, m, xshift, yshift
        integer(4) :: ix, iy
	    complex(fp_kind),device,dimension(n,m) :: in_array, out_array

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                out_array(iy,ix) = in_array(modulo(iy+yshift-1,n)+1,modulo(ix+xshift-1,m)+1)
            endif
        endif
    
	end subroutine cuda_cshift_complex

    
    
    attributes(global) subroutine cuda_make_shift_array(shift_array, shifty, shiftx, n, m)
    
        implicit none
    
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale
        complex(fp_kind),dimension(n,m) :: shift_array
        complex(fp_kind),dimension(n) :: shifty
        complex(fp_kind),dimension(m) :: shiftx

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if (ix <= m) then
                   shift_array(iy,ix) = shiftx(ix) * shifty(iy)   
             endif
        endif   
        
    end subroutine cuda_make_shift_array
     
    
    
    attributes(global) subroutine cuda_complex_exponentiation(aTemp, bTemp,scale, n, m)
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: aTemp, bTemp
        integer(4),value :: n, m, ix, iy
        real(fp_kind),value :: scale

        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
        
        if (iy <= n) then
            if (ix <= m) then
                aTemp(iy,ix) = exp(cmplx(0.0_fp_kind, 1.0_fp_kind) * scale * bTemp(iy,ix))   
            endif
        endif
    
    end subroutine cuda_complex_exponentiation
    
    
    
    attributes(global) subroutine cuda_take_real(in_d, out_d, n, m)   
    
        implicit none
    
        complex(fp_kind),dimension(n,m),device :: in_d
        real(fp_kind),dimension(n,m),device :: out_d
        integer(4),value :: n, m, ix, iy
    
        ix = (blockIdx%y-1)*blockDim%y + threadIdx%y
        iy = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
        if (iy <= n) then
            if(ix <= m) then
                out_d(iy,ix) = real(in_d(iy,ix))   
            endif
        endif
    
    end subroutine 
       

    
end module cuda_array_library
    