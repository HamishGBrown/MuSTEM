module FFTW3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
!      INTEGER FFTW_FORWARD
!      PARAMETER (FFTW_FORWARD=-1)
!      INTEGER FFTW_BACKWARD
!      PARAMETER (FFTW_BACKWARD=+1)
!      INTEGER FFTW_MEASURE
!      PARAMETER (FFTW_MEASURE=0)
!      INTEGER FFTW_EXHAUSTIVE
!      PARAMETER (FFTW_EXHAUSTIVE=8)
!      INTEGER FFTW_PATIENT
!      PARAMETER (FFTW_PATIENT=32)
!      INTEGER FFTW_ESTIMATE
!      PARAMETER (FFTW_ESTIMATE=64)

    interface fft_plan
        module procedure         sfft2_plan,        sfft1_plan,        dfft1_plan,        dfft2_plan
        module procedure inplace_sfft2_plan,inplace_sfft1_plan,inplace_dfft1_plan,inplace_dfft2_plan
    end interface

    interface fft
        module procedure dfft2,sfft2,dfft1,sfft1
    end interface

    interface inplace_fft
        module procedure inplace_dfft2,inplace_dfft1,inplace_sfft1,inplace_sfft2
    end interface

    interface ifft
        module procedure difft2,sifft2,difft1,sifft1
    end interface

    interface inplace_ifft
        module procedure inplace_difft2,inplace_difft1,inplace_sifft1,inplace_sifft2
    end interface
    contains


!Fourier transform planning

    function planning(plan_mode)
        integer*4,intent(in)::plan_mode
        integer(C_INT)::planning

        if(plan_mode==0)then
            planning = FFTW_MEASURE
        elseif(plan_mode==1) then
            planning = FFTW_PATIENT
        elseif(plan_mode==2) then
            planning = FFTW_EXHAUSTIVE
        else
            planning = FFTW_MEASURE
        endif

    end function

    function sfft2_plan(nopiy,nopix,array_in,array_out,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy,nopix
        complex(C_FLOAT_COMPLEX),dimension(nopiy,nopix)::array_in,array_out
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftwf_plan_dft_2d (nopix,nopiy,array_in,array_out,direction,plan_mode_ )

    end function

    function dfft2_plan(nopiy,nopix,array_in,array_out,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy,nopix
        complex(C_DOUBLE_COMPLEX),dimension(nopiy,nopix)::array_in,array_out
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftw_plan_dft_2d (nopix,nopiy,array_in,array_out,direction,plan_mode_ )

    end function

    function sfft1_plan(nopiy,array_in,array_out,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy
        complex(C_FLOAT_COMPLEX),dimension(nopiy)::array_in,array_out
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftwf_plan_dft_1d (nopiy,array_in,array_out,direction,plan_mode_ )

    end function

    function dfft1_plan(nopiy,array_in,array_out,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy
        complex(C_DOUBLE_COMPLEX),dimension(nopiy)::array_in,array_out
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftw_plan_dft_1d (nopiy,array_in,array_out,direction,plan_mode_ )

    end function

    !Inplace planning

    function inplace_sfft2_plan(nopiy,nopix,array,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy,nopix
        complex(C_FLOAT_COMPLEX),dimension(nopiy,nopix)::array
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftwf_plan_dft_2d (nopix,nopiy,array,array,direction,plan_mode_ )

    end function

    function inplace_dfft2_plan(nopiy,nopix,array,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy,nopix
        complex(C_DOUBLE_COMPLEX),dimension(nopiy,nopix)::array
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftw_plan_dft_2d (nopix,nopiy,array,array,direction,plan_mode_ )

    end function

    function inplace_sfft1_plan(nopiy,array,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy
        complex(C_FLOAT_COMPLEX),dimension(nopiy)::array
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftwf_plan_dft_1d (nopiy,array,array,direction,plan_mode_ )

    end function

    function inplace_dfft1_plan(nopiy,array,direction,plan_mode) result(plan)
        integer*4,intent(in)::nopiy
        complex(C_DOUBLE_COMPLEX),dimension(nopiy)::array
        type(C_PTR) :: plan
        integer(C_INT) :: plan_mode_,direction
        integer*4,intent(in)::plan_mode
        optional:: plan_mode

        if(present(plan_mode)) then
            plan_mode_ = planning(plan_mode)
        else
            plan_mode_ = planning(0)
        endif

        plan = fftw_plan_dft_1d (nopiy,array,array,direction,plan_mode_ )

    end function

!Forward Fourier transforms

    function dfft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical::norm
    complex(C_DOUBLE_COMPLEX)  array_in(nopiy,nopix)
    complex(C_DOUBLE_COMPLEX) dfft2(nopiy,nopix)
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,dfft2)
    else
        plan_ = fftw_plan_dft_2d (nopix,nopiy,array_in,dfft2,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,dfft2)
        call fftw_destroy_plan ( plan_ )
    endif
    if(present(norm)) then
        if(norm) dfft2=dfft2/sqrt(float(nopiy*nopix))
    end if
    end function


    function sfft2(nopiy,nopix,array_in,plan,norm)
    !Single precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX)  array_in(nopiy,nopix)
    complex(C_FLOAT_COMPLEX) sfft2(nopiy,nopix)
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,sfft2)
    else
        plan_ = fftwf_plan_dft_2d (nopix,nopiy,array_in,sfft2,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,sfft2)
        call fftwf_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) sfft2=sfft2/sqrt(float(nopiy*nopix))
    end if
    end function

    function dfft1(nopiy,array_in,plan,norm)
    !Double precision 1D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX)  array_in(nopiy)
    complex(C_DOUBLE_COMPLEX) dfft1(nopiy)
    type(C_PTR) :: plan,plan_
    optional :: plan,norm

    if(present(plan))then
        call fftw_execute_dft(plan,array_in,dfft1)
    else
        plan_ = fftw_plan_dft_1d (nopiy,array_in,dfft1,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,dfft1)
        call fftw_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) dfft1 = dfft1/sqrt(float(nopiy))
    end if

    end function

    function sfft1(nopiy,array_in,plan,norm)
    !Single precision 1D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX)  array_in(nopiy)
    complex(C_FLOAT_COMPLEX) sfft1(nopiy)
    type(C_PTR) :: plan,plan_
    optional :: plan,norm

    if(present(plan))then
        call fftwf_execute_dft(plan,array_in,sfft1)
    else
        plan_ = fftwf_plan_dft_1d (nopiy,array_in,sfft1,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,sfft1)
        call fftwf_destroy_plan (plan_)
    endif

    if(present(norm)) then
        if(norm) sfft1 = sfft1/sqrt(float(nopiy))
    end if

    end function

!Inverse Fourier transforms

    function difft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D inverse Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX)  array_in(nopiy,nopix)
    complex(C_DOUBLE_COMPLEX) difft2(nopiy,nopix)
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,difft2)
    else
        plan_ = fftw_plan_dft_2d (nopix,nopiy,array_in,difft2,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,difft2)
        call fftw_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) difft2 = difft2/sqrt(float(nopiy*nopix))
    end if

    end function

    function sifft2(nopiy,nopix,array_in,plan,norm)
    !Single precision 2D inverse Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX)  array_in(nopiy,nopix)
    complex(C_FLOAT_COMPLEX) sifft2(nopiy,nopix)
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,sifft2)
    else
        plan_ = fftwf_plan_dft_2d (nopix,nopiy,array_in,sifft2,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,sifft2)
        call fftwf_destroy_plan ( plan_ )
    endif
    if(present(norm)) then
        if(norm) sifft2 = sifft2/sqrt(float(nopiy*nopix))
    end if
    end function

    function difft1(nopiy,array_in,plan,norm)
    !Double precision 1D inverse Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX)  array_in(nopiy)
    complex(C_DOUBLE_COMPLEX) difft1(nopiy)
    type(C_PTR) :: plan,plan_
    optional :: plan,norm

    if(present(plan))then
        call fftw_execute_dft(plan,array_in,difft1)
    else
        plan_ = fftw_plan_dft_1d (nopiy,array_in,difft1,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,difft1)
        call fftw_destroy_plan ( plan_ )
    endif

    if(present(norm))then
        if(norm) difft1=difft1/sqrt(float(nopiy))
    end if

    end function

    function sifft1(nopiy,array_in,plan,norm)
    !Single precision 1D inverse Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX)  array_in(nopiy)
    complex(C_FLOAT_COMPLEX) sifft1(nopiy)
    type(C_PTR) :: plan,plan_
    optional :: plan,norm

    if(present(plan))then
        call fftwf_execute_dft(plan,array_in,sifft1)
    else
        plan_ = fftwf_plan_dft_1d (nopiy,array_in,sifft1,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,sifft1)
        call fftwf_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) sifft1=sifft1/sqrt(float(nopiy))
    end if

    end function

    subroutine inplace_dfft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX),dimension(nopiy,nopix)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,array_in)
    else

    if(present(norm)) then
        if(norm) array_in=array_in/sqrt(float(nopiy*nopix))
    end if
        plan_ = fftw_plan_dft_2d (nopix,nopiy,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,array_in)
        call fftw_destroy_plan ( plan_ )
    endif

    end subroutine

    subroutine inplace_sfft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX),dimension(nopiy,nopix)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftwf_plan_dft_2d (nopix,nopiy,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,array_in)
        call fftwf_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) array_in = array_in/sqrt(float(nopiy*nopix))
    end if

    end subroutine

    subroutine inplace_dfft1(nopiy,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX),dimension(nopiy)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftw_plan_dft_1d (nopiy,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,array_in)
        call fftw_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) array_in = array_in /sqrt(float(nopiy))
    end if

    end subroutine

    subroutine inplace_sfft1(nopiy,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX),dimension(nopiy)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftwf_plan_dft_1d (nopiy,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,array_in)
        call fftwf_destroy_plan ( plan_ )
    endif

    if(present(norm)) then
        if(norm) array_in=array_in/sqrt(float(nopiy))
    end if
    end subroutine

    subroutine inplace_difft1(nopiy,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX),dimension(nopiy)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftw_plan_dft_1d (nopiy,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,array_in)
        call fftw_destroy_plan ( plan_ )
    endif
    if(present(norm)) then
        if(norm) array_in = array_in/sqrt(float(nopiy))
    endif
    end subroutine

    subroutine inplace_difft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_DOUBLE_COMPLEX),dimension(nopiy,nopix)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftw_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftw_plan_dft_2d (nopix,nopiy,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftw_execute_dft(plan_,array_in,array_in)
        call fftw_destroy_plan ( plan_ )
    endif
    if(present(norm)) then
        if(norm) array_in = array_in/sqrt(float(nopiy*nopix))
    endif
    end subroutine

    subroutine inplace_sifft2(nopiy,nopix,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy,nopix
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX),dimension(nopiy,nopix)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,array_in)
    else
        plan_ = fftwf_plan_dft_2d (nopix,nopiy,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,array_in)
        call fftwf_destroy_plan ( plan_ )
    endif
    if(present(norm)) then
        if(norm) array_in = array_in/sqrt(float(nopiy*nopix))
    endif
    end subroutine

    subroutine inplace_sifft1(nopiy,array_in,plan,norm)
    !Double precision 2D forward Fourier transform
    integer*4,intent(in)::nopiy
    logical,intent(in)::norm
    complex(C_FLOAT_COMPLEX),dimension(nopiy)::array_in
    type(C_PTR) :: plan,plan_
    optional::plan,norm

    if(present(plan)) then
        call fftwf_execute_dft(plan,array_in,array_in)
    else
    if(present(norm)) then
        if(norm) array_in = array_in/sqrt(float(nopiy))
    end if
        plan_ = fftwf_plan_dft_1d (nopiy,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
        call fftwf_execute_dft(plan_,array_in,array_in)
        call fftwf_destroy_plan ( plan_ )
    endif
    end subroutine

end module
