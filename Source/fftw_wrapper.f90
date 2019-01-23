module FFTW3
#ifdef GPUd
use cufft
implicit none
  save

  interface fft1
	  module procedure dfft1d
	  module procedure sfft1d
	end interface fft1

  interface ifft1
	  module procedure dfft1b
	  module procedure sfft1b
	end interface ifft1

      interface fft2
	  module procedure dfft2d
	  module procedure sfft2d
	end interface fft2

      interface ifft2
	  module procedure dfft2b
	  module procedure sfft2b
	end interface ifft2

      interface fft3
	  module procedure dfft3d
	  module procedure sfft3d
	end interface fft3

      interface ifft3
	  module procedure dfft3b
	  module procedure sfft3b
	end interface ifft3


integer*4::a
      contains
  function dfft1d_fun(nopiy,arrayin,norm,plan) result(arrayout)
    integer*4,intent(in) :: nopiy
    complex(16),intent(in) :: arrayin(nopiy)
    logical,intent(in):norm
    integer,intent(in)::plan

    complex(16) :: arrayout(nopiy)
    complex(16),device,dimension(nopiy):: arrayin_d,arrayout_d
    integer :: plan_

    !Plan FFT if no plan passed to function
    if(present(plan)) then
      plan_ = plan
    else
      call cufftplan(plan,nopiy,CUFFT_Z2Z)
    endif

    !Copy data to device
    arrayin_d = arrayin

    !Execute FFTs
    call cufftExec(plan,arrayin_d,arrayout_d,CUFFT_FORWARD)

    !Destroy plans if created
    if(.not. present(plan)) a= cufftDestroy(plan_)

    arrayout= arrayout_d
    if(norm) arrayout = array_out/(sqrt(float(nopiy)))

  end function
!	forward 1D transform
	subroutine dfft1d(nopiyb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb
	integer(4) :: nopiy,nopix,a
	complex(8),dimension(nopiyb) :: array_in
	complex(8),dimension(nopiyb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb) :: array_in_d
	complex(8),device,dimension(nopiyb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d

    array_out=array_out/(dsqrt(dfloat(nopiyb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 1D transform
	subroutine dfft1b(nopiyb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb
	integer(4) :: nopiy,nopix
	complex(8),dimension(nopiyb) :: array_in
	complex(8),dimension(nopiyb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb) :: array_in_d
	complex(8),device,dimension(nopiyb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d

    array_out=array_out/(dsqrt(dfloat(nopiyb)))

	return
	end subroutine






!----------------------------------------------------------------------------------------
!	forward 2D transform
	subroutine dfft2d(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb,nopixb
	integer(4) :: nopiy,nopix
	complex(8),dimension(nopiyb,nopixb) :: array_in
	complex(8),dimension(nopiyb,nopixb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb,nopixb) :: array_in_d
	complex(8),device,dimension(nopiyb,nopixb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopixb,nopiyb,CUFFT_Z2Z)
    !call cufftPlan(plan,nopiyb,nopixb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 2D transform
	subroutine dfft2b(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb,nopixb
	integer(4) :: nopiy,nopix !dummy
	complex(8),dimension(nopiyb,nopixb) :: array_in
	complex(8),dimension(nopiyb,nopixb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb,nopixb) :: array_in_d
	complex(8),device,dimension(nopiyb,nopixb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopixb,nopiyb,CUFFT_Z2Z)
    !call cufftPlan(plan,nopiyb,nopixb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	forward 3D transform
	subroutine dfft3d(nopiyb,nopixb,nopizb,array_in,nopiy,nopix,array_out,nopiya,nopixa)

	implicit none
	integer	plan
	integer(4) :: nopiya,nopixa,nopiza
	integer(4) :: nopiyb,nopixb,nopizb
	integer(4) :: nopiy,nopix,nopiz
	complex(8),dimension(nopiyb,nopixb,nopizb) :: array_in
	complex(8),dimension(nopiyb,nopixb,nopizb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb,nopixb,nopizb) :: array_in_d
	complex(8),device,dimension(nopiyb,nopixb,nopizb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,nopixb,nopizb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb*nopizb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 3D transform
	subroutine dfft3b(nopiyb,nopixb,nopizb,array_in,nopiy,nopix,array_out,nopiya,nopixa)

	implicit none
	integer	plan
	integer(4) :: nopiya,nopixa,nopiza
	integer(4) :: nopiyb,nopixb,nopizb
	integer(4) :: nopiy,nopix,nopiz
	complex(8),dimension(nopiyb,nopixb,nopizb) :: array_in
	complex(8),dimension(nopiyb,nopixb,nopizb) :: array_out

	!device arrays
	complex(8),device,dimension(nopiyb,nopixb,nopizb) :: array_in_d
	complex(8),device,dimension(nopiyb,nopixb,nopizb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,nopixb,nopizb,CUFFT_Z2Z)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb*nopizb)))
	return
	end subroutine


      !----------------------------------------------------------------------------------------------------------------------------------



      !	forward 1D transform
	subroutine sfft1d(nopiyb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb
	integer(4) :: nopiy,nopix
	complex(4),dimension(nopiyb) :: array_in
	complex(4),dimension(nopiyb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb) :: array_in_d
	complex(4),device,dimension(nopiyb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 1D transform
	subroutine sfft1b(nopiyb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb
	integer(4) :: nopiy,nopix
	complex(4),dimension(nopiyb) :: array_in
	complex(4),dimension(nopiyb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb) :: array_in_d
	complex(4),device,dimension(nopiyb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb)))
	return
	end subroutine






!----------------------------------------------------------------------------------------
!	forward 2D transform
	subroutine sfft2d(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb,nopixb
	integer(4) :: nopiy,nopix
	complex(4),dimension(nopiyb,nopixb) :: array_in
	complex(4),dimension(nopiyb,nopixb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb,nopixb) :: array_in_d
	complex(4),device,dimension(nopiyb,nopixb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
    call cufftPlan(plan,nopixb,nopiyb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb*nopixb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 2D transform
	subroutine sfft2b(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	integer	plan
	integer(4) :: nopiyb,nopixb
	integer(4) :: nopiy,nopix
	complex(4),dimension(nopiyb,nopixb) :: array_in
	complex(4),dimension(nopiyb,nopixb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb,nopixb) :: array_in_d
	complex(4),device,dimension(nopiyb,nopixb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopixb,nopiyb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb*nopixb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	forward 3D transform
	subroutine sfft3d(nopiyb,nopixb,nopizb,array_in,nopiy,nopix,array_out,nopiya,nopixa)

	implicit none
	integer	plan
	integer(4) :: nopiya,nopixa,nopiza
	integer(4) :: nopiyb,nopixb,nopizb
	integer(4) :: nopiy,nopix,nopiz
	complex(4),dimension(nopiyb,nopixb,nopizb) :: array_in
	complex(4),dimension(nopiyb,nopixb,nopizb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb,nopixb,nopizb) :: array_in_d
	complex(4),device,dimension(nopiyb,nopixb,nopizb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,nopixb,nopizb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_FORWARD)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb*nopixb*nopizb)))
	return
	end subroutine

!----------------------------------------------------------------------------------------
!	inverse 3D transform
	subroutine sfft3b(nopiyb,nopixb,nopizb,array_in,nopiy,nopix,array_out,nopiya,nopixa)

	implicit none
	integer	plan
	integer(4) :: nopiya,nopixa,nopiza
	integer(4) :: nopiyb,nopixb,nopizb
	integer(4) :: nopiy,nopix,nopiz
	complex(4),dimension(nopiyb,nopixb,nopizb) :: array_in
	complex(4),dimension(nopiyb,nopixb,nopizb) :: array_out

	!device arrays
	complex(4),device,dimension(nopiyb,nopixb,nopizb) :: array_in_d
	complex(4),device,dimension(nopiyb,nopixb,nopizb) :: array_out_d

	!copy data to device
	array_in_d=array_in

	! Initialize the plan
	call cufftPlan(plan,nopiyb,nopixb,nopizb,CUFFT_C2C)

	! Execute FFTs
	call cufftExec(plan,array_in_d,array_out_d,CUFFT_INVERSE)

	! Destroy plans
	a = cufftDestroy(plan)

	! Copy results back to host
	array_out=array_out_d
    array_out=array_out/(sqrt(float(nopiyb*nopixb*nopizb)))
	return
	end subroutine
#else

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

        plan = fftwf_plan_dft_2d (nopiy,nopix,array_in,array_out,direction,plan_mode_ )

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

        plan = fftw_plan_dft_2d (nopiy,nopix,array_in,array_out,direction,plan_mode_ )

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

        plan = fftwf_plan_dft_2d (nopiy,nopix,array,array,direction,plan_mode_ )

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

        plan = fftw_plan_dft_2d (nopiy,nopix,array,array,direction,plan_mode_ )

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
        plan_ = fftw_plan_dft_2d (nopiy,nopix,array_in,dfft2,FFTW_FORWARD,FFTW_ESTIMATE )
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
        plan_ = fftwf_plan_dft_2d (nopiy,nopix,array_in,sfft2,FFTW_FORWARD,FFTW_ESTIMATE )
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
        plan_ = fftw_plan_dft_2d (nopiy,nopix,array_in,difft2,FFTW_BACKWARD,FFTW_ESTIMATE )
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
        plan_ = fftwf_plan_dft_2d (nopiy,nopix,array_in,sifft2,FFTW_BACKWARD,FFTW_ESTIMATE )
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
        plan_ = fftw_plan_dft_2d (nopiy,nopix,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
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
        plan_ = fftwf_plan_dft_2d (nopiy,nopix,array_in,array_in,FFTW_FORWARD,FFTW_ESTIMATE )
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
        plan_ = fftw_plan_dft_2d (nopiy,nopix,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
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
        plan_ = fftwf_plan_dft_2d (nopiy,nopix,array_in,array_in,FFTW_BACKWARD,FFTW_ESTIMATE )
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
#endif
end module
