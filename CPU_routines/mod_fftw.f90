
module cufft_wrapper
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



    contains
          
 
 
    subroutine setup_threading()
    implicit none
      
    integer*4   num_threads,junk
    integer*4   OMP_GET_MAX_THREADS
    external    OMP_GET_MAX_THREADS
    
    !Set the number of threads allowed
    num_threads=OMP_GET_MAX_THREADS()
    num_threads=num_threads/2
    CALL OMP_SET_NUM_THREADS(num_threads)
    call dfftw_init_threads(junk)
    call dfftw_plan_with_nthreads(num_threads)
    write(*,*) '|----------------------------------------------------------------------------|'
    write(*,*) '        The number of threads available is: ',num_threads*2
    write(*,*) '        The number of threads being used is: ',num_threads
    write(*,*) '|----------------------------------------------------------------------------|'
    return
    end subroutine
    
    
    subroutine dfft1b(nopixb,array_in,array_out)

	implicit none
	include 'fftw3.f'
	integer*8	plan
	integer*4	nopixb
	complex*16  array_in(nopixb)
	complex*16  array_out(nopixb)

    call dfftw_plan_dft (plan,nopixb,array_in,array_out,FFTW_BACKWARD,FFTW_ESTIMATE )
	call dfftw_execute ( plan)
	call dfftw_destroy_plan ( plan )
	array_out=array_out/(dsqrt(dfloat(nopixb)))

	return
	end subroutine
	
	subroutine dfft1d(nopixb,array_in,array_out)
	implicit none
	include 'fftw3.f'
	integer*8	plan
	integer*4	nopixb
	complex*16  array_in(nopixb)
	complex*16  array_out(nopixb)

    call dfftw_plan_dft (plan,nopixb,array_in,array_out,FFTW_FORWARD,FFTW_ESTIMATE )
	call dfftw_execute ( plan)
	call dfftw_destroy_plan ( plan )
    array_out=array_out/(dsqrt(dfloat(nopixb)))
    
	return
	end subroutine
	
	subroutine dfft2b(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	include 'fftw3.f'
	integer*8	plan
	integer*4	nopiyb,nopixb
	integer*4   nopiy,nopix !dummy variables
	complex*16  array_in(nopiyb,nopixb)
	complex*16  array_out(nopiyb,nopixb)

    call dfftw_plan_dft_2d (plan,nopiyb,nopixb,array_in,array_out,FFTW_BACKWARD,FFTW_ESTIMATE )
	call dfftw_execute ( plan)
	call dfftw_destroy_plan ( plan )
	array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb)))

	return
	end subroutine

	subroutine dfft2d(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	include 'fftw3.f'
	integer*8	plan
	integer*4	nopiyb,nopixb
	integer*4   nopiy,nopix !dummy variables
	complex*16  array_in(nopiyb,nopixb)
	complex*16  array_out(nopiyb,nopixb)


    call dfftw_plan_dft_2d (plan,nopiyb,nopixb,array_in,array_out,FFTW_FORWARD,FFTW_ESTIMATE )
	call dfftw_execute ( plan)
	call dfftw_destroy_plan ( plan )
    array_out=array_out/(dsqrt(dfloat(nopiyb*nopixb)))
    	
	return
	end subroutine
	
	
	
	subroutine sfft1b(nopixb,array_in,array_out)

	implicit none
	include 'fftw3.f'
	integer*8  plan
	integer*4  nopixb
	complex*8  array_in(nopixb)
	complex*8  array_out(nopixb)

    call sfftw_plan_dft (plan,nopixb,array_in,array_out,FFTW_BACKWARD,FFTW_ESTIMATE )
	call sfftw_execute ( plan)
	call sfftw_destroy_plan ( plan )
	array_out=array_out/(sqrt(float(nopixb)))

	return
	end subroutine
	
	subroutine sfft1d(nopixb,array_in,array_out)
	implicit none
	include 'fftw3.f'
	integer*8	plan
	integer*4	nopixb
	complex*8  array_in(nopixb)
	complex*8  array_out(nopixb)

    call sfftw_plan_dft (plan,nopixb,array_in,array_out,FFTW_FORWARD,FFTW_ESTIMATE )
	call sfftw_execute ( plan)
	call sfftw_destroy_plan ( plan )
    array_out=array_out/(sqrt(float(nopixb)))
    
	return
	end subroutine
	
	subroutine sfft2b(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	include 'fftw3.f'
	integer*8 plan
	integer*4 nopiyb,nopixb
	integer*4 nopiy,nopix !dummy variables
	complex*8 array_in(nopiyb,nopixb)
	complex*8 array_out(nopiyb,nopixb)

    call sfftw_plan_dft_2d (plan,nopiyb,nopixb,array_in,array_out,FFTW_BACKWARD,FFTW_ESTIMATE )
	call sfftw_execute ( plan)
	call sfftw_destroy_plan ( plan )
	array_out=array_out/(sqrt(float(nopiyb*nopixb)))

	return
	end subroutine

	subroutine sfft2d(nopiyb,nopixb,array_in,nopiy,array_out,nopix)

	implicit none
	include 'fftw3.f'
	integer*8 plan
	integer*4 nopiyb,nopixb
	integer*4 nopiy,nopix !dummy variables
	complex*8 array_in(nopiyb,nopixb)
	complex*8 array_out(nopiyb,nopixb)


    call sfftw_plan_dft_2d (plan,nopiyb,nopixb,array_in,array_out,FFTW_FORWARD,FFTW_ESTIMATE )
	call sfftw_execute ( plan)
	call sfftw_destroy_plan ( plan )
    array_out=array_out/(sqrt(float(nopiyb*nopixb)))
    	
	return
	end subroutine
    end module cufft_wrapper


