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

module output

    use m_precision, only: fp_kind
     implicit none

    interface array_from_txt_file
        module procedure array_from_txt_real_1d,array_from_txt_real
    end interface

    interface quad_shift
        module procedure quad_shift_complex,quad_shift_real
    end interface

    interface binary_out
        module procedure binary_out_real, binary_out_complex
    end interface

    character(120) :: output_prefix
    logical :: timing

    contains

    integer*4 function nlines(fnam)
    character*(*),intent(in)::fnam
    integer*4::io
    nlines = 0
    OPEN (50, file = fnam)
    DO
      READ(50,*,iostat=io)
      IF (io/=0) EXIT
      nlines = nlines + 1
    END DO
    CLOSE (50)
    end function

    function array_from_txt_real_1d(fnam,nopiy) result(array_from_txt)

    character*(*),intent(in)::fnam
    integer*4,intent(in)::nopiy

    real(fp_kind)::array_from_txt(nopiy)
    integer*4::i

    open(unit = 50,file = fnam,status='old',action='read')
    do i=1,nopiy
        read(50,*) array_from_txt(i)
    enddo

    close(50)
    end function

    function array_from_txt_real(fnam,nopiy,nopix) result(array_from_txt)

    character*(*),intent(in)::fnam
    integer*4,intent(in)::nopiy,nopix

    real(fp_kind)::array_from_txt(nopiy,nopix)
    integer*4::i,j

    open(unit = 50,file = fnam,status='old',action='read')
    do i=1,nopix
        read(50,*) (array_from_txt(j,i),j=1,nopiy)
    enddo

    close(50)
    end function

    subroutine setup_output_prefix

        use m_user_input, only: get_input

        implicit none
        character(120)::dir,fnam

10      write(6,*) 'Enter the prefix for all outputted filenames:'
        call get_input("Output filename", output_prefix)
        write(*,*)

        if (len_trim(output_prefix).eq.0) goto 10
        call split_filepath(trim(adjustl(output_prefix)),dir,fnam)


        if(len(trim(adjustl(dir)))>0) call system('mkdir '//trim(adjustl(dir)))

        output_prefix = trim(adjustl(output_prefix))

    end subroutine

    subroutine split_filepath(filepath,dir,fnam)
              character(len=*), intent(in) :: filepath
              character(len=*), intent(out):: dir,fnam

              integer*4::j


                j = index(filepath,'/',back=.true.)
                j = max(j,index(filepath,'\',back=.true.))

                if(j>0) then
                    dir = trim(adjustl(filepath(:j-1)))
                    fnam = trim(adjustl(filepath(j:)))
                else
                    dir = ''
                    fnam = trim(adjustl(filepath))
                endif
    end subroutine

    subroutine binary_in(nopiy, nopix, array, filename)
        implicit none

        integer(4) :: nopiy, nopix
        real(fp_kind) :: array(nopiy,nopix)
        character(*) :: filename

        integer :: iostat

        open(unit=19, file=filename, form='unformatted', status='old', convert='big_endian', iostat=iostat)

        if (iostat.ne.0) then
            write(*,*) 'Error reading binary file "', trim(filename), '".'
            write(*,*) 'The program will now halt.'
            pause
            stop

        endif

        read(19, iostat=iostat) array

        if (iostat.ne.0) then
            write(*,*) 'Error reading binary file "', trim(filename), '".'
            write(*,*) 'The program will now halt.'
            pause
            stop
        endif

        close(19)

    end subroutine

    subroutine binary_out_real(nopiy, nopix, array, filename,write_to_screen)

        implicit none

        integer(4) :: nopiy, nopix
        real(fp_kind) :: array(nopiy,nopix)
        character(*) :: filename
        logical,intent(in),optional::write_to_screen

        character(512) :: fnam_out
        character(5) :: ydim,xdim
        integer::ierr

        write(ydim,'(i5)') nopiy
        write(xdim,'(i5)') nopix

        fnam_out = trim(adjustl(filename)) // '_' // trim(adjustl(xdim)) // 'x' // trim(adjustl(ydim)) // '.bin'

        if (present(write_to_screen)) then
            if (write_to_screen) write(*,*) trim(fnam_out)
        else
            write(*,*) trim(fnam_out)
        endif

        open(unit=18, file=fnam_out, form='unformatted', status='unknown', convert='big_endian', iostat=ierr)
        write(18) transpose(array)
        close(18)

        if (ierr.ne.0) then
            write(*,*) '      An error occurred while writing to binary file.'
        endif

    end subroutine

    subroutine binary_out_unwrap(nopiy, nopix, array, filename,write_to_screen,to_bandlimit,nopiyout,nopixout)
        ! Unwrap and output a diffraction pattern

        implicit none

        integer(4) :: nopiy, nopix,nopiy_cbedout,nopix_cbedout
        integer(4),intent(in),optional::nopiyout,nopixout
        real(fp_kind) :: array(nopiy,nopix)
        character(*) :: filename
        real(fp_kind) :: array_unwrapped(nopiy,nopix)
        logical,optional,intent(in)::write_to_screen,to_bandlimit
        logical:: write_to_screen_,to_bandlimit_


        if(present(write_to_screen)) then
            write_to_screen_ = write_to_screen
        else
            write_to_screen_ = .true.
        endif

        if(present(to_bandlimit)) then
            to_bandlimit_ = to_bandlimit
        else
            to_bandlimit_ = .true.
        endif

        if (.not. (present(nopiyout).and.present(nopixout))) then
            nopiy_cbedout = nopiy*2/3
            if(to_bandlimit_) nopiy_cbedout = nopiy_cbedout
            nopix_cbedout = nopix*2/3
            if(to_bandlimit_) nopix_cbedout = nopix_cbedout
        else
            nopiy_cbedout = nopiyout
            nopix_cbedout = nopixout
        endif

        array_unwrapped = cshift(cshift(array, -nopix_cbedout/2,dim = 2),-nopiy_cbedout/2,dim = 1)
        call binary_out(nopiy_cbedout,nopix_cbedout,array_unwrapped(1:nopiy_cbedout,1:nopix_cbedout),filename,write_to_screen_)

    end subroutine

    function quad_shift_real(array,nopiy,nopix)
        real(fp_kind),intent(in)::array(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix

        real(fp_kind)::quad_shift_real(nopiy,nopix)
        integer::Npos_y, Npos_x, Nneg_y, Nneg_x
        integer::shifty, shiftx

        Npos_y = int(floor(float(nopiy)/2))
        Npos_x = int(floor(float(nopix)/2))

        Nneg_y = nopiy - Npos_y - 1
        Nneg_x = nopix - Npos_x - 1

        shifty = -Nneg_y
        shiftx = -Nneg_x

        quad_shift_real = cshift(cshift(array, shifty, 1), shiftx, 2)

    end function

    function quad_shift_complex(array,nopiy,nopix)
        complex(fp_kind),intent(in)::array(nopiy,nopix)
        integer*4,intent(in)::nopiy,nopix

        complex(fp_kind)::quad_shift_complex(nopiy,nopix)
        integer::Npos_y, Npos_x, Nneg_y, Nneg_x
        integer::shifty, shiftx

        Npos_y = int(floor(float(nopiy)/2))
        Npos_x = int(floor(float(nopix)/2))

        Nneg_y = nopiy - Npos_y - 1
        Nneg_x = nopix - Npos_x - 1

        shifty = -Nneg_y
        shiftx = -Nneg_x

        quad_shift_complex = cshift(cshift(array, shifty, 1), shiftx, 2)

    end function



    subroutine printout_1d(image, nopix, filename)

        implicit none

        integer(4) :: nopix
        real(fp_kind) :: image(nopix)
        character(*) :: filename

        integer :: i

        write(*,*) trim(filename)

        open(unit = 19, file = filename, status = 'unknown')
        do i = 1, nopix
            write(19,*) image(i)
        enddo

        close(19)

    end subroutine



    subroutine printout_2d(image, nopiy, nopix, filename)

        implicit none

        integer(4) :: nopix, nopiy
        real(fp_kind) :: image(nopiy,nopix)
        character(*) :: filename

        integer(4) :: ny, nx, iostat

        write(*,*) trim(filename)

        open(unit=18, file=filename, status='unknown', iostat=iostat)

        if (iostat.ne.0) then
            write(*,*) 'Error reading binary file "', trim(filename), '".'
            write(*,*) 'The program will now halt.'
            pause
            stop

        endif

        do ny = 1,nopiy
            do nx = 1, nopix-1
                    write (18,142) image(ny,nx)
142               format(1x,e20.12,$)
            enddo
            write (18,144) image(ny,nopix)
144         format(1x,e20.12)
        enddo

        close(18)

    end subroutine



    subroutine tile_output(in_image, nopiy, nopix, ifacty, ifactx, out_image)

        implicit none

        integer(4) :: nopiy,nopix
        real(fp_kind) :: in_image(nopiy,nopix)
        integer(4) :: ifacty, ifactx
        real(fp_kind) :: out_image(nopiy*ifacty,nopix*ifactx)

        integer(4) :: i, j

        do i = 1, ifacty
        do j = 1, ifactx
            out_image((i-1)*nopiy + 1 : (i-1)*nopiy + nopiy, (j-1)*nopix + 1 : (j-1)*nopix + nopix)  = in_image(1:nopiy,1:nopix)
        enddo
        enddo

    end subroutine



    subroutine interpolate_real2D(a, b)

        use FFTW3

        implicit none

        real(fp_kind)::a(:,:), b(:,:)

        complex(fp_kind),allocatable::a1(:,:), a2(:,:), b1(:,:), b2(:,:)

        integer::nopiy_a, nopix_a, nopiy_b, nopix_b

        integer::Npos_y, Npos_x, Nneg_y, Nneg_x
        integer::shifty, shiftx

        nopiy_a = size(a, 1)
        nopix_a = size(a, 2)
        nopiy_b = size(b, 1)
        nopix_b = size(b, 2)

        allocate(a1(nopiy_a,nopix_a))
        allocate(a2(nopiy_a,nopix_a))

        allocate(b1(nopiy_b,nopix_b))
        allocate(b2(nopiy_b,nopix_b))

        a1 = a

        a2 = fft(nopiy_a, nopix_a, a1,norm=.true.)
        a2 = a2 * sqrt(float(nopiy_a*nopix_a))

        Npos_y = int(floor(float(nopiy_a)/2))
        Npos_x = int(floor(float(nopix_a)/2))

        Nneg_y = nopiy_a - Npos_y - 1
        Nneg_x = nopix_a - Npos_x - 1

        shifty = -Nneg_y
        shiftx = -Nneg_x

        a1 = cshift(a2, shifty, 1)
        a2 = cshift(a1, shiftx, 2)

        b1 = 0.0_fp_kind
        b1(1:nopiy_a,1:nopix_a) = a2

        b2 = cshift(cshift(b1, -shifty, 1), -shiftx, 2)

        b1 = ifft(nopiy_b, nopix_b, b2,norm=.true.)
        b1 = b1 / sqrt(float(nopiy_b*nopix_b))

        b = real(b1) * float(product(shape(b))) / float(product(shape(a)))

    end subroutine




      !----------------------------------------------------------------------
      !   subroutine constructs a filename appending the integer to the input
      !   no extension is added
      !----------------------------------------------------------------------
      subroutine make_fnam_out_noext(output_prefix,fnam_out,num)

      implicit none

      integer(4)   i,itmp,itmp1,num
      character(*) output_prefix,fnam_out
      character*80 fnam_temp,fnam_temp1
      character*80 fnam_temp2,fnam_temp3

      fnam_temp=adjustl(output_prefix)
      write (fnam_temp1,123) num
123   format(i4)

      fnam_temp2=adjustl(fnam_temp1)
      itmp=len_trim(fnam_temp2)
      fnam_temp3=trim(fnam_temp2)
      do i = 1,5
            itmp1=itmp+i
            if(itmp1.lt.6)then
                  fnam_temp3='0'//trim(fnam_temp3)
            endif
      enddo

      fnam_out=trim(fnam_temp)//trim(fnam_temp3)

      return
      end subroutine

      subroutine add_zero_padded_int(output_prefix, fnam_out, num, width)
        implicit none

        character(*) :: output_prefix, fnam_out
        integer :: num, width

        character(10) :: fmt

10      format('(i', i1, '.', i1, ')')

        if (num.ge.0) then
            write(fmt, 10) width, width
        else
            write(fmt, 10) width, width-1
        endif

        write(fnam_out, fmt) num

        fnam_out = trim(output_prefix) // fnam_out

      end subroutine

      subroutine add_zero_padded_int_signed(output_prefix, fnam_out, num, width)
        implicit none

        character(*) :: output_prefix, fnam_out
        integer :: num, width

        character(10) :: fmt

10      format('(sp, i', i1, '.', i1, ')')
        write(fmt, 10) width, width-1

        write(fnam_out, fmt) num

        fnam_out = trim(output_prefix) // fnam_out

      end subroutine

      function defocus_string(defocus,length_df)
            use m_string
            implicit none

            character(:),allocatable :: defocus_string
            integer*4,intent(in),optional::length_df
            real(fp_kind),intent(in)::defocus

            integer*4::length_df_

            if(present(length_df))then
                length_df_=length_df
            else
                length_df_ = ceiling(log10(abs(defocus)))
                if(defocus<0) length_df_ = length_df_+1
            endif

            defocus_string = '_Defocus_'//zero_padded_int(int(defocus),length_df_)//'_Ang'

      end function
      !----------------------------------------------------------------------------
      !subroutine output stem image
      !----------------------------------------------------------------------------


    subroutine output_stem_image(stem_image,fnam_det,defoci)
        use global_variables, only:output_nopiy,output_nopix,zarray,tiley,tilex,interpolation
        use m_string, only: zero_padded_int

        implicit none

        integer(4) :: i_df,i_z,length,lengthdf,nysample,nxsample,n_df,nz
        real(fp_kind),dimension(:,:,:,:) :: stem_image
        real(fp_kind),allocatable :: tiled_image(:,:),montage(:,:)
        real(fp_kind) :: interpolated_image(output_nopiy,output_nopix)
        real(fp_kind),optional:: defoci(size(stem_image,3))
        character(*) :: fnam_det
        character(512) :: fnam_temp,fnam_out

        logical :: many_y, many_x, many_df,many_z,make_montage

        nysample = size(stem_image,1)
        nxsample = size(stem_image,2)
        n_df = size(stem_image,3)
        nz = size(stem_image,4)

        many_y = nysample .gt. 1
        many_x = nxsample .gt. 1
        many_df = n_df .gt. 1
        many_z = nz>1
        make_montage = many_df.and.many_z


        if (make_montage) allocate(montage(output_nopiy*n_df,output_nopix*nz))

        length = ceiling(log10(maxval(zarray)))
        lengthdf = ceiling(log10(maxval(abs(defoci))))
        if(any(defoci<0)) lengthdf = lengthdf+1
        if (many_y .and. many_x) then
            ! y, x, df

            do i_df = 1, n_df
            do i_z = 1, nz
                ! Output STEM image at original sampling

                fnam_out = trim(adjustl(fnam_det))

                if(many_df) fnam_out = trim(adjustl(fnam_out)) //defocus_string(defoci(i_df),lengthdf)

                if(nz>1) fnam_out = trim(adjustl(fnam_out))//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'

                if(.not.interpolation) then
                    call binary_out(nysample, nxsample, stem_image(:,:,i_df,i_z), fnam_out)
                    cycle
                endif

                ! Output STEM image with tiling and interpolation

                if(allocated(tiled_image)) deallocate(tiled_image)
                allocate(tiled_image(tiley*nysample,tilex*nxsample))

                call tile_output(stem_image(:,:,i_df,i_z), nysample, nxsample, tiley, tilex, tiled_image)

                call interpolate_real2D(tiled_image, interpolated_image)

                if(make_montage) montage((i_df-1)*output_nopiy+1:i_df*output_nopiy,&
                                         (i_z-1)*output_nopix+1:i_z*output_nopix) = interpolated_image

                call binary_out(output_nopiy, output_nopix, interpolated_image, fnam_out)
            enddo
            enddo

            if(make_montage.and.interpolation) then
                fnam_out = trim(adjustl(fnam_det)) //'_defocus_thickness_montage'
                 call binary_out(output_nopiy*n_df,output_nopix*nz, montage, fnam_out)
            endif
        elseif (many_y .and. many_df) then
            ! y vs. df
            do i_z = 1, nz
                fnam_temp = trim(adjustl(fnam_det)) // '_Linescan(y)_vs_Defocus'
                if (nz>1) fnam_temp = trim(adjustl(fnam_temp)) //'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
                call binary_out(nysample, n_df, stem_image(:,1,:,i_z), fnam_temp)
            enddo
        elseif (many_y) then
            ! y
            do i_z = 1, nz
            fnam_temp = trim(adjustl(fnam_det)) // '_Linescan(y).txt'
            if (nz>1) fnam_temp = trim(adjustl(fnam_out))//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
            call printout_1d(stem_image(:,1,1,i_z), nysample, fnam_temp)
            enddo
        elseif (many_x .and. many_df) then
            ! x vs. df
            do i_z = 1, nz
            fnam_temp = trim(adjustl(fnam_det)) // '_Linescan(x)_vs_Defocus'
            if (nz>1) fnam_temp = trim(adjustl(fnam_out))//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
            call binary_out(n_df, nxsample, transpose(stem_image(1,:,:,i_z)), fnam_temp)
            enddo
            ! Note: the transpose is so that defocus is passed through to
            ! binary_out() as the fastest varying dimension; binary_out()
            ! then transposes again so that x will be outputted as the fastest varying
            ! dimension and thus will be displayed horizontally in ImageJ

        elseif (many_x) then
            ! x
            do i_z = 1, nz
            fnam_temp = trim(adjustl(fnam_det))
            if (nz>1) fnam_temp = trim(adjustl(fnam_out))//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
            fnam_temp = trim(adjustl(fnam_out))// '_Linescan(x).txt'
            call printout_2d(stem_image(1:1,:,1,i_z), 1, nxsample, fnam_temp)
            ! Note: use printout_2d() to force x horizontal
            enddo
        elseif (many_df) then
            !df
            do i_z = 1, nz
            if (nz>1) fnam_temp = fnam_temp//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
            fnam_temp = trim(adjustl(fnam_det)) // '_DefocusSeries.txt'
            call printout_1d(stem_image(1,1,:,i_z), n_df, fnam_temp)
            enddo

        else
            !Single value
            do i_z = 1,nz
            fnam_temp = trim(adjustl(fnam_det))
            if (nz>1) fnam_temp = trim(adjustl(fnam_out))//'_z='//zero_padded_int(int(zarray(i_z)),length)//'_A'
             fnam_temp = trim(adjustl(fnam_out))// '_Signal.txt'
            call printout_1d(stem_image(:,1,1,I_z), 1, fnam_temp)
            enddo

        endif

    end subroutine

    subroutine binary_out_complex(nopiy_temp,nopix_temp,array,fnam,cutoff,realimag,quadshift)
    implicit none

    integer(4),intent(in):: nopiy_temp,nopix_temp
    complex(fp_kind),intent(in) :: array(nopiy_temp,nopix_temp)
    character*(*),intent(in) :: fnam
    real(fp_kind),intent(in),optional::cutoff
    logical,intent(in),optional::realimag,quadshift

    integer(4) i,j
    real(fp_kind), allocatable :: obj_ret(:,:)
    character*(120) fnam1
    complex(fp_kind)::array_(nopiy_temp,nopix_temp)
    logical::realimag_,mask(nopiy_temp,nopix_temp)

    array_ = array
    if(present(quadshift)) then
        if (quadshift) array_ = quad_shift(array_,nopiy_temp,nopix_temp)
    endif

    realimag_ = .false.
    if (present(realimag)) realimag_ = realimag

    if (realimag_) then
        call binary_out(nopiy_temp,nopix_temp,real(array_),trim(adjustl(fnam))//'_real.bin')
        call binary_out(nopiy_temp,nopix_temp,imag(array_),trim(adjustl(fnam))//'_imag.bin')
    else
        if(allocated(obj_ret)) deallocate(obj_ret)
        allocate(obj_ret(nopiy_temp,nopix_temp))
        obj_ret = real(conjg(array_)*array_,kind=fp_kind)
        if(present(cutoff)) then
            mask = obj_ret/maxval(obj_ret) > cutoff**2
        else
            mask = obj_ret/maxval(obj_ret) > 1e-8
        endif
        fnam1=trim(adjustl(fnam))//'_inten'
        call binary_out(nopiy_temp,nopix_temp,obj_ret,fnam1)
        obj_ret = 0
        forall(i=1:nopiy_temp,j=1:nopix_temp,mask(i,j)) obj_ret(i,j) = atan2(imag(array_(i,j)),real(array_(i,j)))
        fnam1=trim(adjustl(fnam))//'_phase'
        call binary_out(nopiy_temp,nopix_temp,obj_ret,fnam1)
    endif

    return
    end subroutine


    end module
