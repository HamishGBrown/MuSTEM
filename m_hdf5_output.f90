! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! This example shows how to write and read a hyperslab.
!

module m_hdf5_output
use m_precision
use hdf5 ! This module contains all necessary modules
contains
	
	function remove_extension(fnam)
        character*(*),intent(in)::fnam
		character*(*)::remove_extension
        integer*4:: indx
        indx = index(fnam,'.')
		remove_extension = fnam(:indx-1)
    
    end function    
	
     subroutine output_to_hdf5(array,fnam)

     

     IMPLICIT NONE
	
	 character*(*),intent(in)::fnam
	 real(fp_kind),intent(in)::array(:,:,:,:)
	 
	 character*(*)::filename

     ! CHARACTER(LEN=7), PARAMETER :: filename = "sdsf.h5"  ! File name
     CHARACTER(LEN=5), PARAMETER :: dsetname = "Array" ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier
     INTEGER(HID_T) :: memspace      ! memspace identifier

     ! INTEGER(HSIZE_T), DIMENSION(3) :: dimsm = (/7,7,3/) ! Dataset dimensions
                                                         ! ! in memory
     INTEGER(HSIZE_T), DIMENSION(4) :: dimsf != (/5,6/) ! Dataset dimensions.

     INTEGER(HSIZE_T), DIMENSION(2) :: count = (/3,4/)
                                            ! Size of the hyperslab in the file
     INTEGER(HSIZE_T), DIMENSION(2) :: offset = (/1,2/)
                                            !hyperslab offset in the file
     INTEGER(HSIZE_T), DIMENSION(3) :: count_out = (/3,4,1/)
                                            !Size of the hyperslab in memory
     INTEGER(HSIZE_T), DIMENSION(3) :: offset_out = (/3,0,0/)
                                            !hyperslab offset in memory
     INTEGER, DIMENSION(5,6) :: data ! Data to write
     INTEGER, DIMENSION(7,7,3) :: data_out ! Output buffer
     INTEGER :: dsetrank = 4 ! Dataset rank ( in file )
     ! INTEGER :: memrank = 3  ! Dataset rank ( in memory )
     INTEGER :: i, j, k

     INTEGER :: error  ! Error flag
     INTEGER(HSIZE_T), DIMENSION(4) :: data_dims

	
	
	!Get rank and size of input dataset
	data_dims = shape(array)
     !
     ! Initialize FORTRAN interface.
     !
     CALL h5open_f(error)
	  
	 filename = remove_extension(fnam)
     !
     ! Create a new file using default properties.
     !
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

     !
     ! Create the data space for the  dataset.
     !
     CALL h5screate_simple_f(dsetrank, dimsf, dataspace, error)

     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, dsetname, H5T_FLOAT_F, dataspace, &
                      dset_id, error)

     !
     ! Write the dataset.
     !
     CALL h5dwrite_f(dset_id, H5T_FLOAT_F, array, data_dims, error)

     !
     ! Close the dataspace for the dataset.
     !
     CALL h5sclose_f(dataspace, error)

     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, error)

     !
     ! Close the file.
     !
     CALL h5fclose_f(file_id, error)

     !
     ! Close FORTRAN interface.
     !
     CALL h5close_f(error)

     END subroutine
end module