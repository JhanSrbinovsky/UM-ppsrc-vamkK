! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given a netcdf file and a series of dimensions reads teh
!  size of those dimensions and if required check that they
!  are as expected.

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_netcdf_dimreader(        &
  filename,                                 & ! Name of netcdf file
  global_row_length,                        & ! model row length
  global_rows,                              & ! model rows
  input_ndims,                              & ! number of dimensions
  input_dimnames,                           & ! name of dimensions
  file_dimensions,                          & ! dimension lengths in file
  checkdims,                                & ! check rows/columns match (T/F)
  debug)                                      ! Debug flag

USE nudging_control
USE netcdf

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod, ONLY : ereport
USE um_types, ONLY: integer32
IMPLICIT NONE

!************************************************************************
!  Define inout variables

! Character strings for netcdf file
CHARACTER(LEN=*)        :: filename          ! Name of netcdf file
INTEGER, INTENT(IN)  :: global_rows       ! model row length
INTEGER, INTENT(IN)  :: global_row_length ! model rows
INTEGER, INTENT(IN)  :: input_ndims       ! No. dimensions read
CHARACTER(LEN=dimname_length), INTENT(IN)::                       &
  input_dimnames(input_ndims)                  ! Dimension names
INTEGER, INTENT(OUT) ::                                           &
  file_dimensions(input_ndims)            ! output dimensions
LOGICAL, INTENT(IN)  :: checkdims         ! check rows/columns match (T/F)
INTEGER, INTENT(IN)  :: debug             ! debug flag

INTEGER(KIND=integer32) :: file_dimensions_int4(totdims) ! output dimensions
INTEGER(KIND=integer32) :: status                        ! netcdf status flag
INTEGER(KIND=integer32) :: ncid                          ! netcdf file ID
INTEGER(KIND=integer32) :: dimid(totdims)                ! Array storing dim IDs
INTEGER                 :: i                             ! Loop variable

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_NETCDF_DIMREADER',zhook_in,zhook_handle)

!***********************************************************************
!End Header,
! access netcdf file to obtain dimension information

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
  ' Entered NUDGING_NETCDF_DIMREADER routine'
END IF

! Load the netcdf file (read only)
status = nf90_open(filename, nf90_nowrite, ncid)

! If there is an error set errorflag to status, close file and return
IF(status /= 0) THEN
  nmessage = 'Error in loading file' // filename

  CALL ereport('NUDGING_NETCDF_DIMREADER',status,nmessage)
END IF

IF(debug > 15) THEN
  WRITE(OUT,*)                                                         &
  'NUDGING_NETCDF_DIMREADER: Fileopen ',TRIM(filename),                &
   nf90_nowrite, ', ', ncid, ', ', status
END IF

! Loop round dimensions
DO i=1, input_ndims

! Load dimids
  status = nf90_inq_dimid(ncid, input_dimnames(i), dimid(i))

! If there is an error set errorflag to status, close file & return
  IF(status /= 0) THEN
    nmessage = 'Error in loading dimid' // input_dimnames(i)

    CALL ereport('NUDGING_NETCDF_DIMREADER',status,nmessage)
  END IF

! Load dimension lengths
  status = nf90_inquire_dimension(                                    &
           ncid, dimid(i), LEN=file_dimensions_int4(i))

! If there is an error set errorflag to status, close file & return
  IF(status /= 0) THEN
    nmessage = 'Error in loading dimension' // input_dimnames(i)

    CALL ereport('NUDGING_NETCDF_DIMREADER',status,nmessage)
  END IF

! Convert dimension lengths to int *8 for out
  file_dimensions(i) = INT(file_dimensions_int4(i), KIND=8)

  IF(debug > 10) THEN
    WRITE(OUT,*)                                                        &
     ': In NUDGING_NETCDF_DIMREADER: ',                                 &
     'DIMENSION ', input_dimnames(i), ' (', i, '), ',                   &
     'dimid and len =', dimid(i), file_dimensions(i)
  END IF

END DO     ! loop over input dims

! Close netcdf file
status = nf90_close(ncid)

! If there is an error set errorflag to status, close file and return
IF(status /= 0) THEN
  nmessage = 'Error in closing file' // filename

  CALL ereport('NUDGING_NETCDF_DIMREADER',status,nmessage)
END IF

! Can run the routine in dimension checking mode
! If we do this and the dimensions do not match then
! exit with error
IF(checkdims) THEN

! Check that the length of rows are consistent
  IF(file_dimensions(1).ne.global_row_length) THEN
    nmessage = 'Model row length differs from file'
    status = 999

    CALL ereport('NUDGING_NETCDF_DIMRDR',status,nmessage)
  END IF

! Check that the number of rows are consistent
  IF(input_ndims > 1 .AND. file_dimensions(2) /= global_rows)           &
  THEN
    nmessage = 'Model rows differs from file'
    status = 999

    CALL ereport('NUDGING_NETCDF_DIMRDR',status,nmessage)
  END IF

! Check that the length of rows are consistent
  IF(input_ndims > 4 .AND. file_dimensions(5) /= global_row_length)     &
  THEN
    nmessage = 'Model `u` row length differs from file'
    status = 999

    CALL ereport('NUDGING_NETCDF_DIMRDR',status,nmessage)
  END IF

! Check that the number of rows are consistent
  IF(input_ndims > 5 .AND. file_dimensions(6) /= (global_rows-1))       &
  THEN
    nmessage = 'Model `v` rows differs from file'
    status = 999

    CALL ereport('NUDGING_NETCDF_DIMRDR',status,nmessage)
  END IF

END IF ! Check dimensions

!*****************************************************************************
!     Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
  ' Leaving NUDGING_NETCDF_DIMREADER routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_NETCDF_DIMREADER',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_netcdf_dimreader

