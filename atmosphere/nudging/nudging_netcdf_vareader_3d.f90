! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given a netcdf file and variable name return 3D variable

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_VARLOADER.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_netcdf_vareader_3d(   &
  filename,                              &  ! Filename
  varname,                               &  ! Variable name
  dim1,                                  &  ! Dimension 1
  dim2,                                  &  ! Dimension 2
  dim3,                                  &  ! Dimension 3
  variable,                              &  ! Return variable
  debug)                                    ! Debug flag

USE nudging_control
USE netcdf

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

  USE ereport_mod, ONLY : ereport
  USE um_types, ONLY: integer32
  IMPLICIT NONE

! Filename, variable name, dimension names
CHARACTER(LEN=*):: filename             ! name of netcdf file
CHARACTER(LEN=*):: varname              ! name of variable

! Load in the dimensions of the netcdf file
INTEGER, INTENT(IN) :: dim1          ! dimension 1 length
INTEGER, INTENT(IN) :: dim2          ! dimension 2 length
INTEGER, INTENT(IN) :: dim3          ! dimension 3 length

REAL, INTENT(OUT)   :: variable(dim1, dim2, dim3)

INTEGER, INTENT(IN) :: debug       ! debug flag

!*********************************************************************
! Various netcdf flags and identifiiers
INTEGER(KIND=integer32) :: status                ! netcdf status flag
INTEGER(KIND=integer32) :: ncid                  ! netcdf file ID
INTEGER(KIND=integer32) :: varid                 ! variable ID
INTEGER(KIND=integer32) :: numdim                ! number of dimensions
INTEGER(KIND=integer32) :: mydimids (3)          ! array of dimension IDs
INTEGER(KIND=integer32) :: dimlen (3)            ! file dimension lengths
INTEGER                 :: dims (3)              ! expected dimension lengths
INTEGER                 :: i! Loop variable

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_NETCDF_VAREADER_3D',zhook_in,zhook_handle)

!***********************************************************************
! End Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
   ': NUDGING_NETCDF_VAREADER (3D): Entered routine'
END IF

!************************************
! Open netcdf file and obtain the variable ID
status = nf90_open(filename, nf90_nowrite, ncid)

! If there is an error set errorflag to status, close file and return
IF(status /= 0) THEN
  nmessage = 'Error in loading file' // filename

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

!******************************
! Load the variable ID
status = nf90_inq_varid(ncid, varname, varid)

! If there is an error set errorflag to status, close file and return
IF(status /= 0) THEN
  nmessage = 'Error in loading Variable' // varname

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

!*****************************
! Load the number of dimensions
status = nf90_inquire_variable(ncid, varid, ndims=numdim)

! If there is an error set errorflag to status and return
IF(status /= 0) THEN
  nmessage = 'error IN loading '// varname //' num dimensions'

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

! Check that no. dimensions is correct. If not return error
IF(numdim /= 3) THEN
  nmessage = 'Error Number dimensions NE 3' //varname
  status = 999

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

!****************************
! Check that we have the correct dimension lengths

! Obtain length of dimensions from netcdf file
 status = nf90_inquire_variable(ncid, varid, dimids=mydimids)

IF(status /= 0) THEN
  nmessage = 'Error Loading dim IDS' //varname

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

!Construct array for input dimensions
dims = (/dim1, dim2, dim3/)

! oop  over dimensions check input and netcdf dimensions are the same
DO i=1, numdim

! Load the length of dimension from the netcdf file
  status = nf90_inquire_dimension(ncid, mydimids(i), LEN=dimlen(i))

! Check that dimension lengths have been loaded correctly
  IF(status /= 0) THEN
    nmessage = 'Error Loading dimension' //  dim_names(i)

    CALL ereport('NUDGING_NETCDF_VAREADER_4D',status,nmessage)
  END IF

! Check that dimension length is as expected
  IF(dimlen(i) /= dims(i)) THEN
    nmessage = 'Error Dimension '//dim_names(i)                        &
    //' is incorrect length'
    status = 999

    CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
  END IF  ! dimlen equals dims

END DO !dimensions

!***********************************************************************
! If everything checks out then load the variable

! Write out various input to track programm if debug switch requires it
IF(debug > 10) WRITE(OUT,*)                                            &
': NUDGING_NETCDF_VAREADER (3D): about to fill variable ',             &
 varid, ', dimension lengths = ',                                      &
 dim1, dim2, dim3

! Load the Variable
status = nf90_get_var(ncid, varid, variable)

! Check that variable has been loaded correctly
IF(status /= 0) THEN
  nmessage = 'Error Loading Variable ' // varname

  CALL ereport('NUDGING_NETCDF_VAREADER_3D',status,nmessage)
END IF

!********************************************************************
! Close netcdf file
status = nf90_close(ncid)

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
   ' NUDGING_NETCDF_VAREADER (3D): Leaving routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_NETCDF_VAREADER_3D',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_netcdf_vareader_3d

