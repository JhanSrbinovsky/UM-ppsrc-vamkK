
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Dummy routine for converting fields from fieldsfile to GRIB.

!=======================================================================

SUBROUTINE gribify(  UMHdr_in,    &  ! in
                     Gribfile,    &  ! in
                     ErrorStatus )   ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
USE IO_Mod, ONLY:         &
  UM_Header_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN) :: UMHdr_in
CHARACTER(LEN=*), INTENT(IN)     :: gribfile
INTEGER, INTENT(INOUT)           :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "GRIBIFY"

! Local variables
CHARACTER(LEN=80) :: cmessage          ! Error Message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

CMESSAGE = 'Routine should not be callable'
ErrorStatus = 1
CALL EReport(RoutineName,ErrorStatus,CMESSAGE)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Gribify

