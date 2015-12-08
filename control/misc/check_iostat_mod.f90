! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE check_iostat_mod
! Description:
! This module is the interface to a utility subroutine
! that manages status of success/failure of 
! FORTRAN namelist READs within the UM. 
! Failures are handled in a standardised way across the UM. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

CONTAINS

!--------------------------------------------

SUBROUTINE check_iostat(errorstatus, info)

! if IOSTAT is nonzero from a READ then abort UM run.

IMPLICIT NONE  
  
INTEGER,       INTENT(IN) :: errorstatus
CHARACTER(LEN=*), INTENT(IN) :: info

CHARACTER(LEN=*), PARAMETER  :: RoutineName='check_iostat'
CHARACTER(LEN=80)              :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!Commented to prevent premature MPI initialisation by DrHook
!IF (lhook) CALL dr_hook('CHECK_IOSTAT',zhook_in,zhook_handle)

! Report fatal error (ABS(Errorstatus)) for all non-zero errors.
IF (ErrorStatus /= 0) THEN
  CMessage  = ' Error reading ' // info // &
              '. Please check input list against code.'
  CALL ereport (RoutineName, ABS(ErrorStatus), CMessage)
END IF

!Commented to prevent premature MPI initialisation by DrHook
!IF (lhook) CALL dr_hook('CHECK_IOSTAT',zhook_out,zhook_handle)
  
END SUBROUTINE check_iostat

!--------------------------------------------

END MODULE check_iostat_mod
