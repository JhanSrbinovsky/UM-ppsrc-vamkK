
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE EXPAND32B--------------------------------------
!
!    Purpose: Expands from 32 to 64 bit for dump reading routines.
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: Unified Model Documentation Paper No F3
!     ---------------------------------------------------------
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Dump I/O

SUBROUTINE expand32b(length, array, version)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 length,                                                          &
               !IN length of the field to be expanded
 version       !IN model version

REAL                                                              &
 array(length)  !IN/OUT array to be expanded in place

! -------------------------------------------------------------
! Local variables: --------------------------------------------
REAL hold(length)     ! space for expanded array
INTEGER i             ! Loop index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! -------------------------------------------------------------
!   External subroutines called:-------------------------------

IF (lhook) CALL dr_hook('EXPAND32B',zhook_in,zhook_handle)

! DEPENDS ON: expand21
CALL expand21(length,array,hold)
DO i=1,length
  array(i)=hold(i)
END DO

IF (lhook) CALL dr_hook('EXPAND32B',zhook_out,zhook_handle)
RETURN
END SUBROUTINE expand32b
