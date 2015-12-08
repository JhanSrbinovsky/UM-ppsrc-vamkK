! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL  Purpose : Read from ACOBS Files,reformat and place OBS header
!LL            details in COMOBS. The bulk of the required OBS data
!LL            is put into dynamic work array OBS for transmission via
!LL            argument list to GETOBS. OBS is written out to a cache
!LL            file for subsequent reading at later timesteps.
!LL            Thus reread of ACOBS files only required intermittently
!LL            (The routine DAYS does a dd/mm/yy to dayno)
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE days_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE DAYS (KDAY,KMONTH,KYEAR,KDAYS)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER :: kday
      INTEGER :: kmonth
      INTEGER :: kyear
      INTEGER :: kdays
      INTEGER :: ileap

!L           FIND ELAPSED DAYS SINCE 1 JAN 1980.
      INTEGER IDAYS(12,2)


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      DATA IDAYS/0,31,59,90,120,151,181,212,243,273,304,334,            &
     &           0,31,60,91,121,152,182,213,244,274,305,335/

      IF (lhook) CALL dr_hook('DAYS',zhook_in,zhook_handle)

      ILEAP=1
      IF(MOD(KYEAR,4) == 0)ILEAP=2
      KDAYS=(KYEAR-1980)*365+IDAYS(KMONTH,ILEAP)+KDAY+(KYEAR-1977)/4

      IF (lhook) CALL dr_hook('DAYS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DAYS
END MODULE days_mod
