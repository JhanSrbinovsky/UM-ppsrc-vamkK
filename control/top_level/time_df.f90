! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: TIME_DF --------------------------------------------------
!LL
!LL  Purpose: Subroutine to obtain a new model time in days and seconds
!LL           from some reference point, given an increment in days and
!LL           seconds.  Note that the seconds and days increments are
!LL           treated independently so that -ve increments or seconds
!LL           increments larger than the no of seconds in a day are
!LL           handled correctly.
!LL           Forms a service routine for model date/time and internal
!LL           clock purposes, written for 32-bit portability.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level
      SUBROUTINE TIME_DF(DAYS1,SECS1,DEL_DAYS,DEL_SECS,DAYS2,SECS2)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
      INTEGER                                                           &
     &     DAYS1,SECS1,                                                 &
                                   ! IN  - days/seconds (input time)
     &     DEL_DAYS,DEL_SECS,                                           &
                                   ! IN  - days/seconds increments
     &     DAYS2,SECS2             ! OUT - days/seconds (output time)
!
      INTEGER                                                           &
     &     SECS_IN_DAY            ! No of seconds in a day

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      PARAMETER                                                         &
     &    (SECS_IN_DAY=24*60*60)
!
      IF (lhook) CALL dr_hook('TIME_DF',zhook_in,zhook_handle)
      DAYS2 = DAYS1 + DEL_DAYS + DEL_SECS/SECS_IN_DAY
      SECS2 = SECS1 + MOD(DEL_SECS,SECS_IN_DAY)
!
      IF(SECS2 <  0) THEN
         SECS2 = SECS2 + SECS_IN_DAY
         DAYS2 = DAYS2 - 1
      ENDIF
!
      IF(SECS2 >= SECS_IN_DAY) THEN
         SECS2 = SECS2 - SECS_IN_DAY
         DAYS2 = DAYS2 + 1
      ENDIF
!
      IF (lhook) CALL dr_hook('TIME_DF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TIME_DF
