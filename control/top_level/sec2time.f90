! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SEC2TIM0 -------------------------------------------------
!LL
!LL  Purpose: Converts from an integer number of elapsed seconds since
!LL           the model basis time to a calendar date/time, using the
!LL           absolute calendar zero point as a reference.  30-day
!LL           month or standard calendar may be used.
!LL           NB: BASIS_TIME_SECS is the number of seconds from the
!LL           calendar zero point to the basis time for the run, and
!LL           is calculated in INITTIME.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S62
!LL
!LL  Project task: S62
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE SEC2TIME(ELAPSED_DAYS,ELAPSED_SECS                     &
     &,                   BASIS_TIME_DAYS,BASIS_TIME_SECS               &
     &,                   I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND &
     &,                   I_DAY_NUMBER,LCAL360)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE cdaydata_mod, ONLY: days_per_4c, days_per_c, days_per_4y,     &
                              days_per_y, days_in_month, days_to_month
      IMPLICIT NONE
      LOGICAL LCAL360
!
      INTEGER                                                           &
     &     ELAPSED_DAYS,                                                &
                                   ! IN  - elapsed days since basis time
     &     ELAPSED_SECS,                                                &
                                   ! IN  - elapsed secs in part of day
     &     BASIS_TIME_DAYS,                                             &
                                   ! IN  - whole days to basis time
     &     BASIS_TIME_SECS,                                             &
                                   ! IN  - secs in day at basis time
!                                  !       relative to calendar zero
     &     I_SECOND,                                                    &
                                   ! OUT - model time (seconds)
     &     I_MINUTE,                                                    &
                                   ! OUT - model time (minutes)
     &     I_HOUR,                                                      &
                                   ! OUT - model time (hours)
     &     I_DAY,                                                       &
                                   ! OUT - model time (days)
     &     I_MONTH,                                                     &
                                   ! OUT - model time (months)
     &     I_YEAR,                                                      &
                                   ! OUT - model time (years)
     &     I_DAY_NUMBER            ! OUT - model time (day number)

!
!*----------------------------------------------------------------------
!
!  Common blocks
!
!
!  Local variables
!
      INTEGER                                                           &
     &       SECOND      ! number of seconds since calendar zero

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!L----------------------------------------------------------------------
!L 1. Add elapsed time to basis time in days/seconds to get elapsed
!L     since calendar zero, and convert to hours, minutes, seconds and
!L     total days since calendar zero
!L
      IF (lhook) CALL dr_hook('SEC2TIME',zhook_in,zhook_handle)
      SECOND   = BASIS_TIME_SECS+ELAPSED_SECS
      I_DAY    = BASIS_TIME_DAYS+ELAPSED_DAYS+SECOND/86400
      IF (SECOND  >=  0) THEN
        SECOND = MOD(SECOND,86400)
      ELSE
        SECOND = MOD(SECOND,86400)
       ! Correction for fact that rounding was upwards.
        IF (SECOND  /=  0) THEN
          I_DAY  = I_DAY - 1
          SECOND = SECOND + 86400
        END IF
      END IF
      I_HOUR   = MOD(SECOND/3600,24)
      I_MINUTE = MOD(SECOND/60  ,60)
      I_SECOND = MOD(SECOND,60)
!L----------------------------------------------------------------------
!L 2. Convert day number to date
!L
      IF (LCAL360) THEN
!L
!L 2.1 30-day month (360 day year) calendar
!L
      I_YEAR  = I_DAY/360
      I_MONTH = MOD(I_DAY/30,12)+1
      I_DAY   = MOD(I_DAY,30)+1
      I_DAY_NUMBER = I_DAY+30*(I_MONTH-1)

      ELSE
!L
!L 2.2 Gregorian calendar
!L
      I_YEAR = (I_DAY/DAYS_PER_4C)*400
      I_DAY = I_DAY-(I_DAY/DAYS_PER_4C)*DAYS_PER_4C
!L     Catch special case 31 Dec in leap years
      IF (I_DAY == 4*DAYS_PER_C) THEN
        I_YEAR = I_YEAR+400
        I_DAY = DAYS_PER_Y+1
      ELSE
        I_YEAR = I_YEAR+(I_DAY/DAYS_PER_C)*100
        I_DAY = I_DAY-(I_DAY/DAYS_PER_C)*DAYS_PER_C
        I_YEAR = I_YEAR+(I_DAY/DAYS_PER_4Y)*4
        I_DAY = I_DAY-(I_DAY/DAYS_PER_4Y)*DAYS_PER_4Y
        IF (I_DAY == 4*DAYS_PER_Y) THEN
          I_YEAR = I_YEAR+4
          I_DAY = DAYS_PER_Y+1
        ELSE
          I_YEAR = I_YEAR+(I_DAY/DAYS_PER_Y) + 1
          I_DAY = I_DAY-(I_DAY/DAYS_PER_Y)*DAYS_PER_Y + 1
        ENDIF
      ENDIF
      I_DAY_NUMBER = I_DAY
!L     Find month/day from day no in year
      I_MONTH = 1
      DO WHILE ((I_MONTH  <=  12) .and.                                 &
     &  (I_DAY  >   DAYS_IN_MONTH(I_MONTH)))
        I_DAY = I_DAY-DAYS_IN_MONTH(I_MONTH)
        I_MONTH = I_MONTH+1
      ENDDO
!L     Adjust if leap year and after February
      IF (I_MONTH >  2 .AND. MOD(I_YEAR,4) == 0 .AND.                   &
     &    (MOD(I_YEAR,400) == 0 .OR. MOD(I_YEAR,100) /= 0)) THEN
        I_DAY = I_DAY-1
        IF (I_DAY == 0) THEN
          I_MONTH = I_MONTH-1
          I_DAY = DAYS_IN_MONTH(I_MONTH)
          IF (I_MONTH == 2) I_DAY=29
        ENDIF
      ENDIF
      END IF  !  LCAL360
!L----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SEC2TIME',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SEC2TIME
