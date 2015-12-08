! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SETPERLEN---------------------------------------------
!LL
!LL  Purpose:
!LL           Return length of current meaning period using mean-level
!LL           (0, 1, 2 or 3) & current date. Days_in_month is declared
!LL           and provided by comdeck CDAYDATA.
!LL
!LL  Method: where this routine is called at the end of a month, the
!LL          month will already have been incremented, which is why
!LL          i_month-1 is used in setting the period length. Every
!LL          100 years is not leap unless year is divisible by 400.
!LL
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
!
!  Code description:
!   FORTRAN 77 + common extensions also in Fortran 90
!LL
!LL
!LL Programming standard :  UMDP 3 Version 7  (rev. 6/10/94)
!LL
!LL Logical components covered :
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LL-----------------------------------------------------------------
!*L Arguments:------------------------------------------------------
      SUBROUTINE SETPERLEN (MEANLEV,I_MONTH,I_YEAR,PERIODLEN)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE cdaydata_mod, ONLY: days_per_4c, days_per_c, days_per_4y,     &
                              days_per_y, days_in_month, days_to_month
      IMPLICIT NONE

      INTEGER                                                           &
     &  MEANLEV,                                                        &
                         ! IN - Mean level indicator, e.g. 0, 1, 2 or 3
     &  I_MONTH,                                                        &
                         ! IN - model time (month)
     &  I_YEAR,                                                         &
                         ! IN - model time (year)
     &  PERIODLEN        ! OUT - length of current meaning period (days)

!-----------------------------------------------------------------------
! Workspace usage:------------------------------------------------------
! NONE
!-----------------------------------------------------------------------
! External subroutines called:------------------------------------------
! NONE
!*----------------------------------------------------------------------
! Define local variables:-----------------------------------------------

      LOGICAL L_LEAP     ! Leap year indicator

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!  Comdecks:
!-----------------------------------------------------------------------
! End of standard header info

      IF (lhook) CALL dr_hook('SETPERLEN',zhook_in,zhook_handle)
      IF (mod(i_year,4)  ==  0 .AND.                                    &
                                              ! is this a leap year?
     &    (mod(i_year,400)  ==  0 .OR. mod(i_year,100)  /=  0)) then
        L_LEAP = .TRUE.
      ELSE
        L_LEAP = .FALSE.
      END IF

      IF (meanlev  ==  0) then       ! instantaneous data (e.g. part-
        periodlen = 1                ! way through a monthly mean)
      ELSEIF (meanlev  ==  1) then   ! end of monthly mean
        IF (L_LEAP .AND. (i_month  ==  3)) then  ! Is it leap year Feb?
          periodlen = days_in_month(i_month-1) + 1
        ELSE IF (i_month  ==  1) then  ! end of Dec, so can't use
          periodlen = 31               ! days_in_month(i_month-1)
        ELSE
          periodlen = days_in_month(i_month-1)
        END IF
      ELSE IF (meanlev  ==  2) then  ! seasonal mean

!      find season length using current month as pointer
        IF (L_LEAP) then  ! do leap year seasons
          IF (i_month  ==  5) then  ! season=FebMarApr
            periodlen = 90
          ELSE IF ((i_month  ==  3) .or. (i_month  ==  4) .or.          &
     &            (i_month  ==  7) .or. (i_month  ==  12)) then
            periodlen = 91          ! for DJF, JFM, AMJ or SON
          ELSE
            periodlen = 92
          END IF
        ELSE              ! do non-leap year seasons
          IF (i_month  ==  5) then  ! season=FebMarApr
            periodlen = 89
          ELSE IF ((i_month  ==  3) .or. (i_month  ==  4)) then
            periodlen = 90  ! for DJF and JFM
          ELSE IF ((i_month  ==  7) .or. (i_month  ==  12)) then
            periodlen = 91  ! for AMJ and SON
          ELSE
            periodlen = 92  ! for all other seasons
          END IF
        END IF  ! end of IF test of L_LEAP, and end of seasons.
      ELSE IF (meanlev  ==  3) then  ! annual mean

! Bear in mind period 3 may be 366 days if _previous_ year was leap, and
! may not always be 366 days even if current year is a leap year, since
! annual means are often not for calendar years
        IF (L_LEAP .AND. (i_month  >=  3)) then
          periodlen = 366
        ELSE IF (mod(i_year-1,4)  ==  0 .AND.                           &
                                                   ! was last year leap?
     &    (mod(i_year-1,400)  ==  0 .OR. mod(i_year-1,100)  /=  0)) then
          IF (i_month  <=  2) then   !is it no later than end of Jan?
            periodlen = 366          !Then this year had an extra day!
        ELSE
          periodlen = 365
          ENDIF                 !end of Jan
        ELSE
          periodlen = 365       ! Normal year
        ENDIF                   !was last year a leap?
      ELSE                      ! meanlev has unexpected value
        periodlen = 1           ! so set weighting factor=1
        write(6,*)'SETPERLEN: MEANLEV not in allowed range of 0 to 3'
      END IF  ! end of IF tests on meanlev

      IF (lhook) CALL dr_hook('SETPERLEN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SETPERLEN

