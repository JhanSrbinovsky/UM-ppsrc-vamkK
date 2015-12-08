! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!-----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!-----------------------------------------------------------------------
      SUBROUTINE DATE31 (ID, IM, IY, ICD)
!
!     DAYS SINCE 1.1.1900 FROM DAY, MONTH, YEAR
!
        IMPLICIT NONE
!                        INPUT ARGUMENTS
        INTEGER, INTENT(IN) :: ID, IM, IY
!                       OUTPUT ARGUMENTS
        INTEGER, INTENT(OUT) :: ICD
!                       LOCAL VARIABLES
        INTEGER :: IDY, INY
        INTEGER :: K,IYN, days_in_feb
        INTEGER, DIMENSION(12) :: MONTHS
!       external function
        LOGICAL, EXTERNAL :: ISALEAP

! DEPENDS ON: isaleap
        IF (ISALEAP(IY)) THEN
           days_in_feb = 29
        ELSE
           days_in_feb = 28
        ENDIF

        MONTHS = (/31,days_in_feb,31,30,31,30,31,31,30,31,30,31/)

        K = SUM(MONTHS(1:(IM-1)))    ! use array sections and intrinsics

        IDY = K + ID
        INY = IY
        IYN = INY - 1900
        IF (IYN  >   0) THEN
           ICD = IDY + IYN*365 + (IYN-1)/4 - (IYN-1)/100 + (IYN+299)/400
        ELSE
           ICD = IDY + IYN*365 + IYN/4 - IYN/100
        ENDIF


      END SUBROUTINE DATE31
!-----------------------------------------------------------------------
! The JDATE Conversion algorithms are based on the algorithm published
! in a letter to the editor of Communications of the ACM (CACM, volume 1
! number 10, October 1968, p.657) by Henry F. Fliegel and
! Thomas Van Flandern
! This algorithm is valid only for dates from
! 1/3/-4900 G onward when converting from a Julian day number to a date,
! or from 1/3/-4800 when converting from a date to a Julian day number.
! It should be noted that these algorithms are valid only in the
! Gregorian Calendar and the Proleptic Gregorian Calendar (after the
! dates given above). They do not handle dates in the Julian Calendar.
!-----------------------------------------------------------------------
