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
      SUBROUTINE DATE13 (ICD, ID, IM, INY)
!
!     DAY, MONTH, YEAR FROM DAYS SINCE 1.1.1900
!
        IMPLICIT NONE
!                        INPUT ARGUMENTS
        INTEGER, INTENT(IN) :: ICD
!                       OUTPUT ARGUMENTS
        INTEGER, INTENT(OUT) :: ID, IM, INY
!                       LOCAL VARIABLES
        INTEGER :: IDY, IY
        INTEGER :: K,KD,KE,KY,I,K1X, days_in_feb
        INTEGER, DIMENSION(12) :: MONTHS
!       external function
        LOGICAL, EXTERNAL :: ISALEAP

        K = ICD
        KE = 0
        IF (K  >=  366) THEN  ! these allow for the non-leap years 1900
           K = K + 1
           IF (K  >=  73416) THEN         !2100, ...
              K = K + 1
              IF (K  >=  109941) THEN     !2200,
                 K = K + 1
                 IF (K  >=  146466) THEN  !2300 ...
                    K = K + 1
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        IF (K  <=  -36159) THEN   ! and 1800 respectively
           K = K - 1
        ENDIF

        KY = K/1461*4
        KD = K - K/1461*1461
        IF (KD  <   0) THEN
           KD = KD + 1461
           KY = KY - 4
        ENDIF
        KY = KY + 1900
        IF (KD  >   366) THEN
           KD = KD - 1
           KE = KD/365
           KD = KD - KD/365*365
        ENDIF
        IF (KD  ==  0) THEN
           KE = KE - 1
           KD = 365
        ENDIF
        INY = KY + KE
        IDY = KD
        IY = INY

! DEPENDS ON: isaleap
        IF (ISALEAP(IY)) THEN
           days_in_feb = 29
        ELSE
           days_in_feb = 28
        ENDIF

        MONTHS = (/31,days_in_feb,31,30,31,30,31,31,30,31,30,31/)

        K1X = IDY

        DO I=1,12
           K1X = K1X - MONTHS(I)
           IF (K1X  >   0) THEN
              CYCLE
           ELSE
              ID = K1X + MONTHS(I)
              IM = I
              INY = IY

           ENDIF
           EXIT
        END DO

      END SUBROUTINE DATE13
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
