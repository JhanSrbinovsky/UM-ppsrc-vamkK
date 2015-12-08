! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!-----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

      LOGICAL FUNCTION ISALEAP(IY)
!
!     Returns .TRUE. if IY is a Leap year
!     Returns .FALSE. if IY is not a Leap year
!
      IMPLICIT NONE
!     INPUT ARGUMENT
      INTEGER, INTENT(IN) :: IY


      IF (IY/4*4  /=  IY) THEN    ! Divide by 4
         ISALEAP=.FALSE.
      ELSE
        IF (IY/400*400  ==  IY) THEN  ! Century check
           ISALEAP=.TRUE.
        ELSE
          IF (IY/100*100  ==  IY) THEN   ! Century qualifier
             ISALEAP=.FALSE.
          ELSE
            ISALEAP=.TRUE.
          ENDIF
        ENDIF
      ENDIF
      END FUNCTION ISALEAP

!-----------------------------------------------------------------------
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
