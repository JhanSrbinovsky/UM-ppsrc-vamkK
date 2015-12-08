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

      SUBROUTINE  ZPDATE
!
!     Prints version information
!
      IMPLICIT NONE
      PRINT *, '  ZPDATE - F90 fixed format module-free version         &
     &              (Y2K Compliance Checked)'
      PRINT *, '  LAST MODIFIED MONDAY 5th October 1998'
      PRINT *, '        by Stephen Turner (DD)'
      PRINT *, '  Contact Software Engineering Group with any queries.'
      RETURN
      END SUBROUTINE ZPDATE

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

