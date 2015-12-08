! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: add_hours,amend_times
!
! Purpose:  Flux processing routines.
!           add_hours adds an input number of hours
!           to an input date/time (called Clim1)
!           to produce an output date/time (called Valid).
!           amend_times amends date information in Int_Head.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine add_hours(                                             &
! ACLM1TIM start
! Purpose: argument list for variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to CCLM1TIM.
!----------------------------------------------------------------------
     &  Clim1Year, Clim1Month, Clim1Day, Clim1Hour, Clim1Min, Clim1Sec, &
! ACLM1TIM end
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       IAddHrs )

      implicit none

! declaration of argument list

! input time: intent IN
!----------------------------------------------------------------------
! comdeck: CCLM1TIM
! Purpose: declares local variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to ACLM1TIM.
!----------------------------------------------------------------------
! declarations:
! local variables defining a time
      integer Clim1Year, Clim1Month, Clim1Day, Clim1Hour,               &
     &        Clim1Min, Clim1Sec
!----------------------------------------------------------------------
! output time: intent OUT
!----------------------------------------------------------------------
! comdeck: CVALTIM
! Purpose: declares local variables storing a time of validity
!          This deck is linked to AVALTIM.
!----------------------------------------------------------------------
! declarations:
! local variables defining a time of validity
      integer ValidYear, ValidMonth, ValidDay, ValidHour,               &
     &        ValidMin, ValidSec
!----------------------------------------------------------------------

      integer IAddHrs  ! number of hours to add to input time

! no parameters, globals or local arrays used

! declaration of local scalars
      integer CDay  ! century day (of input)
      integer CHour ! century-hour
      integer New_CHour ! new century-hour (after adding input hours)
      integer New_CDay ! new century day

      external date31, date13

!----------------------------------------------------------------------

! 1. Convert the input date to century-day
! DEPENDS ON: date31
      call date31(Clim1Day, Clim1Month, Clim1Year,CDay)

! 2. Calculate century-hour
      CHour = (CDay-1)*24 + Clim1Hour

! 3. Calculate new century-hour
      New_CHour = CHour + IAddHrs

! 4. Calculate new century-day
      New_CDay = 1 + New_Chour/24

! 5. Convert new century-day to new date
! DEPENDS ON: date13
      call date13(New_CDay,ValidDay,ValidMonth,ValidYear)

! 6. Convert rest of values
      ValidHour = New_CHour - 24 * (New_CDay - 1)
      ValidMin = Clim1Min
      ValidSec = Clim1Sec

      return
      END SUBROUTINE add_hours
!----------------------------------------------------------------------
!----------------------------------------------------------------------
