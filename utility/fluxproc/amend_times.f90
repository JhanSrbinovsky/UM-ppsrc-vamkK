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
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine amend_times (                                          &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                      Int_head,Len_IntHd )

      USE lookup_addresses

      IMPLICIT NONE

! declaration of parameters

! declaration of globals used
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
!----------------------------------------------------------------------
! comdeck CVALOFF
! Purpose: stores user inputs relating to validity times and units to
!          use for output of FOAM fluxes
!----------------------------------------------------------------------

      ! maximum number of validity times
      INTEGER,PARAMETER:: MaxTimes = 50

      ! period in hours covered by each of output flux fields
                              ! (6 hours for operational FOAM)
      INTEGER :: ValidityPeriod

      ! number of (validity) times to process
      INTEGER :: NoValidTimes

      ! offset from reference time of validity time of end of flux
      ! period (in hours)
      INTEGER :: IValidOffHr(MaxTimes)

      ! offset from main output unit. Used to output fields for last
                                    ! validity time to separate files
      INTEGER :: IOutUnitOff(MaxTimes)

! control for inserting additional "copies" of lookup tables / fields

      ! For preferred file
      INTEGER:: NoAddTimesPreferred ! number of lookup tables to insert
      INTEGER:: ISrchOffHrPreferred(MaxTimes) ! offset hours to look for
      INTEGER:: INewOffHrPreferred(MaxTimes) ! new offset hours

      ! For previous file
      INTEGER:: NoAddTimesPrevious ! number of lookup tables to insert
      INTEGER:: ISrchOffHrPrevious(MaxTimes) ! offset hour to look for
      INTEGER:: INewOffHrPrevious(MaxTimes) ! new offset hour

      ! For climate file
      INTEGER:: NoAddTimesClimate ! number of lookup tables to insert
      INTEGER:: ISrchOffHrClimate(MaxTimes) ! offset hour to look for
      INTEGER:: INewOffHrClimate(MaxTimes) ! new offset hour

      ! Value to use at land points in final output file (the default
      ! value is rmdi; for testing it is sometimes useful  to set this
      ! value to zero).
      REAL :: output_land_value   ! value at land points

      COMMON / ValOff / ValidityPeriod,                                 &
     &  NoValidTimes, IValidOffHr, IOutUnitOff,                         &
     &  NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,   &
     &  NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,      &
     &  NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,         &
     &  output_land_value

! CVALOFF end

! declaration of argument list
      integer Len_IntHd
      integer Int_Head(Len_IntHd) ! IN/OUT integer part of lookup table

! declarations for validity time

      integer IAddHrs !    used to find validity time
      integer Year1   !    First year in header
      integer Month1  !    First Month in header
      integer Day1    !    First Day in header
      integer Hour1   !    First Hour in header
      integer Min1    !    Always equal to zero
      integer Sec1    !    Always equal to zero

! no other variables used

! declaration of externals
      external add_hours
!----------------------------------------------------------------------

! 1. Set the second time in header
      Int_Head(LBYRD) = ValidYear
      Int_Head(LBMOND) = ValidMonth
      Int_Head(LBDATD) = ValidDay
      Int_Head(LBHRD) = ValidHour

! 2. Set the first time in header

      IAddHrs = 0 - ValidityPeriod

! DEPENDS ON: add_hours
      call add_hours(                                                   &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Year1,Month1,Day1,Hour1,Min1,Sec1,                         &
     &       IAddHrs )

      Int_Head(LBYR)  = Year1
      Int_Head(LBMON) = Month1
      Int_Head(LBDAT) = Day1
      Int_Head(LBHR)  = Hour1

      return
      END SUBROUTINE amend_times
!----------------------------------------------------------------------
