! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: set_climate_times
!
! Purpose: Flux processing routine.
!          Preliminaries for interpolating climate fields in time
!          sets date/times required to extract two climate fields
!          and calculates the weights to give to them
!
! WARNING: does not test ICODE properly yet
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine set_climate_times ( stcode,                            &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
! ACLM1TIM start
! Purpose: argument list for variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to CCLM1TIM.
!----------------------------------------------------------------------
     &  Clim1Year, Clim1Month, Clim1Day, Clim1Hour, Clim1Min, Clim1Sec, &
! ACLM1TIM end
!----------------------------------------------------------------------
! comdeck: ACLM2TIM
! Purpose: argument list for variables storing time of later of two
!          climate fields used to interpolate to the Clim2ity time.
!          This deck is linked to CCLM2TIM.
!----------------------------------------------------------------------
     & Clim2Year, Clim2Month, Clim2Day, Clim2Hour, Clim2Min, Clim2Sec,  &
!----------------------------------------------------------------------
     &     weight1, weight2, icode  )

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      implicit none

! declaration of argument list

      integer stcode ! IN stash code of field being accessed
! validity time of field (intent: IN)
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
! validity time (in lookup headers) of climate fields (intent: OUT)
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
!----------------------------------------------------------------------
! comdeck: CCLM2TIM
! Purpose: declares local variables storing time of later of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to ACLM2TIM.
!----------------------------------------------------------------------
! declarations:
! local variables defining a time
      integer Clim2Year, Clim2Month, Clim2Day, Clim2Hour,               &
     &        Clim2Min, Clim2Sec
!----------------------------------------------------------------------
      real weight1   ! OUT  weight to give to climate field 1
      real weight2   ! OUT  weight to give to climate field 2
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! no parameters

! declaration of globals used
! 5.3      22/10/01     Update input file units. A. Hines
! 6.1      30/07/04     Update input file units for UM6.0. A. Hines
! CUNITNOS defines parameters defining unit numbers and logicals for
! units of fields which are open

      ! for input control files
      INTEGER, PARAMETER :: UnitDbg = 7
      INTEGER, PARAMETER :: UnitHK  = 8
      INTEGER, PARAMETER :: UnitVT  = 9
      INTEGER, PARAMETER :: UnitSlt = 94

      ! for flux fields
      INTEGER, PARAMETER :: UnitPreferred = 10
      INTEGER, PARAMETER :: UnitPrevious  = 11
      INTEGER, PARAMETER :: UnitClimate   = 12

      ! for land / sea masks
      INTEGER, PARAMETER :: UnitNWPlsmt  = 13 ! NWP tracer grid
      INTEGER, PARAMETER :: UnitFOAMlsmt = 15 ! FOAM tracer grid
      INTEGER, PARAMETER :: UnitFOAMlsmu = 16 ! FOAM velocity grid

      ! for output of ancillary file validity times namelist
      ! output ancillary file validity times
      INTEGER, PARAMETER :: UnitVTOut = 18

      ! for main set of output pp files
      ! wind stress & mixing energy
      INTEGER,PARAMETER:: UnitWindsOut     = 20

      INTEGER,PARAMETER:: UnitHeatOut      = 21 ! heat fluxes
      INTEGER,PARAMETER:: UnitMoistureOut  = 25 ! moisture fluxes
      INTEGER,PARAMETER:: UnitSeaIceOut    = 23 ! sea-ice field
      INTEGER,PARAMETER:: UnitReferencesOut= 24 ! reference fields
      INTEGER,PARAMETER:: UnitPressureOut  = 26 ! pressure fields
      INTEGER,PARAMETER:: UnitWindspdOut   = 27 ! wind speed

! for output units
      INTEGER,PARAMETER:: IUnOutLow        = 20
      INTEGER,PARAMETER:: IUnOutHi         = 44

      COMMON / OutUnits / LUnOutOpen, LPreferred, LPrevious, LClimate

      !  logical T => output ! unit is open
      LOGICAL :: LUnOutOpen(IUnOutLow:IUnOutHi)

      LOGICAL :: LPreferred ! T => preferred NWP file is available
      LOGICAL :: LPrevious  ! T => previous NWP file is available
      LOGICAL :: LClimate   ! T => climate file is available
! CUNITNOS end

! no local arrays

! declaration of local scalars
! mid_month_day_# is day number of middle of month determined from
! lookup tables of climate file
      integer mid_month_day_valid  ! for validity time
      integer mid_month_day_clim1  ! for climate field 1
      integer mid_month_day_clim2  ! for climate field 2
      integer Year1     ! year for climate field 1
      integer Year2     ! year for climate field 2
      integer CDay      ! century day
      integer C1Hour    ! century hour of climate field 1
      integer C2Hour    ! century hour of climate field 2
      integer ValHour   ! century hour of validity time

      external climate_month_date, date31
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'set_climate_times'  ! subroutine name for error messages

! 1. calculate the mid-month day of the validity time

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, ValidMonth,                      &
! ACLM1TIM start
! Purpose: argument list for variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to CCLM1TIM.
!----------------------------------------------------------------------
     &  Clim1Year, Clim1Month, Clim1Day, Clim1Hour, Clim1Min, Clim1Sec, &
! ACLM1TIM end
     &       mid_month_day_valid, icode )


! 2. determine the months of the first and second climate fields to use

      if (  ValidDay  >   mid_month_day_valid ) then
        Clim1Month = ValidMonth
      else
        Clim1Month = ValidMonth - 1
      end if

      if ( Clim1Month  ==  0) then
        Clim1Month = 12
      end if

      Clim2Month = Clim1Month + 1

      if ( Clim2Month  ==  13) then
        Clim2Month = 1
      end if

! 3. find mid-month days of the first and second climate months
!    and the full dates in the lookup tables for these fields (one
!    of the main outputsfrom this routine)

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, Clim1Month,                      &
! ACLM1TIM start
! Purpose: argument list for variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to CCLM1TIM.
!----------------------------------------------------------------------
     &  Clim1Year, Clim1Month, Clim1Day, Clim1Hour, Clim1Min, Clim1Sec, &
! ACLM1TIM end
     &       mid_month_day_clim1,icode )

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, Clim2Month,                      &
!----------------------------------------------------------------------
! comdeck: ACLM2TIM
! Purpose: argument list for variables storing time of later of two
!          climate fields used to interpolate to the Clim2ity time.
!          This deck is linked to CCLM2TIM.
!----------------------------------------------------------------------
     & Clim2Year, Clim2Month, Clim2Day, Clim2Hour, Clim2Min, Clim2Sec,  &
!----------------------------------------------------------------------
     &       mid_month_day_clim2,icode )

! 4. find the weights to give two months when interpolating to
!    validity time

! 4.1 find the years for months 1 and 2

      if ( Clim1Month  ==  12 .and. ValidMonth  ==  1 ) then
        Year1 = ValidYear - 1
      else
        Year1 = ValidYear
      end if

      if ( Clim2Month  ==  1 .and. ValidMonth  ==  12 ) then
        Year2 = ValidYear + 1
      else
        Year2 = ValidYear
      end if

! 4.2 find the relative times (in hours) of the three dates

! DEPENDS ON: date31
      call date31(mid_month_day_clim1, Clim1Month, Year1,CDay)
      C1Hour = (CDay-1)*24

! DEPENDS ON: date31
      call date31(mid_month_day_clim2, Clim2Month, Year2,CDay)
      C2Hour = (CDay-1)*24

! DEPENDS ON: date31
      call date31(ValidDay, ValidMonth, ValidYear,CDay)
      ValHour = (CDay-1)*24 + ValidHour

! 4.3 calculate the weights
      weight1 = real( C2Hour - ValHour ) / real( C2Hour - C1Hour )
      weight2 = 1.0 - weight1

9999  continue
      return
      END SUBROUTINE set_climate_times
!----------------------------------------------------------------------
