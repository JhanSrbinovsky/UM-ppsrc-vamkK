! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: climate_month_date
!
! Purpose: Flux processing routine.
!          Calculates the full date of a climate field which matches
!          the input stash code and month number.
!          Also output the day of the middle of the month.
!
! WARNING: This routine contains mid_month_day_valid hard wired
!          as 15
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine climate_month_date( stcode, ValidMonth,                &
! ACLM1TIM start
! Purpose: argument list for variables storing time of earlier of two
!          climate fields used to interpolate to the validity time.
!          This deck is linked to CCLM1TIM.
!----------------------------------------------------------------------
     &  Clim1Year, Clim1Month, Clim1Day, Clim1Hour, Clim1Min, Clim1Sec, &
! ACLM1TIM end
     &       mid_month_day_valid, icode)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                             len_realhd, itgrid, iugrid,                &
                             max_num_fc_times, max_num_clim_times,      &
                             max_num_in_flux, len2_lookuppreferred,     &
                             len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses

      IMPLICIT NONE

! declaration of argument list
      integer stcode     ! IN stash code of field to look for
      integer ValidMonth ! IN month of field to look for
! validity time (in lookup header) of climate field (intent: OUT)
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
      integer mid_month_day_valid ! OUT day number of middle of month
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters

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
!----------------------------------------------------------------------
! comdeck: CLOOKUPS
! Purpose: declares and stores lookup tables for flux input files.
!          Also stores LCyclic, ISeaPt and ILandPt.
!          Some arrays are dimensioned using parameters in PLOOKUPS
!          so *CALL CLOOKUPS must always be preceded by *CALL PLOOKUPS
!----------------------------------------------------------------------
! parameters
      integer ISeaPt    ! value of land / sea mask at a sea point
                        ! ISeaPt = 0 because land / sea mask is
                        ! zero at sea points and 1 at land points
      integer ILandPt   ! ILandPt = 1
      parameter ( ISeaPt = 0, ILandPt = 1 )

! common block:
      common / Lookups /                                                &
     &   Len2_ActualPreferred, Len2_ActualPrevious, Len2_ActualClimate, &
     &   LCyclic, LCyclicO,                                             &
     &   FixHdPreferred, FixHdPrevious, FixHdClimate,                   &
     &   FixHdlsmt,FixHdlsmtO, FixHdlsmuO,                              &
     &   LookupPreferred, LookupPrevious, LookupClimate,                &
     &   LookFldNoPreferred, LookFldNoPrevious, LookFldNoClimate,       &
     &   Lookuplsmt, LookuplsmtO, LookuplsmuO


! actual numbers of fields in files
      integer Len2_ActualPreferred ! for NWP file
      integer Len2_ActualPrevious  ! for NWP file
      integer Len2_ActualClimate   ! for climate file

! additional information on fields read in
      logical LCyclic   ! T => atmosphere grid is cyclic
      logical LCyclicO  ! T => ocean grid is cyclic

! fixed headers
      integer FixHdPreferred(Len_FixHd)  ! for preferred NWP file
      integer FixHdPrevious(Len_FixHd)   ! for previous NWP file
      integer FixHdClimate(Len_FixHd)    ! for climate file
      integer FixHdlsmt(Len_FixHd)       ! for atmos. tracer lsm
      integer FixHdlsmtO(Len_FixHd)      ! for ocean tracer lsm
      integer FixHdlsmuO(Len_FixHd)      ! for ocean velocity lsm

! lookup tables for NWP preferred and previous files and climate files
      integer LookupPreferred(Len1_Lookup, Len2_LookupPreferred)
      integer LookupPrevious(Len1_Lookup, Len2_LookupPrevious)
      integer LookupClimate(Len1_Lookup, Len2_LookupClimate)

! additional lookup entry indicating the field (and lookup table)
! to use to access data with this validity time and stash code.
! This allows the user to specify copies of data with amended dates to
! be used. (See subroutine add_lookups)
      integer LookFldNoPreferred(Len2_LookupPreferred)
      integer LookFldNoPrevious(Len2_LookupPrevious)
      integer LookFldNoClimate(Len2_LookupClimate)

! lookup tables for land / sea masks for 4 grids:
      integer Lookuplsmt(Len1_Lookup)   ! atmosphere tracer
!  Lookuplsmu(Len1_Lookup)   is not used from vn 5.3
      integer LookuplsmtO(Len1_Lookup)  ! FOAM ocean tracer
      integer LookuplsmuO(Len1_Lookup)  ! FOAM ocean velocity
!----------------------------------------------------------------------

! no local arrays

! declaration of local scalars
      logical ItemFound     ! T => item has been found
      integer fld_no        ! number of lookup table of required field
      integer i             ! do loop index for lookup table number
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'climate_month_date'  ! subroutine name for error messages

! 1. Find a match in the climate lookup tables

      ItemFound = .false.
      do i = 1, Len2_ActualClimate
        if ( LookupClimate(LBMON,i)      ==  ValidMonth  .and.          &
     &       LookupClimate(ITEM_CODE,i)  ==  StCode ) then
          ItemFound = .True.
          fld_no = i
          go to 100
        end if
      end do

100   continue

      if ( .not. ItemFound ) then
        icode = 34
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1. unable to find climate field with stcode ', stcode,   &
     &  ' for month ', ValidMonth
        go to 9999
      end if


! 2. Set the date in CCLM1TIM from the lookup table
      Clim1Year  = LookupClimate(LBYR, fld_no)
      Clim1Month = LookupClimate(LBMON, fld_no)
      Clim1Day   = LookupClimate(LBDAT, fld_no)
      Clim1Hour  = LookupClimate(LBHR, fld_no)
      Clim1Min   = LookupClimate(LBMIN, fld_no)
      Clim1Sec   = 0

! 3. Calculate the middle day in the month from the lookup table

!     if ( LookupClimate(LBDAT,  fld_no)  == 
!    #     LookupClimate(LBDATD, fld_no)       ) then
!        mid_month_day_valid = LookupClimate(LBDATD, fld_no)

!     else
!        mid_month_day_valid = 0.5 * (LookupClimate(LBDATD,fld_no) + 1)

!     end if

      mid_month_day_valid = 15
      if ( mid_month_day_valid  <   14 .or.                             &
     &     mid_month_day_valid  >   16       ) then
        icode = 35
        write(UnWarn,*)CErr,CSub,                                       &
     &  ' step 3. Lookup table times for climate fields are strange ',  &
     &  ' mid_month_day_valid = ', mid_month_day_valid
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE climate_month_date
!----------------------------------------------------------------------
