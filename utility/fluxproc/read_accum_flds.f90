! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_accum_flds
!
! Purpose:  Flux processing routine.
!           Reads fields for validity time and validity time minus
!           six hours. These two fields are then manipulated to
!           obtain an accumulation for the six hour period.
!
! Uses:     StCode and to read NWP files;
!           xstcode to read climate fields
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_accum_flds(StCode, IVTOffHr,                      &
     &               ldebug, Int_Head, Real_Head,                       &
     &               ncols, nrows, field,                               &
     &               icode)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                              len_realhd, itgrid, iugrid,               &
                              max_num_fc_times, max_num_clim_times,     &
                              max_num_in_flux, len2_lookuppreferred,    &
                              len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses
      USE Submodel_Mod

      IMPLICIT NONE

! declaration of parameters
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

      real time            ! timescale in seconds for division
!(SAS)**** parameter ( time = 21600 )

! declaration of argument list

! search criteria

!       Uses    StCode to read NWP files
!               stcode to read climate fields
      integer StCode       ! IN StCode value to test

!       Reference date is used with IVTOffHr to define validity
!       time needed
      integer IVTOffHr     ! IN offset from validity time in hours

! debug control variable
      logical ldebug          ! IN T => output debugging info
      logical l_climate_field ! Set to false initially

! lookup tables
      integer Int_Head(Len_IntHd) ! OUT
      real Real_Head(Len_RealHd)  ! OUT

! output field
      integer ncols             ! IN  number of columns
      integer nrows             ! IN  number of rows
      real field(ncols,nrows)   ! OUT field values

! error code
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
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

!----------------------------------------------------------------------
! comdeck: CREFTIM
! Purpose: declares and stores a reference time
!          This deck is linked to AREFTIM.
!----------------------------------------------------------------------

      ! variables define a reference time
      INTEGER :: RefYear
      INTEGER :: RefMonth
      INTEGER :: RefDay
      INTEGER :: RefHour
      INTEGER :: RefMin
      INTEGER :: RefSec
      COMMON /RefTim/RefYear, RefMonth,RefDay,RefHour,RefMin,RefSec

! CREFTIM end
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


! declaration of local arrays
      real fieldVT(ncols,nrows)  ! field at validity time
      real fieldM6(ncols,nrows)  ! field at validity time minus 6 hours
      real fieldint(ncols,nrows) ! intermediate field for calculation

! declaration of local scalars
      real timediv        ! division scale for field (6x3600)
      integer IM6OffHr    ! secondary validity time offset for VT-6

! declaration of logicals
      logical l_preferred_VT   ! OUT test for preferred field at VT
      logical l_preferred_M6   ! OUT test for preferred field at VT-6
      logical l_previous_VT    ! OUT test for previous field at VT
      logical l_previous_M6    ! OUT test for previous field at VT-6

      CHARACTER(LEN=256) cmessage  ! error message



! declaration of externals
      external add_hours, read_one_field, read_climate_field,           &
     &  FieldSub,ScalarMult,check_header


!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_accum_flds'  ! subroutine name for error messages
      l_preferred_M6 = .false.
      l_preferred_VT = .false.
      l_previous_M6 = .false.
      l_previous_VT = .false.
      l_climate_field = .false.

      time = ValidityPeriod * 3600

! 1. calculate validity time minus 6 hours of NWP data required
      IM6OffHr = IVTOffHr - ValidityPeriod
! DEPENDS ON: add_hours
      call add_hours(                                                   &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       IM6OffHr)

!----------------------------------------------------------------------
! 2. Check headers for preferred and previous to see if they exist
!----------------------------------------------------------------------
      if ( LPreferred ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPreferred,                    &
     &                         LookupPreferred,                         &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                         l_preferred_M6)
      endif
      if ( LPrevious ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPrevious,                     &
     &                         LookupPrevious,                          &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                         l_previous_M6)
      endif

! 2.1 Calculate Validity Time and check if VT exists
! DEPENDS ON: add_hours
      call add_hours(                                                   &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       IVTOffHr)
      if ( LPreferred) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPreferred,                    &
     &                         LookupPreferred,                         &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                         l_preferred_VT)
      endif
      if ( LPrevious ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPrevious,                     &
     &                         LookupPrevious,                          &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                         l_previous_VT)
      endif

!----------------------------------------------------------------------
! 3. Read preferred VT&VT-6 if they exist else previous if they do
!----------------------------------------------------------------------
      if ( l_preferred_M6 .and. l_preferred_VT ) then
! DEPENDS ON: read_one_field
        call read_one_field (UnitPreferred, ITEM_CODE, StCode,          &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
     &       icode)

        if ( icode  <=  0) then
! 3.1 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//                                    &
     &     'NWP preferred field (VT) StCode ',                          &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 3.2 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//                                  &
     &     'NWP preferred field (VT) StCode ',                          &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
          l_preferred_VT = .false.
        end if
        icode = 0     ! reset icode

! 3.3 If preferred VT has been read, then read preferred VT-6
        if ( l_preferred_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &          IM6OffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPreferred, ITEM_CODE, StCode,        &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldM6,                                     &
     &       icode)

          if ( icode  <=  0) then
! 3.4 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP preferred field (VT-6) StCode ',                      &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else
! 3.5 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP preferred field (VT-6) StCode ',                      &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_preferred_M6 = .false.
          end if
          icode = 0     ! reset icode
        endif    ! l_preferred_VT
      endif    ! l_preferred_VT / l_preferred_M6

! 3.6 If either preferred VT or preferred M6 has not been read
!     read previous VT and VT-6
      if ( (.not. l_preferred_M6 .or. .not. l_preferred_VT) .and.       &
     &     (      l_previous_M6 .and. l_previous_VT       ) .and.       &
     &            LPrevious ) then
! DEPENDS ON: add_hours
        call add_hours(                                                 &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &        IVTOffHr)
! DEPENDS ON: read_one_field
        call read_one_field (UnitPrevious, ITEM_CODE, StCode,           &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
     &       icode)

        if ( icode  <=  0) then
! 3.7 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//                                    &
     &     'NWP previous field (VT) StCode ',                           &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 3.8 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//                                  &
     &     'NWP previous field (VT) StCode ',                           &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
          l_previous_VT = .false.
        end if
        icode = 0     ! reset icode
        if ( l_previous_VT) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &          IM6OffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPrevious, ITEM_CODE, StCode,         &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldM6,                                     &
     &       icode)

          if ( icode  <=  0) then
! 3.9 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP previous field (VT-6) StCode ',                       &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else
! 3.10 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP previous field (VT-6) StCode ',                       &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_previous_M6 = .false.
          end if
          icode = 0     ! reset icode
          endif   ! l_previous_VT
      endif

! Special case for first validity time in file: if accumulations to
! T+0 are not available, instead then first accumulation is from
! T+0, so only one field required.
!
      if ( IVTOffHr  ==  IValidOffHr(1) ) then
!
! 3.11 If no M6 fields have been read, but preferred VT is available
!      read preferred VT
        if ( (.not. l_preferred_M6 .and. .not. l_previous_M6) .and.     &
     &       l_preferred_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &        IVTOffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPreferred, ITEM_CODE, StCode,        &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
     &       icode)

          if ( icode  <=  0) then

! 3.12 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP preferred field (VT) StCode ',                        &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else

! 3.13 else write warning message and reset icode and validity time
!      for read of previous file
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP preferred field (VT) StCode ',                        &
     &      StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_preferred_VT = .false.
          end if
          icode = 0     ! reset icode

        endif

! 3.14 If no M6 fields have been read, and preferred VT is not
!      available, read previous VT if available
        if ( (.not. l_preferred_M6 .and. .not. l_previous_M6 .and.      &
     &        .not. l_preferred_VT ) .and. l_previous_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &        IVTOffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPrevious, ITEM_CODE, StCode,         &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
     &       icode)

          if ( icode  <=  0) then

! 3.15 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP previous field (VT) StCode ',                         &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else

! 3.16 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP previous field (VT) StCode ',                         &
     &      StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_previous_VT = .false.
          end if
          icode = 0     ! reset icode

        endif

      endif  ! First validity time in file
!----------------------------------------------------------------------
! 4. If there is a preferred field for VT-6 and for VT
!    or previous for VT-6 and VT then do accumulation
!----------------------------------------------------------------------
      if ( (l_preferred_VT .and. l_preferred_M6) .or.                   &
     &     (l_previous_VT .and. l_previous_M6) ) then
! DEPENDS ON: fieldsub
         call FieldSub (ncols,nrows,rmdi,                               &
     &                   fieldVT,fieldM6,                               &
     &                   fieldint,                                      &
     &                   icode,cmessage)

! 4.1 Now divide the result by a scalar using ScalarMult
        timediv = 1.0 / time
! DEPENDS ON: scalarmult
        call ScalarMult (ncols,nrows,rmdi,timediv,                      &
     &                      fieldint,field,                             &
     &                      icode,cmessage)

! 4.2 Write times to integer header
! DEPENDS ON: add_hours
        call add_hours(                                                 &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       IVTOffHr)
! DEPENDS ON: amend_times
        call amend_times (                                              &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                   Int_Head,Len_IntHd )
        goto 9999

      elseif ( l_preferred_VT .OR. l_previous_VT ) then

        timediv = 1.0 / time
! DEPENDS ON: scalarmult
        call ScalarMult (ncols,nrows,rmdi,timediv,                      &
     &                      fieldVT,field,                              &
     &                      icode,cmessage)

! DEPENDS ON: add_hours
        call add_hours(                                                 &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &       IVTOffHr)

! DEPENDS ON: amend_times
        call amend_times (                                              &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                   Int_Head,Len_IntHd )

        goto 9999

      endif   ! test for both fields

!----------------------------------------------------------------------
! 5. Otherwise extract field from climate file if available
!----------------------------------------------------------------------
      if (LClimate) then
! DEPENDS ON: read_climate_field
        call read_climate_field(StCode, IVTOffHr,                       &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
     &           icode)

         if ( icode  <=  0) then
            write(UnStd,*)CStd//CSub//'5. climate field extracted  ',   &
     &       ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            go to 9999
         else

            write(UnWarn,*)CWarn//CSub//                                &
     &       '5. failed to retrieve climate field ',                    &
     &       ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            icode = 0
         end if   ! icode
      end if !  LClimate

!----------------------------------------------------------------------
! 6. If no data has been successfully extracted return an error code
!----------------------------------------------------------------------
      icode = 5
      write(UnErr,*)CErr//CSub//'6. failed to extract any data',        &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr

9999  continue
      return
      END SUBROUTINE read_accum_flds
!----------------------------------------------------------------------
