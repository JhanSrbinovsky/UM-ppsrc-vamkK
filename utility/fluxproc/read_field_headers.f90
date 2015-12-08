! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!
!----------------------------------------------------------------------
! contains routines: read_field_headers, read_field_sizes
!
! Purpose: Flux processing routine.
!          Reads fixed headers and lookup tables of  flux files input
!          to FOAM_Flux_Process (i.e. Preferred or Previous fluxes and
!          climate fluxes)
!          Also works out dimensions ncols, nrowsu, nrowst
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_field_headers (                                   &
! AFLDDIMA start
      ! argument list for dimensions of atmosphere fields
      ! linked to CFLDDIMA.
     &  ncols, nrowst, nrowsu, nrowsv, nrowsuv,                         &
! AFLDDIMA end
                                     icode )

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                              len_realhd, itgrid, iugrid,               &
                              max_num_fc_times, max_num_clim_times,     &
                              max_num_in_flux, len2_lookuppreferred,    &
                              len2_lookupprevious, len2_lookupclimate
      implicit none

! declaration of argument list
!----------------------------------------------------------------------
! comdeck: CFLDDIMA
! Purpose: declares dimensions of fields.
!          This deck is linked to AFLDDIMA and CFLDDIMS
!----------------------------------------------------------------------
! declarations:

! atmosphere grid
      integer ncols     !  number of columns
      integer nrowst    !  number of rows (tracer grid)
      integer nrowsu    !  number of rows (u grid)
      integer nrowsv    !  number of rows (v grid)
      integer nrowsuv   !  number of rows on grid with coincident velys
                        !  (at B grid locations both for B & C grids)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! comdeck: CSELCT
! Purpose: declares and stores information needed for flux selectiong
!----------------------------------------------------------------------
! common block:
      common / CSelect /                                                &
     &   l_B_grid,                                                      &
     &   l_winds_slt,   l_heat_slt,       l_moisture_slt,               &
     &   l_sea_ice_slt, l_references_slt, l_pressure_slt,               &
     &   l_windspd_slt

! debug logical for each selected flux
      logical l_B_grid,                                                 &
     &        l_winds_slt,   l_heat_slt,      l_moisture_slt,           &
     &        l_sea_ice_slt, l_references_slt,l_pressure_slt,           &
     &        l_windspd_slt
!----------------------------------------------------------------------
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
! declaration of logicals
      logical l_climate_field       ! T => Climate Field being used
      integer IROW_NUMBER
      CHARACTER(LEN=80) cmessage

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'read_field_headers' ! subroutine name for error messages


! 0.2 Read StashMaster files
      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)

!----------------------------------------------------------------------
! 1. Read and amend fixed header and lookups of preferred file
!----------------------------------------------------------------------


! 1.0 read headers

      LPreferred = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitPreferred, icode,                        &
     &   Len_FixHd, Len1_Lookup, Len2_LookupPreferred,                  &
     &   Len2_ActualPreferred, FixHdPreferred,                          &
     &   LookupPreferred)

      if (icode  /=  0) then
        LPreferred = .False.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1.0 unable to open and read headers of' //               &
     &  ' preferred flux file '
        icode = 0
      end if

! 1.1 amend headers

      if ( LPreferred ) then


! 1.1.2 amend headers

        l_climate_field = .false.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,  &
     &   l_climate_field, Len1_Lookup, Len2_LookupPreferred,            &
     &   Len2_ActualPreferred,                                          &
     &   LookupPreferred, LookFldNoPreferred, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 1.1 add_lookups failed for preferred file'
          go to 9999
        end if

      end if   ! LPreferred


!----------------------------------------------------------------------
! 2. Read and amend fixed header and lookups of previous file
!----------------------------------------------------------------------

! 2.0 read headers

      LPrevious = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitPrevious, icode,                         &
     &     Len_FixHd, Len1_Lookup, Len2_LookupPrevious,                 &
     &     Len2_ActualPrevious, FixHdPrevious,                          &
     &     LookupPrevious)

      if (icode  /=  0) then
        LPrevious = .False.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 2.0 failed to open and read headers of' //               &
     &  ' previous flux file '
        icode = 0
      end if

! 2.1 amend headers

      if ( LPrevious ) then


! 2.1.2 amend headers

        l_climate_field = .false.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,     &
     &   l_climate_field, Len1_Lookup, Len2_LookupPrevious,             &
     &   Len2_ActualPrevious,                                           &
     &   LookupPrevious, LookFldNoPrevious, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 2.1 add_lookups failed for Previous file'
          icode = icode + 2000
          go to 9999
        end if

      end if   ! LPrevious

!----------------------------------------------------------------------
! 3. Read and amend fixed header and lookups of climate file
!----------------------------------------------------------------------

! 3.0 read headers

      LClimate = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitClimate, icode,                          &
     &   Len_FixHd, Len1_Lookup, Len2_LookupClimate,                    &
     &   Len2_ActualClimate, FixHdClimate,                              &
     &   LookupClimate)

      if (icode  >   0) then
        LClimate = .false.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 3.0 failed to read headers of climate flux file '
        icode = 0
      end if

! 3.1 amend headers

      if ( LClimate ) then


! 3.1.2 amend headers

        l_climate_field = .true.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,        &
     &   l_climate_field, Len1_Lookup, Len2_LookupClimate,              &
     &   Len2_ActualClimate,                                            &
     &   LookupClimate, LookFldNoClimate, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 3.1 add_lookups failed for Climate file'
          icode = icode + 2500
          go to 9999
        end if

      end if   ! LClimate

!----------------------------------------------------------------------
! 4. If no file headers have been read exit with a fatal error
!----------------------------------------------------------------------
      if ( .not. ( LPreferred .or. LPrevious .or. LClimate) ) then
        icode = 16
         write(UnErr,*)CErr,CSub,                                       &
     &  ' step 4. failed to read headers of any flux file '
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_field_headers
!----------------------------------------------------------------------
