! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_select_cntl
!
! Purpose: Flux processing routine.
!          Reads files controlling which fluxes to process
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_select_cntl ( icode )

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      implicit none

! declaration of argument list
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

! No local arrays

! No local scalars

! namelist declaration
      NAMELIST /NmLstSlt/                                               &
     &           l_B_grid,                                              &
     &           l_winds_slt, l_heat_slt, l_moisture_slt,               &
     &           l_sea_ice_slt, l_references_slt, l_pressure_slt,       &
     &           l_windspd_slt

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_select_cntl'  ! subroutine name for error messages

! 1. set default values for variables in NmLstSlt

      l_B_grid = .false.

      l_winds_slt = .false.
      l_heat_slt = .false.
      l_moisture_slt = .false.
      l_sea_ice_slt = .false.
      l_references_slt = .false.
      l_pressure_slt = .false.
      l_windspd_slt = .false.

! 2. read debug control namelist
      read (UnitSlt, NmLstSlt, iostat = icode)

      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 2. unable to read select control namelist'
        icode = 4008
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_select_cntl
!----------------------------------------------------------------------
