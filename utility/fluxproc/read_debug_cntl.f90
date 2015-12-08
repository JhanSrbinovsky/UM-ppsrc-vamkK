! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_debug_cntl
!
! Purpose: Flux processing routine.
!          Reads files controlling diagnostic "debugging" output
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_debug_cntl ( icode )

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
! comdeck: CDEBUG
! Purpose: declares and stores information needed to output debugging
!          diagnostics at user selected points as fields are processed.
!----------------------------------------------------------------------
! common block:
      COMMON / CDebug /                                                 &
     &  NoDbgPts,IColDbg,JRowDbg,l_winds_dbg,l_heat_dbg,l_moisture_dbg, &
     &  l_sea_ice_dbg, l_references_dbg, l_pressure_dbg,l_windspd_dbg

      common / CCDebug / CValues

! declarations:
      integer MaxNoDbgPts   ! parameter: max. number of points to output
      parameter ( MaxNoDbgPts = 20)

! points to output
      integer NoDbgPts      ! actual number of points to output
      integer IColDbg(MaxNoDbgPts)   ! column of each point
      integer JRowDbg(MaxNoDbgPts)   ! row    of each point

! character array for output
      CHARACTER(LEN=11) CValues(MaxNoDbgPts) ! values to write out

! debug logical for each output file
      LOGICAL :: l_winds_dbg
      LOGICAL :: l_heat_dbg
      LOGICAL :: l_moisture_dbg
      LOGICAL :: l_sea_ice_dbg
      LOGICAL :: l_references_dbg
      LOGICAL :: l_pressure_dbg
      LOGICAL :: l_windspd_dbg

! CDEBUG end

! No local arrays

! declaration of local scalars
      integer i  ! loop index

! namelist declaration
      NAMELIST /NmLstDbg/                                               &
     &           NoDbgPts,                                              &
     &           IColDbg, JRowDbg,                                      &
     &           l_winds_dbg, l_heat_dbg, l_moisture_dbg,               &
     &           l_sea_ice_dbg, l_references_dbg, l_pressure_dbg,       &
     &           l_windspd_dbg

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_debug_cntl'  ! subroutine name for error messages

! 1. set default values for variables in NmLstDbg
      NoDbgPts = 0

      do i = 1, MaxNoDbgPts
        IColDbg(i) = 1
        JRowDbg(i) = 1
      end do

      l_winds_dbg = .false.
      l_heat_dbg = .false.
      l_moisture_dbg = .false.
      l_sea_ice_dbg = .false.
      l_references_dbg = .false.
      l_pressure_dbg = .false.
      l_windspd_dbg = .false.

! 2. read debug control namelist
      read (UnitDbg, NmLstDbg, iostat = icode)

      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 2. unable to read debug control namelist'
        icode = 8
        go to 9999
      end if

! 3. open file for debugging output
! DEPENDS ON: open_file
      call Open_file(OutUnitDbg,'Formatted  ','Unknown',icode)
      if (icode  /=  0) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 3. unable to open file for debugging output'
        icode = 9
        go to 9999
      end if

! 4. read and write out contents of namelist
      write(OutUnitDbg, NmLstDbg)

      if ( NoDbgPts  <=  0) then
        write(OutUnitDbg, *) ' no points to output '
      else
        write(OutUnitDbg, *) ' columns of output: '
        write(OutUnitDbg, '(11I11)' ) (IColDbg(i), i=1,NoDbgPts)
        write(OutUnitDbg, *) ' rows of output: '
        write(OutUnitDbg, '(11I11)' ) (JRowDbg(i), i=1,NoDbgPts)
      end if

9999  continue
      return
      END SUBROUTINE read_debug_cntl
!----------------------------------------------------------------------
