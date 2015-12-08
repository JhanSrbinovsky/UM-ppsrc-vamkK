! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines:  open_flux_control_files, open_grid_control_files,
!                     open_control_files
!
! Purpose: Flux processing routine.
!          Opens all control and log files used by FOAM_Flux_Process
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

      subroutine open_control_files ( icode )

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn, cerr, &
                           cstd, csub
      USE chsunits_mod, ONLY : nunits
      implicit none

! declaration of argument list
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
! comdeck: CENVIRON
! Purpose: defines environment variables for units
!          which have to be opened by open_file
!----------------------------------------------------------------------
! declaration of parameters

! declarations of common blocks
      common / LEnviron /    LEnv
      common / CEnviron /    CEnv

! declarations of variables
      integer LEnv(NUnits)      ! lengths of environment variable names
      CHARACTER(LEN=15) CEnv(NUnits) ! names of environment variables
!----------------------------------------------------------------------

! No local arrays

! declaration of local scalars
      integer i  ! do loop index
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'open_cnt_files'  ! subroutine name for error messages

! 0.1 Before opening any files set values in CENVIRON
      do i = 1, NUNITS
        CEnv(i) = ' '
        LEnv(i) = 1
      end do
      CEnv(UnitNWPlsmt)   = 'FFLSMNWPT'
      LEnv(UnitNWPlsmt)   = 9
      CEnv(UnitFOAMlsmt)   = 'FFLSMFOAMT'
      LEnv(UnitFOAMlsmt)   = 15
      CEnv(UnitFOAMlsmu)   = 'FFLSMFOAMU'
      LEnv(UnitFOAMlsmu)   = 16


      CEnv(UnitPreferred) = 'FFPREFERRED'
      LEnv(UnitPreferred) = 11
      CEnv(UnitPrevious)  = 'FFPREVIOUS'
      LEnv(UnitPrevious)  = 10
      CEnv(UnitClimate)   = 'FFCLIMATE'
      LEnv(UnitClimate)   = 9


! 1. open log files

! 1.1 open error, warning and standard log files
! DEPENDS ON: open_file
      call open_file(UnErr, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(6,*)' ERROR opening error log file in FOAM flux '
        write(6,*)' processing. This job will have failed.  '
        write(6,*)CErr,CSub,' step 1.1 unable to open error log '
        icode = 1
        go to 9999
      end if

! 1.2 open warning log file
! DEPENDS ON: open_file
      call open_file(UnWarn, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.2 unable to open warning log '
        icode = 2
        go to 9999
      end if

! 1.3 open standard log file
! DEPENDS ON: open_file
      call open_file(UnStd, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.3 unable to open standard log '
        icode = 3
        go to 9999
      end if

! 2. open control files

! 2.1 open housekeeping file
! DEPENDS ON: open_file
      call open_file(UnitHK, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.1 unable to open housekeeping file'
        icode = 4
        go to 9999
      end if

! 2.2 open Validity time selection file
! DEPENDS ON: open_file
      call open_file(UnitVT, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &   ' step 2.2 unable to open Validity time selection file'
        icode = 5
        go to 9999
      end if

! 2.3 open debug control file
! DEPENDS ON: open_file
      call open_file(UnitDbg, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.3 unable to open debug control file'
        icode = 6
        go to 9999
      end if

! 2.4 open flux selection control file
! DEPENDS ON: open_file
      call open_file(UnitSlt, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.4 unable to open select control file'
        icode = 7
        go to 9999
      end if


9999  continue
      return
      END SUBROUTINE open_control_files
!----------------------------------------------------------------------

