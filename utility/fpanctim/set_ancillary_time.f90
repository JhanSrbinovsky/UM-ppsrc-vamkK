! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains progam:   set_ancillary_time
!
! Purpose: Flux processing routine.
!          Writes namelists which determine the validity times in the
!          fixed headers of the FOAM flux ancillary files
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      PROGRAM set_ancillary_time

      USE UM_Config, ONLY : &
          appInit, &
          exe_set_ancillary_time
      USE ereport_mod, ONLY: ereport
      IMPLICIT NONE

      INTEGER :: me_gc
      INTEGER :: nproc_gc

! declaration of globals

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

! declaration of local scalars

      INTEGER fvhh,fvdd,fvmm,fvyy    ! first validity time
      INTEGER lvhh,lvdd,lvmm,lvyy    ! last validity time
      INTEGER ivhh,ivdd,ivmm,ivyy    ! interval between fields in file

      LOGICAL year360                ! true for 360-day calendar

      INTEGER first_vt_offset        ! offset in hours of first validity
                                     ! time from reference time

      INTEGER last_vt_offset         ! offset in hours of first validity
                                     ! time from reference time

      INTEGER imins, isecs           ! minutes and seconds discarded

      INTEGER icode                  ! error return code; > 0 => error


! declaration of namelists

      NAMELIST /offset_vt/ first_vt_offset,last_vt_offset

      NAMELIST /first_vt/ fvhh,fvdd,fvmm,fvyy
      NAMELIST /last_vt/  lvhh,lvdd,lvmm,lvyy
      NAMELIST /interval/ year360,ivhh,ivdd,ivmm,ivyy

! declaration of externals
      EXTERNAL open_file, add_hours
!----------------------------------------------------------------------

! 0. Preliminaries
      CALL gc_init(' ',me_gc,nproc_gc)
      CALL appInit(exe_set_ancillary_time)
      icode = 0

! 1. Open files for input and output

! 1.1 open housekeeping file
! DEPENDS ON: open_file
      CALL open_file(UnitHK, 'Formatted  ', 'Unknown', icode)
      IF ( icode  /=  0 ) THEN
        WRITE(6,*) ' ERROR: set_ancillary_time: 1.1 ',                  &
     &             ' unable to open housekeeping file'
        go to 9999
      END IF


!  open Validity time selection file
! DEPENDS ON: open_file
      CALL open_file(UnitVT, 'Formatted  ', 'Unknown', icode)
      IF ( icode  /=  0 ) THEN
        WRITE(6,*) ' ERROR: set_ancillary_time: 1.2 ',                  &
     &             ' unable to open Validity time selection file'
        go to 9999
      END IF

! open output file for ancillary file validity times namelist
! DEPENDS ON: open_file
      CALL open_file(UnitVTOut, 'Formatted  ', 'Unknown', icode)
      IF ( icode  /=  0 ) THEN
        WRITE(6,*) ' ERROR: set_ancillary_time: 1.3 ',                  &
     &             ' unable to open Validity times output file'
        go to 9999
      END IF


! 1. Read house keeping file to set reference date
      RefSec = 0
      RefMin = 0
! DEPENDS ON: readhk_fprdhk
      CALL readhk_fprdhk(UnitHK, RefHour, RefDay, RefMonth, RefYear,icode)

! Read offset_vt  namelist
      first_vt_offset =  0
      last_vt_offset =  0

      READ ( UnitVT, offset_vt )
      WRITE ( 6, offset_vt )

! Read interval namelist

      year360=.FALSE.
      ivhh=0
      ivdd=0
      ivmm=0
      ivyy=0

      READ ( UnitVT, interval )
      WRITE ( 6, interval )

! Set first validity time
! DEPENDS ON: add_hours
      CALL add_hours (                                                  &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
     &  fvyy, fvmm, fvdd, fvhh, imins, isecs,                           &
     &  first_vt_offset )

      WRITE (6, first_vt)

! Set last validity time
! DEPENDS ON: add_hours
      CALL add_hours (                                                  &
!----------------------------------------------------------------------
! comdeck: AREFTIM
! Purpose: argument list for variables storing a reference time.
!          This deck is linked to CREFTIM.
!----------------------------------------------------------------------
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec,              &
!----------------------------------------------------------------------
     &  lvyy, lvmm, lvdd, lvhh, imins, isecs,                           &
     &  last_vt_offset )

      WRITE (6, last_vt)

! Write output namelists

      WRITE ( UnitVTOut, first_vt )
      WRITE ( UnitVTOut, interval )
      WRITE ( UnitVTOut, last_vt )

! close files
      CLOSE (UnitHK)
      CLOSE (UnitVT)
      CLOSE (UnitVTOut)

9999  CONTINUE

      END PROGRAM set_ancillary_time
!----------------------------------------------------------------------
