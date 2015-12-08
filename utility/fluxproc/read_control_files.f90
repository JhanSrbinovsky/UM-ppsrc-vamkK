! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_control_files
!
! Purpose: Flux processing routine.
!          Reads all control files used by FOAM_Flux_Process
!          Units added for pressure and windspeed (S. Spall)
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_control_files (icode)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,    &
                           cerr, cstd, csub
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

! declaration of  local arrays
      integer iunit_base(7)   ! main unit numbers for output files

! declaration of local scalars
      integer ivt       ! loop index over validity times
      integer iunit     ! loop index over unit numbers
      integer IAdd      ! value to add to basic unit number
      integer IUnitOpen ! unit number to open

! namelist declaration
      NAMELIST / NamFluxSelect /                                        &
     &  ValidityPeriod,                                                 &
     &  NoValidTimes, IValidOffHr, IOutUnitOff,                         &
     &  NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,   &
     &  NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,      &
     &  NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,         &
     &  output_land_value

! declaration of external subroutines and functions
      external read_debug_cntl, open_file

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_control_files' ! subroutine name for error messages

! 1. Read house keeping file to set reference date
      RefSec = 0
      RefMin = 0
! DEPENDS ON: readhk_fprdhk
      call readhk_fprdhk(UnitHK, RefHour, RefDay, RefMonth, RefYear,icode)

      if ( icode  /=  0 ) then
        write (UnErr,*)CErr,CSub,                                       &
     &   '1. Failed to read housekeeping file'
        goto 9999
      endif
      write(UnStd,*) CStd,CSub,'reference time from housekeeping file:' &
     & , ' RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec = ',      &
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec

! 2. Read debug control file and open debug ouput file
! DEPENDS ON: read_debug_cntl
      call read_debug_cntl ( icode )
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 2.  Failed to read debug control file'
        go to 9999
      end if

! 2.1 Read select control file
! DEPENDS ON: read_select_cntl
      call read_select_cntl ( icode )
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 2.1  Failed to read select control file'
        go to 9999
      end if

! 3. Read validity times control file

! 3.0 Set defaults
      NoAddTimesPreferred = 0
      NoAddTimesPrevious  = 0
      NoAddTimesClimate   = 0

      do ivt = 1, MaxTimes
        IValidOffHr(ivt) = 0
        IOutUnitOff(ivt) = 0
        ISrchOffHrPreferred(ivt) = 0
        INewOffHrPreferred(ivt) = 0
        ISrchOffHrPrevious(ivt) = 0
        INewOffHrPrevious(ivt) = 0
        ISrchOffHrClimate(ivt) = 0
        INewOffHrClimate(ivt) = 0
      end do

      output_land_value = rmdi

! 3.1 read namelist

      read (UnitVT, NamFluxSelect,  iostat = icode)
      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 3.1  Failed to read validity times control file'
        icode = 10
        go to 9999
      end if

      write(UnStd, NamFluxSelect)

! 4. set which units to open for output flux files

      do iunit = IUnOutLow, IUnOutHi
        LUnOutOpen(iunit) = .False.
      end do

      iunit_base(1) = UnitWindsOut
      iunit_base(2) = UnitHeatOut
      iunit_base(3) = UnitMoistureOut
      iunit_base(4) = UnitSeaIceOut
      iunit_base(5) = UnitReferencesOut
      iunit_base(6) = UnitPressureOut
      iunit_base(7) = UnitWindspdOut

      do ivt = 1, NoValidTimes
        IAdd = IOutUnitOff(ivt)
        do iunit = 1, 7
          IUnitOpen = iunit_base(iunit) + IAdd
          if ( IUnitOpen  <   IUnOutLow .or.                            &
     &        IUnitOpen  >   IUnOutHi ) then
            icode = 11
            write(UnErr,*)CErr,CSub,' step 4. Unit number chosen'       &
     &      ,' incorrectly; ivt,iunit =',ivt,iunit
            go to 9999
          else
            LUnOutOpen(IUnitOpen) = .True.
          end if
        end do ! iunit
      end do ! ivt

! 5. open output flux files
      do iunit = IUnOutLow, IUnOutHi
        if ( LUnOutOpen(iunit) ) then

! DEPENDS ON: open_file
          call open_file ( iunit, 'unformatted', 'unknown', icode )
          write(UnStd,*)CStd,CSub, ' step 5.  Opening file ', iunit

          if ( icode  >   0) then
            write(UnErr,*)CErr,CSub,                                    &
     &      ' step 5.  Failed to open output flux file ', iunit
            icode = 12
            go to 9999
          end if  ! icode

        end if ! LUnOutOpen(iunit)
      end do ! iunit

9999  continue
      return
      END SUBROUTINE read_control_files
!----------------------------------------------------------------------
