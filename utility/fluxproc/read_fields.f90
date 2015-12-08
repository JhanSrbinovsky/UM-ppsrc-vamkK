! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: read_fields
!
! Purpose: Flux processing routine.
!          Finds a field according to user's search criteria and
!          returns it and its lookup table by the argument list
!
! Uses StCode to read NWP and climate files
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_fields (StCode, IVTOffHr,                         &
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

      IMPLICIT NONE

! declaration of argument list

      integer StCode       ! IN   stash code of fields to search for

!       Reference date is used with IVTOffHr to define validity
!       time needed
      integer IVTOffHr     ! IN offset from validity time in hours

! debug control variable
      logical ldebug       ! IN T => output debugging info

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

! declaration of local logical

      logical l_climate_field ! Set to false initially

! no local arrays

! declaration of externals
      external add_hours, read_one_field, read_climate_field

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_fields'  ! subroutine name for error messages
      l_climate_field = .false.
! 1. calculate validity time of NWP data required

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

!----------------------------------------------------------------------
! 2. Extract field from preferred file if available
!----------------------------------------------------------------------
      if ( LPreferred ) then
! 2.1 try to read preferred field
! DEPENDS ON: read_one_field
        call read_one_field (UnitPreferred, ITEM_CODE, Stcode,          &
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
     &       ncols, nrows, field,                                       &
     &       icode)

        if ( icode  <=  0) then

! 2.2 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//'NWP preferred field stash code ',  &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
! 2.2.1  Write times to integer header
! DEPENDS ON: amend_times
          call amend_times (                                            &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                   Int_Head,Len_IntHd )
          go to 9999
        else

! 2.3 else write warning message and reset icode
         write(UnWarn,*)CWarn//CSub//'NWP preferred field stash code ', &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
        end if
        icode = 0     ! reset icode

      end if    ! LPreferred

!----------------------------------------------------------------------
! 3. Otherwise extract field from preferred file if available
!----------------------------------------------------------------------
      if ( LPrevious ) then

! 3.1 try to read previous field
! DEPENDS ON: read_one_field
       call read_one_field (UnitPrevious, ITEM_CODE, Stcode,            &
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
     &       ncols, nrows, field,                                       &
     &       icode)


        if ( icode  <=  0) then

! 3.2 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//'NWP previous field stash code ',   &
     &    Stcode, '; IVTOffHr = ', IVTOffHr, ' extracted'
! 3.2.1  Write times to integer header
! DEPENDS ON: amend_times
          call amend_times (                                            &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                   Int_Head,Len_IntHd )
          go to 9999
        else

! 3.3 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//'NWP previous field stash code ', &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
        end if
        icode = 0     ! reset icode

      end if    ! LPrevious

!----------------------------------------------------------------------
! 4. Otherwise extract field from climate file if available
!----------------------------------------------------------------------
      if ( LClimate ) then
! DEPENDS ON: read_climate_field
        call read_climate_field(StCode, IVTOffHr,                       &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
     &           icode)

        if ( icode  <=  0) then
          write(UnStd,*)CStd//CSub//'4. climate field extracted  ',     &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        else

          write(UnWarn,*)CWarn//CSub//                                  &
     &     '4. failed to retrieve climate field ',                      &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          icode = 0
        endif     ! icode
      endif     ! LClimate

!----------------------------------------------------------------------
! 5. If no data has been successfully extracted return an error code
!----------------------------------------------------------------------
      icode = 5
      write(UnErr,*)CErr//CSub//'5. failed to extract any data',        &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr

9999  continue
      return
      END SUBROUTINE read_fields
!----------------------------------------------------------------------
