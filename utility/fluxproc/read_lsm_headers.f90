! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!----------------------------------------------------------------------
! contains routines: read_lsm_headers
!
! Purpose: Flux processing routine.
!          Opens and reads lookup tables for land sea masks used by
!          FOAM_Flux_Process. Also sets LCyclic = T if atmosphere
!          grid has wrap-round points, and LCyclicO = T if
!          ocean grid has wrap-round points.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine read_lsm_headers (                                     &
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
                                   icode)

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
! dimensions of ocean and atmosphere fields
!----------------------------------------------------------------------
! comdeck: CFLDDIMS
! Purpose: declares dimensions of fields.
!          This deck is linked to AFLDDIMS.
!----------------------------------------------------------------------
! declarations:

!atmosphere grid
      integer ncols     !  number of columns
      integer nrowst    !  number of rows (tracer grid)
      integer nrowsu    !  number of rows (u grid)
      integer nrowsv    !  number of rows (v grid)
      integer nrowsuv   !  number of rows on grid with coincident velys
                        !  (at B grid locations both for B & C grids)
!ocean grid
      integer ncolsO    !  number of columns
      integer nrowstO   !  number of rows (tracer grid)
      integer nrowsuO   !  number of rows (velocity grid)
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

! no local arrays

! declaration of local scalars
      integer Len2_Lookup_lsm     ! max 2nd dimension for lsms
      integer Len2_Lookup_Actual  ! actual 2nd dimension for lsms
      integer IROW_NUMBER
      CHARACTER(LEN=80) cmessage

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_lsm_headers'! subroutine name for error messages
      Len2_Lookup_lsm = 1      ! all lsm ancillary files contain 1 field

! 0.1 Read StashMaster files
      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
     &  ICODE,CMESSAGE)

! 1. read atmosphere tracer land / sea mask fixed header and lookup
!    table from an an ancillary file
! DEPENDS ON: read_one_header
      call read_one_header(UnitNWPlsmt, icode,                          &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmt,                     &
     &               Lookuplsmt)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. unable to read NWP tracer land sea mask headers'
        go to 9999
      end if

! 1.1 extract the number of rows and columns from the lookup table
      ncols  = Lookuplsmt(LBNPT)
      nrowst = Lookuplsmt(LBROW)

! 2. Set dimensions for atmosphere velocity grids.
!    Calculations differ for B and C grids.
!    Note that Lookuplsmu and set_lookups_u are no longer used.
      if ( l_B_grid) then
        nrowsu = nrowst - 1
        nrowsv = nrowsu
        nrowsuv = nrowst - 1
      else
        nrowsu = nrowst
        nrowsv = nrowsu - 1
        nrowsuv = nrowst - 1
      end if

! 3. Set LCyclic (T if atmosphere grid has wrap points)
!    if fixhd(4)
      if ( MOD ( FixHdlsmt (4) , 100 )  /=  3 ) then
        LCyclic = .True.
      else
        LCyclic = .False.
      end if

! 4. read ocean tracer land / sea mask lookup table
! DEPENDS ON: read_one_header
      call read_one_header(UnitFOAMlsmt, icode,                         &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmtO,                    &
     &               LookuplsmtO)


      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4. unable to read ocean tracer land sea mask '
        go to 9999
      end if

! 4.1 extract the number of rows and columns from the lookup table
      ncolsO  = LookuplsmtO(LBNPT)
      nrowstO = LookuplsmtO(LBROW)

! 5. read ocean velocity land / sea mask lookup table
! DEPENDS ON: read_one_header
      call read_one_header(UnitFOAMlsmu, icode,                         &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmuO,                    &
     &               LookuplsmuO)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &    ' step 5. unable to read ocean velocity land sea mask '
        go to 9999
      end if

! 5.1 extract the number of rows and columns from the lookup table
      nrowsuO = LookuplsmuO(LBROW)

! 6. Set LCyclicO (T if ocean grid has wrap points)
!    if fixhd(4)
      if ( MOD ( FixHdlsmtO (4) , 100 )  /=  3 ) then
        LCyclicO = .True.
      else
        LCyclicO = .False.
      end if

9999  continue
      return
      END SUBROUTINE read_lsm_headers
!----------------------------------------------------------------------
