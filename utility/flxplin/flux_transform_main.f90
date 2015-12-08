! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! contains progam: Flux_Transform_Main
!
! Purpose: Flux processing routine.
!          Takes pp files of fluxes as input, interpolates them
!          from an atmosphere to an ocean grid and fills in
!          missing data values
!
!    Programming standard : UMDP 3
!
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      Program Flux_Transform_Main

      USE UM_ParVars
      USE UM_Config, ONLY : &
          appInit, &
          exe_flux_transform
      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      USE ppxlook_mod, ONLY: ppxrecs
      USE Submodel_Mod
      IMPLICIT NONE

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

! declaration of local scalars
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
      integer icode  ! error code ; > 0 => fatal error detected
      integer iunit  ! loop counter
      integer ifile  ! loop counter
      integer NoInFiles ! number of input files
      CHARACTER(LEN=80) cmessage
      CHARACTER(LEN=4)  C_NoInFiles  ! Char variable to read env var

      INTEGER :: me_gc
      INTEGER :: nproc_gc

!----------------------------------------------------------------------
! 0. Preliminaries
      CALL gc_init(' ',me_gc,nproc_gc)
      CALL appInit(exe_flux_transform)
      CSub = 'Flux_Transform_Main'  ! subroutine name for err messages

      icode = 0   ! initialise icode

! 0.1 Initialise N_INTERNAL_MODEL/INTERNAL_MODEL_INDEX
      N_INTERNAL_MODEL=2
      INTERNAL_MODEL_INDEX(1)=1    !  Atmos
      INTERNAL_MODEL_INDEX(2)=2    !  Ocean

! 0.2 Read STASHmaster files
      ppxRecs=1
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_A',ICODE,CMESSAGE)
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_O',ICODE,CMESSAGE)

! 0.3 Read NoInFiles

      CALL FORT_GET_ENV('NO_FLUX_FILES',13,C_NoInFiles,4,icode)
      IF (icode  /=  0) THEN
        WRITE(6,*) 'Warning : Environment variable NO_FLUX_FILES ',     &
     &             'has not been set.'
        WRITE(6,*) 'Setting NoInFiles to 6'
        NoInFiles=6
      ELSE
        READ(C_NoInFiles,'(I2)') NoInFiles
        write (6,*) ' '
        write (6,*) 'NO_FLUX_FILES is set to ',NoInFiles
      ENDIF
 ! 1. Open all control and log files
! DEPENDS ON: open_grid_control_files
      call open_grid_control_files( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 1. Failed to open control and log files'
        go to 9999
      end if

! 2. Read debug control file and open debug ouput file
! DEPENDS ON: read_debug_cntl
      call read_debug_cntl( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &' step 2. Failed to read debug control file'
        go to 9999
      end if

! 3. Open and read lookup tables of land-sea masks and find
!    field dimensions
! DEPENDS ON: read_lsm_headers
      call read_lsm_headers(                                            &
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
                            icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 4. Failed to read lookups of lsms'
        go to 9999
      end if

! 4. Open  input and output flux files

! 4.0 Set all logicals showing which files are open to false
      do iunit = IUnOutLow, IUnOutHi
        LUnOutOpen(iunit) = .False.
      end do

! 4.1 Open input flux file  ! should this status be 'old' or 'share' ?
      do ifile = 0, NoInFiles-1
        LUnOutOpen(IUnOutLow+10+ifile)= .True.
! DEPENDS ON: open_file
        call open_file (IUnOutLow+10+ifile,                             &
     &                  'unformatted', 'unknown', icode )
      enddo

! 4.2 Open output flux file
      LUnOutOpen(IUnOutLow+10+NoInFiles)= .True.
! DEPENDS ON: open_file
      call open_file (IUnOutLow+10+NoInFiles,                           &
     &                  'unformatted', 'unknown', icode )

! 5. Do main processing at a lower level
! DEPENDS ON: flux_transform_grid
      call Flux_Transform_Grid(                                         &
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
                               NoInFiles,icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 5. Failed doing main processing'
        go to 9999
      end if

! 6. close files opened in steps 1. - 4.

! DEPENDS ON: close_grid_files
      call close_grid_files

9999  continue
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,Csub,                                        &
     &   'Flux Processing failed with error code = ',icode
        close ( UnErr )
      endif

      END PROGRAM Flux_Transform_Main
!----------------------------------------------------------------------
