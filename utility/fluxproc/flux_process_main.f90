
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!
!----------------------------------------------------------------------
! contains progam: Flux_Process_Main
!
! Purpose: Flux processing routine.
!          Takes fields files of fluxes
!          (a) output from NWP  and
!          (b) NWP or coupled model or observed "climatologies"
!          combines them into fluxes as required by ocean models
!
!          The version of the code to build is determined by which
!          of FLUXPROC, FLXPLPR or FLXPLIN is defined.
!
!          If FLUXPROC is defined, the original executable is built.
!
!          If FLXPLPR (flux parallel processing), the executable to
!          process fluxes to the original NWP grid in parallel is
!          built; interpolation to ocean grid is performed by a
!          separate program.
!
!          If FLXPLIN (flux parallel interpolation), the executable to
!          interpolate fluxes to the ocean grid in parallel is
!          built.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      Program Flux_Process_Main
      USE IO  
      USE UM_ParVars
      USE UM_Config, ONLY : &
          appInit, &
          exe_fluxproc
      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn, &
                           cerr, cstd, csub
      USE ppxlook_mod, ONLY: ppxrecs
      USE Submodel_Mod
      IMPLICIT NONE

      INTEGER :: me_gc
      INTEGER :: nproc_gc

! declaration of parameters

! declaration of globals used

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
      CHARACTER(LEN=80) cmessage

!----------------------------------------------------------------------
! 0. Preliminaries
      CALL gc_init(' ',me_gc,nproc_gc)
      CALL appInit(exe_fluxproc)
      CALL ioInit()

      CSub = 'Flux_Process_Main'  ! subroutine name for error messages

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

 ! 1. Open all control and log files
! DEPENDS ON: open_control_files
      call open_control_files( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 1. Failed to open control and log files'
        go to 9999
      end if

! 2. Read all control files and open output flux files
! DEPENDS ON: read_control_files
      call read_control_files( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &' step 2. Failed to read control files or open output flux files'
        go to 9999
      end if

! 3. Open and read headers (fixed header & lookups) of flux fields
! DEPENDS ON: read_field_headers
      call read_field_headers(                                          &
! AFLDDIMA start
      ! argument list for dimensions of atmosphere fields
      ! linked to CFLDDIMA.
     &  ncols, nrowst, nrowsu, nrowsv, nrowsuv,                         &
! AFLDDIMA end
                              icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 3. Failed to read headers of flux fields'
        go to 9999
      end if

! 4. Open and read lookup tables of land-sea masks and find
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
     &    icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 4. Failed to read lookups of lsms'
        go to 9999
      end if


! 5. Do main processing at a lower level
! DEPENDS ON: flux_process
      call Flux_Process(                                                &
!----------------------------------------------------------------------
! comdeck: AFLDDIMS
! Purpose: argument list for dimensions of fields.
!          This deck is linked to FLDDIMS.
!----------------------------------------------------------------------
     & ncols, nrowst, nrowsu, nrowsv, nrowsuv, ncolsO, nrowstO, nrowsuO,&
!----------------------------------------------------------------------
     &     icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 4. Failed doing main flux processing'
        go to 9999
      end if

! 6. close files opened in steps 1. - 3.

! DEPENDS ON: close_files
      call close_files

9999  continue
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,Csub,                                        &
     &   'Flux Processing failed with error code = ',icode
        close ( UnErr )
      endif

      call ioShutdown()

      END PROGRAM Flux_Process_Main

!----------------------------------------------------------------------
