! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE diagnostics_riv_mod

IMPLICIT NONE

CONTAINS



      Subroutine diagnostics_riv(                                       &
     &                       row_length, rows                           &
     &,                       river_row_length, river_rows              &
     &,                      at_extremity                               &
     &,                      at_extremity_riv                           &
     &,                      RIVEROUT                                   &
     &,                      BOX_OUTFLOW, BOX_INFLOW                    &
! Add inland basin outflow to arguments
     &,                      TWATSTOR,INLANDOUT_RIV                     &

     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork                                                        &
     & )

! Description:
!   Calculates river-related diagnostics (held in STASH section 26).
!
! Method:
!   Each diagnostic is simply copied into the STASHwork array
!   to be passed on to STASH for output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    1    River water storage  (river grid)
!    2    gridbox outflow       (   "     )
!    3    gridbox runoff        (   "     )
!    4    coastal outflow       (ATMOS grid)
!         6    inland basin outflow       (river grid)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE Submodel_Mod
      IMPLICIT NONE

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
!     &                   ! south, east or west of the processor grid
     &,  at_extremity_riv(4) ! Indicates if this processor is at north,
!     &                   ! south, east or west of the processor grid
! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, river_row_length                                                &
                          ! river row length
     &, river_rows        ! river rows

      LOGICAL                                                           &
     & L_RIVERS           ! IN rivers switched on if .TRUE.

      REAL                                                              &
     & RIVEROUT(row_length, rows)                                       &
     &,BOX_OUTFLOW(river_row_length, river_rows)                        &
     &,BOX_INFLOW(river_row_length, river_rows)                         &
     &,TWATSTOR(river_row_length, river_rows)                           &

! Declare inland basin outflow variable
     &,INLANDOUT_RIV(river_row_length, river_rows)


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

! Local variables

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Integer, Parameter :: Sect = 26  !  Section No for RR diagnostics

      CHARACTER(LEN=80)  cmessage
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='diagnostics_riv')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)

! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------
! Section 1.  Initialisation.
! ------------------------------------------------------------------

      IF (lhook) CALL dr_hook('DIAGNOSTICS_RIV',zhook_in,zhook_handle)

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! --------------------------------------------------------------------
! River outflow (at seapoints)
! --------------------------------------------------------------------
! Item 26 004 riverout
      IF(sf(004,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(004,sect,im_index)),riverout,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,004,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 004)"
            goto 9999
         End if
      ENDIF

! ----------------------------------------------------------------------
! River outflow (at seapoints)
! -------------------------------------------------------------------
! Item 26 001 riverout
      IF(sf(001,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(001,sect,im_index)),TWATSTOR,       &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,001,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 001)"
            goto 9999
         End if
      ENDIF
!----------------------------------------------------------------------
! --------------------------------------------------------------------
! River outflow (at seapoints)
! --------------------------------------------------------------------
! Item 26 002 riverout
      IF(sf(002,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(002,sect,im_index)),BOX_OUTFLOW,    &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,002,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 002)"
            goto 9999
         End if
      ENDIF
!-------------------------------------------------------------------
! ------------------------------------------------------------------
! River outflow (at seapoints)
! ------------------------------------------------------------------
! Item 26 003 riverout
      IF(sf(003,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(003,sect,im_index)), BOX_INFLOW,    &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,003,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 003)"
            goto 9999
         End if
      ENDIF
!---------------------------------------------------------------------

! Output inland basin outflow on TRIP grid

! ------------------------------------------------------------------
! Inland basin outflow
! ------------------------------------------------------------------
! Item 26 006 inlandout_riv
      IF(sf(006,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(006,sect,im_index)),                &
     &    INLANDOUT_RIV,                                                &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,006,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 006)"
            goto 9999
         End if
      ENDIF
!---------------------------------------------------------------------

 9999 continue
      If(icode /= 0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      IF (lhook) CALL dr_hook('DIAGNOSTICS_RIV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_riv

END MODULE diagnostics_riv_mod
