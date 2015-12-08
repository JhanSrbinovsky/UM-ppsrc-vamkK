! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top-level control routine for vegetation section
!
! Subroutine Interface:
      SUBROUTINE VEG_CTL(                                               &
     &                   row_length, rows, n_rows                       &
     &,                  global_row_length, global_rows                 &
     &,                  DIM_CS1, DIM_CS2                               &
     &,                  halo_i, halo_j, off_x, off_y, me               &
     &,                  n_proc, n_procx, n_procy                       &
     &,                  g_rows, g_row_length                           &
     &,                  at_extremity                                   &
     &,                  LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL           &
     &,                  A_STEP,ASTEPS_SINCE_TRIFFID                    &
     &,                  PHENOL_PERIOD,TRIFFID_PERIOD                   &
     &,                  L_PHENOL,L_TRIFFID,L_TRIF_EQ                   &
     &,                  ATIMESTEP,FRAC_DISTURB,SATCON                  &
     &,                  G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                &
     &,                  RESP_S_AC,RESP_W_AC                            &
     &,                  CS,FRAC,LAI,CLAY_FRAC,HT                       &
     &,                  CATCH_S,CATCH_T,INFIL_T,Z0_T                   &
      ,                  z0h_t                                          &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &                   STASHwork                                      &
     &                   )


      USE Submodel_Mod 
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE nstypes
      IMPLICIT NONE
!
! Description:  Calls interim control routine VEG_INTCTL
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: vegetation
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!


      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      INTEGER                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows           ! number of rows in a v field

      INTEGER                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        


      INTEGER                                                           &
     & LAND_PTS                                                         &
                                    ! IN Number of land points.
     &,NTILES                                                           &
                                    ! IN Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,A_STEP                                                           &
                                    ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                    ! INOUT Number of atmospheric
!                                   !       timesteps since last call
!                                   !       to TRIFFID.
     &,PHENOL_PERIOD                                                    &
                                    ! IN Phenology period (days).
     &,TRIFFID_PERIOD                                                   &
                                    ! IN TRIFFID period (days).
     &,DIM_CS1, DIM_CS2             ! IN soil carbon dimensions


      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! IN I=LAND_INDEX(L) => the Ith
!                                   !    point in P_FIELD is the Lth
!                                   !    land point.
      LOGICAL                                                           &
     & L_PHENOL                                                         &
                                    ! IN .T. for interactive leaf
!                                   !    phenology.
     &,L_TRIFFID                                                        &
                                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! IN Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_PTS)                                           &
                                    ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,SATCON(LAND_PTS)                                                 &
                                    ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_PTS,NPFT)                                    &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate including phenology.
     &,NPP_AC(LAND_PTS,NPFT)                                            &
                                    ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated wood respiration
!                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_PTS,DIM_CS1)                                      &
                                    ! INOUT Accumulated soil respiration
!                                   !       (kg C/m2).
     &,CS(LAND_PTS,DIM_CS1)                                             &
                               ! INOUT Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                    ! INOUT Fractions of surface types.
     &,LAI(LAND_PTS,NPFT)                                               &
                                    ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_PTS,NPFT)                                                &
                                    ! INOUT Height of plant functional
!                                   !       types (m).
     &,ALBSNC(LAND_PTS)                                                 &
                                    ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_PTS)                                                 &
                                    ! OUT Snow-free albedo.
     &,CATCH_S(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy snow capacity for tiles
!                                   !     (kg/m2).
     &,CATCH_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy capacity for tiles
!                                   !     (kg/m2).
     &,INFIL_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
      ,Z0_T(LAND_PTS,NTILES)       &! OUT Roughness length for tiles (m)
      ,Z0H_T(LAND_PTS,NTILES)       ! OUT Thermal roughness length for 
                                    !     tiles (m)

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
       REAL                                                             &
     &  STASHwork(*)    ! STASH workspace


! Local Variables

      REAL                                                              &
     & C_VEG(LAND_PTS,NPFT)                                             &
                                    ! LOCAL Total carbon content of
!                                   !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                    ! LOCAL Gridbox mean vegetation
!                                   !     carbon (kg C/m2).
     &,CS_TOT(DIM_CS2)                                                  &
                                    ! total soil carbon (kg C/m2).
     &,G_LEAF_PHEN(LAND_PTS,NPFT)                                       &
                                    ! LOCAL Mean leaf turnover rate over
!                                   !     phenology period (/360days).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                    ! LOCAL Carbon Litter
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_PTS)                                               &
                                    ! LOCAL Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! LOCAL Mean leaf turnover rate for
!                                   !     input to PHENOL (/360days).
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! LOCAL LAI of PFTs after phenology.
     &,G_LEAF_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! LOCAL Mean leaf turnover rate for
!                                   !       driving TRIFFID (/360days).
     &,NPP_DR_OUT(LAND_PTS,NPFT)                                        &
                                    ! LOCAL Mean NPP for driving TRIFFID
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! LOCAL Mean wood respiration for
!                                   !       driving TRIFFID
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_PTS,DIM_CS1+1)                                &
                                        ! Mean soil respiration for
                                        ! driving TRIFFID
                                        ! (kg C/m2/360days).
     &,CLAY_FRAC(ROW_LENGTH, ROWS)      ! IN  Clay fraction of soil

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('VEG_CTL',zhook_in,zhook_handle)

! DEPENDS ON: veg_ic
      CALL VEG_IC(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL                  &
     &,           A_STEP,ASTEPS_SINCE_TRIFFID                           &
     &,           PHENOL_PERIOD,TRIFFID_PERIOD                          &
     &,           L_PHENOL,L_TRIFFID,L_TRIF_EQ                          &
     &,           ATIMESTEP,FRAC_DISTURB,SATCON                         &
     &,           G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                       &
     &,           RESP_S_AC,RESP_W_AC                                   &
     &,           CS,FRAC,LAI,HT                                        &
      ,           CATCH_S,CATCH_T,INFIL_T,Z0_T,z0h_t                    &
     &,           C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN        &
     &,           LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT       &
     &,           RESP_S_DR_OUT                                         &
     &                )


! ----------------------------------------------------------------------
! Vegetation output diagnostics
! ----------------------------------------------------------------------

! Check that vegetation diagnostics requested this timestep
      IF (sf(0,19)) THEN

! DEPENDS ON: diagnostics_veg
        CALL diagnostics_veg(                                           &
     &                     row_length, rows, n_rows                     &
     &,                    global_row_length, global_rows               &
     &,                    DIM_CS1, DIM_CS2                             &
     &,                    halo_i, halo_j, off_x, off_y, me             &
     &,                    n_proc, n_procx, n_procy                     &
     &,                    g_rows, g_row_length                         &
     &,                    at_extremity                                 &
     &,                    land_pts                                     &
     &,                    land_index                                   &
     &,                    ntype,npft                                   &
     &,                    c_veg,cv,g_leaf_phen                         &
     &,                    lit_c,lit_c_mn,g_leaf_day                    &
     &,                    lai_phen,g_leaf_dr_out,npp_dr_out            &
     &,                    resp_w_dr_out,resp_s_dr_out,frac_disturb     &
     &,                    frac,lai,ht,cs                               &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &   STASHwork                                                      &
     &     )

      ENDIF

      IF (lhook) CALL dr_hook('VEG_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VEG_CTL
