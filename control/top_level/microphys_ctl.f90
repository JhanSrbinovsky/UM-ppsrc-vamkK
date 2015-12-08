! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine microphys_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

SUBROUTINE microphys_ctl (                                              &

! Parallel variables
  halo_i, halo_j, off_x, off_y, global_row_length,                      &
  at_extremity,                                                         &

! model dimensions.
  row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points,      &
  model_levels, wet_model_levels, bl_levels,                            &
  lspice_dim1,lspice_dim2,lspice_dim3,                                  &
  salt_dim1, salt_dim2, salt_dim3,                                      &
  cdnc_dim1, cdnc_dim2, cdnc_dim3,                                      &

! Model switches
  ltimer, l_rhcpt,                                                      &
  l_dust, l_sulpc_so2, l_sulpc_nh3, l_soot,                             &
  l_biomass, l_ocff, l_nitrate,                                         &
  l_cosp_lsp,                                                           &   

! Model parameters
  rhcrit,                                                               &

! Primary fields passed in
  t_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n,                    &
  cf_n, cfl_n, cff_n,                                                   &
  u_on_p, v_on_p,                                                       &
  snow_depth,                                                           &
  land_sea_mask, ice_fract,                                             &
  p_layer_centres, p_layer_boundaries,                                  &
  rho_r2,                                                               &
  aerosol,                                                              &
  ukca_cdnc,                                                            &
  dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,          &
  so2, nh3, so4_aitken, so4_accu, so4_diss,                             &
  aged_soot, cloud_soot, aged_bmass, cloud_bmass,                       &
  aged_ocff, cloud_ocff, nitr_acc, nitr_diss, biogenic,                 &
  sea_salt_film, sea_salt_jet, arcl,                                    &

! Other fields passed in
  ntml, cumulus,                                                        &
  fland, land_index,                                                    &

! Variables for stochastic physics random parameters
  m_ci,                                                                 &

! diagnostic info

! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
  stashwork4,                                                           &

! SCM diagnostics switches (dummy in full UM)
  nSCMDpkgs, L_SCMDiags,                                                &

! Increment fields passed in/out
  t_inc, q_inc, qcl_inc, qcf_inc,                                       &
  qcf2_inc, qrain_inc, qgraup_inc,                                      &
  cf_inc, cfl_inc, cff_inc,                                             &

! Fields required elsewhere
  ls_rain, ls_snow, micro_tends,                                        &
  cosp_gbx, n_drop_pot,                                                 &
! Field for Rh crit parametrization
  ext_p_layer_centres,                                                  &
  ext_tl,                                                               &
  ext_ql,                                                               &
  ext_qcf,                                                              &
  ext_ice_frac,                                                         &
  ext_land_frac,                                                        &

! error information
  error_code  )

! purpose: Interface to Atmospheric Physics parametrizations.
!          This version interfaces to cloud scheme, boundary layer,
!          hydrology, large scale rain, radiation, convection.
!          Called by atmos_physics1. 

! but Not:
!          gravity wave drag, and optionally energy
!          correction. Diagnostics also missing.

!          IN/OUT etc intents to be added later.

! code description: this code is written to umdp3 programming standards.

    ! Microphysics modules

  USE mphys_diags_mod,       ONLY: l_aggfr_diag, l_point_diag,           &
                                   l_psdep_diag, l_psaut_diag,           &
                                   l_psacw_diag, l_psacr_diag,           &
                                   l_psaci_diag, l_psmlt_diag,           &
                                   l_psmltevp_diag, l_praut_diag,        &
                                   l_pracw_diag, l_prevp_diag,           &
                                   l_pgaut_diag, l_pgacw_diag,           &
                                   l_pgacs_diag, l_pgmlt_diag,           &
                                   l_pifrw_diag, l_piprm_diag,           &
                                   l_pidep_diag, l_piacw_diag,           &
                                   l_piacr_diag, l_pimlt_diag,           &
                                   l_pimltevp_diag, l_pifall_diag,       &
                                   l_psfall_diag, l_prfall_diag,         &
                                   l_pgfall_diag, l_plset_diag,          &
                                   l_plevpset_diag, l_pifrr_diag,        &
                                   frac_agg, mphys_pts,                  &
                                   psdep, psaut,                         &
                                   psacw, psacr, psaci, psmlt,           &
                                   psmltevp, praut, pracw, prevp,        &
                                   pgaut, pgacw, pgacs, pgmlt, pifrw,    &
                                   piprm, pidep, piacw, piacr, pimlt,    &
                                   pimltevp, pifall, psfall, prfall,     &
                                   pgfall, plset, plevpset, pifrr

  USE mphys_bypass_mod,      ONLY: mp_dell

  USE cloud_inputs_mod,      ONLY: l_micro_eros, l_pc2

  USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

    ! General atmosphere modules
  USE arcl_mod,              ONLY: npd_arcl_compnts
  USE atmos_constants_mod,   ONLY: cp, r
  USE conversions_mod,       ONLY: zerodegc
  USE water_constants_mod,   ONLY: lc, lf
  USE dust_parameters_mod,   ONLY: ndiv, krain_dust, ksnow_dust,         &
                                   l_twobin_dust
  USE um_input_control_mod,  ONLY: l_mr_physics1
  USE timestep_mod,          ONLY: timestep

    ! Grid bounds module
  USE atm_fields_bounds_mod, ONLY: qdims, tdims, tdims_s,                &
                                   pdims
  USE level_heights_mod,     ONLY: r_theta_levels, r_rho_levels

    ! Dr Hook modules
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  !$ USE omp_lib

  USE Field_Types
  USE UM_ParParams
    ! COSP modules
  USE cosp_types_mod,        ONLY: cosp_gridbox

  USE diagnostics_lsrain_mod, ONLY: diagnostics_lsrain
  USE ls_ppn_mod, ONLY: ls_ppn
    ! Aerosols modules
  USE mass_calc_mod,      ONLY: mass_calc
  USE nh3dwash_mod,       ONLY: nh3dwash
  USE rainout_intctl_mod, ONLY: rainout_intctl
  USE sl3dwash_mod,       ONLY: sl3dwash
  USE Submodel_Mod 
  IMPLICIT NONE

  
! arguments with intent in. ie: input variables.

! Parallel setup variables
  INTEGER ::                                                            &
    halo_i,                                                             &
                   ! Size of halo in i direction.
    halo_j,                                                             &
                   ! Size of halo in j direction.
    off_x,                                                              &
                   ! Size of small halo in i
    off_y,                                                              &
                   ! Size of small halo in j.
    global_row_length
                   ! number of points on a row

  LOGICAL ::                                                            &
    at_extremity(4)  ! Indicates if this processor is at north,
!                          south, east or west of the processor grid

  LOGICAL ::                                                            &
    l_dust,                                                             &
               ! true for mineral dust calcs
    l_sulpc_so2,                                                        &
                    ! true for sulphur cycle calcs
    l_sulpc_nh3,                                                        &
                    ! true for sulphur/nh3 calcs
    l_soot,                                                             &
                    ! true for soot cycles calcs
    l_biomass,                                                          &
                    ! true for biomass aerosol calcs
    l_ocff,                                                             &
                    ! true for fossil-fuel organic carbon calcs
    l_nitrate
                    ! true for ammonium nitrate aerosol calcs

! Model dimensions
  INTEGER ::                                                            &
    row_length,                                                         &
    rows,                                                               &
    rhc_row_length,                                                     &
                       ! = row_length if L_RHCPT true, 1 otherwise.
    rhc_rows,                                                           &
                       ! = rows       if L_RHCPT true, 1 otherwise.
    n_rows,                                                             &
    land_points,                                                        &
                    ! IN No.of land points being processed, can be 0.
    model_levels,                                                       &
    wet_model_levels,                                                   &
    bl_levels,                                                          &
    lspice_dim1,                                                        &
                      ! Dimensions for 3D diagnostic arrays.
    lspice_dim2,                                                        &
                      ! These are set to 1 in order to save
    lspice_dim3,                                                        &
                      ! memory if the diagnostics are not used.
    salt_dim1,                                                          &

    salt_dim2,                                                          &
                      ! Dimensions for sea-salt aerosol arrays
    salt_dim3

  INTEGER, INTENT(IN) :: cdnc_dim1
  INTEGER, INTENT(IN) :: cdnc_dim2
  INTEGER, INTENT(IN) :: cdnc_dim3
                  ! Dimensions of UKCA cloud drop number concentration

  LOGICAL ::                                                            &
    ltimer,                                                             &
                 ! true then output some timing information
    l_rhcpt
                 ! Switch for 3D diagnosed RHcrit not 1D parameter

! Switch for COSP LS diagnostics
  LOGICAL ::  l_cosp_lsp

!  Additional arguments for SCM diagnostics which are dummy in full UM
  INTEGER ::                                                            &
    nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL ::                                                            &
    L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Model parameters
  REAL ::                                                               &
    rhcrit(1 : qdims%k_end)  ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
! Variable used in stochastic physics random parameters
  REAL ::                                                               &
   m_ci            ! Used to modify ice fall speed

! Primary fields passed in
! Note these are copies of prognostics from atmos_physics1
  REAL, INTENT (inout) ::                                               &
    t_n(      tdims%i_start : tdims%i_end,                              &
              tdims%j_start : tdims%j_end,                              &
                          1 : tdims%k_end ),                            &
    q_n(      qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qcl_n(    qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qcf_n(    qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qcf2_n(   qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qrain_n(  qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qgraup_n( qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cf_n(     qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cfl_n(    qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cff_n(    qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    u_on_p(   pdims%i_start : pdims%i_end,                              &
              pdims%j_start : pdims%j_end,                              &
              pdims%k_start : pdims%k_end ),                            &
    v_on_p(   pdims%i_start : pdims%i_end,                              &
              pdims%j_start : pdims%j_end,                              &
              pdims%k_start : pdims%k_end ),                            &
    snow_depth( tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end ),                          &
    ice_fract(  tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end ),                          &
    fland( land_points )

  INTEGER ::                                                            &
    land_index( land_points )

  LOGICAL ::                                                            &
    land_sea_mask( tdims%i_start : tdims%i_end,                         &
                   tdims%j_start : tdims%j_end )

  REAL ::                                                               &
    p_layer_boundaries( tdims%i_start : tdims%i_end,                    &
                        tdims%j_start : tdims%j_end,                    &
                                    0 : tdims%k_end ),                  &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.

    p_layer_centres(    tdims%i_start : tdims%i_end,                    &
                        tdims%j_start : tdims%j_end,                    &
                                    0 : tdims%k_end ),                  &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.

    rho_r2( tdims_s%i_start : tdims_s%i_end,                            &
            tdims_s%j_start : tdims_s%j_end,                            &
                          1 : tdims_s%k_end)
                               ! Density*earth_radius^2

  INTEGER ::                                                            &
    ntml( tdims%i_start : tdims%i_end,                                  &
          tdims%j_start : tdims%j_end )
                                           ! Height of diagnosed BL top

  LOGICAL ::                                                            &
    cumulus(  tdims%i_start : tdims%i_end,                              &
              tdims%j_start : tdims%j_end )
                                     ! Logical indicator for convection

  REAL, INTENT(INOUT) ::                                                &
                                     ! tracer variables
    aerosol   ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &

    dust_div1(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div1
    dust_div2(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div2
    dust_div3(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div3
    dust_div4(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div4
    dust_div5(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div5
    dust_div6(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
                           !dust mmr in div6
    so2       ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    nh3       ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    so4_aitken( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    so4_accu  ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    so4_diss  ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    aged_soot ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    cloud_soot( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    aged_bmass( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    cloud_bmass(tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    aged_ocff(  tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    cloud_ocff( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    nitr_acc  ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end ),                      &
    nitr_diss ( tdims_s%i_start : tdims_s%i_end,                        &
                tdims_s%j_start : tdims_s%j_end,                        &
                tdims_s%k_start : tdims_s%k_end )

    REAL, INTENT(IN) ::                                                 &
    arcl(       tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end,                            &
                tdims%k_start : tdims%k_end,                            &
                npd_arcl_compnts                ),                      &

    biogenic(   tdims%i_start   : tdims%i_end,                          &
                tdims%j_start   : tdims%j_end,                          &
                              1 : tdims%k_end   ),                      &
 
    sea_salt_film( salt_dim1,                                           &
                   salt_dim2,                                           &
                   salt_dim3 ),                                         &
!          Film-mode sea-salt aerosol number concentration
    sea_salt_jet(  salt_dim1,                                           &
                   salt_dim2,                                           &
                   salt_dim3 )
!          Jet-mode sea-salt aerosol number concentration

  REAL, INTENT(IN) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3) 
!                     CDNC from UKCA

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

! Diagnostics info
  REAL ::                                                               &
   stashwork4(*)     ! STASH workspace


! arguments with intent in/out. ie: input variables changed on output.
  REAL, INTENT (inout) ::                                               &
    t_inc(      tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end,                            &
                            1 : tdims%k_end ),                          &
    q_inc(      qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    qcl_inc(    qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    qcf_inc(    qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    qcf2_inc(   qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    qrain_inc(  qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    qgraup_inc( qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    cf_inc(     qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    cfl_inc(    qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end ),                          &
    cff_inc(    qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end )

  INTEGER ::                                                            &
    error_code

! arguments with intent out. ie: output variables.
  REAL ::                                                               &
    ls_rain( qdims%i_start : qdims%i_end,                               &
             qdims%j_start : qdims%j_end ),                             &
    ls_snow( qdims%i_start : qdims%i_end,                               &
             qdims%j_start : qdims%j_end ),                             &
    micro_tends( qdims%i_start : qdims%i_end,                           &
                 qdims%j_start : qdims%j_end, 2, bl_levels )
                        ! Tendencies from microphys within BL levels
                        ! (TL, K/s; QW, kg/kg/s)
                        ! The number 2 refers to the two categories
                        ! required by b.layer. 1 is TL and 2 is QW

! Gridbox information. Input for COSP
  TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx

! local variables.

  CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'microphys_ctl'

! loop counters
  INTEGER ::                                                            &
    i, j, k,                                                            &
    idiv, jj !loop counter for dust divisions
 
! open mp block size 
  INTEGER :: omp_block 

! Diagnostic switches

! Local variables
  REAL ::                                                               &
   land_frac(tdims%i_start : tdims%i_end,                               &
             tdims%j_start : tdims%j_end)
  REAL ::                                                               &
    rhcpt(rhc_row_length, rhc_rows, 1 : qdims%k_end),                   &
! Critical relative humidity
    ls_rain3d(lspice_dim1,lspice_dim2,lspice_dim3),                     &
! rainfall rate out of each model level for diagonstic
    ls_snow3d(lspice_dim1,lspice_dim2,lspice_dim3),                     &
! snowfall rate out of each model level for diagonstic
    rainfrac3d(lspice_dim1,lspice_dim2,lspice_dim3),                    &
! Rain fraction for diagnostic
    rnout_tracer(lspice_dim1,lspice_dim2),                              &
! Total tracer amount scavenged by rainout (kg m-2)
    lscav_tr(lspice_dim1,lspice_dim2),                                  &
! Total tracer amount scavenged by washout (kg m-2)
    rnout_soot(lspice_dim1,lspice_dim2),                                &
! Total soot amount scavenged by rainout (kg m-2)
    lscav_soot(lspice_dim1,lspice_dim2),                                &
! Total soot amount scavenged by washout (kg m-2)
    rnout_bmass(lspice_dim1,lspice_dim2),                               &
! Total biomass amount scavenged by rainout (kg m-2)
    lscav_bmass(lspice_dim1,lspice_dim2),                               &
! Total biomass amount scavenged by washout (kg m-2)
    rnout_ocff(lspice_dim1, lspice_dim2),                               &
! Total fossil-fuel organic carbon amount scavenged by rainout (kg m-2)
    lscav_ocff(lspice_dim1, lspice_dim2),                               &
! Total fossil-fuel organic carbon amount scavenged by washout (kg m-2)
    rnout_nitrate(lspice_dim1,lspice_dim2),                             &
! Total ammonium nitrate amount scavenged by rainout (kg[N] m-2)
    lscav_nitrate(lspice_dim1,lspice_dim2)
! Total ammonium nitrate amount scavenged by washout (kg[N] m-2)


  REAL ::                                                               &
    t_work(   tdims%i_start : tdims%i_end,                              &
              tdims%j_start : tdims%j_end,                              &
                          1 : tdims%k_end ),                            &
    q_work(   qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qcl_work( qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    qcf_work( qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cf_work(  qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cfl_work( qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    cff_work( qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end ),                            &
    n_drop_3d(qdims%i_start : qdims%i_end,                              &
              qdims%j_start : qdims%j_end,                              &
                          1 : qdims%k_end )
                    ! 3D calculated droplet number

  REAL, INTENT(OUT) ::                                                  &
    n_drop_pot( qdims%i_start : qdims%i_end,                            &
                qdims%j_start : qdims%j_end,                            &
                            1 : qdims%k_end )
!   Potential droplet number (includes values where cloud not present)

  REAL, INTENT(IN) ::                                                   &
    ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,                &
                                           0:qdims%k_end),              &
    ext_tl(0:rhc_row_length+1, 0:rhc_rows+1, 1: qdims%k_end),           &
    ext_ql(0:rhc_row_length+1, 0:rhc_rows+1, 1: qdims%k_end),           &
    ext_qcf(0:rhc_row_length+1,0:rhc_rows+1, 1: qdims%k_end),           &
    ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1),                      &
    ext_land_frac(0:rhc_row_length+1,0:rhc_rows+1)

  REAL ::                                                               &
   dm( tdims%i_start : tdims%i_end,                                     &
       tdims%j_start : tdims%j_end,                                     &
                   1 : tdims%k_end  )   ! mass air p.u. area in lev

      ! Local work arrays for additional microphysics fields if in use
  REAL, DIMENSION (:,:,:), ALLOCATABLE ::                               &
    qcf2_work, qrain_work, qgraup_work
  REAL ::                                                               &
    tinc_np1(   tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end ),                          &
    qinc_np1(   tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end ),                          &
    qclinc_np1( tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end ),                          &
    qcfinc_np1( tdims%i_start : tdims%i_end,                            &
                tdims%j_start : tdims%j_end )

! Local work array for increment diagnostics
  REAL ::                                                               &
    work_3d( tdims%i_start : tdims%i_end,                               &
             tdims%j_start : tdims%j_end,                               &
                         1 : tdims%k_end )

! Scavenged tracers (in column) for diagnostics
  REAL ::                                                               &
    lscav_so2( tdims%i_start : tdims%i_end,                             &
               tdims%j_start : tdims%j_end ),                           &
    lscav_nh3( tdims%i_start : tdims%i_end,                             &
               tdims%j_start : tdims%j_end ),                           &
    lscav_so4ait( tdims%i_start : tdims%i_end,                          &
                  tdims%j_start : tdims%j_end ),                        &
    lscav_so4acc( tdims%i_start : tdims%i_end,                          &
                  tdims%j_start : tdims%j_end ),                        &
    lscav_so4dis( tdims%i_start : tdims%i_end,                          &
                  tdims%j_start : tdims%j_end ),                        &
    lscav_dust_all( tdims%i_start : tdims%i_end,                        &
                    tdims%j_start : tdims%j_end, ndiv ),                &
    dust_all(       tdims%i_start : tdims%i_end,                        &
                    tdims%j_start : tdims%j_end,                        &
                                1 : tdims%k_end, ndiv )

  REAL ::                                                               &
  rainrate,                                                             &
                 !rate corrected for possible -ive precip
  snowrate,                                                             &
                 !rate corrected for possible -ive precip
  rate,                                                                 &

  delta_dust

! Microphysical process rate diagnostics
! Note: These arrays will only increase memory usage and are
!       only referenced if the particular diagnostic is active
!       (i.e. if the associated logical is .true., set from SF
!       STASH flag.

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MICROPHYS_CTL',zhook_in,zhook_handle)

  omp_block = (qdims%j_end - qdims%j_start) + 1

! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------
! Call timer for cloud code
! DEPENDS ON: timer
  IF (ltimer) CALL timer ('AP1M LS Cloud',5)

! Store values in work arrays

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SCHEDULE(STATIC)        &
!$OMP& SHARED(qdims, l_pc2, t_work, q_work, t_n, q_n, cf_n,             &
!$OMP& qcl_work, qcf_work, cf_work, cfl_work, cff_work, qcl_n, qcf_n,   &
!$OMP& cfl_n, cff_n)
  DO k = 1, qdims%k_end

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        IF (l_pc2) THEN
! No diagnostic cloud scheme is required.
!Store temperature and
! vapour content in T_work and q_work

          t_work(i,j,k)   = t_n(i,j,k)
          q_work(i,j,k)   = q_n(i,j,k)

        ELSE  ! l_pc2
! inputs to diagnostic cloud scheme are
! Tl_star, qT_star
! outputs from diag. cloud scheme
! are cloud_fraction, T, q, qcl, qcf
! output T is held in theta_star
! Calculate Tl, store in T_work
! Calculate qT, store in q_work

          t_work(i,j,k) = t_n(i,j,k) -                                  &
                          (lc * qcl_n(i,j,k) ) / cp
          q_work(i,j,k) = q_n(i,j,k) + qcl_n(i,j,k)

        END IF  ! l_pc2

        qcl_work(i,j,k) = qcl_n(i,j,k)
        qcf_work(i,j,k) = qcf_n(i,j,k)
        cf_work(i,j,k)  = cf_n(i,j,k)
        cfl_work(i,j,k) = cfl_n(i,j,k)
        cff_work(i,j,k) = cff_n(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  IF (l_micro_eros) THEN

! DEPENDS ON: pc2_turbulence_ctl
    CALL pc2_turbulence_ctl (                                           &
! Primary fields passed in, updated on exit
     T_work, q_work, qcl_work, cf_work, cfl_work, cff_work,             &
     p_layer_centres(1,1,1),                                            &
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     stashwork4,                                                        &
! SCM diagnostics switches (dummy in full UM)
     nSCMDpkgs, L_SCMDiags,                                             &
! Increment fields passed in, unchanged on exit
     t_inc, q_inc, qcl_inc, cf_inc, cfl_inc)

  END IF  ! l_micro_eros

  IF (l_mcr_qcf2) THEN

        ! Second cloud ice variable in use

    ALLOCATE ( qcf2_work( qdims%i_start : qdims%i_end,                  &
                          qdims%j_start : qdims%j_end,                  &
                                      1 : qdims%k_end ) )

        ! Add qcf and qcf2 to get total ice for cloud scheme
        ! qcf_work and qcf2_work are reset after call to ls_cld

    qcf_work(:,:,:) = qcf_n(:,:,:) + qcf2_n(:,:,:)

  ELSE

    ALLOCATE ( qcf2_work(1,1,1) )

  END IF

  IF (l_mcr_qrain) THEN

        ! Prognostic rain in use

    ALLOCATE ( qrain_work( qdims%i_start : qdims%i_end,                 &
                           qdims%j_start : qdims%j_end,                 &
                                       1 : qdims%k_end ) )

    qrain_work(:,:,:) = qrain_n(:,:,:)

  ELSE

    ALLOCATE ( qrain_work(1,1,1) )

  END IF

  IF (l_mcr_qgraup) THEN

        ! Prognostic graupel in use

    ALLOCATE ( qgraup_work(qdims%i_start : qdims%i_end,                 &
                           qdims%j_start : qdims%j_end,                 &
                                       1 : qdims%k_end ) )

    qgraup_work(:,:,:) = qgraup_n(:,:,:)

  ELSE

    ALLOCATE ( qgraup_work(1,1,1) )

  END IF

! Dry level T values (only required for diagnostics):

  DO k = qdims%k_end + 1, tdims%k_end

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        t_work(i,j,k) = t_n(i,j,k)

      END DO

    END DO

  END DO

  l_aggfr_diag    = sf(100,4)
  l_point_diag    = sf(101,4)
  l_pifrw_diag    = sf(240,4)
  l_piprm_diag    = sf(241,4)
  l_pidep_diag    = sf(243,4)
  l_psdep_diag    = sf(245,4)
  l_piacw_diag    = sf(247,4)
  l_psacw_diag    = sf(248,4)
  l_piacr_diag    = sf(249,4)
  l_psacr_diag    = sf(250,4)
  l_pimltevp_diag = sf(251,4)
  l_psmltevp_diag = sf(252,4)
  l_pimlt_diag    = sf(253,4)
  l_psmlt_diag    = sf(254,4)
  l_psaut_diag    = sf(255,4)
  l_psaci_diag    = sf(256,4)
  l_praut_diag    = sf(257,4)
  l_pracw_diag    = sf(258,4)
  l_prevp_diag    = sf(259,4)
  l_pgaut_diag    = sf(260,4)
  l_pgacw_diag    = sf(261,4)
  l_pgacs_diag    = sf(262,4)
  l_pgmlt_diag    = sf(263,4)
  l_pifall_diag   = sf(265,4)
  l_psfall_diag   = sf(266,4)
  l_prfall_diag   = sf(267,4)
  l_pgfall_diag   = sf(268,4)
  l_plset_diag    = sf(269,4)
  l_plevpset_diag = sf(270,4)
  l_pifrr_diag    = sf(271,4)

!==============================================================
! ALLOCATE arrays for required microphysics diagnostics
!==============================================================

      ! Aggregate Fraction
  IF (l_aggfr_diag) THEN
   ALLOCATE ( frac_agg( qdims%i_start : qdims%i_end,                    &
                        qdims%j_start : qdims%j_end,                    &
                                    1 : qdims%k_end ) )
    frac_agg(:,:,:) = 0.0
  ELSE
    ALLOCATE ( frac_agg(1,1,1) )
  END IF

      ! Microphysics Points Used
  IF (l_point_diag) THEN
   ALLOCATE ( mphys_pts( qdims%i_start : qdims%i_end,                   &
                         qdims%j_start : qdims%j_end,                   &
                                     1 : qdims%k_end ) )
    mphys_pts(:,:,:) = .FALSE.
  ELSE
    ALLOCATE ( mphys_pts(1,1,1) )
  END IF

      ! Homogeneous freezing nucl.
  IF (l_pifrw_diag) THEN
    ALLOCATE ( pifrw( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pifrw(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pifrw(1,1,1) )
  END IF

      ! Homoegenous freezing of rain
  IF (l_pifrr_diag) THEN
    ALLOCATE ( pifrr( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pifrr(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pifrr(1,1,1) )
  END IF

      ! Heterogeneous nucl.
  IF (l_piprm_diag) THEN
    ALLOCATE ( piprm( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    piprm(:,:,:) = 0.0
  ELSE
    ALLOCATE ( piprm(1,1,1) )
  END IF

      ! Deposition of vapour to ice
  IF (l_pidep_diag) THEN
    ALLOCATE ( pidep( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pidep(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pidep(1,1,1) )
  END IF

      ! Deposition of vapour to snow
  IF (l_psdep_diag) THEN
    ALLOCATE ( psdep( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    psdep(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psdep(1,1,1) )
  END IF

      ! Accretion of liq. by ice
  IF (l_piacw_diag) THEN
    ALLOCATE ( piacw( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    piacw(:,:,:) = 0.0
  ELSE
    ALLOCATE ( piacw(1,1,1) )
  END IF

      ! Accretion of liq. by snow
  IF (l_psacw_diag) THEN
    ALLOCATE ( psacw( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    psacw(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psacw(1,1,1) )
  END IF

      ! Collection of rain by ice
  IF (l_piacr_diag) THEN
    ALLOCATE ( piacr( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    piacr(:,:,:) = 0.0
  ELSE
    ALLOCATE ( piacr(1,1,1) )
  END IF

      ! Collection of rain by snow
  IF (l_psacr_diag) THEN
    ALLOCATE ( psacr( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    psacr(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psacr(1,1,1) )
  END IF

      ! Evaporation of melting ice
  IF (l_pimltevp_diag) THEN
    ALLOCATE ( pimltevp( qdims%i_start : qdims%i_end,                   &
                         qdims%j_start : qdims%j_end,                   &
                                     1 : qdims%k_end ) )
    pimltevp(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pimltevp(1,1,1) )
  END IF

      ! Evap. of melting aggregates
  IF (l_psmltevp_diag) THEN
    ALLOCATE ( psmltevp( qdims%i_start : qdims%i_end,                   &
                         qdims%j_start : qdims%j_end,                   &
                                     1 : qdims%k_end ) )
    psmltevp(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psmltevp(1,1,1) )
  END IF

      ! Melting of ice crystals
  IF (l_pimlt_diag) THEN
    ALLOCATE ( pimlt( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pimlt(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pimlt(1,1,1) )
  END IF

      ! Melting of snow aggregates
  IF (l_psmlt_diag) THEN
    ALLOCATE ( psmlt(qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ) )
    psmlt(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psmlt(1,1,1) )
  END IF

      ! Autoconversion of snow
  IF (l_psaut_diag) THEN
    ALLOCATE ( psaut( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    psaut(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psaut(1,1,1) )
  END IF

      ! Collection of ice crystals
  IF (l_psaci_diag) THEN
    ALLOCATE ( psaci( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    psaci(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psaci(1,1,1) )
  END IF

      ! Autoconversion of cloud
  IF (l_praut_diag) THEN
    ALLOCATE ( praut( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    praut(:,:,:) = 0.0
  ELSE
    ALLOCATE ( praut(1,1,1) )
  END IF

      ! Accretion of liq. by rain
  IF (l_pracw_diag) THEN
    ALLOCATE ( pracw( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pracw(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pracw(1,1,1) )
  END IF

      ! Evaporation of rain
  IF (l_prevp_diag) THEN
    ALLOCATE ( prevp( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    prevp(:,:,:) = 0.0
  ELSE
    ALLOCATE ( prevp(1,1,1) )
  END IF

      ! Autoconversion of graupel
  IF (l_pgaut_diag) THEN
    ALLOCATE ( pgaut( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pgaut(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pgaut(1,1,1) )
  END IF

      ! Accretion of liq. by graup
  IF (l_pgacw_diag) THEN
    ALLOCATE ( pgacw( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pgacw(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pgacw(1,1,1) )
  END IF

      ! Collection of snow by graup
  IF (l_pgacs_diag) THEN
    ALLOCATE ( pgacs( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pgacs(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pgacs(1,1,1) )
  END IF

      ! Melting of graupel
  IF (l_pgmlt_diag) THEN
    ALLOCATE ( pgmlt( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    pgmlt(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pgmlt(1,1,1) )
  END IF

      ! Sedimentation of ice crystals
  IF (l_pifall_diag) THEN
    ALLOCATE ( pifall( qdims%i_start : qdims%i_end,                     &
                       qdims%j_start : qdims%j_end,                     &
                                   1 : qdims%k_end ) )
    pifall(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pifall(1,1,1) )
  END IF

      ! Sedimentation of aggregates
  IF (l_psfall_diag) THEN
    ALLOCATE ( psfall( qdims%i_start : qdims%i_end,                     &
                       qdims%j_start : qdims%j_end,                     &
                                   1 : qdims%k_end ) )
    psfall(:,:,:) = 0.0
  ELSE
    ALLOCATE ( psfall(1,1,1) )
  END IF

      ! Sedimentation of rain
  IF (l_prfall_diag) THEN
    ALLOCATE ( prfall( qdims%i_start : qdims%i_end,                     &
                       qdims%j_start : qdims%j_end,                     &
                                   1 : qdims%k_end ) )
    prfall(:,:,:) = 0.0
  ELSE
    ALLOCATE ( prfall(1,1,1) )
  END IF

      ! Sedimentation of graupel
  IF (l_pgfall_diag) THEN
    ALLOCATE ( pgfall( qdims%i_start : qdims%i_end,                     &
                       qdims%j_start : qdims%j_end,                     &
                                   1 : qdims%k_end ) )
    pgfall(:,:,:) = 0.0
  ELSE
    ALLOCATE ( pgfall(1,1,1) )
  END IF

      ! Droplet settling of liquid water
  IF (l_plset_diag) THEN
    ALLOCATE ( plset( qdims%i_start : qdims%i_end,                      &
                      qdims%j_start : qdims%j_end,                      &
                                  1 : qdims%k_end ) )
    plset(:,:,:) = 0.0
  ELSE
    ALLOCATE ( plset(1,1,1) )
  END IF

      ! Evaporated settled droplets
  IF (l_plevpset_diag) THEN
    ALLOCATE ( plevpset( qdims%i_start : qdims%i_end,                   &
                         qdims%j_start : qdims%j_end,                   &
                                     1 : qdims%k_end ) )
    plevpset(:,:,:) = 0.0
  ELSE
    ALLOCATE ( plevpset(1,1,1) )
  END IF

! ----------------------------------------------------------------------
! Section CLD.1.a Calculate diagnostic RHcrit or read as namelist param.
! ----------------------------------------------------------------------

! Lrhcpt_if1:
  IF (l_rhcpt) THEN
!       RHCRIT is 3D diagnosed variable

! DEPENDS ON: ls_calc_rhcrit
    CALL ls_calc_rhcrit( ext_p_layer_centres,                           &
!            Array dimensions
       wet_model_levels, rhc_row_length, rhc_rows,                      &
       global_row_length,                                               &
!            Prognostic Fields
       ext_tl, ext_ql, ext_qcf, ext_land_frac, ext_ice_frac,            &
!            Logical controls
       l_mr_physics1,                                                   &
!            Output
       rhcpt, mp_dell )

  ELSE
!       RHCRIT is 1D Parameter read in from namelist
    DO k = 1, qdims%k_end
      rhcpt(1,1,k) = rhcrit(k)
    END DO
  END IF  ! Lrhcpt_if1

  IF (.NOT. l_pc2) THEN

! ----------------------------------------------------------------------
! Section CLD.1.b Calculate (2A) large-scale cloud fraction/ condensate.
! ----------------------------------------------------------------------

! DEPENDS ON: ls_cld
    CALL ls_cld( p_layer_centres(1,1,1), rhcpt,                         &
                 qdims%k_end, bl_levels,                                &
                 rhc_row_length, rhc_rows,                              &
                 ntml, cumulus, l_mr_physics1, t_work(1,1,1),           &
                 cf_work(1,1,1), q_work(1,1,1), qcf_work(1,1,1),        &
                 qcl_work(1,1,1), cfl_work(1,1,1), cff_work(1,1,1),     &
                 error_code )

  END IF  ! .not. l_pc2

! Call timer for cloud code
! DEPENDS ON: timer
  IF (ltimer) CALL timer ('AP1M LS Cloud',6)

  IF (error_code  ==  0 ) THEN

! ----------------------------------------------------------------------
! Section LSP.1 Call Large Scale precipitation scheme.
! ----------------------------------------------------------------------
! Call timer for LS rain code
! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1M LS Rain ',5)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

! initialise large scale rain and snow to zero.

!$OMP DO SCHEDULE(STATIC)
    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        ls_rain(i,j) = 0.0
        ls_snow(i,j) = 0.0

      END DO

    END DO
!$OMP END DO

! Set global land fraction

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

        land_frac(i,j) = 0.0

      END DO

    END DO
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, land_points

      j = (land_index(k)-1)/ ( tdims%i_end - tdims%i_start + 1 ) + 1

      i = land_index(k) - (j-1)*( tdims%i_end - tdims%i_start + 1 )

      land_frac(i,j) = fland(k)

    END DO
!$OMP END DO

!$OMP END PARALLEL  

    IF (l_mcr_qcf2) THEN  ! Second cloud ice variable in use
        ! Store ice variables in qcf_work and qcf2_work for LS_PPN
      qcf_work(:,:,:)  = qcf_n(:,:,:)
      qcf2_work(:,:,:) = qcf2_n(:,:,:)
    END IF

! theta_star holds current estimate of temperature at new time level.
! 3A
    CALL ls_ppn(                                                        &
                p_layer_boundaries, p_layer_centres(1,1,1),             &
                land_sea_mask,                                          &
   ! primary fields and
   ! cloud fractions
                cf_work, cfl_work, cff_work,                            &
                rhcpt,                                                  &
                lspice_dim1,lspice_dim2,lspice_dim3,                    &
                rho_r2, q_work, qcf_work, qcl_work, t_work,             &
                qcf2_work, qrain_work, qgraup_work,                     &
   ! Wind field and grid size 
   ! for lateral displacement 
   ! of falling ice by shear
                u_on_p, v_on_p,                                         &
   ! Aerosol variables
                sea_salt_film, sea_salt_jet,                            &
                salt_dim1, salt_dim2, salt_dim3,                        &
                ukca_cdnc,                                              &
                cdnc_dim1, cdnc_dim2, cdnc_dim3,                        & 
                biogenic,                                               &
                snow_depth, land_frac,                                  &
                so4_aitken( qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                so4_accu(   qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                so4_diss(   qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                aged_bmass( qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                cloud_bmass(qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                aged_ocff(  qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                cloud_ocff( qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                nitr_acc(   qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                nitr_diss(  qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                aerosol(    qdims%i_start : qdims%i_end,                &
                            qdims%j_start : qdims%j_end,                &
                                        1 : qdims%k_end ),              &
                arcl,                                                   &
    ! Other variables for mphys
                ls_rain, ls_snow,                                       &
                ls_rain3d, ls_snow3d, rainfrac3d,                       &
                n_drop_pot, n_drop_3d,                                  &
                rhc_row_length, rhc_rows,                               &
    ! Variables for stochastic physics random parameters
                m_ci,                                                   &
                error_code )

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)


! Calculate increment fields
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end

      DO j = qdims%j_start, qdims%j_end

        DO i = qdims%i_start, qdims%i_end

          t_inc(i,j,k) = t_work(i,j,k) - t_n(i,j,k) + t_inc(i,j,k)
          q_inc(i,j,k) = q_work(i,j,k) - q_n(i,j,k) + q_inc(i,j,k)
          qcl_inc(i,j,k) = qcl_work(i,j,k) - qcl_n(i,j,k)               &
                         + qcl_inc(i,j,k)
          qcf_inc(i,j,k) = qcf_work(i,j,k) - qcf_n(i,j,k)               &
                         + qcf_inc(i,j,k)
          cf_inc(i,j,k) = cf_work(i,j,k)   - cf_n(i,j,k)                &
                         + cf_inc(i,j,k)
          cfl_inc(i,j,k) = cfl_work(i,j,k)   - cfl_n(i,j,k)             &
                         + cfl_inc(i,j,k)
          cff_inc(i,j,k) = cff_work(i,j,k)   - cff_n(i,j,k)             &
                         + cff_inc(i,j,k)
        END DO

      END DO

    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels

      DO j = qdims%j_start, qdims%j_end

        DO i = qdims%i_start, qdims%i_end

              ! Calculate and save microphys increments to TL and QW,
              ! required in boundary layer.

          micro_tends(i,j,1,k) = t_inc(i,j,k)                           &
                                - (      lc*qcl_inc(i,j,k) )/cp         &
                                - ( (lc+lf)*qcf_inc(i,j,k) )/cp

          micro_tends(i,j,2,k) = q_inc(i,j,k) + qcl_inc(i,j,k)          &
                                              + qcf_inc(i,j,k)

          micro_tends(i,j,1,k) = micro_tends(i,j,1,k)/timestep
          micro_tends(i,j,2,k) = micro_tends(i,j,2,k)/timestep

        END DO

      END DO

    END DO
!$OMP END DO

!$OMP END PARALLEL 

    IF (l_mcr_qcf2)                                                     &
                       ! Second cloud ice variable in use
      qcf2_inc(:,:,:) = qcf2_work(:,:,:) - qcf2_n(:,:,:)                &
                      + qcf2_inc(:,:,:)

    IF (l_mcr_qrain)                                                    &
                        ! Prognostic rain in use
      qrain_inc(:,:,:) = qrain_work(:,:,:) - qrain_n(:,:,:)             &
                       + qrain_inc(:,:,:)

    IF (l_mcr_qgraup)                                                   &
                         ! Prognostic graupel in use
      qgraup_inc(:,:,:) = qgraup_work(:,:,:) - qgraup_n(:,:,:)          &
                        + qgraup_inc(:,:,:)

! Call timer for LS rain code
! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1M LS Rain ',6)
  END IF ! on error code zero

! ----------------------------------------------------------------------
! Section LSP.1.1 Tracer scavenging
! ----------------------------------------------------------------------

! below cloud scavenging of dust

! Scavenging of mineral dust
! (at present this is simplistic scheme based on 4.5 code)

!Initialise tracer scavenging to zero. Moved from ls_ppn, where
! it was being initialised, but not used.

  lscav_so2(:,:)    = 0.0
  lscav_nh3(:,:)    = 0.0
  lscav_so4ait(:,:) = 0.0
  lscav_so4acc(:,:) = 0.0
  lscav_so4dis(:,:) = 0.0

  IF (l_dust) THEN

    CALL mass_calc(                                                     &
   row_length, rows, model_levels, wet_model_levels,                    &
   r_rho_levels(1:row_length,1:rows,1:model_levels),                    &
   r_theta_levels(1:row_length,1:rows,0:model_levels),                  &
   timestep, rho_r2(1:row_length,1:rows,1:model_levels),                &
   q_n, qcl_n, qcf_n,                                                   &
   dm )


!   Put dust arrays together for simplicity

    DO idiv = 1, ndiv

      DO j = tdims%j_start, tdims%j_end

        DO i = tdims%i_start, tdims%i_end

          lscav_dust_all( i, j, idiv ) = 0.0

        END DO

      END DO

    END DO ! idiv

    DO k = 1, tdims%k_end

      DO j = tdims%j_start, tdims%j_end

        DO i = tdims%i_start, tdims%i_end

          dust_all(i,j,k,1)=dust_div1(i,j,k)
          dust_all(i,j,k,2)=dust_div2(i,j,k)
          IF (.NOT.l_twobin_dust) THEN
            dust_all(i,j,k,3)=dust_div3(i,j,k)
            dust_all(i,j,k,4)=dust_div4(i,j,k)
            dust_all(i,j,k,5)=dust_div5(i,j,k)
            dust_all(i,j,k,6)=dust_div6(i,j,k)
          END IF

        END DO

      END DO

    END DO ! k

!$  omp_block=CEILING(REAL(((qdims%j_end-qdims%j_start)+1))/                   &
!$  & omp_get_max_threads())

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, rainrate, snowrate,    &
!$OMP& rate, delta_dust, idiv, jj) SHARED(qdims, ls_rain3d, omp_block,   &
!$OMP& krain_Dust, dust_all, lscav_dust_all, ndiv, timestep, ksnow_dust, &
!$OMP& ls_snow3d, dm) 
    DO jj = qdims%j_start, qdims%j_end, omp_block
      DO k = qdims%k_end, 1, -1

        DO j = jj, MIN(jj+omp_block-1,qdims%j_end)

          DO i = qdims%i_start, qdims%i_end

!            deal with possible -ive precip

            rainrate = MAX( ls_rain3d(i,j,k), 0.0)
            snowrate = MAX( ls_snow3d(i,j,k), 0.0)

!              calc proportion of dust mixing ratio scavenged
!CDIR UNROLL=NDIV

            DO idiv = 1, ndiv

              rate = ( rainrate * krain_dust(idiv) +                      &
                  snowrate * ksnow_dust(idiv) ) *                         &
                  3600.0 * timestep

              delta_dust=dust_all(i,j,k,idiv)*(1.0-1.0/(1.0+ rate ))
!                calc mass of dust removed
              lscav_dust_all(i,j,idiv)=lscav_dust_all(i,j,idiv) +         &
              delta_dust * dm (i,j,k)
!                decrement mixing ratio
              dust_all(i,j,k,idiv)=dust_all(i,j,k,idiv) - delta_dust

            END DO ! idiv

          END DO ! i

        END DO ! j

      END DO ! k
    END DO ! jj
!$OMP END PARALLEL DO

!       put newly calculated values into main arrays

      DO k = 1, tdims%k_end

        DO j = tdims%j_start, tdims%j_end

          DO i = tdims%i_start, tdims%i_end

            dust_div1(i,j,k)=dust_all(i,j,k,1)
            dust_div2(i,j,k)=dust_all(i,j,k,2)
            IF (.NOT.l_twobin_dust) THEN
              dust_div3(i,j,k)=dust_all(i,j,k,3)
              dust_div4(i,j,k)=dust_all(i,j,k,4)
              dust_div5(i,j,k)=dust_all(i,j,k,5)
              dust_div6(i,j,k)=dust_all(i,j,k,6)
            END IF

          END DO

        END DO

      END DO

  END IF ! l_dust


  IF (l_sulpc_so2) THEN

    CALL rainout_intctl(row_length, rows,                               &
                 off_x, off_y,                                          &
                 halo_i, halo_j,                                        &
                 model_levels, wet_model_levels,                        &
                 rho_r2, q_n,                                           &
                 qcf_work, qcl_work,                                    &
                 qcf_n, qcl_n,                                          &
                 ls_rain3d, ls_snow3d,                                  &
                 timestep,                                              &
                 so4_diss(:,:,1:), so4_accu(:,:,1:),                    &
                 rnout_tracer)

    CALL sl3dwash(row_length, rows,                                     &
                  off_x, off_y,                                         &
                  halo_i, halo_j,                                       &
                  model_levels, wet_model_levels,                       &
                  timestep,                                             &
                  rho_r2,                                               &
                  q_n, qcl_n, qcf_n, so2(:,:,1:),                       &
                  ls_rain3d,                                            &
                  lscav_tr)

  END IF ! On test for sulphur cycle

  IF ( (l_sulpc_so2 .OR. l_nitrate) .AND. l_sulpc_nh3) THEN

    CALL nh3dwash(row_length, rows,                                     &
                  off_x, off_y,                                         &
                  halo_i, halo_j,                                       &
                  model_levels, wet_model_levels,                       &
                  timestep,                                             &
                  rho_r2,                                               &
                  q_n, qcl_n, qcf_n, nh3(:,:,1:),                       &
                  ls_rain3d,                                            &
                  lscav_nh3)
  END IF ! On test for (L_sulpc_SO2 .OR. L_nitrate) .AND. L_sulpc_NH3

   IF (l_soot) THEN  ! IF soot modelling is included

    CALL rainout_intctl(row_length, rows,                               &
                 off_x, off_y,                                          &
                 halo_i, halo_j,                                        &
                 model_levels, wet_model_levels,                        &
                 rho_r2, q_n,                                           &
                 qcf_work, qcl_work,                                    &
                 qcf_n, qcl_n,                                          &
                 ls_rain3d, ls_snow3d,                                  &
                 timestep,                                              &
                 cloud_soot(:,:,1:), aged_soot(:,:,1:),                 &
                 rnout_soot)


!       LS washout of soot was neglected at 4.5, although code was
!       allegedly written to perform the calculation. We can't use
!       the same treatment for soot as we use for the S cycle.
!       For now, we will again neglect this process, but note that
!       we may wish to include it in the future.

      DO j = 1, lspice_dim2
        DO i = 1, lspice_dim1
          lscav_soot(i,j)=0.0
        END DO
      END DO

  END IF  ! l_soot

  IF (l_biomass) THEN  ! IF biomass modelling is included

    CALL rainout_intctl(row_length, rows,                               &
                 off_x, off_y,                                          &
                 halo_i, halo_j,                                        &
                 model_levels, wet_model_levels,                        &
                 rho_r2, q_n,                                           &
                 qcf_work, qcl_work,                                    &
                 qcf_n, qcl_n,                                          &
                 ls_rain3d, ls_snow3d,                                  &
                 timestep,                                              &
                 cloud_bmass(:,:,1:), aged_bmass(:,:,1:),               &
                 rnout_bmass)


!       As with soot, LS washout of biomass smoke aerosol is
!       currently neglected. Again, though, we may wish to include
!       it in the future.

    DO j = 1, lspice_dim2
      DO i = 1, lspice_dim1
        lscav_bmass(i,j)=0.0
      END DO
    END DO

  END IF  ! l_biomass

  IF (l_ocff) THEN  ! IF fossil-fuel org carb modelling is included

    CALL rainout_intctl(row_length, rows,                               &
                 off_x, off_y,                                          &
                 halo_i, halo_j,                                        &
                 model_levels, wet_model_levels,                        &
                 rho_r2, q_n,                                           &
                 qcf_work, qcl_work,                                    &
                 qcf_n, qcl_n,                                          &
                 ls_rain3d, ls_snow3d,                                  &
                 timestep,                                              &
                 cloud_ocff(:,:,1:), aged_ocff(:,:,1:),                 &
                 rnout_ocff)


!       As with soot and biomass, LS washout of fossil-fuel organic
!       carbon aerosol is currently neglected. Again, though, we may
!       wish to include it in the future.

    DO j = 1, lspice_dim2
      DO i = 1, lspice_dim1
        lscav_ocff(i,j)=0.0
      END DO
    END DO

  END IF  ! l_ocff

  IF (l_nitrate) THEN  ! IF ammonium nitrate modelling is included

    CALL rainout_intctl(row_length, rows,                               &
                 off_x, off_y,                                          &
                 halo_i, halo_j,                                        &
                 model_levels, wet_model_levels,                        &
                 rho_r2, q_n,                                           &
                 qcf_work, qcl_work,                                    &
                 qcf_n, qcl_n,                                          &
                 ls_rain3d, ls_snow3d,                                  &
                 timestep,                                              &
                 nitr_diss(:,:,1:), nitr_acc(:,:,1:),                   &
                 rnout_nitrate)

!       As with soot, biomass, and OCFF, LS washout of ammonium
!       nitrate is currently neglected. Again, though, we may
!       wish to include it in the future.

    DO j = 1, lspice_dim2
      DO i = 1, lspice_dim1
        lscav_nitrate(i,j)=0.0
      END DO
    END DO

  END IF  ! l_nitrate

! ----------------------------------------------------------------------
! Section LSP.2 Output Diagnostics
! ----------------------------------------------------------------------

! Check that microphysics diagnostics requested this timestep
  IF (error_code  ==  0 .AND. sf(0,4)) THEN
    CALL diagnostics_lsrain(                                            &
                         lspice_dim1,lspice_dim2,lspice_dim3,           &
                         timestep,                                      &
                         at_extremity,                                  &
                         l_dust,                                        &
                         p_layer_centres,                               &
                         t_work, q_work, qcl_work, qcf_work,            &
                         qrain_work, qgraup_work, qcf2_work,            &
                         cf_work, cfl_work, cff_work,                   &
                         t_n, q_n, qcl_n, qcf_n,                        &
                         qrain_n, qgraup_n, qcf2_n,                     &
                         cf_n, cfl_n, cff_n,                            &
                         ls_rain, ls_snow,                              &
                         ls_rain3d,ls_snow3d,rainfrac3d,                &
                         rnout_tracer,lscav_dust_all,lscav_tr,          &
                         lscav_nh3,                                     &
                         rnout_soot, lscav_soot,                        &
                         rnout_bmass, lscav_bmass,                      &
                         rnout_ocff, lscav_ocff,                        &
                         rnout_nitrate, lscav_nitrate,                  &
                         psdep,psaut,psacw,psacr,                       &
                         psaci,psmlt,psmltevp,                          &
                         praut,pracw,prevp,                             &
                         pgaut,pgacw,pgacs,pgmlt,                       &
                         pifrw,pifrr,piprm,pidep,piacw,                 &
                         piacr,pimlt,pimltevp,                          &
                         pifall,psfall,prfall,pgfall,                   &
                         plset, plevpset, n_drop_pot, n_drop_3d,        &
                         frac_agg, mphys_pts,                           &

! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
   stashwork4                                                           &
         )

  END IF   ! on error_code and sf(0,4)
! Copy 3D rain and snow into COSP input structure
! The mass of snow is already contained in the LS ice water variable 
! (aggregates + crystals +[graupel]), so it does not need to be considered here.
  IF (L_cosp_lsp) THEN
    cosp_gbx%rain_ls(:,1:lspice_dim3) = reshape(ls_rain3d, &
              (/lspice_dim1*lspice_dim2,lspice_dim3/))
  END IF

  DEALLOCATE ( plevpset )     ! 33
  DEALLOCATE ( plset )        ! 32
  DEALLOCATE ( pgfall )       ! 31
  DEALLOCATE ( prfall )       ! 30
  DEALLOCATE ( psfall )       ! 29
  DEALLOCATE ( pifall )       ! 28
  DEALLOCATE ( pgmlt )        ! 27
  DEALLOCATE ( pgacs )        ! 26
  DEALLOCATE ( pgacw )        ! 25
  DEALLOCATE ( pgaut )        ! 24
  DEALLOCATE ( prevp )        ! 23
  DEALLOCATE ( pracw )        ! 22
  DEALLOCATE ( praut )        ! 21
  DEALLOCATE ( psaci )        ! 20
  DEALLOCATE ( psaut )        ! 19
  DEALLOCATE ( psmlt )        ! 18
  DEALLOCATE ( pimlt )        ! 17
  DEALLOCATE ( psmltevp )     ! 16
  DEALLOCATE ( pimltevp )     ! 15
  DEALLOCATE ( psacr )        ! 14
  DEALLOCATE ( piacr )        ! 13
  DEALLOCATE ( psacw )        ! 12
  DEALLOCATE ( piacw )        ! 11
  DEALLOCATE ( psdep )        ! 10
  DEALLOCATE ( pidep )        ! 9
  DEALLOCATE ( piprm )        ! 8
  DEALLOCATE ( pifrr )        ! 7
  DEALLOCATE ( pifrw )        ! 6
  DEALLOCATE ( mphys_pts )    ! 5
  DEALLOCATE ( frac_agg )     ! 4
  DEALLOCATE ( qgraup_work )  ! 3
  DEALLOCATE ( qrain_work )   ! 2
  DEALLOCATE ( qcf2_work )    ! 1

! end of routine microphys_ctl
  IF (lhook) CALL dr_hook('MICROPHYS_CTL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE microphys_ctl
