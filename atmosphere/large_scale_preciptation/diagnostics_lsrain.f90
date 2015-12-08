! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: large scale precip
MODULE diagnostics_lsrain_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE diagnostics_lsrain(                                          &
                       lspice_dim1,lspice_dim2,lspice_dim3,             &
                       timestep,                                        &
                       at_extremity,                                    &
                       l_dust,                                          &
                       p_layer_centres,                                 &
                       t, q, qcl, qcf, qrain, qgraup, qcf2,             &
                       cf, cfl, cff,                                    &
                       t_n, q_n, qcl_n, qcf_n, qrain_n, qgraup_n,       &
                       qcf2_n, cf_n, cfl_n, cff_n,                      &
                       ls_rain, ls_snow,                                &
                       ls_rain3d,ls_snow3d,rainfrac3d,                  &
                       rnout_tracer,lscav_dust_all,lscav_tr,            &
                       lscav_nh3,                                       &
                       rnout_soot, lscav_soot,                          &
                       rnout_bmass, lscav_bmass,                        &
                       rnout_ocff, lscav_ocff,                          &
                       rnout_nitrate, lscav_nitrate,                    &
                       psdep,psaut,psacw,psacr,                         &
                       psaci,psmlt,psmltevp,                            &
                       praut,pracw,prevp,                               &
                       pgaut,pgacw,pgacs,pgmlt,                         &
                       pifrw,pifrr,piprm,pidep,piacw,                   &
                       piacr,pimlt,pimltevp,                            &
                       pifall,psfall,prfall,pgfall,                     &
                       plset,plevpset, n_drop_tpr, n_drop_3d,           &
                       frac_agg, mphys_pts,                             &

! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
 stashwork                                                              &
  )

  USE Submodel_Mod
  USE um_input_control_mod,  ONLY: l_mr_physics1

  USE conversions_mod,       ONLY: zerodegc

  USE ac_diagnostics_mod,    ONLY: lsrr, lssr, tinc_ppn

! Grid bounds module
  USE atm_fields_bounds_mod, ONLY: qdims, tdims

! Dr Hook modules
  USE dust_parameters_mod,   ONLY: ndiv
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim

! Purpose:
!          Calculates diagnostics generated from large scale
!          precipitation (UM section 4).

! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the large scale
! precipitation routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing, except where indicated.
! NOTE: Although moisture field diagnostics are available from this
! section (and correspond to equivalent variables calculated at UM4.5
! and earlier), the most appropriate place to extract moisture
! diagnostics is after section 9, which is the end of the moist
! processes calculations during the timestep.

!  Diagnostics currently available: (in order calculated)

! STASH item (all section 4 )
! ------------------------------------------------------------
! 100 aggregate fraction (model levels)
! 101 flag for where microphysics is performed        (model levels)
! 181 temperature increment across ls precip routines (model levels)
! 182 humidity    increment across ls precip routines (model levels)
! 183 qcl         increment across ls precip routines (model levels)
! 184 qcf         increment across ls precip routines (model levels)
! 189 qrain       increment across ls precip routines (model levels)
! 190 qgraup      increment across ls precip routines (model levels)
! 191 qcf2        increment across ls precip routines (model levels)
! 192 cf          increment across ls precip routines (model levels)
! 193 cfl         increment across ls precip routines (model levels)
! 194 cff         increment across ls precip routines (model levels)
! 201 large scale rain amount (kg/m2 per timestep)    (surface)
! 202 large scale snow amount (kg/m2 per timestep)    (surface)
! 203 large scale rainfall rate (kg/m2/s)             (surface)
! 204 large scale snowfall rate (kg/m2/s)             (surface)
!   4 temperature           after ls precip           (model levels)
! 205 cloud water (qcl)     after ls precip           (model levels)
! 206 cloud ice (qcf)       after ls precip           (model levels)
! 207 relative humidity     (percent)                 (model levels)
! 210 cloud droplet number / m3 (from autoconversion) (model levels)
! 211 cloud droplet number / m3 (from aerosols)       (model levels)
! 213 Large-scale rainout of dissolved ammonium nitrate (kg[N]/m2/s)
! 214 Large-scale washout of dissolved ammonium nitrate (kg[N]/m2/s)
!  10 specific humidity (q) after ls precip           (model levels)
! 220 Large scale rainout of soot (kg/m2/s)           (surface)
! 221 Large scale washout of soot (kg/m2/s)           (surface)
! 223 snowfall rate                                   (model levels)
! 224 supercooled liquid water content                (model levels)
! 225 supercooled rainfall rate                       (model levels)
! 227 rain fraction                                   (model levels)
! 228 Large scale rainout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 229 Large scale washout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 237 Large scale rainout of biomass smoke (kg/m2/s)  (surface)
! 238 Large scale washout of biomass smoke (kg/m2/s)  (surface)
! Microphysical process rate diagnostics all on wet_model_levels
! 240 Homogeneous nucleation rate (kg/kg/s)
! 241 Heterogeneous nucleation rate (kg/kg/s)
! 243 Ice deposition rate (kg/kg/s)
! 245 Snow deposition rate (kg/kg/s)
! 247 Ice collection rate of cloud liquid water (riming) (kg/kg/s)
! 248 Snow collection rate of cloud liquid water (riming) (kg/kg/s)
! 249 Ice collection rate of rain (capture) (kg/kg/s)
! 250 Snow collection rate of rain (capture) (kg/kg/s)
! 251 Evaporation rate of melting ice (kg/kg/s)
! 252 Evaporation rate of melting snow (kg/kg/s)
! 253 Melting rate for ice (kg/kg/s)
! 254 Melting rate for snow (kg/kg/s)
! 255 Snow aggregate autoconversion rate (kg/kg/s)
! 256 Snow collection rate of ice (capture) (kg/kg/s)
! 257 Rain autoconversion rate (kg/kg/s)
! 258 Rain collection rate of cloud liquid water (accretion) (kg/kg/s)
! 259 Evaporation rate of rain (kg/kg/s)
! 260 Graupel autoconversion rate  (kg/kg/s)
! 261 Graupel collection rate of cloud water (accretion) (kg/kg/s)
! 262 Graupel collection rate of snow (capture) (kg/kg/s)
! 263 Melting rate for graupel (kg/kg/s)
! 265 Ice crystal sedimentation rate (kg/kg/s)
! 266 Snow sedimentation rate (kg/kg/s)
! 267 Rain sedimentation rate (kg/kg/s)
! 268 Graupel sedimentation rate (kg/kg/s)
! 269 Droplet settling rate of liquid (kg/kg/s)
! 270 Evaporation rate for settling droplets (kg/kg/s)
! 271 Homogeneous freezing of rain (kg/kg/s)


!
! NOTE : The following PC2 diagnostics are part of section 4 but
!        are calculated elsewhere:
! In diagnostics_pc2checks:  141, 142, 143, 130, 131, 144, 132, 133
!                            152, 153, 136, 137, 154, 138, 139
! In pc2_turbulence_control: 281, 282, 283, 292, 293

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

  LOGICAL ::                                                            &
    at_extremity(4),                                                    &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
                         ! -----
    l_dust               ! Switch for mineral dust

! Parameters

  INTEGER, INTENT(IN) ::                                                &
    lspice_dim1,                                                        &
                            ! Dimensions for 3D diagnostic arrays.
    lspice_dim2,                                                        &
                            !  These are set to 1 in order to save
    lspice_dim3         !  memory if the diagnostics are not used.

  REAL, INTENT(IN) ::                                                   &
    timestep


! Primary Arrays used in all models
  REAL, INTENT(IN) ::                                                   &
    p_layer_centres( tdims%i_start : tdims%i_end,                       &
                     tdims%j_start : tdims%j_end,                       &
                                 0 : tdims%k_end ),                     &
    t(               tdims%i_start : tdims%i_end,                       &
                     tdims%j_start : tdims%j_end,                       &
                                 1 : tdims%k_end ),                     &
    q(               qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcl(             qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcf(             qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qrain(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qgraup(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcf2(            qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cf(              qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cfl(             qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cff(             qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &

! Time level n values for increment diagnostics
    t_n  (           tdims%i_start : tdims%i_end,                       &
                     tdims%j_start : tdims%j_end,                       &
                                 1 : tdims%k_end ),                     &
    q_n  (           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcl_n(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcf_n(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qrain_n  (       qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qgraup_n(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    qcf2_n(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cf_n  (          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cfl_n  (         qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    cff_n  (         qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    ls_rain(         qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end ),                     &
    ls_snow(         qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end ),                     &
    ls_snow3d(  lspice_dim1, lspice_dim2, lspice_dim3 ),                &
    rainfrac3d( lspice_dim1, lspice_dim2, lspice_dim3 ),                &
    n_drop_tpr(      qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    n_drop_3d(       qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    frac_agg(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )

LOGICAL, INTENT(IN) :: mphys_pts( qdims%i_start : qdims%i_end,          &
                                  qdims%j_start : qdims%j_end,          &
                                              1 : qdims%k_end )

! Microphysics Process Rate diagnostic arrays
  REAL, INTENT(IN) ::                                                   &
    psdep(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psaut(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psacw(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psacr(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psaci(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psmlt(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psmltevp(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )
  REAL, INTENT(INOUT) ::                                                &
    praut(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pracw(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    prevp(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )
  REAL, INTENT(INOUT) ::                                                &
    pgaut(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pgacw(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pgacs(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pgmlt(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )
  REAL, INTENT(INOUT) ::                                                &
    pifrw(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pifrr(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    piprm(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pidep(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    piacw(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    piacr(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pimlt(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pimltevp(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )
  REAL, INTENT(INOUT) ::                                                &
    pifall(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    psfall(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    prfall(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    pgfall(          qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )
  REAL, INTENT(INOUT) ::                                                &
    plset(           qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end ),                     &
    plevpset(        qdims%i_start : qdims%i_end,                       &
                     qdims%j_start : qdims%j_end,                       &
                                 1 : qdims%k_end )


! Used as input and workspace
  REAL, INTENT(INOUT) ::                                                &
    ls_rain3d(     lspice_dim1, lspice_dim2, lspice_dim3 ),             &
    rnout_tracer(  lspice_dim1, lspice_dim2 ),                          &
    lscav_dust_all(qdims%i_start : qdims%i_end,                         &
                   qdims%j_start : qdims%j_end, ndiv ),                 &
                                             !scavenged mineral dust
    lscav_tr(      lspice_dim1, lspice_dim2 ),                          &
    lscav_nh3(     lspice_dim1, lspice_dim2 ),                          &
    rnout_soot(    lspice_dim1, lspice_dim2 ),                          &
    lscav_soot(    lspice_dim1, lspice_dim2 ),                          &
    rnout_bmass(   lspice_dim1, lspice_dim2 ),                          &
    lscav_bmass(   lspice_dim1, lspice_dim2 ),                          &
    rnout_ocff(    lspice_dim1, lspice_dim2 ),                          &
    lscav_ocff(    lspice_dim1, lspice_dim2 ),                          &
    rnout_nitrate( lspice_dim1, lspice_dim2 ),                          &
    lscav_nitrate( lspice_dim1, lspice_dim2 )

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
  REAL, INTENT(INOUT) ::                                                &
   stashwork(*)     ! STASH workspace for section 4 (LS precip)

! Local variables
  INTEGER ::                                                            &
   i, j, k, ji,                                                         &
       icode                ! Return code  =0 Normal exit  >1 Error

  INTEGER :: sect,item    ! STASH section, item no.s
  PARAMETER (sect = 4) !  for microphysics - large scale rain

  REAL :: work_3d( tdims%i_start : tdims%i_end,                         &
                tdims%j_start : tdims%j_end,                            &
                            1 : tdims%k_end ) ! work space

  CHARACTER(LEN=80):: cmessage

  CHARACTER(LEN=*):: RoutineName 
  PARAMETER ( RoutineName='diagnostics_lsrain')

  INTEGER ::                                                            &
    im_index        ! internal model index

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('DIAGNOSTICS_LSRAIN',zhook_in,zhook_handle)
  icode = 0 ! Initialise error status
  im_index = internal_model_index(atmos_im)

!  Copy diagnostic information to STASHwork for STASH processing

!--------------------------------------------------------------------------
! Item 100: Aggregate Fraction
!--------------------------------------------------------------------------

  IF (icode <= 0 .AND. sf(100,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d( stashwork(si(100,4,im_index)), frac_agg,          &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,100,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,100,                                                &
         icode,cmessage )

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(frac_agg)"
    END IF

  END IF

!--------------------------------------------------------------------------
! Item 101: Flag for where microphysics is performed
!--------------------------------------------------------------------------

  IF (icode <= 0 .AND. sf(101,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d( stashwork(si(101,4,im_index)), mphys_pts,         &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,101,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,101,                                                &
         icode,cmessage )

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(mphys_pts)"
    END IF

  END IF

! increment diagnostics= modified - previous

  item = 181  ! temperature increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end

          work_3d(i,j,k) = t(i,j,k) - t_n(i,j,k)

        END DO ! i
      END DO   ! j
    END DO     ! k

! And set dry level increments to zero explicitly
    DO k = qdims%k_end+1, tdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end

          work_3d(i,j,k) = 0.0

        END DO ! i
      END DO   ! j
    END DO     ! k

    IF (.NOT.ALLOCATED(tinc_ppn)) THEN
      ALLOCATE ( tinc_ppn(qdims%i_end*qdims%j_end,tdims%k_end) )
    END IF

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0,                   &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 181)"//cmessage
    END IF

    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ji = (j-1)*tdims%i_end+i

          tinc_ppn(ji,k) = work_3d(i,j,k)

        END DO ! i
      END DO   ! j
    END DO     ! k

  END IF  !  sf(item,sect)

  item = 182  ! humidity increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = q(i,j,k) - q_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 182)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 183  ! qcl increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = qcl(i,j,k) - qcl_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 183)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 170  ! qcl increment: positive 
  IF (icode <= 0 .AND. sf(item,sect)) THEN 

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MAX(0.0,qcl(i,j,k) - qcl_n(i,j,k)) 
        END DO ! i 
      END DO   ! j 
    END DO     ! k 

! DEPENDS ON: copydiag_3d 
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                & 
         work_3d,                                                       & 
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,   & 
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            & 
         stash_levels,num_stash_levels+1,                               & 
         atmos_im,sect,item,                                            & 
         icode,cmessage) 

     IF (icode >  0) THEN 
       cmessage=": error in copydiag_3d(item 170)"//cmessage 
     END IF 

  END IF  !  sf(item,sect) 

  item = 171  ! qcl increment: negative 
  IF (icode <= 0 .AND. sf(item,sect)) THEN 

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MIN(0.0,qcl(i,j,k) - qcl_n(i,j,k)) 
        END DO ! i 
      END DO   ! j 
    END DO     ! k 

! DEPENDS ON: copydiag_3d 
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                & 
         work_3d,                                                       & 
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,   & 
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            & 
         stash_levels,num_stash_levels+1,                               & 
         atmos_im,sect,item,                                            & 
         icode,cmessage) 

    IF (icode >  0) THEN 
      cmessage=": error in copydiag_3d(item 171)"//cmessage 
    END IF 

  END IF  !  sf(item,sect) 

  item = 184  ! qcf increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = qcf(i,j,k) - qcf_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 184)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 172  ! qcf increment: positive 
  IF (icode <= 0 .AND. sf(item,sect)) THEN 

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MAX(0.0,qcf(i,j,k) - qcf_n(i,j,k)) 
        END DO ! i 
      END DO   ! j 
    END DO     ! k 

! DEPENDS ON: copydiag_3d 
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,   &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage) 

    IF (icode >  0) THEN 
      cmessage=": error in copydiag_3d(item 172)"//cmessage 
    END IF 

  END IF  !  sf(item,sect) 

  item = 173  ! qcf increment: negative 
  IF (icode <= 0 .AND. sf(item,sect)) THEN 

    DO k = 1, qdims%k_end 
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MIN(0.0,qcf(i,j,k) - qcf_n(i,j,k)) 
        END DO ! i 
      END DO   ! j 
    END DO     ! k 

! DEPENDS ON: copydiag_3d 
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,   &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage) 

    IF (icode >  0) THEN 
      cmessage=": error in copydiag_3d(item 173)"//cmessage 
    END IF 

  END IF  !  sf(item,sect) 

  item = 189  ! qrain increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = qrain(i,j,k) - qrain_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 189)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 190  ! qgraup increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = qgraup(i,j,k) - qgraup_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
        work_3d,                                                        &
        qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                  &
        at_extremity,                                                   &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 190)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 191  ! qcf2 increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = qcf2(i,j,k) - qcf2_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 191)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 192  ! cf increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = cf(i,j,k) - cf_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 192)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 193  ! cfl increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = cfl(i,j,k) - cfl_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 193)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 176  ! cfl increment:positive
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MAX(0.0,cfl(i,j,k) - cfl_n(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 176)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 177  ! cfl increment:negative
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MIN(0.0,cfl(i,j,k) - cfl_n(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 177)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 194  ! cff increment
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = cff(i,j,k) - cff_n(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 194)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 178  ! cff increment: positive
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MAX(0.0,cff(i,j,k) - cff_n(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 178)"//cmessage
    END IF

  END IF  !  sf(item,sect)

  item = 179  ! cff increment: negative
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = MIN(0.0,cff(i,j,k) - cff_n(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 179)"//cmessage
    END IF

  END IF  !  sf(item,sect)

! Item 201 Large scale rain

  IF (icode <= 0 .AND. sf(201,4)) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(201,4,im_index)),ls_rain,                &
         qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                 &
         atmos_im,4,201,                                                &
         icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(ls_rain)"
    END IF
! Code to convert rate to amount for a given timestep

    DO i = 1, qdims%i_end*qdims%j_end
      stashwork(si(201,4,im_index)+i-1)=                                &
           stashwork(si(201,4,im_index)+i-1)* timestep
    END DO

  END IF


! Item 202 Large scale snow

  IF (icode <= 0 .AND. sf(202,4)) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(202,4,im_index)),ls_snow,                &
         qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                 &
         atmos_im,4,202,                                                &
         icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(ls_snow)"
    END IF

    DO i = 1, qdims%i_end*qdims%j_end
      stashwork(si(202,4,im_index)+i-1)=                                &
           stashwork(si(202,4,im_index)+i-1)* timestep
    END DO

  END IF


! Item 203 Large scale rain

  IF (icode <= 0 .AND. sf(203,4)) THEN

    IF (.NOT.ALLOCATED(lsrr)) THEN
      ALLOCATE ( lsrr(qdims%i_end*qdims%j_end) )
    END IF

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(203,4,im_index)),ls_rain,                &
         qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                 &
         atmos_im,4,203,                                                &
        icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(ls_rain)"
    END IF

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        ji = (j-1)*qdims%i_end+i
        lsrr(ji) = ls_rain(i,j)
      END DO
    END DO

  END IF


! Item 204 Large scale snow

  IF (icode <= 0 .AND. sf(204,4)) THEN

    IF (.NOT.ALLOCATED(lssr)) THEN
      ALLOCATE ( lssr(qdims%i_end*qdims%j_end) )
    END IF

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(204,4,im_index)),ls_snow,                &
         qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                 &
         atmos_im,4,204,                                                &
         icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(ls_snow)"
    END IF

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        ji = (j-1)*qdims%i_end+i
        lssr(ji) = ls_snow(i,j)
      END DO
    END DO

  END IF



! Items 231-236 mineral dust scavenged by LS precip

  IF (l_dust) THEN

    DO k = 1, ndiv

      IF (icode <= 0 .AND. sf(230+k,4)) THEN

!       Convert "per timestep" diagnostic to "per second":
        DO j = 1, lspice_dim2
          DO i = 1, lspice_dim1
            lscav_dust_all(i, j, k)=lscav_dust_all(i, j, k)/timestep
          END DO
        END DO

! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(230+k,4,im_index)),                  &
           lscav_dust_all(1:qdims%i_end,1:qdims%j_end,k),               &
           qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,               &
           atmos_im,4,230+k,                                            &
           icode,cmessage)

        IF (icode  >   0) THEN
          cmessage=": ERROR IN COPYDIAG(LSCAV_DUST_ALL)"
        END IF

      END IF

    END DO !NDIV

  END IF !L_DUST


! Item 215 NH3 scavenged by LS precip

  IF (icode <= 0 .AND. sf(215,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_nh3(i, j)=lscav_nh3(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(215,4,im_index)),lscav_nh3,              &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,215,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage="microphysics_ctl  : error in copydiag(lscav_nh3)"
    END IF

  END IF

! Item 216 SO2 scavenged by LS precip

  IF (icode <= 0 .AND. sf(216,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_tr(i, j)=lscav_tr(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(216,4,im_index)),lscav_tr,               &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,216,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(lscav_tr)"
    END IF

  END IF

! Item 219 Dissolved SO4 aerosol scavenged by LS precip

  IF (icode <= 0 .AND. sf(219,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        rnout_tracer(i, j)=rnout_tracer(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(219,4,im_index)),rnout_tracer,           &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,219,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(rnout_tracer)"
    END IF

  END IF

! Item 220 Soot scavenged by LS rainout

  IF (icode <= 0 .AND. sf(220,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        rnout_soot(i, j)=rnout_soot(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(220,4,im_index)),rnout_soot,             &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,220,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(rnout_soot)"
    END IF

  END IF

! Item 221 Soot scavenged by LS washout

  IF (icode <= 0 .AND. sf(221,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_soot(i, j)=lscav_soot(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(221,4,im_index)),lscav_soot,             &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,221,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(lscav_soot)"
    END IF

  END IF

! Item 228 Fossil-fuel organic carbon scavenged by LS rainout

  IF (icode <= 0 .AND. sf(228,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        rnout_ocff(i, j)=rnout_ocff(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(228,4,im_index)),rnout_ocff,             &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,228,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(rnout_ocff)"
    END IF

  END IF

! Item 229 Fossil-fuel organic carbon scavenged by LS washout

  IF (icode <= 0 .AND. sf(229,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_ocff(i, j)=lscav_ocff(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(229,4,im_index)),lscav_ocff,             &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,229,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(lscav_ocff)"
    END IF

  END IF

!  Copy T to STASHwork

  IF (icode <= 0 .AND. sf(004,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(004,4,im_index)),t,                   &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0,                   &
         at_extremity,                                                  &
         stlist(1,stindex(1,004,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,004,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(T)"
    END IF

  END IF



!  Copy Cloud water to STASHwork

  IF (icode <= 0 .AND. sf(205,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(205,4,im_index)),qcl,                 &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,205,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,205,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(cloud water)"
    END IF

  END IF

  IF (icode <= 0 .AND. sf(206,4)) THEN

!  Copy Cloud ice to STASHwork

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(206,4,im_index)),qcf,                 &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,206,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,206,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(cloud ice)"
    END IF

  END IF

! -----------------------------
! 207 relative humidity wrt ice (T<0degC) and water (T>0degC) (mdl levs)

  item = 207
  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: qsat_mix
    CALL qsat_mix(work_3d,t,p_layer_centres(1,1,1),                     &
              qdims%i_end*qdims%j_end*qdims%k_end, l_mr_physics1 )
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
          IF (work_3d(i,j,k)  <   0.0) THEN
            work_3d(i,j,k) = 0.0
          END IF

        END DO  ! i
      END DO    ! j
    END DO      ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,     &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 207)"//cmessage
    END IF

  END IF   ! item 207

! -------------------------------
! 208 relative humidity wrt water (on model levels)

  item = 208
  IF (icode  <=  0 .AND. sf(item,sect)) THEN
         ! q saturation is put in work_3d array
! DEPENDS ON: qsat_wat
    CALL qsat_wat(work_3d,t,p_layer_centres(1,1,1),                     &
              qdims%i_end*qdims%j_end*qdims%k_end )
    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
               !  Supersaturation wrt water is limited to =< 100%
          IF (work_3d(i,j,k) > 100.0) THEN
            work_3d(i,j,k) = 100.0
          END IF
               !  Negative humidity also removed from the diagnostic
          IF (work_3d(i,j,k) < 0.0) THEN
            work_3d(i,j,k) = 0.0
          END IF

        END DO  ! i
      END DO    ! j
    END DO      ! k

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)),                &
         work_3d,                                                       &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,     &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d (item 208, rhw)"//cmessage
    END IF

  END IF   ! item 208

!--------------------------------------------------------------------------
! Item 210: 3D droplet number (autoconversion-derived)
!--------------------------------------------------------------------------

! Copy droplet number  to STASHwork

  IF (icode <= 0 .AND. sf(210,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d( stashwork(si(210,4,im_index)), n_drop_3d,         &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,210,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,210,                                                &
         icode,cmessage )

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(n_drop_3d)"
    END IF

  END IF

!--------------------------------------------------------------------------
! Item 211: 3D droplet number (aerosol-derived)
!--------------------------------------------------------------------------

! Copy droplet number  to STASHwork

  IF (icode <= 0 .AND. sf(211,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d( stashwork(si(211,4,im_index)),n_drop_tpr,         &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,211,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,211,                                                &
         icode,cmessage )

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(n_drop_tpr)"
    END IF

  END IF

! Item 213 Ammonium nitrate scavenged by LS rainout
  IF (icode <= 0 .AND. sf(213,4)) THEN
!       Convert "per timestep" diagnostic to "per second":
    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        rnout_nitrate(i, j)=rnout_nitrate(i, j)/timestep
      END DO
    END DO
! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(213,4,im_index)),rnout_nitrate,          &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,213,                                               &
          icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag(rnout_nitrate)"
    END IF
  END IF

! Item 214 Ammonium nitrate scavenged by LS washout
  IF (icode <= 0 .AND. sf(214,4)) THEN
!       Convert "per timestep" diagnostic to "per second":
    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_nitrate(i, j)=lscav_nitrate(i, j)/timestep
      END DO
    END DO
! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(214,4,im_index)),lscav_nitrate,          &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,214,                                               &
          icode,cmessage)
    IF (icode  >   0) THEN
      cmessage=": error in copydiag(lscav_nitrate)"
    END IF

  END IF



!  Copy q  to STASHwork

  IF (icode <= 0 .AND. sf(010,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(010,4,im_index)),q,                   &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,010,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,010,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(q)"
    END IF

  END IF


! Copy ls_rain3d  to STASHwork

  IF (icode <= 0 .AND. sf(222,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(222,4,im_index)),ls_rain3d,           &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,222,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,222,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(ls_rain3d)"
    END IF

  END IF


! Copy ls_snoww3d  to STASHwork

  IF (icode <= 0 .AND. sf(223,4)) THEN

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(223,4,im_index)),ls_snow3d,           &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,223,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,223,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(ls_snow3d)"
    END IF

  END IF


! Need to produce diagnostic 225 before 224 in order to save memory.

  IF(icode <= 0 .AND. sf(225,4)) THEN

!  Supercooled 3D rain content. It is equal to
!  the 3D rainrate at T < 0 and equal to 0 at T > 0
!  Alter the array LS_RAIN3D directly

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF (t(i,j,k)  >=  zerodegc) THEN
! Warm temperatures
            ls_rain3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO

! Copy supercooled rain (now in ls_rain3d)  to STASHwork

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(225,4,im_index)),ls_rain3d,           &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,225,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,225,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(supercooled rain)"
    END IF

  END IF

  IF (icode <= 0 .AND. sf(224,4)) THEN

!  Supercooled liquid water content at TIME LEVEL N. It is equal to
!  the liquid water content at T < 0 and equal to 0 at T > 0
!  Use LS_RAIN3D as the array in order to save memory

    DO k = 1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          IF (t(i,j,k)  <   zerodegc) THEN
! Supercooled temperatures
! Use time level n fields in this diagnostic
            ls_rain3d(i,j,k)=qcl_n(i,j,k)
          ELSE
! Warm temperatures
            ls_rain3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO

! Copy supercooled liquid (now in ls_rain3d)  to STASHwork

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(224,4,im_index)),ls_rain3d,           &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,224,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,224,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(supercooled liq)"
    END IF

  END IF


  IF (icode <= 0 .AND. sf(227,4)) THEN

! Copy rain fraction to stashwork

! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(stashwork(si(227,4,im_index)),rainfrac3d,          &
         qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0,                 &
         at_extremity,                                                  &
         stlist(1,stindex(1,227,4,im_index)),len_stlist,                &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,4,227,                                                &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(rain fraction)"
    END IF

  END IF

! Item 237 Biomass scavenged by LS rainout

  IF (icode <= 0 .AND. sf(237,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        rnout_bmass(i, j)=rnout_bmass(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(237,4,im_index)),rnout_bmass,            &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,237,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(rnout_bmass)"
    END IF

  END IF

! Item 238 Biomass scavenged by LS washout

  IF (icode <= 0 .AND. sf(238,4)) THEN

!       Convert "per timestep" diagnostic to "per second":

    DO i = 1, lspice_dim1
      DO j = 1, lspice_dim2
        lscav_bmass(i, j)=lscav_bmass(i, j)/timestep
      END DO
    END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(238,4,im_index)),lscav_bmass,            &
          qdims%i_end,qdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,4,238,                                               &
          icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(lscav_bmass)"
    END IF

  END IF


        !---------------------------------------------------------------
        ! Homogeneous nucleation
        !---------------------------------------------------------------
  item = 240

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pifrw,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Heterogeneous nucleation
        !---------------------------------------------------------------
  item = 241

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      piprm,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Deposition of ice
        !---------------------------------------------------------------
  item = 243

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pidep,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Deposition of snow aggregates
        !---------------------------------------------------------------
  item = 245

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psdep,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF


        !---------------------------------------------------------------
        ! Ice collection of cloud liquid water (riming)
        !---------------------------------------------------------------
  item = 247

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      piacw,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Snow collection of cloud liquid water (riming)
        !---------------------------------------------------------------
  item = 248

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psacw,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Ice collection of rain (capture)
        !---------------------------------------------------------------
  item = 249

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      piacr,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Snow collection of rain (capture)
        !---------------------------------------------------------------
  item = 250

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psacr,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF
        !---------------------------------------------------------------
        ! Evaporation of melting ice
        !---------------------------------------------------------------
  item = 251

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pimltevp,                                                         &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Evaporation of melting snow
        !---------------------------------------------------------------
  item = 252

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psmltevp,                                                         &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Melting ice
        !---------------------------------------------------------------
  item = 253

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pimlt,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Melting snow
        !---------------------------------------------------------------
  item = 254

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psmlt,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Snow aggregate autoconversion
        !---------------------------------------------------------------
  item = 255

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psaut,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Snow collection of ice (capture)
        !---------------------------------------------------------------
  item = 256

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psaci,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Rain autoconversion
        !---------------------------------------------------------------
  item = 257

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      praut,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Rain collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
  item = 258

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pracw,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Evaporation of rain
        !---------------------------------------------------------------
  item = 259

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      prevp,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Graupel autoconversion
        !---------------------------------------------------------------
  item = 260

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pgaut,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Graupel collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
  item = 261

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pgacw,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Graupel collection of snow (capture)
        !---------------------------------------------------------------
  item =262

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pgacs,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Melting graupel
        !---------------------------------------------------------------
  item = 263

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pgmlt,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Ice crystal sedimentation
        !---------------------------------------------------------------
  item = 265

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pifall,                                                           &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Snow sedimentation
        !---------------------------------------------------------------
  item = 266

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      psfall,                                                           &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Rain sedimentation
        !---------------------------------------------------------------
  item = 267

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      prfall,                                                           &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Graupel sedimentation
        !---------------------------------------------------------------
  item = 268

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pgfall,                                                           &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Droplet settling of liquid
        !---------------------------------------------------------------
  item = 269

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      plset,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Evaporated settled droplets
        !---------------------------------------------------------------
  item = 270

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      plevpset,                                                         &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,*)                                                 &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF

        !---------------------------------------------------------------
        ! Homogeneous nucleation of rain
        !---------------------------------------------------------------
  item = 271

  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (                                                  &
      stashwork(si(item,sect,im_index)),                                &
      pifrr,                                                            &
      qdims%i_end, qdims%j_end, qdims%k_end,0,0,0,0, at_extremity,      &
      stlist(1,stindex(1,item,sect,im_index)),                          &
      len_stlist, stash_levels,num_stash_levels+1,                      &
      atmos_im,sect,item,                                               &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(cmessage,'(A,A,I3)')                                        &
        'ERROR in copydiag_3d for diagnostic ',                         &
        'section 4, item ',item
    END IF
  END IF


! Single point exception handling
  IF (icode /= 0) THEN

    CALL ereport(RoutineName,icode,cmessage)
  END IF

  IF (lhook) CALL dr_hook('DIAGNOSTICS_LSRAIN',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE diagnostics_lsrain
END MODULE diagnostics_lsrain_mod
