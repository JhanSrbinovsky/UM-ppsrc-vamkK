! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_imp_ctl
! *********************************************************************
!
!+ Boundary layer Implicit solver and Large-scale (Area) Cloud Scheme.
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      SUBROUTINE NI_imp_ctl (                                           &

! Parallel variables
     &  global_row_length, global_rows                                  &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, me                    &

! model dimensions.
     &, rhc_row_length, rhc_rows, land_points                           &
     &, ntiles, bl_levels, dst_levels                                   &
     &, dsm_levels, cloud_levels, n_cca_levels, nice, nice_use          &
     &, DIM_CS1, DIM_CS2                                                &
!
! Model switches
     &, model_domain, L_CAL360, L_area_cloud, L_ACF_Cusack              &
     &, L_ACF_Brooks, L_RHCPT, L_emcorr, Ltimer, L_DRY,L_MURK           &
     &, L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_DUST_DIAG                &
     &, L_sulpc_so2, L_sulpc_nh3, L_sulpc_dms, L_soot, L_biomass        &
     &, L_ocff, L_nitrate, L_co2_interactive                            &
     &, L_co2_emits, L_pc2                                              &
     &, NumCycles, CycleNo, lq_mix_bl, L_ukca, L_sice_heatflux          &
     &, L_sice_multilayers                                              &
     &, l_use_cariolle                                                  &
!
! model Parameters
     &, rhcrit, tr_vars, tr_ukca                                        &
     &, co2_mmr                                                         &

! in coordinate information
     &, delta_lambda, delta_phi,lat_rot_NP,long_rot_NP                  &

! in time stepping information.
     &, timestep, val_year, val_day_number, val_hour, val_minute        &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &
     &, f3_at_u                                                         &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork3, STASHwork9                                           &
!
! SCM diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries, rho, rho_wet_theta      &
     &, u, v, w                                                         &
     &, land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels   &
     &, T_conv,q_conv,qcl_conv,qcf_conv                                 &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, smvccl, smvcwt, smvcst, sthf, sthu, sil_orog_land, ho2r2_orog   &
     &, di, ice_fract, di_ncat, ice_fract_ncat, k_sice                  &
     &, u_0, v_0,land_index                                             &
     &, cca,  ccb,  cct,  ccw,  cca_2d, lcbase                          &
      , cca0, ccb0, cct0, ccw0, cca0_2d                                 &
     &, ls_rain, ls_snow, conv_rain, conv_snow, L_scrn, L_plsp          &

! in variables required from BDY_LAYR in IMP_SOLVER
     &, alpha1_sice, ashtf, BQ_GB, BT_GB, dtrdz_charney_grid            &
     &, rdz_charney_grid, dtrdz_u, dtrdz_v, rdz_u, rdz_v                &
     &, cdr10m_u, cdr10m_v, z1_tq                                       &
     &, uStarGBM                                                        &
!ajm    extra variable added
     &, rhokm_u,rhokm_v                                                 &

! in diagnostics (started or from) BDY_LAYR
     &, T_incr_diag_bl, q_incr_diag_bl, qcl_incr_diag_bl                &
     &, qcf_incr_diag_bl, u_incr_diag_bl, v_incr_diag_bl                &
     &, cfl_incr_diag_bl, cff_incr_diag_bl, cf_incr_diag_bl             &
     &, e_sea, fqT, ftl, h_sea, rib_gb                                  &
     &, taux, tauy, vshr, zlcl, zht, dzh                                &
     &, bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
     &, bl_type_7                                                       &
     &, z0m_gb, z0m_eff_gb, z0h_eff_gb                                  &
     &, fme, rhokh                                                      &
     &, TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                          &

! in data required to calculate increment diagnostics
     &, theta_star,q_star                                               &
! in logical for scm surface forcing
     &, L_flux_bc                                                       &

! in data required for tracer mixing :
     &, RHO_ARESIST,ARESIST,RESIST_B,R_B_DUST                           &
     &, KENT, WE_LIM, T_FRAC, ZRZI                                      &
     &, KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                      &
     &, ZHSC,Z_HALF,DUST_FLUX,DUST_EMISS_FRAC                           &
     &, U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE                          &
     &, so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
     &, ocff_hilem, ocff_em                                             &

! IN additional variables for MOSES II
      , TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY                            &
     &, ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,RESFS,Z0HSSI,Z0MSSI         &
     &, CANHC_TILE,FLAKE,WT_EXT_TILE,LW_DOWN,lai_ft,canht_ft            &
     &, SW_TILE,ASHTF_TILE,gc,aresist_tile,resist_b_tile                &
     &, FQT_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT              &
     &, RHOKPM_SICE,Z0H_TILE,Z0M_TILE,CHR1P5M_SICE                      &
     &, FLAND,FLANDG,FLANDG_U,FLANDG_V,TSTAR_SEA,VSHR_LAND,VSHR_SSI     &

! IN additional variables for JULES
     &,RHOKH_MIX_DUMMY,DTSTAR_TILE,DTSTAR,HCONS,EMIS_TILE,EMIS_SOIL     &

! IN MOSES II variables for STASH
     &, GS,GPP,NPP,RESP_P,GPP_FT,NPP_FT,RESP_P_FT,RESP_S                &
     &, RESP_S_TOT,CS                                                   &
     &, RIB_TILE,FSMC,CATCH,G_LEAF                                      &
     &, CO2_EMITS, CO2FLUX                                              &

!IN additional variables for soil moisture nudging scheme
     &, WT_EXT,RA                                                       &

! in/out
! (Note ti and ti_gb are IN only if l_sice_multilayers=T)
     &, t_soil, ti, t_surf, ti_gb                                       &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, cf_latest, cfl_latest, cff_latest                               &
     &, R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen     &
     &, zh,sum_eng_fluxes,sum_moist_flux,rhcpt                          &

! In/Out tracer fields
     &, aerosol, free_tracers                                           &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, DRYDEP2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3          &
     &, soot_new, soot_aged, soot_cld, bmass_new, bmass_agd             &
     &, bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss   &
     &, co2, ozone_tracer                                               &

! in/out fields
     &, ecan, ei, ext, snowmelt                                         &
     &, t1_sd, q1_sd, ntml, cumulus, l_pc2_diag_sh_pts                  &
     &, nbdsc, ntdsc                                                    &
     &, surf_ht_flux_land, snomlt_surf_htf, cH_term                     &

! INOUT additional variables for MOSES II
     &, TSTAR_TILE,FQT_TILE,EPOT_TILE,FTL_TILE                          &
     &, SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR                   &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SICE_CAT,TSTAR_SSI                  &
     &, rib_ssi,taux_land,taux_ssi,tauy_land,tauy_ssi                   &

! OUT additional variables for MOSES II
     &, ESOIL_TILE,ES,EI_TILE                                           &
     &, Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE                       &
     &, SURF_HTF_TILE                                                   &

! error information
     &, Error_code,BL_diag)

! purpose: Interface to boundary layer implicit solver and cloud scheme
!

! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE nstypes
  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, udims_s, udims_l, vdims, vdims_s, vdims_l,                    &
   tdims, tdims_s, tdims_l, qdims, qdims_s, qdims_l,                    &
   pdims, pdims_s, pdims_l, wdims, wdims_s, trdims_ltl

  USE earth_constants_mod, ONLY: g, earth_radius, two_omega

  USE atmos_constants_mod, ONLY: cp, r, repsilon

  USE water_constants_mod, ONLY: lc 

  USE cv_run_mod, ONLY:                                             &
      tice, rad_cloud_decay_opt, l_ccrad, l_3d_cca, l_pc2_diag_sh 

  USE cv_param_mod, ONLY:                                           &
      rad_decay_off

  USE cloud_inputs_mod, ONLY:                                       &
      l_fixbug_pc2_mixph, forced_cu

  USE bl_option_mod, ONLY:                                          &
      Fric_heating, alpha_cd,trweights1,  l_use_bl_diag_term

  USE bl_diags_mod, ONLY:                                           &
      strnewbldiag

  USE swapable_field_mod, ONLY :                                    &
      swapable_field_pointer_type

  USE switches, ONLY:                                               &
      IScrnTDiag

  USE dust_parameters_mod, ONLY: ndiv, ndivh, rhop, drep,           &
      l_twobin_dust

  USE pc2_constants_mod,   ONLY: ls_bl0

  USE level_heights_mod, ONLY:                                      &
      r_theta_levels, r_rho_levels    

  USE turb_diff_mod, ONLY:                                          &
      l_subfilter_vert, l_subfilter_horiz

  USE rad_pcf, ONLY: ip_cloud_mix_max

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  USE Submodel_Mod
  USE Field_Types
  USE UM_ParParams

  USE u_to_p_mod, ONLY: u_to_p
  USE tr_mix_mod, ONLY: tr_mix
  IMPLICIT NONE

! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
     &  global_row_length                                               &
                           ! number of points on a row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, me         ! My processor number

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

      Logical                                                           &
     &  at_extremity(4)                                                 &
!                          Indicates if this processor is at north,
!                          south, east or west of the processor grid
!                          (array index PNorth etc. from parparm.h)
! Switch for  surface flux forcing in SCM
     &, L_flux_bc

! Model dimensions
      Integer                                                           &
     &  rhc_row_length                                                  &
     &, rhc_rows                                                        &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, ntiles                                                          &
                    ! IN No. of land-surface tiles ( MOSES II )
     &, nice                                                            &
                    ! IN No. of sea ice categories
     &, nice_use                                                        &
                    ! IN No. of sea ice categories used fully in surface
                    !    exchange code
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, cloud_levels                                                    &
     &, n_cca_levels
                      ! Number of levels for conv cloud amount :
!                       1 for 2D, nlevs for 3D.
       Integer                                                          &
     &  NumCycles                                                       &
                   ! Number of cycles
     &, CycleNo                                                         &
                   ! Sweep number

     &, DIM_CS1, DIM_CS2
!
! Model switches
      Integer                                                           &
     &  model_domain

      Logical                                                           &
     &  L_CAL360                                                        &
                    ! true if using 360 day calender
     &, L_emcorr                                                        &
                    ! true if energy correction scheme is to be used.
     &, L_dry
                    ! true if model to be run with no moisture

      Logical                                                           &
     &  L_area_cloud                                                    &
                           ! Switch for area cloud fraction param
     &, L_ACF_Cusack                                                    &
                           ! ... to select Cusack
     &, L_ACF_Brooks                                                    &
                           ! ... to select Brooks
     &, L_RHCPT                                                         &
                 ! Switch for 3D diagnosed RHcrit not 1D parameter
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_Murk                                                          &
                          ! Switch for (visibility) aerosol
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! Switch for prognostic mineral dust
     &, L_DUST_DIAG                                                     &
                          ! Switch for diagnostic mineral dust lifting
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for fossil-fuel OC aerosol
     &, L_nitrate                                                       &
                          ! Switch for ammonium nitrate aerosol
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_pc2                                                           &
                          ! Switch for PC2 cloud scheme
     &, L_co2_emits                                                     &
                          ! Switch to include surface emissions
     &, L_ukca                                                          &
                          ! Switch for UKCA sub-model
     &, L_sice_heatflux                                                 &
                          ! Switch for semi-impl sea-ice temperature
     &, l_use_cariolle                                                  &
                          ! Switch for cariolle ozone tracer scheme
     &, L_sice_multilayers
                          ! coupled to multilayers sea ice model
              

! model parameters
      Real                                                              &
     &  timestep

      Real                                                              &
     &  rhcrit(qdims%k_end)                                             &
                                  ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
     &, alpha_tr(bl_levels)  ! Implicit weights for tracers
!                            !  = alpha_cd from RUN_BL namelist
!                            ! or = 1 if TRWEIGHTS1 = ON

      REAL, INTENT(IN)  :: co2_mmr
                             ! Initial or fixed concentration of CO2


      Integer                                                           &
     &  tr_vars                                                         &
                             ! IN number of free tracer variables
     &, tr_ukca              ! IN number of ukca tracer variables

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
! Start blopt8a

! Description:
!   Permissible settings for BL options.


      INTEGER, PARAMETER :: off = 0  ! Switch disabled
      INTEGER, PARAMETER :: on  = 1  ! Switch enabled

!     Options for non-gradient stress following
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!       Brown and Grant (1997), version 2 including a limit on its size

!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)

!     Options for entrainment enhancement in Sc over Cu
      INTEGER, PARAMETER :: Buoyrev_feedback = 1

!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
      INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used as a dynamic criterion in the
!       diagnosis of BL types: version 2 includes changes to
!       cope with BL_LEVELS >> 3km
      INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
!       As 2 but only applied to points diagnosed with Cumulus 
!       and strictly for sea points (fland<0.01, cf 0.5)
      INTEGER, PARAMETER :: DynDiag_Ribased = 4
!       As 3 but also overrides Cumulus diagnosis if 
!          ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
!       Note that here Ri accounts for gradient adjustment by the 
!       non-local scheme.

!     Options for surface exchange
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
      INTEGER, PARAMETER :: Limit_ObukhovL = 3
!       Option under the COR_MO_ITER switch for imposing
!       lower limit on L in very stable conditions.
      INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
      INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

!     Options for convective boundary layers
      INTEGER, PARAMETER ::  UM_std     = 0
      INTEGER, PARAMETER ::  neut_cbl   = 1
      INTEGER, PARAMETER ::  LEM_conven = 2
      INTEGER, PARAMETER ::  LEM_std    = 3

!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a
!     File to set heights for screen diagnostics and options for
!     diagnosis.
!
      REAL Z_OBS_TQ,Z_OBS_WIND
      PARAMETER (                                                       &
     & Z_OBS_TQ = 1.5                                                   &
                         ! Height of screen observations of temperature
!                        ! and humidity.
     &,Z_OBS_WIND = 10.0                                                &
                         ! Height of surface wind observations.
     &)
!
      INTEGER, PARAMETER :: IP_ScrnSurfSim = 0
!                           ! Diagnose the screen temperature using
!                           ! pure surface similarity theory
      INTEGER, PARAMETER :: IP_ScrnDecpl1 = 1
!                           ! Diagnose the screen temperature using 
!                           ! surface similarity theory, but allow 
!                           ! decoupling in very stable conditions
!                           ! based on the quasi-equilibrium radiative
!                           ! solution.
      INTEGER, PARAMETER :: IP_ScrnDecpl2 = 2
!                           ! Diagnose the screen temperature using 
!                           ! including transient effects and radiative
!                           ! cooling

! Local variables

      Integer, Parameter :: sect = 3               ! BL section
      Integer            :: item                   ! stash item code
      Integer            :: im_index               ! internal model

      Character (Len=*), Parameter :: RoutineName='bl_w_mixing'
      Character (Len=80)           :: Cmessage

! Diagnostics info
       Real                                                             &
     & STASHwork3(*)                                                    &
                     ! STASH workspace for section 3 (Boundary Layer)
     &,STASHwork9(*) ! STASH workspace for section 9 (LS Cloud)


! Co-ordinate arrays
      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi


! Trig arrays
      real                                                              &
     &  cos_theta_longitude (tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end)                 &
     &, sin_theta_longitude (tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end)                 &
     &, FV_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,           &
                               tdims_s%j_start:tdims_s%j_end)           &
     &, f3_at_u (udims_s%i_start:udims_s%i_end,                         &
     &                         udims_s%j_start:udims_s%j_end)

! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

!      Real
!     &  p_theta_levels(1-offx :row_length+offx ,
!     &                 1-offy :rows+offy , model_levels)
!     &, exner_rho_levels(1-offx :row_length+offx ,
!     &                   1-offy :rows+offy , model_levels)

! Data arrays
      Real                                                              &
     &  u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
     &      udims_s%k_start:udims_s%k_end)                              &
     &, v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
     &      vdims_s%k_start:vdims_s%k_end)                              &
     &, w(wdims_s%i_start:wdims_s%i_end,wdims_s%j_start:wdims_s%j_end,  &
          wdims_s%k_start:wdims_s%k_end)                                &
     &, rho(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,                              &
     &      pdims_s%k_start:pdims_s%k_end)                              &
     &, rho_wet_theta(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,                        &
                      1:tdims%k_end-1)                                  &
     &, p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,  &
     &    pdims_s%k_start:pdims_s%k_end)                                &
     &, p_layer_centres(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,0:tdims%k_end)        &
     &, p_layer_boundaries(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,0:pdims%k_end)     &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
     &, theta(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
     &, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end)               &
     &, q(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,  &
          qdims_l%k_start:qdims_l%k_end)                                &
     &, qcl(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,&
            qdims_l%k_start:qdims_l%k_end)                              &
     &, qcf(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,&
            qdims_l%k_start:qdims_l%k_end)                              &
     &, T_conv(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_end)                                           &
     &, q_conv(qdims_s%i_start:qdims_s%i_end,                           &
               qdims_s%j_start:qdims_s%j_end,                           &
               qdims_s%k_end)                                           &
     &, qcl_conv(qdims_s%i_start:qdims_s%i_end,                         &
                 qdims_s%j_start:qdims_s%j_end,                         &
                 qdims_s%k_end)                                         &
     &, qcf_conv(qdims_s%i_start:qdims_s%i_end,                         &
                 qdims_s%j_start:qdims_s%j_end,                         &
                 qdims_s%k_end)                                         &
     &, ls_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
     &, ls_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
     &, conv_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
     &, conv_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)


      logical                                                           &
        land_sea_mask(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end)                        &
     &, L_scrn                                                          &
                                 ! Logical to control output
                                 !    of screen level T,Q,QCL,QCF
     &, L_plsp                   ! Logical to control output
                                 !    of Probability of LS Precip

! ancillary arrays and fields required to be saved from timestep to
! timestep.
      Integer                                                           &
     &  land_index (land_points)      ! set from land_sea_mask

      Real                                                              &
     &  u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end)        &
                                ! set to zero
     &, v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)        &
                                ! set to zero
     &, sil_orog_land (land_points)                                     &
                                    ! orog/qrparm.orog.as
     &, ho2r2_orog (land_points)                                        &
                                 ! orog/qrparm.orog.h2root2
     &, smvccl (land_points)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_points)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_points)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_points,dsm_levels)                                    &
                                ! IN Frozen soil moisture content of
                                !     each layer as a fraction of
                                !     saturation.
     &, sthu(land_points,dsm_levels)                                    &
                                ! IN Unfrozen soil moisture content
                                !    of each layer as a fraction of
                                !    saturation.
     &, ice_fract (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
                                     ! ice/qrclim.ice.(month)
     &, di(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
                              ! ice/qrclim.ice_thick.(month)
     &, ice_fract_ncat (pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end, nice)                &
                              !ice fract on categories
     &, di_ncat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                nice)                                                   &
                              ! ice thickness on categories
     &, k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice)
                              ! sea ice effective conductivity in
                              !  sfc layer on categories (W/m2/K)

      Real ::                                                           &
       cca    (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,     &
               n_cca_levels)                                            &
     , ccw    (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,     &
               qdims%k_end)                                             &
     , cca_2d (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)     &
     , cca0   (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,     &
               n_cca_levels)                                            &
     , ccw0   (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,     &
               qdims%k_end)                                             &
     , cca0_2d(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)

      Integer  ::                                                       &
       lcbase (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)     &
      ,ccb    (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)     &
      ,cct    (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)     &
      ,ccb0   (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)     &
      ,cct0   (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)

! in data required to calculate increment diagnostics
      Real                                                              &
     & theta_star(tdims_s%i_start:tdims_s%i_end,                        &
                  tdims_s%j_start:tdims_s%j_end,                        &
                  tdims_s%k_start:tdims_s%k_end)                        &
     &,    q_star(qdims_s%i_start:qdims_s%i_end,                        &
                  qdims_s%j_start:qdims_s%j_end,                        &
                  qdims_s%k_start:qdims_s%k_end)

! in variables passed from BDY_LAYR to IMP_SOLVER
      Real                                                              &
     &  alpha1_sice(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end,nice_use)                 &
     &, ashtf(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use)                       &
     &, BQ_GB(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels)                                                &
     &, BT_GB(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels)                                                &
     &, dtrdz_charney_grid(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,                   &
                           bl_levels)                                   &
     &, rdz_charney_grid(pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end,                     &
                         bl_levels)                                     &
     &, dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,    &
     &                          bl_levels)                              &
     &, dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,    &
     &                            bl_levels)                            &
     &, rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,      &
     &                        2:bl_levels)                              &
     &, rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,      &
     &                          2:bl_levels)                            &
     &, cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
     &, cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
     &, z1_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
      REAL, INTENT(IN) :: uStarGBM(pdims%i_start:pdims%i_end,           &
                                    pdims%j_start:pdims%j_end)
!       ! GBM surface friction velocity for diagnosis of decoupling
!ajm variable added
      Real                                                              &
     &  rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,    &
     &                           bl_levels)                             &
     &, rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,    &
     &                           bl_levels)

! in diagnostics (started or from) BDY_LAYR
      Real                                                              &
            ! output from bdy_layr.
     &  u_incr_diag_bl(udims%i_start:udims%i_end,                       &
          udims%j_start:udims%j_end,udims%k_start:udims%k_end)          &
                              ! u wind      increment for STASH
     &, v_incr_diag_bl(vdims%i_start:vdims%i_end,                       &
          vdims%j_start:vdims%j_end,vdims%k_start:vdims%k_end)          &
                              ! v wind      increment for STASH
     &, T_incr_diag_bl(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,tdims%k_end)           &
     &, q_incr_diag_bl(qdims%i_start:qdims%i_end,                       &
                       qdims%j_start:qdims%j_end,qdims%k_end)           &
     &, qcl_incr_diag_bl(qdims%i_start:qdims%i_end,                     &
                         qdims%j_start:qdims%j_end,qdims%k_end)         &
     &, qcf_incr_diag_bl(qdims%i_start:qdims%i_end,                     &
                         qdims%j_start:qdims%j_end,qdims%k_end)         &
     &, cf_incr_diag_bl(qdims%i_start:qdims%i_end,                      &
                        qdims%j_start:qdims%j_end,qdims%k_end)          &
     &, cfl_incr_diag_bl(qdims%i_start:qdims%i_end,                     &
                         qdims%j_start:qdims%j_end,qdims%k_end)         &
     &, cff_incr_diag_bl(qdims%i_start:qdims%i_end,                     &
                         qdims%j_start:qdims%j_end,qdims%k_end)         &
     &, e_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                 ! needed as diagnostic ?
     &, fqT(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels)                                                  &
                                          ! needed as diagnostic ?
     &, ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels)                                                  &
                                          ! needed as diagnostic
     &, h_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                 ! needed as diagnostic ?
     &, rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                     ! Mean bulk Richardson number for
!                                     lowest layer.
     &, taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
                               bl_levels)                               &
                                           ! needed as diagnostic
     &, tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
                                 bl_levels)                             &
                                             ! needed as diagnostic
     &, vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
     &, zlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
     &, zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                                      ! Max height of turb mixing
     &, dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                                  ! IN inversion thickness
     &, bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
     &, bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
                                  ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
                                  ! IN Effective grid-box roughness
                                  !     lengths for diagnostics
     &, fme(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
     &, rhokh (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               bl_levels)

      REAL, INTENT(INOUT) :: TScrnDcl_SSI(pdims%i_start:pdims%i_end,    &
                                          pdims%j_start:pdims%j_end)
!                           !    Decoupled screen-level temperature
!                           !    over sea or sea-ice
      REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_points,ntiles)
!                           !    Decoupled screen-level temperature
!                           !    over land tiles
      REAL, INTENT(INOUT) :: tStbTrans(pdims%i_start:pdims%i_end,    &
                                          pdims%j_start:pdims%j_end)
!                           !    Time since the transition


      REAL :: stress,stress_1

      Real                                                              &
     &  rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)&
                                ! OUT RHOSTAR*CD_STD*VSHR
                                !     for CLASSIC aerosol scheme
     &, aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
                                ! OUT 1/(CD_STD*VSHR)
                                !     for CLASSIC aerosol scheme
     &, resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                                ! OUT (1/CH-1/(CD_STD)/VSHR
                                !     for CLASSIC aerosol scheme
     &, R_B_DUST(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                 ndiv)                                                  &
                                        ! IN surf layer res for dust
     &, WE_LIM(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)   &
                                    ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)     &
                                    ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)   &
                                    ! IN a fraction of the timestep
     &, WE_LIM_DSC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)&
!                                   ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI_DSC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3) &
                                    ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)&
!                                   ! IN a fraction of the timestep
     &, Z_HALF(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               bl_levels)                                               &
                                    ! IN Z_HALF(*,K) is height of half
!                                   !    level k-1/2.
     &, ZHSC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
                                    ! IN Top of decoupled layer
     &, RA(land_points)                                                 &
                                    ! IN aerodynamic resiatnce (s/m)
     &, WT_EXT(land_points,dsm_levels) ! IN cumulative fract of trans'n

!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag


      Integer                                                           &
     &  KENT(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
                                    ! IN grid-level of SML inversion
     &, KENT_DSC(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! IN grid-level of DSC inversion


! Mineral dust source flux for tracer mixing
      REAL, INTENT(IN) ::                                               &
     &  DUST_FLUX(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  ndiv)                                                 &
     &, DUST_EMISS_FRAC(LAND_POINTS,NTILES)                             &
                                ! OUT fraction of tile can emit dust
     &, U_S_T_TILE(LAND_POINTS,NTILES,NDIVH)                            &
                                           !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_POINTS,NTILES,NDIVH)                        &
                                               !OUT dry soil value
     &, U_S_STD_TILE(LAND_POINTS,NTILES)!OUT friction velocity

! Emissions for tracer mixing
      Real, Intent(In) ::                                               &
     &  so2_hilem (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, so2_em    (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, nh3_em    (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, dms_em    (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, soot_hilem(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, soot_em   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, ocff_hilem(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
     &, ocff_em   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! arguments with intent in/out. ie: input variables changed on output.
      Real                                                              &
     &  T_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
     &, t_soil(land_points,dsm_levels)                                  &
                                       ! slt/qrclim.slt_pm(lev).(month)
     &, ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)    &
                                       ! category sea ice sfc layer temp
                                       ! (IN only if l_sice_multilayers=T)
     &, R_u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,&
     &        udims_s%k_start:udims_s%k_end)                            &
     &, R_v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,&
     &        vdims_s%k_start:vdims_s%k_end)                            &
     &, R_w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,        &
            wdims%k_start:wdims%k_end)                                  &
     &, T_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 tdims%k_end)                                           &
     &, q_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,   &
                 qdims%k_end)                                           &
     &, qcl_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qdims%k_end)                                         &
     &, qcf_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qdims%k_end)                                         &
     &, qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                 tdims%k_end)                                           &
     &, cf_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,  &
                  qdims%k_end)                                          &
     &, cfl_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qdims%k_end)                                         &
     &, cff_latest(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, &
                   qdims%k_end)                                         &
     &, area_cloud_fraction(qdims%i_start:qdims%i_end,                  &
                            qdims%j_start:qdims%j_end,                  &
                            qdims%k_end)                                &
     &, bulk_cloud_fraction(qdims%i_start:qdims%i_end,                  &
                            qdims%j_start:qdims%j_end,                  &
                            qdims%k_end)                                &
     &, cloud_fraction_liquid(qdims%i_start:qdims%i_end,                &
                              qdims%j_start:qdims%j_end,                &
                              qdims%k_end)                              &
     &, cloud_fraction_frozen(qdims%i_start:qdims%i_end,                &
                              qdims%j_start:qdims%j_end,                &
                              qdims%k_end)                              &
     &, zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
     &, sum_eng_fluxes(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end)                       &
     &, sum_moist_flux(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end)

! Tracer variables
      Real, Intent(InOut) ::                                            &
     &  aerosol     (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)                     &
     & ,free_tracers(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
     &               trdims_ltl%k_end, tr_vars)

      Real, Intent(InOut) ::                                            &
     &  DUST_DIV1   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DUST_DIV2   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DUST_DIV3   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DUST_DIV4   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DUST_DIV5   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DUST_DIV6   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     &, DRYDEP2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                NDIV) !dry dep though grav. set.

      Real, Intent(InOut) ::                                            &
     &  so2         (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,dms         (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,so4_aitken  (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,so4_accu    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,so4_diss    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,nh3         (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )

      Real, Intent(InOut) ::                                            &
     &  soot_new    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,soot_aged   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,soot_cld    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,bmass_new   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,bmass_agd   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,bmass_cld   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,ocff_new    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,ocff_aged   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,ocff_cld    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,nitr_acc    (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,nitr_diss   (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,co2         (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )                    &
     & ,ozone_tracer(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end )

      Real                                                              &
             !
     &  ecan(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
                                        !output from sf_evap.
     &, ei(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
                                        !output from sf_evap.
     &, ext(land_points,dsm_levels)                                     &
                                    ! Extraction of water from each
!                                    soil layer (kg/m2/s).
     &, snowmelt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                                        !output from sf_evap.
     &, t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                                ! set to zero initially
     &, q1_sd(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)      &
                                ! set to zero initially
     &, surf_ht_flux_gb(pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end)                      &
                                           !
     &, snomlt_surf_htf(pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end)

      Logical ::                     &
       cumulus (qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)    &
                                   ! *APL bl convection flag
     , l_pc2_diag_sh_pts(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end) 
                ! Carry diagnostic shallow convective information for PC2

      Integer                                                           &
     &  ntml (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
     &, nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
     &, ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

      Real                                                              &
     & cH_term(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,                           &
               pdims_s%k_start:pdims_s%k_end-1)

! IN additional variables for MOSES II

      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                 ! IN Number of tile points.
     &,TILE_INDEX(LAND_POINTS,NTYPE)
!                                ! IN Index of tile points.

      REAL                                                              &
     & TILE_FRAC(land_points,ntiles)                                    &
                                ! IN fractional coverage for each
                                !    surface tile
     &,CANOPY(land_points,ntiles)                                       &
                                   ! IN Surface/canopy water (kg/m2)
     &,ALPHA1(land_points,ntiles)                                       &
                                  ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
     &,FRACA(land_points,ntiles)                                        &
                                   ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
     &,RHOKH_TILE(land_points,ntiles)                                   &
                                      ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
     &,SMC(LAND_POINTS)                                                 &
                                ! IN Soil moisture content in root depth
!                                  (kg/m2).
     &,CHR1P5M(LAND_POINTS,NTILES)                                      &
                                   ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
     &,RESFS(land_points,ntiles)                                        &
                                 ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
      ,Z0HSSI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
      ,Z0MSSI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                ! IN Roughness lengths over sea (m).
      ,Z0M_GB(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                ! IN Gridbox mean roughness length 
!                               !    for momentum (m).
     &,CANHC_TILE(LAND_POINTS,NTILES)                                   &
!                               ! IN Areal heat capacity of canopy
!                               !    for land tiles (J/K/m2).
     &,FLAKE(LAND_POINTS,NTILES)                                        &
                                   ! IN Lake fraction.
     &,WT_EXT_TILE(LAND_POINTS,DSM_LEVELS,NTILES)                       &
!                               ! IN Fraction of evapotranspiration
!                               !    which is extracted from each
!                               !    soil layer by each tile.
     &,LW_DOWN(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                ! IN Surface downward LW radiation
!                               !    (W/m2).
     &,lai_ft(land_points,npft)                                         &
                                   ! IN LAI on vegetated tiles
     &,canht_ft(land_points,npft)                                       &
                                   ! IN CANHT on vegetated tiles
     &,SW_TILE(LAND_POINTS,NTILES)                                      &
                                   ! IN Surface net SW radiation on land
!                               !    tiles (W/m2).
     &,ASHTF_TILE(LAND_POINTS,NTILES)                                   &
!                               ! IN Coefficient to calculate
!                               !    surface heat flux into land
!                               !    tiles.
     &,FQT_ICE(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use)                      &
                                ! IN Surface FQT for sea-ice
     &,FTL_ICE(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use)                      &
                                ! IN Surface FTL for sea-ice
     &,RESFT(LAND_POINTS,NTILES)                                        &
                                   ! IN Total resistance factor.
!                               !    FRACA+(1-FRACA)*RESFS for
!                               !    snow-free land, 1 for snow.
      ,RHOKH_SICE(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
!                               ! IN Surface exchange coefficients
!                               !    for sea and sea-ice
     &,RHOKPM(LAND_POINTS,NTILES)                                       &
                                ! IN Land surface exchange coeff.
     &,RHOKPM_POT(LAND_POINTS,NTILES)                                   &
!                               ! IN Land surface exchange coeff.
!                                    for potential evaporation.
     &,EPOT_TILE(land_points,ntiles)                                    &
!                               ! INOUT surface tile potential
!                               !       evaporation
      ,RHOKPM_SICE(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) &
!                               ! IN Sea-ice surface exchange coeff.
     &,Z0H_TILE(LAND_POINTS,NTILES)                                     &
                                ! IN Tile roughness lengths for heat
!                               !    and moisture (m).
     &,Z0M_TILE(LAND_POINTS,NTILES)                                     &
                                ! IN Tile roughness lengths for
!                               !    momentum.
      ,CHR1P5M_SICE(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)&
!                               ! IN CHR1P5M for sea and sea-ice
!                               !    (leads ignored).
     &,FLAND(LAND_POINTS)                                               &
                                ! IN Land fraction on land tiles.
     &,FLANDG(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end)                            &
                                ! IN Land fraction on all points.
     &,FLANDG_U(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
!                               ! IN Land frac (on U-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,FLANDG_V(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)    &
!                               ! IN Land frac (on V-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,TSTAR_SEA(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
!                               ! IN Open sea sfc temperature (K).
      , vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
!                               ! IN VSHR over land part of gridbox.
      , vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                                ! IN VSHR over sea part of gridbox.
     &, gc(land_points,ntiles)                                          &
                                  !Stomatal conductance to evapn
!                                 !    for land tiles (m/s).
     &, aresist_tile(land_points,ntiles)                                &
!                                  !1/(CD_STD*VSHR) on land tiles
!                                  !for CLASSIC aerosol scheme
     &, resist_b_tile(land_points,ntiles)
!                                  !(1/CH-1/CD_STD)/VSHR on land tiles
!                                  !for CLASSIC aerosol scheme

! Additional variables for JULES

      REAL ::                                                           &
     & RHOKH_MIX_DUMMY(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end)                       &
                                   ! IN  Exchange coeffs for moisture.
     &,DTSTAR_TILE(LAND_POINTS,NTILES)                                  &
                                   ! IN  Change in TSTAR over timestep 
!                                  !     for land tiles
     &,DTSTAR(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use)                       &
                                   ! IN  Change is TSTAR over timestep
!                                  !     for sea-ice
     &,HCONS(LAND_POINTS)                                               &
                                   ! IN  Soil thermal conductivity
!                                  !     including water and ice
     &,EMIS_TILE(LAND_POINTS,NTILES)                                    &
                                   ! IN  Emissivity for land tiles
     &,EMIS_SOIL(LAND_POINTS)
                                   ! IN  Emissivity of underlying soil

! IN MOSES II additional STASH variables
      REAL                                                              &
     & GS(LAND_POINTS)                                                  &
                                ! IN "Stomatal" conductance to
!                               !    evaporation (m/s).
     &,GPP(LAND_POINTS)                                                 &
                                ! IN Gross primary productivity
!                               !    (kg C/m2/s).
     &,NPP(LAND_POINTS)                                                 &
                                ! IN Net primary productivity
!                               !    (kg C/m2/s).
     &,RESP_P(LAND_POINTS)                                              &
                                ! IN Plant respiration (kg C/m2/s).
     &,GPP_FT(LAND_POINTS,NPFT)                                         &
                                ! IN Gross primary productivity
!                               !    on PFTs (kg C/m2/s).
     &,NPP_FT(LAND_POINTS,NPFT)                                         &
                                ! IN Net primary productivity
!                               !    on PFTs (kg C/m2/s).
     &,RESP_P_FT(LAND_POINTS,NPFT)                                      &
                                  !IN Plant respiration on PFTs
!                               !     (kg C/m2/s).
     &,RESP_S(LAND_POINTS,DIM_CS1)                                      &
                                   ! IN Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              & 
                                   ! IN Total soil resp'n (kg C/m2/s).
     &,CS(LAND_POINTS,DIM_CS1)                                          &
                                   ! IN Soil carbon
     &,RIB_TILE(LAND_POINTS,NTILES)                                     &
!                               ! IN RIB for land tiles.
     &,FSMC(LAND_POINTS,NPFT)                                           &
                                ! IN Moisture availability factor.
     &,CATCH(LAND_POINTS,NTILES)                                        &
                                ! IN Surface/canopy water capacity
!                               !    of snow-free land tiles (kg/m2).
     &,G_LEAF(LAND_POINTS,NPFT)                                         &
                                ! IN Leaf turnover rate (/360days).
      ,CO2_EMITS(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                                 !IN CO2 Emissions
      ,CO2FLUX(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                      ! IN ocean CO2 flux
      ,co2_flux_tot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)&
                                      !  total CO2 flux
     &,land_co2(land_points)          !  terrestrial CO2 flux


! INOUT additional variables for MOSES II
      REAL                                                              &
     & TSTAR_TILE(land_points,ntiles)                                   &
                                ! INOUT Surface tile temperature
     &,FQT_TILE(land_points,ntiles)                                     &
                                ! INOUT surface tile moisture flux
     &,FTL_TILE(land_points,ntiles)                                     &
!                               ! INOUT surface tile heat flux
     &,SNOW_TILE(LAND_POINTS,NTILES)                                    &
!                               ! INOUT Snow on tiles (kg/m2).
     &,LE_TILE(LAND_POINTS,NTILES)                                      &
                                   ! INOUT Surface latent heat flux for
!                               !       land tiles (W/m2).
     &,RADNET_SICE(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,nice_use)                  &
!                               ! INOUT Sea-ice surface net radiation.
     &,RADNET_TILE(LAND_POINTS,NTILES)                                  &
!                               ! INOUT Tile surface net radiation.
      ,OLR(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
     &,TSTAR_LAND(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                                   ! INOUT Land mean sfc temperature (K)
     &,TSTAR_SICE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                                   ! OUT Sea-ice sfc temperature (K).
                                   !     (Ice mean over categories)
     &,TSTAR_SICE_CAT(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,nice_use)               &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                                   ! INOUT Sea mean sfc temperature (K).
      ,RIB_SSI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                   ! INOUT Sea mean bulk Richardson
!                                          number for lowest layer.
     &,TAUX_LAND(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                                   ! INOUT W'ly compt of land sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
     &,TAUX_SSI(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                                   ! INOUT W'ly compt of sea sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
     &,TAUY_LAND(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                                   ! INOUT S'ly compt of land sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.
     &,TAUY_SSI(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! INOUT S'ly compt of sea sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.

! OUT additional variables for MOSES II
      REAL                                                              &
     & ESOIL_TILE(land_points,ntiles)                                   &
                                ! OUT Evaporation from bare soil (kg/m2)
      ,ES(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
                                ! OUT Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EI_TILE(LAND_POINTS,NTILES)                                      &
                                   ! OUT EI for land tiles
     &,Q1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! OUT Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! OUT T1P5M over land tiles.
     &,ECAN_TILE(LAND_POINTS,NTILES)                                    &
                                ! OUT ECAN for land tiles
     &,MELT_TILE(LAND_POINTS,NTILES)                                    &
!                               ! OUT Snowmelt on tiles (kg/m2/s).
     &,SURF_HTF_TILE(LAND_POINTS,NTILES)
!                               ! OUT Net downward surface heat flux
!                               !     on tiles (W/m2)

! Local variables for MOSES II only
      REAL                                                              &
        E_SSI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                    ! Evaporation from mean sea
      , EI_SICE(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,nice_use)                     &
                                    ! Output from sf_evap.
     &, SURF_HT_FLUX_LAND(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end)                    &
     &, SURF_HT_FLUX_SICE(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end,nice)               &
                                    ! Category sea ice surface heat flux  
      , FTL_SSI(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
! Surface FTL for mean sea

!local variables for mineral dust
      REAL DUST_ALL(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,                          &
                    1:tdims%k_end,NDIV) !dust mmr
      REAL T_MODELLEVS(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,                       &
                       tdims%k_end) !T on model levs


      Integer                                                           &
     &  Error_code


! local variables.
      REAL                                                              &
     &  DENOM                                                           &
                   ! Denominator in PC2 inhomogeneous ice forcing calc.
     &, Q4                                                              &
                   ! QCF increment in PC2 inhomog.    ice forcing calc.
     &, InCloudIce                                                      &
                   ! Value of QCF/CFF used for calculating sink of CFF
     &, deltacff                                                        &
                   ! Ice cloud fraction increment
     &, w_int(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,      &
              bl_levels)
                   ! w on interior points          

      Real :: cca_at_base(qdims%i_start:qdims%i_end,                    &
                          qdims%j_start:qdims%j_end)
!
! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &, kinvert                                                         &
                                ! vertical index for inverted arrays.
     &, IDIV

      Integer :: i_field  ! counter for swap_bounds_mv

! Diagnostic switches
! a) boundary layer
      Logical                                                           &
     &  su10, sv10, slh, sq1p5, sT1p5, sq_T1p5                          &
     &, simlt, smlt, l_ftl, l_fqw, l_taux, l_tauy

! local variables
!
! Diagnostics controlled by Diagnostic switches

      Real                                                              &
     &  q1p5m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
     &, t1p5m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
     &, rho1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
                                               ! Density at level 1
     &, u10m(udims%i_start:udims%i_end,udims%j_start:udims%j_end)       &
     &, v10m(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)       &
     &, latent_heat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)&
     &, sice_mlt_htf(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end, nice)                   &
                                               !output seaice topmelt
     &, sea_ice_htf(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end, nice)                    &
                                               !output seaice fcondtop
                                               !(downwards conductive flux
                                               ! used to force ice model) 
     &, ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                               !output seaice temp.
                                               ! sfc layer (ice mean)

      Real                                                              &
     &  rhokh_mix (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   bl_levels)

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs              ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! local variables
      Integer                                                           &
     &  nclds      ! Number of radiation cloud levels ( <=  wet levels)
!
!-------Needed for area_cloud formulation-------------------------------
      Integer                                                           &
     &  levels_per_level                                                &
                               ! 3 is hardwired inside ls_arcld
     &, large_levels           ! depends on above and wet_levels
!
!-----------------------------------------------------------------------

      REAL                                                              &
       qn                                                               &
              ! temporary in forced cloud calculation
     , qcl_forced                                                       &
     , cf_forced                                                        &
              ! forced cloud water content and fraction
     , dqcl                                                             &
     , dcfl                                                             &
              ! forced cloud water content and fraction increments
     , z_theta                                                          &
              ! height of theta levels above the surface
     , cf_base                                                          &
     , cf_top                                                           &
              ! forced cloud fraction at cloud base and top
     , zc_depth                                                         &
              ! forced cloud depth
     , rht                                                              &
     , rhc                                                              &
              !       total and critical RH
     , alpha                                                            &
              !       d qsat/dT (kg kg-1 K-1)
     , al                                                               &
              !       1 / (1 + alpha L/cp)  (no units)
     , bs     !       Width of distribution (kg kg-1)

      REAL, PARAMETER :: lcrcp=lc/cp

! Local data arrays
      Real                                                              &
     &  T(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)&
     &, RHCPT(rhc_row_length, rhc_rows, qdims%k_end)

      Real, Target :: ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1)

      Real                                                              &
     &  u_inc_bl(udims_s%i_start:udims_s%i_end,                         &
                  udims_s%j_start:udims_s%j_end,bl_levels)              &
     &, v_inc_bl(vdims_s%i_start:vdims_s%i_end,                         &
                  vdims_s%j_start:vdims_s%j_end,bl_levels)


      Real                                                              &
     &  interp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     & 
                                 ! Workspace in calculation of Ch_term's
     &, temp1, temp2, temp3                                             &
                                 ! Temporary variables in calculation of
                                 ! Ch_term's.
     &, drydep_str(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

      Real, Target ::                                                   &
     &  ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,            &
     &                                         0:qdims%k_end)           &
     &, ext_TL(0:rhc_row_length+1, 0:rhc_rows+1,qdims%k_end)            &
     &, ext_QL(0:rhc_row_length+1, 0:rhc_rows+1,qdims%k_end)            &
     &, ext_QCF(0:rhc_row_length+1,0:rhc_rows+1,qdims%k_end)

       Real ::                                                          &
     &  work2d_1(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end)   &
                                    ! Single-level work array (cloud)
     &, plsp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! Probability of large-scale precip
      REAL, ALLOCATABLE :: f3_at_p(:, :)
                                    ! Coriolis parameter at theta-points

!
! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      Real,Dimension(:,:,:),Allocatable::                               &
     & combined_cloud                                                   &
                             ! Mixed CCA and CF per gridbox
     &,T_earliest                                                       &
     &,q_earliest                                                       &
     &,qcl_earliest                                                     &
     &,qcf_earliest                                                     &
     &,cf_earliest                                                      &
     &,cfl_earliest                                                     &
     &,cff_earliest                                                     &
     &,T_inc_PC2                                                        &
                        !  temperature     increment due to PC2 homog
     &,q_inc_PC2                                                        &
                        !  humidity        increment due to PC2 homog
     &,qcl_inc_PC2                                                      &
                        !  qCL             increment due to PC2 homog
     &,cfl_inc_PC2                                                      &
                        !  cf_liquid       increment due to PC2 homog
     &,bcf_inc_PC2      !  bulk cloud      increment due to PC2 homog
!
      Real, Dimension (:,:,:), Allocatable ::    &
       zeros              & ! Array of zero values
     , TL_force           & ! Forcing of TL by homogenous processes
     , QT_force           & ! Forcing of QT by homogenous processes
     , ccw_cca            & ! Convective cloud water * frac (i.e. gridbox
                            ! mean ccw)
     , cca_3d               ! 3D array of convective cloud frac

      REAL, ALLOCATABLE ::  &
        cca4comb_cld(:,:,:) &! Used to calculate combined cloud
      , ccb4comb_cld(:,:)   &! Used to calculate combined cloud
      , cct4comb_cld(:,:)    ! Used to calculate combined cloud

!
! STASHflag switches for increment diagnostics:
      Logical                                                           &
     & l_u_incr_bl                                                      &
                             ! u wind
     &,l_v_incr_bl                                                      &
                             ! v wind
     &,L_T_incr_bl_lsc                                                  &
                             ! T across BL and LS CLD
     &,L_Tl_incr_bl_lsc                                                 &
                             ! Tl across BL (and LS CLD)
     &,L_q_incr_bl_lsc                                                  &
                             ! Q across BL and LS CLD
     &,L_qtl_incr_bl_lsc                                                &
                             ! QT (q+qCL) across BL (and LS CLD)
     &,L_qcl_incr_bl_lsc                                                &
                             ! qCL across BL and LS CLD
     &,L_qcf_incr_bl_lsc                                                &
                             ! qCF across BL (and LS CLD)
     &,L_apply_diag                                                     &
                             ! flag to determine when to apply
                             ! diagnostics when iterating
     &,L_cf_incr_bl                                                     &
                            ! total cloud fraction
     &,L_cfl_incr_bl                                                    &
                            ! liquid cloud fraction
     &,L_cff_incr_bl        ! frozen cloud fraction

!
! Switch for field calculations to support STASH diagnostics
      Logical                                                           &
     & L_combi_cld           ! combined cloud amount

      Type(swapable_field_pointer_type) :: fields_to_swap(4) ! mv swapbounds

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
! Temporarily create these as other subroutines use them
! can be removed when all Endgame changes are made
      INTEGER ::                                                        &
    halo_i,halo_j,offx,offy,n_rows
!
!- End of Header
! ----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('NI_IMP_CTL',zhook_in,zhook_handle)
      If ( error_code  ==  0) Then
! ----------------------------------------------------------------------
! Section BL.0 Initialisation of variables.
! ----------------------------------------------------------------------
! Temporarily create these as other subroutines use them
! can be removed when all Endgame changes are made
    halo_i=pdims_l%i_end-pdims%i_end
    halo_j=pdims_l%j_end-pdims%j_end
    offx=pdims_s%i_end-pdims%i_end
    offy=pdims_s%j_end-pdims%j_end
    n_rows=vdims%j_end-vdims%j_start+1
!
        IF ( TRWEIGHTS1  ==  ON ) THEN
!         ! Set all implict weights used by tracers to one
!         !  - overweighting not necessary since tracers
!         !    have no feedback on the diffusion coefficients
          DO K=1,BL_LEVELS
            alpha_tr(k) = 1.0
          ENDDO
        ELSE
!        ! Set implict weights used by tracers to those input
          DO K=1,BL_LEVELS
            alpha_tr(k) = alpha_cd(k)
          ENDDO
        ENDIF

! Apply diags at last cycle only
        L_apply_diag = CycleNo == NumCycles

! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
!        !--------------------------------------------------------
!        ! Note that an equivalent block of code exists in routine
!        ! ni_bl_ctl, and needs to be kept consistent.
!        !--------------------------------------------------------
!        ! Windspeed (227, 230) and u, v at 10m on 'B' or 'C' grid
         su10 = (sf(209,3) .OR. sf(225,3) .OR. sf(227,3) .OR.           &
                           sf(230,3) .OR. sf(463,3)) .AND. l_apply_diag
         sv10 = (sf(210,3) .OR. sf(226,3) .OR. sf(227,3) .OR.           &
                           sf(230,3) .OR. sf(463,3)) .AND. l_apply_diag
         slh = sf(234,3) .AND. L_apply_diag
         sq_T1p5 = ( sf(236,3) .OR. sf(237,3) .OR. sf(245,3)            &
              .OR. sf(247,3) .OR. sf(248,3) .OR. sf(250,3) .OR. L_scrn  &
              .OR. sf(341,3) .OR. sf(342,3)                             &
              .OR. sf(253,3) .OR. sf(328,3) .OR. sf(329,3)              &
                    ) .AND. L_apply_diag
         sq1p5 = sq_T1p5 .AND. L_apply_diag
         sT1p5 = sq_T1p5 .AND. L_apply_diag
         ! Sea ice topmelt (single category or multi category)
         simlt = ( sf(235,3) .OR. sf(257,3) ) .AND. L_apply_diag 
         smlt = sf(258,3) .AND. L_apply_diag
         l_ftl = (Fric_heating /= OFF) .OR.                             &
                 ( sf(216,3) .AND. L_apply_diag )
         l_fqw = (Fric_heating /= OFF) .OR.                             &
                 ( sf(222,3) .AND. L_apply_diag )
         l_taux = (Fric_heating /= OFF) .OR.                             &
                  ( ( sf(219,3) .OR. sf(221,3) .OR. sf(463,3) ) .AND.    & 
                    L_apply_diag  )
         l_tauy = (Fric_heating /= OFF) .OR.                             &
                  ( ( sf(220,3) .OR. sf(221,3) .OR. sf(463,3) ) .AND.    & 
                    L_apply_diag  )
         L_u_incr_bl = sf(185,3) .AND. L_apply_diag
         L_v_incr_bl = sf(186,3) .AND. L_apply_diag
         L_T_incr_bl_lsc = (sf(181,9).OR.sf(181,3)) .AND. L_apply_diag
         L_Tl_incr_bl_lsc = sf(189,3) .AND. L_apply_diag
         L_q_incr_bl_lsc = (sf(182,9).OR.sf(182,3)) .AND. L_apply_diag
         L_qtl_incr_bl_lsc = sf(190,3) .AND. L_apply_diag
         L_qcl_incr_bl_lsc =(sf(183,9).OR.sf(183,3).OR.sf(170,3)        &
                             .OR.sf(171,3)).AND. L_apply_diag
         L_qcf_incr_bl_lsc = (sf(184,3).OR.sf(172,3).OR.sf(173,3))      &
                             .AND. L_apply_diag
         L_cf_incr_bl  = sf(192,3) .AND. L_apply_diag
         L_cfl_incr_bl = (sf(193,3).OR.sf(176,3).OR.sf(177,3))          &
                          .AND. L_apply_diag
         L_cff_incr_bl = (sf(194,3).OR.sf(178,3).OR.sf(179,3))          &
                          .AND. L_apply_diag
!
! Flag required for pre-calculation of cloud-related fields needed for
! cloud and vis diagnostics. Not duplicated in ni_bl_ctl.
         l_combi_cld = ( l_plsp .OR.                                    &
                 sf(208,9) .OR. sf(209,9) .OR. sf(210,9) .OR. sf(211,9) &
            .OR. sf(212,9) .OR. sf(213,9) .OR. sf(214,9) .OR. sf(215,9) &
            .OR. sf(216,9) .OR. sf(217,9) .OR. sf(223,9) .OR. sf(231,9) &
            .OR. sf(232,9) .OR. sf(233,9) .OR. SF(234,9)                &
            ) .AND. l_apply_diag
      ! Note this will be affected by the anvil scheme if it is
      ! applied, i.e the cca at the anvil base may have been
      ! scaled by the anvil tower_factor.
      IF (l_3d_cca) THEN
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            IF (lcbase(i,j) /= 0) THEN
              cca_at_base(i,j) = cca(i,j,lcbase(i,j))
            ELSE
              cca_at_base(i,j) = 0.0
            END IF
          END DO
        END DO
      ELSE
        ! CCA is only dimensioned with one level
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            IF (lcbase(i,j) /= 0) THEN
              cca_at_base(i,j) = cca(i,j,1)
            ELSE
              cca_at_base(i,j) = 0.0
            END IF
          END DO
        END DO
      END IF ! l_3d_cca
         

!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! implicit solver for optional output of bl wind increments
!---------------------------------------------------------------------
! We need to store information about the increments of the temperature
! and moisture variables, so copy these to the _earliest variables.
      If (L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR. L_pc2) Then
!
        Allocate ( T_earliest(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,tdims%k_end) )

! Hold initial value of Temperature
        Do k = 1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              T_earliest(i,j,k) = T_latest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

      End if                   ! on STASHflags or PC2
!
      If(L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc                         &
     &   .OR. L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc                   &
     &   .OR. L_qcf_incr_bl_lsc .OR. L_pc2) Then
!
        Allocate ( q_earliest(qdims%i_start:qdims%i_end,                &
                              qdims%j_start:qdims%j_end,qdims%k_end) )
        Allocate ( qcl_earliest(qdims%i_start:qdims%i_end,              &
                                qdims%j_start:qdims%j_end,qdims%k_end) )
        Allocate ( qcf_earliest(qdims%i_start:qdims%i_end,              &
                                qdims%j_start:qdims%j_end,qdims%k_end) )
        Allocate ( cf_earliest(qdims%i_start:qdims%i_end,               &
                               qdims%j_start:qdims%j_end,qdims%k_end) )
        Allocate ( cfl_earliest(qdims%i_start:qdims%i_end,              &
                                qdims%j_start:qdims%j_end,qdims%k_end) )
        Allocate ( cff_earliest(qdims%i_start:qdims%i_end,              &
                                qdims%j_start:qdims%j_end,qdims%k_end) )
!
! Hold initial values of wet parameters
        Do k=1, qdims%k_end
          Do j=qdims%j_start,qdims%j_end
            Do i=qdims%i_start, qdims%i_end
              q_earliest(i,j,k)   = q_latest(i,j,k)
              qcl_earliest(i,j,k) = qcl_latest(i,j,k)
              qcf_earliest(i,j,k) = qcf_latest(i,j,k)
              cf_earliest(i,j,k)  = cf_latest(i,j,k)
              cfl_earliest(i,j,k) = cfl_latest(i,j,k)
              cff_earliest(i,j,k) = cff_latest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

      ELSE
        Allocate ( q_earliest(1,1,1) )
                              
        Allocate ( qcl_earliest(1,1,1) )
                                
        Allocate ( qcf_earliest(1,1,1) )
                                
        Allocate ( cf_earliest(1,1,1) )
                               
        Allocate ( cfl_earliest(1,1,1) )
                                
        Allocate ( cff_earliest(1,1,1) )
!
      End if                  ! on STASHflags or PC2
!
! ----------------------------------------------------------------------
! Section BL.1 Calculate T at old time level.
! Modified to use latest values to avoid time-level inconsistencies
! with cloud data.
! ---------------------------------------------------------------------

        Do k = 1, bl_levels
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do

!  Initialise output arrays to zero

        Do j= pdims%j_start, pdims%j_end
          Do i= pdims%i_start, pdims%i_end
            ecan(i,j)=0.0
            ei(i,j) = 0.0
            snowmelt(i,j)=0.0
          End Do
        End Do

! ----------------------------------------------------------------------
! Section BL.2  Call Implicit solver
!               Call tracer mixing for qcf.
!               Call tracer mixing for other tracers.
! ----------------------------------------------------------------------
!
!       Density of the lowest level is required for the transient screen
!       diagnostic. (Note: rho1 is recalculated in part 6 for
!       visibility, so a minor efficiency saving could be made.)
        If (IScrnTDiag == IP_ScrnDecpl2) Then
          Do j = pdims%j_start, pdims%j_end
            Do i = pdims%i_start, pdims%i_end
              rho1(i,j)=rho(i,j,1)/(r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
            End Do
          End Do
!         Allocate the Coriolis parameter on the p-grid and interpolate
!         from the u-grid.
          ALLOCATE(f3_at_p(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end))
      CALL u_to_p (f3_at_u,                                             & 
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        1,                                              &
                        model_domain,at_extremity, f3_at_p)     

IF (.NOT. l_vatpoles) THEN     
!         This routine does not fill the poles in a global model, so
!         the numbers are filled in by hand.
          If (model_domain  ==  mt_global) Then

            If (at_extremity(PSouth) ) Then
              f3_at_p(:,pdims%j_start) = -two_omega
            End if
            If (at_extremity(PNorth) ) Then
              f3_at_p(:,pdims%j_end) = two_omega
            End if
          End if  
END IF ! vatpoles
          
        Else
          ALLOCATE(f3_at_p(1,1))
        End if


! DEPENDS ON: imps_intct
        CALL IMPS_INTCT(halo_i,halo_j,offx,offy,pdims%i_end,pdims%j_end,&
                        n_rows,global_row_length, proc_row_group,  &
     &                  at_extremity, n_proc, n_procx, n_procy,         &
     &                  neighbour,                                      &
     &                  ntiles, land_points, nice, nice_use,            &

! in values defining model domain, vertical grid of model atmosphere
!    and implicit weights. :
     &                  model_domain, bl_levels,                        &
     &                  alpha_cd,                                       &

! IN U and V momentum fields.
     &                  u, v,                                           &
! IN Non turbulent increments to momentum and thermodynamic variables.
!  (New dynamics only).
     &                  R_u, R_v,                                       &
! in soil/vegetation/land surface data :
     &                  land_sea_mask, land_index,                      &
     &                  dst_levels, dsm_levels,                         &

! in sea/sea-ice data :
     &                  di,ice_fract, di_ncat, ice_fract_ncat,          &
     &                  k_sice, u_0, v_0,                               &

! in cloud data :
! cloud_fraction passed in
     &                  q(qdims_l%i_start,qdims_l%j_start,1),           &
                        qcf(qdims_l%i_start,qdims_l%j_start,1),         &
                        qcl(qdims_l%i_start,qdims_l%j_start,1),         &
                        q_conv, qcf_conv, qcl_conv,                     &
     &                  qcf_latest, qcl_latest, T, T_conv,              &
!
! in everything not covered so far :
                        rho_wet_theta, p_star, timestep,                &
     &                  L_sice_heatflux,L_sice_multilayers,             &

! in variables required in IMP_SOLVER
     &                  alpha1_sice, ashtf, BQ_GB, BT_GB,               &
     &                  dtrdz_charney_grid, rdz_charney_grid,           &
     &                  dtrdz_u, dtrdz_v, rdz_u, rdz_v,                 &
     &                  cdr10m_u, cdr10m_v, z1_tq, zh,                  &
     &                  rhokm_u,rhokm_v,                                &
! in variables for new BL solver                                
     &                  bl_type_1,bl_type_2,bl_type_3,bl_type_4,        &
     &                  bl_type_5,bl_type_6,bl_type_7,                  &

! IN additional variables for MOSES II
     &                  TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY,           &
                        ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,            &
     &                  RESFS,Z0HSSI,Z0MSSI,CANHC_TILE,FLAKE,           &
     &                  WT_EXT_TILE,LW_DOWN,SW_TILE,ASHTF_TILE,         &
     &                  FQT_ICE,FTL_ICE,RESFT,RHOKH_SICE,               &
     &                  RHOKPM,RHOKPM_POT,RHOKPM_SICE,                  &
     &                  Z0H_TILE,Z0M_TILE,CHR1P5M_SICE,                 &
     &                  FLAND,                                          &
     &    FLANDG(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
     &                  FLANDG_U,FLANDG_V,TSTAR_SEA,                    &
! IN additional variables used to calculate atmospheric cooling
!    at the screen level
     &                  L_co2_interactive, co2_mmr, co2(:, :, 1),       &
     &                  rho1, f3_at_p, uStarGBM,                        &

! IN additional variables for JULES
     & RHOKH_MIX_DUMMY,DTSTAR_TILE,DTSTAR,HCONS,EMIS_TILE,EMIS_SOIL,    &

! INOUT data :
! (Note ti and ti_bg are IN only if l_sice_multilayers=T)
     &                  t_soil, ti, ti_gb, t_surf,                      &
     &                  T_latest, Q_latest,                             &

! inout  diagnostic started in bdy_layr not requiring stash flags :
     &                  e_sea, fqT, ftl, h_sea, rhokh,                  &
     &                  taux, tauy,                                     &

! INOUT additional variables for MOSES II
     &                  TSTAR_TILE,FQT_TILE,EPOT_TILE,FTL_TILE,         &
     &                  SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR,  &
     & TSTAR_SICE,TSTAR_SICE_CAT,TSTAR_SSI,                             &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &
     & TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                          &

! out u and v increments
     &                  u_inc_bl, v_inc_bl,                             &

! out  diagnostic not requiring stash flags :
     &                  rhokh_mix,                                      &
     &                  sea_ice_htf, surf_ht_flux_gb,                   &
     &                  SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,            &
     &                  SURF_HTF_TILE,                                  &

! out diagnostic requiring stash flags :
     &                  sice_mlt_htf, snomlt_surf_htf,                  &
     &                  latent_heat, q1p5m, t1p5m, u10m, v10m,          &

! (in) stash flags :-
     &              simlt, smlt, slh, sq1p5, st1p5, su10, sv10,         &
     &              l_ftl, l_fqw, l_taux, l_tauy,                       &

! OUT additional variables for MOSES II
     &                  ESOIL_TILE,ES,EI_TILE,                          &
     &                  Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE,      &
     &                  TSTAR_LAND,E_SSI,EI_SICE,FTL_SSI,               &
     &                  ERROR_CODE,                                     &

! out data required elsewhere in um system :
     &                  ecan, ei, ext, snowmelt,                        &
     &                  BL_diag, lq_mix_bl,                             &
     &                  L_flux_bc, Ltimer )

!       Release space used for the screen diagnostic.
        DEALLOCATE(f3_at_p)


! add boundary layer increment to R_u and R_v
        If ( l_u_incr_bl ) Then
! add boundary layer increment to R_u and R_v
        Do k = 1, bl_levels
          Do j = udims%j_start, udims%j_end
            Do i = udims%i_start, udims%i_end
              u_incr_diag_bl(i,j,k) = u_incr_diag_bl(i,j,k)             &
     &                              + ( u_inc_bl(i,j,k) - R_u(i,j,k) )
              R_u(i,j,k) = u_inc_bl(i,j,k)
            End Do
          End Do
        Enddo
        Else
        Do k = 1, bl_levels
          Do j = udims%j_start, udims%j_end
            Do i = udims%i_start, udims%i_end
              R_u(i,j,k) = u_inc_bl(i,j,k)
            End Do
          End Do
        Enddo
        Endif

        If ( l_v_incr_bl ) Then
        Do k = 1, bl_levels
          Do j = vdims%j_start, vdims%j_end
            Do i = vdims%i_start, vdims%i_end
              v_incr_diag_bl(i,j,k) = v_incr_diag_bl(i,j,k)             &
     &                              + ( v_inc_bl(i,j,k) - R_v(i,j,k) )
              R_v(i,j,k) = v_inc_bl(i,j,k)
            End Do
          End Do
        End Do
        Else
        Do k = 1, bl_levels
          Do j = vdims%j_start, vdims%j_end
            Do i = vdims%i_start, vdims%i_end
              R_v(i,j,k) = v_inc_bl(i,j,k)
            End Do
          End Do
        End Do
        Endif

! Call tr_mix to mix qcf
! output qcf_flux in T to save workspace
! Pass in a zero field for source terms.
        Do j = pdims%j_start,pdims%j_end
          Do i = pdims%i_start,pdims%i_end
            interp(i,j) = 0.
          End Do
        End Do

        CALL  tr_mix (                                                  &
! IN fields
     &        bl_levels, alpha_tr, rhokh_mix(:,:,2:), rhokh_mix(:,:,1)  &
     &       ,dtrdz_charney_grid,interp, interp                         &
     &       ,KENT, WE_LIM, T_FRAC, ZRZI                                &
     &       ,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                &
     &       ,ZH ,ZHSC, Z_HALF                                          &
! INOUT / OUT fields
     &       ,qcf_latest,T,drydep_str                                   &
     &        )
!
!Call tr_mix to mix w in the vertical when the subgrid turbulence
!scheme is activated
!
        If (L_subfilter_vert .OR. L_subfilter_horiz) then

          Do k = 1, bl_levels
            Do j = wdims%j_start,wdims%j_end
              Do i = wdims%i_start,wdims%i_end
                w_int(i,j,k) = w(i,j,k)
              End Do
            End Do
          End Do
!
!rhokm should be used in mixing of w.  rhokh_mix is used instead to avoid
!passing around rhokm on p-points.
!
          CALL  tr_mix (                                                &
! In fields
     &        bl_levels, alpha_cd,rhokh_mix(:,:,2:), rhokh_mix(:,:,1)   &
     &       ,dtrdz_charney_grid,interp, interp                         &
     &       ,KENT, WE_LIM, T_FRAC, ZRZI                                &
     &       ,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                &
     &       ,ZH ,ZHSC, Z_HALF                                          &
! INOUT / OUT fields
     &       ,w_int,T,drydep_str                                        &
     &        )

          Do k = 1, bl_levels
            Do j = wdims%j_start,wdims%j_end
              Do i = wdims%i_start,wdims%i_end
                R_w(i,j,k) = R_w(i,j,k) + (w_int(i,j,k) - w(i,j,k))
              End Do
            End Do
          End Do

        End If     !L_subfilter_vert or L_subfilter_horiz

! apply BL tracer mixing and gravitational settling of dust
! on the last cycle only
      If ( CycleNo == NumCycles ) Then
!
! Gravitational settling of mineral dust
!
        IF (L_DUST) THEN
         
          DO K = 1, tdims%k_end
            DO J = tdims%j_start,tdims%j_end
              DO I = tdims%i_start,tdims%i_end
                T_MODELLEVS(I,J,K)=THETA(I,J,K)*EXNER_THETA_LEVELS(I,J,K)
                DUST_ALL(I,J,K,1)=DUST_DIV1(I,J,K)
                DUST_ALL(I,J,K,2)=DUST_DIV2(I,J,K)
                IF (.NOT.l_twobin_dust) THEN 
                  DUST_ALL(I,J,K,3)=DUST_DIV3(I,J,K)
                  DUST_ALL(I,J,K,4)=DUST_DIV4(I,J,K)
                  DUST_ALL(I,J,K,5)=DUST_DIV5(I,J,K)
                  DUST_ALL(I,J,K,6)=DUST_DIV6(I,J,K)
                END IF
              ENDDO !i
            ENDDO !j
          ENDDO !k


          DO IDIV=1,NDIV
!
! DEPENDS ON: gravsett
            CALL GRAVSETT(                                              &
       DREP(IDIV),RHOP,P_LAYER_CENTRES,P_LAYER_BOUNDARIES,T_MODELLEVS,  &
       DUST_ALL(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                pdims%k_start:pdims%k_end,IDIV),                        & 
       DRYDEP2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
       IDIV))
!
          ENDDO
!
!
          DO K = 1, tdims%k_end
            DO J = tdims%j_start,tdims%j_end
              DO I = tdims%i_start,tdims%i_end
                DUST_DIV1(I,J,K)=DUST_ALL(I,J,K,1)
                DUST_DIV2(I,J,K)=DUST_ALL(I,J,K,2)
                IF(.NOT.l_twobin_dust) THEN
                  DUST_DIV3(I,J,K)=DUST_ALL(I,J,K,3)
                  DUST_DIV4(I,J,K)=DUST_ALL(I,J,K,4)
                  DUST_DIV5(I,J,K)=DUST_ALL(I,J,K,5)
                  DUST_DIV6(I,J,K)=DUST_ALL(I,J,K,6)
                END IF
              ENDDO !ROW_LENGTH
            ENDDO !ROWS
          ENDDO !MODEL_LEVELS
        ENDIF !L_DUST

! Mixing for all non-qcf tracers done in subroutine
        If ( l_bl_tracer_mix .OR. l_sulpc_so2 .OR. l_soot .OR.          &
     &       l_co2_interactive .OR. l_murk .OR. l_biomass .OR.          &
     &       l_dust .OR. l_ocff .OR. L_nitrate) Then

! DEPENDS ON: bl_trmix_dd
          CALL BL_TRMIX_DD(                                             &
! IN arguments
     &         bl_levels,dtrdz_charney_grid,tr_vars, DIM_CS2            &
     &        ,alpha_tr, rhokh_mix, p_star                              &
! IN Control logicals
     &        ,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST                     &
     &        ,L_SULPC_SO2, l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass &
     &        ,l_ocff, l_nitrate, l_co2_interactive                     &
     &        ,l_co2_emits, l_use_cariolle                              &
! IN Emissions fields
     &        ,DUST_FLUX, co2_emits, co2flux, npp, resp_s_tot           &
! IN variables for tr_mix
     &        ,kent, we_lim, t_frac, zrzi                               &
     &        ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc               &
     &        ,zh, zhsc, z_half                                         &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist                                      &
     &,        R_B_DUST,T_SURF, land_points, land_index, ice_fract      &
! IN variables for sresfact
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile                               &
     &,   FLANDG(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
! INOUT Fields to mix
               aerosol(tdims_s%i_start,tdims_s%j_start,1),              &
               free_tracers,                                            &
               ozone_tracer(tdims_s%i_start,tdims_s%j_start,1),         &
               drydep_str, resist_b,                                    &
! Mineral dust
               DUST_DIV1(tdims_s%i_start,tdims_s%j_start,1),            &
               DUST_DIV2(tdims_s%i_start,tdims_s%j_start,1),            &
               DUST_DIV3(tdims_s%i_start,tdims_s%j_start,1),            &
               DUST_DIV4(tdims_s%i_start,tdims_s%j_start,1),            &
               DUST_DIV5(tdims_s%i_start,tdims_s%j_start,1),            &
               DUST_DIV6(tdims_s%i_start,tdims_s%j_start,1),            &
! Sulphur cycle
               so2(tdims_s%i_start,tdims_s%j_start,1),                  &
               dms(tdims_s%i_start,tdims_s%j_start,1),                  &
               so4_aitken(tdims_s%i_start,tdims_s%j_start,1),           &
               so4_accu(tdims_s%i_start,tdims_s%j_start,1),             &
               so4_diss(tdims_s%i_start,tdims_s%j_start,1),             &
               nh3(tdims_s%i_start,tdims_s%j_start,1),                  &
! Soot cycle
               soot_new(tdims_s%i_start,tdims_s%j_start,1),             &
               soot_aged(tdims_s%i_start,tdims_s%j_start,1),            &
               soot_cld(tdims_s%i_start,tdims_s%j_start,1),             &
! Biomass aerosol
               bmass_new(tdims_s%i_start,tdims_s%j_start,1),            &
               bmass_agd(tdims_s%i_start,tdims_s%j_start,1),            &
               bmass_cld(tdims_s%i_start,tdims_s%j_start,1),            &
! Fossil-fuel organic carbon aerosol
               ocff_new(tdims_s%i_start,tdims_s%j_start,1),             &
               ocff_aged(tdims_s%i_start,tdims_s%j_start,1),            &
               ocff_cld(tdims_s%i_start,tdims_s%j_start,1),             &
! Ammonium nitrate aerosol
               nitr_acc(tdims_s%i_start,tdims_s%j_start,1),             &
               nitr_diss(tdims_s%i_start,tdims_s%j_start,1),            &
! Carbon cycle
               co2(tdims_s%i_start,tdims_s%j_start,1),                  &
               co2_flux_tot, land_co2,                                  &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &        STASHwork3,                                               &
! SCM diagnostics (dummy in full UM)
     &        nSCMDpkgs, L_SCMDiags                                     &
     &        )

        End If  ! if tracer mixing required

      End If ! ( CycleNo == NumCycles )
! ----------------------------------------------------------------------
! Section BL.3 Form Ch terms.
! ----------------------------------------------------------------------

! NB: Code assumes that the number of boundary layer levels is at
!     least 2.
        If ( L_use_bl_diag_term ) Then
! Remove factor of r**2 from rho, store in interp

        Do k = 1, bl_levels

          Do j = pdims%j_start,pdims%j_end
            Do i = pdims%i_start,pdims%i_end
              interp(i,j) = ( ( r_theta_levels(i,j,k) -                 &
     &                          r_rho_levels(i,j,k) )                   &
     &                         * rho(i,j,k+1) /                         &
     &                          (r_rho_levels(i,j,k+1) *                &
     &                           r_rho_levels(i,j,k+1))                 &
     &                       + ( r_rho_levels(i,j,k+1) -                &
     &                           r_theta_levels(i,j,k) )                &
     &                         * rho(i,j,k) /                           &
     &                          (r_rho_levels(i,j,k) *                  &
     &                           r_rho_levels(i,j,k))                   &
     &                      ) / (r_rho_levels(i,j,k+1) -                &
     &                           r_rho_levels(i,j,k) )
            End Do
          End Do

          If (k  ==  1) Then
            Do j = pdims%j_start,pdims%j_end
              Do i = pdims%i_start,pdims%i_end

                temp1 = rhokh(i,j,k)

                temp2 = rhokh(i,j,k+1)/ (r_theta_levels(i,j,k+1)        &
     &                                   - r_theta_levels(i,j,k))

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                           / (interp(i,j) *                       &
     &                           (r_rho_levels(i,j,k+1) -               &
     &                            r_rho_levels(i,j,k)))

              End Do
            End Do


          Else If (k  <   bl_levels) Then
            Do j = pdims%j_start,pdims%j_end
              Do i = pdims%i_start,pdims%i_end

                temp1 = rhokh(i,j,k)/ (r_theta_levels(i,j,k)            &
     &                                 - r_theta_levels(i,j,k-1))

                temp2 = rhokh(i,j,k+1)/ (r_theta_levels(i,j,k+1)        &
     &                                   - r_theta_levels(i,j,k))

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                         / (interp(i,j) * (r_rho_levels(i,j,k+1)  &
     &                            - r_rho_levels(i,j,k)))

              End Do
            End Do

          Else ! top boundary layer level

            Do j = pdims%j_start,pdims%j_end
              Do i = pdims%i_start,pdims%i_end

                temp1 = rhokh(i,j,k)/ (r_theta_levels(i,j,k)            &
     &                                 - r_theta_levels(i,j,k-1))

                temp2 = 0.0

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                         / (interp(i,j) * (r_rho_levels(i,j,k+1)  &
     &                             - r_rho_levels(i,j,k)))

              End Do
            End Do

          End If

        End Do
! swop halo values for Ch_term in boundary layer
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   cH_term, pdims%i_end,pdims%j_end,              &
     &                   bl_levels, offx , offy , fld_type_p,  .false.)

! set zero value above Boundary layer
        Do k = bl_levels + 1, pdims_s%k_end - 1
          Do j = pdims_s%j_start, pdims_s%j_end
            Do i = pdims_s%i_start, pdims_s%i_end
              Ch_term(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Endif ! If  L_use_bl_diag_term


! ----------------------------------------------------------------------
! Section BL.4 Convert and calculate theta and q fields from qT and Tl.
! ----------------------------------------------------------------------

        If (.not. L_dry) Then
! If the mixed phase precipitation scheme is used then T and Q are
! required to contain T liquid and Q(vapour+liquid) but at this stage
! will actually contain T liquid ice and Q(vapour+liquid+ice)

! DEPENDS ON: bl_lsp
          CALL bl_lsp( bl_levels, qcf_latest, q_latest, t_latest )

        End If

        If (L_dry) Then
          Do k = bl_levels+1, qdims%k_end
            Do j = qdims%j_start,qdims%j_end
              Do i = qdims%i_start,qdims%i_end
                q_latest(i,j,k) = 0.0
                qcl_latest(i,j,k) = 0.0
                qcf_latest(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If

! Create Tl and qT outside boundary layer
        Do k = bl_levels+1, qdims%k_end
          Do j = qdims%j_start,qdims%j_end
            Do i = qdims%i_start, qdims%i_end
              T_latest(i,j,k) = t_latest(i,j,k) -                       &
     &                        (lc * qcl_latest(i,j,k))                  &
     &                         / cp
              q_latest(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)
            End Do
          End Do
        End Do
!
! Prepare for cloud scheme. Are we using PC2 or not?
!
! zero any negative q_latests
        Do k = 1, qdims%k_end
          Do j = qdims%j_start,qdims%j_end
            Do i = qdims%i_start, qdims%i_end
              If (q_latest(i,j,k)  <   0.) Then
!              print*,' neg qT before BL cld call set to zero',i,j,k
                q_latest(i,j,k) = 0.
              End If
              If (qcf_latest(i,j,k)  <   0.) Then
!              print*,' neg qcf before BL cld call set to zero',i,j,k
                qcf_latest(i,j,k) = 0.
              End If
            End Do
          End Do
        End Do
!
! Calculate diagnostic RHcrit or read as parameter in from namelist
!
! Lrhcpt_if1:
        If (L_RHCPT) Then
!       RHCRIT is 3D diagnosed variable
! Wet_mlev_do1:
          Do k = 1, qdims%k_end
! Rhc_rows_do1:
            Do j = 1, rhc_rows
! Rhc_rowlen_do1:
              Do i = 1, rhc_row_length
                ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)
                ext_TL(i,j,k) = T_latest(i,j,k)
                ext_QL(i,j,k) = q_latest(i,j,k)
                ext_QCF(i,j,k) = qcf_latest(i,j,k)
              End Do ! Rhc_rowlen_do1
            End Do ! Rhc_rows_do1
          End Do ! Wet_mlev_do1
!
! Rhc_rows_do2:
          Do j = 1, rhc_rows
! Rhc_rowlen_do2:
            Do i = 1, rhc_row_length
              ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
              ext_ICE_FRAC(i,j) = ice_fract(i,j)
! extended halo land fraction now passed in from AP2
            End Do ! Rhc_rowlen_do2
          End Do ! Rhc_rows_do2
!
! Synchronize haloes.
!
          i_field = 0
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ext_p_layer_centres(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  qdims%k_end+1
          fields_to_swap(i_field) % rows        =  rhc_rows
          fields_to_swap(i_field) % vector      =  .FALSE.
  
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ext_tl(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  qdims%k_end
          fields_to_swap(i_field) % rows        =  rhc_rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ext_ql(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  qdims%k_end
          fields_to_swap(i_field) % rows        =  rhc_rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ext_qcf(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  qdims%k_end
          fields_to_swap(i_field) % rows        =  rhc_rows
          fields_to_swap(i_field) % vector      =  .FALSE.

! DEPENDS ON: swap_bounds_mv
          CALL swap_bounds_mv( fields_to_swap, i_field,                 &
                               rhc_row_length, 1, 1)

          i_field = 0

          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d    => ext_ice_frac(:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  1
          fields_to_swap(i_field) % rows        =  rhc_rows
          fields_to_swap(i_field) % vector      =  .FALSE.

! DEPENDS ON: swap_bounds_2d_mv
          CALL swap_bounds_2d_mv( fields_to_swap, i_field,              &
                               rhc_row_length, 1, 1)

!
! DEPENDS ON: ls_calc_rhcrit
          CALL ls_calc_rhcrit( ext_p_layer_centres                      &
!              Array dimensions
     &     , qdims%k_end, rhc_row_length, rhc_rows                      &
     &     , global_row_length                                          &
!              Prognostic Fields
     &     , ext_TL, ext_QL, ext_QCF                                    &
     &     , FLANDG(0:rhc_row_length,0:rhc_rows), ext_ICE_FRAC          &
!              Logical control
     &     , lq_mix_bl                                                  &
                           ! Use mixing_ratio
!              Output
     &     , RHCPT, delta_lambda)
!
        Else
!         RHCRIT is 1D Parameter read in from namelist
          Do k = 1, qdims%k_end
            RHCPT(1,1,k) = rhcrit(k)
          End Do
        End if  ! Lrhcpt_if1
!
! Which cloud scheme are we using?
!
        If (L_pc2) then
!
! ----------------------------------------------------------------------
! PC2 cloud scheme
! ----------------------------------------------------------------------
          allocate ( zeros   (qdims%i_start:qdims%i_end,                &
                              qdims%j_start:qdims%j_end,qdims%k_end) )
!
! ----------------------------------------------------------------------
! Inhomogenous forcing of ice
! ----------------------------------------------------------------------
!
! Calculate in-plume ice content (LS) by assuming that they are equal to
! the current values, except when the current value is not defined.
! Also calculate the forcing of ice content Q4F
!
          Do k = 1, qdims%k_end
            Do j = qdims%j_start,qdims%j_end
              Do i = qdims%i_start,qdims%i_end
!
! Calculate Q4. Only perform the calculation if the Q4 is non-zero.
! Since ice is mixed by tracer mixing regardless of whether cumulus
! is present we still need to provide an ice cloud fraction increment
! below cumulus base.
!
              IF (l_fixbug_pc2_mixph) THEN
!
! Calculate change in ice cloud fraction differently depending on
! whether QCF has increased or decreased.
!
                q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k) 

                IF (q4 > 1.0e-12) THEN 
!
! Source. 
!
! Calculate the change in total cloud fraction. 
! Use a weighted (by ice cloud fraction) average of in-cloud ice 
! content and a fixed value to specify the plume ice content. 
! The denominator in the deltaCf calculation is then 
! ((qcf_earliest/cff_earliest)*cffearliest + 
! ls_bl0*(1-cff_earliest) - qcf_earliest 
! and the qcf_earliest terms cancel. 
!
                denom = ls_bl0 * ( 1.0 - cff_earliest(i,j,k) )

                IF ( ABS(denom) > 1.0e-10 ) THEN  

                   denom = q4 / denom 

                   deltacff = (1.0-cff_latest(i,j,k)) * denom  
                   cff_latest(i,j,k) = cff_latest(i,j,k) + deltacff  

! calculate total cf based on minimum overlap
                   cf_latest(i,j,k)  = cf_latest(i,j,k)  + deltacff  

                 ELSE  
                   cf_latest(i,j,k)  = 1.0  
                   cff_latest(i,j,k) = 1.0  
                 END IF

                ELSE IF (q4 < -1.0e-12) THEN 
! 
! Sink.
! 
! Given a reduction in qcf, remove some CFF in order to maintain the  
! same in-cloud IWC.   
! 

                  InCloudIce = qcf_latest(i,j,k) /                      &
                    MAX(cff_latest(i,j,k),1.e-6) 

                  cf_latest(i,j,k) = cf_latest(i,j,k) +                 & 
                    q4/MAX(InCloudIce, ls_bl0) 

                  cff_latest(i,j,k) = cff_latest(i,j,k) +               & 
                    q4/MAX(InCloudIce, ls_bl0) 

                  temp3=cff_latest(i,j,k)

                  IF (cff_latest(i,j,k) > 0.0) THEN 
                    !Prevent very high in-cloud values by increasing CFF.
                    IF ( (qcf_latest(i,j,k) / cff_latest(i,j,k) )       &
                         > 2.0e-3) THEN
                       cff_latest(i,j,k)=qcf_latest(i,j,k)/2.e3
                       temp3=cff_latest(i,j,k)-temp3
                      ! Update total cloud fraction.
                       cf_latest(i,j,k)=cf_latest(i,j,k)+temp3
                    END IF

                  END IF

                END IF

              ELSE
!
! Original code as run in Wilson et al (2008a,b).
! Calculation of change in cloud fraction is consitent with an increase
! in QCF but not consistent with a decrease.
!
                Q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k)
                IF (Q4  /=  0.0) Then
!
! Calculate the change in total cloud fraction.
! Use a weighted (by ice cloud fraction) average of in-cloud ice
! content and a fixed value to specify the plume ice content.
! The denominator in the deltaCf calculation is then
! (qcf_earliest/cff_earliest + ls_bl0*(1-cff_earliest) - qcf_earliest
! and the qcf_earliest terms cancel.
!
                  DENOM = LS_bl0 * (1.0 - cff_earliest(i,j,k))
!
                  IF ( ABS(DENOM)  >   1.0E-10 ) THEN
!
                    DENOM = Q4 / DENOM
                    cf_latest(i,j,k)  = cf_latest(i,j,k)  +             &
     &                           (1.0 - cf_latest(i,j,k))  * DENOM
                    cff_latest(i,j,k) = cff_latest(i,j,k) +             &
     &                           (1.0 - cff_latest(i,j,k)) * DENOM
!
! Otherwise cloud fraction will go to one. In theory, cloud fraction
! can never be reduced by this process.
!
                  ELSE
                    cf_latest(i,j,k)  = 1.0
                    cff_latest(i,j,k) = 1.0
                  END IF
!
                END IF !(Q4  /=  0.0)

                END IF ! l_fixbug_pc2_mixph
!
! Set zero array
                zeros(i,j,k)        = 0.0
              end do
            end do
          end do
!
! ----------------------------------------------------------------------
! Homogenous forcing of the liquid cloud
! ----------------------------------------------------------------------
!
! Calculate forcing in qT and TL. Currently q_latest contains the vapour
! plus liquid content and q_earliest just the initial vapour content
!
          allocate ( QT_force(qdims%i_start:qdims%i_end,                &
                              qdims%j_start:qdims%j_end,qdims%k_end) )
          allocate ( TL_force(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,tdims%k_end) )
          allocate ( T_inc_PC2(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,tdims%k_end) )
          allocate ( q_inc_PC2(qdims%i_start:qdims%i_end,               &
                               qdims%j_start:qdims%j_end,qdims%k_end) )
          allocate ( qcl_inc_PC2(qdims%i_start:qdims%i_end,             &
                                 qdims%j_start:qdims%j_end,qdims%k_end) )
          allocate ( cfl_inc_PC2(qdims%i_start:qdims%i_end,             &
                                 qdims%j_start:qdims%j_end,qdims%k_end) )
          allocate ( bcf_inc_PC2(qdims%i_start:qdims%i_end,             &
                                 qdims%j_start:qdims%j_end,qdims%k_end) )
!
          Do k = 1, qdims%k_end
            Do j = qdims%j_start,qdims%j_end
              Do i = qdims%i_start,qdims%i_end
                qT_force(i,j,k) = ( q_latest(i,j,k)                     &
     &            - (q_earliest(i,j,k) + qcl_earliest(i,j,k)) )
                TL_force(i,j,k) = ( T_latest(i,j,k)                     &
     &            - (t_earliest(i,j,k)- lc * qcl_earliest(i,j,k) / cp) )
              End Do
            End Do
          End Do
!
! Call homogenous forcing routine
!
! DEPENDS ON: pc2_delta_hom_turb
          CALL PC2_DELTA_HOM_TURB(                                      &
! INput variables
     &      p_layer_centres(1,1,1)                                      &
     & ,    timestep                                                    &
! INput variables
     & ,    T_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)   &
     & ,    cf_latest(1,1,1), cfl_latest(1,1,1), cff_latest(1,1,1)      &
     & ,    TL_force(1,1,1),qT_force(1,1,1),zeros(1,1,1),zeros(1,1,1)   &
! OUTput variables
     & ,    T_inc_PC2, q_inc_PC2, qcl_inc_PC2, bcf_inc_PC2, cfl_inc_PC2 &
! INput variables (other quantities)
     & ,    0.0, 0.0, lq_mix_bl)
!
! Diagnostic shallow cloud within PC2
!
          IF ( l_pc2_diag_sh .AND. .NOT. l_ccrad .AND.              &
               Rad_cloud_decay_opt == rad_decay_off ) THEN
            !
            ! The following code is only valid with the follow settings.
            ! 
            ! l_pc2               = .true.
            ! l_pc2_diag_sh     = .true.
            ! l_ccrad           = .false.
            ! Rad_cloud_decay_opt = rad_decay_off

            ALLOCATE ( ccw_cca(qdims%i_start:qdims%i_end,               &
                               qdims%j_start:qdims%j_end,qdims%k_end) )
            ALLOCATE ( cca_3d (qdims%i_start:qdims%i_end,               &
                               qdims%j_start:qdims%j_end,qdims%k_end) )

            Do k=1, qdims%k_end
              Do j=qdims%j_start,qdims%j_end
                Do i=qdims%i_start,qdims%i_end
!
!
! If below cumulus cloud base or we are using a diagnostic
! shallow convective cloud then we simply zero the liquid water
! content instead of using the homogenous BL response. The forcings
! themselves are still applied but to q and T.
! Note that we find, for PC2, that zeroing up to ntml,
! rather than ntml+1, is more physically justified and gives slightly
! better cloud results.
!
                  If ( (l_pc2_diag_sh_pts(i,j) .AND. k <= cct0(i,j)) .OR.     &
                       (cumulus(i,j)           .AND. k <= ntml(i,j)) ) Then
                    T_inc_PC2(i,j,k)   = (-LC/CP) * qcl_earliest(i,j,k)
                    q_inc_PC2(i,j,k)   = qcl_earliest(i,j,k)
                    qcl_inc_PC2(i,j,k) = (-qcl_earliest(i,j,k))
                    cfl_inc_PC2(i,j,k) = (-cfl_earliest(i,j,k))
                    bcf_inc_PC2(i,j,k) = cff_latest(i,j,k)-cf_earliest(i,j,k)
                  End If

                  ! Set convective cloud properties if we are
                  ! using a diagnostic shallow convection for PC2
                  ! Calculate convective properties

                  ! With l_PC2 =true, the cca0 which goes
                  ! to radiation should only contain only
                  ! shallow convection (i.e. no deep or mid-level).
 
                  ! l_pc2_diag_sh and l_ccrad should not be
                  ! true at the same time.
                  IF (l_3d_cca) THEN 

                    ccw_cca(i,j,k) = ccw0(i,j,k)*cca0(i,j,k)
                    cca_3d(i,j,k)  = cca0(i,j,k)

                  ELSE

                    ! Currently with this switch cca is set to be
                    ! a single level field. In which case use cca0_2d
                    ! which should also only contain cca0_2d from
                    ! shallow cloud

                    IF (k <= cct0(i,j) .AND. k >= ccb0(i,j)) THEN
                      ccw_cca(i,j,k) = ccw0(i,j,k) * cca0_2d(i,j)
                      cca_3d(i,j,k)  = cca0_2d(i,j)
                    ELSE
                      ccw_cca(i,j,k) = 0.0
                      cca_3d(i,j,k)  = 0.0
                    END IF
                  END IF


                  If ( l_pc2_diag_sh_pts(i,j) .AND. k <= cct0(i,j)    &
                                .AND. t_earliest(i,j,k) >= tice) Then

                    ! Hand over convective attributes to the large scale 
                    T_inc_PC2(i,j,k)   =  T_inc_PC2(i,j,k)            &
                                        + (LC/CP) *  ccw_cca(i,j,k)
                    q_inc_PC2(i,j,k)   =  q_inc_PC2(i,j,k) -  ccw_cca(i,j,k)  
                    qcl_inc_PC2(i,j,k) =  qcl_inc_PC2(i,j,k)          &
                                        + ccw_cca(i,j,k)
                    cfl_inc_PC2(i,j,k) =  cfl_inc_PC2(i,j,k)          &
                                         + cca_3d(i,j,k)
                    bcf_inc_PC2(i,j,k) =  bcf_inc_PC2(i,j,k)          &
                                    + cca_3d(i,j,k) *(1.0-cff_latest(i,j,k))

                  End if

                  ! Reset the section 0 convective cloud back to zero
                  ! now that increments have been included in PC2
                  ! increments.
                  cca0(i,j,k)  = 0.0 
                  cca0_2d(i,j) = 0.0
                  ccw0(i,j,k)  = 0.0

                End Do  ! i loop
              End Do  ! j
            End Do  ! k

            deallocate (ccw_cca)
            deallocate (cca_3d)

          Else     ! original code

            Do k=1, qdims%k_end
              Do j=qdims%j_start,qdims%j_end
                Do i=qdims%i_start,qdims%i_end
!
!
! If below the LCL in cumulus or unstable (but not shear-dominated) 
! boundary layers then we simply zero the liquid water content 
! instead of using the homogenous BL response. The forcings
! themselves are still applied but to q and T.
                  z_theta = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
                  IF ( ( cumulus(i,j) .AND. k  <=  ntml(i,j) ) .OR.      &
                       ( forced_cu >= 1 .AND. (bl_type_3(i,j) > 0.5      &
                          .OR. bl_type_4(i,j) > 0.5 )                    &
                          .AND. z_theta  <  zlcl(i,j) )  ) THEN
                    T_inc_PC2(i,j,k)   =  (-LC/CP) * qcl_earliest(i,j,k)
                    q_inc_PC2(i,j,k)   =  qcl_earliest(i,j,k)
                    qcl_inc_PC2(i,j,k) =  (-qcl_earliest(i,j,k))
                    cfl_inc_PC2(i,j,k) =  (-cfl_earliest(i,j,k))
                    bcf_inc_PC2(i,j,k) =  cff_latest(i,j,k)              &
                                            -cf_earliest(i,j,k)
                  END IF
                End Do  ! i loop
              End Do  ! j
            End Do  ! k

          END IF    ! test on l_pc2_diag_sh

          DO k = 1, qdims%k_end
            DO j = qdims%j_start,qdims%j_end
              DO i = qdims%i_start,qdims%i_end
!
! Update working version of temperature, moisture and cloud fields with
! increments from the PC2 homogenous response.
!
                T_latest(i,j,k)   = T_earliest(i,j,k) + TL_force(i,j,k) &
     &                             + T_inc_PC2(i,j,k)
                q_latest(i,j,k)   = q_earliest(i,j,k) + qT_force(i,j,k) &
     &                             + q_inc_PC2(i,j,k)
                qcl_latest(i,j,k) = qcl_earliest(i,j,k)                 &
     &                             + qcl_inc_PC2(i,j,k)
                cfl_latest(i,j,k) = cfl_earliest(i,j,k)                 &
     &                             + cfl_inc_PC2(i,j,k)
                cf_latest(i,j,k)  =  cf_earliest(i,j,k)                 &
     &                             + bcf_inc_PC2(i,j,k)
!               qcf_latest and cff_latest are not updated
!
              END DO  ! i loop
            END DO  ! j
          END DO  ! k
!
          DEALLOCATE (zeros)
!------------------------------------------------------------
! Parametrize "forced cumulus clouds" at top of well-mixed BL
!------------------------------------------------------------
        IF (forced_cu == 1) THEN
! Use TL_force to store TL_latest
          DO k = 1, qdims%k_end
            DO j = tdims%j_start,tdims%j_end
              DO i = tdims%i_start,tdims%i_end
                TL_force(i,j,k) = T_latest(i,j,k)                       &
                                  - (lc/cp) * qcl_latest(i,j,k)
              END DO  ! i loop
            END DO  ! j
          END DO  ! k

          DO k = 1, qdims%k_end
! DEPENDS ON: qsat_mix
            CALL qsat_mix(qs(1,1,k),TL_force(1,1,k),                    &
              p_layer_centres(1,1,k),tdims%i_end*tdims%j_end,lq_mix_bl)
          END DO

          DO j = qdims%j_start,qdims%j_end
            DO i = qdims%i_start,qdims%i_end
              DO k = 1, qdims%k_end
                z_theta = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
                zc_depth = zh(i,j)+dzh(i,j)-zlcl(i,j)
                IF ( z_theta >= zlcl(i,j)                               &   
                    .AND. z_theta <= zh(i,j)+dzh(i,j)                   &
                    .AND. zc_depth > 1.0                                &
                    .AND. cfl_latest(i,j,k) < 0.5                       &
                    .AND. bl_type_3(i,j)>0.5                            &
                   ) THEN
                  rht = (q_latest(i,j,k)+qcl_latest(i,j,k))/qs(i,j,k)
                  cf_top  = 0.1
                  ! Make cloud fraction at cloud base a function of 
                  ! the cloud depth:
                  cf_base = MAX( cf_top,                                &
                                 0.3 * MIN( 1.0, zc_depth/300.0 ) )
                  cf_forced = cf_base - (cf_base-cf_top) *              &
                                        (z_theta-zlcl(i,j)) / zc_depth 
                  ! calculate implied RHcrit (between 0.01 and 0.99)
                  rhc = MAX( 0.01, MIN( 0.99,                           &
                                   (rht - SQRT(2.0*cf_forced))/         &
                                   (1.0- SQRT(2.0*cf_forced)) ) )
                  qn = (rht-1.0)/(1.0-rhc)
                  alpha=repsilon*lc*qs(i,j,k)/(r*TL_force(i,j,k)**2)
                  al=1.0/(1.0+lcrcp*alpha)
                  bs=al*(1.0-rhc)*qs(i,j,k)
                  qcl_forced = 0.00001  ! small default cloud water
                  IF (qn > -1.0 .AND. qn < 0.0) THEN
                    qcl_forced = (bs/6.0)*(1.0 + qn)**3
                  ELSE IF (qn < 1.0) THEN
                    qcl_forced = bs*( qn + ((1.0-qn)**3)/6.0 )
                  ELSE
                    qcl_forced = bs*qn
                  END IF
                  IF ( cf_forced > cfl_latest(i,j,k) ) THEN
                    dcfl = cf_forced - cfl_latest(i,j,k)
                    cfl_latest(i,j,k) = cf_forced
                    cf_latest(i,j,k)  = cf_latest(i,j,k) + dcfl
                    dqcl = MAX( qcl_latest(i,j,k), qcl_forced )         &
                            - qcl_latest(i,j,k)
                    qcl_latest(i,j,k) = MAX( qcl_latest(i,j,k),         &
                                             qcl_forced )
                    T_latest(i,j,k)  = T_latest(i,j,k)+(LC/CP) * dqcl
                    q_latest(i,j,k)  = q_latest(i,j,k) - dqcl
                  END IF
                END IF  ! test on z and BLtype
              END DO
            END DO
          END DO
        END IF  ! test on forced_cu eq 1
        DEALLOCATE (QT_force)
        DEALLOCATE (TL_force)
! ----------------------------------------------------------------------
! Copy updated cloud fractions to the in/out variables
! ----------------------------------------------------------------------
!
          If (.not. L_area_cloud) Then

            Do k = 1, qdims%k_end
              Do j = qdims%j_start,qdims%j_end
                Do i = qdims%i_start,qdims%i_end
                  ! For the moment set area cloud fraction
                  ! to the bulk cloud fraction
                  area_cloud_fraction(i,j,k)=cf_latest(i,j,k)
                  ! Ensure it has a value between 0.0 and 1.0  
                  area_cloud_fraction(i,j,k)=&
                       max(min(area_cloud_fraction(i,j,k),1.0),0.0) 
                End Do
              End Do
            End Do

          Else If (L_area_cloud) Then

            If (L_ACF_Brooks) Then

! DEPENDS ON: ls_acf_brooks
              CALL LS_ACF_Brooks (                                      &
                   delta_lambda, delta_phi                              &
                  ,FV_cos_theta_latitude                                &
                  ,bulk_cloud_fraction, cloud_fraction_liquid           &
                  ,cloud_fraction_frozen, cumulus                       &
                  ,area_cloud_fraction )

            End If ! L_ACF_Brooks
!
! call sub-level interpolation parameterisation of cloud area
            If (L_ACF_Cusack) then
!
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
              levels_per_level = 3
              large_levels = ((qdims%k_end - 2)*levels_per_level) + 2
!
! DEPENDS ON: pc2_hom_arcld
              CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,    &
     &         large_levels,levels_per_level,                           &
     &         area_cloud_fraction,T_latest,bulk_cloud_fraction,        &
     &         cloud_fraction_liquid,cloud_fraction_frozen,q_latest,    &
     &         qcl_latest,qcf_latest,                                   &
     &         lq_mix_bl)
            End if ! L_ACF_Cusack
!
          End If ! L_area_cloud
!
! ----------------------------------------------------------------------
! Provide an estimate of convective cloud fraction for visibility
! ----------------------------------------------------------------------
!
          Do j=qdims%j_start,qdims%j_end
            Do i=qdims%i_start,qdims%i_end
              If (cca_at_base(i,j) == 0.0 .AND.                         &
                 (conv_rain(i,j) + conv_snow(i,j)) > 0.0) then
                 ! Convective precipitation exists but no 
                 ! estimate for its cloud fraction. Set it
                 ! to a constant value of 0.2 as an estimate.
                 cca_at_base(i,j) = 0.2
              End If
            End Do
          End Do            

!
! ----------------------------------------------------------------------
!  PC2: End of cloud section
! ----------------------------------------------------------------------
!
        Else   !  L_pc2
!
! ----------------------------------------------------------------------
! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
! in boundary layer, calculate bulk_cloud fields from qT and qcf
! and calculate area_cloud fields.
! ----------------------------------------------------------------------
!
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
          levels_per_level = 3
          large_levels = ((qdims%k_end - 2)*levels_per_level) + 2
!
! DEPENDS ON: ls_arcld
          CALL ls_arcld( p_layer_centres, RHCPT, p_layer_boundaries,    &
     &                 rhc_row_length, rhc_rows, bl_levels,             &
     &                 levels_per_level, large_levels,                  &
     &                 L_area_cloud,L_ACF_Cusack,L_ACF_Brooks,          &
     &                 delta_lambda, delta_phi,                         &
     &                 FV_cos_theta_latitude,                           &
     &                 ntml, cumulus, lq_mix_bl, qcf_latest,            &
     &                 T_latest, q_latest, qcl_latest,                  &
     &                 area_cloud_fraction, bulk_cloud_fraction,        &
     &                 cloud_fraction_liquid, cloud_fraction_frozen,    &
     &                 error_code, me )
!
        End If   ! L_pc2

        If ( L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc ) Then
          Do k=1,tdims%k_end
            Do j=tdims%j_start, tdims%j_end
              Do i=tdims%i_start, tdims%i_end
                T_incr_diag_bl(i,j,k) = T_incr_diag_bl(i,j,k)           &
     &            + ( T_latest(i,j,k) - T_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                   ! on STASHflags

        If ( L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                q_incr_diag_bl(i,j,k) = q_incr_diag_bl(i,j,k)           &
     &            + ( q_latest(i,j,k) - q_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflags

        If ( L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR.               &
     &       L_qtl_incr_bl_lsc ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                qcl_incr_diag_bl(i,j,k) = qcl_incr_diag_bl(i,j,k)       &
     &            + ( qcl_latest(i,j,k) - qcl_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflag

        If ( L_qcf_incr_bl_lsc ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                qcf_incr_diag_bl(i,j,k) = qcf_incr_diag_bl(i,j,k)       &
     &            + ( qcf_latest(i,j,k) - qcf_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflag!

        If ( L_cf_incr_bl ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                cf_incr_diag_bl(i,j,k) = cf_incr_diag_bl(i,j,k)         &
     &            + ( cf_latest(i,j,k) - cf_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif

        If ( L_cfl_incr_bl ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                cfl_incr_diag_bl(i,j,k) = cfl_incr_diag_bl(i,j,k)       &
     &            + ( cfl_latest(i,j,k) - cfl_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif

        If ( L_cff_incr_bl ) Then
          Do k=1,qdims%k_end
            Do j=qdims%j_start, qdims%j_end
              Do i=qdims%i_start, qdims%i_end
                cff_incr_diag_bl(i,j,k) = cff_incr_diag_bl(i,j,k)       &
     &            + ( cff_latest(i,j,k) - cff_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif
!
      End If ! on error code zero

!
! ----------------------------------------------------------------------
! Section BL.4c Combined cloud field calculation for use by visibility
!               (section 3) and cloud scheme (section 9) diagnostics
!               09208 - 09217 and 09223.
! ----------------------------------------------------------------------
!
! L_combi_cld_if1:
      If (error_code  <=  0  .AND.  L_combi_cld .AND. L_apply_diag)  Then
        ! Set the combined cloud area fractions in each gridbox.
        ! Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
      
        ALLOCATE                                            &
          ( combined_cloud (qdims%i_start:qdims%i_end,      &
                            qdims%j_start:qdims%j_end,      &
                            qdims%k_end)                    &
          , cca4comb_cld   (qdims%i_start:qdims%i_end,      &
                            qdims%j_start:qdims%j_end,      &
                            n_cca_levels) &
          , ccb4comb_cld   (qdims%i_start:qdims%i_end,      &
                            qdims%j_start:qdims%j_end)      &
          , cct4comb_cld   (qdims%i_start:qdims%i_end,      &
                            qdims%j_start:qdims%j_end) )

        nclds = MIN(cloud_levels, qdims%k_end)

        !
        ! ***** Code adapted from R2_SET_CLOUD_FIELD. *****
        !
        ! Zero cloud amounts in the upper layers (if necessary).
! Nclds_if1:
        IF (qdims%k_end  >   nclds) Then
! Rad_k_do1:
          DO k = 1, qdims%k_end-nclds
            DO j = qdims%j_start,qdims%j_end
              DO i = qdims%i_start,qdims%i_end
                combined_cloud(i, j, k) = 0.0E+00
              END DO
            END DO
          END DO  ! Rad_k_do1
        END IF  ! Nclds_if1

        ! Use Sec 0 convective cloud
        DO k=1, n_cca_levels
          DO j=qdims%j_start,qdims%j_end
            DO i=qdims%i_start,qdims%i_end
              cca4comb_cld(i,j,k) = cca0(i,j,k)
            END DO
          END DO
        END DO
        DO j=qdims%j_start,qdims%j_end
          DO i=qdims%i_start,qdims%i_end
            ccb4comb_cld(i,j) = ccb0(i,j)
            cct4comb_cld(i,j) = cct0(i,j)
          END DO
        END DO

        ! Calculate combined cloud field
        DO k=qdims%k_end+1-nclds, qdims%k_end

          kinvert = qdims%k_end+1 - k

          IF (l_3d_cca) THEN

            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                combined_cloud(i,j,k) = cca4comb_cld(i,j,kinvert)             &
                                      + (1.0 - cca4comb_cld(i,j,kinvert))     &
                                      * area_cloud_fraction(i,j,kinvert)
              END DO
            END DO

          ELSE

            DO j=qdims%j_start,qdims%j_end
              DO i=qdims%i_start,qdims%i_end
                IF ( (cct4comb_cld(i,j) >= kinvert+1) .AND.                   &
                     (ccb4comb_cld(i,j) <= kinvert) ) THEN
                  combined_cloud(i,j,k) = cca4comb_cld(i,j,1)                 &
                                        + (1.0-cca4comb_cld(i,j,1))           &
                                        * area_cloud_fraction(i,j,kinvert)
                ELSE
                  combined_cloud(i,j,k) = area_cloud_fraction(i,j,kinvert)
                END IF
              END DO
            END DO

          END IF  ! l_3d_cca

        END DO  ! k

      ELSE
        ALLOCATE(combined_cloud(1,1,1))

      END IF  ! L_combi_cld_if1


!
! NB: Combined cloud area fractions in each gridbox set up above.
!     Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!
! L_plsp_if1:
      If ( error_code <= 0  .AND.  L_plsp .AND. L_apply_diag ) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        CALL R2_calc_total_cloud_cover(                                 &
               pdims%i_end*pdims%j_end, qdims%k_end, nclds              &
             , IP_CLOUD_MIX_MAX, combined_cloud(1,1,1), work2d_1        &
             , pdims%i_end*pdims%j_end, qdims%k_end                     &
             )


        DO j=qdims%j_start,qdims%j_end
          DO i=qdims%i_start,qdims%i_end
            IF (cca_at_base(i,j) < 1.0) THEN
              plsp(i,j) = MAX( 0.0                                      &
                             , (work2d_1(i,j)-cca_at_base(i,j))         &
                                / (1.0 - cca_at_base(i,j)) )
            ELSE
              plsp(i,j) = 0.0
            END IF
          END DO
        END DO

      End If  ! L_plsp_if1:

      If ( L_apply_diag ) Then

! ----------------------------------------------------------------------
! Section BL.4d Cloud scheme (section 9) diagnostics.
! ----------------------------------------------------------------------
!
! Check that cloud diagnostics are requested this timestep
! DiagSect09_if1:
      If (error_code  ==  0 .AND. sf(0,9)) Then
!
! DEPENDS ON: diagnostics_lscld
        Call Diagnostics_Lscld(                                         &
     &                       pdims%i_end,pdims%j_end, pdims%k_end       &
     &,                      rhc_row_length, rhc_rows                   &
     &,                      qdims%k_end, bl_levels, cloud_levels       &
     &,                      n_rows, global_row_length, global_rows&
     &,                      halo_i, halo_j, offx , offy , me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length                       &
     &,                      at_extremity, p_layer_centres(1, 1, 1)     &
     &,                      p                                          &
     &,                      T_incr_diag_bl, q_incr_diag_bl             &
     &,                      qcl_incr_diag_bl                           &
     &,                      T_latest, q_latest, qcl_latest, qcf_latest &
     &,                      area_cloud_fraction, bulk_cloud_fraction,  &
                             cfl_latest, cff_latest,                    &
                             rho_wet_theta(tdims%i_start:tdims%i_end,   &
                                           tdims%j_start:tdims%j_end,   &
                                           1:tdims%k_end-1)             &
     &,                      p_star, rhcpt                              &
     &,                      combined_cloud, L_murk, Aerosol, RHcrit    &
     &,                      lq_mix_bl,                                 &
                             delta_lambda, delta_phi,                   &
                             FV_cos_theta_latitude,                     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                             STASHwork9,                                &
                             l_combi_cld                                &
                             )
!
      End If  ! DiagSect09_if1
!
! ----------------------------------------------------------------------
! Section BL.5 Energy correction
! ----------------------------------------------------------------------

      If (L_emcorr .AND. Error_code  ==  0) Then

! Add surface sensible heat flux into diabatic heating
! for use in energy correction procedure.

! DEPENDS ON: flux_diag
        Call flux_diag(ftl, FV_cos_theta_latitude,                      &
     &                 pdims%i_end,pdims%j_end,offx ,offy ,             &
     &                 1.0, sum_eng_fluxes,timestep)
! moisture flux level 1 held in fqt
! Should be total moisture flux from surface to layer 1 ie evaporation
! DEPENDS ON: flux_diag
        Call flux_diag(fqt, FV_cos_theta_latitude,                      &
     &                 pdims%i_end,pdims%j_end,offx ,offy ,             &
     &                 1.0, sum_moist_flux,timestep)

      End If   ! L_emcorr

! ----------------------------------------------------------------------
! Section BL.6 Output Diagnostics requested.
! ----------------------------------------------------------------------
      ! Take rho values on the lowest level for visibility diagnostics
      Do j = pdims%j_start,pdims%j_end
        Do i = pdims%i_start,pdims%i_end
          rho1(i,j)=rho(i,j,1)/(r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
        End Do
      End Do

! diagnostics requested this timestep
      If ( sf(0,3) ) Then
        If ( error_code  ==  0) Then

! DEPENDS ON: diagnostics_bl
          Call diagnostics_bl(                                          &
! IN levels / grids / switches
     &                       bl_levels, land_points, dsm_levels         &
     &,                      DIM_CS1, DIM_CS2, l_murk                   &
     &,                      rhc_row_length, rhc_rows                   &
     &,                      sin_theta_longitude, cos_theta_longitude   &
     &,                      land_index,ntiles,npft,nice,nice_use       &
     &,                      l_dust,l_dust_diag, sq_T1p5                &
! IN fields for diagnostics
     &,                      Aerosol(1:tdims%i_end,1:tdims%j_end,1)     &
     &,                      RHcrit, plsp, cca_at_base                  &
     &,                      ls_rain, ls_snow, conv_rain, conv_snow     &
     &,                      p_star, rhcpt, ntml, cumulus, rho1         &
     &,                      qcf_latest                                 &
     &,                      T_incr_diag_bl, q_incr_diag_bl             &
     &,                      qcl_incr_diag_bl, qcf_incr_diag_bl         &
     &,                      cf_incr_diag_bl                            &
     &,                      cfl_incr_diag_bl, cff_incr_diag_bl         &
     &,                      u_incr_diag_bl,v_incr_diag_bl              &
     &,                      t1p5m, zh, u10m, v10m, q1p5m               &
     &,                      e_sea, h_sea,ei                            &
     &,                      sea_ice_htf, sice_mlt_htf                  &
     &,                      snomlt_surf_htf, zht                       &
     &,                      bl_type_1,bl_type_2,bl_type_3,bl_type_4    &
     &,                      bl_type_5,bl_type_6,bl_type_7              &
     &,                      fqt, ftl, z0m_gb, z0m_eff_gb, z0h_eff_gb   &
     &,                      rib_gb, latent_heat, taux, tauy, fme       &
     &,                      t_soil, surf_ht_flux_gb                    &
     &,                      surf_ht_flux_land,surf_ht_flux_sice        &
     &,                      rib_ssi,ftl_ssi,e_ssi,ei_sice              &
     &,                      vshr_land,vshr_ssi                         &
     &,                      taux_land,taux_ssi,tauy_land,tauy_ssi      &
     &,                      radnet_sice                                &
     &,    flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)  &
     &,                      land_sea_mask                              &
     &,                      sil_orog_land,ho2r2_orog,gs,gpp,npp,resp_p &
     &,                      ecan_tile,esoil_tile,gpp_ft,ftl_tile       &
     &,                      npp_ft,resp_p_ft,resp_s,resp_s_tot,cs      &
     &,                      rib_tile,es,ecan,fsmc,radnet_tile          &
     &,                      tstar_tile,canopy,catch,z0m_tile,g_leaf    &
     &,                      t1p5m_tile,q1p5m_tile,le_tile,ei_tile,olr  &
     &,                      epot_tile,tile_frac                        &
     &,                      co2_flux_tot, land_co2, dust_flux          &     
     &,                      dust_emiss_frac,u_s_t_tile,u_s_t_dry_tile  &
     &,                      u_s_std_tile,drydep2, bl_diag              &
! variables required for soil moisture nudging scheme macro
     &,                      rhokh,resfs,chr1p5m,alpha1,ra,wt_ext       &
     &,                      lai_ft,canht_ft,gc                         &
! MGS extra bl vars for UKCA
     &, rhokh_mix, rho_aresist, aresist, resist_b, r_b_dust             &
     &, dtrdz_charney_grid, kent, we_lim, t_frac, zrzi, kent_dsc        &
     &, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! INOUT stash workspace
     & STASHwork3)

        End If   ! error code zero
      Endif      ! sf(0,3)

!
! Clear up allocatable arrays
      IF (error_code  <=  0  .AND.  L_combi_cld) THEN
        DEALLOCATE ( cct4comb_cld   )
        DEALLOCATE ( ccb4comb_cld   )
        DEALLOCATE ( cca4comb_cld   )
      END IF
      
      End If ! L_apply_diag 
      ! Deallocate array needed for a number of places.
      IF(ALLOCATED(combined_cloud)) THEN
        DEALLOCATE ( combined_cloud )
      END IF
! 
! Deallocate BL_diags on last substep 
      IF (L_apply_diag) THEN
         DEALLOCATE(BL_diag%oblen) 
         DEALLOCATE(BL_diag%ustar) 
         DEALLOCATE(BL_diag%wbsurf) 
         DEALLOCATE(BL_diag%gradrich) 
         DEALLOCATE(BL_diag%wstar) 
         DEALLOCATE(BL_diag%dbdz) 
         DEALLOCATE(BL_diag%dvdzm) 
         DEALLOCATE(BL_diag%rhokm) 
         DEALLOCATE(BL_diag%rhokh) 
         DEALLOCATE(BL_diag%tke) 
         DEALLOCATE(BL_diag%ostressx) 
         DEALLOCATE(BL_diag%ostressy) 
         DEALLOCATE(BL_diag%smltop) 
         DEALLOCATE(BL_diag%dsctop) 
         DEALLOCATE(BL_diag%zhlocal) 
         DEALLOCATE(BL_diag%zhpar) 
         DEALLOCATE(BL_diag%dscbase) 
         DEALLOCATE(BL_diag%cldbase) 
         DEALLOCATE(BL_diag%weparm) 
         DEALLOCATE(BL_diag%weparm_dsc) 
         DEALLOCATE(BL_diag%dTfric) 
         DEALLOCATE(BL_diag%elm3D) 
         DEALLOCATE(BL_diag%elh3D) 
         DEALLOCATE(BL_diag%rhoKmloc) 
         DEALLOCATE(BL_diag%rhoKhloc) 
         DEALLOCATE(BL_diag%rhoKmsurf) 
         DEALLOCATE(BL_diag%rhoKhsurf) 
         DEALLOCATE(BL_diag%rhoKmSc) 
         DEALLOCATE(BL_diag%rhoKhSc) 
         DEALLOCATE(BL_diag%weight1d) 
         DEALLOCATE(BL_diag%sgm_trb) 
         DEALLOCATE(BL_diag%ql_trb) 
         DEALLOCATE(BL_diag%cf_trb) 
         DEALLOCATE(BL_diag%wb_ng) 
         DEALLOCATE(BL_diag%sh) 
         DEALLOCATE(BL_diag%sm) 
         DEALLOCATE(BL_diag%tke_dissp) 
         DEALLOCATE(BL_diag%tke_boy_prod) 
         DEALLOCATE(BL_diag%tke_shr_prod) 
         DEALLOCATE(BL_diag%elm) 
         DEALLOCATE(BL_diag%rhogamq) 
         DEALLOCATE(BL_diag%rhogamt) 
         DEALLOCATE(BL_diag%rhogamv) 
         DEALLOCATE(BL_diag%rhogamu)
      END IF !l_apply_diag
!       
!
      If (L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR. L_PC2) Then
        DEALLOCATE (T_earliest)
      End If
!
      DEALLOCATE (q_earliest)
      DEALLOCATE (qcl_earliest)
      DEALLOCATE (qcf_earliest)
      DEALLOCATE (cf_earliest)
      DEALLOCATE (cfl_earliest)
      DEALLOCATE (cff_earliest)
!
      If (L_PC2) Then
        DEALLOCATE ( T_inc_PC2 )
        DEALLOCATE ( q_inc_PC2 )
        DEALLOCATE ( qcl_inc_PC2 )
        DEALLOCATE ( cfl_inc_PC2 )
        DEALLOCATE ( bcf_inc_PC2 )
      End If

      IF (lhook) CALL dr_hook('NI_IMP_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_imp_ctl
