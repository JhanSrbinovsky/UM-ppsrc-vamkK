! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to hold SCM input variables

MODULE s_main_force

  USE ancil_info, ONLY: nsmax
  USE snow_param, ONLY: rho_snow_const
  USE rad_param,  ONLY: r0


  USE s_maxdim,  ONLY:                                                        &
    mx_mod_lv, mx_wet_lv, mx_st_lv, mx_sm_lv, mx_tr_lv, mx_tr_vars            &
  , dim_cs1, dim_cs2, mx_nsprog, mx_nobs

  USE scm_utils, ONLY:                                                        &
    rmdi, imdi, rw_lng, rw, nmod_lv, nwet_lv, nobs, nbl_lv, nlnd, ntile       &
  , o3_lv, st_lv, sm_lv, npft, ntr_lv, ntr_var, ntype, zhook_in, zhook_out    &
  , jprb, lhook, dr_hook

  USE lake_mod, ONLY:                                                         &
    lake_depth_0, lake_fetch_0, lake_h_mxl_0, lake_shape_0, g_dt_0


  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Main module to holds input variables read in from SCM forcing file
!
! Method:
!   Input from the SCM forcing file sub-namelists are transferred to this
!   module after being read from their separate modules. This module
!   holds variable from the following sub-namelists:
!
!   INDATA; RUNDATA; LOGIC; INJULES; INGWD; INGEOFOR; INOBSFOR; INPROF;
!   RADCLOUD; PHYSWITCH; NC_OBS
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  INTEGER, PARAMETER :: &
    row_length = 1      &
  , rows       = 1


  !---------------------------------------------------------------------------
  ! &INDATA
  !---------------------------------------------------------------------------

  INTEGER ::        &
    year_init       &! Initial year
  , month_init      &! Initial month
  , day_init        &! Initial day
  , hour_init       &! Initial hour
  , min_init        &! Initial minute
  , sec_init         ! Initial second

  ! Variables for tape input
  INTEGER ::        &
    tapeyear_init   &! Initial year
  , tapemonth_init  &! Initial month
  , tapeday_init    &! Initial day
  , tapehour_init   &! Initial hour
  , tapemin_init    &! Initial minute
  , tapesec_init     ! Initial second

  INTEGER ::        &
    salt_dim1       &! Dimensions of sea-salt aerosol arrays
  , salt_dim2       &
  , salt_dim3

  INTEGER, ALLOCATABLE :: &
    soil_type(:)     ! Soil type code: 1) Ice,    2) Fine,
                     !                 3) Medium, 4) Coarse

  REAL, ALLOCATABLE :: &
    gridbox_area (:,:) &! Gridbox area (km2)
  , lat          (:,:) &! Latitude,  read from climate dataset
                        ! if STATS forcing chosen
  , long         (:,:)  ! Longitude, read from climate dataset
                        ! if STATS forcing chosen

  LOGICAL :: &
    gather    ! TRUE if running on the Cray


  !---------------------------------------------------------------------------
  ! &RUNDATA
  !---------------------------------------------------------------------------

  INTEGER ::       &
    change_clim    &! No. of days between changes of
                    ! climatological data (default=10)
  , dump_step      &! No. of timesteps for which mean dump is required
                    ! (default of 0 will cause dump_step to be set such that
                    ! time between dumps is close to 1200s)
  , dump_days(4)   &! No. of days for which mean dump is required (default=1)
  , ndayin         &! No. of whole days in run
  , nminin         &! No. of whole minutes in run
  , nsecin         &! No. of seconds in run
                    ! Note: Total run time = ndayin + nminin + nsecin
  , resdump_days   &! frequency of dumps for restart
  , runno_in       &! Run no. to be read from previous run stored on tape
  , runno_out       ! Run no. to be written to tape

  INTEGER ::       &
    min_trop_level &! First reference model level about 700hpa
  , max_trop_level  ! First reference model level about 50hpa

  INTEGER ::       &
    ntrad1          ! Timestep on which to first call Radiation


  REAL :: timestep  ! Model timestep (s)

  REAL ::          &
    co2start       &! co2 mmr at start of year (kg/kg)
  , co2end         &! co2 mmr at end of year (kg/kg)
  , co2rate         ! co2 mmr rate of change (kg/kg)/year

  CHARACTER(LEN=2) ::  &
    modelID         ! Forcing model identifier

  CHARACTER(LEN=8) ::  &
    exname_in      &! Name of expt. to be read from previous
                    ! run stored on tape up to 6 Characters
  , exname_out      ! Name of expt. to be written to tape up to 6 Characters

  INTEGER, ALLOCATABLE :: &
    nbdsc (:,:)           &! Bottom level of Decoupled StratoCumulus layer
  , ntdsc (:,:)           &! Top level of Decoupled StratoCumulus layer
  , ntml  (:,:)            ! Top level of surface mixed layer

  REAL, ALLOCATABLE :: &
    albsoil(:)         &! Bare soil albedo
  , albobs_sw(:)       &! SW obs/clim SW albedo
  , albobs_vis(:)      &! SW obs/clim VIS albedo
  , albobs_nir(:)       ! SW obs/clim NIR albedo

  REAL, ALLOCATABLE :: &
    so2_em      (:,:)  &! Srf emiss. - Sulphur Dioxide
  , nh3_em      (:,:)  &! Srf emiss. - Ammonia
  , dms_em      (:,:)  &! Srf emiss. - DiMethyl Sulphide
  , soot_em     (:,:)  &! Srf emiss. - Soot (fresh)
  , bmass_em    (:,:)  &! Srf emiss. - Biomass (fresh)
  , ocff_em     (:,:)  &! Srf emiss. - Organic Carbon Fossil Fuel (fresh)
  , soot_hilem  (:,:)  &! High lvl emiss. - Soot (fresh)
  , bmass_hilem (:,:)  &! High lvl emiss. - Biomass (fresh)
  , ocff_hilem  (:,:)  &! High lvl emiss. - Organic Carbon Fossil Fuel (fresh)
  , soot        (:,:)   ! Snow soot content (kg/kg)

  REAL, ALLOCATABLE ::   &
    sum_eng_fluxes (:,:) &
  , sum_moist_flux (:,:) &
  , cclwp          (:,:) &! Condensed water path kg/m2
  , orog           (:,:)  ! Orographic height, (m)

  ! Vertical correlation coeff.for ...
  REAL, ALLOCATABLE ::   &
    cort  (:,:)          &! ... Temperature
  , cord  (:,:)          &! ... Dew pt. depression
  , corvn (:,:)          &! ... Velocity VN
  , corw  (:,:)           ! ... Vertical velocity

  REAL, ALLOCATABLE ::   &
    ozone (:,:,:)         ! Ozone profile (kg/kg)

  ! Sulphur Cycle tracers for wet scavenging.
  REAL, ALLOCATABLE ::   &
    so2        (:,:,:)   &! Sulphur dioxide
  , so4_aitken (:,:,:)   &! Sulphate aerosol - Aitken mode
  , so4_accu   (:,:,:)   &! Sulphate aerosol - Accumulation mode
  , so4_diss   (:,:,:)   &! Sulphate aerosol - Dissolved
  , nh3        (:,:,:)   &! Ammonia
  , dms        (:,:,:)    ! DiMethyl Suphide

  ! Aerosols (Used if l_murk=.true., default .false.)
  REAL, ALLOCATABLE ::   &
    aerosol    (:,:,:)   &! MURK aerosol field (in micro-g/kg)
  , aerosol_em (:,:,:)    ! MURK aerosol emissions (in micro-g/kg/s)

  ! Biomass (MMR)
  REAL, ALLOCATABLE ::   &
    bmass_new  (:,:,:)   &! Biomass - Fresh
  , bmass_aged (:,:,:)   &! Biomass - Aged
  , bmass_cld  (:,:,:)    ! Biomass - In-cloud

  ! Soot (MMR)
  REAL, ALLOCATABLE ::   &
    soot_new   (:,:,:)   &! Soot - Fresh
  , soot_aged  (:,:,:)   &! Soot - Aged
  , soot_cld   (:,:,:)    ! Soot - In-cloud

  ! Organic Carbon Fossil-Fuel
  REAL, ALLOCATABLE ::   &
    ocff_new   (:,:,:)   &! OCFF - Fresh
  , ocff_aged  (:,:,:)   &! OCFF - Aged
  , ocff_cld   (:,:,:)    ! OCFF - In-cloud

  ! Nitrate aerosol
  REAL, ALLOCATABLE::    &
    nitr_acc   (:,:,:)   &! Nitrate - Accumulated
  , nitr_diss  (:,:,:)    ! Nitrate - Dissolved

  ! Carbon cycle
  REAL, ALLOCATABLE ::   &
    co2_emits  (:,:)     &! Srf emiss. - Carbon dioxide (kg/m2/s)
  , co2flux    (:,:)     &
  , co2        (:,:,:)    ! Carbon dioxide

  REAL, ALLOCATABLE ::   &
    dOLR_rts    (:,:)    &! TOA - surface upward LW
  , fland_ctile (:,:)    &! Land fraction on land points. (Coastal tiling)
  , tstar_land  (:,:)    &! Land mean sfc temperature (K)
  , tstar_sea   (:,:)    &! Open sea sfc temperature (K).
  , tstar_sice  (:,:)    &! Sea-ice sfc temperature (K).
  , land_alb    (:,:)    &! Mean land albedo.
  , sice_alb    (:,:)    &! Sea-ice albedo.
  , ddmfx       (:,:)    &! Convective downdraught mass flux at cloud-base
  , zh          (:,:)     ! Height above surface of top of boundary layer (m)

  ! Prognostics for the screen-level temperature
  REAL, ALLOCATABLE ::   &
    tstbtrans(:,:)       &! Time since the last transition to
                          ! stability at the surface
  , tscrndcl_ssi(:,:)    &! Decoupled screen-level temperature
                          ! over sea and sea-ice
  , tscrndcl_tile(:,:)    ! Decoupled screen-level temperature land points



  !---------------------------------------------------------------------------
  ! &LOGIC
  !---------------------------------------------------------------------------
  ! All logical comments refer to a logical state of .TRUE.

  LOGICAL ::        &
    dev_test        &! General logical for development testing
  , test            &! Use detailed sub-timestep diagnostics
  , altdat          &! Use namelist initial profiles of T, q, u and v
  , geoforce        &! Use geostrophic wind forcing
  , geoinit         &! Initialise dump to geostrophic
  , l_geo_centred   &! Use centred discretization to calculate
                     ! increment from geostrophic forcing
  , grafdump_day    &! Output graphical (mean daily values)
  , grafdump_days   &! Output graphical (mean values over DUMP_DAYS)
  , grafdump_step   &! Output graphical (mean values over DUMP_STEP)
  , noforce         &! No large-scale forcing required
  , obs             &! Use observational large-scale forcing
  , obs_surf        &! Use observational surface forcing
  , prindump_day    &! Printout mean daily dump
  , prindump_days   &! Printout mean dump over DUMP_DAYS
  , prindump_obs    &! Printout observational diagnostics every OBS_ timesteps
  , prindump_step   &! Printout mean dump each DUMP_STEP
  , prinstat        &! Printout stats forcing every timestep
  , radcloud_fixed  &! Cloud fixed for radiation
  , stats           &! Use statistical large-scale forcing
  , local_time      &! Use local time for diagnostics instead of GMT
  , l_flux_bc       &! Use prescribed surface fluxes, otherwise T-star forcing
                     ! used. (Sea-point only)
  , ancyc           &! Use annual cycle
                     ! (ie. radiation input then varies throughout year)
  , tapein          &! Read initial data from previous run stored on tape
  , tapeout         &! Store restart information and diagnostic
                     ! output to be stored on tape routines
  , l_spec_z0       &! Use prescibed surface roughness length
  , l_qpos_for       ! Ensure q >= qlimit if large-scale forcing
                     ! causes q < qlimit.

  LOGICAL, ALLOCATABLE :: &
    land_sea_mask (:,:)   &! Land point flag (any type)
  , land_ice_mask (:,:)   &! Land point flag (ice)
  , soil_mask     (:,:)   &! Land point flag (soil)
  , cumulus       (:,:)    ! BL convection flag


  !---------------------------------------------------------------------------
  ! &INJULES
  !---------------------------------------------------------------------------

  INTEGER :: &
    smi_opt  ! Option to define method of initialisation of soil
             ! moisture content: 0: Use smcli, 1: Use fsmc, 2: Use sth

  REAL, ALLOCATABLE :: &
    gs       (:)       &! Stomatal conductance
  , fsmc     (:)       &! Soil Moisture Stress Factor
  , frac_typ (:,:)     &! Fractions of surface types
  , smcli    (:,:)     &! Initial SMC profile in layers (kg/m2)
  , sth      (:,:)      ! Total soil moisture in layers
                        ! as fraction of saturation

  REAL, ALLOCATABLE :: &
    canht    (:,:)     &! Canopy height (m)
  , lai      (:,:)      ! Leaf area index


  REAL, ALLOCATABLE :: &
    z0_tile    (:,:)   &! Tile roughness lengths (m).
  , z0h_tile   (:,:)   &! Tile thermal roughness lengths (m).
  , rgrain     (:,:)   &! Snow grain size (microns)
  , infil_tile (:,:)   &! Max. surface infiltration
  , snow_tile  (:,:)   &! Lying snow on tiles not sure if JULES uses this
  , tstar_tile (:,:)   &! Tile surface temperatures, (K)
  , catch      (:,:)   &! Surf/canopy water capacity
                        ! (snow-free land tiles) (kg/m2)
  , canopy     (:,:)    ! Surf/canopy water
                        ! (snow-free land tiles) (kg/m2)


  ! 2A Triffid variables
  !=====================
  REAL, ALLOCATABLE ::    &
    frac_disturb    (:)   &! Fraction of gridbox in which!
                           ! vegetation is disturbed.
  , npp_ft_acc      (:)   &! Acc. npp_ft
  , resp_w_ft_acc   (:,:) &! Acc. resp_w_ft
  , resp_s_acc      (:,:) &! Acc. resp_s
  , cs              (:,:) &! Soil carbon (kg C/m2)
  , g_leaf_phen_acc (:,:) &! Acc. leaf turnover rate with phenology
  , g_leaf_acc      (:,:)  ! Acc. g_leaf


  REAL, ALLOCATABLE ::    &
    lw_down (:,:)          ! Surface downward LW

  ! FLake lake scheme prognostics
  !===============================
  REAL, ALLOCATABLE ::  &
    lake_depth      (:) &! lake depth (m)
  , lake_fetch      (:) &! typical wind fetch (m)
  , lake_t_mean     (:) &! lake mean Temp (K)
  , lake_t_mxl      (:) &! lake mixed-layer Temp (K)
  , lake_t_ice      (:) &! Temp at ice upper surface (K)
  , lake_h_mxl      (:) &! lake mixed-layer thickness (m)
  , lake_h_ice      (:) &! lake ice thickness (m)
  , lake_shape      (:) &! thermocline shape factor
  , lake_g_dt       (:)  ! lake ht.flx / dT (W m-2 K-1)

  ! Jules scheme soil parameters
  !===============================
  REAL, ALLOCATABLE ::    &
    clapp_levs      (:,:) &
  , sathh_levs      (:,:) &
  , hcap_levs       (:,:) &
  , hcon_levs       (:,:) &
  , satcon_levs     (:,:) &
  , smvccl_levs     (:,:) &
  , smvcwt_levs     (:,:) &
  , smvcst_levs     (:,:)


  !---------------------------------------------------------------------------
  ! &INGWD
  !---------------------------------------------------------------------------
  REAL ::                 &
    sd_orog_land          &! Ancillary fields and fields needed
  , orog_grad_xx_land     &! to be kept from timestep to
  , orog_grad_xy_land     &! timestep- in atmos_physics2 call
  , orog_grad_yy_land


  !---------------------------------------------------------------------------
  ! &INGEOFOR
  !---------------------------------------------------------------------------

  REAL, ALLOCATABLE :: &
    ug(:,:,:,:)        &! Geostrophic U velocity (m/s)
  , vg(:,:,:,:)         ! Geostrophic V velocity (m/s)

  INTEGER ::           &
    ug_opt             &! Specifies u-wind geostrophic forcing format
  , vg_opt              ! Specifies v-wind geostrophic forcing format


  !---------------------------------------------------------------------------
  ! &INOBSFOR
  !---------------------------------------------------------------------------

  ! Relaxation options
  !   0 = No relaxation
  !   1 = Relax back to initial conditions across tau timescale
  !   2 = Relax back to background profile across tau timescale
  !   3 = Relax back to initial across 1 timstep
  !   4 = Relax back to background across 1 timstep

  ! Relaxation option for ...
  INTEGER ::   &
    rlx_t      &! ... Temperature
  , rlx_q      &! ... Spec. humid.
  , rlx_u      &! ... Zonal wind
  , rlx_v      &! ... Merid wind
  , rlx_w       ! ... vert. wind

  REAL ::      &
    obs_pd     &! Period between observational profiles (s)
  , obs_top    &! Height above which observations are replaced
                ! with constants/standard profiles. (m)
  , obs_bot     ! Height below which observations are replaced
                ! with constants/standard profiles. (m)

  ! Pressure threshold (Pa) below which relax ...
  REAL ::      &
    plev_t     &! ... Temperature
  , plev_q     &! ... Spec. humid.
  , plev_u     &! ... Zonal wind
  , plev_v     &! ... Merid wind
  , plev_w      ! ... Vert. wind

  ! Relaxation timescale (s) for ...
  REAL ::      &
    tau_t      &! ... Temperature
  , tau_q      &! ... Spec. humid.
  , tau_u      &! ... Zonal wind
  , tau_v      &! ... Merid wind
  , tau_w       ! ... Vert. wind

  LOGICAL ::   &
    l_vertadv   ! Use interactive vertical advection
                ! NOTE: namelists forcings must be from
                !       HORIZONTAL ADVECTION ONLY IF THIS OPTION IS TRUE

  ! Forcing tendency due to large-scale horizontal/vertical advection of ...
  REAL, ALLOCATABLE ::    &
    q_star(:,:,:,:)       &! ... Spec. Humid, (kg/kg)/day
  , t_inc (:,:,:,:)       &! ... Temperature, (K)/day
  , u_inc (:,:,:,:)       &! ... Zonal wind,  (m/s)/day
  , v_inc (:,:,:,:)       &! ... Merid wind,  (m/s)/day
  , w_inc (:,:,:,:)        ! ... Vert. wind,  (m/s)/day

  ! Observed background fields of ...
  REAL, ALLOCATABLE ::    &
    q_bg  (:,:,:,:)       &! ... Temperature  (K)
  , t_bg  (:,:,:,:)       &! ... Spec. humid. (kg/kg)
  , u_bg  (:,:,:,:)       &! ... Zonal wind   (m/s)
  , v_bg  (:,:,:,:)       &! ... Merid wind   (m/s)
  , w_bg  (:,:,:,:)        ! ... Vert. wind   (m/s)

  ! Surface forcing of ...
  REAL, ALLOCATABLE ::    &
    flux_h       (:,:,:)  &! ... Sensible heat flux  (W/m2)
  , flux_e       (:,:,:)  &! ... Latent heat flux    (W/m2)
  , tstar_forcing(:,:,:)   ! ... Temperature forcing (K)


  !---------------------------------------------------------------------------
  ! &INPROF
  !---------------------------------------------------------------------------

  INTEGER ::              &
    nml_inprof_thetal      ! Switch for interpretation of
                           ! inital theta, qi profiles

  INTEGER, ALLOCATABLE :: &
    iccbi (:,:)           &! Initial Convective cloud base level
  , iccti (:,:)            ! Initial Convective cloud top level


  REAL, ALLOCATABLE ::    &
    canopy_gbi (:)        &! Initial gridbox mean canopy water content (kg/m2)
  , smci       (:)         ! Initial soil moisture content (kg/m2)


  REAL, ALLOCATABLE ::    &
    t_deep_soili (:,:)     ! Initial deep soil temperature profile, (K)

  REAL, ALLOCATABLE ::    &
    ccai(:,:)             &! Initial convective cloud amount (decimal frac.)
  , snodepi(:,:)          &! Initial snow depth (kg/m2)
  , tstari(:,:)           &! Initial surface temperature, (K)
  , z0mseai(:,:)          &! Init. sea srf roughness length (momentum), (m)
  , z0m_scm(:,:)          &! Fixed sea srf roughness length (momentum), (m)
  , z0h_scm(:,:)          &! Fixed sea srf roughness length (heat), (m)
  , sil_orog_land(:,:)    &! Silhouette area of unresolved orography
                           ! per unit horizontal area on land points only.
  , ho2r2_orog(:,:)       &! Standard Deviation of orography equivalent to
                           ! peak to trough height of unresolved orography
                           ! divided by 2SQRT(2) on land points only (m)
  , ice_fract(:,:)        &! Fraction of grid box covered by sea ice
                           ! (decimal fraction)
  , di(:,:)               &! Equivalent thickness of sea-ice (m)
  , u_0(:,:)              &! Comp. of srf current, westerly (m/s)
  , v_0(:,:)               ! Comp. of srf current, southernly (m/s)

  ! JULES
  REAL, ALLOCATABLE ::      &
    i_snowdepth     (:,:)   &! Init. snow depth (m)
  , i_snow_grnd     (:,:)   &! Under canopy snow store (kg/m2)
  , i_rho_snow_grnd (:,:)   &! Init. snow density (kg/m3)
  , nsnow           (:,:)   &! Number of snow levels (Real)
  , i_tsnowlayer    (:,:,:) &! Init. snow layer temperature (K)
  , i_sice          (:,:,:) &! Init. ice in snow pack (kg/m2)
  , i_sliq          (:,:,:) &! Init. liquid in snow pack (kg/m2)
  , i_ds            (:,:,:) &! Init. snow layer thickness (m)
  , i_rgrainl       (:,:,:)  ! Snow grain size in each layer (microns)

  REAL, ALLOCATABLE ::    &
    free_tracers(:,:,:,:)  ! Model tracer fields (kg/kg)

  ! Initial profiles of ...
  REAL, ALLOCATABLE ::    &
    qi     (:,:,:)        &! ... Spec. humid. (kg/kg)
  , theta  (:,:,:)        &! ... Pot. temp.   (K)
  , ui     (:,:,:)        &! ... Zonal  wind  (m/s)
  , vi     (:,:,:)        &! ... Merid. wind  (m/s)
  , wi     (:,:,:)        &! ... Vert   wind  (m/s)
  , p_in   (:,:,:)        &! ... Pressure     (Pa) (rh-levs)
  , w_advi (:,:,:)

  LOGICAL :: &
    kill_interp


  !---------------------------------------------------------------------------
  ! &RADCLOUD
  !---------------------------------------------------------------------------

  INTEGER, ALLOCATABLE ::  &
    iccb_rad(:,:)          &! Fixed convective cloud base level
  , icct_rad(:,:)           ! Fixed convective cloud top level

  REAL, ALLOCATABLE ::     &
    cca_rad    (:,:)       &! Fixed convective cloud amount (fraction)
  , ccwpin_rad (:,:)       &! Fixed convective cloud water path (kg/m2).
  , qcl_rad    (:,:,:)     &! Total cloud water and ice
                            ! content over cloud (kg/kg)
  , qcf_rad    (:,:,:)     &! Set to zero as user will usually
                            ! input combined cloud water and
                            ! ice content over cloud (kg/kg) in QCL_RAD.
  , layer_cloud_rad (:,:,:) ! Layer cloud amount (fraction)


  !---------------------------------------------------------------------------
  ! &PHYSWITCH
  !---------------------------------------------------------------------------

  ! conv_mode determines actions for convection scheme
  !   0 = Run normally
  !   1 = Run for diagnostics every radiation timestep but save dump state
  !       (except CCA,ICCB and ICCT)
  !   2 = Don't run

  INTEGER :: conv_mode

  !---------------------------------------------------------------------------
  ! &NC_obs
  !---------------------------------------------------------------------------

  CHARACTER(LEN=200) :: netCDF_file
  INTEGER        :: obs_t0_prf
  INTEGER        :: source


!=============================================================================

CONTAINS
!-----------------------------------------------------------------------------
  SUBROUTINE set_forcing_defaults

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('SET_FORCING_DEFAULTS',zhook_in,zhook_handle)

    ! &INDATA
    !=========
    gridbox_area (:,:) = 100000.0
    lat          (:,:) = 0.0
    long         (:,:) = 0.0
    soil_type    (:)   = imdi
    year_init          = 1998
    month_init         = 1
    day_init           = 1
    hour_init          = 0
    min_init           = 0
    sec_init           = 0
    tapeyear_init      = 1998
    tapemonth_init     = 1
    tapeday_init       = 1
    tapehour_init      = 0
    tapemin_init       = 0
    tapesec_init       = 0
    salt_dim1          = 1
    salt_dim2          = 1
    salt_dim3          = nmod_lv
    gather             = .FALSE.


    ! &INOBSFOR
    !===========
    l_vertadv = .FALSE.
    obs_pd    = rmdi
    obs_bot   = 20.0
    obs_top   = 40000.0
    rlx_t     = 0
    rlx_q     = 0
    rlx_u     = 2
    rlx_v     = 2
    rlx_w     = 4
    plev_t    = 1.5e5
    plev_q    = 1.5e5
    plev_u    = 1.5e5
    plev_v    = 1.5e5
    plev_w    = 1.5e5
    tau_t     = 3600.0
    tau_q     = 3600.0
    tau_u     = 3600.0
    tau_v     = 3600.0
    tau_w     = 3600.0

    flux_h        (:,:,:) = rmdi
    flux_e        (:,:,:) = rmdi
    tstar_forcing (:,:,:) = rmdi

    t_bg   (:,:,:,:) = rmdi
    q_bg   (:,:,:,:) = rmdi
    u_bg   (:,:,:,:) = rmdi
    v_bg   (:,:,:,:) = rmdi
    w_bg   (:,:,:,:) = rmdi
    u_inc  (:,:,:,:) = 0.0
    v_inc  (:,:,:,:) = 0.0
    w_inc  (:,:,:,:) = 0.0
    t_inc  (:,:,:,:) = 0.0
    q_star (:,:,:,:) = 0.0


    ! &LOGIC
    !===========
    ancyc          = .TRUE.
    altdat         = .TRUE.
    local_time     = .TRUE.
    obs            = .FALSE.
    obs_surf       = .FALSE.
    prindump_obs   = .FALSE.
    stats          = .FALSE.
    noforce        = .FALSE.
    geoforce       = .FALSE.
    geoinit        = .FALSE.
    test           = .FALSE.
    prinstat       = .FALSE.
    radcloud_fixed = .FALSE.
    prindump_day   = .FALSE.
    prindump_days  = .FALSE.
    prindump_step  = .FALSE.
    grafdump_day   = .FALSE.
    grafdump_days  = .FALSE.
    grafdump_step  = .FALSE.
    tapein         = .FALSE.
    tapeout        = .FALSE.
    l_geo_centred  = .FALSE.
    l_flux_bc      = .FALSE.
    l_spec_z0      = .FALSE.
    l_qpos_for     = .TRUE.

    land_ice_mask (:,:) = .FALSE.
    land_sea_mask (:,:) = .FALSE.
    soil_mask     (:,:) = .FALSE.
    cumulus       (:,:) = .FALSE.


    ! &INGEOFOR
    !===========
    ug (:,:,:,:) = 0.0
    vg (:,:,:,:) = 0.0
    ug_opt   = 0
    vg_opt   = 0


    ! &INJULES
    !===========
    smi_opt       = imdi
    lw_down (:,:) = 0.0

    IF (nlnd > 0) THEN
      gs         (:)   = 1.0
      fsmc       (:)   = rmdi
      sth        (:,:) = rmdi
      smcli      (:,:) = rmdi
      canopy     (:,:) = 0.0
      rgrain     (:,:) = r0
      frac_typ   (:,:) = rmdi
      canht      (:,:) = rmdi
      lai        (:,:) = rmdi
      catch      (:,:) = rmdi

      snow_tile  (:,:) = 0.0
      z0_tile    (:,:) = rmdi
      z0h_tile   (:,:) = rmdi
      infil_tile (:,:) = rmdi
      tstar_tile (:,:) = rmdi

      ! 2A Triffid variables
      !=====================
      frac_disturb    (:)   = rmdi
      npp_ft_acc      (:)   = rmdi
      resp_w_ft_acc   (:,:) = rmdi
      resp_s_acc      (:,:) = rmdi
      cs              (:,:) = rmdi
      g_leaf_phen_acc (:,:) = rmdi
      g_leaf_acc      (:,:) = rmdi

      ! JULES
      i_snowdepth       (:,:) = 0.0
      i_snow_grnd       (:,:) = rmdi
      i_rho_snow_grnd   (:,:) = rho_snow_const
      nsnow             (:,:) = 0
      i_tsnowlayer    (:,:,:) = rmdi
      i_sice          (:,:,:) = rmdi
      i_sliq          (:,:,:) = rmdi
      i_ds            (:,:,:) = rmdi
      i_rgrainl       (:,:,:) = rmdi

      ! FLake lake scheme prognostics
      !===============================
      lake_depth  (:) = lake_depth_0
      lake_fetch  (:) = lake_fetch_0
      lake_t_mean (:) = rmdi
      lake_t_mxl  (:) = rmdi
      lake_t_ice  (:) = rmdi
      lake_h_mxl  (:) = lake_h_mxl_0
      lake_h_ice  (:) = rmdi
      lake_shape  (:) = lake_shape_0
      lake_g_dt   (:) = g_dt_0

      clapp_levs  (:,:) = rmdi
      sathh_levs  (:,:) = rmdi
      hcap_levs   (:,:) = rmdi
      hcon_levs   (:,:) = rmdi
      satcon_levs (:,:) = rmdi
      smvccl_levs (:,:) = rmdi
      smvcwt_levs (:,:) = rmdi
      smvcst_levs (:,:) = rmdi

    END IF


    ! &INGWD
    !===========
    sd_orog_land      = 1000.0
    orog_grad_xx_land = 1.0e-3
    orog_grad_xy_land = 1.0e-3
    orog_grad_yy_land = 1.0e-3


    ! &INPROF
    !===========
    nml_inprof_thetal = 0

    theta (:,:,:) = rmdi
    p_in  (:,:,:) = rmdi
    qi    (:,:,:) = rmdi
    ui    (:,:,:) = rmdi
    vi    (:,:,:) = rmdi
    wi    (:,:,:) = 0.0

    IF (nlnd > 0) THEN
      t_deep_soili (:,:) = rmdi
      canopy_gbi   (:)   = rmdi
      smci         (:)   = rmdi
    END IF

    di        (:,:) = 0       ! These model variables are set up as
    ice_fract (:,:) = 0.0     ! constants for a land point.
    snodepi   (:,:) = 0.0
    tstari    (:,:) = 0.0
    z0mseai   (:,:) = 0.0001
    z0m_scm   (:,:) = 0.0
    z0h_scm   (:,:) = 0.0
    ccai      (:,:) = 0.0
    iccbi     (:,:) = 0
    iccti     (:,:) = 0
    u_0       (:,:) = 0.0
    v_0       (:,:) = 0.0

    sil_orog_land (:,:)     = 0.0
    ho2r2_orog    (:,:)     = 0.0
    free_tracers(:,:,:,:)   = 0.0


    ! &RUNDATA
    !===========
    modelID        ='XX'
    dump_step      = 0
    dump_days(:)   = 1
    resdump_days   = 1
    change_clim    = 10
    exname_in      = 'XXXXXXXX'
    exname_out     = 'XXXXXXXX'
    timestep       = 1800.0
    ntrad1         = 1
    ndayin         = imdi
    nminin         = imdi
    nsecin         = imdi
    runno_in       = 0
    runno_out      = 999
    min_trop_level = 0
    max_trop_level = 0
    co2start       = 5.42e-4
    co2end         = 5.42e-4
    co2rate        = 0.0

    cort        (:,:) = 0.9
    cord        (:,:) = 0.9
    corvn       (:,:) = 0.5
    corw        (:,:) = 0.5
    orog        (:,:) = 0.0
    tstar_land  (:,:) = rmdi
    tstar_sea   (:,:) = rmdi
    tstar_sice  (:,:) = rmdi
    land_alb    (:,:) = 0.0
    sice_alb    (:,:) = 0.0
    fland_ctile (:,:) = 0.0
    nbdsc       (:,:) = 0
    ntdsc       (:,:) = 0
    ntml        (:,:) = nbl_lv
    zh          (:,:) = 500.0
    ddmfx       (:,:) = 0.0

    co2_emits   (:,:) = 0.0
    co2flux     (:,:) = 0.0
    cclwp       (:,:) = 0.0
    dolr_rts    (:,:) = 0.0

    so2_em      (:,:) = 0.0
    nh3_em      (:,:) = 0.0
    dms_em      (:,:) = 0.0
    bmass_em    (:,:) = 0.0
    ocff_em     (:,:) = 0.0
    soot_em     (:,:) = 0.0
    soot        (:,:) = 0.0

    aerosol     (:,:,:) = 0.0
    aerosol_em  (:,:,:) = 0.0
    co2         (:,:,:) = 0.0
    nh3         (:,:,:) = 0.0
    dms         (:,:,:) = 0.0
    so2         (:,:,:) = 0.0

    nitr_acc    (:,:,:) = 0.0
    nitr_diss   (:,:,:) = 0.0
    so4_aitken  (:,:,:) = 0.0
    so4_accu    (:,:,:) = 0.0
    so4_diss    (:,:,:) = 0.0

    soot_new    (:,:,:) = 0.0
    bmass_new   (:,:,:) = 0.0
    ocff_new    (:,:,:) = 0.0

    soot_cld    (:,:,:) = 0.0
    bmass_cld   (:,:,:) = 0.0
    ocff_cld    (:,:,:) = 0.0

    soot_aged   (:,:,:) = 0.0
    bmass_aged  (:,:,:) = 0.0
    ocff_aged   (:,:,:) = 0.0

    soot_hilem  (:,:)   = 0.0
    bmass_hilem (:,:)   = 0.0
    ocff_hilem  (:,:)   = 0.0

    ozone       (:,:,:) = rmdi

    sum_eng_fluxes (:,:) = 0.0
    sum_moist_flux (:,:) = 0.0

    IF (nlnd > 0) THEN
      albsoil       (:)   = rmdi
      tscrndcl_tile (:,:) = rmdi
    END IF

    albobs_sw (:)        = rmdi
    albobs_vis(:)        = rmdi
    albobs_nir(:)        = rmdi

    tscrndcl_ssi (:,:) = rmdi
    tstbtrans    (:,:) = rmdi


    ! &RADCLOUD
    !===========
    iccb_rad   (:,:) = 0
    icct_rad   (:,:) = 0
    cca_rad    (:,:) = 0.0
    ccwpin_rad (:,:) = 0.0

    layer_cloud_rad (:,:,:) = 0.0
    qcl_rad         (:,:,:) = 0.0
    qcf_rad         (:,:,:) = 0.0


    ! &PHYSWITCH
    !===========
    conv_mode = 0


    IF (lhook) CALL dr_hook('SET_FORCING_DEFAULTS',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE set_forcing_defaults

!=============================================================================

  SUBROUTINE alloc_forcing

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_FORCING',zhook_in,zhook_handle)

    ! INDATA
    ALLOCATE                        &
      ( soil_type    (nlnd)         &
      , gridbox_area (rw_lng,rw)    &
      , lat          (rw_lng,rw)    &
      , long         (rw_lng,rw) )


    ! RUNDATA
    ALLOCATE                        &
      ( albsoil        (nlnd)       &
      , albobs_sw      (nlnd)       &
      , albobs_vis     (nlnd)       &
      , albobs_nir     (nlnd)       &
      , tscrndcl_tile  (nlnd,ntile) &
      , ntml           (rw_lng,rw)  &
      , nbdsc          (rw_lng,rw)  &
      , ntdsc          (rw_lng,rw)  &
      , cort           (rw_lng,rw)  &
      , cord           (rw_lng,rw)  &
      , corvn          (rw_lng,rw)  &
      , corw           (rw_lng,rw)  &
      , cclwp          (rw_lng,rw)  &
      , orog           (rw_lng,rw)  &
      , co2_emits      (rw_lng,rw)  &
      , co2flux        (rw_lng,rw)  &
      , dolr_rts       (rw_lng,rw)  &
      , fland_ctile    (rw_lng,rw)  &
      , tstar_land     (rw_lng,rw)  &
      , tstar_sea      (rw_lng,rw)  &
      , tstar_sice     (rw_lng,rw)  &
      , land_alb       (rw_lng,rw)  &
      , sice_alb       (rw_lng,rw)  &
      , ddmfx          (rw_lng,rw)  &
      , zh             (rw_lng,rw)  &
      , tstbtrans      (rw_lng,rw)  &
      , tscrndcl_ssi   (rw_lng,rw)  &
      , sum_eng_fluxes (rw_lng,rw)  &
      , sum_moist_flux (rw_lng,rw)  &
      , so2_em         (rw_lng,rw)  &
      , nh3_em         (rw_lng,rw)  &
      , dms_em         (rw_lng,rw)  &
      , soot_em        (rw_lng,rw)  &
      , bmass_em       (rw_lng,rw)  &
      , ocff_em        (rw_lng,rw)  &
      , soot_hilem     (rw_lng,rw)  &
      , bmass_hilem    (rw_lng,rw)  &
      , ocff_hilem     (rw_lng,rw)  &
      , soot           (rw_lng,rw)  &

      , ozone      (rw_lng,rw,o3_lv)     &
      , so2        (rw_lng,rw,nwet_lv)   &
      , so4_aitken (rw_lng,rw,nwet_lv)   &
      , so4_accu   (rw_lng,rw,nwet_lv)   &
      , so4_diss   (rw_lng,rw,nwet_lv)   &
      , nh3        (rw_lng,rw,nmod_lv)   &
      , dms        (rw_lng,rw,nmod_lv)   &
      , aerosol    (rw_lng,rw,nmod_lv)   &
      , aerosol_em (rw_lng,rw,nmod_lv)   &
      , bmass_new  (rw_lng,rw,nmod_lv)   &
      , bmass_aged (rw_lng,rw,nmod_lv)   &
      , bmass_cld  (rw_lng,rw,nmod_lv)   &
      , soot_new   (rw_lng,rw,nmod_lv)   &
      , soot_aged  (rw_lng,rw,nmod_lv)   &
      , soot_cld   (rw_lng,rw,nmod_lv)   &
      , ocff_new   (rw_lng,rw,nmod_lv)   &
      , ocff_aged  (rw_lng,rw,nmod_lv)   &
      , ocff_cld   (rw_lng,rw,nmod_lv)   &
      , nitr_acc   (rw_lng,rw,nmod_lv)   &
      , nitr_diss  (rw_lng,rw,nmod_lv)   &
      , co2        (rw_lng,rw,nmod_lv) )


    ! INPROF
    ALLOCATE                                       &
      ( t_deep_soili    (nlnd,st_lv)               &
      , canopy_gbi      (nlnd)                     &
      , smci            (nlnd)                     &
      , nsnow           (nlnd,ntile)               &
      , i_snowdepth     (nlnd,ntile)               &
      , i_snow_grnd     (nlnd,ntile)               &
      , i_rho_snow_grnd (nlnd,ntile)               &
      , i_tsnowlayer    (nlnd,ntile,nsmax)         &
      , i_sice          (nlnd,ntile,nsmax)         &
      , i_sliq          (nlnd,ntile,nsmax)         &
      , i_ds            (nlnd,ntile,nsmax)         &
      , i_rgrainl       (nlnd,ntile,nsmax)         &
      , iccbi         (rw_lng,rw)                  &
      , iccti         (rw_lng,rw)                  &
      , ccai          (rw_lng,rw)                  &
      , snodepi       (rw_lng,rw)                  &
      , tstari        (rw_lng,rw)                  &
      , z0mseai       (rw_lng,rw)                  &
      , z0m_scm       (rw_lng,rw)                  &
      , z0h_scm       (rw_lng,rw)                  &
      , di            (rw_lng,rw)                  &
      , u_0           (rw_lng,rw)                  &
      , v_0           (rw_lng,rw)                  &
      , sil_orog_land (rw_lng,rw)                  &
      , ho2r2_orog    (rw_lng,rw)                  &
      , ice_fract     (rw_lng,rw)                  &
      , qi            (rw_lng,rw,nwet_lv)          &
      , theta         (rw_lng,rw,nmod_lv)          &
      , ui            (rw_lng,rw,nmod_lv)          &
      , vi            (rw_lng,rw,nmod_lv)          &
      , wi            (rw_lng,rw,0:nmod_lv)        &
      , p_in          (rw_lng,rw,nmod_lv+1)        &
      , w_advi        (rw_lng,rw,0:nmod_lv)        &
      , free_tracers  (rw_lng,rw,ntr_lv,ntr_var) )


    ! INOBSFOR
    ALLOCATE                                       &
      ( q_star        (rw_lng,rw,nobs,nwet_lv)     &
      , q_bg          (rw_lng,rw,nobs,nwet_lv)     &
      , t_inc         (rw_lng,rw,nobs,nmod_lv)     &
      , u_inc         (rw_lng,rw,nobs,nmod_lv)     &
      , v_inc         (rw_lng,rw,nobs,nmod_lv)     &
      , t_bg          (rw_lng,rw,nobs,nmod_lv)     &
      , u_bg          (rw_lng,rw,nobs,nmod_lv)     &
      , v_bg          (rw_lng,rw,nobs,nmod_lv)     &
      , w_inc         (rw_lng,rw,nobs,0:nmod_lv)   &
      , w_bg          (rw_lng,rw,nobs,0:nmod_lv)   &
      , flux_h        (rw_lng,rw,nobs)             &
      , flux_e        (rw_lng,rw,nobs)             &
      , tstar_forcing (rw_lng,rw,nobs) )


    ! LOGIC
    ALLOCATE                                       &
      ( land_sea_mask (rw_lng,rw)                  &
      , land_ice_mask (rw_lng,rw)                  &
      , soil_mask     (rw_lng,rw)                  &
      , cumulus       (rw_lng,rw) )


    ! INJULES
    ALLOCATE                                       &
      ( gs              (nlnd)                     &
      , fsmc            (nlnd)                     &
      , frac_disturb    (nlnd)                     &
      , npp_ft_acc      (nlnd)                     &
      , frac_typ        (nlnd,ntype)               &
      , smcli           (nlnd,sm_lv)               &
      , sth             (nlnd,sm_lv)               &
      , canht           (nlnd,npft)                &
      , lai             (nlnd,npft)                &
      , z0_tile         (nlnd,ntile)               &
      , z0h_tile        (nlnd,ntile)               &
      , rgrain          (nlnd,ntile)               &
      , infil_tile      (nlnd,ntile)               &
      , snow_tile       (nlnd,ntile)               &
      , tstar_tile      (nlnd,ntile)               &
      , catch           (nlnd,ntile)               &
      , canopy          (nlnd,ntile)               &
      , resp_s_acc      (nlnd,dim_cs1)             &
      , cs              (nlnd,dim_cs1)             &
      , resp_w_ft_acc   (nlnd,npft)                &
      , g_leaf_phen_acc (nlnd,npft)                &
      , g_leaf_acc      (nlnd,npft)                &
      , lw_down         (rw_lng,rw)                &
      , lake_depth      (nlnd)                     &
      , lake_fetch      (nlnd)                     &
      , lake_t_mean     (nlnd)                     &
      , lake_t_mxl      (nlnd)                     &
      , lake_t_ice      (nlnd)                     &
      , lake_h_mxl      (nlnd)                     &
      , lake_h_ice      (nlnd)                     &
      , lake_shape      (nlnd)                     &
      , lake_g_dt       (nlnd)                     &
      , clapp_levs      (nlnd,sm_lv)               &
      , sathh_levs      (nlnd,sm_lv)               &
      , hcap_levs       (nlnd,sm_lv)               &
      , hcon_levs       (nlnd,0:sm_lv)             &
      , satcon_levs     (nlnd,0:sm_lv)             &
      , smvccl_levs     (nlnd,sm_lv)               &
      , smvcwt_levs     (nlnd,sm_lv)               &
      , smvcst_levs     (nlnd,sm_lv) )

    ! INGEOFOR
    ALLOCATE                                       &
      ( ug (rw_lng,rw,nobs,nmod_lv)                &
      , vg (rw_lng,rw,nobs,nmod_lv) )


    ! RADCLOUD
    ALLOCATE                                       &
      ( iccb_rad        (rw_lng,rw)                &
      , icct_rad        (rw_lng,rw)                &
      , cca_rad         (rw_lng,rw)                &
      , ccwpin_rad      (rw_lng,rw)                &
      , qcl_rad         (rw_lng,rw,nwet_lv)        &
      , qcf_rad         (rw_lng,rw,nwet_lv)        &
      , layer_cloud_rad (rw_lng,rw,nwet_lv) )

    IF (lhook) CALL dr_hook('ALLOC_FORCING',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_forcing

!=============================================================================

  SUBROUTINE dealloc_forcing

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_FORCING',zhook_in,zhook_handle)


    ! RADCLOUD
    DEALLOCATE                                                                &
      ( layer_cloud_rad,   qcf_rad,          qcl_rad,         ccwpin_rad      &
      , cca_rad,           icct_rad,         iccb_rad )

    ! INGEOFOR
    DEALLOCATE                                                                &
      ( vg,                ug )

    ! INJULES
    DEALLOCATE                                                                &
      ( lw_down,           g_leaf_acc,       g_leaf_phen_acc, resp_w_ft_acc   &
      , cs,                resp_s_acc,       canopy,          catch           &
      , tstar_tile,        snow_tile,        infil_tile,      rgrain          &
      , z0_tile,           z0h_tile,         lai,             canht           &
      , sth,               smcli,            frac_typ,        npp_ft_acc      &
      , frac_disturb,      fsmc,             gs                               &
      , lake_depth,        lake_fetch,       lake_t_mean,     lake_t_mxl      &
      , lake_t_ice,        lake_h_mxl,       lake_h_ice,      lake_shape      &
      , lake_g_dt                                                             &
      , clapp_levs,        sathh_levs,       hcap_levs,       hcon_levs       &
      , satcon_levs,       smvccl_levs,      smvcwt_levs,     smvcst_levs )

    ! LOGIC
    DEALLOCATE                                                                &
      ( cumulus,           soil_mask,        land_ice_mask,   land_sea_mask )

    ! INOBSFOR
    DEALLOCATE                                                                &
      ( tstar_forcing,     flux_e,           flux_h,          w_bg            &
      , w_inc,             v_bg,             u_bg,            t_bg            &
      , v_inc,             u_inc,            t_inc,           q_bg            &
      , q_star )

    ! INPROF
    DEALLOCATE                                                                &
      ( free_tracers,      w_advi,           p_in,            wi              &
      , vi,                ui,               theta,           qi              &
      , ice_fract,         ho2r2_orog,       sil_orog_land,   v_0             &
      , u_0,               di,               z0h_scm,         z0m_scm         &
      , z0mseai,           tstari,           snodepi,         ccai            &
      , iccti,             iccbi,            smci,            canopy_gbi      &
      , t_deep_soili,      i_snowdepth,      i_snow_grnd,     i_rho_snow_grnd &
      , nsnow,             i_tsnowlayer,     i_sice,          i_sliq          &
      , i_ds,              i_rgrainl )


    ! RUNDATA
    DEALLOCATE                                                                &
      ( co2,               nitr_diss,        nitr_acc,        ocff_cld        &
      , ocff_aged,         ocff_new,         soot_cld,        soot_aged       &
      , soot_new,          bmass_cld,        bmass_aged,      bmass_new       &
      , aerosol_em,        aerosol,          dms,             nh3             &
      , so4_diss,          so4_accu,         so4_aitken,      so2             &
      , ozone,             soot,             ocff_hilem,      bmass_hilem     &
      , soot_hilem,        ocff_em,          bmass_em,        soot_em         &
      , dms_em,            nh3_em,           so2_em,          sum_moist_flux  &
      , sum_eng_fluxes,    tscrndcl_ssi,     tstbtrans,       zh              &
      , ddmfx,             sice_alb,         land_alb,        tstar_sice      &
      , tstar_sea,         tstar_land,       fland_ctile,     dOLR_rts        &
      , co2flux,           co2_emits,        orog,            cclwp           &
      , corw,              corvn,            cord,            cort            &
      , ntdsc,             nbdsc,            ntml,            tscrndcl_tile   &
      , albsoil,           albobs_sw,        albobs_vis,      albobs_nir )

    ! INDATA
    DEALLOCATE                                                                &
      ( long,              lat,              gridbox_area,    soil_type )


    IF (lhook) CALL dr_hook('DEALLOC_FORCING',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_forcing

!=============================================================================
END MODULE s_main_force
