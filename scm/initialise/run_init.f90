! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE Run_Init
!
! Purpose:
!   Called by scm_main (Single Column Model main routine) to do
!   initialisations.
!
! Code Description:
!    Language - FORTRAN 90
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================
!     OPTIONS TO SET INITIAL PROFILES
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================

SUBROUTINE run_init                                                           &
  ! (In)
  ( row_length, rows, model_levels, nwet, land_pts, nfor, nbl_levs            &
  , nsoilt_levs, nsoilm_levs, ntiles, nice, nice_use, ntrop, n_cca_lev, ntab  &
  , nprimvars, dayno_init, can_model, ichgf                                   &
  , sec_day, rhcrit, l_mixing_ratio                                           &
  , lcal360, l_snow_albedo                                                    &
  , l_use_sulpc_direct, l_use_seasalt_direct, l_use_soot_direct               &
  , l_use_biogenic, l_use_dust, l_use_bmass_direct, l_use_ocff_direct         &
  , l_use_nitrate_direct, l_use_arclbiom, l_use_arclblck, l_use_arclsslt      &
  , l_use_arclsulp, l_use_arcldust, l_use_arclocff, l_use_arcldlta            &
  , l_murk_rad, l_climat_aerosol, l_use_aod                                   &
  ! (InOut)
  , rho, ti, smcl                                                             &
  ! (Out)
  , iv, iy, idum, iseed, dayno_wint, iccb, icct, cca, cf, cfl, cff, t, q, qcl &
  , qcf, theta_star, pstar, tstar, u, v, w, w_adv, z0msea, zh, t_deep_soil    &
  , canopy_gb, tsi, smc, sthf, sthu, snodep, catch, infil_tile, z0_tile       &
  , z0h_tile_bare, catch_snow, exner_rho_levels, exner_theta_levels, deltap, p&
  , p_theta_levels, rp, rp_theta, ch_tls, ch_qls, ch_uls, ch_vls, ch_wls      &
  , ch_flux_e, ch_flux_h, ch_tstar_forcing                                    &
  , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg, ch_ug, ch_vg                 &
  , flux_e_scm, flux_h_scm, t_inc_scm                                         &
  , u_inc_scm, v_inc_scm, w_inc_scm, q_star_scm, ug_scm, vg_scm               &
  , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                          &
  , tls, qls, uls, vls, wls                                                   &
  , resdump, dap1, dap2, dap3, dab1, dab2, dab3                               &
  , b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt, v_crit, atime            &
  , btime, alfada, dbara, dgrada, pa, tbara, tgrada, tsda, vnbara, vnsda      &
  , vpbara, wbara, wsda, alfadb, dbarb, dgradb, pb, tbarb, tgradb, tsdb       &
  , vnbarb, vnsdb, vpbarb, wbarb, wsdb )


  USE nstypes
  USE conversions_mod,     ONLY: pi
  USE water_constants_mod, ONLY: rho_water
  USE atmos_constants_mod, ONLY: r, cp, kappa, pref
  USE earth_constants_mod, ONLY: g
  
  USE level_heights_mod, ONLY: r_theta_levels

  USE s_main_force, ONLY:                                                     &
     year_init, tapeday_init, soil_type, iccbi, iccti, nml_inprof_thetal      &
   , smi_opt, ndayin, resdump_days, runno_in, runno_out, timestep             &
   , ug, vg, canopy_gbi, ccai, smci, snodepi, t_deep_soili, tstari, ui, vi    &
   , wi, z0mseai, frac_typ, canht, lai, fsmc, smcli, sth, altdat, geoforce    &
   , geoinit, noforce, obs, obs_surf, stats, tapein, tapeout, prindump_obs    &
   , exname_in, exname_out, t_bg, q_bg, u_bg, v_bg, w_bg                      &
   , theta, p_in, flux_e, flux_h, tstar_forcing, qi, t_inc, q_star            &
   , u_inc, v_inc, w_inc, lat, long, rgrain                                   &
   , snow_tile, i_snow_grnd, nsnow, i_rgrainl, i_rho_snow_grnd, i_sice        &
   , i_sliq, i_snowdepth, i_ds, i_tsnowlayer                                  &
   , scm_clapp_levs      => clapp_levs        &
   , scm_sathh_levs      => sathh_levs        &
   , scm_hcap_levs       => hcap_levs         &
   , scm_hcon_levs       => hcon_levs         &
   , scm_satcon_levs     => satcon_levs       &
   , scm_smvccl_levs     => smvccl_levs       &
   , scm_smvcwt_levs     => smvcwt_levs       &
   , scm_smvcst_levs     => smvcst_levs


  USE sw_control_struct
  USE lw_control_struct
  USE spec_sw_lw
  USE mcica_mod
  USE rad_pcf, ONLY: ip_cloud_mcica
  USE switches,   ONLY: l_aggregate

  USE soil_param, ONLY: dzsoil_jules => dzsoil

  ! These parameters are dimensioned as (land_field,sm_levels or 0:sm_levels)
  ! in jules_init.
  USE jules_mod,  ONLY:                      &
      jules_clapp_levs  => clapp_levs        &
    , jules_sathh_levs  => sathh_levs        &
    , jules_hcap_levs   => hcap_levs         &
    , jules_hcon_levs   => hcon_levs         &
    , jules_satcon_levs => satcon_levs       &
    , jules_smvccl_levs => smvccl_levs       &
    , jules_smvcwt_levs => smvcwt_levs       &
    , jules_smvcst_levs => smvcst_levs

  USE scm_utils, ONLY:                                                        &
      dzsoil, old_nml, zhook_in, zhook_out, jprb, lhook, dr_hook

  USE ereport_mod, ONLY : ereport
  USE Control_Max_Sizes
  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(In)
!-----------------------------------------------------------------------------
  INTEGER, INTENT(In) :: &
    row_length           &! Leading x dimension of SCM arrays.
  , rows                 &! Leading y dimension of SCM arrays.
  , model_levels         &! No of levels.
  , nwet                 &! No of model levels in which Q is set.
  , land_pts             &! Number of land points to be processed.
  , nfor                 &! Number terms for observational forcing
  , nbl_levs             &! Number of Boundary layer levels
  , nsoilt_levs          &! Number of soil temperature levels
  , nsoilm_levs          &! Number of soil moisture levels
  , ntiles               &! Number of surface tiles
  , nice                 &! Number of sea-ice categories
  , nice_use             &! Number of sea-ice categories used fully
  , ntrop                &! Max number of levels in the troposphere
  , n_cca_lev            &! No of levels for cca
  , ntab                 &! Dimension of array used in random generator.
  , nprimvars             ! Minimum no. of variables required
                          ! to restart from a dump.



!     Comdecks
!
!     SCM specific :

! Surface scheme.
! Soil Parameters for land surface
! Description:
!   This deck defines the soil parameters used in the SCM.
!
!
!
! Declarations:
  INTEGER, PARAMETER :: &
    nsoilp = 4           ! Number of possible soil parameters

  REAL ::               &
    b_exp_typ(nsoilp)   &! Single Layer :
                         !   (C_EAG in code) Eagleson's exponent for
                         !   calc. sub surf. runoff, P253.4
                         ! Multilayer Hydrology:
                         !   Exponent used in calculation of soil
                         !   water suction and hydraulic conductivity
                         !   (known as B_WAG in HYDROL2A)
                         ! JULES:
                         !   Exponent used in calculation of soil
                         !   water suction and hydraulic conductivity
                         !   (known as B Clapp-Hornberger exponent)

  , sathh_typ(nsoilp)   &! Single layer :
                         !     Dummy
                         ! JULES/Multilayer hydrology :
                         !     Saturated soil water suction

  , satcon_typ(nsoilp)  &! Saturated hydrological conductivity of
                         ! the soil (kg/m2/s)

  , hcap_typ(nsoilp)    &! Soil heat capacity (J/K/m^3)
  , hcon_typ(nsoilp)    &! Soil thermal conductivity (W/m/K)
  , v_crit_typ(nsoilp)  &! Volumetric soil moisture content the critical point
                         ! below this value evaporation falls below its max
                         ! (m^3 per m^3 soil)

  , v_wilt_typ(nsoilp)  &! Volumetric soil moisture content at wilting point
                         ! (m^3/m^3)

  , v_sat_typ(nsoilp)    ! Volumetric soil moisture content at saturation
                         ! (m^3/m^3 soil)


!    (i)  Values of SATCON, V_SAT, B_WAG (B_EXP), SATHH are given by
!         Wageningen for the Maulem/Van Genuchten curves:
!          SATHH = 1 / ALPHA
!          B_WAG = (1-M) / M = 1 / (N-1)
!    (ii) Values of V_CRIT , V_SAT, V_WILT are used which apply
!         for the soil moisture variable defined as the toal soil
!         moisture minus the residuals given by Wageningen.
!    (iii)Soil types ICE, CLAY, LOAM, LOAMY SAND

      COMMON/SOIL_DATA/                                                       &
        b_exp_typ, hcap_typ, hcon_typ, satcon_typ, sathh_typ, v_crit_typ,     &
        v_sat_typ, v_wilt_typ

!-----------------------------------------------------------------------------
                                   ! Unsure if this is still required
                                   ! leave for now
!     Others :
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end


!-----------------------------------------------------------------------------

  INTEGER, INTENT(In) ::     &
    dayno_init               &! Initial day in year
  , can_model                &!
  , ichgf                     ! No. of timesteps between change
                              ! in observational forcing

  INTEGER, INTENT(In) ::     &
    sec_day

  REAL, INTENT(In) ::        &
    rhcrit(nwet)              ! Critical humidity for cloud formation.

! COMMENTS REFER TO TRUE STATUS
!------------------------------
  LOGICAL, INTENT(In) ::     &
    l_mixing_ratio           &! Moisture variables are in mixing ratios
  , lcal360                  &! Use 360 day year
  , l_snow_albedo

  ! Use Direct Effects due to ...
  LOGICAL, INTENT(In) ::     &
    l_use_sulpc_direct       &! ... Sulphate Aerosols
  , l_use_seasalt_direct     &! ... SeaSalt
  , l_use_soot_direct        &! ... Soot
  , l_use_biogenic           &! ... Biogenic Aerosol
  , l_use_dust               &! ... Mineral Dust
  , l_use_bmass_direct       &! ... Biomass
  , l_use_ocff_direct        &! ... Fossil-Fuel Organic Carbon
  , l_use_nitrate_direct      ! ... Nitrate aerosols

  ! From NWP climatology, include Direct Effects due to ...
  LOGICAL, INTENT(In) ::     &
    l_use_arclbiom           &! ... Biomass
  , l_use_arclblck           &! ... Black carbon
  , l_use_arclsslt           &! ... Seasalt
  , l_use_arclsulp           &! ... Sulphate aerosols
  , l_use_arcldust           &! ... Mineral dust
  , l_use_arclocff           &! ... Fossil-Fuel Organic Carbon
  , l_use_arcldlta            ! ... Delta aerosol

  LOGICAL, INTENT(In) ::     &
    l_murk_rad               &! Include Mesoscale Model Aerosols
  , l_climat_aerosol         &! Include Climatological Aerosols
  , l_use_aod                 ! >= 1 of aerosol optical depth diags requested


  INTEGER ::                          &
    land_index      (row_length*rows) &
  , land_ice_index  (row_length*rows) &
  , land_soil_index (row_length*rows)

  REAL, INTENT(Out) ::                       &
    t_bg_scm(row_length,rows,model_levels)   &
  , u_bg_scm(row_length,rows,model_levels)   &
  , v_bg_scm(row_length,rows,model_levels)   &
  , w_bg_scm(row_length,rows,0:model_levels) &
  , q_bg_scm(row_length,rows,nwet)

!-----------------------------------------------------------------------------
! Arguments with INTENT(InOut)
!-----------------------------------------------------------------------------

  REAL, INTENT(InOut) ::                 &
    rho(row_length,rows,model_levels)    &!
  , smcl(land_pts,nsoilm_levs)            ! Soil moisture content
                                          ! in layers (Kg/m2)

  REAL, INTENT(InOut) ::                 &
    ti(row_length,rows,model_levels)      ! Initial temperature profile (K)

!-----------------------------------------------------------------------------
! Arguments with INTENT(Out)
!-----------------------------------------------------------------------------
  ! Random generator variables
  INTEGER, INTENT(Out) :: &
    iv(ntab)              &
  , iy                    &
  , idum                  &
  , iseed                  ! Seed for random number generator

  INTEGER, INTENT(Out) :: &
    dayno_wint             ! Day number relative to winter solstice

  INTEGER, INTENT(Out) :: &
    iccb(row_length,rows) &! Convective cloud base
  , icct(row_length,rows)  ! Convective cloud top

  REAL, INTENT(Out) ::                        &
    cca(row_length,rows,n_cca_lev)            &! Convective cloud amount
  , cf(row_length,rows,nwet)                  &! layer cloud amount
                                               ! (decimal fraction)
  , cfl(row_length,rows,nwet)                 &! liquid layer cloud amount
  , cff(row_length,rows,nwet)                  ! frozen layer cloud amount

  REAL, INTENT(Out) ::                        &
    t(row_length,rows,model_levels)           &! Temperature(K)
  , q(row_length,rows,nwet)                   &! Specific humidity (kg/kg)
  , qcl(row_length,rows,nwet)                 &! Cloud water content(kg/kg)
  , qcf(row_length,rows,nwet)                 &! Cloud ice content (kg/kg)
  , theta_star(row_length,rows,model_levels)  &! Potential temp. increments(K)
  , pstar(row_length,rows)                    &! Surface pressure (Pa)
  , tstar(row_length,rows)                     ! Surface temperature (K)

  REAL, INTENT(Out) ::                        &
    u(row_length,rows,model_levels)           &! Zonal wind (m/s)
  , v(row_length,rows,model_levels)           &! Meridional wind (m/s)
  , w(row_length,rows,0:model_levels)

  REAL, INTENT(Out) ::                        &
    w_adv(row_length,rows,0:model_levels)     &
  , z0msea(row_length,rows)                   &! Sea surface roughness length
  , zh(row_length, rows)                       ! Height above surface of top
                                               ! of boundary layer (m)

  REAL, INTENT(Out) :: &
    t_deep_soil(land_pts,nsoilt_levs)   ! Deep soil temperatures (K)
                                        ! top level not included,=surface

  REAL, INTENT(Out) ::          &
    canopy_gb(land_pts)         &! Canopy water content (kg/m2)
  , tsi(row_length,rows)        &! Temperature of sea-ice
  , smc(land_pts)               &! Soil moisture content(Kg/m^2)
  , sthf(land_pts,nsoilm_levs)  &! Frozen soil moisture content of each
                                 ! layer as a fraction of saturation (kg/m^2)
  , sthu(land_pts,nsoilm_levs)  &! Unfrozen soil moisture content of each
                                 ! layer as a fraction of saturation (kg/m^2)
  , snodep(row_length,rows)      ! Snow depth (kg/m^2)

  REAL, INTENT(Out) ::                 &
    catch(row_length*rows,ntiles)      &! Surface/canopy water capacity of
                                        ! snow-free land tiles (kg/m2)
  , infil_tile(row_length*rows,ntiles) &! Maximum surface infiltration rate
                                        ! for each tile (kg/m2/s)
  , z0_tile(row_length*rows,ntiles)    &! Roughness length for each tile (m)
  , z0h_tile_bare(row_length*rows,ntiles)&! Thermal roughness length 
                                        ! for each tile (m)
  , catch_snow(row_length*rows)         ! Snow capacity for NLT tile (kg/m^2)

  REAL, INTENT(Out) ::                                &
    exner_rho_levels(row_length,rows,model_levels+1)  &
  , exner_theta_levels(row_length,rows,model_levels)  &
  , deltap(row_length,rows,model_levels)              &! Layer Thickness
  , p(row_length,rows,model_levels+1)                 &! Pressure rho levels
  , p_theta_levels(row_length,rows,model_levels)      &! Pressure theta levels
  , rp(row_length,rows,model_levels+1)                &! 1/p on rho levels
  , rp_theta(row_length,rows,model_levels)             ! 1/p on theta levels


  ! Rate of change of increments due to large-scale horizontal and
  ! vertical advection (per second) of ...
  REAL, INTENT(Out) ::                            &
    ch_tls(row_length,rows,nfor-1,model_levels)   &! ... Temp increment
  , ch_qls(row_length,rows,nfor-1,nwet)           &! ... Specific humidity
  , ch_uls(row_length,rows,nfor-1,model_levels)   &! ... Zonal  wind
  , ch_vls(row_length,rows,nfor-1,model_levels)   &! ... Merid. wind
  , ch_wls(row_length,rows,nfor-1,0:model_levels)  ! ... Vert.  wind


  ! Rate of change of background fields for ...
  REAL, INTENT(Out) ::                            &
    ch_t_bg(row_length,rows,nfor-1,model_levels)  &! ... Temperature
  , ch_q_bg(row_length,rows,nfor-1,nwet)          &! ... Specific humidity
  , ch_u_bg(row_length,rows,nfor-1,model_levels)  &! ... Zonal  wind
  , ch_v_bg(row_length,rows,nfor-1,model_levels)  &! ... Merid. wind
  , ch_w_bg(row_length,rows,nfor-1,0:model_levels) ! ... Vert.  wind


  ! Rate of change of geostrophic winds ...
  REAL, INTENT(Out) ::                            &
    ch_ug(row_length,rows,nfor-1,model_levels)    &! ... Zonal  wind
  , ch_vg(row_length,rows,nfor-1,model_levels)     ! ... Merid. wind


  ! Rate of change of forcing of surface ...
  REAL, INTENT(Out) ::                           &
    ch_flux_e(row_length,rows,nfor-1)            &! ... flux_e
  , ch_flux_h(row_length,rows,nfor-1)            &! ... flux_h
  , ch_tstar_forcing(row_length,rows,nfor-1)      ! ... temperature


  REAL, INTENT(Out) ::             &
    flux_e_scm(row_length,rows)    &
  , flux_h_scm(row_length,rows)


  ! Forcing increments (X/timestep)
  REAL, INTENT(Out) ::                          &
    q_star_scm (row_length,rows,nwet)           &! Spec. humid. inc.
                                                 ! (Kg/Kg)/timestep
  , t_inc_scm  (row_length,rows,model_levels)   &! Temp. inc.  (K/timestep)
  , u_inc_scm  (row_length,rows,model_levels)   &! Zonal  wind (m/s)/timestep
  , v_inc_scm  (row_length,rows,model_levels)   &! Merid. wind (m/s)/timestep
  , w_inc_scm  (row_length,rows,0:model_levels)  ! Vert.  wind (m/s)/timestep


  ! Geostrophics wind forcing
  REAL, INTENT(Out) ::                          &
    ug_scm     (row_length,rows,model_levels)   &! Geo. zonal  wind (m/s)
  , vg_scm     (row_length,rows,model_levels)    ! Geo. merid. wind (m/s)


  ! Current LS forcing tendencies (X/day)
  REAL, INTENT(Out) ::                   &
    qls (row_length,rows,nwet)           &! Spec. humid. inc. (Kg/Kg)/day
  , tls (row_length,rows,model_levels)   &! Temp. inc.  (K/day)
  , uls (row_length,rows,model_levels)   &! Zonal  wind (m/s)/day
  , vls (row_length,rows,model_levels)   &! Merid. wind (m/s)/day
  , wls (row_length,rows,0:model_levels)  ! Vert.  wind (m/s)/day


  REAL, INTENT(Out) ::    &
    resdump(row_length,rows,nprimvars)  ! DUMP array of restart variables

!---------------------------------------------------------------------
! Large scale observational forcing
!---------------------------------------------------------------------

! Variables for diagnostic output for observational forcing

  REAL, INTENT(Out) ::                           &! These don't appear
    dap1(row_length,rows,36,model_levels)        &! to do anything
  , dap2(row_length,rows,36,model_levels)        &! except get initialised
  , dap3(row_length,rows,36,nfor-1,model_levels) &
  , dab1(row_length,rows,44)                     &
  , dab2(row_length,rows,44)                     &
  , dab3(row_length,rows,44,nfor-1)

  REAL, INTENT(Out) ::    &
    b_exp  (land_pts)     &! Clapp-Hornberger exponent
  , hcap   (land_pts)     &! Soil heat capacity
  , hcon   (land_pts)     &! Soil thermal conductivity
  , satcon (land_pts)     &! Saturated hydrological conductivity
  , sathh  (land_pts)     &! Saturated soil water suction
  , v_sat  (land_pts)     &! Vol. soil moisture content at saturation
  , v_wilt (land_pts)     &! Vol. soil moisture content at wilting point
  , v_crit (land_pts)      ! Vol. soil moisture content at critical point

!---------------------------------------------------------------------
!     Large scale statistical forcing
!---------------------------------------------------------------------

! Variable for statistical forcing

  REAL, INTENT(Out) ::                   &
    atime                                &! Constants for calculating annual
  , btime                                 ! cycle

  ! Amplitude of seasonal variation of ...
  REAL, INTENT(Out) ::                   &
    alfada(row_length,rows)              &! ... Tuning Factor
  , dbara(row_length,rows,nwet)          &! ... Mean Dew pt. dep. (K)
  , dgrada(row_length,rows,nwet)         &! ... Dew pt. dep. gradient (K/km)
  , pa(row_length,rows, model_levels+1)  &! ... Pressure
  , tbara(row_length,rows,model_levels)  &! ... Temperature (K)
  , tgrada(row_length,rows,model_levels) &! ... Temperature gradient (K/km)
  , tsda(row_length,rows,model_levels)   &! ... SD of temperature (K)
  , vnbara(row_length,rows,model_levels) &! ... Velocity VN (m/s)
  , vnsda(row_length,rows,model_levels)  &! ... SD of velocity VN (m/s)
  , vpbara(row_length,rows,model_levels) &! ... Velocity VP (m/s)
  , wbara(row_length,rows,ntrop)         &! ... Vert. Vel. (mb or HPa/s)
  , wsda(row_length,rows,ntrop)           ! ... SD of Vert. Vel. (mb/s)

  ! Mean of seasonal variation of ...
  REAL, INTENT(Out) ::                   &
    alfadb(row_length,rows)              &! ... Tuning Factor
  , dbarb(row_length,rows,nwet)          &! ... Mean Dew pt. dep (K)
  , dgradb(row_length,rows,nwet)         &! ... Dew pt. dep. gradient (K/km)
  , pb(row_length,rows, model_levels+1)  &! ... Pressure
  , tbarb(row_length,rows,model_levels)  &! ... Temperature (K)
  , tgradb(row_length,rows,model_levels) &! ... Temperature gradient (K/km)
  , tsdb(row_length,rows,model_levels)   &! ... SD of temperature (K)
  , vnbarb(row_length,rows,model_levels) &! ... Velocity VN (m/s)
  , vnsdb(row_length,rows,model_levels)  &! ... SD of velocity VN (m/s)
  , vpbarb(row_length,rows,model_levels) &! ... Velocity VP (m/s)
  , wbarb(row_length,rows,ntrop)         &! ... Vert. Vel. (mb or HPa/s)
  , wsdb(row_length,rows,ntrop)           ! ... SD of Vert. Vel. (mb/s)



!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------
  CHARACTER(len=8), PARAMETER :: routinename = 'run_init'
  CHARACTER(len=80)           :: cmessage ! Error message if Icode >0
  INTEGER                     :: icode

  INTEGER, PARAMETER ::        &
    smc_spec = 0               &! Initial smcl specifed via namelist
  , smc_fsmc = 1               &! Initial smcl calculated via fsmc
  , smc_sth  = 2                ! Initial smcl calculated via sth

  INTEGER ::                   &
    i, j, k, l, m              &! Loop counters
  , runno                      &! Run number of experiment on tape
  , nresdump                   &! Number of restart dumps
  , andayy1, andayy2           &! To calculate andday using time2sec
  , tapeday                    &! Tape year day
  , tapedump_no                &! Number of dumps on tape
  , dummy

  INTEGER ::                   &
    ErrorStatus

  INTEGER ::                   &
    ntml_tmp(row_length,rows)

  INTEGER ::                   &
    tile_pts(ntype)            &! Number of land points which
                                ! include the nth surface type
  , tile_index(land_pts,ntype)  ! Indices of land points which
                                ! include the nth surface type

  INTEGER ::                   &
    file_start                  ! Start point of spectral file name

  REAL ::                  &
    andayy                 &! Number of days in 1 year. (for one year effects)
  , rccb(row_length,rows)  &! Conv. cloud base (real values) for DUMP purposes
  , rcct(row_length,rows)   ! Conv. cloud top  (real values) for DUMP purposes

  REAL ::                  &
    fsmc_min(nsoilp)       &
  , fsmc_max(nsoilp)

  REAL ::                  &
    thl_to_sl              &! Top of mixed layer difference (K)
  , sl_bl                   ! Mixed-layer value for s_L (K)


  REAL ::  &
    factor &
  , tstpfd  ! Timestep fraction of day

  LOGICAL ::               &
    cumulus_tmp(row_length,rows)

  CHARACTER(len=8) ::      &
    exname                  ! Name of experiment on tape

!---------------------------------------------------------------------
!     Define site specific soil parameters and Initialise SWNOCZ
!---------------------------------------------------------------------

  INTEGER ::               &
    path_end

  CHARACTER (LEN=200) ::   &
    mcica_data              ! Path to McICA data file

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('RUN_INIT',zhook_in,zhook_handle)

  tstpfd = timestep/sec_day


  ! Note: jules_init does not currently use land_ice_index
  !       or land_soil_index they are only in the argument list
  !       because an include file was used in jules_init which 
  !       specified them. Jules_init expects them to be dimensioned
  !       as an array of at least size 1
  ! 

  IF (land_pts == 1) land_index(:) = 1
  land_ice_index (:) = 0         ! Not used
  land_soil_index(:) = 0         ! Not used

  ! The routine jules_init is part of the JULES repository. In the
  ! full um it sets the following variables from atm_fields_mod which
  ! the SCM doesn't have. So default values based on soil_type of
  ! the surface of SCM column is used to determine default values
  ! of the following JULEs parameters

  ! The include file s_soilpt which holds initialisation soil data
  ! (via blkdata.F90) is only one level, so all soil levels are
  ! initialised with the same value.

  ! DEPENDS ON: jules_init
  CALL jules_init                                                         &
     ( land_index, land_ice_index, land_soil_index                        &
     , land_pts, ntiles, nsoilm_levs, nice, nice_use, frac_typ            &
     , l_snow_albedo, rgrain, snow_tile, t_deep_soili, i_snow_grnd        &
     , nsnow, i_rgrainl, i_rho_snow_grnd, i_sice, i_sliq, i_snowdepth     &
     , i_ds, i_tsnowlayer )

  IF (land_pts >= 1) THEN
    DO i=1, land_pts
      jules_clapp_levs(i,:)  = b_exp_typ  (soil_type(i))
      jules_sathh_levs(i,:)  = sathh_typ  (soil_type(i))
      jules_hcap_levs(i,:)   = hcap_typ   (soil_type(i))
      jules_hcon_levs(i,:)   = hcon_typ   (soil_type(i))
      jules_satcon_levs(i,:) = satcon_typ (soil_type(i))
      jules_smvccl_levs(i,:) = v_crit_typ (soil_type(i))
      jules_smvcwt_levs(i,:) = v_wilt_typ (soil_type(i))
      jules_smvcst_levs(i,:) = v_sat_typ  (soil_type(i))
    END DO

    ! Jules has set the defaults, now see if SCM user has overrriden them 
    ! in the SCM forcing namelists, test on first element, if unset
    ! it should be rmdi
    IF (scm_clapp_levs  (1,1) /= rmdi)                                   &
                        jules_clapp_levs  (:,:) = scm_clapp_levs  (:,:)

    IF (scm_sathh_levs  (1,1) /= rmdi)                                   &
                        jules_sathh_levs  (:,:) = scm_sathh_levs  (:,:)

    IF (scm_hcap_levs   (1,1) /= rmdi)                                   &
                        jules_hcap_levs   (:,:) = scm_hcap_levs   (:,:)

    IF (scm_hcon_levs   (1,1) /= rmdi)                                   &
                        jules_hcon_levs   (:,:) = scm_hcon_levs   (:,:)

    IF (scm_satcon_levs (1,1) /= rmdi)                                   &
                        jules_satcon_levs (:,:) = scm_satcon_levs (:,:)

    IF (scm_smvccl_levs (1,1) /= rmdi)                                   &
                        jules_smvccl_levs (:,:) = scm_smvccl_levs (:,:)

    IF (scm_smvcwt_levs (1,1) /= rmdi)                                   &
                        jules_smvcwt_levs (:,:) = scm_smvcwt_levs (:,:)

    IF (scm_smvcst_levs (1,1) /= rmdi)                                   &
                        jules_smvcst_levs (:,:) = scm_smvcst_levs (:,:)


    ! Set surface values passed to atmos_phyics2 via argument lists
    DO i=1, land_pts
      b_exp  (i) = jules_clapp_levs  (i,1)
      sathh  (i) = jules_sathh_levs  (i,1)
      hcap   (i) = jules_hcap_levs   (i,1)
      hcon   (i) = jules_hcon_levs   (i,0)
      satcon (i) = jules_satcon_levs (i,0)
      v_crit (i) = jules_smvccl_levs (i,1)
      v_wilt (i) = jules_smvcwt_levs (i,1)
      v_sat  (i) = jules_smvcst_levs (i,1)
    END DO
  END IF ! Test on land pts

!---------------------------------------------------------------------
!     Calculate number of days in year
!---------------------------------------------------------------------

! DEPENDS ON: time2sec
  CALL time2sec (year_init, 1, 1, 0, 0, 0, 0, 0, andayy1, dummy, lcal360)

! DEPENDS ON: time2sec
  CALL time2sec (year_init+1, 1, 1, 0, 0, 0, 0, 0, andayy2, dummy, lcal360)
  andayy = andayy2 - andayy1


  IF (stats) THEN

!---------------------------------------------------------------------
!       Calculate day number relative to winter solstice 
!       (ie day 351 assuming 360 day calendar)
!---------------------------------------------------------------------

    IF (dayno_init  <   351) THEN
      dayno_wint = dayno_init + 9
    ELSE
      dayno_wint = dayno_init - 351
    END IF

!---------------------------------------------------------------------
!       Derive initial data from climate datasets
!---------------------------------------------------------------------
    DO k=1,model_levels+1
      DO j=1, rows
        DO i=1, row_length
          p(i,j,k) = p_in(i,j,k)
        END DO
      END DO
    END DO

! DEPENDS ON: initstat
    CALL initstat                                                             &
      ( row_length, rows, model_levels, nwet, ntrop, andayy, dayno_wint, q, t &
      , lat, long, p_in, pa, pb, alfada, alfadb, tbara, tbarb, tsda, tsdb     &
      , tgrada, tgradb, dbara, dbarb, dgrada, dgradb, vnbara, vnbarb, vnsda   &
      , vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb, atime, btime         &
      , p_theta_levels )

    DO k=1,model_levels
      DO j=1, rows
        DO i=1, row_length
          p_in(i,j,k) = p(i,j,k)
        END DO
      END DO
    END DO
!---------------------------------------------------------------------
!       Initialise random generator
!---------------------------------------------------------------------

! DEPENDS ON: g05cbe
    CALL g05cbe(iseed)
  END IF                     ! stats

!---------------------------------------------------------------------
!     Set initial data from &INPROF
!---------------------------------------------------------------------

  DO j=1, rows
    DO i=1, row_length

      IF (stats .OR. obs .OR. noforce .OR. geoforce) THEN
!  Do for stats as well for now as we have no start data otherwise.
        DO k=1, model_levels
          u(i,j,k) = ui(i,j,k)
          v(i,j,k) = vi(i,j,k)
          w(i,j,k) = wi(i,j,k)
          w_adv(i,j,k) = w(i,j,k)
          t(i,j,k) = ti(i,j,k)
        END DO

        w(i,j,0)     = 0.0
        w_adv(i,j,0) = w(i,j,0)

        DO k = 1, nwet
          q(i,j,k) = qi(i,j,k)
        END DO
      END IF                   ! (obs .or. noforce)

      IF (stats .AND. altdat) THEN
        DO k = 1, model_levels
          t(i,j,k) = ti(i,j,k)
        END DO
        DO k=1, nwet
          q(i,j,k) = qi(i,j,k)
        END DO
      END IF                   ! (stats .and. altdat)

      IF (geoforce .AND. geoinit) THEN
        DO k=1, model_levels
          u(i,j,k) = ug(i,j,1,k)
          v(i,j,k) = vg(i,j,1,k)
        END DO
      END IF                   ! geoforce and geoinit

    END DO
  END DO                     ! i

!---------------------------------------------------------------------
!     Calculate rates of change for large scale observational forcing
!---------------------------------------------------------------------

  IF (obs) THEN

    IF (old_nml) THEN
      ! Reproduce old run format style
      ! Old style only had u_inc, v_inc, and w_inc. and treated
      ! them as background fields or forcing tendencies depending on
      ! whether l_windrlx was .TRUE. or .FALSE.

      ! So make wind _bg = wind _inc fields if relaxation is in old format
      ! so it will reproduce same results in restructured forcing code later
      ! on

      u_bg (:,:,:,:) = u_inc (:,:,:,:)
      v_bg (:,:,:,:) = v_inc (:,:,:,:)
      w_bg (:,:,:,:) = w_inc (:,:,:,:)

      ! This is not the preferred method, users with old namelists
      ! generated from earlier documentation wishing to use
      ! these in the restructured code should really rename u_inc, v_inc
      ! and w_inc to u_bg, v_bg and w_bg with units (m/s)

      ! In addition they should source forcing data for u_inc, v_inc and
      ! w_inc in the units of (m/s)/day

    END IF

    factor = ichgf * timestep

    DO l=1, (nfor-1)
      ch_t_bg(:,:,l,:) = (t_bg(:,:,l+1,:) - t_bg(:,:,l,:)) / factor
      ch_q_bg(:,:,l,:) = (q_bg(:,:,l+1,:) - q_bg(:,:,l,:)) / factor
      ch_u_bg(:,:,l,:) = (u_bg(:,:,l+1,:) - u_bg(:,:,l,:)) / factor
      ch_v_bg(:,:,l,:) = (v_bg(:,:,l+1,:) - v_bg(:,:,l,:)) / factor
      ch_w_bg(:,:,l,:) = (w_bg(:,:,l+1,:) - w_bg(:,:,l,:)) / factor

      ch_tls(:,:,l,:)  = (t_inc  (:,:,l+1,:) - t_inc  (:,:,l,:)) / factor
      ch_qls(:,:,l,:)  = (q_star (:,:,l+1,:) - q_star (:,:,l,:)) / factor
      ch_uls(:,:,l,:)  = (u_inc  (:,:,l+1,:) - u_inc  (:,:,l,:)) / factor
      ch_vls(:,:,l,:)  = (v_inc  (:,:,l+1,:) - v_inc  (:,:,l,:)) / factor
      ch_wls(:,:,l,:)  = (w_inc  (:,:,l+1,:) - w_inc  (:,:,l,:)) / factor
    END DO

    ! Initialise background fields
    t_bg_scm(:,:,:) = t_bg(:,:,1,:)
    q_bg_scm(:,:,:) = q_bg(:,:,1,:)
    u_bg_scm(:,:,:) = u_bg(:,:,1,:)
    v_bg_scm(:,:,:) = v_bg(:,:,1,:)
    w_bg_scm(:,:,:) = w_bg(:,:,1,:)

    ! Initialise LS forcing tendencies
    tls (:,:,:) = t_inc (:,:,1,:)
    qls (:,:,:) = q_star(:,:,1,:)
    uls (:,:,:) = u_inc (:,:,1,:)
    vls (:,:,:) = v_inc (:,:,1,:)
    wls (:,:,:) = w_inc (:,:,1,:)

    ! Set initial increments for timestep
    t_inc_scm (:,:,:) = tstpfd * tls(:,:,:)
    q_star_scm(:,:,:) = tstpfd * qls(:,:,:)
    u_inc_scm (:,:,:) = tstpfd * uls(:,:,:)
    v_inc_scm (:,:,:) = tstpfd * vls(:,:,:)
    w_inc_scm (:,:,:) = tstpfd * wls(:,:,:)


    IF (obs_surf) THEN

      DO l=1, (nfor-1)
        ch_flux_h(:,:,l) = (flux_h(:,:,l+1) - flux_h(:,:,l)) / factor
        ch_flux_e(:,:,l) = (flux_e(:,:,l+1) - flux_e(:,:,l)) / factor
        ch_tstar_forcing(:,:,l) = (tstar_forcing(:,:,l+1)                      &
                                     - tstar_forcing(:,:,l)) / factor
      END DO

      flux_h_scm(:,:) = flux_h(:,:,1)
      flux_e_scm(:,:) = flux_e(:,:,1)
      tstar(:,:)      = tstar_forcing(:,:,1)

    END IF  ! obs_surf

  END IF ! obs


  IF (geoforce) THEN
    DO k = 1, model_levels
      DO l = 1, nfor - 1
        DO j = 1, rows
          DO i = 1, row_length
            ch_ug(i,j,l,k)  = (ug(i,j,l+1,k) - ug(i,j,l,k))         &
                            /             (ichgf * timestep)
            ch_vg(i,j,l,k)  = (vg(i,j,l+1,k) - vg(i,j,l,k))         &
                            /             (ichgf * timestep)
          END DO
        END DO
      END DO
    END DO

    ug_scm (:,:,:) = ug (:,:,1,:)
    vg_scm (:,:,:) = vg (:,:,1,:)

  END IF


  IF (tapein) THEN

!---------------------------------------------------------------------
!       Read tape data if required to carry on from previous run
!---------------------------------------------------------------------

    READ (50) exname,runno,tapedump_no
    WRITE (6,'(A17,/,A30,A6,/,A30,i4,/,A30,i4)')                          &
      ' from tape header',                                                &
      ' expt. name is                ',exname,                            &
      ' run no. is                   ',runno,                             &
      ' no. of dumps on tape are     ',tapedump_no

!---------------------------------------------------------------------
!       Check for correct data set
!---------------------------------------------------------------------

    IF (exname  ==  exname_in .AND. runno  ==  runno_in) THEN

!---------------------------------------------------------------------
!         Look for correct day - tapeday_init input in namelist
!         INDATA.
!---------------------------------------------------------------------

      DO i=1, tapedump_no
        IF (stats) THEN
          READ (50) tapeday, resdump, iv, iy, idum
        ELSE IF (obs) THEN
          READ (50) tapeday, resdump
        END IF
        WRITE (6,'(A16,I4,A22,I4)')                                       &
          ' tape year day= ', tapeday,                                    &
          '      start year day= ', tapeday_init

        IF (tapeday  ==  (tapeday_init-1) ) THEN
          GOTO 9999
        END IF

!           If the end of the tape is reached and the specified day
!           tapeday_init not found o/p error message and stop run.

        IF (i  ==  tapedump_no) THEN
          icode = 520
          WRITE (6,'(A39,A6,A4,I3/)')                             &
            ' Initial day not found on data set m20.',exname,'.run',runno

          
          CALL ereport(routinename, icode, cmessage)
        END IF
      END DO
9999  CLOSE (50)

!---------------------------------------------------------------------
!         Read initial data from tape in DUMP format.
!---------------------------------------------------------------------

! DEPENDS ON: dumpinit
      CALL dumpinit                                                           &
        ( row_length, rows, nprimvars, land_pts, model_levels, nwet, nbl_levs &
        , nsoilt_levs, nsoilm_levs, ntrop, n_cca_lev, resdump, u, v, w, t     &
        , theta, q, qcl, qcf, cf, p, rho, t_deep_soil, smc                    &
        , canopy_gb, snodep, tstar, zh, z0msea, cca, rccb, rcct, smcl )

!         Sometimes (when initial wind in dump is arbitrary) we want
!         to reset the wind to geostrophic. To do this set geoinit to
!         true in logic namelist

      DO j=1, rows
        DO i=1, row_length
          IF (geoinit .AND. geoforce) THEN
            DO k=1, model_levels
              u(i,j,k) = ug(i,j,1,k)
              v(i,j,k) = vg(i,j,1,k)
            END DO
          END IF
          iccb(i,j) = INT(rccb(i,j))
          icct(i,j) = INT(rcct(i,j))
          tsi(i,j) = tstar(i,j)
        END DO
      END DO

!---------------------------------------------------------------------
!         G05CGE restores the state of the basic generator
!         routine G05DDE following the call to G05CFE at the end of
!         each day.
!---------------------------------------------------------------------

! DEPENDS ON: g05cge
      CALL g05cge(idum,iv,iy)
    ELSE
      icode = 521
      WRITE (6,'(A21,A6,A4,I3,A10,/)')                                  &
        'Initial data set m20', exname_in,'.run',runno_in,' not found'

      
      CALL ereport(routinename, icode, cmessage)
    END IF                   ! (exname  ==  exname_in
                          ! .and. runno  ==  runno_in)


  ELSE                      ! not tapein
!---------------------------------------------------------------------
!       Set initial values if no tape data to be used
!---------------------------------------------------------------------
!

    IF (land_pts > 0) THEN
!-----------------------------------------------------------------------
! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
      CALL tilepts (land_pts, frac_typ, tile_pts, tile_index)
!-----------------------------------------------------------------------
! Initialise tiled and gridbox mean vegetation parameters
!-----------------------------------------------------------------------
      WRITE(6,*) 'RUN_INIT: CALLING SPARM'
! DEPENDS ON: sparm
      CALL sparm                                                  &
        ( land_pts, ntiles, can_model, l_aggregate                &
        , tile_pts, tile_index, frac_typ, canht                   &
        , lai, satcon, catch_snow, catch, infil_tile, z0_tile     &
        , z0h_tile_bare )

    END IF

    DO k=1, model_levels+1
      DO j=1, rows
        DO i=1, row_length
          p(i,j,k) = p_in(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          cca(i,j,k) = ccai(i,j)
        END DO
      END DO
    END DO

    DO j=1, rows
      DO i=1, row_length
        tstar  (i,j) = tstari  (i,j)
        tsi    (i,j) = tstar   (i,j)
        z0msea (i,j) = z0mseai (i,j)
        iccb   (i,j) = iccbi   (i,j)
        icct   (i,j) = iccti   (i,j)
        snodep (i,j) = snodepi (i,j) ! Includes snow depth on land/sea
      END DO
    END DO


    IF (land_pts > 0) THEN

    !--------------------------------------------------------------
    ! Initialise the canopy water
    !--------------------------------------------------------------
      DO i=1, land_pts
        IF (canopy_gbi(i) == RMDI) THEN
          icode = 522
          WRITE (6,'(5(TR1,A52/))')                                    &
            '====================================================',    &
            '| Initial Gridbox Mean Canopy Water (canopy_gbi)   |',    &
            '| must be specified for each land point in         |',    &
            '| namelist.                                        |',    &
            '===================================================='

          
          CALL ereport(routinename, icode, cmessage)

        ELSE
          canopy_gb(i) = canopy_gbi(i)
        END IF
      END DO

    !--------------------------------------------------------------
    ! Initialise the soil moisture content (SMC)
    !--------------------------------------------------------------
      DO i=1, land_pts
        IF (smci(i) == RMDI) THEN
          icode =523
          WRITE (6,'(4(TR1,A52/))')                                  &
            '====================================================',  &
            '| Initial Soil Moisture Content (smci) must be     |',  &
            '| specified for each land point in namelist.       |',  &
            '===================================================='
          
          CALL ereport(routinename, icode, cmessage)
        ELSE
          smc(i) = smci(i)
        END IF
      END DO

    !--------------------------------------------------------------
    ! Initialise the deep soil temperature
    !--------------------------------------------------------------
      DO k=1, nsoilm_levs
        DO i=1, land_pts
          IF (t_deep_soili(i,k) == RMDI) THEN
            icode=524
            WRITE (6,'(2(TR1,A52/),(TR1,A15,I2,A18,T53,A1/)'//         &
                     ',2(TR1,A52/))')                                  &
              '====================================================',  &
              '| Initial Deep Soil Temperature (t_deep_soili)     |',  &
              '| must specify ',nsoilt_levs,' soil temperatures','|',  &
              '| for each land point in namelist.                 |',  &
              '===================================================='
            
            CALL ereport(routinename, icode, cmessage)
          ELSE
            t_deep_soil(i,k) = t_deep_soili(i,k)
          END IF
        END DO
      END DO               ! nsoilm_levs

    !--------------------------------------------------------------
    ! Initialise the soil moisture in root zone layers (SMCL)
    !--------------------------------------------------------------

      SELECT CASE (smi_opt)
        CASE(smc_spec)
        ! Soil moisture (smcli) to be specified by user namelist.

        ! Check for missing data
          DO k=1, nsoilm_levs
            DO i=1, land_pts
              IF (smcli(i,k) == RMDI) THEN
                icode=525
                WRITE (6,'(2(TR1,A52/),'//                                &
                         '  (TR1,A25,I2,A21,T53,A1/),'//                  &
                         ' 2(TR1,A52/))')                                 &
                  '====================================================', &
                  '| Initial Soil Moisture Content in layers (smcli)  |', &
                  '| must specify values on ', nsoilm_levs,               &
                                             ' soil moisture levels','|', &
                  '| for each land point in namelist.                 |', &
                  '===================================================='
                
                CALL ereport(routinename, icode, cmessage)
              ELSE
                smcl(i,k) = smcli(i,k)
              END IF
            END DO
          END DO

        CASE(smc_fsmc)
          ! Calculate max & min limits of FSMC
          DO i=1, nsoilp
            fsmc_min(i) = -  v_wilt_typ(i)                                    &
                          / (v_crit_typ(i) - v_wilt_typ(i))
            fsmc_max(i) =   (v_sat_typ(i)  - v_wilt_typ(i))                   &
                          / (v_crit_typ(i) - v_wilt_typ(i))
          END DO

          ! Soil moisture to be initialised by soil stress factor
          ! (FSMC) from namelist.


          DO i=1, land_pts

            ! Check for missing data
            IF (fsmc(i) == RMDI) THEN
              icode=526
              WRITE (6,'(4(TR1,A52/))')                                  &
                '====================================================',  &
                '| Soil moisture stress factor (fsmc) must be       |',  &
                '| specified for each land point in namelist.       |',  &
                '===================================================='
              
              CALL ereport(routinename, icode, cmessage)
            END IF

            ! Check fsmc is within range for given soil type
            IF ((fsmc(i) < fsmc_min(soil_type(i))) .OR.                  &
                (fsmc(i) > fsmc_max(soil_type(i)))) THEN
              icode=527
              WRITE (6,'(5(TR1,A52/),'//                                 &
                       ' 4(TR1,A19,F6.3,A10,F6.3,T53,A1/),'//            &
                       ' 2(TR1,A52/))')                                  &
                '====================================================',  &
                '| Soil moisture stress factor (fsmc) is outside    |',  &
                '| acceptable range. For the specified soil types,  |',  &
                '| fsmc should be in ranges:                        |',  &
                '|                                                  |',  &
                '|  1) Ice        : ',fsmc_min(1),' < fsmc < ',          &
                                      fsmc_max(1),                 '|',  &
                '|  2) Clay       : ',fsmc_min(2),' < fsmc < ',          &
                                      fsmc_max(2),                 '|',  &
                '|  3) Loam       : ',fsmc_min(3),' < fsmc < ',          &
                                      fsmc_max(3),                 '|',  &
                '|  4) Loamy Sand : ',fsmc_min(4),' < fsmc < ',          &
                                      fsmc_max(4),                 '|',  &
                '|                                                  |',  &
                '===================================================='
              
              CALL ereport(routinename, icode, cmessage)
            END IF
          END DO

            DO k=1, nsoilm_levs
              DO i=1, land_pts
                smcl(i,k) = rho_water*dzsoil_jules(k)                    &
                          * ( fsmc(i)*v_crit(i) + (1-fsmc(i))*v_wilt(i) )
              END DO
            END DO

        CASE(smc_sth)
          ! Soil moisture is to be initialised by total soil
          ! moisture as a fraction of saturation (STH).

          DO k=1, nsoilm_levs
            DO i=1, land_pts

              ! Check for missing data
              IF (sth(i,k) == RMDI) THEN
                icode=528
                WRITE (6,'(2(TR1,A52/),'//                                 &
                         '  (TR1,A37,I2,A5,T53,A1/),'//                    &
                         ' 2(TR1,A52/))')                                  &
                  '====================================================',  &
                  '| Total Soil Moisture (sth) as a fraction of       |',  &
                  '| saturation, must specify values on ', nsoilm_levs,    &
                                                             ' soil','|',  &
                  '| moisture levels for each land point in namelist. |',  &
                  '===================================================='
                
                CALL ereport(routinename, icode, cmessage)
              END IF

              ! Check sth range
              IF ((sth(i,k) < 0.0) .OR.                                    &
                  (sth(i,k) > 1.0)) THEN
                icode=529
                WRITE (6,'(5(TR1,A52/))')                                  &
                  '====================================================',  &
                  '| Specified values for Total Soil Moisture (sth)   |',  &
                  '| as a fraction of saturation, must be within the  |',  &
                  '| range: 0.0 ~ 1.0                                 |',  &
                  '===================================================='
                
                CALL ereport(routinename, icode, cmessage)
              END IF

                smcl(i,k) = sth(i,k)*rho_water*v_sat(i) * dzsoil_jules(k)

            END DO
          END DO

        CASE(IMDI)
          icode=530
          WRITE (6,'(5(TR1,A52/))')                                  &
            '====================================================',  &
            '| An initialization option for soil moisture       |',  &
            '| content in layers (smi_opt) must be specified    |',  &
            '| when running over land points                    |',  &
            '===================================================='

          
          CALL ereport(routinename, icode, cmessage)
      END SELECT
    END IF                    ! land_pts

  END IF                     ! tapein

!---------------------------------------------------------------------
!     Initialise the frozen and unfrozen soil moisture.
!---------------------------------------------------------------------
  ! DEPENDS ON: freeze_soil_jules
  CALL freeze_soil_jules                                                    &
    ( land_pts, nsoilm_levs, jules_clapp_levs, dzsoil_jules                 &
    , jules_sathh_levs, smcl, t_deep_soil, jules_smvcst_levs                &
    , sthu, sthf )

!---------------------------------------------------------------------
!     Calculate pressure, exner_theta_levels and pstar
!---------------------------------------------------------------------

! DEPENDS ON: calc_press
  CALL calc_press                                                           &
! In
    ( model_levels, nwet, rows, row_length, p, theta, q                     &
    , .TRUE., .TRUE.                                                        &
! InOut
    , rho                                                                   &
! Out
    , exner_theta_levels, exner_rho_levels, p_theta_levels, rp, rp_theta    &
    , pstar )


!---------------------------------------------------------------------
!     Calculate DELTAP
!---------------------------------------------------------------------

  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        deltap(i,j,k) = p(i,j,k+1)-p(i,j,k)
      END DO
    END DO
  END DO
!---------------------------------------------------------------------
!     Initialise cloud water (QCL,QCF)
!---------------------------------------------------------------------

  IF (.NOT. tapein) THEN
    IF (NML_inprof_thetal == 0) THEN

      ! Standard, to input theta and qi as theta and qv
      DO k=1, nwet
! DEPENDS ON: initqlcf
        CALL initqlcf                                                         &
           ( p_theta_levels(:,:,k), rhcrit, q(:,:,k), t(:,:,k)                &
           , model_levels, row_length, rows, cf(:,:,k)                        &
           , qcf(:,:,k), qcl(:,:,k), nbl_levs, k )
      END DO
      DO i=1, row_length
        DO j=1, rows
          DO  k=1, nwet
            IF ( qcl(i,j,k)+qcf(i,j,k) > 0.0 ) THEN
              ! use simple water-content-dependent split 
              cfl(i,j,k) = cf(i,j,k)*qcl(i,j,k)/(qcl(i,j,k)+qcf(i,j,k))
              cff(i,j,k) = 1.0 - cfl(i,j,k)
            ELSE
              cfl(i,j,k) = 0.0
              cff(i,j,k) = 0.0
            END IF
          END DO  ! k
        END DO  ! j
      END DO  ! i

    ELSE
      !-------------------------------------------------------
      ! Alternative, to input theta and qi as theta_l and qt
      ! where theta_l = (T -L*qcl/cp)/pi = theta - L*qcl/(cp*pi)
      !-------------------------------------------------------
      ! Set up fields ready to pass to ls_cld

      DO j=1, rows
        DO i=1, row_length
          ! Just initialise assuming NOT convection
          cumulus_tmp(i,j) = .FALSE.

          DO k=1, nwet
            qcf(i,j,k) = 0.0
          END DO

          !-------------------------------------------------------
          ! Impose well-mixed SL(=T -L*qcl/cp+g*z/cp) where input
          ! theta_l is well-mixed in the BL
          ! because SL is the BL's conserved thermodynamic variable
          !-------------------------------------------------------
          IF (NML_inprof_thetal == 2) THEN

            DO k=1, nwet
              ! Convert "T" back to "theta"=theta_l from namelist "theta"
              t(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
            END DO

            ! Find where t (currently =theta_l from namelist inprof)
            ! stops being well-mixed.
            k=2
            DO WHILE ( k < nwet .AND. &
                       ABS(t(i,j,k)-t(i,j,1)) < 0.01)
              k=k+1
            END DO

            ntml_tmp(i,j) = k-1

            ! mixed layer s_L = mixed_layer_theta_l * pi(surf)
            sl_bl = t(i,j,1)*(pstar(i,j)/pref)**kappa

            DO k=1, ntml_tmp(i,j)
              ! Save theta_l*pi
              thl_to_sl = t(i,j,k)*exner_theta_levels(i,j,k)

              ! Convert s_L to T_L and
              ! set mixed layer T_L = mixed layer s_L - gz/cp
              t(i,j,k) = sl_bl                                                &
                       - g*(r_theta_levels(i,j,k)-r_theta_levels(i,j,0))/cp

              ! Calc change on converting from mixed theta_l to mixed s_L
              thl_to_sl = t(i,j,k) - thl_to_sl
            END DO

            ! For rest of profile, adjust by same amount as top of mixed
            ! layer, in order to preserve inversion jump
            DO k=ntml_tmp(i,j)+1, nwet
              t(i,j,k) = t(i,j,k)*exner_theta_levels(i,j,k) + thl_to_sl
            END DO

          END IF  ! test on NML_inprof_thetal=2
        END DO
      END DO

      ! Arguments passed to ls_cld:
      !   t = T_L = T-L*qcl/cp
      !   q = q_w = q+qcl
! DEPENDS ON: ls_cld
      CALL ls_cld                                                             &
        ( p_theta_levels(1,1,1), rhcrit, nwet, nbl_levs, row_length, rows     &
        , ntml_tmp, cumulus_tmp, l_mixing_ratio, t, cf, q, qcf, qcl, cfl, cff &
        , ErrorStatus )

      DO k=1, nwet
        DO j=1, rows
          DO i=1, row_length
            ! Reset initial profiles
            ti(i,j,k) = t(i,j,k)
            qi(i,j,k) = q(i,j,k)
          END DO
        END DO
      END DO

    END IF
  END IF

!---------------------------------------------------------------------
!       Convert temperature to potential temperature and t_inc to
!       theta_star
!---------------------------------------------------------------------

! Set T to Theta and t_inc to theta_star
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        theta(i,j,k)      = t(i,j,k)                                          &
                          / exner_theta_levels(i,j,k)

  ! Note: This replicates old code behaviour. This has been identified and
  !       points to a bug that may have always been present.
  !       This bug means that forcings on the 1st timestep are a
  !       factor of 86400/timestep larger that subsequent timesteps.
  !       This is to be fixed at UM8.0, theta_star may be affected in a
  !       similar manner by this bug.

        theta_star(i,j,k) = t_inc(i,j,nfor,k)                                 &
                          / exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

!---------------------------------------------------------------------
!     Zero diagnostics for observational forcing
!---------------------------------------------------------------------

  IF (obs .AND. prindump_obs) THEN
    DO i=1, row_length
      DO m=1, rows
        DO k=1, model_levels
          DO j=1, 36
            dap1(i,m,j,k) = 0.0
            dap2(i,m,j,k) = 0.0
            DO  l = 1, nfor-1
              dap3(i,m,j,k,l) = 0.0
            END DO
          END DO
        END DO
        DO j=1, 44
          dab1(i,m,j) = 0.0
          dab2(i,m,j) = 0.0
          DO k=1, nfor-1
            dab3(i,m,j,k) = 0.0
          END DO
        END DO
      END DO                   ! i
    END DO                    ! m
  END IF                     ! obs


!---------------------------------------------------------------------
!     Radiation - LWLKIN in UM code deck LWTRAN, SWLKIN in UM code
!     SWTRAN.
!     Initialise the LW and SW spectral band tables - standard
!     Radiation code.
!     For Edwards-Slingo radiation code 3C use R2_LW_SPECIN and
!     R2_SW_SPECIN to pick up the spectral files as specified in the
!     control structures (rather than the old method of using units
!     57 (SW) and 80 (LW) ).
!---------------------------------------------------------------------
!     Longwave
  ErrorStatus = 0

  DO j=1,n_lwcall
    file_start=SCAN(lw_control(j)%spectral_file,'/',.TRUE.)

    SELECT CASE ( lw_control(j)%spectral_file(file_start+1:file_start+3) )

    CASE ('ses','SES')
! DEPENDS ON: ses_inilw
      CALL ses_inilw                                                          &
        ( ErrorStatus, Cmessage                                               &
        , lw_control(j)%spectral_file                                         &
        , lw_control(j)%l_ch4, lw_control(j)%l_n2o                            &
        , lw_control(j)%l_cfc11, lw_control(j)%l_cfc12                        &
        , lw_control(j)%l_cfc113, lw_control(j)%l_cfc114                      &
        , lw_control(j)%l_hcfc22, lw_control(j)%l_hfc125                      &
        , lw_control(j)%l_hfc134a                                             &
        , l_climat_aerosol                                                    &
        , l_use_dust, l_use_arcldust                                          &
        , l_use_sulpc_direct, l_use_arclsulp                                  &
        , l_use_soot_direct, l_use_arclblck                                   &
        , l_use_bmass_direct, l_use_arclbiom                                  &
        , l_use_seasalt_direct, l_use_arclsslt                                &
        , l_use_ocff_direct, l_use_arclocff                                   &
        , l_use_biogenic, l_use_arcldlta                                      &
        , l_use_nitrate_direct                                                &
        , l_murk_rad, l_use_aod                                               &
        , lw_control(j)%l_gas, lw_control(j)%l_continuum                      &
        , lw_control(j)%l_drop, lw_control(j)%l_aerosol                       &
        , lw_control(j)%l_ice                                                 &
        , lw_control(j)%i_solar_src, lw_spectrum(j) )

    CASE DEFAULT
! DEPENDS ON: r2_lw_specin
      CALL r2_lw_specin                                                       &
        ( ErrorStatus, Cmessage, lw_control(j)%spectral_file                  &
        , lw_control(j)%l_ch4, lw_control(j)%l_n2o, lw_control(j)%l_cfc11     &
        , lw_control(j)%l_cfc12, lw_control(j)%l_cfc113                       &
        , lw_control(j)%l_hcfc22, lw_control(j)%l_hfc125                      &
        , lw_control(j)%l_hfc134A, l_climat_aerosol, l_use_dust               &
        , l_use_arcldust, l_use_sulpc_direct, l_use_arclsulp                  &
        , l_use_soot_direct, l_use_arclblck, l_use_bmass_direct               &
        , l_use_arclbiom, l_use_seasalt_direct, l_use_arclsslt                &
        , l_use_ocff_direct, l_use_arclocff, l_use_biogenic, l_use_arcldlta   &
        , l_use_nitrate_direct, l_murk_rad, l_use_aod, lw_control(j)%l_gas    &
        , lw_control(j)%l_continuum, lw_control(j)%l_drop                     &
        , lw_control(j)%l_aerosol, lw_control(j)%l_ice, lw_spectrum(j) )
    END SELECT
  END DO

  IF (ErrorStatus  /=  0) THEN
    
    CALL ereport(routinename, errorstatus, cmessage)
  END IF
!     Shortwave
  ErrorStatus = 0

  DO j=1, n_swcall
    file_start=SCAN(sw_control(j)%spectral_file,'/',.TRUE.)

    SELECT CASE ( sw_control(j)%spectral_file(file_start+1:file_start+3) )

    CASE ('ses','SES')
! DEPENDS ON: ses_inisw
      CALL ses_inisw                                                          &
        ( ErrorStatus, Cmessage                                               &
        , sw_control(j)%spectral_file                                         &
        , sw_control(j)%l_o2                                                  &
        , sw_control(j)%l_ch4                                                 &
        , sw_control(j)%l_n2o                                                 &
        , l_climat_aerosol                                                    &
        , l_use_dust, l_use_arcldust                                          &
        , l_use_sulpc_direct, l_use_arclsulp                                  &
        , l_use_soot_direct, l_use_arclblck                                   &
        , l_use_bmass_direct, l_use_arclbiom                                  &
        , l_use_seasalt_direct, l_use_arclsslt                                &
        , l_use_ocff_direct, l_use_arclocff                                   &
        , l_use_biogenic, l_use_arcldlta                                      &
        , l_use_nitrate_direct                                                &
        , l_murk_rad                                                          &
        , sw_control(j)%l_rayleigh, sw_control(j)%l_gas                       &
        , sw_control(j)%l_continuum, sw_control(j)%l_drop                     &
        , sw_control(j)%l_aerosol, sw_control(j)%l_ice                        &
        , sw_control(j)%i_solar_src, sw_spectrum(j) )

    CASE DEFAULT
! DEPENDS ON: r2_sw_specin
      CALL r2_sw_specin                                                       &
        ( ErrorStatus, Cmessage, sw_control(j)%spectral_file                  &
        , sw_control(j)%l_o2, l_climat_aerosol, l_use_dust, l_use_arcldust    &
        , l_use_sulpc_direct, l_use_arclsulp, l_use_soot_direct               &
        , l_use_arclblck, l_use_bmass_direct, l_use_arclbiom                  &
        , l_use_seasalt_direct, l_use_arclsslt, l_use_ocff_direct             &
        , l_use_arclocff, l_use_biogenic, l_use_arcldlta                      &
        , l_use_nitrate_direct, l_murk_rad, sw_control(j)%l_rayleigh          &
        , sw_control(j)%l_gas, sw_control(j)%l_continuum                      &
        , sw_control(j)%l_drop, sw_control(j)%l_aerosol, sw_control(j)%l_ice  &
        , sw_spectrum(j) )
    END SELECT
  END DO

  IF (ErrorStatus  /=  0) THEN

    CALL ereport(routinename, errorstatus, cmessage)
  END IF

  IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                          &
      (lw_control(1)%i_cloud == ip_cloud_mcica)) THEN

    path_end=SCAN(sw_control(1)%spectral_file,'/',.TRUE.)
    mcica_data(:path_end)=sw_control(1)%spectral_file(:path_end)
    mcica_data(path_end+1:)='mcica_data'

    WRITE (6,'(/,A,A)') 'Reading McICA data file: ',mcica_data
    CALL read_mcica_data(mcica_data)
  END IF

!---------------------------------------------------------------------
!     Write restart dump information to tape
!---------------------------------------------------------------------

  IF (tapeout) THEN
    nresdump = INT( ndayin / resdump_days)
    IF (MOD(ndayin, resdump_days)  /=  0) THEN
      nresdump = nresdump + 1
    END IF
    WRITE (55) exname_out,runno_out,nresdump
  END IF

  IF (lhook) CALL dr_hook('RUN_INIT',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE run_init

