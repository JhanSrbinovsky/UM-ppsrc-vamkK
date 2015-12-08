! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! ---------------------------------------------------------------------
! Purpose: Interface level to tr_mix for tracer variables that
!          require boundary layer mixing and/or dry deposition
!          with MOSES II surface scheme.
!
!          Called by BL_TRACER_INTCTL.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
SUBROUTINE bl_trmix_dd(                                                 &
! IN arguments
        bl_levels, dtrdz_charney_grid, tr_vars, dim_cs2,                &
        alpha_cd, rhokh_mix,  p_star,                                   &
! IN Control logicals
        l_murk_advect,l_bl_tracer_mix,l_dust,                           &
        l_sulpc_so2,l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass,        &
        l_ocff, l_nitrate, l_co2_interactive,                           &
        l_co2_emits, l_use_cariolle,                                    &
! IN Emissions fields
        dust_flux, co2_emits, co2flux, npp, resp_s,                     &
! IN variables needed for tr_mix
        kent, we_lim, t_frac, zrzi,                                     &
        kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh, zhsc, z_uv,     &
! IN for dry deposition of tracers
        rho_aresist, aresist,                                           &
        r_b_dust,tstar, land_points, land_index, ice_frac,              &
! IN variables for sresfact
        ntype, ntiles, tile_pts, tile_index, tile_frac,                 &
        canopy, catch, snow_tile, gc,                                   &
        aresist_tile, resist_b_tile, flandg,                            &
! INOUT Fields to mix
        murk, free_tracers, ozone_tracer, drydep_str, resist_b,         &
! INOUT Mineral Dust
        dust_div1,dust_div2,dust_div3,                                  &
        dust_div4,dust_div5,dust_div6,                                  &
! INOUT Sulphur cycle
        so2, dms, so4_aitken, so4_accu, so4_diss, nh3,                  &
! INOUT Soot cycle
        soot_new, soot_aged, soot_cld,                                  &
! INOUT Biomass aerosol
        bmass_new, bmass_agd, bmass_cld,                                &
! INOUT Fossil-fuel organic carbon aerosol
        ocff_new, ocff_agd, ocff_cld,                                   &
! INOUT Ammonium nitrate aerosol
        nitr_acc, nitr_diss,                                            &
! INOUT Carbon cycle
        co2, co2_flux_tot, land_co2,                                    &
! STASH related variables
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
        stashwork,                                                      &
! SCM diagnostics (dummy in full UM)
        nSCMDpkgs, L_SCMDiags                                           &
        )

  USE dust_parameters_mod, ONLY: ndiv, rhop, drep, l_twobin_dust
  USE atm_fields_bounds_mod, ONLY:                                      &
   pdims, tdims, tdims_s, trdims_ltl
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  USE c_bm_bdy_mod, ONLY: resb_freshbmass, resb_agedbmass,              &
                          resb_bmassincloud, ress_bmass
  USE c_nitrate_bdy_mod, ONLY: resb_nitr_acc, resb_nitr_diss,           &
                               ress_nitr_acc, ress_nitr_diss
  USE c_ocff_bdy_mod, ONLY: resb_freshocff, resb_agedocff,              &
                            resb_ocffincloud, ress_ocff
  USE c_st_bdy_mod, ONLY: resb_freshsoot, resb_agedsoot,                &
                          resb_sootincloud, ress_soot
  USE c_sulbdy_mod, ONLY: resb_so2, resb_nh3, resb_so4_ait,             &
                          resb_so4_acc, resb_so4_dis, ress_so2,         &
                          ress_nh3, ress_so4_ait, ress_so4_acc,         &
                          ress_so4_dis, cond_lim, r_snow, asnow
  USE proc_info_mod, ONLY: at_extremity
  USE tr_mix_mod, ONLY: tr_mix
  USE missing_data_mod, ONLY: rmdi
  USE Submodel_Mod

  IMPLICIT NONE

! Model dimensions
  INTEGER, INTENT(IN) ::                                                &
    bl_levels,                                                          &
    tr_vars,                                                            &
                          ! number of free tracer variables
    dim_cs2         ! soil carbon dimension

! Model switches
  LOGICAL, INTENT(IN) ::                                                &
    l_murk_advect,                                                      &
                          ! Switch for advecting aerosol
    l_bl_tracer_mix,                                                    &
                          ! Switch for BL mixing of free tracers
    l_dust,                                                             &
                          ! switch for mineral dust
    l_sulpc_so2,                                                        &
                          ! Switch for Sulphur Cycle
    l_sulpc_nh3,                                                        &
                          ! NH3 included in Sulphur Cycle
    l_sulpc_dms,                                                        &
                          ! DMS included in Sulphur Cycle
    l_soot,                                                             &
                          ! Switch for Soot Cycle
    l_biomass,                                                          &
                          ! Switch for Biomass aerosol
    l_ocff,                                                             &
                          ! Switch for Fossil-fuel OC aerosol
    l_nitrate,                                                          &
                          ! Switch for ammonium nitrate aerosol
    l_co2_interactive,                                                  &
                          ! Switch for interactive CO2
    l_co2_emits,                                                        &
                          ! Switch for anthro CO2 emissions
    l_use_cariolle
                          ! Switch for cariolle ozone tracer scheme
! Model parameters
  REAL, INTENT(IN) ::                                                   &
    alpha_cd(bl_levels),                                                &
    rhokh_mix (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               bl_levels),                                              &
    we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                      !rho*entrainment rate implied by
                                      !placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                      !(z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                      !a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                      !rho*entrainment rate implied by
                                      !  placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                      !(z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                      !a fraction of the timestep
    z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
           bl_levels),                                                  &
                                      !Z_uv(*,K) is height of half
                                      !    level k-1/2.
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                      !Top of decoupled layer
    zh  (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                      !Top of surface mixed layer
    dtrdz_charney_grid(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,bl_levels)
                                      ! For tracer mixing

  INTEGER, INTENT(IN) ::                                                &
    kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                      !grid-level of SML inversion
    kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                      !grid-level of DSC inversion
    land_points,                                                        &
                                      !No.of land points being processed
    land_index(land_points),                                            &
                                      !set from land_sea_mask
! For sresfact
    ntype,                                                              &
                                      !No. of tile types
    ntiles,                                                             &
                                      !No.of land-surface tiles
    tile_pts(ntype),                                                    &
                                      !No.of tile points
    tile_index(land_points,ntype) !Index of tile points

! Required fields (input)
  REAL, INTENT(IN) ::                                                   &
    p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Emissions fields (input)
  REAL, INTENT(IN) ::                                                   &
   dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv)
                                     ! dust emission flux

  REAL, INTENT(IN) ::                                                   &
    co2_emits (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                        ! anthro CO2 emissions
                                        !  (Kg CO2 m-2 s-1)

    co2flux   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                        ! ocean  CO2 emissions
                                        !  (Kg CO2 m-2 s-1)

! For dry deposition of tracers (input)
  REAL, INTENT(IN) ::                                                   &
    rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                       !RHOSTAR*CD_STD*VSHR 
                                       !for CLASSIC aerosol scheme
    aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                       !1/(CD_STD*VSHR)
    r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv), &
                                          !surf layer res for dust
    tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                               ! temperature
    ice_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                       !sea_ice fractn (ice/qrclim.ice)
! For sresfact
    tile_frac(land_points,ntiles),                                      &
                                       !fractional coverage for each
!                                       surface tile
    canopy(land_points,ntiles),                                         &
                                       !Surface/canopy water (kg/m2)
    catch(land_points,ntiles),                                          &
                                       !Surface/canopy water capacity
!                                       of snow-free land tiles (kg/m2)
    snow_tile(land_points,ntiles),                                      &
                                       !Snow on tiles (kg/m2).
    gc(land_points,ntiles),                                             &
                                       !Stomatal conductance to evapn
!                                       for land tiles (m/s)
    aresist_tile(land_points,ntiles),                                   &
                                       ! 1/(CD_STD*VSHR) on land tiles
                                       !for CLASSIC aerosol scheme
    resist_b_tile(land_points,ntiles),                                  &
                                     !(1/CH-1/CD_STD)/VSHR on land tiles
                                       !for CLASSIC aerosol scheme
    flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                     !land fraction of grid box
!                                       (0 if all sea, 1 if all land)
    npp(land_points),                                                   &
                                       ! net primary productivity
                                       !  (Kg C m-2 s-1)
    resp_s(dim_cs2)
                                       ! soil respiration
                                       !  (Kg C m-2 s-1)
! Tracer fields for mixing
  REAL, INTENT(INOUT) ::                                                &
    murk        (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    free_tracers(tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end,                         &
                 trdims_ltl%k_end, tr_vars)

  REAL, INTENT(INOUT) ::                                                &
    drydep_str(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
    resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                       !(1/CH-1/CD_STD)/VSHR
    co2_flux_tot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                       ! total CO2 flux
                                       !  (Kg CO2 m-2 s-1)
    land_co2(land_points)             ! terrestrial CO2 flux
                                       !  (Kg CO2 m-2 s-1)

  REAL, INTENT(INOUT) ::                                                &
    dust_div1   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dust_div2   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dust_div3   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dust_div4   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dust_div5   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dust_div6   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end)
  REAL, INTENT(INOUT) ::                                                &
    so2         (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    dms         (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    so4_aitken  (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    so4_accu    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    so4_diss    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    nh3         (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end)
  REAL, INTENT(INOUT) ::                                                &
    soot_new    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    soot_aged   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    soot_cld    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    bmass_new   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    bmass_agd   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    bmass_cld   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    ocff_new    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    ocff_agd    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    ocff_cld    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    nitr_acc    (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    nitr_diss   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    co2         (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end),         &
    ozone_tracer(tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end, tdims_s%k_end)

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER nSCMDpkgs              ! No of SCM diagnostics packages
  LOGICAL L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages


  REAL, INTENT(INOUT) ::  stashwork(*)

! Argument comdeck declaration
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
!*L------------------COMDECK CCARBON------------------------------------
! Purpose: declares variables and parameters for the carbon cycle
! History:
! version  date         change
! 5.5      26/02/03     add M_CARBON. C Jones.
!----------------------------------------------------------------------
!carbon cycle and vegetation parameters
      REAL                                                              &
     & M_CO2                                                            &
                                  ! molecular weight of CO2
     &,M_AIR                                                            &
                                  ! molecular weight of dry air
     &,M_CARBON                                                         &
                                  ! molecular weight of carbon
     &,EPSILON                                                          &
                                  ! Ratio of molecular weights of water
!                                 !  and dry air.
     &,EPCO2                                                            &
                                  ! Ratio of molecular weights of CO2
!                                 !  and dry air.
     &,EPO2                                                             &
                                  ! Ratio of molecular weights of O2
!                                 !  and dry air.
     &,CO2CONV_A2O                                                      &
                                  ! conversion factor for atmos to
!                                 !  ocean passing of CO2 (mmr to ppmv)
     &,CO2CONV_O2A                ! conversion factor for ocean to
!                                 !  atmos passing of CO2 flux
!                                 !  (mol C/m2/yr to Kg CO2/m2/s)

      PARAMETER (M_AIR=28.966, EPCO2=1.5194, M_CO2=M_AIR*EPCO2,         &
     &           M_CARBON = 12.0, EPSILON = 0.62198, EPO2 = 1.106)

      PARAMETER (CO2CONV_A2O = M_AIR * 1E6 / M_CO2,                     &
     &           CO2CONV_O2A = M_CO2 * 1e-3 / (360.0 * 24.0 * 3600.0))
!*----------------------------------------------------------------------

! Local variables
  INTEGER, PARAMETER :: sect = 3               ! BL section
  INTEGER            :: item                   ! stash item code
  INTEGER            :: im_index               ! internal model
  INTEGER            :: n_tracer               ! looper
  INTEGER            :: errorstatus            ! error reporting
  INTEGER            :: idiv            ! loop counter, dust divs
  INTEGER            :: lev1            ! =1 no. levs for vgrav calc
  CHARACTER (LEN=2)  :: cdiv            ! character string for idiv

  REAL        :: zero_field(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end) ! field of 0s
  REAL        :: tr_flux( pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end,bl_levels)
                                                       ! Fluxes returned
                                                       ! from mixing

  CHARACTER (LEN=*), PARAMETER :: RoutineName='bl_tracer_mixing'
  CHARACTER (LEN=80)           :: cmessage

  INTEGER :: i, j, k, l                !loop counters

  REAL ::                                                               &
    work_2d(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                        !2d work array for stash
    str_resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                        !Rb for S Cycle dry dep
    str_resist_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                        !Rs for S Cycle dry dep
    res_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                        !Ra/(Ra+Rb+Rs) for dry dep
   res_factor_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),&
                                        !Ra/(Ra+Rb+Rs) mean over land
    resist_s(pdims%i_end*pdims%j_end),                                  &
                                        !stomatal resistance
    vstokes1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 !dust settling velocity, lowest layer
    drydep_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                ndiv),                                                  &
                                          ! dust deposition flux
    work1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels,ndiv),                                              &
                                              ! workspace
    work2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
    work3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
    snow_f                          !calculated snow fraction

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
!
!---------------------------------------------------------------------
! Initial setup
!---------------------------------------------------------------------
  IF (lhook) CALL dr_hook('BL_TRMIX_DD',zhook_in,zhook_handle)
  zero_field(:,:) = 0.0
  errorstatus = 0
  im_index = internal_model_index(atmos_im)
!  Initialise RES_FACTOR array to zero to prevent use of unset values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      res_factor(i,j) = 0.0
    END DO
  END DO

!---------------------------------------------------------------------
! Mixing for Aerosol
!---------------------------------------------------------------------
  IF ( l_murk_advect ) THEN
    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rhokh_mix(:,:,1),        &
         dtrdz_charney_grid, zero_field, zero_field,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         murk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              1:bl_levels ),tr_flux,drydep_str                          &
         )


item = 129
    IF ( sf(item,sect) ) THEN      ! diagnostic flux required
! DEPENDS ON: copydiag_3d
      CALL copydiag_3d( stashwork(si(item,sect,im_index)),              &
          tr_flux,                                                      &
          pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,      &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,           &
          stash_levels,num_stash_levels+1,                              &
          atmos_im,sect,item,                                           &
          errorstatus ,cmessage)
    END IF

    IF (errorstatus /=0) THEN
      CALL ereport( RoutineName, errorstatus, cmessage )
    END IF

  END IF

!------------------------------------------------------------------
! Mixing for mineral dust, including dry deposition
!------------------------------------------------------------------

  IF (l_dust) THEN

    lev1 = 1

    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work1(i,j,k,1)=dust_div1(i,j,k)
          work1(i,j,k,2)=dust_div2(i,j,k)
          IF(.NOT.l_twobin_dust) THEN
            work1(i,j,k,3)=dust_div3(i,j,k)
            work1(i,j,k,4)=dust_div4(i,j,k)
            work1(i,j,k,5)=dust_div5(i,j,k)
            work1(i,j,k,6)=dust_div6(i,j,k)
          END IF
        END DO !i
      END DO !j
    END DO !BL_LEVELS

    DO idiv = 1, ndiv

! DEPENDS ON: vgrav
      CALL vgrav(                                                       &
    lev1,drep(idiv),                                                    &
    rhop,p_star,tstar,                                                  &
    vstokes1,work2,work3                                                &
    )
!     Do not allow tracer to settle directly from level 1
!     but treat this in parallel with tracer mixing here. 

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j)=aresist(i,j) * ( vstokes1(i,j) + 1./          &
               (aresist(i,j)+r_b_dust(i,j,idiv)+                        &
               aresist(i,j)*r_b_dust(i,j,idiv)*vstokes1(i,j)) )
        END DO !i
      END DO !j

      CALL tr_mix(                                                      &
! IN fields
            bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,          &
            dtrdz_charney_grid,                                         &
            dust_flux(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end, idiv), res_factor,     &
            kent, we_lim, t_frac, zrzi,                                 &
            kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv, &
! INOUT / OUT fields
            work1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  1:bl_levels,idiv ),tr_flux,drydep_str                 &
            )


      item = 440 + idiv      !dust dry dep flux, from lowest layer
      IF ( sf(item,sect) ) THEN

! Change sign of dry dep flux (otherwise negative)
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            work_2d(i,j) = -drydep_str(i,j)
          END DO
        END DO

! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,        &
         pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                  &
         atmos_im,sect,item,                                            &
         errorstatus,cmessage)

        IF (errorstatus /=0) THEN
          CALL ereport( RoutineName, errorstatus, cmessage )
        END IF

      END IF !stashflag
    END DO !NDIV

    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          dust_div1(i,j,k)=work1(i,j,k,1)
          dust_div2(i,j,k)=work1(i,j,k,2)
          IF(.NOT.l_twobin_dust) THEN 
            dust_div3(i,j,k)=work1(i,j,k,3)
            dust_div4(i,j,k)=work1(i,j,k,4)
            dust_div5(i,j,k)=work1(i,j,k,5)
            dust_div6(i,j,k)=work1(i,j,k,6)
          END IF
        END DO !i
      END DO !j
    END DO !BL_LEVELS

  END IF !L_DUST
!------------------------------------------------------------------
! Mixing for Sulphur Cycle variables, including dry deposition
!------------------------------------------------------------------

  IF (l_sulpc_so2) THEN

!  Calculate resistance factors for species to be dry deposited.

! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.

!   For RESIST_B check to eliminate possible negative values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (resist_b(i,j)  <   0.0)     THEN
          resist_b(i,j) = 0.0
        END IF
      END DO
    END DO

!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.

!   Initialise 'stomatal' resistance array to zero:

    DO k = 1, pdims%i_end*pdims%j_end
      resist_s(k)=0.0
    END DO

! CODE FOR SO2

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_so2*resist_b(i,j)
          str_resist_s(i,j) = ress_so2*resist_s(k)

!  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
!  (0 is acceptable over sea).
          IF (ice_frac(i,j)  >   0.0 )  THEN
            str_resist_s(i,j)=r_snow
          END IF

!  Calculate RES_FACTOR for SO2 (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for SO2 to calculate correct RES_FACTOR_LAND values
! if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                   ntiles, tile_index, tile_pts, .TRUE.,                &
                   canopy, catch, gc, tile_frac, snow_tile,             &
                   aresist, aresist_tile, resb_so2, ress_so2,           &
                   resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for SO2 combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix SO2   (Note emiss added in AERO_CTL via TRSRCE call)

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields 
         so2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             1:bl_levels ), tr_flux, drydep_str                         &
         )

! Write diagnostics to STASH

    item = 270                      !dry dep flux SO2
    IF ( sf(item,sect) ) THEN

! Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF

! CODE FOR NH3 (if present)

    IF (l_sulpc_nh3) THEN

!  Calculate RES_FACTOR for NH3 to allow dry deposition in same way
!  as for SO2 (including code for snow and ice)

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          k=(j-1)*pdims%i_end + i

          IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
            str_resist_b(i,j) = resb_nh3*resist_b(i,j)
            str_resist_s(i,j) = ress_nh3*resist_s(k)

!  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
!  (0 is acceptable over sea).
            IF (ice_frac(i,j)  >   0.0 )  THEN
              str_resist_s(i,j)=r_snow
            END IF

!  Calculate RES_FACTOR for NH3  (only correct for sea points)

            res_factor(i,j) = aresist(i,j) /                            &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

          END IF

        END DO
      END DO

! Call SRESFACT for NH3 to calculate correct RES_FACTOR_LAND values
! if land present:
      IF (land_points >  0) THEN

! DEPENDS ON: sresfact
        CALL sresfact (land_points, land_index,                         &
                   ntiles, tile_index, tile_pts, .TRUE.,                &
                   canopy, catch, gc, tile_frac, snow_tile,             &
                   aresist, aresist_tile, resb_nh3, ress_nh3,           &
                   resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for NH3 combining sea and land values
        DO j = pdims%j_start,pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +       &
                           flandg(i,j)*res_factor_land(i,j)
          END DO
        END DO

      END IF

! Mix NH3   (Note emiss added in AERO_CTL via TRSRCE call)

      CALL tr_mix(                                                      &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT/OUT fields
         nh3(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,      &
             1:bl_levels ), tr_flux, drydep_str                         &
         )

! Write diagnostics to STASH

      item = 300                      !dry dep flux NH3
      IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            work_2d (i,j) = -drydep_str(i,j)
          END DO
        END DO

! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,        &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

        IF (errorstatus /=0) THEN
          CALL ereport( RoutineName, errorstatus, cmessage )
        END IF

      END IF

    END IF                        ! End If nh3 condition

! CODE FOR SO4_AIT

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end+ i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_so4_ait*resist_b(i,j)
          str_resist_s(i,j) = ress_so4_ait*resist_s(k)

!  Calculate RES_FACTOR for SO4_AIT  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for SO4_AIT to calculate correct RES_FACTOR_LAND values
! if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_so4_ait, ress_so4_ait,     &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for SO4_AIT combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix so4_aitken

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         so4_aitken(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
         )

! Write diagnostics to STASH

    item = 271                      !dry dep flux SO4_AIT
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF

! CODE FOR SO4_ACC

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_so4_acc*resist_b(i,j)
          str_resist_s(i,j) = ress_so4_acc*resist_s(k)

!  Calculate RES_FACTOR for SO4_ACC  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
! if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_so4_acc, ress_so4_acc,     &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for SO4_ACC combining sea and land values
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix so4_accu

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         so4_accu(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 272                      !dry dep flux SO4_ACC
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF

! CODE FOR SO4_DIS

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_so4_dis*resist_b(i,j)
          str_resist_s(i,j) = ress_so4_dis*resist_s(k)

!  Calculate RES_FACTOR for SO4_DIS  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
! if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_so4_dis, ress_so4_dis,     &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for SO4_DIS combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix so4_diss

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         so4_diss(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 273                      !dry dep flux SO4_DIS
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF

! CODE FOR DMS (if present)

    IF (l_sulpc_dms) THEN

! (Note no dry dep for DMS, and emiss added in AERO_CTL via TRSRCE call)

! Mix DMS

      CALL tr_mix(                                                      &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),       &
         dtrdz_charney_grid, zero_field, zero_field,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         dms(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             1:bl_levels ), tr_flux, drydep_str                         &
         )

    END IF                          ! End If dms condition

  END IF                            ! END IF l_sulpc_so2

!---------------------------------------------------------------------
! Mixing for Carbon cycle
!---------------------------------------------------------------------
  IF ( l_co2_interactive ) THEN

    DO i = pdims%i_start, pdims%i_end
      DO j = pdims%j_start, pdims%j_end
        co2_flux_tot(i,j)=0.0

!  (i) CO2 emissions from ancillary file.
        IF(l_co2_emits) THEN
          IF ( co2_emits(i,j)  /=  rmdi ) THEN
            co2_flux_tot(i,j)=co2_flux_tot(i,j) + co2_emits(i,j)
          END IF
        END IF

!  (ii) CO2 flux from ocean. (+ve implies air to sea)
        IF ( co2flux(i,j)  /=  rmdi ) THEN
          co2_flux_tot(i,j)=co2_flux_tot(i,j) - co2flux(i,j)
        END IF

      END DO
    END DO

!  (iii) CO2 flux from land processes. (+ve implies biosphere to atmos)
    DO l = 1, land_points
      j = (land_index(l)-1)/pdims%i_end + 1
      i =  land_index(l) - (j-1)*pdims%i_end
      land_co2(l) = (resp_s(l) - npp(l)) * m_co2 / m_carbon
      co2_flux_tot(i,j)=co2_flux_tot(i,j) + land_co2(l)
    END DO

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),       &
         dtrdz_charney_grid, co2_flux_tot, zero_field,                  &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         co2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,      &
             1:bl_levels ), tr_flux, drydep_str                         &
         )

  END IF   ! C-cycle

!------------------------------------------------------------------
! Mixing for Soot variables, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of soot are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of soot.

  IF (l_soot)  THEN   ! If soot modelling is included

!  Calculate resistance factors for species to be dry deposited.

! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.

!   For RESIST_B check to eliminate possible negative values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (resist_b(i,j)  <   0.0)     THEN
          resist_b(i,j) = 0.0
        END IF
      END DO
    END DO

!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.

!   Initialise 'stomatal' resistance array to zero:

    DO k = 1, pdims%i_end*pdims%j_end
      resist_s(k)=0.0
    END DO

! CODE FOR FRESH SOOT

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_freshsoot*resist_b(i,j)
          str_resist_s(i,j) = ress_soot*resist_s(k)

!  Calculate RES_FACTOR for fresh soot  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for fresh soot to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_freshsoot, ress_soot,      &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for fresh soot combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix fresh soot

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         soot_new(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 301                      !dry dep flux fresh soot
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR AGED SOOT

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_agedsoot*resist_b(i,j)
          str_resist_s(i,j) = ress_soot*resist_s(k)

!  Calculate RES_FACTOR for aged soot  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for aged soot to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_agedsoot, ress_soot,       &
                 resist_b_tile,  res_factor_land)

! Recalculate RES_FACTOR for aged soot combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix aged soot

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         soot_aged(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

! Write diagnostics to STASH

    item = 302                      !dry dep flux aged soot
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR CLOUD SOOT

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_sootincloud*resist_b(i,j)
          str_resist_s(i,j) = ress_soot*resist_s(k)

!  Calculate RES_FACTOR for cloud soot  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for cloud soot to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_sootincloud, ress_soot,    &
                 resist_b_tile,  res_factor_land)

! Recalculate RES_FACTOR for cloud soot combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix cloud soot

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         soot_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 303                   !dry (occult) dep flux cloud soot
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

  END IF  ! L_SOOT Soot modelling included

! END of soot section

!------------------------------------------------------------------
! Mixing for Biomass aerosol, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of biomass aerosol are dealt with in aero_ctl.
! Here we consider just boundary layer mixing and dry deposition.

  IF (l_biomass)  THEN   ! If biomass aerosol is included

!  Calculate resistance factors for species to be dry deposited.

! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.

!   For RESIST_B check to eliminate possible negative values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (resist_b(i,j)  <   0.0)     THEN
          resist_b(i,j) = 0.0
        END IF
      END DO
    END DO

!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.

!   Initialise 'stomatal' resistance array to zero:

    DO k = 1, pdims%i_end*pdims%j_end
      resist_s(k)=0.0
    END DO

! CODE FOR FRESH BIOMASS

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_freshbmass*resist_b(i,j)
          str_resist_s(i,j) = ress_bmass*resist_s(k)

!  Calculate RES_FACTOR for fresh biomass  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for fresh biomass to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_freshbmass, ress_bmass,    &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for fresh biomass combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix fresh biomass smoke

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         bmass_new(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,&
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

! Write diagnostics to STASH

    item = 396                      !dry dep flux fresh biomass
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR AGED BIOMASS

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_agedbmass*resist_b(i,j)
          str_resist_s(i,j) = ress_bmass*resist_s(k)

!  Calculate RES_FACTOR for aged biomass  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for aged biomass to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_agedbmass, ress_bmass,     &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for aged smoke combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix aged smoke

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,  zh ,zhsc, z_uv,   &
! INOUT / OUT fields
         bmass_agd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

! Write diagnostics to STASH

    item = 397                      !dry dep flux aged biomass
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR CLOUD BIOMASS

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_bmassincloud*resist_b(i,j)
          str_resist_s(i,j) = ress_bmass*resist_s(k)

!  Calculate RES_FACTOR for cloud smoke  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for cloud smoke to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points >  0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_bmassincloud, ress_bmass,  &
                 resist_b_tile,  res_factor_land)

! Recalculate RES_FACTOR for cloud smoke combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix cloud smoke

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         bmass_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

! Write diagnostics to STASH

    item = 398                   !dry (occult) dep flux cloud smoke
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

  END IF  ! L_BIOMASS Biomass smoke modelling included

! END of biomass section

!------------------------------------------------------------------
! Mixing for Fossil-fuel organic carbon, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of ocff are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of ocff.

  IF (l_ocff)  THEN   ! If ocff modelling is included

!  Calculate resistance factors for species to be dry deposited.

! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.

!   For RESIST_B check to eliminate possible negative values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (resist_b(i,j)  <  0.0)     THEN
          resist_b(i,j) = 0.0
        END IF
      END DO
    END DO

!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.

!   Initialise 'stomatal' resistance array to zero:

    DO k = 1, pdims%i_end*pdims%j_end
      resist_s(k)=0.0
    END DO

! CODE FOR FRESH OCFF

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) < 1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_freshocff*resist_b(i,j)
          str_resist_s(i,j) = ress_ocff*resist_s(k)

!  Calculate RES_FACTOR for fresh ocff  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for fresh ocff to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points > 0) THEN

      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_freshocff, ress_ocff,      &
                 resist_b_tile,  res_factor_land)

! Recalculate RES_FACTOR for fresh ocff combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix fresh ocff

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         ocff_new(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 407                        !dry dep flux fresh ocff
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR AGED OCFF

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) < 1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_agedocff*resist_b(i,j)
          str_resist_s(i,j) = ress_ocff*resist_s(k)

!  Calculate RES_FACTOR for aged ocff  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for aged ocff to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points > 0) THEN

      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_agedocff, ress_ocff,       &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for aged ocff combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix aged OCFF

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         ocff_agd(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 408                  ! dry dep flux aged ocff
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR CLOUD OCFF

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) < 1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_ocffincloud*resist_b(i,j)
          str_resist_s(i,j) = ress_ocff*resist_s(k)

!  Calculate RES_FACTOR for cloud ocff  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for cloud ocff to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points > 0) THEN

      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_ocffincloud, ress_ocff,    &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for cloud ocff combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix cloud ocff

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         ocff_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 409                 !MSR dry (occult) dep flux cloud ocff
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

  END IF  ! L_ocff ocff modelling included

! END of fossil-fuel organic carbon section

!------------------------------------------------------------------
! Mixing for ammonium nitrate, including dry deposition
!------------------------------------------------------------------

! Note:  There are no direct emissions of ammonium nitrate.  Here
! we consider just boundary layer mixing and dry deposition of the
! two modes of ammonium nitrate.

  IF (l_nitrate) THEN  ! If ammonium nitrate modelling is included

!  Calculate resistance factors for species to be dry deposited.

! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.

!   For RESIST_B check to eliminate possible negative values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (resist_b(i,j)  <  0.0)     THEN
          resist_b(i,j) = 0.0
        END IF
      END DO
    END DO

!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.

!   Initialise 'stomatal' resistance array to zero:

    DO k = 1, pdims%i_end*pdims%j_end
      resist_s(k)=0.0
    END DO

! CODE FOR ACCUMULATION-MODE AMMONIUM NITRATE

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) < 1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_nitr_acc * resist_b(i,j)
          str_resist_s(i,j) = ress_nitr_acc * resist_s(k)

!  Calculate RES_FACTOR for accumulation nitrate (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for accumulation nitrate to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points > 0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_nitr_acc, ress_nitr_acc,   &
                 resist_b_tile,  res_factor_land)

! Recalculate RES_FACTOR for accumulation nitrate combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix accumulation nitrate

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         nitr_acc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

! Write diagnostics to STASH

    item = 274             !dry dep flux accumulation nitrate
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN

        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

! CODE FOR DISSOLVED AMMONIUM NITRATE

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) < 1.0) THEN
!  For grid boxes containing sea:
          str_resist_b(i,j) = resb_nitr_diss * resist_b(i,j)
          str_resist_s(i,j) = ress_nitr_diss * resist_s(k)

!  Calculate RES_FACTOR for dissolved nitrate (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
             (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

! Call SRESFACT for dissolved nitrate to calculate correct RES_FACTOR_LAND
! values if land present:
    IF (land_points > 0) THEN

! DEPENDS ON: sresfact
      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, .FALSE.,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_nitr_diss, ress_nitr_diss, &
                 resist_b_tile, res_factor_land)

! Recalculate RES_FACTOR for dissolved nitrate combining sea and land values

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                           flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

! Mix dissolved nitrate

    CALL tr_mix(                                                        &
! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,    &
! INOUT / OUT fields
         nitr_diss(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

! Write diagnostics to STASH

    item = 275             !dry dep flux dissolved nitrate
    IF ( sf(item,sect) ) THEN

!  Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO

! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
          pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,                 &
          atmos_im,sect,item,                                           &
          errorstatus,cmessage)

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END IF  ! SF(item,sect)

  END IF   ! L_nitrate

! END of ammonium nitrate section.

!---------------------------------------------------------------------
! Mixing for Cariolle Ozone tracer
!---------------------------------------------------------------------
  IF ( l_use_cariolle ) THEN
    CALL tr_mix(                                                        &
! IN fields
        bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),        &
        dtrdz_charney_grid, zero_field, zero_field,                     &
        kent, we_lim, t_frac, zrzi,                                     &
        kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,     &
! INOUT / OUT fields
        ozone_tracer(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,1:bl_levels ),           &
        tr_flux, drydep_str                                             &
        )

    item = 480
    IF ( sf(item,sect) ) THEN      ! diagnostic flux required
! DEPENDS ON: copydiag_3d
      CALL copydiag_3d( stashwork(si(item,sect,im_index)),              &
         tr_flux,                                                       &
         pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,       &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         errorstatus ,cmessage)
    END IF

    IF (errorstatus /=0) THEN
      CALL ereport( RoutineName, errorstatus, cmessage )
    END IF

  END IF

!---------------------------------------------------------------------
! Mixing for Free Tracers
! Will mimic 4.5 here, but note that free tracer levels exist only
! at the top of the model so this may be incorrect if
! model_levels /= tr_levels!
!---------------------------------------------------------------------

  IF ( l_bl_tracer_mix ) THEN

    IF (trdims_ltl%k_end /= pdims%k_end) THEN
      WRITE(cmessage,*) 'Tracer mixing gives unexpected results',       &
                        'when model_levels /= tr_levels'
      errorstatus = -10

      CALL ereport( RoutineName, errorstatus, cmessage )
    END IF

    DO n_tracer = 1, tr_vars
      CALL tr_mix(                                                      &
! IN fields
           bl_levels,alpha_cd,rhokh_mix(:,:,2:), rhokh_mix(:,:,1),      &
           dtrdz_charney_grid, zero_field, zero_field,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zh ,zhsc, z_uv,  &
! INOUT / OUT fields
           free_tracers(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,1:bl_levels,n_tracer),&
           tr_flux, drydep_str                                          &
           )


      item = 99 + n_tracer
      IF ( sf(item,sect) )  THEN     ! diagnostic flux required
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d( stashwork(si(item,sect,im_index)),            &
            tr_flux,                                                    &
            pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,    &
            stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
            stash_levels,num_stash_levels+1,                            &
            atmos_im,sect,item,                                         &
            errorstatus, cmessage)
      END IF

      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF

    END DO
  END IF

  IF (lhook) CALL dr_hook('BL_TRMIX_DD',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bl_trmix_dd
