! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Mid-level convection scheme

MODULE mid_conv_6a_mod

IMPLICIT NONE

!
! Description:
!   Mid level convection scheme
!
!   Called by GLUE_CONV.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CONTAINS

SUBROUTINE  mid_conv_6a(nbl,nlev,ntra,n_cca_lev,npnts,trlev,      &
                       bland, w_max, exner_layer_centres,         &
                       exner_layer_boundaries,                    &
                       l_calc_dxek, l_q_interact,                 &
                       l_tracer,midtrig,ntml,ntpar,freeze_lev,    &
                       pstar,p_layer_centres,p_layer_boundaries,  &
                       z_theta, z_rho,                            &
                       r_theta, r_rho,                            &
                       rho_theta, rho,                            &
                       r2rho_th, r2rho,                           & 
                       dr_across_th, dr_across_rh,                &
                       q,th,                                      &
                       timestep,u,v,recip_pstar,qse,              &
                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
                       qcl,tracer,cape_out,cclwp,ccw,cca,         &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                       dqclbydt,dthbydt,                          &
                       dubydt,dvbydt,dtrabydt,                    &
                       detrain_up,detrain_dwn,                    &
                       entrain_up,entrain_dwn,                    &
                       iccb,icct,lcca,                            &
                       lcbase,lctop,rain,snow,rain_3d,snow_3d,    &
                       up_flux,up_flux_half,                      &
                       dwn_flux,tcw,l_mid_all,cca_2d,             &
                       uw_mid,vw_mid, cfl_limited, error_point    &
                       )


USE water_constants_mod, ONLY: lc, lf

USE cv_run_mod, ONLY:                                             &
    l_mom, l_eman_dd, cape_opt, cape_min,                         &
    w_cape_limit, cape_timescale, mid_cmt_opt,                    &
    mid_cnv_pmin, cca2d_md_opt, cca_md_knob, ccw_md_knob,         &
    l_anvil, cnv_wat_load_opt, l_cv_conserve_check, l_3d_cca

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, srf_precip,                            &
    a_land, a_sea, b_land, b_sea, xsbmin,                         &
    thpixs_mid, qpixs_mid, thpixs_deep, qpixs_deep,               &
    mparb, c_mid, d_mid, wcape_fac, delthst

USE cv_dependent_switch_mod, ONLY:                                &
    md_on, mdet_md_on, md_ent_on, md_new_termc

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                     &
    flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,         &
    flg_mf_midlev

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: cp, c_virtual

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE PrintStatus_mod
USE lift_par_6a_mod
USE convec2_6a_mod
USE water_loading_mod
USE cor_engy_6a_mod

IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:

INTEGER, INTENT(IN) :: &
  nbl                  & ! No. of boundary layer levels
 ,nlev                 & ! No. of model layers
 ,ntra                 & ! No. of tracer fields
 ,n_cca_lev            & ! No. of convective cloud amount levels (1 for 2D,
                         ! nlevs for 3D)
 ,npnts                & ! No. of deep convection points
 ,trlev                  ! No. of model levels on which tracers are included

LOGICAL, INTENT(IN) :: bland(npnts) ! Land/sea mask

REAL, INTENT(IN)    ::                 &
  exner_layer_centres(npnts,0:nlev)    & ! Exner
 ,exner_layer_boundaries(npnts,0:nlev)   ! Exner at half level above
                                         ! exner_layer_centres

LOGICAL, INTENT(IN) ::  &
  l_calc_dxek           & ! Switch for calculation of condensate increment
 ,l_q_interact          & ! Switch allows overwriting parcel variables when
                          ! calculating condensate incr.
 ,l_tracer                ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: &
  midtrig(npnts)       & ! Lowest trigger levelfor convection
 ,ntml(npnts)          & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(npnts)         & ! Top level of initial parcel ascent in BL scheme 
                         ! defined relative to theta,q grid
 ,freeze_lev(npnts)      ! freezing level

REAL, INTENT(IN)    ::             &
  pstar(npnts)                     & ! Surface pressure (Pa)
 ,p_layer_centres(npnts,0:nlev)    & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:nlev)   ! Pressure at half level above
                                     ! p_layer_centres (Pa)

REAL, INTENT(IN) :: z_theta(npnts,nlev)      ! height of theta levels (m)
REAL, INTENT(IN) :: z_rho(npnts,nlev)        ! height of rho levels (m)
REAL, INTENT(IN) :: r_theta(npnts,0:nlev)    ! radius of theta levels (m)
REAL, INTENT(IN) :: r_rho(npnts,nlev)        ! radius of rho levels (m)
REAL, INTENT(IN) :: rho_theta(npnts,nlev)    ! density for theta lev (kg/m3)
REAL, INTENT(IN) :: rho(npnts,nlev)          ! density for rho lev (kg/m3)
REAL, INTENT(IN) :: r2rho_th(npnts,nlev)     ! radius**2 density for 
                                             ! theta lev (kg/m)
REAL, INTENT(IN) :: r2rho(npnts,nlev)        ! radius**2 density for 
                                             ! rho lev (kg/m)
REAL, INTENT(IN) :: dr_across_th(npnts,nlev) ! thickness of theta levels (m)
REAL, INTENT(IN) :: dr_across_rh(npnts,nlev) ! thickness of rho levels (m)

REAL, INTENT(IN)    ::  &
  q(npnts,nlev)         & ! Model mixing ratio (kg/kg)
 ,th(npnts,nlev)        & ! Model potential temperature (K)
 ,timestep              & ! Model timestep (s)
 ,u(npnts,nlev)         & ! Model U field (m/s)
 ,v(npnts,nlev)         & ! Model V field (m/s)
 ,w_max(npnts)          & ! max w in column
                          ! for use in scale dependent cape timescale
 ,recip_pstar(npnts)    & ! Reciprocal of pstar array
 ,qse(npnts,nlev)         ! Saturation mixing ratio of
                          ! cloud environment (kg/kg)

! Arguments with intent INOUT:

REAL, INTENT(INOUT) ::   &
  bulk_cf(npnts,nlev)    & ! Bulk total cloud volume ( )
 ,cf_frozen(npnts,nlev)  & ! Frozen water cloud volume ( )
 ,cf_liquid(npnts,nlev)  & ! Liq water cloud volume ( )
 ,qcf(npnts,nlev)        & ! Ice condensate mix ratio (kg/kg)
 ,qcl(npnts,nlev)          ! Liq condensate mix ratio (kg/kg)

REAL, INTENT(INOUT) :: tracer(npnts,trlev,ntra) !Model tracer fields (kg/kg)

REAL, INTENT(INOUT) :: tcw(npnts) ! Total condensed water(kg/m2/s)

! Arguments with intent OUT:

REAL, INTENT(OUT) ::   &
  cape_out(npnts)      & ! Saved convective available
                         ! potential energy for diagnostic output (J/kg)
 ,cclwp(npnts)         & ! Condensed water path (kg/m^2)
 ,ccw(npnts,nlev)      & ! Convective cloud liquid water
                         ! on model levels (kg/kg)
 ,cca(npnts,n_cca_lev) & ! Convective cloud amount on model levels (0-1)
 ,dbcfbydt(npnts,nlev) & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(npnts,nlev) & ! Increments to ice cloud volume due to convection(/s)
 ,dcflbydt(npnts,nlev) & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(npnts,nlev)   & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(npnts,nlev) & ! Increments to ice condensate due to convection
                         ! (kg/kg/s)
 ,dqclbydt(npnts,nlev) & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(npnts,nlev)  & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(npnts,nlev+1) & ! Increments to U due to CMT (m/s2)
 ,dvbydt(npnts,nlev+1)   ! Increments to V due to CMT (m/s2)

REAL, INTENT(OUT) ::        &
  dtrabydt(npnts,nlev,ntra)   !Increment to tracer convection (kg/kg/s)

REAL, INTENT(OUT) ::       &
  detrain_up(npnts,nlev)   & ! Fractional detrainment rate into updraughts
                             ! (Pa/s)
 ,detrain_dwn(npnts,nlev)  & ! Fractional detrainment rate into downdraughts 
                             ! (Pa/s)
 ,entrain_up(npnts,nlev)   & ! Fractional entrainment rate into updraughts
                             ! (Pa/s)
 ,entrain_dwn(npnts,nlev)    ! Fractional entrainment rate into downdraughts
                             ! (Pa/s)

INTEGER, INTENT(OUT) :: &
  iccb(npnts)           & ! Convective cloud base level
 ,icct(npnts)             ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(npnts) ! Lowest conv. cloud amt. (%)

INTEGER, INTENT(OUT) :: &
  lcbase(npnts)         & ! Lowest conv. cloud base level
 ,lctop(npnts)            ! Lowest conv. cloud top level

REAL, INTENT(OUT) ::  &
  rain(npnts)         & ! Surface convective rainfall (kg/m2/s)
 ,snow(npnts)         & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(npnts,nlev) & ! Convective rainfall flux (kg/m2/s)
 ,snow_3d(npnts,nlev) & ! Convective snowfall flux (kg/m2/s)
 ,up_flux(npnts,nlev) & ! Updraught mass flux (Pa/s)
 ,up_flux_half(npnts,nlev) & ! Updraught mass flux on half levels(Pa/s)
 ,dwn_flux(npnts,nlev)  ! Downdraught mass flux (Pa/s)

LOGICAL, INTENT(OUT) :: &
  l_mid_all(npnts)        ! Points where mid level convection
                          ! occurs at some level
REAL,INTENT(INOUT) :: &
  cca_2d(npnts)         !2D convective cloud amount(%)

! CMT diagnostics
REAL, INTENT(OUT) ::      &
  uw_mid(npnts,nlev)      & ! U component of stress from mid-level convection
                            ! (kg/m/s2)
 ,vw_mid(npnts,nlev)        ! V component of stress from mid-level convection
                            ! (kg/m/s2)
REAL, INTENT(OUT) ::      &
  cfl_limited(npnts)        ! Inidicator of CFL limited mid-level convection

INTEGER, INTENT(OUT) ::   &
  error_point               ! 0 no problem, > 0 problem point
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

REAL :: zk(npnts)               !heights for use in calc

REAL :: zkp12(npnts)            !of moist static energy

REAL :: zkp1(npnts)

REAL :: entrain_coef(npnts)     ! entrainment coefficients unset

INTEGER :: index1(npnts),index2(npnts)

INTEGER :: ncposs               ! No. of points which may convect

INTEGER :: nconv                ! No. of convecting points

REAL :: amdetk(npnts)           ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

REAL :: cape(npnts)             ! Convective available potential
                                ! energy (J/kg)

REAL :: dcpbydt(npnts)          ! Rate of change of cape (J/kg/s)

REAL :: depth(npnts)            ! Depth of convective cloud (m)

REAL :: delexkp1(npnts)         ! Difference in exner ratio
                                ! across layer k+1

REAL :: eminds(npnts)           ! Minimum buoyancy for convection
                                ! to initiate from level k
                                ! (Kelvin)

REAL :: ekp14(npnts)            ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(Dimensionless)

REAL :: ekp34(npnts)            ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekm14(npnts)            ! Entrainment coefficients at
                                ! level k-1+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: exk(npnts)              ! Exner ratio at layer k

REAL :: exkp1(npnts)            ! Exner ratio at layer k+1

REAL :: flxmax(npnts)           ! Maximum initial convective
                                ! mass flux (Pa/s)

REAL :: flx_init(npnts)         ! Initial mass flux at cloud base (Pa/s)

REAL :: flx_init_new(npnts)     ! flx_init scaled to destroy cape
                                ! over timescale cape_timescale (Pa/s)

REAL :: flxmax_init(npnts)      ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer in Pa/s)

REAL :: max_cfl(npnts)          ! Max cfl ratio over a convecting
                                ! layer

REAL :: precip(npnts,nlev)      ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

REAL :: pk(npnts)               ! Pressure at midpoint of layer k (Pa)

REAL :: pkp1(npnts)             ! Pressure at midpoint of layer k+1 (Pa)

REAL :: delpk(npnts)            ! Pressure difference over layer k (Pa)

REAL :: delpkp1(npnts)          ! Pressure difference over layer k+1 (Pa)

REAL :: delpkp12(npnts)         ! Pressure difference between
                                ! layers k and k+1 (Pa)

REAL :: delp_uv_k(npnts)        ! Pressure difference across uv layer k (Pa)

REAL :: delp_uv_kp1(npnts)      ! Pressure difference across uv layer k+1 (Pa)

REAL :: rhum(npnts)             ! Dummy relative humidity
                                ! (only used on shallow points)

LOGICAL :: bgmk(npnts)          ! Mask for points where parcel in
                                ! layer k is saturated
LOGICAL :: blatent(npnts)       ! Mask for points where latent heat has 
                                ! been released

LOGICAL :: bwater(npnts,2:nlev) ! Mask for points at which
                                ! condensate is liquid

LOGICAL :: bwk(npnts)           ! Mask for liquid condensate on k
LOGICAL :: bwkp1(npnts)         ! Mask for liquid condensate on k+1

LOGICAL :: bterm(npnts)         ! Mask for points which have
                                ! stopped convecting

LOGICAL :: bconv(npnts)         ! Mask for points at which
                                ! convection is occurring

LOGICAL :: bcposs(npnts)        ! Mask for points passing
                                ! initial stability test

LOGICAL :: binit(npnts)         ! Mask for points initiating

! Parcel variables

REAL :: qpi(npnts)              ! Initial parcel mixing ratio (kg/kg)

REAL :: qp(npnts,nlev)          ! Parcel mixing ratio (kg/kg)

REAL :: thpi(npnts)             ! Initial parcel potential temp.(K)

REAL :: thp(npnts,nlev)         ! Parcel potential temp (K)

REAL :: up(npnts,nlev)          ! Parcel U (m/s)

REAL :: vp(npnts,nlev)          ! Parcel V (m/s)

REAL :: trap(npnts,nlev,ntra)   ! Tracer content of parcel (kg/kg)

REAL :: expi(npnts)             ! Initial parcel exner pressure

REAL :: flx(npnts,nlev)         ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(npnts,nlev)    ! Minmum parcel buoyancy excess

REAL :: thpixs_v(npnts,nlev)    ! Theta parcel excess (K)

REAL :: qpixs_v(npnts,nlev)     ! Q parcel excess(kg/kg)

REAL :: dpmin(npnts)            ! work array for parcel excess cal

REAL :: qclp(npnts,nlev)         ! Parcel liquid condensate mixing
                                 ! ratio in layer k (kg/kg)

REAL :: qcfp(npnts,nlev)         ! Parcel frozen condensate mixing
                                 ! ratio in layer k (kg/kg)

! Parameters

REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed
REAL, PARAMETER :: minflx = TINY(flx_init_new)  ! minimum allowable
                                                ! initial mass flux

! CMT variables

INTEGER :: kterm(npnts)         ! Level index for termination of

REAL :: eflux_u_ud(npnts)       ! Vertical eddy flux of momentum
                                ! due to UD at top of layer (Pa m/s2)

REAL :: eflux_v_ud(npnts)       ! Vertical eddy flux of momentum
                                ! due to UD at bottom of layer (Pa m/s2)

LOGICAL :: l_mom_gk             ! true if Gregory-Kershaw CMT

! Cape scaling/closure variables

INTEGER :: start_lev(npnts)     !Level at which convection
                                !initiates

INTEGER :: start_lev_c2(npnts)   ! PC2 Compressed convection
                                ! initiation level

INTEGER :: det_lev(npnts)       ! Level at which split final
                                ! detrainment last occurred

INTEGER :: nterm                ! No. of points where conv.
                                ! has terminated

INTEGER :: index_nterm(npnts)   ! Index for points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage

REAL :: scale_f(npnts)          ! scaling factor

REAL :: cape_ts_new(npnts)      ! Used as variable in RH-based
                                ! closure

REAL :: relh(npnts)             ! RH integral (average when
                                ! convection terminates)

REAL :: dptot(npnts)            ! Delta P integral

! Downdraught scheme variables

INTEGER :: npossdd              ! Max. no. of downdraughts
                                ! possible

INTEGER :: nnodd                ! No. of downdraughts not possible

INTEGER :: index_possdd(npnts)  ! Index of downdraughts possible

INTEGER :: index_nodd(npnts)    ! Index of downdraughts not
                                ! possible

REAL :: deltap_cld              ! Pressure thickness of convective
                                ! cloud (Pa)

! Arrays required by Emanuel downdraughts

INTEGER ::          &
  kterm_mid(npnts)  & ! termination level of highest mid level
, kterm_max


! Local compressed arrays

LOGICAL ::          &
  bgmkp1_c(npnts)   & ! Mask for points where parcel in layer k+1 is saturated
 ,bgmkp1_c2(npnts)  & ! Mask for points where parcel in layer k+1 is saturated
 ,bwk_c(npnts)      & ! bwater mask in layer k
 ,bwk_c2(npnts)     & ! bwater mask in layer k
 ,bwkp1_c(npnts)    & ! bwater mask in layer k+1
 ,bwkp1_c2(npnts)     ! bwater mask in layer k+1

LOGICAL :: blatent_c2(npnts)    ! Mask for points where latent heat has 
                                ! been released

REAL :: deltak_c2(npnts)        ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

REAL :: dqek_c2(npnts)          ! Increment to q due to
                                ! convection in layer k (kg/kg)

REAL :: dqekp1_c2(npnts)        ! Increment to q due to
                                ! convection in layer k+1 (kg/kg)

REAL :: dthek_c2(npnts)         ! Increment to potential temp.
                                ! due to convection in layer k

REAL :: dthekp1_c2(npnts)       ! Increment to potential temp.
                                ! due to convection in layer k+1

REAL :: dtraek_c2(npnts,ntra)   ! Increment to model tracer due
                                ! to conv. at level k (kg/kg/s)

REAL :: dtraekp1_c2(npnts,ntra) ! Increment to model tracer due
                                ! to conv. at level k+1 (kg/kg/s)

REAL :: duek_c2(npnts)          ! Increment to model U in layer k
                                ! due to CMT (m/s2)

REAL :: duekp1_c2(npnts)        ! Increment to model U in layer
                                ! k+1 due to CMT (m/s2)

REAL :: dvek_c2(npnts)          ! Increment to model V in layer k
                                ! due to CMT (m/s2)

REAL :: dvekp1_c2(npnts)        ! Increment to model V in layer
                                ! k+1 due to CMT (m/s2)

REAL :: flxk_c(npnts), flxk_c2(npnts) !Parcel mass flux in layer k
                                ! (Pa/s)

REAL :: flxkp12_c2(npnts)       ! Half level mass flux (Pa/s)

REAL :: flxkp1_c2(npnts)         ! Parcel mass flux in layer k+1
                                ! (Pa/s)

REAL :: prekp1_c2(npnts)        ! Precip. from parcel as it rises
                                ! from layer k to k+1 (kg/m2/s)

REAL ::           &
  qpk_c(npnts)    & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpk_c2(npnts)   & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpk(npnts)      & ! ad. entrain.
 ,qpkp1_c(npnts)  & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qpkp1_c2(npnts) & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qek_c(npnts)    & ! Env. mixing ratio in layer k (kg/kg)
 ,qek_c2(npnts)   & ! Env. mixing ratio in layer k (kg/kg)
 ,qek(npnts)      & ! for ad entrain.
 ,qekp1_c(npnts)  & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1_c2(npnts) & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1(npnts)    & ! for ad entrain.
 ,qsek_c2(npnts)  & ! Saturation mixing ratio of cld. env. in layer k (kg/kg)
 ,qsek(npnts)     & ! for ad entrain.
 ,qsekp1_c2(npnts)& ! Saturation mixing ratio of cld. env. in layer k+1 (kg/kg)
 ,qsekp1(npnts)     !for ad entrain.

REAL ::            &
  thek_c(npnts)    & ! Env. potential temp in layer k (K)
 ,thek_c2(npnts)   & ! Env. potential temp in layer k (K)
 ,thek(npnts)      & ! for ad entrain.
 ,thekp1_c(npnts)  & ! Env. potential temp in layer k (K)
 ,thekp1_c2(npnts) & ! Env. potential temp in layer k (K)
 ,thekp1(npnts)    & ! for ad entrain.
 ,thpk_c(npnts)    & ! Parcel potential temp in layer k (K)
 ,thpk_c2(npnts)   & ! Parcel potential temp in layer k (K)
 ,thpk(npnts)      & ! for ad entrain.
 ,thpkp1_c(npnts)  & ! Parcel potential temp in layer k (K)
 ,thpkp1_c2(npnts)   ! Parcel potential temp in layer k (K)

REAL :: traek_c2(npnts,ntra)  ! Env. Tracer content in layer k (kg/kg)
REAL :: traekp1_c2(npnts,ntra)! Env. Tracer content in layer k+1 (kg/kg)
REAL :: trapk_c2(npnts,ntra)  ! Parcel Tracer content in layer k (kg/kg)
REAL :: trapkp1_c2(npnts,ntra)! Parcel Tracer content in layer k+1 (kg/kg)

REAL :: rbuoyk_c(npnts), rbuoyk_c2(npnts)       ! Par. buoyancy at k (K)
REAL :: rbuoykp1_c(npnts),rbuoykp1_c2(npnts)    ! Par. buoyancy at k+1 (K)

REAL :: watldek_c(npnts), watldek_c2(npnts)     ! Env. water loading
                                                ! in layer k (kg/kg)
REAL :: watldpk_c(npnts), watldpk_c2(npnts)     ! Par. water loading
                                                ! in layer k (kg/kg)
REAL :: watldekp1_c(npnts), watldekp1_c2(npnts) ! Env. water loading
                                                ! in layer k+1 (kg/kg)
REAL :: watldpkp1_c(npnts), watldpkp1_c2(npnts) ! Par. water loading
                                                ! layer k+1 (kg/kg)

REAL :: Qlkp1_c(npnts),   &   ! Amount of condensation to liquid water 
        Qlkp1_c2(npnts)       ! in the parcel (kg/kg)
REAL :: Qfkp1_c(npnts),   &   ! Amount of deposition to ice water
        Qfkp1_c2(npnts)       ! in the parcel (kg/kg)
REAL :: Frezkp1_c(npnts), &   ! Amount of freezing from liquid
        Frezkp1_c2(npnts)     ! to frozen water in the parcel (kg/kg)

REAL :: uek_c2(npnts)    ! Env. U in layer k (m/s)
REAL :: uekp1_c2(npnts)  ! Env. U in layer k+1 (m/s)
REAL :: vek_c2(npnts)    ! Env. V in layer k (m/s)
REAL :: vekp1_c2(npnts)  ! Env. V in layer k+1 (m/s)
REAL :: upk_c2(npnts)    ! Parcel U in layer k (m/s)
REAL :: upkp1_c2(npnts)  ! Parcel U in layer k+1 (m/s)
REAL :: vpk_c2(npnts)    ! Parcel V in layer k (m/s)
REAL :: vpkp1_c2(npnts)  ! Parcel V in layer k+1 (m/s)

REAL ::              &
  qclek_c(npnts)     & ! Environment liquid condensate mixing ratio
 ,qclek_c2(npnts)    & ! in layer k (kg/kg)
 ,qclekp1_c(npnts)   & ! Environment liquid condensate mixing ratio
 ,qclekp1_c2(npnts)  & ! in layer k+1 (kg/kg)
 ,qcfek_c(npnts)     & ! Environment frozen condensate mixing ratio in
 ,qcfek_c2(npnts)    & ! layer k (kg/kg)
 ,qcfekp1_c(npnts)   & ! Environment frozen condensate mixing ratio in
 ,qcfekp1_c2(npnts)  & ! layer k+1 (kg/kg)
 ,qclpk_c(npnts)     & ! Parcel liquid condensate mixing ratio in
 ,qclpk_c2(npnts)    & ! layer k (kg/kg)
 ,qclpkp1_c(npnts)   & ! Parcel liquid condensate mixing ratio in
 ,qclpkp1_c2(npnts)  & ! layer k+1 (kg/kg)
 ,qcfpk_c(npnts)     & ! Parcel frozen condensate mixing ratio in
 ,qcfpk_c2(npnts)    & ! layer k (kg/kg) 
 ,qcfpkp1_c(npnts)   & ! Parcel frozen condensate mixing ratio in
 ,qcfpkp1_c2(npnts)    !  layer k+1 (kg/kg)

REAL :: cflek_c2(npnts),cflekp1_c2(npnts)
                                ! Environment liquid water
                                ! cloud volume ( )

REAL :: cffek_c2(npnts),cffekp1_c2(npnts)
                                ! Environment frozen water
                                ! cloud volume ( )

REAL :: bcfek_c2(npnts),bcfekp1_c2(npnts)
                                ! Environment bulk total
                                ! cloud volume ( )

REAL :: dqclek_c2(npnts),dqclekp1_c2(npnts)
                                ! Environment increments
                                ! to liquid condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dqcfek_c2(npnts),dqcfekp1_c2(npnts)
                                ! Environment increments
                                ! to frozen condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dcflek_c2(npnts),dcflekp1_c2(npnts)
                                ! Environment increments
                                ! to liquid water cloud volume due
                                ! to convection (/s)

REAL :: dcffek_c2(npnts),dcffekp1_c2(npnts)
                                ! Environment increments
                                ! to frozen water cloud volume due
                                ! to convection (/s)

REAL :: dbcfek_c2(npnts),dbcfekp1_c2(npnts)
                                ! Environment increments
                                ! to bulk total cloud volume due
                                ! to convection (/s)

REAL :: amdetk_c2(npnts)
LOGICAL :: bgmk_c2(npnts)
LOGICAL :: bland_c2(npnts)
LOGICAL :: blowst_c2(npnts)
LOGICAL :: bterm_c2(npnts)
REAL :: cape_c2(npnts)
REAL :: cca_2d_c2(npnts)
REAL :: cclwp_c2(npnts)
REAL :: ccwkp1_c2(npnts)
REAL :: dcpbydt_c2(npnts)
REAL :: delpk_c2(npnts)
REAL :: delpkp12_c2(npnts)
REAL :: delpkp1_c2(npnts)
REAL :: delp_uv_k_c2(npnts)
REAL :: delp_uv_kp1_c2(npnts)
REAL :: depth_c2(npnts)
REAL :: dptot_c2(npnts)
REAL :: eflux_u_ud_c2(npnts)
REAL :: eflux_v_ud_c2(npnts)
REAL :: ekp14_c(npnts),ekp14_c2(npnts)
REAL :: ekp34_c(npnts),ekp34_c2(npnts)
REAL :: eminds_c(npnts)
REAL :: exk_c(npnts), exk_c2(npnts)
REAL :: exkp1_c(npnts),exkp1_c2(npnts)
REAL :: expi_c2(npnts)
INTEGER :: icct_c2(npnts)
INTEGER :: iccb_c2(npnts)
INTEGER :: lctop_c2(npnts)
INTEGER :: lcbase_c2(npnts)
REAL :: lcca_c2(npnts)
REAL :: max_cfl_c2(npnts)
REAL :: pk_c(npnts),pk_c2(npnts)
REAL :: pkp1_c(npnts),pkp1_c2(npnts)
REAL :: pstar_c2(npnts)
REAL :: qpi_c2(npnts)
REAL :: relh_c2(npnts)
REAL :: rbuoy_p_here_c2(npnts)
REAL :: the_here_c2(npnts)
REAL :: thp_here_c2(npnts)
REAL :: qe_here_c2(npnts)
REAL :: qp_here_c2(npnts)
REAL :: tcw_c2(npnts)
REAL :: thpi_c2(npnts)
REAL :: xsbmin_v_c2(npnts)
REAL :: qsat_lcl(npnts)   ! Not used
REAL ::        &
  ekp14_plus1  & ! (1+epk14)
 ,repss          ! 1/((1+epk14)*(1+epk34)) 


!===============================================================
! CCRad Variables local variables
!===============================================================
INTEGER          :: dum_iccb    ! Holding variables for
INTEGER          :: dum_icct    ! identifying cld base and tops
                                ! on gridpoints with more than 1
                                ! (possibly 2) mid-level
                                ! convective events
                                !
INTEGER          :: n_mdcld     ! Number of gridpoints with
                                ! single mid-level cloud banks
                                !
INTEGER          :: n_mdcld_mult   ! Number of gridpoints with
                                   ! multiple mid-level cloud
                                   ! banks
                                   !

!===============================================================
! Allocatable arrays, because we do not know how many gridpoints
! will contain mid-level until we test for them.
!===============================================================

INTEGER                            :: dum1(npnts)
INTEGER                            :: dum2(npnts)

INTEGER, ALLOCATABLE :: mdcldi(:)
                                ! INDICES in full array of
                                ! gridpoints single mid-level
                                ! cloud banks
INTEGER, ALLOCATABLE :: mdcldi_mult(:)
                                ! index in full array of
                                ! gridpoints multiple mid-level
                                ! cloud banks

!===============================================================
! Compressed arrays for gridpoints which have one bank of
! mid-level cloud and multiple banks of mid-level cloud.
! Requires the use of compressed allocatable arrays because the
! location and number of gridpoints with mid-level has not been
! diagnosed yet.
!===============================================================

!===============================================================
! For multiple mid-level cloud
!===============================================================
INTEGER, ALLOCATABLE ::  iccb_md_c(:)
INTEGER, ALLOCATABLE ::  icct_md_c(:)
INTEGER, ALLOCATABLE ::  freeze_lev_md_c(:)
REAL, ALLOCATABLE ::  cca_2d_md_c(:)
REAL, ALLOCATABLE ::  cca_md_c(:,:)
REAL, ALLOCATABLE ::  ccw_md_c(:,:)
REAL, ALLOCATABLE ::  z_theta_md_c(:,:)
REAL, ALLOCATABLE ::  z_rho_md_c(:,:)
REAL, ALLOCATABLE ::  p_lyr_bnds_md_c(:,:)

!===============================================================
! For multiple mid-level cloud
!===============================================================
INTEGER, ALLOCATABLE ::  iccb_md_mult_c(:)
INTEGER, ALLOCATABLE ::  icct_md_mult_c(:)
INTEGER, ALLOCATABLE ::  freeze_lev_md_mult_c(:)
REAL, ALLOCATABLE ::  cca_2d_md_mult_c(:)
REAL, ALLOCATABLE ::  cca_md_mult_c(:,:)
REAL, ALLOCATABLE ::  ccw_md_mult_c(:,:)
REAL, ALLOCATABLE ::  z_theta_md_mult_c(:,:)
REAL, ALLOCATABLE ::  z_rho_md_mult_c(:,:)
REAL, ALLOCATABLE ::  p_lyr_bnds_md_mult_c(:,:)

!===============================================================
! End CCRad Variables local variables
!===============================================================

!   required by check on -ve q

REAL :: qminincolumn(npnts)     ! Minimum value for q in column
                                ! (kg/kg)
REAL :: temp1(npnts)            ! work array
REAL :: temp2(npnts)            ! work array
REAL :: temp3(npnts)            ! work array

REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

INTEGER ::        &
 nmid             & ! total number of points where mid-level
                    ! convection occurs
,nmax_layers        ! Maximum number of allow mid-level layers

REAL ::      &
  rh_test    &   ! critical RH value for CAPE_OPT=6
 ,rh_fac         ! factor for CAPE_OPT=6 calculation

! Loop counters

INTEGER :: i,j,k,ktra,kt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



IF (lhook) CALL dr_hook('MID_CONV_6A',zhook_in,zhook_handle)
error_point = 0

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

DO i = 1,npnts
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  l_mid_all(i) = .FALSE.
  blatent(i)   = .FALSE.
END DO

!-----------------------------------------------------------------------
! Decide whether Gregory-Kershaw CMT scheme is to be used if l_mom is true.
! The Gregory-Kershaw scheme requires calculations in the main plume ascent
! loop whereas the alternative diffusive scheme is called after the plume
! calculation.

IF (mid_cmt_opt == 1 ) THEN
  l_mom_gk = .FALSE.       ! Use diffusive scheme
ELSE
  l_mom_gk = l_mom         ! Use Gregory-Kershaw scheme
END IF

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!-----------------------------------------------------------------------
DO i = 1,npnts
  kterm(i)       = 0
  kterm_mid(i)   = 0
END DO


!intialise parcel values over all levels
DO k = 1, nlev
  DO i = 1, npnts
    qp(i,k)     = 0.0
    thp(i,k)    = 0.0
    qclp(i,k)   = 0.0
    qcfp(i,k)   = 0.0
    flx(i,k)    = 0.0
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
  END DO
END DO

IF (l_mom_gk) THEN
  DO k=1,nlev
    DO i = 1,npnts
      up(i,k) = 0.0
      vp(i,k) = 0.0
    END DO
  END DO
END IF

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k=1,nlev
      DO i = 1,npnts
        trap(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF

DO k=1,nlev
  DO i=1, npnts
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
  END DO
END DO

IF (l_mom) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
  DO k = 1,nlev
    DO i = 1,npnts
      uw_mid(i,k) = 0.0
      vw_mid(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,npnts
        dtrabydt(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
IF (flg_up_flx .OR. flg_mf_midlev) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux(i,k)      = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux_half(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dwn_flux(i,k)     = 0.0
    END DO
  END DO
END IF
! Now required in all cases
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_up(i,k)   = 0.0
    END DO
  END DO
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_up(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_dwn(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_dwn(i,k)  = 0.0
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
! Zeroes values to remove mid-level dependence on shallow and deep
! cloud
DO i = 1,npnts
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
END DO

DO k= 1,n_cca_lev
  DO i= 1,npnts
    cca(i,k) = 0.0
  END DO
END DO


DO i = 1,npnts
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling and closure calculations
!-----------------------------------------------------------------------
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)         = 0.0
  cape_out(i)     = 0.0
  dcpbydt(i)      = 0.0
  max_cfl(i)      = 0.0
  det_lev(i)      = 0
  start_lev(i)    = 0
  relh(i)         = 0.0
  dptot(i)        = 0.0
  cfl_limited(i)  = 0.0

!-----------------------------------------------------------------------
! 2.6  Initialise eddy flux arrays for updraught
!-----------------------------------------------------------------------
  eflux_u_ud(i)   = 0.0
  eflux_v_ud(i)   = 0.0

!-----------------------------------------------------------------------
! 2.7  Initialise surface precipitation arrays
!-----------------------------------------------------------------------
  rain(i)         = 0.0
  snow(i)         = 0.0
END DO

DO i = 1,npnts
  entrain_coef(i) = -99.0      ! not set
END DO

!initialise ekm14 for unused adaptive entrainment
DO i =1, npnts
  ekm14(i)    = 0.0
END DO

!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate xsbmin_v, thpixs_v and qpixs_v constants based on layer
! thickness (Pa).
DO k = 1,nlev-1
  DO i = 1,npnts
    dpmin(i) = MIN( ((p_layer_centres(i,k) -                    &
                         p_layer_centres(i,k+1))/5000.0),1.0)

    xsbmin_v(i,k) = dpmin(i) *0.2

    IF (midtrig(i)  ==  ntml(i)                                 &
                     .OR. midtrig(i)  ==  ntml(i)+1) THEN

      thpixs_v(i,k) = dpmin(i) * thpixs_deep
      qpixs_v(i,k)  = qpixs_deep
    ELSE
      thpixs_v(i,k) = dpmin(i) * thpixs_mid
      qpixs_v(i,k)  = qpixs_mid
    END IF
  END DO
END DO  ! nlev

! Set bwater=.true. on points where water will condense rather than
! ice.

! DEPENDS ON: flag_wet
CALL flag_wet(npnts,npnts,nlev,th,exner_layer_centres,bwater)


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

DO k = 2,nlev-1
!-----------------------------------------------------------------------
! Initialise environment variables
! NB These variable are only used by layer_cn.
!-----------------------------------------------------------------------
  DO i = 1,npnts
    thek(i)   = th(i,k)
    qek(i)    = q(i,k)
    qsek(i)   = qse(i,k)
    thekp1(i) = th(i,k+1)
    qekp1(i)  = q(i,k+1)
    qsekp1(i) = qse(i,k+1)
    bwk(i)    = bwater(i,k)
    bwkp1(i)  = bwater(i,k+1)
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    zk(i)     = z_theta(i,k)
    zkp12(i)  = z_rho(i,k+1)
    zkp1(i)   = z_theta(i,k+1)
    rhum(i)   = q(i,k) / qse(i,k)
  END DO

!-----------------------------------------------------------------------
! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
!-----------------------------------------------------------------------
  DO i = 1,npnts
    IF ( .NOT. bconv(i)) THEN
      expi(i)     = exner_layer_centres(i,k)
      bgmk(i)     = .FALSE.
      depth(i)    = 0.0
      thpi(i)     = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi(i)
      qpi(i)      = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi(i)
      IF (l_q_interact) THEN
        qclp(i,k) = qcl(i,k)
        qcfp(i,k) = qcf(i,k)
      ELSE
        qclp(i,k) = 0.0
        qcfp(i,k) = 0.0
      END IF
      IF (l_mom_gk) THEN
        up(i,k)   = u(i,k)
        vp(i,k)   = v(i,k)
      END IF
    END IF  !not bconv
  END DO
  IF (l_tracer) THEN
    DO ktra=1,ntra
      DO i = 1,npnts
        IF ( .NOT. bconv(i)) THEN
          trap(i,k,ktra)  = tracer(i,k,ktra)
        END IF  !not bconv
      END DO    ! npnts
    END DO
  END IF

! Initialise binit i.e. set as no convection initialised in layer
! at start of this levels calculations
  DO i = 1,npnts
    binit(i) = .FALSE.
  END DO

!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
  CALL layer_cn(k,npnts,nlev,                                     &
                   mdet_md_on, md_ent_on,                         &
                   ntml,ntpar,                                    &
                   .FALSE.,.FALSE.,.FALSE.,                       &
                   bconv,bwk,bwkp1,                               &
                   exner_layer_boundaries,                        &
                   exner_layer_centres,                           &
                   p_layer_boundaries,p_layer_centres,            &
                   recip_pstar,entrain_coef,rhum,                 &
                   zk, zkp12, zkp1,                               &
                   thek, qek,qsek, thekp1,qekp1,qsekp1,           &
                   thpk,qpk ,qsat_lcl, ekm14,                     &
                   pkp1,delpkp1,exkp1,                            &
                   pk,delpk,delpkp12,exk,delexkp1,                &
                   delp_uv_k, delp_uv_kp1,                        &
                   ekp14,ekp34,amdetk)


! Set ekm14 for next pass through loop
  DO i = 1, npnts
    ekm14(i) = ekp14(i)
  END DO

  DO i = 1,npnts
    ! Maximum initial convective mass flux
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

    ! Minimum buoyancy for convection to start from level k
    eminds(i) = mparb * delpkp12(i) * recip_pstar(i)
  END DO

!-----------------------------------------------------------------------
! Initial test to check if convection is possible in layer k
!-----------------------------------------------------------------------
! convection is possible from layer k to k+1 if
! - the point was convecting (bconv = .T.) and did not terminate
!   in the previous layer  OR
! - or k > ntml+1 and the stability is low enough at k+1
!  (At levels above nbl this is the same as midtrig(i)=k)

  DO i = 1,npnts
    bcposs(i) = bconv(i) .OR.                                     &
                 (p_layer_centres(i,k) > mid_cnv_pmin .AND.       &
                 k > ntml(i) + 1 .AND.                            &
                 (( th(i,k) * (1.0 + c_virtual * q(i,k)) -        &
                 th(i,k+1) * (1.0 + c_virtual * q(i,k+1)) +       &
                 delthst + MAX(0.0,(q(i,k)-qse(i,k+1))) *         &
                 (lc/(cp * exkp1(i)))) > 0.)) .OR.                &
                 (k == ntml(i) .AND. ntml(i) == nbl-1)

  END DO  ! npnts

! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
  ncposs = 0
  DO i = 1,npnts
    IF (bcposs(i)) THEN
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    END IF
  END DO

!-----------------------------------------------------------------------
! Compress points where convection may occur
! NB This process is used to update some single level variables that are 
! defined on level k and kp1 by using the full field variables.
! NB The order in which the variables are compressed 
! is the same as the argument list for LIFT_PAR
!-----------------------------------------------------------------------
  IF (ncposs  >   0) THEN
    DO i = 1,ncposs
      !INTENT(IN) for lift_par/water_loading
      thek_c(i)     = th(index1(i),k)
      thekp1_c(i)   = th(index1(i),k+1)
      qek_c(i)      = q(index1(i),k)
      qekp1_c(i)    = q(index1(i),k+1)
      qclek_c(i)    = qcl(index1(i),k)
      qcfek_c(i)    = qcf(index1(i),k)
      qclekp1_c(i)  = qcl(index1(i),k+1)
      qcfekp1_c(i)  = qcf(index1(i),k+1)
      pk_c(i)       = pk(index1(i))
      pkp1_c(i)     = pkp1(index1(i))
      exkp1_c(i)    = exkp1(index1(i))
      thpk_c(i)     = thp(index1(i),k)
      qpk_c(i)      = qp(index1(i),k)
      qclpk_c(i)    = qclp(index1(i),k)
      qcfpk_c(i)    = qcfp(index1(i),k)
      ekp14_c(i)    = ekp14(index1(i))
      ekp34_c(i)    = ekp34(index1(i))
      bwk_c(i)      = bwater(index1(i),k)
      bwkp1_c(i)    = bwater(index1(i),k+1)
      !INTENT(IN) for water_loading only
      exk_c(i)      = exk(index1(i))
    END DO
  END IF  ! ncposs>0

!-----------------------------------------------------------------------
! 3.2  Lift parcel from layer k to layer k+1
!-----------------------------------------------------------------------
  IF (ncposs > 0) THEN
  
    CALL lift_par_6a(ncposs, thek_c, thekp1_c,                      &
                qek_c, qekp1_c, qclek_c, qcfek_c,                   &
                qclekp1_c, qcfekp1_c,                               &
                pk_c, pkp1_c, exkp1_c,                              &
                thpk_c, qpk_c, qclpk_c, qcfpk_c,                    &
                ekp14_c, ekp34_c,                                   &
                l_q_interact, bwk_c, bwkp1_c,                       &
                !Out
                bgmkp1_c, thpkp1_c, qpkp1_c,                        &
                qclpkp1_c, qcfpkp1_c,                               &
                Qlkp1_c, Qfkp1_c, Frezkp1_c)

!-----------------------------------------------------------------------
! Calculate the water loading for level k and k+1
!-----------------------------------------------------------------------
    CALL water_loading(ncposs, pk_c, exk_c, thek_c, thpk_c,         &
                     qclek_c, qcfek_c, qclpk_c, qcfpk_c,            &
                     watldek_c, watldpk_c)

    CALL water_loading(ncposs, pkp1_c, exkp1_c, thekp1_c, thpkp1_c, &
                     qclekp1_c, qcfekp1_c, qclpkp1_c, qcfpkp1_c,    &
                     watldekp1_c, watldpkp1_c)

! NEC compiler directive
!CDIR NODEP

!-----------------------------------------------------------------------
! Test if convection is starting from layer k
!-----------------------------------------------------------------------
    DO i = 1,ncposs ! Loop over points which may convect

! Calculate buoyancy (virt. potential temp.) of parcel in layer k and k+1
      rbuoyk_c(i)   = thpk_c(i) * (1.0 + c_virtual *qpk_c(i)        &
                    - watldpk_c(i))                                 &
                    - thek_c(i) * (1.0 + c_virtual *qek_c(i)        &
                    - watldek_c(i))                                  

      rbuoykp1_c(i) = thpkp1_c(i) * (1.0 + c_virtual *qpkp1_c(i)    &
                    - watldpkp1_c(i))                               &
                    - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i)    &
                    - watldekp1_c(i))                                  

      ! Allow parcel to convect at midtrig or above if
      !1. it has not convected in the past (bconv=.F.,bterm=.F.) AND
      !2. k >= midtrig AND
      !3. the buoyancy of parcel at k+1 > min buoyancy required for
      !   convection + XSBMIN
      ! Flag convecting points with logical array bconv

      IF ( k > det_lev(index1(i)) .AND.                            &
        .NOT. bconv(index1(i)) .AND. .NOT. bterm(index1(i)) .AND. &
       ( k  >   midtrig(index1(i)).OR.                              &
       (k == ntml(index1(i)).AND.ntml(index1(i)) == nbl-1)) .AND.   &
             rbuoykp1_c(i)  >   (eminds(index1(i)) + xsbmin))THEN
        bconv(index1(i)) = .TRUE.
        start_lev(index1(i)) = k
        binit(index1(i)) = .TRUE.
        l_mid_all(index1(i)) = .TRUE. ! mid level convection
      END IF

      ! Set parcel mass flux (UM Documentation paper 27, section 1.5)
      ! If mass flux out of the initial layer is greater than the mass flux
      ! of the layer over the timestep then limit mass flux to mass of layer.

      IF (bconv(index1(i)).AND.k == start_lev(index1(i))) THEN
        flxk_c(i) = 1.0e-3*pstar(index1(i)) * (d_mid + c_mid *      &
                      pstar(index1(i)) * ((rbuoykp1_c(i) - xsbmin)     &
                       / delpkp12(index1(i))))
        IF (flxk_c(i)  >   flxmax(index1(i))) THEN
          flxk_c(i) = flxmax(index1(i))
        END IF

        ! Write compressed mass flux back to full array
        flx(index1(i),k) = flxk_c(i)
        
        ! Apply the parcel perturbation increments
        dthbydt(index1(i),k) = dthbydt(index1(i),k)              &
                               - flxk_c(i)/delpk(index1(i))      &
                               *(thpk_c(i)-thek_c(i))
        dqbydt(index1(i),k)  = dqbydt(index1(i),k)               &
                               - flxk_c(i)/delpk(index1(i))      &
                               *(qpk_c(i)-qek_c(i))
        
        ! Set mixing detrainment rate at first convecting level to zero
        amdetk(index1(i))      = 0.0

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))    = flxk_c(i)
        flxmax_init(index1(i)) = flxmax(index1(i))

      END IF

    END DO  !ncposs loop
  END IF  !ncposs>0

! Calculate number of points which are convecting  (nconv)
! set compression indices (index2).
  nconv = 0
  DO i = 1,ncposs
    IF (bconv(index1(i))) THEN
      nconv         = nconv + 1
      index2(nconv) = i
    END IF
  END DO

!-----------------------------------------------------------------------
! Second compression to form arrays of length nconv to be passed
! into CONVEC2.
! NB This process is used to update some single level variables that are 
! defined on level k and kp1 by using the full field variables.
! NB The order in which the variables are compressed 
! is the same as the argument list for CONVEC2.
!-----------------------------------------------------------------------
  IF (nconv  >   0) THEN
    !Compression for INTENT(IN)
    DO i = 1,nconv
      start_lev_c2(i)   = start_lev(index1(index2(i)))
      pstar_c2(i)       = pstar(index1(index2(i)))
      pk_c2(i)          = pk_c(index2(i))
      pkp1_c2(i)        = pkp1_c(index2(i))
      delpk_c2(i)       = delpk(index1(index2(i)))
      delpkp1_c2(i)     = delpkp1(index1(index2(i)))
      delpkp12_c2(i)    = delpkp12(index1(index2(i)))
      delp_uv_k_c2(i)   = delp_uv_k(index1(index2(i)))
      delp_uv_kp1_c2(i) = delp_uv_kp1(index1(index2(i)))
      exk_c2(i)         = exk_c(index2(i))
      exkp1_c2(i)       = exkp1_c(index2(i))
      thek_c2(i)        = thek_c(index2(i))
      thekp1_c2(i)      = thekp1_c(index2(i))
      qek_c2(i)         = qek_c(index2(i))
      qekp1_c2(i)       = qekp1_c(index2(i))
      qclek_c2(i)       = qclek_c(index2(i))
      qclekp1_c2(i)     = qclekp1_c(index2(i))
      qcfek_c2(i)       = qcfek_c(index2(i))
      qcfekp1_c2(i)     = qcfekp1_c(index2(i))
      qsek_c2(i)        = qse(index1(index2(i)),k)
      qsekp1_c2(i)      = qse(index1(index2(i)),k+1)
      cflek_c2(i)       = cf_liquid(index1(index2(i)),k)
      cflekp1_c2(i)     = cf_liquid(index1(index2(i)),k+1)
      cffek_c2(i)       = cf_frozen(index1(index2(i)),k)
      cffekp1_c2(i)     = cf_frozen(index1(index2(i)),k+1)
      bcfek_c2(i)       = bulk_cf(index1(index2(i)),k)
      bcfekp1_c2(i)     = bulk_cf(index1(index2(i)),k+1)
      thpk_c2(i)        = thpk_c(index2(i))
      qpk_c2(i)         = qpk_c(index2(i))
      qclpk_c2(i)       = qclpk_c(index2(i))
      qcfpk_c2(i)       = qcfpk_c(index2(i))
      thpi_c2(i)        = thpi(index1(index2(i)))
      qpi_c2(i)         = qpi(index1(index2(i)))
      expi_c2(i)        = expi(index1(index2(i)))
      rbuoyk_c2(i)      = rbuoyk_c(index2(i))
      rbuoykp1_c2(i)    = rbuoykp1_c(index2(i))
      xsbmin_v_c2(i)    = xsbmin_v(index1(index2(i)),k)
      watldek_c2(i)     = watldekp1_c(index2(i))
      watldekp1_c2(i)   = watldekp1_c(index2(i))
      watldpk_c2(i)     = watldpkp1_c(index2(i))
      watldpkp1_c2(i)   = watldpkp1_c(index2(i))
      Qlkp1_c2(i)       = Qlkp1_c(index2(i))
      Qfkp1_c2(i)       = Qfkp1_c(index2(i))
      Frezkp1_c2(i)     = Frezkp1_c(index2(i))      
      ekp14_c2(i)       = ekp14_c(index2(i))
      ekp34_c2(i)       = ekp34_c(index2(i))
      amdetk_c2(i)      = amdetk(index1(index2(i)))
      flxk_c2(i)        = flx(index1(index2(i)),k)
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        uek_c2(i)       = u(index1(index2(i)),k)
        uekp1_c2(i)     = u(index1(index2(i)),k+1)
        vek_c2(i)       = v(index1(index2(i)),k)
        vekp1_c2(i)     = v(index1(index2(i)),k+1)
        upk_c2(i)       = up(index1(index2(i)),k)
        vpk_c2(i)       = vp(index1(index2(i)),k)
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          traek_c2(i,ktra)   = tracer(index1(index2(i)),k,ktra)
          traekp1_c2(i,ktra) = tracer(index1(index2(i)),k+1,ktra)
          trapk_c2(i,ktra)   = trap(index1(index2(i)),k,ktra)
        END DO
      END DO
    END IF
    DO i = 1,nconv
      bgmk_c2(i)        = bgmk(index1(index2(i)))
      bgmkp1_c2(i)      = bgmkp1_c(index2(i))
      bwk_c2(i)         = bwk_c(index2(i))
      bwkp1_c2(i)       = bwkp1_c(index2(i))
      blowst_c2(i)      = binit(index1(index2(i)))
      bland_c2(i)       = bland(index1(index2(i)))
    END DO
    !Compression for INTENT(INOUT)
    DO i = 1,nconv
      lcbase_c2(i)      = lcbase(index1(index2(i)))
      lctop_c2(i)       = lctop(index1(index2(i)))
      thpkp1_c2(i)      = thpkp1_c(index2(i))
      qpkp1_c2(i)       = qpkp1_c(index2(i))
      qclpkp1_c2(i)     = qclpkp1_c(index2(i))
      qcfpkp1_c2(i)     = qcfpkp1_c(index2(i))
      dthek_c2(i)       = dthbydt(index1(index2(i)),k)
      dqek_c2(i)        = dqbydt(index1(index2(i)),k)
      dqclek_c2(i)      = dqclbydt(index1(index2(i)),k)
      dqcfek_c2(i)      = dqcfbydt(index1(index2(i)),k)
      tcw_c2(i)         = tcw(index1(index2(i)))
      depth_c2(i)       = depth(index1(index2(i)))
      cclwp_c2(i)       = cclwp(index1(index2(i)))
      lcca_c2(i)        = lcca(index1(index2(i)))
      cape_c2(i)        = cape(index1(index2(i)))
      dcpbydt_c2(i)     = dcpbydt(index1(index2(i)))
      relh_c2(i)        = relh(index1(index2(i)))
      dptot_c2(i)       = dptot(index1(index2(i)))
      max_cfl_c2(i)     = max_cfl(index1(index2(i))) 
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv 
        eflux_u_ud_c2(i)= eflux_u_ud(index1(index2(i)))
        eflux_v_ud_c2(i)= eflux_v_ud(index1(index2(i)))
        duek_c2(i)      = dubydt(index1(index2(i)),k)
        dvek_c2(i)      = dvbydt(index1(index2(i)),k)
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
        END DO
      END DO
    END IF
    DO i = 1,nconv
      bterm_c2(i)       = .FALSE.
      blatent_c2(i)     = blatent(index1(index2(i)))
    END DO
    !Compression for INTENT(OUT)
    ! Most INTENT(OUT) variables do not need to be initialised because they
    ! are always set in CONVEC2 or if they are not set in CONVEC2 then they
    ! are not used. However, several of the cloud variables do need to be
    ! initialised.
    DO i = 1,nconv
      iccb_c2(i)        = iccb(index1(index2(i)))
      icct_c2(i)        = icct(index1(index2(i)))
      cca_2d_c2(i)      = cca_2d(index1(index2(i)))
    END DO

!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------
    
    CALL convec2_6a  (k, nconv, npnts, nlev, ntra, md_on, md_new_termc,      &
                      start_lev_c2,                                          &
                      pstar_c2, pk_c2, pkp1_c2, delpk_c2,                    &
                      delpkp1_c2, delpkp12_c2, delp_uv_k_c2, delp_uv_kp1_c2, &
                      exk_c2, exkp1_c2,                                      &
                      thek_c2, thekp1_c2, qek_c2, qekp1_c2,                  &
                      qclek_c2, qclekp1_c2, qcfek_c2, qcfekp1_c2,            &
                      qsek_c2, qsekp1_c2,                                    &
                      cflek_c2, cflekp1_c2,  cffek_c2,  cffekp1_c2,          &
                      bcfek_c2,  bcfekp1_c2,                                 &
                      thpk_c2, qpk_c2, qclpk_c2, qcfpk_c2,                   &
                      thpi_c2, qpi_c2, expi_c2,                              &
                      rbuoyk_c2, rbuoykp1_c2, xsbmin_v_c2,                   &
                      watldek_c2, watldekp1_c2, watldpk_c2, watldpkp1_c2,    &
                      Qlkp1_c2, Qfkp1_c2, Frezkp1_c2,                        &
                      ekp14_c2, ekp34_c2, amdetk_c2, flxk_c2,                &
                      uek_c2, uekp1_c2, vek_c2, vekp1_c2,                    &
                      upk_c2, vpk_c2,                                        &
                      traek_c2, traekp1_c2, trapk_c2,                        &
                      l_q_interact, l_mom_gk, l_tracer,                      &
                      bgmk_c2, bgmkp1_c2, bwk_c2,                            &
                      bwkp1_c2, blowst_c2, bland_c2,                         &
                      ! In/out
                      lcbase_c2, lctop_c2,                                   &
                      thpkp1_c2, qpkp1_c2, qclpkp1_c2, qcfpkp1_c2,           &
                      dthek_c2, dqek_c2, dqclek_c2, dqcfek_c2,               &
                      tcw_c2, depth_c2, cclwp_c2, lcca_c2,                   &
                      cape_c2, dcpbydt_c2, relh_c2, dptot_c2, max_cfl_c2,    &
                      eflux_u_ud_c2, eflux_v_ud_c2,                          &
                      duek_c2, dvek_c2,                                      &
                      dtraek_c2,                                             &
                      bterm_c2, blatent_c2,                                  &
                      ! Out
                      iccb_c2, icct_c2,                                      &
                      dcflek_c2, dcffek_c2, dbcfek_c2,                       &
                      dthekp1_c2, dqekp1_c2, dqclekp1_c2, dqcfekp1_c2,       &
                      dcflekp1_c2, dcffekp1_c2, dbcfekp1_c2,                 &
                      prekp1_c2, deltak_c2, flxkp12_c2, flxkp1_c2,           &
                      cca_2d_c2, ccwkp1_c2,                                  &
                      upkp1_c2, vpkp1_c2,                                    &
                      duekp1_c2, dvekp1_c2,                                  &
                      trapkp1_c2,                                            &
                      dtraekp1_c2)

  END IF ! nconv > 0

!-----------------------------------------------------------------------
! Decompression of compressed variables coming out of
! of CONVEC2.
! NB The order in which the variables are decompressed 
! is the same as the the argument list for CONVEC2.
!-----------------------------------------------------------------------
  DO i = 1,npnts
    depth(i)      = 0.0
    bgmk(i)       = .FALSE.
    bterm(i)      = .FALSE.
  END DO

  IF (nconv  >   0) THEN
    !Decompression for INTENT(IN)
    DO i = 1,nconv
      bgmk(index1(index2(i)))         = bgmkp1_c2(i)
    END DO
    !Decompression for INTENT(INOUT)
    DO i = 1,nconv
      lcbase(index1(index2(i)))       = lcbase_c2(i)
      lctop(index1(index2(i)))        = lctop_c2(i)
      thp(index1(index2(i)),k+1)      = thpkp1_c2(i)
      qp(index1(index2(i)),k+1)       = qpkp1_c2(i)
      qclp(index1(index2(i)),k+1)     = qclpkp1_c2(i)
      qcfp(index1(index2(i)),k+1)     = qcfpkp1_c2(i)
      dthbydt(index1(index2(i)),k)    = dthek_c2(i)
      dqbydt(index1(index2(i)),k)     = dqek_c2(i)
      dqclbydt(index1(index2(i)),k)   = dqclek_c2(i)
      dqcfbydt(index1(index2(i)),k)   = dqcfek_c2(i)
      tcw(index1(index2(i)))          = tcw_c2(i)
      depth(index1(index2(i)))        = depth_c2(i)            
      cclwp(index1(index2(i)))        = cclwp_c2(i)
      lcca(index1(index2(i)))         = lcca_c2(i)
      cape(index1(index2(i)))         = cape_c2(i)
      dcpbydt(index1(index2(i)))      = dcpbydt_c2(i)
      relh(index1(index2(i)))         = relh_c2(i)
      dptot(index1(index2(i)))        = dptot_c2(i)
      max_cfl(index1(index2(i)))      = max_cfl_c2(i)
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
        eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
        dubydt(index1(index2(i)),k)   = duek_c2(i)
        dvbydt(index1(index2(i)),k)   = dvek_c2(i)
      END DO
    END IF
    IF  (l_tracer) THEN
      DO i = 1,nconv
        DO ktra = 1,ntra
          dtrabydt(index1(index2(i)),k,ktra)    = dtraek_c2(i,ktra)
        END DO
      END DO
    END IF
    DO i = 1,nconv
      bterm(index1(index2(i)))        = bterm_c2(i)    
      blatent(index1(index2(i)))      = blatent_c2(i)
    END DO
    !Decompression for INTENT(OUT)
    DO i = 1,nconv
      iccb(index1(index2(i)))         = iccb_c2(i)
      icct(index1(index2(i)))         = icct_c2(i)
      dcflbydt(index1(index2(i)),k)   = dcflek_c2(i)
      dcffbydt(index1(index2(i)),k)   = dcffek_c2(i)
      dbcfbydt(index1(index2(i)),k)   = dbcfek_c2(i)
      dthbydt(index1(index2(i)),k+1)  = dthekp1_c2(i)
      dqbydt(index1(index2(i)),k+1)   = dqekp1_c2(i)
      dqclbydt(index1(index2(i)),k+1) = dqclekp1_c2(i)
      dqcfbydt(index1(index2(i)),k+1) = dqcfekp1_c2(i)
      dcflbydt(index1(index2(i)),k+1) = dcflekp1_c2(i)
      dcffbydt(index1(index2(i)),k+1) = dcffekp1_c2(i)
      dbcfbydt(index1(index2(i)),k+1) = dbcfekp1_c2(i)
      precip(index1(index2(i)),k+1)   = prekp1_c2(i)
      !dummy                          = deltak_c2(i)   
      !dummy                          = flxkp12_c2(i)
      flx(index1(index2(i)),k+1)      = flxkp1_c2(i)
      cca_2d(index1(index2(i)))       = cca_2d_c2(i)
      ccw(index1(index2(i)),k+1)      = ccwkp1_c2(i)
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        up(index1(index2(i)),k+1)     = upkp1_c2(i)
        vp(index1(index2(i)),k+1)     = vpkp1_c2(i)
        dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
        dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
      END DO
    END IF
    IF  (l_tracer) THEN
      DO i = 1,nconv
        DO ktra = 1,ntra
          trap(index1(index2(i)),k+1,ktra)      = trapkp1_c2(i,ktra)
          dtrabydt(index1(index2(i)),k+1,ktra)  = dtraekp1_c2(i,ktra)
        END DO
      END DO
    END IF

  END IF     ! nconv >0

!-----------------------------------------------------------------------
! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
!      removed by convection over timescale cape_timescale.
!-----------------------------------------------------------------------

! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (npnts) with index_nterm

  nterm = 0
  DO i = 1,npnts
    IF (bterm(i)) THEN
      nterm = nterm + 1
      index_nterm(nterm) = i
    END IF
  END DO

  IF (nterm >  0) THEN

    SELECT CASE (cape_opt)

    ! Default 4a convection scheme - RH-based CAPE closure
    CASE(0)

! NEC compiler directive
!CDIR NODEP
      DO j = 1,nterm
        i = index_nterm(j)
        IF (dcpbydt(i)  >   0.0) THEN
          cape_ts_new(i) =                                              &
                  MIN(MAX(900.0*(1.0 - relh(i)/dptot(i))/0.1,60.0)      &
                                                      ,cape_timescale)
        END IF
      END DO

    ! Modified 4a convection scheme - RH-based CAPE closure
    ! timescale limited to timestep
    CASE(1)

! NEC compiler directive
!CDIR NODEP
      DO j = 1,nterm
        i = index_nterm(j)
        IF (dcpbydt(i)  >   0.0) THEN
          ! Only use RH cape if thickness of convective layer is greater 
          ! than 150hPa
          IF (p_layer_centres(i,start_lev(i))-                    &
               p_layer_centres(i,k)  >   15000.0)  THEN

            cape_ts_new(i) =                                            &
                    MIN(MAX(cape_timescale*(1.0-relh(i)/dptot(i))/0.4   &
                  ,timestep)                                            &
                  ,cape_timescale)

          ELSE
            cape_ts_new(i) = cape_timescale
          END IF
        END IF  ! dcpbydt > 0
      END DO

    ! Fixed cape timescale
    CASE(2)

! NEC compiler directive
!CDIR NODEP
      DO j = 1,nterm
        i = index_nterm(j)
        IF (dcpbydt(i)  >   0.0) THEN
          cape_ts_new(i) =  cape_timescale
        END IF
      END DO

    ! w based cape closure; if w_cape_limit > 1000. reverts to cape_timescale
    CASE(3)

      ! Initialise array to cape timescale and then alter as required.
      cape_ts_new(:) =  cape_timescale

      IF ( w_cape_limit < 1000.0 ) THEN
      !  This section includes test on w_max

!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            ! new denominator introduced at vn6.6
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) =   cape_timescale * w_cape_limit/           &
                      (w_cape_limit+ (w_max(i)-w_cape_limit)*wcape_fac)
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF !  w_max(i) > w_cape_limit
          END IF  ! dcpbydt > 0
        END DO
      END IF  ! w_cape_limit

    ! Option 4 - a w based cape closure; Grid-box area scaled CAPE closure
    CASE(4)

      IF ( w_cape_limit < 1000.0 ) THEN
      !  This section includes test on w_max

!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            ! new denominator introduced at vn6.6
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
            ELSE
              cape_ts_new(i) = cape_timescale * cape(i) / cape_min +    &
                               cape_timescale * EXP(-cape(i) / cape_min)
            END IF !  w_max(i) > w_cape_limit
          END IF  ! dcpbydt > 0
        END DO
      ELSE
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN

            cape_ts_new(i) = cape_timescale * cape(i) / cape_min +      &
                             cape_timescale * EXP( - cape(i) /cape_min)

          END IF  ! dcpbydt > 0
        END DO
      END IF  ! w_cape_limit

    ! Option 5 - a w based cape closure; Unavailable from UMUI
    CASE(5)

      IF ( w_cape_limit < 1000.0 ) THEN
      !  This section includes test on w_max

!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
            ELSE
              IF ( relh(i) / dptot(i) >= 0.75 ) THEN
                cape_ts_new(i) = cape_timescale *                        &
                                    ( 0.2373 / (relh(i) / dptot(i))**5)
              ELSE
                cape_ts_new(i) = cape_timescale
              END IF ! relh(i) / dptot(i) >= 0.75
            END IF !  w_max(i) > w_cape_limit
          END IF  ! dcpbydt > 0
        END DO
      ELSE
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            IF ( relh(i) / dptot(i) >= 0.75 ) THEN
              cape_ts_new(i) = cape_timescale *                         &
                                  ( 0.2373 / (relh(i) / dptot(i))**5)
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF ! relh(i) / dptot(i) >= 0.75
          END IF  ! dcpbydt > 0
        END DO
      END IF  ! w_cape_limit

    ! RH and w based CAPE option
    ! Expects a sensible w_cape_limit or will do nothing
    CASE(6)

      rh_test = 0.60         ! critical RH value
      rh_fac  = 1.0/ (1.0 - rh_test)

      IF ( w_cape_limit < 1000.0 ) THEN
      !  This section includes test on w_max
!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            ! work out any reduction in cape_timescale due to RH
            ! linearly falls to 1/2 given cape_timescale for RH above rh_test
            IF ( relh(i) / dptot(i) >= rh_test ) THEN
              cape_ts_new(i) = cape_timescale*0.5*                        &
                               (1.0+(1.0-(relh(i) / dptot(i)))*rh_fac)
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF
            ! Further reduction if w_max above critical value
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) = cape_ts_new(i) * w_cape_limit/        &
                       (w_cape_limit + (w_max(i)-w_cape_limit)*wcape_fac)

            END IF
            ! Limit CAPE timescale to convective timestep
            cape_ts_new(i) = MAX(cape_ts_new(i), timestep)

          END IF  ! dcpbydt > 0
        END DO
      END IF  ! w_cape_limit

    END SELECT        ! cape_opt

    ! Use new cape timescale calculated above

! NEC compiler directive
!CDIR NODEP
    DO j = 1,nterm
      i = index_nterm(j)
      IF (dcpbydt(i)  >   0.0) THEN
        flx_init_new(i) = flx_init(i)*cape(i)/(cape_ts_new(i)*dcpbydt(i))

        IF (flx_init_new(i)  >   flxmax_init(i)) THEN
          flx_init_new(i) = flxmax_init(i)
        END IF

        ! Scale max_cfl with cape scale

        max_cfl(i) = max_cfl(i) * flx_init_new(i) / flx_init(i)
      ELSE
        flx_init_new(i) = flx_init(i)
      END IF  ! dcpbydt > 0

    END DO

! Work out scaled mass flux needed to keep cfl ratio below limit.
! L_CAPE assumed to be true.

    DO j = 1,nterm
      i = index_nterm(j)
      max_cfl(i) = max_cfl(i) * timestep

      IF (max_cfl(i)  >   cfl_limit) THEN
        flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
        cfl_limited(i) = 1.0
      ELSE
        flx_init_new(i) = flx_init_new(i)
      END IF

      IF (flx_init_new(i)  >   flxmax_init(i)) THEN
        flx_init_new(i) = flxmax_init(i)
      END IF
      max_cfl(i) = 0.0

! Scale cloud fraction

      IF (flx_init_new(i)  >   0.0) THEN
        scale_f(i) = flx_init_new(i) / flx_init(i)
        cca_2d(i) = cca_2d(i) + 0.06 * LOG(scale_f(i))

        ! Check scaled cloud fraction not smaller than minimum value
        ! (2.0E-5) or greater than unity.

        cca_2d(i) = MAX(2.0e-5,cca_2d(i))
        IF (cca_2d(i)  >   1.0) THEN
          cca_2d(i) = 1.0
        END IF
      END IF

    END DO  ! nterm

!-----------------------------------------------------------------------
! Check for false/true convection and reset variables appropriately
!-----------------------------------------------------------------------
    DO j = 1,nterm
      i = index_nterm(j)
      IF ( (flx_init_new(i)  <=  minflx)                    &
                  .OR. ( (icct(i)-iccb(i)) <= 3 )           &
                  .OR. ( .NOT. blatent(i) ) ) THEN
        ! False convection diagnosed if:
        ! - the new initial mass flux is less than zero 
        ! - or the convecting layer it too thin 
        ! - or the parcel has never released latent heat
        ! 3d variables are reset below by setting scale_f to zero.
        flx_init(i)     = 0.0
        flx_init_new(i) = 0.0
        scale_f(i)      = 0.0
        cca_2d(i)       = 0.0
        iccb(i)         = 0
        icct(i)         = 0
        tcw(i)          = 0.0
        cclwp(i)        = 0.0
        lcca(i)         = 0.0
        lctop(i)        = 0
        lcbase(i)       = 0
        l_mid_all(i)    = .FALSE.
        kterm_mid(i)    = 0
        kterm(i)        = 0
        det_lev(i)      = 0
      ELSE
        ! True convection
        kterm_mid(i)    = k
        kterm(i)        = k
        det_lev(i)      = k+1
      END IF
    END DO  ! nterm

!-----------------------------------------------------------------------
! Apply closure and cfl scaling
!-----------------------------------------------------------------------
    DO kt = 2, k+1
      DO j = 1,nterm
        i = index_nterm(j)
        IF (kt  >=  start_lev(i)) THEN
          dthbydt(i,kt)   = dthbydt(i,kt)  * scale_f(i)
          dqbydt(i,kt)    = dqbydt(i,kt)   * scale_f(i)
          dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
          dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
          dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
          dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
          dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
          IF (l_mom_gk) THEN
            dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
            dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
          END IF
          IF (l_tracer) THEN
            DO ktra = 1,ntra
              dtrabydt(i,kt,ktra) =dtrabydt(i,kt,ktra)*scale_f(i)
            END DO
          END IF

          flx(i,kt)    = flx(i,kt)    * scale_f(i)
          precip(i,kt) = precip(i,kt) * scale_f(i)

        END IF ! kt>start_lev and flx_init_new >0
      END DO  ! j loop
    END DO  ! kt loop

!-----------------------------------------------------------------------
! 3.5  Downdraft calculation - on all points where convection is
!      terminating.
!
!      Subroutine DD_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

    IF (l_eman_dd) THEN

      ! save termination level

      DO j=1,nterm
        i = index_nterm(j)
        kterm_mid(i) = kterm(i)
      END DO

    ELSE        ! original down draughts

      npossdd = 0
      nnodd = 0
      DO j = 1,nterm
        i = index_nterm(j)
        tempnum = 0.0
        IF (iccb(i)  >   0) THEN
          deltap_cld = p_layer_centres(i,iccb(i))                   &
                               - p_layer_centres(i,k)
          DO kt = iccb(i), k+1
            tempnum = tempnum + precip(i,kt)
          END DO
        ELSE
          deltap_cld = 0.0
        END IF

! Downdraughts possible if pressure thickness of convective
! cloud (deltap_cld) is greater than 15000m, the point is saturated
! and the precip. in the layer is greater than a threashold
! value (1E-12).


        IF (deltap_cld  >   15000.0 .AND. bgmk(i) .AND.             &
                                     tempnum  >   1e-12) THEN
          npossdd = npossdd + 1
          index_possdd(npossdd) = i
        ELSE
          nnodd = nnodd + 1
          index_nodd(nnodd) = i
        END IF
      END DO  ! nterm loop

      ! If some downdraughts are possible (npossdd > 0) then call
      ! downdraught code

      IF (npossdd  >   0) THEN

! DEPENDS ON: dd_call_6a
        CALL dd_call_6a(npnts, npossdd, k, nlev, trlev, ntra        &
            ,iccb, icct, index_possdd                               &
            ,l_tracer                                               &
            ,bwater(1,2)                                            &
            ,exner_layer_centres, exner_layer_boundaries            &
            ,p_layer_centres, p_layer_boundaries, pstar             &
            ,recip_pstar, timestep , cca_2d                         &
            ,thp(1,1), qp(1,1), th(1,1), q(1,1),qse                 &
            ,trap, tracer, flx(1,1)                                 &
            ,precip(1,1), dthbydt(1,1), dqbydt(1,1), dtrabydt       &
            ,rain, snow ,rain_3d, snow_3d                           &
            ,dwn_flux, entrain_dwn, detrain_dwn)

      END IF


      ! Surface precipitation calculation for points where downdraught not
      ! possible.

      IF (nnodd  >   0) THEN
! DEPENDS ON: evap_bcb_nodd
        CALL evap_bcb_nodd(npnts,nnodd,k,iccb,index_nodd,         &
                       bwater(1,2),timestep,                      &
                       pstar,p_layer_centres,p_layer_boundaries,  &
                       exner_layer_centres,exner_layer_boundaries,&
                       th, q, qse, cca_2d,                        &
                       precip, dthbydt, dqbydt,                   &
                       rain, snow, rain_3d, snow_3d)

      END IF

    END IF  ! Emanuel test


! If convection has terminated write cape to diagnostic output
! variable (cape_out). Note if more than one lot of mid-level
! convection or if mid above deep or shallow then stores highest
! convection results only.
! Zero integrals as convection terminated.
! NB convection may start again at the
! same locations at higher levels in the atmosphere.
! NEC compiler directive
!CDIR NODEP
    DO j = 1,nterm
      i=index_nterm(j)
      cape_out(i)   = cape(i)
      dcpbydt(i)    = 0.0
      cape(i)       = 0.0
      bconv(i)      = .FALSE.
      bterm(i)      = .FALSE.     ! reset bterm to false
      start_lev(i)  = 0
      relh(i)       = 0.0      ! needed for RH CAPE
      dptot(i)      = 0.0      !
      bgmk(i)       = .FALSE.
      blatent(i)    = .FALSE.
    END DO

  END IF        ! nterm >0

!-----------------------------------------------------------------------
! Write out entrainment, detrainment and half-level mass flux diagnostics.
! They will be scaled by the full level mass flux outside 
! of the level loop
!-----------------------------------------------------------------------
! Calculate fractional entrainment rate for level k.
  IF (flg_entr_up .OR. mid_cmt_opt == 1) THEN
    DO i = 1,nconv
      entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i))        &
               * (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)  &
               * (1.0 + ekp14_c2(i)))
    END DO
  END IF

! Calculate fractional detrainment rate for level k
  IF (flg_detr_up) THEN
    DO i = 1,nconv
      detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)              &
                      + deltak_c2(i) * (1.0 - amdetk_c2(i)))
    END DO
  END IF

! Calculate the half level mass flux for level k
! Only the scaling factor between full level and half levels is calculated
! here. This is scaled by the full level mass flux outside the level loop
  IF (flg_up_flx_half) THEN
    DO i =1,nconv
      up_flux_half(index1(index2(i)),k) = (1.0 - deltak_c2(i))      &
                * (1.0 - amdetk_c2(i)) * (1.0 + ekp14_c2(i))
    END DO
  END IF

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------
END DO

!-----------------------------------------------------------------------
! Write out updraught massflux diagnostics and scale the 
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
IF (flg_up_flx .OR. flg_mf_midlev) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux(i,k) = flx(i,k)
    END DO
  END DO
END IF

IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF

! Has mid-level convection gone to the top of the model ?

DO i=1,npnts
  IF (kterm_mid(i) == nlev-1 .AND. l_mid_all(i) ) THEN
    IF (printstatus >= prstatus_normal) THEN
      WRITE(6,'(a41,i8)') 'WARNING: mid_conv has reached nlev-1 at ',i 
! Extra info which may prove useful to debugging problem
      WRITE(6,'(a7,g16.8)')    ' rain  ',rain(i)
      WRITE(6,'(a7,g16.8)')    ' snow  ',snow(i)
      WRITE(6,'(a7,400g16.8)') ' theta ',(th(i,k),      k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' q     ',(q(i,k),       k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' qcl   ',(qcl(i,k),     k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' qcf   ',(qcf(i,k),     k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' mf    ',(flx(i,k),     k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' dq    ',(dqbydt(i,k),  k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' dqcl  ',(dqclbydt(i,k),k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' dqcf  ',(dqcfbydt(i,k),k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' dth   ',(dthbydt(i,k), k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' du    ',(dubydt(i,k),  k=1,nlev)
      WRITE(6,'(a7,400g16.8)') ' dv    ',(dvbydt(i,k),  k=1,nlev)
    END IF

    error_point = i   
    IF (lhook) CALL dr_hook('MID_CONV_6A',zhook_out,zhook_handle)
    RETURN        
  END IF
END DO


! (Was copied out of scaling if test to ensure these limits at all
! times, not just when cca_2d is scaled)

! Check scaled cloud fraction not smaller than minimum value
! (2.0E-5) or greater than unity.

DO i=1, npnts
  IF (cca_2d(i) > 0.0) THEN
    cca_2d(i) = MAX(2.0e-5, cca_2d(i))
    IF (cca_2d(i) > 1.0) THEN
      cca_2d(i) = 1.0
    END IF
  END IF
END DO

!-----------------------------------------------------------------------
! How many columns have mid-level convection
!-----------------------------------------------------------------------

nmid = 0
DO i=1,npnts
  IF (l_mid_all(i)) THEN
    nmid = nmid +1
    index1(nmid) = i
  END IF
END DO



!-----------------------------------------------------------------------
! Emanuel down draughts - one call for all points with mid-level conv
!-----------------------------------------------------------------------

IF (l_eman_dd) THEN

! how many points have mid level convection?
! What is the maximum termination level ?

   kterm_max = 2

   DO i=1,npnts
     IF (l_mid_all(i)) THEN

        IF (kterm_mid(i) >  kterm_max) THEN
           kterm_max = kterm_mid(i)
        END IF

     END IF
   END DO

   ! Call routine to compress to just those points with convection
   ! and call Emanuel downdraughts then expand back to full grid.

! DEPENDS ON: eman_cex
   CALL eman_cex(npnts, nmid, kterm_max, nlev, trlev ,ntra        &
,                  kterm_mid, l_mid_all, l_tracer                 &
,                  exner_layer_centres,exner_layer_boundaries     &
,                  p_layer_centres, p_layer_boundaries            &
,                  timestep, th, q, qse, tracer ,precip           &
,                  dthbydt, dqbydt, dtrabydt                      &
,                  rain, snow, dwn_flux)

END IF

!-----------------------------------------------------------------------
! 4.0  Diffusive Convective momentum tranport option
!-----------------------------------------------------------------------

IF (l_mom .AND. mid_cmt_opt == 1) THEN

! Maximum number of mid-level layers, required for holding cloud bases
! and top of each layer. Note each layer must be at least 2 levels and have at
! least one level between. The bottom model level is not used and don't expect
! convection in any stratospheric levels so trying estimating by dividing by 6.
! This should give a number bigger than required but not excessive.

  nmax_layers = nlev/6

  ! Only call if there are points with mid_level convection.
  IF (nmid > 0) THEN

! DEPENDS ON: mid_conv_dif_cmt 
    CALL mid_conv_dif_cmt(npnts, nmid, nlev, nmax_layers,         &
                       index1, l_mid_all,                         &
                       timestep,                                  &
                       u, v, r_theta, r_rho,                      &
                       z_theta, z_rho, rho, rho_theta,            &
                       p_layer_boundaries,                        &
                       flx, entrain_up,                           &
                       dubydt, dvbydt, uw_mid, vw_mid )

  END IF  ! test on mid points
END IF  ! test on mid_cmt_opt 1

! ---------------------------------------------------------------------
! 5.0 Energy correction calculation - involves integral over
!     whole column. Located here it will given a different result
!     to being called in glue.
!     Removed call as not correct for new dynamics
!     UM documentation paper 27 - section 12.
! ---------------------------------------------------------------------
! First work out which points convection has occurred at .



IF (nmid >  0) THEN

  IF (l_cv_conserve_check) THEN
      CALL cor_engy_6a(npnts,nmid,nlev,index1,r2rho_th,                      &
                    dr_across_th,exner_layer_centres,p_layer_boundaries,     &
                    dqbydt,dqclbydt,dqcfbydt,rain,snow,                      &
                    dthbydt)
  END IF

!-----------------------------------------------------------------------
! 6.0  Correct negative/very small humidities
!-----------------------------------------------------------------------
! only check columns where convection has occurred.

  DO j = 1,nmid
    i=index1(j)
    qminincolumn(j) = q(i,nlev)
  END DO
  DO k = 1,nlev-1
    DO j = 1,nmid
      i=index1(j)
      IF (q(i,k)  <   qminincolumn(j)) THEN
        qminincolumn(j) = q(i,k)
      END IF
    END DO
  END DO

  ! Ensure Q does not go below global allowed minimum (QMIN)

  DO j = 1,nmid
    qminincolumn(j)=MAX(qmin,qminincolumn(j))
  END DO

  ! Apply an artificial upwards flux from k-1 level to ensure Q
  ! remains above minimum value in the column.

  DO k = nlev,2,-1
    DO j = 1,nmid
      i=index1(j)
      IF (dqbydt(i,k) /= 0.0) THEN
        temp1(j)=q(i,k) + dqbydt(i,k) * timestep

        IF (temp1(j)  <   qminincolumn(j)) THEN

          dqbydt(i,k-1) = dqbydt(i,k-1) -                         &
              ((qminincolumn(j) - q(i,k)) / timestep-dqbydt(i,k)) &
               * (r2rho_th(i,k)*dr_across_th(i,k))                &
               / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

          dqbydt(i,k) = (qminincolumn(j) - q(i,k)) / timestep
        END IF
      END IF
    END DO ! nmid loop
  END DO  ! nlev

  ! check no negative q at bottom
  k=1
    DO j = 1,nmid
      i=index1(j)
      temp1(j)=q(i,k) + dqbydt(i,k) * timestep

      IF (temp1(j)  <   0.0 .AND.                                       &
          printstatus >= prstatus_normal) THEN
        WRITE(6,'(a19,i6,a9,g26.18,a7,g26.18)')                         &
          ' negative q mid, i:', i,' q after ',temp1(j),' dq/dt ',dqbydt(i,k)
        WRITE(6,'(a9,400g16.8)') ' q inc   ',(dqbydt(i,k),  k=1,nlev)
        WRITE(6,'(a9,400g16.8)') ' qcl inc ',(dqclbydt(i,k),k=1,nlev)
        WRITE(6,'(a9,400g16.8)') ' qcf inc ',(dqcfbydt(i,k),k=1,nlev)
      END IF
    END DO ! nmid loop


END IF    ! test on nmid

!-----------------------------------------------------------------------
! 7.0  3d convective cloud
!-----------------------------------------------------------------------
! Initialise output array

DO k = 1,nlev
  DO i = 1,npnts
    cca(i,k) = 0.0
  END DO
END DO

!-----------------------------------------------------------------------------
! 7.1 Calculate CCA fow MID-CONVECTION levels only
!-----------------------------------------------------------------------------

n_mdcld      = 0    ! Initialise no of points with mid-cld
n_mdcld_mult = 0    ! Initialise no of points with multiple mid-cld

DO i=1, npnts
  IF ((l_mid_all(i)) .AND. (lcbase(i) /= 0)) THEN

    ! Mid-level cloud present
    ! Create index of gridpoints with single/multiple mid-level
    ! cloud

    IF (lcbase(i) /= iccb(i)) THEN

      n_mdcld_mult       = n_mdcld_mult + 1
      dum2(n_mdcld_mult) = i

    ELSE IF (lcbase(i) == iccb(i)) THEN

      n_mdcld       = n_mdcld + 1
      dum1(n_mdcld) = i
    END IF
  END IF
END DO

!---------------------------------------------------------------
! Calculate CCA_2D of Mid Cloud
!---------------------------------------------------------------
SELECT CASE (cca2d_md_opt)
  CASE(srf_precip)
    DO i=1, npnts
      IF ((l_mid_all(i)) .AND. (lcbase(i) /= 0)) THEN
        !-----------------------------------------------------------
        ! CCA_2D based on surface precipitation rate from mid
        ! cloud
        !-----------------------------------------------------------

        ! NOTE: at present the a_ and b_ parameters for LAND and SEA
        !       are equal, so there will be no difference between
        !       land and sea points.

        IF ((rain(i) + snow(i)) > 0.0) THEN

          IF (bland(i)) THEN
            ! Land point
            tempnum = a_land + b_land                                 &
                    * ALOG(86400.0 * (rain(i)+snow(i)))
          ELSE
            ! Sea point
            tempnum = a_sea + b_sea                                   &
                    * ALOG(86400.0 * (rain(i)+snow(i)))
          END IF

          cca_2d(i) = MAX(2.0e-5, tempnum)


          ! Grab lowest cca value before any tuning occurs
          ! This will overwrite lcca in ni_conv_ctl only if neither
          ! shallow or deep have occured.  This is under a switch in
          ! the 4a scheme.
          !
          ! NOTE: Downdraughts & Evaporation still being fed cca_2d
          !       derived from TCW, This issue may require further
          !       investigation.
          lcca(i) = cca_2d(i)

        END IF

      END IF      ! mid-level cloud present
    END DO      ! npnts

  CASE(total_condensed_water)
      ! cca_2d left unchanged from code, which is based on
      ! TCW (Total Condensed Water) (This is a rate)

END SELECT      ! cca2d_md_opt

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv

!-------------------------------------------------------------------
! 7.2 Apply CCA_3D
!-------------------------------------------------------------------
! A value of cca_2d has been calculated, though we need to know
! which levels to apply it to. So we identify which layers require
! CCA based on layers (above shallow/deep, if they have occurred)
! with non-zero ccw.
!-------------------------------------------------------------------

!===================================================================
! Single mid-level cloud bank gridpoints
!===================================================================

IF (n_mdcld > 0) THEN

  ! Resize compressed arrays for single Mid-level events

  ALLOCATE (mdcldi(n_mdcld))

  ! Now we know the number and indexes of gridpoints with
  ! single/multiple mid-level cloud, we can now assign the
  ! compressed mid-level arrays to pass to the anvil scheme.

  ALLOCATE (iccb_md_c       (n_mdcld       ))
  ALLOCATE (icct_md_c       (n_mdcld       ))
  ALLOCATE (freeze_lev_md_c (n_mdcld       ))
  ALLOCATE (cca_2d_md_c     (n_mdcld       ))
  ALLOCATE (cca_md_c        (n_mdcld,  n_cca_lev))
  ALLOCATE (ccw_md_c        (n_mdcld,  nlev))
  ALLOCATE (z_theta_md_c    (n_mdcld,  nlev))
  ALLOCATE (z_rho_md_c      (n_mdcld,  nlev))
  ALLOCATE (p_lyr_bnds_md_c (n_mdcld,0:nlev))

  DO i=1, n_mdcld
    mdcldi          (i) = dum1(i)
    iccb_md_c       (i) = iccb       (mdcldi(i))
    icct_md_c       (i) = icct       (mdcldi(i))
    cca_2d_md_c     (i) = cca_2d     (mdcldi(i))
    freeze_lev_md_c (i) = freeze_lev (mdcldi(i))
  END DO

  DO k=1, nlev
    DO i=1, n_mdcld
      ccw_md_c        (i,k) = ccw     (mdcldi(i),k)
      z_theta_md_c    (i,k) = z_theta (mdcldi(i),k)
      z_rho_md_c      (i,k) = z_rho   (mdcldi(i),k)
      p_lyr_bnds_md_c (i,k) =                        &
                    p_layer_boundaries(mdcldi(i),k)
    END DO
  END DO

  DO k=1, n_cca_lev
    DO i=1, n_mdcld
      cca_md_c(i,k) = 0.0
    END DO
  END DO

  !-----------------------------------------------------------------
  ! End Resizing compressed arrays for single Mid-level events.
  ! There is only one mid-level cloud bank if a mid-event has
  ! occurred, so does not require checking to see if there are
  ! multiple mid-level cloud banks.
  !-----------------------------------------------------------------

  IF (l_anvil) THEN
    CALL calc_3d_cca                                                  &
      ( n_mdcld, n_mdcld, nlev, n_cca_lev, nbl, iccb_md_c, icct_md_c  &
      , p_lyr_bnds_md_c, freeze_lev_md_c, cca_2d_md_c, cca_md_c       &
      , z_theta_md_c, z_rho_md_c )
  ELSE

    !---------------------------------------------------------------
    ! Do not use anvil scheme and apply cca_2d_md_c to all cloud
    ! levels between mid-level base and top
    !---------------------------------------------------------------
    DO k=1, n_cca_lev
      DO i=1, n_mdcld
        IF ( k >= iccb_md_c(i) .AND. &
             k <= icct_md_c(i) ) THEN
          cca_md_c(i,k) = cca_2d_md_c(i)
        ELSE
          cca_md_c(i,k) = 0.0
        END IF
      END DO      ! i (n_mdcld)
    END DO      ! k (n_cca_lev)
  END IF     ! l_anvil

  !-----------------------------------------------------------------
  ! Merge cca_md/ccw_md to full cca array and scale ccw_md_c by
  ! ccw_md_knob
  !-----------------------------------------------------------------

  DO k=1, n_cca_lev
    DO i=1, n_mdcld
      cca(mdcldi(i),k) = cca(mdcldi(i),k) + cca_md_c(i,k)
    END DO      ! i (n_mdcld)
  END DO      ! k (n_cca_lev)

  ! l_dcpl_cld4pc2 is set to true for 5a scheme, so
  ! CCRad tuning knobs applied in glue_conv

  ! Deallocate compressed single mid-level cloud arrays

  DEALLOCATE (mdcldi)
  DEALLOCATE (iccb_md_c)
  DEALLOCATE (icct_md_c)
  DEALLOCATE (freeze_lev_md_c)
  DEALLOCATE (cca_2d_md_c)
  DEALLOCATE (cca_md_c)
  DEALLOCATE (ccw_md_c)
  DEALLOCATE (z_theta_md_c)
  DEALLOCATE (z_rho_md_c)
  DEALLOCATE (p_lyr_bnds_md_c)

END IF ! n_mdcld > 0

!===================================================================
! Multiple mid-level cloud bank gridpoints
!===================================================================

IF (n_mdcld_mult > 0) THEN
  ! Resize index array for gridpoints with multiple Mid-level cloud
  ! banks and apply indices

  ALLOCATE(mdcldi_mult(n_mdcld_mult))

  !-----------------------------------------------------------------
  ! Allocate arrays with multiple-mid level cloud gridpoints
  !-----------------------------------------------------------------
  ALLOCATE (iccb_md_mult_c       (n_mdcld_mult       ))
  ALLOCATE (icct_md_mult_c       (n_mdcld_mult       ))
  ALLOCATE (freeze_lev_md_mult_c (n_mdcld_mult       ))
  ALLOCATE (cca_2d_md_mult_c     (n_mdcld_mult       ))

  ! Lets allocate them with levels first since it will then be continuous in
  ! memory later on when we call calc_3d_cca.
  ALLOCATE (cca_md_mult_c        ( n_cca_lev, n_mdcld_mult))
  ALLOCATE (ccw_md_mult_c        (  nlev, n_mdcld_mult))
  ALLOCATE (z_theta_md_mult_c    (  nlev, n_mdcld_mult))
  ALLOCATE (z_rho_md_mult_c      (  nlev, n_mdcld_mult))
  ALLOCATE (p_lyr_bnds_md_mult_c (0:nlev, n_mdcld_mult))

  DO i=1, n_mdcld_mult
    mdcldi_mult          (i) = dum2(i)
    iccb_md_mult_c       (i) = 0
    icct_md_mult_c       (i) = 0
    freeze_lev_md_mult_c (i) = freeze_lev (mdcldi_mult(i))
    cca_2d_md_mult_c     (i) = cca_2d     (mdcldi_mult(i))

    DO k=1, nlev
      ccw_md_mult_c       (k,i) = ccw     (mdcldi_mult(i),k)
      z_theta_md_mult_c   (k,i) = z_theta (mdcldi_mult(i),k)
      z_rho_md_mult_c     (k,i) = z_rho   (mdcldi_mult(i),k)
      p_lyr_bnds_md_mult_c(k,i) =                             &
                        p_layer_boundaries(mdcldi_mult(i),k)
    END DO

    DO k=1, n_cca_lev
      cca_md_mult_c(k,i) = 0.0
    END DO

  ! For multiple mid-level cloud banks, increment up model levels
  ! through compressed ccw_md_c array. Enter CAL3DCCA on locating
  ! cloud/bases of multiple clouds.

    DO k=2, nlev

      !-------------------------------------------------------------
      ! Check for cloud base
      !-------------------------------------------------------------
      IF ((ccw_md_mult_c(k,i) >  0.0) .AND.                         &
          (iccb_md_mult_c(i)  == 0)   .AND.                         &
          (icct_md_mult_c(i)  == 0)) THEN
        iccb_md_mult_c(i) = k
      END IF

      !-------------------------------------------------------------
      ! Check for cloud top
      !-------------------------------------------------------------
      IF ((ccw_md_mult_c(k,i)   <= 0.0) .AND.                       &
          (ccw_md_mult_c(k-1,i) >  0.0) .AND.                       &
          (iccb_md_mult_c(i)    /= 0)   .AND.                       &
          (icct_md_mult_c(i)    == 0)) THEN
        icct_md_mult_c(i) = k-1
      END IF

      !-------------------------------------------------------------
      ! Check for for anvil if both a cloud base and top found
      !-------------------------------------------------------------
      IF (iccb_md_mult_c(i) /= 0 .AND.                               &
          icct_md_mult_c(i) /= 0) THEN

        ! Apply CCA to vertically continuous cloud
        IF (l_anvil) THEN
          ! Since we are performing this on each individual point we can
          ! pass down the array with level information in the first rank
          ! instead of the second rank.  This reduces the amount of temporary
          ! data which is required otherwise it has to step through each
          ! array by number of points to get the all level information for it.
          CALL calc_3d_cca                                           &
            ( 1, 1, nlev, n_cca_lev, nbl, iccb_md_mult_c(i)          &
            , icct_md_mult_c(i), p_lyr_bnds_md_mult_c(:,i)           &
            , freeze_lev_md_mult_c(i), cca_2d_md_mult_c(i)           &
            , cca_md_mult_c(:,i), z_theta_md_mult_c(:,i)             &
            , z_rho_md_mult_c(:,i) )
        ELSE

          ! Copy cca_2d_md_c to all levels from cloud base to cloud
          ! top of located cloud

          DO j=iccb_md_mult_c(i), icct_md_mult_c(i)
            cca_md_mult_c(j,i) = cca_2d_md_mult_c(i)
          END DO      ! j
        END IF      ! l_anvil

        iccb_md_mult_c(i) = 0
        icct_md_mult_c(i) = 0

      END IF      ! Test on iccb_md_c(i) and icct_md_c(i)
    END DO      ! k (nlev)

  !-----------------------------------------------------------------
  ! Merge cca_md/ccw_md to cca full array and scale ccw_md_c by
  ! ccw_md_knob
  !-----------------------------------------------------------------

    DO k=1, n_cca_lev
      cca(mdcldi_mult(i),k) = cca(mdcldi_mult(i),k)                   &
                            + cca_md_mult_c(k,i)
    END DO      ! k (n_cca_lev)
  END DO      ! i (n_mdcld_mult)

  ! l_dcpl_cld4pc2 is set to true for 5a scheme, so
  ! CCRad tuning knobs applied in glue_conv

  !-----------------------------------------------------------------
  ! Deallocate arrays
  !-----------------------------------------------------------------
  DEALLOCATE (mdcldi_mult)
  DEALLOCATE (iccb_md_mult_c)
  DEALLOCATE (icct_md_mult_c)
  DEALLOCATE (freeze_lev_md_mult_c)
  DEALLOCATE (cca_2d_md_mult_c)
  DEALLOCATE (cca_md_mult_c)
  DEALLOCATE (ccw_md_mult_c)
  DEALLOCATE (z_theta_md_mult_c)
  DEALLOCATE (z_rho_md_mult_c)
  DEALLOCATE (p_lyr_bnds_md_mult_c)

END IF      ! n_mdcld_mult

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('MID_CONV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE mid_conv_6a
END MODULE mid_conv_6a_mod
