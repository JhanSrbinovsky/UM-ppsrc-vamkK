! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Mid-level convection scheme
!
SUBROUTINE  mid_conv_4a(nbl,nlev,ntra,n_cca_lev,npnts,trlev,         &
                       bland, w_max, exner_layer_centres,         &
                       exner_layer_boundaries,                    &
                       l_calc_dxek, l_q_interact,                 &
                       l_tracer,midtrig,ntml,ntpar,freeze_lev,    &
                       pstar,p_layer_centres,p_layer_boundaries,  &
                       r_theta,r_rho,                             &
                       z_theta,z_rho,rho,rho_theta,               &
                       q,q1_sd,t1_sd,th,                          &
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
                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out, &
                       uw_mid,vw_mid, cfl_limited, error_point    &
                       )

!------------------------------------------------------------------------
! Purpose:
!   Mid level convection scheme - works on all points.
!
!   Called by GLUE_CONV.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP 3 v8.3 programming standards
!------------------------------------------------------------------------

USE water_constants_mod, ONLY: lc, lf

USE cv_run_mod, ONLY:                                             &
    l_mom, l_eman_dd, cape_opt, cape_min,                         &
    w_cape_limit, cape_timescale, mid_cmt_opt, mid_cnv_pmin,      &
    cca2d_md_opt, cca_md_knob, ccw_md_knob, l_anvil,              &
    l_dcpl_cld4pc2, cnv_wat_load_opt, l_ccrad

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, srf_precip,                            &
    a_land, a_sea, b_land, b_sea, xsbmin,                         &
    thpixs_mid, qpixs_mid, thpixs_deep, qpixs_deep,               &
    mparb, c_mid, d_mid, wcape_fac, delthst

USE cv_dependent_switch_mod, ONLY:                                &
    md_on, mdet_md_on, md_ent_on, md_sdet_on, md_new_termc,       &
    cape_ts_w

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                     &
    flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn

USE atmos_constants_mod, ONLY: cp, c_virtual

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE PrintStatus_mod
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

REAL, INTENT(IN)    :: exner_layer_centres(npnts,0:nlev) !Exner

REAL, INTENT(IN)    :: exner_layer_boundaries(npnts,0:nlev)
                                ! Exner at half level above
                                ! exner_layer_centres

LOGICAL, INTENT(IN) :: l_calc_dxek ! Switch for calculation of
                                   ! condensate increment

LOGICAL, INTENT(IN) :: l_q_interact ! Switch allows overwriting
                                    ! parcel variables when
                                    ! calculating condensate incr.

LOGICAL, INTENT(IN) :: l_tracer ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: midtrig(npnts) ! Lowest trigger level
                                      ! for convection

INTEGER, INTENT(IN) :: ntml(npnts) ! Top level of surface mixed
                                   ! layer defined relative to
                                   ! theta,q grid

INTEGER, INTENT(IN) :: ntpar(npnts) ! Top level of initial parcel
                                    ! ascent in BL scheme defined
                                    ! relative to theta,q grid

INTEGER, INTENT(IN) :: freeze_lev(npnts) ! freezing level

REAL, INTENT(IN)    :: pstar(npnts) ! Surface pressure (Pa)

REAL, INTENT(IN)    :: p_layer_centres(npnts,0:nlev) ! Pressure (Pa)

REAL, INTENT(IN)    :: p_layer_boundaries(npnts,0:nlev) ! Pressure
                                                        ! at half level above
                                                        ! p_layer_centres (Pa)

REAL, INTENT(IN)    ::     &
  z_theta(npnts,nlev)      &  ! height of theta levels (m)
, z_rho(npnts,nlev)        &  ! height of rho levels (m)
, r_theta(npnts,0:nlev)    &  ! radius of theta levels (m)
, r_rho(npnts,nlev)        &  ! radius of rho levels (m)
, rho(npnts,nlev)          &  ! Density on rho levels (kg/m3)
, rho_theta(npnts,nlev)       ! Density on theta levels (kg/m3)

REAL, INTENT(IN)    :: q(npnts,nlev) ! Model mixing ratio (kg/kg)

REAL, INTENT(IN)    :: q1_sd(npnts) ! Standard deviation of
                                    ! turbulent flucts. of layer 1 q (kg/kg)

REAL, INTENT(IN)    :: t1_sd(npnts) ! Standard deviation of
                                    ! turbulent flucts. of layer 1 temp. (K)

REAL, INTENT(IN)    :: th(npnts,nlev) !Model potential temperature (K)

REAL, INTENT(IN)    :: timestep ! Model timestep (s)

REAL, INTENT(IN)    :: u(npnts,nlev) !Model U field (m/s)

REAL, INTENT(IN)    :: v(npnts,nlev) !Model V field (m/s)

REAL, INTENT(IN) :: w_max(npnts)    !  max w in column
                                    ! for use in scale dependent
                                    ! cape timescale

REAL, INTENT(IN) :: recip_pstar(npnts) ! Reciprocal of pstar array

REAL,INTENT(IN) :: qse(npnts,nlev) ! Saturation mixing ratio of
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

REAL, INTENT(OUT) :: cape_out(npnts) ! Saved convective available
                                ! potential energy for diagnostic
                                ! output (J/kg)

REAL, INTENT(OUT) :: cclwp(npnts)
                                ! Condensed water path (kg/m^2)

REAL, INTENT(OUT) :: ccw(npnts,nlev)
                                ! Convective cloud liquid water
                                ! on model levels (kg/kg)

REAL, INTENT(OUT) :: cca(npnts,n_cca_lev)
                                ! Convective cloud amount
                                ! on model levels (0-1)

REAL, INTENT(OUT) :: dbcfbydt(npnts,nlev) ! Increments to
                                ! total cld volume due to
                                ! convection(/s)

REAL, INTENT(OUT) :: dcffbydt(npnts,nlev) ! Increments to ice
                                ! cloud volume due to convection (/s)

REAL, INTENT(OUT) :: dcflbydt(npnts,nlev) ! Increments to liq
                                ! cloud volume due to convection (/s)

REAL, INTENT(OUT) :: dqbydt(npnts,nlev) ! Increments to q due to
                                ! convection (kg/kg/s)

REAL, INTENT(OUT) :: dqcfbydt(npnts,nlev) ! Increments to ice
                                ! condensate due to convection
                                ! (kg/kg/s)

REAL, INTENT(OUT) :: dqclbydt(npnts,nlev) ! Increments to liq
                                ! condensate due to convection
                                ! (kg/kg/s)

REAL, INTENT(OUT) :: dthbydt(npnts,nlev) ! Increments to potential
                                ! temp. due to convection (K/s)

REAL, INTENT(OUT) :: dubydt(npnts,nlev+1) ! Increments to U due
                                ! to CMT (m/s2)

REAL, INTENT(OUT) :: dvbydt(npnts,nlev+1) ! Increments to V due
                                ! to CMT (m/s2)

REAL, INTENT(OUT) :: dtrabydt(npnts,nlev,ntra) !Increment to
                                ! tracer convection (kg/kg/s)

REAL, INTENT(OUT) :: detrain_up(npnts,nlev) ! Fractional
                                ! detrainment rate into updraughts (Pa/s)

REAL, INTENT(OUT) :: detrain_dwn(npnts,nlev) ! Fractional
                                ! detrainment rate into downdraughts (Pa/s)

REAL, INTENT(OUT) :: entrain_up(npnts,nlev) ! Fractional 
                                ! entrainment rate into updraughts (Pa/s)

REAL, INTENT(OUT) :: entrain_dwn(npnts,nlev) ! Fractional
                                ! entrainment rate into downdraughts (Pa/s)

INTEGER, INTENT(OUT) :: iccb(npnts) ! Convective cloud base level

INTEGER, INTENT(OUT) :: icct(npnts) ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(npnts) ! Lowest conv. cloud amt. (%)

INTEGER, INTENT(OUT) :: lcbase(npnts) ! Lowest conv. cloud base level

INTEGER, INTENT(OUT) :: lctop(npnts) ! Lowest conv. cloud top level

REAL, INTENT(OUT) :: rain(npnts) ! Surface convective rainfall (kg/m2/s)

REAL, INTENT(OUT) :: snow(npnts) ! Surface convective snowfall (kg/m2/s)

REAL, INTENT(OUT) :: rain_3d(npnts,nlev) ! Convective rainfall flux (kg/m2/s)

REAL, INTENT(OUT) :: snow_3d(npnts,nlev) ! Convective snowfall flux (kg/m2/s)

REAL, INTENT(OUT) :: up_flux(npnts,nlev) ! Updraught mass flux (Pa/s)

REAL, INTENT(OUT) :: up_flux_half(npnts,nlev)
                                ! Updraught mass flux
                                ! on rho levels (Pa/s)

REAL, INTENT(OUT) :: dwn_flux(npnts,nlev) ! Downdraught mass flux (Pa/s)

LOGICAL, INTENT(OUT) :: l_mid_all(npnts)
                                ! Points where mid level convection
                                ! occurs at some level
REAL,INTENT(INOUT) :: cca_2d(npnts) !2D convective cloud amount(%)

! Adaptive detrainment output variables

REAL, INTENT(INOUT) ::      &
  rbuoy_p_out(npnts,nlev) & !buoyancy excess
 ,the_out(npnts,nlev)     & !th_E in parcel routine
 ,thp_out(npnts,nlev)     & !th_P in parcel routine
 ,qe_out(npnts,nlev)      & !q_E in parcel routine
 ,qp_out(npnts,nlev)        !q_P in parcel routine

! CMT diagnostics - INOUT as initialised in glue_conv
REAL, INTENT(INOUT) ::    &
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
! Adaptive detrainment output variables

REAL :: rbuoy_p_here(npnts)       !buoyancy excess

REAL :: the_here(npnts)         !th_E in parcel routine

REAL :: thp_here(npnts)         !th_P in parcel routine

REAL :: qe_here(npnts)          !q_E in parcel routine

REAL :: qp_here(npnts)          !q_P in parcel routine

REAL :: rbuoy_p_old(npnts)      !buoyancy excess on previous
                                !level

REAL :: zk(npnts)               !heights for use in calc

REAL :: zkp12(npnts)            !of moist static energy

REAL :: zkp1(npnts)

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

REAL :: dqsthk(npnts)           ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k (kg/kg/K)

REAL :: dqsthkp1(npnts)         ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k+1 (kg/kg/K)

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
                                ! the initial convecting layer
                                ! in Pa/s)

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

LOGICAL :: l_shallow(npnts)     ! Dummy variable (=.F.)

LOGICAL :: l_mid(npnts)         ! Dummy variable (=.T.)

LOGICAL :: cumulus(npnts)       ! Dummy variable (=.F.)

LOGICAL :: bgmk(npnts)          ! Mask for points where parcel in
                                ! layer k is saturated

LOGICAL :: bwater(npnts,2:nlev) ! Mask for points at which
                                ! condensate is liquid

LOGICAL :: bwk(npnts)            !mask for liquid condensate on k
LOGICAL :: bwkp1(npnts)          !mask for liquid condensate on k+1

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

REAL :: expi(npnts)
                                ! Initial parcel exner pressure

REAL :: xpk(npnts,nlev)         ! Parcel cloud water (kg/kg)

REAL :: flx(npnts,nlev)         ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(npnts,nlev)    ! Minmum parcel buoyancy excess

REAL :: thpixs_v(npnts,nlev)    ! Theta parcel excess (K)

REAL :: qpixs_v(npnts,nlev)     ! Q parcel excess(kg/kg)

REAL :: dpmin(npnts)            ! work array for parcel excess cal


! PC2

REAL :: qclp(npnts,nlev)         ! Parcel liquid condensate mixing
                                 ! ratio in layer k (kg/kg)

REAL :: qcfp(npnts,nlev)         ! Parcel frozen condensate mixing
                                 ! ratio in layer k (kg/kg)

! Parameters


REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed


! CMT variables


INTEGER :: kterm(npnts)         ! Level index for termination of

REAL :: eflux_u_ud(npnts)       ! Vertical eddy flux of momentum
                                ! due to UD at top of layer
                                ! (Pa m/s2)

REAL :: eflux_v_ud(npnts)       ! Vertical eddy flux of momentum
                                ! due to UD at bottom of layer
                                ! (Pa m/s2)

REAL :: flxkp12(nlev,npnts)     ! Mass flux on half level (Pa/s)

LOGICAL :: l_mom_gk             ! true if Gregory-Kershaw CMT

! Cape scaling/closure variables

INTEGER :: start_lev(npnts)     !Level at which convection
                                !initiates

INTEGER :: start_lev3c(npnts)   ! PC2 Compressed convection
                                ! initiation level

INTEGER :: det_lev(npnts)       ! Level at which split final
                                ! detrainment last occurred

INTEGER :: nterm                ! No. of points where conv.
                                ! has terminated

INTEGER :: index_nterm(npnts)   ! Index for points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage

REAL :: scale_f(npnts)          ! scaling factor

REAL :: dthef(npnts)            ! Theta increment from convection
                                ! in model level at which split
                                ! final detrainment last occurred
                                ! (K/s)

REAL :: dqf(npnts)              ! Specific humidity increment
                                ! from convection in model level
                                ! at which split final detrainment
                                ! last occurred (kg/kg/s)

REAL :: dqclf(npnts)            ! As dqf but for qcl (kg/kg/s)
REAL :: dqcff(npnts)            ! As dqf but for qcf (kg/kg/s)
REAL :: dcflf(npnts)            ! As dqf but for cfl (/s)
REAL :: dcfff(npnts)            ! As dqf but for cff (/s)
REAL :: dbcff(npnts)            ! As dqf but for bcf (/s)

REAL :: duef(npnts)             ! As for dthef but for U increments (m/s2)

REAL :: dvef(npnts)             ! As for dthef but for V increments (m/s2)

REAL :: dtraef(npnts,ntra)      ! As for dthef but for tracer
                                ! increments (kg/kg/s)

REAL :: cape_ts_new(npnts)      ! Used as variable in RH-based
                                ! closure

REAL :: relh(npnts)             ! RH integral (average when
                                ! convection terminates)

REAL :: dptot(npnts)            ! Delta P integral


! Downdraught scheme variables

INTEGER :: npossdd              ! Max. no. of downdraughts possible

INTEGER :: nnodd                ! No. of downdraughts not possible

INTEGER :: index_possdd(npnts)  ! Index of downdraughts possible

INTEGER :: index_nodd(npnts)    ! Index of downdraughts not possible

REAL :: deltap_cld              ! Pressure thickness of convective
                                ! cloud (Pa)

! Local compressed arrays

LOGICAL :: bconv_c2(npnts)

LOGICAL :: bgmkp1_c(npnts), bgmkp1_c2(npnts) ! Mask for points
                                ! where parcel in layer k+1
                                ! is saturated

LOGICAL :: bwk_c(npnts), bwk_c2(npnts) ! bwater mask in layer k

LOGICAL :: bwkp1_c(npnts), bwkp1_c2(npnts) ! bwater mask in layer k+1

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

REAL :: flxk_c(npnts), flxk_c2(npnts) !Parcel mass flux in layer k (Pa/s)

REAL :: flxkp12_c2(npnts)       ! Half level mass flux (Pa/s)

REAL :: prekp1_c2(npnts)        ! Precip. from parcel as it rises
                                ! from layer k to k+1 (kg/m2/s)

REAL :: qpk_c(npnts), qpk_c2(npnts) ! Parcel mixing ratio in layer k(kg/kg)

REAL :: qpk(npnts)              !ad. entrain.

REAL :: qpkp1_c(npnts), qpkp1_c2(npnts) ! Parcel mixing ratio layer k +1 (kg/kg)

REAL :: qek_c(npnts), qek_c2(npnts) ! Env. mixing ratio in layer k (kg/kg)

REAL :: qek(npnts)             !for ad entrain.

REAL :: qekp1_c(npnts), qekp1_c2(npnts) ! Env. mixing ratio in
                                ! layer k+1 (kgkg-1)

REAL :: qekp1(npnts)             !for ad entrain.

REAL :: qsek_c2(npnts)          ! Saturation mixing ratio of
                                ! cld. env. in layer k (kg/kg)

REAL :: qsek(npnts)             !for ad entrain.

REAL :: qsekp1_c(npnts), qsekp1_c2(npnts) ! Saturation mixing
                                ! ratio of cld. env. in layer k+1 (kg/kg)

REAL :: qsekp1(npnts)             !for ad entrain.

REAL :: thek_c(npnts), thek_c2(npnts) ! Env. potential temp in layer k (K)

REAL :: thek(npnts)             !for ad entrain.

REAL :: thekp1_c(npnts), thekp1_c2(npnts) ! Env. potential temp i in layer k (K)

REAL :: thekp1(npnts)             !for ad entrain.

REAL :: thpk_c(npnts), thpk_c2(npnts) ! Parcel potential temp in layer k (K)

REAL :: thpk(npnts)             !for ad entrain.

REAL :: thpkp1_c(npnts), thpkp1_c2(npnts)! Parcel potential temp in layer k (K)

REAL :: traek_c(npnts,ntra), traek_c2(npnts,ntra) ! Tracer content
                                ! cld. env. in layer k (kgkg-1)

REAL :: traekp1_c(npnts,ntra), traekp1_c2(npnts,ntra) ! Tracer
                                ! content of cloud env.
                                ! in layer k+1 (kg/kg)

REAL :: trapk_c(npnts,ntra), trapk_c2(npnts,ntra) ! Tracer cont.
                                ! of parcel in layer k (kg/kg)

REAL :: trapkp1_c(npnts,ntra), trapkp1_c2(npnts,ntra) ! Tracer
                            ! cont.of parcel in layer k+1 (kg/kg)

REAL :: rbuoy_c(npnts), rbuoy_c2(npnts)! Buoyancy of parcel at k+1 (Kelvin)

REAL :: uek_c(npnts), uek_c2(npnts) ! Model U field on layer k (m/s)

REAL :: uekp1_c(npnts), uekp1_c2(npnts)! Model U field on layer k+1 (m/s)

REAL :: vek_c(npnts), vek_c2(npnts) ! Model V field on layer k (m/s)

REAL :: vekp1_c(npnts), vekp1_c2(npnts) ! Model V field on layer k+1 (m/s)

REAL :: upk_c(npnts), upk_c2(npnts) ! Parcel U in layer k
                                ! after entrainment (m/s)


REAL :: upkp1_c(npnts), upkp1_c2(npnts) ! Parcel U in layer k+1
                                ! after entrainment (m/s)

REAL :: vpk_c(npnts), vpk_c2(npnts) ! Parcel V in layer k
                                ! after entrainment (m/s)

REAL :: vpkp1_c(npnts), vpkp1_c2(npnts) ! Parcel V in layer k+1
                                ! after entrainment (m/s)

REAL :: xsqkp1_c(npnts), xsqkp1_c2(npnts) ! Excess water vapour
                                ! in parcel at k+1 (kg/kg)

REAL :: qclek_c(npnts), qclek_c2(npnts) ! Environment liquid
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qclekp1_c(npnts), qclekp1_c2(npnts) ! Environment liquid
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qcfek_c(npnts), qcfek_c2(npnts) ! Environment frozen
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qcfekp1_c(npnts), qcfekp1_c2(npnts) ! Environment frozen
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qclpk_c(npnts), qclpk_c2(npnts) ! Parcel liquid
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qclpkp1_c(npnts), qclpkp1_c2(npnts) ! Parcel liquid
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qcfpk_c(npnts), qcfpk_c2(npnts) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qcfpkp1_c(npnts), qcfpkp1_c2(npnts) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

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
REAL :: ccw_c2(npnts)
LOGICAL :: cumulus_c2(npnts)
REAL :: dcpbydt_c2(npnts)
REAL :: delexkp1_c2(npnts)
REAL :: delpk_c2(npnts)
REAL :: delpkp1_c2(npnts)
REAL :: delp_uv_k_c2(npnts)
REAL :: delp_uv_kp1_c2(npnts)
REAL :: depth_c2(npnts)
REAL :: dptot_c2(npnts)
REAL :: dqsthkp1_c2(npnts)
REAL :: dqsthk_c2(npnts)
REAL :: eflux_u_ud_c2(npnts)
REAL :: eflux_v_ud_c2(npnts)
REAL :: ekp14_c(npnts),ekp14_c2(npnts)
REAL :: ekp34_c(npnts),ekp34_c2(npnts)
REAL :: eminds_c(npnts)
REAL :: exk_c2(npnts)
REAL :: exkp1_c(npnts),exkp1_c2(npnts)
REAL :: expi_c2(npnts)
INTEGER :: icct_c2(npnts)
INTEGER :: iccb_c2(npnts)
INTEGER :: lctop_c2(npnts)
INTEGER :: lcbase_c2(npnts)
REAL :: lcca_c2(npnts)
LOGICAL :: l_shallow_c2(npnts)
LOGICAL :: l_mid_c2(npnts)
REAL :: max_cfl_c2(npnts)
REAL :: pk_c(npnts),pk_c2(npnts)
REAL :: pkp1_c(npnts),pkp1_c2(npnts)
REAL :: pstar_c2(npnts)
REAL :: q1_sd_c2(npnts)
REAL :: qpi_c2(npnts)
REAL :: qpixs_v_c2(npnts)
REAL :: rbuoy_p_here_c2(npnts)
REAL :: the_here_c2(npnts)
REAL :: thp_here_c2(npnts)
REAL :: qe_here_c2(npnts)
REAL :: qp_here_c2(npnts)
REAL :: rbuoy_p_old_c2(npnts)
REAL :: relh_c2(npnts)
REAL :: tcw_c2(npnts)
REAL :: thpi_c2(npnts)
REAL :: thpixs_v_c2(npnts)
REAL :: t1_sd_c2(npnts)
REAL :: xpk_c(npnts),xpk_c2(npnts)
REAL :: xsbmin_v_c2(npnts)

INTEGER ::        &
 nmid             & ! total number of points where mid-level
                    ! convection occurs
,nmax_layers        ! Maximum number of allow mid-level layers


! CCRad Variables local variables
INTEGER ::     &
  dum_iccb     &! Holding variables for identifying cld base and tops
, dum_icct     &! on gridpoints with more than 1 (possibly 2) mid-level
                ! convective events

, n_mdcld      &! No. of points with single mid-level cloud layers
, n_mdcld_mult  ! No. of points with multiple mid-level cloud layers


! Allocatable arrays, because we do not know how many gridpoints
! will contain mid-level until we test for them.
INTEGER ::     &
  dum1(npnts)  &
, dum2(npnts)

INTEGER, ALLOCATABLE :: &
  mdcldi(:)       &! Indices in full array of gridpoints with
                   ! single mid-level cloud layers
, mdcldi_mult(:)   ! Indices in full array of gridpoints
                   ! multiple mid-level cloud layers


! Compressed arrays for gridpoints which have one bank of
! mid-level cloud and multiple banks of mid-level cloud.
! Requires the use of compressed allocatable arrays because the
! location and number of gridpoints with mid-level has not been
! diagnosed yet.
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

! End CCRad Variables local variables

! Arrays required by Emanuel downdraughts

INTEGER ::            &
  kterm_mid(npnts)    & ! termination level of highest mid level
, kterm_max

REAL ::      &
  rh_test    &   ! critical RH value for CAPE_OPT=6
 ,rh_fac         ! factor for CAPE_OPT=6 calculation

! Loop counters

INTEGER :: i,j,k,ktra,kt,kk

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('MID_CONV_4A',zhook_in,zhook_handle)

error_point = 0
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

!initialise SCM diagnostics


! Initialise logicals

DO i = 1,npnts
!        blowst(i)    = .true.
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  cumulus(i)   = .FALSE.
  l_shallow(i) = .FALSE.
  l_mid(i)     = .TRUE.
  l_mid_all(i) = .FALSE.

! Zero bottom level mass flux array as never set by levels loop.

  flx(i,1) = 0.0
  kterm_mid(i) = 0  ! initialise kterm_mid array
END DO

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------

! Re-calculate XSBMIN and THPIXS constants based on layer
! thickness (Pa).

DO k = 1,nlev-1
  DO i = 1,npnts
    dpmin(i) = MIN( ((p_layer_centres(i,k) -                    &
                         p_layer_centres(i,k+1))/5000.),1.0)

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
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)

! DEPENDS ON: flag_wet
CALL flag_wet(npnts,npnts,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!-----------------------------------------------------------------------

DO k=1, nlev
  DO i=1, npnts
    precip(i,k)  = 0.0
    ccw(i,k)     = 0.0
    dthbydt(i,k) = 0.0
    dqbydt(i,k)  = 0.0
  END DO
END DO

DO k=1, n_cca_lev
  DO i=1, npnts
    cca(i,k) = 0.0
  END DO
END DO

IF (l_mom) THEN

  l_mom_gk = .TRUE.        ! use Gregory-Kershaw CMT

  DO k = 1,nlev+1
    DO i = 1,npnts
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

IF (l_tracer) THEN
  DO ktra = 1,ntra
    ! No need to initialise dtrabydt as done in glue_conv
    DO i = 1,npnts
      dtraef(i,ktra) = 0.0
    END DO
  END DO
END IF  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------

IF (flg_up_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux(i,k) = 0.0
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
      dwn_flux(i,k) = 0.0
    END DO
  END DO
END IF
! Now required in all cases
DO k = 1,nlev
  DO i = 1,npnts
    entrain_up(i,k) = 0.0
  END DO
END DO

IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------


IF (l_ccrad) THEN
  ! Zeroes values to remove mid-level dependence on shallow and deep
  ! cloud
  DO i = 1,npnts
    iccb   (i) = 0
    icct   (i) = 0
    tcw    (i) = 0.0
    cca_2d (i) = 0.0
  END DO
END IF      ! l_ccrad

DO i = 1,npnts
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
END DO

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics (now in glue)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling and closure calculations
!-----------------------------------------------------------------------

DO i = 1,npnts
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)      = 0.0
  cape_out(i)  = 0.0
  dcpbydt(i)   = 0.0
  max_cfl(i)   = 0.0
  det_lev(i)   = 0
  start_lev(i) = 0
  dthef(i)      = 0.0
  dqf(i)        = 0.0
  dqclf(i)      = 0.0
  dqcff(i)      = 0.0
  dcflf(i)      = 0.0
  dcfff(i)      = 0.0
  dbcff(i)      = 0.0
  duef(i)       = 0.0
  dvef(i)       = 0.0
  relh(i)       = 0.0
  dptot(i)      = 0.0
  cfl_limited(i) = 0.0
!-----------------------------------------------------------------------
! 2.6  Initialise eddy flux arrays for updraught
!-----------------------------------------------------------------------

  eflux_u_ud(i) = 0.0
  eflux_v_ud(i) = 0.0

!-----------------------------------------------------------------------
! 2.7  Initialise surface precipitation arrays
!-----------------------------------------------------------------------

  rain(i) = 0.0
  snow(i) = 0.0

END DO

!-----------------------------------------------------------------------
! 2.8 Initialise PC2 arrays
!-----------------------------------------------------------------------

DO k = 1,nlev
  DO i = 1,npnts
    dqclbydt(i,k)  = 0.0
    dqcfbydt(i,k)  = 0.0
    dbcfbydt(i,k)  = 0.0
    dcflbydt(i,k)  = 0.0
    dcffbydt(i,k)  = 0.0
  END DO
END DO

! Set rhum to an arbitrary value since it is not used in LAYER_CN
! for mid-level convection. Therefore value not level dependent.

DO i = 1,npnts
  rhum(i) = 0.0
END DO


!initialise ekm14
DO i =1, npnts
  ekm14(i) =0.0
END DO

!Initialise adaptive entrainment variables
!Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1

DO i = 1, npnts
  thek(i)=th(i,1)
  qek(i)=q(i,1)
  qsek(i)=qse(i,1)
  thpk(i)=thp(i,1)
  qpk(i)=qp(i,1)
  zk(i) = z_theta(i,1)
END DO
DO i = 1, npnts
  thekp1(i)=th(i,2)
  qekp1(i)=q(i,2)
  qsekp1(i)=qse(i,2)
  bwk(i)=bwater(i,2)
  bwkp1(i)=bwater(i,2)
  zkp12(i)=z_rho(i, 2)
  zkp1(i)=z_theta(i, 2)
END DO

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

DO k = 2,nlev-1

  !  Initialise adaptive diagnostics for this level

  DO i = 1,npnts
    rbuoy_p_here(i) =0.0
    the_here(i) = th(i,k)
    thp_here(i) = th(i,k)
    qe_here(i) = q(i,k)
    qp_here(i) = q(i,k)
    rbuoy_p_here_c2(i) =0.0
    the_here_c2(i) = 0.0
    thp_here_c2(i) = 0.0
    qe_here_c2(i) = 0.0
    qp_here_c2(i) = 0.0
  END DO
  !Initialise adaptive entrainment variables
  DO i = 1, npnts
    thek(i)=th(i,k)
    qek(i)=q(i,k)
    qsek(i)=qse(i,k)
    thekp1(i)=th(i,k+1)
    qekp1(i)=q(i,k+1)
    qsekp1(i)=qse(i,k+1)
    bwk(i)=bwater(i,k)
    bwkp1(i)=bwater(i,k+1)
   !Note that unlike p_layer_boundaries, where k indexing is offset
   !by one compared to the dynamics numbering, z retains the numbering
   !convention for dynamics variables i.e. for theta levels, k->k
   !and for rho levels k+1/2 -> k+1
   !check this with Rachel
    zk(i) = z_theta(i,k)
    zkp12(i)=z_rho(i, k+1)
    zkp1(i)=z_theta(i, k+1)
  END DO
  IF(k  ==  2) THEN      !set to environmental values
    DO i = 1,npnts
      thpk(i)=th(i,2)
      qpk(i)=q(i,2)
    END DO
  ELSE
    DO i = 1,npnts
      thpk(i)=thp(i,k)
      qpk(i)=qp(i,k)
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

! DEPENDS ON: layer_cn_4a
  CALL layer_cn_4a(k,npnts,npnts,nlev,                               &
                   exner_layer_boundaries,                           &
                   exner_layer_centres,                              &
                   p_layer_boundaries,                               &
                   p_layer_centres,                                  &
                   pstar,pk,pkp1,delpk,delpkp1,                      &
                   delpkp12,ekp14,ekp34,amdetk,exk,exkp1,            &
                   delexkp1,delp_uv_k,delp_uv_kp1,recip_pstar,       &
                   rhum,l_shallow,ntml,ntpar                         &
                   , cumulus, mdet_md_on, md_ent_on                  &
                   , bconv                                           &
                   ,thek,qek, qsek, thekp1,qekp1,qsekp1              &
                   ,thpk, qpk                                        &
                   ,bwk,bwkp1,ekm14                                  &
                   ,zk, zkp12, zkp1                                  &
                   )

! Set ekm14 for next pass through loop
  DO i = 1, npnts
    ekm14(i) = ekp14(i)
  END DO

! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)

  IF (k == 2) THEN
! DEPENDS ON: dqs_dth_4a
    CALL dqs_dth_4a(dqsthk,k,th(1,k),qse(1,k),exk,npnts)
  ELSE
    DO i = 1,npnts
      dqsthk(i) = dqsthkp1(i)
    END DO
  END IF

! DEPENDS ON: dqs_dth
  CALL dqs_dth_4a(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,npnts)


! Set other grid dependent constants

  DO i = 1,npnts

    ! Maximum initial convective mass flux

    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

    ! Minimum buoyancy for convection to start from level k

    eminds(i) = mparb * delpkp12(i) * recip_pstar(i)

  END DO

! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
! qclp and qcfp zero at this point but will add an initial
! value at cloud base

  DO i = 1,npnts
    IF ( .NOT. bconv(i)) THEN
      expi(i)  = exk(i)
      xpk(i,k)  = 0.0
      qclp(i,k) = 0.0
      qcfp(i,k) = 0.0
      flx(i,k) = 0.0
      bgmk(i)  = .FALSE.
      depth(i) = 0.0
      thpi(i)  = th(i,k) + thpixs_v(i,k)
      thp(i,k)  = thpi(i)
      qpi(i)   = q(i,k) + qpixs_v(i,k)
      qp(i,k)   = qpi(i)
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


! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k > ntml+1 and the stability is low enough at k+1
!  (At levels above nbl this is the same as midtrig(i)=k)

  DO i = 1,npnts
    bcposs(i) = bconv(i) .OR.                                     &
                (p_layer_centres(i,k) > mid_cnv_pmin .AND.        &
                k > ntml(i) + 1 .AND.                             &
                (( th(i,k) * (1.0 + c_virtual * q(i,k)) -         &
                th(i,k+1) * (1.0 + c_virtual * q(i,k+1)) +        &
                delthst + MAX(0.0,(q(i,k)-qse(i,k+1))) *          &
                (lc/(cp * exkp1(i)))) > 0.0)) .OR.                &
                (k == ntml(i) .AND. ntml(i) == nbl-1)

  END DO  ! npnts

! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)

  ncposs = 0
  DO i = 1,npnts
    IF (bcposs(i)) THEN
      ncposs = ncposs + 1
      index1(ncposs) = i
    END IF
  END DO

! Compress points where convection may occur

  IF (ncposs  >   0) THEN
    DO i = 1,ncposs
      thek_c(i)   = th(index1(i),k)
      thekp1_c(i) = th(index1(i),k+1)
      qek_c(i)    = q(index1(i),k)
      qekp1_c(i)  = q(index1(i),k+1)
      qsekp1_c(i) = qse(index1(i),k+1)
      thpk_c(i)   = thp(index1(i),k)
      qpk_c(i)    = qp(index1(i),k)
      xpk_c(i)    = xpk(index1(i),k)
      bwk_c(i)    = bwater(index1(i),k)
      bwkp1_c(i)  = bwater(index1(i),k+1)
      pk_c(i)     = pk(index1(i))
      pkp1_c(i)   = pkp1(index1(i))
      ekp14_c(i)  = ekp14(index1(i))
      ekp34_c(i)  = ekp34(index1(i))
      exkp1_c(i)  = exkp1(index1(i))
      bgmkp1_c(i) = bgmk(index1(i)) ! bgmk into lift_par,
      eminds_c(i) = eminds(index1(i))
    END DO
    IF (l_q_interact .OR. cnv_wat_load_opt == 1) THEN
      DO i = 1,ncposs
        ! PC2 variables or needed for water loading
        qclekp1_c(i)  = qcl(index1(i),k+1)
        qcfekp1_c(i)  = qcf(index1(i),k+1)
      END DO
    END IF
    IF (l_q_interact) THEN
      DO i = 1,ncposs
        !   PC2 variables
        qclek_c(i)    = qcl(index1(i),k)
        qcfek_c(i)    = qcf(index1(i),k)
        qclpk_c(i)    = qclp(index1(i),k)
        qcfpk_c(i)    = qcfp(index1(i),k)
      END DO
    END IF
    IF (l_mom_gk) THEN
      DO i = 1,ncposs
        uek_c(i)    = u(index1(i),k)
        uekp1_c(i)  = u(index1(i),k+1)
        vek_c(i)    = v(index1(i),k)
        vekp1_c(i)  = v(index1(i),k+1)
        upk_c(i)    = up(index1(i),k)
        vpk_c(i)    = vp(index1(i),k)
       END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,ncposs
          traek_c(i,ktra) = tracer(index1(i),k,ktra)
          traekp1_c(i,ktra) = tracer(index1(i),k+1,ktra)
          trapk_c(i,ktra) = trap(index1(i),k,ktra)
        END DO
      END DO
    END IF
  END IF  ! ncposs>0

!-----------------------------------------------------------------------
! 3.2  Lift parcel from layer k to layer k+1
!-----------------------------------------------------------------------

! DEPENDS ON: lift_par_4a
  CALL lift_par_4a(ncposs,npnts,thpkp1_c,qpkp1_c,xsqkp1_c,        &
           bgmkp1_c,bwkp1_c,bwk_c,thpk_c,qpk_c,xpk_c,thekp1_c,    &
           qekp1_c,thek_c,qek_c,qsekp1_c,                         &
           qclpkp1_c,qclpk_c,qclekp1_c,qclek_c,l_q_interact,      &
           qcfpkp1_c,qcfpk_c,qcfekp1_c,qcfek_c,                   &
           pk_c,pkp1_c,exkp1_c,ekp14_c,ekp34_c,l_mom_gk,upkp1_c,  &
           vpkp1_c,upk_c,vpk_c,uek_c,uekp1_c,vek_c,vekp1_c,       &
           l_tracer,ntra,trapkp1_c,trapk_c,traekp1_c,             &
           traek_c)


! Loop over points which may convect (ncposs)

! NEC compiler directive
!CDIR NODEP

  DO i = 1,ncposs

    ! Calculate buoyancy (virt. potential temp.) of parcel in layer k+1

    IF (cnv_wat_load_opt == 0) THEN
      rbuoy_c(i) = thpkp1_c(i) * (1.0 + c_virtual *qpkp1_c(i))      &
                 - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i))
    ELSE IF (cnv_wat_load_opt == 1) THEN
      ! Include water loading using. The best guess for the parcel's qcl+qcf 
      ! level k+1 is the total parcel condensate at level k.
      rbuoy_c(i) = thpkp1_c(i) * (1.0 + c_virtual *qpkp1_c(i)       &
                 - xpk_c(i))                                        &
                 - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i)       &
                 - qclekp1_c(i) - qcfekp1_c(i))
    END IF

! Allow parcel to convect at midtrig or above if
!1. it has not convected in the past (bconv=.F.,bterm=.F.) AND
!2. k >= midtrig AND
!3. the buoyancy of parcel at k+1 > min buoyancy required for
!   convection + XSBMIN
! Flag convecting points with logical array bconv

    IF ( .NOT. bconv(index1(i)) .AND. .NOT. bterm(index1(i)) .AND.&
     ( k  >   midtrig(index1(i)).OR.                              &
     (k == ntml(index1(i)).AND.ntml(index1(i)) == nbl-1)) .AND.   &
           rbuoy_c(i)  >   (eminds_c(i) + xsbmin))THEN
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
                    pstar(index1(i)) * ((rbuoy_c(i) - xsbmin)     &
                     / delpkp12(index1(i))))
      IF (flxk_c(i)  >   flxmax(index1(i))) THEN
        flxk_c(i) = flxmax(index1(i))
      END IF
      IF (flg_up_flx) THEN
        up_flux(index1(i),k)=flxk_c(i)
      END IF

! Write compressed mass flux back to full array

      flx(index1(i),k) = flxk_c(i)

! At ntml set mixing detrainment rate equal to zero
! Store diagnostics linked to initial convective mass flux for
! calculation of final closure.

      amdetk(index1(i))      = 0.0
      flx_init(index1(i))    = flxk_c(i)
      flxmax_init(index1(i)) = flxmax(index1(i))

      IF (l_q_interact) THEN

!           Initialize QCLP(*,K) and QCFP(*,K) at start level and
!           perform the Parcel Lift again for such points: needed
!           for PC2.   Duplicates code from
!           Convection Parcel Lifting Scheme, LIFPAR.
!           Use non-compressed arrays for variables not updated
!           in or since LIFT_PAR
        qclpk_c(i)  = qcl(index1(i),k) +                          &
                     (1.0 + ekp14(index1(i))) *                   &
                     (qcl(index1(i),k+1) - qcl(index1(i),k))
        qclpkp1_c(i) = ( qclpk_c(i) +                             &
               (ekp14(index1(i)) * qcl(index1(i),k)) +            &
               (ekp34(index1(i)) *  (1.0 + ekp14(index1(i)))      &
                 * qcl(index1(i),k+1)) )/                         &
          ( (1.0 + ekp14(index1(i))) * (1.0 + ekp34(index1(i))) )

        qcfpk_c(i) = qcf(index1(i),k) +                           &
                    (1.0 + ekp14(index1(i))) *                    &
                    (qcf(index1(i),k+1) - qcf(index1(i),k))
        qcfpkp1_c(i) = ( qcfpk_c(i) +                             &
               (ekp14(index1(i)) * qcf(index1(i),k)) +            &
               (ekp34(index1(i)) * (1.0 + ekp14(index1(i)))       &
                 * qcf(index1(i),k+1)) )/                         &
          ( (1.0 + ekp14(index1(i))) * (1.0 + ekp34(index1(i))) )

! Working through the algebra of the above equations, we see that
!
!   ekp14(index1(i))     = 0 => qclpkp1_c(i) = qcl(index1(i),k+1)
!   ekp34(index1(i))     = 0 => qclpkp1_c(i) = qcl(index1(i),k+1)
!   qcl  (index1(i),k)   = 0 => qclpkp1_c(i) = qcl(index1(i),k+1)
!   qcl  (index1(i),k+1) = 0 => qclpkp1_c(i) = 0
!
! and similarly for the qcf variables. In finite arithmetic, however, these
! relationships break down. To avoid problems later on, we enforce them
! directly.

        IF (ekp14(index1(i))   == 0.0 .OR. &
            ekp34(index1(i))   == 0.0 .OR. &
            qcl  (index1(i),k) == 0.0) THEN
          qclpkp1_c(i) = qcl(index1(i),k+1)
        END IF
        IF (qcl(index1(i),k+1) == 0.0) qclpkp1_c(i) = 0.0

        IF (ekp14(index1(i))   == 0.0 .OR. &
            ekp34(index1(i))   == 0.0 .OR. &
            qcf  (index1(i),k) == 0.0) THEN
          qcfpkp1_c(i) = qcf(index1(i),k+1)
        END IF
        IF (qcf(index1(i),k+1) == 0.0) qcfpkp1_c(i) = 0.0

      END IF ! l_q_interact

    END IF

  END DO  !ncposs loop

! L_q_interact_if0:
  IF (l_q_interact) THEN
! ----------------------------------------------------------------------
!     Follow-on calculation from QCLP(*,K) and QCFP(*,K) initialization.
!       Duplicates code from Convection Parcel Lifting Scheme, LIFPAR.
! ----------------------------------------------------------------------
! If test for convection different from dp/sh conv.
! Use bconv and k eq start_lev as above

    DO i=1, ncposs
      IF (bconv(index1(i)).AND.k == start_lev(index1(i))) THEN
! ----------------------------------------------------------------------
!       CURRENTLY MIXED PHASE PARCEL IS FORBIDDEN. MELT OR FREEZE THE
!       ENTRAINED LAYER CLOUD AND ADJUST PARCEL TEMPERATURE ACCORDINGLY.
! ----------------------------------------------------------------------
 
        IF (bwater(index1(i),k+1) .AND. qcfpkp1_c(i)  >   0.0) THEN
          qclpkp1_c(i) = qclpkp1_c(i) + qcfpkp1_c(i)
          thpkp1_c(i)  = ( thpkp1_c(i) - ( qcfpkp1_c(i) * lf /      &
                                   (cp * exkp1(index1(i))) ) )
          qcfpkp1_c(i) = 0.0
        ELSE IF (.NOT. bwater(index1(i),k+1) .AND.                  &
                 qclpkp1_c(i)  >   0.0) THEN
          qcfpkp1_c(i) = qclpkp1_c(i) + qcfpkp1_c(i)
          thpkp1_c(i)  = ( thpkp1_c(i) + ( qclpkp1_c(i) * lf /      &
                                   (cp * exkp1(index1(i))) ) )
          qclpkp1_c(i) = 0.0
        END IF
      END IF
    END DO
  END IF  ! L_q_interact_if0

! Calculate number of points which are convecting  (nconv)
! set compression indices (index2).

  nconv = 0
  DO i = 1,ncposs
    IF (bconv(index1(i))) THEN
      nconv = nconv + 1
      index2(nconv) = i
    END IF
  END DO

! Second compression to form arrays of length nconv to be passed
! into CONVEC2
! Input variables to CONVEC2

  IF (nconv  >   0) THEN
    DO i = 1,nconv
      thek_c2(i)   = thek_c(index2(i))
      thekp1_c2(i) = thekp1_c(index2(i))
      qek_c2(i)    = qek_c(index2(i))
      qekp1_c2(i)  = qekp1_c(index2(i))
      dqsthkp1_c2(i) = dqsthkp1(index1(index2(i)))
      qsekp1_c2(i)   = qsekp1_c(index2(i))
      pstar_c2(i)    = pstar(index1(index2(i)))
      thpkp1_c2(i)   = thpkp1_c(index2(i))
      qpkp1_c2(i)    = qpkp1_c(index2(i))
      xsqkp1_c2(i) = xsqkp1_c(index2(i))
      rbuoy_c2(i)  = rbuoy_c(index2(i))
      qsek_c2(i)   = qse(index1(index2(i)),k)
      dqsthk_c2(i) = dqsthk(index1(index2(i)))
      thpi_c2(i)   = thpi(index1(index2(i)))
      qpi_c2(i)    = qpi(index1(index2(i)))
      expi_c2(i)   = expi(index1(index2(i)))
      bconv_c2(i)  = bconv(index1(index2(i)))
      bwk_c2(i)    = bwk_c(index2(i))
      bwkp1_c2(i)  = bwkp1_c(index2(i))
      bgmkp1_c2(i) = bgmkp1_c(index2(i))
      bland_c2(i)  = bland(index1(index2(i)))
      blowst_c2(i) = binit(index1(index2(i)))
      l_shallow_c2(i) = .FALSE.
      l_mid_c2(i)     = .TRUE.
      cumulus_c2(i)   = .FALSE.
      ekp14_c2(i)  = ekp14_c(index2(i))
      ekp34_c2(i)  = ekp34_c(index2(i))
      amdetk_c2(i) = amdetk(index1(index2(i)))
      pk_c2(i)     = pk_c(index2(i))
      pkp1_c2(i)   = pkp1_c(index2(i))
      exk_c2(i)    = exk(index1(index2(i)))
      exkp1_c2(i)  = exkp1_c(index2(i))
      delexkp1_c2(i)    = delexkp1(index1(index2(i)))
      delpk_c2(i)       = delpk(index1(index2(i)))
      delpkp1_c2(i)     = delpkp1(index1(index2(i)))
      delp_uv_k_c2(i)   = delp_uv_k(index1(index2(i)))
      delp_uv_kp1_c2(i) = delp_uv_kp1(index1(index2(i)))
      t1_sd_c2(i)     = t1_sd(index1(index2(i)))
      q1_sd_c2(i)     = q1_sd(index1(index2(i)))
      xsbmin_v_c2(i)  = xsbmin_v(index1(index2(i)),k)
    END DO
    IF (l_q_interact) THEN
      DO i=1,nconv
        !           PC2 variables
        qclek_c2(i)   = qclek_c(index2(i))
        qclekp1_c2(i) = qclekp1_c(index2(i))
        qcfek_c2(i)   = qcfek_c(index2(i))
        qcfekp1_c2(i) = qcfekp1_c(index2(i))
        qclpk_c2(i)   = qclpk_c(index2(i))
        qcfpk_c2(i)   = qcfpk_c(index2(i))
        !           (PC2) Compress input cloud fields
        cflek_c2(i)   = cf_liquid(index1(index2(i)),k)
        cflekp1_c2(i) = cf_liquid(index1(index2(i)),k+1)
        cffek_c2(i)   = cf_frozen(index1(index2(i)),k)
        cffekp1_c2(i) = cf_frozen(index1(index2(i)),k+1)
        bcfek_c2(i)   = bulk_cf(index1(index2(i)),k)
        bcfekp1_c2(i) = bulk_cf(index1(index2(i)),k+1)
        !           Compress convective base indicator
        start_lev3c(i) = start_lev(index1(index2(i)))
      END DO
    END IF
    IF (l_mom_gk) THEN
      DO i=1,nconv
        uek_c2(i)    = uek_c(index2(i))
        uekp1_c2(i)  = uekp1_c(index2(i))
        vek_c2(i)    = vek_c(index2(i))
        vekp1_c2(i)  = vekp1_c(index2(i))
        upkp1_c2(i)  = upkp1_c(index2(i))
        vpkp1_c2(i)  = vpkp1_c(index2(i))
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          traek_c2(i,ktra)   = traek_c(index2(i),ktra)
          traekp1_c2(i,ktra) = traekp1_c(index2(i),ktra)
          trapkp1_c2(i,ktra) = trapkp1_c(index2(i),ktra)
        END DO
      END DO
    END IF

! Input/output variables to/from CONVEC2

    DO i = 1,nconv
      thpk_c2(i)   = thpk_c(index2(i))
      qpk_c2(i)    = qpk_c(index2(i))
      xpk_c2(i)    = xpk_c(index2(i))
      flxk_c2(i)   = flx(index1(index2(i)),k)
      bgmk_c2(i)   = bgmk(index1(index2(i)))
      bterm_c2(i)  = .FALSE.
      dthek_c2(i)  = dthbydt(index1(index2(i)),k)
      dqek_c2(i)   = dqbydt(index1(index2(i)),k)
      dthekp1_c2(i)  = dthbydt(index1(index2(i)),k+1)
      dqekp1_c2(i)   = dqbydt(index1(index2(i)),k+1)
      tcw_c2(i)     = tcw(index1(index2(i)))
      depth_c2(i)   = depth(index1(index2(i)))
      cclwp_c2(i)   = cclwp(index1(index2(i)))
      cape_c2(i)    = cape(index1(index2(i)))
      dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
      relh_c2(i)    = relh(index1(index2(i)))
      rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
      the_here_c2(i)=the_here(index1(index2(i)))
      thp_here_c2(i)=thp_here(index1(index2(i)))
      qe_here_c2(i)=qe_here(index1(index2(i)))
      qp_here_c2(i)=qp_here(index1(index2(i)))
      rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
      dptot_c2(i)   = dptot(index1(index2(i)))
      thpixs_v_c2(i) = thpixs_v(index1(index2(i)),k)
      qpixs_v_c2(i)  = qpixs_v(index1(index2(i)),k)
      iccb_c2(i)     = iccb(index1(index2(i)))
      icct_c2(i)     = icct(index1(index2(i)))
      cca_2d_c2(i)   = cca_2d(index1(index2(i)))
      lcca_c2(i)     = lcca(index1(index2(i)))
      lcbase_c2(i)   = lcbase(index1(index2(i)))
      lctop_c2(i)    = lctop(index1(index2(i)))
    END DO
    IF (l_q_interact) THEN
      DO i = 1,nconv
        !           PC2 variables
        qclpkp1_c2(i) = qclpkp1_c(index2(i))
        qcfpkp1_c2(i) = qcfpkp1_c(index2(i))
        !           Compress increment fields
        dqclek_c2(i)   = dqclbydt(index1(index2(i)),k)
        dqclekp1_c2(i) = dqclbydt(index1(index2(i)),k+1)
        dqcfek_c2(i)   = dqcfbydt(index1(index2(i)),k)
        dqcfekp1_c2(i) = dqcfbydt(index1(index2(i)),k+1)
        dcflek_c2(i)   = dcflbydt(index1(index2(i)),k)
        dcflekp1_c2(i) = dcflbydt(index1(index2(i)),k+1)
        dcffek_c2(i)   = dcffbydt(index1(index2(i)),k)
        dcffekp1_c2(i) = dcffbydt(index1(index2(i)),k+1)
        dbcfek_c2(i)   = dbcfbydt(index1(index2(i)),k)
        dbcfekp1_c2(i) = dbcfbydt(index1(index2(i)),k+1)
      END DO
    END IF
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        upk_c2(i)    = upk_c(index2(i))
        vpk_c2(i)    = vpk_c(index2(i))
        duek_c2(i)   = dubydt(index1(index2(i)),k)
        dvek_c2(i)   = dvbydt(index1(index2(i)),k)
        duekp1_c2(i)   = dubydt(index1(index2(i)),k+1)
        dvekp1_c2(i)   = dvbydt(index1(index2(i)),k+1)
        eflux_u_ud_c2(i) = eflux_u_ud(index1(index2(i)))
        eflux_v_ud_c2(i) = eflux_v_ud(index1(index2(i)))
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          trapk_c2(i,ktra) = trapk_c(index2(i),ktra)
          dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
          dtraekp1_c2(i,ktra) = dtrabydt(index1(index2(i)),k+1,ktra)
        END DO
      END DO
    END IF

!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2_4a5a
    CALL convec2_4a5a(nconv,npnts,nlev,ntra,k,md_on,md_sdet_on,md_new_termc,&
               start_lev3c,                                           &
               l_tracer,l_mom_gk,l_q_interact,l_calc_dxek,            &
               l_shallow_c2,l_mid_c2,cumulus_c2,                      &
               bwk_c2,bwkp1_c2,bgmkp1_c2,bland_c2,blowst_c2,          &
               timestep,                                              &
               thek_c2,thekp1_c2,qek_c2,qekp1_c2,qclek_c2,qclekp1_c2, &
               qcfek_c2,qcfekp1_c2,cflek_c2,cflekp1_c2,               &
               cffek_c2,cffekp1_c2,bcfek_c2,bcfekp1_c2,               &
               uek_c2,uekp1_c2,vek_c2,vekp1_c2,                       &
               traek_c2,traekp1_c2,trapkp1_c2,                        &
               qsekp1_c2,dqsthkp1_c2,pstar_c2,thpkp1_c2,qpkp1_c2,     &
               upkp1_c2,vpkp1_c2,xsqkp1_c2,rbuoy_c2,qsek_c2,dqsthk_c2,&
               thpi_c2,qpi_c2,expi_c2,ekp14_c2,ekp34_c2,amdetk_c2,    &
               pk_c2,pkp1_c2,exk_c2,exkp1_c2,delexkp1_c2,delpk_c2,    &
               delpkp1_c2,delp_uv_k_c2, delp_uv_kp1_c2,               &
               t1_sd_c2,q1_sd_c2,thpixs_v_c2,qpixs_v_c2,xsbmin_v_c2,  &
               rbuoy_p_old_c2,                                        &
             ! In/out
               bconv_c2,bgmk_c2,                                      &
               thpk_c2,qpk_c2,qclpk_c2,qcfpk_c2,qclpkp1_c2,qcfpkp1_c2,&
               upk_c2,vpk_c2,trapk_c2,xpk_c2,flxk_c2,                 &
               dthek_c2,dqek_c2,dqclek_c2,dqcfek_c2,dcflek_c2,        &
               dcffek_c2,dbcfek_c2,duek_c2,dvek_c2,dtraek_c2,         &
               tcw_c2,depth_c2,cclwp_c2,cape_c2,dcpbydt_c2,           &
               eflux_u_ud_c2,eflux_v_ud_c2,                           &
             ! Out
               iccb_c2,icct_c2,lcbase_c2,lctop_c2,                    &
               bterm_c2,                                              &
               prekp1_c2, dthekp1_c2,dqekp1_c2,                       &
               dqclekp1_c2,dqcfekp1_c2,dcflekp1_c2,dcffekp1_c2,       &
               dbcfekp1_c2,duekp1_c2,dvekp1_c2,dtraekp1_c2,           &
               cca_2d_c2,ccw_c2,lcca_c2,deltak_c2,flxkp12_c2,         &
               max_cfl_c2,relh_c2,dptot_c2,rbuoy_p_here_c2,           &
               the_here_c2,thp_here_c2,qe_here_c2,qp_here_c2)

! Calculate fractional entrainment rate for level k.
! If convection has terminated (bterm=.T.) then set
! fractional entrainment rate for k+1 to zero.
   
    DO i = 1,nconv
      entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *     &
                 (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i) &
                * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
      IF (bterm_c2(i)) THEN
        entrain_up(index1(index2(i)),k+1) = 0.0
      END IF
    END DO

! Calculate fractional detrainment rate for level k
! (and k+1 if bterm=.T.)

    IF (flg_detr_up) THEN
      DO i = 1,nconv
        detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)            &
                    + deltak_c2(i) * (1.0 - amdetk_c2(i)))        &
                  * flx(index1(index2(i)),k)
        IF (bterm_c2(i)) THEN
          detrain_up(index1(index2(i)),k+1) =                       &
                 -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
        END IF
      END DO
    END IF
  END IF !if nconv>0

! Write CONVEC2 compressed output arrays back to full fields

  DO i = 1,npnts
    thp(i,k+1)    = 0.0
    qp(i,k+1)     = 0.0
    xpk(i,k+1)    = 0.0
    flx(i,k+1)    = 0.0
    depth(i)      = 0.0
    precip(i,k+1) = 0.0
    qclp(i,k+1)   = 0.0
    qcfp(i,k+1)   = 0.0
    bgmk(i)       = .FALSE.
    bterm(i)      = .FALSE.
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO i = 1,npnts
        trap(i,k+1,ktra) = 0.0
      END DO
    END DO
  END IF

  IF (l_mom_gk) THEN
    DO i = 1,npnts
      up(i,k+1) = 0.0
      vp(i,k+1) = 0.0
    END DO
  END IF

  IF (nconv  >   0) THEN
    DO i = 1,nconv
      thp(index1(index2(i)),k+1) = thpkp1_c2(i)
      qp(index1(index2(i)),k+1)  = qpkp1_c2(i)
      xpk(index1(index2(i)),k+1) = xpk_c2(i)
      flx(index1(index2(i)),k+1) = flxk_c2(i)
      depth(index1(index2(i)))   = depth_c2(i)
      precip(index1(index2(i)),k+1) = prekp1_c2(i)
      bgmk(index1(index2(i)))    = bgmk_c2(i)
      bterm(index1(index2(i)))   = bterm_c2(i)
      dthbydt(index1(index2(i)),k) = dthek_c2(i)
      dqbydt(index1(index2(i)),k)  = dqek_c2(i)
      dthbydt(index1(index2(i)),k+1) = dthekp1_c2(i)
      dqbydt(index1(index2(i)),k+1)  = dqekp1_c2(i)
      cca_2d(index1(index2(i)))  = cca_2d_c2(i)
      tcw(index1(index2(i)))     = tcw_c2(i)
      iccb(index1(index2(i)))    = iccb_c2(i)
      icct(index1(index2(i)))    = icct_c2(i)
      cclwp(index1(index2(i)))   = cclwp_c2(i)
      lcca(index1(index2(i)))    = lcca_c2(i)
      lcbase(index1(index2(i)))  = lcbase_c2(i)
      lctop(index1(index2(i)))   = lctop_c2(i)
      ccw(index1(index2(i)),k+1) = ccw_c2(i)
      cape(index1(index2(i)))    = cape_c2(i)
      dcpbydt(index1(index2(i))) = dcpbydt_c2(i)
      rbuoy_p_here(index1(index2(i))) = rbuoy_p_here_c2(i)
      the_here(index1(index2(i))) = the_here_c2(i)
      thp_here(index1(index2(i))) = thp_here_c2(i)
      qe_here(index1(index2(i)))  = qe_here_c2(i)
      qp_here(index1(index2(i)))  = qp_here_c2(i)
    END DO
    DO i = 1,nconv
      max_cfl(index1(index2(i))) =                                  &
                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
      relh(index1(index2(i)))  = relh_c2(i)
      dptot(index1(index2(i))) = dptot_c2(i)
    END DO
    IF (l_q_interact) THEN
      DO i = 1,nconv
        qclp(index1(index2(i)),k+1)   = qclpkp1_c2(i)
        qcfp(index1(index2(i)),k+1)   = qcfpkp1_c2(i)
        dqclbydt(index1(index2(i)),k)   = dqclek_c2(i)
        dqcfbydt(index1(index2(i)),k)   = dqcfek_c2(i)
        dcflbydt(index1(index2(i)),k)   = dcflek_c2(i)
        dcffbydt(index1(index2(i)),k)   = dcffek_c2(i)
        dbcfbydt(index1(index2(i)),k)   = dbcfek_c2(i)
        dqclbydt(index1(index2(i)),k+1) = dqclekp1_c2(i)
        dqcfbydt(index1(index2(i)),k+1) = dqcfekp1_c2(i)
        dcflbydt(index1(index2(i)),k+1) = dcflekp1_c2(i)
        dcffbydt(index1(index2(i)),k+1) = dcffekp1_c2(i)
        dbcfbydt(index1(index2(i)),k+1) = dbcfekp1_c2(i)
      END DO
    END IF
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        up(index1(index2(i)),k+1) = upk_c2(i)
        vp(index1(index2(i)),k+1) = vpk_c2(i)
        dubydt(index1(index2(i)),k) = duek_c2(i)
        dvbydt(index1(index2(i)),k) = dvek_c2(i)
        dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
        dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
        eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
        eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
        flxkp12(k,index1(index2(i)))  = flxkp12_c2(i)
      END DO
    END IF
    IF  (l_tracer) THEN
      DO i = 1,nconv
        DO ktra = 1,ntra
          trap(index1(index2(i)),k+1,ktra)     = trapk_c2(i,ktra)
          dtrabydt(index1(index2(i)),k,ktra)   = dtraek_c2(i,ktra)
          dtrabydt(index1(index2(i)),k+1,ktra) = dtraekp1_c2(i,ktra)
        END DO
      END DO
    END IF
    IF (flg_up_flx) THEN
      DO i = 1,nconv
        up_flux(index1(index2(i)),k+1) = flxk_c2(i)
      END DO
    END IF
    IF (flg_up_flx_half) THEN
      DO i = 1,nconv
        up_flux_half(index1(index2(i)),k+1) =  flxkp12_c2(i)
      END DO
    END IF
  END IF     ! nconv >0

!   Write adaptive diagnostics for this level to full array for output


!Write rbuoy for previous level for next pass through loop
  DO i = 1,npnts
    rbuoy_p_old(i) = rbuoy_p_here(i)
  END DO

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
          !  than 150hPa
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
              cape_ts_new(i) = cape_ts_w * w_cape_limit/               &
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
              cape_ts_new(i) = cape_timescale * w_cape_limit/w_max(i)
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
                cape_ts_new(i) = cape_timescale *                       &
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

    DO j = 1,nterm
      i = index_nterm(j)
      max_cfl(i) = max_cfl(i) * timestep

      IF (max_cfl(i)  >   cfl_limit) THEN
        flx_init_new(i) = flx_init_new(i) * cfl_limit           &
                              / max_cfl(i)
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

! Carry out cape and cfl scaling

    DO kt = 2, k+1
! NEC compiler directive
!CDIR NODEP
      DO j = 1,nterm
      i = index_nterm(j)
        IF (kt  >=  start_lev(i) .AND. flx_init_new(i)            &
                                            >   0.0) THEN

          IF (kt  ==  det_lev(i)) THEN
            dthbydt(i,kt)  = ((dthbydt(i,kt) - dthef(i))          &
                                       * scale_f(i)) + dthef(i)
            dqbydt(i,kt)   = ((dqbydt(i,kt) - dqf(i))             &
                                       * scale_f(i)) + dqf(i)

            IF (l_q_interact) THEN  ! PC2
              dqclbydt(i,kt)   = ((dqclbydt(i,kt) - dqclf(i))     &
                                       * scale_f(i)) + dqclf(i)
              dqcfbydt(i,kt)   = ((dqcfbydt(i,kt) - dqcff(i))     &
                                       * scale_f(i)) + dqcff(i)
              dcflbydt(i,kt)   = ((dcflbydt(i,kt) - dcflf(i))     &
                                       * scale_f(i)) + dcflf(i)
              dcffbydt(i,kt)   = ((dcffbydt(i,kt) - dcfff(i))     &
                                       * scale_f(i)) + dcfff(i)
              dbcfbydt(i,kt)   = ((dbcfbydt(i,kt) - dbcff(i))     &
                                       * scale_f(i)) + dbcff(i)
            END IF

            IF (l_mom_gk) THEN
              dubydt(i,kt) = ((dubydt(i,kt) - duef(i))            &
                                       * scale_f(i)) + duef(i)
              dvbydt(i,kt) = ((dvbydt(i,kt) - dvef(i))            &
                                       * scale_f(i)) + dvef(i)

            END IF

            IF (l_tracer) THEN
              DO ktra = 1,ntra
                dtrabydt(i,kt,ktra) = ((dtrabydt(i,kt,ktra) -     &
                                   dtraef(i,ktra)) * scale_f(i))  &
                                        + dtraef(i,ktra)
              END DO
            END IF
          ELSE
            dthbydt(i,kt) = dthbydt(i,kt) * scale_f(i)
            dqbydt(i,kt)  = dqbydt(i,kt) * scale_f(i)

            IF (l_q_interact) THEN  ! PC2
              dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
              dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
              dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
              dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
              dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
            END IF

            IF (l_mom_gk) THEN
              dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
              dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
              IF (kt <  k+1) THEN
                flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
              END IF
            END IF
            IF (l_tracer) THEN
              DO ktra = 1,ntra
                dtrabydt(i,kt,ktra) =dtrabydt(i,kt,ktra)*scale_f(i)
              END DO
            END IF

          END IF !kt=det_lev
          flx(i,kt)    = flx(i,kt) * scale_f(i)
          precip(i,kt) = precip(i,kt) * scale_f(i)

          IF (flg_up_flx) THEN
            up_flux(i,kt) = flx(i,kt)
          END IF
          IF (kt <  k+1) THEN
            IF (flg_up_flx_half) THEN
              up_flux_half(i,kt+1) = flxkp12(kt,i)
            END IF
          END IF

          entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)

          IF (flg_detr_up) THEN
            detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
          END IF

        END IF ! kt>start_lev and flx_init_new >0
      END DO  ! j loop
    END DO  ! kt loop

! Set cape scaling parameters for next level

! NEC compiler directive
!CDIR NODEP
    DO j = 1,nterm
      i = index_nterm(j)
      dthef(i) = dthbydt(i,k+1)
      dqf(i)   = dqbydt(i,k+1)
      IF (l_q_interact) THEN  ! PC2
        dqclf(i)   = dqclbydt(i,k+1)
        dqcff(i)   = dqcfbydt(i,k+1)
        dcflf(i)   = dcflbydt(i,k+1)
        dcfff(i)   = dcffbydt(i,k+1)
        dbcff(i)   = dbcfbydt(i,k+1)
      END IF
      det_lev(i)= k+1
      IF (l_mom_gk) THEN
        duef(i)   = dubydt(i,k+1)
        dvef(i)   = dvbydt(i,k+1)
      END IF
    END DO  ! nterm loop
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO j = 1,nterm
          i = index_nterm(j)
          dtraef(i,ktra) = dtrabydt(i,k+1,ktra)
        END DO  ! nterm loop
      END DO
    END IF

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
        kterm_mid(i) = k
      END DO

    ELSE        ! original down draughts

      npossdd = 0
      nnodd = 0
      DO j = 1,nterm
        i = index_nterm(j)
        tempnum = 0.0
        IF (iccb(i)  >   0) THEN
          deltap_cld = p_layer_centres(i,iccb(i))                     &
                               - p_layer_centres(i,k)
          DO kt = iccb(i), k+1
            tempnum = tempnum + precip(i,kt)
          END DO
        ELSE
          deltap_cld = 0.0
        END IF

! Downdraughts possible if pressure thickness of convective
! cloud (deltap_cld) is greater than 15000m, the point is saturated
! and the precip. in the layer is greater than a threshold
! value (1E-12).

        IF (deltap_cld  >   15000.0 .AND. bgmk(i) .AND.               &
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

! DEPENDS ON: dd_call_4a5a
        CALL dd_call_4a5a(npnts, npossdd, k, nlev, trlev, ntra      &
            ,iccb, icct, index_possdd                               &
            ,l_tracer                                               &
            ,bwater(1,2)                                            &
            ,exner_layer_centres, exner_layer_boundaries            &
            ,p_layer_centres, p_layer_boundaries, pstar             &
            ,recip_pstar, timestep , cca_2d                         &
            ,thp(1,1), qp(1,1), th(1,1), q(1,1), qse                &
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

! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).

  ! Prevent model from wiping out single layer cloud. Problems with the
  ! radiation scheme may well have occured if iccb == icct, however the
  ! radiation scheme expects the cloud base and top boundaries not
  ! layers as specified by iccb and icct in the convection scheme so
  ! that even with single layer clouds iccb and icct that are passed to
  ! the radiation scheme should not be equal to each other.

    IF (.NOT. l_ccrad) THEN
      DO j=1, nterm
        i = index_nterm(j)

        IF (iccb(i) == icct(i)) THEN
          iccb   (i) = 0
          icct   (i) = 0
          cca_2d (i) = 0.0
          tcw    (i) = 0.0
          cclwp  (i) = 0.0
        END IF

        IF (lcbase(i) == lctop(i)) THEN
          lcbase (i) = 0
          lctop  (i) = 0
          lcca   (i) = 0.0
        END IF

      END DO
    END IF      ! l_ccrad

! If convection has terminated write cape to diagnostic output
! variable (cape_out). Note if more than one lot of mid-level
! convection or if mid above deep or shallow then stores highest
! convection results only.
! Set kterm array which holds the level index for termination
! of convection.
! Zero integrals as convection terminated. Note may convect again at
! same locations at higher levels in the atmosphere.

! NEC compiler directive
!CDIR NODEP
    DO j = 1,nterm
      i=index_nterm(j)
      cape_out(i) = cape(i)
      dcpbydt(i) = 0.0
      cape(i) = 0.0
      kterm(i) = k
      kterm_mid(i) = k
      bconv(i) = .FALSE.
      bterm(i) = .FALSE.     ! reset bterm to false
      start_lev(i) = 0
      relh(i)     = 0.0      ! needed for RH CAPE
      dptot(i)    = 0.0      !
    END DO

  END IF        ! nterm >0

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

END DO

! Has mid-level convection gone to the top of the model ?

DO i=1,npnts
  IF (kterm_mid(i) == nlev-1 .AND. l_mid_all(i)) THEN
    IF (printstatus >= prstatus_normal) THEN
      write(6,'(a33,i8)') 'WARNING: mid_conv reached top at ',i 
    END IF
    error_point = i  
    IF (lhook) CALL dr_hook('MID_CONV_4A',zhook_out,zhook_handle)
    RETURN        
  END IF
END DO

IF (l_ccrad) THEN

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

END IF      ! l_ccrad

!-----------------------------------------------------------------------
! Emanuel down draughts - one call for all points with mid-level conv
!-----------------------------------------------------------------------

IF (l_eman_dd) THEN

! how many points have mid level convection?
! What is the maximum termination level ?

   kterm_max = 2
   nmid = 0
   DO i=1,npnts
     IF (l_mid_all(i)) THEN
        nmid = nmid +1
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
,                  timestep                                       &
,                  th, q, qse, tracer ,precip                     &
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

  nmid = 0
  DO i=1,npnts
    IF (l_mid_all(i)) THEN
      nmid = nmid +1
      index1(nmid) = i
    END IF
  END DO
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

!-----------------------------------------------------------------------------
! 5.0 Calculate CCA
!-----------------------------------------------------------------------------

IF (l_ccrad) THEN

  n_mdcld      = 0    ! Initialise no of points with mid-cld
  n_mdcld_mult = 0    ! Initialise no of points with multiple mid-cld

  DO i=1, npnts
    IF ((l_mid_all(i)) .AND. (lcbase(i) /= 0)) THEN

      ! Mid-level cloud present
      ! Create index of gridpoints with single/multiple mid-level cloud
      IF (lcbase(i) /= iccb(i)) THEN

        n_mdcld_mult       = n_mdcld_mult + 1
        dum2(n_mdcld_mult) = i

      ELSE IF (lcbase(i) == iccb(i)) THEN

        n_mdcld       = n_mdcld + 1
        dum1(n_mdcld) = i
      END IF

    END IF
  END DO


  !-----------------------------------------------------------------
  ! Calculate CCA_2D
  !-----------------------------------------------------------------
  SELECT CASE (cca2d_md_opt)
    CASE(srf_precip)
      DO i=1, npnts
        IF ((l_mid_all(i)) .AND. (lcbase(i) /= 0)) THEN
         !-------------------------------------------------------------
         ! CCA_2D based on surface precipitation rate from mid cloud
         !-------------------------------------------------------------
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
            ! shallow or deep have occured
            !
            ! NOTE: Downdraughts & Evaporation still being fed cca_2d
            !       derived from TCW, This issue may require further
            !       investigation.
            lcca(i) = cca_2d(i)
          END IF

        END IF      ! mid-level cloud present
      END DO      ! npnts

    CASE(total_condensed_water)
      ! cca_2d_md left unchanged from code, which is based on
      ! TCW (Total Condensed Water) (This is a rate)

  END SELECT      ! cca2d_md_opt


  IF (.NOT. l_dcpl_cld4pc2) THEN
    !-------------------------------------------------------------------
    ! Apply Mid-leveld CCA Tuning factor
    !-------------------------------------------------------------------
    ! Because of decoupling move up a level so scaling is only applied
    ! to section 0 variables, section 5 variables continue to be
    ! scaled here to maintain bit-comparison, would ideally be moved
    ! up to Glue_Conv in future
    DO i=1, npnts
      ! Apply Mid-Level CCA Tuning factor
      !----------------------------------------
      cca_2d(i) = cca_md_knob*cca_2d(i)

      ! Make sure cca_2d_md is within limits
      !----------------------------------------
      cca_2d(i) = MAX(0.0,    cca_2d(i))
      cca_2d(i) = MIN(1.0e+0, cca_2d(i))
    END DO      ! i (npnts)
  END IF

  !---------------------------------------------------------------------
  ! 5.1 Apply CCA_3D
  !---------------------------------------------------------------------
  ! A value of cca_2d has been calculated, though we need to know which
  ! levels to apply it to. So we identify which layers require CCA based
  ! on layers (above shallow/deep, if they have occurred) with non-zero
  ! ccw.
  !---------------------------------------------------------------------

  !=====================================================================
  ! Single mid-level cloud bank gridpoints
  !=====================================================================

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

    !-------------------------------------------------------------------
    ! End Resizing compressed arrays for single Mid-level events. There
    ! is only one mid-level cloud bank if a mid-event has occurred, so
    ! does not require checking to see if there are multiple mid-level
    ! cloud banks.
    !-------------------------------------------------------------------


    IF (l_anvil) THEN
! DEPENDS ON: calc_3d_cca
      CALL calc_3d_cca                                                  &
        ( n_mdcld, n_mdcld, nlev, n_cca_lev, nbl, iccb_md_c, icct_md_c  &
        , p_lyr_bnds_md_c, freeze_lev_md_c, cca_2d_md_c, cca_md_c       &
        , z_theta_md_c, z_rho_md_c )

    ELSE
      !-----------------------------------------------------------------
      ! Do not use anvil scheme and apply cca_2d_md_c to all cloud
      ! levels between mid-level base and top
      !-----------------------------------------------------------------
      DO k=1, n_cca_lev
        DO i=1, n_mdcld
          IF ( k >= iccb_md_c(i) .AND. &
               k <= icct_md_c(i) ) THEN
            cca_md_c(i,k) = cca_2d_md_c(i)
          ELSE
            cca_md_c(i,k) = 0.0
          END IF
        END DO      ! i (n_mdcld)
      END DO      ! k (nlev)

    END IF     ! l_anvil

    !-------------------------------------------------------------------
    ! Merge cca_md/ccw_md to full cca array and scale ccw_md_c by
    ! ccw_md_knob
    !-------------------------------------------------------------------
    DO k=1, n_cca_lev
      DO i=1, n_mdcld
        cca(mdcldi(i),k) = cca(mdcldi(i),k) + cca_md_c(i,k)
      END DO      ! i (n_mdcld)
    END DO      ! k (nlev)


    IF (.NOT. l_dcpl_cld4pc2) THEN
      DO k=1, nlev
        DO i=1, n_mdcld
          ccw(mdcldi(i),k) = ccw_md_c(i,k)*ccw_md_knob
        END DO      ! i (n_mdcld)
      END DO      ! k (nlev)
    END IF

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

  !=====================================================================
  ! Multiple mid-level cloud bank gridpoints
  !=====================================================================

  IF (n_mdcld_mult > 0) THEN
    ! Resize index array for gridpoints with multiple Mid-level cloud
    ! banks and apply indices

    ALLOCATE(mdcldi_mult(n_mdcld_mult))

    !-------------------------------------------------------------------
    ! Allocate arrays with multiple-mid level cloud gridpoints
    !-------------------------------------------------------------------
    ALLOCATE (iccb_md_mult_c       (n_mdcld_mult       ))
    ALLOCATE (icct_md_mult_c       (n_mdcld_mult       ))
    ALLOCATE (freeze_lev_md_mult_c (n_mdcld_mult       ))
    ALLOCATE (cca_2d_md_mult_c     (n_mdcld_mult       ))

    ! Lets allocate them with levels first since it will then be continuous in
    ! memory later on when we call calc_3d_cca.
    ALLOCATE (cca_md_mult_c        (  n_cca_lev,n_mdcld_mult))
    ALLOCATE (ccw_md_mult_c        (  nlev     ,n_mdcld_mult))
    ALLOCATE (z_theta_md_mult_c    (  nlev     ,n_mdcld_mult))
    ALLOCATE (z_rho_md_mult_c      (  nlev     ,n_mdcld_mult))
    ALLOCATE (p_lyr_bnds_md_mult_c (0:nlev     ,n_mdcld_mult))

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
        p_lyr_bnds_md_mult_c(k,i) =                                     &
                          p_layer_boundaries(mdcldi_mult(i),k)
      END DO

      DO k=1, n_cca_lev
        cca_md_mult_c(k,i) = 0.0
      END DO

    ! For multiple mid-level cloud banks, increment up model levels
    ! through compressed ccw_md_c array. Enter CAL3DCCA on locating
    ! cloud/bases of multiple clouds.

      DO k=2, nlev

        !---------------------------------------------------------------
        ! Check for cloud base
        !---------------------------------------------------------------
        IF ((ccw_md_mult_c(k,i) >  0.0) .AND.                           &
            (iccb_md_mult_c(i)  == 0)   .AND.                           &
            (icct_md_mult_c(i)  == 0)) THEN
          iccb_md_mult_c(i) = k
        END IF

        !---------------------------------------------------------------
        ! Check for cloud top
        !---------------------------------------------------------------
        IF ((ccw_md_mult_c(k,i)   <= 0.0) .AND.                         &
            (ccw_md_mult_c(k-1,i) >  0.0) .AND.                         &
            (iccb_md_mult_c(i)    /= 0)   .AND.                         &
            (icct_md_mult_c(i)    == 0)) THEN
          icct_md_mult_c(i) = k-1
        END IF

        !---------------------------------------------------------------
        ! Check for for anvil if both a cloud base and top found
        !---------------------------------------------------------------
        IF (iccb_md_mult_c(i) /= 0 .AND.                                &
            icct_md_mult_c(i) /= 0) THEN

          ! Apply CCA to vertically continuous cloud
          IF (l_anvil) THEN
            ! Since we are performing this on each individual point we can
            ! pass down the array with level information in the first rank
            ! instead of the second rank.  This reduces the amount of temporary
            ! data which is required otherwise it has to step through each
            ! array by number of points to get the all level information for it.
! DEPENDS ON: calc_3d_cca
            CALL calc_3d_cca                                            &
              ( 1, 1, nlev, n_cca_lev,nbl, iccb_md_mult_c(i)            &
              , icct_md_mult_c(i), p_lyr_bnds_md_mult_c(:,i)            &
              , freeze_lev_md_mult_c(i), cca_2d_md_mult_c(i)            &
              , cca_md_mult_c(:,i), z_theta_md_mult_c(:,i)              &
              , z_rho_md_mult_c(:,i) )

          ELSE

            ! Copy cca_2d_md_c to all levels from
            ! cloud base to cloud top of located cloud
            DO j=iccb_md_mult_c(i), icct_md_mult_c(i)
              cca_md_mult_c(j,i) = cca_2d_md_mult_c(i)
            END DO      ! j

          END IF      ! l_anvil

          iccb_md_mult_c(i) = 0
          icct_md_mult_c(i) = 0

        END IF      ! Test on iccb_md_c(i) and icct_md_c(i)
    
        IF (.NOT. l_dcpl_cld4pc2) THEN
          ccw(mdcldi_mult(i),k) = ccw_md_mult_c(k,i)*ccw_md_knob
        END IF

      END DO      ! k (nlev)


    !-------------------------------------------------------------------
    ! Merge cca_md/ccw_md to cca full array and scale ccw_md_c by
    ! ccw_md_knob
    !-------------------------------------------------------------------
      DO k=1, n_cca_lev
        cca(mdcldi_mult(i),k) = cca(mdcldi_mult(i),k)                   &
                              + cca_md_mult_c(k,i)
      END DO

    END DO



    !-------------------------------------------------------------------
    ! Deallocate arrays
    !-------------------------------------------------------------------
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

END IF      ! l_ccrad

!-----------------------------------------------------------------------
! 6.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('MID_CONV_4A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE mid_conv_4a
