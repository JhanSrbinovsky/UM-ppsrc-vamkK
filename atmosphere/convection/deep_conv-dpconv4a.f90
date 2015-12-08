! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Deep convection scheme
!
SUBROUTINE deep_conv_4a(nbl,nlev,ntra,n_cca_lev,n_dp,trlev,       &
                       bland, delthvu, exner_layer_centres,       &
                       exner_layer_boundaries,                    &
                       l_calc_dxek, l_q_interact,                 &
                       l_tracer,ntml,ntpar,pstar,                 &
                       p_layer_centres, p_layer_boundaries,       &
                       z_theta, z_rho,r_theta,r_rho,              &
                       rho_theta,rho,r2rho_th,r2rho,              &
                       dr_across_th,dr_across_rh,                 &
                       q,q1_sd,t1_sd,th,                          &
                       timestep, u, v, uw0, vw0,                  &
                       w_max, wstar, zlcl_uv,                     &
                       freeze_lev,recip_pstar,qse,                &
                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
                       qcl,tracer,cape_out,cclwp,ccw,cca,         &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                       dqclbydt,dthbydt,                          &
                       dubydt,dvbydt,dtrabydt,                    &
                       detrain_up,detrain_dwn,                    &
                       entrain_up,entrain_dwn,                    &
                       iccb,icct,lcca,                            &
                       lcbase,lctop,rain,snow,rain_3d,snow_3d,    &
                       up_flux, up_flux_half,                     &
                       dwn_flux,uw_deep,vw_deep,kterm,tcw,cca_2d, &
                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out, &
                       ind_cape_reduced,cape_ts_used,cfl_limited, &
                       ind_deep,error_point)


! Purpose:
!   Deep convection scheme - works on points diagnosed as deep in
!   subroutine CONV_DIAG.

!   Called by GLUE_CONV.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP 3 v8.3 programming standards

USE atmos_constants_mod, ONLY:                                    &
    r, cp, kappa, pref, repsilon, c_virtual

USE cv_run_mod, ONLY:                                             &
    l_mom, l4a_kterm, l_eman_dd, cape_opt, cape_min,              &
    w_cape_limit, cape_timescale, deep_cmt_opt, bl_cnv_mix,       &
    limit_pert_opt, cca2d_dp_opt, cca_dp_knob, ccw_dp_knob,       &
    l_anvil, l_dcpl_cld4pc2, cnv_wat_load_opt, l_ccrad

USE cv_dependent_switch_mod, ONLY:                                &
    dp_on, mdet_dp_on, dp_sdet_on, dp_ent_on, dp_new_termc,       &
    cape_ts_w

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, srf_precip,                            &
    a_land, a_sea, b_land, b_sea,                                 &
    thpixs_deep, qpixs_deep, c_mass, wcape_fac,                   &
    max_dp_thpert, min_dp_thpert, max_dp_qpert_fac

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx,  flg_up_flx_half, flg_dwn_flx,                    &
    flg_entr_up, flg_detr_up, flg_detr_up, flg_detr_dwn,          &
    flg_entr_dwn, flg_uw_dp, flg_vw_dp


USE earth_constants_mod, ONLY: g


USE water_constants_mod, ONLY: lc, lf, tm

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE PrintStatus_mod
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:


INTEGER, INTENT(IN) :: nbl      ! No. of boundary layer levels

INTEGER, INTENT(IN) :: nlev     ! No. of model layers

INTEGER, INTENT(IN) :: ntra     ! No. of tracer fields

INTEGER, INTENT(IN) :: n_cca_lev! No. of convective cloud
                                ! amount levels (1 for 2D,
                                ! nlevs for 3D)

INTEGER, INTENT(IN) :: n_dp     ! No. of deep convection points

INTEGER, INTENT(IN) :: trlev    ! No. of model levels on which
                                ! tracers are included

LOGICAL, INTENT(IN) :: bland(n_dp) ! Land/sea mask

REAL, INTENT(IN)    :: delthvu(n_dp) ! a measure of CAPE used to cal wcld

REAL, INTENT(IN)    :: exner_layer_centres(n_dp,0:nlev) !Exner

REAL, INTENT(IN)    :: exner_layer_boundaries(n_dp,0:nlev)
                                ! Exner at half level above
                                ! exner_layer_centres

LOGICAL, INTENT(IN) :: l_calc_dxek ! Switch for calculation of
                                   ! condensate increment

LOGICAL, INTENT(IN) :: l_q_interact ! Switch allows overwriting
                                    ! parcel variables when
                                    ! calculating condensate incr.

LOGICAL, INTENT(IN) :: l_tracer ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: ntml(n_dp) ! Top level of surface mixed
                                  ! layer defined relative to
                                  ! theta,q grid

INTEGER, INTENT(IN) :: ntpar(n_dp) ! Top level of initial parcel
                                   ! ascent in BL scheme defined
                                   ! relative to theta,q grid

REAL, INTENT(IN)    :: pstar(n_dp) ! Surface pressure (Pa)

REAL, INTENT(IN)    :: p_layer_centres(n_dp,0:nlev) ! Pressure (Pa)

REAL, INTENT(IN)    :: p_layer_boundaries(n_dp,0:nlev) ! Pressure
                                                       ! at half level above
                                                       ! p_layer_centres (Pa)

REAL, INTENT(IN)    ::    &
  z_theta(n_dp,nlev)      & ! height of theta levels (m)
, z_rho(n_dp,nlev)        & ! height of rho levels (m)
, r_theta(n_dp,0:nlev)    & ! radius of theta levels (m)
, r_rho(n_dp,nlev)        & ! radius of rho levels (m)
, rho_theta(n_dp,nlev)    & ! density for theta lev (kg/m3)
, rho(n_dp,nlev)          & ! density for rho lev (kg/m3)
, r2rho_th(n_dp,nlev)     & ! radius**2 density for theta lev (kg/m)
, r2rho(n_dp,nlev)        & ! radius**2 density for rho lev (kg/m)
, dr_across_th(n_dp,nlev) & ! thickness of theta levels (m)
, dr_across_rh(n_dp,nlev)   ! thickness of rho levels (m)

REAL, INTENT(IN)    :: q(n_dp,nlev) ! Model mixing ratio (kg/kg)

REAL, INTENT(IN)    :: q1_sd(n_dp) ! Standard deviation of
                                   ! turbulent flucts. of layer 1 q (kg/kg)

REAL, INTENT(IN)    :: t1_sd(n_dp) ! Standard deviation of
                                   ! turbulent flucts. of layer 1 temp. (K)

REAL, INTENT(IN)    :: th(n_dp,nlev) !Model potential temperature (K)

REAL, INTENT(IN)    :: timestep      ! Model timestep (s)

REAL, INTENT(IN)    :: u(n_dp,nlev)  ! Model U field (m/s)

REAL, INTENT(IN)    :: v(n_dp,nlev)  ! Model V field (m/s)

REAL, INTENT(IN)    :: uw0(n_dp)     ! U-comp of surface stress (N/m2)

REAL, INTENT(IN)    :: vw0(n_dp)     ! V-comp of surface stress (N/m2)

REAL, INTENT(IN)    :: wstar(n_dp) ! Convective velocity scale (m/s)

REAL, INTENT(IN)    :: w_max(n_dp) ! max w in column
                                   !for use in scale dependent cape timescale

REAL, INTENT(IN)    :: zlcl_uv(n_dp) !Lifting condensation level
                                     ! defined for the uv grid (m)

INTEGER, INTENT(IN) :: freeze_lev(n_dp) ! Level index for
                                        ! freezing level

REAL, INTENT(IN) :: recip_pstar(n_dp) ! Reciprocal of pstar array

REAL, INTENT(IN) :: qse(n_dp,nlev) ! Saturation mixing ratio of
                                   ! cloud environment (kg/kg)

! Arguments with intent INOUT:

REAL, INTENT(INOUT) :: bulk_cf(n_dp,nlev) ! Bulk total cloud volume ( )

REAL, INTENT(INOUT) :: cf_frozen(n_dp,nlev) ! Frozen water cloud volume ( )

REAL, INTENT(INOUT) :: cf_liquid(n_dp,nlev) ! Liq water cloud volume ( )

REAL, INTENT(INOUT) :: qcf(n_dp,nlev) ! Ice condensate mix ratio (kg/kg)

REAL, INTENT(INOUT) :: qcl(n_dp,nlev) ! Liq condensate mix ratio (kg/kg)

REAL, INTENT(INOUT) :: tracer(n_dp,trlev,ntra) !Model tracer fields (kg/kg)

! Arguments with intent OUT:

REAL, INTENT(OUT) :: cape_out(n_dp) ! Saved convective available
                                    ! potential energy for diagnostic
                                    ! output (J/kg)

REAL, INTENT(OUT) :: cclwp(n_dp)    ! Condensed water path (k/m2)

REAL, INTENT(OUT) :: ccw(n_dp,nlev)
                                ! Convective cloud liquid
                                ! water on model levels (kg/kg)

REAL, INTENT(OUT) :: cca(n_dp,n_cca_lev)
                                ! Convective cloud amount
                                ! on model levels (0-1)

REAL, INTENT(OUT) :: dbcfbydt(n_dp,nlev) ! Increments to
                                         ! total cld volume due to
                                         ! convection(/s)

REAL, INTENT(OUT) :: dcffbydt(n_dp,nlev) ! Increments to ice
                                         ! cloud volume due to convection (/s)

REAL, INTENT(OUT) :: dcflbydt(n_dp,nlev) ! Increments to liq
                                         ! cloud volume due to convection (/s)

REAL, INTENT(OUT) :: dqbydt(n_dp,nlev) ! Increments to q due to
                                       ! convection (kg/kg/s)

REAL, INTENT(OUT) :: dqcfbydt(n_dp,nlev) ! Increments to ice
                                         ! condensate due to convection
                                         ! (kg/kg/s)

REAL, INTENT(OUT) :: dqclbydt(n_dp,nlev) ! Increments to liq
                                         ! condensate due to convection
                                         ! (kg/kg/s)

REAL, INTENT(OUT) :: dthbydt(n_dp,nlev) ! Increments to potential
                                        ! temp. due to convection (K/s)

REAL, INTENT(OUT) :: dubydt(n_dp,nlev+1) ! Increments to U due
                                         ! to CMT (m/s2)

REAL, INTENT(OUT) :: dvbydt(n_dp,nlev+1) ! Increments to V due
                                         ! to CMT (m/s2)

REAL, INTENT(OUT) :: dtrabydt(n_dp,nlev,ntra) !Increment to tracer
                                              ! due to convection (kg/kg/s)

REAL, INTENT(OUT) :: detrain_up(n_dp,nlev) ! Fractional
                                           ! detrainment rate into updraughts
                                           ! (Pa/s)

REAL, INTENT(OUT) :: detrain_dwn(n_dp,nlev) ! Fractional
                                            ! detrainment rate into
                                            ! downdraughts (Pa/s)


REAL, INTENT(OUT) :: entrain_up(n_dp,nlev) ! Fractional
                                           ! entrainment rate into updraughts
                                           ! (Pa/s)

REAL, INTENT(OUT) :: entrain_dwn(n_dp,nlev) ! Fractional
                                            ! entrainment rate into
                                            ! downdraughts (Pa/s)

INTEGER, INTENT(OUT) :: iccb(n_dp) ! Convective cloud base level

INTEGER, INTENT(OUT) :: icct(n_dp) ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(n_dp) ! Lowest conv. cloud amt. (%)

INTEGER, INTENT(OUT) :: lcbase(n_dp) ! Lowest conv. cloud base level

INTEGER, INTENT(OUT) :: lctop(n_dp) ! Lowest conv. cloud top level

REAL, INTENT(OUT) :: rain(n_dp) ! Surface convective rainfall (kg/m2/s)

REAL, INTENT(OUT) :: snow(n_dp) ! Surface convective snowfall (kg/m2/s)

REAL, INTENT(OUT) :: rain_3d(n_dp,nlev) ! Convective rainfall flux (kg/m2/s)

REAL, INTENT(OUT) :: snow_3d(n_dp,nlev) ! Convective snowfall flux (kg/m2/s)

REAL, INTENT(OUT) :: up_flux(n_dp,nlev) ! Updraught mass flux (Pa/s)

REAL, INTENT(OUT) :: up_flux_half(n_dp,nlev) ! Updraught mass flux
                                             ! on rho levels (Pa/s)

REAL, INTENT(OUT) :: dwn_flux(n_dp,nlev) ! Downdraught mass flux (Pa/s)

REAL, INTENT(OUT) :: uw_deep(n_dp,nlev) ! X-comp. of stress
                                        ! from deep convection (kg/m/s2)

REAL, INTENT(OUT) :: vw_deep(n_dp,nlev) ! Y-comp. of stress
                                        ! from deep convection (kg/m/s2)

INTEGER, INTENT(OUT) :: kterm(n_dp) ! Level at which deep
                                    ! convection terminates,
                                    ! required by mid level scheme
REAL, INTENT(OUT) :: tcw(n_dp)   ! Total condensed water(kg/m2/s)
                                 ! required by mid-level CCA cal.

REAL, INTENT(OUT) :: cca_2d(n_dp) ! 2D convective cloud amount (%)

! Adaptive detrainment output variables

REAL, INTENT(OUT) :: rbuoy_p_out(n_dp,nlev)  !buoyancy excess

REAL, INTENT(OUT) :: the_out(n_dp,nlev)    !th_E in parcel routine

REAL, INTENT(OUT) :: thp_out(n_dp,nlev)    !th_P in parcel routine

REAL, INTENT(OUT) :: qe_out(n_dp,nlev)     !q_E in parcel routine

REAL, INTENT(OUT) :: qp_out(n_dp,nlev)     !q_P in parcel routine

REAL, INTENT(OUT) :: ind_cape_reduced(n_dp) ! 1 - if CAPE reduced
                                            ! applies to CAPE option 3
REAL, INTENT(OUT) :: cape_ts_used(n_dp)     ! CAPE timescale used for
                                            ! deep convection (s)
REAL, INTENT(OUT) :: cfl_limited(n_dp)      ! Indicator of CFL limited conv

REAL, INTENT(OUT) :: ind_deep(n_dp)         ! 1.0 if real deep convection

INTEGER, INTENT(OUT) :: error_point     ! 0 no error
                                        ! location of problem deep point
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Adaptive detrainment output variables

REAL :: rbuoy_p_here(n_dp)     !buoyancy excess

REAL :: the_here(n_dp)         !th_E in parcel routine

REAL :: thp_here(n_dp)         !th_P in parcel routine

REAL :: qe_here(n_dp)          !q_E in parcel routine

REAL :: qp_here(n_dp)          !q_P in parcel routine

REAL :: rbuoy_p_old(n_dp)      !buoyancy excess on previous level

REAL :: zk(n_dp)               !heights for use in calc
REAL :: zkp12(n_dp)            !of moist static energy
REAL :: zkp1(n_dp)

INTEGER :: index1(n_dp),index2(n_dp)

INTEGER :: ncposs               ! No. of points which may convect

INTEGER :: nconv                ! No. of convecting points

REAL :: amdetk(n_dp)            ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

REAL :: b_calc                  ! Coefficient in thpert calc.

REAL :: c_calc                  ! Coefficient in thpert calc.

REAL :: cape(n_dp)              ! Convective available potential
                                ! energy (J/kg)

REAL :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

REAL :: dcpbydt(n_dp)           ! Rate of change of cape (J/kg/s)

REAL :: depth(n_dp)             ! Depth of convective cloud (m)

REAL :: delexkp1(n_dp)          ! Difference in exner ratio
                                ! across layer k+1

REAL :: dqsthk(n_dp)            ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k (kg/kg/K)

REAL :: dqsthkp1(n_dp)          ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k+1 (kg/kg/K)

REAL :: ekp14(n_dp)             ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekp34(n_dp)             ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekm14(n_dp)             ! Entrainment coefficients at
                                ! level k-1+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: exk(n_dp)               ! Exner ratio at layer k

REAL :: exkp1(n_dp)             ! Exner ratio at layer k+1

REAL :: flxmax(n_dp)            ! Maximum initial convective
                                ! mass flux (Pa/s)

REAL :: flx_init(n_dp)          ! Initial mass flux at cloud base
                                ! (Pa/s)

REAL :: flx_init_new(n_dp)      ! flx_init scaled to destroy cape
                                ! over timescale cape_timescale (Pa/s)

REAL :: flxmax_init(n_dp)       ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer
                                ! in Pa/s)

REAL :: max_cfl(n_dp)           ! Max cfl ratio over a convecting
                                ! layer

REAL :: p_lcl(n_dp)             ! Pressure at LCL (Pa)

REAL :: precip(n_dp,nlev)       ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

REAL :: pk(n_dp)                ! Pressure at midpoint of layer
                                ! k (Pa)

REAL :: pkp1(n_dp)              ! Pressure at midpoint of layer
                                ! k+1 (Pa)

REAL :: delpk(n_dp)             ! Pressure difference over layer
                                ! k (Pa)

REAL :: delpkp1(n_dp)           ! Pressure difference over layer
                                ! k+1 (Pa)

REAL :: delpkp12(n_dp)          ! Pressure difference between
                                ! layers k and k+1 (Pa)

REAL :: delp_uv_k(n_dp)         ! Pressure difference across uv
                                ! layer k (Pa)

REAL :: delp_uv_kp1(n_dp)       ! Pressure difference across uv
                                ! layer k+1 (Pa)

REAL :: q_lcl(n_dp)             ! Mixing ratio at LCL (kg/kg)

REAL :: qse_lcl(n_dp)           ! Saturated q at LCL (kg/kg)

REAL :: rhum(n_dp)              ! Dummy relative humidity
                                ! (only used on shallow points)

REAL :: t_lcl(n_dp)             ! Temperature at LCL (K)

REAL :: th_lcl(n_dp)            ! Theta at LCL (K)

REAL :: thv_pert                ! Theta_v parcel perturbation (K)

REAL :: thpert(n_dp)            ! Theta parcel perturbation (K)

REAL :: qpert(n_dp)             ! q parcel perturbation (kg/kg)

INTEGER :: start_lev3c(n_dp)    ! PC2 Compressed convection
                                ! initiation level

LOGICAL :: l_shallow(n_dp)      ! Dummy variable (=.F.)

LOGICAL :: l_mid(n_dp)          ! Dummy variable (=.F.)

LOGICAL :: cumulus(n_dp)        ! Dummy variable (=.T.)

LOGICAL :: bgmk(n_dp)           ! Mask for points where parcel in
                                ! layer k is saturated

LOGICAL :: bwater(n_dp,2:nlev)  ! Mask for points at which
                                ! condensate is liquid
LOGICAL :: bwk(n_dp)            !mask for liquid condensate on k

LOGICAL :: bwkp1(n_dp)          !mask for liquid condensate on k+1

LOGICAL :: blowst(n_dp)         ! Dummy variable indicating low
                                ! enough stability for convection
                                ! to occur

LOGICAL :: bterm(n_dp)          ! Mask for points which have
                                ! stopped convecting

LOGICAL :: bconv(n_dp)          ! Mask for points at which
                                ! convection is occurring

LOGICAL :: bcposs(n_dp)         ! Mask for points passing
                                ! initial stability test

! Parcel variables

REAL :: qpi(n_dp)               ! Initial parcel mixing ratio (kg/kg)

REAL :: qp(n_dp,nlev)           ! Parcel mixing ratio (kg/kg)

REAL :: thpi(n_dp)              ! Initial parcel potential temp.(K)

REAL :: thp(n_dp,nlev)          ! Parcel potential temp (K)

REAL :: trap(n_dp,nlev,ntra)    ! Tracer content of parcel (kg/kg)

REAL :: expi(n_dp)              ! Initial parcel exner pressure

REAL :: xpk(n_dp,nlev)          ! Parcel cloud water (kg/kg)

REAL :: flx(n_dp,nlev)          ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(n_dp,nlev)     ! Minmum parcel buoyancy excess

REAL :: thpixs_v(n_dp,nlev)     ! Theta parcel excess (K)

REAL :: qpixs_v(n_dp,nlev)      ! Q parcel excess(kg/kg)

! PC2

REAL :: qclp(n_dp,nlev)         ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

REAL :: qcfp(n_dp,nlev)         ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)

! Parameters


REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed


! CMT variables  - those used depend on scheme

! Required by Gregory-Kershaw scheme operating in plume calculation

INTEGER ::         &
 nstart(n_dp)        ! Level for start of plume

REAL ::            &
 eflux_u_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     !  top of layer (Pa m/s)
,eflux_v_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     ! bottom of layer (Pa m/s)
,up(n_dp,nlev)     & ! Parcel U (m/s)
,vp(n_dp,nlev)     & ! Parcel V (m/s)
,zsurf(n_dp)         ! Height of start of plume = 0.1*zlcl

! Required by Turbulence base scheme called after plume calculation

INTEGER ::             &
 nlcl_uv(n_dp)         & ! Level index for LCL
,ntop_uv(n_dp)         & ! Level index for top of layer
,n_0degc(n_dp)         & ! Level index for zero degrees
,cu_term(n_dp)         & ! Indicies for CMT subs
,cu_tend(n_dp)           ! Indicies for CMT subs

REAL ::                &
 mass_dwn(nlev,n_dp)   & ! Downdraught mass flux (Pa/s)
,p_uv(nlev,n_dp)       & ! Pressure of model level (Pa)
,phalf_uv(nlev,n_dp)   & ! Pressure of half level (Pa)
,plcl_uv(n_dp)         & ! Pressure at LCL (Pa)
,ptop_uv(n_dp)         & ! Pressure at top of cloud layer (Pa)
,p_0degc_uv(n_dp)      & ! Pressure of zero degree level (Pa)
,rho_uv(nlev,n_dp)     & ! Density on uv level (kg/m3)
,visc(nlev,n_dp)       & ! CMT eddy viscosity (m2/s)
,uw(nlev,n_dp)         & ! U- comp stress profile (N/m2)
                         ! (units vary through calls)
,vw(nlev,n_dp)         & ! V-comp stress profile (N/m2)
,uw_base(nlev,n_dp)    & ! Cloud base U stress (N/m2)
,vw_base(nlev,n_dp)    & ! Cloud base V stress (N/m2)
,ue_p(nlev,n_dp)       & ! Environment U profile (m/s)
,ve_p(nlev,n_dp)         ! Environment V profile (m/s)

REAL :: exk_temp                ! Temporary exner

! Required by all version of CMT

REAL :: flxkp12(nlev,n_dp)      ! Mass flux on half level (Pa/s)

REAL :: mb(n_dp)                ! Cloud base mass flux (Pa/s)

LOGICAL :: l_mom_gk             ! true if Gregory-Kershaw CMT required

! Cape scaling/closure variables

INTEGER :: det_lev(n_dp)        ! Level at which split final
                                ! detrainment last occurred

INTEGER :: nterm                ! No. of points where conv.
                                ! has terminated

INTEGER :: index_nterm(n_dp)    ! Index for points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage

REAL :: scale_f(n_dp)           ! scale factor

REAL :: cape_ts_new(n_dp)       ! Used as variable in RH-based
                                ! closure

REAL :: relh(n_dp)              ! RH integral (average when
                                ! convection terminates)

REAL :: dptot(n_dp)             ! Delta P integral

! Downdraught scheme variables

INTEGER :: npossdd              ! Max. no. of downdraughts
                                ! possible

INTEGER :: nnodd                ! No. of downdraughts not possible

INTEGER :: index_possdd(n_dp)   ! Index of downdraughts possible

INTEGER :: index_nodd(n_dp)     ! Index of downdraughts not
                                ! possible
INTEGER :: kmax_term            ! maximum termination level + 1

REAL :: deltap_cld              ! Pressure thickness of convective
                                ! cloud (Pa)

! Local compressed arrays

LOGICAL :: bconv_c2(n_dp)

LOGICAL :: bgmkp1_c(n_dp), bgmkp1_c2(n_dp) ! Mask for points
                                           ! where parcel in layer k+1
                                           ! is saturated

LOGICAL :: bwk_c(n_dp), bwk_c2(n_dp) ! bwater mask in layer k

LOGICAL :: bwkp1_c(n_dp), bwkp1_c2(n_dp) ! bwater mask in layer k+1

REAL :: deltak_c2(n_dp)         ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

REAL :: dqek_c2(n_dp)           ! Increment to q due to
                                ! convection in layer k (kg/kg)

REAL :: dqekp1_c2(n_dp)         ! Increment to q due to
                                ! convection in layer k+1 (kg/kg)

REAL :: dthek_c2(n_dp)          ! Increment to potential temp.
                                ! due to convection in layer k

REAL :: dthekp1_c2(n_dp)        ! Increment to potential temp.
                                ! due to convection in layer k+1

REAL :: dtraek_c2(n_dp,ntra)    ! Increment to model tracer due
                                ! to conv. at level k (kg/kg/s)

REAL :: dtraekp1_c2(n_dp,ntra)  ! Increment to model tracer due
                                ! to conv. at level k+1 (kg/kg/s)

REAL :: duek_c2(n_dp)           ! Increment to model U in layer k
                                ! due to CMT (m/s2)

REAL :: duekp1_c2(n_dp)         ! Increment to model U in layer
                                ! k+1 due to CMT (m/s2)

REAL :: dvek_c2(n_dp)           ! Increment to model V in layer k

REAL :: dvekp1_c2(n_dp)         ! Increment to model V in layer
                                ! k+1 due to CMT (m/s2)

REAL :: flxk_c(n_dp), flxk_c2(n_dp) !Parcel mass flux in layer k (Pa/s)

REAL :: flxkp12_c2(n_dp)        ! Half level mass flux (Pa/s)

REAL :: prekp1_c2(n_dp)         ! Precip. from parcel as it rises
                                ! from layer k to k+1 (kg/m2/s)

REAL :: qpk_c(n_dp), qpk_c2(n_dp) ! Parcel mixing ratio in
                                ! layer k(kg/kg)
REAL :: qpk(n_dp)

REAL :: qpkp1_c(n_dp), qpkp1_c2(n_dp) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

REAL :: qek_c(n_dp), qek_c2(n_dp) ! Env. mixing ratio in
                                  ! layer k (kg/kg)
REAL :: qek(n_dp)

REAL :: qekp1_c(n_dp), qekp1_c2(n_dp) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)
REAL :: qekp1(n_dp)

REAL :: qsek_c2(n_dp)           ! Saturation mixing ratio of
                                ! cld. env. in layer k (kg/kg)
REAL :: qsek(n_dp)

REAL :: qsekp1_c(n_dp), qsekp1_c2(n_dp) ! Saturation mixing ratio
                                        ! of cld. env. in layer k+1 (kg/kg)
REAL :: qsekp1(n_dp)

REAL :: thek_c(n_dp), thek_c2(n_dp) ! Env. potential temp in layer k (K)

REAL :: thek(n_dp)

REAL :: thekp1_c(n_dp), thekp1_c2(n_dp) ! Env. potential temp i in layer k (K)

REAL :: thekp1(n_dp)

REAL :: thpk_c(n_dp), thpk_c2(n_dp) ! Parcel potential temp in layer k (K)

REAL :: thpk(n_dp)

REAL :: thpkp1_c(n_dp), thpkp1_c2(n_dp)! Parcel potential temp in layer k (K)

REAL :: thpkp1(n_dp)

REAL :: traek_c(n_dp,ntra), traek_c2(n_dp,ntra) ! Tracer content
                                                ! cld. env. in layer k (kgkg-1)

REAL :: traekp1_c(n_dp,ntra), traekp1_c2(n_dp,ntra) ! Tracer
                                                    ! content of cloud env.
                                                    ! in layer k+1 (kg/kg)

REAL :: trapk_c(n_dp,ntra), trapk_c2(n_dp,ntra) ! Tracer cont.
                                                ! of parcel in layer k (kg/kg)

REAL :: trapkp1_c(n_dp,ntra), trapkp1_c2(n_dp,ntra) ! Tracer cont.
                                       ! of parcel in layer k+1 (kg/kg)

REAL :: rbuoy_c(n_dp), rbuoy_c2(n_dp) ! Buoyancy of parcel at k+1 (Kelvin)

REAL :: uek_c(n_dp), uek_c2(n_dp) ! Model U field on layer k (m/s)

REAL :: uekp1_c(n_dp), uekp1_c2(n_dp)! Model U field on layer k+1 (m/s)

REAL :: vek_c(n_dp), vek_c2(n_dp) ! Model V field on layer k (m/s)

REAL :: vekp1_c(n_dp), vekp1_c2(n_dp) ! Model V field on layer k+1 (m/s)

REAL :: upk_c(n_dp), upk_c2(n_dp) ! Parcel U in layer k
                                  ! after entrainment (m/s)

REAL :: upkp1_c(n_dp), upkp1_c2(n_dp) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

REAL :: vpk_c(n_dp), vpk_c2(n_dp) ! Parcel V in layer k
                                  ! after entrainment (m/s)

REAL :: vpkp1_c(n_dp), vpkp1_c2(n_dp) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

REAL :: xsqkp1_c(n_dp), xsqkp1_c2(n_dp) ! Excess water vapour
                                        ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

REAL :: qclek_c(n_dp), qclek_c2(n_dp) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qclekp1_c(n_dp), qclekp1_c2(n_dp) ! Environment liquid
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qcfek_c(n_dp), qcfek_c2(n_dp) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qcfekp1_c(n_dp), qcfekp1_c2(n_dp) ! Environment frozen
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qclpk_c(n_dp), qclpk_c2(n_dp) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qclpkp1_c(n_dp), qclpkp1_c2(n_dp) ! Parcel liquid
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qcfpk_c(n_dp), qcfpk_c2(n_dp) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qcfpkp1_c(n_dp), qcfpkp1_c2(n_dp) ! Parcel frozen
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: cflek_c2(n_dp),cflekp1_c2(n_dp) ! Environment liquid water
                                        ! cloud volume ( )

REAL :: cffek_c2(n_dp),cffekp1_c2(n_dp) ! Environment frozen water
                                        ! cloud volume ( )

REAL :: bcfek_c2(n_dp),bcfekp1_c2(n_dp) ! Environment bulk total
                                        ! cloud volume ( )

REAL :: dqclek_c2(n_dp),dqclekp1_c2(n_dp) ! Environment increments
                                          ! to liquid condensate mixing
                                          ! ratio to convection (kg/kg/s)

REAL :: dqcfek_c2(n_dp),dqcfekp1_c2(n_dp) ! Environment increments
                                          ! to frozen condensate mixing
                                          ! ratio to convection (kg/kg/s)

REAL :: dcflek_c2(n_dp),dcflekp1_c2(n_dp) ! Environment increments
                                          ! to liquid water cloud volume due
                                          ! to convection (/s)

REAL :: dcffek_c2(n_dp),dcffekp1_c2(n_dp) ! Environment increments
                                          ! to frozen water cloud volume due
                                          ! to convection (/s)

REAL :: dbcfek_c2(n_dp),dbcfekp1_c2(n_dp) ! Environment increments
                                          ! to bulk total cloud volume due
                                          ! to convection (/s)

REAL :: amdetk_c2(n_dp)
LOGICAL :: bgmk_c2(n_dp)
LOGICAL :: bland_c2(n_dp)
LOGICAL :: blowst_c2(n_dp)
LOGICAL :: bterm_c2(n_dp)
REAL :: cape_c2(n_dp)
REAL :: cca_2d_c2(n_dp)
REAL :: cclwp_c2(n_dp)
REAL :: ccw_c2(n_dp)
LOGICAL :: cumulus_c(n_dp), cumulus_c2(n_dp)
REAL :: dcpbydt_c2(n_dp)
REAL :: delexkp1_c2(n_dp)
REAL :: delpk_c2(n_dp)
REAL :: delpkp1_c2(n_dp)
REAL :: delp_uv_k_c2(n_dp)
REAL :: delp_uv_kp1_c2(n_dp)
REAL :: depth_c2(n_dp)
REAL :: dptot_c2(n_dp)
REAL :: dqsthkp1_c2(n_dp)
REAL :: dqsthk_c2(n_dp)
REAL :: eflux_u_ud_c2(n_dp)
REAL :: eflux_v_ud_c2(n_dp)
REAL :: ekp14_c(n_dp),ekp14_c2(n_dp)
REAL :: ekp34_c(n_dp),ekp34_c2(n_dp)
REAL :: exk_c2(n_dp)
REAL :: exkp1_c(n_dp),exkp1_c2(n_dp)
REAL :: expi_c2(n_dp)
INTEGER :: icct_c2(n_dp)
INTEGER :: iccb_c2(n_dp)
INTEGER :: lctop_c2(n_dp)
INTEGER :: lcbase_c2(n_dp)
REAL :: lcca_c2(n_dp)
LOGICAL :: l_shallow_c2(n_dp)
LOGICAL :: l_mid_c2(n_dp)
REAL :: max_cfl_c2(n_dp)
REAL :: pk_c(n_dp),pk_c2(n_dp)
REAL :: pkp1_c(n_dp),pkp1_c2(n_dp)
REAL :: pstar_c2(n_dp)
REAL :: q1_sd_c2(n_dp)
REAL :: qpi_c2(n_dp)
REAL :: qpixs_v_c2(n_dp)
REAL :: relh_c2(n_dp)
REAL :: rbuoy_p_here_c2(n_dp)
REAL :: the_here_c2(n_dp)
REAL :: thp_here_c2(n_dp)
REAL :: qe_here_c2(n_dp)
REAL :: qp_here_c2(n_dp)
REAL :: rbuoy_p_old_c2(n_dp)
REAL :: tcw_c2(n_dp)
REAL :: thpi_c2(n_dp)
REAL :: thpixs_v_c2(n_dp)
REAL :: t1_sd_c2(n_dp)
REAL :: xpk_c(n_dp),xpk_c2(n_dp)
REAL :: xsbmin_v_c2(n_dp)
LOGICAL :: b_nodd(n_dp)   ! points with no downdraught
LOGICAL :: b_dd(n_dp)     ! points with downdraught on termination

REAL ::      &
  rh_test    &   ! critical RH value for CAPE_OPT=6
 ,rh_fac         ! factor for CAPE_OPT=6 calculation

REAL ::         &
  wcld(n_dp)    &  ! Convective veloicty scale
 ,zlcl(n_dp)       ! lifting condensation level


! Loop counters

INTEGER :: i,j,k,ktra,kt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('DEEP_CONV_4A',zhook_in,zhook_handle)

! Initialise error_point

error_point=0     

!initialise SCM diagnostics


DO k = 1,nlev
  DO i = 1,n_dp
    rbuoy_p_out(i,k)=0.0
    the_out(i,k)=th(i,k)
    thp_out(i,k)=th(i,k)
    qe_out(i,k)=q(i,k)
    qp_out(i,k)=q(i,k)
  END DO
END DO



! Initialise logicals


DO i = 1,n_dp
  blowst(i)    = .TRUE.
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  cumulus(i)   = .TRUE.
  l_shallow(i) = .FALSE.
  l_mid(i)     = .FALSE.
  b_nodd(i)    = .FALSE.
  b_dd(i)      = .FALSE.
END DO

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays (now in glue)
!-----------------------------------------------------------------------

! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)

DO k = 1,nlev-1
  DO i = 1,n_dp
    xsbmin_v(i,k) = MIN( ((p_layer_centres(i,k) -               &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = MIN( ((p_layer_centres(i,k) -               &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_deep

    qpixs_v(i,k)  = qpixs_deep
  END DO
END DO  ! nlev

! Calculate cloud base mass flux

DO i = 1,n_dp
  mb(i) = c_mass * wstar(i)
END DO

! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.


DO i = 1,n_dp
  k=ntml(i)
  p_lcl(i)  = p_layer_boundaries(i,k)
  th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
  t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
  q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
END DO

! Calculate saturation mixing ratio at LCL

! DEPENDS ON: qsat
CALL qsat(qse_lcl,t_lcl,p_lcl,n_dp)


!----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!----------------------------------------------------------------------
IF (l_mom) THEN
  SELECT CASE (deep_cmt_opt)

  CASE(2)         ! Gregory-Kershaw CMT

    ! need level near surface for initial parcel U & V values
    ! zsurf = 0.1*z_lcl

    DO i = 1,n_dp
      zsurf(i)  = 0.1*z_rho(i,ntml(i))
    END DO
    DO k=nlev-1,1,-1
      DO i = 1,n_dp
       IF (zsurf(i) <= z_theta(i,k)) THEN
         nstart(i) = k
       END IF
      END DO
    END DO
    l_mom_gk = .TRUE.

  CASE DEFAULT    ! (0/1/5) Alan Grant's eddy viscosity based CMT

! Note: In terms of array indices p and phalf follow the convention
!       used in the boundary layer scheme. phalf(k,*) refers to the
!       lower boundary of uv layer k. This follows the convention for
!       um UM4.5 and before

!       Also note that p_layer_boundaries(0) and p_layer_centres(0)
!       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.

!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)


    k=1
    DO i = 1,n_dp
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      nlcl_uv(i)    = ntml(i) + 1
      n_0degc(i)    = freeze_lev(i)
      kterm(i)=0
    END DO

    DO i = 1,n_dp
      DO k = 2,nlev
        p_uv(k,i)     = p_layer_boundaries(i,k-1)
        phalf_uv(k,i) = p_layer_centres(i,k-1)
        ue_p(k,i)     = u(i,k)
        ve_p(k,i)     = v(i,k)
        exk_temp      = (p_uv(k,i)/pref)**0.288
        rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp *        &
                        (th(i,k-1) + th(i,k)))
      END DO
      plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
      p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
      rho_uv(1,i)     = rho_uv(2,i)
    END DO

    l_mom_gk = .FALSE.

  END SELECT      ! test on deep_cmt_opt

ELSE

! Initialise variable
  l_mom_gk = .FALSE.

END IF     !L_mom

!-----------------------------------------------------------------------
! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! Re-set th and q xs's at ntml

IF (limit_pert_opt == 1 .OR. limit_pert_opt == 2) THEN
  DO i = 1,n_dp
    IF (t_lcl(i) >  tm) THEN
      dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                    / (r * t_lcl(i) * t_lcl(i))
    ELSE
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                    / (r * t_lcl(i) * t_lcl(i))
    END IF

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                    + c_virtual * qse_lcl(i)

    thv_pert = -0.5 * (th(i,ntml(i)+1)                            &
                    * (1.0+c_virtual * q(i,ntml(i)+1))            &
                    - th(i,ntml(i)) * (1.0 + c_virtual            &
                    * q(i,ntml(i)))) + 0.5

    c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i)                &
                    - q_lcl(i)) - thv_pert

    thpert(i) = MAX(MIN(-c_calc / b_calc, max_dp_thpert),         &
                          min_dp_thpert)
                          ! ignore term in thpert**2

    thpixs_v(i,ntml(i)) = thpert(i)

    qpert(i)  = MAX(MIN(qse_lcl(i) + ((p_lcl(i) / pref)           &
                          **kappa) * thpert(i) * dq_sat_env       &
                          - q_lcl(i),                             &
                          max_dp_qpert_fac * qse_lcl(i)),0.0)

    qpixs_v(i,ntml(i))  = qpert(i)

  END DO ! n_dp

ELSE IF (limit_pert_opt == 0) THEN

  DO i = 1,n_dp
    IF (t_lcl(i) >  tm) THEN
      dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                    / (r * t_lcl(i) * t_lcl(i))
    ELSE
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                    / (r * t_lcl(i) * t_lcl(i))
    END IF

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                    + c_virtual * qse_lcl(i)

    thv_pert = -0.5 * (th(i,ntml(i)+1)                            &
                    * (1.0+c_virtual * q(i,ntml(i)+1))            &
                    -th(i,ntml(i)) * (1.0 + c_virtual             &
                    * q(i,ntml(i)))) + 0.5

    c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i)                &
                    - q_lcl(i)) - thv_pert

    thpert(i) = -c_calc / b_calc   ! ignore term in thpert**2

    thpixs_v(i,ntml(i)) = thpert(i)

    qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)                   &
                           **kappa) * thpert(i) * dq_sat_env      &
                             - q_lcl(i)

    qpixs_v(i,ntml(i))  = qpert(i)

  END DO ! n_dp

END IF


! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)


! DEPENDS ON: flag_wet
CALL flag_wet(n_dp,n_dp,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

DO k=1, nlev
  DO i=1, n_dp
    precip(i,k)   = 0.0
    ccw(i,k)      = 0.0
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
  END DO
END DO

DO k=1, n_cca_lev
  DO i=1, n_dp
    cca(i,k) = 0.0
  END DO
END DO


IF (l_mom) THEN
  DO k = 1,nlev+1
    DO i = 1,n_dp
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

! No need to initialise dtrabydt as done in glue_conv

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
DO i = 1,n_dp
  ind_cape_reduced(i) = 0.0
  cape_ts_used(i) = 0.0
  cfl_limited(i)      = 0.0
  ind_deep(i)     = 1.0  ! not setup to detect false deep convection
END DO

IF (flg_up_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux_half(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      dwn_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      entrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      detrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      entrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      detrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF

IF (l_mom) THEN
  IF (flg_uw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        uw_deep(i,k) = 0.0
      END DO
    END DO
  END IF
  IF (flg_vw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        vw_deep(i,k) = 0.0
      END DO
    END DO
  END IF
END IF  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------

DO i = 1,n_dp
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
END DO

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling and closure calculations
!-----------------------------------------------------------------------

DO i = 1,n_dp
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)      = 0.0
  cape_out(i)  = 0.0
  dcpbydt(i)   = 0.0
  max_cfl(i)   = 0.0
  det_lev(i)   = 0
  relh(i)       = 0.0
  dptot(i)      = 0.0

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


! Set rhum to an arbitrary value since it is not used in LAYER_CN
! for deep convection. Moved outside level loop as values unchanged.


!        rhum(i) = 0.0
  rhum(i) = q(i,1)/qse(i,1)

END DO



! set SCM adaptive diagnostics for level k = 1

DO i = 1,n_dp
  rbuoy_p_out(i,1) = 0.0
  the_out(i,1) = th(i,1)
  thp_out(i,1) = th(i,1)
  qe_out(i,1) = q(i,1)
  qp_out(i,1) = q(i,1)
END DO
!Also, initialise rbuoy_p_old, i.e. for previous level, to 0.0
!for level1
DO i = 1,n_dp
  rbuoy_p_old(i) = 0.0
END DO
!initialise ekm14
DO i =1, n_dp
   ekm14(i) =0.0
END DO

!Initialise adaptive entrainment variables
!intitilaise to level 2 'cos that's where parcel lift starts from
DO i = 1, n_dp
  thek(i)=th(i,2)
  qek(i)=q(i,2)
  qsek(i)=qse(i,2)
  thekp1(i)=th(i,2)
  qekp1(i)=q(i,2)
  qsekp1(i)=qse(i,2)
  thpk(i)=th(i,2)
  qpk(i)=q(i,2)
  bwk(i)=bwater(i,2)
  bwkp1(i)=bwater(i,2)
!Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1
  zk(i) = z_theta(i,2)
  zkp12(i)=z_rho(i, 3)
  zkp1(i)=z_theta(i, 3)
END DO

!intialise parcel values over all levels to make sure we have
!non-garbage values at points which don't convect (do not intialise
!inside loop, 'cos some values are set at the end of each pass at k+1
!and must not be overwritten at k at the start of the next pass).
DO k = 2, nlev -1
  DO i = 1, n_dp
    qp(i,k) = q(i,k)
    thp(i,k) = th(i,k)
  END DO
END DO

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
DO k = 2,nlev-1

  DO i = 1, n_dp
    rhum(i) = q(i,k)/qse(i,k)
  END DO


!  Initialise adaptive diagnostics for this level

  DO i = 1,n_dp
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
  DO i = 1, n_dp
    thek(i)=th(i,k)
    qek(i)=q(i,k)
    qsek(i)=qse(i,k)
    thekp1(i)=th(i,k+1)
    qekp1(i)=q(i,k+1)
    qsekp1(i)=qse(i,k+1)
    !note for all levels above k=2, thp(i, k) for this current pass
    !will have been set as thp(i, k+1) at end of previous pass through
    !loop, and similarly for qp (where point has not convected before
    !will be set to environmental value)
    IF(k  ==  2) THEN      !set to environmental values
      thpk(i)=th(i,2)
      qpk(i)=q(i,2)
    ELSE                   !value from previous pass through loop
      thpk(i)=thp(i,k)
      qpk(i)=qp(i,k)
    END IF
    bwk(i)=bwater(i,k)
    bwkp1(i)=bwater(i,k+1)
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    zk(i) = z_theta(i,k)
    zkp12(i)=z_rho(i, k+1)
    zkp1(i)=z_theta(i, k+1)
  END DO

! Set initial parcel properties for thp, qp if convection
! is not occurring at level k

  DO i = 1,n_dp
  ! not convecting and not convected in column before
  !NB Calc for thp, qp have been moved to before call to LAYER_CN for
  !adaptive mod,
  !and settings for thpk, qpk added, but should not change other calcs
  !which are carried out as before after call to layer_cn.

    IF ( .NOT. bconv(i).AND.det_lev(i) == 0) THEN
      thpi(i)  = th(i,k) + thpixs_v(i,k)
      thp(i,k) = thpi(i)
      thpk(i)=thp(i,k)         !set thpk for first convecting
                               !level
      qpi(i)   = q(i,k) + qpixs_v(i,k)
      qp(i,k)  = qpi(i)
      qpk(i) = qp(i,k)        !set qpk for first convecting level
    END IF
  END DO  ! n_dp

!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn_4a
  CALL layer_cn_4a(k,n_dp,n_dp,nlev,                              &
                exner_layer_boundaries,                           &
                exner_layer_centres,                              &
                p_layer_boundaries,                               &
                p_layer_centres,                                  &
                pstar,pk,pkp1,delpk,delpkp1,                      &
                delpkp12,ekp14,ekp34,amdetk,exk,exkp1,            &
                delexkp1,delp_uv_k,delp_uv_kp1,recip_pstar,       &
                rhum,l_shallow,ntml,ntpar                         &
                ,cumulus, mdet_dp_on, dp_ent_on                   &
                ,bconv                                            &
                ,thek, qek,qsek, thekp1,qekp1,qsekp1              &
                ,thpk, qpk                                        &
                ,bwk, bwkp1, ekm14                                &
                ,zk, zkp12, zkp1                                  &
                 )

! Set ekm14 for next pass through loop
   DO i = 1, n_dp
     ekm14(i) = ekp14(i)
   END DO



! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)


  IF (k == 2) THEN
! DEPENDS ON: dqs_dth_4a
    CALL dqs_dth_4a(dqsthk,k,th(1,k),qse(1,k),exk,n_dp)
  ELSE
    DO i = 1,n_dp
      dqsthk(i) = dqsthkp1(i)
    END DO
  END IF

! DEPENDS ON: dqs_dth
  CALL dqs_dth_4a(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_dp)


! Set other grid dependent constants


  DO i = 1,n_dp


! Maximum initial convective mass flux


    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  END DO  ! n_dp


! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k

  DO i = 1,n_dp
! not convecting and not convected in column before
! PC2 qclp and qcfp zero at this point but will add an initial
! value at cloud base
!NB init. for thp, qp has been moved to before call to LAYER_CN
!for adaptive mod,
!and settings for thpk, qpk added, but should not change other calcs.

    IF ( .NOT. bconv(i).AND.det_lev(i) == 0) THEN
      expi(i)  = exk(i)
      xpk(i,k)  = 0.0
      qclp(i,k) = 0.0
      qcfp(i,k) = 0.0
      flx(i,k) = 0.0
      bgmk(i)  = .FALSE.
      depth(i) = 0.0

      IF (l_mom_gk) THEN  ! Gregory Kershaw CMT

! Set initial parcel values at cloud base to values of near surface winds
          up(i,k)   = u(i,nstart(i))
          vp(i,k)   = v(i,nstart(i))

      END IF
    END IF
  END DO  ! n_dp

  IF (l_tracer) THEN
    DO ktra=1,ntra
      DO i = 1,n_dp
        IF ( .NOT. bconv(i)) THEN
          trap(i,k,ktra)  = tracer(i,k,ktra)
        END IF  !not bconv
      END DO
    END DO
  END IF


! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k = ntml

  DO i = 1,n_dp

    bcposs(i) = bconv(i) .OR. k  ==  ntml(i)

  END DO  ! n_dp


! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)

  ncposs = 0
  DO i = 1,n_dp
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
      cumulus_c(i) = .TRUE.         ! bgmkp1_c out
      uek_c(i)    = u(index1(i),k)
      uekp1_c(i)  = u(index1(i),k+1)
      vek_c(i)    = v(index1(i),k)
      vekp1_c(i)  = v(index1(i),k+1)
      upk_c(i)    = up(index1(i),k)
      vpk_c(i)    = vp(index1(i),k)
    END DO
    IF (l_q_interact .OR. cnv_wat_load_opt == 1) THEN
      DO i = 1,ncposs
!           PC2 variables or needed for water loading
        qclekp1_c(i)  = qcl(index1(i),k+1)
        qcfekp1_c(i)  = qcf(index1(i),k+1)
      END DO
    END IF
    IF (l_q_interact) THEN
      DO i = 1,ncposs
!           PC2 variables
        qclek_c(i)    = qcl(index1(i),k)
        qcfek_c(i)    = qcf(index1(i),k)
        qclpk_c(i)    = qclp(index1(i),k)
        qcfpk_c(i)    = qcfp(index1(i),k)
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,ncposs
          trapk_c(i,ktra) = trap(index1(i),k,ktra)
          traek_c(i,ktra) = tracer(index1(i),k,ktra)
          traekp1_c(i,ktra) = tracer(index1(i),k+1,ktra)
        END DO
      END DO
    END IF
  END IF  ! ncposs>0

!-----------------------------------------------------------------------
! 3.2  Lift parcel from layer k to layer k+1
!-----------------------------------------------------------------------

! DEPENDS ON: lift_par_4a
  CALL lift_par_4a(ncposs,n_dp,thpkp1_c,qpkp1_c,xsqkp1_c,         &
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

! Allow parcel to convect from ntml.
! Flag convecting points with logical array bconv

    IF (k  ==  ntml(index1(i))) THEN
      bconv(index1(i)) = .TRUE.  ! convection active
      blowst(index1(i)) = .TRUE. ! convection initialised in layer

! Set parcel mass flux (UM Documentation paper 27, section 1.5)
! If mass flux out of the initial layer is greater than the mass flux
! of the layer over the timestep then limit mass flux to mass of layer.

      flxk_c(i) = mb(index1(i)) * g *                             &
                      p_layer_centres(index1(i),k) / (r *         &
                        thpk_c(i) * (p_layer_centres(index1(i),k) &
                           / pref)**kappa)

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

    ELSE     ! k=ntml test
      blowst(index1(i))=.FALSE. ! not initialised in this layer
    END IF   ! k=ntml test

! Reset threshold for forced detrainment
! to the initial (positive or negative) buoyancy (limit positive buoy.
! threshold to XSBMIN fn(delta P)), ONLY for first 5 levels of lift

    IF (k  >=  ntml(index1(i)) .AND.                              &
        k  <=  ntml(index1(i)) + 4) THEN

      xsbmin_v(index1(i),k) = MIN ( xsbmin_v(index1(i),k),        &
      - 0.5 * ( th(index1(i),ntml(index1(i))+1) *                 &
      (1.0 + c_virtual * q(index1(i),ntml(index1(i)) + 1))        &
      - th(index1(i),ntml(index1(i))) * (1.0 +                    &
      c_virtual * q(index1(i),ntml(index1(i))))  ) + 0.5 )

    END IF

  END DO  !ncposs

! L_q_interact_if0:
  IF (l_q_interact) THEN
! ----------------------------------------------------------------------
!     Follow-on calculation from QCLP(*,K) and QCFP(*,K) initialization.
!       Duplicates code from Convection Parcel Lifting Scheme, LIFPAR.
! ----------------------------------------------------------------------
! ntml is level from which convection starts [dp/sh conv only]

    DO i=1, ncposs
      IF (k  ==  ntml(index1(i))) THEN
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
      uek_c2(i)    = uek_c(index2(i))
      uekp1_c2(i)  = uekp1_c(index2(i))
      vek_c2(i)    = vek_c(index2(i))
      vekp1_c2(i)  = vekp1_c(index2(i))
      dqsthkp1_c2(i) = dqsthkp1(index1(index2(i)))
      qsekp1_c2(i)   = qsekp1_c(index2(i))
      pstar_c2(i)    = pstar(index1(index2(i)))
      thpkp1_c2(i)   = thpkp1_c(index2(i))
      qpkp1_c2(i)    = qpkp1_c(index2(i))
      upkp1_c2(i)    = upkp1_c(index2(i))
      vpkp1_c2(i)    = vpkp1_c(index2(i))
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
      blowst_c2(i) = blowst(index1(index2(i)))
      l_shallow_c2(i) = .FALSE.
      l_mid_c2(i)     = .FALSE.
      cumulus_c2(i)   = .TRUE.
      ekp14_c2(i)  = ekp14_c(index2(i))
      ekp34_c2(i)  = ekp34_c(index2(i))
      amdetk_c2(i)    = amdetk(index1(index2(i)))
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
    END DO
    IF (l_q_interact) THEN
      DO i = 1,nconv
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
        start_lev3c(i) = ntml(index1(index2(i)))
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
    END DO
    IF (l_q_interact) THEN
      DO i = 1,nconv
!           PC2 variables
        qclpkp1_c2(i) = qclpkp1_c(index2(i))
        qcfpkp1_c2(i) = qcfpkp1_c(index2(i))
!           Compress increment fields
        dqclek_c2(i) = dqclbydt(index1(index2(i)),k)
        dqclekp1_c2(i) = dqclbydt(index1(index2(i)),k+1)
        dqcfek_c2(i) = dqcfbydt(index1(index2(i)),k)
        dqcfekp1_c2(i) = dqcfbydt(index1(index2(i)),k+1)
        dcflek_c2(i) = dcflbydt(index1(index2(i)),k)
        dcflekp1_c2(i) = dcflbydt(index1(index2(i)),k+1)
        dcffek_c2(i) = dcffbydt(index1(index2(i)),k)
        dcffekp1_c2(i) = dcffbydt(index1(index2(i)),k+1)
        dbcfek_c2(i) = dbcfbydt(index1(index2(i)),k)
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
    DO i = 1,nconv
      tcw_c2(i)    = tcw(index1(index2(i)))
      depth_c2(i)  = depth(index1(index2(i)))
      cclwp_c2(i)  = cclwp(index1(index2(i)))
      cape_c2(i)    = cape(index1(index2(i)))
      dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
      relh_c2(i)    = relh(index1(index2(i)))
      dptot_c2(i)   = dptot(index1(index2(i)))
      rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
      ! note rbuoy_p_old has been set at end of last pass through loop
      ! and needs to be compressed to be fed into convec2 and hence parcel
      rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
      the_here_c2(i)=the_here(index1(index2(i)))
      thp_here_c2(i)=thp_here(index1(index2(i)))
      qe_here_c2(i)=qe_here(index1(index2(i)))
      qp_here_c2(i)=qp_here(index1(index2(i)))
      thpixs_v_c2(i) = thpixs_v(index1(index2(i)),k)
      qpixs_v_c2(i)  = qpixs_v(index1(index2(i)),k)
      xsbmin_v_c2(i) = xsbmin_v(index1(index2(i)),k)
      iccb_c2(i)     = iccb(index1(index2(i)))
      icct_c2(i)     = icct(index1(index2(i)))
      cca_2d_c2(i)   = cca_2d(index1(index2(i)))
      lcca_c2(i)     = lcca(index1(index2(i)))
      lcbase_c2(i)   = lcbase(index1(index2(i)))
      lctop_c2(i)    = lctop(index1(index2(i)))
    END DO
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          trapk_c2(i,ktra) = trapk_c(index2(i),ktra)
          dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
          dtraekp1_c2(i,ktra) = dtrabydt(index1(index2(i)),k+1,ktra)
        END DO
      END DO
    END IF
  END IF  ! nconv>0

!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.

!      Subroutine CONVEC2

!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------


! DEPENDS ON: convec2_4a5a
  CALL convec2_4a5a(nconv,n_dp,nlev,ntra,k,dp_on, dp_sdet_on, dp_new_termc,&
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


  IF (nconv  >   0) THEN
    IF (flg_entr_up) THEN
      DO i = 1,nconv
        entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *    &
             (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)      &
           * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
        IF (bterm_c2(i)) THEN
          entrain_up(index1(index2(i)),k+1) = 0.0
        END IF
      END DO
    END IF

! Calculate fractional detrainment rate for level k
! (and k+1 if bterm=.T.)

    IF (flg_detr_up) THEN
      DO i = 1,nconv
        detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)            &
                 + deltak_c2(i) * (1.0 - amdetk_c2(i)))             &
                      * flx(index1(index2(i)),k)
        IF (bterm_c2(i)) THEN
          detrain_up(index1(index2(i)),k+1) =                       &
             -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
        END IF
      END DO
    END IF
  END IF

! Write CONVEC2 compressed output arrays back to full fields

  DO i = 1,n_dp
    !qp, thp already initialised outside main loop
    xpk(i,k+1)    = 0.0
    flx(i,k+1)    = 0.0
    depth(i)      = 0.0
    precip(i,k+1) = 0.0
    qclp(i,k+1)   = 0.0
    qcfp(i,k+1)   = 0.0
  END DO
  DO i = 1,n_dp
    bgmk(i)       = .FALSE.
    bterm(i)      = .FALSE.
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO i = 1,n_dp
        trap(i,k+1,ktra) = 0.0
      END DO
    END DO
  END IF

  IF (l_mom_gk) THEN
    DO i = 1,n_dp
      up(i,k+1) = 0.0
      vp(i,k+1) = 0.0
    END DO
  END IF

  IF (nconv  >   0)THEN
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
       max_cfl(index1(index2(i))) =                             &
                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
       relh(index1(index2(i)))  = relh_c2(i)
       dptot(index1(index2(i))) = dptot_c2(i)
     END DO

     IF (l_q_interact) THEN
       DO i = 1,nconv
!          PC2 variables
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
     IF (l_mom_gk) THEN ! Gregory Kershaw CMT

       DO i = 1,nconv
         up(index1(index2(i)),k+1) = upk_c2(i)
         vp(index1(index2(i)),k+1) = vpk_c2(i)
       END DO

       DO i = 1,nconv
         dubydt(index1(index2(i)),k) = duek_c2(i)
         dvbydt(index1(index2(i)),k) = dvek_c2(i)
         dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
         dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
         eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
         eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
       END DO
     END IF  ! L_mom_gk

     IF (l_mom) THEN    ! needed for all versions
       DO i = 1,nconv
         flxkp12(k,index1(index2(i)))  = flxkp12_c2(i)
       END DO
     END IF  ! L_mom

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
       DO i = 1 ,nconv
         up_flux(index1(index2(i)),k+1) = flxk_c2(i)
       END DO
     END IF
! Note that because of numbering system,
! p_layer_centres(k) which lies on level k+1/2 equates to
! z_rho(k+1), hence up_flux_half on k+1/2 will be numbered
! as k+1 in final output
     IF (flg_up_flx_half) THEN
       DO i = 1 ,nconv
         up_flux_half(index1(index2(i)),k+1) = flxkp12_c2(i)
       END DO
     END IF
 
  END IF

!   Write adaptive diagnostics for this level to full array for output

  DO i = 1,n_dp
    rbuoy_p_out(i,k) = rbuoy_p_here(i)
    the_out(i,k) = the_here(i)
    thp_out(i,k) = thp_here(i)
    qe_out(i,k) = qe_here(i)
    qp_out(i,k) = qp_here(i)
  END DO

!   Write rbuoy for this level to rbuoy_p_old for previous level
!   for use in parcel

  DO i = 1,n_dp
    rbuoy_p_old(i) = rbuoy_p_here(i)
  END DO

!-----------------------------------------------------------------------
! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
!      removed by convection over timescale cape_timescale.
!-----------------------------------------------------------------------

! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_dp) with index_nterm


  nterm = 0
  DO i = 1,n_dp
    IF (bterm(i)) THEN
      nterm = nterm + 1
      index_nterm(nterm) = i
    END IF
  END DO

! Note nterm is often small so vectorisation of the loop over nterm
! will not always result in faster code.

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
          IF (cape_ts_new(i) < cape_timescale) THEN
            ind_cape_reduced(i) = 1.0
          END IF
          cape_ts_used(i) = cape_ts_new(i)
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
          cape_ts_new(i) =                                              &
                   MIN(MAX(cape_timescale*(1.0-relh(i)/dptot(i))/0.4    &
                  ,timestep)                                            &
                  ,cape_timescale)

          IF (cape_ts_new(i) < cape_timescale) THEN
            ind_cape_reduced(i) = 1.0
          END IF
          cape_ts_used(i) = cape_ts_new(i)
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
          cape_ts_used(i) = cape_ts_new(i)
        END IF
      END DO

! w based cape closure; if w_cape_limit > 1000. reverts to cape_timescale
    CASE(3)

      ! Initialise araay to cape timescale and then alter as required.
      cape_ts_new(:) =  cape_timescale

      IF ( w_cape_limit < 1000.0 ) THEN
      !  This section includes test on w_max

!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            ! new denominator introduced at vn6.6
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) =   cape_ts_w * w_cape_limit/           &
                      (w_cape_limit+ (w_max(i)-w_cape_limit)*wcape_fac)
              ! set indicator that CAPE reduced
              ind_cape_reduced(i) = 1.0
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF !  w_max(i) > w_cape_limit
            cape_ts_used(i) = cape_ts_new(i)
          END IF  ! dcpbydt > 0
        END DO
      END IF  ! w_cape_limit

! Option 4 - a w based cape closure; Grid-box area scaled CAPE closure
    CASE(4)

      IF ( w_cape_limit < 1000.0 ) THEN
      ! This section includes test on w_max

!CDIR NODEP
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
            ELSE
              cape_ts_new(i) = cape_timescale * cape(i) / cape_min +    &
                               cape_timescale * EXP(-cape(i) / cape_min)
            END IF !  w_max(i) > w_cape_limit
            cape_ts_used(i) = cape_ts_new(i)
          END IF  ! dcpbydt > 0
        END DO
      ELSE
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN

            cape_ts_new(i) = cape_timescale * cape(i) / cape_min +      &
                             cape_timescale * EXP( - cape(i) /cape_min)

            cape_ts_used(i) = cape_ts_new(i)
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
              ind_cape_reduced(i) = 1.0
            ELSE
              IF ( relh(i) / dptot(i) >= 0.75 ) THEN
                cape_ts_new(i) = cape_timescale *                       &
                                 ( 0.2373 / (relh(i) / dptot(i))**5)
                ind_cape_reduced(i) = 1.0
              ELSE
                cape_ts_new(i) = cape_timescale
              END IF ! relh(i) / dptot(i) >= 0.75
            END IF !  w_max(i) > w_cape_limit
            cape_ts_used(i) = cape_ts_new(i)
          END IF  ! dcpbydt > 0
        END DO
      ELSE
        DO j = 1, nterm
          i = index_nterm(j)
          IF ( dcpbydt(i) > 0.0 ) THEN
            IF ( relh(i) / dptot(i) >= 0.75 ) THEN
              cape_ts_new(i) = cape_timescale *                          &
                                  ( 0.2373 / (relh(i) / dptot(i))**5)
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF ! relh(i) / dptot(i) >= 0.75
            cape_ts_used(i) = cape_ts_new(i)
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
              ind_cape_reduced(i) =1.0
            ELSE
              cape_ts_new(i) = cape_timescale
            END IF
            ! Further reduction if w_max above critical value
            IF ( w_max(i) > w_cape_limit ) THEN
              cape_ts_new(i) = cape_ts_new(i) * w_cape_limit/         &
                       (w_cape_limit + (w_max(i)-w_cape_limit)*wcape_fac)
              ind_cape_reduced(i) =1.0
            END IF
            ! Limit CAPE timescale to convective timestep
            cape_ts_new(i) = MAX(cape_ts_new(i), timestep)
            cape_ts_used(i) = cape_ts_new(i)

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

      ! If convection has terminated write cape to diagnostic output
      ! variable (cape_out).

      cape_out(i) = cape(i)
      cape(i)     = 0.0
      dcpbydt(i)  = 0.0

    END DO


! Work out scaled mass flux needed to keep cfl ratio below limit.
! Note CAPE closure assumed

    DO j = 1,nterm
      i = index_nterm(j)
      max_cfl(i) = max_cfl(i) * timestep

      IF (max_cfl(i)  >   cfl_limit) THEN
        flx_init_new(i) = flx_init_new(i) * cfl_limit             &
                                       / max_cfl(i)
        cfl_limited(i) = 1.0       ! deep convection is CFL limited
      ELSE
        flx_init_new(i) = flx_init_new(i)        
      END IF

      IF (flx_init_new(i)  >   flxmax_init(i)) THEN
        flx_init_new(i) = flxmax_init(i)
      END IF
      max_cfl(i) = 0.0
    END DO      ! j (nterm)


    !
    ! Scale cloud fraction
    !
    IF (l_ccrad) THEN

      DO j = 1,nterm
        i = index_nterm(j)

        IF (flx_init_new(i) > 0.0) THEN
          scale_f(i) = flx_init_new(i) / flx_init(i)
          cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

          ! set the flx_init to the new value to provide the real initial mass
          ! flux under all conditions
          flx_init(i) = flx_init_new(i)

        END IF

        ! Check scaled cloud fraction not smaller than minimum value
        ! (2.0E-5) or greater than unity.
        !
        ! (Was moved out of scaling if test to ensure these limits
        ! at all times, not just when cca_2d is scaled)

        cca_2d(i) = MAX(2.0e-5, cca_2d(i))
        cca_2d(i) = MIN(1.0e+0, cca_2d(i))

      END DO      ! j (nterm)

    ELSE      ! original

      DO j = 1,nterm
        i = index_nterm(j)
        IF (flx_init_new(i) > 0.0) THEN
          scale_f(i) = flx_init_new(i) / flx_init(i)
          cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

          !
          ! Check scaled cloud fraction not smaller than minimum value
          ! (2.0E-5) or greater than unity.
          !

          cca_2d(i) = MAX(2.0e-5, cca_2d(i))

          IF (cca_2d(i) > 1.0) THEN
            cca_2d(i) = 1.0
          END IF

          ! set the flx_init to the new value to provide the real initial mass
          ! flux under all conditions
          flx_init(i) = flx_init_new(i)

        END IF
      END DO      ! j (nterm)

    END IF      ! l_ccrad



! Carry out cape and cfl scaling


    DO kt = 2, k+1
      DO j = 1,nterm
        i = index_nterm(j)
        IF (kt  >=  ntml(i) .AND. flx_init_new(i) >   0.0) THEN

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
          END IF
          IF (l_mom) THEN     ! required for all versions
            IF (kt <  k+1) THEN
              flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
            END IF
          END IF
          IF (l_tracer) THEN
            DO ktra = 1,ntra
             dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            END DO
          END IF

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
          IF (flg_entr_up) THEN
            entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)
          END IF
          IF (flg_detr_up) THEN
            detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
          END IF


        END IF !kt>ntml and flx_init_new>0
      END DO  ! j loop
    END DO  ! kt loop

! Set cape scaling parameters for next level

    DO j = 1,nterm
      i = index_nterm(j)
      det_lev(i)= k+1
    END DO  ! nterm loop

!-----------------------------------------------------------------------
! 3.5  Downdraft calculation - on all points where convection is
!      terminating.

!      Subroutine DD_CALL

!      UM Documentation Paper 27, part 2

!-----------------------------------------------------------------------

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
! Set logical for use later

      IF (deltap_cld  >   15000.0 .AND. bgmk(i) .AND.               &
                                tempnum  >   1e-12) THEN
        b_dd(i) = .TRUE.
      ELSE
        b_nodd(i) = .TRUE.
      END IF
    END DO  ! nterm

! If convection has terminated write cape to diagnostic output
! variable (cape_out).
! Set kterm array which holds the level index for termination
! of convection.


    DO j = 1,nterm
      i=index_nterm(j)
      IF (k >= ntml(i)) THEN
        kterm(i) = k
      END IF
      bconv(i) = .FALSE.
    END DO
  END IF ! nterm > 0

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

END DO



!-----------------------------------------------------------------------
! 4.0 Down draughts  - now 2 options
!                      original code
!                      Emanuel down draught code
! Note the level at which deep convection terminates has been stored
! in the above updraught loop as Kterm.
!-----------------------------------------------------------------------

IF (l_eman_dd) THEN


! Work out maximum termination level
  kmax_term = 2
  DO i = 1,n_dp
    IF(kterm(i) >  kmax_term) THEN
      kmax_term = kterm(i)
    END IF
  END DO

! DEPENDS ON: eman_dd
  CALL eman_dd (n_dp,kmax_term,nlev,trlev,ntra,kterm,l_tracer     &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries        &
,                      timestep                                   &
,                      th, q, qse, tracer, precip                 &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow ,dwn_flux                       &
                     )

ELSE         ! original down draught code
!-----------------------------------------------------------------------
! 4.1  Downdraft calculation - on all points where convection is
!      terminating.

!      Subroutine DD_ALL_CALL

!      UM Documentation Paper 27, part 2

!-----------------------------------------------------------------------

  npossdd = 0
  DO i = 1,n_dp
    IF (b_dd(i)) THEN
      npossdd = npossdd +1
      index_possdd(npossdd) = i
    END IF
  END DO

  IF (npossdd  >   0) THEN

! Work out maximum termination level
    kmax_term = 2
    DO i = 1,npossdd
      IF(kterm(index_possdd(i)) >  kmax_term) THEN
        kmax_term = kterm(index_possdd(i))
      END IF
    END DO

! DEPENDS ON: dd_all_call_4a5a
    CALL dd_all_call_4a5a (n_dp,npossdd,kmax_term,nlev,trlev,ntra &
,                      kterm, iccb, icct, index_possdd            &
,                      l_tracer, bwater(1,2)                      &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      recip_pstar,timestep , cca_2d              &
,                      thp, qp, th, q, qse, trap,tracer, flx,precip &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow, rain_3d, snow_3d, dwn_flux     &
,                      entrain_dwn, detrain_dwn)


  END IF

!-----------------------------------------------------------------------
! 4.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
  nnodd = 0
  DO i = 1,n_dp
    IF (b_nodd(i)) THEN
      nnodd = nnodd +1
      index_nodd(nnodd) = i
    END IF
  END DO

  IF (nnodd  >   0) THEN

! Work out maximum termination level
    kmax_term = 2
    DO i = 1,nnodd
      IF(kterm(index_nodd(i)) >  kmax_term) THEN
        kmax_term = kterm(index_nodd(i))
      END IF
    END DO
! Only add 1 if kmax_term is less than model levels
! (which should be true).
    IF (kmax_term  <  nlev ) THEN
      kmax_term = kmax_term + 1
    END IF

! Surface precipitation calculation

! DEPENDS ON: evap_bcb_nodd_all
    CALL evap_bcb_nodd_all(n_dp,nnodd,kmax_term,kterm             &
,                      iccb, index_nodd, bwater(1,2)              &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      timestep , cca_2d, th, q, qse, precip      &
,                      dthbydt, dqbydt, rain, snow                &
,                      rain_3d, snow_3d)
  END IF


END IF        ! test on down draught type


!-----------------------------------------------------------------------
! 4.3 Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!-----------------------------------------------------------------------

IF (.NOT. l_ccrad) THEN
  DO i=1, n_dp

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

!-----------------------------------------------------------------------
! 5.0  Convective Momentum Transport (if L_mom = .true.)
!------------------------------ -----------------------------------------

IF (l_mom) THEN
  SELECT CASE (deep_cmt_opt)
  CASE (2)       ! Gregory-Kershaw deep CMT

    ! Do nothing here as calculated in above parcel ascent

  CASE DEFAULT   ! (0/1/5) Alan Grants Eddy viscosity CMT

    IF (l4a_kterm) THEN

    ! altered to use kterm instead of ntpar
      DO i = 1,n_dp
        IF (kterm(i) >=  nlcl_uv(i)) THEN
          ntop_uv(i)    = kterm(i) + 1
        ELSE     ! case where deep convection fails
    ! I think in this case the cloud base mass flux will be zero so
    ! there will be no CMT. (The value will not matter)
          ntop_uv(i)    = ntpar(i) + 1
        END IF
        ptop_uv(i)    = phalf_uv(ntop_uv(i),i)
      END DO
    ELSE          ! original code
      DO i = 1,n_dp
        ntop_uv(i)    = ntpar(i) + 1
        ptop_uv(i)    = phalf_uv(ntop_uv(i),i)
      END DO
    END IF       ! test on l4a_kterm

    nterm = 0

! Set cloud base mass flux equal to mass flux at half level below
! the LCL.

    DO i = 1,n_dp
      IF (kterm(i)  >=  nlcl_uv(i) .AND. kterm(i) < nlev-1) THEN
        nterm = nterm + 1
        cu_term(nterm) = i
        cu_tend(nterm) = i
        mb(i) = flxkp12(nlcl_uv(i),i)
        DO j = 1,nlev
          flxkp12(j,i) = 0.0
        END DO
      ELSE IF (kterm(i) == nlev-1) THEN
      ! Problem deep convection has gone to the top of the model
      ! Return location of deep problem point

        IF (printstatus >= prstatus_normal) THEN
          write(6,*) 'Deep point ',i,' kterm ',kterm(i),' nlcl ',nlcl_uv(i) 
          write(6,*) ' theta ',(th(i,k),k=1,nlev)
          write(6,*) ' q ',(q(i,k),k=1,nlev)
          write(6,*) ' qcl ',(qcl(i,k),k=1,nlev)
          write(6,*) ' qcf ',(qcf(i,k),k=1,nlev)

        END IF
        error_point = i   
        IF (lhook) CALL dr_hook('DEEP_CONV_4A',zhook_out,zhook_handle)
        RETURN        

      END IF

! initialise output arrays as lower level subroutines don't set all
! values

      DO j = 1,nlev
        uw(j,i)=0.0
        vw(j,i)=0.0
        visc(j,i)=0.0
      END DO

    END DO

    IF (nterm  >   0) THEN
! DEPENDS ON: cmt_mass
      CALL cmt_mass(n_dp, n_dp, nlev, nterm, cu_term,                 &
              kterm, cu_tend, n_0degc, nlcl_uv, ntop_uv,              &
              mb, p_0degc_uv, plcl_uv, ptop_uv, phalf_uv, p_uv,       &
            ! Output arguments
              flxkp12 ,mass_dwn, visc)

! DEPENDS ON: deep_grad_stress
      CALL deep_grad_stress(n_dp,n_dp,n_dp,nlev,nlcl_uv,ntop_uv,         &
                            nterm,cu_term,cu_tend,                       &
                            ue_p,ve_p,visc,phalf_uv,p_uv,rho_uv,timestep,&
                            ! Output
                            uw,vw)

! DEPENDS ON: deep_ngrad_stress
      CALL deep_ngrad_stress(n_dp,n_dp,n_dp,nterm,nlev,               &
                             nlcl_uv,ntop_uv,cu_term,cu_tend,cu_tend, &
                             pstar,uw0,vw0,zlcl_uv,ue_p,ve_p,visc,    &
                             flxkp12,p_uv,phalf_uv,rho_uv,timestep,   &
                             ! Input/output
                             uw,vw,                                   &
                             ! Output
                             uw_base,vw_base,uw_deep,vw_deep)

! DEPENDS ON: deep_cmt_incr
      CALL deep_cmt_incr(n_dp,n_dp,n_dp,nlev,nterm,                  &
                         nlcl_uv,ntop_uv,cu_term,cu_tend,            &
                         zlcl_uv,phalf_uv,p_uv,rho_uv,               &
                         uw_base,vw_base,uw,vw,                      &
                         ! Output
                         dubydt,dvbydt)

    END IF  ! nterm > 0

  CASE (3,4)       ! new turbulence scheme

    nterm = 0   ! count of number of deep points which actually convected
    DO i = 1,n_dp
      IF (kterm(i)  >=  nlcl_uv(i)) THEN
        nterm = nterm + 1
        cu_term(nterm) = i

! Use CAPE scaled mass flux as initial mass flux rather than CRM derived value
! to be consistent with thermodynamic part of convection.
        mb(i) = flxkp12(nlcl_uv(i),i)
        zlcl(i) = z_rho(i,ntml(i))

! Cloud velocity scale - derived from CRM simulations
!             wcld = (C_mass*wstar*CAPE)**(1/3)

        wcld(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i)) &
              * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
      END IF
    END DO

    IF (nterm > 0) THEN
! DEPENDS ON: deep_turb_cmt
      CALL deep_turb_cmt (n_dp, nterm, nlev,deep_cmt_opt,          &
                    ntml, kterm,cu_term,freeze_lev,                &
                    timestep,                                      &
                    uw0, vw0, mb, wcld, wstar ,zlcl,               &
                    flx,                                           &
                    r_rho, r_theta, z_rho, z_theta,                &
                    rho, rho_theta,                                &
                    r2rho, r2rho_th, dr_across_th, dr_across_rh,   &
                    u, v,                                          &
                    dubydt, dvbydt, uw_deep, vw_deep)

    END IF           ! nterm > 0

  END SELECT       ! deep_cmt_opt

END IF  ! L_mom

!-----------------------------------------------------------------------
! 6.0  Mix the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.
!-----------------------------------------------------------------------
IF (bl_cnv_mix == 1) THEN
  DO i = 1,n_dp
    index1(i) = i
  END DO
! DEPENDS ON: mix_ipert_4a5a
  CALL mix_ipert_4a5a(n_dp, n_dp,nlev, nbl, ntml, index1,             &
                      p_layer_boundaries,                             &
                      exner_layer_centres, dthbydt, dqbydt, flx_init, &
                      thpert, qpert)
END IF

!-----------------------------------------------------------------------
! 7.0  CCRad - Calculate CCA for Deep levels only
!-----------------------------------------------------------------------

IF (l_ccrad) THEN

  !-----------------------------------------------------------------
  ! 7.1 Calculate CCA_2D of Deep Cloud
  !-----------------------------------------------------------------
  SELECT CASE (cca2d_dp_opt)
    CASE(srf_precip)
      DO i=1, n_dp
        IF (iccb(i) /= 0) THEN ! Deep convection was successful

          ! Use determination of CCA_2D based on surface precip. rate
          ! from deep cloud.

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
            ! NOTE: Downdraughts & Evaporation still being fed cca_2d
            !       derived from TCW, This issue may require further
            !       investigation.
            lcca(i) = cca_2d(i)

          END IF

        END IF      ! iccb
      END DO      ! i (n_dp)

    CASE(total_condensed_water)
      ! cca_2d_dp left unchanged from code, which is based on
      ! TCW (Total Condensed Water) (This is a rate)

  END SELECT


  IF (.NOT. l_dcpl_cld4pc2) THEN
    !-------------------------------------------------------------------
    ! Apply Deep CCA/CCW Tuning factors
    !-------------------------------------------------------------------
    ! Because of decoupling move up a level so scaling is only applied
    ! to section 0 variables, section 5 variables continue to be
    ! scaled here to maintain bit-comparison, would ideally be moved
    ! up to Glue_Conv in future
    DO i=1, n_dp
      cca_2d(i) = cca_2d(i)*cca_dp_knob

      ! Make sure cca_2d_dp within limits
      cca_2d(i) = MAX(0.0,    cca_2d(i))
      cca_2d(i) = MIN(1.0e+0, cca_2d(i))
    END DO

    DO k=1, nlev
      DO i=1, n_dp
        ccw(i,k) = ccw(i,k)*ccw_dp_knob
      END DO
    END DO

  END IF


  !---------------------------------------------------------------------
  ! 7.3 Apply CCA_2D to 3d cloud profile
  !---------------------------------------------------------------------

  IF (l_anvil) THEN
! DEPENDS ON: CALC_3D_CCA
    CALL calc_3d_cca                                                    &
       ( n_dp, n_dp, nlev, n_cca_lev, nbl, iccb, icct                   &
       , p_layer_boundaries, freeze_lev, cca_2d, cca, z_theta, z_rho )

    ! NOTE: iccb, icct are layer centres (theta levels) at this
    !       point.

  ELSE

    ! Apply cca_2d to all levels from deep base to deep top
    DO k=1, n_cca_lev
      DO i=1, n_dp
        IF ( k >= iccb(i) .AND. &
             k <= icct(i) ) THEN
          cca(i,k) = cca_2d(i)
        END IF
      END DO
    END DO

  END IF ! l_anvil

END IF ! l_ccrad
!-----------------------------------------------------------------------
! 8.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEEP_CONV_4A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE deep_conv_4a
