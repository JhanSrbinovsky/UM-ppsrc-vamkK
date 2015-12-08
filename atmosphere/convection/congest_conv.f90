! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Congestus convection scheme
!

SUBROUTINE congest_conv(nbl,nlev,ntra,n_cca_lev,n_cg,trlev,       &
                       bland,delthvu,exner_layer_centres,         &
                       exner_layer_boundaries,                    &
                       l_calc_dxek,l_q_interact,                  &
                       l_tracer,ntml,ntpar,                       &
                       pstar,p_layer_centres,                     &
                       p_layer_boundaries,z_theta,z_rho,          &
                       r_theta,r_rho,                             &
                       r2rho_th, r2rho, dr_across_th,dr_across_rh,&
                       q,q1_sd,t1_sd,th,timestep,u,v,uw0,vw0,     &
                       wstar,wthvs,entrain_coef,                  &
                       zlcl_uv,ztop_uv,freeze_lev,                &
                       recip_pstar,qse,                           &
                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
                       qcl,tracer,cape_out,cclwp,ccw,cca,         &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                       dqclbydt,dthbydt,                          &
                       dubydt,dvbydt,dtrabydt,                    &
                       detrain_up,detrain_dwn,                    &
                       entrain_up,entrain_dwn,                    &
                       iccb,icct,lcca,                            &
                       lcbase,lctop,rain,snow,                    &
                       rain_3d,snow_3d, up_flux, up_flux_half,    &
                       dwn_flux,uw_shall,vw_shall,kterm,          &
                       tcw,cca_2d,                                &
                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out)

!-----------------------------------------------------------------------
!
! Purpose:
! Congestus convection scheme - works on points diagnosed as congestus
! in  subroutine CONV_DIAG. (At present based on shallow routine.)
!
!   Called by GLUE_CONV.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP 3 v8.3 programming standards
!-----------------------------------------------------------------------

USE atmos_constants_mod, ONLY:                                    &
    r, cp , pref, kappa, repsilon, c_virtual

USE cv_run_mod, ONLY:                                             &
    l_mom, bl_cnv_mix, icvdiag, l_cv_conserve_check,              &
    cca2d_sh_opt, cca_sh_knob, ccw_sh_knob, cnv_wat_load_opt,     &
    l_ccrad

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, grant_lock,                            &
    thpixs_shallow, qpixs_shallow, c_mass

USE cv_dependent_switch_mod, ONLY:                                &
    cg_on, mdet_cg_on, cg_ent_on, cg_sdet_on, cg_new_termc

USE cv_stash_flg_mod, ONLY:                                               &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx, flg_entr_up, flg_entr_dwn,  &
    flg_detr_up, flg_detr_dwn, flg_uw_shall, flg_vw_shall

USE earth_constants_mod, ONLY: g

USE water_constants_mod, ONLY: lc, lf, tm

USE yomhook, ONLY  : lhook, dr_hook
USE parkind1, ONLY : jprb, jpim
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
 ,n_cg                 & ! No. of congestus convection points
 ,trlev                  ! No. of model levels on which tracers are included

LOGICAL, INTENT(IN) :: bland(n_cg) ! Land/sea mask

REAL, INTENT(IN)    :: delthvu(n_cg) !Integral of undilute parcel
                                     ! buoyancy over convective cloud
                                     ! layer (Kelvin m)

REAL, INTENT(IN)    ::                &
  exner_layer_centres(n_cg,0:nlev)    & ! Exner
 ,exner_layer_boundaries(n_cg,0:nlev)   ! Exner at half level above
                                        ! exner_layer_centres

LOGICAL, INTENT(IN) :: &
  l_calc_dxek          & ! Switch for calculation of condensate increment
 ,l_q_interact         & ! Switch allows overwriting parcel variables
                         ! when calculating condensate incr.
 ,l_tracer               ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: &
  ntml(n_cg)           & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(n_cg)            ! Top level of initial parcel ascent in BL scheme
                         ! defined relative to theta, q grid

REAL, INTENT(IN)    ::            &
  pstar(n_cg)                     & ! Surface pressure (Pa)
 ,p_layer_centres(n_cg,0:nlev)    & ! Pressure (Pa)
 ,p_layer_boundaries(n_cg,0:nlev)   ! Pressure at half level above
                                    ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

REAL, INTENT(IN)    ::    &
  z_theta(n_cg,nlev)      & ! height of theta levels (m)
, z_rho(n_cg,nlev)        & ! height of rho levels (m)
, r_theta(n_cg,0:nlev)    & ! radius of theta levels (m)
, r_rho(n_cg,nlev)        & ! radius of rho levels (m)
, r2rho_th(n_cg,nlev)     & ! r2*rho theta levels (kg/m)
, r2rho(n_cg,nlev)        & ! r2*rho rho levels (kg/m)
, dr_across_th(n_cg,nlev) & ! dr across theta levels (m)
, dr_across_rh(n_cg,nlev)   ! dr across rho levels (m)

REAL, INTENT(IN)  ::  &
  q(n_cg,nlev)        & ! Model mixing ratio (kg/kg)
 ,q1_sd(n_cg)         & ! Standard deviation of turbulent flucts. of layer 1 q
                        ! (kg/kg)
 ,t1_sd(n_cg)         & ! Standard deviation of turbulent flucts. of layer 1
                        ! temp. (K)
 ,th(n_cg,nlev)       & ! Model potential temperature (K)
 ,timestep            & ! Model timestep (s)
 ,u(n_cg,nlev)        & ! Model U field (m/s)
 ,v(n_cg,nlev)        & ! Model V field (m/s)
 ,uw0(n_cg)           & ! U-comp of surface stress (N/m2)
 ,vw0(n_cg)           & ! V-comp of surface stress (N/m2)
 ,wstar(n_cg)         & ! Convective velocity scale (m/s)
 ,wthvs(n_cg)         & ! Surface flux of THV (Pa m/s2)
 ,entrain_coef(n_cg)  & ! entrianment coefficients
 ,zlcl_uv(n_cg)       & ! Lifting condensation level defined for the uv grid (m)
 ,ztop_uv(n_cg)         ! Top of cloud layer defined for the uv grid (m)

INTEGER, INTENT(IN) :: freeze_lev(n_cg) ! Level index for freezing level

REAL, INTENT(IN) :: recip_pstar(n_cg)  ! Reciprocal of pstar array

REAL, INTENT(IN) :: qse(n_cg,nlev) ! Saturation mixing ratio of
                                   ! cloud environment (kg/kg)

! Arguments with intent INOUT:

REAL, INTENT(INOUT) ::   &
  bulk_cf(n_cg,nlev)     & ! Bulk total cloud volume ( )
 ,cf_frozen(n_cg,nlev)   & ! Frozen water cloud volume ( )
 ,cf_liquid(n_cg,nlev)   & ! Liq water cloud volume ( )
 ,qcf(n_cg,nlev)         & ! Ice condensate mix ratio (kg/kg)
 ,qcl(n_cg,nlev)         & ! Liq condensate mix ratio (kg/kg)
 ,tracer(n_cg,trlev,ntra)  ! Model tracer fields(kg/kg)

! Arguments with intent OUT:

REAL, INTENT(OUT) :: &
  cape_out(n_cg)     & ! Saved convective available potential energy for
                       ! diagnostic output (J/kg)
 ,cclwp(n_cg)        & ! Condensed water path (k/m2)
 ,ccw(n_cg,nlev)     & ! Convective cloud liquid water on model levels (g/kg)
 ,cca(n_cg,nlev)       ! Convective cloud amount on model levels (fraction)

REAL, INTENT(OUT) ::  &
  dbcfbydt(n_cg,nlev) & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_cg,nlev) & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_cg,nlev) & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_cg,nlev)   & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_cg,nlev) & ! Increments to ice condensate due to conv (kg/kg/s)
 ,dqclbydt(n_cg,nlev) & ! Increments to liq condensate due to conv (kg/kg/s)
 ,dthbydt(n_cg,nlev)  & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_cg,nlev+1) & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_cg,nlev+1)   ! Increments to V due to CMT (m/s2)

REAL, INTENT(OUT) ::       &
  dtrabydt(n_cg,nlev,ntra)   !Increment to tracer due to convection (kg/kg/s)

REAL, INTENT(OUT) ::     &
  detrain_up(n_cg,nlev)  & ! Fractional detrainment rate into updraughts (Pa/s)
 ,detrain_dwn(n_cg,nlev) & ! Fractional detrainment rate into 
                           ! downdraughts (Pa/s)
 ,entrain_up(n_cg,nlev)  & ! Fractional entrainment rate into updraughts (Pa/s)
 ,entrain_dwn(n_cg,nlev)   ! Fractional entrainment rate into
                           ! downdraughts (Pa/s)

INTEGER, INTENT(OUT) :: iccb(n_cg) ! Convective cloud base level 

INTEGER, INTENT(OUT) :: icct(n_cg) ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(n_cg) ! Lowest conv. cloud amt. (%)


INTEGER, INTENT(OUT) :: lcbase(n_cg) ! Lowest conv. cloud base level

INTEGER, INTENT(OUT) :: lctop(n_cg) ! Lowest conv. cloud top level 

REAL, INTENT(OUT) ::   &
  rain(n_cg)           & ! Surface convective rainfall (kg/m2/s)
 ,snow(n_cg)           & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(n_cg,nlev)   & ! rainfall flux (kg/m2/s)
 ,snow_3d(n_cg,nlev)     ! snowfall flux (kg/m2/s)

REAL, INTENT(OUT) ::     &
  up_flux(n_cg,nlev)     & ! Updraught mass flux (Pa/s)
 ,up_flux_half(n_cg,nlev)& ! Updraught mass flux (Pa/s)
 ,dwn_flux(n_cg,nlev)    & ! Downdraught mass flux (Pa/s)
 ,uw_shall(n_cg,nlev)    & ! X-comp. of stress from shallow convection (kg/m/s2)
 ,vw_shall(n_cg,nlev)      ! Y-comp. of stress from shallow convection (kg/m/s2)

INTEGER, INTENT(OUT) :: kterm(n_cg) ! termination level

REAL, INTENT(OUT) ::   &
  tcw(n_cg)            & ! Total condensed water(kg/m2/s)
 ,cca_2d(n_cg)           ! 2D convective cloud amount (%)

! Adaptive detrainment output variables

REAL, INTENT(OUT) ::      &
  rbuoy_p_out(n_cg,nlev)  & !buoyancy excess
 ,the_out(n_cg,nlev)      & !th_E in parcel routine
 ,thp_out(n_cg,nlev)      & !th_P in parcel routine
 ,qe_out(n_cg,nlev)       & !q_E in parcel routine
 ,qp_out(n_cg,nlev)         !q_P in parcel routine

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Adaptive detrainment output variables

REAL ::               &
  rbuoy_p_here(n_cg)  & !buoyancy excess
 ,the_here(n_cg)      & !th_E in parcel routine
 ,thp_here(n_cg)      & !th_P in parcel routine
 ,qe_here(n_cg)       & !q_E in parcel routine
 ,qp_here(n_cg)       & !q_P in parcel routine
 ,rbuoy_p_old(n_cg)     !buoyancy excess on previous level

REAL :: zk(n_cg)               !heights for use in calc
REAL :: zkp12(n_cg)            !of moist static energy
REAL :: zkp1(n_cg)


INTEGER :: index1(n_cg),index2(n_cg)

INTEGER :: ncposs               ! No. of points which may convect

INTEGER :: nconv                ! No. of convecting points

REAL :: amdetk(n_cg)            ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

REAL :: b_calc                  ! Coefficient in thpert calc.

REAL :: c_calc                  ! Coefficient in thpert calc.

REAL :: cape(n_cg)              ! Convective available potential
                                ! energy (J/kg)

REAL :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

REAL :: dcpbydt(n_cg)           ! Rate of change of cape (J/kg/s)

REAL :: depth(n_cg)             ! Depth of convective cloud (m)

REAL :: delexkp1(n_cg)          ! Difference in exner ratio
                                ! across layer k+1

REAL :: dqsthk(n_cg)            ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k (kg/kg/K)

REAL :: dqsthkp1(n_cg)          ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k+1 (kg/kg/K)


REAL ::        &
  ekp14(n_cg)  & ! Entrainment coefficients at level k+1/4 multiplied by
                 ! appropriate layer thickness (dimensionless)
 ,ekp34(n_cg)  & ! Entrainment coefficients at level k+3/4 multiplied by
                 ! appropriate layer thickness (dimensionless)
 ,ekm14(n_cg)    ! Entrainment coefficients at level k-1+1/4 multiplied by
                 ! appropriate layer thickness (dimensionless)

REAL ::              &
  exk(n_cg)          & ! Exner ratio at layer k
 ,exkp1(n_cg)        & ! Exner ratio at layer k+1
 ,flxmax(n_cg)       & ! Maximum initial convective mass flux (Pa/s)
 ,flx_init(n_cg)     & ! Initial mass flux at cloud base (Pa/s)
 ,flx_init_new(n_cg) & ! flx_init scaled (Pa/s)
 ,flxmax_init(n_cg)  & ! Maximum possible initial mass flux (limited to the
                       ! mass in the initial convecting layer in Pa/s)
 ,max_cfl(n_cg)      & ! Max cfl ratio over a convecting layer
 ,p_lcl(n_cg)        & ! Pressure at LCL (Pa)
 ,precip(n_cg,nlev)    ! Amount of precip from each layer (kg/m/s)

REAL ::              &
  pk(n_cg)           & ! Pressure at midpoint of layer k (Pa)
 ,pkp1(n_cg)         & ! Pressure at midpoint of layer k+1 (Pa)
 ,delpk(n_cg)        & ! Pressure difference over layer k (Pa) 
 ,delpkp1(n_cg)      & ! Pressure difference over layer k+1 (Pa)
 ,delpkp12(n_cg)     & ! Pressure difference between layers k and k+1 (Pa)
 ,delp_uv_k(n_cg)    & ! Pressure difference across uv layer k (Pa)
 ,delp_uv_kp1(n_cg)    ! Pressure difference across uv layer k+1 (Pa)

REAL ::              &
  q_lcl(n_cg)        & ! Mixing ratio at LCL (kg/kg)
 ,qse_lcl(n_cg)      & ! Saturated q at LCL (kg/kg)
 ,rhum(n_cg)         & ! Dummy relative humidity (only used on shallow points)
 ,t_lcl(n_cg)        & ! Temperature at LCL (K)
 ,th_lcl(n_cg)       & ! Theta at LCL (K)
 ,thv_pert(n_cg)     & ! Theta_v parcel pertubation (K)
 ,thpert(n_cg)       & ! Theta parcel pertubation (K)
 ,qpert(n_cg)          ! q parcel pertubation (kg/kg)

INTEGER :: start_lev3c(n_cg)    ! Compressed convection
                                ! initiation level

REAL :: wsc(n_cg)               ! Convective velocity scale (m/s)

REAL :: wsc_o_mb(n_cg)          ! Convective velocity scale /mb

LOGICAL ::           &
  l_shallow(n_cg)    & ! Dummy variable (=.T.)
 ,l_mid(n_cg)        & ! Dummy variable (=.F.)
 ,cumulus(n_cg)        ! Dummy variable (=.T.)

LOGICAL ::           &
  bgmk(n_cg)         & ! Mask for points where parcel in layer k is saturated
 ,bwater(n_cg,2:nlev)& ! Mask for points at which condensate is liquid
 ,bwk(n_cg)          & ! mask for liquid condensate on k
 ,bwkp1(n_cg)        & ! mask for liquid condensate on k+1
 ,blowst(n_cg)       & ! Dummy variable indicating low enough stability for
                       ! convection to occur
 ,bterm(n_cg)        & ! Mask for points which have stopped convecting
 ,bconv(n_cg)        & ! Mask for points at which convection is occurring
 ,bcposs(n_cg)         ! Mask for points passing initial stability test

! Parcel variables

REAL :: qpi(n_cg)               ! Initial parcel mixing ratio (kg/kg)

REAL :: qp(n_cg,nlev)           ! Parcel mixing ratio (kg/kg)

REAL :: thpi(n_cg)              ! Initial parcel potential temp.(K)

REAL :: thp(n_cg,nlev)          ! Parcel potential temp (K)

REAL :: up(n_cg,nlev)           ! Parcel U (m/s)

REAL :: vp(n_cg,nlev)           ! Parcel V  (m/s)

REAL :: trap(n_cg,nlev,ntra)    ! Tracer content of parcel (kg/kg)

REAL :: expi(n_cg)              ! Initial parcel exner pressure

REAL :: xpk(n_cg,nlev)          ! Parcel cloud water (kg/kg)

REAL :: flx(n_cg,nlev)          ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(n_cg,nlev)     ! Minmum parcel buoyancy excess

REAL :: thpixs_v(n_cg,nlev)     ! Theta parcel excess (K)

REAL :: qpixs_v(n_cg,nlev)      ! Q parcel excess(kg/kg)


! PC2

REAL ::           &
  qclp(n_cg,nlev) & ! Parcel liquid condensated mixing ratio in layer k (kg/kg)
 ,qcfp(n_cg,nlev)   ! Parcel frozen condensated mixing ratio in layer k (kg/kg)

! Parameters

REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed


! CMT variables

INTEGER ::           &
  nlcl_uv(n_cg)      & ! Level index for LCL
 ,ntop_uv(n_cg)      & ! Level index for top of layer
 ,n_0degc(n_cg)      & ! Level index for zero degrees
 ,cu_term(n_cg)      & ! Indicies for CMT subs
 ,cu_tend(n_cg)        !

REAL ::                &
  exk_temp             & ! Temporary exner
 ,eflux_u_ud(n_cg)     & ! Vertical eddy flux of momentum due to UD at
                         ! top of layer (Pa m/s2)
 ,eflux_v_ud(n_cg)     & ! Vertical eddy flux of momentum due to UD at
                         ! bottom of layer (Pa m/s2)
 ,flxkp12(n_cg,nlev)   & ! Mass flux on half level (Pa/s)
 ,mb(n_cg)               ! Cloud base mass flux (Pa/s)

REAL :: p_uv(nlev,n_cg)         ! Pressure of model level (Pa)

REAL :: phalf_uv(nlev,n_cg)     ! Pressure of half level (Pa)

REAL :: plcl_uv(n_cg)           ! Pressure at LCL (Pa)

REAL :: ptop_uv(n_cg)           ! Pressure at top of cloud layer (Pa)

REAL :: p_0degc_uv(n_cg)        ! Pressure of zero degree level (Pa)

REAL :: rho_uv(nlev,n_cg)       ! Density on uv level (kg/m3)

REAL :: uw(nlev,n_cg)           ! U- comp stress profile (N/m2)
                                ! (units change through calls)

REAL :: ue_p(nlev,n_cg)         ! Environment U profile (m/s)

REAL :: vw(nlev,n_cg)           ! V-comp stress profile (N/m2)

REAL :: ve_p(nlev,n_cg)         ! Environment V profile (m/s)

REAL :: zcld(n_cg)              ! Depth of cloud layer (m)

LOGICAL :: l_mom_gk     ! Set to true if using Gregory-Kershaw CMT scheme


! CFL scaling variables

INTEGER ::           &
  det_lev(n_cg)      & ! Level at which split final detrainment last occurred
 ,nterm              & ! No. of points where conv. has terminated
 ,index_nterm(n_cg)    ! Index for points where conv. has terminated

REAL ::               &
  tempnum             & ! Temporary variable for storage
 ,scale_f(n_cg)         ! store scaling factor

! Downdraught scheme variables

INTEGER ::            &
  nnodd               & ! No. of downdraughts not possible
 ,index_nodd(n_cg)    & ! Index of downdraughts not possible
 ,npossdd             & ! No. downdraughts possible
 ,index_possdd(n_cg)  & ! Index of downdraughts
 ,kmax_term             ! maximum termination level + 1

REAL :: deltap_cld      ! pressure thickness of convective cloud (Pa)

! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag

INTEGER :: ntpar_max          ! max ntpar value

! parameters etc for qmin checks

REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

REAL :: qminincolumn(n_cg)     ! Minimum value for q in column
                               ! (kg/kg)
REAL :: temp1(n_cg)            ! work array

! Local compressed arrays

LOGICAL ::        &
  bconv_c2(n_cg)  & ! 
 ,bgmkp1_c(n_cg)  & ! Mask for points where parcel in layer k+1 is saturated
 ,bgmkp1_c2(n_cg) & ! Mask for points where parcel in layer k+1 is saturated
 ,bwk_c(n_cg)     & ! bwater mask in layer k
 ,bwk_c2(n_cg)    & ! bwater mask in layer k
 ,bwkp1_c(n_cg)   & ! bwater mask in layer k+1
 ,bwkp1_c2(n_cg)    ! bwater mask in layer k+1

REAL :: deltak_c2(n_cg)         ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

REAL ::            &
  dqek_c2(n_cg)    & ! Increment to q due to convection in layer k (kg/kg)
 ,dqekp1_c2(n_cg)  & ! Increment to q due to convection in layer k+1 (kg/kg)
 ,dthek_c2(n_cg)   & ! Increment to potential temp. due to conv in layer k
 ,dthekp1_c2(n_cg)   ! Increment to potential temp. due to conv in layer k+1

REAL ::                &
  dtraek_c2(n_cg,ntra) & ! Incr to model tracer due to conv. at lev k (kg/kg/s)
 ,dtraekp1_c2(n_cg,ntra) ! Incr to model tracer due to conv. at lev k+1(kg/kg/s)

REAL ::             &
  duek_c2(n_cg)     & ! Increment to model U in layer k due to CMT (m/s2)
 ,duekp1_c2(n_cg)   & ! Increment to model U in layer k+1 due to CMT (m/s2)
 ,dvek_c2(n_cg)     & ! Increment to model V in layer k due to CMT (m/s2)
 ,dvekp1_c2(n_cg)     ! Increment to model V in layer k+1 due to CMT (m/s2)

REAL ::           &
  flxk_c(n_cg)    & ! Parcel mass flux in layer k (Pa/s)
 ,flxk_c2(n_cg)   & ! Parcel mass flux in layer k (Pa/s)
 ,flxkp12_c2(n_cg)  ! Half level mass flux (Pa/s)

REAL ::           &
  prekp1_c2(n_cg) ! Precip from parcel as it rises from layer k to k+1 (kg/m2/s)

REAL ::           &
  qpk_c(n_cg)     & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpk_c2(n_cg)    & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpk(n_cg)       & ! Parcel mixing ratio in layer k (kg/kg)
 ,qpkp1_c(n_cg)   & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qpkp1_c2(n_cg)  & ! Parcel mixing ratio in layer k+1 (kg/kg)
 ,qek_c(n_cg)     & ! Env. mixing ratio in layer k (kg/kg)
 ,qek_c2(n_cg)    & ! Env. mixing ratio in layer k (kg/kg)
 ,qek(n_cg)       & ! Env. mixing ratio in layer k (kg/kg)
 ,qekp1_c(n_cg)   & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1_c2(n_cg)  & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qekp1(n_cg)     & ! Env. mixing ratio in layer k+1 (kgkg-1)
 ,qsek_c2(n_cg)   & ! Saturation mixing ratio of cld. env. in layer k (kg/kg)
 ,qsek(n_cg)      & ! Saturation mixing ratio of cld. env. in layer k (kg/kg)
 ,qsekp1_c(n_cg)  & ! Saturation mixing ratio of cld. env. in layer k+1 (kg/kg)
 ,qsekp1_c2(n_cg) & ! Saturation mixing ratio of cld. env. in layer k+1 (kg/kg)
 ,qsekp1(n_cg)      ! Saturation mixing ratio of cld. env. in layer k+1 (kg/kg)

REAL ::           &
  thek_c(n_cg)    & ! Env. potential temp in layer k (K)
 ,thek_c2(n_cg)   & ! Env. potential temp in layer k (K)
 ,thek(n_cg)      & ! Env. potential temp in layer k (K)
 ,thekp1_c(n_cg)  & ! Env. potential temp i in layer k+1 (K)
 ,thekp1_c2(n_cg) & ! Env. potential temp i in layer k+1 (K)
 ,thekp1(n_cg)    & ! Env. potential temp i in layer k+1 (K)
 ,thpk_c(n_cg)    & ! Parcel potential temp in layer k (K)
 ,thpk_c2(n_cg)   & ! Parcel potential temp in layer k (K)
 ,thpk(n_cg)      & ! Parcel potential temp in layer k (K)
 ,thpkp1_c(n_cg)  & ! Parcel potential temp in layer k+1 (K)
 ,thpkp1_c2(n_cg) & ! Parcel potential temp in layer k+1 (K)
 ,thpkp1(n_cg)      ! Parcel potential temp in layer k+1 (K)

REAL ::                 &
  traek_c(n_cg,ntra)    & ! Tracer content cld. env. in layer k (kg/kg)
 ,traek_c2(n_cg,ntra)   & ! Tracer content cld. env. in layer k (kg/kg)
 ,traekp1_c(n_cg,ntra)  & ! Tracer content of cloud env. in layer k+1 (kg/kg)
 ,traekp1_c2(n_cg,ntra) & ! Tracer content of cloud env. in layer k+1 (kg/kg)
 ,trapk_c(n_cg,ntra)    & ! Tracer cont.of parcel in layer k (kg/kg)
 ,trapk_c2(n_cg,ntra)   & ! Tracer cont.of parcel in layer k (kg/kg)
 ,trapkp1_c(n_cg,ntra)  & ! Tracer cont. of parcel in layer k+1 (kg/kg)
 ,trapkp1_c2(n_cg,ntra)   ! Tracer cont. of parcel in layer k+1 (kg/kg)

REAL :: rbuoy_c(n_cg), rbuoy_c2(n_cg) ! Buoyancy of parcel at k+1 (Kelvin)

REAL :: uek_c(n_cg), uek_c2(n_cg) ! Model U field on layer k (m/s)

REAL :: uekp1_c(n_cg), uekp1_c2(n_cg)! Model U field on layer k+1 (m/s)


REAL :: vek_c(n_cg), vek_c2(n_cg)     ! Model V field on layer k (m/s)

REAL :: vekp1_c(n_cg), vekp1_c2(n_cg) ! Model V field on layer k+1 (m/s)

REAL :: upk_c(n_cg), upk_c2(n_cg) ! Parcel U in layer k
                                  ! after entrainment (m/s)


REAL :: upkp1_c(n_cg), upkp1_c2(n_cg) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

REAL :: vpk_c(n_cg), vpk_c2(n_cg) ! Parcel V in layer k
                                  ! after entrainment (m/s)

REAL :: vpkp1_c(n_cg), vpkp1_c2(n_cg) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

REAL :: xsqkp1_c(n_cg), xsqkp1_c2(n_cg) ! Excess water vapour
                                        ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

REAL ::            &
  qclek_c(n_cg)    & ! Env liquid condensate mixing ratio in layer k (kg/kg)
 ,qclek_c2(n_cg)   & ! Env liquid condensate mixing ratio in layer k (kg/kg)
 ,qclekp1_c(n_cg)  & ! Env liquid condensate mixing ratio in layer k+1 (kg/kg)
 ,qclekp1_c2(n_cg) & ! Env liquid condensate mixing ratio in layer k+1 (kg/kg)
 ,qcfek_c(n_cg)    & ! Env frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfek_c2(n_cg)   & ! Env frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfekp1_c(n_cg)  & ! Env frozen condensate mixing ratio in layer k+1 (kg/kg)
 ,qcfekp1_c2(n_cg) & ! Env frozen condensate mixing ratio in layer k+1 (kg/kg)
 ,qclpk_c(n_cg)    & ! Parcel liquid condensate mixing ratio in layer k (kg/kg)
 ,qclpk_c2(n_cg)   & ! Parcel liquid condensate mixing ratio in layer k (kg/kg)
 ,qclpkp1_c(n_cg)  & ! Parcel liquid condensate mixing ratio in layer k+1(kg/kg)
 ,qclpkp1_c2(n_cg) & ! Parcel liquid condensate mixing ratio in layer k+1(kg/kg)
 ,qcfpk_c(n_cg)    & ! Parcel frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfpk_c2(n_cg)   & ! Parcel frozen condensate mixing ratio in layer k (kg/kg)
 ,qcfpkp1_c(n_cg)  & ! Parcel frozen condensate mixing ratio in layer k+1(kg/kg)
 ,qcfpkp1_c2(n_cg)   ! Parcel frozen condensate mixing ratio in layer k+1(kg/kg)

REAL ::            &
  cflek_c2(n_cg)   & ! Environment liquid water cloud volume ( )
 ,cflekp1_c2(n_cg) & ! Environment liquid water cloud volume ( )
 ,cffek_c2(n_cg)   & ! Environment frozen water cloud volume ( )
 ,cffekp1_c2(n_cg) & ! Environment frozen water cloud volume ( )
 ,bcfek_c2(n_cg)   & ! Environment bulk total cloud volume ( )
 ,bcfekp1_c2(n_cg)   ! Environment bulk total cloud volume ( )

REAL ::            &
  dqclek_c2(n_cg)  & ! Environment increments to liquid condensate mixing
 ,dqclekp1_c2(n_cg)& ! ratio to convection (kg/kg/s)
 ,dqcfek_c2(n_cg)  & ! Environment increments to frozen condensate mixing
 ,dqcfekp1_c2(n_cg)  !  ratio to convection (kg/kg/s)

REAL ::             &
  dcflek_c2(n_cg)   & ! Environment increments to liquid water cloud 
 ,dcflekp1_c2(n_cg) & ! volume due to convection (/s)
 ,dcffek_c2(n_cg)   & ! Environment increments to frozen water cloud
 ,dcffekp1_c2(n_cg) & ! volume due to convection (/s)
 ,dbcfek_c2(n_cg)   & ! Environment increments to bulk total cloud
 ,dbcfekp1_c2(n_cg)   ! volume due to convection (/s)

REAL :: amdetk_c2(n_cg)
LOGICAL :: bgmk_c2(n_cg)
LOGICAL :: bland_c2(n_cg)
LOGICAL :: blowst_c2(n_cg)
LOGICAL :: bterm_c2(n_cg)
REAL :: cape_c2(n_cg)
REAL :: cca_2d_c2(n_cg)
REAL :: cclwp_c2(n_cg)
REAL :: ccw_c2(n_cg)
LOGICAL :: cumulus_c(n_cg), cumulus_c2(n_cg)
REAL :: dcpbydt_c2(n_cg)
REAL :: delexkp1_c2(n_cg)
REAL :: delpk_c2(n_cg)
REAL :: delpkp1_c2(n_cg)
REAL :: delp_uv_k_c2(n_cg)
REAL :: delp_uv_kp1_c2(n_cg)
REAL :: depth_c2(n_cg)
REAL :: dptot_c2(n_cg)
REAL :: dqsthkp1_c2(n_cg)
REAL :: dqsthk_c2(n_cg)
REAL :: eflux_u_ud_c2(n_cg)
REAL :: eflux_v_ud_c2(n_cg)
REAL :: ekp14_c(n_cg),ekp14_c2(n_cg)
REAL :: ekp34_c(n_cg),ekp34_c2(n_cg)
REAL :: exk_c2(n_cg)
REAL :: exkp1_c(n_cg),exkp1_c2(n_cg)
REAL :: expi_c2(n_cg)
INTEGER :: icct_c2(n_cg)
INTEGER :: iccb_c2(n_cg)
INTEGER :: lctop_c2(n_cg)
INTEGER :: lcbase_c2(n_cg)
REAL :: lcca_c2(n_cg)
LOGICAL :: l_shallow_c2(n_cg)
LOGICAL :: l_mid_c2(n_cg)
REAL :: max_cfl_c2(n_cg)
REAL :: pk_c(n_cg),pk_c2(n_cg)
REAL :: pkp1_c(n_cg),pkp1_c2(n_cg)
REAL :: pstar_c2(n_cg)
REAL :: q1_sd_c2(n_cg)
REAL :: qpi_c2(n_cg)
REAL :: qpixs_v_c2(n_cg)
REAL :: relh_c2(n_cg)
REAL :: rbuoy_p_here_c2(n_cg)
REAL :: the_here_c2(n_cg)
REAL :: thp_here_c2(n_cg)
REAL :: qe_here_c2(n_cg)
REAL :: qp_here_c2(n_cg)
REAL :: rbuoy_p_old_c2(n_cg)
REAL :: tcw_c2(n_cg)
REAL :: thpi_c2(n_cg)
REAL :: thpixs_v_c2(n_cg)
REAL :: t1_sd_c2(n_cg)
REAL :: xpk_c(n_cg),xpk_c2(n_cg)
REAL :: xsbmin_v_c2(n_cg)
REAL :: qsat_lcl(n_cg)    ! not used
LOGICAL :: b_nodd(n_cg)   ! points with no downdraught
LOGICAL :: b_dd(n_cg)     ! points with downdraught on termination

!===============================================================
! CCRad Variables local variables As SHALLOW
!===============================================================

REAL   :: overlap_fac(n_cg)  ! Factor designed to improve
                             ! shallow Cu cover by allowing
                             ! for non-vertical clouds.

REAL   :: zpr         ! method (BL Fluxes) only if
                      !   l_ccrad = T .AND. cca2d_sh_opt  = 1

INTEGER, PARAMETER   :: cca2d_total_condensed_water = 0
INTEGER, PARAMETER   :: cca2d_grant_lock            = 1 ! Shallow cnv
INTEGER, PARAMETER   :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv


!===============================================================
! End CCRad Variables local variables
!===============================================================

! Loop counters

INTEGER :: i,i2,j,k,ktra,kt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



IF (lhook) CALL dr_hook('CONGEST_CONV',zhook_in,zhook_handle)

!initialise SCM diagnostics
DO k = 1,nlev
  DO i = 1,n_cg
    rbuoy_p_out(i,k)=0.0
    the_out(i,k)=th(i,k)
    thp_out(i,k)=th(i,k)
    qe_out(i,k)=q(i,k)
    qp_out(i,k)=q(i,k)
  END DO
END DO

! Initialise logicals

DO i = 1,n_cg
  blowst(i)    = .TRUE.
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  cumulus(i)   = .TRUE.
  l_shallow(i) = .TRUE.
  l_mid(i)     = .FALSE.
  b_nodd(i)    = .FALSE.
  b_dd(i)      = .FALSE.
END DO

l_mom_gk = .FALSE.    ! not using Gregory-Kershaw scheme

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------
! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)

DO k = 1,nlev-1
  DO i = 1,n_cg
    xsbmin_v(i,k) = MIN( ((p_layer_centres(i,k) -                 &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = MIN( ((p_layer_centres(i,k) -                 &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_shallow

    qpixs_v(i,k)  = qpixs_shallow
  END DO
END DO  ! nlev

! Calculate convective velocity scale and cloud base mass flux

DO i = 1,n_cg
  wsc(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i))   &
             * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
  mb(i)  = c_mass * wstar(i)
  zcld(i) = ztop_uv(i) - zlcl_uv(i)
  wsc_o_mb(i) = wsc(i)/mb(i)
END DO

! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.


DO i = 1,n_cg
  k =ntml(i)
  p_lcl(i)  = p_layer_boundaries(i,k)
  th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
  t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
  q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
END DO

! Calculate saturation mixing ratio at LCL


! DEPENDS ON: qsat_mix
CALL qsat_mix(qse_lcl,t_lcl,p_lcl,n_cg,.FALSE.)


! Note: In terms of array indices p and phalf follow the convention
!       used in the boundary layer scheme. phalf(k,*) refers to the
!       lower boundary of uv layer k. This follows the convention for
!       um UM4.5 and before
!
!       Also note that p_layer_boundaries(0) and p_layer_centres(0)
!       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.
!
!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)

! Calculate theta and q pertubations (pertubation is based on
! environment buoyancy gradient)
! Reset th and q xs's at ntml

DO i = 1,n_cg
  IF (t_lcl(i) >  tm) THEN
    dq_sat_env = repsilon * lc * qse_lcl(i) / (r * t_lcl(i) * t_lcl(i))
  ELSE
    dq_sat_env = repsilon * (lc+lf) * qse_lcl(i) / (r * t_lcl(i) * t_lcl(i))
  END IF

  b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0              &
               + c_virtual * qse_lcl(i)

  thv_pert(i) = -0.17 * wthvs(i) / mb(i) +                        &
              ( th(i,ntml(i)+1) *(1.0 + c_virtual*q(i,ntml(i)+1)) &
              - th(i,ntml(i))  * (1.0 + c_virtual*q(i,ntml(i))) )


  c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i)) - thv_pert(i)

  thpert(i)   = -c_calc / b_calc  !ignore term in THPERT**2

  thpixs_v(i,ntml(i)) = thpert(i)

  qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)**kappa)             &
                             * thpert(i) * dq_sat_env - q_lcl(i)

  qpixs_v(i,ntml(i))  = qpert(i)

END DO !n_cg


! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)


! DEPENDS ON: flag_wet
CALL flag_wet(n_cg,n_cg,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

DO k = 1,nlev
  DO i = 1,n_cg
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
    xpk(i,k)    = 0.0
    dthbydt(i,k) = 0.0
    dqbydt(i,k)  = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
  END DO
END DO

IF (l_mom) THEN
  DO k = 1,nlev+1
    DO i = 1,n_cg
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,n_cg
        dtrabydt(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF  ! L_tracer


!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------

IF (flg_up_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      up_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      up_flux_half(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      dwn_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      entrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      detrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      entrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_cg
      detrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF

IF (l_mom) THEN
  IF (flg_uw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_cg
        uw_shall(i,k) = 0.0
      END DO
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_cg
        vw_shall(i,k) = 0.0
      END DO
    END DO
  END IF
END IF  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
DO i = 1,n_cg
  cca_2d(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
END DO

DO k=1, n_cca_lev
  DO i=1, n_cg
    cca(i,k) = 0.0
  END DO
END DO

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics - done in glue routine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling calculations
!-----------------------------------------------------------------------

DO i = 1,n_cg
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)      = 0.0
  cape_out(i)  = 0.0
  dcpbydt(i)   = 0.0
  max_cfl(i)   = 0.0
  det_lev(i)   = 0

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


! set SCM adaptive diagnostics for level k = 1

DO i = 1,n_cg
  rbuoy_p_out(i,1) = 0.0
  the_out(i,1) = th(i,1)
  thp_out(i,1) = th(i,1)
  qe_out(i,1) = q(i,1)
  qp_out(i,1) = q(i,1)
END DO
! Also, initialise rbuoy_p_old, i.e. for previous level,
! to 0.0 for level1
DO i = 1,n_cg
  rbuoy_p_old(i) = 0.0
END DO
!initialise ekm14
DO i =1, n_cg
  ekm14(i) =0.0
END DO
!Initialise adaptive entrainment variables
!intitilaise to level 2 'cos that's where parcel lift starts from
DO i = 1, n_cg
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
  DO i = 1, n_cg
    qp(i,k) = q(i,k)
    thp(i,k) = th(i,k)
  END DO
END DO


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

IF (cg_on == 1) THEN   ! Adaptive forced detrainment on
                       ! No limit on convection top
  ntpar_max = nlev-3   ! What is a sensible value to have here?

ELSE                   ! Top limited
  ntpar_max=0
  DO i = 1,n_cg
    IF (ntpar(i) >  ntpar_max) THEN
      ntpar_max=ntpar(i)
    END IF
  END DO
END IF


! Add a level (needed to ensure tracers give same answers on different
! CPU configurations). Not sure why.
ntpar_max=ntpar_max+1

! This test should not really be required as don't expect any
! shallow convection to reach model top unless a very funny set of
! model levels (very shallow atmosphere).

IF (ntpar_max == nlev) THEN
  ntpar_max=nlev-1
END IF


!      Do k = 2,nlev-1       ! original level loop

DO k = 2,ntpar_max

! Set relative humidity in layer k (rhum)

  DO i = 1,n_cg
    rhum(i) = q(i,k) / qse(i,k)
  END DO

!  Initialise adaptive diagnostics for this level

  DO i = 1,n_cg
    rbuoy_p_here(i) =0.0
    the_here(i) = th(i,k)
    qe_here(i) = q(i,k)
    rbuoy_p_here_c2(i) =0.0
    the_here_c2(i) = 0.0
    thp_here_c2(i) = 0.0
    qe_here_c2(i) = 0.0
    qp_here_c2(i) = 0.0
  END DO

!Initialise adaptive entrainment variables
  DO i = 1, n_cg
    thek(i)=th(i,k)
    qek(i)=q(i,k)
    qsek(i)=qse(i,k)
    thekp1(i)=th(i,k+1)
    qekp1(i)=q(i,k+1)
    qsekp1(i)=qse(i,k+1)
    IF(k  ==  2) THEN      !set to environmental values
      thpk(i)=th(i,2)
      qpk(i)=q(i,2)
  !only if bconv is true do we have non-zero value for thp(i,k)
  !in this iteration of the main loop - so set to environment
  !value
    ELSE IF (thp(i,k)  <   1) THEN
      thpk(i)=th(i,2)
      qpk(i)=q(i,2)
    ELSE
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

  DO i = 1,n_cg
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
      qpk(i) = qp(i,k)         !set qpk for first convecting level
    END IF
  END DO  ! n_cg


!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
  CALL layer_cn(k,n_cg,nlev                                       &
,               mdet_cg_on, cg_ent_on                             &
,               ntml,ntpar                                        &
,               .FALSE.,.TRUE.,.FALSE.                            &
,               bconv,bwk,bwkp1                                   &
,               exner_layer_boundaries                            &
,               exner_layer_centres                               &
,               p_layer_boundaries,p_layer_centres                &
,               recip_pstar,entrain_coef,rhum                     &
,               zk, zkp12, zkp1                                   &
,               thek, qek,qsek, thekp1,qekp1,qsekp1               &
,               thpk,qpk ,qsat_lcl, ekm14                         &
,               pkp1,delpkp1,exkp1                                &
,               pk,delpk,delpkp12,exk,delexkp1                    &
,               delp_uv_k, delp_uv_kp1                            &
,               ekp14,ekp34,amdetk)


! Set ekm14 for next pass through loop
  DO i = 1, n_cg
    ekm14(i) = ekp14(i)
  END DO


! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)


  IF (k == 2) THEN
! DEPENDS ON: dqs_dth
    CALL dqs_dth(dqsthk,k,th(1,k),qse(1,k),exk,n_cg)
  ELSE
    DO i = 1,n_cg
      dqsthk(i) = dqsthkp1(i)
    END DO
  END IF

! DEPENDS ON: dqs_dth
  CALL dqs_dth(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_cg)


! Set other grid dependent constants


! Maximum initial convective mass flux

  DO i = 1,n_cg

    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

  END DO

! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k


  DO i = 1,n_cg
! not convecting and not convected in column before
! PC2 qclp and qcfp zero at this point but will add an initial
! value at cloud base
    IF ( .NOT. bconv(i).AND.det_lev(i) == 0) THEN
      expi(i)  = exk(i)
      xpk(i,k) = 0.0
      qclp(i,k) = 0.0
      qcfp(i,k) = 0.0
      flx(i,k) = 0.0
      bgmk(i)  = .FALSE.
      depth(i) = 0.0
      thpi(i)  = th(i,k) + thpixs_v(i,k)
      thp(i,k) = thpi(i)
      qpi(i)   = q(i,k) + qpixs_v(i,k)
      qp(i,k)  = qpi(i)
      IF (l_mom_gk) THEN
        up(i,k) = u(i,k)
        vp(i,k) = v(i,k)
      END IF
    END IF
  END DO  ! n_cg
  IF (l_tracer) THEN
    DO ktra=1,ntra
      DO i = 1,n_cg
        IF ( .NOT. bconv(i)) THEN
          trap(i,k,ktra)  = tracer(i,k,ktra)
        END IF  !not bconv
      END DO
    END DO
  END IF


! Scale entrainment coefficients with cloud base mass flux
! and convective velocity scale

  DO i = 1,n_cg

! Leaving commented out code as may use again
!          If (k  >=  ntml(i)) then

    ! If original entrainment coefficient then scale as before else
    ! leave as calculated by layer_cn.

!            If (entrain_coef(i) < 0.0 .or. icvdiag == 7 ) then
!              ekp14(i)  = ekp14(i) * wsc_o_mb(i)
!              ekp34(i)  = ekp34(i) * wsc_o_mb(i)
!              amdetk(i) = amdetk(i) * wsc_o_mb(i)
!            End If
!          End If

! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k = ntml

    bcposs(i) = bconv(i) .OR. k  ==  ntml(i)

  END DO  ! n_cg

! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)

  ncposs = 0
  DO i = 1,n_cg
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
    IF (l_q_interact) THEN
      DO i = 1,ncposs
!           PC2 variables
        qclek_c(i)    = qcl(index1(i),k)
        qclekp1_c(i)  = qcl(index1(i),k+1)
        qcfek_c(i)    = qcf(index1(i),k)
        qcfekp1_c(i)  = qcf(index1(i),k+1)
        qclpk_c(i)    = qclp(index1(i),k)
        qcfpk_c(i)    = qcfp(index1(i),k)
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

! DEPENDS ON: lift_par_5a
  CALL lift_par_5a(ncposs,n_cg,thpkp1_c,qpkp1_c,xsqkp1_c,         &
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
      bconv(index1(i)) = .TRUE.
!          End If


! Set parcel mass flux
! UM Documentation paper 27, section 1.5

!          If (bconv(index1(i)).and.k == ntml(index1(i))) then
! as k=ntml
      flxk_c(i) = mb(index1(i)) * g * p_layer_centres(index1(i),k)&
                  / (r * thpk_c(i) *                              &
                     (p_layer_centres(index1(i),k)/ pref)**kappa)


! Write compressed mass flux back to full array

      flx(index1(i),k) = flxk_c(i)

! At ntml set mixing detrainment rate equal to zero
! Store diagnostics linked to initial convective mass flux for
! calculation of final closure.

      amdetk(index1(i))      = 0.0
      flx_init(index1(i))    = flxk_c(i)
      flxmax_init(index1(i)) = flxmax(index1(i))

      blowst(index1(i)) = .TRUE.

      IF (l_q_interact) THEN

!           Initialize QCLP(*,K) and QCFP(*,K) at start level and
!           perform the Parcel Lift again for such points: needed
!           for PC2.   Duplicates code from
!           Convection Parcel Lifting Scheme, LIFPAR.
!           Assume "if k eq ntml" equiv "if bwork(i,3) in convec4a
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

    ELSE
      blowst(index1(i)) = .FALSE.   ! not initial layer

    END IF



! Reset threashold for forced detrainment to the initial (negative)
! buoyancy


    xsbmin_v(index1(i),k) = thv_pert(index1(i))

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
      blowst_c2(i) = blowst(index1(index2(i)))
      l_shallow_c2(i) = .TRUE.
      l_mid_c2(i)     = .FALSE.
      cumulus_c2(i)   = .TRUE.
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
      uek_c2(i)    = uek_c(index2(i))
      uekp1_c2(i)  = uekp1_c(index2(i))
      vek_c2(i)    = vek_c(index2(i))
      vekp1_c2(i)  = vekp1_c(index2(i))
      upkp1_c2(i)    = upkp1_c(index2(i))
      vpkp1_c2(i)    = vpkp1_c(index2(i))
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
      upk_c2(i)    = upk_c(index2(i))
      vpk_c2(i)    = vpk_c(index2(i))
      xpk_c2(i)    = xpk_c(index2(i))
      flxk_c2(i)   = flx(index1(index2(i)),k)
      bgmk_c2(i)   = bgmk(index1(index2(i)))
      bterm_c2(i)  = .FALSE.
      dthek_c2(i)  = dthbydt(index1(index2(i)),k)
      dqek_c2(i)   = dqbydt(index1(index2(i)),k)
      dthekp1_c2(i)  = dthbydt(index1(index2(i)),k+1)
      dqekp1_c2(i)   = dqbydt(index1(index2(i)),k+1)
      tcw_c2(i)    = tcw(index1(index2(i)))
      depth_c2(i)  = depth(index1(index2(i)))
      cclwp_c2(i)  = cclwp(index1(index2(i)))
      cape_c2(i)    = cape(index1(index2(i)))
      dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
      relh_c2(i)    = 0.0 ! dummy variable
      dptot_c2(i)   = 0.0 ! dummy variable
      rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
      the_here_c2(i)=the_here(index1(index2(i)))
      thp_here_c2(i)=thp_here(index1(index2(i)))
      qe_here_c2(i)=qe_here(index1(index2(i)))
      qp_here_c2(i)=qp_here(index1(index2(i)))
      rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
      eflux_u_ud_c2(i) = eflux_u_ud(index1(index2(i)))
      eflux_v_ud_c2(i) = eflux_v_ud(index1(index2(i)))
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
        duek_c2(i)   = dubydt(index1(index2(i)),k)
        dvek_c2(i)   = dvbydt(index1(index2(i)),k)
        duekp1_c2(i)   = dubydt(index1(index2(i)),k+1)
        dvekp1_c2(i)   = dvbydt(index1(index2(i)),k+1)
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

! Force congestus convection to stop at ntpar unless using adaptive forced
! detrainment.

    IF (cg_on == 0) THEN
! Part of BTERM comes from values in array before convec2.
! Original code set bterm_c2 according to this test before call to
! convec2.
      DO i = 1,nconv
        IF (k  ==  ntpar(index1(index2(i))))THEN
          bterm_c2(i) = .TRUE.
        END IF
      END DO

    END IF
!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.

!      Subroutine CONVEC2

!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2_4a5a
    CALL convec2_4a5a(nconv,n_cg,nlev,ntra,k,cg_on, cg_sdet_on, cg_new_termc,&
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


    IF (flg_entr_up) THEN
      DO i = 1,nconv
        entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *      &
                    (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i) &
                   * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
        IF (bterm_c2(i)) THEN
          entrain_up(index1(index2(i)),k+1) = 0.0
        END IF
      END DO
    END IF


! Calculate fractional detrainment rate for level k
!(and k+1 if bterm=.T.)


    IF (flg_detr_up) THEN
      DO i = 1,nconv
        detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)              &
                        + deltak_c2(i) * (1.0 - amdetk_c2(i)))        &
                      * flx(index1(index2(i)),k)
        IF (bterm_c2(i)) THEN
          detrain_up(index1(index2(i)),k+1) =                        &
                 -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
        END IF
      END DO
    END IF

  END IF   ! nconv > 0

! Write CONVEC2 compressed output arrays back to full fields

  DO i = 1,n_cg
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
      DO i = 1,n_cg
        trap(i,k+1,ktra) = 0.0
      END DO
    END DO
  END IF

  IF (l_mom_gk) THEN
    DO i = 1,n_cg
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
    END DO
    DO i = 1,nconv
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
      qe_here(index1(index2(i))) = qe_here_c2(i)
      qp_here(index1(index2(i))) = qp_here_c2(i)
    END DO
    DO i = 1,nconv
      max_cfl(index1(index2(i))) =                                  &
                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
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
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        dubydt(index1(index2(i)),k) = duek_c2(i)
        dvbydt(index1(index2(i)),k) = dvek_c2(i)
        dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
        dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
        eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
        eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
      END DO
    END IF
    IF (l_mom) THEN      ! require whatever scheme used
      DO i = 1,nconv
        flxkp12(index1(index2(i)),k)  = flxkp12_c2(i)
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
  END IF      ! nconv > 0

!   Write adaptive diagnostics for this level to full array for output

  DO i = 1,n_cg
    rbuoy_p_out(i,k) = rbuoy_p_here(i)
    the_out(i,k) = the_here(i)
    thp_out(i,k) = thp_here(i)
    qe_out(i,k) = qe_here(i)
    qp_out(i,k) = qp_here(i)
  END DO
!   Write rbuoy for this level to rbuoy_p_old for previous level
!   for use in parcel

  DO i = 1,n_cg
    rbuoy_p_old(i) = rbuoy_p_here(i)
  END DO


!-----------------------------------------------------------------------
! 3.4  CFL scaling
!-----------------------------------------------------------------------

! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_cg) with index_nterm


  nterm = 0
  DO i = 1,n_cg
    IF (bterm(i)) THEN
      nterm = nterm + 1
      index_nterm(nterm) = i
    END IF
  END DO

  IF (nterm >  0) THEN


! Work out scaled mass flux needed to keep cfl ratio below limit.
! Note L_CAPE not applied to shallow convection

    DO j = 1,nterm
      i = index_nterm(j)

      max_cfl(i) = max_cfl(i) * timestep
      IF (max_cfl(i)  >   cfl_limit) THEN
        flx_init_new(i) = flx_init(i) * cfl_limit / max_cfl(i)
      ELSE
        flx_init_new(i) = flx_init(i)
      END IF

      IF (flx_init_new(i)  >   flxmax_init(i)) THEN
        flx_init_new(i) = flxmax_init(i)
      END IF
      max_cfl(i) = 0.0
    END DO      ! j (nterm)

! Scale cloud fraction

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
      !
        cca_2d(i) = MAX(2.0e-5, cca_2d(i))
        cca_2d(i) = MIN(1.0e+0, cca_2d(i))

      END DO      ! j (nterm)

    ELSE ! original code

      DO j = 1,nterm
        i = index_nterm(j)
        IF (flx_init_new(i)  >   0.0) THEN
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
            IF (kt <  k+1) THEN
              flxkp12(i,kt) = flxkp12(i,kt) * scale_f(i)
            END IF
          END IF
          IF (l_tracer) THEN
            DO ktra = 1,ntra
              dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            END DO
          END IF

          flx(i,kt)    = flx(i,kt) *  scale_f(i)
          precip(i,kt) = precip(i,kt) *  scale_f(i)

          IF (flg_entr_up) THEN
            entrain_up(i,kt) = entrain_up(i,kt) *  scale_f(i)
          END IF
          IF (flg_detr_up) THEN
            detrain_up(i,kt) = detrain_up(i,kt) *  scale_f(i)
          END IF

        END IF !kt >ntml and flx_init_new >0
      END DO  ! j loop
    END DO  ! kt loop

! Set final detrainment level (but not used).

    DO j = 1,nterm
      i = index_nterm(j)
      det_lev(i)= k+1
    END DO  ! nterm loop

!-----------------------------------------------------------------------
! 3.5  Downdraught calculation - on all points where convection is
!      terminating. Downdraughts are possible for some deeper shallow
!      convection.
!
!      Subroutine DD_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

    npossdd=0
    nnodd = 0

    DO i = 1,nterm
      i2=index_nterm(i)
      tempnum=0.0
      IF(iccb(i2) >  0) THEN
        deltap_cld=p_layer_centres(i2,iccb(i2))-p_layer_centres(i2,k)

        DO kt=iccb(i2),k+1
          tempnum=tempnum+precip(i2,kt)
        END DO
      ELSE
        deltap_cld = 0.0
      END IF

! Set logicals for use later

      IF (deltap_cld >  15000.0.AND.bgmk(i2)                        &
                              .AND.tempnum >  1e-12) THEN
        b_dd(i2) = .TRUE.
      ELSE
        b_nodd(i2) = .TRUE.
      END IF
    END DO  ! nterm loop


! If convection has terminated write cape to diagnostic output
! variable (cape_out).
! Set kterm array which holds the level index for termination
! of convection.

    DO j = 1,nterm
      i=index_nterm(j)
      cape_out(i) = cape(i)
      dcpbydt(i) = 0.0
      cape(i) = 0.0
      kterm(i) = k
      bconv(i) = .FALSE.
    END DO

  END IF  ! nterm > 0

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

END DO
! Mass flux diagnostics
IF (flg_up_flx) THEN
  DO k=1,nlev
    DO i = 1,n_cg
      up_flux(i,k) = flx(i,k)
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k=1,nlev-1
    DO i = 1,n_cg
      up_flux_half(i,k+1) = flxkp12(i,k)
    END DO
  END DO
END IF


!-----------------------------------------------------------------------
! 4.0 All shallow convection will terminate at some level. This level
!     has been stored in the main level loop.
!     The convection will either have a down draught or none will be
!     possible.
!-----------------------------------------------------------------------
! 4.1  Downdraft calculation - on all points where convection is
!      terminating.

!      Subroutine DD_ALL_CALL

!      UM Documentation Paper 27, part 2

!-----------------------------------------------------------------------

npossdd = 0
DO i = 1,n_cg
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
  CALL dd_all_call_4a5a (n_cg,npossdd,kmax_term,nlev,trlev,ntra   &
,                      kterm, iccb, icct, index_possdd, l_tracer  &
,                      bwater(1,2)                                &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      recip_pstar,timestep , cca_2d              &
,                      thp, qp, th, q, qse, trap,tracer, flx,precip &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow, rain_3d, snow_3d               &
,                      dwn_flux, entrain_dwn, detrain_dwn)

END IF

!-----------------------------------------------------------------------
! 4.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
nnodd = 0
DO i = 1,n_cg

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

! Only add 1 if kmax_term is less than model levels (which should be
! true).
  IF (kmax_term  <  nlev ) THEN
    kmax_term = kmax_term + 1
  END IF


! Surface precipitation calculation

! DEPENDS ON: evap_bcb_nodd_all
  CALL evap_bcb_nodd_all(n_cg,nnodd,kmax_term,kterm               &
                     , iccb,index_nodd,bwater(1,2)                &
                     , exner_layer_centres,exner_layer_boundaries &
                     , p_layer_centres,p_layer_boundaries,pstar   &
                     , timestep, cca_2d, th, q, qse, precip       &
                     , dthbydt,dqbydt                             &
                     , rain, snow, rain_3d, snow_3d)

END IF

! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).

IF (.NOT. l_ccrad) THEN
  DO i = 1,n_cg
    IF (iccb(i)  ==  icct(i)) THEN
      iccb(i)   = 0
      icct(i)   = 0
      cca_2d(i) = 0.0
      tcw(i)    = 0.0
      cclwp(i)  = 0.0
    END IF
    IF (lcbase(i)  ==  lctop(i)) THEN
      lcbase(i) = 0
      lctop(i)  = 0
      lcca(i)   = 0.0
    END IF
  END DO
END IF      ! l_ccrad

!-----------------------------------------------------------------------
! 5.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

IF (l_mom) THEN


! Initialize arrays required for Convective Momentum Transport(CMT)

  k=1
  DO i = 1,n_cg
    p_uv(k,i)     = p_layer_boundaries(i,k-1)
    phalf_uv(k,i) = p_layer_centres(i,k-1)
    ue_p(k,i)     = u(i,k)
    ve_p(k,i)     = v(i,k)
  END DO

  DO i = 1,n_cg
    nlcl_uv(i)    = ntml(i) + 1
    n_0degc(i)    = freeze_lev(i)
  END DO

  DO i = 1,n_cg
    DO k = 2,nlev
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      exk_temp      = (p_uv(k,i)/pref)**kappa
      rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp *           &
                      (th(i,k-1) + th(i,k)))
    END DO
    plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
    p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
    rho_uv(1,i)     = rho_uv(2,i)
  END DO

! altered to use kterm instead of ntpar
  DO i=1,n_cg
    IF (kterm(i) >= nlcl_uv(i)) THEN
      ntop_uv(i) = kterm(i) +1
    ELSE     ! case where congestus convection fails
! I think in this case the mass flux will be zero so no CMT
      ntop_uv(i) = ntpar(i) + 1
    END IF
    ptop_uv(i) = phalf_uv(ntop_uv(i),i)
  END DO

! Calculate CMT for required points

  nterm = 0

  DO i = 1, n_cg
    nterm = nterm + 1
    cu_term(nterm) = i
    cu_tend(nterm) = i
  END DO

! Note using shallow CMT assumptions but may be using top from kterm
! May not be a sensible choice.

  IF (nterm  >   0) THEN

! DEPENDS ON: shallow_grad_stress
    CALL shallow_grad_stress(n_cg,n_cg,nterm,nlev,cu_term,        &
                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,   &
                             plcl_uv,ptop_uv,p_uv,phalf_uv,       &
                             rho_uv,ue_p,ve_p,timestep,           &
                             ! IN
                             uw,vw)

! DEPENDS ON: shallow_base_stress
    CALL shallow_base_stress(n_cg,n_cg,n_cg,nlev,nterm,cu_term,   &
                             cu_tend,nlcl_uv,ntop_uv,mb,wsc,      &
                             zlcl_uv,zcld,uw0,vw0,plcl_uv,        &
                             ptop_uv,ue_p,ve_p,phalf_uv,p_uv,     &
                             rho_uv,timestep,flg_uw_shall,        &
                             flg_vw_shall,                        &
                             ! INOUT
                             uw,vw,                               &
                             ! OUT
                             uw_shall,vw_shall)


! DEPENDS ON: shallow_cmt_incr
    CALL shallow_cmt_incr(n_cg,n_cg,n_cg,nlev,nterm,cu_term,      &
                          cu_tend,nlcl_uv,ntop_uv,uw,vw,phalf_uv, &
                          rho_uv,zlcl_uv,                         &
                          !OUT
                          dubydt,dvbydt)

  END IF  ! nterm>0
END IF ! L_mom

!-----------------------------------------------------------------------
! 6.0  Energy correction calculation - removed as old code not correct
!     for new dynamics grid ( attempts to correct this give problems).
!     UM documentation paper 27 - section 12.
!-----------------------------------------------------------------------
DO i = 1,n_cg
  index1(i) = i
END DO

!-----------------------------------------------------------------------
! 7.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.

DO i = 1,n_cg
  qminincolumn(i) = q(i,nlev)
END DO
DO k = 1,nlev-1
  DO i = 1,n_cg
    IF (q(i,k)  <   qminincolumn(i)) THEN
      qminincolumn(i) = q(i,k)
    END IF
  END DO
END DO


! Ensure Q does not go below global allowed minimum (QMIN)

DO i = 1,n_cg
  qminincolumn(i)=MAX(qmin,qminincolumn(i))
END DO


! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.

DO k = nlev,2,-1
  DO i = 1,n_cg
    IF (dqbydt(i,k) /= 0.0) THEN
      temp1(i)=q(i,k) + dqbydt(i,k) * timestep
      IF (temp1(i)  <   qminincolumn(i)) THEN

        dqbydt(i,k-1) = dqbydt(i,k-1) -                         &
            ((qminincolumn(i) - q(i,k)) / timestep-dqbydt(i,k)) &
             * (r2rho_th(i,k)*dr_across_th(i,k))                &
            / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

        dqbydt(i,k) = (qminincolumn(i) - q(i,k)) / timestep
      END IF
    END IF
  END DO ! n_s loop
END DO  ! nlev

! check negative q

k=1
  DO i = 1,n_cg
    temp1(i)=q(i,k) + dqbydt(i,k) * timestep
    IF (temp1(i)  <   qminincolumn(i) .AND.                            &
        printstatus >= prstatus_normal) THEN
      WRITE(6,'(a20,i6,a9,g26.18,a7,g26.18)') ' negative q cong, i:',  &
            i,' q after ',temp1(i),' dq/dt ',dqbydt(i,k)
    END IF
  END DO ! n_cg loop

!-----------------------------------------------------------------------
! 8.0  Mixing of the convective increments in the boundary
!      layer.
!-----------------------------------------------------------------------

IF (bl_cnv_mix == 1) THEN

!      Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.

! DEPENDS ON: mix_ipert_4a5a
  CALL mix_ipert_4a5a(n_cg, nlev, nbl, ntml, p_layer_boundaries,  &
               exner_layer_centres, dthbydt, dqbydt, flx_init,    &
               thpert, qpert)
ELSE

!      Mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.

! DEPENDS ON: mix_inc
  CALL mix_inc(n_cg,n_cg,n_cg,nlev,nbl,ntml,                      &
             dthbydt,dqbydt,dubydt,dvbydt,l_tracer,ntra,dtrabydt, &
             p_layer_boundaries,p_layer_centres,index1)

END IF
!-----------------------------------------------------------------------
! 9.0 Moisture and Energy correction calculation 
!     appropriate for the new dynamics grid and vertical coordinate.
!     UM documentation paper 27 - section 12.
!-----------------------------------------------------------------------

IF (l_cv_conserve_check) THEN

! DEPENDS ON: cor_engy_5a
  CALL cor_engy_5a(n_cg,n_cg,nlev,index1, r_theta, r_rho                     &
                    ,r2rho_th, r2rho, dr_across_th, dr_across_rh             &
                    ,exner_layer_centres, th, u, v                           &
                    ,dubydt, dvbydt, dqclbydt, dqcfbydt                      &
                    ,rain,snow,dqbydt,dthbydt)
END IF
!-----------------------------------------------------------------------
! 9.0  Calculate convective cloud amount on model levels - no anvils
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------------
! 9.1 CCRad - Calculate CCA fow shallow levels only
!-----------------------------------------------------------------------------

IF (l_ccrad) THEN

  DO i=1, n_cg

    overlap_fac(i) = 0.0

    IF (iccb(i) /= 0) THEN ! Shallow convection occured

    !---------------------------------------------------------------
    ! Grant and Lock (2004) LES show mb/wsc nicely scales the cloud
    ! fraction profiles but not the TCA.  Also the UM overlap
    ! assumption in radiation is maximal.  This implies significant
    ! underestimate of TCA. So include a further parametrization of
    ! Cu "overlap", based again on the LES of Grant and Lock (2004).
    ! This increases cca_2d proportional to the ratio of the cloud
    ! to sub-cloud layer depths.  In order to preserve the grid-box
    ! cloud water, ccw will be divided by the same factor.
    !---------------------------------------------------------------

      overlap_fac(i) = 2.0*zcld(i) / z_rho(i,ntml(i)+1)

    END IF     ! iccb
  END DO     ! n_cg

  !---------------------------------------------------------------
  ! 9.11 Calculate CCA   - use shallow values
  !---------------------------------------------------------------
  SELECT CASE (cca2d_sh_opt)
    CASE(grant_lock)

      DO i=1, n_cg
        IF (iccb(i) /= 0) THEN ! Shallow convection occured
          cca_2d(i) = 2.0*mb(i)/wsc(i)

          ! Grab lowest cca value before any tuning occurs
          ! This will overwrite lcca in ni_conv_ctl only if neither
          ! shallow or deep have occured.  This is under a switch in
          ! the 4a scheme.
          !
          ! NOTE: Downdraughts & Evaporation still being fed cca_2d
          !       derived from TCW, This issue may require further
          !       investigation.
          lcca(i) = cca_2d(i)

        END IF     ! iccb
      END DO     ! n_cg

    CASE(total_condensed_water)
        ! cca_2d is left unchanged from that calculated in the
        ! code, which is based on TCW (Total Condensed Water)
        ! (TCW is a rate)

  END SELECT


  DO i=1, n_cg

    overlap_fac(i) = MAX( 0.5, overlap_fac(i) )
    overlap_fac(i) = MIN( 5.0, overlap_fac(i) )

    IF (overlap_fac(i)*cca_2d(i) > 0.99) THEN
      overlap_fac(i) = 0.99/cca_2d(i)
    END IF
  END DO      ! i (n_sh)


  !-------------------------------------------------------------------
  ! 9.12 Fill cca with cca_2d where non-zero ccw
  !-------------------------------------------------------------------
  DO k=1, nlev
    DO i=1, n_cg
      IF (iccb(i) /= 0) THEN ! Shallow convection occured

        IF (ccw(i,k) > 0.0) THEN
          zpr = (z_rho(i,k) - z_rho(i,ntml(i)+1)) / zcld(i)

          ! Apply Shape-function
          !
          ! Apply overlap_fac to cca, also preserving grid-box water
          ! by dividing ccw by overlap_fac, at least at cloud-base

          ccw(i,k)  = ccw(i,k)/overlap_fac(i)
          zpr       = MIN(1.0,zpr)
          cca(i,k)  = overlap_fac(i)*cca_2d(i)                    &
                         * 0.25*( 1.0 + 3.0*EXP(-5.0*zpr) )

        END IF       ! ccw
      END IF       ! iccb
    END DO       ! i (n_cg)
  END DO       ! k (nlev)

ELSE        ! Non CCRAD option

  DO k=1, n_cca_lev
    DO i=1, n_cg
      IF (k >= iccb(i) .AND. k < icct(i)) THEN
        cca(i,k) = cca_2d(i)
      END IF
    END DO
  END DO


END IF      ! l_ccrad


!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONGEST_CONV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE congest_conv
