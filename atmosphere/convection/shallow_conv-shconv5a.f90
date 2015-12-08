! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Shallow convection scheme
!
SUBROUTINE shallow_conv_5a(nbl,nlev,ntra,n_cca_lev,n_sh,trlev,    &
                       bland,delthvu,exner_layer_centres,         &
                       exner_layer_boundaries, l_calc_dxek,       &
                       l_q_interact, l_tracer, ntml, ntpar,       &
                       pstar,p_layer_centres,                     &
                       p_layer_boundaries,z_theta,z_rho,          &
                       r_theta,r_rho,r2rho_th,r2rho,              &
                       dr_across_th,dr_across_rh,                 &
                       q,q1_sd,t1_sd,th,timestep,u,v,w,uw0,vw0,   &
                       wstar,wthvs,entrain_coef,delta_smag,       &
                       zlcl_uv,ztop_uv,freeze_lev,recip_pstar,qse,&
                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
                       qcl,tracer,cape_out,cclwp,ccw,cca,         &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                       dqclbydt,dthbydt,                          &
                       dubydt,dvbydt,dtrabydt,                    &
                       detrain_up,detrain_dwn,                    &
                       entrain_up,entrain_dwn,                    &
                       iccb,icct,lcca,                            &
                       lcbase,lctop,rain,snow,                    &
                       rain_3d, snow_3d, up_flux, up_flux_half,   &
                       dwn_flux,uw_shall,vw_shall,tcw,cca_2d,     &
                       kterm, ind_shall,                          &
                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out  &
                       )

! Purpose:
!   Shallow convection scheme - works on points diagnosed as shallow in
!   subroutine CONV_DIAG.
!
!   Called by GLUE_CONV.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP 3 v8.3 programming standards

USE atmos_constants_mod, ONLY:                                    &
    cp, r, kappa, pref, repsilon, c_virtual


USE water_constants_mod, ONLY: lc, lf, tm

USE cv_run_mod, ONLY:                                             &
    l_mom, l_safe_conv, l_cv_conserve_check,                      &
    sh_pert_opt, bl_cnv_mix, icvdiag,                             &
    cca2d_sh_opt, cca_sh_knob, ccw_sh_knob, limit_pert_opt,       &
    cnv_wat_load_opt, l_ccrad

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, grant_lock,                            &
    thpixs_shallow, qpixs_shallow, c_mass,                        &
    max_sh_thpert, min_sh_thpert, max_sh_qpert_fac, beta_cu

USE cv_dependent_switch_mod, ONLY:                                &
    sh_on, mdet_sh_on, sh_ent_on, sh_sdet_on, sh_new_termc, sh_grey

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                     &
    flg_entr_up, flg_detr_up, flg_entr_dwn, flg_detr_dwn,         &
    flg_uw_shall, flg_vw_shall, flg_mf_shall

USE bl_option_mod, ONLY:                                          &
    Kprof_cu, off

USE earth_constants_mod, ONLY: g

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
 ,n_sh                 & ! No. of shallow convection points
 ,trlev                  ! No. of model levels on which tracers are included

LOGICAL, INTENT(IN) :: bland(n_sh) ! Land/sea mask

REAL, INTENT(IN)    :: delthvu(n_sh) !Integral of undilute parcel
                                     ! buoyancy over convective cloud
                                     ! layer (Kelvin m)

REAL, INTENT(IN)    ::               &
  exner_layer_centres(n_sh,0:nlev)   & ! Exner
 ,exner_layer_boundaries(n_sh,0:nlev)  ! Exner at half level above
                                       ! exner_layer_centres

LOGICAL, INTENT(IN) :: &
  l_calc_dxek          & ! Switch for calculation of condensate increment
 ,l_q_interact         & ! Switch allows overwriting parcel variables when
                         ! calculating condensate incr.
 ,l_tracer               ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: &
  ntml(n_sh)           & ! Top level of surface mixed layer defined relative to
                         ! theta,q grid
 ,ntpar(n_sh)            ! Top level of initial parcel ascent in BL scheme 
                         ! defined relative to theta,q grid

REAL, INTENT(IN)    ::          & 
  pstar(n_sh)                   & ! Surface pressure (Pa)
 ,p_layer_centres(n_sh,0:nlev)  & ! Pressure (Pa)
 ,p_layer_boundaries(n_sh,0:nlev) ! Pressure at half level above
                                  ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

REAL, INTENT(IN)    ::      &
  z_theta(n_sh,nlev)        & ! height of theta levels (m)
 ,z_rho(n_sh,nlev)          & ! height of rho levels (m)
 ,r_theta(n_sh,0:nlev)      & ! r on theta levels (m)
 ,r_rho(n_sh,nlev)          & ! r on rho levels (m)
 ,r2rho_th(n_sh,nlev)       & ! r2*rho theta levels (kg/m)
 ,r2rho(n_sh,nlev)          & ! r2*rho rho levels (kg/m)
 ,dr_across_th(n_sh,nlev)   & ! thickness of theta levels (m)
 ,dr_across_rh(n_sh,nlev)     ! thickness of rho levels (m)

REAL, INTENT(IN)    :: q(n_sh,nlev) ! Model mixing ratio (kg/kg)

REAL, INTENT(IN)    :: &
  q1_sd(n_sh)          & ! Standard deviation of turbulent flucts. of layer 1 q
                         ! (kg/kg)
 ,t1_sd(n_sh)            ! Standard deviation of turbulent flucts. of layer 1
                         ! temp. (K)

REAL, INTENT(IN)    :: th(n_sh,nlev) ! Model potential temperature (K)

REAL, INTENT(IN)    :: timestep    ! Model timestep (s)

REAL, INTENT(IN)    :: &
  u(n_sh,nlev)         & ! Model U field (m/s)
 ,v(n_sh,nlev)         & ! Model V field (m/s)
 ,w(n_sh,nlev)           ! Model W field (m/s)

REAL, INTENT(IN)    :: &
  uw0(n_sh)            & ! U-comp of surface stress (N/m2)
 ,vw0(n_sh)            & ! V-comp of surface stress (N/m2)
 ,wstar(n_sh)          & ! Convective velocity scale (m/s)
 ,wthvs(n_sh)          & ! Surface flux of THV (Pa m/s2)
 ,entrain_coef(n_sh)   & ! Entrainment coefficients
 ,delta_smag(n_sh)     & ! grid size (m)
 ,zlcl_uv(n_sh)        & ! Lifting condensation level defined for the uv 
                         ! grid (m)
 ,ztop_uv(n_sh)          ! Top of cloud layer defined for the uv grid (m)

INTEGER, INTENT(IN) :: freeze_lev(n_sh) ! Level index for freezing level

REAL, INTENT(IN) ::    &
  recip_pstar(n_sh)    & ! Reciprocal of pstar array
 ,qse(n_sh,nlev)         ! Saturation mixing ratio of cloud environment (kg/kg)

! Arguments with intent INOUT:

REAL, INTENT(INOUT) ::   &
  bulk_cf(n_sh,nlev)     & ! Bulk total cloud volume ( )
 ,cf_frozen(n_sh,nlev)   & ! Frozen water cloud volume ( )
 ,cf_liquid(n_sh,nlev)   & ! Liq water cloud volume ( )
 ,qcf(n_sh,nlev)         & ! Ice condensate mix ratio (kg/kg)
 ,qcl(n_sh,nlev)         & ! Liq condensate mix ratio (kg/kg)
 ,tracer(n_sh,trlev,ntra)  ! Model tracer fields (kg/kg)

! Arguments with intent OUT:

REAL, INTENT(OUT) :: &
  cape_out(n_sh)     & ! Saved convective available potential energy for
                       ! diagnostic output (J/kg)
 ,cclwp(n_sh)        & ! Condensed water path (kg/m2)
 ,ccw(n_sh,nlev)     & ! Convective cloud liquid water on model levels (kg/kg)
 ,cca(n_sh,n_cca_lev)  ! Convective cloud amount on model levels (fraction)

REAL, INTENT(OUT) ::   &
  dbcfbydt(n_sh,nlev)  & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_sh,nlev)  & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_sh,nlev)  & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_sh,nlev)    & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_sh,nlev)  & ! Increments to ice condensate due to convection
                         ! (kg/kg/s)
 ,dqclbydt(n_sh,nlev)  & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(n_sh,nlev)   & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_sh,nlev+1)  & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_sh,nlev+1)    ! Increments to V due to CMT (m/s2)

REAL, INTENT(OUT) ::      &
  dtrabydt(n_sh,nlev,ntra)  ! Increment to tracer due to convection (kg/kg/s)


REAL, INTENT(OUT) ::     &
  detrain_up(n_sh,nlev)  & ! Fractional detrainment rate into updraughts (Pa/s)
 ,detrain_dwn(n_sh,nlev) & ! Fractional detrainment rate into downdraughts 
                           ! (Pa/s)
 ,entrain_up(n_sh,nlev)  & ! Fractional entrainment rate into updraughts (Pa/s)
 ,entrain_dwn(n_sh,nlev)   ! Fractional entrainment rate into downdraughts
                           ! (Pa/s)

INTEGER, INTENT(OUT) :: &
  iccb(n_sh)            & ! Convective cloud base level
 ,icct(n_sh)              ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(n_sh) ! Lowest conv. cloud amt. (%)

INTEGER, INTENT(OUT) :: &
  lcbase(n_sh)          & ! Lowest conv. cloud base level 
 ,lctop(n_sh)             ! Lowest conv. cloud top level 

REAL, INTENT(OUT) :: rain(n_sh) ! Surface convective rainfall (kg/m2/s)

REAL, INTENT(OUT) :: snow(n_sh) ! Surface convective snowfall (kg/m2/s)

REAL, INTENT(OUT) :: rain_3d(n_sh,nlev) ! Convective rainfall flux (kg/m2/s)

REAL, INTENT(OUT) :: snow_3d(n_sh,nlev) ! Convective snowfall flux (kg/m2/s)

REAL, INTENT(OUT) ::       &   
  up_flux(n_sh,nlev)       & ! Updraught mass flux (Pa/s)
 ,up_flux_half(n_sh,nlev)  & ! Updraught mass flux (Pa/s)
 ,dwn_flux(n_sh,nlev)        ! Downdraught mass flux (Pa/s)

REAL, INTENT(OUT) :: uw_shall(n_sh,nlev) ! X-comp. of stress
                                         ! from shallow convection (kg/m/s2)

REAL, INTENT(OUT) :: vw_shall(n_sh,nlev) ! Y-comp. of stress
                                         ! from shallow convection (kg/m/s2)

REAL, INTENT(OUT) :: tcw(n_sh)  ! Total condensed water(kg/m2/s)

REAL, INTENT(OUT) :: cca_2d(n_sh) ! 2D convective cloud amount (%)

INTEGER, INTENT(OUT) :: kterm(n_sh) ! termination level for shallow
                                    ! convection
REAL, INTENT(OUT) ::      &
  ind_shall(n_sh)            ! 1.0 if real shallow convection else 0.0

! Adaptive detrainment output variables

REAL, INTENT(OUT) ::      &
  rbuoy_p_out(n_sh,nlev)  & ! buoyancy excess
 ,the_out(n_sh,nlev)      & ! th_E in parcel routine
 ,thp_out(n_sh,nlev)      & ! th_P in parcel routine
 ,qe_out(n_sh,nlev)       & ! q_E in parcel routine
 ,qp_out(n_sh,nlev)         ! q_P in parcel routine

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Adaptive detrainment output variables

REAL ::               &
  rbuoy_p_here(n_sh)  & ! buoyancy excess
 ,the_here(n_sh)      & ! th_E in parcel routine
 ,thp_here(n_sh)      & ! th_P in parcel routine
 ,qe_here(n_sh)       & ! q_E in parcel routine
 ,qp_here(n_sh)       & ! q_P in parcel routine
 ,rbuoy_p_old(n_sh)     ! buoyancy excess from previous k

REAL :: zk(n_sh)                ! Heights for use in calc

REAL :: zkp12(n_sh)             ! of moist static energy

REAL :: zkp1(n_sh)

INTEGER :: index1(n_sh),index2(n_sh)

INTEGER :: ncposs               ! No. of points which may convect

INTEGER :: nconv                ! No. of convecting points

REAL :: amdetk(n_sh)            ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

REAL :: b_calc                  ! Coefficient in thpert calc.

REAL :: c_calc                  ! Coefficient in thpert calc.

REAL :: cape(n_sh)              ! Convective available potential
                                ! energy (J/kg)

REAL :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

REAL :: dcpbydt(n_sh)           ! Rate of change of cape (J/kg/s)

REAL :: depth(n_sh)             ! Depth of convective cloud (m)

REAL :: delexkp1(n_sh)          ! Difference in exner ratio
                                ! across layer k+1

REAL :: dqsthk(n_sh)            ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k (kg/kg/K)

REAL :: dqsthkp1(n_sh)          ! Gradient of saturation mixing
                                ! ratio of cloud environment with
                                ! theta in layer k+1 (kg/kg/K)


REAL ::               &
  ekp14(n_sh)         & ! Entrainment coefficients at level k+1/4 multiplied
                        ! by appropriate layer thickness (dimensionless)
 ,ekp34(n_sh)         & ! Entrainment coefficients at level k+3/4 multiplied
                        ! by appropriate layer thickness (dimensionless)
 ,ekm14(n_sh)           ! Entrainment coefficients at level k-1+1/4 multiplied
                        ! by appropriate layer thickness (dimensionless)

REAL :: exk(n_sh)               ! Exner ratio at layer k

REAL :: exkp1(n_sh)             ! Exner ratio at layer k+1

REAL ::                  & 
  flxmax(n_sh)           & ! Maximum initial convective mass flux (Pa/s)
 ,flx_init(n_sh)         & ! Initial mass flux at cloud base (Pa/s)
 ,flx_init_new(n_sh)     & ! flx_init scaled (Pa/s)
 ,flx_init_term(n_sh)    & ! flx_init copy on termination
 ,flxmax_init(n_sh)      & ! Maximum possible initial mass flux (limited to 
                           ! the mass inthe initial convecting layer in Pa/s)
 ,max_cfl(n_sh)            ! Max cfl ratio over a convecting layer

REAL :: p_lcl(n_sh)             ! Pressure at LCL (Pa)

REAL :: precip(n_sh,nlev)       ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

REAL :: pk(n_sh)                ! Pressure at midpoint of layer
                                ! k (Pa)

REAL :: pkp1(n_sh)              ! Pressure at midpoint of layer
                                ! k+1 (Pa)

REAL :: delpk(n_sh)             ! Pressure difference over layer
                                ! k (Pa)

REAL :: delpkp1(n_sh)           ! Pressure difference over layer
                                ! k+1 (Pa)

REAL :: delpkp12(n_sh)          ! Pressure difference between
                                ! layers k and k+1 (Pa)

REAL :: delp_uv_k(n_sh)         ! Pressure difference across uv
                                ! layer k (Pa)

REAL :: delp_uv_kp1(n_sh)       ! Pressure difference across uv
                                ! layer k+1 (Pa)

REAL :: q_lcl(n_sh)             ! Mixing ratio at LCL (kg/kg)

REAL :: qse_lcl(n_sh)           ! Saturated q at LCL (kg/kg)

REAL :: rhum(n_sh)              ! Dummy relative humidity
                                ! (only used on shallow points)

REAL :: t_lcl(n_sh)             ! Temperature at LCL (K)

REAL :: th_lcl(n_sh)            ! Theta at LCL (K)

REAL :: dthv_ma                 ! Moist adiabtic change in thv
                                ! from ntml to ntml+1 (K)

REAL :: thv_pert(n_sh)          ! Theta_v parcel pertubation (K)

REAL :: thpert(n_sh)            ! Theta parcel pertubation (K)

REAL :: qpert(n_sh)             ! q parcel pertubation (kg/kg)

REAL :: rho_k                   ! density on level k


INTEGER :: start_lev3c(n_sh)    ! Compressed convection
                                ! initiation level

REAL :: wsc(n_sh)               ! Convective velocity scale (m/s)

REAL :: wsc_o_mb(n_sh)          ! Convective velocity scale /mb

LOGICAL :: l_shallow(n_sh)      ! Dummy variable (=.T.)

LOGICAL :: l_mid(n_sh)          ! Dummy variable (=.F.)

LOGICAL :: cumulus(n_sh)        ! Dummy variable (=.T.)

LOGICAL ::             &
  bgmk(n_sh)           & ! Mask for points where parcel in layer k is saturated
 ,bgmk_term(n_sh)      & ! Mask at termination
 ,bwater(n_sh,2:nlev)  & ! Mask for points at which condensate is liquid
 ,bwk(n_sh)            & ! Mask for liquid condensate on k
 ,bwkp1(n_sh)          & ! Mask for liquid condensate on k+1
 ,blowst(n_sh)         & ! Dummy variable indicating low enough stability
                         ! for convection to occur
 ,bterm(n_sh)          & ! Mask for points which have stopped convecting
 ,bconv(n_sh)          & ! Mask for points at which convection is occurring
 ,bcposs(n_sh)           ! Mask for points passing initial stability test

! Parcel variables

REAL :: qpi(n_sh)               ! Initial parcel mixing ratio
                                !(kg/kg)

REAL :: qp(n_sh,nlev)           ! Parcel mixing ratio (kg/kg)

REAL :: thpi(n_sh)              ! Initial parcel potential temp.
                                !(K)

REAL :: thp(n_sh,nlev)          ! Parcel potential temp (K)

REAL :: up(n_sh,nlev)           ! Parcel U (m/s)

REAL :: vp(n_sh,nlev)           ! Parcel V  (m/s)

REAL :: trap(n_sh,nlev,ntra)    ! Tracer content of parcel
                                ! (kg/kg)

REAL :: expi(n_sh)              ! Initial parcel exner pressure

REAL :: xpk(n_sh,nlev)          ! Parcel cloud water (kg/kg)

REAL :: flx(n_sh,nlev)          ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(n_sh,nlev)     ! Minmum parcel buoyancy excess

REAL :: thpixs_v(n_sh,nlev)     ! Theta parcel excess (K)

REAL :: qpixs_v(n_sh,nlev)      ! Q parcel excess(kg/kg)


! PC2

REAL :: qclp(n_sh,nlev)         ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

REAL :: qcfp(n_sh,nlev)         ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)


! Parameters

REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed


! CMT variables

INTEGER :: nlcl_uv(n_sh)        ! Level index for LCL

INTEGER :: ntop_uv(n_sh)        ! Level index for top of layer

INTEGER :: n_0degc(n_sh)        ! Level index for zero degrees

INTEGER :: cu_term(n_sh),cu_tend(n_sh) !Indicies for CMT subs

REAL :: exk_temp                ! Temporary exner

REAL :: eflux_u_ud(n_sh)        ! Vertical eddy flux of momentum
                                ! due to UD at top of layer
                                ! (Pa m/s2)

REAL :: eflux_v_ud(n_sh)        ! Vertical eddy flux of momentum
                                ! due to UD at bottom of layer
                                ! (Pa m/s2)

REAL :: flxkp12(n_sh,nlev)      ! Mass flux on half level (Pa/s)

REAL :: mb(n_sh)                ! Cloud base mass flux (Pa/s)

REAL :: p_uv(nlev,n_sh)         ! Pressure of model level (Pa)

REAL :: phalf_uv(nlev,n_sh)     ! Pressure of half level (Pa)

REAL :: plcl_uv(n_sh)           ! Pressure at LCL (Pa)

REAL :: ptop_uv(n_sh)           ! Pressure at top of cloud layer
                                ! (Pa)

REAL :: p_0degc_uv(n_sh)        ! Pressure of zero degree level
                                ! (Pa)

REAL :: rho_uv(nlev,n_sh)       ! Density on uv level (kg/m3)

REAL :: uw(nlev,n_sh)           ! U- comp stress profile (N/m2)
                                ! (units change through calls)

REAL :: ue_p(nlev,n_sh)         ! Environment U profile (m/s)

REAL :: vw(nlev,n_sh)           ! V-comp stress profile (N/m2)

REAL :: ve_p(nlev,n_sh)         ! Environment V profile (m/s)

REAL :: zcld(n_sh)              ! Depth of cloud layer (m)

LOGICAL :: l_mom_gk             ! true if Gregory-Kershaw CMT


! CFL scaling variables


INTEGER :: det_lev(n_sh)        ! Level at which split final
                                ! detrainment last occurred

INTEGER :: nterm                ! No. of points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage

REAL :: weight_param            ! Weighting factor

REAL :: scale_f(n_sh)           ! store scaling factor

REAL :: cca_2d_term(n_sh)       ! store 2d CCA on termination

! Downdraught scheme variables


INTEGER :: nnodd                ! No. of downdraughts not possible

INTEGER :: index_nodd(n_sh)     ! Index of downdraughts not
                                ! possible
INTEGER :: npossdd              ! No. downdraughts possible

INTEGER :: index_possdd(n_sh)   ! Index of downdraughts

INTEGER :: kmax_term            ! maximum termination level + 1

REAL :: deltap_cld              ! pressure thickness of convective
                                ! cloud (Pa)


! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag

INTEGER :: ntpar_max          ! max ntpar value


! parameters ect for qmin checks

REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

REAL :: qminincolumn(n_sh)     ! Minimum value for q in column (kg/kg)

REAL :: temp1(n_sh)            ! work array

! Local compressed arrays


LOGICAL :: bconv_c2(n_sh)

LOGICAL :: bgmkp1_c(n_sh), bgmkp1_c2(n_sh) ! Mask for points
                                ! where parcel in layer k+1
                                ! is saturated

LOGICAL :: bwk_c(n_sh), bwk_c2(n_sh) ! bwater mask in layer k

LOGICAL :: bwkp1_c(n_sh), bwkp1_c2(n_sh) ! bwater mask in layer k+1

REAL :: deltak_c2(n_sh)         ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

REAL :: dqek_c2(n_sh)           ! Increment to q due to
                                ! convection in layer k (kg/kg)

REAL :: dqekp1_c2(n_sh)         ! Increment to q due to
                                ! convection in layer k+1 (kg/kg)

REAL :: dthek_c2(n_sh)          ! Increment to potential temp.
                                ! due to convection in layer k

REAL :: dthekp1_c2(n_sh)        ! Increment to potential temp.
                                ! due to convection in layer k+1

REAL :: dtraek_c2(n_sh,ntra)    ! Increment to model tracer due
                                ! to conv. at level k (kg/kg/s)

REAL :: dtraekp1_c2(n_sh,ntra)  ! Increment to model tracer due
                                ! to conv. at level k+1 (kg/kg/s)

REAL :: duek_c2(n_sh)           ! Increment to model U in layer k
                                ! due to CMT (m/s2)

REAL :: duekp1_c2(n_sh)         ! Increment to model U in layer
                                ! k+1 due to CMT (m/s2)

REAL :: dvek_c2(n_sh)           ! Increment to model V in layer k

REAL :: dvekp1_c2(n_sh)         ! Increment to model V in layer
                                ! k+1 due to CMT (m/s2)

REAL :: flxk_c(n_sh), flxk_c2(n_sh) !Parcel mass flux in layer k
                                ! (Pa/s)

REAL :: flxkp12_c2(n_sh)        ! Half level mass flux (Pa/s)

REAL :: prekp1_c2(n_sh)         ! Precip. from parcel as it rises
                                ! from layer k to k+1 (kg/m2/s)

REAL :: qpk_c(n_sh), qpk_c2(n_sh) ! Parcel mixing ratio in
                                ! layer k(kg/kg)
REAL :: qpk(n_sh)               !ad. entrain.

REAL :: qpkp1_c(n_sh), qpkp1_c2(n_sh) ! Parcel mixing ratio
                                ! in layer k+1 (kg/kg)

REAL :: qek_c(n_sh), qek_c2(n_sh) ! Env. mixing ratio in
                                ! layer k (kg/kg)
REAL :: qek(n_sh)               !ad. entrain.

REAL :: qekp1_c(n_sh), qekp1_c2(n_sh) ! Env. mixing ratio in
                                ! layer k+1 (kgkg-1)
REAL :: qekp1(n_sh)               !ad. entrain.

REAL :: qsek_c2(n_sh)           ! Saturation mixing ratio of
                                ! cld. env. in layer k (kg/kg)
REAL :: qsek(n_sh)              !ad. entrain.

REAL :: qsekp1_c(n_sh), qsekp1_c2(n_sh) ! Saturation mixing ratio
                                ! of cld. env. in layer k+1
                                ! (kg/kg)
REAL :: qsekp1(n_sh)            !ad. entrain.

REAL :: thek_c(n_sh), thek_c2(n_sh) ! Env. potential temp
                                ! in layer k (K)
REAL :: thek(n_sh)              !ad. entrain.

REAL :: thekp1_c(n_sh), thekp1_c2(n_sh) ! Env. potential temp i
                                ! in layer k (K)
REAL :: thekp1(n_sh)            !ad. entrain.

REAL :: thpk_c(n_sh), thpk_c2(n_sh) ! Parcel potential temp
                                ! in layer k (K)
REAL :: thpk(n_sh)              !ad. entrain.

REAL :: thpkp1_c(n_sh), thpkp1_c2(n_sh)! Parcel potential temp
                                ! in layer k (K)

REAL :: traek_c(n_sh,ntra), traek_c2(n_sh,ntra) ! Tracer content
                                ! cld. env. in layer k (kgkg-1)

REAL :: traekp1_c(n_sh,ntra), traekp1_c2(n_sh,ntra) ! Tracer
                                ! content of cloud env.
                                ! in layer k+1 (kg/kg)

REAL :: trapk_c(n_sh,ntra), trapk_c2(n_sh,ntra) ! Tracer cont.
                                ! of parcel in layer k (kg/kg)

REAL :: trapkp1_c(n_sh,ntra), trapkp1_c2(n_sh,ntra) ! Tracer cont.
                                ! of parcel in layer k+1 (kg/kg)

REAL :: rbuoy_c(n_sh), rbuoy_c2(n_sh) ! Buoyancy of parcel at k+1 (Kelvin)

REAL :: uek_c(n_sh), uek_c2(n_sh) ! Model U field on layer k (m/s)

REAL :: uekp1_c(n_sh), uekp1_c2(n_sh)! Model U field on layer k+1 (m/s)

REAL :: vek_c(n_sh), vek_c2(n_sh) ! Model V field on layer k (m/s)

REAL :: vekp1_c(n_sh), vekp1_c2(n_sh) ! Model V field on layer k+1 (m/s)

REAL :: upk_c(n_sh), upk_c2(n_sh) ! Parcel U in layer k
                                  ! after entrainment (m/s)


REAL :: upkp1_c(n_sh), upkp1_c2(n_sh) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

REAL :: vpk_c(n_sh), vpk_c2(n_sh) ! Parcel V in layer k
                                  ! after entrainment (m/s)

REAL :: vpkp1_c(n_sh), vpkp1_c2(n_sh) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

REAL :: xsqkp1_c(n_sh), xsqkp1_c2(n_sh) ! Excess water vapour
                                        ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

REAL :: qclek_c(n_sh), qclek_c2(n_sh) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qclekp1_c(n_sh), qclekp1_c2(n_sh) ! Environment liquid
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qcfek_c(n_sh), qcfek_c2(n_sh) ! Environment frozen
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qcfekp1_c(n_sh), qcfekp1_c2(n_sh) ! Environment frozen
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qclpk_c(n_sh), qclpk_c2(n_sh) ! Parcel liquid
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qclpkp1_c(n_sh), qclpkp1_c2(n_sh) ! Parcel liquid
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: qcfpk_c(n_sh), qcfpk_c2(n_sh) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qcfpkp1_c(n_sh), qcfpkp1_c2(n_sh) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: cflek_c2(n_sh),cflekp1_c2(n_sh) ! Environment liquid water
                                ! cloud volume ( )

REAL :: cffek_c2(n_sh),cffekp1_c2(n_sh) ! Environment frozen water
                                ! cloud volume ( )

REAL :: bcfek_c2(n_sh),bcfekp1_c2(n_sh) ! Environment bulk total
                                ! cloud volume ( )

REAL :: dqclek_c2(n_sh),dqclekp1_c2(n_sh) ! Environment increments
                                ! to liquid condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dqcfek_c2(n_sh),dqcfekp1_c2(n_sh) ! Environment increments
                                ! to frozen condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dcflek_c2(n_sh),dcflekp1_c2(n_sh) ! Environment increments
                                ! to liquid water cloud volume due
                                ! to convection (/s)

REAL :: dcffek_c2(n_sh),dcffekp1_c2(n_sh) ! Environment increments
                                ! to frozen water cloud volume due
                                ! to convection (/s)

REAL :: dbcfek_c2(n_sh),dbcfekp1_c2(n_sh) ! Environment increments
                                ! to bulk total cloud volume due
                                ! to convection (/s)

REAL :: amdetk_c2(n_sh)
LOGICAL :: bgmk_c2(n_sh)
LOGICAL :: bland_c2(n_sh)
LOGICAL :: blowst_c2(n_sh)
LOGICAL :: bterm_c2(n_sh)
REAL :: cape_c2(n_sh)
REAL :: cca_2d_c2(n_sh)
REAL :: cclwp_c2(n_sh)
REAL :: ccw_c2(n_sh)
LOGICAL :: cumulus_c(n_sh), cumulus_c2(n_sh)
REAL :: dcpbydt_c2(n_sh)
REAL :: delexkp1_c2(n_sh)
REAL :: delpk_c2(n_sh)
REAL :: delpkp1_c2(n_sh)
REAL :: delp_uv_k_c2(n_sh)
REAL :: delp_uv_kp1_c2(n_sh)
REAL :: depth_c2(n_sh)
REAL :: dptot_c2(n_sh)
REAL :: dqsthkp1_c2(n_sh)
REAL :: dqsthk_c2(n_sh)
REAL :: eflux_u_ud_c2(n_sh)
REAL :: eflux_v_ud_c2(n_sh)
REAL :: ekp14_c(n_sh),ekp14_c2(n_sh)
REAL :: ekp34_c(n_sh),ekp34_c2(n_sh)
REAL :: exk_c2(n_sh)
REAL :: exkp1_c(n_sh),exkp1_c2(n_sh)
REAL :: expi_c2(n_sh)
INTEGER :: icct_c2(n_sh)
INTEGER :: iccb_c2(n_sh)
INTEGER :: lctop_c2(n_sh)
INTEGER :: lcbase_c2(n_sh)
REAL :: lcca_c2(n_sh)
LOGICAL :: l_shallow_c2(n_sh)
LOGICAL :: l_mid_c2(n_sh)
REAL :: max_cfl_c2(n_sh)
REAL :: pk_c(n_sh),pk_c2(n_sh)
REAL :: pkp1_c(n_sh),pkp1_c2(n_sh)
REAL :: pstar_c2(n_sh)
REAL :: q1_sd_c2(n_sh)
REAL :: qpi_c2(n_sh)
REAL :: qpixs_v_c2(n_sh)
REAL :: relh_c2(n_sh)
REAL :: rbuoy_p_here_c2(n_sh)
REAL :: the_here_c2(n_sh)
REAL :: thp_here_c2(n_sh)
REAL :: qe_here_c2(n_sh)
REAL :: qp_here_c2(n_sh)
REAL :: rbuoy_p_old_c2(n_sh)
REAL :: tcw_c2(n_sh)
REAL :: thpi_c2(n_sh)
REAL :: thpixs_v_c2(n_sh)
REAL :: t1_sd_c2(n_sh)
REAL :: xpk_c(n_sh),xpk_c2(n_sh)
REAL :: xsbmin_v_c2(n_sh)
REAL :: qsat_lcl(n_sh)         ! not used
LOGICAL :: b_nodd(n_sh)   ! points with no downdraught
LOGICAL :: b_dd(n_sh)     ! points with downdraught on termination

INTEGER :: n_real_sh   ! real shallow ascents i.e. terminate above ntml

!===============================================================
! CCRad Variables local variables
!===============================================================

REAL   :: overlap_fac(n_sh)  ! Factor designed to improve
                             ! shallow Cu cover by allowing
                             ! for non-vertical clouds.

REAL   :: zpr         ! method (BL Fluxes) only if
                      !   l_ccrad = T .AND. cca2d_sh_opt  = 1

!===============================================================
! End CCRad Variables local variables
!===============================================================

! Loop counters

INTEGER :: i,i2,j,k,ktra,kt










INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle








IF (lhook) CALL dr_hook('SHALLOW_CONV_5A',zhook_in,zhook_handle)


!initialise SCM diagnostics

DO k = 1,nlev
  DO i = 1,n_sh
    rbuoy_p_out(i,k)=0.0
    the_out(i,k)=th(i,k)
    thp_out(i,k)=th(i,k)
    qe_out(i,k)=q(i,k)
    qp_out(i,k)=q(i,k)
  END DO
END DO

! Initialise logicals

DO i = 1,n_sh
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


!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------
! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)


DO k = 1,nlev-1
  DO i = 1,n_sh
    xsbmin_v(i,k) = MIN( ((p_layer_centres(i,k) -                 &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = MIN( ((p_layer_centres(i,k) -                 &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_shallow

    qpixs_v(i,k)  = qpixs_shallow
  END DO
END DO  ! nlev

! Calculate convective velocity scale and cloud base mass flux

DO i = 1,n_sh
  wsc(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i)) &
               * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
  mb(i)  = c_mass * wstar(i)
  zcld(i) = ztop_uv(i) - zlcl_uv(i)
  wsc_o_mb(i) = wsc(i)/mb(i)

  kterm(i) = 0
END DO

! Define the LCL

IF ( sh_pert_opt == 0) THEN

!         Define the LCL at the half level above ntml. Find
!         environmental T at p_lcl by approximating theta there with
!         th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is
!         tunable.  Similarly for q.

  DO i = 1,n_sh
    k =ntml(i)
    p_lcl(i)  = p_layer_boundaries(i,k)
    th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
    t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
    q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
  END DO

ELSE  ! Sh_pert_Opt = 1

!         Define the LCL at the half level above ntml. Find
!         environmental T at p_lcl by approximating theta there with
!         th(i,k) Similarly for q.

  DO i = 1,n_sh
    k =ntml(i)
    p_lcl(i)  = p_layer_boundaries(i,k)
    th_lcl(i) = th(i,k)
    t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
    q_lcl(i)  = q(i,k)
  END DO

END IF

! Calculate saturation mixing ratio at LCL

! DEPENDS ON: qsat_mix
CALL qsat_mix(qse_lcl,t_lcl,p_lcl,n_sh,.FALSE.)

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

l_mom_gk = .FALSE.       ! not Gregory-Kershaw CMT
                         ! Shallow code uses turbulence based CMT

IF (l_mom) THEN

! Initialize arrays required for Convective Momentum Transport(CMT)

  k=1
  DO i = 1,n_sh
    p_uv(k,i)     = p_layer_boundaries(i,k-1)
    phalf_uv(k,i) = p_layer_centres(i,k-1)
    ue_p(k,i)     = u(i,k)
    ve_p(k,i)     = v(i,k)
  END DO

  DO i = 1,n_sh
    nlcl_uv(i)    = ntml(i) + 1
    ntop_uv(i)    = ntpar(i) + 1
    n_0degc(i)    = freeze_lev(i)
  END DO

  DO i = 1,n_sh
    DO k = 2,nlev
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      exk_temp      = (p_uv(k,i)/pref)**kappa
      rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp * (th(i,k-1) + th(i,k)))
    END DO
    plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
    ptop_uv(i)      = phalf_uv(ntop_uv(i),i)
    p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
    rho_uv(1,i)     = rho_uv(2,i)
  END DO
END IF     !L_mom

! Calculate theta and q pertubations (pertubation is based on
! environment buoyancy gradient)
! Reset th and q xs's at ntml

IF ( sh_pert_opt == 0) THEN
  DO i = 1,n_sh

    k = ntml(i)

    IF (t_lcl(i) >  tm) THEN
      dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                      / (r * t_lcl(i) * t_lcl(i))

      ! Estimate of moist adiabatic lapse rate

      dthv_ma    = ( (lc/cp) - (1.+c_virtual)*th(i,k) )*          &
                   dq_sat_env*(g/cp)/(1.0+(lc/cp)*dq_sat_env)
    ELSE
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                      / (r * t_lcl(i) * t_lcl(i))

      ! Estimate of moist adiabatic lapse rate (in K/m)

      dthv_ma    = ( ((lc+lf)/cp) - (1.0+c_virtual)*th(i,k) )*    &
                   dq_sat_env*(g/cp)/(1.0+((lc+lf)/cp)*dq_sat_env)
    END IF

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                        + c_virtual * qse_lcl(i)

    ! Calculate theta_v perturbation:

    IF (Kprof_cu == off) THEN
      thv_pert(i) = -0.17 * wthvs(i) / mb(i)                      &
                    + (th(i,k+1) * (1.0 + c_virtual               &
                    * q(i,k+1)) - th(i,k)                         &
                    * (1.0 + c_virtual * q(i,k)))
    ELSE
      ! "entrainment flux" at LCL given by BL scheme
      thv_pert(i) =  (th(i,k+1) * (1.0 + c_virtual                 &
                    * q(i,k+1)) - th(i,k)                          &
                    * (1.0 + c_virtual * q(i,k)))
    END IF

    c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i))    &
                         - thv_pert(i)

    thpert(i) = -c_calc / b_calc  ! ignore term in thpert**2

    thpixs_v(i,k) = thpert(i)

    qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)                   &
                           **kappa) * thpert(i) * dq_sat_env      &
                           - q_lcl(i)

    qpixs_v(i,ntml(i))  = qpert(i)

  END DO !n_sh

ELSE ! Sh_pert_opt = 1
  IF (limit_pert_opt == 1 .OR. limit_pert_opt == 2) THEN
    DO i = 1,n_sh

      k = ntml(i)

      IF (t_lcl(i) >  tm) THEN
        dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                        / (r * t_lcl(i) * t_lcl(i))

!             Estimate of moist adiabatic lapse rate

        dthv_ma    = ( (lc/cp) - (1.+c_virtual)*th(i,k) )*          &
                     dq_sat_env*(g/cp)/(1.+(lc/cp)*dq_sat_env)
      ELSE
        dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                        / (r * t_lcl(i) * t_lcl(i))

!             Estimate of moist adiabatic lapse rate (in K/m)

        dthv_ma    = ( ((lc+lf)/cp) - (1.+c_virtual)*th(i,k) )*     &
                     dq_sat_env*(g/cp)/(1.+((lc+lf)/cp)*dq_sat_env)
      END IF

      b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                          + c_virtual * qse_lcl(i)


!           Calculate theta_v perturbation:

!           First convert moist adiabatic lapse rate to thv difference
!           between levels k and k+1

      rho_k = p_layer_centres(i,k) /                              &
             (r * th(i,k) * (p_layer_centres(i,k)/ pref)**kappa)
      dthv_ma = -dthv_ma*                                         &
             (p_layer_centres(i,k+1)-p_layer_centres(i,k)) /      &
             (rho_k*g)

!           Make perturbation relative to a target lapse rate (namely
!           0.6*dthv_ma, which is approximately what is seen in LES)
      IF (Kprof_cu == off) THEN
        thv_pert(i) = -0.17 * wthvs(i) / mb(i) +  0.6*dthv_ma     &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))     &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

      ELSE
        ! "entrainment flux" at LCL done by BL scheme
        thv_pert(i) = 0.6*dthv_ma                                 &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))     &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

      END IF

!           limit thv_pert to physically sensible values

      thv_pert(i) = MAX(MIN(thv_pert(i), max_sh_thpert),          &
                        min_sh_thpert)

      c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i))  &
               - thv_pert(i)

      thpert(i) = MAX(MIN(-c_calc / b_calc, max_sh_thpert),       &
                      min_sh_thpert)  ! ignore term in thpert**2

      thpixs_v(i,k) = thpert(i)

      qpert(i)  = MAX(MIN(qse_lcl(i) + ((p_lcl(i) / pref)         &
                        **kappa) * thpert(i) * dq_sat_env         &
                        - q_lcl(i),                               &
                        max_sh_qpert_fac * qse_lcl(i)),0.0)

      qpixs_v(i,ntml(i))  = qpert(i)

    END DO !n_sh

  ELSE IF (limit_pert_opt == 0) THEN
    DO i = 1,n_sh

      k = ntml(i)

      IF (t_lcl(i) >  tm) THEN
        dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                        / (r * t_lcl(i) * t_lcl(i))

!             Estimate of moist adiabatic lapse rate

        dthv_ma    = ( (lc/cp) - (1.0+c_virtual)*th(i,k) )*         &
                     dq_sat_env*(g/cp)/(1.0+(lc/cp)*dq_sat_env)
      ELSE
        dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                        / (r * t_lcl(i) * t_lcl(i))

!             Estimate of moist adiabatic lapse rate (in K/m)

        dthv_ma    = ( ((lc+lf)/cp) - (1.0+c_virtual)*th(i,k) )*    &
                     dq_sat_env*(g/cp)/(1.0+((lc+lf)/cp)*dq_sat_env)
      END IF

      b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                          + c_virtual * qse_lcl(i)


!           Calculate theta_v perturbation:

!           First convert moist adiabatic lapse rate to thv difference
!           between levels k and k+1

      rho_k = p_layer_centres(i,k) /                              &
             (r * th(i,k) * (p_layer_centres(i,k)/ pref)**kappa)
      dthv_ma = -dthv_ma*                                         &
             (p_layer_centres(i,k+1)-p_layer_centres(i,k)) /      &
             (rho_k*g)

!           Make perturbation relative to a target lapse rate (namely
!           0.6*dthv_ma, which is approximately what is seen in LES)

      IF (Kprof_cu == off) THEN
        thv_pert(i) = -0.17 * wthvs(i) / mb(i) +  0.6*dthv_ma     &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))     &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

      ELSE
        ! "entrainment flux" at LCL done by BL scheme
        thv_pert(i) = 0.6*dthv_ma                                 &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))     &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

      END IF

      c_calc   = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i))  &
                         - thv_pert(i)

      thpert(i) = -c_calc / b_calc  ! ignore term in thpert**2

      thpixs_v(i,k) = thpert(i)

      qpert(i)  = qse_lcl(i) + ((p_lcl(i) / pref)                 &
                             **kappa) * thpert(i) * dq_sat_env    &
                             - q_lcl(i)

      qpixs_v(i,ntml(i))  = qpert(i)

    END DO !n_sh

  END IF

END IF

! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)

! DEPENDS ON: flag_wet
CALL flag_wet(n_sh,n_sh,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

DO k = 1,nlev
  DO i = 1,n_sh
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
    flxkp12(i,k)  = 0.0
  END DO
END DO

IF (l_mom) THEN
  DO k = 1,nlev+1
    DO i = 1,n_sh
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

! No need to initialise dtrabydt as done in glue_conv

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------

IF (flg_up_flx .OR. flg_mf_shall) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      up_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      up_flux_half(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      dwn_flux(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      entrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      detrain_up(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
     entrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      detrain_dwn(i,k) = 0.0
    END DO
  END DO
END IF

IF (l_mom) THEN
  IF (flg_uw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        uw_shall(i,k) = 0.0
      END DO
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        vw_shall(i,k) = 0.0
      END DO
    END DO
  END IF
END IF  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
DO i = 1,n_sh
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
  DO i=1, n_sh
    cca(i,k) = 0.0
  END DO
END DO

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics - done in glue routine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling calculations
!-----------------------------------------------------------------------

DO i = 1,n_sh
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)      = 0.0
  cape_out(i)  = 0.0
  dcpbydt(i)   = 0.0
  max_cfl(i)   = 0.0
  det_lev(i)   = 0
  cca_2d_term(i) = 0.0
  flx_init_term(i) = 0.0

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

DO i = 1,n_sh
  rbuoy_p_here(i) = 0.0
  the_here(i) = th(i,1)
  thp_here(i) = th(i,1)
  qe_here(i) = q(i,1)
  qp_here(i) = q(i,1)
  rbuoy_p_old(i) = 0.0
END DO
!initialise ekm14
DO i =1, n_sh
  ekm14(i) =0.0
END DO
!Initialise adaptive entrainment variables
DO i = 1, n_sh
  thek(i)=th(i,1)
  qek(i)=q(i,1)
  qsek(i)=qse(i,1)
  thekp1(i)=th(i,2)
  qekp1(i)=q(i,2)
  qsekp1(i)=qse(i,2)
  thpk(i)=thp(i,1)
  qpk(i)=qp(i,1)
  bwk(i)=bwater(i,2)
  bwkp1(i)=bwater(i,2)
  !Note that unlike p_layer_boundaries, where k indexing is offset
  !by one compared to the dynamics numbering, z retains the numbering
  !convention for dynamics variables i.e. for theta levels, k->k
  !and for rho levels k+1/2 -> k+1
  zk(i) = z_theta(i,1)
  zkp12(i)=z_rho(i, 2)
  zkp1(i)=z_theta(i, 2)
END DO

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

IF (sh_on == 1 .OR. sh_grey == 1) THEN   
                      ! Adaptive forced detrainment or grey shallow param
                      ! No limit on convection top
  ntpar_max = nlev-3  ! What is a sensible value to have here?

ELSE                  ! Top limited

  ntpar_max=0
  DO i = 1,n_sh
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

! initialize SCM diagnostics for this pass through the loop
!NB do not re-initialise rbuoy_p_old 'cos value from k-1 needed
  DO i = 1,n_sh
    rbuoy_p_here(i) = 0.0
    the_here(i) = th(i,k)
    thp_here(i) = th(i,k)
    qe_here(i) = q(i,k)
    qp_here(i) = q(i,k)
    rbuoy_p_here_c2(i) = 0.0
    the_here_c2(i) = 0.0
    thp_here_c2(i) = 0.0
    qe_here_c2(i) = 0.0
    qp_here_c2(i) = 0.0
  END DO
  !Initialise adaptive entrainment variables
  DO i = 1, n_sh
    thek(i)=th(i,k)
    qek(i)=q(i,k)
    qsek(i)=qse(i,k)
    thekp1(i)=th(i,k+1)
    qekp1(i)=q(i,k+1)
    qsekp1(i)=qse(i,k+1)
    IF(k  ==  2) THEN      !set to environmental values
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

  ! Set relative humidity in layer k (rhum)

  DO i = 1,n_sh
    rhum(i) = q(i,k) / qse(i,k)
  END DO

!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
  CALL layer_cn(k,n_sh,nlev                                       &
,                   mdet_sh_on, sh_ent_on                         &
,                   ntml,ntpar                                    &
,                   .TRUE.,.FALSE.,.FALSE.                        &
,                   bconv,bwk,bwkp1                               &
,                   exner_layer_boundaries                        &
,                   exner_layer_centres                           &
,                   p_layer_boundaries,p_layer_centres            &
,                   recip_pstar,entrain_coef,rhum                 &
,                   zk, zkp12, zkp1                               &
,                   thek, qek,qsek, thekp1,qekp1,qsekp1           &
,                   thpk,qpk , qsat_lcl, ekm14                    &
,                   pkp1,delpkp1,exkp1                            &
,                   pk,delpk,delpkp12,exk,delexkp1                &
,                   delp_uv_k, delp_uv_kp1                        &
,                   ekp14,ekp34,amdetk)

! Set ekm14 for next pass through loop
  DO i = 1, n_sh
    ekm14(i) = ekp14(i)
  END DO

! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)

  IF (k == 2) THEN
! DEPENDS ON: dqs_dth
    CALL dqs_dth(dqsthk,k,th(1,k),qse(1,k),exk,n_sh)
  ELSE
    DO i = 1,n_sh
      dqsthk(i) = dqsthkp1(i)
    END DO
  END IF

! DEPENDS ON: dqs_dth
  CALL dqs_dth(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_sh)

! Set other grid dependent constants

  DO i = 1,n_sh

    ! Maximum initial convective mass flux

    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

  END DO

! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k

  DO i = 1,n_sh
! not convecting and not convected in column before
! PC2 qclp and qcfp zero at this point but will add an initial
! value at cloud base
    IF ( .NOT. bconv(i).AND.det_lev(i) == 0) THEN
      expi(i)  = exk(i)
      xpk(i,k) = 0.0
      qclp(i,k) = 0.0
      qcfp(i,k) = 0.0
      flx(i,k) = 0.0
      bgmk(i)       = .FALSE.
      bgmk_term(i)  = .FALSE.
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
  END DO  ! n_sh
  IF (l_tracer) THEN
    DO ktra=1,ntra
      DO i = 1,n_sh
        IF ( .NOT. bconv(i)) THEN
           trap(i,k,ktra)  = tracer(i,k,ktra)
        END IF  !not bconv
      END DO
    END DO
  END IF


! Scale entrainment coefficients with cloud base mass flux
! and convective velocity scale

  DO i = 1,n_sh

    IF (k  >=  ntml(i)) THEN
      ! If original entrainment coefficient then scale as before else leave
      ! as calculated by layer_cn.
      IF (entrain_coef(i) < 0.0 ) THEN
        ekp14(i)  = ekp14(i) * wsc_o_mb(i)
        ekp34(i)  = ekp34(i) * wsc_o_mb(i)
        amdetk(i) = amdetk(i) * wsc_o_mb(i)
      END IF
    END IF

    ! Carry out initial test to see if convection is possible from layer
    ! k to k+1. Set bcposs = .T. if
    ! 1. the point was convecting (bconv = .T.) and did not terminate
    ! in the previous layer  OR
    ! 2. k = ntml

    bcposs(i) = bconv(i) .OR. k  ==  ntml(i)

  END DO  ! n_sh

! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)

  ncposs = 0
  DO i = 1,n_sh
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
  CALL lift_par_5a(ncposs,n_sh,thpkp1_c,qpkp1_c,xsqkp1_c,         &
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

! Set parcel mass flux
! UM Documentation paper 27, section 1.5

      flxk_c(i) = mb(index1(i)) * g *                             &
                      p_layer_centres(index1(i),k) / (r *         &
                        thpk_c(i) * (p_layer_centres(index1(i),k) &
                          / pref)**kappa)

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

! Reset threashold for forced detrainment to the initial
! (potentially negative) buoyancy

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
    END DO
    DO i = 1,nconv
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
      rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
      the_here_c2(i)=the_here(index1(index2(i)))
      thp_here_c2(i)=thp_here(index1(index2(i)))
      qe_here_c2(i)=qe_here(index1(index2(i)))
      qp_here_c2(i)=qp_here(index1(index2(i)))
      rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
      dptot_c2(i)   = 0.0 ! dummy variable
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

! Force shallow convection to stop at ntpar unless using adaptive forced
! detrainment or shallow grey zone param

    IF (sh_on == 0 .AND. sh_grey == 0) THEN

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
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2_4a5a
    CALL convec2_4a5a(nconv,n_sh,nlev,ntra,k,sh_on, sh_sdet_on, sh_new_termc,&
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
                   (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)  &
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
          detrain_up(index1(index2(i)),k+1) =                         &
                 -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
        END IF
      END DO
    END IF

  END IF   ! nconv > 0

! Write CONVEC2 compressed output arrays back to full fields

  DO i = 1,n_sh
    thp(i,k+1)    = 0.0
    qp(i,k+1)     = 0.0
    xpk(i,k+1)    = 0.0
    flx(i,k+1)    = 0.0
    depth(i)      = 0.0
    precip(i,k+1) = 0.0
    qclp(i,k+1)   = 0.0
    qcfp(i,k+1)   = 0.0
  END DO
  DO i = 1,n_sh
    bgmk(i)       = .FALSE.
    bterm(i)      = .FALSE.
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO i = 1,n_sh
        trap(i,k+1,ktra) = 0.0
      END DO
    END DO
  END IF

  IF (l_mom_gk) THEN
    DO i = 1,n_sh
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
    IF (l_mom) THEN    ! needed for all versions
      DO i = 1,nconv
        flxkp12(index1(index2(i)),k)  = flxkp12_c2(i)
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
  END IF      ! nconv > 0

!   Write adaptive diagnostics for this level to full array for output

  DO i = 1,n_sh
    rbuoy_p_out(i,k) = rbuoy_p_here(i)
    the_out(i,k) = the_here(i)
    thp_out(i,k) = thp_here(i)
    qe_out(i,k) = qe_here(i)
    qp_out(i,k) = qp_here(i)
  END DO

!  write rbuoy_here to rbuoy_p_old for next pass through loop
  DO i = 1,n_sh
    rbuoy_p_old(i) = rbuoy_p_here(i)
  END DO

!-----------------------------------------------------------------------
! 3.4  Store information at termination for use later
!-----------------------------------------------------------------------

  DO i = 1,n_sh
    IF (bterm(i)) THEN
      cca_2d_term(i)  = cca_2d(i)
      ! Set final detrainment level (but not used)
      det_lev(i)= k+1
      flx_init_term(i) = flx_init(i)
      bgmk_term(i) = bgmk(i)
      ! If convection has terminated write cape to diagnostic output
      ! variable (cape_out).
      cape_out(i) = cape(i)
      dcpbydt(i) = 0.0
      cape(i) = 0.0
      ! Set kterm array which holds the level index for termination
      ! of convection.
      kterm(i) = k
      bconv(i) = .FALSE.
    END IF
  END DO

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

END DO

!-----------------------------------------------------------------------
! 4.0 Work out scaled mass flux
! First scale by grey zone parametrization, if requested
! Then to keep cfl ratio below limit.
! Note L_CAPE not applied to shallow convection
!-----------------------------------------------------------------------

DO i = 1,n_sh
  IF (kterm(i)  >=  ntml(i)) THEN   ! actual shallow convection

    flx_init_new(i) = flx_init_term(i)
    IF (sh_grey == 1) THEN
      !---------------------------------------------------------------
      ! Calculate Honnert et al style subgrid weighting as a function
      ! of the cloud-top height to grid size ratio
      !---------------------------------------------------------------
      k = kterm(i)
      weight_param = 1.0 -                                            &
             TANH( beta_cu*z_theta(i,k)/delta_smag(i)) *              &
             MAX( 0.0, (4.0-delta_smag(i)/z_theta(i,k)) )/4.0
      ! Scale flx_init using parametrization weighting
      flx_init_new(i) = weight_param*flx_init_new(i)
    END IF

    max_cfl(i) = max_cfl(i) * timestep

    IF (max_cfl(i)  >   cfl_limit) THEN
      flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
    END IF

    IF (flx_init_new(i)  >   flxmax_init(i)) THEN
      flx_init_new(i) = flxmax_init(i)
    END IF
    max_cfl(i) = 0.0

    IF (flx_init_new(i) > 0.0) THEN
      scale_f(i) = flx_init_new(i) / flx_init_term(i)
      ! set flx_init to the new value to provide the real initial mass
      ! flux in all conditions
      flx_init(i) = flx_init_new(i)   ! needed for later
    END IF
  END IF
END DO

kmax_term = 2

DO i = 1,n_sh
  IF(kterm(i)+1 >  kmax_term) THEN
    kmax_term = kterm(i)+1
  END IF
END DO

IF (kmax_term > nlev) THEN
  kmax_term = nlev
END IF

IF (l_safe_conv) THEN

  ! Ensure scale_f not a very small number otherwise may get an almost zero
  ! mass flux which may cause problems in CMT or downdraught calculation 

  DO i = 1,n_sh
    ! No proper shallow convection
    IF (scale_f(i) < 1.e-6 .OR. kterm(i) <= ntml(i)) THEN 
!      write(6,*) ' Problem shallow : i ',cape_out(i),scale_f(i),flx_init(i),&
!                   thpert(i),qpert(i),thv_pert(i),ntml(i),ntpar(i),kterm(i)
!      write(6,*) ' dtheta : ',(dthbydt(i,k),k=1,ntpar(i)+2)
!      write(6,*) ' dq     : ',(dqbydt(i,k),k=1,ntpar(i)+2)
!      write(6,*) ' dqcl   : ',(dqclbydt(i,k),k=1,ntpar(i)+2)
!      write(6,*) ' dqcf   : ',(dqcfbydt(i,k),k=1,ntpar(i)+2)
!      write(6,*) ' buoy   : ',(rbuoy_p_out(i,k),k=1,ntpar(i)+2)
!      write(6,*) ' precip   : ',(precip(i,k),k=1,ntpar(i)+2)
      scale_f(i)  = 0.0
      flx_init(i) = 0.0
      flx_init_new(i) = 0.0
      cape_out(i) = 0.0         ! reset to zero as failed deep
      ind_shall(i) = 0.0  ! real shallow event indicator
    ELSE
      ind_shall(i) = 1.0  ! real shallow event indicator
    END IF
  END DO    

ELSE 
  ! Original code no check that shallow events are real
  DO i = 1,n_sh
    ind_shall(i) = 1.0  ! All shallow as 1.0
  END DO    
END IF  ! l_safe_conv

!-----------------------------------------------------------------------
! 4.1  Carry out closure scaling and cfl scaling on updraught values
!-----------------------------------------------------------------------
IF (l_safe_conv) THEN
  ! Note if scale_f(i) = 0.0 then all increments are zeroed - no shallow
  ! convection occurs.

  DO kt = 2, kmax_term
    DO i = 1,n_sh
  
      !   Note any failed deep convection has increments scaled by zero.
      IF (kt  >=  ntml(i) .AND. kt <= kterm(i)+1 ) THEN

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
            IF (kt <  kterm(i)+1) THEN
              flxkp12(i,kt) = flxkp12(i,kt) * scale_f(i)
            END IF
          END IF
          IF (l_tracer) THEN
            DO ktra = 1,ntra
              dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            END DO
          END IF
          flx(i,kt)    = flx(i,kt) * scale_f(i)
          precip(i,kt) = precip(i,kt) * scale_f(i)

          IF (flg_up_flx .OR. flg_mf_shall) THEN
            up_flux(i,kt) = flx(i,kt)
          END IF
          IF (flg_entr_up) THEN
            entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)
          END IF
          IF (flg_detr_up) THEN
            detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
          END IF

      END IF !kt>ntml and flx_init_new>0
    END DO  ! i loop
  END DO  ! kt loop

  ! Indicate no shallow convection occurred to other processes by resetting 
  ! kterm after setting existing increments to zero through the above scaling.
  DO i = 1,n_sh
    IF (scale_f(i) < 1.e-6) THEN
      kterm(i) = 0 
    END IF
  END DO    


ELSE    ! original unsafe code where failing shallow increments remain 

  DO kt = 2, kmax_term
    DO i = 1,n_sh
  
      IF (kt  >=  ntml(i) .AND. kt <= kterm(i)+1 .AND.                       &
                                       flx_init_new(i) >   0.0 ) THEN
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
            IF (kt <  kterm(i)+1) THEN
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
          IF (flg_up_flx .OR. flg_mf_shall) THEN        
            up_flux(i,kt) = flx(i,kt)
          END IF
          IF (flg_entr_up) THEN
            entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)
          END IF
          IF (flg_detr_up) THEN
            detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
          END IF
 
      END IF !kt>ntml and flx_init_new>0
    END DO  ! i loop
  END DO  ! kt loop

END IF

! Mass flux on half levels 
IF (flg_up_flx_half) THEN
  DO kt = 1, nlev-1
    DO i = 1,n_sh
      up_flux_half(i,kt+1) = flxkp12(i,kt)
    END DO    
  END DO    
END IF

! Scale cloud fraction

IF (l_ccrad) THEN

  DO i = 1,n_sh
    cca_2d(i) = cca_2d_term(i)
    IF (flx_init_new(i) > 0.0) THEN
      cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))
    END IF

    ! Check scaled cloud fraction not smaller than minimum value
    ! (2.0E-5) or greater than unity.
    !
    ! (Was moved out of scaling if test to ensure these limits
    ! at all times, not just when cca_2d is scaled)

    cca_2d(i) = MAX(2.0e-5, cca_2d(i))
    cca_2d(i) = MIN(1.0e+0, cca_2d(i))
  END DO      ! i

ELSE      ! original

  DO i = 1,n_sh
    cca_2d(i) = cca_2d_term(i)
    IF (flx_init_new(i) > 0.0) THEN
      cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

      ! Check scaled cloud fraction not smaller than minimum value
      ! (2.0E-5) or greater than unity.

      cca_2d(i) = MAX(2.0e-5, cca_2d(i))
      IF (cca_2d(i) > 1.0) THEN
        cca_2d(i) = 1.0
      END IF
    END IF
  END DO      ! i

END IF      ! l_ccrad


!-----------------------------------------------------------------------
! 5.0 All shallow convection will terminate at some level. This level
!     has been stored in the main level loop.
!     The convection will either have a down draught or none will be
!     possible.
!-----------------------------------------------------------------------

DO i = 1,n_sh
  IF (kterm(i) >= ntml(i)) THEN
    k = kterm(i)
    tempnum=0.0    ! stores total column precipitation from convection      
    IF (iccb(i) >  0) THEN
      deltap_cld=p_layer_centres(i,iccb(i)) -p_layer_centres(i,k)
      DO kt=iccb(i),k+1
        tempnum = tempnum + precip(i,kt)
      END DO
    ELSE
      deltap_cld = 0.0
    END IF

    ! Set logicals for use later
    ! Downdraughts allowed in the cloud depth > 15km
    ! and the total column precipitation is > 1e-12.
    IF (deltap_cld >  15000.0 .AND. bgmk_term(i) .AND.tempnum >  1e-12) THEN
      b_dd(i) = .TRUE.
    ELSE  
    ! If no downdraught then go through evaporation below cloud base only
      b_nodd(i) = .TRUE.
    END IF
  END IF
END DO  ! n_sh

IF (l_safe_conv) THEN
  ! Overrule b_nodd if shallow convection not real for the column
  ! so that neither the downdraught nor the evaporation below cloud base is
  ! called for failed shallow points.
   
  DO i = 1,n_sh    
    IF (b_nodd(i) .AND. scale_f(i) < 1.0e-6) THEN
      b_nodd(i) = .FALSE.
    END IF
  END DO
  ! Find out which points have real shallow convection
  n_real_sh=0
  DO i = 1,n_sh
    IF (scale_f(i) >= 1.0e-6) THEN
      n_real_sh=n_real_sh+1
      index1(n_real_sh) = i
    END IF
  END DO
ELSE
  n_real_sh = n_sh
  DO i = 1,n_sh
    index1(i) = i
  END DO
END IF

!-----------------------------------------------------------------------
! 5.1  Downdraft calculation - on all points where convection is
!      terminating.
!
!      Subroutine DD_ALL_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

npossdd = 0
DO i = 1,n_sh
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
  CALL dd_all_call_4a5a(n_sh,npossdd,kmax_term,nlev,trlev,ntra    &
                  , kterm,iccb,icct,index_possdd,l_tracer         &
                  , bwater(1,2),exner_layer_centres               &
                  , exner_layer_boundaries,p_layer_centres        &
                  , p_layer_boundaries,pstar,recip_pstar,timestep &
                  , cca_2d,thp,qp,th,q,qse,trap,tracer,flx,precip &
                  , dthbydt,dqbydt,dtrabydt,rain,snow,rain_3d     &
                  , snow_3d,dwn_flux,entrain_dwn,detrain_dwn)

END IF

!-----------------------------------------------------------------------
! 5.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
nnodd = 0
DO i = 1,n_sh

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
  CALL evap_bcb_nodd_all(n_sh,nnodd,kmax_term,kterm               &
,                      iccb, index_nodd, bwater(1,2)              &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      timestep , cca_2d, th, q, qse, precip      &
,                      dthbydt, dqbydt                            &
,                      rain, snow, rain_3d, snow_3d)
END IF

!-----------------------------------------------------------------------
! 6.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

IF (l_mom) THEN
  nterm = 0

  DO i = 1, n_sh
    ! Only do CMT for points where shallow convection has occurred
    ! Will not do for cases of failed shallow
    IF (ind_shall(i) == 1.0) THEN
      nterm = nterm + 1
      cu_term(nterm) = i
      cu_tend(nterm) = i
    END IF
  END DO

  IF (nterm  >   0) THEN

! DEPENDS ON: shallow_grad_stress
    CALL shallow_grad_stress(n_sh,n_sh,nterm,nlev,cu_term,        &
                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,   &
                             plcl_uv,ptop_uv,p_uv,phalf_uv,       &
                             rho_uv,ue_p,ve_p,timestep,           &
                             ! IN
                             uw,vw)

! DEPENDS ON: shallow_base_stress
    CALL shallow_base_stress(n_sh,n_sh,n_sh,nlev,nterm,cu_term,   &
                             cu_tend,nlcl_uv,ntop_uv,mb,wsc,      &
                             zlcl_uv,zcld,uw0,vw0,plcl_uv,        &
                             ptop_uv,ue_p,ve_p,phalf_uv,p_uv,     &
                             rho_uv,timestep,                     &
                             ! INOUT
                             uw,vw,                               &
                             ! OUT
                             uw_shall,vw_shall)

! DEPENDS ON: shallow_cmt_incr
    CALL shallow_cmt_incr(n_sh,n_sh,n_sh,nlev,nterm,cu_term,      &
                          cu_tend,nlcl_uv,ntop_uv,uw,vw,phalf_uv, &
                          rho_uv,zlcl_uv,                         &
                          !OUT
                          dubydt,dvbydt)

  END IF  ! nterm>0
END IF ! L_mom


!-----------------------------------------------------------------------
! 7.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.

DO i = 1,n_sh
  qminincolumn(i) = q(i,nlev)
END DO
DO k = 1,nlev-1
  DO i = 1,n_sh
    IF (q(i,k)  <   qminincolumn(i)) THEN
      qminincolumn(i) = q(i,k)
    END IF
  END DO
END DO

! Ensure Q does not go below global allowed minimum (QMIN)

DO i = 1,n_sh
  qminincolumn(i)=MAX(qmin,qminincolumn(i))
END DO

! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.

DO k = nlev,2,-1
  DO i = 1,n_sh
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
  DO i = 1,n_sh
    temp1(i)=q(i,k) + dqbydt(i,k) * timestep
    IF (temp1(i)  <   qminincolumn(i) .AND.                            &
        printstatus >= prstatus_normal ) THEN
      WRITE(6,'(a21,i6,a9,g26.18,a7,g26.18)') ' negative q shall, i:',  &
            i,' q after ',temp1(i),' dq/dt ',dqbydt(i,k)
    END IF
  END DO ! n_sh loop
!-----------------------------------------------------------------------
! 8.0  Mixing of the convective increments in the boundary
!      layer.
!-----------------------------------------------------------------------

IF (bl_cnv_mix == 1) THEN

!      Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.

! DEPENDS ON: mix_ipert_4a5a
  CALL mix_ipert_4a5a(n_sh, n_real_sh, nlev, nbl, ntml, index1,        &
                      p_layer_boundaries,                              &
                      exner_layer_centres, dthbydt, dqbydt, flx_init,  &
                      thpert, qpert)

ELSE

!      Mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.

! DEPENDS ON: mix_inc
  CALL mix_inc (n_sh,n_sh,n_sh,nlev,nbl,ntml,                       &
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
  CALL cor_engy_5a(n_sh,n_real_sh,nlev,index1, r_theta, r_rho                &
                    ,r2rho_th, r2rho, dr_across_th, dr_across_rh             &
                    ,exner_layer_centres, th, u, v                           &
                    ,dubydt, dvbydt, dqclbydt, dqcfbydt                      &
                    ,rain,snow,dqbydt,dthbydt)
END IF

!-----------------------------------------------------------------------
! 10.0  Calculate convective cloud amount on model levels - no anvils
!-----------------------------------------------------------------------
! Initialise output array

DO k = 1,nlev
  DO i = 1,n_sh
    cca(i,k) = 0.0
  END DO
END DO

!-----------------------------------------------------------------------
! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!-----------------------------------------------------------------------

IF (.NOT. l_ccrad) THEN
  DO i=1, n_sh
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

!-------------------------------------------------------------------------
! 10.1 CCRad - Calculate CCA fow shallow levels only
!-------------------------------------------------------------------------

IF (l_ccrad) THEN

  DO i=1, n_sh

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

      overlap_fac(i) = 2.0 * (  z_rho(i,ntpar(i)+1)                 &
                              - z_rho(i, ntml(i)+1) )               &
                     / z_rho(i,ntml(i)+1)

    END IF     ! iccb
  END DO     ! n_sh

  !---------------------------------------------------------------
  ! 10.11 Calculate CCA
  !---------------------------------------------------------------
  SELECT CASE (cca2d_sh_opt)
    CASE(grant_lock)

      DO i=1, n_sh
        IF (iccb(i) /= 0) THEN ! Shallow convection occured
          tempnum = 2.0*mb(i)/wsc(i)

          cca_2d(i) = MAX(2.0e-5, tempnum)

          ! Will be used by NAME, grab lowest cca_2d before any
          ! Tuning knobs applied
          lcca(i) = cca_2d(i)
        END IF     ! iccb
      END DO     ! n_sh


    CASE(total_condensed_water)
        ! cca_2d is left unchanged from that calculated in the
        ! code, which is based on TCW (Total Condensed Water)
        ! (TCW is a rate)

  END SELECT


  DO i=1, n_sh

    overlap_fac(i) = MAX( 0.5, overlap_fac(i) )
    overlap_fac(i) = MIN( 5.0, overlap_fac(i) )

    IF (overlap_fac(i)*cca_2d(i) > 0.99) THEN
      overlap_fac(i) = 0.99/cca_2d(i)
    END IF
  END DO      ! i (n_sh)


  !-------------------------------------------------------------------
  ! 10.12 Fill cca with cca_2d where non-zero ccw
  !-------------------------------------------------------------------
  DO k=1, nlev
    DO i=1, n_sh
      IF (iccb(i) /= 0) THEN ! Shallow convection occured

        IF (ccw(i,k) > 0.0) THEN

          zpr = (z_rho(i,k)          - z_rho(i,ntml(i)+1))            &
              / (z_rho(i,ntpar(i)+1) - z_rho(i,ntml(i)+1))

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
    END DO       ! i (n_sh)
  END DO       ! k (nlev)

ELSE        ! Non CCRAD option

  DO k=1, n_cca_lev
    DO i=1, n_sh
      IF (k >= iccb(i) .AND. k < icct(i)) THEN
        cca(i,k) = cca_2d(i)
      END IF
    END DO
  END DO

END IF      ! l_ccrad

!-----------------------------------------------------------------------
! 11.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SHALLOW_CONV_5A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE shallow_conv_5a
