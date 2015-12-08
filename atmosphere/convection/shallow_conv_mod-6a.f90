! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Shallow convection scheme

MODULE shallow_conv_6a_mod

IMPLICIT NONE

!
! Description:
!   Shallow convection scheme
!   Works only on points diagnosed as shallow in subroutine CONV_DIAG
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

SUBROUTINE shallow_conv_6a(nbl,nlev,ntra,n_cca_lev,n_sh,trlev,    &
                       bland,delthvu,exner_layer_centres,         &
                       exner_layer_boundaries, l_calc_dxek,       &
                       l_q_interact, l_tracer, ntml, ntpar,       &
                       pstar,p_layer_centres,                     &
                       p_layer_boundaries,                        &
                       z_theta, z_rho,                            &
                       r_theta, r_rho,                            &
                       rho_theta, rho,                            &
                       r2rho_th, r2rho,                           & 
                       dr_across_th, dr_across_rh,                &
                       q,th,timestep,u,v,w,uw0,vw0,               &
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
                       kterm,                                     &
                       ind_shallow)
                       
USE atmos_constants_mod, ONLY:                                    &
    cp, r, kappa, pref, repsilon, c_virtual


USE water_constants_mod, ONLY: lc, lf, tm

USE cv_run_mod, ONLY:                                             &
    l_mom, sh_pert_opt, bl_cnv_mix, icvdiag,                      &
    cca2d_sh_opt, cca_sh_knob, ccw_sh_knob, limit_pert_opt,       &
    cnv_wat_load_opt, l_cv_conserve_check


USE cv_param_mod, ONLY:                                           &
    total_condensed_water, grant_lock,                            &
    thpixs_shallow, qpixs_shallow, c_mass,                        &
    max_sh_thpert, min_sh_thpert, max_sh_qpert_fac, beta_cu

USE cv_dependent_switch_mod, ONLY:                                &
    sh_on, mdet_sh_on, sh_ent_on, sh_new_termc, sh_grey

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
USE lift_par_6a_mod
USE convec2_6a_mod
USE water_loading_mod
USE cor_engy_6a_mod
USE mix_ipert_6a_mod

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

REAL, INTENT(IN) :: z_theta(n_sh,nlev)      ! height of theta levels (m)
REAL, INTENT(IN) :: z_rho(n_sh,nlev)        ! height of rho levels (m)
REAL, INTENT(IN) :: r_theta(n_sh,0:nlev)    ! radius of theta levels (m)
REAL, INTENT(IN) :: r_rho(n_sh,nlev)        ! radius of rho levels (m)
REAL, INTENT(IN) :: rho_theta(n_sh,nlev)    ! density for theta lev (kg/m3)
REAL, INTENT(IN) :: rho(n_sh,nlev)          ! density for rho lev (kg/m3)
REAL, INTENT(IN) :: r2rho_th(n_sh,nlev)     ! radius**2 density for 
                                            ! theta lev (kg/m)
REAL, INTENT(IN) :: r2rho(n_sh,nlev)        ! radius**2 density for 
                                            ! rho lev (kg/m)
REAL, INTENT(IN) :: dr_across_th(n_sh,nlev) ! thickness of theta levels (m)
REAL, INTENT(IN) :: dr_across_rh(n_sh,nlev) ! thickness of rho levels (m)

REAL, INTENT(IN)    :: q(n_sh,nlev) ! Model mixing ratio (kg/kg)

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

REAL, INTENT(OUT) :: up_flux(n_sh,nlev) ! Updraught mass flux (Pa/s)

REAL, INTENT(OUT) :: up_flux_half(n_sh,nlev) ! Updraught mass flux on 
                                             !half levels(Pa/s)
                                             
REAL, INTENT(OUT) :: dwn_flux(n_sh,nlev) ! Downdraught mass flux (Pa/s)

REAL, INTENT(OUT) :: uw_shall(n_sh,nlev) ! X-comp. of stress
                                         ! from shallow convection (kg/m/s2)

REAL, INTENT(OUT) :: vw_shall(n_sh,nlev) ! Y-comp. of stress
                                         ! from shallow convection (kg/m/s2)

REAL, INTENT(OUT) :: tcw(n_sh)  ! Total condensed water(kg/m2/s)

REAL, INTENT(OUT) :: cca_2d(n_sh) ! 2D convective cloud amount (%)

REAL, INTENT(OUT) :: ind_shallow(n_sh)  ! 1.0 if real shallow convection 
                                        ! else 0.0

INTEGER, INTENT(OUT) :: kterm(n_sh) ! termination level for shallow
                                    ! convection
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

REAL :: zk(n_sh)               !heights for use in calc

REAL :: zkp12(n_sh)            !of moist static energy

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

REAL :: ekp14(n_sh)             ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekp34(n_sh)             ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekm14(n_sh)             ! Entrainment coefficients at
                                ! level k-1+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: exk(n_sh)               ! Exner ratio at layer k

REAL :: exkp1(n_sh)             ! Exner ratio at layer k+1

REAL :: flxmax(n_sh)            ! Maximum initial convective
                                ! mass flux (Pa/s)

REAL :: flx_init(n_sh)          ! Initial mass flux at cloud base
                                ! (Pa/s)

REAL :: flx_init_new(n_sh)      ! flx_init scaled (Pa/s)

REAL :: flxmax_init(n_sh)       ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer
                                ! in Pa/s)

REAL :: max_cfl(n_sh)           ! Max cfl ratio over a convecting
                                ! layer

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


INTEGER :: start_lev_c2(n_sh)    ! Compressed convection
                                ! initiation level

REAL :: wsc(n_sh)               ! Convective velocity scale (m/s)

REAL :: wsc_o_mb(n_sh)          ! Convective velocity scale /mb

LOGICAL :: bgmk(n_sh)           ! Mask for points where parcel in
                                ! layer k is saturated

LOGICAL :: blatent(n_sh)        ! Mask for points where latent heat has 
                                ! been released

LOGICAL :: bwater(n_sh,2:nlev)  ! Mask for points at which
                                ! condensate is liquid

LOGICAL :: bwk(n_sh)            !mask for liquid condensate on k

LOGICAL :: bwkp1(n_sh)          !mask for liquid condensate on k+1

LOGICAL :: blowst(n_sh)         ! Dummy variable indicating low
                                ! enough stability for convection
                                ! to occur

LOGICAL :: bterm(n_sh)          ! Mask for points which have
                                ! stopped convecting

LOGICAL :: bconv(n_sh)          ! Mask for points at which
                                ! convection is occurring

LOGICAL :: bcposs(n_sh)         ! Mask for points passing
                                ! initial stability test


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

REAL :: flx(n_sh,nlev)          ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(n_sh,nlev)     ! Minmum parcel buoyancy excess

REAL :: thpixs_v(n_sh,nlev)     ! Theta parcel excess (K)

REAL :: qpixs_v(n_sh,nlev)      ! Q parcel excess(kg/kg)

REAL :: qclp(n_sh,nlev)         ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

REAL :: qcfp(n_sh,nlev)         ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)


! Parameters

REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed
REAL, PARAMETER :: minflx = TINY(flx_init_new)  ! minimum allowable
                                                ! initial mass flux

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

INTEGER :: index_nterm(n_sh)    ! Index for points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage

REAL :: scale_f(n_sh)           ! store scaling factor

REAL :: weight_param            ! Weighting factor

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

INTEGER :: ntpar_max            ! max ntpar value


! parameters ect for qmin checks

REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

REAL :: qminincolumn(n_sh)     ! Minimum value for q in column (kg/kg)

REAL :: temp1(n_sh)            ! work array

! Local compressed arrays

LOGICAL :: bgmkp1_c(n_sh), bgmkp1_c2(n_sh) ! Mask for points
                                ! where parcel in layer k+1
                                ! is saturated

LOGICAL :: bwk_c(n_sh), bwk_c2(n_sh) ! bwater mask in layer k

LOGICAL :: bwkp1_c(n_sh), bwkp1_c2(n_sh) ! bwater mask in layer k+1

LOGICAL :: blatent_c2(n_sh)     ! Mask for points where latent heat has 
                                ! been released

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

REAL :: flxkp1_c2(n_sh)         ! Parcel mass flux in layer k+1
                                ! (Pa/s)

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

REAL :: qsekp1_c2(n_sh)         ! Saturation mixing ratio
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

REAL :: traek_c2(n_sh,ntra)  ! Env. Tracer content in layer k (kg/kg)
REAL :: traekp1_c2(n_sh,ntra)! Env. Tracer content in layer k+1 (kg/kg)
REAL :: trapk_c2(n_sh,ntra)  ! Parcel Tracer content in layer k (kg/kg)
REAL :: trapkp1_c2(n_sh,ntra)! Parcel Tracer content in layer k+1 (kg/kg)

REAL :: rbuoyk_c(n_sh), rbuoyk_c2(n_sh)       ! Par. buoyancy at k (K)
REAL :: rbuoykp1_c(n_sh),rbuoykp1_c2(n_sh)    ! Par. buoyancy at k+1 (K)

REAL :: watldek_c(n_sh), watldek_c2(n_sh)     ! Env. water loading
                                              ! in layer k (kg/kg)
REAL :: watldpk_c(n_sh), watldpk_c2(n_sh)     ! Par. water loading
                                              ! in layer k (kg/kg)
REAL :: watldekp1_c(n_sh), watldekp1_c2(n_sh) ! Env. water loading
                                              ! in layer k+1 (kg/kg)
REAL :: watldpkp1_c(n_sh), watldpkp1_c2(n_sh) ! Par. water loading
                                              ! in layer k+1 (kg/kg)

REAL :: Qlkp1_c(n_sh),   &   ! Amount of condensation to liquid water 
        Qlkp1_c2(n_sh)       ! in the parcel (kg/kg)
REAL :: Qfkp1_c(n_sh),   &   ! Amount of deposition to ice water
        Qfkp1_c2(n_sh)       ! in the parcel (kg/kg)
REAL :: Frezkp1_c(n_sh), &   ! Amount of freezing from liquid
        Frezkp1_c2(n_sh)     ! to ice water in the parcel (kg/kg)

REAL :: uek_c2(n_sh)    ! Env. U in layer k (m/s)
REAL :: uekp1_c2(n_sh)  ! Env. U in layer k+1 (m/s)
REAL :: vek_c2(n_sh)    ! Env. V in layer k (m/s)
REAL :: vekp1_c2(n_sh)  ! Env. V in layer k+1 (m/s)
REAL :: upk_c2(n_sh)    ! Parcel U in layer k (m/s)
REAL :: upkp1_c2(n_sh)  ! Parcel U in layer k+1 (m/s)
REAL :: vpk_c2(n_sh)    ! Parcel V in layer k (m/s)
REAL :: vpkp1_c2(n_sh)  ! Parcel V in layer k+1 (m/s)

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
REAL :: ccwkp1_c2(n_sh)
REAL :: dcpbydt_c2(n_sh)
REAL :: delpk_c2(n_sh)
REAL :: delpkp12_c2(n_sh)
REAL :: delpkp1_c2(n_sh)
REAL :: delp_uv_k_c2(n_sh)
REAL :: delp_uv_kp1_c2(n_sh)
REAL :: depth_c2(n_sh)
REAL :: dptot_c2(n_sh)
REAL :: eflux_u_ud_c2(n_sh)
REAL :: eflux_v_ud_c2(n_sh)
REAL :: ekp14_c(n_sh),ekp14_c2(n_sh)
REAL :: ekp34_c(n_sh),ekp34_c2(n_sh)
REAL :: exk_c(n_sh), exk_c2(n_sh)
REAL :: exkp1_c(n_sh),exkp1_c2(n_sh)
REAL :: expi_c2(n_sh)
INTEGER :: icct_c2(n_sh)
INTEGER :: iccb_c2(n_sh)
INTEGER :: lctop_c2(n_sh)
INTEGER :: lcbase_c2(n_sh)
REAL :: lcca_c2(n_sh)
REAL :: max_cfl_c2(n_sh)
REAL :: pk_c(n_sh),pk_c2(n_sh)
REAL :: pkp1_c(n_sh),pkp1_c2(n_sh)
REAL :: pstar_c2(n_sh)
REAL :: qpi_c2(n_sh)
REAL :: relh_c2(n_sh)
REAL :: rbuoy_p_here_c2(n_sh)
REAL :: the_here_c2(n_sh)
REAL :: thp_here_c2(n_sh)
REAL :: qe_here_c2(n_sh)
REAL :: qp_here_c2(n_sh)
REAL :: tcw_c2(n_sh)
REAL :: thpi_c2(n_sh)
REAL :: xsbmin_v_c2(n_sh)
REAL :: qsat_lcl(n_sh)         ! not used
LOGICAL :: b_nodd(n_sh)   ! points with no downdraught
LOGICAL :: b_dd(n_sh)     ! points with downdraught on termination

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






IF (lhook) CALL dr_hook('SHALLOW_CONV_6A',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

DO i = 1,n_sh
  blowst(i)    = .TRUE.
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  b_nodd(i)    = .FALSE.
  b_dd(i)      = .FALSE.
  blatent(i)   = .FALSE.
END DO

DO i =1, n_sh
  kterm(i)        = 0
  ind_shallow(i)  = 0.0
END DO

l_mom_gk = .FALSE.       ! not Gregory-Kershaw CMT
                         ! Shallow code uses turbulence based CMT

!-----------------------------------------------------------------------
! 2.1  Initialise parcel properties and increment arrays
!-----------------------------------------------------------------------

!intialise parcel values over all levels
DO k = 1, nlev
  DO i = 1, n_sh
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
    DO i = 1,n_sh
      up(i,k) = 0.0
      vp(i,k) = 0.0
    END DO
  END DO
END IF

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k=1,nlev
      DO i = 1,n_sh
        trap(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF

DO k = 1,nlev
  DO i = 1,n_sh
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
    DO i = 1,n_sh
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,n_sh
        dtrabydt(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
IF (flg_up_flx .OR. flg_mf_shall) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      up_flux(i,k)      = 0.0
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
      dwn_flux(i,k)     = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      entrain_up(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      detrain_up(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
     entrain_dwn(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      detrain_dwn(i,k)  = 0.0
    END DO
  END DO
END IF
IF (l_mom) THEN
  IF (flg_uw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        uw_shall(i,k)   = 0.0
      END DO
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        vw_shall(i,k)   = 0.0
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

DO k = 1,n_cca_lev
  DO i = 1,n_sh
    cca(i,k) = 0.0
  END DO
END DO

DO i = 1,n_sh
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling calculations
!-----------------------------------------------------------------------
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)         = 0.0
  cape_out(i)     = 0.0
  dcpbydt(i)      = 0.0
  max_cfl(i)      = 0.0
  det_lev(i)      = 0

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


!initialise ekm14 for unused adaptive entrainment
DO i =1, n_sh
  ekm14(i) =0.0
END DO

!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
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

! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! Reset thpixs and qpixs at ntml
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
      thv_pert(i) =  (th(i,k+1) * (1.0 + c_virtual                &
                    * q(i,k+1)) - th(i,k)                         &
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

  DO i = 1,n_sh

    k = ntml(i)

    IF (t_lcl(i) >  tm) THEN
      dq_sat_env = repsilon * lc * qse_lcl(i)                     &
                      / (r * t_lcl(i) * t_lcl(i))

!         Estimate of moist adiabatic lapse rate

      dthv_ma    = ( (lc/cp) - (1.+c_virtual)*th(i,k) )*          &
                   dq_sat_env*(g/cp)/(1.+(lc/cp)*dq_sat_env)
    ELSE
      dq_sat_env = repsilon * (lc+lf) * qse_lcl(i)                &
                      / (r * t_lcl(i) * t_lcl(i))

!         Estimate of moist adiabatic lapse rate (in K/m)

      dthv_ma    = ( ((lc+lf)/cp) - (1.+c_virtual)*th(i,k) )*     &
                   dq_sat_env*(g/cp)/(1.+((lc+lf)/cp)*dq_sat_env)
    END IF

    b_calc   = t_lcl(i) * c_virtual * dq_sat_env + 1.0            &
                        + c_virtual * qse_lcl(i)


!       Calculate theta_v perturbation:

!       First convert moist adiabatic lapse rate to thv difference
!       between levels k and k+1

    rho_k   = p_layer_centres(i,k) /                              &
            (r * th(i,k) * (p_layer_centres(i,k)/ pref)**kappa)
    dthv_ma = -dthv_ma*                     &
            (p_layer_centres(i,k+1)-p_layer_centres(i,k)) /       &
            (rho_k*g)

!       Make perturbation relative to a target lapse rate (namely
!       0.6*dthv_ma, which is approximately what is seen in LES)

    IF (Kprof_cu == off) THEN
      thv_pert(i) = -0.17 * wthvs(i) / mb(i) +  0.6*dthv_ma       &
                - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))         &
                - th(i,k)  *(1.0 + c_virtual*q(i,k))  )
    ELSE
      ! "entrainment flux" at LCL done by BL scheme
      thv_pert(i) = 0.6*dthv_ma                                   &
                    - (  th(i,k+1)*(1.0 + c_virtual*q(i,k+1))     &
                    - th(i,k)  *(1.0 + c_virtual*q(i,k))  )

    END IF

!       limit thv_pert to physically sensible values

    thv_pert(i) = MAX(MIN(thv_pert(i), max_sh_thpert),            &
                  min_sh_thpert)

    c_calc      = th_lcl(i) * c_virtual * (qse_lcl(i) - q_lcl(i)) &
                - thv_pert(i)

    thpert(i)   = MAX(MIN(-c_calc / b_calc, max_sh_thpert),       &
                  min_sh_thpert)  ! ignore term in thpert**2

    thpixs_v(i,k) = thpert(i)

    qpert(i)    = MAX(MIN(qse_lcl(i) + ((p_lcl(i) / pref)         &
                **kappa) * thpert(i) * dq_sat_env                 &
                - q_lcl(i),                                       &
                  max_sh_qpert_fac * qse_lcl(i)),0.0)

    qpixs_v(i,ntml(i))  = qpert(i)

  END DO !n_sh

END IF

! Set bwater=.true. on points where water will condense rather than
! ice.
! DEPENDS ON: flag_wet
CALL flag_wet(n_sh,n_sh,nlev,th,exner_layer_centres,bwater)


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

IF (sh_on == 1 .OR. sh_grey == 1) THEN
                     ! Adaptive forced detrainment or grey shallow param
                     ! No limit on convection top
  ntpar_max = nlev-3 ! What is a sensible value to have here?

ELSE                 ! Top limited

  ntpar_max=0
  DO i = 1,n_sh
    IF (ntpar(i)+1 >  ntpar_max) THEN
      ntpar_max=ntpar(i)+1
    END IF
  END DO
! Ensure that ntpar_max does not exceed nlev-1
  ntpar_max = MIN(ntpar_max, nlev-1)
END IF

DO k = 2,ntpar_max  !loop over model levels

!-----------------------------------------------------------------------
! Initialise environment variables
! NB These variable are only used by layer_cn.
!-----------------------------------------------------------------------
  DO i = 1,n_sh
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
! Initialise parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k and has not convected in column before
!----------------------------------------------------------------------- 
  DO i = 1,n_sh
    IF ( .NOT. bconv(i) .AND. det_lev(i) == 0) THEN
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


!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
  CALL layer_cn(k,n_sh,nlev,                                      &
                   mdet_sh_on, sh_ent_on,                         &
                   ntml,ntpar,                                    &
                   .TRUE.,.FALSE.,.FALSE.,                        &
                   bconv,bwk,bwkp1,                               &
                   exner_layer_boundaries,                        &
                   exner_layer_centres,                           &
                   p_layer_boundaries,p_layer_centres,            &
                   recip_pstar,entrain_coef,rhum,                 &
                   zk, zkp12, zkp1,                               &
                   thek, qek,qsek, thekp1,qekp1,qsekp1,           &
                   thpk,qpk , qsat_lcl, ekm14,                    &
                   pkp1,delpkp1,exkp1,                            &
                   pk,delpk,delpkp12,exk,delexkp1,                &
                   delp_uv_k, delp_uv_kp1,                        &
                   ekp14,ekp34,amdetk)

! Set ekm14 for next pass through loop
  DO i = 1, n_sh
    ekm14(i) = ekp14(i)
  END DO

! Maximum initial convective mass flux
  DO i = 1,n_sh
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  END DO

! Scale entrainment coefficients with cloud base mass flux
! and convective velocity scale
  DO i = 1,n_sh
    IF (k  >=  ntml(i)) THEN
      ! If original entrainment coefficient then scale as before else leave
      ! as calculated by layer_cn.
      IF (entrain_coef(i) < 0.0 ) THEN
        ekp14(i)  = ekp14(i)  * wsc_o_mb(i)
        ekp34(i)  = ekp34(i)  * wsc_o_mb(i)
        amdetk(i) = amdetk(i) * wsc_o_mb(i)
      END IF
    END IF
  END DO  ! n_sh
  
!-----------------------------------------------------------------------
! Initial test to check if convection is possible in layer k
!-----------------------------------------------------------------------
! Convection is possible if
! - the point was convecting (bconv = .T.) and did not terminate
!   in the previous layer
! - or if at the top level of the surface mixed layer (k = ntml)
  DO i = 1,n_sh
    bcposs(i) = bconv(i) .OR. k  ==  ntml(i)
  END DO  ! n_sh

! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
  ncposs = 0
  DO i = 1,n_sh
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

! Allow parcel to convect from ntml.
      IF (k  ==  ntml(index1(i))) THEN
        bconv(index1(i))  = .TRUE.  ! convection active
        blowst(index1(i)) = .TRUE.  ! convection initialised in layer

! Set parcel mass flux
        flxk_c(i)         = mb(index1(i)) * g                       &
                          * p_layer_centres(index1(i),k)            &
                          / ( r * thpk_c(i)                         &
                          * (p_layer_centres(index1(i),k)/ pref)**kappa )

! Write compressed mass flux back to full array
        flx(index1(i),k)  = flxk_c(i)

! Store diagnostics linked to initial convective mass flux for
! calculation of final closure.
        flx_init(index1(i))    = flxk_c(i)
        flxmax_init(index1(i)) = flxmax(index1(i))

! Set mixing detrainment rate at first convecting level to zero
        amdetk(index1(i)) = 0.0

      ELSE
        blowst(index1(i)) = .FALSE. ! convection not initialise in layer
      END IF

! Reset threshold for forced detrainment to the initial
! (potentially negative) buoyancy
      xsbmin_v(index1(i),k) = thv_pert(index1(i))

    END DO  !ncposs
  END IF    !ncposs>0
  
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
      start_lev_c2(i)   = ntml(index1(index2(i)))
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
      blowst_c2(i)      = blowst(index1(index2(i)))
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
      relh_c2(i)        = 0.0 ! dummy variable
      dptot_c2(i)       = 0.0 ! dummy variable
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

! Force shallow convection to stop at the parcel top from conv_diag (ntpar)
! unless using adaptive forced detrainment or grez-zone shallow.
    IF (sh_on == 0 .AND. sh_grey == 0) THEN   
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

    CALL convec2_6a  (k, nconv, n_sh, nlev, ntra, sh_on, sh_new_termc,       &
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
  DO i = 1,n_sh
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
      !dummy                          = relh_c2(i)
      !dummy                          = dptot_c2(i)
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

  END IF      ! nconv > 0

!-----------------------------------------------------------------------
! 3.4  CFL scaling
!-----------------------------------------------------------------------
! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_sh) with index_nterm

  nterm = 0
  DO i = 1,n_sh
    IF (bterm(i)) THEN
      nterm = nterm + 1
      index_nterm(nterm) = i
    END IF
  END DO

  IF (nterm >  0) THEN

! Work out scaled mass flux needed to keep cfl ratio below limit.
! First scale by grey zone parametrization, if requested
! Then to keep cfl ratio below limit.
! Note L_CAPE not applied to shallow convection

    DO j = 1,nterm
      i = index_nterm(j)

      flx_init_new(i) = flx_init(i)
      IF (sh_grey == 1) THEN
        !-----------------------------------------------------------
        ! Calculate Honnert et al style subgrid weighting as a 
        ! function of the cloud-top height to grid size ratio
        !-----------------------------------------------------------
        weight_param = 1.0 -                                         &
              TANH( beta_cu*z_theta(i,k)/delta_smag(i)) *            &
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
    END DO      ! j (nterm)

    !
    ! Scale cloud fraction
    !
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
        kterm(i)        = 0
      ELSE
        ! True convection
        kterm(i)        = k
        ind_shallow(i)  = 1.0
      END IF
    END DO  ! nterm

!-----------------------------------------------------------------------
! Apply cfl scaling
!-----------------------------------------------------------------------
    DO kt = 2, k+1
      DO j = 1,nterm
        i = index_nterm(j)
        IF (kt  >=  ntml(i)) THEN
          dthbydt(i,kt)   = dthbydt(i,kt)  * scale_f(i)
          dqbydt(i,kt)    = dqbydt(i,kt)   * scale_f(i)
          dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
          dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
          dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
          dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
          dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
          IF (l_mom) THEN
            dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
            dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
          END IF
          IF (l_tracer) THEN
            DO ktra = 1,ntra
              dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
            END DO
          END IF

          flx(i,kt)    = flx(i,kt)    *  scale_f(i)
          precip(i,kt) = precip(i,kt) *  scale_f(i)

        END IF !kt >ntml and flx_init_new >0
      END DO  ! j loop
    END DO  ! kt loop


    DO j = 1,nterm
      i = index_nterm(j)
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
        deltap_cld=p_layer_centres(i2,iccb(i2)) -p_layer_centres(i2,k)
        DO kt=iccb(i2),k+1
          tempnum=tempnum+precip(i2,kt)
        END DO
      ELSE
        deltap_cld = 0.0
      END IF

! Set logical to determine if downdraughts are allowed or not.
      IF (deltap_cld >  15000.0.AND.bgmk(i2) .AND.tempnum >  1e-12) THEN
        b_dd(i2) = .TRUE.
      ELSE
        b_nodd(i2) = .TRUE.
      END IF
    END DO  ! nterm loop

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

    DO j = 1,nterm
      i=index_nterm(j)
      cape_out(i) = cape(i)
      dcpbydt(i) = 0.0
      cape(i) = 0.0
      bconv(i) = .FALSE.
      det_lev(i)= k+1 ! Set final detrainment level (but not used).
    END DO

  END IF  ! nterm > 0

!-----------------------------------------------------------------------
! Write out entrainment, detrainment and half-level mass flux diagnostics.
! They will be scaled by the full level mass flux outside 
! of the level loop
!-----------------------------------------------------------------------
! Calculate fractional entrainment rate for level k.
  IF (flg_entr_up) THEN
    DO i = 1,nconv
      entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i))        &
                * (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i) &
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
IF (flg_up_flx .OR. flg_mf_shall) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      up_flux(i,k) = flx(i,k)
    END DO
  END DO
END IF

IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_sh
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF


!-----------------------------------------------------------------------
! 4.0  Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer 
!-----------------------------------------------------------------------

CALL mix_ipert_6a(n_sh, nlev, nbl, ntml, p_layer_boundaries,       &
                  exner_layer_centres, dthbydt, dqbydt, flx_init,  &
                  thpert, qpert)

!-----------------------------------------------------------------------
! 5.0 All shallow convection will terminate at some level. This level
!     has been stored in the main level loop.
!     The convection will either have a down draught or none will be
!     possible.
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

! DEPENDS ON: dd_all_call_6a
  CALL dd_all_call_6a(n_sh,npossdd,kmax_term,nlev,trlev,ntra      &
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
    nterm = nterm + 1
    cu_term(nterm) = i
    cu_tend(nterm) = i
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
! 7.0  Energy correction calculation - removed as old code not correct
!     for new dynmaics grid ( attempts to correct this give problems).
!     UM documentation paper 27 - section 12.
!-----------------------------------------------------------------------
DO i = 1,n_sh
  index1(i) = i
END DO

  IF (l_cv_conserve_check) THEN
      CALL cor_engy_6a(n_sh,n_sh,nlev,index1,r2rho_th,                       &
                       dr_across_th,exner_layer_centres,p_layer_boundaries,  &
                       dqbydt,dqclbydt,dqcfbydt,rain,snow,                   &
                       dthbydt)
  END IF
!-----------------------------------------------------------------------
! 8.0  Correct negative/very small humidities
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
! 9.0  Calculate convective cloud amount on model levels - no anvils
!-----------------------------------------------------------------------
! Initialise output array

DO k = 1,nlev
  DO i = 1,n_sh
    cca(i,k) = 0.0
  END DO
END DO

!-------------------------------------------------------------------------
! 9.1 CCRad - Calculate CCA fow shallow levels only
!-------------------------------------------------------------------------

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
! 9.11 Calculate CCA
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
! 9.12 Fill cca with cca_2d where non-zero ccw
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

!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SHALLOW_CONV_6A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE shallow_conv_6a
END MODULE shallow_conv_6a_mod
