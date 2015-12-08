! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE glue_conv_4a_mod

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  IMPLICIT NONE

CONTAINS
!
!  Gather-scatter routine for deep and shallow convection points
!
SUBROUTINE glue_conv_4a(np_field,npnts,nlev,nbl,call_number,seg_num  &
,                 th,q,qcl,qcf, qrain, qgraup, qcf2               &
,                 cf_liquid,cf_frozen,bulk_cf,pstar               &
,                 bland,u,v,w                                     &
,                 tracer, dthbydt, dqbydt, dqclbydt, dqcfbydt     &
,                 dcflbydt, dcffbydt, dbcfbydt, dubydt, dvbydt    &
,                 rain, snow, rain_3d, snow_3d                    &
,              cca0, iccb0, icct0, cclwp0, ccw0, lcbase0, cca0_2d &
, lctop, lcca, cca,  iccb,  icct,  cclwp,  ccw,  lcbase,  cca_2d  &
,                 freeze_lev, deep_cfl_limited, mid_cfl_limited   &
,                 l_mid_all,kterm_deep,kterm_shall                &
,                 precip_deep,precip_shall,precip_mid,precip_cong &
,                 wstar_dn,wstar_up,mb1,mb2,kterm_congest         &
,                 n_cumulus,uw0,vw0,w_max,zlcl,zlcl_uv,ztop_uv    &
,                 entrain_coef,deep_flag,past_precip,past_conv_ht &
,                 cape_out                                        &
,                 n_dp, n_cg, n_sh, n_md                          &
,                 r_rho,r_theta,rho,rho_theta,delta_smag          &
,                 exner_layer_boundaries                          &
,                 exner_layer_centres                             &
,                 p_layer_boundaries                              &
,                 p_layer_centres                                 &
,                 z_theta, z_rho                                  &
,                 timestep,t1_sd,q1_sd                            &
,                 ntml,ntpar,conv_type,l_shallow_bl               &
,                 l_pc2_diag_sh_pts, l_congestus, l_mid           &
,                 cumulus_bl,wstar,wthvs,delthvu_bl,ql_ad         &
,                 qsat_lcl,ftl,fqt                                &
,                 l_tracer,ntra,trlev,n_cca_lev                   &
,                 l_mixing_ratio, l_mcr_qrain, l_mcr_qgraup       & 
,                 l_mcr_qcf2, l_calc_dxek,l_q_interact            &
,                 up_flux_half, up_flux, dwn_flux                 &
,                 entrain_up,detrain_up,entrain_dwn,detrain_dwn   &
,                 uw_deep,vw_deep,uw_shall,vw_shall,uw_mid,vw_mid &
,         wqt_flux_sh,wthetal_flux_sh,wthetav_flux_sh,wql_flux_sh &
,          mf_deep,mf_congest,mf_shall,mf_midlev                  &
,          dt_deep,dt_congest,dt_shall,dt_midlev                  &
,          dq_deep,dq_congest,dq_shall,dq_midlev                  &
,          du_deep,du_congest,du_shall,du_midlev                  &
,          dv_deep,dv_congest,dv_shall,dv_midlev                  &
,          ind_cape_reduced, cape_ts_used, ind_deep, ind_shall    &
                 )

! Purpose:
!   Gather-scatter routine for deep and shallow convection points.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v6 programming standards

USE cv_run_mod,  ONLY:                                            &
    l_mom,        adapt,        termconv,                         &
    cca2d_sh_opt, cca2d_md_opt, cca2d_dp_opt,                     &
    cca_sh_knob,  cca_md_knob,  cca_dp_knob,                      &
    ccw_sh_knob,  ccw_md_knob,  ccw_dp_knob,                      &
    bl_cnv_mix,   l_anvil,      l_dcpl_cld4pc2,                   &
    l_ccrad,      l_3d_cca,     l_pc2_diag_sh

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                     &
    flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,         &
    flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall,             &
    flg_uw_mid, flg_vw_mid

USE earth_constants_mod, ONLY: g

USE water_constants_mod, ONLY: tm

USE ereport_mod, ONLY : ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE


!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:


INTEGER, INTENT(IN) :: np_field ! No of points in a full field
                                !(note all multi-dimensional fields
                                !  passed to routine MUST be
                                !   dimensioned with np_field)

INTEGER, INTENT(IN) :: npnts    ! No. of points in segment

INTEGER, INTENT(IN) :: nlev     ! No. of model layers

INTEGER, INTENT(IN) :: ntml(npnts) ! Top level of surface mixed
                                   ! layer defined relative to
                                   ! theta,q grid

INTEGER, INTENT(IN) :: ntpar(npnts) ! Top level of initial parcel
                                    ! ascent in BL scheme defined
                                    ! relative to theta,q grid

INTEGER, INTENT(IN) :: n_cca_lev! No. of convective cloud
                                ! amount levels (1 for 2D,
                                ! nlevs for 3D)

INTEGER, INTENT(IN) :: nbl      ! No. of boundary layer levels

INTEGER, INTENT(IN) :: call_number ! Current convection sweep

INTEGER, INTENT(IN) :: seg_num  ! segment number 

INTEGER, INTENT(IN) :: ntra     ! No. of tracer fields

INTEGER, INTENT(IN) :: trlev    ! No. of model levels on which
                                ! tracers are included

INTEGER, INTENT(IN) :: n_dp     ! Number of deep points

INTEGER, INTENT(IN) :: n_cg     ! Number of congestus points

INTEGER, INTENT(IN) :: n_sh     ! Number of shallow points

INTEGER, INTENT(IN) :: n_md     ! Number of mid points

REAL, INTENT(IN) :: delthvu_bl(npnts) !Integral of undilute parcel
                                      ! buoyancy over convective cloud
                                      ! layer (Kelvin m)

REAL, INTENT(IN) :: &
  ql_ad(npnts)      & ! adiabatic liquid water content at inversion (kg/kg)
, qsat_lcl(npnts)   & ! qsat at cloud base (kg/kg) (NOT USED 4A)
, ftl(npnts)        & ! Surface sensible heat flux from BL
                      ! (W/m2) i.e. cp*rho*w'tl'
, fqt(npnts)          ! Total water flux from surface (kg/m2/s)
                      ! i.e. rho*w'qT'

REAL, INTENT(IN) :: exner_layer_centres(np_field,0:nlev) !Exner

REAL, INTENT(IN) :: exner_layer_boundaries(np_field,0:nlev)!Exner
                                                           ! at half level above
                                                           ! exner_layer_centres

REAL, INTENT(IN) :: pstar(npnts) ! Surface pressure (Pa)

REAL, INTENT(IN) :: p_layer_centres(np_field,0:nlev) !Pressure(Pa)


REAL, INTENT(IN) :: p_layer_boundaries(np_field,0:nlev) ! Pressure
                                                        ! at half level above
                                                        ! p_layer_centres (Pa)

REAL, INTENT(IN) :: z_theta(np_field,nlev) ! height of theta
                                           ! levels above surface (m)

REAL, INTENT(IN) :: z_rho(np_field,nlev)    ! height of rho levels
                                            ! above surface (m)


REAL, INTENT(IN) :: t1_sd(np_field) ! Standard deviation of
                                    ! turbulent fluctuations of
                                    ! layer 1 temp. (K)

REAL, INTENT(IN) :: q1_sd(np_field) ! Standard deviation of
                                    ! turbulent fluctuations of
                                    ! layer 1 q (kg/kg)

REAL, INTENT(IN) :: th(np_field,nlev) ! Model potential temperature (K)

REAL, INTENT(IN) :: q(np_field,nlev) ! Model water vapour (kg/kg)

REAL, INTENT(IN) :: timestep    ! Model timestep (s)

REAL, INTENT(IN) :: uw0(np_field) ! U-comp of surface stress(N/m2)

REAL, INTENT(IN) :: vw0(np_field) ! V-comp of surface stress(N/m2)

REAL, INTENT(IN) :: w_max(np_field) ! max w in column

REAL, INTENT(IN) :: u(np_field,nlev) !Model U field (m/s)

REAL, INTENT(IN) :: v(np_field,nlev) !Model V field (m/s)

REAL, INTENT(IN) :: w(np_field,nlev) !Model W field (m/s)

REAL, INTENT(IN) :: wstar(npnts) ! Convective velocity scale (m/s)

REAL, INTENT(IN) :: wthvs(npnts) ! Surface flux of THV  (Pa m/s2)

REAL, INTENT(IN) :: zlcl(npnts) ! Lifting condensation level accurate
                                ! height (m) (NOT USED)

REAL, INTENT(IN) :: zlcl_uv(npnts) ! Lifting condensation level
                                   ! defined for the uv grid (m)

REAL, INTENT(IN) :: ztop_uv(npnts) ! Top of cloud layer
                                   ! defined for the uv grid (m)

REAL, INTENT(IN) :: r_rho(np_field,nlev)      ! radius rho lev (m)
REAL, INTENT(IN) :: r_theta(np_field,0:nlev)  ! radius th lev (m)
REAL, INTENT(IN) :: rho(np_field,nlev)        ! density kg/m3
REAL, INTENT(IN) :: rho_theta(np_field,nlev)  ! th lev kg/m3
REAL, INTENT(IN) :: delta_smag(np_field)      ! grid size used in smag (m)

LOGICAL, INTENT(IN) :: l_tracer ! Switch for inclusion of tracers

LOGICAL, INTENT(IN) :: &
  l_mixing_ratio       & ! true - input moisture mixing ratio
                         ! false - input moisture specific humidity
 ,l_mcr_qrain          & ! true - prognostic rain
 ,l_mcr_qgraup         & ! true - prognostic graupel
 ,l_mcr_qcf2           & ! true - prognostic qcf2 
 ,l_calc_dxek          & ! Switch for calculation of condensate increment
 ,l_q_interact           ! Switch allows overwriting parcel variables when
                         ! calculating condensate incr.

LOGICAL, INTENT(IN) :: l_shallow_bl(npnts) ! Shallow cumulus indicator

LOGICAL, INTENT(IN) :: l_pc2_diag_sh_pts(npnts) ! Carry
                                                ! diagnostic shallow convective
                                                ! information for PC2

LOGICAL, INTENT(IN) :: l_congestus(npnts) ! congestus cumulus

LOGICAL, INTENT(IN) :: l_mid(npnts)       ! possible mid-level cnv

LOGICAL, INTENT(IN) :: cumulus_bl(npnts)  ! Cumulus indicator

LOGICAL, INTENT(IN) :: bland(npnts) ! Land/sea mask

! NOT USED 4A
REAL, INTENT(IN) ::   &
 entrain_coef(npnts)  &  ! entrainmentceofficients
,deep_flag(npnts)     &  ! history of deep convection
,past_precip(npnts)   &  ! history of convective precip
,past_conv_ht(npnts)     ! history of convective height


! Arguments with intent INOUT:

! NOTE - All moist variables passed down to this routine are :
! specific humidities if l_mixing_ratio = .false.
! mixing ratios       if l_mixing_ratio = .true.

REAL, INTENT(INOUT) ::      &
  qcl(np_field,nlev)        & ! Liq condensate (kg/kg)
 ,qcf(np_field,nlev)        & ! Ice condensate (kg/kg)
 ,qrain(np_field,nlev)      & ! rain (kg/kg)
 ,qgraup(np_field,nlev)     & ! graupel (kg/kg)
 ,qcf2(np_field,nlev)       & ! 2nd ice type (kg/kg)
 ,cf_liquid(np_field,nlev)  & ! Liq water cloud volume (fraction)
 ,cf_frozen(np_field,nlev)  & ! Frozen water cloud volume (fraction?)
 ,bulk_cf(np_field,nlev)      ! Bulk total cloud volume ( )


REAL, INTENT(INOUT) :: tracer(np_field,trlev,ntra) !Model tracer
                                                   ! fields (kg/kg)


! Arguments with intent OUT:


REAL, INTENT(OUT) :: dqclbydt(np_field,nlev) ! Increments to liq
                                ! condensate due to convection
                                ! (kg/kg/s)

REAL, INTENT(OUT) :: dqcfbydt(np_field,nlev) ! Increments to ice
                                ! condensate due to convection
                                ! (kg/kg/s)

REAL, INTENT(OUT) :: dcflbydt(np_field,nlev) ! Increments to liq
                                ! cloud volume due to convection
                                ! (/s)

REAL, INTENT(OUT) :: dcffbydt(np_field,nlev) ! Increments to ice
                                ! cloud volume due to convection
                                ! (/s)

REAL, INTENT(OUT) :: dbcfbydt(np_field,nlev) ! Increments to
                                ! total cld volume due to
                                ! convection(/s)

REAL, INTENT(OUT) :: dthbydt(np_field,nlev) ! Increments to
                         ! potential temp. due to convection (K/s)

REAL, INTENT(OUT) :: dqbydt(np_field,nlev) ! Increments to q due
                                           ! to convection (kg/kg/s)

REAL, INTENT(OUT) :: dubydt(np_field,nlev+1) ! Increments to U due
                                             ! to CMT (m/s2)

REAL, INTENT(OUT) :: dvbydt(np_field,nlev+1) ! Increments to V due
                                             ! to CMT (m/s2)

REAL, INTENT(OUT) :: rain(npnts) ! Surface convective rainfall (kg/m2/s)

REAL, INTENT(OUT) :: snow(npnts) ! Surface convective snowfall (kg/m2/s)

REAL, INTENT(OUT) :: rain_3d(np_field,nlev) ! convective rainfall flux (kg/m2/s)

REAL, INTENT(OUT) :: snow_3d(np_field,nlev) ! convective snowfall flux (kg/m2/s)


! Section 5 Convective Cloud properties
REAL, INTENT(OUT) ::       &
  cca(np_field,n_cca_lev)  &! Cnv. cld amount (0-1)
, ccw(np_field,nlev)       &! Cnv. in-cld water (kg/kg)
, cclwp(npnts)             &! Cond. water path (kg/m^2)
, lcca(npnts)               ! Cnv. cld amount at base of lowest
                            ! cloud layer (no anvil). (0-1)

INTEGER, INTENT(OUT) ::    &
  iccb(npnts)              &! Cnv. cld base level (Highest Cld layer)
, icct(npnts)              &! Cnv. cld top  level (Highest Cld layer)
, lcbase(npnts)            &! Cnv. cld base level (Lowest  Cld layer)
, lctop(npnts)              ! Cnv. cld top  level (Lowest  Cld layer)


! Section 0 Convective Cloud properties for Radiative impacts
REAL, INTENT(OUT) ::       &
  cca0(np_field,n_cca_lev) &! Cnv. cld amount (0-1)
, ccw0(np_field,nlev)      &! Cnv. in-cld water (kg/kg)
, cclwp0(npnts)             ! Cond. water path (kg/m^2)

INTEGER, INTENT(OUT) ::    &
  iccb0(npnts)             &! Cnv. cld base level (Highest Cld layer)
, icct0(npnts)             &! Cnv. cld top  level (Highest Cld layer)
, lcbase0(npnts)            ! Cnv. cld base level (Lowest  Cld layer)

INTEGER, INTENT(OUT) :: freeze_lev(npnts) !index for freezing lev

INTEGER, INTENT(OUT) :: kterm_deep(npnts) ! index deep conv
INTEGER, INTENT(OUT) :: kterm_shall(npnts) ! termination level for shallow

LOGICAL, INTENT(OUT) :: l_mid_all(npnts)  ! on exit true if mid level
                                          ! convection has triggered

REAL, INTENT(OUT) ::      &
  deep_cfl_limited(npnts) & ! indicator for cfl limited deep conv
 ,mid_cfl_limited(npnts)    ! indicator for cfl limited mid-levelconv

REAL, INTENT(OUT) :: precip_deep(npnts)  ! deep precip (kg/m2/s)

REAL, INTENT(OUT) :: precip_shall(npnts) ! shallow precip(kg/m2/s)

REAL, INTENT(OUT) :: precip_mid(npnts)   ! mid precip (kg/m2/s)
REAL, INTENT(OUT) :: precip_cong(npnts)  !congestus precip(kg/m2/s)

REAL, INTENT(OUT) :: up_flux(np_field,nlev) ! Updraught mass flux (Pa/s)

REAL, INTENT(OUT) :: up_flux_half(np_field,nlev)
                                ! Updraught mass flux
                                ! on half levels (Pa/s)

REAL, INTENT(OUT) :: dwn_flux(np_field,nlev) ! Downdraught mass
                                ! flux (Pa/s)

REAL, INTENT(OUT) :: entrain_up(np_field,nlev) ! Fractional
                                ! entrainment rate into updraughts
                                ! (Pa/s)

REAL, INTENT(OUT) :: detrain_up(np_field,nlev) ! Fractional
                                ! detrainment rate into updraughts
                                ! (Pa/s)

REAL, INTENT(OUT) :: entrain_dwn(np_field,nlev) ! Fractional
                                ! entrainment rate into
                                ! downdraughts (Pa/s)

REAL, INTENT(OUT) :: detrain_dwn(np_field,nlev) ! Fractional
                                ! detrainment rate into
                                ! downdraughts (Pa/s)

! Diagnostics relating to momentum fluxes

REAL, INTENT(OUT) ::     &
 uw_deep(np_field,nlev)  & ! X-comp. of stress from deep convection (kg/m/s2)
,vw_deep(np_field,nlev)  & ! Y-comp. of stress from deep convection (kg/m/s2)
,uw_shall(np_field,nlev) & ! X-comp. of stress from shallow convection (kg/m/s2)
,vw_shall(np_field,nlev) & ! Y-comp. of stress from shallow convection (kg/m/s2)
,uw_mid(np_field,nlev)   & ! U comp of stress from mid convection (kg/m/s2)
,vw_mid(np_field,nlev)     ! V comp of stress from mid convection (kg/m/s2)

REAL, INTENT(OUT) :: cape_out(npnts) ! Saved convective available
                                     ! potential energy for diagnostic
                                     ! output (Jkg-1)

! Flux diagnostics - only available from turbulence based schemes
! (not actually used in this routine - OUT shouldnt be used)

REAL              ::                                              &
  wqt_flux_sh(np_field,nlev)                                      &
                                 ! w'qt' flux
, wql_flux_sh(np_field,nlev)                                      &
                                ! w'ql' flux
, wthetal_flux_sh(np_field,nlev)                                  &
                                ! w'thetal' flux
, wthetav_flux_sh(np_field,nlev)                                  &
                                ! w'thetav' flux
, wstar_dn(npnts)                                                 &
, wstar_up(npnts)                                                 &
, mb1(npnts)                                                      &
, mb2(npnts)                                                      &
, mf_deep(np_field,nlev),mf_congest(np_field,nlev)                &
, mf_shall(np_field,nlev),mf_midlev(np_field,nlev)                &
, dt_deep(np_field,nlev),dt_congest(np_field,nlev)                &
, dt_shall(np_field,nlev),dt_midlev(np_field,nlev)                &
, dq_deep(np_field,nlev),dq_congest(np_field,nlev)                &
, dq_shall(np_field,nlev),dq_midlev(np_field,nlev)                &
, du_deep(np_field,nlev),du_congest(np_field,nlev)                &
, du_shall(np_field,nlev),du_midlev(np_field,nlev)                &
, dv_deep(np_field,nlev),dv_congest(np_field,nlev)                &
, dv_shall(np_field,nlev),dv_midlev(np_field,nlev)

REAL, INTENT(OUT) ::         &
 ind_cape_reduced(np_field)  & ! indicates cape timescale reduced
,cape_ts_used(np_field)      & ! cape timescale used for deep (s)
,ind_deep(np_field)          & ! 1.0 if deep convection 0.0 otherwise
,ind_shall(np_field)           ! 1.0 if shallow convection 0.0 otherwise

! This is not used in this routine (shouldnt use OUT)
INTEGER           ::                                         &
  kterm_congest(npnts)   ! termination level for congestus

!-----------------------------------------------------------------------
! Redundant arguments
!-----------------------------------------------------------------------

INTEGER, INTENT(IN) :: n_cumulus
INTEGER, INTENT(IN) :: conv_type(npnts)

!-----------------------------------------------------------------------
! local required in move mixing to glue

REAL :: dtrabydt(npnts,nlev,ntra) ! Increment to tracer due to
                                  ! convection (kg/kg/s)


CHARACTER (LEN=12), PARAMETER ::  routinename = 'glue_conv_4a'
CHARACTER (LEN=80) ::                                             &
  cmessage                   ! error message
INTEGER :: errorstatus       ! error status

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to DEEP convection scheme.
! Arrays are identified by underscore DP (_dp) and are of length
! n_dp where n_dp is the number of points diagnosed as deep in the
! boundary layer diagnosis routine. For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

INTEGER :: error_point      ! location of problem deep point

INTEGER :: dpi(n_dp)        ! index for deep points in full grid

INTEGER :: ntml_dp(n_dp)

INTEGER :: ntpar_dp(n_dp)

LOGICAL :: bland_dp(n_dp)
 
REAL ::                                                            &
  pstar_dp(n_dp)                                                   &
, recip_pstar_dp(n_dp)                                             &
, t1_sd_dp(n_dp)                                                   &
, q1_sd_dp(n_dp)                                                   &
, uw0_dp(n_dp)                                                     &
, vw0_dp(n_dp)                                                     &
, zlcl_uv_dp(n_dp)                                                 &
, wstar_dp(n_dp)                                                   &
, delthvu_dp(n_dp)

REAL ::                                                            &
  p_layer_centres_dp(n_dp,0:nlev)                                  &
, p_layer_boundaries_dp(n_dp,0:nlev)                               &
, exner_layer_centres_dp(n_dp,0:nlev)                              &
, exner_layer_boundaries_dp(n_dp,0:nlev)                           &
, r2rho_dp(n_dp,nlev)                                              &
, r2rho_th_dp(n_dp,nlev)                                           &
, rho_dp(n_dp,nlev)                                                &
, rho_theta_dp(n_dp,nlev)                                          &
, dr_across_th_dp(n_dp,nlev)                                       &
, dr_across_rh_dp(n_dp,nlev)                                       &
, z_theta_dp(n_dp,nlev)                                            &
, z_rho_dp(n_dp,nlev)                                              &
, r_theta_dp(n_dp,0:nlev)                                          &
, r_rho_dp(n_dp,nlev)

REAL :: w_max_dp(n_dp)

REAL :: u_dp(n_dp,nlev)

REAL :: v_dp(n_dp,nlev)

REAL :: th_dp(n_dp,nlev)

REAL :: q_dp(n_dp,nlev)

REAL :: tracer_dp(n_dp,trlev,ntra)

REAL :: dthbydt_dp(n_dp,nlev)

REAL :: dqbydt_dp(n_dp,nlev)

REAL :: dubydt_dp(n_dp,nlev+1)

REAL :: dvbydt_dp(n_dp,nlev+1)

REAL :: rain_dp(n_dp)

REAL :: snow_dp(n_dp)

REAL :: rain_3d_dp(n_dp,nlev)

REAL :: snow_3d_dp(n_dp,nlev)

REAL :: tcw_dp(n_dp)

INTEGER ::                                                        &
  iccb_dp(n_dp)                                                   &
, icct_dp(n_dp)                                                   &
, lcbase_dp(n_dp)                                                 &
, lctop_dp(n_dp)                                                  &
, freeze_lev_dp(n_dp)

REAL :: cclwp_dp(n_dp)

REAL :: cca_dp(n_dp,n_cca_lev)

REAL :: ccw_dp(n_dp,nlev)

REAL :: lcca_dp(n_dp)

REAL :: up_flux_dp(n_dp,nlev)

REAL :: up_flux_half_dp(n_dp,nlev)

REAL :: dwn_flux_dp(n_dp,nlev)

REAL :: entrain_up_dp(n_dp,nlev)

REAL :: detrain_up_dp(n_dp,nlev)

REAL :: entrain_dwn_dp(n_dp,nlev)

REAL :: detrain_dwn_dp(n_dp,nlev)

REAL :: uw_deep_dp(n_dp,nlev)

REAL :: vw_deep_dp(n_dp,nlev)

REAL :: cape_out_dp(n_dp)

REAL :: qse_dp(n_dp,nlev)

!PC2
REAL ::                                                           &
  qcl_dp(n_dp,nlev)                                               &
, qcf_dp(n_dp,nlev)                                               &
, cf_liquid_dp(n_dp,nlev)                                         &
, cf_frozen_dp(n_dp,nlev)                                         &
, bulk_cf_dp(n_dp,nlev)                                           &
, dqclbydt_dp(n_dp,nlev)                                          &
, dqcfbydt_dp(n_dp,nlev)                                          &
, dcflbydt_dp(n_dp,nlev)                                          &
, dcffbydt_dp(n_dp,nlev)                                          &
, dbcfbydt_dp(n_dp,nlev)

REAL :: dtrabydt_dp(n_dp,nlev,ntra) ! Increment to tracer due to
                                    ! convection (kg/kg/s)
INTEGER :: kterm_dp(n_dp)    ! required by mid level scheme

REAL :: cca_2d_dp(n_dp)      ! required by mid level scheme

REAL :: ind_cape_reduced_dp(n_dp)
REAL :: cape_ts_used_dp(n_dp)
REAL :: cfl_limited_dp(n_dp)
REAL :: ind_deep_dp(n_dp)

!  compressed SCM diagnostics for adaptive modset

REAL :: rbuoy_p_out_dp(n_dp,nlev)

REAL :: the_out_dp(n_dp,nlev)

REAL :: thp_out_dp(n_dp,nlev)

REAL :: qe_out_dp(n_dp,nlev)

REAL :: qp_out_dp(n_dp,nlev)
!temporary local diagnostics
INTEGER :: thp_out_copy_count
INTEGER :: thp_out_md_count
INTEGER :: th_count

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to SHALLOW convection scheme.
! Arrays are identified by underscore SH (_sh) and are of length
! n_sh where n_sh is the number of points diagnosed as shallow in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

INTEGER :: shi(n_sh)      ! index for shallow points in full grid

INTEGER :: ntml_sh(n_sh)

INTEGER :: ntpar_sh(n_sh)

REAL :: pstar_sh(n_sh)

REAL :: recip_pstar_sh(n_sh)

REAL :: p_layer_centres_sh(n_sh,0:nlev)

REAL :: p_layer_boundaries_sh(n_sh,0:nlev)

LOGICAL :: bland_sh(n_sh)

REAL :: exner_layer_centres_sh(n_sh,0:nlev)

REAL :: exner_layer_boundaries_sh(n_sh,0:nlev)

REAL ::                        &
  z_theta_sh(n_sh,nlev)        &
, z_rho_sh(n_sh,nlev)          &
, r2rho_th_sh(n_sh,nlev)       &
, dr_across_th_sh(n_sh,nlev)

REAL :: delthvu_sh(n_sh)

REAL :: t1_sd_sh(n_sh)

REAL :: q1_sd_sh(n_sh)

REAL :: uw0_sh(n_sh)

REAL :: vw0_sh(n_sh)

REAL :: wstar_sh(n_sh)

REAL :: wthvs_sh(n_sh)

REAL :: zlcl_uv_sh(n_sh)

REAL :: ztop_uv_sh(n_sh)

REAL :: u_sh(n_sh,nlev)

REAL :: v_sh(n_sh,nlev)

REAL :: th_sh(n_sh,nlev)

REAL :: q_sh(n_sh,nlev)

REAL :: tracer_sh(n_sh,trlev,ntra)

REAL :: dthbydt_sh(n_sh,nlev)

REAL :: dqbydt_sh(n_sh,nlev)

REAL :: dubydt_sh(n_sh,nlev+1)

REAL :: dvbydt_sh(n_sh,nlev+1)

REAL :: rain_sh(n_sh)

REAL :: snow_sh(n_sh)

REAL :: rain_3d_sh(n_sh,nlev)

REAL :: snow_3d_sh(n_sh,nlev)

REAL :: tcw_sh(n_sh)

INTEGER :: iccb_sh(n_sh)

INTEGER :: icct_sh(n_sh)

REAL :: cclwp_sh(n_sh)

REAL :: cca_sh(n_sh,n_cca_lev)

REAL :: ccw_sh(n_sh,nlev)

REAL :: lcca_sh(n_sh)

INTEGER :: lcbase_sh(n_sh)

INTEGER :: lctop_sh(n_sh)

INTEGER :: freeze_lev_sh(n_sh)

REAL :: up_flux_sh(n_sh,nlev)

REAL :: up_flux_half_sh(n_sh,nlev)

REAL :: dwn_flux_sh(n_sh,nlev)

REAL :: entrain_up_sh(n_sh,nlev)

REAL :: detrain_up_sh(n_sh,nlev)

REAL :: entrain_dwn_sh(n_sh,nlev)

REAL :: detrain_dwn_sh(n_sh,nlev)

REAL :: uw_shall_sh(n_sh,nlev)

REAL :: vw_shall_sh(n_sh,nlev)

REAL :: cape_out_sh(n_sh)

REAL :: qse_sh(n_sh,nlev)

REAL :: qcl_sh(n_sh,nlev)

REAL :: qcf_sh(n_sh,nlev)

REAL :: cf_liquid_sh(n_sh,nlev)

REAL :: cf_frozen_sh(n_sh,nlev)

REAL :: bulk_cf_sh(n_sh,nlev)

REAL :: dqclbydt_sh(n_sh,nlev)

REAL :: dqcfbydt_sh(n_sh,nlev)

REAL :: dcflbydt_sh(n_sh,nlev)

REAL :: dcffbydt_sh(n_sh,nlev)

REAL :: dbcfbydt_sh(n_sh,nlev)
REAL :: dtrabydt_sh(n_sh,nlev,ntra) ! Increment to tracer due to
                                    ! convection (kg/kg/s)
REAL :: cca_2d_sh(n_sh)        ! required by mid level scheme

!  compressed SCM diagnostics for adaptive modset

REAL :: rbuoy_p_out_sh(n_sh,nlev)

REAL :: the_out_sh(n_sh,nlev)

REAL :: thp_out_sh(n_sh,nlev)

REAL :: qe_out_sh(n_sh,nlev)

REAL :: qp_out_sh(n_sh,nlev)

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to MID-LEVEL convection scheme.
! Arrays are identified by underscore MD (_md) and are of length
! npnts where npnts is the total number of points in the grid (since
! mid-level convection may occur on any point previously diagnosed as
! shallow or deep).
!-----------------------------------------------------------------------

!  compressed SCM diagnostics for adaptive modset

REAL :: rbuoy_p_out_md(n_md,nlev)

REAL :: the_out_md(n_md,nlev)

REAL :: thp_out_md(n_md,nlev)

REAL :: qe_out_md(n_md,nlev)

REAL :: qp_out_md(n_md,nlev)

INTEGER :: midtrig(n_md)   ! Level at which mid level convection
                           ! may start
INTEGER :: mdi(n_md)       ! index for mid points in full grid

REAL :: p_layer_centres_md(n_md,0:nlev)

REAL :: p_layer_boundaries_md(n_md,0:nlev)

REAL :: exner_layer_centres_md(n_md,0:nlev)

REAL :: exner_layer_boundaries_md(n_md,0:nlev)

REAL :: rho_theta_md(n_md,nlev), rho_md(n_md,nlev)

REAL :: z_theta_md(n_md,nlev), z_rho_md(n_md,nlev)
REAL :: r_theta_md(n_md,0:nlev), r_rho_md(n_md,nlev)

REAL :: w_max_md(n_md)

REAL :: u_md(n_md,nlev)

REAL :: v_md(n_md,nlev)

REAL :: th_md(n_md,nlev)

REAL :: q_md(n_md,nlev)

REAL :: tracer_md(n_md,trlev,ntra)

REAL :: dthbydt_md(n_md,nlev)

REAL :: dqbydt_md(n_md,nlev)

REAL :: dubydt_md(n_md,nlev+1)

REAL :: dvbydt_md(n_md,nlev+1)

REAL :: rain_md(n_md)

REAL :: snow_md(n_md)

REAL :: rain_3d_md(n_md,nlev)

REAL :: snow_3d_md(n_md,nlev)

REAL :: tcw_md(n_md)

INTEGER :: iccb_md(n_md)

INTEGER :: icct_md(n_md)

REAL :: cclwp_md(n_md)

REAL :: cca_md(n_md,n_cca_lev)

REAL :: ccw_md(n_md,nlev)

REAL :: lcca_md(n_md)

INTEGER :: lcbase_md(n_md)

INTEGER :: lctop_md(n_md)

INTEGER :: freeze_lev_md(n_md)

REAL :: t1_sd_md(n_md)

REAL :: q1_sd_md(n_md)

REAL :: pstar_md(n_md)

REAL :: recip_pstar_md(n_md)

LOGICAL :: bland_md(n_md)

INTEGER :: ntpar_md(n_md)

INTEGER :: ntml_md(n_md)

REAL :: up_flux_md(n_md,nlev)

REAL :: up_flux_half_md(n_md,nlev)

REAL :: dwn_flux_md(n_md,nlev)

REAL :: entrain_up_md(n_md,nlev)

REAL :: detrain_up_md(n_md,nlev)

REAL :: entrain_dwn_md(n_md,nlev)

REAL :: detrain_dwn_md(n_md,nlev)

REAL :: cape_out_md(n_md)

REAL :: cfl_limited_md(n_md)

REAL :: qcl_md(n_md,nlev)

REAL :: qcf_md(n_md,nlev)

REAL :: cf_liquid_md(n_md,nlev)

REAL :: cf_frozen_md(n_md,nlev)

REAL :: bulk_cf_md(n_md,nlev)

REAL :: dqclbydt_md(n_md,nlev)

REAL :: dqcfbydt_md(n_md,nlev)

REAL :: dcflbydt_md(n_md,nlev)

REAL :: dcffbydt_md(n_md,nlev)

REAL :: dbcfbydt_md(n_md,nlev)
REAL :: dtrabydt_md(n_md,nlev,ntra) ! Increment to tracer due to
                                    ! convection (kg/kg/s)
REAL :: qse_md(n_md,nlev)      ! Saturation mixing ratio of
                               ! cloud environment (kg/kg)
REAL :: cca_2d_md(n_md)        ! required by mid level scheme

LOGICAL :: l_mid_md(n_md)     ! true if mid level convection
                              ! compressed version of l_mid_all


!-----------------------------------------------------------------------
! Mixing ratio / specific humidity conversions
!-----------------------------------------------------------------------
!  The full model holds q, qcl and qcf, all of which should be used in the
!  conversions and evaluation of qt. (High resolution versions can also
!  hold other moist variables.) The current mass flux scheme without
!  PC2 assume qt=q i.e. it ignores the presence of qcl and qcf.
!  This assumption is also true of the new shallow turbulence based scheme. 
! 
!  Let m indicate mixing ratio and q specific then
!
!  Now taking account of all moist variables
!      qt = qv + qcl + qcf + (qrain + qgraup + qcf2)
!      mt = mv + mcl + mcf + (mrain + mgraup + mcf2)
!
! Note there could be problems if qcf2 is present as in this case qcf holds
! the snow like ice and qcf2 the ice crystals so what is "seen" and incremented
! by convection (i.e. qcf) could be inappropriate.
!
!  Not PC2
!      dqcl = 0.0,      dqcf = 0.0
!      dmcl = 0.0,      dmcf = 0.0
!
!  Before scheme
!  mv(t) = qv(t)/(1-qt(t))          or qv(t) = mv(t)/(1+mt(t))
!
!  after scheme
!  mv(t+dt) = qv(t+dt)/(1-qt(t+dt)) or qv(t+dt) = mv(t+dt)/(1+mt(t+dt))
!
!  where
!  dmv= mv(t+dt) -  mv(t)           or dqv  =  qv(t+dt)- qv(t)
!
! How to convert increments
!
!  dmv = (dqv + qv*dqt - qt*dqv)/[(1-qt)(1-qt-dqt)]
!
!   or
!
!  dqv = (dmv + mt*dmv - mv*dmt)/[(1+mt)(1+mt+dmt)]
!
! if mt = mv  then dqv = dmv /[(1+mt)(1+mt+dmt)]
!
! Conversion of qsat
! ------------------
! rsat - mixing ratio saturation value
!
! qsat=rsat/(1+rsat)      and rsat=qsat(1-qsat)
!
!-----------------------------------------------------------------------
!  Arrays required for mixing ratio to specific conversions

REAL :: qt(npnts,nlev)    ! total water variable either mixing
                          ! ratio or specific humidity (kg/kg)
                          ! qt = q + qcl + qcf + (qrain+qgraup+qcf2)
                          ! ( ) extra moist variables not in most runs. 
REAL ::                &
  dqt                  &  ! total moisture increment per timestep
, dqtt                 &  ! total moisture increment per s
, denom                   ! 1./denominator

!-----------------------------------------------------------------------
!  Arrays required by 3d CCA calculation

REAL :: cca_2d(npnts)        ! 2d convective cloud (section 5)
REAL :: cca0_2d(npnts)       ! 2d convective cloud (section 0)

REAL :: tcw(npnts)           ! total column water on all points

! Saturation mixing ratio calulation and 1/pstar

REAL :: recip_pstar(npnts)      ! Reciprocal of pstar array

REAL :: qse(npnts,nlev)         ! Saturation specific humidity of
                                ! cloud environment (kg/kg)

REAL :: pt(npnts)               ! Temporary store for P in calc.
                                ! of sat. mixing ratio. (Pa)
REAL :: tt(npnts)               ! Temporary store for T in calc.
                                ! of saturation mixing ratio. (K)
REAL :: ttkm1(npnts)            ! Temporary store for T in layer
                                ! k-1 for use in freezing level
                                ! calc. for anvil. (K)


!  compressed SCM diagnostics for adaptive modset

REAL :: rbuoy_p_out(npnts,nlev)

REAL :: the_out(npnts,nlev)

REAL :: thp_out(npnts,nlev)

REAL :: qe_out(npnts,nlev)

REAL :: qp_out(npnts,nlev)

! copies so we can see if mid convection modifies the profile

REAL :: rbuoy_p_out_copy(npnts,nlev)

REAL :: the_out_copy(npnts,nlev)

REAL :: thp_out_copy(npnts,nlev)

REAL :: qe_out_copy(npnts,nlev)

REAL :: qp_out_copy(npnts,nlev)

!  arrays required by grid calculations

REAL ::                    &
 r2rho_th(npnts,nlev)      & ! radius**2 density theta lev (kg/m)
,r2rho(npnts,nlev)         & ! radius**2 density rho lev (kg/m)
,dr_across_th(npnts,nlev)  & ! thickness of theta levels (m)
,dr_across_rh(npnts,nlev)    ! thickness of rho levels (m)

!  arrays required by energy correction

INTEGER :: index1(npnts)
INTEGER :: nconv_all         ! Total number of points convecting

!   required by check on -ve q

REAL :: qminincolumn(npnts)     ! Minimum value for q in column
                                ! (kg/kg)
REAL :: temp1(npnts)            ! work array

! array required by tracers at end

REAL :: limited_step(npnts)     ! Reduced step size for tracer
                                ! mixing

REAL :: step_test(npnts)        ! Work array used in reducing step


! Parameters


REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

REAL, PARAMETER :: safety_margin = 1.0e-100 ! Small number used in
                                ! tracer step reduction




! Loop counters


INTEGER :: i,j,k,ktra,idp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



IF (lhook) CALL dr_hook('GLUE_CONV_4A',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise output variables
! Note a lot of the variables are initialised in ni_conv_ctl so don't
! need to be iniitalised here.
!-----------------------------------------------------------------------

DO i=1, npnts
  ind_cape_reduced(i) = 0.0
  cape_ts_used(i)     = 0.0
  ind_deep(i)         = 0.0
  ind_shall(i)        = 0.0
  kterm_shall(i)      = ntpar(i)   ! set to same as ntpar
END DO


DO k = 1,nlev
  DO i = 1,npnts
    dqbydt(i,k)              = 0.0
    dthbydt(i,k)             = 0.0
    dqclbydt(i,k)            = 0.0
    dqcfbydt(i,k)            = 0.0
    dcflbydt(i,k)            = 0.0
    dcffbydt(i,k)            = 0.0
    dbcfbydt(i,k)            = 0.0
    ccw(i,k)                 = 0.0
  END DO
END DO

! Initialisation of qt used for mixing/specific conversions
! Taking account of all possible moist variables even if not use in the 
! convection scheme

DO k=1,nlev
  DO i = 1,npnts
    qt(i,k) = q(i,k) + qcl(i,k) + qcf(i,k)
  END DO
END DO

IF (l_mcr_qrain) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qrain(i,k)
    END DO
  END DO
END IF
IF (l_mcr_qgraup) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qgraup(i,k)
    END DO
  END DO
END IF
IF (l_mcr_qcf2) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qcf2(i,k)
    END DO
  END DO
END IF

IF (l_mom) THEN
  DO k=1,nlev+1
    DO i = 1,npnts
      dubydt(i,k)              = 0.0
      dvbydt(i,k)              = 0.0
    END DO
  END DO
END IF
IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,npnts
        dtrabydt(i,k,ktra)  = 0.0   ! tracer increments
      END DO
    END DO
  END DO
END IF

! Initialise diagnostics if requested

IF (flg_up_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux(i,k)             = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux_half(i,k)             = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dwn_flux(i,k)            = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_up(i,k)          = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_up(i,k)          = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_dwn(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_dwn(i,k)         = 0.0
    END DO
  END DO
END IF

! Extra variable initialisation for safety

DO i = 1,npnts
  tcw(i)      = 0.0
END DO

IF (flg_uw_dp) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_deep(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_vw_dp) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_deep(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_uw_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_shall(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_vw_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_shall(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_uw_mid) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_mid(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_vw_mid) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_mid(i,k) = 0.0
    END DO
  END DO
END IF


! Required to get same convective cloud as old scheme
! Need to pass values from deep and shallow to mid level.
! Currently conversion of 2d to 3d makes use of total 2d from
! shallow/deep and mid in a column. If in future the conversion
! does not work on a column total this could be removed.

DO i=1,n_md
  tcw_md(i) = 0.0
  iccb_md(i) = 0
  icct_md(i) = 0
  cca_2d_md(i) = 0.0
END DO

!-----------------------------------------------------------------------
!  grid information
!-----------------------------------------------------------------------
! Calculate quantities involving density etc for future use

DO k=1,nlev
  DO i=1,npnts
    dr_across_rh(i,k) = r_theta(i,k) - r_theta(i,k-1)
    r2rho(i,k)        = r_rho(i,k)*r_rho(i,k)*rho(i,k)
  END DO
END DO

! rho_theta only set for nlev-1

k=1     ! bottom theta level thicker
  DO i=1,npnts
    dr_across_th(i,k) = r_rho(i,k+1) - r_theta(i,0)
    r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
  END DO

DO k=2,nlev-1
  DO i=1,npnts
    dr_across_th(i,k) = r_rho(i,k+1) - r_rho(i,k)
    r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
  END DO
END DO

k=nlev     ! top layer  (hope not used ?
           !             assume density as layer below)
  DO i=1,npnts
    dr_across_th(i,k) = r_theta(i,nlev) - r_rho(i,k)
    r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k-1)
  END DO

!-----------------------------------------------------------------------
! 1.0 Section to calculate fields used by all convection types.
! Note cheaper to do here than in individual routines.
! Create saturation mixing ratio arrays  & calculate freeze level
!-----------------------------------------------------------------------

! Calculate 1/pstar and initialize freeze_lev array.

DO i = 1,npnts
  recip_pstar(i)=1.0 / pstar(i)
  freeze_lev(i) = 1
END DO

! Loop over levels

DO k = 1,nlev

! Find freezing level

  IF (k  ==  1) THEN
    DO i = 1,npnts
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)

! Commented out as initialisation sets freeze_lev to 1.
! Code left incase altered in future.

!            If (tt(i)  <   TM) then
!              freeze_lev(i) = k
!            End If

    END DO
  ELSE
    DO i = 1,npnts
      ttkm1(i) = tt(i)
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)
      IF (tt(i)  <   tm .AND. ttkm1(i)  >=  tm) THEN
        IF (freeze_lev(i) == 1) THEN
          freeze_lev(i) = k
        END IF
      END IF
    END DO
  END IF


! Calculate saturation specific humidity as mass flux schemes all expect
! specific input

! DEPENDS ON: qsat_mix
  CALL qsat_mix(qse(1,k),tt,pt,npnts,.false.)

END DO  ! nlev

!-----------------------------------------------------------------------
! 1.0 DEEP Convection
! 1.1 Compress input variable arrays for deep convection scheme to
!     length n_dp (deep points only)
!-----------------------------------------------------------------------
IF (n_dp  >   0) THEN
  j = 0
  DO i = 1,npnts
    IF (cumulus_bl(i).AND..NOT.l_shallow_bl(i)) THEN
      j                        = j+1
      dpi(j)                   = i
    END IF
  END DO

  ! Initialise output variables
  DO k=1, nlev
    DO i=1, n_dp
      rain_3d_dp(i,k) = 0.0
      snow_3d_dp(i,k) = 0.0
    END DO
  END DO
! In only variables

  DO j=1,n_dp
    bland_dp(j)           = bland(dpi(j))
    ntml_dp(j)            = ntml(dpi(j))
    ntpar_dp(j)           = ntpar(dpi(j))
    pstar_dp(j)           = pstar(dpi(j))
    recip_pstar_dp(j)     = recip_pstar(dpi(j))
    q1_sd_dp(j)           = q1_sd(dpi(j))
    t1_sd_dp(j)           = t1_sd(dpi(j))
    uw0_dp(j)             = uw0(dpi(j))
    vw0_dp(j)             = vw0(dpi(j))
    w_max_dp(j)           = w_max(dpi(j))
    wstar_dp(j)           = wstar(dpi(j))
    zlcl_uv_dp(j)         = zlcl_uv(dpi(j))
    freeze_lev_dp(j)      = freeze_lev(dpi(j))
    delthvu_dp(j)         = delthvu_bl(dpi(j))
  END DO

  DO k = 0,nlev
    DO j=1,n_dp
      p_layer_centres_dp(j,k)    = p_layer_centres(dpi(j),k)
      p_layer_boundaries_dp(j,k) = p_layer_boundaries(dpi(j),k)

      exner_layer_centres_dp(j,k)    = exner_layer_centres(dpi(j),k)
      exner_layer_boundaries_dp(j,k) = exner_layer_boundaries(dpi(j),k)
      r_theta_dp(j,k)  = r_theta(dpi(j),k)
    END DO
  END DO

  DO k = 1,nlev
    DO j=1,n_dp
      u_dp(j,k)           = u(dpi(j),k)
      v_dp(j,k)           = v(dpi(j),k)
      th_dp(j,k)          = th(dpi(j),k)
      z_theta_dp(j,k)     = z_theta(dpi(j),k)
      z_rho_dp(j,k)       = z_rho(dpi(j),k)
      r_rho_dp(j,k)       = r_rho(dpi(j),k)
      r2rho_th_dp(j,k)     = r2rho_th(dpi(j),k)
      r2rho_dp(j,k)        = r2rho(dpi(j),k)
      rho_theta_dp(j,k)    = rho_theta(dpi(j),k)
      rho_dp(j,k)          = rho(dpi(j),k)
      dr_across_th_dp(j,k) = dr_across_th(dpi(j),k)
      dr_across_rh_dp(j,k) = dr_across_rh(dpi(j),k)
    END DO
  END DO

! G-R mass flux scheme requires input of specific humidity

  IF (l_mixing_ratio) THEN

    DO k = 1,nlev
      DO j = 1,n_dp
        q_dp(j,k)   = q(dpi(j),k)  /(1.0+qt(dpi(j),k))
        qse_dp(j,k) = qse(dpi(j),k)  ! copy specific value
        qcl_dp(j,k) = qcl(dpi(j),k)/(1.0+qt(dpi(j),k))
        qcf_dp(j,k) = qcf(dpi(j),k)/(1.0+qt(dpi(j),k))
      END DO
    END DO

  ELSE   ! Input is specific humidity therefore no problems

    DO k = 1,nlev
      DO j = 1,n_dp
        q_dp(j,k)   = q(dpi(j),k)
        qse_dp(j,k) = qse(dpi(j),k) ! copy specific value
        qcl_dp(j,k) = qcl(dpi(j),k)
        qcf_dp(j,k) = qcf(dpi(j),k)
      END DO
    END DO

  END IF ! (l_mixing_ratio)


! In/out variables


  DO k = 1,nlev
    DO j=1,n_dp
      cf_liquid_dp(j,k)   = cf_liquid(dpi(j),k)
      cf_frozen_dp(j,k)   = cf_frozen(dpi(j),k)
      bulk_cf_dp(j,k)     = bulk_cf(dpi(j),k)
    END DO
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,trlev
        DO j=1,n_dp
          tracer_dp(j,k,ktra)  = tracer(dpi(j),k,ktra)
        END DO
      END DO
      DO k = 1,nlev
        DO j=1,n_dp
          dtrabydt_dp(j,k,ktra)  = 0.0
        END DO
      END DO
    END DO
  END IF


!-----------------------------------------------------------------------
! 1.2 Call deep convection code
!-----------------------------------------------------------------------
! DEPENDS ON: deep_conv_4a
  CALL deep_conv_4a(                                              &
                 !IN
                 nbl,nlev,ntra,n_cca_lev,n_dp,trlev,              &
                 bland_dp, delthvu_dp, exner_layer_centres_dp,    &
                 exner_layer_boundaries_dp,                       &
                 l_calc_dxek, l_q_interact,                       &
                 l_tracer, ntml_dp, ntpar_dp,                     &
                 pstar_dp,p_layer_centres_dp,                     &
                 p_layer_boundaries_dp,                           &
                 z_theta_dp, z_rho_dp, r_theta_dp, r_rho_dp,      &
                 rho_theta_dp,rho_dp,r2rho_th_dp,r2rho_dp,        &
                 dr_across_th_dp,dr_across_rh_dp,                 &
                 q_dp,q1_sd_dp,                                   &
                 t1_sd_dp,th_dp,timestep,u_dp,v_dp,               &
                 uw0_dp, vw0_dp, w_max_dp, wstar_dp,              &
                 zlcl_uv_dp ,freeze_lev_dp,                       &
                 recip_pstar_dp,qse_dp,                           &
                 !INOUT
                 bulk_cf_dp,cf_frozen_dp,cf_liquid_dp,            &
                 qcf_dp,qcl_dp,tracer_dp,                         &
                 !OUT
                 cape_out_dp,cclwp_dp,ccw_dp,cca_dp,              &
                 dbcfbydt_dp,dcffbydt_dp,dcflbydt_dp,             &
                 dqbydt_dp,dqcfbydt_dp,dqclbydt_dp,               &
                 dthbydt_dp,dubydt_dp,dvbydt_dp,dtrabydt_dp,      &
                 detrain_up_dp,detrain_dwn_dp,entrain_up_dp,      &
                 entrain_dwn_dp,                                  &
                 iccb_dp,icct_dp,lcca_dp,lcbase_dp,lctop_dp,      &
                 rain_dp, snow_dp, rain_3d_dp, snow_3d_dp,        &
                 up_flux_dp, up_flux_half_dp,                     &
                 dwn_flux_dp,uw_deep_dp,vw_deep_dp,kterm_dp,      &
                 tcw_dp,cca_2d_dp,                                &
                 rbuoy_p_out_dp,the_out_dp,thp_out_dp,            &
                 qe_out_dp,qp_out_dp,ind_cape_reduced_dp,         &
                 cape_ts_used_dp,cfl_limited_dp,ind_deep_dp,      &
                 error_point                                      &
                 )
  IF (error_point /= 0) THEN
     errorstatus = 2   ! will cause model to fail 
     write(cmessage,'(a37,i12,a8,i3,a9,i2)')                      &
     'Deep conv went to model top at point ',                     &
      dpi(error_point),' in seg ',seg_num,' on call ',call_number

      CALL ereport(routinename, errorstatus, cmessage )

  END IF   

!-----------------------------------------------------------------------
! 1.3 Write data from deep convection points to full arrays
!-----------------------------------------------------------------------

  DO i = 1,n_dp
    cape_out(dpi(i))           = cape_out_dp(i)
    cclwp(dpi(i))              = cclwp_dp(i)
    iccb(dpi(i))               = iccb_dp(i)
    icct(dpi(i))               = icct_dp(i)

    lcca(dpi(i))               = lcca_dp(i)
    lcbase(dpi(i))             = lcbase_dp(i)
    lctop(dpi(i))              = lctop_dp(i)

    rain(dpi(i))               = rain_dp(i)
    snow(dpi(i))               = snow_dp(i)
    precip_deep(dpi(i))        = rain_dp(i) + snow_dp(i)
    kterm_deep(dpi(i))         = kterm_dp(i)
    ind_cape_reduced(dpi(i))   = ind_cape_reduced_dp(i)
    cape_ts_used(dpi(i))       = cape_ts_used_dp(i)
    deep_cfl_limited(dpi(i))   = cfl_limited_dp(i)
    ind_deep(dpi(i))           = ind_deep_dp(i)

    tcw   (dpi(i)) = tcw_dp   (i)
    cca_2d(dpi(i)) = cca_2d_dp(i)
  END DO

  IF (l_mom) THEN
    DO k=1,nlev+1
      DO i = 1,n_dp
        dubydt(dpi(i),k)         = dubydt_dp(i,k)
        dvbydt(dpi(i),k)         = dvbydt_dp(i,k)
      END DO
    END DO
  END IF

  DO k = 1,nlev
    DO i = 1,n_dp
      dthbydt(dpi(i),k)        = dthbydt_dp(i,k)
      dcflbydt(dpi(i),k)       = dcflbydt_dp(i,k)
      dcffbydt(dpi(i),k)       = dcffbydt_dp(i,k)
      dbcfbydt(dpi(i),k)       = dbcfbydt_dp(i,k)
      ccw(dpi(i),k)            = ccw_dp(i,k)
    END DO
  END DO

  DO k=1, n_cca_lev
    DO i=1, n_dp
      cca(dpi(i),k) = cca_dp(i,k)
    END DO
  END DO


  IF (l_ccrad) THEN

    IF (l_dcpl_cld4pc2) THEN
    ! Then scaling is done by CCRad, even if PC2 = .true.
    ! To zero CCA (section 0) for PC2 use the CCRAD knobs
    ! cca_knobs should be set to 0.0 in UMUI for PC2

    ! Use of Anvils will affect both section 0/5 diags if CCRad knobs
    ! are not set to 0.0, this may required further attention in future
      DO k=1, nlev
        DO i=1, n_dp
          ccw0(dpi(i),k) = ccw_dp(i,k) * ccw_dp_knob
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_dp
          cca0(dpi(i),k) = cca_dp(i,k) * cca_dp_knob
        END DO
      END DO

    ELSE

      DO k=1, nlev
        DO i=1, n_dp
          ccw0(dpi(i),k) = ccw_dp(i,k)
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_dp
          cca0(dpi(i),k) = cca_dp(i,k)
        END DO
      END DO
    END IF

    DO i=1, n_dp
      cclwp0 (dpi(i)) = cclwp  (dpi(i))
      iccb0  (dpi(i)) = iccb   (dpi(i))
      icct0  (dpi(i)) = icct   (dpi(i))
      lcbase0(dpi(i)) = lcbase (dpi(i))
!         cca0_2d(dpi(i)) = cca0_2d(dpi(i))
    END DO

  ELSE ! Non-CCRad Code

    IF (l_q_interact) THEN
    ! PC2 Code (PC2 has no deep contribution)
    ELSE
    ! Non-PC2 Code
      DO i=1, n_dp
        cclwp0 (dpi(i)) = cclwp (dpi(i))
        iccb0  (dpi(i)) = iccb  (dpi(i))
        icct0  (dpi(i)) = icct  (dpi(i))
        lcbase0(dpi(i)) = lcbase(dpi(i))
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_dp
          cca0(dpi(i),k) = cca_dp(i,k)
        END DO
      END DO
    END IF

  END IF        ! l_ccrad




! G-R outputs specific humidity

  IF (l_mixing_ratio) THEN  ! Requires conversion

    IF (l_q_interact) THEN  ! PC2

      DO k = 1,nlev
        DO i = 1,n_dp
          dqtt = dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
          dqt  = dqtt*timestep
          denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
          dqbydt(dpi(i),k)   =  denom *                               &
               ( dqbydt_dp(i,k)*(1.0-qt(dpi(i),k))+q_dp(i,k)*dqtt )
          dqclbydt(dpi(i),k) =  denom *                               &
             ( dqclbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcl_dp(i,k)*dqtt )
          dqcfbydt(dpi(i),k) =  denom *                               &
             ( dqcfbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcf_dp(i,k)*dqtt )
        END DO
      END DO

    ELSE                    ! Not PC2
                          ! No qcl and qcf increments anyway
      DO k = 1,nlev
        DO i = 1,n_dp
          dqt   = dqbydt_dp(i,k)*timestep
          denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
          dqbydt(dpi(i),k)    =  dqbydt_dp(i,k) *denom
          dqclbydt(dpi(i),k)  = 0.0
          dqcfbydt(dpi(i),k)  = 0.0
        END DO
      END DO

    END IF                  ! test on PC2

  ELSE        ! output is specific humidity therefore no problems

    DO k = 1,nlev
      DO i = 1,n_dp
        dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
        dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
        dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
      END DO
    END DO

  END IF      ! Test on l_mixing_ratio


  IF (flg_up_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        up_flux(dpi(i),k)        = up_flux_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_up_flx_half) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        up_flux_half(dpi(i),k)        = up_flux_half_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_dwn_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        dwn_flux(dpi(i),k)       = dwn_flux_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        entrain_up(dpi(i),k)     = entrain_up_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        detrain_up(dpi(i),k)     = detrain_up_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        entrain_dwn(dpi(i),k)    = entrain_dwn_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        detrain_dwn(dpi(i),k)    = detrain_dwn_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_uw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        uw_deep(dpi(i),k)        = uw_deep_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_vw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        vw_deep(dpi(i),k)        = vw_deep_dp(i,k)
      END DO
    END DO
  END IF

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,nlev
        DO i = 1, n_dp
          dtrabydt(dpi(i),k,ktra)  = dtrabydt_dp(i,k,ktra)
        END DO
      END DO
    END DO
  END IF

! Merge dp 3d rain & snow profiles

  DO k=1, nlev
    DO i=1, n_dp
      rain_3d(dpi(i),k) = rain_3d_dp(i,k)
      snow_3d(dpi(i),k) = snow_3d_dp(i,k)
    END DO
  END DO

END IF !n_dp > 0


!-----------------------------------------------------------------------
! 2.0 SHALLOW convection
! 2.1 Compress input variable arrays for shallow convection scheme to
!     length n_sh (shallow points only)
!-----------------------------------------------------------------------
IF (n_sh  >   0) THEN
  j = 0
  DO i = 1,npnts
    IF (cumulus_bl(i).AND.l_shallow_bl(i)) THEN
      j                        = j+1
      shi(j)                   = i
    END IF
  END DO

! initialise 3d rain varibles

  DO k=1, nlev
    DO i=1, n_sh
      rain_3d_sh(i,k) = 0.0
      snow_3d_sh(i,k) = 0.0
    END DO
  END DO

! In only variables

  DO j = 1,n_sh
    bland_sh(j)           = bland(shi(j))
    delthvu_sh(j)         = delthvu_bl(shi(j))
    ntml_sh(j)            = ntml(shi(j))
    ntpar_sh(j)           = ntpar(shi(j))
    pstar_sh(j)           = pstar(shi(j))
    recip_pstar_sh(j)     = recip_pstar(shi(j))
    q1_sd_sh(j)           = q1_sd(shi(j))
    t1_sd_sh(j)           = t1_sd(shi(j))
    uw0_sh(j)             = uw0(shi(j))
    vw0_sh(j)             = vw0(shi(j))
    wstar_sh(j)           = wstar(shi(j))
    wthvs_sh(j)           = wthvs(shi(j))
    zlcl_uv_sh(j)         = zlcl_uv(shi(j))
    ztop_uv_sh(j)         = ztop_uv(shi(j))
    freeze_lev_sh(j)      = freeze_lev(shi(j))
  END DO

  DO k = 0,nlev
    DO j = 1,n_sh
      p_layer_centres_sh(j,k)        = p_layer_centres(shi(j),k)
      p_layer_boundaries_sh(j,k)     = p_layer_boundaries(shi(j),k)
      exner_layer_centres_sh(j,k)    = exner_layer_centres(shi(j),k)
      exner_layer_boundaries_sh(j,k) = exner_layer_boundaries(shi(j),k)
    END DO
  END DO

  DO k = 1,nlev
    DO j = 1,n_sh
      u_sh(j,k)           = u(shi(j),k)
      v_sh(j,k)           = v(shi(j),k)
      th_sh(j,k)          = th(shi(j),k)
      z_theta_sh(j,k)     = z_theta(shi(j),k)
      z_rho_sh(j,k)       = z_rho(shi(j),k)
    END DO
  END DO

! Moisture - scheme requires input of specific humidity

  IF (l_mixing_ratio) THEN
!  Conversion required
    DO k = 1,nlev
      DO j = 1,n_sh
        q_sh(j,k)   = q(shi(j),k)  /(1.0+qt(shi(j),k))
        qse_sh(j,k) = qse(shi(j),k)  ! copy specific value
        qcl_sh(j,k) = qcl(shi(j),k)/(1.0+qt(shi(j),k))
        qcf_sh(j,k) = qcf(shi(j),k)/(1.0+qt(shi(j),k))
      END DO
    END DO

  ELSE     ! Input is specific humidity therefore no problems

    DO k = 1,nlev
      DO j = 1,n_sh
        q_sh(j,k)   = q(shi(j),k)
        qse_sh(j,k) = qse(shi(j),k) ! copy specific value
        qcl_sh(j,k) = qcl(shi(j),k)
        qcf_sh(j,k) = qcf(shi(j),k)
      END DO
    END DO

  END IF   ! Test on l_mixing_ratio

! In/out variables


  DO k = 1,nlev
    DO j = 1,n_sh
      cf_liquid_sh(j,k)   = cf_liquid(shi(j),k)
      cf_frozen_sh(j,k)   = cf_frozen(shi(j),k)
      bulk_cf_sh(j,k)     = bulk_cf(shi(j),k)
    END DO
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,trlev
        DO j = 1,n_sh
          tracer_sh(j,k,ktra)  = tracer(shi(j),k,ktra)
        END DO
      END DO
      DO k = 1,nlev
        DO j = 1,n_sh
          dtrabydt_sh(j,k,ktra)  = 0.0
        END DO
      END DO
    END DO
  END IF


!-----------------------------------------------------------------------
! 2.2 Call shallow convection code
!-----------------------------------------------------------------------
! DEPENDS ON: shallow_conv_4a
  CALL shallow_conv_4a(                                           &
                     !IN
                     nbl,nlev,ntra,n_cca_lev,n_sh,trlev,          &
                     bland_sh,delthvu_sh,exner_layer_centres_sh,  &
                     exner_layer_boundaries_sh,                   &
                     l_calc_dxek,                                 &
                     l_q_interact, l_tracer, ntml_sh, ntpar_sh,   &
                     pstar_sh,p_layer_centres_sh,                 &
                     p_layer_boundaries_sh,z_theta_sh,z_rho_sh,   &
                     r2rho_th_sh, dr_across_th_sh,                &
                     q_sh,q1_sd_sh,t1_sd_sh,th_sh,timestep,       &
                     u_sh,v_sh,uw0_sh,vw0_sh,wstar_sh,wthvs_sh,   &
                     zlcl_uv_sh,ztop_uv_sh,freeze_lev_sh,         &
                     recip_pstar_sh,qse_sh,                       &
                     !INOUT
                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh,        &
                     qcf_sh,qcl_sh,tracer_sh,                     &
                     !OUT
                     cape_out_sh,cclwp_sh,ccw_sh,cca_sh,          &
                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh,         &
                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh,           &
                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh,  &
                     detrain_up_sh,detrain_dwn_sh,entrain_up_sh,  &
                     entrain_dwn_sh,                              &
                     iccb_sh,icct_sh,                             &
                     lcca_sh,lcbase_sh,lctop_sh,                  &
                     rain_sh, snow_sh, rain_3d_sh, snow_3d_sh,    &
                     up_flux_sh, up_flux_half_sh,                 &
                     dwn_flux_sh,uw_shall_sh,vw_shall_sh,         &
                     tcw_sh,cca_2d_sh                             &
                     ,rbuoy_p_out_sh,the_out_sh,thp_out_sh        &
                     ,qe_out_sh,qp_out_sh                         &
                     )

!-----------------------------------------------------------------------
! 2.3 Write data from shallow convection points to full arrays
!-----------------------------------------------------------------------

  DO i = 1,n_sh
    cape_out(shi(i))           = cape_out_sh(i)
    cclwp(shi(i))              = cclwp_sh(i)
    iccb(shi(i))               = iccb_sh(i)
    icct(shi(i))               = icct_sh(i)

    lcca(shi(i))               = lcca_sh(i)
    lcbase(shi(i))             = lcbase_sh(i)
    lctop(shi(i))              = lctop_sh(i)
    rain(shi(i))               = rain_sh(i)
    snow(shi(i))               = snow_sh(i)
    precip_shall(shi(i))       = rain_sh(i) + snow_sh(i)

    tcw   (shi(i)) = tcw_sh   (i)
    cca_2d(shi(i)) = cca_2d_sh(i)
    ind_shall(shi(i)) = 1.0       ! all shall points 
  END DO

  IF (l_mom) THEN
    DO k=1,nlev+1
      DO i = 1,n_sh
        dubydt(shi(i),k)         = dubydt_sh(i,k)
        dvbydt(shi(i),k)         = dvbydt_sh(i,k)
      END DO
    END DO
  END IF

  DO k = 1,nlev
    DO i = 1,n_sh
      dthbydt(shi(i),k)        = dthbydt_sh(i,k)
      dcflbydt(shi(i),k)       = dcflbydt_sh(i,k)
      dcffbydt(shi(i),k)       = dcffbydt_sh(i,k)
      dbcfbydt(shi(i),k)       = dbcfbydt_sh(i,k)
      ccw(shi(i),k)            = ccw_sh(i,k)
    END DO
  END DO

  DO k=1, n_cca_lev
    DO i=1, n_sh
      cca(shi(i),k) = cca_sh(i,k)
    END DO
  END DO


  IF (l_ccrad) THEN

    DO i=1, n_sh
      iccb0  (shi(i)) = iccb  (shi(i))
      icct0  (shi(i)) = icct  (shi(i))
      cclwp0 (shi(i)) = cclwp (shi(i))
      lcbase0(shi(i)) = lcbase(shi(i))
    END DO

    IF (l_dcpl_cld4pc2) THEN

      DO k=1, nlev
        DO i=1, n_sh
          ccw0(shi(i),k) = ccw_sh(i,k) * ccw_sh_knob
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_sh
          cca0(shi(i),k) = cca_sh(i,k) * cca_sh_knob
        END DO
      END DO

    ELSE

      DO k=1, nlev
        DO i=1, n_sh
          ccw0(shi(i),k) = ccw_sh(i,k)
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_sh
          cca0(shi(i),k) = cca_sh(i,k)
        END DO
      END DO

    END IF  ! l_dcpl_cld4pc2

  ELSE ! Non-ccrad code

    IF ( l_q_interact .AND. (.NOT. l_pc2_diag_sh) ) THEN
    ! Do nothing as PC2 only keeps shallow component
    ! if l_pc2_diag_sh=.TRUE.
    ELSE
    ! Either PC2=.FALSE. OR PC2=.TRUE. and it wants
    ! Shallow cloud

      DO i=1, n_sh
        iccb0  (shi(i)) = iccb  (shi(i))
        icct0  (shi(i)) = icct  (shi(i))
        cclwp0 (shi(i)) = cclwp (shi(i))
        lcbase0(shi(i)) = lcbase(shi(i))
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_sh
          cca0(shi(i),k) = cca_sh(i,k)
        END DO
      END DO

!         ! Although ccw0 is not passed up to atm_step
!         ! because (l_ccrad=.false.), shallow points
!         ! are kept because they may needed in ni_imp_ctl
!         ! by l_pc2_diag_sh. This needs, needs looking into/possibly
!         ! fixing.

!         DO k=1, nlev
!           DO i=1, n_sh
!             ccw0(shi(i),k) = ccw_sh(i,k)
!           END DO
!         END DO
!
!         DO i=1, n_sh
!           cca0_2d(shi(i)) = cca_2d_sh(i)
!         END DO

    END IF

  END IF ! l_ccrad




! G-R requires input of specific humidity

  IF (l_mixing_ratio) THEN ! requires conversion

    IF (l_q_interact) THEN  ! PC2

      DO k = 1,nlev
        DO i = 1,n_sh
          dqtt  = dqbydt_sh(i,k)+dqclbydt_sh(i,k)+dqcfbydt_sh(i,k)
          dqt   = dqtt*timestep
          denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
          dqbydt(shi(i),k)   = denom *                                   &
               ( dqbydt_sh(i,k)*(1.0-qt(shi(i),k))+q_sh(i,k)*dqtt )
          dqclbydt(shi(i),k) = denom *                                   &
             ( dqclbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcl_sh(i,k)*dqtt )
          dqcfbydt(shi(i),k) = denom *                                   &
             ( dqcfbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcf_sh(i,k)*dqtt )
        END DO
      END DO

    ELSE                    ! Not PC2
                         ! No qcl and qcf increments anyway
      DO k = 1,nlev
        DO i = 1,n_sh
          dqt   = dqbydt_sh(i,k)*timestep
          denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
          dqbydt(shi(i),k)    = dqbydt_sh(i,k) *denom
          dqclbydt(shi(i),k)  = 0.0
          dqcfbydt(shi(i),k)  = 0.0
        END DO
      END DO

    END IF                  ! end test on PC2

  ELSE   ! output is specific humidity therefore no problems

    DO k = 1,nlev
      DO i = 1,n_sh
        dqbydt(shi(i),k)   = dqbydt_sh(i,k)
        dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
        dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
      END DO
    END DO

  END IF      ! Test on l_mixing_ratio

  IF (flg_up_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        up_flux(shi(i),k)        = up_flux_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_up_flx_half) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        up_flux_half(shi(i),k)   = up_flux_half_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_dwn_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        dwn_flux(shi(i),k)       = dwn_flux_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        entrain_up(shi(i),k)     = entrain_up_sh(i,k)
      END DO
    END DO 
  END IF
  IF (flg_detr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        detrain_up(shi(i),k)     = detrain_up_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        entrain_dwn(shi(i),k)    = entrain_dwn_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        detrain_dwn(shi(i),k)    = detrain_dwn_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_uw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        uw_shall(shi(i),k)       = uw_shall_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k = 1,nlev
      DO i = 1,n_sh
        vw_shall(shi(i),k)       = vw_shall_sh(i,k)
      END DO
    END DO
  END IF

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,nlev
        DO i = 1,n_sh
          dtrabydt(shi(i),k,ktra)  = dtrabydt_sh(i,k,ktra)
        END DO
      END DO
    END DO
  END IF

! Merge sh 3d rain & snow profiles

  DO k=1,nlev
    DO i=1,n_sh
      rain_3d(shi(i),k) = rain_3d_sh(i,k)
      snow_3d(shi(i),k) = snow_3d_sh(i,k)
    END DO
  END DO

END IF     ! n_sh > 0



!-----------------------------------------------------------------------
! 3.0 MID-LEVEL Convection
! 3.1 Set lowest level that mid level convection can trigger
!-----------------------------------------------------------------------
IF (n_md > 0) THEN

  DO k=1, nlev
    DO i=1, n_md
      rain_3d_md(i,k) = 0.0
      snow_3d_md(i,k) = 0.0
    END DO
  END DO

  j = 0
  DO i = 1,npnts
    IF (l_mid(i)) THEN
      j      = j+1
      mdi(j) = i
    END IF
  END DO


  idp=0
  DO i = 1,n_md
    IF (.NOT. cumulus_bl(mdi(i))) THEN
      midtrig(i) = ntml(mdi(i)) + 1
      IF (ntml(mdi(i))  ==  nbl) THEN
        midtrig(i) = ntml(mdi(i))
      END IF
    ELSE

!  Cumulus points ie deep or shallow convection has occurred

      midtrig(i) = ntpar(mdi(i)) + 2


! NTPAR has a maximum value which can be less than the top level for
! deep. Deep convection may terminate above this. Need to prevent mid
! level convection occurring at the same levels as deep.

      IF(.NOT.l_shallow_bl(mdi(i))) THEN   ! deep points
        IF (kterm_deep(mdi(i)) >  ntpar(mdi(i))+2) THEN
          midtrig(i) = kterm_deep(mdi(i))
        END IF
      END IF  ! deep points
    END IF
  END DO

!-----------------------------------------------------------------------
! 3.2 Copy all input arrays to arrays ending in _md for passing to
!     mid-level scheme
!-----------------------------------------------------------------------

! In only variables

  DO j=1,n_md
    bland_md(j)           = bland(mdi(j))
    ntml_md(j)            = ntml(mdi(j))
    ntpar_md(j)           = ntpar(mdi(j))
    pstar_md(j)           = pstar(mdi(j))
    recip_pstar_md(j)     = recip_pstar(mdi(j))
    q1_sd_md(j)           = q1_sd(mdi(j))
    t1_sd_md(j)           = t1_sd(mdi(j))
    w_max_md(j)           = w_max(mdi(j))
    freeze_lev_md(j)      = freeze_lev(mdi(j))
  END DO

  IF (.NOT. l_ccrad) THEN
    DO i=1, n_md
      tcw_md(i)    = tcw(mdi(i))
      cca_2d_md(i) = cca_2d(mdi(i))
    END DO
  END IF

  DO i=1, n_md
    iccb_md(i) = iccb(mdi(i))
    icct_md(i) = icct(mdi(i))
  END DO

  DO k = 0,nlev
    DO i = 1,n_md
      p_layer_centres_md(i,k)      = p_layer_centres(mdi(i),k)
      p_layer_boundaries_md(i,k)   = p_layer_boundaries(mdi(i),k)
      exner_layer_centres_md(i,k)  = exner_layer_centres(mdi(i),k)
      exner_layer_boundaries_md(i,k) = exner_layer_boundaries(mdi(i),k)
      r_theta_md(i,k)              = r_theta(mdi(i),k)
    END DO
  END DO
  DO k = 1,nlev
    DO i = 1,n_md
      u_md(i,k)           = u(mdi(i),k)
      v_md(i,k)           = v(mdi(i),k)
      th_md(i,k)          = th(mdi(i),k)
      rho_md(i,k)         = rho(mdi(i),k)
      rho_theta_md(i,k)   = rho_theta(mdi(i),k)
      z_theta_md(i,k)     = z_theta(mdi(i),k)
      z_rho_md(i,k)       = z_rho(mdi(i),k)
      r_rho_md(i,k)       = r_rho(mdi(i),k)
      cf_liquid_md(i,k)   = cf_liquid(mdi(i),k)
      cf_frozen_md(i,k)   = cf_frozen(mdi(i),k)
      bulk_cf_md(i,k)     = bulk_cf(mdi(i),k)
    END DO
  END DO


! Moisture - scheme requires input of specific humidity

  IF (l_mixing_ratio) THEN
  !  Conversion required
    DO k = 1,nlev
      DO i = 1,n_md
        q_md(i,k)   = q(mdi(i),k)  /(1.0+qt(mdi(i),k))
        qse_md(i,k) = qse(mdi(i),k)  ! copy specific value
        qcl_md(i,k) = qcl(mdi(i),k)/(1.0+qt(mdi(i),k))
        qcf_md(i,k) = qcf(mdi(i),k)/(1.0+qt(mdi(i),k))
      END DO
    END DO

  ELSE       ! input is specific humidity therefore no problems

    DO k = 1,nlev
      DO i = 1,n_md
        q_md(i,k)   = q(mdi(i),k)
        qse_md(i,k) = qse(mdi(i),k)  ! copy specific value
        qcl_md(i,k) = qcl(mdi(i),k)
        qcf_md(i,k) = qcf(mdi(i),k)
      END DO
    END DO

  END IF     ! Test on l_mixing_ratio

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,trlev
        DO i = 1,n_md
          tracer_md(i,k,ktra)  = tracer(mdi(i),k,ktra)
        END DO
      END DO
      DO k = 1,nlev
        DO i = 1,n_md
          dtrabydt_md(i,k,ktra)  = 0.0
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
! 3.3 Call mid-level convection code
!-----------------------------------------------------------------------

! DEPENDS ON: mid_conv_4a
  CALL mid_conv_4a(                                               &
                !IN
                nbl,nlev,ntra,n_cca_lev,n_md,trlev,               &
                bland_md,w_max_md,exner_layer_centres_md,         &
                exner_layer_boundaries_md,                        &
                l_calc_dxek, l_q_interact,                        &
                l_tracer, midtrig, ntml_md, ntpar_md,             &
                freeze_lev_md,pstar_md, p_layer_centres_md,       &
                p_layer_boundaries_md,                            &
                r_theta_md, r_rho_md,                             &
                z_theta_md, z_rho_md, rho_md, rho_theta_md,       &
                q_md, q1_sd_md, t1_sd_md,                         &
                th_md,timestep,u_md,v_md,recip_pstar_md,qse_md,   &
                !INOUT
                bulk_cf_md,cf_frozen_md,cf_liquid_md,             &
                qcf_md,qcl_md,tracer_md,                          &
                !OUT
                cape_out_md,cclwp_md,ccw_md,cca_md,               &
                dbcfbydt_md,dcffbydt_md,dcflbydt_md,              &
                dqbydt_md,dqcfbydt_md,dqclbydt_md,                &
                dthbydt_md,dubydt_md,dvbydt_md,dtrabydt_md,       &
                detrain_up_md,detrain_dwn_md,entrain_up_md,       &
                entrain_dwn_md,                                   &
                iccb_md,icct_md,                                  &
                lcca_md,lcbase_md,lctop_md,                       &
                rain_md, snow_md, rain_3d_md, snow_3d_md,         &
                up_flux_md, up_flux_half_md,                      &
                dwn_flux_md,tcw_md,l_mid_md,cca_2d_md,            &
                rbuoy_p_out_md,the_out_md,thp_out_md,             &
                qe_out_md,qp_out_md,                              &
                uw_mid, vw_mid, cfl_limited_md, error_point       &
                )

  IF (error_point /= 0) THEN
    errorstatus = 3   ! will cause model to fail 
    write(cmessage,'(a47,i12,a16,i2)')                                  &
    'Mid conv went to the top of the model at point ',                  &
    mdi(error_point),' in seg on call ',call_number

    CALL ereport(routinename, errorstatus, cmessage )

  END IF

!-----------------------------------------------------------------------
! 3.4 Write data from mid-level convection to full arrays
!-----------------------------------------------------------------------

  DO i = 1,n_md

! Cloud variables - only overwrite deep or shallow values if
! iccb_md  and icct_md > 0

    IF (iccb_md(i)  >   0 .AND. icct_md(i)  >   0) THEN
      iccb(mdi(i))  = iccb_md(i)
      icct(mdi(i))  = icct_md(i)
    END IF

! Overwrite lowest cloud values only if cumulus = .F.

    IF (.NOT. cumulus_bl(mdi(i))) THEN
      lcca(mdi(i))   = lcca_md(i)
      lcbase(mdi(i)) = lcbase_md(i)
      lctop(mdi(i))  = lctop_md(i)
    END IF

! Write remaining data to full arrays

    cape_out(mdi(i))   = cape_out(mdi(i)) + cape_out_md(i)
    cclwp(mdi(i))      = cclwp(mdi(i))    + cclwp_md(i)
    rain(mdi(i))       = rain(mdi(i))     + rain_md(i)
    snow(mdi(i))       = snow(mdi(i))     + snow_md(i)
    precip_mid(mdi(i)) = rain_md(i)       + snow_md(i)
    mid_cfl_limited(mdi(i)) = cfl_limited_md(i)

    l_mid_all(mdi(i))  = l_mid_md(i)

  !=====================================================================
  ! NOTE: At this point if l_ccrad = T, then cca_2d_md is ONLY equal
  !       to that from mid-level cloud on a given grid point.  The
  !       original code I.E. l_ccrad = F means that cca_2d_md will
  !       include that from sh/dp aswell.
  !=====================================================================

  END DO


! Merge md 3d rain & snow profiles
  DO k=1,nlev
    DO i=1,n_md
      rain_3d(mdi(i),k) = rain_3d(mdi(i),k) + rain_3d_md(i,k)
      snow_3d(mdi(i),k) = snow_3d(mdi(i),k) + snow_3d_md(i,k)
    END DO
  END DO


  IF (l_mom) THEN
    DO k=1,nlev+1
      DO i = 1,n_md
        dubydt(mdi(i),k) = dubydt(mdi(i),k) + dubydt_md(i,k)
        dvbydt(mdi(i),k) = dvbydt(mdi(i),k) + dvbydt_md(i,k)
      END DO
    END DO
  END IF

  DO k = 1,nlev
    DO i = 1,n_md
      dthbydt(mdi(i),k)  = dthbydt(mdi(i),k)  + dthbydt_md(i,k)
      dcflbydt(mdi(i),k) = dcflbydt(mdi(i),k) + dcflbydt_md(i,k)
      dcffbydt(mdi(i),k) = dcffbydt(mdi(i),k) + dcffbydt_md(i,k)
      dbcfbydt(mdi(i),k) = dbcfbydt(mdi(i),k) + dbcfbydt_md(i,k)
      ccw(mdi(i),k)      = ccw(mdi(i),k)      + ccw_md(i,k)
    END DO
  END DO


  DO k=1, n_cca_lev
    DO i=1, n_md
      cca(mdi(i),k) = cca(mdi(i),k) + cca_md(i,k)
    END DO
  END DO


  IF (l_ccrad) THEN

    DO i=1, n_md
      cclwp0 (mdi(i)) = cclwp (mdi(i))
      iccb0  (mdi(i)) = iccb  (mdi(i))
      icct0  (mdi(i)) = icct  (mdi(i))
      lcbase0(mdi(i)) = lcbase(mdi(i))
    END DO

    IF (l_dcpl_cld4pc2) THEN
    ! Then scaling is done by CCRad, even if PC2 = .true.
    ! To zero CCA (section 0) for PC2 use the CCRAD knobs
    ! cca_knobs should be set to 0.0 in umui for PC2

    ! Anvils if selected will affect both section 0/5 diags
    ! if ccrads are not set to 0.0, this may required further
    ! attention in future
      DO k=1, nlev
        DO i=1, n_md
          ccw0(mdi(i),k) = ccw0(mdi(i),k) + ccw_md(i,k) * ccw_md_knob
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_md
          cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k) * cca_md_knob
        END DO
      END DO

    ELSE

      DO k=1, nlev
        DO i=1, n_md
          ccw0(mdi(i),k) = ccw0(mdi(i),k) + ccw_md(i,k)
        END DO
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_md
          cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k)
        END DO
      END DO

    END IF ! l_dcpl_cld4pc2


  ! If l_ccrad=T add on cca_2d_sh,cca_2d_dp to cca_2d.
  ! Contributions from cca_2d_sh/cca_2d_dp have already
  ! been applied to cca_sh/cca_dp

  ! However, contributions from cca_2d_sh/cca_2d_dp need
  ! to be included in the cca_2d diagnostic so as not to
  ! upset any Downstream products that use it.
  ! (at this cca_2d point only holds cca_2d_md)

  ! cca_2d_md entered mid_conv as an empty array,
  ! i.e. no contribution from shallow/deep.
  ! Add any contribution from mid level to preserve
  ! CCA_2d diagnostic.

    DO i=1, n_md
      cca_2d(mdi(i)) = cca_2d(mdi(i)) + cca_2d_md(i)
    END DO

  ELSE ! Non-CCRad

    IF (l_q_interact) THEN
    ! PC2 Code (PC2 has no mid-level contribution)
    ELSE
    ! Non-PC2 Code
      DO i=1, n_md
        cclwp0 (mdi(i)) = cclwp (mdi(i))
        iccb0  (mdi(i)) = iccb  (mdi(i))
        icct0  (mdi(i)) = icct  (mdi(i))
        lcbase0(mdi(i)) = lcbase(mdi(i))
      END DO

      DO k=1, n_cca_lev
        DO i=1, n_md
          cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k)
        END DO
      END DO
    END IF

  ! In non-ccrad configs, the cca_2d_md array is not empty when
  ! entering mid-conv; it holds contributions from deep_conv/
  ! shallow_conv, so we can just copy cca_2d_md back to cca_2d
    DO i=1, n_md
      cca_2d(mdi(i)) = cca_2d_md(i)
    END DO

  END IF ! l_ccrad





! G-R outputs specific humidity

  IF (l_mixing_ratio) THEN  ! Requires conversion

    IF (l_q_interact) THEN  ! PC2

      DO k = 1,nlev
        DO i = 1,n_md
          dqtt= dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
          dqt = dqtt*timestep
          denom = 1.0/((1.0-qt(mdi(i),k))*(1.0-qt(mdi(i),k)-dqt))
          dqbydt(mdi(i),k)   =  dqbydt(mdi(i),k) + denom *           &
                     ( dqbydt_md(i,k)*(1.0-qt(mdi(i),k))+ q(mdi(i),k)*dqtt )
          dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + denom *          &
                     ( dqclbydt_md(i,k)*(1.0-qt(mdi(i),k))+ qcl(mdi(i),k)*dqtt )
          dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) +  denom *         &
                     ( dqcfbydt_md(i,k)*(1.0-qt(mdi(i),k))+ qcf(mdi(i),k)*dqtt )
        END DO
      END DO

    ELSE                    ! Not PC2
                            ! No qcl and qcf increments anyway
      DO k = 1,nlev
        DO i = 1,n_md
          dqt   = dqbydt_md(i,k)*timestep
          denom = 1.0/((1.0-qt(mdi(i),k))*(1.0-qt(mdi(i),k)-dqt))
          dqbydt(mdi(i),k)    = dqbydt(mdi(i),k) + dqbydt_md(i,k)*denom
          dqclbydt(mdi(i),k)  = 0.0
          dqcfbydt(mdi(i),k)  = 0.0
        END DO
      END DO

    END IF                  ! test on PC2

  ELSE        ! output is specific humidity therefore no problems

    DO k = 1,nlev
      DO i = 1,n_md
        dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + dqbydt_md(i,k)
        dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + dqclbydt_md(i,k)
        dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + dqcfbydt_md(i,k)
      END DO
    END DO

  END IF      ! Test on l_mixing_ratio


  IF (flg_up_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        up_flux(mdi(i),k)     = up_flux(mdi(i),k) + up_flux_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_up_flx_half) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        up_flux_half(mdi(i),k)= up_flux_half(mdi(i),k) + up_flux_half_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_dwn_flx) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        dwn_flux(mdi(i),k)     = dwn_flux(mdi(i),k) + dwn_flux_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        entrain_up(mdi(i),k)   = entrain_up(mdi(i),k) + entrain_up_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_up) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        detrain_up(mdi(i),k)   = detrain_up(mdi(i),k) + detrain_up_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        entrain_dwn(mdi(i),k)  = entrain_dwn(mdi(i),k) + entrain_dwn_md(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_dwn) THEN
    DO k = 1,nlev
      DO i = 1,n_md
        detrain_dwn(mdi(i),k)  = detrain_dwn(mdi(i),k) + detrain_dwn_md(i,k)
      END DO
    END DO
  END IF
  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,nlev
        DO i = 1,n_md
          dtrabydt(mdi(i),k,ktra)  = dtrabydt(mdi(i),k,ktra)      &
                                    + dtrabydt_md(i,k,ktra)
        END DO
      END DO
    END DO
  END IF

END IF ! (n_md > 0) Test


! ---------------------------------------------------------------------
!write adaptive scm diagnostics
!commented out because current version causes problems with
!substepping - code still available for future development work.
! ---------------------------------------------------------------------
!#if defined(SCMA)
!      sname='rbuoy_p_out'
!      lname='buoyancy excess as used in Parcel'
!      units='K '
!      call SCMoutput(rbuoy_p_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='the_out'
!      lname='th_E as used in Parcel'
!      units='K '
!      call SCMoutput(the_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='thp_out'
!      lname='th_P as used in Parcel'
!      units='K '
!      call SCMoutput(the_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='qe_out'
!      lname='q_E as used in Parcel'
!      units='kg/kg'
!      call SCMoutput(qe_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='qp_out'
!      lname='q_P as used in Parcel'
!      units='kg/kg '
!      call SCMoutput(qp_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')

!#endif

! ---------------------------------------------------------------------
! 4.0 Energy correction calculation - involves integral over
!     whole column therefore this should not be placed in the
!     separate calls to deep, shallow and mid. Some locations will
!     have both shallow and mid level or deep and mid level convection

!     UM documentation paper 27 - section 12.
! ---------------------------------------------------------------------
! First work out which points convection has occurred at .

nconv_all=0
DO i = 1,npnts
  IF (cumulus_bl(i).OR.l_mid_all(i)) THEN
    nconv_all = nconv_all + 1
    index1(nconv_all) = i
  END IF
END DO

IF (nconv_all >  0) THEN

! DEPENDS ON: cor_engy_4a
  CALL cor_engy_4a(np_field,npnts,nconv_all,nlev,dthbydt,dqbydt,snow &
                  ,exner_layer_centres,p_layer_boundaries,index1)

!-----------------------------------------------------------------------
! 5.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.


  DO j = 1,nconv_all
    i=index1(j)
    qminincolumn(j) = q(i,nlev)
  END DO
  DO k = 1,nlev-1
    DO j = 1,nconv_all
      i=index1(j)
      IF (q(i,k)  <   qminincolumn(j)) THEN
        qminincolumn(j) = q(i,k)
      END IF
    END DO
  END DO


! Ensure Q does not go below global allowed minimum (QMIN)

  DO j = 1,nconv_all
    qminincolumn(j)=MAX(qmin,qminincolumn(j))
  END DO


! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.


  DO k = nlev,2,-1
! NEC SX6 compiler directive so that next loop vectorized
!CDIR NODEP
    DO j = 1,nconv_all
      i=index1(j)
      temp1(j)=q(i,k) + dqbydt(i,k) * timestep

      IF (temp1(j)  <   qminincolumn(j)) THEN

        dqbydt(i,k-1) = dqbydt(i,k-1) -                           &
              ((qminincolumn(j) - q(i,k)) / timestep-dqbydt(i,k)) &
               * ( - (p_layer_boundaries(i,k) -                   &
                      p_layer_boundaries(i,k-1)))                 &
               / ( - (p_layer_boundaries(i,k-1) -                 &
                      p_layer_boundaries(i,k-2)))

        dqbydt(i,k) = (qminincolumn(j) - q(i,k)) / timestep
      END IF
    END DO ! nconv_all loop
  END DO  ! nlev

!-----------------------------------------------------------------------
! 6.0  Subroutine MIX_INC mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.
!-----------------------------------------------------------------------
  IF (bl_cnv_mix == 0) THEN
! DEPENDS ON: mix_inc
    CALL mix_inc(np_field,npnts,nconv_all,nlev,nbl,ntml,          &
               dthbydt,dqbydt,dubydt,dvbydt,l_tracer,ntra,        &
               dtrabydt,p_layer_boundaries,                       &
               p_layer_centres,index1)


  END IF
END IF    ! test on nconv_all

!-----------------------------------------------------------------------
! 7.0  Calculate convective cloud amount on model levels if
!      l_3d_cca = .true.
!-----------------------------------------------------------------------
IF (.NOT. l_ccrad) THEN

  DO k=1, n_cca_lev
    DO i=1, npnts
      cca(i,k) = 0.0
    END DO
  END DO

  IF (l_3d_cca) THEN

    IF (l_anvil) THEN

! DEPENDS ON: calc_3d_cca
      CALL calc_3d_cca                                            &
         ( np_field, npnts, nlev, n_cca_lev, nbl, iccb, icct      &
         , p_layer_boundaries, freeze_lev, cca_2d, cca, z_theta   &
         , z_rho )

    ELSE

      DO k=1, n_cca_lev
        DO i=1, npnts
          IF ( k >= iccb(i) .AND. &
               k <= icct(i) ) THEN
            cca(i,k) = cca_2d(i)
          END IF
        END DO
      END DO

    END IF

  ELSE ! CCA is a single level field and l_ccrad is false

    DO i=1, npnts
      cca(i,1) = cca_2d(i)
    END DO

  END IF ! l_3d_cca


!-----------------------------------------------------------------------------
! Section 0 cloud properties
!-----------------------------------------------------------------------------
  IF (l_q_interact) THEN

    IF (l_pc2_diag_sh) THEN
      DO k=1, n_cca_lev
        DO i=1, n_sh
          IF ( k >= iccb_sh(i) .AND. &
               k < icct_sh(i) ) THEN
            cca0(shi(i),k) = cca_2d_sh(i)
          END IF
        END DO
      END DO
    ELSE
      DO k=1, n_cca_lev
        DO i=1, npnts
          cca0(i,k) = 0.0
        END DO
      END DO
    END IF

  ELSE ! NOT l_q_interact

    DO k=1, n_cca_lev
      DO i=1, npnts
        cca0(i,k) = cca(i,k)
      END DO
    END DO

  END IF

END IF ! .NOT. CCRad

!-----------------------------------------------------------------------
! 8.0  Update tracer field
!      More efficient to do here rather than in subroutines.
!      This seems to be an expensive bit of code on NEC.
!      Changed to operate on just convective points.
!-----------------------------------------------------------------------

IF (l_tracer.AND.(nconv_all >  0)) THEN


! Adjust timestep to prevent any negative values invading the tracer
! fields (adjusted timestep is a function of geograhical  location and
! tracer type.

  DO ktra=1,ntra
! initialise step to timestep
    DO i = 1,nconv_all
      limited_step(i) =timestep
    END DO

    DO k=1,nlev
      DO j = 1,nconv_all
        i=index1(j)
! negative increments  may have a problem
        IF (dtrabydt(i,k,ktra)  <   0.0 ) THEN

          step_test(j) = (0.9999 * ABS(tracer(i,k,ktra))) /       &
                    (ABS(dtrabydt(i,k,ktra)) + safety_margin)

          IF (step_test(j)   <   limited_step(j) ) THEN
! then increment is bigger than tracer and timestep needs to be reduced
            limited_step (j) = step_test(j)
          END IF

        END IF
      END DO
    END DO

! Update tracer field using limited_step.

    DO k = 1,nlev
! NEC compiler directive
!CDIR NODEP
      DO j = 1,nconv_all
        i=index1(j)
        tracer(i,k,ktra) = tracer(i,k,ktra) +                       &
                         dtrabydt(i,k,ktra) * limited_step(j)
      END DO
    END DO

  END DO  ! ktra loop
END IF   ! L_tracer

!-----------------------------------------------------------------------


! Check no CCA/CCA is greater than 1.0
DO k=1, n_cca_lev
  DO i=1, npnts
    cca0(i,k) = MAX(0.0, cca0(i,k))
    cca(i,k)  = MAX(0.0, cca(i,k))

    IF (cca0(i,k) >= 1.0) THEN
      cca0(i,k) = MIN(0.99e+0, cca0(i,k))
    END IF
    IF (cca(i,k) >= 1.0) THEN
      cca(i,k) = MIN(0.99e+0, cca(i,k))
    END IF

  END DO
END DO


!-------------------------------------------------------

IF (lhook) CALL dr_hook('GLUE_CONV_4A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE glue_conv_4a

    END MODULE glue_conv_4a_mod
