! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE conv_diag_4a_mod

  USE UM_ParParams
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

CONTAINS

!
! To diagnose convective occurrence and type
!
      SUBROUTINE conv_diag_4a(                                         &

! IN values defining field dimensions and subset to be processed :
      row_length, rows                                                 &

! IN values defining vertical grid of model atmosphere :
     , bl_levels, model_levels, wet_model_levels                       &
     , land_points                                                     &
     , p, P_theta_lev, exner_rho                                       &
     , rho_only, rho_theta, z_full, z_half                             &

! IN Model switches
     , l_mixing_ratio, l_ctile, l_extra_call                           &
     , no_cumulus                                                      &

! IN Cloud data :
     , qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :
     , pstar, q, theta, exner_theta_levels, u_p, v_p, u_0_p, v_0_p     &
     , tstar_land, tstar_sea, tstar_sice, z0msea                       &
     , L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm          &
     , tstar, land_mask, flandg, ice_fract, timestep                   &
     , w_copy, w_max                                                   &
     , deep_flag, past_precip, past_conv_ht                            &

! SCM Diagnostics (dummy values in full UM)
     , nSCMDpkgs, L_SCMDiags                                           &

! OUT data required elsewhere in UM system :
     , zh,zhpar,dzh,z_lcl,z_lcl_uv,delthvu,ql_ad,ntml,ntpar,NLCL       &
     , cumulus, L_shallow,l_congestus, l_congestus2, conv_type, cin    &
     , CAPE, wstar, wthvs, entrain_coef, qsat_lcl                      &
     , ERROR)


! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                       &
  pdims, pdims_s, tdims_s, tdims, qdims, wdims

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m) 
 ,r_rho_levels                  ! Radii on rho levels (m)

! Trig information
USE trignometric_mod,  ONLY:                                           &
  sin_theta_longitude, cos_theta_longitude

USE cv_derived_constants_mod, ONLY:                                    &
  ls, lsrcp, lcrcp, gamma_dry, ra2 

USE cv_diag_param_mod, ONLY:                                           &
  a_plume, b_plume, max_t_grad, sc_cftol,                              &
  a_bolton, b_bolton, c_bolton, d_bolton

USE cv_param_mod, ONLY:                                                &
  max_diag_thpert

USE cv_run_mod, ONLY:                                                  &
  icvdiag, tv1_sd_opt, limit_pert_opt

    USE surf_param, ONLY : z0sice, z0h_z0m_sice
    USE c_rough,    ONLY : z0hsea
    USE water_constants_mod, ONLY: lc, lf, tm

USE earth_constants_mod, ONLY: g, earth_radius

USE atmos_constants_mod, ONLY:                                         &
   vkman, cp, kappa, r, repsilon, c_virtual, recip_kappa

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   To diagnose convective occurrence and type
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine Arguments
!
! Arguments with intent IN:
!
! (a) Defining horizontal grid and subset thereof to be processed.

INTEGER, INTENT(IN) ::  &
  row_length            & ! Local number of points on a row
 ,rows                    ! Local number of rows in a theta field

! (b) Defining vertical grid of model atmosphere.

INTEGER, INTENT(IN) :: &
  bl_levels            & ! Max. no. of "boundary" levels allowed.
 ,model_levels         & ! number of model levels
 ,wet_model_levels     & ! number of wet model levels
 ,land_points            ! number of land points

REAL, INTENT(IN) ::                          &
  p(pdims_s%i_start:pdims_s%i_end,           & ! pressure  on rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,           &
    pdims_s%k_start:pdims_s%k_end)           &
 ,p_theta_lev(tdims%i_start:tdims%i_end,     & ! P on theta lev (Pa)
              tdims%j_start:tdims%j_end,     &
                          1:tdims%k_end)     &
 ,exner_rho(pdims_s%i_start:pdims_s%i_end,   & ! Exner on rho level
            pdims_s%j_start:pdims_s%j_end,   & !
            pdims_s%k_start:pdims_s%k_end)   &
 ,rho_only(row_length,rows,1:tdims%k_end)    & ! density (kg/m3)
 ,rho_theta(row_length,rows,1:tdims%k_end-1) & ! rho th lev (kg/m3)
 ,z_full(row_length,rows,1:tdims%k_end)      & ! height th lev (m)
 ,z_half(row_length,rows,1:tdims%k_end)        ! height rho lev (m)


LOGICAL,INTENT(IN) ::   &
  l_mixing_ratio        & ! true moisture input as mixing ratios
                          ! false moisture input as specific humidity
 ,l_ctile               & ! true if coastal tiling
 ,l_flux_bc             & ! true if SCM using specified surface fluxes
 ,l_spec_z0             & ! true if roughness length has been specified
 ,l_extra_call            ! true this is an additional call to conv_diag
                          ! within a timestep

LOGICAL,INTENT(IN) ::   &
  no_cumulus(row_length,rows)   ! Points overruled by BL     

! (c) Cloud data.

REAL, INTENT(IN) ::                          &
  qcf(qdims%i_start:qdims%i_end,             & ! Cloud ice (kg per kg air)
      qdims%j_start:qdims%j_end,             &
                  1:qdims%k_end)             & 
 ,qcl(qdims%i_start:qdims%i_end,             & !Cloud liquid water (kg/kg air)
      qdims%j_start:qdims%j_end,             &
                  1:qdims%k_end)             & 
 ,cloud_fraction(qdims%i_start:qdims%i_end,  & !  Cloud fraction
                 qdims%j_start:qdims%j_end,  &
                             1:qdims%k_end)

! (d) Atmospheric + any other data not covered so far, incl control.

REAL, INTENT(IN) ::                             &
  pstar(row_length, rows)                       & ! Surface pressure (Pa)
 ,q(qdims%i_start:qdims%i_end,                  & ! water vapour (kg/kg) 
    qdims%j_start:qdims%j_end,                  &
                1:qdims%k_end)                  & 
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,              &
                    1:tdims%k_end)              &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)  
REAL, INTENT(IN) ::            &
  u_p(row_length, rows)        & ! U(1) on P-grid.
 ,v_p(row_length, rows)        & ! V(1) on P-grid.
 ,u_0_p(row_length,rows)       & ! W'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,v_0_p(row_length,rows)       & ! S'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,flux_e(row_length,rows)      & ! Specified surface
                                 !    latent heat flux (W/m^2)
 ,flux_h(row_length,rows)      & ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
 ,z0msea(row_length,rows)      & ! Sea roughness length for momentum (m)
 ,z0m_scm(row_length,rows)     & ! Namelist input z0m (if >0)
 ,z0h_scm(row_length,rows)       ! Namelist input z0h (if >0)

REAL, INTENT(IN) :: &
  tstar_land(row_length, rows) & ! Surface T on land
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows)   ! Surface T on sea-ice 

REAL, INTENT(INOUT) ::         &
  tstar(row_length,rows)         ! Surface temperature (K)
                                 ! NOTE only inout for SCM

LOGICAL,INTENT(IN) ::          &
  land_mask(row_length, rows)   ! T If land, F Elsewhere.

REAL, INTENT(IN) ::            &
  flandg(pdims_s%i_start:pdims_s%i_end, & ! Land fraction of gridbox
         pdims_s%j_start:pdims_s%j_end) & ! on all points
 ,ice_fract(row_length,rows)              ! fraction of sea that has ice

REAL, INTENT(IN) ::                  &
  timestep                           & ! timestep (seconds).
 ,w_copy(wdims%i_start:wdims%i_end,  & ! vertical velocity W (m/s)
         wdims%j_start:wdims%j_end,  &  
                     0:wdims%k_end)  & ! Not exact match to module values
 ,w_max(row_length,rows)               ! Column max vertical velocity (m/s)

REAL, INTENT(IN) ::              &
  deep_flag(row_length,rows)     & ! 0-1.0, 1 if deep last time step
 ,past_precip(row_length,rows)   & ! convective precip rate last step
                                   ! or a decayed value.
 ,past_conv_ht(row_length,rows)    ! convective height (m)

     
! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::  &
  nSCMDpkgs               ! No of diagnostics packages

LOGICAL,INTENT(IN) ::   &
  L_SCMDiags(nSCMDpkgs)   ! Logicals for diagnostics packages

REAL, INTENT(INOUT) ::           & 
  zh(row_length,rows)              ! Height above surface of top
                                   !  of boundary layer (metres).

REAL               ::            & ! Not used (OUT shouldnt be used)
 qsat_lcl(row_length,rows)       & ! qsat at cloud base (kg/kg) 
                                   ! (not used in routine)
 ,entrain_coef(row_length,rows)    ! Entrainment coefficient
                                   ! (not used in routine)

REAL, INTENT(OUT) ::             &
  zhpar(row_length,rows)         & ! Height of max parcel ascent (m)
 ,dzh(row_length,rows)           & ! Dummy here (used from 5A)
 ,z_lcl(row_length,rows)         & ! Height of lifting condensation 
                                   ! level  (m)
 ,z_lcl_uv(row_length,rows)      & ! Height of lifting condensation
                                   ! level on uv grid (m)
 ,delthvu(row_length,rows)       & ! Integral of undilute parcel buoyancy
                                   ! over convective cloud layer
                                   ! (for convection scheme)
 ,ql_ad(row_length,rows)         & ! adiabatic liquid water content at 
                                   ! inversion or cloud top (kg/kg)
 ,cape(row_length, rows)         & ! CAPE from parcel ascent (m2/s2)
 ,cin(row_length, rows)            ! CIN from parcel ascent (m2/s2)

INTEGER, INTENT(OUT) ::      &
  ntml(row_length,rows)      & ! Number of model levels in the
                               ! turbulently mixed layer.
 ,ntpar(row_length,rows)     & ! Max levels for parcel ascent
 ,nlcl(row_length,rows)        ! No. of model layers below the
                               ! lifting condensation level.

LOGICAL,INTENT(OUT) ::         &
  cumulus(row_length,rows)     & ! Logical indicator for convection
 ,l_shallow(row_length,rows)   & ! Logical indicator for shallow Cu
 ,l_congestus(row_length,rows) & ! Logical indicator for congestus Cu
 ,l_congestus2(row_length,rows)  ! Logical ind 2 for congestus Cu

! Required for extra call to conv_diag as values from original BL call
! may return zero values based on initial timestep profiles.
REAL, INTENT(INOUT) ::    &
  wstar(row_length,rows)  & ! Convective sub-cloud velocity scale (m/s)
 ,wthvs(row_length,rows)    ! surface flux of wthv (K m/s)
INTEGER, INTENT(INOUT) :: &
  error                     ! 0 - no error in this routine

!
! Redundant arguments
! -------------------
! Convective type array ::
INTEGER, INTENT(IN) :: &
  conv_type(row_length, rows)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'conv_diag_4a'
CHARACTER (len=80) ::  cmessage       ! error message

INTEGER ::  &
  i,j       & ! LOCAL Loop counter (horizontal field index).
 ,ii        & ! Local compressed array counter.
 ,k         & ! LOCAL Loop counter (vertical level index).
 ,mbl       & ! Maximum number of model layers allowed in the
!             ! mixing layer; set to bl_levels-1.
 ,nunstable   ! total number of unstable points

! uncompressed arrays - all points

REAL ::                         &
  qs_star(row_length, rows)     & ! Saturated sp humidity at surface
 ,qs_star_sice(row_length, rows)& ! Saturated sp humidity at sea-ice surface
 ,fb_surf(row_length, rows)     & ! Change in theta_v from surface
                                  ! to layer 1 (note diff from BL)
 ,dqsdt(row_length,rows)        & ! d(QSAT)/dT 
 ,z0(row_length,rows)           & ! roughness length (m)
 ,z0m_land(row_length,rows)     & ! roughness length for momentum over land (m)
 ,z0h_land(row_length,rows)     & ! roughness length for heat over land (m)
 ,z0m_sea(row_length,rows)      & ! roughness length for momentum over sea(m)
 ,z0h_sea(row_length,rows)        ! roughness length for heat over sea(m)

! Used in calculation to decide on unstable points 

REAL ::           &
  theta1          &  ! Potential temperature in layer 1
 ,ushear          &  ! U wind shear from level 1 to surface
 ,vshear          &  ! V wind shear from level 1 to surface
 ,wshr1           &  ! magnitude of surface wind shear
 ,wshr2           &  ! (wshr1)**2
 ,rhostar         &  ! surface air density
 ,theta_star      &  ! theta at surface
 ,wthvbar         &  ! surface buoyancy flux
 ,cd              &  ! bulk transfer coefficient for momentum
 ,ch              &  ! bulk transfer coefficient for heat
 ,Rib             &  ! Ri for surface exchange
 ,ustar           &  ! surface friction velocity
 ,ustar_eff       &  ! effective surface friction velocity
 ,w_m             &  ! 
 ,w_s_cubed       &  !
 ,wstar_tmp       &  ! convective velocity scale
 ,ustar2_sea      &  ! ustar^2 over sea
 ,ustar2_land     &  ! ustar^2 over land
 ,wthvbar_sea     &  ! surface buoyancy flux over sea
 ,wthvbar_land       ! surface buoyancy flux

! compressed arrays store only values for unstable columns

INTEGER ::                   &
  index_i(row_length*rows)   & ! column number of unstable points
 ,index_j(row_length*rows)     ! row number of unstable points

! compressed copys for unstable points - names as original array plus _c

INTEGER ::                   &
  ntml_c(row_length*rows)    &
 ,nlcl_c(row_length*rows) 

REAL ::                                                            &
  q_c(row_length*rows, 1:tdims%k_end)                              &
 ,qcl_c(row_length*rows, 1:tdims%k_end)                            &
 ,qcf_c(row_length*rows, 1:tdims%k_end)                            &
 ,z_full_c(row_length*rows, 1:tdims%k_end)                         &
 ,z_half_c(row_length*rows, 1:tdims%k_end)                         &
 ,exner_theta_levels_c(row_length*rows, 1:tdims%k_end)             &
 ,P_theta_lev_c(row_length*rows, 1:tdims%k_end)                    &
 ,r_theta_levels_c( row_length *rows, 1:tdims%k_end)               &
 ,z_lcl_c(row_length*rows)                                         &
 ,zh_c(row_length*rows)                                            &
 ,delthvu_c(row_length*rows)                                       &
 ,cape_c(row_length*rows)                                          &
 ,cin_c(row_length*rows) 

! Arrays only used for unstable calculations

LOGICAL ::                   &
  topbl(row_length*rows)     & ! Flag set when top of boundary layer
                               ! is reached.
 ,topprof(row_length*rows)   & ! Flag set when top of ascent
                               ! is reached.
 ,above_lcl(row_length*rows) & ! Flag set when parcel above LCL.
 ,shmin(row_length*rows)       ! Flag for finding min in parcel
                               ! buoyancy below 3km (for shallow Cu)

REAL ::                                   &
  t(row_length*rows, 1:tdims%k_end)       & ! temperature (from theta)
 ,tl(row_length*rows, 1:tdims%k_end)      & ! Ice/liquid water temperature,
                                            ! but replaced by T in LS_CLD.
 ,qw(row_length*rows, 1:qdims%k_end)      & ! Total water content
 ,svl(row_length*rows, 1:qdims%k_end)     & ! Liquid/frozen water virtual
                                            ! static energy over CP.
 ,env_svl(row_length*rows, 1:qdims%k_end) & ! Density (virtual) static energy
                                            ! over CP for layer.
 ,par_svl(row_length*rows, 1:qdims%k_end) & ! Density (virtual) static energy
                                            ! over CP of parcel for level.
 ,qs(row_length*rows, 1:qdims%k_end)        ! Saturated sp humidity at pressure
                                            ! and temperature of sucessive 
                                            ! levels.
REAL ::                            &
  t_lcl(row_length*rows)           & ! Temperature at lifting condensation level
 ,p_lcl(row_length*rows)           & ! Pressure at lifting condensation level.
 ,sl_plume(row_length*rows)        & ! Liquid/frozen water static energy
                                     ! over CP for a plume rising without
                                     ! dilution from level 1.
 ,qw_plume(row_length*rows)        & ! QW for a plume rising without
                                     ! dilution from level 1.
 ,th_ref(row_length*rows)          & ! theta - reference value
 ,t_ref(row_length*rows)           & ! t - reference value
 ,qsat_lev(row_length*rows)        & ! qsat for reference temperature
 ,qsat_env(row_length*rows)        & ! qsat for environment temperature
 ,dt_dens_parc_t(row_length*rows)  & ! t_dens_parc-t_dens_env at ntpar
 ,dt_dens_parc_tmin(row_length*rows) & ! t_dens_parc-t_dens_env at kshmin
 ,thv_pert(row_length*rows)        & ! threshold thv of parcel
 ,dtv_min(row_length*rows)         & ! min TV of parcel in cld layer
 ,tv1_sd(row_length*rows)          & ! Approx to standard dev of level
                                     ! 1 virtual temperature (K).
 ,env_svl_km1(row_length*rows)     & ! Density (virtual) static energy
                                     ! over CP for last layer considered.
 ,par_svl_km1(row_length*rows)     & ! Density (virtual) static energy over
                                     ! CP of parcel at last level considered.
! Added for improved parcel top - mainly used for finding an ascent
! capped by an inversion.
 ,dt_dens_parc_t2(row_length*rows) & ! 2nd copy of Dt_dens_parc_T
 ,dtv_min2(row_length*rows)        & ! 2nd copy min TV of parcel
 ,delthvu2(row_length*rows)        & ! 2nd copy
 ,zh2(row_length*rows)             & ! 2nd copy
 ,max_buoy(row_length*rows)          ! max parcel buoyancy 

! arrays holding various key model levels

INTEGER ::                       &
  kshmin(row_length*rows)        & ! Position of buoyancy minimum above
                                   ! topbl (used for shallow Cu diag)
 ,kcucheck(row_length*rows)      & ! Position of level just below 2.5km
                                   ! (used for gradient check to
                                   !  diagnose convection)
 ,k_plume(row_length*rows)       & ! start level for surface-driven plume
 ,k_max(row_length*rows)         & ! level of max parcel buoyancy
 ,nlcl_min(row_length*rows)      & ! min level of LCL
 ,k_neutral(row_length*rows)     & ! level of neutral parcel buoyancy
 ,k_inv(row_length*rows)         & ! level from inversion testing
 ,freeze_lev(row_length*rows)      ! freezing level

! parcel calculation

REAL ::                                       &
  t_parc(row_length*rows, 1:qdims%k_end)      & ! Temperature of parcel.
 ,t_dens_parc(row_length*rows, 1:qdims%k_end) & ! Density potential temperature
                                                ! of parcel.
 ,t_dens_env(row_length*rows, 1:qdims%k_end)  & ! Density potential temperature 
                                                ! of environment.
 ,denv_bydz(row_length*rows, 1:qdims%k_end)   & ! Gradient of density potential
                                                ! temperature in the environment
 ,dpar_bydz(row_length*rows, 1:qdims%k_end)   & ! Gradient of density potential
                                                ! temperature of the parcel.
 ,buoyancy(row_length*rows, 1:qdims%k_end)    & ! undilute parcel buoyancy (K)
 ,buoyancy_dil(row_length*rows, 1:qdims%k_end)& ! dilute parcel buoyancy (K)
 ,th_par_km1(row_length*rows)                   ! parcel theta at level below

REAL ::           &
  z_surf          &  ! approx height of top of surface layer
 ,vap_press       &  ! Vapour pressure.
 ,grad_cld        &  ! SVL gradient in layer above LCL.
 ,grad_sub        &  ! SVL gradient in layer below LCL.
 ,q_liq_env       &  ! Condensed water content of environment.
 ,dq_sat_env      &  ! DQSAT/DT for environment
 ,lrcp_const      &  ! lc or lc+lf over cp
 ,lrcp_const_env  &  ! lc or lc+lf over cp
 ,lrcp_const_parc &  ! lc or lc+lf over cp
 ,l_const         &  ! lc or lc+lf
 ,l_const_env     &  ! lc or lc+lf
 ,grcp            &  ! g/cp
 ,dz              &  ! layer depth
 ,z_pr            &  ! used in estimating th_ref at next level
 ,th_par          &  ! theta value for parcel
 ,inc             &  ! CIN/CAPE increment for layer
 ,dq_sat_par      &  ! dqsat/dT for parcel
 ,dq_sat_par_dil  &  ! dqsat/dT for dilute parcel
 ,temp_parc       &  ! average temperature of parcel after entrainment
 ,q_vap_parc      &  ! Vapour content of undilute parcel
 ,q_liq_parc      &  ! Liquid water content of undilute parcel
 ,qcl_parc        &  ! parcel qcl - dilute parcel cal
 ,qcf_parc        &  ! parcel qcf - dilute parcel cal 
 ,factor          &  ! multiplying factor 
 ,denv_bydz2      &  ! Gradient of density potential
                     ! temperature in the environment.
 ,dpar_bydz2         ! Gradient of density potential
                     ! temperature of the parcel.
     
! required for average w calculation

REAL  ::                                     &
  dmass_theta(row_length*rows,1:tdims%k_end)   ! r**2rho*dr on theta levels

REAL  ::                   &
  w_avg(row_length*rows)   & ! mean w over layer (m/s)
 ,w_avg2(row_length*rows)  & ! mean w over layer (m/s)
 ,mass(row_length*rows)      ! mass for column

! Arrays added for dilute parcel calculation
      
REAL  ::                                         &
  entrain_fraction(row_length*rows,1:tdims%k_end)& ! fraction of environmental
                                                   ! air to mix with parcel
 ,qw_parc(row_length*rows,1:tdims%k_end)         & ! parcel total water
 ,sl_parc(row_length*rows,1:tdims%k_end)         & ! parcel SL
 ,ql_parc(row_length*rows,1:tdims%k_end)         & ! parcel water
 ,t_parc_dil(row_length*rows,1:tdims%k_end)      & ! dilute parcel temeperature
 ,t_dens_parc_dil(row_length*rows,1:tdims%k_end) & ! dilute parcel t_dens
 ,ql_parc_dil(row_length*rows,1:tdims%k_end)       ! dilute parcel liquid water

REAL  ::                         &
  th_ref_dil(row_length*rows)    & ! dilute parcel ref potential temperature
 ,th_par_km_dil(row_length*rows) & ! dilute parcel ref potential temperature 2nd
 ,t_ref_dil(row_length*rows)     & ! dilute parcel reference temperature 
 ,qsat_lev_dil(row_length*rows)    ! qsat for dilute parcel

REAL ::                          &
     tv1_sd_temp(row_length, rows) ! temp store to allow parallel omp
! Arrays added for extra conv_diag calls
INTEGER ::                     &
  ntml_copy(row_length,rows)


LOGICAL ::     &
  l_keep_water & ! if true keeps water loading in plume
                 ! false removed if water exceeds 1g/kg
 ,l_wtest        ! do w test

!
! Model constants:
!



      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
! 
! svl = (Tl+gz/cp)*(1+(1/repsilon-1)qt)   
!     ~ (Tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T) 
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation 
!  e = rp/(repsilon+r)   - accurate    
!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CONV_DIAG_4A',zhook_in,zhook_handle)
      ERROR = 0

!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

      IF ( bl_levels <  1 .OR. rows <  1 .OR. tdims%k_end <  1 .OR.    &
            qdims%k_end <  1 ) THEN
        ERROR = 1
        goto 9999

      END IF

!-----------------------------------------------------------------------
! 1.1a copy ntml values where do not want to overwrite
!-----------------------------------------------------------------------
      If (l_extra_call) THEN
        Do j=1, rows
          Do i=1,row_length
            If (no_cumulus(i,j)) THEN
              ntml_copy(i,j) = ntml(i,j)
            END IF
          END DO
        END DO
      END IF

!-----------------------------------------------------------------------
! 1.2 initialisation of output arrays
!-----------------------------------------------------------------------
      Do j=1,rows
        Do i=1,row_length

          cumulus(i,j)     = .false.
          L_shallow(i,j)   = .false.
          L_congestus(i,j) = .false.
          L_congestus2(i,j) = .false.
          ntml(i,j) = 1
          nlcl(i,j) = 1
          ntpar(i,j) = 1
          delthvu(i,j)  = 0.0
          ql_ad(i,j)    = 0.0
          CAPE(i,j)     = 0.0
          CIN(i,j)      = 0.0
          zhpar(i,j)    = 0.0
          dzh(i,j)      =-9.9E9     ! initialise to large and negative 
          z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)

! Set LCL for UV grid to level below z_lcl (note that (at vn5.1) UV
! levels are half levels (below) relative to P-grid. Consistent with
! BL scheme's treatment of staggered vertical grid.)

          z_lcl_uv(i,j)= z_full(i,j,nlcl(i,j))

        END DO
      END DO
!

!-----------------------------------------------------------------------
! set variables
!-----------------------------------------------------------------------
!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL  = bl_levels - 1
!-----------------------------------------------------------------------
! 1.4 Set appropriate roughness length
!-----------------------------------------------------------------------
      IF (tv1_sd_opt == 2) THEN


!       ! using grid-box mean in coastal points, hence need land and sea 
!       ! separately
        DO J=1, rows
        DO I=1, row_length
          ! Land: assume z0m=0.1m and z0h=z0m/10.
          z0m_land(i,j) = 0.1
          z0h_land(i,j) = 0.01
          ! Sea: use parametrized values
          !      z0h_sea updated later for low wind speed limit
          z0m_sea(i,j) = z0msea(i,j)
          z0h_sea(i,j) = MAX( 2.56e-9/z0msea(i,j), 7.0e-08 )
        END DO ! I
        END DO ! J
        IF ( L_spec_z0 ) THEN
          ! Code to use z0mh_scm if namelist specifies it
          DO J=1, rows
          DO I=1, row_length
            If ( z0m_SCM(i,j)  >   0.0 ) THEN
              z0m_sea(i,j)  = z0m_SCM(i,j)
              z0m_land(i,j) = z0m_SCM(i,j)
            END IF
            If ( z0h_SCM(i,j)  >   0.0 ) THEN
              z0h_sea(i,j)  = z0h_SCM(i,j)
              z0h_land(i,j) = z0h_SCM(i,j)
            END IF
          END DO ! I
          END DO ! J
        END IF
      END IF

      DO J=1, rows
      DO I=1, row_length
        If (land_mask(I,j)) THEN
!       ! Approximate z0 for land as 0.1.
          z0(i,j) = 0.1
        Else
          z0(i,j) = z0hsea
        END IF
      END DO ! I
      END DO ! J
 
      If ( L_spec_z0 ) THEN

        ! Code to use z0h_scm if Namelist specifies it
        DO J=1, rows
        DO I=1, row_length
          If ( z0h_SCM(i,j)  >   0.0 ) THEN
            z0(i,j) = z0h_SCM(i,j)
          END IF ! z0h_scm
        END DO ! I
        END DO ! J
      END IF

!-----------------------------------------------------------------------
! 1.5 Surface buoyancy flux and calculation of unstable points
!-----------------------------------------------------------------------

      IF ( .NOT. L_flux_bc) THEN        ! used by most UM runs

       IF (tv1_sd_opt == 2) THEN
        
!-----------------------------------------------------------------------
! Calculate the surface buoyancy flux
! new method includes stability dependence and area mean for coastal 
! and sea-ice points
!-----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
        CALL qsat_mix(qs_star,tstar_sea,pstar,row_length*rows,        &
                     l_mixing_ratio)
! DEPENDS ON: qsat_mix
        CALL qsat_mix(qs_star_sice,tstar_sice,pstar,row_length*rows,  &
                     l_mixing_ratio)
        grcp = g/cp

        ii = 0
        DO j=1,rows
        DO i=1,row_length

         theta1 = theta(i,j,1)
         rhostar = pstar(i,j) / ( R*tstar(i,j) )

         Ushear = U_P(i,j) - U_0_P(i,j)
         Vshear = V_P(i,j) - V_0_P(i,j)
         wshr2 = MAX (1.0E-6 , Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)
!        !-------------------------------------------------------------
!        ! Sea 
!        !-------------------------------------------------------------
         wthvbar_sea = 0.0
         ustar2_sea   = 0.0
         IF ( flandg(i,j) < 0.99 ) THEN
!          ! Include a crude stability dependence
           rib = - ( g / tstar_sea(i,j) ) * &
            (tstar_sea(i,j) - ( theta1*exner_theta_levels(i,j,1)       &
                    +grcp*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
           cd = ( vkman / LOG(z_full(i,j,1)/z0m_sea(i,j)) )**2
           cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
           ustar2_sea = cd * wshr2
           IF ( .NOT.L_spec_z0 ) THEN
             ! include low wind speed limit (larger z0)
             z0h_sea(i,j) = MAX( 2.52e-6/(SQRT(ustar2_sea)+1.0e-05),   &
                                 z0h_sea(i,j) ) 
           END IF
           ch = vkman**2 / ( LOG(z_full(i,j,1)/z0m_sea(i,j)) *         &
                             LOG(z_full(i,j,1)/z0h_sea(i,j)))
           ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
           wthvbar_sea = ch * wshr1 * ( tstar_sea(i,j) -               &
             ( theta1*exner_theta_levels(i,j,1) + grcp*z_full(i,j,1) ) &
               + 0.61*theta1*(qs_star(i,j)-q(i,j,1)) )
!          !-------------------------------------------------------------
!          ! Sea-ice
!          !-------------------------------------------------------------
           IF ( ice_fract(i,j) > 0.01 ) THEN
!            ! Include a crude stability dependence
             rib = - ( g / tstar_sice(i,j) ) * &
              (tstar_sice(i,j) - ( theta1 * exner_theta_levels(i,j,1)  &
                    +grcp*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
             cd = ( vkman / LOG(z_full(i,j,1)/z0sice) )**2
             cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
             ustar2_sea = (1.0-ice_fract(i,j)) * ustar2_sea +          &
                               ice_fract(i,j) * cd*wshr2
             ch = vkman**2 / ( LOG(z_full(i,j,1)/z0sice) *             &
                               LOG(z_full(i,j,1)/                      &
                                  (z0sice*z0h_z0m_sice)) )
             ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
             wthvbar_sea = (1.0-ice_fract(i,j)) * wthvbar_sea +        &
                         ice_fract(i,j) * ch*wshr1*( tstar_sice(i,j) - &
              ( theta1*exner_theta_levels(i,j,1) + grcp*z_full(i,j,1) )&
                + 0.61*theta1*(qs_star_sice(i,j)-q(i,j,1)) )
           END IF
         END IF
!        !-------------------------------------------------------------
!        ! Land
!        !-------------------------------------------------------------
         wthvbar_land = 0.0
         ustar2_land   = 0.0
         IF ( flandg(i,j) > 0.01 ) THEN
!          Include a crude stability dependence
           rib = - ( g / tstar_land(i,j) ) * &
            (tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1)  &
             +grcp*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
           cd = ( vkman / LOG(z_full(i,j,1)/z0m_land(i,j)) ) ** 2
           cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
           ustar2_land = cd * wshr2
           ch = vkman**2 / ( LOG(z_full(i,j,1)/z0m_land(i,j)) *      &
                             LOG(z_full(i,j,1)/z0h_land(i,j) ) )
           ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
           wthvbar_land = ch * wshr1 * &
            ( tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1) &
             +grcp*z_full(i,j,1) ) )
         END IF
!        !-------------------------------------------------------------
!        ! Combine to cell average values and then take sqrt for ustar
!        !-------------------------------------------------------------
          IF ( flandg(i,j) < 0.01 ) THEN
            ustar = ustar2_sea
            wthvbar = wthvbar_sea
          ELSE IF ( flandg(i,j) > 0.99 ) THEN
            ustar = ustar2_land
            wthvbar = wthvbar_land
          ELSE
!           ! Take area-weighted mean
            ustar = (1.0 - flandg(i,j) ) * ustar2_sea +                &
                     flandg(i,j) * ustar2_land 
            wthvbar = (1.0 - flandg(i,j) ) * wthvbar_sea +             &
                     flandg(i,j) * wthvbar_land
          END IF
!
          ustar = SQRT(ustar)
          fb_surf(i,j) = g * wthvbar /                                 &
             ( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

          IF (fb_surf(i,j)  >   0.0) THEN
            ii= ii+1
            w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
            w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
            tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
          END IF ! (fb_surf > 0.0)

        END DO ! I
        END DO ! J

       ELSE  ! tv1_sd_opt /= 2
! old (neutral stability) method

! DEPENDS ON: qsat_mix
       Call qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, theta1, rhostar,           &
!$OMP& Ushear, Vshear, wshr2, wshr1, cd, wthvbar, ustar,                &
!$OMP& w_s_cubed, w_m)

!$OMP SINGLE
       ii = 0
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
        DO j=1,rows
        DO i=1,row_length

!-----------------------------------------------------------------------
! Calculate the surface buoyancy flux
! Uses approximation for unstable CD as 1.5*neutral value (defined
! Garratt p54) and that CH=CD. Approximate z0 for land as 0.1.
!-----------------------------------------------------------------------

         theta1 = theta(i,j,1)
         rhostar = pstar(i,j) / ( R*tstar(i,j) )

         Ushear = U_P(i,j) - U_0_P(i,j)
         Vshear = V_P(i,j) - V_0_P(i,j)
         wshr2 = MAX (1.0E-6 , Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)
         cd  = 1.5 * ( vkman/LOG(z_full(i,j,1)/z0(i,j)) )**2

         IF (land_mask(i,j)) THEN         ! land
            wthvbar = wshr1 * CD                                          &
             * ( TSTAR(I,j)*((100000.0/PSTAR(I,j))**kappa) - theta1 )
         ELSE
            wthvbar = wshr1 * CD                                          &
             * ( TSTAR(I,j)*((100000.0/PSTAR(I,j))**kappa) - theta1       &
               + 0.61*theta1*(QS_STAR(I,j)-Q(I,j,1)) )
         END IF

         ustar = SQRT(cd  * wshr2)
         fb_surf(i,j) = G * wthvbar /                                   &
                             ( rhostar * theta1*(1.0+0.61*Q(i,j,1)) )

         IF (fb_surf(i,j)  >   0.0) THEN

          IF (tv1_sd_opt == 0) THEN 
! old method assumes the BL depth ZH=300.            
           
           tv1_sd_temp(i,j) = 1.93 * wthvbar / ( rhostar * ( 75.0 *          &
                         fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
       
          ELSE IF (tv1_sd_opt == 1) THEN
! improved method uses BL depth from the previous timestep
           w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
           w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
           tv1_sd_temp(i,j) =  1.93 * wthvbar/( rhostar * w_m )   

          END IF
         END IF ! (fb_surf > 0.0)

        END DO ! I
        END DO ! J
!$OMP END DO

!$OMP SINGLE
        ii = 0
        DO j=1,rows
        DO i=1,row_length
         IF (fb_surf(i,j)  >   0.0) THEN
           ii = ii + 1
           tv1_sd(ii) = tv1_sd_temp(i,j)
         END IF
        END DO ! I
        END DO ! J
!$OMP END SINGLE

!$OMP END PARALLEL
      END IF


      ELSE ! if L_flux_bc       (used by some SCM runs)

!       !------------------------------------------
!       ! If specified surface fluxes are required
!       !------------------------------------------
        Do J=1, rows
        Do I=1, row_length
          ! For taylor expansion about T0=SL(K=1)
          tstar(I,J) = theta(i,j,1) * exner_theta_levels(i,j,1)         &
                                             +gamma_dry*z_full(I,J,1)
        END DO
        END DO

! DEPENDS ON: qsat_mix
        Call qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

! dqsat/dT - approximation same expression whether specific humidity or 
!            mixing ratio.

        DO J=1, rows
        DO I=1, row_length
          dqsdt(I,J) = (repsilon * LC * qs_star(I,J))                    &
                       / ( R * tstar(I,J) * tstar(I,J) )
        END DO
        END DO

       If (icvdiag <= 1 .OR. icvdiag >= 5 ) THEN

        ii = 0

        Do j=1,rows
        Do i=1,row_length

         Ushear = U_P(I,j) - U_0_P(I,j)
         Vshear = V_P(I,j) - V_0_P(I,j)
!        ! Need to have a higher minimum wind speed limit with
!        ! specified fluxes in order not to generate huge TSTAR
         wshr2 = MAX (0.1, Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)

         ! Calculate WTHV from namelist flux_h and flux_e (in W/m2)
         wthvbar = ((flux_h(I,j)/CP)+0.61*(flux_e(I,j)/LC))             &
                 * ((100000.0/PSTAR(I,j))**kappa)

         cd  = 1.5 * ( vkman/LOG(z_full(I,j,1)/z0(i,j)) )**2

         theta1 = theta(I,j,1)

         ! Taylor expansion for qsat(T*) about SL(k=1)
         TSTAR(I,j) = ( theta1 + (wthvbar/(wshr1*CD))                   &
                    -   0.61*theta1                                     &
                    *   (QS_STAR(I,j)-Q(I,j,1)-DQSDT(I,j)*TSTAR(I,j)))  &
                    /   ( (100000.0/PSTAR(I,j))**kappa +                &
                             0.61*theta1*DQSDT(I,j) )

         RHOSTAR = PSTAR(I,j) / ( R*TSTAR(I,j) )

         USTAR = SQRT(CD * wshr2)
         FB_SURF(I,j) = G * wthvbar /                                   &
                             ( RHOSTAR * theta1*(1.0+0.61*Q(I,j,1)) )


         If (fb_surf(i,j)  >   0.0) then
          ii= ii+1
          If (tv1_sd_opt == 0) then 
! old method assumes the BL depth ZH=300.
           tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *           &
                         fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
 
          Else 
! improved method uses BL depth from the previous timestep
           w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
           w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
           tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
          End If
         End If ! (fb_surf > 0.0)

        END DO ! I
        END DO ! J
       
       Else    ! (icvdiag = 2-4)

        ii=0
        Do j=1,rows
        Do i=1,row_length

         Ushear = U_P(I,j) - U_0_P(I,j)
         Vshear = V_P(I,j) - V_0_P(I,j)
!        ! Need to have a higher minimum wind speed limit with
!        ! specified fluxes in order not to generate huge TSTAR
         wshr2 = MAX (0.1, Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)

         ! Calculate WTHV from namelist flux_h and flux_e (in W/m2)
         wthvbar = ((flux_h(I,j)/CP)+0.61*(flux_e(I,j)/LC))             &
                 * ((100000.0/PSTAR(I,j))**kappa)

         CD = 1.5 * ( VKMAN/log(z_full(I,j,1)/z0(i,j)) )**2

         theta1 = theta(I,j,1)
! new 
!         WSTAR_tmp = ( ZH(I,j) * G * 
!                   WTHVBAR / ( THETA1*(1.0+0.61*Q(I,j,1)) ) )**(1./3.)

         !  Use "effective" ustar, allowing for convective eddies 
!         USTAR_EFF = SQRT( CD*WSHR2 + BETA*BETA*WSTAR_tmp*WSTAR_tmp )

         ! Taylor expansion for qsat(T*) about SL(k=1)
!new         TSTAR(I,j) = ( theta1 + (wthvbar/(ustar_eff1*CD))          &
         TSTAR(I,j) = ( theta1 + (wthvbar/(wshr1*CD))                   &
                    -   0.61*theta1                                     &
                    *   (QS_STAR(I,j)-Q(I,j,1)-DQSDT(I,j)*TSTAR(I,j)))  &
                    /   ( (100000.0/PSTAR(I,j))**kappa +                &
                             0.61*theta1*DQSDT(I,j) )

         RHOSTAR = PSTAR(I,j) / ( R*TSTAR(I,j) )

         USTAR = SQRT(CD * wshr2)
         FB_SURF(I,j) = G * wthvbar /                                   &
                             ( RHOSTAR * theta1*(1.0+0.61*Q(I,j,1)) )

         If (fb_surf(i,j)  >   0.0) then
          ii= ii+1
          If (tv1_sd_opt == 0) then 
! old method assumes the BL depth ZH=300.
           tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *           &
                         fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
          Else 
! improved method uses BL depth from the previous timestep
           w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
           w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
           tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
          End If
         End If ! (fb_surf > 0.0)

        END DO ! I
        END DO ! J
  
       END IF  ! test on icvdiag

      END IF !  L_flux_bc

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, ii)
!-----------------------------------------------------------------------
! 2.0 Decide on unstable points  ( fb_surf > 0.0 )
!     Only work on these points for the rest of the calculations.
!-----------------------------------------------------------------------

! this confirms indirection indices as injective
! changing indexing algorithm may invalidate this
!$OMP MASTER
      nunstable = 0           ! total number of unstable points

      Do j=1,rows
      Do i=1,row_length
        If ( fb_surf(i,j)  >   0.0 ) THEN
          nunstable = nunstable + 1
          index_i(nunstable) = i
          index_j(nunstable) = j
        END IF 
      END DO
      END DO
!$OMP END MASTER
!$OMP BARRIER

!-----------------------------------------------------------------------
! 2.1 initialise just unstable points
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable

          shmin(ii)   = .false.
          topbl(ii)   = .false.
          topprof(ii) = .false.
          kshmin(ii)     = 1
          kcucheck(ii)   = 1
          k_max(ii)      = 1
          k_neutral(ii)  = 1
          k_inv(ii)      = 1
          freeze_lev(ii) = 1
          dtv_min(ii)   = 0.0
          dtv_min2(ii)  = 0.0
          delthvu2(ii)  = 0.0
          max_buoy(ii)  = 0.0
          ntml_c(ii)    = 1
          nlcl_c(ii)    = 1
          delthvu_c(ii) = 0.0
          cape_c (ii)   = 0.0
          cin_c  (ii)   = 0.0

      END DO
!$OMP END DO

! compress boundary layer depth
!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii)  = zh(i,j)
        z_lcl_c(ii) = z_half(i,j,nlcl_c(ii)+1)
      END DO
!$OMP END DO




!-----------------------------------------------------------------------
! 2.1 Calculate mass of layers 
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      DO k=1, tdims%k_end-1
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          dmass_theta(ii,k) = rho_theta(i,j,k)*                         &
                                  (r_theta_levels(i,j,k)/earth_radius)  &
                                 *(r_theta_levels(i,j,k)/earth_radius)  &
                    *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
        END DO
      END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 3.0 Calculate various quantities required by the parcel calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 3.1 Section BL.1 Calculate T at old time level.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, tdims%k_end
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          T(ii,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

            ! initialise t_parc at all points
            ! added for safety of qsat_mix calls later

          t_parc(ii,k) = t(ii,k)
          q_c(ii,k)   = q(i,J,K)
          qcl_c(ii,k) = qcl(i,J,K)
          qcf_c(ii,k) = qcf(i,J,K)
          z_full_c(ii,k) = z_full(i,j,k)  
          z_half_c(ii,k) = z_half(i,j,k)  
          exner_theta_levels_c(ii,k) = exner_theta_levels(i,j,k)
          P_theta_lev_c(ii,k) = P_theta_lev(i,j,k)
          r_theta_levels_c(ii,k) =r_theta_levels (i,j,k)
        END DO
      END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 3.2 Calculate total water content, QW and Liquid water temperature, TL
!     Definitions for Qw and Tl the same whether q, qcl, qcf  are 
!     specific humidities  or mixing ratio.
!     
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      DO k=1,qdims%k_end
        Do ii=1, nunstable

          QW(ii,K) = q_c(ii,K) + qcl_c(ii,k) + qcf_c(ii,k)
                                ! BL doc P243.10
          TL(ii,K) = T(ii,K) - LCRCP*qcl_c(ii,k) - LSRCP*qcf_c(ii,k)
                                ! BL doc P243.9

!
! Calculate SVL: conserved variable  a form of moist static energy /cp 
!       svl = (Tl+gz/cp)*(1+(1/repsilon-1)qt)  - specific humidity
!

          SVL(ii,K) = ( TL(ii,K) + gamma_dry * z_full_c(ii,K) )             &
                                     * ( 1.0 + c_virtual*QW(ii,K) )

        END DO
      END DO

!$OMP END DO

!$OMP END PARALLEL 



!-----------------------------------------------------------------------
! 4.0 Parcel calculation - various options available 
!-----------------------------------------------------------------------
!  icvdiag option         explaination
!     0              original undilute parcel calculation    
!     1              undilute parcel calculation (improved version).
!     2 and above    dilute parcel calculation
!        2           1/z entrainment rates ?
!        3            0.55/z
!=======================================================================
! Choice of diagnosis method
!=======================================================================
      IF (ICVDIAG == 0) THEN   ! old code Still required at present
!=======================================================================
!-----------------------------------------------------------------------
!! 0.1 Set TOPBL to .FALSE. and calculate boundary layer top using
!!     a non-local, plume method.
! ----------------------------------------------------------------------
!
      Do ii=1, nunstable

! Initialise CAPE and CIN
        CAPE_c(ii) = 0.0
        CIN_c(ii)  = 0.0

        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent if unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*ZH
!-----------------------------------------------------------------------
        z_surf = 0.1 * zh_c(ii)

        Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
                    SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

            k_plume(ii) = k_plume(ii) + 1

        END DO
      END DO       ! ii loop 

      If (limit_pert_opt == 2) Then
        Do ii=1, nunstable
          sl_plume(ii) = TL(ii,k_plume(ii))                              &
                              + gamma_dry * z_full_c(ii,k_plume(ii))
          thv_pert(ii) = min( max( a_plume,                              &
                          min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) &
                          ), max_diag_thpert )
          qw_plume(ii) = QW(ii,k_plume(ii))
        End Do      ! ii loop 
      Else If (limit_pert_opt == 0 .or. limit_pert_opt == 1) Then
        Do ii=1, nunstable
          sl_plume(ii) = TL(ii,k_plume(ii))                              &
                              + gamma_dry * z_full_c(ii,k_plume(ii))
          thv_pert(ii) = max( a_plume,                                   &
                          min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ))
          qw_plume(ii) = QW(ii,k_plume(ii))
        End Do      ! ii loop 
      End If
      
!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/repsilon   q specific humidity
!   vapour pressure e ~ qp/(repsilon+q)   q mixing ratio

      If (l_mixing_ratio) THEN
        Do ii=1, nunstable
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
                                       P_theta_lev_c(ii,k_plume(ii))   &
                            / (repsilon+Q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
           i = index_i(ii)   
           j = index_j(ii)   
           P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      Else       ! not l_mixing_ratio
! expression for specific humidity
        Do ii=1, nunstable
          vap_press = Q_c(ii,k_plume(ii)) *                            &
               P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*repsilon )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      END IF ! test on l_mixing_ratio  

!-----------------------------------------------------------------------
! calculate parcel water by linearising qsat about the environmental
! temperature
! Note: the following calculation is for parcel in level 1:
!
        k = 1 
! DEPENDS ON: qsat_mix
        Call QSAT_MIX(QS(1,K),T(1,K),P_theta_lev_c(1,K),               &
                                           nunstable,l_mixing_ratio)

        Do ii=1, nunstable

         IF(T(ii,K) >  TM) THEN
           DQ_SAT_ENV = repsilon*LC*QS(ii,K)/(R*T(ii,K)**2)

           Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -               &
             DQ_SAT_ENV*( SL_plume(ii)-gamma_dry*z_full_c(ii,K)-T(ii,K) )   &
                                   ) / (1.0+LCRCP*DQ_SAT_ENV) )

           Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*    &
              ( TL(ii,K) - T(ii,K) ) ) / (1.0+LCRCP*DQ_SAT_ENV) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
           Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)         &
                          - Q_LIQ_ENV
           T_PARC(ii,k)=                                               &
                  SL_plume(ii)-gamma_dry*z_full_c(ii,K)+LCRCP*Q_LIQ_PARC
         Else
           DQ_SAT_ENV=repsilon*LS*QS(ii,K)/(R*T(ii,K)**2)

           Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -          &
             DQ_SAT_ENV*( SL_plume(ii)-gamma_dry*z_full_c(ii,K)-T(ii,K) )   &
                                   ) / (1.0+LSRCP*DQ_SAT_ENV) )

           Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*    &
               ( TL(ii,K) - T(ii,K) ) ) / (1.0+LSRCP*DQ_SAT_ENV) )

! add on difference in environment's ql between RH_CRIT and RH_CRIT=1

           Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)         &
                           - Q_LIQ_ENV

           T_PARC(ii,k) =                                              &
                  SL_plume(ii)-gamma_dry*z_full_c(ii,K)+LSRCP*Q_LIQ_PARC
         END IF     ! test on TM

         Q_VAP_PARC=QW_plume(ii)-Q_LIQ_PARC

         T_DENS_PARC(ii,k) =                                           &
               T_PARC(ii,k) *(1.0+C_VIRTUAL*Q_VAP_PARC-Q_LIQ_PARC)
         T_DENS_ENV(ii,k)=T(ii,K)*                                     &
                      (1.0+C_VIRTUAL*Q_c(ii,K)-QCL_c(ii,K)-QCF_c(ii,K))

         ENV_SVL_KM1(ii) = T_DENS_ENV(ii,k)+gamma_dry*(z_full_c(ii,K))
         PAR_SVL_KM1(ii) = T_DENS_PARC(ii,k)+gamma_dry*(z_full_c(ii,K))
       END DO         ! ii loop


!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      do j=1,rows
        do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        END DO
      END DO
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
      END DO
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
      DO k = 2,qdims%k_end

        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
           IF ( P_LCL(ii)  <   P(I,j,K) ) THEN
             NLCL(I,j) = K-1
             NLCL_c(ii) = K-1
             Z_LCL(I,j) = Z_HALF_c(ii,NLCL(I,j)+1)
             Z_LCL_UV(I,j) = z_full_c(ii,NLCL(I,j))
           END IF
        END DO
      END DO

!-----------------------------------------------------------------------
!! 0.3 Now compare plume s_VL with each model layer s_VL in turn to
!!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------
!
      DO  k = 2,qdims%k_end

! DEPENDS ON: qsat_mix
        Call QSAT_MIX(QS(1,K),T(1,K),P_theta_lev_c(1,K),               &
                                             nunstable,l_mixing_ratio)

!CDIR nodep
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   

!-----------------------------------------------------------------------
! Only perform parcel ascent if unstable
!-----------------------------------------------------------------------

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                              &
               .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1
!
!-----------------------------------------------------------------------
! Set flag to true when level BELOW is above the lcl (in case P_LCL is
! in the upper half of layer NLCL+1). Implies ABOVE_LCL at NLCL+3.
!-----------------------------------------------------------------------
          IF ( k > 2 .and. P_LCL(ii)  <   P(I,j,K-1) ) THEN
            ABOVE_LCL(ii)=.FALSE.
          Else IF ( k == 2 .and. P_LCL(ii)  <   Pstar(I,j) ) THEN
            ABOVE_LCL(ii)=.FALSE.
          Else
            ABOVE_LCL(ii)=.TRUE.
          END IF
!
!         !-----------------------------------------------------------
!         ! calculate parcel potential temperature by linearising
!         ! q_sat about the environmental temperature.
!         !-----------------------------------------------------------
!
          IF(T(ii,K) >  TM) THEN
            DQ_SAT_ENV=repsilon*LC*QS(ii,K)/(R*T(ii,K)**2)
            Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -              &
              DQ_SAT_ENV*( SL_plume(ii)-gamma_dry*z_full_c(ii,K)-T(ii,K) )  &
                                   ) / (1.0+LCRCP*DQ_SAT_ENV) )
            Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*   &
               ( TL(ii,K) - T(ii,K) ) ) / (1.0+LCRCP*DQ_SAT_ENV) )
! add on the difference in the environment's ql as calculated by the
! partial condensation scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
            Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)        &
                           - Q_LIQ_ENV
            T_PARC(ii,k) =                                             &
                   SL_plume(ii)-gamma_dry*z_full_c(ii,K)+LCRCP*Q_LIQ_PARC
          Else       
            DQ_SAT_ENV=repsilon*LS*QS(ii,K)/(R*T(ii,K)**2)
            Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -              &
              DQ_SAT_ENV*( SL_plume(ii)-gamma_dry*z_full_c(ii,K)-T(ii,K) )  &
                                   ) / (1.0+LSRCP*DQ_SAT_ENV) )
            Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*   &
               ( TL(ii,K) - T(ii,K) ) ) / (1.0+LSRCP*DQ_SAT_ENV) )
! add on difference in environment's ql between RH_CRIT and RH_CRIT=1
            Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)        &
                           - Q_LIQ_ENV
            T_PARC(ii,k) =                                             &
                   SL_plume(ii)-gamma_dry*z_full_c(ii,K)+LSRCP*Q_LIQ_PARC
          END IF       ! test on TM
          Q_VAP_PARC=QW_plume(ii)-Q_LIQ_PARC
!
          T_DENS_PARC(ii,k)=T_PARC(ii,k)*                              &
                             (1.0+C_VIRTUAL*Q_VAP_PARC-Q_LIQ_PARC)
          T_DENS_ENV(ii,k)=T(ii,K)*                                    &
                     (1.0+C_VIRTUAL*Q_c(ii,K)-QCL_c(ii,K)-QCF_C(ii,K))

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------
          DPAR_BYDZ2 = (T_DENS_PARC(ii,k) + gamma_dry*z_full_c(ii,K) -      &
                       PAR_SVL_KM1(ii)) /                                   &
                    (z_full_c(ii,K) - z_full_c(ii,K-1))
          DENV_BYDZ2 = (T_DENS_ENV(ii,k) + gamma_dry*z_full_c(ii,K) -       &
                       ENV_SVL_KM1(ii)) /                                   &
                    (z_full_c(ii,K) - z_full_c(ii,K-1))

          IF ( TOPBL(ii) .AND. ( DENV_BYDZ2  <   DPAR_BYDZ2 .OR.       &
                     PAR_SVL_KM1(ii) - ENV_SVL_KM1(ii)  <=  0.0 )      &
                                    .AND. .NOT. SHMIN(ii) ) THEN
              SHMIN(ii) = .TRUE.
              DT_DENS_PARC_TMIN(ii) = PAR_SVL_KM1(ii) -                &
                                              ENV_SVL_KM1(ii)
              KSHMIN(ii) = K-1
          END IF
          IF ( .NOT.TOPBL(ii) .AND. K  >   K_plume(ii) .AND.           &
        (  ( T_DENS_PARC(ii,k)-T_DENS_ENV(ii,k)  <=  - THV_PERT(ii) )  &
!
!                      plume non buoyant
!
          .OR. (ABOVE_LCL(ii) .AND. (DENV_BYDZ2  >   1.25*DPAR_BYDZ2)) &
!
!                      or environmental virtual temperature gradient
!                      significantly larger than parcel gradient
!                      above lifting condensation level
!
                 .OR. (k  >   qdims%k_end-1)                           &
!                      or reached top of model
               )                                                       &
               ) THEN
!
            TOPBL(ii) = .TRUE.
            ZH(I,j) = Z_HALF_c(ii,K)
            NTML(i,j) = K-1
            DT_DENS_PARC_T(ii) = PAR_SVL_KM1(ii) - ENV_SVL_KM1(ii)
            IF ( DELTHVU(I,j)  >   0.0) THEN
! compensate for any negative buoyancy of parcel in cloud layer
              DELTHVU(I,j) = DELTHVU(I,j) - DTV_MIN(ii) *              &
                                          ( ZH(I,j) - Z_LCL(I,j) )
            END IF

! Added to estimate undilute CAPE and CIN. Note for this option the CAPE
! may not always be that of an ascent reaching the level of neutral buoyancy.

            inc = g * (t_dens_parc(ii,k) - t_dens_env(ii,k))           &
                  *(z_half_c(ii,k+1) - z_half_c(ii,k))/t_dens_env(ii,k)

            IF (inc < 0.0) THEN
              CIN_c(ii) = CIN_c(ii) + inc
            ELSE     ! CAPE holds only postive part       
              CAPE_c(ii) = CAPE_c(ii) + inc           
            END IF          

          END IF           ! test on topbl

          ENV_SVL_KM1(ii) = T_DENS_ENV(ii,k) + gamma_dry*z_full_c(ii,K)
          PAR_SVL_KM1(ii) = T_DENS_PARC(ii,k) + gamma_dry*z_full_C(ii,K)
!
          IF (K  >   NLCL(I,j) .AND. .NOT. TOPBL(ii)) THEN
            DTV_MIN(ii) = MIN( DTV_MIN(ii),                            &
                           (T_DENS_PARC(ii,k)-T_DENS_ENV(ii,k)) *      &
                           ( (100000.0/P_theta_lev_c(ii,K))**kappa ) )
            DELTHVU(I,j) = DELTHVU(I,j) + (T_DENS_PARC(ii,k)           &
                                   -T_DENS_ENV(ii,k)) *                &
                           ( (100000.0/P_theta_lev_c(ii,K))**kappa ) * &
                           ( Z_HALF_c(ii,K+1) - Z_HALF_c(ii,K) )
          END IF
!

        END DO   ! ii loop
      END DO     ! k loop
!-----------------------------------------------------------------------
!! 0.4 Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        ZHPAR(I,j) = ZH(I,j)
        NTPAR(I,j) = NTML(I,j)
      END DO
      END DO

!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
      Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     if lifting condensation level lower than parcel ascent, and is
!     within BL_LEVELS, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any if LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
        IF ( NTML(I,j)-NLCL(I,j)  >=  2                                &
                              .AND. NLCL(I,j)  >   K_plume(ii)         &
                                    .AND. NLCL(I,j)  <   MBL-1 ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = ZH
!     For cumulus top of mixed layer = ZLCL
!     NTML >= MBL indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------
          IF (NTML(I,j)  >=  MBL) THEN
            CUMULUS(I,j) = .TRUE.
          Else

! TEMPORARY: SHOULD REALLY BE DONE WITH SCALE HEIGHT
            IF (NTML(I,j)  >   KCUCHECK(ii)                            &
                   .AND. NLCL(I,j)  <=  KCUCHECK(ii)-2) THEN

              GRAD_CLD =  ABS( QW(ii,KCUCHECK(ii)) -                   &
                                        QW(ii,NLCL(I,j)) ) /           &
                ( z_full_C(ii,KCUCHECK(ii)) - z_full_c(ii,NLCL(I,j)) )
            Else
              GRAD_CLD =  ABS( QW(ii,NTML(I,j)) -                      &
                                        QW(ii,NLCL(I,j)) ) /           &
                    ( z_full_c(ii,NTML(I,j)) - z_full_C(ii,NLCL(I,j)) )
            END IF

            GRAD_SUB =  ABS( QW(ii,NLCL(I,j)) -                        &
                                      QW(ii,K_plume(ii)) ) /           &
                 ( z_full_c(ii,NLCL(I,j)) - z_full_c(ii,K_plume(ii)) )

            IF (GRAD_CLD  >   1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!-----------------------------------------------------------------------

! TEMPORARY: SHOULD REALLY BE DONE WITH SCALE HEIGHT
              IF (NTML(I,j)  <=  KCUCHECK(ii)) THEN
              GRAD_CLD =  ABS( QW(ii,NTML(I,j)-1) - QW(ii,NLCL(I,j)) ) &
                 /( z_full_c(ii,NTML(I,j)-1) - z_full_c(ii,NLCL(I,j)) )
              END IF

              IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                CUMULUS(I,j) = .TRUE.
              END IF
            Else    

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identified as
! well-mixed)

! First check using level below (recalculate GRAD_SUB)

              IF (NLCL(I,j) - K_plume(ii)  >=  2) THEN
                 GRAD_SUB =  ABS( QW(ii,NLCL(I,j)-1) -                 &
                                      QW(ii,K_plume(ii)) ) /           &
                 ( z_full_c(ii,NLCL(I,j)-1) - z_full_c(ii,K_plume(ii)) )

                 IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
                   CUMULUS(I,j) =.TRUE.
                 END IF

              END IF

! If still diagnosing well-mixed, check using level above
! (recalculate GRAD_CLD)

              IF (.NOT. CUMULUS(I,j) ) THEN

               IF (NTML(I,j)  >   KCUCHECK(ii)                         &
                   .AND. NLCL(I,j)  <=  KCUCHECK(ii)-2) THEN

                GRAD_CLD =  ABS( QW(ii,KCUCHECK(ii)) -                 &
                                        QW(ii,NLCL(I,j)+1) ) /         &
               ( z_full_c(ii,KCUCHECK(ii)) - z_full_c(ii,NLCL(I,j)+1) )
               Else
                GRAD_CLD =  ABS( QW(ii,NTML(I,j)) -                    &
                                        QW(ii,NLCL(I,j)+1) ) /         &
                  ( z_full_c(ii,NTML(I,j)) - z_full_c(ii,NLCL(I,j)+1) )
               END IF

               IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
                 CUMULUS(I,j) =.TRUE.
               END IF

              END IF    ! not cumulus
            END IF      ! cloud gradient test
          END IF        ! ntml > mbl test
        END IF        !  level test

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!-----------------------------------------------------------------------
        K=NLCL(I,j)
        IF ( LAND_MASK(I,j) .AND. CUMULUS(I,j) .AND.                   &
                                       NTPAR(I,j)  <   MBL ) THEN
          DO WHILE ( K  <=  NTPAR(I,j) .AND. CLOUD_FRACTION(I,j,K)     &
                                                  >=  SC_CFTOL )
            K = K + 1
          END DO
          IF (K  ==  NTPAR(I,j)+1) CUMULUS(I,j) = .FALSE.
        END IF


        IF ( CUMULUS(I,j) ) THEN

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
          IF ( SHMIN(ii) ) THEN

            IF ( w_copy(I,j,BL_LEVELS)  <   0.0 .AND.                  &
                 (z_full_c(ii,NTPAR(I,j))  <=  2500.0 .OR.             &
                  T(ii,NTPAR(I,j))  >=  TM)                            &
             .AND. (z_full_c(ii,KSHMIN(ii)) - z_full_c(ii,NTPAR(I,j))) &
                <=  1.25*(ZHPAR(I,j) - Z_LCL(I,j)) .AND.               &
               DT_DENS_PARC_TMIN(ii)  <=  0.55*DT_DENS_PARC_T(ii) )    &
            THEN

              L_SHALLOW(I,j) = .TRUE.
            END IF

          END IF

!-----------------------------------------------------------------------
!      Set mixed layer depth to Z_LCL
!-----------------------------------------------------------------------
          IF (P_LCL(ii)  <   (P_theta_lev(I,j,NLCL(I,j)+1))) THEN
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set Z_LCL to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
             NLCL(I,j) = NLCL(I,j)+1
             Z_LCL(I,j) = Z_HALF(I,j,NLCL(I,j)+1)
             Z_LCL_UV(I,j)= z_full(I,j,NLCL(I,j))
          END IF
          ZH(I,j) = Z_LCL(I,j)
          NTML(I,j) = NLCL(I,j)

!      If CUMULUS has been diagnosed but DELTHVU is negative, reset
!      CUMULUS and L_SHALLOW to FALSE but leave ZH and NTML at LCL

            IF (DELTHVU(I,j)  <=  0.0) THEN

              CUMULUS(I,j)   = .FALSE.
              L_SHALLOW(I,j) = .FALSE.

            END IF

        Else

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within BL_LEVELS
!-----------------------------------------------------------------------
          IF (NTML(I,j)  >   MBL) THEN
            NTML(I,j) = MBL
            NTPAR(I,j) = MBL
            ZH(I,j) = Z_HALF(I,j,MBL+1)
            ZHPAR(I,j) = ZH(I,j)
          END IF
          IF (NLCL(I,j)  >   MBL) THEN
              NLCL(I,j) = MBL
              Z_LCL(I,j) = ZH(I,j)
              Z_LCL_UV(I,j)=z_full(I,j,MBL-1)
          END IF

        END IF       ! cumulus test
      END DO

! Expand CAPE and CIN values back to full arrays
      Do ii=1,nunstable
        i = index_i(ii)
        j = index_j(ii)
        CAPE(i,j) = CAPE_c(ii)
        CIN(i,j)  = CIN_c(ii)
      END DO



!=======================================================================
! Use  new code  to calculate both parcel and ntpar differently
! This code is taken from 5A scheme
!=======================================================================
      Else if (icvdiag == 1 .OR. icvdiag == 5 ) then


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(z_surf, ii, lrcp_const, l_const, &
!$OMP& dq_sat_env, q_liq_parc, q_liq_env, q_vap_parc, dz, z_pr, i, j,   &
!$OMP& th_par, grad_cld, grad_sub, k, vap_press, lrcp_const_parc, inc)

!$OMP DO SCHEDULE(DYNAMIC,1)
      Do ii=1, nunstable

        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
         z_surf = 0.1 * zh_c(ii)

         Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
                    SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

            k_plume(ii) = k_plume(ii) + 1

         END DO
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
      DO ii=1, nunstable

        nlcl_min(ii) = 2

!-----------------------------------------------------------------------
! Convection scheme requires NLCL at least 2
! Also require ZLCL > 150m (approx nlcl=2 for G3 levels)
!-----------------------------------------------------------------------
        K=3
        DO WHILE ( z_half_c(ii,k) < 150.0 .AND. K < MBL ) 
          K=K+1
        END DO
        nlcl_min(ii) = K-1

      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do ii=1, nunstable

         sl_plume(ii) = TL(ii,k_plume(ii))                              &
                              + gamma_dry * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                   &
                        min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                                &
                              /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

      END DO
!$OMP END DO

!
!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/repsilon   q specific humidity
!   vapour pressure e ~ qp/(repsilon+q)   q mixing ratio

      If (l_mixing_ratio) THEN
!$OMP DO SCHEDULE(STATIC)
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
                                       P_theta_lev_c(ii,k_plume(ii))   &
                            / (repsilon+Q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
           P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO
!$OMP END DO

      Else
! expression for specific humidity
!$OMP DO SCHEDULE(STATIC)
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          vap_press = Q_c(ii,k_plume(ii)) *                            &
               P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*repsilon )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
            P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO
!$OMP END DO

      END IF ! test on l_mixing_ratio  
!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
!$OMP DO SCHEDULE(STATIC)
      do j=1,rows
        do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        END DO
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
        zh2(ii)  = zh_c(ii)
      END DO
!$OMP END DO

!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  k = 2,qdims%k_end

!$OMP DO SCHEDULE(STATIC)
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If ( P_LCL(ii)  <   P(i,j,K) ) THEN
! compressed copies
            nlcl_c(ii) = K-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = K-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          END IF     ! test on p_lcl
        END DO       ! ii loop
!$OMP END DO

      END DO         ! k loop

!-----------------------------------------------------------------------
! 4.0 Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Note initial parcel conditions different from those used in the G-R
! mass flux convection scheme.
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

       DO  K = 1,qdims%k_end

! Require t_ref on all point for qsat call
!$OMP DO SCHEDULE(STATIC)
         Do ii=1, nunstable
           t_ref(ii) = th_ref(ii)*exner_theta_levels_c(ii,k)
         END DO
!$OMP END DO

!$OMP BARRIER
!$OMP SINGLE
! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)               &
                                          ,nunstable,l_mixing_ratio)
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
         Do ii=1, nunstable
           If(T_ref(ii) >  TM) THEN
             lrcp_const = lcrcp
             l_const    = lc
           Else
             lrcp_const = lsrcp
             l_const    = ls
           END IF

! dqsat/dT - same whether q specific humidity or mixing ratio

           dq_sat_env = repsilon*l_const*qsat_lev(ii)/(R*T_ref(ii)**2)

           q_liq_parc = max( 0.0, ( qw_plume(ii) - QSat_lev(ii)             &
             -dq_sat_env*( sl_plume(ii)-gamma_dry*z_full_c(ii,K)-T_ref(ii) )&
                                   ) / (1.0+lrcp_const*dq_sat_env) )
           q_liq_env  = max( 0.0, ( qw(ii,K) - QSat_lev(ii)            &
              -dq_sat_env*( TL(ii,K)               - T_ref(ii) )       &
                                   ) / (1.0+Lrcp_const*dq_sat_env) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
           q_liq_parc = q_liq_parc + qcl_c(ii,k)                            &
                                  + qcf_c(ii,k)- q_liq_env
           T_parc(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,K)               &
                                        +lrcp_const*q_liq_parc

! May need to recalculate if T_parc is > Tm and T_ref < Tm

           If (T_ref(ii) <= TM.and.T_parc(ii,k) >  TM) THEN

! recalculate using corrected latent heats
             lrcp_const_parc = lcrcp
             q_liq_parc = max( 0.0, ( qw_plume(ii) - Qsat_lev(ii)           &
            -dq_sat_env*( sl_plume(ii)-gamma_dry*z_full_c(ii,K)-T_ref(ii) ) &
                           ) / (1.0+lrcp_const_parc*dq_sat_env) )
             q_liq_parc = q_liq_parc + qcl_c(ii,k)                     &
                                     + qcf_c(ii,k)- q_liq_env

! revised at parcel calculation

             T_parc(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,K)             &
                                        +lrcp_const_parc*q_liq_parc

           END IF

           q_vap_parc=qw_plume(ii)-q_liq_parc
!
           t_dens_parc(ii,k)=T_parc(ii,k)*                             &
                          (1.0+c_virtual*q_vap_parc-q_liq_parc)


           t_dens_env(ii,k)=T(ii,K)*                                   &
                       (1.0+c_virtual*Q_c(ii,K)-qcl_c(ii,k)-qcf_c(ii,k))

           buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

           env_svl(ii,k) = t_dens_env(ii,k)  + gamma_dry*z_full_c(ii,K)
           par_svl(ii,k) = t_dens_parc(ii,k) + gamma_dry*z_full_c(ii,K)

           If (k >= 2) THEN

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------

             dz = z_full_c(ii,K) - z_full_c(ii,K-1)

             dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
             denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

           END IF   ! test on k

! calculate t_ref for next level
           IF (k > 1 .AND. k < qdims%k_end-1) THEN
             z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                  &
                               /(z_full_c(ii,k)-z_full_c(ii,k-1))
             th_par = T_parc(ii,k)/exner_theta_levels_c(ii,k)
             th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

! Check sensible value otherwise set to previous reference value
! Problems can occur near top of model where calculation are nolonger 
! important.
             If (th_ref(ii) < 0.0) THEN
               th_ref(ii) = th_par_km1(ii)
             END IF
             If (th_par > 0.0) THEN   
               th_par_km1(ii) = th_par
             END IF
  
          END IF

        END DO    ! ii loop
!$OMP END DO
      END DO      ! level loop

!-----------------------------------------------------------------------
! tests on parcel ascent
!-----------------------------------------------------------------------
!   Now compare plume s_VL with each model layer s_VL in turn to
!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------

      Do  k = 2,qdims%k_end

!-----------------------------------------------------------------------
! Only perform tests if parcel ascent If unstable
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
        Do ii=1,nunstable

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                              &
               .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1

! freezing level

          If (t(ii,k) <  TM.and.t(ii,k-1) >= TM) THEN  
            If (freeze_lev(ii) == 1) THEN
              freeze_lev(ii) = k    
            END IF
          END IF

!-----------------------------------------------------------------------
! Set flag to true when level below is at least one level above the lcl
! and above the lcl transition zone
! Code implies ABOVE_LCL at NLCL+3 or greater.

          If (k-1 >  nlcl_c(ii)+1                                      &
                       .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN
            above_lcl(ii)=.true.
          Else
            above_lcl(ii)=.false.
          END IF

!-----------------------------------------------------------------------
! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
!-----------------------------------------------------------------------
! Not reached LNB continue testing

          If ( .not.topprof(ii).and.k >  k_plume(ii) )THEN

            If (buoyancy(ii,k) >  max_buoy(ii)) THEN
              max_buoy(ii) = buoyancy(ii,k)
              k_max(ii)    = k  
            END IF 

! Is parcel still buoyant ?

            If ( (buoyancy(ii,k)  <=  - thv_pert(ii))                  &
!                      or reached top of model
                 .OR. (k  >   qdims%k_end-1)  ) THEN

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,K)

! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

              If ( delthvu_c(ii)  >   0.0) THEN
! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
                                      ( z_half_c(ii,K) - z_lcl_c(ii) )
              END IF                                                     
            END IF
          END IF

!-----------------------------------------------------------------------
! Tests applied once found top of parcel ascent.
! Aim - to establish if the ascent has an inversion above the top
!       i.e. the ascent may indicate shallow /congestus convection.
! Sets indicator shmin = .true. if conditions met and stops testing.
!
! Conditions are ;
! either  denv/dz(k) < dpar/dz(k)
!   or    par_svl(k-1) -env_svl(k-1) <= 0.0
!
!-----------------------------------------------------------------------

          If ( topbl(ii) .and. ( denv_bydz(ii,k)  <   dpar_bydz(ii,k)  &
                      .OR. buoyancy(ii,k-1)  <=  0.0 )                 &
                                    .and. .not. shmin(ii) ) THEN
            shmin(ii) = .TRUE.
            Dt_dens_parc_TMIN(ii) = buoyancy(ii,k-1)
            kshmin(ii) = K-1
          END IF

!-----------------------------------------------------------------------
! Tests applied to find parcel top
!
!-----------------------------------------------------------------------

          If ( .not.topbl(ii) .and. K  >   k_plume(ii) .and.           &
       (  ( buoyancy(ii,k) <=  - thv_pert(ii)).OR.                     &
!
!                      plume non buoyant
!
       (above_lcl(ii).and.(denv_bydz(ii,k) >  1.25*dpar_bydz(ii,k)))   &
!
!                      or environmental virtual temperature gradient
!                      signIficantly larger than parcel gradient
!                      above lIfting condensation level
!
                 .OR. (k  >   qdims%k_end-1)                           &
!                      or reached top of model
               )                                                       &
               ) THEN
!
            topbl(ii) = .TRUE.
            zh2(ii) = z_half_c(ii,K)
            k_inv(ii) = K-1

            Dt_dens_parc_T2(ii) = buoyancy(ii,k-1)
            If ( delthvu2(ii)  >   0.0) THEN
! compensate for any negative buoyancy of parcel in cloud layer
              delthvu2(ii) = delthvu2(ii) - dtv_min2(ii) *             &
                                    ( z_half_c(ii,k) - z_lcl_c(ii) )
            END IF
          END IF          ! test on .not.topbl

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          If (k > nlcl_c(ii) .and. k < qdims%k_end ) THEN

            inc = g  * buoyancy(ii,k)                                  &
                * (z_half_c(ii,K+1) - z_half_c(ii,K))/t_dens_env(ii,k)

            !---------------------------------------------------------- 
            ! If not reached an inversion or level of neutral buoyancy
            !---------------------------------------------------------- 

            If (.not. topbl(ii)) THEN

              dtv_min2(ii) = MIN( dtv_min2(ii),                        &
                             buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

              delthvu2(ii) = delthvu2(ii) +                            &
                    buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K)) & 
                              /exner_theta_levels_c(ii,k)
            END IF

            !---------------------------------------------------------- 
            ! If not reached level of neutral buoyancy (top of ascent)
            !---------------------------------------------------------- 

            If (.not. topprof(ii)) THEN
                
              ! Note only calculating CIN and CAPE from ascents reaching
              ! level of neutral buoyancy. This may not always correspond 
              ! to the diagnosed top for the convection scheme.

              IF (inc <  0.0) THEN
                CIN_c(ii)  = CIN_c(ii) + inc
              ELSE      ! CAPE holds only positive part
                CAPE_c(ii) = CAPE_c(ii) + inc
              END IF

              dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

              delthvu_c(ii) = delthvu_c(ii) +                         &
                  buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K))  & 
                              /exner_theta_levels_c(ii,k)
 
            END IF    ! test on topprof

          END IF

!-----------------------------------------------------------------------
        END DO   ! ii loop
!$OMP END DO
      END DO     ! level loop

!-----------------------------------------------------------------------
! Average vertical velocity over a layer  - required for shallow
!   convection test.
!-----------------------------------------------------------------------
! Layer from top in cloud w to value 1500km above cloud top?
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      END DO
!$OMP END DO

      DO k=1,tdims%k_end-1
!$OMP DO SCHEDULE(STATIC)
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= k_neutral(ii).and.                                   &
             z_full_c(ii,k) <= (z_half_c(ii,k_neutral(ii)+1)+1500.)) THEN

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)

          END IF
        END DO
!$OMP END DO
      END DO

!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable
        if (mass(ii)  >  0.0 ) THEN
          w_avg(ii) = w_avg(ii)/mass(ii)
        endif
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable
        w_avg2(ii) = 0.0
        mass(ii)   = 0.0
      END DO
!$OMP END DO

      DO k=1,tdims%k_end-1

!$OMP DO SCHEDULE(STATIC)
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= k_inv(ii) .and.                                     &
             z_full_c(ii,k) <= (z_half_c(ii,k_inv(ii)+1)+1500.)) THEN

            mass(ii)   = mass(ii) + dmass_theta(ii,k)
            w_avg2(ii) = w_avg2(ii) + w_copy(i,j,k)*dmass_theta(ii,k)
          END IF
        END DO
!$OMP END DO

      END DO

!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable
        If (mass(ii)  >  0.0 ) THEN
          w_avg2(ii) = w_avg2(ii)/mass(ii)
        END IF
      END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Default parcel top properties are assumed to be those when the
! ascent reaches the level of neutral buoyancy. These may not be those
! required in the case of shallow convection.
! Shallow convection requires the possible identifcation of an inversion
! at the top of the ascent. This may not be detected by the LNB test.
! The gradient tests are designed to detect the shallow top.
!-----------------------------------------------------------------------
! Modify top based on topbl test if ascent is likely to be shallow

!$OMP DO SCHEDULE(STATIC)
      Do ii=1,nunstable

        If (shmin(ii) ) THEN    ! found an inversion    
! points where k_inv not the same as k_neutral and level below freezing
! may be shallow or congestus or deep

          If (k_inv(ii) == k_neutral(ii)) THEN
!  Both methods give same answer for top level leave shmin set
            ntml_c(ii) = k_neutral(ii)

! Inversion top lower than level of neutral buoyancy.
! Check also, either below freezing level or less than 2500m for shallow 
! convection.
          Else if ((k_inv(ii) <  freeze_lev(ii) .or.                   &
                           z_full_C(ii,k_inv(ii)+1)  <=  2500.0 )      &
                   .and. k_inv(ii) <  k_neutral(ii) )THEN     


            If ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,k_inv(ii)))    &
               <=  1.25*(z_half_c(ii,k_inv(ii)+1) - z_lcl_c(ii)).and.  &
               (dt_dens_parc_tmin(ii)  <=  0.55*dt_dens_parc_t2(ii))   &
               .and.     (w_avg2(ii)  <   0.0)  ) THEN

! May be shallow or congestus
! set values to those found from inversion testing
               ntml_c(ii)  = k_inv(ii)
               delthvu_c(ii) = delthvu2(ii)
               zh_c(ii)    = zh2(ii)
               w_avg(ii)   = w_avg2(ii)
               dt_dens_parc_t(ii) = dt_dens_parc_t2(ii)

            Else   ! Assume not shallow or congestus
               ntml_c(ii) = k_neutral(ii)
               shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
            Endif
          Else   ! Assume deep  and therefore top LNB
              ntml_c(ii) = k_neutral(ii)
              shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
          Endif

        Else    !  No inversion found  i.e. shmin=false
          ntml_c(ii) = k_neutral(ii)
        Endif   ! shmin test

      END DO    ! ii loop
!$OMP END DO

! Is there any need to recalculate w_avg  for any points ?
! Assuming not but may need to check this.
!-----------------------------------------------------------------------
!!     Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
! Expand back up to full arrays 

!$OMP DO SCHEDULE(STATIC)
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do j = 1,rows
        Do i = 1, row_length
          zhpar(i,j) = zh(i,j)
        END DO
      END DO
!$OMP END DO
!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------

! loop over only unstable 
      If (icvdiag == 1) THEN

!$OMP SINGLE

!CDIR NODEP
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any If LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! 4A code  test ntml-nlcl >= 2 (5A code requires 3 cloud levels)

        If ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
                              .and. nlcl(i,j)  >   k_plume(ii)          &
                                    .and. nlcl(i,j)  <   MBL-1 ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     is above boundary layer indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (zhpar(I,j) >= 3000.0) THEN
            cumulus(i,j) = .TRUE.
          Else


! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine
! Assume moisture gradient tests stay the same whether specific or 
! mixing ratio though imply slightly different moisture gradients.

            If (ntml(i,j)  >   kcucheck(ii)                            &
                   .and. nlcl(i,j)  <=  kcucheck(ii)-2) THEN

              grad_cld =  ABS( QW(ii,kcucheck(ii)) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            Else
              grad_cld =  ABS( QW(ii,ntml_c(ii)) -                      &
                                        QW(ii,nlcl_c(ii)) ) /           &
                    ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            END IF

            grad_sub =  ABS( QW(ii,nlcl_c(ii)) -                        &
                                      QW(ii,k_plume(ii)) ) /            &
                 ( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            If (grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and THEN a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml_c(ii)  <=  kcucheck(ii)) THEN
              grad_cld =  ABS( QW(ii,ntml_c(ii)-1) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                 ( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              END IF

              If ( grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              END IF
            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl_c(ii) - k_plume(ii)  >=  2) THEN
                 grad_sub =  ABS( QW(ii,nlcl(i,j)-1) -                 &
                                      QW(ii,k_plume(ii)) ) /           &
                ( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 If ( grad_cld  >   1.10*grad_sub) THEN
                   cumulus(i,j) =.TRUE.
                 END IF

              END IF

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) THEN

               If (ntml_c(ii)  >   kcucheck(ii)                        &
                   .and. nlcl_c(ii)  <=  kcucheck(ii)-2) THEN

                grad_cld =  ABS( QW(ii,kcucheck(ii)) -                 &
                                        QW(ii,nlcl_c(ii)+1) ) /        &
               ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               Else
                grad_cld =  ABS( QW(ii,ntml_c(ii)) -                   &
                                        QW(ii,nlcl_c(ii)+1) ) /        &
                  ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               END IF

               If ( grad_cld  >   1.10*grad_sub) THEN
                 cumulus(i,j) =.TRUE.
               END IF

              END IF
            END IF
          END IF
        END IF
      END DO       ! ii loop 

!$OMP END SINGLE

! loop over only unstable 
      Else If (icvdiag == 5) THEN

!$OMP DO SCHEDULE(STATIC)      
!CDIR NODEP
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any If LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! 4A code  test ntml-nlcl >= 2 (5A code requires 3 cloud levels)

        If ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
          .and. nlcl(i,j)  <   MBL-1 ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     is above boundary layer indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (zhpar(I,j) >= 3000.0) then
            cumulus(i,j) = .TRUE.
            !nlcl is not permitted to be less than 2 if Cu is diagnosed.
            nlcl(i,j)    = max(2, nlcl(i,j))
            z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
            z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

          Else if (nlcl(i,j) > k_plume(ii)) then


! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine
! Assume moisture gradient tests stay the same whether specific or 
! mixing ratio though imply slightly different moisture gradients.

            If (ntml(i,j)  >   kcucheck(ii)                             &
                   .and. nlcl(i,j)  <=  kcucheck(ii)-2) then

              grad_cld =  ABS( QW(ii,kcucheck(ii)) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            Else
              grad_cld =  ABS( QW(ii,ntml_c(ii)) -                      &
                                        QW(ii,nlcl_c(ii)) ) /           &
                    ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            End If

            grad_sub =  ABS( QW(ii,nlcl_c(ii)) -                        &
                                      QW(ii,k_plume(ii)) ) /            &
                 ( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            If (grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and then a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml_c(ii)  <=  kcucheck(ii)) then
              grad_cld =  ABS( QW(ii,ntml_c(ii)-1) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                 ( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              End If

              If ( grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              End If
            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl_c(ii) - k_plume(ii)  >=  2) then
                 grad_sub =  ABS( QW(ii,nlcl(i,j)-1) -                 &
                                      QW(ii,k_plume(ii)) ) /           &
                ( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 If ( grad_cld  >   1.10*grad_sub) then
                   cumulus(i,j) =.TRUE.
                 End If

              End If

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) then

               If (ntml_c(ii)  >   kcucheck(ii)                        &
                   .and. nlcl_c(ii)  <=  kcucheck(ii)-2) then

                grad_cld =  ABS( QW(ii,kcucheck(ii)) -                 &
                                        QW(ii,nlcl_c(ii)+1) ) /        &
               ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               Else
                grad_cld =  ABS( QW(ii,ntml_c(ii)) -                   &
                                        QW(ii,nlcl_c(ii)+1) ) /        &
                  ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               End If

               If ( grad_cld  >   1.10*grad_sub) then
                 cumulus(i,j) =.TRUE.
               End If

              End If
            End If
          End If
        End If
      End Do       ! ii loop 
!$OMP END DO
      END IF     ! icvdiag
!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are done on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land.
!-----------------------------------------------------------------------
!$OMP DO
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        K=nlcl(i,j)

        If ( land_mask(i,j) .and. cumulus(i,j) .and.                    &
                                       ntpar(i,j)  <   MBL ) THEN
          Do while( K  <=  ntpar(i,j) .and. cloud_fraction(i,j,K)       &
                                                  >=  SC_CFTOL )
            K = K + 1
          END DO
          If (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        END IF

      END DO       ! ii loop 
!$OMP END DO

      !-----------------------------------------------------------------
      !      Final checks on the LCL
      !-----------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      DO ii=1,nunstable
        i = index_i(ii) 
        j = index_j(ii) 
        If ( cumulus(i,j) ) then
          If (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) then
            !-----------------------------------------------------------
            ! If LCL is diagnosed in the upper half of the layer set 
            ! z_lcl to the height of the upper layer interface
            ! (in code above LCL is always set to the lower interface).
            !-----------------------------------------------------------
            nlcl(i,j) = nlcl(i,j)+1
          End If

          !---------------------------------------------------------
          ! nlcl is not permitted to be less than nlcl_min
          !---------------------------------------------------------
          nlcl(i,j)    = max(nlcl_min(ii), nlcl(i,j))
          IF ( nlcl(i,j) >= ntpar(i,j)-1 ) THEN
            ! Cloud layer too shallow so rediagnose as well-mixed 
            cumulus(i,j)   = .false.
            L_shallow(i,j) = .false.
          END IF

          z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
          z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

        END IF

      END DO       ! ii loop 
!$OMP END DO

!-----------------------------------------------------------------------
! Original shallow diagnosis no congestus diagnosis
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
!CDIR NODEP
        Do ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          If ( cumulus(i,j) ) THEN

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
! Conditions for shallow convection
!    wadv < 0.0            (descending air)
!   top of parcel ascent < 2500. or T (top of parcel) > TM
!   height of min buoyancy (above Bl) - height of parcel top T level
!      <1.25(height parcel top - z lifting condensation level)
!   t_dens_parc -t_dens_env at kshmin <0.55t_dens_parc -t_dens_env
!   at ntpar
!
!   The last 2 conditions are looking for a strong inversion at the top
!   of the shallow cumulus.
!-----------------------------------------------------------------------
           If ( shmin(ii) ) THEN

            If ( w_avg(ii)  <   0.0 .and.                              &
              (z_full_c(ii,ntpar(i,j))  <=  2500.0 .OR.                &
                  T(ii,ntpar(i,j))  >=  TM)                            &
             .and. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar(i,j))) &
                <=  1.25*(zhpar(i,j) - z_lcl_c(ii)) .and.              &
               Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) )    &
            THEN

              L_shallow(i,j) = .TRUE.
! may be problem with ntpar diagnosis for deep if wadv test sets
! L_shallow  false

            END IF

           END IF       ! test on shmin

!-----------------------------------------------------------------------
!      Set mixed layer depth to z_lcl
!-----------------------------------------------------------------------
           zh(i,j) = z_lcl(i,j)
           ntml(i,j) = nlcl(i,j)

!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to FALSE but leave zh and ntml at LCL

           If (delthvu(i,j)  <=  0.0) THEN

             cumulus(i,j)   = .false.
             L_shallow(i,j) = .false.

           END IF

        Else      ! not cumulus

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within bl_levels
!-----------------------------------------------------------------------
          If (ntml(i,j)  >   MBL) THEN
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          END IF
          If (nlcl(i,j)  >   MBL) THEN
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          END IF

        END IF        ! test on cumulus

      END DO          ! ii loop  
!$OMP END DO

!$OMP END PARALLEL

!=======================================================================
! Use  new corrected code  - when every one uses this or dilute options
!                            THEN options 1 & 5 should be removed as 
!                            these contain a bug.
!                            Currently correcting option 6 introduced at UM7.5
!=======================================================================
    ELSE IF (icvdiag == 6) THEN

      DO ii=1, nunstable

        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
         z_surf = 0.1 * zh_c(ii)

         DO WHILE( z_full_c(ii,k_plume(ii))  <   z_surf .AND.           &
!                   ! not reached z_surf
                    SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

            k_plume(ii) = k_plume(ii) + 1

         END DO
      END DO

      DO ii=1, nunstable

        nlcl_min(ii) = 2

!-----------------------------------------------------------------------
! Convection scheme requires NLCL at least 2
! Also require ZLCL > 150m (approx nlcl=2 for G3 levels)
!-----------------------------------------------------------------------
        K=3
        DO WHILE ( z_half_c(ii,k) < 150.0 .AND. K < MBL ) 
          K=K+1
        END DO
        nlcl_min(ii) = K-1

      END DO

      DO ii=1, nunstable

         sl_plume(ii) = TL(ii,k_plume(ii))                              &
                              + gamma_dry * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                   &
                        min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                                &
                              /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

      END DO

!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/repsilon   q specific humidity
!   vapour pressure e ~ qp/(repsilon+q)   q mixing ratio

      IF (l_mixing_ratio) THEN
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
                                       P_theta_lev_c(ii,k_plume(ii))   &
                            / (repsilon+Q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
           P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      ELSE
! expression for specific humidity
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          vap_press = Q_c(ii,k_plume(ii)) *                            &
              P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*repsilon )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
                           (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                         - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                     ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
            P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      END IF ! test on l_mixing_ratio  
!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      DO j=1,rows
        Do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        END DO
      END DO
      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
        zh2(ii)  = zh_c(ii)
      END DO
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  k = 2,qdims%k_end
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If ( P_LCL(ii)  <   P(i,j,K) ) THEN
! compressed copies
            nlcl_c(ii) = K-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = K-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          END IF     ! test on p_lcl
        END DO       ! ii loop
      END DO         ! k loop

!-----------------------------------------------------------------------
! 4.0 Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Note initial parcel conditions different from those used in the G-R
! mass flux convection scheme.
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

       DO  k = 1,qdims%k_end

         ! Require t_ref on all point for qsat call
         DO ii=1, nunstable
           t_ref(ii) = th_ref(ii)*exner_theta_levels_c(ii,k)
         END DO

! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)               &
                                          ,nunstable,l_mixing_ratio)

! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_env,t(1,k),p_theta_lev_c(1,k)               &
                                          ,nunstable,l_mixing_ratio)
         Do ii=1, nunstable
           IF(T_ref(ii) >  TM) THEN
             lrcp_const = lcrcp
             l_const    = lc
           ELSE
             lrcp_const = lsrcp
             l_const    = ls
           END IF
           IF(T(ii,k) >  TM) THEN
             lrcp_const_env = lcrcp
             l_const_env    = lc
           ELSE
             lrcp_const_env = lsrcp
             l_const_env    = ls
           END IF

           ! dqsat/dT - same whether q specific humidity or mixing ratio

           dq_sat_env = repsilon*l_const_env*qsat_lev(ii)/(R*T_ref(ii)**2)

           dq_sat_par = repsilon*l_const*qsat_env(ii)/(R*T(ii,k)**2)

           q_liq_parc = max( 0.0, ( qw_plume(ii) - qsat_lev(ii)              &
                            -dq_sat_par*( sl_plume(ii)                       &
                                       -gamma_dry*z_full_c(ii,K)-T_ref(ii) ) &
                                  ) / (1.0+lrcp_const*dq_sat_par) )

           q_liq_env  = max( 0.0, ( qw(ii,K) - qsat_env(ii)                 &
                  -dq_sat_env*( TL(ii,K)               - T(ii,k) )          &
                                  ) / (1.0+lrcp_const_env*dq_sat_env) )

         ! add on the difference in the environment's ql as calculated by the
         ! UM cloud scheme (using some RH_CRIT value) and what it
         ! would be If RH_CRIT=1. This THEN imitates partial condensation
         ! in the parcel.

           q_liq_parc = q_liq_parc + qcl_c(ii,k) + qcf_c(ii,k)- q_liq_env
           T_parc(ii,k)=sl_plume(ii)                                        &
                         -gamma_dry*z_full_c(ii,K)+lrcp_const*q_liq_parc

          ! May need to recalculate if T_parc is > Tm and T_ref < Tm

           IF (T_ref(ii) <= TM .AND. T_parc(ii,k) >  TM) THEN

             ! recalculate using corrected latent heats
             lrcp_const_parc = lcrcp
             q_liq_parc = max( 0.0, ( qw_plume(ii) - Qsat_lev(ii)           &
                                  -dq_sat_par*( sl_plume(ii)                &
                                      -gamma_dry*z_full_c(ii,K)-T_ref(ii) ) &
                                ) / (1.0+lrcp_const_parc*dq_sat_par) )
             q_liq_parc = q_liq_parc + qcl_c(ii,k)+ qcf_c(ii,k)- q_liq_env

             ! revised at parcel calculation

             T_parc(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,K)             &
                                        +lrcp_const_parc*q_liq_parc

           END IF

           q_vap_parc=qw_plume(ii)-q_liq_parc

           t_dens_parc(ii,k)=T_parc(ii,k)*(1.0+c_virtual*q_vap_parc-q_liq_parc)

           t_dens_env(ii,k) =T(ii,K)*                                      &
                             (1.0+c_virtual*Q_c(ii,K)-qcl_c(ii,k)-qcf_c(ii,k))

           buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

           env_svl(ii,k) = t_dens_env(ii,k)  + gamma_dry*z_full_c(ii,K)
           par_svl(ii,k) = t_dens_parc(ii,k) + gamma_dry*z_full_c(ii,K)

           IF (k >= 2) THEN

          !-------------------------------------------------------------
          ! Find vertical gradients in parcel and environment SVL
          ! (using values from level below (i.e. K-1)).
          !-------------------------------------------------------------

             dz = z_full_c(ii,K) - z_full_c(ii,K-1)

             dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
             denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

           END IF   ! test on k

          ! calculate t_ref for next level
           IF (k >  1 .and. k <   qdims%k_end-1) THEN
             z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                  &
                               /(z_full_c(ii,k)-z_full_c(ii,K-1))
             th_par = T_parc(ii,k)/exner_theta_levels_c(ii,k)
             th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

          ! Check sensible value otherwise set to previous reference value
          ! Problems can occur near top of model where calculation are nolonger 
          ! important.
             IF (th_ref(ii) < 0.0) THEN
               th_ref(ii) = th_par_km1(ii)
             END IF
             IF (th_par > 0.0) THEN
               th_par_km1(ii) = th_par
             END IF
  
           END IF

        END DO    ! ii loop
      END DO      ! level loop

      !-----------------------------------------------------------------------
      ! tests on parcel ascent
      !-----------------------------------------------------------------------
      !   Now compare plume s_VL with each model layer s_VL in turn to
      !     find the first time that plume has negative buoyancy.
      !-----------------------------------------------------------------------

      DO  k = 2,qdims%k_end

      !-----------------------------------------------------------------------
      ! Only perform tests if parcel ascent If unstable
      !-----------------------------------------------------------------------
        DO ii=1,nunstable

        !  Find level just below 2.5km (for use in Cu diagnosis)

          IF ( z_full_c(ii,K)  >   2500.0                              &
               .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1

        ! freezing level

          IF (t(ii,k) <  TM.and.t(ii,k-1) >= TM) THEN  
            If (freeze_lev(ii) == 1) THEN
              freeze_lev(ii) = k    
            END IF
          END IF

      !-----------------------------------------------------------------------
      ! Set flag to true when level below is at least one level above the lcl
      ! and above the lcl transition zone
      ! Code implies ABOVE_LCL at NLCL+3 or greater.

          IF (k-1 >  nlcl_c(ii)+1                                      &
                      .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN
            above_lcl(ii)=.true.
          ELSE
            above_lcl(ii)=.false.
          END IF

      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

          IF ( .not.topprof(ii).and.k >  k_plume(ii) )THEN

            If (buoyancy(ii,k) >  max_buoy(ii)) THEN
              max_buoy(ii) = buoyancy(ii,k)
              k_max(ii)    = k  
            END IF 

            ! Is parcel still buoyant ?

            IF ( (buoyancy(ii,k)  <=  - thv_pert(ii))                  &
            !    or reached top of model
                .OR. (k  >   qdims%k_end-1)  ) THEN

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,K)

             ! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

              IF ( delthvu_c(ii)  >   0.0) THEN
             ! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
                                    ( z_half_c(ii,K) - z_lcl_c(ii) )
              END IF                                                     
            END IF
          END IF

!-----------------------------------------------------------------------
! Tests applied once found top of parcel ascent.
! Aim - to establish if the ascent has an inversion above the top
!       i.e. the ascent may indicate shallow /congestus convection.
! Sets indicator shmin = .true. if conditions met and stops testing.
!
! Conditions are ;
! either  denv/dz(k) < dpar/dz(k)
!   or    par_svl(k-1) -env_svl(k-1) <= 0.0
!
!-----------------------------------------------------------------------

          IF ( topbl(ii) .AND. ( denv_bydz(ii,k)  <   dpar_bydz(ii,k)  &
                      .OR. buoyancy(ii,k-1)  <=  0.0 )                 &
                                    .AND. .NOT. shmin(ii) ) THEN
            shmin(ii) = .TRUE.
            Dt_dens_parc_TMIN(ii) = buoyancy(ii,k-1)
            kshmin(ii) = K-1
          END IF

        !-----------------------------------------------------------------
        ! Tests applied to find parcel top
        !-----------------------------------------------------------------

          IF ( .not.topbl(ii) .and. K  >   k_plume(ii) .and.           &
       (  ( buoyancy(ii,k) <=  - thv_pert(ii)).OR.                     &
        !
        !                      plume non buoyant
        !
       (above_lcl(ii).and.(denv_bydz(ii,k) >  1.25*dpar_bydz(ii,k)))   &

        !       or environmental virtual temperature gradient
        !       significantly larger than parcel gradient
        !       above lifting condensation level

                 .OR. (k  >   qdims%k_end-1)                           &
        !        or reached top of model
               )                                                       &
               ) THEN

            topbl(ii) = .TRUE.
            zh2(ii) = z_half_c(ii,K)
            k_inv(ii) = K-1

            Dt_dens_parc_T2(ii) = buoyancy(ii,k-1)
            IF ( delthvu2(ii)  >   0.0) THEN
        ! compensate for any negative buoyancy of parcel in cloud layer
              delthvu2(ii) = delthvu2(ii) - dtv_min2(ii) *             &
                                    ( z_half_c(ii,k) - z_lcl_c(ii) )
            END IF
          END IF          ! test on .not.topbl

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          IF (k > nlcl_c(ii) .AND. k < qdims%k_end ) THEN

            inc = g  * buoyancy(ii,k)                             &
                * (z_half_c(ii,K+1) - z_half_c(ii,K))/t_dens_env(ii,k)

            !---------------------------------------------------------- 
            ! If not reached an inversion or level of neutral buoyancy
            !---------------------------------------------------------- 

            IF (.NOT. topbl(ii)) THEN

              dtv_min2(ii) = MIN( dtv_min2(ii),                        &
                             buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

              delthvu2(ii) = delthvu2(ii) +                            &
                    buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K)) & 
                              /exner_theta_levels_c(ii,k)
            END IF

            !---------------------------------------------------------- 
            ! If not reached level of neutral buoyancy (top of ascent)
            !---------------------------------------------------------- 

            IF (.NOT. topprof(ii)) THEN
                
              ! Note only calculating CIN and CAPE from ascents reaching
              ! level of neutral buoyancy. This may not always correspond 
              ! to the diagnosed top for the convection scheme.

              IF (inc <  0.0) THEN
                CIN_c(ii)  = CIN_c(ii) + inc
              ELSE
                CAPE_c(ii) = CAPE_c(ii) + inc
              END IF 
              dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

              delthvu_c(ii) = delthvu_c(ii) +                         &
                  buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K))  & 
                             /exner_theta_levels_c(ii,k)
 
            END IF    ! test on topprof

          END IF

       !-----------------------------------------------------------------------
        END DO   ! ii loop

      END DO     ! level loop

     !-----------------------------------------------------------------------
     ! Average vertical velocity over a layer  - required for shallow
     !   convection test.
     !-----------------------------------------------------------------------
     ! Layer from top in cloud w to value 1500km above cloud top?
     !-----------------------------------------------------------------------
      DO ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      END DO
      DO k=1,tdims%k_end-1
        DO ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          IF (k >= k_neutral(ii).AND.                                     &
            z_full_c(ii,k) <= (z_half_c(ii,k_neutral(ii)+1)+1500.)) THEN

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)

          END IF
        END DO
      END DO
      DO ii=1,nunstable
        IF (mass(ii)  >  0.0 ) THEN
          w_avg(ii) = w_avg(ii)/mass(ii)
        endif
      END DO

      DO ii=1,nunstable
        w_avg2(ii) = 0.0
        mass(ii)   = 0.0
      END DO
      DO k=1,tdims%k_end-1
        DO ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= k_inv(ii) .AND.                                     &
             z_full_c(ii,k) <= (z_half_c(ii,k_inv(ii)+1)+1500.)) THEN

            mass(ii)   = mass(ii) + dmass_theta(ii,k)
            w_avg2(ii) = w_avg2(ii) + w_copy(i,j,k)*dmass_theta(ii,k)
          END IF
        END DO
      END DO
      DO ii=1,nunstable
        IF (mass(ii)  >  0.0 ) THEN
          w_avg2(ii) = w_avg2(ii)/mass(ii)
        END IF
      END DO

    !-----------------------------------------------------------------------
    ! Default parcel top properties are assumed to be those when the
    ! ascent reaches the level of neutral buoyancy. These may not be those
    ! required in the case of shallow convection.
    ! Shallow convection requires the possible identifcation of an inversion
    ! at the top of the ascent. This may not be detected by the LNB test.
    ! The gradient tests are designed to detect the shallow top.
    !-----------------------------------------------------------------------
    ! Modify top based on topbl test if ascent is likely to be shallow

      DO ii=1,nunstable

        IF (shmin(ii) ) THEN    ! found an inversion    
    ! points where k_inv not the same as k_neutral and level below freezing
    ! may be shallow or congestus or deep

          IF (k_inv(ii) == k_neutral(ii)) THEN
    !  Both methods give same answer for top level leave shmin set
            ntml_c(ii) = k_neutral(ii)

    ! Inversion top lower than level of neutral buoyancy.
    ! Check also, either below freezing level or less than 2500m for shallow 
    ! convection.
          ELSE IF ((k_inv(ii) <  freeze_lev(ii) .OR.                   &
                           z_full_C(ii,k_inv(ii)+1)  <=  2500.0 )      &
                   .AND. k_inv(ii) <  k_neutral(ii) )THEN     


            IF ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,k_inv(ii)))    &
               <=  1.25*(z_half_c(ii,k_inv(ii)+1) - z_lcl_c(ii)).AND.  &
               (dt_dens_parc_tmin(ii)  <=  0.55*dt_dens_parc_t2(ii))   &
               .AND.     (w_avg2(ii)  <   0.0)  ) THEN

    ! May be shallow or congestus
    ! set values to those found from inversion testing
               ntml_c(ii)  = k_inv(ii)
               delthvu_c(ii) = delthvu2(ii)
               zh_c(ii)    = zh2(ii)
               w_avg(ii)   = w_avg2(ii)
               dt_dens_parc_t(ii) = dt_dens_parc_t2(ii)

            ELSE   ! Assume not shallow or congestus
               ntml_c(ii) = k_neutral(ii)
               shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
            END IF
          ELSE   ! Assume deep  and therefore top LNB
              ntml_c(ii) = k_neutral(ii)
              shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
          END IF

        ELSE    !  No inversion found  i.e. shmin=false
          ntml_c(ii) = k_neutral(ii)
        END IF   ! shmin test

      END DO    ! ii loop

     ! Is there any need to recalculate w_avg  for any points ?
     ! Assuming not but may need to check this.
     !-----------------------------------------------------------------------
     !     Save parcel ascent top: this will be used to allow mixing and
     !     entrainment into decoupled Sc of single layer thickness when it
     !     occurs above Cu.
     !-----------------------------------------------------------------------
     ! Expand back up to full arrays 

      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      END DO

      DO j = 1,rows
        DO i = 1, row_length
          zhpar(i,j) = zh(i,j)
        END DO
      END DO
     !-----------------------------------------------------------------------
     !     Test height derived above against lifting condensation level
     !-----------------------------------------------------------------------
     ! loop over only unstable 

      
!CDIR NODEP
      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
      !-----------------------------------------------------------------------
      !     Check lifting condensation level against height of parcel ascent.
      !     If lifting condensation level lower than parcel ascent, and is
      !     within bl_levels, then decide on type of cloudy layer.
      !     Gradient tests and cumulus parametriztion require a minimum number 
      !     of grid-levels between LCL and top of parcel ascent, otherwise 
      !     define as stratocumulus.  To avoid resolution sensitivity, also 
      !     require this depth to be more than 400m, ie cumulus clouds are 
      !     deeper than 400m.  Note this depth is physically plausible and 
      !     close to the two grid-level requirement in the original L38 
      !     implementation.
      !-----------------------------------------------------------------------
      ! 4A code  test ntml-nlcl >= 2 (5A code requires 3 cloud levels)

        IF ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
          .and. zhpar(i,j)-z_lcl(i,j) >= 400.                           &
          .and. nlcl(i,j)  <   MBL-1 ) THEN

      !-----------------------------------------------------------------------
      !     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
      !     For stratocumulus top of mixed layer = zh
      !     For cumulus top of mixed layer = ZLCL
      !     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
      !     above 3000m indicates convection.
      !     Diagnosis is done by comparing gradients
      !-----------------------------------------------------------------------

          IF (zhpar(i,j) >= 3000.0) THEN
            cumulus(i,j) = .TRUE.
            !nlcl is not permitted to be less than 2 if Cu is diagnosed.
            nlcl(i,j)    = max(2, nlcl(i,j))
            nlcl_c(ii)   = nlcl(i,j)
            z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
            z_lcl_c(ii)  = z_lcl(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

          ELSE IF (nlcl(i,j) > k_plume(ii)) THEN


      ! Current test is against a height of ~<2.5km
      ! This could be replaced by a scale height if a suitable method
      ! for determining a sensible height was possible from profile/cumulus
      ! depth information available in this routine
      ! Assume moisture gradient tests stay the same whether specific or 
      ! mixing ratio though imply slightly different moisture gradients.

            IF (ntml(i,j)  >   kcucheck(ii)                             &
                  .AND. nlcl(i,j)  <=  kcucheck(ii)-2) THEN

              grad_cld =  ABS( QW(ii,kcucheck(ii)) -                    &
                                       QW(ii,nlcl_c(ii)) ) /            &
               ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            ELSE
              grad_cld =  ABS( QW(ii,ntml_c(ii)) -                      &
                                       QW(ii,nlcl_c(ii)) ) /            &
                   ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            END IF

            grad_sub =  ABS( QW(ii,nlcl_c(ii)) -                        &
                                     QW(ii,k_plume(ii)) ) /             &
                ( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            IF (grad_cld  >   1.10*grad_sub) THEN
      !-----------------------------------------------------------------------
      !     Not well mixed, however it is possible that the depth of a well
      !     mixed boundary layer has increased but not yet been mixed yet so
      !     test gradient from next level down.
      !     Note typical cumulus profiles are expected to have a fairly  
      !     uniform q profile from the surface to the cloud base and THEN a
      !     decreasing profile of q above this in the cloud. Typical the
      !     decreasing gradient from the cloud base to 2.5km will be the
      !     order of > 1.10 the below cloud value.
      !-----------------------------------------------------------------------

      ! test against a height ~ 2.5km

              IF (ntml_c(ii)  <=  kcucheck(ii)) THEN
              grad_cld =  ABS( QW(ii,ntml_c(ii)-1) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                ( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              END IF

              IF ( grad_cld  >   1.10*grad_sub) THEN
      !-----------------------------------------------------------------------
      !      Diagnose a cumulus layer
      !-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              END IF
            ELSE

      ! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
      ! and not yet been mixed (so could have been erroneously identIfied as
      ! well-mixed)

      ! First check using level below (recalculate grad_sub)

              IF (nlcl_c(ii) - k_plume(ii)  >=  2) THEN
                 grad_sub =  ABS( QW(ii,nlcl(i,j)-1) -                 &
                                      QW(ii,k_plume(ii)) ) /           &
               ( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 IF ( grad_cld  >   1.10*grad_sub) THEN
                   cumulus(i,j) =.TRUE.
                 END IF

              END IF

      ! If still diagnosing well-mixed, check using level above
      ! (recalculate grad_cld)

              IF (.NOT. cumulus(i,j) ) THEN

               IF (ntml_c(ii)  >   kcucheck(ii)                        &
                  .AND. nlcl_c(ii)  <=  kcucheck(ii)-2) THEN

                grad_cld =  ABS( QW(ii,kcucheck(ii)) -                 &
                                       QW(ii,nlcl_c(ii)+1) ) /         &
              ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               ELSE
                grad_cld =  ABS( QW(ii,ntml_c(ii)) -                   &
                                       QW(ii,nlcl_c(ii)+1) ) /        &
                 ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               END IF

               IF ( grad_cld  >   1.10*grad_sub) THEN
                 cumulus(i,j) =.TRUE.
               END IF

              END IF
            END IF
          END IF
        END IF
      END DO       ! ii loop 

      !-----------------------------------------------------------------------
      !      Check that a cumulus layer has not been erroneously diagnosed in   
      !      a deep cloudy region
      !      As the above checks are done on the total water rather than q it
      !      is possible the conditions can be met in areas where the level of  
      !      prognostic qcl or qcf is high. The type of mistake is only
      !      thought to occur over land.
      !-----------------------------------------------------------------------
      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        K=nlcl(i,j)

        IF ( land_mask(i,j) .AND. cumulus(i,j) .AND.                    &
                                      ntpar(i,j)  <   MBL ) THEN
          DO WHILE( K  <=  ntpar(i,j) .AND. cloud_fraction(i,j,K)       &
                                                 >=  SC_CFTOL )
            K = K + 1
          END DO
          IF (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        END IF

      END DO       ! ii loop 
      !-----------------------------------------------------------------
      !      Final checks on the LCL
      !-----------------------------------------------------------------
      DO ii=1,nunstable
        i = index_i(ii) 
        j = index_j(ii) 
        If ( cumulus(i,j) ) then
          If (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) then
            !-----------------------------------------------------------
            ! If LCL is diagnosed in the upper half of the layer set 
            ! z_lcl to the height of the upper layer interface
            ! (in code above LCL is always set to the lower interface).
            !-----------------------------------------------------------
            nlcl(i,j) = nlcl(i,j)+1
          End If

          !---------------------------------------------------------
          ! nlcl is not permitted to be less than nlcl_min
          !---------------------------------------------------------
          nlcl(i,j)    = max(nlcl_min(ii), nlcl(i,j))
          IF ( nlcl(i,j) >= ntpar(i,j)-1 ) THEN
            ! Cloud layer too shallow so rediagnose as well-mixed 
            cumulus(i,j)   = .false.
            L_shallow(i,j) = .false.
          END IF

          z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
          z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

        END IF

      END DO       ! ii loop 
      !-----------------------------------------------------------------------
      ! Original shallow diagnosis no congestus diagnosis
      !-----------------------------------------------------------------------
!CDIR NODEP
        DO ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          IF ( cumulus(i,j) ) THEN

      !-----------------------------------------------------------------------
      !       If cumulus has been diagnosed, determine whether it is shallow
      !       or deep convection
      !-----------------------------------------------------------------------
      ! Conditions for shallow convection
      !    wadv < 0.0            (descending air)
      !   top of parcel ascent < 2500. or T (top of parcel) > TM
      !   height of min buoyancy (above Bl) - height of parcel top T level
      !      <1.25(height parcel top - z lifting condensation level)
      !   t_dens_parc -t_dens_env at kshmin <0.55t_dens_parc -t_dens_env
      !   at ntpar
      !
      !   The last 2 conditions are looking for a strong inversion at the top
      !   of the shallow cumulus.
      !-----------------------------------------------------------------------
           IF ( shmin(ii) ) THEN

            IF ( w_avg(ii)  <   0.0 .AND.                              &
             (z_full_c(ii,ntpar(i,j))  <=  2500.0 .OR.                &
                 T(ii,ntpar(i,j))  >=  TM)                            &
            .AND. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar(i,j))) &
               <=  1.25*(zhpar(i,j) - z_lcl_c(ii)) .AND.              &
              Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) )    &
           THEN

              L_shallow(i,j) = .TRUE.
     ! may be problem with ntpar diagnosis for deep if wadv test sets
     ! L_shallow  false

            END IF

           END IF       ! test on shmin

     !-----------------------------------------------------------------------
     !      Set mixed layer depth to z_lcl
     !-----------------------------------------------------------------------
           zh(i,j) = z_lcl(i,j)
           ntml(i,j) = nlcl(i,j)

     !      If cumulus has been diagnosed but delthvu is negative, reset
     !      cumulus and L_shallow to FALSE but leave zh and ntml at LCL

           IF (delthvu(i,j)  <=  0.0) THEN

             cumulus(i,j)   = .false.
             L_shallow(i,j) = .false.

           END IF

        ELSE      ! not cumulus

     !-----------------------------------------------------------------------
     !      If not cumulus, reset parameters to within bl_levels
     !-----------------------------------------------------------------------
          IF (ntml(i,j)  >   MBL) THEN
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          END IF
          IF (nlcl(i,j)  >   MBL) THEN
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          END IF

        END IF        ! test on cumulus

      END DO          ! ii loop  

    ELSE IF (icvdiag == 7) THEN
    l_wtest  = .true.       ! w test

      DO ii=1, nunstable

        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
        z_surf = 0.1 * zh_c(ii)

        DO WHILE( z_full_c(ii,k_plume(ii))  <   z_surf .AND.           &
!                   ! not reached z_surf
                   SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

          k_plume(ii) = k_plume(ii) + 1

        END DO
      END DO

      DO ii=1, nunstable

        nlcl_min(ii) = 2

!-----------------------------------------------------------------------
! Convection scheme requires NLCL at least 2
! Also require ZLCL > 150m (approx nlcl=2 for G3 levels)
!-----------------------------------------------------------------------
        k=3
        DO WHILE ( z_half_c(ii,k) < 150.0 .AND. k < MBL ) 
          k=k+1
        END DO
        nlcl_min(ii) = k-1

      END DO

      DO ii=1, nunstable

         sl_plume(ii) = TL(ii,k_plume(ii))                              &
                              + gamma_dry * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                   &
                        min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                                &
                              /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

      END DO

!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/repsilon   q specific humidity
!   vapour pressure e ~ qp/(repsilon+q)   q mixing ratio

      IF (l_mixing_ratio) THEN
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
                                       P_theta_lev_c(ii,k_plume(ii))   &
                            / (repsilon+Q_c(ii,k_plume(ii)) )
          IF (vap_press  >   0.0) THEN
            T_LCL(ii) = a_bolton + b_bolton/                           &
                            (c_bolton*LOG(T(ii,k_plume(ii)))           &
                                          - LOG(vap_press) - d_bolton )
            P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          ELSE
            P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      ELSE
! expression for specific humidity
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          vap_press = Q_c(ii,k_plume(ii)) *                            &
              P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*repsilon )
          IF (vap_press  >   0.0) THEN
            T_LCL(ii) = a_bolton + b_bolton/                           &
                           (c_bolton*LOG(T(ii,k_plume(ii)))            &
                                         - LOG(vap_press) - d_bolton )
            P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                &
                     ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          ELSE
            P_LCL(ii) = pstar(i,j)
          END IF
 
        END DO

      END IF ! test on l_mixing_ratio  
!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      DO j=1,rows
        DO i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        END DO
      END DO
      DO ii=1, nunstable
        i        = index_i(ii)   
        j        = index_j(ii)   
        zh_c(ii) = zh(i,j)   
        zh2(ii)  = zh_c(ii)
      END DO
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  k = 2,wet_model_levels
        DO ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          IF ( P_LCL(ii)  <   P(i,j,k) ) THEN
! compressed copies
            nlcl_c(ii) = k-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = k-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          END IF     ! test on p_lcl
        END DO       ! ii loop
      END DO         ! k loop

!-----------------------------------------------------------------------
! 4.0 Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Note initial parcel conditions different from those used in the G-R
! mass flux convection scheme.
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

       DO  k = 1,wet_model_levels

         ! Require t_ref on all point for qsat call
         DO ii=1, nunstable
           t_ref(ii) = th_ref(ii)*exner_theta_levels_c(ii,k)
         END DO

! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)               &
                                          ,nunstable,l_mixing_ratio)

! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_env,t(1,k),p_theta_lev_c(1,k)               &
                                          ,nunstable,l_mixing_ratio)
         DO ii=1, nunstable
           IF(T_ref(ii) >  TM) THEN
             lrcp_const = lcrcp
             l_const    = lc
           ELSE
             lrcp_const = lsrcp
             l_const    = ls
           END IF
           IF(T(ii,k) >  TM) THEN
             lrcp_const_env = lcrcp
             l_const_env    = lc
           ELSE
             lrcp_const_env = lsrcp
             l_const_env    = ls
           END IF

           ! dqsat/dT - same whether q specific humidity or mixing ratio

           dq_sat_env = repsilon*l_const_env*qsat_lev(ii)/(R*T_ref(ii)**2)

           dq_sat_par = repsilon*l_const*qsat_env(ii)/(R*T(ii,k)**2)

           q_liq_parc = max( 0.0, ( qw_plume(ii) - qsat_lev(ii)              &
                        -dq_sat_par*( sl_plume(ii)                           &
                                       -gamma_dry*z_full_c(ii,k)-T_ref(ii) ) &
                                  ) / (1.0+lrcp_const*dq_sat_par) )

           q_liq_env  = max( 0.0, ( qw(ii,k) - qsat_env(ii)                 &
                  -dq_sat_env*( TL(ii,k)               - T(ii,k) )          &
                                  ) / (1.0+lrcp_const_env*dq_sat_env) )

           ! add on the difference in the environment's ql as calculated by the
           ! UM cloud scheme (using some RH_CRIT value) and what it
           ! would be If RH_CRIT=1. This THEN imitates partial condensation
           ! in the parcel.

           q_liq_parc = q_liq_parc + qcl_c(ii,k) + qcf_c(ii,k)- q_liq_env
           T_parc(ii,k)=sl_plume(ii)                                          &
                        -gamma_dry*z_full_c(ii,k)+lrcp_const*q_liq_parc

           ! May need to recalculate if T_parc is > Tm and T_ref < Tm

           IF (T_ref(ii) <= TM .AND. T_parc(ii,k) >  TM) THEN

             ! recalculate using corrected latent heats
             lrcp_const_parc = lcrcp
             q_liq_parc = max( 0.0, ( qw_plume(ii) - Qsat_lev(ii)             &
                          -dq_sat_par*( sl_plume(ii)                          &
                                     -gamma_dry*z_full_c(ii,k)-T_ref(ii) )    &
                                ) / (1.0+lrcp_const_parc*dq_sat_par) )
             q_liq_parc = q_liq_parc + qcl_c(ii,k)+ qcf_c(ii,k)- q_liq_env

             ! revised at parcel calculation

             T_parc(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,k)             &
                                        +lrcp_const_parc*q_liq_parc

           END IF

           q_vap_parc=qw_plume(ii)-q_liq_parc

           t_dens_parc(ii,k)=T_parc(ii,k)*(1.0+c_virtual*q_vap_parc-q_liq_parc)

           t_dens_env(ii,k) =T(ii,k)*                                      &
                             (1.0+c_virtual*Q_c(ii,k)-qcl_c(ii,k)-qcf_c(ii,k))

           buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

           env_svl(ii,k) = t_dens_env(ii,k)  + gamma_dry*z_full_c(ii,k)
           par_svl(ii,k) = t_dens_parc(ii,k) + gamma_dry*z_full_c(ii,k)

           IF (k >= 2) THEN

          !-------------------------------------------------------------
          ! Find vertical gradients in parcel and environment SVL
          ! (using values from level below (i.e. k-1)).
          !-------------------------------------------------------------

             dz = z_full_c(ii,k) - z_full_c(ii,k-1)

             dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
             denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

           END IF   ! test on k

          ! calculate t_ref for next level
           IF (k >  1 .and. k <   wet_model_levels-1) THEN
             z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                  &
                               /(z_full_c(ii,k)-z_full_c(ii,k-1))
             th_par = T_parc(ii,k)/exner_theta_levels_c(ii,k)
             th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

          ! Check sensible value otherwise set to previous reference value
          ! Problems can occur near top of model where calculation are nolonger 
          ! important.
             IF (th_ref(ii) < 0.0) THEN
               th_ref(ii) = th_par_km1(ii)
             END IF
             IF (th_par > 0.0) THEN
               th_par_km1(ii) = th_par
             END IF
  
           END IF

        END DO    ! ii loop
      END DO      ! level loop

      !-----------------------------------------------------------------------
      ! tests on parcel ascent
      !-----------------------------------------------------------------------
      !   Now compare plume s_VL with each model layer s_VL in turn to
      !     find the first time that plume has negative buoyancy.
      !-----------------------------------------------------------------------

      DO  k = 2,wet_model_levels

      !-----------------------------------------------------------------------
      ! Only perform tests if parcel ascent If unstable
      !-----------------------------------------------------------------------
        DO ii=1,nunstable

        !  Find level just below 2.5km (for use in Cu diagnosis)

          IF ( z_full_c(ii,k)  >   2500.0                              &
               .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = k-1

        ! freezing level

          IF (t(ii,k) <  TM.and.t(ii,k-1) >= TM) THEN  
            If (freeze_lev(ii) == 1) THEN
              freeze_lev(ii) = k    
            END IF
          END IF

      !-----------------------------------------------------------------------
      ! Set flag to true when level below is at least one level above the lcl
      ! and above the lcl transition zone
      ! Code implies ABOVE_LCL at NLCL+3 or greater.

          IF (k-1 >  nlcl_c(ii)+1                                      &
                      .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN
            above_lcl(ii)=.true.
          ELSE
            above_lcl(ii)=.false.
          END IF

      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

          IF ( .not.topprof(ii).and.k >  k_plume(ii) )THEN

            IF (buoyancy(ii,k) >  max_buoy(ii)) THEN
              max_buoy(ii) = buoyancy(ii,k)
              k_max(ii)    = k  
            END IF 

            ! Is parcel still buoyant ?

            IF ( (buoyancy(ii,k)  <=  - thv_pert(ii))                  &
            !    or reached top of model
                .OR. (k  >   wet_model_levels-1)  ) THEN

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,k)

             ! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

              IF ( delthvu_c(ii)  >   0.0) THEN
             ! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
                                    ( z_half_c(ii,k) - z_lcl_c(ii) )
              END IF                                                     
            END IF
          END IF


!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          IF (k > nlcl_c(ii) .AND. k < wet_model_levels ) THEN

            inc = g  * buoyancy(ii,k)                             &
                * (z_half_c(ii,k+1) - z_half_c(ii,k))/t_dens_env(ii,k)

            !---------------------------------------------------------- 
            ! If not reached level of neutral buoyancy (top of ascent)
            !---------------------------------------------------------- 

            IF (.NOT. topprof(ii)) THEN
                
              ! Note only calculating CIN and CAPE from ascents reaching
              ! level of neutral buoyancy. This may not always correspond 
              ! to the diagnosed top for the convection scheme.

              IF (inc <  0.0) THEN
                CIN_c(ii)  = CIN_c(ii) + inc
              ELSE
                CAPE_c(ii) = CAPE_c(ii) + inc
              END IF 
              dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

              delthvu_c(ii) = delthvu_c(ii) +                         &
                  buoyancy(ii,k)*(z_half_c(ii,k+1) - z_half_c(ii,k))  & 
                             /exner_theta_levels_c(ii,k)
 
            END IF    ! test on topprof

          END IF

       !-----------------------------------------------------------------------
        END DO   ! ii loop

      END DO     ! level loop

     !-----------------------------------------------------------------------
     ! Average vertical velocity over a layer  - required for shallow
     !   convection test.
     !-----------------------------------------------------------------------
     ! Layer from top in cloud w to value 1500km above cloud top?
     !-----------------------------------------------------------------------
      DO ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      END DO
      DO k=1,model_levels-1
        DO ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          IF (k >= k_neutral(ii).AND.                                     &
            z_full_c(ii,k) <= (z_half_c(ii,k_neutral(ii)+1)+1500.)) THEN

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)

          END IF
        END DO
      END DO
      DO ii=1,nunstable
        IF (mass(ii)  >  0.0 ) THEN
          w_avg(ii) = w_avg(ii)/mass(ii)
        endif
      END DO


     ! Expand back up to full arrays 

      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntml_c(ii) = k_neutral(ii)
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      END DO

      DO j = 1,rows
        DO i = 1, row_length
          zhpar(i,j) = zh(i,j)
        END DO
      END DO
     !-----------------------------------------------------------------------
     !     Test height derived above against lifting condensation level
     !-----------------------------------------------------------------------
     ! loop over only unstable 

      
!CDIR NODEP
      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
      !-----------------------------------------------------------------------
      !     Check lifting condensation level against height of parcel ascent.
      !     If lifting condensation level lower than parcel ascent, and is
      !     within bl_levels, then decide on type of cloudy layer.
      !     Gradient tests and cumulus parametriztion require a minimum number 
      !     of grid-levels between LCL and top of parcel ascent, otherwise 
      !     define as stratocumulus.  To avoid resolution sensitivity, also 
      !     require this depth to be more than 400m, ie cumulus clouds are 
      !     deeper than 400m.  Note this depth is physically plausible and 
      !     close to the two grid-level requirement in the original L38 
      !     implementation.
      !-----------------------------------------------------------------------
      ! 4A code  test ntml-nlcl >= 2 (5A code requires 3 cloud levels)

        IF ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
          .and. zhpar(i,j)-z_lcl(i,j) >= 400.                           &
          .and. nlcl(i,j)  <   MBL-1 ) THEN

      !-----------------------------------------------------------------------
      !     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
      !     For stratocumulus top of mixed layer = zh
      !     For cumulus top of mixed layer = ZLCL
      !     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
      !     above 3000m indicates convection.
      !     Diagnosis is done by comparing gradients
      !-----------------------------------------------------------------------

          IF (zhpar(i,j) >= 3000.0) THEN
            cumulus(i,j) = .TRUE.
            !nlcl is not permitted to be less than 2 if Cu is diagnosed.
            nlcl(i,j)    = max(2, nlcl(i,j))
            nlcl_c(ii)   = nlcl(i,j)
            z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
            z_lcl_c(ii)  = z_lcl(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

          ELSE IF (nlcl(i,j) > k_plume(ii)) THEN


      ! Current test is against a height of ~<2.5km
      ! This could be replaced by a scale height if a suitable method
      ! for determining a sensible height was possible from profile/cumulus
      ! depth information available in this routine
      ! Assume moisture gradient tests stay the same whether specific or 
      ! mixing ratio though imply slightly different moisture gradients.

            IF (ntml(i,j)  >   kcucheck(ii)                             &
                  .AND. nlcl(i,j)  <=  kcucheck(ii)-2) THEN

              grad_cld =  ABS( QW(ii,kcucheck(ii)) -                    &
                                       QW(ii,nlcl_c(ii)) ) /            &
               ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            ELSE
              grad_cld =  ABS( QW(ii,ntml_c(ii)) -                      &
                                       QW(ii,nlcl_c(ii)) ) /            &
                   ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            END IF

            grad_sub =  ABS( QW(ii,nlcl_c(ii)) -                        &
                                     QW(ii,k_plume(ii)) ) /             &
                ( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            IF (grad_cld  >   1.10*grad_sub) THEN
      !-----------------------------------------------------------------------
      !     Not well mixed, however it is possible that the depth of a well
      !     mixed boundary layer has increased but not yet been mixed yet so
      !     test gradient from next level down.
      !     Note typical cumulus profiles are expected to have a fairly  
      !     uniform q profile from the surface to the cloud base and THEN a
      !     decreasing profile of q above this in the cloud. Typical the
      !     decreasing gradient from the cloud base to 2.5km will be the
      !     order of > 1.10 the below cloud value.
      !-----------------------------------------------------------------------

      ! test against a height ~ 2.5km

              IF (ntml_c(ii)  <=  kcucheck(ii)) THEN
              grad_cld =  ABS( QW(ii,ntml_c(ii)-1) -                    &
                                        QW(ii,nlcl_c(ii)) ) /           &
                ( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              END IF

              IF ( grad_cld  >   1.10*grad_sub) THEN
      !-----------------------------------------------------------------------
      !      Diagnose a cumulus layer
      !-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              END IF
            ELSE

      ! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
      ! and not yet been mixed (so could have been erroneously identIfied as
      ! well-mixed)

      ! First check using level below (recalculate grad_sub)

              IF (nlcl_c(ii) - k_plume(ii)  >=  2) THEN
                 grad_sub =  ABS( QW(ii,nlcl(i,j)-1) -                 &
                                      QW(ii,k_plume(ii)) ) /           &
               ( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 IF ( grad_cld  >   1.10*grad_sub) THEN
                   cumulus(i,j) =.TRUE.
                 END IF

              END IF

      ! If still diagnosing well-mixed, check using level above
      ! (recalculate grad_cld)

              IF (.NOT. cumulus(i,j) ) THEN

                IF (ntml_c(ii)  >   kcucheck(ii)                       &
                  .AND. nlcl_c(ii)  <=  kcucheck(ii)-2) THEN

                  grad_cld =  ABS( QW(ii,kcucheck(ii)) -               &
                               QW(ii,nlcl_c(ii)+1) ) /                 &
                              ( z_full_c(ii,kcucheck(ii)) -            &
                                z_full_c(ii,nlcl_c(ii)+1) )
                ELSE
                  grad_cld =  ABS( QW(ii,ntml_c(ii)) -                 &
                               QW(ii,nlcl_c(ii)+1) ) /                 &
                              ( z_full_c(ii,ntml_c(ii)) -              &
                                z_full_c(ii,nlcl_c(ii)+1) )
                END IF

                IF ( grad_cld  >   1.10*grad_sub) THEN
                  cumulus(i,j) =.TRUE.
                END IF

              END IF ! NOT cumulus
            END IF   ! grad_cld test
          END IF     ! nlcl > k_plume
        END IF       ! cloudy boundary layer
      END DO         ! ii loop 

      !-----------------------------------------------------------------------
      !      Check that a cumulus layer has not been erroneously diagnosed in   
      !      a deep cloudy region
      !      As the above checks are done on the total water rather than q it
      !      is possible the conditions can be met in areas where the level of  
      !      prognostic qcl or qcf is high. The type of mistake is only
      !      thought to occur over land.
      !-----------------------------------------------------------------------
      DO ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        k=nlcl(i,j)

        IF ( land_mask(i,j) .AND. cumulus(i,j) .AND.                    &
                                      ntpar(i,j)  <   MBL ) THEN
          DO WHILE( k  <=  ntpar(i,j) .AND. cloud_fraction(i,j,k)       &
                                                 >=  SC_CFTOL )
            k = k + 1
          END DO
          IF (k  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        END IF

      END DO       ! ii loop 
      !-----------------------------------------------------------------
      !      Final checks on the LCL
      !-----------------------------------------------------------------
      DO ii=1,nunstable
        i = index_i(ii) 
        j = index_j(ii) 
        IF ( cumulus(i,j) ) THEN
          IF (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) THEN
            !-----------------------------------------------------------
            ! If LCL is diagnosed in the upper half of the layer set 
            ! z_lcl to the height of the upper layer interface
            ! (in code above LCL is always set to the lower interface).
            !-----------------------------------------------------------
            nlcl(i,j) = nlcl(i,j)+1
          END IF

          !---------------------------------------------------------
          ! nlcl is not permitted to be less than nlcl_min
          !---------------------------------------------------------
          nlcl(i,j)    = max(nlcl_min(ii), nlcl(i,j))
          IF ( nlcl(i,j) >= ntpar(i,j)-1 ) THEN
            ! Cloud layer too shallow so rediagnose as well-mixed 
            cumulus(i,j)   = .false.
            L_shallow(i,j) = .false.
          END IF

          z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
          z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))

        END IF

      END DO       ! ii loop 
      !-----------------------------------------------------------------------
      ! Original shallow diagnosis no congestus diagnosis
      !-----------------------------------------------------------------------
!CDIR NODEP
        DO ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          IF ( cumulus(i,j) ) THEN

      !-----------------------------------------------------------------------
      !       If cumulus has been diagnosed, determine whether it is shallow
      !       or deep convection
      !-----------------------------------------------------------------------
      ! Conditions for shallow convection
      !    wadv < 0.0            (descending air)
      !   top of parcel ascent < 2500. or T (top of parcel) > TM
      !-----------------------------------------------------------------------
            IF ( ( z_full_c(ii,ntpar(i,j))  <=  2500.0 ) .OR.         &
                 ( T(ii,ntpar(i,j))  >=  TM ) ) THEN
              IF (l_wtest) THEN
                ! Only shallow if there is significant descent
                ! 0.02m/s seems reasonable but a bit arbitrary
                ! may be better to set via namelist.
                If (w_avg(ii) < 0.02) THEN
                  L_shallow(i,j) = .TRUE.
                END IF
              ELSE
                L_shallow(i,j) = .TRUE.
              END IF  
     ! may be problem with ntpar diagnosis for deep if wadv test sets
     ! L_shallow  false

            END IF

     !-----------------------------------------------------------------------
     !      Set mixed layer depth to z_lcl
     !-----------------------------------------------------------------------
           zh(i,j) = z_lcl(i,j)
           ntml(i,j) = nlcl(i,j)

     !      If cumulus has been diagnosed but delthvu is negative, reset
     !      cumulus and L_shallow to FALSE but leave zh and ntml at LCL

           IF (delthvu(i,j)  <=  0.0) THEN

             cumulus(i,j)   = .false.
             L_shallow(i,j) = .false.

           END IF

        ELSE      ! not cumulus

     !-----------------------------------------------------------------------
     !      If not cumulus, reset parameters to within bl_levels
     !-----------------------------------------------------------------------
          IF (ntml(i,j)  >   MBL) THEN
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          END IF
          IF (nlcl(i,j)  >   MBL) THEN
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          END IF

        END IF        ! test on cumulus

      END DO          ! ii loop  
!=======================================================================
! End of option 7: simplified undilute diagnosis.
!=======================================================================


!=======================================================================
! Option 2. dilute parcel - deep entrainment rates ? 
!=======================================================================
! Calculation of LCL unchanged
! Parcel ascent now dilute
! 
!=======================================================================

    ELSE IF (icvdiag >= 2 .AND. icvdiag < 4) THEN

        l_keep_water  = .false. ! water loading not kept, reduced if >1g/kg
        l_wtest  = .true.       ! w test
     
!-----------------------------------------------------------------------
! Set entrainment fractions
!-----------------------------------------------------------------------
      If (icvdiag == 2) THEN

!  0.55/z entrainment rates eg Jakob & Siebesma   e*dz

        DO k= 1,qdims%k_end-1
          Do ii=1, nunstable
            i = index_i(ii)   
            j = index_j(ii)   
            entrain_fraction(ii,k) = 0.55*                             &
                       (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))     &
                       /z_full_c(ii,k)
          END DO
        END DO 

      Else If (icvdiag ==3) THEN

! 1/z entrainment rates - similar to many CRM results for deep and 
!                         shallow convection.

        DO k= 1,qdims%k_end-1
          Do ii=1, nunstable
            i = index_i(ii)   
            j = index_j(ii)   
            entrain_fraction(ii,k) = 1.0 *                             &
     &                 (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))     &
     &                 /z_full_c(ii,k)

          END DO
        END DO 

      END IF

      k=qdims%k_end
        Do ii=1, nunstable
            entrain_fraction(ii,k) = 0.0
        END DO

!-----------------------------------------------------------------------
! Work out initial parcel properties and LCL
!-----------------------------------------------------------------------
      Do ii=1, nunstable

        k_plume(ii) = 1
!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
        z_surf = 0.1 * zh_c(ii)

        Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
     &              SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

          k_plume(ii) = k_plume(ii) + 1

        END DO
      END DO

      Do ii=1, nunstable
         sl_plume(ii) = TL(ii,k_plume(ii))                             &
     &                        + gamma_dry * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                  &
     &                 min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                               &
     &                        /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

! 2nd set for undilute calculation
         th_ref_dil(ii)    = th_ref(ii)
         th_par_km_dil(ii) = th_ref_dil(ii)

      END DO

!-----------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/repsilon   q specific humidity
!   vapour pressure e ~ qp/(repsilon+q)   q mixing ratio

      If (l_mixing_ratio) THEN      ! expression for mixing ratio

        Do ii=1, nunstable
          vap_press = 0.01*q_c(ii,k_plume(ii)) *                       &
     &                                 P_theta_lev_c(ii,k_plume(ii))   &
     &                      / (repsilon+q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) THEN
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          END IF
        END DO

      Else                           ! expression for specific humidity

        Do ii=1, nunstable
          vap_press = q_c(ii,k_plume(ii)) *                            &
     &                P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*repsilon )

          If (vap_press  >   0.0) THEN
            T_LCL(ii) = a_bolton + b_bolton/                           &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

            P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                &
                      ( T_LCL(ii) / T(ii,k_plume(ii)) )**recip_kappa
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          END IF
        END DO

      END IF ! test on l_mixing_ratio  

!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      Do j=1,rows
        Do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        END DO
      END DO
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
      END DO


!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  k = 2,qdims%k_end
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If ( P_LCL(ii)  <   P(i,j,K) ) THEN
! compressed copies
            nlcl_c(ii) = K-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = K-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          END IF
        END DO       !  ii loop 
      END DO         !  k  loop

!-----------------------------------------------------------------------
! Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted up
! Dilute parcel ascent - mix in environmental air above lifting 
! condensation level
! 
!  Initial testing want tparc_undilute and t_parc (dilute)
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

      DO  k = 1,qdims%k_end

! Require t_ref on all point for qsat call
        Do ii=1, nunstable
          t_ref(ii)     = th_ref(ii)*exner_theta_levels_c(ii,k)
          t_ref_dil(ii) = th_ref_dil(ii)*exner_theta_levels_c(ii,k)
        END DO

! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)                &
                                             ,nunstable,l_mixing_ratio)

! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_env,t(1,k),p_theta_lev_c(1,k)                &
                                             ,nunstable,l_mixing_ratio)
! Second set
! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_lev_dil,t_ref_dil,p_theta_lev_c(1,k)        &
                                             ,nunstable,l_mixing_ratio)

        Do ii=1, nunstable

! Undilute parcel calculation required as a reference

          If(T_ref(ii) >  TM) THEN
            lrcp_const = lcrcp
            l_const    = lc
          Else
            lrcp_const = lsrcp
            l_const    = ls
          END IF

          If(T(ii,k) >  TM) THEN
            lrcp_const_env = lcrcp
            l_const_env    = lc
          Else
            lrcp_const_env = lsrcp
            l_const_env    = ls
          END IF

          dq_sat_env = repsilon*l_const_env*qsat_env(ii)/(R*T(ii,k)**2)
          dq_sat_par = repsilon*l_const*qsat_lev(ii)/(R*T_ref(ii)**2)

          q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii)              &
                                 -dq_sat_par*( sl_plume(ii)                 &
                                       -gamma_dry*z_full_c(ii,K)-T_ref(ii) )&
                                       ) / (1.0+lrcp_const*dq_sat_par) )

          q_liq_env  = MAX( 0.0, ( qw(ii,K) - qsat_env(ii)                  &
                     -dq_sat_env*( TL(ii,K)               - T(ii,k) )       &
                                       ) / (1.0+Lrcp_const*dq_sat_env) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
          ql_parc(ii,k) = q_liq_parc + qcl_c(ii,k)                          &
                                         + qcf_c(ii,k)- q_liq_env
          T_PARC(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,K)                &
                                        +lrcp_const*q_liq_parc

! May need to recalculate if T_parc is > Tm and T_ref < Tm

          If (T_ref(ii) <= TM.and.T_parc(ii,k) >  TM) THEN

! recalculate using corrected latent heats
            lrcp_const_parc = lcrcp

            q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii)           &
                              -dq_sat_par*( sl_plume(ii)                   &
                                     -gamma_dry*z_full_c(ii,K)-T_ref(ii) ) &
                                    ) / (1.0+lrcp_const_parc*dq_sat_par) )

            ql_parc(ii,k) = q_liq_parc + qcl_c(ii,k)                       &
                                         + qcf_c(ii,k)- q_liq_env
! revised at parcel calculation

            T_PARC(ii,k)=sl_plume(ii)-gamma_dry*z_full_c(ii,K)             &
                                        +lrcp_const_parc*ql_parc(ii,k)

          END IF

          q_vap_parc=qw_plume(ii)-ql_parc(ii,k)
!          
          t_dens_parc(ii,k)=T_PARC(ii,k)*                              &
                          (1.0+c_virtual*q_vap_parc-ql_parc(ii,k))


! calculate t_ref for next level
          IF (k >  1 .AND. k <   qdims%k_end-1) THEN
            z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                   &
                               /(z_full_c(ii,k)-z_full_c(ii,K-1))
            th_par = t_parc(ii,k)/exner_theta_levels_c(ii,k)
            th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

! Check sensible value otherwise set to previous reference value
! Problems can occur near top of model where calculation are nolonger 
! important.
            If (th_ref(ii) < 0.0) THEN
              th_ref(ii) = th_par_km1(ii)
            END IF
            If (th_par > 0.0) THEN   
              th_par_km1(ii) = th_par
            END IF
          END IF

! dilute parcel same as undilute parcel

          If (k <= nlcl_c(ii)) THEN

            sl_parc(ii,k) = sl_plume(ii)
            qw_parc(ii,k) = qw_plume(ii)
            t_parc_dil(ii,k)  = t_parc(ii,k)
            ql_parc_dil(ii,k) = ql_parc(ii,k)         
            t_dens_parc_dil(ii,k)= t_dens_parc(ii,k)
            th_ref_dil(ii) = th_ref(ii)
            th_par_km_dil(ii) = th_par_km1(ii)

          Else

! Dilute  parcel ascents now required

            If(t_ref_dil(ii) >  TM) THEN
              lrcp_const = lcrcp
              l_const    = lc
            Else
              lrcp_const = lsrcp
              l_const    = ls
            END IF

!-----------------------------------------------------------------------
! Dilute parcel
!-----------------------------------------------------------------------
! Mix in entrain_fraction from environmental air from level below and 
! raise this to current level.
! Assume mix in fraction of mass from environment.
! Estimate parcel properties after mixing air from environment with 
! parcel. Temperature given approximately by average

            temp_parc = (t_parc_dil(ii,k-1)                            &
                                + entrain_fraction(ii,k)*t(ii,k-1))    &
                         /(1.+entrain_fraction(ii,k))

     
            qw_parc(ii,k) = (qw_parc(ii,k-1) +                         &
                               entrain_fraction(ii,k)*qw(ii,k-1))      & 
                          /(1.+entrain_fraction(ii,k)) 

            qcl_parc = (ql_parc_dil(ii,k-1)   +                        &    
                               entrain_fraction(ii,k)*qcl_c(ii,k-1))   & 
                          /(1.+entrain_fraction(ii,k)) 

            qcf_parc = (0.0     +                                      &    
                               entrain_fraction(ii,k)*qcf_c(ii,k-1))   & 
                          /(1.+entrain_fraction(ii,k)) 

! All condensed water either ice or liquid based on t_ref ?
            sl_parc(ii,k) = temp_parc - lrcp_const*(qcl_parc+qcf_parc) &
                                 +gamma_dry*z_full_c(ii,k-1) 

            dq_sat_par_dil = repsilon*l_const*qsat_lev_dil(ii)          &
                                                  /(R*t_ref_dil(ii)**2)

            q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)         &
                            -dq_sat_par_dil*( sl_parc(ii,k)                   &
                                    -gamma_dry*z_full_c(ii,K)-t_ref_dil(ii))  &
                                     ) / (1.0+lrcp_const*dq_sat_par_dil) )

!
! add on the dIfference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
            ql_parc_dil(ii,k) = q_liq_parc + qcl_c(ii,k)                    &
                                           + qcf_c(ii,k)- q_liq_env

            t_parc_dil(ii,k)=sl_parc(ii,k)-gamma_dry*z_full_c(ii,K)         &
                                        +lrcp_const*ql_parc_dil(ii,k)

! May need to recalculate if T_parc is > Tm and T_ref < Tm

            If (t_ref_dil(ii) <= TM.and.t_parc_dil(ii,k) >  TM) THEN

! recalculate using corrected latent heats
               lrcp_const_parc = lcrcp

               q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)      &
                               -dq_sat_par_dil*(sl_parc(ii,k)                 &
                                 -gamma_dry*z_full_c(ii,K)-t_ref_dil(ii))     &
                                   ) / (1.0+lrcp_const_parc*dq_sat_par_dil) )

               ql_parc_dil(ii,k) = q_liq_parc + qcl_c(ii,k)                  &
                                    + qcf_c(ii,k)- q_liq_env

! revised at parcel calculation

              t_parc_dil(ii,k)=sl_parc(ii,k)-gamma_dry*z_full_c(ii,K)       &
                                     +lrcp_const_parc*ql_parc_dil(ii,k)

            END IF   ! test on t_ref

            q_vap_parc=qw_parc(ii,k)-ql_parc_dil(ii,k)

            If (.not.l_keep_water) THEN
!  water removed from parcel after condesation
              If (ql_parc_dil(ii,k).gt.0.001) THEN
                ql_parc_dil(ii,k) = 0.001
                qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
              END IF
            END IF
            t_dens_parc_dil(ii,k)=t_parc_dil(ii,k)*                    &
                          (1.0+c_virtual*q_vap_parc-ql_parc_dil(ii,k))

! calculate dilute t_ref for next level
            IF (k >  1 .AND. k <   qdims%k_end-1) THEN
              z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                 &
                               /(z_full_c(ii,k)-z_full_c(ii,K-1))

              th_par = t_parc_dil(ii,k)/exner_theta_levels_c(ii,k)
              th_ref_dil(ii) = th_par*(1.+z_pr) - th_par_km_dil(ii)*z_pr
! Check new reference sensible
              If (th_ref_dil(ii) < 0.0) THEN
                th_ref_dil(ii) = th_par_km_dil(ii)
              END IF
              If (th_par > 0.0) THEN   
                th_par_km_dil(ii) = th_par
              END IF
            END IF      ! k level test

          END IF   ! test on LCL
 
          t_dens_env(ii,k)=T(ii,K)*                                    &
                        (1.0+c_virtual*Q_c(ii,K)-qcl_c(ii,k)-qcf_c(ii,k))

          buoyancy(ii,k)     = t_dens_parc(ii,k)     - t_dens_env(ii,k)
          buoyancy_dil(ii,k) = t_dens_parc_dil(ii,k) - t_dens_env(ii,k)

          env_svl(ii,k) = t_dens_env(ii,k)      + gamma_dry*z_full_c(ii,K)
          par_svl(ii,k) = t_dens_parc_dil(ii,k) + gamma_dry*z_full_c(ii,K)

          If (k >= 2) THEN

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------

            dz = z_full_c(ii,K) - z_full_c(ii,K-1)

            dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
            denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

          END IF   ! test on k

        END DO    ! ii loop
      END DO      ! level loop

!-----------------------------------------------------------------------
! tests on parcel buoyancy
!-----------------------------------------------------------------------
!   Now compare plume s_VL with each model layer s_VL in turn to
!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------

      DO  k = 2,qdims%k_end

        Do ii=1,nunstable
!-----------------------------------------------------------------------
! Only perform tests if parcel ascent If unstable
!-----------------------------------------------------------------------

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                             &
               .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1

! freezing level

          If (t(ii,k) <  TM.and.t(ii,k-1) >= TM) THEN  
            If (freeze_lev(ii) == 1) THEN
              freeze_lev(ii) = k    
            END IF
          END IF

!-----------------------------------------------------------------------
! No flag for above_lcl required. Reduce thv_pert by a factor dependent
! on height relative to LCL.

          If (k-1 >  nlcl_c(ii)+1                                        &
                         .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN

! decrease thv_pert by exp(-(z-zlcl)/2000.)

            factor = exp( (z_lcl_c(ii)-z_full_c(ii,k))*1.e-3)

! set to zero if z-zlcl >1000?

            If ((z_full_c(ii,k)-z_lcl_c(ii)).gt. 1000.) THEN
              factor =0.0
            END IF 

          Else
            factor= 1.0  
          END IF

!-----------------------------------------------------------------------
! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
!-----------------------------------------------------------------------
! Not reached LNB continue testing

          If ( .not.topprof(ii).and.k >  k_plume(ii) )THEN
            If (buoyancy_dil(ii,k) >  max_buoy(ii)) THEN
              max_buoy(ii) = buoyancy_dil(ii,k)
              k_max(ii)    = k  
            END IF 

! Is parcel still buoyant ?

            If ( (buoyancy_dil(ii,k)  <=  - thv_pert(ii)*factor)       &
!                      or reached top of model
                 .OR. (k  >   qdims%k_end-1)  ) THEN

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,K)

! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy_dil(ii,k-1)

              If ( delthvu_c(ii)  >   0.0) THEN
! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
                                      ( z_half_c(ii,K) - z_lcl_C(ii) )
              END IF                                                     
            END IF
          END IF

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          IF ( k > nlcl_c(ii) .AND. k < qdims%k_end ) THEN

! Not reached top of ascent
            If (.not. topprof(ii)) THEN
              dtv_min(ii) = MIN( dtv_min(ii),                          &
                             buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

! undilute value
              delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)*          &
                           ( z_half_c(ii,K+1) - z_half_c(ii,K) )       &
                              /exner_theta_levels_c(ii,k)

! calculation of CIN and CAPE from profiles - undilute values 
              inc =g * buoyancy(ii,k)                                  &
                * (z_half_c(ii,K+1) - z_half_c(ii,K))/t_dens_env(ii,k)

              IF (inc <  0.0) THEN
                CIN_c(ii)  = CIN_c(ii) + inc
              ELSE        ! CAPE holds only positive part
                CAPE_c(ii) = CAPE_c(ii) + inc
              END IF


            END IF    ! test on topprof

          END IF      ! test on level k

!-----------------------------------------------------------------------

        END DO   ! ii loop
      END DO

!-----------------------------------------------------------------------
! Default parcel top properties are assumed to be those when the
! ascent reaches the level of neutral buoyancy. These may not be those
! required in the case of shallow convection.
! Shallow convection requires the possible identifcation of an inversion
! at the top of the ascent. This may not be detected by the LNB test.
! The gradient tests are designed to detect the shallow top.
!-----------------------------------------------------------------------

      Do ii=1,nunstable
        ntml_c(ii) = k_neutral(ii)
      END DO
!-----------------------------------------------------------------------
!!     Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
! Expand back up to full arrays 

      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      END DO
! set zhpar(i,j)
      Do j = 1,rows
        Do i = 1, row_length
          zhpar(i,j) = zh(i,j)
        END DO
      END DO

!-----------------------------------------------------------------------
! Average vertical velocity over a layer  - required for shallow
!   convection test.
!-----------------------------------------------------------------------
! Layer from top in cloud w to value 1500km above cloud top?
!-----------------------------------------------------------------------
      If (l_wtest) THEN
      Do ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      END DO
      DO k=1,tdims%k_end-1
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= ntml_c(ii).and.                                     &
             z_full_c(ii,k) <= (z_half_c(ii,ntml_c(ii)+1)+1500.)) THEN

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)
          END IF
        END DO
      END DO
      Do ii=1,nunstable
        If (mass(ii)  >  0.0 ) THEN
          w_avg(ii) = w_avg(ii)/mass(ii)
        END IF
      END DO

      END IF     ! test on w
!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
!CDIR NODEP
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any If LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! 4A scheme  ntpar-ntml >=2  (must be 2 cloud levels, 5A code requires 3)

        If ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
                              .and. nlcl(i,j)  >   k_plume(ii)          &
                                    .and. nlcl(i,j)  <   MBL-1 ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     is above boundary layer indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (zhpar(I,j) >= 3000.0) THEN
            cumulus(i,j) = .TRUE.
          Else

! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine

            If (ntml(i,j)  >   kcucheck(ii)                             &
                   .and. nlcl(i,j)  <=  kcucheck(ii)-2) THEN

              grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl_c(ii)) ) &
                /( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            Else
              grad_cld = ABS( QW(ii,ntml_c(ii)) - QW(ii,nlcl_c(ii)) )   &
                /( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            END IF

            grad_sub   =  ABS( QW(ii,nlcl_c(ii)) - QW(ii,k_plume(ii)) ) &
                 /( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            If (grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and then a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml_c(ii)  <=  kcucheck(ii)) THEN
              grad_cld = ABS( QW(ii,ntml_c(ii)-1) - QW(ii,nlcl_c(ii)) ) &
                 /( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              END IF

              If ( grad_cld  >   1.10*grad_sub) THEN
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              END IF

            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl_c(ii) - k_plume(ii)  >=  2) THEN

                 grad_sub = ABS( QW(ii,nlcl(i,j)-1) - QW(ii,k_plume(ii)) ) &
                 /( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 If ( grad_cld  >   1.10*grad_sub) THEN
                   cumulus(i,j) =.TRUE.
                 END IF

              END IF

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) THEN

               If (ntml_c(ii)  >   kcucheck(ii)                             &
                   .and. nlcl_c(ii)  <=  kcucheck(ii)-2) THEN

                grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl_c(ii)+1) ) &
                 /( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               Else
                grad_cld = ABS( QW(ii,ntml_c(ii)) - QW(ii,nlcl_c(ii)+1) )   &
                  /( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               END IF

               If ( grad_cld  >   1.10*grad_sub) THEN
                 cumulus(i,j) =.TRUE.
               END IF

              END IF   ! not cumulus
            END IF     ! test on cloud gradient
          END IF       ! test on cloud top height
        END IF         ! tests on nlcl
      END DO           ! ii loop 

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are done on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land.
!-----------------------------------------------------------------------
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        K=nlcl(i,j)

        If ( Land_MASK(i,j) .and. cumulus(i,j) .and.                    &
                                       ntpar(i,j)  <   MBL ) THEN
          Do while( K  <=  ntpar(i,j) .and. cloud_fraction(i,j,K)       &
                                                  >=  SC_CFTOL )
            K = K + 1
          END DO
          If (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        END IF
      END DO       ! ii loop 

!-----------------------------------------------------------------------
! Original shallow diagnosis no congestus diagnosis
!-----------------------------------------------------------------------
!CDIR NODEP
        Do ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          If ( cumulus(i,j) ) THEN

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
! Conditions for shallow convection 
! 
!   top of parcel ascent < 2500. or T (top of parcel) > TM
!-----------------------------------------------------------------------

            If ( z_full_c(ii,ntpar(i,j))  <=  2500.0 .OR.               &
                   T(ii,ntpar(i,j))  >=  TM ) THEN

              If (l_wtest) THEN
! Only shallow if descending air
                If (w_avg(ii).lt.0.0) THEN
                  L_shallow(i,j) = .TRUE.
                END IF
                  
              Else
                L_shallow(i,j) = .TRUE.
              END IF
            END IF

!-----------------------------------------------------------------------
!      Set mixed layer depth to z_lcl
!-----------------------------------------------------------------------
          If (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) THEN
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set z_lcl to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
            nlcl(i,j)    = nlcl(i,j)+1
            z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
            z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))
          END IF
          zh(i,j)   = z_lcl(i,j)
          ntml(i,j) = nlcl(i,j)


!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to FALSE but leave zh and NTML at LCL
!      Need undilute CAPE to be > 0

          If (delthvu(i,j)<=  0.0 ) THEN

            cumulus(i,j)   = .false.
            L_shallow(i,j) = .false.

          END IF

        Else      ! not cumulus

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within bl_levels
!-----------------------------------------------------------------------
          If (ntml(i,j)  >   MBL) THEN
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          END IF
          If (nlcl(i,j)  >   MBL) THEN
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          END IF

        END IF        ! test on cumulus

      END DO          ! ii loop  


!=======================================================================
! Option ?. New diagnosis - yet to be written  
!=======================================================================

      Else 
        error = 2     ! fatal return code
        ! Statement assumes icvdiag is in the range that format statement
        ! covers. 
        write(cmessage,'(a40,i6)')                                   &
              'Convection diagnosis option not allowed ',icvdiag
        CALL ereport(routinename,error,cmessage)

!=======================================================================
!    End of choice of diagnosis code
!=======================================================================
      END IF

!=======================================================================
! Extra calculation required if conv_diag called more than once per
! model physics timestep
!=======================================================================
      IF (l_extra_call) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j)
        DO j=1, rows
          DO i=1,row_length
            IF (fb_surf(i,j) > 0.0) THEN
              wstar(i,j) = (zh(i,j)*fb_surf(i,j))**(1.0/3.0)
              wthvs(i,j) = fb_surf(i,j)*theta(i,j,1)                        &
                                     *exner_theta_levels(i,j,1)/g
            ELSE
              wstar(i,j) = 0.0
              wthvs(i,j) = 0.0
            END IF
! BL overruled first sweep diagnosis of cumulus for a good reason
! So prevent subsequent sweeps from diagnosing convection
            IF (no_cumulus(i,j)) THEN
              cumulus(i,j)   = .false.
              L_shallow(i,j) = .false.

! Copy back orignal ntml value
              ntml(i,j) = ntml_copy(i,j)

            END IF
          END DO
        END DO
!$OMP END PARALLEL DO

      END IF        ! test on l_extra_call

!=======================================================================
!-----------------------------------------------------------------------
 9999  CONTINUE  ! Branch for error exit.

       IF (lhook) CALL dr_hook('CONV_DIAG_4A',zhook_out,zhook_handle)
       RETURN
      END SUBROUTINE CONV_DIAG_4A
    END MODULE conv_diag_4a_mod
