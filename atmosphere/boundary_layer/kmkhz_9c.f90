
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE KMKHZ  --------------------------------------------------
!
!  Purpose: To calculate the non-local turbulent mixing
!           coefficients KM and KH
!
!  Programming standard: UMDP3 vn8.4
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE kmkhz (                                                      &
! IN levels/switches
 bl_levels,lq_mix_bl,BL_diag,nSCMDpkgs,L_SCMDiags,                      &
! IN fields
 p,rho_tq,rho_uv,rho_dry_tq,t,q,qcl,qcf,cf,qw,tl,dzl,rdz,z_tq,z_uv,     &
 rad_hr,micro_tends,                                                    &
 bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,               &
 v_s,fb_surf,rhostar_gb,ntpar,nlcl,zh_prev,                             &
 zhpar,z_lcl,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,                     &
! INOUT fields
 ftl,fqw,zh,cumulus,ntml,w,etadot,t1_sd,q1_sd,                          &
! OUT fields
 rhokm,rhokh,rhokm_top,rhokh_top,zhsc,                                  &
 unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,                        &
 ntdsc,nbdsc,f_ngstress, grad_t_adj, grad_q_adj,                        &
 rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb,                   &
 tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc,                            &
 kent, we_lim, t_frac_tr, zrzi_tr,                                      &
 kent_dsc, we_lim_dsc, t_frac_dsc_tr, zrzi_dsc_tr                       &
 )

  USE atm_fields_bounds_mod, ONLY: pdims, tdims, rkmdims
  USE timestep_mod, ONLY: timestep
  USE missing_data_mod, ONLY: rmdi
  USE level_heights_mod, ONLY: eta_theta_levels 
  USE water_constants_mod, ONLY: lc, lf, tm
  USE atmos_constants_mod, ONLY: cp, r, repsilon, c_virtual
  USE cv_run_mod, ONLY: l_param_conv
  USE bl_diags_mod, ONLY : strnewbldiag
  USE bl_option_mod, ONLY:                                              &
      entr_enhance_by_cu, Buoyrev_feedback, subs_couple_fix, on,        &
      relax_sc_over_cu, a_grad_adj, max_t_grad, flux_grad, Locketal2000,&
      HoltBov1993, LockWhelan2006, entr_smooth_dec
  USE earth_constants_mod, ONLY: g
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! IN arguments
  INTEGER, INTENT(IN) ::                                                &
   bl_levels
                              ! IN No. of atmospheric levels for
                              !    which boundary layer fluxes are
                              !    calculated.
  LOGICAL, INTENT(IN) ::                                                &
   lq_mix_bl              ! IN True if using mixing ratios

!     Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
    nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
    L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

  REAL, INTENT(IN) ::                                                   &
   zmaxb_for_dsc,                                                       &
   zmaxt_for_dsc
                              ! IN Max heights to look for DSC cloud
                              !    base and top

  INTEGER, INTENT(IN) ::                                                &
   ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Top level of parcel ascent.
                              !    Used in convection scheme.
                              !    NOTE: CAN BE > BL_LEVELS-1
   nlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
         ! IN No. of model layers below the lifting condensation level. 
         !    lifting condensation level.

  LOGICAL, INTENT(IN) ::                                                &
   l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN Flag to indicate shallow
                              !    convection (only for A05_4A)

  REAL, INTENT(IN) ::                                                   &
   bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN A buoyancy parameter for clear air
                              !    on p,T,q-levels (full levels).
   bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN A buoyancy parameter for clear air
                              !    on p,T,q-levels (full levels).
   bqm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN A buoyancy parameter for clear air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN A buoyancy parameter for clear air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   bqm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                              ! IN A buoyancy parameter for cloudy air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                              ! IN A buoyancy parameter for cloudy air
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                              ! IN Saturated lapse rate factor
                              !    on p,T,q-levels (full levels).
   a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                              ! IN Saturated lapse rate factor
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                              ! IN Saturated lapse rate factor
                              !    on intermediate levels (half levels):
                              !    (*,K) elements are k+1/2 values.
   dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                              ! IN Partial derivative of QSAT w.r.t.
                              !    temperature.
   p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! IN P(*,K) is pressure at full level k.
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN Total water content (kg per kg air).
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN Liquid/frozen water temperature (K).
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! IN Temperature (K).
   qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN Cloud ice (kg per kg air)
   qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN Cloud liquid water
   q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! IN specific humidity
   cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN Cloud fractions for boundary levs.
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                              ! IN Z_tq(*,K) is the height of the
                              !    k-th full level above the surface.
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                              ! IN Z_uv(*,K) is the height of level
                              !       k-1/2 above the surface (m).
   dzl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN Layer depths (m).  DZL(,K) is the
                              !    distance from layer boundary K-1/2
                              !    to layer boundary K+1/2.  For K=1
                              !    the lower boundary is the surface.
   rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! IN Reciprocal of distance between
                              !    full levels (m-1).  1/RDZ(,K) is
                              !    the vertical distance from level
                              !    K-1 to level K, except that for
                              !    K=1 it is the height of the
                              !    lowest atmospheric full level.
   rho_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                              ! IN density on UV (ie. rho) levels,
                              !    used in RHOKH so dry density if
                              !    Lq_mix_bl is true
   rho_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                              ! IN density on TQ (ie. theta) levels,
                              !    used in RHOKM so wet density
   rho_dry_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              bl_levels)
                              ! IN density on TQ (ie. theta) levels,
                              !    used in non-turb flux integration
                              !    so dry density if Lq_mix_bl is true

  REAL, INTENT(IN) ::                                                   &
   v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! IN Surface friction velocity (m/s)
   fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! IN Surface buoyancy flux over density
                              !       (m^2/s^3).
   rhostar_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! IN Surface air density in kg/m3
   z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Height of lifting condensation
                              !    level.
   zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Height of top of NTPAR
                              !    NOTE: CAN BE ABOVE BL_LEVELS-1
   zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN boundary layer height (m) from
                              !    previous timestep

  REAL, INTENT(IN) ::                                                   &
   rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2,bl_levels),                                                 &
                              ! IN (LW,SW) radiative heating rates (K/s)
   micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               2, bl_levels)
                              ! IN Tendencies from microphysics
                              !    (TL, K/s; QW, kg/kg/s)

! INOUT arrays
  INTEGER, INTENT(INOUT) ::                                             &
   ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
      ! INOUT Number of model levels in the turbulently mixed layer.

  LOGICAL, INTENT(INOUT) ::                                             &
   cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! INOUT Flag for Cu in the bl

  REAL, INTENT(INOUT) ::                                                &
   zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! INOUT Boundary layer height (m).
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
   fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! INOUT "Explicit" fluxes of TL and QW
                              !       (rho*Km/s, rho*m/s)
                              !       IN:  level 1 (surface flux)
                              !       OUT: entrainment-level flux
   t1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! INOUT Standard Deviation of level 1
                              !    temperature (K).
   q1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! INOUT Standard Deviation of level 1
                              !    specific humidity (kg/kg).
   w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:bl_levels),  &
                              ! INOUT Vertical velocity (m/s)
   etadot(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          0:bl_levels)
                              ! INOUT d(ETA)/dt

! OUT arrays
  INTEGER, INTENT(OUT) ::                                               &
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! OUT Top level for turb mixing in
                              !       cloud layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! OUT Bottom level of any decoupled
                              !       turbulently mixed Sc layer
   sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! OUT Flags for whether discontinuous
   dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! OUT   inversions are diagnosed
   kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! OUT Grid-levels of SML and DSC
   kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! OUT   inversions (for tracer mixing)

  LOGICAL, INTENT(OUT) ::                                               &
   unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! OUT Flag to indicate an unstable
                              !     surface layer.
   dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! OUT Flag set if decoupled
                              !     stratocumulus layer found
   coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! OUT Flag to indicate Sc layer weakly
                              !     decoupled (implies mixing at SML
                              !     top is through K profiles rather
                              !     than entrainment parametrization)

  REAL, INTENT(OUT) ::                                                  &
   rhokm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                              ! OUT Non-local turbulent mixing
                              !     coefficient for momentum.
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels),                                                  &
                              ! OUT Non-local turbulent mixing
                              !     coefficient for scalars.
   rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
         2:bl_levels),                                                  &
                              ! OUT Top-down turbulent mixing
                              !     coefficient for momentum.
   rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
         2:bl_levels),                                                  &
                              ! OUT Top-down turbulent mixing
                              !     coefficient for scalars.
   f_ngstress(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              2:bl_levels),                                             &
                              ! OUT dimensionless function for
                              !     non-gradient stresses
   rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels),                                                  &
   rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
           2:bl_levels)
                              ! OUT f2 and fsc term shape profiles
                              !       multiplied by rho

  REAL, DIMENSION(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  bl_levels+1), INTENT(OUT) ::                          &
    ft_nt,                                                              &
                              ! OUT Non-turbulent heat and moisture
    fq_nt                 !       fluxes (rho*Km/s, rho*m/s)

  REAL, INTENT(OUT) ::                                                  &
   tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! OUT Total heat fluxes at inversions
   tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              !     (rho*Km/s)
   totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! OUT Total moisture fluxes at
   totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              !      inversions (rho*m/s)
   ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! OUT Non-turbulent heat and moisture
   fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              !      flux at the base of the DSC layer
   grad_t_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! OUT Temperature gradient adjustment
                              !     for non-local mixing in unstable
                              !     turbulent boundary layer.
   grad_q_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! OUT Humidity gradient adjustment
                              !     for non-local mixing in unstable
                              !     turbulent boundary layer.
   zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! OUT Cloud layer height (m).

      ! The following are used in tracer mixing.
      ! At 8B 3 elements were used - here just (i,j,2) is used.
  REAL, INTENT(OUT) ::                                                  &
   we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),       &
                              ! OUT rho*entrainment rate implied by
                              !     placing of subsidence
   zrzi_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                              ! OUT (z-z_base)/(z_i-z_base)
   t_frac_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                              ! OUT a fraction of the timestep
   we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),   &
                              ! OUT rho*entrainment rate implied by
                              !     placing of subsidence
   zrzi_dsc_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                              ! OUT (z-z_base)/(z_i-z_base)
   t_frac_dsc_tr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)
                              ! OUT a fraction of the timestep

!----------------------------------------------------------------------
!    Local and other symbolic constants :-





  REAL :: etar,grcp,lcrcp,lfrcp,ls,lsrcp
  PARAMETER (                                                           &
   etar=1.0/(1.0-repsilon),                                             &
                                ! Used in buoyancy parameter BETAC.
   grcp=g/cp,                                                           &
                                ! Adiabatic lapse rate.
   lcrcp=lc/cp,                                                         &
                                ! Latent heat of condensation / CP.
   lfrcp=lf/cp,                                                         &
                                ! Latent heat of fusion / CP.
   ls=lc+lf,                                                            &
                                ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                                ! Latent heat of sublimation / CP.
  )

  REAL :: a_plume,b_plume,a_ga_hb93,a_ga_lw06,max_svl_grad,dfsw_frac,   &
          sc_cftol,ct_resid,svl_coup,svl_coup_max,dec_svl_grad,fgf
  PARAMETER (                                                           &
   a_plume=0.2,                                                         &
   b_plume=3.26,                                                        &
   a_ga_hb93=7.2,                                                       &
   a_ga_lw06=10.0,                                                      &
   max_svl_grad=1.0e-3,                                                 &
                            ! maximum SVL gradient in a mixed layer
   dec_svl_grad=1.0e-3,                                                 &
                            ! SVL gradient required for weak decoupling
   sc_cftol=0.1,                                                        &
                            ! CF required for a Sc layer to be diagnosed
   ct_resid=200.,                                                       &
                            ! Parcel cloud-top residence time (in s)
   svl_coup_max=1.0,                                                    &
   svl_coup=0.5,                                                        &
                            ! Parameters controlling positioning of
                            ! surface-driven entrainment
   dfsw_frac = 0.35,                                                    &
                            ! Fraction of SW flux difference to 
                            ! contribute to net flux difference across 
                            ! cloud top
   fgf=0.0)                 ! Adiabatic gradient factor for ice

!  Define local storage.

!  (a) Workspace.

  LOGICAL, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) ::                      &
   cloud_base,                                                          &
                          ! Flag set when cloud base is reached.
   dsc_save           ! Copy of DSC needed to indicate
                          !   decoupling diagnosed in EXCF_NL

  REAL, DIMENSION(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end,bl_levels) ::               &
   qs,                                                                  &
                          ! Saturated sp humidity at pressure
                          !   and temperature of sucessive levels.
   cfl,                                                                 &
                          ! Liquid cloud fraction.
   cff,                                                                 &
                          ! Frozen cloud fraction.
   dqcldz,                                                              &
                          ! Vertical gradient of in-cloud liquid cloud
                          !   water in a well-mixed layer.
   dqcfdz,                                                              &
                          ! Vertical gradient of in-cloud frozen cloud
                          !   water in a well-mixed layer.
   sls_inc,                                                             &
                          ! SL and QW increments due to large-scale
   qls_inc,                                                             &
                          !    vertical advection (K s^-1, s^-1)
   df_over_cp,                                                          &
                          ! Radiative flux change over layer / c_P
   dflw_over_cp,                                                        &
                          ! LW radiative flux change over layer / c_P
   dfsw_over_cp,                                                        &
                          ! SW radiative flux change over layer / c_P
   svl,                                                                 &
                          ! Liquid/frozen water virtual temperature / CP
   sl,                                                                  &
                          ! TL + G*Z/CP (K)
   z_top,                                                               &
                          ! Z_TOP(*,K) is the height of
                          !   level k+1/2 above the surface.
   w_grad                 ! Gradient of w

  REAL, DIMENSION(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end,2:bl_levels) ::             &
  db_ga_dry,                                                            &
                          ! Out-of-cloud (DRY) and in-cloud buoyancy
  db_noga_dry,                                                          &
                          !   jumps used in flux integral calculation.
  db_ga_cld,                                                            &
                          !   GA terms include gradient adjustment
  db_noga_cld        !   arising from non-gradient fluxes. (m/s2)

!-----------------------------------------------------------------------
! The following fluxes, flux changes are in units of rho*Km/s for heat
! and rho*m/s for humidity
!-----------------------------------------------------------------------
  REAL ::                                                               &
   dfmic (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels,2),                                                 &
                                           ! Flux changes from microphys
   dfsubs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels,2),                                                 &
                                           !   and subsidence
   frad (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
          bl_levels+1),                                                 &
                                           ! Fluxes from net radiation,
   frad_lw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels+1),                                                &
                                           !   LW,
   frad_sw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels+1),                                                &
                                           !   SW,
   fmic (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
          bl_levels+1,2),                                               &
                                           !   microphys and subsidence;
   fsubs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
          bl_levels+1,2)!   for T, Q separately.

  REAL, DIMENSION(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end) ::                         &
   bflux_surf,                                                          &
                          ! Buoyancy flux at the surface.
   bflux_surf_sat,                                                      &
                          ! Saturated-air surface buoyancy flux
   db_top,                                                              &
                          ! Buoyancy jump at the top of the BL
   db_dsct,                                                             &
                          ! Buoyancy jump at the DSC layer top
   df_top_over_cp,                                                      &
                          ! Radiative flux change at cloud top / c_P
   df_dsct_over_cp,                                                     &
                          ! Radiative flux change at DSC top / CP
   svl_plume,                                                           &
                          ! SVL, SL and QW for a plume rising without
   sl_plume,                                                            &
                          !   dilution from level 1.
   qw_plume,                                                            &
                          !
   env_svl_km1,                                                         &
                          ! Density potential temperature layer K-1
   qcl_ic_top,                                                          &
   qcf_ic_top,                                                          &
                          ! In-cloud liquid and frozen water contents
                          !   at the top of the model layer
   bt_top,                                                              &
                          ! Buoyancy parameter at the top of the SML
   bt_dsct,                                                             &
                          !   and DSC
   btt_top,                                                             &
                          ! In-cloud buoyancy param at the top of the BL
   btt_dsct,                                                            &
                          !   and DSC
   btc_top,                                                             &
                          ! Cloud fraction weighted buoyancy parameter
   btc_dsct,                                                            &
                          !   at the top of the SML and DSC
   db_top_cld,                                                          &
                          ! In-cloud buoyancy jump at the top of the BL
   db_dsct_cld,                                                         &
                          !   and DSC
   cld_factor,                                                          &
                          ! Fraction of grid box potentially giving
   cld_factor_dsc,                                                      &
                          !   evaporative entrainment, for SML and DSC
   chi_s_top,                                                           &
                          ! Mixing fraction of just saturated mixture
   chi_s_dsct,                                                          &
                          !   at top of the SML and DSC layer
   zeta_s,                                                              &
                          ! Non-cloudy fraction of mixing layer for
                          !   surface forced entrainment term.
   zeta_r,                                                              &
                          ! Non-cloudy fraction of mixing layer for
                          !   cloud-top radiative cooling entrainment
   zeta_r_dsc,                                                          &
                          !   term in SML and DSC layers
   zc,                                                                  &
                          ! Cloud depth (not cloud fraction weighted).
   zc_dsc,                                                              &
                          !   for SML and DSC layer (m)
   z_cld,                                                               &
                          ! Cloud fraction weighted depth of cloud.
   z_cld_dsc,                                                           &
                          !   for SML and DSC layers (m)
   dscdepth,                                                            &
                          ! Depth of DSC layer (m)
   d_siems,                                                             &
   d_siems_dsc,                                                         &
                          ! Siems (1990) et al. cloud-top entrainment
                          ! instability parm for SML and DSC inversions
   br_fback,                                                            &
   br_fback_dsc,                                                        &
                          ! Weight for degree of buoyancy reversal
                          ! feedback for SML and DSC inversions
   tv1_sd,                                                              &
                          ! Standard Deviation of level 1 Tv
   ft_nt_zh,                                                            &
                          ! FT_NT at ZH
   ft_nt_zhsc,                                                          &
                          ! FT_NT at ZHSC
   fq_nt_zh,                                                            &
                          ! FQ_NT at ZH
   fq_nt_zhsc,                                                          &
                          ! FQ_NT at ZHSC
   df_inv_sml,                                                          &
                          ! Radiative flux divergences
   df_inv_dsc,                                                          &
                          !   over inversion grid-level
   cf_sml,                                                              &
                          ! cloud fraction of SML
   cf_dsc,                                                              &
                          ! cloud fraction of DSC layer
   z_cf_base,                                                           &
                          ! cloud base height from cloud scheme
   z_ctop,                                                              &
                          ! cloud top height
   dqw_sml,                                                             &
                          ! QW and SL changes across SML disc inv
   dsl_sml,                                                             &
                          !
   dqw_dsc,                                                             &
                          ! QW and SL changes across DSC disc inv
   dsl_dsc,                                                             &
                          !
   rhokh_surf_ent,                                                      &
                          ! SML surf-driven entrainment KH
   rhokh_top_ent,                                                       &
                          ! SML top-driven entrainment KH
   rhokh_dsct_ent,                                                      &
                          ! DSC top-driven entrainment KH
   zdsc_base,                                                           &
                          ! Height of base of K_top in DSC
   we_parm,                                                             &
                          ! Parametrised entrainment rates (m/s)
   we_dsc_parm,                                                         &
                          !   for surf and DSC layers
   we_rho,                                                              &
                          ! rho*entrainment rate
   we_rho_dsc,                                                          &
                          ! rho*entrainment rate for DSC
   w_ls,                                                                &
                          ! large-scale (subs) velocity
   w_ls_dsc,                                                            &
                          !   at subgrid inversion heights
   zh_np1,                                                              &
                          ! estimate of ZH at end of timestep
   zhsc_np1,                                                            &
                          ! estimate of ZHSC at end of timestep
   zh_frac,                                                             &
                          ! (ZH-ZHALF)/DZ
   zhsc_frac,                                                           &
                          ! (ZHSC-ZHALF)/DZ
   zrzi,                                                                &
                          ! (z-z_base)/(z_i-z_base)
   zrzi_dsc,                                                            &
                          ! (z-z_base)/(z_i-z_base)
   t_frac,                                                              &
   t_frac_dsc,                                                          &
                          ! Fraction of timestep inversion is above
                          !   entr.t flux-level for SML and DSC layers
   svl_diff_frac          ! 1 - svl difference between ntdsc and ntml 
                          ! divided by svl coupling threshold
  INTEGER, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) ::                      &
   k_cloud_top,                                                         &
                          ! Level number of top of b.l. cloud.
   k_cloud_dsct,                                                        &
                          ! Level number of top of dec. cloud.
   ntml_save,                                                           &
                          ! Copy of NTML
   ntml_prev,                                                           &
                          ! NTML from previous timestep
   k_plume,                                                             &
                          ! Start grid-level for surface-driven plume
   k_cbase            ! grid-level above cloud-base

  INTEGER ::                                                            &
   w_nonmono(pdims%i_start:pdims%i_end,                                 &
             pdims%j_start:pdims%j_end,bl_levels)
                          ! 0/1 flag for w being non-monotonic

! NEC vectorization
  INTEGER, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) ::                      &
   k_level,                                                             &
                          ! array to store level selection
   k_cff              ! level counter for CFF

!  (b) Scalars.

  REAL ::                                                               &
   virt_factor,                                                         &
                         ! Temporary in calculation of buoyancy
                         ! parameters.
   dsldz,                                                               &
                         ! Vertical gradient of SL in a well-mixed
                         ! layer.
   dqw,                                                                 &
                         ! Total water content change across layer
                         ! interface.
   dsl,                                                                 &
                         ! Liquid/ice static energy change across
                         ! layer interface.
   dsl_ga,                                                              &
                         ! As DSL but inc gradient adjustment
   dqw_ga,                                                              &
                         ! As DQW but inc gradient adjustment
   dqcl,                                                                &
                         ! Cloud liquid water change across layer
                         ! interface.
   dqcf,                                                                &
                         ! Cloud frozen water change across layer
                         ! interface.
   q_vap_parc,                                                          &
                         ! Vapour content of parcel
   q_liq_parc,                                                          &
                         ! Condensed water content of parcel
   q_liq_env,                                                           &
                         ! Condensed water content of environment
   t_parc,                                                              &
                         ! Temperature of parcel
   t_dens_parc,                                                         &
                         ! Density potential temperature of parcel
   t_dens_env,                                                          &
                         ! Density potential temperature of
                         ! environment
   denv_bydz,                                                           &
                         ! Gradient of density potential
                         ! temperature in environment
   dpar_bydz,                                                           &
                         ! Gradient of density potential
                         ! temperature of parcel
   rho_dz,                                                              &
                         ! rho*dz
   r_d_eta,                                                             &
                         ! 1/(eta(k+1)-eta(k))
   svl_lapse,                                                           &
                      ! Lapse rate of SVL above inversion (K/m)
   sl_lapse,                                                            &
                      ! Lapse rate of SL above inversion (K/m)
   qw_lapse,                                                            &
                      ! Lapse rate of QW above inversion (kg/kg/m)
   svl_lapse_base,                                                      &
                      ! Lapse rate of SVL above inversion (K/m)
   svl_top,                                                             &
                      ! s_VL at half level above inversion (K)
   dsvl_top,                                                            &
                      ! s_VL jump across inversion grid layer (K)
   tothf_efl,                                                           &
                      ! total heat flux at entrainment flux grid-level
   totqf_efl,                                                           &
                      ! Total QW flux at entrainment flux grid-level
   ml_tend,                                                             &
                      ! mixed layer tendency (d/dt)
   fa_tend,                                                             &
                      ! free atmospheric tendency (d/dt)
   inv_tend,                                                            &
                      ! limit on inversion grid-level tendency (d/dt)
   dflw_inv,                                                            &
                      ! temporary in LW rad divergence calculation
   dfsw_inv,                                                            &
                      ! temporary in SW rad divergence calculation
   dz_disc_min,                                                         &
                      ! smallest allowed DZ_DISC
   db_disc,                                                             &
                      ! Temporary disc inversion buoyancy jump
   w_s_ent,                                                             &
                      ! numerical (subsidence) entrainment rate
   dz_disc,                                                             &
                      ! height of ZH below Z_uv(NTML+2)
   z_surf,                                                              &
                      ! approx height of top of surface layer
   quad_a,                                                              &
                      ! term `a' in quadratic solver for DZ_DISC
   quad_bm,                                                             &
                      ! term `-b'in quadratic solver for DZ_DISC
   quad_c,                                                              &
                      ! term `c' in quadratic solver for DZ_DISC
   w_m,                                                                 &
                      ! scaling velocity for momentum
   w_h,                                                                 &
                      ! scaling velocity for heat
   wstar3,                                                              &
                      ! cube of convective velocity scale
   w_s_cubed,                                                           &
                      ! convective velocity scale
   z_cbase,                                                             &
                      ! cloud base height (m)
   zdsc_cbase,                                                          &
                      ! DSC cloud base height (m)
   cf_for_wb,                                                           &
                      ! CF for use in wb calculation for decoupling
   cfl_ml,cff_ml,                                                       &
                      ! liquid and frozen mixed layer fractions
   dfsw_top,                                                            &
                      ! SW radiative flux change assoc with cloud-top
   wb_test,                                                             &
                      ! test wb (m2/s-3)
   ratio,                                                               &
                      ! temporary ratio
   c_ws,                                                                &
                      ! Empirical constant multiplying Wstar
   pr_neut,                                                             &
                      ! Neutral Prandtl number
   cu_depth_fac,                                                        &
                      ! 0 to 1 factor related to cumulus depth
   w_curv,                                                              &
                      ! curvature of w
   w_curv_nm,                                                           &
   w_del_nm,                                                            &
                      ! terms used to find where w is non-monotonic
   svl_diff           ! svl difference between ntdsc and ntml 

  INTEGER ::                                                            &
   i,                                                                   &
                   ! Loop counter (horizontal field index).
   j,                                                                   &
                   ! Offset counter in certain I loops.
   k,                                                                   &
                   ! Loop counter (vertical level index).
   kl,                                                                  &
                   ! K
   kp2,                                                                 &
                   ! K+2
   kp,km,                                                               &
   kmax,                                                                &
                   ! level of maximum
   k_rad_lim   ! limit on levels within which to search for
                   !   the max LW radiative cooling

  LOGICAL ::                                                            &
   monotonic_inv,                                                       &
                   ! Flag that inversion grid-level properties are
                   ! monotonic, otherwise can't do subgrid
                   ! profile reconstruction (_DISC_INV flags)
   moisten
                   ! Indicator of whether inversion grid-level should
                   !   moisten this timestep (or dry)

  !Variables for cache-blocking
  INTEGER            :: jj          ! Block index

  INTEGER            :: jblock      ! Not a parameter. Its value needs
                                    ! to be flexible at runtime.

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('KMKHZ',zhook_in,zhook_handle)

!Set value of jblock. Platform dependent.
  jblock = pdims%j_end

!Start OpenMP parallel region

!$OMP  PARALLEL DEFAULT(SHARED)                                         &
!$OMP& PRIVATE (i, j, k, kl, rho_dz, kp, km, r_d_eta, dflw_inv,         &
!$OMP& dfsw_inv, kp2, dz_disc_min, dz_disc, dsvl_top, dfsw_top,         &
!$OMP& sl_lapse, qw_lapse, wstar3, c_ws, w_m, pr_neut, w_h, dsldz,      &
!$OMP& virt_factor, z_cbase, dqw, dsl, dsl_ga, dqw_ga, cf_for_wb,       &
!$OMP& zdsc_cbase, cfl_ml, cff_ml, kmax, db_disc,                       &
!$OMP& dqcl, dqcf, k_rad_lim, cu_depth_fac,                             &
!$OMP& t_parc, q_liq_parc, q_liq_env, q_vap_parc,                       &
!$OMP& t_dens_parc, t_dens_env, dpar_bydz, denv_bydz, quad_a, quad_bm,  &
!$OMP& quad_c, svl_lapse, svl_lapse_base, z_surf, monotonic_inv,        &
!$OMP& w_curv_nm,w_del_nm,w_curv,jj,svl_diff)
!-----------------------------------------------------------------------
! Index to subroutine KMKHZ8C

! 1. Set up local variables, etc
! 2. Look for decoupled cloudy mixed-layer above SML top
! 3. Diagnose a discontinuous inversion structure.
! 4. Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 5. Calculate uniform mixed-layer cloud fractions and cloud depths
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
! 7. Calculate inputs for the top of b.l. entrainment parametrization
! 8. Calculate the radiative flux change across cloud top
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
! 10.Call subroutine EXCF_NL
!      - calculates parametrized entrainment rate, K profiles and
!        non-gradient flux/stress functions
! 11.Calculate "explicit" entrainment fluxes of SL and QW.

!-----------------------------------------------------------------------
! 1. Set up local variables, etc
!-----------------------------------------------------------------------
! 1.1 Calculate Z_TOP (top of levels) and NTML from previous timestep
!-----------------------------------------------------------------------

!$OMP SINGLE
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      ntml_prev(i,j) = 1
    END DO
  END DO

  DO k = 1, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        z_top(i,j,k) = z_uv(i,j,k+1)
          !------------------------------------------------------------
          !find NTML from previous TS (for accurate gradient adjustment
          !of profiles - also note that NTML LE BL_LEVELS-1)
          !------------------------------------------------------------
        IF ( zh_prev(i,j)  >=  z_uv(i,j,k+1) ) ntml_prev(i,j)=k
      END DO
    END DO
  END DO
!$OMP END SINGLE nowait

  k = bl_levels
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      z_top(i,j,k) = z_uv(i,j,k) + dzl(i,j,k)
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 1.2 Calculate SVL: conserved buoyancy-like variable
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        sl(i,j,k)  = tl(i,j,k) + grcp * z_tq(i,j,k)
        svl(i,j,k) = sl(i,j,k) * ( 1.0 + c_virtual*qw(i,j,k) )
      END DO
    END DO
  END DO
!$OMP END DO nowait

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs(1,1,k),t(1,1,k),p(1,1,k),                          &
                  pdims%i_end*pdims%j_end,lq_mix_bl)
  END DO
!$OMP END DO

!--------------------------------------------------------------------
! 1.3 Integrate non-turbulent increments to give flux profiles:
!     FT_NT, FQ_NT  are the flux profiles from non-turbulent processes
!                  (consisting of radiative FRAD, subsidence FSUBS and
!                   microphysical FMIC fluxes)
!--------------------------------------------------------------------
! For heat, units of rho * Km/s
! For humidity, units of rho * m/s
!----------------------------------
IF (subs_couple_fix == ON) THEN

  DO k=1, bl_levels
    km = MAX( 1, k-1 )
    kp = MIN( bl_levels, k+1 )
!$OMP DO SCHEDULE(STATIC)
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        w_grad(i,j,k) = (w(i,j,k)-w(i,j,km))*rdz(i,j,k)
        w_curv_nm = w(i,j,kp)-2.0*w(i,j,k)+w(i,j,km)
        w_del_nm  = w(i,j,kp)-w(i,j,km)
        w_nonmono(i,j,k) = 0
        IF ( ABS(w_curv_nm) > ABS(w_del_nm) ) w_nonmono(i,j,k) = 1
        sls_inc(i,j,k) = 0.0
        qls_inc(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO
  END DO

!$OMP DO SCHEDULE(STATIC)
  DO jj=pdims%j_start,pdims%j_end,jblock
    DO k = 2,bl_levels-1
      DO j=jj,MIN((jj+jblock)-1,pdims%j_end)
        DO i=pdims%i_start,pdims%i_end

          IF ( etadot(i,j,k)  < - TINY(1.0) .AND.                       &
               etadot(i,j,k-1)< - TINY(1.0) ) THEN
!           !-----------------------------------------------------------
!           ! Only needed in subsidence regions
!           ! Also don't attempt coupling with dynamics if w has
!           ! significant vertical structure
!           !-----------------------------------------------------------
            w_curv = (w_grad(i,j,k+1)-w_grad(i,j,k))/dzl(i,j,k)
            IF ( ABS(w_curv) > 1.e-6 .AND. w_nonmono(i,j,k) == 1 ) THEN
               ! large curvature at a turning point
              sls_inc(i,j,k-1) = 0.0  ! need to make sure increments in
              qls_inc(i,j,k-1) = 0.0  ! level below are also set to zero
              etadot(i,j,k-1) = 0.0
              etadot(i,j,k)   = 0.0
              etadot(i,j,k+1) = 0.0
              w(i,j,k-1) = 0.0
              w(i,j,k)   = 0.0
              w(i,j,k+1) = 0.0
            ELSE
              kp = k+1
              km = kp-1
              r_d_eta = 1.0 /(eta_theta_levels(kp)-eta_theta_levels(km))
              sls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                &
     &                                * ( sl(i,j,kp) - sl(i,j,km) )
              qls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                &
     &                                * ( qw(i,j,kp) - qw(i,j,km) )
            END IF  ! safe to calculate increments
          END IF

        END DO
      END DO
    END DO
  END DO !jj
!$OMP END DO 

ELSE  ! subs_couple_fix off

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end

        IF ( etadot(i,j,k)  <=  0.0 ) THEN
          ! only needed in subsidence regions
          kp = MIN( bl_levels, k+1 )
          km = kp-1
          r_d_eta = 1.0 /( eta_theta_levels(kp) - eta_theta_levels(km) )
          sls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                    &
                                          * ( sl(i,j,kp) - sl(i,j,km) )
          qls_inc(i,j,k) = - etadot(i,j,k) * r_d_eta                    &
                                          * ( qw(i,j,kp) - qw(i,j,km) )
        ELSE
          sls_inc(i,j,k) = 0.0
          qls_inc(i,j,k) = 0.0
        END IF

      END DO
    END DO
  END DO
!$OMP END DO nowait

END IF  ! test on subs_couple_fix

!$OMP DO SCHEDULE(STATIC)
  DO k = 1,bl_levels
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end

        rho_dz = rho_dry_tq(i,j,k) * dzl(i,j,k)

        dflw_over_cp(i,j,k) = - rad_hr(i,j,1,k) * rho_dz
        dfsw_over_cp(i,j,k) = - rad_hr(i,j,2,k) * rho_dz
        df_over_cp(i,j,k)   = dflw_over_cp(i,j,k) + dfsw_over_cp(i,j,k)

        dfmic(i,j,k,1)  = - micro_tends(i,j,1,k) * rho_dz
        dfmic(i,j,k,2)  = - micro_tends(i,j,2,k) * rho_dz
        dfsubs(i,j,k,1) = - sls_inc(i,j,k)       * rho_dz
        dfsubs(i,j,k,2) = - qls_inc(i,j,k)       * rho_dz

      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
        ! Non-turbulent fluxes are defined relative to the surface
        ! so set them to zero at the surface
      frad(i,j,1)  = 0.0
      frad_lw(i,j,1) = 0.0
      frad_sw(i,j,1) = 0.0
      fsubs(i,j,1,1) = 0.0 ! for heat
      fsubs(i,j,1,2) = 0.0 ! for humidity
      fmic(i,j,1,1) = 0.0  ! for heat
      fmic(i,j,1,2) = 0.0  ! for humidity
    END DO
  END DO
!$OMP END DO

!$OMP SINGLE
  DO k = 2, bl_levels+1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        frad(i,j,k)   = frad(i,j,k-1)    + df_over_cp(i,j,k-1)
        frad_lw(i,j,k)= frad_lw(i,j,k-1) + dflw_over_cp(i,j,k-1)
        frad_sw(i,j,k)= frad_sw(i,j,k-1) + dfsw_over_cp(i,j,k-1)
        fsubs(i,j,k,1)= fsubs(i,j,k-1,1) + dfsubs(i,j,k-1,1)
        fsubs(i,j,k,2)= fsubs(i,j,k-1,2) + dfsubs(i,j,k-1,2)
        fmic(i,j,k,1) = fmic(i,j,k-1,1)  + dfmic(i,j,k-1,1)
        fmic(i,j,k,2) = fmic(i,j,k-1,2)  + dfmic(i,j,k-1,2)
      END DO
    END DO
  END DO
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels+1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        ft_nt(i,j,k) = frad(i,j,k) + fmic(i,j,k,1) + fsubs(i,j,k,1)
        fq_nt(i,j,k) =               fmic(i,j,k,2) + fsubs(i,j,k,2)
      END DO
    END DO
  END DO
!$OMP END DO nowait

!-----------------------------------------------------------------------
! 1.4 Set UNSTABLE flag and find first level above surface layer
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      unstable(i,j) = (fb_surf(i,j) >  0.0)
      k_plume(i,j)  = -1
    END DO
  END DO
!$OMP END DO

      !------------------------------------------------------------
      ! Find grid-level above top of surface layer, taken
      ! to be at a height, z_surf, given by:
      !       Z_SURF = 0.1*ZH_PREV
      ! Use ZH_prev since that will have determined the shape
      ! of the time-level n profiles.
      !------------------------------------------------------------

!$OMP SINGLE
  DO k = 1, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end

        IF ( unstable(i,j) ) THEN

          z_surf = 0.1 * zh_prev(i,j)
          IF ( z_tq(i,j,k) >= z_surf .AND. k_plume(i,j) == -1 ) THEN
                 !reached z_surf
            k_plume(i,j)=k
          END IF
          IF ( svl(i,j,k+1) >= svl(i,j,k)                               &
                  .AND. k_plume(i,j) == -1 ) THEN
                 !reached inversion
            k_plume(i,j)=k
          END IF

        END IF
      END DO
    END DO
  END DO

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF (k_plume(i,j) == -1) k_plume(i,j)=1
    END DO
  END DO
!$OMP END SINGLE nowait

!-----------------------------------------------------------------------
! 2.  Look for decoupled cloudy mixed-layer above SML top
!-----------------------------------------------------------------------
! 2.1  (IF NOT CUMULUS: starting from level 3 and below 2.5km):
!      find cloud-base above SML inversion, ie. above NTML+1,
!      then cloud-top (ie. CF < SC_CFTOL)
!      and finally check that cloud is well-mixed.
!-----------------------------------------------------------------------
!      Initialise variables


!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      cloud_base(i,j) = .FALSE.
      dsc(i,j) = .FALSE.
      coupled(i,j) = .FALSE.
      zhsc(i,j)    = 0.0
      ntdsc(i,j)   = 0
    END DO
  END DO
!$OMP END DO

!$OMP SINGLE
  DO k = 3, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end

!----------------------------------------------------------------------
!..Find cloud-base (where cloud here means CF > SC_CFTOL)
!----------------------------------------------------------------------

        IF ( .NOT. cumulus(i,j) .AND.                                   &
             z_tq(i,j,k) < zmaxb_for_dsc .AND.                          &
             k  >   ntml(i,j)+1 .AND. cf(i,j,k)  >   sc_cftol           &
                              .AND. .NOT.cloud_base(i,j)                &
!                                  not yet found cloud-base
                              .AND. .NOT.dsc(i,j) ) THEN
!                                  not yet found a Sc layer
          cloud_base(i,j) = .TRUE.
        END IF
        IF ( cloud_base(i,j) .AND. .NOT.dsc(i,j) .AND.                  &
!                  found cloud-base but not yet reached cloud-top
             cf(i,j,k+1) < sc_cftol .AND.                               &
             z_tq(i,j,k) < zmaxt_for_dsc                                &
!                  got to cloud-top below ZMAXT_FOR_DSC
           ) THEN
          cloud_base(i,j) = .FALSE.         ! reset CLOUD_BASE
            !-----------------------------------------------------------
            ! Look to see if at least top of cloud is well mixed:
            ! test SVL-gradient for top 2 pairs of levels, in case
            ! cloud top extends into the inversion.
            ! Parcel descent in Section 4.0 below will determine depth
            ! of mixed layer.
            !----------------------------------------------------------
          IF ( (svl(i,j,k)-svl(i,j,k-1))                                &
                       /(z_tq(i,j,k)-z_tq(i,j,k-1))                     &
                                              <   max_svl_grad ) THEN
            dsc(i,j) = .TRUE.
            ntdsc(i,j) = k
            zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
          ELSE IF ( (svl(i,j,k-1)-svl(i,j,k-2))                         &
                        /(z_tq(i,j,k-1)-z_tq(i,j,k-2))                  &
                                            <   max_svl_grad ) THEN
              !---------------------------------------------------------
              ! Well-mixed layer with top at k-1 or k.  Check whether
              ! there is a buoyancy inversion between levels k-1 and k
              ! in a manner similar to the surface-driven plume: compare
              ! the buoyancy gradient between levels K-1 and K for an
              ! undiluted parcel and the environment
              !---------------------------------------------------------
            sl_plume(i,j) = tl(i,j,k-1) + grcp * z_tq(i,j,k-1)
            qw_plume(i,j) = qw(i,j,k-1)
! -------------------------------------------------------------------
! calculate parcel water by linearising qsat about the environmental
! temperature.
! -------------------------------------------------------------------
            IF(t(i,j,k) >  tm) THEN
              q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - qs(i,j,k) -      &
                dqsdt(i,j,k)*                                           &
                ( sl_plume(i,j)-grcp*z_tq(i,j,k)-t(i,j,k) )             &
                                     ) *a_qs(i,j,k) )
              q_liq_env = MAX( 0.0, ( qw(i,j,k) - qs(i,j,k)             &
                          -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )        &
                                     ) *a_qs(i,j,k) )
! add on the difference in the environment's ql as calculated by the
! partial condensation scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
              q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)         &
                             - q_liq_env
              t_parc = sl_plume(i,j) - grcp * z_tq(i,j,k) +             &
                               lcrcp*q_liq_parc
            ELSE
              q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - qs(i,j,k) -      &
                dqsdt(i,j,k)*                                           &
                  ( sl_plume(i,j)-grcp*z_tq(i,j,k)-t(i,j,k) )           &
                                     ) *a_qs(i,j,k) )
              q_liq_env = MAX( 0.0, ( qw(i,j,k) - qs(i,j,k)             &
                 -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )                 &
                                     ) *a_qs(i,j,k) )
! add on difference in environment's ql between RH_CRIT and RH_CRIT=1
              q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)         &
                             - q_liq_env
              t_parc = sl_plume(i,j) - grcp * z_tq(i,j,k) +             &
                               lsrcp*q_liq_parc
            END IF
            q_vap_parc=qw_plume(i,j)-q_liq_parc

            t_dens_parc=t_parc*(1.0+c_virtual*q_vap_parc-q_liq_parc)
            t_dens_env=t(i,j,k)*                                        &
                       (1.0+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))
! find vertical gradients in parcel and environment SVL (using values
! from level below (K-1))
            env_svl_km1(i,j) = t(i,j,k-1) * ( 1.0+c_virtual*q(i,j,k-1)  &
                 -qcl(i,j,k-1)-qcf(i,j,k-1) ) + grcp*z_tq(i,j,k-1)
            dpar_bydz=(t_dens_parc+grcp*z_tq(i,j,k)-                    &
                        env_svl_km1(i,j)) /                             &
                    (z_tq(i,j,k)-z_tq(i,j,k-1))
            denv_bydz=(t_dens_env+grcp*z_tq(i,j,k)-                     &
                        env_svl_km1(i,j))/                              &
                    (z_tq(i,j,k)-z_tq(i,j,k-1))

            IF ( denv_bydz >  1.25*dpar_bydz ) THEN
                ! there is an inversion between levels K-1 and K
              IF ( k  >=  ntml(i,j)+3 ) THEN
                  ! if NTDSC EQ NTML+1 then assume we're looking
                  ! at the same inversion and so don't set DSC
                ntdsc(i,j) = k-1
                zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
                dsc(i,j) = .TRUE.
              END IF
            ELSE
                ! no inversion between levels K-1 and K, assume there
                ! is an inversion between K and K+1 because of CF change
              ntdsc(i,j) = k
              zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
              dsc(i,j) = .TRUE.
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
!$OMP END SINGLE

!-----------------------------------------------------------------------
! 2.2 If the layer to ZHPAR is a cumulus layer capped by cloud and
!       an inversion, declare this layer a decoupled cloud layer and
!       set ZHSC and NTDSC accordingly.
!-----------------------------------------------------------------------

  IF (relax_sc_over_cu == on) THEN
!
! Diagnosed simply if significant cloud fraction at ZHPAR
! below the height threshold zmaxt_for_dsc
!
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k = ntpar(i,j)
      IF ( cumulus(i,j) .AND. k < bl_levels  ) THEN
                ! cumulus layer within BL_LEVELS 
        IF ( z_tq(i,j,k) < zmaxt_for_dsc .AND.                          &
                ! cloud top below zmaxt_for_dsc
             ( MAX( cf(i,j,k-1),cf(i,j,k),cf(i,j,k+1) ) > sc_cftol )    &
                ! cloud-top sufficiently cloudy
            ) THEN
          dsc(i,j)  = .TRUE.
          zhsc(i,j) = zhpar(i,j)
          ntdsc(i,j)= ntpar(i,j)
        END IF
      END IF
    END DO
  END DO
!$OMP END DO nowait

  ELSE
!
! Original code, only diagnosed if shallow cu or not l_param_conv
!
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( (l_param_conv .AND.                                              &
              l_shallow(i,j) .AND. ntpar(i,j)  <   bl_levels )          &
                ! shallow cumulus layer within BL_LEVELS (A05_4A)
           .OR. (.NOT. l_param_conv .AND.                                   &
                cumulus(i,j) .AND. ntpar(i,j)  <   bl_levels ) ) THEN
                ! cumulus layer and inversion found
        IF ( cf(i,j,ntpar(i,j))  >   sc_cftol  .OR.                     &
             cf(i,j,ntpar(i,j)+1)  >   sc_cftol ) THEN
             ! cloudy
          dsc(i,j)  = .TRUE.
          zhsc(i,j) = zhpar(i,j)
          ntdsc(i,j)= ntpar(i,j)
        END IF
      END IF
    END DO
  END DO
!$OMP END DO nowait

 END IF  ! test on relax_sc_over_cu

!-----------------------------------------------------------------------
! 2.3 Calculate the radiative flux changes across cloud top for the
!      stratocumulus layer and thence a first guess for the top-down
!      mixing depth of this layer, DSCDEPTH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_cloud_dsct(i,j) = 0
      df_dsct_over_cp(i,j) = 0.0
      df_inv_dsc(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP SINGLE
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
          !-------------------------------------------------------------
          ! Find the layer with the greatest LW radiative flux jump in
          ! the upper half of the boundary layer and assume that this
          ! marks the top of the DSC layer.
          ! Necessary as radiation is not usually called every timestep.
          !-------------------------------------------------------------
          ! Limit the search to above the SML.
        k_rad_lim = ntml(i,j)+2

        IF ( dsc(i,j) .AND. k >= k_rad_lim .AND. k <= ntdsc(i,j)+2      &
              .AND. z_tq(i,j,k) > 0.5*zhsc(i,j)                         &
              .AND. dflw_over_cp(i,j,k) > df_dsct_over_cp(i,j) ) THEN
          k_cloud_dsct(i,j) = k
            ! Set K_CLOUD_DSCT to the level below if its DF is greater
            ! than half the maximum.  DF in level K_CLOUD_DSCT+1 is then
            ! included as DF_INV_DSC below.
          IF (dflw_over_cp(i,j,k-1)  >   0.5*dflw_over_cp(i,j,k))       &
             k_cloud_dsct(i,j) = k-1
          df_dsct_over_cp(i,j) = dflw_over_cp(i,j,k)
        END IF

      END DO
    END DO
  END DO

      !-----------------------------------------------------------------
      !  Find bottom grid-level (K_LEVEL) for cloud-top radiative flux
      !  divergence: higher of base of LW radiatively cooled layer,
      !  ZH and 0.5*ZHSC, since cooling must be in upper part of layer
      !  in order to generate turbulence.
      !-----------------------------------------------------------------

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_level(i,j) = k_cloud_dsct(i,j)
      IF ( k_cloud_dsct(i,j)  >   1 ) THEN
        k_rad_lim = ntml(i,j)+1
        k=k_cloud_dsct(i,j)-1
        kl=MAX(1,k)  ! only to avoid out-of-bounds compiler warning
        DO WHILE ( k  >   k_rad_lim                                     &
                  .AND. dflw_over_cp(i,j,kl)  >   0.0                   &
                  .AND. z_tq(i,j,kl)  >   0.5*zhsc(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=MAX(1,k)
        END DO
      END IF
    END DO
  END DO
!$OMP END SINGLE

      !-----------------------------------------------------------------
      ! Calculate LW and SW flux divergences and combine into
      ! cloud-top turbulence forcing.
      ! Need to account for radiative divergence in cloud in inversion
      ! grid-level, DF_INV_DSC. Assume DF_OVER_CP(K_cloud_dsct+2) is
      ! representative of clear-air rad divergence and so subtract this
      ! `clear-air' part from the grid-level divergence.
      !-----------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF ( k_cloud_dsct(i,j)  >   0 ) THEN
        dflw_inv = 0.0
        dfsw_inv = 0.0
        IF ( k_cloud_dsct(i,j) <  bl_levels ) THEN
          k = k_cloud_dsct(i,j)+1
          IF ( k  <   bl_levels ) THEN
            dflw_inv = dflw_over_cp(i,j,k)                              &
                       - dflw_over_cp(i,j,k+1)                          &
                              * dzl(i,j,k)/dzl(i,j,k+1)
            dfsw_inv = dfsw_over_cp(i,j,k)                              &
                       - dfsw_over_cp(i,j,k+1)                          &
                              * dzl(i,j,k)/dzl(i,j,k+1)
          ELSE
            dflw_inv = dflw_over_cp(i,j,k)
            dfsw_inv = dfsw_over_cp(i,j,k)
          END IF
          dflw_inv = MAX( dflw_inv, 0.0 )
          dfsw_inv = MIN( dfsw_inv, 0.0 )
        END IF
        df_inv_dsc(i,j) = dflw_inv + dfsw_inv

        df_dsct_over_cp(i,j) = frad_lw(i,j,k_cloud_dsct(i,j)+1)         &
                             - frad_lw(i,j,k_level(i,j))                &
                             + dflw_inv

        dfsw_top = frad_sw(i,j,k_cloud_dsct(i,j)+1)                     &
                 - frad_sw(i,j,k_level(i,j))                            &
                 + dfsw_inv

          !-----------------------------------------------------------
          ! Combine SW and LW cloud-top divergences into a net
          ! divergence by estimating SW flux divergence at a given
          ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
          ! Empirically (from LEM data) a reasonable fit is found
          ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = dfsw_frac
          !-----------------------------------------------------------
        df_dsct_over_cp(i,j) = MAX( 0.0,                                &
                      df_dsct_over_cp(i,j) + dfsw_frac * dfsw_top )
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2.4 Set NBDSC, the bottom level of the DSC layer.
!     Note that this will only be used to give an estimate of the layer
!     depth, DSCDEPTH, used to calculate the entrainment
!     rate (the dependence is only weak), and that a more accurate
!     algorithm is subsequently used to determine the depth over which
!     the top-down mixing profiles will be applied.  If DSC is FALSE,
!     DSCDEPTH = 0.  The plume descent here uses a radiative
!     perturbation to the cloud-layer SVL (use level NTDSC-1 in case
!     SVL is not yet well-mixed to NTDSC), based roughly
!     on a typical cloud-top residence time.  If the plume does not sink
!     and the cloud is decoupled from the surface (ie. above Alan
!     Grant's ZH), then it is assumed to be stable, ie. St rather than
!     Sc, and no mixing or entrainment is applied to it.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      nbdsc(i,j) = ntdsc(i,j)+1
      IF (dsc(i,j)) THEN
! The depth of the radiatively-cooled layer tends to be less than O(50m)
! and so RAD_HR will be an underestimate of the cooling tendency there.
! Compensate by multiplying by DZL/50. (~4)
! Recall that DF_OVER_CP(I,j,K) = RAD_HR * RHO_DRY_TQ * DZL
! Thus use cloud-top radiative forcing as follows:

        k = ntdsc(i,j)
        rho_dz = rho_dry_tq(i,j,k) * dzl(i,j,k)
        svl_plume(i,j)=svl(i,j,k-1)                                     &
           - ct_resid * dzl(i,j,k)*df_dsct_over_cp(i,j) / ( 50.*rho_dz )

      ELSE
        svl_plume(i,j)=0.0
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP SINGLE
  DO k = bl_levels-1, 1, -1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        IF ( k <  ntdsc(i,j) .AND. svl_plume(i,j)  <   svl(i,j,k) ) THEN
          nbdsc(i,j) = k+1     ! marks lowest level within ML
        END IF
      END DO
    END DO
  END DO
!$OMP END SINGLE

!----------------------------------------------------------------------
! 2.5 Tidy up variables associated with decoupled layer
!       NOTE that NTDSC GE 3 if non-zero
!----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
! Note that ZHSC-Z_UV(NTML+2) may = 0, so this test comes first!
      IF (cumulus(i,j) .AND. dsc(i,j))                                  &
                       nbdsc(i,j) = MAX( nbdsc(i,j), ntml(i,j)+2 )
      IF ( ntdsc(i,j)  >=  1 ) THEN
        IF ( nbdsc(i,j) <   ntdsc(i,j)+1 ) THEN
          dscdepth(i,j) =                                               &
                      z_uv(i,j,ntdsc(i,j)+1) - z_uv(i,j,nbdsc(i,j))
        ELSE
         !----------------------------------------------------------
         ! Indicates a layer of zero depth
         !----------------------------------------------------------
          IF (ntdsc(i,j) == ntpar(i,j)) THEN
            !----------------------------------------------------------
            ! Indicates a Sc layer at the top of Cu: force mixing
            ! over single layer.
            !----------------------------------------------------------
            dscdepth(i,j) = dzl(i,j,ntdsc(i,j))
          ELSE
            dsc(i,j)=.FALSE.
            ntdsc(i,j)=0
            zhsc(i,j)=0.0
            df_dsct_over_cp(i,j) = 0.0
            k_cloud_dsct(i,j) = 0
            df_inv_dsc(i,j)   = 0.0
            dscdepth(i,j) = 0.0
          END IF
        END IF
      ELSE  ! ntdsc eq 0, just to make sure!
        dscdepth(i,j)=0.0
        dsc(i,j)=.FALSE.
        zhsc(i,j)=0.0
        df_dsct_over_cp(i,j) = 0.0
        k_cloud_dsct(i,j) = 0
        df_inv_dsc(i,j)   = 0.0
      END IF
    END DO
  END DO
!$OMP END DO

!----------------------------------------------------------------------
!2.6 If decoupled cloud-layer found test to see if it is, in fact,
!  only weakly decoupled from the surface mixed-layer:
!  if SVL difference between NTML and NTDSC is less than svl_coup (in K)
!  then assume there is still some coupling.  This will mean that
!  the surface-driven entrainment term will be applied at ZHSC, no
!  subgrid inversion or entrainment will be calculated for ZH and
!  ZHSC will be the length scale used in the entrainment inputs.
!  Note that for CUMULUS "surface-driven entrainment" will be done
!  by the convection scheme.
!----------------------------------------------------------------------
IF (entr_smooth_dec == on) THEN
  !-------------------------------------------------------------
  ! entr_smooth_dec:
  ! OFF - original method
  ! ON  - taper off surface terms to zero for svl_diff between 
  !       svl_coup and svl_coup_max; also ignore cumulus diags
  !-------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      coupled(i,j)       = .FALSE.
      svl_diff_frac(i,j) = 0.0   ! Fully coupled by default
      IF ( dsc(i,j) ) THEN
        !-------------------------------------------------------------
        ! Calculate cloud to surface mixed layer SVL difference
        ! - avoid ntdsc as can be within base of inversion
        !-------------------------------------------------------------
        svl_diff           = 0.0
        IF ( ntdsc(i,j) >= 2 )                                          &
                  svl_diff = svl(i,j,ntdsc(i,j)-1) - svl(i,j,ntml(i,j))
        IF ( svl_diff  < svl_coup_max ) THEN
          coupled(i,j) = .TRUE.
          svl_diff_frac(i,j) = 1.0 - MAX( 0.0,                          &
                         (svl_diff-svl_coup)/(svl_coup_max-svl_coup) )
                         ! to give 1 for svl_diff<svl_coup and
                         ! decrease linearly to 0 at svl_coup_max
          dscdepth(i,j) = zhsc(i,j)
        END IF
      END IF  ! dsc test
    END DO
  END DO
!$OMP END DO

ELSE  ! entr_smooth_dec test

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      coupled(i,j)       = .FALSE.
      svl_diff_frac(i,j) = 0.0   ! Fully coupled by default
      IF ( dsc(i,j) .AND. .NOT.cumulus(i,j) ) THEN
        !-----------------------------------------------------------
        ! Note this IF test structure is required because if DSC is
        ! false then NTDSC = 0 and cannot be used to index SVL.
        !-----------------------------------------------------------
        IF ( svl(i,j,ntdsc(i,j)) - svl(i,j,ntml(i,j)) < svl_coup )      &
               coupled(i,j) = .TRUE.
      END IF
    END DO
  END DO
!$OMP END DO

END IF  ! entr_smooth_dec test
!-----------------------------------------------------------------------
! 3. Diagnose a discontinuous inversion structure:
!    - to this point in the code, ZH and ZHSC mark the half-level at
!      the base of the inversion
!    - now they will be interpolated into the level above assuming
!      SVL(NTML+1) is a volume average over a subgrid discontinuous
!      inversion structure
!    - discontinuous jumps of SL and QW (and thence buoyancy) can be
!      calculated and used to determine the entrainment rate
!    - parametrized grid-level fluxes at NTML,NTDSC can then be made
!      consistent with this assumed inversion structure
!-----------------------------------------------------------------------
       ! If any `problems' are encountered with this interpolation of ZH
       ! (such as ZH diagnosed at or below Z_UV(NTML+1)), then NTML
       ! is lowered a level and ZH is set fractionally below what has
       ! become Z_UV(NTML+2).  This distance is such that for a net
       ! dZH/dt of 1.E-4 m/s, ZH will be diagnosed as spending at least
       ! half the timestep in level NTML+2, leaving the growth only
       ! marginally affected.  Conversely, it allows a subsiding
       ! inversion to fall more readily.

  dz_disc_min = 0.5 * timestep * 1.e-4

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      sml_disc_inv(i,j) = 0  ! initialise flags to indicate whether a
      dsc_disc_inv(i,j) = 0  ! discontinuous inversion is diagnosed
!..First interpolate to find ZH

      k = ntml(i,j)
!..by default, keep ZH at the half-level where it was diagnosed
!..initially and use grid-level jumps

      dsl_sml(i,j) = sl(i,j,k+1) - sl(i,j,k)
      dqw_sml(i,j) = qw(i,j,k+1) - qw(i,j,k)

      IF ( .NOT.cumulus(i,j) .AND. .NOT.coupled(i,j) .AND.              &
                    k  >   1 .AND. k  <=  bl_levels-2 ) THEN

          !  Require SVL and SL to be monotonically increasing
          !  and QW to be simply monotonic
        monotonic_inv = ( svl(i,j,k+2) > svl(i,j,k+1).AND.              &
                          svl(i,j,k+1) > svl(i,j,k) )                   &
                  .AND. ( sl(i,j,k+2) > sl(i,j,k+1) .AND.               &
                          sl(i,j,k+1) > sl(i,j,k) )                     &
                  .AND. ( ( qw(i,j,k+2) > qw(i,j,k+1) .AND.             &
                            qw(i,j,k+1) > qw(i,j,k) )                   &
                      .OR.( qw(i,j,k+2) < qw(i,j,k+1) .AND.             &
                            qw(i,j,k+1) < qw(i,j,k) ) )

        IF ( monotonic_inv ) THEN

          IF ( k  <=  bl_levels-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis so best just to ignore lapse)
            svl_lapse = MAX(0.0,                                        &
                  ( svl(i,j,k+3) - svl(i,j,k+2) ) * rdz(i,j,k+3)  )
            IF ( svl_lapse  >                                           &
                  ( svl(i,j,k+2) - svl(i,j,k+1) ) * rdz(i,j,k+2) )      &
                  svl_lapse = 0.0
          ELSE
            svl_lapse = 0.0
          END IF
          IF ( k  >=  k_plume(i,j)+2 ) THEN
                ! Use mean mixed layer gradient (if resolved) to allow
                ! for stablisation by gradient-adjustment
                ! Ignore level K in case inversion is dropping
            svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,k_plume(i,j)) )/    &
                          (z_tq(i,j,k-1)-z_tq(i,j,k_plume(i,j)))
            svl_lapse_base = MAX( 0.0, svl_lapse_base )
          ELSE
            svl_lapse_base = 0.0
          END IF

          quad_a  = 0.5*( svl_lapse - svl_lapse_base )
          quad_bm = svl(i,j,k+2) - svl(i,j,k)                           &
              - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )             &
              - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +          &
                                                      dzl(i,j,k+1) )
          quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -          &
              svl_lapse_base * (                                        &
                z_uv(i,j,k+1)-z_tq(i,j,k) + 0.5*dzl(i,j,k+1) ) )

          IF ( quad_bm  >   0.0 ) THEN
            IF ( quad_c  <=  0.0) THEN
                  ! SVL extrapolated from K to K+1 is greater than
                  ! the level K+1 value - inversion needs to rise so
                  ! place it as high as possible
              dz_disc = dz_disc_min
            ELSE IF ( quad_bm*quad_bm  >=  4.0*quad_a*quad_c ) THEN
                  ! solve equation for DZ_DISC...
              IF ( quad_a  /=  0.0 ) THEN
                    !   ...quadratic if QUAD_A NE 0
                dz_disc = ( quad_bm - SQRT( quad_bm*quad_bm             &
                                         - 4.0*quad_a*quad_c )          &
                                ) / (2.0*quad_a)
              ELSE
                    !   ...linear if QUAD_A EQ 0
                dz_disc = quad_c / quad_bm
              END IF
            ELSE
              dz_disc = 99999.9  ! large dummy value
            END IF

            IF ( dz_disc  >   0.9 * dzl(i,j,k+1) ) THEN
                ! ZH diagnosed very close to or below Z_UV(K+1):
              IF ( svl(i,j,k)-svl(i,j,k-1)  >   0.0) THEN
                    ! top of ML stably stratified so lower NTML but
                    ! set ZH only fractionally (DZ_DISC_MIN)
                    ! below the top of the inversion level.
                ntml(i,j) = ntml(i,j) - 1
                k=ntml(i,j)
                dz_disc = dz_disc_min
              ELSE
                    ! top of ML well-mixed so don't lower the inversion
                    ! level but set ZH just (DZ_DISC_MIN) above the
                    ! half-level to allow the inversion to subside if
                    ! necessary.
                dz_disc = dzl(i,j,k+1) - dz_disc_min
              END IF
            END IF

          ELSE
!.. ignoring lapse rates
            dsvl_top = svl(i,j,k+2) - svl(i,j,k)
            dz_disc = dzl(i,j,k+1) *                                    &
                            (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
          END IF

          zh(i,j) = z_uv(i,j,k+2) - dz_disc
          sml_disc_inv(i,j) = 1 ! set flag to indicate disc inv found

!-----------------------------------------------------------
!..Calculate SML inversion discontinuous jumps of SL and QW
!-----------------------------------------------------------
            ! Allow for lapse rate above inversion, if known
          dz_disc = z_tq(i,j,k+2) - zh(i,j)
          sl_lapse = 0.0
          qw_lapse = 0.0
          IF ( k  <=  bl_levels-3 ) THEN
            sl_lapse = MAX( 0.0,                                        &
               ( sl(i,j,k+3) - sl(i,j,k+2) )*rdz(i,j,k+3) )
            qw_lapse = MIN( 0.0,                                        &
               ( qw(i,j,k+3) - qw(i,j,k+2) )*rdz(i,j,k+3) )
          END IF
            !-----------------
            ! First SL jump
            !-----------------
            ! Only reduce 2 level jump by at most half
          dsl_sml(i,j) = sl(i,j,k+2) - sl(i,j,k)
          dsl_sml(i,j) = dsl_sml(i,j) -                                 &
                 MIN( 0.5*dsl_sml(i,j), sl_lapse*dz_disc )
            !-----------------
            ! Next QW jump
            !-----------------
          IF ( qw(i,j,k+2) < qw(i,j,k+1) .AND.                          &
               qw(i,j,k+1) < qw(i,j,k) ) THEN
              ! QW monotonically decreasing across inversion
              ! Only allow for QW lapse rate if both it and the
              ! 2 grid-level jump are negative (expected sign)
            dqw_sml(i,j) = qw(i,j,k+2) - qw(i,j,k)
            IF ( dqw_sml(i,j) < 0.0 ) THEN
              dqw_sml(i,j) = dqw_sml(i,j) -                             &
                MAX( 0.5*dqw_sml(i,j), qw_lapse*dz_disc )
            END IF
          ELSE IF ( qw(i,j,k+2) > qw(i,j,k+1) .AND.                     &
                    qw(i,j,k+1) > qw(i,j,k) ) THEN
              ! QW monotonically increasing across inversion
              ! Suggests something unusual is going so not clear how
              ! to proceed, so currently leaving DQW as 2 level jump
            dqw_sml(i,j) = qw(i,j,k+2) - qw(i,j,k)
          END IF

        END IF  ! Monotonic inversion
      END IF ! not cumulus and not at top of bl_levels
!-----------------------------------------------------------------------
!..Second interpolate to find ZHSC
!-----------------------------------------------------------------------
      IF ( dsc(i,j) ) THEN
        k = ntdsc(i,j)
!..by default, keep ZHSC at the half-level where it was diagnosed
!..initially and use grid-level jumps
        dsl_dsc(i,j) = sl(i,j,k+1) - sl(i,j,k)
        dqw_dsc(i,j) = qw(i,j,k+1) - qw(i,j,k)
        IF ( k  <=  bl_levels-2 ) THEN

          !  Require SVL and SL to be monotonically increasing
          !  and QW to be simply monotonic
          monotonic_inv = ( svl(i,j,k+2) > svl(i,j,k+1).AND.            &
                            svl(i,j,k+1) > svl(i,j,k) )                 &
                    .AND. ( sl(i,j,k+2) > sl(i,j,k+1) .AND.             &
                            sl(i,j,k+1) > sl(i,j,k) )                   &
                    .AND. ( ( qw(i,j,k+2) > qw(i,j,k+1) .AND.           &
                              qw(i,j,k+1) > qw(i,j,k) )                 &
                        .OR.( qw(i,j,k+2) < qw(i,j,k+1) .AND.           &
                              qw(i,j,k+1) < qw(i,j,k) ) )

          IF ( monotonic_inv ) THEN

            IF ( k  <=  bl_levels-3 ) THEN
! need to test for K+1 to K+2 gradient in case profile is concave
! (would mess up the inversion diagnosis so best just to ignore)
              svl_lapse = MAX(0.0,                                      &
                    ( svl(i,j,k+3) - svl(i,j,k+2) )*rdz(i,j,k+3) )
              IF ( svl_lapse  >                                         &
                    ( svl(i,j,k+2) - svl(i,j,k+1) )*rdz(i,j,k+2) )      &
                    svl_lapse = 0.0
            ELSE
              svl_lapse = 0.0
            END IF
            IF ( k  >=  nbdsc(i,j)+2 ) THEN
                ! Use mean mixed layer gradient (if resolved) to allow
                ! for stablisation by gradient-adjustment
                ! Ignore level K in case inversion is dropping
              svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,nbdsc(i,j)) )/    &
                            (z_tq(i,j,k-1)-z_tq(i,j,nbdsc(i,j)))
              svl_lapse_base = MAX( 0.0, svl_lapse_base )
            ELSE
              svl_lapse_base = 0.0
            END IF

            quad_a  = 0.5*( svl_lapse - svl_lapse_base )
            quad_bm = svl(i,j,k+2) - svl(i,j,k)                         &
                 - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )          &
                 - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +       &
                                                        dzl(i,j,k+1) )
            quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -        &
                 svl_lapse_base * (                                     &
                  z_uv(i,j,k+1)-z_tq(i,j,k) + 0.5*dzl(i,j,k+1) ) )

            IF ( quad_bm  >   0.0 ) THEN
              IF ( quad_c  <=  0.0) THEN
                  ! SVL extrapolated from K to K+1 is greater than
                  ! the level K+1 value - inversion needs to rise
                dz_disc = dz_disc_min
              ELSE IF ( quad_bm*quad_bm  >=  4.0*quad_a*quad_c ) THEN
                  ! solve equation for DZ_DISC...
                IF ( quad_a  /=  0.0 ) THEN
                    !   ...quadratic if QUAD_A NE 0
                  dz_disc = ( quad_bm - SQRT( quad_bm*quad_bm           &
                                           - 4.0*quad_a*quad_c )        &
                                  ) / (2.0*quad_a)
                ELSE
                    !   ...linear if QUAD_A EQ 0
                  dz_disc = quad_c / quad_bm
                END IF
              ELSE
                dz_disc = 99999.9  ! large dummy value
              END IF

              IF ( dz_disc  >   0.9 * dzl(i,j,k+1) ) THEN
                IF ( ntdsc(i,j) == 2 ) THEN
                  dz_disc = dzl(i,j,k+1)
                ELSE
                ! ZHSC diagnosed very close to or below Z_UV(K+1):
                  IF ( svl(i,j,k)-svl(i,j,k-1)  >   0.0) THEN
                    ! top of ML stably stratified so lower NTDSC but
                    ! set ZHSC only fractionally (DZ_DISC_MIN)
                    ! below the top of the inversion level.
                    ntdsc(i,j) = ntdsc(i,j) - 1
                    k=ntdsc(i,j)
                    dz_disc = dz_disc_min
                    dscdepth(i,j) = dscdepth(i,j) - dzl(i,j,k+1)
                    ! Note that all but DZ_DISC_MIN of this layer will
                    ! be added back on to DSCDEPTH a few lines below
                  ELSE
                    ! top of ML well-mixed so don't lower the inversion
                    ! level but set ZHSC just (DZ_DISC_MIN) above the
                    ! half-level to allow the inversion to subside if
                    ! necessary.
                    dz_disc = dzl(i,j,k+1) - dz_disc_min
                  END IF
                END IF
              END IF

            ELSE  ! QUAD_BM le 0
!.. ignoring lapse rates
              dsvl_top = svl(i,j,k+2) - svl(i,j,k)
              dz_disc = dzl(i,j,k+1) *                                  &
                              (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
            END IF

            zhsc(i,j) = z_uv(i,j,k+2) - dz_disc
            dscdepth(i,j) = dscdepth(i,j) + zhsc(i,j) - z_uv(i,j,k+1)
            dsc_disc_inv(i,j) = 1  ! set flag to indicate disc inv found

!-----------------------------------------------------------
!..Calculate DSC inversion discontinuous jumps of SL and QW
!-----------------------------------------------------------
            ! Allow for lapse rate above inversion, if known
            dz_disc = z_tq(i,j,k+2) - zhsc(i,j)
            sl_lapse = 0.0
            qw_lapse = 0.0
            IF ( k  <=  bl_levels-3 ) THEN
              sl_lapse = MAX( 0.0,                                      &
                 ( sl(i,j,k+3) - sl(i,j,k+2) )*rdz(i,j,k+3) )
              qw_lapse = MIN( 0.0,                                      &
                 ( qw(i,j,k+3) - qw(i,j,k+2) )*rdz(i,j,k+3) )
            END IF
            !-----------------
            ! First SL jump
            !-----------------
            ! Only reduce 2 level jump by at most half
            dsl_dsc(i,j) = sl(i,j,k+2) - sl(i,j,k)
            dsl_dsc(i,j) = dsl_dsc(i,j) -                               &
                 MIN( 0.5*dsl_dsc(i,j), sl_lapse*dz_disc )
            !-----------------
            ! Next QW jump
            !-----------------
            IF ( qw(i,j,k+2) < qw(i,j,k+1) .AND.                        &
                 qw(i,j,k+1) < qw(i,j,k) ) THEN
              ! QW monotonically decreasing across inversion
              ! Only allow for QW lapse rate if both it and the
              ! 2 grid-level jump are negative (expected sign)
              dqw_dsc(i,j) = qw(i,j,k+2) - qw(i,j,k)
              IF ( dqw_dsc(i,j) < 0.0 ) THEN
                dqw_dsc(i,j) = dqw_dsc(i,j) -                           &
                  MAX( 0.5*dqw_dsc(i,j), qw_lapse*dz_disc )
              END IF
            ELSE IF ( qw(i,j,k+2) > qw(i,j,k+1) .AND.                   &
                      qw(i,j,k+1) > qw(i,j,k) ) THEN
              ! QW monotonically increasing across inversion
              ! Suggests something unusual is going so not clear how
              ! to proceed, so currently leaving DQW as 2 level jump
              dqw_dsc(i,j) = qw(i,j,k+2) - qw(i,j,k)
            END IF

          END IF  ! monotonic inversion
        END IF ! test on K LT BL_LEVELS-2
      END IF ! test on DSC
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 4.  Calculate the within-layer vertical gradients of cloud liquid
!      and frozen water
! 4.1 Calculate gradient adjustment terms
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      virt_factor = 1.0 + c_virtual*q(i,j,1) - qcl(i,j,1) - qcf(i,j,1)
      grad_t_adj(i,j) = 0.0
      grad_q_adj(i,j) = 0.0
      IF ( unstable(i,j) ) THEN
          ! Here this is an estimate of the gradient adjustment applied
          ! the previous timestep (assumes T1_SD has not changed much,
          ! which in turn assumes the surface fluxes have not)
        IF (flux_grad  ==  Locketal2000) THEN
          grad_t_adj(i,j) = MIN( max_t_grad ,                           &
                           a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
          grad_q_adj(i,j) = 0.0
        ELSE IF (flux_grad  ==  HoltBov1993) THEN
            ! Use constants from Holtslag and Boville (1993)
            ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
            ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
          wstar3 = fb_surf(i,j) * zh_prev(i,j)
          c_ws = 0.6
          w_m =( v_s(i,j)**3 + c_ws*wstar3 )**(1.0/3.0)

          grad_t_adj(i,j) = a_ga_hb93*(wstar3**(1.0/3.0))*ftl(i,j,1)    &
                            / ( rhostar_gb(i,j)*w_m*w_m*zh_prev(i,j) )
!            GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FQW(I,j,1)
!     &                       / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH_PREV(I,j) )
          grad_q_adj(i,j) = 0.0
        ELSE IF (flux_grad  ==  LockWhelan2006) THEN
            ! Use constants from LockWhelan2006
            ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
            ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
          wstar3  = fb_surf(i,j) * zh_prev(i,j)
          c_ws    = 0.42   !  = 0.75^3
          pr_neut = 0.75
          w_h = ( ( v_s(i,j)**3+c_ws*wstar3 )**(1.0/3.0) )/ pr_neut

          grad_t_adj(i,j) = a_ga_lw06 * ftl(i,j,1)                      &
                             / ( rhostar_gb(i,j)*w_h*zh_prev(i,j) )
          grad_q_adj(i,j) = a_ga_lw06 * fqw(i,j,1)                      &
                             / ( rhostar_gb(i,j)*w_h*zh_prev(i,j) )
        END IF
      END IF  ! test on UNSTABLE

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!! 4.2 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the current layer
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end

        IF (k  <=  ntml_prev(i,j)) THEN
          dsldz = -grcp + grad_t_adj(i,j)
        ELSE
          dsldz = -grcp
        END IF

        virt_factor = 1.0 + c_virtual*q(i,j,k) - qcl(i,j,k) -           &
                            qcf(i,j,k)

        dqcldz(i,j,k) = -( dsldz*dqsdt(i,j,k)                           &
                       + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) )         &
                        / ( 1.0 + lcrcp*dqsdt(i,j,k) )
        dqcfdz(i,j,k) = -( dsldz*dqsdt(i,j,k)                           &
                       + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) ) * fgf   &
                        / ( 1.0 + lsrcp*dqsdt(i,j,k) )

          ! limit calculation to greater than a small cloud fraction
        IF ( qcl(i,j,k) + qcf(i,j,k)  >   0.0                           &
             .AND. cf(i,j,k)  >   1.e-3 ) THEN
          cfl(i,j,k) = cf(i,j,k) * qcl(i,j,k) /                         &
                       ( qcl(i,j,k) + qcf(i,j,k) )
          cff(i,j,k) = cf(i,j,k) * qcf(i,j,k) /                         &
                       ( qcl(i,j,k) + qcf(i,j,k) )
        ELSE
          cfl(i,j,k) = 0.0
          cff(i,j,k) = 0.0
        END IF

      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 5.  Calculate uniform mixed-layer cloud fraction and thence
!        estimate Sc layer cloud depth (not cloud fraction weighted).
!        (If DSC=.FALSE. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
! First the SML
!---------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      cloud_base(i,j)= .FALSE.
      zc(i,j)        = 0.0
      k_cbase(i,j)   = 0
      z_cf_base(i,j) = zh(i,j)
      z_ctop(i,j)    = zh(i,j)
        ! Use a single CF for whole mixed-layer (more realistic).
        ! Include NTML+1 if a subgrid inversion has been diagnosed
      IF ( coupled(i,j) .OR. cumulus(i,j) .OR. ntml(i,j) == 1 ) THEN
        cf_sml(i,j)=0.0
      ELSE
        k = ntml(i,j)
        cf_sml(i,j) = MAX( cf(i,j,k), cf(i,j,k-1) )
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! First find cloud-base as seen by the cloud scheme
! [K_LEVEL=first level below NTML with CF<SC_CFTOL and
!  height Z_CF_BASE=half-level above]
! to use as first guess or lower limit
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_level(i,j) = ntml(i,j)
      IF ( cf_sml(i,j)  >   sc_cftol ) THEN
        k_level(i,j) = ntml(i,j)-1
        DO WHILE ( cf(i,j,MAX(k_level(i,j),1))  >   sc_cftol           &
                   .AND. k_level(i,j)  >=  2 )
          k_level(i,j) = k_level(i,j) - 1
        END DO
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( cf_sml(i,j)  >   sc_cftol ) THEN
        IF ( k_level(i,j)  ==  1 .AND.                                 &
             cf(i,j,MAX(k_level(i,j),1))  >   sc_cftol) THEN
          z_cf_base(i,j) = 0.0
        ELSE
          z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
        END IF
        zc(i,j) = z_ctop(i,j) - z_cf_base(i,j)
      END IF
    END DO
  END DO
!$OMP END DO

!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------

!$OMP SINGLE
  DO k = bl_levels, 1, -1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        IF ( .NOT.cloud_base(i,j) .AND. k  <=  ntml(i,j) .AND.          &
             cf_sml(i,j)  >   sc_cftol ) THEN
               ! within cloudy boundary layer
          IF ( k  ==  1) THEN
            cloud_base(i,j) = .TRUE.
          ELSE
            IF ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .TRUE.
          END IF
          k_cbase(i,j) = k
        END IF
      END DO
    END DO
  END DO
!$OMP END SINGLE

!Initialise K_CFF = lowest level with ice cloud

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_cff(i,j) = k_cbase(i,j)
      IF (k_cff(i,j) > 1) THEN
        DO WHILE ( cff(i,j,k_cff(i,j))  >   sc_cftol                    &
                      .AND. k_cff(i,j)  >   1 )
          k_cff(i,j) = k_cff(i,j) - 1
        END DO
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

          !--------------------------------------------------
          ! Use adiabatic qcl gradient to estimate cloud-base
          ! from in-cloud qcl in level K_CBASE
          ! If k_cbase = 0 then it hasn't been initialised
          !--------------------------------------------------

      IF ( cloud_base(i,j) .AND. k_cbase(i,j) /= 0 ) THEN
        z_cbase = z_tq(i,j,k_cbase(i,j)) -                              &
                  qcl(i,j,k_cbase(i,j)) /                               &
                  ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
        IF ( dqcfdz(i,j,k_cbase(i,j))  >   0.0 ) THEN
          z_cbase = MIN( z_cbase, z_tq(i,j,k_cbase(i,j)) -              &
                 qcf(i,j,k_cbase(i,j)) /                                &
                 ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )      &
                       )
        ELSE
              !---------------------------------------------------------
              ! No adiabatic QCF gradient so find lowest level, K_CFF,
              ! with CFF>SC_CFTOL and assume cloud-base within that leve
              !---------------------------------------------------------
          IF ( cff(i,j,k_cff(i,j))  <=  sc_cftol .AND.                  &
                        k_cff(i,j)  <   k_cbase(i,j) )                  &
               k_cff(i,j) = k_cff(i,j) + 1
                   ! will want to raise K_CFF back up one level unless
                   ! level 1 is cloudy or no sig frozen cloud at all
          z_cbase = MIN( z_cbase, z_top(i,j,k_cff(i,j)) -               &
                    dzl(i,j,k_cff(i,j))                                 &
                  * cff(i,j,k_cff(i,j))/cf(i,j,k_cff(i,j)) )
        END IF
            !------------------------------------------------------
            ! use cloud-base as seen by cloud scheme as lower limit
            ! and base of level NTML+1 as upper limit
            !------------------------------------------------------
        z_cbase = MIN( z_uv(i,j,ntml(i,j)+1),                           &
                       MAX( z_cf_base(i,j), z_cbase) )

        zc(i,j) = z_ctop(i,j) - z_cbase
      END IF

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Second DSC layer
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      cloud_base(i,j) = .FALSE.
      zc_dsc(i,j) = 0.0
      cf_dsc(i,j) = 0.0
      k_cbase(i,j) = 0
      z_cf_base(i,j) = zhsc(i,j)
      z_ctop(i,j)    = zhsc(i,j)

      IF ( dsc(i,j) ) THEN
        k = ntdsc(i,j)
        cf_dsc(i,j) = MAX( cf(i,j,k), cf(i,j,k-1) )
      END IF
    END DO
  END DO
!$OMP END DO

!-------------------------------------------------------------
! Find cloud-base as seen by cloud scheme, Z_CF_BASE,
! to use as first guess or lower limit and find cloud top.
!-------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_level(i,j) = ntdsc(i,j)
      IF ( cf_dsc(i,j)  >   sc_cftol ) THEN
          ! assume level NTDSC is cloudy so start from NTDSC-1
        k_level(i,j) = MAX( 2, ntdsc(i,j) - 1 )
        DO WHILE ( cf(i,j,k_level(i,j))  >   sc_cftol           &
                 .AND. k_level(i,j)  >=  2 )
          k_level(i,j) = k_level(i,j) - 1
        END DO
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( cf_dsc(i,j)  >   sc_cftol ) THEN
        IF ( k_level(i,j)  ==  1 .AND.                                  &
             cf(i,j,MAX(k_level(i,j),1))  >   sc_cftol) THEN
          z_cf_base(i,j) = 0.0
        ELSE
          z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
        END IF
        zc_dsc(i,j) = z_ctop(i,j) - z_cf_base(i,j)   ! first guess
      END IF
    END DO
  END DO
!$OMP END DO

!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------

!$OMP SINGLE
  DO k = bl_levels, 1, -1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        IF ( .NOT.cloud_base(i,j) .AND. k  <=  ntdsc(i,j) .AND.         &
             cf_dsc(i,j)  >   sc_cftol ) THEN
               ! within cloudy boundary layer
          IF ( k  ==  1) THEN
            cloud_base(i,j) = .TRUE.
          ELSE
            IF ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .TRUE.
          END IF
          k_cbase(i,j) = k
        END IF
      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END SINGLE

! Initialise K_CFF

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_cff(i,j) = k_cbase(i,j)
      IF (k_cff(i,j) > 1) THEN
        DO WHILE ( cff(i,j,k_cff(i,j))  >   sc_cftol                    &
                      .AND. k_cff(i,j)  >   1)
          k_cff(i,j) = k_cff(i,j) - 1
        END DO
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

        !--------------------------------------------------
        ! use adiabatic qcl gradient to estimate cloud-base
        ! from in-cloud qcl in level K_CBASE
        !--------------------------------------------------
      IF ( cloud_base(i,j) .AND. k_cbase(i,j)  /= 0 ) THEN
        z_cbase = z_tq(i,j,k_cbase(i,j)) -                              &
                  qcl(i,j,k_cbase(i,j)) /                               &
                  ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
        IF ( dqcfdz(i,j,k_cbase(i,j))  >   0.0 ) THEN
          z_cbase = MIN( z_cbase, z_tq(i,j,k_cbase(i,j)) -              &
                qcf(i,j,k_cbase(i,j)) /                                 &
                ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )       &
                       )
        ELSE
            !----------------------------------------------------------
            ! No adiabatic QCF gradient so find lowest level, K_CFF,
            ! with CFF>SC_CFTOL and assume cloud-base within that level
            !----------------------------------------------------------
          IF ( cff(i,j,k_cff(i,j))  <=  sc_cftol .AND.                  &
                        k_cff(i,j)  <   k_cbase(i,j) )                  &
               k_cff(i,j) = k_cff(i,j) + 1
                 ! will want to raise K_CFF back up one level unless
                 ! level 1 is cloudy or no sig frozen cloud at all
          z_cbase = MIN( z_cbase, z_top(i,j,k_cff(i,j)) -               &
                    dzl(i,j,k_cff(i,j))                                 &
                   * cff(i,j,k_cff(i,j))/cf(i,j,k_cff(i,j)) )
        END IF
          !------------------------------------------------------
          ! use cloud-base as seen by cloud scheme as lower limit
          ! and base of level NTDSC+1 as upper limit
          !------------------------------------------------------
        z_cbase = MIN( z_uv(i,j,ntdsc(i,j)+1),                          &
                       MAX( z_cf_base(i,j) , z_cbase) )

        zc_dsc(i,j) = z_ctop(i,j) - z_cbase
      END IF

    END DO !I
  END DO !J
!$OMP END DO

      !-----------------------------------------------------------------
      !  Layer cloud depth cannot be > the layer depth itself.
      !-----------------------------------------------------------------

!$OMP SINGLE
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      zc_dsc(i,j) = MIN( zc_dsc(i,j), dscdepth(i,j) )
    END DO
  END DO
!$OMP END SINGLE

!-----------------------------------------------------------------------
! 6. Calculate buoyancy flux factor used in the diagnosis of decoupling
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        dqw = qw(i,j,k) - qw(i,j,k-1)
        dsl = sl(i,j,k) - sl(i,j,k-1)
        IF (k  <=  ntml_prev(i,j)) THEN
          dsl_ga = dsl - grad_t_adj(i,j)/rdz(i,j,k)
          dqw_ga = dqw - grad_q_adj(i,j)/rdz(i,j,k)
        ELSE
          dsl_ga = dsl
          dqw_ga = dqw
        END IF
          !----------------------------------------------------------
          ! CF_FOR_WB is uniform `bl' CF for use within cloud layers
          !----------------------------------------------------------
        cf_for_wb = 0.0
        z_cbase = zh(i,j)-zc(i,j)
        zdsc_cbase = zhsc(i,j)-zc_dsc(i,j)
        IF ( z_tq(i,j,k)  <=  zh(i,j) .AND.                             &
             z_tq(i,j,k)  >=  z_cbase) cf_for_wb = cf_sml(i,j)
        IF ( z_tq(i,j,k)  <=  zhsc(i,j) .AND.                           &
             z_tq(i,j,k)  >=  zdsc_cbase) cf_for_wb = cf_dsc(i,j)
          !----------------------------------------------------------
          ! WB = -K_SURF*(DB/DZ - gamma_buoy) - K_TOP*DB/DZ
          ! This is integrated in EXCF_NL, iterating the K profiles.
          ! Here the relevant integrated DB/DZ factors are calculated
          !----------------------------------------------------------
        db_ga_dry(i,j,k) = - g *                                        &
                 ( btm(i,j,k-1)*dsl_ga + bqm(i,j,k-1)*dqw_ga )
        db_noga_dry(i,j,k)  = - g *                                     &
                 ( btm(i,j,k-1)*dsl + bqm(i,j,k-1)*dqw )
        db_ga_cld(i,j,k) = - g *                                        &
                 ( btm_cld(i,j,k-1)*dsl_ga + bqm_cld(i,j,k-1)*dqw_ga )
        db_noga_cld(i,j,k)  = - g *                                     &
                 ( btm_cld(i,j,k-1)*dsl + bqm_cld(i,j,k-1)*dqw )
          !-------------------------------------------------------
          ! Weight cloud layer factors with cloud fraction
          !-------------------------------------------------------
        db_ga_cld(i,j,k) = db_ga_dry(i,j,k)*(1.0-cf_for_wb) +           &
                           db_ga_cld(i,j,k)*cf_for_wb
        db_noga_cld(i,j,k)  = db_noga_dry(i,j,k)*(1.0-cf_for_wb) +      &
                              db_noga_cld(i,j,k)*cf_for_wb
      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 7. Calculate inputs for the top of b.l. entrainment parametrization
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      zeta_r_dsc(i,j) = 0.0
      chi_s_dsct(i,j) = 0.0
      d_siems_dsc(i,j) = 0.0
      br_fback_dsc(i,j)= 0.0
      cld_factor_dsc(i,j) = 0.0
      bt_dsct(i,j) = 0.0
      btt_dsct(i,j) = 0.0
      btc_dsct(i,j) = 0.0
      db_dsct(i,j) = 0.0
      db_dsct_cld(i,j) = 0.0
      chi_s_top(i,j) = 0.0
      d_siems(i,j) = 0.0
      br_fback(i,j)= 0.0
      cld_factor(i,j) = 0.0
      bt_top(i,j) = 0.0
      btt_top(i,j) = 0.0
      btc_top(i,j) = 0.0
      db_top(i,j) = 0.0
      db_top_cld(i,j) = 0.0    ! default required if COUPLED
      z_cld(i,j) = 0.0
      z_cld_dsc(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 7.1 Calculate surface buoyancy flux
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

        ! use mixed-layer average of buoyancy parameters
      bflux_surf(i,j) = 0.5 * g * (                                     &
           (btm(i,j,1)+btm(i,j,ntml(i,j)))*ftl(i,j,1) +                 &
           (bqm(i,j,1)+bqm(i,j,ntml(i,j)))*fqw(i,j,1) )

      IF ( bflux_surf(i,j)  >   0.0 ) THEN
        bflux_surf_sat(i,j) = 0.5 * g * (                               &
           (btm_cld(i,j,1)+btm_cld(i,j,ntml(i,j)))*ftl(i,j,1) +         &
           (bqm_cld(i,j,1)+bqm_cld(i,j,ntml(i,j)))*fqw(i,j,1) )
        IF ( coupled(i,j) ) bflux_surf_sat(i,j) = 0.5 * g * (           &
           (btm_cld(i,j,1)+btm_cld(i,j,ntdsc(i,j)))*ftl(i,j,1) +        &
           (bqm_cld(i,j,1)+bqm_cld(i,j,ntdsc(i,j)))*fqw(i,j,1) )
      ELSE
        bflux_surf_sat(i,j) = 0.0
      END IF

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 7.2 Calculation of cloud fraction weighted thickness of
!     cloud in the SML and DSC layer (or to the surface if COUPLED)
!-----------------------------------------------------------------------

!$OMP SINGLE
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        IF ( k  <=  ntml(i,j)+1 ) THEN
          z_cld(i,j) = z_cld(i,j) +                                     &
                     cf(i,j,k) * 0.5 * dzl(i,j,k) +                     &
                      MIN( cfl(i,j,k) * 0.5 * dzl(i,j,k) ,              &
                              qcl(i,j,k) / dqcldz(i,j,k) )
          IF ( dqcfdz(i,j,k)  >   0.0) THEN
            z_cld(i,j) = z_cld(i,j) +                                   &
                      MIN( cff(i,j,k) * 0.5 * dzl(i,j,k) ,              &
                              qcf(i,j,k) / dqcfdz(i,j,k) )
          ELSE
            z_cld(i,j) = z_cld(i,j) + cff(i,j,k) * 0.5 * dzl(i,j,k)
          END IF
        END IF

        IF ( dsc(i,j) .AND. k <= ntdsc(i,j)+1 .AND.                     &
             ( coupled(i,j) .OR.                                        &
                   z_top(i,j,k) >= zhsc(i,j)-zc_dsc(i,j) ) ) THEN
          z_cld_dsc(i,j) = z_cld_dsc(i,j) +                             &
                     cf(i,j,k) * 0.5 * dzl(i,j,k) +                     &
                      MIN( cfl(i,j,k) * 0.5 * dzl(i,j,k) ,              &
                              qcl(i,j,k) / dqcldz(i,j,k) )
          IF ( dqcfdz(i,j,k)  >   0.0) THEN
            z_cld_dsc(i,j) = z_cld_dsc(i,j) +                           &
                      MIN( cff(i,j,k) * 0.5 * dzl(i,j,k) ,              &
                              qcf(i,j,k) / dqcfdz(i,j,k) )
          ELSE
            z_cld_dsc(i,j) = z_cld_dsc(i,j) +                           &
                                    cff(i,j,k) * 0.5 * dzl(i,j,k)
          END IF
        END IF
      END DO
    END DO
  END DO
!$OMP END SINGLE

!-----------------------------------------------------------------------
! 7.3 Calculate the buoyancy jumps across the inversions
!----------------------------------------------------------------------
!..Where appropriate, replace grid-level jumps with jumps calculated
!..assuming a discontinuous subgrid inversion structure.
!----------------------------------------------------------------------

!$OMP SINGLE
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
        !--------------------------
        ! First the SML
        !--------------------------
      qcl_ic_top(i,j) = 0.0
      qcf_ic_top(i,j) = 0.0
      cfl_ml = 0.0
      cff_ml = 0.0

      km = ntml(i,j)
      IF ( sml_disc_inv(i,j) == 1 ) THEN
          !-----------------------------------------------------
          ! Extrapolate water contents from level with max CF,
          ! out of NTML and NTML-1 (ie near top of SML),
          ! to the top of the mixed layer
          !-----------------------------------------------------
        IF (cf_sml(i,j)  >   0.0) THEN
          kmax = km
          IF (cf(i,j,km-1) > cf(i,j,km)) kmax = km-1

          cfl_ml = cf_sml(i,j)*cfl(i,j,kmax)                            &
                                 /(cfl(i,j,kmax)+cff(i,j,kmax)+1.e-10)
          cff_ml = cf_sml(i,j)*cff(i,j,kmax)                            &
                                 /(cfl(i,j,kmax)+cff(i,j,kmax)+1.e-10)

          IF (cfl_ml > 0.01) qcl_ic_top(i,j) = qcl(i,j,kmax)/cfl_ml     &
                           + ( zh(i,j)-z_tq(i,j,km) )*dqcldz(i,j,km)
          IF (cff_ml > 0.01) qcf_ic_top(i,j) = qcf(i,j,kmax)/cff_ml     &
                           + ( zh(i,j)-z_tq(i,j,km) )*dqcfdz(i,j,km)
        END IF

        dqw = dqw_sml(i,j)
        dsl = dsl_sml(i,j)
          ! ignore any cloud above the inversion
        dqcl = - cfl_ml*qcl_ic_top(i,j)
        dqcf = - cff_ml*qcf_ic_top(i,j)

        db_disc = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +             &
                 (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +        &
                 (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf   )

        IF ( db_disc > 0.03  ) THEN
            ! Diagnosed inversion statically stable and at least ~1K
          db_top(i,j) = db_disc
        ELSE
            ! Diagnosed inversion statically UNstable
            ! Reset flag to use entrainment K (rather than fluxes)
          sml_disc_inv(i,j) = 0
          zh(i,j) = z_uv(i,j,ntml(i,j)+1)
        END IF
      END IF  ! disc inversion diagnosed

      IF ( sml_disc_inv(i,j) == 0 ) THEN
           ! Calculate using simple grid-level differences
        kp = km+1
        dqw = qw(i,j,kp) - qw(i,j,km)
        dsl = sl(i,j,kp) - sl(i,j,km)
        qcl_ic_top(i,j) = qcl(i,j,km)
        qcf_ic_top(i,j) = qcf(i,j,km)
        dqcl = qcl(i,j,kp) - qcl_ic_top(i,j)
        dqcf = qcf(i,j,kp) - qcf_ic_top(i,j)
        db_top(i,j) = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +         &
                  (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +       &
                  (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf )
      END IF  ! no disc inversion diagnosed

      db_top_cld(i,j) = g * ( btm_cld(i,j,km)*dsl                       &
                            + bqm_cld(i,j,km)*dqw )
      chi_s_top(i,j) = -qcl_ic_top(i,j) /                               &
                          (a_qsm(i,j,km)*dqw - a_dqsdtm(i,j,km)*dsl)
      chi_s_top(i,j) = MAX( 0.0, MIN( chi_s_top(i,j), 1.) )

      IF ( db_top(i,j)  <   0.003 ) THEN
          ! Diagnosed inversion statically unstable:
          ! ensure DB>0 so that entrainment is non-zero and
          ! instability can be removed.
        db_top(i,j) = 0.003
        db_top_cld(i,j) = 0.0  ! set buoyancy reversal
        chi_s_top(i,j) = 0.0   ! term to zero
      END IF

    END DO
  END DO
!$OMP END SINGLE

        !--------------------------
        ! Then the DSC layer
        !--------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( dsc(i,j) ) THEN

        qcl_ic_top(i,j) = 0.0
        qcf_ic_top(i,j) = 0.0
        cfl_ml = 0.0
        cff_ml = 0.0

        km = ntdsc(i,j)
        IF ( dsc_disc_inv(i,j) == 1 ) THEN
           !-----------------------------------------------------
           ! Extrapolate water contents from level with max CF,
           ! out of NTDSC and NTDSC-1 (ie near top of DSC),
           ! to the top of the mixed layer
           !-----------------------------------------------------
          IF (cf_dsc(i,j) > 0.0) THEN
            kmax = km
            IF (cf(i,j,km-1) > cf(i,j,km)) kmax = km-1

            cfl_ml = cf_dsc(i,j)*cfl(i,j,kmax)                          &
                                  /(cfl(i,j,kmax)+cff(i,j,kmax)+1.e-10)
            cff_ml = cf_dsc(i,j)*cff(i,j,kmax)                          &
                                  /(cfl(i,j,kmax)+cff(i,j,kmax)+1.e-10)
            IF (cfl_ml > 0.01) qcl_ic_top(i,j) = qcl(i,j,kmax)/cfl_ml   &
                          + ( zhsc(i,j)-z_tq(i,j,km) )*dqcldz(i,j,km)
            IF (cff_ml > 0.01) qcf_ic_top(i,j) = qcf(i,j,kmax)/cff_ml   &
                          + ( zhsc(i,j)-z_tq(i,j,km) )*dqcfdz(i,j,km)

          END IF

          dqw = dqw_dsc(i,j)
          dsl = dsl_dsc(i,j)
           ! ignore any cloud above the inversion
          dqcl = - cfl_ml*qcl_ic_top(i,j)
          dqcf = - cff_ml*qcf_ic_top(i,j)

          db_disc = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +           &
                  (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +       &
                  (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf   )

          IF ( db_disc > 0.03  ) THEN
              ! Diagnosed inversion statically stable
            db_dsct(i,j) = db_disc
          ELSE
              ! Diagnosed inversion statically UNstable
              ! Reset flag to use entrainment K (rather than fluxes)
            dsc_disc_inv(i,j) = 0
            zhsc(i,j) = z_uv(i,j,ntdsc(i,j)+1)
          END IF
        END IF  ! disc inversion diagnosed

        IF ( dsc_disc_inv(i,j) == 0 ) THEN
           ! Calculate using simple grid-level differences
          kp = km+1
          dqw = qw(i,j,kp) - qw(i,j,km)
          dsl = sl(i,j,kp) - sl(i,j,km)
          qcl_ic_top(i,j) = qcl(i,j,km)
          qcf_ic_top(i,j) = qcf(i,j,km)
          dqcl = qcl(i,j,kp) - qcl_ic_top(i,j)
          dqcf = qcf(i,j,kp) - qcf_ic_top(i,j)
          db_dsct(i,j) = g * ( btm(i,j,km)*dsl + bqm(i,j,km)*dqw +      &
                    (lcrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcl +     &
                    (lsrcp*btm(i,j,km) - etar*bqm(i,j,km)) * dqcf )
        END IF  ! no disc inversion diagnosed

        db_dsct_cld(i,j) = g * ( btm_cld(i,j,km)*dsl                    &
                               + bqm_cld(i,j,km)*dqw )
        chi_s_dsct(i,j) = -qcl_ic_top(i,j) /                            &
                             (a_qsm(i,j,km)*dqw - a_dqsdtm(i,j,km)*dsl)
        chi_s_dsct(i,j) = MAX( 0.0, MIN( chi_s_dsct(i,j), 1.) )

        IF ( db_dsct(i,j) < 0.003 ) THEN
           ! Diagnosed inversion statically unstable:
           ! ensure DB>0 so that entrainment is non-zero and
           ! instability can be removed.
          db_dsct(i,j) = 0.003
          db_dsct_cld(i,j) = 0.0  ! set buoyancy reversal
          chi_s_dsct(i,j) = 0.0   ! term to zero
        END IF
      END IF  ! test on DSC

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 7.3 Calculation of other SML and DSC inputs to entr param.
!     If COUPLED then SML are not used as no "entrainment" is then
!     applied at ZH.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
          !------------------------------------------------------
          ! Calculation of SML inputs.
          !------------------------------------------------------
      k = ntml(i,j)
      kp2=MIN(k+1+sml_disc_inv(i,j),bl_levels)
      cld_factor(i,j) = MAX( 0.0 , cf_sml(i,j)-cf(i,j,kp2) )
      bt_top(i,j)  = g * btm(i,j,k)
      btt_top(i,j) = g * btm_cld(i,j,k)
      btc_top(i,j) = btt_top(i,j)
          !---------------------------------------------------
          ! Calculation of DSC inputs
          !---------------------------------------------------
      IF (dsc(i,j)) THEN
        k = ntdsc(i,j)
        kp2=MIN(k+1+dsc_disc_inv(i,j),bl_levels)
        cld_factor_dsc(i,j) = MAX( 0.0 , cf_dsc(i,j)-cf(i,j,kp2) )
        bt_dsct(i,j)  = g * btm(i,j,k)
        btt_dsct(i,j) = g * btm_cld(i,j,k)
        btc_dsct(i,j) = btt_dsct(i,j)
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 7.4 Next those terms which depend on the presence of buoyancy reversal
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      z_cld(i,j) = MIN( z_cld(i,j), zh(i,j) )
      z_cld_dsc(i,j) = MIN( z_cld_dsc(i,j), zhsc(i,j) )
        !---------------------------------------------------------------
        ! First the surface mixed layer.
        !---------------------------------------------------------------
      IF ( coupled(i,j) ) THEN
        zeta_s(i,j) = 1.0 - z_cld_dsc(i,j) / zhsc(i,j)
        zeta_r(i,j) = 1.0 - zc_dsc(i,j) / zhsc(i,j)
      ELSE
        zeta_s(i,j) = 1.0 - z_cld(i,j) / zh(i,j)
        zeta_r(i,j) = 1.0 - zc(i,j) / zh(i,j)
      END IF

      IF (db_top_cld(i,j)  >=  0.0) THEN
          !--------------------------------------------------
          ! i.e. no buoyancy reversal (or default if COUPLED)
          !--------------------------------------------------
        db_top_cld(i,j) = 0.0
        d_siems(i,j) = 0.0
        br_fback(i,j)= 0.0
      ELSE
          !----------------------------
          ! IF (DB_TOP_CLD(I,j)  <   0.0)
          ! i.e. buoyancy reversal
          !----------------------------
        db_top_cld(i,j) = -db_top_cld(i,j) * cld_factor(i,j)
        d_siems(i,j) = MAX( 0.0,                                        &
             chi_s_top(i,j) * db_top_cld(i,j) / (db_top(i,j)+1.e-14) )
          ! Linear feedback dependence for D<0.1
        br_fback(i,j)= MIN( 1.0, 10.0*d_siems(i,j) )
        zeta_r(i,j)  = zeta_r(i,j) + (1.0-zeta_r(i,j))*br_fback(i,j)
      END IF
        !---------------------------------------------------------------
        ! Now the decoupled Sc layer (DSC).
        !---------------------------------------------------------------
      IF (dsc(i,j)) THEN
        IF ( coupled(i,j) ) THEN
          zeta_r_dsc(i,j) = 1.0 - zc_dsc(i,j) / zhsc(i,j)
        ELSE
          zeta_r_dsc(i,j) = 1.0 - zc_dsc(i,j) / dscdepth(i,j)
        END IF

        IF (db_dsct_cld(i,j)  >=  0.0) THEN
            !----------------------------
            ! i.e. no buoyancy reversal
            !----------------------------
          db_dsct_cld(i,j) = 0.0
          d_siems_dsc(i,j) = 0.0
          br_fback_dsc(i,j)= 0.0
        ELSE
            !----------------------------
            ! IF (DB_DSCT_CLD(I,j)  <   0.0)
            ! i.e. buoyancy reversal
            !----------------------------
          db_dsct_cld(i,j) = -db_dsct_cld(i,j) * cld_factor_dsc(i,j)
          d_siems_dsc(i,j) = MAX( 0.0, chi_s_dsct(i,j)                  &
                          * db_dsct_cld(i,j) / (db_dsct(i,j)+1.e-14) )
            ! Linear feedback dependence for D<0.1
          br_fback_dsc(i,j)= MIN( 1.0, 10.0*d_siems_dsc(i,j) )

          IF ( entr_enhance_by_cu == Buoyrev_feedback                   &
               .AND. cumulus(i,j)                                       &
               .AND. d_siems_dsc(i,j) < 0.1                             &
               .AND. d_siems_dsc(i,j) > 1.e-14 ) THEN
              ! Assume mixing from cumulus can enhance the
              ! buoyancy reversal feedback in regime 0<D<0.1.
              ! Make enhancement dependent on Cu depth:
              !       Cu_depth_fac->0 below 400m, 1 above 1000m
            cu_depth_fac = 0.5*( 1.0+                                   &
                      TANH( ((zhpar(i,j)-zh(i,j))-700.0)/100.) )
              ! BR_FBACK = unchanged for Cu<400m, ->1 for Cu>1000.
            br_fback_dsc(i,j) = cu_depth_fac +                          &
                               (1.0-cu_depth_fac)*br_fback_dsc(i,j)
          END IF

          zeta_r_dsc(i,j) = zeta_r_dsc(i,j) +                           &
                           (1.0-zeta_r_dsc(i,j))*br_fback_dsc(i,j)

        END IF
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 8. Calculate the radiative flux change across cloud top for mixed-
!    layer to ZH.  Restrict search for maximum divergence to below
!    NTML+2.  This may introduce errors if NTML changes a lot during
!    the radiative timestep but can't be helped.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_cloud_top(i,j) = 0
      df_top_over_cp(i,j) = 0.0
      df_inv_sml(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP SINGLE
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
          ! restrict search to `close' to ZH
        k_rad_lim = ntml(i,j)+1
          !-------------------------------------------------------------
          ! Find the layer below K_RAD_LIM with the greatest LW
          ! radiative flux jump in the upper half of the BL
          ! and assume that this is the top of the SML.
          !-------------------------------------------------------------
        IF (dflw_over_cp(i,j,k)  >   df_top_over_cp(i,j)                &
                    .AND. k  <=  k_rad_lim                              &
                    .AND. z_tq(i,j,k)  >   0.5*zh(i,j) ) THEN
          k_cloud_top(i,j) = k
          IF ( k >  1 ) THEN
              ! Set K_CLOUD_TOP to the level below if its DF is
              ! greater than half the maximum.  DF in level
              ! K_CLOUD_TOP+1 is then included as DF_INV_SML below.
            IF (dflw_over_cp(i,j,k-1)  >   0.5*dflw_over_cp(i,j,k))     &
              k_cloud_top(i,j) = k-1
          END IF
          df_top_over_cp(i,j) = dflw_over_cp(i,j,k)
        END IF

      END DO
    END DO
  END DO
!$OMP END SINGLE

      !-----------------------------------------------------------------
      !  Find bottom grid-level (K_LEVEL) for cloud-top radiative fux
      !  divergence: higher of base of LW radiatively cooled layer,
      !  0.5*ZH, since cooling must be in upper part of layer
      !-----------------------------------------------------------------

!$OMP SINGLE
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_level(i,j) = k_cloud_top(i,j)
      IF ( k_cloud_top(i,j)  >   1 ) THEN
        k_rad_lim = 1
        k=k_cloud_top(i,j)-1
        kl=MAX(1,k)  ! only to avoid out-of-bounds compiler warning
        DO WHILE ( k  >   k_rad_lim                                     &
                  .AND. dflw_over_cp(i,j,kl)  >   0.0                   &
                  .AND. z_tq(i,j,kl)  >   0.5*zh(i,j) )
          k_level(i,j) = k
          k = k-1
          kl=MAX(1,k)

        END DO
      END IF
    END DO
  END DO
!$OMP END SINGLE

      !-----------------------------------------------------------------
      ! Calculate LW and SW flux divergences and combine into
      ! cloud-top turbulence forcing.
      ! Need to account for radiative divergences in cloud in inversion
      ! grid-level, DF_INV_SML. Assume DF_OVER_CP(K_cloud_top+2) is
      ! representative of clear-air rad divergence and so subtract this
      ! `clear-air' part from the grid-level divergence.
      !-----------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF ( k_cloud_top(i,j)  >   0 ) THEN
        dflw_inv = 0.0
        dfsw_inv = 0.0
        IF ( k_cloud_top(i,j) <  bl_levels ) THEN
          k = k_cloud_top(i,j)+1
          IF ( k  <   bl_levels ) THEN
            dflw_inv = dflw_over_cp(i,j,k)                              &
                       - dflw_over_cp(i,j,k+1)                          &
                              * dzl(i,j,k)/dzl(i,j,k+1)
            dfsw_inv = dfsw_over_cp(i,j,k)                              &
                       - dfsw_over_cp(i,j,k+1)                          &
                              * dzl(i,j,k)/dzl(i,j,k+1)
          ELSE
            dflw_inv = dflw_over_cp(i,j,k)
            dfsw_inv = dfsw_over_cp(i,j,k)
          END IF
          dflw_inv = MAX( dflw_inv, 0.0 )
          dfsw_inv = MIN( dfsw_inv, 0.0 )
        END IF
        df_inv_sml(i,j) = dflw_inv + dfsw_inv

        df_top_over_cp(i,j) = frad_lw(i,j,k_cloud_top(i,j)+1)           &
                             - frad_lw(i,j,k_level(i,j))                &
                             + dflw_inv

        dfsw_top = frad_sw(i,j,k_cloud_top(i,j)+1)                      &
                 - frad_sw(i,j,k_level(i,j))                            &
                 + dfsw_inv

          !-----------------------------------------------------------
          ! Combine SW and LW cloud-top divergences into a net
          ! divergence by estimating SW flux divergence at a given
          ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
          ! Empirically (from LEM data) a reasonable fit is found
          ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = dfsw_frac
          !-----------------------------------------------------------
        df_top_over_cp(i,j) = MAX( 0.0,                                 &
                      df_top_over_cp(i,j) + dfsw_frac * dfsw_top )
      END IF
    END DO
  END DO
!$OMP END DO

! ------------------------------------------------------------------
! 9. Calculate the non-turbulent fluxes at the layer boundaries.
!  - the radiative flux at the inversion allows for an estimate
!    of the FA flux divergence within the inversion grid-level
!  - because the radiative time-step is usually longer the radiative
!    cloud-top grid-level (K_CLOUD_TOP) is allowed to differ from
!    the actual one (NTML)
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
!  - if no rad cooling was identified in layer, need to set
!    K_CLOUD_TOP and K_CLOUD_DSCT to top level in mixed layer
! ------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( k_cloud_top(i,j)  ==  0 ) k_cloud_top(i,j) = ntml(i,j)

      ft_nt_zh(i,j)   = frad(i,j,k_cloud_top(i,j)+1)                    &
                           + df_inv_sml(i,j)
      ft_nt_zh(i,j)   = ft_nt_zh(i,j)   + fmic(i,j,ntml(i,j)+2,1)       &
                                        + fsubs(i,j,ntml(i,j),1)
      fq_nt_zh(i,j)   = fmic(i,j,ntml(i,j)+2,2)                         &
                      + fsubs(i,j,ntml(i,j),2)

      ft_nt_zhsc(i,j) = 0.0
      ft_nt_zhsc(i,j) = 0.0
      fq_nt_zhsc(i,j) = 0.0
      IF ( dsc(i,j) ) THEN
        IF ( k_cloud_dsct(i,j)  ==  0 ) k_cloud_dsct(i,j) = ntdsc(i,j)
        ft_nt_zhsc(i,j) = frad(i,j,k_cloud_dsct(i,j)+1)                 &
                             + df_inv_dsc(i,j)
        ft_nt_zhsc(i,j) = ft_nt_zhsc(i,j) + fmic(i,j,ntdsc(i,j)+2,1)    &
                                          + fsubs(i,j,ntdsc(i,j),1)
        fq_nt_zhsc(i,j) = fmic(i,j,ntdsc(i,j)+2,2)                      &
                        + fsubs(i,j,ntdsc(i,j),2)
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 10.  Subroutine EXCF_NL.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      ntml_save(i,j) = ntml(i,j)  ! needed to identify changes
      dsc_save(i,j)  = dsc(i,j)   !      in excf_nl
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

! DEPENDS ON: excf_nl
  CALL excf_nl (                                                        &
! IN levels/switches
   bl_levels,nSCMDpkgs,L_SCMDiags,                                      &
! IN fields
   rdz,z_uv,z_tq,rho_uv,rho_tq,rhostar_gb,v_s,fb_surf,zhpar,            &
   btm,bqm,btm_cld,bqm_cld,cf_sml,cf_dsc,                               &
   bflux_surf,bflux_surf_sat,zeta_s,svl_diff_frac,                      &
   df_top_over_cp,zeta_r,bt_top,btt_top,btc_top,                        &
   db_top,db_top_cld,chi_s_top,br_fback,                                &
   df_dsct_over_cp,zeta_r_dsc,bt_dsct,btt_dsct,                         &
   db_dsct,db_dsct_cld,chi_s_dsct,br_fback_dsc,                         &
   db_ga_dry,db_noga_dry,db_ga_cld,db_noga_cld,                         &
   dsl_sml,dqw_sml,dsl_dsc,dqw_dsc,ft_nt,fq_nt,                         &
! INOUT fields
   dsc,cumulus,coupled,ntml,zh,zhsc,dscdepth,ntdsc,zc,zc_dsc,           &
   ft_nt_zh, ft_nt_zhsc, fq_nt_zh, fq_nt_zhsc,                          &
! OUT fields
   rhokm, rhokh, rhokm_top, rhokh_top,                                  &
   rhokh_top_ent, rhokh_dsct_ent, rhokh_surf_ent,                       &
   rhof2,rhofsc,f_ngstress,zdsc_base,nbdsc                              &
  )

!$OMP  PARALLEL DEFAULT(SHARED)                                         &
!$OMP& PRIVATE (i, j, k, wstar3, w_m, w_h, ml_tend, fa_tend, inv_tend,  &
!$OMP& tothf_efl, moisten, totqf_efl, w_s_ent, w_s_cubed)

!-----------------------------------------------------------------------
!-adjust SML/DSC properties depending on diagnoses in EXCF_NL
! Note that the non-turbulent fluxes at inversions will have been
! swapped in EXCF_NL (ie. FT/Q_NT_ZH/ZHSC)
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( dsc(i,j) .AND. .NOT.dsc_save(i,j) ) THEN
!..decoupling diagnosed in EXCF_NL - change parameters around
        dsl_dsc(i,j) = dsl_sml(i,j)
        dqw_dsc(i,j) = dqw_sml(i,j)
        db_dsct(i,j) = db_top(i,j)  ! copy diagnostics across
        zc_dsc(i,j)  = zc(i,j)      !
        dsc_disc_inv(i,j) = sml_disc_inv(i,j)
        sml_disc_inv(i,j) = 0
        dsl_sml(i,j) = 0.0
        dqw_sml(i,j) = 0.0
        df_inv_dsc(i,j) = df_inv_sml(i,j)
        df_inv_sml(i,j) = 0.0
      END IF
      IF ( .NOT.dsc(i,j) .AND. dsc_save(i,j) ) THEN
!..decoupled layer removed in EXCF_NL; either...
        IF ( ntml_save(i,j)  ==  ntml(i,j) ) THEN
          !...had no turbulence forcing
          dsl_dsc(i,j) = 0.0
          dqw_dsc(i,j) = 0.0
          dsc_disc_inv(i,j) = 0
          df_inv_dsc(i,j) = 0.0
        ELSE
          !...recoupled with surface layer
          dsl_sml(i,j) = dsl_dsc(i,j)
          dqw_sml(i,j) = dqw_dsc(i,j)
          dsl_dsc(i,j) = 0.0
          dqw_dsc(i,j) = 0.0
          sml_disc_inv(i,j) = dsc_disc_inv(i,j)
          dsc_disc_inv(i,j) = 0
          df_inv_sml(i,j) = df_inv_dsc(i,j)
          df_inv_dsc(i,j) = 0.0
        END IF
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 11.  Calculate "explicit" entrainment fluxes of SL and QW.
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the DSC base
! ------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      ft_nt_dscb(i,j) = ft_nt(i,j,1)
      IF ( nbdsc(i,j)  >   1 ) THEN
        k = nbdsc(i,j)  ! NBDSC marks the lowest flux-level
                           !    within the DSC layer
                           ! Interpolate non-turb flux to base
                           !    of DSC layer:
        ft_nt_dscb(i,j) = ft_nt(i,j,k-1) +                              &
                  (ft_nt(i,j,k)-ft_nt(i,j,k-1))                         &
                 *(zdsc_base(i,j)-z_uv(i,j,k-1))/dzl(i,j,k-1)
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!..Specify entrainment fluxes at NTML+1 and NTDSC+1 directly through FTL
!..and FQW (and set the entrainment RHOKH to zero).
!..The turbulent flux at the subgrid ZH is given by the entrainment
!..parametrization ( = -w_e*'jump across inversion').  Together with
!..the non-turbulent flux profile (rad+microphys+subs), this gives the
!..total flux at the subgrid ZH.  The linear total flux profile is then
!..interpolated onto the half-level below ZH (the entrainment flux
!..grid-level).  The total flux divergence across the inversion grid
!..level is then checked for consistency with the entrainment/subsidence
!..balance, eg. a falling inversion should warm the inversion grid-level
!..Finally, the non-turbulent flux is sutracted off to give the required
!..turbulent component at the entrainment flux grid-level.
!..For momentum, given the horizontal interpolation required, together
!..with the lack of accuracy in assuming a discontinuous inversion,
!..entrainment continues to be specified using the specified RHOKM.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        ftl(i,j,k) = 0.0
        fqw(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!..First the surface-based mixed layer (if entraining and a
!..discontinuous inversion structure was diagnosed -
!..ie. the inversion is well-defined)
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      zh_np1(i,j)   = 0.0
      t_frac(i,j)   = 0.0
      zrzi(i,j)     = 0.0
      we_rho(i,j)   = 0.0
      tothf_zh(i,j) = 0.0
      totqf_zh(i,j) = 0.0
      k=ntml(i,j)+1

        ! Only RHOKH_ENT is passed out of EXCFNL so recalculate WE:
      we_parm(i,j) = rdz(i,j,k)*                                        &
                         ( rhokh_top_ent(i,j)+rhokh_surf_ent(i,j) )     &
                                                    / rho_uv(i,j,k)

      IF ( sml_disc_inv(i,j)  ==  1 .AND. .NOT.coupled(i,j) .AND.       &
           (rhokh_top_ent(i,j)+rhokh_surf_ent(i,j))  >   0.0 ) THEN

!-----------------------------------------------------------------------
!..Calculate ZH at end of timestep, ZH_NP1
!-----------------------------------------------------------------------
!..linearly interpolate vertical velocity to ZH
        IF ( zh(i,j)  >=  z_tq(i,j,k) ) THEN
          w_ls(i,j) = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) )              &
                      * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
        ELSE
          w_ls(i,j) = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) )              &
                      * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
        END IF
        w_ls(i,j) = MIN ( w_ls(i,j), 0.0 )
          ! only interested in subsidence

        zh_np1(i,j) = zh(i,j) +                                         &
                        timestep * ( we_parm(i,j) + w_ls(i,j) )
        zh_np1(i,j) = MAX( zh_np1(i,j), z_uv(i,j,k-1) )
        IF ( zh_np1(i,j)  >   z_top(i,j,k+1) ) THEN
            ! limit ZH and W_e (and therefore the entraiment fluxes)
            ! because the inversion cannot rise more than one level
            ! in a timestep.
          zh_np1(i,j) = z_top(i,j,k+1)
          we_parm(i,j) =                                                &
                  (z_top(i,j,k+1) - zh(i,j))/timestep - w_ls(i,j)
        END IF
!-----------------------------------------------------------------------
!..Decide on which grid-level to apply entrainment flux
!-----------------------------------------------------------------------
        IF ( zh_np1(i,j)  >   z_uv(i,j,ntml(i,j)+2) ) THEN
            ! ZH risen above level K+1 so specify appropriate flux
            ! at this level and raise NTML by one (this means
            ! gradient-adjustment is also applied at half-level
            ! old_NTML+1).  Note KH profiles should already be
            ! calculated at level NTML+1 because ZH is above this level.
          ntml(i,j) = ntml(i,j) + 1
          k=ntml(i,j)+1

            ! T_FRAC is fraction of timestep inversion is above
            ! the entrainment flux grid-level (at Z_UV(K))
          t_frac(i,j) = (zh_np1(i,j)-z_uv(i,j,k)) /                     &
                        (zh_np1(i,j)-zh(i,j))
            ! ZH_FRAC is the timestep-average fraction of mixed layer
            ! air in the inversion grid-level, level NTML+1
          zh_frac(i,j) = 0.5*t_frac(i,j)*(zh_np1(i,j)-z_uv(i,j,k) )     &
                         / dzl(i,j,k)

        ELSE IF ( zh_np1(i,j)  >=  z_uv(i,j,ntml(i,j)+1) ) THEN
            ! ZH always between half-levels NTML+1 and NTML+2

          t_frac(i,j) = 1.0
          zh_frac(i,j) = ( 0.5*(zh(i,j)+zh_np1(i,j)) - z_uv(i,j,k) )    &
                         / dzl(i,j,k)

        ELSE
            ! ZH falls below half-level NTML+1
            ! Keep implicit (diffusive) entrainment but apply
            ! at the level below
          IF (ntml(i,j)  >=  2) THEN     ! ftl(k=1) is surface flux!
            ntml(i,j) = ntml(i,j) - 1
            k=ntml(i,j)+1
            rhokh_top(i,j,k+1) = 0.0   ! also need to remove diffusion
            rhokh(i,j,k+1)     = 0.0   ! at old entrainment grid-level
          END IF
          t_frac(i,j) = 0.0
          zh_frac(i,j) = 0.0
          sml_disc_inv(i,j) = 0

        END IF  ! test on where to apply entrainment flux

        we_rho(i,j) = rho_uv(i,j,k) * we_parm(i,j)
        zrzi(i,j)   = z_uv(i,j,k)*2.0/(zh(i,j)+zh_np1(i,j))

      END IF   ! test on SML_DISC_INV, etc
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the surface and the parametrized flux at the inversion
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

        ! Entrainment flux applied to level NTML+1 which is the
        ! flux-level above the top of the SML
      k=ntml(i,j)+1

      IF ( t_frac(i,j)  >   0.0 ) THEN

        rhokh_top(i,j,k) = 0.0   ! apply entrainment explicitly
        rhokh(i,j,k) = 0.0       !      "

        tothf_zh(i,j) = - we_rho(i,j)*dsl_sml(i,j) + ft_nt_zh(i,j)
          ! Linearly interpolate to entrainment flux grid-level
        tothf_efl = ft_nt(i,j,1) + ftl(i,j,1) +                         &
                   ( tothf_zh(i,j)-ft_nt(i,j,1)-ftl(i,j,1) )*zrzi(i,j)
          ! Ensure total heat flux gradient in inversion grid-level is
          ! consistent with inversion rising (ie. implies cooling in
          ! level K relative to the mixed layer) or falling
          ! (implies warming)

        ml_tend = -( tothf_zh(i,j)-ft_nt(i,j,1)-ftl(i,j,1) ) / zh(i,j)
        fa_tend = 0.0
        IF ( k+1  <=  bl_levels )                                       &
            fa_tend = - ( ft_nt(i,j,k+2) - ft_nt(i,j,k+1) )             &
                        / dzl(i,j,k+1)
        inv_tend =       zh_frac(i,j) * ml_tend                         &
                 + (1.0-zh_frac(i,j)) * fa_tend
        IF (we_parm(i,j)+w_ls(i,j)  >=  0.0) THEN
            ! Inversion moving up so inversion level should cool
            ! Ensure it does cool relative to ML
          tothf_efl = MIN( tothf_efl,                                   &
                           ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
            ! Ensure inversion level won't end up colder than
            ! NTML by end of timestep.
            ! Set INV_TEND to max allowable cooling rate, also
            ! allowing for change in ML_TEND arising from this change
            ! to TOTHF_EFL:
          inv_tend = (sl(i,j,k-1)-sl(i,j,k))/timestep                   &
                      + (ft_nt(i,j,1)+ftl(i,j,1))/z_uv(i,j,k)
          tothf_efl = MAX( tothf_efl,                                   &
                           (ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k))         &
                                /(1.+ dzl(i,j,k)/z_uv(i,j,k)) )
        ELSE  ! WE_PARM+W_LS < 0
            ! Ensure inversion level does warm relative to ML
          tothf_efl = MAX( tothf_efl,                                   &
                           ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        END IF
          ! Turbulent entrainment flux is then the residual of the total
          ! flux and the net flux from other processes
        ftl(i,j,k) =  t_frac(i,j) * ( tothf_efl - ft_nt(i,j,k) )

      ELSE   ! NOT specifying entrainment flux but KH
          ! Include entrainment KH in K-profiles, if greater
          ! (for COUPLED layers these will be zero)
        rhokh_top(i,j,k) = MAX( rhokh_top(i,j,k), rhokh_top_ent(i,j) )
        rhokh(i,j,k)     = MAX( rhokh(i,j,k), rhokh_surf_ent(i,j) )

      END IF  ! test on T_FRAC gt 0

    END DO
  END DO
!$OMP END DO

!-------------------------------------------------
!..Second the decoupled mixed layer, if entraining
!-------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      zhsc_np1(i,j)   = 0.0
      t_frac_dsc(i,j) = 0.0
      zrzi_dsc(i,j)   = 0.0
      we_rho_dsc(i,j) = 0.0
      tothf_zhsc(i,j) = 0.0
      totqf_zhsc(i,j) = 0.0

      k=ntdsc(i,j)+1
      we_dsc_parm(i,j) = rdz(i,j,k)*rhokh_dsct_ent(i,j)                 &
                                                 / rho_uv(i,j,k)

      IF ( dsc_disc_inv(i,j)  ==  1                                     &
                    .AND. rhokh_dsct_ent(i,j)  >   0.0 ) THEN

!-----------------------------------------------------------------------
!..Calculate ZHSC at end of timestep, ZHSC_NP1
!-----------------------------------------------------------------------
!..interpolate vertical velocity to ZH
        IF ( zhsc(i,j)  >=  z_tq(i,j,k) ) THEN
          w_ls_dsc(i,j) = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) ) *        &
                           (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
        ELSE
          w_ls_dsc(i,j) = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) ) *        &
                           (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
        END IF
        w_ls_dsc(i,j) = MIN ( w_ls_dsc(i,j), 0.0 )
          ! only interested in subsidence

        zhsc_np1(i,j) = zhsc(i,j) +                                     &
              timestep * ( we_dsc_parm(i,j) + w_ls_dsc(i,j) )
        zhsc_np1(i,j) = MAX( zhsc_np1(i,j), z_uv(i,j,k-1) )
        IF ( zhsc_np1(i,j)  >   z_top(i,j,k+1) ) THEN
            ! limit ZHSC and W_e (and therefore the entrainment fluxes)
            ! because the inversion cannot rise more than one level
            ! in a timestep.
          zhsc_np1(i,j) = z_top(i,j,k+1)
          we_dsc_parm(i,j) =                                            &
             (z_top(i,j,k+1) - zhsc(i,j))/timestep - w_ls_dsc(i,j)
        END IF
!-----------------------------------------------------------------------
!..Decide on which grid-level to apply entrainment flux
!-----------------------------------------------------------------------
        IF ( zhsc_np1(i,j)  >   z_uv(i,j,ntdsc(i,j)+2) ) THEN
            ! ZHSC risen above level K+1 so specify appropriate
            ! flux at this level and raise NTDSC by one

          ntdsc(i,j) = ntdsc(i,j) + 1
          k = ntdsc(i,j)+1
          t_frac_dsc(i,j) = (zhsc_np1(i,j)-z_uv(i,j,k)) /               &
                            (zhsc_np1(i,j)-zhsc(i,j))

          zhsc_frac(i,j) = 0.5*t_frac_dsc(i,j)*                         &
                           ( zhsc_np1(i,j)-z_uv(i,j,k) )/ dzl(i,j,k)

        ELSE IF ( zhsc_np1(i,j)  >   z_uv(i,j,ntdsc(i,j)+1) ) THEN
            ! ZHSC always between half-levels NTDSC+1 and NTDSC+2

          t_frac_dsc(i,j) = 1.0
          zhsc_frac(i,j) = ( 0.5*(zhsc(i,j)+zhsc_np1(i,j))              &
                                         - z_uv(i,j,k) )/ dzl(i,j,k)

        ELSE
            ! ZHSC falls below half-level NTDSC+1
            ! Keep implicit (diffusive) entrainment but apply
            ! at the level below
          ntdsc(i,j) = ntdsc(i,j) - 1  ! could reduce NTDSC to 1
          k = ntdsc(i,j)+1
          rhokh_top(i,j,k+1) = 0.0
          rhokh(i,j,k+1)     = 0.0

          t_frac_dsc(i,j)   = 0.0
          zhsc_frac(i,j)    = 0.0
          dsc_disc_inv(i,j) = 0

        END IF  ! test on where to apply entrainment flux

        we_rho_dsc(i,j) = rho_uv(i,j,k) * we_dsc_parm(i,j)
          ! for z'/z_i' assume height of DSC base is fixed in time
        zrzi_dsc(i,j) =( z_uv(i,j,k)-(zhsc(i,j)-dscdepth(i,j)) )        &
                      /( dscdepth(i,j)+0.5*(zhsc_np1(i,j)-zhsc(i,j)) )

      END IF   ! test on DSC_DISC_INV, etc
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!..Linearly interpolate between the known total (turb+rad+subs+micro)
!..flux at the DSC base and the parametrized flux at the inversion
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

        ! Entrainment flux applied to level NTDSC+1 which is the
        ! flux-level above the top of the DSC layer
      k=ntdsc(i,j)+1

      IF ( t_frac_dsc(i,j)  >   0.0 ) THEN


        rhokh_top(i,j,k) = 0.0   ! apply entrainment explicitly
        rhokh(i,j,k)     = 0.0   !      "

        tothf_zhsc(i,j) = - we_rho_dsc(i,j)*dsl_dsc(i,j)                &
                                + ft_nt_zhsc(i,j)
        tothf_efl = ft_nt_dscb(i,j) +                                   &
                    ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )*zrzi_dsc(i,j)
          ! Ensure total heat flux gradient in inversion grid-level is
          ! consistent with inversion rising (implies cooling in
          ! level K, relative to mixed layer) or falling
          ! (implies warming)
        ml_tend = - ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )/ dscdepth(i,j)
        fa_tend = 0.0
        IF ( k+1  <=  bl_levels )                                       &
            fa_tend = - ( ft_nt(i,j,k+2) - ft_nt(i,j,k+1) )             &
                        / dzl(i,j,k+1)
        inv_tend =       zhsc_frac(i,j) * ml_tend                       &
                 + (1.0-zhsc_frac(i,j)) * fa_tend

        IF (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  0.0) THEN
            ! Inversion moving up so inversion level should cool
            ! Ensure it does cool relative to ML
          tothf_efl = MIN( tothf_efl,                                   &
                           ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
            ! Ensure inversion level won't end up colder than
            ! NTDSC by end of timestep.
          inv_tend = (sl(i,j,k-1)-sl(i,j,k))/timestep                   &
                     + ft_nt_dscb(i,j)/dscdepth(i,j)
          tothf_efl = MAX( tothf_efl,                                   &
                        (ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k))            &
                         /(1.+ dzl(i,j,k)/dscdepth(i,j))   )
        ELSE   ! WE_DSC_PARM+W_LS_DSC < 0
          tothf_efl = MAX( tothf_efl,                                   &
                           ft_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
        END IF
          ! Turbulent entrainment flux is then the residual of the total
          ! flux and the net flux from other processes
        ftl(i,j,k) = t_frac_dsc(i,j) * ( tothf_efl - ft_nt(i,j,k) )

      ELSE IF ( dsc(i,j) ) THEN

          ! Not specifying entrainment flux but KH
          ! Include entrainment KH in K-profile, if greater
        rhokh_top(i,j,k) = MAX( rhokh_top(i,j,k),rhokh_dsct_ent(i,j) )

      END IF  ! IF NOT DSC

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Specify QW entrainment fluxes
!-----------------------------------------------------------------------
! Calculate the non-turbulent fluxes at the layer boundaries
!  - the subsidence flux at the inversion is taken from the
!    flux grid-level below it (assumes the divergence across
!    the inversion is physically above the BL)
!  - the microphysical flux at the inversion is taken from the
!    flux grid-level just above it (assumes the divergence across
!    the inversion grid-level is physically within the BL)
! ------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      fq_nt_dscb(i,j) = fq_nt(i,j,1)
      IF ( nbdsc(i,j)  >   1 ) THEN
        k = nbdsc(i,j)  ! NBDSC marks the lowest flux-level
                           !    within the DSC layer
                           ! Interpolate non-turb flux to base
                           !    of DSC layer:
        fq_nt_dscb(i,j) = fq_nt(i,j,k-1) +                              &
                  (fq_nt(i,j,k)-fq_nt(i,j,k-1))                         &
                 *(zdsc_base(i,j)-z_uv(i,j,k-1))/dzl(i,j,k-1)
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Calculate grid-level QW fluxes at inversion, ensuring the turbulent,
! microphysical and subsidence fluxes are correctly coupled.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k = ntml(i,j)+1
      IF ( t_frac(i,j)  >   0.0 ) THEN
          ! Calculate total (turb+micro+subs) QW flux at subgrid
          ! inversion height
        totqf_zh(i,j) = - we_rho(i,j)*dqw_sml(i,j) + fq_nt_zh(i,j)
          ! Interpolate to entrainment flux-level below
        totqf_efl = fq_nt(i,j,1) + fqw(i,j,1) +  zrzi(i,j) *            &
                         ( totqf_zh(i,j) - fq_nt(i,j,1) - fqw(i,j,1) )
          ! Need to ensure the total QW flux gradient in inversion
          ! grid-level is consistent with inversion rising or falling.
          ! If QW(K) is drier than mixed layer then inversion rising
          ! implies moistening in level K relative to mixed layer
          ! while falling would imply relative drying of level K.
          ! If QW(K) is moister than ML then want opposite tendencies.
        ml_tend = - ( totqf_zh(i,j)-fq_nt(i,j,1)-fqw(i,j,1) ) /zh(i,j)
        fa_tend = 0.0
        IF ( k+1  <=  bl_levels )                                       &
          fa_tend = - ( fq_nt(i,j,k+2)-fq_nt(i,j,k+1) )                 &
                      / dzl(i,j,k+1)
        inv_tend =       zh_frac(i,j) * ml_tend                         &
                 + (1.0-zh_frac(i,j)) * fa_tend

        IF (we_parm(i,j)+w_ls(i,j) >=  0.0) THEN
            ! inversion moving up so inversion will moisten/dry
            ! depending on relative QW in level below
          moisten = ( qw(i,j,k) <= qw(i,j,k-1) )
        ELSE
            ! inversion moving down so inversion will moisten/dry
            ! depending on relative QW in level above
          moisten = ( qw(i,j,k) <= qw(i,j,k+1) )
        END IF

        IF ( moisten ) THEN
            ! Ensure inversion level does moisten relative to ML
          totqf_efl = MAX( totqf_efl,                                   &
                          fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          IF (we_parm(i,j)+w_ls(i,j)  >=  0.0) THEN
              ! Ensure inversion level won't end up more moist than
              ! NTML by end of timestep.
              ! Set INV_TEND to max allowable moistening rate, also
              ! allowing for change in ML_TEND arising from this change
              ! to TOTQF_EFL:
            inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                 &
                        + (fq_nt(i,j,1)+fqw(i,j,1))/z_uv(i,j,k)
            totqf_efl = MIN( totqf_efl,                                 &
                          (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))          &
                           /(1.+ dzl(i,j,k)/z_uv(i,j,k))   )
          END IF
        ELSE
          totqf_efl = MIN( totqf_efl,                                   &
                          fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          IF (we_parm(i,j)+w_ls(i,j)  >=  0.0) THEN
              ! Ensure inversion level won't end up drier than
              ! NTML by end of timestep.
              ! Set INV_TEND to max allowable drying rate:
            inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                 &
                        + (fq_nt(i,j,1)+fqw(i,j,1))/z_uv(i,j,k)
            totqf_efl = MAX( totqf_efl,                                 &
                          (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))          &
                           /(1.+ dzl(i,j,k)/z_uv(i,j,k))   )
          END IF
        END IF
        fqw(i,j,k) = t_frac(i,j) *                                      &
                         ( totqf_efl - fq_nt(i,j,k) )
      END IF

    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Now decoupled layer
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF ( t_frac_dsc(i,j)  >   0.0 ) THEN

        k = ntdsc(i,j)+1

          ! Calculate total (turb+micro) QW flux at subgrid inversion
        totqf_zhsc(i,j) = - we_rho_dsc(i,j)*dqw_dsc(i,j)                &
                            + fq_nt_zhsc(i,j)
          ! Interpolate to entrainment flux-level
        totqf_efl = fq_nt_dscb(i,j) +                                   &
                  ( totqf_zhsc(i,j) - fq_nt_dscb(i,j) )*zrzi_dsc(i,j)

        ml_tend = - ( totqf_zhsc(i,j)-fq_nt_dscb(i,j) )/dscdepth(i,j)
        fa_tend = 0.0
        IF ( k+1  <=  bl_levels )                                       &
           fa_tend = - ( fq_nt(i,j,k+2)-fq_nt(i,j,k+1) )                &
                       / dzl(i,j,k+1)
        inv_tend =       zhsc_frac(i,j) * ml_tend                       &
                 + (1.0-zhsc_frac(i,j)) * fa_tend

        IF (we_dsc_parm(i,j)+w_ls_dsc(i,j) >=  0.0) THEN
            ! inversion moving up so inversion will moisten/dry
            ! depending on relative QW in level below
          moisten = ( qw(i,j,k) <= qw(i,j,k-1) )
        ELSE
            ! inversion moving down so inversion will moisten/dry
            ! depending on relative QW in level above
          moisten = ( qw(i,j,k) <= qw(i,j,k+1) )
        END IF

        IF ( moisten ) THEN
          totqf_efl = MAX( totqf_efl,                                   &
                            fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          IF (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  0.0) THEN
              ! Ensure inversion level won't end up more moist than
              ! NTDSC by end of timestep.
            inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                 &
                       + fq_nt_dscb(i,j)/dscdepth(i,j)
            totqf_efl = MIN( totqf_efl,                                 &
                          (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))          &
                           /(1.+ dzl(i,j,k)/dscdepth(i,j))   )
          END IF
        ELSE
          totqf_efl = MIN( totqf_efl,                                   &
                            fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k) )
          IF (we_dsc_parm(i,j)+w_ls_dsc(i,j)  >=  0.0) THEN
              ! Ensure inversion level won't end up drier than
              ! NTDSC by end of timestep.
              ! Set INV_TEND to max allowable drying rate:
            inv_tend = (qw(i,j,k-1)-qw(i,j,k))/timestep                 &
                       + fq_nt_dscb(i,j)/dscdepth(i,j)
            totqf_efl = MAX( totqf_efl,                                 &
                          (fq_nt(i,j,k+1)+inv_tend*dzl(i,j,k))          &
                           /(1.+ dzl(i,j,k)/dscdepth(i,j))   )
          END IF
        END IF
        fqw(i,j,k) = t_frac_dsc(i,j) * ( totqf_efl - fq_nt(i,j,k) )
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Calculate effective entrainment (ie. reduced to allow for subsidence
! increments in the ML) for use in tracer mixing.  Take theta_l as a
! representative scalar field since jump should always be the same sign
! and code therefore simpler.
! In this version the inversion fluxes are only implemented at one
! grid-level so only one element of these 3D arrays is used.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, 3
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        we_lim(i,j,k)    = 0.0
        t_frac_tr(i,j,k) = 0.0
        zrzi_tr(i,j,k)   = 0.0
        we_lim_dsc(i,j,k)    = 0.0
        t_frac_dsc_tr(i,j,k) = 0.0
        zrzi_dsc_tr(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      kent(i,j)   = ntml(i,j)+1
      t_frac_tr(i,j,2) = t_frac(i,j)
      zrzi_tr(i,j,2) = zrzi(i,j)
      IF ( t_frac(i,j)  >   0.0 ) THEN
        w_s_ent = 0.0
        k = ntml(i,j)
        IF ( dsl_sml(i,j)  /=  0.0 ) w_s_ent =                          &
            MIN( 0.0, -sls_inc(i,j,k) * dzl(i,j,k) /dsl_sml(i,j) )
          ! Only allow w_e to be reduced to zero!
        we_lim(i,j,2) = rho_uv(i,j,k+1) *                               &
                          MAX( 0.0, we_parm(i,j) + w_s_ent )
      END IF
      kent_dsc(i,j)   = ntdsc(i,j)+1
      t_frac_dsc_tr(i,j,2) = t_frac_dsc(i,j)
      zrzi_dsc_tr(i,j,2) = zrzi_dsc(i,j)
      IF ( t_frac_dsc(i,j)  >   0.0 ) THEN
        w_s_ent = 0.0
        k = ntdsc(i,j)
        IF ( dsl_dsc(i,j)  /=  0.0 ) w_s_ent =                          &
            MIN( 0.0, -sls_inc(i,j,k) * dzl(i,j,k) /dsl_dsc(i,j) )
          ! Only allow w_e to be reduced to zero!
        we_lim_dsc(i,j,2) = rho_uv(i,j,k) *                             &
                          MAX( 0.0, we_dsc_parm(i,j) + w_s_ent )
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 12. Update standard deviations and gradient adjustment to use this
!     timestep's ZH (code from SF_EXCH)
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( unstable(i,j) ) THEN
        IF (flux_grad  ==  Locketal2000) THEN
          w_s_cubed = 0.25 * zh(i,j) * fb_surf(i,j)
          IF (w_s_cubed  >   0.0) THEN
            w_m  =                                                      &
           ( w_s_cubed + v_s(i,j) * v_s(i,j) * v_s(i,j) ) ** (1.0/3.0)
            t1_sd(i,j) = 1.93 * ftl(i,j,1) / (rhostar_gb(i,j) * w_m)
            q1_sd(i,j) = 1.93 * fqw(i,j,1) / (rhostar_gb(i,j) * w_m)
            tv1_sd(i,j) = t(i,j,1) *                                    &
              ( 1.0 + c_virtual*q(i,j,1) - qcl(i,j,1) - qcf(i,j,1) ) *  &
              ( bt(i,j,1)*t1_sd(i,j) + bq(i,j,1)*q1_sd(i,j) )
            t1_sd(i,j) = MAX ( 0.0 , t1_sd(i,j) )
            q1_sd(i,j) = MAX ( 0.0 , q1_sd(i,j) )
            IF (tv1_sd(i,j)  <=  0.0) THEN
              tv1_sd(i,j) = 0.0
              t1_sd(i,j) = 0.0
              q1_sd(i,j) = 0.0
            END IF
          END IF
          grad_t_adj(i,j) = MIN( max_t_grad ,                           &
                           a_grad_adj * t1_sd(i,j) / zh(i,j) )
          grad_q_adj(i,j) = 0.0
        ELSE IF (flux_grad  ==  HoltBov1993) THEN
            ! Use constants from Holtslag and Boville (1993)
            ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
            ! Neut limit GAMMA_TH = 7.2*wstar*FTL1/(ustar^2*zh)
          wstar3 = fb_surf(i,j) * zh(i,j)
          w_m =( v_s(i,j)**3 + 0.6*wstar3 )**(1.0/3.0)

          grad_t_adj(i,j) = a_ga_hb93*(wstar3**(1.0/3.0))*ftl(i,j,1)    &
                            / ( rhostar_gb(i,j)*w_m*w_m*zh(i,j) )
!            GRAD_Q_ADJ(I,j) = A_GA_HB93*(WSTAR3**(1.0/3.0))*FQW(I,j,1)
!     &                       / ( RHOSTAR_GB(I,j)*W_M*W_M*ZH(I,j) )
           ! Set q term to zero for same empirical reasons as Lock et al
          grad_q_adj(i,j) = 0.0
        ELSE IF (flux_grad  ==  LockWhelan2006) THEN
            ! Use constants LockWhelan2006
            ! Conv limit GAMMA_TH = 10 *FTL1/(wstar*zh)
            ! Neut limit GAMMA_TH = 7.5*FTL1/(ustar*zh)
          wstar3 = fb_surf(i,j) * zh(i,j)
          w_h =( ((4.0/3.0)*v_s(i,j))**3 + wstar3 )**(1.0/3.0)

          grad_t_adj(i,j) = a_ga_lw06 * ftl(i,j,1)                      &
                             / ( rhostar_gb(i,j)*w_h*zh(i,j) )
          grad_q_adj(i,j) = a_ga_lw06 * fqw(i,j,1)                      &
                             / ( rhostar_gb(i,j)*w_h*zh(i,j) )
        END IF
      END IF  ! test on UNSTABLE
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!- Save diagnostics
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( dsc(i,j) ) THEN
        IF (BL_diag%l_dscbase) THEN
          BL_diag%dscbase(i,j)= zhsc(i,j)-dscdepth(i,j)
        END IF
        IF (BL_diag%l_cldbase) THEN
          BL_diag%cldbase(i,j)= zhsc(i,j)-zc_dsc(i,j)
        END IF
        IF (BL_diag%l_weparm_dsc) THEN
          BL_diag%weparm_dsc(i,j)= we_dsc_parm(i,j)
        END IF
      ELSE
        IF (BL_diag%l_dscbase) THEN
          BL_diag%dscbase(i,j)= rmdi
        END IF
        IF (BL_diag%l_cldbase) THEN
          BL_diag%cldbase(i,j)= zh(i,j)-zc(i,j)
        END IF
        IF (BL_diag%l_weparm_dsc) THEN
          BL_diag%weparm_dsc(i,j)= we_parm(i,j)
        END IF
      END IF
      IF (BL_diag%l_weparm) THEN
        BL_diag%weparm(i,j)= we_parm(i,j)
      END IF
    END DO
  END DO
!$OMP END DO


!End of OpenMP parallel region
!$OMP END PARALLEL

  IF (lhook) CALL dr_hook('KMKHZ',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE kmkhz
