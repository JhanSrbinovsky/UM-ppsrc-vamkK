
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BDY_EXPL2----------------------------------------------
!
!  Purpose: Calculate the explicit turbulent fluxes of heat, moisture
!           and momentum between atmospheric levels
!           within the boundary layer, and/or the effects of these
!           fluxes on the primary model variables.
!
! Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 24.
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE bdy_expl2 (                                                  &
! IN values defining vertical grid of model atmosphere :
 bl_levels,p_theta_levels,land_pts,land_index,                          &
! IN U, V and W momentum fields.
 u_p,v_p,u_0_p,v_0_p,                                                   &
! IN variables for TKE scheme
 pstar,p_half,                                                          &
! IN from other part of explicit boundary layer code
 rho_uv,rho_tq,rho_dry_tq,dzl_charney,rdz,rdz_charney_grid,             &
 z_tq,z_uv,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt,  &
 recip_l_mo_sea,flandg,rib_gb, sil_orog_land, z0m_eff_gb,               &
! IN cloud/moisture data :
 cf,q,qcf,qcl,t,qw,tl,                                                  &
! IN everything not covered so far :
 rad_hr,micro_tends,fb_surf,u_s,h_blend_orog,                           &
 lq_mix_bl,zh_prev,nlcl,zhpar,z_lcl,ho2r2_orog,sd_orog,                 &
! SCM Diagnostics (dummy values in full UM) & stash diagnostics
 nSCMDpkgs,L_SCMDiags,BL_diag,                                          &
! INOUT variables
 zh,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot,t1_sd,q1_sd,   &
! INOUT variables on TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! OUT new variables for message passing
 tau_fd_x, tau_fd_y, f_ngstress, rhogamu, rhogamv,                      &
! OUT Diagnostic not requiring STASH flags :
 zht,shallowc,cu_over_orog,                                             &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,           &
 bl_type_7,                                                             &
! OUT data required for tracer mixing :
 ntml, kent, we_lim, t_frac, zrzi,                                      &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                            &
! OUT data required elsewhere in UM system :
 zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0                                   &
 )

  USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
  USE atm_fields_bounds_mod, ONLY: pdims, tdims, wdims, qdims, tdims_l, &
      pdims_s
  USE stochastic_physics_run_mod, ONLY: l_rp2, par_mezcla
  USE bl_diags_mod, ONLY: strnewbldiag
  USE bl_option_mod, ONLY: off, max_t_grad, a_grad_adj, sg_orog_mixing, &
      h_scale, t_drain, formdrag, explicit_stress, iDynDiag, DynDiag_ZL,&
      DynDiag_ZL_corrn, DynDiag_ZL_CuOnly,  DynDiag_Ribased,            &
      RiCrit_sharp, zhloc_depth_fac, non_local_bl, on, l_full_lambdas,  &
      nl_bl_levels, local_fa, free_trop_layers, to_sharp_across_1km,    &
      sbl_op, Equilibrium_SBL
  USE cv_run_mod, ONLY:                                                 &
      l_param_conv
  USE turb_diff_ctl_mod, ONLY:                                          &
      visc_m, visc_h, rneutml_sq, max_diff, delta_smag
  USE turb_diff_mod, ONLY:                                              &
      l_subfilter_vert, l_subfilter_horiz, l_subfilter_blend,           &
      turb_startlev_vert, turb_endlev_vert
  USE earth_constants_mod, ONLY: g    
  USE atmos_constants_mod, ONLY: vkman, cp    
  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
!$ USE omp_lib

  IMPLICIT NONE

!  Inputs :-
  INTEGER, INTENT(IN) ::                                                &
   land_pts,                                                            &
                               ! No.of land points in whole grid.
   bl_levels,                                                           &
                               ! IN Max. no. of "boundary" levels
   nlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! IN No of levels to LCL

!     Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

  REAL, INTENT(IN) ::                                                   &
    p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   bl_levels+1),                                        &
   rho_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1),                                                 &
                                   ! IN density on UV (ie. rho) levels;
                                   !    used in RHOKH so dry density if
                                   !    Lq_mix_bl is true
   rho_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                   ! IN density on TQ (ie. theta) levels;
                                   !    used in RHOKM so wet density
   rho_dry_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              bl_levels),                                               &
                                   ! IN density on TQ (ie. theta) levels;
                                   !    used in non-turb flux integration
                                   !    so dry density if Lq_mix_bl is true
   dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               bl_levels),                                              &
                                   ! IN DZL(,K) is depth in m of theta
                                   !    level K, i.e. distance from
                                   !    boundary K-1/2 to boundary K+1/2
   rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! IN RDZ(,1) is the reciprocal of
                                   !    the height of level 1, i.e. of
                                   !    the middle of layer 1.  For
                                   !    K > 1, RDZ(,K) is the
                                   !    reciprocal of the vertical
                                   !    distance from level K-1 to
                                   !    level K.
   rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels),                                         &
                                   ! IN RDZ(,1) is the reciprocal of
                                   !       the height of level 1,
                                   !       i.e. of the middle of layer 1
                                   !       For K > 1, RDZ(,K) is the
                                   !       reciprocal of the vertical
                                   !       distance from level K-1 to
                                   !       level K.
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                   ! IN Z_tq(*,K) is height of full
                                   !    level k.
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                                    ! OUT Z_uv(*,K) is height of half
                                    ! level k-1/2.
   rhostar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! IN Surface air density
   u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! IN U on P-grid.
   v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! IN V on P-grid.
   bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                   ! IN A buoyancy parameter for clear
                                   !    air on p,T,q-levels
                                   !    (full levels).
   bq(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),   &
                                   ! IN A buoyancy parameter for clear
                                   !    air on p,T,q-levels
                                   !    (full levels).
   bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                   ! IN A buoyancy parameter for cloudy
                                   !    air on p,T,q-levels
                                   !    (full levels).
   bq_cld(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,          &
          bl_levels),                                                   &
                                   ! IN A buoyancy parameter for cloudy
                                   !    air on p,T,q-levels
                                   !    (full levels).
   bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                   ! IN A grid-box mean buoyancy param
                                   ! on p,T,q-levels (full levels).
   bq_gb(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),&
                                   ! IN A grid-box mean buoyancy param
                                   ! on p,T,q-levels (full levels).
   a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                   ! IN Saturated lapse rate factor
                                   !    on p,T,q-levels (full levels).
   a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                                   ! IN Saturated lapse rate factor
                                   !    on p,T,q-levels (full levels).
   dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                   ! IN Derivative of q_SAT w.r.t. T

  REAL, INTENT(IN) ::                                                   &
   recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end), &
                                   ! IN Reciprocal of the surface
                                   !    Obukhov length over sea (m^-1).
   flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                   ! IN Land fraction on all tiles
  p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! P_HALF(*,K) is pressure at half
                                ! level k-1/2.
   pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! IN Surface pressure (Pascals).

! (f) Atmospheric + any other data not covered so far, incl control.

  REAL, INTENT(IN) ::                                                   &
   rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2,bl_levels),                                                 &
                                    ! IN (LW,SW) rad heating rate (K/s)
    micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                2, bl_levels),                                          &
                           ! Tendencies from microphys within BL levels
                           ! (TL, K/s; QW, kg/kg/s)
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                    ! IN Surface flux buoyancy over
                                    ! density (m^2/s^3)

   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                    ! IN Surface friction velocity
                                    !    (m/s)
   h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),   &
                                    ! IN Blending height used as part
                                    ! of effective roughness scheme
   zh_prev(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                    ! IN boundary layer height from
                                    !    previous timestep
   rib_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                 ! IN  Bulk Richardson number for lowest
                                 ! layer
   sil_orog_land(land_pts),                                             &
                                 ! IN Silhouette area of unresolved
                                 ! orography per unit horizontal area
   zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Height of top of initial
                                 !     parcel ascent
   z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! IN Height of LCL

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
   nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
   L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
   lq_mix_bl

  REAL, INTENT(IN) ::                                                   &
   u_0_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! IN W'ly component of surface
!                                       current (m/s). P grid
   v_0_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! IN S'ly component of surface
!                                       current (m/s). P grid
   ho2r2_orog(land_pts),                                                &
                                   ! IN peak to trough height of 
!                                       unresolved orography
!                                       on land points only (m)
   sd_orog(land_pts),                                                   &
                                   ! IN Standard Deviation of unresolved 
!                                       orography on land points only (m)
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN Effective grid-box roughness
!                                 length for momentum

  INTEGER, INTENT(IN) ::                                                &
   land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.
! (e) Cloud data.
  REAL, INTENT(IN) ::                                                   &
   cf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, bl_levels),  &
                                          ! IN Cloud fraction (decimal).
   qcf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                     ! IN Cloud ice (kg per kg air)
   qcl(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                     ! IN Cloud liquid water
   q(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),    &
                                     ! IN specific humidity
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                     ! IN temperature
   qw(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end, bl_levels),  &
                                   ! IN Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                                   ! IN Ice/liquid water temperature

! INOUT variables
  REAL, INTENT(INOUT) ::                                                &
   zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                                   ! INOUT Height above surface of top
                                   !       of boundary layer (metres).
   fqw(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                   ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                   ! INOUT Exchange coeffs for moisture.
    w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,0:bl_levels), &
    etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,         &
           0:bl_levels),                                                &
   t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                    ! INOUT Standard deviation of
                                    ! turbulent fluctuations of layer 1
                                    ! temperature; for use in
                                    ! initiating convection.
   q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                    ! INOUT Standard deviation of turbulent
                                    !    fluctuations of layer 1
                                    !    humidity; for use in initiating
                                    !    convection.

  REAL, INTENT(INOUT) ::                                                &
   rhokm(pdims_s%i_start:pdims_s%i_end,                                 &
         pdims_s%j_start:pdims_s%j_end ,bl_levels)
!            Exchange coefficients for momentum on P-grid

! INOUT but not used: variables used in the 1A version (TKE-based schemes)
  REAL, INTENT(INOUT) ::                                                &
    e_trb(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,  &
        bl_levels),                                                     &
    tsq_trb(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,&
        bl_levels),                                                     &
    qsq_trb(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,&
        bl_levels),                                                     &
    cov_trb(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,&
        bl_levels),                                                     &
    zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

  LOGICAL, INTENT(INOUT) ::                                             &
   cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! INOUT Logical switch for trade Cu
   l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   ! INOUT Flag to indicate shallow
                                   !     convection

  INTEGER, INTENT(INOUT) ::                                             &
   ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT Top level of initial parcel
                                 !  ascent. Used in convection scheme.

!  Outputs :-
!  (a) Calculated anyway (use STASH space from higher level) :-
 REAL, INTENT(OUT) ::                                                   &
  f_ngstress(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,2:bl_levels),                &
  rhogamu(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
                    ! Counter gradient terms for u
                    ! defined at theta level K-1
  rhogamv(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
                    ! Counter gradient terms for v
                    ! defined at theta level K-1
  tau_fd_x(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end, &
           bl_levels),                                                  &
  tau_fd_y(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end, &
           bl_levels),                                                  &
   bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if stable
                                   !     b.l. diagnosed, 0.0 otherwise.
   bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if Sc over
                                   !     stable surface layer diagnosed,
                                   !     0.0 otherwise.
   bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if well
                                   !     mixed b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer (not over
                                   !     cumulus) diagnosed,
                                   !     0.0 otherwise.
   bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer over cumulus
                                   !     diagnosed, 0.0 otherwise.
   bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if a
                                   !     cumulus capped b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! OUT Indicator set to 1.0 if a
                                   !     Shear-dominated unstable b.l.
                                   !     diagnosed, 0.0 otherwise.

  REAL, INTENT(OUT) ::                                                  &
    zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                   ! OUT Max height of turb mixing
    wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT Convective velocity scale (m/s)
    wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT surface flux of thv (Km/s)
    shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! OUT Shallow Cu diagnostic
                                   !   Indicator set to 1.0 if shallow,
                                   !   0.0 if not shallow or not cumulus
    cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                   ! OUT Indicator for cumulus
                                   !     over steep orography
                                   !   Indicator set to 1.0 if true,
                                   !   0.0 if false. Exclusive.
    we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT a fraction of the timestep
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! OUT Top of decoupled layer

  INTEGER, INTENT(OUT) ::                                               &
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Bottom level of any decoupled
                                   !     turbulently-mixed Sc layer.
    kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! OUT grid-level of SML inversion
    kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                  ! OUT grid-level of DSC inversion
   ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! OUT Number of model levels in the
                                 !     surface-based turbulently mixed
                                 !     layer.

!-2 Genuinely output, needed by other atmospheric routines :-
  REAL, INTENT(OUT) ::                                                  &
    uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                             ! OUT U-component of surface wind stress
                             !     on P-grid
    vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                             ! OUT V-component of surface wind stress
                             !     on P-grid
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-


! Derived local parameters.
  REAL :: lcrcp,ls,lsrcp
  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )
! Parameters also passed to EX_COEF
! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.
  INTEGER ::    k_log_layr
  PARAMETER (k_log_layr=2)
!-----------------------------------------------------------------------
! Constant in TKE diagnostic (K_M = C_TKE * l_m * e^0.5)
  REAL ::      c_tke
  PARAMETER (c_tke=0.5)
!-----------------------------------------------------------------------
!  Workspace :-
  REAL ::                                                               &
   l_int                     ! Length scale for TKE diagnostic

  REAL, ALLOCATABLE ::                                                  &
   visc_bl_h(:,:,:),                                                    &
                        ! OUT: lambda^2*S*FH (on BL_LEVELS)
   visc_bl_h_rho(:,:,:)    ! visc_BL_h on rho levels

  REAL ::                                                               &
   a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                ! Saturated lapse rate factor
                                ! on intermediate levels (half levels).
   a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! Saturated lapse rate factor
                                ! on intermediate levels (half levels).
   bqm(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                ! A buoyancy parameter for clear air
                                ! on intermediate levels (half levels).
   bqm_cld(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,         &
           bl_levels),                                                  &
                                ! A buoyancy parameter for cloudy air
                                ! on intermediate levels (half levels).
   btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                ! A buoyancy parameter for clear air
                                ! on intermediate levels (half levels).
   btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                                ! A buoyancy parameter for cloudy air
                                ! on intermediate levels (half levels).
   cfm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                ! Estimate of cloud fraction
                                ! on intermediate levels (half levels).
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2:bl_levels),                                                   &
                                ! Buoyancy gradient across layer
                                !  interface.
   dbdz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                                ! Buoyancy gradient across layer
                                !  interface, inc gradient adjustment
   dvdzm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels),                                                  &
                                ! Modulus of wind shear.
   ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels), &
                                ! Local Richardson number.
   ri_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                                ! Local Richardson number, inc grad adj
   deltap_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             bl_levels),                                                &
!                                 Difference in pressure between levels
!                                 on UV points
   grad_q_adj(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end),     &
                                ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   grad_t_adj(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   rhokhz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                                ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
   rhokh_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels),                                              &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of heat and
                                ! moisture.
   rhokh_th(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                ! local scheme rhokh on th-levels, 
                                ! index k held on th-level(k-1),  
                                ! same as rhokm
   rhokmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                                ! Non-local turbulent mixing
!                                 coefficient for momentum.
   rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels),                                              &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of momentum.
   weight_1dbl(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end ,bl_levels),                   &
                                ! Weighting applied to 1D BL scheme 
                                ! to blend with Smagorinsky scheme,
                                ! index k held on theta level (k-1)
   weight_1dbl_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,bl_levels),                &
                                ! weight_1dbl interpolated to rho levels
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
                                ! Mixing length for momentum
   elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
   elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                                ! Mixing length for heat (m), 
                                ! held on theta and rho levels, resp.
   fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! stability function for momentum transport
                                ! level 1 value is dummy
   fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! stability function for heat and moisture.
                                ! level 1 value is dummy
   sigma_h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! Standard deviation of subgrid
                                ! orography for sg mixing options (m)

  REAL, DIMENSION (:,:,:), ALLOCATABLE ::                               &
    visc_h_rho                  ! visc_h on rho levels

      ! Terms for non-gradient flux parametrization
      !  (=0 unless using 8C code with FLUX_GRAD=LockWhelan2006)
  REAL ::                                                               &
    ft_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1),                                                 &
                                ! Non-turbulent heat and moisture flux
    fq_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1)          !  (on rho levels, surface flux(K=1)=0)
  REAL ::                                                               &
    rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          2:bl_levels),                                                 &
                                ! f2 and fsc term shape profiles
    rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           2:bl_levels)

  REAL ::                                                               &
    tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                ! Total heat fluxes at inversions
    tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                !
    totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                ! Total moisture fluxes at inversions
    totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                !
    ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                ! Non-turbulent heat and moisture flux
    fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                !    at the base of the DSC layer.

  REAL ::                                                               &
  zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                ! Height above surface of top of
                                !  boundary layer (metres) as
                                !  determined from the local
                                !  Richardson number profile.
  riout(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                ! Gradient Richardson number
                                ! for SCM output
                                ! RIOUT(K) is on theta-level K
  dtldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2:bl_levels),                                                   &
                                ! TL+gz/cp gradient between
                                ! levels K and K-1
  dtldz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                                ! TL+gz/cp gradient between
                                ! levels K and K-1, inc gradient adjust
  dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels)
                                ! QW gradient between

  INTEGER ::                                                            &
   ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                   ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the local
!                                    Richardson number profile.
   ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the parcel ascent.
   sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! Flags for whether discontinuous
   dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! inversions are diagnosed

  LOGICAL ::                                                            &
   unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! Logical switch for unstable
                                 !    surface layer.
   dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! Flag set if decoupled
                                 ! stratocumulus layer found
   coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! Flag to indicate Sc layer weakly
                                 ! coupled to surface (ie weakly
                                 ! decoupled)
   dynamic_bl_diag(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),&
                                 ! Flag to indicate the dynamic 
                                 ! diagnosis (iDynDiag) has 
                                 ! determined the BL type
   topbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! Flag for having reached
                                 ! the top of the turbulently mixed
                                 ! layer.

!  Local scalars :-
  REAL ::                                                               &
   dzu,                                                                 &
              ! Westerly wind shear between levels K+1 and K.
   dzv,                                                                 &
              ! Southerly wind shear between levels K+1 and K.
   lambda_min,                                                          &
              ! Min value of length scale LAMBDA.
   lambdah,                                                             &
              ! Asymptotic mixing length for turbulent transport
              ! of heat/moisture.
   vkz,                                                                 &
              ! Temporary in calculation of ELH.
   f_log,                                                               &
              ! Temporary in calculation of logarithmic correction
   zmaxb_for_dsc,                                                       &
   zmaxt_for_dsc
              ! Max heights to look for DSC cloud base and top

  REAL ::                                                               &
    weight1,                                                            &
    weight2,                                                            &
    weight3,                                                            &
    z_scale,                                                            &
               ! scaling with height
    zpr,                                                                &
               ! z/sigma_h
    slope,                                                              &
               ! subgrid orographic slope
    grcp,                                                               &
               ! G/CP
    dtldzm,                                                             &
               ! TL+gz/cp gradient interpolated to Z_TQ
    dqwdzm     ! QW gradient interpolated to Z_TQ

  INTEGER ::                                                            &
   i,j,                                                                 &
                  ! LOCAL Loop counter (horizontal field index).
   k,ient,                                                              &
                  ! LOCAL Loop counter (vertical level index).
   kp,km,                                                               &
                  ! K+/-1,
   l,                                                                   &
                  ! LOCAL Loop counter for land points
   ntop       ! NTPAR restricted below BL_LEVELS

  INTEGER ::                                                            &
   omp_block,                                                           &
                   ! for open mp blocking 
   jj              
                   ! for indexing over open mp block

  REAL, PARAMETER :: max_abs_obkhov = 1.0e6
                   ! Maximum permitted magnitude of the Obukhov
                   ! length (m).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BDY_EXPL2',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1) Set various diagnostics and switches
!-----------------------------------------------------------------------
  IF ( nl_bl_levels == off .OR.                                         &
       nl_bl_levels > bl_levels ) nl_bl_levels = bl_levels
                                      ! if unset or too large

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
        dynamic_bl_diag(i,j) = .FALSE.
    END DO
  END DO

!------------------------------------------------------------------
!  Initialize weighting applied to 1d BL scheme 
!  (used to blend with 3D Smagorinsky scheme)
!------------------------------------------------------------------
  DO k=1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        weight_1dbl(i,j,k) = 1.0
        weight_1dbl_rho(i,j,k) = 1.0
      END DO
    END DO
  END DO

  grcp = g/cp
!-----------------------------------------------------------------------
! Set surface scaling diagnostics
!-----------------------------------------------------------------------
 ! Obukhov length
  IF (BL_diag%l_oblen) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
!       Limit the magnitude of the Obukhov length to avoid
!       problems with packing.
        BL_diag%oblen(i,j)= u_s(i,j)*u_s(i,j)*u_s(i,j)
        IF ( BL_diag%oblen(i,j) <                                       &
             max_abs_obkhov*ABS(vkman*fb_surf(i,j)) ) THEN
          BL_diag%oblen(i,j)=-BL_diag%oblen(i,j)/(vkman*fb_surf(i,j))
        ELSE
          BL_diag%oblen(i,j)=-SIGN(max_abs_obkhov, fb_surf(i,j))
        ENDIF
      END DO
    END DO
  END IF
! Ustar
  IF (BL_diag%l_ustar) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%ustar(i,j)=u_s(i,j)
      END DO
    END DO
  END IF
! Surface buoyancy flux
  IF (BL_diag%l_wbsurf) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%wbsurf(i,j)=fb_surf(i,j)
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.  Interpolate BT and BQ to half levels and calculate Ri
!-----------------------------------------------------------------------
! DEPENDS ON: btq_int
  CALL btq_int (                                                        &
! IN levels
   bl_levels,                                                           &
! IN fields
   z_tq,z_uv,bq,bt,bq_cld,bt_cld,a_qs,a_dqsdt,cf,                       &
! OUT fields
   bqm,btm,bqm_cld,btm_cld,a_qsm,a_dqsdtm,cfm                           &
    )
!-----------------------------------------------------------------------
! Calculate lapse rates
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1,  weight2, weight3,&
!$OMP& dtldzm, dqwdzm, zpr, dzv, dzu, l, slope)

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      grad_t_adj(i,j) = MIN( max_t_grad,                                &
                             a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dtldz(i,j,k)    = ( tl(i,j,k) - tl(i,j,k-1) )                   &
                                *rdz_charney_grid(i,j,k) + grcp
        dtldz_ga(i,j,k) = dtldz(i,j,k)
        IF ( z_tq(i,j,k) <= zh_prev(i,j) ) THEN
          dtldz_ga(i,j,k) = dtldz_ga(i,j,k) - grad_t_adj(i,j)
        END IF
        dqwdz(i,j,k)    = ( qw(i,j,k) - qw(i,j,k-1) )                   &
                               * rdz_charney_grid(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO


! Local Ri-based calculation of RHOKM and RHOKH:
! Calculate `buoyancy' gradient, DBDZ, on theta-levels
! NOTE: DBDZ(K) is on theta-level K-1

!$OMP DO SCHEDULE(STATIC)
  DO k = 3, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        weight1 = r_rho_levels(i,j,k) -                                 &
                  r_rho_levels(i,j,k-1)
        weight2 = r_theta_levels(i,j,k-1)-                              &
                  r_rho_levels(i,j,k-1)
        weight3 = r_rho_levels(i,j,k) -                                 &
                  r_theta_levels(i,j,k-1)
        dtldzm = weight2 * dtldz(i,j,k)                                 &
               + weight3 * dtldz(i,j,k-1)
        dqwdzm = weight2 * dqwdz(i,j,k)                                 &
               + weight3 * dqwdz(i,j,k-1)
        dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dtldzm +                       &
                          bq_gb(i,j,k-1)*dqwdzm )/weight1
!       ! Now with gradient adjustment
        dtldzm = weight2 * dtldz_ga(i,j,k)                              &
               + weight3 * dtldz_ga(i,j,k-1)
        dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dtldzm +                    &
                             bq_gb(i,j,k-1)*dqwdzm )/weight1
      END DO
    END DO
  END DO
!$OMP END DO

  k = 2
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dtldz(i,j,k) +                   &
                        bq_gb(i,j,k-1)*dqwdz(i,j,k) )
      dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dtldz_ga(i,j,k) +             &
                           bq_gb(i,j,k-1)*dqwdz(i,j,k) )
    END DO
  END DO
!$OMP END DO
!--------------------------------------------------
! Calculate modulus of shear on theta-levels
! dvdzm(k) is on theta-level(k-1)
!--------------------------------------------------
IF (l_subfilter_blend) THEN
! On entry, visc_m is 3D shear(k) on theta-level(k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dvdzm(i,j,k) = visc_m(i,j,k-1)
      END DO
    END DO
  END DO
!$OMP END DO

ELSE

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dzu = u_p(i,j,k) - u_p(i,j,k-1)
        dzv = v_p(i,j,k) - v_p(i,j,k-1)
        dvdzm(i,j,k) = MAX ( 1.0e-12 ,                                  &
                 SQRT(dzu*dzu + dzv*dzv) * rdz(i,j,k)  )
      END DO
    END DO
  END DO
!$OMP END DO

END IF
!-----------------------------------------------------------------------
! 2.1 Orographic enhancement of subgrid mixing
!-----------------------------------------------------------------------
!  Set-up 2D array for standard deviation of subgrid orography.  
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      sigma_h(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO l = 1, land_pts
    j=(land_index(l)-1)/pdims%i_end + 1
    i=land_index(l) - (j-1)*pdims%i_end
    sigma_h(i,j) = MIN( sd_orog(l), 300.0 )
  END DO
!$OMP END DO
!-----------------------------------------------------------------------
!  Enhance resolved shear through unresolved subgrid drainage flows.
!-----------------------------------------------------------------------
  IF (sg_orog_mixing >= 2) THEN

!$OMP DO SCHEDULE(STATIC)
    DO k=2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          IF (sigma_h(i,j) > 1.0 ) THEN
            zpr = z_tq(i,j,k-1)/sigma_h(i,j)
            ! Height dependence, to reduce effect to zero with height
            !   gives z_scale~[1,0.95,0.5,0] at zpr=[0,0.6,1,1.7]
            weight1 = 0.5*( 1.0 - TANH(4.0*(zpr-1.0) ) )

            ! Take slope ~ sd/h_scale for small sd; 
            !            tends to 0.2 for large sd
            slope = 1.0 / SQRT( 25.0 + (h_scale/sigma_h(i,j))**2 )

            dvdzm(i,j,k) = MAX ( dvdzm(i,j,k),                          &
                                 weight1*slope*t_drain*dbdz(i,j,k) )

            IF (k==2 .AND. BL_diag%L_dvdzm)                             &
              BL_diag%dvdzm(i,j,1)=weight1*slope*t_drain*dbdz(i,j,k)

          END IF
        END DO
      END DO
    END DO
!$OMP END DO
  END IF     ! sg_orog_mixing

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ri(i,j,k)    = dbdz(i,j,k)    / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
        ri_ga(i,j,k) = dbdz_ga(i,j,k) / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL 
  IF (BL_diag%l_gradrich) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%gradrich(i,j,k)=ri(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_dbdz) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%dbdz(i,j,k)=dbdz(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_dvdzm) THEN
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%dvdzm(i,j,k)=dvdzm(i,j,k)
        END DO
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
! 3.  Orographic formdrag - distributed drag option
!-----------------------------------------------------------------------
  IF (formdrag ==  explicit_stress) THEN
!------------------------------------------------------------------
!      Set stresses to zero
!------------------------------------------------------------------
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          tau_fd_x(i,j,k) = 0.0
          tau_fd_y(i,j,k) = 0.0
        END DO
      END DO
    END DO
!------------------------------------------------------------------
!      Calculate stress profiles
!------------------------------------------------------------------
! DEPENDS ON: fm_drag
    CALL fm_drag (                                                      &
! IN levels
      land_pts, land_index, bl_levels,                                  &
! IN fields
      u_p, v_p, rho_tq, z_uv, z_tq, z0m_eff_gb, zh_prev, rib_gb,        &
      sil_orog_land,                                                    &
! OUT fields
      tau_fd_x(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               1:bl_levels),                                            &
      tau_fd_y(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               1:bl_levels)                                             &
      )
!------------------------------------------------------------------
!      Orographic stress diagnostics
!------------------------------------------------------------------
    IF (BL_diag%l_ostressx) THEN
      DO k = 1, bl_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
              BL_diag%ostressx(i,j,k)=tau_fd_x(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF (BL_diag%l_ostressy) THEN
      DO k = 1, bl_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
              BL_diag%ostressy(i,j,k)=tau_fd_y(i,j,k)
          END DO
        END DO
      END DO
    END IF

  END IF
!-----------------------------------------------------------------------
! 4. Apply dynamic diagnosis of shear-driven layers.
!-----------------------------------------------------------------------
!       In cases where the parcel ascent continues right through the
!       boundary layer, we diagnose cumulus only where the surface
!       buoyancy flux is sufficiently unstable, testing the ratio of
!       the depth of the inversion to the Obukhov length. A value of
!       1 -- 2 is reasonable for this test and 1.6 is selected, but
!       no great precision is attached to this value. Since this is
!       of importance mainly at sea-points, to avoid complications
!       with coastal tiling, the scheme operates only at points
!       where the land fraction is below 0.5.

  omp_block = pdims%j_end
!$ omp_block = pdims%j_end/omp_get_max_threads()

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, jj, ntop, z_scale)
  IF (iDynDiag == DynDiag_ZL) THEN

! Original version - causes spuriously deep boundary layers if
! BL_LEVELS is >> 3km

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF ( flandg(i,j) < 0.5 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          IF ( -z_uv(i,j,ntop+1) * recip_l_mo_sea(i,j) < 1.6 ) THEN
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
            ! Override the provisional cumulus diagnosis if the
            ! actual surface buoyancy flux indicates stability.
          IF ( fb_surf(i,j) < 0.0 ) THEN
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = 1
            zh(i,j)        = z_uv(i,j,2)
          END IF
        END IF

      END DO
    END DO
!$OMP END DO

  ELSE IF (iDynDiag == DynDiag_ZL_corrn) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

          !------------------------------------------------------------
          ! As original code above, except restrict depth scale to <3km
          ! and, if near-neutral, ignore the BL depth diagnosed by the
          ! adiabatic parcel (ie ZH, NTML) completely.
          !------------------------------------------------------------
        IF ( flandg(i,j) < 0.5 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
          IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 ) THEN
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = 1
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
        END IF
          !-----------------------------------------------------
          ! Override the provisional cumulus diagnosis if the
          ! actual surface buoyancy flux indicates stability,
          ! irrespective of whether over land or sea.
          !-----------------------------------------------------
        IF ( fb_surf(i,j) < 0.0 ) THEN
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          ntml(i,j)      = 1
          zh(i,j)        = z_uv(i,j,2)
        END IF

      END DO
    END DO
!$OMP END DO

  ELSE IF (iDynDiag == DynDiag_ZL_CuOnly) THEN

!$OMP DO SCHEDULE(STATIC) 
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

          !-----------------------------------------------------
          ! Override the provisional cumulus diagnosis if the
          ! actual surface buoyancy flux indicates stability,
          ! irrespective of whether over land or sea.
          !-----------------------------------------------------
        IF ( fb_surf(i,j) < 0.0 ) THEN
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          ntml(i,j)      = 1
          zh(i,j)        = z_uv(i,j,2)
        END IF
          !------------------------------------------------------------
          ! As DynDiag_ZL_corrn but only affects cumulus points and 
          ! points that are entirely sea.  Note that DynDiag_ZL_corrn
          ! has been found to switch off non-local mixing in 
          ! stratocumulus, where surface fluxes are typically small
          !------------------------------------------------------------
        IF ( cumulus(i,j) .AND. flandg(i,j) < 0.01 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
          IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 ) THEN
            dynamic_bl_diag(i,j) = .TRUE.
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
        END IF

      END DO
    END DO
!$OMP END DO

  ELSE IF (iDynDiag == DynDiag_Ribased ) THEN
    !------------------------------------------------------------
    ! As DynDiag_ZL_CuOnly but also allow ZH(Ri) to overrule the 
    ! Cumulus diagnosis
    !------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        topbl(i,j)           = .FALSE.
      END DO
    END DO
!$OMP END DO
    !---------------------------------------------------------------
    !  Loop over levels to find Ri > RiCrit_sharp (=0.25) to find 
    !  level to which Ri really is close to neutral
    !---------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO jj =pdims%j_start, pdims%j_end, omp_block 
      DO k = 2, bl_levels
        DO j = jj, MIN(jj+omp_block-1,pdims%j_end)
          DO i = pdims%i_start, pdims%i_end
            IF ( .NOT.topbl(i,j) .AND.                                    &
              (ri_ga(i,j,k) >  RiCrit_sharp .OR. k >  bl_levels-1) )  THEN
              topbl(i,j) = .TRUE.
              zh_local(i,j) = z_uv(i,j,k)
            END IF
          END DO  ! Loop over points
        END DO  ! Loop over points
      END DO  ! Loop over levels
    END DO
!$OMP END DO 
    !---------------------------------------------------------------
    !  Overrule Cumulus flag where close to neutral BL
    !---------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        !-----------------------------------------------------
        ! Override the provisional cumulus diagnosis if the
        ! actual surface buoyancy flux indicates stability,
        ! irrespective of whether over land or sea.
        !-----------------------------------------------------
        IF ( fb_surf(i,j) < 0.0 ) THEN
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          ntml(i,j)      = 1
          zh(i,j)        = z_uv(i,j,2)
        END IF

        IF ( cumulus(i,j) .AND. flandg(i,j) < 0.01 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
          IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 .OR.                   &
                  ! - ZH/L indicates BL close to neutral
                zh_local(i,j) > zh(i,j)+zhloc_depth_fac*(z_scale-zh(i,j))&
                  ! ZH(Ri>RiCrit) more than zhloc_depth_fac up the 
                  ! cloud layer, indicating significant shear disruption
              ) THEN
            dynamic_bl_diag(i,j) = .TRUE.
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
        END IF

      END DO
    END DO
!$OMP END DO


  END IF  ! tests on iDynDiag


!$OMP END PARALLEL 
!-----------------------------------------------------------------------
! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------
! 5.1  Calculate the non-local terms and diffusion coefficients
!-----------------------------------------------------------------------
! Set NTML_NL to NTML as passed in from initial diagnosis routine
!-----------------------------------------------------------------------
  zmaxb_for_dsc = 2500.0
  zmaxt_for_dsc = 3000.0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ntml_nl(i,j) = ntml(i,j)
    END DO
  END DO

  IF (nl_bl_levels < bl_levels) THEN
        ! Set to huge value to make if-test in KMKHZ redundent
    zmaxb_for_dsc = 1.0e10
    zmaxt_for_dsc = zmaxb_for_dsc
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF ( ntml_nl(i,j) > nl_bl_levels-1 ) THEN
          ntml_nl(i,j) = nl_bl_levels-1
          zh(i,j)      = z_uv(i,j,ntml_nl(i,j)+1)
        END IF
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
! Initialise non-local K and fluxes to zero; necessary for levels
! above NL_BL_LEVELS
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ftl(i,j,k) = 0.0
        fqw(i,j,k) = 0.0
        rhokmz(i,j,k) = 0.0
        rhokhz(i,j,k) = 0.0
        rhokm_top(i,j,k) = 0.0
        rhokh_top(i,j,k) = 0.0
        f_ngstress(i,j,k) = 0.0
      END DO
    END DO
  END DO
      ! Initialise Lock-Whelan non-gradient terms to zero
      ! Only calculated in KMKHZ8C
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ft_nt(i,j,k)  = 0.0
        fq_nt(i,j,k)  = 0.0
      END DO
    END DO
  END DO
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rhof2(i,j,k)  = 0.0
        rhofsc(i,j,k) = 0.0
      END DO
    END DO
  END DO
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      tothf_zh(i,j)   = 0.0
      tothf_zhsc(i,j) = 0.0
      totqf_zh(i,j)   = 0.0
      totqf_zhsc(i,j) = 0.0
      ft_nt_dscb(i,j) = 0.0
      fq_nt_dscb(i,j) = 0.0
    END DO
  END DO

  IF (l_subfilter_vert .AND. .NOT. l_subfilter_blend) THEN
    non_local_bl = off
  END IF

  IF (non_local_bl == on) THEN

! DEPENDS ON: kmkhz
    CALL kmkhz (                                                        &
! IN levels/switches
       nl_bl_levels,lq_mix_bl,BL_diag, nSCMDpkgs,L_SCMDiags,            &
! IN fields
       p_theta_levels,rho_tq,rho_uv,rho_dry_tq,t,q,qcl,qcf,cf,qw,tl,    &
       dzl_charney,rdz_charney_grid,z_tq,z_uv,                          &
       rad_hr,micro_tends,                                              &
       bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,         &
       u_s,fb_surf,rhostar,ntpar,nlcl,zh_prev,                          &
       zhpar,z_lcl,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,               &
! INOUT fields
       ftl,fqw,zh,cumulus,ntml_nl,w,etadot,t1_sd,q1_sd,                 &
! OUT fields
       rhokmz,rhokhz,rhokm_top,rhokh_top,zhsc,                          &
       unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,                  &
       ntdsc,nbdsc,                                                     &
       f_ngstress(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  2:nl_bl_levels),                                      &
       grad_t_adj, grad_q_adj,                                          &
       rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb,             &
       tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc,                      &
       kent, we_lim, t_frac, zrzi,                                      &
       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                       &
       )

  ELSE   ! not NON_LOCAL_BL

         !-------------------------------------------------------------
         ! Set all variables from the non-local scheme to zero or "off"
         !  - reset all fluxes and K's arising from the non-local scheme
         !-------------------------------------------------------------

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
          ! surface mixed layer
        unstable(i,j) = (fb_surf(i,j) >  0.0)
        cumulus(i,j) = .FALSE.
        l_shallow(i,j) = .FALSE.
        sml_disc_inv(i,j) = 0
        ntpar(i,j)   = 0
        ntml_nl(i,j) = -1    ! to ensure correct diagnostics
        zh(i,j)      = 0.0
        grad_t_adj(i,j) = 0.0
        grad_q_adj(i,j) = 0.0
          ! decoupled mixed layer
        dsc(i,j)     = .FALSE.
        dsc_disc_inv(i,j) = 0
        ntdsc(i,j)   = 0
        nbdsc(i,j)   = 0
        zhsc(i,j)    = 0.0
        coupled(i,j) = .FALSE.
          ! entrainment variables for non-local tracer mixing
        kent(i,j) = 2
        kent_dsc(i,j) = 2
        DO ient = 1, 3
          t_frac(i,j,ient) = 0.0
          zrzi(i,j,ient)   = 0.0
          we_lim(i,j,ient) = 0.0
          t_frac_dsc(i,j,ient) = 0.0
          zrzi_dsc(i,j,ient)   = 0.0
          we_lim_dsc(i,j,ient) = 0.0
        END DO
      END DO
    END DO

  END IF  ! test on NON_LOCAL_BL

!-----------------------------------------------------------------------
! 5.1  Call local coeff calculation for levels 2 to bl_levels
!-----------------------------------------------------------------------
! DEPENDS ON: ex_coef
  CALL ex_coef (                                                        &
! IN levels/logicals
   bl_levels,k_log_layr,lq_mix_bl,l_subfilter_vert,l_subfilter_horiz,   &
   l_subfilter_blend,rneutml_sq,delta_smag,nSCMDpkgs,L_SCMDiags,BL_diag,&
! IN fields
   sigma_h,flandg,dbdz,dvdzm,ri,rho_tq,z_uv,z_tq,z0m_eff_gb,            &
   h_blend_orog,ntpar,ntml_nl,ntdsc,nbdsc,u_p,v_p,u_s,fb_surf,qw,tl,    &
! IN/OUT fields
   cumulus,weight_1dbl,                                                 &
! OUT fields
   lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm,rhokh_th,       &
   fm_3d,fh_3d                                                          &
   )


! interpolate rhokh_th to rho levels 2 to bl_levels

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j, k,     &
!$OMP& weight1, weight2, weight3, lambdah, z_scale, vkz, f_log)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        weight1 = r_theta_levels(i,j,k) -                               &
                r_theta_levels(i,j, k-1)
        weight2 = r_theta_levels(i,j,k) -                               &
                r_rho_levels(i,j,k)
        weight3 = r_rho_levels(i,j,k) -                                 &
                r_theta_levels(i,j,k-1)
        IF ( k  ==  bl_levels ) THEN
                ! assume rhokh_th(BL_LEVELS+1) is zero
          rhokh(i,j,k) = ( weight2/weight1 ) * rhokh_th(i,j,k)
          IF (l_subfilter_blend) weight_1dbl_rho(i,j,k) =               &
                                 (weight2/weight1) * weight_1dbl(i,j,k)
        ELSE
          rhokh(i,j,k) = (weight3/weight1) * rhokh_th(i,j,k+1)          &
                       + (weight2/weight1) * rhokh_th(i,j,k)
          IF (l_subfilter_blend) weight_1dbl_rho(i,j,k) =               &
                                (weight3/weight1)*weight_1dbl(i,j,k+1)  &
                              + (weight2/weight1)*weight_1dbl(i,j,k)
        END IF

        IF (local_fa == free_trop_layers) THEN
          ! elh already included in rhokh_th so no need to calculate 
          ! here, but interpolate elh separately for diagnostic
          IF (BL_diag%L_elh3D) THEN
            IF ( k  ==  bl_levels ) THEN
              ! assume rhokh_th(BL_LEVELS+1) is zero
              elh_rho(i,j,k) = ( weight2/weight1 ) * elh(i,j,k)
            ELSE
              elh_rho(i,j,k) =                                          &
                weight3/weight1 *                                       &
                        elh(i,j,k+1)                                    &
               +weight2/weight1 *                                       &
                        elh(i,j,k)
            END IF
            BL_diag%elh3D(i,j,k)=elh_rho(i,j,k)
          END IF
        ELSE
          IF ((sbl_op/=Equilibrium_SBL).OR.(fb_surf(i,j) >  0.0)) THEN
!-------------------------------------------------------------------
!  Include mixing length, ELH, in RHOKH.
!  Code moved from EX_COEF to avoid interpolation
!-------------------------------------------------------------------
          IF ( k >= ntml_local(i,j)+2 .AND. l_full_lambdas .AND.        &
               local_fa == to_sharp_across_1km ) THEN
            ! Assuming only LOCAL_FA = "to_sharp_across_1km" option
            ! will have L_FULL_LAMBDAS.  
            ! If other LOCAL_FA options are coded here then
            ! changes must be included in section 2.1 of ex_coef
            IF (l_rp2) THEN
              lambdah = MAX ( lambda_min , par_mezcla*zh_local(i,j) )
            ELSE
              lambdah = MAX ( lambda_min , 0.15*zh_local(i,j) )
            END IF
            z_scale = 1000.0
            weight1 = 0.5*( 1.0 -                                       &
                        TANH(3.*((z_uv(i,j,k)/z_scale )-1.0) ) )
            lambdah = lambdah * weight1                                 &
                         + lambda_min*( 1.0 -  weight1)
!           ! no need to do log profile correction as klog_layr eq 2
            vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
            elh_rho(i,j,k) = vkz / (1.0 + vkz/lambdah )
          END IF
! Reinstate UKV drainage flow bug here, where lambdah was not enhanced 
! as intended (and as was done in ex_coef)!
          IF (sg_orog_mixing == 3) THEN
            IF (l_rp2) THEN
              lambdah = MAX ( lambda_min , par_mezcla*zh_local(i,j) )
            ELSE
              lambdah = MAX ( lambda_min , 0.15*zh_local(i,j) )
            END IF
            IF (k >= ntml_local(i,j)+2 .AND. .NOT.l_full_lambdas) THEN
              lambdah = lambda_min
            END IF
            IF (k <= k_log_layr) THEN
              vkz   = vkman * ( z_tq(i,j,k) - z_tq(i,j,k-1) )
              f_log = LOG( ( z_tq(i,j,k) + z0m_eff_gb(i,j)   ) /        &
                           ( z_tq(i,j,k-1) + z0m_eff_gb(i,j) ) )
              elh_rho(i,j,k) = vkz / ( f_log + vkz / lambdah )
            ELSE
              vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
              elh_rho(i,j,k) = vkz / (1.0 + vkz/lambdah )
            END IF
          END IF
! End of UKV bug!

          IF (BL_diag%L_elh3D) BL_diag%elh3D(i,j,k)=elh_rho(i,j,k)

          rhokh(i,j,k) = elh_rho(i,j,k) * rhokh(i,j,k)

         END IF   ! test on sbl_op
        END IF   ! test on local_fa = free_trop_layers

            ! Finally multiply RHOKH by dry density
        IF (lq_mix_bl)                                                  &
           rhokh(i,j,k) = rho_uv(i,j,k) * rhokh(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO 

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
        !----------------------------------------------------------
        ! Use local NTML if significantly higher (to allow for
        ! local averaging) than the non-local or if the non-local
        ! is on the ground (=1)
        !----------------------------------------------------------
      IF ( .NOT.cumulus(i,j) .AND.                                      &
                ( ntml_local(i,j)  >   ntml_nl(i,j)+1                   &
                  .OR. ntml_nl(i,j)  ==  1 )            ) THEN
        ntml(i,j) = ntml_local(i,j)
        sml_disc_inv(i,j) = 0   ! reset flag for subgrid inversion
      ELSE
        ntml(i,j) = ntml_nl(i,j)
      END IF
        !----------------------------------------------------------
        ! If local NTML is higher than NTDSC then ignore DSC layer
        ! for diagnostics but keep mixing associated with it
        !----------------------------------------------------------
      IF ( ntml_local(i,j)  >   ntdsc(i,j)+1 ) THEN
        dsc_disc_inv(i,j) = 0
        ntdsc(i,j) = 0
        nbdsc(i,j) = 0
        zhsc(i,j)  = 0.0
        dsc(i,j)   = .FALSE.
        coupled(i,j) = .FALSE.
      END IF
    END DO
  END DO

! Calculate max of two coeffs
! DEPENDS ON: kmkh
  CALL kmkh (                                                           &
! IN data
   bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                              &
   ntml,cumulus,ntdsc,dsc,sml_disc_inv,dsc_disc_inv,                    &
   weight_1dbl, weight_1dbl_rho,                                        &
! INOUT data
   rhokm,rhokh,rhokmz(1,1,2),rhokhz(1,1,2),                             &
   rhokm_top(1,1,2),rhokh_top(1,1,2)                                    &
   )

  IF (BL_diag%l_weight1d) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%weight1d(i,j,k)=weight_1dbl(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_rhokm) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhokm(i,j,k)=rhokm(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_rhokh) THEN
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          BL_diag%rhokh(i,j,k)=rhokh(i,j,k)
        END DO
      END DO
    END DO
  END IF

! Calculation of TKE diagnostic.
! Stored on theta-levels with TKE(K) on theta-level(k-1),
! consistent with RHOKM(K), RI(K), etc.
! The K=1 value could be set to a diagnosed surface value (eg as a
! function of ustar, wstar) but is currently just set to zero

  IF (BL_diag%l_tke) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%tke(i,j,1) = 0.0
      END DO
    END DO

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(z_tq, zh, zhsc, &
!$OMP& z0m_eff_gb, BL_diag, rho_tq, rhokm, bl_levels, pdims,             &
!$OMP& fb_surf, elm, dsc)  PRIVATE(i, j, k, vkz, l_int)
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          vkz = vkman * ( z_tq(i,j,k-1) + z0m_eff_gb(i,j) )

! neutral or stable bl scale (from ex_coef)
          l_int= elm(i,j,k)

          IF ( fb_surf(i,j) > 0.0 ) THEN
! convective bl length scale

            IF (z_tq(i,j,k-1) < zh(i,j)) THEN
              l_int= MAX( l_int,                                        &
                          (zh(i,j)-z_tq(i,j,k-1)) /                     &
                          ( 1.0 +  (zh(i,j)-z_tq(i,j,k-1))/vkz ) )
            END IF      ! z < zh

          END IF   ! FB_SURF(I,j) > 0.0

! decoupled cloud layer
          IF ( dsc(i,j)) THEN
            IF ( z_tq(i,j,k-1) < zhsc(i,j) .AND.                        &
                 z_tq(i,j,k-1) > BL_diag%dscbase(i,j) ) THEN
                 ! Take max of SML and DSC length scales
              l_int= MAX( l_int,                                        &
                          (zhsc(i,j)-z_tq(i,j,k-1)) /                   &
                  (1.0 +  (zhsc(i,j)-z_tq(i,j,k-1))/                    &
                          (z_tq(i,j,k-1)-BL_diag%dscbase(i,j)) )  )

            END IF
          END IF

            ! Don't let l_int get very small
          l_int = MAX( l_int, 5.0 )

            ! TKE diagnostic - taking 5 m2/s2 as a suitable maximum
            !  - large values typically generated when RHOKM_local is
            !    large (unstable Ri) but (currently neutral) l_int is
            !    small
          BL_diag%tke(i,j,k)=  MIN( 5.0,                                &
                ( rhokm(i,j,k) / (rho_tq(i,j,k-1)*c_tke*l_int) )**2 )

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF    ! BL_diag%L_tke

  IF (l_subfilter_blend) THEN
!   ! Blended diffusion coefficients now held in rhokm and rhokh
!   ! so copy to visc_m,h for horizontal diffusion too.
!   ! Need to interpolate rhokh back to theta levels for visc_h

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1, weight2,      &
!$OMP& weight3)

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          visc_m(i,j,k) = rhokm(i,j,k+1)/rho_tq(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels-1
      DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        weight1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
        weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
        weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
        visc_h(i,j,k) = ( weight2 * (rhokh(i,j,k+1)/rho_uv(i,j,k+1))   &
                        + weight3 * (rhokh(i,j,k)  /rho_uv(i,j,k)  ) ) &
                                         / weight1
      END DO
      END DO
    END DO
!$OMP END DO
    k = 1
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_h(i,j,k) = rhokh_th(i,j,k+1)/rho_tq(i,j,k)
      END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL 
    
  ELSE IF (l_subfilter_horiz .OR. l_subfilter_vert) THEN

    ! visc_m,h on IN are just S and visc_m,h(k) are co-located with w(k)
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          visc_m(i,j,k) = visc_m(i,j,k)*rneutml_sq(i,j,k)
          visc_h(i,j,k) = visc_h(i,j,k)*rneutml_sq(i,j,k)
       END DO
      END DO
    END DO

    DO k = 1, bl_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          ! stability functions are indexed with Ri, fm(k) on w(k-1)
          visc_h(i,j,k) = visc_h(i,j,k)*fh_3d(i,j,k+1)
          visc_m(i,j,k) = visc_m(i,j,k)*fm_3d(i,j,k+1)
! APL why apply this cap here for implicit vertical diffusion?
! (also applied in atm_step_phys_init for horiz diffn, that actually needs it)?
          visc_h(i,j,k) = MIN(visc_h(i,j,k),max_diff(i,j))
          visc_m(i,j,k) = MIN(visc_m(i,j,k),max_diff(i,j))
        END DO
      END DO
    END DO

! visc_m and visc _h are now lambda^2*S*FM and lambda^2*S*FH

  IF (l_subfilter_vert) THEN

! visc_h_rho(k) is held on rho(k), same as BL's rhokh
 ALLOCATE (visc_h_rho(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end, bl_levels))

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1, weight2,      &
!$OMP& weight3)

!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          weight1 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)
          IF ( k  ==  bl_levels ) THEN
! assume visc_h(bl_levels) is zero (Ri and thence f_h not defined)
            visc_h_rho(i,j,k) = ( weight2/weight1 ) * visc_h(i,j,k-1)
          ELSE
            visc_h_rho(i,j,k) = ( weight3/weight1 ) * visc_h(i,j,k)      &
                              + ( weight2/weight1 ) * visc_h(i,j,k-1)
          END IF
        END DO
      END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL 

! Overwrite the diffusion coefficients from the local BL scheme
!(RHOKM and RHOKH) with those obtained from the Smagorinsky scheme

    DO k = 2, bl_levels
      IF (k >= turb_startlev_vert .AND.                                 &
          k <= turb_endlev_vert) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = visc_m(i,j,k-1)*rho_tq(i,j,k-1)
            rhokh(i,j,k) = visc_h_rho(i,j,k)*rho_uv(i,j,k)
          END DO
        END DO
      ELSE
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = 0.0
            rhokh(i,j,k) = 0.0
          END DO
        END DO
      END IF
    END DO

    DEALLOCATE (visc_h_rho)

  END IF ! L_subfilter_vert
  END IF ! L_subfilter_horiz or L_subfilter_vert or L_subfilter_blend

!-----------------------------------------------------------------------
! Diagnose boundary layer type.
!      Seven different types are considered:
!      1 - Stable b.l.
!      2 - Stratocumulus over a stable surface layer.
!      3 - Well mixed buoyancy-driven b.l. (possibly with stratocumulus)
!      4 - Decoupled stratocumulus (not over cumulus).
!      5 - Decoupled stratocumulus over cumulus.
!      6 - Cumulus capped b.l.
!      7 - Shear-dominated unstable b.l.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the depth diagnostics

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

!$OMP DO SCHEDULE(STATIC) 
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! Top of surface mixed layer (Ksurf profile)
      IF (BL_diag%l_smltop) THEN
        BL_diag%smltop(i,j)=zh(i,j)
      END IF
! Top of decoupled stratocu layer
      IF (BL_diag%l_dsctop) THEN
        BL_diag%dsctop(i,j)=zhsc(i,j)
      END IF
! Height of diagnosis parcel top
      IF (BL_diag%l_zhpar) THEN
        BL_diag%zhpar(i,j)=zhpar(i,j)
      END IF
! Max height of BL turbulent mixing 
      zht(i,j) = MAX( zh(i,j) , zhsc(i,j) )
      IF ( ntml(i,j)  >   ntml_nl(i,j) ) THEN
        ! Higher local K allowed so reset ZH, ZHT diagnostics
        zh(i,j)  = MAX( zh(i,j) , zh_local(i,j) )
        zht(i,j) = MAX( zht(i,j), zh_local(i,j) )
      END IF

      bl_type_1(i,j) = 0.0
      bl_type_2(i,j) = 0.0
      bl_type_3(i,j) = 0.0
      bl_type_4(i,j) = 0.0
      bl_type_5(i,j) = 0.0
      bl_type_6(i,j) = 0.0
      bl_type_7(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (dynamic_bl_diag(i,j)) THEN
!       ! shear-dominated, via iDynDiag option
        bl_type_7(i,j) = 1.0
      ELSE
        IF (.NOT.unstable(i,j) .AND. .NOT.dsc(i,j) .AND.                &
            .NOT.cumulus(i,j)) THEN
!         ! Stable b.l.
          bl_type_1(i,j) = 1.0
        ELSE IF (.NOT.unstable(i,j) .AND. dsc(i,j) .AND.                &
                 .NOT.cumulus(i,j)) THEN
!         ! Stratocumulus over a stable surface layer
          bl_type_2(i,j) = 1.0
        ELSE IF (unstable(i,j) .AND. .NOT.cumulus(i,j) .AND.            &
                .NOT.dsc(i,j) ) THEN
!         ! Well mixed b.l. (possibly with stratocumulus)
          IF ( ntml(i,j)  >   ntml_nl(i,j) ) THEN
            ! shear-dominated - currently identified
            ! by local NTML overriding non-local
            bl_type_7(i,j) = 1.0
          ELSE
            ! buoyancy-dominated
            bl_type_3(i,j) = 1.0
          END IF
        ELSE IF (unstable(i,j) .AND. dsc(i,j) .AND.                     &
                                        .NOT.cumulus(i,j)) THEN
!         ! Decoupled stratocumulus (not over cumulus)
          bl_type_4(i,j) = 1.0
        ELSE IF (dsc(i,j) .AND. cumulus(i,j)) THEN
!         ! Decoupled stratocumulus over cumulus
          bl_type_5(i,j) = 1.0
        ELSE IF (.NOT.dsc(i,j) .AND. cumulus(i,j)) THEN
!         ! Cumulus capped b.l.
          bl_type_6(i,j) = 1.0
        END IF
      END IF

    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL
!-----------------------------------------------------------------------
! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------
! DEPENDS ON: ex_flux_tq
  CALL ex_flux_tq (                                                     &
! IN levels etc
    bl_levels,nSCMDpkgs,L_SCMDiags,                                     &
! IN fields
    tl,qw,rdz_charney_grid, rhokh, rhokhz(1,1,2),grad_t_adj,grad_q_adj, &
    rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb, tothf_zh,      &
    tothf_zhsc, totqf_zh, totqf_zhsc, weight_1dbl_rho,                  &
    ntml_nl, ntdsc, nbdsc,                                              &
! INOUT fields
    ftl,fqw                                                             &
    )


  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        f_ngstress(i,j,k) = weight_1dbl(i,j,k) * f_ngstress(i,j,k)
      END DO
    END DO
  END DO


!-----------------------------------------------------------------------
! 5.6.1 Calculate explicit surface fluxes of U and V on
!       P-grid for convection scheme
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      uw0(i,j) = -rhokm(i,j,1) *                                        &
                        ( u_p(i,j,1) - u_0_p(i,j) )
      vw0(i,j) = -rhokm(i,j,1) *                                        &
                        ( v_p(i,j,1) - v_0_p(i,j) )
    END DO
  END DO
!-----------------------------------------------------------------------
! 5.7 Set NTML to max number of turbulently mixed layers
!      Calculate quantities to pass to convection scheme.
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      wstar(i,j) = 0.0
      wthvs(i,j) = 0.0
      cu_over_orog(i,j) = 0.0
      IF ( cumulus(i,j) ) THEN
        IF ( fb_surf(i,j)  >   0.0 ) THEN
          wstar(i,j) = ( zh(i,j)*fb_surf(i,j) )**(1.0/3.0)
          wthvs(i,j) = fb_surf(i,j) / ( g * bt(i,j,1) )
        END IF
        wstar(i,j) = MAX( 0.1, wstar(i,j) )
        IF (.NOT. l_param_conv) THEN
          ntml(i,j) = MAX( 2, ntml_nl(i,j) - 1 )
        END IF
      ELSE
        ntml(i,j) = MAX( ntml_nl(i,j) , ntdsc(i,j) )
      END IF
! Limit explicitly calculated surface stresses
! to a physically plausible level.
      IF ( uw0(i,j)  >=  5.0 ) THEN
        uw0(i,j) =  5.0
      ELSE IF ( uw0(i,j)  <=  -5.0 ) THEN
        uw0(i,j) = -5.0
      END IF
      IF ( vw0(i,j)  >=  5.0 ) THEN
        vw0(i,j) =  5.0
      ELSE IF ( vw0(i,j)  <=  -5.0 ) THEN
        vw0(i,j) = -5.0
      END IF
      IF (BL_diag%l_wstar .AND. (fb_surf(i,j) >0.0))  THEN
        BL_diag%wstar(i,j)= (zh(i,j)*fb_surf(i,j))**(1.0/3.0)
      END IF
    END DO
  END DO

  IF (l_param_conv) THEN

! Check for CUMULUS having been diagnosed over steep orography.
! Reset to false but keep NTML at NLCL (though decrease by 2 so that
! coupling between BL and convection scheme can be maintained).
! Reset type diagnostics.

    DO l = 1, land_pts
      j=(land_index(l)-1)/pdims%i_end + 1
      i=land_index(l) - (j-1)*pdims%i_end
      IF (cumulus(i,j) .AND. ho2r2_orog(l)  >   900.0) THEN
        cumulus(i,j) = .FALSE.
        l_shallow(i,j) = .FALSE.
        bl_type_5(i,j) = 0.0
        bl_type_6(i,j) = 0.0
        cu_over_orog(i,j) = 1.0
        IF (ntml(i,j)  >=  3) ntml(i,j) = ntml(i,j) - 2
      END IF
    END DO

! Check that CUMULUS and L_SHALLOW are still consistent

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF ( .NOT. cumulus(i,j) ) l_shallow(i,j) = .FALSE.
      END DO
    END DO

  END IF    ! (l_param_conv)
!-----------------------------------------------------------------------
!     Set shallow convection diagnostic: 1.0 if L_SHALLOW (and CUMULUS)
!                                        0.0 if .NOT. CUMULUS
!-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( cumulus(i,j) .AND. l_shallow(i,j) ) THEN
        shallowc(i,j) = 1.0
      ELSE
        shallowc(i,j) = 0.0
      END IF
    END DO
  END DO

  IF (lhook) CALL dr_hook('BDY_EXPL2',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_expl2
