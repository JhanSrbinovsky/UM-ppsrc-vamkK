
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE EXCF_NL ----------------------------------------------
!
!  Purpose: To calculate non-local exchange coefficients,
!           entrainment parametrization and non-gradient flux terms.
!
!  Programming standard: UMDP3
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE excf_nl (                                                    &
! IN levels/switches
 bl_levels,nSCMDpkgs,L_SCMDiags,                                        &
! IN fields
 rdz,z_uv,z_tq,rho_uv,rho_tq,rhostar_gb,v_s,fb_surf,zhpar,              &
 btm,bqm,btm_cld,bqm_cld,cf_sml,cf_dsc,                                 &
 bflux_surf,bflux_surf_sat,zeta_s,svl_diff_frac,                        &
 df_top_over_cp,zeta_r,bt_top,btt_top,btc_top,                          &
 db_top,db_top_cld,chi_s_top,br_fback,                                  &
 df_dsct_over_cp,zeta_r_dsc,bt_dsct,btt_dsct,                           &
 db_dsct,db_dsct_cld,chi_s_dsct,br_fback_dsc,                           &
 db_ga_dry,db_noga_dry,db_ga_cld,db_noga_cld,                           &
 dsl_sml,dqw_sml,dsl_dsc,dqw_dsc,ft_nt, fq_nt,                          &
! INOUT fields
 dsc,cumulus,coupled,ntml,zh,zhsc,dscdepth,ntdsc,zc,zc_dsc,             &
 ft_nt_zh, ft_nt_zhsc, fq_nt_zh, fq_nt_zhsc,                            &
! OUT fields
 rhokm, rhokh, rhokm_top, rhokh_top,                                    &
 rhokh_top_ent, rhokh_dsct_ent, rhokh_surf_ent,                         &
 rhof2,rhofsc,f_ngstress,zdsc_base,nbdsc                                &
)
  USE atm_fields_bounds_mod, ONLY: pdims, rkmdims
  USE atmos_constants_mod, ONLY: vkman
  USE bl_option_mod, ONLY:                                              &
      on, off, dec_thres_cloud, ng_stress,                              &
      BrownGrant97, BrownGrant97_limited, flux_grad, Locketal2000,      &
      HoltBov1993, LockWhelan2006, entr_smooth_dec,                     &
      Kprof_cu, klcl_entr, max_cu_depth
  USE earth_constants_mod, ONLY: g
  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE stochastic_physics_run_mod, ONLY: l_rp2 , a_ent_1_rp, g1_rp
  IMPLICIT NONE

! IN fields
  INTEGER, INTENT(IN) ::                                                &
   bl_levels
                   ! IN maximum number of boundary layer levels

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
    nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
    L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  bl_levels), INTENT(IN) ::                             &
   rdz,                                                                 &
                                ! IN Reciprocal of distance between
                                !    T,q-levels (m^-1). 1/RDZ(,K) is
                                !    the vertical distance from level
                                !    K-1 to level K, except that for
                                !    K=1 it is just the height of the
                                !    lowest atmospheric level.
   z_uv,                                                                &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_UV(*,K) is the height
                                !    of the k-th u,v-level (half level
                                !    k-1/2) above the surface
   z_tq,                                                                &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, Z_TQ(*,K) is the height
                                !    of the k-th T,q-level (full level
                                !    k) above the surface
   rho_uv,                                                              &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_UV(*,K) is the
                                !    density at the k-th u,v-level
                                !    above the surface
   rho_tq,                                                              &
                                ! IN For a vertically staggered grid
                                !    with a u,v-level first above the
                                !    surface, RHO_TQ(*,K) is the
                                !    density of the k-th T,q-level
                                !    above the surface
   bqm,                                                                 &
                                ! IN Buoyancy parameters for clear and
   btm,                                                                 &
                                !    cloudy air on half levels
   bqm_cld,                                                             &
                                !    (*,K) elements are k+1/2 values
   btm_cld                  !

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
        INTENT(IN) ::                                                   &
   rhostar_gb,                                                          &
                                ! IN Surface density (kg/m3)
   v_s,                                                                 &
                                ! IN Surface friction velocity (m/s).
   fb_surf,                                                             &
                                ! IN Buoyancy flux at the surface over
                                !    density (m^2/s^3).
   zhpar,                                                               &
                                ! IN Height of top of NTPAR
                                !    NOTE: CAN BE ABOVE BL_LEVELS-1
   bflux_surf,                                                          &
                                ! IN Surface buoyancy flux (kg/m/s^3).
   bflux_surf_sat,                                                      &
                                ! IN Saturated-air surface buoyancy
                                !    flux.
   db_top,                                                              &
                                ! IN Buoyancy jump across the top of
                                !    the SML (m/s^2).
   df_top_over_cp,                                                      &
                                ! IN Radiative flux change at cloud top
                                !    divided by c_P (K.kg/m^2/s).
   bt_top,                                                              &
                                ! IN Buoyancy parameter at the top of
                                !    the b.l. (m/s^2/K).
   btt_top,                                                             &
                                ! IN In-cloud buoyancy parameter at
                                !    the top of the b.l. (m/s^2/K).
   btc_top,                                                             &
                                ! IN Cloud fraction weighted buoyancy
                                !    parameter at the top of the b.l.
   db_top_cld,                                                          &
                                ! IN In-cloud buoyancy jump at the
                                !    top of the b.l. (m/s^2).
   chi_s_top,                                                           &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the b.l.
   zeta_s,                                                              &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for surface forced
                                !    entrainment term.
   zeta_r,                                                              &
                                ! IN Non-cloudy fraction of mixing
                                !    layer for cloud top radiative
                                !    cooling entrainment term.
   db_dsct,                                                             &
                                ! IN Buoyancy jump across the top of
                                !    the DSC layer (m/s^2).
   df_dsct_over_cp,                                                     &
                                ! IN Radiative flux change at DSC top
                                !    divided by c_P (K.kg/m^2/s).
   bt_dsct,                                                             &
                                ! IN Buoyancy parameters at the top of
   btt_dsct,                                                            &
                                !    the DSC layer (m/s^2/K)
   db_dsct_cld,                                                         &
                                ! IN In-cloud buoyancy jump at the
                                !    top of the DSC (m/s^2).
   chi_s_dsct,                                                          &
                                ! IN Mixing fraction of just saturated
                                !    mixture at top of the DSC
   zeta_r_dsc,                                                          &
                                ! IN Non-cloudy fraction of DSC
                                !    for cloud top radiative
                                !    cooling entrainment term.
   br_fback,                                                            &
   br_fback_dsc,                                                        &
                                ! IN Weights for degree of buoyancy
                                !    reversal feedback
   dqw_sml,                                                             &
                                ! IN QW change across SML disc inv
   dsl_sml,                                                             &
                                ! IN SL change across SML disc inv
   dqw_dsc,                                                             &
                                ! IN QW change across DSC disc inv
   dsl_dsc,                                                             &
                                ! IN SL change across DSC disc inv
   cf_sml,                                                              &
                                ! IN cloud fraction of SML
   cf_dsc,                                                              &
                                ! IN cloud fraction of DSC layer
   svl_diff_frac                ! IN Fractional svl decoupling difference 

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  2:bl_levels), INTENT(IN) ::                           &
   db_ga_dry,                                                           &
                                ! IN Cloudy and cloud-free buoyancy
   db_noga_dry,                                                         &
                                !    jumps for flux integral
   db_ga_cld,                                                           &
                                !    calculation (m/s2):
   db_noga_cld

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  bl_levels+1), INTENT(IN) ::                           &
   ft_nt,                                                               &
                                ! IN Non-turbulent heat (rho*Km/s) and
   fq_nt                    !      moisture (rho*m/s) fluxes

! INOUT fields
  LOGICAL, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end), INTENT(INOUT) ::       &
   coupled,                                                             &
                                ! INOUT Flag to indicate Sc layer
                                !       weakly decoupled
   cumulus,                                                             &
                                ! INOUT Flag for cumulus
   dsc                      ! INOUT Flag set if decoupled stratocu
                                !       layer found.

  INTEGER, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end), INTENT(INOUT) ::       &
   ntml,                                                                &
                                ! INOUT  Number of turbulently mixed
                                !        layers.
   ntdsc                    ! INOUT  Top level of any decoupled
                                !        turbulently mixed Sc layer

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
        INTENT(INOUT) ::                                                &
   zc,                                                                  &
                                ! INOUT Cloud depth (not cloud fraction
                                !    weighted) (m).
   zc_dsc,                                                              &
                                ! INOUT Cloud depth (not cloud fraction
                                !    weighted) (m).
   zhsc,                                                                &
                                ! INOUT Cloud-layer height (m)
   zh,                                                                  &
                                ! INOUT Boundary layer height (m)
   dscdepth,                                                            &
                                ! INOUT Decoupled cloud-layer depth (m)
   ft_nt_zh,                                                            &
                                ! INOUT Non-turbulent heat (rho*Km/s)
   ft_nt_zhsc,                                                          &
                                !        and moisture (rho*m/s) fluxes
   fq_nt_zh,                                                            &
                                !        evaluated at the SML and DSC
   fq_nt_zhsc               !        inversions

! OUT fields
  INTEGER, INTENT(OUT) ::                                               &
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
         ! OUT Bottom level of any decoupled turbulently mixed Sc layer.

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  2:bl_levels), INTENT(OUT) ::                          &
   rhokm,                                                               &
                                ! OUT Layer k-1 - to - layer k
                                !     turbulent mixing coefficient
                                !     for momentum (kg/m/s).
   rhokh,                                                               &
                                ! OUT Layer k-1 - to - layer k
                                !     turbulent mixing coefficient
                                !     for heat and moisture (kg/m/s).
   rhokm_top,                                                           &
                                ! OUT exchange coefficient for
                                !     momentum due to top-down mixing
   rhokh_top,                                                           &
                                ! OUT exchange coefficient for
                                !     heat and moisture due to top-down
                                !     mixing
   f_ngstress,                                                          &
                                ! OUT dimensionless function for
                                !     non-gradient stresses
   rhof2,                                                               &
                                ! OUT f2 and fsc term shape profiles
   rhofsc                   !       multiplied by rho

  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
        INTENT(OUT) ::                                                  &
   zdsc_base,                                                           &
                                ! OUT Height of base of K_top in DSC
   rhokh_surf_ent,                                                      &
                                ! OUT SML surface-driven entrainment KH
   rhokh_top_ent,                                                       &
                                ! OUT SML top-driven entrainment KH
   rhokh_dsct_ent           ! OUT DSC top-driven entrainment KH
!  ---------------------------------------------------------------------
!    Local and other symbolic constants :-





  REAL :: a_ent_1,a_ent_2,c_t,a_ent_shr,dec_thres_clear
  REAL :: g1
  INTEGER :: n_steps
  PARAMETER (                                                           &
   a_ent_2=0.056,                                                       &
                                ! Entrainment parameter.
   c_t=1.0,                                                             &
                                ! Parameter in Zilitinkevich term.
   dec_thres_clear=1.0,                                                 &
                                ! Decoupling threshold for cloud-free
                                ! boundary layers (larger makes
                                ! decoupling less likely)
   n_steps=3                                                            &
                                ! Number of steps through the mixed
                                ! layer per sweep
  )

  REAL :: s_m,a_ngs
  PARAMETER (                                                           &
   s_m   = 1.0,                                                         &
                                ! empirical parameters in
   a_ngs = 2.7                                                          &
                                ! non-gradient stresses
  )
!*
!  Define local storage.
!  (a) Workspace.

  INTEGER ::                                                            &
   ksurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
             ! First Theta-level above surface layer well-mixed SC layer
  LOGICAL, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) ::                      &
   scbase,                                                              &
                                ! Flag to signal base of CML reached
   test_well_mixed,                                                     &
                                ! Flag to test wb integration
                                ! for a well-mixed layer
   ksurf_iterate,                                                       &
                                ! Flag to perform iteration to
                                ! find top of Ksurf
   ktop_iterate             ! Flag to perform iteration to
                                ! find base of Ktop
  REAL ::                                                               &
   kh_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           2:bl_levels)
                                ! Shape factor for non-local
                                ! turbulent mixing coefficient
  REAL, DIMENSION(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)::&
   w_m_top,                                                             &
                      ! Turbulent velocity scale for momentum
                      !   evaluated at the top of the b.l.
   w_h_top,                                                             &
                      ! Turbulent velocity scale for scalars
                      !   evaluated at the top of the b.l.
   prandtl_top,                                                         &
                      ! Turbulent Prandtl number
                      !   evaluated at the top of the b.l.
   kh_top_factor,                                                       &
                      ! Factors to ensure K_H and K_M profiles are
   km_top_factor,                                                       &
                      !   continuous at ZH and ZHSC
   kh_sct_factor,                                                       &
                      !               "
   km_sct_factor,                                                       &
                      !               "
   kh_dsct_factor,                                                      &
                      !               "
   km_dsct_factor,                                                      &
                      !               "
   v_top,                                                               &
                      ! Velocity scale for top-down convection
   v_top_dsc,                                                           &
                      !
   v_sum,                                                               &
                      ! total velocity scale
   v_sum_dsc,                                                           &
                      !
   zsml_top,                                                            &
                      ! Height of top of surf-driven K in SML
   zsml_base,                                                           &
                      ! Height of base of top-driven K in SML
   Kprof_cu_top,                                                        &
                      ! Height of top of K profile in cumulus
   rhokh_lcl,                                                           &
                      ! rho*KH at the LCL in cumulus
   cu_depth_scale,                                                      &
                      ! Depth scale for decay of K profile 
                      ! above LCL in cumulus (in m)
   v_surf,                                                              &
                      ! Velocity scale for surface-up conve
   scdepth,                                                             &
                      ! Depth of top-driven mixing in SML
   z_inv,                                                               &
                      ! inversion height (top of K profiles)
   wb_surf_int,                                                         &
                      ! Estimate of wb integrated over surface layer
   wb_dzrad_int,                                                        &
                      ! Estimate of wb integrated over cloud-top region
   dzrad,                                                               &
                      ! Depth of cloud-top (radiatively cooled) region
   v_ktop,                                                              &
                      ! velocity scale for K_top profile
   v_ksum,                                                              &
                      ! total velocity scale
   z_cbase,                                                             &
                      ! cloud base height
   wb_ratio,                                                            &
                      ! WBN_INT/WBP_INT
   dec_thres,                                                           &
                      ! Local decoupling threshold
   wbp_int,                                                             &
                      ! Positive part of buoyancy flux integral
   wbn_int,                                                             &
                      ! Negative part of buoyancy flux integral
   zinv_pr,                                                             &
                      ! Height of layer top above surface
   khtop,                                                               &
                      ! temporary KH_top in wb integration
   khsurf,                                                              &
                      ! temporary KH_surf in wb integration
   zwb0,                                                                &
                      ! height at which wb assumed to go to zero
   z_top_lim,                                                           &
                      ! upper height limit on K profile
   z_bot_lim,                                                           &
                      ! lower height limit on K profile
   z_inc,                                                               &
                      ! Step size (m)
   rho_we,                                                              &
                      ! rho*param.d entrainment rates...
   rho_we_sml,                                                          &
                      !  ...for surf and DSC layers (kg/m2/s)
   rho_we_dsc,                                                          &
                      !
   cf_ml,                                                               &
                      ! Mixed layer cloud fraction
   df_ctop,                                                             &
                      ! Cloud-top radiative flux divergence
   dqw,                                                                 &
                      ! QW jump across inversion
   dsl,                                                                 &
                      ! SL jump across inversion
   ft_nt_dscb,                                                          &
                      ! Non-turbulent heat flux at DSC base
   fq_nt_dscb,                                                          &
                      ! Non-turbulent moisture flux at DSC base
   tothf_zi,                                                            &
                      ! Total heat and moisture fluxes at
   totqf_zi       !   the inversion height, Zi

!  (b) Scalars.

  REAL ::                                                               &
   Prandtl,                                                             &
                    ! Turbulent Prandtl number.
   pr_neut,                                                             &
                    ! Neutral limit for Prandtl number
   pr_conv,                                                             &
                    ! Convective limit for Prandtl number
   zk_uv,                                                               &
                    ! Height above surface of u,v-level.
   zk_tq,                                                               &
                    ! Height above surface of T,q-level.
   wstar3,                                                              &
                    ! Cube of free-convective velocity scale
   c_ws,                                                                &
                    ! Empirical constant multiplying Wstar
   w_s_cubed_uv,                                                        &
                    ! WSTAR for u,v-level
   w_s_cubed_tq,                                                        &
                    !   and T,q-level
   w_m_uv,                                                              &
                    ! Turbulent velocity scale for momentum: u,v-level
   w_m_tq,                                                              &
                    !   and T,q-level
   w_h_uv,                                                              &
                    ! Turbulent velocity scale for scalars: u,v-level
   w_h_tq,                                                              &
                    !   and T,q-level
   w_m_hb_3,                                                            &
                    ! Cube of W_M, as given by Holtslag and Boville, 93
   w_m_neut,                                                            &
                    ! Neutral limit for W_M
   sf_term,                                                             &
                    ! Surface flux term for entrainment parametrization.
   sf_shear_term,                                                       &
                    ! Surface shear term for entrainment paramn.
   ir_term,                                                             &
                    ! Indirect radiative term for entrainment paramn.
   dr_term,                                                             &
                    ! Direct radiative term for entrainment paramn.
   evap_term,                                                           &
                    ! Evaporative term in entrainment parametrization.
   zil_corr,                                                            &
                    ! Zilitinkevich correction term in entrn. paramn.
   zeta_s_fac,                                                          &
                    ! Factor involving ZETA_S.
   zeta_r_sq,                                                           &
                    ! ZETA_R squared.
   zr,                                                                  &
                    ! Ratio ZC/ZH.
   z_pr,                                                                &
                    ! Height above surface layer
   zh_pr,                                                               &
                    ! Height of layer top above surface
   z_ratio,                                                             &
                    ! Ratio of heights
   zcml_base,                                                           &
                    ! Height of base of cloud mixed layer
   rhokh_ent,                                                           &
                    ! entrainment eddy viscosity
   frac_top,                                                            &
                    ! Fraction of turbulent mixing driven from the top
   factor,                                                              &
                    ! Temporary scalar
   alpha_t,                                                             &
                    ! Parametrized fraction of cloud-top
                    ! radiative cooling within the inversion
   dz_inv,                                                              &
                    ! Parametrizzed inversion thickness (m)
   l_rad,                                                               &
                    ! Estimate of e-folding radiative flux
                    ! decay depth (assumed >= 25m)
   wb_cld,                                                              &
                     ! Cloud layer buoyancy flux
   wb_scld,                                                             &
                     ! Sub-cloud layer buoyancy flux
   cld_frac,                                                            &
                     ! Vertical fraction of layer containing cloud
   zb_ktop,                                                             &
                     ! height of base of K_top profile
   db_ratio,                                                            &
                     ! Temporary in ZWB0 calculation
   gamma_wbs,                                                           &
                     ! Surface layer wb gradient
   wsl_dzrad_int,                                                       &
                     ! Estimate of wsl and wqw integrated over
   wqw_dzrad_int,                                                       &
                     !   the cloud-top region
   wslng,                                                               &
                     ! Non-gradient part of SL flux
   wqwng,                                                               &
                     ! Non-gradient part of QW flux
   f2, fsc       ! Shape functions for non-gradient fluxes

  INTEGER ::                                                            &
   i,j,                                                                 &
                    ! Loop counter (horizontal field index).
   k,                                                                   &
                    ! Loop counter (vertical level index).
   n_sweep,                                                             &
                    ! sweep counter
   ns           ! step counter

! 2D arrays for optimisation

  INTEGER, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) ::                      &
   ntop,                                                                &
                     ! top level of surf-driven K profile
   ntml_new,                                                            &
                     ! temporary in NTML calculation
   kwb0          ! level at which wb assumed to go to zero

  INTEGER, DIMENSION(pdims%i_end*pdims%j_end) ::                        &
   up            ! indicator of upward/downward sweep

  LOGICAL, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) :: status_ntml

! Array introduced to calculate kwb0
  LOGICAL, DIMENSION(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end) :: kstatus

! Variables for vector compression

  INTEGER :: ij_len
  INTEGER :: ic
  INTEGER :: c_len
  LOGICAL,DIMENSION(pdims%i_end*pdims%j_end)  :: to_do
  INTEGER,DIMENSION(pdims%i_end*pdims%j_end)  :: ind_todo
  INTEGER :: c_len_i
  LOGICAL,DIMENSION(pdims%i_end*pdims%j_end)  :: todo_inner
  INTEGER,DIMENSION(pdims%i_end*pdims%j_end)  :: ind_todo_i

  INTEGER :: i1, j1, l

  INTEGER, PARAMETER :: jblock = 4  !Cache blocking - block size
  INTEGER            :: jj          !Cache blocking - loop index 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('EXCF_NL',zhook_in,zhook_handle)
!
!----------------------------------------------------------------------- 
!     Set values of A_ENT_1, G1 and A_ENT_SHR 

  IF (L_RP2) THEN 
    a_ent_1 = a_ent_1_rp       ! Entrainment parameter 
    g1 = g1_rp                 ! Velocity scale parameter 
  ELSE 
    a_ent_1 = 0.23             ! Entrainment parameter 
    g1 = 0.85                  ! Velocity scale parameter 
  END IF 

  a_ent_shr = 5.0 * a_ent_1 / 0.23   ! Entrainment parameter. 

!-----------------------------------------------------------------------
! Index to subroutine EXCFNL8C

! 0. Calculate top-of-b.l. velocity scales and Prandtl number.
! 1. Calculate the top-of-b.l. entrainment parametrization
! 2. Estimate the depths of top-down and surface-up mixing.
!   2.1 First test for well-mixed boundary layer
!   2.2 Iterate to find top of surface-driven mixing, ZSML_TOP,
!   2.3 Iterate to find the base of the top-driven K profile, ZDSC_BASE.
! 3. Calculate factors required to ensure that the K profiles are
!    continuous at the inversion
! 4. Calculate height dependent turbulent transport coefficients
!    within the mixing layers.

!-----------------------------------------------------------------------
! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.
!-----------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(SHARED)                                       &
!$OMP& PRIVATE(i, j, k, jj, c_ws, wstar3, pr_neut, pr_conv, w_m_neut, &
!$OMP& zeta_s_fac, sf_term, sf_shear_term, zeta_r_sq, ir_term, zr,    &
!$OMP& evap_term, dz_inv, l_rad, alpha_t, dr_term, zil_corr,          &
!$OMP& rhokh_ent, frac_top, zh_pr, factor,                            &
!$OMP& wsl_dzrad_int, wqw_dzrad_int, db_ratio, zb_ktop, f2, fsc,      &
!$OMP& z_ratio, z_pr, wslng, wqwng, wb_scld, wb_cld, cld_frac, l,     &
!$OMP& j1, i1, ic,  w_m_hb_3, zk_uv, zk_tq, Prandtl, w_h_uv, w_h_tq,  &
!$OMP& w_m_uv, w_m_tq, w_s_cubed_tq, w_s_cubed_uv, gamma_wbs)

!cdir collapse
!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      rhokh_surf_ent(i,j) = 0.0
      rhokh_top_ent(i,j) = 0.0
      rhokh_dsct_ent(i,j) = 0.0
      v_top(i,j) = 0.0
      v_surf(i,j)= 0.0
      v_sum(i,j) = 0.0
      v_top_dsc(i,j) = 0.0
      v_sum_dsc(i,j) = 0.0
      rho_we(i,j)     = 0.0
      rho_we_sml(i,j) = 0.0
      rho_we_dsc(i,j) = 0.0

      IF (fb_surf(i,j)  >=  0.0) THEN

          ! Free-convective velocity scale cubed

        IF (coupled(i,j)) THEN
          wstar3 = zhsc(i,j) * fb_surf(i,j)
        ELSE
          wstar3 =   zh(i,j) * fb_surf(i,j)
        END IF

        IF (flux_grad  ==  Locketal2000) THEN

            ! Turbulent velocity scale for momentum

          c_ws = 0.25
          w_m_top(i,j) = (v_s(i,j)*v_s(i,j)*v_s(i,j) +                  &
                          c_ws*wstar3)**(1.0/3.0)

            ! Turbulent Prandtl number and velocity scale for scalars
            ! gives 0.375<Pr<0.75 for convective to neutral conditions
          prandtl_top(i,j) = 0.75 *                                     &
                       ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +          &
                           (1.0/25.0)*wstar3*w_m_top(i,j) ) /           &
                       ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +          &
                           (2.0/25.0)*wstar3*w_m_top(i,j) )
          w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
        ELSE IF (flux_grad  ==  HoltBov1993) THEN
          c_ws = 0.6
          w_m_top(i,j) = (v_s(i,j)*v_s(i,j)*v_s(i,j) +                  &
                          c_ws*wstar3)**(1.0/3.0)
            ! Using Lock et al interpolation but with
            ! HB93 range of 0.6<Pr<1
          pr_neut = 1.0
          pr_conv = 0.6
          prandtl_top(i,j) = pr_neut *                                  &
                ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                 &
                                  (1.0/25.0)*wstar3*w_m_top(i,j) ) /    &
                ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                 &
                    (pr_neut/(25.0*pr_conv))*wstar3*w_m_top(i,j) )
          w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
        ELSE IF (flux_grad  ==  LockWhelan2006) THEN
          pr_neut = 0.75
          pr_conv = 0.6
          c_ws    = 0.42   !  ~ 0.75^3
            ! Slightly contrived notation since really we know W_H_TOP
            ! and Prandtl range but this makes similarity to
            ! Lock et al clearer (possibly!)
          w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) +                     &
                       c_ws*wstar3 )**(1.0/3.0)
          w_h_top(i,j) = w_m_neut / pr_neut

          prandtl_top(i,j) = pr_neut *                                  &
                ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                 &
                    (1.0/25.0)*wstar3*w_m_neut ) /                      &
                ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                 &
                    (pr_neut/(25.0*pr_conv))*wstar3*w_m_neut )

          w_m_top(i,j) = w_h_top(i,j) * prandtl_top(i,j)
        END IF

      ELSE
        w_m_top(i,j) = v_s(i,j)
        prandtl_top(i,j) = 0.75
        w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 1. Calculate the top-of-b.l. entrainment parametrization
!-----------------------------------------------------------------------
! 1.1 Initialise 3D arrays
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
!cdir collapse
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        rhokh(i,j,k) = 0.0
        rhokm(i,j,k) = 0.0
        rhokh_top(i,j,k) = 0.0
        rhokm_top(i,j,k) = 0.0
        rhof2(i,j,k)  = 0.0
        rhofsc(i,j,k) = 0.0
        f_ngstress(i,j,k) = 0.0
        kh_surf(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO

!cdir collapse
!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
!-----------------------------------------------------------------------
! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
!     and store b.l. top quantities for later use.
!-----------------------------------------------------------------------
!      FIRST the top of the SML (if not coupled)
!-----------------------------------------------
      k = ntml(i,j)+1
      IF ( .NOT.coupled(i,j) .AND. fb_surf(i,j)  >=  0.0 ) THEN
                ! Correction (controlled by STOPWE_SBL in 6B scheme):
                !   require FB_SURF>0 too
            !-----------------------------------------------------------
            ! Calculate the surface buoyancy flux term
            !-----------------------------------------------------------
        zeta_s_fac = (1.0 - zeta_s(i,j)) * (1.0 - zeta_s(i,j))
        sf_term = a_ent_1 * MAX ( 0.0 ,                                 &
                            ( (1.0 - zeta_s_fac) * bflux_surf(i,j)      &
                              + zeta_s_fac * bflux_surf_sat(i,j) ) )
            !-----------------------------------------------------------
            ! Calculate the surface shear term
            !-----------------------------------------------------------
        sf_shear_term =  a_ent_shr * v_s(i,j) * v_s(i,j) * v_s(i,j)     &
                        * rho_uv(i,j,k)  / zh(i,j)
            !-----------------------------------------------------------
            ! Calculate the indirect radiative term
            !-----------------------------------------------------------
        zeta_r_sq = zeta_r(i,j)*zeta_r(i,j)
        ir_term = ( bt_top(i,j)*zeta_r_sq +                             &
                    btt_top(i,j)*(1.0-zeta_r_sq) )                      &
                  * a_ent_1 * df_top_over_cp(i,j)
            !-----------------------------------------------------------
            ! Calculate the evaporative term
            !-----------------------------------------------------------
        IF ( db_top(i,j)  >   0.0) THEN
          zr = SQRT( zc(i,j) / zh(i,j) )
          evap_term = a_ent_2 * rho_uv(i,j,k)                           &
                    * chi_s_top(i,j) * chi_s_top(i,j)                   &
                    * zr * zr * zr * db_top_cld(i,j)                    &
                    * SQRT( zh(i,j) * db_top(i,j) )
        ELSE
          evap_term = 0.0
        END IF
            !-----------------------------------------------------------
            ! Combine forcing terms to calculate the representative
            ! velocity scales
            !-----------------------------------------------------------
        v_sum(i,j) = ( (sf_term + sf_shear_term +                       &
                        ir_term + evap_term)                            &
                     * zh(i,j) /(a_ent_1*rho_uv(i,j,k)) )**(1.0/3.0)
        v_top(i,j) = ( (ir_term+evap_term) * zh(i,j)                    &
                             / (a_ent_1*rho_uv(i,j,k)) )**(1.0/3.0)
        v_surf(i,j) = ( (sf_term) * zh(i,j)                             &
                             / (a_ent_1*rho_uv(i,j,k)) )**(1.0/3.0)
            !-----------------------------------------------------------
            ! Calculate the direct radiative term
            !  can only calculate for DB_TOP > 0
            !-----------------------------------------------------------
        IF ( db_top(i,j)  >   0.0) THEN
          dz_inv  = MIN( v_sum(i,j)*v_sum(i,j) / db_top(i,j) ,100.0 )
          l_rad   = 15.0 * MAX( 1.0 , 200./(zc(i,j)+1.0e-14) )
          alpha_t = 1.0 - EXP(-0.5*dz_inv/l_rad)
             ! Make enhancement due to buoyancy reversal feedback
          alpha_t = alpha_t + br_fback(i,j)*(1.0-alpha_t)
          dr_term = btc_top(i,j) * alpha_t * df_top_over_cp(i,j)
             !----------------------------------------------------------
             ! Combine terms to calculate the entrainment
             ! mixing coefficients
             !----------------------------------------------------------
          zil_corr = c_t * ( (sf_term + sf_shear_term +                 &
                              ir_term + evap_term) /                    &
                      (rho_uv(i,j,k) * SQRT(zh(i,j))) )**(2.0/3.0)

          rho_we_sml(i,j) = (sf_term + sf_shear_term                    &
                      +      ir_term + evap_term + dr_term)             &
                         / ( db_top(i,j) + zil_corr )

          rhokh_ent = rho_we_sml(i,j)/ rdz(i,j,k)

          frac_top = v_top(i,j) / ( v_top(i,j)+w_h_top(i,j)+1.0e-14 )
          rhokh_surf_ent(i,j) = rhokh_ent * ( 1.0 - frac_top )
          rhokh_top_ent(i,j) = rhokh_ent * frac_top

          rhokm(i,j,k) = prandtl_top(i,j) * rhokh_surf_ent(i,j)         &
                       * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))       &
                       * rho_tq(i,j,k-1) / rho_uv(i,j,k)
          rhokm_top(i,j,k) = 0.75 * rhokh_top_ent(i,j)                  &
                       * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))       &
                       * rho_tq(i,j,k-1) / rho_uv(i,j,k)

        END IF    ! test on DB_TOP GT 0
      END IF
!----------------------------------------------------------------
!      THEN the top of the DSC (if coupled use ZHSC length-scale)
!----------------------------------------------------------------
      IF ( ntdsc(i,j)  >   0 ) THEN
        k = ntdsc(i,j)+1
        IF (coupled(i,j)) THEN
              !--------------------------------------------------------
              ! Calculate the surface buoyancy flux term
              !--------------------------------------------------------
          zeta_s_fac = (1.0 - zeta_s(i,j)) * (1.0 - zeta_s(i,j))
          sf_term = a_ent_1 * MAX ( 0.0 ,                               &
                            ( (1.0 - zeta_s_fac) * bflux_surf(i,j)      &
                              + zeta_s_fac * bflux_surf_sat(i,j) ) )
              !--------------------------------------------------------
              ! Calculate the surface shear term
              !--------------------------------------------------------
          IF (fb_surf(i,j)  >=  0.0) THEN
            sf_shear_term = a_ent_shr * v_s(i,j)*v_s(i,j)*v_s(i,j)      &
                            * rho_uv(i,j,k)  / zhsc(i,j)
          ELSE
            sf_shear_term = 0.0
          END IF
          v_surf(i,j) = ( (sf_term) * zhsc(i,j)                         &
                              / (a_ent_1*rho_uv(i,j,k)) )**(1.0/3.0)
          IF (entr_smooth_dec == on) THEN
            ! taper surface terms to zero depending on svl_diff_frac
            sf_term =       sf_term       * svl_diff_frac(i,j)
            sf_shear_term = sf_shear_term * svl_diff_frac(i,j)
          END IF
        ELSE
          sf_term = 0.0
          sf_shear_term = 0.0
        END IF
          !-----------------------------------------------------------
          ! Calculate the indirect radiative term
          !-----------------------------------------------------------
        zeta_r_sq = zeta_r_dsc(i,j)*zeta_r_dsc(i,j)
        ir_term = ( bt_dsct(i,j)*zeta_r_sq +                            &
                      btt_dsct(i,j)*(1.0-zeta_r_sq) )                   &
                      * a_ent_1 * df_dsct_over_cp(i,j)
          !-----------------------------------------------------------
          ! Calculate the evaporative term
          !-----------------------------------------------------------
        IF (db_dsct(i,j)  >   0.0) THEN
          zr = SQRT( zc_dsc(i,j) / dscdepth(i,j) )
          evap_term = a_ent_2 * rho_uv(i,j,k)                           &
                    * chi_s_dsct(i,j) * chi_s_dsct(i,j)                 &
                    * zr * zr * zr * db_dsct_cld(i,j)                   &
                    * SQRT( dscdepth(i,j) * db_dsct(i,j) )
        ELSE
          evap_term = 0.0
        END IF
          !-----------------------------------------------------------
          ! Combine forcing terms to calculate the representative
          ! velocity scales
          !-----------------------------------------------------------
        v_sum_dsc(i,j) = ( (sf_term + sf_shear_term +                   &
                              ir_term + evap_term)                      &
                  * dscdepth(i,j) / (a_ent_1*rho_uv(i,j,k)) )**(1./3.)
        v_top_dsc(i,j) =( (ir_term + evap_term) * dscdepth(i,j) /       &
                             (a_ent_1*rho_uv(i,j,k)) )**(1.0/3.0)
          !-----------------------------------------------------------
          ! Calculate the direct radiative term
          !-----------------------------------------------------------
        IF (db_dsct(i,j)  >   0.0) THEN
          dz_inv  = MIN( v_sum_dsc(i,j)*v_sum_dsc(i,j)/db_dsct(i,j),    &
                         100.0 )
          l_rad   = 15.0 * MAX( 1.0 , 200./(zc_dsc(i,j)+1.0) )
          alpha_t = 1.0 - EXP(-0.5*dz_inv/l_rad)
             ! Make enhancement due to buoyancy reversal feedback
          alpha_t = alpha_t + br_fback_dsc(i,j)*(1.0-alpha_t)
          dr_term = btc_top(i,j) * alpha_t * df_dsct_over_cp(i,j)
             !----------------------------------------------------------
             ! Finally combine terms to calculate the entrainment
             ! rate and mixing coefficients
             !----------------------------------------------------------
          zil_corr = c_t * ( (sf_term + sf_shear_term +                 &
                              ir_term + evap_term) /                    &
                   (rho_uv(i,j,k) * SQRT(dscdepth(i,j))) )**(2.0/3.0)
          rho_we_dsc(i,j) = ( sf_term + sf_shear_term                   &
                  +           ir_term + evap_term + dr_term )           &
                             / ( db_dsct(i,j) + zil_corr )
          rhokh_dsct_ent(i,j) = rho_we_dsc(i,j)/ rdz(i,j,k)
          rhokm_top(i,j,k) = 0.75 * rhokh_dsct_ent(i,j)                 &
                       * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))       &
                       * rho_tq(i,j,k-1) / rho_uv(i,j,k)
        END IF   ! test on DB_DSCT gt 0
      END IF
    END DO
  END DO
!$OMP END DO

!  If there is no turbulence generation in DSC layer, ignore it.

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( v_top_dsc(i,j)  <=  0.0 ) THEN
        dsc(i,j) = .FALSE.
        ntdsc(i,j) = 0
        zhsc(i,j) = 0.0
        zc_dsc(i,j) = 0.0
        dscdepth(i,j) = 0.0
        coupled(i,j) = .FALSE.
      END IF
    END DO
  END DO
!$OMP END DO

! ----------------------------------------------------------------------
! 2.0 Estimate the depths of top-down and surface-up mixing.
!     These amount to diagnoses of recoupling and decoupling.
!     The K_H profiles are applied over layers such that the ratio
!        WBN_INT/WBP_INT = DEC_THRES (parameter),
!     where WBN_INT and WBP_INT are the magnitudes of the integrals of
!     the negative and positive parts, respectively, of the resulting
!     buoyancy flux profile (given by - KH * DB_FOR_FLUX).
! ----------------------------------------------------------------------
! 2.1 First test for well-mixed boundary layer
!     (ie. both KH profiles extending from cloud-top to the surface).
!     If the parcel ascent diagnosed:
!        DSC    - test for well-mixed up to ZHSC = recoupling
!        no DSC - test for well-mixed up to ZH   = decoupling
! -----------------------------------------------------------
! Default settings

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      zsml_top(i,j)  = zh(i,j)
      zsml_base(i,j) = 0.1 * zh(i,j)
      zdsc_base(i,j) = 0.1 * zhsc(i,j)
      dec_thres(i,j) = dec_thres_cloud  ! use cloudy by default
      z_inv(i,j) = 0.0    ! inversion height (top of K profiles)
      zwb0(i,j)  = 0.0    ! height at which WB goes to zero
      wbp_int(i,j) = 0.0
      wbn_int(i,j) = 0.0
      wb_surf_int(i,j) = 0.0
      wb_dzrad_int(i,j) = 0.0
      dzrad(i,j) = 100.0
      tothf_zi(i,j) = 0.0
      totqf_zi(i,j) = 0.0
      kstatus(i,j)= .TRUE.
      kwb0(i,j)  = 2
      ntop(i,j)  = -1
      ksurf(i,j) = 1
    END DO
  END DO
!$OMP END DO

! Find KSURF, the first theta-level above the surface layer

!$OMP DO SCHEDULE(STATIC)
  DO jj = pdims%j_start, pdims%j_end, jblock
    DO k = 2, bl_levels
!cdir collapse
      DO j = jj, MIN(jj+jblock-1, pdims%j_end)
        DO i = pdims%i_start,pdims%i_end
          IF ( z_tq(i,j,k-1)  <   0.1*zh(i,j) ) ksurf(i,j) = k
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

! Find kwb0, level with lowest positive cloud-free buoyancy gradient

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, jblock
  DO k = 2, bl_levels
    DO j = jj, MIN(jj+jblock-1, pdims%j_end)
      DO i = pdims%i_start,pdims%i_end
        IF (kstatus(i,j)) THEN
          IF ( (db_ga_dry(i,j,k) <=  0.0) .OR.                          &
               (k >= ntml(i,j)) ) THEN
            kstatus(i,j)=.FALSE.
            kwb0(i,j)=k
          END IF
        END IF
      END DO
    END DO
  END DO
END DO
!$OMP END DO

! Set flags for iterating wb integral to calculate depth of mixing,
! one each for KSURF and K_TOP.  Note these will be updated depending on
! what happpens on testing for a well-mixed layer in section 2.2.

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      test_well_mixed(i,j) = .FALSE.
      ksurf_iterate(i,j)= .FALSE.
      ktop_iterate(i,j) = .FALSE.
      IF ( ntdsc(i,j)  >   2 ) THEN
        ktop_iterate(i,j) = .TRUE.
      END IF
    END DO ! I
  END DO ! J
!$OMP END DO

!cdir collapse
!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( bflux_surf(i,j)  >   0.0) THEN
          ! can only be coupled to an unstable SML
        IF ( ntdsc(i,j)  >   2 ) THEN
            ! well-resolved DSC layer (6B DECFIX included)
            ! ...test for recoupling with SML
          test_well_mixed(i,j) = .TRUE.
          z_inv(i,j)  = zhsc(i,j)
          z_cbase(i,j)= z_inv(i,j) - zc_dsc(i,j)
          cf_ml(i,j)  = cf_dsc(i,j)
          v_ktop(i,j) = v_top_dsc(i,j)
          v_ksum(i,j) = v_sum_dsc(i,j)
          df_ctop(i,j)= df_dsct_over_cp(i,j)
          rho_we(i,j) = rho_we_dsc(i,j)
          dsl(i,j)    = dsl_dsc(i,j)
          dqw(i,j)    = dqw_dsc(i,j)
          tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zhsc(i,j)
          totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zhsc(i,j)
          ntop(i,j)   = ntdsc(i,j) - 1
            ! assuming wb goes to zero by the lowest of ZH
            ! or cloud-base, but above surface layer
        ELSE IF ( .NOT.dsc(i,j) .AND. .NOT.cumulus(i,j) .AND.           &
                  ntml(i,j)  >   2) THEN
            ! well-resolved SML       (6B DECFIX included)
            ! ...test for decoupling
            ! Note: code can only deal with one DSC layer at a time so
            ! can't decouple SML if a DSC layer already exists.
            ! ---------------------------------------------------------
            ! 6B DECFIX switch included:
            ! If the BL layer is cloud-free then use a less restrictive
            ! threshold - ideally, the parcel ascent would have
            ! found the correct BL top in this case but this test is
            ! kept to keep negative buoyancy fluxes under control
            ! (ie. DEC_THRES_CLEAR=1 ensures wbn_int < |wbp_int|)
          IF (zc(i,j)  ==  0.0) dec_thres(i,j) = dec_thres_clear
            ! ---------------------------------------------------------
          test_well_mixed(i,j) = .TRUE.
          z_inv(i,j)  = zh(i,j)
          z_cbase(i,j)= z_inv(i,j) - zc(i,j)
          cf_ml(i,j)  = cf_sml(i,j)
          v_ktop(i,j) = v_top(i,j)
          v_ksum(i,j) = v_sum(i,j)
          df_ctop(i,j)= df_top_over_cp(i,j)
          rho_we(i,j) = rho_we_sml(i,j)
          dsl(i,j)    = dsl_sml(i,j)
          dqw(i,j)    = dqw_sml(i,j)
          tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zh(i,j)
          totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zh(i,j)
          ntop(i,j)   = ntml(i,j) - 1
        END IF
      END IF
    END DO ! I
  END DO ! J
!$OMP END DO

! ----------------------------------------------------------------------
! 2.1.1 Estimate wb integral over radiatively cooled cloud-top region,
!       from Z_INV to Z_INV-DZRAD.
!       DZRAD taken to be constant (100m) for simplicity, but also
!       integration depth taken to extend down to at least
!       Z_TQ(NTML-1) but without going below cloud-base.
!       This code also invoked for cloud-free cases
!           - does this matter?
!           - only ignoring top grid level of integration which
!             shouldn't be important/relevant for cloud-free layers.
! ----------------------------------------------------------------------
!cdir collapse
!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( test_well_mixed(i,j) ) THEN
        dzrad(i,j) = 100.0
        IF ( ntop(i,j)  >   2 ) THEN
          DO WHILE ( z_tq(i,j,ntop(i,j))  >   z_inv(i,j)-dzrad(i,j)     &
               .AND. z_tq(i,j,ntop(i,j)-1)  >   z_inv(i,j)-z_cbase(i,j))
            ntop(i,j) = ntop(i,j) - 1
            IF (ntop(i,j) == 2 ) THEN
              EXIT
            END IF
          END DO
        END IF
        dzrad(i,j) = z_inv(i,j) - z_tq(i,j,ntop(i,j))

        wsl_dzrad_int = dzrad(i,j) *                                    &
                        ( 0.66*df_ctop(i,j) - rho_we(i,j)*dsl(i,j) )
        wqw_dzrad_int = - dzrad(i,j) * rho_we(i,j) * dqw(i,j)

        wb_dzrad_int(i,j) = wsl_dzrad_int * (                           &
                              (1.0-cf_ml(i,j))*btm(i,j,ntop(i,j)+1) +   &
                                cf_ml(i,j)*btm_cld(i,j,ntop(i,j)+1) )   &
                          + wqw_dzrad_int * (                           &
                              (1.0-cf_ml(i,j))*bqm(i,j,ntop(i,j)+1) +   &
                                cf_ml(i,j)*bqm_cld(i,j,ntop(i,j)+1) )
        wb_dzrad_int(i,j) = g * wb_dzrad_int(i,j)
        wb_dzrad_int(i,j) = MAX( 0.0, wb_dzrad_int(i,j) )

          ! Include WB_DZRAD_INT in WBP_INT as it set to be >0
        wbp_int(i,j) = wbp_int(i,j) + wb_dzrad_int(i,j)

      ELSE

        wb_dzrad_int(i,j) = -1.0  ! To identify not calculated

      END IF
    END DO ! I
  END DO ! J
!$OMP END DO

! For WB diagnostics, convert integrated WB to uniform profile
! ----------------------------------------------------------------------
! 2.1.2 Estimate wb integral over surface layer
!       (and up to next theta-level, namely Z_TQ(KSURF) )
!       assuming a linear profile going to zero at ZWB0
! ----------------------------------------------------------------------
!cdir collapse

!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( test_well_mixed(i,j) ) THEN
        IF ( kwb0(i,j)  ==  ntml(i,j) ) THEN
          zwb0(i,j) = zh(i,j)
        ELSE IF ( kwb0(i,j)  ==  2 ) THEN
          zwb0(i,j) = z_uv(i,j,2)
        ELSE
          k=kwb0(i,j)
            ! now DB_GA_DRY(K) LE 0 and DB_GA_DRY(K-1) GT 0
            ! so interpolate:
          db_ratio = db_ga_dry(i,j,k-1)                                 &
                  / ( db_ga_dry(i,j,k-1) - db_ga_dry(i,j,k) )
          db_ratio = MAX( 0.0, db_ratio )  ! trap for rounding error
          zwb0(i,j)=z_uv(i,j,k-1) +                                     &
                    db_ratio * (z_uv(i,j,k)-z_uv(i,j,k-1))
        END IF
        wb_surf_int(i,j) = bflux_surf(i,j) * z_tq(i,j,ksurf(i,j)) *     &
                      ( 1.0 - z_tq(i,j,ksurf(i,j))/(2.0*zwb0(i,j)))
        wb_surf_int(i,j) = MAX( 1.0e-14, wb_surf_int(i,j) )
      ELSE
          ! only include surface layer contribution for unstable mixing
        wb_surf_int(i,j) = 1.0e-14
      END IF

      wbp_int(i,j) = wbp_int(i,j) + wb_surf_int(i,j) ! must be >0

        ! Save surface and bl-top layer integral for diagnostics
    END DO ! I
  END DO ! J
!$OMP END DO

! ----------------------------------------------------------------------
! 2.1.3 Loop over well-mixed boundary layer integrating WB
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC)
DO jj = pdims%j_start, pdims%j_end, jblock
  DO k = 2, bl_levels-1
    DO j = jj, MIN(jj+jblock-1, pdims%j_end)
      DO i = pdims%i_start,pdims%i_end

        IF ( test_well_mixed(i,j) ) THEN
          ! ----------------------------------------------
          ! worth testing layer as well-mixed to cloud-top
          ! ----------------------------------------------
          zb_ktop = 0.1*z_inv(i,j)
          zinv_pr(i,j) = z_inv(i,j) - zb_ktop
          ! DB(K)is the K to K-1 difference and already
          ! integrated up to KSURF, so start this loop at KSURF+1
          IF ( (k >= ksurf(i,j)+1) .AND. (k <= ntop(i,j)) ) THEN

            khtop(i,j) = 0.0
            khsurf(i,j)= 0.0
            f2         = 0.0
            fsc        = 0.0

            z_pr = z_uv(i,j,k) - zb_ktop
            IF (z_pr  >   0.0 .AND. z_pr  <   zinv_pr(i,j)) THEN
              z_ratio = z_pr/zinv_pr(i,j)

              IF (flux_grad  ==  LockWhelan2006) THEN
                khtop(i,j) = 3.6 * vkman * rho_uv(i,j,k) * v_ktop(i,j)  &
                                 * zinv_pr(i,j) * (z_ratio**3)        &
                                 * ( (1.0-z_ratio)*(1.0-z_ratio) )
                f2 = rho_uv(i,j,k) * 0.5 * z_ratio                      &
                                         * 2.0**( z_ratio**4 )
                IF ( v_ksum(i,j)  >   0.0 ) THEN
                  fsc = rho_uv(i,j,k) * 3.5 * (v_ktop(i,j)/v_ksum(i,j)) &
                        * (z_ratio**3) * (1.0-z_ratio)
                END IF
              ELSE
                khtop(i,j) = g1 * vkman * rho_uv(i,j,k) * v_ktop(i,j)   &
                                * (( 1.0 - z_ratio )**0.8)              &
                                * z_pr * z_ratio
              END IF
            END IF

            z_pr = z_uv(i,j,k)
            IF ( z_pr  <   z_inv(i,j)) THEN
              !--------------------------------
              ! include surface-driven profile
              !--------------------------------
              khsurf(i,j) = vkman * rho_uv(i,j,k) *                     &
                       w_h_top(i,j)*z_pr*( 1.0 - z_pr/z_inv(i,j) )      &
                                        *( 1.0 - z_pr/z_inv(i,j) )
            END IF

            IF (flux_grad  ==  LockWhelan2006) THEN
              wslng = (f2+fsc)*tothf_zi(i,j) - ft_nt(i,j,k)
              wqwng = (f2+fsc)*totqf_zi(i,j) - fq_nt(i,j,k)
            END IF

            IF ( z_tq(i,j,k)  <=  z_cbase(i,j) ) THEN
              ! Completely below cloud-base so use cloud-free formula
              wb_scld = khsurf(i,j) * db_ga_dry(i,j,k) +                &
                        khtop(i,j) * db_noga_dry(i,j,k)
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_scld = wb_scld + ( g/rdz(i,j,k) ) *                  &
                  ( btm(i,j,k-1)*wslng + bqm(i,j,k-1)*wqwng )
              END IF
              wb_cld  = 0.0
            ELSE IF (z_tq(i,j,k-1)  >=  z_cbase(i,j)) THEN
              ! Completely above cloud-base so use cloudy formula
              wb_cld = ( khsurf(i,j) * db_ga_cld(i,j,k) +               &
                         khtop(i,j)  * db_noga_cld(i,j,k) )
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_cld = wb_cld + ( g/rdz(i,j,k) ) * (                  &
                         ( btm(i,j,k-1)*(1.0-cf_ml(i,j)) +              &
                           btm_cld(i,j,k-1)*cf_ml(i,j) )*wslng +        &
                         ( bqm(i,j,k-1)*(1.0-cf_ml(i,j)) +              &
                           bqm_cld(i,j,k-1)*cf_ml(i,j) )*wqwng )
              END IF
              wb_scld = 0.0
            ELSE
              ! cloud-base within this integration range
              ! so treat cloud and sub-cloud layer wb separately
              wb_scld = khsurf(i,j) * db_ga_dry(i,j,k) +                &
                        khtop(i,j) * db_noga_dry(i,j,k)
              wb_cld = ( khsurf(i,j) * db_ga_cld(i,j,k) +               &
                         khtop(i,j)  * db_noga_cld(i,j,k) )
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_scld = wb_scld + ( g/rdz(i,j,k) ) *                  &
                  ( btm(i,j,k-1)*wslng + bqm(i,j,k-1)*wqwng )
                wb_cld = wb_cld + ( g/rdz(i,j,k) ) * (                  &
                         ( btm(i,j,k-1)*(1.0-cf_ml(i,j)) +              &
                           btm_cld(i,j,k-1)*cf_ml(i,j) )*wslng +        &
                         ( bqm(i,j,k-1)*(1.0-cf_ml(i,j)) +              &
                           bqm_cld(i,j,k-1)*cf_ml(i,j) )*wqwng )
              END IF
              cld_frac = (z_tq(i,j,k)-z_cbase(i,j))                     &
                        /(z_tq(i,j,k)-z_tq(i,j,k-1))
              wb_cld  = cld_frac * wb_cld
              wb_scld = (1.0-cld_frac) * wb_scld
            END IF
            IF (wb_cld  >=  0.0) THEN
              wbp_int(i,j) = wbp_int(i,j) + wb_cld
            ELSE
              wbn_int(i,j) = wbn_int(i,j) - wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i,j) = wbp_int(i,j) + wb_scld
            ELSE
              wbn_int(i,j) = wbn_int(i,j) - wb_scld
            END IF

          END IF ! K

        END IF ! TEST_WELL_MIXED
      END DO ! I
    END DO ! J
  END DO ! K
END DO
!$OMP END DO

! ----------------------------------------------------------------------
! 2.1.4 Test WB_Ratio to see if well-mixed layer allowed
!-----------------------------------------------------------------------
!cdir collapse

!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( test_well_mixed(i,j) ) THEN

        wb_ratio(i,j) = wbn_int(i,j)/wbp_int(i,j)

        IF ( wb_ratio(i,j)  <=  dec_thres(i,j) ) THEN
            ! No need to test depth of mixing any further as
            ! well-mixed layer buoyancy flux integral criteria.
            ! SML will simply stay well-mixed (and so use defaults)
          ksurf_iterate(i,j)= .FALSE.
          ktop_iterate(i,j) = .FALSE.
          IF ( dsc(i,j) ) THEN
              ! Recouple DSC with SML:
              ! move surface driven entrainment
              ! RHOKH(z_i) = rho * w_e * DZL and w_e ~ 1/DB_TOP, so:
            IF ( db_top(i,j) >  0.0 .AND. db_dsct(i,j) >  0.01 ) THEN
                                          ! can't calc Zil. term
              rhokh_surf_ent(i,j) = rhokh_surf_ent(i,j) *               &
                     ( rho_uv(i,j,ntdsc(i,j)+1) * db_top(i,j) *         &
                       rdz(i,j,ntml(i,j)+1) ) /                         &
                     ( rho_uv(i,j,ntml(i,j)+1) * db_dsct(i,j) *         &
                                          rdz(i,j,ntdsc(i,j)+1) )
              rhokm(i,j,ntdsc(i,j)+1) = rhokm(i,j,ntml(i,j)+1) *        &
                     ( rho_tq(i,j,ntdsc(i,j)) * db_top(i,j) *           &
                       rdz(i,j,ntml(i,j)+1) ) /                         &
                     ( rho_tq(i,j,ntml(i,j)) * db_dsct(i,j) *           &
                                        rdz(i,j,ntdsc(i,j)+1) )
            END IF
              ! redesignate top-driven entrainment at ZHSC
              ! (ignore that calculated at ZH)
            rhokh_top_ent(i,j) = rhokh_dsct_ent(i,j)
            zh(i,j) = zhsc(i,j)
            ntml(i,j) = ntdsc(i,j)
            v_top(i,j) = v_top_dsc(i,j)
            zsml_base(i,j) = 0.1 * zh(i,j)
            zc(i,j) = zc_dsc(i,j)
            zhsc(i,j) = 0.0
            ntdsc(i,j) = 0
            v_top_dsc(i,j) = 0.0
            zdsc_base(i,j) = 0.0
            zc_dsc(i,j)    = 0.0
            ft_nt_zh(i,j)   = ft_nt_zhsc(i,j)
            ft_nt_zhsc(i,j) = 0.0
            fq_nt_zh(i,j)   = fq_nt_zhsc(i,j)
            fq_nt_zhsc(i,j) = 0.0
            dsc(i,j) = .FALSE.
            cumulus(i,j) = .FALSE.
            coupled(i,j) = .FALSE.
          END IF  ! recoupled DSC layer
        ELSE   ! buoyancy flux threshold violated
            !---------------------------------
            ! Extent of mixing must be reduced
            !---------------------------------
          IF ( .NOT.cumulus(i,j) ) ksurf_iterate(i,j) = .TRUE.
          ktop_iterate(i,j)  = .TRUE.
          IF (.NOT.dsc(i,j)) THEN
              ! Set up a `COUPLED' decoupled layer,
              !   implies no explicit `entrainment' at ZH.
              ! Note a new ZH (and thence NTML) will be calculated by
              ! wb integral iteration.
            IF (cumulus(i,j)) zk_uv=SQRT(zh(i,j)-1000000.)
                              ! APLTEST: shouldn't ever happen!
            dsc(i,j) = .TRUE.
            coupled(i,j) = .TRUE.
            ntdsc(i,j) = ntml(i,j)
            zhsc(i,j) = zh(i,j)
            zc_dsc(i,j) = zc(i,j)
            v_top_dsc(i,j) = v_top(i,j)
            v_sum_dsc(i,j) = v_sum(i,j)
            ft_nt_zhsc(i,j) = ft_nt_zh(i,j)
            fq_nt_zhsc(i,j) = fq_nt_zh(i,j)
              ! put all entrainment into RHOKH_TOP
            rhokh_dsct_ent(i,j) = rhokh_top_ent(i,j)                    &
                                + rhokh_surf_ent(i,j)
            rhokh_top_ent(i,j) = 0.0
            rhokh_surf_ent(i,j) = 0.0
            rhokm_top(i,j,ntml(i,j)+1) = rhokm_top(i,j,ntml(i,j)+1)     &
                                       + rhokm(i,j,ntml(i,j)+1)
            rhokm(i,j,ntml(i,j)+1) = 0.0
          END IF
        END IF   ! test on WB_RATIO LE DEC_THRES
      END IF   ! testing for well-mixed layer (TEST_WELL_MIXED)

    END DO ! I
  END DO ! J
!$OMP END DO

! ----------------------------------------------------------------------
! 2.2 Start iteration to find top of surface-driven mixing, ZSML_TOP,
!     within predetermined maximum and minimum height limits.
!     The solution is the height that gives WB_RATIO = DEC_THRES.
!     Procedure used makes 3 sweeps (up, down and up again), using
!     progressively smaller increments (Z_INC), each time stopping when
!     the buoyancy flux threshold or the height limits are reached.
!--------------------------------------------------------------------
!     If boundary layer is stable then ignore surface driven mixing.
!--------------------------------------------------------------------
!cdir collapse

!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF ( ksurf_iterate(i,j) ) THEN
          !-----------------------------------------------------------
          ! Mixing must extend just above surface layer
          ! (not clear precisely how to define this here: for now use
          !  KSURF calculated from ZH)
          !-----------------------------------------------------------
          ! limit K-surf to below cloud-base
        z_bot_lim(i,j)=z_uv(i,j,ksurf(i,j)+1)                           &
               + 0.1 * (z_uv(i,j,ksurf(i,j)+2)-z_uv(i,j,ksurf(i,j)+1))
          ! limit K-surf to below cloud-top radiatively cooled layer
        z_top_lim(i,j)=MAX( z_bot_lim(i,j), zhsc(i,j)- dzrad(i,j) )

        z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)
          !-----------------------------------------------------
          ! Initial increment to ZSML_TOP found by dividing
          ! up depth of layer within which it is allowed:
          ! Start with ZSML_TOP at lower limit and work upwards
          !-----------------------------------------------------
        z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j))                      &
                     / FLOAT(n_steps)
        zsml_top(i,j) = z_bot_lim(i,j)

        wb_ratio(i,j) = dec_thres(i,j) - 1.0 ! to be < DEC_THRES

      END IF ! KSURF_ITERATE

    END DO
  END DO
!$OMP END DO

!$OMP SINGLE 
  ij_len=pdims%i_end*pdims%j_end
  DO i = 1, ij_len
    to_do(i)    = .FALSE.
    ind_todo(i) = i
    up(i)       = 1
  END DO

  c_len=ij_len

  l=0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      l=l+1
      IF (ksurf_iterate(i,j)) THEN
        to_do(l)=.TRUE.
      END IF
    END DO
  END DO
!$OMP END SINGLE

  DO n_sweep = 1, 3


!$OMP MASTER

        ! Compress to_do and ind_todo (will have new length c_len)
! DEPENDS ON: excfnl_cci

    CALL excfnl_cci(c_len, to_do, ind_todo)

        ! Restart inner interation with the points of outer
    c_len_i = c_len
    todo_inner(1:c_len_i) = to_do(1:c_len_i)
    ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

!$OMP END MASTER
!$OMP BARRIER

    DO ns = 1, n_steps

!$OMP MASTER

          ! Calculate active elements and compress
! DEPENDS ON: excfnl_compin
      CALL excfnl_compin(up, wb_ratio, dec_thres, 1,                    &
                         c_len_i, ind_todo_i, todo_inner)

!$OMP END MASTER
!$OMP BARRIER

!cdir nodep
!$OMP DO SCHEDULE(STATIC)
      DO ic = 1, c_len_i
        j1=(ind_todo_i(ic)-1)/pdims%i_end+1
        i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

        zsml_top(i1,j1)=zsml_top(i1,j1)+z_inc(i1,j1)

            ! assume wb goes to zero at ZSML_TOP
        wb_surf_int(i1,j1) =                                            &
             bflux_surf(i1,j1) * z_tq(i1,j1,ksurf(i1,j1)) *             &
             ( 1.0 - z_tq(i1,j1,ksurf(i1,j1))/                          &
                                        (2.0*zsml_top(i1,j1)) )
        wb_surf_int(i1,j1) = MAX(1.0e-14,wb_surf_int(i1,j1))
            ! Note: WB_DZRAD_INT not included as K_SURF restricted
            !       to below zi-dzrad
        wbp_int(i1,j1) = wb_surf_int(i1,j1)  ! must be > 0
        wbn_int(i1,j1) = 0.0

        z_inv(i1,j1) = zsml_top(i1,j1)

      END DO ! ic c_len_i
!$OMP END DO

!..Integrate buoyancy flux profile given this ZSML_TOP

!$OMP DO SCHEDULE(STATIC)
      DO jj = 1, c_len_i, jblock
      DO k = 2, bl_levels-1
!cdir nodep
        DO ic = jj, MIN(jj+jblock-1, c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          IF ( k  >=  ksurf(i1,j1)+1 .AND.                              &
               k  <=  ntop(i1,j1) ) THEN

            z_pr = z_uv(i1,j1,k)
            IF (z_pr  <   z_inv(i1,j1)) THEN
              kh_surf(i1,j1,k) = w_h_top(i1,j1)                         &
                                    *z_pr*(1.0-z_pr/z_inv(i1,j1) )      &
                                         *(1.0-z_pr/z_inv(i1,j1) )
            ELSE
              kh_surf(i1,j1,k) = 0.0
            END IF
                !-----------------------------------------------------
                ! No F2 or FSc terms here because we're only
                ! considering effects driven from the surface
                !-----------------------------------------------------
            IF (z_cbase(i1,j1)  >   z_tq(i1,j1,k)) THEN
                  ! cloud-base above this range so use dry WB
              wb_scld= kh_surf(i1,j1,k) * db_ga_dry(i1,j1,k)
              wb_cld = 0.0
            ELSE IF (z_cbase(i1,j1)  <   z_tq(i1,j1,k-1)) THEN
                  ! cloud-base below this range so use cloudy WB
              wb_cld = kh_surf(i1,j1,k) * db_ga_cld(i1,j1,k)
              wb_scld=0.0
            ELSE
                  ! cloud-base within this integration range
                  ! so treat cloud and sub-cloud layer wb separately
              cld_frac = (z_tq(i1,j1,k)-z_cbase(i1,j1))                 &
                        /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac                                        &
                         * kh_surf(i1,j1,k)*db_ga_cld(i1,j1,k)
              wb_scld = (1.0-cld_frac)                                  &
                         * kh_surf(i1,j1,k)*db_ga_dry(i1,j1,k)
            END IF
            IF (wb_cld  >=  0.0) THEN
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_cld
            ELSE
              wbn_int(i1,j1)= wbn_int(i1,j1) - wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_scld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)- wb_scld
            END IF

          END IF ! K

        END DO ! ic c_len_i
      END DO ! K

      DO ic = jj, MIN(jj+jblock-1, c_len_i)
        j1=(ind_todo_i(ic)-1)/pdims%i_end+1
        i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
        wb_ratio(i1,j1) = wbn_int(i1,j1)/wbp_int(i1,j1)
      END DO ! ic c_len_i
   END DO
!$OMP END DO

    END DO  ! loop stepping up through ML (N_steps)

!cdir nodep
!$OMP DO SCHEDULE(STATIC)
    DO ic = 1, c_len
      l=ind_todo(ic)
      j1=(l-1)/pdims%i_end+1
      i1=l-(j1-1)*pdims%i_end

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to calculate WB for ZSML_TOP at a current Z_INC
      z_inc(i1,j1)= z_inc(i1,j1)/FLOAT(n_steps+1)

      IF ((up(l) == 1.AND.wb_ratio(i1,j1) >= dec_thres(i1,j1)).OR.      &
                 ! hit thres while working up
          (up(l) == 0.AND.wb_ratio(i1,j1) <= dec_thres(i1,j1))) THEN
                 ! hit thres while working down
        up(l) = 1-up(l)   ! change direction of sweep
        z_inc(i1,j1)= - z_inc(i1,j1)
      ELSE IF (zsml_top(i1,j1) >= z_top_lim(i1,j1)-1.0) THEN
                 ! hit upper height limit (give-or-take 1m) without
                 ! reaching threshold
        to_do(ic)=.FALSE.
        zsml_top(i1,j1) = z_top_lim(i1,j1)
      ELSE IF (zsml_top(i1,j1)  <=  z_bot_lim(i1,j1)+ 1.0) THEN
                 ! hit lower height limit (give-or-take 1m) without
                 ! reaching threshold
        to_do(ic)=.FALSE.
        zsml_top(i1,j1) = z_bot_lim(i1,j1)
      END IF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

    END DO ! c_len
!$OMP END DO

  END DO ! n_sweep

!$OMP DO SCHEDULE(DYNAMIC)
DO jj = pdims%j_start, pdims%j_end, jblock
  DO j = jj, MIN(jj+jblock-1, pdims%j_end)
    DO i = pdims%i_start,pdims%i_end
      ntml_new(i,j) = 2
      status_ntml(i,j)=.TRUE.
    END DO
  END DO

  DO k = 2, bl_levels-2
!cdir collapse
    DO j = jj, MIN(jj+jblock-1, pdims%j_end)
      DO i = pdims%i_start,pdims%i_end
        IF ( ksurf_iterate(i,j) .AND. status_ntml(i,j) ) THEN
            ! -------------
            ! find new NTML
            ! -------------
          IF  (z_uv(i,j,k+1)  <   zsml_top(i,j)) THEN
            ntml_new(i,j) = k+1
          ELSE
            status_ntml(i,j)=.FALSE.
          END IF
            ! --------------------------------------------------------
            ! Rounding error previously found to give
            !      ZSML_TOP > Z_TOP_LIM = ZHSC
            ! Test on ZSML_TOP hitting thresholds consequently changed
            ! but also include the following failsafe tests here.
            ! --------------------------------------------------------
          ntml(i,j) = MIN( ntdsc(i,j), ntml_new(i,j)-1 )
          zh(i,j)   = MIN(  zhsc(i,j), zsml_top(i,j) )

        END IF  ! KSURF_ITERATE true

      END DO
    END DO
  END DO
END DO
!$OMP END DO

! ----------------------------------------------------------------------
! 2.3 Now repeat the above procedure to find the base of the
!     top-driven K profile, ZDSC_BASE.
! ----------------------------------------------------------------------
!cdir collapse
!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF ( ktop_iterate(i,j) ) THEN

          ! Lower limit on base of DSC layer
        z_bot_lim(i,j) = 0.1 * zh(i,j)
          ! Upper limit on base of DSC layer
        z_top_lim(i,j) = zhsc(i,j) - dzrad(i,j)
          ! If cumulus limit base of top-driven mixing to above ZH
        IF ( cumulus(i,j) ) THEN
          z_bot_lim(i,j) = MIN( zh(i,j), z_top_lim(i,j) )
        END IF

        z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)

!..Divide up depth of layer within which ZDSC_BASE is allowed
        z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j))                      &
                       /FLOAT(n_steps)
        zdsc_base(i,j) = z_bot_lim(i,j)
                           ! will start at Z_BOT_LIM+Z_INC

        wb_ratio(i,j) = dec_thres(i,j) + 1.0 ! to be > DEC_THRES

      END IF ! KTOP_ITERATE

    END DO

    DO i = pdims%i_start,pdims%i_end
      IF ( ktop_iterate(i,j) .AND. wb_dzrad_int(i,j)  <   0.0 ) THEN
          !-------------------------------------------------------------
          ! Estimation of wb integral over radiatively cooled cloud-top
          ! region not yet performed (ie. DSC over stable surface) so
          ! do it now.
          !-------------------------------------------------------------
        z_inv(i,j)  = zhsc(i,j)
        z_cbase(i,j)= z_inv(i,j) - zc_dsc(i,j)
        cf_ml(i,j)  = cf_dsc(i,j)
        df_ctop(i,j)= df_dsct_over_cp(i,j)
        rho_we(i,j) = rho_we_dsc(i,j)
        dsl(i,j)    = dsl_dsc(i,j)
        dqw(i,j)    = dqw_dsc(i,j)
        tothf_zi(i,j)= - rho_we(i,j)*dsl(i,j) + ft_nt_zhsc(i,j)
        totqf_zi(i,j)= - rho_we(i,j)*dqw(i,j) + fq_nt_zhsc(i,j)
        ntop(i,j)   = ntdsc(i,j) - 1
        dzrad(i,j)  = 100.0
        IF ( ntop(i,j)  >   1 ) THEN
          DO WHILE ( z_tq(i,j,ntop(i,j))  >   z_inv(i,j)-dzrad(i,j)     &
               .AND. z_tq(i,j,ntop(i,j)-1)  >   z_inv(i,j)-z_cbase(i,j))
            ntop(i,j) = ntop(i,j) - 1
            IF (ntop(i,j) == 1 ) THEN
              EXIT
            END IF
          END DO
        END IF
        dzrad(i,j) = z_inv(i,j) - z_tq(i,j,ntop(i,j))

        wsl_dzrad_int = dzrad(i,j) *                                    &
                        ( 0.66*df_ctop(i,j) - rho_we(i,j)*dsl(i,j) )
        wqw_dzrad_int = - dzrad(i,j) * rho_we(i,j) * dqw(i,j)

        wb_dzrad_int(i,j) = wsl_dzrad_int * (                           &
                              (1.0-cf_ml(i,j))*btm(i,j,ntop(i,j)+1) +   &
                                cf_ml(i,j)*btm_cld(i,j,ntop(i,j)+1) )   &
                          + wqw_dzrad_int * (                           &
                              (1.0-cf_ml(i,j))*bqm(i,j,ntop(i,j)+1) +   &
                                cf_ml(i,j)*bqm_cld(i,j,ntop(i,j)+1) )
        wb_dzrad_int(i,j) = g * wb_dzrad_int(i,j)
        wb_dzrad_int(i,j) = MAX( 0.0, wb_dzrad_int(i,j) )
      END IF

      wb_dzrad_int(i,j) = MAX( 1.0e-14, wb_dzrad_int(i,j) )

    END DO ! I
  END DO ! J
!$OMP END DO

!$OMP MASTER

  ij_len=pdims%i_end*pdims%j_end
  DO i = 1, ij_len
    to_do(i)    = .FALSE.
    ind_todo(i) = i
    up(i)     = 1
  END DO

  c_len=ij_len

  l=0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      l=l+1
      IF (ktop_iterate(i,j)) THEN
        to_do(l)=.TRUE.
      END IF
    END DO
  END DO

!$OMP END MASTER
!$OMP BARRIER

  DO n_sweep = 1, 3

!$OMP MASTER

        ! Compress to_do and ind_todo (will have new length c_len)
! DEPENDS ON: excfnl_cci
    CALL excfnl_cci(c_len, to_do, ind_todo)

        ! Restart inner interation with the points of outer
    c_len_i = c_len
    todo_inner(1:c_len_i) = to_do(1:c_len_i)
    ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)

!$OMP END MASTER
!$OMP BARRIER

    DO ns = 1, n_steps

!$OMP MASTER

          ! Calculate active elements and compress
! DEPENDS ON: excfnl_compin
      CALL excfnl_compin(up, wb_ratio, dec_thres, 2,                    &
                         c_len_i, ind_todo_i, todo_inner)
!$OMP END MASTER
!$OMP BARRIER

!cdir nodep
!$OMP DO SCHEDULE(STATIC)
      DO ic = 1, c_len_i
        j1=(ind_todo_i(ic)-1)/pdims%i_end+1
        i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

        zdsc_base(i1,j1) = zdsc_base(i1,j1)+z_inc(i1,j1)
        scbase(i1,j1)    = .TRUE.
            ! Flag for NBDSC found (only needed for LockWhelan2006)
        IF (flux_grad  ==  LockWhelan2006) scbase(i1,j1) = .FALSE.
        ft_nt_dscb(i1,j1)= 0.0
        fq_nt_dscb(i1,j1)= 0.0

        wbn_int(i1,j1) = 0.0
        wbp_int(i1,j1) = wb_dzrad_int(i1,j1)
        IF ( ksurf_iterate(i1,j1) .AND.                                 &
             zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) THEN
              ! only include surface flux if K_SURF is included
              ! in the wb calculation and K profiles overlap
          wbp_int(i1,j1) = wbp_int(i1,j1) + wb_surf_int(i1,j1)
          wbn_int(i1,j1) = 0.0
        END IF
        zinv_pr(i1,j1) = zhsc(i1,j1)-zdsc_base(i1,j1)

      END DO ! ic c_len_i
!$OMP END DO

!..Integrate buoyancy flux profile given this ZDSC_BASE

!$OMP DO SCHEDULE(STATIC)
      DO jj = 1, c_len_i, jblock
      DO k = 2, bl_levels-1
!cdir nodep
        DO ic = jj, MIN(jj+jblock-1, c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          IF ((k >= ksurf(i1,j1)+1).AND.(k <= ntop(i1,j1))) THEN

            khtop(i1,j1) = 0.0
            f2           = 0.0
            fsc          = 0.0
            IF (.NOT. scbase(i1,j1) ) THEN
              ft_nt_dscb(i1,j1) = ft_nt(i1,j1,k)
              fq_nt_dscb(i1,j1) = fq_nt(i1,j1,k)
            END IF
            z_pr = z_uv(i1,j1,k) - zdsc_base(i1,j1)

            IF (z_pr >   0.0 .AND.z_pr <  zinv_pr(i1,j1)) THEN

              IF (.NOT. scbase(i1,j1) ) THEN
                scbase(i1,j1) = .TRUE.
                z_ratio = (zdsc_base(i1,j1)-z_uv(i1,j1,k-1))            &
                         /(z_uv(i1,j1,k)-z_uv(i1,j1,k-1))
                ft_nt_dscb(i1,j1) = ft_nt(i1,j1,k-1) +                  &
                      (ft_nt(i1,j1,k)-ft_nt(i1,j1,k-1))*z_ratio
                fq_nt_dscb(i1,j1) = fq_nt(i1,j1,k-1) +                  &
                      (fq_nt(i1,j1,k)-fq_nt(i1,j1,k-1))*z_ratio
              END IF

              z_ratio = z_pr/zinv_pr(i1,j1)

              IF (flux_grad  ==  LockWhelan2006) THEN
                khtop(i1,j1) = 3.6 * vkman * rho_uv(i1,j1,k)            &
                             * v_top_dsc(i1,j1) * zinv_pr(i1,j1)        &
                  * (z_ratio**3) * ( (1.0-z_ratio)*(1.0-z_ratio) )
                f2 = rho_uv(i1,j1,k) * 0.5 * z_ratio                    &
                                     * 2.0**( z_ratio**4 )
                IF ( v_sum_dsc(i1,j1)  >   0.0 ) THEN
                  fsc = 3.5 * rho_uv(i1,j1,k)                           &
                           * (v_top_dsc(i1,j1)/v_sum_dsc(i1,j1))        &
                           * (z_ratio**3) * (1.0-z_ratio)
                END IF
              ELSE
                khtop(i1,j1) = g1 * vkman * rho_uv(i1,j1,k)             &
                       * v_top_dsc(i1,j1) * (( 1.0 - z_ratio )**0.8)    &
                       * z_pr * z_ratio
              END IF

            END IF

            khsurf(i1,j1) = 0.0
            IF ( zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) THEN
                  ! only include K_surf if profiles overlap
                  ! otherwise layers are independent
              khsurf(i1,j1) = kh_surf(i1,j1,k)
            END IF

            IF (flux_grad  ==  LockWhelan2006) THEN
              wslng = (f2+fsc)*(tothf_zi(i1,j1)-ft_nt_dscb(i1,j1))      &
                          - ( ft_nt(i1,j1,k)-ft_nt_dscb(i1,j1) )
              wqwng = (f2+fsc)*(totqf_zi(i1,j1)-fq_nt_dscb(i1,j1))      &
                          - ( fq_nt(i1,j1,k)-fq_nt_dscb(i1,j1) )
            END IF

            IF ( z_tq(i1,j1,k)  <=  z_cbase(i1,j1) ) THEN
                  ! Completely below cloud-base so use cloud-free form
              wb_scld = khsurf(i1,j1)* db_ga_dry(i1,j1,k) +             &
                        khtop(i1,j1) * db_noga_dry(i1,j1,k)
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_scld = wb_scld + ( g/rdz(i1,j1,k) ) *                &
                   ( btm(i1,j1,k-1)*wslng + bqm(i1,j1,k-1)*wqwng )
              END IF
              wb_cld  = 0.0
            ELSE IF (z_tq(i1,j1,k-1)  >=  z_cbase(i1,j1)) THEN
                  ! Completely above cloud-base so use cloudy formula
              wb_cld = ( khsurf(i1,j1) * db_ga_cld(i1,j1,k) +           &
                         khtop(i1,j1)  * db_noga_cld(i1,j1,k) )
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_cld = wb_cld + ( g/rdz(i1,j1,k) ) * (                &
                     ( btm(i1,j1,k-1)*(1.0-cf_ml(i1,j1)) +              &
                       btm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wslng +        &
                     ( bqm(i1,j1,k-1)*(1.0-cf_ml(i1,j1)) +              &
                       bqm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wqwng )
              END IF
              wb_scld = 0.0
            ELSE
                  ! cloud-base within this integration range
                  ! so treat cloud and sub-cloud layer wb separately
              wb_scld = khsurf(i1,j1) * db_ga_dry(i1,j1,k) +            &
                        khtop(i1,j1) * db_noga_dry(i1,j1,k)
              wb_cld = ( khsurf(i1,j1) * db_ga_cld(i1,j1,k) +           &
                         khtop(i1,j1)  * db_noga_cld(i1,j1,k) )
              IF (flux_grad  ==  LockWhelan2006) THEN
                wb_scld = wb_scld + ( g/rdz(i1,j1,k) ) *                &
                   ( btm(i1,j1,k-1)*wslng + bqm(i1,j1,k-1)*wqwng )
                wb_cld = wb_cld + ( g/rdz(i1,j1,k) ) * (                &
                     ( btm(i1,j1,k-1)*(1.0-cf_ml(i1,j1)) +              &
                       btm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wslng +        &
                     ( bqm(i1,j1,k-1)*(1.0-cf_ml(i1,j1)) +              &
                       bqm_cld(i1,j1,k-1)*cf_ml(i1,j1) )*wqwng )
              END IF
              cld_frac = (z_tq(i1,j1,k)-z_cbase(i1,j1))                 &
                        /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac * wb_cld
              wb_scld = (1.0-cld_frac) * wb_scld
            END IF
            IF (wb_cld  >=  0.0) THEN
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_cld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_scld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_scld
            END IF

          END IF ! K
        END DO ! ic c_len_i
      END DO ! K
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      DO ic = 1, c_len_i
        j1=(ind_todo_i(ic)-1)/pdims%i_end+1
        i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
        wb_ratio(i1,j1)=wbn_int(i1,j1)/wbp_int(i1,j1)
      END DO ! ic c_len_i
!$OMP END DO

    END DO  ! loop stepping up through ML

!cdir nodep
!$OMP DO SCHEDULE(STATIC)
    DO ic = 1, c_len
      l=ind_todo(ic)
      j1=(l-1)/pdims%i_end+1
      i1=l-(j1-1)*pdims%i_end

!..sub-divide current Z_INC into one more part than there will be steps
!..as there is no need to recalculate WB at a current Z_INC
      z_inc(i1,j1)= z_inc(i1,j1)/FLOAT(n_steps+1)

      IF (                                                              &
         (up(l) == 1 .AND.wb_ratio(i1,j1) <= dec_thres(i1,j1)).OR.      &
                ! hit thres while working up
         (up(l) == 0 .AND.wb_ratio(i1,j1) >= dec_thres(i1,j1))) THEN
                ! hit thres while working down
        up(l) = 1-up(l)   ! change direction of sweep
        z_inc(i1,j1)=- z_inc(i1,j1)
      ELSE IF ( zdsc_base(i1,j1) >= z_top_lim(i1,j1)-1.0 .OR.           &
                zdsc_base(i1,j1) <=  z_bot_lim(i1,j1)+1.0 ) THEN
            ! hit height limits (give-or-take 1m) without
            ! reaching threshold
        to_do(ic)=.FALSE.
      END IF

!..Note that if the threshold has not been passed then the next sweep
!..continues in the same direction (but with reduced increment).

    END DO ! c_len
!$OMP END DO

  END DO  ! loop over sweeps

! Convert integrated WB to profiles of WB itself for diagnostics
! ----------------------------------------------------------------------
! 2.4 Set depth of cloud-top driven mixing in SML when there is a DSC
!     layer above (eg. fog under Sc) to be the SML layer depth
!     and, if option selected, the top of any K profile above the LCL 
!     in cumulus to be the lower of the parcel top, the DSC base and 
!     500m above the LCL but at least 1.1*z_lcl
! ----------------------------------------------------------------------
IF (entr_smooth_dec == on) THEN

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( cumulus(i,j) .OR.                                            &
           ( dsc(i,j) .AND. zdsc_base(i,j) < zh(i,j) ) ) THEN
          ! ignore SML `cloud-top' driven mixing
        zsml_base(i,j) = zh(i,j)
        v_top(i,j)     = 0.0
      ELSE
        zsml_base(i,j) = 0.1*zh(i,j)
      END IF
    END DO  ! loop over j
  END DO  ! loop over I
!$OMP END DO

ELSE ! entr_smooth_dec off

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      IF ( cumulus(i,j) .OR. coupled(i,j) ) THEN
          ! ignore SML `cloud-top' driven mixing
        zsml_base(i,j) = zh(i,j)
        v_top(i,j)     = 0.0
      ELSE
        zsml_base(i,j) = 0.1*zh(i,j)
      END IF
    END DO  ! loop over j
  END DO  ! loop over I
!$OMP END DO

END IF  ! test on entr_smooth_dec

IF ( Kprof_cu /= off ) THEN

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      ! Start with constant depth (input from namelist)
      Kprof_cu_top(i,j) = zh(i,j)+max_cu_depth
      ! Keep top of Kprof for cu below DSC base:
      IF (dsc(i,j)) THEN
        Kprof_cu_top(i,j) = MIN( Kprof_cu_top(i,j), zdsc_base(i,j) )
      END IF
      ! ...but at least 1.1*z_lcl:
      Kprof_cu_top(i,j) = MAX( 1.1*zh(i,j), Kprof_cu_top(i,j) )
      ! ...but no higher than parcel top:
      Kprof_cu_top(i,j) = MIN( zhpar(i,j), Kprof_cu_top(i,j) )

      cu_depth_scale(i,j) = (Kprof_cu_top(i,j)-zh(i,j))/3.0

      rhokh_lcl(i,j) = 0.0
      IF (kprof_cu == klcl_entr) THEN
        ! Use the BL entrainment parametrization as calculated above
        rhokh_lcl(i,j) = MIN( rhokh_surf_ent(i,j), 5.0)
      END IF  ! no other options coded yet
    END DO  ! loop over j
  END DO  ! loop over I
!$OMP END DO

END IF  ! test on Kprof_cu
!-----------------------------------------------------------------------
! 3.  Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment
!     level is approached.
!-----------------------------------------------------------------------
!cdir collapse

!$OMP DO SCHEDULE(DYNAMIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k=ntml(i,j)+1
      kh_top_factor(i,j) = MAX( 0.7 , 1.0 - SQRT(                       &
               rhokh_surf_ent(i,j) /                                    &
                     ( rho_uv(i,j,k)*w_h_top(i,j)*vkman*zh(i,j) ) ) )
      km_top_factor(i,j) = MAX( 0.7 , 1.0 - SQRT( rhokm(i,j,k) /        &
                 ( rho_tq(i,j,k-1)*w_m_top(i,j)*vkman*zh(i,j) ) ) )
      scdepth(i,j) = zh(i,j) - zsml_base(i,j)
      factor = g1 * rho_uv(i,j,k) * v_top(i,j) *vkman *scdepth(i,j)
      IF ( factor  >   0.0) THEN
        kh_sct_factor(i,j) = 1.0 -                                      &
                             ( rhokh_top_ent(i,j) / factor )**1.25
                                                          ! 1.25=1/0.8
      ELSE
        kh_sct_factor(i,j) = 1.0
      END IF
      factor = g1 * rho_tq(i,j,k-1) * v_top(i,j) *                      &
                      vkman * scdepth(i,j) * 0.75
      IF ( factor  >   0.0) THEN
        km_sct_factor(i,j) = 1.0 -                                      &
                             ( rhokm_top(i,j,k) / factor )**1.25
                                                          ! 1.25=1/0.8
      ELSE
        km_sct_factor(i,j) = 1.0
      END IF

      IF (ntdsc(i,j)  >   0) THEN
        !-------------------------------------------------------------
        ! Set up factors to ensure K profile continuity at ZHSC;
        ! no need to limit size of factor as precise shape of top-down
        ! mixing profile not important.
        ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus
        ! layer exists, i.e. NTDSC > 0.
        !-------------------------------------------------------------
        k=ntdsc(i,j)+1
        dscdepth(i,j) = zhsc(i,j) - zdsc_base(i,j)
        factor = g1*rho_uv(i,j,k)*v_top_dsc(i,j)*vkman*dscdepth(i,j)
        IF ( factor  >   0.0) THEN
          kh_dsct_factor(i,j) = 1.0 -                                   &
                              ( rhokh_dsct_ent(i,j) / factor )**1.25
                                                          ! 1.25=1/0.8
        ELSE
          kh_dsct_factor(i,j) = 1.0
        END IF

        factor = 0.75 * g1 * rho_tq(i,j,k-1) * v_top_dsc(i,j) *         &
                             vkman * dscdepth(i,j)
        IF ( factor  >   0.0) THEN
          km_dsct_factor(i,j) = 1.0 -                                   &
                             ( rhokm_top(i,j,k) / factor )**1.25
                                                          ! 1.25=1/0.8
        ELSE
          km_dsct_factor(i,j) = 1.0
        END IF
      END IF
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 4.  Calculate height dependent turbulent
!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------

! Reset identifiers of base of decoupled layer mixing

!$OMP MASTER

!cdir collapse
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      scbase(i,j) = .FALSE.
      nbdsc(i,j)  = 0
    END DO
  END DO

!$OMP END MASTER
!$OMP BARRIER

!-------------------------------------------------------------
! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles
! for the surface mixed layer.
! This is a variation on an up-side-down version of the cubic
! surface-forced profiles below.  Implement between at least
! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
! Note this may well include NTML+1: entrainment fluxes will
! be dealt with in KMKHZ.
!-------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, jblock
  DO k = 2, bl_levels
!cdir collapse
    DO j = jj, MIN(jj+jblock-1, pdims%j_end)
      DO i = pdims%i_start,pdims%i_end

!           Calculate the height of u,v-level above the surface
! *APL: z0m removed from z in K(z)
        zk_uv = z_uv(i,j,k)

!           Calculate the height of T,q-level above the surface

        zk_tq = z_tq(i,j,k-1)

        IF ( zk_uv  <   zh(i,j) .AND.                                   &
             zk_uv  >   zsml_base(i,j) ) THEN
          z_pr  = zk_uv - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          z_ratio = z_pr/zh_pr
          IF (flux_grad  ==  LockWhelan2006) THEN

            rhokh_top(i,j,k) = 3.6*vkman * rho_uv(i,j,k) * v_top(i,j)   &
                               * zh_pr * (z_ratio**3)                 &
                               * (( 1.0 - z_ratio )**2)

            IF ( .NOT.coupled(i,j) ) THEN
              rhof2(i,j,k)  = rho_uv(i,j,k) * 0.5 * z_ratio             &
                                          * 2.0**( z_ratio**4 )
              IF ( v_sum(i,j)  >   0.0 ) THEN
                rhofsc(i,j,k) = 3.5 * rho_uv(i,j,k)                     &
                                * (v_top(i,j)/v_sum(i,j))               &
                                * (z_ratio**3) * (1.0-z_ratio)
              END IF
            END IF

          ELSE  ! Not LockWhelan2006

            rhokh_top(i,j,k) = rho_uv(i,j,k) * v_top(i,j) * g1 *        &
              vkman * ( ( 1.0 - kh_sct_factor(i,j)*z_ratio )**0.8 )     &
                                             * z_pr * z_ratio
          END IF

        END IF
          !-------------------------------------------------------
          !   For LockWhelan2006, KM_TOP could be changed to match
          !   the shape of KH_TOP.  This has not been done on the
          !   grounds that the change in shape arises with the
          !   inclusion of the other non-gradient terms.
          !-------------------------------------------------------
        IF ( zk_tq  <   zh(i,j) .AND.                                   &
             zk_tq  >   zsml_base(i,j) ) THEN
          z_pr = zk_tq - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          rhokm_top(i,j,k) = 0.75 * rho_tq(i,j,k-1) * v_top(i,j) *      &
                g1 * vkman *                                            &
                ( ( 1.0 - km_sct_factor(i,j)*z_pr/zh_pr )**0.8 )        &
                                           * z_pr * z_pr / zh_pr
                                                      ! PRANDTL=0.75
        END IF
          !-------------------------------------------------------------
          ! Add contribution to top-down mixing coefficient
          ! profiles for decoupled stratocumulus layers when
          ! one exists
          !-------------------------------------------------------------
        IF ( zk_uv  <   zhsc(i,j) .AND.                                 &
                zk_uv  >   zdsc_base(i,j) ) THEN
          IF (.NOT. scbase(i,j) ) THEN
            scbase(i,j) = .TRUE.
              ! identifies lowest layer below which there is mixing
            nbdsc(i,j) = k
          END IF
            !-----------------------------------------------------------
            ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing
            ! profiles and add to any generated in the surface mixing
            ! layer.
            ! This is a variation on an up-side-down version of the
            ! cubic surface-forced profiles above.  Implement between
            ! at least the top of the `surface layer' (at Z=0.1*ZH) and
            ! ZHSC.
            !-----------------------------------------------------------
          z_pr = zk_uv - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          z_ratio = z_pr/zh_pr


          IF (flux_grad  ==  LockWhelan2006) THEN

            rhokh_top(i,j,k) = 3.6*vkman * rho_uv(i,j,k)                &
                             * v_top_dsc(i,j) * zh_pr * (z_ratio**3)  &
                             * (( 1.0 - z_ratio )**2)

            rhof2(i,j,k)  = rho_uv(i,j,k) * 0.5 * z_ratio               &
                                          * 2.0**( z_ratio**4 )
            IF ( v_sum_dsc(i,j)  >   0.0 ) THEN
              rhofsc(i,j,k) = 3.5 * rho_uv(i,j,k)                       &
                                * (v_top_dsc(i,j)/v_sum_dsc(i,j))       &
                                * (z_ratio**3) * (1.0-z_ratio)
            END IF

          ELSE  ! Not LockWhelan2006

            rhokh_top(i,j,k) = rhokh_top(i,j,k) +                       &
               rho_uv(i,j,k)*v_top_dsc(i,j)*g1*vkman*                   &
                  ( ( 1.0 - kh_dsct_factor(i,j)*z_ratio )**0.8 )        &
                                             * z_pr * z_ratio
          END IF
        END IF
          !-------------------------------------------------------------
          ! Now momentum
          !-------------------------------------------------------------
        IF ( zk_tq  <   zhsc(i,j) .AND.                                 &
             zk_tq  >   zdsc_base(i,j) ) THEN
          z_pr = zk_tq - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          rhokm_top(i,j,k) = rhokm_top(i,j,k) +                         &
             0.75*rho_tq(i,j,k-1)*v_top_dsc(i,j)*g1*vkman*              &
                ( ( 1.0 - km_dsct_factor(i,j)*z_pr/zh_pr )**0.8 )       &
                                        * z_pr * z_pr / zh_pr
        END IF

      END DO
    END DO
  END DO
END DO
!$OMP END DO

!----------------------------------------------------
! Now K_SURF profiles
!----------------------------------------------------
  IF (flux_grad  ==  LockWhelan2006) THEN
!----------------------------------------------------
! Lock and Whelan formulation
!----------------------------------------------------

    c_ws = 0.42     ! ~ PR_NEUT^3 by design
    pr_neut = 0.75
    pr_conv = 0.6

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, jblock
    DO k = 2, bl_levels
!cdir collapse
      DO j = jj, MIN(jj+jblock-1, pdims%j_end)
        DO i = pdims%i_start,pdims%i_end

          zk_uv = z_uv(i,j,k)
          zk_tq = z_tq(i,j,k-1)

          IF (fb_surf(i,j)  >=  0.0) THEN

!           Calculate the free-convective scaling velocity at z(k)

            IF (coupled(i,j)) THEN  !  coupled
              wstar3 = zhsc(i,j) * fb_surf(i,j)
            ELSE
              wstar3 = zh(i,j) * fb_surf(i,j)
            END IF

            IF (zk_uv  <=  0.1*zh(i,j)) THEN
!             Surface layer calculation
              w_s_cubed_uv = 10.*c_ws * zk_uv * fb_surf(i,j)
            ELSE
!             Outer layer calculation
              w_s_cubed_uv = c_ws * wstar3
            END IF

            IF (zk_tq  <=  0.1*zh(i,j)) THEN
!             Surface layer calculation
              w_s_cubed_tq = 10.*c_ws * zk_tq * fb_surf(i,j)
            ELSE
!             Outer layer calculation
              w_s_cubed_tq = c_ws * wstar3
            END IF

!           Turbulent velocity scale for scalars

            w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_uv )    &
                                   **(1.0/3.0)
            w_h_uv = w_m_neut/pr_neut

            ! Also calc on TQ levels for W_M_TQ
            w_m_neut = ( v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq )    &
                                   **(1.0/3.0)
            w_h_tq = w_m_neut/pr_neut

!           Turbulent Prandtl number and velocity scale for scalars

            Prandtl = pr_neut*                                          &
            ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                     &
             (1.0/(c_ws*25.0))*w_s_cubed_tq*w_m_neut ) /                &
            ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                     &
             (1.0/(c_ws*25.0))*(pr_neut/pr_conv)*w_s_cubed_tq*w_m_neut )

            w_m_tq = Prandtl * w_h_tq

            IF ( zk_uv  <   zh(i,j) ) THEN
              !---------------------------------------------------------
              ! Calculate RHOKH(w_h,z/z_h)
              !---------------------------------------------------------

              rhokh(i,j,k) = rho_uv(i,j,k) * w_h_uv * vkman * zk_uv *   &
                                    ( 1.0 - ( zk_uv / zh(i,j) ) ) *     &
                                    ( 1.0 - ( zk_uv / zh(i,j) ) )

            END IF
            IF ( zk_tq  <   zh(i,j) ) THEN
              !---------------------------------------------------------
              ! Calculate RHOKM(w_m,z/z_h)
              !---------------------------------------------------------

              rhokm(i,j,k) = rho_tq(i,j,k-1) * w_m_tq * vkman * zk_tq * &
                                    ( 1.0 - ( zk_tq / zh(i,j) ) ) *     &
                                    ( 1.0 - ( zk_tq / zh(i,j) ) )

            END IF
          END IF
        END DO
      END DO
    END DO
 END DO
!$OMP END DO

  ELSE
!----------------------------------------------------
! Lock et al and Holtstalg and Boville formulations
!----------------------------------------------------

      ! Default to Lock et al
    c_ws = 0.25
    pr_neut = 0.75
    pr_conv = 0.375
    IF (flux_grad  ==  HoltBov1993) THEN
      c_ws = 0.6
      pr_neut = 1.0
      pr_conv = 0.6
    END IF

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, jblock
    DO k = 2, bl_levels
!cdir collapse
      DO j = jj, MIN(jj+jblock-1, pdims%j_end)
        DO i = pdims%i_start,pdims%i_end

!         Calculate the height of u,v-level above the surface
! *APL: z0m removed from z in K(z)
          zk_uv = z_uv(i,j,k)

!         Calculate the height of T,q-level above the surface

          zk_tq = z_tq(i,j,k-1)

          IF (fb_surf(i,j)  >=  0.0) THEN

!           Calculate the free-convective scaling velocity at z(k)

            IF (zk_uv  <=  0.1*zh(i,j)) THEN

!             Surface layer calculation

              w_s_cubed_uv = 10.*c_ws * zk_uv * fb_surf(i,j)
            ELSE

!             Outer layer calculation

              IF (coupled(i,j)) THEN  !  coupled and cloudy
                w_s_cubed_uv = c_ws * zhsc(i,j) * fb_surf(i,j)
              ELSE
                w_s_cubed_uv = c_ws * zh(i,j) * fb_surf(i,j)
              END IF
            END IF

            IF (zk_tq  <=  0.1*zh(i,j)) THEN

!             Surface layer calculation

              w_s_cubed_tq = 10.*c_ws * zk_tq * fb_surf(i,j)
            ELSE

!             Outer layer calculation

              IF (coupled(i,j)) THEN  !  coupled and cloudy
                w_s_cubed_tq = c_ws * zhsc(i,j) * fb_surf(i,j)
              ELSE
                w_s_cubed_tq = c_ws * zh(i,j) * fb_surf(i,j)
              END IF
            END IF

!           Turbulent velocity scale for momentum

            w_m_uv = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_uv)        &
                                   **(1.0/3.0)

            w_m_tq = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq)        &
                                   **(1.0/3.0)

!           Turbulent Prandtl number and velocity scale for scalars

            Prandtl = pr_neut*( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +   &
               (1.0/(c_ws*25.0))*w_s_cubed_uv*w_m_uv ) /                &
                              ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +   &
               (1.0/(c_ws*25.0))*(pr_neut/pr_conv)*w_s_cubed_uv*w_m_uv )
            w_h_uv = w_m_uv / Prandtl

            IF ( zk_uv  <   zh(i,j) ) THEN
              !---------------------------------------------------------
              ! Calculate RHOKH(w_h,z/z_h)
              !---------------------------------------------------------

              rhokh(i,j,k) = rho_uv(i,j,k) * w_h_uv * vkman * zk_uv *   &
                  ( 1.0 - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) ) *  &
                  ( 1.0 - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) )

            ELSE IF ( Kprof_cu /= off .AND. cumulus(i,j) .AND.          &
                      zk_uv < Kprof_cu_top(i,j) ) THEN
!             ! Exponential decay from ZH but tends to zero 
!             ! at Kprof_cu_top
              rhokh(i,j,k) = rhokh_lcl(i,j) *                           &
                       EXP(-(zk_uv-zh(i,j))/cu_depth_scale(i,j)) *      &
                       (1.0-(zk_uv-zh(i,j))/(Kprof_cu_top(i,j)-zh(i,j)))
            END IF
            IF ( zk_tq  <   zh(i,j) ) THEN
              !---------------------------------------------------------
              ! Calculate RHOKM(w_m,z/z_h)
              !---------------------------------------------------------

              rhokm(i,j,k) = rho_tq(i,j,k-1) * w_m_tq * vkman * zk_tq * &
                  ( 1.0 - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) ) *  &
                  ( 1.0 - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) )

            ELSE IF ( Kprof_cu /= off .AND. cumulus(i,j) .AND.          &
                      zk_tq < Kprof_cu_top(i,j) ) THEN
!             ! Exponential decay from ZH but tends to zero
!             ! at Kprof_cu_top
              rhokm(i,j,k) = prandtl_top(i,j) * rhokh_lcl(i,j) *        &
                       EXP(-(zk_tq-zh(i,j))/cu_depth_scale(i,j)) *      &
                       (1.0-(zk_tq-zh(i,j))/(Kprof_cu_top(i,j)-zh(i,j)))
            END IF
          END IF
        END DO
      END DO
    END DO
 END DO
!$OMP END DO

  END IF  ! Test on Flux_grad

  IF ( ng_stress  ==  BrownGrant97 .OR.                                 &
       ng_stress  ==  BrownGrant97_limited ) THEN

!$OMP DO SCHEDULE(STATIC)
  DO jj = pdims%j_start, pdims%j_end, jblock
    DO k = 2, bl_levels
!cdir collapse
      DO j = jj, MIN(jj+jblock-1, pdims%j_end)
        DO i = pdims%i_start,pdims%i_end
          zk_tq = z_tq(i,j,k-1)   ! stresses are calc on theta-levs
          IF ( fb_surf(i,j)  >   0.0 .AND. zk_tq  <   zh(i,j) ) THEN
              !---------------------------------------------------------
              ! Calculate non-gradient stress function
              ! (Brown and Grant 1997)
              ! Shape function chosen such that non-gradient stress
              ! goes to zero at 0.1*ZH and ZH
              !---------------------------------------------------------
            IF ( zk_tq  >   0.1*zh(i,j) ) THEN
              z_pr = zk_tq - 0.1*zh(i,j)
              zh_pr = 0.9*zh(i,j)
                !
                ! Outer layer calculation
                !
              IF (coupled(i,j)) THEN  !  coupled and cloudy
                wstar3 = zhsc(i,j) * fb_surf(i,j)
              ELSE
                wstar3 =   zh(i,j) * fb_surf(i,j)
              END IF

                ! Use the Holtslag and Boville velocity scale for
                ! non-gradient stress stability dependence, as in BG97
              w_m_hb_3 = v_s(i,j)*v_s(i,j)*v_s(i,j) + 0.6*wstar3
              f_ngstress(i,j,k) = ( rho_tq(i,j,k-1)/rhostar_gb(i,j) )   &
                * s_m * ( a_ngs * wstar3 / w_m_hb_3 )                   &
                   * ( z_pr / zh_pr ) * ( 1.0 -  ( z_pr / zh_pr ) ) *   &
                                        ( 1.0 -  ( z_pr / zh_pr ) )
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL


  IF (lhook) CALL dr_hook('EXCF_NL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE excf_nl
