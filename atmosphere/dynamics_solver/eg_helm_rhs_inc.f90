! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_helm_rhs_inc_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_helm_rhs_inc(R_u, R_v, R_w, R_p, R_theta, R_rho,           &
                                  R_etadot, rho_divu_out,                      &
                                  R_m_v_d, R_m_cl_d, R_m_cf_d,                 &
                                  R_m_r_d, R_m_gr_d, R_m_cf2_d,                &
                                  u, v, w, exner, exner_star, rho,             &
                                  thetav, etadot,m_v, m_cl, m_cf, m_r, m_gr,   &
                                  m_cf2, R_u_d, R_v_d, R_w_d, R_theta_d,       &
                                  R_rho_d, Ih, row_length, rows, n_rows,       &
                                  model_levels, model_domain, inner_it_val,    &
                                  S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf,    &
                                  S_m_cf2,S_m_r,S_m_gr, RHS,psi_w_surf,        &
                                  psi_w_lid)
 
      USE timestep_mod,      ONLY : timestep
      USE proc_info_mod,     ONLY : at_extremity, gc_proc_row_group
      USE level_heights_mod, ONLY : eta_theta_levels
      USE um_parvars,        ONLY : offx, offy, halo_i, halo_j
      USE eg_alpha_mod,      ONLY : alpha_u,     alpha_v,    alpha_w,          &
                                    alpha_theta, alpha_rho, alpha_p

      USE yomhook,           ONLY : lhook, dr_hook
      USE parkind1,          ONLY : jprb, jpim
      USE eg_vert_damp_mod,  ONLY : mu_w
      USE eg_helmholtz_mod
      USE eg_sisl_init_mod
      USE eg_v_at_poles_mod
      USE integrity_mod
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE UM_ParParams
      USE metric_terms_mod
      USE update_moisture_mod
      USE eg_dxout_mod

      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE gravity_mod
      USE eg_parameters_mod,   ONLY: pole_consts
      USE atmos_constants_mod, ONLY: p_zero, R
      USE proc_info_mod,       ONLY: n_proc, mype=>me
      
      USE PrintStatus_mod
      USE eg_swap_bounds_mod

      IMPLICIT NONE

!
! Description: Code to calculate the star variables for use
!              in the Helmholtz problem and updating of the
!              time level (n+1) fields.
!
! Method: ENDGame formulation version 1.01
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! In

      INTEGER, INTENT(IN) ::  row_length, rows, n_rows, model_levels,         &
                              model_domain, inner_it_val

      REAL                                                                    &
      psi_w_surf(row_length,rows),                                            &
      psi_w_lid (row_length,rows)

! Time step, hydrostatic flag  and off-centre weights
      
      REAL, INTENT(IN) :: Ih

! sources at departure points

      REAL, INTENT(IN) ::    R_u_d(udims_s%i_start:udims_s%i_end,     &
                                   udims_s%j_start:udims_s%j_end,     &
                                   udims_s%k_start:udims_s%k_end)

      REAL, INTENT(IN) ::    R_v_d(vdims_s%i_start:vdims_s%i_end,     &
                                   vdims_s%j_start:vdims_s%j_end,     &
                                   vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(IN) ::      R_w_d(wdims%i_start:wdims%i_end,       &
                                     wdims%j_start:wdims%j_end,       &
                                     wdims%k_start:wdims%k_end)

      REAL, INTENT(IN) ::     R_theta_d(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_rho_d(pdims_s%i_start:tdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)


      REAL, INTENT(IN) ::   R_m_v_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cl_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cf_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::   R_m_r_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_gr_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) :: R_m_cf2_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)
     
! Output 

       REAL :: RHS(pdims_s%i_start:pdims_s%i_end,                     & 
                   pdims_s%j_start:pdims_s%j_end,                     &
                   pdims_s%k_start:pdims_s%k_end)


      REAL, INTENT(OUT)   :: R_u(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,       &
                                 udims_s%k_start:udims_s%k_end)

      REAL, INTENT(OUT)   :: R_v(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,       &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(OUT)   ::   R_w(wdims%i_start:wdims%i_end,         &
                                   wdims%j_start:wdims%j_end,         &
                                   wdims%k_start:wdims%k_end)

      REAL, INTENT(OUT)   ::  R_theta(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT)   :: R_rho(pdims_s%i_start:wdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)


      REAL, INTENT(OUT)   ::  R_p(pdims_s%i_start:wdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(OUT)   :: R_etadot(wdims%i_start:wdims%i_end,      &
                                      wdims%j_start:wdims%j_end,      &
                                      wdims%k_start:wdims%k_end)  
       
      REAL, INTENT(OUT) ::                                            &
                     rho_divu_out(pdims_s%i_start:wdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)


! time level (n+1) fields

      REAL, INTENT(INOUT) ::     u(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,         &
                                 udims_s%k_start:udims_s%k_end)

      REAL, INTENT(INOUT) ::     v(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,         &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(INOUT) ::     w(wdims_s%i_start:wdims_s%i_end,       &
                                 wdims_s%j_start:wdims_s%j_end,         &
                                 wdims_s%k_start:wdims_s%k_end)

      REAL, INTENT(INOUT) ::     exner(pdims_s%i_start:pdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,     &
                                     pdims_s%k_start:pdims_s%k_end+1)

      REAL, INTENT(INOUT) :: exner_star(pdims_s%i_start:pdims_s%i_end,  &
                                      pdims_s%j_start:pdims_s%j_end)
 
      REAL, INTENT(INOUT) ::     rho(pdims_s%i_start:pdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,       &
                                   pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(INOUT) ::     thetav(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,    &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::      m_v(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::     m_cl(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::     m_cf(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::      m_r(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::     m_gr(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::    m_cf2(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)


      REAL, INTENT(INOUT) ::   etadot(wdims_s%i_start:wdims_s%i_end,  &
                                      wdims_s%j_start:wdims_s%j_end,  &
                                      wdims_s%k_start:wdims_s%k_end)

!     Fast Physics source terms                 
      REAL                :: S_u(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,       &
                                 udims_s%k_start:udims_s%k_end)

      REAL                :: S_v(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,       &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL                ::   S_w(wdims%i_start:wdims%i_end,         &
                                   wdims%j_start:wdims%j_end,         &
                                   wdims%k_start:wdims%k_end)

      REAL                :: S_thetav(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL              ::    S_m_v(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_cl(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_cf(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::    S_m_r(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_gr(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::  S_m_cf2(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)
          
! Local variables

      INTEGER :: i, j, k
      LOGICAL :: L_call_from_solver
      REAL    :: T1, T3, T4, D1, a_p

      REAL    :: rdxi1, rdxi2, rdxi3
      REAL    :: u_at_w, v_at_w
      
      REAL    :: R_p_err,R_rho_err,R_etadot_err,R_w_err,R_theta_err,  &
                 R_u_err, R_v_err
                    
      INTEGER :: istat              


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_HELM_RHS_STAR',zhook_in,zhook_handle)

      IF (integrity_test) THEN
! check all input fields

        CALL check_hash_m(                                            &
                  exner,           SIZE(exner),           'pinp1',    & !20
                  u,               SIZE(u),               'u_np1',    &
                  v,               SIZE(v),               'v_np1',    &
                  rho,             SIZE(rho),             'r_np1',    &
                  m_v,             SIZE(m_v),             'mvnp1',    &
                  m_cl,            SIZE(m_cl),            'mclp1',    &
                  m_cf,            SIZE(m_cf),            'mcfp1',    &
                  m_r,             SIZE(m_r),             'mrnp1',    &
                  m_gr,            SIZE(m_gr),            'mgrp1',    & !30
                  m_cf2,           SIZE(m_cf2),           'mcf21',    &
                  thetav,          SIZE(thetav),          'tvnp1',    &
                  w,               SIZE(w),               'w_np1',    &
                  etadot,          SIZE(etadot),          'ednp1',    &
                  S_u,             SIZE(S_u),             's_u__',    &
                  S_v,             SIZE(S_v),             's_v__',    &
                  S_w,             SIZE(S_w),             's_w__',    &
                  S_thetav,        SIZE(S_thetav),        's_the',    &
                  S_m_v,           SIZE(S_m_v),           's_m_v',    &
                  S_m_cl,          SIZE(S_m_cl),          's_mcl',    & !40
                  S_m_cf,          SIZE(S_m_cf),          's_mcf',    &
                  S_m_r,           SIZE(S_m_r),           's_m_r',    &
                  S_m_gr,          SIZE(S_m_gr),          's_mgr',    &
                  S_m_cf2,         SIZE(S_m_cf2),         'smcf2')
        CALL check_hash_m(                                            &
                  exner_ref_pro,   SIZE(exner_ref_pro),   'piref',    &
                  thetav_ref_pro,  SIZE(thetav_ref_pro),  'tvref',    &
                  rho_ref_pro,     SIZE(rho_ref_pro),     'r_ref',    &
                  thetav_ref_eta,  SIZE(thetav_ref_eta),  'tvree',    &
                  rho_ref_eta,     SIZE(rho_ref_eta),     'rreta',    &
                  xi1_p,           SIZE(xi1_p),           'xi1_p',    &
                  xi1_u,           SIZE(xi1_u),           'xi1_u',    &
                  xi2_p,           SIZE(xi2_p),           'xi2_p',    &
                  xi2_v,           SIZE(xi2_v),           'xi2_v',    &
                  phi_at_p,        SIZE(phi_at_p),        'phi_p',    &
                  phi_at_u,        SIZE(phi_at_u),        'phi_u',    &
                  phi_at_v,        SIZE(phi_at_v),        'phi_v',    &
                  phi_at_eta,      SIZE(phi_at_eta),      'phiet',    & !60
                  intw_u2p,        SIZE(intw_u2p),        'iwu2p',    &
                  intw_v2p,        SIZE(intw_v2p),        'iwv2u',    &
                  intw_p2u,        SIZE(intw_p2u),        'iwp2u',    &
                  intw_p2v,        SIZE(intw_p2v),        'iwp2v',    &
                  intw_rho2w,      SIZE(intw_rho2w),      'iwr2w',    &
                  intw_w2rho,      SIZE(intw_w2rho),      'iww2r',    &
                  deta_xi3,        SIZE(deta_xi3),        'dexi3',    &
                  deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t',    &
                  deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u',    &
                  deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v',    & !70
                  dxi1_xi3,        SIZE(dxi1_xi3),        'dx1x3',    &
                  dxi2_xi3,        SIZE(dxi2_xi3),        'dx2x3',    &
                  h1_p,            SIZE(h1_p),            'h1_p_',    &
                  h1_xi1_u,        SIZE(h1_xi1_u),        'h1x1u',    &
                  h1_xi2_v,        SIZE(h1_xi2_v),        'h1x2v',    &
                  h1_p_eta,        SIZE(h1_p_eta),        'h1pet',    &
                  h2_p,            SIZE(h2_p),            'h2_p_',    &
                  h2_xi1_u,        SIZE(h2_xi1_u),        'h2x1u',    &
                  h2_xi2_v,        SIZE(h2_xi2_v),        'h2xiv',    &
                  h2_p_eta,        SIZE(h2_p_eta),        'h2pet',    & !80
                  h3_p,            SIZE(h3_p),            'h3_p_',    &
                  h3_xi1_u,        SIZE(h3_xi1_u),        'h3x1u',    &
                  h3_xi2_v,        SIZE(h3_xi2_v),        'h3x2v',    &
                  h3_p_eta,        SIZE(h3_p_eta),        'h3pet',    &
                  f1_star,         SIZE(f1_star),         'f1sta',    &
                  f2_star,         SIZE(f2_star),         'f2sta',    &
                  f3_star,         SIZE(f3_star),         'f3sta',    &
                  f1_comp,         SIZE(f1_comp),         'f1cmp',    &
                  f2_comp,         SIZE(f2_comp),         'f2cmp',    &
                  f3_comp,         SIZE(f3_comp),         'f3cmp',    & !90
                  g_theta,         SIZE(g_theta),         'gthet',    &
                  g_rho,           SIZE(g_rho),           'g_rho',    &
                  Csxi1_p,         SIZE(Csxi1_p),         'cxi1p',    &
                  Csxi1_u,         SIZE(Csxi1_u),         'cxi1u',    &
                  Csxi2_p,         SIZE(Csxi2_p),         'cxi2p',    &
                  Csxi2_v,         SIZE(Csxi2_v),         'cxi2v',    &
                  Snxi1_p,         SIZE(Snxi1_p),         'sxi1p',    &
                  Snxi1_u,         SIZE(Snxi1_u),         'sxi1u',    &
                  Snxi2_p,         SIZE(Snxi2_p),         'sxi2p',    &
                  Snxi2_v,         SIZE(Snxi2_v),         'sxi2v')      !100

      END IF
!
! Add the fast physics source terms
!

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,udims,R_u, S_u,R_u_d,u,      &
!$OMP&            vdims,R_v,S_v,R_v_d,v)
      DO k = 1, model_levels

        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            R_u(i,j,k) = S_u(i,j,k) + R_u_d(i,j,k) - 2.0*u(i,j,k)
          END DO
        END DO

        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            R_v(i,j,k) = S_v(i,j,k) + R_v_d(i,j,k) - 2.0*v(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO PRIVATE(i,j) SHARED( model_levels,pdims,R_w,S_w,    &
!$OMP&            R_w_d,Ih,mu_w,w)
      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            R_w(i,j,k) = S_w(i,j,k) + R_w_d(i,j,k)                    &
                                - (2.0*Ih + mu_w(i,j,k))*w(i,j,k)
          END DO
        END DO
      END DO  
!$OMP END PARALLEL DO   
      
      IF( inner_it_val == 1 ) THEN
!$OMP PARALLEL DO PRIVATE(i,j) SHARED( model_levels, pdims,R_theta,   &
!$OMP&            S_thetav,R_theta_d,thetav)
        DO k = 0, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              R_theta (i,j,k) = S_thetav(i,j,k) + R_theta_d(i,j,k)    &
                                 - thetav(i,j,k)                 
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO  
      ELSE
!$OMP PARALLEL DO SHARED(R_theta)
        DO k = 0, model_levels
          R_theta (:,:,k) = 0.0
        END DO
!$OMP END PARALLEL DO
      END IF
                      
      IF( inner_it_val == 1 ) THEN

        CALL update_moisture(m_v, m_cl, m_cf, m_r,m_gr, m_cf2,        &
                             R_m_v_d, R_m_cl_d, R_m_cf_d,             &
                             R_m_r_d, R_m_gr_d, R_m_cf2_d,            &
                             S_m_v, S_m_cl, S_m_cf,                   &
                             S_m_r, S_m_gr, S_m_cf2)
      END IF


! Adjust arrival point values to account for extra terms
 
!------------------------------------------------------------------------------
! Compute Psi_w_surf to be used in eg_sisl_init to compute bottom R_w
!------------------------------------------------------------------------------
      a_p = 1.0/(alpha_w*timestep)
      k = 0
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          psi_w_surf(i,j) = ( Ih*(                                    &
                                (h3_p_eta(i,j,0)/h1_p_eta(i,j,0))*    &
                            dxi1_xi3(i,j,0)*(intw_u2p(i,1)*u(i-1,j,1) &
                                       +intw_u2p(i,2)*u(i,j,1))   +   &
                       (h3_p_eta(i,j,0)/h2_p_eta(i,j,0))*             &
                       dxi2_xi3(i,j,0)*(intw_v2p(j,1)*v(i,j-1,1)      &
                                        +intw_v2p(j,2)*v(i,j,1))) -   &
                                      R_w_d(i,j,0) )*a_p

          psi_w_lid(i,j) = -R_w_d(i,j,model_levels)*a_p
        END DO
      END DO

      L_call_from_solver = .TRUE.

      CALL eg_sisl_init(                                              &
               row_length, rows, n_rows, model_levels,                &
               .TRUE.,L_call_from_solver,.FALSE., Ih,g_theta,         &
               u, v, w,  thetav, rho,                                 &
               m_v, m_cl, m_cf, m_r, m_gr, m_cf2,                     &
               exner, exner_star, R_u, R_v, R_w, R_theta, R_rho,      &
               S_m_v, S_m_cl, S_m_cf, S_m_r, S_m_gr, S_m_cf2, etadot, &
               psi_w_surf, psi_w_lid)
               
      
! Compute R_etadot 
!$OMP PARALLEL DO PRIVATE(i,j,u_at_w,v_at_w) SHARED(model_levels,     &
!$OMP&            pdims,intw_u2p,u,intw_rho2w,intw_v2p,               &
!$OMP&            v,R_etadot,HM_etadot,etadot,w,dxi1_xi3,h3_p_eta,    &
!$OMP&            h1_p_eta,dxi2_xi3,h2_p_eta)
      DO k = 1, model_levels-1
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end

            u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +   &
                                       intw_u2p(i,2)*u(i,j,k+1) ) +   &
                     intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +     &
                                       intw_u2p(i,2)*u(i,j,k) )

            v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +   &
                                       intw_v2p(j,2)*v(i,j,k+1) ) +   &
                     intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +     &
                                       intw_v2p(j,2)*v(i,j,k) )

            R_etadot(i,j,k) = - HM_etadot(i,j,k)*etadot(i,j,k)        &
                              +  w(i,j,k)                             &
                              -  u_at_w*dxi1_xi3(i,j,k)               &
                                *h3_p_eta(i,j,k)/h1_p_eta(i,j,k)      &
                              -  v_at_w*dxi2_xi3(i,j,k)               &
                                 *h3_p_eta(i,j,k)/h2_p_eta(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO      
      k = 0
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          u_at_w =  intw_u2p(i,1)*u(i-1,j,k+1) +                      &
                    intw_u2p(i,2)*u(i,j,k+1) 

          v_at_w = intw_v2p(j,1)*v(i,j-1,k+1) +                       &
                   intw_v2p(j,2)*v(i,j,k+1)

          R_etadot(i,j,k) =   w(i,j,k)                                &
                             - u_at_w*dxi1_xi3(i,j,k)                 &
                              *h3_p_eta(i,j,k)/h1_p_eta(i,j,k)        &
                             - v_at_w*dxi2_xi3(i,j,k)                 &
                              *h3_p_eta(i,j,k)/h2_p_eta(i,j,k)
        END DO
      END DO
      k = model_levels
      R_etadot(:,:,k) = 0.0
     
! Compute R_pi
!$OMP PARALLEL DO PRIVATE(i,j,T1) SHARED(model_levels,pdims,          &
!$OMP&            intw_w2rho,thetav,R_p,rho,exner,Hm_pp)
      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end  
            T1 = intw_w2rho(k,1)*thetav(i,j,k)                        &
               + intw_w2rho(k,2)*thetav(i,j,k-1) 
            R_p(i,j,k) = - rho(i,j,k) + p_zero/(R*T1)                 &
                          * exner(i,j,k)**Hm_pp
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO      
     
! Compute R_rho 
      IF( inner_it_val == 1 ) THEN
!$OMP PARALLEL DO PRIVATE(i,j,T1) SHARED(model_levels,pdims,          &
!$OMP&            rho_divu_out,rho,R_rho,rho_ref_pro)
        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
           ! Compute (rho-rho_ref)*div(u) (R_rho contains rho - rho*div(u)
              rho_divu_out(i,j,k) = (rho(i,j,k) - R_rho(i,j,k))       &
                              * (1.0 - rho_ref_pro(i,j,k)/rho(i,j,k))
           ! first inner loop use R_rho as usual   
              R_rho(i,j,k) = R_rho(i,j,k) + R_rho_d(i,j,k)            &
                             - 2.0*rho(i,j,k)
            END DO
          END DO
        END DO 
!$OMP END PARALLEL DO      
      ELSE
        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end             
              R_rho(i,j,k) = (rho(i,j,k) - R_rho(i,j,k))              &
                           * (1.0 - rho_ref_pro(i,j,k)/rho(i,j,k))    &
                           - rho_divu_out(i,j,k)
            END DO
          END DO
        END DO
      END IF
      
! Output initial errors
      IF( PrintStatus >= PrStatus_Diag) THEN
        R_p_err      = -1.0
        R_rho_err    = -1.0
        R_etadot_err = -1.0
        R_w_err      = -1.0
        R_theta_err  = -1.0
        R_u_err      = -1.0
        R_v_err      = -1.0

        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              R_p_err = MAX(R_p_err,ABS(R_p(i,j,k)/exner(i,j,k)))
              R_rho_err = MAX(R_rho_err,ABS(R_rho(i,j,k)/rho(i,j,k)))
              R_etadot_err = MAX(R_etadot_err,ABS(R_etadot(i,j,k)))
              R_w_err = MAX(R_w_err,ABS(R_w(i,j,k)))
              R_theta_err = MAX(R_theta_err,ABS(R_theta(i,j,k)        &
                                                /thetav(i,j,k)))
            END DO
          END DO
          DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
              R_u_err = MAX(R_u_err,ABS(R_u(i,j,k)))         
            END DO
          END DO

          DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
              R_v_err = MAX(R_v_err,ABS(R_v(i,j,k)))         
            END DO
          END DO
        END DO 

        CALL gc_rmax(1,n_proc,istat,R_p_err)
        CALL gc_rmax(1,n_proc,istat,R_rho_err)
        CALL gc_rmax(1,n_proc,istat,R_etadot_err)
        CALL gc_rmax(1,n_proc,istat,R_w_err)
        CALL gc_rmax(1,n_proc,istat,R_theta_err)
        CALL gc_rmax(1,n_proc,istat,R_u_err)
        CALL gc_rmax(1,n_proc,istat,R_v_err)

        IF( mype == 0 ) THEN
          WRITE(6,fmt='(A)')                                          &
                '_________________________________________________'     
          WRITE(6,fmt='(A,E25.5)') '    Error in R_u      ',          &
                                   R_u_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_v      ',          &
                                   R_v_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_w      ',          &
                                   R_w_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_theta  ',          &
                                   R_theta_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_etadot ',          &
                                   R_etadot_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_rho    ',          &
                                   R_rho_err
          WRITE(6,fmt='(A,E25.5)') '    Error in R_p      ',          &
                                   R_p_err
          WRITE(6,fmt='(A)')                                          &
                '_________________________________________________'
        END IF
      END IF
!==== Elimination procedure =============================================

!$OMP PARALLEL DO PRIVATE(i,j,T1) SHARED(model_levels,pdims,R_w,Ih,   &
!$OMP&            mu_w,R_etadot,HM_b,HM_w,R_theta,intw_w2rho,         &
!$OMP&            thetav_ref_pro,RHS,rho_ref_pro,R_p,R_rho)
      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end  
!           Eliminate w' from w equation            
            R_w(i,j,k) = R_w(i,j,k) +                                 &
                         (Ih+mu_w(i,j,k))                             &
                           *R_etadot(i,j,k)                           &
                          - HM_b(i,j,k)*HM_w(i,j,k)*R_theta(i,j,k)
            
            T1 = intw_w2rho(k,1)*R_theta(i,j,k) /thetav_ref_pro(i,j,k)&
               + intw_w2rho(k,2)*R_theta(i,j,k-1)                     &
                    /thetav_ref_pro(i,j,k-1)

!           Eliminate theta' & rho' from Equation of state                
            RHS(i,j,k)   = rho_ref_pro(i,j,k)*(R_p(i,j,k) - T1)       &
                            - R_rho(i,j,k)  
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO    


! Fix R_v at the poles

      IF ( model_domain == mt_global ) THEN
        IF ( at_extremity(PSouth) ) THEN
          CALL eg_v_at_poles(R_u,R_v, 1.0, udims%j_start,             &
                             vdims%j_start,                           &
                             udims_s,vdims_s)
        END IF

        IF ( at_extremity(PNorth) ) THEN
          CALL eg_v_at_poles(R_u,R_v,-1.0, udims%j_end, vdims%j_end,  &
                             udims_s,vdims_s)
        END IF
      END IF

      CALL eg_swap_bounds(R_u,udims_s,fld_type_u,.true.)

      CALL eg_swap_bounds(R_v,vdims_s,fld_type_v,.true.)

! Now build RHS
!$OMP PARALLEL DO PRIVATE(i,j,T3,T4,D1,rdxi1,rdxi2,rdxi3)
      DO k = 1, model_levels
        rdxi3 = 1.0/(eta_theta_levels(k) - eta_theta_levels(k-1))
        DO j = pdims%j_start, pdims%j_end
          rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
          DO i = pdims%i_start, pdims%i_end
            rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))
                                      
            RHS(i,j,k) =  RHS(i,j,k)                                  &
                         +HM_vol(i,j,k)*(                             &
                             rdxi1*( Hm_rhox(i,j,k)  *R_u(i,j,k)      &
                                    -Hm_rhox(i-1,j,k)*R_u(i-1,j,k) )  &
                            +rdxi2*( Hm_rhoy(i,j,k)  *R_v(i,j,k)      &
                                    -Hm_rhoy(i,j-1,k)*R_v(i,j-1,k) )  &
                                            )

!           Add in D1 terms - T3 = 0 at top, T4 = 0 at bottom
            T3 = R_w(i,j,k)  *Hlm_Ck(i,j,k)  
            T4 = R_w(i,j,k-1)*Hlm_Ck(i,j,k-1)
               
            D1 = HM_vol(i,j,k)*(Hm_rhoz(i,j,k)*T3                     &
                                 - Hm_rhoz(i,j,k-1)*T4)*rdxi3
               
            D1 = D1 + rho_ref_pro(i,j,k)*(                            &
                                  intw_w2rho(k,1)*Hm_theta(i,j,k)     &
                                    /thetav_ref_pro(i,j,k)*T3         &
                                + intw_w2rho(k,2)*Hm_theta(i,j,k-1)   &
                                     /thetav_ref_pro(i,j,k-1)*T4      &
                                            )
                                   
            RHS(i,j,k) = (RHS(i,j,k) + D1)*Hlm_Lp(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO    

      IF (lhook) CALL dr_hook('EG_HELM_RHS_STAR',zhook_out,zhook_handle)

      END SUBROUTINE eg_helm_rhs_inc
      END MODULE eg_helm_rhs_inc_mod
