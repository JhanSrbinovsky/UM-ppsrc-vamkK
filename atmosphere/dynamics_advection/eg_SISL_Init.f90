! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_init_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sisl_init(                                              &
         row_length,rows,n_rows,model_levels, l_inc_solver,           &
         l_call_from_solver, l_call_from_f1sp, ih, g_theta,u, v, w,   &
         thetav, rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2,              &
         exner,exner_star, r_u, r_v, r_w, r_theta, r_rho, r_m_v,      &
         r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2, etadot,              &
         psi_w_surf,psi_w_lid )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod,         ONLY : model_domain
USE level_heights_mod,     ONLY : xi3_at_u => r_at_u,                 &
                                  xi3_at_v => r_at_v,                 &
                                  xi3_at_theta=>r_theta_levels,       &
                                  xi3_at_rho=>r_rho_levels,           &
                                  eta_theta_levels, eta_rho_levels
USE timestep_mod,          ONLY : timestep
USE eg_helmholtz_mod,      ONLY : hm_vol
USE eg_sisl_init_uvw_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE metric_terms_mod
USE coriolis_mod, ONLY : f1_star,f2_star,f3_star
USE um_parvars,   ONLY : halo_i, halo_j, offx, offy
USE eg_alpha_mod, ONLY : alpha_rho
USE eg_parameters_mod, ONLY : l_rho_av_zz

IMPLICIT NONE
!
! Description: computes time level n arrival point quantities:
!              Ru, Rv, Rw, Rtheta, Rrho, Rm_v, Rm_cl, Rm_cf
!  
!
! Method: ENDGame formulation version 1.01,
!         section 11 (solution procedure), paragraph 1.
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

! Loop index bounds for arrays defined on p, u, v points respectively

LOGICAL, INTENT(IN) :: l_inc_solver,                                  &
                       l_call_from_solver,l_call_from_f1sp

! SI time weights   & hydrosctatic switch

REAL, INTENT(IN) :: ih


! Gravity arrays

REAL, INTENT(IN) ::                                                   &
  g_theta (1-offx:row_length+offx, 1-offy:rows+offy,                  &
           0:model_levels)

REAL, INTENT(IN) ::                                                   &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),           &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),         &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),          &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),      &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy),                &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels ),         &
  m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels )

REAL, INTENT(IN) ::                                                   &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(INOUT) ::                                                &
  r_u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels)          &
, r_v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels)        &
, r_w(row_length,rows,0:model_levels)                                 &
, r_theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)     &
, r_rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)         &
, r_m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)       &
, r_m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)      &
, r_m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)      &
, r_m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)       &
, r_m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)      &
, r_m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)


! Local variables

INTEGER :: i,j,k

REAL :: beta_dt,  d_xi1_term, d_xi2_term, deta_term, del_rho

REAL :: rdxi1, rdxi2, rdxi3

REAL :: rho_av, rho_ref_term, rho_switch

REAL                                                                    &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SISL_INIT',zhook_in,zhook_handle)

IF (integrity_test) THEN
  IF(l_call_from_f1sp.or.l_call_from_solver) THEN
    CALL check_hash_m(                                                &
                  exner,           SIZE(exner),           'pinp1',    & 
                  u,               SIZE(u),               'u_np1',    &
                  v,               SIZE(v),               'v_np1')
    IF(l_call_from_f1sp) THEN
      CALL check_hash_m(                                              &
                  etadot,          SIZE(etadot),          'ednp1')
    ENDIF
  ELSE
    CALL check_hash_m(                                                &
                  exner,           SIZE(exner),           'pi___',    & 
                  u,               SIZE(u),               'u____',    &
                  v,               SIZE(v),               'v____',    &
                  etadot,          SIZE(etadot),          'ed___')
  ENDIF

  CALL check_hash_m(                                                  &
                  hm_vol,          SIZE(hm_vol),          'hmvol',    &
                  exner_ref_pro,   SIZE(exner_ref_pro),   'piref',    &
                  xi1_u,           SIZE(xi1_u),           'xi1_u',    &
                  xi2_v,           SIZE(xi2_v),           'xi2_v',    &
                  intw_u2p,        SIZE(intw_u2p),        'iwu2p',    &
                  intw_v2p,        SIZE(intw_v2p),        'iwv2u',    &
                  deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t',    &
                  deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u',    &
                  deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v',    &
                  eta_theta_levels,SIZE(eta_theta_levels),'etatl',    &
                  h1_xi2_v,        SIZE(h1_xi2_v),        'h1x2v',    &
                  h1_p,            SIZE(h1_p),            'h1_p_',    &
                  h1_p_eta,        SIZE(h1_p_eta),        'h1pet',    &
                  h2_xi1_u,        SIZE(h2_xi1_u),        'h2x1u',    &
                  h2_p,            SIZE(h2_p),            'h2_p_',    &
                  h2_p_eta,        SIZE(h2_p_eta),        'h2pet',    &
                  h3_xi1_u,        SIZE(h3_xi1_u),        'h3x1u',    &
                  h3_xi2_v,        SIZE(h3_xi2_v),        'h3x2v',    &
                  h3_p_eta,        SIZE(h3_p_eta),        'h3pet')
END IF

CALL eg_sisl_init_uvw(                                                &
         row_length, rows, n_rows, model_levels, ih, g_theta,         &
         u, v, w, thetav, rho, m_v, m_cl, m_cf, m_r,                  &
         m_gr, m_cf2, exner, exner_star, r_u, r_v, r_w ,              &
         l_call_from_solver, l_call_from_f1sp,psi_w_surf,psi_w_lid)

not_from_solver: IF(.NOT. l_call_from_solver ) THEN

!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED (model_levels,pdims,r_m_v,    &
!$OMP&            r_m_cl,r_m_cf,r_m_r,r_m_gr,r_m_cf2,m_v,m_cl,        &
!$OMP&            m_r,m_cf,m_gr,m_cf2) DEFAULT(NONE) SCHEDULE(STATIC)
DO k=0, model_levels
   DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
         r_m_v(i,j,k)   = r_m_v(i,j,k)   + m_v(i,j,k)
         r_m_cl(i,j,k)  = r_m_cl(i,j,k)  + m_cl(i,j,k)
         r_m_cf(i,j,k)  = r_m_cf(i,j,k)  + m_cf(i,j,k)
         r_m_r(i,j,k)   = r_m_r(i,j,k)   + m_r(i,j,k)
         r_m_gr(i,j,k)  = r_m_gr(i,j,k)  + m_gr(i,j,k)
         r_m_cf2(i,j,k) = r_m_cf2(i,j,k) + m_cf2(i,j,k)
      END DO
   END DO
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! Compute R_theta term
!-----------------------------------------------------------------------

    DO k=0, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          r_theta(i,j,k) = r_theta(i,j,k) + thetav(i,j,k)
        END DO
      END DO
    END DO
END IF not_from_solver


IF( (.NOT. l_call_from_f1sp) .AND.                                        &
    (l_inc_solver .OR. (.NOT.l_call_from_solver)) ) THEN
!-----------------------------------------------------------------------
! Compute R_rho term
!-----------------------------------------------------------------------
! add -beta*Dt*div(u) term to R_rho - see eqn (5.15) in formulation doc

    del_rho = 1.0
    
    rho_switch = 1.0
    IF ( l_rho_av_zz ) rho_switch = 0.0

    IF(l_call_from_solver) THEN
      beta_dt  = del_rho*alpha_rho*timestep
    ELSE
      beta_dt  = del_rho*(1.0 - alpha_rho)*timestep
    END IF

    k = 1
      rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
      DO j=pdims%j_start, pdims%j_end
        rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
        DO i=pdims%i_start, pdims%i_end
          rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u(i,j,k) -                   &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u(i-1,j,k) )*rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v(i,j,k) -                   &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v(i,j-1,k) )*rdxi2

          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                        h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*        &
                        etadot(i,j,k) )*rdxi3

          r_rho(i,j,k) = rho(i,j,k)*( 1.0 -                           &
                         beta_dt*                                     &
                           (d_xi1_term + d_xi2_term + deta_term)      &
                          /(h1_p(i,j,k)*h2_p(i,j,k)*                  &
                                     h3_p(i,j,k)*deta_xi3(i,j,k) ))
        END DO
      END DO

!$OMP PARALLEL DO PRIVATE(i,j,k,rdxi1,rdxi2,rdxi3,d_xi1_term,d_xi2_term,&
!$OMP&            deta_term,rho_av,rho_ref_term) SHARED (model_levels,&
!$OMP&            eta_theta_levels,pdims,xi2_v,xi1_u,                 &
!$OMP&            deta_xi3_u,u,h2_xi1_u,h3_xi1_u,                     &
!$OMP&            h3_xi2_v,v,h1_xi2_v, intw_rho2w,intw_w2rho,         &
!$OMP&            deta_xi3_v,h1_p_eta,h2_p_eta,h3_p_eta,              &
!$OMP&            deta_xi3_theta,etadot,r_rho,rho,beta_dt,h1_p,h2_p,  &
!$OMP&        h3_p,deta_xi3,rho_switch) DEFAULT(NONE) SCHEDULE(STATIC)
    DO k=2, model_levels-1
      rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
      DO j=pdims%j_start, pdims%j_end
        rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
        DO i=pdims%i_start, pdims%i_end
          rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u(i,j,k) -                   &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u(i-1,j,k) )*rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v(i,j,k) -                   &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v(i,j-1,k) )*rdxi2

          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                    h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*            &
                    etadot(i,j,k) -                                   &
                       h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*           &
                       h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                       etadot(i,j,k-1) )*rdxi3

          rho_av = intw_w2rho(k,1)*(intw_rho2w(k,1)*rho(i,j,k+1)      &
                                   +intw_rho2w(k,2)*rho(i,j,k)    )   &
                 + intw_w2rho(k,2)*(intw_rho2w(k-1,1)*rho(i,j,k)      &
                                   +intw_rho2w(k-1,2)*rho(i,j,k-1)) 

          rho_ref_term = rho_switch + (1.0 - rho_switch)*             &
                         rho_av/rho(i,j,k)

          r_rho(i,j,k) = rho(i,j,k)*( 1.0  -                          &
                         beta_dt*rho_ref_term*                        &
                          (d_xi1_term + d_xi2_term + deta_term)       &
                         /( h1_p(i,j,k)*h2_p(i,j,k)*                  &
                                     h3_p(i,j,k)*deta_xi3(i,j,k) ))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    k = model_levels
    rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
      DO j=pdims%j_start, pdims%j_end
        rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
        DO i=pdims%i_start, pdims%i_end
          rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u(i,j,k) -                   &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u(i-1,j,k) )*rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v(i,j,k) -                   &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v(i,j-1,k) )*rdxi2

          deta_term = ( - h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*        &
                      h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*      &
                      etadot(i,j,k-1) )*rdxi3

          r_rho(i,j,k) = rho(i,j,k)*( 1.0 -                           &
                         beta_dt*                                     &
                          (d_xi1_term + d_xi2_term + deta_term)       &
                         /( h1_p(i,j,k)*h2_p(i,j,k)*                  &
                                      h3_p(i,j,k)*deta_xi3(i,j,k) ))
        END DO
      END DO
ENDIF

IF (integrity_test.AND..NOT.l_call_from_solver                        &
                  .AND..NOT.l_call_from_f1sp) THEN
  CALL update_hash_m(                                                 &
                  R_theta,         SIZE(R_theta),         'R_t__',    &
                  R_rho,           SIZE(R_rho),           'R_r__',    &
                  R_m_v,           SIZE(R_m_v),           'R_mv_',    &
                  R_m_cl,          SIZE(R_m_cl),          'Rmcl_',    &
                  R_m_cf,          SIZE(R_m_cf),          'Rmcf_',    &
                  R_m_r,           SIZE(R_m_r),           'Rmr__',    &
                  R_m_gr,          SIZE(R_m_gr),          'Rmgr_',    &
                  R_m_cf2,         SIZE(R_m_cf2),         'Rmc2_')
  CALL update_hash_m(                                                 &
                  R_u,             SIZE(R_u),             'R_u__',    & 
                  R_v,             SIZE(R_v),             'R_v__',    &
                  R_w,             SIZE(R_w),             'R_w__')
END IF

IF (lhook) CALL dr_hook('EG_SISL_INIT',zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_init
END MODULE eg_sisl_init_mod
