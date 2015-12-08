! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_f1sp_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_f1sp(theta_star, m_star,mcl_star,mcf_star               &
                   ,mcf2_star,mrain_star,mgraup_star                  &
                   ,u,v,w,un,vn,wn,exner,exner_star                   &
                   ,rho, thetav, etadot                               &
                   ,m_v, m_cl, m_cf, m_r, m_gr, m_cf2,ih              &
                   ,l_slice,row_length, rows, n_rows, model_levels    &
                   ,L_mcr_qcf2,L_mcr_qrain                            &
                   ,L_mcr_qgraup,psi_w_surf, psi_w_lid)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: recip_epsilon
USE timestep_mod,      ONLY : timestep
USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,       &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v
USE proc_info_mod,     ONLY : model_domain,                           &
                              global_row_length,at_extremity,         &
                              gc_proc_row_group, me, n_proc
USE eg_helmholtz_mod
USE eg_sisl_init_mod
USE eg_v_at_poles_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE atm_fields_mod, ONLY : 
USE PrintStatus_mod
USE UM_ParParams
USE metric_terms_mod
USE fields_rhs_mod, ONLY : r_u_p2,r_v_p2,r_w_p2,r_u_d, r_v_d, r_w_d,  &
                           r_theta_d  ,r_m_v_d, r_m_cl_d, r_m_cf_d,   &
                           r_m_cf2_d,r_m_r_d,r_m_gr_d

USE eg_alpha_mod
USE coriolis_mod
USE eg_parameters_mod, ONLY : pole_consts 
USE gravity_mod

USE um_parvars, ONLY : halo_i, halo_j, offx, offy 
USE eg_swap_bounds_mod


IMPLICIT NONE
!
! Description: Code to calculate the predictor fields for calling
!              atmos_physics2
!  
!              In ND atm_step, R_u, R_v_p2 and R_w are always increments.
!
!              theta_star and q_star are full fields EXCEPT 
!              immediately after Atmos_physics1 which returns
!              increments. 
!  
!              The same should apply for ENDGame.
!
! Method:
!         See ENDGame Formulation version 3.03 Section 7.10
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




LOGICAL, INTENT(IN) ::  l_slice,L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup
INTEGER, INTENT(IN) ::  row_length, rows, n_rows, model_levels

! Time step, hydrostatic flag  and off-centre weights

REAL, INTENT(IN) :: ih

REAL  del_rho

REAL                                                                  &
  psi_w_surf(row_length,rows),                                        &
  psi_w_lid (row_length,rows)


REAL, INTENT(OUT) ::                                                  &
 theta_star(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL, INTENT(INOUT) , DIMENSION (1-offx:row_length+offx,              &
                                 1-offy:rows+offy, 0:model_levels)::  &
  m_star, mcl_star, mcf_star,mcf2_star,mrain_star,mgraup_star

REAL, INTENT(IN) ::                                                   &
  u(-offx:row_length+offx-1,1-offy:rows+offy,model_levels)            &
, v(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels)          &
, w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)           &
, un(-offx:row_length+offx-1,1-offy:rows+offy,model_levels)           &
, vn(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels)         &
, wn(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)          &
, exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels)         &
, rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)           &
, thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)


REAL, DIMENSION(1-offx:row_length+offx,1-offy:rows+offy,              &
                                            0:model_levels),          &
      INTENT(IN) :: m_v, m_cl, m_cf, m_r, m_gr, m_cf2

REAL, INTENT(INOUT) ::                                                &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL                                                                  &
  r_theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)     &
, r_rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)         &
, exner_star(1-offx:row_length+offx, 1-offy:rows+offy) 

REAL                                                                  &
  r_m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)       &
, r_m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)      &
, r_m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels) 

REAL, DIMENSION(1-offx:row_length+offx,1-offy:rows+offy,              &
                                            0:model_levels)::         &
     r_m_v, r_m_cl, r_m_cf


INTEGER :: i, j, k
LOGICAL :: l_call_from_solver
REAL    :: rdeta

REAL max_r_u   ! diagnostic, min/max of r_u/r_v_p2
REAL max_r_v
REAL max_r_w
REAL max_r_theta
REAL max_m_star
REAL max_mcl_star
REAL max_mcf_star
REAL max_mcf2_star
REAL max_mrain_star
REAL max_mgraup_star

REAL min_r_u
REAL min_r_v
REAL min_r_w
REAL min_r_theta
REAL min_m_star
REAL min_mcl_star
REAL min_mcf_star
REAL min_mcf2_star
REAL min_mrain_star
REAL min_mgraup_star

REAL max_u   ! diagnostic, min/max of r_u/r_v
REAL max_v
REAL max_w
REAL max_thetav
REAL min_u
REAL min_v
REAL min_w
REAL min_thetav

INTEGER ierr


IF (lhook) CALL dr_hook('EG_F1SP',zhook_in,zhook_handle)

del_rho = 0.0

IF ( .NOT. l_slice ) del_rho = 1.0

r_u_p2       = 0.0       ! Initialize to zero
r_v_p2       = 0.0
r_w_p2       = 0.0

! The following arrays do not need initialization they can be deleted
r_rho     = 0.0

r_theta   = 0.0
r_m_v     = 0.0
r_m_cl    = 0.0
r_m_cf    = 0.0
r_m_r     = 0.0
r_m_gr    = 0.0
r_m_cf2   = 0.0


! In ND atm_step, R_u_p2, R_v and R_w are always increments.
! theta_star and q_star are full fields EXCEPT immediately after
! Atmos_physics1 which returns increments. These  are input into
! SL_thermo. On exit from sl_thermo, theta_star and q_star are full
! fields.

! Adjust arrival point values to account for extra terms

!------------------------------------------------------------------------------
! Compute Psi_w_surf to be used in eg_sisl_init to compute bottom R_w
!------------------------------------------------------------------------------

k = 0
DO j=pdims%j_start, pdims%j_end
   DO i=pdims%i_start, pdims%i_end

      psi_w_surf(i,j) = ( ih*(                                        &
                          (h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*          &
                          dxi1_xi3(i,j,k)*(intw_u2p(i,1)*u(i-1,j,k+1) &
                                           +intw_u2p(i,2)*u(i,j,k+1))+&
                          (h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*          &
                          dxi2_xi3(i,j,k)*(intw_v2p(j,1)*v(i,j-1,k+1) &
                                          +intw_v2p(j,2)*v(i,j,k+1)))-&
                                r_w_d(i,j,k) )/(alpha_w*timestep)

      psi_w_lid(i,j) = -R_w_d(i,j,model_levels)/(alpha_w*timestep)

   END DO
END DO

  l_call_from_solver = .TRUE.

! here we want alphas. In SISL init betas are used. In order to get alphas in SISL init we need to call it with betas!

CALL eg_sisl_init(                                                    &
         row_length,rows,n_rows,model_levels,.FALSE.,                 &
         l_call_from_solver,.TRUE., ih,g_theta, u, v, w,              &
         thetav, rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2, exner,       &
         exner_star,r_u_p2,r_v_p2,r_w_p2,r_theta,r_rho,r_m_v, r_m_cl, &
         r_m_cf, r_m_r, r_m_gr, r_m_cf2, etadot,                      &
         psi_w_surf,psi_w_lid )

DO k=1, model_levels
   DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end

!        Eq. 7.63: u^F1SP = R^n_u+alpha_u \Delta t Psi^(n+1)_u
!        Eq. 7.1 :                alpha_u \Delta t Psi^(n+1)_u = -R^n_u + u + alpha_u \Delta t S^u
! however, note that SISL Init computes the RHS term as in Eq. 7.4, i.e
!        u + beta_u  \Delta t Psi^(n+1)_u + beta_u \Delta t S^u = R^n_u
! therefore, with beta==alpha
!            alpha_u \Delta t Psi^(n+1)_u                       = R^n_u - u -beta_u \Delta t S^u
!        

! note further that here S==0, as it is the yet unknown fast physics source term.
!
! also note that, when comparing with the documentation:
!            \alpha^S_u == 1 and R_u_p2 as returned from
!
! and physics already contains the Delta t


!        this is u^F1SP:
         r_u_p2(i,j,k) = r_u_d(i,j,k) + r_u_p2(i,j,k) - u(i,j,k)
!        since here we want the tendency, we have to subtract un:
         r_u_p2(i,j,k) = r_u_p2(i,j,k) - un(i,j,k)

      END DO
   END DO

   DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end

         r_v_p2(i,j,k) = r_v_p2(i,j,k) - v(i,j,k) - vn(i,j,k) +       &
                         r_v_d(i,j,k)

      END DO
   END DO
END DO

CALL eg_swap_bounds(r_u_p2,udims_s,fld_type_u,.TRUE.)
CALL eg_swap_bounds(r_v_p2,vdims_s,fld_type_v,.TRUE.)
IF( model_domain == mt_LAM ) THEN
!   DEPENDS ON: fill_external_halos
  CALL fill_external_halos(r_u_p2,row_length,rows,model_levels,offx,offy)
!   DEPENDS ON: fill_external_halos
  CALL fill_external_halos(r_v_p2,row_length,n_rows,model_levels,offx,offy)  
END IF

DO k=1, model_levels-1
   DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end

        r_w_p2(i,j,k) = r_w_p2(i,j,k) - ih*(w(i,j,k)+wn(i,j,k))       &
                       +r_w_d(i,j,k)

      END DO
   END DO
END DO

DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end

    r_w_p2(i,j,0)            = w(i,j,0)-wn(i,j,0)
    r_w_p2(i,j,model_levels) = 0.

  END DO
END DO

! Fix R_v at the poles
IF( model_domain == mt_global ) THEN
   IF( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(r_u_p2,r_v_p2, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)

   END IF

   IF( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(r_u_p2,r_v_p2,-1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)
   END IF

END IF


DO k=1, model_levels-1
  rdeta = alpha_theta*timestep/(eta_rho_levels(k+1) - eta_rho_levels(k))
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end

!     theta predictor:

!     Eq 7.66 of formulation. Note that version 2 does not use the 
!                             full u.grad(theta^star0), only the 
!                             vertical component!

      theta_star(i,j,k) = R_theta_d(i,j,k) + thetav_ref_pro(i,j,k)

      theta_star(i,j,k) = theta_star(i,j,k) - rdeta*etadot(i,j,k)     &
                         *(thetav_ref_eta(i,j,k+1)-thetav_ref_eta(i,j,k))
 
      m_star  (i,j,k) = r_m_v_d (i,j,k)

! convert to from virtual dry to potential temperature
      theta_star(i,j,k) = theta_star(i,j,k) / ( 1.+                   &
                                     m_star  (i,j,k) * recip_epsilon)

      mcl_star(i,j,k) = r_m_cl_d(i,j,k)
      mcf_star(i,j,k) = r_m_cf_d(i,j,k)
    END DO
  END DO
END DO

DO j=pdims%j_start, pdims%j_end
   DO i=pdims%i_start, pdims%i_end
      k = 0
      theta_star(i,j,k) = R_theta_d(i,j,k) + thetav_ref_pro(i,j,k)      
      m_star  (i,j,k)   = r_m_v_d (i,j,k)
      theta_star(i,j,k) = theta_star(i,j,k) / ( 1.+                   &
                                     m_star  (i,j,k) * recip_epsilon)
      mcl_star(i,j,k)   = r_m_cl_d(i,j,k)
      mcf_star(i,j,k)   = r_m_cf_d(i,j,k)

      k = model_levels
      theta_star(i,j,k) = R_theta_d(i,j,k) + thetav_ref_pro(i,j,k)      
      m_star  (i,j,k)   = r_m_v_d (i,j,k)
      theta_star(i,j,k) = theta_star(i,j,k) / ( 1.+                   &
                                     m_star  (i,j,k) * recip_epsilon)
      mcl_star(i,j,k)   = r_m_cl_d(i,j,k)
      mcf_star(i,j,k)   = r_m_cf_d(i,j,k)
   END DO
END DO

IF(L_mcr_qcf2)   THEN
   DO k = 0, model_levels
      DO j=pdims%j_start, pdims%j_end
         DO i=pdims%i_start, pdims%i_end
            mcf2_star  (i,j,k) = r_m_cf2_d(i,j,k)
         END DO
      END DO
   END DO
END IF

IF(L_mcr_qrain)  THEN
   DO k = 0, model_levels
      DO j=pdims%j_start, pdims%j_end
         DO i=pdims%i_start, pdims%i_end
            mrain_star (i,j,k) = r_m_r_d  (i,j,k)
         END DO
      END DO
   END DO
END IF

IF(L_mcr_qgraup) THEN
   DO k = 0, model_levels
      DO j=pdims%j_start, pdims%j_end
         DO i=pdims%i_start, pdims%i_end
            mgraup_star(i,j,k) = r_m_gr_d (i,j,k)
         END DO
      END DO
   END DO
END IF


IF ( PrintStatus == PrStatus_Diag) THEN
  CALL eg_swap_bounds(m_star,tdims_s,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(mcl_star,tdims_s,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(mcf_star,tdims_s,fld_type_p,.FALSE.)
  IF(L_mcr_qcf2) THEN
    CALL eg_swap_bounds(mcf2_star,tdims_s,fld_type_p,.FALSE.)
  END IF
  IF(L_mcr_qrain) THEN
    CALL eg_swap_bounds(mrain_star,tdims_s,fld_type_p,.FALSE.)
  END IF
  IF(L_mcr_qgraup) THEN
    CALL eg_swap_bounds(mgraup_star,tdims_s,fld_type_p,.FALSE.)
  END IF

  max_r_u     = MAXVAL(R_u_p2+un)
  min_r_u     = MINVAL(R_u_p2+un)
  max_r_v     = MAXVAL(R_v_p2+vn)
  min_r_v     = MINVAL(R_v_p2+vn)
  max_r_w     = MAXVAL(R_w_p2+wn(1:row_length,1:rows,:))
  min_r_w     = MINVAL(R_w_p2+wn(1:row_length,1:rows,:))
  max_r_theta = MAXVAL(theta_star(1:row_length,1:rows,:))
  min_r_theta = MINVAL(theta_star(1:row_length,1:rows,:))
  max_u       = MAXVAL(u)
  min_u       = MINVAL(u)
  max_v       = MAXVAL(v)
  min_v       = MINVAL(v)
  max_w       = MAXVAL(w)
  min_w       = MINVAL(w)
  max_thetav  = MAXVAL(thetav)
  min_thetav  = MINVAL(thetav)
  max_m_star   = MAXVAL(m_star)
  max_mcl_star = MAXVAL(mcl_star)
  max_mcf_star = MAXVAL(mcf_star)
  min_m_star   = MINVAL(m_star)
  min_mcl_star = MINVAL(mcl_star)
  min_mcf_star = MINVAL(mcf_star)

  IF(L_mcr_qcf2) THEN
    max_mcf2_star    = MAXVAL(mcf2_star)
    min_mcf2_star    = MINVAL(mcf2_star)
  END IF

  IF(L_mcr_qrain) THEN
    max_mrain_star  = MAXVAL(mrain_star)
    min_mrain_star  = MINVAL(mrain_star)
  END IF

  IF(L_mcr_qgraup) THEN
    max_mgraup_star = MAXVAL(mgraup_star)
    min_mgraup_star = MINVAL(mgraup_star)
  END IF

  CALL gc_rmax(1,n_proc,ierr,max_r_u)
  CALL gc_rmax(1,n_proc,ierr,max_r_v)
  CALL gc_rmax(1,n_proc,ierr,max_r_w)
  CALL gc_rmax(1,n_proc,ierr,max_r_theta)
  CALL gc_rmin(1,n_proc,ierr,min_r_u)
  CALL gc_rmin(1,n_proc,ierr,min_r_v)
  CALL gc_rmin(1,n_proc,ierr,min_r_w)
  CALL gc_rmin(1,n_proc,ierr,min_r_theta)
  CALL gc_rmax(1,n_proc,ierr,max_u)
  CALL gc_rmax(1,n_proc,ierr,max_v)
  CALL gc_rmax(1,n_proc,ierr,max_w)
  CALL gc_rmax(1,n_proc,ierr,max_thetav)
  CALL gc_rmin(1,n_proc,ierr,min_u)
  CALL gc_rmin(1,n_proc,ierr,min_v)
  CALL gc_rmin(1,n_proc,ierr,min_w)
  CALL gc_rmin(1,n_proc,ierr,min_thetav)

  CALL gc_rmax(1,n_proc,ierr,max_m_star)
  CALL gc_rmax(1,n_proc,ierr,max_mcl_star)
  CALL gc_rmax(1,n_proc,ierr,max_mcf_star)
  CALL gc_rmin(1,n_proc,ierr,min_m_star)
  CALL gc_rmin(1,n_proc,ierr,min_mcl_star)
  CALL gc_rmin(1,n_proc,ierr,min_mcf_star)


  IF(L_mcr_qcf2) THEN
    CALL gc_rmax(1,n_proc,ierr,max_mcf2_star)
    CALL gc_rmin(1,n_proc,ierr,min_mcf2_star)
  END IF

  IF(L_mcr_qrain) THEN
    CALL gc_rmax(1,n_proc,ierr,max_mrain_star)
    CALL gc_rmin(1,n_proc,ierr,min_mrain_star)
  END IF

  IF(L_mcr_qgraup) THEN
    CALL gc_rmax(1,n_proc,ierr,max_mgraup_star)
    CALL gc_rmin(1,n_proc,ierr,min_mgraup_star)
  END IF

  IF ( me == 0 ) THEN

    WRITE(6,fmt='(A)') '=============================================='
    WRITE(6,fmt='(A)') 'F1SP predictor (max, min):'
    WRITE(6,fmt='(A,2E15.5)') 'U^f1sp          :',max_r_u,min_r_u
    WRITE(6,fmt='(A,2E15.5)') 'V^f1sp          :',max_r_v,min_r_v
    WRITE(6,fmt='(A,2E15.5)') 'W^f1sp          :',max_r_w,min_r_w
    WRITE(6,fmt='(A,2E15.5)') 'theta^f1sp      :',max_r_theta,min_r_theta
    WRITE(6,fmt='(A,2E15.5)') 'm_star^f1sp     :',max_m_star,min_m_star
    WRITE(6,fmt='(A,2E15.5)') 'mcl_star^f1sp   :',max_mcl_star,min_mcl_star
    WRITE(6,fmt='(A,2E15.5)') 'mcf_star^f1sp   :',max_mcf_star,min_mcf_star

    IF(L_mcr_qcf2) WRITE(6,fmt='(A,2E15.5)')   'mcf2_star^f1sp  :',max_mcf2_star,min_mcf2_star 
    IF(L_mcr_qrain) WRITE(6,fmt='(A,2E15.5)')  'mrain_star^f1sp :',max_mrain_star,min_mrain_star
    IF(L_mcr_qgraup) WRITE(6,fmt='(A,2E15.5)') 'mgraup_star^f1sp:',max_mgraup_star,min_mgraup_star

  END IF
 

  max_r_u     = MAXVAL(R_u_p2)
  min_r_u     = MINVAL(R_u_p2)
  max_r_v     = MAXVAL(R_v_p2)
  min_r_v     = MINVAL(R_v_p2)
  max_r_w     = MAXVAL(R_w_p2)
  min_r_w     = MINVAL(R_w_p2)

  CALL gc_rmax(1,n_proc,ierr,max_r_u)
  CALL gc_rmax(1,n_proc,ierr,max_r_v)
  CALL gc_rmax(1,n_proc,ierr,max_r_w)
  CALL gc_rmin(1,n_proc,ierr,min_r_u)
  CALL gc_rmin(1,n_proc,ierr,min_r_v)
  CALL gc_rmin(1,n_proc,ierr,min_r_w)


  IF ( me == 0 ) THEN

    WRITE(6,fmt='(A,2E15.5)') 'U               :',max_u,min_u
    WRITE(6,fmt='(A,2E15.5)') 'V               :',max_v,min_v
    WRITE(6,fmt='(A,2E15.5)') 'W               :',max_w,min_w
    WRITE(6,fmt='(A,2E15.5)') 'theta           :',max_thetav,min_thetav
    WRITE(6,fmt='(A,2E15.5)') 'R_u             :',max_r_u,min_r_u
    WRITE(6,fmt='(A,2E15.5)') 'R_v             :',max_r_v,min_r_v
    WRITE(6,fmt='(A,2E15.5)') 'R_w             :',max_r_w,min_r_w
    write(6,fmt='(A)') '=============================================='

  END IF

END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(  psi_w_surf,  SIZE(psi_w_surf),  'psiws')

IF (lhook) CALL dr_hook('EG_F1SP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_f1sp
END MODULE eg_f1sp_mod
