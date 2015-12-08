! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_helm_rhs_star_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_helm_rhs_star(R_u,R_v,R_w,R_etadot,R_p, R_theta, R_rho,    &
                                  R_m_v_d, R_m_cl_d, R_m_cf_d,                 &
                                  R_m_r_d, R_m_gr_d, R_m_cf2_d,                &
                                  u, v, w, exner, exner_star,                  &
                                  exner_prime, rho, thetav, etadot, m_v, m_cl, &
                                  m_cf, m_r, m_gr, m_cf2, R_w_d,  del_rho, Ih, &
                                  timestep,  alpha_w, row_length, rows,        &
                                  n_rows, model_levels,inner_it_val,S_u,S_v,   &
                                  S_w,S_thetav,S_m_v,S_m_cl,S_m_cf,            &
                                  S_m_cf2,S_m_r,S_m_gr,psi_w_surf, psi_w_lid)

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
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

      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE gravity_mod
      USE eg_parameters_mod, ONLY : pole_consts,l_rho_av_zz

      USE UM_parvars, ONLY:  offx,offy,halo_i,halo_j,at_extremity
      USE proc_info_mod,  ONLY:  model_domain

      USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels
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
      INTEGER, INTENT(IN) ::  row_length, rows, n_rows, model_levels
! Time step, hydrostatic flag  and off-centre weights

      REAL, INTENT(IN) :: del_rho, Ih, timestep, alpha_w

! sources at departure points

      REAL, INTENT(IN) ::                                                      &
       R_w_d(row_length,rows,0:model_levels)
     
      REAL                                                                     &
      psi_w_surf(row_length,rows),                                             &
      psi_w_lid (row_length,rows)

! Output 

      REAL, INTENT(OUT) ::                                                     &
       R_u(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),             &
       R_v(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),           &
       R_w(row_length,rows,0:model_levels),                                    &
       R_etadot(row_length,rows,0:model_levels),                               &
       R_theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
       R_rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
       R_p(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


       REAL, INTENT(INOUT) ::                                                  &
        exner_star(1-offx:row_length+offx, 1-offy:rows+offy)


! time level (n+1) fields
      REAL, INTENT(IN) ::                                                      &
       u(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),               &
       v(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),             &
       w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),              &
       exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),          &
       exner_prime(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
       rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),              &
       thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

      REAL, INTENT(INOUT) ::                                                   &
           m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
          m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
          m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
           m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
          m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
         m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

      REAL, INTENT(INOUT) ::                                                   &
        etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)


      REAL, INTENT(IN) ::                                                      &
           R_m_v_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
          R_m_cl_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
          R_m_cf_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
           R_m_r_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
          R_m_gr_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
         R_m_cf2_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

      INTEGER, INTENT(IN) :: inner_it_val

! Local variables


      INTEGER :: i, j, k, ktmp
      LOGICAL :: L_call_from_solver
      REAL    :: ugrad, vgrad, wgrad, T1, T2, P1

      REAL  dummy(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

      REAL    :: rdxi1, rdxi2, rdxi3

!     Fast Physics source terms                 
      REAL                                                                     &
             S_u( -offx:row_length-1+offx,1-offy:  rows  +offy,model_levels),  &
             S_v(1-offx:row_length  +offx, -offy:n_rows-1+offy,model_levels),  &
             S_w(row_length,rows,0:model_levels),                              &
        S_thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
           S_m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
           S_m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
         S_m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)
         
      REAL :: rho_av, rho_ref_term, rho_switch   

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_HELM_RHS_STAR',zhook_in,zhook_handle)

      IF (INTEGRITY_TEST) THEN
!       check all input fields

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
      END IF
      
      rho_switch = 1.0
      IF ( l_rho_av_zz ) rho_switch = 0.0    
!
! Add the fast physics source terms
!

!$OMP PARALLEL DO PRIVATE(i,j)
      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               R_u(i,j,k) = S_u(i,j,k)
            END DO
         END DO

         DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
               R_v(i,j,k) = S_v(i,j,k)
            END DO
         END DO

         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
               R_w(i,j,k) = S_w(i,j,k)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

     R_w(:,:,0) = S_w(:,:,0)
     R_theta (:,:,:) = S_thetav(:,:,:)

      IF( inner_it_val == 1 ) THEN

         CALL update_moisture(m_v, m_cl, m_cf, m_r, m_gr, m_cf2,              &
                              R_m_v_d, R_m_cl_d, R_m_cf_d,                    &
                              R_m_r_d, R_m_gr_d, R_m_cf2_d,                   &
                              S_m_v, S_m_cl, S_m_cf,                          &
                              S_m_r, S_m_gr,S_m_cf2)
      ENDIF


! Adjust arrival point values to account for extra terms
 
!------------------------------------------------------------------------------
! Compute Psi_w_surf to be used in eg_sisl_init to compute bottom R_w
!------------------------------------------------------------------------------

      k = 0
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            psi_w_surf(i,j) = ( Ih*(                                           &
                                (h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*             &
                                dxi1_xi3(i,j,k)*(intw_u2p(i,1)*u(i-1,j,k+1)    &
                                                 +intw_u2p(i,2)*u(i,j,k+1)) +  &
                                (h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*             &
                                dxi2_xi3(i,j,k)*(intw_v2p(j,1)*v(i,j-1,k+1)    &
                                                 +intw_v2p(j,2)*v(i,j,k+1))) - &
                                      R_w_d(i,j,k) )/(alpha_w*timestep)

            psi_w_lid(i,j) = -R_w_d(i,j,model_levels)/(alpha_w*timestep)

         END DO
      END DO

      L_call_from_solver = .true.

      CALL eg_sisl_init(                                                       &
               row_length, rows, n_rows, model_levels,  .FALSE.,               &
               L_call_from_solver,.FALSE., Ih, g_theta, u, v, w,               &
               thetav, rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2, exner,          &
               exner_star, R_u, R_v, R_w, R_theta, R_rho, S_m_v, S_m_cl,       &
               S_m_cf, S_m_r, S_m_gr, S_m_cf2, dummy,psi_w_surf, psi_w_lid )

! Calculate "tilde" quantities for implicit treatment of Coriolis

! R_u_star
      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               rdxi1      = 1.0/(xi1_p(i+1)-xi1_p(i))
               R_u(i,j,k) = R_u(i,j,k) -u(i,j,k) +                             &
                            HM_u(i,j,k)*(                                      &
                               exner_prime(i+1,j,k)*deta_xi3(i+1,j,k) -        &
                               exner_prime(i,j,k)*deta_xi3(i,j,k))             &
                            *rdxi1

            END DO
         END DO
      END DO

! R_v_star
      DO k = 1, model_levels
         DO j = vdims%j_start, vdims%j_end
            rdxi2 = 1.0/(xi2_p(j+1)-xi2_p(j))
            DO i = vdims%i_start, vdims%i_end
               R_v(i,j,k) = R_v(i,j,k) -v(i,j,k) +                             &
                            HM_v(i,j,k)*(                                      &
                               exner_prime(i,j+1,k)*deta_xi3(i,j+1,k) -        &
                               exner_prime(i,j,k)*deta_xi3(i,j,k))             &
                            *rdxi2

            END DO
         END DO
      END DO

! R_w_star,  R_etadot_star, R_thetav_star & R_p_star

      DO k = 1, model_levels-1
         rdxi3 = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k))
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end

                R_w(i,j,k) = R_w(i,j,k) - Ih*w(i,j,k) +                        &
                             HM_w(i,j,k)*(                                     &
                              (thetav(i,j,k)/thetav_ref_pro(i,j,k)-1.0) *      &
                              (exner_ref_pro(i,j,k+1)-exner_ref_pro(i,j,k)) +  &
                              (exner_prime(i,j,k+1) - exner_prime(i,j,k)) )    &
                            *rdxi3

! Now etadot: use method in eg_SISL_Init() for consistency and better accuracy

                R_etadot(i,j,k) =                                              &
                          -(h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*dxi1_xi3(i,j,k)*  &
                                   ( intw_rho2w(k,1)*(                         &
                                         intw_u2p(i,1)*u(i-1,j,k+1)            &
                                       + intw_u2p(i,2)*u(i,j,k+1) )            &
                                    +intw_rho2w(k,2)*(                         &
                                         intw_u2p(i,1)*u(i-1,j,k)              &
                                       + intw_u2p(i,2)*u(i,j,k) ) )            &
                          -(h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*dxi2_xi3(i,j,k)*  &
                                   ( intw_rho2w(k,1)*(                         &
                                         intw_v2p(j,1)*v(i,j-1,k+1)            &
                                       + intw_v2p(j,2)*v(i,j,k+1) )            &
                                    +intw_rho2w(k,2)*(                         &
                                         intw_v2p(j,1)*v(i,j-1,k)              &
                                       + intw_v2p(j,2)*v(i,j,k) ) )
            END DO
         END DO
      END DO

! Surface BC stored in R_etadot

      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            k = 0
            R_etadot(i,j,k) =                                                  &
                          -(h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*dxi1_xi3(i,j,k)*  &
                                   ( intw_u2p(i,1)*u(i-1,j,k+1) +              &
                                              intw_u2p(i,2)*u(i,j,k+1) )       &
                          -(h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*dxi2_xi3(i,j,k)*  &
                                   ( intw_v2p(j,1)*v(i,j-1,k+1) +              &
                                              intw_v2p(j,2)*v(i,j,k+1) )
           k = model_levels
           R_etadot(i,j,k) = 0.0
           END DO
      END DO   

      DO k = 2, model_levels-1
         rdxi3 = 1.0/(eta_theta_levels(k) - eta_theta_levels(k-1))
         DO j = pdims%j_start, pdims%j_end
            rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
            DO i = pdims%i_start, pdims%i_end
               rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

! exner_prime / exner_ref
               P1 = exner_prime(i,j,k) / exner_ref_pro(i,j,k)

! Ratio of vertically averaged thetav_ref and thetav_np1
               T1 = thetav_ref_eta(i,j,k) /                                    &
                     ( intw_w2rho(k,1) * thetav(i,j,k) +                       &
                       intw_w2rho(k,2) * thetav(i,j,k-1) )

! Vertically averaged ratio of thetav and thetav_ref
               T2 = intw_w2rho(k,1) * ( thetav(i,j,k) /                       &
                                        thetav_ref_pro(i,j,k) ) +             &
                    intw_w2rho(k,2) * ( thetav(i,j,k-1) /                     &
                                        thetav_ref_pro(i,j,k-1) )

! R^star_pi: see eqn 9.25
               R_p(i,j,k) = HM_p(i,j,k) - T1 * (1.0 + P1)**HM_pp    +         &
                             HM_pp * P1 - (T2-1.0)


! Now do density

               T2    = h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)       &
                                      *u(i,j,k)

               T1    = h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*deta_xi3_u(i-1,j,k) &
                                      *u(i-1,j,k)

               ugrad = (T2 - T1)*rdxi1


               T2    = h3_xi2_v(i,j,k)*h1_xi2_v(i,j,k)*deta_xi3_v(i,j,k)       &
                                      *v(i,j,k)

               T1    = h3_xi2_v(i,j-1,k)*h1_xi2_v(i,j-1,k)*deta_xi3_v(i,j-1,k) &
                                      *v(i,j-1,k)

               vgrad = (T2 - T1)*rdxi2


               T2    = h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*h3_p_eta(i,j,k)         &
                                  *deta_xi3_theta(i,j,k)*etadot(i,j,k)

               T1    = h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*h3_p_eta(i,j,k-1)   &
                                  *deta_xi3_theta(i,j,k-1)*etadot(i,j,k-1)

               wgrad = (T2 - T1)*rdxi3
               
               rho_av = intw_w2rho(k,1)                                        &
                      *(intw_rho2w(k,1)*(rho(i,j,k+1)-rho_ref_pro(i,j,k+1))    &
                       +intw_rho2w(k,2)*(rho(i,j,k)  -rho_ref_pro(i,j,k))  )   &
                    + intw_w2rho(k,2)                                          &
                      *(intw_rho2w(k-1,1)*(rho(i,j,k)  -rho_ref_pro(i,j,k))    &
                       +intw_rho2w(k-1,2)*(rho(i,j,k-1)-rho_ref_pro(i,j,k-1))) 
               
               rho_ref_term = rho_switch*(rho(i,j,k)-rho_ref_pro(i,j,k))       &
                            + (1.0 - rho_switch)*rho_av
               
               R_rho(i,j,k) = -del_rho*HM_vol(i,j,k)*                          &
                               rho_ref_term *                                  &
                               ( ugrad + vgrad + wgrad )

            END DO
         END DO
      END DO   
! Surface/top boundaries
      k = 1
      DO ktmp = 1,2
         IF ( ktmp == 2 ) k = model_levels
         rdxi3 = 1.0/(eta_theta_levels(k) - eta_theta_levels(k-1))
         DO j = pdims%j_start, pdims%j_end
            rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
            DO i = pdims%i_start, pdims%i_end
               rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

! exner_prime / exner_ref
               P1 = exner_prime(i,j,k) / exner_ref_pro(i,j,k)

! Ratio of vertically averaged thetav_ref and thetav_np1
               T1 = thetav_ref_eta(i,j,k) /                                    &
                     ( intw_w2rho(k,1) * thetav(i,j,k) +                       &
                       intw_w2rho(k,2) * thetav(i,j,k-1) )

! Vertically averaged ratio of thetav and thetav_ref
               T2 = intw_w2rho(k,1) * ( thetav(i,j,k) /                       &
                                        thetav_ref_pro(i,j,k) ) +             &
                    intw_w2rho(k,2) * ( thetav(i,j,k-1) /                     &
                                        thetav_ref_pro(i,j,k-1) )

! R^star_pi: see eqn 9.25
               R_p(i,j,k) = HM_p(i,j,k) - T1 * (1.0 + P1)**HM_pp    +         &
                             HM_pp * P1 - (T2-1.0)


! Now do density

               T2    = h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)       &
                                      *u(i,j,k)

               T1    = h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*deta_xi3_u(i-1,j,k) &
                                      *u(i-1,j,k)

               ugrad = (T2 - T1)*rdxi1


               T2    = h3_xi2_v(i,j,k)*h1_xi2_v(i,j,k)*deta_xi3_v(i,j,k)       &
                                      *v(i,j,k)

               T1    = h3_xi2_v(i,j-1,k)*h1_xi2_v(i,j-1,k)*deta_xi3_v(i,j-1,k) &
                                      *v(i,j-1,k)

               vgrad = (T2 - T1)*rdxi2


               T2    = h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*h3_p_eta(i,j,k)         &
                                  *deta_xi3_theta(i,j,k)*etadot(i,j,k)

               T1    = h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*h3_p_eta(i,j,k-1)   &
                                  *deta_xi3_theta(i,j,k-1)*etadot(i,j,k-1)

               wgrad = (T2 - T1)*rdxi3
               
               rho_av = rho(i,j,k)-rho_ref_pro(i,j,k)
               
               rho_ref_term = rho_switch*(rho(i,j,k)-rho_ref_pro(i,j,k))       &
                            + (1.0 - rho_switch)*rho_av
               
               R_rho(i,j,k) = -del_rho*HM_vol(i,j,k)*                          &
                               rho_ref_term *                                  &
                               ( ugrad + vgrad + wgrad )

            END DO
         END DO
      END DO


! Fix R_v at the poles

      IF( model_domain == mt_global ) THEN
         IF( at_extremity(PSouth) ) THEN
            CALL eg_v_at_poles(R_u,R_v, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)
         END IF

         IF( at_extremity(PNorth) ) THEN
            CALL eg_v_at_poles(R_u,R_v,-1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)
         END IF
      END IF

      CALL eg_swap_bounds(R_u,udims_s,fld_type_u,.true.)
      CALL eg_swap_bounds(R_v,vdims_s,fld_type_v,.true.)
      
      IF (lhook) CALL dr_hook('EG_HELM_RHS_STAR',zhook_out,zhook_handle)

      END SUBROUTINE eg_helm_rhs_star
      END MODULE eg_helm_rhs_star_mod
