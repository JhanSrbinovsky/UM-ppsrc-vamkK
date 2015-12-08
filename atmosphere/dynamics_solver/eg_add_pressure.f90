! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE eg_add_pressure_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_add_pressure(exner_prime, exner_np1, exner_star_np1,       &
                 R_v_South, R_v_North,                                         &
                 u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,         &
                 m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,               &
                 m_cf2_np1,                                                    &
                 R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,                      &
                 R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a,                      &
                 R_etadot, R_p_a,                                              &
                 alpha_w, timestep,  l_eliminate_rho,                          &
                 Ih, offx, offy, n_rows, row_length, rows, model_levels)


      USE atmos_constants_mod, ONLY : p_zero, kappa, R
      USE eg_vert_damp_mod, ONLY : mu_w
      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE eg_helmholtz_mod
      USE eg_coriolis_star_mod
      USE eg_v_at_poles_mod
      USE eg_calc_p_star_mod
      USE proc_info_mod, ONLY: global_row_length,                              &
                               at_extremity,                                   &
                               gc_proc_row_group,                              &
                               model_domain,                                   &
                               me
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE UM_ParParams
      USE metric_terms_mod
      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE gravity_mod
      USE eg_parameters_mod, ONLY : pole_consts
      USE eg_swap_bounds_mod
      USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels

      IMPLICIT NONE

!
! Description: 
!         Update values to obtain the timelevel n+1 fields
!
!
! Method: ENDGame formulation version 3.02
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions

      Integer,                          Intent(IN)    :: offx, offy
      Integer,                          Intent(IN)    :: row_length
      Integer,                          Intent(IN)    :: rows
      Integer,                          Intent(IN)    :: model_levels
      Integer,                          Intent(IN)    :: n_rows
      Real,                             Intent(IN)    :: Ih, alpha_w, timestep


! Fields at timelevel n+1

      Real, Intent(OUT)   ::                                                   &
            u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),      &
            v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),    &
            w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
            thetav_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),&
            rho_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels),     &
            etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                &
                       0:model_levels),                                        &
            exner_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

      Real, Intent(IN)    ::                                                   &
              m_v_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
             m_cl_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
             m_cf_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
              m_r_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
             m_gr_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
            m_cf2_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)
                       
      Real :: exner_star_np1(1-offx:row_length+offx,1-offy:rows+offy)

! RHS terms from time level n

      Real, Intent(IN)    ::                                                   &
           R_u_d(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),       &
           R_v_d(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),     &
           R_w_d(row_length,rows,0:model_levels),                              &
       R_theta_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
         R_rho_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      Real,  Intent(IN)   ::                                                   &
         R_v_North(1-offx:row_length+offx,model_levels),                       &
         R_v_South(1-offx:row_length+offx,model_levels)

! estimates of RHS terms at time level n+1

      Real, Intent(IN)    ::                                                   &
       R_u_a(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),           &
       R_v_a(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),         &
       R_w_a(row_length,rows,0:model_levels),                                  &
       R_theta_a(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
       R_rho_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels),          &
       R_p_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

       Real :: R_etadot(row_length,rows,0:model_levels)

! Current estimate of the pressure perturbation

      Real ::                                                                  &
        exner_prime(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
        
      Real :: w_surf(1-offx:row_length+offx,1-offy:rows+offy)
      Real :: w_lid(1-offx:row_length+offx,1-offy:rows+offy)


! Local variables

      Real    :: D, P1, T1
      REAL    :: rdxi1, rdxi2, rdxi3
      Integer :: i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      LOGICAL l_eliminate_rho

      IF (lhook) CALL dr_hook('EG_ADD_PRESSURE',zhook_in,zhook_handle)

! Update to get etadot and hence w

!$OMP PARALLEL DO PRIVATE(i,j,rdxi3) SHARED(model_levels,eta_rho_levels,pdims, &
!$OMP&                               etadot_np1,Hlm_Ck,R_w_a,R_w_d,Ih,         &
!$OMP&                               mu_w,R_etadot,HM_w,exner_prime,R_theta_a, &
!$OMP&                               R_theta_d,exner_ref_pro,thetav_ref_pro,   &
!$OMP&                               HM_etadot,w_np1)
      DO k = 1, model_levels-1
         rdxi3 = 1.0/(eta_rho_levels(k+1)-eta_rho_levels(k))
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end

! Etadot is calculated using eqn(9.43) of EG2.02

               etadot_np1(i,j,k) = Hlm_Ck(i,j,k) * ( R_w_a(i,j,k) +            &
                           R_w_d(i,j,k) + (Ih*(1.0)+mu_w(i,j,k))  &
                           *R_etadot(i,j,k) -                                  &
                           HM_w(i,j,k)*                                        &
                        ((exner_prime(i,j,k+1)-exner_prime(i,j,k))   &
                        + (R_theta_a(i,j,k)+R_theta_d(i,j,k))*                 &
                          (exner_ref_pro(i,j,k+1)-exner_ref_pro(i,j,k))        &
                         /thetav_ref_pro(i,j,k) ) *rdxi3 )

! w is calculated using egn(9.30) of EG2.02

               w_np1(i,j,k) = HM_etadot(i,j,k)*etadot_np1(i,j,k)               &
                                          - R_etadot(i,j,k)

            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

! Surface value used for BCs

      D = 1.0/(alpha_w*timestep)
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
            k = 0
            w_np1(i,j,k) = -R_etadot(i,j,k)
            k = model_levels
            w_np1(i,j,k) = 0.0
            etadot_np1(i,j,k) = 0.0
            w_surf(i,j)   = D*( R_w_d(i,j,0) + Ih*R_etadot(i,j,0) )
                           
            w_lid(i,j)    = D*R_w_d(i,j,model_levels)
         END DO
      END DO


! Now update the u velocity using eqn(9.26) of EG2.02

!$OMP PARALLEL DO PRIVATE(i,j,rdxi1,rdxi2) SHARED(model_levels,udims,xi1_p,    &
!$OMP&            u_np1,R_u_d,R_u_a,HM_u,exner_prime,deta_xi3,                 &
!$OMP&            vdims,xi2_p,R_v_d,R_v_a,HM_v,v_np1)        
      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               rdxi1 = 1.0/(xi1_p(i+1)-xi1_p(i))

               u_np1(i,j,k) = R_u_d(i,j,k) + R_u_a(i,j,k) - HM_u(i,j,k) *      &
               (exner_prime(i+1,j,k)*deta_xi3(i+1,j,k)                         &
                -exner_prime(i,j,k)*deta_xi3(i,j,k) )                          &
               *rdxi1

            END DO
        END DO

! Now update the v velocity using eqn(9.27) of EG2.02

         DO j = vdims%j_start, vdims%j_end
             rdxi2 = 1.0/(xi2_p(j+1)-xi2_p(j))
             DO i = vdims%i_start, vdims%i_end
                v_np1(i,j,k) = R_v_d(i,j,k) + R_v_a(i,j,k) - HM_v(i,j,k) *     &
                (exner_prime(i,j+1,k)*deta_xi3(i,j+1,k)                        &
                 -exner_prime(i,j,k)*deta_xi3(i,j,k) ) *rdxi2
              
             END DO
         END DO
      END DO
!$OMP END PARALLEL DO


! Update thetav using eqn(9.29) of EG2.02

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,pdims,thetav_np1,R_theta_d, &
!$OMP&            R_theta_a,HM_theta,etadot_np1,thetav_ref_pro) 
      DO k = 1, model_levels-1
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
               thetav_np1(i,j,k) = R_theta_d(i,j,k) + R_theta_a(i,j,k)         &
                                  - HM_theta(i,j,k)*etadot_np1(i,j,k)          &
                                  + thetav_ref_pro(i,j,k)

             END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
            k = 0
            thetav_np1(i,j,k) = R_theta_d(i,j,k) + R_theta_a(i,j,k)            &
                              + thetav_ref_pro(i,j,k)

            k = model_levels
            thetav_np1(i,j,k) = R_theta_d(i,j,k) + R_theta_a(i,j,k)            &
                              + thetav_ref_pro(i,j,k)

         END DO
      END DO

! Update the density - but first need to fix u and v haloes

      IF( model_domain == mt_global ) THEN
         IF( at_extremity(PSouth) ) THEN

            CALL eg_v_at_poles(u_np1,v_np1, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)

         END IF

         IF( at_extremity(PNorth) ) THEN

            CALL eg_v_at_poles(u_np1,v_np1, -1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)
         END IF
      END IF

      CALL eg_swap_bounds(u_np1,udims_s,fld_type_u,.true.)
      CALL eg_swap_bounds(v_np1,vdims_s,fld_type_v,.true.)

! Update density using eqn(9.32) of EG2.02 and, finally, pressure.

      IF ( l_eliminate_rho ) THEN
!$OMP PARALLEL DO PRIVATE(i,j,rdxi3,rdxi2,rdxi1, P1, T1) SHARED(model_levels,  &
!$OMP&            eta_theta_levels,pdims,xi2_v,xi1_u,Hm_pp,exner_prime,        &
!$OMP&            exner_ref_pro,intw_w2rho,thetav_np1,thetav_ref_pro,rho_np1,  &
!$OMP&            rho_ref_pro,Hm_p,R_p_a,exner_np1)
        DO k = 1, model_levels
          rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
          DO j = pdims%j_start, pdims%j_end
            rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
            DO i = pdims%i_start, pdims%i_end
              rdxi1          = 1.0/( xi1_u(i) - xi1_u(i-1) )

              P1             = Hm_pp*exner_prime(i,j,k)/exner_ref_pro(i,j,k)
              T1             = intw_w2rho(k,1)*( thetav_np1(i,j,k)             &
                                                  /thetav_ref_pro(i,j,k) )     &
                                +intw_w2rho(k,2)*( thetav_np1(i,j,k-1)         &
                                                  /thetav_ref_pro(i,j,k-1) )   &
                                 -1.0
              rho_np1(i,j,k) = rho_ref_pro(i,j,k)*( 1.0 -                      &
                                   1.0/Hm_p(i,j,k)*( R_p_a(i,j,k) - P1 + T1 ))


              exner_np1(i,j,k) = exner_ref_pro(i,j,k)+exner_prime(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO PRIVATE(i,j,rdxi3,rdxi2,rdxi1) SHARED(model_levels,          &
!$OMP&            eta_theta_levels,pdims,xi2_v,xi1_u,rho_np1,R_rho_d,          &
!$OMP&            R_rho_a,rho_ref_pro,HM_vol,HM_rhox,u_np1,HM_rhoy,            &
!$OMP&            v_np1,HM_rhoz,etadot_np1,exner_np1,exner_ref_pro,exner_prime)
        DO k = 1, model_levels
          rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
          DO j = pdims%j_start, pdims%j_end
            rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
            DO i = pdims%i_start, pdims%i_end
              rdxi1          = 1.0/( xi1_u(i) - xi1_u(i-1) )

                rho_np1(i,j,k) = R_rho_d(i,j,k) + R_rho_a(i,j,k)               &
                                  + rho_ref_pro(i,j,k)                         &
                                  - HM_vol(i,j,k)*(                            &
                                    ( HM_rhox(i,j,k)*u_np1(i,j,k)              &
                                    - HM_rhox(i-1,j,k)*u_np1(i-1,j,k) )        &
                                  *rdxi1                                       &
                                  + ( HM_rhoy(i,j,k)*v_np1(i,j,k)              &
                                    - HM_rhoy(i,j-1,k)*v_np1(i,j-1,k) )        &
                                  *rdxi2                                       &
                                  + (HM_rhoz(i,j,k)*etadot_np1(i,j,k)          &
                                    - HM_rhoz(i,j,k-1)*etadot_np1(i,j,k-1))    &
                                  *rdxi3 )

              exner_np1(i,j,k) = exner_ref_pro(i,j,k)+exner_prime(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

!----------------------------------------------------------------------------
! Compute surface pressure
!----------------------------------------------------------------------------

     CALL EG_Calc_P_star(                                                      &
                       model_levels, row_length, rows, exner_np1,              &
                       thetav_np1,                                             &
                       m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,         &
                       g_theta, exner_star_np1, w_surf, w_lid)
        
      CALL eg_swap_bounds(w_np1,wdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(etadot_np1,wdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(thetav_np1,tdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(rho_np1,pdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(exner_np1,wdims_s,fld_type_p,.false.)


! IF cyclic LAM fix boundary pressure to balance v-forcing terms

!     IF( model_domain == mt_cyclic_lam ) THEN
!        DO k = 1, model_levels
!           DO i = pdims%i_start, pdims%i_end

! North boundary
!              j = vdims%j_end
!              D = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
!              exner_np1(i,j+1,k) = exner_np1(i,j,k) + R_v_North(i,k)/D

! South boundary
!              j = vdims%j_start
!              D = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
!              exner_np1(i,j,k) = exner_np1(i,j+1,k) - R_v_South(i,k)/D

!           END DO
!        END DO
!     END IF


!-------------------------------------------------------------------------
! Compute "starred" Coriolis terms at new timelevel
!-------------------------------------------------------------------------

      CALL eg_coriolis_star(rho_np1, m_v_np1, m_cl_np1, m_cf_np1,     &
                            m_r_np1, m_gr_np1, m_cf2_np1)

      IF (lhook) CALL dr_hook('EG_ADD_PRESSURE',zhook_out,zhook_handle)

      END SUBROUTINE eg_add_pressure
      END MODULE eg_add_pressure_mod
