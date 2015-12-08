! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE eg_add_pressure_inc_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_add_pressure_inc(exner_prime, exner_np1, exner_star_np1,   &
                 R_v_South, R_v_North,                                         &
                 u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,         &
                 m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,               &
                 m_cf2_np1, R_w_d, R_u, R_v, R_w, R_theta, R_etadot, R_p,      &
                 Ih,n_rows, row_length, rows, model_levels)

      USE um_parvars,          ONLY : offx, offy, halo_i, halo_j
      USE level_heights_mod,   ONLY : eta_theta_levels, eta_rho_levels,        &
                                      xi3_at_theta=>r_theta_levels,            &
                                      xi3_at_rho=>r_rho_levels
      USE timestep_mod,        ONLY : timestep
      USE eg_alpha_mod,        ONLY : alpha_w
      USE atmos_constants_mod, ONLY : p_zero, kappa, R
      USE eg_vert_damp_mod,    ONLY : mu_w
      USE yomhook,             ONLY : lhook, dr_hook
      USE parkind1,            ONLY : jprb, jpim
      USE eg_helmholtz_mod
      USE eg_coriolis_star_mod
      USE eg_v_at_poles_mod
      USE eg_calc_p_star_mod
      USE proc_info_mod,       ONLY : global_row_length,                       &
                                      at_extremity,                            &
                                      gc_proc_row_group,                       &
                                      model_domain
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE UM_ParParams
      USE metric_terms_mod
      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE gravity_mod
      USE eg_parameters_mod,   ONLY : pole_consts
      USE eg_swap_bounds_mod

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

      INTEGER,       INTENT(IN)    :: row_length
      INTEGER,       INTENT(IN)    :: rows
      INTEGER,       INTENT(IN)    :: model_levels
      INTEGER,       INTENT(IN)    :: n_rows
      REAL,          INTENT(IN)    :: Ih
   
! Fields at timelevel n+1

      REAL, INTENT(OUT) :: u_np1(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,       &
                                 udims_s%k_start:udims_s%k_end)

      REAL, INTENT(OUT) :: v_np1(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,       &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(OUT) :: w_np1(wdims_s%i_start:wdims_s%i_end,       &
                                 wdims_s%j_start:wdims_s%j_end,       &
                                 wdims_s%k_start:wdims_s%k_end)

      REAL, INTENT(OUT) :: thetav_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) :: rho_np1(pdims_s%i_start:pdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(OUT) :: etadot_np1(wdims_s%i_start:wdims_s%i_end,  &
                                      wdims_s%j_start:wdims_s%j_end,  &
                                      wdims_s%k_start:wdims_s%k_end)

      REAL, INTENT(OUT) :: exner_np1(pdims_s%i_start:pdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,   &
                                     pdims_s%k_start:pdims_s%k_end+1)

      REAL, INTENT(IN) ::     m_v_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::    m_cl_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::    m_cf_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::     m_r_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::    m_gr_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::   m_cf2_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)
                       
      REAL ::          exner_star_np1(pdims_s%i_start:pdims_s%i_end,  &
                                      pdims_s%j_start:pdims_s%j_end)

! RHS terms from time level n

      REAL, INTENT(IN)    :: R_w_d(wdims%i_start:wdims%i_end,         &
                                   wdims%j_start:wdims%j_end,         &
                                   wdims%k_start:wdims%k_end)

      REAL,  INTENT(IN)   :: R_v_North(vdims_s%i_start:vdims_s%i_end, &
                                       vdims_s%k_start:vdims_s%k_end)

      REAL,  INTENT(IN)   :: R_v_South(vdims_s%i_start:vdims_s%i_end, &
                                       vdims_s%k_start:vdims_s%k_end)

! estimates of RHS terms at time level n+1

      REAL, INTENT(INOUT) :: R_u(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,       &
                                 udims_s%k_start:udims_s%k_end)

      REAL, INTENT(INOUT) :: R_v(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,       &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(INOUT) ::   R_w(wdims%i_start:wdims%i_end,         &
                                   wdims%j_start:wdims%j_end,         &
                                   wdims%k_start:wdims%k_end)

      REAL, INTENT(INOUT) ::  R_theta(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::  R_p(pdims_s%i_start:wdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(INOUT) :: R_etadot(wdims%i_start:wdims%i_end,      &
                                      wdims%j_start:wdims%j_end,      &
                                      wdims%k_start:wdims%k_end)  

! Current estimate of the pressure perturbation

      REAL :: exner_prime(pdims_s%i_start:wdims_s%i_end,              &
                          pdims_s%j_start:pdims_s%j_end,              &
                          pdims_s%k_start:pdims_s%k_end)
        
      REAL :: w_surf(wdims_s%i_start:wdims_s%i_end,       &
                     wdims_s%j_start:wdims_s%j_end)

      REAL :: w_lid(wdims_s%i_start:wdims_s%i_end,       &
                    wdims_s%j_start:wdims_s%j_end)

! Local variables

      REAL    :: D, kp, p0
      REAL    :: rdxi1, rdxi2, rdxi3
      REAL    :: T2
      INTEGER :: i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
              
      REAL :: field_inc               
      
      IF (lhook) CALL dr_hook('EG_ADD_PRESSURE',zhook_in,zhook_handle)
      

! Back Subs to get u and v
!$OMP PARALLEL DO PRIVATE(i,j,rdxi1,rdxi2,field_inc)
      DO k = 1, model_levels
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            rdxi1 = 1.0/(xi1_p(i+1)-xi1_p(i))
            field_inc = R_u(i,j,k) - rdxi1*HM_u(i,j,k)*(              &
                         exner_prime(i+1,j,k)*deta_xi3(i+1,j,k)       &
                          -exner_prime(i,j,k)*deta_xi3(i,j,k) )
               
            u_np1(i,j,k) = u_np1(i,j,k) + field_inc               
          END DO
        END DO

        DO j = vdims%j_start, vdims%j_end
          rdxi2 = 1.0/(xi2_p(j+1)-xi2_p(j))
          DO i = vdims%i_start, vdims%i_end
            field_inc = R_v(i,j,k)                                    &
                         - rdxi2*HM_v(i,j,k)*(                        &
                             exner_prime(i,j+1,k)*deta_xi3(i,j+1,k)   &
                            -exner_prime(i,j,k)*deta_xi3(i,j,k) )

            v_np1(i,j,k) = v_np1(i,j,k) + field_inc                
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
      
! Back Subs to get etadot and theta      
!$OMP PARALLEL DO PRIVATE(i,j,rdxi3,field_inc)
      DO k = 1, model_levels-1
        rdxi3 = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k))
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
              
            field_inc = R_w(i,j,k)                                    &
                        - rdxi3*HM_w(i,j,k)*(                         &
                          exner_prime(i,j,k+1)-exner_prime(i,j,k))   
                                
            field_inc = Hlm_Ck(i,j,k)*field_inc
               
            etadot_np1(i,j,k) = etadot_np1(i,j,k) + field_inc
               
!           Alternative (to that below) w back subs not
!           enforcing R_etadot = 0

            w_np1(i,j,k) = w_np1(i,j,k) + Hm_etadot(i,j,k)*field_inc  &
                            - R_etadot(i,j,k)
               
            R_theta(i,j,k)    = R_theta(i,j,k)                        &
                              - Hm_theta(i,j,k)*field_inc

            thetav_np1(i,j,k) = R_theta(i,j,k) + thetav_np1(i,j,k)

          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    
      k = model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          thetav_np1(i,j,0) = R_theta(i,j,0) + thetav_np1(i,j,0)
          thetav_np1(i,j,k) = R_theta(i,j,k) + thetav_np1(i,j,k)       
        END DO
      END DO

      IF( model_domain == mt_global ) THEN
        IF( at_extremity(PSouth) ) THEN
          CALL eg_v_at_poles(u_np1,v_np1, 1.0, udims%j_start,         &
                             vdims%j_start,                           &
                             udims_s,vdims_s)
        END IF

        IF( at_extremity(PNorth) ) THEN
          CALL eg_v_at_poles(u_np1,v_np1,-1.0, udims%j_end,           &
                             vdims%j_end,                             &
                             udims_s,vdims_s)
        END IF
      END IF

      CALL eg_swap_bounds(u_np1,udims_s,fld_type_u,.true.)
      CALL eg_swap_bounds(v_np1,vdims_s,fld_type_v,.true.)

! Surface value used for BCs
     etadot_np1(:,:,0)            = 0.0
     etadot_np1(:,:,model_levels) = 0.0
     w_np1(:,:,model_levels)      = 0.0

     D = 1.0/(alpha_w*timestep)
     DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end

! w back subs not enforcing R_etadot = 0
         w_np1(i,j,0) = w_np1(i,j,0) - R_etadot(i,j,0)
                            
         w_surf(i,j)   = D*( R_w_d(i,j,0) - Ih*w_np1(i,j,0) )
         w_lid(i,j)    = D*R_w_d(i,j,model_levels)
       END DO
     END DO
      
!------------------------------------------------------------------------------

! Update density using eqn(9.32) of EG2.02 and, finally, pressure.

     kp = 1.0/HM_pp
     P0 = R/p_zero
!$OMP PARALLEL DO PRIVATE(i,j,T2)
     DO k = 1, model_levels
       DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

! R_theta now contains theta'
           T2 = intw_w2rho(k,1)*R_theta(i,j,k)                        &
                                   /thetav_ref_pro(i,j,k)             &
                    +intw_w2rho(k,2)*R_theta(i,j,k-1)                 &
                                   /thetav_ref_pro(i,j,k-1)

           rho_np1(i,j,k) = rho_np1(i,j,k)                            &
                                +rho_ref_pro(i,j,k)*(                 &
                        Hm_pp*exner_prime(i,j,k)/exner_ref_pro(i,j,k) &
                                 -T2  + R_p(i,j,k) )
                                 
!          Use exner_prime to update exner and eqn of state to update rho
           exner_np1(i,j,k) = exner_np1(i,j,k) + exner_prime(i,j,k)
                                  
         END DO
       END DO
     END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------------
! Compute surface pressure
!----------------------------------------------------------------------------

     CALL EG_Calc_P_star(                                             &
                       model_levels,  row_length, rows, exner_np1,    &
                       thetav_np1,                                    &
                       m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,&
                       g_theta, exner_star_np1, w_surf, w_lid)
        
      CALL eg_swap_bounds(w_np1,wdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(etadot_np1,wdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(thetav_np1,tdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(rho_np1,pdims_s,fld_type_p,.false.)

      CALL eg_swap_bounds(exner_np1,wdims_s,fld_type_p,.false.)

! IF cyclic LAM fix boundary pressure to balance v-forcing terms

      IF( model_domain == mt_cyclic_lam ) THEN
        DO k = 1, model_levels
          DO i = pdims%i_start, pdims%i_end

! North boundary
            j = vdims%j_end
            D = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
            exner_np1(i,j+1,k) = exner_np1(i,j,k) + R_v_North(i,k)/D

! South boundary
            j = vdims%j_start
            D = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
            exner_np1(i,j,k) = exner_np1(i,j+1,k) - R_v_South(i,k)/D

          END DO
        END DO
      END IF

!-------------------------------------------------------------------------
! Compute "starred" Coriolis terms at new timelevel
!-------------------------------------------------------------------------

      CALL eg_coriolis_star(rho_np1, m_v_np1, m_cl_np1, m_cf_np1,     &
                            m_r_np1, m_gr_np1, m_cf2_np1)

      IF (lhook) CALL dr_hook('EG_ADD_PRESSURE',zhook_out,zhook_handle)

      END SUBROUTINE eg_add_pressure_inc
      END MODULE eg_add_pressure_inc_mod
