! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE eg_helm_fixd_rhs_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_helm_fixd_rhs(RHS,eta_rho_levels,row_length, rows,&
                                  model_levels,offx, offy)

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE eg_helmholtz_mod
      USE integrity_mod
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE Field_Types
      USE metric_terms_mod
      USE fields_rhs_mod

      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE eg_swap_bounds_mod

      IMPLICIT NONE

!
! Description: Code to calculate the Fixed RHS terms
!              in the Helmholtz problem
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

      INTEGER,  INTENT(IN)    :: offx, offy
      INTEGER,  INTENT(IN)    :: row_length, rows,  model_levels

      REAL, INTENT(IN) ::   eta_rho_levels(model_levels)

      REAL, INTENT(OUT)        ::                                              &
        RHS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


! Local temporary variables

      REAL    :: Rk, Rkm1(1-offx:row_length+offx,1-offy:rows+offy)
      REAL    :: rdxi1, rdxi2, rdxi3
      INTEGER :: i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('EG_HELM_FIXD_RHS',zhook_in,zhook_handle)


      IF (integrity_test) THEN
!     check all input fields
        CALL check_hash_m(                                            &
                  R_u_d,           SIZE(R_u_d),           'R_u_d',    & !1
                  R_v_d,           SIZE(R_v_d),           'R_v_d',    &
                  R_w_d,           SIZE(R_w_d),           'R_w_d',    &
                  R_theta_d,       SIZE(R_theta_d),       'R_t_d',    &
                  R_rho_d,         SIZE(R_rho_d),         'R_r_d',    &
                  R_p_p_d,         SIZE(R_p_p_d),         'Rpp_d')

          CALL eg_helmholtz_check_integrity()

          CALL  check_hash_m(                                         &
                   deta_xi3,        SIZE(deta_xi3),        'dexi3',   &
                   deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t',   &
                   deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u',   &
                   deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v',   &
                   dxi1_xi3,        SIZE(dxi1_xi3),        'dx1x3',   &
                   dxi2_xi3,        SIZE(dxi2_xi3),        'dx2x3')

        END IF

!CALL eg_swap_bounds( r_u_d, udims_s,fld_type_u,.TRUE.)
!CALL eg_swap_bounds( r_v_d, vdims_s,fld_type_v,.TRUE.)

! Calculate RHS^n eqn(9.45) of EG3.01

      DO k = 1, model_levels
         DO j = pdims%j_start, pdims%j_end
            rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
            DO i = pdims%i_start, pdims%i_end
              rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

               RHS(i,j,k) = -rho_ref_pro(i,j,k)*(                           &
                             intw_w2rho(k,1)*R_theta_d(i,j,k)               &
                                            /thetav_ref_pro(i,j,k)          &
                            +intw_w2rho(k,2)*R_theta_d(i,j,k-1)             &
                                            /thetav_ref_pro(i,j,k-1) )      &
                            /HM_p(i,j,k) - R_rho_d(i,j,k)


               RHS(i,j,k) = RHS(i,j,k) + HM_vol(i,j,k)*(                    &
                                       ( HM_rhox(i,j,k)*R_u_d(i,j,k)        &
                                       - HM_rhox(i-1,j,k)*R_u_d(i-1,j,k)    &
                                       )*rdxi1                              &
                                      +( HM_rhoy(i,j,k)*R_v_d(i,j,k)        &
                                       - HM_rhoy(i,j-1,k)*R_v_d(i,j-1,k)    &
                                       )*rdxi2   )
            END DO
         END DO
      END DO


! These are the terms arising from D_1 in (9.47)

      Rkm1(:,:) = 0.0
      DO k = 1, model_levels-1
         rdxi3 = 1.0/( eta_rho_levels(k+1) - eta_rho_levels(k))
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end

! Use same definintion of D_1(X) as in LHS :
!            D_1(X) = E_k.X_k + F_k.X_(k-1)

! NB. HM_vol is subsumed into Hlm_Ek and Hlm_Fk, see eg_set_helm_lhs.

               Rk = Hlm_Ck(i,j,k) *                                            &
                    ( R_w_d(i,j,k) - HM_w(i,j,k) * R_theta_d(i,j,k) *          &
                      ( exner_ref_pro(i,j,k+1) - exner_ref_pro(i,j,k) )        &
                      *rdxi3/thetav_ref_pro(i,j,k) )

               RHS(i,j,k) = RHS(i,j,k) + Hlm_Ek(i,j,k)*Rk                      &
                                       + Hlm_Fk(i,j,k)*Rkm1(i,j)
               Rkm1(i,j)  = Rk

            END DO
         END DO
      END DO

! Fix up top 

      k = model_levels
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
            RHS(i,j,k) = RHS(i,j,k) + Hlm_Fk(i,j,k)*Rkm1(i,j)
         END DO
      END DO

      IF (integrity_test) CALL update_hash_m(RHS,SIZE(RHS),'rn___')

      IF (lhook) CALL dr_hook('eg_helm_fixd_rhs',zhook_out,zhook_handle)

      END SUBROUTINE eg_helm_fixd_rhs
      END MODULE eg_helm_fixd_rhs_mod
