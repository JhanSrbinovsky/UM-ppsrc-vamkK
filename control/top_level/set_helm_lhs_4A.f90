! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE  eg_set_helm_lhs_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_set_helm_lhs(                                              &
       row_length, rows, n_rows, model_levels, ih)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE eg_vert_damp_mod,      ONLY : mu_w
USE level_heights_mod,     ONLY : eta_theta_levels, eta_rho_levels
USE eg_helmholtz_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE metric_terms_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE UM_ParVars
USE proc_info_mod,     ONLY : model_domain

IMPLICIT NONE
!
! Description:calculates the constant
!              part of the Helmholtz coefficient matrix.
!  
!
! Method: Appendix F, ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments



INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

REAL, INTENT(IN) :: ih

! Local variables

REAL    :: c, T, rdxi1, rdxi2, rdxi3
INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SET_HELM_LHS',zhook_in,zhook_handle)


! Precalculate some useful constants for reuse in the
! RHS of the Helmholtz equation
! In (9.40) we  D1_k = E_k*P_k + F_k*P_(k-1)
! NOTE : (a) We've absorbed the common HM_vol factor into the definitions
!            of Ek and Fk.
!        (b) Version of model using (1/rho)*grad(p) had A_k and B_k
!            coefficients for a more comlicated D2 operator.
IF (integrity_test) hlm_ck(:,:,model_levels) = 0.



!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,C,rdxi3)                  &
!$OMP&         SHARED(eta_rho_levels,                                   &
!$OMP&            HM_etadot,Ih,mu_w,hm_w,hm_theta,hm_b,hlm_ck,          &
!$OMP&            model_levels,pdims) SCHEDULE(STATIC)
DO k = 1, model_levels-1
   rdxi3 = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k))
   DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

         c = HM_etadot(i,j,k)*(Ih+mu_w(i,j,k))                         &
               -hm_w(i,j,k)*hm_theta(i,j,k)*hm_b(i,j,k)

         hlm_ck(i,j,k) = 1.0/c

      END DO
   END DO
END DO
!$OMP END PARALLEL DO

hlm_ck(:,:,0)            = 0.0
hlm_ck(:,:,model_levels) = 0.0

!$OMP PARALLEL DO PRIVATE(i,j,k, C,T,rdxi3) SHARED(model_levels,         &
!$OMP&           eta_theta_levels,pdims,rho_ref_pro,hm_p,hm_vol,      &
!$OMP&           hlm_ek,hm_rhoz,intw_w2rho,hm_theta,thetav_ref_pro,   &
!$OMP&           hlm_fk) DEFAULT(NONE) SCHEDULE(STATIC)
   DO k = 2, model_levels-1
      rdxi3 = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            C             = rho_ref_pro(i,j,k)/hm_p(i,j,k)
            T             = hm_vol(i,j,k)*rdxi3

            hlm_ek(i,j,k) = hm_rhoz(i,j,k)*T                          &
                           +C*intw_w2rho(k,1)*hm_theta(i,j,k)         &
                             /thetav_ref_pro(i,j,k)

            hlm_fk(i,j,k) = -hm_rhoz(i,j,k-1)*T                       &
                            +C*intw_w2rho(k,2)*hm_theta(i,j,k-1)      &
                              /thetav_ref_pro(i,j,k-1)

         END DO
      END DO
   END DO
!$OMP END PARALLEL DO

! Top and bottom values for Ek and Fk

   k = model_levels
   T = (eta_rho_levels(1)-eta_theta_levels(0))
   C = (eta_theta_levels(k)-eta_rho_levels(k))
   rdxi3 = 1.0/( eta_theta_levels(1)-eta_theta_levels(0) )
   rdxi1 = 1.0/( eta_theta_levels(k)-eta_theta_levels(k-1) )
!$OMP PARALLEL DO PRIVATE(i,j) SHARED(pdims,hlm_ek,hm_rhoz,hm_vol,      &
!$OMP&            rho_ref_pro,hm_theta,T,hm_p,thetav_ref_pro,rdxi3,   &
!$OMP&            hlm_fk,C,rdxi1,k) DEFAULT(NONE) SCHEDULE(STATIC)

   DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
         hlm_ek(i,j,1)  = ( hm_rhoz(i,j,1)*hm_vol(i,j,1) +            &
                            rho_ref_pro(i,j,1) * hm_theta(i,j,1) *    &
                          T/(hm_p(i,j,1)*thetav_ref_pro(i,j,1)))*rdxi3
         hlm_fk(i,j,1) = 0.0

         hlm_ek(i,j,k) = 0.0
         hlm_fk(i,j,k) = (-hm_rhoz(i,j,k-1)*hm_vol(i,j,k) +           &
                          rho_ref_pro(i,j,k) *hm_theta(i,j,k-1) *     &
                         C/(hm_p(i,j,k)*thetav_ref_pro(i,j,k-1)))*rdxi1

      END DO
   END DO
!$OMP END PARALLEL DO

! Calculate the constant part of the Helmholtz matrix
! Simplest way is to build it in pieces!

!$OMP PARALLEL DO PRIVATE(i,j,k,rdxi1,rdxi2) SHARED(model_levels,pdims, &
!$OMP&            xi2_v,xi1_u,rho_ref_pro,hm_pp,hm_p,                 &
!$OMP&            exner_ref_pro,hlm_le,hm_rhox,hm_u,xi1_p,hm_vol,     &
!$OMP&            hlm_ln,hm_rhoy,hm_v,xi2_p,hlm_lp,                   &
!$OMP&          deta_xi3,hlm_ls,hlm_lw) DEFAULT(NONE) SCHEDULE(STATIC)
   DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
            rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
         DO i=pdims%i_start, pdims%i_end
            rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

! Main diagonal entry

           hlm_lp(i,j,k) = -rho_ref_pro(i,j,k)*hm_pp/                 &
                        (hm_p(i,j,k)*exner_ref_pro(i,j,k))

! Now do "e" and "w" points with update at "p"

           hlm_le(i,j,k) = hm_rhox(i,j,k)*hm_u(i,j,k)/                &
                           (xi1_p(i+1) - xi1_p(i))
           hlm_lw(i,j,k) = hm_rhox(i-1,j,k)*hm_u(i-1,j,k)/            &
                           (xi1_p(i) - xi1_p(i-1))

           hlm_le(i,j,k) = hm_vol(i,j,k)*hlm_le(i,j,k)*rdxi1
           hlm_lw(i,j,k) = hm_vol(i,j,k)*hlm_lw(i,j,k)*rdxi1

! Now do "n" and "s" points with update at "p"

           hlm_ln(i,j,k) = hm_rhoy(i,j,k)*hm_v(i,j,k)/                &
                           (xi2_p(j+1) - xi2_p(j))
           hlm_ls(i,j,k) = hm_rhoy(i,j-1,k)*hm_v(i,j-1,k)/            &
                           (xi2_p(j) - xi2_p(j-1))

           hlm_ln(i,j,k) = hm_vol(i,j,k)*hlm_ln(i,j,k)*rdxi2
           hlm_ls(i,j,k) = hm_vol(i,j,k)*hlm_ls(i,j,k)*rdxi2

           hlm_lp(i,j,k)=hlm_lp(i,j,k)-((hlm_ln(i,j,k)+hlm_ls(i,j,k))+&
                                        (hlm_le(i,j,k)+hlm_lw(i,j,k)))&
                                       *deta_xi3(i,j,k)

           hlm_le(i,j,k) = hlm_le(i,j,k)*deta_xi3(i+1,j,k)
           hlm_lw(i,j,k) = hlm_lw(i,j,k)*deta_xi3(i-1,j,k)
           hlm_ln(i,j,k) = hlm_ln(i,j,k)*deta_xi3(i,j+1,k)
           hlm_ls(i,j,k) = hlm_ls(i,j,k)*deta_xi3(i,j-1,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO


! Now do "u" and "d" points with final update of "p"
! The top and bottom levels are done separately

!$OMP PARALLEL DO PRIVATE(i,j,k,rdxi3,rdxi2) SHARED(model_levels,pdims, &
!$OMP&            xi2_v,xi1_u,rho_ref_pro,hm_pp,hm_p,                 &
!$OMP&            exner_ref_pro,hlm_le,hm_rhox,hm_u,xi1_p,hm_vol,     &
!$OMP&            hlm_ln,hm_rhoy,hm_v,xi2_p,hlm_lp,                   &
!$OMP&            deta_xi3,hlm_ls,hlm_lw,hlm_ld,hm_w,hlm_lu,          &
!$OMP&            eta_rho_levels,hlm_ek,hlm_fk,hlm_ck) DEFAULT(NONE)  &
!$OMP&        SCHEDULE(STATIC)
DO k = 1, model_levels
 
   IF( k == 1 ) THEN
      rdxi2 = 1.0/(eta_rho_levels(k))
      rdxi3 = 1.0/(eta_rho_levels(k+1)-eta_rho_levels(k))
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            hlm_lu(i,j,k) = hlm_ek(i,j,k)*hlm_ck(i,j,k)*hm_w(i,j,k)*rdxi3

            hlm_lp(i,j,k) = hlm_lp(i,j,k) - hlm_lu(i,j,k)
            hlm_ld(i,j,k) = 0.0

         END DO
      END DO
   ELSE IF( k == model_levels ) THEN
      rdxi2 = 1.0/(eta_rho_levels(k)-eta_rho_levels(k-1))
      rdxi3 = 1.0/(1.-eta_rho_levels(k))
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            hlm_ld(i,j,k) = -hlm_fk(i,j,k)*hlm_ck(i,j,k-1)*hm_w(i,j,k-1)*rdxi2

            hlm_lp(i,j,k) = hlm_lp(i,j,k) - hlm_ld(i,j,k)
            hlm_lu(i,j,k) = 0.0

         END DO
      END DO
   ELSE
      rdxi2 = 1.0/(eta_rho_levels(k)-eta_rho_levels(k-1))
      rdxi3 = 1.0/(eta_rho_levels(k+1)-eta_rho_levels(k))
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            hlm_lu(i,j,k) = hlm_ek(i,j,k)*hlm_ck(i,j,k)*hm_w(i,j,k)*rdxi3
            hlm_ld(i,j,k) =-hlm_fk(i,j,k)*hlm_ck(i,j,k-1)*hm_w(i,j,k-1)*rdxi2

            hlm_lp(i,j,k) = hlm_lp(i,j,k) - (hlm_lu(i,j,k) + hlm_ld(i,j,k))

         END DO
      END DO
   END IF
! Now rescale to make the diagonal unity
! NB. Rescaling is applied to RHS in eg_helm_var_rhs.
   DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
         hlm_lp(i,j,k) = 1.0/hlm_lp(i,j,k)
         hlm_lu(i,j,k) = hlm_lu(i,j,k)*hlm_lp(i,j,k)
         hlm_ld(i,j,k) = hlm_ld(i,j,k)*hlm_lp(i,j,k)
         hlm_le(i,j,k) = hlm_le(i,j,k)*hlm_lp(i,j,k)
         hlm_lw(i,j,k) = hlm_lw(i,j,k)*hlm_lp(i,j,k)
         hlm_ln(i,j,k) = hlm_ln(i,j,k)*hlm_lp(i,j,k)
         hlm_ls(i,j,k) = hlm_ls(i,j,k)*hlm_lp(i,j,k)
      END DO
   END DO
END DO
!$OMP END PARALLEL DO

! Dirichelt BCs for a LAM

IF( model_domain == mt_LAM ) THEN
   IF (neighbour(pwest)  ==  nodomain) THEN
      HLM_Lp(1,:,:) = 1.0
      HLM_Ln(1,:,:) = 0.0
      HLM_Ls(1,:,:) = 0.0
      HLM_Le(1,:,:) = 0.0
      HLM_Lw(1,:,:) = 0.0
      HLM_Lu(1,:,:) = 0.0
      HLM_Ld(1,:,:) = 0.0
   END IF
   IF (neighbour(peast)  ==  nodomain) THEN
      i = pdims%i_end - 1     ! one less p point in LAM's
      HLM_Lp(i,:,:) = 1.0
      HLM_Ln(i,:,:) = 0.0
      HLM_Ls(i,:,:) = 0.0
      HLM_Le(i,:,:) = 0.0
      HLM_Lw(i,:,:) = 0.0
      HLM_Lu(i,:,:) = 0.0
      HLM_Ld(i,:,:) = 0.0

      i = i + 1               ! Force last column to input
      HLM_Lp(i,:,:) = 1.0
      HLM_Ln(i,:,:) = 0.0
      HLM_Ls(i,:,:) = 0.0
      HLM_Le(i,:,:) = 0.0
      HLM_Lw(i,:,:) = 0.0
      HLM_Lu(i,:,:) = 0.0
      HLM_Ld(i,:,:) = 0.0
   END IF
   IF (neighbour(psouth)  ==  nodomain) THEN
      HLM_Lp(:,1,:) = 1.0
      HLM_Ln(:,1,:) = 0.0
      HLM_Ls(:,1,:) = 0.0
      HLM_Le(:,1,:) = 0.0
      HLM_Lw(:,1,:) = 0.0
      HLM_Lu(:,1,:) = 0.0
      HLM_Ld(:,1,:) = 0.0
   END IF
   IF (neighbour(pnorth)  ==  nodomain) THEN
      j = pdims%j_end
      HLM_Lp(:,j,:) = 1.0
      HLM_Ln(:,j,:) = 0.0
      HLM_Ls(:,j,:) = 0.0
      HLM_Le(:,j,:) = 0.0
      HLM_Lw(:,j,:) = 0.0
      HLM_Lu(:,j,:) = 0.0
      HLM_Ld(:,j,:) = 0.0
   END IF
END IF

! Build tridiagonal factorization for SOR routine

k = 1
DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end
      Hu_k(i,j,k) =  hlm_lu(i,j,k)
   END DO
END DO

DO k = 2, model_levels
   DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
         Hd_k(i,j,k) = 1.0/(1.0 - Hlm_Ld(i,j,k)*Hu_k(i,j,k-1))
         Hu_k(i,j,k) = Hlm_Lu(i,j,k)*Hd_k(i,j,k)
      END DO
   END DO
END DO

IF (lhook) CALL dr_hook('EG_SET_HELM_LHS',zhook_out,zhook_handle)

END SUBROUTINE eg_set_helm_lhs
END MODULE eg_set_helm_lhs_mod
