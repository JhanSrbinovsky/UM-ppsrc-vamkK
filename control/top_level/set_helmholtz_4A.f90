! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_set_helmholtz_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_set_helmholtz(                                      &
       row_length, rows, n_rows, model_levels, halo_i, halo_j,    &
       offx, offy, l_slice, l_test_tracer,                        &
       l_cartesian, l_shallow, l_inc_solver,                      &
       f1_comp, f2_comp, f3_comp, ih)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE proc_info_mod,       ONLY : model_domain
USE timestep_mod,        ONLY : timestep
USE atmos_constants_mod, ONLY : r, cp, kappa, p_zero
USE level_heights_mod,   ONLY : eta_theta_levels,                     &
                                eta_rho_levels,                       &
                                xi3_at_u_w=>r_at_u_w,                 &
                                xi3_at_v_w=>r_at_v_w,                 &
                                xi3_at_theta => r_theta_levels,       &
                                xi3_at_rho=>r_rho_levels

USE eg_helmholtz_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE eg_alpha_mod
USE Field_Types
USE atm_fields_bounds_mod
USE metric_terms_mod
USE eg_swap_bounds_mod
USE domain_params

IMPLICIT NONE
!
! Description: calculate constants required for the
!              solution to the  Helmholtz problem.
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

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,    &
                       halo_i, halo_j, offx, offy
!
LOGICAL, INTENT(IN) :: l_slice, l_test_tracer, l_cartesian, l_shallow


LOGICAL, INTENT(IN) :: l_inc_solver

! L_test_tracer allows passive tracer testing with slice

!
! SI time weights

REAL, INTENT(IN) ::  ih

! Coriolis parameters

REAL, INTENT(IN)  ::                                              &
  f1_comp (row_length,rows),                                      &
  f2_comp (row_length,rows),                                      &
  f3_comp (row_length,rows)

INTEGER :: i, j, k

REAL :: xi3_theta_avg1, xi3_theta_avg2, dx, dy, dz, deta, tmp
REAL :: rdeta, rdxi1, rdxi2, T1, T2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

LOGICAL, SAVE :: first_call = .TRUE.

IF (lhook) CALL dr_hook('EG_SET_HELMHOLTZ',zhook_in,zhook_handle)

! Compute thetav_ref_eta

! multiplying by one provides a different answer on the IBM!
  DO k=1, model_levels
    thetav_ref_eta(:,:,k) = intw_w2rho(k,1)*thetav_ref_pro(:,:,k)   &
                           +intw_w2rho(k,2)*thetav_ref_pro(:,:,k-1) 
  END DO

! This is HM_pi_pi
hm_pp = (1.0-kappa)/kappa

!------------------------------------------------------------------------!
! Compute HM_theta.
!------------------------------------------------------------------------!

DO k=1, model_levels-1
  rdeta = alpha_theta*timestep/(eta_rho_levels(k+1)-eta_rho_levels(k))
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end

        T1 = intw_w2rho(k,1)*thetav_ref_pro(i,j,k)                                   &
            +intw_w2rho(k,2)*thetav_ref_pro(i,j,k-1)

        T2 = intw_w2rho(k+1,1)*thetav_ref_pro(i,j,k+1)                               &
            +intw_w2rho(k+1,2)*thetav_ref_pro(i,j,k)

        hm_theta(i,j,k) = (T2 - T1) *rdeta
    END DO
  END DO
END DO

hm_theta(:,:,0)            = 0.0
hm_theta(:,:,model_levels) = 0.0

!------------------------------------------------------------------------!
! Compute HM_u, HM_v, HM_rhox, HM_rhoy, HM_p, HM_pp
!------------------------------------------------------------------------!

DO k=1, model_levels
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end

       tmp = h1_xi1_u(i,j,k)*deta_xi3_u(i,j,k)

       t2     =  intw_p2u(i,1)*thetav_ref_eta(i  ,j,k)                &
                +intw_p2u(i,2)*thetav_ref_eta(i+1,j,k)

       hm_u(i,j,k) = alpha_u*timestep*cp*t2/tmp                    

       t2     =  intw_p2u(i,1)*rho_ref_pro(i  ,j,k)                   &
                +intw_p2u(i,2)*rho_ref_pro(i+1,j,k)

       hm_rhox(i,j,k) = h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*                &
                        deta_xi3_u(i,j,k)*t2                 

    END DO
  END DO

  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
            
       tmp    = h2_xi2_v(i,j,k)*deta_xi3_v(i,j,k)

       t2     =  intw_p2v(j,1)*thetav_ref_eta(i,j  ,k)                &
                +intw_p2v(j,2)*thetav_ref_eta(i,j+1,k)

       hm_v(i,j,k) = alpha_v*timestep*cp*t2/tmp    
       
       t2     =  intw_p2v(j,1)*rho_ref_pro(i,j  ,k)                  &
                +intw_p2v(j,2)*rho_ref_pro(i,j+1,k)

       hm_rhoy(i,j,k) = h3_xi2_v(i,j,k)*h1_xi2_v(i,j,k)*               &
                        deta_xi3_v(i,j,k)*t2                              
    END DO
  END DO

! This is HM_pi_rho
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
                                         
      hm_p(i,j,k) = r*rho_ref_pro(i,j,k)*thetav_ref_eta(i,j,k)        &
                    /(p_zero*exner_ref_pro(i,j,k)**hm_pp)
    END DO
  END DO
END DO

! Bottom boundary

DO j=pdims%j_start, pdims%j_end
   DO i=pdims%i_start, pdims%i_end
      hm_p(i,j,0) = r*rho_ref_pro(i,j,1)*thetav_ref_pro(i,j,0)        &
                   /(p_zero*exner_ref_pro(i,j,1)**hm_pp)
   END DO
END DO

! cjs 110609 Should be 1.0
!      HM_p(:,:,:) = 1.0


!------------------------------------------------------------------------!
! Compute rho_ref_eta needed by HM_rhoz below.
!------------------------------------------------------------------------!

! cjs 170909 Upper and lower boundary values were missing
! Use the assumption: rho_ref_pro does not change below k=1/2 and
!                     above k=N-1/2.

rho_ref_eta(:,:,0) = rho_ref_pro(:,:,1)

DO k=1, model_levels-1
  rho_ref_eta(:,:,k) = intw_rho2w(k,1)*rho_ref_pro(:,:,k+1) +     &
                       intw_rho2w(k,2)*rho_ref_pro(:,:,k)
END DO

rho_ref_eta(:,:,model_levels) = rho_ref_pro(:,:,model_levels)

!------------------------------------------------------------------------!
! Compute HM_rhoz
!------------------------------------------------------------------------!

! eta_dot=0 at k=0 and k=N so HM_rhoz term not really needed there

hm_rhoz(:,:,0) = 0.0
hm_rhoz(:,:,model_levels) = 0.0

DO k=1, model_levels-1
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
       hm_rhoz(i,j,k) = h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*          &
                        h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*    &
                        rho_ref_eta(i,j,k)
    END DO
  END DO
END DO

!------------------------------------------------------------------------!

! Compute HM_w, HM_etadot, HM_vol (HM_V)
! NOTE: HM_w, HM_etadot not need to be defined on k=model_levels
!       see eqns (9.20), (9.22)

!------------------------------------------------------------------------!

DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      hm_w(i,j,k) = alpha_w*timestep*cp*thetav_ref_pro(i,j,k)     &
                    /(h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k))

      hm_etadot(i,j,k) = h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)
    END DO
  END DO
END DO

DO k=1, model_levels-1
  rdeta = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k))
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      hm_b(i,j,k) = (exner_ref_pro(i,j,k+1)-exner_ref_pro(i,j,k)) &
                    * rdeta/thetav_ref_pro(i,j,k)
    END DO
  END DO
END DO
k = model_levels
DO j=pdims%j_start, pdims%j_end
   DO i=pdims%i_start, pdims%i_end
      hm_b(i,j,0) = 0.0 !hm_b(i,j,1)
      hm_b(i,j,k) = 0.0 !hm_b(i,j,k-1)
   END DO
END DO

IF ( l_slice .OR. l_test_tracer ) THEN

  DO k=1, model_levels
    dz = eta_theta_levels(k)-eta_theta_levels(k-1)
    DO j=pdims%j_start, pdims%j_end
      dy = xi2_v(j)-xi2_v(j-1)
      DO i=pdims%i_start, pdims%i_end
        dx  = xi1_u(i)-xi1_u(i-1)
        tmp = h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)
        hm_vol(i,j,k)  = alpha_rho*timestep/tmp
        ec_vol(i,j,k)  = tmp*dx*dy*dz
        ec_area(i,j,k) = h1_p(i,j,k)*h2_p(i,j,k)*dx*dy
      END DO
    END DO
  END DO

ELSE
  ec_area = 0.0
  DO k=1, model_levels
    dz = eta_theta_levels(k)-eta_theta_levels(k-1)
    DO j=pdims%j_start, pdims%j_end
      dy = xi2_v(j)-xi2_v(j-1)
      DO i=pdims%i_start, pdims%i_end
        dx = xi1_u(i)-xi1_u(i-1)
        tmp = h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)
        hm_vol(i,j,k) = alpha_rho*timestep/tmp
        ec_vol(i,j,k) = tmp*dx*dy*dz
      END DO
    END DO
  END DO

END IF

CALL eg_swap_bounds(hm_u,udims_s,fld_type_u,.FALSE.)

CALL eg_swap_bounds(hm_v,vdims_s,fld_type_v,.FALSE.)

CALL eg_swap_bounds(hm_w,wdims_s,fld_type_p,.FALSE.)

IF (first_call)                                                       &
  CALL eg_swap_bounds(hm_etadot,wdims_s,fld_type_p,.FALSE.)

IF (alpha_changed .OR. first_call)                                    &
  CALL eg_swap_bounds(hm_vol,pdims_s,fld_type_p,.FALSE.)

CALL eg_swap_bounds(hm_theta,tdims_s,fld_type_p,.FALSE.)

CALL eg_swap_bounds(hm_p,tdims_s,fld_type_p,.FALSE.)

CALL eg_swap_bounds(hm_rhox,udims_s,fld_type_u,.FALSE.)

CALL eg_swap_bounds(hm_rhoy,vdims_s,fld_type_v,.FALSE.)

CALL eg_swap_bounds(hm_rhoz,wdims_s,fld_type_p,.FALSE.)

IF (integrity_test) THEN

  CALL eg_helmholtz_update_integrity()

  CALL update_hash_m(deta_xi3,        SIZE(deta_xi3),        'dexi3', &
                     deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t', &
                     deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u', &
                     deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v', &
                     dxi1_xi3,        SIZE(dxi1_xi3),        'dx1x3', &
                     dxi2_xi3,        SIZE(dxi2_xi3),        'dx2x3')
END IF

first_call = .FALSE.

IF (lhook) CALL dr_hook('EG_SET_HELMHOLTZ',zhook_out,zhook_handle)

END SUBROUTINE eg_set_helmholtz
END MODULE eg_set_helmholtz_mod
