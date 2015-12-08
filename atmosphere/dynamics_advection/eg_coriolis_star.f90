! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_coriolis_star_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_coriolis_star(rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)


USE parkind1,          ONLY: jpim, jprb       !DrHook
USE yomhook,           ONLY: lhook, dr_hook   !DrHook
USE level_heights_mod, ONLY: eta_theta_levels
USE proc_info_mod,     ONLY: model_domain
USE domain_params,     ONLY: mt_global

USE proc_info_mod,     ONLY : model_domain
USE integrity_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE Field_Types
USE metric_terms_mod
USE coriolis_mod
USE eg_swap_bounds_mod
USE dynamics_input_mod, ONLY : l_simple_coriolis

IMPLICIT NONE
!
! Description:
!  
!
! Method:
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

! Arrays

REAL    rho(pdims_s%i_start:pdims_s%i_end,                            &
            pdims_s%j_start:pdims_s%j_end,                            &
            pdims_s%k_start:pdims_s%k_end),                           &
        m_v(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_cl(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_cf(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
        m_r(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_gr(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
      m_cf2(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end)

! Local
INTEGER i, j, k
REAL    total_rho, rdeta

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_CORIOLIS_STAR',zhook_in,zhook_handle)

IF( .NOT. l_simple_coriolis ) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k,rdeta,total_rho) DEFAULT(NONE)        &
!$OMP& SHARED(pdims_s,eta_theta_levels,pdims,rho,intw_w2rho,m_v,      &
!$OMP& m_cl,m_r,m_cf,m_gr,m_cf2,f1_comp,h1_p,xi1_u,f2_comp,           &
!$OMP& h2_p,xi2_v,f3_comp,h3_p,deta_xi3,f1_star,f2_star,f3_star)      &
!$OMP& SCHEDULE(STATIC)
  DO k = pdims_s%k_start,pdims_s%k_end
    rdeta = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
         total_rho = rho(i,j,k)*( 1.0                                 &
                    + intw_w2rho(k,1)*( m_v(i,j,k) + m_cl(i,j,k)      &
                                      + m_r(i,j,k) + m_cf(i,j,k)      &
                                      + m_gr(i,j,k)+ m_cf2(i,j,k))    &
                    + intw_w2rho(k,2)*(m_v(i,j,k-1)+m_cl(i,j,k-1)     &
                                     + m_r(i,j,k-1)+m_cf(i,j,k-1)     &
                                     + m_gr(i,j,k-1)+m_cf2(i,j,k-1) ))
         total_rho = 1.0/total_rho

         f1_star(i,j,k) = total_rho*f1_comp(i,j)                      &
                            /(h1_p(i,j,k)*(xi1_u(i)-xi1_u(i-1)) )
         f2_star(i,j,k) = total_rho*f2_comp(i,j)                      &
                            /(h2_p(i,j,k)*(xi2_v(j)-xi2_v(j-1)) )
         f3_star(i,j,k) = rdeta*total_rho*f3_comp(i,j)                &
                         /(h3_p(i,j,k)*deta_xi3(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                    f1_star,         SIZE(f1_star),         'f1sta',  &
                    f2_star,         SIZE(f2_star),         'f2sta',  &
                    f3_star,         SIZE(f3_star),         'f3sta')

IF (lhook) CALL dr_hook('EG_CORIOLIS_STAR',zhook_out,zhook_handle)

END SUBROUTINE eg_coriolis_star
END MODULE eg_coriolis_star_mod
