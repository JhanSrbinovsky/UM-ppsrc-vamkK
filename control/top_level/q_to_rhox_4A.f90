! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_q_to_rhox_mod
IMPLICIT NONE

CONTAINS
SUBROUTINE eg_q_to_rhox( row_length, rows, model_levels,          &
             wet_levels, l_dry, halo_i, halo_j, offx, offy,       &
             intw_w2rho,q, qcl, qcf, rho, m_v, m_cl, m_cf )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod

IMPLICIT NONE
!
! Description: compute rho_v, rho_cl, rho_cf from specific humidity
!              fields
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER, INTENT(IN) :: row_length, rows, model_levels,          &
                       wet_levels, halo_i, halo_j, offx, offy

LOGICAL, INTENT(IN) :: l_dry

REAL, INTENT(IN) :: intw_w2rho(model_levels,2)

REAL, INTENT(IN) ::                                               &
  q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        0:model_levels)                                           &
, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
        0:model_levels)                                           &
, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
        0:model_levels)

REAL                                                              &
  rho(1-offx:row_length+offx,1-offy:rows+offy,                    &
      model_levels )                                              &
, m_v(1-offx:row_length+offx,1-offy:rows+offy,                    &
      0:model_levels )                                            &
, m_cl(1-offx:row_length+offx,1-offy:rows+offy,                   &
       0:model_levels )                                           &
, m_cf(1-offx:row_length+offx,1-offy:rows+offy,                   &
       0:model_levels )

! Local Variables

INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_Q_TO_RHOX',zhook_in,zhook_handle)

! Convert density to dry density and then qX to rhoX

  DO k = 0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        m_v(i,j,k)  = 0.0
        m_cl(i,j,k) = 0.0
        m_cf(i,j,k) = 0.0
      END DO
    END DO
  END DO

IF (lhook) CALL dr_hook('EG_Q_TO_RHOX',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_q_to_rhox
END MODULE eg_q_to_rhox_mod
