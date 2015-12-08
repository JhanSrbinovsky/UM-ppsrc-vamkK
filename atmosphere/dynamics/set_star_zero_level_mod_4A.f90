! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! set the zeroth level for physics _star variables
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics

MODULE set_star_zero_level_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE set_star_zero_level(                                       &
             theta_star,                                              &
             q_star,                                                  &
             qcl_star,                                                &
             qcf_star,                                                &
             cf_star,                                                 &
             cfl_star,                                                &
             cff_star,                                                &
             qcf2_star,                                               &
             qrain_star,                                              &
             qgraup_star,                                             &
             L_mcr_qgraup,                                            &
             L_mcr_qrain,                                             &
             L_mcr_qcf2)      

USE atm_fields_bounds_mod
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

REAL, INTENT (INOUT) :: theta_star                                    &
                                 (tdims_s%i_start:tdims_s%i_end,      &
                                  tdims_s%j_start:tdims_s%j_end,      &
                                  tdims_s%k_start:tdims_s%k_end) 
REAL, INTENT (INOUT) :: q_star                                        &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)      
REAL, INTENT (INOUT) :: qcl_star                                      &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)      
REAL, INTENT (INOUT) :: qcf_star                                      &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)
REAL, INTENT (INOUT) :: qcf2_star                                     &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end) 
REAL, INTENT (INOUT) :: qrain_star                                    &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end) 
REAL, INTENT (INOUT) :: qgraup_star                                   &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end) 
REAL, INTENT (INOUT) :: cf_star                                       &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)       
REAL, INTENT (INOUT) :: cfl_star                                      &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)
REAL, INTENT (INOUT) :: cff_star                                      &
                                 (qdims_s%i_start:qdims_s%i_end,      &
                                  qdims_s%j_start:qdims_s%j_end,      &
                                  qdims_s%k_start:qdims_s%k_end)

LOGICAL, INTENT (IN) ::  l_mcr_qgraup
LOGICAL, INTENT (IN) ::  l_mcr_qrain
LOGICAL, INTENT (IN) ::  l_mcr_qcf2

INTEGER i,j
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('set_star_zero_level',zhook_in,zhook_handle)


DO j=tdims%j_start, tdims%j_end
  DO i=tdims%i_start, tdims%i_end
    theta_star(i,j,0) = theta_star(i,j,1)
  END DO
END DO

DO j=qdims%j_start, qdims%j_end
  DO i=qdims%i_start, qdims%i_end
    q_star(i,j,0)   = q_star(i,j,1)
    qcl_star(i,j,0) = qcl_star(i,j,1)
    qcf_star(i,j,0) = qcf_star(i,j,1)
    cf_star(i,j,0)  = cf_star(i,j,1)
    cfl_star(i,j,0) = cfl_star(i,j,1)
    cff_star(i,j,0) = cff_star(i,j,1)
  END DO
END DO


IF (l_mcr_qcf2)   qcf2_star  (:, :, 0) = qcf2_star(:, :, 1)
IF (l_mcr_qrain)  qrain_star (:, :, 0) = qrain_star(:, :, 1)
IF (l_mcr_qgraup) qgraup_star(:, :, 0) = qgraup_star(:, :, 1)

IF (lhook) CALL dr_hook('set_star_zero_level',zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_star_zero_level
END MODULE set_star_zero_level_mod
