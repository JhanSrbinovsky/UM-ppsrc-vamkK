! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_forcing_mod

IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics


CONTAINS
SUBROUTINE eg_idl_forcing(                                        &
! in data fields.
  u, v, theta,  exner_theta_levels, exner,                        &
! in/out
  theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star,  &
  qgraup_star, r_u, r_v,                                          &
  l_expl_horz_drag, L_HeldSuarez, L_HeldSuarez1_drag,             &
  L_HeldSuarez2_drag,                                             &
! error information
  error_code  )

USE eg_explicit_horz_drag_mod, ONLY: eg_explicit_horz_drag
USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod
USE held_suarez_mod,           ONLY: eg_held_suarez
USE mphys_inputs_mod,          ONLY: l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup

IMPLICIT NONE
!
! Description: Initialises star state for ENDGame when there is no
!              physics. Optionally applies Held-Suarez test of
!              dynamical core (BAMS 75, 1825-1830).
!  
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


! Data arrays
REAL, INTENT (INOUT) ::                                               &
  theta(             tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end )

REAL, INTENT (INOUT) ::                                               &
  theta_star(        tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  q_star(            tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  qcl_star(          tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  qcf_star(          tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  qcf2_star(         tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  qrain_star       ( tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  qgraup_star(       tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  r_u(               udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  r_v(               vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end,                   &
                     vdims_s%k_start:vdims_s%k_end),                  &
  u(                 udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  v(                 vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end ,                  &
                     vdims_s%k_start:vdims_s%k_end),                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  exner(             pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1)

INTEGER error_code

LOGICAL :: l_expl_horz_drag, L_HeldSuarez, L_HeldSuarez2_drag,        &
           l_heldsuarez1_drag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Variables for Held Suarez test case:
INTEGER :: i, j,k
REAL    :: temp1, temp2, SuHe_sigma_cutoff, friction_level,           &
           base_frictional_timescale, SuHe_pole_equ_deltaT,           &
           SuHe_static_stab, SuHe_level_weight,                       &
           SuHe_newtonian_timescale_ka,                               &
           SuHe_newtonian_timescale_ks, newtonian_timescale, sigma

REAL    ::                                                            &
           theta_eq(tdims_s%i_start:tdims_s%i_end,                    &
                    tdims_s%j_start:tdims_s%j_end,                    &
                    tdims_s%k_start:tdims_s%k_end)


! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_IDL_FORCING',zhook_in,zhook_handle)

! initial NULL change atmos_physics1 creates tendencies. So if 
! nothing happens all tendencies are zero.
theta_star  = 0.
q_star      = 0.
qcl_star    = 0.
qcf_star    = 0.
IF (l_mcr_qcf2)   qcf2_star   = 0.
IF (l_mcr_qrain)  qrain_star  = 0.
IF (l_mcr_qgraup) qgraup_star = 0.

IF (l_expl_horz_drag) CALL eg_explicit_horz_drag(r_u,r_v,u,v)
IF (L_HeldSuarez)     CALL eg_held_suarez(                            &
                                   ! in data fields.
                                 u, v, theta,  exner_theta_levels,    &
                                 exner,                               &
                                 ! in/out
                                 theta_star, r_u, r_v,                &
                                 l_expl_horz_drag, L_HeldSuarez,      &
                                 L_HeldSuarez1_drag,                  &
                                 L_HeldSuarez2_drag,                  &
                                 ! error information
                                 error_code  )



IF (lhook) CALL dr_hook('EG_IDL_FORCING',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_forcing


!======================================================================

SUBROUTINE eg_idl_forcing2(                                           &

! in data fields.
  u, v,exner_theta_levels, exner,                                     &
! in/out
  r_u, r_v,                                                           &
  L_HeldSuarez1_drag, L_HeldSuarez2_drag,                             &
! error information
  error_code  )

USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod
USE held_suarez_mod,           ONLY: eg_held_suarez2

IMPLICIT NONE
!
! Description: Initialises star state for ENDGame when there is no
!              physics. Optionally applies Held-Suarez test of
!              dynamical core (BAMS 75, 1825-1830).
!  
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


REAL, INTENT (INOUT) ::                                               &
  r_u(               udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  r_v(               vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end,                   &
                     vdims_s%k_start:vdims_s%k_end),                  &
  u(                 udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  v(                 vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end ,                  &
                     vdims_s%k_start:vdims_s%k_end),                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  exner(             pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1)

INTEGER error_code
LOGICAL l_heldsuarez1_drag,l_heldsuarez2_drag


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_IDL_FORCING2',zhook_in,zhook_handle)

IF (L_HeldSuarez2_drag) CALL eg_held_suarez2(                         &
                                 ! in data fields.
                                 u, v,exner_theta_levels, exner,      &
                                 ! in/out
                                 r_u, r_v,                            &
                                 L_HeldSuarez1_drag,                  &
                                 L_HeldSuarez2_drag,                  &
                                 ! error information
                                 error_code  )

IF (lhook) CALL dr_hook('EG_IDL_FORCING2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_forcing2

END MODULE eg_idl_forcing_mod
