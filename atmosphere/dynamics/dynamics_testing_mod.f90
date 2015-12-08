! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input switches/settings and the check_run_dyntest 
!   routine for logic checking the selected settings. Used by the
!   dynamics to enable testing and some idealised running. 
!   Each switch would usually not be used by the everyday user
!   and is defaulted to the common setting.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics
!
! Method:
!   Switches are initialised to false and read in from the
!   namelists. The module may then be used directly where the switches
!   are needed within the dynamics code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dynamics_testing_mod

USE missing_data_mod, ONLY: rmdi, imdi

USE gcr_input_mod, ONLY:                                         &
    GCR_Diagnostics, GCR_its_avg_step

USE var_end_mod, ONLY:                                           &
    lambda_p_end,phi_p_end,lambda_u_end,phi_v_end,               &
    dlambda_p_end,dphi_p_end,dlambda_u_end,dphi_v_end

USE var_input_mod, ONLY:                                         &
    lam_var, phi_var,var_ratio,lam_ratio,phi_ratio,              &
    phi_frac, lam_frac

IMPLICIT NONE

! the following general users need not alter.
! power users may like to alter these via the namelist run_dyntest

LOGICAL :: L_Physics           = .FALSE.  ! T: physics to be included
LOGICAL :: L_dynamics_only     = .FALSE.  ! T: run without physics
LOGICAL :: L_Run_With_Physics2 = .FALSE.  ! T: physics2 to be included
LOGICAL :: L_exclude_Physics2  = .FALSE.  ! T: physics2 to be excluded
LOGICAL :: L_perturb_IC_theta  = .FALSE.  ! T: perturb theta on ts1
LOGICAL :: L_Backwards         = .FALSE.  ! F Integrate backwards without 
                                          ! physics
LOGICAL :: L_dry               = .FALSE.  ! T run with no moisture  
LOGICAL :: L_adjust_wet        = .FALSE.  ! T perform simple dry&moist 
                                          ! adjustment (ND only)
LOGICAL :: L_idealised_data    = .FALSE.  ! T run idealised problem
LOGICAL :: L_free_slip         = .FALSE.  ! Use free slip lower boundary 
                                          ! condition 
LOGICAL :: L_hydrostatic_EG    = .FALSE.  ! choose to run hydrostatic in EG.

! trapping
LOGICAL :: L_trap_uv    = .FALSE.  ! .true. trap excessive u,v  (ND only)
LOGICAL :: L_trap_w     = .FALSE.  ! .true. trap excessive w    (ND only)
LOGICAL :: L_trap_theta = .FALSE.  ! .true. trap diabatic heating/cooling

INTEGER :: trap_option = imdi     ! trapping and printing options (ND only)
                                  ! 0 = reset, no prints
                                  ! 1 = reset + print message
                                  ! 2 = no reset but print details max/min

INTEGER :: Cw_test_lev = imdi     ! lowest level to test for max courant number

REAL :: max_thinc = rmdi          ! max diabatic heating/cooling for trapping
REAL :: Cw_max    = rmdi          ! max Courant number for w for trapping 
REAL :: uv_limit  = rmdi          ! max horizontal wind for trapping   


NAMELIST/RUN_Dyntest/                                                  &
         L_dynamics_only, L_Backwards, L_hydrostatic_EG,               &
         L_exclude_Physics2, L_perturb_IC_theta,                       &
         L_dry, L_adjust_wet, L_free_slip,                             &
         L_trap_uv, L_trap_w, L_trap_theta,                            &
         Cw_max, uv_limit, max_thinc, trap_option, Cw_test_lev,        &
         GCR_Diagnostics, GCR_its_avg_step,                            &
         lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,           &
         dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,       &     
         lam_var, phi_var,                                             &
         var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac
         
CONTAINS

SUBROUTINE check_run_dyntest()

! Description:
!   Subroutine to apply logic checks and set variables based on the 
!   options selected in the run_dyntest namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE dynamics_input_mod, ONLY: l_endgame, Ih

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_DYNTEST',zhook_in,zhook_handle)

! ---------------------------------------------
! Dynamics logic control - New Dynamics and ENDGAME
! ---------------------------------------------
IF (l_dynamics_only) THEN
    l_physics           = .FALSE.
ELSE
    l_physics           = .TRUE.
END IF

IF (l_exclude_physics2) THEN
    l_run_with_physics2 = .FALSE.
ELSE
    l_run_with_physics2 = .TRUE.
END IF
! ---------------------------------------------
! Dynamics logic control - ENDGAME Only
! ---------------------------------------------
IF (l_endgame) THEN
  IF (l_hydrostatic_EG) THEN
    Ih =0.0
  ELSE
    Ih =1.0
  END IF
END IF

IF (lhook) CALL dr_hook('CHECK_RUN_DYNTEST',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_dyntest

END MODULE dynamics_testing_mod
