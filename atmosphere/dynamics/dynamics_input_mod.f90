! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input switches/settings as used by the
!   semi_lagrangian code and the check_run_dyn routine for logic 
!   checking the selected settings.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

! Method:
!   Switches are initialised to false and read in from the
!   namelists. The module may then be used directly where the switches
!   are needed within the dynamics code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dynamics_input_mod

USE missing_data_mod, ONLY: rmdi, imdi

USE gcr_input_mod, ONLY:                                         &
    L_GCR_cycle_opt,GCR_zero_init_guess,GCR_use_residual_tol,    &
    GCR_adi_add_full_soln,GCR_use_tol_abs,L_gcr_fast_x,          &
    GCR_max_iterations,GCR_diagnostics,GCR_precon_option,        &
    GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,    &
    GCR_tol,GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_tol2,      &
    GCR_ADI_pseudo_timestep,G_term_tol,GCR_its_avg_step

USE lbc_input_mod, ONLY:                                         &
    L_LBC_balance,L_lbc_new,L_transparent,n_rims_to_do,L_lbc_old
    
USE var_input_mod, ONLY:                                         &
    L_regular   
    
USE eg_alpha_ramp_mod, ONLY:                                     &
     alpha_relax_type   
       
USE um_input_control_mod, ONLY: model_domain, l_endgame   
USE domain_params,        ONLY: mt_global   
       
IMPLICIT NONE

! ------------------------------------------------
! Parameters not read in by the run_dyn namelist.
! ------------------------------------------------


! T: reset polar values to mean value every polar_reset_timesteps
LOGICAL, PARAMETER :: L_polar_reset   = .FALSE.
LOGICAL, PARAMETER :: L_interp_depart = .FALSE. 
                                      ! interpolate to find u,v departure pts
LOGICAL, PARAMETER :: L_qwaterload    = .TRUE.  
                                      ! true if using adding waterloading terms
LOGICAL, PARAMETER :: L_fint_theta    = .FALSE. ! true:  fully-interpolating 
                                      ! semi-lagrangian theta advection will be
                                      ! used false: standard non-interpolating 
                                      ! in the vertical

! interval for resetting polar to mean
INTEGER, PARAMETER :: polar_reset_timesteps = 1

REAL, PARAMETER :: tol_sc_fact  = 1.0      ! linear solver scaling factor
REAL, PARAMETER :: T_surf_ref   = 290.0    ! SI Reference surface temperature 
REAL, PARAMETER :: p_surf_ref   = 1.0e5    ! SI Reference surface pressure 
REAL, PARAMETER :: extrp_weight = 1.5      ! time interpolation weight for first
                                           ! SL iteration.
 
! -------------------------------
! Items read in by run_dyn namelist.
! -------------------------------
    
LOGICAL :: L_mix_ratio    = .FALSE. ! use MR code 
LOGICAL :: L_new_tdisc    = .FALSE. ! activate new implicit time scheme
LOGICAL :: L_thmono_fixed = .FALSE. ! for non-interpolating theta - 
                                    ! use corrected code or not

INTEGER :: IntRand_seed = imdi    ! integer seed for random number generation.
INTEGER :: NumCycles    = imdi    !  Parameters for dynamics-physics cycling
INTEGER :: i_ND_solver_vn = imdi  ! dynamics solver version
                                  ! 1 - 2A scheme
                                  ! 2 - 2B scheme

!Time weights for integration scheme
REAL :: alpha_1    = rmdi   
REAL :: alpha_2    = rmdi
REAL :: alpha_3    = rmdi   
REAL :: alpha_4    = rmdi
REAL :: alpha_1_2  = rmdi   
REAL :: alpha_2_2  = rmdi
REAL :: alpha_3_2  = rmdi   
REAL :: alpha_4_2  = rmdi


! ENDGAME only namelist items
LOGICAL :: l_simple_coriolis    =.FALSE. ! option to use a simplified
                                         ! representation of the Coriolis terms
LOGICAL :: L_fix_mass           =.FALSE. ! Enforce global mass conservation by 
                                          ! rescaling density at end of timestep
                                          ! should only be used in global models
LOGICAL :: L_eg_dry_static_adj  =.FALSE. ! Enforces static stability on 
                                          ! input theta field 

LOGICAL :: L_inc_solver         =.FALSE. ! Flag for incremental solver
LOGICAL :: L_conserv_smooth_lap =.FALSE. ! to use the option of redistributing 
                                          ! mass based on only the Laplacian in
                                          ! the mass conservation routine

LOGICAL :: L_accel_convergence  =.FALSE. ! accelerated convergence switches
LOGICAL :: L_init_Fnm1          =.FALSE. ! Fnm1 initialisation
LOGICAL :: l_impl_horz_drag     =.FALSE. ! Implicit horizontal drag switch
LOGICAL :: l_fast_vert_adjust   =.FALSE. ! Switch to use simplified 
                                         ! adjust_vert_bound code
LOGICAL :: l_sl_bc_correction   =.FALSE. ! Flag for using corrected treatment 
                                         ! of upper and lower boundaries in
                                         ! semi-Lagrangian advection scheme.

LOGICAL :: l_filter_dump        =.FALSE. ! Option to filter initial fields 
                                          ! during first timestep. 
                                       ! Usually done when initialising 
                                       ! global model with either ND or EC data.

INTEGER :: eg_vert_damp_profile = imdi    ! Vertical damping flag
REAL :: eta_s                   = rmdi    ! Height (in eta) above
                                          ! which to apply damping
                              
REAL :: Ih                      = rmdi    ! Hydrostatic switch 
                                          !(1=nonhydrostatic; 0=hydrostatic)

LOGICAL :: L_RK_dps             =.TRUE.  ! Runge-Kutta departure point scheme
LOGICAL :: L_eliminate_rho      =.TRUE.  ! Diagnostic rho
LOGICAL :: l_expl_horz_drag     =.TRUE.  ! Explicit horizontal drag switch ??
INTEGER :: InnIts               = 2       ! Number of inner iterations  
REAL :: eg_vert_damp_coeff      = 0.05  

! End of ENDGAME only

! ----------------------
! namelist run_dyn
! ----------------------

      NAMELIST/RUN_Dyn/                                                &
       L_GCR_cycle_opt, L_mix_ratio, L_new_tdisc,  L_thmono_fixed,     & 
       IntRand_seed, L_transparent,                                    &
       n_rims_to_do, i_ND_solver_vn,                                   &
       alpha_1,alpha_2,alpha_3,alpha_4,                                &
       alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
       NumCycles, GCR_use_residual_Tol, GCR_adi_add_full_soln,         &
       GCR_max_iterations,GCR_precon_option,                           &
       GCR_tol, GCR_tol2,GCR_ADI_pseudo_timestep,                      &
! ENDGAME
       L_fix_mass,                                                     &
       alpha_relax_type,                                               &
       eg_vert_damp_profile, eta_s,                                    &
       l_filter_dump 

CONTAINS

SUBROUTINE check_run_dyn()

! Description:
!   Subroutine to apply logic checks and set variables based on the 
!   options selected in the run_dyn namelist.

USE precon_constants_mod, ONLY:                                   &
          vert_plus_xyz_ADI_precon,xyz_ADI_precon,                &
          vert_plus_xz_ADI_precon,xz_ADI_precon
! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_DYN',zhook_in,zhook_handle)

! ---------------------------------------------
! Dynamics logic control - New Dynamics and ENDGAME
! ---------------------------------------------
IF (GCR_use_residual_Tol) THEN
    GCR_use_tol_abs=.FALSE.
    GCR_tol_res = GCR_tol
    GCR_tol_res2 = GCR_tol2
ELSE     
    GCR_use_tol_abs=.TRUE.
    GCR_tol_abs = GCR_tol
    GCR_tol_abs2 = GCR_tol2    
END IF

! ---------------------------------------------
! Dynamics logic control - ENDGame only
! ---------------------------------------------
IF (model_domain /= mt_global .AND. l_endgame) THEN 
    l_fix_mass = .FALSE.  
END IF  
! ---------------------------------------------
! Dynamics logic control - New Dynamics only
! ---------------------------------------------
IF (.NOT. l_endgame) THEN
  IF (GCR_precon_option == vert_plus_xyz_ADI_precon   .OR.      &
       GCR_precon_option == xyz_ADI_precon            .OR.      &
       GCR_precon_option == vert_plus_xz_ADI_precon   .OR.      &
       GCR_precon_option == xz_ADI_precon   ) THEN
    GCR_n_ADI_pseudo_timesteps = 1
  END IF
END IF


IF (lhook) CALL dr_hook('CHECK_RUN_DYN',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_dyn

END MODULE  dynamics_input_mod
