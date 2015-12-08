! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for SL dynamics.

! Description:
!   Module containing input switches/settings as used by the
!   semi_lagrangian code and the check_run_sl routine
!   for logic checking the selected settings.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

! Method:
!   Switches are initialised to false and read in from the
!   UMUI. The module may then be used directly where the switches
!   are needed within the dynamics semi_lagrangian code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE sl_input_mod

USE missing_data_mod, ONLY: rmdi, imdi

USE var_look_mod, ONLY:                    &
    max_look, halo_lam ,halo_phi,          &
    recip_dlam,recip_dphi,                 &
    look_lam,look_phi

USE comms_2d_mod, ONLY: L_2dcomm

IMPLICIT NONE

! -------------------------------
! Parameters not read in by namelist.
! -------------------------------

INTEGER, PARAMETER :: Number_SL_choices = 4
INTEGER, PARAMETER :: Theta_SL = 1
INTEGER, PARAMETER :: moist_SL = 2
INTEGER, PARAMETER :: Wind_SL = 3
INTEGER, PARAMETER :: rho_SL = 4

INTEGER, PARAMETER :: unselected = 0
INTEGER, PARAMETER :: standard = 1
INTEGER, PARAMETER :: accurate = 2


! -------------------------------
! Items not read in by namelist.
! -------------------------------

LOGICAL :: L_conserv(Number_SL_choices) = .FALSE.  
                                    ! T if conservation to be enforced.
LOGICAL :: L_moist_nonhydro_conserve    = .FALSE.    
                                    ! improved nonhydrostatic conservation  
                                    ! of q, qcl and qcf

LOGICAL :: L_mono(Number_SL_choices) = .FALSE.    
                                    ! T if monotonicity to be enforced.
LOGICAL :: L_high(Number_SL_choices) = .FALSE.     
                                    ! T if high order interpolation scheme
                                          ! to be used.
LOGICAL :: L_Ritchie_high = .FALSE. ! T if high order scheme to be
                                    ! used in Ritchie routine.
LOGICAL :: L_Ritchie_mono = .FALSE. ! T if monotone scheme to be used
                                    ! in Ritchie routine.

INTEGER :: thmono_levels = 0   ! levels for monotone fully-interp theta

INTEGER :: Depart_scheme = 1   ! which departure point scheme to use (Ritchie)

!-----------------------------------------------------------
! run_sl namelist items
! ----------------------------------------------------------

LOGICAL :: L_conserve_tracers = .FALSE. ! Run tracer advection code with 
                                        ! conservation on 
LOGICAL :: L_2d_sl_geometry = .FALSE.   ! T then only perform vector 
                                        ! co-ordinate geometry in 2d
LOGICAL :: L_sl_halo_reprod = .FALSE.   ! T: sl code bit repoducible with 
                                        ! any sensible halo size

! select what type of moisture conservation is required.
INTEGER :: moisture_conservation = imdi

! a code saying which high monotone scheme to use.
INTEGER :: high_order_scheme(Number_SL_choices) = imdi

! a code saying which  monotone scheme to use.
INTEGER :: monotone_scheme(Number_SL_choices) = imdi

INTEGER :: Ritchie_high_order_scheme = imdi ! which high order scheme to use
INTEGER :: Ritchie_monotone_scheme   = imdi ! which monotone scheme to use
INTEGER :: Depart_order              = imdi ! how many iterations/term 
                                            ! to use for Ritchie scheme
INTEGER :: interp_vertical_search_tol= imdi ! number of levels either side 
                                            ! of default level to search.

INTEGER :: Instability_diagnostics   = imdi ! >0 if wanted, 0 otherwise

REAL :: thmono_height = rmdi             ! top for monotone fully-interp theta

! ----------------------
! Namelist
! ----------------------

      NAMELIST/RUN_SL/ thmono_height, L_2d_sl_geometry,                &
       L_sl_halo_reprod, L_2dcomm,                                     &
       high_order_scheme, monotone_scheme,                             &
       Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
       Depart_order,interp_vertical_search_tol,                        &
       L_conserve_tracers,                                             &
       Instability_diagnostics,moisture_conservation

CONTAINS

SUBROUTINE check_run_sl()

! Description:
!   Subroutine to apply logic checks and setup variables based on the 
!   options selected in the run_sl namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_SL',zhook_in,zhook_handle)

!-------------------------------------------------------
! SL input setup
!-------------------------------------------------------
! moisture conservation setup
!-------------------------------------------------------
IF (moisture_conservation==accurate) THEN
  L_conserv(moist_SL) = .TRUE.
  L_moist_nonhydro_conserve    = .TRUE.
ELSE IF (moisture_conservation==standard) THEN
  L_conserv(moist_SL) = .TRUE.
  L_moist_nonhydro_conserve    = .FALSE.
END IF
!-------------------------------------------------------
! monotone scheme setup
!-------------------------------------------------------      
IF (monotone_scheme(Theta_SL)/=unselected) THEN
  L_mono(Theta_SL) = .TRUE.
END IF
IF (monotone_scheme(moist_SL)/=unselected) THEN
  L_mono(moist_SL) = .TRUE.
END IF
IF (monotone_scheme(Wind_SL)/=unselected) THEN   
  L_mono(Wind_SL) = .TRUE.
END IF
IF (monotone_scheme(rho_SL)/=unselected) THEN
  L_mono(rho_SL) = .TRUE.         
END IF
!-------------------------------------------------------
! high order scheme setup
!-------------------------------------------------------      
IF (high_order_scheme(Theta_SL)/=unselected) THEN
  L_high(Theta_SL) = .TRUE.
END IF
IF (high_order_scheme(moist_SL)/=unselected) THEN
  L_high(moist_SL) = .TRUE.
END IF
IF (high_order_scheme(Wind_SL)/=unselected) THEN   
  L_high(Wind_SL) = .TRUE.
END IF
IF (high_order_scheme(rho_SL)/=unselected) THEN
  L_high(rho_SL) = .TRUE.
END IF
!-------------------------------------------------------
! Ritchie scheme setup
!-------------------------------------------------------      
IF (Ritchie_high_order_scheme/=unselected) THEN
  L_Ritchie_high = .TRUE.
END IF
IF (Ritchie_monotone_scheme/=unselected) THEN
  L_Ritchie_mono = .TRUE.
END IF

IF (lhook) CALL dr_hook('CHECK_RUN_SL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_sl

END MODULE  sl_input_mod
