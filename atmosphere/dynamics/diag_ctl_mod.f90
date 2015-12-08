! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Semi-Lagrangian parameters.
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: dynamics

      MODULE diag_ctl_mod

! Description:
!          Diagnostic printing control parameters
!
! Method:
!         The required parameters are read in by NAMELIST
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
 
      IMPLICIT NONE
      
      LOGICAL :: L_print_w
      LOGICAL :: L_print_div
      LOGICAL :: L_diag_print_ops     ! diagnostic prints for ops
      LOGICAL :: L_print_pe     ! print diagnostics on all pe's if true
      LOGICAL :: L_print_shear        ! wind diagnostic prints
      LOGICAL :: L_print_max_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_L2norms  ! l2norm diagnostic prints
      LOGICAL :: L_diag_L2helm   ! l2norm diagnostic prints from solver
      LOGICAL :: L_diag_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_noise    ! w diagnostic prints
      LOGICAL :: L_flush6        ! if T then flush buffers on failure
      LOGICAL :: L_diag_print  ! Print diagnostics
      LOGICAL :: L_print_lapse ! Print lapse_rate diagnostics
      LOGICAL :: L_print_wmax ! Print max w diagnostic
      LOGICAL :: L_print_theta1 ! Print level 1 theta diagnostic

      INTEGER :: Instability_diagnostics  ! >0 if wanted, 0 otherwise
      INTEGER :: print_step    ! To control diagnostic printing interval
      INTEGER :: diag_interval ! diagnostic printing sampling frequency
      INTEGER :: norm_lev_start ! start level for norm diagnostics
      INTEGER :: norm_lev_end   ! end level for norm diagnostics
      INTEGER :: first_norm_print ! first timestep for norm printing
      INTEGER :: dom_w_in    ! define start for block printing
      INTEGER :: dom_e_in    ! define end for block printing
      INTEGER :: dom_s_in    ! define start for block printing
      INTEGER :: dom_n_in    ! define end for block printing
      INTEGER :: blockx_in    ! define size for block printing
      INTEGER :: blocky_in    ! define size for  block printing
      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta

      INTEGER, ALLOCATABLE :: time_w_max(:) ! Timestep of max w
      INTEGER, ALLOCATABLE :: time_div_max(:) ! Timestep of max div
      INTEGER, ALLOCATABLE :: time_div_min(:) ! Timestep of min div
      INTEGER, ALLOCATABLE :: time_lapse_min(:) ! Timestep of min
      INTEGER, ALLOCATABLE :: time_max_shear(:) !Timestep max shear
      INTEGER, ALLOCATABLE :: time_max_wind(:) ! Timestep of max wind
      INTEGER, ALLOCATABLE :: time_KE_max(:) ! Timestep of max KE
      INTEGER, ALLOCATABLE :: time_KE_min(:) ! Timestep of min KE
      INTEGER, ALLOCATABLE :: time_noise_max(:) ! Timestep of max

      REAL :: w_print_limit ! w Threshold for diagnostic printing
      REAL :: w_conv_limit  ! w Threshold for limiting convection
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL, ALLOCATABLE :: max_w_run(:)  ! Max w at a level
      REAL, ALLOCATABLE :: max_div_run(:) ! Max divergence at a level
      REAL, ALLOCATABLE :: min_div_run(:) ! Min divergence at a level
      REAL, ALLOCATABLE :: min_lapse_run(:) ! Min dtheta/dz at a level
      REAL, ALLOCATABLE :: max_shear_run(:) ! Max shear at a level
      REAL, ALLOCATABLE :: max_wind_run(:) ! Max wind at a level
      REAL, ALLOCATABLE :: max_KE_run(:)   ! Max KE at a level
      REAL, ALLOCATABLE :: min_KE_run(:)   ! Min KE at a level
      REAL, ALLOCATABLE :: max_noise_run(:) ! Max noise at a level

      END MODULE diag_ctl_mod
