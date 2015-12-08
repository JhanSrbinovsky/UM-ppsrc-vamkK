! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ PE_Helmholtz solver parameters.

MODULE gcr_input_mod

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE

! Description:
!          PE_Helmholtz parameters constant for run          
!          set by NAMELIST input          
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: dynamics_solver
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

!----------------------------------------------
! the following are not read in by a namelist and fixed
!----------------------------------------------

! Speed up solver convergence when dynamics-physics cycling is used.
! not used in any jobs as FAST is slow in real jobs.
! true use 'faster' non reproducible code
LOGICAL, PARAMETER :: L_gcr_fast_x        = .FALSE. 
LOGICAL, PARAMETER :: GCR_zero_init_guess = .TRUE. ! True if initial guess zero

INTEGER, PARAMETER ::  max_numcycles = 5  ! Max number of dynamics cycles
INTEGER, PARAMETER ::  GCR_Restart_value  = 2  
                                 ! No. of iterations before restarts

! Switch controlling diagnostic output.                               
INTEGER, PARAMETER ::  no_GCR_diags = 0        ! 0 = none
INTEGER, PARAMETER ::  initial_and_final = 1   ! 1 = initial and final residuals
INTEGER, PARAMETER ::  ini_fin_itr= 2          ! 2 = all
INTEGER, PARAMETER ::  iteration_count= 3      ! 3 = iteration count processing

REAL, PARAMETER    ::  G_term_tol  = 0.900     ! tolerance for vertical G term

!----------------------------------------------
! the following are not read in by a namelist but 
! are set dependent upon othe namelist inputs.
!----------------------------------------------
LOGICAL  ::  GCR_use_tol_abs  = .FALSE.  
INTEGER  ::  GCR_n_ADI_pseudo_timesteps = imdi
                                 ! Number of ADI pseudo timesteps to perform.
                                 
INTEGER  ::  GCR_its_switch(max_numcycles) = imdi ! Iterations analysis switch
INTEGER  ::  GCR_max_its(max_numcycles)    = imdi ! Max its this period
INTEGER  ::  GCR_min_its(max_numcycles)    = imdi ! Min its this period
INTEGER  ::  GCR_max_time(max_numcycles)   = imdi ! max timestep number
INTEGER  ::  GCR_min_time(max_numcycles)   = imdi ! min timestep number
INTEGER  ::  GCR_sum_its(max_numcycles)    = imdi ! Sum its over period 
                                 
REAL  ::  GCR_tol_res = rmdi
REAL  ::  GCR_tol_res2= rmdi
REAL  ::  GCR_tol_abs = rmdi
REAL  ::  GCR_tol_abs2= rmdi

!----------------------------------------------
! the following are read in by a namelist:
! either run_dyn or run_dyntest
!----------------------------------------------
LOGICAL  :: L_GCR_cycle_opt       = .FALSE.  ! Speed up solver convergence when 
                                             !  dynamics-physics cycling
LOGICAL  :: GCR_use_residual_tol  = .FALSE. 
LOGICAL  :: GCR_adi_add_full_soln = .FALSE.  ! true then use full equation on
                                             ! RHS on second and subsequent 
                                             ! ADI timesteps
                                                                    
INTEGER  ::  GCR_max_iterations   = imdi
INTEGER  ::  GCR_precon_option    = imdi
INTEGER  ::  GCR_Diagnostics      = initial_and_final  ! Initial and final norms
INTEGER  ::  GCR_its_avg_step(3)  = imdi ! Iterations analysis step now

REAL  ::  GCR_tol                 = rmdi
REAL  ::  GCR_tol2                = rmdi
REAL  ::  GCR_ADI_pseudo_timestep = rmdi


END MODULE gcr_input_mod
