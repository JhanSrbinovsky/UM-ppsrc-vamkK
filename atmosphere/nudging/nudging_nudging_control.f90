! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Main control routine for the nudging

!  Loads variable, interpolates it onto model levels
!  Performs the nudging and calculates diagnostics

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_nudging_control(  &
  varname,                          &     ! Variable name
  grid_type,                        &     ! Grid type
  field_type,                       &     ! Field type
  global_row_length,                &     ! (Global) row length
  global_rows,                      &     ! (Global) rows
  file_levels,                      &     ! File levels
  proc_row_length_min,              &     ! Min. column
  proc_row_length_max,              &     ! Max. column
  proc_rows_min,                    &     ! Min. row
  proc_rows_max,                    &     ! Max. row
  model_timestep,                   &     ! Model timestep size
  sin_theta_latitude_min,           &     ! sine of min latitude
  sin_theta_latitude_max,           &     ! sine of max latitude
  base_lambda,                      &     ! minimum longitude in model
  delta_lambda,                     &     ! longitude increment in model
  base_phi,                         &     ! minimum latitude in model
  delta_phi,                        &     ! latitude increment in model
  base_lambda_ana,                  &     ! minimum longitude in data
  delta_lambda_ana,                 &     ! longitude increment in data
  base_phi_ana,                     &     ! minimum latitude in data
  delta_phi_ana,                    &     ! latitude increment in data
  variable_file_data,               &     ! data on first timestep
  variable_file_data2,              &     ! data on second timestep
  logp_surf_data,                   &     ! pressure on first timestep
  logp_surf_data2,                  &     ! pressure on second timestep
  frac,                             &     ! Fraction between 1 & 2
  model_levels,                     &     ! Model levels
  model_pressure_model_levels,      &     ! Model pressure
  model_trop_pressure,              &     ! Tropopause Pressure
  model_variable_model_levels,      &     ! Variable nudged
  diag_var_data,                    &     ! Diag. 1 (analyses value)
  diag_var_model,                   &     ! Diag. 2 (model value)
  diag_tend_nudging,                &     ! Diag. 3 (nudging tendency)
  diag_tend_model,                  &     ! Diag. 4 (model tendency)
  diag_relax,                       &     ! Diag. 5 (relaxation paremeter)
  debug                             &     ! Debug flag
)

USE nudging_control                       ! use standard nudging switches

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)  :: varname                 ! Variable name
INTEGER, INTENT(IN)       :: grid_type               ! Grid type
LOGICAL, INTENT(IN)       :: field_type              ! Field type
INTEGER, INTENT(IN)       :: global_row_length       ! (global) row length
INTEGER, INTENT(IN)       :: global_rows             ! (Global) rows
INTEGER, INTENT(IN)       :: file_levels             ! File levels
INTEGER, INTENT(IN)       :: proc_row_length_min     ! min column in grid
INTEGER, INTENT(IN)       :: proc_row_length_max     ! max column in grid
INTEGER, INTENT(IN)       :: proc_rows_min           ! min row in grid
INTEGER, INTENT(IN)       :: proc_rows_max           ! max row in grid
REAL, INTENT(IN)          :: model_timestep          ! size of model timestep
REAL, INTENT(IN)          :: sin_theta_latitude_min  ! sin min lat (T grid)
REAL, INTENT(IN)          :: sin_theta_latitude_max  ! sin max lat (T grid)
REAL, INTENT(IN)          :: frac                    ! Fraction between 1 & 2
INTEGER, INTENT(IN)       :: model_levels            ! number of model levels
REAL, INTENT(IN)  :: base_lambda             ! minimum longitude in the model
REAL, INTENT(IN)  :: delta_lambda            ! longitude increment in model
REAL, INTENT(IN)  :: base_phi                ! minimum latitude in model
REAL, INTENT(IN)  :: delta_phi               ! latitude increment in model
REAL, INTENT(IN)  :: base_lambda_ana         ! base longitude in the data
REAL, INTENT(IN)  :: delta_lambda_ana        ! longitude increment in the data
REAL, INTENT(IN)  :: base_phi_ana            ! base latitude in the data
REAL, INTENT(IN)  :: delta_phi_ana           ! latitude increment in the data

REAL, INTENT(IN)          :: variable_file_data                       &
  (global_row_length, global_rows, file_levels)

REAL, INTENT(IN)          :: variable_file_data2                      &
  (global_row_length, global_rows, file_levels)

! Analyses log surface pressure
REAL, INTENT(IN)        :: logp_surf_data                             &
  (global_row_length, global_rows)

! Analyses log surface pressure
REAL, INTENT(IN)        :: logp_surf_data2                            &
  (global_row_length, global_rows)

! Model pressure
REAL, INTENT(IN) :: model_pressure_model_levels(                      &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! tropopause pressure 
REAL, INTENT(IN)  :: model_trop_pressure(                              &   
 proc_row_length_min:proc_row_length_max,                              &   
 proc_rows_min:proc_rows_max)  

! Nudged variable
REAL, INTENT(INOUT) :: model_variable_model_levels(                   &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! Nudging Diagnostics (1-5)
 REAL, INTENT(OUT) :: diag_var_data(                                   &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

 REAL, INTENT(INOUT) :: diag_var_model(                                &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

 REAL, INTENT(OUT) :: diag_tend_nudging(                               &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

 REAL, INTENT(OUT) :: diag_tend_model(                                 &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

 REAL, INTENT(OUT) :: diag_relax(                                      &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

 INTEGER, INTENT(IN) :: debug          ! Debug flag

!--------------------
! Other variables (non IO)

! Analyses on model levels
REAL    :: data_variable_model_levels (                                &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

! Model variable before nudging
REAL    :: last_variable_model_levels (                                &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

INTEGER :: tropopause_level(                                           & 
 proc_row_length_min:proc_row_length_max,                              &      
 proc_rows_min:proc_rows_max)                     ! Relaxation Parameter 

! Relaxation Parameter
REAL    :: relax_par (                                                 &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max, 1:model_levels )

! Analyses surface Pressure
REAL        :: data_logp_surf_level                                    &
  (global_row_length, global_rows)

! Analyses surface Pressure
REAL        :: data_pressure_surf_level                                &
  (global_row_length, global_rows)

!*******************************************************************
! End of Header

! Standard subroutine entry comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
   ': NUDGING_NUDGING_CTL: Entering routine'
END IF

! If we are using hybrid pressure levels load surface pressure
IF(ndg_analysis_source == 0) THEN

! Interpolate log surface pressure onto timestep
! DEPENDS ON: nudging_interpolation_two_d
  CALL nudging_interpolation_two_d(        &
    global_row_length,                     &  ! Row length
    global_rows,                           &  ! Rows
    logp_surf_data,                        &  ! Variable at start
    logp_surf_data2,                       &  ! Variable at end
    frac,                                  &  ! Frac. between start & end
    data_logp_surf_level,                  &  ! Interpolated variable
    0,                                     &  ! Linear Interpolation
    no_debug)                                 ! Debug flag

! convert lnP to P
  data_pressure_surf_level=EXP(data_logp_surf_level)

END IF

! If in debug mode write out typical relaxation parameter
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
   ': NUDGING_NUDGING_CTL: About to load on pres lev',                 &
   ' for "', varname, '" is: ',                                        &
   data_pressure_surf_level(dbg_x,dbg_y)
END IF

!  Load data interpolated on to this model step
! DEPENDS ON: nudging_varloader_pressure_lev
CALL nudging_varloader_pressure_lev( &
  varname,                           &    ! Variable name
  global_row_length,                 &    ! (global) Row length
  global_rows,                       &    ! (Global) rows
  file_levels,                       &    ! File levels
  proc_row_length_min,               &    ! Min. column in the PE
  proc_row_length_max,               &    ! Max. column in the PE
  proc_rows_min,                     &    ! Min. row in the PE
  proc_rows_max,                     &    ! Max. row in the PE
  base_lambda,                       &    ! minimum longitude in model
  delta_lambda,                      &    ! longitude increment in model
  base_phi,                          &    ! minimum latitude in model
  delta_phi,                         &    ! latitude increment in model
  base_lambda_ana,                   &    ! minimum longitude in data
  delta_lambda_ana,                  &    ! longitude increment in data
  base_phi_ana,                      &    ! minimum latitude in data
  delta_phi_ana,                     &    ! latitude increment in data
  variable_file_data,                &    ! Variable on first time step
  variable_file_data2,               &    ! Variable on second timestep
  frac,                              &    ! Fraction between  1 & 2
  model_levels,                      &    ! Model levels
  model_pressure_model_levels,       &    ! Model pressure
  data_pressure_surf_level,          &    ! Surface pressure (in analyses)
  data_variable_model_levels,        &    ! Variable on model levels
  debug)                                  ! Debug flag

! If in debug mode write out typical interpolated value
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
   ': NUDGING_NUDGING_CTL: Typical interpolated value',                &
   ' for "', varname, '" is: ',                                        &
   data_variable_model_levels                                          &
     (dbg_x,dbg_y,dbg_z)
END IF

! Calculate tropopause level from tropopause pressure 
! DEPENDS ON: nudging_calc_tropopause_level      
CALL nudging_calc_tropopause_level( & 
  proc_row_length_min,              & 
  proc_row_length_max,              & 
  proc_rows_min,                    & 
  proc_rows_max,                    & 
  model_levels,                     & 
  model_pressure_model_levels,      & 
  model_trop_pressure,              & 
  tropopause_level,                 & 
  no_debug) 

! Load the relaxation parameter for variable
! DEPENDS ON: nudging_call_relax
CALL nudging_call_relax(            &
  proc_row_length_min,              & 
  proc_row_length_max,              & 
  proc_rows_min,                    &
  proc_rows_max,                    &
  model_levels,                     &      ! Model levels
  varname,                          &      ! Variable name
  model_timestep,                   &      ! Timestep size
  tropopause_level,                 &      ! Tropopause level
  relax_par,                        &      ! Relaxation parameter
  no_debug)                                ! Debug flag

! If in debug mode write out typical relaxation parameter
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
   ': NUDGING_NUDGING_CTL: Typical relaxation parameter',              &
   ' for "', varname, '" is: ',                                        &
   relax_par(dbg_x,dbg_y,dbg_z), relax_par(dbg_x,dbg_y,dbg_z+5)
END IF

IF(debug > 5) THEN
  WRITE(OUT, *)                                                        &
   ': NUDGING_NUDGING_CTL: Typical UNUDGED value',                     &
   ' for "', varname, '" is: ',                                        &
   model_variable_model_levels                                         &
  (dbg_x,dbg_y,dbg_z)
END IF

! set the pre nudged model variable for diagnostic calculations
last_variable_model_levels = model_variable_model_levels

! Call function to nudge model towards data
! DEPENDS ON: nudging_nudging_three_d
CALL nudging_nudging_three_d(     &
  model_levels,                   &     ! Model levels
  proc_row_length_min,            &     ! Min. column in the PE
  proc_row_length_max,            &     ! Max. column in the PE
  proc_rows_min,                  &     ! Min. row in the PE
  proc_rows_max,                  &     ! Max. row in the PE
  model_variable_model_levels,    &     ! Model variable
  data_variable_model_levels,     &     ! analysis variable
  relax_par,                      &     ! size of timestep
  no_debug)                             ! Debug flag

IF(debug > 5) THEN
  WRITE(OUT, *)                                                        &
    ': NUDGING_NUDGING_CTL: Typical NUDGED value',                     &
    ' for "', varname, '" is: ',                                       &
    model_variable_model_levels                                        &
     (dbg_x,dbg_y,dbg_z)
END IF

! Set nudging diagnostics
! DEPENDS ON: nudging_calc_diags
CALL nudging_calc_diags(       &
  model_levels,                &          ! Model levels
  proc_row_length_min,         &          ! Min. column in the PE
  proc_row_length_max,         &          ! Max. column in the PE
  proc_rows_min,               &          ! Min. row in the PE
  proc_rows_max,               &          ! Max. row in the PE
  last_variable_model_levels,  &          ! Model variable on last timestep
  model_variable_model_levels, &          ! Model variable
  data_variable_model_levels,  &          ! Analyses variable
  relax_par,                   &          ! Relaxation parameter (2D)
  diag_var_data,               &          ! Diag 1 (variable in analyses)
  diag_var_model,              &          ! Diag 2 (variable in model)
  diag_tend_nudging,           &          ! Diag 3 (nudging tend)
  diag_tend_model,             &          ! Diag 4 (model tend)
  diag_relax,                  &          ! Diag 5 (relax par (3D))
  no_debug)                               ! Debug flag

! Standard subroutine exit comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
   ': NUDGING_NUDGING_CTL: Leaving routine'
END IF

RETURN
END SUBROUTINE nudging_nudging_control

