! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Load a variable from files, interpolate in time & height
!  and return on model levels & timestep

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_NUDGING_CONTROL.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_varloader_pressure_lev( &
  varname,                                 & ! Variable name
  global_row_length,                       & ! Row length
  global_rows,                             & ! Rows
  file_levels,                             & ! Levels in file
  proc_row_length_min,                     & ! Min. column in the PE
  proc_row_length_max,                     & ! Max. column in the PE
  proc_rows_min,                           & ! Min. row in the PE
  proc_rows_max,                           & ! Max. row in the PE
  base_lambda,                             & ! minimum longitude in model
  delta_lambda,                            & ! longitude increment in model
  base_phi,                                & ! minimum latitude in model
  delta_phi,                               & ! latitude increment in model
  base_lambda_ana,                         & ! minimum longitude in data
  delta_lambda_ana,                        & ! longitude increment in data
  base_phi_ana,                            & ! minimum latitude in data
  delta_phi_ana,                           & ! latitude increment in data
  variable_file_data,                      & ! variable at first timestep
  variable_file_data2,                     & ! variable at second timestep
  frac,                                    & ! fraction between 1 & 2
  model_levels,                            & ! model levels
  pressure_model_levels,                   & ! model pressure
  pressure_surf_level,                     & ! surface pressure
  variable_model_levels,                   & ! variable on model levels
  debug)                                     ! Debug flag

USE nudging_control               ! Use standard nudging switches
USE nudging_ecmwf_60level_def     ! Use ECMWF model elvel information
USE nudging_jra_plevel_def        ! Use JRA pressure level
USE atmos_constants_mod, ONLY: kappa,pref   ! Atmospheric constants

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)  :: varname              ! Variable name
INTEGER, INTENT(IN)       :: global_row_length    ! Row length
INTEGER, INTENT(IN)       :: global_rows ! Rows
INTEGER, INTENT(IN)       :: file_levels          ! Levels in file
INTEGER, INTENT(IN)       :: proc_row_length_min  ! Min. column in the PE
INTEGER, INTENT(IN)       :: proc_row_length_max  ! Max. column in the PE
INTEGER, INTENT(IN)       :: proc_rows_min        ! Min. row in the PE
INTEGER, INTENT(IN)       :: proc_rows_max        ! Max. row in the PE
INTEGER, INTENT(IN)       :: frac                 ! Frac between start & end
INTEGER, INTENT(IN)       :: model_levels         ! Model levels
REAL, INTENT(IN)          :: base_lambda          ! MInimum longitude in model
REAL, INTENT(IN)          :: delta_lambda         ! longitude increment in model
REAL, INTENT(IN)          :: base_phi             ! Minimum latitude in model
REAL, INTENT(IN)          :: delta_phi            ! Latitude increment in model
REAL, INTENT(IN)          :: base_lambda_ana      ! Minimum longitude in data
REAL, INTENT(IN)          :: delta_lambda_ana     ! Longitude increment in data
REAL, INTENT(IN)          :: base_phi_ana         ! Minimum latitude in data
REAL, INTENT(IN)          :: delta_phi_ana        ! Latitude increment in data

! Temporary Arrays used to transfer data from netcdf files to interpolation
REAL, INTENT(IN)          :: variable_file_data                        &
  (global_row_length, global_rows, file_levels)

REAL, INTENT(IN)          :: variable_file_data2                       &
  (global_row_length, global_rows, file_levels)

! Model pressure on model levels
REAL, INTENT(IN) :: pressure_model_levels(                             &
  proc_row_length_min:proc_row_length_max,                             &
  proc_rows_min:proc_rows_max, model_levels)

! Surface pressure from the data
REAL, INTENT(IN) :: pressure_surf_level                                &
 (global_row_length, global_rows)

! Output variable interpolated in time and on to model levels
REAL, INTENT(OUT) :: variable_model_levels(                            &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

INTEGER, INTENT(IN) :: debug               ! Debug flag

! Variable interpolated in time, but not in height
REAL :: variable_file_levels                                           &
  (global_row_length, global_rows, file_levels)

! Used in flipping analyses from top down to bottom up
! Legacy of sawitching from data obtained at the BADC to data from BDAN
REAL :: variable_file_levels_temp                                     &
  (global_row_length, global_rows, file_levels)

! Variable on ECMWF hybrid p-levels and model grid
REAL             :: variable_interp(                                  &
  proc_row_length_min:proc_row_length_max,                            &
  proc_rows_min:proc_rows_max, 1:file_levels)

! Pressure on ECMWF hybrid p-levels and anlysis grid
REAL             :: pressure_ecmwf_levels                              &
  (global_row_length, global_rows, file_levels)

! Pressure on ECMWF hybrid p-levels and model grid
REAL             :: pressure_interp(                                   &
  proc_row_length_min:proc_row_length_max,                             &
  proc_rows_min:proc_rows_max, 1:file_levels)

INTEGER :: i, j, k, l, a, b         ! Loop variables
INTEGER :: dummy

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_VARLOADER_PRESSURE_LEV',zhook_in,zhook_handle)

! End of Header
!******************************************************************************

! Standard Subroutine Entry Comment
IF(debug > 10) WRITE(OUT, *)                                           &
 ': NUDGING_VARLOADER_PRESSURE_LEV: Entering Routine'

! Initialise data to zero
variable_file_levels = 0

!Drop a line to let us know where we are
IF(debug > 10) WRITE(OUT, *)                                          &
 ': NUDGING_VARLOADER_PRESSURE_LEV: About to Interpolate'

! Interpolate between the two input data files onto the model timestep
! DEPENDS ON: nudging_interpolation_three_d
CALL nudging_interpolation_three_d(  &
  global_row_length,                 & ! Row length
  global_rows,                       & ! Rows
  file_levels,                       & ! File levels
  variable_file_data(:,:,:),         & ! Variable at start
  variable_file_data2(:,:,:),        & ! Variable at end
  frac,                              & ! Frac between start & end
  variable_file_levels,              & ! Interpolated variable
  0,                                 & ! Linear interpolation
  debug)                               ! Debug flag

! If we know the data is in top down format flip it
! This is a legacy of BADC data which is bottom up
! Converted to BDAN data which is top down, so have to flip
IF(data_source == 0) THEN
  DO k=1, file_levels
    l = (file_levels+1) - k
    variable_file_levels_temp(:,:,l) = variable_file_levels(:,:,k)
  END DO
  DO k=1, file_levels
    variable_file_levels(:,:,k) = variable_file_levels_temp(:,:,k)
  END DO
END IF

! Depending on level structure of analyses source choose a method
! for converting from analyses levels to model levels
SELECT CASE(ndg_analysis_source)

! ECMWF hybrid pressure levels
  CASE(0)

! Calculate levels on basis
    SELECT CASE(file_levels)

      CASE(60)

! Calculate the ECMWF hybrid pressure levels
! DEPENDS ON: nudging_calc_ecmwf_60levs
        CALL nudging_calc_ecmwf_60levs(   &
          global_row_length,              & ! Length of model row
          global_rows,                    & ! Number of model rows
          file_levels,                    & ! Number of file levels
          pressure_surf_level,            & ! Surface pressure
          pressure_ecmwf_levels,          & ! Return hybrid p levels
          no_debug)                         ! Debug flag

      CASE(91)

! Calculate the ECMWF hybrid pressure levels
! DEPENDS ON: nudging_calc_ecmwf_91levs
        CALL nudging_calc_ecmwf_91levs(  &
          global_row_length,             & ! Length of model row
          global_rows,                   & ! Number of model rows
          file_levels,                   & ! Number of file levels
          pressure_surf_level,           & ! Surface pressure
          pressure_ecmwf_levels,         & ! Return hybrid p levels
          no_debug)                        ! Debug flag

      CASE DEFAULT
        nmessage='Unknown number of levels in ECMWF Data'
        dummy = 999

        CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                  &
                     dummy,nmessage)
    END SELECT     ! file levels

! Convert temperature to potential temperature
    IF(varname == temp_name)                                             &
      variable_file_levels = variable_file_levels*                       &
                            (pref/pressure_ecmwf_levels)**kappa

! Drop a line to let us know where we are
    IF(debug > 10) WRITE(OUT, *)                                         &
      ': NUDGING_VARLOADER_PRESSURE_LEV: Fin Convert Theta,',            &
      'typical "',varname,'" value equals ',                             &
      variable_file_levels(dbg_x,dbg_y,dbg_z),                           &
      pref, pressure_ecmwf_levels(dbg_x,dbg_y,dbg_z), kappa

! Convert from analysis to model grid
!DEPENDS ON: NUDGING_ANA2MOD_GRID
     CALL nudging_ana2mod_grid(       &
       global_row_length,             &  ! Global length of rows (in data)
       global_rows,                   &  ! Global numver of rows (in data)
       file_levels,                   &  ! Number of levels in data
       proc_row_length_min,           &  ! min column in model pe
       proc_row_length_max,           &  ! max column in model pe
       proc_rows_min,                 &  ! min row in model pe
       proc_rows_max,                 &  ! max row in model pe
       variable_file_levels,          &  ! Variable on data grid
       variable_interp,               &  ! Variable on model grid
       base_lambda,                   &  ! Min longitude in model
       delta_lambda,                  &  ! Longitude increment in model
       base_phi,                      &  ! Min latitude in model
       delta_phi,                     &  ! Latitude increment in model
       base_lambda_ana,               &  ! Min longitude in data
       delta_lambda_ana,              &  ! Longitude increment in data
       base_phi_ana,                  &  ! Min latitude in data
       delta_phi_ana,                 &  ! Latitude increment in data
       debug)                            ! Debug flag

!DEPENDS ON: NUDGING_ANA2MOD_GRID
     CALL nudging_ana2mod_grid(       &
       global_row_length,             &  ! Global length of rows (in data)
       global_rows,                   &  ! Global numver of rows (in data)
       file_levels,                   &  ! Number of levels in data
       proc_row_length_min,           &  ! min column in model pe
       proc_row_length_max,           &  ! max column in model pe
       proc_rows_min,                 &  ! min row in model pe
       proc_rows_max,                 &  ! max row in model pe
       pressure_ecmwf_levels,         &  ! Pressure on data grid
       pressure_interp,               &  ! Pressure on model grid
       base_lambda,                   &  ! Min longitude in model
       delta_lambda,                  &  ! Longitude increment in model
       base_phi,                      &  ! Min latitude in model
       delta_phi,                     &  ! Latitude increment in model
       base_lambda_ana,               &  ! Min longitude in data
       delta_lambda_ana,              &  ! Longitude increment in data
       base_phi_ana,                  &  ! Min latitude in data
       delta_phi_ana,                 &  ! Latitude increment in data
       debug)                            ! Debug flag

! Interpolate the analyses ftom ECMWF hybrid pressure levels
! to UM hybrid height levels
!DEPENDS ON: NUDGING_ECMWF_TO_MOD
     CALL nudging_ecmwf_to_mod(       &
       global_row_length,             &  ! Row length
       global_rows,                   &  ! Rows
       model_levels,                  &  ! ECMWF hybrid levels
       file_levels,                   &  ! Model levels
       proc_row_length_min,           &  ! Min. column in the PE
       proc_row_length_max,           &  ! Max. column in the PE
       proc_rows_min,                 &  ! Min. row in the PE
       proc_rows_max,                 &  ! Max. row in the PE
       pressure_model_levels,         &  ! Pressure on model levels
       pressure_interp,               &  ! pressure on ECMWF hybrid levels
       variable_interp,               &  ! Variable on ecmwf levels
       variable_model_levels,         &  ! Variable on  model levels
       debug)                            ! Miscellanea

!variable_model_levels = pressure_interp

! ECMWF fixed pressure levels
  CASE(1)

    IF(file_levels /= ecmwf_nplevs) THEN
      nmessage = 'Number of levels in file differ from standard'
      dummy = 999

      CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV', dummy, nmessage)
    END IF

! Convert temperature to potential temperature
    IF(varname == temp_name) THEN
      DO k=1, file_levels
        variable_file_levels(:,:,k) =   variable_file_levels(:,:,k)*     &
                                      (pref/data_pressure(k))**kappa
      END DO
    END IF

! Interpolate analyses onto UM grid

! Drop a line to let us know where we are
    IF(debug > 10) WRITE(OUT, *)                                        &
    ': NUDGING_VARLOADER_PRESSURE_LEV: Finished Interpolation,',        &
    'Typical "t" value equals ',variable_file_levels(dbg_x,dbg_y,dbg_z)

! Interpolate analyses from ECMWF fixed p-levels to model levels
! DEPENDS ON: nudging_pres_to_mod
    CALL nudging_pres_to_mod(        &
      global_row_length,             & ! Row length
      global_rows,                   & ! Rows
      model_levels,                  & ! ECMWF fixed p levels
      file_levels,                   & ! model levels
      proc_row_length_min,           & ! Min. column in the PE
      proc_row_length_max,           & ! Max. column in the PE
      proc_rows_min,                 & ! Min. row in the PE
      proc_rows_max,                 & ! Max. row in the PE
      pressure_model_levels,         & ! Pressure on model levels
      data_pressure,                 & ! ECMWF fixed P levels
      variable_file_levels,          & ! Variable on ecmwf levels
      variable_model_levels,         & ! Variable on model levels
      debug)                           ! Debug flag

! If variable is on UM hybrid height levels then trivial
  CASE(2)

! Convert from analysis to model grid
!DEPENDS ON: NUDGING_ANA2MOD_GRID
    CALL nudging_ana2mod_grid(       &
     global_row_length,             &  ! Global length of rows (in data)
     global_rows,                   &  ! Global numver of rows (in data)
     file_levels,                   &  ! Number of levels in data
     proc_row_length_min,           &  ! min column in model pe
     proc_row_length_max,           &  ! max column in model pe
     proc_rows_min,                 &  ! min row in model pe
     proc_rows_max,                 &  ! max row in model pe
     variable_file_levels,          &  ! Variable on data grid
     variable_interp,               &  ! Variable on model grid
     base_lambda,                   &  ! Min longitude in model
     delta_lambda,                  &  ! Longitude increment in model
     base_phi,                      &  ! Min latitude in model
     delta_phi,                     &  ! Latitude increment in model
     base_lambda_ana,               &  ! Min longitude in data
     delta_lambda_ana,              &  ! Longitude increment in data
     base_phi_ana,                  &  ! Min latitude in data
     delta_phi_ana,                 &  ! Latitude increment in data
     debug)                            ! Debug flag


    variable_model_levels(:,:,:) = 0.

! As on same levels no interpolation  required
    DO j= proc_rows_min, proc_rows_max
      DO i= proc_row_length_min, proc_row_length_max
        variable_model_levels(i,j, ndg_lev_bottom:ndg_lev_top) =          &
              variable_interp(i,j, ndg_lev_bottom:ndg_lev_top)
      END DO
    END DO

! Convert temperature to potential temperature
    IF(varname == temp_name) THEN
      DO k=ndg_lev_bottom, ndg_lev_top
        DO j = proc_rows_min, proc_rows_max
          DO i = proc_row_length_min, proc_row_length_max
            variable_model_levels(i,j,k) =                                   &
             variable_model_levels(i,j,k)*                                   &
             (pref/pressure_model_levels(i,j,k))**kappa
          END DO
        END DO
      END DO
    END IF

! Variable on JRA fixed Pressure levels
  CASE(3)

    IF(file_levels /= jra_nplevs) THEN
      nmessage = 'Number of levels in file differ from standard'
      dummy = 999

      CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                     &
                 dummy, nmessage)
    END IF

! Convert temperature to potential temperature
    IF(varname == temp_name) THEN
      DO k=1, file_levels
        variable_file_levels(:,:,k) =                                    &
         variable_file_levels(:,:,k)*                                    &
         (pref/jra_pressure(k))**kappa
      END DO
    END IF

! Interpolate analyses onto UM grid

! Drop a line to let us know where we are
    IF(debug > 10) WRITE(OUT, *)                                        &
     ': NUDGING_VARLOADER_PRESSURE_LEV: Finished Interpolation,',       &
     'Typical "t" value equals ',variable_file_levels(dbg_x,dbg_y,dbg_z)

! Interpolate analyses from JRA fixed P levels to model levels
! DEPENDS ON: nudging_pres_to_mod
    CALL nudging_pres_to_mod(          &
      global_row_length,               & ! Row length
      global_rows,                     & ! Rows
      model_levels,                    & ! JRA levels
      file_levels,                     & ! model levels
      proc_row_length_min,             & ! Min. column in the PE
      proc_row_length_max,             & ! Max. column in the PE
      proc_rows_min,                   & ! Min. row in the PE
      proc_rows_max,                   & ! Max. row in the PE
      pressure_model_levels,           & ! Pressure on model levels
      jra_pressure,                    & ! JRA pressure levels
      variable_file_levels,            & ! Variable on JRA levels
      variable_model_levels,           & ! Variable on model levels
      debug)                             ! Debug flag

  CASE DEFAULT

    nmessage = 'Unkown Analysis data source'
    dummy = 999

    CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                       &
                 dummy, nmessage)

END SELECT         ! Pressure level scheme

! Standard Subroutine Exit Comment
IF(debug > 10) WRITE(OUT, *)                                           &
 ': NUDGING_VARLOADER_PRESSURE_LEV: Leaving Routine'

IF (lhook) CALL dr_hook('NUDGING_VARLOADER_PRESSURE_LEV',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_varloader_pressure_lev

