! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Calculates diagnostics for the nudging

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_NUDGING_CONTROL.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_calc_diags(        &
  model_levels,                       & ! NUmber of levels in the Data
  proc_row_length_min,                & ! Minimum column in the PE
  proc_row_length_max,                & ! Maximum column in the PE
  proc_rows_min,                      & ! Minimum row in the PE
  proc_rows_max,                      & ! Maximum row in the PE
  model_variable_model_levels,        & ! Input model temperature
  model_nudged_var_mod_levs,          & ! Output model temperature
  data_variable_model_levels,         & ! Input data temperature
  relax_par,                          & ! size of timestep
  var_diag1,                          & ! Diagnostic array
  var_diag2,                          & ! Diagnostics
  var_diag3,                          & ! Diagnostics
  var_diag4,                          & ! Diagnostics
  var_diag5,                          & ! Diagnostics
  debug)                                ! Mis variables

USE nudging_control                     ! Access Nudging Diagnsocts
USE nudging_d1_defs                     ! Access D1 array information

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN)       :: model_levels    ! Number of levels in nectdf file
INTEGER, INTENT(IN)  :: proc_row_length_min  ! min column in grid
INTEGER, INTENT(IN)  :: proc_row_length_max  ! max column in grid
INTEGER, INTENT(IN)  :: proc_rows_min        ! min row in grid
INTEGER, INTENT(IN)  :: proc_rows_max        ! max row in grid

! Model variable before nudging
REAL, INTENT(IN)   :: model_variable_model_levels(              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

! Model variable after nudging
REAL, INTENT(IN)   :: model_nudged_var_mod_levs(                &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

! Data variable on model levels
REAL, INTENT(IN)      :: data_variable_model_levels (           &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max, 1:model_levels )

! Relaxation parameter (3D)
REAL, INTENT(IN)      :: relax_par (                            &
  proc_row_length_min:proc_row_length_max,                      & 
  proc_rows_min:proc_rows_max,1:model_levels)

! Output diagnostic 1 (data variable on model levels)
REAL, INTENT(INOUT)  :: var_diag1(                              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels)

! Output diagnostic 2 (model variable on model levels)
REAL, INTENT(INOUT)  :: var_diag2(                              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

! Output diagnostic 3 (tendency due to nudging)
REAL, INTENT(INOUT)  :: var_diag3(                              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

! Output diagnostic 4 (tendency due to model)
REAL, INTENT(INOUT)  :: var_diag4(                              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

! Output diagnostic 5 (3D relax par)
REAL, INTENT(INOUT)  :: var_diag5(                              &
  proc_row_length_min:proc_row_length_max,                      &
  proc_rows_min:proc_rows_max,1:model_levels )

INTEGER, INTENT(IN)  :: debug                ! debug flag
INTEGER              :: i, j, k              ! Loop variables
REAL                 :: lastvar              ! Holder in tendency calc
INTEGER, PARAMETER   :: n_nudge_diags = 10   ! Max number of diagnostics

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_CALC_DIAGS',zhook_in,zhook_handle)

!***********************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 5) THEN
  WRITE(OUT,*)                                                    &
  ' NUDGING_CALC_DIAGS: Entered Routine'
END IF

! Loop over the chosen regions to calc diagnostics

DO i=(proc_row_length_min), (proc_row_length_max)
  DO j=(proc_rows_min), (proc_rows_max)
    DO k=1, model_levels

! This is the value of the variable on the last timestep
      lastvar = var_diag2(i,j,k)

! 1st diagnostic = data on model levels
      IF(n_nudge_diags >= 1) THEN
        var_diag1(i,j,k) =                                         &
         data_variable_model_levels(i,j,k)
      END IF

! 2nd diagnostic = nudged variable
      IF(n_nudge_diags >= 2) THEN
        var_diag2(i,j,k) =                                         &
         model_nudged_var_mod_levs(i,j,k)
      END IF

! 3rd diagnostic = tendency due to nudging
      IF(n_nudge_diags >= 3) THEN
        var_diag3(i,j,k) =                                          &
         model_nudged_var_mod_levs(i,j,k)                           &
        - model_variable_model_levels(i,j,k)
      END IF

! 4th diagnostic = tendency due to all other factors
      IF(n_nudge_diags >= 4) THEN
        var_diag4(i,j,k) =                                           &
         model_variable_model_levels(i,j,k)                          &
        - lastvar
      END IF

! 5th diagnostic = relaxation parameter
      IF(n_nudge_diags >= 5) THEN
        var_diag5(i,j,k) =                                           &
         relax_par(i,j,k)
      END IF

    END DO      ! Loop over levels
  END DO        ! Loop over proc_rows
END DO          ! Loop over proc_row_length

! Standard Subroutine Exit Comment
IF(debug > 5) THEN
  WRITE(OUT,*)                                                       &
  ' NUDGING_CALC_DIAGS: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_CALC_DIAGS',zhook_out,zhook_handle)

END SUBROUTINE nudging_calc_diags

