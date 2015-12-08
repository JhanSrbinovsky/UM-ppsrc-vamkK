! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Converts data from analysis to model grid.

!  Part of the Nudged model (see nudging_main.F90)

!  Called from

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_ana2mod_grid( &
  global_row_length,             &  ! Length of analysis rows
  global_rows,                   &  ! Number of analysis rows
  file_levels,                   &  ! Number of analysis levels
  proc_row_length_min,           &  ! Minimum  column in local grid
  proc_row_length_max,           &  ! Maximum column in local grid
  proc_rows_min,                 &  ! Minimum row in local grid
  proc_rows_max,                 &  ! Maximum row in local grid
  variable_file_levels,          &  ! Variable on analysis grid
  variable_interp,               &  ! Variable on model grid
  base_lambda,                   &  ! base longitude for model
  delta_lambda,                  &  ! longitude increment for model
  base_phi,                      &  ! base latitude for model
  delta_phi,                     &  ! latitude increment for model
  base_lambda_ana,               &  ! base longitude for data
  delta_lambda_ana,              &  ! longitude increment for data
  base_phi_ana,                  &  ! base latitude for data
  delta_phi_ana,                 &  ! latitude increment for data
  debug)                            ! Debug flag


USE conversions_mod, ONLY: pi
USE nudging_control                 ! Use standard nudging switches

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook


IMPLICIT NONE


INTEGER, INTENT(IN) :: global_row_length    ! Length of rows in data grid
INTEGER, INTENT(IN) :: global_rows          ! Number of rows in data grid
INTEGER, INTENT(IN) :: file_levels          ! Levels in data
INTEGER, INTENT(IN) :: proc_row_length_min  ! Min. column in the PE
INTEGER, INTENT(IN) :: proc_row_length_max  ! Max. column in the PE
INTEGER, INTENT(IN) :: proc_rows_min        ! Min. row in the PE
INTEGER, INTENT(IN) :: proc_rows_max        ! Max. row in the PE

REAL, INTENT(IN)    :: base_lambda          ! base longitude for model
REAL, INTENT(IN)    :: delta_lambda         ! longitude increment for model
REAL, INTENT(IN)    :: base_phi             ! base latitude for model
REAL, INTENT(IN)    :: delta_phi            ! latitude increment for model

REAL, INTENT(IN)    :: base_lambda_ana      ! base longitude for data
REAL, INTENT(IN)    :: delta_lambda_ana     ! longitude increment for the data
REAL, INTENT(IN)    :: base_phi_ana         ! base latitude for data
REAL, INTENT(IN)    :: delta_phi_ana        ! latitude increment for data

INTEGER, INTENT(IN) :: debug                ! debug flag

! Variable on data grid
REAL, INTENT(IN) :: variable_file_levels                           &
  (global_row_length, global_rows, file_levels)

! Variable on ECMWF hybrid p-levels and model grid
REAL, INTENT(OUT) :: variable_interp(                              &
  proc_row_length_min:proc_row_length_max,                         &
  proc_rows_min:proc_rows_max,1:file_levels)

!**************************************
! End of I/O

REAL    :: base_lambda_anarad           ! base longitude in the data (radians)
REAL    :: delta_lambda_anarad          ! longitude eincrement in the data (radians)
REAL    :: base_phi_anarad              ! base latitude in teh data (radians)
REAL    :: delta_phi_anarad             ! latitude increment in the data (radians)
REAL    :: lambda_anarad                ! longitude in the data (radians)
REAL    :: phi_anarad                   ! latitude in teh data (radians)

REAL    :: lambda                       ! longitude in the model
REAL    :: phi                          ! latitude in the model
REAL    :: step1                        ! step for interpolation
REAL    :: step2                        ! step for interpolation

REAL    :: phi_ana_max                  ! Maximum longitude in the analysis
REAL    :: lambda_ana_max               ! Maximum latitude in the analysis
REAL    :: base_lambda_loc              ! Local base longitude
REAL    :: base_phi_loc                 ! Local base latitude

INTEGER :: i, j, k, ii, jj, kk, dummy   ! loop variables


INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_ANA2MOD_GRID',zhook_in,zhook_handle)


!**********************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) WRITE(OUT, *)                                           &
 'NUDGING_ANA2MOD_GRID: Entering Routine'

! Initialise output variable to 0
variable_interp(:,:,:) = 0.

! Convert angles from degrees to radians
! NB do it this way to ensure same value of PI as UM
base_lambda_anarad  = base_lambda_ana*(pi/180.)
delta_lambda_anarad = delta_lambda_ana*(pi/180.)
base_phi_anarad     = base_phi_ana*(pi/180.)
delta_phi_anarad    = delta_phi_ana*(pi/180.)

! Note starting coordinates of analysis for future use
lambda_anarad = base_lambda_anarad
phi_anarad    = base_phi_anarad

! Calculate local starting points
base_lambda_loc = base_lambda + (proc_row_length_min -1 )*delta_lambda
base_phi_loc = base_phi + (proc_rows_min -1 )*delta_phi

! Drop a line to let us know where we are
IF(debug > 15) WRITE(OUT,'(A,4(A,2e10.4),A)')                          &
  'NUDGING_ANA2MOD_GRID: Various Co-ordinates: ',                      &
  ', A:', lambda_anarad,  base_lambda_loc,                             &
  ', B:', delta_lambda_anarad,  delta_lambda,                          &
  ', C:', phi_anarad, base_phi_loc,                                    &
  ', D:', delta_phi_anarad, delta_phi,                                 &
  '.'

!**********************************************************************
! Locate strating point in teh analyses of phi
ii = 1

! Calculate maximum value of phi
phi_ana_max =  base_phi_anarad +                                       &
         (global_rows-1)*delta_phi_anarad

! Loop over phi finding strating point
! NB ensure we don't pass the maximum value
DO WHILE((phi_anarad+delta_phi_anarad).le.base_phi_loc.AND.            &
         (phi_anarad+delta_phi_anarad).le.phi_ana_max)

! whilst still not reached phi or end keep incrementing
  phi_anarad   = phi_anarad + delta_phi_anarad
  ii = ii+1
END DO

! set model phi to base phi
phi = base_phi_loc

! Loop over all rows (phi)
DO i=proc_rows_min, proc_rows_max

! Loop over phi finidng new starting point
! NB ensure we don't pass the maximum value
  DO WHILE((phi_anarad+delta_phi_anarad) <= phi .AND.                  &
         (phi_anarad+delta_phi_anarad) <= phi_ana_max)

! whilst still not reached phi or end keep incrementing
    phi_anarad   = phi_anarad + delta_phi_anarad
    ii = ii+1
  END DO

! Evaluate stride used for interpolation
  step1 = (phi-phi_anarad)/delta_phi_anarad

  IF (debug > 15) WRITE(OUT, '(A,2i4.4,3e10.4)')                       &
    ': NUDGING_ANA2MOD_GRID: PHI ', i, ii,  phi_anarad,                &
    base_phi_loc, phi, step1

!**********************************************
! Loop over all columns (lambda)

! Locate strating point in teh analyses of lambda
  jj = 1

! Calculate maximum value of phi
  lambda_ana_max =  base_lambda_anarad +                               &
         (global_row_length-1)*delta_lambda_anarad

! reset analysis lambda from previous loop
  lambda_anarad  = base_lambda_anarad

! Loop over phi finidng starting point
! NB ensure we don't pass the maximum value
  DO WHILE((lambda_anarad+delta_lambda_anarad) <= base_lambda_loc .AND. &
           (lambda_anarad+delta_lambda_anarad) <= lambda_ana_max)

! whilst still not reached phi or end keep incrementing
    lambda_anarad   = lambda_anarad + delta_lambda_anarad
    jj = jj+1
  END DO

! set model phi to base phi
  lambda = base_lambda_loc

! Loop over all rows (phi)
  DO j=proc_row_length_min, proc_row_length_max

! Loop over phi finidng strating point
! NB ensure we don't pass the maximum value
    DO WHILE((lambda_anarad+delta_lambda_anarad) <= lambda .AND.       &
           (lambda_anarad+delta_lambda_anarad) <= lambda_ana_max)

! whilst still not reached phi or end keep incrementing
      lambda_anarad   = lambda_anarad + delta_lambda_anarad
      jj = jj+1
    END DO

! Evaluate stride used for interpolation
    step2 = (lambda-lambda_anarad)
    step2 = step2/delta_lambda_anarad

    IF (debug > 15) WRITE(OUT, '(A,4i4.4,5e10.4)')                     &
      ' NUDGING_ANA2MOD_GRID: Lambda ', i, ii, j, jj, lambda_anarad,   &
      base_lambda, lambda, step2,  variable_file_levels(jj,   ii,  1)

! Loop over model levels
    DO k=1, file_levels

! DEPENDS ON: NUDGING_BI_INTERPOLATION_ZERO_D
      CALL nudging_bi_interpolation_zero_d(  &
        variable_file_levels(jj,   ii,  k),  &  ! Variable at x_0, y_0
        variable_file_levels(jj,   ii+1,k),  &  ! Variable at x_0, y_1
        variable_file_levels(jj+1, ii,  k),  &  ! Variable at x_1, y_0
        variable_file_levels(jj+1, ii+1,k),  &  ! Variable at x_1, y_1
        step2,   &  ! fraction between x_0 and x_1
        step1,                               &  ! fraction between y_0 and y_1
        variable_interp(j,i,k),   &  ! intrepolated variable
        0,    &  ! (Bi) Linear interpolation
        no_debug)               ! Debug flag

    END DO

! increment model lambda
    lambda = lambda + delta_lambda

  END DO      ! Loop over proc_row_length

! Increment model phi
  phi = phi + delta_phi

END DO        ! Loop over proc_rows

! Standard Subroutine Exit Comment
IF(debug > 10) WRITE(OUT, *)                                           &
 ' NUDGING_ANA2MOD_GRID: Leaving Routine'

IF (lhook) CALL dr_hook('NUDGING_ANA2MOD_GRID',zhook_out,zhook_handle)

END SUBROUTINE nudging_ana2mod_grid

