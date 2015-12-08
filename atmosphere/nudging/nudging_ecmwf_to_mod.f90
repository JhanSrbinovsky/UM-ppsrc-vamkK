! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  INterpolate variable from ECMWF pressure levels to
!  UM model levels linearly in lnP

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_VARLOADER.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------

SUBROUTINE nudging_ecmwf_to_mod(     &
  global_row_length,                 &  ! Length of model row
  global_rows,                       &  ! Number of model rows
  model_levels,                      &  ! Number of model levels
  file_levels,                       &  ! Number of file levels
  proc_row_length_min,               &  ! Minimum column in the PE
  proc_row_length_max,               &  ! Maximum column in the PE
  proc_rows_min,                     &  ! Minimum row in the PE
  proc_rows_max,                     &  ! Maximum row in the PE
  model_pressure_model_levels,       &  ! Pressure on model levels
  data_pressure_ecmwf_levels,        &  ! Surface pressure
  data_variable_ecmwf_levels,        &  ! Variable on ecmwf levels (In)
  data_variable_model_levels,        &  ! Variable on model levels (Out)
  debug)                                ! Miscellanea

USE nudging_control              ! Use standard nudging switches
USE nudging_ecmwf_60level_def    ! ECMWF pressure levels definition module

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: global_row_length      ! model row length
INTEGER, INTENT(IN) :: global_rows            ! model rows
INTEGER, INTENT(IN) :: model_levels           ! Model levels
INTEGER, INTENT(IN) :: file_levels            ! Ecmwf levels

INTEGER, INTENT(IN) :: proc_row_length_min   ! Min. column in PE
INTEGER, INTENT(IN) :: proc_row_length_max   ! Max. column in PE
INTEGER, INTENT(IN) :: proc_rows_min         ! Min. row in PE
INTEGER, INTENT(IN) :: proc_rows_max         ! Max. row in PE

! Pressure on model levels
REAL, INTENT(IN) :: model_pressure_model_levels(                      &
  proc_row_length_min:proc_row_length_max,                            &
  proc_rows_min:proc_rows_max, model_levels)

! Pressure on ECMWF levels
REAL, INTENT(IN) :: data_pressure_ecmwf_levels(                        &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,file_levels)

! Variable on ecmwf levels
REAL, INTENT(IN) :: data_variable_ecmwf_levels(                        &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,file_levels)

! Variable interpolated on to model levels
REAL, INTENT(OUT) :: data_variable_model_levels(                       &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,model_levels)

! Various debugging, error functions etc.
INTEGER, INTENT(IN) :: debug           ! Debug flag

INTEGER             :: i, j, k, l, m   ! Loop counters
REAL                :: step            ! Frac. for interp.

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_ECMWF_TO_MOD',zhook_in,zhook_handle)

!*****************************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
  ': NUDGING_ECMWF_TO_MOD: Entering Routine'
END IF

! Loop round all latitudes and all longitudes
DO i= proc_rows_min, proc_rows_max
  DO j= proc_row_length_min, proc_row_length_max

! Initialise the data loop
    l=1
    m=l+1

! Increment L until P levels straddle 1st model P level
    DO WHILE(model_pressure_model_levels(j,i,1) <                      &
            data_pressure_ecmwf_levels(j,i,m).AND.                     &
            m < file_levels)
      l= l+1
      m= l+1
    END DO

! Loop over model levels interpolating on to these levels
    DO k=1, model_levels

! Calculate interpolation amount
      step = LOG(model_pressure_model_levels(j,i,k)) -                 &
             LOG(data_pressure_ecmwf_levels (j,i,l))
      step = step/(LOG(data_pressure_ecmwf_levels(j,i,(l+1)))          &
            - LOG(data_pressure_ecmwf_levels(j,i,l)))

! If in debug mode provide some more information
      IF(debug > 15 .AND.                                              &
        i == proc_row_length_min .AND. j == proc_rows_min) THEN
        WRITE(OUT,*)                                                   &
          'NUDGING_ECMWF_TO_MOD: Press calc uses ',                    &
           model_pressure_model_levels(j,i,k),                         &
           data_pressure_ecmwf_levels (j,i,l),                         &
           data_pressure_ecmwf_levels(j,i,m),                          &
           step, k, l
      END IF

!****************************************************
! Limit large extrapolations if this option is chosen

      IF(l_extrapcutoff) THEN
        IF(step > 1) THEN

! If in debuf mode let us know that we are extrapolating
          IF(debug > 15) THEN
            WRITE(OUT,*)                                               &
              'NUDGING_ECMWF_TO_MOD: Extrapolation',                   &
              step, i, j, k
          END IF

! Cutoff extrapolations at 2
          IF(step > 2.) THEN
            step=2.
          END IF

        ELSE IF(step < 0) THEN

! If in debuf mode let us know that we are extrapolating
          IF(debug > 15) THEN
            WRITE(OUT,*)                                                &
            'NUDGING_ECMWF_TO_MOD: Extrapolation',                      &
             step, i, j, k
          END IF

! Cutoff extrapolations at -1
          IF(step < -1.) THEN
            step = -1.
          END IF

        END IF   ! Step 
      END IF     ! L_extrapcuttof

!***********************************************************
! Interpolate the analyses on to the model levels
! DEPENDS ON: nudging_interpolation_zero_d
      CALL nudging_interpolation_zero_d(    &
        data_variable_ecmwf_levels(j,i,l),  &     ! Variable at start
        data_variable_ecmwf_levels(j,i,m),  &     ! Variable at end
        step,                               &     ! Frac between start and end
        data_variable_model_levels(j,i,k),  &     ! Interpolated variable
        0,                                  &     ! Linear interpolation
        no_debug)                                 ! Debug flag

! If in debug mode print information about the interpolation
      IF(debug > 15 .AND.                                              &
        i == proc_row_length_min .AND. j == proc_rows_min) THEN
        WRITE(OUT,*)                                                   &
          'NUDGING_ECMWF_TO_MOD: Variable calc uses ',                 &
           data_variable_model_levels(j,i,k),                          &
           data_variable_ecmwf_levels(j,i,l),                          &
           data_variable_ecmwf_levels(j,i,m),                          &
           step, k, l
      END IF

! If the pressure of the next ECMWF is greater than the next
! model level then increment to the ECMWF level
! Also check that we have not exceeded the number of ECMWF or model
! levels
      IF(k < model_levels) THEN
        m=l+1
        DO WHILE(model_pressure_model_levels(j,i,(k+1)) <           &
                data_pressure_ecmwf_levels(j,i,m) .AND.             &
                 m < file_levels)
          l=l+1
          m=l+1
        END  DO
      END IF

    END DO      ! Loop over levels
  END DO        ! Loop over proc_row_length
END DO          ! Loop over proc_rows

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT, *)  'PE',                                                 &
  ' NUDGING_ECMWF_TO_MOD: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_ECMWF_TO_MOD',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_ecmwf_to_mod

