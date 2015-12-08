! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Converts variable from fixed pressure levels to model levels

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_VARLOADER

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_pres_to_mod(       &
  global_row_length,                  & ! Length of model row
  global_rows,                        & ! Number of model rows
  model_levels,                       & ! Number of model levels
  file_levels,                        & ! Number of file levels
  proc_row_length_min,                & ! Minimum column in the PE
  proc_row_length_max,                & ! Maximum column in the PE
  proc_rows_min,                      & ! Minimum row in the PE
  proc_rows_max,                      & ! Maximum row in the PE
  model_pressure_model_levels,        & ! Pressure on model levels
  data_pressure,                      & ! Pressure on data levels
  data_variable_pressure_levels,      & ! Variable on ecmwf pressure levels
  data_variable_model_levels,         & ! Variable interpolated on to model levels
  debug)                                ! Miscellanea

USE nudging_control                     ! Use standard nudging switches

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: global_row_length    ! No. of pts in a row
INTEGER, INTENT(IN) :: global_rows          ! No of rows
INTEGER, INTENT(IN) :: model_levels         ! No. of model levels
INTEGER, INTENT(IN) :: file_levels          ! No. of model levels
INTEGER, INTENT(IN) :: proc_row_length_min  ! Min. column in the PE
INTEGER, INTENT(IN) :: proc_row_length_max  ! Max. column in the PE
INTEGER, INTENT(IN) :: proc_rows_min        ! Min. row in the PE
INTEGER, INTENT(IN) :: proc_rows_max        ! Max. row in the PE

! model pressure
REAL, INTENT(IN) :: model_pressure_model_levels(                       &
  proc_row_length_min:proc_row_length_max,                             &
  proc_rows_min:proc_rows_max,1:model_levels)

! data pressure on data levels
REAL, INTENT(IN):: data_pressure                                       &
                  (file_levels)

! data variable on data levels
REAL, INTENT(IN):: data_variable_pressure_levels                       &
 (global_row_length, global_rows, file_levels)

! data variable on model levels
REAL, INTENT(OUT) :: data_variable_model_levels(                       &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

INTEGER, INTENT(IN) :: debug          ! Debug flag
INTEGER             :: i, j, k, l, m  !  Loop counters
REAL                :: step           ! Frac for interpolation

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_PRES_TO_MOD',zhook_in,zhook_handle)

!*******************************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
 ': NUDGING_PRES_TO_MOD: Entering Routine'
END IF

! Loop round all latitudes and all longitudes
DO i= proc_rows_min, proc_rows_max
  DO j= proc_row_length_min, proc_row_length_max

! Initialise the data loop
    l=1
    m=l+1
! Increment L until pressure levels straddle 1st model level pressure
    DO WHILE(model_pressure_model_levels(j,i,1) <                     &
       data_pressure(l+1) .AND.                            &
       m < file_levels .AND. data_pressure(l+1) > 0)
      l= l+1
      m= l+1
    END DO

! Loop over model levels interpolating on to these levels
    DO k=1, model_levels

! Calculate step for intrepolation
      step = LOG(model_pressure_model_levels(j,i,k)) -                 &
             LOG(data_pressure(l))
      step = step/(LOG(data_pressure(l+1))                             &
             - LOG(data_pressure(l)))

! If in debug mode provide some information on the interpolation
      IF(debug > 15 .AND.                                              &
        i == proc_row_length_min.AND.j == proc_rows_min) THEN
        WRITE(OUT,*)                                                   &
         'NUDGING_PRES_TO_MOD: Press calc uses ',                      &
          model_pressure_model_levels(j,i,k),                          &
          data_pressure(l),                                            &
          data_pressure(l+1),                                          &
          step, k, l
      END IF

!******************************************************************
! Limit large extrapolations if this option is chosen

      IF(l_extrapcutoff) THEN
        IF(step > 1.) THEN

! If in debuf mode let us know that we are extrapolating
          IF(debug > 15) THEN
            WRITE(OUT,*)                                               &
              'NUDGING_PRES_TO_MOD: Extrapolation',                    &
               step, i, j, k
          END IF

! Limit extrapolation to 2
          IF(step > 2.) THEN
            step=2.
          END IF

        ELSE IF(step < 0) THEN

! If in debuf mode let us know that we are extrapolating
          IF(debug > 15) THEN
            WRITE(OUT,*)                                                 &
             'NUDGING_PRES_TO_MOD: Extrapolation',                       &
              step, i, j, k
          END IF

! Limit extraoplation to -1
          IF(step < -1.) THEN
            step=-1.
          END IF

        END IF     ! step
      END IF       ! l_extrapcutoff

!*****************************************************************
! Do the interpolation
! DEPENDS ON: nudging_interpolation_zero_d
      CALL nudging_interpolation_zero_d(          &
        data_variable_pressure_levels(j,i,(l)),   &   ! Starting variable
        data_variable_pressure_levels(j,i,(l+1)), &   ! Finishing variable
        step,                                     &   ! Fraction between 1 & 2
        data_variable_model_levels(j,i,k),        &   ! Inerpolated variable
        0,                                        &   ! Linear interp
        no_debug)                                      ! Debug flag

! If in debug mode provide information about what we're doing
      IF(debug > 15 .AND.                                              &
        i == proc_row_length_min.AND.j == proc_rows_min) THEN
        WRITE(OUT,*)                                                   &
          'NUDGING_PRES_TO_MOD: Press calc uses ',                     &
           data_variable_model_levels(j,i,k),                          &
           data_variable_pressure_levels (j,i,l),                      &
           data_variable_pressure_levels(j,i, (l+1)),                  &
           step, k, l
      END IF

! If the pressure of next ECMWF gt than  next model level
! increment the ECMWF level
! Check not exceeded the number of ECMWF levels
      m=l+1
      DO WHILE(model_pressure_model_levels(j,i,(k+1)) <               &
                data_pressure(m).AND.                                 &
                m < file_levels.AND.data_pressure(l+1) > 0)
        l=l+1
        m=l+1
      END DO

    END DO     ! levels
  END DO       ! rows
END DO         ! row_length

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
 ': NUDGING_PRES_TO_MOD: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_PRES_TO_MOD',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_pres_to_mod
