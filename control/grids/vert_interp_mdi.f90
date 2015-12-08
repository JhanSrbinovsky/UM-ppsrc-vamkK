! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine vert_interp_mdi
MODULE vert_interp_mdi_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE vert_interp_mdi(                                       &
                           data_in, row_length, data_rows,        &
                           data_levels, desired_r,                &
                           halo_x1, halo_y1,                      &
                           halo_x2, halo_y2,                      &
                           r_at_data, interp_order,               &
                           mdi, data_out )

! Purpose:
!          Performs vertical interpolation of a field to a
!          desired r surface given the r value at each of
!          the data points. Where the desired surface is below/above
!          the bottom/top data point a missing data indicator is
!          returned.
!          r can be any quantity that is monotonic increasing with
!          model level, such as height, or potential temperature
!          (assuming that the model is statically stable).
!          Interpolation order can be specified via the argument
!          interp_order and may be linear, cubic or quintic.

! Method:
!          Lagrange interpolation.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90

USE interpor_mod, ONLY: &
  interp_order_linear,                  &
  interp_order_cubic

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) :: row_length    ! number of points on a row
INTEGER, INTENT(IN) :: data_rows     ! number of rows of data
INTEGER, INTENT(IN) :: data_levels   ! number of levels of data
INTEGER, INTENT(IN) :: interp_order  ! 1 = linear, 3 = cubic, 5=quintic
INTEGER, INTENT(IN) :: halo_x1
INTEGER, INTENT(IN) :: halo_y1
INTEGER, INTENT(IN) :: halo_x2
INTEGER, INTENT(IN) :: halo_y2

REAL, INTENT(IN) :: desired_r        ! desired value to which data should be
                                     ! interpolated to.
REAL, INTENT(IN) :: mdi              ! missing data indicator.

REAL, INTENT(IN) :: data_in (1-halo_x1:row_length+halo_x1,              &
                             1-halo_y1:data_rows+halo_y1, data_levels)
REAL, INTENT(IN) :: r_at_data (1-halo_x2:row_length+halo_x2,            &
                               1-halo_y2:data_rows+halo_y2, data_levels)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) :: data_out (row_length, data_rows)

! Local variables

INTEGER :: i,j,k ! loopers
INTEGER :: last
INTEGER :: level_below(row_length, data_rows)

REAL :: r_here
REAL :: r_here_plus
REAL :: r_here_plus2
REAL :: r_here_minus
REAL :: r_here_plus3
REAL :: r_here_minus2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('VERT_INTERP_MDI',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
! ----------------------------------------------------------------------

! For scalar platforms we can take advantage of the fact that fields
! tend to vary smoothly by storing the last level we found and searching
! up or down from this point until the right level is found. This is not as
! simple to follow, but much quicker.
last = 1
DO j = 1, data_rows
  DO i = 1, row_length

    level_below(i,j) = 0

    IF ( r_at_data(i,j,last) <= desired_r ) THEN
      DO k = last+1, data_levels - 1
        IF ( r_at_data(i,j,k) > desired_r ) THEN
          level_below(i,j) = k-1
          EXIT
        END IF
      END DO
    ELSE
      DO k = last-1, 1, -1
        IF ( r_at_data(i,j,k) <= desired_r ) THEN
          level_below(i,j) = k
          EXIT
        END IF
      END DO
    END IF

    last = MAX( level_below(i,j), 1)

  END DO
END DO

! if requested level is above top of model, set to zero, which will
! be converted to missing data indicator.
DO j = 1, data_rows
  DO i = 1, row_length
    IF ( desired_r  >   r_at_data(i,j,data_levels) ) THEN
      level_below(i,j) = 0
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

DO j = 1, data_rows
  DO i = 1, row_length

    IF (level_below(i,j)  ==  0) THEN
      data_out(i,j) = mdi

    ELSE IF (level_below(i,j)  ==  1 .OR.                         &
             level_below(i,j)  ==  data_levels - 1                &
             .OR. interp_order  ==  interp_order_linear ) THEN
! linearly interpolate
      data_out (i,j) = ( (desired_r -                             &
                          r_at_data(i,j,level_below(i,j)) )       &
                          * data_in (i,j,level_below(i,j)+1)      &
                         -(desired_r -                            &
                           r_at_data(i,j,level_below(i,j)+1)) *   &
                           data_in (i,j,level_below(i,j)) ) /     &
                        ( r_at_data(i,j,level_below(i,j)+1) -     &
                          r_at_data(i,j,level_below(i,j)) )

    ELSE IF (level_below(i,j)  ==  2 .OR.                         &
             level_below(i,j)  ==  data_levels - 2                &
             .OR. interp_order  ==  interp_order_cubic ) THEN

! cubicly interpolate

      r_here_minus = r_at_data(i,j,level_below(i,j)-1)
      r_here = r_at_data(i,j,level_below(i,j))
      r_here_plus = r_at_data(i,j,level_below(i,j)+1)
      r_here_plus2 = r_at_data(i,j,level_below(i,j)+2)

      data_out (i,j) = ( (desired_r - r_here) *                   &
                                (desired_r - r_here_plus )*       &
                                (desired_r - r_here_plus2 ) ) /   &
                              ( (r_here_minus - r_here) *         &
                                (r_here_minus - r_here_plus )*    &
                                (r_here_minus - r_here_plus2 ) ) *&
                              data_in (i,j,level_below(i,j)-1) +  &
                              ( (desired_r - r_here_minus) *      &
                                (desired_r - r_here_plus )*       &
                                (desired_r - r_here_plus2 ) ) /   &
                              ( (r_here - r_here_minus) *         &
                                (r_here - r_here_plus )*          &
                                (r_here - r_here_plus2 ) ) *      &
                              data_in (i,j,level_below(i,j)) +    &
                              ( (desired_r - r_here_minus) *      &
                                (desired_r - r_here )*            &
                                (desired_r - r_here_plus2 ) ) /   &
                              ( (r_here_plus - r_here_minus) *    &
                                (r_here_plus - r_here )*          &
                                (r_here_plus - r_here_plus2 ) ) * &
                              data_in (i,j,level_below(i,j)+1) +  &
                              ( (desired_r - r_here_minus) *      &
                                (desired_r - r_here )*            &
                                (desired_r - r_here_plus ) ) /    &
                              ( (r_here_plus2 - r_here_minus) *   &
                                (r_here_plus2 - r_here )*         &
                                (r_here_plus2 - r_here_plus ) ) * &
                              data_in (i,j,level_below(i,j)+2)

    ELSE
! interpolate quinticly

      r_here_minus2 = r_at_data(i,j,level_below(i,j)-2)
      r_here_minus = r_at_data(i,j,level_below(i,j)-1)
      r_here = r_at_data(i,j,level_below(i,j))
      r_here_plus = r_at_data(i,j,level_below(i,j)+1)
      r_here_plus2 = r_at_data(i,j,level_below(i,j)+2)
      r_here_plus3 = r_at_data(i,j,level_below(i,j)+3)

      data_out (i,j) = ((desired_r - r_here_minus) *              &
                              (desired_r - r_here )*              &
                              (desired_r - r_here_plus )*         &
                              (desired_r - r_here_plus2 )*        &
                              (desired_r - r_here_plus3 ))/       &
                            ( (r_here_minus2 - r_here_minus) *    &
                              (r_here_minus2 - r_here )*          &
                              (r_here_minus2 - r_here_plus )*     &
                              (r_here_minus2 - r_here_plus2 )*    &
                              (r_here_minus2 - r_here_plus3 ) ) * &
                              data_in (i,j,level_below(i,j)-2) +  &
                            ((desired_r - r_here_minus2) *        &
                              (desired_r - r_here )*              &
                              (desired_r - r_here_plus )*         &
                              (desired_r - r_here_plus2 )*        &
                              (desired_r - r_here_plus3 ))/       &
                            ( (r_here_minus - r_here_minus2) *    &
                              (r_here_minus - r_here )*           &
                              (r_here_minus - r_here_plus )*      &
                              (r_here_minus - r_here_plus2 )*     &
                              (r_here_minus - r_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)-1) +  &
                            ((desired_r - r_here_minus2) *        &
                              (desired_r - r_here_minus )*        &
                              (desired_r - r_here_plus )*         &
                              (desired_r - r_here_plus2 )*        &
                              (desired_r - r_here_plus3 ))/       &
                            ( (r_here - r_here_minus2) *          &
                              (r_here - r_here_minus )*           &
                              (r_here - r_here_plus )*            &
                              (r_here - r_here_plus2 )*           &
                              (r_here - r_here_plus3 ) ) *        &
                              data_in (i,j,level_below(i,j)) +    &
                            ((desired_r - r_here_minus2) *        &
                              (desired_r - r_here_minus )*        &
                              (desired_r - r_here )*              &
                              (desired_r - r_here_plus2 )*        &
                              (desired_r - r_here_plus3 ))/       &
                            ( (r_here_plus - r_here_minus2) *     &
                              (r_here_plus - r_here_minus )*      &
                              (r_here_plus - r_here )*            &
                              (r_here_plus - r_here_plus2 )*      &
                              (r_here_plus - r_here_plus3 ) ) *   &
                              data_in (i,j,level_below(i,j)+1) +  &
                            ((desired_r - r_here_minus2) *        &
                              (desired_r - r_here_minus )*        &
                              (desired_r - r_here )*              &
                              (desired_r - r_here_plus )*         &
                              (desired_r - r_here_plus3 ))/       &
                            ( (r_here_plus2 - r_here_minus2) *    &
                              (r_here_plus2 - r_here_minus )*     &
                              (r_here_plus2 - r_here )*           &
                              (r_here_plus2 - r_here_plus )*      &
                              (r_here_plus2 - r_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)+2) +  &
                            ((desired_r - r_here_minus2) *        &
                              (desired_r - r_here_minus )*        &
                              (desired_r - r_here )*              &
                              (desired_r - r_here_plus )*         &
                              (desired_r - r_here_plus2 ))/       &
                            ( (r_here_plus3 - r_here_minus2) *    &
                              (r_here_plus3 - r_here_minus )*     &
                              (r_here_plus3 - r_here )*           &
                              (r_here_plus3 - r_here_plus )*      &
                              (r_here_plus3 - r_here_plus2 ) ) *  &
                              data_in (i,j,level_below(i,j)+3)

    END IF

  END DO
END DO

! end of routine

IF (lhook) CALL dr_hook('VERT_INTERP_MDI',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vert_interp_mdi
END MODULE vert_interp_mdi_mod
