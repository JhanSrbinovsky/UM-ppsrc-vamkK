! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
MODULE vert_interp_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE vert_interp                                            &
                      (data_in, data_points,                      &
                       data_levels, desired_r,                    &
                       r_at_data, interp_order,                   &
                       data_out )

! Purpose:
!          Performs vertical interpolation of a field to a
!          desired r surface given the r value at each of
!          the data points. Where the desired surface is below/above
!          the top/bottom point, then a value is found by using
!          linear extrapolation from the nearest 2 points, except when
!          interp_order=2 is chosen, when values from the top/bottom
!          level are taken.
!          r can be any quantity that is monotonic increasing with
!          model level, such as height, or potential temperature
!          (assuming that the model is statically stable).

! Method:
!          Lagrange interpolation.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90

USE interpor_mod, ONLY: &
  interp_order_linear,                  &
  interp_order_linear_noex,             &
  interp_order_cubic

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) :: data_points    ! number of rows of data
INTEGER, INTENT(IN) :: data_levels    ! number of levels of data
INTEGER, INTENT(IN) :: interp_order   ! 1 = linear, 3 = cubic, 5=quintic
                                      ! 2 = linear with no extrapolation at top
                                      !     or bottom of model.

REAL, INTENT(IN) :: data_in (:, :)    ! input data field
REAL, INTENT(IN) :: r_at_data (:, :)  ! r surface
REAL, INTENT(IN) :: desired_r(:)      ! desired value to
                                      ! which data should
                                      ! be interpolated to

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) :: data_out (:)     ! output data

! Local variables
INTEGER :: j,k       ! loopers

INTEGER :: last      ! value at last point

INTEGER :: level_below(data_points)

REAL :: r_here
REAL :: r_here_plus
REAL :: r_here_plus2
REAL :: r_here_minus
REAL :: r_here_plus3
REAL :: r_here_minus2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('VERT_INTERP',zhook_in,zhook_handle)


! For scalar platforms we can take advantage of the fact that fields
! tend to vary smoothly by storing the last level we found and searching
! up or down from this point until the right level is found. This is not as
! simple to follow, but much quicker.
last=1
DO j=1, data_points

  level_below(j) = data_levels

  IF (r_at_data(j,last) <= desired_r(j)) THEN
    DO k=last+1, data_levels
      IF (r_at_data(j,k) > desired_r(j)) THEN
        level_below(j) = k-1
        EXIT
      END IF
    END DO
  ELSE
    DO k = last-1, 1, -1
      IF (r_at_data(j,k) <= desired_r(j)) THEN
        level_below(j) = k
        EXIT
      END IF
    END DO
  END IF

  last = MAX(level_below(j),1)
END DO

IF (interp_order  /=  interp_order_linear_noex) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                  &
!$OMP&         SHARED(data_points, desired_r, r_at_data,          &
!$OMP&                data_out, data_in, data_levels)             &
!$OMP&         PRIVATE(j)
DO j = 1, data_points

! If requested level is above top of model, do linear
! extrapolation using data on top and second top levels.
    IF ( desired_r(j)  >   r_at_data(j,data_levels) ) THEN
      data_out(j) = data_in(j,data_levels) + (desired_r(j)        &
      - r_at_data(j,data_levels)) * (data_in(j,data_levels)       &
      - data_in(j,data_levels-1))/(r_at_data(j,data_levels)       &
      - r_at_data(j,data_levels-1))
    ELSE IF (desired_r(j) == r_at_data(j,data_levels) ) THEN
      data_out(j) = data_in(j,data_levels)
    END IF

! If requested level is below bottom of model, do linear
! extrapolation using data on first and second levels.
    IF ( desired_r(j)  <   r_at_data(j,1) ) THEN
      data_out(j) = data_in(j,1) + (desired_r(j)                  &
      - r_at_data(j,1)) * (data_in(j,1)                           &
      - data_in(j,2))/(r_at_data(j,1) - r_at_data(j,2))
    ELSE IF (desired_r(j) == r_at_data(j,1) ) THEN
      data_out(j) = data_in(j,1)
    END IF

END DO
!$OMP END PARALLEL DO

ELSE ! No linear extrapolation at top or bottom

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                  &
!$OMP&         SHARED(data_points, desired_r, r_at_data, data_in, &
!$OMP&                data_out, data_levels)                      &
!$OMP&         PRIVATE(j)
  DO j = 1, data_points

!         Top : Set to top input data

    IF ( desired_r(j)  >=  r_at_data(j,data_levels) ) THEN
      data_out(j) = data_in(j,data_levels)
    END IF

!       Bottom : Set to bottom input data

    IF ( desired_r(j)  <=  r_at_data(j,1) ) THEN
      data_out(j) = data_in(j,1)
    END IF

  END DO
!$OMP END PARALLEL DO

END IF

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP&         SHARED(level_below, data_levels, data_out, desired_r,    &
!$OMP&                r_at_data, data_in, data_points, interp_order)    &
!$OMP&         PRIVATE(j, r_here, r_here_minus, r_here_plus,            &
!$OMP&                 r_here_plus2, r_here_plus3, r_here_minus2)
DO j = 1, data_points

    IF (level_below(j)  ==  0.OR.                                 &
        level_below(j)  ==  data_levels) THEN

       ! Output data has already been calculated - do nothing...
        data_out (j) = data_out (j)

    ELSE IF ( ABS( desired_r(j) - r_at_data(j, level_below(j)))   &
                         < EPSILON( desired_r(j))) THEN

        ! Levels are essentially the same, copy in to out.
        data_out(j) = data_in(j,level_below(j))

    ELSE IF (level_below(j)  ==  1 .OR.                           &
        level_below(j)  ==  data_levels - 1                       &
        .OR. interp_order  ==  interp_order_linear                &
        .OR. interp_order  ==  interp_order_linear_noex) THEN
! linearly interpolate
      data_out (j) = ( (desired_r(j) -                            &
                          r_at_data(j,level_below(j)) )           &
                          * data_in (j,level_below(j)+1)          &
                         -(desired_r(j) -                         &
                           r_at_data(j,level_below(j)+1)) *       &
                           data_in (j,level_below(j)) ) /         &
                        ( r_at_data(j,level_below(j)+1) -         &
                          r_at_data(j,level_below(j)) )

    ELSE IF (level_below(j)  ==  2 .OR.                           &
             level_below(j)  ==  data_levels - 2                  &
             .OR. interp_order  ==  interp_order_cubic ) THEN

! cubicly interpolate

      r_here_minus = r_at_data(j,level_below(j)-1)
      r_here = r_at_data(j,level_below(j))
      r_here_plus = r_at_data(j,level_below(j)+1)
      r_here_plus2 = r_at_data(j,level_below(j)+2)

      data_out (j) = ( (desired_r(j) - r_here) *                  &
                           (desired_r(j) - r_here_plus )*         &
                           (desired_r(j) - r_here_plus2 ) ) /     &
                         ( (r_here_minus - r_here) *              &
                           (r_here_minus - r_here_plus )*         &
                           (r_here_minus - r_here_plus2 ) ) *     &
                         data_in (j,level_below(j)-1) +           &
                         ( (desired_r(j) - r_here_minus) *        &
                           (desired_r(j) - r_here_plus )*         &
                           (desired_r(j) - r_here_plus2 ) ) /     &
                         ( (r_here - r_here_minus) *              &
                           (r_here - r_here_plus )*               &
                           (r_here - r_here_plus2 ) ) *           &
                         data_in (j,level_below(j)) +             &
                         ( (desired_r(j) - r_here_minus) *        &
                           (desired_r(j) - r_here )*              &
                           (desired_r(j) - r_here_plus2 ) ) /     &
                         ( (r_here_plus - r_here_minus) *         &
                           (r_here_plus - r_here )*               &
                           (r_here_plus - r_here_plus2 ) ) *      &
                         data_in (j,level_below(j)+1) +           &
                         ( (desired_r(j) - r_here_minus) *        &
                           (desired_r(j) - r_here )*              &
                           (desired_r(j) - r_here_plus ) ) /      &
                         ( (r_here_plus2 - r_here_minus) *        &
                           (r_here_plus2 - r_here )*              &
                           (r_here_plus2 - r_here_plus ) ) *      &
                         data_in (j,level_below(j)+2)

    ELSE
! interpolate quinticly

      r_here_minus2 = r_at_data(j,level_below(j)-2)
      r_here_minus = r_at_data(j,level_below(j)-1)
      r_here = r_at_data(j,level_below(j))
      r_here_plus = r_at_data(j,level_below(j)+1)
      r_here_plus2 = r_at_data(j,level_below(j)+2)
      r_here_plus3 = r_at_data(j,level_below(j)+3)

      data_out (j) = ((desired_r(j) - r_here_minus) *             &
                         (desired_r(j) - r_here )*                &
                         (desired_r(j) - r_here_plus )*           &
                         (desired_r(j) - r_here_plus2 )*          &
                         (desired_r(j) - r_here_plus3 ))/         &
                       ( (r_here_minus2 - r_here_minus) *         &
                         (r_here_minus2 - r_here )*               &
                         (r_here_minus2 - r_here_plus )*          &
                         (r_here_minus2 - r_here_plus2 )*         &
                         (r_here_minus2 - r_here_plus3 ) ) *      &
                         data_in (j,level_below(j)-2) +           &
                       ((desired_r(j) - r_here_minus2) *          &
                         (desired_r(j) - r_here )*                &
                         (desired_r(j) - r_here_plus )*           &
                         (desired_r(j) - r_here_plus2 )*          &
                         (desired_r(j) - r_here_plus3 ))/         &
                       ( (r_here_minus - r_here_minus2) *         &
                         (r_here_minus - r_here )*                &
                         (r_here_minus - r_here_plus )*           &
                         (r_here_minus - r_here_plus2 )*          &
                         (r_here_minus - r_here_plus3 ) ) *       &
                         data_in (j,level_below(j)-1) +           &
                       ((desired_r(j) - r_here_minus2) *          &
                         (desired_r(j) - r_here_minus )*          &
                         (desired_r(j) - r_here_plus )*           &
                         (desired_r(j) - r_here_plus2 )*          &
                         (desired_r(j) - r_here_plus3 ))/         &
                       ( (r_here - r_here_minus2) *               &
                         (r_here - r_here_minus )*                &
                         (r_here - r_here_plus )*                 &
                         (r_here - r_here_plus2 )*                &
                         (r_here - r_here_plus3 ) ) *             &
                         data_in (j,level_below(j)) +             &
                       ((desired_r(j) - r_here_minus2) *          &
                         (desired_r(j) - r_here_minus )*          &
                         (desired_r(j) - r_here )*                &
                         (desired_r(j) - r_here_plus2 )*          &
                         (desired_r(j) - r_here_plus3 ))/         &
                       ( (r_here_plus - r_here_minus2) *          &
                         (r_here_plus - r_here_minus )*           &
                         (r_here_plus - r_here )*                 &
                         (r_here_plus - r_here_plus2 )*           &
                         (r_here_plus - r_here_plus3 ) ) *        &
                         data_in (j,level_below(j)+1) +           &
                       ((desired_r(j) - r_here_minus2) *          &
                         (desired_r(j) - r_here_minus )*          &
                         (desired_r(j) - r_here )*                &
                         (desired_r(j) - r_here_plus )*           &
                         (desired_r(j) - r_here_plus3 ))/         &
                       ( (r_here_plus2 - r_here_minus2) *         &
                         (r_here_plus2 - r_here_minus )*          &
                         (r_here_plus2 - r_here )*                &
                         (r_here_plus2 - r_here_plus )*           &
                         (r_here_plus2 - r_here_plus3 ) ) *       &
                         data_in (j,level_below(j)+2) +           &
                       ((desired_r(j) - r_here_minus2) *          &
                         (desired_r(j) - r_here_minus )*          &
                         (desired_r(j) - r_here )*                &
                         (desired_r(j) - r_here_plus )*           &
                         (desired_r(j) - r_here_plus2 ))/         &
                       ( (r_here_plus3 - r_here_minus2) *         &
                         (r_here_plus3 - r_here_minus )*          &
                         (r_here_plus3 - r_here )*                &
                         (r_here_plus3 - r_here_plus )*           &
                         (r_here_plus3 - r_here_plus2 ) ) *       &
                         data_in (j,level_below(j)+3)

    END IF

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('VERT_INTERP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vert_interp
END MODULE vert_interp_mod
