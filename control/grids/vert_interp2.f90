! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine vert_interp2
MODULE vert_interp2_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE vert_interp2(                                        &
                         data_in, row_length, data_rows,        &
                         data_levels, desired_p,                &
                         halo_x1, halo_y1,                      &
                         halo_x2, halo_y2,                      &
                         p_at_data, interp_order,               &
                         data_out )

! Purpose:
!          Performs vertical interpolation of a field to a
!          desired p surface given the value of p at each of
!          the data points. Where the desired surface is below/above
!          the bottom/top data point the lowest/highest level data
!          value is returned.
!          p can be any quantity that is monotonic decreasing with
!          model level, such as pressure, or Exner pressure.
!          It is usual, and recommended, that the p surface be exner.
!          This routine can  do linear, cubic or quintic interpolation,
!          as specified by the interp_order argument.

! Method:
!          Lagrange interpolation.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

USE interpor_mod, ONLY: &
  interp_order_linear,                  &
  interp_order_cubic

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT (IN) :: row_length     ! number of points on a row
INTEGER, INTENT (IN) :: data_rows      ! number of rows of data
INTEGER, INTENT (IN) :: data_levels    ! number of levels of data
INTEGER, INTENT (IN) :: interp_order   ! 1 = linear, 3 = cubic, 5=quintic
INTEGER, INTENT (IN) :: halo_x1
INTEGER, INTENT (IN) :: halo_x2
INTEGER, INTENT (IN) :: halo_y1
INTEGER, INTENT (IN) :: halo_y2

REAL, INTENT (IN) :: desired_p         ! desired value to which data
                                       ! should be  interpolated to.

REAL, INTENT (IN) :: data_in (1-halo_x1:row_length+halo_x1,       &
                              1-halo_y1:data_rows+halo_y1, data_levels)
REAL, INTENT (IN) :: p_at_data (1-halo_x2:row_length+halo_x2,     &
                                1-halo_y2:data_rows+halo_y2, data_levels)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT (OUT) :: data_out (row_length, data_rows)

! Local variables

INTEGER :: i, j, k   ! Loopers

INTEGER :: level_below(row_length, data_rows)

INTEGER :: last  ! level from last loop

REAL :: p_here
REAL :: p_here_plus
REAL :: p_here_plus2
REAL :: p_here_minus
REAL :: p_here_plus3
REAL :: p_here_minus2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('VERT_INTERP2',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!   (if requested level is above top of model, set to -1, which will
!    use value from top level)
! ----------------------------------------------------------------------

! For scalar platforms we can take advantage of the fact that fields
! tend to vary smoothly by storing the last level we found and searching
! up or down from this point until the right level is found. This is not as
! simple to follow, but much quicker.

last=1
DO j = 1, data_rows
  DO i = 1, row_length

    ! Searching up. Default is -1 (i.e. top level)
    level_below(i,j) = -1
    IF (p_at_data(i,j,last) >= desired_p) THEN
      DO k = last+1, data_levels
        IF (p_at_data(i,j,k) < desired_p) THEN
          level_below(i,j) = k-1
          EXIT
        END IF
      END DO

    ELSE

      ! Searching down. Default is 0 (i.e. bottom level)
      level_below(i,j) = 0
      DO k = last-1, 1, -1
        IF (p_at_data(i,j,k) >= desired_p) THEN
          level_below(i,j) = k
          EXIT
        END IF
      END DO

    END IF

    last = MAX(level_below(i,j), 1)
  END DO
END DO


! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

DO j = 1, data_rows
  DO i = 1, row_length

    IF (level_below(i,j)  ==  -1) THEN
      data_out(i,j) = data_in(i,j,data_levels)

    ELSE IF (level_below(i,j)  ==  0) THEN
      data_out(i,j) = data_in(i,j,1)

    ELSE IF (level_below(i,j)  ==  1 .OR.                         &
             level_below(i,j)  ==  data_levels - 1                &
             .OR. interp_order  ==  interp_order_linear ) THEN
! linearly interpolate
      data_out (i,j) = ( (desired_p -                             &
                          p_at_data(i,j,level_below(i,j)) )       &
                          * data_in (i,j,level_below(i,j)+1)      &
                         -(desired_p -                            &
                           p_at_data(i,j,level_below(i,j)+1)) *   &
                           data_in (i,j,level_below(i,j)) ) /     &
                        ( p_at_data(i,j,level_below(i,j)+1) -     &
                          p_at_data(i,j,level_below(i,j)) )

    ELSE IF (level_below(i,j)  ==  2 .OR.                         &
             level_below(i,j)  ==  data_levels - 2                &
             .OR. interp_order  ==  interp_order_cubic ) THEN
! cubicly interpolate

      p_here_minus = p_at_data(i,j,level_below(i,j)-1)
      p_here = p_at_data(i,j,level_below(i,j))
      p_here_plus = p_at_data(i,j,level_below(i,j)+1)
      p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)

      data_out (i,j) = ( (desired_p - p_here) *                   &
                                (desired_p - p_here_plus )*       &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here_minus - p_here) *         &
                                (p_here_minus - p_here_plus )*    &
                                (p_here_minus - p_here_plus2 ) ) *&
                              data_in (i,j,level_below(i,j)-1) +  &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here_plus )*       &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here - p_here_minus) *         &
                                (p_here - p_here_plus )*          &
                                (p_here - p_here_plus2 ) ) *      &
                              data_in (i,j,level_below(i,j)) +    &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here )*            &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here_plus - p_here_minus) *    &
                                (p_here_plus - p_here )*          &
                                (p_here_plus - p_here_plus2 ) ) * &
                              data_in (i,j,level_below(i,j)+1) +  &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here )*            &
                                (desired_p - p_here_plus ) ) /    &
                              ( (p_here_plus2 - p_here_minus) *   &
                                (p_here_plus2 - p_here )*         &
                                (p_here_plus2 - p_here_plus ) ) * &
                              data_in (i,j,level_below(i,j)+2)


    ELSE
! interpolate quinticly

      p_here_minus2 = p_at_data(i,j,level_below(i,j)-2)
      p_here_minus = p_at_data(i,j,level_below(i,j)-1)
      p_here = p_at_data(i,j,level_below(i,j))
      p_here_plus = p_at_data(i,j,level_below(i,j)+1)
      p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)
      p_here_plus3 = p_at_data(i,j,level_below(i,j)+3)

      data_out (i,j) = ((desired_p - p_here_minus) *              &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_minus2 - p_here_minus) *    &
                              (p_here_minus2 - p_here )*          &
                              (p_here_minus2 - p_here_plus )*     &
                              (p_here_minus2 - p_here_plus2 )*    &
                              (p_here_minus2 - p_here_plus3 ) ) * &
                              data_in (i,j,level_below(i,j)-2) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_minus - p_here_minus2) *    &
                              (p_here_minus - p_here )*           &
                              (p_here_minus - p_here_plus )*      &
                              (p_here_minus - p_here_plus2 )*     &
                              (p_here_minus - p_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)-1) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here - p_here_minus2) *          &
                              (p_here - p_here_minus )*           &
                              (p_here - p_here_plus )*            &
                              (p_here - p_here_plus2 )*           &
                              (p_here - p_here_plus3 ) ) *        &
                              data_in (i,j,level_below(i,j)) +    &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_plus - p_here_minus2) *     &
                              (p_here_plus - p_here_minus )*      &
                              (p_here_plus - p_here )*            &
                              (p_here_plus - p_here_plus2 )*      &
                              (p_here_plus - p_here_plus3 ) ) *   &
                              data_in (i,j,level_below(i,j)+1) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_plus2 - p_here_minus2) *    &
                              (p_here_plus2 - p_here_minus )*     &
                              (p_here_plus2 - p_here )*           &
                              (p_here_plus2 - p_here_plus )*      &
                              (p_here_plus2 - p_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)+2) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 ))/       &
                            ( (p_here_plus3 - p_here_minus2) *    &
                              (p_here_plus3 - p_here_minus )*     &
                              (p_here_plus3 - p_here )*           &
                              (p_here_plus3 - p_here_plus )*      &
                              (p_here_plus3 - p_here_plus2 ) ) *  &
                              data_in (i,j,level_below(i,j)+3)

    END IF

  END DO
END DO

! end of routine

IF (lhook) CALL dr_hook('VERT_INTERP2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vert_interp2
END MODULE vert_interp2_mod
