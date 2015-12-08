! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_errfunc---------------------------------------------
!
!  Purpose: To fastly calculate values of error function in the MY
!           model.
!           The calculation is based on the expansion up to
!           the 13th order.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_errfunc(nn, x, y)

  USE conversions_mod, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nn       ! size of array

  REAL, INTENT(IN)    :: x(nn)    ! input array

  REAL, INTENT(OUT)   :: y(nn)    ! output array

! Local Variables
  INTEGER             :: i        ! Loop index

  REAL, SAVE ::                                                         &
     c01,                                                               &
          ! expansion coefficient of x
     c03,                                                               &
          ! expansion coefficient of x**3
     c05,                                                               &
          ! expansion coefficient of x**5
     c07,                                                               &
          ! expansion coefficient of x**7
     c09,                                                               &
          ! expansion coefficient of x**9
     c11,                                                               &
          ! expansion coefficient of x**11
     c13,                                                               &
          ! expansion coefficient of x**13
     factor
          ! common factor to all the coefficients

  REAL ::                                                               &
     x02,                                                               &
         ! x powered by 2
     x04,                                                               &
         ! x powered by 4
     x06,                                                               &
         ! x powered by 6
     x08,                                                               &
         ! x powered by 8
     x10,                                                               &
         ! x powered by 10
     x12
         ! x powered by 12

  LOGICAL, SAVE       :: first = .TRUE.
                                  ! flag to indication first run

  REAL, PARAMETER ::                                                    &
     erfmax = 1.0
         ! upper limit of the value to avoid it outside domain

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_ERRFUNC',zhook_in,zhook_handle)

  IF (first) THEN
    factor =  2.0 / SQRT(pi)
    c01 = factor * 1.0
    c03 = factor * 1.0 /    3.0
    c05 = factor * 1.0 /   10.0
    c07 = factor * 1.0 /   42.0
    c09 = factor * 1.0 /  216.0
    c11 = factor * 1.0 / 1320.0
    c13 = factor * 1.0 / 9360.0
    first = .FALSE.
  END IF
  DO i = 1, nn
    x02 = x(i) * x(i)
    x04 = x02 * x02
    x06 = x04 * x02
    x08 = x06 * x02
    x10 = x08 * x02
    x12 = x10 * x02
    y(i) = x(i) * (                                                     &
          + c01                                                         &
          - c03 * x02                                                   &
          + c05 * x04                                                   &
          - c07 * x06                                                   &
          + c09 * x08                                                   &
          - c11 * x10                                                   &
          + c13 * x12)
    IF(x(i) > 0) THEN
      y(i) = MIN(y(i), erfmax)
    ELSE
      y(i) = MAX(y(i), -erfmax)
    END IF
  END DO
  IF (lhook) CALL dr_hook('MYM_ERRFUNC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_errfunc
