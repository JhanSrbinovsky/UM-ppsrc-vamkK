! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Chebyshev polynomial function

REAL FUNCTION ChebPoly (x, n) ! in, in

! Method:
!
!   Use a recursion relation.
!
!   The Chebyshev polynomial T_n of order n is given by
!
!     T_0(x) = 1
!     T_1(x) = x
!
!     T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x),  (n > 1)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Function arguments:

REAL,    INTENT(IN) :: x
INTEGER, INTENT(IN) :: n  ! Order (non-negative).

! Local variables:

INTEGER :: i

REAL :: LastButOne_val
REAL :: Last_val
REAL :: Latest_val

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header --------------------------------------------------------------

IF (lhook) CALL dr_hook('CHEBPOLY',zhook_in,zhook_handle)

IF (n == 0) ChebPoly = 1.0
IF (n == 1) ChebPoly = x

IF (n > 1) THEN

  LastButOne_val = 1.0
  Last_val       = x

  DO i = 2, n
    Latest_val     = 2.0 * x * Last_val - LastButOne_val
    LastButOne_val = Last_val
    Last_val       = Latest_val
  END DO

  ChebPoly = Latest_val

END IF

IF (lhook) CALL dr_hook('CHEBPOLY',zhook_out,zhook_handle)
RETURN

END FUNCTION ChebPoly
