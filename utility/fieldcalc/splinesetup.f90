! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height

!=======================================================================

SUBROUTINE SplineSetup (n, x, y, b, c, d)

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.  Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: n
REAL,    INTENT(IN) :: x(n), y(n)
REAL,    INTENT(OUT) :: b(n), c(n), d(n)

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SplineSetup"

! Local variables:
INTEGER :: i
REAL :: t
REAL :: deltax(n-1), deltay(n-1), dydx(n-1)
!-----------------------------------------------------------------------
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!    for  x(i) < x < x(i+1)
!  input..
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!  output..
!    b, c, d  = arrays of spline coefficients as defined above.
!  using  p  to denote differentiation,
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!  the accompanying function subprogram SplineEval can be used
!  to evaluate the spline.
!-----------------------------------------------------------------------

IF ( n >= 2 ) THEN

  deltax(1:n-1) = x(2:n) - x(1:n-1)
  deltay(1:n-1) = y(2:n) - y(1:n-1)
  dydx  (1:n-1) = deltay(1:n-1)/deltax(1:n-1)

  IF ( n == 2 ) THEN
    b(1:2) = dydx(1)
    c(1:2) = 0.
    d(1:2) = 0.
  ELSE

   !--------------------------------------------------------------------
   ! Tridiagonal system:  b = diagonal, d = offdiagonal, c = RH side.
   ! Third derivs at x(1) and x(n) obtained from divided differences
   !--------------------------------------------------------------------
    b(1)     = -deltax(1)
    b(2:n-1) =  deltax(1:n-2) + deltax(2:n-1)
    b(n)     = -deltax(n-1)

    c(2:n-1) =  dydx(2:n-1)-dydx(1:n-2)
    IF ( n /= 3 ) THEN
      c(1) = c(  3)/b(  3) - c(  2)/b(  2)
      c(n) = c(n-2)/b(n-2) - c(n-1)/b(n-1)
      c(1) = c(1)*b(1)**2/(b(  3)-b(1))
      c(n) = c(n)*b(n)**2/(b(n-2)-b(n))
    ELSE
      c(1) = 0.
      c(n) = 0.
    ENDIF
    b(2:n-1) = 2.*b(2:n-1)
    !----------------------------------------------  forward elimination
    DO i = 2, n
      t=deltax(i-1)/b(i-1)
      b(i) = b(i) - t*deltax(i-1)
      c(i) = c(i) - t*c(i-1)
    ENDDO
    !------------------------------------------------  back substitution
    c(n) = c(n)/b(n)
    DO i = n-1, 1, -1
      c(i) = (c(i) - deltax(i)*c(i+1))/b(i) !2nd derivative/ 6
    ENDDO
!------------------------------ c(i) is now the sigma(i) of the text
!---------------------------------------- compute polynomial coeff.s
!-----------------------------------------------------------------------
    d(1:n-1) =(c(2:n)-c(1:n-1))/deltax(1:n-1) !3rd derivative/ 6
    d(n)     = d(n-1)                         !

    b(1:n-1) = dydx(1:n-1) - deltax(1:n-1)*(c(2:n) + 2.*c(1:n-1))
    b(n)     = dydx(n-1)   + deltax(n-1)  *(c(n-1) + 2.*c(n)    )

    c(1:n)   = 3.*c(1:n)  ! 2nd derivative/2

  ENDIF
ENDIF

END SUBROUTINE SplineSetup

!=======================================================================
