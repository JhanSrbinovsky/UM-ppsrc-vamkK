! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height


!=======================================================================

SUBROUTINE SplineSetupV(m, n, kmax, x, y, b, c, d, MaxWLev)

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
! 6.2     19/12/05 New subroutine. Vectorised version of existing
!                  SplineSetup subroutine. J-C Rioual (NEC)
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: m, n,kmax
REAL,    INTENT(IN) :: x(m,n,kmax), y(m,n,kmax)
REAL,    INTENT(OUT) :: b(m,n,kmax), c(m,n,kmax), d(m,n,kmax)

INTEGER, DIMENSION(m,n), INTENT(IN) :: MaxWLev

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SplineSetupV"

! Local variables:
INTEGER :: i, j, k
REAL :: t(m,n)
REAL :: deltax(m,n,kmax-1), deltay(m,n,kmax-1), dydx(m,n,kmax-1)

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

IF ( kmax >= 2 ) THEN

  do k=1, kmax-1
     do j=1, n
        do i=1, m
           if (maxwlev(i,j) /=0 ) then
              deltax(i,j,k)=x(i,j,k+1)-x(i,j,k)
              deltay(i,j,k)=y(i,j,k+1)-y(i,j,k)
              dydx(i,j,k)=deltay(i,j,k)/deltax(i,j,k)
           end if
        end do
     end do
  end do

  IF ( kmax == 2 ) THEN
    do j=1, n
        do i=1, m
           if (maxwlev(i,j) /=0 ) then
              b(i,j,1) = dydx(i,j,1)
              b(i,j,2) = dydx(i,j,1)
              c(i,j,1:2) = 0.
              d(i,j,1:2) = 0.
           end if
        end do
    end do

  ELSE
   !--------------------------------------------------------------------
   ! Tridiagonal system:  b = diagonal, d = offdiagonal, c = RH side.
   ! Third derivs at x(1) and x(n) obtained from divided differences
   !--------------------------------------------------------------------
    do j=1, n
       do i=1, m
          if (maxwlev(i,j) /=0 ) then
             b(i,j,1)     = -deltax(i,j,1)
             b(i,j,kmax)  = -deltax(i,j,kmax-1)
          end if
       end do
    end do

    do k=2, kmax-1
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
                b(i,j,k) =  deltax(i,j,k-1) + deltax(i,j,k)
                c(i,j,k) =  dydx(i,j,k)-dydx(i,j,k-1)
             end if
          end do
       end do
    end do

    do j=1, n
       do i=1, m
          if (maxwlev(i,j) /=0 ) then
             if ( kmax /= 3 ) then
                c(i,j,1) = c(i,j,3)/b(i,j,3) - c(i,j,2)/b(i,j,2)
                c(i,j,kmax) = c(i,j,kmax-2)/b(i,j,kmax-2) - &
                            & c(i,j,kmax-1)/b(i,j,kmax-1)
                c(i,j,1) = c(i,j,1)*b(i,j,1)**2/(b(i,j,3)-b(i,j,1))
                c(i,j,kmax) = c(i,j,kmax)*b(i,j,kmax)**2 &
                            & /(b(i,j,kmax-2)-b(i,j,kmax))
             else
                c(i,j,1) = 0.
                c(i,j,n) = 0.
             endif
          end if
       end do
    end do

    do k=2, kmax-1
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
                b(i,j,k) = 2.*b(i,j,k)
             end if
          end do
       end do
    end do

    !----------------------------------------------  forward elimination
    DO k = 2, kmax
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
               t(i,j)=deltax(i,j,k-1)/b(i,j,k-1)
               b(i,j,k) = b(i,j,k) - t(i,j)*deltax(i,j,k-1)
               c(i,j,k) = c(i,j,k) - t(i,j)*c(i,j,k-1)
             end if
          end do
       end do
    end do

    !------------------------------------------------  back substitution

    do j=1, n
       do i=1, m
         if (maxwlev(i,j) /=0 ) then
            c(i,j,kmax) = c(i,j,kmax)/b(i,j,kmax)
         end if
       end do
    end do

    do k = kmax-1, 1, -1
      do j=1, n
         do i=1, m
            if (maxwlev(i,j) /=0 ) then
               c(i,j,k) = (c(i,j,k) - deltax(i,j,k)*c(i,j,k+1)) &
                         & /b(i,j,k) !2nd derivative/ 6
            end if
         end do
       end do
    end do

!------------------------------ c(:,:,i) is now the sigma(i) of the text
!---------------------------------------- compute polynomial coeff.s
!-----------------------------------------------------------------------
    do k=1, kmax-1
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
                d(i,j,k) =(c(i,j,k+1)-c(i,j,k)) &
                         & /deltax(i,j,k) !3rd derivative/ 6
             end if
          end do
       end do
    end do

    do j=1, n
       do i=1, m
          if (maxwlev(i,j) /=0 ) then
             d(i,j,kmax)     = d(i,j,kmax-1)                         !
          end if
       end do
    end do

    do k=1, kmax-1
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
                b(i,j,k)=dydx(i,j,k)-deltax(i,j,k)* &
                        & (c(i,j,k+1)+2.*c(i,j,k))
             end if
          end do
        end do
    end do

    do j=1, n
       do i=1, m
          if (maxwlev(i,j) /=0 ) then
             b(i,j,kmax)=dydx(i,j,kmax-1)  +deltax(i,j,kmax-1) &
            &  *(c(i,j,kmax-1)+2.*c(i,j,kmax))
          end if
       end do
    end do

    do k=1,kmax
       do j=1, n
          do i=1, m
             if (maxwlev(i,j) /=0 ) then
               c(i,j,k)   = 3.*c(i,j,k)  ! 2nd derivative/2
             end if
          end do
       end do
    end do

  ENDIF
ENDIF

END SUBROUTINE SplineSetupV

!=======================================================================
