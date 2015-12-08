! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   solves tridiagonal matrix -  convection scheme

SUBROUTINE tridiag_all(n,nvec,nvec_len,a,b,c,r,u)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!------------------------------------------------------------------------
! Description:
!  Solves the equations A.X = Y,  where A is a tridiagonal matrix
!          for several matrices A at once.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  n                    & ! maximum size of vectors X and Y
 ,nvec                 & ! Number of vectors/matrices to solve
 ,nvec_len(nvec)         ! length of each vector.

REAL, INTENT(IN) :: &
  a(nvec,n)         & ! Components of tridiagonal matrix
 ,b(nvec,n)         & ! Components of tridiagonal matrix
 ,c(nvec,n)         & ! Components of tridiagonal matrix
 ,r(nvec,n)           ! vector Y, R.H.S of linear equation

REAL, INTENT(OUT) :: &
  u(nvec,n)            ! solution vectors


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER ::     &
  i,j            ! loop counter

REAL ::        &
  gam(nvec,n)  & ! work array
 ,bet(nvec)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG_ALL',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
DO i=1,nvec
  bet(i) = b(i,1)
  u(i,1) = r(i,1)/bet(i)
END DO

DO j=2,n
  DO i=1, nvec
    IF (j <= nvec_len(i)) THEN
      gam(i,j) = c(i,j-1)/bet(i)
      bet(i)   = b(i,j) - a(i,j)*gam(i,j)
!for info !    if (bet(i) == 0.0) stop   ! was in original code
      u(i,j)   = (r(i,j) - a(i,j)*u(i,j-1))/bet(i)
    END IF
  END DO
END DO

DO j=n-1,1,-1
  DO i=1, nvec
    IF (j <= (nvec_len(i)-1)) THEN
      u(i,j) = u(i,j) - gam(i,j+1)*u(i,j+1)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG_ALL',zhook_out,zhook_handle)

!-----------------------------------------------------------------------
RETURN
END SUBROUTINE tridiag_all
