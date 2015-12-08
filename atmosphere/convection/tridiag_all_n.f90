! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  solves tridiagonal matrix -  convection scheme
!
SUBROUTINE tridiag_all_n(N,nvec,A,B,C,R,U)
!
! Purpose: Solves the equations A.X = Y,  where A is a tridiagnol matrix
!
!          for several matrices A at once.
!          Version where all vetcors same length.
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) :: &
  N                    & ! maximum size of vectors X and Y
 ,nvec                   ! Number of vectors/matrices to solve


REAL, INTENT(IN) ::    &
  A(nvec,N)            & ! Components of tridiagonal matrix
 ,B(nvec,N)            & ! Components of tridiagonal matrix
 ,C(nvec,N)            & ! Components of tridiagonal matrix
 ,R(nvec,N)              ! vector Y, R.H.S of linear equation

REAL, INTENT(OUT) ::   &
  U(nvec,N)              ! solution vectors


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: &
  i,j            ! loop counters

REAL ::         &
  gam(nvec,n)   & ! work array
 ,bet(nvec)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG_ALL_N',zhook_in,zhook_handle)

  DO i=1,nvec
    bet(i) = B(i,1)
    U(i,1) = R(i,1)/bet(i)
  END DO

  DO j=2,N
    DO i=1, nvec

      gam(i,j) = C(i,j-1)/bet(i)
      bet(i)   = B(i,j) - A(i,j)*gam(i,j)
      u(i,j)   = (R(i,j) - A(i,j)*U(i,j-1))/bet(i)

    END DO
  END DO

  DO j=N-1,1,-1
    DO i=1, nvec

      U(i,j) = U(i,j) - Gam(i,j+1)*U(i,j+1)

    END DO
  END DO

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG_ALL_N',zhook_out,zhook_handle)
RETURN
END SUBROUTINE tridiag_all_n
