! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Solves the equation A.X=Y where A is a tridiagonal matrix
!

SUBROUTINE tridiag(a,b,c,r,u,n)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------
! Description:
!   Solves the equation A.X=Y where A is a tridiagonal matrix 
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  n                      ! size of vectors


REAL, INTENT(IN)    :: & 
  a(n)                 & ! Components of tridiagonal matrix
 ,b(n)                 & !
 ,c(n)                 & !
 ,r(n)                   ! RHS of linear equation

REAL, INTENT(OUT) ::   &
  u(n)                   ! Solution vector

! Local variables

INTEGER ::       &
  j                ! Loop counter

REAL ::          &
  gam(n)         & ! Dimension of gam, should be the same as n
 ,bet    

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG',zhook_in,zhook_handle)

!------------------------------------------------------------------------

bet=b(1)
u(1)=r(1)/bet
DO j=2,n
  gam(j)=c(j-1)/bet
  bet=b(j)-a(j)*gam(j)
  u(j)=(r(j)-a(j)*u(j-1))/bet
END DO
DO j=n-1,1,-1
  u(j)=u(j)-gam(j+1)*u(j+1)
END DO
!------------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDIAG',zhook_out,zhook_handle)

!------------------------------------------------------------------------

RETURN
END SUBROUTINE tridiag
