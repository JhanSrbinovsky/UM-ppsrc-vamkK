! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Solves system of 1D tridiagonal equations

FUNCTION eg_vert_trisolve(a,b,c,d,KK,tri_idx)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

!
! Description:
!                Used in initialisation of baroclinic test case.  
!
! Method:
!  
! Tri-diagonal solve by back use of Thomas algorithm
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER, INTENT(IN)     :: KK, tri_idx
REAL                    :: eg_vert_trisolve(tri_idx:KK) 
REAL, INTENT(IN)        :: b(tri_idx:KK), c(tri_idx:KK-1), d(tri_idx:KK)
REAL, INTENT(IN)        :: a(tri_idx+1:KK) 
REAL                    :: d_new(tri_idx:KK)
REAL                    :: c_new(tri_idx:KK-1)
INTEGER                 :: k
REAL                    :: denom

IF (lhook) CALL dr_hook('EG_VERT_TRISOLVE',zhook_in,zhook_handle)

c_new(tri_idx) = c(tri_idx)/b(tri_idx)
d_new(tri_idx) = d(tri_idx)/b(tri_idx)

DO k = tri_idx+1,KK-1
  denom = b(k)-c_new(k-1)*a(k)
  c_new(k) = c(k)/denom
  d_new(k) = (d(k)-d_new(k-1)*a(k))/denom
  IF (ABS(b(k)) < ABS(a(k) + c(k)) ) THEN
    WRITE(6,*) 'WARNING - Matrix is not diagonally dominant'
    WRITE(6,*) k,b(k),a(k)+c(k)
  END IF
END DO
denom = b(KK)-c_new(KK-1)*a(KK)
d_new(KK) = (d(KK)-d_new(KK-1)*a(KK))/denom

eg_vert_trisolve(KK) = d_new(KK)
DO k = KK-1,tri_idx,-1
 eg_vert_trisolve(k) = d_new(k)-c_new(k)*eg_vert_trisolve(k+1)
END DO

IF (lhook) CALL dr_hook('EG_VERT_TRISOLVE',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_vert_trisolve
