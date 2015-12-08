! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

 REAL FUNCTION eps(m,n)

! use recurrence relation for normalized associated Legendre
! functions i.e.
! eps(m,n)*alf(m,n,x)= x*alf(m,n-1,x)-eps(m,n-1)*alf(m,n-2,x)
! where eps(m,n)=  sqrt((n*n-m*m)/(4*n*n-1.0))

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

 IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

 INTEGER :: m,n

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('EPS',zhook_in,zhook_handle)
 eps= SQRT((n*n-m*m)/(4*n*n-1.0))

 IF (lhook) CALL dr_hook('EPS',zhook_out,zhook_handle)
 RETURN
 END FUNCTION eps

