! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     computes the normalized associated Legendre function
!

REAL FUNCTION alf(m,n,x)

USE conversions_mod, ONLY: pi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER :: m,n,nn,i
REAL    :: x,somx2,fact,pmm,pmmp1,pnn,eps

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ALF',zhook_in,zhook_handle)
IF (m < 0 .OR. m > n .OR. ABS(x) > 1.0) PRINT*,'bad arguments'
!     calculate alf(m,m,x) (pmm) directly
pmm= 1.0
IF (m > 0) THEN
 somx2= SQRT( (1.0-x)*(1.0+x) )
 fact= 1.0
 DO i=1,m
   pmm= pmm*somx2*SQRT(fact/(fact+1))
   fact= fact +2.0
 END DO
END IF
pmm= SQRT(m+0.5)*pmm

! use recurrence relation for normalized associated Legendre
! functions i.e.
! eps(m,n)*alf(m,n,x)= x*alf(m,n-1,x)-eps(m,n-1)*alf(m,n-2,x)
! where eps(m,n)=  sqrt((n*n-m*m)/(4*n*n-1.0))


IF (n == m) THEN
 alf= pmm
ELSE
!DEPENDS ON:eps
 pmmp1= x*pmm/eps(m,m+1)
 IF (n == m+1) THEN
  alf= pmmp1
 ELSE
  DO nn= m+2,n
!DEPENDS ON:eps
   pnn= (x*pmmp1-eps(m,nn-1)*pmm)/eps(m,nn)
   pmm= pmmp1
   pmmp1= pnn
  END DO
  alf= pnn
 END IF
END IF

!   enforce spherical harmonic normalization where
!   longitudinal integral gives 2*pi factor
alf= alf/SQRT(2*pi)

IF (lhook) CALL dr_hook('ALF',zhook_out,zhook_handle)
RETURN
END FUNCTION alf


