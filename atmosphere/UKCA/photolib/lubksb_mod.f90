! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE lubksb_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     For matrix inversion

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CONTAINS
SUBROUTINE lubksb(a,n,np,indx,b)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: np
INTEGER, INTENT(IN) :: indx(n)
REAL, INTENT(IN)    :: a(np,np)
REAL, INTENT(INOUT) :: b(n)

! Local variables
INTEGER :: i
INTEGER :: ii
INTEGER :: j
INTEGER :: ll

REAL :: sum

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('LUBKSB',zhook_in,zhook_handle)
ii=0
DO i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  IF (ii /= 0) THEN
    DO j=ii,i-1
      sum=sum-a(i,j)*b(j)
    END DO
  ELSE IF (sum /= 0.) THEN
    ii=i
  END IF
  b(i)=sum
END DO
DO i=n,1,-1
  sum=b(i)
  IF (i < n) THEN
    DO j=i+1,n
      sum=sum-a(i,j)*b(j)
    END DO
  END IF
  b(i)=sum/a(i,i)
END DO
IF (lhook) CALL dr_hook('LUBKSB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lubksb
END MODULE lubksb_mod
