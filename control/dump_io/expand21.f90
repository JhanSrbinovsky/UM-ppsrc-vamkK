! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE EXPAND21
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Dump I/O

!  Purpose: Unpacks IEEE 32-bit data into IEEE 64-bit data.


SUBROUTINE expand21(n, in, out)
!--expands input array 'in' from 32-bit into 'out' in 64-bit

!  n       the number of floating point words to convert
!  in      the input array of 32-bit numbers
!  out     the output array of 64-bit numbers

!--

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types
IMPLICIT NONE


! Argument variables
INTEGER (KIND=integer64), INTENT(IN) :: n
REAL (KIND=real32), INTENT(IN)       :: in(1:n)
REAL (KIND=real64), INTENT(OUT)      :: out(1:n)

! Local integer variables
INTEGER (KIND=integer64) :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('EXPAND21',zhook_in,zhook_handle)
DO i = 1,n
  OUT(i) = IN(i)
END DO

IF (lhook) CALL dr_hook('EXPAND21',zhook_out,zhook_handle)
RETURN
END SUBROUTINE expand21
