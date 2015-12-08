! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Gets the cpu time from the system

! Function Interface:
REAL FUNCTION get_cpu_time()


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!   There used to be no Fortran standard for calculating CPU time.
!   This routine was therefore an in interface to propriatory library
!   functions, but now uses the Fortran 95 intrinsic.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Timer

! Code Description:
!   Language: Fortran 95

REAL :: time

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('GET_CPU_TIME',zhook_in,zhook_handle)
CALL CPU_TIME(time)
get_cpu_time = time

IF (lhook) CALL dr_hook('GET_CPU_TIME',zhook_out,zhook_handle)
RETURN
END FUNCTION get_cpu_time
