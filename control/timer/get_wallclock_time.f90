! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------------

!  Gets the elapsed time from the system

! Function Interface:
REAL FUNCTION get_wallclock_time()

!$ USE omp_lib

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!   The system function SYSTEM_CLOCK is used to return the numbers of
!   seconds which have elapsed.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Timer

! Code Description:
!   Language: FORTRAN 90


!- End of header

! Local variables
INTEGER :: tid, tmax          ! threading variables

INTEGER, ALLOCATABLE, SAVE :: start_count(:)
INTEGER, ALLOCATABLE, SAVE :: old_count(:)
REAL, ALLOCATABLE, SAVE    :: rollover(:)
REAL, ALLOCATABLE, SAVE    :: oneover_count_rate(:)

INTEGER       :: my_count, count_rate, count_max, elapsed_count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle







! Intrinsic procedures called:
INTRINSIC SYSTEM_CLOCK




IF (lhook) CALL dr_hook('GET_WALLCLOCK_TIME',zhook_in,zhook_handle)


tid = 0
!$ tid = omp_get_thread_num()
!$OMP CRITICAL(gwt_alloc)
IF (.NOT. ALLOCATED(start_count)) THEN
  tmax = 1
!$ tmax = omp_get_max_threads()
  ALLOCATE(start_count(0:tmax-1))
  start_count(:) = -1
  ALLOCATE(old_count(0:tmax-1))
  old_count(:) = 0
  ALLOCATE(rollover(0:tmax-1))
  rollover(:) = 0.0
  ALLOCATE(oneover_count_rate(0:tmax-1))
  oneover_count_rate(:) = 0.0
END IF
!$OMP END CRITICAL(gwt_alloc)

CALL SYSTEM_CLOCK(COUNT=my_count, COUNT_RATE=count_rate, &
                  COUNT_MAX=count_max)

IF ((old_count(tid)  <   start_count(tid)) .AND.                       &
    ((my_count <   old_count(tid)) .OR.                                &
     (my_count >   start_count(tid)))) THEN
  IF (count_rate /= 0) THEN
    rollover(tid)=rollover(tid)+(REAL(count_max)/REAL(count_rate))
  END IF
END IF

IF (start_count(tid)  ==  -1) THEN
  start_count(tid) = my_count
  IF (count_rate /= 0) THEN
    oneover_count_rate(tid)=1.0/REAL(count_rate)
  ELSE
    oneover_count_rate(tid) = 1.0
  END IF
END IF

elapsed_count = my_count - start_count(tid)

IF (elapsed_count  <   0) elapsed_count = elapsed_count + count_max

get_wallclock_time = rollover(tid)+                                    &
                    (REAL(elapsed_count)*oneover_count_rate(tid))
old_count(tid) = my_count

IF (lhook) CALL dr_hook('GET_WALLCLOCK_TIME',zhook_out,zhook_handle)
RETURN
END FUNCTION get_wallclock_time
