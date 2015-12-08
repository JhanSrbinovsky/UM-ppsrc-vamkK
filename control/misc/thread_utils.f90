! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owner's HTML page
! This file belongs in section: C70

MODULE thread_utils
  USE um_types, ONLY : integer64
!$ USE OMP_LIB




  IMPLICIT NONE

  TYPE threadLock
    LOGICAL                          :: inUse = .FALSE.
!$  INTEGER(KIND=omp_lock_kind)      :: lock
    INTEGER(KIND=integer64)          :: owner = -1
  END TYPE threadLock

  TYPE(threadLock),                                                            &
      SAVE, POINTER, PRIVATE         :: locks(:) => NULL()
  INTEGER, PARAMETER, PRIVATE        :: chunk=50
  INTEGER, PARAMETER, PRIVATE        :: noOwner=-1
  INTEGER                            :: numLocks

CONTAINS

  FUNCTION newLock()                                                           &



    RESULT(R)

    IMPLICIT NONE
    TYPE(threadLock), POINTER :: temp(:)
    INTEGER(KIND=integer64)   :: R

    IF ( .NOT.ASSOCIATED(locks) ) ALLOCATE(locks(chunk))

    DO R=1,SIZE(locks)
      IF (.NOT.locks(R)%inUse) EXIT
    END DO
    
    IF (R==SIZE(locks)+1) THEN
      ! Increase storage
      ALLOCATE (temp (SIZE(locks) + chunk ))
      temp(1:SIZE(locks))=locks(1:SIZE(locks))
      DEALLOCATE (locks)
      locks=>temp
      NULLIFY (temp)
    END IF

!$  CALL omp_init_lock( locks(R)%lock )
    locks(R)%inUse=.TRUE.
    locks(R)%owner=noOwner

  END FUNCTION newLock

  FUNCTION releaseLock(l)                                                      &



    RESULT(R)

    IMPLICIT NONE
    INTEGER(KIND=integer64), INTENT(IN) :: l
    INTEGER(KIND=integer64)             :: R

!$  CALL omp_destroy_lock( locks(l)%lock)
    locks(l)%inUse=.FALSE.
    locks(l)%owner=noOwner

    R=0

  END FUNCTION releaseLock

  FUNCTION Lock(l)                                                             &



    RESULT(R)
    IMPLICIT NONE
    INTEGER(KIND=integer64), INTENT(IN) :: l
    INTEGER(KIND=integer64)             :: tid
    INTEGER(KIND=integer64)             :: R

    tid=0
!$  tid=omp_get_thread_num()   
!$  CALL omp_set_lock(locks(l)%lock)
    locks(l)%owner=tid

    R = 0

  END FUNCTION Lock

  FUNCTION Unlock(l)                                                           &



    RESULT(R)

    IMPLICIT NONE
    INTEGER(KIND=integer64), INTENT(IN) :: l
    INTEGER(KIND=integer64)             :: R
    INTEGER(KIND=integer64)             :: tid

    tid=0
!$  tid=omp_get_thread_num()
!$  CALL omp_unset_lock(locks(l)%lock)
    locks(l)%owner=noOwner
    R = 0
  END FUNCTION Unlock

  FUNCTION threadFlush()                                                       &



    RESULT(R)
    IMPLICIT NONE
    INTEGER(KIND=integer64)             :: R
!$OMP FLUSH
    R = 0
  END FUNCTION threadFlush

  FUNCTION threadID()                                                          &



    RESULT(R)
    IMPLICIT NONE
    INTEGER(KIND=integer64)             :: R
    R = 0
!$  R = omp_get_thread_num()
  END FUNCTION threadID

END MODULE thread_utils
