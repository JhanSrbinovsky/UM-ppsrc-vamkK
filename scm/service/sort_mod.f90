!
! MODULE SCM_UTILS--------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************

MODULE sort_mod

!-------------------------------------------------------------------------------
! Description:
!   SCM Sort : Sorts data points (x,y) to be monotonically increasing in x
!
! Method:
!   Sorting is done using bubble sorting on x
!-------------------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sort(x,y)

    USE um_types
    USE parkind1, ONLY: jpim, jprb
    USE yomhook,  ONLY: lhook, dr_hook

    IMPLICIT NONE

    REAL, INTENT(InOut) :: x(:),y(:)
    REAL    :: temp_x, temp_y
    INTEGER :: j, k, num
    LOGICAL :: sorted

    ! Dr Hook
    !==============================
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('SORT',zhook_in,zhook_handle)
    !==============================

    num    = SIZE(x)
    sorted = .FALSE.
    k=0

    DO WHILE (.NOT. sorted)

      sorted = .TRUE.
      k = k+1

      DO j=1, num-k
        IF ( x(j) > x(j+1) ) THEN

          temp_x = x(j)
          temp_y = y(j)

          x(j) = x(j+1)
          y(j) = y(j+1)

          x(j+1) = temp_x
          y(j+1) = temp_y
          sorted = .FALSE.

        END IF
      END DO
    END DO

    IF (lhook) CALL dr_hook('SORT',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE sort

!-----------------------------------------------------------------------------
END MODULE sort_mod
