! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE cp2_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE cp2(rd1, rd2, nx, ny, jmax)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     copy file
!
      integer nx, ny, i, j, jmax
      real rd1(nx, ny), rd2(nx, ny)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CP2',zhook_in,zhook_handle)
      do j = 1, jmax
        do i = 1, nx
          rd2(i,j) = rd1(i,j)
        end do
      end do
      IF (lhook) CALL dr_hook('CP2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE cp2
END MODULE cp2_mod
