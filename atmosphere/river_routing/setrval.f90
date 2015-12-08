! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE setrval_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE setrval(nx, ny, rdim, rdat, jmax)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     set uniform value of rdat into rdim(i,j)
!
      integer nx, ny, i, j, jmax
      real rdim(nx, ny), rdat

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SETRVAL',zhook_in,zhook_handle)
      do i = 1, nx
        do j = 1, jmax
          rdim(i,j) = rdat
        end do
      end do
      IF (lhook) CALL dr_hook('SETRVAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE setrval
END MODULE setrval_mod
