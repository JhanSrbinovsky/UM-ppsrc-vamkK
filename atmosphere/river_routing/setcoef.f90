! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE setcoef_mod

IMPLICIT NONE

CONTAINS



      SUBROUTINE setcoef(nx, ny, rlen, rvel, ratmed, rc, jmax)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     calculate coefficient 'c' in the routing model
!     basically by c = u/d, but actual d is assumed to d*ratmed
!     considering the meandering of the river.
!
!     c has dimension of [1/s]
!     sea grids: c = 0.0
!
      integer nx, ny, i, j, jmax
      real rlen(nx, ny), rvel(nx, ny), ratmed, rc(nx, ny)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('SETCOEF',zhook_in,zhook_handle)
      do i= 1, nx
        do j = 1, ny
          if (rlen(i,j) >  0.0) then
            rc(i, j) = rvel(i,j) / (rlen(i, j) * ratmed)
          else
            rc(i, j) = 0.0
          end if
        end do
      end do
      IF (lhook) CALL dr_hook('SETCOEF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE setcoef
END MODULE setcoef_mod
