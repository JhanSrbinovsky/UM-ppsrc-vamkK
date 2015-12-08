! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE setinit_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE setinit(nx, ny, area, rinit, sto, jmax)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     assume water equivalent to rinit [mm] over 1 grid box is stored
!     as in the initial stage
!     [mm] x [m^2] = [10^(-3) m^3] ===> [kg]
!
      integer nx, ny, i, j, jmax
      real area(nx, ny), sto(nx, ny), rinit

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('SETINIT',zhook_in,zhook_handle)
      do j = 1, jmax
        do i = 1, nx
          sto(i,j) = rinit * area(i,j)
        end do
      end do
      IF (lhook) CALL dr_hook('SETINIT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE setinit
END MODULE setinit_mod
