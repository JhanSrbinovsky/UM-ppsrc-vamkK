! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
!
!
MODULE mmd2kgs_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE mmd2kgs(rin, igrcn, din, area, nx, ny, rmiss, jmax)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     convert rin [mm/day] fo din [kg] using area [m^2]
!     [mm/day] x [m^2] = [10^(-3) m^3/day]
!                       ===> [kg/day] / (3600*24) --> [kg/s]
!
      INTEGER nx, ny, i, j, igrcn(nx, ny), jmax
      REAL rin(nx, ny)
      REAL din(nx, ny), area(nx, ny), rc, rmiss

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      DATA rc/86400.0/

      IF (lhook) CALL dr_hook('MMD2KGS',zhook_in,zhook_handle)

      DO j = 1, jmax
        DO i = 1, nx

           IF ((rin(i,j) /= rmiss)) THEN
            din(i,j) = rin(i,j) * area(i,j) / rc
          ELSE
            din(i,j) = 0.0
          END IF
        END DO
      END DO
      IF (lhook) CALL dr_hook('MMD2KGS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE mmd2kgs
END MODULE mmd2kgs_mod
