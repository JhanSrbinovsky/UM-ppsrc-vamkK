! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE setarea_mod

IMPLICIT NONE

CONTAINS



      SUBROUTINE setarea(nx, ny, area, jmax, offset_ny)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     set area [m^2] of each grid box
!
      INTEGER, INTENT(IN) :: nx, ny, jmax, offset_ny
      REAL, INTENT(OUT) :: area(nx, ny)
      REAL atmp
      INTEGER j_offset
      REAL arealat1
      INTEGER i, j 
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      j_offset = offset_ny 


      IF (lhook) CALL dr_hook('SETAREA',zhook_in,zhook_handle)
      DO j = 1, jmax
        IF ((j+j_offset) <= 90) THEN
! DEPENDS ON: arealat1
          atmp = arealat1(INT(ABS(91.0-(j+j_offset))))*1.0E06
        ELSE
! DEPENDS ON: arealat1
          atmp = arealat1(INT(ABS((j+j_offset)-90.0)))*1.0E06
        END IF

        DO i = 1, nx
          area(i,j) = atmp
        END DO
      END DO
      IF (lhook) CALL dr_hook('SETAREA',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE setarea

! ******************************COPYRIGHT******************************
END MODULE setarea_mod
