! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE setlen_mod

IMPLICIT NONE

CONTAINS




      SUBROUTINE setlen (nx, ny, igrcn, inextx, inexty, rlen            &
           , jmax, rmiss, offset_nx, offset_ny)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!     length from (i, j) to the destination in [m]
!
!     river mouth : distance to 1 grid north
!     sea         : 0.0
!
      INTEGER nx, ny, jmax, offset_nx, offset_ny
      INTEGER inextx(nx, ny), inexty(nx, ny), igrcn(nx, ny), i, j
      REAL getlat0, getlon0, rx, ry, rx2, ry2
      REAL rlen(nx, ny), rmin, rmax, givelen, rmiss
      INTEGER i_offset, j_offset

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('SETLEN',zhook_in,zhook_handle)
      rx2 = 0.0
      ry2 = 0.0
      rmin = 1.0E10
      rmax = 0.0

      j_offset = offset_ny
      i_offset = offset_nx

!
      DO j = 1, jmax
! DEPENDS ON: getlat0
        ry = getlat0(j+j_offset, 180)
        DO i = 1, nx
! DEPENDS ON: getlon0
          rx = getlon0(i+i_offset, 360)
!
          IF ((igrcn(i,j) >= 1).AND.(igrcn(i,j) <= 8)) THEN
! DEPENDS ON: getlon0
            rx2 = getlon0(inextx(i,j)+i_offset, 360)
! DEPENDS ON: getlat0
            ry2 = getlat0(inexty(i,j)+j_offset, 180)
! DEPENDS ON: givelen
            rlen(i, j) = givelen(rx, ry, rx2, ry2) * 1000.0
          ELSE IF (igrcn(i,j) == 9) THEN
! DEPENDS ON: getlat0
            ry2 = getlat0((j+j_offset)+1, 180)
! DEPENDS ON: givelen
            rlen(i, j) = givelen(rx, ry, rx, ry2) * 1000.0
          ELSE IF (igrcn(i,j) == 10) THEN
! DEPENDS ON: getlat0
            ry2 = getlat0((j+j_offset)+1, 180)
! DEPENDS ON: givelen
            rlen(i, j) = givelen(rx, ry, rx, ry2) * 1000.0
          ELSE
            rlen(i,j) = 0.0
          END IF
          IF (igrcn(i,j) /= rmiss) THEN
            IF (rlen(i,j) <  rmin) rmin = rlen(i,j)
            IF (rlen(i,j) >  rmax) rmax = rlen(i,j)
          END IF
!
        END DO
      END DO

      WRITE (6,'(A,F10.2,2X,F8.2)') 'setlen: proc local rmax, rmin = ', &
           rmax, rmin
      IF (lhook) CALL dr_hook('SETLEN',zhook_out,zhook_handle)
      RETURN
!
      END SUBROUTINE setlen
END MODULE setlen_mod
