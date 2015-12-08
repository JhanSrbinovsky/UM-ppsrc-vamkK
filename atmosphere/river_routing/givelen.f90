! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION givelen(rx1, ry1, rx2, ry2)

!     give the length (km) between (rx1, ry1) to (rx2, ry2)
!     sphere approximation is applied.
!     see page 621 of Rika-Nenpyo (1995)
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing


      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE

      REAL givelat, givelon, giverade, rlat, dx, dy, re
      REAL rx1, rx2, ry1, ry2, dlon

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('GIVELEN',zhook_in,zhook_handle)

      dlon = abs(rx1 - rx2)
      if (dlon >= 180.0) dlon = abs(360.0 - dlon)

      if (rx1 == rx2) then
        rlat = (ry1+ry2) / 2.0
! DEPENDS ON: givelat
        givelen = givelat(rlat) * abs(ry1-ry2)
      else if (ry1 == ry2) then
        rlat = ry1
! DEPENDS ON: givelon
        givelen = givelon(rlat) * dlon

      else
        rlat = (ry1+ry2) / 2.0
! DEPENDS ON: giverade
        re = giverade(rlat)
! DEPENDS ON: givelon
        dx = givelon(rlat) * dlon / re
! DEPENDS ON: givelat
        dy = givelat(rlat) * abs(ry1 - ry2) / re
!        write (6, '(5f15.7)') dx*re, dy*re, dx, dy, re
        givelen = acos(cos(dx)*cos(dy)) * re
      end if

      IF (lhook) CALL dr_hook('GIVELEN',zhook_out,zhook_handle)

      END FUNCTION givelen
