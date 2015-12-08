! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing

      REAL FUNCTION givelat(rlat)

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE

      REAL rlat
      ! These should really be from the constants module.
      REAL, PARAMETER :: e2 = 0.006694470
      REAL, PARAMETER :: pi = 3.14159265358979
      REAL, PARAMETER :: ra = 6378.136

      REAL sind
      EXTERNAL sind

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('GIVELAT',zhook_in,zhook_handle)

      givelat = pi / 180.0 * ra * ( 1 - e2)                             &
! DEPENDS ON: sind
     &     / SQRT((1 - e2*sind(rlat)*sind(rlat))**3)

      IF (lhook) CALL dr_hook('GIVELAT',zhook_out,zhook_handle)

      END FUNCTION givelat
