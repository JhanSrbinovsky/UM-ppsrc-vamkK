! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing

      REAL FUNCTION givelon(rlat)

!     give the length (km) of 1 degree longitude  (dx)
!     at the latitude of rlat(degree)
!     see page 621 of Rika-Nenpyo (1995)

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE

      REAL rlat 
      
! These should really be from the constants module.
      REAL, PARAMETER :: e2 = 0.006694470
      REAL, PARAMETER :: pi = 3.14159265358979
      REAL, PARAMETER :: ra = 6378.136

      REAL sind, cosd
      EXTERNAL sind, cosd

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('GIVELON',zhook_in,zhook_handle)

! DEPENDS ON: cosd
      givelon = pi / 180.0 * ra * cosd(rlat)                            &
! DEPENDS ON: sind
          / SQRT(1 - e2*sind(rlat)*sind(rlat))

      IF (lhook) CALL dr_hook('GIVELON',zhook_out,zhook_handle)

      END FUNCTION givelon
