! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
MODULE parcons_mod

IMPLICIT NONE
!*    *PARAMETER* OF GLOBAL CONSTANTS.
!
      REAL, PARAMETER :: G = 9.806
      REAL, PARAMETER :: PI = 3.14159265358978
      REAL, PARAMETER :: CIRC = 40000000.0
      REAL, PARAMETER :: ZPI  = 2.0*PI
      REAL, PARAMETER :: RAD  = PI/180.0
      REAL, PARAMETER :: DEG  = 180.0/PI
      REAL, PARAMETER :: R    = CIRC/ZPI
!
!*     VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     *G*         REAL      ACCELERATION OF GRAVITY.
!     *PI*        REAL      PI.
!     *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
!     *RAD*       REAL      PI / 180.
!     *DEG*       REAL      180. / PI.
!     *ZPI*       REAL      2. * PI.
!     *R*         REAL      EARTH RADIUS        (METRES).
!

END MODULE parcons_mod
