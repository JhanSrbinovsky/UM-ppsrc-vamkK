! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE pole_bearing(row_length, rows, lat_rot_NP, &
    long_rot_NP, bear_rot_NP)

  USE trignometric_mod, ONLY: true_latitude, true_longitude
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Description:
!   Calculates bearing of Grid North (from True North) for gridpoints
!   in a Local Area Model rotated grid.
!
! Method:
!   Spherical trig.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows         ! grid size
  REAL, INTENT(IN) ::  &
       lat_rot_NP,     &        ! Real latitude and longitude
       long_rot_NP              !  of 'pseudo' N pole in radians.

  REAL, DIMENSION(row_length,rows), INTENT(OUT) :: &
       bear_rot_NP              ! Bearing of 'pseudo' N pole (rads)

! Local variables
  REAL :: long_diff(row_length,rows)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! End of header

  IF (lhook) CALL dr_hook('POLE_BEARING',zhook_in,zhook_handle)

! Calculate bearing of Grid North
  long_diff = true_longitude - long_rot_NP

  bear_rot_NP = atan2( -cos(lat_rot_NP)*sin(long_diff),  &
        cos(true_latitude)*sin(lat_rot_NP) -                  &
        sin(true_latitude)*cos(lat_rot_NP)*cos(long_diff) )

  IF (lhook) CALL dr_hook('POLE_BEARING',zhook_out,zhook_handle)

END SUBROUTINE pole_bearing
