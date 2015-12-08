! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculation of horizon angle along a 1D slice

SUBROUTINE horizon_1d(orog, nd_points, n_points, horiz_limit,           &
                      sep_ang, horiz_ang1, horiz_ang2)

  USE conversions_mod, ONLY: pi
  USE earth_constants_mod, ONLY: earth_radius
  IMPLICIT NONE

! Description:
!   Determines horizon angles in the forwards and backwards directions
!
! Method:
!   Uses a method similar to that outlined in Dozier and Frew (1990
!   IEEE Trans. Geosci. Remote Sens., 28, 963-969) adapted for use
!   within a limited horizon and allowing for the curvature of the
!   Earth.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

! Subroutine arguments

  INTEGER, INTENT(IN) :: nd_points, n_points, horiz_limit

  REAL, INTENT(IN) :: orog(1-horiz_limit:nd_points+horiz_limit)
  REAL, INTENT(IN) :: sep_ang(1-horiz_limit:nd_points+horiz_limit-1)

  REAL, INTENT(OUT) :: horiz_ang1(nd_points), horiz_ang2(nd_points)

! Local variables

  INTEGER :: i,k
  INTEGER :: horiz_point(1-horiz_limit:n_points+horiz_limit)
  REAL :: ang, h_ang

! End of header

  horiz_ang1 = 0.0
  horiz_point = -horiz_limit

! First proceed forwards along the strip calculating horizon angles in
! the backwards direction
  DO i=1,n_points
    k=1
!   Set maximum horizon angle to the horizontal
    horiz_ang1(i) = pi/2.0
!   Initially each point is its own horizon
    horiz_point(i) = i
    DO
!     Loop backwards from each point until we reach the horizon limit
      IF (k > horiz_limit) EXIT
!     Only consider points with higher orography
      IF ( orog(i-k) > orog(horiz_point(i)) ) THEN
        ang = SUM(sep_ang(i-k:i-1))
!       Horizon angle considers curvature of the earth
        h_ang = ATAN2( (earth_radius+orog(i-k))*SIN(ang),               &
          orog(i-k)*COS(ang) - orog(i) - earth_radius*(1.0-COS(ang)) )
        IF (h_ang < horiz_ang1(i)) THEN
!         Record the new horizon angle and the point which forms it
          horiz_ang1(i) = h_ang
          horiz_point(i) = i-k
        END IF
!       If the horizon angle to this point is smaller (zenith to horizon)
!       than the further horizon angle viewed from this point then we
!       need not consider any further points.
        IF (i-k >= 1 .AND. h_ang < horiz_ang1(i-k)) EXIT
      END IF
!     If this point is its own horizon we can stop
      IF (horiz_point(i-k) == i-k) THEN
        EXIT
!     If this point's own horizon is within our horizon limit then we
!     can jump straight to the further horizon point
      ELSE IF (horiz_point(i-k) >= i-horiz_limit) THEN
        k=i-horiz_point(i-k)
!     Otherwise we look to the next point along
      ELSE
        k=k+1
      END IF
    END DO
  END DO

  horiz_ang2 = 0.0
  horiz_point = n_points+horiz_limit+1

! Now proceed backwards along the strip calculating horizon angles in
! the forwards direction
  DO i=n_points,1,-1
    k=1
    horiz_ang2(i) = pi/2.0
    horiz_point(i) = i
    DO
      IF (k > horiz_limit) EXIT
      IF ( orog(i+k) > orog(horiz_point(i)) ) THEN
        ang = SUM(sep_ang(i:i+k-1))
        h_ang = ATAN2( (earth_radius+orog(i+k))*SIN(ang),               &
          orog(i+k)*COS(ang) - orog(i) - earth_radius*(1.0-COS(ang)) )
        IF (h_ang < horiz_ang2(i)) THEN
          horiz_ang2(i) = h_ang
          horiz_point(i) = i+k
        END IF
        IF (i+k <= n_points .AND. h_ang < horiz_ang2(i+k)) EXIT
      END IF
      IF (horiz_point(i+k) == i+k) THEN
        EXIT
      ELSE IF (horiz_point(i+k) <= i+horiz_limit) THEN
        k=horiz_point(i+k)-i
      ELSE
        k=k+1
      END IF
    END DO
  END DO

END SUBROUTINE horizon_1d
