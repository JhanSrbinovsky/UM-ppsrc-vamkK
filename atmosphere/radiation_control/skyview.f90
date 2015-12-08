! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculation of sky-view factor

SUBROUTINE skyview(row_length, rows, bear_rot_NP,                       &
                   phi, delta_phi, delta_lambda, orog)

  USE solinc_data, ONLY:                                                &
    slope_aspect, slope_angle, horiz_ang, horiz_aspect, horiz_limit,    &
    sky, n_horiz_ang
  USE conversions_mod, ONLY: pi
  IMPLICIT NONE

! Description:
!   Determines the sky-view factor for the inclined surface.
!
! Method:
!   Uses the model orography and 16 lines of sight aligned with the model
!   grid together with the pre-determined surface slope. Sky-view factor
!   is determined using spherical trigonometry based on the view from
!   the inclined plane. This is similar to equation (9) from Lai et al.
!   2010 J. Geohys. Res., 115, D01104 and equivalent to eqn (7b) from
!   Dozier and Frew 1990 IEEE Trans. Geosci. Remote Sens., 28, 963-969. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

! Subroutine arguments

  INTEGER, INTENT(IN) :: row_length, rows

  REAL, INTENT(IN) :: bear_rot_NP(row_length,rows)

  REAL, INTENT(IN) :: phi(1-horiz_limit:rows+horiz_limit)
  REAL, INTENT(IN) :: delta_phi(1-horiz_limit:rows+horiz_limit)
  REAL, INTENT(IN) :: delta_lambda(1-horiz_limit:row_length+horiz_limit)
  REAL, INTENT(IN) :: orog(1-horiz_limit:row_length+horiz_limit,        &
                           1-horiz_limit:rows+horiz_limit)

! Local constants

  INTEGER, PARAMETER :: ip_N   = 1
  INTEGER, PARAMETER :: ip_NE  = 2
  INTEGER, PARAMETER :: ip_E   = 3
  INTEGER, PARAMETER :: ip_SE  = 4
  INTEGER, PARAMETER :: ip_S   = 5
  INTEGER, PARAMETER :: ip_SW  = 6
  INTEGER, PARAMETER :: ip_W   = 7
  INTEGER, PARAMETER :: ip_NW  = 8

  INTEGER, PARAMETER :: ip_NNE = 9
  INTEGER, PARAMETER :: ip_ENE = 10
  INTEGER, PARAMETER :: ip_ESE = 11
  INTEGER, PARAMETER :: ip_SSE = 12
  INTEGER, PARAMETER :: ip_SSW = 13
  INTEGER, PARAMETER :: ip_WSW = 14
  INTEGER, PARAMETER :: ip_WNW = 15
  INTEGER, PARAMETER :: ip_NNW = 16

  REAL, PARAMETER :: pi_over_2 = pi/2.0

! Local variables

  INTEGER :: i,ii,j,jj,k,kk,gj
  INTEGER :: nd_points

  REAL :: sinlat(1-horiz_limit:rows+horiz_limit)
  REAL :: coslat(1-horiz_limit:rows+horiz_limit)
  REAL :: cos_delta_lambda(1-horiz_limit:row_length+horiz_limit)
  REAL :: cos_two_delta_lambda(1-horiz_limit:row_length+horiz_limit)
  REAL :: max_ang

  REAL :: sep_ang(1-horiz_limit:2*MAX(row_length,rows)+horiz_limit-1)
  REAL :: orog_1d(1-horiz_limit:2*MAX(row_length,rows)+horiz_limit)
  REAL :: horiz_ang1(2*MAX(row_length,rows))
  REAL :: horiz_ang2(2*MAX(row_length,rows))

  REAL :: horiz_aspect_rot(row_length,rows,n_horiz_ang)
  REAL :: tilt_ang(row_length,rows,n_horiz_ang)
  REAL :: delta_ang(row_length,rows)

! End of header

  ALLOCATE(horiz_ang(row_length,rows,n_horiz_ang))
  ALLOCATE(horiz_aspect(row_length,rows,n_horiz_ang))

  nd_points=2*MAX(row_length,rows)

! Calculate arrays used more than once
  sinlat=SIN(phi)
  coslat=COS(phi)
  cos_delta_lambda = COS(delta_lambda)
  cos_two_delta_lambda = COS(2.0*delta_lambda)

! Max separation angle to catch errors where 1D slice wraps around grid
  max_ang = pi_over_2/REAL(horiz_limit)

! Assume grid directions approximate to great circles from the local
! bearings. Horizon angles:
!  8   1   2
!   \  |  /
!    \ | /
! 7---   ---3
!    / | !   /  |  !  6   5   4



  DO j=1,rows
    DO i=1,row_length
      delta_ang(i,j)=ATAN(delta_lambda(i)*coslat(j)/delta_phi(j))
    END DO
  END DO

  horiz_aspect_rot(:,:,ip_N)  = 0.0
  horiz_aspect_rot(:,:,ip_NE) = delta_ang
  horiz_aspect_rot(:,:,ip_E)  = pi_over_2
  horiz_aspect_rot(:,:,ip_SE) = pi-delta_ang
  horiz_aspect_rot(:,:,ip_S)  = pi
  horiz_aspect_rot(:,:,ip_SW) = pi+delta_ang
  horiz_aspect_rot(:,:,ip_W)  = pi*1.5
  horiz_aspect_rot(:,:,ip_NW) = pi*2.-delta_ang

  DO j=1,rows
    DO i=1,row_length
      delta_ang(i,j)=ATAN(0.5*delta_lambda(i)*coslat(j)/delta_phi(j))
    END DO
  END DO

  horiz_aspect_rot(:,:,ip_NNE) = delta_ang
  horiz_aspect_rot(:,:,ip_SSE) = pi-delta_ang
  horiz_aspect_rot(:,:,ip_SSW) = pi+delta_ang
  horiz_aspect_rot(:,:,ip_NNW) = pi*2.-delta_ang

  DO j=1,rows
    DO i=1,row_length
      delta_ang(i,j)=ATAN(2.0*delta_lambda(i)*coslat(j)/delta_phi(j))
    END DO
  END DO

  horiz_aspect_rot(:,:,ip_ENE) = delta_ang
  horiz_aspect_rot(:,:,ip_ESE) = pi-delta_ang
  horiz_aspect_rot(:,:,ip_WSW) = pi+delta_ang
  horiz_aspect_rot(:,:,ip_WNW) = pi*2.-delta_ang

! Add bearing of 'pseudo' N pole so aspects are relative to
! true North.
  DO k=1,n_horiz_ang
    DO j=1,rows
      DO i=1,row_length
        horiz_aspect(i,j,k) =                                           &
          MODULO(horiz_aspect_rot(i,j,k) + bear_rot_NP(i,j), pi*2.)
      END DO
    END DO
  END DO


! Calculate horizon angles along 1D slices.

! Horizon angles W and E
! ----------------------
  DO j=1,rows
    sep_ang(1-horiz_limit:row_length+horiz_limit)=delta_lambda*coslat(j)
    orog_1d(1-horiz_limit:row_length+horiz_limit)=orog(:,j)

! DEPENDS ON: horizon_1d
    CALL horizon_1d(orog_1d, nd_points, row_length, horiz_limit,        &
                    sep_ang, horiz_ang1, horiz_ang2)

    horiz_ang(:,j,ip_W) = horiz_ang1(1:row_length)
    horiz_ang(:,j,ip_E) = horiz_ang2(1:row_length)
  END DO

! Horizon angles S and N
! ----------------------
  DO i=1,row_length
    sep_ang(1-horiz_limit:rows+horiz_limit) = delta_phi
    orog_1d(1-horiz_limit:rows+horiz_limit) = orog(i,:)

    CALL horizon_1d(orog_1d, nd_points, rows, horiz_limit,              &
                    sep_ang, horiz_ang1, horiz_ang2)

    horiz_ang(i,:,ip_S) = horiz_ang1(1:rows)
    horiz_ang(i,:,ip_N) = horiz_ang2(1:rows)
  END DO


! Slices along the diagonals use a sheared array.

! Horizon angles NW and SE
! ------------------------
  DO j=1-horiz_limit,rows+horiz_limit
    i=1-horiz_limit
    orog_1d(i) = orog(i,j)
    kk=j
!   Map the diagonals onto 1D slices, wrapping around the grid to
!   minimise the number of loops required
    DO i=2-horiz_limit,row_length+horiz_limit
      k = MODULO(j-i,rows+2*horiz_limit)+1-horiz_limit
      sep_ang(i-1)=MIN(ACOS(sinlat(k)*sinlat(kk)+coslat(k)*coslat(kk)   &
                       *cos_delta_lambda(i-1) ), max_ang)
      orog_1d(i) = orog(i,k)
      kk=k
    END DO

    CALL horizon_1d(orog_1d, nd_points, row_length, horiz_limit,        &
                    sep_ang, horiz_ang1, horiz_ang2)

!   Map the slices back on to the 2D grid
    DO i=1,row_length
      k = MODULO(j-i,rows+2*horiz_limit)+1-horiz_limit
      IF (k >= 1 .AND. k <= rows) THEN
        horiz_ang(i,k,ip_NW) = horiz_ang1(i)
        horiz_ang(i,k,ip_SE) = horiz_ang2(i)
      END IF
    END DO
  END DO

! Horizon angles SW and NE
! ------------------------
  DO j=1-horiz_limit,rows+horiz_limit
    i=1-horiz_limit
    orog_1d(i) = orog(i,j)
    kk=j
    DO i=2-horiz_limit,row_length+horiz_limit
      k = MODULO(j+i-2*(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      sep_ang(i-1)=MIN(ACOS(sinlat(k)*sinlat(kk)+coslat(k)*coslat(kk)   &
                       *cos_delta_lambda(i-1) ), max_ang)
      orog_1d(i) = orog(i,k)
      kk=k
    END DO

    CALL horizon_1d(orog_1d, nd_points, row_length, horiz_limit,        &
                    sep_ang, horiz_ang1, horiz_ang2)

    DO i=1,row_length
      k = MODULO(j+i-2*(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      IF (k >= 1 .AND. k <= rows) THEN
        horiz_ang(i,k,ip_SW) = horiz_ang1(i)
        horiz_ang(i,k,ip_NE) = horiz_ang2(i)
      END IF
    END DO
  END DO


! The intermediate diagonals require a shear of 2 grid-boxes, with
! points in-between given by the mean of the 2 points bisected by the
! line of sight.

! Horizon angles NNW and SSE
! --------------------------
  DO j=1-horiz_limit,rows+horiz_limit
    i=1-horiz_limit/2
    k=MODULO(j-2*i+(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
    orog_1d(2*i) = orog(i,k)
    kk=k
!   1D slices are created by looping along the slowest changing
!   direction, setting up both the main grid-points and the previous
!   bisected point on each loop
    DO i=2-horiz_limit/2,row_length+horiz_limit/2
      k=MODULO(j-2*i+(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
!     The separation angle between main points is calculated and then
!     halved to give the angle to the bisected points
      sep_ang(2*i-2)=MIN(ACOS(sinlat(k)*sinlat(kk)+coslat(k)*coslat(kk) &
                         *cos_delta_lambda(i-1) ) * 0.5, max_ang)
      sep_ang(2*i-1)= sep_ang(2*i-2)
      orog_1d(2*i) = orog(i,k)
      jj=MODULO(k+1-(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
!     The orographic height of the bisected point is interpolated
      orog_1d(2*i-1) = (orog(i,jj)+orog(i-1,jj))*0.5
      kk=k
    END DO
    ! First point in array is beyond the horizon so use arbitrary values
    orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
    sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

    CALL horizon_1d(orog_1d, nd_points, 2*row_length, horiz_limit,      &
                    sep_ang, horiz_ang1, horiz_ang2)

!   The horizon angles calculated for the main grid-points are then
!   mapped back on to the 2D grid
    DO i=1,row_length
      k=MODULO(j-2*i+(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      IF (k >= 1 .AND. k <= rows) THEN
        horiz_ang(i,k,ip_NNW) = horiz_ang1(2*i)
        horiz_ang(i,k,ip_SSE) = horiz_ang2(2*i)
      END IF
    END DO
  END DO

! Horizon angles SSW and NNE
! --------------------------
  DO j=1-horiz_limit,rows+horiz_limit
    i=1-horiz_limit/2
    k=MODULO(j+2*i-3*(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
    orog_1d(2*i) = orog(i,k)
    kk=k
    DO i=2-horiz_limit/2,row_length+horiz_limit/2
      k=MODULO(j+2*i-3*(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      sep_ang(2*i-2)=MIN(ACOS(sinlat(k)*sinlat(kk)+coslat(k)*coslat(kk) &
                         *cos_delta_lambda(i-1) ) * 0.5, max_ang)
      sep_ang(2*i-1)= sep_ang(2*i-2)
      orog_1d(2*i) = orog(i,k)
      jj=MODULO(k-1-(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      orog_1d(2*i-1) = (orog(i,jj)+orog(i-1,jj))*0.5
      kk=k
    END DO
    orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
    sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

    CALL horizon_1d(orog_1d, nd_points, 2*row_length, horiz_limit,      &
                    sep_ang, horiz_ang1, horiz_ang2)

    DO i=1,row_length
      k=MODULO(j+2*i-3*(1-horiz_limit),rows+2*horiz_limit)+1-horiz_limit
      IF (k >= 1 .AND. k <= rows) THEN
        horiz_ang(i,k,ip_SSW) = horiz_ang1(2*i)
        horiz_ang(i,k,ip_NNE) = horiz_ang2(2*i)
      END IF
    END DO
  END DO


! Horizon angles ESE and WNW
! --------------------------
  DO i=1-horiz_limit,row_length+horiz_limit
    j=1-horiz_limit/2
    k = MODULO(i-2*j+(1-horiz_limit),row_length+2*horiz_limit)          &
      + 1-horiz_limit
    orog_1d(2*j) = orog(k,j)
    kk=k
    DO j=2-horiz_limit/2,rows+horiz_limit/2
      k = MODULO(i-2*j+(1-horiz_limit),row_length+2*horiz_limit)        &
        + 1-horiz_limit
      sep_ang(2*j-2)=MIN(ACOS(sinlat(j)*sinlat(j-1)                     &
        +coslat(j)*coslat(j-1)*cos_two_delta_lambda(kk) )*0.5, max_ang)
      sep_ang(2*j-1)= sep_ang(2*j-2)
      orog_1d(2*j) = orog(k,j)
      ii = MODULO(k+1-(1-horiz_limit),row_length+2*horiz_limit)         &
         + 1-horiz_limit
      orog_1d(2*j-1) = (orog(ii,j)+orog(ii,j-1))*0.5
      kk=k
    END DO
    orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
    sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

    CALL horizon_1d(orog_1d, nd_points, 2*rows, horiz_limit,            &
                    sep_ang, horiz_ang1, horiz_ang2)

    DO j=1,rows
      k = MODULO(i-2*j+(1-horiz_limit),row_length+2*horiz_limit)        &
        + 1-horiz_limit
      IF (k >= 1 .AND. k <= row_length) THEN
        horiz_ang(k,j,ip_ESE) = horiz_ang1(2*j)
        horiz_ang(k,j,ip_WNW) = horiz_ang2(2*j)
      END IF
    END DO
  END DO

! Horizon angles WSW and ENE
! --------------------------
  DO i=1-horiz_limit,row_length+horiz_limit
    j=1-horiz_limit/2
    k = MODULO(i+2*j-3*(1-horiz_limit),row_length+2*horiz_limit)        &
      + 1-horiz_limit
    orog_1d(2*j) = orog(k,j)
    kk=k
    DO j=2-horiz_limit/2,rows+horiz_limit/2
      k = MODULO(i+2*j-3*(1-horiz_limit),row_length+2*horiz_limit)      &
        + 1-horiz_limit
      sep_ang(2*j-2)=MIN(ACOS(sinlat(j)*sinlat(j-1)                     &
        +coslat(j)*coslat(j-1)*cos_two_delta_lambda(kk) )*0.5, max_ang)
      sep_ang(2*j-1)= sep_ang(2*j-2)
      orog_1d(2*j) = orog(k,j)
      ii = MODULO(k-1-(1-horiz_limit),row_length+2*horiz_limit)         &
         + 1-horiz_limit
      orog_1d(2*j-1) = (orog(ii,j)+orog(ii,j-1))*0.5
      kk=k
    END DO
    orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
    sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

    CALL horizon_1d(orog_1d, nd_points, 2*rows, horiz_limit,            &
                    sep_ang, horiz_ang1, horiz_ang2)

    DO j=1,rows
      k = MODULO(i+2*j-3*(1-horiz_limit),row_length+2*horiz_limit)      &
        + 1-horiz_limit
      IF (k >= 1 .AND. k <= row_length) THEN
        horiz_ang(k,j,ip_WSW) = horiz_ang1(2*j)
        horiz_ang(k,j,ip_ENE) = horiz_ang2(2*j)
      END IF
    END DO
  END DO


! The horizon angle must include self-shadowing from the slope of the
! grid-box in each direction.

  DO k=1,n_horiz_ang
    DO j=1,rows
      DO i=1,row_length
        tilt_ang(i,j,k) = ATAN2( 1.0, -TAN(slope_angle(i,j))            &
          *COS(slope_aspect(i,j)-horiz_aspect(i,j,k)) )
        horiz_ang(i,j,k) = MIN( horiz_ang(i,j,k), tilt_ang(i,j,k) )
      END DO
    END DO
  END DO

! Now calculate the sky-view factor considering the view from the
! inclined plane.

  ALLOCATE(sky(row_length,rows))
  DO j=1,rows
    DO i=1,row_length
      sky(i,j)=SUM(SIN(horiz_ang(i,j,:)+pi_over_2-tilt_ang(i,j,:))**2)  &
        /(REAL(n_horiz_ang)*COS(slope_angle(i,j)))
    END DO
  END DO

END SUBROUTINE skyview
