! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Correction for the angle of solar incidence on sloping terrain.
!
! Description:
!   Calculate the orographic correction to be applied to the
!   Direct SW flux at the surface using the angle of solar
!   incidence on sloping terrain.
!   Uses data from global data module 'solinc_data'.
!
! Method:
!   Spherical trigonometry.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.

SUBROUTINE solinc(row_length,rows,cos_zenith_angle,cosz_beg,cosz_end)

  USE conversions_mod, ONLY: pi
  USE solinc_data, ONLY:                                                      &
    sol_bearing, orog_corr, slope_aspect, slope_angle,                        &
    horiz_ang, horiz_aspect, n_horiz_ang, l_skyview
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Subroutine arguments

  INTEGER, INTENT(IN) :: &
       row_length, rows         ! grid size
  REAL, DIMENSION(row_length,rows), INTENT(IN) :: &
       cos_zenith_angle, &      ! Cosine of the solar zenith angle
       cosz_beg, cosz_end       ! Cos at beginning and end of timestep

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! Local variables

  INTEGER :: i,j
  INTEGER :: max_sub(row_length,rows)
  INTEGER :: min_sub(row_length,rows)

  REAL :: rel_horiz_aspect(row_length,rows,n_horiz_ang)
  REAL :: cosz_hor(row_length,rows)
  REAL :: cosz_lit(row_length,rows)
  REAL :: lit_frac(row_length,rows)

!- End of header

  IF (lhook) CALL dr_hook('SOLINC',zhook_in,zhook_handle)

! orog_corr is equal to the cosine of the angle between the incoming
! solar insolation and the normal to the mean slope, divided by
! (cosine of the solar zenith angle x cosine of the slope angle).

  WHERE (cos_zenith_angle /= 0.0 .AND. slope_angle /= 0.0)

       orog_corr = 1.0 + TAN( ACOS(cos_zenith_angle) ) * &
          TAN(slope_angle) * COS( sol_bearing - slope_aspect )

       orog_corr = MAX(orog_corr,EPSILON(orog_corr))

  ELSEWHERE
       orog_corr = 1.0
  END WHERE

  IF (l_skyview) THEN
!   Calculate cosine of terrain horizon angle
    rel_horiz_aspect = MODULO( horiz_aspect - &
      SPREAD(sol_bearing,3,n_horiz_ang), pi*2.)
    min_sub=MINLOC(rel_horiz_aspect, DIM=3)
    max_sub=MAXLOC(rel_horiz_aspect, DIM=3)
    DO j=1,rows
      DO i=1,row_length
        cosz_hor(i,j) = COS( (horiz_ang(i,j,min_sub(i,j))                     &
          *(pi*2. - rel_horiz_aspect(i,j,max_sub(i,j)))                       &
          +horiz_ang(i,j,max_sub(i,j))*rel_horiz_aspect(i,j,min_sub(i,j)))    &
          /(pi*2. - rel_horiz_aspect(i,j,max_sub(i,j))                        &
          +rel_horiz_aspect(i,j,min_sub(i,j))) )
      END DO
    END DO

!   Adjust orographic correction for times when sun is below horizon
!   (assumes change of zenith angle is roughly linear over the timestep)
    WHERE (cos_zenith_angle /= 0.0)
      WHERE (cosz_hor > cosz_beg .AND. cosz_hor > cosz_end)
!       Sun is below the horizon
        orog_corr = EPSILON(orog_corr)
      ELSEWHERE (cosz_hor > cosz_beg .AND. cosz_hor < cosz_end)
!       Sun rises during timestep
        lit_frac=(cosz_end-cosz_hor)/(cosz_end-cosz_beg)
        cosz_lit=(cosz_end+cosz_hor)*0.5
        orog_corr = lit_frac*(cosz_lit+SIN(ACOS(cosz_lit))*TAN(slope_angle)*  &
          COS(sol_bearing-slope_aspect))/cos_zenith_angle    
        orog_corr = MAX(orog_corr,EPSILON(orog_corr))
      ELSEWHERE (cosz_hor < cosz_beg .AND. cosz_hor > cosz_end)
!       Sun sets during timestep
        lit_frac=(cosz_beg-cosz_hor)/(cosz_beg-cosz_end)
        cosz_lit=(cosz_beg+cosz_hor)*0.5
        orog_corr = lit_frac*(cosz_lit+SIN(ACOS(cosz_lit))*TAN(slope_angle)*  &
          COS(sol_bearing-slope_aspect))/cos_zenith_angle
        orog_corr = MAX(orog_corr,EPSILON(orog_corr))
      END WHERE
    END WHERE
  END IF

  IF (lhook) CALL dr_hook('SOLINC',zhook_out,zhook_handle)

END SUBROUTINE solinc
