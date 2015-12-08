! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine T_vert_interp_to_p
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
MODULE t_vert_interp_to_p_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE t_vert_interp_to_p(                                    &
                           t, theta, row_length, rows,            &
                           model_levels, desired_p,               &
                           off_x, off_y,                          &
                            halo_i, halo_j,                       &
                           p_theta_levels,                        &
                           boundary_layer_levels,                 &
                           exner_theta_levels,                    &
                           t_out )

! Purpose:
!          Performs vertical interpolation of temperature to a
!          desired p surface assuming that T varies linearly with
!          geopotential height, which itself is assumed to vary
!          linearly with Exner pressure. Hence T varies linearly with
!          exner pressure. Where the desired surface is above/below
!          the top/bottom data point extrapolation is done
!          respectively: isothermal or
!          from the first level above a fixed height (200m)

! Method:
!          Simple linear interpolation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2programming standards.

USE atm_fields_bounds_mod, ONLY: tdims_s, wdims_l

USE earth_constants_mod, ONLY: g
USE atmos_constants_mod, ONLY: kappa, p_zero, r, lapse
USE level_heights_mod, ONLY:  r_theta_levels


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) ::   row_length    ! number of points on a row
INTEGER, INTENT(IN) ::   rows          ! number of rows of data
INTEGER, INTENT(IN) ::   model_levels  ! number of levels of data
INTEGER, INTENT(IN) ::   boundary_layer_levels
INTEGER, INTENT(IN) ::   off_x
INTEGER, INTENT(IN) ::   off_y         ! halo sizes
INTEGER, INTENT(IN) ::   halo_i
INTEGER, INTENT(IN) ::   halo_j        ! large halo sizes


REAL, INTENT(IN) ::   desired_p

REAL, INTENT(IN) ::  t (row_length, rows, model_levels)
REAL, INTENT(IN) ::  theta (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) ::  exner_theta_levels                           &
                           (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) ::  p_theta_levels                               &
                           (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)


! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) ::  t_out (row_length, rows)

! Local variables

INTEGER ::  i,j,k  
INTEGER ::  level_extrap (row_length, rows)
INTEGER ::  level_below(row_length, rows)

REAL    ::  desired_exner 
REAL    ::  power  
REAL    ::  extrap_height  !height to determine level_extrap

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 0. Initialize arrays
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook('T_VERT_INTERP_TO_P',zhook_in,zhook_handle)
DO j = 1, rows
  DO i = 1, row_length
    level_below(i,j) = 0
    level_extrap(i,j) = 1
  END DO
END DO

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!         Use nearest level below extrap_height if extrapolation needed
!         200m is chosed as realistic physical value for atmosphere
!         If first theta level is above this first level is used
! ----------------------------------------------------------------------
extrap_height = 200 !m
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      IF ( (r_theta_levels(i,j,k)-r_theta_levels(i,j,0) )  <=     &
            extrap_height ) THEN
        level_extrap(i,j) = k
      END IF
    END DO
  END DO
END DO


! change desired pressure into exner equivalent.

desired_exner = (desired_p/p_zero) ** kappa


DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      IF ( exner_theta_levels(i,j,k)  >=  desired_exner ) THEN
        level_below(i,j) = k
      END IF
    END DO
  END DO
END DO

! if requested level is above top of model, set level_below to -1,
! which will be converted to use isothermal extrapolation

DO j = 1, rows
  DO i = 1, row_length
    IF ( desired_exner  <                                         &
         exner_theta_levels(i,j,model_levels) ) THEN
      level_below(i,j) = -1
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

power = lapse * r / g
DO j = 1, rows
  DO i = 1, row_length

    IF (level_below(i,j)  ==  -1) THEN
! isothermal extrapolation above top level
      t_out(i,j) = t(i,j,model_levels)

    ELSE IF (level_below(i,j)  ==  0) THEN
! extrapolate
      t_out(i,j) = t(i,j,level_extrap(i,j)) *                     &
                   (desired_p/                                    &
                 p_theta_levels(i,j,level_extrap(i,j)))           &
                  ** power

    ELSE
! linear interpolation
      t_out(i,j) =                                                &
                  ( ( desired_exner -                             &
                     exner_theta_levels(i,j,level_below(i,j)) )   &
                          * t(i,j,level_below(i,j)+1)             &
                     -( desired_exner -                           &
                 exner_theta_levels(i,j,level_below(i,j)+1) ) *   &
                           t(i,j,level_below(i,j)) ) /            &
                 ( exner_theta_levels(i,j,level_below(i,j)+1) -   &
                     exner_theta_levels(i,j,level_below(i,j)) )

    END IF

  END DO
END DO

! end of routine

IF (lhook) CALL dr_hook('T_VERT_INTERP_TO_P',zhook_out,zhook_handle)
RETURN
END SUBROUTINE t_vert_interp_to_p

END MODULE t_vert_interp_to_p_mod
