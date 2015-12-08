! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine vert_h_onto_p
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
MODULE vert_h_onto_p_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE vert_h_onto_p(                                         &
                           data_in, row_length, data_rows,        &
                           data_levels, desired_p,                &
                           p_theta_levels,                        &
                           theta, exner_theta_levels,             &
                           exner_rho_levels,                      &
                           boundary_layer_levels,                 &
                           off_x, off_y, halo_i, halo_j,          &
                           p_at_data, interp_order,               &
                           data_out )

! Purpose:
!          Performs linear interpolation in exner of the height
!          field onto a pressure surface. Where surface is below
!          bottom of model, data isd extrapolated as in the UM.

! Method:

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

USE atm_fields_bounds_mod, ONLY: pdims_s, tdims_s, wdims_l, pdims_l

USE level_heights_mod, ONLY:                                      &
                r_theta_levels, r_rho_levels

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: r, cp, kappa, lapse, p_zero

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) ::   row_length    ! number of points on a row
INTEGER, INTENT(IN) ::   data_rows     ! number of rows of data
INTEGER, INTENT(IN) ::   data_levels   ! number of levels of data
INTEGER, INTENT(IN) ::   boundary_layer_levels
INTEGER, INTENT(IN) ::   interp_order ! 1 = linear, 3 = cubic, 5=quintic
INTEGER, INTENT(IN) ::   off_x
INTEGER, INTENT(IN) ::   off_y         ! halo sizes
INTEGER, INTENT(IN) ::   halo_i
INTEGER, INTENT(IN) ::   halo_j        ! large halo sizes

REAL, INTENT(IN) ::   desired_p

REAL, INTENT(IN) ::  data_in (row_length, data_rows, data_levels)
REAL, INTENT(IN) ::  p_at_data (pdims_s%i_start:pdims_s%i_end,    &
                                pdims_s%j_start:pdims_s%j_end,    &
                                             data_levels) 
REAL, INTENT(IN) ::  p_theta_levels(tdims_s%i_start:tdims_s%i_end,&
                                    tdims_s%j_start:tdims_s%j_end,&
                                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) ::  theta         (tdims_s%i_start:tdims_s%i_end,&
                                    tdims_s%j_start:tdims_s%j_end,&
                                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) :: exner_theta_levels                            &
                                   (tdims_s%i_start:tdims_s%i_end,&
                                    tdims_s%j_start:tdims_s%j_end,&
                                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) :: exner_rho_levels                              &
                            (pdims_s%i_start:pdims_s%i_end,       &
                             pdims_s%j_start:pdims_s%j_end,       &
                             pdims_s%k_start:pdims_s%k_end + 1) 
                                             
! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) :: data_out (row_length, data_rows)

! Local variables

INTEGER :: i,j,k 
INTEGER :: level_extrap (row_length, data_rows)
INTEGER :: level_below(row_length, data_rows)

REAL    :: power 
REAL    :: t_ref_level_1 
REAL    :: extrap_height
REAL    :: desired_exner
REAL    ::   mdi              ! missing data indicator.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!            Use first level above extrap_height if extrapolation needed
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook('VERT_H_ONTO_P',zhook_in,zhook_handle)
extrap_height = 2000 !m
DO k = 1, data_levels
  DO j = 1, data_rows
    DO i = 1, row_length
      IF ( (r_theta_levels(i,j,k)-r_theta_levels(i,j,0) )  <=     &
            extrap_height ) THEN
        level_extrap(i,j) = k
      END IF
    END DO
  END DO
END DO

DO j = 1, data_rows
  DO i = 1, row_length
    level_below(i,j) = 0
  END DO
END DO

DO k = 1, data_levels - 1
  DO j = 1, data_rows
    DO i = 1, row_length
      IF ( p_at_data(i,j,k)  >=  desired_p ) THEN
        level_below(i,j) = k
      END IF
    END DO
  END DO
END DO

! if requested level is above top of model, set to data_levels+1,
! which will be converted to missing data indicator.

DO j = 1, data_rows
  DO i = 1, row_length
    IF ( desired_p  <   p_at_data(i,j,data_levels) ) THEN
      level_below(i,j) = data_levels+1
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

! change desired pressure into exner equivalent.

desired_exner = (desired_p/p_zero) ** kappa

! required to calculate values below bottom of model
power = (r * lapse) / g

DO j = 1, data_rows
  DO i = 1, row_length

    IF (level_below(i,j)  ==  0) THEN

      t_ref_level_1= lapse *                                      &
           ( r_theta_levels(i,j,level_extrap(i,j)) -              &
             r_rho_levels(i,j,1)  ) /                             &
           (1.- ( p_theta_levels(i,j,level_extrap(i,j)) /         &
                  p_at_data(i,j,1) )**power )

      data_out(i,j) = data_in(i,j,1) +                            &
                     ( t_ref_level_1 / lapse)      *              &
                     (1. - (desired_p/p_at_data(i,j,1))**power)

    ELSE IF (level_below(i,j)  ==  data_levels+1) THEN
      data_out(i,j) = data_in(i,j,data_levels) -                  &
            cp*(desired_exner -exner_rho_levels(i,j,data_levels)) &
            *theta(i,j,data_levels)/g

    ELSE
! linearly interpolate
      data_out (i,j) = ( ( desired_exner -                        &
                       exner_rho_levels(i,j,level_below(i,j)) )   &
                          * data_in (i,j,level_below(i,j)+1)      &
                     -( desired_exner -                           &
                    exner_rho_levels(i,j,level_below(i,j)+1) ) *  &
                           data_in (i,j,level_below(i,j)) ) /     &
                   ( exner_rho_levels(i,j,level_below(i,j)+1) -   &
                       exner_rho_levels(i,j,level_below(i,j)) )

    END IF
  END DO
END DO

! end of routine

IF (lhook) CALL dr_hook('VERT_H_ONTO_P',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vert_h_onto_p
END MODULE vert_h_onto_p_mod
