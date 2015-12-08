! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_diff_q_varres

SUBROUTINE turb_diff_q_varres(                                    &
                       moist,                                     &
                       r_theta_levels, r_rho_levels,              &
                       sec_theta_latitude,                        &
                       cos_v_latitude, sec_v_latitude,            &
                       off_x, off_y, halo_i, halo_j,              &
                       recip_dlamp, recip_dphip,                  &
                       recip_dlamu, recip_dphiv,                  &
                       timestep,                                  &
                       rows, n_rows, row_length,                  &
                       model_levels, wet_levels, levels,          &
                       coeff_u, coeff_v,                &
                       moist_star)

! Purpose:
!          Turbulent diffusion of a theta field
!          Based on conservative diffusion operator.
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParParams
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                   ! number of point on a row.
, rows                                                            &
                   ! number of rows.
, n_rows                                                          &
                   ! number of v rows.
, model_levels                                                    &
                   ! number of model levels.
, wet_levels                                                      &
                 ! number of wet model levels.
, halo_i                                                          &
                ! Size of halo in i direction.
, halo_j                                                          &
                ! Size of halo in j direction.
, off_x                                                           &
                ! Size of small halo in i
, off_y                                                           &
                ! Size of small halo in j.
, levels     ! number of levels to process

! VarRes horizontal co-ordinate information

REAL, INTENT(IN) ::                                               &
  recip_dlamp(1-halo_i:row_length+halo_i)                         &
, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
, recip_dlamu(1-halo_i:row_length+halo_i)                         &
, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
, timestep

! vertical co-ordinate arrays.
REAL, INTENT(IN) ::                                               &
  r_theta_levels (1-halo_i:row_length+halo_i,                     &
                  1-halo_j:rows+halo_j, 0:model_levels)           &
, r_rho_levels (1-halo_i:row_length+halo_i,                       &
                1-halo_j:rows+halo_j, model_levels)

!  trigonometric functions
REAL, INTENT(IN) ::                                               &
  sec_theta_latitude(1-off_x:row_length+off_x,                    &
                     1-off_y:rows+off_y)                          &
, cos_v_latitude(1-off_x:row_length+off_x,                        &
                     1-off_y:n_rows+off_y)                        &
, sec_v_latitude(1-off_x:row_length+off_x,                        &
                     1-off_y:n_rows+off_y)                        &
!  diffusion coefficient,
, coeff_u(0:row_length+1, 0:  rows+1, levels)                     &
, coeff_v(0:row_length+1, 0:n_rows+1, levels)

REAL, INTENT(IN) ::                                               &
  moist (1-halo_i:row_length+halo_i,                              &
         1-halo_j:rows+halo_j, wet_levels )

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

REAL, INTENT(INOUT) ::                                            &
  moist_star (1-off_x:row_length+off_x,                           &
              1-off_y:rows+off_y, wet_levels )

! Local Variables.

INTEGER ::                                                        &
  i, j, k     ! Loop indices

! Local arrays

REAL ::                                                           &
  delta_z(1-off_x:row_length+off_x,                               &
          1-off_y:rows+off_y)                                     &
, dt_over_r_squared_delz(1-off_x:row_length+off_x,                &
                       1-off_y:rows+off_y)                        &
, temp(1-off_x:row_length+off_x, 1-off_y:rows)                    &
, lambda_term(row_length, rows)                                   &
, phi_term(row_length, rows)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set values and calculate delta_z
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TURB_DIFF_Q_VARRES',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, lambda_term,     &
!$OMP& temp, phi_term, delta_z, dt_over_r_squared_delz)

! calculate D(r)/D(eta) about theta levels
!$OMP DO SCHEDULE(STATIC)
DO k = 1, levels

  IF( k < model_levels )THEN
    DO j = 1-off_y, rows+off_y
      DO i = 1-off_x, row_length+off_x

!              delta_z(i,j,k) = 1.0
        delta_z(i,j) = r_rho_levels(i,j,k+1) -                    &
                           r_rho_levels(i,j,k )
        dt_over_r_squared_delz(i,j) = timestep /                  &
             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) *    &
                                              delta_z(i,j) )
      END DO
    END DO

  ELSE ! k = model_levels

! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
    DO j = 1-off_y, rows+off_y
      DO i = 1-off_x, row_length+off_x
        delta_z(i,j) = 1.0
        dt_over_r_squared_delz(i,j) = timestep /                  &
             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )
      END DO
    END DO
  ENDIF !  k < model_levels

! ----------------------------------------------------------------------
! Section 2.0  Horizontal Diffusion
! ----------------------------------------------------------------------
!
! "conservative diffusion" diffuses dzQ rather than Q.  Use 
! dz*(dQ/d_lambda) instead of 1/d_lambda * (dzQ) (and similarly for phi)
! to avoid instabilities ocurring in the presence of orography due to 
! volume differences between neighbouring gridboxes.
!


! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
  DO j = 1, rows
    DO i = 1-off_x, row_length
      temp(i,j) = (moist(i+1,j,k) - moist(i,j,k) )                &
                    * recip_dlamp(i)                              &
                    * delta_z(i,j) * coeff_u(i,j,k)               &
                    * sec_v_latitude(i,j) * sec_v_latitude(i,j)
    END DO

    DO i = 1, row_length
      lambda_term(i,j) = recip_dlamu(i-1) *                       &
                          ( temp(i,j) - temp(i-1,j))
    END DO
  END DO


! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! Section 2.3  Add diffusion terms to physics increment moist_star.
! ----------------------------------------------------------------------

  DO j = 0, rows
    DO i = 1, row_length
      temp(i,j) = (moist(i,j+1,k)  -  moist(i,j,k) )              &
                    * recip_dphip(i,j)                            &
                    * delta_z(i,j) * coeff_v(i,j,k)               &
                    * cos_v_latitude(i,j)
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
      phi_term(i,j) = (temp(i,j) - temp(i,j-1)) *                 &
                        recip_dphiv(i,j-1) *                      &
                        sec_theta_latitude(i,j)
      moist_star(i,j,k) = moist_star(i,j,k) +                     &
                          dt_over_r_squared_delz(i,j) *           &
                         ( lambda_term(i,j) + phi_term(i,j) )
    END DO
  END DO

END DO ! k = 1, levels
!$OMP END DO

!$OMP END PARALLEL

! DEPENDS ON: swap_bounds
CALL swap_bounds(                                                 &
                 moist_star, row_length, rows, levels,            &
                 off_x, off_y, fld_type_p, .false.)

IF (lhook) CALL dr_hook('TURB_DIFF_Q_VARRES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE turb_diff_q_VarRes

