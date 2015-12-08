! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_diff_v_VarRes

SUBROUTINE turb_diff_v_varres(                                    &
                       v, r_at_v,                                 &
                       r_theta_levels, r_rho_levels,              &
                       sec_v_latitude,                            &
                       cos_theta_latitude, sec_theta_latitude,    &
                       off_x, off_y, halo_i, halo_j,              &
                       recip_dlamp, recip_dphip,                  &
                       recip_dlamu, recip_dphiv,                  &
                       timestep,                                  &
                       rows, n_rows, row_length,                  &
                       model_levels, levels,                      &
                       coeff_theta, coeff_centre, R_v )

! Purpose:
!          Turbulent diffusion of a v-type field
!          Based on conservative diffusion operator.
!          3d diffusion coefficient
!          Active wherever diff_coeff > 0
!          If you need to switch off over steep slopes
!          then set diff_coeff = 0 at required points before entry.
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
                   ! number of v-rows
, model_levels                                                    &
                   ! number of model levels.
, halo_i                                                          &
                   ! Size of halo in i direction.
, halo_j                                                          &
                   ! Size of halo in j direction.
, off_x                                                           &
                   ! Size of small halo in i
, off_y                                                           &
                   ! Size of small halo in j.
, levels           ! levels to process

! VarRes horizontal co-ordinate information
REAL, INTENT(IN) ::                                               &
  recip_dlamp(1-halo_i:row_length+halo_i)                         &
, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
, recip_dlamu(1-halo_i:row_length+halo_i)                         &
, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
, timestep

! vertical co-ordinate arrays and input arrays
REAL, INTENT(IN) ::                                               &
  r_at_v (1-halo_i:row_length+halo_i,                             &
          1-halo_j:n_rows+halo_j, model_levels)                   &
, r_rho_levels (1-halo_i:row_length+halo_i,                       &
                1-halo_j:rows+halo_j, model_levels)               &
, r_theta_levels (1-halo_i:row_length+halo_i,                     &
                  1-halo_j:rows+halo_j, 0:model_levels)           &
, v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     model_levels)

!  trigonometric functions
REAL, INTENT(IN) ::                                               &
  sec_v_latitude (1-off_x:row_length+off_x,                       &
                  1-off_y:n_rows+off_y)                           &
, cos_theta_latitude (1-off_x:row_length+off_x,                   &
                      1-off_y:rows+off_y)                         &
, sec_theta_latitude (1-off_x:row_length+off_x,                   &
                      1-off_y:rows+off_y)                         &
!  diffusion coefficients common to both u and v
!  needed at theta points and at centre of grid cell
, coeff_theta (0:row_length+1, 0:rows+1, levels)                  &
, coeff_centre(0:row_length+1, 0:rows+1, levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

REAL, INTENT(INOUT) ::                                            &
  R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
       model_levels)

! Local Variables.

INTEGER ::                                                        &
  i, j, k           ! Loop indices

! Local arrays

REAL ::                                                           &
  r_theta_at_v(1-off_x:row_length+off_x,                          &
               1-off_y:n_rows+off_y, 0:levels)                    &
, delta_z(1-off_x:row_length+off_x,                               &
          1-off_y:n_rows+off_y)                                   &
, recip_r_squared_delz(1-off_x:row_length+off_x,                  &
                       1-off_y:n_rows+off_y)                      &
, temp(1-off_x:row_length, n_rows+ 1)                             &
, diff_lambda(1-off_x:row_length, n_rows)                         &
, lambda_term(row_length, n_rows)                                 &
, phi_term(row_length, n_rows)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set r field  and calculate delta_z
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TURB_DIFF_V_VARRES',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, phi_term, temp,  &
!$OMP& lambda_term, delta_z, recip_r_squared_delz)

!$OMP DO SCHEDULE(STATIC)
DO k = 0, levels
  DO j = 1-off_y, n_rows+off_y
    DO i = 1-off_x, row_length+off_x
      r_theta_at_v(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +             &
                                  r_theta_levels(i  ,j,k) )
    END DO
  END DO
END DO ! k = 0, levels
!$OMP END DO 

!$OMP MASTER

! call swap_bounds to set halo points

! DEPENDS ON: swap_bounds
CALL swap_bounds                                                  &
                (r_theta_at_v,                                    &
                 row_length, n_rows, levels+ 1,                   &
                 off_x, off_y, fld_type_v, .false.)

!$OMP END MASTER
!$OMP BARRIER 

! calculate D(r)/D(eta) about theta levels
!$OMP DO SCHEDULE(STATIC)
DO k = 1, levels

  IF ( k < model_levels ) THEN
    DO j = 1-off_y, n_rows+off_y
      DO i = 1-off_x, row_length+off_x
        delta_z(i,j) = r_theta_at_v(i,j,k) -                      &
                           r_theta_at_v(i,j,k-1)
        recip_r_squared_delz(i,j) = 1.0 / ( r_at_v(i,j,k) *       &
                                               r_at_v(i,j,k) *    &
                                               delta_z(i,j) )
      END DO
    END DO
  ELSE ! ( k = model_levels)
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
    DO j = 1-off_y, n_rows+off_y
      DO i = 1-off_x, row_length+off_x
        delta_z(i,j) = 1.0
        recip_r_squared_delz(i,j) = 1.0 / ( r_at_v(i,j,k) *      &
                                               r_at_v(i,j,k) )
      END DO
    END DO
  ENDIF ! k < model_levels


! ----------------------------------------------------------------------
! Section 2.0 Diffusion
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

  DO j = 1, n_rows
    DO i = 1-off_x, row_length
      temp(i,j) = (v(i+1,j,k) - v(i,j,k) )                        &
               * recip_dlamp(i)                                   &
               * delta_z(i,j) * coeff_centre(i,j,k)               &
               * sec_theta_latitude(i,j) * sec_theta_latitude(i,j)
    END DO
    DO i = 1, row_length
      lambda_term(i,j) = recip_dlamu(i-1) *                       &
                           ( temp(i,j) - temp(i-1,j) )
    END DO
  END DO

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! Section 2.3  Add diffusion terms to physics increment R_v.
! ----------------------------------------------------------------------

  DO j = 1, n_rows + 1
    DO i = 1, row_length
      temp(i,j) = (v(i,j,k) - v(i,j-1,k) )                        &
                    * recip_dphiv(i,j-1)                          &
                    * delta_z(i,j-1) * coeff_theta(i,j,k)         &
                    * cos_theta_latitude(i,j)
    END DO
  END DO
  
  DO j = 1, n_rows
    DO i = 1, row_length
      phi_term(i,j) = ( temp(i,j+1) - temp(i,j) ) *               &
                       recip_dphip(i,j) *                         &
                       sec_v_latitude(i,j)
      R_v(i,j,k) = R_v(i,j,k) + timestep *                        &
                       recip_r_squared_delz(i,j) *                &
                     ( lambda_term(i,j) + phi_term(i,j) )
    END DO
  END DO

END DO  ! k = 1, levels
!$OMP END DO

!$OMP END PARALLEL 

! DEPENDS ON: swap_bounds
CALL swap_bounds(                                                 &
                  R_v, row_length, n_rows, levels,                &
                  off_x, off_y, fld_type_v, .true.)

IF (lhook) CALL dr_hook('TURB_DIFF_V_VARRES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE turb_diff_v_varres

