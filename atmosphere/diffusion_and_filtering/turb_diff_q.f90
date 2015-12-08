! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_diff_q

      Subroutine turb_diff_q(                                           &
                             moist,                                     &
                             r_theta_levels, r_rho_levels,              &
                             FV_sec_theta_latitude,                     &
                             cos_v_latitude, sec_v_latitude,            &
                             off_x, off_y, halo_i, halo_j,              &
                             delta_lambda, delta_phi, timestep,         &
                             rows, n_rows, row_length,                  &
                             model_levels, wet_levels, levels,          &
                             coeff_u, coeff_v, moist_star)

      USE turb_diff_ctl_mod, ONLY: turb_ustart, turb_uend,              &
                                   turb_phi_st, turb_phi_end

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

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_levels                                                      &
                       ! number of wet model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, me                                                              &
                  ! processor id
     &, levels     ! number of levels to process

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, timestep

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           !  trigonometric functions
     &  FV_sec_theta_latitude(1-off_x:row_length+off_x,                 &
     &                        1-off_y:rows+off_y)                       &
     &, cos_v_latitude(1-off_x:row_length+off_x,                        &
     &                     1-off_y:n_rows+off_y)                        &
     &, sec_v_latitude(1-off_x:row_length+off_x,                        &
     &                     1-off_y:n_rows+off_y)                        &
!  diffusion coefficient , only 1 side needs halo
      , coeff_u(0:row_length+1, 0:rows+1, levels)                       &
      , coeff_v(0:row_length+1, 0:n_rows+1, levels)

      Real                                                              &
     &  moist (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels )

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  moist_star (1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y, wet_levels )

! Local Variables.

      Integer                                                           &
     &  i, j, k     ! Loop indices

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi

! Local arrays

      Real                                                              &
     &  delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, levels )                            &
     &, dt_over_r_squared_delz(1-off_x:row_length+off_x,                &
     &                       1-off_y:rows+off_y, levels )               &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows)                    &
     &, lambda_term(row_length, rows)                                   &
     &, phi_term(row_length, rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set values and calculate delta_z
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TURB_DIFF_Q',zhook_in,zhook_handle)
      recip_delta_lambda = 1.0 / delta_lambda
      recip_delta_phi = 1.0 / delta_phi

! calculate D(r)/D(eta) about theta levels
      Do k = 1, levels
        if( k < model_levels )then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = r_rho_levels(i,j,k+1) -                  &
     &                           r_rho_levels(i,j,k )
              dt_over_r_squared_delz(i,j,k) = timestep /                &
     &             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) *    &
     &                                              delta_z(i,j,k) )
            End Do
          End Do
        else ! k = model_levels
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = 1.0
              dt_over_r_squared_delz(i,j,k) = timestep /                &
     &             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )
            End Do
          End Do
        endif !  k < model_levels
      End Do   !  k = 1, levels

! ----------------------------------------------------------------------
! Section 2.0  Horizontal Diffusion
! ----------------------------------------------------------------------
!
! "conservative diffusion" diffuses dzQ rather than Q.  Use 
! dz*(dQ/d_lambda) instead of 1/d_lambda * (dzQ) (and similarly for phi)
! to avoid instabilities ocurring in the presence of orography due to 
! volume differences between neighbouring gridboxes.
!

      Do k = 1, levels

! set lambda_term = 0 so no lambda diffusion on polar filter rows
       lambda_term(:,:) = 0.0

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
        Do j = turb_ustart, turb_uend
          Do i = 1-off_x, row_length
            temp(i,j) = (moist(i+1,j,k) - moist(i,j,k) )                &
     &                    * recip_delta_lambda                          &
     &                    * delta_z(i,j,k) * coeff_u(i,j,k)             &
     &                    * sec_v_latitude(i,j) * sec_v_latitude(i,j)
          End Do

          Do i = 1, row_length
            lambda_term(i,j) = recip_delta_lambda *                     &
     &                          ( temp(i,j) - temp(i-1,j))
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

        Do j = turb_phi_st - 1, turb_phi_end
          Do i = 1, row_length
            temp(i,j) = (moist(i,j+1,k)  -  moist(i,j,k) )              &
     &                    * recip_delta_phi                             &
     &                    * delta_z(i,j,k) * coeff_v(i,j,k)             &
     &                    * cos_v_latitude(i,j)
          End Do
        End Do

        Do j = turb_phi_st, turb_phi_end
          Do i = 1, row_length
            phi_term(i,j) = (temp(i,j) - temp(i,j-1)) * recip_delta_phi*&
     &                        FV_sec_theta_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do j = turb_phi_st, turb_phi_end
          Do i = 1, row_length
              moist_star(i,j,k) = moist_star(i,j,k) +                   &
     &                            dt_over_r_squared_delz(i,j,k) *       &
     &                         ( lambda_term(i,j) + phi_term(i,j) )
          End Do
        End Do

      END DO ! k = 1, levels

! DEPENDS ON: swap_bounds
      call Swap_Bounds(                                                 &
     &                 moist_star, row_length, rows, levels,            &
     &                 off_x, off_y, fld_type_p, .false.)

      IF (lhook) CALL dr_hook('TURB_DIFF_Q',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE turb_diff_q

