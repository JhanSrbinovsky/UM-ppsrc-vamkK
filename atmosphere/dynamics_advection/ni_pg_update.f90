! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
        SUBROUTINE NI_pg_update(                                        &
     &                         q, qcl, qcf,                             &
     &                         qcf2, qrain, qgraup,                     &
     &                         q_star, qcl_star, qcf_star,              &
     &                         qcf2_star, qrain_star, qgraup_star,      &
     &                         mix, mix_cl, mix_cf, mix_cf2,            &
     &                         mix_rain, mix_graup, mix_star,           &
     &                         mix_cl_star, mix_cf_star, mix_cf2_star,  &
     &                         mix_rain_star, mix_graup_star,           &
     &                         L_mcr_cf2, L_mcr_rain, L_mcr_graup,      &
     &                         theta_star, theta, exner,                &
     &                         r_theta_levels, r_rho_levels,            &
     &                         sec_theta_latitude,                      &
     &                         delta_lambda, delta_phi, timestep,       &
     &                         recip_dlambda_p, recip_dphi_p,           &
     &                         wt_lambda_p, wt_phi_p,                   &
     &                         alpha_3, alpha_4,                        &
     &                         row_length, rows, n_rows, model_levels,  &
     &                         wet_levels,                              &
     &                         model_domain,                            &
     &                         first_constant_r_rho_level,              &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         at_extremity,                            &
     &                         CycleNo, L_new_tdisc,                    &
     &                         L_qwaterload,                            &
     &                         R_u, R_v, R_w, L_mix_ratio, L_regular )

! Purpose: Interface routine to PG_update
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
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
                   ! number of model levels where moisture is held
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_mix_ratio                                                     &
     &, L_regular                                                       &
                     ! Variable resolution if false
     &, L_mcr_cf2, L_mcr_rain, L_mcr_graup                              &
     &, L_new_tdisc

      Logical :: L_qwaterload           ! add waterloading terms


      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, timestep

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  recip_dlambda_p(1-halo_i:row_length+halo_i)                     &
     &, recip_dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)  &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

      Real                                                              &
           ! time-weighting parameters, see WP 154.
     &  alpha_3                                                         &
     &, alpha_4

      Real                                                              &
           ! trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real, Intent (InOut) ::                                           &
                              ! primary model variables
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &       wet_levels)                                                &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &      wet_levels)                                                 &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &      wet_levels)                                                 &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &      wet_levels)                                                 &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &        wet_levels)                                               &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &         wet_levels)                                              &
     &, q_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &            wet_levels)                                           &
     &, qcl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &           wet_levels)                                            &
     &, qcf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &           wet_levels)                                            &
     &, qcf2_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &            wet_levels)                                           &
     &, qrain_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &             wet_levels)                                          &
     &, qgraup_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &              wet_levels)                                         &
     &, mix (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_levels)                                                &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &      wet_levels)                                                 &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &      wet_levels)                                                 &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &      wet_levels)                                                 &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &        wet_levels)                                               &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &         wet_levels)                                              &
     &, mix_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &            wet_levels)                                           &
     &, mix_cl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           wet_levels)                                            &
     &, mix_cf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           wet_levels)                                            &
     &, mix_cf2_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &            wet_levels)                                           &
     &, mix_rain_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
     &             wet_levels)                                          &
     &, mix_graup_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &              wet_levels)                                         &
     &, theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

! Arguments with Intent INOUT. ie: Input variables changed on output.

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
           ! See WP154 and WP162 for definitions.
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &        model_levels)                                             &
     &, R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &       model_levels)                                              &
     &, R_w (row_length, rows, model_levels)

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

!    Local arrays

      Real                                                              &
     &  moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         wet_levels)                                              &
     &, moist_star (1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y,wet_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle




      IF (lhook) CALL dr_hook('NI_PG_UPDATE',zhook_in,zhook_handle)
      If ( L_mix_ratio ) then  ! moist terms are mixing ratios

        moist = 1.0 + mix + mix_cl + mix_cf
        moist_star = 1.0 + mix_star + mix_cl_star + mix_cf_star
        If(L_mcr_cf2)then
          moist      = moist + mix_cf2
          moist_star = moist_star + mix_cf2_star
        endif
        If(L_mcr_rain)then
          moist      = moist + mix_rain
          moist_star = moist_star + mix_rain_star
        endif
        If(L_mcr_graup)then
          moist      = moist + mix_graup
          moist_star = moist_star + mix_graup_star
        endif

! DEPENDS ON: pg_update
        Call pg_update(                                                 &
     &                 moist, moist_star,                               &
     &                 mix, mix_star, L_mix_ratio,                      &
     &                 theta_star, theta, exner,                        &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 sec_theta_latitude,                              &
     &                 delta_lambda, delta_phi,                         &
     &                 recip_dlambda_p, recip_dphi_p,                   &
     &                 wt_lambda_p, wt_phi_p,                           &
     &                 timestep, alpha_3, alpha_4,                      &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 model_domain,                                    &
     &                 first_constant_r_rho_level,                      &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 L_regular, at_extremity,                         &
     &                 CycleNo, L_new_tdisc,                            &
     &                 R_u, R_v, R_w)

      else  ! L_mix_ratio is false moist terms are specific quantities

      ! Include condensate/hydrometeors in theta_v calculation  
        If (L_qwaterload) Then  
          moist      = qcl + qcf  
          moist_star = qcl_star + qcf_star  
  
          If (L_mcr_cf2) Then  
            moist      = moist + qcf2  
            moist_star = moist_star + qcf2_star  
          End If  
          If (L_mcr_rain) Then  
            moist      = moist + qrain  
            moist_star = moist_star + qrain_star  
          End If  
          If (L_mcr_graup)  Then  
            moist      = moist + qgraup  
            moist_star = moist_star + qgraup_star  
          End If  
        Else !  L_qwaterload = false
          moist = 0.0  
          moist_star = 0.0  
        End If !  L_qwaterload
        
! DEPENDS ON: pg_update
        Call pg_update(                                                 &
     &                 moist, moist_star,                               &
     &                 q, q_star, L_mix_ratio,                          &
     &                 theta_star, theta, exner,                        &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 sec_theta_latitude,                              &
     &                 delta_lambda, delta_phi,                         &
     &                 recip_dlambda_p, recip_dphi_p,                   &
     &                 wt_lambda_p, wt_phi_p,                           &
     &                 timestep, alpha_3, alpha_4,                      &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 model_domain,                                    &
     &                 first_constant_r_rho_level,                      &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 L_regular, at_extremity,                         &
     &                 CycleNo, L_new_tdisc,                            &
     &                 R_u, R_v, R_w)

      endif   !L_mix_ratio

      IF (lhook) CALL dr_hook('NI_PG_UPDATE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_pg_update
