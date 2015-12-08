! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
        SUBROUTINE NI_Update_Rho(                                       &
     &                      rho, rho_n, rho_np1, inc_rho,               &
     &                      u, v, w, R_u, R_v, R_w,                     &
     &                      q, qcl, qcf, qcf2, qrain, qgraup,           &
     &                      q_star, qcl_star, qcf_star, qcf2_star,      &
     &                      qrain_star, qgraup_star,                    &
     &                      q_np1, qcl_np1, qcf_np1, qcf2_np1,          &
     &                      qrain_np1, qgraup_np1,                      &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_v_star, mix_cl_star, mix_cf_star,       &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      L_mcr_cf2, L_mcr_rain, L_mcr_graup,         &
     &                      timestep, CycleNo, NumCycles,               &
     &                      rows, n_rows, row_length,                   &
     &                      model_levels, wet_model_levels,             &
     &                      model_domain, first_constant_r_rho_level,   &
     &                      alpha_1, alpha_2, rims_to_do,               &
     &                      nproc, gc_proc_row_group,                   &
     &                      L_regular, at_extremity, global_row_length, &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      cos_v_latitude, delta_lambda, delta_phi,    &
     &                      dlambda_p, dphi_p,                          &
     &                      recip_dlambda_u, recip_dphi_v,              &
     &                      wt_lambda_p, wt_lambda_u,                   &
     &                      wt_phi_p, wt_phi_v,                         &
     &                      r_theta_levels, r_rho_levels,               &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      FV_sec_theta_latitude,                      &
     &                      wet_to_dry_n, wet_to_dry_np1,               &
     &                      L_do_increment, L_new_tdisc,                &
     &                      L_mix_ratio, L_dry )

! Purpose: Interface routine to update_rho
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held.
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, offx                                                            &
     &, offy                                                            &
     &, rims_to_do                                                      &
                        ! rim size of lbc weights = 1
     &, nproc                                                           &
                    ! Total number of processors
     &, gc_proc_row_group                                               &
                          ! Group id for processors on the same row
     &, global_row_length                                               &
     &, NumCycles                                                       &
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                    ! false if variable resolution
     &, L_new_tdisc                                                     &
     &, L_mcr_cf2, L_mcr_rain, L_mcr_graup

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  timestep                                                        &
     &, alpha_1                                                         &
     &, alpha_2

      Real                                                              &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, rho_n (1-offx:row_length+offx,                                  &
     &         1-offy:rows+offy, model_levels)

      Real, Intent(In) ::                                               &
     &  q_np1 (1-offx:row_length+offx, 1-offy:rows+offy,                &
     &          wet_model_levels)                                       &
     &, qcl_np1 (1-offx:row_length+offx, 1-offy:rows+offy,              &
     &          wet_model_levels)                                       &
     &, qcf_np1 (1-offx:row_length+offx, 1-offy:rows+offy,              &
     &          wet_model_levels)                                       &
     &, qcf2_np1 (1-offx:row_length+offx, 1-offy:rows+offy,             &
     &          wet_model_levels)                                       &
     &, qrain_np1 (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &          wet_model_levels)                                       &
     &, qgraup_np1 (1-offx:row_length+offx, 1-offy:rows+offy,           &
     &          wet_model_levels)

      Real, Intent(Out) ::                                              &
     &  inc_rho(1-offx:row_length+offx,                                 &
     &             1-offy:rows+offy,model_levels)                       &
     &, wet_to_dry_n (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &          model_levels)

      Real                                                              &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Real                                                              &
           ! trigonometric arrays.
     &  FV_sec_theta_latitude (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy)                        &
     &, cos_v_latitude (1-offx:row_length+offx,                         &
     &                  1-offy:n_rows+offy)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)                                   &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, recip_dlambda_u(1-halo_i : row_length + halo_i)                 &
     &, recip_dphi_v( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : n_rows+halo_j )                        &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Logical :: L_do_increment
      Logical :: L_mix_ratio
      Logical :: L_dry

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  rho (1-offx:row_length+offx, 1-offy:rows+offy,                  &
     &       model_levels)

      Real, Intent (InOut) ::                                           &
                              ! primary model variables
     &  rho_np1 (1-offx:row_length+offx, 1-offy:rows+offy,              &
     &       model_levels)                                              &
     &, q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     wet_model_levels)                                            &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_model_levels)                                            &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_model_levels)                                            &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &     wet_model_levels)                                            &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_model_levels)                                            &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &     wet_model_levels)                                            &
     &, q_star (1-offx:row_length+offx, 1-offy:rows+offy,               &
     &          wet_model_levels)                                       &
     &, qcl_star (1-offx:row_length+offx, 1-offy:rows+offy,             &
     &          wet_model_levels)                                       &
     &, qcf_star (1-offx:row_length+offx, 1-offy:rows+offy,             &
     &          wet_model_levels)                                       &
     &, qcf2_star (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &          wet_model_levels)                                       &
     &, qrain_star (1-offx:row_length+offx, 1-offy:rows+offy,           &
     &          wet_model_levels)                                       &
     &, qgraup_star (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_model_levels)                                       &
     &, mix_v (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cf(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cf2(1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, wet_model_levels)                 &
     &, mix_rain(1-halo_i:row_length+halo_i,                            &
     &           1-halo_j:rows+halo_j, wet_model_levels)                &
     &, mix_graup(1-halo_i:row_length+halo_i,                           &
     &            1-halo_j:rows+halo_j, wet_model_levels)               &
     &, mix_v_star(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy,wet_model_levels)                   &
     &, mix_cl_star(1-offx:row_length+offx,                             &
     &              1-offy:rows+offy,wet_model_levels)                  &
     &, mix_cf_star(1-offx:row_length+offx,                             &
     &              1-offy:rows+offy,wet_model_levels)                  &
     &, mix_cf2_star(1-offx:row_length+offx,                            &
     &               1-offy:rows+offy,wet_model_levels)                 &
     &, mix_rain_star(1-offx:row_length+offx,                           &
     &                1-offy:rows+offy,wet_model_levels)                &
     &, mix_graup_star(1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy,wet_model_levels)               &
     &, mix_v_np1 (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &               wet_model_levels)                                  &
     &, mix_cl_np1 (1-offx:row_length+offx, 1-offy:rows+offy,           &
     &               wet_model_levels)                                  &
     &, mix_cf_np1 (1-offx:row_length+offx, 1-offy:rows+offy,           &
     &               wet_model_levels)                                  &
     &, mix_cf2_np1 (1-offx:row_length+offx,                            &
     &               1-offy:rows+offy,wet_model_levels)                 &
     &, mix_rain_np1 (1-offx:row_length+offx,                           &
     &               1-offy:rows+offy,wet_model_levels)                 &
     &, mix_graup_np1 (1-offx:row_length+offx,                          &
     &               1-offy:rows+offy,wet_model_levels)

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

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('NI_UPDATE_RHO',zhook_in,zhook_handle)
        If (L_mix_ratio)then     
  !   put mixing ratios quantities in arg list

! DEPENDS ON: update_rho
          Call Update_Rho(                                              &
     &                 timestep, NumCycles, CycleNo,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_model_levels, model_domain,    &
     &                 first_constant_r_rho_level,                      &
     &                 alpha_1, alpha_2, rims_to_do,                    &
     &                 nproc, gc_proc_row_group,                        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 dlambda_p, dphi_p,                               &
     &                 recip_dlambda_u, recip_dphi_v,                   &
     &                 wt_lambda_p, wt_lambda_u,                        &
     &                 wt_phi_p, wt_phi_v,                              &
     &                 R_u, R_v, R_w,                                   &
     &                 u, v, w, rho, rho_np1,                           &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 mix_v_np1, mix_cl_np1, mix_cf_np1,               &
     &                 mix_cf2_np1, mix_rain_np1, mix_graup_np1,        &
     &                 L_mcr_cf2, L_mcr_rain, L_mcr_graup,              &
     &                 r_theta_levels, r_rho_levels, eta_theta_levels,  &
     &                 eta_rho_levels, FV_sec_theta_latitude,           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 rho_n, inc_rho, L_do_increment,                  &
     &                 L_new_tdisc, L_mix_ratio, L_dry)

      else  ! L_mix_ratio F put specific quantities in arg list

! DEPENDS ON: update_rho
        Call Update_Rho(                                                &
     &                 timestep, NumCycles, CycleNo,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_model_levels, model_domain,    &
     &                 first_constant_r_rho_level,                      &
     &                 alpha_1, alpha_2, rims_to_do,                    &
     &                 nproc, gc_proc_row_group,                        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 dlambda_p, dphi_p,                               &
     &                 recip_dlambda_u, recip_dphi_v,                   &
     &                 wt_lambda_p, wt_lambda_u,                        &
     &                 wt_phi_p, wt_phi_v,                              &
     &                 R_u, R_v, R_w,                                   &
     &                 u, v, w, rho, rho_np1,                           &
     &                 q, qcl, qcf, qcf2, qrain, qgraup,                &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 q_np1, qcl_np1, qcf_np1,                         &
     &                 qcf2_np1, qrain_np1, qgraup_np1,                 &
     &                 L_mcr_cf2, L_mcr_rain, L_mcr_graup,              &
     &                 r_theta_levels, r_rho_levels, eta_theta_levels,  &
     &                 eta_rho_levels, FV_sec_theta_latitude,           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 rho_n, inc_rho, L_do_increment,                  &
     &                 L_new_tdisc, L_mix_ratio, L_dry)

      endif  !  L_mix_ratio

      IF (lhook) CALL dr_hook('NI_UPDATE_RHO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_Update_Rho
