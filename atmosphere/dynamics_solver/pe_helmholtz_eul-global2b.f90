! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine PE_Helmholtz_global_2B
      SUBROUTINE PE_Helmholtz_global_2B(                                &
     &                       u, v, w, r_theta_levels,                   &
     &                       r_rho_levels, p, rho, rho_np1,             &
     &                       theta, theta_star, theta_np1,              &
     &                       q, q_star, q_np1,                          &
     &                       qcl, qcf, qcf2, qrain, qgraup,             &
     &                       qcl_star, qcf_star,                        &
     &                       qcf2_star, qrain_star, qgraup_star,        &
     &                       qcl_np1 , qcf_np1 ,                        &
     &                       qcf2_np1 , qrain_np1 , qgraup_np1 ,        &
     &                       L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,     &
     &                       L_qwaterload,                              &
     &                       rho_Km, cH, G_term_tol,                    &
     &                       exner_rho_levels, frictional_timescale,    &
     &                       sec_theta_latitude, cos_v_latitude,        &
     &                       FV_cos_theta_latitude,                     &
     &                       FV_sec_theta_latitude,                     &
     &                       f3_at_u, f3_at_v,                          &
     &                       timestep, timestep_number,                 &
     &                       row_length, rows, n_rows,                  &
     &                       model_levels, wet_model_levels,            &
     &                       boundary_layer_levels,                     &
     &                       delta_lambda, delta_phi,                   &
     &                       lambda_p, phi_p, lambda_u, phi_v,          &
     &                       dlambda_p, dphi_p, dlambda_u, dphi_v,      &
     &                       recip_dlamp, recip_dphip,                  &
     &                       recip_dlamu, recip_dphiv,                  &
     &                       wt_lambda_p, wt_phi_p,                     &
     &                       wt_lambda_u, wt_phi_v,                     &
     &                       GCR_max_iterations, GCR_diagnostics,       &
     &                       GCR_its_switch, GCR_its_avg_step,          &
     &                       GCR_max_its, GCR_min_its, GCR_sum_its,     &
     &                       GCR_max_time, GCR_min_time,                &
     &                       GCR_tol_res, GCR_tol_abs,                  &
     &                       GCR_use_tol_abs, GCR_zero_init_guess,      &
     &                       GCR_use_residual_Tol,                      &
     &                       GCR_adi_add_full_soln, L_gcr_fast_x,       &
     &                       GCR_precon_option, GCR_ADI_Pseudo_timestep,&
     &                       GCR_n_ADI_pseudo_timesteps,                &
     &                       eta_theta_levels, eta_rho_levels,          &
     &                       alpha_1, alpha_2, alpha_3, alpha_4,        &
     &                       alpha_Cd,                                  &
     &                       model_domain, L_physics,                   &
     &                       GCR_Restart_value,                         &
     &                       first_constant_r_rho_level,                &
     &                       first_constant_r_rho_level_m1,             &
     &                       R_u, R_v, R_w, exner_prime,                &
     &                       dtheta_dr_term,                            &
     &                       me, n_proc, n_procx, n_procy,              &
     &                       halo_i, halo_j, l_datastart,               &
     &                       L_regular, at_extremity, off_x, off_y,     &
     &                       proc_row_group, proc_col_group,            &
     &                       global_row_length, global_rows,            &
     &                       g_rows, g_row_length,                      &
     &                       CycleNo, L_new_tdisc, L_fint_theta,        &
     &                       ldump                                      &
     &                       )
      USE diag_print_mod

! Purpose:
!          Sets up and solves Primitive equation Helmholtz equation.
!          Calculates increments to u, v, w and exner.
!
! Method:
!          Is described in ;
!
!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE atmos_constants_mod, ONLY:  &
          cp, r, kappa, recip_epsilon

       USE global_2d_sums_mod, ONLY: &
           global_2d_sums

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

! Parameters required for dimensioning some of the arguments

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                         ! number of model levels where moisture is held
     &, boundary_layer_levels                                           &
                              ! number of boundary layer levels.
     &, timestep_number

      Integer                                                           &
     &  me                                                              &
                     ! My processor number
     &, n_proc                                                          &
                     ! Total number of processors
     &, n_procx                                                         &
                     ! Number of processors in longitude
     &, n_procy                                                         &
                     ! Number of processors in latitude
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_rows                                                          &
                   ! Local number of rows in a v field
     &, off_x, off_y                                                    &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, CycleNo

      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                    !  true for regular resolution
     &, L_new_tdisc                                                     &
     &, L_fint_theta    ! true:  fully-interpolating semi-lagrangian
                        !        theta advection will be used
                     ! false: standard non-interpolating in the vertical

      Logical  :: L_mcr_qcf2             ! is qcf2 present
      Logical  :: L_mcr_qrain            ! is qrain present
      Logical  :: L_mcr_qgraup           ! is qgraup present
      Logical  :: L_qwaterload           ! add waterloading terms

      Integer                                                           &
     &  first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1 ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Logical                                                           &
     &  L_Physics   ! true if physics wanted
 
      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, timestep                                                        &
     &, G_term_tol       ! tolerance for vertical G term

!  VarRes horizontal co-ordinate spacing etc.
      Real                                                              &
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dlambda_u(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, dphi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)      &
     &, recip_dlamp(1-halo_i : row_length + halo_i)                     &
     &, recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

! time-weighting coefficients. See WP 154 for details.
      Real                                                              &
     &  alpha_1                                                         &
     &, alpha_2                                                         &
     &, alpha_3                                                         &
     &, alpha_4                                                         &
     &, alpha_Cd(boundary_layer_levels)

! trigonometric functions
      Real                                                              &
     &  sec_theta_latitude(1-off_x:row_length+off_x,1-off_y:rows+off_y) &
     &, cos_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
             ! Finite Volume cosine
     &, FV_sec_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)
             ! and secant arrays

      Real                                                              &
           ! components of coriolis force.
     &  f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, f3_at_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

      Real                                                              &
           ! primary model variables
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)  &
     &, v (1-off_x:row_length+off_x,1-off_y:n_rows+off_y, model_levels) &
     &, w (1-off_x:row_length+off_x,1-off_y:rows+off_y, 0:model_levels) &
     &, p (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)   &
     &, rho (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels) &
     &, rho_np1 (1-off_x:row_length+off_x,1-off_y:rows+off_y,           &
     &           model_levels)                                          &
     &, theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels) &
     &, theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, theta_np1(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
     &            model_levels)                                         &
     &, q (1-halo_i:row_length+halo_i,                                  &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, qcl (1-halo_i:row_length+halo_i,                                &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, qcf (1-halo_i:row_length+halo_i,                                &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, qcf2 (1-halo_i:row_length+halo_i,                               &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, qrain (1-halo_i:row_length+halo_i,                              &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, qgraup (1-halo_i:row_length+halo_i,                             &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, q_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          wet_model_levels)                                       &
     &, qcl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, qcf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, qcf2_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_model_levels)                                       &
     &, qrain_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, qgraup_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          wet_model_levels)                                       &
     &, q_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          wet_model_levels)                                       &
     &, qcl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, qcf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, qcf2_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, qrain_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_model_levels)                                       &
     &, qgraup_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)

      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)&
     &, R_v (1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels)&
     &, R_w (row_length, rows, model_levels)                            &
     &, exner_rho_levels (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
     &                    model_levels)                                 &
     &, frictional_timescale(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! GCR(k) arguments

      ! True if this is a dumping period
      Logical, Intent(In) :: ldump

      Logical                                                           &
     &  GCR_use_tol_abs                                                 &
     &, GCR_zero_init_guess                                             &
                             ! True if initial guess to solution is
                             ! zero.
     &, GCR_use_residual_Tol                                            &
     &, GCR_adi_add_full_soln                                           &
                              ! true then use full equation on RHS
                            ! on second and subsequent ADI timesteps
     &, L_gcr_fast_x        ! true then user faster non reproducible
                            !      code

      Real                                                              &
     &  GCR_tol_res                                                     &
     &, GCR_tol_abs                                                     &
     &, GCR_ADI_pseudo_timestep

      Integer                                                           &
     &  GCR_Restart_value                                               &
                           ! After how many iterations do we restart
     &, GCR_Diagnostics                                                 &
                           !
     &, GCR_max_iterations                                              &
     &, GCR_its_switch                                                  &
                            ! Iterations analysis switch
     &, GCR_its_avg_step(3)                                             &
                               ! Iterations analysis step now
     &, GCR_max_its                                                     &
                            ! Max iterations this period
     &, GCR_min_its                                                     &
                            ! Min iterations this period
     &, GCR_max_time                                                    &
                            ! Timestep number for max GCR its
     &, GCR_min_time                                                    &
                            ! Timestep number for min GCR its
     &, GCR_sum_its                                                     &
                           ! Sum iterations over test period

     &, GCR_precon_option                                               &
     &, GCR_n_ADI_pseudo_timesteps

! Physics arrays
      Real                                                              &
     &  rho_Km(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         0:boundary_layer_levels-1)                               &
     &, cH(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels-1)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  exner_prime (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               model_levels)                                      &
     &, dtheta_dr_term (row_length, rows, model_levels)

! Local Variables.

      Logical                                                           &
     &  L_poles      !  switch for use in etadot_calc, Global = T

      Integer                                                           &
     &  i, j, k                                                         &
                   ! Loop indices
     &, i_start, i_stop, i_stop_u                                       &
     &, j_start, j_stop, j_stop_v                                       &
     &, j_begin, j_end                                                  &
     &, solver_row_length                                               &
     &, number_dif_horiz_points                                         &
     &, info

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi                                                 &
     &, recip_timestep

      Real                                                              &
     &  weight1                                                         &
     &, weight2                                                         &
     &, recip_Cp                                                        &
     &, temp1                                                           &
     &, temp2                                                           &
     &, temp3                                                           &
     &, temp4                                                           &
     &, Ju                                                              &
     &, Jv                                                              &
     &, G_term

       Real                                                             &
                  !  VarRes
     &  wk_temp(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

! Helmholtz equation coefficients
      Real                                                              &
     & HM_Cxx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Czz (1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C2n(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C3 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C4 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C5 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &, HM_Cxz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_RHS(row_length, rows, model_levels)

      Real                                                              &
     &  Au (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &, Av (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &, Fu (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &, Fv (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &, recip_G_term (row_length, rows, model_levels-1)                 &
     &, K_term (row_length, rows, model_levels-1)

! 3-d work arrays
      Real                                                              &
     &  thetav_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, rescale (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &       wet_model_levels)                                          &
     &, moist_star (1-off_x:row_length+off_x,                           &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, moist_np1  (1-off_x:row_length+off_x,                           &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, weight_upper(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          &
     &, weight_lower(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          &
     &, recip_1p_cd_lambda( -off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, recip_1p_cd_phi(1-off_x:row_length+off_x,                       &
     &                   -off_y:rows+off_y, model_levels)

      Real,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     &  dry_density                                                     &
     &, dexprdr                                                         &
     &, thetav_star_rho                                                 &
     &, etadot                                                          &
     &, etadot_star                                                     &
     &, etadot_star2                                                    &
     &, v_star_minus_v                                                  &
     &, solver_lbc                                                      &
     &, L_solver_lbc

! 2-d work arrays for storing interpolated fields.
      Real                                                              &
     &  interp1(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp2(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp3(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp4(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp5(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp6(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, u_mask( -off_x:row_length+off_x, 1-off_y:rows+off_y)            &
     &, v_mask(1-off_x:row_length+off_x,  -off_y:rows+off_y)

       Real                                                             &
     &  sum_s(model_levels)                                             &
     &, sum_n(model_levels)                                             &
     &, temp_s(row_length,model_levels)                                 &
     &, temp_n(row_length,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! External Routines:
      External                                                          &
     &  GCR_k, Etadot_Calc

! ----------------------------------------------------------------------
! Section 1.   Calculate Residual in equation of state.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PE_HELMHOLTZ_GLOBAL_2B',zhook_in,zhook_handle)
      
      ALLOCATE ( dry_density(1-off_x:row_length+off_x,                  &
     &                   1-off_y:rows+off_y, model_levels) )
      ALLOCATE ( thetav_star_rho(1-off_x:row_length+off_x,              &
     &                   1-off_y:rows+off_y, model_levels) )

      ! Include condensate/hydrometeors in theta_v calculation
      If (L_qwaterload) Then
        moist      = qcl + qcf
        moist_star = qcl_star + qcf_star

        If (L_mcr_qcf2)then
          moist      = moist + qcf2
          moist_star = moist_star + qcf2_star
        End If
        If (L_mcr_qrain)then
          moist      = moist + qrain
          moist_star = moist_star + qrain_star
        End If
        If (L_mcr_qgraup) Then
          moist      = moist + qgraup
          moist_star = moist_star + qgraup_star
        End If
      Else
        moist = 0.0
        moist_star = 0.0
      End If 

! Number of points over which solver works
      solver_row_length = global_row_length
      L_poles = .true.
      number_dif_horiz_points = (global_rows-2) * global_row_length + 2

      i_start = 1
      i_stop = row_length
      i_stop_u = row_length
      j_start = 1
      j_stop = rows
      j_stop_v = n_rows
      j_begin = 1
      j_end = rows

      if (at_extremity(PSouth)) j_begin = 2
      if (at_extremity(PNorth)) j_end = rows - 1

      recip_delta_lambda = 1. / delta_lambda
      recip_delta_phi = 1. / delta_phi
      recip_timestep = 1. / timestep

! calculate mask which decides whether a v/u value substituted into
! the u/v equation is a boundary condition or not.
! Only 0 for LAMs on boundaries when using original lbc algorithm.
      u_mask = 1.0
      v_mask = 1.0

! ----------------------------------------------------------------------

! set up vertical interpolation weights between theta levels and rho
! levels
      Do k = 1, model_levels
        Do j = 0, rows + 1
          Do i = 0, row_length + 1
            weight_upper(i,j,k) = ( r_rho_levels(i,j,k) -               &
     &                              r_theta_levels(i,j,k-1) ) /         &
     &                            ( r_theta_levels(i,j,k) -             &
     &                              r_theta_levels(i,j,k-1) )
            weight_lower(i,j,k) = 1.0 - weight_upper(i,j,k)
          End Do
        End Do
      End Do  !  k = 1, model_levels

! Calculate thetav at theta levs, store in HM_Cxx2 as tempry. workspace
      weight1 = recip_epsilon
      Do k = 1, wet_model_levels
        Do j = 0, rows + 1
          Do i = 0, row_length + 1
            temp1 = ( 1. +  ( (weight1 - 1.)* q(i,j,k) ) - moist(i,j,k))
            dry_density(i,j,k) = ( 1. + ((weight1 - 1.)* q_star(i,j,k)) &
     &                    - moist_star(i,j,k) )
            HM_Cxx2(i,j,k) = theta(i,j,k) * temp1
            thetav_star(i,j,k) = theta_star(i,j,k) * dry_density(i,j,k)
          End Do
        End Do
      End Do  !  k = 1, wet_model_levels

! Interpolate thetav to rho levels, theta to rho levels, and calculate
! Residual in equation of state.
! theta at rho levels is stored in HM_Cyy2
! Residual in equation of state is stored in HM_RHS.

      recip_Cp = 1. / Cp

      k = 1
      Do j = 0, rows + 1
         Do i = 0, row_length + 1
! theta at rho levels
            HM_Cyy2(i,j,k) = theta(i,j,k)
! thetav at rho levels
            HM_Cxy2(i,j,k) = HM_Cxx2(i,j,k)
            thetav_star_rho(i,j,k) = thetav_star(i,j,k)
        End Do
      End Do
! set RHS at all points
! LAMs do not use boundaries but RHS is used to set first guess in GCR_k
      If (L_qwaterload) Then

        Do j = 1, rows
          Do i = 1, row_length
            HM_RHS(i,j,k) = kappa * rho(i,j,k) * exner_rho_levels(i,j,k)  &
     &                      *  thetav_star_rho(i,j,k)                     &
     &                      - r_rho_levels(i,j,k) * r_rho_levels(i,j,k)   &
     &                      * p(i,j,k) * recip_Cp
          End Do
        End Do

      Else !  L_qwaterload  FALSE

        Do j = 1, rows
          Do i = 1, row_length
            HM_RHS(i,j,k) = kappa * rho(i,j,k) * exner_rho_levels(i,j,k) *&
     &                         ( thetav_star_rho(i,j,k) + HM_Cxy2(i,j,k) *&
     &                   ( (1. - q(i,j,k)) / (1.-q_star(i,j,k)) - 1. ) ) -&
     &                         r_rho_levels(i,j,k) * r_rho_levels(i,j,k) *&
     &                                    p(i,j,k) * recip_Cp
          End Do
        End Do

      End If   !  L_qwaterload

      Do k = 2, model_levels

! (1-q)/(1-q_star) - 1 with q terms averaged to rho levels

          Do j = 1, rows
            Do i = 1, row_length
              interp3(i,j) = ( 1.-(weight_upper(i,j,k) * q(i,j,k) +     &
     &                             weight_lower(i,j,k)* q(i,j,k-1)) ) / &
     &                      ( 1.-(weight_upper(i,j,k) * q_star(i,j,k) + &
     &                      weight_lower(i,j,k) * q_star(i,j,k-1) ) ) - &
     &                        1.0
            End Do
          End Do

          Do j = 0, rows + 1
            Do i = 0, row_length + 1
              HM_Cyy2(i,j,k) = weight_upper(i,j,k) * theta(i,j,k) +     &
     &                         weight_lower(i,j,k) * theta(i,j,k-1)
! thetav at rho levels
              HM_Cxy2(i,j,k) = weight_upper(i,j,k) * HM_Cxx2(i,j,k) +   &
     &                         weight_lower(i,j,k) * HM_Cxx2(i,j,k-1)
              thetav_star_rho(i,j,k) = weight_upper(i,j,k) *            &
     &                                 thetav_star(i,j,k) +             &
     &                                 weight_lower(i,j,k) *            &
     &                                 thetav_star(i,j,k-1)
            End Do
          End Do

        If (L_qwaterload) Then
! set RHS at all points
! LAMs do not use boundaries but RHS is used to set first guess in GCR_k
          Do j = 1, rows
            Do i = 1, row_length
              HM_RHS(i,j,k) = kappa * rho(i,j,k)                        &
     &                        * exner_rho_levels(i,j,k)                 &
     &                    *  thetav_star_rho(i,j,k)                     &
     &                    - r_rho_levels(i,j,k) * r_rho_levels(i,j,k) * &
     &                      p(i,j,k) * recip_Cp
            End Do
          End Do

        Else  !  L_qwaterload  FALSE

          Do j = 1, rows
            Do i = 1, row_length
              HM_RHS(i,j,k) = kappa * rho(i,j,k) *                      &
     &                                        exner_rho_levels(i,j,k) * &
     &                      ( thetav_star_rho(i,j,k) + HM_Cxy2(i,j,k) * &
     &                                                 interp3(i,j) ) - &
     &                      r_rho_levels(i,j,k) * r_rho_levels(i,j,k) * &
     &                                 p(i,j,k) * recip_Cp
            End Do
          End Do

        End If ! L_qwaterload

      End Do  !  k = 2, model_levels

! Due to changes in u, v discretization theta_np1 is required
! instead of theta_star in H terms Cxy2, Cyx2, Cxx2, Cyy2, Cxp, Cyp
! Obtain corresponding thetav_star.
! Note that thetav_star in HM_RHS, initialised earlier,
! has intentionally not been replaced as it should not change.
      If ( L_new_tdisc .and. CycleNo > 1 ) Then

      ! Include condensate/hydrometeors in theta_v calculation
      If (L_qwaterload) Then
        moist_np1      = qcl_np1 + qcf_np1

        If (L_mcr_qcf2)then
          moist_np1      = moist_np1 + qcf2_np1
        End If
        If (L_mcr_qrain)then
          moist_np1      = moist_np1 + qrain_np1
        End If
        If (L_mcr_qgraup) Then
          moist_np1  = moist_np1 + qgraup_np1
        End If
      Else
        moist_np1 = 0.0
      End If

        k=1
          Do j = j_start - 1, j_stop + 1
            Do i = i_start - 1, i_stop + 1
              thetav_star(i,j,k) = theta_np1(i,j,k)                     &
     &                             * ( 1.+(weight1-1.)*q_np1(i,j,k)     &
     &                             - moist_np1(i,j,k) )
              thetav_star_rho(i,j,k) = thetav_star(i,j,k)
            End Do
          End Do

        Do k = 2, wet_model_levels
          Do j = j_start - 1, j_stop + 1
            Do i = i_start - 1, i_stop + 1
              thetav_star(i,j,k) = theta_np1(i,j,k)                     &
     &                             * ( 1.+(weight1-1.)*q_np1(i,j,k)     &
     &                             - moist_np1(i,j,k) )
              thetav_star_rho(i,j,k) = weight_upper(i,j,k) *            &
     &                                 thetav_star(i,j,k) +             &
     &                                 weight_lower(i,j,k) *            &
     &                                 thetav_star(i,j,k-1)
            End Do
          End Do
        End Do

      End If

! ----------------------------------------------------------------------
! Section 2.   Set up K, G, Au, Av, Fu and Fv terms of WP 154.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.1  Form K and G terms.
! ----------------------------------------------------------------------

! calculate G and K terms
      temp1 = Cp * alpha_2 * alpha_4 * timestep * timestep
      temp2 = alpha_4 * Cp * timestep
      temp4 = 1. / (G_term_tol - 1.)

      If ( .NOT. L_fint_theta ) Then

       If ( ( .NOT. L_new_tdisc ) .OR. CycleNo == 1 ) Then

        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            dtheta_dr_term(i,j,k) = ( theta_star(i,j,k+1) -             &
     &                                theta_star(i,j,k) ) /             &
     &                          ( r_theta_levels(i,j,k+1) -             &
     &                            r_theta_levels(i,j,k) )
          End Do
        End Do
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              dtheta_dr_term(i,j,k) = ( theta_star(i,j,k+1) -           &
     &                                  theta_star(i,j,k-1) ) /         &
     &                              ( r_theta_levels(i,j,k+1) -         &
     &                                r_theta_levels(i,j,k-1) )
            End Do
          End Do
        End Do ! k = 2, model_levels - 1

       Else

        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            dtheta_dr_term(i,j,k) =                                     &
     &                (theta_np1(i,j,k+1) - theta_np1(i,j,k)) /         &
     &                  (r_theta_levels(i,j,k+1) -                      &
     &                   r_theta_levels(i,j,k))
          End Do
        End Do
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              dtheta_dr_term(i,j,k) =                                   &
     &                (theta_np1(i,j,k+1) - theta_np1(i,j,k-1)) /       &
     &                 (r_theta_levels(i,j,k+1) -                       &
     &                  r_theta_levels(i,j,k-1))
            End Do
          End Do
        End Do

       End If ! If ( .NOT. L_new_tdisc ) .OR. CycleNo == 1

      Else ! L_fint_theta = true :

! 0 dtheta_dr_term if fully-interpolating SL theta advection is used
        Do k = 1, model_levels-1
          Do j = 1, rows
            Do i = 1, row_length
              dtheta_dr_term(i,j,k) = 0.0
            End Do
          End Do
        End Do

      End If ! If L_fint_theta

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
! dry_density converts theta to thetav
            temp3 =  - temp1 * (exner_rho_levels(i,j,k+1) -             &
     &                          exner_rho_levels(i,j,k)) /              &
     &                         ( r_rho_levels(i,j,k+1) -                &
     &                           r_rho_levels(i,j,k)) *                 &
     &                          dtheta_dr_term(i,j,k) *                 &
     &                          dry_density(i,j,k)
            G_term = 1. + temp3 / ( 1. + cH(i,j,k) * timestep )

            If ( G_term < G_term_tol ) Then
              cH(i,j,k) = ( temp3 * temp4 - 1.) * recip_timestep
              G_term = G_term_tol
            End If

            recip_G_term(i,j,k) = 1. / G_term
            K_term(i,j,k) = temp2 * thetav_star(i,j,k) *                &
     &                             recip_G_term(i,j,k)
            dtheta_dr_term(i,j,k) = dtheta_dr_term(i,j,k) /             &
     &                              ( 1. + cH(i,j,k) * timestep )
          End Do
        End Do
      End Do ! k = 1, model_levels - 1

! ----------------------------------------------------------------------
! Section 2.2  Form Cd terms.
!              Code assumes that the number of boundary layer levels
!              is at least 2.
! ----------------------------------------------------------------------

      Do k = 1, model_levels

          Do j = j_start - 1, j_stop + 1
            Do i = i_start - 2, i_stop + 1
              recip_1p_cd_lambda(i,j,k) = 1.0
            End Do
          End Do
          Do j = j_start - 2, j_stop + 1
            Do i = i_start - 1, i_stop + 1
              recip_1p_cd_phi(i,j,k) = 1.0
            End Do
          End Do


! ----------------------------------------------------------------------
! Section 2.3  Form Ju and hence Au and Fu terms.
! ----------------------------------------------------------------------

! calculate Ju and hence Au and Fu
! temp1 holds averaged term

        temp3 = alpha_3 * alpha_3 * timestep * timestep

! set whole domain values (needed to calculate or set u_star)
         Do j = j_start-1, j_stop+1
           Do i = i_start - 1, i_stop 
            temp1 = .25 * (v_mask(i+1,j) * recip_1p_cd_phi(i+1,j,k)     &
     &                  + v_mask(i,j) * recip_1p_cd_phi(i,j,k)          &
     &                  + v_mask(i+1,j-1) * recip_1p_cd_phi(i+1,j-1,k)  &
     &                  + v_mask(i,j-1)* recip_1p_cd_phi(i,j-1,k))
            Ju = 1. / recip_1p_cd_lambda(i,j,k)  + temp3 * temp1 *      &
     &            f3_at_u(i,j) * f3_at_u(i,j)

            Au(i,j,k) = 1. / Ju
            Fu(i,j,k) = alpha_3 * timestep * f3_at_u(i,j) / Ju
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.4  Form Jv and hence Av and Fv terms.
! ----------------------------------------------------------------------

! calculate Jv and hence Av and Fv
! temp1 holds averaged  term

        temp4 = alpha_3 * alpha_3 * timestep * timestep

! set whole domain values (needed to calculate or set v_star)
         Do j = j_start - 1, j_stop
           Do i = i_start - 1, i_stop + 1 
            temp1 = .25 * (u_mask(i-1,j) * recip_1p_cd_lambda(i-1,j,k)  &
     &                  + u_mask(i,j) * recip_1p_cd_lambda(i,j,k)       &
     &               + u_mask(i-1,j+1) * recip_1p_cd_lambda(i-1,j+1,k)  &
     &                  + u_mask(i,j+1)* recip_1p_cd_lambda(i,j+1,k))

            Jv = 1. / recip_1p_cd_phi(i,j,k) + temp4 * temp1 *          &
     &              f3_at_v(i,j) * f3_at_v(i,j)
            Av(i,j,k) = 1. / Jv
            Fv(i,j,k) = alpha_3 * timestep * f3_at_v(i,j) / Jv
          End Do
        End Do

      End Do  !  k = 1, model_levels

! ----------------------------------------------------------------------
! Section 3    Set up Helmholtz equation coefficients and
!              right-hand-side.
! ----------------------------------------------------------------------

      ALLOCATE ( etadot(row_length, rows, 0:model_levels) )
      ALLOCATE ( etadot_star(row_length, rows, 0:model_levels) )
      ALLOCATE ( etadot_star2(row_length, rows, 0:model_levels) )
      ALLOCATE ( v_star_minus_v(1-off_x:row_length+off_x,               &
     &                   1-off_y:n_rows+off_y, model_levels) )
! ----------------------------------------------------------------------
! Section 3.1  Re-scale Helmholtz equation right-hand-side and
!              set C3 and C4 coefficients.
! ----------------------------------------------------------------------
! at the moment HM_RHS holds residual in equation of state,
! HM_Cxy2 holds thetav at rho levels, HM_Cyy2 holds theta at rho levels.
! Store dry rho in dry_density

! store d(r)/d(eta) about rho levels in HM_Czz
! In CycleNo > 1 dry_density needs to be defined in terms of rho^(1)
! for the improved rho-discretization to take place. Rho^(1) is an
! alpha1 weighted t-average between rho and rho_np1. dry_density
! is used in Helmholtz coefficients Cxx1, Cyy1, Czz, C5.

      Do k = 1, model_levels
        Do j = j_start - 1, j_stop + 1
          Do i = i_start - 1, i_stop + 1
            HM_Czz(i,j,k) = (r_theta_levels(i,j,k) -                    &
     &                       r_theta_levels(i,j,k-1)) /                 &
     &                     (eta_theta_levels(k) - eta_theta_levels(k-1))
            rescale(i,j,k) = recip_timestep * HM_Czz(i,j,k) *           &
     &                                           rho(i,j,k)
            End Do
          End Do
      End Do  ! k = 1, model_levels

      k = 1

      If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then

        Do j = j_start, j_stop
           Do i = i_start - 1, i_stop
            interp6(i,j) = dry_density(i,j,k) ! cajm2
          End Do
        End Do
        Do j = j_start - 1, j_stop + 1
           Do i = i_start - 1, i_stop + 1
            dry_density(i,j,k) = rho(i,j,k) *                           &
     &                           (1. - q(i,j,k) - moist(i,j,k))
          End Do
        End Do

      Else  !  CycleNo > 1 .AND. L_new_tdisc

        Do j = j_begin - 1, j_end + 1
          Do i = i_start-1, i_stop + 1
            interp6(i,j) = dry_density(i,j,k) ! cajm2
            dry_density(i,j,k)=(1. - alpha_1) * rho(i,j,k) *            &
     &                         (1. - q(i,j,k) - moist(i,j,k) ) +        &
     &                                alpha_1 * rho_np1(i,j,k) *        &
     &                         (1. - q_np1(i,j,k) - moist_np1(i,j,k) )
          End Do
        End Do

      End If  !  CycleNo == 1 .OR. ( .NOT. L_new_tdisc )

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            temp2 = 1./(1. - q_star(i,j,k)-moist_star(i,j,k))
            temp3 = HM_Czz(i,j,k)                                       &
     &              /(HM_Cxy2(i,j,k) * timestep * kappa *               &
     &                exner_rho_levels(i,j,k) * temp2)
            HM_RHS(i,j,k) = - HM_RHS(i,j,k) * temp3
            HM_C4(i,j,k) = ( r_rho_levels(i,j,k) *                      &
     &                       r_rho_levels(i,j,k) * p(i,j,k) /           &
     &                       (R * exner_rho_levels(i,j,k)) -            &
     &                       kappa * HM_Cxy2(i,j,k) * rho(i,j,k) )      &
     &                     * temp3
            HM_C3(i,j,k) = HM_Czz(i,j,k) * rho(i,j,k)*interp6(i,j)      &
     &                     / (temp2 * HM_Cxy2(i,j,k))
          End Do
        End Do

      Do k = 2, model_levels

        If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp6(i,j) =  dry_density(i,j,k) ! cajm2
                interp2(i,j) = weight_upper(i,j,k) *                    &
     &                        (1./ (1. - q_star(i,j,k)                  &
     &                        - moist_star(i,j,k) ) ) +                 &
     &                         weight_lower(i,j,k) *                    &
     &                        (1./ (1. - q_star(i,j,k-1)                &
     &                        - moist_star(i,j,k-1) ) )
            End Do
          End Do
          Do j = j_start - 1, j_stop + 1
            Do i = i_start - 1, i_stop + 1
                dry_density(i,j,k) = rho(i,j,k) *                       &
     &                           ( weight_upper(i,j,k)                  &
     &                           * (1. - q(i,j,k) - moist(i,j,k) ) +    &
     &                           weight_lower(i,j,k) *                  &
     &                           (1. - q(i,j,k-1) - moist(i,j,k-1) ) )
            End Do
          End Do
        Else  !  CycleNo > 1 .AND. L_new_tdisc

            Do j = j_begin - 1, j_end + 1
              Do i = i_start - 1, i_stop + 1
                interp6(i,j) =  dry_density(i,j,k) ! cajm2
                interp2(i,j) =                                          &
     &               weight_upper(i,j,k) *                              &
     &               (1./ (1. - q_star(i,j,k) - moist_star(i,j,k) ) )   &
     &               + weight_lower(i,j,k) *                            &
     &               (1./ (1. - q_star(i,j,k-1) - moist_star(i,j,k-1)))
                dry_density(i,j,k) =                                    &
     &                     (1.-alpha_1) * rho(i,j,k) *                  &
     &                     ( weight_upper(i,j,k) *                      &
     &                     (1. - q(i,j,k) - moist(i,j,k) ) +            &
     &                     weight_lower(i,j,k) *                        &
     &                     (1. - q(i,j,k-1)-moist(i,j,k-1)) ) +         &
     &                     alpha_1 * rho_np1(i,j,k) *                   &
     &                     ( weight_upper(i,j,k) *                      &
     &                     (1. - q_np1(i,j,k)-moist_np1(i,j,k))+        &
     &                     weight_lower(i,j,k) *                        &
     &                     (1. - q_np1(i,j,k-1)-moist_np1(i,j,k-1)))
              End Do
            End Do

        End If  !  CycleNo == 1 .OR. ( .NOT. L_new_tdisc )

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            temp3 = HM_Czz(i,j,k)                                       &
     &              /(HM_Cxy2(i,j,k) * timestep * kappa *               &
     &                exner_rho_levels(i,j,k) * interp2(i,j))
            HM_RHS(i,j,k) = - HM_RHS(i,j,k) * temp3
            HM_C4(i,j,k) = ( r_rho_levels(i,j,k) *                      &
     &                       r_rho_levels(i,j,k) * p(i,j,k) /           &
     &                       (R * exner_rho_levels(i,j,k)) -            &
     &                       kappa * HM_Cxy2(i,j,k) * rho(i,j,k) )      &
     &                     * temp3
            HM_C3(i,j,k) = HM_Czz(i,j,k) * rho(i,j,k)*interp6(i,j)      &
     &                     / (interp2(i,j) * HM_Cxy2(i,j,k))
          End Do
        End Do
      End Do  !  k = 2, model_levels

! Calculate dry rho at theta levels
      Do k = 1, model_levels - 1
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            HM_C5(i,j,k) = ( dry_density(i,j,k) *                       &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      dry_density(i,j,k+1) *                      &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) )
          End Do
        End Do
      End Do !  k = 1, model_levels - 1

! ----------------------------------------------------------------------
! Section 3.2  Lambda direction coefficients (Cxx1,Cxx2,etc).
! ----------------------------------------------------------------------
! Initialise interp6 to zero at poles
      If(at_extremity(PSouth))then
        Do i = 0, row_length + 1
          interp6(i,1) = 0.0
        End Do
      End If
      If(at_extremity(PNorth))then
        Do i = 0, row_length + 1
          interp6(i,rows) = 0.0
        End Do
      End If

! calculate HM_Cxx1, HM_Cxx2 and HM_Cx terms.
! Calculate terms requiring averaging to u points.
! average r and store in interp2 as 1/r
! average dry rho * d(r)/d(eta) and store in interp3
! average R_v and store in interp4

      Do k = 1, model_levels

          temp1 = recip_delta_lambda * 2.
          Do j = j_start - 1, j_stop + 1
            Do i = i_start - 1, i_stop
              interp2(i,j) = temp1 / (r_rho_levels(i+1,j,k) +           &
     &                              r_rho_levels(i,j,k))
            End Do
          End Do

          Do j = j_start, j_stop
            Do i = i_start - 1, i_stop
              interp3(i,j) = (dry_density(i+1,j,k) * HM_Czz(i+1,j,k) +  &
     &                        dry_density(i,j,k) * HM_Czz(i,j,k)) * .5
              HM_Cxx1(i,j,k) = interp3(i,j) * interp2(i,j)
            End Do
          End Do

! calculate HM_Cxx2 and HM_Cx terms.

        temp1 = alpha_1 * alpha_3 * timestep
        temp2 = Cp * 0.5
        Do j = j_start - 1, j_stop + 1
          Do i = i_start - 1, i_stop
            HM_Cyx2(i,j,k) = temp2 * (thetav_star_rho(i,j,k) +          &
     &                                thetav_star_rho(i+1,j,k) ) *      &
     &                       interp2(i,j)                               &
     &                       * FV_sec_theta_latitude(i,j)
            HM_Cxx2(i,j,k) = temp1 * Au(i,j,k) * HM_Cyx2(i,j,k)
          End Do
        End Do

        Do j = j_start, j_stop
          Do i = i_start - 1, i_stop
            HM_Cxy1(i,j,k) = temp1 * Fu(i,j,k)
          End Do
        End Do

! Form u_star, store in interp6

         Do j = j_begin, j_end
          Do i = i_start-1, i_stop
            interp6(i,j) =u(i,j,k)+ alpha_1 *( Au(i,j,k) * R_u(i,j,k)   &
     &                   + Fu(i,j,k) * .25 *                            &
     &                    (R_v(i+1,j,k) * recip_1p_cd_phi(i+1,j,k)      &
     &                   + R_v(i,j,k) * recip_1p_cd_phi(i,j,k)          &
     &                   + R_v(i+1,j-1,k) * recip_1p_cd_phi(i+1,j-1,k)  &
     &                   + R_v(i,j-1,k) * recip_1p_cd_phi(i,j-1,k) ))
          End Do
        End Do

! Modify right-hand-side with respect to the terms similar to the
! HM_Cxx2 and HM_Cx coefficients

         Do j = j_begin, j_end
           Do i = i_start, i_stop
            HM_RHS(i,j,k) = HM_RHS(i,j,k) + sec_theta_latitude(i,j) *   &
     &                        (interp3(i,j) * interp6(i,j)              &
     &                         * interp2(i,j) -                         &
     &                         interp3(i-1,j) * interp6(i-1,j)          &
     &                         * interp2(i-1,j))
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 3.3  Phi direction coefficients (Cyy1,etc) and add terms to
!              right hand side.
! ----------------------------------------------------------------------
! calculate HM_Cyy1, HM_Cyy2 and HM_Cyx1, HM_Cxy2 terms.

! Calculate terms requiring averaging to v points.
! average r and store in interp2 as 1/r
! average dry rho * d(r)/d(eta) and store in interp3
! average R_u in both lamda and phi and store in interp4

        temp1 = recip_delta_phi * 2
        Do j = j_start - 1, j_stop
          Do i = i_start - 1, i_stop + 1
            interp2(i,j) = temp1 / (r_rho_levels(i,j+1,k) +             &
     &                              r_rho_levels(i,j,k))
          End Do
        End Do
        Do j = j_start - 1, j_stop
          Do i = i_start, i_stop
            interp3(i,j) = (dry_density(i,j+1,k) * HM_Czz(i,j+1,k) +    &
     &                      dry_density(i,j,k) * HM_Czz(i,j,k)) * .5
            HM_Cyy1(i,j,k) = interp3(i,j) * interp2(i,j)                &
     &                       * cos_v_latitude(i,j)
          End Do
        End Do


        temp1 = alpha_1 * alpha_3 * timestep
        temp2 = Cp * 0.5
        Do j = j_start - 1, j_stop
          Do i = i_start - 1, i_stop + 1
            HM_Cxy2(i,j,k) = temp2 * (thetav_star_rho(i,j+1,k) +        &
     &                                thetav_star_rho(i,j,k) ) *        &
     &                       interp2(i,j)
            HM_Cyy2(i,j,k) = temp1 * Av(i,j,k) * HM_Cxy2(i,j,k)
          End Do
        End Do

        Do j = j_start - 1, j_stop
           Do i = i_start, i_stop
            HM_Cyx1(i,j,k) = temp1 * Fv(i,j,k)
          End Do
        End Do

! Form v_star, store in HM_Cz

        Do j = j_start-1, j_stop
          Do i = i_start, i_stop
            HM_Cz(i,j,k)= v(i,j,k)+ alpha_1 * (Av(i,j,k) * R_v(i,j,k)   &
     &                   - Fv(i,j,k) * 0.25 *                           &
     &                (R_u(i-1,j+1,k) * recip_1p_cd_lambda(i-1,j+1,k)   &
     &               + R_u(i,j+1,k) * recip_1p_cd_lambda(i,j+1,k)       &
     &               + R_u(i-1,j,k) * recip_1p_cd_lambda(i-1,j,k)       &
     &               + R_u(i,j,k)* recip_1p_cd_lambda(i,j,k)) )
          End Do
        End Do

! Modify right-hand-side with respect to the terms similar to the
! HM_Cyy2 coefficient

        Do j = j_begin, j_end
          Do i = i_start, i_stop
            HM_RHS(i,j,k) = HM_RHS(i,j,k)                               &
     &                    +  FV_sec_theta_latitude(i,j) *               &
     &                      ( HM_Cz(i,j,k) * interp3(i,j) *             &
     &                        interp2(i,j) * cos_v_latitude(i,j) -      &
     &                        HM_Cz(i,j-1,k) * interp3(i,j-1) *         &
     &                        interp2(i,j-1) * cos_v_latitude(i,j-1) )
          End Do
        End Do

! Modify right hand side at the poles
! save terms at poles for summing.

          If(at_extremity(PSouth))then
            Do i = 1, row_length
              temp_s(i,k) = HM_Cz(i,1,k) * interp3(i,1) *               &
     &                      interp2(i,1)
            End Do
          End If
          If(at_extremity(PNorth))then
            Do i = 1, row_length
              temp_n(i,k) = HM_Cz(i,n_rows,k) * interp3(i,n_rows) *     &
     &                      interp2(i,n_rows)
            End Do
          End If

! subtract u from u_star (stored in interp6) (store in HM_Czz)
!     loop changed to agree with old code  Do j = j_start, j_stop
          Do j = j_begin, j_end
            Do i = i_start - 1, i_stop
              HM_Czz(i,j,k) = interp6(i,j) - u(i,j,k)
            End Do
          End Do

! subtract v from v_star (stored in v_star_minus_v)
            Do j = j_start - 1, j_stop
              Do i = i_start, i_stop
                v_star_minus_v(i,j,k) = HM_Cz(i,j,k) - v(i,j,k)
              End Do
            End Do

      End Do  !  k = 1, model_levels

! sum terms at poles and add on
        If( at_extremity(PSouth)) Then

          CALL global_2d_sums(temp_s, row_length, 1, 0, 0, model_levels, &
                              sum_s, proc_row_group)

            Do k = 1, model_levels
              sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) *        &
     &                   cos_v_latitude(1,1) / global_row_length
              Do i = 1, row_length
                HM_RHS(i,1,k) = HM_RHS(i,1,k) + sum_s(k)
              End Do
            End Do  ! k = 1, model_levels

        End If ! at_extremity(PSouth)

        If( at_extremity(PNorth)) Then

          CALL global_2d_sums(temp_n, row_length, 1, 0, 0, model_levels, &
                              sum_n, proc_row_group)

            Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) *     &
     &                     cos_v_latitude(1,n_rows) / global_row_length
              Do i = 1, row_length
                HM_RHS(i,rows,k) = HM_RHS(i,rows,k) - sum_n(k)
              End Do
            End Do  ! k = 1, model_levels

        End If ! at_extremity(PNorth)

! ----------------------------------------------------------------------
! Section 3.4  r direction coefficients (C5,Cz,Czz,Cxz,Cyz) and
!              add terms to right hand side.
! ----------------------------------------------------------------------

! Calculate etadot_star, the known part of the increment to etadot.
! HM_Czz currently holds u_star-u
! v_star_minus_v currently holds v_star-v
! store w term in etadot

        Do k = 1, model_levels-1
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              etadot_star2(i,j,k) = alpha_2 * R_w(i,j,k) *              &
     &                            recip_G_term(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: etadot_calc
      Call Etadot_Calc (r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  HM_Czz, v_star_minus_v, etadot_star2,           &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  dlambda_p, dphi_p,                              &
     &                  wt_lambda_p, wt_lambda_u,                       &
     &                  wt_phi_p, wt_phi_v,                             &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  proc_row_group, at_extremity, global_row_length,&
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  off_x, off_y, 0, 0, 0, 0,                       &
     &                  i_start, i_stop,                                &
     &                  j_start, j_stop, j_begin, j_end,                &
     &                  L_regular, L_poles,                             &
     &                  temp_s, temp_n, etadot_star)

! Calculate etadot.

! DEPENDS ON: etadot_calc
      Call Etadot_Calc (r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u, v, w,                                        &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  dlambda_p, dphi_p,                              &
     &                  wt_lambda_p, wt_lambda_u,                       &
     &                  wt_phi_p, wt_phi_v,                             &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  proc_row_group, at_extremity, global_row_length,&
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  off_x, off_y, off_x, off_y, 0, 0,               &
     &                  i_start, i_stop,                                &
     &                  j_start, j_stop, j_begin, j_end,                &
     &                  L_regular, L_poles,                             &
     &                  temp_s, temp_n, etadot)

! HM_Czz currently holds u_star-u
! v_star_minus_v currently holds v_star-v

      Do k = 1, model_levels

! calculate terms required for vertical averaging and differencing
! store the terms for differencing in interp1&2 where interp1 is the
! value at the level above and interp2 the value at the level below.
! store the terms for averaging in interp3&4 where interp3 is the
! value at the level above and interp4 the value at the level below.

! set the coefficients HM_Czz and HM_Cz

        If ( k == 1 ) Then

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp2(i,j) = 0.
            End Do
          End Do

          Do j = j_start, j_stop
            Do i = i_start, i_stop
! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
              interp5(i,j) = alpha_2 * temp2
              interp1(i,j) = HM_C5(i,j,k) *                             &
     &                       ( etadot(i,j,k) + etadot_star(i,j,k) ) /   &
     &                         temp2

            End Do
          End Do

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp4(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)
              interp3(i,j) = interp4(i,j)
              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = HM_C5(i,j,k) * temp1
              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)
            End Do
          End Do

        Else If ( k == model_levels ) Then

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp2(i,j) = interp1(i,j)
              interp1(i,j) = 0.
              interp4(i,j) = interp3(i,j)
              interp3(i,j) = 0.
              HM_Cz(i,j,k) = 0.
              HM_Czz(i,j,k) = 0.
            End Do
          End Do

        Else If ( k == model_levels - 1 ) Then

          Do j = j_start, j_stop
            Do i = i_start, i_stop
! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
              interp2(i,j) = interp1(i,j)
              interp1(i,j) = HM_C5(i,j,k) *                             &
     &                       ( etadot(i,j,k) +  etadot_star(i,j,k) ) /  &
     &                         temp2
              interp5(i,j) = alpha_2 * temp2

            End Do
          End Do

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp4(i,j) = interp3(i,j)
              interp3(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)
              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = HM_C5(i,j,k) * temp1
              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)
            End Do
          End Do

        Else  !  1 < k < model_levels - 1

          Do j = j_start, j_stop
            Do i = i_start, i_stop
! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
              interp2(i,j) = interp1(i,j)
              interp1(i,j) = HM_C5(i,j,k) *                             &
     &                       ( etadot(i,j,k) +  etadot_star(i,j,k) ) /  &
     &                         temp2
              interp5(i,j) = alpha_2 * temp2
            End Do
          End Do

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              interp4(i,j) = interp3(i,j)
              interp3(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)
              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = HM_C5(i,j,k) * temp1
              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)
            End Do
          End Do

        End If  ! k = 1

! Modify right-hand-side

        temp1 = 1. / (eta_theta_levels(k) - eta_theta_levels(k-1))

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            HM_RHS(i,j,k) = HM_RHS(i,j,k) +                             &
     &                      (interp1(i,j) - interp2(i,j)) * temp1 +     &
     &                                               HM_C3(i,j,k) *     &
     &                        (weight_upper(i,j,k) * interp3(i,j) +     &
     &                         weight_lower(i,j,k) * interp4(i,j) )
          End Do
        End Do

      End Do  !  k = 1, model_levels

      temp1 = recip_delta_lambda * 2.
      temp2 = recip_delta_phi * 2.0

      Do k = 1, first_constant_r_rho_level - 1
! Calculate HM_Cxz and HM_Cyz
! Calculate dr/d lambda

        Do j = j_start - 1, j_stop + 1
          Do i = i_start - 1, i_stop
            HM_Cxp(i,j,k) =  ( r_rho_levels(i+1,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) * 2.0 /            &
     &                       (thetav_star_rho(i+1,j,k) +                &
     &                        thetav_star_rho(i,j,k) )
          End Do
        End Do
        Do j = j_start, j_stop
          Do i = i_start - 1, i_stop
            HM_Cxz(i,j,k) =  ( r_theta_levels(i+1,j,k) -                &
     &                         r_theta_levels(i,j,k) )                  &
     &                          * temp1                                 &
     &                          * sec_theta_latitude(i,j)               &
     &                          / ( r_theta_levels(i+1,j,k) +           &
     &                              r_theta_levels(i,j,k) )
          End Do
        End Do

! Calculate dr/d phi
        Do j = j_start-1, j_stop
          Do i = i_start - 1, i_stop + 1
            HM_Cyp(i,j,k) = (r_rho_levels(i,j+1,k) -                    &
     &                       r_rho_levels(i,j,k) ) * 2.0 /              &
     &                      (thetav_star_rho(i,j+1,k) +                 &
     &                       thetav_star_rho(i,j,k) )
          End Do
        End Do

        Do j = j_start-1, j_stop
          Do i = i_start, i_stop
            HM_Cyz(i,j,k) = (r_theta_levels(i,j+1,k) -                  &
     &                         r_theta_levels(i,j,k) )                  &
     &                         * temp2 /                                &
     &                        (r_theta_levels(i,j+1,k) +                &
     &                         r_theta_levels(i,j,k) )
          End Do
        End Do

      End Do ! k = 1, first_constant_r_rho_level - 1

      DEALLOCATE ( dry_density )
      DEALLOCATE ( thetav_star_rho )
      DEALLOCATE ( etadot )
      DEALLOCATE ( etadot_star )
      DEALLOCATE ( etadot_star2 )
      DEALLOCATE ( v_star_minus_v )

! ----------------------------------------------------------------------
! Section 4.   Call elliptic solver to solve Helmholtz equation.
! ----------------------------------------------------------------------

! re-scale HM_cyx2 and HM_cxy2 to incorporate friction in solver
      Do k = 1, boundary_layer_levels

        Do j = j_begin, j_end
          Do i = i_start, i_stop
            HM_Cyx2(i,j,k) = HM_Cyx2(i,j,k)                             &
     &                       * recip_1p_cd_lambda(i,j,k)
          End Do
        End Do

        Do j = j_begin - 1, j_end
          Do i = i_start, i_stop
            HM_Cxy2(i,j,k) = HM_Cxy2(i,j,k)                             &
     &                       * recip_1p_cd_phi(i,j,k)
          End Do
        End Do
      End Do

      If(at_extremity(PSouth))then
        Do k = 1, model_levels
          Do i = 0, row_length+1
            thetav_star(i,0,k) = 0.
          End Do
        End Do
      End If
      If(at_extremity(PNorth))then
        Do k = 1, model_levels
          Do i = 0, row_length+1
            thetav_star(i,rows+1,k) = 0.
          End Do
        End Do
      End If

      Do k = 1, first_constant_r_rho_level - 1
        Do j= j_start - 1, j_stop + 1
          Do i= i_start - 1, i_stop + 1
            HM_C2n(i,j,k) = thetav_star(i,j,k) /                        &
     &                     (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
          End Do
        End Do
      End Do ! k = 1, first_constant_r_rho_level - 1

! initialise parts of coefficients not set so far to zero.

        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1-off_x, row_length+off_x
              HM_Cxx1(i,1,k) = 0.
              HM_Cyx2(i,1,k) = 0.
              HM_Cxy1(i,1,k) = 0.
              HM_Cxx2(i,1,k) = 0.
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1-off_x, row_length+off_x
              HM_Cxz(i,1,k) = 0.0
              HM_Cxp(i,1,k) = 0.0
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1-off_x, row_length+off_x
              HM_Cxx1(i,rows,k) = 0.
              HM_Cxx2(i,rows,k) = 0.
              HM_Cxy1(i,rows,k) = 0.
              HM_Cyx2(i,rows,k) = 0.
              HM_Cyy1(i,rows,k) = 0.
              HM_Cyy2(i,rows,k) = 0.
              HM_Cxy2(i,rows,k) = 0.
              HM_Cyx1(i,rows,k) = 0.
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1-off_x, row_length+off_x
              HM_Cxz(i,rows,k) = 0.0
              HM_Cxp(i,rows,k) = 0.0
            End Do
          End Do
        End If

      if ( L_print_L2helm ) then

      if( norm_start_lev == norm_end_lev ) then
        do k = 1, model_levels

! DEPENDS ON: helm_l2_print
        Call helm_l2_print(                                             &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,               &
     &                     HM_Cxp, HM_Cyp, HM_C2n,                      &
     &                     HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
     &                     HM_C3, HM_C4, HM_C5, HM_RHS,                 &
     &                     row_length, rows, model_levels,              &
     &                     first_constant_r_rho_level_m1,               &
     &                     k, k, off_x, off_y,                          &
     &                     .false., .false., 1, L_print_allpe, me )
        end do ! k = 1, model_levels
      else ! norms from norm_start_lev to norm_end_lev
! DEPENDS ON: helm_l2_print
        Call helm_l2_print(                                             &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,               &
     &                     HM_Cxp, HM_Cyp, HM_C2n,                      &
     &                     HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
     &                     HM_C3, HM_C4, HM_C5, HM_RHS,                 &
     &                     row_length, rows, model_levels,              &
     &                     first_constant_r_rho_level_m1,               &
     &                     norm_start_lev, norm_end_lev, off_x, off_y,  &
     &                     .false., .false., 1, L_print_allpe, me )
      end if   ! norm_start_lev == norm_end_lev
! DEPENDS ON: um_fort_flush
        if ( L_flush ) call UM_FORT_FLUSH(6,info)

      endif ! L_print_L2helm

! Call GCR(k) to solve problem.
! DEPENDS ON: gcr_k_2B
      CALL GCR_k_2B(                                                    &
     &                 row_length, rows, n_rows, model_levels,          &
     &                 model_domain, timestep_number, CycleNo,          &
     &                 GCR_Diagnostics, GCR_max_its, GCR_min_its,       &
     &                 GCR_sum_its, GCR_its_switch, GCR_its_avg_step,   &
     &                 GCR_max_time, GCR_min_time,                      &
     &                 GCR_max_iterations, GCR_tol_res,                 &
     &                 GCR_Restart_value, GCR_tol_abs,                  &
     &                 GCR_use_tol_abs,                                 &
     &                 GCR_zero_init_guess,                             &
     &                 GCR_use_residual_Tol,                            &
     &                 GCR_precon_option, GCR_ADI_Pseudo_timestep,      &
     &                 GCR_n_ADI_pseudo_timesteps,                      &
     &                 GCR_adi_add_full_soln, L_gcr_fast_x, L_regular,  &
     &                 first_constant_r_rho_level,                      &
     &                 first_constant_r_rho_level_m1,                   &
     &                 lambda_p, phi_p, lambda_u, phi_v,                &
     &                 dlambda_p, dphi_p, dlambda_u, dphi_v,            &
     &                 recip_dlamp, recip_dphip,                        &
     &                 recip_dlamu, recip_dphiv,                        &
     &                 wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,    &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 FV_sec_theta_latitude, FV_cos_theta_latitude,    &
     &                 HM_Cxx1, HM_Cxx2,                                &
     &                 HM_Cyy1, HM_Cyy2,                                &
     &                 HM_Czz, HM_Cz, HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,   &
     &                 HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,              &
     &                 HM_C2n, HM_C3, HM_C4, HM_C5,                     &
     &                 HM_RHS, rescale,                                 &
     &                 weight_upper, weight_lower, exner_prime,         &
     &                 off_x, off_y, me, n_proc, n_procx, n_procy,      &
     &                 at_extremity, l_datastart,                       &
     &                 proc_row_group, proc_col_group,                  &
     &                 global_row_length, global_rows,                  &
     &                 number_dif_horiz_points, halo_i, halo_j,         &
     &                 i_start, i_stop, j_start, j_stop,                &
     &                 j_begin, j_end, solver_row_length,               &
     &                 g_rows, g_row_length,                            &
     &                 ldump                                            &
     &                 )

! DEPENDS ON: swap_bounds
      call swap_bounds(                                                 &
     &                 exner_prime, row_length, rows, model_levels,     &
     &                 off_x, off_y, fld_type_p, .false.)

! un-re-scale HM_cyx2 and HM_cxy2 to remove friction term
      Do k = 1, boundary_layer_levels
        Do j = j_begin, j_end + 1
          Do i = i_start - 1, i_stop
            HM_Cyx2(i,j,k) = HM_Cyx2(i,j,k) / recip_1p_cd_lambda(i,j,k)
          End Do
        End Do
        Do j = j_begin - 1,j_end
          Do i = i_start, i_stop_u + 1
            HM_Cxy2(i,j,k) = HM_Cxy2(i,j,k) / recip_1p_cd_phi(i,j,k)
          End Do
        End Do
      End Do  ! k = 1, boundary_layer_levels

! ----------------------------------------------------------------------
! Section 5.   Calculate increments to u,v and w.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 5.1  Calculate increments to w.
! ----------------------------------------------------------------------

      ALLOCATE ( dexprdr(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels - 1) )

      Do k = 1, model_levels - 1

        Do j = j_start - 1, j_stop + 1
          Do i = i_start - 1, i_stop + 1
            dexprdr(i,j,k) = (exner_prime(i,j,k+1) - exner_prime(i,j,k))&
     &                       / ( r_rho_levels(i,j,k+1) -                &
     &                           r_rho_levels(i,j,k) )
          End Do
        End Do

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            R_w(i,j,k) = R_w(i,j,k) * recip_G_term(i,j,k) -             &
     &                      K_term(i,j,k) * dexprdr(i,j,k)
          End Do
        End Do

      End Do ! k = 1, model_levels - 1

! ----------------------------------------------------------------------
! Section 5.2  Calculate increments to u and v.
! ----------------------------------------------------------------------

        Do  k = 1, model_levels
! calculate first derivatives times Cp theta /r.
          If ( k >= first_constant_r_rho_level ) Then

            Do j = j_begin, j_end + 1
              Do i = i_start - 1, i_stop
                interp5(i,j) =  (exner_prime(i+1,j,k) -                 &
     &                           exner_prime(i,j,k)) *                  &
     &                           HM_Cyx2(i,j,k)
              End Do
            End Do

            Do j = j_begin - 1, j_end
              Do i = i_start, i_stop_u + 1
                interp6(i,j) =  (exner_prime(i,j+1,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                          * HM_Cxy2(i,j,k)
              End Do
            End Do

          Else ! k < first_constant_r_rho_level

! form average vertical pressure derivative
            If ( k > 1 ) Then

              Do j = j_begin - 1, j_end + 1
                Do i = i_start - 1, i_stop + 1
                  interp1(i,j) = weight_upper(i,j,k) *                  &
     &                           thetav_star(i,j,k) * dexprdr(i,j,k)    &
     &                         + weight_lower(i,j,k)                    &
     &                           * thetav_star(i,j,k-1)                 &
     &                           * dexprdr(i,j,k-1)
                End Do
              End Do

            Else  !  k = 1

              Do j = j_begin - 1, j_end + 1
                Do i = i_start - 1, i_stop + 1
                  interp1(i,j) = weight_upper(i,j,k) *                  &
     &                           thetav_star(i,j,k) * dexprdr(i,j,k)
                End Do
              End Do
            End If  ! k > 1

            Do j = j_begin, j_end + 1
              Do i = i_start - 1, i_stop
                interp5(i,j) = ((exner_prime(i+1,j,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                           - HM_Cxp(i,j,k) * .5*                  &
     &                             (interp1(i,j)+interp1(i+1,j)) )      &
     &                          * HM_Cyx2(i,j,k)
              End Do
            End Do

            Do j = j_begin - 1, j_end
              Do i = i_start, i_stop_u + 1
                interp6(i,j) = ((exner_prime(i,j+1,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                           - HM_Cyp(i,j,k) * .5*                  &
     &                             (interp1(i,j)+interp1(i,j+1)) )      &
     &                           * HM_Cxy2(i,j,k)
              End Do
            End Do

        End If ! k >= first_constant_r_rho_level

! zero u derivative at pole
        If (at_extremity(PSouth)) Then
           Do i = 0, row_length
             interp5(i,1) = 0.
           End Do
        End If

        If (at_extremity(PNorth)) Then
           Do i = 0, row_length
             interp5(i,rows) = 0.
           End Do
        End If

! v increment
! first calculate lambda derivative of exner_prime, store in interp1
! calculate R_u at phi points, store in interp2
! Final calculation of v_prime is after u increment calculation.
        Do j = j_start, j_stop_v + 1
          Do i = i_start, i_stop
              interp2(i,j) = .5 * (R_u(i-1,j,k) *                       &
     &                             recip_1p_cd_lambda(i-1,j,k)          &
     &                           + R_u(i,j,k) *                         &
     &                             recip_1p_cd_lambda(i,j,k))
              interp1(i,j) = .5 * (interp5(i-1,j) *                     &
     &                             recip_1p_cd_lambda(i-1,j,k)          &
     &                           + interp5(i,j) *                       &
     &                             recip_1p_cd_lambda(i,j,k))
            End Do
          End Do

! u increment
! (not updated on LAM boundaries)

        temp1 = alpha_3 * timestep

        Do j = j_begin, j_end
          Do i = i_start, i_stop_u
              R_u(i,j,k) = Au(i,j,k) *                                  &
     &                         (R_u(i,j,k) - temp1 * interp5(i,j))      &
     &                   + Fu(i,j,k) * .25 * (                          &
     &                     R_v(i,j,k) * recip_1p_cd_phi(i,j,k)          &
     &                   + R_v(i+1,j,k) * recip_1p_cd_phi(i+1,j,k)      &
     &                   + R_v(i,j-1,k) * recip_1p_cd_phi(i,j-1,k)      &
     &                   + R_v(i+1,j-1,k)* recip_1p_cd_phi(i+1,j-1,k)   &
     &        - temp1 * (interp6(i,j) * recip_1p_cd_phi(i,j,k)          &
     &                 + interp6(i+1,j) * recip_1p_cd_phi(i+1,j,k)      &
     &                 + interp6(i,j-1) * recip_1p_cd_phi(i,j-1,k)      &
     &                 + interp6(i+1,j-1)* recip_1p_cd_phi(i+1,j-1,k))  &
     &                  )
            End Do
          End Do

! zero u derivative at pole
        If (at_extremity(PSouth)) Then
           Do i = 1, row_length
             R_u(i,1,k) = 0.
           End Do
        End If

        If (at_extremity(PNorth)) Then
           Do i = 1, row_length
             R_u(i,rows,k) = 0.
           End Do
        End If

        Do j = j_start, j_stop_v
          Do i = i_start, i_stop
              R_v(i,j,k) = Av(i,j,k) *                                  &
     &                         (R_v(i,j,k) - temp1 * interp6(i,j) )     &
     &                       - Fv(i,j,k) *                              &
     &                         .5 * (interp2(i,j) + interp2(i,j+1)      &
     &                               - temp1 *                          &
     &                           (interp1(i,j) + interp1(i,j+1) ) )
            End Do
          End Do

      End Do  ! k = 1, model_levels

      DEALLOCATE ( dexprdr )

      IF (lhook) CALL dr_hook('PE_HELMHOLTZ_GLOBAL_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PE_Helmholtz_global_2B
