! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine PE_Helmholtz_eul_mix
      Subroutine PE_Helmholtz_eul_mix(                                  &
     &                       u, v, w, r_theta_levels,                   &
     &                       r_rho_levels, p, rho, rho_np1,             &
     &                       theta, theta_star, theta_np1,              &
     &                       mix_v, mix_cl, mix_cf,                     &
     &                       mix_cf2, mix_rain, mix_graup,              &
     &                       mix_v_star, mix_cl_star, mix_cf_star,      &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                       mix_v_np1, mix_cl_np1, mix_cf_np1,         &
     &                       mix_cf2_np1, mix_rain_np1, mix_graup_np1,  &
     &                       L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,     &
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
     &                       dtheta_dr_term, exner_lbcs,                &
     &                       LENRIM, LBC_SIZE, LBC_START,               &
     &                       RIMWIDTH, RIMWEIGHTS,                      &
     &                       me, n_proc, n_procx, n_procy,              &
     &                       halo_i, halo_j, l_datastart,               &
     &                       L_regular,at_extremity,                    &
     &                       rims_to_do, off_x, off_y,                  &
     &                       proc_row_group, proc_col_group,            &
     &                       global_row_length, global_rows,            &
     &                       g_rows, g_row_length,                      &
     &                       CycleNo, L_new_tdisc, L_fint_theta,        &
     &                       L_lbc_new, L_fixed_lbcs, L_transparent,    &
     &                       ldump)
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

      Use swapable_field_mod, Only: &
          swapable_field_pointer_type

      USE atmos_constants_mod, ONLY:  &
          cp, r, kappa, repsilon, recip_epsilon

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
     &, n_rows                                                          &
                   ! Local number of rows in a v field
     &, off_x,off_y                                                     &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, proc_col_group                                                  &
     &, rims_to_do                                                      &
                           ! zone where lbc rim wts = 1
     &, RIMWIDTH                                                        &
                          ! IN : Size of boundary region
     &, LENRIM(Nfld_max,NHalo_max)                                      &
                          ! IN : Size of single level of LBC
     &, LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
                          ! IN : Size of a side of a LBC
     &, LBC_START(4,Nfld_max,NHalo_max)                                 &
                          ! IN : Start of a side in a LBC
     &, CycleNo                                                        

      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        

      Logical                                                           &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup

      Integer                                                           &
     &  neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                    !  true for regular resolution
     &, L_lbc_new                                                       &
                    !  not active in this subroutine
     &, L_fixed_lbcs                                                    &
                    !  true for fixed lbcs
     &, L_transparent                                                   &
                    !  true for transparent lbcs (idealised runs)
     &, L_new_tdisc                                                     &
     &, L_fint_theta
                     ! true:  fully-interpolating semi-lagrangian
                     !        theta advection will be used
                     ! false: standard non-interpolating in the vertical

      REAL                                                              &
     &  RIMWEIGHTS(RIMWIDTH)  ! IN : weight to apply to LBC

! Only needed as in argument list: Not used
      Real, Intent (In) ::                                              &
     &  exner_lbcs(LENRIM(fld_type_p,halo_type_extended), model_levels)

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
     &  L_Physics  ! true if physics wanted

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
     &, theta_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &              model_levels)                                       &
     &, mix_v  (1-halo_i:row_length+halo_i,                             &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, mix_cf (1-halo_i:row_length+halo_i,                             &
     &     1-halo_j:rows+halo_j, wet_model_levels)                      &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &       wet_model_levels)                                          &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &       wet_model_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &       wet_model_levels)                                          &
     &, mix_v_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, mix_cl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          wet_model_levels)                                       &
     &, mix_cf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          wet_model_levels)                                       &
     &, mix_cf2_star (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_star (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_star (1-off_x:row_length+off_x,                       &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_v_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_model_levels)                                       &
     &, mix_cl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, mix_cf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, mix_cf2_np1 (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_np1 (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_np1 (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)


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
     &  L_poles    !  switch for use in etadot_calc, Global = T

      Integer                                                           &
     &  i, j, k                                                         &
                   ! Loop indices
     &, number_dif_horiz_points                                         &
     &, j0, j1, j1p1, info, count(model_levels)                         &
     &, i_start, i_end, i_end_u                                         &
     &, j_start, j_end, j_one

      Integer :: i_field  ! counter for swappable fields
     
      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi                                                 &
     &, recip_timestep

      Real                                                              &
     &  weight1                                                         &
     &, temp1                                                           &
     &, temp2                                                           &
     &, temp3                                                           &
     &, temp4                                                           &
     &, Ju                                                              &
     &, Jv                                                              &
     &, G_term

       Real                                                             &
              !  VarRes
     &  wk_temp(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)       &
     &, weight2
! Helmholtz equation coefficients
      Real, Target ::                                                   &
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
     &,HM_C3 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C4 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C5 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_RHS(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

! Note HM_C5 redundant storage, dry_rho_theta could be passed into
! GCR_k and save this work-space

! terms in WP154, stored for use in calculating increments after
! solving Helmholtz equation.
      Real                                                              &
     &  Au (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)  &
     &, Av (1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels) &
     &, Fu (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)  &
     &, Fv (1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels) &
     &, recip_G_term (row_length, rows, model_levels-1)                 &
     &, K_term (row_length, rows, model_levels-1)

! 3-d work arrays
      Real                                                              &
     &  thetav_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, rescale (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, weight_upper(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          &
     &, weight_lower(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          &
     &, recip_1p_cd_lambda(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, recip_1p_cd_phi(1-off_x:row_length+off_x,                       &
     &                  1-off_y:rows+off_y, model_levels)

      Real,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     & dry_density                                                      &
     &,dry_rho_theta                                                    &
     &,thetav_star_rho                                                  &
     &,etadot                                                           &
     &,etadot_star                                                      &
     &,etadot_star2

! 2-d work arrays for storing interpolated fields.
      Real                                                              &
     &  interp1(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp2(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp3(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp4(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp5(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, interp6(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, u_mask(1-off_x:row_length+off_x, 1-off_y:rows+off_y)            &
     &, v_mask(1-off_x:row_length+off_x, 1-off_y:rows+off_y)            &
     &, u_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

       Real                                                             &
     &  sum_s(model_levels)                                             &
     &, sum_n(model_levels)                                             &
     &, temp_s(row_length,model_levels)                                 &
     &, temp_n(row_length,model_levels)

      Real :: rhs_cut(row_length,rows,model_levels)
      real                                                              &
     &  convert_to_v
      Real,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     &  convert_to_vstar                                                &
     &, theta_v                                                         &
     &, theta_v_rho                                                     &
     &, dr_deta                                                         &
     &, moist                                                           &
     &, moist_star, moist_np1

      Real, Allocatable, Target ::                                      &
     &  u_star_minus_u(:,:,:)                                           &
     &, v_star_minus_v(:,:,:)

       Type(swapable_field_pointer_type) :: fields_to_swap(15)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines:
      External                                                          &
     &  GCR_k, Etadot_Calc

! ----------------------------------------------------------------------
! Section 1.   Calculate Residual in equation of state.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('PE_HELMHOLTZ_EUL_MIX',zhook_in,zhook_handle)

      ALLOCATE ( dry_density(1-off_x:row_length+off_x,                  &
     &                   1-off_y:rows+off_y, model_levels) )
      ALLOCATE ( thetav_star_rho(1-off_x:row_length+off_x,              &
     &                   1-off_y:rows+off_y, model_levels) )
      ALLOCATE (moist (1-halo_i:row_length+halo_i,                      &
     &     1-halo_j:rows+halo_j, wet_model_levels) )
      ALLOCATE (moist_star (1-off_x:row_length+off_x,                   &
     &        1-off_y:rows+off_y,  wet_model_levels) )
      If ( CycleNo > 1 .and. L_new_tdisc ) Then
        ALLOCATE (moist_np1 (1-off_x:row_length+off_x,                  &
     &          1-off_y:rows+off_y,  wet_model_levels) )
      End If

! Number of points over which solver works
      L_poles = .false.      
      If (model_domain  ==  mt_global) Then
       L_poles = .true.
       number_dif_horiz_points = (global_rows-2) * global_row_length + 2
      Else If (model_domain  ==  mt_lam) Then
       number_dif_horiz_points = (global_rows - 2 * rims_to_do) *       &
     &                           (global_row_length - 2 * rims_to_do)
      Else If (model_domain  ==  mt_cyclic_lam) Then
       number_dif_horiz_points = (global_rows-2) * global_row_length
      Else If (model_domain  ==  mt_bi_cyclic_lam) Then
        number_dif_horiz_points= global_rows * global_row_length
      End If

      moist     =mix_v     +mix_cl     +mix_cf
      moist_star=mix_v_star+mix_cl_star+mix_cf_star

      If(L_mcr_qcf2)then
        moist      = moist + mix_cf2
        moist_star = moist_star + mix_cf2_star
      endif
      If(L_mcr_qrain)then
        moist      = moist + mix_rain
        moist_star = moist_star + mix_rain_star
      endif
      If(L_mcr_qgraup)then
        moist      = moist + mix_graup
        moist_star = moist_star + mix_graup_star
      endif

      j0 = 1
      j1 = rows
      j1p1 = rows + 1
      j_one = 1
      If (model_domain  /=  mt_bi_cyclic_lam) then
      If (at_extremity(PSouth)) Then
        j0 = 2
        j_one = 2
      End If
      If (at_extremity(PNorth)) Then
        j1 = rows-1
        j1p1 = j1
      End If
      Endif
      i_start = 1
      i_end = row_length
      i_end_u = row_length
      j_start = 1
      j_end = n_rows
      If (model_domain  ==  mt_LAM) Then
        If (at_extremity(PWest)) Then
          i_start = 2
        End If
        If (at_extremity(PEast)) Then
          i_end = row_length-1
          i_end_u = row_length-2
        End If
        If (at_extremity(PSouth)) Then
          j_start = 2
          j_one = 1
        End If
        If (at_extremity(PNorth)) Then
          j_end = n_rows-1
          j1p1 = rows
        End If
      End If
!zg
      If (model_domain  ==  mt_cyclic_LAM) Then
        If (at_extremity(PSouth)) Then
          j_start = 2
        End If
        If (at_extremity(PNorth)) Then
          j_end = n_rows-1
        End If
      End If
      recip_delta_lambda = 1. / delta_lambda
      recip_delta_phi = 1. / delta_phi
      recip_timestep = 1. / timestep

! calculate mask which decides whether a v value substituted into
! the u equation is a boundary condition or not. Only active in LAM
! calculate mask which decides whether a u value substituted into
! the v equation is a boundary condition or not. Only active in LAM
! Only 0 in LAM case on boundaries.
      Do j = j0-1, j1+1
        Do i = 0, row_length
          u_mask(i,j) = 1.0
        End Do
      End Do
      Do j = j0-1, j1
        Do i = 0, row_length+1
          v_mask(i,j) = 1.0
        End Do
      End Do

! Zero weights at LAM boundary points
      If (model_domain  ==  mt_LAM) Then
        If (at_extremity(PSouth)) Then
! u
          Do i = 0, row_length
            u_mask(i,1) = 0.0
          End Do
! v
          Do i = 0, row_length+1
            v_mask(i,1) = 0.0
          End Do
        End If
        If (at_extremity(PNorth)) Then
! u
          Do i = 0, row_length
            u_mask(i,rows) = 0.0
          End Do
! v
          Do i = 0, row_length+1
            v_mask(i,n_rows) = 0.0
          End Do
        End If
        If (at_extremity(PWest)) Then
! u
          Do j = j0-1, j1+1
            u_mask(0,j) = 0.0
            u_mask(1,j) = 0.0
          End Do
! v
          Do j = j0-1, j1
            v_mask(0,j) = 0.0
            v_mask(1,j) = 0.0
          End Do
        End If
        If (at_extremity(PEast)) Then
! u
          Do j = j0-1, j1+1
            u_mask(row_length,j) = 0.0
            u_mask(row_length-1,j) = 0.0
          End Do
! v
          Do j = j0-1, j1
            v_mask(row_length,j) = 0.0
            v_mask(row_length+1,j) = 0.0
          End Do
        End If
      End If ! on model_domain = mt_LAM

! ----------------------------------------------------------------------
       ALLOCATE ( convert_to_vstar(1-off_x:row_length+off_x,            &
     &                    1-off_y:rows+off_y, model_levels) )
       ALLOCATE ( theta_v(1-off_x:row_length+off_x,                     &
     &                    1-off_y:rows+off_y, model_levels) )

      weight1 = recip_epsilon
! set up vertical interpolation weights between theta levels and rho
! levels
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,convert_to_v)                &
!$OMP& SHARED(model_levels,j0,j1,row_length,weight_upper,weight_lower,  &
!$OMP&        r_rho_levels,r_theta_levels,wet_model_levels,weight1,     &
!$OMP&        mix_v_star,moist_star,theta,theta_star,convert_to_vstar,  &
!$OMP&        mix_v,moist,theta_v,thetav_star)
!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            weight_upper(i,j,k) = (r_rho_levels(i,j,k)                  &
     &                             - r_theta_levels(i,j,k-1) )          &
     &                           / (r_theta_levels(i,j,k)               &
     &                              - r_theta_levels(i,j,k-1) )
            weight_lower(i,j,k) = 1.0 - weight_upper(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO NOWAIT

! Calculate thetav at theta levs
!$OMP DO SCHEDULE(STATIC)
      Do k = 1, wet_model_levels
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            convert_to_v = ( 1. +  weight1 * mix_v(i,j,k) )/            &
     &                                (1. + moist(i,j,k) )
            convert_to_vstar(i,j,k) = (1. + weight1 * mix_v_star(i,j,k))&
     &                              / (1. + moist_star(i,j,k))
            theta_v(i,j,k) = theta(i,j,k) * convert_to_v
            thetav_star(i,j,k) = theta_star(i,j,k) *                    &
     &                           convert_to_vstar(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO
!$OMP END PARALLEL
      Do k = 1 + wet_model_levels, model_levels
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            theta_v(i,j,k) = theta(i,j,k)
            thetav_star(i,j,k) = theta_star(i,j,k)
            convert_to_vstar(i,j,k) = 1.0
          End Do
        End Do
      End Do

! Interpolate thetav to rho levels, theta to rho levels, and calculate
! Residual in equation of state.
! theta at rho levels is stored in HM_Cyy2
! Residual in equation of state is stored in HM_RHS.

       ALLOCATE ( theta_v_rho(1-off_x:row_length+off_x,                 &
     &                    1-off_y:rows+off_y, model_levels) )

      temp1 = 1. / Cp

      k = 1
      Do j = j0-1, j1+1
        Do i = 0, row_length+1
! thetav at rho levels
          theta_v_rho(i,j,k) = theta_v(i,j,k)
          thetav_star_rho(i,j,k) = thetav_star(i,j,k)
! RHS term.
          HM_RHS(i,j,k) = kappa * rho(i,j,k) * exner_rho_levels(i,j,k)  &
     &                    *  thetav_star_rho(i,j,k)                     &
      
     &                  - r_rho_levels(i,j,k) * r_rho_levels(i,j,k) *   &
     &                    p(i,j,k) * temp1

        End Do
      End Do

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)         &
!$OMP& SHARED(model_levels,j0,j1,row_length,theta_v_rho,weight_upper,   &
!$OMP&        theta_v,weight_lower,thetav_star_rho,thetav_star,HM_RHS,  &
!$OMP&        rho,exner_rho_levels,r_rho_levels,p,temp1)
      Do k = 2, model_levels
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
! thetav at rho levels
              theta_v_rho(i,j,k) = weight_upper(i,j,k) * theta_v(i,j,k) &
     &                     +   weight_lower(i,j,k) * theta_v(i,j,k-1)
              thetav_star_rho(i,j,k) = weight_upper(i,j,k) *            &
     &                                 thetav_star(i,j,k) +             &
     &                                 weight_lower(i,j,k) *            &
     &                                 thetav_star(i,j,k-1)
! RHS term.
              HM_RHS(i,j,k) = kappa * rho(i,j,k)                        &
     &                        * exner_rho_levels(i,j,k)                 &
     &                    *  thetav_star_rho(i,j,k)                     &
     &                    - r_rho_levels(i,j,k) * r_rho_levels(i,j,k) * &
     &                      p(i,j,k) * temp1
            End Do
          End Do
      End Do
!$OMP END PARALLEL DO
      DEALLOCATE ( theta_v)

! Due to changes in u, v discretization theta_np1 is required
! instead of theta_star in H terms Cxy2, Cyx2, Cxx2, Cyy2, Cxp, Cyp
! Obtain corresponding thetav_star.
! Note that thetav_star in HM_RHS, initialised earlier,
! has intentionally not been replaced as it should not change.
      If ( L_new_tdisc .and. CycleNo > 1 ) Then

        moist_np1 = mix_v_np1 + mix_cl_np1 + mix_cf_np1

        If (L_mcr_qcf2) moist_np1 = moist_np1 + mix_cf2_np1
        If (L_mcr_qrain) moist_np1 = moist_np1 + mix_rain_np1
        If (L_mcr_qgraup) moist_np1 = moist_np1 + mix_graup_np1

        k=1
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              thetav_star(i,j,k) = theta_np1(i,j,k) *                   &
     &            (1.+mix_v_np1(i,j,k)/repsilon)/(1.+moist_np1(i,j,k))
              thetav_star_rho(i,j,k) = thetav_star(i,j,k)
            End Do
          End Do

        Do k = 2, wet_model_levels
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              thetav_star(i,j,k) = theta_np1(i,j,k) *                   &
     &            (1.+mix_v_np1(i,j,k)/repsilon)/(1.+moist_np1(i,j,k))
              thetav_star_rho(i,j,k) = weight_upper(i,j,k) *            &
     &                                 thetav_star(i,j,k) +             &
     &                                 weight_lower(i,j,k) *            &
     &                                 thetav_star(i,j,k-1)
            End Do
          End Do
        End Do

        Do k = 1 + wet_model_levels, model_levels
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              thetav_star(i,j,k) = theta_np1(i,j,k)
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
            dtheta_dr_term(i,j,k) =                                     &
     &                (theta_star(i,j,k+1) - theta_star(i,j,k)) /       &
     &                  (r_theta_levels(i,j,k+1) -                      &
     &                   r_theta_levels(i,j,k))
          End Do
        End Do
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)         &
!$OMP& SHARED(model_levels,rows,row_length,dtheta_dr_term,theta_star,   &
!$OMP&        r_theta_levels)
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              dtheta_dr_term(i,j,k) =                                   &
     &                (theta_star(i,j,k+1) - theta_star(i,j,k-1)) /     &
     &                 (r_theta_levels(i,j,k+1) -                       &
     &                  r_theta_levels(i,j,k-1))
            End Do
          End Do
        End Do
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& PRIVATE(i,j,k,temp3,g_term)                                      &
!$OMP& SHARED(model_levels,rows,row_length,temp1,exner_rho_levels,      &
!$OMP&        r_rho_levels, cH, timestep,G_term_tol,temp4,              &
!$OMP&        recip_timestep,recip_G_term,K_term,temp2,thetav_star,     &
!$OMP&        dtheta_dr_term,convert_to_vstar)
      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
! dry_density converts theta to thetav
            temp3 =  - temp1 * (exner_rho_levels(i,j,k+1) -             &
     &                          exner_rho_levels(i,j,k)) /              &
     &                         (r_rho_levels(i,j,k+1) -                 &
     &                          r_rho_levels(i,j,k)) *                  &
     &                          dtheta_dr_term(i,j,k)                   &
     &                    * convert_to_vstar(i,j,k)
            G_term = 1. + temp3 / ( 1. + cH(i,j,k) * timestep )

            If ( G_term  <   G_term_tol ) Then
              cH(i,j,k) = ( temp3 * temp4 - 1.) * recip_timestep
              G_term = G_term_tol
            End If

            recip_G_term(i,j,k) = 1. / G_term
            K_term(i,j,k) = temp2 * thetav_star(i,j,k) *                &
     &                      recip_G_term(i,j,k)
            dtheta_dr_term(i,j,k) = dtheta_dr_term(i,j,k)               &
     &                  / (1. + cH(i,j,k) * timestep)
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO

      Do k = 1, model_levels

! ----------------------------------------------------------------------
! Section 2.2  Form Cd terms.
!              Code assumes that the number of boundary layer levels
!              is at least 2.
! ----------------------------------------------------------------------

        If (k  <=  boundary_layer_levels) Then

          If (L_physics ) Then

! Remove factor of r**2 from rho, store in interp5

            Do j = j0-1, j1+1
              Do i = 0, row_length + 1
                interp5(i,j) = rho(i,j,k) / (r_rho_levels(i,j,k) *      &
     &                                       r_rho_levels(i,j,k) )
              End Do
            End Do

            If (k  ==  1) Then
              Do j = j0-1, j1+1
                Do i = 0, row_length

                  temp1 = rho_Km(i,j,k-1) + rho_Km(i+1,j,k-1)

                  temp2 = rho_Km(i,j,k)/                                &
     &              (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) +     &
     &              rho_Km(i+1,j,k)/                                    &
     &              ( r_rho_levels(i+1,j,k+1) - r_rho_levels(i+1,j,k))

                  interp1(i,j) = alpha_Cd(k) * ( temp1 + temp2 )        &
     &              / (interp5(i,j) * (r_theta_levels(i,j,k) -          &
     &              r_theta_levels(i,j,k-1))                            &
     &              + interp5(i+1,j) * (r_theta_levels(i+1,j,k)         &
     &              - r_theta_levels(i+1,j,k-1)))

                End Do
              End Do

            Else If (k  <   boundary_layer_levels) Then
              Do j = j0-1, j1+1
                Do i = 0, row_length

                  temp1 = rho_Km(i,j,k-1)/                              &
     &              (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)) +     &
     &              rho_Km(i+1,j,k-1)/                                  &
     &              ( r_rho_levels(i+1,j,k) - r_rho_levels(i+1,j,k-1))

                  temp2 = rho_Km(i,j,k)/                                &
     &              (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) +     &
     &              rho_Km(i+1,j,k)/                                    &
     &              ( r_rho_levels(i+1,j,k+1) - r_rho_levels(i+1,j,k))

                  interp1(i,j) = alpha_Cd(k) * ( temp1 + temp2 )        &
     &              / (interp5(i,j) * (r_theta_levels(i,j,k) -          &
     &              r_theta_levels(i,j,k-1))                            &
     &              + interp5(i+1,j) * (r_theta_levels(i+1,j,k)         &
     &              - r_theta_levels(i+1,j,k-1)))

                End Do
              End Do

            Else ! top boundary layer level

              Do j = j0-1, j1+1
                Do i = 0, row_length
                  temp1 = rho_Km(i,j,k-1)/                              &
     &            (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)) +       &
     &            rho_Km(i+1,j,k-1)/                                    &
     &            (r_rho_levels(i+1,j,k) - r_rho_levels(i+1,j,k-1))

                  temp2 = 0.0

                  interp1(i,j) = alpha_Cd(k) * ( temp1 + temp2 )        &
     &            / (interp5(i,j) * (r_theta_levels(i,j,k) -            &
     &            r_theta_levels(i,j,k-1))                              &
     &            + interp5(i+1,j) * (r_theta_levels(i+1,j,k)           &
     &            - r_theta_levels(i+1,j,k-1)))
                End Do
              End Do

            End If

! store Cd at v points in interp2

            If (k  ==  1) Then

              Do j = j0-1, j1
                Do i = 0, row_length+1

                  temp1 = rho_Km(i,j,k-1) + rho_Km(i,j+1,k-1)

                  temp2 = rho_Km(i,j,k)/                                &
     &            (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) +       &
     &            rho_Km(i,j+1,k)/                                      &
     &            ( r_rho_levels(i,j+1,k+1)-r_rho_levels(i,j+1,k))

                  interp2(i,j) = alpha_Cd(k) * ( temp1 + temp2 )        &
     &            / (interp5(i,j) * (r_theta_levels(i,j,k) -            &
     &            r_theta_levels(i,j,k-1)) + interp5(i,j+1) *           &
     &            (r_theta_levels(i,j+1,k)-r_theta_levels(i,j+1,k-1)) )
                End Do
              End Do

            Else If (k  <   boundary_layer_levels) Then

              Do j = j0-1, j1
                Do i = 0, row_length+1

                  temp1 = rho_Km(i,j,k-1)/                              &
     &              (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)) +     &
     &              rho_Km(i,j+1,k-1)/                                  &
     &              ( r_rho_levels(i,j+1,k)-r_rho_levels(i,j+1,k-1))

                  temp2 = rho_Km(i,j,k)/                                &
     &            (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) +       &
     &            rho_Km(i,j+1,k)/                                      &
     &            ( r_rho_levels(i,j+1,k+1)-r_rho_levels(i,j+1,k))

                  interp2(i,j) = alpha_Cd(k) * ( temp1 + temp2 )        &
     &            / (interp5(i,j) * (r_theta_levels(i,j,k) -            &
     &            r_theta_levels(i,j,k-1)) + interp5(i,j+1) *           &
     &            (r_theta_levels(i,j+1,k)-r_theta_levels(i,j+1,k-1)) )
                End Do
              End Do

            Else

              Do j = j0-1, j1
                Do i = 0, row_length+1
                  temp1 = rho_Km(i,j,k-1)/                              &
     &            (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)) +       &
     &            rho_Km(i,j+1,k-1)/                                    &
     &            (r_rho_levels(i,j+1,k)-r_rho_levels(i,j+1,k-1))

                  temp2 = 0.0

                  interp2(i,j) = alpha_Cd (k)* ( temp1 + temp2 )        &
     &            / (interp5(i,j) * (r_theta_levels(i,j,k) -            &
     &            r_theta_levels(i,j,k-1)) + interp5(i,j+1) *           &
     &            (r_theta_levels(i,j+1,k)-r_theta_levels(i,j+1,k-1)) )
                End Do

              End Do
            End If

          Else     ! If NOT L_Physics

! store Cd at u points in interp1

            Do j = j0-1, j1+1
              Do i = 0, row_length
                interp1(i,j) = frictional_timescale(k)
              End Do
            End Do

! store Cd at v points in interp2

            Do j = j0-1, j1
              Do i = 0, row_length+1
                interp2(i,j) = frictional_timescale(k)
              End Do
            End Do

          End If

! New for 2.6
! Zero cd_at_u, cd_at_v at LAM boundary points
          If (model_domain  ==  mt_LAM) Then
            If (at_extremity(PSouth)) Then
! u
              Do i = 0, row_length
                interp1(i,1) = 0.0
              End Do
! v
              Do i = 0, row_length+1
                interp2(i,1) = 0.0
              End Do
            End If
            If (at_extremity(PNorth)) Then
! u
              Do i = 0, row_length
                interp1(i,rows) = 0.0
              End Do
! v
              Do i = 0, row_length+1
                interp2(i,n_rows) = 0.0
              End Do
            End If
            If (at_extremity(PWest)) Then
! u
              Do j = j0-1, j1+1
                interp1(0,j) = 0.0
                interp1(1,j) = 0.0
              End Do
! v
              Do j = j0-1, j1
                interp2(0,j) = 0.0
                interp2(1,j) = 0.0
              End Do
            End If
            If (at_extremity(PEast)) Then
! u
              Do j = j0-1, j1+1
                interp1(row_length,j) = 0.0
                interp1(row_length-1,j) = 0.0
              End Do
! v
              Do j = j0-1, j1
                interp2(row_length,j) = 0.0
                interp2(row_length+1,j) = 0.0
              End Do
            End If
          End If ! on model_domain = mt_LAM

          Do j = j0-1, j1+1
            Do i = 0, row_length
              recip_1p_cd_lambda(i,j,k) = 1.0                           &
     &                              / (1.0 + interp1(i,j) * timestep)
            End Do
          End Do
          Do j = j0-1, j1
            Do i = 0, row_length+1
              recip_1p_cd_phi(i,j,k) = 1.0                              &
     &                              / (1.0 + interp2(i,j) * timestep)
            End Do
          End Do
        Else ! above BL
          Do j = j0-1, j1+1
            Do i = 0, row_length
              recip_1p_cd_lambda(i,j,k) = 1.0
            End Do
          End Do
          Do j = j0-1, j1
            Do i = 0, row_length+1
              recip_1p_cd_phi(i,j,k) = 1.0
            End Do
          End Do
        End If ! on inside BL

! ----------------------------------------------------------------------
! Section 2.3  Form Ju and hence Au and Fu terms.
! ----------------------------------------------------------------------

! changed for 2.6

! calculate Ju and hence Au and Fu
! temp1 holds averaged term

        temp3 = alpha_3 * alpha_3 * timestep * timestep

        Do j = j0, j1
          Do i = 0, row_length
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

        Do j = j0-1, j1
          Do i = 1, row_length
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

! End loop over model_levels
      End Do


! ----------------------------------------------------------------------
! Section 3    Set up Helmholtz equation coefficients and
!              right-hand-side.
! ----------------------------------------------------------------------

      ALLOCATE ( dry_rho_theta(row_length, rows, model_levels) )
      ALLOCATE ( etadot(1-off_x:row_length+off_x,                       &
     &                   1-off_y:rows+off_y, 0:model_levels) )
      ALLOCATE ( etadot_star(1-off_x:row_length+off_x,                  &
     &                   1-off_y:rows+off_y, 0:model_levels) )
      ALLOCATE ( etadot_star2(1-off_x:row_length+off_x,                 &
     &                   1-off_y:rows+off_y, 0:model_levels) )
      ALLOCATE ( v_star_minus_v(1-off_x:row_length+off_x,               &
     &                   1-off_y:n_rows+off_y, model_levels) )
      ALLOCATE ( dr_deta(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, 0:model_levels) )
! ----------------------------------------------------------------------
! Section 3.1  Re-scale Helmholtz equation right-hand-side and
!              set C3 and C4 coefficients.
! ----------------------------------------------------------------------
! at the moment HM_RHS holds residual in equation of state,
! Store dry rho in dry_density

! calculate d(r)/d(eta) about rho levels
! In CycleNo > 1 dry_density needs to be defined in terms of rho^(1)
! for the improved rho-discretization to take place. Rho^(1) is an
! alpha1 weighted t-average between rho and rho_np1. dry_density
! is used in Helmholtz coefficients Cxx1, Cyy1, Czz, C5.

      Do k = 1, model_levels
        Do j = j0-off_y , j1+off_y
          Do i = 0, row_length+1
            dr_deta(i,j,k) = (r_theta_levels(i,j,k) -                   &
     &                      r_theta_levels(i,j,k-1)) /                  &
     &                     (eta_theta_levels(k) -                       &
     &                      eta_theta_levels(k-1))
          End Do
        End Do
      End Do


        Do k = 1, model_levels
          Do j = j0-off_y , j1+off_y
            Do i = 0, row_length+1
              rescale(i,j,k) = recip_timestep * dr_deta(i,j,k)          &
     &                     * rho(i,j,k)
            End Do
          End Do
        End Do

      If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then
        k = 1
        Do j = j0-off_y , j1+off_y
          Do i = 0, row_length+1
            dry_density(i,j,k) = rho(i,j,k) / (1. + moist(i,j,k))
          End Do
        End Do
      Else
        k = 1
        Do j = j0-off_y , j1+off_y
          Do i = 0, row_length+1
            dry_density(i,j,k)=(1.-alpha_1)*rho(i,j,k)/(1.+moist(i,j,k))&
     &                    + alpha_1*rho_np1(i,j,k)/(1.+moist_np1(i,j,k))
          End Do
        End Do
      End If

        Do j = 1, rows
          Do i = 1, row_length
            temp2 = 1.+ moist_star(i,j,k)
            temp3 = dr_deta(i,j,k)                                      &
     &              /(theta_v_rho(i,j,k) * timestep * kappa *           &
     &                exner_rho_levels(i,j,k) * temp2)
            HM_RHS(i,j,k) = - HM_RHS(i,j,k) * temp3
            HM_C4(i,j,k) = ( r_rho_levels(i,j,k) *                      &
     &                       r_rho_levels(i,j,k) * p(i,j,k) /           &
     &                       (R * exner_rho_levels(i,j,k)) -            &
     &                       kappa * theta_v_rho(i,j,k) * rho(i,j,k) )  &
     &                     * temp3
            HM_C3(i,j,k) = dr_deta(i,j,k) * rho(i,j,k)                  &
     &                     / (temp2 * theta_v_rho(i,j,k))
          End Do
        End Do

      Do k = 2, model_levels
        If (k  <=  wet_model_levels) Then
          If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = weight_upper(i,j,k) *                    &
     &                        (1. + moist_star(i,j,k) ) +               &
     &                         weight_lower(i,j,k) *                    &
     &                        (1. + moist_star(i,j,k-1) )
                dry_density(i,j,k) = rho(i,j,k) /                       &
     &                         ( weight_upper(i,j,k)                    &
     &                           * (1. + moist(i,j,k)) +                &
     &                           weight_lower(i,j,k) *                  &
     &                         (1. + moist(i,j,k-1)) )
              End Do
            End Do
          Else
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = weight_upper(i,j,k) *                    &
     &                        (1. + moist_star(i,j,k) ) +               &
     &                         weight_lower(i,j,k) *                    &
     &                        (1. + moist_star(i,j,k-1) )
                dry_density(i,j,k) = (1.-alpha_1)*rho(i,j,k) /          &
     &                 ( weight_upper(i,j,k)*(1.+moist(i,j,k))          &
     &                  + weight_lower(i,j,k)*(1.+moist(i,j,k-1)) )     &
     &               + alpha_1 * rho_np1(i,j,k) /                       &
     &                 ( weight_upper(i,j,k)*(1.+moist_np1(i,j,k))      &
     &                  + weight_lower(i,j,k)*(1.+moist_np1(i,j,k-1)) )
              End Do
            End Do
          End If
        Else If (k  ==  wet_model_levels+1) Then
          If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = ( weight_upper(i,j,k)  +                 &
     &                           weight_lower(i,j,k) *                  &
     &                         (1. + moist_star(i,j,k-1) ) )
                dry_density(i,j,k) = rho(i,j,k) /                       &
     &                       ( weight_upper(i,j,k) +                    &
     &                         weight_lower(i,j,k) *                    &
     &                        (1. + moist(i,j,k-1) ) )
              End Do
            End Do
          Else
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = ( weight_upper(i,j,k)  +                 &
     &                           weight_lower(i,j,k) *                  &
     &                         (1. + moist_star(i,j,k-1) ) )
                dry_density(i,j,k) = (1.-alpha_1)*rho(i,j,k) /          &
     &                  ( weight_upper(i,j,k) +                         &
     &                     weight_lower(i,j,k)*(1.+ moist(i,j,k-1)) )   &
     &                + alpha_1*rho_np1(i,j,k) /                        &
     &                  ( weight_upper(i,j,k) +                         &
     &                     weight_lower(i,j,k)*(1.-moist_np1(i,j,k-1)) )
              End Do
            End Do
          End If
        Else
          If ( CycleNo == 1 .OR. ( .NOT. L_new_tdisc ) ) Then
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = 1.
                dry_density(i,j,k) = rho(i,j,k)
              End Do
            End Do
          Else
            Do j = j0-1, j1+1
              Do i = 0, row_length+1
                interp2(i,j) = 1.
                dry_density(i,j,k) = (1.-alpha_1)*rho(i,j,k)            &
     &                             +  alpha_1*rho_np1(i,j,k)
              End Do
            End Do
          End If
        End If

        Do j = 1, rows
          Do i = 1, row_length
            temp3 = dr_deta(i,j,k)                                      &
     &              /(theta_v_rho(i,j,k) * timestep * kappa *           &
     &                exner_rho_levels(i,j,k) * interp2(i,j))
            HM_RHS(i,j,k) = - HM_RHS(i,j,k) * temp3
            HM_C4(i,j,k) = ( r_rho_levels(i,j,k) *                      &
     &                       r_rho_levels(i,j,k) * p(i,j,k) /           &
     &                       (R * exner_rho_levels(i,j,k)) -            &
     &                       kappa * theta_v_rho(i,j,k) * rho(i,j,k) )  &
     &                     * temp3
            HM_C3(i,j,k) = dr_deta(i,j,k) * rho(i,j,k)                  &
     &                     / (interp2(i,j) * theta_v_rho(i,j,k))
          End Do
        End Do
      End Do

      DEALLOCATE (theta_v_rho)

      do k=1,model_levels
! (m_star-m)/(1+m_star)(1+m)  with m terms averaged to rho levels
        If (k == 1)then
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              interp3(i,j) =   (moist_star(i,j,k)-moist(i,j,k)) /       &
     &         (  (1 + moist_star(i,j,k))*(1+moist(i,j,k))  )
            End Do
          End Do

        elseIf ( k  <=  wet_model_levels) Then
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              interp3(i,j) =                                            &
     &  ( weight_upper(i,j,k) *(moist_star(i,j,k)-moist(i,j,k)) +       &
     &    weight_lower(i,j,k) *(moist_star(i,j,k-1)-moist(i,j,k-1)) ) / &

     &  ( (1. +  (weight_upper(i,j,k)*moist_star(i,j,k) +               &
     &            weight_lower(i,j,k)*moist_star(i,j,k-1) ) ) *         &
     &    (1. +  (weight_upper(i,j,k)*moist(i,j,k) +                    &
     &            weight_lower(i,j,k)*moist(i,j,k-1) ) )  )


            End Do
          End Do
        Else If ( k  ==  wet_model_levels+1) Then
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              interp3(i,j) =                                            &
     &  (  weight_lower(i,j,k) *(moist_star(i,j,k-1)-moist(i,j,k-1)) ) /&

     &  ( (1. +  ( weight_lower(i,j,k)*moist_star(i,j,k-1) ) ) *        &
     &    (1. +  ( weight_lower(i,j,k)*moist(i,j,k-1)      ) )  )
            End Do
          End Do
        Else
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              interp3(i,j) = 0.0
            End Do
          End Do
        End If

          Do j = j0-1, j1+1
            Do i = 0, row_length+1
        HM_RHS(i,j,k) = HM_RHS(i,j,k)  -                                &
     &    dr_deta(i,j,k) *rho(i,j,k)*recip_timestep*interp3(i,j)
            End Do
          End Do

      end do
      DEALLOCATE (moist)
      DEALLOCATE (moist_star)
      If ( CycleNo > 1 .and. L_new_tdisc) DEALLOCATE (moist_np1)

! Calculate dry rho at theta levels
      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
! ------------------------------------------------------
! switch weights for lin interp
! GTG 200799
! ------------------------------------------------------
            dry_rho_theta(i,j,k) = ( dry_density(i,j,k) *               &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      dry_density(i,j,k+1) *                      &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) )
          End Do
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 3.2  Lambda direction coefficients (Cxx1,Cxx2,etc).
! ----------------------------------------------------------------------

! Initialise u_star to zero at poles
      If (model_domain  /=  mt_bi_cyclic_lam) then
      If(at_extremity(PSouth))then
        Do i = 1, row_length
          u_star(i,1) = 0.0
        End Do
      End If
      If(at_extremity(PNorth))then
        Do i = 1, row_length
          u_star(i,rows) = 0.0
        End Do
      End If
      Endif

      ALLOCATE (u_star_minus_u(1-off_x:row_length+off_x,                &
     &                    1-off_y:rows+off_y, model_levels) )

      Do k = 1, model_levels

! Calculate terms requiring averaging to u points.
! average r and store in interp2 as 1/r
! average dry rho * d(r)/d(eta) and store in interp3
! average R_v and store in interp4

        temp1 = recip_delta_lambda * 2.
        Do j = j0, j1
          Do i = 0, row_length
            interp2(i,j) = temp1 / (r_rho_levels(i+1,j,k) +             &
     &                              r_rho_levels(i,j,k))
            interp3(i,j) = (dry_density(i+1,j,k) * dr_deta(i+1,j,k) +   &
     &                      dry_density(i,j,k) * dr_deta(i,j,k)) * .5
            HM_Cxx1(i,j,k) = interp3(i,j) * interp2(i,j)
          End Do
        End Do

! calculate HM_Cxx1, HM_Cxx2 and HM_Cx terms.

        temp1 = alpha_1 * alpha_3 * timestep
        temp2 = Cp * 0.5
        Do j = j0, j1
          Do i = 1, row_length
            HM_Cyx2(i,j,k) = temp2 * (thetav_star_rho(i,j,k) +          &
     &                                thetav_star_rho(i+1,j,k) ) *      &
     &                       interp2(i,j)                               &
     &                       * FV_sec_theta_latitude(i,j)
            HM_Cxx2(i,j,k) = temp1 * Au(i,j,k) * HM_Cyx2(i,j,k)
          End Do
        End Do

        Do j = j0, j1
          Do i = 1, row_length
            HM_Cxy1(i,j,k) = temp1 * Fu(i,j,k)
          End Do
        End Do

! Form u_star


        Do j = j0 ,j1
          Do i = 0, row_length
            u_star(i,j) =u(i,j,k)+ alpha_1 *( Au(i,j,k) * R_u(i,j,k)    &
     &                   + Fu(i,j,k) * .25 *                            &
     &                    (R_v(i+1,j,k) * recip_1p_cd_phi(i+1,j,k)      &
     &                   + R_v(i,j,k) * recip_1p_cd_phi(i,j,k)          &
     &                   + R_v(i+1,j-1,k) * recip_1p_cd_phi(i+1,j-1,k)  &
     &                   + R_v(i,j-1,k) * recip_1p_cd_phi(i,j-1,k) ))
          End Do
        End Do
! set u + alpha_1 u' correctly since u' is known and held in R_u
        If (model_domain  ==  mt_LAM) Then
          If(at_extremity(PNorth)) Then
            Do i = 1, row_length
              u_star(i,rows) = u(i,rows,k)+ alpha_1 *  R_u(i,rows,k)
            End Do
          End If
          If(at_extremity(PSouth)) Then
            Do i = 1, row_length
              u_star(i,1) = u(i,1,k)+ alpha_1 *  R_u(i,1,k)
            End Do
          End If
          If(at_extremity(PWest)) Then
            Do j = j0 ,j1
              u_star(1,j) = u(1,j,k)+ alpha_1 *  R_u(1,j,k)
            End Do
          End If
          If(at_extremity(PEast)) Then
            Do j = j0 ,j1
              u_star(row_length,j) = u(row_length,j,k)                  &
     &                               + alpha_1 * R_u(row_length,j,k)
              u_star(row_length-1,j) = u(row_length-1,j,k)              &
     &                               + alpha_1 * R_u(row_length-1,j,k)
            End Do
          End If
        End If

! Modify right-hand-side with respect to the terms similar to the
! HM_Cxx2 and HM_Cx coefficients

        Do j = j0, j1
          Do i = i_start, i_end
            HM_RHS(i,j,k) = HM_RHS(i,j,k) + sec_theta_latitude(i,j) *   &
     &                        (interp3(i,j) * u_star(i,j)               &
     &                         * interp2(i,j) -                         &
     &                         interp3(i-1,j) * u_star(i-1,j)           &
     &                         * interp2(i-1,j))
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 3.3  Phi direction coefficients (Cyy1,etc) and add terms to
!              right hand side.
! ----------------------------------------------------------------------

! Calculate terms requiring averaging to v points.
! average r and store in interp2 as 1/r
! average dry rho * d(r)/d(eta) and store in interp3
! average R_u in both lamda and phi and store in interp4

        temp1 = recip_delta_phi * 2
        Do j = j0-1, j1
          Do i = 1, row_length
            interp2(i,j) = temp1 / (r_rho_levels(i,j+1,k) +             &
     &                              r_rho_levels(i,j,k))
            interp3(i,j) = (dry_density(i,j+1,k) * dr_deta(i,j+1,k) +   &
     &                      dry_density(i,j,k) * dr_deta(i,j,k)) * .5
            HM_Cyy1(i,j,k) = interp3(i,j) * interp2(i,j)                &
     &                       * cos_v_latitude(i,j)
          End Do
        End Do

! calculate HM_Cyy1, HM_Cyy2 and HM_Cyx1, HM_Cxy2 terms.

        temp1 = alpha_1 * alpha_3 * timestep
        temp2 = Cp * 0.5
        Do j = 1,j1
          Do i = 1, row_length
            HM_Cxy2(i,j,k) = temp2 * (thetav_star_rho(i,j+1,k) +        &
     &                                thetav_star_rho(i,j,k) ) *        &
     &                       interp2(i,j)
            HM_Cyy2(i,j,k) = temp1 * Av(i,j,k) * HM_Cxy2(i,j,k)
          End Do
        End Do

        Do j = 1, j1
          Do i = 1, row_length
            HM_Cyx1(i,j,k) = temp1 * Fv(i,j,k)
          End Do
        End Do

! Form v_star
! store in HM_Cz


! subtract v_star (stored in v_star_minus_v)

        Do j = j0-1 ,j1
          Do i = 1, row_length
            HM_Cz(i,j,k)= v(i,j,k)+ alpha_1 * (Av(i,j,k) * R_v(i,j,k)   &
     &                   - Fv(i,j,k) * 0.25 *                           &
     &                (R_u(i-1,j+1,k) * recip_1p_cd_lambda(i-1,j+1,k)   &
     &               + R_u(i,j+1,k) * recip_1p_cd_lambda(i,j+1,k)       &
     &               + R_u(i-1,j,k) * recip_1p_cd_lambda(i-1,j,k)       &
     &               + R_u(i,j,k)* recip_1p_cd_lambda(i,j,k)) )
          End Do
        End Do

! set v + alpha_1 v' correctly since v' is known and held in R_v
        If (model_domain  ==  mt_LAM) Then
          If(at_extremity(PNorth)) Then
            Do i = 1, row_length
              HM_Cz(i,n_rows,k) = v(i,n_rows,k)                         &
     &                          + alpha_1 * R_v(i,n_rows,k)
            End Do
          End If
          If(at_extremity(PSouth)) Then
            Do i = 1, row_length
              HM_Cz(i,1,k) = v(i,1,k) + alpha_1 * R_v(i,1,k)
            End Do
          End If
          If(at_extremity(PWest)) Then
            Do j = j0 ,j1
              HM_Cz(1,j,k) = v(1,j,k)+ alpha_1 *  R_v(1,j,k)
            End Do
          End If
          If(at_extremity(PEast)) Then
            Do j = j0 ,j1
              HM_Cz(row_length,j,k) = v(row_length,j,k)                 &
     &                               + alpha_1 * R_v(row_length,j,k)
            End Do
          End If
        End If

! Modify right-hand-side with respect to the terms similar to the
! HM_Cyy2 coefficient

        Do j = j0, j1
          Do i = i_start, i_end
            HM_RHS(i,j,k) = HM_RHS(i,j,k)                               &
     &                    +  FV_sec_theta_latitude(i,j) *               &
     &                      ( HM_Cz(i,j,k) * interp3(i,j) *             &
     &                        interp2(i,j) * cos_v_latitude(i,j) -      &
     &                        HM_Cz(i,j-1,k) * interp3(i,j-1) *         &
     &                        interp2(i,j-1) * cos_v_latitude(i,j-1) )
          End Do
        End Do

        If (model_domain  ==  mt_Global) Then
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

        End If

! subtract u from u_star  (store in u_star_minus_u)
        If (model_domain  ==  mt_Global) Then
          Do j = j0 ,j1
            Do i = 1, row_length
              u_star_minus_u(i,j,k) = u_star(i,j) - u(i,j,k)
            End Do
          End Do
        Else If (model_domain  ==  mt_LAM .or.                          &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          Do j = 1, rows
            Do i = 1, row_length
              u_star_minus_u(i,j,k) = u_star(i,j) - u(i,j,k)
            End Do
          End Do
        End If

! subtract v from v_star (stored in v_star_minus_v)
        Do j = 1 ,n_rows
          Do i = 1, row_length
            v_star_minus_v(i,j,k) = HM_Cz(i,j,k) - v(i,j,k)
          End Do
        End Do


! end loop over levels
      End Do

      DEALLOCATE ( dr_deta)

      If (model_domain  ==  mt_Global) Then
! sum terms at poles and add on

        If( at_extremity(PSouth)) Then

          CALL global_2d_sums(temp_s, row_length, 1, 0, 0, model_levels, &
                              sum_s, proc_row_group)

          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) *          &
     &                 cos_v_latitude(1,1) / global_row_length

            Do i = 1, row_length
              HM_RHS(i,1,k) = HM_RHS(i,1,k) + sum_s(k)
            End Do
          End Do
        End If
        If( at_extremity(PNorth)) Then

          CALL global_2d_sums(temp_n, row_length, 1, 0, 0, model_levels, &
                              sum_n, proc_row_group)

          Do k = 1, model_levels
            sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) *       &
     &                 cos_v_latitude(1,n_rows) / global_row_length

            Do i = 1, row_length
              HM_RHS(i,rows,k) = HM_RHS(i,rows,k) - sum_n(k)
            End Do
          End Do
        End If

      End If

! ----------------------------------------------------------------------
! Section 3.4  r direction coefficients (C5,Cz,Czz,Cxz,Cyz) and
!              add terms to right hand side.
! ----------------------------------------------------------------------

! Calculate etadot_star, the known part of the increment to etadot.
! u_star_minus_u currently holds u_star-u
! v_star_minus_v currently holds v_star-v
! store w term in etadot

        Do k = 1, model_levels-1
          Do j = 1, rows
            Do i = 1, row_length
              etadot_star2(i,j,k) = alpha_2 * R_w(i,j,k) *              &
     &                            recip_G_term(i,j,k)
            End Do
          End Do
        End Do


      i_field = 0

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => u_star_minus_u(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_u
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .TRUE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => v_star_minus_v(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_v
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  n_rows
      fields_to_swap(i_field) % vector       =  .TRUE.

! DEPENDS ON: swap_bounds_mv
      CALL SWAP_BOUNDS_MV(fields_to_swap, i_field,                      &
                          row_length, off_x, off_y)

! DEPENDS ON: etadot_calc
      Call Etadot_Calc (r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u_star_minus_u, v_star_minus_v, etadot_star2,   &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  lambda_p, phi_p, lambda_u, phi_v,               &
     &                  wt_lambda_u, wt_phi_v,                          &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  proc_row_group, at_extremity, global_row_length,&
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  off_x, off_y, off_x, off_y, off_x, off_y,       &
     &                  1, row_length, 1, rows, j0, j1,                 &
     &                  L_regular, .true.,                              &
     &                  temp_s, temp_n, etadot_star)

       DEALLOCATE (u_star_minus_u)
! Calculate etadot.

! DEPENDS ON: etadot_calc
      Call Etadot_Calc (r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u, v, w,                                        &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  lambda_p, phi_p, lambda_u, phi_v,               &
     &                  wt_lambda_u, wt_phi_v,                          &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  proc_row_group, at_extremity, global_row_length,&
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  off_x, off_y, off_x, off_y, off_x, off_y,       &
     &                  1, row_length, 1, rows, j0, j1,                 &
     &                  L_regular, L_poles,                             &
     &                  temp_s, temp_n, etadot)


! v_star_minus_v currently holds v_star-v

      Do k = 1, model_levels

! calculate terms required for vertical averaging and differencing
! store the terms for differencing in interp1&2 where interp1 is the
! value at the level above and interp2 the value at the level below.
! store the terms for averaging in interp3&4 where interp3 is the
! value at the level above and interp4 the value at the level below.

! set the coefficients HM_Czz and HM_Cz

        If (k  ==  1) Then
          Do j = 1, rows
            Do i = 1, row_length
              interp2(i,j) = 0.
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))

              interp5(i,j) = alpha_2 * temp2
              interp1(i,j) = dry_rho_theta(i,j,k) *                     &
     &                       ( etadot(i,j,k) + etadot_star(i,j,k) )     &
     &                       / temp2

            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              interp4(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)&
     &                        *convert_to_vstar(i,j,k)
              interp3(i,j) = interp4(i,j)

              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = dry_rho_theta(i,j,k) * temp1

              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)              &
     &                       * convert_to_vstar(i,j,k)

            End Do
          End Do

        Else If (k  ==  model_levels) Then

          Do j = 1, rows
            Do i = 1, row_length

              interp2(i,j) = interp1(i,j)
              interp1(i,j) = 0.

              interp4(i,j) = interp3(i,j)
              interp3(i,j) = 0.

              HM_Cz(i,j,k) = 0.
              HM_Czz(i,j,k) = 0.
            End Do
          End Do


        Else If (k  ==  model_levels - 1) Then

          Do j = 1, rows
            Do i = 1, row_length

! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))

              interp2(i,j) = interp1(i,j)
              interp1(i,j) = dry_rho_theta(i,j,k) *                     &
     &                       ( etadot(i,j,k) +  etadot_star(i,j,k) )    &
     &                       / temp2

              interp5(i,j) = alpha_2 * temp2

            End Do
          End Do

          Do j = 1, rows
            Do i = 1, row_length
              interp4(i,j) = interp3(i,j)
              interp3(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)&
     &                        *convert_to_vstar(i,j,k)

              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = dry_rho_theta(i,j,k) * temp1

              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)              &
     &                    *  convert_to_vstar(i,j,k)

            End Do
          End Do

        Else

          Do j = 1, rows
            Do i = 1, row_length

! calculate d(eta)/d(r) at theta levels
              temp2 = (eta_rho_levels(k+1) - eta_rho_levels(k)) /       &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))

              interp2(i,j) = interp1(i,j)
              interp1(i,j) = dry_rho_theta(i,j,k) *                     &
     &                       ( etadot(i,j,k) +  etadot_star(i,j,k) )    &
     &                       / temp2

              interp5(i,j) = alpha_2 * temp2
            End Do
          End Do

          Do j = 1, rows
            Do i = 1, row_length

              interp4(i,j) = interp3(i,j)
              interp3(i,j) = etadot_star2(i,j,k) * dtheta_dr_term(i,j,k)&
     &                        *convert_to_vstar(i,j,k)

              temp1 = interp5(i,j) * K_term(i,j,k)
              HM_Czz(i,j,k) = dry_rho_theta(i,j,k) * temp1

              HM_Cz(i,j,k) = temp1 * dtheta_dr_term(i,j,k)              &
     &                      *convert_to_vstar(i,j,k)

            End Do
          End Do

        End If

        If ( k  <   first_constant_r_rho_level) Then
          Do j = 1, rows
            Do i = 1, row_length
              HM_C5(i,j,k) = dry_rho_theta(i,j,k)
            End Do
          End Do
        End If

! Modify right-hand-side

        temp1 = 1. / (eta_theta_levels(k) -                             &
     &                eta_theta_levels(k-1))
        Do j = 1, rows
          Do i = 1, row_length
            HM_RHS(i,j,k) = HM_RHS(i,j,k) +                             &
     &                      (interp1(i,j) - interp2(i,j)) * temp1       &
     &                    + HM_C3(i,j,k) *                              &
     &                      (weight_upper(i,j,k) * interp3(i,j) +       &
     &                       weight_lower(i,j,k) * interp4(i,j) )
          End Do
        End Do

! end loop over model_levels
      End Do

      temp1 = recip_delta_lambda * 2.
      temp2 = recip_delta_phi * 2.0
      Do k = 1, first_constant_r_rho_level - 1
! Calculate HM_Cxz and HM_Cyz
! Calculate dr/d lambda
        Do j = j0, j1
          Do i = 1, row_length
            HM_Cxp(i,j,k) =  ( r_rho_levels(i+1,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) * 2.0 /            &
     &                       (thetav_star_rho(i+1,j,k) +                &
     &                        thetav_star_rho(i,j,k) )
            HM_Cxz(i,j,k) =  ( r_theta_levels(i+1,j,k) -                &
     &                         r_theta_levels(i,j,k) )                  &
     &                          * temp1                                 &
     &                          * sec_theta_latitude(i,j)               &
     &                          / ( r_theta_levels(i+1,j,k) +           &
     &                              r_theta_levels(i,j,k) )
          End Do
        End Do

! Calculate dr/d phi
        Do j = 1, n_rows
          Do i = 1, row_length
            HM_Cyp(i,j,k) = (r_rho_levels(i,j+1,k) -                    &
     &                       r_rho_levels(i,j,k) ) * 2.0 /              &
     &                      (thetav_star_rho(i,j+1,k) +                 &
     &                       thetav_star_rho(i,j,k) )
            HM_Cyz(i,j,k) = (r_theta_levels(i,j+1,k) -                  &
     &                         r_theta_levels(i,j,k) )                  &
     &                         * temp2 /                                &
     &                        (r_theta_levels(i,j+1,k) +                &
     &                         r_theta_levels(i,j,k) )
          End Do
        End Do

      End Do

        DEALLOCATE ( dry_density )
        DEALLOCATE ( dry_rho_theta )
        DEALLOCATE ( thetav_star_rho )
        DEALLOCATE ( etadot )
        DEALLOCATE ( etadot_star )
        DEALLOCATE ( etadot_star2 )
        DEALLOCATE ( v_star_minus_v )
        DEALLOCATE ( convert_to_vstar )

! ----------------------------------------------------------------------
! Section 4.   Call elliptic solver to solve Helmholtz equation.
! ----------------------------------------------------------------------

! re-scale HM_cyx2 and HM_cxy2 to incorporate friction in solver
      Do k = 1, boundary_layer_levels
        Do j = j0, j1
          Do i = 1, row_length
            HM_Cyx2(i,j,k) = HM_Cyx2(i,j,k)                             &
     &                       * recip_1p_cd_lambda(i,j,k)
          End Do
        End Do
        Do j = 1,j1
          Do i = 1, row_length
            HM_Cxy2(i,j,k) = HM_Cxy2(i,j,k)                             &
     &                       * recip_1p_cd_phi(i,j,k)
          End Do
        End Do
      End Do

      If (model_domain  /=  mt_bi_cyclic_lam) then
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
      Endif
      If (model_domain  ==  mt_Global) Then
! initialise parts of coefficients not set so far to zero.

        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,1,k) = 0.
              HM_Cyx2(i,1,k) = 0.
              HM_Cxy1(i,1,k) = 0.
              HM_Cxx2(i,1,k) = 0.
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,1,k) = 0.0
              HM_Cxp(i,1,k) = 0.0
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,rows,k) = 0.
              HM_Cxx2(i,rows,k) = 0.
              HM_Cyy1(i,rows,k) = 0.
              HM_Cyy2(i,rows,k) = 0.
              HM_Cxy1(i,rows,k) = 0.
              HM_Cxy2(i,rows,k) = 0.
              HM_Cyx1(i,rows,k) = 0.
              HM_Cyx2(i,rows,k) = 0.
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,rows,k) = 0.0
              HM_Cxp(i,rows,k) = 0.0
              HM_Cyz(i,rows,k) = 0.0
              HM_Cyp(i,rows,k) = 0.0
            End Do
          End Do
        End If

      Else If (model_domain  ==  mt_LAM) Then
! Set coefficients and residuals on boundaries
        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,1,k) = 0.0
              HM_Cxx2(i,1,k) = 0.0
              HM_Cxy1(i,1,k) = 0.0
              HM_Cyx2(i,1,k) = 0.0
              HM_Cyy1(i,1,k) = 0.0
              HM_Cyy2(i,1,k) = 0.0
              HM_Cxy2(i,1,k) = 0.0
              HM_Cyx1(i,1,k) = 0.0
              HM_RHS(i,1,k) = 0.0
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,1,k) = 0.0
              HM_Cxp(i,1,k) = 0.0
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,rows,k) = 0.0
              HM_Cxx2(i,rows,k) = 0.0
              HM_Cxy1(i,rows,k) = 0.
              HM_Cyx2(i,rows,k) = 0.
              HM_Cyy1(i,n_rows,k) = 0.
              HM_Cyy2(i,n_rows,k) = 0.
              HM_Cxy2(i,n_rows,k) = 0.
              HM_Cyx1(i,n_rows,k) = 0.
              HM_Cyy1(i,rows,k) = 0.
              HM_Cyy2(i,rows,k) = 0.
              HM_Cxy2(i,rows,k) = 0.
              HM_Cyx1(i,rows,k) = 0.
              HM_RHS(i,rows,k) = 0.0
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,rows,k) = 0.0
              HM_Cxp(i,rows,k) = 0.0
              HM_Cyz(i,rows,k) = 0.0
              HM_Cyp(i,rows,k) = 0.0
            End Do
          End Do
        End If
        If(at_extremity(PWest))then
          Do k = 1, model_levels
            Do j = 1, rows
              HM_Cxx1(1,j,k) = 0.0
              HM_Cxx2(1,j,k) = 0.0
              HM_Cxy1(1,j,k) = 0.0
              HM_Cyx2(1,j,k) = 0.0
              HM_Cyy1(1,j,k) = 0.0
              HM_Cyy2(1,j,k) = 0.0
              HM_Cxy2(1,j,k) = 0.0
              HM_Cyx1(1,j,k) = 0.0
              HM_RHS(1,j,k) = 0.0
           End Do
          End Do
        End If
        If(at_extremity(PEast))then
          Do k = 1, model_levels
            Do j = 1, rows
              HM_Cxx1(row_length,j,k) = 0.0
              HM_Cxx2(row_length,j,k) = 0.0
              HM_Cxy1(row_length,j,k) = 0.0
              HM_Cyx2(row_length,j,k) = 0.0

              HM_Cxx1(row_length-1,j,k) = 0.0
              HM_Cxx2(row_length-1,j,k) = 0.0
              HM_Cxy1(row_length-1,j,k) = 0.0
              HM_Cyx2(row_length-1,j,k) = 0.0

              HM_Cyy1(row_length,j,k) = 0.0
              HM_Cyy2(row_length,j,k) = 0.0
              HM_Cxy2(row_length,j,k) = 0.0
              HM_Cyx1(row_length,j,k) = 0.0
              HM_RHS(row_length,j,k) = 0.0
            End Do
          End Do
        End If

      Else If (model_domain  ==  mt_cyclic_LAM) Then
        If (at_extremity(PSouth)) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,1,k) = HM_Cxx1(i,2,k)
              HM_Cxx2(i,1,k) = HM_Cxx2(i,2,k)
              HM_RHS(i,1,k) = HM_RHS(i,2,k)
              HM_Cxy1(i,1,k) = HM_Cxy1(i,2,k)
              HM_Cyx1(i,1,k) = HM_Cyx1(i,2,k)
              HM_Cyx2(i,1,k) = HM_Cyx2(i,2,k)
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,1,k) = HM_Cxz(i,2,k)
              HM_Cxp(i,1,k) = HM_Cxp(i,2,k)
            End Do
          End Do
        End If
        If (at_extremity(PNorth)) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              HM_Cxx1(i,rows,k) = HM_Cxx1(i,rows-1,k)
              HM_Cxx2(i,rows,k) = HM_Cxx2(i,rows-1,k)
              HM_RHS(i,rows,k)  = HM_RHS(i,rows-1,k)
              HM_Cxy1(i,rows,k) = HM_Cxy1(i,rows-1,k)
              HM_Cxy2(i,rows,k) = HM_Cxy2(i,rows-1,k)
              HM_Cyx1(i,rows,k) = HM_Cyx1(i,rows-1,k)
              HM_Cyx2(i,rows,k) = HM_Cyx2(i,rows-1,k)
! set Cyy1, Cyy2 coefficients which lie outside domain
              HM_Cyy1(i,rows,k) = 0.
              HM_Cyy2(i,rows,k) = 0.
            End Do
          End Do
          Do k = 1, first_constant_r_rho_level_m1
            Do i = 1, row_length
              HM_Cxz(i,rows,k) = HM_Cxz(i,rows-1,k)
              HM_Cxp(i,rows,k) = HM_Cxp(i,rows-1,k)
              HM_Cyz(i,rows,k) = 0.0
              HM_Cyp(i,rows,k) = 0.0
            End Do
          End Do
        End If
      End If

      i_field = 0

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Czz(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cz(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_RHS(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxx1(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxx2(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxy1(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxy2(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxp(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =                   &
                                first_constant_r_rho_level_m1
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cxz(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =                   &
                                first_constant_r_rho_level_m1
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyy1(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyy2(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyx1(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyx2(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =  model_levels
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyp(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =                   &
                                first_constant_r_rho_level_m1
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field        => HM_Cyz(:,:,:)
      fields_to_swap(i_field) % field_type   =  fld_type_p
      fields_to_swap(i_field) % levels       =                   &
                                first_constant_r_rho_level_m1
      fields_to_swap(i_field) % rows         =  rows
      fields_to_swap(i_field) % vector       =  .FALSE.

! DEPENDS ON: swap_bounds_mv
      CALL SWAP_BOUNDS_MV(fields_to_swap, i_field,                      &
                          row_length, off_x, off_y)

      IF ( L_print_L2helm ) THEN

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              rhs_cut(i,j,k) = HM_RHS(i,j,k)
            END DO
          END DO
        END DO

      IF( norm_start_lev == norm_end_lev ) THEN

        DO k = 1, model_levels
! DEPENDS ON: helm_l2_print
        CALL helm_l2_print(                                             &
                           HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                           HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,               &
                           HM_Cxp, HM_Cyp, thetav_star,                 &
                           HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
                           HM_C3, HM_C4, HM_C5, rhs_cut,                &
                           row_length, rows, model_levels,              &
                           first_constant_r_rho_level_m1,               &
                           k, k, off_x, off_y,                          &
                           .false., .false., 1, L_print_allpe, me )
        END DO ! k = 1, model_levels

      ELSE ! norms from norm_start_lev to norm_end_lev

! DEPENDS ON: helm_l2_print
        CALL helm_l2_print(                                             &
                           HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                           HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,               &
                           HM_Cxp, HM_Cyp, thetav_star,                 &
                           HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
                           HM_C3, HM_C4, HM_C5, rhs_cut,                &
                           row_length, rows, model_levels,              &
                           first_constant_r_rho_level_m1,               &
                           norm_start_lev, norm_end_lev, off_x, off_y,  &
                           .false., .false., 1, L_print_allpe, me )
      END IF   ! norm_start_lev == norm_end_lev
! DEPENDS ON: um_fort_flush
        IF ( L_flush ) CALL um_fort_flush(6,info)

      END IF ! L_print_L2helm

! Call GCR(k) to solve problem.
! DEPENDS ON: gcr_k
      Call GCR_k(                                                       &
     &                 rows, row_length, model_levels,                  &
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
     &                 GCR_adi_add_full_soln, L_gcr_fast_x,             &
     &                 first_constant_r_rho_level,                      &
     &                 first_constant_r_rho_level_m1,                   &
     &                 delta_lambda, delta_phi,                         &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 FV_sec_theta_latitude, cos_v_latitude,           &
     &                 FV_cos_theta_latitude,                           &
     &                 HM_Cxx1, HM_Cxx2,                                &
     &                 HM_Cyy1, HM_Cyy2,                                &
     &                 HM_Czz, HM_Cz, HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,   &
     &                 HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,              &
     &                 thetav_star, HM_C3, HM_C4, HM_C5,                &
     &                 HM_RHS, rescale,                                 &
     &                 weight_upper, weight_lower, exner_prime,         &
     &                 off_x, off_y, me, n_rows, n_proc,                &
     &                 n_procx, n_procy,                                &
     &                 at_extremity, neighbour, l_datastart,            &
     &                 proc_row_group, global_row_length,               &
     &                 proc_col_group,                                  &
     &                 number_dif_horiz_points, halo_i, halo_j,         &
     &                 global_rows,                                     &
     &                 g_rows, g_row_length,                            &
     &                 ldump                                            &
     &                 )

! DEPENDS ON: swap_bounds
      call swap_bounds(                                                 &
     &                 exner_prime, row_length, rows, model_levels,     &
     &                 off_x, off_y, fld_type_p, .false.)

! un-re-scale HM_cyx2 and HM_cxy2 to remove friction term
      Do k = 1, boundary_layer_levels
        Do j = j0, j1+1
          Do i = 0, row_length
            HM_Cyx2(i,j,k) = HM_Cyx2(i,j,k)                             &
     &                       / recip_1p_cd_lambda(i,j,k)
          End Do
        End Do
        Do j = j0-1,j1
          Do i = 1, row_length+1
            HM_Cxy2(i,j,k) = HM_Cxy2(i,j,k)                             &
     &                       / recip_1p_cd_phi(i,j,k)
          End Do
        End Do
      End Do


! ----------------------------------------------------------------------
! Section 5.   Calculate increments to u,v and w.
! ----------------------------------------------------------------------

      ALLOCATE ( etadot(1-off_x:row_length+off_x,                       &
     &                   1-off_y:rows+off_y, model_levels) )

! ----------------------------------------------------------------------
! Section 5.1  Calculate increments to w.
! ----------------------------------------------------------------------

      Do k = 1, model_levels - 1
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            etadot(i,j,k) = (exner_prime(i,j,k+1) -                     &
     &                       exner_prime(i,j,k) )                       &
     &                    / (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do
      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) Then
          Do j = 1, rows
            Do i = 1, row_length
              R_w(i,j,k) = R_w(i,j,k) * recip_G_term(i,j,k)             &
     &                   - K_term(i,j,k) * etadot(i,j,k)
            End Do
          End Do
        Else If (model_domain  ==  mt_cyclic_LAM) Then
          Do j = j0, j1
            Do i = 1, row_length
              R_w(i,j,k) = R_w(i,j,k) * recip_G_term(i,j,k)             &
     &                   - K_term(i,j,k) * etadot(i,j,k)
            End Do
          End Do
        Else If (model_domain  ==  mt_LAM) Then
! (not updated on LAM boundaries)
          Do j = j0, j1
            Do i = i_start, i_end
              R_w(i,j,k) = R_w(i,j,k) * recip_G_term(i,j,k)             &
     &                   - K_term(i,j,k) * etadot(i,j,k)
            End Do
          End Do
        End If
      End Do

! ----------------------------------------------------------------------
! Section 5.2  Calculate increments to u and v.
! ----------------------------------------------------------------------

        Do  k = 1, model_levels
! calculate first derivatives times Cp theta /r.
          If ( k  >=  first_constant_r_rho_level ) Then
            Do j = j0, j1+1
              Do i = 0, row_length
                interp5(i,j) =  (exner_prime(i+1,j,k) -                 &
     &                           exner_prime(i,j,k)) *                  &
     &                           HM_Cyx2(i,j,k)
              End Do
            End Do

            Do j = j0-1, j1
              Do i = 1, row_length+1
                interp6(i,j) =  (exner_prime(i,j+1,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                          * HM_Cxy2(i,j,k)
              End Do
            End Do

          Else

! form average vertical pressure derivative
            If (k  >   1 ) Then
              Do j = j0-1, j1+1
                Do i = 0, row_length+1
                  interp1(i,j) = weight_upper(i,j,k) *                  &
     &                           thetav_star(i,j,k) * etadot(i,j,k)     &
     &                         + weight_lower(i,j,k)                    &
     &                           * thetav_star(i,j,k-1)                 &
     &                           * etadot(i,j,k-1)
                End Do
              End Do
            Else
              Do j = j0-1, j1+1
                Do i = 0, row_length+1
                  interp1(i,j) = weight_upper(i,j,k) *                  &
     &                           thetav_star(i,j,k) * etadot(i,j,k)
                End Do
              End Do
            End If

            Do j = j_one, j1p1
              Do i = 0, row_length
                interp5(i,j) = ((exner_prime(i+1,j,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                           - HM_Cxp(i,j,k) * .5*                  &
     &                             (interp1(i,j)+interp1(i+1,j)) )      &
     &                          * HM_Cyx2(i,j,k)
              End Do
            End Do

            Do j = j0-1, j1
              Do i = 1, row_length+1
                interp6(i,j) = ((exner_prime(i,j+1,k) -                 &
     &                           exner_prime(i,j,k))                    &
     &                           - HM_Cyp(i,j,k) * .5*                  &
     &                             (interp1(i,j)+interp1(i,j+1)) )      &
     &                           * HM_Cxy2(i,j,k)
              End Do
            End Do

          End If

          If (model_domain  ==  mt_Global) Then
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

          Else If (model_domain  ==  mt_LAM) Then
! zero derivatives at boundaries
            If (at_extremity(PEast)) Then
              Do j = 1-off_y, rows+off_y
                interp5(row_length,j) = 0.
                interp5(row_length-1,j) = 0.
              End Do
            End If
            If (at_extremity(PWest)) Then
              Do j = 1-off_y, rows+off_y
                interp5(1,j) = 0.
              End Do
            End If
          End If
          If (model_domain  ==  mt_LAM .or.                             &
     &        model_domain  ==  mt_cyclic_LAM) Then
            If (at_extremity(PSouth)) Then
              Do i = 1-off_x, row_length+off_x
                interp6(i,1) = 0.
!                interp5(i,1) = 0.
              End Do
            End If
            If (at_extremity(PNorth)) Then
              Do i = 1-off_x, row_length+off_x
!                interp5(i,rows) = 0.
                interp6(i,rows-1) = 0.
              End Do
            End If
          End If

! v increment
! first calculate lambda derivative of exner_prime, store in interp1
! calculate R_u at phi points, store in interp2
! Final calculation of v_prime is after u increment calculation.
          Do j = j_start, j_end+1
            Do i = i_start, i_end
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
          Do j = j0, j1
            Do i = i_start, i_end_u
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

          If (model_domain  ==  mt_Global) Then
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

          End If

          Do j = j_start, j_end
            Do i = i_start, i_end
              R_v(i,j,k) = Av(i,j,k) *                                  &
     &                         (R_v(i,j,k) - temp1 * interp6(i,j) )     &
     &                       - Fv(i,j,k) *                              &
     &                         .5 * (interp2(i,j) + interp2(i,j+1)      &
     &                               - temp1 *                          &
     &                           (interp1(i,j) + interp1(i,j+1) ) )
            End Do
          End Do

! end loop over levels
        End Do

        DEALLOCATE ( etadot )

! End of routine.
      IF (lhook) CALL dr_hook('PE_HELMHOLTZ_EUL_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PE_Helmholtz_eul_mix
