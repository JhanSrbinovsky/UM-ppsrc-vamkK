! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
        SUBROUTINE NI_PE_Helmholtz(                                     &
     &                      u, v, w, r_theta_levels,                    &
     &                      r_rho_levels, p, rho, rho_np1, theta,       &
     &                      theta_star, theta_np1, q, qcl, qcf,         &
     &                      qcf2, qrain, qgraup,                        &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_v_star, mix_cl_star, mix_cf_star,       &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      q_star, q_np1, qcl_star, qcf_star,          &
     &                      qcf2_star, qrain_star, qgraup_star,         &
     &                      qcl_np1, qcf_np1,                           &
     &                      qcf2_np1, qrain_np1, qgraup_np1,            &
     &                      rho_Km, cH, G_term_tol,                     &
     &                      exner_rho_levels, exner_theta_levels,       &
     &                      frictional_timescale,                       &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      FV_cos_theta_latitude,                      &
     &                      FV_sec_theta_latitude,                      &
     &                      f3_at_u, f3_at_v,                           &
     &                      timestep, timestep_number,                  &
     &                      row_length, rows, n_rows,                   &
     &                      model_levels, wet_levels, bl_levels,        &
     &                      delta_lambda, delta_phi,                    &
     &                      lambda_p, phi_p, lambda_u, phi_v,           &
     &                      dlambda_p, dphi_p, dlambda_u, dphi_v,       &
     &                      recip_dlamp, recip_dphip,                   &
     &                      recip_dlamu, recip_dphiv,                   &
     &                      wt_lambda_p, wt_phi_p,                      &
     &                      wt_lambda_u, wt_phi_v,                      &
     &                      GCR_max_iterations, GCR_diagnostics,        &
     &                 GCR_its_switch, GCR_its_avg_step,                &
     &                 GCR_max_its, GCR_min_its, GCR_sum_its,           &
     &                 GCR_max_time, GCR_min_time,                      &
     &                      GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
     &                      GCR_zero_init_guess,                        &
     &                      GCR_use_residual_Tol,                       &
     &                      GCR_adi_add_full_soln, L_gcr_fast_x,        &
     &                      GCR_precon_option, GCR_ADI_Pseudo_timestep, &
     &                      GCR_n_ADI_pseudo_timesteps,                 &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      alpha_1, alpha_2, alpha_3, alpha_4,         &
     &                      alpha_Cd,                                   &
     &                      model_domain, L_physics,                    &
     &                      GCR_Restart_value,                          &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
     &                      exner_lbcs, LENRIM, LBC_SIZE, LBC_START,    &
     &                      RIMWIDTH, RIMWEIGHTS,                       &
!!! parallel variables
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      halo_i, halo_j, datastart,                  &
     &                      L_regular, at_extremity,                    &
     &                      rims_to_do, off_x, off_y,                   &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      global_row_length, global_rows,             &
     &                      g_rows, g_row_length, CycleNo,              &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      L_mix_ratio, L_fint_theta, L_new_tdisc,     &
     &                      L_qwaterload, L_lbc_new, L_fixed_lbcs,      &
     &                      L_transparent, ldump                        &
     &                      )

! Purpose: Interface routine to PE_Helmholtz
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE dynamics_input_mod, ONLY :  i_ND_solver_vn   ! i.e. A10_2A/2B switch
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE Field_Types
      USE UM_ParParams
      IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.


      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_levels                                                      &
                   ! number of model levels where moisture is held
     &, bl_levels ! number of boundary layer levels.

      Integer                                                           &
     &  mype                                                            &
                       ! My processor number
     &, nproc                                                           &
                    ! Total number of processors
     &, nproc_x                                                         &
                     ! Number of processors in longitude
     &, nproc_y                                                         &
                     ! Number of processors in latitude
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, datastart(3)                                                    &
                     ! First gridpoints held by this processor.
     &, gc_proc_row_group                                               &
                          ! Group id for processors on the same row
     &, n_rows                                                          &
                   ! Local number of rows in a v field
     &, off_x,off_y                                                     &
     &, global_row_length                                               &
     &, gc_proc_col_group                                               &
     &, global_rows                                                     &
     &, rims_to_do                                                      &
                           ! zone where rim wts = 1
     &, RIMWIDTH                                                        &
                          ! IN : Size of boundary region
     &, LENRIM(Nfld_max,NHalo_max)                                      &
                          ! IN : Size of single level of LBC
     &, LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
                          ! IN : Size of a side of a LBC
     &, LBC_START(4,Nfld_max,NHalo_max)                                 &
                          ! IN : Start of a side in a LBC
     &, norm_lev_start                                                  &
                         !  start level for L2norm diagnostics
     &, norm_lev_end                                                    &
                         !  end level for L2norm diagnostics
     &, CycleNo

      Integer                                                           &
     &  g_rows(0:nproc-1)                                               &
     &, g_row_length(0:nproc-1)

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                          ! true for regular resolution
     &, L_lbc_new                                                       &
     &, L_fixed_lbcs                                                    &
     &, L_transparent                                                   &
     &, L_new_tdisc                                                     &
     &, L_fint_theta
                     ! true:  fully-interpolating semi-lagrangian
                     !        theta advection will be used
                     ! false: standard non-interpolating in the vertical

      Integer                                                           &
     &  first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1 ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, timestep_number

      Logical                                                           &
     &  L_Physics   ! true if physics wanted


      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, timestep                                                        &
     &, G_term_tol       ! tolerance for vertical G term

      Real                                                              &
           !VarRes horizontal co-ordinate information etc.
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
     &, alpha_Cd(bl_levels)

! trigonometric functions
      Real                                                              &
     &  cos_theta_latitude(1-off_x:row_length+off_x,1-off_y:rows+off_y) &
     &, sec_theta_latitude(1-off_x:row_length+off_x,1-off_y:rows+off_y) &
     &, cos_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, sec_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
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

      REAL                                                              &
     &  RIMWEIGHTS(RIMWIDTH)  ! IN : weight to apply to LBC

      REAL                                                              &
     & exner_lbcs(LENRIM(fld_type_p,halo_type_extended),MODEL_LEVELS+1)

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
     &              model_levels)

      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)&
     &, R_v (1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels)&
     &, R_w (row_length, rows, model_levels)                            &
     &, exner_rho_levels (1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
     &                    model_levels)                                 &
     &, exner_theta_levels (1-off_x:row_length+off_x,1-off_y:rows+off_y,&
     &                      model_levels)                               &
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
     &, GCR_max_iterations                                              &
     &, GCR_precon_option                                               &
     &, GCR_n_ADI_pseudo_timesteps

! Physics arrays
      Real                                                              &
     &  rho_Km(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         0:bl_levels-1)                                           &
     &, cH(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels-1)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  exner_prime (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               model_levels)                                      &
     &, dtheta_dr_term (row_length, rows, model_levels)


      Logical                                                           &
     & L_mix_ratio                                                      &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup
      Logical :: L_qwaterload               ! add waterloading terms

      Real, Intent (InOut) ::                                           &
                              ! primary model variables
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     wet_levels)                                                  &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_levels)                                                  &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_levels)                                                  &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &     wet_levels)                                                  &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_levels)                                                  &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &     wet_levels)                                                  &
     &, q_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          wet_levels)                                             &
     &, q_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          wet_levels)                                             &
     &, qcl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_levels)                                             &
     &, qcf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_levels)                                             &
     &, qcf2_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_levels)                                             &
     &, qrain_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_levels)                                             &
     &, qgraup_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_levels)                                             &
     &, qcl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_levels)                                             &
     &, qcf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_levels)                                             &
     &, qcf2_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_levels)                                             &
     &, qrain_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_levels)                                             &
     &, qgraup_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          wet_levels)                                             &
     &, mix_v_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_levels)                                             &
     &, mix_cl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_levels)                                             &
     &, mix_cf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_levels)                                             &
     &, mix_cf2_np1 (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_levels)                     &
     &, mix_rain_np1 (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_levels)                     &
     &, mix_graup_np1 (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_levels)                     &
     &, mix_v (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cf(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cf2 (1-halo_i:row_length+halo_i,                            &
     &           1-halo_j:rows+halo_j, wet_levels)                      &
     &, mix_rain(1-halo_i:row_length+halo_i,                            &
     &           1-halo_j:rows+halo_j, wet_levels)                      &
     &, mix_graup(1-halo_i:row_length+halo_i,                           &
     &            1-halo_j:rows+halo_j, wet_levels)                     &
     &, mix_v_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y,wet_levels)                       &
     &, mix_cl_star(1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y,wet_levels)                      &
     &, mix_cf_star(1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y,wet_levels)                      &
     &, mix_cf2_star(1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_levels)                     &
     &, mix_rain_star(1-off_x:row_length+off_x,                         &
     &                1-off_y:rows+off_y,wet_levels)                    &
     &, mix_graup_star(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y,wet_levels)

      Integer                                                           &
     &  info

      INTEGER, PARAMETER :: i_solver_vn_10A_2A = 1
      INTEGER, PARAMETER :: i_solver_vn_10A_2B = 2
      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'NI_PE_HELMHOLTZ'

      INTEGER            :: ErrorStatus
      CHARACTER (LEN=256):: Cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      

      IF (lhook) CALL dr_hook('NI_PE_HELMHOLTZ',zhook_in,zhook_handle)

      SELECT CASE ( i_ND_solver_vn )

        CASE ( i_solver_vn_10A_2A )      ! A10_2A branch

        IF(L_mix_ratio)THEN
! DEPENDS ON: pe_helmholtz_eul_mix
          CALL PE_Helmholtz_eul_mix(                                      &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho, rho_np1,               &
                             theta, theta_star, theta_np1,               &
                             mix_v, mix_cl, mix_cf,                      &
                             mix_cf2, mix_rain, mix_graup,               &
                             mix_v_star, mix_cl_star, mix_cf_star,       &
                             mix_cf2_star, mix_rain_star, mix_graup_star,&
                             mix_v_np1, mix_cl_np1, mix_cf_np1,          &
                             mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime,                 &
                             dtheta_dr_term, exner_lbcs,                 &
                             LENRIM, LBC_SIZE, LBC_START,                &
                             RIMWIDTH, RIMWEIGHTS,                       &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity,                    &
                             rims_to_do, off_x, off_y,                   &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             L_lbc_new, L_fixed_lbcs, L_transparent,     &
                             ldump )
        
        ELSE IF ( model_domain  ==  mt_global) THEN
        
! DEPENDS ON: PE_Helmholtz_global
          CALL PE_Helmholtz_global(                                       &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho,rho_np1,theta,          &
                             theta_star, theta_np1, q, q_star, q_np1,    &
                             qcl, qcf, qcf2, qrain, qgraup,              &
                             qcl_star, qcf_star,                         &
                             qcf2_star, qrain_star, qgraup_star,         &
                             qcl_np1 , qcf_np1 ,                         &
                             qcf2_np1 , qrain_np1 , qgraup_np1 ,         &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             L_qwaterload,                               &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity, off_x, off_y,      &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             ldump )
        
        ELSE !   L_mix_ratio = false AND NOT Global
        
! DEPENDS ON: pe_helmholtz_eul
          CALL PE_Helmholtz_eul(                                          &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho,rho_np1,theta,          &
                             theta_star, theta_np1, q, q_star, q_np1,    &
                             qcl, qcf, qcf2, qrain, qgraup,              &
                             qcl_star, qcf_star,                         &
                             qcf2_star, qrain_star, qgraup_star,         &
                             qcl_np1 , qcf_np1 ,                         &
                             qcf2_np1 , qrain_np1 , qgraup_np1 ,         &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             L_qwaterload,                               &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
                             exner_lbcs, LENRIM, LBC_SIZE, LBC_START,    &
                             RIMWIDTH, RIMWEIGHTS,                       &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity,                    &
                             rims_to_do, off_x, off_y,                   &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             L_lbc_new, L_fixed_lbcs, L_transparent,     &
                             ldump )
        END IF
        
        CASE ( i_solver_vn_10A_2B )      ! A10_2B branch

        IF(L_mix_ratio)THEN
! DEPENDS ON: pe_helmholtz_eul_mix_2B
          CALL PE_Helmholtz_eul_mix_2B(                                   &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho, rho_np1,               &
                             theta, theta_star, theta_np1,               &
                             mix_v, mix_cl, mix_cf,                      &
                             mix_cf2, mix_rain, mix_graup,               &
                             mix_v_star, mix_cl_star, mix_cf_star,       &
                             mix_cf2_star, mix_rain_star, mix_graup_star,&
                             mix_v_np1, mix_cl_np1, mix_cf_np1,          &
                             mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime,                 &
                             dtheta_dr_term, exner_lbcs,                 &
                             LENRIM, LBC_SIZE, LBC_START,                &
                             RIMWIDTH, RIMWEIGHTS,                       &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity,                    &
                             rims_to_do, off_x, off_y,                   &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             L_lbc_new, L_fixed_lbcs, L_transparent,     &
                             ldump )
        
        ELSE IF ( model_domain  ==  mt_global) THEN
        
! DEPENDS ON: PE_Helmholtz_global_2B
          CALL PE_Helmholtz_global_2B(                                    &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho,rho_np1,theta,          &
                             theta_star, theta_np1, q, q_star, q_np1,    &
                             qcl, qcf, qcf2, qrain, qgraup,              &
                             qcl_star, qcf_star,                         &
                             qcf2_star, qrain_star, qgraup_star,         &
                             qcl_np1 , qcf_np1 ,                         &
                             qcf2_np1 , qrain_np1 , qgraup_np1 ,         &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             L_qwaterload,                               &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity, off_x, off_y,      &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             ldump )
        
        ELSE !   L_mix_ratio = false AND NOT Global
        
! DEPENDS ON: pe_helmholtz_eul_2B
          CALL PE_Helmholtz_eul_2B(                                       &
                             u,v,w, r_theta_levels,                      &
                             r_rho_levels, p,rho,rho_np1,theta,          &
                             theta_star, theta_np1, q, q_star, q_np1,    &
                             qcl, qcf, qcf2, qrain, qgraup,              &
                             qcl_star, qcf_star,                         &
                             qcf2_star, qrain_star, qgraup_star,         &
                             qcl_np1 , qcf_np1 ,                         &
                             qcf2_np1 , qrain_np1 , qgraup_np1 ,         &
                             L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
                             L_qwaterload,                               &
                             rho_Km, cH, G_term_tol,                     &
                             exner_rho_levels, frictional_timescale,     &
                             sec_theta_latitude, cos_v_latitude,         &
                             FV_cos_theta_latitude,                      &
                             FV_sec_theta_latitude,                      &
                             f3_at_u, f3_at_v,                           &
                             timestep, timestep_number,                  &
                             row_length, rows, n_rows,                   &
                             model_levels, wet_levels,                   &
                             bl_levels,                                  &
                             delta_lambda, delta_phi,                    &
                             lambda_p, phi_p, lambda_u, phi_v,           &
                             dlambda_p, dphi_p, dlambda_u, dphi_v,       &
                             recip_dlamp, recip_dphip,                   &
                             recip_dlamu, recip_dphiv,                   &
                             wt_lambda_p, wt_phi_p,                      &
                             wt_lambda_u, wt_phi_v,                      &
                             GCR_max_iterations, GCR_diagnostics,        &
                        GCR_its_switch, GCR_its_avg_step,                &
                        GCR_max_its, GCR_min_its, GCR_sum_its,           &
                        GCR_max_time, GCR_min_time,                      &
                             GCR_tol_res, GCR_tol_abs, GCR_use_tol_abs,  &
                             GCR_zero_init_guess,                        &
                             GCR_use_residual_Tol,                       &
                             GCR_adi_add_full_soln, L_gcr_fast_x,        &
                             GCR_precon_option, GCR_ADI_Pseudo_timestep, &
                             GCR_n_ADI_pseudo_timesteps,                 &
                             eta_theta_levels, eta_rho_levels,           &
                             alpha_1, alpha_2, alpha_3, alpha_4,         &
                             alpha_Cd,                                   &
                             model_domain, L_physics,                    &
                             GCR_Restart_value,                          &
                             first_constant_r_rho_level,                 &
                             first_constant_r_rho_level_m1,              &
                             R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
                             exner_lbcs, LENRIM, LBC_SIZE, LBC_START,    &
                             RIMWIDTH, RIMWEIGHTS,                       &
!!! parallel variables
                             mype, nproc, nproc_x, nproc_y,              &
                             halo_i, halo_j, datastart,                  &
                             L_regular, at_extremity,                    &
                             rims_to_do, off_x, off_y,                   &
                             gc_proc_row_group, gc_proc_col_group,       &
                             global_row_length, global_rows,             &
                             g_rows, g_row_length,                       &
                             CycleNo, L_new_tdisc, L_fint_theta,         &
                             L_lbc_new, L_fixed_lbcs, L_transparent,     &
                             ldump )
        END IF

        CASE DEFAULT

          ErrorStatus = 10
          WRITE (Cmessage,'(A,I6,A)') 'Dynamics solver option ',        &
                             i_ND_solver_vn,' not recognised.'

          CALL Ereport ( RoutineName, ErrorStatus, Cmessage)

        END SELECT

      IF (lhook) CALL dr_hook('NI_PE_HELMHOLTZ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_PE_Helmholtz
