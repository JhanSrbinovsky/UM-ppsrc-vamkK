! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_SL_Thermo
!
      Subroutine NI_SL_Thermo(                                          &
     &                         moisture_array_size,                     &
     &                         theta, q, qcl, qcf,                      &
     &                         qcf2, qrain, qgraup,                     &
     &                         mix_v, mix_cl, mix_cf,                   &
     &                         mix_cf2, mix_rain, mix_graup,            &
     &                         cf_pc2, cfl_pc2, cff_pc2,                &
     &                         e_trb, tsq, qsq, cov,                    &
     &                         q_star, qcl_star, qcf_star,              &
     &                         qcf2_star, qrain_star, qgraup_star,      &
     &                         mix_v_star, mix_cf_star, mix_cl_star,    &
     &                         mix_cf2_star, mix_rain_star,             &
     &                         mix_graup_star, cf_pc2_star,             &
     &                         cfl_pc2_star, cff_pc2_star,              &
     &                         exner_star, theta_star, theta_np1,       &
     &                         w, w_adv, u_adv, v_adv,                  &
     &                         exner_theta_levels,                      &
     &                         p_star, p, p_theta_levels, rho,          &
     &                         eta_rho_levels, eta_theta_levels,        &
     &                         r_rho_levels, r_theta_levels,            &
     &                         row_length, rows, n_rows,                &
     &                         model_levels, wet_model_levels,bl_levels,&
     &                         alpha_2, check_bottom_levels,            &
     &                         interp_vertical_search_tol,              &
     &                         first_constant_r_rho_level,              &
     &                         delta_lambda, delta_phi,                 &
     &                         glambda_p, phi_p, glambda_u, phi_v,      &
     &                         gdlambda_p, dphi_p, gdlambda_u, dphi_v,  &
     &                         grecip_dlamp, recip_dphip, grecip_dlamu, &
     &                         recip_dphiv, wt_lambda_p, wt_phi_p,      &
     &                         wt_lambda_u, wt_phi_v, lambda_p_rm,      &
     &                         lambda_p_rp, lambda_u_rm, lambda_u_rp,   &
     &                         phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,  &
     &                         recip_lambda_p_m, recip_lambda_p_0,      &
     &                         recip_lambda_p_p, recip_lambda_p_p2,     &
     &                         recip_lambda_u_m, recip_lambda_u_0,      &
     &                         recip_lambda_u_p, recip_lambda_u_p2,     &
     &                         recip_phi_p_m, recip_phi_p_0,            &
     &                         recip_phi_p_p, recip_phi_p_p2,           &
     &                         recip_phi_v_m, recip_phi_v_0,            &
     &                         recip_phi_v_p, recip_phi_v_p2,           &
     &                         base_lambda, base_phi,                   &
     &                         lambda_p_end, phi_p_end,                 &
     &                         dlambda_p_end, dphi_p_end, dphi_v_end,   &
     &                         recip_dlam, recip_dphi, max_look,        &
     &                         look_lam, look_phi, halo_lam, halo_phi,  &
     &                         timestep,  FV_cos_theta_latitude,        &
     &                           cos_theta_latitude, sec_theta_latitude,&
     &                           sin_theta_latitude, tan_theta_latitude,&
     &                           cos_v_latitude, sec_v_latitude,        &
     &                           sin_v_latitude, tan_v_latitude,        &
     &                           sin_theta_longitude,                   &
     &                           cos_theta_longitude,                   &
     &                           LAM_max_cfl, theta_lbcs,               &
     &                           rimwidth, rimweights,                  &
     &                           lenrim, lbc_size, lbc_start,           &
     &                           r_at_u, r_at_v,                        &
     &                           me, n_proc, n_procx, n_procy,          &
     &                           off_x, off_y, halo_i, halo_j,          &
     &                           datastart, g_i_pe, g_j_pe, l_2dcomm,   &
     &                           size_2dcomm, group_2dcomm,             &
     &                           max_comm_size, at_extremity,           &
     &                           global_row_length, global_rows,        &
     &                           proc_row_group, proc_col_group,        &
     &                           proc_all_group,                        &
     &                           Depart_scheme, Depart_order,           &
     &                           high_order_scheme_theta,               &
     &                           monotone_scheme_theta,                 &
     &                           high_order_scheme_moist,               &
     &                           monotone_scheme_moist,                 &
     &                           L_ritchie_high, L_ritchie_mono,        &
     &                           ritchie_high_order_scheme,             &
     &                           ritchie_monotone_scheme,               &
     &                           model_domain, L_high_theta,            &
     &                           L_mono_theta, thmono_levels,           &
     &                           L_thmono_fixed,                        &
     &                           L_high_moist, L_mono_moist,            &
     &                           L_conserv_moist, L_pc2,                &
     &                           L_2d_sl_geometry, L_mix_ratio,         &
     &                           L_mcr_cf2, L_mcr_rain, L_mcr_graup,    &
     &                           L_regular, L_sl_halo_reprod,           &
     &                           L_fint_theta, L_lbc_old,               &
     &                           L_new_tdisc, CycleNo,                  &
     &                           depart_r_theta,                        &
     &                           depart_lambda, depart_phi, depart_r_w, &
     &                           Error_Code )


! Purpose: Interface routine to SL_Thermo
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
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                     ! number of points on a row
     &, rows                                                            &
                     ! number of rows.
     &, n_rows                                                          &
                     ! number of v rows.
     &, model_levels                                                    &
                     ! Number of model levels.
     &, wet_model_levels                                                &
                         ! Number of model levels where moisture held
     &, bl_levels                                                       &
                         ! Number of model levels of boundary layer
     &, me                                                              &
                     ! My processor number
     &, n_proc                                                          &
                     ! Total number of processors
     &, n_procx                                                         &
                     ! Number of processors in longitude
     &, n_procy                                                         &
                     ! Number of processors in latitude
     &, halo_i                                                          &
                     ! Size of large halo in i.
     &, halo_j                                                          &
                     ! Size of large halo in j.
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, thmono_levels                                                   &
     &, datastart(3)                                                    &
                     ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, proc_all_group                                                  &
                       ! Group id for all processors
     &, max_look                                                        &
                            ! max size of look-up arrays for searches
     &, global_row_length                                               &
                            ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of points in a column
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                             ! processor on my processor-row
                             ! holding a given value in i direction
     &, g_j_pe(1-halo_j:global_rows      +halo_j)                       &
                             ! processor on my processor-column
                             ! holding a given value in j direction
     &, rimwidth                                                        &
                      ! Width of boundaries in LBCs
     &, rimweights(rimwidth)                                            &
                      ! Weights to apply to the LBCs
     &, lenrim                                                          &
                      ! Size of single level of LBC data
     &, lbc_size(4)                                                     &
                      ! Size of each side of LBC data
     &, lbc_start(4)                                                    &
                      ! Start of each side in LBC data
     &, moisture_array_size
      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: max_comm_size ! error check size for comms on demand

      Integer                                                           &
     &  first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, CycleNo

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_regular                                                       &
     &, L_new_tdisc                                                     &
     &, L_lbc_old                                                       &
                         ! false for new lbc treatment
     &, L_fint_theta                                                    &
                     ! true:  fully-interpolating semi-lagrangian
                     !        theta advection will be used
                     ! false: standard non-interpolating in the vertical
     &, L_thmono_fixed  
!                    ! true to use corrected intelligent limiter code
!                    ! false to use code operational up to vn7.5

      Integer                                                           &
     &  LAM_max_cfl(2)

      Integer                                                           &
     &  high_order_scheme_theta                                         &
                                 ! a code saying which high order
                           ! scheme to use for theta.
     &, high_order_scheme_moist                                         &
                                 ! a code saying which high order
                           ! scheme to use for moist variables.
     &, monotone_scheme_theta                                           &
                              ! a code saying which monotone
                           ! scheme to use for theta.
     &, monotone_scheme_moist                                           &
                              ! a code saying which monotone
                           ! scheme to use for moist variables.
     &, Depart_scheme                                                   &
                        ! code saying which departure point scheme to
                        ! use.
     &, Depart_order                                                    &
                        ! for the chosen departure point scheme how
                        ! many iterations/terms to use.
     &, ritchie_high_order_scheme                                       &
                                  ! code choosing high order
                          ! interpolation scheme used in ritchie routine
     &, ritchie_monotone_scheme                                         &
                                  ! code choosing monotone
                          ! interpolation scheme used in ritchie routine
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Logical                                                           &
     &  L_high_theta                                                    &
                       ! True, if high order interpolation required
                       !       for theta.
     &, L_mono_theta                                                    &
                       ! True, if interpolation required to be monotone
                       !       for theta.
     &, L_high_moist                                                    &
                       ! True, if high order interpolation required
                       !       for moist variables.
     &, L_mono_moist                                                    &
                       ! True, if interpolation required to be monotone
                       !       for moist variables.
     &, L_conserv_moist                                                 &
                        ! True, if interpolation to be monotone and
                       !       conservative for moist variables.
     &, L_pc2                                                           &
                       ! True, if clouds need advecting for PC2 scheme
     &, L_mix_ratio                                                     &
     &, L_mcr_cf2                                                       &
                       ! True if second ice variable is in use
     &, L_mcr_rain                                                      &
                       ! True if prognostic rain is in use
     &, L_mcr_graup    ! True if prognostic graupel is in use

      Logical                                                           &
     &  L_ritchie_high                                                  &
                       ! True if high order interpolation scheme to be
                      ! used in ritchie scheme
     &, L_ritchie_mono                                                  &
                       ! True if monotone interpolation scheme to be
                      ! used in ritchie scheme
     &, L_2d_sl_geometry

      Real                                                              &
     &  delta_lambda                                                    &
                      ! holds spacing between points in the i
                      ! direction for the input data field.
     &, delta_phi                                                       &
                      ! holds spacing between points in the j
                      ! direction for the input data field.
     &, base_lambda                                                     &
     &, base_phi                                                        &
     &, lambda_p_end                                                    &
     &, phi_p_end                                                       &
     &, dlambda_p_end                                                   &
     &, dphi_p_end                                                      &
     &, dphi_v_end                                                      &
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, timestep                                                        &
     &, alpha_2

! look-up table halos
      Integer                                                           &
     &  halo_lam                                                        &
     &, halo_phi

!VarRes horizontal co-ordinate look-up table
      Integer                                                           &
     &  look_lam(1-halo_lam : max_look-halo_lam)                        &
     &, look_phi(1-halo_phi : max_look-halo_phi)

!VarRes horizontal co-ordinate information etc.
      Real                                                              &
     &  glambda_p( 1-halo_i : global_row_length+halo_i)                 &
     &, glambda_u( 1-halo_i : global_row_length+halo_i)                 &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )                           &
     &, gdlambda_p(1-halo_i : global_row_length+halo_i)                 &
     &, gdlambda_u(1-halo_i : global_row_length+halo_i)                 &
     &, dphi_p  ( 1-halo_i : row_length + halo_i                        &
     &,           1-halo_j : rows+halo_j )                              &
     &, dphi_v  ( 1-halo_i : row_length + halo_i                        &
     &,           1-halo_j : n_rows+halo_j )                            &
     &, grecip_dlamp(1-halo_i : global_row_length + halo_i)             &
     &, grecip_dlamu(1-halo_i : global_row_length + halo_i)             &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, recip_dphiv( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)    &
     &, lambda_p_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_p_rp (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rp (1-halo_i : row_length + halo_i)                    &
     &, phi_p_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_p_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_v_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, phi_v_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, recip_lambda_p_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_phi_p_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_v_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )

! 1st iteration estimate for theta^n+1
      Real                                                              &
     &  theta_np1(1-off_x:row_length+off_x,                             &
     &             1-off_y:rows+off_y, model_levels)

      Real, Intent (InOut) ::                                           &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     wet_model_levels)                                            &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_model_levels)                                         &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         wet_model_levels)                                        &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, cf_pc2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &           wet_model_levels)                                      &
     &, cfl_pc2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_model_levels)                                      &
     &, cff_pc2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_model_levels)                                      &
     &, e_trb (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         bl_levels)                                               &
     &, tsq   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         bl_levels)                                               &
     &, qsq   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         bl_levels)                                               &
     &, cov   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         bl_levels)


      Real                                                              &
     &  eta_rho_levels(model_levels)                                    &
     &, eta_theta_levels(0:model_levels)

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, u_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, v_adv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &         model_levels)                                            &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         0:model_levels)                                          &
     &, w (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     0:model_levels)                                              &
     &, rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Real                                                              &
     &  p_star (row_length, rows)                                       &
     &, p (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)  &
     &, p_theta_levels (1-off_x:row_length+off_x,                       &
     &                  1-off_y:rows+off_y, model_levels)               &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)

      Real                                                              &
                      ! Trig functions.
     &  FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sin_theta_latitude (row_length, rows)                           &
     &, tan_theta_latitude (row_length, rows)                           &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, tan_v_latitude (row_length, n_rows)                             &
     &, cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

      Real                                                              &
     &  theta_lbcs(lenrim,model_levels)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  depart_lambda (row_length, rows, model_levels)                  &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_phi (row_length, rows, model_levels)                     &
                                                    ! Phi Co-ordinate
                                                     ! of co-ordinate of
                                                      ! departure point.
     &, depart_r_w (row_length, rows, model_levels)                     &
                                                      ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.
     &, exner_star(1-off_x:row_length+off_x,                            &
                                                      ! Departure value
     &             1-off_y:rows+off_y, model_levels)                    &
                                                      ! of exner.
     &, depart_r_theta (row_length, rows, model_levels)   ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.

! Star variables hold physics increments on input, latest values on
! output
      Real, Intent (InOut) ::                                           &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, model_levels)                    &
     &, q_star(1-off_x:row_length+off_x,                                &
     &         1-off_y:rows+off_y, wet_model_levels)                    &
     &, qcl_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, qcf_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, qcf2_star(1-off_x:row_length+off_x,                             &
     &            1-off_y:rows+off_y, wet_model_levels)                 &
     &, qrain_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, qgraup_star(1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y, wet_model_levels)               &
     &, cf_pc2_star (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)              &
     &, cfl_pc2_star(1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)              &
     &, cff_pc2_star(1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)

      Real, Intent (InOut) ::                                           &
     &  mix_v (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cf(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_cf2 (1-halo_i:row_length+halo_i,                            &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_rain(1-halo_i:row_length+halo_i,                            &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_graup(1-halo_i:row_length+halo_i,                           &
     &         1-halo_j:rows+halo_j, wet_model_levels)                  &
     &, mix_v_star  (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_cl_star (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_cf_star (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_cf2_star  (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_star (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_star (1-off_x:row_length+off_x,                       &
     &               1-off_y:rows+off_y,wet_model_levels)

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local Variables.

      Real                                                              &
     &  mono_mass(row_length, rows, model_levels)

! scalars
      Real                                                              &
     &  drkp1, drk

      Integer                                                           &
     &  i, j, k   ! loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1  Mixing ratio set up first
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('NI_SL_THERMO',zhook_in,zhook_handle)
      If ( L_mix_ratio ) then

! ----------------------------------------------------------------------
! Section 2  Set delta_p at data points, ensure positive.
! ----------------------------------------------------------------------

        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
            drk   = r_theta_levels(i,j,k)   - r_theta_levels(i,j,k-1)
            mono_mass(i,j,k) = rho(i,j,k+1) * drkp1 + rho(i,j,k) * drk
            mono_mass(i,j,k) = mono_mass(i,j,k)/r_theta_levels(i,j,k)**2
          End Do
        End Do
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
              drk   = r_rho_levels(i,j,k)     - r_theta_levels(i,j,k-1)
              mono_mass(i,j,k) = rho(i,j,k+1) * drkp1 + rho(i,j,k) * drk
              mono_mass(i,j,k) = mono_mass(i,j,k)/                      &
     &                                         r_theta_levels(i,j,k)**2
            End Do
          End Do
        End Do  !  k = 2, model_levels - 1
        k = model_levels
        Do j = 1, rows
          Do i = 1, row_length
            drk   = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            mono_mass(i,j,k) = rho(i,j,k) * drk
            mono_mass(i,j,k) = mono_mass(i,j,k)/r_theta_levels(i,j,k)**2
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 3   Call SL_Thermo
! ----------------------------------------------------------------------

! DEPENDS ON: sl_thermo
        Call SL_Thermo(                                                 &
     &                 moisture_array_size, timestep,                   &
     &                 theta, mix_v, mix_cl, mix_cf,                    &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 cf_pc2, cfl_pc2, cff_pc2,                        &
     &                 e_trb, tsq, qsq, cov,                            &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 cf_pc2_star, cfl_pc2_star, cff_pc2_star,         &
     &                 exner_star, theta_star, theta_np1,               &
     &                 u_adv, v_adv, w_adv, w,                          &
     &                 eta_rho_levels, eta_theta_levels,                &
     &                 r_rho_levels, r_theta_levels,                    &
     &                 exner_theta_levels, mono_mass,                   &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_model_levels, bl_levels,       &
     &                 alpha_2, check_bottom_levels,                    &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, glambda_u, phi_v,              &
     &                 gdlambda_p, dphi_p, gdlambda_u, dphi_v,          &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 grecip_dlamu, recip_dphiv,                       &
     &                 wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,    &
     &                 lambda_p_rm, lambda_p_rp,                        &
     &                 lambda_u_rm, lambda_u_rp,                        &
     &                 phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,          &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_lambda_u_m, recip_lambda_u_0,              &
     &                 recip_lambda_u_p, recip_lambda_u_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 recip_phi_v_m, recip_phi_v_0,                    &
     &                 recip_phi_v_p, recip_phi_v_p2,                   &
     &                 Base_lambda, base_phi,                           &
     &                 lambda_p_end, phi_p_end,                         &
     &                 dlambda_p_end, dphi_p_end, dphi_v_end,           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 cos_theta_latitude, sec_theta_latitude,          &
     &                 sin_theta_latitude, tan_theta_latitude,          &
     &                 cos_v_latitude, sec_v_latitude,                  &
     &                 sin_v_latitude, tan_v_latitude,                  &
     &                 sin_theta_longitude,                             &
     &                 cos_theta_longitude,                             &
     &                 LAM_max_cfl, theta_lbcs,                         &
     &                 rimwidth, rimweights,                            &
     &                 lenrim, lbc_size, lbc_start,                     &
     &                 r_at_u, r_at_v,                                  &
     &                 me, n_proc, n_procx, n_procy,                    &
     &                 datastart, g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,&
     &                 group_2dcomm, max_comm_size, at_extremity,       &
     &                 global_row_length, global_rows,                  &
     &                 proc_row_group, proc_col_group, proc_all_group,  &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 Depart_scheme,                                   &
     &                 Depart_order,                                    &
     &                 high_order_scheme_theta,                         &
     &                 monotone_scheme_theta,                           &
     &                 high_order_scheme_moist,                         &
     &                 monotone_scheme_moist,                           &
     &                 L_Ritchie_high,                                  &
     &                 L_Ritchie_mono,                                  &
     &                 Ritchie_high_order_scheme,                       &
     &                 Ritchie_monotone_scheme,                         &
     &                 model_domain, L_high_theta,                      &
     &                 L_mono_theta, thmono_levels,                     &
     &                 L_thmono_fixed,                                  &
     &                 L_high_moist, L_mono_moist,                      &
     &                 L_conserv_moist, L_pc2,                          &
     &                 L_2d_sl_geometry, L_sl_halo_reprod,              &
     &                 L_regular, L_fint_theta,                         &
     &                 L_mcr_cf2, L_mcr_rain, L_mcr_graup,              &
     &                 L_lbc_old, L_new_tdisc, CycleNo,                 &
     &                 depart_r_theta,                                  &
     &                 depart_lambda, depart_phi, depart_r_w,           &
     &                 Error_Code)

      else   ! L_mix_ratio=.false.

! ----------------------------------------------------------------------
! Section 4    Specific quantities
!              Set delta_p at data points, ensure positive.
! ----------------------------------------------------------------------

        k = 1
        Do j = 1, rows
          Do i = 1, row_length
            mono_mass(i,j,k) = p_star(i,j) - p(i,j,2)
          End Do
        End Do
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              mono_mass(i,j,k) = p(i,j,k) - p(i,j,k+1)
            End Do
          End Do
        End Do
        k = model_levels
        Do j = 1, rows
          Do i = 1, row_length
            mono_mass(i,j,k) = p(i,j,k) - p_theta_levels(i,j,k)
          End Do
        End Do

! DEPENDS ON: sl_thermo
        Call SL_Thermo(                                                 &
     &                 moisture_array_size, timestep,                   &
     &                 theta, q, qcl, qcf,                              &
     &                 qcf2, qrain, qgraup,                             &
     &                 cf_pc2, cfl_pc2, cff_pc2,                        &
     &                 e_trb, tsq, qsq, cov,                            &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 cf_pc2_star, cfl_pc2_star, cff_pc2_star,         &
     &                 exner_star, theta_star, theta_np1,               &
     &                 u_adv, v_adv, w_adv, w,                          &
     &                 eta_rho_levels, eta_theta_levels,                &
     &                 r_rho_levels, r_theta_levels,                    &
     &                 exner_theta_levels, mono_mass,                   &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_model_levels, bl_levels,       &
     &                 alpha_2, check_bottom_levels,                    &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, glambda_u, phi_v,              &
     &                 gdlambda_p, dphi_p, gdlambda_u, dphi_v,          &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 grecip_dlamu, recip_dphiv,                       &
     &                 wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,    &
     &                 lambda_p_rm, lambda_p_rp,                        &
     &                 lambda_u_rm, lambda_u_rp,                        &
     &                 phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,          &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_lambda_u_m, recip_lambda_u_0,              &
     &                 recip_lambda_u_p, recip_lambda_u_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 recip_phi_v_m, recip_phi_v_0,                    &
     &                 recip_phi_v_p, recip_phi_v_p2,                   &
     &                 Base_lambda, base_phi,                           &
     &                 lambda_p_end, phi_p_end,                         &
     &                 dlambda_p_end, dphi_p_end, dphi_v_end,           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 cos_theta_latitude, sec_theta_latitude,          &
     &                 sin_theta_latitude, tan_theta_latitude,          &
     &                 cos_v_latitude, sec_v_latitude,                  &
     &                 sin_v_latitude, tan_v_latitude,                  &
     &                 sin_theta_longitude,                             &
     &                 cos_theta_longitude,                             &
     &                 LAM_max_cfl, theta_lbcs,                         &
     &                 rimwidth, rimweights,                            &
     &                 lenrim, lbc_size, lbc_start,                     &
     &                 r_at_u, r_at_v,                                  &
     &                 me, n_proc, n_procx, n_procy,                    &
     &                 datastart, g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,&
     &                 group_2dcomm, max_comm_size, at_extremity,       &
     &                 global_row_length, global_rows,                  &
     &                 proc_row_group, proc_col_group, proc_all_group,  &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 Depart_scheme,                                   &
     &                 Depart_order,                                    &
     &                 high_order_scheme_theta,                         &
     &                 monotone_scheme_theta,                           &
     &                 high_order_scheme_moist,                         &
     &                 monotone_scheme_moist,                           &
     &                 L_Ritchie_high,                                  &
     &                 L_Ritchie_mono,                                  &
     &                 Ritchie_high_order_scheme,                       &
     &                 Ritchie_monotone_scheme,                         &
     &                 model_domain, L_high_theta,                      &
     &                 L_mono_theta, thmono_levels,                     &
     &                 L_thmono_fixed,                                  &
     &                 L_high_moist, L_mono_moist,                      &
     &                 L_conserv_moist, L_pc2,                          &
     &                 L_2d_sl_geometry, L_sl_halo_reprod,              &
     &                 L_regular, L_fint_theta,                         &
     &                 L_mcr_cf2, L_mcr_rain, L_mcr_graup,              &
     &                 L_lbc_old, L_new_tdisc, CycleNo,                 &
     &                 depart_r_theta,                                  &
     &                 depart_lambda, depart_phi, depart_r_w,           &
     &                 Error_Code)
  
      endif  ! L_mix_ratio

! End of routine.
      IF (lhook) CALL dr_hook('NI_SL_THERMO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_SL_Thermo

