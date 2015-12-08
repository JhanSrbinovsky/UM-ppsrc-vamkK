! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_SL_Full_Wind
!
      Subroutine NI_SL_Full_Wind(                                       &
     &                         u, u_np1, v, v_np1, w, w_np1,            &
     &                         u_adv, v_adv, w_adv,                     &
     &                         theta, theta_np1, exner,                 &
     &                         q, qcl, qcf, qcf2, qrain, qgraup,        &
     &                         q_np1, qcl_np1, qcf_np1, qcf2_np1,       &
     &                         qrain_np1, qgraup_np1,                   &
     &                         mix_v, mix_cl, mix_cf,                   &
     &                         mix_v_np1, mix_cl_np1, mix_cf_np1,       &
     &                         mix_cf2, mix_rain, mix_graup,            &
     &                         mix_cf2_np1, mix_rain_np1, mix_graup_np1,&
     &                         L_mix_ratio,                             &
     &                         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,   &
     &                         depart_lambda_w, depart_phi_w,           &
     &                         depart_r_w,                              &
     &                         r_theta_levels, r_rho_levels,            &
     &                         eta_theta_levels, eta_rho_levels,        &
     &                         cos_theta_latitude, sec_theta_latitude,  &
     &                         sin_theta_latitude, cos_v_latitude,      &
     &                         sec_v_latitude, sin_v_latitude,          &
     &                         tan_theta_latitude, tan_v_latitude,      &
     &                         cos_theta_longitude,                     &
     &                         sin_theta_longitude,                     &
     &                         f1_at_v, f2_at_u, f3_at_u, f3_at_v,      &
     &                         delta_lambda, delta_phi, timestep,       &
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
     &                         alpha_3, alpha_4, LAM_max_cfl,           &
     &                         n_Y_arrays, n_Yw_arrays,                 &
     &                         n_Yd_arrays, n_Ydw_arrays,               &
     &                         u_lbcs, v_lbcs, w_lbcs,                  &
     &                         LENRIM, LBC_SIZE, LBC_START, RIMWIDTH,   &
     &                         model_domain, row_length, rows, n_rows,  &
     &                         model_levels, wet_model_levels,          &
     &                         Depart_scheme, Depart_order,             &
     &                         high_order_scheme, monotone_scheme,      &
     &                         L_trivial_trigs,                         &
     &                         L_high, L_mono, L_conserv,               &
     &                         L_Robert_high, L_Robert_mono,            &
     &                         Robert_high_order_scheme,                &
     &                         Robert_monotone_scheme,                  &
     &                         first_constant_r_rho_level,              &
     &                         check_bottom_levels,                     &
     &                         interp_vertical_search_tol,              &
     &                         r_at_u, r_at_v,                          &
     &                         me, n_proc, n_procx, n_procy,            &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         global_row_length, global_rows,          &
     &                         l_datastart, at_extremity, g_i_pe,       &
     &                         g_j_pe, l_2dcomm, size_2dcomm,           &
     &                         group_2dcomm, proc_row_group,            &
     &                         proc_col_group, proc_all_group,          &
     &                         L_2d_sl_geometry, L_sl_halo_reprod,      &
     &                         L_free_slip, L_regular,                  &
     &                         L_qwaterload, L_interp_depart,           &
     &                         L_lbc_old, L_new_tdisc, CycleNo,         &
     &                         R_u, R_v, R_w, Error_code )


! Purpose:
!          Performs semi-Lagrangian advection of the full vector wind.
!          This is the first advection step in the integration scheme.
!          The second step is dealt with by SL_Full_Wind2.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
!          and
!
!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE atmos_constants_mod, ONLY: recip_epsilon

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Parameters required for dimensioning some arguments
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
     &, wet_model_levels                                                &
                         ! number of model levels where moisture is held
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, me                                                              &
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
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same row
     &, proc_all_group                                                  &
                       ! Group id for all processors
     &, max_look                                                        &
                          ! max size of look-up arrays for searches
     &, global_row_length                                               &
                          ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of points in a column
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                                                  ! processor on my
                    ! processor-row holding a given value in i direction
     &, g_j_pe(1-halo_j:global_rows      +halo_j)                       &
                                                  ! processor on my
                 ! processor-column holding a given value in j direction
     &, LAM_max_cfl(2)                                                  &
     &, n_Y_arrays                                                      &
                         ! = 1 for global, 3 for LAM
     &, n_Yw_arrays                                                     &
                         ! = 1 for global, 2 for LAM
     &, n_Yd_arrays                                                     &
                         ! = 1 for global, 3 for LAM
     &, n_Ydw_arrays                                                    &
                         ! = 1 for global, 2 for LAM
     &, rimwidth                                                        &
                 ! numbr of rows/cols that use lbcs in their calculation
     &, LENRIM(Nfld_max,NHalo_max)                                      &
     &, LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
     &, LBC_START(4,Nfld_max,NHalo_max)                                 &
     &, CycleNo
      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_lbc_old                                                       &
                         ! false for new lbc treatment
     &, L_new_tdisc      ! if true activate time discretization
                         ! for iterative scheme


      Integer                                                           &
     &  high_order_scheme                                               &
                           ! a code saying which high order scheme to
                           ! use. 1 = tensor tri-cubic lagrange order
                           ! (j,i,k) no other options available at
                           ! present.
     &, monotone_scheme                                                 &
                        ! a code saying which monotone scheme to use.
                        ! 1 = tri-linear order (j,i,k)
                        ! no other options available at present.
     &, Depart_scheme                                                   &
                        ! code saying which departure point scheme to
                        ! use.
     &, Depart_order                                                    &
                        ! for the chosen departure point scheme how
                        ! many iterations/terms to use.
     &, Robert_high_order_scheme                                        &
                                 ! code choosing high order
                          ! interpolation scheme used in Robert routine
     &, Robert_monotone_scheme                                          &
                                 ! code choosing monotone
                          ! interpolation scheme used in Robert routine
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      Logical                                                           &
     &  L_2d_sl_geometry                                                &
                         ! True, then only perform vector co-ordinate
                         !       geometry in 2d.
     &, L_free_slip                                                     &
                         ! True, free-slip lower boundary condition
     &, L_high                                                          &
                       ! True, if high order interpolation required.
     &, L_mono                                                          &
                       ! True, if interpolation required to be monotone.
     &, L_conserv                                                       &
                       ! True, if interpolation to be monotone and
                       !       conservative.
     &, L_mix_ratio                                                     &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                           &
     &,  L_regular                                                      &
                     !  False if variable res. Default .true.
     &, L_interp_depart

      Logical :: L_qwaterload             ! add waterloading terms
      Logical :: L_trivial_trigs ! True if trivial_trigs (Cartesian grid)

      Logical                                                           &
     &  L_Robert_high                                                   &
                      ! True if high order interpolation scheme to be
                      ! used in Robert scheme
     &, L_Robert_mono ! True if monotone interpolation scheme to be
                      ! used in Robert scheme

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, timestep                                                        &
     &, base_phi                                                        &
     &, base_lambda                                                     &
     &, lambda_p_end                                                    &
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, phi_p_end                                                       &
     &, dlambda_p_end                                                   &
     &, dphi_p_end                                                      &
     &, dphi_v_end

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

      Real                                                              &
           ! time-weighting parameters, see WP 154.
     &  alpha_3                                                         &
     &, alpha_4

      Real                                                              &
           ! trigonometric functions
     &  cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sin_theta_latitude (row_length, rows)                           &
     &, tan_theta_latitude (row_length, rows)                           &
     &, cos_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, sec_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, tan_v_latitude (row_length, n_rows)                             &
     &, cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

      Real                                                              &
           ! components of coriolis force.
     &  f1_at_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)        &
     &, f2_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, f3_at_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

      Real, Intent (InOut) ::                                           &
                              ! primary model variables
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, u_np1(1-off_x:row_length+off_x,                                 &
     &        1-off_y:rows+off_y,model_levels)                          &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, v_np1(1-off_x:row_length+off_x,                                 &
     &        1-off_y:n_rows+off_y, model_levels)                       &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels) &
     &, w_np1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        0:model_levels)                                           &
     &, u_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, v_adv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &         model_levels)                                            &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         0:model_levels)                                          &
     &, q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     wet_model_levels)                                            &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_model_levels)                                            &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_model_levels)                                            &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     wet_model_levels)                                            &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_model_levels)                                            &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_model_levels)                                            &
     &, q_np1 (1-off_x:row_length+off_x,                                &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, qcl_np1 (1-off_x:row_length+off_x,                              &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, qcf_np1 (1-off_x:row_length+off_x,                              &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, qcf2_np1 (1-off_x:row_length+off_x,                             &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, qrain_np1 (1-off_x:row_length+off_x,                            &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, qgraup_np1 (1-off_x:row_length+off_x,                           &
     &     1-off_y:rows+off_y, wet_model_levels)                        &
     &, theta_np1 (1-off_x:row_length+off_x,                            &
     &         1-off_y:rows+off_y, model_levels)                        &
     &, theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &          wet_model_levels)                                       &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &          wet_model_levels)                                       &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &          wet_model_levels)                                       &
     &, mix_v_np1 (1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_cl_np1 (1-off_x:row_length+off_x,                           &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_cf_np1 (1-off_x:row_length+off_x,                           &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_cf2_np1 (1-off_x:row_length+off_x,                          &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_rain_np1 (1-off_x:row_length+off_x,                         &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_graup_np1 (1-off_x:row_length+off_x,                        &
     &             1-off_y:rows+off_y, wet_model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)                                           &
     &, eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)

      Real                                                              &
           ! co-ordinates of departure points for w.
     &  depart_lambda_w(row_length, rows, model_levels)                 &
     &, depart_phi_w(row_length, rows, model_levels)                    &
     &, depart_r_w(row_length, rows, model_levels)

      REAL                                                              &
     &  u_lbcs(LENRIM(fld_type_u,halo_type_extended),MODEL_LEVELS)      &
     &, v_lbcs(LENRIM(fld_type_v,halo_type_extended),MODEL_LEVELS)      &
     &, w_lbcs(LENRIM(fld_type_p,halo_type_extended),MODEL_LEVELS)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent INOUT. ie: Input variables changed on output.

! Arguments with Intent OUT. ie: Output variables.

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Real                                                              &
           ! See WP154 and WP162 for definitions.
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &         model_levels)                                            &
     &, R_w (row_length, rows, model_levels)

! Local Variables

      Real                                                              &
     &  thetav (row_length+1, rows+1, model_levels)                     &
     &, moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &          wet_model_levels)                                       &
     &, moist_np1 (1-off_x:row_length+off_x,                            &
     &     1-off_y:rows+off_y, wet_model_levels)   

      Real, Dimension (:,:,:), ALLOCATABLE :: thetav_np1

      Integer                                                           &
     &  i, j, k                                                         &
                  ! loop counters
     &, j1        ! loop bound


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('NI_SL_FULL_WIND',zhook_in,zhook_handle)

      j1 = rows
      If (model_domain == mt_Global) Then
        If (at_extremity(PNorth)) j1 = rows-1
      End If

      IF(L_mix_ratio)then

        moist = 1 + mix_v + mix_cl + mix_cf
        If(L_mcr_qcf2)then
          moist = moist + mix_cf2
        endif
        If(L_mcr_qrain)then
          moist = moist + mix_rain
        endif
        If(L_mcr_qgraup)then
          moist = moist + mix_graup
        endif

        Do k = 1, wet_model_levels
          Do j = 1, j1 + 1
            Do i = 1, row_length + 1
              thetav(i,j,k) = theta(i,j,k) *                            &
     &                       ( 1. +   recip_epsilon * mix_v(i,j,k) ) /  &
     &                                                moist(i,j,k)
            End Do
          End Do
        End Do ! k = 1, wet_model_levels
!
! For iterative scheme, compute thetav^(1), i.e. thetav at tn+1, as well
!
        If ( CycleNo > 1 .AND. L_new_tdisc ) Then

! Use do loops as arrays have different halo sizes
          Do k = 1, wet_model_levels
            Do j = 1-off_y, rows+off_y
              Do i = 1-off_x, row_length+off_x
                moist(i,j,k) = 1. + mix_v_np1(i,j,k) + mix_cl_np1(i,j,k)&
     &                            + mix_cf_np1(i,j,k)
              End Do
            End Do
          End Do
          If ( L_mcr_qcf2 ) Then
            Do k = 1, wet_model_levels
              Do j = 1-off_y, rows+off_y
                Do i = 1-off_x, row_length+off_x
                  moist(i,j,k) = moist(i,j,k) + mix_cf2_np1(i,j,k)
                End Do
              End Do
            End Do
          End If
          If ( L_mcr_qrain ) Then
            Do k = 1, wet_model_levels
              Do j = 1-off_y, rows+off_y
                Do i = 1-off_x, row_length+off_x
                  moist(i,j,k) = moist(i,j,k) + mix_rain_np1(i,j,k)
                End Do
              End Do
            End Do
          End If
          If ( L_mcr_qgraup ) Then
            Do k = 1, wet_model_levels
              Do j = 1-off_y, rows+off_y
                Do i = 1-off_x, row_length+off_x
                  moist(i,j,k) = moist(i,j,k) + mix_graup_np1(i,j,k)
                End Do
              End Do
            End Do
          End If

          Allocate( thetav_np1 (row_length+1, rows+1, model_levels) )

          Do k = 1, wet_model_levels
            Do j = 1, j1 + 1
              Do i = 1, row_length + 1
                thetav_np1(i,j,k) = theta_np1(i,j,k) *                  &
     &                       ( 1. + recip_epsilon*mix_v_np1(i,j,k) ) /  &
     &                                                moist(i,j,k)
              End Do
            End Do
          End Do ! k = 1, wet_model_levels

        Else

          Allocate( thetav_np1 (1,1,1) )

        End If


      else  ! L_mix_ratio=.false.

      ! Include condensate/hydrometeors in theta_v calculation
        If (L_qwaterload) Then
          moist= qcl + qcf
          If (L_mcr_qcf2) Then
            moist = moist + qcf2
          End If
          If (L_mcr_qrain) Then
            moist = moist + qrain
          End If
          If (L_mcr_qgraup) Then
            moist = moist + qgraup
          End If
        Else
          moist = 0.0
        End If

        If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
          Allocate( thetav_np1 (1,1,1) )
          Do k = 1, wet_model_levels
            Do j = 1, j1 + 1
              Do i = 1, row_length + 1
                thetav(i,j,k) = theta(i,j,k) *                          &
     &          ( 1.0 + ((recip_epsilon-1.0)*q(i,j,k)) - moist(i,j,k) ) 
              End Do
            End Do
          End Do ! k = 1, wet_model_levels
        Else
          Allocate( thetav_np1 (row_length+1, rows+1, model_levels) )
! For iterative scheme, compute thetav^(1), i.e. thetav at tn+1, as well
        If (L_qwaterload) Then
          moist_np1= qcl_np1 + qcf_np1
          If (L_mcr_qcf2) Then
            moist_np1 = moist_np1 + qcf2_np1
          End If
          If (L_mcr_qrain) Then
            moist_np1 = moist_np1 + qrain_np1
          End If
          If (L_mcr_qgraup) Then
            moist_np1 = moist_np1 + qgraup_np1
          End If
        Else
          moist_np1 = 0.0
        End If
          Do k = 1, wet_model_levels
            Do j = 1, j1+1
              Do i = 1, row_length+1
                thetav(i,j,k) = theta(i,j,k) *                          &
     &          ( 1.0 + ((recip_epsilon-1.0)*q(i,j,k)) - moist(i,j,k) ) 
                thetav_np1(i,j,k) = theta_np1(i,j,k) *                  &
     &          ( 1.0 + ((recip_epsilon-1.0)*q_np1(i,j,k))              &
     &                                             - moist_np1(i,j,k) ) 
              End Do
            End Do
          End Do ! k = 1, wet_model_levels
        End If

      endif  ! L_mix_ratio

      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
        Do k = 1 + wet_model_levels, model_levels
          Do j = 1, j1 + 1
            Do i = 1, row_length + 1
              thetav(i,j,k) = theta(i,j,k)
            End Do
          End Do
        End Do
      Else
        Do k = 1 + wet_model_levels, model_levels
          Do j = 1, j1 + 1
            Do i = 1, row_length + 1
              thetav(i,j,k) = theta(i,j,k)
              thetav_np1(i,j,k) = theta_np1(i,j,k)
            End Do
          End Do
        End Do
      End If

! DEPENDS ON: sl_full_wind
        Call SL_Full_wind(                                              &
     &                    u, u_np1, v, v_np1, w, w_np1,                 &
     &                    u_adv, v_adv, w_adv,                          &
     &                    thetav, thetav_np1, exner,                    &
     &                    depart_lambda_w, depart_phi_w,                &
     &                    depart_r_w,                                   &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    eta_theta_levels, eta_rho_levels,             &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    sin_theta_latitude, cos_v_latitude,           &
     &                    sec_v_latitude, sin_v_latitude,               &
     &                    tan_theta_latitude, tan_v_latitude,           &
     &                    cos_theta_longitude,                          &
     &                    sin_theta_longitude,                          &
     &                    f1_at_v, f2_at_u, f3_at_u, f3_at_v,           &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    glambda_p, phi_p, glambda_u, phi_v,           &
     &                    gdlambda_p, dphi_p, gdlambda_u, dphi_v,       &
     &                    grecip_dlamp, recip_dphip, grecip_dlamu,      &
     &                    recip_dphiv, wt_lambda_p, wt_phi_p,           &
     &                    wt_lambda_u, wt_phi_v, lambda_p_rm,           &
     &                    lambda_p_rp, lambda_u_rm, lambda_u_rp,        &
     &                    phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,       &
     &                    recip_lambda_p_m, recip_lambda_p_0,           &
     &                    recip_lambda_p_p, recip_lambda_p_p2,          &
     &                    recip_lambda_u_m, recip_lambda_u_0,           &
     &                    recip_lambda_u_p, recip_lambda_u_p2,          &
     &                    recip_phi_p_m, recip_phi_p_0,                 &
     &                    recip_phi_p_p, recip_phi_p_p2,                &
     &                    recip_phi_v_m, recip_phi_v_0,                 &
     &                    recip_phi_v_p, recip_phi_v_p2,                &
     &                    base_lambda, base_phi,                        &
     &                    lambda_p_end, phi_p_end, dlambda_p_end,       &
     &                    dphi_p_end, dphi_v_end,                       &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    alpha_3, alpha_4, LAM_max_cfl,                &
     &                    n_Y_arrays, n_Yw_arrays,                      &
     &                    n_Yd_arrays, n_Ydw_arrays,                    &
     &                    u_lbcs,v_lbcs,w_lbcs,                         &
     &                    LENRIM,LBC_SIZE,LBC_START,RIMWIDTH,           &
     &                    model_domain,                                 &
     &                    row_length, rows, n_rows, model_levels,       &
     &                    Depart_scheme, Depart_order,                  &
     &                    high_order_scheme, monotone_scheme,           &
     &                    L_trivial_trigs,                              &
     &                    L_high, L_mono, L_conserv,                    &
     &                    L_Robert_high, L_Robert_mono,                 &
     &                    Robert_high_order_scheme,                     &
     &                    Robert_monotone_scheme,                       &
     &                    first_constant_r_rho_level,                   &
     &                    check_bottom_levels,                          &
     &                    interp_vertical_search_tol,                   &
     &                    r_at_u, r_at_v,                               &
     &                    me, n_proc, n_procx, n_procy,                 &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    global_row_length, global_rows,               &
     &                    l_datastart, at_extremity, g_i_pe, g_j_pe,    &
     &                    l_2dcomm, size_2dcomm, group_2dcomm,          &
     &                    proc_row_group, proc_col_group,               &
     &                    proc_all_group,                               &
     &                    L_sl_halo_reprod, L_2d_sl_geometry,           &
     &                    L_free_slip, L_regular, L_interp_depart,      &
     &                    L_lbc_old, L_new_tdisc, CycleNo,              &
     &                    R_u, R_v, R_w, Error_code)

      Deallocate ( thetav_np1 )
! End of routine.
      IF (lhook) CALL dr_hook('NI_SL_FULL_WIND',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_SL_Full_Wind

