! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Thermo

      Subroutine SL_Thermo(                                             &
     &                     moisture_array_size, timestep,               &
     &                     theta, vap, cl, cf, cf2, rain, graup,        &
     &                     cf_pc2, cfl_pc2, cff_pc2,                    &
     &                     e_trb, tsq, qsq, cov,                        &
     &                     vap_star, cl_star, cf_star,                  &
     &                     cf2_star, rain_star, graup_star,             &
     &                     cf_pc2_star, cfl_pc2_star, cff_pc2_star,     &
     &                     exner_star, theta_star, theta_np1,           &
     &                     u_adv, v_adv, w_adv, w,                      &
     &                     eta_rho_levels, eta_theta_levels,            &
     &                     r_rho_levels, r_theta_levels,                &
     &                     exner_theta_levels, mono_mass,               &
     &                     row_length, rows, n_rows,                    &
     &                     model_levels, wet_model_levels, bl_levels,   &
     &                     alpha_2, check_bottom_levels,                &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     delta_lambda, delta_phi,                     &
     &                     glambda_p, phi_p, glambda_u, phi_v,          &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamp, recip_dphip,                   &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,&
     &                     lambda_p_rm, lambda_p_rp,                    &
     &                     lambda_u_rm, lambda_u_rp,                    &
     &                     phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,      &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_lambda_u_m, recip_lambda_u_0,          &
     &                     recip_lambda_u_p, recip_lambda_u_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     recip_phi_v_m, recip_phi_v_0,                &
     &                     recip_phi_v_p, recip_phi_v_p2,               &
     &                     Base_lambda, base_phi,                       &
     &                     lambda_p_end, phi_p_end,                     &
     &                     dlambda_p_end, dphi_p_end, dphi_v_end,       &
     &                     recip_dlam, recip_dphi, max_look,            &
     &                     look_lam, look_phi, halo_lam, halo_phi,      &
     &                           FV_cos_theta_latitude,                 &
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
! roar's bit
     &                           me, n_proc, n_procx, n_procy,          &
     &                           l_datastart, g_i_pe, g_j_pe,           &
     &                           l_2dcomm, size_2dcomm,                 &
     &                           group_2dcomm, max_comm_size,           &
     &                           at_extremity,                          &
     &                           global_row_length, global_rows,        &
     &                           proc_row_group, proc_col_group,        &
     &                           proc_all_group,                        &
     &                           off_x, off_y, halo_i, halo_j,          &
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
     &                           L_2d_sl_geometry, L_sl_halo_reprod,    &
     &                           L_regular, L_fint_theta,               &
     &                           L_mcr_cf2, L_mcr_rain, L_mcr_graup,    &
     &                           L_lbc_old, L_new_tdisc, CycleNo,       &
     &                           depart_r_theta,                        &
     &                           depart_lambda_w, depart_phi_w,         &
     &                           depart_r_w, Error_Code)

! Purpose:
!          Performs semi-Lagrangian advection of theta, moisture,
!          cl and cf. Also outputs the departure point for use in 
!          calculating w later.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE mym_option_mod, ONLY:                                         &
                            bdy_tke, l_adv_turb_field, tke_levels,      &
                            high_order_scheme_adv_turb,                 &
                            monotone_scheme_adv_turb,                   &
                            l_high_adv_turb,                            &
                            l_mono_adv_turb,                            &
                            l_conserv_adv_turb

      USE cloud_inputs_mod, ONLY: l_fixbug_pc2_mixph
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
     &, l_datastart(3)                                                  &
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
      INTEGER :: size_int      ! error check size for comms on demand
      INTEGER :: size_int_mult ! error check size for comms on demand
      INTEGER :: max_comm_size ! error check size for comms on demand

      Integer                                                           &
     &  first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, CycleNo

      Integer, Parameter ::                                             &
     &  etrb_array_size = 4

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_new_tdisc                                                     &
     &, L_regular                                                       &
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
     &, L_pc2          ! True, if clouds need advecting for PC2 scheme

      Logical                                                           &
     &  L_ritchie_high                                                  &
                       ! True if high order interpolation scheme to be
                      ! used in ritchie scheme
     &, L_ritchie_mono                                                  &
                       ! True if monotone interpolation scheme to be
                      ! used in ritchie scheme
     &, L_2d_sl_geometry                                                &
     &, L_mcr_cf2                                                       &
                       ! True if second ice variable is in use
     &, L_mcr_rain                                                      &
                       ! True if prognostic rain is in use
     &, L_mcr_graup    ! True if prognostic graupel is in use

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

      Real                                                              &
     &  theta_np1(1-off_x:row_length+off_x,                             &
     &             1-off_y:rows+off_y, model_levels)

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

      Real, Intent (InOut) ::                                           &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, vap (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &       wet_model_levels)                                          &
     &, cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &       wet_model_levels)                                          &
     &, cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &        wet_model_levels)                                         &
     &, rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &         wet_model_levels)                                        &
     &, graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
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
     &     0:model_levels)

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Real                                                              &
     &  exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, mono_mass(row_length, rows, model_levels)

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
     &  depart_lambda_w (row_length, rows, model_levels)                &
                                                        ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_phi_w (row_length, rows, model_levels)                   &
                                                      ! Phi Co-ordinate
                                                     ! of co-ordinate of
                                                      ! departure point.
     &, depart_r_w (row_length, rows, model_levels)                     &
                                                      ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_r_theta (row_length, rows, model_levels)                 &
                                                      ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.
     &, exner_star(1-off_x:row_length+off_x,                            &
                                                      ! Departure value
     &             1-off_y:rows+off_y, model_levels)  ! of exner.

! Star variables hold physics increments on input, latest values on
! output
      Real, Intent (InOut) ::                                           &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, model_levels)                    &
     &, vap_star(1-off_x:row_length+off_x,                              &
     &         1-off_y:rows+off_y, wet_model_levels)                    &
     &, cl_star(1-off_x:row_length+off_x,                               &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, cf_star(1-off_x:row_length+off_x,                               &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, cf2_star(1-off_x:row_length+off_x,                              &
     &            1-off_y:rows+off_y, wet_model_levels)                 &
     &, rain_star(1-off_x:row_length+off_x,                             &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, graup_star(1-off_x:row_length+off_x,                            &
     &              1-off_y:rows+off_y, wet_model_levels)               &
     &, cf_pc2_star (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)              &
     &, cfl_pc2_star(1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)              &
     &, cff_pc2_star(1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y, wet_model_levels)

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                         ! Loop indices
     &, min_k, temp                                                     &
     &, n_inputs                                                        &
     &, n_ext_fields    ! number of ext_data fields required, 1 for
                        ! global model, 2 for LAM.

      Integer                                                           &
     &  type           ! a code saying what points the grid points are
                       ! on. For theta this code is 4.

      Logical                                                           &
     &  L_vector       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.

      Real                                                              &
     &  dummy, dummy3D(1,1,1)

      Real                                                              &
     &  del_theta_NI, del_theta_3D
!                      ! temporary variables used in corrected
!                      ! intelligent limiter code.

      Logical                                                           &
     &  L_do_halos                                                      &
                          ! update the halos?
     &, L_do_boundaries   ! update the boundaries?

! arrays

      Real                                                              &
     &  depart_lambda (row_length, rows, model_levels)                  &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_phi (row_length, rows, model_levels)                     &
                                                    ! Phi Co-ordinate of
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_r (row_length, rows, model_levels)     ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi    (row_length, rows, model_levels)

      Real                                                              &
     &  work(row_length, rows, model_levels)

      Real                                                              &
     &  theta_sl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &            model_levels)

      Real                                                              &
     &  super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,    &
     &           wet_model_levels,moisture_array_size)                  &
     &, data_out_super(row_length, rows,                                &
     &           wet_model_levels,moisture_array_size)                  &
     &, super_array_etrb(1-halo_i:row_length+halo_i,                    &
     &                   1-halo_j:rows+halo_j,                          &
     &                   tke_levels - 1,etrb_array_size)                &
     &, data_out_super_etrb(row_length, rows,                           &
     &                   tke_levels - 1,etrb_array_size)
      integer count
      Real, Dimension (:,:,:), Allocatable :: work2
      Real, Dimension (:,:,:), Allocatable :: theta_sl2

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Find Departure Point for w points.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_THERMO',zhook_in,zhook_handle)
      type = 3  ! w points.
      Error_Code=0

! DEPENDS ON: departure_point
      Call Departure_Point(                                             &
     &                     type, timestep, u_adv, v_adv, w_adv,         &
     &                     eta_rho_levels, eta_theta_levels,            &
     &                     r_rho_levels, r_theta_levels,                &
     &                     r_at_u, r_at_v,                              &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     delta_lambda, delta_phi,                     &
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
     &                     cos_theta_latitude, sec_theta_latitude,      &
     &                     sin_theta_latitude, tan_theta_latitude,      &
     &                     cos_v_latitude, sec_v_latitude,              &
     &                     sin_v_latitude, tan_v_latitude,              &
     &                     sin_theta_longitude,                         &
     &                     cos_theta_longitude,                         &
     &                     model_domain, Depart_scheme,                 &
     &                     Depart_order, L_ritchie_high, L_ritchie_mono,&
     &                     ritchie_high_order_scheme,                   &
     &                     ritchie_monotone_scheme,                     &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     L_2d_sl_geometry, L_regular,                 &
     &                     first_constant_r_rho_level,                  &
     &                     LAM_max_cfl, rows,                           &

! Roar's bit
     &                     me, n_proc, n_procx, n_procy,                &
     &                     off_x, off_y, halo_i, halo_j, l_datastart,   &
     &                     global_row_length, global_rows, g_i_pe,      &
     &                     g_j_pe, l_2dcomm, size_2dcomm,               &
     &                     group_2dcomm,                                &
     &                     proc_row_group,proc_col_group,               &
     &                     proc_all_group, at_extremity,                &
     &                     L_sl_halo_reprod,                            &

     &                     depart_lambda_w, depart_phi_w, depart_r_w)

! ----------------------------------------------------------------------
! Section 2.0  Calculate estimate of theta at new time level.
!              Copy w departure point locations.
!              Limiting trajectory near ground so that lowest value
!              is at the bottom data level is performed inside
!              calc_non_int_sl_theta.
! ----------------------------------------------------------------------

      If (Error_Code  ==  0 ) Then
! copy departure point locations.

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              depart_lambda(i,j,k) = depart_lambda_w(i,j,k)
              depart_phi(i,j,k) = depart_phi_w(i,j,k)
              depart_r(i,j,k) = depart_r_w(i,j,k)
            End Do
          End Do
        End Do

        If ( thmono_levels > 0 .and. ( .not. L_fint_theta ) ) Then
          ALLOCATE ( theta_sl2(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j, model_levels) )
          ALLOCATE ( work2(row_length, rows, model_levels) )
        End If ! thmono_levels > 0 .and. ( .not. L_fint_theta )

! copy physics theta_increment into theta_sl. Copy exner into exner_sl
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_sl(i,j,k) = theta_star(i,j,k)
            End Do
          End Do
          If ( thmono_levels > 0 .and. ( .not. L_fint_theta ) ) Then
            Do j = 1, rows
              Do i = 1, row_length
                theta_sl2(i,j,k) = theta(i,j,k) + theta_star(i,j,k)
              endDo
            endDo
          End If ! thmono_levels > 0
        End Do

! ----------------------------------------------------------------------
! Section 2.1  Perform non-interpolating in the vertical semi-Lagrangian
!              advection of theta.
! ----------------------------------------------------------------------

        If (model_domain  ==  mt_LAM) Then
           n_ext_fields = 2
        Else
           n_ext_fields = 1
        End If

! DEPENDS ON: calc_index
        Call Calc_Index(                                                &
     &                      row_length, rows, model_levels,             &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular, depart_lambda, depart_phi,       &
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

!  L_mono_theta = .false. in Calc_non_int_sl_theta
!  so advection of theta_phys incs can be monotone if desired
!  Note (MD): Even when fully interpolating theta advection is requested
!             the call to the following subroutine is needed to shift 
!             these departure points which are below the ground.  
! DEPENDS ON: calc_non_int_sl_theta
        Call calc_non_int_sl_theta(                                     &
     &                          w, theta, theta_lbcs,                   &
     &                          r_theta_levels, r_rho_levels,           &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          row_length, rows, model_levels,         &
     &                          rimwidth, rimweights,                   &
     &                          lenrim, lbc_size, lbc_start,            &
     &                          n_ext_fields,                           &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          model_domain, timestep, alpha_2,        &
     &                          high_order_scheme_theta,                &
     &                          monotone_scheme_theta,                  &
     &                          L_high_theta, .false.,                  &
     &                          L_sl_halo_reprod, L_lbc_old,            &
     &                          L_regular, L_new_tdisc, CycleNo,        &
     &                          depart_r, me, n_proc, n_procx, n_procy, &
     &                          off_x, off_y, halo_i, halo_j,           &
     &                          l_datastart,  at_extremity,             &
     &                          global_row_length, global_rows,         &
     &                  proc_row_group, proc_col_group, proc_all_group, &
     &                          g_i_pe, g_j_pe, l_2dcomm,               &
     &                        size_2dcomm, group_2dcomm, max_comm_size, &
     &                          theta_star, theta_np1,                  &
     &                          Error_Code)

!     store theta departure points for use in tracer code
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              depart_r_theta(i,j,k) = depart_r(i,j,k)
            End Do
          End Do
        End Do

      End If  ! Error_Code == 0

! ----------------------------------------------------------------------
! Section 2.2  Interpolate physics increment
! ----------------------------------------------------------------------

      If ( Error_Code == 0 ) Then

        If ( .not. L_fint_theta ) Then

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &          theta_sl, row_length, rows,                             &
     &          model_levels, halo_i, halo_j, fld_type_p, .false. )

          If ( model_domain == mt_LAM ) Then

           If ( L_lbc_old ) Then

            L_do_halos=.TRUE.
            L_do_boundaries=.TRUE.

! DEPENDS ON: zero_lateral_boundaries
            CALL ZERO_LATERAL_BOUNDARIES(                               &
     &                 ROW_LENGTH,ROWS,halo_i,halo_j,MODEL_LEVELS,      &
     &                 fld_type_p,theta_sl,                             &
     &                 1, AT_EXTREMITY,                                 &
     &                 L_do_boundaries,L_do_halos)

           End If ! L_lbc_old

          End If ! model_domain == mt_LAM

          If ( thmono_levels > 0 ) Then

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                     theta_sl2, row_length, rows, model_levels,   &
     &                     halo_i, halo_j, fld_type_p, .false. )

            If ( model_domain == mt_LAM ) Then

             If ( L_lbc_old ) Then

! DEPENDS ON: set_lateral_boundaries
              CALL SET_LATERAL_BOUNDARIES(                              &
     &                                    ROW_LENGTH, ROWS,             &
     &                                    halo_i, halo_j, MODEL_LEVELS, &
     &                                    fld_type_p,                   &
     &                                   theta_sl2(1-halo_i,1-halo_j,1),&
     &                                    LENRIM, LBC_SIZE, LBC_START,  &
     &                                    HALO_I, HALO_J, THETA_LBCS,   &
     &                                    RIMWIDTH, RIMWIDTH,           &
     &                                    RIMWEIGHTS, AT_EXTREMITY,     &
     &                                    L_do_boundaries, L_do_halos)
            End If ! L_lbc_old

            End If ! model_domain == mt_LAM

          End If ! thmono_levels > 0

        Else ! L_fint_theta = true :

          Do k = 1, model_levels 
            Do j = 1, rows 
              Do i = 1, row_length
! theta_sl = theta_n + theta_inc after phys1 = theta after phys1
                theta_sl(i,j,k) = theta_sl(i,j,k) + theta(i,j,k)
              End Do
            End Do
          End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &        theta_sl, row_length, rows,                               &
     &        model_levels, halo_i, halo_j, fld_type_p, .false. )

          If ( model_domain == mt_LAM ) Then

           If ( L_lbc_old ) Then

! Set data on edge processor edge haloes from lateral boundary data 
            L_do_halos = .TRUE.           
            L_do_boundaries = .FALSE.      

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &                               ROW_LENGTH, ROWS, halo_i, halo_j,  &
     &                               MODEL_LEVELS, fld_type_p, theta_sl,&
     &                               LENRIM, LBC_SIZE, LBC_START,       &
     &                               halo_i, halo_j, theta_lbcs,        &
     &                               RIMWIDTH, RIMWIDTH,                &
     &                               RIMWEIGHTS, AT_EXTREMITY,          &
     &                               L_do_boundaries, L_do_halos )

           End If !  L_lbc_old

          End If ! model_domain == mt_LAM

        End If ! .not. L_fint_theta

        L_vector = .false.
        n_inputs = 1
        if(l_2dcomm)then
          size_int = row_length * rows * model_levels
        else
          size_int = model_levels * global_row_length
        endif

! DEPENDS ON: interpolation
        Call Interpolation(                                             &
     &                     theta_sl, dummy3D, dummy3D,                  &
     &                     eta_theta_levels(1),                         &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     dummy,                                       &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     type, n_inputs, check_bottom_levels,         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, model_levels,              &
     &                     rows,                                        &
     &                     row_length, rows, model_levels,              &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_p_rm, lambda_p_rp,                    &
     &                     phi_p_rm, phi_p_rp,                          &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme_theta,                     &
     &                     monotone_scheme_theta,                       &
     &                     FV_cos_theta_latitude, L_regular,            &
     &                     L_vector, model_domain, L_high_theta,        &
     &                     L_mono_theta, .false.,                       &
     &                     depart_r, depart_lambda, depart_phi,         &

! Roar's bit
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     g_j_pe, l_2dcomm, size_2dcomm,               &
     &                     group_2dcomm, size_int, proc_all_group,      &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     0, 0, L_sl_halo_reprod,                      &
     &                     off_x,off_y,                                 &

     &                     work, dummy3D, dummy3D, Error_Code)

! Apply monotone limiter only when a non-interpolating in the vertical
! scheme is used        
        If ( .not. L_fint_theta ) Then

          If ( thmono_levels > 0 ) Then
            if(l_2dcomm)then
              size_int = row_length * rows * thmono_levels
            else
              size_int = thmono_levels * global_row_length
            endif
!  Do a full-interpolation of theta
! DEPENDS ON: interpolation
            Call Interpolation(                                         &
     &                       theta_sl2, dummy3D, dummy3D,               &
     &                       eta_theta_levels(1),                       &
     &                       r_theta_levels(1-halo_i,1-halo_j,1),       &
     &                       dummy,                                     &
     &                       r_theta_levels(1-halo_i,1-halo_j,1),       &
     &                       type, n_inputs, check_bottom_levels,       &
     &                       interp_vertical_search_tol,                &
     &                       first_constant_r_rho_level,                &
     &                       row_length, rows, model_levels,            &
     &                       rows,                                      &
     &                       row_length, rows, thmono_levels,           &
     &                       delta_lambda, delta_phi,                   &
     &                       base_lambda, base_phi,                     &
     &                       glambda_u, phi_v,                          &
     &                       gdlambda_p, dphi_p, gdlambda_u, dphi_v,    &
     &                       grecip_dlamu, recip_dphiv,                 &
     &                       lambda_p_rm, lambda_p_rp,                  &
     &                       phi_p_rm, phi_p_rp,                        &
     &                       recip_lambda_p_m, recip_lambda_p_0,        &
     &                       recip_lambda_p_p, recip_lambda_p_p2,       &
     &                       recip_phi_p_m, recip_phi_p_0,              &
     &                       recip_phi_p_p, recip_phi_p_p2,             &
     &                       i_out, j_out,                              &
     &                       weight_lambda, weight_phi,                 &
     &                       high_order_scheme_theta,                   &
     &                       monotone_scheme_theta,                     &
     &                       FV_cos_theta_latitude, L_regular,          &
     &                       L_vector, model_domain, L_high_theta,      &
     &                       .true. , .false.,                          &
     &                       depart_r, depart_lambda, depart_phi,       &
     &                       me, n_proc, n_procx, n_procy,              &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                       l_datastart, at_extremity, g_i_pe,         &
     &                       g_j_pe, l_2dcomm, size_2dcomm,             &
     &                       group_2dcomm, size_int, proc_all_group,    &
     &                       proc_row_group, proc_col_group, 1,         &
     &                       0, 0, L_sl_halo_reprod,                    &
     &                       off_x,off_y,                               &
     &                       work2, dummy3D, dummy3D, Error_Code)

          End If ! If  thmono_levels > 0
! add on interpolated physics increment
! Limit non-interpolated value by monotone value at lowest thmono_levels
! vn7.6 - corrected code introduced using L_thmono_fixed=true
!
          If ( L_thmono_fixed ) Then          
            Do k = 1, thmono_levels
              Do j = 1, rows
                Do i = 1, row_length
                  work(i,j,k) = theta_star(i,j,k) + work(i,j,k)
                  del_theta_NI = work(i,j,k)  - theta(i,j,k)                 
                  del_theta_3D = work2(i,j,k) - theta(i,j,k)
                  If ( del_theta_NI*del_theta_3D .lt. 0.0) Then
!                   ! Nigel not keen on line below !
                    theta_star(i,j,k) = theta(i,j,k) 
                  Else If (abs(del_theta_NI) .gt. abs(del_theta_3D))    & 
     &                                                              Then
                    theta_star(i,j,k) = work2(i,j,k)
                  Else
                    theta_star(i,j,k) = work(i,j,k)
                  End If
                End Do
              End Do
            End Do !  k = 1, thmono_levels
          Else ! L_thmono_fixed = false :          
            Do k = 1, thmono_levels
              Do j = 1, rows
                Do i = 1, row_length
                  work(i,j,k) = theta_star(i,j,k) + work(i,j,k)
!  Limit non-interpolated value by monotone fully-interpolated value
!  NB could do this the other way round
                  If( abs(work(i,j,k) - theta(i,j,k)) <=                &
     &                abs(work2(i,j,k) - theta(i,j,k)) )Then
                    theta_star(i,j,k) = work(i,j,k)
                  Else  ! monotone increment is smaller
                    theta_star(i,j,k) = work2(i,j,k)
                  End If
                End Do
              End Do
            End Do !  k = 1, thmono_levels
          End If  !  L_thmono_fixed

          Do k = thmono_levels+1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = theta_star(i,j,k) + work(i,j,k)
              End Do
            End Do
          End Do

          If ( thmono_levels > 0 ) Then
            DEALLOCATE ( theta_sl2 )
            DEALLOCATE ( work2 )
          End If ! thmono_levels > 0

        Else ! L_fint_theta = true :

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = work(i,j,k)
              End Do
            End Do
          End Do

        End If ! L_fint_theta

      End If ! Error_code = 0

! ----------------------------------------------------------------------
! Section 3.0  Calculate vap, cl and cf at departure point.
!              Trajectory limit near ground already done in section 2.
!              Trajectory limit at top required if
!              model_levels > wet_model_levels.
! ----------------------------------------------------------------------

      If (Error_Code  ==  0 ) Then

! Perform trajectory limiter near top of model if
! model_levels > wet_model_levels
        If (model_levels  >   wet_model_levels ) Then

! is height field constant on this level,
! if so much cheaper code can be used.
          If (wet_model_levels  >=  first_constant_r_rho_level) Then

! assume once one level has no data above top then no lower level
! has either.
            min_k = max(1,wet_model_levels-interp_vertical_search_tol)
            k = wet_model_levels
            temp =1
            Do while (k  >=  min_k .and. temp  >   0 )
              temp = 0
              Do j = 1, rows
                Do i = 1, row_length
                  If (depart_r(i,j,k)  >                                &
     &                r_theta_levels(1,1,wet_model_levels) ) Then
                    depart_r(i,j,k) =                                   &
     &                         r_theta_levels(1,1,wet_model_levels)
                    temp = temp + 1
                  End If
                End Do
              End Do
              k = k - 1
            End Do

          Else

! calculate trajectory limit at departure point

! assume once one level has no data above top then no lower level
! has either.
            min_k = max(1,wet_model_levels-interp_vertical_search_tol)
            k = wet_model_levels
            temp=1
            Do while (k  >=  min_k .and. temp  >   0 )
              temp = 0

! DEPENDS ON: bi_linear_h
              Call bi_linear_h(                                         &
     &                      r_theta_levels(1-halo_i,1-halo_j,           &
     &                                     wet_model_levels),           &
     &                      depart_lambda(1,1,k),                       &
     &                      depart_phi(1,1,k),                          &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows,                           &
     &                      i_out(1,1,k), j_out(1,1,k),                 &
     &                      weight_lambda(1,1,k), weight_phi(1,1,k),    &
     &                      model_domain, me, n_procx, n_procy,n_proc,  &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      1, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      work)

              Do j = 1, rows
                Do i = 1, row_length

                  If (depart_r(i,j,k) >   work(i,j,1) ) Then
! move trajectory down to highest level of data.
                    depart_r(i,j,k) = work(i,j,1)
                    temp = temp + 1
                  End If

                End Do
              End Do
              call gc_isum (1,n_proc,error_code,temp)
              k = k - 1
            End Do

          End If

        End If

! ----------------------------------------------------------------------
! Section 3.1  Create new variable for advection, vap(n) + vap_inc(phys)
! ----------------------------------------------------------------------

        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_array(i,j,k,1) = vap(i,j,k) + vap_star(i,j,k)
              super_array(i,j,k,2) = cl(i,j,k) + cl_star(i,j,k)
              super_array(i,j,k,3) = cf(i,j,k) + cf_star(i,j,k)
            End Do
          End Do
        End Do
        count=3

        IF (L_pc2) THEN
          IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
            DO k = 1, wet_model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  super_array(i,j,k,count+1) = cfl_pc2(i,j,k) +         &
                                               cfl_pc2_star (i,j,k) +   &
                                               cff_pc2(i,j,k) +         &
                                               cff_pc2_star (i,j,k) -   &
                                               cf_pc2(i,j,k) -          &
                                               cf_pc2_star (i,j,k)
                  super_array(i,j,k,count+2) = cfl_pc2(i,j,k) +         &
                                               cfl_pc2_star (i,j,k)
                  super_array(i,j,k,count+3) = cff_pc2(i,j,k) +         &
                                               cff_pc2_star (i,j,k)
                  super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
                END DO
              END DO
            END DO
          ELSE
! original method, advect the bulk cloud fraction
            DO k = 1, wet_model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  super_array(i,j,k,count+1) = cf_pc2(i,j,k) +          &
                                               cf_pc2_star (i,j,k)
                  super_array(i,j,k,count+2) = cfl_pc2(i,j,k) +         &
                                               cfl_pc2_star (i,j,k)
                  super_array(i,j,k,count+3) = cff_pc2(i,j,k) +         &
                                               cff_pc2_star (i,j,k)
                  super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
                END DO
              END DO
            END DO
          END IF ! l_fixbug_pc2_mixph
          count=count+4
        END IF ! L_pc2

        If (L_mcr_cf2) Then   ! Second cloud ice variable in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
              super_array(i,j,k,count+1) = cf2(i,j,k) + cf2_star(i,j,k)
             End Do
            End Do
          End Do
          count=count+1
        End If

        If (L_mcr_rain) Then   ! Prognostic rain variable in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_array(i,j,k,count+1) = rain(i,j,k) +              &
     &                                       rain_star(i,j,k)
             End Do
            End Do
          End Do
          count=count+1
        End If

        If (L_mcr_graup) Then  ! Prognostic graupel in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
              super_array(i,j,k,count+1) = graup(i,j,k) +               &
     &                                     graup_star(i,j,k)
             End Do
            End Do
          End Do
        End If

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &       super_array, row_length, rows,                             &
     &       moisture_array_size*wet_model_levels,                      &
     &       halo_i, halo_j, fld_type_p, .false. )


       If ( (model_domain == mt_LAM .and. L_lbc_old) .or.               &
     &       model_domain == mt_cyclic_LAM ) Then
! overwrite X_sl with X values on LAM boundaries
          If (at_extremity(PSouth)) Then
            Do k = 1, wet_model_levels
              Do j = 1-halo_j, 1
                Do i = 1-halo_i, row_length + halo_i
                  super_array(i,j,k,1) = vap(i,j,k)
                  super_array(i,j,k,2) = cl(i,j,k)
                  super_array(i,j,k,3) = cf(i,j,k)
                End Do
              End Do
            End Do
            count=3

            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, 1
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl_pc2(i,j,k) +     &
                                                   cff_pc2(i,j,k) -     &
                                                   cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, 1
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              END IF ! l_fixbug_pc2_mixph
              count=count+4
            endif !l_pc2

            if(l_mcr_cf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, 1
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_cf2

            if(l_mcr_rain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, 1
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_rain

            if(l_mcr_graup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, 1
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_graup

          endif  ! at_extremity PSouth

          If (at_extremity(PNorth)) Then
            Do k = 1, wet_model_levels
              Do j = rows, rows+halo_j
                Do i = 1-halo_i, row_length + halo_i
                  super_array(i,j,k,1) = vap(i,j,k)
                  super_array(i,j,k,2) = cl(i,j,k)
                  super_array(i,j,k,3) = cf(i,j,k)
                End Do
              End Do
            End Do
            count=3

            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = rows, rows+halo_j
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl_pc2(i,j,k) +     &
                                                   cff_pc2(i,j,k) -     &
                                                   cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO 
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = rows, rows+halo_j
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              END IF ! l_fixbug_pc2_mixph
              count=count+4
            endif !l_pc2

            if(l_mcr_cf2)then
              Do k = 1, wet_model_levels
                Do j = rows, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_cf2

            if(l_mcr_rain)then
              Do k = 1, wet_model_levels
                Do j = rows, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_rain

            if(l_mcr_graup)then
              Do k = 1, wet_model_levels
                Do j = rows, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_graup

          endif  ! at_extremity PNorth

        End If  ! model_domain (mt_lam/old lbcs or mt_cyclic lam)

        If ( model_domain == mt_LAM ) Then
        
         If (L_lbc_old) Then

          If (at_extremity(PWest)) Then
            Do k = 1, wet_model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = 1-halo_i, 1
                  super_array(i,j,k,1) = vap(i,j,k)
                  super_array(i,j,k,2) = cl(i,j,k)
                  super_array(i,j,k,3) = cf(i,j,k)
                End Do
              End Do
            End Do
            count=3

            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = 1-halo_i, 1
                      super_array(i,j,k,count+1) = cfl_pc2(i,j,k) +     &
                                                   cff_pc2(i,j,k) -     &
                                                   cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO                  
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = 1-halo_i, 1
                      super_array(i,j,k,count+1) = cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              END IF
              count=count+4
            END IF !l_pc2

            if(l_mcr_cf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, 1
                    super_array(i,j,k,count+1) = cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_cf2

            if(l_mcr_rain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, 1
                    super_array(i,j,k,count+1) = rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_rain

            if(l_mcr_graup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, 1
                    super_array(i,j,k,count+1) = graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_graup

          endif  ! at_extremity PWest

          If (at_extremity(PEast)) Then

            Do k = 1, wet_model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = row_length, row_length + halo_i
                  super_array(i,j,k,1) = vap(i,j,k)
                  super_array(i,j,k,2) = cl(i,j,k)
                  super_array(i,j,k,3) = cf(i,j,k)
                End Do
              End Do
            End Do
            count=3

            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = row_length, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl_pc2(i,j,k) +     &
                                                   cff_pc2(i,j,k) -     &
                                                   cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = row_length, row_length + halo_i
                      super_array(i,j,k,count+1) = cf_pc2(i,j,k)
                      super_array(i,j,k,count+2) = cfl_pc2(i,j,k)
                      super_array(i,j,k,count+3) = cff_pc2(i,j,k)
                    END DO
                  END DO
                END DO
              END IF !l_fixbug_pc2_mixph
              count=count+4
            END IF !l_pc2

            if(l_mcr_cf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length, row_length + halo_i
                    super_array(i,j,k,count+1) = cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_cf2

            if(l_mcr_rain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length, row_length + halo_i
                    super_array(i,j,k,count+1) = rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_rain

            if(l_mcr_graup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length, row_length + halo_i
                    super_array(i,j,k,count+1) = graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_graup

          Endif  ! at_extremity PEast

         Else   !  new_lbcs

! DEPENDS ON: fill_external_halos
            call fill_external_halos(                                   &
     &                               super_array, row_length, rows,     &
     &                             moisture_array_size*wet_model_levels,&
     &                               halo_i, halo_j)

         Endif  !  L_lbc_old

        Endif   ! model_domain == mt_lam

          If (l_pc2) then
! DEPENDS ON: fill_external_halos
            call fill_external_halos(                                   &
     &                               super_array(1-halo_i,1-halo_j,1,7),&
     &                               row_length, rows,                  &
     &                               wet_model_levels,halo_i,halo_j)
          Endif  !  l_pc2

! ----------------------------------------------------------------------
! Section 3.2  Call interpolation routine.
! ----------------------------------------------------------------------

        L_vector = .false.
        if(l_2dcomm)then
          size_int_mult = row_length * rows * wet_model_levels
        else
          size_int_mult = wet_model_levels * global_row_length
        endif

! DEPENDS ON: interpolation_multi
        Call Interpolation_multi(                                       &
     &                     super_array,                                 &
     &                     eta_theta_levels(1),                         &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     mono_mass,                                   &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     moisture_array_size,                         &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, wet_model_levels,          &
     &                     rows,                                        &
     &                     row_length, rows, wet_model_levels,          &
     &                     delta_lambda, delta_phi,                     &
     &                     gdlambda_u, dphi_v,                          &
     &                     lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,&
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme_moist,                     &
     &                     monotone_scheme_moist,                       &
     &                     FV_cos_theta_latitude, L_regular,            &
     &                     model_domain, L_high_moist,                  &
     &                     L_mono_moist, L_conserv_moist,               &
     &                     depart_r, depart_lambda, depart_phi,         &

! Roar's bit
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j, global_row_length,           &
     &                     global_rows, row_length, rows, n_rows,       &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     g_j_pe, l_2dcomm, size_2dcomm,               &
     &                     group_2dcomm, size_int_mult, proc_all_group, &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     0, 0, L_sl_halo_reprod,                      &
     &                     off_x, off_y,                                &
     &                     data_out_super, Error_Code )
!
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              vap_star(i,j,k) = data_out_super(i,j,k,1)
              cl_star(i,j,k) = data_out_super(i,j,k,2)
              cf_star(i,j,k) = data_out_super(i,j,k,3)
            End Do
          End Do
        End Do
        count=3

        IF (l_pc2) THEN
          IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
            DO k = 1, wet_model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  cf_pc2_star(i,j,k) = data_out_super(i,j,k,count+2)    &
                                      +data_out_super(i,j,k,count+3)    &
                                      -data_out_super(i,j,k,count+1)
                  cfl_pc2_star(i,j,k) = data_out_super(i,j,k,count+2)
                  cff_pc2_star(i,j,k) = data_out_super(i,j,k,count+3)
                  exner_star(i,j,k)   = data_out_super(i,j,k,count+4)
                END DO
              END DO
            END DO
          ELSE
! original method, advect the bulk cloud fraction
            DO k = 1, wet_model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  cf_pc2_star(i,j,k)  = data_out_super(i,j,k,count+1)
                  cfl_pc2_star(i,j,k) = data_out_super(i,j,k,count+2)
                  cff_pc2_star(i,j,k) = data_out_super(i,j,k,count+3)
                  exner_star(i,j,k)   = data_out_super(i,j,k,count+4)
                END DO
              END DO
            END DO
          END IF ! l_fixbug_pc2_mixph
          count=count+4
        END IF  ! l_pc2

        If (l_mcr_cf2) Then
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf2_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_cf2

        If (l_mcr_rain) Then
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                rain_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_rain

        If (l_mcr_graup) Then
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                graup_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_graup

      End If  ! error_code=0

! ----------------------------------------------------------------------
! Section 4.0  Calculate turbulent variables at departure point.
!              Trajectory limit near ground already done in section 2.
!              Trajectory limit at top required if
!              model_levels > tke_levels - 1.
! ----------------------------------------------------------------------

      IF (Error_Code  ==  0 .AND.                                       &
                            bdy_tke > 0 .AND. l_adv_turb_field) THEN

! Perform trajectory limiter near top of model if
! model_levels > tke_levels - 1
        IF (model_levels  >   tke_levels-1 ) THEN

! is height field constant on this level,
! if so much cheaper code can be used.
          IF (tke_levels-1  >=  first_constant_r_rho_level) THEN

! assume once one level has no data above top then no lower level
! has either.
            min_k = MAX(1,tke_levels-1-interp_vertical_search_tol)
            k = tke_levels-1
            temp =1
            DO WHILE (k  >=  min_k .AND. temp  >   0 )
              temp = 0
              DO j = 1, rows
                DO i = 1, row_length
                  IF (depart_r(i,j,k)  >                                &
                      r_theta_levels(1,1,tke_levels-1) ) THEN
                    depart_r(i,j,k) =                                   &
                               r_theta_levels(1,1,tke_levels-1)
                    temp = temp + 1
                  END IF
                END DO
              END DO
              k = k - 1
            END DO

          ELSE

! calculate trajectory limit at departure point

! assume once one level has no data above top then no lower level
! has either.
            min_k = MAX(1,tke_levels-1-interp_vertical_search_tol)
            k = tke_levels - 1
            temp=1
            DO WHILE (k  >=  min_k .AND. temp  >   0 )
              temp = 0

! DEPENDS ON: bi_linear_h
              CALL bi_linear_h(                                         &
                            r_theta_levels(1-halo_i,1-halo_j,           &
                                           tke_levels - 1),             &
                            depart_lambda(1,1,k),                       &
                            depart_phi(1,1,k),                          &
                            row_length, rows, 1,                        &
                            row_length, rows, 1,                        &
                            row_length, rows,                           &
                            i_out(1,1,k), j_out(1,1,k),                 &
                            weight_lambda(1,1,k), weight_phi(1,1,k),    &
                            model_domain, me, n_procx, n_procy,n_proc,  &
                            halo_i, halo_j, l_datastart,                &
                            global_row_length, g_i_pe, at_extremity,    &
                            global_rows      , g_j_pe, l_2dcomm,        &
                            size_2dcomm, group_2dcomm,                  &
                            1, proc_row_group,                          &
                            L_sl_halo_reprod, L_regular,                &
                            work)

              DO j = 1, rows
                DO i = 1, row_length

                  IF (depart_r(i,j,k) >   work(i,j,1) ) THEN
! move trajectory down to highest level of data.
                    depart_r(i,j,k) = work(i,j,1)
                    temp = temp + 1
                  END IF

                END DO
              END DO
              CALL gc_isum (1,n_proc,error_code,temp)
              k = k - 1
            END DO

          END IF

        END IF
        

! ----------------------------------------------------------------------
! Section 4.1  Put values into super arrays
! ----------------------------------------------------------------------

!       e_trb(:,:,K) denotes valus on the (K-1)th theta level.
        DO k = 1, tke_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              super_array_etrb(i,j,k,1) = e_trb(i,j,k + 1)
            END DO
          END DO
        END DO
        count = 1

!       (:,:,K) denotes valus on the (K-1)th theta level.
        IF (bdy_tke == 3) THEN
          DO k = 1, tke_levels - 1
            DO j = 1, rows
              DO i = 1, row_length
                super_array_etrb(i,j,k,2) = tsq(i,j,k + 1)
                super_array_etrb(i,j,k,3) = qsq(i,j,k + 1)
                super_array_etrb(i,j,k,4) = cov(i,j,k + 1)
              END DO
            END DO
          END DO
          count = 4
        END IF

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
             super_array_etrb, row_length, rows,                        &
             count*(tke_levels - 1),                                    &
             halo_i, halo_j, fld_type_p, .FALSE. )

       IF ( (model_domain == mt_LAM .AND. L_lbc_old) .OR.               &
             model_domain == mt_cyclic_LAM ) THEN
! overwrite X_sl with X values on LAM boundaries
          IF (at_extremity(PSouth)) THEN
            DO k = 1, tke_levels - 1
              DO j = 1-halo_j, 1
                DO i = 1-halo_i, row_length + halo_i
                  super_array_etrb(i,j,k,1) = e_trb(i,j,k + 1)
                END DO
              END DO
            END DO
            IF (bdy_tke == 3) THEN
              DO k = 1, tke_levels - 1
                DO j = 1-halo_j, 1
                  DO i = 1-halo_i, row_length + halo_i
                    super_array_etrb(i,j,k,2) = tsq(i,j,k + 1)
                    super_array_etrb(i,j,k,3) = qsq(i,j,k + 1)
                    super_array_etrb(i,j,k,4) = cov(i,j,k + 1)
                  END DO
                END DO
              END DO
            END IF

          END IF  ! at_extremity PSouth

          IF (at_extremity(PNorth)) THEN
            DO k = 1, tke_levels - 1
              DO j = rows, rows+halo_j
                DO i = 1-halo_i, row_length + halo_i
                  super_array_etrb(i,j,k,1) = e_trb(i,j,k + 1)
                END DO
              END DO
            END DO
            IF (bdy_tke == 3) THEN
              DO k = 1, tke_levels - 1
                DO j = rows, rows+halo_j
                  DO i = 1-halo_i, row_length + halo_i
                    super_array_etrb(i,j,k,2) = tsq(i,j,k + 1)
                    super_array_etrb(i,j,k,3) = qsq(i,j,k + 1)
                    super_array_etrb(i,j,k,4) = cov(i,j,k + 1)
                  END DO
                END DO
              END DO
            END IF
          END IF  ! at_extremity PNorth

        END IF  ! model_domain (mt_lam or mt_cyclic lam)

        IF (model_domain == mt_LAM) THEN

         IF (L_lbc_old) THEN

          IF (at_extremity(PWest)) THEN
            DO k = 1, tke_levels - 1
              DO j = 1-halo_j, rows+halo_j
                DO i = 1-halo_i, 1
                  super_array_etrb(i,j,k,1) = e_trb(i,j,k + 1)
                END DO
              END DO
            END DO
            IF (bdy_tke == 3) THEN
              DO k = 1, tke_levels - 1
                DO j = 1-halo_j, rows+halo_j
                  DO i = 1-halo_i, 1
                    super_array_etrb(i,j,k,2) = tsq(i,j,k + 1)
                    super_array_etrb(i,j,k,3) = qsq(i,j,k + 1)
                    super_array_etrb(i,j,k,4) = cov(i,j,k + 1)
                  END DO
                END DO
              END DO
            END IF
          END IF  ! at_extremity PWest

          IF (at_extremity(PEast)) THEN

            DO k = 1, tke_levels - 1
              DO j = 1-halo_j, rows+halo_j
                DO i = row_length, row_length + halo_i
                  super_array_etrb(i,j,k,1) = e_trb(i,j,k + 1)
                END DO
              END DO
            END DO
            IF (bdy_tke == 3) THEN
              DO k = 1, tke_levels - 1
                DO j = 1-halo_j, rows+halo_j
                  DO i = row_length, row_length + halo_i
                    super_array_etrb(i,j,k,2) = tsq(i,j,k + 1)
                    super_array_etrb(i,j,k,3) = qsq(i,j,k + 1)
                    super_array_etrb(i,j,k,4) = cov(i,j,k + 1)
                  END DO
                END DO
              END DO
            END IF

          END IF  ! at_extremity PEast

         ELSE   !  new_lbcs

! DEPENDS ON: fill_external_halos
            CALL fill_external_halos(                                   &
                                super_array_etrb, row_length, rows,     &
                                etrb_array_size*(tke_levels-1),         &
                                     halo_i, halo_j)

         END IF  !  L_lbc_old


        END IF   ! model_domain   mt_lam

! ----------------------------------------------------------------------
! Section 4.2  Call interpolation routine.
! ----------------------------------------------------------------------

        L_vector = .FALSE.

! DEPENDS ON: interpolation_multi
        CALL Interpolation_multi(                                       &
                           super_array_etrb,                            &
                           eta_theta_levels(1),                         &
                           r_theta_levels(1-halo_i,1-halo_j,1),         &
                           mono_mass,                                   &
                           r_theta_levels(1-halo_i,1-halo_j,1),         &
                           count,                                       &
                           check_bottom_levels,                         &
                           interp_vertical_search_tol,                  &
                           first_constant_r_rho_level,                  &
                           row_length, rows, tke_levels - 1,            &
                           rows,                                        &
                           row_length, rows, tke_levels - 1,            &
                           delta_lambda, delta_phi,                     &
                           gdlambda_u, dphi_v,                          &
                           lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,&
                           recip_lambda_p_m, recip_lambda_p_0,          &
                           recip_lambda_p_p, recip_lambda_p_p2,         &
                           recip_phi_p_m, recip_phi_p_0,                &
                           recip_phi_p_p, recip_phi_p_p2,               &
                           i_out, j_out,                                &
                           weight_lambda, weight_phi,                   &
                           high_order_scheme_adv_turb,                  &
                           monotone_scheme_adv_turb,                    &
                           FV_cos_theta_latitude, L_regular,            &
                           model_domain, l_high_adv_turb,               &
                           l_mono_adv_turb, l_conserv_adv_turb,         &
                           depart_r, depart_lambda, depart_phi,         &

! Roar's bit
                           me, n_proc, n_procx, n_procy,                &
                           halo_i, halo_j, global_row_length,           &
                           global_rows, row_length, rows, n_rows,       &
                           l_datastart, at_extremity, g_i_pe,           &
                           g_j_pe,l_2dcomm,size_2dcomm,                 &
                           group_2dcomm,size_int, proc_all_group,       &
                           proc_row_group, proc_col_group, 1,           &
                           0, 0, L_sl_halo_reprod,                      &
                           off_x, off_y,                                &
                           data_out_super_etrb, Error_Code )
!
        DO k = 1, tke_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              e_trb(i,j,k + 1) = data_out_super_etrb(i,j,k,1)
            END DO
          END DO
        END DO
        IF (bdy_tke == 3) THEN
          DO k = 1, tke_levels - 1
            DO j = 1, rows
              DO i = 1, row_length
                tsq(i,j,k + 1) = data_out_super_etrb(i,j,k,2)
                qsq(i,j,k + 1) = data_out_super_etrb(i,j,k,3)
                cov(i,j,k + 1) = data_out_super_etrb(i,j,k,4)
              END DO
            END DO
          END DO
        END IF
      END IF ! error_code=0 .AND. ! bdy_tke > 0 .AND. l_adv_turb_field

! End of routine.
      IF (lhook) CALL dr_hook('SL_THERMO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SL_Thermo

