! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Full_Wind
!

      Subroutine SL_Full_Wind(                                          &
     &                         u, u_np1, v, v_np1, w, w_np1,            &
     &                         u_adv, v_adv, w_adv,                     &
     &                         thetav, thetav_np1, exner,               &
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
     &                         lambda_p_end, phi_p_end, dlambda_p_end,  &
     &                         dphi_p_end, dphi_v_end,                  &
     &                         recip_dlam, recip_dphi, max_look,        &
     &                         look_lam, look_phi, halo_lam, halo_phi,  &
     &                         alpha_3, alpha_4, LAM_max_cfl,           &
     &                         n_Y_arrays, n_Yw_arrays,                 &
     &                         n_Yd_arrays, n_Ydw_arrays,               &
     &                         u_lbcs,v_lbcs,w_lbcs,                    &
     &                         LENRIM,LBC_SIZE,LBC_START,RIMWIDTH,      &
     &                         model_domain,                            &
     &                         row_length, rows, n_rows, model_levels,  &
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
     &                         L_sl_halo_reprod, L_2d_sl_geometry,      &
     &                         L_free_slip, L_regular, L_interp_depart, &
     &                         L_lbc_old, L_new_tdisc, CycleNo,         &
     &                         R_u, R_v, R_w, Error_code)

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
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE earth_constants_mod, ONLY: g

      Use swapable_field_mod, Only: &
          swapable_field_pointer_type

      USE atmos_constants_mod, ONLY: cp
      
      USE conversions_mod, ONLY: pi

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
     &, L_interp_depart                                                 &
     &, L_regular

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
     &, base_lambda                                                     &
     &, base_phi                                                        &
     &, lambda_p_end                                                    &
     &, phi_p_end                                                       &
     &, dlambda_p_end                                                   &
     &, dphi_p_end                                                      &
     &, dphi_v_end                                                      &
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, timestep

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

      Real                                                              &
           ! primary model variables
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, u_np1(1-off_x:row_length+off_x,1-off_y:rows+off_y,              &
     &        model_levels)                                             &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, v_np1(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,           &
     &        model_levels)                                             &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels) &
     &, w_np1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        0:model_levels)                                           &
     &, u_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, v_adv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &         model_levels)                                            &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         0:model_levels)                                          &
     &, thetav (row_length+1, rows+1, model_levels)                     &
     &, thetav_np1 (row_length+1, rows+1, model_levels)                 &
     &, exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

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


! Local Variables.


      Integer                                                           &
     &  i, j, k                                                         &
                  ! Loop indices
     &, j0, j1                                                          &
                  ! derived from loop boundaries.
     &, gi

      Integer :: i_field     ! counter

      Real                                                              &
     &  weight1                                                         &
     &, weight2                                                         &
     &, weight3                                                         &
     &, term1                                                           &
     &, term2                                                           &
     &, recip_delta_lambda                                              &
     &, recip_delta_phi

      Real, Target ::                                                   &
     &  Yu (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      model_levels, n_Y_arrays)                                   &
     &, Yv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,         &
     &      model_levels, n_Y_arrays)                                   &
     &, Yw (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      0:model_levels, n_Yw_arrays) 

      Real                                                              &
     &  work(row_length+1, rows+1, model_levels)                        &
     &, thetav_rho_levels (row_length+1, rows+1, model_levels)          &
     &, ave_vpg(row_length+1, rows+1, first_constant_r_rho_level)

      Real                                                              &
     &  interp1 (row_length, rows)                                      &
     &, interp2 (row_length, rows)                                      &
     &, interp3 (row_length, rows)                                      &
     &, interp_wk (1:row_length+1, 0:rows+1)

      Real,Target ::                                                    &
     &  interp_work1 (1-off_x:row_length+off_x,                         &
     &                1-off_y:rows+off_y, model_levels)                 &
     &, interp_work2 (1-off_x:row_length+off_x,                         &
     &                1-off_y:rows+off_y, model_levels)                 &
     &, interp_work3 (1-off_x:row_length+off_x,                         &
     &                1-off_y:rows+off_y, model_levels) 
    
      Real ::                                                           &
     &  interp_work4 (row_length, rows,   model_levels)                 &
     &, interp_work5 (row_length, rows,   model_levels)                 &
     &, interp_work6 (row_length, rows,   model_levels)                 &
     &, interp_work7 (row_length, n_rows, model_levels)                 &
     &, interp_work8 (row_length, n_rows, model_levels)                 &
     &, interp_work9 (row_length, n_rows, model_levels)                 &
     &, work4 (row_length, rows)                                        &
     &, work5 (row_length, n_rows)                                      &
     &, top,bottom

      Integer                                                           &
     &  i_out (row_length, rows, 1)                                     &
     &, j_out (row_length, rows, 1)

      Real                                                              &
     &  weight_lambda (row_length, rows, 1)                             &
     &, weight_phi    (row_length, rows, 1)

      integer temp

      real min_phi,max_phi,min_lambda,max_lambda
      real min_phiv,max_phiv,min_lambdav,max_lambdav
      real lam_bnd, phi_bnd
      integer extra_i, extra_j

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)

      Logical                                                           &
     &  L_do_halos                                                      &
                          ! zero the halos?
     &, L_do_boundaries   ! zero the boundaries?

      REAL                                                              &
     &  rimweights(RIMWIDTH) ! dummy array to pass to
                             ! SET_LATERAL_BOUNDARIES

      Type(swapable_field_pointer_type) :: fields_to_swap(8)
      Real :: domain_size_x
      Real :: domain_size_y

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines:
      External                                                          &
     &  SL_Vector_u, SL_Vector_v, SL_Vector_w

! ----------------------------------------------------------------------
! Section 0.   Set r at u and v points. Set u at poles.
!              Calculate thetav and thetav on rho levels.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_FULL_WIND',zhook_in,zhook_handle)

      j0 = 1
      j1 = rows
      If (model_domain  ==  mt_Global) Then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      End If

      if ( L_regular) then
        lam_bnd = delta_lambda
        phi_bnd = delta_phi
        recip_delta_lambda = 1. / delta_lambda
        recip_delta_phi = 1. / delta_phi
      else  ! variable resolution
        lam_bnd = dlambda_p_end
        phi_bnd = dphi_p_end
      endif ! L_regular

      if(l_INTERP_DEPART) then
! Set Lateral boundary limits.
        If (model_domain  ==  mt_Global) Then
          If (L_sl_halo_reprod) Then
            min_lambda = 0.
            max_lambda = 2 * pi
          Else
          min_lambda = (2-halo_i) * lam_bnd
          max_lambda = 2 * pi + (halo_i-2) * lam_bnd
          End If
          min_phi = -pi*.5
          max_phi = pi*.5
          min_phiv=min_phi
          max_phiv=max_phi
          min_lambdav=min_lambda
          max_lambdav=max_lambda
        Else If (model_domain  ==  mt_LAM) Then
! allow for extra area defined by CFL parameter
          extra_i = (LAM_max_cfl(1) - 1)
          extra_j = (LAM_max_cfl(2) - 1)
! Check the LAM in the horizontal for each type.
! check horizontal in the LAM at u points.
          min_lambda = (0.5 - extra_i) * lam_bnd
          min_phi = base_phi - extra_j *  phi_bnd
          if ( L_regular) then
            max_lambda = (global_row_length-1.5+extra_i) * delta_lambda
            max_phi = base_phi + (global_rows-1+extra_j) * delta_phi
          else  ! variable resolution
            max_lambda = lambda_p_end - (0.5-extra_i) * dlambda_p_end - &
     &                    base_lambda
            max_phi = phi_p_end + extra_j * dphi_p_end
          endif ! L_regular
! check horizontal in the LAM at v points.
          min_lambdav = - extra_i *  lam_bnd
          min_phiv = base_phi+ (0.5-extra_j) * phi_bnd
          if ( L_regular) then
            max_lambdav = (global_row_length-1+extra_i) * delta_lambda
            max_phiv = base_phi + (global_rows-1.5+extra_j) * delta_phi
          else  ! variable resolution
            max_lambdav = lambda_p_end + extra_i * dlambda_p_end -      &
     &                     base_lambda
            max_phiv = phi_p_end -(0.5-extra_j)* dphi_v_end
          endif ! L_regular

        Else If (model_domain  ==  mt_cyclic_LAM) Then
! check horizontal in the periodic in x LAM.

          min_phi = base_phi
          min_phiv = base_phi+0.5 * phi_bnd
          min_lambda = (2-halo_i) * lam_bnd
          if ( L_regular) then
            max_lambda = (global_row_length+halo_i-2) * delta_lambda
            max_phi = base_phi + delta_phi * (global_rows-1)
            max_phiv = base_phi + delta_phi * (global_rows-1.5)
                    domain_size_x=global_row_length* delta_lambda 
          else  ! variable resolution
            max_lambda = lambda_p_end + (halo_i-1) * dlambda_p_end -    &
     &                   base_lambda
            domain_size_x=glambda_p(global_row_length+1)-glambda_p(1)
            max_phi = phi_p_end
            max_phiv = phi_p_end - 0.5* dphi_p_end
          endif ! L_regular
          min_lambdav=min_lambda
          max_lambdav=max_lambda
        Elseif (model_domain  ==  mt_bi_cyclic_LAM) then
          min_lambda = (2-halo_i) * lam_bnd
          min_phi = base_phi+(2-halo_j) * phi_bnd
          if ( L_regular) then
            max_lambda = (global_row_length + halo_i-2) * delta_lambda
            max_phi = base_phi + (global_rows+halo_j-2) * delta_phi
            domain_size_x=global_row_length* delta_lambda 
            domain_size_y=global_rows* delta_phi
          else  ! variable resolution
            max_lambda = lambda_p_end + (halo_i-1) * dlambda_p_end      &
     &                  -  base_lambda
            max_phi = phi_p_end + (halo_j-1) * dphi_p_end
            domain_size_x=glambda_p(global_row_length+1)-glambda_p(1)
            domain_size_y=phi_p_end-base_phi+phi_bnd 
          endif ! L_regular
          min_phiv=min_phi
          max_phiv=max_phi
          min_lambdav=min_lambda
          max_lambdav=max_lambda

        End If
      EndIf  ! L_interp_depart
      
! calculate vertical derivative at all points.
! store in Yw to save re-calculation.

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1, weight2,       &
!$OMP& interp1, interp2, interp3, interp_wk, term1, term2, weight3, gi)

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels - 1
        Do j = 1, j1+1
          Do i = 1, row_length+1
            Yw(i,j,k,1) = Cp * thetav(i,j,k) *                          &
     &                 (exner(i,j,k+1) - exner(i,j,k)) /                &
     &                 (r_rho_levels(i,j,k+1) -                         &
     &                  r_rho_levels(i,j,k))
          End Do
        End Do
      End Do
!$OMP END DO 

      k = 0

!$OMP DO SCHEDULE(STATIC)
        Do j = 1, j1+1
          Do i = 1, row_length+1
            Yw(i,j,k,1) = - g
          End Do
        End Do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, first_constant_r_rho_level - 1
! ave_vpg term at rho levels
        Do j = 1, j1+1
          Do i = 1, row_length+1
            ave_vpg(i,j,k) = ((r_rho_levels(i,j,k)                      &
     &                         - r_theta_levels(i,j,k-1))               &
     &                         * Yw(i,j,k,1) +                          &
     &                        (r_theta_levels(i,j,k)                    &
     &                         - r_rho_levels(i,j,k))                   &
     &                         * Yw(i,j,k-1,1) )                        &
     &                      / ( r_theta_levels(i,j,k) -                 &
     &                          r_theta_levels(i,j,k-1) )
          End Do
        End Do
      End Do
!$OMP END DO

      k = 1
!$OMP DO SCHEDULE(STATIC)
        Do j = 1, j1+1
          Do i = 1, row_length+1
! w at rho levels
            work(i,j,k) = ( r_rho_levels(i,j,k) -                       &
     &                      r_theta_levels(i,j,k-1) ) * w(i,j,k) /      &
     &                    ( r_theta_levels(i,j,k) -                     &
     &                      r_theta_levels(i,j,k-1) )

! thetav at rho levels
            thetav_rho_levels(i,j,k) = thetav(i,j,k)

          End Do
        End Do
!$OMP END DO
 
!$OMP DO SCHEDULE(STATIC)
      Do k = 2, model_levels - 1
        Do j = 1, j1+1
          Do i = 1, row_length+1
! w at rho levels
            weight1 = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1))/  &
     &                (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
            weight2 = 1.0 - weight1
            work(i,j,k) = weight1 * w(i,j,k) +  weight2 * w(i,j,k-1)
! thetav at rho levels
            thetav_rho_levels(i,j,k) = weight2 * thetav(i,j,k-1) +      &
     &                                 weight1 * thetav(i,j,k)

          End Do
        End Do
      End Do
!$OMP END DO

! ----------------------------------------------------
! cjs 081299 Removed isothermal layer. Since w is
!            always zero on top boundary, the following bit of
!            code could be removed and the above k
!            loop extended to k = model_levels.
! ----------------------------------------------------

      k = model_levels
!$OMP DO SCHEDULE(STATIC)
        Do j = 1, j1+1
          Do i = 1, row_length+1
            weight1 = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)) / &
     &                (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
            weight2 = 1.0 - weight1

! w at rho levels
            work(i,j,k) = weight2 * w(i,j,k-1)
            thetav_rho_levels(i,j,k) = (weight2 * thetav(i,j,k-1) +     &
     &                                  weight1 * thetav(i,j,k) )
          End Do
        End Do
!$OMP END DO

! ----------------------------------------------------------------------
! Section 1.  Calculate Yu at u points, Yv at v points and Yw at w
!             points, using the notation of working paper 162.
!             Set R_u, R_v and R_w to hold the value at the arrival
!             point. These equate, almost, to Xu, Xv and Xw. See 162 for
!             details of the implementation.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1  Calculate Yu at u points. R_u to X1 as defined in WP162.
! ----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels

! Calculate pressure gradient term and store in interp2.

        if ( L_regular) then
        If ( k >= first_constant_r_rho_level ) Then
          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) =  Cp * (thetav_rho_levels(i,j,k) +          &
     &                              thetav_rho_levels(i+1,j,k))         &
     &                     * (exner(i+1,j,k) - exner(i,j,k)) *          &
     &                          recip_delta_lambda                      &
     &                     * sec_theta_latitude(i,j)/                   &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i+1,j,k))
            End Do
          End Do

        Else ! k < first_constant_r_rho_level

          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) = ( Cp * (thetav_rho_levels(i,j,k) +         &
     &                               thetav_rho_levels(i+1,j,k))        &
     &                         * (exner(i+1,j,k) - exner(i,j,k))        &
     &                         - (r_rho_levels(i+1,j,k) -               &
     &                            r_rho_levels(i,j,k) ) *               &
     &                         (ave_vpg(i+1,j,k) + ave_vpg(i,j,k)) )    &
     &                     * recip_delta_lambda                         &
     &                     * sec_theta_latitude(i,j)/                   &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i+1,j,k))
            End Do
          End Do

        EndIf ! k >= first_constant_r_rho_level

! Store w at u points in interp1, v in interp3

        Do j = j0,j1
          Do i = 1, row_length
            interp1(i,j) = .5 * (work(i,j,k) + work(i+1,j,k))
            interp3(i,j) = .25 * (v(i,j,k) + v(i+1,j,k) +               &
     &                            v(i,j-1,k) + v(i+1,j-1,k))
          End Do
        End Do

        else  ! variable resolution
        If ( k >= first_constant_r_rho_level ) Then
          Do j = j0, j1
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              interp2(i,j) = Cp * ( thetav_rho_levels(i,j,k) *          &
     &                                             wt_lambda_p(i+1) +   &
     &                                   thetav_rho_levels(i+1,j,k) *   &
     &                                   (1.0 - wt_lambda_p(i+1)) ) *   &
     &                            ( exner(i+1,j,k) - exner(i,j,k) ) *   &
     &                   grecip_dlamp(gi) * sec_theta_latitude(i,j) /   &
     &                    ( r_rho_levels(i,j,k) * wt_lambda_p(i+1) +    &
     &             (1.0 - wt_lambda_p(i+1)) * r_rho_levels(i+1,j,k) )
            End Do
          End Do

        Else ! k < first_constant_r_rho_level
          Do j = j0, j1
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              interp2(i,j) = ( Cp * ( thetav_rho_levels(i,j,k) *        &
     &                                             wt_lambda_p(i+1) +   &
     &                                   thetav_rho_levels(i+1,j,k) *   &
     &                                   (1.0 - wt_lambda_p(i+1)) ) *   &
     &                                (exner(i+1,j,k) - exner(i,j,k)) - &
     &                  (r_rho_levels(i+1,j,k) - r_rho_levels(i,j,k)) * &
     &                            ( ave_vpg(i,j,k) * wt_lambda_p(i+1) + &
     &                                               ave_vpg(i+1,j,k) * &
     &                                   (1.0 - wt_lambda_p(i+1)) ) ) * &
     &                     grecip_dlamp(gi) * sec_theta_latitude(i,j) / &
     &                       ( r_rho_levels(i,j,k) * wt_lambda_p(i+1) + &
     &               (1.0 - wt_lambda_p(i+1)) * r_rho_levels(i+1,j,k) )
            End Do
          End Do

        EndIf ! k >= first_constant_r_rho_level

! Store w at u points in interp1,
! v in interp3 and rho in interp4

        Do j = j0,j1
          Do i = 1, row_length + 1
             interp_wk(i,j) = wt_phi_v(i,j) * v(i,j-1,k) +              &
     &                       (1.0 - wt_phi_v(i,j)) * v(i,j,k)
          End Do
        End Do

        Do j = j0,j1
          Do i = 1, row_length
            interp3(i,j) = wt_lambda_p(i+1) * interp_wk(i,j) +          &
     &                     (1.0 - wt_lambda_p(i+1)) * interp_wk(i+1,j)
            interp1(i,j) = wt_lambda_p(i+1) * work(i,j,k) +             &
     &                        (1.0 - wt_lambda_p(i+1)) * work(i+1,j,k)
          End Do
        End Do

      endif ! L_regular

! Form Yu, R_u from the calculated terms.
! Add R_u to Yu as this is the term to be interpolated to the
! departure point and R_u on input holds the physics increment
        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          weight2 = 1. - alpha_4
          weight3 = 1. - alpha_3
          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
          Do j = j0, j1
            Do i = 1, row_length
              term1 = interp3(i,j) * f3_at_u(i,j)                       &
     &              - interp2(i,j)
              term2 = interp1(i,j) * f2_at_u(i,j)

              Yu(i,j,k,1) = (term1 * weight3 - term2 * weight2 )        &
     &                    * timestep + u(i,j,k) + R_u(i,j,k)
              R_u(i,j,k) = (term1 * alpha_3 - term2 * alpha_4 )         &
     &                   * timestep - u(i,j,k)
            End Do
          End Do
          Else ! CycleNo == 1
            Do j = j0, j1
              Do i = 1, row_length
                term1 = interp3(i,j) * f3_at_u(i,j)                     &
     &                - interp2(i,j)
                term2 = interp1(i,j) * f2_at_u(i,j)
                Yu(i,j,k,1) = (term1 * weight3 - term2 * weight2 )      &
     &                      * timestep + u(i,j,k) + R_u(i,j,k)
              End Do
            End Do
          End If
        Else If (model_domain  ==  mt_LAM ) Then
          Do j = j0, j1
            Do i = 1, row_length
              term1 = interp3(i,j) * f3_at_u(i,j)                       &
     &              - interp2(i,j)
              term2 = interp1(i,j) * f2_at_u(i,j)

              Yu(i,j,k,1) = term1 * timestep
              Yu(i,j,k,2) = term2 * timestep
            End Do
          End Do

        End If  !  model_domain == mt_Global etc.

      End Do  !  k = 1, model_levels
!$OMP END DO

! northern and southern boundary values of Yu are set after section 1.2

! ----------------------------------------------------------------------
! Section 1.2  Calculate Yv at v points. R_v to X2 as defined in WP162.
! ----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels

! Now interpolate fields in horizontal and calculate horizontal
! derivatives.
! Store pressure gradient term at v points in interp2, u in interp3

        if ( L_regular) then

        If ( k < first_constant_r_rho_level ) Then
          Do j = 1, j1
            Do i = 1, row_length
              interp2(i,j) = ( Cp * (thetav_rho_levels(i,j+1,k) +       &
     &                               thetav_rho_levels(i,j,k))          &
     &                         * (exner(i,j+1,k) - exner(i,j,k))        &
     &                         - (r_rho_levels(i,j+1,k) -               &
     &                            r_rho_levels(i,j,k) ) *               &
     &                         (ave_vpg(i,j+1,k) + ave_vpg(i,j,k)) )    &
     &                     * recip_delta_phi /                          &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i,j+1,k))
            End Do
          End Do

        Else ! k >= first_constant_r_rho_level

          Do j = 1, j1
            Do i = 1, row_length
              interp2(i,j) =  Cp * (thetav_rho_levels(i,j,k) +          &
     &                              thetav_rho_levels(i,j+1,k))         &
     &                     * (exner(i,j+1,k) - exner(i,j,k)) *          &
     &                          recip_delta_phi /                       &
     &                       (r_rho_levels(i,j+1,k) +                   &
     &                        r_rho_levels(i,j,k))
            End Do
          End Do

        EndIf ! k < first_constant_r_rho_level

        Do j = 1, j1
          Do i = 1, row_length
              interp3(i,j) = .25 * (u(i,j,k) + u(i-1,j,k) +             &
     &                              u(i,j+1,k) + u(i-1,j+1,k))
          End Do
        End Do

        else  ! variable resolution
        If ( k < first_constant_r_rho_level ) Then
          Do j = 1, j1
            Do i = 1, row_length
              interp2(i,j) = ( Cp * (thetav_rho_levels(i,j,k) *         &
     &                                          wt_phi_p(i,j+1) +       &
     &                               thetav_rho_levels(i,j+1,k) *       &
     &                                 (1.0 -wt_phi_p(i,j+1)) ) *       &
     &                          (exner(i,j+1,k) - exner(i,j,k)) -       &
     &                  (r_rho_levels(i,j+1,k) - r_rho_levels(i,j,k)) * &
     &                              (ave_vpg(i,j,k) * wt_phi_p(i,j+1) + &
     &                                               ave_vpg(i,j+1,k) * &
     &                                     (1.0 -wt_phi_p(i,j+1)) ) ) * &
     &                                               recip_dphip(i,j) / &
     &                        ( r_rho_levels(i,j,k) * wt_phi_p(i,j+1) + &
     &                                          r_rho_levels(i,j+1,k) * &
     &                                         (1.0 -wt_phi_p(i,j+1)) )
            End Do
          End Do

        Else ! k >= first_constant_r_rho_level
          Do j = 1, j1
            Do i = 1, row_length
              interp2(i,j) = Cp * (thetav_rho_levels(i,j,k) *           &
     &                                              wt_phi_p(i,j+1) +   &
     &                                   thetav_rho_levels(i,j+1,k) *   &
     &                                    (1.0 - wt_phi_p(i,j+1)) ) *   &
     &                              (exner(i,j+1,k) - exner(i,j,k)) *   &
     &                                             recip_dphip(i,j) /   &
     &                     ( r_rho_levels(i,j,k) * wt_phi_p(i,j+1) +    &
     &                                       r_rho_levels(i,j+1,k) *    &
     &                                    (1.0 - wt_phi_p(i,j+1)) )
            End Do
          End Do

        EndIf ! k < first_constant_r_rho_level

        Do j = 1, j1+1
          Do i = 1, row_length
            interp_wk(i,j) = wt_lambda_u(i) * u(i-1,j,k) +              &
     &                             (1.0 - wt_lambda_u(i)) * u(i,j,k)
          End Do
        End Do

        Do j = 1, j1
          Do i = 1, row_length
            interp1(i,j) = wt_phi_p(i,j+1) * work(i,j,k) +              &
     &                    (1.0 - wt_phi_p(i,j+1)) * work(i,j+1,k)
            interp3(i,j) = wt_phi_p(i,j+1) * interp_wk(i,j) +           &
     &                     (1.0 - wt_phi_p(i,j+1)) * interp_wk(i,j+1)
          End Do
        End Do

      endif ! L_regular

! Form Yv, R_v from the calculated terms.
! Add R_v to Yv as this is the term to be interpolated to the
! departure point and R_v on input holds the physics increment
        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          weight2 = (1. - alpha_4)
          weight3 = (1. - alpha_3)

          if ( L_regular) then

          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
          Do j = 1, n_rows
            Do i = 1, row_length
              term1 = interp3(i,j) * f3_at_v(i,j)                       &
     &              + interp2(i,j)
              term2 = f1_at_v(i,j) * 0.5*(work(i,j,k) + work(i,j+1,k))

              Yv(i,j,k,1) = - (term1 * weight3 - term2 * weight2)       &
     &                      * timestep + v(i,j,k) + R_v(i,j,k)
              R_v(i,j,k) = - (term1 * alpha_3 - term2 * alpha_4)        &
     &                     * timestep - v(i,j,k)
            End Do
          End Do
          Else
            Do j = 1, n_rows
              Do i = 1, row_length
                term1 = interp3(i,j) * f3_at_v(i,j)                     &
     &                + interp2(i,j)
                term2 = f1_at_v(i,j) * 0.5*(work(i,j,k)+work(i,j+1,k))
                Yv(i,j,k,1) = - (term1 * weight3 - term2 * weight2)     &
     &                      * timestep + v(i,j,k) + R_v(i,j,k)
              End Do
            End Do
          End If

          else  ! variable resolution

          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

          Do j = 1, n_rows
            Do i = 1, row_length
                term1 = interp3(i,j) * f3_at_v(i,j) + interp2(i,j)
                term2 = f1_at_v(i,j) * interp1(i,j)
                Yv(i,j,k,1) = - (term1 * weight3 - term2 * weight2) *   &
     &                           timestep + v(i,j,k) + R_v(i,j,k)
                R_v(i,j,k) = - (term1 * alpha_3 - term2 * alpha_4) *    &
     &                          timestep - v(i,j,k)

            End Do
          End Do

          Else !  CycleNo > 1

            Do j = 1, n_rows
              Do i = 1, row_length
                term1 = interp3(i,j) * f3_at_v(i,j)                     &
     &                + interp2(i,j)
                term2 = f1_at_v(i,j) * interp1(i,j)
                Yv(i,j,k,1) = - (term1 * weight3 - term2 * weight2) *   &
     &                        timestep + v(i,j,k) + R_v(i,j,k)
              End Do
            End Do

          End If ! CycleNo == 1

          endif ! L_regular


        Else If (model_domain == mt_LAM ) Then

        if ( L_regular) then

          Do j = 1, n_rows
            Do i = 1, row_length
              term1 = interp3(i,j) * f3_at_v(i,j) + interp2(i,j)
              term2 = f1_at_v(i,j)*0.5 * (work(i,j,k) + work(i,j+1,k))
              Yv(i,j,k,1) = term1 * timestep
              Yv(i,j,k,2) = term2 * timestep
            End Do
          End Do

        else  ! variable resolution

          Do j = 1, n_rows
            Do i = 1, row_length
              term1 = interp3(i,j) * f3_at_v(i,j) + interp2(i,j)
              term2 = f1_at_v(i,j) * interp1(i,j)
              Yv(i,j,k,1) = term1 * timestep
              Yv(i,j,k,2) = term2 * timestep
            End Do
          End Do

        endif ! L_regular

        End If  !   model_domain == mt_Global

      End Do   !  k = 1, model_levels
!$OMP END DO 

!$OMP MASTER

! set northern and southern boundary values of Yu using Yv to define
! the cartesian polar values, as for the u wind field at the poles
      If (model_domain  ==  mt_Global) Then
! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       Yv,                                        &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       halo_i, halo_j, global_row_length,         &
     &                       proc_row_group, at_extremity)

        If (at_extremity(PSouth)) then
        if ( L_regular) then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              Yu (i,1,k,1) = - mag_vector_sp(k) * sin ( (gi-.5)*        &
     &                       delta_lambda - dir_vector_sp(k))
            End Do
          End Do
        else  ! variable resolution
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              Yu (i,1,k,1) = - mag_vector_sp(k) * sin ( glambda_u(gi) - &
     &                               base_lambda - dir_vector_sp(k))
            End Do
          End Do
        endif ! L_regular
        End If
        If (at_extremity(PNorth)) then
        if ( L_regular) then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              Yu (i,rows,k,1) = mag_vector_np(k) *                      &
     &                         sin ( (gi-.5)*delta_lambda -             &
     &                                dir_vector_np(k))
            End Do
          End Do
        else  ! variable resolution
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              Yu (i,rows,k,1) = mag_vector_np(k) * sin ( glambda_u(gi) -&
     &                               base_lambda - dir_vector_np(k))
            End Do
          End Do
        endif ! L_regular
        End If
      End If

!$OMP END MASTER
!$OMP BARRIER

! ----------------------------------------------------------------------
! Section 1.3  Calculate Yw at w points. R_w to X3 as defined in WP162.
! ----------------------------------------------------------------------

! first calculate f2 * u at w points and store in interp1.
! and f1 * v at w points and store in interp2.

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels - 1
! calculate f2 * u at w points except at poles
      if ( L_regular) then
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = .5 * (u(i,j,k+1) * f2_at_u(i,j) +                   &
     &                    u(i-1,j,k+1) * f2_at_u(i-1,j))
            term2 = .5 * (u(i,j,k) * f2_at_u(i,j) +                     &
     &                    u(i-1,j,k) * f2_at_u(i-1,j))
! interpolate in vertical
            interp1(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do
      else  ! variable resolution
! calculate f2 * u at w points except at poles
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = wt_lambda_u(i) * u(i-1,j,k+1) * f2_at_u(i-1,j) +    &
     &              (1.0 - wt_lambda_u(i)) * u(i,j,k+1) * f2_at_u(i,j)
            term2 = wt_lambda_u(i) * u(i-1,j,k) * f2_at_u(i-1,j) +      &
     &              (1.0 - wt_lambda_u(i)) * u(i,j,k) * f2_at_u(i,j)
! interpolate in vertical
            interp1(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do

      endif ! L_regular

!  set f2 * u, f1 * v fields at pole to zero.

        If (model_domain  ==  mt_Global) Then
          If (at_extremity(PSouth)) then
            Do i = 1, row_length
              interp1(i,1) = 0.
              interp2(i,1) = 0.
            End Do
          End If
          If (at_extremity(PNorth)) then
            Do i = 1, row_length
              interp1(i,rows) = 0.
              interp2(i,rows) = 0.
            End Do
          End If
        End If

      if ( L_regular) then

! calculate f1 * v at w points except at poles
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = .5 * (v(i,j,k+1) * f1_at_v(i,j) +                   &
     &                    v(i,j-1,k+1) * f1_at_v(i,j-1))
            term2 = .5 * (v(i,j,k) * f1_at_v(i,j) +                     &
     &                    v(i,j-1,k) * f1_at_v(i,j-1))
! interpolate in vertical
            interp2(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do

      else  ! variable resolution
! calculate f1 * v at w points except at poles
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = wt_phi_v(i,j) * v(i,j-1,k+1) * f1_at_v(i,j-1) +     &
     &              (1.0 - wt_phi_v(i,j)) * v(i,j,k+1) * f1_at_v(i,j)
            term2 = wt_phi_v(i,j) * v(i,j-1,k) * f1_at_v(i,j-1) +       &
     &              (1.0 - wt_phi_v(i,j)) * v(i,j,k) * f1_at_v(i,j)
! interpolate in vertical
            interp2(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do

      endif ! L_regular


! Now calculate Yw term.
        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
! calculate Yw for wet levels
          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
            Do j = 1, rows
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &              - g                                                 &
     &              - Yw(i,j,k,1)
                R_w (i,j,k) = term1 * alpha_4 * timestep - w(i,j,k)
                Yw(i,j,k,1) = term1*(1.-alpha_4) * timestep + w(i,j,k)
              End Do
            End Do

          Else ! If ( CycleNo == 1 .OR. .NOT. L_new_tdisc )
! store temporarily Yw (vpg) in interp_work1 to be used later
            Do j = 1, rows
              Do i = 1, row_length
                interp_work1(i,j,k) = Yw(i,j,k,1)
              End Do
            End Do
            Do j = 1, rows
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &                - g - Yw(i,j,k,1)
                Yw(i,j,k,1) = term1*(1.-alpha_4)*timestep + w(i,j,k)
              End Do
            End Do
          End If
        Else ! model_domain == mt_LAM
          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
! calculate Yw for wet levels
            Do j = 1, rows
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &              - g                                                 &
     &              - Yw(i,j,k,1)

                Yw (i,j,k,1) = term1 * timestep
              End Do
            End Do

          Else ! CycleNo > 1 .and. L_new_tdisc

            Do j = 1, rows
              Do i = 1, row_length
                interp_work1(i,j,k) = Yw(i,j,k,1)
              End Do
            End Do
            Do j = 1, rows
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &              - g                                                 &
     &              - Yw(i,j,k,1)
                Yw (i,j,k,1) = term1 * timestep
              End Do
            End Do

          End If ! CycleNo == 1 .OR. .NOT. L_new_tdisc

        End if ! model_domain

      End Do  !  k = 1, model_levels - 1
!$OMP END DO

! set values at top and bottom.
      if ( L_free_slip ) then

!$OMP DO SCHEDULE(STATIC)
        Do j = 1, rows
          Do i = 1, row_length
            Yw(i,j,0,1) = w(i,j,0)
            Yw(i,j,model_levels,1) = 0.
          End Do
        End Do
!$OMP END DO

      else   !  usual no-slip lower boundary condition

!$OMP DO SCHEDULE(STATIC)
      Do j = 1, rows
        Do i = 1, row_length
          Yw(i,j,0,1) = 0.
          Yw(i,j,model_levels,1) = 0.
        End Do
      End Do
!$OMP END DO

      End if   !  L_free_slip

      If (model_domain  ==  mt_LAM) Then
! Copy u, v and w into Yu(v)(w) arrays

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              Yu(i,j,k,3) = u(i,j,k) + R_u(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              Yv(i,j,k,3) = v(i,j,k) + R_v(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 0, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              Yw(i,j,k,2) = w(i,j,k)
            End Do
          End Do
        End Do  !  k = 0, model_levels
!$OMP END DO

      End If  !  model_domain  ==  mt_LAM

! ----------------------------------------------------------------------
! Section 1.4  Calculate R_u,R_v,R_w to X1, X2, X3 as defined in WP162
!              but for the iterative scheme. Yu, Yv, Yw have been
!              computed at sections 1.1, 1.2, 1.3.
! ----------------------------------------------------------------------
!
! Compute implicit part of the SISL scheme using
! available timelevel n+1 estimates from previous sweep.
! Currently only works for Global model
!
      If ( CycleNo > 1 .and. L_new_tdisc ) Then
!
! use theta^(1) for implicit vpg term in Ru, Rv
!
        k = 1
!$OMP DO SCHEDULE(STATIC)
        Do j = 1, j1+1
          Do i = 1, row_length+1
            weight1 = ( r_rho_levels(i,j,k) -                           &
     &                  r_theta_levels(i,j,k-1) ) /                     &
     &                ( r_theta_levels(i,j,k) -                         &
     &                  r_theta_levels(i,j,k-1) )
            weight2 = 1.0 - weight1
!
! Form ave_vpg(i,j,1) = weight1*Yw(i,j,1)+weight2*Yw(i,j,0)
! avoiding recomputing Yw (needed in w momentum below).
! Yw(,,0,1) is already defined = -g.
!
            ave_vpg(i,j,k) = weight1 * Cp * thetav_np1(i,j,k) *         &
     &                      ( exner(i,j,k+1) - exner(i,j,k) ) /         &
     &                      ( r_rho_levels(i,j,k+1) -                   &
     &                        r_rho_levels(i,j,k) ) +                   &
     &                       weight2*(-g)
          End Do
        End Do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 2, first_constant_r_rho_level - 1
          Do j = 1, j1+1
            Do i = 1, row_length+1
              weight1 = ( r_rho_levels(i,j,k) -                         &
     &                    r_theta_levels(i,j,k-1) ) /                   &
     &                  ( r_theta_levels(i,j,k) -                       &
     &                    r_theta_levels(i,j,k-1) )
              weight2 = 1.0 - weight1
!
! Form ave_vpg(i,j,k) = weight1*Yw(i,j,k)+weight2*Yw(i,j,k-1)
! avoiding recomputing Yw (needed in sub SL_vector_u below).
!
              ave_vpg(i,j,k) = weight1 * Cp * thetav_np1(i,j,k) *       &
     &                        ( exner(i,j,k+1) - exner(i,j,k) ) /       &
     &                        ( r_rho_levels(i,j,k+1) -                 &
     &                          r_rho_levels(i,j,k) ) +                 &
     &                         weight2 * Cp * thetav_np1(i,j,k-1) *     &
     &                        ( exner(i,j,k) - exner(i,j,k-1) ) /       &
     &                        ( r_rho_levels(i,j,k) -                   &
     &                          r_rho_levels(i,j,k-1) )
            End Do
          End Do
        End Do
!$OMP END DO
        k = 1

!$OMP DO SCHEDULE(STATIC)
          Do j = 1, j1+1
            Do i = 1, row_length+1
! w at rho levels
              work(i,j,k) = (r_rho_levels(i,j,k)                        &
     &                    - r_theta_levels(i,j,k-1) ) * w_np1(i,j,k) /  &
     &                     ( r_theta_levels(i,j,k)                      &
     &                    -  r_theta_levels(i,j,k-1) )
! thetav at rho levels
              thetav_rho_levels(i,j,k) = thetav_np1(i,j,k)
            End Do
          End Do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 2, model_levels - 1
          Do j = 1, j1+1
            Do i = 1, row_length+1
! w at rho levels
            weight1 = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1))/  &
     &                (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
            weight2 = 1.0 - weight1
            work(i,j,k) = weight1*w_np1(i,j,k)+weight2*w_np1(i,j,k-1)
! thetav at rho levels
            thetav_rho_levels(i,j,k) = weight2 * thetav_np1(i,j,k-1) +  &
     &                                 weight1 * thetav_np1(i,j,k)
            End Do
          End Do
        End Do  !  k = 2, model_levels - 1
!$OMP END DO

        k = model_levels

!$OMP DO SCHEDULE(STATIC)
          Do j = 1, j1+1
            Do i = 1, row_length+1
            weight1 = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)) / &
     &                (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
            weight2 = 1.0 - weight1
! w at rho levels
            work(i,j,k) = weight2 * w_np1(i,j,k-1)
            thetav_rho_levels(i,j,k) = (weight2 * thetav_np1(i,j,k-1) + &
     &                                  weight1 * thetav_np1(i,j,k) )
            End Do
          End Do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels

! Calculate pressure gradient term and store in interp2.
          If ( k >= first_constant_r_rho_level ) Then
            Do j = j0, j1
              Do i = 1, row_length
                interp2(i,j) =  Cp * (thetav_rho_levels(i,j,k) +        &
     &                                thetav_rho_levels(i+1,j,k))       &
     &                       * (exner(i+1,j,k) - exner(i,j,k)) *        &
     &                            recip_delta_lambda                    &
     &                       * sec_theta_latitude(i,j)/                 &
     &                         (r_rho_levels(i,j,k) +                   &
     &                          r_rho_levels(i+1,j,k))
              End Do
            End Do
          Else
            Do j = j0, j1
              Do i = 1, row_length
                interp2(i,j) = ( Cp * (thetav_rho_levels(i,j,k) +       &
     &                                 thetav_rho_levels(i+1,j,k))      &
     &                           * (exner(i+1,j,k) - exner(i,j,k))      &
     &                           - (r_rho_levels(i+1,j,k) -             &
     &                              r_rho_levels(i,j,k) ) *             &
     &                           (ave_vpg(i+1,j,k) + ave_vpg(i,j,k)) )  &
     &                       * recip_delta_lambda                       &
     &                       * sec_theta_latitude(i,j)/                 &
     &                         (r_rho_levels(i,j,k) +                   &
     &                          r_rho_levels(i+1,j,k))
              End Do
            End Do
          End If

! Form R_u from the calculated terms.
          Do j = j0, j1
            Do i = 1, row_length
              term1 = .25 * ( v(i,j,k) + v(i+1,j,k) +                   &
     &                        v(i,j-1,k) + v(i+1,j-1,k) )               &
     &                    * f3_at_u(i,j) - interp2(i,j)
              term2 = .5*(work(i,j,k) + work(i+1,j,k)) * f2_at_u(i,j)
              R_u(i,j,k) = ( term1 * alpha_3 - term2 * alpha_4 )        &
     &                   *   timestep - u(i,j,k)
            End Do
          End Do
!
          If (model_domain == mt_LAM .and. L_lbc_old) Then
! boundary area uses alphas = 1
! southern area
            If (at_extremity(PSouth)) Then
              Do j = 1, halo_j
                Do i = 1, row_length
                  term1 = .25 * ( v(i,j,k) + v(i+1,j,k) +               &
     &                            v(i,j-1,k) + v(i+1,j-1,k) )           &
     &                        * f3_at_u(i,j) - interp2(i,j)
                  term2 = .5*(work(i,j,k) + work(i+1,j,k))*f2_at_u(i,j)
                  R_u(i,j,k) = ( term1 - term2 )*timestep - u(i,j,k)
                End Do
              End Do
            End If
! western area
            If (at_extremity(PWest))  Then
              Do j = 1, rows
                Do i = 1, halo_i
                  term1 = .25 * ( v(i,j,k) + v(i+1,j,k) +               &
     &                            v(i,j-1,k) + v(i+1,j-1,k) )           &
     &                        * f3_at_u(i,j) - interp2(i,j)
                  term2 = .5*(work(i,j,k) + work(i+1,j,k))*f2_at_u(i,j)
                  R_u(i,j,k) = ( term1 - term2 )*timestep - u(i,j,k)
                End Do
              End Do
            End If
! eastern area
            If (at_extremity(PEast))  Then
              Do j = 1, rows
                Do i = row_length - halo_i, row_length
                  term1 = .25 * ( v(i,j,k) + v(i+1,j,k) +               &
     &                            v(i,j-1,k) + v(i+1,j-1,k) )           &
     &                        * f3_at_u(i,j) - interp2(i,j)
                  term2 = .5*(work(i,j,k) + work(i+1,j,k))*f2_at_u(i,j)
                  R_u(i,j,k) = ( term1 - term2 )*timestep - u(i,j,k)
                End Do
              End Do
            End If
! northern area
            If (at_extremity(PNorth))  Then
              Do j = rows-halo_j+1, rows
                Do i = 1, row_length
                  term1 = .25 * ( v(i,j,k) + v(i+1,j,k) +               &
     &                            v(i,j-1,k) + v(i+1,j-1,k) )           &
     &                        * f3_at_u(i,j) - interp2(i,j)
                  term2 = .5*(work(i,j,k) + work(i+1,j,k))*f2_at_u(i,j)
                  R_u(i,j,k) = ( term1 - term2 )*timestep - u(i,j,k)
                End Do
              End Do
            End If

          End If ! mt_LAM .and. L_lbc_old

        End Do  !  k = 1, model_levels
!$OMP END DO

! Compute R_v term
!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
! Now interpolate fields in horizontal and calculate horizontal
! derivatives.
! Store pressure gradient term at v points in interp2,
! u in interp3 and rho in interp4

          If ( k < first_constant_r_rho_level ) Then
            Do j = 1, j1
              Do i = 1, row_length
                interp2(i,j) = ( Cp * (thetav_rho_levels(i,j+1,k) +     &
     &                               thetav_rho_levels(i,j,k))          &
     &                         * (exner(i,j+1,k) - exner(i,j,k))        &
     &                         - (r_rho_levels(i,j+1,k) -               &
     &                            r_rho_levels(i,j,k) ) *               &
     &                         (ave_vpg(i,j+1,k) + ave_vpg(i,j,k)) )    &
     &                     * recip_delta_phi /                          &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i,j+1,k))
              End Do
            End Do
          Else
            Do j = 1, j1
              Do i = 1, row_length
                interp2(i,j) =  Cp * (thetav_rho_levels(i,j,k) +        &
     &                                thetav_rho_levels(i,j+1,k))       &
     &                       * (exner(i,j+1,k) - exner(i,j,k))          &
     &                       * recip_delta_phi /                        &
     &                         (r_rho_levels(i,j+1,k) +                 &
     &                          r_rho_levels(i,j,k))
              End Do
            End Do
          End If

          Do j = 1, n_rows
            Do i = 1, row_length
              term1 = .25 * (u(i,j,k) + u(i-1,j,k) +                    &
     &                      u(i,j+1,k) + u(i-1,j+1,k))* f3_at_v(i,j)    &
     &              + interp2(i,j)
              term2 = f1_at_v(i,j) * 0.5*(work(i,j,k) + work(i,j+1,k))
              R_v(i,j,k) = - (term1 * alpha_3 - term2 * alpha_4)        &
     &                     * timestep - v(i,j,k)
            End Do
          End Do

          If (model_domain  ==  mt_LAM .and. L_lbc_old ) Then
! boundary area uses alphas = 1
! southern area
            If (at_extremity(PSouth)) Then
              Do j = 1, halo_j
                Do i = 1, row_length
                  term1 = .25 * (u(i,j,k) + u(i-1,j,k) +                &
     &                          u(i,j+1,k) + u(i-1,j+1,k))*f3_at_v(i,j) &
     &                  + interp2(i,j)
                  term2 = f1_at_v(i,j) * 0.5*(work(i,j,k)+work(i,j+1,k))
                  R_v(i,j,k) = -( term1 - term2 )*timestep - v(i,j,k)
                End Do
              End Do
            End If
! western area
            If (at_extremity(PWest))  Then
              Do j = 1, n_rows
                Do i = 1, halo_i
                  term1 = .25 * (u(i,j,k) + u(i-1,j,k) +                &
     &                          u(i,j+1,k) + u(i-1,j+1,k))*f3_at_v(i,j) &
     &                  + interp2(i,j)
                  term2 = f1_at_v(i,j) * 0.5*(work(i,j,k)+work(i,j+1,k))
                  R_v(i,j,k) = -( term1 - term2 )*timestep - v(i,j,k)
                End Do
              End Do
            End If
! eastern area
            If (at_extremity(PEast))  Then
              Do j = 1, n_rows
                Do i = row_length-halo_i+1, row_length
                  term1 = .25 * (u(i,j,k) + u(i-1,j,k) +                &
     &                          u(i,j+1,k) + u(i-1,j+1,k))*f3_at_v(i,j) &
     &                  + interp2(i,j)
                  term2 = f1_at_v(i,j) * 0.5*(work(i,j,k)+work(i,j+1,k))
                  R_v(i,j,k) = -( term1 - term2 )*timestep - v(i,j,k)
                End Do
              End Do
            End If
! northern area
            If (at_extremity(PNorth))  Then
              Do j = n_rows-halo_j, n_rows
                Do i = 1, row_length
                  term1 = .25 * (u(i,j,k) + u(i-1,j,k) +                &
     &                          u(i,j+1,k) + u(i-1,j+1,k))*f3_at_v(i,j) &
     &                  + interp2(i,j)
                  term2 = f1_at_v(i,j) * 0.5*(work(i,j,k)+work(i,j+1,k))
                  R_v(i,j,k) = -( term1 - term2 )*timestep - v(i,j,k)
                End Do
              End Do
            End If

          End If ! mt_LAM .and. L_lbc_old

        End Do  !  k = 1, model_levels
!$OMP END DO

!-------------------------------------------------------------
! Compute R_w term
!-------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels - 1

! calculate f2 * u at w points except at poles
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = .5 * (u_np1(i,j,k+1) * f2_at_u(i,j) +               &
     &                    u_np1(i-1,j,k+1) * f2_at_u(i-1,j))
            term2 = .5 * (u_np1(i,j,k) * f2_at_u(i,j) +                 &
     &                    u_np1(i-1,j,k) * f2_at_u(i-1,j))
! interpolate in vertical
            interp1(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do

!  set f2 * u, f1 * v fields at pole to zero.

        If (model_domain == mt_Global) Then
          If (at_extremity(PSouth)) then
            Do i = 1, row_length
              interp1(i,1) = 0.
              interp2(i,1) = 0.
            End Do
          End If
          If (at_extremity(PNorth)) then
            Do i = 1, row_length
              interp1(i,rows) = 0.
              interp2(i,rows) = 0.
            End Do
          End If
        End If

! calculate f1 * v at w points except at poles
        Do j = j0, j1
          Do i = 1, row_length
! interpolate in horizontal
            term1 = .5 * (v_np1(i,j,k+1) * f1_at_v(i,j) +               &
     &                    v_np1(i,j-1,k+1) * f1_at_v(i,j-1))
            term2 = .5 * (v_np1(i,j,k) * f1_at_v(i,j) +                 &
     &                    v_np1(i,j-1,k) * f1_at_v(i,j-1))
! interpolate in vertical
            interp2(i,j) = ((r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)) * term1 +             &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k)) * term2 )/          &
     &                      (r_rho_levels(i,j,k+1) -                    &
     &                       r_rho_levels(i,j,k))
          End Do
        End Do

! Now calculate Rw term.
! Note that Rw is not recomputed using theta_np1.
! This is because of the coupling of w and theta eqn.
! Instead, theta^(1) term enters the Kterm which
! gives w'=w^{n+1}-w^n
!
        Do j = 1, rows
          Do i = 1, row_length
            term1 = interp1(i,j) - interp2(i,j)                         &
     &            - g -  interp_work1(i,j,k)
            R_w (i,j,k) = term1 * alpha_4 * timestep - w(i,j,k)
          End Do
        End Do

        If (model_domain == mt_LAM .and. L_lbc_old ) Then
! boundary area uses alphas = 1
! southern area
          If (at_extremity(PSouth)) Then
            Do j = 1, halo_j
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &                - g -  interp_work1(i,j,k)
                R_w (i,j,k) = term1 * timestep - w(i,j,k)
              End Do
            End Do
          End If
! western area
          If (at_extremity(PWest))  Then
            Do j = 1, rows
              Do i = 1, halo_i
                term1 = interp1(i,j) - interp2(i,j)                     &
     &                - g -  interp_work1(i,j,k)
                R_w (i,j,k) = term1 * timestep - w(i,j,k)
              End Do
            End Do
          End If
! eastern area
          If (at_extremity(PEast))  Then
            Do j = 1, rows
              Do i = row_length-halo_i+1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &                - g -  interp_work1(i,j,k)
                R_w (i,j,k) = term1 * timestep - w(i,j,k)
              End Do
            End Do
          End If
! northern area
          If (at_extremity(PNorth))  Then
            Do j =  rows-halo_j+1, rows
              Do i = 1, row_length
                term1 = interp1(i,j) - interp2(i,j)                     &
     &                - g -  interp_work1(i,j,k)
                R_w (i,j,k) = term1 * timestep - w(i,j,k)
              End Do
            End Do
          End If

        End If ! model_domain == mt_LAM .and. L_lbc_old

      End Do  !  k = 1, model_levels - 1
!$OMP END DO


      End If ! CycleNo

!$OMP END PARALLEL 
! ----------------------------------------------------------------------
! Section 2.   Call SL_Vector_u to calculate R_u term.
! ----------------------------------------------------------------------

      i_field = 0
      Do k = 1, n_Y_arrays
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => Yu(:,:,:,k)
        fields_to_swap(i_field) % field_type  =  fld_type_u
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => Yv(:,:,:,k)
        fields_to_swap(i_field) % field_type  =  fld_type_v
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  n_rows
        fields_to_swap(i_field) % vector      =  .TRUE.
      End Do

      Do k = 1, n_Yw_arrays
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => Yw(:,:,:,k)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels+1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
      End Do

! DEPENDS ON: swap_bounds_mv
      Call SWAP_BOUNDS_MV( fields_to_swap, i_field, row_length,         &
                           halo_i, halo_j )

      If (model_domain == mt_LAM) Then

        L_do_boundaries=.FALSE.
        L_do_halos=.TRUE.

! Add halo values to the Yu on the edge processors
! Set to zero for forcing terms but to lateral boundary values for u
        Do k = 1, n_Y_arrays - 1

! DEPENDS ON: zero_lateral_boundaries
            CALL ZERO_LATERAL_BOUNDARIES(                               &
     &        ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_u,    &
     &        Yu(1-halo_i,1-halo_j,1,k), 1, AT_EXTREMITY,               &
     &        L_do_boundaries, L_do_halos)

! DEPENDS ON: zero_lateral_boundaries
            CALL ZERO_LATERAL_BOUNDARIES(                               &
     &        ROW_LENGTH,N_ROWS,halo_i,halo_j,MODEL_LEVELS,fld_type_v,  &
     &        Yv(1-halo_i,1-halo_j,1,k), 1, AT_EXTREMITY,               &
     &        L_do_boundaries, L_do_halos)

        End Do  !  k = 1, n_Y_arrays - 1

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,ROWS,halo_i,halo_j,MODEL_LEVELS,fld_type_u,    &
     &        Yu(1-halo_i,1-halo_j,1,3),                                &
     &        LENRIM(fld_type_u,halo_type_extended),                    &
     &        LBC_SIZE(1,fld_type_u,halo_type_extended),                &
     &        LBC_START(1,fld_type_u,halo_type_extended),               &
     &        halo_i,halo_j,u_lbcs,                                     &
     &        RIMWIDTH,1,rimweights,AT_EXTREMITY,                       &
     &        L_do_boundaries, L_do_halos)

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,N_ROWS,halo_i,halo_j,MODEL_LEVELS,fld_type_v,  &
     &        Yv(1-halo_i,1-halo_j,1,3),                                &
     &        LENRIM(fld_type_v,halo_type_extended),                    &
     &        LBC_SIZE(1,fld_type_v,halo_type_extended),                &
     &        LBC_START(1,fld_type_v,halo_type_extended),               &
     &        halo_i,halo_j,v_lbcs,                                     &
     &        RIMWIDTH,1,rimweights,AT_EXTREMITY,                       &
     &        L_do_boundaries, L_do_halos)

! DEPENDS ON: zero_lateral_boundaries
            CALL ZERO_LATERAL_BOUNDARIES(                               &
     &        ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,  &
     &        Yw(1-halo_i,1-halo_j,0,1),1,AT_EXTREMITY,                 &
     &        L_do_boundaries,L_do_halos)

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,ROWS,halo_i,halo_j,MODEL_LEVELS,fld_type_p,    &
     &        Yw(1-halo_i,1-halo_j,0,2),                                &
     &        LENRIM(fld_type_p,halo_type_extended),                    &
     &        LBC_SIZE(1,fld_type_p,halo_type_extended),                &
     &        LBC_START(1,fld_type_p,halo_type_extended),               &
     &        halo_i,halo_j,w_lbcs,                                     &
     &        RIMWIDTH, 1, rimweights, AT_EXTREMITY,                    &
     &        L_do_boundaries, L_do_halos)

      End If  !  model_domain == mt_LAM


      IF(L_INTERP_DEPART) then

!!!!! interpolate departure points for w points to u,v points
!! interpolate from w (theta levels) on to rho levels

        Do k = 2, model_levels
          Do j = 1, rows
            Do i = 1, row_length

          weight1 = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k) )/     &
     &              (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1) )
          weight3 = 1.0 - weight1

              interp_work1 (i,j,k) =                                    &
     &           (weight1  * depart_lambda_w(i,j,k-1) +                 &
     &            weight3  * depart_lambda_w(i,j,k)   )
              interp_work2 (i,j,k) =                                    &
     &           (weight1  * depart_phi_w(i,j,k-1) +                    &
     &            weight3  * depart_phi_w(i,j,k)   )
              interp_work3 (i,j,k) =                                    &
     &           (weight1  * depart_r_w(i,j,k-1) +                      &
     &            weight3  * depart_r_w(i,j,k)   )
              if(interp_work3 (i,j,k)  <  r_rho_levels(i,j,1))then
                interp_work3 (i,j,k)= r_rho_levels(i,j,1)
              endif
              if(interp_work3 (i,j,k)  >                                &
     &             r_rho_levels(i,j,model_levels))then
                interp_work3 (i,j,k)= r_rho_levels(i,j,model_levels)
              endif

            End Do
          End Do
        End do

! DEPENDS ON: calc_index
              Call Calc_Index(                                          &
     &                      row_length, rows, 1,                        &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular,                                  &
     &                      depart_lambda_w(1,1,1), depart_phi_w(1,1,1),&
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

!!!    need to set up arrays on level 1
! DEPENDS ON: bi_linear_h
              Call bi_linear_h(                                         &
     &                      r_theta_levels(1-halo_i,1-halo_j,0),        &
     &                      depart_lambda_w(1,1,1),                     &
     &                      depart_phi_w(1,1,1),                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain, me, n_procx, n_procy,n_proc,  &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      1, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      work4)
!! work4 holds surface at departure point coordinate

        k = 1
          Do j = 1, rows
            Do i = 1, row_length

              interp_work1 (i,j,k) = depart_lambda_w(i,j,k)
              interp_work2 (i,j,k) = depart_phi_w(i,j,k)
              interp_work3 (i,j,k) = depart_r_w(i,j,k) -                &
     &                   r_theta_levels(i,j,1) + r_rho_levels(i,j,1)

              if(interp_work3 (i,j,k)  <  work4(i,j))then
                interp_work3 (i,j,k)= work4(i,j)
              endif

            End Do
          End Do

        i_field = 0

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => interp_work1(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => interp_work2(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => interp_work3(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

! DEPENDS ON: swap_bounds_mv
        Call SWAP_BOUNDS_MV( fields_to_swap, i_field, row_length,       &
                           off_x, off_y)

!!  section below needed for 0, 2pi neighbouring points

        if(model_domain == mt_global)then
        if(at_extremity(PEast))then
        Do k = 1, model_levels
          Do j = 1, rows
            i =  row_length  + off_x
            if (abs(interp_work1(i,j,k)-interp_work1(i-1,j,k))          &
     &           >  pi)then
              interp_work1(i,j,k) =interp_work1(i,j,k) +2.0*pi
            endif
          End Do
        End Do
        endif
        if(at_extremity(Pwest))then
        Do k = 1, model_levels
          Do j = 1, rows
            i =  1- off_x
            if (abs(interp_work1(i,j,k)-interp_work1(i+1,j,k))          &
     &           >  Pi)then
              interp_work1(i,j,k) =interp_work1(i,j,k) -2.0*pi
            endif
          End Do
        End Do
        endif
        endif !model_domain

! use andy whites formula to find departure points at u points

        Do k = 1, model_levels
          Do j = 1, rows       !ajm
            Do i = 1, row_length
              interp_work6 (i,j,k)= .5*(interp_work3 (i,j,k) +          &
     &                                  interp_work3 (i+1,j,k) )

            top= cos (interp_work2 (i,j,k))*sin(interp_work1 (i,j,k))   &
     &      + cos (interp_work2 (i+1,j,k))*sin(interp_work1 (i+1,j,k))
            bottom= cos (interp_work2 (i,j,k))*cos(interp_work1 (i,j,k))&
     &    + cos (interp_work2 (i+1,j,k))*cos(interp_work1 (i+1,j,k))

              interp_work4(i,j,k)=atan2(top,bottom)
              if(interp_work4(i,j,k) <   -delta_lambda/2.)              &
     &        interp_work4(i,j,k)=interp_work4(i,j,k)+2.*pi

              interp_work5 (i,j,k)= asin(                               &
     &  (sin(interp_work2 (i,j,k)) + sin(interp_work2 (i+1,j,k))) /     &
     &  (sqrt(2.) * sqrt( 1.+                                           &
     &    sin(interp_work2 (i,j,k))*sin(interp_work2 (i+1,j,k))         &
     & + cos(interp_work2 (i,j,k))*cos(interp_work2 (i+1,j,k)) *        &
     & cos( interp_work1 (i+1,j,k) -interp_work1 (i,j,k) ) ) ) )
            End Do
          End Do
        End do

!    the code below is very inefficient
        if(model_domain == mt_global)then
        if(at_extremity(PEast))then
        Do k = 1, model_levels
          Do j = 1, rows
            do i =  row_length-halo_i,row_length
            if (interp_work4(i,j,k) <  0.0) then
              interp_work4(i,j,k) =interp_work4(i,j,k) +2.0*pi
            endif
            End Do
          End Do
        End Do
        endif
        if(at_extremity(Pwest))then
        Do k = 1, model_levels
          Do j = 1, rows
            do i =  1,halo_i
            if (interp_work4(i,j,k) >  pi)then
              interp_work4(i,j,k) =interp_work4(i,j,k) -2.0*pi
            endif
          End Do
          End Do
        End Do
        endif

!  end of inefficient code

         if(at_extremity(PSouth))then
         j=1
        Do k = 1, model_levels
         Do i = 1, row_length
          interp_work4 (i,j,k)= interp_work1 (i,j,k) +delta_lambda/2.
          interp_work5 (i,j,k)= interp_work2 (i,j,k)
          interp_work6 (i,j,k)= interp_work3 (i,j,k)
         End do
         End do
          endif
         if(at_extremity(PNorth))then
         j=rows
        Do k = 1, model_levels
         Do i = 1, row_length
          interp_work4 (i,j,k)= interp_work1 (i,j,k) +delta_lambda/2.
          interp_work5 (i,j,k)= interp_work2 (i,j,k)
          interp_work6 (i,j,k)= interp_work3 (i,j,k)
         End do
        End do
          endif
        endif   ! model_domain

!!!!! AJM mes fix 27/05/04
         if(model_domain == mt_LAM)then
         if (at_extremity(PEast))then
          do k=1,model_levels
          do j=1,rows
          interp_work4 (row_length,j,k)=                                &
     &            interp_work1 (row_length,j,k) +delta_lambda/2.
          End do
          End do
         endif
         endif
!! AJM mes fix ends 27/05/04

! DEPENDS ON: check_sl_domain
         Call check_sl_domain(                                          &
     &                         model_domain, interp_work5,interp_work4, &
     &                         row_length, rows, model_levels,          &
     &                         domain_size_x, domain_size_y,            &
     &                         max_lambda, min_lambda,                  &
     &                         max_phi, min_phi )

! move points below bottom rho level up to bottom rho level
           k=1
           temp=1
            Do while ( temp  >   0 )
              temp = 0

! DEPENDS ON: calc_index
              Call Calc_Index(                                          &
     &                      row_length, rows, 1,                        &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular,                                  &
     &                      interp_work4(1,1,k), interp_work5(1,1,k),   &
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

! DEPENDS ON: bi_linear_h
              Call bi_linear_h(                                         &
     &                      r_rho_levels(1-halo_i,1-halo_j,1),          &
     &                      interp_work4(1,1,k),                        &
     &                      interp_work5(1,1,k),                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain, me, n_procx, n_procy,n_proc,  &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      1, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      work4)

              Do j = 1, rows
                Do i = 1, row_length
                  If (interp_work6(i,j,k) <   work4(i,j) ) Then
! move trajectory up to lowest level of data.
                    interp_work6(i,j,k) = work4(i,j)
                    temp = temp + 1
                  End If

                End Do     ! i loop
              End Do       ! j loop
              k = k + 1
              call gc_isum (1,n_proc,error_code,temp)
            End Do         ! do while
      ENDIF  ! L_INTERP_DEPART

! DEPENDS ON: sl_vector_u
      Call SL_Vector_u(                                                 &
     &                  u, u_adv, v_adv, w_adv, Yu, Yv, Yw,             &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  r_at_u, r_at_v,                                 &
     &                  interp_work4,interp_work5,interp_work6,         &
     &                  L_regular, L_interp_depart,                     &
     &                  cos_theta_latitude, sec_theta_latitude,         &
     &                  sin_theta_latitude, cos_v_latitude,             &
     &                  sec_v_latitude, sin_v_latitude,                 &
     &                  tan_theta_latitude, tan_v_latitude,             &
     &                  cos_theta_longitude,                            &
     &                  sin_theta_longitude,                            &
     &                  delta_lambda, delta_phi,                        &
     &                  glambda_p, phi_p, glambda_u, phi_v,             &
     &                  gdlambda_p, dphi_p, gdlambda_u, dphi_v,         &
     &                  grecip_dlamp, recip_dphip,                      &
     &                  grecip_dlamu, recip_dphiv, wt_lambda_p,         &
     &                  wt_phi_p, wt_lambda_u, wt_phi_v,                &
     &                  lambda_p_rm, lambda_p_rp,                       &
     &                  lambda_u_rm, lambda_u_rp,                       &
     &                  phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,         &
     &                  recip_lambda_p_m, recip_lambda_p_0,             &
     &                  recip_lambda_p_p, recip_lambda_p_p2,            &
     &                  recip_lambda_u_m, recip_lambda_u_0,             &
     &                  recip_lambda_u_p, recip_lambda_u_p2,            &
     &                  recip_phi_p_m, recip_phi_p_0,                   &
     &                  recip_phi_p_p, recip_phi_p_p2,                  &
     &                  recip_phi_v_m, recip_phi_v_0,                   &
     &                  recip_phi_v_p, recip_phi_v_p2,                  &
     &                  Base_lambda, base_phi,                          &
     &                  lambda_p_end, phi_p_end,                        &
     &                  dlambda_p_end, dphi_p_end, dphi_v_end,          &
     &                  recip_dlam, recip_dphi, max_look,               &
     &                  look_lam, look_phi, halo_lam, halo_phi,         &
     &                  n_Y_arrays, n_Yd_arrays, n_Ydw_arrays,          &
     &                  timestep, alpha_3, alpha_4, LAM_max_cfl,        &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  Depart_scheme, Depart_order,                    &
     &                  high_order_scheme, monotone_scheme,             &
     &                  model_domain, L_2d_sl_geometry, L_high, L_mono, &
     &                  L_conserv, L_Robert_high, L_Robert_mono,        &
     &                  Robert_high_order_scheme,                       &
     &                  Robert_monotone_scheme,                         &
     &                  check_bottom_levels,                            &
     &                  interp_vertical_search_tol,                     &
     &                  first_constant_r_rho_level,                     &
     &                  L_trivial_trigs,                                &
     &                  me, n_proc, n_procx, n_procy,                   &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  global_row_length, global_rows,                 &
     &                  l_datastart, at_extremity, g_i_pe, g_j_pe,      &
     &                  l_2dcomm, size_2dcomm,                          &
     &                  group_2dcomm, proc_row_group, proc_col_group,   &
     &                  proc_all_group,                                 &
     &                  L_sl_halo_reprod, L_lbc_old,                    &
     &                  L_new_tdisc, CycleNo, R_u,                      &
     &                  Error_code)

! ----------------------------------------------------------------------
! Section 3.   Call SL_Vector_v to calculate R_v term.
! ----------------------------------------------------------------------

      If (Error_code  ==  0) Then
      IF(L_INTERP_DEPART) then
!! may change 7,8,9 to 4,5,6

! use andy whites formula to find departure points at v points

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
            top= cos (interp_work2 (i,j,k))*sin(interp_work1 (i,j,k))   &
     &      + cos (interp_work2 (i,j+1,k))*sin(interp_work1 (i,j+1,k))
            bottom= cos (interp_work2 (i,j,k))*cos(interp_work1 (i,j,k))&
     &    + cos (interp_work2 (i,j+1,k))*cos(interp_work1 (i,j+1,k))

              interp_work7(i,j,k)=atan2(top,bottom)
              if(interp_work7(i,j,k) <   -delta_lambda/2.)              &
     &        interp_work7(i,j,k)=interp_work7(i,j,k)+2.*pi

              interp_work8 (i,j,k)= asin(                               &
     &  (sin(interp_work2 (i,j,k)) + sin(interp_work2 (i,j+1,k))) /     &
     &  (sqrt(2.) * sqrt( 1.+                                           &
     &    sin(interp_work2 (i,j,k))*sin(interp_work2 (i,j+1,k))         &
     & + cos(interp_work2 (i,j,k))*cos(interp_work2 (i,j+1,k)) *        &
     & cos( interp_work1 (i,j+1,k) -interp_work1 (i,j,k) ) ) ) )

              interp_work9 (i,j,k)= .5*(interp_work3 (i,j,k) +          &
     &                                  interp_work3 (i,j+1,k))
            End Do
          End Do
        End do

!    the code below is very inefficient
        if(model_domain == mt_global)then
        if(at_extremity(PEast))then
        Do k = 1, model_levels
          Do j = 1, n_rows
            do i =  row_length-halo_i,row_length
            if (interp_work7(i,j,k) <  0.0) then
              interp_work7(i,j,k) =interp_work7(i,j,k) +2.0*pi
            endif
            End Do
          End Do
        End Do
        endif
        if(at_extremity(Pwest))then
        Do k = 1, model_levels
          Do j = 1, n_rows
            do i =  1,halo_i
            if (interp_work7(i,j,k) >  pi)then
              interp_work7(i,j,k) =interp_work7(i,j,k) -2.0*pi
            endif
          End Do
          End Do
        End Do
        endif
        endif !model_domain
!  end of inefficient code

!!!!! AJM mes fix 27/05/04
         if(model_domain == mt_LAM)then
         if (at_extremity(PNorth))then
          do k=1,model_levels
          do i=1,row_length
          interp_work8 (i,n_rows,k)=                                    &
     &            interp_work2 (i,n_rows,k) +delta_phi/2.
          End do
          End do
         endif
         endif
!! AJM mes fix ends 27/05/04

! DEPENDS ON: check_sl_domain
          Call check_sl_domain(                                         &
     &                         model_domain, interp_work8,interp_work7, &
     &                         row_length, n_rows, model_levels,        &
     &                         domain_size_x, domain_size_y,            &
     &                         max_lambdav, min_lambdav,                &
     &                         max_phiv, min_phiv )

! don't need to move pole values for v, are directly calculated

! move points below bottom rho level up to bottom rho level
           k=1
           temp=1
            Do while ( temp  >   0 )
              temp = 0

! DEPENDS ON: calc_index
              Call Calc_Index(                                          &
     &                      row_length, n_rows, 1,                      &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular,                                  &
     &                      interp_work7(1,1,k), interp_work8(1,1,k),   &
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

! DEPENDS ON: bi_linear_h
              Call bi_linear_h(                                         &
     &                      r_rho_levels(1-halo_i,1-halo_j,1),          &
     &                      interp_work7(1,1,k),                        &
     &                      interp_work8(1,1,k),                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, n_rows, 1,                      &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain, me, n_procx, n_procy,n_proc,  &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      2, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      work5)

              Do j = 1, n_rows
                Do i = 1, row_length
                  If (interp_work9(i,j,k) <   work5(i,j) ) Then
! move trajectory up to lowest level of data.
                    interp_work9(i,j,k) = work5(i,j)
                    temp = temp + 1
                  End If

                End Do     ! i loop
              End Do       ! j loop
              k = k + 1
              call gc_isum (1,n_proc,error_code,temp)
            End Do         ! do while
      ENDIF  ! L_INTERP_DEPART


! DEPENDS ON: sl_vector_v
        Call SL_Vector_v(                                               &
     &                    v, u_adv, v_adv, w_adv, Yu, Yv, Yw,           &
     &                    eta_theta_levels, eta_rho_levels,             &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    r_at_u, r_at_v,                               &
     &                    interp_work7, interp_work8, interp_work9,     &
     &                    L_regular, L_interp_depart,                   &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    sin_theta_latitude, cos_v_latitude,           &
     &                    sec_v_latitude, sin_v_latitude,               &
     &                    tan_theta_latitude, tan_v_latitude,           &
     &                    cos_theta_longitude,                          &
     &                    sin_theta_longitude,                          &
     &                    delta_lambda, delta_phi,                      &
     &                    glambda_p, phi_p, glambda_u, phi_v,           &
     &                    gdlambda_p, dphi_p, gdlambda_u, dphi_v,       &
     &                    grecip_dlamp, recip_dphip,                    &
     &                    grecip_dlamu, recip_dphiv, wt_lambda_p,       &
     &                    wt_phi_p, wt_lambda_u, wt_phi_v,              &
     &                    lambda_p_rm, lambda_p_rp,                     &
     &                    lambda_u_rm, lambda_u_rp,                     &
     &                    phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,       &
     &                    recip_lambda_p_m, recip_lambda_p_0,           &
     &                    recip_lambda_p_p, recip_lambda_p_p2,          &
     &                    recip_lambda_u_m, recip_lambda_u_0,           &
     &                    recip_lambda_u_p, recip_lambda_u_p2,          &
     &                    recip_phi_p_m, recip_phi_p_0,                 &
     &                    recip_phi_p_p, recip_phi_p_p2,                &
     &                    recip_phi_v_m, recip_phi_v_0,                 &
     &                    recip_phi_v_p, recip_phi_v_p2,                &
     &                    Base_lambda, base_phi,                        &
     &                    lambda_p_end, phi_p_end, dlambda_p_end,       &
     &                    dphi_p_end, dphi_v_end,                       &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    n_Y_arrays, n_Yd_arrays, n_Ydw_arrays,        &
     &                    timestep, alpha_3, alpha_4, LAM_max_cfl,      &
     &                    row_length, rows, n_rows, model_levels,       &
     &                    Depart_scheme, Depart_order,                  &
     &                    high_order_scheme, monotone_scheme,           &
     &                    model_domain, L_2d_sl_geometry,               &
     &                    L_high, L_mono,                               &
     &                    L_conserv, L_Robert_high, L_Robert_mono,      &
     &                    Robert_high_order_scheme,                     &
     &                    Robert_monotone_scheme,                       &
     &                    check_bottom_levels,                          &
     &                    interp_vertical_search_tol,                   &
     &                    first_constant_r_rho_level,                   &
     &                    L_trivial_trigs,                              &
     &                    me, n_proc, n_procx, n_procy,                 &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    global_row_length, global_rows,               &
     &                    l_datastart, at_extremity, g_i_pe, g_j_pe,    &
     &                    l_2dcomm, size_2dcomm,                        &
     &                    group_2dcomm, proc_row_group, proc_col_group, &
     &                    proc_all_group,                               &
     &                    L_sl_halo_reprod, L_lbc_old,                  &
     &                    L_new_tdisc, CycleNo, R_v,                    &
     &                    Error_code)

      End If

! ----------------------------------------------------------------------
! Section 4.   Call SL_Vector_w to calculate R_w term.
! ----------------------------------------------------------------------

      If (Error_code  ==  0) Then

! Call SL_Vector_w to calculate advection of vertical velocity.

! DEPENDS ON: sl_vector_w
        Call SL_Vector_w(                                               &
     &                    w, Yu, Yv, Yw,                                &
     &                    eta_theta_levels, eta_rho_levels,             &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    r_at_u, r_at_v,                               &
     &                    depart_lambda_w, depart_phi_w,                &
     &                    depart_r_w,                                   &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    sin_theta_latitude,                           &
     &                    cos_theta_longitude,                          &
     &                    sin_theta_longitude,                          &
     &                    delta_lambda, delta_phi,                      &
     &                    glambda_p, phi_p, glambda_u, phi_v,           &
     &                    gdlambda_p, dphi_p, gdlambda_u, dphi_v,       &
     &                    grecip_dlamp, recip_dphip,                    &
     &                    grecip_dlamu, recip_dphiv,                    &
     &                    lambda_p_rm, lambda_p_rp,                     &
     &                    lambda_u_rm, lambda_u_rp,                     &
     &                    phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,       &
     &                    recip_lambda_p_m, recip_lambda_p_0,           &
     &                    recip_lambda_p_p, recip_lambda_p_p2,          &
     &                    recip_lambda_u_m, recip_lambda_u_0,           &
     &                    recip_lambda_u_p, recip_lambda_u_p2,          &
     &                    recip_phi_p_m, recip_phi_p_0,                 &
     &                    recip_phi_p_p, recip_phi_p_p2,                &
     &                    recip_phi_v_m, recip_phi_v_0,                 &
     &                    recip_phi_v_p, recip_phi_v_p2,                &
     &                    Base_lambda, base_phi,                        &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &     
     &                    n_Y_arrays, n_Yw_arrays,                      &
     &                    n_Yd_arrays, n_Ydw_arrays,                    &
     &                    timestep, alpha_3, alpha_4,                   &
     &                    row_length, rows, n_rows, model_levels,       &
     &                    high_order_scheme, monotone_scheme,           &
     &                    model_domain, L_2d_sl_geometry, L_high,       &
     &                    L_mono, L_conserv, check_bottom_levels,       &
     &                    interp_vertical_search_tol,                   &
     &                    first_constant_r_rho_level,                   &
     &                    L_trivial_trigs,                              &
     &                    me, n_proc, n_procx, n_procy,                 &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    global_row_length, global_rows,               &
     &                    l_datastart, at_extremity, g_i_pe, g_j_pe,    &
     &                    l_2dcomm, size_2dcomm,                        &
     &                    group_2dcomm, proc_row_group, proc_col_group, &
     &                    proc_all_group, L_regular,                    &
     &                    L_sl_halo_reprod, L_lbc_old,                  &
     &                    L_new_tdisc, CycleNo, R_w,                    &
     &                    Error_code)

      End If

! End of routine.
      IF (lhook) CALL dr_hook('SL_FULL_WIND',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SL_Full_Wind

