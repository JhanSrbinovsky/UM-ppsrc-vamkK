! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine SL_Vector_w


SUBROUTINE sl_vector_w(                                                        &
  w, yu, yv, yw,                                                               &
  eta_theta_levels, eta_rho_levels,                                            &
  r_theta_levels, r_rho_levels,                                                &
  r_at_u, r_at_v,                                                              &
  depart_lambda_w, depart_phi_w,                                               &
  depart_r_w,                                                                  &
  cos_theta_latitude, sec_theta_latitude,                                      &
  sin_theta_latitude,                                                          &
  cos_theta_longitude,                                                         &
  sin_theta_longitude,                                                         &
  delta_lambda, delta_phi,                                                     &
  glambda_p, phi_p, glambda_u, phi_v,                                          &
  gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                      &
  grecip_dlamp, recip_dphip,                                                   &
  grecip_dlamu, recip_dphiv,                                                   &
  lambda_p_rm, lambda_p_rp,                                                    &
  lambda_u_rm, lambda_u_rp,                                                    &
  phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,                                      &
  recip_lambda_p_m, recip_lambda_p_0,                                          &
  recip_lambda_p_p, recip_lambda_p_p2,                                         &
  recip_lambda_u_m, recip_lambda_u_0,                                          &
  recip_lambda_u_p, recip_lambda_u_p2,                                         &
  recip_phi_p_m, recip_phi_p_0,                                                &
  recip_phi_p_p, recip_phi_p_p2,                                               &
  recip_phi_v_m, recip_phi_v_0,                                                &
  recip_phi_v_p, recip_phi_v_p2,                                               &
  base_lambda, base_phi,                                                       &
  recip_dlam, recip_dphi, max_look,                                            &
  look_lam, look_phi, halo_lam, halo_phi,                                      &
  n_y_arrays, n_yw_arrays,                                                     &
  n_yd_arrays, n_ydw_arrays,                                                   &
  timestep, alpha_3, alpha_4,                                                  &
  row_length, rows, n_rows, model_levels,                                      &
  high_order_scheme, monotone_scheme,                                          &
  model_domain, l_2d_sl_geometry,                                              &
  l_high, l_mono, l_conserv,                                                   &
  check_bottom_levels,                                                         &
  interp_vertical_search_tol,                                                  &
  first_constant_r_rho_level,                                                  &
  l_trivial_trigs,                                                             &
  me, n_proc, n_procx, n_procy,                                                &
  off_x, off_y, halo_i, halo_j,                                                &
  global_row_length, global_rows,                                              &
  l_datastart,  at_extremity,                                                  &
  g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,                                       &
  group_2dcomm, proc_row_group,                                                &
  proc_col_group, proc_all_group,                                              &
  l_regular, l_sl_halo_reprod, l_lbc_old,                                      &
  l_new_tdisc, cycleno, r_w, error_code)

! Purpose:
!          Performs vector semi-Lagrangian integration of values at
!          w points given the forcing functions for the first advection
!          step.

! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson

!          and

!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.



USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER                                                                        &
  row_length                                                                   &
                                ! number of point on a row.
  , rows                                                                       &
                                ! number of rows.
  , n_rows                                                                     &
                                ! number of v rows.
  , model_levels                                                               &
                                ! number of model levels.
  , me                                                                         &
                                ! My processor number
  , n_proc                                                                     &
                                ! Total number of processors
  , n_procx                                                                    &
                                ! Number of processors in longitude
  , n_procy                                                                    &
                                ! Number of processors in latitude
  , halo_i                                                                     &
                                ! Size of halo in i.
  , halo_j                                                                     &
                                ! Size of halo in j.
  , off_x                                                                      &
                                ! Size of small halo in i
  , off_y                                                                      &
                                ! Size of small halo in j.
  , l_datastart(3)                                                             &
                                ! First gridpoints held by this processor.
  , proc_row_group                                                             &
                                ! Group id for processors on the same row
  , proc_col_group                                                             &
                                ! Group id for processors on the same column
  , proc_all_group                                                             &
                                ! Group id for all processors
  , max_look                                                                   &
                                ! max size of look-up arrays for searches
  , global_row_length                                                          &
                                ! global number of points on a row
  , global_rows                                                                &
                                ! global number of rows
  , g_i_pe(1-halo_i:global_row_length+halo_i)                                  &
                                ! processor on my
                                ! processor-row holding a given value in i direction
  , g_j_pe(1-halo_j:global_rows+halo_j)                                        &
                                ! processor on my
                                ! processor-column holding a given value in j direction
  , n_y_arrays                                                                 &
                                ! = 1 for global, 3 for LAM
  , n_yw_arrays                                                                &
                                ! = 1 for global, 2 for LAM
  , n_yd_arrays                                                                &
                                ! = 1 for global, 3 for LAM
  , n_ydw_arrays                                                               &
                                ! = 1 for global, 2 for LAM
  , cycleno

LOGICAL                                                                        &
  l_sl_halo_reprod                                                             &
                                ! if true then sl code bit repoducible with
                                ! any sensible halo size
  , l_lbc_old                                                                  &
                                ! false for new lbc treatment
  , l_new_tdisc

LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
INTEGER :: size_int      ! error check size for comms on demand

INTEGER                                                                        &
  high_order_scheme                                                            &
                                ! a code saying which high order scheme to
                                ! use. 1 = tensor tri-cubic lagrange order
                                ! (j,i,k) no other options available at
                                ! present.
  , monotone_scheme                                                            &
                                ! a code saying which monotone scheme to use.
                                ! 1 = tri-linear order (j,i,k)
                                ! no other options available at present.
  , interp_vertical_search_tol                                                 &
                                ! used in interpolation code.
  , check_bottom_levels ! used in interpolation code, and is
! the number of levels to check to see
! if the departure point lies inside the
! orography.

INTEGER                                                                        &
  model_domain                                                                 &
                                ! holds integer code for model domain
  , first_constant_r_rho_level

LOGICAL                                                                        &
  l_2d_sl_geometry                                                             &
                                ! True, then only perform vector co-ordinate
                                !       geometry in 2d.
  , l_high                                                                     &
                                ! True, if high order interpolation required.
  , l_mono                                                                     &
                                ! True, if interpolation required to be monotone.
  , l_conserv                                                                  &
                                ! True, if interpolation to be monotone and
                                !       conservative.
  , l_regular

LOGICAL :: l_trivial_trigs ! True if trivial_trigs (Cartesian grid)

REAL                                                                           &
  delta_lambda                                                                 &
                                ! grid-length in lambda direction
  , delta_phi                                                                  &
                                ! grid-length in phi direction
  , base_lambda                                                                &
  , base_phi                                                                   &
  , recip_dlam                                                                 &
  , recip_dphi                                                                 &
  , alpha_3                                                                    &
  , alpha_4                                                                    &
  , timestep

! look-up table halos
INTEGER                                                                        &
  halo_lam                                                                     &
  , halo_phi

!VarRes horizontal co-ordinate look-up table
INTEGER                                                                        &
  look_lam(1-halo_lam : max_look-halo_lam)                                     &
  , look_phi(1-halo_phi : max_look-halo_phi)

!VarRes horizontal co-ordinate information etc.
REAL                                                                           &
  glambda_p( 1-halo_i : global_row_length+halo_i)                              &
  , glambda_u( 1-halo_i : global_row_length+halo_i)                            &
  , phi_p    ( 1-halo_i : row_length + halo_i                                  &
  ,            1-halo_j : rows + halo_j )                                      &
  , phi_v    ( 1-halo_i : row_length + halo_i                                  &
  ,            1-halo_j : n_rows+halo_j )                                      &
  , gdlambda_p(1-halo_i : global_row_length+halo_i)                            &
  , gdlambda_u(1-halo_i : global_row_length+halo_i)                            &
  , dphi_p  ( 1-halo_i : row_length + halo_i                                   &
  ,           1-halo_j : rows+halo_j )                                         &
  , dphi_v  ( 1-halo_i : row_length + halo_i                                   &
  ,           1-halo_j : n_rows+halo_j )                                       &
  , grecip_dlamp(1-halo_i : global_row_length + halo_i)                        &
  , grecip_dlamu(1-halo_i : global_row_length + halo_i)                        &
  , recip_dphip( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows+halo_j )                                      &
  , recip_dphiv( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : n_rows+halo_j )                                    &
  , lambda_p_rm (1-halo_i : row_length + halo_i)                               &
  , lambda_p_rp (1-halo_i : row_length + halo_i)                               &
  , lambda_u_rm (1-halo_i : row_length + halo_i)                               &
  , lambda_u_rp (1-halo_i : row_length + halo_i)                               &
  , phi_p_rm   ( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows+halo_j )                                      &
  , phi_p_rp   ( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : rows+halo_j )                                      &
  , phi_v_rm   ( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : n_rows+halo_j )                                    &
  , phi_v_rp   ( 1-halo_i : row_length + halo_i                                &
  ,              1-halo_j : n_rows+halo_j )                                    &
  , recip_lambda_p_m(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_p_0(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_p_p(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_p_p2(1-halo_i : row_length+halo_i)                            &
  , recip_lambda_u_m(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_u_0(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_u_p(1-halo_i : row_length+halo_i)                             &
  , recip_lambda_u_p2(1-halo_i : row_length+halo_i)                            &
  , recip_phi_p_m ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : rows+halo_j )                                   &
  , recip_phi_p_0 ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : rows+halo_j )                                   &
  , recip_phi_p_p ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : rows+halo_j )                                   &
  , recip_phi_p_p2( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : rows+halo_j )                                   &
  , recip_phi_v_m ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : n_rows+halo_j )                                 &
  , recip_phi_v_0 ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : n_rows+halo_j )                                 &
  , recip_phi_v_p ( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : n_rows+halo_j )                                 &
  , recip_phi_v_p2( 1-halo_i : row_length + halo_i                             &
  ,                 1-halo_j : n_rows+halo_j )

REAL                                                                           &
                                ! trigonometric functions
  cos_theta_latitude (1-off_x:row_length+off_x,                                &
  1-off_y:rows+off_y)                                                          &
  , sec_theta_latitude (1-off_x:row_length+off_x,                              &
  1-off_y:rows+off_y)                                                          &
  , sin_theta_latitude (row_length, rows)                                      &
  , cos_theta_longitude (row_length, rows)                                     &
  , sin_theta_longitude (row_length, rows)                                     &
  , cos_u_longitude (row_length, rows)                                         &
  , sin_u_longitude (row_length, rows)

REAL                                                                           &
                                ! primary model variables
  w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels)

REAL                                                                           &
                                ! co-ordinates of departure points for w.
  depart_lambda_w(row_length, rows, model_levels)                              &
  , depart_phi_w(row_length, rows, model_levels)                               &
  , depart_r_w(row_length, rows, model_levels)

REAL                                                                           &
  yu (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                        &
  model_levels, n_y_arrays)                                                    &
  , yv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,                    &
  model_levels, n_y_arrays)                                                    &
  , yw (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                      &
  0:model_levels, n_yw_arrays)

REAL                                                                           &
                                ! model level arrays
  eta_theta_levels(0:model_levels)                                             &
  , eta_rho_levels(model_levels)

REAL                                                                           &
                                ! vertical co-ordinate arrays
  r_theta_levels (1-halo_i:row_length+halo_i,                                  &
  1-halo_j:rows+halo_j, 0:model_levels)                                        &
  , r_rho_levels (1-halo_i:row_length+halo_i,                                  &
  1-halo_j:rows+halo_j, model_levels)                                          &
  , r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                  &
  model_levels)                                                                &
  , r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,                &
  model_levels)

INTEGER                                                                        &
  i_out (row_length, rows, model_levels)                                       &
  , j_out (row_length, rows, model_levels)

REAL                                                                           &
  weight_lambda (row_length, rows, model_levels)                               &
  , weight_phi    (row_length, rows, model_levels)

INTEGER                                                                        &
  neighbour(4)         ! Array with the Ids of the four neighbours
! in the horizontal plane

LOGICAL                                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
! south, east or west of the processor grid

! Parameters

! Arguments with Intent INOUT. ie: Input variables changed on output.

REAL                                                                           &
                                ! See WP154 and WP162 for definitions.
  r_w (row_length, rows, model_levels)

! Arguments with Intent OUT. ie: Output variables.

INTEGER                                                                        &
  error_code     ! Non-zero on exit if error detected.

! Local Variables.

INTEGER                                                                        &
  i, j, k, gi                                                                  &
                                ! Loop indices
  , imin, imax, jmin, jmax

INTEGER                                                                        &
  itype                                                                        &
                                ! a code saying what points the grid points are
                                ! on.
  , number_of_inputs !the number of fields to interpolate in any one
! call to the interpolation routine.

LOGICAL                                                                        &
  l_vector       ! True, if data is a horizontal vector component,
! False, then data is a scalar.

REAL                                                                           &
  cos_lambda_ad                                                                &
                                ! cosine(lambda_a - lambda_d)
  , sin_lambda_ad                                                              &
                                ! sine(lambda_a - lambda_d)
  , cos_phi_d                                                                  &
  , sin_phi_d                                                                  &
  , trig_u                                                                     &
  , trig_v                                                                     &
  , trig_w                                                                     &
  , dummy                                                                      &
  , delta_r        ! dummy array

REAL                                                                           &
  yu_d (row_length, rows, model_levels-1, n_yd_arrays)                         &
  , yv_d (row_length, rows, model_levels-1, n_yd_arrays)                       &
  , yw_d (row_length, rows, model_levels-1, n_ydw_arrays)                      &
  , lambda_a(row_length)

REAL                                                                           &
  tmp1          ! temporary variable to handle int -> real

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle



! ----------------------------------------------------------------------
! Section 1.  Call interpolation routine for each term.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SL_VECTOR_W',zhook_in,zhook_handle)

IF (l_2dcomm) THEN
  size_int=row_length * rows *(model_levels-1)
ELSE
  size_int=global_row_length*(model_levels-1)
END IF

! ----------------------------------------------------------------------
! Section 1.1 Call interpolation routine for Yu.
! ----------------------------------------------------------------------

! DEPENDS ON: calc_index
CALL calc_index(                                                               &
  row_length, rows, model_levels-1,                                            &
  delta_lambda, delta_phi,                                                     &
  base_lambda, base_phi,                                                       &
  glambda_p, phi_p, grecip_dlamp, recip_dphip,                                 &
  recip_dlam, recip_dphi, max_look,                                            &
  look_lam, look_phi, halo_lam, halo_phi,                                      &
  l_regular, depart_lambda_w, depart_phi_w,                                    &
  halo_i, halo_j,                                                              &
  global_row_length,                                                           &
  row_length, rows, l_datastart,                                               &
  i_out, j_out,                                                                &
  weight_lambda, weight_phi)


IF (error_code  ==  0 ) THEN

  l_vector = .TRUE.

  IF (.NOT. l_2d_sl_geometry .AND.                                             &
    (model_domain  ==  mt_global .OR.                                          &
    model_domain  ==  mt_cyclic_lam .OR.                                       &
    model_domain  ==  mt_bi_cyclic_lam) ) THEN

    itype = 1 ! data is at u points.

    number_of_inputs = 1

    ! DEPENDS ON: interpolation
    CALL interpolation(                                                        &
      yu, dummy, dummy,                                                        &
      eta_rho_levels,                                                          &
      r_at_u, delta_r, r_rho_levels, itype,                                    &
      number_of_inputs,                                                        &
      check_bottom_levels,                                                     &
      interp_vertical_search_tol,                                              &
      first_constant_r_rho_level,                                              &
      row_length, rows, model_levels,                                          &
      rows,                                                                    &
      row_length, rows, model_levels-1,                                        &
      delta_lambda, delta_phi,                                                 &
      base_lambda, base_phi,                                                   &
      glambda_u, phi_v,                                                        &
      gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                  &
      grecip_dlamu, recip_dphiv,                                               &
      lambda_u_rm, lambda_u_rp,                                                &
      phi_p_rm, phi_p_rp,                                                      &
      recip_lambda_u_m, recip_lambda_u_0,                                      &
      recip_lambda_u_p, recip_lambda_u_p2,                                     &
      recip_phi_p_m, recip_phi_p_0,                                            &
      recip_phi_p_p, recip_phi_p_p2,                                           &
      i_out, j_out,                                                            &
      weight_lambda, weight_phi,                                               &
      high_order_scheme, monotone_scheme,                                      &
      cos_theta_latitude, l_regular,                                           &
      l_vector, model_domain, l_high, l_mono,                                  &
      .FALSE.,                                                                 &
      depart_r_w, depart_lambda_w, depart_phi_w,                               &
      me, n_proc, n_procx, n_procy,                                            &
      halo_i, halo_j,                                                          &
      global_row_length, global_rows,                                          &
      row_length, rows, n_rows, rows,                                          &
      l_datastart, at_extremity, g_i_pe,                                       &
      g_j_pe, l_2dcomm, size_2dcomm,                                           &
      group_2dcomm, size_int, proc_all_group,                                  &
      proc_row_group, proc_col_group, 1, 0, 0,                                 &
      l_sl_halo_reprod, off_x,off_y,                                           &
      yu_d, dummy, dummy, error_code)

    ! ----------------------------------------------------------------------
    ! Section 2.2 Call interpolation routine for Yv.
    ! ----------------------------------------------------------------------

    itype = 2 ! data at v points.

    number_of_inputs = 1

    ! DEPENDS ON: interpolation
    CALL interpolation(                                                        &
      yv, dummy, dummy,                                                        &
      eta_rho_levels,                                                          &
      r_at_v, delta_r, r_rho_levels, itype,                                    &
      number_of_inputs,                                                        &
      check_bottom_levels,                                                     &
      interp_vertical_search_tol,                                              &
      first_constant_r_rho_level,                                              &
      row_length, n_rows, model_levels,                                        &
      rows,                                                                    &
      row_length, rows, model_levels-1,                                        &
      delta_lambda, delta_phi,                                                 &
      base_lambda, base_phi,                                                   &
      glambda_u, phi_v,                                                        &
      gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                  &
      grecip_dlamu, recip_dphiv,                                               &
      lambda_p_rm, lambda_p_rp,                                                &
      phi_v_rm, phi_v_rp,                                                      &
      recip_lambda_p_m, recip_lambda_p_0,                                      &
      recip_lambda_p_p, recip_lambda_p_p2,                                     &
      recip_phi_v_m, recip_phi_v_0,                                            &
      recip_phi_v_p, recip_phi_v_p2,                                           &
      i_out, j_out,                                                            &
      weight_lambda, weight_phi,                                               &
      high_order_scheme, monotone_scheme,                                      &
      cos_theta_latitude, l_regular,                                           &
      l_vector, model_domain, l_high, l_mono,                                  &
      .FALSE.,                                                                 &
      depart_r_w, depart_lambda_w, depart_phi_w,                               &
      me, n_proc, n_procx, n_procy,                                            &
      halo_i, halo_j,                                                          &
      global_row_length, global_rows,                                          &
      row_length, rows, n_rows, n_rows,                                        &
      l_datastart, at_extremity, g_i_pe,                                       &
      g_j_pe, l_2dcomm, size_2dcomm,                                           &
      group_2dcomm, size_int, proc_all_group,                                  &
      proc_row_group, proc_col_group, 1, 0, 0,                                 &
      l_sl_halo_reprod, off_x,off_y,                                           &
      yv_d, dummy, dummy, error_code)

  END IF ! on 2d geometry and domain

END IF !  Error_Code == 0

! ----------------------------------------------------------------------
! Section 2.3 Call interpolation routine for Yw.
! ----------------------------------------------------------------------

IF (error_code  ==  0 ) THEN

  itype = 3 ! data at w points.
  l_vector = .FALSE.

  IF (model_domain  ==  mt_global .OR.                                         &
    model_domain  ==  mt_cyclic_lam .OR.                                       &
    model_domain  ==  mt_bi_cyclic_lam) THEN

    number_of_inputs = 1

    ! DEPENDS ON: interpolation
    CALL interpolation(                                                        &
      yw, dummy, dummy,                                                        &
      eta_theta_levels,                                                        &
      r_theta_levels, delta_r, r_theta_levels,                                 &
      itype, number_of_inputs,                                                 &
      check_bottom_levels,                                                     &
      interp_vertical_search_tol,                                              &
      first_constant_r_rho_level+1,                                            &
      row_length, rows, model_levels+1,                                        &
      rows,                                                                    &
      row_length, rows, model_levels-1,                                        &
      delta_lambda, delta_phi,                                                 &
      base_lambda, base_phi,                                                   &
      glambda_u, phi_v,                                                        &
      gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                  &
      grecip_dlamu, recip_dphiv,                                               &
      lambda_p_rm, lambda_p_rp,                                                &
      phi_p_rm, phi_p_rp,                                                      &
      recip_lambda_p_m, recip_lambda_p_0,                                      &
      recip_lambda_p_p, recip_lambda_p_p2,                                     &
      recip_phi_p_m, recip_phi_p_0,                                            &
      recip_phi_p_p, recip_phi_p_p2,                                           &
      i_out, j_out,                                                            &
      weight_lambda, weight_phi,                                               &
      high_order_scheme, monotone_scheme,                                      &
      cos_theta_latitude, l_regular,                                           &
      l_vector, model_domain, l_high, l_mono,                                  &
      l_conserv,                                                               &
      depart_r_w, depart_lambda_w, depart_phi_w,                               &
      me, n_proc, n_procx, n_procy,                                            &
      halo_i, halo_j,                                                          &
      global_row_length, global_rows,                                          &
      row_length, rows, n_rows, rows,                                          &
      l_datastart, at_extremity, g_i_pe,                                       &
      g_j_pe, l_2dcomm, size_2dcomm,                                           &
      group_2dcomm, size_int, proc_all_group,                                  &
      proc_row_group, proc_col_group, 1, 0, 0,                                 &
      l_sl_halo_reprod, off_x,off_y,                                           &
      yw_d, dummy, dummy, error_code)

  ELSE

    number_of_inputs = 2

    ! DEPENDS ON: interpolation
    CALL interpolation(                                                        &
      yw(1-halo_i,1-halo_j,0,2),                                               &
      yw(1-halo_i,1-halo_j,0,1), dummy,                                        &
      eta_theta_levels,                                                        &
      r_theta_levels, delta_r, r_theta_levels,                                 &
      itype, number_of_inputs,                                                 &
      check_bottom_levels,                                                     &
      interp_vertical_search_tol,                                              &
      first_constant_r_rho_level+1,                                            &
      row_length, rows, model_levels+1,                                        &
      rows,                                                                    &
      row_length, rows, model_levels-1,                                        &
      delta_lambda, delta_phi,                                                 &
      base_lambda, base_phi,                                                   &
      glambda_u, phi_v,                                                        &
      gdlambda_p, dphi_p, gdlambda_u, dphi_v,                                  &
      grecip_dlamu, recip_dphiv,                                               &
      lambda_p_rm, lambda_p_rp,                                                &
      phi_p_rm, phi_p_rp,                                                      &
      recip_lambda_p_m, recip_lambda_p_0,                                      &
      recip_lambda_p_p, recip_lambda_p_p2,                                     &
      recip_phi_p_m, recip_phi_p_0,                                            &
      recip_phi_p_p, recip_phi_p_p2,                                           &
      i_out, j_out,                                                            &
      weight_lambda, weight_phi,                                               &
      high_order_scheme, monotone_scheme,                                      &
      cos_theta_latitude, l_regular,                                           &
      l_vector, model_domain, l_high, l_mono,                                  &
      l_conserv,                                                               &
      depart_r_w, depart_lambda_w, depart_phi_w,                               &
      me, n_proc, n_procx, n_procy,                                            &
      halo_i, halo_j,                                                          &
      global_row_length, global_rows,                                          &
      row_length, rows, n_rows, rows,                                          &
      l_datastart, at_extremity, g_i_pe,                                       &
      g_j_pe, l_2dcomm, size_2dcomm,                                           &
      group_2dcomm, size_int, proc_all_group,                                  &
      proc_row_group, proc_col_group, 1, 0, 0,                                 &
      l_sl_halo_reprod, off_x,off_y,                                           &
      yw_d(1,1,1,1), yw_d(1,1,1,2),                                            &
      dummy, error_code)

  END IF  ! model domain

END IF   ! Error_Code == 0

! ----------------------------------------------------------------------
! Section 3.  Form R_w from the interpolated terms using
!             expressions given in WP 162.
! ----------------------------------------------------------------------

IF ( l_regular ) THEN
  DO i = 1, row_length
    gi = l_datastart(1) + i - 1
    tmp1 = float(gi-1)
    lambda_a(i) = tmp1 * delta_lambda
  END DO
ELSE ! variable resolution
  DO i = 1, row_length
    gi = l_datastart(1) + i - 1
    lambda_a(i) = glambda_p(gi) - base_lambda
  END DO
END IF ! L_regular

IF (error_code  ==  0 ) THEN

  IF (model_domain  ==  mt_global .OR.                                         &
    model_domain  ==  mt_cyclic_lam .OR.                                       &
    model_domain  ==  mt_bi_cyclic_lam) THEN

    IF (l_2d_sl_geometry .OR. l_trivial_trigs) THEN

      DO k = 1, model_levels-1
        DO j = 1, rows
          DO i = 1, row_length
            r_w(i,j,k) =  yw_d(i,j,k,1) + r_w(i,j,k)
          END DO
        END DO
      END DO

    ELSE

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(model_levels, rows,   &
!$OMP& row_length, lambda_a, depart_lambda_w, depart_phi_w, cos_theta_latitude,&
!$OMP& sin_theta_latitude, r_w, yu_d, yv_d, yw_d) PRIVATE(cos_lambda_ad,       &
!$OMP& sin_lambda_ad, cos_phi_d, sin_phi_d, trig_u, trig_w, trig_v, i, j, k)
      DO k = 1, model_levels-1
        DO j = 1, rows
          DO i = 1, row_length
            cos_lambda_ad = COS( lambda_a(i) -                                 &
              depart_lambda_w(i,j,k) )
            sin_lambda_ad = SIN( lambda_a(i) -                                 &
              depart_lambda_w(i,j,k) )
            cos_phi_d = COS ( depart_phi_w(i,j,k))
            sin_phi_d = SIN ( depart_phi_w(i,j,k))
            trig_u = sin_lambda_ad * cos_theta_latitude(i,j)
            trig_v = sin_theta_latitude(i,j) * cos_phi_d -                     &
              cos_theta_latitude(i,j) * sin_phi_d *                            &
              cos_lambda_ad
            trig_w = sin_theta_latitude(i,j) * sin_phi_d +                     &
              cos_theta_latitude(i,j) * cos_phi_d *                            &
              cos_lambda_ad
            r_w(i,j,k) = trig_u * yu_d(i,j,k,1)                                &
              + trig_v * yv_d(i,j,k,1)                                         &
              + trig_w * yw_d(i,j,k,1)                                         &
              + r_w(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO 

    END IF

  ELSE IF (model_domain  ==  mt_lam ) THEN
    ! NB: code is the same whether 2d geometry option on or off

    IF ( cycleno == 1 .OR. .NOT. l_new_tdisc ) THEN

      ! internal area uses real alphas
      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            r_w(i,j,k) = (yw_d(i,j,k,1) +                                      &
              yw_d(i,j,k,2) * (1. - alpha_4) )                                 &
              + alpha_4 * yw(i,j,k,1)                                          &
              - w(i,j,k)
          END DO
        END DO
      END DO

      IF ( l_lbc_old ) THEN
        ! boundary area uses alphas = 1
        ! southern area
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels - 1
            DO j = 1, halo_j
              DO i = 1, row_length
                r_w(i,j,k) =  yw_d(i,j,k,1)                                    &
                  + yw(i,j,k,1)                                                &
                  - w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Western Area
        IF (at_extremity(pwest)) THEN
          DO k = 1, model_levels - 1
            DO j = 1, rows
              DO i = 1, halo_i
                r_w(i,j,k) =  yw_d(i,j,k,1)                                    &
                  + yw(i,j,k,1)                                                &
                  - w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Eastern Area
        IF (at_extremity(peast)) THEN
          DO k = 1, model_levels - 1
            DO j = 1, rows
              DO i = row_length-halo_i+1, row_length
                r_w(i,j,k) =  yw_d(i,j,k,1)                                    &
                  + yw(i,j,k,1)                                                &
                  - w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Northern Area
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels - 1
            DO j = rows-halo_j+1, rows
              DO i = 1, row_length
                r_w(i,j,k) =  yw_d(i,j,k,1)                                    &
                  + yw(i,j,k,1)                                                &
                  - w(i,j,k)
              END DO
            END DO
          END DO
        END IF

      END IF ! L_lbc_old

    ELSE ! if CycleNo > 1 and L_new_tdisc

      ! First do the boundaries to set appropriate loop indices
      ! thus avoiding overwriting R_u or using temp storage

      ! Loop indices for R_u update when processor is away from boundaries

      imin = 1
      imax = row_length
      jmin = 1
      jmax = row_length

      IF ( l_lbc_old ) THEN
        ! boundary area uses alphas = 1
        ! southern area
        IF (at_extremity(psouth)) THEN
          DO k = 1, model_levels - 1
            DO j = 1, halo_j
              DO i = 1, row_length
                r_w(i,j,k) =  yw_d(i,j,k,1) + r_w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Western Area
        IF (at_extremity(pwest)) THEN
          DO k = 1, model_levels - 1
            ! use jmin to avoid overwriting S boundary
            DO j = jmin, rows
              DO i = 1, halo_i
                r_w(i,j,k) =  yw_d(i,j,k,1) + r_w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Eastern Area
        IF (at_extremity(peast)) THEN
          DO k = 1, model_levels - 1
            ! use jmin to avoid overwriting S boundary
            DO j = jmin, rows
              DO i = row_length-halo_i+1, row_length
                r_w(i,j,k) =  yw_d(i,j,k,1) + r_w(i,j,k)
              END DO
            END DO
          END DO
        END IF
        ! Northern Area
        IF (at_extremity(pnorth)) THEN
          DO k = 1, model_levels - 1
            DO j = rows-halo_j+1, rows
              ! use imin, imax to avoid overwriting E and W boundary
              DO i = imin, imax
                r_w(i,j,k) =  yw_d(i,j,k,1) + r_w(i,j,k)
              END DO
            END DO
          END DO
        END IF

      END IF ! L_lbc_old

      ! internal area uses real alphas
      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            r_w(i,j,k) = (yw_d(i,j,k,1) +                                      &
              yw_d(i,j,k,2) * (1. - alpha_4) ) +                               &
              r_w(i,j,k)
          END DO
        END DO
      END DO

    END IF ! CycleNo ==1 ...

  END IF ! on model domain

  k = model_levels
  DO j = 1, rows
    DO i = 1, row_length
      r_w(i,j,k) = 0.
    END DO
  END DO

END IF

! End of routine.
IF (lhook) CALL dr_hook('SL_VECTOR_W',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sl_vector_w

