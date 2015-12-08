! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Departure_Point.
!

      Subroutine Departure_Point(                                       &
     &                          type, timestep, u_adv, v_adv, w_adv,    &
     &                          eta_rho_levels, eta_theta_levels,       &
     &                          r_rho_levels, r_theta_levels,           &
     &                          r_at_u, r_at_v,                         &
     &                          row_length, rows, n_rows, model_levels, &
     &                          delta_lambda, delta_phi,                &
     &                          glambda_p, phi_p, glambda_u, phi_v,     &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamp, recip_dphip,              &
     &                          grecip_dlamu, recip_dphiv,              &
     &                    wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v, &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          lambda_u_rm, lambda_u_rp,               &
     &                          phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp, &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_lambda_u_m, recip_lambda_u_0,     &
     &                          recip_lambda_u_p, recip_lambda_u_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          recip_phi_v_m, recip_phi_v_0,           &
     &                          recip_phi_v_p, recip_phi_v_p2,          &
     &                          base_lambda, base_phi,                  &
     &                          lambda_p_end, phi_p_end, dlambda_p_end, &
     &                          dphi_p_end, dphi_v_end,                 &
     &                          recip_dlam, recip_dphi, max_look,       &
     &                          look_lam, look_phi, halo_lam, halo_phi, &
     &                          cos_theta_latitude, sec_theta_latitude, &
     &                          sin_theta_latitude, tan_theta_latitude, &
     &                          cos_v_latitude, sec_v_latitude,         &
     &                          sin_v_latitude, tan_v_latitude,         &
     &                          sin_theta_longitude,                    &
     &                          cos_theta_longitude,                    &
     &                          model_domain, Depart_scheme,            &
     &                          Depart_order, L_high, L_mono,           &
     &                          high_order_scheme, monotone_scheme,     &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          L_2d_sl_geometry, L_regular,            &
     &                          first_constant_r_rho_level,             &
     &                          LAM_max_cfl, rows_depart,               &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          off_x, off_y, halo_i, halo_j, datastart,&
     &                          global_row_length, global_rows,         &
     &                          g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,  &
     &                          group_2dcomm, proc_row_group,           &
     &                          proc_col_group, proc_all_group,         &
     &                          at_extremity,                           &
     &                          L_sl_halo_reprod,                       &
     &                          depart_lambda, depart_phi, depart_r)

! Purpose:
!          Finds departure point for trajectory defined by the
!          given wind field.
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


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                       ! Dimension of u_adv array in i direction.
     &, rows                                                            &
                       ! Dimension of u_adv array in j direction.
     &, n_rows                                                          &
                       ! Dimension of v_adv array in j direction.
     &, model_levels                                                    &
                       ! Dimension of u_adv array in k direction.
     &, rows_depart                                                     &
                       ! Dimension of depart arrays in j direction.
     &, me                                                              &
                       ! My processor number
     &, n_proc                                                          &
                       ! Total number of processors
     &, n_procx                                                         &
                       ! Number of processors in longitude
     &, n_procy                                                         &
                       ! Number of processors in latitude
     &, halo_i                                                          &
                       ! Size of halo in i direction.
     &, halo_j                                                          &
                       ! Size of halo in j direction.
     &, off_x                                                           &
     &, off_y                                                           &
     &, datastart(3)                                                    &
                       ! First gridpoints held by this processor.
     &, max_look                                                        &
                       ! max size of look-up arrays for searches
     &, global_row_length                                               &
                            ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of points in a column
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, proc_all_group                                                  &
                       ! Group id for all processors
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                                                  ! processor on my
                   ! processor-row holding a given value in i direction
     &, g_j_pe(1-halo_j:global_rows      +halo_j)                       &
                                                  ! processor on my
                ! processor-column holding a given value in j direction
     &, LAM_max_cfl(2)

      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_2d_sl_geometry                                                &
     &, L_regular      ! False if variable resolution


      Integer                                                           &
     &  type                                                            &
                    ! Defines via an integer code the nature of the
                    ! grid for which the departure points are required.
                    ! The codes are given in
                    ! terms of a primary variable that would be held
                    ! at that point and are u=1, v=2, w=3, theta=4,
                    ! rho=5.
     &, Depart_scheme                                                   &
                      ! code saying which departure point scheme to
                      ! use.
     &, Depart_order                                                    &
                      ! for the chosen departure point scheme how
                      ! many iterations/terms to use.
     &, high_order_scheme                                               &
                          ! code choosing high order interpolation
                          ! scheme used in Robert routine
     &, monotone_scheme                                                 &
                          ! code choosing monotone interpolation
                          ! scheme used in Robert routine
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels                                             &
                            ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.
     &, first_constant_r_rho_level

      Logical                                                           &
     &  L_high                                                          &
                      ! True if high order interpolation scheme to be
                      ! used in Robert scheme
     &, L_mono        ! True if monotone interpolation scheme to be
                      ! used in Robert scheme

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  delta_lambda                                                    &
                         ! holds spacing between points in the i
                         ! direction.
     &, delta_phi                                                       &
                         ! holds spacing between points in the j
                         ! direction.
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
       Integer                                                          &
     &  halo_lam                                                        &
     &, halo_phi

!VarRes horizontal co-ordinate look-up table
      Integer                                                           &
     &  look_lam(1-halo_lam : max_look-halo_lam)                        &
     &, look_phi(1-halo_phi : max_look-halo_phi)

!VarRes horizontal co-ordinate information
     Real                                                               &
     &  glambda_p( 1-halo_i : global_row_length+halo_i )                &
     &, glambda_u( 1-halo_i : global_row_length+halo_i )                &
     &, phi_p    ( 1-halo_i : row_length + halo_i,                      &
     &             1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i,                      &
     &             1-halo_j : n_rows+halo_j )                           &
     &, gdlambda_p( 1-halo_i : global_row_length+halo_i)                &
     &, gdlambda_u( 1-halo_i : global_row_length+halo_i)                &
     &, dphi_p   ( 1-halo_i : row_length + halo_i,                      &
     &             1-halo_j : rows + halo_j )                           &
     &, dphi_v   ( 1-halo_i : row_length + halo_i,                      &
     &             1-halo_j : n_rows+halo_j )                           &
     &, grecip_dlamp(1-halo_i : global_row_length + halo_i)             &
     &, grecip_dlamu(1-halo_i : global_row_length + halo_i)             &
     &, recip_dphip( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : rows + halo_j )                         &
     &, recip_dphiv( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : n_rows+halo_j )                         &
     &, wt_lambda_p(1-halo_i : row_length+halo_i)                       &
     &, wt_phi_p   (1-halo_i : row_length+halo_i, 1-halo_j:rows+halo_j) &
     &, wt_lambda_u(1-halo_i : row_length+halo_i)                       &
     &, wt_phi_v( 1-halo_i : row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, lambda_p_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_p_rp (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rp (1-halo_i : row_length + halo_i)                    &
     &, phi_p_rm   ( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : rows + halo_j )                         &
     &, phi_p_rp   ( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : rows + halo_j )                         &
     &, phi_v_rm   ( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : n_rows+halo_j )                         &
     &, phi_v_rp   ( 1-halo_i : row_length + halo_i,                    &
     &               1-halo_j : n_rows+halo_j )                         &
     &, recip_lambda_p_m (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_p_0 (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_p_p (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_m (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_0 (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_p (1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_phi_p_m   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : rows + halo_j )                    &
     &, recip_phi_p_0   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : rows + halo_j )                    &
     &, recip_phi_p_p   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : rows + halo_j )                    &
     &, recip_phi_p_p2  ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : rows + halo_j )                    &
     &, recip_phi_v_m   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : n_rows+halo_j )                    &
     &, recip_phi_v_0   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : n_rows+halo_j )                    &
     &, recip_phi_v_p   ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : n_rows+halo_j )                    &
     &, recip_phi_v_p2  ( 1-halo_i : row_length + halo_i,               &
     &                    1-halo_j : n_rows+halo_j )

      Real                                                              &
     &  eta_rho_levels(model_levels)                                    &
     &, eta_theta_levels(0:model_levels)

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)                                           &
     &, u_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, v_adv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &         model_levels)                                            &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         0:model_levels)

      Real                                                              &
                      ! Trig functions.
     &  cos_theta_latitude (1-off_x:row_length+off_x,                   &
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

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
                      ! Departure point co-ordinates.
     &  depart_lambda (row_length, rows_depart, model_levels)           &
     &, depart_phi (row_length, rows_depart, model_levels)              &
     &, depart_r (row_length, rows_depart, model_levels)

      Integer :: ErrorStatus    

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Local Variables.

! scalars
! arrays

! External Routines:
! subroutines
      External Ritchie

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Error trap Un-supported options.
! ----------------------------------------------------------------------

! Check grid type is recognised.

      IF (lhook) CALL dr_hook('DEPARTURE_POINT',zhook_in,zhook_handle)
      If ( type  >   6 .or. type  <   1) Then
        ErrorStatus = 10

        Call Ereport("Departure_point ", ErrorStatus,                   &
     &               "calculation given unknown grid type" )
      End If

! Check scheme choice is valid.

      If (depart_scheme <= 0 .or. depart_scheme > 1 ) Then
        ErrorStatus = 10

        Call Ereport("Departure_point ", ErrorStatus,                   &
     &               "calculation given unknown scheme choice" )
      End If

! Check scheme order choice is sensible.

      If (depart_order < 0 .or. depart_order > 20 ) Then
        ErrorStatus = 10

        Call Ereport("Departure_point ", ErrorStatus,                   &
     &               "Unrealistic value for scheme order" )
      End If

! ----------------------------------------------------------------------
! Section 2. Call appropriate Departure point scheme.
! ----------------------------------------------------------------------

      If (Depart_scheme  ==  1) Then
! use scheme based on Ritchie 1994

! DEPENDS ON: ritchie
          Call Ritchie(                                                 &
     &                  type, timestep, u_adv, v_adv, w_adv,            &
     &                  eta_rho_levels, eta_theta_levels,               &
     &                  r_rho_levels, r_theta_levels,                   &
     &                  r_at_u, r_at_v,                                 &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  glambda_p, phi_p, glambda_u, phi_v,             &
     &                  gdlambda_p, dphi_p, gdlambda_u, dphi_v,         &
     &                  grecip_dlamp, recip_dphip,                      &
     &                  grecip_dlamu, recip_dphiv,                      &
     &                  wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,   &
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
     &                  base_lambda, base_phi,                          &
     &                  lambda_p_end, phi_p_end, dlambda_p_end,         &
     &                  dphi_p_end, dphi_v_end,                         &
     &                  recip_dlam, recip_dphi, max_look,               &
     &                  look_lam, look_phi, halo_lam, halo_phi,         &
     &                  cos_theta_latitude, sec_theta_latitude,         &
     &                  sin_theta_latitude, tan_theta_latitude,         &
     &                  cos_v_latitude, sec_v_latitude,                 &
     &                  sin_v_latitude, tan_v_latitude,                 &
     &                  sin_theta_longitude,                            &
     &                  cos_theta_longitude,                            &
     &                  Depart_order, rows_depart, model_domain,        &
     &                  L_high, L_mono, high_order_scheme,              &
     &                  monotone_scheme,                                &
     &                  check_bottom_levels,                            &
     &                  interp_vertical_search_tol,                     &
     &                  L_2d_sl_geometry, first_constant_r_rho_level,   &
     &                  L_regular, LAM_max_cfl,                         &
     &                  me, n_proc, n_procx, n_procy,                   &
     &                  halo_i, halo_j, datastart,                      &
     &                  global_row_length, global_rows,                 &
     &                  g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,          &
     &                  group_2dcomm, proc_row_group, proc_col_group,   &
     &                  proc_all_group, at_extremity,                   &
     &                  L_sl_halo_reprod, off_x, off_y,                 &
     &                  depart_lambda, depart_phi, depart_r)

      End If

! End of routine.
      IF (lhook) CALL dr_hook('DEPARTURE_POINT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Departure_Point

