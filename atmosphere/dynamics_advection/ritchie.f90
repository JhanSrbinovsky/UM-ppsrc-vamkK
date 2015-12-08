! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Ritchie.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Dynamics Advection

      Subroutine Ritchie(                                               &
     &                    type, timestep, u_adv, v_adv, w_adv,          &
     &                    eta_rho_levels, eta_theta_levels,             &
     &                    r_rho_levels, r_theta_levels,                 &
     &                    r_at_u, r_at_v,                               &
     &                    row_length, rows, n_rows, model_levels,       &
     &                    delta_lambda, delta_phi,                      &
     &                    glambda_p, phi_p, glambda_u, phi_v,           &
     &                    gdlambda_p, dphi_p, gdlambda_u, dphi_v,       &
     &                    grecip_dlamp, recip_dphip,                    &
     &                    grecip_dlamu, recip_dphiv,                    &
     &                    wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v, &
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
     &                    base_lambda, base_phi,                        &
     &                    lambda_p_end, phi_p_end, dlambda_p_end,       &
     &                    dphi_p_end, dphi_v_end,                       &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    sin_theta_latitude, tan_theta_latitude,       &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    sin_v_latitude, tan_v_latitude,               &
     &                    sin_theta_longitude,                          &
     &                    cos_theta_longitude,                          &
     &                    Depart_order, rows_depart, model_domain,      &
     &                    L_high, L_mono, high_order_scheme,            &
     &                    monotone_scheme,                              &
     &                    check_bottom_levels,                          &
     &                    interp_vertical_search_tol,                   &
     &                    L_2d_sl_geometry, first_constant_r_rho_level, &
     &                    L_regular, LAM_max_cfl,                       &
     &                    me, n_proc, n_procx, n_procy,                 &
     &                    halo_i, halo_j, l_datastart,                  &
     &                    global_row_length, global_rows,               &
     &                    g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,        &
     &                    group_2dcomm, proc_row_group,                 &
     &                    proc_col_group, proc_all_group, at_extremity, &
     &                    L_sl_halo_reprod, off_x, off_y,               &
     &                    depart_lambda, depart_phi, depart_r)

! Purpose:
!          Finds departure point on the sphere for trajectory defined
!          by the given wind field using method of Ritchie.
!
! Method:
!          Is described in the proposed semi-lagrangian advection
!          scheme for the semi-implicit Unified Model integration
!          scheme. FR Division working paper No. 162.
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
USE VECTLIB_MOD, ONLY :                                                 &
      ASIN_V

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Input variables.

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
     &, l_datastart(3)                                                  &
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
    &, g_j_pe(1-halo_j:global_rows+halo_j)                       &
                                                  ! processor on my
               ! processor-column holding a given value in j direction
     &, LAM_max_cfl(2)

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_regular    ! False if variable resolution

      LOGICAL :: l_2dcomm
      INTEGER :: size_2dcomm
      INTEGER :: size_int
      INTEGER :: group_2dcomm

      Integer                                                           &
     &  type                                                            &
                    ! Defines via an integer code the nature of the
                    ! grid for which the departure points are required.
                    ! The codes are given in
                    ! terms of a primary variable that would be held
                    ! at that point and are u=1, v=2, w=3, theta=4.
     &, Depart_order                                                    &
                      ! Number of iterations to be used in finding
                      ! midpoint velocities.
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

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Logical                                                           &
     &  L_high                                                          &
                      ! True if high order interpolation scheme to be
                      ! used in Robert scheme
     &, L_mono                                                          &
                      ! True if monotone interpolation scheme to be
                      ! used in Robert scheme
     &, L_2d_sl_geometry

      Real                                                              &
     &  delta_lambda                                                    &
                       ! holds spacing between points in the i
                       ! direction.
     &, delta_phi                                                       &
                       ! holds spacing between points in the j
                       ! direction.
     &, Base_lambda                                                     &
     &, Base_phi                                                        &
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

!VarRes horizontal co-ordinate information 
      Real                                                              &
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
     &, recip_lambda_p_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_p(1-halo_i : row_length+halo_i)                  &
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


! Output variables.

      Real                                                              &
                      ! Departure point co-ordinates.
     &  depart_lambda (row_length, rows_depart, model_levels)           &
     &, depart_phi (row_length, rows_depart, model_levels)              &
     &, depart_r (row_length, rows_depart, model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k ,it                                                     &
                    ! loop counters
     &, Error_code                                                      &
     &, temp, extra_i, extra_j                                          &
     &, rlat                                                            &
             ! Use auxiliary coordinates polewards of rlat in
             ! the global model.
     &, ref_row, min_k                                                  &
     &, j0, j1, j0r, j1r, j0_SP, j1_SP, j0_NP, j1_NP                    &
                  ! derived from loop boundaries.
     &, gi, gj
      Real tmp1, tmp2

      Real                                                              &
     &  weight1                                                         &
                ! weights for vertical interpolation.
     &, weight2                                                         &
                !
     &, weight3                                                         &
                !
     &, depart_lambda_dashed                                            &
                             ! departure point on auxiliary co-ordinate
                             ! system
     &, cos_depart_phi                                                  &
                       ! trig functions for the departure point
     &, cos_depart_lambda                                               &
     &, sin_depart_phi                                                  &
     &, sin_depart_lambda                                               &
     &, exp1                                                            &
             ! used to simplify long expressions
     &, exp2                                                            &
             !
     &, G_comp                                                          &
               ! components of transformation matrix for the wind
     &, S_comp                                                          &
               ! on to the auxiliary co-ordinate system.
     &, u_dashed                                                        &
                 ! winds in auxiliary co-ordinate system.
     &, v_dashed                                                        &
                 !
     &, max_lambda                                                      &
     &, max_phi                                                         &
     &, min_lambda                                                      &
     &, min_phi                                                         &
     &, lam_bnd                                                         &
     &, phi_bnd                                                         &
     &, dummy

      Logical                                                           &
              ! used for the interpolation routine.
     &  L_vector                                                        &
     &, L_conserv

      Real                                                              &
                      !  work arrays to hold interpolated wind
                      ! fields. Names are generic as which two fields
                      ! are interpolated varies with type.
     &  interp_work1 (0:row_length+1, rows+1, model_levels)             &
                                                            ! work
     &, interp_work2 (row_length, 0:rows+1, model_levels)               &
                                                          ! arrays
                                                          ! for winds.
     &, interpw_wind1 (row_length, rows_depart, model_levels)           &
                                                              ! winds
     &, interpw_wind2 (row_length, rows_depart, model_levels)           &
                                                              ! interp-
     &, interpw_wind3 (row_length, rows_depart, model_levels)           &
                                                              ! olated.
     &, interpw_wind (0:row_length+1, 0:rows_depart+1)                  &
     &, lambda_a(row_length)                                            &
     &, phi_a(row_length, rows_depart)                                  &
     &, depart_phi_dashed(row_length, rows_depart, model_levels) !
                             ! departure point on auxiliary co-ordinate
                             ! system.
      Real                                                              &
     &  r_d_s_e (row_length, rows_depart, check_bottom_levels)

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi    (row_length, rows, model_levels)

      Real                                                              &
                      ! magnitude and direction of polar vector wind
     &  mag_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, dir_vector_sp (model_levels)

      Real :: domain_size_x 
      Real :: domain_size_y

      real :: tmp(row_length)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Functions: None

! ---------------------------------------------------------------------
! Section 0     Initialisation.
! ---------------------------------------------------------------------

      IF (lhook) CALL dr_hook('RITCHIE',zhook_in,zhook_handle)
      If ( L_regular ) Then
        lam_bnd = delta_lambda
        phi_bnd = delta_phi
        Do j = 1, rows_depart
          gj = l_datastart(2) + j - 1
          tmp2 = float(gj-1)
          Do i = 1, row_length
            phi_a(i,j) = tmp2 * delta_phi + Base_phi
          endDo
        endDo
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          tmp1 = float(gi) - 0.5
          lambda_a(i) = tmp1 * delta_lambda
        endDo
      else  ! variable resolution
        lam_bnd = dlambda_p_end
        phi_bnd = dphi_p_end
        Do j = 1, rows_depart
          Do i = 1, row_length
            phi_a(i,j) = phi_p(i,j)
          endDo
        endDo
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          lambda_a(i) = glambda_u(gi) - Base_lambda
        endDo
      endIf ! L_regular

! set latitude. Latitudes polewards of this value must take into
! account spherical geometry in the global model. Value is hard
! wired to 80 degrees.

      If (model_domain  ==  mt_Global) then
! replace delta_phi with dphi_p_end, NEED to CHECK ON THIS
        rlat = nint( Pi / ( 18. * phi_bnd ) )
        If(rlat <  1) rlat=1
      Else
        rlat=0
      End if

! Set Lateral boundary limits.
      If ( model_domain == mt_Global ) Then
        If (L_sl_halo_reprod) Then
          min_lambda = 0.
          max_lambda = 2 * pi
        Else
          min_lambda = (2-halo_i) * lam_bnd
          max_lambda = 2 * pi + (halo_i-2) * lam_bnd
        End If

        min_phi = -pi*.5
        max_phi = pi*.5

      Else If (model_domain == mt_LAM) Then
! allow for extra area defined by CFL parameter
        extra_i = (LAM_max_cfl(1) - 1)
        extra_j = (LAM_max_cfl(2) - 1)
! Check the LAM in the horizontal for each type.
        If(Type == 1)Then
! check horizontal in the LAM at u points.
          min_lambda = (0.5 - extra_i) * lam_bnd
          min_phi = base_phi - extra_j * phi_bnd
          If ( L_regular ) Then
            max_lambda = (global_row_length-1.5+extra_i) * delta_lambda
            max_phi = base_phi + (global_rows-1+extra_j) * delta_phi
          else  ! variable resolution
            max_lambda = lambda_p_end - (0.5-extra_i) * lam_bnd -       &
     &                   base_lambda
            max_phi = phi_p_end + extra_j * phi_bnd
          end If ! L_regular

        Else if(type == 2)then
! check horizontal in the LAM at v points.
          min_lambda = - extra_i * lam_bnd
          min_phi = base_phi+ (0.5-extra_j) * phi_bnd
          If ( L_regular ) Then
            max_lambda = (global_row_length-1+extra_i) * delta_lambda
            max_phi = base_phi + (global_rows-1.5+extra_j) * delta_phi
          else  ! variable resolution
            max_lambda = lambda_p_end + extra_i * lam_bnd - base_lambda
            max_phi = phi_p_end -(0.5-extra_j)* dphi_v_end
          end If ! L_regular

        Else
! check horizontal in the LAM at w points.
          min_lambda = - extra_i * lam_bnd
          min_phi = base_phi - extra_j * phi_bnd
          If ( L_regular ) Then
            max_lambda = (global_row_length-1+extra_i) * delta_lambda
            max_phi = base_phi + (global_rows-1+extra_j) * delta_phi
          else  ! variable resolution
            max_lambda = lambda_p_end + extra_i * lam_bnd - base_lambda
            max_phi = phi_p_end + extra_j * phi_bnd
          end If ! L_regular

        End If

      Else If (model_domain == mt_cyclic_LAM) Then
! check horizontal in the periodic in x LAM.

        If(Type == 1)Then
          min_phi = base_phi
          If ( L_regular ) Then
            max_phi = base_phi + delta_phi * (global_rows - 1)
          else  ! variable resolution
            max_phi = phi_p_end
          end If ! L_regular

        Else if(type == 2)then
! check horizontal in the LAM at v points.
          min_phi = base_phi + 0.5 * phi_bnd
          If ( L_regular ) Then
            max_phi = base_phi + delta_phi * (global_rows - 1.5)
          else  ! variable resolution
            max_phi = phi_p_end - 0.5 * phi_bnd
          end If ! L_regular

        Else
! check horizontal in the LAM at w points.
          min_phi = base_phi
          If ( L_regular ) Then
            max_phi = base_phi + delta_phi * (global_rows - 1)
          else  ! variable resolution
            max_phi = phi_p_end
          end If ! L_regular

        End If

        min_lambda = (2-halo_i) * lam_bnd
        If ( L_regular ) Then
          max_lambda = (global_row_length + halo_i-2) * delta_lambda
          domain_size_x=global_row_length* delta_lambda
        else  ! variable resolution
          max_lambda = lambda_p_end + (halo_i-1) * lam_bnd - base_lambda
          domain_size_x=glambda_p(global_row_length+1) - glambda_p(1)
        end If ! L_regular

      Elseif (model_domain == mt_bi_cyclic_LAM) then

        min_lambda = (2-halo_i) * lam_bnd
        min_phi = base_phi+(2-halo_j)*phi_bnd
        If ( L_regular ) Then
          max_phi = base_phi + (global_rows + halo_j-2) * delta_phi
          max_lambda = (global_row_length + halo_i-2) * delta_lambda
          domain_size_x=global_row_length* delta_lambda 
          domain_size_y=global_rows* delta_phi
        else  ! variable resolution
          max_lambda = lambda_p_end + (halo_i-1) * lam_bnd - base_lambda
          max_phi = phi_p_end + (halo_j-1) * phi_bnd
          domain_size_x=glambda_p(global_row_length+1) - glambda_p(1)
          domain_size_y=phi_p_end - base_phi  + phi_bnd
        end If ! L_regular

      End If  !  model_domain == mt_Global

! initialise depart_phi_dashed

      Do k = 1, model_levels
        Do j = 1, rows_depart
          Do i = 1, row_length
            depart_phi_dashed(i,j,k)= 0.
          End Do
        End Do
      End Do

      If (model_domain  ==  mt_Global .and. type  /=  2) Then
!  Calculate vector wind at poles. Not required for v points.

! DEPENDS ON: polar_vector_wind_n
          Call Polar_Vector_Wind_n(                                     &
     &                       v_adv,                                     &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       halo_i, halo_j, global_row_length,         &
     &                       proc_row_group, at_extremity)

      End If

      If (model_domain  ==  mt_Global) Then
! Initialise variables used to exclude the poles from the computation
! and also initialise variables which determine how many rows on
! a processor are inside the polar caps.

        j0 = 1
        j1 = rows_depart
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows_depart-1

! South Pole

        j0_SP = 1
        If (at_extremity(PSouth) .and. type  /=  2) Then
          j0_SP = j0_SP + 1
        End If

! No polar rows on this processor
        If( j0_SP + l_datastart(2) - 1  >   rlat) Then
! set j1_SP so no polar points done.
          j1_SP = j0_SP - 1
        Else
! set number of rows of south pole calculation
          j1_SP = rlat - (l_datastart(2)-1)
! check to see this does not exceed total possible on this processor.
! 36proc bug       If (j1_SP  >   l_datastart(2) + rows_depart - 1) Then
          If (j1_SP  >   rows_depart ) Then
            j1_SP = rows_depart
          End If
        End If

! North Pole.
        j1_NP = rows_depart
        If (at_extremity(PNorth) .and. type  /=  2) then
          j1_NP = j1_NP - 1
        End If

        If (type  ==  2) Then
          ref_row = global_rows - rlat
        Else
          ref_row = global_rows - rlat + 1
        End If

! No polar rows on this processor
        If( j1_NP + l_datastart(2) - 1  <   ref_row) Then
! set j0_NP so no polar points done.
          j0_NP = j1_NP + 1
        Else
! set number of rows of south pole calculation
          j0_NP =  ref_row - l_datastart(2) + 1
! check to see this does not exceed total possible on this processor.
          If (j0_NP  <   1) Then
            j0_NP = 1
          End If
        End If

! set remainder
        j0r = j1_SP + 1
        j1r = j0_NP - 1

      Else
! Limited area code
        j0 = 1
        j1 = rows_depart
        j0_SP = 1
        j1_SP = 0
        j0_NP = 1
        j1_NP = 0
        j0r = 1
        j1r = rows_depart
      End If

      If ( type == 1 ) then

! ---------------------------------------------------------------------
! Section 1.    Calculate departure point for a variable at a u point.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Section 1.1   Interpolate velocities to u points.
! ---------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& weight1,weight2,weight3,interpw_wind)                            &
!$OMP& SHARED(model_levels,rows_depart,row_length,u_adv,rows,           &
!$OMP&  r_theta_levels,r_rho_levels,w_adv,interp_work1,interpw_wind3,   &
!$OMP&  j0,j1,v_adv, wt_lambda_p,interpw_wind2,wt_phi_v,L_regular,      &
!$OMP&  interpw_wind1)
        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind1(i,j,k)=u_adv(i,j,k)
            End Do
          End Do
!       End Do

! interpolate w on to rho levels

!       Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length+1

              weight1 = r_theta_levels(i,j,k)
              weight2 = r_theta_levels(i,j,k-1)
              weight3 = r_rho_levels(i,j,k)

              interp_work1 (i,j,k) = ((weight1 - weight3)               &
     &                                 *w_adv(i,j,k-1) +                &
     &                               (weight3 - weight2)                &
     &                                 *w_adv(i,j,k) )                  &
     &                               /(weight1 - weight2)
            End Do
          End Do
!       End do

! perform horizontal interpolation to obtain w at u points.

        If ( L_regular ) Then
!       Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind3(i,j,k)= .5*(interp_work1 (i,j,k) +          &
     &             interp_work1 (i+1,j,k) )
            End Do
          End Do

! interpolate v onto the u points. They lie on the same level.

          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind2(i,j,k)= .25 * (v_adv(i,j,k)                 &
     &                              + v_adv(i+1,j,k) + v_adv(i,j-1,k)   &
     &                              + v_adv(i+1,j-1,k))
            End do
          End Do
!       End Do ! k = 1, model_levels

        else  ! variable resolution

!       Do k = 1, model_levels
! interpolation to obtain w at u points.
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind3(i,j,k)= wt_lambda_p(i+1) *                  &
     &                                      interp_work1(i,j,k) +       &
     &                                   (1.0 - wt_lambda_p(i+1)) *     &
     &                                      interp_work1(i+1,j,k)
            End Do
          End Do
! interpolate v onto the u points.
          Do j = j0-1, j1
            Do i = 1, row_length
              interpw_wind(i,j) = wt_lambda_p(i+1) * v_adv(i,j,k) +     &
     &                                   (1.0 - wt_lambda_p(i+1)) *     &
     &                                             v_adv(i+1,j,k)
            End do
          End Do
          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind2(i,j,k)= wt_phi_v(i,j) * interpw_wind(i,j-1) &
     &                                        + (1.0 - wt_phi_v(i,j)) * &
     &                                            interpw_wind(i,j)
            End do
          End Do
!       End Do  ! k = 1, model_levels

        end If ! L_regular
        End Do  ! k = 1, model_levels
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Section 1.2   Estimate velocities using depart_order iterations.
!               Do not use this method at the poles.
! ---------------------------------------------------------------------
        Do it=1,depart_order

          weight1 = timestep * timestep / 24.
          weight2 = timestep * timestep / 8.

          If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.

                depart_r(i,j,k)=r_at_u(i,j,k) - 0.5 *                   &
     &                               timestep*interpw_wind3(i,j,k)

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_at_u(i,j,k)
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2*r_at_u(i,j,k))
              End Do
            End Do
            End Do
          Else

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_u,        &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_theta_latitude,weight1,interpw_wind2,depart_phi,phi_a,      &
!$OMP&  tan_theta_latitude,weight2)
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.

                depart_r(i,j,k)=r_at_u(i,j,k) - 0.5 *                   &
     &                               timestep*interpw_wind3(i,j,k)

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_at_u(i,j,k) *                        &
     &                         (1.0 + weight1*                          &
     &                          (interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) *                 &
     &                           (sec_theta_latitude(i,j)*              &
     &                            sec_theta_latitude(i,j) - 1.)         &
     &                           - interpw_wind2(i,j,k) *               &
     &                             interpw_wind2(i,j,k) ) /             &
     &                           (r_at_u(i,j,k)*r_at_u(i,j,k)) )
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2*r_at_u(i,j,k))                       &
     &                        + tan_theta_latitude(i,j) * weight2 *     &
     &                           interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) /                 &
     &                           (r_at_u(i,j,k)*r_at_u(i,j,k))
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO 
          End If  ! model_domain = mt_LAM and 3d geometry

! For the north polar area

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,                               &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2,tmp)                 &
!$OMP& SHARED(model_levels,j0_NP,j1_NP,row_length,depart_r,r_at_u,      &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

! calculate the  midpoint departure points in the new co-ordinate
! system.

                depart_r(i,j,k)=r_at_u(i,j,k) -                         &
     &                              (timestep*interpw_wind3(i,j,k))* .5


                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                                *r_at_u(i,j,k))

                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_at_u(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)

                exp1=cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                *cos_theta_latitude(i,j) -                        &
     &                 sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                          + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

          If (at_extremity(PNorth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of r, lambda, phi at the poles use t=0 values. (will
! not be used)

            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,rows,k)= r_at_u(i,rows,k)
                depart_lambda(i,rows,k) = lambda_a(i)
                depart_phi(i,rows,k) = phi_a(i,rows-1)
              End Do
            End Do

          End If

! For the south polar area

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_SP, j1_SP,row_length,depart_r,r_at_u,     &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

! Estimate the departure points in the new co-ordinate system.

                depart_r(i,j,k)=r_at_u(i,j,k) -                         &
     &                              (timestep*interpw_wind3(i,j,k)) *.5
                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                                *r_at_u(i,j,k))
                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                   /(2*r_at_u(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)
                exp1=cos_depart_phi*sin_depart_lambda
                exp2=  cos_depart_phi*cos_depart_lambda                 &
     &                 *cos_theta_latitude(i,j) -                       &
     &                 sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                         + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

          If (at_extremity(PSouth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of r, lambda, phi at the poles use t=0 values. (will
! not be used)

            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,1,k)= r_at_u(i,1,k)
                depart_lambda(i,1,k) = lambda_a(i)
                depart_phi(i,1,k)= Base_phi
              End Do
            End Do

          End If

! ---------------------------------------------------------------------
! Check the values are in range before they are interpolated.

! DEPENDS ON: check_sl_domain
          Call check_sl_domain(                                         &
     &                         model_domain, depart_phi, depart_lambda, &
     &                         row_length, rows_depart, model_levels,   &
     &                         domain_size_x,domain_size_y,             &
     &                         max_lambda, min_lambda,                  &
     &                         max_phi, min_phi )

! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, rows, model_levels,               &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction
! DEPENDS ON: bi_linear_h
            Call Bi_Linear_H (r_rho_levels, depart_lambda, depart_phi,  &
     &                        row_length, rows, model_levels,           &
     &                        row_length, rows_depart,                  &
     &                        check_bottom_levels,                      &
     &                        row_length, rows,                         &
     &                        i_out, j_out,                             &
     &                        weight_lambda, weight_phi,                &
     &                        model_domain,                             &
     &                        me, n_procx, n_procy,n_proc,              &
     &                        halo_i, halo_j, l_datastart,              &
     &                        global_row_length, g_i_pe, at_extremity,  &
     &                        global_rows      , g_j_pe, l_2dcomm,      &
     &                        size_2dcomm, group_2dcomm,                &
     &                        0, proc_all_group,                        &
     &                        L_sl_halo_reprod, L_regular,              &
     &                        r_d_s_e )

          Do k = 1, check_bottom_levels
            Do j = 1, rows_depart
              Do i = 1, row_length

                If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                  depart_r(i,j,k) = r_d_s_e(i,j,k)
                End If

              End Do
            End Do
          End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
          temp=1
          min_k = max(1, model_levels - interp_vertical_search_tol)
          k = model_levels
          Do while (k  >=  min_k .and. temp  >   0 )
            temp = 0
            Do j = 1, rows_depart
              Do i = 1, row_length
                If (depart_r(i,j,k) >   r_at_u(i,j,model_levels) ) Then
                  depart_r(i,j,k) = r_at_u(i,j,model_levels)
                  temp = temp + 1
                End If
              End Do
            End Do
            k = k - 1
          End Do

! --------------------------------------------------------------------
! Interpolate data to the new points depart_lambda,depart_phi,depart_r
! using tri-linear interpolation.

          L_vector=.True.
          L_conserv=.False.
          if(l_2dcomm)then
            size_int=row_length* rows* model_levels
          else
            size_int=global_row_length* model_levels
          endif

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          u_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_u, dummy, r_rho_levels, 1, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, rows, model_levels,         &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_u_rm, lambda_u_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_u_m, recip_lambda_u_0,     &
     &                          recip_lambda_u_p, recip_lambda_u_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind1, dummy, dummy, Error_Code)

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          v_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_v, dummy, r_rho_levels, 2, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, n_rows, model_levels,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_v_rm, phi_v_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_v_m, recip_phi_v_0,           &
     &                          recip_phi_v_p, recip_phi_v_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_v_latitude, L_regular, L_vector,    &
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, n_rows,       &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind2, dummy, dummy, Error_Code)

          L_vector=.False.

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          w_adv, dummy, dummy,                    &
     &                          eta_theta_levels,                       &
     &                          r_theta_levels, dummy, r_theta_levels,  &
     &                          3, 1, check_bottom_levels,              &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level+1,           &
     &                          row_length, rows, model_levels+1,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind3, dummy, dummy, Error_Code)

! calculate the winds for those points on the auxiliary grid.
! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,cos_depart_phi,           &
!$OMP& G_comp,S_comp,u_dashed,v_dashed) SCHEDULE(STATIC)                &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_phi_dashed,   &
!$OMP&  depart_phi,cos_theta_latitude,sin_theta_latitude,depart_lambda,&
!$OMP&  lambda_a,interpw_wind1,interpw_wind2,j0_SP, j1_SP)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                  +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)*&
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                  /  cos_depart_phi

                S_comp=  sin_theta_latitude(i,j)*                       &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                  / cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                   - S_comp * interpw_wind2(i,j,k)

                v_dashed =  S_comp*interpw_wind1(i,j,k)                 &
     &                    +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
!         End Do

! For the south polar area.

!         Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                 +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)* &
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                 /  cos_depart_phi

                S_comp=  sin_theta_latitude(i,j)*                       &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                 /  cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                     - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                    +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
          End Do
!$OMP END PARALLEL DO

        End Do ! end loop over iterations

! ---------------------------------------------------------------------
! Section 1.3 :End of iteration step calculate final values.
!              Poles excluded.
! ---------------------------------------------------------------------

        weight1 = timestep * timestep * timestep / 8.
        weight2 = 2./3.

        If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_u,         &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=  r_at_u(i,j,k) -                         &
     &                           (timestep*interpw_wind3(i,j,k))
              depart_lambda(i,j,k) = lambda_a(i) -                      &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_at_u(i,j,k)
              depart_phi(i,j,k) = phi_a(i,j) -                          &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_at_u(i,j,k)
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_u,         &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2,               &
!$OMP&  tan_theta_latitude,weight2,weight1)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=  r_at_u(i,j,k) -                         &
     &                           (timestep*interpw_wind3(i,j,k))
              depart_lambda(i,j,k) = lambda_a(i) -                      &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_at_u(i,j,k) *                   &
     &                              (1.0 - tan_theta_latitude(i,j)      &
     &                               *interpw_wind2(i,j,k) * timestep   &
     &                               /(2. * r_at_u(i,j,k)))
              depart_phi(i,j,k) = phi_a(i,j) -                          &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_at_u(i,j,k)                       &
     &                            + (sec_theta_latitude(i,j) *          &
     &                               sec_theta_latitude(i,j) - weight2) &
     &                              *weight1 * interpw_wind2(i,j,k)     &
     &                              *interpw_wind1(i,j,k)               &
     &                              *interpw_wind1(i,j,k)               &
     &                              / (r_at_u(i,j,k) *                  &
     &                                 r_at_u(i,j,k) * r_at_u(i,j,k) )
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        End If ! model_domain  ==  mt_LAM .and. 3d_sl_geometry

! For the north polar area

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_r,r_at_u,     &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi,j0_SP, j1_SP) SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j0_NP, j1_NP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)= r_at_u(i,j,k) -                          &
     &                        (timestep*interpw_wind3(i,j,k))
              depart_lambda_dashed= -                                   &
     &                           (timestep*interpw_wind1(i,j,k))        &
     &                            /(cos(depart_phi_dashed(i,j,k))       &
     &                            *r_at_u(i,j,k))

              depart_phi_dashed(i,j,k)= -                               &
     &                                 (timestep*interpw_wind2(i,j,k))  &
     &                                 /(r_at_u(i,j,k))

! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)

              exp1=cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &             *cos_theta_latitude(i,j) -                           &
     &              sin_depart_phi*sin_theta_latitude(i,j)

              depart_lambda(i,j,k) = lambda_a(i)                        &
     &                               + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
!       End Do

! For the south polar area.

!       Do k = 1, model_levels
          Do j = j0_SP, j1_SP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_at_u(i,j,k) -                           &
     &                         (timestep*interpw_wind3(i,j,k))

              depart_lambda_dashed= -                                   &
     &                             (timestep*interpw_wind1(i,j,k))      &
     &                             /(cos(depart_phi_dashed(i,j,k))      &
     &                             *r_at_u(i,j,k))

              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                /(r_at_u(i,j,k))


! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)
              exp1= cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &              *cos_theta_latitude(i,j) -                          &
     &              sin_depart_phi*sin_theta_latitude(i,j)

              depart_lambda(i,j,k) = lambda_a(i)                        &
     &                             + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
        End Do
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Section 1.4 :  Method for the Northern and Southern boundaries.
! ---------------------------------------------------------------------
        If(model_domain  ==  mt_Global) then
! The Poles

          If (at_extremity(PSouth)) then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,1,k)= Pi + dir_vector_sp(k)
                depart_phi(i,1,k)= - Pi* .5  +                          &
     &                          ( mag_vector_sp(k)*timestep)            &
     &                           /r_at_u(i,1,k)
                depart_r(i,1,k)=r_at_u(i,1,k) -                         &
     &                        (timestep*interpw_wind3(i,1,k))
              End Do
            End Do
          End If
          If (at_extremity(PNorth)) then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,rows,k)= dir_vector_np(k)
                depart_phi(i,rows,k)= pi*.5- (mag_vector_np(k)*timestep)&
     &                              / r_at_u(i,rows,k)
                depart_r(i,rows,k)= r_at_u(i,rows,k) -                  &
     &                           (timestep*interpw_wind3(i,rows,k))

              End Do
            End Do
          End If

        End If

      Else if (type == 2) then

! ---------------------------------------------------------------------
! Section 2      Calculate departure point for a variable at a v point.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Section 2.1    Interpolate velocities to v points.
! ---------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& weight1,weight2,weight3,interpw_wind)                            &
!$OMP& SHARED(model_levels,rows_depart,row_length,u_adv,                &
!$OMP&  r_theta_levels,r_rho_levels,w_adv,interp_work1,interpw_wind3,   &
!$OMP&        v_adv, wt_lambda_u,interpw_wind2,wt_phi_p,L_regular,      &
!$OMP&  interpw_wind1)
        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind2(i,j,k)=v_adv(i,j,k)
            End Do
          End Do
!       End Do

! interpolate w on r_rho_levels

!        Do k = 1, model_levels
          Do j = 1, rows_depart+1
            Do i = 1, row_length

              weight1 = r_theta_levels(i,j,k)
              weight2 = r_theta_levels(i,j,k-1)
              weight3 = r_rho_levels(i,j,k)

              interp_work1 (i,j,k) = ((weight1 - weight3)               &
     &                                 *w_adv(i,j,k-1) +                &
     &                               (weight3 - weight2)                &
     &                                 *w_adv(i,j,k) )                  &
     &                               /(weight1 - weight2)

            End Do
          End Do
!       End do

! horizontal interpolation to v points

        If ( L_regular ) Then
!       Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind3 (i,j,k)= .5*(interp_work1 (i,j,k) +         &
     &                               interp_work1 (i,j+1,k))
            End Do
          End Do

! interpolate u onto the v points they lie on the same level
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind1(i,j,k)= .25 * (u_adv(i,j,k)                 &
     &                               + u_adv(i-1,j,k) + u_adv(i,j+1,k)  &
     &                               + u_adv(i-1,j+1,k))
            End do
          End Do
!       End Do    !  k = 1, model_levels
        else  ! variable resolution
!       Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind3(i,j,k) = wt_phi_p(i,j+1) *                  &
     &                                           interp_work1 (i,j,k) + &
     &                                        (1.0 - wt_phi_p(i,j+1)) * &
     &                                          interp_work1 (i,j+1,k)
            End Do
          End Do
! interpolate u onto the v points they lie on the same level
          Do j = 1, rows_depart + 1
            Do i = 1, row_length
              interpw_wind(i,j) = wt_lambda_u(i) * u_adv(i-1,j,k) +     &
     &                            (1.0 - wt_lambda_u(i)) * u_adv(i,j,k)
            End do
          End Do
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind1(i,j,k) = wt_phi_p(i,j+1) * interpw_wind(i,j)&
     &                                       + (1.0 - wt_phi_p(i,j+1)) *&
     &                                               interpw_wind(i,j+1)
            End do
          End Do
!       End Do    !  k = 1, model_levels
        end If ! L_regular
        End Do    !  k = 1, model_levels
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Section 2.2   Estimate velocities using depart_order iterations.
! ---------------------------------------------------------------------
! lambda_a(i) and phi_a(j) need re-defining
      If ( L_regular ) Then
        Do j = j0r, j1r
          gj = l_datastart(2) + j - 1
          tmp2 = float(gj)-0.5
          Do i = 1, row_length
            phi_a(i,j) = tmp2 * delta_phi + Base_phi
          endDo
        endDo
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          tmp1 = float(gi-1)
          lambda_a(i) = tmp1 * delta_lambda
        endDo
      else  ! variable resolution
        Do j = j0r, j1r
          Do i = 1, row_length
            phi_a(i,j) = phi_v(i,j)
          endDo
        endDo
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          lambda_a(i) = glambda_p(gi) - Base_lambda
        endDo
      endIf ! L_regular
        Do it=1,depart_order
          weight1 = timestep * timestep / 24.
          weight2 = timestep * timestep / 8.

          If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_v,        &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_v_latitude,        interpw_wind2,depart_phi,phi_a) 
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_r(i,j,k)=r_at_v(i,j,k) -                         &
     &                               (timestep*interpw_wind3(i,j,k))*.5
                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_v_latitude(i,j)* 0.5 /             &
     &                           r_at_v(i,j,k)
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2*r_at_v(i,j,k))
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO
          Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_v,        &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_v_latitude,weight1,interpw_wind2,depart_phi,phi_a,      &
!$OMP&  tan_v_latitude,weight2)

            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_r(i,j,k)=r_at_v(i,j,k) -                         &
     &                               (timestep*interpw_wind3(i,j,k))*.5
              depart_lambda(i,j,k) = lambda_a(i) -                      &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_v_latitude(i,j)* 0.5 /             &
     &                           r_at_v(i,j,k) *                        &
     &                         (1.0 + weight1*                          &
     &                          (interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) *                 &
     &                           (sec_v_latitude(i,j)*                  &
     &                            sec_v_latitude(i,j) - 1.)             &
     &                           - interpw_wind2(i,j,k) *               &
     &                             interpw_wind2(i,j,k) ) /             &
     &                           (r_at_v(i,j,k)*r_at_v(i,j,k)) )
              depart_phi(i,j,k) = phi_a(i,j) -                          &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2*r_at_v(i,j,k))                       &
     &                        + tan_v_latitude(i,j) * weight2 *         &
     &                           interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) /                 &
     &                           (r_at_v(i,j,k)*r_at_v(i,j,k))
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO
          End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_NP,j1_NP,row_length,depart_r,r_at_v,      &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_v_latitude,sin_v_latitude,lambda_a,           &
!$OMP&  depart_lambda,depart_phi,j0_SP,j1_SP) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

! calculate the  midpoint departure points in the new co-ordinate
! system.

                depart_r(i,j,k)=r_at_v(i,j,k) -                         &
     &                              (timestep*interpw_wind3(i,j,k))* .5

                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_at_v(i,j,k))

                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_at_v(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)

                exp1=cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &               *cos_v_latitude(i,j) -                             &
     &               sin_depart_phi*sin_v_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                         + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_v_latitude(i,j)+                    &
     &                        sin_depart_phi*cos_v_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
!         End Do

! For the south polar area.

!         Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

                depart_r(i,j,k)=r_at_v(i,j,k) -                         &
     &                              (timestep*interpw_wind3(i,j,k)) *.5

                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_at_v(i,j,k))

                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_at_v(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)
                exp1=cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                 *cos_v_latitude(i,j) -                           &
     &                 sin_depart_phi*sin_v_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                        + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_v_latitude(i,j)+                    &
     &                        sin_depart_phi*cos_v_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Check the values are in range before they are interpolated.

! DEPENDS ON: check_sl_domain
          Call check_sl_domain(                                         &
     &                         model_domain, depart_phi, depart_lambda, &
     &                         row_length, rows_depart, model_levels,   &
     &                         domain_size_x,domain_size_y,             &
     &                         max_lambda, min_lambda,                  &
     &                         max_phi, min_phi )

! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, n_rows, model_levels,             &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction

! DEPENDS ON: bi_linear_h
            Call Bi_Linear_H (r_rho_levels,                             &
     &                        depart_lambda, depart_phi,                &
     &                        row_length, rows, model_levels,           &
     &                        row_length, rows_depart,                  &
     &                        check_bottom_levels,                      &
     &                        row_length, rows,                         &
     &                        i_out, j_out,                             &
     &                        weight_lambda, weight_phi,                &
     &                        model_domain,                             &
     &                        me, n_procx, n_procy,n_proc,              &
     &                        halo_i, halo_j, l_datastart,              &
     &                        global_row_length, g_i_pe, at_extremity,  &
     &                        global_rows      , g_j_pe, l_2dcomm,      &
     &                        size_2dcomm, group_2dcomm,                &
     &                        2, proc_all_group,                        &
     &                        L_sl_halo_reprod, L_regular,              &
     &                        r_d_s_e)

          Do k = 1, check_bottom_levels
            Do j = 1, rows_depart
              Do i = 1, row_length

                If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                  depart_r(i,j,k) = r_d_s_e(i,j,k)
                End If

              End Do
            End Do
          End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
          temp=1
          min_k = max(1, model_levels - interp_vertical_search_tol)
          k = model_levels
          Do while (k  >=  min_k .and. temp  >   0 )
            temp = 0
            Do j = 1, rows_depart
              Do i = 1, row_length
                If (depart_r(i,j,k) >   r_at_v(i,j,model_levels) ) Then
                  depart_r(i,j,k) = r_at_v(i,j,model_levels)
                  temp = temp + 1
                End If
              End Do
            End Do
            k = k - 1
          End Do

! --------------------------------------------------------------------
! Interpolate data to the new points depart_lambda,depart_phi,r
! direction

          L_vector=.True.
          L_conserv=.False.
          if(l_2dcomm)then
            size_int=row_length* n_rows* model_levels
          else
            size_int=global_row_length* model_levels
          endif

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          u_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_u, dummy, r_rho_levels, 1, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, rows, model_levels,         &
     &                          rows,                                   &
     &                          row_length, n_rows, model_levels,       &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_u_rm, lambda_u_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_u_m, recip_lambda_u_0,     &
     &                          recip_lambda_u_p, recip_lambda_u_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 2,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind1, dummy, dummy, Error_Code)

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          v_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_v, dummy, r_rho_levels, 2, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, n_rows, model_levels,       &
     &                          rows,                                   &
     &                          row_length, n_rows, model_levels,       &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_v_rm, phi_v_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_v_m, recip_phi_v_0,           &
     &                          recip_phi_v_p, recip_phi_v_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_v_latitude, L_regular, L_vector,    &
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, n_rows,       &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 2,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind2, dummy, dummy, Error_Code)

          L_vector=.False.

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          w_adv, dummy, dummy,                    &
     &                          eta_theta_levels,                       &
     &                          r_theta_levels, dummy, r_theta_levels,  &
     &                          3, 1, check_bottom_levels,              &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level+1,           &
     &                          row_length, rows, model_levels+1,       &
     &                          rows,                                   &
     &                          row_length, n_rows, model_levels,       &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 2,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind3, dummy, dummy, Error_Code)

! calculate the winds for those points on the auxiliary grid.
! For the north polar area

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,cos_depart_phi,           &
!$OMP& G_comp,S_comp,u_dashed,v_dashed) SCHEDULE(STATIC)                &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_phi_dashed,   &
!$OMP&  depart_phi,cos_v_latitude,sin_v_latitude,depart_lambda,&
!$OMP&  lambda_a,interpw_wind1,interpw_wind2,j0_SP, j1_SP)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_v_latitude(i,j) +  &
     &                  sin(depart_phi(i,j,k))*sin_v_latitude(i,j)*     &
     &                   cos(depart_lambda(i,j,k) - lambda_a(i)))       &
     &                  /  cos_depart_phi

                S_comp= sin_v_latitude(i,j)*                            &
     &                   sin(depart_lambda(i,j,k) - lambda_a(i))        &
     &                  /  cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                   - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                   +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
!         End Do

! For the south polar area.

!         Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_v_latitude(i,j) +  &
     &                  sin(depart_phi(i,j,k))*sin_v_latitude(i,j)*     &
     &                   cos(depart_lambda(i,j,k) - lambda_a(i)))       &
     &                  /  cos_depart_phi

                S_comp= sin_v_latitude(i,j)*                            &
     &                   sin(depart_lambda(i,j,k) - lambda_a(i))        &
     &                 / cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                   - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                   +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
          End Do
!$OMP END PARALLEL DO

! End loop over depart_order iterations.
        End Do

! ---------------------------------------------------------------------
! Section 2.3 :End of iteration step calculate final departure point
!              values.
! ---------------------------------------------------------------------
        weight1 = timestep * timestep * timestep / 8.
        weight2 = 2./3.

        If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_v,         &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_v_latitude,phi_a,depart_phi,interpw_wind2)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)= r_at_v(i,j,k) -                          &
     &                         (timestep*interpw_wind3(i,j,k))
              depart_lambda(i,j,k) = lambda_a(i) -                      &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_v_latitude (i,j)               &
     &                               /r_at_v(i,j,k)
              depart_phi(i,j,k) = phi_a(i,j) -                          &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_at_v(i,j,k)
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_at_v,         &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_v_latitude,phi_a,depart_phi,interpw_wind2,                   &
!$OMP&  tan_v_latitude,weight2,weight1)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)= r_at_v(i,j,k) -                          &
     &                         (timestep*interpw_wind3(i,j,k))
              depart_lambda(i,j,k) = lambda_a(i) -                      &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_v_latitude (i,j)               &
     &                               /r_at_v(i,j,k) *                   &
     &                              (1.0 - tan_v_latitude(i,j)          &
     &                               *interpw_wind2(i,j,k) * timestep   &
     &                               /(2. * r_at_v(i,j,k)))
              depart_phi(i,j,k) = phi_a(i,j) -                          &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_at_v(i,j,k)                       &
     &                            + (sec_v_latitude(i,j) *              &
     &                               sec_v_latitude(i,j) - weight2)     &
     &                              *weight1 * interpw_wind2(i,j,k)     &
     &                              *interpw_wind1(i,j,k)               &
     &                              *interpw_wind1(i,j,k)               &
     &                              / (r_at_v(i,j,k) *                  &
     &                                 r_at_v(i,j,k) * r_at_v(i,j,k) )
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda, tmp,                          &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_r,r_at_v,     &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_v_latitude,sin_v_latitude,lambda_a,           &
!$OMP&  depart_lambda,depart_phi,j0_SP, j1_SP) SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j0_NP, j1_NP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_at_v(i,j,k) -                           &
     &                         (timestep*interpw_wind3(i,j,k))

              depart_lambda_dashed= -                                   &
     &                             (timestep*interpw_wind1(i,j,k))      &
     &                             /(cos(depart_phi_dashed(i,j,k))      &
     &                             *r_at_v(i,j,k))


              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                 /(r_at_v(i,j,k))

! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)

                 exp1=cos_depart_phi*sin_depart_lambda
                 exp2= cos_depart_phi*cos_depart_lambda                 &
     &              *cos_v_latitude(i,j) -                              &
     &               sin_depart_phi*sin_v_latitude(i,j)

              depart_lambda(i,j,k) = lambda_a(i)                        &
     &                             + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_v_latitude(i,j)+                    &
     &                        sin_depart_phi*cos_v_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
!       End Do

! For the south polar area.

!       Do k = 1, model_levels
          Do j = j0_SP, j1_SP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)= r_at_v(i,j,k) -                          &
     &                        (timestep*interpw_wind3(i,j,k))

              depart_lambda_dashed= -                                   &
     &                          (timestep*interpw_wind1(i,j,k))         &
     &                          /(cos(depart_phi_dashed(i,j,k))         &
     &                          *r_at_v(i,j,k))


              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                /(r_at_v(i,j,k))

! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)
                 exp1=cos_depart_phi*sin_depart_lambda
                 exp2= cos_depart_phi*cos_depart_lambda                 &
     &              *cos_v_latitude(i,j) -                              &
     &              sin_depart_phi*sin_v_latitude(i,j)

              depart_lambda(i,j,k) = lambda_a(i)                        &
     &                            + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_v_latitude(i,j)+                    &
     &                        sin_depart_phi*cos_v_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
        End Do
!$OMP END PARALLEL DO

      Else if(type == 3.or.type == 4) then

! ---------------------------------------------------------------------
! Section 3      Calculate departure point for a variable at a w or
!                theta point.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Section 3.1:Interpolate velocities to w/theta gridpoints.
! ---------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,rows_depart,row_length,w_adv,                &
!$OMP&  interpw_wind3)
        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              interpw_wind3(i,j,k)=w_adv(i,j,k)
            End Do
          End Do
        End Do
!$OMP END PARALLEL DO

! Interpolate u,v on to the same level as w.

        If ( L_regular ) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& weight1,weight2,weight3)                                         &
!$OMP& SHARED(model_levels,row_length,r_at_u,r_at_v,j0,j1,              &
!$OMP&  r_theta_levels,u_adv,interp_work1, v_adv, interp_work2)
        Do k = 1, model_levels-1
!  calculate u on theta levels, store in interp_work1, except at poles.

          Do j = j0, j1
            Do i = 0, row_length
              weight1 = r_at_u(i,j,k+1)
              weight2 = r_at_u(i,j,k)
              weight3 = .5 * (r_theta_levels(i,j,k) +                   &
     &                        r_theta_levels(i+1,j,k))
              interp_work1 (i,j,k) = ((weight1 - weight3)               &
     &                                 *u_adv(i,j,k) +                  &
     &                               (weight3 - weight2)                &
     &                                 *u_adv(i,j,k+1) )                &
     &                               /(weight1 - weight2)
            End Do
          End Do

!  calculate v on theta levels, store in interp_work2
          Do j = j0-1, j1
            Do i = 1, row_length
              weight1 = r_at_v(i,j,k+1)
              weight2 = r_at_v(i,j,k)
              weight3 = .5 * (r_theta_levels(i,j,k) +                   &
     &                      r_theta_levels(i,j+1,k))
              interp_work2 (i,j,k) = ((weight1 - weight3)               &
     &                                 *v_adv(i,j,k) +                  &
     &                               (weight3 - weight2)                &
     &                                 *v_adv(i,j,k+1) )                &
     &                               /(weight1 - weight2)
            End Do
          End Do
        End Do  !  k = 1, model_levels
!$OMP END PARALLEL DO
        else  ! variable resolution
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& weight1,interpw_wind)                                            &
!$OMP& SHARED(model_levels,row_length,j0,j1,wt_lambda_p,r_at_u,u_adv,   &
!$OMP& r_theta_levels,interp_work1,v_adv,interp_work2,wt_phi_p,R_at_v)  
        Do k = 1, model_levels - 1
          Do j = j0, j1
            Do i = 0, row_length
              interpw_wind(i,j) = wt_lambda_p(i+1) *                    &
     &                                     r_theta_levels(i,j,k) +      &
     &                                  (1.0 - wt_lambda_p(i+1)) *      &
     &                                      r_theta_levels(i+1,j,k)
              weight1 = ( r_at_u(i,j,k+1) - interpw_wind(i,j) ) /       &
     &                  ( r_at_u(i,j,k+1) - r_at_u(i,j,k) )
              interp_work1 (i,j,k) =  weight1 * u_adv(i,j,k) +          &
     &                                (1.0 - weight1) * u_adv(i,j,k+1)
            End Do
          End Do
!  calculate v on theta levels, store in interp_work2
          Do j = j0-1, j1
            Do i = 1, row_length
              interpw_wind(i,j) = wt_phi_p(i,j+1) *                     &
     &                                    r_theta_levels(i,j,k) +       &
     &                                  (1.0 - wt_phi_p(i,j+1)) *       &
     &                                    r_theta_levels(i,j+1,k)
              weight1 = ( r_at_v(i,j,k+1) - interpw_wind(i,j) ) /       &
     &                  ( r_at_v(i,j,k+1) - r_at_v(i,j,k))
              interp_work2 (i,j,k) =  weight1 * v_adv(i,j,k) +          &
     &                                (1.0 - weight1) * v_adv(i,j,k+1)
            End Do
          End Do

        End Do  !  k = 1, model_levels - 1
!$OMP END PARALLEL DO
        end If ! L_regular

! At the top of the model winds remain the same

        k = model_levels
! u field
        Do j = j0, j1
          Do i = 0, row_length
            interp_work1(i,j,k) = u_adv(i,j,k)
          End Do
        End Do
! v field
        Do j = j0-1, j1
          Do i = 1, row_length
            interp_work2(i,j,k) = v_adv(i,j,k)
          End Do
        End Do

! Average the interpolated u,v fields to obtain u,v at grid points
! holding w values

        If ( L_regular ) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0, j1,row_length,interpw_wind1,             &
!$OMP&   interp_work1,interpw_wind2,interp_work2)
        Do k = 1, model_levels
          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind1(i,j,k)=(interp_work1(i-1,j,k)+              &
     &                             interp_work1(i,j,k))* .5
              interpw_wind2(i,j,k)=(interp_work2(i,j-1,k)+              &
     &                              interp_work2(i,j,k))* .5
            End Do
          End Do
        End Do  !  k = 1, model_levels
!$OMP END PARALLEL DO
        else  ! variable resolution
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0, j1,row_length,interpw_wind1,             &
!$OMP&   interp_work1,interpw_wind2,interp_work2,wt_lambda_u,wt_phi_v)
        Do k = 1, model_levels
          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind1(i,j,k) = wt_lambda_u(i) *                   &
     &                                   interp_work1(i-1,j,k) +        &
     &                                  (1.0 - wt_lambda_u(i)) *        &
     &                                   interp_work1(i,j,k)
              interpw_wind2(i,j,k) = wt_phi_v(i,j) *                    &
     &                                  interp_work2(i,j-1,k) +         &
     &                                  (1.0 - wt_phi_v(i,j)) *         &
     &                                  interp_work2(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels
!$OMP END PARALLEL DO
        end If ! L_regular

! ---------------------------------------------------------------------
! Section 3.2 :  Estimate velocities using depart_order number of
!                iterations. Do not use this method for the poles.
! ---------------------------------------------------------------------

! lambda_a(i) must be re-defined
      If ( L_regular ) Then
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          tmp1 = float(gi-1)
          lambda_a(i) = tmp1 * delta_lambda
        endDo
      else  ! variable resolution
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          lambda_a(i) = glambda_p(gi) - Base_lambda
        endDo
      endIf ! L_regular
        Do it=1,depart_order
          weight1 = timestep * timestep / 24.
          weight2 = timestep * timestep / 8.
          If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_theta_levels,&
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_theta_latitude,interpw_wind2,depart_phi,phi_a)
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_theta_levels(i,j,k)
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2.0*r_theta_levels(i,j,k))
                depart_r(i,j,k)=r_theta_levels(i,j,k) -                 &
     &                              (timestep*interpw_wind3(i,j,k))*.5
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO
          Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_theta_levels,&
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_theta_latitude,weight1,interpw_wind2,depart_phi,phi_a,      &
!$OMP&  tan_theta_latitude,weight2)
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_theta_levels(i,j,k) *                &
     &                         (1.0 + weight1*                          &
     &                          (interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) *                 &
     &                           (sec_theta_latitude(i,j)*              &
     &                            sec_theta_latitude(i,j) - 1.)         &
     &                           - interpw_wind2(i,j,k) *               &
     &                             interpw_wind2(i,j,k) ) /             &
     &                           (r_theta_levels(i,j,k)                 &
     &                            *r_theta_levels(i,j,k)) )
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2.0*r_theta_levels(i,j,k))             &
     &                        + tan_theta_latitude(i,j) * weight2 *     &
     &                           interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) /                 &
     &                           (r_theta_levels(i,j,k)                 &
     &                            *r_theta_levels(i,j,k))
                depart_r(i,j,k)=r_theta_levels(i,j,k) -                 &
     &                              (timestep*interpw_wind3(i,j,k))*.5
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO
          End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,                               &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2,tmp)                 &
!$OMP& SHARED(model_levels,j0_NP,j1_NP,row_length,depart_r,             &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi,r_theta_levels) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

! calculate the  midpoint departure points in the new co-ordinate
! system.

                depart_r(i,j,k)=r_theta_levels(i,j,k) -                 &
     &                             (timestep*interpw_wind3(i,j,k))* .5


                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_theta_levels(i,j,k))


                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_theta_levels(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)

                exp1= cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                 *cos_theta_latitude(i,j) -                       &
     &                 sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                          + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))

            End Do
          End Do
!$OMP END PARALLEL DO


          If (at_extremity(PNorth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of depart_r, theta, phi at the poles,set to point 1 values.
! (will not be used to obtain final solution.)

            If ( L_regular ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,rows,k)=r_theta_levels(i,rows,k)
                depart_lambda(i,rows,k)= delta_lambda
                depart_phi(i,rows,k)= (l_datastart(2) - 1 +             &
     &               rows-1)*delta_phi + Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            else  ! variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                depart_r(i,rows,k) = r_theta_levels(i,rows,k)
                depart_lambda(i,rows,k) = gdlambda_p(gi)
                depart_phi(i,rows,k )= phi_p(i,rows-1)
              End Do
            End Do  !  k = 1, model_levels
            end If ! L_regular

          End If

! For the south polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_SP, j1_SP,row_length,depart_r,            &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi,r_theta_levels) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

                depart_r(i,j,k)=r_theta_levels(i,j,k) -                 &
     &                              (timestep*interpw_wind3(i,j,k)) *.5


                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_theta_levels(i,j,k))

                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_theta_levels(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)
                exp1= cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                *cos_theta_latitude(i,j) -                        &
     &                sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                         + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

          If (at_extremity(PSouth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of depart_r, theta, phi at the poles,set to point 1 values.
! (will not be used to obtain fianl solution.)

            If ( L_regular ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,1,k)=r_theta_levels(i,1,k)
                depart_lambda(i,1,k)= delta_lambda
                depart_phi(i,1,k)= Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            else  ! variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                depart_r(i,1,k) = r_theta_levels(i,1,k)
                depart_lambda(i,1,k) = gdlambda_p(gi)
                depart_phi(i,1,k) = Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            end If ! L_regular

          End If

! --------------------------------------------------------------------
! Check the values are in range before they are interpolated.

! DEPENDS ON: check_sl_domain
          Call check_sl_domain(                                         &
     &                         model_domain, depart_phi, depart_lambda, &
     &                         row_length, rows_depart, model_levels,   &
     &                         domain_size_x,domain_size_y,             &
     &                         max_lambda, min_lambda,                  &
     &                         max_phi, min_phi )

! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, rows, model_levels,               &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

! this if statement is in purely for the check on depart_r being inside
! the domain

          If (type  ==  3) Then
            temp = 0
          Else
            temp = 1
          End If

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction

! DEPENDS ON: bi_linear_h
          Call Bi_Linear_H (r_theta_levels(1-halo_i,1-halo_j,temp),     &
     &                      depart_lambda, depart_phi,                  &
     &                      row_length, rows, model_levels,             &
     &                      row_length, rows_depart,                    &
     &                      check_bottom_levels,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain,                               &
     &                      me, n_procx, n_procy,n_proc,                &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      0, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      r_d_s_e)

          Do k = 1, check_bottom_levels
            Do j = 1, rows_depart
              Do i = 1, row_length

                If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                  depart_r(i,j,k) = r_d_s_e(i,j,k)
                End If

              End Do
            End Do
          End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
          temp=1
          min_k = max(1, model_levels - interp_vertical_search_tol)
          k = model_levels
          Do while (k  >=  min_k .and. temp  >   0 )
            temp = 0
            Do j = 1, rows_depart
              Do i = 1, row_length
                If (depart_r(i,j,k) >                                   &
     &                   r_theta_levels(i,j,model_levels) ) Then
                  depart_r(i,j,k) = r_theta_levels(i,j,model_levels)
                  temp = temp + 1
                End If
              End Do
            End Do
            k = k - 1
          End Do

! Create departure height for u,v interpolations limited to
! r rho levels(,,model_levels) at top
! strictly should check levels below top as well.
          k = model_levels
            Do j = 1, rows_depart
              Do i = 1, row_length
                depart_r(i,j,k) = r_rho_levels(i,j,model_levels)
              End Do
            End Do

! ---------------------------------------------------------------------
! Interpolate data to the new points depart_lambda,depart_phi,depart_r
! direction

          L_vector=.True.
          L_conserv=.False.
          if(l_2dcomm)then
            size_int=row_length* rows* model_levels
          else
            size_int=global_row_length* model_levels
          endif

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          u_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_u, dummy, r_rho_levels, 1, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, rows, model_levels,         &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_u_rm, lambda_u_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_u_m, recip_lambda_u_0,     &
     &                          recip_lambda_u_p, recip_lambda_u_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind1, dummy, dummy, Error_Code)

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          v_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_v, dummy, r_rho_levels, 2, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, n_rows, model_levels,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_v_rm, phi_v_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_v_m, recip_phi_v_0,           &
     &                          recip_phi_v_p, recip_phi_v_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_v_latitude, L_regular, L_vector,    &
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, n_rows,       &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind2, dummy, dummy, Error_Code)

! Reset departure height for w interpolations as
! r theta levels(,,model_levels) at top
          k = model_levels
            Do j = 1, rows_depart
              Do i = 1, row_length
                depart_r(i,j,k) = r_theta_levels(i,j,model_levels)
              End Do
            End Do

          L_vector=.False.

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          w_adv, dummy, dummy,                    &
     &                          eta_theta_levels,                       &
     &                          r_theta_levels, dummy, r_theta_levels,  &
     &                          3, 1, check_bottom_levels,              &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level+1,           &
     &                          row_length, rows, model_levels+1,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind3, dummy, dummy, Error_Code)

! calculate the winds for those points on the auxiliary grid.
! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,cos_depart_phi,           &
!$OMP& G_comp,S_comp,u_dashed,v_dashed) SCHEDULE(STATIC)                &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_phi_dashed,   &
!$OMP&  depart_phi,cos_theta_latitude,sin_theta_latitude,depart_lambda, &
!$OMP&  lambda_a,interpw_wind1,interpw_wind2,j0_SP, j1_SP)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                  +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)*&
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                  /  cos_depart_phi

                S_comp= sin_theta_latitude(i,j)*                        &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                 / cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                     - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                     +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
!         End Do

! For the south polar area.

!         Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                  +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)*&
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                  /  cos_depart_phi

                S_comp= sin_theta_latitude(i,j)*                        &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                 / cos_depart_phi



                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                   - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                    +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
          End Do
!$OMP END PARALLEL DO

        End Do ! end loop over iterations

! ---------------------------------------------------------------------
! Section 3.3 :End of iteration step calculate final values.
!              Poles excluded.
! ---------------------------------------------------------------------
        weight1 = timestep * timestep * timestep / 8.
        weight2 = 2./3.

        If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_Theta_levels, &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=r_theta_levels(i,j,k) -                   &
     &                         (timestep*interpw_wind3(i,j,k))
                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_theta_levels(i,j,k)

                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_theta_levels(i,j,k)
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_theta_levels, &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2,               &
!$OMP&  tan_theta_latitude,weight2,weight1)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=r_theta_levels(i,j,k) -                   &
     &                         (timestep*interpw_wind3(i,j,k))
                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_theta_levels(i,j,k) *           &
     &                              (1.0 - tan_theta_latitude(i,j)      &
     &                               *interpw_wind2(i,j,k) * timestep   &
     &                               /(2. * r_theta_levels(i,j,k)))

                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_theta_levels(i,j,k)               &
     &                            + (sec_theta_latitude(i,j) *          &
     &                               sec_theta_latitude(i,j) - weight2) &
     &                              *weight1 * interpw_wind2(i,j,k)     &
     &                              *interpw_wind1(i,j,k)               &
     &                              *interpw_wind1(i,j,k)               &
     &                              / (r_theta_levels(i,j,k) *          &
     &                                 r_theta_levels(i,j,k) *          &
     &                                 r_theta_levels(i,j,k) )
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2) SCHEDULE(STATIC)    &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_r,            &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi,j0_SP, j1_SP,r_theta_levels)
        Do k = 1, model_levels
          Do j = j0_NP, j1_NP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_theta_levels(i,j,k) -                   &
     &                         (timestep*interpw_wind3(i,j,k))
              depart_lambda_dashed= -                                   &
     &                             (timestep*interpw_wind1(i,j,k))      &
     &                             /(cos(depart_phi_dashed(i,j,k))      &
     &                             *r_theta_levels(i,j,k))
              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                 /(r_theta_levels(i,j,k))

! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)

              exp1=cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &              *cos_theta_latitude(i,j) -                          &
     &              sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                            + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
!       End Do

! For the south polar area.

!       Do k = 1, model_levels
          Do j = j0_SP, j1_SP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_theta_levels(i,j,k) -                   &
     &                         (timestep*interpw_wind3(i,j,k))

              depart_lambda_dashed= -                                   &
     &                            (timestep*interpw_wind1(i,j,k))       &
     &                            /(cos(depart_phi_dashed(i,j,k))       &
     &                            *r_theta_levels(i,j,k))

              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                /(r_theta_levels(i,j,k))


! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)
              exp1=cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &              *cos_theta_latitude(i,j) -                          &
     &              sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                             + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
        End Do
!$OMP END PARALLEL DO 

! ---------------------------------------------------------------------
! Section 3.4 :  Method for the poles. Equn for r.
! ---------------------------------------------------------------------
        If (model_domain  ==  mt_Global) Then
          If (at_extremity(PSouth)) then
           Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,1,k)= Pi + dir_vector_sp(k)
                depart_phi(i,1,k)= - Pi* .5  +                          &
     &                         ( mag_vector_sp(k)*timestep)             &
     &                         /r_theta_levels(i,1,k)
                depart_r(i,1,k)= r_theta_levels(i,1,k) -                &
     &                         (timestep*w_adv(i,1,k))
              End Do
            End Do

          End If
          If (at_extremity(PNorth)) then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,rows,k)= dir_vector_np(k)
                depart_phi(i,rows,k)= pi* .5 -                          &
     &                              ( mag_vector_np(k)*timestep)        &
     &                           / r_theta_levels(i,rows,k)
                depart_r(i,rows,k)= r_theta_levels(i,rows,k) -          &
     &                           (timestep*w_adv(i,rows,k))

              End Do
            End Do
          End If
        End If

      Else if(type == 5) then

! ---------------------------------------------------------------------
! Section 4      Calculate departure point for a variable at a rho
!                point.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Section 4.1:Interpolate velocities to w/theta gridpoints.
! ---------------------------------------------------------------------

! interpolate w on to rho levels

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& weight1,weight2,weight3)                                         &
!$OMP& SHARED(model_levels,rows_depart,row_length,u_adv,rows,           &
!$OMP&  r_theta_levels,r_rho_levels,w_adv,interp_work1,interpw_wind3,   &
!$OMP&  j0,j1,v_adv, wt_lambda_u,interpw_wind2,wt_phi_v,L_regular,      &
!$OMP&  interpw_wind1)

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length

              weight1 = r_theta_levels(i,j,k)
              weight2 = r_theta_levels(i,j,k-1)
              weight3 = r_rho_levels(i,j,k)

              interpw_wind3 (i,j,k) = ((weight1 - weight3)              &
     &                                 *w_adv(i,j,k-1) +                &
     &                               (weight3 - weight2)                &
     &                                 *w_adv(i,j,k) )                  &
     &                               /(weight1 - weight2)
            End Do
          End Do
!       End do

! Interpolate u,v on to the same level as w.

! Average the  u,v fields to obtain u,v at grid points
! holding rho values

        If ( L_regular ) Then
!       Do k = 1, model_levels
          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind1(i,j,k)=(u_adv(i-1,j,k)+                     &
     &                             u_adv(i,j,k))* .5
              interpw_wind2(i,j,k)=(v_adv(i,j-1,k)+                     &
     &                              v_adv(i,j,k))* .5
            End Do
          End Do
!       End Do  !  k = 1, model_levels
        else  ! variable resolution
!       Do k = 1, model_levels
          Do j = j0, j1
            Do i = 1, row_length
              interpw_wind1(i,j,k) = wt_lambda_u(i) * u_adv(i-1,j,k) +  &
     &                                        (1.0 - wt_lambda_u(i)) *  &
     &                                                  u_adv(i,j,k)
              interpw_wind2(i,j,k)= wt_phi_v(i,j) * v_adv(i,j-1,k) +    &
     &                            (1.0 - wt_phi_v(i,j)) * v_adv(i,j,k)
            End Do
          End Do
!       End Do  !  k = 1, model_levels
        end If ! L_regular
        End Do  !  k = 1, model_levels
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Section 4.2 :  Estimate velocities using depart_order number of
!                iterations. Do not use this method for the poles.
! ---------------------------------------------------------------------

! lambda_a(i) must be re-defined
      If ( L_regular ) Then
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          tmp1 = float(gi-1)
          lambda_a(i) = tmp1 * delta_lambda
        endDo
      else  ! variable resolution
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          lambda_a(i) = glambda_p(gi) - Base_lambda
        endDo
      endIf ! L_regular
        Do it=1,depart_order
          weight1 = timestep * timestep / 24.
          weight2 = timestep * timestep / 8.
          If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_rho_levels(i,j,k)
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2.0*r_rho_levels(i,j,k))
                depart_r(i,j,k)=r_rho_levels(i,j,k) -                   &
     &                              (timestep*interpw_wind3(i,j,k))*.5
              End Do
            End Do
            End Do
          Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(model_levels,j0r, j1r,row_length,depart_r,r_rho_levels,  &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP&  sec_theta_latitude,weight1,interpw_wind2,depart_phi,phi_a,      &
!$OMP&  tan_theta_latitude,weight2)
            Do k = 1, model_levels
            Do j = j0r, j1r
              Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
! allow for spherical geometry near the poles.

                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                          timestep*interpw_wind1(i,j,k)           &
     &                          *sec_theta_latitude(i,j)* 0.5 /         &
     &                           r_rho_levels(i,j,k) *                  &
     &                         (1.0 + weight1*                          &
     &                          (interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) *                 &
     &                           (sec_theta_latitude(i,j)*              &
     &                            sec_theta_latitude(i,j) - 1.)         &
     &                           - interpw_wind2(i,j,k) *               &
     &                             interpw_wind2(i,j,k) ) /             &
     &                           (r_rho_levels(i,j,k)                   &
     &                            *r_rho_levels(i,j,k)) )
                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                         (timestep*interpw_wind2(i,j,k))          &
     &                         /(2.0*r_rho_levels(i,j,k))               &
     &                        + tan_theta_latitude(i,j) * weight2 *     &
     &                           interpw_wind1(i,j,k) *                 &
     &                           interpw_wind1(i,j,k) /                 &
     &                           (r_rho_levels(i,j,k)                   &
     &                            *r_rho_levels(i,j,k))
                depart_r(i,j,k)=r_rho_levels(i,j,k) -                   &
     &                              (timestep*interpw_wind3(i,j,k))*.5
              End Do
            End Do
            End Do
!$OMP END PARALLEL DO
          End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,                               &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2,tmp)                 &
!$OMP& SHARED(model_levels,j0_NP,j1_NP,row_length,depart_r,r_rho_levels,&
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

! calculate the  midpoint departure points in the new co-ordinate
! system.

                depart_r(i,j,k)=r_rho_levels(i,j,k) -                   &
     &                             (timestep*interpw_wind3(i,j,k))* .5


                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_rho_levels(i,j,k))


                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_rho_levels(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)

                exp1= cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                 *cos_theta_latitude(i,j) -                       &
     &                 sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                          + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

          If (at_extremity(PNorth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of depart_r, theta, phi at the poles,set to point 1 values.
! (will not be used to obtain final solution.)

            If ( L_regular ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,rows,k)=r_rho_levels(i,rows,k)
                depart_lambda(i,rows,k)= delta_lambda
                depart_phi(i,rows,k)= (l_datastart(2) - 1 +             &
     &               rows-1)*delta_phi + Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            else  ! variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                depart_r(i,rows,k) = r_rho_levels(i,rows,k)
                depart_lambda(i,rows,k) = gdlambda_p(gi)
                depart_phi(i,rows,k) = phi_p(i,rows-1) 
              End Do
            End Do  !  k = 1, model_levels
            end If ! L_regular

          End If

! For the south polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2)                     &
!$OMP& SHARED(model_levels,j0_SP,j1_SP,row_length,depart_r,r_rho_levels,&
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi) SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

                depart_r(i,j,k)=r_rho_levels(i,j,k) -                   &
     &                              (timestep*interpw_wind3(i,j,k)) *.5


                depart_lambda_dashed= -                                 &
     &                               (timestep*interpw_wind1(i,j,k))    &
     &                               /(2*cos(depart_phi_dashed(i,j,k))  &
     &                               *r_rho_levels(i,j,k))

                depart_phi_dashed(i,j,k)= -                             &
     &                                  (timestep*interpw_wind2(i,j,k)) &
     &                                  /(2*r_rho_levels(i,j,k))

! convert back to the old co-ordinate system

                cos_depart_phi= cos(depart_phi_dashed(i,j,k))
                cos_depart_lambda= cos(depart_lambda_dashed)
                sin_depart_phi= sin(depart_phi_dashed(i,j,k))
                sin_depart_lambda= sin(depart_lambda_dashed)
                exp1= cos_depart_phi*sin_depart_lambda
                exp2= cos_depart_phi*cos_depart_lambda                  &
     &                *cos_theta_latitude(i,j) -                        &
     &                sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                         + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
            End Do
          End Do
!$OMP END PARALLEL DO

          If (at_extremity(PSouth) .and.                                &
     &        model_domain  ==  mt_Global) then
! set values of depart_r, theta, phi at the poles,set to point 1 values.
! (will not be used to obtain fianl solution.)

            If ( L_regular ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_r(i,1,k)=r_rho_levels(i,1,k)
                depart_lambda(i,1,k)= delta_lambda
                depart_phi(i,1,k)= Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            else  ! variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                depart_r(i,1,k) = r_rho_levels(i,1,k)
                depart_lambda(i,1,k) = gdlambda_p(gi)
                depart_phi(i,1,k) = Base_Phi
              End Do
            End Do  !  k = 1, model_levels
            end If ! L_regular

          End If

! --------------------------------------------------------------------
! Check the values are in range before they are interpolated.

! DEPENDS ON: check_sl_domain
          Call check_sl_domain(                                         &
     &                         model_domain, depart_phi, depart_lambda, &
     &                         row_length, rows_depart, model_levels,   &
     &                         domain_size_x,domain_size_y,             &
     &                         max_lambda, min_lambda,                  &
     &                         max_phi, min_phi )

! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, rows, model_levels,               &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

! this if statement is in purely for the check on depart_r being inside
! the domain

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction

! DEPENDS ON: bi_linear_h
          Call Bi_Linear_H (r_rho_levels(1-halo_i,1-halo_j,1),          &
     &                      depart_lambda, depart_phi,                  &
     &                      row_length, rows, model_levels,             &
     &                      row_length, rows_depart,                    &
     &                      check_bottom_levels,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain,                               &
     &                      me, n_procx, n_procy,n_proc,                &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe,  l_2dcomm,       &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      0, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      r_d_s_e)

          Do k = 1, check_bottom_levels
            Do j = 1, rows_depart
              Do i = 1, row_length

                If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                  depart_r(i,j,k) = r_d_s_e(i,j,k)
                End If

              End Do
            End Do
          End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
          temp=1
          min_k = max(1, model_levels - interp_vertical_search_tol)
          k = model_levels
          Do while (k  >=  min_k .and. temp  >   0 )
            temp = 0
            Do j = 1, rows_depart
              Do i = 1, row_length
                If (depart_r(i,j,k) >                                   &
     &                   r_rho_levels(i,j,model_levels) ) Then
                  depart_r(i,j,k) = r_rho_levels(i,j,model_levels)
                  temp = temp + 1
                End If
              End Do
            End Do
            k = k - 1
          End Do

! ---------------------------------------------------------------------
! Interpolate data to the new points depart_lambda,depart_phi,depart_r
! direction

          L_vector=.True.
          L_conserv=.False.
          if(l_2dcomm)then
            size_int=row_length* rows* model_levels
          else
            size_int=global_row_length* model_levels
          endif

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          u_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_u, dummy, r_rho_levels, 1, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, rows, model_levels,         &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_u_rm, lambda_u_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_u_m, recip_lambda_u_0,     &
     &                          recip_lambda_u_p, recip_lambda_u_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind1, dummy, dummy, Error_Code)

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          v_adv, dummy, dummy,                    &
     &                          eta_rho_levels,                         &
     &                          r_at_v, dummy, r_rho_levels, 2, 1,      &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level,             &
     &                          row_length, n_rows, model_levels,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_v_rm, phi_v_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_v_m, recip_phi_v_0,           &
     &                          recip_phi_v_p, recip_phi_v_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, n_rows,       &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind2, dummy, dummy, Error_Code)


          L_vector=.False.

! DEPENDS ON: interpolation
          Call  Interpolation(                                          &
     &                          w_adv, dummy, dummy,                    &
     &                          eta_theta_levels,                       &
     &                          r_theta_levels, dummy, r_theta_levels,  &
     &                          3, 1, check_bottom_levels,              &
     &                          interp_vertical_search_tol,             &
     &                          first_constant_r_rho_level+1,           &
     &                          row_length, rows, model_levels+1,       &
     &                          rows,                                   &
     &                          row_length, rows, model_levels,         &
     &                          delta_lambda, delta_phi,                &
     &                          base_lambda, base_phi,                  &
     &                          glambda_u, phi_v,                       &
     &                          gdlambda_p, dphi_p, gdlambda_u, dphi_v, &
     &                          grecip_dlamu, recip_dphiv,              &
     &                          lambda_p_rm, lambda_p_rp,               &
     &                          phi_p_rm, phi_p_rp,                     &
     &                          recip_lambda_p_m, recip_lambda_p_0,     &
     &                          recip_lambda_p_p, recip_lambda_p_p2,    &
     &                          recip_phi_p_m, recip_phi_p_0,           &
     &                          recip_phi_p_p, recip_phi_p_p2,          &
     &                          i_out, j_out,                           &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_theta_latitude, L_regular, L_vector,&
     &                          model_domain, L_high, L_mono, L_conserv,&
     &                          depart_r, depart_lambda, depart_phi,    &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j,                         &
     &                          global_row_length, global_rows,         &
     &                          row_length, rows, n_rows, rows,         &
     &                          l_datastart, at_extremity, g_i_pe,      &
     &                          g_j_pe, l_2dcomm, size_2dcomm,          &
     &                          group_2dcomm, size_int, proc_all_group, &
     &                          proc_row_group, proc_col_group, 0,      &
     &                          0, 0, L_sl_halo_reprod, off_x, off_y,   &
     &                          interpw_wind3, dummy, dummy, Error_Code)

! calculate the winds for those points on the auxiliary grid.
! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,cos_depart_phi,           &
!$OMP& G_comp,S_comp,u_dashed,v_dashed) SCHEDULE(STATIC)                &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_phi_dashed,   &
!$OMP&  depart_phi,cos_theta_latitude,sin_theta_latitude,depart_lambda, &
!$OMP&  lambda_a,interpw_wind1,interpw_wind2,j0_SP, j1_SP)
          Do k = 1, model_levels
            Do j = j0_NP, j1_NP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                  +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)*&
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                  /  cos_depart_phi

                S_comp= sin_theta_latitude(i,j)*                        &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                 / cos_depart_phi

                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                     - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                     +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
!         End Do

! For the south polar area.

!         Do k = 1, model_levels
            Do j = j0_SP, j1_SP
              Do i = 1, row_length

                cos_depart_phi=cos(depart_phi_dashed(i,j,k))

                G_comp= ( cos(depart_phi(i,j,k))*cos_theta_latitude(i,j)&
     &                  +sin(depart_phi(i,j,k))*sin_theta_latitude(i,j)*&
     &                 cos(depart_lambda(i,j,k) - lambda_a(i)))         &
     &                  /  cos_depart_phi

                S_comp= sin_theta_latitude(i,j)*                        &
     &                 sin(depart_lambda(i,j,k) - lambda_a(i))          &
     &                 / cos_depart_phi



                u_dashed = G_comp*interpw_wind1(i,j,k)                  &
     &                   - S_comp * interpw_wind2(i,j,k)

                v_dashed = S_comp*interpw_wind1(i,j,k)                  &
     &                    +  G_comp* interpw_wind2(i,j,k)

                interpw_wind1(i,j,k)=u_dashed
                interpw_wind2(i,j,k)=v_dashed

              End Do
            End Do
          End Do
!$OMP END PARALLEL DO

        End Do ! end loop over iterations

! ---------------------------------------------------------------------
! Section 4.3 :End of iteration step calculate final values.
!              Poles excluded.
! ---------------------------------------------------------------------
        weight1 = timestep * timestep * timestep / 8.
        weight2 = 2./3.

        If (model_domain == mt_LAM .and. .not. L_2d_sl_geometry) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_rho_levels,   &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=r_rho_levels(i,j,k) -                     &
     &                         (timestep*interpw_wind3(i,j,k))
                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_rho_levels(i,j,k)

                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_rho_levels(i,j,k)
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        Else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(model_levels,j0r, j1r,row_length,depart_r,r_rho_levels,   &
!$OMP&  timestep,interpw_wind3, depart_lambda,lambda_a,interpw_wind1,   &
!$OMP& sec_theta_latitude,phi_a,depart_phi,interpw_wind2,               &
!$OMP&  tan_theta_latitude,weight2,weight1)
          Do k = 1, model_levels
          Do j = j0r, j1r
            Do i = 1, row_length

! obtain estimate of lambda,phi and r in spherical co-ordinates.
!  Not near the poles.

              depart_r(i,j,k)=r_rho_levels(i,j,k) -                     &
     &                         (timestep*interpw_wind3(i,j,k))
                depart_lambda(i,j,k) = lambda_a(i) -                    &
     &                              timestep*interpw_wind1(i,j,k)       &
     &                              *sec_theta_latitude (i,j)           &
     &                               /r_rho_levels(i,j,k) *             &
     &                              (1.0 - tan_theta_latitude(i,j)      &
     &                               *interpw_wind2(i,j,k) * timestep   &
     &                               /(2. * r_rho_levels(i,j,k)))

                depart_phi(i,j,k) = phi_a(i,j) -                        &
     &                             (timestep*interpw_wind2(i,j,k))      &
     &                             /r_rho_levels(i,j,k)                 &
     &                            + (sec_theta_latitude(i,j) *          &
     &                               sec_theta_latitude(i,j) - weight2) &
     &                              *weight1 * interpw_wind2(i,j,k)     &
     &                              *interpw_wind1(i,j,k)               &
     &                              *interpw_wind1(i,j,k)               &
     &                              / (r_rho_levels(i,j,k) *            &
     &                                 r_rho_levels(i,j,k) *            &
     &                                 r_rho_levels(i,j,k) )
            End Do
          End Do
          End Do
!$OMP END PARALLEL DO
        End If  !   mt_LAM   + 3d_sl_geometry

! For the north polar area.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,depart_lambda_dashed,     &
!$OMP&  cos_depart_phi,cos_depart_lambda,tmp,                           &
!$OMP&  sin_depart_phi,sin_depart_lambda,exp1,exp2) SCHEDULE(STATIC)    &
!$OMP& SHARED(model_levels,j0_NP, j1_NP,row_length,depart_r,            &
!$OMP&  timestep,interpw_wind3,interpw_wind1,depart_phi_dashed,         &
!$OMP&  interpw_wind2,cos_theta_latitude,sin_theta_latitude,lambda_a,   &
!$OMP&  depart_lambda,depart_phi,j0_SP, j1_SP,r_rho_levels)
        Do k = 1, model_levels
          Do j = j0_NP, j1_NP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_rho_levels(i,j,k) -                     &
     &                         (timestep*interpw_wind3(i,j,k))
              depart_lambda_dashed= -                                   &
     &                             (timestep*interpw_wind1(i,j,k))      &
     &                             /(cos(depart_phi_dashed(i,j,k))      &
     &                             *r_rho_levels(i,j,k))
              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                 /(r_rho_levels(i,j,k))

! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)

              exp1=cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &              *cos_theta_latitude(i,j) -                          &
     &              sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                            + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
!       End Do

! For the south polar area.

!       Do k = 1, model_levels
          Do j = j0_SP, j1_SP
            Do i = 1, row_length

! calculate the departure points in the new co-ordinate system

              depart_r(i,j,k)=r_rho_levels(i,j,k) -                     &
     &                         (timestep*interpw_wind3(i,j,k))

              depart_lambda_dashed= -                                   &
     &                            (timestep*interpw_wind1(i,j,k))       &
     &                            /(cos(depart_phi_dashed(i,j,k))       &
     &                            *r_rho_levels(i,j,k))

              depart_phi_dashed(i,j,k)= -                               &
     &                                (timestep*interpw_wind2(i,j,k))   &
     &                                /(r_rho_levels(i,j,k))


! convert back to the old co-ordinate system

              cos_depart_phi= cos(depart_phi_dashed(i,j,k))
              cos_depart_lambda= cos(depart_lambda_dashed)
              sin_depart_phi= sin(depart_phi_dashed(i,j,k))
              sin_depart_lambda= sin(depart_lambda_dashed)
              exp1=cos_depart_phi*sin_depart_lambda
              exp2= cos_depart_phi*cos_depart_lambda                    &
     &              *cos_theta_latitude(i,j) -                          &
     &              sin_depart_phi*sin_theta_latitude(i,j)

                depart_lambda(i,j,k) = lambda_a(i)                      &
     &                             + Atan2(exp1,exp2)

                tmp(i)=  cos_depart_phi * cos_depart_lambda             &
     &                        * sin_theta_latitude(i,j)+                &
     &                        sin_depart_phi*cos_theta_latitude(i,j)
              End Do
              call asin_v(row_length, tmp, depart_phi(1,j,k))
          End Do
        End Do
!$OMP END PARALLEL DO

! ---------------------------------------------------------------------
! Section 4.4 :  Method for the poles. Equn for r.
! ---------------------------------------------------------------------
        If (model_domain  ==  mt_Global) Then
          If (at_extremity(PSouth)) then
           Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,1,k)= Pi + dir_vector_sp(k)
                depart_phi(i,1,k)= - Pi* .5  +                          &
     &                         ( mag_vector_sp(k)*timestep)             &
     &                         /r_rho_levels(i,1,k)
                depart_r(i,1,k)= r_rho_levels(i,1,k) -                  &
     &                         (timestep*w_adv(i,1,k))
              End Do
            End Do

          End If
          If (at_extremity(PNorth)) then
            Do k = 1, model_levels
              Do i = 1, row_length
                depart_lambda(i,rows,k)= dir_vector_np(k)
                depart_phi(i,rows,k)= pi* .5 -                          &
     &                              ( mag_vector_np(k)*timestep)        &
     &                           / r_rho_levels(i,rows,k)
                depart_r(i,rows,k)= r_rho_levels(i,rows,k) -            &
     &                           (timestep*w_adv(i,rows,k))

              End Do
            End Do
          End If
        End If
! End conditional over type to perform the departure point calculation
! for.
      End if

! ---------------------------------------------------------------------
! Section 5. Check the final values for the departure point lie inside
!            the data domain.
! ---------------------------------------------------------------------

! Check the values are in range before they are interpolated.

! DEPENDS ON: check_sl_domain
      Call check_sl_domain(                                             &
     &                     model_domain, depart_phi, depart_lambda,     &
     &                     row_length, rows_depart, model_levels,       &
     &                     domain_size_x,domain_size_y,                 &
     &                     max_lambda, min_lambda,                      &
     &                     max_phi, min_phi )

! DEPENDS ON: calc_index
      Call Calc_Index(                                                  &
     &                    row_length, rows_depart, check_bottom_levels, &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

      If (type  ==  3 .or. type  ==  4) Then
! this if statement is in purely for the check on depart_r being inside
! the domain

        If (type  ==  3) Then
          temp = 0
        Else
          temp = 1
        End If

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction

! DEPENDS ON: bi_linear_h
          Call Bi_Linear_H (r_theta_levels(1-halo_i,1-halo_j,temp),     &
     &                      depart_lambda, depart_phi,                  &
     &                      row_length, rows, model_levels,             &
     &                      row_length, rows_depart,                    &
     &                      check_bottom_levels,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain,                               &
     &                      me, n_procx, n_procy,n_proc,                &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      1, proc_all_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      r_d_s_e)

        Do k = 1, check_bottom_levels
          Do j = 1, rows_depart
            Do i = 1, row_length

              If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                depart_r(i,j,k) = r_d_s_e(i,j,k)
              End If

            End Do
          End Do
        End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
        temp=1
        min_k = max(1, model_levels - interp_vertical_search_tol)
        k = model_levels
        Do while (k  >=  min_k .and. temp  >   0 )
          temp = 0
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_r(i,j,k) >                                     &
     &                 r_theta_levels(i,j,model_levels) ) Then
                depart_r(i,j,k) =                                       &
     &                        r_theta_levels(i,j,model_levels)
                temp = temp + 1
              End If
            End Do
          End Do
          k = k - 1
        End Do

      Else

! check not below bottom data boundary.
! Limit search to bottom few levels depending on user set parameter

! Calculate r at the bottom level at the departure point
! also in just the lambda direction and in just the phi direction

! set pole handling flag, don't check poles for u points,
! check all points for v and rho
          If (type  ==  1) Then
            temp = 0
          Else If (type  ==  5) Then
            temp = 1
          Else
            temp = 2
          End If

! DEPENDS ON: bi_linear_h
          Call Bi_Linear_H (r_rho_levels,                               &
     &                      depart_lambda, depart_phi,                  &
     &                      row_length, rows, model_levels,             &
     &                      row_length, rows_depart,                    &
     &                      check_bottom_levels,                        &
     &                      row_length, rows,                           &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi,                  &
     &                      model_domain,                               &
     &                      me, n_procx, n_procy,n_proc,                &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      global_rows      , g_j_pe, l_2dcomm,        &
     &                      size_2dcomm, group_2dcomm,                  &
     &                      temp, proc_all_group,                       &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      r_d_s_e)

        Do k = 1, check_bottom_levels
          Do j = 1, rows_depart
            Do i = 1, row_length

              If (depart_r(i,j,k) <   r_d_s_e(i,j,k) ) Then
! move trajectory up to lowest level of data.
                depart_r(i,j,k) = r_d_s_e(i,j,k)
              End If

            End Do
          End Do
        End Do

! check not above top data boundary.
! Limit search to top few levels depending on user set parameter
! assume once one level has no data above top then no lower level
! has either.
        temp=1
        min_k = max(1, model_levels - interp_vertical_search_tol)
        k = model_levels
        Do while (k  >=  min_k .and. temp  >   0 )
          temp = 0
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_r(i,j,k) >   r_rho_levels(i,j,model_levels) )  &
     &          Then
                depart_r(i,j,k) = r_rho_levels(i,j,model_levels)
                temp = temp + 1
              End If
            End Do
          End Do
          k = k - 1
        End Do

      End If

! ---------------------------------------------------------------------
! End of routine.
      IF (lhook) CALL dr_hook('RITCHIE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Ritchie

