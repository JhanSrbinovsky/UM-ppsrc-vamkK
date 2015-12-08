! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_moist_nonhydro_conserve_mix
!
      Subroutine SL_moist_nonhydro_conserve_mix(                        &
     &                           moisture_array_size,                   &
     &                           mix_v,mix_cl,mix_cf, timestep,         &
     &                           mix_cf2,mix_rain,mix_graup,cf,cff,cfl, &
     &                           L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,   &
     &                           L_pc2,                                 &
     &                           mix_v_inter,mix_cl_inter,mix_cf_inter, &
     &                           u_adv, v_adv, w_adv,                   &
     &                           eta_rho_levels, eta_theta_levels,      &
     &                           r_rho_levels, r_theta_levels,          &
     &                           exner_theta_levels, exner_rho_levels,  &
     &                           rho_n, rho_np1,                        &
     &                           row_length, rows, n_rows,              &
     &                           model_levels, wet_model_levels,        &
     &                           delta_lambda, delta_phi,               &
     &                           base_lambda, base_phi,                 &
     &                           glambda_p, phi_p,                      &
     &                           grecip_dlamp, recip_dphip,             &
     &                           lambda_rm,lambda_rp, phi_rm, phi_rp,   &
     &                           recip_lambda_m, recip_lambda_0,        &
     &                           recip_lambda_p, recip_lambda_p2,       &
     &                           recip_phi_m, recip_phi_0,              &
     &                           recip_phi_p, recip_phi_p2,             &
     &                           recip_dlam, recip_dphi, max_look,      &
     &                           look_lam, look_phi, halo_lam, halo_phi,&
     &                           FV_cos_theta_latitude,                 &
     &                           cos_theta_latitude, sec_theta_latitude,&
     &                           sin_theta_latitude, tan_theta_latitude,&
     &                           cos_v_latitude, sec_v_latitude,        &
     &                           sin_v_latitude, tan_v_latitude,        &
     &                           sin_theta_longitude,                   &
     &                           cos_theta_longitude,                   &
     &                           sin_u_longitude, cos_u_longitude,      &
     &                           wet_to_dry_n,wet_to_dry_np1,           &
     &                         depart_lambda, depart_phi,               &
     &                         depart_r,                                &
! LAM bit
     &                           LAM_max_cfl, L_regular, rims_to_do,    &

! roar's bit
     &                           r_at_u, r_at_v,                        &
     &                           me, n_proc, n_procx, n_procy,          &
     &                           halo_i, halo_j, l_datastart,           &
     &                           g_i_pe, g_j_pe, l_2dcomm, size_2dcomm, &
     &                           group_2dcomm, at_extremity,            &
     &                           global_row_length, global_rows,        &
     &                           proc_row_group, proc_col_group,        &
     &                           proc_all_group, off_x, off_y,          &
     &                           L_sl_halo_reprod,                      &

     &                           Pi, Depart_scheme, Depart_order,       &
     &                           high_order_scheme_moist,               &
     &                           monotone_scheme_moist,                 &
     &                           L_ritchie_high, L_ritchie_mono,        &
     &                           ritchie_high_order_scheme,             &
     &                           ritchie_monotone_scheme,               &
     &                           model_domain, L_high_moist,            &
     &                           L_mono_moist, L_conserv_moist,         &
     &                           alpha_2, check_bottom_levels,          &
     &                           interp_vertical_search_tol,            &
     &                           first_constant_r_rho_level,            &
     &                           L_2d_sl_geometry,                      &
     &                           mix_v_star,mix_cl_star,mix_cf_star,    &
     &                           mix_cf2_star,mix_rain_star,            &
     &                           mix_graup_star,                        &
     &                           cf_star,cff_star,cfl_star,             &
     &                           exner_star,                            &
     &                           Error_Code)

! Purpose:
!          Performs semi-Lagrangian advection of tracers
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
     &, max_look                                                        &
                       ! max size of look-up arrays for searches
     &, l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, proc_all_group                                                  &
                       ! Group id for all processors
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
     &, moisture_array_size

      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: size_int_qc_mult  ! error check size for comms on demand

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Logical                                                           &
     &  L_sl_halo_reprod ! if true then sl code bit repoducible with
                         ! any sensible halo size

      Integer                                                           &
     &  LAM_max_cfl(2)                                                  &
     &, rims_to_do

      Integer                                                           &
     &  high_order_scheme_moist                                         &
                                 ! a code saying which high order
                           ! scheme to use for moist variables.
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
     &  L_high_moist                                                    &
                       ! True, if high order interpolation required
                       !       for moist variables.
     &, L_mono_moist                                                    &
                       ! True, if interpolation required to be monotone
                       !       for moist variables.
     &, L_conserv_moist                                                 &
                        ! True, if interpolation to be monotone and
                       !       conservative for moist variables.
     &, L_regular                                                       &
                    !  False if variable resolution
     &, L_mcr_qcf2                                                      &
                   ! True if second ice variable is in use
     &, L_mcr_qrain                                                     &
                    ! True if prognostic rain is in use
     &, L_mcr_qgraup                                                    &
                     ! True if prognostic graupel is in use
     &, L_pc2          ! True, if clouds need advecting for PC2 scheme

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
     &, Pi                                                              &
                      ! 3.1415etc
     &, base_lambda                                                     &
     &, base_phi                                                        &
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
     &  look_lam(1-halo_lam:max_look-halo_lam)                          &
     &, look_phi(1-halo_phi:max_look-halo_phi)

!  VarRes horizontal co-ordinate spacing etc.
       Real                                                             &
     &  glambda_p(1-halo_i : global_row_length + halo_i)                &
     &, phi_p( 1-halo_i : row_length + halo_i                           &
     &,        1-halo_j : rows + halo_j )                               &
     &, grecip_dlamp(1-halo_i : global_row_length + halo_i)             &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, lambda_rm (1-halo_i : row_length + halo_i)                      &
     &, lambda_rp (1-halo_i : row_length + halo_i)                      &
     &, phi_rm    ( 1-halo_i : row_length + halo_i                      &
     &,             1-halo_j : rows + halo_j )                          &
     &, phi_rp    ( 1-halo_i : row_length + halo_i                      &
     &,             1-halo_j : rows + halo_j )                          &
     &, recip_lambda_m(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_0(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_p(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_p2(1-halo_i : row_length + halo_i)                 &
     &, recip_phi_m( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_0( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_p( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_p2( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : rows + halo_j )

      Real                                                              &
     &  mix_v (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_model_levels)                                            &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &       wet_model_levels)                                          &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &       wet_model_levels)                                          &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &        wet_model_levels)                                         &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &         wet_model_levels)                                        &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &          wet_model_levels)                                       &
     &, cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &     wet_model_levels)                                            &
     &, cfl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, cff (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, mix_v_inter (1-off_x:row_length+off_x,                          &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_cl_inter (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y ,wet_model_levels)              &
     &, mix_cf_inter (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y ,wet_model_levels)              &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)

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
     &         0:model_levels)

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Real                                                              &
     &  rho_n   (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, rho_np1 (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels)

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
     &, sin_theta_longitude (row_length, rows)                          &
     &, cos_u_longitude (row_length, rows)                              &
     &, sin_u_longitude (row_length, rows)

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


! Arguments with Intent OUT. ie: Output variables.


! Star variables hold physics increments on input, latest values on
! output
      Real                                                              &
     &  mix_v_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_cl_star(1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y, wet_model_levels)               &
     &, mix_cf_star(1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y, wet_model_levels)               &
     &, mix_cf2_star(1-off_x:row_length+off_x,                          &
     &            1-off_y:rows+off_y, wet_model_levels)                 &
     &, mix_rain_star(1-off_x:row_length+off_x,                         &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, mix_graup_star(1-off_x:row_length+off_x,                        &
     &              1-off_y:rows+off_y, wet_model_levels)               &
     &, cf_star (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, cfl_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, cff_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_model_levels)                  &
     &, exner_star(1-off_x:row_length+off_x,                            &
                                                      ! Departure value
     &             1-off_y:rows+off_y, model_levels)  ! of exner.

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid



! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, k1                                                     &
                        ! Loop indices
     &, min_k, temp                                                     &
     &, check_levels                                                    &
     &, n_inputs                                                        &
     &, n_ext_fields    ! number of ext_data fields required, 1 for
                        ! global model, 2 for LAM.

      Logical                                                           &
     &  L_vector       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.

      Real                                                              &
     &  dummy

      Logical                                                           &
     &  L_do_halos                                                      &
                          ! update the halos?
     &, L_do_boundaries   ! update the boundaries?

! arrays

      Real                                                              &
     &  work(row_length, rows, model_levels)                            &
     &, work_np1(row_length, rows, model_levels)                        &
     &, drk, drkp1

      Real                                                              &
     &  super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,    &
     &           wet_model_levels,moisture_array_size)                  &
     &, data_out_super(row_length, rows,                                &
     &           wet_model_levels,moisture_array_size)
      integer count

      Integer                                                           &
     &  i_out (row_length, rows, wet_model_levels)                      &
     &, j_out (row_length, rows, wet_model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, wet_model_levels)              &
     &, weight_phi    (row_length, rows, wet_model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines:
! subroutines
      External Departure_Point, Interpolation_qcon_multi, Bi_Linear_H
      External Calc_Index

! Functions: None

! ----------------------------------------------------------------------
!  Section 0.    Set up type variable
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_MOIST_NONHYDRO_CONSERVE_MIX',zhook_in,zhook_handle)
! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, rows, wet_model_levels,           &
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

! ----------------------------------------------------------------------
! Section 1.0  Calculate q, qcl and qcf at departure point.
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
              Call bi_linear_h                                          &
     &                     (r_theta_levels(1-halo_i,1-halo_j,           &
     &                                     wet_model_levels),           &
     &                      depart_lambda(1,1,k),                       &
     &                      depart_phi(1,1,k),                          &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows,                           &
     &                      i_out(1,1,k), j_out(1,1,k),                 &
     &                      weight_lambda(1,1,k), weight_phi(1,1,k),    &
     &                      model_domain,                               &
     &                      me, n_procx, n_procy,n_proc,                &
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


!qcon block start
! ----------------------------------------------------------------------
! Section 2  Set appropriate weighted rho*r*r*delta_r at data points
! and at time t and t+deltat.
! ----------------------------------------------------------------------

! Weights come from rewriting: Sum_k(qbarr * rho*r*r*dr) =
!                              Sum_k(q * weighted average of rho*r*r*dr)

! Store in work (for current timestep n) and work_np1 (for next
! timestep n+1) -  'work' array therefore reused

! It is assumed that rho_n holds r-squared scaled current value of rho
! and that rho_np1 holds r-squared scaled value of rho at next timestep

        k = 1
! Note it is assumed that q(0) = q(1) and so q(0) contribution has
! been absorbed into that of q(1), hence different form of drk.
! This is not essential part of alogorithm and could be changed.
        Do j = 1, rows
          Do i = 1, row_length
            drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
            drk   = r_theta_levels(i,j,k)   - r_theta_levels(i,j,k-1)
          work(i,j,k)     = wet_to_dry_n(i,j,k+1)*rho_n(i,j,k+1)*drkp1  &
     &                   +  wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
        work_np1(i,j,k) = wet_to_dry_np1(i,j,k+1)*rho_np1(i,j,k+1)*drkp1&
     &                  + wet_to_dry_np1(i,j,k)  *rho_np1(i,j,k)  *drk
          End Do
        End Do

        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
              drk   = r_rho_levels(i,j,k)     - r_theta_levels(i,j,k-1)
          work(i,j,k)     = wet_to_dry_n(i,j,k+1)*rho_n(i,j,k+1)*drkp1  &
     &                    + wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
       work_np1(i,j,k) = wet_to_dry_np1(i,j,k+1)*rho_np1(i,j,k+1)*drkp1 &
     &                 + wet_to_dry_np1(i,j,k)  *rho_np1(i,j,k)  *drk
            End Do
          End Do
        End Do

        k = model_levels
        Do j = 1, rows
          Do i = 1, row_length
              drk   = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              work(i,j,k)   = wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
            work_np1(i,j,k) = wet_to_dry_np1(i,j,k)*rho_np1(i,j,k)*drk
          End Do
        End Do

!qcon block ends
! ----------------------------------------------------------------------
! Section 3  Create new variable for advection, q(n) +q_inc(phys)
! ----------------------------------------------------------------------

        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_array(i,j,k,1) = mix_v(i,j,k) + mix_v_star(i,j,k)
              super_array(i,j,k,2) = mix_cl(i,j,k) + mix_cl_star(i,j,k)
              super_array(i,j,k,3) = mix_cf(i,j,k) + mix_cf_star(i,j,k)
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
                  super_array(i,j,k,count+1) = cfl(i,j,k) +             &
                                               cfl_star (i,j,k) +       &
                                               cff(i,j,k) +             &
                                               cff_star (i,j,k) -       &
                                               cf(i,j,k) -              &
                                               cf_star (i,j,k)
                  super_array(i,j,k,count+2) = cfl(i,j,k) + cfl_star (i,j,k)
                  super_array(i,j,k,count+3) = cff(i,j,k) + cff_star (i,j,k)
                  super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
                END DO
              END DO
            END DO
          ELSE
! original method, advect the bulk cloud fraction
            DO k = 1, wet_model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  super_array(i,j,k,count+1) = cf(i,j,k) + cf_star (i,j,k)
                  super_array(i,j,k,count+2) = cfl(i,j,k) + cfl_star (i,j,k)
                  super_array(i,j,k,count+3) = cff(i,j,k) + cff_star (i,j,k)
                  super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
                END DO
              END DO
            END DO
          END IF ! l_fixbug_pc2_mixph
          count=count+4
        END IF  ! L_pc2

        If (L_mcr_qcf2) Then   ! Second cloud ice variable in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_array(i,j,k,count+1)= mix_cf2(i,j,k) +            &
     &                                      mix_cf2_star(i,j,k)
             End Do
            End Do
          End Do
          count=count+1
        End If

        If (L_mcr_qrain) Then   ! Prognostic rain variable in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_array(i,j,k,count+1)=                             &
     &                 mix_rain(i,j,k) + mix_rain_star(i,j,k)
             End Do
            End Do
          End Do
          count=count+1
        End If

        If (L_mcr_qgraup) Then  ! Prognostic graupel in use
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_array(i,j,k,count+1)=                             &
     &                 mix_graup(i,j,k) + mix_graup_star(i,j,k)
             End Do
            End Do
          End Do
        End If

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &       super_array, row_length, rows,                             &
     &       moisture_array_size*wet_model_levels,                      &
     &       halo_i, halo_j, fld_type_p, .false. )

        If (model_domain == mt_LAM .or.                                 &
     &      model_domain == mt_cyclic_LAM) Then
! overwrite X_sl with X values on LAM boundaries
          If (at_extremity(PSouth)) Then
            Do k = 1, wet_model_levels
              Do j = 1-halo_j, rims_to_do
                Do i = 1-halo_i, row_length + halo_i
                  super_array(i,j,k,1) = mix_v(i,j,k)
                  super_array(i,j,k,2) = mix_cl(i,j,k)
                  super_array(i,j,k,3) = mix_cf(i,j,k)
                End Do
              End Do
            End Do
            count=3
            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rims_to_do
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl(i,j,k) +         &
                                                   cff(i,j,k) -         &
                                                   cf(i,j,k)

                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rims_to_do
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              END IF ! l_fixbug_pc2_mixph
              count=count+4
            END IF !l_pc2
            
            if(l_mcr_qcf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rims_to_do
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qcf2
            if(l_mcr_qrain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rims_to_do
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_rain(i,j,k)
                 End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qrain
            if(l_mcr_qgraup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rims_to_do
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qgraup
          endif  ! at_extremity Psouth
          If (at_extremity(PNorth)) Then
            Do k = 1, wet_model_levels
              Do j = rows - rims_to_do + 1, rows+halo_j
                Do i = 1-halo_i, row_length + halo_i
                  super_array(i,j,k,1) = mix_v(i,j,k)
                  super_array(i,j,k,2) = mix_cl(i,j,k)
                  super_array(i,j,k,3) = mix_cf(i,j,k)
                End Do
              End Do
            End Do
            count=3
            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = rows - rims_to_do + 1, rows+halo_j
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl(i,j,k) +           &
                                                   cff(i,j,k) -           &
                                                   cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = rows - rims_to_do + 1, rows+halo_j
                    DO i = 1-halo_i, row_length + halo_i
                      super_array(i,j,k,count+1) = cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              END IF ! l_fixbug_pc2_mixph
              count=count+4
            END IF !l_pc2
            if(l_mcr_qcf2)then
              Do k = 1, wet_model_levels
                Do j = rows - rims_to_do + 1, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qcf2
            if(l_mcr_qrain)then
              Do k = 1, wet_model_levels
                Do j = rows - rims_to_do + 1, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qrain
            if(l_mcr_qgraup)then
              Do k = 1, wet_model_levels
                Do j = rows - rims_to_do + 1, rows+halo_j
                  Do i = 1-halo_i, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qgraup
          endif  ! at_extremity Pnorth
        End If  ! model_domain (mt_lam or mt_ew_cyclic lam)

        If (model_domain == mt_LAM) Then
          If (at_extremity(PWest)) Then
            Do k = 1, wet_model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = 1-halo_i, rims_to_do
                  super_array(i,j,k,1) = mix_v(i,j,k)
                  super_array(i,j,k,2) = mix_cl(i,j,k)
                  super_array(i,j,k,3) = mix_cf(i,j,k)
                End Do
              End Do
            End Do
            count=3
            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = 1-halo_i, rims_to_do
                      super_array(i,j,k,count+1) = cfl(i,j,k) +         &
                                                   cff(i,j,k) -         &
                                                   cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = 1-halo_i, rims_to_do
                      super_array(i,j,k,count+1) = cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              END IF ! l_fixbug_pc2_mixph
              count=count+4
            END IF !l_pc2
            if(l_mcr_qcf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, rims_to_do
                    super_array(i,j,k,count+1) = mix_cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qcf2
            if(l_mcr_qrain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, rims_to_do
                    super_array(i,j,k,count+1) = mix_rain(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qrain
            if(l_mcr_qgraup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = 1-halo_i, rims_to_do
                    super_array(i,j,k,count+1) = mix_graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qgraup
          endif  ! at_extremity PWest
          If (at_extremity(PEast)) Then
            Do k = 1, wet_model_levels
              Do j = 1-halo_j, rows+halo_j
                Do i = row_length - rims_to_do + 1, row_length + halo_i
                  super_array(i,j,k,1) = mix_v(i,j,k)
                  super_array(i,j,k,2) = mix_cl(i,j,k)
                  super_array(i,j,k,3) = mix_cf(i,j,k)
                End Do
              End Do
            End Do
            count=3
            IF (L_pc2) THEN
              IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = row_length - rims_to_do +1, row_length + halo_i
                      super_array(i,j,k,count+1) = cfl(i,j,k) +         &
                                                   cff(i,j,k) -         &
                                                   cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              ELSE
! original method, advect the bulk cloud fraction
                DO k = 1, wet_model_levels
                  DO j = 1-halo_j, rows+halo_j
                    DO i = row_length - rims_to_do +1, row_length + halo_i
                      super_array(i,j,k,count+1) = cf(i,j,k)
                      super_array(i,j,k,count+2) = cfl(i,j,k)
                      super_array(i,j,k,count+3) = cff(i,j,k)
                    END DO
                  END DO
                END DO
              END IF !   l_fixbug_pc2_mixph
              count=count+4
            END IF !l_pc2
            if(l_mcr_qcf2)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length - rims_to_do +1, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_cf2(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qcf2
            if(l_mcr_qrain)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length - rims_to_do +1, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_rain(i,j,k)
                  End Do
                End Do
              End Do
            count=count+1
            endif !l_mcr_qrain
            if(l_mcr_qgraup)then
              Do k = 1, wet_model_levels
                Do j = 1-halo_j, rows+halo_j
                  Do i = row_length - rims_to_do +1, row_length + halo_i
                    super_array(i,j,k,count+1) = mix_graup(i,j,k)
                  End Do
                End Do
              End Do
              count=count+1
            endif !l_mcr_qgraup
          endif  ! at_extremity Peast
        endif   ! model_domain   mt_lam

        if(l_pc2)then
! DEPENDS ON: fill_external_halos
          call fill_external_halos(                                     &
     &                         super_array(1-halo_i,1-halo_j,1,7),      &
     &                         row_length, rows,                        &
     &                         wet_model_levels,halo_i,halo_j)
        endif

! ----------------------------------------------------------------------
! Section 4  Call interpolation routine.
! ----------------------------------------------------------------------

        L_vector = .false.
        if(l_2dcomm)then
          size_int_qc_mult = row_length * rows * wet_model_levels
        else
          size_int_qc_mult = wet_model_levels * global_row_length
        endif

! DEPENDS ON: interpolation_qcon_multi
        Call Interpolation_qcon_multi(                                  &
     &                     super_array,                                 &
     &                     eta_theta_levels(1),                         &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     work, work_np1,                              &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     moisture_array_size,                         &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, wet_model_levels,          &
     &                     rows,                                        &
     &                     row_length, rows, wet_model_levels,          &
     &                     glambda_p, phi_p,                            &
     &                     lambda_rm, lambda_rp, phi_rm, phi_rp,        &
     &                     recip_lambda_m, recip_lambda_0,              &
     &                     recip_lambda_p, recip_lambda_p2,             &
     &                     recip_phi_m, recip_phi_0,                    &
     &                     recip_phi_p, recip_phi_p2,                   &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme_moist,                     &
     &                     monotone_scheme_moist,                       &
     &                     FV_cos_theta_latitude, L_regular,            &
     &                     L_vector, model_domain, L_high_moist,        &
     &                     L_mono_moist, L_conserv_moist,               &
     &                     depart_r, depart_lambda, depart_phi,         &

! Roar's bit
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j, global_row_length,           &
     &                     global_rows, row_length, rows,               &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     g_j_pe, l_2dcomm, size_2dcomm,               &
     &                     group_2dcomm, size_int_qc_mult,              &
     &                     proc_all_group,                              &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     0,0, L_sl_halo_reprod,                       &
     &                     off_x,off_y,                                 &

     &                     data_out_super,Error_Code)

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
             mix_v_star(i,j,k) = data_out_super(i,j,k,1)
             mix_cl_star(i,j,k) = data_out_super(i,j,k,2)
             mix_cf_star(i,j,k) = data_out_super(i,j,k,3)
            End Do
          End Do
        End Do
        count=3
        IF (l_pc2) THEN
          IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
            DO k = 1, model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  cf_star(i,j,k) = data_out_super(i,j,k,count+2)       &
                                  +data_out_super(i,j,k,count+3)       &
                                    -data_out_super(i,j,k,count+1)
                  cfl_star(i,j,k) = data_out_super(i,j,k,count+2)
                  cff_star(i,j,k) = data_out_super(i,j,k,count+3)
                  exner_star(i,j,k) = data_out_super(i,j,k,count+4)
                END DO
              END DO
            END DO
          ELSE
! original method, advect the bulk cloud fraction
            DO k = 1, model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  cf_star(i,j,k) = data_out_super(i,j,k,count+1)
                  cfl_star(i,j,k) = data_out_super(i,j,k,count+2)
                  cff_star(i,j,k) = data_out_super(i,j,k,count+3)
                  exner_star(i,j,k) = data_out_super(i,j,k,count+4)
                END DO
              END DO
            END DO
          END IF ! l_fixbug_pc2_mixph
          count=count+4
        END IF  ! l_pc2
        If (l_mcr_qcf2) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_cf2_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_qcf2
        If (l_mcr_qrain) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_rain_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_qrain
        If (l_mcr_qgraup) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_graup_star(i,j,k) = data_out_super(i,j,k,count+1)
              End Do
            End Do
          End Do
          count=count+1
        endif  ! l_mcr_qgraup

      End If

! End of routine.
      IF (lhook) CALL dr_hook('SL_MOIST_NONHYDRO_CONSERVE_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SL_moist_nonhydro_conserve_mix

