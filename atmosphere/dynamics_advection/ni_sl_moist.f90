! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_SL_MOIST
!
      Subroutine NI_SL_MOIST(                                           &
     &                       moisture_array_size,                       &
     &                       q, qcl, qcf, cf, cff, cfl,                 &
     &                       qcf2, qrain, qgraup,                       &
     &                       q_phys1, qcl_phys1, qcf_phys1,             &
     &                       cf_phys1, cff_phys1,cfl_phys1,             &
     &                       qcf2_phys1, qrain_phys1, qgraup_phys1,     &
     &                       mix_v, mix_cl, mix_cf,                     &
     &                       mix_cf2, mix_rain, mix_graup,              &
     &                       mix_v_phys1, mix_cl_phys1, mix_cf_phys1,   &
     &                       mix_cf2_phys1, mix_rain_phys1,             &
     &                       mix_graup_phys1,                           &
     &                       L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,     &
     &                       L_pc2, L_mix_ratio,                        &
     &                           eta_theta_levels,                      &
     &                           r_rho_levels, r_theta_levels,          &
     &                           exner_theta_levels,                    &
     &                           rho_n, rho_np1,                        &
     &                           depart_lambda, depart_phi,             &
     &                           depart_r_theta,                        &
     &                           row_length, rows, n_rows,              &
     &                           model_levels, wet_levels,              &
     &                           delta_lambda, delta_phi,               &
     &                           base_lambda, base_phi,                 &
     &                           glambda_p, phi_p,                      &
     &                           grecip_dlamp, recip_dphip,             &
     &                           lambda_rm, lambda_rp, phi_rm, phi_rp,  &
     &                           recip_lambda_m, recip_lambda_0,        &
     &                           recip_lambda_p, recip_lambda_p2,       &
     &                           recip_phi_m, recip_phi_0,              &
     &                           recip_phi_p, recip_phi_p2,             &
     &                           recip_dlam, recip_dphi, max_look,      &
     &                           look_lam, look_phi, halo_lam, halo_phi,&
     &                           FV_cos_theta_latitude,                 &
     &                           wet_to_dry_n, wet_to_dry_np1,          &
     &                           L_regular, rims_to_do,                 &
     &                           mype, nproc, nproc_x, nproc_y,         &
     &                           halo_i, halo_j, datastart,             &
     &                           g_i_pe, g_j_pe, l_2dcomm, size_2dcomm, &
     &                           group_2dcomm, at_extremity,            &
     &                           global_row_length, global_rows,        &
     &                           gc_proc_row_group, gc_proc_col_group,  &
     &                           gc_proc_all_group, offx, offy,         &
     &                           L_sl_halo_reprod,                      &
     &                           high_order_scheme_moist,               &
     &                           monotone_scheme_moist,                 &
     &                           model_domain,  L_high_moist,           &
     &                           L_mono_moist, L_conserv_moist,         &
     &                           check_bottom_levels,                   &
     &                           interp_vertical_search_tol,            &
     &                           first_constant_r_rho_level,            &
     &                           Errorstatus )


! Purpose: Top level routine to SL_Moist
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
     &, wet_levels                                                      &
                   ! Number of model levels where moisture held
     &, mype                                                            &
                       ! My processor number
     &, nproc                                                           &
                    ! Total number of processors
     &, nproc_x                                                         &
                     ! Number of processors in longitude
     &, nproc_y                                                         &
                     ! Number of processors in latitude
     &, halo_i                                                          &
                     ! Size of large halo in i.
     &, halo_j                                                          &
                     ! Size of large halo in j.
     &, offx                                                            &
                    ! Size of small halo in i
     &, offy                                                            &
                    ! Size of small halo in j.
     &, max_look                                                        &
                       ! max size of look-up arrays for searches
     &, datastart(3)                                                    &
                     ! First gridpoints held by this processor.
     &, rims_to_do                                                      &
     &, gc_proc_row_group                                               &
                          ! Group id for processors on the same row
     &, gc_proc_col_group                                               &
                          ! Group id for processors on the same column
     &, gc_proc_all_group                                               &
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

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Logical                                                           &
     &  L_sl_halo_reprod ! if true then sl code bit repoducible with
                         ! any sensible halo size
      Integer                                                           &
     &  high_order_scheme_moist                                         &
                                 ! a code saying which high order
                           ! scheme to use for moist variables.
     &, monotone_scheme_moist                                           &
                              ! a code saying which monotone
                           ! scheme to use for moist variables.
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
     &, L_mix_ratio                                                     &
     &, L_regular                                                       &
     &, L_mcr_qcf2                                                      &
     &, L_mcr_qrain                                                     &
     &, L_mcr_qgraup                                                    &
     &, L_pc2

      Real                                                              &
     &  delta_lambda                                                    &
                      ! holds spacing between points in the i
                      ! direction for the input data field.
     &, delta_phi                                                       &
                      ! holds spacing between points in the j
                      ! direction for the input data field.
     &, base_lambda                                                     &
     &, base_phi                                                        &
     &, recip_dlam                                                      &
     &, recip_dphi

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


      Real, Intent (InOut) ::                                           &
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     wet_levels)                                                  &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_levels)                                                &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_levels)                                                &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &       wet_levels)                                                &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &       wet_levels)                                                &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &       wet_levels)                                                &
     &, cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &       wet_levels)                                                &
     &, cfl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_levels)                                                &
     &, cff (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_levels)

      Real                                                              &
     & eta_theta_levels(0:model_levels)

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
     &  rho_n (1-offx:row_length+offx,                                  &
     &                  1-offy:rows+offy, model_levels)                 &
     &, rho_np1 (1-offx:row_length+offx,                                &
     &                  1-offy:rows+offy, model_levels)                 &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)              &
     &, wet_to_dry_n (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &          model_levels)

      Real                                                              &
                      ! Trig functions.
     &  FV_cos_theta_latitude (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy)


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
     &, exner_star(1-offx:row_length+offx,                              &
                                                    ! Departure value
     &             1-offy:rows+offy, model_levels)                      &
                                                    ! of exner.
     &, depart_r_theta (row_length, rows, model_levels)   ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.

! Star variables hold physics increments on input, latest values on
! output
      Real, Intent (InOut) ::                                           &
     &  q_phys1(1-offx:row_length+offx,                                 &
     &          1-offy:rows+offy, wet_levels)                           &
     &, qcl_phys1(1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, wet_levels)                         &
     &, qcf_phys1(1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, wet_levels)                         &
     &, qcf2_phys1(1-offx:row_length+offx,                              &
     &          1-offy:rows+offy, wet_levels)                           &
     &, qrain_phys1(1-offx:row_length+offx,                             &
     &            1-offy:rows+offy, wet_levels)                         &
     &, qgraup_phys1(1-offx:row_length+offx,                            &
     &            1-offy:rows+offy, wet_levels)                         &
     &, cf_phys1(1-offx:row_length+offx,                                &
     &            1-offy:rows+offy, wet_levels)                         &
     &, cff_phys1(1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, wet_levels)                         &
     &, cfl_phys1(1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, wet_levels)

      Real, Intent (InOut) ::                                           &
     &  mix_v (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cf(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels)                        &
     &, mix_cf2(1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, wet_levels)                       &
     &, mix_rain(1-halo_i:row_length+halo_i,                            &
     &           1-halo_j:rows+halo_j, wet_levels)                      &
     &, mix_graup(1-halo_i:row_length+halo_i,                           &
     &            1-halo_j:rows+halo_j, wet_levels)                     &
     &, mix_v_phys1  (1-offx:row_length+offx,                           &
     &               1-offy:rows+offy,wet_levels)                       &
     &, mix_cl_phys1 (1-offx:row_length+offx,                           &
     &               1-offy:rows+offy,wet_levels)                       &
     &, mix_cf_phys1 (1-offx:row_length+offx,                           &
     &               1-offy:rows+offy,wet_levels)                       &
     &, mix_cf2_phys1  (1-offx:row_length+offx,                         &
     &               1-offy:rows+offy,wet_levels)                       &
     &, mix_rain_phys1 (1-offx:row_length+offx,                         &
     &               1-offy:rows+offy,wet_levels)                       &
     &, mix_graup_phys1 (1-offx:row_length+offx,                        &
     &               1-offy:rows+offy,wet_levels)

      Integer                                                           &
     &  Errorstatus    ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('NI_SL_MOIST',zhook_in,zhook_handle)
      IF(l_mix_ratio)then

! DEPENDS ON: sl_moist_nonhydro_conserve
        Call SL_moist_nonhydro_conserve(                                &
     &                 moisture_array_size,                             &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 cf, cff, cfl,                                    &
     &                 mix_v_phys1, mix_cl_phys1, mix_cf_phys1,         &
     &                 mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1,  &
     &                 cf_phys1, cff_phys1, cfl_phys1,                  &
     &                 exner_star,                                      &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,    &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels,  r_theta_levels,                   &
     &                 exner_theta_levels,                              &
     &                 rho_n, rho_np1,                                  &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 delta_lambda, delta_phi,                         &
     &                 base_lambda, base_phi,                           &
     &                 glambda_p, phi_p,                                &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 lambda_rm, lambda_rp, phi_rm, phi_rp,            &
     &                 recip_lambda_m, recip_lambda_0,                  &
     &                 recip_lambda_p, recip_lambda_p2,                 &
     &                 recip_phi_m, recip_phi_0,                        &
     &                 recip_phi_p, recip_phi_p2,                       &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi,                              &
     &                 halo_lam, halo_phi,                              &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 L_regular, rims_to_do,                           &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,           &
     &                 group_2dcomm, at_extremity,                      &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group, gc_proc_col_group,            &
     &                 gc_proc_all_group, offx, offy,                   &
     &                 L_sl_halo_reprod,                                &
     &                 high_order_scheme_moist,                         &
     &                 monotone_scheme_moist,                           &
     &                 model_domain,   L_high_moist,                    &
     &                 L_mono_moist, L_conserv_moist,                   &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 ErrorStatus)

      else  ! specific humidity

! DEPENDS ON: sl_moist_nonhydro_conserve
        Call SL_moist_nonhydro_conserve(                                &
     &                 moisture_array_size,                             &
     &                 q, qcl, qcf,                                     &
     &                 qcf2, qrain, qgraup,                             &
     &                 cf, cff, cfl,                                    &
     &                 q_phys1, qcl_phys1, qcf_phys1,                   &
     &                 qcf2_phys1, qrain_phys1, qgraup_phys1,           &
     &                 cf_phys1, cff_phys1, cfl_phys1,                  &
     &                 exner_star,                                      &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,    &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels, r_theta_levels,                    &
     &                 exner_theta_levels,                              &
     &                 rho_n, rho_np1,                                  &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 delta_lambda, delta_phi,                         &
     &                 base_lambda, base_phi,                           &
     &                 glambda_p, phi_p,                                &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 lambda_rm, lambda_rp, phi_rm, phi_rp,            &
     &                 recip_lambda_m, recip_lambda_0,                  &
     &                 recip_lambda_p, recip_lambda_p2,                 &
     &                 recip_phi_m, recip_phi_0,                        &
     &                 recip_phi_p, recip_phi_p2,                       &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi,                              &
     &                 halo_lam, halo_phi,                              &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 L_regular, rims_to_do,                           &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,           &
     &                 group_2dcomm, at_extremity,                      &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group, gc_proc_col_group,            &
     &                 gc_proc_all_group, offx, offy,                   &
     &                 L_sl_halo_reprod,                                &
     &                 high_order_scheme_moist,                         &
     &                 monotone_scheme_moist,                           &
     &                 model_domain,   L_high_moist,                    &
     &                 L_mono_moist, L_conserv_moist,                   &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 ErrorStatus)

          endif     !  L_mix_ratio

      IF (lhook) CALL dr_hook('NI_SL_MOIST',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_SL_MOIST

