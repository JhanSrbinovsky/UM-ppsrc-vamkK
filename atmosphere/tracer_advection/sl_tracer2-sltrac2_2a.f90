! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_tracer2
!

      Subroutine SL_tracer2(                                            &
     &                           super_array_size,                      &
     &                           super_array, super_tracer_phys2,       &
     &                           eta_theta_levels,                      &
     &                           r_rho_levels, r_theta_levels,          &
     &                           rho_n, rho_np1,                        &
     &                           row_length, rows, model_levels,        &
     &                           bl_levels, delta_lambda, delta_phi,    &
     &                           glambda_p, phi_p,                      &
     &                           grecip_dlamp, recip_dphip,             &
     &                           lambda_p_rm, lambda_p_rp,              &
     &                           phi_p_rm, phi_p_rp,                    &
     &                           recip_lambda_p_m, recip_lambda_p_0,    &
     &                           recip_lambda_p_p, recip_lambda_p_p2,   &
     &                           recip_phi_p_m, recip_phi_p_0,          &
     &                           recip_phi_p_p, recip_phi_p_p2,         &
     &                           base_lambda, base_phi,                 &
     &                           recip_dlam, recip_dphi, max_look,      &
     &                           look_lam, look_phi,                    &
     &                           halo_lam, halo_phi,                    &
     &                           FV_cos_theta_latitude,                 &
     &                           wet_to_dry_n, wet_to_dry_np1,          &
! LAM bit
     &                           me, n_proc, n_procx, n_procy,          &
     &                           halo_i, halo_j, l_datastart,           &
     &                           g_i_pe, g_j_pe, l_2dcomm, size_2dcomm, &
     &                           group_2dcomm, at_extremity,            &
     &                           global_row_length,                     &
     &                           global_rows,                           &
     &                           gc_all_proc_group,                     &
     &                           proc_row_group,                        &
     &                           proc_col_group, off_x, off_y,          &
     &                           L_regular, L_sl_halo_reprod,           &
     &                           high_order_scheme_moist,               &
     &                           monotone_scheme_moist,                 &
     &                           model_domain, L_high_moist,            &
     &                           L_mono_moist, L_conserv_moist,         &
     &                           check_bottom_levels,                   &
     &                           interp_vertical_search_tol,            &
     &                           first_constant_r_rho_level,            &
     &                           depart_lambda, depart_phi,             &
     &                           depart_r,                              &
     &                           CO2, L_CO2_interactive,                &
     &                           murk, L_murk_advect,                   &
     &                           soot_new, soot_agd, soot_cld, L_soot,  &
     &                           bmass_new, bmass_agd, bmass_cld,       &
     &                           L_biomass,                             &
     &                           ocff_new, ocff_agd, ocff_cld, l_ocff,  &
     &                           DUST_DIV1,DUST_DIV2,DUST_DIV3,         &
     &                           DUST_DIV4,DUST_DIV5,DUST_DIV6,         &
     &                           L_DUST,                                &
     &                           so2, so4_aitken, so4_accu,             &
     &                           so4_diss, nh3, dms,                    &
     &                           L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms, &
     &                           nitr_acc, nitr_diss, L_nitrate,        &
     &                           tracers, tr_levels, tr_vars,           &
     &                           tracer_phys1, tracer_phys2,            &
     &                           tr_lbc_vars,tracer_lbcs,               &
     &                           A_max_trvars,A_tr_active_lbc_index,    &
     &                           tracers_ukca,tr_ukca,                  &
     &                           ukca_tracer_phys1, ukca_tracer_phys2,  &
     &                           tr_lbc_ukca,ukca_tracer_lbcs,          &
     &                           A_max_ukcavars,                        &
     &                           UKCA_tr_active_lbc_index,              &
     &                           rimwidth,rimweights,                   &
     &                           lenrim,lbc_size,lbc_start,             &
     &                           i_start, i_end, j_start, j_end,        &
     &                           L_USE_CARIOLLE, OZONE_TRACER,          &
     &                           l_qpos_diag_pr, qpos_diag_limit,       &
     &                           q_pos_tracer_method,                   &
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
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE dust_parameters_mod, ONLY: l_twobin_dust
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
     &, model_levels                                                    &
                     ! Number of model levels.
     &, bl_levels                                                       &
                     ! Number of boundary layer levels
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
     &, l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, max_look                                                        &
                             ! max size of look-up arrays for searches
     &, gc_all_proc_group                                               &
     &, global_row_length                                               &
                            ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of rows
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                             ! processor on my processor-row
                             ! holding a given value in i direction
     &, g_j_pe(1-halo_j:global_rows      +halo_j)                       &
                             ! processor on my processor-column
                             ! holding a given value in j direction
     &, rimwidth                                                        &
                       ! Width of boundaries in LBCs     
     &, lenrim                                                          &
                       ! Size of single level of LBC data
     &, lbc_size(4)                                                     &
                       ! Size of each side of LBC data   
     &, lbc_start(4)                                                    
                       ! Start of each side in LBC data

      REAL rimweights(rimwidth)
                       ! Weights to apply to the LBCs    

      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: size_int_qc_mult ! error check size for comms on demand

      Integer                                                           &
     &  tr_levels                                                       &
     &, tr_vars                                                         &
     &, tr_ukca                                                         &
     &, tr_lbc_vars                                                     &
                       ! No. of active tracer lbcs
     &, tr_lbc_ukca                                                     &
                       ! No. of active UKCA tracer lbcs
     &, A_max_trvars                                                    &
                       ! Max number of tracer vars
     &, A_max_ukcavars                                                  &
                       ! Max number of UKCA tracer vars
     &, A_tr_active_lbc_index(A_max_trvars)                             &
                       ! Active tracer lbcs index
     &, UKCA_tr_active_lbc_index(A_max_ukcavars)                       
                       ! Active UKCA tracer lbcs index

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

      Integer   i_start, i_end, j_start, j_end

      Integer, Intent(In) :: q_pos_tracer_method ! which qpos method to use

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
     &, L_regular

      Logical, Intent(In) :: l_qpos_diag_pr  ! True for qpos diagnostic prints

      Real                                                              &
     &  delta_lambda                                                    &
                      ! holds spacing between points in the i
                      ! direction for the input data field.
     &, delta_phi                                                       &
                      ! holds spacing between points in the j
                      ! direction for the input data field.
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, base_lambda                                                     &
     &, base_phi

      Real, Intent(In) :: qpos_diag_limit  ! lower limit for qpos diagnostics

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
     &  glambda_p(1-halo_i : global_row_length+halo_i)                  &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, grecip_dlamp(1-halo_i : global_row_length+halo_i)               &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, lambda_p_rm(1-halo_i : row_length+halo_i)                       &
     &, lambda_p_rp(1-halo_i : row_length+halo_i)                       &
     &, phi_p_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_p_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, recip_lambda_p_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_phi_p_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )

      Real, Intent(In) ::                                               &
     &  tracer_lbcs(lenrim,tr_levels,tr_lbc_vars) ! Tracer lbcs
      Real, Intent(In) ::                                               &
     &  ukca_tracer_lbcs(lenrim,tr_levels,tr_lbc_ukca) 
                                                  ! UKCA Tracer lbcs

      Logical, Intent(In) ::                                            &
     &  L_CO2_interactive                                               &
     &, L_murk_advect                                                   &
     &, L_Soot                                                          &
     &, L_biomass                                                       &
     &, L_ocff                                                          &
     &, L_DUST                                                          &
     &, L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms                           &
     &, L_USE_CARIOLLE                                                  &
     &, L_nitrate
     

      Real, Intent(InOut) ::                                            &
     & CO2 (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &,murk(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &,soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,soot_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         model_levels)                                            &
     &,bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         model_levels)                                            &
     &,nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &         model_levels)                                            &
     &, DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, so2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &           model_levels)                                          &
     &, so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, nh3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, dms(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, tracers(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          tr_levels,tr_vars)                                      &
     &, tracers_ukca(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          tr_levels,tr_ukca)                                      &
! Add cariolle specific parameters for ozone tracer     
     &, OZONE_TRACER(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, eta_theta_levels(0:model_levels)

      Real                                                              &
     &  rho_n   (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, rho_np1 (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)

      Real                                                              &
                      ! Trig functions.
     &  FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

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

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, l                                                      &
                       ! Loop indices
     &, temp                                                            &
     &, count                                                           &
     &, tr_start                                                        &
     &, lbc_index       ! pointer to tracer for active tracer lbcs

      Logical                                                           &
     &  L_vector                                                        &
                       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.
     &, L_do_halos                                                      &
                                  ! update the halos?
     &, L_do_boundaries           ! update the boundaries?

! arrays

      Real                                                              &
     &  work(row_length, rows, model_levels)                            &
     &, work_np1(row_length, rows, model_levels)                        &
     &, drk, drkp1                                                      &
     &, work4(row_length, rows)

      integer super_array_size, array_size_count, super_array_size_qpos

      Real :: super_array(1-halo_i:row_length+halo_i,                   &
                          1-halo_j:rows+halo_j,                         &
                          model_levels, super_array_size)     
      Real :: super_tracer_phys2(row_length,                            &
                                 rows,                                  &
                                 model_levels, super_array_size)         
      Real :: tracer_phys1(1-halo_i:row_length+halo_i,                  &
                           1-halo_j:rows+halo_j,                        &
                           tr_levels, tr_vars)   
      Real :: tracer_phys2(row_length, rows, tr_levels, tr_vars) 
      Real :: ukca_tracer_phys1(1-halo_i:row_length+halo_i,             &
                                1-halo_j:rows+halo_j,                   &
                                tr_levels, tr_ukca)   
      Real :: ukca_tracer_phys2(row_length, rows, tr_levels, tr_ukca)
      Real :: data_out_super(row_length, rows, model_levels, super_array_size)

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi    (row_length, rows, model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
!  Section 0.    Initialise array_size_count
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_TRACER2',zhook_in,zhook_handle)
      array_size_count=0
!    qpos is not called for murk 
      super_array_size_qpos=super_array_size
      if (L_murk_advect)super_array_size_qpos=super_array_size_qpos-1

!qcon block start
! ----------------------------------------------------------------------
! Section 3.1  Set appropriate weighted rho*r*r*delta_r at data points
! and at time t and t+deltat.
! ----------------------------------------------------------------------

! Weights come from rewriting: Sum_k(qbarr * rho*r*r*dr) =
!                              Sum_k(q * weighted average of rho*r*r*dr)

! Store in work (for current timestep n) and work_np1 (for next
! timestep n+1) -  'work' array therefore reused

! It is assumed that rho_n holds r-squared scaled current value of rho
! and that rho_np1 holds r-squared scaled value of rho at next timestep


! Note it is assumed that q(0) = q(1) and so q(0) contribution has
! been absorbed into that of q(1), hence different form of drk.
! This is not essential part of alogorithm and could be changed.
   
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(drkp1, drk, i, j, k)

        k = 1

!$OMP DO SCHEDULE(STATIC)
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
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
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
!$OMP END DO

        k = model_levels

!$OMP DO SCHEDULE(STATIC)
        Do j = 1, rows
          Do i = 1, row_length
              drk   = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              work(i,j,k)   = wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
            work_np1(i,j,k) = wet_to_dry_np1(i,j,k)*rho_np1(i,j,k)*drk
          End Do
        End Do
!$OMP END DO

!$OMP END PARALLEL 

!qcon block ends

! ----------------------------------------------------------------------
! Section 1.0  do interpolation of tracers
! ----------------------------------------------------------------------

        L_vector = .false.

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

! super array set up in first tr_set_phys

      if(super_array_size >= 1)then

      if(l_2dcomm) then
        size_int_qc_mult =  row_length * rows * model_levels
      else
        size_int_qc_mult = model_levels * global_row_length
      end if

! DEPENDS ON: interpolation_qcon_multi
         Call Interpolation_qcon_multi(                                 &
     &                     super_array,                                 &
     &                     eta_theta_levels(1),                         &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     work, work_np1,                              &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     super_array_size, check_bottom_levels,       &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, model_levels,              &
     &                     rows,                                        &
     &                     row_length, rows, model_levels,              &
     &                     glambda_p, phi_p,                            &
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
     &                     L_vector, model_domain, L_high_moist,        &
     &                     L_mono_moist, L_conserv_moist,               &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j, global_row_length,           &
     &                     global_rows, row_length, rows,               &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     g_j_pe, l_2dcomm, size_2dcomm,               &
     &                     group_2dcomm, size_int_qc_mult,              &
     &                     gc_all_proc_group,                           &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     0, 0, L_sl_halo_reprod,                      &
     &                     off_x, off_y,                                &
     &                     data_out_super, Error_Code )
      endif

!  add on atmos_physics2 increment to readvected tracer
      Do l=1,super_array_size
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              super_tracer_phys2(i,j,k,l) = data_out_super(i,j,k,l) +   &
     &                  super_tracer_phys2(i,j,k,l)
            End Do
          End Do
        End Do
      End Do

      if (super_array_size_qpos >0)then
!    check tracer values are all non-negative
! DEPENDS ON: q_pos_ctl
        call Q_Pos_Ctl(                                                 &
                         super_tracer_phys2, row_length, rows,          &
                         model_levels*super_array_size_qpos,            &
                         super_array_size_qpos, bl_levels,              &
                         global_row_length, global_rows,                &
                         me, n_proc, 0,0,                               &
                         model_domain,                                  &
                         halo_type_no_halo, q_pos_tracer_method,  0.0,  &
                         l_qpos_diag_pr, qpos_diag_limit,               &
                         'TRACER call from sl_tracer2')
      endif      ! super_array_size_qpos >0

! ----------------------------------------------------------------------
! section 2.0  set end of timestep tracers
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.1  carbon cycle.
! ----------------------------------------------------------------------
      array_size_count=0
      if(l_CO2_interactive)then

        array_size_count=array_size_count +1

        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              co2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      endif  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 2.2  Soot cycle.
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)
      If (l_Soot) Then

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              soot_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              soot_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              soot_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
!$OMP END DO 

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER

      endif  ! l_soot

! ----------------------------------------------------------------------
! Section 2.3  Biomass aerosol.
! ----------------------------------------------------------------------
      If (l_Biomass) Then

!$OMP MASTER 
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              bmass_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              bmass_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              bmass_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER
      endif  ! l_biomass

! ----------------------------------------------------------------------
! Section 2.4  sulphur cycle.
! ----------------------------------------------------------------------
      If (l_Sulpc_so2) Then

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              so4_aitken(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              so4_accu(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              so4_diss(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
              so2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+3)
            End Do
          End Do
        End Do
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +3
!$OMP END MASTER 
!$OMP BARRIER

        if(L_sulpc_nh3)then

!$OMP MASTER 
          array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = i_start,i_end
                nh3(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
!$OMP END DO
        endif

        if(L_sulpc_dms)then

!$OMP MASTER
          array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = i_start,i_end
                dms(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
!$OMP END DO
        endif

      endif  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 2.5  mineral dust.
! ----------------------------------------------------------------------
      IF (L_DUST) THEN


!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO K = 1, MODEL_LEVELS
          Do j = j_start,j_end
            Do i = i_start,i_end
              DUST_DIV1(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              DUST_DIV2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              IF(.NOT.l_twobin_dust) THEN
                DUST_DIV3(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
                DUST_DIV4(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+3)
                DUST_DIV5(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+4)
                DUST_DIV6(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+5)
              END IF 
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        IF (l_twobin_dust) THEN
          array_size_count=array_size_count + 1
        ELSE
          array_size_count=array_size_count + 5
        END IF
!$OMP END MASTER 
!$OMP BARRIER
      ENDIF  ! L_DUST

! ----------------------------------------------------------------------
! New addition  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      If (L_OCFF) Then
!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              ocff_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              ocff_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              ocff_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER

      endif  ! l_ocff
!$OMP END PARALLEL 

! ----------------------------------------------------------------------
! Section 2.5.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (L_USE_CARIOLLE) THEN
        array_size_count=array_size_count +1
        DO K = 1, MODEL_LEVELS
          Do j = j_start,j_end
            Do i = i_start,i_end
       OZONE_TRACER(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      ENDIF  ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! Section 2.5.2  Ammonium nitrate aerosol
! ----------------------------------------------------------------------
      If (L_nitrate) Then

        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              nitr_acc(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              nitr_diss(i,j,k)= super_tracer_phys2(i,j,k,array_size_count+1)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +1
      endif  ! l_ocff

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE ^^^
!  and in the same location/order in sl_tracer1 and tr_set_phys
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_vars>0 ) THEN

        do count=1,tr_vars
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            Do j = j_start,j_end
              Do i = i_start,i_end
              tracers(i,j,k,count) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do
      End IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 2.7  UKCA tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_ukca> 0 ) THEN

        do count=1,tr_ukca
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            Do j = j_start,j_end
              Do i = i_start,i_end
              tracers_ukca(i,j,k,count) =                               &
     &           super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do

      End IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 2.99  Murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      IF (L_Murk_advect) then

        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              murk(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      End IF  ! L_Murk_advect

      IF (lhook) CALL dr_hook('SL_TRACER2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SL_tracer2
