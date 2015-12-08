! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_NI_init

      Subroutine IDL_NI_init(                                           &
     & model_domain, row_length, rows, n_rows, model_levels,wet_levels  &
     &,TR_VARS,TR_UKCA,TR_LEVELS,bl_levels,first_constant_r_rho_level   &
     &,cos_theta_latitude, sec_theta_latitude, f3_at_u, f3_at_v         &
     &,timestep, first_atmstep_call, L_regular                          &
     &,delta_x, delta_y, delta_lambda, delta_phi, base_phi, base_lambda &
     &,lat_rot_NP_in, long_rot_NP_in                                    &
     &,lambda_p, phi_p, lambda_u, phi_v, lambda_p_end, phi_p_end        &
     &,r_theta_levels, r_rho_levels, r_at_u, r_at_v, z_orog_print       &
     &,eta_theta_levels, eta_rho_levels                                 &
! Multi-processor
     &,offx,offy,halo_i,halo_j, mype, nproc, at_extremity, datastart    &
     &,gc_all_proc_group, global_row_length, global_rows                &
     &,g_rows, g_row_length, nproc_x                                    &
! Primary fields
     &,theta, rho, exner_theta_levels, exner_rho_levels                 &
     &,p, p_theta_levels, pstar, q, qcl, qcf, qcf2, qrain, qgraup       &
     &,cf, cfl, cff, u, v, w, u_adv, v_adv, w_adv                       &
! LAM lateral boundary data
     &,rimwidth, rimweights, lenrim, lbc_size, lbc_start                &
     &,theta_lbc, theta_lbc_tend, exner_lbc, exner_lbc_tend             &
     &,rho_lbc, rho_lbc_tend, q_lbc, q_lbc_tend                         &
     &,qcl_lbc, qcl_lbc_tend, qcf_lbc, qcf_lbc_tend                     &
     &,qcf2_lbc, qcf2_lbc_tend, qrain_lbc, qrain_lbc_tend               &
     &,qgraup_lbc, qgraup_lbc_tend, cf_bulk_lbc, cf_bulk_lbc_tend       &
     &,cf_liquid_lbc, cf_liquid_lbc_tend                                &
     &,cf_frozen_lbc, cf_frozen_lbc_tend                                &
     &,u_lbc, u_lbc_tend, v_lbc, v_lbc_tend, w_lbc, w_lbc_tend          &
     &,u_adv_lbc, u_adv_lbc_tend, v_adv_lbc, v_adv_lbc_tend             &
     &,w_adv_lbc, w_adv_lbc_tend                                        &
! Grid info for idealised
     &,height_domain_in, height_domain, big_layers                      &
     &,transit_layers, mod_layers, surface_type, p_surface              &
! Profile settings
     &,                  tprofile_number, qprofile_number               &
     &,                  uvprofile_number                               &
     &,                  Brunt_Vaisala                                  &
     &,                  theta_surface, dtheta_dz1, height_dz1          &
     &,                  u_in, v_in, height_u_in, ujet_lat, ujet_width  &
     &,                  u_ramp_start, u_ramp_end, f_plane, r_plane     &
     &,                  q1, num_profile_data                           &
     &,                  zprofile_data, tprofile_data, qprofile_data    &
     &,                  num_uvprofile_data, z_uvprofile_data           &
     &,                  uprofile_data, vprofile_data                   &
! Forcing settings
     &,                  max_model_levels                               &
     &,                  max_num_profile_data, max_num_force_times      &
     &,                  tforce_option, qforce_option, uvforce_option   &
     &,                  num_tforce_levels, num_tforce_times            &
     &,                  num_qforce_levels, num_qforce_times            &
     &,                  num_uvforce_levels, num_uvforce_times          &
     &,                  z_tforce_data, tforce_data                     &
     &,                  z_qforce_data, qforce_data                     &
     &,                  z_uvforce_data,uforce_data, vforce_data        &
     &,                  tforce_data_modlev, qforce_data_modlev         &
     &,                  uforce_data_modlev, vforce_data_modlev         &
! Dynamical core settings
     &,                  SuHe_pole_equ_deltaT, SuHe_static_stab         &
     &,                  base_frictional_timescale                      &
     &,                  frictional_timescale, SuHe_sigma_cutoff        &
     &,                  SuHe_level_weight, L_SH_Williamson             &
!  Horizontal function parameters
     &,                  t_horizfn_number, uv_horizfn_number            &
     &, t_horizfn_data, L_perturb_t, perturb_magnitude_t                &
     &, L_perturb_q, perturb_magnitude_q, L_perturb_correlate_tq        &
     &, L_perturb_correlate_vert, L_perturb_correlate_time              &
     &, perturb_type, perturb_height                                    &
!  Profiles for fixed lbcs and sponge zones
     &,                  u_ref, v_ref, theta_ref, exner_ref, rho_ref    &
     &,                  q_ref                                          &
     &,L_fix_orog_hgt_lbc, orog_hgt_lbc                                 &
     &,zprofile_orog, idl_interp_option, hf                             &
!  Options
     &,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2                     &
     &,L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc     &
     &,L_constant_dz, L_rotating, L_fixed_lbcs, L_polar_wind_zero       &
     &,L_wind_balance, L_rotate_winds, L_pressure_balance, L_physics    &
     &,L_dry, L_sponge                                                  &
     &,L_trivial_trigs, L_perturb, L_code_test, L_cyclone, L_baroclinic &
     &,h_print, timestep_number, h_o_actual, grow_steps                 &
     &,h_o_per_step, h_o, grid_number,  grid_flat                       &
     &,first_theta_height, thin_theta_height, big_factor, mag           &
     &,lambda_fraction, phi_fraction, half_width_x, half_width_y        &
     &,plat_size_x, plat_size_y, Witch_power                            &
     &,idl_max_num_bubbles, idl_bubble_option, idl_bubble_max           &
     &,idl_bubble_height, idl_bubble_xoffset, idl_bubble_yoffset        &
     &,idl_bubble_width, idl_bubble_depth, L_idl_bubble_saturate        &
! Other fields
     &,                  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts       &
     &,                  nproc_y, gc_proc_row_group, n_cca_lev          &
     &,IdlSurfFluxSeaOption,IdlSurfFluxSeaParams,L_flux_bc,flux_h,flux_e&
     &,L_spec_z0, z0m_scm, z0h_scm, roughlen_z0m, roughlen_z0h          &
     &,i_hour, i_minute, i_second                                       &
     &,                  problem_number, rad_hr, orography              &
     &,tstar_tile, ntiles, land_field, land_index                       &
     &,cumulus, nbdsc, ntdsc, cca, ccb, cct, cclwp, tstar               &
     &,land_sea_mask, SW_incs, LW_incs, t1_sd, q1_sd, zh                &
     &,area_cloud_fraction, ti, z0msea, ntml, u_0, v_0, u_0_p, v_0_p)

      USE earth_constants_mod, ONLY: g, earth_radius

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParParams
      USE problem_mod, ONLY: standard, monsoon, dynamical_core,         &
                             idealised_problem, standard_namelist
      IMPLICIT NONE
!
! Description:
!   this routine does the initialisation for idealised test runs.
!
! Method: take code down one level from ATMSTEP
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!


!  starting level for theta/q variables (= 0/1 for VATPOLES/ND) 
!INTEGER, PARAMETER :: stlev = 0
INTEGER, PARAMETER :: stlev = 1

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! Local number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_levels                                                      &
                         ! number of model levels where moisture
                         ! variables are held
     &, bl_levels                                                       &
                         ! number of  boundary_layer_levels
     &, n_cca_lev                                                       &
     &, first_constant_r_rho_level                                      &
     &, halo_i                                                          &
                         ! Size of halo in i direction.
     &, halo_j                                                          &
                         ! Size of halo in j direction.
     &, offx                                                            &
                         ! Size of small halo in i
     &, offy                                                            &
                         ! Size of small halo in j.
     &, TR_VARS                                                         &
     &, TR_UKCA                                                         &
     &, TR_LEVELS

      Integer :: ntiles      ! Number of tiles in MOSESII
      Integer :: land_field  ! Number of land points
      Integer :: land_index(land_field)
      Real :: tstar_tile(land_field,ntiles)

      Logical, Intent(In) :: first_atmstep_call

      Real                                                              &
     &  timestep                                                        &
     &, height_domain                                                   &
     &, orog_hgt_lbc                                                    &
     &, zprofile_orog, hf                                               &
     &, dtheta_dz1(3)                                                   &
                        !Allows different values of dtheta_dz to be set
     &, height_dz1(2)                                                   &
                        ! at different heights specified by the
                        ! height_dz variable.
     &, theta_ref(model_levels)                                         &
                                 !theta profile for use in sponge & lbcs
     &, exner_ref(model_levels + 1)                                     &
                                     ! Exner profile for use in lbcs
     &, rho_ref(model_levels)                                           &
                                ! rho profile for use in lbcs
     &, u_ref(model_levels)                                             &
                              ! u profile for use in lbcs
     &, v_ref(model_levels)   ! u_adv profile for use in lbcs

      ! Idealised perturbation settings
      Logical, Intent(In) :: L_perturb_t
      Logical, Intent(In) :: L_perturb_q
      Logical, Intent(In) :: L_perturb_correlate_tq
      Logical, Intent(In) :: L_perturb_correlate_vert
      Logical, Intent(In) :: L_perturb_correlate_time
      Real,    Intent(In) :: perturb_magnitude_t
      Real,    Intent(In) :: perturb_magnitude_q
      Real,    Intent(In) :: perturb_height(2)
      Integer, Intent(In) :: perturb_type

      ! Variables for idealised surface fluxes
      Integer, Intent(In) :: IdlSurfFluxSeaOption ! Surface flux option
      Real, Intent(In)    :: IdlSurfFluxSeaParams(10) ! Idl flux params
      Real, Intent(InOut) :: flux_h(row_length, rows) ! Idealised
      Real, Intent(InOut) :: flux_e(row_length, rows) ! surface fluxes
      Logical, Intent(InOut) :: L_flux_bc ! Switch for specified fluxes

      ! Variables for specifying roughness length
      Logical, Intent(In) :: L_spec_z0
      Real, Intent(InOut) :: z0m_scm(row_length, rows) 
!                                          ! SCM specified z0m (m)
      Real, Intent(InOut) :: z0h_scm(row_length, rows) 
!                                          ! SCM specified z0h (m)
      Real, Intent(In)    :: roughlen_z0m
      Real, Intent(In)    :: roughlen_z0h

      ! Time variables
      Integer, Intent(In) :: i_hour    ! Model time (UTC, hour)
      Integer, Intent(In) :: i_minute  ! Model time (UTC, minute)
      Integer, Intent(In) :: i_second  ! Model time (UTC, second)

      Integer                                                           &
     &  mype                                                            &
                     ! My processor number
     &, nproc                                                           &
                  ! Total number of processors
     &, nproc_y                                                         &
                    ! number of processors in N-S direction
     &, gc_all_proc_group                                               &
                          ! Group identifier for all processors.
     &, gc_proc_row_group                                               &
                          ! Group identifier for all processors.
     &, datastart(2)                                                    &
                           ! First gridpoints held by this processor
     &, global_row_length                                               &
                             ! global number of points on a row
     &, global_rows          ! global number of rows

      Integer                                                           &
     &  g_rows(0:nproc-1)                                               &
     &, g_row_length(0:nproc-1)                                         &
     &, nproc_x

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, base_phi                                                        &
     &, base_lambda                                                     &
     &, lat_rot_NP_in                                                   &
     &, long_rot_NP_in                                                  &
     &, lambda_p_end                                                    &
     &, phi_p_end

      Real                                                              &
           ! VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )

      Logical, Intent(In) :: L_regular


      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, r_at_u(1-halo_i:row_length+halo_i,                              &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, r_at_v(1-halo_i:row_length+halo_i,                              &
     &               1-halo_j:n_rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, cos_theta_latitude(1-offx:row_length+offx,                      &
     &                      1-offy:rows+offy)                           &
     &, sec_theta_latitude(1-offx:row_length+offx,                      &
     &                      1-offy:rows+offy)                           &
     &, f3_at_u (1-offx:row_length+offx,                                &
     &                      1-offy:rows+offy)                           &
     &, f3_at_v (1-offx:row_length+offx,                                &
     &                      1-offy:n_rows+offy)                         &
     &, z_orog_print(0:model_levels)


! Suarez Held variables
      Real                                                              &
     &  SuHe_pole_equ_deltaT                                            &
     &, SuHe_static_stab                                                &
     &, SuHe_sigma_cutoff                                               &
     &, base_frictional_timescale                                       &
     &, frictional_timescale(model_levels)                              &
     &, SuHe_level_weight(model_levels)

      Logical                                                           &
     &  L_SH_williamson

      Real                                                              &
     &  orography(row_length,rows)

      Real, Intent (InOut) ::                                           &
     &  rho(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, exner_theta_levels(1-offx:row_length+offx, 1-offy:rows+offy,    &
     &        model_levels)                                             &
     &, exner_rho_levels(1-offx:row_length+offx, 1-offy:rows+offy,      &
     &        model_levels+1)                                           &
     &, p_theta_levels(1-offx:row_length+offx,  1-offy:rows+offy,       &
     &        model_levels)                                             &
     &, p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &    model_levels+1)                                               &
     &, pstar(row_length, rows)                                         &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &          stlev:model_levels)                                           &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      stlev:wet_levels)                                                 &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        stlev:wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        stlev:wet_levels)                                               &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &      stlev:wet_levels)                                                 &
     &, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &      stlev:wet_levels)                                                 &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &      stlev:wet_levels)                                                 &
     &, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &      stlev:wet_levels)                                                 &
     &, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        stlev:wet_levels)                                               &
     &, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        stlev:wet_levels)                                               &
     &, u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                             &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &        model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        0:model_levels)

!set-ancil below here
      Logical                                                           &
     &  cumulus(row_length, rows)                                       &
                                  ! bl convection flag
     &, land_sea_mask(row_length, rows)

      Real                                                              &
     &  SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, ti(row_length, rows)                                            &
                              ! set equal to tstar
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, rad_hr(row_length, rows, bl_levels,2)   ! BL rad heating rates

! Convection
      Real                                                              &
     &  cca (row_length, rows, n_cca_lev)                               &
     &, cclwp(row_length, rows) ! condensed water path (KG/M**2)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

      Integer                                                           &
     &  ntml(row_length, rows)                                          &
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)

! Diagnostic variables

      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_levels)

      Real                                                              &
     &  tstar(row_length, rows)

      Real                                                              &
     &  u_0(row_length, rows)                                           &
                                ! set to zero
     &, v_0(row_length, n_rows)                                         &
                                ! set to zero
     &, u_0_p(row_length, rows)                                         &
                                  ! set to zero
     &, v_0_p(row_length, rows) ! set to zero

      Real                                                              &
     &  z0msea (row_length, rows) ! veg/qrparm.veg.rough

! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10
! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11

      ! Stash indexing arrays
      Integer len_a_ixsts
      Integer len_a_spsts
      Integer a_ixsts(len_a_ixsts)     ! stash index array
      Real    a_spsts(len_a_spsts)     ! atmos stash array

! Idealised  variables

      INTEGER :: max_model_levels     ! max no. of model levels
      INTEGER :: max_num_profile_data ! max no. levels in forcing data
      INTEGER :: max_num_force_times  ! max no. times in forcing data
      INTEGER :: tforce_option        ! Theta forcing option
      INTEGER :: qforce_option        ! Moisture forcing option
      INTEGER :: uvforce_option       ! Horizontal wind forcing option
      INTEGER :: num_tforce_levels    ! No. levels in T forcing data
      INTEGER :: num_tforce_times     ! No. times in T forcing data
      INTEGER :: num_qforce_levels    ! No. levels for Q forcing data
      INTEGER :: num_qforce_times     ! No. times in Q forcing data
      INTEGER :: num_uvforce_levels   ! No. levels for U,V forcing data
      INTEGER :: num_uvforce_times    ! No. times in U,V forcing data

      REAL :: h_o
      REAL :: h_print
      REAL :: h_o_actual   ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain_in
      REAL :: delta_x           ! Resolution at equator
      REAL :: delta_y           ! Resolution at equator
      REAL :: big_factor
      REAL :: mag
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)        ! Input values of zonal u
      REAL :: v_in(4)        ! Input values of southerly wind v
      REAL :: height_u_in(3) ! heights specified by height_u_in variable
      REAL :: u_ramp_start   ! ramping starting latitude for u
      REAL :: u_ramp_end     ! ramping ending latitude for u
      REAL :: ujet_lat     ! To specify centre lat (degrees) of jet core
      REAL :: ujet_width     ! To specify width (degrees) of jet
      REAL :: t_horizfn_data(10) ! Data values describng horizontal t fn
      REAL :: q1             ! Allows different values of moisture
      REAL :: f_plane        ! fixed latitude Coriolis term
      REAL :: ff_plane
      REAL :: r_plane        ! reference latitude for row 1 (bottom row)

      ! Namelist profile data
      REAL :: zprofile_data(max_num_profile_data)    ! heights for t,q
      REAL :: tprofile_data(max_num_profile_data)    ! theta profile
      REAL :: qprofile_data(max_num_profile_data)    ! humidity profile
      REAL :: z_uvprofile_data(max_num_profile_data) ! heights for u,v
      REAL :: uprofile_data(max_num_profile_data)    ! u-wind profile
      REAL :: vprofile_data(max_num_profile_data)    ! v-wind profile

      REAL :: cool_rate
      REAL :: q_ref(model_levels)
      REAL :: force_time_interval ! Time interval between forcing data
      ! Heights of forcing data for T,Q,UV
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      ! Forcing data arrays for T,Q,U,V
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      ! Forcing data interpolated onto model levels for T,Q,U,V
      REAL :: tforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: qforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: uforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: vforce_data_modlev(max_model_levels, max_num_force_times)

      INTEGER :: surface_type        ! idealised orography type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number    ! temperature profile option
      INTEGER :: qprofile_number    ! moisture profile option
      INTEGER :: uvprofile_number   ! wind profile option
      INTEGER :: num_profile_data   ! no of values in z,q,tprofile_data
      INTEGER :: num_uvprofile_data   ! no of values in u,vprofile_data
      INTEGER :: t_horizfn_number   ! horiz function no. for temperature
      INTEGER :: uv_horizfn_number  ! horiz function no. for wind field
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: sponge_width
      INTEGER :: sponge_width_theta
      INTEGER :: sponge_depth
      INTEGER :: timestep_number
      INTEGER :: problem_number
      INTEGER :: idl_interp_option  ! Profile interpolation option

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_u  ! Sets a constant value of u from u_in(1)
      LOGICAL :: L_constant_dz ! Sets const dtheta_dz from dtheta_dz1(1)
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_fixed_lbcs    ! Set fixed lateral boundary conditions
      LOGICAL :: L_pressure_balance  ! Geostrophically balance pressures
      LOGICAL :: L_wind_balance  ! Geostrophically balance initial winds
      LOGICAL :: L_rotate_winds  ! rotate input u,v (true) to LAM u,v
      LOGICAL :: L_polar_wind_zero   ! set u=0 on polar row
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating      ! .true. for Earth's rotation
      LOGICAL :: L_perturb       ! add random perturb. to surface theta
      LOGICAL :: L_code_test     ! User switch for testing code
      LOGICAL :: L_cyclone       ! true if cyclone simulation
      LOGICAL :: L_baroclinic    ! true if baroclinic wave simulation
      LOGICAL :: L_force
      LOGICAL :: L_physics       ! physics switch
      LOGICAL :: L_dry           ! moisture switch
      LOGICAL :: L_sponge        ! sponge switch
      ! Bubble perturbation options
      INTEGER :: idl_max_num_bubbles
      INTEGER :: idl_bubble_option(idl_max_num_bubbles)
      REAL :: idl_bubble_max(idl_max_num_bubbles)
      REAL :: idl_bubble_height(idl_max_num_bubbles)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_width(idl_max_num_bubbles)
      REAL :: idl_bubble_depth(idl_max_num_bubbles)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fix_orog_hgt_lbc

      LOGICAL :: L_mcr_qcf2    ! true if using 2nd prognostic cloud ice
      LOGICAL :: L_mcr_qrain   ! true if using prognostic rain
      LOGICAL :: L_mcr_qgraup  ! true if using prognostic graupel
      LOGICAL :: L_pc2         ! true if using PC2 cloud scheme

! Lateral boundary variables

      LOGICAL, Intent (In) ::                                           &
     &  L_mcr_qcf2_lbc                                                  &
                          ! true if prognostic 2nd cloud ice in lbcs
     &, L_mcr_qrain_lbc                                                 &
                          ! true if prognostic rain in lbcs
     &, L_mcr_qgraup_lbc                                                &
                          ! true if prognostic graupel in lbcs
     &, L_pc2_lbc         ! true if prognostic cloud fracs in lbcs

      INTEGER, Intent (In) ::                                           &
     &  rimwidth                                                        &
                             ! Size of boundary region
     &, lenrim(Nfld_max,NHalo_max)                                      &
                             ! Size of single level of LBC
     &, lbc_size(4,Nfld_max,NHalo_max)                                  &
                             ! Size of a side of a LBC
     &, lbc_start(4,Nfld_max,NHalo_max)
                             ! Start of a side in a LBC

      REAL, Intent (In) ::                                              &
     &  rimweights(rimwidth) ! Weight to apply to LBC

      Real, Intent (InOut) ::                                           &
     &  u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
     &        model_levels)                                             &
                                    ! in/out : u lbc
     &, u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
     &             model_levels)                                        &
                                    ! in : u lbc tendency
     &, v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
     &        model_levels)                                             &
                                    ! in/out : v lbc
     &, v_lbc_tend(lenrim(fld_type_v,halo_type_extended),               &
     &             model_levels)                                        &
                                    ! in : v lbc tendency
     &, w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
     &        0:model_levels)                                           &
                                    ! in/out : v lbc
     &, w_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
     &             0:model_levels)                                      &
                                    ! in : v lbc tendency
     &, rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
     &          model_levels)                                           &
                                    ! in/out : rho lbc
     &, rho_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
     &               model_levels)                                      &
                                    ! in : rho lbc tendency
     &, theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
     &            model_levels)                                         &
                                    ! in/out : theta lbc
     &, theta_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
     &                 model_levels)                                    &
                                    ! in : theta lbc tendency
     &, q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
     &        wet_levels)                                               &
                                    ! in/out : q lbc
     &, q_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
     &             wet_levels)                                          &
                                    ! in : q lbc tendency
     &, qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
     &          wet_levels)                                             &
                                    ! in/out : qcl lbc
     &, qcl_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
     &               wet_levels)                                        &
                                    ! in : qcl lbc tendency
     &, qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
     &          wet_levels)                                             &
                                    ! in/out : qcl lbc
     &, qcf_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
     &               wet_levels)                                        &
                                    ! in : qcl lbc tendency
     &, qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
     &          wet_levels)                                             &
                                    ! in/out : qcf2 lbc
     &, qcf2_lbc_tend(lenrim(fld_type_p,halo_type_extended),            &
     &               wet_levels)                                        &
                                    ! in : qcf2 lbc tendency
     &, qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
     &          wet_levels)                                             &
                                    ! in/out : qrain lbc
     &, qrain_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
     &               wet_levels)                                        &
                                    ! in : qrain lbc tendency
     &, qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
     &          wet_levels)                                             &
                                    ! in/out : qgraup lbc
     &, qgraup_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
     &               wet_levels)                                        &
                                    ! in : qgraup lbc tendency
     &, cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
     &          wet_levels)                                             &
                                    ! in/out : cf_bulk lbc
     &, cf_bulk_lbc_tend(lenrim(fld_type_p,halo_type_extended),         &
     &               wet_levels)                                        &
                                    ! in : cf_bulk lbc tendency
     &, cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
     &          wet_levels)                                             &
                                    ! in/out : cf_liquid lbc
     &, cf_liquid_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
     &               wet_levels)                                        &
                                    ! in : cf_liquid lbc tendency
     &, cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
     &          wet_levels)                                             &
                                    ! in/out : cf_frozen lbc
     &, cf_frozen_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
     &               wet_levels)                                        &
                                    ! in : cf_frozen lbc tendency
     &, exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
     &            model_levels+1)                                       &
                                    ! in/out : exner lbc
     &, exner_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
     &                 model_levels+1)                                  &
                                         ! in : exner lbc tendenc
     &, u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
     &            model_levels)                                         &
                                    ! in/out : u_adv lbc
     &, u_adv_lbc_tend(lenrim(fld_type_u,halo_type_extended),           &
     &                 model_levels)                                    &
                                    ! in : u_adv lbc tendency
     &, v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
     &            model_levels)                                         &
                                    ! in/out : v_adv lbc
     &, v_adv_lbc_tend(lenrim(fld_type_v,halo_type_extended),           &
     &                 model_levels)                                    &
                                    ! in : v_adv lbc tendency
     &, w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
     &            0:model_levels)                                       &
                                    ! in/out : w lbc
     &, w_adv_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
     &                 0:model_levels) ! in : w lbc tendency


      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error

      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if return code >0

! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='IDL_NI_INIT')

      ! Work arrays with extended halos (_eh)
      ! Needed so that external halo values in LAMS can be set correctly
      ! for the lateral boundary arrays.
      REAL                                                              &
     &  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           stlev:model_levels)                                          &
     &, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &           model_levels)                                          &
     &, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                 &
     &           1-halo_j:rows+halo_j, model_levels+1)

      ! Loop counters
      INTEGER i,j ,ij, k, k2      ! Loop counters

      REAL z_at_theta   ! height of theta level in forcing data interp
      REAL z_at_rho     ! height of rho level in forcing data interp
      REAL weight       ! weight used in forcing data interpolation
      REAL, ALLOCATABLE, DIMENSION(:,:), SAVE :: orog_per_step

      ! Variables for idealised surface fluxes
      REAL :: hrl,xfact

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_NI_INIT',zhook_in,zhook_handle)
      If (first_atmstep_call) Then

        If (mype == 0) Then
          Write (6,*) ' =============================================', &
     &                '=========================='
          Write (6,*) ' '
          Write (6,*) '                         UM IDEALISED SETTINGS'
          Write (6,*) ' '
          Write (6,*) ' =============================================', &
     &                '=========================='
        End If
      End If

!----------------------------------------------------------------------
!
!                Start of section for Timestep = 1
!    (i.e. first timestep of a new run, but NOT a continuation run)
!
!----------------------------------------------------------------------

      If (timestep_number == 1 ) Then

        If (mype == 0) Then
          IF ( L_rotating ) THEN
            Write (6,*) ' '
            Write (6,*) ' Planet is rotating '
          ELSE
            Write (6,*) ' '
            Write (6,*) ' ***   NO PLANETARY ROTATION   *** '
          END IF ! L_rotating
          Write (6,*) ' '
          Write (6,*) ' VERTICAL GRID AND OROGRAPHY '
        End If

! h_o_actual needs to be set>0 to ensure correct surface type generated
!  when not setting h_o in namelist
        h_o_actual = 10.0     ! any value > 0.0
        if ( grow_steps  >   1)then
          allocate (orog_per_step(row_length, rows))

          If (mype == 0) Then
            Write (6,*) '   Growing Orography '
          End If

          if ( surface_type   ==  surface_mask .or.                     &
     &          surface_type   ==  surface_dump ) then
! store orography amount per timestep
            Do j = 1, rows
              Do i = 1, row_length
                orog_per_step(i,j) = orography(i,j ) / grow_steps
! reset orography for first timestep
                orography(i,j )  = orog_per_step(i,j)
              End do
            End do
          else    ! growing idealised data
            h_o_per_step =  h_o / grow_steps
            h_o_actual =  h_o_per_step
            If (mype == 0) Then
              Write (6,*) '   Hill/mountain grows over ',grow_steps,    &
     &                    ' timesteps'
              Write (6,*) '   Final hill/mountain height ',h_o,' metres'
            End If   !(mype == 0)
          endif ! surface_type   ==  surface_mask .or.
!                   surface_type   ==  surface_dump
        else  ! grow_steps  <=  1
          if ( surface_type   ==  surface_mask .or.                     &
     &          surface_type   ==  surface_dump ) then
           if(mype  ==  0)then
              print*,'Orography fixed from start  '
            Endif   !(mype  ==  0)
          else  ! idealised orography
            h_o_per_step = 0.0
            h_o_actual = h_o
            if(mype  ==  0)then
              print*,'Hill/mountain fixed max height ',h_o,' metres'
            Endif   !(mype  ==  0)
           endif ! surface_type   ==  surface_mask .or.
!                    surface_type   ==  surface_dump
       endif  ! grow_steps  >   1

        if ( grid_number  /=  vert_dump ) then  ! if change grid

!   Calculate eta_levels (normalized vertical grid)
!    Only needs to be done at start
! DEPENDS ON: idl_calc_eta_levels
          Call IDL_Calc_eta_levels(                                     &
     &                      row_length, rows, model_levels              &
     &,                     mype, nproc                                 &
     &,                     first_constant_r_rho_level                  &
     &,                     eta_theta_levels, eta_rho_levels            &
!  Grid information
     &,                     grid_number, height_domain                  &
     &,                     first_theta_height, thin_theta_height       &
     &,                     big_layers, transit_layers, mod_layers      &
     &,                     big_factor, mag                             &
     &,                     L_code_test)

!    height_domain set for idealised problem may have been changed
!    in Calc_eta_levels. Anyway, dump value of model top needs changing
          height_domain_in = height_domain

        end if        ! grid_number  /=  vert_dump

!   Calculate required orography
!  Use h_print to pass in h_o_actual since overwritten if real orography
!                                        (surface_mask or surface_dump)
        h_print = h_o_actual

! DEPENDS ON: idl_surface_setup
        Call IDL_Surface_setup(                                         &
     &                    row_length, rows, model_levels                &
     &,                   global_row_length, global_rows                &
     &,                   halo_i, halo_j                                &
     &,                   mype, nproc, at_extremity, model_domain       &
     &,                   datastart, gc_all_proc_group                  &
     &,                   delta_lambda, delta_phi, Base_phi             &
     &,                   n_rows, base_lambda                           &
!  VarRes Grid Spacing
     &,                   lambda_p, phi_p, lambda_u, phi_v              &
     &,                   lambda_p_end, phi_p_end                       &
     &,                   L_regular                                     &
     &,                   delta_x, delta_y                              &
     &,                   orography, r_theta_levels                     &
     &,                   L_fix_orog_hgt_lbc, orog_hgt_lbc, rimwidth    &
     &,                   surface_type                                  &
     &,                   h_print, lambda_fraction, phi_fraction        &
     &,                   half_width_x, half_width_y                    &
     &,                   plat_size_x, plat_size_y, Witch_power         &
     &,                   L_code_test)

        if ( grow_steps  >   1)then
! For growing real orography send final max value into Generate_grid
!  so grid information over orography can be printed
          if ( surface_type   ==  surface_mask .or.                     &
     &         surface_type   ==  surface_dump ) then
            h_print = h_print * grow_steps
          endif   ! surface_type   ==  surface_mask .or.
!                   surface_type   ==  surface_dump
        endif    !  grow_steps  >   1

!   Now generate r_theta_levels and r_rho_levels
! DEPENDS ON: idl_generate_grid
        Call IDL_Generate_grid(                                         &
     &                   row_length, rows, n_rows, model_levels         &
     &,                  bl_levels                                      &
     &,                  first_constant_r_rho_level                     &
     &,                  halo_i, halo_j, mype                           &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  r_at_u, r_at_v                                 &
     &,                  eta_theta_levels, eta_rho_levels               &
     &,                  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts       &
     &,                  grid_number, grid_flat                         &
     &,                  height_domain_in                               &
     &,                  h_print, z_orog_print                          &
     &,                  L_code_test)
         print*,'Eta_grid generated problem_number ',problem_number

! Copy single halo arrays into extended halo arrays
        Do k = 1, model_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
              theta_eh(i,j,k) = theta(i,j,k)
              rho_eh(i,j,k)   = rho(i,j,k)
              exner_rho_levels_eh(i,j,k)   = exner_rho_levels(i,j,k)
            End Do
          End Do
        End Do
        k = model_levels + 1
        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
            exner_rho_levels_eh(i,j,k)   = exner_rho_levels(i,j,k)
          End Do
        End Do

!   Generate intitial profiles for fields
!   For grid_number = vert_dump only prints out grid information
! DEPENDS ON: idl_initial_data
        Call IDL_Initial_data(                                          &
     &                   model_domain, row_length, rows, n_rows         &
     &,                  model_levels, wet_levels                       &
     &,                  TR_VARS, TR_UKCA, TR_LEVELS, bl_levels         &
     &,                  first_constant_r_rho_level                     &
     &,                  cos_theta_latitude, sec_theta_latitude         &
     &,                  f3_at_u, f3_at_v, timestep                     &
     &,                  delta_x, delta_y                               &
     &,                  offx, offy, halo_i, halo_j                     &
     &,                  mype, nproc, at_extremity                      &
     &,                  datastart, gc_all_proc_group                   &
     &,                  global_row_length, global_rows                 &
     &,                  delta_lambda, delta_phi                        &
!  VarRes Grid Spacing
     &,                  lambda_p, phi_p, lambda_u, phi_v               &
     &,                  L_regular                                      &
     &,                  base_phi, base_lambda                          &
     &,                  lat_rot_NP_in, long_rot_NP_in                  &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  r_at_u, r_at_v, z_orog_print                   &
     &,                  eta_theta_levels, eta_rho_levels               &
     &,                  theta_eh, rho_eh, exner_rho_levels_eh          &
     &,                  q, qcl, qcf, qcf2, qrain, qgraup               &
     &,                  u_adv, v_adv, w_adv                            &
     &,                  g_rows, g_row_length                           &
     &,                  nproc_x, nproc_y                               &
!  Grid information
     &,                  height_domain_in                               &
     &,                  big_layers, transit_layers, mod_layers         &
     &,                  surface_type, p_surface                        &
! Profile settings
     &,                  tprofile_number, qprofile_number               &
     &,                  uvprofile_number                               &
     &,                  Brunt_Vaisala                                  &
     &,                  theta_surface, dtheta_dz1, height_dz1          &
     &,                  u_in, v_in, height_u_in, ujet_lat, ujet_width  &
     &,                  u_ramp_start, u_ramp_end, f_plane, r_plane     &
     &,                  q1, max_num_profile_data, num_profile_data     &
     &,                  zprofile_data, tprofile_data, qprofile_data    &
     &,                  num_uvprofile_data, z_uvprofile_data           &
     &,                  uprofile_data, vprofile_data                   &
     &,                  zprofile_orog, idl_interp_option, hf           &
! Dynamical core settings
     &,                  SuHe_pole_equ_deltaT, SuHe_static_stab         &
     &,                  base_frictional_timescale                      &
     &,                  frictional_timescale, SuHe_sigma_cutoff        &
     &,                  SuHe_level_weight, L_SH_Williamson             &
!  Horizontal function parameters
     &,                  t_horizfn_number, uv_horizfn_number            &
     &,                  t_horizfn_data                                 &
     &,                  L_perturb_t, perturb_magnitude_t               &
     &,                  L_perturb_q, perturb_magnitude_q               &
     &,                  L_perturb_correlate_tq                         &
     &,                  L_perturb_correlate_vert                       &
     &,                  L_perturb_correlate_time                       &
     &,                  perturb_type, perturb_height                   &
!  Profiles for fixed lbcs and sponge zones
     &,                  u_ref, v_ref, theta_ref, exner_ref, rho_ref    &
     &,                  q_ref                                          &
!  Options
     &,                  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup          &
     &,                  L_constant_dz, L_rotating                      &
     &,                  L_fixed_lbcs, L_polar_wind_zero                &
     &,                  L_wind_balance, L_rotate_winds                 &
     &,                  L_pressure_balance, L_physics, L_dry, L_sponge &
     &,                  L_cyclone,  L_baroclinic                       &
     &,                  idl_max_num_bubbles, idl_bubble_option         &
     &,                  idl_bubble_max, idl_bubble_height              &
     &,                  idl_bubble_xoffset, idl_bubble_yoffset         &
     &,                  idl_bubble_width, idl_bubble_depth             &
     &,                  L_idl_bubble_saturate                          &
     &,                  L_trivial_trigs, L_code_test)

        ! Copy idealised data from extended halo work arrays into
        ! model fields and set other fields.

! DEPENDS ON: idl_ideal_set_fields
        Call IDL_ideal_set_fields(                                      &
     &              model_domain, row_length, rows, n_rows              &
     &,             model_levels, wet_levels                            &
     &,             offx, offy, halo_i, halo_j                          &
     &,             mype, nproc, at_extremity                           &
     &,             r_theta_levels, r_rho_levels                        &
     &,             theta_eh, rho_eh, exner_rho_levels_eh               &
     &,             u,v,w                                               &
     &,             u_adv, v_adv, w_adv                                 &
     &,             theta                                               &
     &,             q, qcl, qcf, rho                                    &
     &,             exner_rho_levels                                    &
     &,             exner_theta_levels                                  &
     &,             p, p_theta_levels, pstar                            &
     &,             .true.)

        If(model_domain  ==  mt_global)then
! DEPENDS ON: polar_reset_mean
          Call Polar_Reset_Mean(                                        &
     &                      exner_rho_levels,rho,theta,w,               &
     &                      q,qcl,qcf,                                  &
     &                      cf,cfl,cff,                                 &
     &                      row_length, rows, model_levels,             &
     &                      wet_levels, global_row_length,              &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_y, gc_proc_row_group,          &
     &                      at_extremity)
        endif

          ! Set up ancillary data fields if physics is on
        If (L_Physics) Then

! DEPENDS ON: idl_set_ancil
          Call IDL_Set_Ancil(                                           &
     &               row_length, rows, n_rows, model_levels             &
     &,              wet_levels, n_cca_lev                              &
     &,              halo_i, halo_j                                     &
     &,              cumulus, nbdsc, ntdsc                              &
     &,              cca, ccb, cct, cclwp                               &
     &,              tstar, land_sea_mask                               &
     &,              SW_incs, LW_incs                                   &
     &,              t1_sd, q1_sd, zh                                   &
     &,              area_cloud_fraction,  cf,cfl,cff                   &
     &,              ti, z0msea, ntml, u_0, v_0, u_0_p, v_0_p           &
     &,              theta_eh, tstar_tile, ntiles                       &
     &,              land_field, land_index                             &
     &,              L_code_test, rad_hr, theta_surface, bl_levels )

          End If ! on (L_Physics)

      End If ! on timestep = 1
!----------------------------------------------------------------------
!                End of section for Timestep = 1
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!
!              Start of section for first call to ATMSTEP
!         (i.e. first timestep of a new run OR continuation run)
!
!----------------------------------------------------------------------

      If (first_atmstep_call) Then

        !-----------------------------------------------------------
        ! Interpolate theta forcing data to model levels if required
        !-----------------------------------------------------------

        If (tforce_option  >   0) Then

          ! Check to make sure the namelist profile data extends
          ! to the top of the model.
          If (eta_theta_levels(model_levels)*height_domain_in           &
     &         >   z_tforce_data(num_tforce_levels)) Then
            Write(Cmessage,*)                                           &
     &        'Idealised namelist forcing profile data (T)'             &
     &        //'does not extend to the top of the model.'              &
     &        //'Please modify the namelist data.'
            ErrorStatus = 1

            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          ! Interpolate theta from namelist profile to model levels
          Do k = 1, model_levels
            z_at_theta = eta_theta_levels(k) * height_domain_in
            Do k2 = 1, num_tforce_levels-1
              Do j = 1, num_tforce_times

                  If (z_at_theta  >   z_tforce_data(k2) .and.           &
     &                z_at_theta  <=  z_tforce_data(k2+1)) Then

                    weight = (z_at_theta - z_tforce_data(k2))           &
     &                     /(z_tforce_data(k2+1) - z_tforce_data(k2))

                    tforce_data_modlev(k,j) = tforce_data(k2,j)         &
     &               + weight*(tforce_data(k2+1,j) - tforce_data(k2,j))
                  End If

              End Do
            End Do
          End Do

        End If

        !-----------------------------------------------------------
        ! Interpolate q forcing data to model levels if required
        !-----------------------------------------------------------

        If (qforce_option  >   0) Then

          ! Check to make sure the namelist profile data extends
          ! to the top of the model.
          If (eta_theta_levels(model_levels)*height_domain_in           &
     &         >   z_qforce_data(num_qforce_levels)) Then
            Write(Cmessage,*)                                           &
     &        'Idealised namelist forcing profile data (q)'             &
     &        //'does not extend to the top of the model.'              &
     &        //'Please modify the namelist data.'
            ErrorStatus = 1

            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          ! Interpolate humidity from namelist profile to model levels
          Do k = 1, model_levels
            z_at_theta = eta_theta_levels(k) * height_domain_in
            Do k2 = 1, num_qforce_levels-1
              Do j = 1, num_qforce_times

                  If (z_at_theta  >   z_qforce_data(k2) .and.           &
     &                z_at_theta  <=  z_qforce_data(k2+1)) Then

                    weight = (z_at_theta - z_qforce_data(k2))           &
     &                       /(z_qforce_data(k2+1) - z_qforce_data(k2))

                    qforce_data_modlev(k,j) = qforce_data(k2,j)         &
     &               + weight*(qforce_data(k2+1,j) - qforce_data(k2,j))
                  End If

              End Do
            End Do
          End Do

        End If

        !-----------------------------------------------------------
        ! Interpolate u,v forcing data to model levels if required
        !-----------------------------------------------------------

        If (uvforce_option  >   0) Then

          ! Check to make sure the namelist profile data extends
          ! to the top of the model.
          If (eta_theta_levels(model_levels)*height_domain_in           &
     &         >   z_uvforce_data(num_uvforce_levels)) Then
            Write(Cmessage,*)                                           &
     &        'Idealised namelist forcing profile data (uv)'            &
     &        //'does not extend to the top of the model.'              &
     &        //'Please modify the namelist data.'
            ErrorStatus = 1

            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          ! Interpolate wind from namelist profile to model levels
          Do k = 1, model_levels
            z_at_rho = eta_rho_levels(k) * height_domain_in
            Do k2 = 1, num_uvforce_levels-1
              Do j = 1, num_uvforce_times

                  If (z_at_rho  >   z_uvforce_data(k2) .and.            &
     &                z_at_rho  <=  z_uvforce_data(k2+1)) Then

                    weight = (z_at_rho - z_uvforce_data(k2))            &
     &                     /(z_uvforce_data(k2+1) - z_uvforce_data(k2))

                    uforce_data_modlev(k,j) = uforce_data(k2,j)         &
     &               + weight*(uforce_data(k2+1,j) - uforce_data(k2,j))
                    vforce_data_modlev(k,j) = vforce_data(k2,j)         &
     &               + weight*(vforce_data(k2+1,j) - vforce_data(k2,j))
                  End If

              End Do
            End Do
          End Do

        End If

        !-----------------------------------------------------
        !          Set lbcs for idealised LAM
        !------------------------------------------------------

        If ((model_domain == mt_lam) .and. L_fixed_lbcs) Then

          ! Copy initial model field boundary data into LAM lbcs
          ! Set theta_lbc, exner_lbc and rho_lbc from extended halo
          ! work arrays. Set u_lbc, v_lbc, w_lbc from extended halo
          ! u_adv, v_adv and w_adv. If L_fixed_lbcs=.true. then
          ! LAM lbcs will be fixed for the whole run

! DEPENDS ON: idl_fix_lam_lbcs
          Call IDL_Fix_Lam_LBCs(                                        &
     &     row_length,rows,n_rows,model_levels,wet_levels               &
     &,    halo_i,halo_j,at_extremity                                   &
     &,    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2                 &
     &,    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc &
     &,    rimwidth, rimweights, lenrim, lbc_size, lbc_start            &
     &,    theta_lbc, q_lbc, qcl_lbc, qcf_lbc                           &
     &,    qcf2_lbc, qrain_lbc, qgraup_lbc                              &
     &,    cf_bulk_lbc, cf_liquid_lbc, cf_frozen_lbc                    &
     &,    rho_lbc, exner_lbc                                           &
     &,    u_lbc, v_lbc, w_lbc, u_adv_lbc, v_adv_lbc, w_adv_lbc         &
     &,    theta_eh, rho_eh, exner_rho_levels_eh                        &
     &,    q, qcl, qcf, qcf2, qrain, qgraup                             &
     &,    cf, cfl, cff, u_adv, v_adv, w_adv                            &
     &     )

          ! Copy LBC data into LBC_TEND. Both remain fixed
          ! throughout the run if no idealised lbc forcing

! DEPENDS ON: idl_copy_lbc_tend
          Call IDL_COPY_LBC_TEND(                                       &
     &     lenrim                                                       &
     &,    model_levels,wet_levels                                      &
     &,    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2                 &
     &,    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc &
     &,    u_lbc, u_lbc_tend, v_lbc, v_lbc_tend, w_lbc, w_lbc_tend      &
     &,    rho_lbc, rho_lbc_tend, theta_lbc, theta_lbc_tend             &
     &,    q_lbc,q_lbc_tend, qcl_lbc,qcl_lbc_tend, qcf_lbc,qcf_lbc_tend &
     &,    qcf2_lbc, qcf2_lbc_tend, qrain_lbc, qrain_lbc_tend           &
     &,    qgraup_lbc, qgraup_lbc_tend, cf_bulk_lbc, cf_bulk_lbc_tend   &
     &,    cf_liquid_lbc, cf_liquid_lbc_tend                            &
     &,    cf_frozen_lbc, cf_frozen_lbc_tend                            &
     &,    exner_lbc, exner_lbc_tend, u_adv_lbc, u_adv_lbc_tend         &
     &,    v_adv_lbc, v_adv_lbc_tend, w_adv_lbc, w_adv_lbc_tend         &
     &     )

        End If ! on mt_lam and L_fixed_lbcs

      End If ! on first_atmstep_call = .true.
!----------------------------------------------------------------------
!              End of section for first call to ATMSTEP
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!
!              Start of section for timestep > 0
!
!----------------------------------------------------------------------
      If ( timestep_number > 0) Then

        !--------------------------------------------------------------
        !
        !                   Growing Orography
        !
        !--------------------------------------------------------------
        If ( timestep_number <= grow_steps ) Then

!  If growing orography then need to recalculate surface and
!   r vertical coordinate values (but eta values remain fixed)
          if ( surface_type   ==  surface_mask .or.                     &
     &         surface_type   ==  surface_dump ) then
            h_o_actual = h_print
! grow orography by  orog_per_step
            Do j = 1, rows
              Do i = 1, row_length
                ij = i + (j-1)*row_length
                orography(i,j) = orography(i,j ) +                      &
     &                                    orog_per_step(i,j)
              End do
            End do
            if(mype  ==  0)then
              print*,'Growing real orography '
            Endif   !(mype  ==  0)
          else  !  growing idealised orography
            h_o_actual = h_o_actual + h_o_per_step
            if(mype  ==  0)then
               print*,'Hill/mountain growing ',h_o_actual ,' metres'
            Endif   !(mype  ==  0)
          endif !  surface_type   ==  surface_mask .or.
!                    surface_type   ==  surface_dump

! DEPENDS ON: idl_surface_setup
          Call IDL_Surface_setup(                                       &
     &                    row_length, rows, model_levels                &
     &,                   global_row_length, global_rows                &
     &,                   halo_i, halo_j                                &
     &,                   mype, nproc, at_extremity, model_domain       &
     &,                   datastart, gc_all_proc_group                  &
     &,                   delta_lambda, delta_phi, Base_phi             &
     &,                   n_rows, base_lambda                           &
!  VarRes Grid Spacing
     &,                   lambda_p, phi_p, lambda_u, phi_v              &
     &,                   lambda_p_end, phi_p_end                       &
     &,                   L_regular                                     &
     &,                   delta_x, delta_y                              &
     &,                   orography, r_theta_levels                     &
     &,                   L_fix_orog_hgt_lbc, orog_hgt_lbc, rimwidth    &
     &,                   surface_type                                  &
     &,                   h_o_actual, lambda_fraction, phi_fraction     &
     &,                   half_width_x, half_width_y                    &
     &,                   plat_size_x, plat_size_y, Witch_power         &
     &,                   L_code_test)

!   Regenerate r_theta_levels and r_rho_levels
! DEPENDS ON: idl_generate_grid
          Call IDL_Generate_grid(                                       &
     &                    row_length, rows, n_rows, model_levels        &
     &,                   bl_levels                                     &
     &,                   first_constant_r_rho_level                    &
     &,                   halo_i, halo_j, mype                          &
     &,                   r_theta_levels, r_rho_levels                  &
     &,                   r_at_u, r_at_v                                &
     &,                   eta_theta_levels, eta_rho_levels              &
     &,                   a_ixsts,len_a_ixsts, a_spsts,len_a_spsts      &
     &,                   grid_number, grid_flat                        &
     &,                   height_domain_in                              &
     &,                   h_o_actual, z_orog_print                      &
     &,                   L_code_test)

          if(L_fixed_lbcs)then
            if(mype  ==  0)then
              print*,'WARNING  WARNING   Growing orography '
              print*,'It is possible that implied lbc orography is '
              print*,'inconsistent with interior orography due to '
              print*,'asymptotic behaviour in definition (e.g. Witch)'
              print*,'Call routine to make fixed lbcs consistent'
            Endif   !(mype  ==  0)
          endif     ! L_fixed_lbcs

        End If ! on (timestep_number <= grow_steps)


        !--------------------------------------------------------------
        !
        !                 Idealised Surface Fluxes
        !
        !--------------------------------------------------------------
        ! Overwrite sea surface fluxes with idealised fluxes
        ! if idealised options turned on (IdlSurfFluxSeaOption  /=  0)

        !---------------------------------
        ! Option 1: Zero surface fluxes
        !---------------------------------
        If (IdlSurfFluxSeaOption == 1) Then

          If (mype == 0 .and. first_atmstep_call) Then
            Write(6,*) ' '
            Write(6,*) '  Specifing zero sea surface heat fluxes'
            Write(6,*) ' '
          End If

          L_flux_bc = .True.

          Do j = 1, rows
            Do i = 1, row_length
              flux_h(i,j) = 0.0
              flux_e(i,j) = 0.0
            End Do
          End Do

        End If ! on (IdlSurfFluxSeaOption  ==  1)

        !---------------------------------
        ! Option 2: Diurnal cycle
        ! (positive surface fluxes during the day, zero at night)
        !---------------------------------
        If (IdlSurfFluxSeaOption == 2) Then

          If (mype == 0 .and. first_atmstep_call) Then
            Write(Unit=6,Fmt=*) ' SURFACE FLUXES'
            Write(Unit=6,Fmt=*)                                         &
     &        '  Specifing sea surface flux diurnal cycle'
            Write(Unit=6,Fmt='(A37,F7.2)')                              &
     &        '   Maximum sensible heat flux (W/m2):',                  &
     &                    IdlSurfFluxSeaParams(1)
            Write(Unit=6,Fmt='(A37,F7.2)')                              &
     &        '   Maximum latent heat flux (W/m2):  ',                  &
     &                    IdlSurfFluxSeaParams(2)
            Write(Unit=6,Fmt='(A37,F7.2)')                              &
     &        '   Time (UTC) of max flux (hours):   ',                  &
     &                    IdlSurfFluxSeaParams(3)
            Write(Unit=6,Fmt='(A37,F7.2)')                              &
     &        '   Length of the day (hours):        ',                  &
     &                    IdlSurfFluxSeaParams(4)
            Write(Unit=6,Fmt=*) ' '
          End If

          L_flux_bc = .True.

          ! IdlSurfFluxSeaParams(1) = max sensible heat flux
          ! IdlSurfFluxSeaParams(2) = max latent heat flux
          ! IdlSurfFluxSeaParams(3) = Time (UTC) of max flux (hours)
          ! IdlSurfFluxSeaParams(4) = Length of day (hours)

          ! Calculate current time (UTC) in hours
          hrl   = ((i_hour*60.0 + i_minute)*60.0 + i_second)/3600.0

          ! Set up diurnally varying function
          xfact  = COS( Pi/2.*(IdlSurfFluxSeaParams(3)-hrl)/            &
     &                        (IdlSurfFluxSeaParams(4)/2.)  )

          ! Limit fluxes to being positive (upward)
          If (xfact <= 0.) xfact=0.

          ! Set diurnally varying fluxes
          Do j = 1, rows
            Do i = 1, row_length
              flux_h(i,j) = IdlSurfFluxSeaParams(1)*xfact**1.5
              flux_e(i,j) = IdlSurfFluxSeaParams(2)*xfact**1.3
            End Do
          End Do

        End If ! on (IdlSurfFluxSeaOption  ==  2)
        
        !---------------------------------
        ! Option 3: Fixed surface fluxes
        !---------------------------------
        If (IdlSurfFluxSeaOption == 3) Then

          If (mype == 0 .and. first_atmstep_call) Then
            Write(6,*) ' '
            Write(6,*) '  Specifing fixed sea surface heat fluxes'
          End If

          L_flux_bc = .True.

          Do j = 1, rows
            Do i = 1, row_length
              flux_h(i,j) = IdlSurfFluxSeaParams(1)
              flux_e(i,j) = IdlSurfFluxSeaParams(2)
             End Do
          End Do

          If (mype == 0 .and. first_atmstep_call) Then
            Write(6,*) ' '
            Write(6,*) '   Sensible heat flux (Wm-2) = ',               &
     &                     IdlSurfFluxSeaParams(1)
            Write(6,*) '   Latent heat flux (Wm-2)   = ',               &
     &                     IdlSurfFluxSeaParams(2)
            Write(6,*) ' '
          End If

        End If ! on (IdlSurfFluxSeaOption == 3)

        !--------------------------------------------------------------
        !
        !                 Fixing roughness length
        !
        !--------------------------------------------------------------

        If (L_spec_z0) then

          Do j = 1, rows
            Do i = 1, row_length
              z0m_scm(i,j) = roughlen_z0m 
              z0h_scm(i,j) = roughlen_z0h 
            End Do
          End Do

        End If    !L_spec_z0

      End If      ! timestep_number > 0
!----------------------------------------------------------------------
!              End of section for timestep > 0
!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_NI_INIT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_NI_init
      
