! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Initial_data_3a

      Subroutine IDL_Initial_data_3a(                                   &
     &                      model_domain, row_length, rows, n_rows      &
     &,                     model_levels, wet_model_levels              &
     &,                     TR_VARS, TR_UKCA, TR_LEVELS                 &
     &,                     boundary_layer_levels                       &
     &,                     first_constant_rho_level                    &
     &,                     cos_theta_latitude, sec_theta_latitude      &
     &,                     f3_at_u, f3_at_v, timestep                  &
     &,                     delta_x, delta_y                            &
     &,                     off_x, off_y, halo_i, halo_j                &
     &,                     me, n_proc, at_extremity                    &
     &,                     l_datastart, all_proc_group                 &
     &,                     global_row_length, global_rows              &
     &,                     delta_lambda, delta_phi                     &
!  VarRes Grid Spacing
     &,                    lambda_p, phi_p, lambda_u, phi_v             &
     &,                    L_regular                                    &
     &,                     Base_phi, base_lambda                       &
     &,                     lat_rot_NP_in, long_rot_NP_in               &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     r_at_u, r_at_v, z_orog_print                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta_eh, rho_eh, exner_rho_levels_eh       &
     &,                     q, qcl, qcf, qcf2, qrain, qgraup            &
     &,                     u_adv, v_adv, w_adv                         &
     &,                     g_rows, g_row_length                        &
     &,                     nproc_x, nproc_y                            &
!  Grid information
     &,                  height_domain                                  &
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

! Purpose:
!          Sets up initial data for idealised problems.
!          Initial temperature profile at all points is set equal to
!          the equilibrium tmeperature profile at the equator.
!

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE earth_constants_mod, ONLY: g, earth_radius


      USE atmos_constants_mod, ONLY:                                    &
          r, cp, recip_epsilon, p_zero, recip_kappa 

      USE conversions_mod, ONLY: pi_over_180 
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParParams
      USE qprofile_mod, ONLY: qp_dry, qp_qsat, qp_namelist_rh, qp_namelist, &
                              qp_dump
      USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv,        &
                              tp_bv_isoth, tp_dyn_core, tp_dyn_core_lam,    &
                              tp_namelist, tp_dump
      USE uvhoriz_mod, ONLY: uv_horiz_const, uv_horiz_ramp,                 &
                             uv_horiz_balance, uv_horiz_deform,             &
                             uv_vert_const, uv_vert_interp,                 &
                             uv_vert_namelist, uv_vert_dump
      IMPLICIT NONE


!  starting level for theta/q variables (= 0/1 for VATPOLES/ND) 
INTEGER, PARAMETER :: stlev = 0
!INTEGER, PARAMETER :: stlev = 1

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
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, boundary_layer_levels                                           &
                                ! number of  boundary_layer_levels
     &, first_constant_rho_level                                        &
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, off_x                                                           &
                             ! Size of small halo in i
     &, off_y                                                           &
                             ! Size of small halo in j.
     &, TR_VARS                                                         &
     &, TR_UKCA                                                         &
     &, TR_LEVELS

      Integer, Intent (In) :: max_num_profile_data

      Real                                                              &
     &  timestep                                                        &
     &, theta_surface                                                   &
     &, p_surface                                                       &
     &, height_domain                                                   &
     &, zprofile_orog                                                   &
                       ! Height of orog of initial/forcing profile
     &, hf                                                              &
     &, u_in(4)                                                         &
                   ! Input values of zonal u
     &, v_in(4)                                                         &
                   ! Input values of southerly wind v
     &, ujet_lat                                                        &
                   ! To specify centre latitude (degrees) of jet core
     &, ujet_width                                                      &
                   ! To specify width (degrees) of jet
     &, height_u_in(3)                                                  &
                        ! heights specified by the height_u_in variable
     &, dtheta_dz1(3)                                                   &
                        !Allows different values of dtheta_dz to be set
     &, height_dz1(2)                                                   &
                        ! at different heights specified by the
                        ! height_dz variable.
     &, Brunt_Vaisala                                                   &
     &, t_horizfn_data(10)                                              &
                            ! Data values describing horizontal t fn
     &, q1                                                              &
                           ! Allows different values of moisture
     &, zprofile_data(max_num_profile_data)                             &
                                               ! heights for t,q
     &, tprofile_data(max_num_profile_data)                             &
                                               ! theta profile
     &, qprofile_data(max_num_profile_data)                             &
                                               ! humidity profile
     &, z_uvprofile_data(max_num_profile_data)                          &
                                               ! heights for u,v
     &, uprofile_data(max_num_profile_data)                             &
                                               ! u-wind profile
     &, vprofile_data(max_num_profile_data)                             &
                                               ! v-wind profile
     &, theta_ref(model_levels)                                         &
                                 !theta profile for use in sponge & lbcs
     &, exner_ref(model_levels + 1)                                     &
                                     ! Exner profile for use in lbcs
     &, rho_ref(model_levels)                                           &
                                ! rho profile for use in lbcs
     &, u_ref(model_levels)                                             &
                              ! u profile for use in lbcs
     &, v_ref(model_levels)                                             &
                              ! u_adv profile for use in lbcs
     &, q_ref(wet_model_levels)                                         &
                                  ! q profile for use in lbcs
     &, delta_x, delta_y                                                &
                              ! Resolution at equator
     &, u_ramp_start                                                    &
                        ! ramping starting latitude for u
     &, u_ramp_end                                                      &
                        ! ramping ending latitude for u
     &, r_plane                                                         &
                      ! reference latitude for row 1 (bottom row)
     &, f_plane                                                         &
                      ! fixed latitude Coriolis term
     &, idl_bubble_max(idl_max_num_bubbles)                             &
     &, idl_bubble_height(idl_max_num_bubbles)                          &
     &, idl_bubble_xoffset(idl_max_num_bubbles)                         &
     &, idl_bubble_yoffset(idl_max_num_bubbles)                         &
     &, idl_bubble_width(idl_max_num_bubbles)                           &
     &, idl_bubble_depth(idl_max_num_bubbles)

      Integer                                                           &
     &  tprofile_number                                                 &
                           ! temperature profile option
     &, qprofile_number                                                 &
                           ! moisture profile option
     &, uvprofile_number                                                &
                           ! u,v wind profile option
     &, num_profile_data                                                &
                           ! number of values in z,q,tprofile_data
     &, num_uvprofile_data                                              &
                             ! number of values in u,vprofile_data
     &, t_horizfn_number                                                &
                           ! horizontal function no. for temperature
     &, uv_horizfn_number                                               &
                           ! horizontal function no. for wind field
     &, idl_interp_option                                               &
                           ! Profile interpolation option
     &, big_layers, transit_layers, mod_layers                          &
     &, surface_type                                                    &
                           ! idealised orography type
     &, idl_max_num_bubbles                                             &
                            ! Maximum number of bubbles allowed
     &, idl_bubble_option(idl_max_num_bubbles)  ! Bubble option number

      Logical                                                           &
     &  L_constant_dz                                                   &
                       ! Sets constant dtheta_dz from dtheta_dz1(1)
     &, L_fixed_lbcs                                                    &
                       ! Set fixed lateral boundary conditions
     &, L_wind_balance                                                  &
                       ! Geostrophically balance initial winds
     &, L_pressure_balance                                              &
                           ! Geostrophically balance pressures
     &, L_rotate_winds                                                  &
                         ! rotate input u,v (true) to LAM u,v
     &, L_polar_wind_zero                                               &
                            ! set u=0 on polar row
     &, L_rotating                                                      &
                       ! Planet rotation
     &, L_physics                                                       &
                       ! physics switch
     &, L_dry                                                           &
                       ! moisture switch
     &, L_sponge                                                        &
                       ! sponge switch
     &, L_trivial_trigs                                                 &
                           !  makes grid Cartesian if .true.
     &, L_code_test                                                     &
                       ! user switch
     &, L_mcr_qcf2                                                      &
                     ! true if using second prognostic cloud ice
     &, L_mcr_qrain                                                     &
                     ! true if using prognostic rain
     &, L_mcr_qgraup                                                    &
                     ! true if using prognostic graupel
     &, L_cyclone                                                       &
                     ! true if cyclone simulation
     &, L_baroclinic                                                    &
                     ! true if baroclinic wave simulation
     &, L_idl_bubble_saturate(idl_max_num_bubbles) ! Saturate bubble

      ! Idealised perturbation settings
      Logical, Intent(In) :: L_perturb_t
      Logical, Intent(In) :: L_perturb_q
      Logical, Intent(In) :: L_perturb_correlate_tq
      Logical, Intent(In) :: L_perturb_correlate_vert
      Logical, Intent(In) :: L_perturb_correlate_time
      Real,    Intent(In) :: perturb_magnitude_t
      Real,    Intent(In) :: perturb_magnitude_q
 
      Integer, Intent(In) :: perturb_type

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc                                                          &
                   ! Total number of processors
     &, all_proc_group                                                  &
                       ! Group identifier for all processors.
     &, l_datastart(2)                                                  &
                             ! First gridpoints held by this processor
     &, global_row_length                                               &
                             ! global number of points on a row
     &, global_rows          ! global number of rows

      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, nproc_x, nproc_y                                                &
     &, ICODE,d,info

      Real  rbuf(1-halo_i:global_row_length+halo_i,                     &
     &           1-halo_j:global_rows+halo_j)

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
     &, long_rot_NP_in

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
     &, cos_theta_latitude(1-off_x:row_length+off_x,                    &
     &                      1-off_y:rows+off_y)                         &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                      1-off_y:rows+off_y)                         &
     &, f3_at_u (1-off_x:row_length+off_x,                              &
     &                      1-off_y:rows+off_y)                         &
     &, f3_at_v (1-off_x:row_length+off_x,                              &
     &                      1-off_y:n_rows+off_y)                       &
     &, z_orog_print(0:model_levels)

      Real, Intent (InOut) ::                                           &
     &  u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                             &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &        model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        0:model_levels)                                           &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    stlev:wet_model_levels)                                             &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      stlev:wet_model_levels)                                           &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      stlev:wet_model_levels)                                           &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &      stlev:wet_model_levels)                                           &
     &, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &      stlev:wet_model_levels)                                           &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &      stlev:wet_model_levels)

      ! Work arrays with extended halos (_eh)
      ! Needed so that external halo values in LAMS can be set correctly
      ! for the lateral boundary arrays.
      REAL, Intent (InOut) ::                                           &
     &  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           stlev:model_levels)                                          &
     &, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &           model_levels)                                          &
     &, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                 &
     &           1-halo_j:rows+halo_j, model_levels+1)


     Real,    Intent(InOut) :: perturb_height(2)

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

! local variables
      Integer                                                           &
     &  i, j, k, k2, fcrl                                               &
     &, gi, gj                                                          &
     &, j_start, j_end                                                  &
     &, u_field_size                                                    &
     &, row_min                                                         &
     &, v_field_size                                                    &
     &, length

      Real                                                              &
        recip_delta_lambda                                              &
     &, recip_delta_phi                                                 &
     &, BV_squared_over_g                                               &
     &, foverRT                                                         &
     &, goverRT                                                         &
     &, goverCpT                                                        &
     &, delta_z                                                         &
     &, x,y                                                             &
     &, temp                                                            &
     &, temp1                                                           &
     &, temp2                                                           &
     &, rtemp                                                           &
     &, weight                                                          &
     &, p_test                                                          &
     &, cos_ramp_start                                                  &
     &, cos_ramp_end                                                    &
     &, z_at_theta                                                      &
     &, z_at_u                                                          &
                ! height above mean sea level at u points
     &, z_at_v                                                          &
                ! height above mean sea level at v points
     &, theta_mid, theta_mid2                                           &
     &, sigma_ref(model_levels)                                         &
     &, sigma_to_kappa(model_levels)

! Williamson variables1. 1. parameters.

      Real                                                              &
     &  p_d                                                             &
     &, p_pl                                                            &
     &, lapse_rate_d                                                    &
     &, lapse_rate_i                                                    &
     &, delta_phi_0                                                     &
     &, A_will                                                          &
     &, phi_0                                                           &
     &, T_0

! 2. derived variables
      Real                                                              &
     &  p_i                                                             &
     &, p_eq                                                            &
     &, power_d                                                         &
     &, power_i                                                         &
     &, latitude                                                        &
     &, p_lim                                                           &
     &, minusgoverRTref

       Real                                                             &
     &  eta_model, z_at_orog, hs                                        &
     &, eta_profile(num_profile_data)

      Real                                                              &
     &  exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)          &
     &, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work3(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work4(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work5(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, work6(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, work7(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, work8(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, rot_coeffu1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, rot_coeffu2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, rot_coeffv1(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, rot_coeffv2(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, dtheta_dz(model_levels)                                         &
     &, exner_theta_ref(model_levels)                                   &
     &, latitude_theta(rows)                                            &
     &, p_ref(model_levels + 1)                                         &
                                 ! p profile for lbc calc for rho
     &, disti_u(row_length,rows)                                        &
                                   ! i-distance fn on u-grid
     &, distj_v(row_length,n_rows)                                      &
                                   ! j-distance fn on v-grid
     &, thetafn_j(rows)            ! 1-D horizontal theta variation fn


      ! Error reporting
      Character (Len=*),  Parameter :: RoutineName='idl_initial_data_3a'
      Character (Len=256)           :: Cmessage
      Integer                       :: ErrorStatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! ----------------------------------------------------------------------
! Section 0.  Initialise Data fields.
!             Set reference vertical grid
!             Set eta_theta_levels
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('IDL_INITIAL_DATA_3A',zhook_in,zhook_handle)
!     fcrl shorthand for first_constant_rho_level
      fcrl = first_constant_rho_level

      if(me  ==  0)then
        print*,'timestep = ',timestep,' seconds'
        print*,'levels = ',model_levels
        print*,'Equatorial East-West gridlength = ',delta_x,' metres '  &
     &        ,  0.001*delta_x,' kilometres'
        print*,'Equatorial North-South gridlength = ',delta_y,' metres' &
     &        ,  0.001*delta_y,' kilometres'
      Endif   ! me  ==  0

      if ( tprofile_number  /=  tp_dump ) then   ! set idealised data

      if(me  ==  0)then
        print*,' Height_domain (top theta level) = '                    &
     &        ,  height_domain,' metres'
        print*,' first_constant_rho_level set to ', fcrl
        print*,' big_layers = ',big_layers
        print*,' transit_layers = ',transit_layers
      if(L_trivial_trigs .or. .not. L_pressure_balance)then
        temp1 = SQRT(u_in(1) * u_in(1) + v_in(1) * v_in(1))
        print*,' Westerly component = ', u_in(1) ,' m/s '
        print*,' Southerly component = ', v_in(1) ,' m/s '
        if( u_in(1)  >=  v_in(1))then
          print*,' Courant number for constant flow = '                 &
     &          ,    temp1 * timestep / delta_x
        else
          print*,' Courant number for constant flow = '                 &
     &          ,    temp1 * timestep / delta_y
        endif  ! u_in(1) >= v_in(1)
      else ! L_pressure_balance = .true.
        print*,' Westerly component = ', u_in(1) ,' m/s '
        print*,' Courant number at equator = '                          &
     &          ,   u_in(1) * timestep / delta_x
        print*,'Note:  Pressure balance in spherical geometry'
        print*,' Southerly wind component v will be set to 0'
      endif  ! L_trivial_trigs .or. .not. L_pressure_balance
        print*,' Surface temperature = ', theta_surface ,' K'
        print*,' Surface pressure = ', p_surface ,' Pa'
        print*,' Brunt-Vaisala frequency = ', Brunt_Vaisala             &
     &        ,' per second'
      Endif   !(me  ==  0)

      BV_squared_over_g =  Brunt_Vaisala * Brunt_Vaisala / g
      goverRT = g/(R*theta_surface)
      minusgoverRTref = - g / ( R * 250.0)

      recip_delta_lambda = 1./delta_lambda
      recip_delta_phi = 1./delta_phi

! ----------------------------------------------------------------------
! Section 1.  Set temperature profiles and other data
!            Based on user choice of profile number
! ----------------------------------------------------------------------

      If (me == 0) Then
        Write (6,*) ' '
        Write (6,*) ' TEMPERATURE PROFILE '
      End If

! Copy surface value into level 0
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          theta_eh(i,j,0) = theta_surface
        END DO
      END DO


!-----------------------------------------------------------------
! Section 1.1  constant static stability  dtheta_dz
!-------------------------------------------------------------------

      If (tprofile_number  ==  tp_dthetadz) Then

! DEPENDS ON: idl_tprofile1
        Call IDL_tprofile1(                                             &
     &                   row_length, rows, model_levels                 &
     &,                  me, halo_i, halo_j                             &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                  exner_theta_levels                             &
     &,                  theta_ref, exner_ref                           &
     &,                  height_domain, p_surface                       &
     &,                  theta_surface, dtheta_dz1, height_dz1          &
     &,                  L_constant_dz                                  &
     &,                  L_code_test)


!-----------------------------------------------------------------
! Section 1.2  ISOTHERMAL   temperature constant
!        (constant) temperature =  theta_surface
!-------------------------------------------------------------------

      elseif(tprofile_number  ==  tp_isothermal) then

! DEPENDS ON: idl_tprofile2
        Call IDL_tprofile2(                                             &
     &                   row_length, rows, model_levels                 &
     &,                  me, halo_i, halo_j                             &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                  exner_theta_levels                             &
     &,                  theta_ref, exner_ref                           &
     &,                  height_domain                                  &
     &,                  p_surface, theta_surface                       &
     &,                  L_code_test)

!-----------------------------------------------------------------
! Section 1.3  constant Brunt-Vaisala frequency
!-------------------------------------------------------------------

      elseif(tprofile_number  ==  tp_BruntV) then

! DEPENDS ON: idl_tprofile3
        Call IDL_tprofile3(                                             &
     &                   row_length, rows, model_levels                 &
     &,                  me, halo_i, halo_j                             &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                  exner_theta_levels                             &
     &,                  theta_ref, exner_ref                           &
     &,                  height_domain                                  &
     &,                  p_surface, theta_surface, Brunt_Vaisala        &
     &,                  L_code_test)


!-----------------------------------------------------------------
! Section 1.4  constant Brunt-Vaisala frequency
!        ** PLUS isothermal big_layers ***
!------------------------------------------------------------------

      elseif(tprofile_number  ==  tp_BV_isoth) then

! DEPENDS ON: idl_tprofile4
        Call IDL_tprofile4(                                             &
     &                   row_length, rows, model_levels                 &
     &,                  me, halo_i, halo_j                             &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                  exner_theta_levels                             &
     &,                  theta_ref, exner_ref                           &
     &,                  height_domain, big_layers                      &
     &,                  p_surface, theta_surface, Brunt_Vaisala        &
     &,                  L_code_test)

!-----------------------------------------------------------------
! Section 1.5  constant potential temperature (isentropic)
!------------------------------------------------------------------

      elseif(tprofile_number  ==  5) then

! DEPENDS ON: idl_tprofile5
        Call IDL_tprofile5(                                             &
                           row_length, rows, model_levels,              &
                           me, halo_i, halo_j,                          &
                           r_theta_levels, r_rho_levels,                &
                           eta_theta_levels, eta_rho_levels,            &
                           theta_eh(1-halo_i,1-halo_j,1),               &
                           exner_rho_levels_eh,                         &
                           exner_theta_levels,                          &
                           theta_ref, exner_ref,                        &
                           height_domain, p_surface, theta_surface,     &
                           L_code_test)

!-----------------------------------------------------------------
! Section 1.6  Zonally symmetric initial data based on dynamical core
!              set up (Held-Suarez)
!-------------------------------------------------------------------

      Else If (tprofile_number  ==  tp_dyn_core) Then

! Global option only; normal haloes and variables
! DEPENDS ON: idl_tprofile6
        Call IDL_tprofile6(                                             & 
     &                   row_length, rows, model_levels                 &
     &,                  me, off_x, off_y, halo_i, halo_j               &
     &,                  cos_theta_latitude                             &
     &,                  r_theta_levels, r_rho_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &
     &,                  theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                  exner_theta_levels                             &
     &,                  theta_ref, exner_ref                           &
     &,                  height_domain, big_layers                      &
     &,                  p_surface, theta_surface                       &
     &,                  SuHe_pole_equ_deltaT, SuHe_static_stab         &
     &,                  L_SH_Williamson                                &
     &,                  L_code_test)

!-----------------------------------------------------------------
! Section 1.7  Initial data based on dynamical core set up (Held-Suarez)
!              Symmetric on rotated grid (for LAM configurations)
!-------------------------------------------------------------------

      Else If (tprofile_number  ==  tp_dyn_core_lam) Then

        If (model_domain  /=  mt_lam) Then
          ErrorStatus = 1
          WRITE(cmessage,'(A,I8)')'ERROR: tProfile allowed for LAM ' // &
          'domain ONLY; tProfile_number = ',tprofile_number
          CALL ereport(RoutineName, ErrorStatus, cmessage)
        Endif

! DEPENDS ON: idl_tprofile7
        Call IDL_tprofile7(                                             &
     &                      row_length, rows, model_levels              &
     &,                     me, l_datastart, halo_i, halo_j             &
     &,                     delta_phi                                   &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                     exner_theta_levels                          &
     &,                     theta_ref, exner_ref                        &
     &,                     height_domain, big_layers                   &
     &,                     p_surface, theta_surface, r_plane           &
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     L_SH_Williamson                             &
     &,                     L_code_test)
     

!-----------------------------------------------------------------
! Section 1.9  set profile from namelist data
!
!-------------------------------------------------------------------

      Else If (tprofile_number  ==  tp_namelist) Then

! DEPENDS ON: idl_tprofile9
        Call IDL_tprofile9(                                             &
     &                      row_length, rows, model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
                        ,theta_eh(1-halo_i,1-halo_j,1)                  &
     &,                  exner_rho_levels_eh                            &
     &,                     exner_theta_levels                          &
     &,                     theta_ref, exner_ref                        &
     &,                     height_domain, big_layers, surface_type     &
     &,                     zprofile_data, tprofile_data                &
     &,                     p_surface, theta_surface, num_profile_data  &
     &,                     zprofile_orog, idl_interp_option, hf        &
     &,                     L_code_test)

!-----------------------------------------------------------------
! Profile choice not known
!
!-------------------------------------------------------------------

      else
        cmessage = '**  No profile choice made ***'
        ErrorStatus = 1
        CALL ereport(RoutineName, ErrorStatus, cmessage)
      endif           ! end if for setting vertical profile

!*******************************************************************

      If (tprofile_number  ==  tp_dyn_core) Then
!  Only small halo set for dynamical core profiles
!  Fill large haloes for rho calculation later
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(exner_rho_levels_eh, row_length, rows,           &
     &                 model_levels+1,                                  &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(theta_eh, row_length, rows, model_levels+stlev,  &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)
!  Exner_theta_levels has large halo for this routine only
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(Exner_theta_levels, row_length, rows,            &
     &                 model_levels, halo_i, halo_j, fld_type_p,.FALSE.)
      End If !tprofile_number  ==  tp_dyncore

!-----------------------------------------------------------------
! Section 2  Set pressure fields
!            If balanced pressure then will be overwritten later
!-------------------------------------------------------------------

      p_test = 1.0

!    Check for silly pressure at top
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            If (exner_rho_levels_eh(i,j,model_levels) < p_test) Then
              p_test =  exner_rho_levels_eh(i,j,model_levels)
            endif
          end do
        end do

      if(p_test  <   0.0) then
        cmessage = 'Pressure negative at model top, cannot run'
        ErrorStatus = 1
        CALL ereport(RoutineName, ErrorStatus, cmessage)
      endif   !(p_test  <   0.0)

!-----------------------------------------------------------------
! Section 3  Set moisture fields
!-------------------------------------------------------------------

! DEPENDS ON: idl_set_initial_humidity
        Call idl_set_initial_humidity(                                  &
                            height_domain                               &
      ,                     row_length, rows                            &
      ,                     model_levels, wet_model_levels              &
      ,                     me, halo_i, halo_j                          &
      ,                     qprofile_number, L_dry, q1                  &
      ,                     max_num_profile_data, num_profile_data      &
      ,                     zprofile_data, qprofile_data                &
      ,                     zprofile_orog, idl_interp_option, hf        &
      ,                     r_theta_levels, eta_theta_levels            &
      ,                     q(1-halo_i,1-halo_j,1)                      &
      ,                     theta_eh(1-halo_i,1-halo_j,1)               &
      ,                     exner_theta_levels                          &
      ,                     q_ref, theta_ref, exner_ref                 &
      ,                     L_code_test)


        Write (6,*) '   All condensate/hydrometeor variables ',         &
     &              '(qcl,qcf...) will be set to zero'

! Set additional microphysics variables to zero if in use
      If (L_mcr_qcf2)   qcf2(:,:,:)   = 0.0
      If (L_mcr_qrain)  qrain(:,:,:)  = 0.0
      If (L_mcr_qgraup) qgraup(:,:,:) = 0.0

!-------------------------------------------------------------------
!
!          Calculate 3D field of hydrostatic exner pressure
!      from potential temperature (theta) and water vapour (q)
!
!-------------------------------------------------------------------

       ! Exner pressure is calculated in the tprofile routines
       ! assuming a dry atmosphere. Exner pressure is recalculated
       ! here to include the humidity field.
       ! At present, only call if this is not a dry run to
       ! retain bit comparison with dry runs in previous versions

       If (.not. L_dry) Then

! DEPENDS ON: idl_calc_exner
         Call IDL_calc_exner(                                           &
     &                      row_length, rows, model_levels              &
     &,                     wet_model_levels, me, halo_i, halo_j        &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta_eh(1-halo_i,1-halo_j,1)               &
     &,                     q(1-halo_i,1-halo_j,1)                      &
     &,                     exner_rho_levels_eh                         &
     &,                     exner_theta_levels                          &
     &,                     height_domain                               &
     &,                     theta_ref, q_ref, exner_ref                 &
     &,                     p_surface, theta_surface                    &
     &,                     L_code_test)

       End If


!-------------------------------------------------------------------
!
!           Generate potential temperature anomaly (bubble)
!                 with option of making it saturated
!
!-------------------------------------------------------------------

      If (idl_bubble_option(1) > 0) Then

! DEPENDS ON: idl_initialise_bubble
        Call IDL_initialise_bubble(                                     &
                            row_length, rows, halo_i, halo_j            &
      ,                     model_levels, wet_model_levels              &
      ,                     delta_lambda, delta_phi, base_phi           &
      ,                     r_theta_levels                              &
      ,                     p_zero, me, l_datastart                     &
      ,                     global_row_length, global_rows              &
      ,                     idl_max_num_bubbles, idl_bubble_option      &
      ,                     idl_bubble_max, idl_bubble_height           &
      ,                     idl_bubble_xoffset, idl_bubble_yoffset      &
      ,                     idl_bubble_width, idl_bubble_depth          &
      ,                     L_idl_bubble_saturate, L_trivial_trigs      &
      ,                     exner_theta_levels                          &
      ,                     theta_eh(1-halo_i,1-halo_j,1)               &
      ,                     q(1-halo_i,1-halo_j,1)                      &
        )

      End If ! idl_bubble_option

!-------------------------------------------------------------------
!
!           Generate random number perturbations and add to
!      potential temperature (theta) and water vapour (q) fields
!
!-------------------------------------------------------------------

      If (L_perturb_t .or. L_perturb_q) Then

! DEPENDS ON: idl_random_perturb
        Call IDL_random_perturb(                                        &
     &                      row_length, rows, halo_i, halo_j            &
     &,                     model_levels, wet_model_levels              &
     &,                     eta_theta_levels, height_domain             &
     &,                     me, n_proc, all_proc_group                  &
     &,                     global_row_length, global_rows              &
     &,                     g_rows, g_row_length                        &
     &,                     L_perturb_t, perturb_magnitude_t            &
     &,                     L_perturb_q, perturb_magnitude_q            &
     &,                     L_perturb_correlate_tq                      &
     &,                     L_perturb_correlate_vert                    &
     &,                     L_perturb_correlate_time                    &
     &,                     perturb_type, perturb_height                &
     &,                     theta_eh(1-halo_i,1-halo_j,1)               &
     &,                     q(1-halo_i,1-halo_j,1)                      &
     &                      )

      End If ! L_perturb

!----------------------------------------------------------------------
! Section 5  Set u,v wind fields
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Section 5.1  Set u,v wind fields in the vertical
!            (constant horizontally on height levels)
!
!----------------------------------------------------------------------

      If (me == 0) Then
        Write (6,Fmt=*) ' '
        Write (6,Fmt=*) ' HORIZONTAL WIND  '
      End If
      !-----------------------------------------------------------------
      ! Section 5.1.1  u,v wind constant with height
      !-----------------------------------------------------------------
      If (uvprofile_number  ==  uv_vert_const) Then

        If (me == 0) Then
          Write (6,Fmt='(A29,F6.2,A6,F6.2,A4)')                         &
     &     '   Constant with height, u = ',u_in(1),                     &
     &                            ', v = ',v_in(1),' m/s'
        End If

        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_in(1)
            End Do
          End Do
          Do j = 1-halo_j, n_rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              v_adv(i,j,k) = v_in(1)
            End Do
          End Do
          u_ref(k) = u_in(1)
          v_ref(k) = v_in(1)
        End Do ! k= 1, model_levels

      !----------------------------------------------------------------
      ! Section 5.1.2  u,v wind interpolated in height between 4 values
      !----------------------------------------------------------------
      Else If (uvprofile_number  ==  uv_vert_interp) Then

        If (me == 0) Then
          Write (6,Fmt=*)                                               &
     &      '   Linearly interpolating in height between:'
          Write (6,Fmt='(A7,F6.2,A6,F6.2,A17)')                         &
     &      '   u = ',u_in(1),', v = ',v_in(1),' m/s at sea level'
          Write (6,Fmt='(A7,F6.2,A6,F6.2,A15,F8.2,A2)')                 &
     &      '   u = ',u_in(2),', v = ',v_in(2),' m/s at height ',       &
     &          height_u_in(1),' m'
          Write (6,Fmt='(A7,F6.2,A6,F6.2,A15,F8.2,A2)')                 &
     &      '   u = ',u_in(3),', v = ',v_in(3),' m/s at height ',       &
     &          height_u_in(2),' m'
          Write (6,Fmt='(A7,F6.2,A6,F6.2,A15,F8.2,A37)')                &
     &      '   u = ',u_in(4),', v = ',v_in(4),' m/s at height ',       &
     &          height_u_in(3),' m and constant above to top of model'
        End If

        Do k= 1, model_levels

          ! Set u-component
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              z_at_u = r_at_u(i,j,k) - Earth_radius
              If (z_at_u  <   height_u_in(1)) Then
                u_adv(i,j,k) = u_in(1) + (u_in(2) - u_in(1)) *          &
     &          z_at_u / height_u_in(1)
              Else If (z_at_u  <   height_u_in(2)) Then
                u_adv(i,j,k) = u_in(2) + (u_in(3) - u_in(2)) *          &
     &          (z_at_u - height_u_in(1)) /                             &
     &                              ( height_u_in(2) - height_u_in(1))
              Else If (z_at_u  <   height_u_in(3)) Then
                u_adv(i,j,k) = u_in(3) + (u_in(4) - u_in(3)) *          &
     &          (z_at_u - height_u_in(2)) /                             &
     &                             ( height_u_in(3) - height_u_in(2))
              Else
                u_adv(i,j,k) = u_in(4)
              End If
            End Do
          End Do

          ! Set v-component
          Do j = 1-halo_j, n_rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              z_at_v = r_at_v(i,j,k) - Earth_radius
              If (z_at_v  <   height_u_in(1)) Then
                v_adv(i,j,k) = v_in(1) + (v_in(2) - v_in(1)) *          &
     &          z_at_v / height_u_in(1)
              Else If (z_at_v  <   height_u_in(2)) Then
                v_adv(i,j,k) = v_in(2) + (v_in(3) - v_in(2)) *          &
     &          (z_at_v - height_u_in(1)) /                             &
     &                              ( height_u_in(2) - height_u_in(1))
              Else If (z_at_v  <   height_u_in(3)) Then
                v_adv(i,j,k) = v_in(3) + (v_in(4) - v_in(3)) *          &
     &          (z_at_v - height_u_in(2)) /                             &
     &                             ( height_u_in(3) - height_u_in(2))
              Else
                v_adv(i,j,k) = v_in(4)
              End If
            End Do
          End Do

          ! Set reference (no orography) u,v profiles
          temp = eta_rho_levels(k) * height_domain
          v_ref(k) = v_in(1)
          If (temp  <   height_u_in(1)) Then
            u_ref(k) = u_in(1) + (u_in(2) - u_in(1)) *                  &
     &         temp / height_u_in(1)
            v_ref(k) = v_in(1) + (v_in(2) - v_in(1)) *                  &
     &         temp / height_u_in(1)
          Else If (temp  <   height_u_in(2)) Then
            u_ref(k) = u_in(2) + (u_in(3) - u_in(2)) *                  &
     &         (temp - height_u_in(1)) /                                &
     &                             ( height_u_in(2) - height_u_in(1))
            v_ref(k) = v_in(2) + (v_in(3) - v_in(2)) *                  &
     &         (temp - height_u_in(1)) /                                &
     &                             ( height_u_in(2) - height_u_in(1))
          Else If (temp  <   height_u_in(3)) Then
            u_ref(k) = u_in(3) + (u_in(4) - u_in(3)) *                  &
     &         (temp - height_u_in(2)) /                                &
     &                             ( height_u_in(3) - height_u_in(2))
            v_ref(k) = v_in(3) + (v_in(4) - v_in(3)) *                  &
     &         (temp - height_u_in(2)) /                                &
     &                             ( height_u_in(3) - height_u_in(2))
          Else
            u_ref(k) = u_in(4)
            v_ref(k) = v_in(4)
          End If

        End Do ! on k=1,model_levels

      !----------------------------------------------------------------
      ! Section 5.1.3  u,v wind interpolated in hght from namelist data
      !----------------------------------------------------------------
      Else If (uvprofile_number  ==  uv_vert_namelist) Then

        If (me == 0) Then
          Write (6,*) '   Setting horizontal winds from namelist ',     &
     &      'profile (see initial profile print out below)'
        End If
        ! Check to make sure the namelist profile data extends
        ! to the top of the model.
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            If (r_rho_levels(i,j,model_levels) - Earth_radius           &
     &           >   z_uvprofile_data(num_uvprofile_data)) Then
              Write(Cmessage,*)                                         &
     &          'Idealised namelist u,v vertical profile data'          &
     &          //'does not extend to the top of the model.'            &
     &          //'Please modify the namelist data.'
              ErrorStatus = 1

              Call Ereport( RoutineName, ErrorStatus, Cmessage )
            End If
          End Do
        End Do

        ! Set u-component
        Do k = 1, model_levels
          Do k2 = 1, num_uvprofile_data-1
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i

                z_at_u = r_at_u(i,j,k) - Earth_radius

                If (z_at_u > z_uvprofile_data(k2) .and.                 &
     &              z_at_u <= z_uvprofile_data(k2+1)) Then

                  weight = (z_at_u - z_uvprofile_data(k2))              &
     &                 /(z_uvprofile_data(k2+1) - z_uvprofile_data(k2))

                  u_adv(i,j,k) = uprofile_data(k2) + weight*            &
     &                   (uprofile_data(k2+1) - uprofile_data(k2))
                End If

              End Do
            End Do
          End Do
          ! Set reference profile
          u_ref(k) = u_adv(1,1,k)
        End Do

        ! Set v-component
        Do k = 1, model_levels
          Do k2 = 1, num_uvprofile_data-1
            Do j = 1-halo_j, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i

                z_at_v = r_at_v(i,j,k) - Earth_radius
                If (z_at_v > z_uvprofile_data(k2) .and.                 &
     &              z_at_v <= z_uvprofile_data(k2+1)) Then

                  weight = (z_at_v - z_uvprofile_data(k2))              &
     &                 /(z_uvprofile_data(k2+1) - z_uvprofile_data(k2))

                  v_adv(i,j,k) = vprofile_data(k2) + weight*            &
     &                   (vprofile_data(k2+1) - vprofile_data(k2))
                End If

              End Do
            End Do
          End Do
          ! Set reference profile
          v_ref(k) = v_adv(1,1,k)
        End Do

! Set up u_ref, v_ref for damping layers
        Do k = 1, model_levels
          z_at_u = eta_rho_levels(k)*height_domain
          Do k2 = 1, num_profile_data-1
            If (z_at_u  >   zprofile_data(k2) .and.                     &
     &        z_at_u  <=  zprofile_data(k2+1)) Then

              weight = (z_at_u - zprofile_data(k2))                     &
     &                 /(zprofile_data(k2+1) - zprofile_data(k2))
              u_ref(k) = uprofile_data(k2) + weight*                    &
     &                   (uprofile_data(k2+1) - uprofile_data(k2))
              v_ref(k) = vprofile_data(k2) + weight*                    &
     &                   (vprofile_data(k2+1) - vprofile_data(k2))
            End If
          End Do  ! k2 loop
        End Do    ! k loop

      !-----------------------------------------------------------------
      ! Alternative interpolation of initial profile over orography
      !-----------------------------------------------------------------
      !
      !  idl_interp_option = 1: constant on height levels
      !                         (default above, no need to modify)
      !  idl_interp_option = 2: hybrid height everywhere up to a
      !                         specified height "hf" (the same as model
      !                         levels are defined). hs = height_domain
      !  idl_interp_option = 3: as option 2 but only when orography is
      !                         less than input profile orography.
      !                         hs = zprofile_orog
      !
      ! Sets up an eta coord for each level of the initial profile data
      ! and an eta coordinate for each model column, and interpolates
      ! in eta space if the model level height is less than "hf" and
      ! the model orography height is less than "hs"
      !-----------------------------------------------------------------

      If (idl_interp_option == 2 .or. idl_interp_option == 3) Then

        If (idl_interp_option == 2) hs = height_domain
        If (idl_interp_option == 3) hs = zprofile_orog

        ! Set up an eta coord for each level of the initial profile data
        eta_profile(1) = 0.0
        Do k2 = 2, num_profile_data
          eta_profile(k2) =  (z_uvprofile_data(k2) - zprofile_orog)     &
     &                      /(hf - zprofile_orog)
        End Do

        ! Interpolate in eta space and overwrite u where appropriate
        Do k = 1, model_levels
          Do k2 = 1, num_profile_data - 1
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                z_at_u = r_at_u(i,j,k) - Earth_radius
                z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
                eta_model = (z_at_u - z_at_orog)/(hf - z_at_orog)

                If ( (z_at_orog <= hs ) .and. (z_at_u < hf) .and.       &
     &               (eta_model > eta_profile(k2))  .and.               &
     &               (eta_model <= eta_profile(k2+1)) ) Then

                  weight = (eta_model - eta_profile(k2))                &
     &                     /(eta_profile(k2+1) - eta_profile(k2))
                  u_adv(i,j,k) = uprofile_data(k2) + weight*            &
     &                      (uprofile_data(k2+1) - uprofile_data(k2))
                End If
              End Do
            End Do
          End Do
        End Do

        ! Interpolate in eta space and overwrite v where appropriate
        Do k = 1, model_levels
          Do k2 = 1, num_profile_data - 1
            Do j = 1-halo_j, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                z_at_v = r_at_v(i,j,k) - Earth_radius
                z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
                eta_model = (z_at_v - z_at_orog)/(hf - z_at_orog)

                If ( (z_at_orog <= hs ) .and. (z_at_v < hf) .and.       &
     &               (eta_model > eta_profile(k2))  .and.               &
     &               (eta_model <= eta_profile(k2+1)) ) Then

                  weight = (eta_model - eta_profile(k2))                &
     &                     /(eta_profile(k2+1) - eta_profile(k2))
                  v_adv(i,j,k) = vprofile_data(k2) + weight*            &
     &                      (vprofile_data(k2+1) - vprofile_data(k2))
                End If
              End Do
            End Do
          End Do
        End Do

      End If ! on idl_interp_option = 2 or 3

      !----------------------------------------------------------------
      ! Section 5.1.4  u,v wind profile taken from input dump
      !----------------------------------------------------------------
      Else If (uvprofile_number  ==  uv_vert_dump) Then

        ! Do nothing
        If (me == 0) Then
          Write (6,*) '   Using u,v wind field from the dump'
        End If

      End If    ! on uvprofile_number


!----------------------------------------------------------------------
! Section 5.2  Modify u,v wind fields in the horizontal
!             (i.e. on height levels)
!
!----------------------------------------------------------------------

      !----------------------------------------------------------------
      ! Section 5.2.1  u,v wind horizontally constant (on hght levels)
      !----------------------------------------------------------------
      If (uv_horizfn_number  ==  uv_horiz_const ) Then

        ! Horizontally constant wind field
        ! so the wind field set in the previous section is unchanged.

      !----------------------------------------------------------------
      ! Section 5.2.2  u,v wind ramped down towards poles
      !----------------------------------------------------------------
      Else If (uv_horizfn_number  ==  uv_horiz_ramp) Then

        If (me  ==  0) Then
          Write(6,*) 'IDL_Initial_Data_3A: '                            &
             //' u,v wind ramped down towards poles '
          Write(6,*) 'u_ramp_start = ', u_ramp_start
          Write(6,*) 'u_ramp_end = ', u_ramp_end
        End If   ! me  ==  0  

        ! Ramp u down to zero

        !  starting from prescribed latitude (degrees) u_ramp_start
        cos_ramp_start = cos(Pi_over_180 * u_ramp_start)

        !  ending at prescribed latitude  (degrees) u_ramp_end
        cos_ramp_end   = cos(Pi_over_180 * u_ramp_end)

        ! Calculate latitude of each row except N pole 
        ! For v-at-poles n_rows=rows except at N pole
        row_min = MIN( rows, n_rows )
        Do j = 1-halo_j, row_min+halo_j

          If (cos_theta_latitude(1,j) <= cos_ramp_end) Then
 
            Do k = 1, model_levels
              Do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = 0.0
                v_adv(i,j,k) = 0.0
              End Do
            End Do  !  k = 1, model_levels

          Else If (cos_theta_latitude(1,j) <= cos_ramp_start) Then

            weight = (cos_theta_latitude(1,j) - cos_ramp_end) /         &
                              (cos_ramp_start - cos_ramp_end)
            Do k = 1, model_levels
              Do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = weight * u_adv(i,j,k)
                v_adv(i,j,k) = weight * v_adv(i,j,k)
              End Do
            End Do  !  k = 1, model_levels
            
          End If

        End Do ! j = 1-halo_j, row_min+halo_j

        ! North pole will not have been done
        ! North pole halo for v, v_adv only 
        DO j = n_rows, n_rows + halo_j 
          If (cos_theta_latitude(1,j) <= cos_ramp_end) Then
 
            Do k = 1, model_levels
              Do i = 1-halo_i, row_length+halo_i
                v_adv(i,j,k) = 0.0
              End Do
            End Do  !  k = 1, model_levels

          Else If (cos_theta_latitude(1,j) <= cos_ramp_start) Then

            weight = (cos_theta_latitude(1,j) - cos_ramp_end) /         &
                              (cos_ramp_start - cos_ramp_end)
            Do k = 1, model_levels
              Do i = 1-halo_i, row_length+halo_i
                v_adv(i,j,k) = weight * v_adv(i,j,k)
              End Do
            End Do  !  k = 1, model_levels
            
          End If
        END DO ! j = n_rows, n_rows + halo_j 

      !----------------------------------------------------------------
      ! Section 5.2.2  u,v balanced from input dump pressure field
      !----------------------------------------------------------------
      Else If (uv_horizfn_number  ==  uv_horiz_balance) Then

        If (me  ==  0) Then
          Write(6,*) 'IDL_Initial_Data_3A: Wind field derived from '    &
     &       //'geostrophic balance with input pressure field. '        &
     &       //'OPTION NOT YET AVAILABLE.'
!          Write(6,*) 'IDL_Initial_Data_3A: Wind field derived from '
!     &       //'geostrophic balance with input pressure field.'
!          Write(6,*) 'Input wind field will be overwritten'
        End If   !(me  ==  0)

!       Call IDL_geostrophy (not yet active)

      !----------------------------------------------------------------
      ! Section 5.2.3  u,v deformation field (LAM only)
      !----------------------------------------------------------------
      Else If (uv_horizfn_number  ==  uv_horiz_deform) Then

        If (me  ==  0) Then
          Print*,'** Deformation wind field **'
        End If   !(me  ==  0)

        If (model_domain  /=  mt_lam) Then
          Write(Cmessage,*)                                             &
     &          'ERROR: Deformation field only allowed for LAM domain'
          ErrorStatus = 1

          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

! DEPENDS ON: idl_setup_def_front
        Call IDL_setup_def_front(                                       &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     global_row_length, global_rows              &
     &,                     l_datastart                                 &
     &,                     halo_i, halo_j, off_x, off_y                &
     &,                     delta_lambda, delta_phi                     &
     &,                     lambda_p, phi_p, lambda_u, phi_v            &
     &,                     L_regular                                   &
     &,                     base_phi                                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     q, theta_eh, rho_eh, u_adv, v_adv           &
     &,                     exner_rho_levels_eh                         &
     &,                     exner_theta_levels                          &
     &,                     f3_at_u                                     &
     &,                     u_ref, t_horizfn_data                       &
     &,                     height_domain                               &
     &,                     L_trivial_trigs                             &
     &                      )

      End If  ! on uv_horizfn_number

!----------------------------------------------------------------
!   A routine to initialise a limited area domain
!   with a warm core vortex, analogous
!   to a spinning down cyclone or an idealised hurricane
!----------------------------------------------------------------

      IF (L_cyclone) THEN
! DEPENDS ON: idl_cyclone
        Call IDL_CYCLONE(                                               &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     global_row_length, global_rows              &
     &,                     l_datastart                                 &
     &,                     halo_i, halo_j, off_x, off_y                &
     &,                     delta_lambda, delta_phi                     &
     &,                     base_phi                                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     q, theta_eh, u_adv, v_adv                   &
     &,                     exner_rho_levels_eh                         &
     &,                     exner_theta_levels                          &
     &,                     f3_at_u                                     &
     &,                     t_horizfn_data                              &
     &,                     theta_ref                                   &
     &,                     L_trivial_trigs                             &
     &                      )
      END IF

!----------------------------------------------------------------
!  A routine to initialise a limited area domain
!  with a baroclinic jet and tropopause perturbation
!  in order to simulate baroclinic instability.
!----------------------------------------------------------------

      IF (L_baroclinic) THEN
! DEPENDS ON: idl_baroclinic
        Call IDL_BAROCLINIC(                                            &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     global_row_length, global_rows              &
     &,                     l_datastart                                 &
     &,                     halo_i, halo_j, off_x, off_y                &
     &,                     delta_lambda, delta_phi                     &
     &,                     base_phi                                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     q, theta_eh, u_adv, v_adv                   &
     &,                     exner_rho_levels_eh                         &
     &,                     exner_theta_levels                          &
     &,                     f3_at_u                                     &
     &,                     t_horizfn_data                              &
     &,                     theta_ref                                   &
     &,                     L_trivial_trigs                             &
     &                      )
      END IF


!----------------------------------------------------------------------
! Section 5.3  Other u,v, wind options
!
!----------------------------------------------------------------------

      !----------------------------------------------------------------
      ! Section 5.3.1  Set polar points to zero for global
      !----------------------------------------------------------------

      If (L_polar_wind_zero .and. (model_domain  ==  mt_global)) Then

        If (me  ==  0) Then
          Write(6,*) 'IDL_Initial_Data_3A: '                            &
             //' Set polar points to zero for global '
          Write(6,*) 'L_polar_wind_zero = ', L_polar_wind_zero
        End If   ! me  ==  0  

        If (at_extremity(PSouth)) Then
          Do k= 1,model_levels
            Do j = 1-halo_j, 1
              Do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = 0.0
                v_adv(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If        !at_extremity(PSouth)

        If (at_extremity(PNorth)) Then
          Do k= 1,model_levels
            Do j = rows, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = 0.0
              End Do
            End Do
            Do j = n_rows, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                v_adv(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If     !at_extremity(PNorth)

      End If     !  L_polar_wind_zero

      !----------------------------------------------------------------
      ! Section 5.3.2  Rotate winds for LAM equatorial lat-long grid
      !----------------------------------------------------------------
      If (L_rotate_winds) Then

! rotate winds for LAM equatorial latitude-longitude since prescribed
!  winds are relative to standard latitude-longitude grid
! Code below assumes that on a level u = constant1, v = constant2
! Hence we can assume that these are B-grid values on reg. lat-lon grid

        If (model_domain  ==  mt_global) Then

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot rotate winds'       &
     &      //' (L_rotate_winds=.true.) if global model.'
          End If  !(me  ==  0)

        Else If (uv_horizfn_number  ==  uv_horiz_balance) Then

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot rotate winds'       &
     &      //' (L_rotate_winds=.true.) if wind is balanced from'       &
     &      //' input pressure field'                                   &
     &      //' (i.e. uv_horizfn_number=uv_horiz_balance).'
          End If  !(me  ==  0)

        Else If (uv_horizfn_number  ==  uv_horiz_deform) Then

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot rotate winds'       &
     &      //' (L_rotate_winds=.true.) if deformation wind field'      &
     &      //' (uv_horizfn_number=uv_horiz_deform) is chosen.'
          End If  !(me  ==  0)

        Else

          If (me  ==  0) Then
            Print*,'** Input u field is the true westerly wind **'
            Print*,'** Therefore apply rotation to derive  **'
            Print*,'** u and v components on the rotated LAM grid **'
          End If   !(me  ==  0)

          u_field_size = (row_length + 2*halo_i) * (rows + 2*halo_j)
          v_field_size = (row_length + 2*halo_i) * (n_rows + 2*halo_j)

! calculate lat/longitude for u points on equatorial C-grid
! longitude in work1, latitude in work2
          If (L_regular) Then
          Do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            Do i = 1-halo_i, row_length+halo_i
              gi = l_datastart(1) + i - 1
              work1(i,j) = (base_lambda + (gi-1) * delta_lambda)        &
     &                         / Pi_over_180
              work2(i,j) = (base_phi + (gj-.5) * delta_phi)             &
     &                         / Pi_over_180
            End Do
          End Do
          Else
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                work1(i,j) =  lambda_u(i)/ Pi_over_180
                work2(i,j) =  phi_p(i,j)/ Pi_over_180
              End Do
            End Do
          End If   ! L_regular

! true_longitude in work3, true_latitude in work4
! DEPENDS ON: eqtoll
          Call eqtoll(work2, work1                                      &
     &,               work4, work3                                      &
     &,               lat_rot_NP_in, long_rot_NP_in                     &
     &,               u_field_size )

! Calculate rotation coefficients for u-points
! true_longitude in work3, rot_longitude in work4
! DEPENDS ON: w_coeff
          Call w_coeff(rot_coeffu1, rot_coeffu2,                        &
     &                 work3, work1,                                    &
     &                 lat_rot_NP_in, long_rot_NP_in,                   &
     &                  u_field_size )

! calculate lat/longitude for v points on equatorial C-grid
! longitude in work5, latitude in work6
          If (L_regular) Then
          Do j = 1-halo_j, n_rows+halo_j
            gj = l_datastart(2) + j - 1
            Do i = 1-halo_i, row_length+halo_i
              gi = l_datastart(1) + i - 1
              work5(i,j) = (base_lambda + (gi-.5) * delta_lambda)       &
     &                         / Pi_over_180
              work6(i,j) = (base_phi + (gj-1) * delta_phi)              &
     &                         / Pi_over_180
            End Do
          End Do
          Else
            Do j = 1-halo_j, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                work5(i,j) =  lambda_p(i)/ Pi_over_180
                work6(i,j) =  phi_v(i,j)/ Pi_over_180
              End Do
            End Do
          End If   ! L_regular

! true_longitude in work7, true_latitude in work8
! DEPENDS ON: eqtoll
          Call eqtoll(work6, work5                                      &
     &,               work8, work7                                      &
     &,               lat_rot_NP_in, long_rot_NP_in                     &
     &,               v_field_size )

! Calculate rotation coefficients for v-points
! true_longitude in work7, rot_longitude in work5
! DEPENDS ON: w_coeff
          Call w_coeff(rot_coeffv1, rot_coeffv2,                        &
     &               work7, work5,                                      &
     &              lat_rot_NP_in, long_rot_NP_in,                      &
     &               v_field_size )

          DO  k = 1, model_levels

! Copy  u  at v points into work5 array
! Copy  v  at u points into work1 array
! Can just copy since u=constant1 and v=constant2 on a level
            Do j = 1-halo_j, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                work5(i,j) = u_adv(i,j,k)
                work1(i,j) = v_adv(i,j,k)
              End Do
            End Do
! Need v  at last row so just copy from penultimate row
            j = rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              work1(i,j) = v_adv(i,j-1,k)
            End Do

! Find rotated u-component
! DEPENDS ON: w_lltoeq
            Call W_LLTOEQ(rot_coeffu1,rot_coeffu2,                      &
     &         u_adv(1-halo_i,1-halo_j,k), work1,                       &
     &         work3, work4, u_field_size, .TRUE.)
!  u field returned in work3 - work4 is v field - not needed

! Copy into u_adv
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = work3(i,j)
              End Do
            End Do

! Find rotated v-component
! DEPENDS ON: w_lltoeq
            Call W_LLTOEQ(rot_coeffv1,rot_coeffv2,                      &
     &         work5, v_adv(1-halo_i,1-halo_j,k),                       &
     &         work7, work8, v_field_size, .TRUE.)
! Non-halo v field returned in work8 - work7 is u field - not needed

! Copy into v_adv
            Do j = 1-halo_j, n_rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                v_adv(i,j,k) = work8(i,j)
              End Do
            End Do

          End Do  ! k = 1, model_levels

        End If   ! on (uv_horizfn_number)

      End If   ! L_rotate_winds


      !----------------------------------------------------------------
      ! Section 5.3.3  Balance pressure from wind field
      !----------------------------------------------------------------
      ! SCANIA balance pressure field balanced from u-field
      ! NB  v must = 0 for this spherical geometry version
      !----------------------------------------------------------------
      If (L_pressure_balance) Then

        If ( L_rotating) Then

          If (me  ==  0) Then
            Print*,'** Pressure field balanced from wind field **'
            If (model_domain  ==  mt_global) Then
              Print*,'**  u is constant or function of latitude only**'
            Else ! Various LAM configurations
              Print*,'** u is constant or function of latitude on '
              Print*,'   ROTATED LAM grid'
            EndIf     !  model_domain  ==  mt_global
          End If   !(me  ==  0)

! DEPENDS ON: idl_pr_balance
          Call IDL_pr_balance(                                          &
     &                      model_domain, row_length, rows, n_rows      &
     &,                     model_levels                                &
     &,                     delta_x, delta_y, halo_i, halo_j            &
     &,                     me, l_datastart                             &
     &,                     delta_lambda, delta_phi                     &
     &,                     lambda_p, phi_p                             &
     &,                     L_regular                                   &
     &,                     Base_phi, base_lambda                       &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     theta_eh(1-halo_i,1-halo_j,1)               &
     &,                     exner_rho_levels_eh, u_adv, v_adv           &
     &,                     theta_ref, exner_ref                        &
! Profile settings
     &,                     height_domain, theta_surface, Brunt_Vaisala &
     &,                     u_in, v_in, u_ramp_start, u_ramp_end        &
     &,                     ujet_lat, ujet_width                        &
     &,                     f_plane, r_plane                            &
!  Options
     &,                     L_trivial_trigs, L_code_test)

        Else If (uv_horizfn_number  ==  uv_horiz_balance) Then

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot balance pressure field'&
            //' (L_rotate_winds=.true.) if wind is balanced from'       &
            //' the input pressure field'                               &
            //' (i.e. uv_horizfn_number=uv_horiz_balance).'
          End If  !(me  ==  0)

        Else If (uv_horizfn_number  ==  uv_horiz_deform) Then

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot balance pressure field'&
            //' (L_rotate_winds=.true.) if deformation wind field'      &
            //' (uv_horizfn_number=uv_horiz_deform) is chosen.'
          End If  !(me  ==  0)

        Else   !  L_rotating is .FALSE.

          If (me  ==  0) Then
            Write(6,*) 'IDL_Initial_Data_3A: Cannot balance pressure field'&
            //' with winds if rotation is off (L_rotating=.false.)'
          End If  !(me  ==  0)

        End If  !  L_rotating

! Set reference profiles
        DO k = 1, model_levels
          u_ref(k) = u_in(1)
          v_ref(k) = v_in(1)
        END DO

      End If  !  L_pressure_balance

! ----------------------------------------------------------------------
! Section 6. Calculate rho from equation of state and multiply by
!            r*r.
! ----------------------------------------------------------------------
      If (me  ==  0) Then
        Write (6,*) ' '
        Write (6,*) ' AIR DENSITY   '
        Write (6,*) '   Calculating density of air from ',              &
     &                 'Equation of State'
      End If   !(me  ==  0)

      do k = 1, model_levels + 1
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
           If (exner_rho_levels_eh(i,j,k) <= 0.0) Then
        Print*,' pressure negative on processor ',me,'  level ',k
        Print*,' pressure negative at i= ',i,'  j= ',j
        WRITE(6,'(A)')' Is your height_domain too high ?'
        Print*,'Your profile results in an atmosphere of limited extent'
        Print*,' Run is STOPPING at start of '
        Print*,' Section 6 of Subroutine IDL_Initial_data_3A'
           ErrorStatus = 1
           cmessage = "Negative pressure detected"
           CALL ereport(RoutineName, ErrorStatus, cmessage)
           End If ! exner_rho_levels_eh(i,j,k) <= 0.0
          end do
        end do
      end do

      do k = 1, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            ! Store pressure temporarily in rho_eh
            rho_eh(i,j,k) = p_zero * exner_rho_levels_eh(i,j,k)         &
     &                       ** recip_Kappa
          end do
        end do
      end do

      do k = 1, model_levels + 1
        p_ref(k) = p_zero * exner_ref(k) ** recip_Kappa
      end do

! Note: rho_eh contains pressure initially and is then overwritten
!       by the air density
      k = 1
      Do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          rho_eh(i,j,k) = r_rho_levels(i,j,k) * r_rho_levels(i,j,k) *   &
     &                     rho_eh(i,j,k) / (R * theta_eh(i,j,k) *       &
     &                     (1. +  (recip_epsilon -1.) * q(i,j,k))       &
     &                      * exner_rho_levels_eh(i,j,k))
        End Do
      End Do
      rtemp = eta_rho_levels(k) * height_domain + Earth_radius
      rho_ref(k) = rtemp * rtemp * p_ref(k) /                           &
     &                   (R * theta_ref(k) * exner_ref(k))

      Do k = 2, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)) /  &
     &               (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))

              If (k-1 > wet_model_levels) Then
                temp = weight * theta_eh(i,j,k) +                       &
     &                  (1.0 - weight) * theta_eh(i,j,k-1)
              Else If (k > wet_model_levels) Then
                temp = weight * theta_eh(i,j,k) +                       &
     &                  (1.0 - weight) * theta_eh(i,j,k-1) *            &
     &                  (1. +  (recip_epsilon -1.) * q(i,j,k-1))
              Else
                temp = weight * theta_eh(i,j,k) *                       &
     &                  (1. +  (recip_epsilon -1.) * q(i,j,k))          &
     &                  + (1.0 - weight) * theta_eh(i,j,k-1) *          &
     &                  (1. +  (recip_epsilon -1.) * q(i,j,k-1))
              End If  ! k > wet_model_levels

              ! Calculate density of air * r^2
              rho_eh(i,j,k) = r_rho_levels(i,j,k)*r_rho_levels(i,j,k) * &
     &                         rho_eh(i,j,k) /                          &
     &                        (R * temp * exner_rho_levels_eh(i,j,k))

            End Do
          End Do
          weight = (eta_rho_levels(k) - eta_theta_levels(k-1)) /        &
     &               (eta_theta_levels(k) - eta_theta_levels(k-1))
          temp =    weight  * theta_ref(k) +                            &
     &                (1.0 - weight) * theta_ref(k-1)
          rtemp = eta_rho_levels(k) * height_domain + Earth_radius
          rho_ref(k) = rtemp * rtemp * p_ref(k) /                       &
     &                     (R * temp * exner_ref(k))
        End Do   !  k = 2, model_levels


! Initialise density for baroclinic and cyclone expts
! to ensure background surface pressure is
! close to 1000 hPa

      IF (L_baroclinic .OR. L_cyclone) THEN
        DO k = 1, model_levels + 1
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              rho_eh(i,j,k) =  rho_ref(k)
            END DO
          END DO
        END DO
      END IF


! ----------------------------------------------------------------------
      else    !    tprofile_number  ==  tp_dump

        if(me  ==  0)then
          print*,'tprofile_number = ',tprofile_number
          print*,'Vertical profile unchanged - as in input dump '
        Endif   !(me  ==  0)

      endif    ! tprofile_number  /=  tp_dump

!-----------------------------------------------------------------
! Output selected diagnostics
!
!-------------------------------------------------------------------

      If (me  ==  0) Then
        Write (6,*) ' '
        Write (6,*) ' INITIAL PROFILE'
        Write(Unit=6,Fmt=*) ' '
        Write(Unit=6,Fmt=*)                                             &
     &    ' Theta-Lev  Height   Delta_z   Theta  WaterVapour',          &
     &    '  Height   Delta_z'
        Write(Unit=6,Fmt=*)                                             &
     &    '           (for sea level pt)                    ',          &
     &    ' (for highest peak)'
        Write(Unit=6,Fmt=*)                                             &
     &    '              (m)      (m)      (K)     (g/kg)',             &
     &    '       (m)      (m)'
        Do k = model_levels, 1, -1
          temp = (eta_theta_levels(k) -  eta_theta_levels(k-1))         &
     &                                * height_domain
          temp1 = z_orog_print(k) - z_orog_print(k-1)
          Write(Unit=6,                                                 &
     &          Fmt='(4X,I3,4X,F8.2,2X,F7.2,2X,F7.2,2X,F8.5,3X,         &
     &                F8.2,2X,F7.2)')                                   &
     &     k, eta_theta_levels(k)*height_domain , temp                  &
     &,    theta_ref(k), q_ref(k)*1000.,z_orog_print(k),temp1
        End Do
      ! Write out surface theta
        k=0
        Write(Unit=6,                                                   &
     &        Fmt='(4X,I3,4X,F8.2,11X,F7.2,13X,F8.2)')                  &
     &     k, eta_theta_levels(k)*height_domain, theta_surface,         &
     &     z_orog_print(k)

      ! Write out rho-level profile
        Write(Unit=6,Fmt=*) ' '
        Write(Unit=6,Fmt=*)                                             &
     &    '  Rho-Lev   Height   Pressure  U-wind  V-wind  Density'
        Write(Unit=6,Fmt=*)                                             &
     &    '              (m)      (hPa)    (m/s)   (m/s)  (kg/kg)'
        k = model_levels
        temp = (eta_theta_levels(k) - eta_rho_levels(k) +               &
     &            eta_theta_levels(k)) * height_domain
        k = model_levels + 1
        temp1 = p_zero * exner_ref(k) ** recip_kappa
        Write(Unit=6,Fmt='(4X,I3,4X,F8.2,3X,F7.2,2X,A20)')              &
     &   k, temp, temp1/100., '[Nominal rho-height]'
        Do k = model_levels , 1, -1
          temp1 = p_zero * exner_ref(k) ** recip_kappa
          rtemp = eta_rho_levels(k) * height_domain + Earth_radius
          temp = rho_ref(k)/(rtemp * rtemp)
          Write(Unit=6,Fmt= &
              '(4X,I3,4X,F8.2,3X,F7.2,2X,F6.2,2X,F6.2,3X,F6.4)')        &
     &     k, eta_rho_levels(k)*height_domain, temp1/100.               &
     &,    u_ref(k), v_ref(k), temp
        End Do

        if(first_constant_rho_level  >   model_levels)then
          print*,' ****** EXTREME  DANGER  ********'
          print*,' Choose lower level for first_constant_rho_level'
          print*,' since you have set it ABOVE model top !!!'
        endif    !first_constant_rho_level  >   model_levels

        if( 2.0 * z_orog_print(0)   >                                   &
     &   r_rho_levels(1,1,first_constant_rho_level) - Earth_radius)then
          print*,' *** WARNING        WARNING      WARNING  ****'
          print*,' You may be flattening the levels too quickly  '
          print*,' Run will fail if you are using QUADRATIC flattening'
          print*,' Choose higher level for first_constant_rho_level'
          print*,' Current value is level ',first_constant_rho_level
        endif   ! flattening test

        Write (6,*) ' ====================== End of Idealised ',        &
     &              'Settings ======================='
        Write (6,*) ' '
      End If        ! (me  ==  0)

      IF (lhook) CALL dr_hook('IDL_INITIAL_DATA_3A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Initial_data_3A
!  End subroutine IDL_Initial_data_3A

