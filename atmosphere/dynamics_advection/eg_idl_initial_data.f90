! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_idl_initial_data(                                       &
                      model_domain, row_length, rows, n_rows          &
,                     model_levels, wet_model_levels                  &
,                     tr_vars,tr_levels, boundary_layer_levels        &
,                     first_constant_rho_level                        &
,                     cos_theta_latitude, sec_theta_latitude          &
,                     f3_at_u, f3_at_v, timestep                      &
,                     delta_x, delta_y                                &
,                     off_x, off_y, halo_i, halo_j                    &
,                     me, n_proc, at_extremity                        &
,                     l_datastart, all_proc_group                     &
,                     global_row_length, global_rows                  &
,                     delta_lambda, delta_phi                         &
!  VarRes Grid Spacing
,                    lambda_p, phi_p, lambda_u, phi_v                 &
,                    l_regular                                        &
,                     base_phi, base_lambda                           &
,                     lat_rot_np_in, long_rot_np_in                   &
,                     r_theta_levels, r_rho_levels                    &
,                     r_at_u, r_at_v, z_orog_print                    &
,                     eta_theta_levels, eta_rho_levels                &
,                     theta_eh, rho_eh, exner_rho_levels_eh           &
,                     q, qcl, qcf, qcf2, qrain, qgraup                &
,                     u_adv, v_adv, w_adv                             &
,                     g_rows, g_row_length                            &
,                     nproc_x, nproc_y                                &
!  Grid information
,                  height_domain                                      &
,                  big_layers, transit_layers                         &
,                  surface_type, p_surface                            &
! Profile settings
,                  tprofile_number, qprofile_number                   &
,                  uvprofile_number                                   &
,                  brunt_vaisala                                      &
,                  theta_surface, dtheta_dz1, height_dz1              &
,                  u_in, v_in, height_u_in, ujet_lat, ujet_width      &
,                  u_ramp_start, u_ramp_end, f_plane, r_plane         &
,                  q1, max_num_profile_data, num_profile_data         &
,                  zprofile_data, tprofile_data, qprofile_data        &
,                  num_uvprofile_data, z_uvprofile_data               &
,                  uprofile_data, vprofile_data                       &
,                  zprofile_orog, idl_interp_option, hf               &
! Dynamical core settings
,                  suhe_pole_equ_deltat, suhe_static_stab             &
,                  base_frictional_timescale                          &
,                  frictional_timescale, suhe_sigma_cutoff            &
,                  suhe_level_weight, l_sh_williamson                 &
!  Horizontal function parameters
,                  t_horizfn_number, uv_horizfn_number                &
,                  t_horizfn_data                                     &
,                  l_perturb_t, perturb_magnitude_t                   &
,                  l_perturb_q, perturb_magnitude_q                   &
,                  l_perturb_correlate_tq                             &
,                  l_perturb_correlate_vert                           &
,                  l_perturb_correlate_time                           &
,                  perturb_type, perturb_height                       &
!  Profiles for fixed lbcs and sponge zones
,                  u_ref, v_ref, theta_ref, exner_ref, rho_ref        &
!  Options
,                  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup              &
,                  l_constant_dz, l_rotating                          &
,                  l_fixed_lbcs, l_polar_wind_zero                    &
,                  l_wind_balance, l_rotate_winds                     &
,                  l_pressure_balance, l_physics, l_dry, l_sponge     &
,                  l_cyclone,  l_baroclinic                           &
,                  idl_max_num_bubbles, idl_bubble_option             &
,                  idl_bubble_max, idl_bubble_height                  &
,                  idl_bubble_xoffset, idl_bubble_yoffset             &
,                  idl_bubble_width, idl_bubble_depth                 &
,                  l_idl_bubble_saturate                              &
,                  l_cartesian, l_code_test)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atmos_constants_mod, ONLY : r,epsln=>repsilon,cp,p_zero
USE earth_constants_mod, ONLY : earth_radius, g 
USE conversions_mod, ONLY: pi, pi_over_180

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
USE eqtoll_mod, ONLY: eqtoll
USE w_coeff_mod, ONLY: w_coeff
USE w_lltoeq_mod, ONLY: w_lltoeq
IMPLICIT NONE
!
! Description:
!          Sets up initial data for idealised problems.
!          Initial temperature profile at all points is set equal to
!          the equilibrium tmeperature profile at the equator.
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


INTEGER                                                               &
  model_domain                                                        &
                   ! holds integer code for model domain
, row_length                                                          &
                   ! number of points on a row
, rows                                                                &
                   ! number of rows in a theta field
, n_rows                                                              &
                       ! Local number of rows in a v field
, model_levels                                                        &
                   ! number of model levels
, wet_model_levels                                                    &
                   ! number of model levels where moisture
                   ! variables are held
, boundary_layer_levels                                               &
                          ! number of  boundary_layer_levels
, first_constant_rho_level                                            &
, halo_i                                                              &
                       ! Size of halo in i direction.
, halo_j                                                              &
                       ! Size of halo in j direction.
, off_x                                                               &
                       ! Size of small halo in i
, off_y                                                               &
                       ! Size of small halo in j.
, tr_vars                                                             &
, tr_levels

INTEGER, INTENT (IN) :: max_num_profile_data

REAL                                                                  &
  timestep                                                            &
, theta_surface                                                       &
, p_surface                                                           &
, height_domain                                                       &
, zprofile_orog                                                       &
                 ! Height of orog of initial/forcing profile
, hf                                                                  &
, u_in(4)                                                             &
             ! Input values of zonal u
, v_in(4)                                                             &
             ! Input values of southerly wind v
, ujet_lat                                                            &
             ! To specify centre latitude (degrees) of jet core
, ujet_width                                                          &
             ! To specify width (degrees) of jet
, height_u_in(3)                                                      &
                  ! heights specified by the height_u_in variable
, dtheta_dz1(3)                                                       &
                  !Allows different values of dtheta_dz to be set
, height_dz1(2)                                                       &
                  ! at different heights specified by the
                  ! height_dz variable.
, brunt_vaisala                                                       &
, t_horizfn_data(10)                                                  &
                      ! Data values describing horizontal t fn
, q1                                                                  &
                     ! Allows different values of moisture
, zprofile_data(max_num_profile_data)                                 &
                                         ! heights for t,q
, tprofile_data(max_num_profile_data)                                 &
                                         ! theta profile
, qprofile_data(max_num_profile_data)                                 &
                                         ! humidity profile
, z_uvprofile_data(max_num_profile_data)                              &
                                         ! heights for u,v
, uprofile_data(max_num_profile_data)                                 &
                                         ! u-wind profile
, vprofile_data(max_num_profile_data)                                 &
                                         ! v-wind profile
, theta_ref(model_levels)                                             &
                           !theta profile for use in sponge     & lbcs
, exner_ref(model_levels + 1)                                         &
                               ! Exner profile for use in lbcs
, rho_ref(model_levels)                                               &
                          ! rho profile for use in lbcs
, u_ref(model_levels)                                                 &
                        ! u profile for use in lbcs
, v_ref(model_levels)                                                 &
                        ! u_adv profile for use in lbcs
, q_ref(wet_model_levels)                                             &
                            ! q profile for use in lbcs
, delta_x, delta_y                                                    &
                        ! Resolution at equator
, u_ramp_start                                                        &
                  ! ramping starting latitude for u
, u_ramp_end                                                          &
                  ! ramping ending latitude for u
, r_plane                                                             &
                ! reference latitude for row 1 (bottom row)
, f_plane                                                             &
                ! fixed latitude Coriolis term
, idl_bubble_max(idl_max_num_bubbles)                                 &
, idl_bubble_height(idl_max_num_bubbles)                              &
, idl_bubble_xoffset(idl_max_num_bubbles)                             &
, idl_bubble_yoffset(idl_max_num_bubbles)                             &
, idl_bubble_width(idl_max_num_bubbles)                               &
, idl_bubble_depth(idl_max_num_bubbles)

INTEGER                                                               &
  tprofile_number                                                     &
                     ! temperature profile option
, qprofile_number                                                     &
                     ! moisture profile option
, uvprofile_number                                                    &
                     ! u,v wind profile option
, num_profile_data                                                    &
                     ! number of values in z,q,tprofile_data
, num_uvprofile_data                                                  &
                       ! number of values in u,vprofile_data
, t_horizfn_number                                                    &
                     ! horizontal function no. for temperature
, uv_horizfn_number                                                   &
                     ! horizontal function no. for wind field
, idl_interp_option                                                   &
                     ! Profile interpolation option
, big_layers                                                          &
, transit_layers                                                      &
, surface_type                                                        &
                     ! idealised orography type
, idl_max_num_bubbles                                                 &
                      ! Maximum number of bubbles allowed
, idl_bubble_option(idl_max_num_bubbles)  ! Bubble option number

LOGICAL                                                               &
  l_constant_dz                                                       &
                 ! Sets constant dtheta_dz from dtheta_dz1(1)
, l_fixed_lbcs                                                        &
                 ! Set fixed lateral boundary conditions
, l_wind_balance                                                      &
                 ! Geostrophically balance initial winds
, l_pressure_balance                                                  &
                     ! Geostrophically balance pressures
, l_rotate_winds                                                      &
                   ! rotate input u,v (true) to LAM u,v
, l_polar_wind_zero                                                   &
                      ! set u=0 on polar row
, l_rotating                                                          &
                 ! Planet rotation
, l_physics                                                           &
                 ! physics switch
, l_dry                                                               &
                 ! moisture switch
, l_sponge                                                            &
                 ! sponge switch
, l_cartesian                                                         &
                     !  makes grid Cartesian if .true.
, l_code_test                                                         &
                 ! user switch
, l_mcr_qcf2                                                          &
               ! true if using second prognostic cloud ice
, l_mcr_qrain                                                         &
               ! true if using prognostic rain
, l_mcr_qgraup                                                        &
               ! true if using prognostic graupel
, l_cyclone                                                           &
               ! true if cyclone simulation
, l_baroclinic                                                        &
               ! true if baroclinic wave simulation
, l_idl_bubble_saturate(idl_max_num_bubbles) ! Saturate bubble

! Idealised perturbation settings
LOGICAL, INTENT(IN) :: l_perturb_t
LOGICAL, INTENT(IN) :: l_perturb_q
LOGICAL, INTENT(IN) :: l_perturb_correlate_tq
LOGICAL, INTENT(IN) :: l_perturb_correlate_vert
LOGICAL, INTENT(IN) :: l_perturb_correlate_time
REAL,    INTENT(IN) :: perturb_magnitude_t
REAL,    INTENT(IN) :: perturb_magnitude_q
REAL,    INTENT(IN) :: perturb_height(2)
INTEGER, INTENT(IN) :: perturb_type

INTEGER                                                               &
  me                                                                  &
             ! My processor number
, n_proc                                                              &
             ! Total number of processors
, all_proc_group                                                      &
                 ! Group identifier for all processors.
, l_datastart(2)                                                      &
                       ! First gridpoints held by this processor
, global_row_length                                                   &
                       ! global number of points on a row
, global_rows          ! global number of rows

INTEGER                                                               &
  g_rows(0:n_proc-1)                                                  &
, g_row_length(0:n_proc-1)                                            &
, nproc_x, nproc_y                                                    &
, icode,d,info

REAL  rbuf(1-halo_i:global_row_length+halo_i,                         &
           1-halo_j:global_rows+halo_j)

LOGICAL                                                               &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Include physical constants


REAL                                                                  &
     ! horizontal co-ordinate information
  delta_lambda                                                        &
, delta_phi                                                           &
, base_phi                                                            &
, base_lambda                                                         &
, lat_rot_np_in                                                       &
, long_rot_np_in

REAL                                                                  &
     ! VarRes horizontal co-ordinate information
  lambda_p(1-halo_i:row_length+halo_i)                                &
, phi_p(1-halo_j:rows+halo_j)                                         &
, lambda_u(1-halo_i:row_length+halo_i)                                &
, phi_v(1-halo_j:n_rows+halo_j)

LOGICAL, INTENT(IN) :: l_regular

REAL                                                                  &
     ! vertical co-ordinate information
  r_theta_levels(1-halo_i:row_length+halo_i,                          &
                 1-halo_j:rows+halo_j,0:model_levels)                 &
, r_rho_levels(1-halo_i:row_length+halo_i,                            &
               1-halo_j:rows+halo_j, model_levels)                    &
, r_at_u(1-halo_i:row_length+halo_i,                                  &
               1-halo_j:rows+halo_j, model_levels)                    &
, r_at_v(1-halo_i:row_length+halo_i,                                  &
               1-halo_j:n_rows+halo_j, model_levels)                  &
, eta_theta_levels(0:model_levels)                                    &
, eta_rho_levels(model_levels)                                        &
, cos_theta_latitude(1-off_x:row_length+off_x,                        &
                      1-off_y:rows+off_y)                             &
, sec_theta_latitude(1-off_x:row_length+off_x,                        &
                      1-off_y:rows+off_y)                             &
, f3_at_u (1-off_x:row_length+off_x,                                  &
                      1-off_y:rows+off_y)                             &
, f3_at_v (1-off_x:row_length+off_x,                                  &
                      1-off_y:n_rows+off_y)                           &
, z_orog_print(0:model_levels)

REAL, INTENT (INOUT) ::                                               &
  u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        model_levels)                                                 &
, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,           &
        model_levels)                                                 &
, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        0:model_levels)                                               &
, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
    0:model_levels)                                                   &
, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
      0:model_levels)                                                 &
, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
      0:model_levels)                                                 &
, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
      0:model_levels)                                                 &
, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
      0:model_levels)                                                 &
, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
      0:model_levels)

! Work arrays with extended halos (_eh)
! Needed so that external halo values in LAMS can be set correctly
! for the lateral boundary arrays.
REAL, INTENT (INOUT) ::                                               &
  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
           model_levels)                                              &
, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
           model_levels)                                              &
, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                     &
           1-halo_j:rows+halo_j, model_levels+1)

! Suarez Held variables
REAL                                                                  &
  suhe_pole_equ_deltat                                                &
, suhe_static_stab                                                    &
, suhe_sigma_cutoff                                                   &
, base_frictional_timescale                                           &
, frictional_timescale(model_levels)                                  &
, suhe_level_weight(model_levels)

LOGICAL                                                               &
  l_sh_williamson

! local variables
INTEGER                                                               &
  i                                                                   &
, j                                                                   &
, k                                                                   &
, k2                                                                  &
, fcrl                                                                &
, gi, gj                                                              &
, j_start, j_end                                                      &
, u_field_size                                                        &
, v_field_size                                                        &
, length

REAL                                                                  &
  recip_epsln                                                         &
, recip_delta_lambda                                                  &
, recip_delta_phi                                                     &
, bv_squared_over_g                                                   &
, foverrt                                                             &
, goverrt                                                             &
, govercpt                                                            &
, delta_z                                                             &
, x,y                                                                 &
, temp                                                                &
, temp1                                                               &
, temp2                                                               &
, rtemp                                                               &
, weight                                                              &
, p_test                                                              &
, cos_ramp_start                                                      &
, cos_ramp_end                                                        &
, z_at_theta                                                          &
, z_at_u                                                              &
          ! height above mean sea level at u points
, z_at_v                                                              &
          ! height above mean sea level at v points
, theta_mid, theta_mid2                                               &
, sigma_ref(model_levels)                                             &
, sigma_to_kappa(model_levels)

! Williamson variables1. 1. parameters.

REAL                                                                  &
  p_d                                                                 &
, p_pl                                                                &
, lapse_rate_d                                                        &
, lapse_rate_i                                                        &
, delta_phi_0                                                         &
, a_will                                                              &
, phi_0                                                               &
, t_0

! 2. derived variables
REAL                                                                  &
  p_i                                                                 &
, p_eq                                                                &
, power_d                                                             &
, power_i                                                             &
, latitude                                                            &
, p_lim                                                               &
, minusgoverrtref

 REAL                                                                 &
  eta_model, z_at_orog, hs                                            &
, eta_profile(num_profile_data)

REAL                                                                  &
  exner_theta_levels(1-halo_i:row_length+halo_i                       &
,                    1-halo_j:rows+halo_j, model_levels)              &
, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)             &
, work2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)             &
, work3(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)             &
, work4(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)             &
, work5(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)           &
, work6(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)           &
, work7(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)           &
, work8(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)           &
, rot_coeffu1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)       &
, rot_coeffu2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)       &
, rot_coeffv1(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)     &
, rot_coeffv2(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)     &
, dtheta_dz(model_levels)                                             &
, exner_theta_ref(model_levels)                                       &
, latitude_theta(rows)                                                &
, p_ref(model_levels + 1)                                             &
                           ! p profile for lbc calc for rho
, disti_u(row_length,rows)                                            &
                             ! i-distance fn on u-grid
, distj_v(row_length,n_rows)                                          &
                             ! j-distance fn on v-grid
, thetafn_j(rows)            ! 1-D horizontal theta variation fn


! Error reporting
CHARACTER (LEN=*),  PARAMETER :: routinename='idl_initial_data'
CHARACTER (LEN=256)           :: cmessage
INTEGER                       :: errorstatus


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_INITIAL_DATA',zhook_in,zhook_handle)

! ---------------------------------------------------------------------
! Section 0.  Initialise Data fields.
!             Set reference vertical grid
!             Set eta_theta_levels
! ---------------------------------------------------------------------


!     fcrl shorthand for first_constant_rho_level

fcrl = first_constant_rho_level

IF(me  ==  0)THEN
  WRITE(6,FMT='(A,E16.8,A)')'timestep = ',timestep,' seconds'
  WRITE(6,FMT='(A,I4)')'levels = ',model_levels
  WRITE(6,FMT='(A,E16.8,A,E16.8,A)')'equatorial east-west gridlength = '  &
             ,delta_x,' metres '      &
        ,  0.001*delta_x,' kilometres'
  WRITE(6,FMT='(A,E16.8,A,E16.8,A)')'equatorial north-south gridlength = '&
             ,delta_y,' metres'     &
        ,  0.001*delta_y,' kilometres'
END IF   ! me  ==  0

IF ( tprofile_number  /=  tp_dump ) THEN   ! set idealised data

  IF(me  ==  0)THEN
    WRITE(6,FMT='(A,E16.8,A)')' Height_domain (top theta level) = '     &
          ,  height_domain,' metres'
    WRITE(6,FMT='(A,I4)')' first_constant_rho_level set to ', fcrl
    WRITE(6,FMT='(A,I4)')' big_layers = ',big_layers
    WRITE(6,FMT='(A,I4)')' transit_layers = ',transit_layers

    IF(l_cartesian .OR. .NOT. l_pressure_balance)THEN
      temp1 = SQRT(u_in(1) * u_in(1) + v_in(1) * v_in(1))
      WRITE(6,FMT='(A,E16.8,A)')' westerly component = ', u_in(1) ,' m/s '
      WRITE(6,FMT='(A,E16.8,A)')' southerly component = ', v_in(1) ,' m/s '

      IF( u_in(1)  >=  v_in(1))THEN
        WRITE(6,FMT='(A,E16.8)')' Courant number for constant flow = '  &
            ,    temp1 * timestep / delta_x
      ELSE
        WRITE(6,FMT='(A,E16.8)')' Courant number for constant flow = '      &
            ,    temp1 * timestep / delta_y
      END IF  ! u_in(1) >= v_in(1)

    ELSE ! L_pressure_balance = .true.
      WRITE(6,FMT='(A,E16.8,A)')' westerly component = ', u_in(1) ,' m/s '
      WRITE(6,FMT='(A,E16.8)')' Courant number at equator = '               &
            ,   u_in(1) * timestep / delta_x
      WRITE(6,FMT='(A)')'Note:  Pressure balance in spherical geometry'
      WRITE(6,FMT='(A)')' Southerly wind component v will be set to 0'

    END IF  ! L_Cartesian .or. .not. L_pressure_balance

    WRITE(6,FMT='(A,E16.8,A)')' surface temperature = ', theta_surface ,' K'
    WRITE(6,FMT='(A,E16.8,A)')' surface pressure = ', p_surface ,' Pa'
    WRITE(6,FMT='(A,E16.8,A)')' Brunt-Vaisala frequency = ', brunt_vaisala  &
        ,' per second'

  END IF   !(me  ==  0)

  bv_squared_over_g =  brunt_vaisala * brunt_vaisala / g
  goverrt = g/(r*theta_surface)
  minusgoverrtref = - g / ( r * 250.0)

  recip_delta_lambda = 1./delta_lambda
  recip_delta_phi = 1./delta_phi
  recip_epsln = 1./ epsln

!----------------------------------------------------------------------
! Section 5  Set u,v wind fields
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Section 5.1  Set u,v wind fields in the vertical
!            (constant horizontally on height levels)

!----------------------------------------------------------------------

  IF (me == 0) THEN
    WRITE (6,FMT=*) ' '
    WRITE (6,FMT=*) ' HORIZONTAL WIND  '
  END IF

!----------------------------------------------------------------------
! Section 5.1.1  u,v wind constant with height
!----------------------------------------------------------------------

  IF (uvprofile_number  ==  uv_vert_const) THEN

    IF (me == 0) THEN
      WRITE (6,FMT='(A29,F6.2,A6,F6.2,A4)')                             &
            '   Constant with height, u = ',u_in(1),                    &
                                   ', v = ',v_in(1),' m/s'
    END IF

    DO k = 1, model_levels
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          u_adv(i,j,k) = u_in(1)
        END DO
      END DO
      DO j = 1-halo_j, n_rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          v_adv(i,j,k) = v_in(1)
        END DO
      END DO
      u_ref(k) = u_in(1)
      v_ref(k) = v_in(1)
    END DO ! k= 1, model_levels

!----------------------------------------------------------------------
! Section 5.1.2  u,v wind interpolated in height between 4 values
!----------------------------------------------------------------------

  ELSE IF (uvprofile_number  ==  uv_vert_interp) THEN

    IF (me == 0) THEN
      WRITE (6,FMT=*)                                                   &
        '   Linearly interpolating in height between:'
      WRITE (6,FMT='(A7,F6.2,A6,F6.2,A17)')                             &
        '   u = ',u_in(1),', v = ',v_in(1),' m/s at sea level'
      WRITE (6,FMT='(A7,F6.2,A6,F6.2,A15,F8.2,A2)')                     &
        '   u = ',u_in(2),', v = ',v_in(2),' m/s at height ',           &
            height_u_in(1),' m'
      WRITE (6,FMT='(A7,F6.2,A6,F6.2,A15,F8.2,A2)')                     &
        '   u = ',u_in(3),', v = ',v_in(3),' m/s at height ',           &
            height_u_in(2),' m'
      WRITE (6,FMT='(A7,F6.2,A6,F6.2,A15,F8.2,A37)')                    &
        '   u = ',u_in(4),', v = ',v_in(4),' m/s at height ',           &
          height_u_in(3),' m and constant above to top of model'
    END IF

    DO k = 1, model_levels

! Set u-component
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          z_at_u = r_at_u(i,j,k) - earth_radius
          IF (z_at_u  <   height_u_in(1)) THEN
            u_adv(i,j,k) = u_in(1) + (u_in(2) - u_in(1)) *              &
            z_at_u / height_u_in(1)
          ELSE IF (z_at_u  <   height_u_in(2)) THEN
            u_adv(i,j,k) = u_in(2) + (u_in(3) - u_in(2)) *              &
            (z_at_u - height_u_in(1)) /                                 &
                              ( height_u_in(2) - height_u_in(1))
          ELSE IF (z_at_u  <   height_u_in(3)) THEN
            u_adv(i,j,k) = u_in(3) + (u_in(4) - u_in(3)) *              &
            (z_at_u - height_u_in(2)) /                                 &
                             ( height_u_in(3) - height_u_in(2))
          ELSE
            u_adv(i,j,k) = u_in(4)
          END IF
        END DO
      END DO

! Set v-component
      DO j = 1-halo_j, n_rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          z_at_v = r_at_v(i,j,k) - earth_radius
          IF (z_at_v  <   height_u_in(1)) THEN
            v_adv(i,j,k) = v_in(1) + (v_in(2) - v_in(1)) *              &
            z_at_v / height_u_in(1)
          ELSE IF (z_at_v  <   height_u_in(2)) THEN
            v_adv(i,j,k) = v_in(2) + (v_in(3) - v_in(2)) *              &
                           (z_at_v - height_u_in(1)) /                  &
                           ( height_u_in(2) - height_u_in(1))
          ELSE IF (z_at_v  <   height_u_in(3)) THEN
            v_adv(i,j,k) = v_in(3) + (v_in(4) - v_in(3)) *              &
                           (z_at_v - height_u_in(2)) /                  &
                           ( height_u_in(3) - height_u_in(2))
          ELSE
            v_adv(i,j,k) = v_in(4)
          END IF
        END DO
      END DO

! Set reference (no orography) u,v profiles
      temp = eta_rho_levels(k) * height_domain
      v_ref(k) = v_in(1)
      IF (temp  <   height_u_in(1)) THEN
        u_ref(k) = u_in(1) + (u_in(2) - u_in(1)) *                      &
                   temp / height_u_in(1)
        v_ref(k) = v_in(1) + (v_in(2) - v_in(1)) *                      &
                   temp / height_u_in(1)
      ELSE IF (temp  <   height_u_in(2)) THEN
        u_ref(k) = u_in(2) + (u_in(3) - u_in(2)) *                      &
                  (temp - height_u_in(1)) /                             &
                  ( height_u_in(2) - height_u_in(1))
        v_ref(k) = v_in(2) + (v_in(3) - v_in(2)) *                      &
                   (temp - height_u_in(1)) /                            &
                   ( height_u_in(2) - height_u_in(1))
      ELSE IF (temp  <   height_u_in(3)) THEN
        u_ref(k) = u_in(3) + (u_in(4) - u_in(3)) *                      &
                   (temp - height_u_in(2)) /                            &
                   ( height_u_in(3) - height_u_in(2))
        v_ref(k) = v_in(3) + (v_in(4) - v_in(3)) *                      &
                  (temp - height_u_in(2)) /                             &
                  ( height_u_in(3) - height_u_in(2))
      ELSE
        u_ref(k) = u_in(4)
        v_ref(k) = v_in(4)
      END IF

    END DO ! on k=1,model_levels

!----------------------------------------------------------------
! Section 5.1.3  u,v wind interpolated in hght from namelist data
!----------------------------------------------------------------

  ELSE IF (uvprofile_number  ==  uv_vert_namelist) THEN

    IF (me == 0) THEN
      WRITE (6,*) '   Setting horizontal winds from namelist ',         &
        'profile (see initial profile print out below)'
    END IF

! Check to make sure the namelist profile data extends
! to the top of the model.
    DO j = 1-halo_j, rows+halo_j
      DO i = 1-halo_i, row_length+halo_i
        IF (r_rho_levels(i,j,model_levels) - earth_radius               &
           >   z_uvprofile_data(num_uvprofile_data)) THEN
          WRITE(cmessage,*)                                             &
              'Idealised namelist u,v vertical profile data'            &
              //'does not extend to the top of the model.'              &
              //'Please modify the namelist data.'
          errorstatus = 1

          CALL ereport( routinename, errorstatus, cmessage )
        END IF
      END DO
    END DO

! Set u-component
    DO k = 1, model_levels
      DO k2 = 1, num_uvprofile_data-1
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i

            z_at_u = r_at_u(i,j,k) - earth_radius

            IF (z_at_u > z_uvprofile_data(k2) .AND.                     &
                z_at_u <= z_uvprofile_data(k2+1)) THEN

              weight = (z_at_u - z_uvprofile_data(k2))                  &
                      /(z_uvprofile_data(k2+1) - z_uvprofile_data(k2))

              u_adv(i,j,k) = uprofile_data(k2) + weight*                &
                            (uprofile_data(k2+1) - uprofile_data(k2))
            END IF

          END DO
        END DO
      END DO

! Set reference profile
      u_ref(k) = u_adv(1,1,k)
    END DO

! Set v-component
    DO k = 1, model_levels
      DO k2 = 1, num_uvprofile_data-1
        DO j = 1-halo_j, n_rows+halo_j
          DO i = 1-halo_i, row_length+halo_i

            z_at_v = r_at_v(i,j,k) - earth_radius
            IF (z_at_v > z_uvprofile_data(k2) .AND.                     &
                z_at_v <= z_uvprofile_data(k2+1)) THEN

              weight = (z_at_v - z_uvprofile_data(k2))                  &
                      /(z_uvprofile_data(k2+1) - z_uvprofile_data(k2))

              v_adv(i,j,k) = vprofile_data(k2) + weight*                &
                            (vprofile_data(k2+1) - vprofile_data(k2))
            END IF

          END DO
        END DO
      END DO

! Set reference profile
      v_ref(k) = v_adv(1,1,k)
    END DO

!-----------------------------------------------------------------
! Alternative interpolation of initial profile over orography
!-----------------------------------------------------------------

!  idl_interp_option = 1: constant on height levels
!                         (default above, no need to modify)
!  idl_interp_option = 2: hybrid height everywhere up to a
!                         specified height "hf" (the same as model
!                         levels are defined). hs = height_domain
!  idl_interp_option = 3: as option 2 but only when orography is
!                         less than input profile orography.
!                         hs = zprofile_orog

! Sets up an eta coord for each level of the initial profile data
! and an eta coordinate for each model column, and interpolates
! in eta space if the model level height is less than "hf" and
! the model orography height is less than "hs"
!-----------------------------------------------------------------

    IF (idl_interp_option == 2 .OR. idl_interp_option == 3) THEN

      IF (idl_interp_option == 2) hs = height_domain
      IF (idl_interp_option == 3) hs = zprofile_orog

! Set up an eta coord for each level of the initial profile data
      eta_profile(1) = 0.0
      DO k2 = 2, num_profile_data
        eta_profile(k2) = (z_uvprofile_data(k2) - zprofile_orog)       &
                          /(hf - zprofile_orog)
      END DO

! Interpolate in eta space and overwrite u where appropriate
      DO k = 1, model_levels
        DO k2 = 1, num_profile_data - 1
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              z_at_u = r_at_u(i,j,k) - earth_radius
              z_at_orog  = r_theta_levels(i,j,0) - earth_radius
              eta_model = (z_at_u - z_at_orog)/(hf - z_at_orog)

              IF ( (z_at_orog <= hs ) .AND. (z_at_u < hf) .AND.        &
                   (eta_model > eta_profile(k2))  .AND.                &
                   (eta_model <= eta_profile(k2+1)) ) THEN

                weight = (eta_model - eta_profile(k2))                 &
                        /(eta_profile(k2+1) - eta_profile(k2))
                u_adv(i,j,k) = uprofile_data(k2) + weight*             &
                              (uprofile_data(k2+1) - uprofile_data(k2))
              END IF
            END DO
          END DO
        END DO
      END DO

! Interpolate in eta space and overwrite v where appropriate

      DO k = 1, model_levels
        DO k2 = 1, num_profile_data - 1
          DO j = 1-halo_j, n_rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              z_at_v = r_at_v(i,j,k) - earth_radius
              z_at_orog  = r_theta_levels(i,j,0) - earth_radius
              eta_model = (z_at_v - z_at_orog)/(hf - z_at_orog)

              IF ( (z_at_orog <= hs ) .AND. (z_at_v < hf) .AND.        &
                   (eta_model > eta_profile(k2))  .AND.                &
                   (eta_model <= eta_profile(k2+1)) ) THEN

                weight = (eta_model - eta_profile(k2))                 &
                        /(eta_profile(k2+1) - eta_profile(k2))
                v_adv(i,j,k) = vprofile_data(k2) + weight*             &
                              (vprofile_data(k2+1) - vprofile_data(k2))
              END IF
            END DO
          END DO
        END DO
      END DO

    END IF ! on idl_interp_option = 2 or 3

!----------------------------------------------------------------
! Section 5.1.4  u,v wind profile taken from input dump
!----------------------------------------------------------------

  ELSE IF (uvprofile_number  ==  uv_vert_dump) THEN

! Do nothing

    IF (me == 0) THEN
      WRITE (6,*) '   Using u,v wind field from the dump'
    END IF

  END IF    ! on uvprofile_number

!     maybe it is generally not necessary since any other idealised
!     case will most likely initialise it. Assuming that if we do
!     not read uv from the dump we also do not want w from the dump:

  IF (uvprofile_number  /=  uv_vert_dump) THEN
    IF (me == 0) WRITE(0,*) 'initialising w_adv=0'
    w_adv = 0.
  END IF

!----------------------------------------------------------------------
! Section 5.2  Modify u,v wind fields in the horizontal
!             (i.e. on height levels)

!----------------------------------------------------------------------

!----------------------------------------------------------------
! Section 5.2.1  u,v wind horizontally constant (on hght levels)
!----------------------------------------------------------------

  IF (uv_horizfn_number  ==  uv_horiz_const ) THEN

! Horizontally constant wind field
! so the wind field set in the previous section is unchanged.

!----------------------------------------------------------------
! Section 5.2.2  u,v wind ramped down towards poles
!----------------------------------------------------------------

  ELSE IF (uv_horizfn_number  ==  uv_horiz_ramp) THEN

! Ramp u down to zero
! Global domain only so only small halo region done
!  starting from prescribed latitude (degrees) u_ramp_start

    cos_ramp_start = COS(pi_over_180 * u_ramp_start)

!  ending at prescribed latitude  (degrees) u_ramp_end
    cos_ramp_end   = COS(pi_over_180 * u_ramp_end)

! Calculate latitude of each row

    DO j = 1-off_y,n_rows+off_y
      DO i = 1-off_x,row_length+off_x
        IF (cos_theta_latitude(i,j)  <=  cos_ramp_end) THEN
          DO k = 1, model_levels
            u_adv(i,j,k) = 0.0
            v_adv(i,j,k) = 0.0
          END DO
        ELSE IF (cos_theta_latitude(i,j)  <=  cos_ramp_start) THEN
          weight = (cos_theta_latitude(i,j) - cos_ramp_end) /           &
                   (cos_ramp_start - cos_ramp_end)
          DO k = 1, model_levels
            u_adv(i,j,k) = weight * u_adv(i,j,k)
            v_adv(i,j,k) = weight * v_adv(i,j,k)
          END DO
        END IF
      END DO
    END DO

! North pole halo for u, u_adv only (because rows=n_rows+1)

    j = rows+1
    DO i = 1-off_x,row_length+off_x
      IF (cos_theta_latitude(i,j)  <=  cos_ramp_end) THEN
        DO k = 1, model_levels
          u_adv(i,j,k) = 0.0
        END DO
      ELSE IF (cos_theta_latitude(i,j)  <=  cos_ramp_start)THEN
        weight = (cos_theta_latitude(i,j) - cos_ramp_end) /             &
                 (cos_ramp_start - cos_ramp_end)
        DO k = 1, model_levels
          u_adv(i,j,k) = weight * u_adv(i,j,k)
        END DO
      END IF
    END DO

!----------------------------------------------------------------
! Section 5.2.2  u,v balanced from input dump pressure field
!----------------------------------------------------------------

  ELSE IF (uv_horizfn_number  ==  uv_horiz_balance) THEN

    IF (me  ==  0) THEN
      WRITE(6,*) 'IDL_Initial_Data: Wind field derived from '           &
         //'geostrophic balance with input pressure field. '            &
         //'OPTION NOT YET AVAILABLE.'
    END IF   !(me  ==  0)

!       Call IDL_geostrophy (not yet active)

!----------------------------------------------------------------
! Section 5.2.3  u,v deformation field (LAM only)
!----------------------------------------------------------------

  ELSE IF (uv_horizfn_number  ==  uv_horiz_deform) THEN

    IF (me  ==  0) THEN
      WRITE(6,FMT='(A)')'** Deformation wind field **'
    END IF   !(me  ==  0)

    IF (model_domain  /=  mt_lam) THEN
      WRITE(cmessage,*)                                               &
          'ERROR: Deformation field only allowed for LAM domain'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
    END IF

    WRITE(cmessage,*)                                                 &
          'ERROR: Deformation field not eg-comp'
    errorstatus = 1
    CALL ereport( routinename, errorstatus, cmessage )

  END IF  ! on uv_horizfn_number

!----------------------------------------------------------------
!   A routine to initialise a limited area domain
!   with a warm core vortex, analogous
!   to a spinning down cyclone or an idealised hurricane
!----------------------------------------------------------------

  IF (l_cyclone) THEN
    WRITE(cmessage,*)                                                  &
         'ERROR: cyclone field not eg-comp'
    errorstatus = 1
    CALL ereport( routinename, errorstatus, cmessage )
  END IF

!----------------------------------------------------------------
!  A routine to initialise a limited area domain
!  with a baroclinic jet and tropopause perturbation
!  in order to simulate baroclinic instability.
!----------------------------------------------------------------

  IF (l_baroclinic) THEN
! DEPENDS ON: idl_baroclinic
    WRITE(cmessage,*)                                                    &
          'ERROR: barloclinic       not eg-comp'
    errorstatus = 1
    CALL ereport( routinename, errorstatus, cmessage )
  END IF

!----------------------------------------------------------------------
! Section 5.3  Other u,v, wind options
!----------------------------------------------------------------------

!----------------------------------------------------------------
! Section 5.3.1  Set polar points to zero for global
!----------------------------------------------------------------

  IF (l_polar_wind_zero .AND. (model_domain  ==  mt_global)) THEN

    IF (at_extremity(psouth)) THEN
      DO k= 1,model_levels
        DO j = 1-halo_j, 1
          DO i = 1-halo_i, row_length+halo_i
            u_adv(i,j,k) = 0.0
            v_adv(i,j,k) = 0.0
          END DO
        END DO
      END DO
    END IF        !at_extremity(PSouth)

    IF (at_extremity(pnorth)) THEN
      DO k= 1,model_levels
        DO j = rows, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            u_adv(i,j,k) = 0.0
          END DO
        END DO
        DO j = n_rows, n_rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            v_adv(i,j,k) = 0.0
          END DO
        END DO
      END DO
    END IF     !at_extremity(PNorth)

  END IF     !  L_polar_wind_zero

!----------------------------------------------------------------
! Section 5.3.2  Rotate winds for LAM equatorial lat-long grid
!----------------------------------------------------------------

  IF (l_rotate_winds) THEN

! rotate winds for LAM equatorial latitude-longitude since prescribed
!  winds are relative to standard latitude-longitude grid
! Code below assumes that on a level u = constant1, v = constant2
! Hence we can assume that these are B-grid values on reg. lat-lon grid

    IF (model_domain  ==  mt_global) THEN

      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot rotate winds'              &
          //' (L_rotate_winds=.true.) if global model.'
      END IF  !(me  ==  0)

    ELSE IF (uv_horizfn_number  ==  uv_horiz_balance) THEN

      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot rotate winds'              &
          //' (L_rotate_winds=.true.) if wind is balanced from'         &
          //' input pressure field'                                     &
          //' (i.e. uv_horizfn_number=uv_horiz_balance).'
      END IF  !(me  ==  0)

    ELSE IF (uv_horizfn_number  ==  uv_horiz_deform) THEN

      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot rotate winds'              &
        //' (L_rotate_winds=.true.) if deformation wind field'          &
        //' (uv_horizfn_number=uv_horiz_deform) is chosen.'
      END IF  !(me  ==  0)

    ELSE

      IF (me  ==  0) THEN
        WRITE(6,FMT='(A)')'** Input u field is the true westerly wind **'
        WRITE(6,FMT='(A)')'** Therefore apply rotation to derive  **'
        WRITE(6,FMT='(A)')'** u and v components on the rotated LAM grid **'
      END IF   !(me  ==  0)

      u_field_size = (row_length + 2*halo_i) * (rows + 2*halo_j)
      v_field_size = (row_length + 2*halo_i) * (n_rows + 2*halo_j)

! calculate lat/longitude for u points on equatorial C-grid
! longitude in work1, latitude in work2

      IF (l_regular) THEN
        DO j = 1-halo_j, rows+halo_j
          gj = l_datastart(2) + j - 1
          DO i = 1-halo_i, row_length+halo_i
            gi = l_datastart(1) + i - 1
            work1(i,j) = (base_lambda + (gi-.5) * delta_lambda)           &
                         /pi_over_180
            work2(i,j) = (base_phi + (gj-1) * delta_phi)                  &
                         /pi_over_180
          END DO
        END DO
      ELSE
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            work1(i,j) =  lambda_u (i)/ pi_over_180
            work2(i,j) =  phi_p (j)/ pi_over_180
          END DO
        END DO
      END IF   ! L_regular

! true_longitude in work3, true_latitude in work4
      CALL eqtoll(work2, work1,                                         &
                  work4, work3,                                         &
                  lat_rot_np_in, long_rot_np_in,                        &
                  u_field_size )

! Calculate rotation coefficients for u-points
! true_longitude in work3, rot_longitude in work4
      CALL w_coeff(rot_coeffu1, rot_coeffu2,                            &
                   work3, work1,                                        &
                   lat_rot_np_in, long_rot_np_in,                       &
                    u_field_size )

! calculate lat/longitude for v points on equatorial C-grid
! longitude in work5, latitude in work6
      IF (l_regular) THEN
        DO j = 1-halo_j, n_rows+halo_j
          gj = l_datastart(2) + j - 1
          DO i = 1-halo_i, row_length+halo_i
            gi = l_datastart(1) + i - 1
            work5(i,j) = (base_lambda + (gi-1) * delta_lambda)            &
                       / pi_over_180
            work6(i,j) = (base_phi + (gj-.5) * delta_phi)                 &
                       / pi_over_180
          END DO
        END DO
      ELSE
        DO j = 1-halo_j, n_rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            work5(i,j) =  lambda_p (i)/ pi_over_180
            work6(i,j) =  phi_v (j)/ pi_over_180
          END DO
        END DO
      END IF   ! L_regular

! true_longitude in work7, true_latitude in work8
      CALL eqtoll(work6, work5,                                         &
                  work8, work7,                                         &
                  lat_rot_np_in, long_rot_np_in,                        &
                  v_field_size )

! Calculate rotation coefficients for v-points
! true_longitude in work7, rot_longitude in work5
      CALL w_coeff(rot_coeffv1, rot_coeffv2,                            &
                 work7, work5,                                          &
                lat_rot_np_in, long_rot_np_in,                          &
                 v_field_size )

      DO  k = 1, model_levels

! Copy  u  at v points into work5 array
! Copy  v  at u points into work1 array
! Can just copy since u=constant1 and v=constant2 on a level
        DO j = 1-halo_j, n_rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            work5(i,j) = u_adv(i,j,k)
            work1(i,j) = v_adv(i,j,k)
          END DO
        END DO
! Need v  at last row so just copy from penultimate row
        j = rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          work1(i,j) = v_adv(i,j-1,k)
        END DO

! Find rotated u-component
        CALL w_lltoeq(rot_coeffu1,rot_coeffu2,                          &
           u_adv(1-halo_i,1-halo_j,k), work1,                           &
           work3, work4, u_field_size, .TRUE.)
!  u field returned in work3 - work4 is v field - not needed

! Copy into u_adv
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            u_adv(i,j,k) = work3(i,j)
          END DO
        END DO

! Find rotated v-component
        CALL w_lltoeq(rot_coeffv1,rot_coeffv2,                          &
           work5, v_adv(1-halo_i,1-halo_j,k),                           &
           work7, work8, v_field_size, .TRUE.)
! Non-halo v field returned in work8 - work7 is u field - not needed

! Copy into v_adv
        DO j = 1-halo_j, n_rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            v_adv(i,j,k) = work8(i,j)
          END DO
        END DO

      END DO  ! k = 1, model_levels

    END IF   ! on (uv_horizfn_number)

  END IF   ! L_rotate_winds


!----------------------------------------------------------------
! Section 5.3.3  Balance pressure from wind field
!----------------------------------------------------------------
! SCANIA balance pressure field balanced from u-field
! NB  v must = 0 for this spherical geometry version
!----------------------------------------------------------------

  IF (l_pressure_balance) THEN

    IF (.NOT. l_rotating) THEN
      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot balance pressure field'    &
          //' with winds if rotation is off (L_rotating=.false.)'
      END IF  !(me  ==  0)

    ELSE IF (uv_horizfn_number  ==  uv_horiz_balance) THEN
      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot balance pressure field'    &
          //' (L_rotate_winds=.true.) if wind is balanced from'         &
          //' the input pressure field'                                 &
          //' (i.e. uv_horizfn_number=uv_horiz_balance).'
      END IF  !(me  ==  0)

    ELSE IF (uv_horizfn_number  ==  uv_horiz_deform) THEN

      IF (me  ==  0) THEN
        WRITE(6,*) 'IDL_Initial_Data: Cannot balance pressure field'    &
          //' (L_rotate_winds=.true.) if deformation wind field'        &
          //' (uv_horizfn_number=uv_horiz_deform) is chosen.'
      END IF  !(me  ==  0)

    ELSE

      IF (me  ==  0) THEN
        WRITE(6,FMT='(A)')'** Pressure field balanced from wind field **'

        IF (model_domain  ==  mt_global) THEN
          WRITE(6,FMT='(A)')'**  u is constant or function of latitude only**'
        ELSE ! Various LAM configurations
          WRITE(6,FMT='(A)')'** u is constant or function of latitude on '
          WRITE(6,FMT='(A)')'   ROTATED LAM grid'
        END IF     !  model_domain  ==  mt_global

      END IF   !(me  ==  0)

! DEPENDS ON: idl_pr_balance
      CALL idl_pr_balance(                                              &
                        model_domain, row_length, rows, n_rows,         &
                        model_levels,                                   &
                        delta_x, delta_y, halo_i, halo_j,               &
                        me, l_datastart,                                &
                        delta_lambda, delta_phi,                        &
                        lambda_p, phi_p,                                &
                        l_regular,                                      &
                        base_phi, base_lambda,                          &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        theta_eh, exner_rho_levels_eh, u_adv, v_adv,    &
                        theta_ref, exner_ref,                           &
! Profile settings
                        height_domain, theta_surface, brunt_vaisala,    &
                        u_in, v_in, u_ramp_start, u_ramp_end,           &
                        ujet_lat, ujet_width,                           &
                        f_plane, r_plane,                               &
!  Options
                        l_cartesian, l_code_test)

    END IF

  END IF  !  L_pressure_balance


! ----------------------------------------------------------------------

ELSE    !    tprofile_number  ==  tp_dump

  IF(me  ==  0)THEN
    WRITE(6,FMT='(A,I2)')'tprofile_number = ',tprofile_number
    WRITE(6,FMT='(A)')'Vertical profile unchanged - as in input dump '
  END IF   !(me  ==  0)

END IF    ! tprofile_number  /=  tp_dump

!-----------------------------------------------------------------
! Output selected diagnostics

!-------------------------------------------------------------------

IF (me  ==  0) THEN

  IF(first_constant_rho_level  >   model_levels)THEN
    WRITE(cmessage,*)                                                 &
     'ERROR: lower level for first_constant_rho_level above model top'
    errorstatus = 1

    CALL ereport( routinename, errorstatus, cmessage )
  END IF    !first_constant_rho_level  >   model_levels

  IF ( l_cartesian ) THEN
    IF( 2.0 * z_orog_print(0)   >                                       &
     r_rho_levels(1,1,first_constant_rho_level))THEN
      WRITE(cmessage,*)                                                 &
       'ERROR: flattening the levels too quickly'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
     END IF  
  ELSE  
    IF( 2.0 * z_orog_print(0)   >                                       &
     r_rho_levels(1,1,first_constant_rho_level) - earth_radius)THEN
      WRITE(cmessage,*)                                                 &
       'ERROR: flattening the levels too quickly'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
     END IF
  END IF   ! flattening test

  WRITE (6,*) ' ====================== End of Idealised ',            &
              'Settings ======================='
  WRITE (6,*) ' '
END IF        ! (me  ==  0)

IF (lhook) CALL dr_hook('EG_IDL_INITIAL_DATA',zhook_out,zhook_handle)

END SUBROUTINE eg_idl_initial_data
