! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_ni_init_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_ni_init(  pi                                        &
,model_domain, row_length, rows, n_rows, model_levels,wet_levels      &
,tr_vars, tr_levels, bl_levels, first_constant_r_rho_level            &
,cos_theta_latitude, sec_theta_latitude, f3_at_u, f3_at_v             &
,timestep, first_atmstep_call, l_regular                              &
,delta_x, delta_y, delta_lambda, delta_phi, base_phi, base_lambda     &
,lat_rot_np_in, long_rot_np_in, phi_p, phi_v                          &
,r_theta_levels, r_rho_levels, r_at_u, r_at_v, r_at_u_w, r_at_v_w     &
,z_orog_print, eta_theta_levels, eta_rho_levels                       &
! Multi-processor
,offx,offy,halo_i,halo_j, mype, nproc, at_extremity, datastart        &
,gc_all_proc_group, global_row_length, global_rows                    &
,g_rows, g_row_length, nproc_x                                        &
! Primary fields
,theta, rho,   exner_rho_levels, q, qcl, qcf, qcf2, qrain, qgraup     &
,cf, cfl, cff, u, v, w, u_adv, v_adv, w_adv                           &
! LAM lateral boundary data
! Grid info for idealised
,height_domain_in, height_domain, big_layers                          &
,transit_layers, surface_type, p_surface                              &
! Profile settings
,tprofile_number, qprofile_number, uvprofile_number                   &
,brunt_vaisala, theta_surface, dtheta_dz1, height_dz1                 &
,u_in, v_in, height_u_in, ujet_lat, ujet_width, u_ramp_start          &
,u_ramp_end, f_plane, r_plane, q1, num_profile_data                   &
,zprofile_data, tprofile_data, qprofile_data, num_uvprofile_data      &
,z_uvprofile_data, uprofile_data, vprofile_data                       &
! Forcing settings
,max_model_levels, max_num_profile_data, max_num_force_times          &
,tforce_option, qforce_option, uvforce_option, num_tforce_levels      &
,num_tforce_times, num_qforce_levels, num_qforce_times                &
,num_uvforce_levels, num_uvforce_times, z_tforce_data                 &
,tforce_data, z_qforce_data, qforce_data, z_uvforce_data              &
,uforce_data, vforce_data, tforce_data_modlev                         &
,qforce_data_modlev, uforce_data_modlev, vforce_data_modlev           &
! Dynamical core settings
,suhe_pole_equ_deltat, suhe_static_stab                               &
,base_frictional_timescale, frictional_timescale                      &
,suhe_sigma_cutoff, suhe_level_weight, l_sh_williamson                &
!  Horizontal function parameters
,t_horizfn_number, uv_horizfn_number, t_horizfn_data                  &
,l_perturb_t, perturb_magnitude_t,l_perturb_q                         &
,perturb_magnitude_q, l_perturb_correlate_tq                          &
,l_perturb_correlate_vert, l_perturb_correlate_time                   &
,perturb_type, perturb_height                                         &
!  Profiles for fixed lbcs and sponge zones
,u_ref, v_ref, theta_ref, exner_ref, rho_ref, zprofile_orog,          &
idl_interp_option, hf                                                 &
!  Options
,l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                                &
,l_constant_dz, l_rotating, l_fixed_lbcs, l_polar_wind_zero           &
,l_wind_balance, l_rotate_winds, l_pressure_balance, l_physics        &
,l_dry, l_sponge,l_cartesian,  l_code_test,l_cyclone, l_baroclinic    &
,h_print, timestep_number, h_o_actual, grow_steps                     &
,h_o_per_step, h_o, grid_number,  grid_flat                           &
,first_theta_height, thin_theta_height, big_factor, vert_grid_ratio   &
,lambda_fraction, phi_fraction, half_width_x, half_width_y            &
,idl_max_num_bubbles, idl_bubble_option, idl_bubble_max               &
,idl_bubble_height, idl_bubble_xoffset, idl_bubble_yoffset            &
,idl_bubble_width, idl_bubble_depth, l_idl_bubble_saturate            &
! Other fields
,a_ixsts,len_a_ixsts, a_spsts,len_a_spsts, nproc_y                    &
,gc_proc_row_group,  idlsurffluxseaoption                             &
,idlsurffluxseaparams, l_flux_bc, flux_h, flux_e                      &
,i_hour, i_minute, i_second, problem_number,  orography )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE water_constants_mod
USE eg_idl_generate_grid_mod
USE atm_fields_bounds_mod
USE ereport_mod, ONLY : ereport
USE UM_ParParams
USE problem_mod, ONLY: standard, monsoon, dynamical_core,             &
                       idealised_problem, standard_namelist
USE eg_idl_calc_eta_levels_mod
IMPLICIT NONE
!
! Description:
!   this routine does the initialisation for idealised test runs.
!   Copy of original idl_ni_init with changes for ENDGame
!  
!
! Method: Currently the subroutine generates orography and height
!         fields. Initial conditions need to be calculated at a
!         later stage in eg_SISL_setcon() when other model
!         parameters are specified.
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
, wet_levels                                                          &
                   ! number of model levels where moisture
                   ! variables are held
, bl_levels                                                           &
                   ! number of  boundary_layer_levels
, first_constant_r_rho_level                                          &
, halo_i                                                              &
                   ! Size of halo in i direction.
, halo_j                                                              &
                   ! Size of halo in j direction.
, offx                                                                &
                   ! Size of small halo in i
, offy                                                                &
                   ! Size of small halo in j.
, tr_vars                                                             &
, tr_levels



LOGICAL, INTENT(IN) :: first_atmstep_call

REAL                                                                  &
  timestep                                                            &
, height_domain                                                       &
, zprofile_orog, hf                                                   &
, dtheta_dz1(3)                                                       &
                  !Allows different values of dtheta_dz to be set
, height_dz1(2)                                                       &
                  ! at different heights specified by the
                  ! height_dz variable.
, theta_ref(model_levels)                                             &
                           !theta profile for use in sponge     & lbcs
, exner_ref(model_levels + 1)                                         &
                               ! Exner profile for use in lbcs
, rho_ref(model_levels)                                               &
                          ! rho profile for use in lbcs
, u_ref(model_levels)                                                 &
                        ! u profile for use in lbcs
, v_ref(model_levels)   ! u_adv profile for use in lbcs

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

! Variables for idealised surface fluxes
INTEGER, INTENT(IN) :: idlsurffluxseaoption ! Surface flux option
REAL, INTENT(IN)    :: idlsurffluxseaparams(10) ! Idl flux params
REAL, INTENT(INOUT) :: flux_h(row_length, rows) ! Idealised
REAL, INTENT(INOUT) :: flux_e(row_length, rows) ! surface fluxes
LOGICAL, INTENT(INOUT) :: l_flux_bc ! Switch for specified fluxes

! Time variables
INTEGER, INTENT(IN) :: i_hour    ! Model time (UTC, hour)
INTEGER, INTENT(IN) :: i_minute  ! Model time (UTC, minute)
INTEGER, INTENT(IN) :: i_second  ! Model time (UTC, second)

INTEGER                                                               &
  mype                                                                &
               ! My processor number
, nproc                                                               &
            ! Total number of processors
, nproc_y                                                             &
              ! number of processors in N-S direction
, gc_all_proc_group                                                   &
                    ! Group identifier for all processors.
, gc_proc_row_group                                                   &
                    ! Group identifier for all processors.
, datastart(2)                                                        &
                     ! First gridpoints held by this processor
, global_row_length                                                   &
                       ! global number of points on a row
, global_rows          ! global number of rows

INTEGER                                                               &
  g_rows(0:nproc-1)                                                   &
, g_row_length(0:nproc-1)                                             &
, nproc_x

LOGICAL                                                               &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Include physical constants
REAL                                                                  &
     ! physical constants
 pi    


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
, r_at_u_w(1-halo_i:row_length+halo_i,                                &
               1-halo_j:rows+halo_j, 0:model_levels)                  &
, r_at_v_w(1-halo_i:row_length+halo_i,                                &
               1-halo_j:n_rows+halo_j, 0:model_levels)                &
, eta_theta_levels(0:model_levels)                                    &
, eta_rho_levels(model_levels)                                        &
, cos_theta_latitude(1-offx:row_length+offx,                          &
                      1-offy:rows+offy)                               &
, sec_theta_latitude(1-offx:row_length+offx,                          &
                      1-offy:rows+offy)                               &
, f3_at_u (1-offx:row_length+offx,                                    &
                      1-offy:rows+offy)                               &
, f3_at_v (1-offx:row_length+offx,                                    &
                      1-offy:n_rows+offy)                             &
, z_orog_print(0:model_levels)


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

REAL                                                                  &
  orography(row_length,rows)

REAL, INTENT (INOUT) ::                                               &
  rho(1-offx:row_length+offx, 1-offy:rows+offy,                       &
        model_levels)                                                 &
, exner_rho_levels(1-offx:row_length+offx, 1-offy:rows+offy,          &
        model_levels+1)                                               &
, theta(1-offx:row_length+offx, 1-offy:rows+offy,                     &
          model_levels)                                               &
, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
      0:model_levels)                                                 &
, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
        0:model_levels)                                               &
, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
        0:model_levels)                                               &
, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
      0:model_levels)                                                 &
, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
      0:model_levels)                                                 &
, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
      0:model_levels)                                                 &
, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                &
      wet_levels)                                                     &
, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
        wet_levels)                                                   &
, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
        wet_levels)                                                   &
, u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)           &
, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)         &
, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)         &
, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        model_levels)                                                 &
, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,           &
        model_levels)                                                 &
, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        0:model_levels)


! Diagnostic variables


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
INTEGER len_a_ixsts
INTEGER len_a_spsts
INTEGER a_ixsts(len_a_ixsts)     ! stash index array
REAL    a_spsts(len_a_spsts)     ! atmos stash array

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
REAL :: height_domain_in
REAL :: delta_x           ! Resolution at equator
REAL :: delta_y           ! Resolution at equator
REAL :: big_factor
REAL :: vert_grid_ratio
REAL :: first_theta_height
REAL :: thin_theta_height
REAL :: p_surface
REAL :: theta_surface
REAL :: brunt_vaisala
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

REAL :: r_plane        ! reference latitude for row 1 (bottom row)

! Namelist profile data
REAL :: zprofile_data(max_num_profile_data)    ! heights for t,q
REAL :: tprofile_data(max_num_profile_data)    ! theta profile
REAL :: qprofile_data(max_num_profile_data)    ! humidity profile
REAL :: z_uvprofile_data(max_num_profile_data) ! heights for u,v
REAL :: uprofile_data(max_num_profile_data)    ! u-wind profile
REAL :: vprofile_data(max_num_profile_data)    ! v-wind profile

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
INTEGER :: timestep_number
INTEGER :: problem_number
INTEGER :: idl_interp_option  ! Profile interpolation option

LOGICAL :: l_constant_dz ! Sets const dtheta_dz from dtheta_dz1(1)
LOGICAL :: l_cartesian   ! new  cartesian flag (replaces l_trivial_trigs!)
LOGICAL :: l_fixed_lbcs    ! Set fixed lateral boundary conditions
LOGICAL :: l_pressure_balance  ! Geostrophically balance pressures
LOGICAL :: l_wind_balance  ! Geostrophically balance initial winds
LOGICAL :: l_rotate_winds  ! rotate input u,v (true) to LAM u,v
LOGICAL :: l_polar_wind_zero   ! set u=0 on polar row
LOGICAL :: l_rotating      ! .true. for Earth's rotation
LOGICAL :: l_code_test     ! User switch for testing code
LOGICAL :: l_cyclone       ! true if cyclone simulation
LOGICAL :: l_baroclinic    ! true if baroclinic wave simulation
LOGICAL :: l_physics       ! physics switch
LOGICAL :: l_dry           ! moisture switch
LOGICAL :: l_sponge        ! sponge switch
! Bubble perturbation options
INTEGER :: idl_max_num_bubbles
INTEGER :: idl_bubble_option(idl_max_num_bubbles)
REAL :: idl_bubble_max(idl_max_num_bubbles)
REAL :: idl_bubble_height(idl_max_num_bubbles)
REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
REAL :: idl_bubble_width(idl_max_num_bubbles)
REAL :: idl_bubble_depth(idl_max_num_bubbles)
LOGICAL :: l_idl_bubble_saturate(idl_max_num_bubbles)

LOGICAL :: l_mcr_qcf2    ! true if using 2nd prognostic cloud ice
LOGICAL :: l_mcr_qrain   ! true if using prognostic rain
LOGICAL :: l_mcr_qgraup  ! true if using prognostic graupel


INTEGER                                                               &
 errorstatus      ! Return code : 0 Normal Exit : >0 Error

CHARACTER(LEN=256)                                                         &
 cmessage         ! Error message if return code >0

! Local variables:

CHARACTER(LEN=*) routinename
PARAMETER (   routinename='IDL_NI_INIT')

! Work arrays with extended halos (_eh)
! Needed so that external halo values in LAMS can be set correctly
! for the lateral boundary arrays.

REAL                                                                  &
  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
           model_levels)                                              &
, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
           model_levels)                                              &
, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                     &
           1-halo_j:rows+halo_j, model_levels+1)


! Loop counters
INTEGER i,                                                            &
        j ,                                                           &
        ij,                                                           &
        k,                                                            &
        k2      ! Loop counters

REAL z_at_theta   ! height of theta level in forcing data interp
REAL z_at_rho     ! height of rho level in forcing data interp
REAL weight       ! weight used in forcing data interpolation
REAL, ALLOCATABLE, SAVE :: orog_per_step(:,:)

! Variables for idealised surface fluxes
REAL :: hrl,                                                          &
        xfact

INTEGER ierr

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_NI_INIT',zhook_in,zhook_handle)


! ---------------------------------------------------------------------

IF (first_atmstep_call) THEN

  IF (mype == 0) THEN
    WRITE (6,*) ' =============================================',     &
                '=========================='
    WRITE (6,*) ' '
    WRITE (6,*) '                         UM IDEALISED SETTINGS'
    WRITE (6,*) ' '
    WRITE (6,*) ' =============================================',     &
                '=========================='
    WRITE (6,*) ' '
    WRITE (6,*) ' VERTICAL GRID AND OROGRAPHY '
  END IF

! h_o_actual needs to be set>0 to ensure correct surface type generated
!  when not setting h_o in namelist
  h_o_actual = 10.0     ! any value > 0.0
  IF ( grow_steps  >   1)THEN
    IF (.not.ALLOCATED(orog_per_step)) THEN
      ALLOCATE (orog_per_step(row_length, rows),STAT=ierr)
      IF (ierr /= 0) THEN
         CALL ereport( 'EG_IDL_NI_INIT', ierr, 'Allocation failed' )
      END IF
    END IF

    IF (mype == 0) THEN
      WRITE (6,*) '   Growing Orography '
    END IF

    IF ( surface_type   ==  surface_mask .OR.                         &
          surface_type   ==  surface_dump ) THEN
! store orography amount per timestep
      DO j = 1, rows
        DO i = 1, row_length
          orog_per_step(i,j) = orography(i,j ) / grow_steps
! reset orography for first timestep
          orography(i,j )  = orog_per_step(i,j)
        END DO
      END DO
    ELSE    ! growing idealised data
      h_o_per_step =  h_o / grow_steps
      h_o_actual =  h_o_per_step
      IF (mype == 0) THEN
        WRITE (6,fmt='(A,I5,A)') '   Hill/mountain grows over ',       &
                    grow_steps,                                       &
                    ' timesteps'
        WRITE (6,fmt='(A,E16.8,A)') '   final hill/mountain height ',     &
                                h_o,' metres'
      END IF   !(mype == 0)
    END IF ! surface_type   ==  surface_mask .or.
!                   surface_type   ==  surface_dump
  ELSE  ! grow_steps  <=  1
    IF ( surface_type   ==  surface_mask .OR.                         &
          surface_type   ==  surface_dump ) THEN
     IF(mype  ==  0)THEN
        WRITE(6,*)'Orography fixed from start  '
      END IF   !(mype  ==  0)
    ELSE  ! idealised orography
      h_o_per_step = 0.0
      h_o_actual = h_o
      IF(mype  ==  0)THEN
        WRITE(6,fmt='(A,E16.8,A)')'hill/mountain fixed MAX height '       &
               ,h_o,' metres'
      END IF   !(mype  ==  0)
     END IF ! surface_type   ==  surface_mask .or.
!                    surface_type   ==  surface_dump
 END IF  ! grow_steps  >   1

  IF ( grid_number  /=  vert_dump  .AND.       &
       grid_number  /=  vert_idl_um_grid ) THEN  ! if change grid

!   Calculate eta_levels (normalized vertical grid)
!    Only needs to be done at start
    CALL eg_idl_calc_eta_levels(                                      &
                      model_levels,mype,first_constant_r_rho_level,   &
                      eta_theta_levels, eta_rho_levels,               &
!  Grid information
                      grid_number, height_domain,                     &
                      first_theta_height, thin_theta_height,          &
                      big_layers, transit_layers, big_factor,         &
                      vert_grid_ratio)

!    height_domain set for idealised problem may have been changed
!    in Calc_eta_levels. Anyway, dump value of model top needs changing
    height_domain_in = height_domain
  END IF        ! grid_number  /=  vert_dump

!   Calculate required orography
!  Use h_print to pass in h_o_actual since overwritten if real orography
!                                        (surface_mask or surface_dump)
  h_print = h_o_actual
! DEPENDS ON: eg_idl_surface_setup
    CALL eg_idl_surface_setup(                                        &
                    row_length, rows, model_levels                    &
,                   halo_i, halo_j                                    &
,                   mype, nproc, at_extremity, model_domain           &
,                   n_rows                                            &
,                   orography, r_theta_levels                         &
,                   r_at_u_w(1-halo_i,1-halo_j,0)                     &
,                   r_at_v_w(1-halo_i,1-halo_j,0)                     &
,                   surface_type                                      &
,                   h_print, lambda_fraction, phi_fraction            &
,                   half_width_x, half_width_y)

  IF ( grow_steps  >   1)THEN
! For growing real orography send final max value into Generate_grid
!  so grid information over orography can be printed
    IF ( surface_type   ==  surface_mask .OR.                         &
         surface_type   ==  surface_dump ) THEN
      h_print = h_print * grow_steps
    END IF   ! surface_type   ==  surface_mask .or.
!                   surface_type   ==  surface_dump
  END IF    !  grow_steps  >   1
!   Now generate r_theta_levels and r_rho_levels

    CALL eg_idl_generate_grid(                                        &
                   first_constant_r_rho_level                         &
,                  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts           &
,                  grid_number, grid_flat,height_domain_in            &
,                  h_print, z_orog_print)

  IF(mype  ==  0) THEN
    WRITE(6,fmt='(A,I5)')'Eta_grid generated problem_number ',         &
                         problem_number
  END IF

! Copy single halo arrays into extended halo arrays
  DO k = 1, model_levels
    DO j = 1-offy, rows+offy
      DO i = 1-offx, row_length+offx
        theta_eh(i,j,k) = theta(i,j,k)
        rho_eh(i,j,k)   = rho(i,j,k)
        exner_rho_levels_eh(i,j,k)   = exner_rho_levels(i,j,k)
      END DO
    END DO
  END DO
  k = model_levels + 1
  DO j = 1-offy, rows+offy
    DO i = 1-offx, row_length+offx
      exner_rho_levels_eh(i,j,k)   = exner_rho_levels(i,j,k)
    END DO
  END DO


  IF (timestep_number == 1) THEN
!   Generate intitial profiles for fields
!   For grid_number = vert_dump only prints out grid information
! DEPENDS ON: eg_idl_initial_data
    CALL eg_idl_initial_data(                                         &
                   model_domain, row_length, rows, n_rows             &
,                  model_levels, wet_levels                           &
,                  tr_vars,tr_levels, bl_levels                       &
,                  first_constant_r_rho_level                         &
,                  cos_theta_latitude, sec_theta_latitude             &
,                  f3_at_u, f3_at_v, timestep                         &
,                  delta_x, delta_y                                   &
,                  offx, offy, halo_i, halo_j                         &
,                  mype, nproc, at_extremity                          &
,                  datastart, gc_all_proc_group                       &
,                  global_row_length, global_rows                     &
,                  delta_lambda, delta_phi                            &
!  VarRes Grid Spacing
,                  lambda_p, phi_p, lambda_u, phi_v                   &
,                  l_regular                                          &
,                  base_phi, base_lambda                              &
,                  lat_rot_np_in, long_rot_np_in                      &
,                  r_theta_levels, r_rho_levels                       &
,                  r_at_u, r_at_v, z_orog_print                       &
,                  eta_theta_levels, eta_rho_levels                   &
,                  theta_eh, rho_eh, exner_rho_levels_eh              &
,                  q, qcl, qcf, qcf2, qrain, qgraup                   &
,                  u_adv, v_adv, w_adv                                &
,                  g_rows, g_row_length                               &
,                  nproc_x, nproc_y                                   &
!  Grid information
,                  height_domain_in                                   &
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

  ! Copy idealised data from extended halo work arrays into
  ! model fields and set other fields.



    u(1:row_length,1:rows,:)   = u_adv(1:row_length,1:rows,:)
    v(1:row_length,1:n_rows,:) = v_adv(1:row_length,1:n_rows,:)
    w(1:row_length,1:rows,:)   = w_adv(1:row_length,1:rows,:)

    IF(model_domain  ==  mt_global)THEN
! DEPENDS ON: polar_reset_mean
      CALL polar_reset_mean(                                          &
                      exner_rho_levels,rho,theta,w,                   &
                      q,qcl,qcf,                                      &
                      cf,cfl,cff,                                     &
                      row_length, rows, model_levels,                 &
                      wet_levels, global_row_length,                  &
                      offx, offy, halo_i, halo_j,                     &
                      nproc, nproc_y, gc_proc_row_group,              &
                      at_extremity)
    END IF

  END IF ! timestep_number == 1

END IF ! on first_atmstep_call
!----------------------------------------------------------------------
!                End of section for Timestep = 1
!----------------------------------------------------------------------

!----------------------------------------------------------------------

!              Start of section for first call to ATMSTEP
!         (i.e. first timestep of a new run OR continuation run)

!----------------------------------------------------------------------

IF (first_atmstep_call) THEN

  !-----------------------------------------------------------
  ! Interpolate theta forcing data to model levels if required
  !-----------------------------------------------------------

  IF (tforce_option  >   0) THEN

    ! Check to make sure the namelist profile data extends
    ! to the top of the model.
    IF (eta_theta_levels(model_levels)*height_domain                  &
         >   z_tforce_data(num_tforce_levels)) THEN
      WRITE(cmessage,*)                                               &
        'Idealised namelist forcing profile data (T)'                 &
        //'does not extend to the top of the model.'                  &
        //'Please modify the namelist data.'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
    END IF

    ! Interpolate theta from namelist profile to model levels
    DO k = 1, model_levels
      z_at_theta = eta_theta_levels(k) * height_domain
      DO k2 = 1, num_tforce_levels-1
        DO j = 1, num_tforce_times

            IF (z_at_theta  >   z_tforce_data(k2) .AND.               &
                z_at_theta  <=  z_tforce_data(k2+1)) THEN

              weight = (z_at_theta - z_tforce_data(k2))               &
                     /(z_tforce_data(k2+1) - z_tforce_data(k2))

              tforce_data_modlev(k,j) = tforce_data(k2,j)             &
               + weight*(tforce_data(k2+1,j) - tforce_data(k2,j))
            END IF

        END DO
      END DO
    END DO

  END IF

  !-----------------------------------------------------------
  ! Interpolate q forcing data to model levels if required
  !-----------------------------------------------------------

  IF (qforce_option  >   0) THEN

    ! Check to make sure the namelist profile data extends
    ! to the top of the model.
    IF (eta_theta_levels(model_levels)*height_domain                  &
         >   z_qforce_data(num_qforce_levels)) THEN
      WRITE(cmessage,*)                                               &
        'Idealised namelist forcing profile data (q)'                 &
        //'does not extend to the top of the model.'                  &
        //'Please modify the namelist data.'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
    END IF

    ! Interpolate humidity from namelist profile to model levels
    DO k = 1, model_levels
      z_at_theta = eta_theta_levels(k) * height_domain
      DO k2 = 1, num_qforce_levels-1
        DO j = 1, num_qforce_times

            IF (z_at_theta  >   z_qforce_data(k2) .AND.               &
                z_at_theta  <=  z_qforce_data(k2+1)) THEN

              weight = (z_at_theta - z_qforce_data(k2))               &
                       /(z_qforce_data(k2+1) - z_qforce_data(k2))

              qforce_data_modlev(k,j) = qforce_data(k2,j)             &
               + weight*(qforce_data(k2+1,j) - qforce_data(k2,j))
            END IF

        END DO
      END DO
    END DO

  END IF

  !-----------------------------------------------------------
  ! Interpolate u,v forcing data to model levels if required
  !-----------------------------------------------------------

  IF (uvforce_option  >   0) THEN

    ! Check to make sure the namelist profile data extends
    ! to the top of the model.
    IF (eta_theta_levels(model_levels)*height_domain                  &
         >   z_uvforce_data(num_uvforce_levels)) THEN
      WRITE(cmessage,*)                                               &
        'Idealised namelist forcing profile data (uv)'                &
        //'does not extend to the top of the model.'                  &
        //'Please modify the namelist data.'
      errorstatus = 1

      CALL ereport( routinename, errorstatus, cmessage )
    END IF

    ! Interpolate wind from namelist profile to model levels
    DO k = 1, model_levels
      z_at_rho = eta_rho_levels(k) * height_domain
      DO k2 = 1, num_uvforce_levels-1
        DO j = 1, num_uvforce_times

            IF (z_at_rho  >   z_uvforce_data(k2) .AND.                &
                z_at_rho  <=  z_uvforce_data(k2+1)) THEN

              weight = (z_at_rho - z_uvforce_data(k2))                &
                     /(z_uvforce_data(k2+1) - z_uvforce_data(k2))

              uforce_data_modlev(k,j) = uforce_data(k2,j)             &
               + weight*(uforce_data(k2+1,j) - uforce_data(k2,j))
              vforce_data_modlev(k,j) = vforce_data(k2,j)             &
               + weight*(vforce_data(k2+1,j) - vforce_data(k2,j))
            END IF

        END DO
      END DO
    END DO

  END IF
!----------------------------------------------------------------------
! Section to cal idl_fix_lam_lbcs deleted
!----------------------------------------------------------------------
END IF ! on first_atmstep_call = .true.
!----------------------------------------------------------------------
!              End of section for first call to ATMSTEP
!----------------------------------------------------------------------

!----------------------------------------------------------------------

!              Start of section for timestep > 0

!----------------------------------------------------------------------
IF ( timestep_number > 0) THEN

  !--------------------------------------------------------------
  !
  !                   Growing Orography
  !
  !--------------------------------------------------------------
  IF ( timestep_number <= grow_steps ) THEN

!  If growing orography then need to recalculate surface and
!   r vertical coordinate values (but eta values remain fixed)
    IF ( surface_type   ==  surface_mask .OR.                         &
         surface_type   ==  surface_dump ) THEN
      h_o_actual = h_print
! grow orography by  orog_per_step
      DO j = 1, rows
        DO i = 1, row_length
          ij = i + (j-1)*row_length
          orography(i,j) = orography(i,j ) +                          &
                                    orog_per_step(i,j)
        END DO
      END DO
      IF(mype  ==  0)THEN
        WRITE(6,*)'Growing real orography '
      END IF   !(mype  ==  0)
    ELSE  !  growing idealised orography
      h_o_actual = h_o_actual + h_o_per_step
      IF(mype  ==  0)THEN
         WRITE(6,fmt='(A,E16.8,A)')'hill/mountain growing ',h_o_actual    &
                               ,' metres'
      END IF   !(mype  ==  0)
    END IF !  surface_type   ==  surface_mask .or.
!                    surface_type   ==  surface_dump

! DEPENDS ON: eg_idl_surface_setup
      CALL eg_idl_surface_setup(                                      &
                    row_length, rows, model_levels                    &
,                   halo_i, halo_j                                    &
,                   mype, nproc, at_extremity, model_domain           &
,                   n_rows                                            &
,                   orography, r_theta_levels                         &
,                   r_at_u_w(1-halo_i,1-halo_j,0)                     &
,                   r_at_v_w(1-halo_i,1-halo_j,0)                     &
,                   surface_type                                      &
,                   h_o_actual, lambda_fraction, phi_fraction         &
,                   half_width_x, half_width_y)

!   Regenerate r_theta_levels and r_rho_levels

       CALL eg_idl_generate_grid(                                     &
                    first_constant_r_rho_level                        &
,                   a_ixsts,len_a_ixsts, a_spsts,len_a_spsts          &
,                   grid_number, grid_flat,height_domain_in           &
,                   h_o_actual, z_orog_print)

    IF(l_fixed_lbcs)THEN
      IF(mype  ==  0)THEN
        CALL ereport( 'EG_IDL_NI_INIT', 1,                            &
                  'Call routine to make fixed lbcs consistent' )
      END IF   !(mype  ==  0)
    END IF     ! L_fixed_lbcs

  END IF ! on (timestep_number <= grow_steps)


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
  IF (idlsurffluxseaoption == 1) THEN

    IF (mype == 0 .AND. first_atmstep_call) THEN
      WRITE(6,*) ' '
      WRITE(6,*) '  Specifing zero sea surface heat fluxes'
      WRITE(6,*) ' '
    END IF

    l_flux_bc = .TRUE.

    DO j = 1, rows
      DO i = 1, row_length
        flux_h(i,j) = 0.0
        flux_e(i,j) = 0.0
      END DO
    END DO

  END IF ! on (IdlSurfFluxSeaOption  ==  1)

  !---------------------------------
  ! Option 2: Diurnal cycle
  ! (positive surface fluxes during the day, zero at night)
  !---------------------------------
  IF (idlsurffluxseaoption == 2) THEN

    IF (mype == 0 .AND. first_atmstep_call) THEN
      WRITE(UNIT=6,FMT=*) ' SURFACE FLUXES'
      WRITE(UNIT=6,FMT=*)                                             &
        '  Specifing sea surface flux diurnal cycle'
      WRITE(UNIT=6,FMT='(A37,F7.2)')                                  &
        '   Maximum sensible heat flux (W/m2):',                      &
                    idlsurffluxseaparams(1)
      WRITE(UNIT=6,FMT='(A37,F7.2)')                                  &
        '   Maximum latent heat flux (W/m2):  ',                      &
                    idlsurffluxseaparams(2)
      WRITE(UNIT=6,FMT='(A37,F7.2)')                                  &
        '   Time (UTC) of max flux (hours):   ',                      &
                    idlsurffluxseaparams(3)
      WRITE(UNIT=6,FMT='(A37,F7.2)')                                  &
        '   Length of the day (hours):        ',                      &
                    idlsurffluxseaparams(4)
      WRITE(UNIT=6,FMT=*) ' '
    END IF

    l_flux_bc = .TRUE.

    ! IdlSurfFluxSeaParams(1) = max sensible heat flux
    ! IdlSurfFluxSeaParams(2) = max latent heat flux
    ! IdlSurfFluxSeaParams(3) = Time (UTC) of max flux (hours)
    ! IdlSurfFluxSeaParams(4) = Length of day (hours)

    ! Calculate current time (UTC) in hours
    hrl   = ((i_hour*60.0 + i_minute)*60.0 + i_second)/3600.0

    ! Set up diurnally varying function
    xfact  = COS( pi/2.*(idlsurffluxseaparams(3)-hrl)/                &
                        (idlsurffluxseaparams(4)/2.)  )

    ! Limit fluxes to being positive (upward)
    IF (xfact <= 0.) xfact=0.

    ! Set diurnally varying fluxes
    DO j = 1, rows
      DO i = 1, row_length
        flux_h(i,j) = idlsurffluxseaparams(1)*xfact**1.5
        flux_e(i,j) = idlsurffluxseaparams(2)*xfact**1.3
      END DO
    END DO

  END IF ! on (IdlSurfFluxSeaOption  ==  2)

 END IF      ! timestep_number > 0
!----------------------------------------------------------------------
!              End of section for timestep > 0
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EG_IDL_NI_INIT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_ni_init
END MODULE eg_idl_ni_init_mod
