! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables and allocatable arrays which are local to atm_step 
!  (and below)
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

MODULE atm_step_local

IMPLICIT NONE

LOGICAL :: l_rad_step       ! T :activate radiation this timestep
LOGICAL :: l_rad_step_diag  ! T :activate fast radiation this timestep 
LOGICAL :: l_rad_step_prog  ! T :activate slow radiation this timestep
LOGICAL :: l_do_rims

LOGICAL :: l_do_inc_vels
LOGICAL :: l_do_inc_rho     ! flag for rho incr diagnostic

INTEGER :: cycleno ! Number of cycles (iterations) for iterative SISL

INTEGER :: tr_size ! size of a single tracer - FLUME variable

LOGICAL :: l_tracer ! T if *any* tracer variables present

LOGICAL :: l_print_l2norms ! diagnostic printing of l2norms

INTEGER :: first_constant_r_rho_level    ! 1st rho level with r constant
INTEGER :: first_constant_r_rho_level_m1 ! Max (1,first_constant_r_rho_level)

INTEGER :: i_start
INTEGER :: i_stop 
INTEGER :: j_start
INTEGER :: j_stop
INTEGER :: j_begin
INTEGER :: j_end

INTEGER :: lambda_start ! pointer for start of lambda_p/lambda_u on this pe

! used in interpolation code - no. of levels to check for whether the departure 
! point lies inside the orography.
INTEGER :: check_bottom_levels 

INTEGER :: rhc_row_length
INTEGER :: rhc_rows

!  Dummy variables for SCM Diagnostics,
!  for passing down to atmos_physics1 and 2 routines
INTEGER, PARAMETER :: nscmdpkgs = 12
LOGICAL :: l_scmdiags(nscmdpkgs)
LOGICAL :: l_flux_bc    ! T if prescribed surface fluxes to be used

! time-stepping weights. Values may be different
! at different cycles (UMUI controlled).
REAL :: alpha1, alpha2, alpha3, alpha4

REAL :: increment_factor ! For calculating value of LBC at next TL

! Variables required for call to SET_LATERAL_BOUNDARIES
INTEGER :: lbc_size        ! size of a single level of LBC
LOGICAL :: l_do_halos      ! update the halos?
LOGICAL :: l_do_boundaries ! update the boundaries?

! -------------------------------------------------

REAL, DIMENSION(:,:,:), ALLOCATABLE:: inc_u, inc_v, inc_w, inc_t, inc_rho 
REAL, DIMENSION(:,:,:), ALLOCATABLE:: inc_q, inc_qcl, inc_qcf, inc_cf, inc_cfl
REAL, DIMENSION(:,:,:), ALLOCATABLE:: inc_cff, inc_qrain, inc_qgraup, inc_qcf2

INTEGER :: nd_o3  ! Total size of ozone array supplied

! Code to do with tropopause diagnostics from O3_to_3D
! Declare logicals for Stash Diagnostics used in call to O3_to_3D
LOGICAL :: l_o3_trop_level   ! stash code 2,280
LOGICAL :: l_o3_trop_height  ! stash code 2,281
LOGICAL :: l_t_trop_level    ! stash code 2,282
LOGICAL :: l_t_trop_height   ! stash code 2,283

INTEGER :: land_pts_trif ! For dimensioning variables in NI_bl_ctl
INTEGER :: npft_trif     ! Depending on whether TRIFFID is in use 
 
REAL, POINTER :: cs(:)  => NULL() ! soil carbon content & accum soil respiration 
REAL, POINTER :: rsa(:) => NULL() ! (target varies according to l_TRIFFID)

INTEGER :: dim_cs1 ! soil C dimension: 1 for single, 4 for RothC
INTEGER :: dim_cs2 ! soil C dimension: 1 for single, LAND_FIELD for RothC

INTEGER :: co2_dim_len ! For dimension 3-D CO2 field to be passed
INTEGER :: co2_dim_row !     to NI_bl_ctl
INTEGER :: co2_dim_lev !     and NI_rad_ctl

! Array dimensions for sea-salt aerosol
INTEGER :: salt_dim1
INTEGER :: salt_dim2 
INTEGER :: salt_dim3

! Array dimensions for Aero_Ctl
INTEGER :: aero_dim1
INTEGER :: aero_dim2 
INTEGER :: aero_dim3

! Declare the tropopause variables output for O3_to_3D as allocatable
REAL, DIMENSION (:,:), ALLOCATABLE:: o3_trop_level
REAL, DIMENSION (:,:), ALLOCATABLE:: t_trop_level
REAL, DIMENSION (:,:), ALLOCATABLE:: o3_trop_height
REAL, DIMENSION (:,:), ALLOCATABLE:: t_trop_height

!3d ozone (expanded from zonal) for radiation
REAL, DIMENSION (:,:,:), ALLOCATABLE:: ozone3D

! Local Arrays to store microphysics fields
REAL, ALLOCATABLE, TARGET :: qcf2_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: qrain_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: qgraup_star(:,:,:)

! Local arrays for when using mixing ratios
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v, mix_cl, mix_cf
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v_inter, mix_cl_inter, mix_cf_inter
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf2, mix_rain, mix_graup  
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf2_inter, mix_rain_inter
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_graup_inter
REAL, DIMENSION (:,:,:), ALLOCATABLE :: q_store, qcl_store, qcf_store  
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf2_store, qrain_store, qgraup_store

!    Local arrays for when using mixing ratios
REAL, ALLOCATABLE, TARGET :: mix_v_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cl_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cf_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cf2_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_rain_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_graup_star(:,:,:)

!  stashworki = stashwork for section i
REAL, DIMENSION (:), ALLOCATABLE:: stashwork1,stashwork2,stashwork3
REAL, DIMENSION (:), ALLOCATABLE:: stashwork4,stashwork5,stashwork6
REAL, DIMENSION (:), ALLOCATABLE:: stashwork8,stashwork9,stashwork12
REAL, DIMENSION (:), ALLOCATABLE:: stashwork13,stashwork14,stashwork17
REAL, DIMENSION (:), ALLOCATABLE:: stashwork19,stashwork26,stashwork30
REAL, DIMENSION (:), ALLOCATABLE:: stashwork10,stashwork18,stashwork35
REAL, DIMENSION (:), ALLOCATABLE:: stashwork39

! Local dynamic arrays for phys1 and phys2 increments for moisture:
REAL, DIMENSION (:,:,:), ALLOCATABLE :: q_phys1, qcl_phys1, qcf_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: q_phys2, qcl_phys2, qcf_phys2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf2_phys1, qrain_phys1, qgraup_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf2_phys2, qrain_phys2, qgraup_phys2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v_phys1, mix_cl_phys1, mix_cf_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v_phys2, mix_cl_phys2, mix_cf_phys2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf2_phys1, mix_rain_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf2_phys2, mix_rain_phys2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_graup_phys1, mix_graup_phys2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cf_phys1, cfl_phys1, cff_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cf_phys2, cfl_phys2, cff_phys2

! local dynamic arrays for PC2
REAL, DIMENSION (:,:,:), ALLOCATABLE :: t_inc_pres, q_inc_pres, qcl_inc_pres
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf_inc_pres, cf_inc_pres 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cfl_inc_pres, cff_inc_pres 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: t_dini, q_dini, qcl_dini
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf_dini, cf_dini
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cfl_dini, cff_dini, rhts 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qtts, tlts, ptts

! Extra variables needed for cycling
! Vars ending in _phys1 are copies holding the value the original
! variable had after exiting phys1.
! Vars ending in _np1 are tn+1 estimates holding the value the original
! variable had at the end of the last cycle (provided that CycleNo>1).
! obtained from the last
! cycle when CycleNo>1.

REAL, DIMENSION (:,:,:), ALLOCATABLE :: r_u_phys1, r_v_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: thetastar_phys1, qstar_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qclstar_phys1, qcfstar_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf2_star_phys1, qrain_star_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qgraup_star_phys1, ti_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cca_phys1, area_cld_frac_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: bulk_cld_frac_phys1, bulk_cld_liq_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: bulk_cld_fr_phys1

REAL, SAVE, ALLOCATABLE, TARGET :: u_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: v_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: w_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: theta_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: rho_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: q_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: qcf_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: qcl_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: qcf2_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: qrain_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: qgraup_np1(:,:,:)

REAL, SAVE, ALLOCATABLE, TARGET :: etadot_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: thetav_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: exner_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: exner_surf_np1(:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_v_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cl_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cf_np1(:,:,:)                     
REAL, SAVE, ALLOCATABLE, TARGET :: m_r_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_gr_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cf2_np1(:,:,:)


REAL, DIMENSION (:,:), ALLOCATABLE:: z0msea_phys1, zh_phys1, ddmfx_phys1
REAL, DIMENSION (:,:), ALLOCATABLE:: deep_flag_phys1, past_precip_phys1,&
                                     past_conv_ht_phys1
REAL, DIMENSION (:,:), ALLOCATABLE:: t_land_ctile_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE:: t_sice_ctile_phys1
REAL, DIMENSION (:,:), ALLOCATABLE:: t_surf_phys1, t_sf_tile_phys1
REAL, DIMENSION (:,:), ALLOCATABLE:: snow_tile_phys1, dolr_phys1
REAL, DIMENSION (:,:), ALLOCATABLE:: rho_lbc_real_tend

INTEGER, DIMENSION (:,:), ALLOCATABLE :: ccb_phys1, cct_phys1

! Additional variables needed for cycling when mixing ratios are used

REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v_star_phys1, mix_cl_star_phys1 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf_star_phys1, mix_cf2_star_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_rain_star_phys1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_graup_star_phys1 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_v_np1, mix_cl_np1, mix_cf_np1 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_cf2_np1, mix_rain_np1
REAL, DIMENSION (:,:,:), ALLOCATABLE :: mix_graup_np1

! Variables required for ice category selection
REAL, POINTER :: p_ti(:,:,:)        => NULL()
REAL, POINTER :: p_ice_fract(:,:,:) => NULL()
REAL, POINTER :: p_ice_thick(:,:,:) => NULL()
REAL, POINTER :: p_tstar_sice(:,:,:) => NULL()
REAL, POINTER :: p_snodep_sice(:,:,:) => NULL()
REAL, POINTER :: p_ice_thick_rad(:,:,:) => NULL()
REAL, POINTER :: p_ice_fract_rad(:,:,:) => NULL()

! increment diagnostics:
REAL, DIMENSION (:,:,:), ALLOCATABLE :: u_incr_diagnostic, v_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: t_incr_diagnostic 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: q_incr_diagnostic 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcl_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf_incr_diagnostic 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qrain_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qgraup_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: qcf2_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cf_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cfl_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: cff_incr_diagnostic
REAL, DIMENSION (:,:,:), ALLOCATABLE :: w_incr_diagnostic

REAL, DIMENSION (:,:), ALLOCATABLE :: w_local_mask

! Allocatable arrays for use in AC_CTL call
REAL, DIMENSION (:,:,:), ALLOCATABLE :: work_q, work_qcl, work_qcf

! Workspace defined as allocatable arrays, since they each communicate
! fields between near adjacent calls and only need to use memory for
! a subset of the total routine.
REAL, DIMENSION (:,:,:), ALLOCATABLE :: exner_prime ! soln to helmholtz solver
REAL, DIMENSION (:,:,:), ALLOCATABLE :: dtheta_dr_term 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: depart_lambda, depart_phi 
REAL, DIMENSION (:,:,:), ALLOCATABLE :: depart_r_theta, depart_r_w

REAL, ALLOCATABLE :: rhcpt(:,:,:)

! ----------------------------------------------------------------

LOGICAL ::  gather   ! Convert to sea_ice points within BL code
PARAMETER ( gather = .true. ) ! (was l_compress_seaice)

INTEGER :: i, j, k, l   ! loop counters
INTEGER :: ij           ! 2D array index for offx/y variables
INTEGER :: lbc_pt       ! pointer to LBC array
INTEGER :: ij_u         ! 2D array index for U offx/y variables
INTEGER :: ij_v         ! 2D array index for V offx/y variables
INTEGER :: ji           ! 2D array index for halo_i/j variables
INTEGER :: icount       ! diagnostic count
INTEGER :: j_ti, j_ice_fract, j_ice_thick

INTEGER :: n_y_arrays    ! = 1 for global, 3 for LAM
INTEGER :: n_yw_arrays   ! = 1 for global, 2 for LAM
INTEGER :: n_yd_arrays   ! = 1 for global, 3 for LAM
INTEGER :: n_ydw_arrays  ! = 1 for global, 2 for LAM
!
INTEGER :: i_field       ! fields in multivariate swapbounds

INTEGER :: info          ! icode return from UM_FORT_FLUSH

REAL  ::   constant 
REAL  ::   h_print     ! needed for printing idealised orography

LOGICAL :: l_update_lbcs, l_apply_lbcs, l_balance, gcr_zero_guess_it

INTEGER :: itemp
INTEGER :: gi

! Oxidant mass-mixing ratios and concentrations, for use in sulphur
! cycle.
REAL, DIMENSION(:,:,:), ALLOCATABLE :: o3_mmr, hno3_mmr, h2o2_mmr
REAL, DIMENSION(:,:,:), ALLOCATABLE :: oh_conc, ho2_conc  

LOGICAL :: l_physics_store

! Local variables for using the aerosol climatology for NWP
     
! Internal array of mass-mixing ratios
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: arcl
      
! Number of requested species within the climatology
INTEGER n_arcl_species

! Corresponding number of requested components
INTEGER n_arcl_compnts

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
REAL, DIMENSION (:,:,:), ALLOCATABLE:: &
      cf_bulk_nohalo, cf_liquid_nohalo, cf_frozen_nohalo


LOGICAL, SAVE :: first_atmstep_call = .TRUE.

END MODULE atm_step_local

