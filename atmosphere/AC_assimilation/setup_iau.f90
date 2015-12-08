! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ IAU set-up routine

MODULE setup_iau_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE Setup_IAU

! Description:
!
!   Basically, set up a data structure for each of the IAU increment files
!   containing the IAU weights etc., and trap any setup errors.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE IAU_mod, ONLY :                 &
    MaxFileNameLen,                 &
    MaxNum_IAU_incs,                &
    IAU_unit,                       &
    IAU_NumFldCodes,                &
    IAU_FldCodes,                   &
    IAU_FldDescs,                   &
    TimeID_Once,                    &
    TimeID_Constant,                &
    CallFreq_Never,                 &
    CallFreq_EveryCallIncsToAdd,    &
    CallFreq_EndInc1Window,         &
    CallFreq_LastIAUCall,           &
    IAU_incs,                       &
    IAU_LocFldLens,                 &
    q_min,                          &
    IAU_FirstCallTS,                &
    IAU_LastCallTS,                 &
    Num_IAU_incs,                   &
    L_IAU_SpecifyInc1Filter,        &
    IAU_FilterType,                 &
    IAU_StartMin,                   &
    IAU_EndMin,                     &
    IAU_ApexMin,                    &
    IAU_Cutoff_period,              &
    IAU_SBE_period,                 &
    L_IAU_IncDiags,                 &
    L_IAU_CalcExnerIncs,            &
    L_IAU_CalcThetaIncs,            &
    L_IAU_CalcRhoIncs,              &
    L_IAU_IncTStar,                 &
    L_IAU_UseSfctempPerts,          &
    L_IAU_UseSoilmPerts,            &
    L_IAU_ResetPoles,               &
    L_IAU_RemoveSS,                 &
    L_IAU_ApplyQTCorrections,       &
    L_IAU_Add_qT_prime_to_q,        &
    L_IAU_IncrementIce,             &
    L_IAU_ScaleCloud,               &
    L_IAU_LimitUpperThetaIncs,      &
    IAU_LimitUpperThetaIncs_pBound, &
    IAU_LimitUpperThetaIncs_maxInc, &
    IAU_QLimitsCallFreq,            &
    L_IAU_QLimitsDiags,             &
    L_IAU_RmNonTropQIncs,           &
    IAU_qsatScale,                  &
    IAU_qsatScale_maxLev,           &
    L_IAU_LimitIncOp,               &
    IAU_qCLThreshScale,             &
    IAU_trop_min_RH,                &
    IAU_nonTrop_max_q,              &
    IAU_nonTrop_min_q,              &
    IAU_nonTrop_max_RH,             &
    IAU_trop_min_p,                 &
    IAU_trop_max_PV,                &
    IAU_nonTrop_max_p,              &
    L_IAU_IgnoreTopLevThetaIncs

USE um_input_control_mod, ONLY:     &
    lcal360,                        &
    model_domain
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE IO, ONLY: &
    setpos,&
    buffin,&
    file_open,&
    file_close
USE model_file
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE UM_ParVars
USE Control_Max_Sizes
USE domain_params
USE dust_parameters_mod, ONLY: l_twobin_dust
USE lookup_addresses
USE turb_diff_mod, ONLY: l_qpos, qlimit

USE calc_tfiltwts_mod, ONLY: calc_tfiltwts
USE Submodel_Mod

USE chsunits_mod, ONLY : nunits

IMPLICIT NONE


! DEPENDS ON: read_flh
! DEPENDS ON: ioerror
! DEPENDS ON: time2sec

! Common blocks:

! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!
! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! This file belongs in section: Top Level

!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:

!------------------   End of Physics   ---------------------------------

      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(model_levels_max) ! Timestep of max w
      INTEGER :: time_div_max(model_levels_max) ! Timestep of max div
      INTEGER :: time_div_min(model_levels_max) ! Timestep of min div
      INTEGER :: time_lapse_min(model_levels_max) ! Timestep of min
      INTEGER :: time_max_shear(model_levels_max) !Timestep max shear
      INTEGER :: time_max_wind(model_levels_max) ! Timestep of max wind
      INTEGER :: time_KE_max(model_levels_max) ! Timestep of max KE
      INTEGER :: time_KE_min(model_levels_max) ! Timestep of min KE
      INTEGER :: time_noise_max(model_levels_max) ! Timestep of max

      REAL:: frictional_timescale(model_levels_max) ! For idealised case
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(0:model_levels_max) ! Max w at a level
      REAL :: max_div_run(model_levels_max) ! Max divergence at a level
      REAL :: min_div_run(model_levels_max) ! Min divergence at a level
      REAL :: min_lapse_run(model_levels_max) ! Min dtheta/dz at a level
      REAL :: max_shear_run(model_levels_max) ! Max shear at a level
      REAL :: max_wind_run(model_levels_max) ! Max wind at a level
      REAL :: max_KE_run(model_levels_max)   ! Max KE at a level
      REAL :: min_KE_run(model_levels_max)   ! Min KE at a level
      REAL :: max_noise_run(model_levels_max) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA
!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
        rpemax, rpemin, ipesum, rpesum,                                 &
        max_w_run, min_theta1_run, dtheta1_run,                         &
        max_div_run, min_div_run, min_lapse_run,                        &
        max_shear_run, max_wind_run,                                    &
        max_noise_run, max_KE_run, min_KE_run,                          &
        time_KE_max, time_KE_min,                                       &
        time_w_max, time_div_max, time_div_min, time_lapse_min,         &
        time_max_shear, time_max_wind,                                  &
        time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(model_levels_max)
      REAL :: friction_level(model_levels_max)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                              &
       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
       SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
       base_frictional_timescale, SuHe_sigma_cutoff,                   &
       L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
       SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: vert_grid_ratio
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(model_levels_max)
      REAL :: rho_ref(model_levels_max)
      REAL :: exner_ref(model_levels_max + 1)
      REAL :: q_ref(model_levels_max)
      REAL :: u_ref(model_levels_max)
      REAL :: v_ref(model_levels_max)
      REAL :: z_orog_print(0:model_levels_max)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: qforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: uforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: vforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind

! ENDGAME
      REAL :: T_surface
      REAL :: Eccentricity
      ! Following two variables used only if L_rotate_grid=.true.
      REAL :: grid_NP_lon ! Longitude (degrees) of grid's north pole
      REAL :: grid_NP_lat ! Latitude (degrees) of grid's north pole
      REAL :: AA_jet_u0   ! See QJRMS 133,1605--1614
      REAL :: AA_jet_A    !
      REAL :: theta_pert
      REAL :: ring_height
      REAL :: angular_velocity ! Planet's angular velocity
      REAL :: T0_P, T0_E ! deep atmosphere baroclinic wave surface temperatures
      INTEGER :: Trefer_number
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: AA_jet_m  ! See QJRMS 133,1605--1614
      INTEGER :: AA_jet_n  !
      INTEGER :: chain_number ! Run continuation number
      LOGICAL :: L_rotate_grid    ! .true. for rotating North pole of grid
      LOGICAL :: L_baro_Perturbed ! Used for baroclinic test to specify
                                  ! pert or steady
      LOGICAL :: L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,  &
                 L_HeldSuarez2_drag,                                        &
                 L_baro_inst, L_isothermal, L_exact_profile, L_balanced,    &
                 L_solid_body
      LOGICAL :: L_deep_baro_inst ! deep atmosphere baroclinic wave switch          


      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type
      INTEGER :: b_const, k_const ! deep atmosphere baroclinic wave parameters

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                              &
       h_o, h_o_actual, h_o_per_step,                                  &
       lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
       Witch_power, plat_size_x, plat_size_y,                          &
       height_domain, delta_x, delta_y, big_factor, mag, vert_grid_ratio, &
       first_theta_height, thin_theta_height, p_surface,               &
       theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
       u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
       ujet_lat, ujet_width,                                           &
       t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
       u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
       z_orog_print, grow_steps,                                       &
       surface_type, grid_number, grid_flat,                           &
       tprofile_number, qprofile_number, uvprofile_number,             &
       num_profile_data, num_uvprofile_data,                           &
       tforce_option, qforce_option, uvforce_option,                   &
       num_tforce_levels, num_tforce_times,                            &
       num_qforce_levels, num_qforce_times,                            &
       num_uvforce_levels, num_uvforce_times,                          &
       L_pforce, pforce_option, num_pforce_times,                      &
       first_constant_r_rho_level_new,                                 &
       big_layers, transit_layers, mod_layers,                         &
       zprofile_data, tprofile_data, qprofile_data,                    &
       z_uvprofile_data, uprofile_data, vprofile_data,                 &
       tforce_time_interval, qforce_time_interval,                     &
       uvforce_time_interval, newtonian_timescale,                     &
       z_tforce_data, z_qforce_data, z_uvforce_data,                   &
       tforce_data, qforce_data, uforce_data, vforce_data,             &
       tforce_data_modlev, qforce_data_modlev,                         &
       uforce_data_modlev, vforce_data_modlev,                         &
       pforce_time_interval, p_surface_data,                           &
       L_initialise_data,                                              &
       L_perturb_t, perturb_magnitude_t,                               &
       L_perturb_q, perturb_magnitude_q,                               &
       L_perturb_correlate_tq,                                         &
       L_perturb_correlate_vert,                                       &
       L_perturb_correlate_time,                                       &
       perturb_type, perturb_height,                                   &
       L_constant_dz, L_polar_wind_zero,                               &
       L_wind_balance, L_rotate_winds,                                 &
       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
       L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
       L_pressure_balance, L_vert_Coriolis,                            &
       cool_rate, L_force, L_force_lbc,                                &
       zprofile_orog, idl_interp_option, hf,                           &
       L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
       L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
       idl_bubble_option, idl_bubble_max                               &
      , idl_bubble_height, idl_bubble_width, idl_bubble_depth          &
      , idl_bubble_xoffset,idl_bubble_yoffset                          &
      , L_idl_bubble_saturate,                                         &
       L_rotating, L_fixed_lbcs, L_code_test,                          &
       L_baroclinic, L_cyclone,                                        &
       L_damp, L_geo_for, L_bomex,                                     &
       DMPTIM, HDMP, ZDMP,                                             &
       u_geo, v_geo,                                                   &
!ENDGAME
       T_surface, chain_number,                                        &
       Trefer_number,                                                  &
       tstep_plot_frequency, tstep_plot_start, Eccentricity,           &
       L_rotate_grid, grid_NP_lon, grid_NP_lat,                        &
       AA_jet_m, AA_jet_n, AA_jet_u0, AA_jet_A, L_baro_Perturbed,      &
       L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,       &
       L_HeldSuarez2_drag,                                             &
       L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const, k_const,    &      
       ring_height, theta_pert, L_isothermal,                          &
       L_exact_profile, L_balanced, L_solid_body, angular_velocity
! CRUNTIMC end
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'Setup_IAU'

! Local variables:

INTEGER :: i
INTEGER :: ICode
INTEGER :: TS_len_secs
INTEGER :: Code
INTEGER :: Level
INTEGER :: LenIO
INTEGER :: LenBase
INTEGER :: FieldNum
INTEGER :: Lookup_len
INTEGER :: IncNum
INTEGER :: LocFldLen
INTEGER :: pack_code
INTEGER :: flds_len
INTEGER :: TS_num
INTEGER :: TimeId
INTEGER :: VAR_StartSec
INTEGER :: VAR_EndSec
INTEGER :: UM_StartSec
INTEGER :: UM_EndSec
INTEGER :: ElapsedDays
INTEGER :: ElapsedSecs
INTEGER :: DTSecs
INTEGER :: DTMins
INTEGER :: VTSecs
INTEGER :: VTMins
INTEGER :: NumWeights
INTEGER :: WeightNum
INTEGER :: FirstNonZeroWeightIndex
INTEGER :: LastNonZeroWeightIndex
INTEGER :: Sec
INTEGER :: IAUDustBinNum
INTEGER :: IAUNumDustBins

INTEGER :: VAR_StartSecSave(MaxNum_IAU_incs)
INTEGER :: VAR_EndSecSave  (MaxNum_IAU_incs)
INTEGER :: UM_StartSecSave (MaxNum_IAU_incs)
INTEGER :: UM_EndSecSave   (MaxNum_IAU_incs)

REAL :: A_IO

REAL, ALLOCATABLE :: Weights(:)

LOGICAL :: contains_new_qT_code
LOGICAL :: contains_old_qT_code
LOGICAL :: AllPacked 
LOGICAL :: Using_q_qCL_or_qCF
LOGICAL :: Using_qT
LOGICAL :: IAUHasDustBin(6)

CHARACTER(LEN=2)   :: IncNumStr
CHARACTER(LEN=9)   :: VAR_StartStr
CHARACTER(LEN=9)   :: VAR_EndStr
CHARACTER(LEN=9)   :: UM_StartStr
CHARACTER(LEN=9)   :: UM_EndStr
CHARACTER(LEN=9)   :: TimeStr
CHARACTER(LEN=9)   :: IAUNumDustBinsStr
CHARACTER(LEN=336) :: CMessage

CHARACTER(LEN=60)  :: UsableFlds(MaxNum_IAU_incs)

CHARACTER(MaxFileNameLen) :: BaseFilePath
CHARACTER(MaxFileNameLen) :: FilePath

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

!-------------------------------------------------------------------------------
! [1]: Check namelist variables.
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('SETUP_IAU',zhook_in,zhook_handle)
IF (Num_IAU_incs > MaxNum_IAU_incs) THEN
  ICode = 1
  WRITE (CMessage,'(A,I6,A,I6,A)')                                          &
    'Number of IAU increment files (', Num_IAU_incs, ') exceeds maximum (', &
    MaxNum_IAU_incs, ')'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

IF (IAU_QLimitsCallFreq /= CallFreq_Never              .AND. &
    IAU_QLimitsCallFreq /= CallFreq_EveryCallIncsToAdd .AND. &
    IAU_QLimitsCallFreq /= CallFreq_EndInc1Window      .AND. &
    IAU_QLimitsCallFreq /= CallFreq_LastIAUCall) THEN
  ICode = 1
  WRITE (CMessage,'(A,I6,A)') &
    'Value of IAU_QLimitsCallFreq (', IAU_QLimitsCallFreq, ') unsupported'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

IF (L_IAU_RmNonTropQIncs .AND. IAU_QLimitsCallFreq == CallFreq_Never) THEN
  ICode = 1
  CMessage = 'L_IAU_RmNonTropQIncs set, but not calling QLimits'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

!-------------------------------------------------------------------------------
! [2]: Set up q_min, depending on whether q_pos routine is being used.
!-------------------------------------------------------------------------------

IF (L_QPOS) THEN
  q_min = qlimit
ELSE
  q_min = 1.0E-8
END IF

!-------------------------------------------------------------------------------
! [3]: Get local lengths of fields that may be read from the increment files.
!-------------------------------------------------------------------------------

DO i = 1, IAU_NumFldCodes

  Code = IAU_FldCodes(i)

  ! Default:
  IAU_LocFldLens(i) = theta_field_size

  ! Exceptions:
  IF (Code == 2) IAU_LocFldLens(i) = u_field_size ! u
  IF (Code == 3) IAU_LocFldLens(i) = v_field_size ! v

  IF (Code == 4 .AND. L_IAU_CalcThetaIncs)          &
    IAU_LocFldLens(i) = 0 ! Don't read theta

  IF (Code == 9 .AND. .NOT. L_IAU_UseSoilmPerts)    &
    IAU_LocFldLens(i) = 0 ! Don't read soilm perts

  IF (Code == 24 .AND. .NOT. L_IAU_UseSfctempPerts) &
    IAU_LocFldLens(i) = 0 ! Don't read sfctemp perts

  IF (Code == 253 .AND. L_IAU_CalcRhoIncs)          &
    IAU_LocFldLens(i) = 0 ! Don't read rho

  IF (Code == 255 .AND. L_IAU_CalcExnerIncs)        &
    IAU_LocFldLens(i) = 0 ! Don't read exner

  IF (Code == 407 .AND. .NOT.L_IAU_CalcExnerIncs)   &
    IAU_LocFldLens(i) = 0 ! Don't read p

END DO

!-------------------------------------------------------------------------------
! [4]: Print IAU details common to all increments.
!-------------------------------------------------------------------------------

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (6,'(A)') ''
  WRITE (6,'(A)') '<><><><><><><><><><><><><><><><><><><><><><><><><>'
  WRITE (6,'(A)') 'Start of Incremental Analysis Update (IAU) details'
  WRITE (6,'(A)') ''
  WRITE (6,'(A,L6)') 'Write increment diagnostics? ',L_IAU_IncDiags
  WRITE (6,'(A,L6)') 'Apply qT corrections?        ',L_IAU_ApplyQTCorrections
  WRITE (6,'(A,L6)') 'Replace q with q+qTprime?      ', &
                                                 L_IAU_Add_qT_prime_to_q
  WRITE (6,'(A,L6)') 'Calculate ice cloud incs?    ',L_IAU_IncrementIce
  WRITE (6,'(A,L6)') 'Scale cloud incs?            ',L_IAU_ScaleCloud
  WRITE (6,'(A,L6)') 'Limit Var_DiagCloud qCL incs?',L_IAU_LimitIncOp
  IF (L_IAU_LimitIncOp) THEN
    WRITE (6,'(A,ES10.4)') &
      ' - Scaling factor for limit to qCL incs: ', IAU_qCLThreshScale
  END IF
  WRITE (6,'(A,L6)') 'Calculate exner incs?        ',L_IAU_CalcExnerIncs
  WRITE (6,'(A,L6)') 'Calculate theta incs?        ',L_IAU_CalcThetaIncs
  WRITE (6,'(A,L6)') 'Calculate rho   incs?        ',L_IAU_CalcRhoIncs
  WRITE (6,'(A,L6)') 'Ignore top-lev theta incs?   ',L_IAU_IgnoreTopLevThetaIncs
  WRITE (6,'(A,L6)') 'Remove supersaturation?      ',L_IAU_RemoveSS
  IF (L_IAU_RemoveSS) THEN
    IF (IAU_qsatScale_maxLev > 0) THEN
      WRITE (6,'(A,ES10.4)') &
      ' - Scaling for saturation limit:              ', IAU_qsatScale
    END IF
    WRITE (6,'(A,I4)') &
      ' - Max level for scaling of saturation limit: ', IAU_qsatScale_maxLev
  END IF
  WRITE (6,'(A,L6)') 'Increment TStar/TSoil?       ', L_IAU_IncTStar
  IF (Model_domain == mt_global) THEN
    WRITE (6,'(A,L6)') 'Reset polar rows?            ', L_IAU_ResetPoles
  END IF
  WRITE (6,'(A,L6)') 'Limit upper-level theta incs?', L_IAU_LimitUpperThetaIncs
  IF (L_IAU_LimitUpperThetaIncs) THEN
    WRITE (6,'(A,ES10.4)') &
      ' - Pressure boundary: ', IAU_LimitUpperThetaIncs_pBound
    WRITE (6,'(A,ES10.4)') &
      ' - Max abs increment: ', IAU_LimitUpperThetaIncs_maxInc
  END IF
  SELECT CASE (IAU_QLimitsCallFreq)
    CASE (CallFreq_Never)
      WRITE (6,'(A)') 'QLimits call frequency:       Never'
    CASE (CallFreq_EveryCallIncsToAdd)
      WRITE (6,'(A)') 'QLimits call frequency:       Every call on which '// &
                                                'there are incs to add'
    CASE (CallFreq_EndInc1Window)
      WRITE (6,'(A)') 'QLimits call frequency:       End of Inc1 time window'
    CASE (CallFreq_LastIAUCall)
      WRITE (6,'(A)') 'QLimits call frequency:       Last IAU call'
  END SELECT
  IF (IAU_QLimitsCallFreq /= CallFreq_Never) THEN
    WRITE (6,'(A,L6)') &
       '- Write diagnostics?                   ',  L_IAU_QLimitsDiags
    WRITE (6,'(A,L6)') &
       '- Remove non-trop q increments?        ',  L_IAU_RmNonTropQIncs
    WRITE (6,'(A,ES10.4)') &
      ' - Lower limit to apply to trop RH:     ', IAU_trop_min_RH
    WRITE (6,'(A,ES10.4)') &
      ' - Upper limit to apply to non-trop q:  ', IAU_nonTrop_max_q
    WRITE (6,'(A,ES10.4)') &
      ' - Lower limit to apply to non-trop q:  ', IAU_nonTrop_min_q
    WRITE (6,'(A,ES10.4)') &
      ' - Upper limit to apply to non-trop RH: ', IAU_nonTrop_max_RH
    WRITE (6,'(A,ES10.4)') &
      ' - Minimum tropospheric pressure (hPa): ', IAU_trop_min_p    * 100.0
    WRITE (6,'(A,ES10.4)') &
      ' - Maximum tropospheric ABS(PV) (PVU):  ', IAU_trop_max_PV   * 1.0E06
    WRITE (6,'(A,ES10.4)') &
      ' - Maximum non-trop     pressure (hPa): ', IAU_nonTrop_max_p * 100.0
  END IF
END IF ! (PrintStatus >= PrStatus_Normal)

!-------------------------------------------------------------------------------
! [5]: Loop over increment files, setting up a data structure for each one.
!-------------------------------------------------------------------------------

ALLOCATE (IAU_incs(Num_IAU_incs))

! Get base filepath for increments:
CALL FORT_GET_ENV ( FT_ENVIRON  (IAU_unit), & ! in
                    LEN_FT_ENVIR(IAU_unit), & ! in
                    BaseFilePath,           & ! out
                    MaxFileNameLen,         & ! in
                    ICode )                   ! out

IF (ICode > 0) THEN
  CMessage = 'Error reading env variable '// FT_ENVIRON(IAU_unit)
  CALL EReport (RoutineName, ICode, CMessage)
END IF

LenBase = LEN_TRIM(BaseFilePath)

! Timestep length in seconds:
TS_len_secs = SECS_PER_STEPIM(A_IM)

! Loop over increment files:
DO IncNum = 1, Num_IAU_incs

  ! Get increment number as a string:
  WRITE (IncNumStr,'(I2.2)') IncNum

  ! Initialise flags:
  IAU_incs(IncNum) % FieldsToAdd           = .FALSE.
  IAU_incs(IncNum) % Contains_q_qCL_or_qCF = .FALSE.
  IAU_incs(IncNum) % Contains_qT           = .FALSE.
  IAU_incs(IncNum) % FieldsStored          = .FALSE.

  !-----------------------------------------------------------------------------
  ! [5.1]: Read fixed header and lookup headers.
  !-----------------------------------------------------------------------------

  ! Construct filepath, matching the method used in VarProg_AnalysePF:
  IF (IncNum == 1) THEN
    FilePath = BaseFilePath
  ELSE
    i = MIN(LenBase, MaxFileNameLen - 2)
    ! Append increment number:
    WRITE (FilePath(i+1:i+2),'(I2.2)') IncNum
  END IF
  IAU_incs(IncNum) % FilePath =  FilePath

  ! Open file:
  CALL MODEL_FILE_OPEN (IAU_unit, FilePath, MaxFileNameLen, 0, 1, ICode)
  IF (ICode > 0) THEN
    CMessage(1:80) = 'Error opening IAU increment no.'// IncNumStr //':'
    CMessage(81:)  = TRIM(FilePath)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Read in fixed-length header:
  CALL READ_FLH (IAU_unit, IAU_incs(IncNum) % FixHd, Len_FixHd, &
                 ICode, CMessage(161:))
  IF (ICode > 0) THEN
    CMessage(1:80)   = 'Error reading fixed header from increment no.'// &
                       IncNumStr
    CMessage(81:160) = 'Message from READ_FLD:'
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  IAU_incs(IncNum) % Len1Lookup = IAU_incs(IncNum) % FixHd(151)
  IAU_incs(IncNum) % Len2Lookup = IAU_incs(IncNum) % FixHd(152)

  ! Read in lookup tables:
  IF (IAU_incs(IncNum) % Len2Lookup /= 0) THEN

    ALLOCATE (IAU_incs(IncNum) % Lookup(IAU_incs(IncNum) % Len1Lookup,  &
                                        IAU_incs(IncNum) % Len2Lookup))

    CALL SETPOS (IAU_unit, IAU_incs(IncNum) % FixHd(150)-1, ICode)
    IF (ICode > 0) THEN
      CMessage = 'SETPOS error moving to start of lookup tables for '// &
                 'increment no.'// IncNumStr
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    Lookup_len = IAU_incs(IncNum) % Len1Lookup &
               * IAU_incs(IncNum) % Len2Lookup

    CALL BUFFIN ( IAU_unit,                       & 
                  IAU_incs(IncNum) % Lookup, & 
                  Lookup_len,                     & 
                  LenIO,                          & 
                  A_IO )                             
    IF ( (A_IO  /= -1.0      ) .OR.   &
         (LenIO /= Lookup_len) ) THEN
      CALL IOERROR ('Buffer in of lookups from IAU increment', &
                    A_IO, LenIO, Lookup_len)
      ICode    = NINT(A_IO) + 1
      CMessage = 'Error reading lookup tables from increment no.'// IncNumStr
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

  END IF

  ! Close file:
  CALL MODEL_FILE_CLOSE (IAU_unit, FilePath, MaxFileNameLen, 1, 0, ICode)
  IF (ICode > 0) THEN
    CMessage(1:80) = 'Error closing IAU increment no.'// IncNumStr //':'
    CMessage(81:)  = TRIM(FilePath)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Cycle if there are no fields in the increment:
  IF (IAU_incs(IncNum) % Len2Lookup == 0) CYCLE

  !-----------------------------------------------------------------------------
  ! [5.2]: Get length of array required to hold data, and determine whether
  !        the data is all 32-bit packed, in which case the data will be
  !        stored in packed form.
  !-----------------------------------------------------------------------------

  AllPacked = .TRUE.
  flds_len  = 0

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup

    Code = IAU_incs(IncNum) % Lookup(ITEM_CODE, FieldNum)

    ! Get local field length:
    LocFldLen = 0
    DO i = 1, IAU_NumFldCodes
      IF (IAU_FldCodes(i) == Code) LocFldLen = IAU_LocFldLens(i)
    END DO

    flds_len = flds_len + LocFldLen

    pack_code = MOD(IAU_incs(IncNum) % Lookup(LBPACK, FieldNum), 10)
    IF (pack_code /= 2) AllPacked = .FALSE.

  END DO

  ! Cycle if there are no relevant fields in the file:
  IF (flds_len == 0) THEN
    DEALLOCATE (IAU_incs(IncNum) % Lookup)
    CYCLE
  ELSE
    IAU_incs(IncNum) % FieldsToAdd = .TRUE.
  END IF

  IAU_incs(IncNum) % StorePacked = AllPacked
  IAU_incs(IncNum) % flds_len    = flds_len

  !-----------------------------------------------------------------------------
  ! [5.3]: See whether there are increments to qT, q, qCL, qCF, dust1-6.
  !-----------------------------------------------------------------------------

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code  = IAU_incs(IncNum) % Lookup(ITEM_CODE, FieldNum)
    Level = IAU_incs(IncNum) % Lookup(LBLEV,     FieldNum)
    IF (Level <= wet_levels) THEN
      IF      (Code == 10 .OR. Code == 12 .OR. Code == 254) THEN
        IAU_incs(IncNum) % Contains_q_qCL_or_qCF = .TRUE.
      ELSE IF (Code == 16207 .OR. Code == 18001) THEN
        IAU_incs(IncNum) % Contains_qT           = .TRUE.
      END IF
    END IF
  END DO

  ! Abort if file contains qT fields with both the new and old STASH codes:
  contains_new_qT_code = .FALSE.
  contains_old_qT_code = .FALSE.
  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code  = IAU_incs(IncNum) % Lookup(ITEM_CODE, FieldNum)
    Level = IAU_incs(IncNum) % Lookup(LBLEV,     FieldNum)
    IF (Level <= wet_levels) THEN
      IF (Code == 16207) contains_new_qT_code = .TRUE.
      IF (Code == 18001) contains_old_qT_code = .TRUE.
    END IF
  END DO
  IF (contains_new_qT_code .AND. contains_old_qT_code) THEN
    ICode = 1
    CMessage = 'Old qT code (18001) being used with new one (16207) in '// &
               'increment no.'// IncNumStr
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Abort if the file contains the wrong number of dust bins for the model
  ! being run (2 bin or 6 bin)
  IAUNumDustBins = 0
  DO i = 1, 6
    IAUHasDustBin(i) = .FALSE.
  END DO
  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code = IAU_incs(IncNum) % Lookup(ITEM_CODE, FieldNum)
    ! All dust bins have codes in this range:
    IF (Code >= 431 .AND. Code <= 436) THEN
      IAUDustBinNum = Code-430
      IF (.NOT. IAUHasDustBin(IAUDustBinNum)) THEN
        IAUHasDustBin(IAUDustBinNum) = .TRUE.
        IAUNumDustBins = IAUNumDustBins + 1
      END IF
    END IF
  END DO
  WRITE (IAUNumDustBinsStr,'(I9)') IAUNumDustBins
  IF (IAUNumDustBins > 0) THEN
    IF (l_twobin_dust .AND. IAUNumDustBins /= 2) THEN
      ICode = 1
      CMessage = 'Number of dust bins in IAU is ' // IAUNumDustBinsStr  // &
                 ', but the dust scheme is being run with 2 size bins.'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
    IF (.NOT. l_twobin_dust .AND. IAUNumDustBins /= 6) THEN
      ICode = 1
      CMessage = 'Number of dust bins in IAU is ' // IAUNumDustBinsStr  // &
                 ', but the dust scheme is being run with 6 size bins.'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
  END IF

  !-----------------------------------------------------------------------------
  ! [5.4]: Get list of usable fields.
  !-----------------------------------------------------------------------------

  IF (PrintStatus >= PrStatus_Normal) THEN
    UsableFlds(IncNum) = ''
    DO i = 1, IAU_NumFldCodes
      IF (ANY(IAU_incs(IncNum) % Lookup(ITEM_CODE,:) == IAU_FldCodes(i))) &
        UsableFlds(IncNum) = TRIM(UsableFlds(IncNum)) //','//             &
                             TRIM(IAU_FldDescs(i))
    END DO
    UsableFlds(IncNum) = UsableFlds(IncNum)(2:)
  END IF

  !-----------------------------------------------------------------------------
  ! [5.5]: Calculate filter details.
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [5.5.1]: Get increment data and validity times relative to basis time.
  !-----------------------------------------------------------------------------

  CALL TIME2SEC ( IAU_incs(IncNum) % FixHd(21), & ! in
                  IAU_incs(IncNum) % FixHd(22), & ! in
                  IAU_incs(IncNum) % FixHd(23), & ! in
                  IAU_incs(IncNum) % FixHd(24), & ! in
                  IAU_incs(IncNum) % FixHd(25), & ! in
                  0,                            & ! in
                  BASIS_TIME_DAYS,              & ! in
                  BASIS_TIME_SECS,              & ! in
                  ElapsedDays,                  & ! out
                  ElapsedSecs,                  & ! out
                  LCAL360 )                       ! in

  DTSecs = ElapsedSecs + ElapsedDays * 86400
  DTMins = DTSecs / 60

  CALL TIME2SEC ( IAU_incs(IncNum) % FixHd(28), & ! in
                  IAU_incs(IncNum) % FixHd(29), & ! in
                  IAU_incs(IncNum) % FixHd(30), & ! in
                  IAU_incs(IncNum) % FixHd(31), & ! in
                  IAU_incs(IncNum) % FixHd(32), & ! in
                  0,                            & ! in
                  BASIS_TIME_DAYS,              & ! in
                  BASIS_TIME_SECS,              & ! in
                  ElapsedDays,                  & ! out
                  ElapsedSecs,                  & ! out
                  LCAL360 )                       ! in

  VTSecs = ElapsedSecs + ElapsedDays * 86400
  VTMins = VTSecs / 60

  IF (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter) THEN

    !---------------------------------------------------------------------------
    ! [5.5.2]: Filter details are specified via the IAU namelist.
    !---------------------------------------------------------------------------

    IAU_incs(1) % TendencyIncs = .FALSE.

    IF ( (VTMins < IAU_StartMin) .OR. &
         (VTMins > IAU_EndMin) ) THEN
      ICode = 1
      CMessage = 'Validity time of first IAU increment outside IAU period'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    !---------------------------------------------------------------------------
    ! [5.5.2.1]: Calculate weights covering the whole filter span.
    !---------------------------------------------------------------------------

    NumWeights = (IAU_EndMin - IAU_StartMin)*60 / TS_len_secs + 1

    ALLOCATE (Weights(NumWeights))

    CALL Calc_TFiltWts ( IAU_StartMin*60,    & ! in
                         IAU_EndMin*60,      & ! in
                         TS_len_secs,        & ! in
                         IAU_ApexMin*60,     & ! in
                         IAU_Cutoff_period,  & ! in
                         IAU_SBE_period,     & ! in
                         IAU_FilterType,     & ! in
                         NumWeights,         & ! in
                         Weights )             ! out

    !---------------------------------------------------------------------------
    ! [5.5.2.2]: Eliminate zero weights at the start and end of the span.
    !---------------------------------------------------------------------------

    FirstNonZeroWeightIndex = 0
    LastNonZeroWeightIndex  = 0

    DO i = 1, NumWeights
      IF (Weights(i) /= 0.0) THEN
        IF (FirstNonZeroWeightIndex == 0) FirstNonZeroWeightIndex = i
        LastNonZeroWeightIndex = i
      END IF
    END DO

    NumWeights = LastNonZeroWeightIndex - FirstNonZeroWeightIndex + 1

    ALLOCATE (IAU_incs(1) % Weights(NumWeights))

    IAU_incs(1) % Weights(1:NumWeights) = Weights(FirstNonZeroWeightIndex: &
                                                  LastNonZeroWeightIndex)

    DEALLOCATE (Weights)

    IAU_incs(1) % InsertionStartTS = IAU_StartMin * 60 / TS_len_secs &
                                   + FirstNonZeroWeightIndex - 1
    IAU_incs(1) % InsertionEndTS   = IAU_StartMin * 60 / TS_len_secs &
                                   + LastNonZeroWeightIndex  - 1

  ELSE

    !---------------------------------------------------------------------------
    ! [5.5.3]: Filter details to be deduced from information written to the
    !          fixed header by VAR.
    !
    ! The key parameter is the "TimeID", supplied as FixHd(10).
    !
    ! There are two main possibilities:
    !
    ! 1. TimeID = TimeID_Once
    !    Add the whole increment at the UM timestep nearest to the validity
    !    time supplied in the fixed header.
    !
    ! 2. TimeID /= TimeID_Once
    !    Treat the increments as (per second) tendencies for a time window
    !    starting at the validity time and ending at the data time supplied in
    !    the fixed header; the period over which the tendency will have been
    !    applied within VAR.
    !
    !    The start and end points of this window are adjusted so that they
    !    coincide with UM timesteps, and a correction factor calculated to
    !    account for the change in window length.
    !
    !    The shape of the weighting function to be used over the window is
    !    determined by the specific value of TimeID. Currently the only choice
    !    is uniform weights, but further choices may be added in the future.
    !---------------------------------------------------------------------------

    ! TimeID supplied by VAR:
    TimeID = IAU_incs(IncNum) % FixHd(10)

    IF (TimeID /= TimeID_Once .AND. &
        TimeID /= TimeID_Constant) THEN
      ICode = 1
      WRITE (CMessage,'(A,I6,A)') &
        'Unsupported TimeID (', TimeID, ') for increment no.'// IncNumStr
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    VAR_StartSec = VTSecs
    VAR_EndSec   = DTSecs

    IF (VAR_EndSec < VAR_StartSec) THEN
      ICode = 1
      CMessage = 'VAR end time before start time for increment no.'// &
                 IncNumStr
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    ! Adjust start and end times to coincide with UM timesteps:
    UM_StartSec = TS_len_secs * (VAR_StartSec / TS_len_secs)
    UM_EndSec   = TS_len_secs * (VAR_EndSec   / TS_len_secs)
    IF ((VAR_StartSec - UM_StartSec) * 2 > TS_len_secs) THEN
      UM_StartSec = UM_StartSec + TS_len_secs
    END IF
    IF ((VAR_EndSec   - UM_EndSec)   * 2 > TS_len_secs) THEN
      UM_EndSec   = UM_EndSec   + TS_len_secs
    END IF

    ! Save start and end times for printout step:
    IF (PrintStatus >= PrStatus_Normal) THEN
      VAR_StartSecSave(IncNum) = VAR_StartSec
      VAR_EndSecSave  (IncNum) = VAR_EndSec
      UM_StartSecSave (IncNum) = UM_StartSec
      UM_EndSecSave   (IncNum) = UM_EndSec
    END IF

    IF (TimeID == TimeID_Once) THEN

      IAU_incs(IncNum) % TendencyIncs = .FALSE.

      IF (VAR_StartSec /= VAR_EndSec) THEN
        ICode = 1
        CMessage = 'TimeID=TimeID_Once, but VT and DT do not match for '// &
                   'increment no.'// IncNumStr
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      ALLOCATE (IAU_incs(IncNum) % Weights(1))

      IAU_incs(IncNum) % Weights(1)       = 1.0
      IAU_incs(IncNum) % InsertionStartTS = UM_StartSec / TS_len_secs
      IAU_incs(IncNum) % InsertionEndTS   = UM_EndSec   / TS_len_secs

    ELSE

      IAU_incs(IncNum) % TendencyIncs = .TRUE.

      IF (UM_StartSec == UM_EndSec) THEN
        ICode = 1
        CMessage = 'TimeID/=TimeID_Once, but zero period in UM for '// &
                   'increment no.'// IncNumStr
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      ! Tendencies applied at the end of each timestep within the window:
      IAU_incs(IncNum) % InsertionStartTS = (UM_StartSec / TS_len_secs) + 1
      IAU_incs(IncNum) % InsertionEndTS   =  UM_EndSec   / TS_len_secs

      NumWeights = IAU_incs(IncNum) % InsertionEndTS   &
                 - IAU_incs(IncNum) % InsertionStartTS &
                 + 1

      ALLOCATE (IAU_incs(IncNum) % Weights(NumWeights))

      ! Calculate filter weights:
      SELECT CASE (TimeID)
        CASE (TimeID_Constant)
          IAU_incs(IncNum) % Weights(:) = 1.0
      END SELECT

      ! Scale weights to convert from (per second) tendencies to
      ! instantaneous increments. (This includes a correction factor to
      ! account for the change of window moving from VAR to the UM.)
      IAU_incs(IncNum) % Weights(:) = IAU_incs(IncNum) % Weights(:)   &
                                    * REAL(TS_len_secs)               &
                                    * REAL(VAR_EndSec - VAR_StartSec) &
                                    / REAL(UM_EndSec  - UM_StartSec)

    END IF ! (TimeID == TimeID_Once)

  END IF ! (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter)

  !-----------------------------------------------------------------------------
  ! [5.5.4]: Check that insertion period does not start before basis time.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % InsertionStartTS < 0) THEN
    ICode = 1
    CMessage = 'Insertion period starts before basis time basis time for '// &
               'increment no.'// IncNumStr
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  !-----------------------------------------------------------------------------
  ! [5.6]: Set StoreFields flag, indicating whether the increments should be
  !        stored in memory during the increment's insertion period.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % InsertionEndTS   > &
      IAU_incs(IncNum) % InsertionStartTS) THEN
    IAU_incs(IncNum) % StoreFields = .TRUE.
  ELSE
    IAU_incs(IncNum) % StoreFields = .FALSE.
  END IF

  !-----------------------------------------------------------------------------
  ! [5.7]: Update timestep numbers for first and last IAU call.
  !-----------------------------------------------------------------------------

  IF (IAU_FirstCallTS == -1) THEN
    IAU_FirstCallTS = IAU_incs(IncNum) % InsertionStartTS
    IAU_LastCallTS  = IAU_incs(IncNum) % InsertionEndTS
  ELSE
    IAU_FirstCallTS = MIN(IAU_FirstCallTS, &
                          IAU_incs(IncNum) % InsertionStartTS)
    IAU_LastCallTS  = MAX(IAU_LastCallTS,  &
                          IAU_incs(IncNum) % InsertionEndTS)
  END IF

END DO ! (IncNum)

!-------------------------------------------------------------------------------
! [6]: Apply some restrictions to the way humidity and cloud are updated.
!
!      There are currently two methods for obtaining humidity and cloud
!      increments:
!        1. Diagnose them from qT increments in the IAU files.
!        2. Update them directly from the q, qCL and/or qCF increments in
!           the IAU files, and then apply further modifications if using the
!           PC2 cloud scheme.
!      For simplicity, we currently demand that only one of these methods is
!      used on an individual IAU timestep.
!-------------------------------------------------------------------------------

! Check that there are no increment files containing q, qCL or qCF if they
! also contain qT:
DO IncNum = 1, Num_IAU_incs

  IF (IAU_incs(IncNum) % Contains_q_qCL_or_qCF .AND. &
      IAU_incs(IncNum) % Contains_qT) THEN
    ICode = 1
    WRITE (CMessage(1:80),'(A,I2.2,A)') &
      'Increment no.', IncNum, ' contains q, qCL or qCF, but also qT.'
    CMessage(81:) = 'This combination is currently unsupported.'
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

END DO

! Check that q, qCL or qCF increments will not be used together with qT
! increments on the same IAU timestep:
DO TS_num = IAU_FirstCallTS, IAU_LastCallTS

  Using_q_qCL_or_qCF = .FALSE.
  Using_qT           = .FALSE.

  DO IncNum = 1, Num_IAU_incs
    IF (IAU_incs(IncNum) % FieldsToAdd) THEN
      IF (TS_num >= IAU_incs(IncNum) % InsertionStartTS .AND. &
          TS_num <= IAU_incs(IncNum) % InsertionEndTS) THEN
        IF (IAU_incs(IncNum) % Contains_q_qCL_or_qCF) &
          Using_q_qCL_or_qCF = .TRUE.
        IF (IAU_incs(IncNum) % Contains_qT)           &
          Using_qT           = .TRUE.
      END IF
    END IF
  END DO

  IF (Using_q_qCL_or_qCF .AND. Using_qT) THEN
    ICode = 1
    WRITE(CMessage(1:80),'(A,I6)')                                        &
      'qT incs being used together with q, qCL or qCF incs on timestep ', &
      TS_num
    CMessage(81:) = 'This combination is currently unsupported.'
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

END DO

! Print warning if asking to diagnose ice cloud increments from qT
! increments, but no qT increments available:
IF (L_IAU_IncrementIce .AND. .NOT.ANY(IAU_incs(:) % Contains_qT)) THEN
  ICode = -1
  CMessage = 'L_IAU_IncrementIce set, but no qT increments available'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

!-------------------------------------------------------------------------------
! [7]: Print details for each increment.
!-------------------------------------------------------------------------------

IF (PrintStatus >= PrStatus_Normal) THEN

  WRITE (6,'(A)') ''
  WRITE (6,'(A)') '(Times given below are relative to the model basis time.)'
  WRITE (6,'(A)') ''

  DO IncNum = 1, Num_IAU_incs

    ! Get increment number as a string:
    WRITE (IncNumStr,'(I2.2)') IncNum

    WRITE (6,'(A)') '** Details for increment no.'// IncNumStr //':'

    IF (.NOT.IAU_incs(IncNum) % FieldsToAdd) THEN
      WRITE (6,'(A)') '   No usable fields'
      CYCLE
    ELSE
      WRITE (6,'(A)') '   Usable fields: '// TRIM(UsableFlds(IncNum))
    END IF

    IF (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter) THEN

      IF (IAU_FilterType == 1) THEN
        WRITE (6,'(A)') ' Filter type: Uniform'
      ELSE IF (IAU_FilterType == 2) THEN
        WRITE (6,'(A)') ' Filter type: Triangular'
        WRITE (6,'(A,I8)') ' Apex minute:           ', IAU_ApexMin
      ELSE IF (IAU_FilterType == 3) THEN
        WRITE (6,'(A)') ' Filter type: Lanczos windowed'
        WRITE (6,'(A,ES10.4,A)') ' Cut-off period: ', IAU_Cutoff_period, ' hrs'
      ELSE IF (IAU_FilterType == 4) THEN
        WRITE (6,'(A)') ' Filter type: Dolph'
        WRITE (6,'(A,ES10.4,A)')          &
                             ' Stop band edge period: ', IAU_SBE_period, ' hrs'
        WRITE (6,'(A,I8)') ' Insertion start min:   ', IAU_StartMin
        WRITE (6,'(A,I8)') ' Insertion end   min:   ', IAU_EndMin
      END IF

    ELSE

      VAR_StartSec = VAR_StartSecSave(IncNum)
      VAR_EndSec   = VAR_EndSecSave  (IncNum)
      UM_StartSec  = UM_StartSecSave (IncNum)
      UM_EndSec    = UM_EndSecSave   (IncNum)

      ! Time strings for start and end of time windows:
      WRITE (VAR_StartStr,'(I3,":",I2.2,":",I2.2)') &
        VAR_StartSec/3600, MOD(VAR_StartSec/60, 60), MOD(VAR_StartSec, 60) 
      WRITE (VAR_EndStr,  '(I3,":",I2.2,":",I2.2)') &
        VAR_EndSec  /3600, MOD(VAR_EndSec  /60, 60), MOD(VAR_EndSec,   60) 
      WRITE (UM_StartStr, '(I3,":",I2.2,":",I2.2)') &
        UM_StartSec /3600, MOD(UM_StartSec /60, 60), MOD(UM_StartSec,  60) 
      WRITE (UM_EndStr,   '(I3,":",I2.2,":",I2.2)') &
        UM_EndSec   /3600, MOD(UM_EndSec   /60, 60), MOD(UM_EndSec,    60) 

      WRITE (6,'(A,I8)') '   TimeID:       ', IAU_incs(IncNum) % FixHd(10)
      WRITE (6,'(A,I8)') '   Tendency?:    ', IAU_incs(IncNum) % TendencyIncs
      IF (IAU_incs(IncNum) % TendencyIncs) THEN
        WRITE (6,'(A)') '   VAR window (hr:mm:ss): ', &
                    VAR_StartStr, ' to ', VAR_EndStr
        WRITE (6,'(A)') '   UM  window (hr:mm:ss): ', &
                    UM_StartStr,  ' to ', UM_EndStr
      ELSE
        WRITE (6,'(A)') '   Increment VT   (hr:mm:ss): ', VAR_StartStr
        WRITE (6,'(A)') '   Insertion time (hr:mm:ss): ', UM_StartStr
      END IF

    END IF

    NumWeights = SIZE(IAU_incs(IncNum) % Weights)

    IF (NumWeights > 1) THEN
      WRITE (6,'(A)') '   IAU Weights:  hr:mm:ss  Weight'
      DO WeightNum = 1, NumWeights
        Sec = (IAU_incs(IncNum) % InsertionStartTS + WeightNum - 1) &
            * TS_len_secs
        WRITE (TimeStr,'(I3,":",I2.2,":",I2.2)')  &
          Sec/3600, MOD(Sec/60, 60), MOD(Sec, 60)
        WRITE (6,'(A,A,ES18.10)') &
          '                 ', TimeStr, IAU_incs(IncNum) % Weights(WeightNum)
      END DO
    END IF

  END DO ! (IncNum)

  WRITE (6,'(A)') ''
  WRITE (6,'(A)') 'End of Incremental Analysis Update (IAU) details'
  WRITE (6,'(A)') '<><><><><><><><><><><><><><><><><><><><><><><><>'
  WRITE (6,'(A)') ''

END IF ! (PrintStatus >= PrStatus_Normal)
IF (lhook) CALL dr_hook('SETUP_IAU',zhook_out,zhook_handle)
RETURN


END SUBROUTINE Setup_IAU

END MODULE setup_iau_mod
