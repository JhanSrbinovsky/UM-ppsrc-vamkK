! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Main routine for IAU scheme

SUBROUTINE IAU (                                  &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                 l_mixing_ratio,                  & ! in
                 u,      v,          w,           & ! inout
                 u_adv,  v_adv,      w_adv,       & ! inout
                 theta,  rho,        murk,        & ! inout
                 q,      qCL,        qCF,         & ! inout
                 TStar,  TStar_tile, TSoil, smcl, & ! inout
                 tstar_anom,                      & ! inout
                 p_star, p,                       & ! inout
                 p_theta_levels,                  & ! inout
                 exner,                           & ! inout
                 exner_theta_levels,              & ! inout
                 snow_depth,                      & ! inout
                 area_cloud_fraction,             & ! inout
                 bulk_cloud_fraction,             & ! inout
                 cloud_fraction_liquid,           & ! inout
                 cloud_fraction_frozen,           & ! inout
                 dust_div1, dust_div2, dust_div3, & ! inout
                 dust_div4, dust_div5, dust_div6, & ! inout
                 ozone_tracer )                     ! inout

! Description:
!
!   Main routine for the Incremental Analysis Update (IAU) scheme, as described
!   in UMDP 31.
!
!   Model fields can be updated in two ways:
!
!     1. Directly, via corresponding increment fields in a series of IAU
!        increment files.
!     2. Indirectly, via the increments to other variables.
!
!   Currently, direct incrementing is supported for the following fields:
!
!     u, v, w, theta, rho, exner, p, q, qCL, qCF, murk (aerosol), ozone,
!     sfctemp, soilm, dust_div1-6
!
!   theta, rho and exner may also be updated indirectly via the following
!   relationships:
!
!     theta - hydrostatic equation
!     rho   - equation of state
!     exner - expression in terms of p
!
!   These three options are activated via namelist switches. When one of these
!   is activated, direct updates to the corresponding field are switched off.
!
!   As well as the above set of fields, the IAU increment files may also
!   contain increments to the total specific humidity qT. When qT increments
!   are available, there are two options: (1) use the routine Var_DiagCloud to
!   diagnose increments to q, qCL and (optionally) qCF, or (2) add qT' to q. 
!   Direct and indirect increments to q, qCL or qCF are not allowed on
!   the same timestep.
!
!   When the PC2 cloud scheme is in use, direct updates to q, qCL or qCF
!   are followed up by adjustments to q, qCL and theta to ensure compatibilty
!   with the scheme.
!
!   There is also an option to pass the level-one temperature increments on to
!   the surface temperature fields (TStar and TStar_tile) at land points, and
!   the top-level soil temperature (TSoil(1)) field.
!
!   Sfctemp and soilm are used as initial perturbations in ensemble forecasts.
!   Sfctemp is passed to TStar (at all points) and to TStar_tile and TSoil(1)
!   over land. Soilm is passed to all levels of smcl.
!   These fields are automatically processed if present in the increment file.
!
!   Increments to u, v, and w are always passed on to their advected
!   counterparts u_adv, v_adv and w_adv.
!   
!   At the end of the routine two options are available to limit humidity and
!   cloud. The first option is to remove supersaturation wrt water. The second
!   is to apply limits to tropospheric and non-troposphere humidities via the
!   routine QLimits, which also removes cloud outside the troposphere. The
!   frequency with which QLimits is called is controlled via a namelist
!   variable. For example, it can be called on every IAU call for which
!   increments are to be added, or only at the end of the IAU insertion period.
!
!   Finally, in order to minimise problems with the spin-down in precipitation
!   following the addition of an analysis increment, options are available to
!   modify the output of the cloud incrementing operator (Var_DiagCloud).
!   Also, the removal of supersaturation noted above can be done with respect
!   to a scaled value of qsat.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE Submodel_Mod

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod, ONLY :   &
    pdims, pdims_s,                 &
    qdims,          qdims_l,        &
    tdims, tdims_s,                 &
    udims, udims_s, udims_l,        &
    vdims, vdims_s, vdims_l,        &
    wdims, wdims_s, wdims_l

USE atmos_constants_mod, ONLY :     &
    kappa,                          &
    cp,                             &
    c_virtual,                      &
    p_zero,                         &
    r,                              &
    recip_kappa

USE Control_Max_Sizes, ONLY :       &
    max_121_rows,                   &
    max_look,                       &
    max_n_intf_a,                   &
    max_req_thpv_levs,              &
    max_sponge_width,               &
    max_updiff_levels

USE conversions_mod, ONLY :         &
    pi

USE dyn_coriolis_mod, ONLY :        &
    f3_at_v

USE earth_constants_mod, ONLY :     &
    g

USE ereport_mod, ONLY :             &
    ereport

USE IAU_mod, ONLY :                 &
    MaxFileNameLen,                 &
    IAU_unit,                       &
    IAU_NumFldCodes,                &
    IAU_FldCodes,                   &
    CallFreq_Never,                 &
    CallFreq_EveryCallIncsToAdd,    &
    CallFreq_EndInc1Window,         &
    CallFreq_LastIAUCall,           &
    IAU_incs,                       &
    IAU_LocFldLens,                 &
    q_min,                          &
    IAU_LastCallTS,                 &
    Num_IAU_incs,                   &
    L_IAU_CalcExnerIncs,            &
    L_IAU_CalcThetaIncs,            &
    L_IAU_CalcRhoIncs,              &
    L_IAU_IncTStar,                 &
    L_IAU_IncTStar_tile,            &
    L_IAU_IncTSurf_Snow,            &
    L_IAU_IncTLake,                 &
    L_IAU_RemoveSS,                 &
    L_IAU_ApplyQTCorrections,       &
    L_IAU_Add_qT_prime_to_q,        &
    L_IAU_IncrementIce,             &
    L_IAU_ScaleCloud,               &
    L_IAU_LimitUpperThetaIncs,      &
    L_IAU_SetOzoneMin,              &
    IAU_LimitUpperThetaIncs_pBound, &
    IAU_LimitUpperThetaIncs_maxInc, &
    IAU_QLimitsCallFreq,            &
    L_IAU_QLimitsDiags,             &
    L_IAU_RmNonTropQIncs,           &
    IAU_trop_min_RH,                &
    IAU_nonTrop_max_q,              &
    IAU_nonTrop_min_q,              &
    IAU_nonTrop_max_RH,             &
    IAU_trop_min_p,                 &
    IAU_trop_max_PV,                &
    IAU_nonTrop_max_p,              &
    L_IAU_IgnoreTopLevThetaIncs,    &
    IAU_qsatScale,                  &
    IAU_qsatScale_maxLev,           &
    L_IAU_LimitIncOp,               &
    IAU_qCLThreshScale

USE IO, ONLY :                      &
    File_Close,                     &
    File_Open

USE model_file

USE level_heights_mod, ONLY :       &
    r_theta_levels,                 &
    r_rho_levels,                   &
    r_at_u,                         &
    r_at_v

USE cloud_inputs_mod, ONLY :        &
    rhcrit,                         &
    l_cld_area,                     &
    l_pc2                           


USE ancilcta_namelist_mod, ONLY:    &
    l_sstanom

USE nstypes, ONLY :                 &
    lake

USE parkind1, ONLY :                &
    jprb,                           &
    jpim

USE PrintStatus_mod, ONLY :         &
    PrintStatus,                    &
    PrStatus_Normal

USE swapable_field_mod, ONLY :      &
    swapable_field_pointer_type

USE trignometric_mod, ONLY :        &
    sec_v_latitude,                 &
    tan_v_latitude,                 &
    sec_theta_latitude

USE UM_ParVars, ONLY :              &
    fld_type_p,                     &
    fld_type_u,                     &
    fld_type_v,                     &
    model_levels_max,               &
    halo_i,                         &
    halo_j,                         &
    offx,                           &
    offy

USE visbty_constants_mod, ONLY :    &
    aero0,                          &
    aeromax

USE water_constants_mod, ONLY :     &
    lc,                             &
    tm

USE yomhook, ONLY:                  &
    lhook,                          &
    dr_hook

USE lookup_addresses

USE qlimits_mod, ONLY: qlimits
USE readiaufield_mod, ONLY: readiaufield
USE var_diagcloud_mod, ONLY: var_diagcloud
USE um_input_control_mod,  ONLY:    &
    model_domain
IMPLICIT NONE

! DEPENDS ON: qsat_wat
! DEPENDS ON: expand21
! DEPENDS ON: pack21
! DEPENDS ON: calc_p_star
! DEPENDS ON: calc_exner_at_theta
! DEPENDS ON: calc_p_from_exner
! DEPENDS ON: pc2_assim
! DEPENDS ON: swap_bounds_mv
! DEPENDS ON: calc_pv_at_theta

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
!-----------------Start of TYPCONA------------------------------------

! Constants for routines independent of resolution.
! Requires Use Control_Max_Sizes for MAX_REQ_THPV_LEVS in cconsts.h
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end

! typcona.h originally contained constants for the atmosphere.
! Almost all of these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, RAD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block did not correspond to constants specified in 
! argcona.h, so it remains here, even though argcona.h has been deleted.
! cderv_trig and CDERIVED in cconsts.h should be moved to modules too.

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end
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

! Subroutine arguments:

LOGICAL, INTENT(IN) :: l_mixing_ratio ! Use mixing ratio code?

REAL,    INTENT(INOUT), TARGET ::           &
  u      ( udims_s%i_start : udims_s%i_end, &
           udims_s%j_start : udims_s%j_end, &
           udims_s%k_start : udims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  v      ( vdims_s%i_start : vdims_s%i_end, &
           vdims_s%j_start : vdims_s%j_end, &
           vdims_s%k_start : vdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  w      ( wdims_s%i_start : wdims_s%i_end, &
           wdims_s%j_start : wdims_s%j_end, &
           wdims_s%k_start : wdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  u_adv  ( udims_l%i_start : udims_l%i_end, &
           udims_l%j_start : udims_l%j_end, &
           udims_l%k_start : udims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  v_adv  ( vdims_l%i_start : vdims_l%i_end, &
           vdims_l%j_start : vdims_l%j_end, &
           vdims_l%k_start : vdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  w_adv  ( wdims_l%i_start : wdims_l%i_end, &
           wdims_l%j_start : wdims_l%j_end, &
           wdims_l%k_start : wdims_l%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  theta  ( tdims_s%i_start : tdims_s%i_end, &
           tdims_s%j_start : tdims_s%j_end, &
           tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  rho    ( pdims_s%i_start : pdims_s%i_end, &
           pdims_s%j_start : pdims_s%j_end, &
           pdims_s%k_start : pdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  murk   ( tdims_s%i_start : tdims_s%i_end, &
           tdims_s%j_start : tdims_s%j_end, &
           tdims_s%k_start : tdims_s%k_end )

! Note that level-zero humidity fields are not currently updated for ENDGame.
! We may add support for this later.
REAL,    INTENT(INOUT) ::                   &
  q      ( qdims_l%i_start : qdims_l%i_end, &
           qdims_l%j_start : qdims_l%j_end, &
           qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qCL    ( qdims_l%i_start : qdims_l%i_end, &
           qdims_l%j_start : qdims_l%j_end, &
           qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qCF    ( qdims_l%i_start : qdims_l%i_end, &
           qdims_l%j_start : qdims_l%j_end, &
           qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) :: TStar      ( theta_field_size )
REAL,    INTENT(INOUT) :: tstar_anom ( theta_field_size )

REAL,    INTENT(INOUT) :: TStar_tile ( land_field, ntiles )

REAL,    INTENT(INOUT) :: TSoil      ( land_field, st_levels ), &
                          smcl       ( land_field, sm_levels )

REAL,    INTENT(INOUT) :: p_star     ( pdims_s%i_start : pdims_s%i_end,  &
                                       pdims_s%j_start : pdims_s%j_end )

REAL,    INTENT(INOUT) ::                                  &
  p                     ( pdims_s%i_start : pdims_s%i_end, &
                          pdims_s%j_start : pdims_s%j_end, &
                          pdims_s%k_start : pdims_s%k_end + 1 )

REAL,    INTENT(INOUT) ::                                  &
  p_theta_levels        ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                  &
  exner                 ( pdims_s%i_start : pdims_s%i_end, &
                          pdims_s%j_start : pdims_s%j_end, &
                          pdims_s%k_start : pdims_s%k_end + 1 )

REAL,    INTENT(INOUT) ::                                  &
  exner_theta_levels    ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) :: snow_depth ( theta_field_size )

REAL,    INTENT(INOUT) :: area_cloud_fraction          &
                        ( qdims%i_start : qdims%i_end, &
                          qdims%j_start : qdims%j_end, &
                                      1 : qdims%k_end )

REAL,    INTENT(INOUT) ::                                  &
  bulk_cloud_fraction   ( qdims_l%i_start : qdims_l%i_end, &
                          qdims_l%j_start : qdims_l%j_end, &
                          qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  cloud_fraction_liquid ( qdims_l%i_start : qdims_l%i_end, &
                          qdims_l%j_start : qdims_l%j_end, &
                          qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  cloud_fraction_frozen ( qdims_l%i_start : qdims_l%i_end, &
                          qdims_l%j_start : qdims_l%j_end, &
                          qdims_l%k_start : qdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  dust_div1             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div2             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div3             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div4             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div5             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div6             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                  &
  ozone_tracer          ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'IAU'

REAL, PARAMETER :: CC_tol = 1.0E-12 ! qCL/qCF tolerance for cloud clearing
REAL, PARAMETER :: oz_min = 1.0E-8  ! Minimum value allowed for ozone
                                    ! after addition of ozone increments

REAL, PARAMETER :: p_zero_recip = 1.0/p_zero
REAL, PARAMETER :: exner_exp    = 1.0/kappa - 1.0

! PC2 options.
! Seek advice from the PC2 team before altering these parameters from .TRUE.;
! you may need to have put in place large amounts of extra code first.
LOGICAL, PARAMETER :: L_pc2_cond = .TRUE. ! Do condensation and liquid cloud
                                          ! fraction changes?
LOGICAL, PARAMETER :: L_pc2_cfl  = .TRUE. ! Do liquid cloud fraction changes?
                                          ! This requires that the
                                          ! condensation as a result of
                                          ! assimilation is calculated either
                                          ! directly or via the estimate
                                          ! selected by setting l_pc2_cond to
                                          ! .TRUE..
                                          ! (Note that one must not run with
                                          ! condensation changes on but the
                                          ! liquid cloud fraction changes off.)
LOGICAL, PARAMETER :: L_pc2_cff  = .TRUE. ! Do ice cloud fraction changes?
                                          ! This requires that the ice
                                          ! increment from assimilation is
                                          ! calculated directly.

! Local variables:

INTEGER :: i, j, k
INTEGER :: TS_len_secs
INTEGER :: Sec
INTEGER :: IncNum
INTEGER :: FieldNum
INTEGER :: WeightNum
INTEGER :: tile_num
INTEGER :: field_size
INTEGER :: flds_StartAddr
INTEGER :: s_addr
INTEGER :: e_addr
INTEGER :: addr
INTEGER :: Code
INTEGER :: LocFldLen
INTEGER :: buflen      ! global field size read in
INTEGER :: ICode
INTEGER :: i_field     ! Counter for swap_bounds

REAL :: pBound
REAL :: maxInc
REAL :: delta_r
REAL :: thetaV_in_int
REAL :: thetaV_new_int
REAL :: Weight
REAL :: Weight_1
REAL :: Weight_2
REAL :: Term_1
REAL :: Term_2

REAL :: t1_in     ( tdims%i_start : tdims%i_end, &
                    tdims%j_start : tdims%j_end )
REAL :: t1_inc    ( theta_field_size )   ! work array for theta(1) perts
REAL :: p_tmp     ( pdims%i_start : pdims%i_end, &
                    pdims%j_start : pdims%j_end )
REAL :: q_sat_wat ( qdims%i_start : qdims%i_end, &
                    qdims%j_start : qdims%j_end )

! Array for holding a single IAU increment field:
REAL, ALLOCATABLE :: IAUIncFld ( : )

! Work arrays required for diagnosing humidity/cloud incs from qT incs:
REAL, ALLOCATABLE :: p_theta_levels_old        (:,:,:)
REAL, ALLOCATABLE :: qT                        (:,:,:)
REAL, ALLOCATABLE :: qT_plus                   (:,:,:)
REAL, ALLOCATABLE :: qSatW                     (:)
REAL, ALLOCATABLE :: VarTemp                   (:)
REAL, ALLOCATABLE :: AreaDivBulk               (:,:,:)
REAL, ALLOCATABLE :: cloud_fraction_liquid_0   (:)
REAL, ALLOCATABLE :: cloud_fraction_liquid_plus(:)
REAL, ALLOCATABLE :: qCL_0                     (:)
REAL, ALLOCATABLE :: qCL_plus                  (:)
REAL, ALLOCATABLE :: Var_qT                    (:)
REAL, ALLOCATABLE :: initialScaling            (:,:,:)
REAL, ALLOCATABLE :: scalingCoeff              (:,:,:)
REAL, ALLOCATABLE :: qcl_Inc                   (:,:,:)
REAL, ALLOCATABLE :: Cl_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Cf_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Cb_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Var_BGT                   (:)
REAL, ALLOCATABLE :: cloud_fraction_frozen_0   (:)
REAL, ALLOCATABLE :: cloud_fraction_frozen_plus(:)
REAL, ALLOCATABLE :: qCF_0                     (:)
REAL, ALLOCATABLE :: qCF_max                   (:)
REAL, ALLOCATABLE :: qCF_plus                  (:)

! Work arrays required for PC2 calculations:
REAL, ALLOCATABLE :: q_work   (:,:,:)
REAL, ALLOCATABLE :: qcl_work (:,:,:)
REAL, ALLOCATABLE :: qcf_in   (:,:,:)
REAL, ALLOCATABLE :: t_work   (:,:,:)
REAL, ALLOCATABLE :: p_work   (:,:,:)
REAL, ALLOCATABLE :: delta_q  (:,:,:)
REAL, ALLOCATABLE :: delta_qcl(:,:,:)
REAL, ALLOCATABLE :: delta_qcf(:,:,:)
REAL, ALLOCATABLE :: delta_t  (:,:,:)
REAL, ALLOCATABLE :: delta_p  (:,:,:)

! Work arrays required for call to QLimits:
REAL, ALLOCATABLE :: q_in       (:,:,:)
REAL, ALLOCATABLE :: pv_at_theta(:,:,:)

! Work arrays for sfctemp and soilm perturbations
REAL,ALLOCATABLE :: t2_inc    (:)
REAL,ALLOCATABLE :: s1_inc    (:)

! Other large work arrays to be allocated on demand:
REAL, ALLOCATABLE :: exner_in  (:,:,:)
REAL, ALLOCATABLE :: thetaV_in (:,:,:)
REAL, ALLOCATABLE :: thetaV_new(:,:,:)

TYPE(swapable_field_pointer_type) :: fields_to_swap(4)

LOGICAL :: IncsThisTS
LOGICAL :: qTIncsThisTS
LOGICAL :: CallQLimitsThisTS
LOGICAL :: FileOpen
LOGICAL :: ReadIAUField_FirstCallThisInc

CHARACTER(LEN=2)   :: IncNumStr
CHARACTER(LEN=9)   :: TimeStr
CHARACTER(LEN=336) :: CMessage

CHARACTER(MaxFileNameLen) :: FilePath

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook('IAU',zhook_in,zhook_handle)

! Array for holding a single IAU increment field:
IF (l_vatpoles) THEN
  ALLOCATE ( IAUIncFld ( v_field_size ) )
ELSE
  ALLOCATE ( IAUIncFld ( theta_field_size ) )
END IF ! vatpoles

! Timestep length in seconds:
TS_len_secs = SECS_PER_STEPIM(A_IM)

IF (PrintStatus >= PrStatus_Normal) THEN
  Sec = STEPim(a_im) * TS_len_secs
  WRITE (TimeStr,'(I3,":",I2.2,":",I2.2)')  &
    Sec/3600, MOD(Sec/60, 60), MOD(Sec, 60)
  WRITE (6,*) 'IAU: Entering IAU. Time since basis time (hr:mm:ss): ', TimeStr
END IF

!-------------------------------------------------------------------------------
! [1]: Return if there are no increments to add on this timestep.
!-------------------------------------------------------------------------------

IncsThisTS = .FALSE.
DO IncNum = 1, Num_IAU_incs
  IF (IAU_incs(IncNum) % FieldsToAdd) THEN
    IF (STEPim(a_im) >= IAU_incs(IncNum) % InsertionStartTS .AND. &
        STEPim(a_im) <= IAU_incs(IncNum) % InsertionEndTS) THEN
      WeightNum = STEPim(a_im) - IAU_incs(IncNum) % InsertionStartTS + 1
      IF (IAU_incs(IncNum) % Weights(WeightNum) /= 0.0) THEN
        IncsThisTS = .TRUE.
        EXIT
      END IF
    END IF
  END IF
END DO

IF (.NOT.IncsThisTS) THEN
  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE (6,*) 'IAU: No increments to add on this timestep'
    WRITE (6,*) 'IAU: Exiting IAU. Time since basis time (hr:mm:ss): ', TimeStr
  END IF
  RETURN
END IF

!-------------------------------------------------------------------------------
! [2]: If exner increments are to be calculated from p increments, make sure p
!      is consistent with exner.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcExnerIncs) THEN

  p(1:row_length,1:rows,:) = exner(1:row_length,1:rows,:)**recip_kappa &
                           * p_zero

END IF

!-------------------------------------------------------------------------------
! [3]: Determine whether there are qT increments to be added on this timestep.
!-------------------------------------------------------------------------------

qTIncsThisTS = .FALSE.
DO IncNum = 1, Num_IAU_incs
  IF (IAU_incs(IncNum) % Contains_qT) THEN
    IF (STEPim(a_im) >= IAU_incs(IncNum) % InsertionStartTS .AND. &
        STEPim(a_im) <= IAU_incs(IncNum) % InsertionEndTS) THEN
      qTIncsThisTS = .TRUE.
    END IF
  END IF
END DO

!-------------------------------------------------------------------------------
! [4]: Determine whether QLimits will be called on this timestep.
!-------------------------------------------------------------------------------

SELECT CASE (IAU_QLimitsCallFreq)
  CASE (CallFreq_Never)
    CallQLimitsThisTS = .FALSE.
  CASE (CallFreq_EveryCallIncsToAdd)
    CallQLimitsThisTS = .TRUE.
  CASE (CallFreq_EndInc1Window)
    CallQLimitsThisTS = ( IAU_incs(1) % FieldsToAdd .AND. &
                          STEPim(a_im) == IAU_incs(1) % InsertionEndTS )
  CASE (CallFreq_LastIAUCall)
    CallQLimitsThisTS = ( STEPim(a_im) == IAU_LastCallTS )
END SELECT

!-------------------------------------------------------------------------------
! [5]: Save fields required for calculating indirect field increments.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! [5.1]: Exner and thetaV.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

  ALLOCATE (exner_in (pdims%i_start:pdims%i_end, &
                      pdims%j_start:pdims%j_end, &
                      pdims%k_start:pdims%k_end+1))
  ALLOCATE (thetaV_in(tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                      tdims%k_start:tdims%k_end)  )

  exner_in (:,:,:)             = exner(pdims%i_start:pdims%i_end, &
                                       pdims%j_start:pdims%j_end,: )

  thetaV_in(:,:,qdims%k_start:qdims%k_end)  =                     &
                                 theta(tdims%i_start:tdims%i_end, &
                                       tdims%j_start:tdims%j_end, &
                                       qdims%k_start:qdims%k_end) &
                               * ( 1.0 + C_Virtual                &
                               *   q  (qdims%i_start:qdims%i_end, &
                                       qdims%j_start:qdims%j_end, &
                                       qdims%k_start:qdims%k_end) )

  IF (qdims%k_end < tdims%k_end) THEN
    thetaV_in(:,:,qdims%k_end+1: ) = theta(tdims%i_start:tdims%i_end, &
                                           tdims%j_start:tdims%j_end, &
                                           qdims%k_end+1: )
  END IF

END IF

!-------------------------------------------------------------------------------
! [5.2]: Level 1 temperature.
!-------------------------------------------------------------------------------

! Level 1 temperature:
IF (L_IAU_IncTStar) THEN

  t1_in(:,:) = exner_theta_levels(tdims%i_start:tdims%i_end,    &
                                  tdims%j_start:tdims%j_end, 1) &
             * theta             (tdims%i_start:tdims%i_end,    &
                                  tdims%j_start:tdims%j_end, 1)

END IF

!-------------------------------------------------------------------------------
! [5.3]: Fields required for processing qT increments.
!-------------------------------------------------------------------------------

IF (qTIncsThisTS) THEN

  ! NOTE : The zeroth theta level is being excluded in the initial
  !        implementation of ENDGAME IAU (UM 8.2) - to be reconsidered later.
  !        Moisture and temperature arrays therefore start from level 1.
  ALLOCATE (qT     (qdims%i_start:qdims%i_end, &
                    qdims%j_start:qdims%j_end, &
                                1:qdims%k_end))
  ALLOCATE (qT_plus(qdims%i_start:qdims%i_end, &
                    qdims%j_start:qdims%j_end, &
                                1:qdims%k_end))

  qT(:,:,:) = q  (qdims%i_start:qdims%i_end, &
                  qdims%j_start:qdims%j_end, &
                              1:qdims%k_end) &
            + qCL(qdims%i_start:qdims%i_end, &
                  qdims%j_start:qdims%j_end, &
                              1:qdims%k_end)

  IF (L_IAU_IncrementIce .AND. L_IAU_ApplyQTCorrections)   &
    qT(:,:,:) = qT(:,:,:) + qCF(qdims%i_start:qdims%i_end, &
                                qdims%j_start:qdims%j_end, &
                                            1:qdims%k_end)

  ! Initialise qT_plus:
  qT_plus(:,:,:) = qT(:,:,:)

  IF (.NOT.L_IAU_Add_qT_prime_to_q) THEN

    field_size = qdims%i_end * qdims%j_end * qdims%k_end

    ALLOCATE (p_theta_levels_old(qdims%i_start:qdims%i_end, &
                                 qdims%j_start:qdims%j_end, &
                                             1:qdims%k_end))

    ALLOCATE (qSatW  (field_size))
    ALLOCATE (VarTemp(field_size))

    p_theta_levels_old(:,:,:) = p_theta_levels(qdims%i_start:qdims%i_end, &
                                               qdims%j_start:qdims%j_end, &
                                                           1:qdims%k_end)

    VarTemp(:) = RESHAPE ( exner_theta_levels(qdims%i_start:qdims%i_end, &
                                              qdims%j_start:qdims%j_end, &
                                                          1:qdims%k_end) &
                         * theta             (qdims%i_start:qdims%i_end, &
                                              qdims%j_start:qdims%j_end, &
                                                          1:qdims%k_end) &
                         - qCL               (qdims%i_start:qdims%i_end, &
                                              qdims%j_start:qdims%j_end, &
                                                          1:qdims%k_end) &
                         * Lc/Cp                                         &
                         , (/field_size/) )

    CALL QSAT_WAT (qSatW,                              &
        VarTemp,                                       &
        RESHAPE(p_theta_levels                         &
                          (qdims%i_start:qdims%i_end,  &
                           qdims%j_start:qdims%j_end,  &
                                       1:qdims%k_end), &
                (/field_size/)), field_size)

    IF (L_IAU_IncrementIce) THEN
      ALLOCATE (Var_BGT(field_size))
      Var_BGT(:) = RESHAPE ( exner_theta_levels(qdims%i_start:qdims%i_end, &
                                                qdims%j_start:qdims%j_end, &
                                                            1:qdims%k_end) &
                           * theta             (qdims%i_start:qdims%i_end, &
                                                qdims%j_start:qdims%j_end, &
                                                            1:qdims%k_end) &
                           , (/field_size/) )
    END IF

  END IF ! (.NOT.L_IAU_Add_qT_prime_to_q)

END IF ! (qTIncsThisTS)

!-------------------------------------------------------------------------------
! [5.4]: Fields required for PC2 calculations.
!-------------------------------------------------------------------------------

IF (.NOT.qTIncsThisTS .AND. &
    L_pc2 .AND. (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff)) THEN

  ALLOCATE (q_work  (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end))
  ALLOCATE (qcl_work(qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end))
  ALLOCATE (qcf_in  (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end))
  ALLOCATE (p_work  (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end))
  ALLOCATE (t_work  (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end))

  ! Calculate increments to vapour, cloud, temperature and pressure:
  DO k =             1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i= qdims%i_start, qdims%i_end
        q_work  (i,j,k) = q  (i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k)
        qcf_in  (i,j,k) = qcf(i,j,k)
        p_work  (i,j,k) = p_theta_levels    (i,j,k)
        t_work  (i,j,k) = exner_theta_levels(i,j,k) * theta(i,j,k)
      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [5.5]: Fields required for QLimits call.
!-------------------------------------------------------------------------------

IF (CallQLimitsThisTS) THEN

  ALLOCATE (q_in(qdims%i_start:qdims%i_end, &
                 qdims%j_start:qdims%j_end, &
                             1:qdims%k_end))

  q_in(:,:,:) = q(qdims%i_start:qdims%i_end, &
                  qdims%j_start:qdims%j_end, 1: )

END IF

!-------------------------------------------------------------------------------
! [6]: Loop over IAU increment files, and perform direct field updates.
!-------------------------------------------------------------------------------

DO IncNum = 1, Num_IAU_incs

  ! Get increment number as a string:
  WRITE (IncNumStr,'(I2.2)') IncNum

  ! Cycle if there are no increments to add on this timestep;
  IF (          .NOT.IAU_incs(IncNum) % FieldsToAdd      .OR. &
      STEPim(a_im) < IAU_incs(IncNum) % InsertionStartTS .OR. &
      STEPim(a_im) > IAU_incs(IncNum) % InsertionEndTS) CYCLE

  ! Get and report weight:
  WeightNum = STEPim(a_im) - IAU_incs(IncNum) % InsertionStartTS + 1
  Weight    = IAU_incs(IncNum) % Weights(WeightNum)
  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE (6,'(A,ES17.10,A)')                              &
      ' IAU: Increments to be added with weight ', Weight, &
      ' for inc no.'// IncNumStr
  END IF

  FileOpen = .FALSE.

  !-----------------------------------------------------------------------------
  ! [6.1]: If required, allocate storage space for increments.
  !-----------------------------------------------------------------------------

  IF (     IAU_incs(IncNum) % StoreFields .AND. &
      .NOT.IAU_incs(IncNum) % FieldsStored) THEN

    IF (IAU_incs(IncNum) % StorePacked) THEN
      ALLOCATE (IAU_incs(IncNum) % flds32bit(IAU_incs(IncNum) % flds_len))
    ELSE
      ALLOCATE (IAU_incs(IncNum) % flds     (IAU_incs(IncNum) % flds_len))
    END IF

  END IF

  !-----------------------------------------------------------------------------
  ! [6.2]: Loop over fields in the increment file.
  !-----------------------------------------------------------------------------

  flds_StartAddr = 1 ! Field start address in flds or flds32bit array

  ReadIAUField_FirstCallThisInc = .TRUE.

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup

    Code    = IAU_incs(IncNum) % Lookup(ITEM_CODE, FieldNum)
    k       = IAU_incs(IncNum) % Lookup(LBLEV,     FieldNum)
    buflen  = IAU_incs(IncNum) % Lookup(LBLREC,    FieldNum)

    ! Get local field length:
    LocFldLen = 0
    DO i = 1, IAU_NumFldCodes
      IF (IAU_FldCodes(i) == Code) LocFldLen = IAU_LocFldLens(i)
    END DO

    IF (LocFldLen == 0) CYCLE ! Unused field

    !---------------------------------------------------------------------------
    ! [6.2.1]: Obtain increment field, and save for next call if necessary.
    !---------------------------------------------------------------------------

    IF (IAU_incs(IncNum) % FieldsStored) THEN ! Field already in memory

      s_addr = flds_StartAddr
      e_addr = flds_StartAddr + LocFldLen - 1

      IF (IAU_incs(IncNum) % StorePacked) THEN
        ! Unpack field:
        CALL expand21 ( LocFldLen,                                   & ! in
                        IAU_incs(IncNum) % flds32bit(s_addr:e_addr), & ! in
                        IAUIncFld(1:LocFldLen) )                       ! out
      ELSE
        IAUIncFld(1:LocFldLen) = IAU_incs(IncNum) % flds(s_addr:e_addr)
      END IF

    ELSE ! Read field from file

      FilePath = IAU_incs(IncNum) % FilePath

      IF (.NOT.FileOpen) THEN
        CALL MODEL_FILE_OPEN (IAU_unit, FilePath, MaxFileNameLen, 0, 1, ICode)
        IF (ICode > 0) THEN
          CMessage(1:80) = 'Error opening IAU increment no.'// IncNumStr //':'
          CMessage(81:)  = TRIM(FilePath)
          CALL EReport (RoutineName, ICode, CMessage)
        END IF
        FileOpen = .TRUE.
      END IF

      CALL ReadIAUField (                                &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                          IAU_incs(IncNum) % FixHd,      & ! in
                          IAU_incs(IncNum) % Len1Lookup, & ! in
                          IAU_incs(IncNum) % Len2Lookup, & ! in
                          IAU_incs(IncNum) % Lookup,     & ! inout
                          IncNum,                        & ! in
                          FieldNum,                      & ! in
                          LocFldLen,                     & ! in
                          ReadIAUField_FirstCallThisInc, & ! in
                          IAUIncFld(1:LocFldLen) )         ! out

      ReadIAUField_FirstCallThisInc = .FALSE.

      IF (IAU_incs(IncNum) % StoreFields) THEN

        s_addr = flds_StartAddr
        e_addr = flds_StartAddr + LocFldLen - 1

        IF (IAU_incs(IncNum) % StorePacked) THEN
          CALL pack21 ( LocFldLen,              &                     ! in
                        IAUIncFld(1:LocFldLen), &                     ! in
                        IAU_incs(IncNum) % flds32bit(s_addr:e_addr) ) ! out
        ELSE
          IAU_incs(IncNum) % flds(s_addr:e_addr) = IAUIncFld(1:LocFldLen) 
        END IF

      END IF

    END IF

    !---------------------------------------------------------------------------
    ! [6.2.2]: Add the increment to the relevant model field, applying limits
    !          where necessary.
    !---------------------------------------------------------------------------

    addr = 1

    ! u and u_adv:
    IF (Code == 2 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          u    (i,j,k) = u    (i,j,k) + Weight * IAUIncFld(addr)
          u_adv(i,j,k) = u_adv(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! v and v_adv:
    IF (Code == 3 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          v    (i,j,k) = v    (i,j,k) + Weight * IAUIncFld(addr)
          v_adv(i,j,k) = v_adv(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF


    ! w and w_adv:
    IF ( Code == 150 .AND. &
         (k == 9999 .OR. (k >= 1 .AND. k <= model_levels)) ) THEN
      IF (k == 9999) k = 0
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          w    (i,j,k) = w    (i,j,k) + Weight * IAUIncFld(addr)
          w_adv(i,j,k) = w_adv(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! theta:
    IF (Code == 4 .AND. k >= 1 .AND. k <= model_levels) THEN
      ! ** Want to remove the following extra test at some point:
      IF (.NOT.(L_IAU_IgnoreTopLevThetaIncs .AND. k == model_levels)) THEN

        ! Apply limit to size of upper-level theta increments?
        IF (L_IAU_LimitUpperThetaIncs) THEN

          ! Shorten variable names:
          pBound = IAU_LimitUpperThetaIncs_pBound
          maxInc = IAU_LimitUpperThetaIncs_maxInc

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end

              IF (p(i,j,k) < pBound .AND. &
                  ABS(Weight * IAUIncFld(addr)) > maxInc) THEN
                WRITE(6,*) 'IAU: Theta inc (', Weight * IAUIncFld(addr), &
                           ',) restricted at level ', k, ' pressure ',   &
                           p(i,j,k), ' for inc no.'// IncNumStr

                IF (Weight * IAUIncFld(addr) < 0) THEN
                   theta(i,j,k) = theta(i,j,k) - maxInc
                ELSE
                   theta(i,j,k) = theta(i,j,k) + maxInc
                END IF
              ELSE
                theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(addr)
              END IF
              addr = addr + 1
            END DO
          END DO

        ELSE

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(addr)
              addr = addr + 1
            END DO
          END DO

        END IF ! (L_IAU_LimitUpperThetaIncs)

      END IF
    END IF

    ! rho:
    IF (Code == 253 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho(i,j,k) = rho(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! exner:
    IF (Code == 255 .AND. k >= 1 .AND. k <= model_levels+1) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          exner(i,j,k) = exner(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! p:
    IF (Code == 407 .AND. k >= 1 .AND. k <= model_levels+1) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          p(i,j,k) = p(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! q:
    IF (Code == 10 .AND. k >= 1 .AND. k <= wet_levels &
                                .AND. k <  model_levels) THEN
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q(i,j,k) = q(i,j,k) + Weight * IAUIncFld(addr)
          ! Apply lower limit:
          IF (q(i,j,k) < q_min) q(i,j,k) = q_min
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qCL:
    IF (Code == 254 .AND. k >= 1 .AND. k <= wet_levels) THEN
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qCL(i,j,k) = qCL(i,j,k) + Weight * IAUIncFld(addr)
          ! Check for non-positive qCLs:
          IF (qCL(i,j,k) <= 0.0) THEN
            qCL(i,j,k) = 0.0
            cloud_fraction_liquid(i,j,k) = 0.0
            IF (qCF(i,j,k) <= CC_tol) THEN
              area_cloud_fraction(i,j,k) = 0.0
              bulk_cloud_fraction(i,j,k) = 0.0
            END IF
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qCF:
    IF (Code == 12 .AND. k >= 1 .AND. k <= wet_levels) THEN
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qCF(i,j,k) = qCF(i,j,k) + Weight * IAUIncFld(addr)
          ! Check for non-positive qCFs:
          IF (qCF(i,j,k) <= 0.0) THEN
            qCF(i,j,k) = 0.0
            cloud_fraction_frozen(i,j,k) = 0.0
            IF (qCL(i,j,k) <= CC_tol) THEN
              area_cloud_fraction(i,j,k) = 0.0
              bulk_cloud_fraction(i,j,k) = 0.0
            END IF
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qT:
    IF ((Code == 16207 .OR. Code == 18001) .AND. &
        k >= 1 .AND. k <= wet_levels) THEN
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          qT_plus(i,j,k) = qT_plus(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
          IF (L_IAU_ApplyQTCorrections) THEN
            ! Apply lower limit:
            IF (qT_plus(i,j,k) < q_min) qT_plus(i,j,k) = q_min
          END IF
        END DO
      END DO
    END IF

    ! Aerosol:
    IF (Code == 90 .AND. k >= 1 .AND. k <= bl_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          murk(i,j,k) = murk(i,j,k) + Weight * IAUIncFld(addr)
          ! Apply upper (AEROMAX) and lower (AERO0) limits:
          murk(i,j,k) = MIN(MAX(AERO0,murk(i,j,k)),AEROMAX)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! dust_div1:
    IF (Code == 431 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div1(i,j,k) = dust_div1(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div1(i,j,k) <= 0.0) THEN
            dust_div1(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div2:
    IF (Code == 432 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div2(i,j,k) = dust_div2(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div2(i,j,k) <= 0.0) THEN
            dust_div2(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div3:
    IF (Code == 433 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div3(i,j,k) = dust_div3(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div3(i,j,k) <= 0.0) THEN
            dust_div3(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div4:
    IF (Code == 434 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div4(i,j,k) = dust_div4(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div4(i,j,k) <= 0.0) THEN
            dust_div4(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div5:
    IF (Code == 435 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div5(i,j,k) = dust_div5(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div5(i,j,k) <= 0.0) THEN
            dust_div5(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div6:
    IF (Code == 436 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div6(i,j,k) = dust_div6(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div6(i,j,k) <= 0.0) THEN
            dust_div6(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! Ozone tracer:                              
    IF (Code == 480 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ozone_tracer(i,j,k) = ozone_tracer(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive ozone tracer:
          IF (ozone_tracer(i,j,k) <= 0.0) THEN
            ozone_tracer(i,j,k) = 0.0
            IF(L_IAU_SetOzoneMin) ozone_tracer(i,j,k) = oz_min
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! SfcTemp: SST (sea) & TStar/TSoil (land)
    IF (Code == 24 .AND. k == 9999) THEN

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (6,*) 'IAU: Adding surface temp perts to TStar + TSoil(1)'
      END IF
      ALLOCATE (t2_inc(theta_field_size))

      DO i = 1, theta_field_size
        t2_inc(i) = Weight * IAUIncFld(addr)
        addr = addr + 1
      END DO

      ! Add perts to TStar (Land and Sea)
      DO i = 1, theta_field_size
        TStar(i) = TStar(i) + t2_inc(i)
      END DO
      ! Add perts to tstar_anom if SST anomaly field is being used,
      ! otherwise anomaly update in REPLANCA will overwrite pert values
      IF (l_sstanom) THEN
        DO i = 1, theta_field_size
          tstar_anom(i) = tstar_anom(i) + t2_inc(i)
        END DO
      END IF

      ! Assuming an uncertainty in the snow depth, we perturb TSoil and
      ! TStar_tiles up to double the depth of the threshold used for
      ! analysis incr: i.e. 0.1 kg/m**2 instead of 0.05 kg/m**2

      ! Add sfctemp perts to TSoil (land)
      DO i = 1, land_field

        IF ( snow_depth(land_index(i)) <= 0.1 ) THEN
          TSoil(i, 1) = TSoil(i, 1) + t2_inc(land_index(i))
        END IF

      END DO

      ! Also add increments to TStar_tile:
      DO tile_num = 1, ntiles
        DO i = 1, land_field

          IF ( snow_depth(land_index(i)) <= 0.1 ) THEN
            TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &
                                   + t2_inc(land_index(i))
          ENDIF

        END DO
      END DO

      DEALLOCATE (t2_inc)

    END IF

    ! Soil Moisture (sm_levels)
    IF (Code == 9 .AND. k >= 1 .AND. k <= sm_levels) THEN

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (6,*) 'IAU: Adding soil moisture perts to smcl'
      END IF

      ! Check whether data on land-points or full 2-D grid
      field_size = global_row_length * global_rows
      IF (buflen == field_size) THEN
        ALLOCATE (s1_inc(theta_field_size))

        DO i = 1, theta_field_size
          s1_inc(i) = Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO

        DO i = 1, land_field
          smcl(i,k) = smcl(i,k) + s1_inc(land_index(i))
        END DO

        IF (ALLOCATED(s1_inc)) DEALLOCATE (s1_inc)

      ELSE
        ! Land-points only
        DO i = 1, land_field
          smcl(i,k) = smcl(i,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO

      END IF


    END IF

    !---------------------------------------------------------------------------
    ! [6.2.3]: Update start address in flds or flds32bit array.
    !---------------------------------------------------------------------------
                                                                        
    flds_StartAddr = flds_StartAddr + LocFldLen

  END DO ! FieldNum

  !-----------------------------------------------------------------------------
  ! [6.3]: Set FieldsStored flag.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % StoreFields) IAU_incs(IncNum) % FieldsStored = .TRUE.

  !-----------------------------------------------------------------------------
  ! [6.4]: Close increment file.
  !-----------------------------------------------------------------------------

  IF (FileOpen) THEN
    CALL MODEL_FILE_CLOSE (IAU_unit, FilePath, MaxFileNameLen, 1, 0, ICode)
    IF (ICode > 0) THEN
      CMessage(1:80) = 'Error closing IAU increment no.'// IncNumStr //':'
      CMessage(81:)  = TRIM(FilePath)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
    FileOpen = .FALSE.
  END IF

END DO ! IncNum

!-------------------------------------------------------------------------------
! [7]: Perform indirect field updates.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! [7.1]: exner/p.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcExnerIncs) THEN
  ! Calculate exner from p:
  exner(1:row_length,1:rows,:) = ( p(1:row_length,1:rows,:) &
                                 * p_zero_recip             &
                                 )**kappa
ELSE                          ! Calculate p from exner
  ! Calculate p from exner:
  p(1:row_length,1:rows,:) =                                &
            p_zero * exner(1:row_length,1:rows,:)**recip_kappa
END IF

! Update pstar:
CALL Calc_P_star ( r_theta_levels, r_rho_levels,   & ! in
                   p, rho,                         & ! in
                   row_length, rows, model_levels, & ! in
                   offx, offy, halo_i, halo_j,     & ! in
                   p_star )                          ! out

! Update exner on theta levels:
CALL Calc_Exner_at_theta ( r_theta_levels, r_rho_levels,   & ! in
                           exner,                          & ! in
                           row_length, rows, model_levels, & ! in
                           offx, offy, halo_i, halo_j,     & ! in
                           exner_theta_levels,             & ! out
                           .FALSE. )                         ! out

! Get p on theta levels from exner on theta levels:
CALL Calc_P_from_Exner ( p_theta_levels,                 & ! out
                         row_length, rows,               & ! in
                         tdims_s%k_end-tdims_s%k_start+1,& ! in
                         offx, offy,                     & ! in
                         exner_theta_levels,.FALSE. )      ! in

!-------------------------------------------------------------------------------
! [7.2]: thetaV.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

  ALLOCATE (thetaV_new(tdims%i_start:tdims%i_end, &
                       tdims%j_start:tdims%j_end, &
                                   1:tdims%k_end))

  ! Add hydrostatic thetaV increments:
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        IF (k /= model_levels) THEN
          delta_r = r_rho_levels(i,j,k+1) &
                  - r_rho_levels(i,j,k  )
        ELSE
          delta_r = 2.0                     &
                  * ( r_theta_levels(i,j,k) &
                    - r_rho_levels  (i,j,k) &
                    )
        END IF

        thetaV_new(i,j,k) = thetaV_in(i,j,k)       &
                          - g/cp * delta_r         &
                          * ( 1.0                  &
                            / ( exner    (i,j,k+1) &
                              - exner    (i,j,k  ) &
                            )                      &
                            - 1.0                  &
                            / ( exner_in(i,j,k+1)  &
                              - exner_in(i,j,k  )  &
                            ) )

      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [7.3]: rho.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcRhoIncs) THEN

  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        ! Interpolate old and new thetaVs to rho points:
        IF (k /= 1) THEN

          Weight_1 = r_rho_levels  (i,j,k  ) &
                   - r_theta_levels(i,j,k-1)
          Weight_2 = r_theta_levels(i,j,k  ) &
                   - r_rho_levels  (i,j,k  )

          thetaV_in_int  = ( Weight_1 * thetaV_in (i,j,k  ) &
                           + Weight_2 * thetaV_in (i,j,k-1) &
                           )                                &
                         / ( Weight_1 + Weight_2 )

          thetaV_new_int = ( Weight_1 * thetaV_new(i,j,k  ) &
                           + Weight_2 * thetaV_new(i,j,k-1) &
                           )                                &
                         / ( Weight_1 + Weight_2 )

        ELSE

          thetaV_in_int  = thetaV_in (i,j,k)
          thetaV_new_int = thetaV_new(i,j,k)

        END IF

        Term_1 = exner   (i,j,k)**exner_exp / thetaV_new_int
        Term_2 = exner_in(i,j,k)**exner_exp / thetaV_in_int

        rho(i,j,k) = rho(i,j,k)             &
                   + ( Term_1 - Term_2 )    &
                   * r_rho_levels(i,j,k)**2 &
                   * p_zero / (kappa * cp)

      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [7.4]: theta.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs) THEN

  theta(qdims%i_start:qdims%i_end,          &
        qdims%j_start:qdims%j_end,          &
                    1:qdims%k_end) =        &
             thetaV_new(:,:,1:qdims%k_end)  &
          / ( 1.0 + C_Virtual               &
             * q(qdims%i_start:qdims%i_end, &
                 qdims%j_start:qdims%j_end, &
                             1:qdims%k_end) &
             )

  theta(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
                      qdims%k_end+1:) = thetaV_new(:,:,qdims%k_end+1:)

END IF

!-------------------------------------------------------------------------------
! [7.5]: Deallocate work arrays.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN
  DEALLOCATE (exner_in )
  DEALLOCATE (thetaV_in)
  DEALLOCATE (thetaV_new)
END IF

!-------------------------------------------------------------------------------
! [7.6]: Humidity and cloud updates via qT increments.
!-------------------------------------------------------------------------------

IF (qTIncsThisTS .AND. .NOT.L_IAU_Add_qT_prime_to_q) THEN

  ALLOCATE (AreaDivBulk(qdims%i_start:qdims%i_end, &
                        qdims%j_start:qdims%j_end, &
                                    1:qdims%k_end))

  ALLOCATE (cloud_fraction_liquid_0   (field_size))
  ALLOCATE (cloud_fraction_liquid_plus(field_size))
  ALLOCATE (qCL_0                     (field_size))
  ALLOCATE (qCL_plus                  (field_size))
  ALLOCATE (Var_qT                    (field_size))

  !----------------------------------------------------------------------------
  ! [7.6.1]: obtain initial and incremented values from Var scheme
  !-----------------------------------------------------------------------------

  ! Use Var cloud diagnostic to obtain  initial T, from qT & Tl
  Var_qT  = RESHAPE(qT(qdims%i_start:qdims%i_end, &
                       qdims%j_start:qdims%j_end, &
                                   1:qdims%k_end), (/field_size/))

  CALL Var_DiagCloud(                                     &
      field_size,                                         &
      RESHAPE(p_theta_levels_old, (/field_size/)),        &
      RESHAPE(SPREAD(SPREAD(                              &
              rhcrit(1:wet_levels),1,rows),1,row_length), &
              (/field_size/)),                            &
      L_IAU_IncrementIce,                                 &
      Var_qT,                                             &
      CMessage,                                           &
      ICode,                                              &
      TL=VarTemp)

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  ! (VarTemp now holds initial T values)

  IF (L_IAU_IncrementIce) THEN

    ALLOCATE (cloud_fraction_frozen_0   (field_size))
    ALLOCATE (cloud_fraction_frozen_plus(field_size))
    ALLOCATE (qCF_0                     (field_size))
    ALLOCATE (qCF_max                   (field_size))
    ALLOCATE (qCF_plus                  (field_size))

    ! Use Var cloud diagnostic to obtain  initial qCL,qCF,Cl & Cf
    qcl_0 = RESHAPE(qCL (qdims%i_start:qdims%i_end, &
                         qdims%j_start:qdims%j_end, &
                                     1:qdims%k_end), (/field_size/))
    qcf_0 = RESHAPE(qCF (qdims%i_start:qdims%i_end, &
                         qdims%j_start:qdims%j_end, &
                                     1:qdims%k_end), (/field_size/))

    CALL Var_DiagCloud(                                     &
        field_size,                                         &
        RESHAPE(p_theta_levels_old, (/field_size/)),        &
        RESHAPE(SPREAD(SPREAD(                              &
                rhcrit(1:wet_levels),1,rows),1,row_length), &
                (/field_size/)),                            &
        L_IAU_IncrementIce,                                 &
        Var_qT,                                             &
        CMessage,                                           &
        ICode,                                              &
        CL    = cloud_fraction_liquid_0,                    &
        qCL   = qCL_0,                                      &
        BGqCL = RESHAPE(qCL (qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                        (/field_size/)),                    &
        qCF   = qCF_0,                                      &
        CF    = cloud_fraction_frozen_0,                    &
        BGqCF = RESHAPE(qCF (qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                        (/field_size/)),                    &
        BGT   = Var_BGT,                                    &
        T     = VarTemp)

  ELSE

    ! Use Var cloud diagnostic to obtain  initial qCL & Cl

    CALL Var_DiagCloud(                                     &
        field_size,                                         &
        RESHAPE(p_theta_levels_old, (/field_size/)),        &
        RESHAPE(SPREAD(SPREAD(                              &
                rhcrit(1:wet_levels),1,rows),1,row_length), &
                (/field_size/)),                            &
        L_IAU_IncrementIce,                                 &
        Var_qT,                                             &
        CMessage,                                           &
        ICode,                                              &
        CL    = cloud_fraction_liquid_0,                    &
        qCL   = qCL_0,                                      &
        T     = VarTemp)

  END IF

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  ! Assign VarTemp to hold incremented T values:
  IF (L_IAU_IncrementIce) THEN

    ! Assign Var_BGT to store the T associated with qcl_0
    ! Use qCL_Plus as a temporary variable
    ! This should be done in the qcl / Cl version too but
    !   will go in later as a bug-fix so that bit-comparison
    !   can be maintained when the code for qcf incrementing is
    !   implemented but switched off
    qCL_Plus = VarTemp
    VarTemp  = VarTemp + (                              &
         RESHAPE(                                       &
         exner_theta_levels(qdims%i_start:qdims%i_end,  &
                            qdims%j_start:qdims%j_end,  &
                                        1:qdims%k_end)  &
         * theta(qdims%i_start:qdims%i_end,             &
                            qdims%j_start:qdims%j_end,  &
                                        1:qdims%k_end), &
        (/field_size/)) - Var_BGT)
    Var_BGT  = qCL_Plus

  ELSE

    VarTemp = RESHAPE(                                 &
         exner_theta_levels(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
         * theta(qdims%i_start:qdims%i_end,            &
                 qdims%j_start:qdims%j_end,            &
                             1:qdims%k_end),           &
        (/field_size/))

  END IF

  IF (L_IAU_IncrementIce) THEN

    ! Calculate maximum qcf values (for use later)
    qCL_Plus = 0.0
    qCF_Plus = 0.0
    Var_qT = RESHAPE(qT_plus(qdims%i_start:qdims%i_end,  &
                             qdims%j_start:qdims%j_end,  &
                                         1:qdims%k_end)  &
                     +  qCF (qdims%i_start:qdims%i_end,  &
                             qdims%j_start:qdims%j_end,  &
                                         1:qdims%k_end), &
                     (/field_size/))
    qCF_max = Var_qT

    CALL Var_DiagCloud(                                     &
        field_size,                                         &
        RESHAPE(p_theta_levels                              &
                (qdims%i_start:qdims%i_end,                 &
                 qdims%j_start:qdims%j_end,                 &
                             1:qdims%k_end),                &
                (/field_size/)),                            &
        RESHAPE(SPREAD(SPREAD(                              &
                rhcrit(1:wet_levels),1,rows),1,row_length), &
                (/field_size/)),                            &
        L_IAU_IncrementIce,                                 &
        Var_qT,                                             &
        CMessage,                                           &
        ICode,                                              &
        qCL    = qCL_Plus,                                  &
        BGqCL  = RESHAPE(qCL(qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                         (/field_size/)),                   &
        qCF    = qCF_Plus,                                  &
        qCFmax = qCF_max,                                   &
        BGqCF  = RESHAPE(qCF(qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                         (/field_size/)),                   &
        BGT    = Var_BGT,                                   &
        T      = VarTemp)

    IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

    qCL_Plus  = qCL_0
    qCF_Plus  = qCF_0
    Var_qT = RESHAPE(qT     (qdims%i_start:qdims%i_end,  &
                             qdims%j_start:qdims%j_end,  &
                                         1:qdims%k_end)  &
                     +  qCF (qdims%i_start:qdims%i_end,  &
                             qdims%j_start:qdims%j_end,  &
                                         1:qdims%k_end), &
                     (/field_size/))

    WHERE (qCF_max <= 0.0)

      qcf_max = Var_qT

    ENDWHERE

  END IF

  Var_qT = RESHAPE(qT_plus(qdims%i_start:qdims%i_end, &
                        qdims%j_start:qdims%j_end,    &
                                    1:qdims%k_end),   &
                   (/field_size/))

  IF (L_IAU_IncrementIce) THEN

    ! Use Var cloud diagnostic to obtain  final qCL,qCF,Cl & Cf
    CALL Var_DiagCloud(                                     &
        field_size,                                         &
        RESHAPE(p_theta_levels                              &
                (qdims%i_start:qdims%i_end,                 &
                 qdims%j_start:qdims%j_end,                 &
                             1:qdims%k_end),                &
                (/field_size/)),                            &
        RESHAPE(SPREAD(SPREAD(                              &
                rhcrit(1:wet_levels),1,rows),1,row_length), &
                (/field_size/)),                            &
        L_IAU_IncrementIce,                                 &
        Var_qT,                                             &
        CMessage,                                           &
        ICode,                                              &
        qCL    = qCL_Plus,                                  &
        BGqCL  = RESHAPE(qCL(qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                         (/field_size/)),                   &
        qCF    = qCF_Plus,                                  &
        BGqCF  = RESHAPE(qCF(qdims%i_start:qdims%i_end,     &
                             qdims%j_start:qdims%j_end,     &
                                         1:qdims%k_end),    &
                         (/field_size/)),                   &
        qCFmax = qCF_max,                                   &
        CL     = cloud_fraction_liquid_plus,                &
        CF     = cloud_fraction_frozen_plus,                &
        BGT    = Var_BGT,                                   &
        T      = VarTemp)

  ELSE

    ! Use Var cloud diagnostic to obtain  final qCL & Cl
    CALL Var_DiagCloud(                                     &
        field_size,                                         &
        RESHAPE(p_theta_levels                              &
                (qdims%i_start:qdims%i_end,                 &
                 qdims%j_start:qdims%j_end,                 &
                             1:qdims%k_end),                &
                (/field_size/)),                            &
        RESHAPE(SPREAD(SPREAD(                              &
                rhcrit(1:wet_levels),1,rows),1,row_length), &
                (/field_size/)),                            &
        L_IAU_IncrementIce,                                 &
        Var_qT,                                             &
        CMessage,                                           &
        ICode,                                              &
        CL     = cloud_fraction_liquid_plus,                &
        qCL    = qCL_Plus,                                  &
        T      = VarTemp)

  END IF

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  !-----------------------------------------------------------------------------
  ! [7.6.2]: apply increments
  !-----------------------------------------------------------------------------

  ! qcl

  ALLOCATE (qcl_Inc(qdims%i_start:qdims%i_end, &
                    qdims%j_start:qdims%j_end, &
                                1:qdims%k_end))

  IF (L_IAU_ScaleCloud) THEN

    ! Scale qcl increments to be physical values
    ALLOCATE (initialScaling(qdims%i_start:qdims%i_end, &
                             qdims%j_start:qdims%j_end, &
                                         1:qdims%k_end))
    ALLOCATE (scalingCoeff  (qdims%i_start:qdims%i_end, &
                             qdims%j_start:qdims%j_end, &
                                         1:qdims%k_end))
    WHERE(RESHAPE(qCL_0,(/row_length,rows,wet_levels/)) >      &
          qCL(qdims%i_start:qdims%i_end,                       &
              qdims%j_start:qdims%j_end,                       &
                          1:qdims%k_end) .AND.                 &
      RESHAPE(qCL_plus - qCL_0,(/row_length,rows,wet_levels/)) &
      < 0.0)
      qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                        (/row_length, rows, wet_levels/))
      WHERE (qCL(qdims%i_start:qdims%i_end, &
                 qdims%j_start:qdims%j_end, &
                             1:qdims%k_end) > 0.0)
        initialScaling = -qCL(qdims%i_start:qdims%i_end,  &
                              qdims%j_start:qdims%j_end,  &
                                          1:qdims%k_end)* &
            (1.0 - exp(qcl_Inc /                          &
                       qCL(qdims%i_start:qdims%i_end,     &
                           qdims%j_start:qdims%j_end,     &
                                       1:qdims%k_end)))
        scalingCoeff = -qCL(qdims%i_start:qdims%i_end,           &
                            qdims%j_start:qdims%j_end,           &
                                        1:qdims%k_end)*          &
        (1.0-exp(-RESHAPE(qCL_0,(/row_length,rows,wet_levels/))/ &
        qCL(qdims%i_start:qdims%i_end,                           &
            qdims%j_start:qdims%j_end,                           &
                        1:qdims%k_end)))
        ! The following calculation should remain >0 so no error
        ! trap needed
        scalingCoeff =                                      &
            (scalingCoeff +                                 &
             qCL(qdims%i_start:qdims%i_end,                 &
                 qdims%j_start:qdims%j_end,                 &
                             1:qdims%k_end)) /              &
            (scalingCoeff +                                 &
             RESHAPE(qCL_0,(/row_length,rows,wet_levels/)))
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      qcl_Inc = initialScaling - &
                (initialScaling - qcl_Inc) * scalingCoeff
    ELSEWHERE
      qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                        (/row_length, rows, wet_levels/))
    END WHERE

  ELSE

    qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                      (/row_length, rows, wet_levels/))

  END IF

  ! limit qcl increment:
  IF (L_IAU_LimitIncOp) THEN
    WHERE (ABS(qcl_Inc(qdims%i_start:qdims%i_end,   &
                       qdims%j_start:qdims%j_end,   &
                       1:qdims%k_end)) > &
      IAU_qclthreshscale * ABS(qT_plus(:,:,:) - qT(:,:,:)))
      qcl_Inc(qdims%i_start:qdims%i_end,   &
              qdims%j_start:qdims%j_end,   &
              1:qdims%k_end) =             &
      SIGN(IAU_qclthreshscale*ABS(qT_plus(:,:,:)-qT(:,:,:)),&
                                  qcl_inc(qdims%i_start:qdims%i_end,  &
                                           qdims%j_start:qdims%j_end, &
                                            1:qdims%k_end)            &
                                   ) 
    END WHERE
  ENDIF

  ! Update:
  qCL(qdims%i_start:qdims%i_end,   &
      qdims%j_start:qdims%j_end,   &
                  1:qdims%k_end) = &
  qCL(qdims%i_start:qdims%i_end,   &
      qdims%j_start:qdims%j_end,   &
                  1:qdims%k_end) + qcl_Inc(:,:,:)

  ! qcf
  IF (L_IAU_IncrementIce) THEN

    qCF(qdims%i_start:qdims%i_end,   &
        qdims%j_start:qdims%j_end,   &
                    1:qdims%k_end) = &
    qCF(qdims%i_start:qdims%i_end,   &
        qdims%j_start:qdims%j_end,   &
                    1:qdims%k_end) + &
        RESHAPE(qCF_plus - qCF_0,    &
                (/row_length, rows, wet_levels/))

  END IF
  
  ! Store input values of area / bulk cloud fractions
  IF (L_CLD_AREA) THEN
    
    WHERE (bulk_cloud_fraction(qdims%i_start:qdims%i_end, &
                               qdims%j_start:qdims%j_end, &
                                           1:qdims%k_end) &
           <= 0.0)

      AreaDivBulk = 1.0

    ELSEWHERE

      AreaDivBulk =                                      &
        area_cloud_fraction(qdims%i_start:qdims%i_end,   &
                            qdims%j_start:qdims%j_end,   &
                                        1:qdims%k_end) / &
        bulk_cloud_fraction(qdims%i_start:qdims%i_end,   &
                            qdims%j_start:qdims%j_end,   &
                                        1:qdims%k_end)

    ENDWHERE

  END IF

  ! Cl
  IF (L_IAU_ScaleCloud) THEN

    ! Scale Cl increments to be physical values
    ALLOCATE (Cl_Inc(qdims%i_start:qdims%i_end,  &
                     qdims%j_start:qdims%j_end,  &
                                 1:qdims%k_end))
    ALLOCATE (Cf_Inc(qdims%i_start:qdims%i_end,  &
                     qdims%j_start:qdims%j_end,  &
                                 1:qdims%k_end))
    WHERE(RESHAPE(cloud_fraction_liquid_0,             &
                  (/row_length,rows,wet_levels/)) >    &
      cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
      .AND. RESHAPE(cloud_fraction_liquid_plus -       &
                    cloud_fraction_liquid_0,           &
                    (/row_length, rows, wet_levels/)) < 0.0)
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, wet_levels/))
      WHERE (cloud_fraction_liquid                 &
                       (qdims%i_start:qdims%i_end, &
                        qdims%j_start:qdims%j_end, &
                                    1:qdims%k_end) > 0.0)
        initialScaling = - cloud_fraction_liquid       &
                           (qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
                         * (1.0 - exp(Cl_Inc /         &
                            cloud_fraction_liquid      &
                          (qdims%i_start:qdims%i_end,  &
                           qdims%j_start:qdims%j_end,  &
                                       1:qdims%k_end)))
        scalingCoeff =                                      &
            - cloud_fraction_liquid                         &
                       (qdims%i_start:qdims%i_end,          &
                        qdims%j_start:qdims%j_end,          &
                                    1:qdims%k_end) *        &
              (1.0 - exp(-RESHAPE(cloud_fraction_liquid_0,  &
                          (/row_length,rows,wet_levels/)) / &
              cloud_fraction_liquid                         &
                       (qdims%i_start:qdims%i_end,          &
                        qdims%j_start:qdims%j_end,          &
                                    1:qdims%k_end)))
        ! The following calculation should remain >0 so
        ! no error trap needed
        scalingCoeff =                                       &
            (scalingCoeff +                                  &
             cloud_fraction_liquid                           &
                       (qdims%i_start:qdims%i_end,           &
                        qdims%j_start:qdims%j_end,           &
                                    1:qdims%k_end)) /        &
            (scalingCoeff + RESHAPE(cloud_fraction_liquid_0, &
                            (/row_length,rows,wet_levels/)))
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cl_Inc = initialScaling - (initialScaling - Cl_Inc) * &
               scalingCoeff
    ELSEWHERE(RESHAPE(cloud_fraction_liquid_0,          &
                      (/row_length,rows,wet_levels/)) < &
              cloud_fraction_liquid                     &
              (qdims%i_start:qdims%i_end,               &
               qdims%j_start:qdims%j_end,               &
                           1:qdims%k_end) .AND.         &
              RESHAPE(cloud_fraction_liquid_plus -      &
                      cloud_fraction_liquid_0,          &
                      (/row_length, rows, wet_levels/)) > 0.0)
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, wet_levels/))
      WHERE ((1.0 - cloud_fraction_liquid       &
                    (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end)) > 0.0)
        initialScaling = (1.0 - cloud_fraction_liquid   &
                          (qdims%i_start:qdims%i_end,   &
                           qdims%j_start:qdims%j_end,   &
                                       1:qdims%k_end))* &
                         (1.0 - exp(-Cl_Inc / (1.0 -    &
                           cloud_fraction_liquid        &
                         (qdims%i_start:qdims%i_end,    &
                          qdims%j_start:qdims%j_end,    &
                                      1:qdims%k_end))))
        scalingCoeff = (1.0 - cloud_fraction_liquid              &
            (qdims%i_start:qdims%i_end,                          &
             qdims%j_start:qdims%j_end,                          &
                         1:qdims%k_end)) *                       &
            (1.0 - exp(-(1.0 - RESHAPE(cloud_fraction_liquid_0,  &
                               (/row_length,rows,wet_levels/)))/ &
            (1.0 - cloud_fraction_liquid                         &
                   (qdims%i_start:qdims%i_end,                   &
                    qdims%j_start:qdims%j_end,                   &
                                1:qdims%k_end))))
        scalingCoeff =                                &
            ((1.0 - cloud_fraction_liquid             &
                    (qdims%i_start:qdims%i_end,       &
                     qdims%j_start:qdims%j_end,       &
                                 1:qdims%k_end)) -    &
             scalingCoeff) / ((1.0 -                  &
             RESHAPE(cloud_fraction_liquid_0,         &
                     (/row_length,rows,wet_levels/))) &
             - scalingCoeff)
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cl_Inc = initialScaling + (Cl_Inc - initialScaling) * &
               scalingCoeff
    ELSEWHERE
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, wet_levels/))
    END WHERE
    ! Reuse cloud_fraction_liquid_0 to store initial values for
    ! use in Cb incrementing
    cloud_fraction_liquid_0 = RESHAPE(cloud_fraction_liquid &
                         (qdims%i_start:qdims%i_end,                  &
                          qdims%j_start:qdims%j_end,                  &
                                      1:qdims%k_end), (/field_size/))
    cloud_fraction_liquid(qdims%i_start:qdims%i_end,     &
                          qdims%j_start:qdims%j_end,     &
                                      1:qdims%k_end) =   &
        cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
        + Cl_Inc

  ELSE

    cloud_fraction_liquid(qdims%i_start:qdims%i_end,     &
                          qdims%j_start:qdims%j_end,     &
                                      1:qdims%k_end) =   &
        cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
        + RESHAPE(cloud_fraction_liquid_plus -           &
                  cloud_fraction_liquid_0,               &
                  (/row_length, rows, wet_levels/))

  END IF

  ! Set Cl <= 1.0
  WHERE (                                              &
      cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
      > 1.0)

    cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                          qdims%j_start:qdims%j_end, &
                                      1:qdims%k_end) &
        = 1.0

  ENDWHERE

  ! Set qcl & Cl >= 0.0
  WHERE (                                              &
      cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
      <= 0.0                                           &
      .OR. qCL(qdims%i_start:qdims%i_end,              &
               qdims%j_start:qdims%j_end,              &
                           1:qdims%k_end) <= CC_tol)

    qCL(qdims%i_start:qdims%i_end, &
        qdims%j_start:qdims%j_end, &
                    1:qdims%k_end) = 0.0
    cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                          qdims%j_start:qdims%j_end, &
                                      1:qdims%k_end) &
        = 0.0

  ENDWHERE

  IF (L_IAU_IncrementIce) THEN

    ! Cf
    IF (L_IAU_ScaleCloud) THEN

      ! Scale Cf increments to be physical values
      WHERE(RESHAPE(cloud_fraction_frozen_0,             &
                    (/row_length,rows,wet_levels/)) >    &
        cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
        .AND. RESHAPE(cloud_fraction_frozen_plus -       &
                      cloud_fraction_frozen_0,           &
                      (/row_length,rows,wet_levels/)) < 0.0)
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,wet_levels/))
        WHERE (cloud_fraction_frozen       &
               (qdims%i_start:qdims%i_end, &
                qdims%j_start:qdims%j_end, &
                            1:qdims%k_end) > 0.0)
          initialScaling = - cloud_fraction_frozen       &
                             (qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
                           * (1.0 - exp(Cf_Inc /         &
                              cloud_fraction_frozen      &
                            (qdims%i_start:qdims%i_end,  &
                             qdims%j_start:qdims%j_end,  &
                                         1:qdims%k_end)))
          scalingCoeff =                                      &
              - cloud_fraction_frozen                         &
                (qdims%i_start:qdims%i_end,                   &
                 qdims%j_start:qdims%j_end,                   &
                             1:qdims%k_end) *                 &
                (1.0 - exp(-RESHAPE(cloud_fraction_frozen_0,  &
                            (/row_length,rows,wet_levels/)) / &
                cloud_fraction_frozen                         &
                (qdims%i_start:qdims%i_end,                   &
                 qdims%j_start:qdims%j_end,                   &
                             1:qdims%k_end)))
          ! The following calculation should remain >0 so
          ! no error trap needed
          scalingCoeff =                                       &
              (scalingCoeff +                                  &
               cloud_fraction_frozen                           &
               (qdims%i_start:qdims%i_end,                     &
                qdims%j_start:qdims%j_end,                     &
                            1:qdims%k_end)) /                  &
              (scalingCoeff + RESHAPE(cloud_fraction_frozen_0, &
                              (/row_length,rows,wet_levels/)))
        ELSEWHERE
          initialScaling = 0.0
          scalingCoeff = 0.0
        END WHERE
        Cf_Inc = initialScaling - (initialScaling - Cf_Inc) * &
                 scalingCoeff
      ELSEWHERE(RESHAPE(cloud_fraction_frozen_0,          &
                        (/row_length,rows,wet_levels/)) < &
                cloud_fraction_frozen                     &
                (qdims%i_start:qdims%i_end,               &
                 qdims%j_start:qdims%j_end,               &
                             1:qdims%k_end) .AND.         &
                RESHAPE(cloud_fraction_frozen_plus -      &
                        cloud_fraction_frozen_0,          &
                        (/row_length,rows,wet_levels/)) > 0.0)
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,wet_levels/))
        WHERE ((1.0 - cloud_fraction_frozen       &
                      (qdims%i_start:qdims%i_end, &
                       qdims%j_start:qdims%j_end, &
                                   1:qdims%k_end)) > 0.0)
          initialScaling = (1.0 - cloud_fraction_frozen &
                      (qdims%i_start:qdims%i_end,       &
                       qdims%j_start:qdims%j_end,       &
                                   1:qdims%k_end))*     &
                           (1.0 - exp(-Cf_Inc / (1.0 -  &
                             cloud_fraction_frozen      &
                      (qdims%i_start:qdims%i_end,       &
                       qdims%j_start:qdims%j_end,       &
                                   1:qdims%k_end))))
          scalingCoeff = (1.0 - cloud_fraction_frozen            &
              (qdims%i_start:qdims%i_end,                        &
               qdims%j_start:qdims%j_end,                        &
                           1:qdims%k_end)) *                     &
              (1.0 - exp(-(1.0 -RESHAPE(cloud_fraction_frozen_0, &
                           (/row_length,rows,wet_levels/))) /    &
              (1.0 - cloud_fraction_frozen                       &
                     (qdims%i_start:qdims%i_end,                 &
                      qdims%j_start:qdims%j_end,                 &
                                  1:qdims%k_end))))
          scalingCoeff =                             &
              ((1.0 - cloud_fraction_frozen          &
                      (qdims%i_start:qdims%i_end,    &
                       qdims%j_start:qdims%j_end,    &
                                   1:qdims%k_end)) - &
               scalingCoeff) / ((1.0 -          &
               RESHAPE(cloud_fraction_frozen_0, &
               (/row_length,rows,wet_levels/))) - scalingCoeff)
        ELSEWHERE
          initialScaling = 0.0
          scalingCoeff = 0.0
        END WHERE
        Cf_Inc = initialScaling + (Cf_Inc - initialScaling) * &
                 scalingCoeff
      ELSEWHERE
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,wet_levels/))
      END WHERE
      ! Reuse cloud_fraction_frozen_0 to store initial values for
      ! use in Cb incrementing
      cloud_fraction_frozen_0 = RESHAPE(cloud_fraction_frozen &
          (qdims%i_start:qdims%i_end,                         &
           qdims%j_start:qdims%j_end,                         &
                       1:qdims%k_end), (/field_size/))
      cloud_fraction_frozen(qdims%i_start:qdims%i_end,   &
                            qdims%j_start:qdims%j_end,   &
                                        1:qdims%k_end) = &
        cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
        qdims%j_start:qdims%j_end,                       &
                    1:qdims%k_end)                       &
        + Cf_Inc

    ELSE

    cloud_fraction_frozen(qdims%i_start:qdims%i_end,     &
                          qdims%j_start:qdims%j_end,     &
                                      1:qdims%k_end) =   &
        cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
        + RESHAPE(cloud_fraction_frozen_plus -           &
         cloud_fraction_frozen_0,                        &
         (/row_length, rows, wet_levels/))

    END IF

    ! Set Cf <= 1.0
    WHERE (                                            &
      cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
        > 1.0)

      cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) = 1.0

    ENDWHERE

    ! Set qcf & Cf >= 0.0
    WHERE (                                              &
        cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
        <= 0.0                                           &
        .OR. qCF(qdims%i_start:qdims%i_end,              &
                 qdims%j_start:qdims%j_end,              &
                             1:qdims%k_end) <= CC_tol)

      qCF(qdims%i_start:qdims%i_end, &
          qdims%j_start:qdims%j_end, &
                      1:qdims%k_end) = 0.0
      cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
          = 0.0

    ENDWHERE

  END IF

  ! Set other cloud variables to zero if qCL & qCF <= CC_tol
  WHERE (   qCL(qdims%i_start:qdims%i_end,           &
                qdims%j_start:qdims%j_end,           &
                            1:qdims%k_end) <= CC_tol &
      .AND. qCF(qdims%i_start:qdims%i_end,           &
                qdims%j_start:qdims%j_end,           &
                            1:qdims%k_end) <= CC_tol)

    qCF(qdims%i_start:qdims%i_end, &
        qdims%j_start:qdims%j_end, &
                    1:qdims%k_end) = 0.0
    cloud_fraction_frozen(qdims%i_start:qdims%i_end, &
                          qdims%j_start:qdims%j_end, &
                                      1:qdims%k_end) &
        = 0.0
    area_cloud_fraction(qdims%i_start:qdims%i_end, &
                        qdims%j_start:qdims%j_end, &
                                    1:qdims%k_end) &
        = 0.0
    bulk_cloud_fraction(qdims%i_start:qdims%i_end, &
                        qdims%j_start:qdims%j_end, &
                                    1:qdims%k_end) &
        = 0.0

  ENDWHERE

  IF (L_IAU_IncrementIce) THEN

    WHERE (   qCL(qdims%i_start:qdims%i_end,           &
                  qdims%j_start:qdims%j_end,           &
                              1:qdims%k_end) <= CC_tol &
        .AND. qCF(qdims%i_start:qdims%i_end,           &
                  qdims%j_start:qdims%j_end,           &
                              1:qdims%k_end) <= CC_tol)

      qCL(qdims%i_start:qdims%i_end,               &
                        qdims%j_start:qdims%j_end, &
                                    1:qdims%k_end) = 0.0
      cloud_fraction_liquid(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) = 0.0

    ENDWHERE

  END IF

  ! Apply stored area / bulk cloud fractions to output
  IF (L_IAU_ScaleCloud) THEN

    ! Scale Cb increments to be physical values
    ! Use Cb_inc and Cf_Inc as temporary storage for incremented
    ! and unincremented MIN(Cl+Cf, 1.0)
    ALLOCATE (Cb_Inc (qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))
    Cb_Inc = MIN(cloud_fraction_liquid               &
                       (qdims%i_start:qdims%i_end,   &
                        qdims%j_start:qdims%j_end,   &
                                    1:qdims%k_end) + &
                 cloud_fraction_frozen               &
                       (qdims%i_start:qdims%i_end,   &
                        qdims%j_start:qdims%j_end,   &
                                    1:qdims%k_end), 1.0)
    IF (L_IAU_IncrementIce) THEN
      Cf_Inc = RESHAPE(MIN(cloud_fraction_liquid_0 +     &
                           cloud_fraction_frozen_0,1.0), &
                       (/row_length,rows,wet_levels/))
    ELSE
      Cf_Inc = MIN(RESHAPE(cloud_fraction_liquid_0,      &
                       (/row_length,rows,wet_levels/)) + &
                       cloud_fraction_frozen             &
                       (qdims%i_start:qdims%i_end,       &
                        qdims%j_start:qdims%j_end,       &
                                    1:qdims%k_end), 1.0)
    END IF
    WHERE(Cf_Inc >                                       &
          bulk_cloud_fraction(qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
          .AND. Cb_Inc - Cf_Inc < 0.0)
      Cb_Inc = Cb_Inc - Cf_Inc
      WHERE (bulk_cloud_fraction         &
             (qdims%i_start:qdims%i_end, &
              qdims%j_start:qdims%j_end, &
                          1:qdims%k_end) > 0.0)
        initialScaling = - bulk_cloud_fraction           &
                           (qdims%i_start:qdims%i_end,   &
                            qdims%j_start:qdims%j_end,   &
                                        1:qdims%k_end) * &
                         (1.0 - exp(Cb_Inc /             &
                                    bulk_cloud_fraction  &
                           (qdims%i_start:qdims%i_end,   &
                            qdims%j_start:qdims%j_end,   &
                                        1:qdims%k_end)))
        scalingCoeff  = - bulk_cloud_fraction            &
                             (qdims%i_start:qdims%i_end, &
                              qdims%j_start:qdims%j_end, &
                                          1:qdims%k_end) &
                           * (1.0 - exp(-Cf_Inc /        &
                             bulk_cloud_fraction         &
                          (qdims%i_start:qdims%i_end,    &
                           qdims%j_start:qdims%j_end,    &
                                       1:qdims%k_end)))
        ! The following calculation should remain >0 so no error
        !trap needed
        scalingCoeff =                                    &
         (scalingCoeff +                                  &
          bulk_cloud_fraction(qdims%i_start:qdims%i_end,  &
                              qdims%j_start:qdims%j_end,  &
                                          1:qdims%k_end)) &
          / (scalingCoeff + Cf_Inc  )
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cb_Inc = initialScaling - (initialScaling - Cb_Inc) * &
                                scalingCoeff
    ELSEWHERE(Cf_Inc <                          &
              bulk_cloud_fraction               &
              (qdims%i_start:qdims%i_end,       &
               qdims%j_start:qdims%j_end,       &
                           1:qdims%k_end) .AND. &
              Cb_Inc - Cf_Inc > 0.0)
      Cb_Inc = Cb_Inc - Cf_Inc
      WHERE ((1.0 - bulk_cloud_fraction         &
                    (qdims%i_start:qdims%i_end, &
                     qdims%j_start:qdims%j_end, &
                                 1:qdims%k_end)) > 0.0)
        initialScaling = (1.0 - bulk_cloud_fraction             &
                          (qdims%i_start:qdims%i_end,           &
                           qdims%j_start:qdims%j_end,           &
                                       1:qdims%k_end))          &
                       * (1.0 - exp(-Cb_Inc / (1.0 -            &
                                            bulk_cloud_fraction &
                         (qdims%i_start:qdims%i_end,            &
                          qdims%j_start:qdims%j_end,            &
                                      1:qdims%k_end))))
        scalingCoeff = (1.0 - bulk_cloud_fraction             &
                          (qdims%i_start:qdims%i_end,         &
                           qdims%j_start:qdims%j_end,         &
                                       1:qdims%k_end))        &
                        * (1.0 - exp(-(1.0 - Cf_Inc) / (1.0 - &
                          bulk_cloud_fraction                 &
                         (qdims%i_start:qdims%i_end,          &
                          qdims%j_start:qdims%j_end,          &
                                      1:qdims%k_end))))
        scalingCoeff =                             &
            ((1.0 - bulk_cloud_fraction            &
                    (qdims%i_start:qdims%i_end,    &
                     qdims%j_start:qdims%j_end,    &
                                 1:qdims%k_end)) - &
             scalingCoeff) /                       &
            ((1.0 - Cf_Inc  ) - scalingCoeff)
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cb_Inc = initialScaling + (Cb_Inc - initialScaling) * &
                                scalingCoeff
    ELSEWHERE
      Cb_Inc = Cb_Inc - Cf_Inc
    END WHERE
    bulk_cloud_fraction(qdims%i_start:qdims%i_end,     &
                        qdims%j_start:qdims%j_end,     &
                                    1:qdims%k_end) =   &
        bulk_cloud_fraction(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end) &
        + Cb_Inc

  ELSE

    bulk_cloud_fraction(qdims%i_start:qdims%i_end,      &
                        qdims%j_start:qdims%j_end,      &
                                    1:qdims%k_end) =    &
       MIN(1.0,MAX(0.0,                                 &
        cloud_fraction_liquid                           &
                            (qdims%i_start:qdims%i_end, &
                             qdims%j_start:qdims%j_end, &
                                         1:qdims%k_end) &
        +                                               &
        cloud_fraction_frozen                           &
                            (qdims%i_start:qdims%i_end, &
                             qdims%j_start:qdims%j_end, &
                                         1:qdims%k_end)))

  END IF

  IF (.NOT. L_CLD_AREA) THEN

    area_cloud_fraction  (qdims%i_start:qdims%i_end,   &
                          qdims%j_start:qdims%j_end,   &
                                      1:qdims%k_end) = &
      bulk_cloud_fraction(qdims%i_start:qdims%i_end,   &
                          qdims%j_start:qdims%j_end,   &
                                      1:qdims%k_end)

  ELSE

    area_cloud_fraction(qdims%i_start:qdims%i_end,     &
                        qdims%j_start:qdims%j_end,     &
                                    1:qdims%k_end) =   &
        MIN(1.0,                                       &
        AreaDivBulk *                                  &
        bulk_cloud_fraction(qdims%i_start:qdims%i_end, &
                            qdims%j_start:qdims%j_end, &
                                        1:qdims%k_end))

  END IF

  ! q
  IF (L_IAU_IncrementIce) THEN

    IF (L_IAU_ScaleCloud) THEN

      q(qdims%i_start:qdims%i_end,                    &
        qdims%j_start:qdims%j_end,                    &
                    1:qdims%k_end) =                  &
          q(qdims%i_start:qdims%i_end,                &
            qdims%j_start:qdims%j_end,                &
                        1:qdims%k_end)     +          &
          (qT_plus - qT)                            - &
          qcl_Inc -                                   &
          RESHAPE(qCF_plus - qCF_0,                   &
                  (/row_length, rows, wet_levels/))

    ELSE

      q(qdims%i_start:qdims%i_end,                    &
        qdims%j_start:qdims%j_end,                    &
                    1:qdims%k_end) =                  &
          q(qdims%i_start:qdims%i_end,                &
            qdims%j_start:qdims%j_end,                &
                        1:qdims%k_end)     +          &
          (qT_plus - qT)                            - &
          RESHAPE(qCL_plus - qCL_0,                   &
                  (/row_length, rows, wet_levels/)) - &
          RESHAPE(qCF_plus - qCF_0,                   &
                  (/row_length, rows, wet_levels/))

    END IF

  ELSE

    IF (L_IAU_ScaleCloud) THEN

      q(qdims%i_start:qdims%i_end,                    &
        qdims%j_start:qdims%j_end,                    &
                    1:qdims%k_end) =                  &
          q(qdims%i_start:qdims%i_end,                &
            qdims%j_start:qdims%j_end,                &
                        1:qdims%k_end)     +          &
          (qT_plus - qT)                            - &
          qcl_Inc

    ELSE

      q(qdims%i_start:qdims%i_end,                    &
        qdims%j_start:qdims%j_end,                    &
                    1:qdims%k_end) =                  &
          q(qdims%i_start:qdims%i_end,                &
            qdims%j_start:qdims%j_end,                &
                        1:qdims%k_end)     +          &
          (qT_plus - qT)                            - &
          RESHAPE(qCL_plus - qCL_0,                   &
                  (/row_length, rows, wet_levels/))

    END IF

  END IF

  ! Apply lower limit
  WHERE (q(qdims%i_start:qdims%i_end, &
           qdims%j_start:qdims%j_end, &
                       1:qdims%k_end) < q_min)

    q(qdims%i_start:qdims%i_end, &
      qdims%j_start:qdims%j_end, &
                  1:qdims%k_end) = q_min

  ENDWHERE

  !-----------------------------------------------------------------------------
  ! [7.6.3]: Deallocate work arrays.
  !-----------------------------------------------------------------------------

  DEALLOCATE (AreaDivBulk)
  DEALLOCATE (cloud_fraction_liquid_0)
  DEALLOCATE (cloud_fraction_liquid_plus)
  DEALLOCATE (p_theta_levels_old)
  DEALLOCATE (qCL_0)
  DEALLOCATE (qCL_plus)
  DEALLOCATE (qSatW)
  DEALLOCATE (Var_qT)
  DEALLOCATE (VarTemp)

  IF (L_IAU_IncrementIce) THEN
    DEALLOCATE (cloud_fraction_frozen_0)
    DEALLOCATE (cloud_fraction_frozen_plus)
    DEALLOCATE (qCF_0)
    DEALLOCATE (qCF_max)
    DEALLOCATE (qCF_plus)
    DEALLOCATE (Var_BGT)
  END IF

  IF (L_IAU_ScaleCloud) THEN
    DEALLOCATE (initialScaling)
    DEALLOCATE (scalingCoeff)
    DEALLOCATE (Cl_Inc)
    DEALLOCATE (Cf_Inc)
    DEALLOCATE (Cb_Inc)
  END IF

  DEALLOCATE (qcl_Inc)

END IF ! (qTIncsThisTS .AND. .NOT.L_IAU_Add_qT_prime_to_q)

IF (qTIncsThisTS) THEN

  ! Update q using qT_prime when Var_DiagCloud is bypassed:
  IF (L_IAU_Add_qT_prime_to_q) THEN
    q(qdims%i_start:qdims%i_end,       &
      qdims%j_start:qdims%j_end,       &
      1:qdims%k_end) =                 &
    qT_plus(qdims%i_start:qdims%i_end, &
            qdims%j_start:qdims%j_end, &
            1:qdims%k_end)  -          &
    qcl(qdims%i_start:qdims%i_end,     &
        qdims%j_start:qdims%j_end,     &
        1:qdims%k_end)

    IF (L_IAU_IncrementIce .AND. L_IAU_ApplyQTCorrections) THEN
      q(qdims%i_start:qdims%i_end, &
        qdims%j_start:qdims%j_end, &
        1:qdims%k_end) =           &
      q(qdims%i_start:qdims%i_end, &
        qdims%j_start:qdims%j_end, &
        1:qdims%k_end) -           &
      qcf(qdims%i_start:qdims%i_end, &
          qdims%j_start:qdims%j_end, &
          1:qdims%k_end)
    END IF
    
  END IF

  DEALLOCATE (qT)
  DEALLOCATE (qT_plus)

ENDIF

!-------------------------------------------------------------------------------
! [7.7]: PC2 updates to water vapour, cloud and temperature.
!-------------------------------------------------------------------------------

IF (.NOT.qTIncsThisTS .AND. &
   L_pc2 .AND. (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff)) THEN

  ALLOCATE (delta_q  (qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))
  ALLOCATE (delta_qcl(qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))
  ALLOCATE (delta_qcf(qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))
  ALLOCATE (delta_p  (qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))
  ALLOCATE (delta_t  (qdims%i_start:qdims%i_end, &
                      qdims%j_start:qdims%j_end, &
                                  1:qdims%k_end))

  ! Calculate increments to vapour, cloud, pressure and temperature.
  ! (The _work arrays currently contain values on entry to the routine.)
  DO k =             1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i=  qdims%i_start, qdims%i_end
        delta_q  (i,j,k) = q  (i,j,k)                - q_work  (i,j,k)
        delta_qcl(i,j,k) = qcl(i,j,k)                - qcl_work(i,j,k)
        delta_qcf(i,j,k) = qcf(i,j,k)                - qcf_in  (i,j,k)
        delta_p  (i,j,k) = p_theta_levels    (i,j,k) - p_work  (i,j,k)
        delta_t  (i,j,k) = exner_theta_levels(i,j,k) &
                         * theta(i,j,k)              - t_work  (i,j,k)
      END DO
    END DO
  END DO

  CALL pc2_assim ( REAL(TS_len_secs),                            & ! in
                   l_pc2_cfl,                                    & ! in
                   l_pc2_cff,                                    & ! in
                   l_mixing_ratio,                               & ! in
                   t_work,                                       & ! inout
                   bulk_cloud_fraction  (1:row_length,1:rows,:), & ! inout
                   cloud_fraction_liquid(1:row_length,1:rows,:), & ! inout
                   cloud_fraction_frozen(1:row_length,1:rows,:), & ! inout
                   q_work,                                       & ! inout
                   qcl_work,                                     & ! inout
                   qcf_in,                                       & ! in
                   p_work,                                       & ! in
                   delta_t,                                      & ! in
                   delta_q,                                      & ! in
                   delta_qcl,                                    & ! in
                   delta_qcf,                                    & ! in
                   delta_p                                       & ! in
                  )

  ! q_work, qcl_work and t_work have now been updated by the forcing and
  ! condensation terms together. cloud_fraction_liquid (and _frozen and _bulk)
  ! has been updated by the condensation. qcf_work and p_work are not updated.

  ! Now copy q_work, qcl_work and t_work variables back into the INOUT
  ! variables if we haven't a confident direct estimate of condensation from
  ! another source:
  IF (L_pc2_cond) THEN
    DO k =             1, qdims%k_end
      DO j = qdims%j_start, qdims%j_end
        DO i= qdims%i_start, qdims%i_end
          q    (i,j,k) = q_work  (i,j,k)
          qcl  (i,j,k) = qcl_work(i,j,k)
          ! Remember the INOUT variable is theta, not temperature:
          theta(i,j,k) = t_work  (i,j,k) / exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DEALLOCATE (q_work)
  DEALLOCATE (qcl_work)
  DEALLOCATE (qcf_in)
  DEALLOCATE (p_work)
  DEALLOCATE (t_work)
  DEALLOCATE (delta_q)
  DEALLOCATE (delta_qcl)
  DEALLOCATE (delta_qcf)
  DEALLOCATE (delta_p)
  DEALLOCATE (delta_t)

END IF ! (.NOT.qTIncsThisTS .AND. &
       !  L_pc2 .AND. (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff))

!-------------------------------------------------------------------------------
! [7.8]: TSoil(1), TStar and TStar_tile.
!-------------------------------------------------------------------------------

IF (L_IAU_IncTStar) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE (6,*) 'IAU: Adding t1 incs to TStar and TSoil(1)'
  END IF

  ! Get increment to level-one temperature:
  addr = 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t1_inc(addr) = exner_theta_levels(i,j,1) &
                   * theta             (i,j,1) &
                   - t1_in(i,j)
      addr = addr + 1
    END DO
  END DO

  ! Add increments to TSoil(1) and TStar where
  ! snow_depth <= 0.05 kg/m**2 or explicitly requested:
  DO i = 1, land_field

    IF ( snow_depth(land_index(i)) <= 0.05 ) THEN
      TSoil(i,1) = TSoil(i,1) + t1_inc(land_index(i))
      TStar(land_index(i)) =  TStar(land_index(i)) &
                           + t1_inc(land_index(i))
    ELSE IF (L_IAU_IncTSurf_Snow) THEN
!     Limit the increment not to allow melting.
      IF (TSoil(i,1) < tm) TSoil(i,1) = MIN(TSoil(i,1) &
                           + t1_inc(land_index(i)), tm)
      IF (TStar(land_index(i)) < tm) TStar(land_index(i)) = &
                           MIN(TStar(land_index(i))         &
                           + t1_inc(land_index(i)), tm)
    END IF

  END DO

  ! If using a tile scheme, also add increments to TStar_tile:
  IF (L_IAU_IncTStar_tile) THEN

    DO tile_num = 1, ntiles
      IF ( (tile_num /= lake) .OR. L_IAU_IncTLake ) THEN
        DO i = 1, land_field

          IF ( snow_depth(land_index(i)) <= 0.05 ) THEN
            TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &
                                   + t1_inc(land_index(i))
          ELSE IF (L_IAU_IncTSurf_Snow) THEN
            IF (TStar_tile(i,tile_num) < tm) TStar_tile(i,tile_num) = &
                             MIN(TStar_tile(i,tile_num)               &
                             + t1_inc(land_index(i)), tm)
          END IF

        END DO
      END IF
    END DO

  END IF

END IF ! (L_IAU_IncTStar)

!-------------------------------------------------------------------------------
! [8]: Remove supersaturation wrt water.
!-------------------------------------------------------------------------------

IF (L_IAU_RemoveSS) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE (6,*) 'IAU: Removing supersaturation wrt water'
  END IF

  DO k = 1, qdims%k_end

    ! Get temperature and pressure:
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ! Use t1_in work array to hold temp on this level:
        t1_in(i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
        p_tmp(i,j) = p_theta_levels    (i,j,k)
      END DO
    END DO

    ! Calculate saturated specific humidity wrt water:
    CALL QSAT_WAT (q_sat_wat, t1_in, p_tmp, row_length*rows)

    ! Limit q to q_sat, or a fraction of qsat:
    IF (k > IAU_qsatScale_maxLev) THEN
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q(i,j,k) = MIN(q(i,j,k), q_sat_wat(i,j)) 
        END DO
      END DO
    ELSE
      DO j = qdims%j_start, qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          q(i,j,k) = MIN(q(i,j,k), IAU_qsatScale * q_sat_wat(i,j)) 
        END DO
      END DO
    END IF

  END DO

END IF ! L_IAU_RemoveSS

!-------------------------------------------------------------------------------
! [9]: Apply bounds to tropospheric and non-tropospheric humidities.
!-------------------------------------------------------------------------------

IF (CallQLimitsThisTS) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE (6,*) 'IAU: Calling QLimits'
  END IF

  ! Need to update halos ready for calculation of pv_at_theta:
  i_field = 1
  fields_to_swap(i_field) % field      => u(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_u
  fields_to_swap(i_field) % levels     =  udims%k_end - udims%k_start + 1 
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => v(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_v
  fields_to_swap(i_field) % levels     =  vdims%k_end - vdims%k_start + 1 
  fields_to_swap(i_field) % rows       =  n_rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => theta(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_p
  fields_to_swap(i_field) % levels     =  tdims%k_end - tdims%k_start + 1 
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => rho(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_p
  fields_to_swap(i_field) % levels     =  pdims%k_end - pdims%k_start + 1 
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .FALSE.
         
  CALL swap_bounds_mv (fields_to_swap, i_field, row_length, offx, offy)
 
  ALLOCATE (pv_at_theta(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end))

  CALL Calc_PV_at_theta ( u, v, theta, rho,             & ! in
                          r_theta_levels, r_rho_levels, & ! in
                          r_at_u, r_at_v,               & ! in
                          sec_v_latitude,               & ! in
                          tan_v_latitude,               & ! in
                          sec_theta_latitude,           & ! in
                          f3_at_v,                      & ! in
                          delta_lambda, delta_phi,      & ! in
                          Model_domain,                 & ! in
                          pv_at_theta )                   ! out

  CALL QLimits ( L_IAU_QLimitsDiags,    & ! in
                 L_IAU_RmNonTropQIncs,  & ! in
                 IAU_trop_min_RH,       & ! in
                 IAU_nonTrop_max_q,     & ! in
                 IAU_nonTrop_min_q,     & ! in
                 IAU_nonTrop_max_RH,    & ! in
                 IAU_trop_min_p,        & ! in
                 IAU_trop_max_PV,       & ! in
                 IAU_nonTrop_max_p,     & ! in
                 q_in,                  & ! in
                 pv_at_theta,           & ! in
                 theta,                 & ! in
                 p,                     & ! in
                 p_theta_levels,        & ! in
                 exner_theta_levels,    & ! in
                 q,                     & ! inout
                 qCL,                   & ! inout
                 qCF,                   & ! inout
                 area_cloud_fraction,   & ! inout
                 bulk_cloud_fraction,   & ! inout
                 cloud_fraction_liquid, & ! inout
                 cloud_fraction_frozen )  ! inout

  DEALLOCATE (pv_at_theta)
  DEALLOCATE (q_in)

END IF ! (CallQLimitsThisTS)

DEALLOCATE(IAUIncFld)

!-------------------------------------------------------------------------------
! [10]: Loop through IAU structures and deallocate arrays that are no longer
!       needed.
!-------------------------------------------------------------------------------

DO IncNum = 1, Num_IAU_incs

  IF (STEPim(a_im) == IAU_incs(IncNum) % InsertionEndTS .AND. &
                      IAU_incs(IncNum) % FieldsToAdd) THEN

    DEALLOCATE (IAU_incs(IncNum) % Lookup)
    DEALLOCATE (IAU_incs(IncNum) % Weights)

    IF (IAU_incs(IncNum) % FieldsStored) THEN
      IF (IAU_incs(IncNum) % StorePacked) THEN
        DEALLOCATE (IAU_incs(IncNum) % flds32bit)
      ELSE
        DEALLOCATE (IAU_incs(IncNum) % flds)
      END IF
      IAU_incs(IncNum) % FieldsStored = .FALSE.
    END IF

  END IF

END DO

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (6,*) 'IAU: Exiting IAU. Time since basis time (hr:mm:ss): ', TimeStr
END IF

IF (lhook) CALL dr_hook('IAU',zhook_out,zhook_handle)
RETURN

END SUBROUTINE IAU

