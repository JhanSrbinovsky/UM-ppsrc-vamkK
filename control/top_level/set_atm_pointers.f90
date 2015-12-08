! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SET_ATM_POINTERS ---------------------------------------
!LL
!LL  Set pointers for primary atmosphere fields
!LL
!LL Programming Standard: Unified Model DP NO. 3, Version 8.3
!LL
!LL  Purpose:   Sets integer pointers to atmospheric
!LL             variables from STASHIN addresses.
!LL
!LL  External documentation: UMDP NO. C4 
!LL
!LLEND-------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

      SUBROUTINE SET_ATM_POINTERS(                                      &
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ARGPTRA start
      ! Pointers for ATMOSPHERE model variables. Configuration dependent.

      ! Addresses in D1 array of primary variables held in primary and
      ! secondary space

      ! Data variables stored in primary space.
        JU, JV, JW, JRHO, JTHETA, JQ, JQCL, JQCF,                       &
        JEXNER_RHO_LEVELS, JU_ADV, JV_ADV, JW_ADV,                      &
      ! Data variables stored in secondary space.
        JP, JP_THETA_LEVELS,  JEXNER_THETA_LEVELS,                      &
      ! Cloud Fields
        JCCA, JCF_AREA, JCF_BULK, JCF_LIQUID, JCF_FROZEN,               &
      ! Soil Ancillary Fields
        J_DEEP_SOIL_TEMP,  JSMCL, JSTHU, JSTHF,                         &
      ! Radiation increments
        JSW_INCS, JLW_INCS,                                             &
      ! Ozone
        JOZONE,                                                         &
      ! Tracers and Aerosols
        JTRACER, JMURK, JMURK_SOURCE,                                   &
        JSO2, JDMS, JSO4_AITKEN, JSO4_ACCU, JSO4_DISS, JH2O2,           &
        JNH3, JSOOT_NEW, JSOOT_AGD, JSOOT_CLD, JSO2_NATEM,              &
        JOH, JHO2, JH2O2_LIMIT,JO3_CHEM, JHadCM2_SO4, JCO2,             &
      ! User Ancillary fields
        JUSER_MULT1, JUSER_MULT2, JUSER_MULT3, JUSER_MULT4,             &
        JUSER_MULT5, JUSER_MULT6, JUSER_MULT7, JUSER_MULT8,             &
        JUSER_MULT9, JUSER_MULT10, JUSER_MULT11, JUSER_MULT12,          &
        JUSER_MULT13, JUSER_MULT14, JUSER_MULT15, JUSER_MULT16,         &
        JUSER_MULT17, JUSER_MULT18, JUSER_MULT19, JUSER_MULT20,         &
      ! Lateral Boundary Conditions and tendencies
        JOROG_LBC,JU_LBC,JU_LBC_TEND,JV_LBC,JV_LBC_TEND,JW_LBC,         &
        JW_LBC_TEND,JRHO_LBC,JRHO_LBC_TEND,JTHETA_LBC,JTHETA_LBC_TEND,  &
        JQ_LBC,JQ_LBC_TEND,JQCL_LBC, JQCL_LBC_TEND,JQCF_LBC,            &
        JQCF_LBC_TEND,JEXNER_LBC,JEXNER_LBC_TEND,JU_ADV_LBC,            &
        JU_ADV_LBC_TEND,JV_ADV_LBC,JV_ADV_LBC_TEND,JW_ADV_LBC,          &
        JW_ADV_LBC_TEND,JTRACER_LBC,JTRACER_LBC_TEND,                   &
        JTR_UKCA_LBC,JTR_UKCA_LBC_TEND,                                 &
        JSO2_LBC,JSO2_LBC_TEND,JDMS_LBC,JDMS_LBC_TEND,JSO4_AITKEN_LBC,  &
        JSO4_AITKEN_LBC_TEND,JSO4_ACCU_LBC,JSO4_ACCU_LBC_TEND,          &
        JSO4_DISS_LBC,JSO4_DISS_LBC_TEND,                               &
        JNH3_LBC,JNH3_LBC_TEND,JSOOT_NEW_LBC,JSOOT_NEW_LBC_TEND,        &
        JSOOT_AGD_LBC,JSOOT_AGD_LBC_TEND,JSOOT_CLD_LBC,                 &
        JSOOT_CLD_LBC_TEND,                                             &
! Tropopause-based Ozone
        JTPPSOZONE,                                                     &
      ! Biomass aerosol
        JBMASS_NEW, JBMASS_AGD, JBMASS_CLD,                             &
        JBMASS_NEW_LBC, JBMASS_NEW_LBC_TEND,                            &
        JBMASS_AGD_LBC, JBMASS_AGD_LBC_TEND,                            &
        JBMASS_CLD_LBC, JBMASS_CLD_LBC_TEND,                            &
      ! Additional microphysics fields and lbcs
        JQCF2,JQRAIN,JQGRAUP,JQCF2_LBC,JQCF2_LBC_TEND,JQRAIN_LBC,       &
        JQRAIN_LBC_TEND,JQGRAUP_LBC,JQGRAUP_LBC_TEND,JDUST_DIV1,        &
        JDUST_DIV2,JDUST_DIV3,                                          &
        JDUST_DIV4,JDUST_DIV5,JDUST_DIV6,                               &
        JDUST_DIV1_LBC,JDUST_DIV1_LBC_TEND,JDUST_DIV2_LBC,              &
        JDUST_DIV2_LBC_TEND, JDUST_DIV3_LBC,JDUST_DIV3_LBC_TEND,        &
        JDUST_DIV4_LBC,JDUST_DIV4_LBC_TEND,JDUST_DIV5_LBC,              &
        JDUST_DIV5_LBC_TEND,JDUST_DIV6_LBC,JDUST_DIV6_LBC_TEND,         &
        JCF_BULK_LBC,JCF_BULK_LBC_TEND,JCF_LIQUID_LBC,                  &
        JCF_LIQUID_LBC_TEND,JCF_FROZEN_LBC,JCF_FROZEN_LBC_TEND,         &
! Pointer for direct PAR flux
        JDIRPAR,                                                        &
        JTR_UKCA, JMURK_LBC, JMURK_LBC_TEND,                            &
! Pointers for UKCA oxidant fields
        JOH_UKCA, JHO2_UKCA, JH2O2_UKCA, JO3_UKCA,                      &
! Convective Cloud Fields
        JLCBASE, JCCW_RAD,                                              &
! Ozone tracer and cariolle parameters
        JOZONE_TRACER,JO3_PROD_LOSS,JO3_P_L_VMR,JO3_VMR,JO3_P_L_TEMP,   &
        JO3_TEMP,JO3_P_L_COLO3,JO3_COLO3,                               &
! Pointers for Aerosol climatologies
        JARCLBIOG_BG, JARCLBIOM_FR, JARCLBIOM_AG, JARCLBIOM_IC,         &
        JARCLBLCK_FR, JARCLBLCK_AG, JARCLSSLT_FI, JARCLSSLT_JT,         &
        JARCLSULP_AC, JARCLSULP_AK, JARCLSULP_DI, JARCLDUST_B1,         &
        JARCLDUST_B2, JARCLDUST_B3, JARCLDUST_B4, JARCLDUST_B5,         & 
        JARCLDUST_B6, JARCLOCFF_FR, JARCLOCFF_AG, JARCLOCFF_IC,         &
        JARCLDLTA_DL,                                                   &
! Fossil-fuel organic carbon aerosol
        JOCFF_NEW, JOCFF_AGD, JOCFF_CLD,                                &
        JOCFF_NEW_LBC,JOCFF_NEW_LBC_TEND,JOCFF_AGD_LBC,                 &
        JOCFF_AGD_LBC_TEND,JOCFF_CLD_LBC,JOCFF_CLD_LBC_TEND,            &
! Ammonium nitrate aerosol
        JHNO3_UKCA, JNITR_ACC, JNITR_DISS,                              &
        JNITR_ACC_LBC, JNITR_ACC_LBC_TEND, JNITR_DISS_LBC,              & 
        JNITR_DISS_LBC_TEND,                                            &
! TKE based turbulence scheme
        JE_TRB, JTSQ_TRB,                                               &
        JQSQ_TRB, JCOV_TRB, JZHPAR_SHCU,                                &
! ENDGame
        JDRYRHO,JETADOT,JTHETAV,JPSIWS,JPSIWL,JMV,JMCL,JMCF,JMCF2,      &
        JMRAIN,JMGRAUP,JEXNERSURF,                                      &
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
        JTSOIL_TILE, JSMCL_TILE, JSTHF_TILE, JSNOW_DEPTH3L,             &
        JSNOW_MASS3L, JSNOW_TMP3L, JSNOW_RHO3L, JSNOW_RHO1L, JSNOW_AGE, & 
        JSNOW_FLG3L,                                                    &
!}CABLE         
! ARGPTRA end
                        ICODE,CMESSAGE)

      use cable_data_mod, only : cable_set_atm_pointers
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ukca_d1_defs, ONLY: ukca_item_sulpc, ukca_sect
      USE clmchfcg_scenario_mod, ONLY: nsulpat
      USE atm_fields_bounds_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE um_input_control_mod, ONLY:                                  &
           l_sulpc_online_oxidants
      USE rad_input_mod, ONLY: lexpand_ozone, lexpand_tpps_ozone
      USE lookup_addresses
      USE cv_run_mod, ONLY: l_3d_cca, l_ccrad
      USE ukca_option_mod, ONLY: l_ukca
      USE Submodel_Mod


      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE
!L
!*L Arguments
!L
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
! TYPDUMA needs TYPSIZE included first
!L --------------- Dump headers (atmosphere)-------------
      INTEGER :: A_FIXHD(LEN_FIXHD)    ! fixed length header
      INTEGER :: A_INTHD(A_LEN_INTHD)  ! integer header
      INTEGER :: A_CFI1(A_LEN_CFI1+1)  ! compress field index
      INTEGER :: A_CFI2(A_LEN_CFI2+1)  ! compress field index
      INTEGER :: A_CFI3(A_LEN_CFI3+1)  ! compress field index

      REAL::A_REALHD(A_LEN_REALHD)                    ! real header
      REAL::A_LEVDEPC(A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1)! level  dep const
      REAL::A_ROWDEPC(A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1)! row    dep const
      REAL::A_COLDEPC(A_LEN1_COLDEPC*A_LEN2_COLDEPC+1)! column dep const
      REAL::A_FLDDEPC(A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1)! field  dep const
      REAL::A_EXTCNST(A_LEN_EXTCNST+1)                ! extra constants
      REAL::A_DUMPHIST(LEN_DUMPHIST+1)                ! temp hist file

      ! Meaningful parameter names for integer constants header:
! ----------------------- include file: IHEADAPM -----------------------
! Description: Meaningful parameter names to index A_INTHD array in
!              atmosphere dump, ie INTEGER CONSTANTS, and reduce magic
!              numbers in code.
!
      INTEGER,PARAMETER:: ih_a_step          = 1  ! Timestep no.
      INTEGER,PARAMETER:: ih_rowlength       = 6  ! No. of points E-W
      INTEGER,PARAMETER:: ih_rows            = 7  ! No. of points N-S

      ! No. of model levels (0=surface)
      INTEGER,PARAMETER:: ih_model_levels    = 8

      ! No. of model levels with moisture
      INTEGER,PARAMETER:: ih_wet_levels      = 9

      ! No. of deep soil temperature levels
      INTEGER,PARAMETER:: ih_soilT_levels    = 10

      INTEGER,PARAMETER:: ih_cloud_levels    = 11 ! No. of cloud levels
      INTEGER,PARAMETER:: ih_tracer_levels   = 12 ! No. of tracer levels

      ! No. of boundary layer levels
      INTEGER,PARAMETER:: ih_boundary_levels = 13
      INTEGER,PARAMETER:: ih_N_types         = 15 ! No. of field types

       ! Height generation method
      INTEGER,PARAMETER:: ih_height_gen      = 17

      ! First rho level at which height is constant
      INTEGER,PARAMETER:: ih_1_c_rho_level   = 24

      INTEGER,PARAMETER:: ih_land_points     = 25 ! No. of land points
      INTEGER,PARAMETER:: ih_ozone_levels    = 26 ! No. of ozone levels

      ! No. of deep soil moisture levels
      INTEGER,PARAMETER:: ih_soilQ_levels    = 28

      ! Number of convective cloud levels
      INTEGER,PARAMETER:: ih_convect_levels  = 34
      INTEGER,PARAMETER:: ih_rad_step        = 35 ! Radiation timestep
      INTEGER,PARAMETER:: ih_AMIP_flag       = 36 ! Flag for AMIP run
      INTEGER,PARAMETER:: ih_AMIP_year       = 37 ! First AMIP year
      INTEGER,PARAMETER:: ih_AMIP_month      = 38 ! First AMIP month
      INTEGER,PARAMETER:: ih_AMIP_day        = 49 ! First AMIP day
      INTEGER,PARAMETER:: ih_ozone_month     = 40 ! Current ozone month
      INTEGER,PARAMETER:: ih_SH_zonal_quad   = 41 ! L_SH_zonal_quadratics
      INTEGER,PARAMETER:: ih_SH_zonal_begin  = 42 ! SH_zonal_begin
      INTEGER,PARAMETER:: ih_SH_zonal_period = 43 ! SH_zonal_period
      INTEGER,PARAMETER:: ih_SH_level_weight = 44 ! SuHe_level_weight
      INTEGER,PARAMETER:: ih_SH_sigma_cutoff = 45 ! SuHe_sigma_cutoff
      INTEGER,PARAMETER:: ih_friction_time   = 46 ! frictional_timescale

! IHEADAPM end
      ! Meaningful parameter names for real constants header:
! ----------------------- include file: RHEADAPM -----------------------
! Description: Meaningful parameter names to index A_REALHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.

      ! East-West   grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaEW         = 1

      ! North-South grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaNS         = 2

      ! Latitude  of first p point in degrees
      INTEGER,PARAMETER:: rh_baselat         = 3

      ! Longitude of first p point in degrees
      INTEGER,PARAMETER:: rh_baselong        = 4

      ! Latitude  of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlat          = 5

      ! Longitude of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlong         = 6

      ! Height of top theta level (m)
      INTEGER,PARAMETER:: rh_z_top_theta     =16

      ! total moisture of the atmosphere
      INTEGER,PARAMETER:: rh_tot_m_init      =18

      ! total mass of atmosphere
      INTEGER,PARAMETER:: rh_tot_mass_init   =19

      ! total energy of atmosphere
      INTEGER,PARAMETER:: rh_tot_energy_init =20

      ! energy correction = energy drift
      INTEGER,PARAMETER:: rh_energy_corr     =21

! RHEADAPM end
      ! Meaningful parameter names for fixed header:
! ----------------------- include file: FHEADAPM -----------------------
! Description: Meaningful parameter names to index A_FIXHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.
 
      ! Start of Row Dependent Constant
      INTEGER,PARAMETER:: fh_RowDepCStart   = 115

      ! Start of Col Dependent Constant
      INTEGER,PARAMETER:: fh_ColDepCStart   = 120

! FHEADAPM end
      ! PP headers

      INTEGER :: A_LOOKUP(LEN1_LOOKUP,A_LEN2_LOOKUP) ! lookup heads
      INTEGER :: A_MPP_LOOKUP(MPP_LEN1_LOOKUP,A_LEN2_LOOKUP)
      INTEGER :: a_ixsts(len_a_ixsts)     ! stash index array

      REAL    :: a_spsts(len_a_spsts)     ! atmos stash array
! TYPDUMA end
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end
! TYPPTRA
! include TYPSIZE first
! Pointers for ATMOSPHERE model variables. Configuration dependent.
! Addresses in D1 array of model variables held in primary and
! secondary space.

      ! 1: Array Variables (dimensions are resolution dependent.)

      ! 1.1: Data variables stored in primary space.

      INTEGER :: JEXNERSURF                          ! surface exner
      INTEGER :: JDRYRHO(pdims%k_start:pdims%k_end)  ! Dry Density
      INTEGER :: JETADOT(wdims%k_start:wdims%k_end)  ! Etadot
      INTEGER :: JTHETAV(tdims%k_start:tdims%k_end)  ! Potential virtual
                                                     ! dry temperature
      INTEGER :: JPSIWS
      INTEGER :: JPSIWL  
      INTEGER :: JMV     (qdims%k_start:qdims%k_end) ! humidity mixing ratio
      INTEGER :: JMCL    (qdims%k_start:qdims%k_end) ! humidity mixing ratio
      INTEGER :: JMCF    (qdims%k_start:qdims%k_end) ! humidity mixing ratio
      INTEGER :: JMCF2   (qdims%k_start:qdims%k_end) ! humidity mixing ratio
      INTEGER :: JMRAIN  (qdims%k_start:qdims%k_end) ! humidity mixing ratio
      INTEGER :: JMGRAUP (qdims%k_start:qdims%k_end) ! humidity mixing ratio


      INTEGER :: JU(udims%k_start:udims%k_end)    ! u component of wind
      INTEGER :: JV(vdims%k_start:vdims%k_end)    ! v component of wind
      INTEGER :: JW(wdims_s%k_start:wdims_s%k_end)  ! w component of wind
      INTEGER :: JRHO(pdims%k_start:pdims%k_end)    ! Density
      INTEGER :: JTHETA(tdims%k_start:tdims%k_end)  ! Potential temperature
      INTEGER :: JQ     (qdims%k_start:qdims%k_end) ! Specific humidity
      INTEGER :: JQCL   (qdims%k_start:qdims%k_end) ! qcl
      INTEGER :: JQCF   (qdims%k_start:qdims%k_end) ! qcf
      INTEGER :: JQCF2  (qdims%k_start:qdims%k_end) ! second ice
      INTEGER :: JQRAIN (qdims%k_start:qdims%k_end) ! rain
      INTEGER :: JQGRAUP(qdims%k_start:qdims%k_end) ! graupel

      INTEGER :: JE_TRB  (tdims%k_start:tdims%k_end) ! TKE
      INTEGER :: JTSQ_TRB(tdims%k_start:tdims%k_end) ! Self covariance of thetal'
      INTEGER :: JQSQ_TRB(tdims%k_start:tdims%k_end) ! Self coveriance of qw'
      INTEGER :: JCOV_TRB(tdims%k_start:tdims%k_end) ! Correlation between thetal' and qw'
      INTEGER :: JZHPAR_SHCU            ! height of mixed layer used
                                        ! to evaluate the non-grad buoy flux

      ! Exner pressure on rho levels
      INTEGER :: JEXNER_RHO_LEVELS(pdims%k_start:pdims%k_end+1)

      INTEGER :: JU_ADV(udims%k_start:udims%k_end) ! Advective u component of wind
      INTEGER :: JV_ADV(vdims%k_start:vdims%k_end) ! Advective v component of wind
      INTEGER :: JW_ADV(wdims_s%k_start:wdims_s%k_end) ! Advective w component of wind

      ! 1.2: Data variables stored in secondary space.
      INTEGER :: JP(pdims%k_start:pdims%k_end+1)  ! Pressure on rho le

      ! Pressure on theta levels
      INTEGER :: JP_THETA_LEVELS(tdims%k_start:tdims%k_end)

      ! Exner pressure on theta levels
      INTEGER :: JEXNER_THETA_LEVELS(tdims%k_start:tdims%k_end)

      ! 1.3: Cloud Fields
      INTEGER :: JCCW_RAD(qdims%k_end)               ! CCW profile to radiation
      INTEGER :: JCCA    (N_CCA_LEV)                 ! Convective cloud amount
                                                     ! n_cca_lev set in dervsize 
      INTEGER :: JCF_AREA  (              qdims%k_end) ! Area Cloud Fraction
      INTEGER :: JCF_BULK  (qdims%k_start:qdims%k_end) ! Bulk Cloud Fraction
      INTEGER :: JCF_LIQUID(qdims%k_start:qdims%k_end) ! Liquid cloud fraction
      INTEGER :: JCF_FROZEN(qdims%k_start:qdims%k_end) ! Frozen cloud fraction

      ! 1.4: Soil Ancillary fields
      INTEGER :: J_DEEP_SOIL_TEMP(ST_LEVELS)   ! Deep soil temperature

      INTEGER :: JSMCL(SM_LEVELS)   ! Soil moisture content in layers
      INTEGER :: JSTHU(SM_LEVELS)   ! Unfrozen soil moisture fraction
      INTEGER :: JSTHF(SM_LEVELS)   ! Frozen soil moisture fraction

      ! 1.5: Radiation Increments 
      !     (best not to use EG module refs here)
      INTEGER :: JSW_INCS(0:model_levels+1)    ! SW radiation increments
      INTEGER :: JLW_INCS(0:model_levels)      ! LW radiation increments
! PAR radiation increment
      INTEGER :: JDIRPAR

      ! 1.6: Ozone
      INTEGER :: JOZONE(o3dims2%k_start:o3dims2%k_end)   ! Ozone
!  tropopause-based ozone
      INTEGER :: JTPPSOZONE(tpps_ozone_levels)

      ! 1.7: Tracer and aerosol fields
      INTEGER :: JTRACER (trdims_xstl%k_start:trdims_xstl%k_end,TR_VARS+1)  ! Tracers
      INTEGER :: JTR_UKCA(trdims_xstl%k_start:trdims_xstl%k_end,TR_UKCA+1)  ! UKCA Tracers

      INTEGER :: JMURK_SOURCE(tdims%k_start:tdims%k_end) ! multi-level murk source
      INTEGER :: JMURK       (tdims%k_start:tdims%k_end) ! multi-level murk content

      INTEGER :: JDUST_DIV1(tdims%k_start:tdims%k_end)   ! dust mmr, division 1
      INTEGER :: JDUST_DIV2(tdims%k_start:tdims%k_end)   ! dust mmr, division 2
      INTEGER :: JDUST_DIV3(tdims%k_start:tdims%k_end)   ! dust mmr, division 3
      INTEGER :: JDUST_DIV4(tdims%k_start:tdims%k_end)   ! dust mmr, division 4
      INTEGER :: JDUST_DIV5(tdims%k_start:tdims%k_end)   ! dust mmr, division 5
      INTEGER :: JDUST_DIV6(tdims%k_start:tdims%k_end)   ! dust mmr, division 6

      INTEGER :: JSO2       (tdims%k_start:tdims%k_end) ! sulphur dioxide gas
      INTEGER :: JDMS       (tdims%k_start:tdims%k_end) ! dimethyl sulphide gas
      INTEGER :: JSO4_AITKEN(tdims%k_start:tdims%k_end) ! Aitken mode sulphate aer
      INTEGER :: JSO4_ACCU  (tdims%k_start:tdims%k_end) ! accumulation mode sulpha
      INTEGER :: JSO4_DISS  (tdims%k_start:tdims%k_end) ! dissloved  sulphate aero
      INTEGER :: JH2O2      (tdims%k_start:tdims%k_end) ! hydrogen peroxide mmr
      INTEGER :: JNH3       (tdims%k_start:tdims%k_end) ! ammonia gas mmr

      INTEGER :: JSOOT_NEW(tdims%k_start:tdims%k_end)       ! fresh soot mmr
      INTEGER :: JSOOT_AGD(tdims%k_start:tdims%k_end)       ! aged soot mmr
      INTEGER :: JSOOT_CLD(tdims%k_start:tdims%k_end)       ! soot in cloud mmr

      INTEGER :: JBMASS_NEW(tdims%k_start:tdims%k_end)      ! fresh biomass mmr
      INTEGER :: JBMASS_AGD(tdims%k_start:tdims%k_end)      ! aged biomass mmr
      INTEGER :: JBMASS_CLD(tdims%k_start:tdims%k_end)      ! cloud biomass mmr

      INTEGER :: JOCFF_NEW(tdims%k_start:tdims%k_end)       ! fresh OCFF mmr
      INTEGER :: JOCFF_AGD(tdims%k_start:tdims%k_end)       ! aged OCFF mmr
      INTEGER :: JOCFF_CLD(tdims%k_start:tdims%k_end)       ! OCFF in cloud mmr

      INTEGER :: JSO2_NATEM (tdims%k_start:tdims%k_end) ! natural SO2 emissions
      INTEGER :: JOH        (tdims%k_start:tdims%k_end) ! hydroxyl radical ancilla
      INTEGER :: JHO2       (tdims%k_start:tdims%k_end) ! hydrogen dioxide ancilla
      INTEGER :: JH2O2_LIMIT(tdims%k_start:tdims%k_end) ! limiting H2O2 ancillary
      INTEGER :: JO3_CHEM   (tdims%k_start:tdims%k_end) ! ozone for chemistry anci
      INTEGER :: JHadCM2_SO4(2)                         ! HadCM2 sulphate loading
      ! JHadCM2_SO4: Should really be NSULPAT (= 2) but to use NSULPAT,
      !  must use clmchfcg_scenario_mod before every include of TYPPTR

      INTEGER :: JCO2      (tdims%k_start:tdims%k_end)      ! 3D CO2 FIELD

      INTEGER :: JOH_UKCA  (tdims%k_start:tdims%k_end)      ! OH MMR from UKCA
      INTEGER :: JHO2_UKCA (tdims%k_start:tdims%k_end)      ! HO2 MMR from UKCA
      INTEGER :: JH2O2_UKCA(tdims%k_start:tdims%k_end)      ! H2O2 MMR from UKCA
      INTEGER :: JO3_UKCA  (tdims%k_start:tdims%k_end)      ! O3 MMR from UKCA
      INTEGER :: JHNO3_UKCA(tdims%k_start:tdims%k_end)      ! HNO3 MMR from UKCA

      INTEGER :: JOZONE_TRACER(tdims%k_start:tdims%k_end)   ! Prognostic O3 Tracer(Cariol)
      INTEGER :: JO3_PROD_LOSS(tdims%k_start:tdims%k_end)   ! Cariol O3 Prod-Loss (P-L)
      INTEGER :: JO3_P_L_VMR  (tdims%k_start:tdims%k_end)   ! Cariol O3 P-L wrt VMR 
      INTEGER :: JO3_VMR      (tdims%k_start:tdims%k_end)   ! Cariol O3 Vol Mix Ratio-VMR
      INTEGER :: JO3_P_L_TEMP (tdims%k_start:tdims%k_end)   ! Cariol O3 P-L wrt temp 
      INTEGER :: JO3_TEMP     (tdims%k_start:tdims%k_end)   ! Cariol O3 temp  
      INTEGER :: JO3_P_L_COLO3(tdims%k_start:tdims%k_end)   ! Cariol O3 P-L wrt colO3 
      INTEGER :: JO3_COLO3    (tdims%k_start:tdims%k_end)   ! Cariol O3 column (colO3) 

      INTEGER :: JARCLBIOG_BG(tdims%k_start:tdims%k_end)    ! Biogenic aerosol climatology 
      INTEGER :: JARCLBIOM_FR(tdims%k_start:tdims%k_end)    ! Biomass burning (fresh) aerosol clim
      INTEGER :: JARCLBIOM_AG(tdims%k_start:tdims%k_end)    ! Biomass burning (aged) aerosol clim
      INTEGER :: JARCLBIOM_IC(tdims%k_start:tdims%k_end)    ! Biomass burning (in-cloud) aerosol clim

      INTEGER :: JARCLBLCK_FR(tdims%k_start:tdims%k_end)    ! Black carbon (fresh) aerosol clim
      INTEGER :: JARCLBLCK_AG(tdims%k_start:tdims%k_end)    ! Black carbon (aged) aerosol clim
      INTEGER :: JARCLSSLT_FI(tdims%k_start:tdims%k_end)    ! Sea salt (film mode) aerosol clim 
      INTEGER :: JARCLSSLT_JT(tdims%k_start:tdims%k_end)    ! Sea salt (jet mode) aerosol clim

      INTEGER :: JARCLSULP_AC(tdims%k_start:tdims%k_end)    ! Sulphate (accumulation mode) aero clim
      INTEGER :: JARCLSULP_AK(tdims%k_start:tdims%k_end)    ! Sulphate (Aitken mode) aerosol clim 
      INTEGER :: JARCLSULP_DI(tdims%k_start:tdims%k_end)    ! Sulphate (dissolved) aerosol clim

      INTEGER :: JARCLDUST_B1(tdims%k_start:tdims%k_end)    ! Dust (bin 1) aerosol climatology 
      INTEGER :: JARCLDUST_B2(tdims%k_start:tdims%k_end)    ! Dust (bin 2) aerosol climatology 
      INTEGER :: JARCLDUST_B3(tdims%k_start:tdims%k_end)    ! Dust (bin 3) aerosol climatology 
      INTEGER :: JARCLDUST_B4(tdims%k_start:tdims%k_end)    ! Dust (bin 4) aerosol climatology 
      INTEGER :: JARCLDUST_B5(tdims%k_start:tdims%k_end)    ! Dust (bin 5) aerosol climatology 
      INTEGER :: JARCLDUST_B6(tdims%k_start:tdims%k_end)    ! Dust (bin 6) aerosol climatology 

      INTEGER :: JARCLOCFF_FR(tdims%k_start:tdims%k_end)    ! Org carbon fossil fuel (fresh) aero clim
      INTEGER :: JARCLOCFF_AG(tdims%k_start:tdims%k_end)    ! Org carbon fossil fuel (aged) aero clim
      INTEGER :: JARCLOCFF_IC(tdims%k_start:tdims%k_end)    ! Org carbon fossil fuel (in-cloud) aero clim
      INTEGER :: JARCLDLTA_DL(tdims%k_start:tdims%k_end)    ! Delta aerosol climatology
      INTEGER :: JNITR_ACC(tdims%k_start:tdims%k_end)       ! Accumulation nitrate aerosol
      INTEGER :: JNITR_DISS(tdims%k_start:tdims%k_end)      ! Dissolved nitrate aerosol
!

! 1.8: Multi-level user ancillary fields
      INTEGER :: JUSER_MULT1(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT2(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT3(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT4(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT5(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT6(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT7(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT8(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT9(MODEL_LEVELS)     ! multi-level user ancilla
      INTEGER :: JUSER_MULT10(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT11(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT12(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT13(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT14(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT15(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT16(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT17(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT18(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT19(MODEL_LEVELS)    ! multi-level user ancilla
      INTEGER :: JUSER_MULT20(MODEL_LEVELS)    ! multi-level user ancilla

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      INTEGER :: JOROG_LBC                       ! Orography LBC
      INTEGER :: JU_LBC                          ! U LBC
      INTEGER :: JV_LBC                          ! V LBC
      INTEGER :: JW_LBC                          ! W LBC
      INTEGER :: JRHO_LBC                        ! RHO LBC
      INTEGER :: JTHETA_LBC                      ! Theta LBC
      INTEGER :: JQ_LBC                          ! Q LBC
      INTEGER :: JQCL_LBC                        ! QCL LBC
      INTEGER :: JQCF_LBC                        ! QCF LBC
      INTEGER :: JQCF2_LBC                       ! 2nd Ice LBC
      INTEGER :: JQRAIN_LBC                      ! Rain LBC
      INTEGER :: JQGRAUP_LBC                     ! Graupel LBC
      INTEGER :: JCF_BULK_LBC                    ! CF_BULK LBC
      INTEGER :: JCF_LIQUID_LBC                  ! CF_LIQUID_LBC
      INTEGER :: JCF_FROZEN_LBC                  ! CF_FROZEN_LBC
      INTEGER :: JEXNER_LBC                      ! Exner LBC
      INTEGER :: JU_ADV_LBC                      ! U_ADV LBC
      INTEGER :: JV_ADV_LBC                      ! V_ADV LBC
      INTEGER :: JW_ADV_LBC                      ! W_ADV LBC
      INTEGER :: JMURK_LBC                       ! Murk aerosol LBC
      INTEGER :: JTRACER_LBC(TR_LBC_VARS+1)      ! Tracer LBCs
      INTEGER :: JTR_UKCA_LBC(TR_LBC_UKCA+1)     ! UKCA Tracer LBCs
      INTEGER :: JDUST_DIV1_LBC                  ! DUST_DIV1 LBC
      INTEGER :: JDUST_DIV2_LBC                  ! DUST_DIV2 LBC
      INTEGER :: JDUST_DIV3_LBC                  ! DUST_DIV3 LBC
      INTEGER :: JDUST_DIV4_LBC                  ! DUST_DIV4 LBC
      INTEGER :: JDUST_DIV5_LBC                  ! DUST_DIV5 LBC 
      INTEGER :: JDUST_DIV6_LBC                  ! DUST_DIV6 LBC
      INTEGER :: JSO2_LBC                        ! SO2 LBC
      INTEGER :: JDMS_LBC                        ! DMS LBC
      INTEGER :: JSO4_AITKEN_LBC                 ! SO4_AITKEN LBC
      INTEGER :: JSO4_ACCU_LBC                   ! SO4_ACCU LBC
      INTEGER :: JSO4_DISS_LBC                   ! SO4_DISS_LBC
      INTEGER :: JNH3_LBC                        ! NH3 LBC
      INTEGER :: JSOOT_NEW_LBC                   ! SOOT_NEW LBC
      INTEGER :: JSOOT_AGD_LBC                   ! SOOT_AGD LBC
      INTEGER :: JSOOT_CLD_LBC                   ! SOOT_CLD LBC
      INTEGER :: JBMASS_NEW_LBC                  ! BMASS_NEW LBC
      INTEGER :: JBMASS_AGD_LBC                  ! BMASS_AGD LBC
      INTEGER :: JBMASS_CLD_LBC                  ! BMASS_CLD LBC
      INTEGER :: JOCFF_NEW_LBC                   ! OCFF_NEW LBC
      INTEGER :: JOCFF_AGD_LBC                   ! OCFF_AGD LBC
      INTEGER :: JOCFF_CLD_LBC                   ! OCFF_CLD LBC
      INTEGER :: JNITR_ACC_LBC                   ! NITR_ACC_LBC
      INTEGER :: JNITR_DISS_LBC                  ! NITR_DISS_LBC
      ! Lateral Boundary Condition tendencies
      INTEGER :: JU_LBC_TEND                     ! U LBC  tendencies
      INTEGER :: JV_LBC_TEND                     ! V LBC tendencies
      INTEGER :: JW_LBC_TEND                     ! W LBC tendencies
      INTEGER :: JRHO_LBC_TEND                   ! RHO LBC tendencies
      INTEGER :: JTHETA_LBC_TEND                 ! Theta LBC tendencies
      INTEGER :: JQ_LBC_TEND                     ! Q LBC tendencies
      INTEGER :: JQCL_LBC_TEND                   ! QCL LBC tendencies
      INTEGER :: JQCF_LBC_TEND                   ! QCF LBC tendencies
      INTEGER :: JQCF2_LBC_TEND                  ! 2nd Ice
      INTEGER :: JQRAIN_LBC_TEND                 ! Rain LBC tendencies
      INTEGER :: JQGRAUP_LBC_TEND                ! Graupel
      INTEGER :: JCF_BULK_LBC_TEND               ! CF_BULK LBC tend'cies
      INTEGER :: JCF_LIQUID_LBC_TEND             ! CF_LIQUID_LBC t'cies
      INTEGER :: JCF_FROZEN_LBC_TEND             ! CF_FROZEN_LBC t'cies
      INTEGER :: JEXNER_LBC_TEND                 ! Exner LBC tendencies
      INTEGER :: JU_ADV_LBC_TEND                 ! U_ADV LBC tendencies
      INTEGER :: JV_ADV_LBC_TEND                 ! V_ADV LBC tendencies
      INTEGER :: JW_ADV_LBC_TEND                 ! W_ADV LBCtendencies
      INTEGER :: JMURK_LBC_TEND                  ! Murk aerosol LBC tend
      INTEGER :: JTRACER_LBC_TEND(TR_LBC_VARS+1) ! Tracer LBC tendencies
      INTEGER :: JTR_UKCA_LBC_TEND(TR_LBC_UKCA+1)! UKCA Tracer LBC tend
      INTEGER :: JDUST_DIV1_LBC_TEND             ! DUST_DIV1 LBC tend
      INTEGER :: JDUST_DIV2_LBC_TEND             ! DUST_DIV2 LBC tend
      INTEGER :: JDUST_DIV3_LBC_TEND             ! DUST_DIV3 LBC tend
      INTEGER :: JDUST_DIV4_LBC_TEND             ! DUST_DIV4 LBC tend
      INTEGER :: JDUST_DIV5_LBC_TEND             ! DUST_DIV5 LBC tend
      INTEGER :: JDUST_DIV6_LBC_TEND             ! DUST_DIV6 LBC tend
      INTEGER :: JSO2_LBC_TEND                   ! SO2 LBC tend
      INTEGER :: JDMS_LBC_TEND                   ! DMS LBC tend
      INTEGER :: JSO4_AITKEN_LBC_TEND            ! SO4_AITKEN LBC tend
      INTEGER :: JSO4_ACCU_LBC_TEND              ! SO4_ACCU LBC tend
      INTEGER :: JSO4_DISS_LBC_TEND              ! SO4_DISS_LBC tend
      INTEGER :: JNH3_LBC_TEND                   ! NH3 LBC tend
      INTEGER :: JSOOT_NEW_LBC_TEND              ! SOOT_NEW LBC tend
      INTEGER :: JSOOT_AGD_LBC_TEND              ! SOOT_AGD LBC tend
      INTEGER :: JSOOT_CLD_LBC_TEND              ! SOOT_CLD LBC tend
      INTEGER :: JBMASS_NEW_LBC_TEND             ! BMASS_NEW LBC tend
      INTEGER :: JBMASS_AGD_LBC_TEND             ! BMASS_AGD LBC tend
      INTEGER :: JBMASS_CLD_LBC_TEND             ! BMASS_CLD LBC tend
      INTEGER :: JOCFF_NEW_LBC_TEND              ! OCFF_NEW LBC tend
      INTEGER :: JOCFF_AGD_LBC_TEND              ! OCFF_AGD LBC tend
      INTEGER :: JOCFF_CLD_LBC_TEND              ! OCFF_CLD LBC tend
      INTEGER :: JNITR_ACC_LBC_TEND              ! NITR_ACC_LBC tend
      INTEGER :: JNITR_DISS_LBC_TEND             ! NITR_DISS_LBC tend

      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      INTEGER :: JTSTAR         ! Surface temperature
      INTEGER :: JLAND          ! Land sea mask
      INTEGER :: JTSTAR_ANOM    ! Surface temperature anomaly
!   2.15: Fields for coastal tiling
      INTEGER :: JFRAC_LAND  ! Land fraction in grid box
      INTEGER :: JTSTAR_LAND ! Land surface temperature
      INTEGER :: JTSTAR_SEA  ! Sea surface temperature
      INTEGER :: JTSTAR_SICE ! Sea-ice surface temperature
      INTEGER :: JTSTAR_SICE_CAT ! Sea-ice surface temperature on categories
! Set pointers for sea-ice and land albedos
      INTEGER :: JSICE_ALB                   ! Sea-ice albedo
      INTEGER :: JLAND_ALB                   ! Mean land albedo

      ! 2.2: Data variables stored in secondary space.

      INTEGER :: JPSTAR          ! Surface pressure

      ! 2.3: Cloud fields
      INTEGER :: JLCBASE         ! Lowest Convective cloud base
      INTEGER :: JCCB            ! Convective cloud base
      INTEGER :: JCCT            ! Convective cloud top

      INTEGER :: JCCLWP          ! Convective cloud liquid water path

      INTEGER :: JDEEPFLAG       ! Flag for history of deep convection
      INTEGER :: JPASTPRECIP     ! Past convective precipitation
      INTEGER :: JPASTCONVHT     ! Past convective height

      ! 2.4: Boundary layer fields

      INTEGER :: JZH                         ! Boundary layer depth

      INTEGER :: jddmfx ! Convective downdraught 
!                       ! mass-flux at cloud-base 

      ! Standard deviation of turbulent fluctuations of layer 1
                                    ! temperature
      INTEGER :: JT1_SD

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      INTEGER :: JQ1_SD

      ! Decoupled screen-level temperature
      INTEGER :: JTScrnDcl_TILE
      INTEGER :: JTScrnDcl_SSI
      INTEGER :: JtStbTrans

      ! Number of model levels in the  turbulently mixed layer
      INTEGER :: JNTML

      ! Top level for turb mixing in any decoupled Sc layer
      INTEGER :: JNTDSC

      ! Bottom level for turb mixing in any decoupled Sc layer
      INTEGER :: JNBDSC

      INTEGER :: JCUMULUS      ! Boundary layer convection flag

      ! 2.4: Soil Ancillary fields

      INTEGER :: JSAT_SOILW_SUCTION          ! Saturated soil water sucti
      INTEGER :: JTHERM_CAP     ! Thermal capacity
      INTEGER :: JTHERM_COND    ! Thermal conductivity
      INTEGER :: JVOL_SMC_CRIT  ! Vol smc at critical point
      INTEGER :: JVOL_SMC_WILT  ! Vol smc at wilting point
      INTEGER :: JVOL_SMC_SAT   ! Vol smc at saturation
      INTEGER :: JSAT_SOIL_COND ! Saturated soil conductivity
      INTEGER :: JCLAPP_HORN    ! Clapp-Hornberger B coefficient

      ! 2.5: Other surface fields
      INTEGER :: JCANOPY_WATER  ! Canopy Water
      INTEGER :: JZ0            ! Roughness length; sea points on first timestep
      INTEGER :: JGS            ! Gridbox mean canopy conductance

      ! 2.6: Orographic Ancillary fields

      INTEGER :: JOROG          ! Orographic height
      INTEGER :: JOROG_SD       ! Standard Deviation of orography
      INTEGER :: JOROG_SIL      ! Silhouette area of orography
      INTEGER :: JOROG_HO2      ! Peak to trough height/(2*sqrt2)
      INTEGER :: JOROG_GRAD_X   ! Orographic gradient x
      INTEGER :: JOROG_GRAD_Y   ! Orographic gradient y
      INTEGER :: JOROG_UNFILT   ! Unfiltered orographic height
      INTEGER :: JOROG_GRAD_XX  ! Orographic gradient xx
      INTEGER :: JOROG_GRAD_XY  ! Orographic gradient xy
      INTEGER :: JOROG_GRAD_YY  ! Orographic gradient yy

      ! 2.7: Sea/Sea Ice

      INTEGER :: JU_SEA         ! Surface current (u component)
      INTEGER :: JV_SEA         ! Surface current (v component)
      INTEGER :: JU_0_P         ! Surace  current (u) on p-grid
      INTEGER :: JV_0_P         ! Surface current (v) on p-grid
      INTEGER :: JICE_FRACTION  ! Sea ice fraction
      INTEGER :: JICE_THICKNESS ! Sea ice depth
      INTEGER :: JTI            ! Sea ice temperature
      INTEGER :: JICE_FRACT_CAT ! Sea ice fraction on categories
      INTEGER :: JICE_THICK_CAT ! Sea ice thickness on categories
      INTEGER :: JTI_CAT        ! Sea ice temperature on categories
      INTEGER :: JICE_K_CAT     ! Sea ice effect conductivity on categories

      ! 2.8: Snow

      INTEGER :: JSNODEP        ! Snow depth on land
      INTEGER :: JSNODEP_SEA    ! Snow depth on sea ice
      INTEGER :: JSNODEP_SEA_CAT ! Snow depth on sea ice catagories
      INTEGER :: JCATCH_SNOW    ! Coniferous canopy
!                               ! snow capacity
      INTEGER :: JSNOW_GRND     ! Snow below canopy
      INTEGER :: JSNSOOT        ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

      INTEGER :: JSOIL_CLAY                    ! soil clay fraction
      INTEGER :: JSOIL_SILT                    ! soil silt fraction
      INTEGER :: JSOIL_SAND                    ! soil sand fraction
      INTEGER :: JDUST_MREL1                   ! soil rel mass, div 1
      INTEGER :: JDUST_MREL2                   ! soil rel mass, div 2
      INTEGER :: JDUST_MREL3                   ! soil rel mass, div 3
      INTEGER :: JDUST_MREL4                   ! soil rel mass, div 4
      INTEGER :: JDUST_MREL5                   ! soil rel mass, div 5
      INTEGER :: JDUST_MREL6                   ! soil rel mass, div 6


      INTEGER :: JSO2_EM        ! sulphur dioxide emission
      INTEGER :: JDMS_EM        ! dimethyl sulphide emission
      INTEGER :: JSO2_HILEM     ! high level SO2 emissions
      INTEGER :: JNH3_EM        ! ammonia gas surface emiss
      INTEGER :: JSOOT_EM       ! fresh soot surface emissions
      INTEGER :: JSOOT_HILEM    ! fresh soot high lev emissions
      INTEGER :: JBMASS_EM       ! fresh bmass surface emissions
      INTEGER :: JBMASS_HILEM    ! fresh bmass high lev emissions
      INTEGER :: JOCFF_EM        ! fresh OCFF surface emissions
      INTEGER :: JOCFF_HILEM     ! fresh OCFF high-level emissions
      INTEGER :: JDMS_CONC       ! seawater dimethyl sulphide conc.
      INTEGER :: JDMS_OFLUX      ! DMS flux from ocean model

      ! 2.10: User ancillary fields
      INTEGER :: JUSER_ANC1         ! user ancillary field 1
      INTEGER :: JUSER_ANC2                  ! user ancillary field 2
      INTEGER :: JUSER_ANC3                  ! user ancillary field 3
      INTEGER :: JUSER_ANC4                  ! user ancillary field 4
      INTEGER :: JUSER_ANC5                  ! user ancillary field 5
      INTEGER :: JUSER_ANC6                  ! user ancillary field 6
      INTEGER :: JUSER_ANC7                  ! user ancillary field 7
      INTEGER :: JUSER_ANC8                  ! user ancillary field 8
      INTEGER :: JUSER_ANC9                  ! user ancillary field 9
      INTEGER :: JUSER_ANC10                 !user ancillary field 10
      INTEGER :: JUSER_ANC11                 ! user ancillary field 11
      INTEGER :: JUSER_ANC12                 ! user ancillary field 12
      INTEGER :: JUSER_ANC13                 ! user ancillary field 13
      INTEGER :: JUSER_ANC14                 ! user ancillary field 14
      INTEGER :: JUSER_ANC15                 ! user ancillary field 15
      INTEGER :: JUSER_ANC16                 ! user ancillary field 16
      INTEGER :: JUSER_ANC17                 ! user ancillary field 17
      INTEGER :: JUSER_ANC18                 ! user ancillary field 18
      INTEGER :: JUSER_ANC19                 ! user ancillary field 19
      INTEGER :: JUSER_ANC20                 ! user ancillary field 20

      !   2.11: Store arrays for energy correction calculation
      INTEGER :: JNET_FLUX                   ! Net energy flux
      INTEGER :: JNET_MFLUX                  ! Net moisture flux

      !   2.12: Tiled Vegetation and Triffid fields
      INTEGER :: JFRAC_TYP        ! Fractions of surface type
      INTEGER :: JFRAC_CON1       ! Fractions of surface type
      INTEGER :: JFRAC_CON2       ! Fractions of surface type
      INTEGER :: JFRAC_CON3       ! Fractions of surface type
      INTEGER :: JFRAC_CON4       ! Fractions of surface type
      INTEGER :: JFRAC_CON5       ! Fractions of surface type
      INTEGER :: JFRAC_CON6       ! Fractions of surface type
      INTEGER :: JFRAC_CON7       ! Fractions of surface type
      INTEGER :: JFRAC_CON8       ! Fractions of surface type
      INTEGER :: JFRAC_CON9       ! Fractions of surface type
      INTEGER :: JLAI_PFT         ! LAI of plant functional types
      INTEGER :: JCANHT_PFT       ! Canopy hght of plant func types
      INTEGER :: JDISTURB         ! Disturbed fraction of vegetation
      INTEGER :: jsoil_alb        ! Snow-free albedo of bare soil
      INTEGER :: jobs_alb_sw      ! Observed snow-free SW albedo 
      INTEGER :: jobs_alb_vis     ! Observed snow-free VIS albedo 
      INTEGER :: jobs_alb_nir     ! Observed snow-free NIR albedo 
      INTEGER :: JSOIL_CARB       ! Soil carbon content
      INTEGER :: JSOIL_CARB1      ! Soil carbon content DPM
      INTEGER :: JSOIL_CARB2      ! Soil carbon content RPM
      INTEGER :: JSOIL_CARB3      ! Soil carbon content BIO
      INTEGER :: JSOIL_CARB4      ! Soil carbon content HUM
      INTEGER :: JNPP_PFT_ACC     ! Accumulated NPP on PFTs
      INTEGER :: JG_LF_PFT_ACC    ! Accum. leaf turnover rate PFTs
      INTEGER :: JG_PHLF_PFT_ACC  ! Accumulated phenological leaf
                                    ! turnover rate on PFTs
      INTEGER :: JRSP_W_PFT_ACC   ! Accum. wood respiration on PFTs
      INTEGER :: JRSP_S_ACC       ! Accumulated soil respiration
      INTEGER :: JRSP_S_ACC1      ! Accumulated soil respiration DPM
      INTEGER :: JRSP_S_ACC2      ! Accumulated soil respiration RPM
      INTEGER :: JRSP_S_ACC3      ! Accumulated soil respiration BIO
      INTEGER :: JRSP_S_ACC4      ! Accumulated soil respiration HUM
      INTEGER :: JCAN_WATER_TILE  ! Canopy water content on tiles
      INTEGER :: JCATCH_TILE      ! Canopy capacity on tiles
      INTEGER :: JINFIL_TILE      ! Max infiltration rate on tiles
      INTEGER :: JRGRAIN_TILE     ! Snow grain size on tiles
      INTEGER :: JSNODEP_TILE     ! Snow depth on tiles
      INTEGER :: JTSTAR_TILE      ! Surface temperature on tiles
      INTEGER :: JZ0_TILE         ! Surface roughness on tiles
      INTEGER :: JZ0H_TILE        ! Surface thermal roughness on tiles
      INTEGER :: JDOLR            ! TOA - surface upward LW at
                                  ! radiation timestep
      INTEGER :: JLW_DOWN         ! Surface downward LW at
                                  ! radiation timestep
      INTEGER :: JSW_TILE         ! Surface net SW on land tiles at
                                  ! radiation timestep
      INTEGER :: jurbhgt          ! Building height
      INTEGER :: jurbhwr          ! Urban H/W ratio
      INTEGER :: jurbwrr          ! Width ratio
      INTEGER :: jurbdisp         ! Displacement height
      INTEGER :: jurbztm          !
      INTEGER :: jurbalbwl        ! Wall albedo
      INTEGER :: jurbalbrd        ! Road albedo
      INTEGER :: jurbemisw        ! Wall emmissivity
      INTEGER :: jurbemisr        ! Road emmissivity

      !   2.13: Slab Model

!   2.14: Carbon cycle fields
      INTEGER :: J_CO2FLUX   ! Ocean CO2 flux (Kg CO2/m2/s1)
      INTEGER :: J_CO2_EMITS ! Surface CO2 emissions (Kg CO2/m2/s1)

!   2.15: Fields carried forward from previous version.
!         May not be required
      INTEGER :: JSURF_RESIST_NIT  ! Surface resistance on
                                    ! non-ice tiles
      INTEGER :: JROOT_DEPTH_NIT   ! Root depth on non-ice tiles
      INTEGER :: JZ0V_TYP          ! Vegetative Roughness length on
                                    ! tiles
      INTEGER :: JTSNOW            ! Snow surface layer temperature
      INTEGER :: JICE_EDGE
      INTEGER :: JOROG_TENDENCY    ! Orographic tendencies
      INTEGER :: JOROG_SD_TENDENCY ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
      INTEGER :: JETATHETA
      INTEGER :: JETARHO
      INTEGER :: JRHCRIT
      INTEGER :: JSOIL_THICKNESS
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      INTEGER :: Jzseak_theta ! zsea(k) on theta levels
      INTEGER :: JCk_theta    ! C(k)    on theta levels
      INTEGER :: Jzseak_rho   ! zsea(k) on rho levels
      INTEGER :: JCk_rho      ! C(k)    on rho levels
      ! Addresses in Row and Col  dependent constants array.
      INTEGER :: JLAMBDA_INPUT_P
      INTEGER :: JLAMBDA_INPUT_U
      INTEGER :: JPHI_INPUT_P
      INTEGER :: JPHI_INPUT_V
!   2.16: Fields for large-scale hydrology scheme.

      INTEGER :: JTI_MEAN          !Mean topographic index
      INTEGER :: JTI_SIG           !Standard dev. in topographic index
      INTEGER :: JFEXP             !Exponential decay in soil
!                                  ! saturated conductivity
      INTEGER :: JGAMMA_INT        !Integrated gamma function
      INTEGER :: JWATER_TABLE      !Water table depth
      INTEGER :: JFSFC_SAT         !Surface saturation fraction
      INTEGER :: JF_WETLAND        !Wetland fraction
      INTEGER :: JSTHZW           !Soil moist fract. in deep-zw layer.
      INTEGER :: JA_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JC_FSAT          !Fitting parameter for Fsat in LSH.
      INTEGER :: JA_FWET          !Fitting parameter for Fwet in LSH.
      INTEGER :: JC_FWET          !Fitting parameter for Fwet in LSH.


!   2.17: Fields for River routing.
      INTEGER :: JRIV_SEQUENCE   ! River sequence
      INTEGER :: JRIV_DIRECTION  ! River direction
      INTEGER :: JRIV_STORAGE    ! River water storage
      INTEGER :: JTOT_SURFROFF   ! Accumulated surface runoff
      INTEGER :: JTOT_SUBROFF    !     "       sub-surface runoff
      INTEGER :: JRIV_INLANDATM       ! inland basin outflow
! Field for water conservation due to lake evaporation 
      INTEGER :: JACC_LAKE_EVAP  ! Accumulated lake evaporation 
! Fields for grid-to-grid river routing (river routing 2A)
      INTEGER :: JRIV_IAREA      ! Drainage area
      INTEGER :: JRIV_SLOPE      ! Grid-cell slope
      INTEGER :: JRIV_FLOWOBS1   ! Initial values of flow
      INTEGER :: JRIV_INEXT      ! Flow direction (x)
      INTEGER :: JRIV_JNEXT      ! Flow direction (y)
      INTEGER :: JRIV_LAND       ! Land-type (land/river/sea)
      INTEGER :: JRIV_SUBSTORE   ! Subsurface storage
      INTEGER :: JRIV_SURFSTORE  ! Surface storage
      INTEGER :: JRIV_FLOWIN     ! Surface inflow
      INTEGER :: JRIV_BFLOWIN    ! Subsurface inflow

      INTEGER :: JC_SOLAR 
      INTEGER :: JC_BLUE 
      INTEGER :: JC_LONGWAVE 
      INTEGER :: JC_TAUX 
      INTEGER :: JC_TAUY 
      INTEGER :: JC_W10 
      INTEGER :: JC_SENSIBLE 
      INTEGER :: JC_SUBLIM 
      INTEGER :: JC_EVAP 
      INTEGER :: JC_FCONDTOPN 
      INTEGER :: JC_TOPMELTN 
      INTEGER :: JC_LSRAIN 
      INTEGER :: JC_LSSNOW 
      INTEGER :: JC_CVRAIN 
      INTEGER :: JC_CVSNOW 
      INTEGER :: JC_RIVEROUT
      INTEGER :: JC_CALVING

!   2.18: JULES variables
      INTEGER :: JSNOWDEPTH      ! Snow depth on ground on tiles (m)
      INTEGER :: JRHO_SNOW_GRND  ! Snowpack bulk density (kg/m3)
      INTEGER :: JNSNOW          ! Number of snow layers on ground on tiles
      INTEGER :: JDS             ! Snow layer thickness (m)
      INTEGER :: JSICE           ! Snow layer ice mass on tiles (Kg/m2)
      INTEGER :: JSLIQ           ! Snow layer liquid mass on tiles (Kg/m2)
      INTEGER :: JTSNOWLAYER     ! Snow layer temperature (K)
      INTEGER :: JRHO_SNOW       ! Snow layer densities (kg/m3)
      INTEGER :: JRGRAINL        ! Snow layer grain size on tiles (microns)
!  FLake lake scheme
      INTEGER :: JLAKE_DEPTH
      INTEGER :: JLAKE_FETCH
      INTEGER :: JLAKE_T_MEAN
      INTEGER :: JLAKE_T_MXL
      INTEGER :: JLAKE_T_ICE
      INTEGER :: JLAKE_H_MXL
      INTEGER :: JLAKE_H_ICE
      INTEGER :: JLAKE_SHAPE
      INTEGER :: JLAKE_G_DT
!{CABLE: introduced tiled prognostics 
!    Fields for CABLE
      INTEGER :: JTSOIL_TILE(SM_LEVELS)  ! Tiled soil temperature
      INTEGER :: JSMCL_TILE(SM_LEVELS)   ! Tiled soil moisture content in layers
      INTEGER :: JSTHF_TILE(SM_LEVELS)   ! Tiled frozen soil moisture fraction
      INTEGER :: JSNOW_DEPTH3L(3)        ! Tiled snow depth
      INTEGER :: JSNOW_MASS3L(3)         ! Tiled snow mass
      INTEGER :: JSNOW_TMP3L(3)          ! Tiled snow temperature
      INTEGER :: JSNOW_RHO3L(3)          ! Tiled snow density
      INTEGER :: JSNOW_RHO1L             ! Tiled mean snow density
      INTEGER :: JSNOW_AGE               ! Tiled snow age
      INTEGER :: JSNOW_FLG3L             ! Flag for use of 3 level snow scheme
!}CABLE

! Addresses in D1 array of primary variables: scalars
      COMMON/CARGPT_ATMOS/                                              &
! Data variables
        JTSTAR, JLAND, JPSTAR, JTSTAR_ANOM,                             &
        JFRAC_LAND, JTSTAR_LAND, JTSTAR_SEA, JTSTAR_SICE,               &
        JTSTAR_SICE_CAT, JSICE_ALB, JLAND_ALB,                          &
      ! Cloud fields
        JCCB, JCCT, JCCLWP,JDEEPFLAG,JPASTPRECIP,JPASTCONVHT,           &
      ! Boundary layer fields
        JZH, jddmfx, JT1_SD, JQ1_SD, JTScrnDcl_TILE, JTScrnDcl_SSI,     &
        JtStbTrans, JNTML, JNTDSC, JNBDSC, JCUMULUS,                    &
      ! Soil Ancillary fields
        JSAT_SOILW_SUCTION, JTHERM_CAP, JTHERM_COND,                    &
        JVOL_SMC_CRIT, JVOL_SMC_WILT, JVOL_SMC_SAT,                     &
        JSAT_SOIL_COND, JCLAPP_HORN,                                    &
      ! Other surface fields
        JCANOPY_WATER, JZ0, JGS,                                        &
      ! Orographic Ancillary fields
        JOROG, JOROG_SD, JOROG_SIL, JOROG_HO2,                          &
        JOROG_GRAD_X, JOROG_GRAD_Y, JOROG_UNFILT,                       &
        JOROG_GRAD_XX, JOROG_GRAD_XY, JOROG_GRAD_YY,                    &
      ! Sea/Sea Ice
        JU_SEA, JV_SEA, JU_0_P, JV_0_P,                                 &
        JICE_FRACTION, JICE_THICKNESS, JTI,                             &
        JICE_FRACT_CAT, JICE_THICK_CAT, JTI_CAT, JICE_K_CAT,            &
      ! Snow
        JSNODEP, JSNODEP_SEA, JSNSOOT, JCATCH_SNOW, JSNOW_GRND,         &
        JSNODEP_SEA_CAT,                                                &
      ! Aerosol emission fields,including mineral dust parent soil props
        JSOIL_CLAY, JSOIL_SILT, JSOIL_SAND,JDUST_MREL1,JDUST_MREL2,     &
        JDUST_MREL3,JDUST_MREL4,JDUST_MREL5,JDUST_MREL6,                &
        JSO2_EM, JDMS_EM, JSO2_HILEM, JNH3_EM, JSOOT_EM, JSOOT_HILEM,   &
        JBMASS_EM, JBMASS_HILEM, JOCFF_EM, JOCFF_HILEM, JDMS_CONC,      &
        JDMS_OFLUX,                                                     &
      ! User ancillary fields
        JUSER_ANC1,  JUSER_ANC2, JUSER_ANC3, JUSER_ANC4,                &
        JUSER_ANC5,  JUSER_ANC6, JUSER_ANC7, JUSER_ANC8,                &
        JUSER_ANC9,  JUSER_ANC10, JUSER_ANC11, JUSER_ANC12,             &
        JUSER_ANC13,  JUSER_ANC14, JUSER_ANC15, JUSER_ANC16,            &
        JUSER_ANC17,  JUSER_ANC18, JUSER_ANC19, JUSER_ANC20,            &
      ! Pointers for ATMOSPHERE model constants. Scalars only.
        JETATHETA, JETARHO, JRHCRIT, JSOIL_THICKNESS,                   &
        Jzseak_theta,Jck_theta,Jzseak_rho,Jck_rho,                      &
      ! pointers for input variable grid info.
        JLAMBDA_INPUT_P,  JLAMBDA_INPUT_U,                              &
        JPHI_INPUT_P, JPHI_INPUT_V,                                     &
      ! pointers for tiled vegetation and triffid
        JFRAC_TYP, JLAI_PFT, JCANHT_PFT, JDISTURB, jsoil_alb,           &
        jobs_alb_sw, jobs_alb_vis, jobs_alb_nir,                        & 
        JFRAC_CON1, JFRAC_CON2, JFRAC_CON3, JFRAC_CON4, JFRAC_CON5,     &
        JFRAC_CON6, JFRAC_CON7, JFRAC_CON8, JFRAC_CON9,                 &
        JSOIL_CARB,                                                     &
        JSOIL_CARB1, JSOIL_CARB2, JSOIL_CARB3, JSOIL_CARB4,             &
        JNPP_PFT_ACC, JG_LF_PFT_ACC, JG_PHLF_PFT_ACC,                   &
        JRSP_W_PFT_ACC, JRSP_S_ACC,                                     &
        JRSP_S_ACC1, JRSP_S_ACC2, JRSP_S_ACC3,                          &
        JRSP_S_ACC4, JCAN_WATER_TILE, JCATCH_TILE,                      &
        JINFIL_TILE, JRGRAIN_TILE, JSNODEP_TILE,                        &
        JTSTAR_TILE, JZ0_TILE, JZ0H_TILE,                               &
        JDOLR, JLW_DOWN, JSW_TILE,JNET_FLUX,JNET_MFLUX,                 &
      ! Pointers for MORUSES - new two-tile urban scheme
        jurbhgt, jurbhwr, jurbwrr, jurbdisp, jurbztm,                   &
        jurbalbwl, jurbalbrd, jurbemisw, jurbemisr,                     &
      ! pointers for carbon cycle
        J_CO2FLUX,J_CO2_EMITS,                                          &
      ! pointers for large-scale hydrology
        JFEXP, JTI_MEAN, JTI_SIG, JGAMMA_INT,                           &
        JWATER_TABLE, JFSFC_SAT, JF_WETLAND, JSTHZW,                    &
        JA_FSAT,      JC_FSAT,   JA_FWET,    JC_FWET,                   &
      ! pointers for river routing
        JRIV_SEQUENCE, JRIV_DIRECTION, JRIV_STORAGE,                    &
        JTOT_SURFROFF, JTOT_SUBROFF, JRIV_INLANDATM,                    & 
        JACC_LAKE_EVAP                                                  &
      ! pointers for coupling fields
       , JC_SOLAR, JC_BLUE, JC_LONGWAVE, JC_TAUX                        &
       , JC_TAUY, JC_W10, JC_SENSIBLE, JC_SUBLIM                        &
       , JC_EVAP, JC_FCONDTOPN, JC_TOPMELTN, JC_LSRAIN                  & 
       , JC_LSSNOW, JC_CVRAIN, JC_CVSNOW, JC_RIVEROUT, JC_CALVING       &
      ! pointers for JULES
       , JSNOWDEPTH, JRHO_SNOW_GRND                                     &
       , JNSNOW                                                         &
       , JDS, JSICE, JSLIQ, JTSNOWLAYER, JRHO_SNOW, JRGRAINL            &
      ! FLake lake scheme
       , JLAKE_DEPTH, JLAKE_FETCH, JLAKE_T_MEAN, JLAKE_T_MXL            &
       , JLAKE_T_ICE, JLAKE_H_MXL, JLAKE_H_ICE,  JLAKE_SHAPE            &
       , JLAKE_G_DT

! TYPPTRA end
      INTEGER                                                           &
          ICODE                  ! OUT: Error return code
!
      CHARACTER(LEN=80)                                                      &
          CMESSAGE               ! OUT: Error return message
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1
      !First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150
      !Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
      !Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is -1.
      ! Similarly for A_UKCA_INDEX.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      ! A_TR_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: A_TR_StashItem(A_MAX_TRVARS)

      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)
      ! UKCA_tr_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: UKCA_tr_StashItem(A_MAX_UKCAVARS) 

      ! A_TR_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: A_TR_LBC_StashItem(A_MAX_TRVARS) 
      INTEGER :: A_TR_active_lbc_index(A_MAX_TRVARS) 

      ! UKCA_tr_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: UKCA_tr_LBC_StashItem(A_MAX_UKCAVARS) 
      INTEGER :: UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)

      COMMON/ATRACER/A_TR_INDEX, A_TR_StashItem,                        &
     &               A_TR_LBC_StashItem, A_TR_active_lbc_index,         &
     &               A_UKCA_INDEX, UKCA_tr_StashItem,                   &
     &               UKCA_tr_LBC_StashItem, UKCA_tr_active_lbc_index

! CTRACERA end
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

! local variables

      INTEGER                                                           &
              IVAR,                                                     &
                                 ! Loop counts
              JVAR,                                                     &
                                 ! Loop counts
              IFLD,                                                     &
              LEV                                                       &
             ,im_ident                                                  &
                            !  Internal Model Identifier
             ,im_index                                                  &
                            !  Internal Model Index in Stash arrays
             ,Sect_No       !  Stash section number

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     Set to atmosphere internal model
      IF (lhook) CALL dr_hook('SET_ATM_POINTERS',zhook_in,zhook_handle)
      im_ident  = atmos_im
      im_index  = internal_model_index(im_ident)
      Sect_No   = 0


! Set pointers for atmospheric primary variables from STASH :-

      JU(udims%k_start)                = SI(  2,Sect_No,im_index)
      JV(vdims%k_start)                = SI(  3,Sect_No,im_index)
      JTHETA(tdims%k_start)            = SI(  4,Sect_No,im_index)
      JQ    (qdims%k_start)            = SI( 10,Sect_No,im_index)
      JQCF  (qdims%k_start)            = SI( 12,Sect_No,im_index)
      JTSTAR                           = SI( 24,Sect_No,im_index)
      JLAND                            = SI( 30,Sect_No,im_index)
      JOROG                            = SI( 33,Sect_No,im_index)
      JW(wdims_s%k_start)              = SI(150,Sect_No,im_index)
      JRHO   (pdims%k_start)           = SI(253,Sect_No,im_index)
      JQCL   (qdims%k_start)           = SI(254,Sect_No,im_index)
      JEXNER_RHO_LEVELS(pdims%k_start) = SI(255,Sect_No,im_index)
      JQCF2  (qdims%k_start)           = SI(271,Sect_No,im_index)
      JQRAIN (qdims%k_start)           = SI(272,Sect_No,im_index)
      JQGRAUP(qdims%k_start)           = SI(273,Sect_No,im_index)

!     ENDGame prognostics
      JTHETAV(tdims%k_start)           = SI(388,Sect_No,im_index)
      JDRYRHO(pdims%k_start)           = SI(389,Sect_No,im_index)
      JETADOT(tdims%k_start)           = SI(387,Sect_No,im_index)
      JPSIWS                           = SI(390,Sect_No,im_index)
      JPSIWL                           = SI(397,Sect_No,im_index)
      JEXNERSURF                       = SI(398,Sect_No,im_index)
      JMV(qdims%k_start)               = SI(391,Sect_No,im_index)
      JMCL(qdims%k_start)              = SI(392,Sect_No,im_index)
      JMCF(qdims%k_start)              = SI(393,Sect_No,im_index)
      JMCF2(qdims%k_start)             = SI(396,Sect_No,im_index)
      JMRAIN(qdims%k_start)            = SI(394,Sect_No,im_index)
      JMGRAUP(qdims%k_start)           = SI(395,Sect_No,im_index)

!     for TKE based turbulence scheme
      JE_TRB  (tdims%k_start)          = SI(70,Sect_No,im_index)
      JTSQ_TRB(tdims%k_start)          = SI(71,Sect_No,im_index)
      JQSQ_TRB(tdims%k_start)          = SI(72,Sect_No,im_index)
      JCOV_TRB(tdims%k_start)          = SI(73,Sect_No,im_index)
      JZHPAR_SHCU                      = SI(74,Sect_No,im_index)



! Set extra pointers for coastal tiling.
      JFRAC_LAND     = SI(505,Sect_No,im_index)
      JTSTAR_LAND    = SI(506,Sect_No,im_index)
      JTSTAR_SEA     = SI(507,Sect_No,im_index)
      JTSTAR_SICE    = SI(508,Sect_No,im_index)
      JTSTAR_SICE_CAT= SI(441,Sect_No,im_index)
! Set pointers for seaice and land albedos
      JSICE_ALB      = SI(509,Sect_No,im_index)
      JLAND_ALB      = SI(510,Sect_No,im_index)
!
! Set extra pointers for large-scale hydrology.
      JTI_MEAN       = SI(274,Sect_No,im_index)
      JTI_SIG        = SI(275,Sect_No,im_index)
      JFEXP          = SI(276,Sect_No,im_index)
      JGAMMA_INT     = SI(277,Sect_No,im_index)
      JWATER_TABLE   = SI(278,Sect_No,im_index)
      JFSFC_SAT      = SI(279,Sect_No,im_index)
      JF_WETLAND     = SI(280,Sect_No,im_index)
      JSTHZW         = SI(281,Sect_No,im_index)
      JA_FSAT        = SI(282,Sect_No,im_index)
      JC_FSAT        = SI(283,Sect_No,im_index)
      JA_FWET        = SI(284,Sect_No,im_index)
      JC_FWET        = SI(285,Sect_No,im_index)
      JACC_LAKE_EVAP = SI(290,Sect_No,im_index) 


      DO LEV= 1, wdims_s%k_end
        JW(LEV)=JW(LEV-1)+theta_off_size
      END DO

      DO LEV= udims%k_start+1, udims%k_end
        JU(LEV)    =JU(LEV-1)     + u_off_size
      END DO

      DO LEV= vdims%k_start+1, vdims%k_end
        JV(LEV)    =JV(LEV-1)     + v_off_size
      END DO

      DO LEV= pdims%k_start+1, pdims%k_end
        JRHO(LEV)    =JRHO(LEV-1)      + theta_off_size
        JDRYRHO(LEV) =JDRYRHO(LEV-1)   + theta_off_size
      END DO

      DO LEV= tdims%k_start+1, tdims%k_end
        JTHETA(LEV) =JTHETA(LEV-1)  + theta_off_size
        JTHETAV(LEV)=JTHETAV(LEV-1) + theta_off_size
        JETADOT(LEV)=JETADOT(LEV-1) + theta_off_size
      END DO

      DO LEV= pdims%k_start+1, pdims%k_end+1
        JEXNER_RHO_LEVELS(LEV)=JEXNER_RHO_LEVELS(LEV-1)+theta_off_size
      END DO

      DO LEV= qdims%k_start+1, qdims%k_end
        JQ(LEV)    =JQ(LEV-1)     + theta_halo_size
        JQCL(LEV)  =JQCL(LEV-1)   + theta_halo_size
        JQCF(LEV)  =JQCF(LEV-1)   + theta_halo_size
        JQCF2(LEV)   =JQCF2(LEV-1)   + theta_halo_size
        JQRAIN(LEV)  =JQRAIN(LEV-1)  + theta_halo_size
        JQGRAUP(LEV) =JQGRAUP(LEV-1) + theta_halo_size
        JMV(LEV)     =JMV(LEV-1)     + theta_off_size
        JMCL(LEV)    =JMCL(LEV-1)    + theta_off_size
        JMCF(LEV)    =JMCF(LEV-1)    + theta_off_size
        JMCF2(LEV)   =JMCF2(LEV-1)   + theta_off_size
        JMRAIN(LEV)  =JMRAIN(LEV-1)  + theta_off_size
        JMGRAUP(LEV) =JMGRAUP(LEV-1) + theta_off_size
      END DO

      DO LEV= tdims%k_start+1, tdims%k_end
        JE_TRB(LEV) = JE_TRB(LEV-1) + theta_halo_size
        JTSQ_TRB(LEV) = JTSQ_TRB(LEV-1) + theta_halo_size
        JQSQ_TRB(LEV) = JQSQ_TRB(LEV-1) + theta_halo_size
        JCOV_TRB(LEV) = JCOV_TRB(LEV-1) + theta_halo_size
      ENDDO
! 
!     Convective downdraught mass-flux at cloud base
      jddmfx         = SI(493,Sect_No,im_index) 

      ! Pointers required to save coupling fields as
      ! prognostics when employing OASIS as the coupler.
      ! In non-coupled models these pointers will
      ! simply end up with a value of 1.
      JC_SOLAR = SI(171,Sect_No,im_index)
      JC_BLUE =  SI(172,Sect_No,im_index)
      JC_LONGWAVE =  SI(174,Sect_No,im_index)
      JC_TAUX =  SI(176,Sect_No,im_index)
      JC_TAUY =  SI(177,Sect_No,im_index)
      JC_SENSIBLE =  SI(179,Sect_No,im_index)
      IF (NICE_USE .EQ. 1) THEN
        JC_SUBLIM =  SI(180,Sect_No,im_index)  ! Single cat field
      ELSE
        JC_SUBLIM =  SI(182,Sect_No,im_index)  ! Multi cat field
      ENDIF
      JC_EVAP =  SI(181,Sect_No,im_index)
      JC_FCONDTOPN =  SI(184,Sect_No,im_index)
      JC_TOPMELTN =  SI(185,Sect_No,im_index)
      JC_LSRAIN =  SI(186,Sect_No,im_index)
      JC_LSSNOW =  SI(187,Sect_No,im_index)
      JC_CVRAIN =  SI(188,Sect_No,im_index)
      JC_CVSNOW =  SI(189,Sect_No,im_index)
      JC_CALVING =  SI(190,Sect_No,im_index)
      JC_W10 =  SI(191,Sect_No,im_index)
      JC_RIVEROUT =  SI(192,Sect_No,im_index)


! Set pointers for optional atmospheric primary variables.
        JZH                         = SI( 25,Sect_No,im_index)
        JU_ADV(udims%k_start)       = SI(256,Sect_No,im_index)
        JV_ADV(vdims%k_start)       = SI(257,Sect_No,im_index)
        JW_ADV(wdims_s%k_start)     = SI(258,Sect_No,im_index)
        JNTML                       = SI(259,Sect_No,im_index)
        JNBDSC                      = SI(260,Sect_No,im_index)
        JNTDSC                      = SI(261,Sect_No,im_index)
        JCUMULUS                    = SI(262,Sect_No,im_index)
        JT1_SD                      = SI(263,Sect_No,im_index)
        JQ1_SD                      = SI(264,Sect_No,im_index)
        JCF_AREA  (1)               = SI(265,Sect_No,im_index)
        JCF_BULK  (qdims%k_start)   = SI(266,Sect_No,im_index)
        JCF_LIQUID(qdims%k_start)   = SI(267,Sect_No,im_index)
        JCF_FROZEN(qdims%k_start)   = SI(268,Sect_No,im_index)

        IF (L_3D_CCA .OR. L_CCRAD) THEN
          JCCA(1) = SI(211,Sect_No,im_index)
          DO LEV=2, N_CCA_LEV      ! n_cca_lev set in dervsize 
            JCCA(LEV)=JCCA(LEV-1)+THETA_FIELD_SIZE
          ENDDO
        ELSE
          JCCA(1)      = SI( 13,Sect_No,im_index)
        ENDIF

        JCCB           = SI( 14,Sect_No,im_index)
        JCCT           = SI( 15,Sect_No,im_index)
        JCCLWP         = SI( 16,Sect_No,im_index)
        jdeepflag      = SI(342,Sect_No,im_index)
        jpastprecip    = SI(343,Sect_No,im_index)
        jpastconvht    = SI(344,Sect_No,im_index)
        JLCBASE        = SI( 21,Sect_No,im_index)
        JCANOPY_WATER  = SI( 22,Sect_No,im_index)
        JCCW_RAD(1)    = SI(212,Sect_No,im_index)

        DO LEV= 2, qdims%k_end
          JCCW_RAD(LEV) =JCCW_RAD(LEV-1)+THETA_FIELD_SIZE
          JCF_AREA(LEV) =JCF_AREA(LEV-1)+THETA_FIELD_SIZE
        ENDDO

      DO LEV= wdims_s%k_start+1, wdims_s%k_end
        JW_ADV(LEV)=JW_ADV(LEV-1)+theta_halo_size
      END DO

      DO LEV= udims%k_start+1, udims%k_end
        JU_ADV(LEV)=JU_ADV(LEV-1)+u_halo_size
      END DO

      DO LEV= vdims%k_start+1, vdims%k_end
        JV_ADV(LEV)=JV_ADV(LEV-1)+v_halo_size
      END DO

      DO LEV= qdims%k_start+1, qdims%k_end
        JCF_BULK(LEV)  =JCF_BULK(LEV-1)+THETA_halo_SIZE
        JCF_LIQUID(LEV)=JCF_LIQUID(LEV-1)+THETA_halo_SIZE
        JCF_FROZEN(LEV)=JCF_FROZEN(LEV-1)+THETA_halo_SIZE
      END DO

!  Set pointers for secondary fields in D1
      JEXNER_THETA_LEVELS(tdims%k_start) = SI(406,Sect_No,im_index)
      JP                 (pdims%k_start) = SI(407,Sect_No,im_index)
      JP_THETA_LEVELS    (tdims%k_start) = SI(408,Sect_No,im_index)
      JPSTAR                             = SI(409,Sect_No,im_index)
      JSW_INCS(0)                        = SI(410,Sect_No,im_index)
      JLW_INCS(0)                        = SI(411,Sect_No,im_index)

! This code should work for both ND and V-AT-POLES
      DO LEV= 1, model_levels+1 
        JSW_INCS(LEV)=JSW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO
      DO LEV= 1, model_levels
        JLW_INCS(LEV)=JLW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO

! Direct PAR flux 
      JDIRPAR               = SI(460,Sect_no,im_index)

      DO LEV= tdims%k_start+1, tdims%k_end
        JEXNER_THETA_LEVELS(LEV)=JEXNER_THETA_LEVELS(LEV-1)+            &
                                   theta_off_size
        JP_THETA_LEVELS(LEV)    =JP_THETA_LEVELS(LEV-1)+theta_off_size
      END DO

      DO LEV= pdims%k_start+1, pdims%k_end+1
        JP(LEV)                 =JP(LEV-1)+theta_off_size
      END DO

!  Set pointers for ancillary fields in D1 from STASH
!     Soil fields
      JSMCL(1)            = SI(  9,Sect_No,im_index)
      J_DEEP_SOIL_TEMP(1) = SI( 20,Sect_No,im_index)
      JVOL_SMC_WILT       = SI( 40,Sect_No,im_index)
      JVOL_SMC_CRIT       = SI( 41,Sect_No,im_index)
      JVOL_SMC_SAT        = SI( 43,Sect_No,im_index)
      JSAT_SOIL_COND      = SI( 44,Sect_No,im_index)
      JTHERM_CAP          = SI( 46,Sect_No,im_index)
      JTHERM_COND         = SI( 47,Sect_No,im_index)
      JSAT_SOILW_SUCTION  = SI( 48,Sect_No,im_index)
      JCLAPP_HORN         = SI(207,Sect_No,im_index)
      JSTHU(1)            = SI(214,Sect_No,im_index)
      JSTHF(1)            = SI(215,Sect_No,im_index)

      DO LEV=2,ST_LEVELS
        J_DEEP_SOIL_TEMP(LEV)=J_DEEP_SOIL_TEMP(LEV-1)+LAND_FIELD
      ENDDO

      DO LEV=2,SM_LEVELS
        JSMCL(LEV)=JSMCL(LEV-1)+LAND_FIELD
        JSTHU(LEV)=JSTHU(LEV-1)+LAND_FIELD
        JSTHF(LEV)=JSTHF(LEV-1)+LAND_FIELD
      END DO

!     Other surface fields
      JZ0          = SI( 26,Sect_No,im_index) ! roughness length (used for sea)
      JGS          = SI(213,Sect_No,im_index) ! stomatal conductance

! Orography fields
      JOROG_SIL      = SI(17,Sect_No,im_index)   ! Silhouette area
      JOROG_HO2      = SI(18,Sect_No,im_index)   ! Peak to trough ht.
      JOROG_SD       = SI(34,Sect_No,im_index)
      JOROG_GRAD_X   = SI( 5,Sect_No,im_index)
      JOROG_GRAD_Y   = SI( 6,Sect_No,im_index)
      JOROG_UNFILT   = SI( 7,Sect_No,im_index)
      JOROG_GRAD_XX  = SI(35,Sect_No,im_index)
      JOROG_GRAD_XY  = SI(36,Sect_No,im_index)
      JOROG_GRAD_YY  = SI(37,Sect_No,im_index)

! Sea/Sea Ice fields
      JU_SEA         = SI( 28,Sect_No,im_index)
      JV_SEA         = SI( 29,Sect_No,im_index)
      JICE_FRACTION  = SI( 31,Sect_No,im_index)
      JICE_THICKNESS = SI( 32,Sect_No,im_index)
      JTI            = SI( 49,Sect_No,im_index)
      JICE_FRACT_CAT = SI(413,Sect_No,im_index)
      JICE_THICK_CAT = SI(414,Sect_No,im_index)
      JTI_CAT        = SI(415,Sect_No,im_index)
      JICE_K_CAT     = SI(440,Sect_No,im_index)
      JU_0_P         = SI(269,Sect_No,im_index)
      JV_0_P         = SI(270,Sect_No,im_index)

! Snow fields
      JSNODEP        = SI( 23,Sect_No,im_index) ! Snow depth over land
      JSNODEP_SEA    = SI( 95,Sect_No,im_index) ! Snow depth on sea ice
      JSNODEP_SEA_CAT= SI(416,Sect_No,im_index) ! Snow depth on ice cats
      JSNSOOT        = SI(221,Sect_No,im_index) ! Snow soot content
      JCATCH_SNOW    = SI(241,Sect_No,im_index)
      JSNOW_GRND     = SI(242,Sect_No,im_index)

! Decoupled screen temperatures
      JTScrnDcl_TILE = SI(490,Sect_No,im_index) ! Decoupled screen-level
                                                ! temperature on tiles
      JTScrnDcl_SSI  = SI(491,Sect_No,im_index) ! Decoupled screen-level
                                                ! temperature on sea/s.ice
      JtStbTrans     = SI(492,Sect_No,im_index) ! Time since the transition

! Ozone
      JOZONE(o3dims2%k_start)     = SI(60,Sect_No,im_index)
! Check for zonal ozone and calculate pointers accordingly
      LEXPAND_OZONE=.FALSE.
      IF (A_LOOKUP(LBNPT,PPINDEX(60,im_index)) == 1) THEN
        LEXPAND_OZONE = .TRUE.
      ENDIF

      DO LEV= o3dims2%k_start+1, o3dims2%k_end
        IF(LEXPAND_OZONE) THEN
!         Ozone held as zonal averages, i.e. one value per row
          JOZONE(LEV)=JOZONE(LEV-1)+ROWS
        ELSE
          JOZONE(LEV)=JOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

!! Tropopause-based Ozone
      IF (tpps_ozone_levels >  0) THEN
        JTPPSOZONE(1)     = SI(341,Sect_No,im_index)

        !Check for zonal tpps_ozone and calculate pointers accordingly
        LEXPAND_TPPS_OZONE=.FALSE.
        IF (A_LOOKUP(LBNPT,PPINDEX(341,im_index)) == 1) THEN
          LEXPAND_TPPS_OZONE = .TRUE.
        ENDIF
      END IF

      DO LEV=2,TPPS_OZONE_LEVELS
        IF(LEXPAND_TPPS_OZONE) THEN
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+ROWS
        ELSE
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

! Add prognostic ozone tracer and cariolle parameters to section 0
! 
        JOZONE_TRACER = SI(480,sect_no,im_index)
        JO3_PROD_LOSS = SI(481,sect_no,im_index)
        JO3_P_L_VMR   = SI(482,sect_no,im_index)
        JO3_VMR       = SI(483,sect_no,im_index)
        JO3_P_L_TEMP  = SI(484,sect_no,im_index)
        JO3_TEMP      = SI(485,sect_no,im_index)
        JO3_P_L_COLO3 = SI(486,sect_no,im_index)
        JO3_COLO3     = SI(487,sect_no,im_index)

        DO  LEV = tdims%k_start+1, tdims%k_end
          JOZONE_TRACER(LEV) = JOZONE_TRACER(LEV-1) + THETA_OFF_SIZE
          JO3_PROD_LOSS(LEV) = JO3_PROD_LOSS(LEV-1) + ROWS
          JO3_P_L_VMR(LEV)   = JO3_P_L_VMR(LEV-1) + ROWS
          JO3_VMR(LEV)       = JO3_VMR(LEV-1) + ROWS
          JO3_P_L_TEMP(LEV)  = JO3_P_L_TEMP(LEV-1) + ROWS 
          JO3_TEMP(LEV)      = JO3_TEMP(LEV-1) + ROWS
          JO3_P_L_COLO3(LEV) = JO3_P_L_COLO3(LEV-1) + ROWS
          JO3_COLO3(LEV)     = JO3_COLO3(LEV-1) + ROWS
        END DO


! STOCHEM fields - removed

! Add sources and aerosol ancillaries

      JMURK_SOURCE(tdims%k_start)= SI(57,Sect_No,im_index) !Murk source
      JSO2_EM                    = SI(58,Sect_No,im_index) !Sulphur dioxide emiss.
      JDMS_EM                    = SI(59,Sect_No,im_index) !Dimethyl sulphide emiss.
      JMURK       (tdims%k_start)= SI(90,Sect_No,im_index) !Murk concentration

! Add for Sulphur Cycle
      JSO2       (tdims%k_start)= SI(101,Sect_No,im_index) !Sulphur dioxide gas
      JDMS       (tdims%k_start)= SI(102,Sect_No,im_index) !Dimethyl sulphide gas
      JSO4_AITKEN(tdims%k_start)= SI(103,Sect_No,im_index) !Aitken mode SO4 aerosol
      JSO4_ACCU  (tdims%k_start)= SI(104,Sect_No,im_index) !Accumulation mode SO4 aer
      JSO4_DISS  (tdims%k_start)= SI(105,Sect_No,im_index) !Dissolved SO4 aerosol
      JH2O2      (tdims%k_start)= SI(106,Sect_No,im_index) !Hydrogen peroxide mmr
      JNH3       (tdims%k_start)= SI(107,Sect_No,im_index) !Ammonia gas
      JSOOT_NEW  (tdims%k_start)= SI(108,Sect_No,im_index) !Fresh soot
      JSOOT_AGD  (tdims%k_start)= SI(109,Sect_No,im_index) !Aged soot
      JSOOT_CLD  (tdims%k_start)= SI(110,Sect_No,im_index) !Soot in cloud
      JBMASS_NEW (tdims%k_start)= SI(111,Sect_No,im_index) !Fresh biomass smoke
      JBMASS_AGD (tdims%k_start)= SI(112,Sect_No,im_index) !Aged biomass smoke
      JBMASS_CLD (tdims%k_start)= SI(113,Sect_No,im_index) !Cloud biomass smoke
      JOCFF_NEW  (tdims%k_start)= SI(114,Sect_No,im_index) !Fresh ocff
      JOCFF_AGD  (tdims%k_start)= SI(115,Sect_No,im_index) !Aged socff
      JOCFF_CLD  (tdims%k_start)= SI(116,Sect_No,im_index) !Ocff in cloud
      JSO2_NATEM (tdims%k_start)= SI(121,Sect_No,im_index) !Natural SO2 emissions
      JOH        (tdims%k_start)= SI(122,Sect_No,im_index) !OH 3_D ancillary
      JHO2       (tdims%k_start)= SI(123,Sect_No,im_index) !HO2 3_D ancillary
      JH2O2_LIMIT(tdims%k_start)= SI(124,Sect_No,im_index) !H2O2 LIMIT 3_D ancillary
      JO3_CHEM   (tdims%k_start)= SI(125,Sect_No,im_index) !O3 for chemistry 3_D anc
      JSO2_HILEM    =SI(126,Sect_No,im_index)  !High level SO2 emissions
      JNH3_EM       =SI(127,Sect_No,im_index)  !Ammonia surface emiss
      JSOOT_EM      =SI(128,Sect_No,im_index)  !Fresh soot surf emiss
      JSOOT_HILEM   =SI(129,Sect_No,im_index)  !Fresh soot high emiss
      JBMASS_EM     =SI(130,Sect_No,im_index)  !Fresh bmass surf emiss
      JBMASS_HILEM  =SI(131,Sect_No,im_index)  !Elevated bmass emiss
      JDMS_CONC     =SI(132,Sect_No,im_index)  !DMS conc in seawater
      JDMS_OFLUX    =SI(133,Sect_No,im_index)  !DMS flux from ocean
      JOCFF_EM      =SI(134,Sect_No,im_index)  !Fresh OCFF surf emiss
      JOCFF_HILEM   =SI(135,Sect_No,im_index)  !Fresh OCFF high emiss

! Aerosol climatologies
      JARCLBIOG_BG  =SI(351,Sect_No,im_index)  ! Biogenic aerosol climatology
      JARCLBIOM_FR  =SI(352,Sect_No,im_index)  ! Biomass burning (fresh) aerosol clim
      JARCLBIOM_AG  =SI(353,Sect_No,im_index)  ! Biomass burning (aged) aerosol clim
      JARCLBIOM_IC  =SI(354,Sect_No,im_index)  ! Biomass burning (in-cloud) aerosol clim
      JARCLBLCK_FR  =SI(355,Sect_No,im_index)  ! Black carbon (fresh) aerosol clim
      JARCLBLCK_AG  =SI(356,Sect_No,im_index)  ! Black carbon (aged) aerosol clim
      JARCLSSLT_FI  =SI(357,Sect_No,im_index)  ! Sea salt (film mode) aerosol clim 
      JARCLSSLT_JT  =SI(358,Sect_No,im_index)  ! Sea salt (jet mode) aerosol clim
      JARCLSULP_AC  =SI(359,Sect_No,im_index)  ! Sulphate (accumulation mode) aero clim
      JARCLSULP_AK  =SI(360,Sect_No,im_index)  ! Sulphate (Aitken mode) aerosol clim 
      JARCLSULP_DI  =SI(361,Sect_No,im_index)  ! Sulphate (dissolved) aerosol clim
      JARCLDUST_B1  =SI(362,Sect_No,im_index)  ! Dust (bin 1) aerosol climatology 
      JARCLDUST_B2  =SI(363,Sect_No,im_index)  ! Dust (bin 2) aerosol climatology 
      JARCLDUST_B3  =SI(364,Sect_No,im_index)  ! Dust (bin 3) aerosol climatology 
      JARCLDUST_B4  =SI(365,Sect_No,im_index)  ! Dust (bin 4) aerosol climatology 
      JARCLDUST_B5  =SI(366,Sect_No,im_index)  ! Dust (bin 5) aerosol climatology 
      JARCLDUST_B6  =SI(367,Sect_No,im_index)  ! Dust (bin 6) aerosol climatology 
      JARCLOCFF_FR  =SI(368,Sect_No,im_index)  ! Org carbon fossil fuel (fresh) aero clim
      JARCLOCFF_AG  =SI(369,Sect_No,im_index)  ! Org carbon fossil fuel (aged) aero clim
      JARCLOCFF_IC  =SI(370,Sect_No,im_index)  ! Org carbon fossil fuel (in-cloud) aero clim
      JARCLDLTA_DL  =SI(371,Sect_No,im_index)  ! Delta aerosol climatology
      
      DO LEV= tdims%k_start+1, tdims%k_end
        JARCLBIOG_BG(LEV) = JARCLBIOG_BG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_FR(LEV) = JARCLBIOM_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_AG(LEV) = JARCLBIOM_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_IC(LEV) = JARCLBIOM_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_FR(LEV) = JARCLBLCK_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_AG(LEV) = JARCLBLCK_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_FI(LEV) = JARCLSSLT_FI(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_JT(LEV) = JARCLSSLT_JT(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AC(LEV) = JARCLSULP_AC(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AK(LEV) = JARCLSULP_AK(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_DI(LEV) = JARCLSULP_DI(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B1(LEV) = JARCLDUST_B1(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B2(LEV) = JARCLDUST_B2(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B3(LEV) = JARCLDUST_B3(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B4(LEV) = JARCLDUST_B4(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B5(LEV) = JARCLDUST_B5(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B6(LEV) = JARCLDUST_B6(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_FR(LEV) = JARCLOCFF_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_AG(LEV) = JARCLOCFF_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_IC(LEV) = JARCLOCFF_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLDLTA_DL(LEV) = JARCLDLTA_DL(LEV-1)+THETA_FIELD_SIZE
      END DO

! Mineral dust scheme

      JSOIL_CLAY   =SI(418,Sect_No,im_index)  ! soil clay fraction
      JSOIL_SILT   =SI(419,Sect_No,im_index)  ! soil silt fraction
      JSOIL_SAND   =SI(420,Sect_No,im_index)  ! soil sand fraction

      JDUST_MREL1=SI(421,Sect_No,im_index) !relative soil mass in div1
      JDUST_MREL2=SI(422,Sect_No,im_index) !relative soil mass in div2
      JDUST_MREL3=SI(423,Sect_No,im_index) !relative soil mass in div3
      JDUST_MREL4=SI(424,Sect_No,im_index) !relative soil mass in div4
      JDUST_MREL5=SI(425,Sect_No,im_index) !relative soil mass in div5
      JDUST_MREL6=SI(426,Sect_No,im_index) !relative soil mass in div6


      JDUST_DIV1(tdims%k_start)=SI(431,Sect_No,im_index)  ! dust mmr, division 1
      JDUST_DIV2(tdims%k_start)=SI(432,Sect_No,im_index)  ! dust mmr, division 2
      JDUST_DIV3(tdims%k_start)=SI(433,Sect_No,im_index)  ! dust mmr, division 3
      JDUST_DIV4(tdims%k_start)=SI(434,Sect_No,im_index)  ! dust mmr, division 4
      JDUST_DIV5(tdims%k_start)=SI(435,Sect_No,im_index)  ! dust mmr, division 5
      JDUST_DIV6(tdims%k_start)=SI(436,Sect_No,im_index)  ! dust mmr, division 6

      DO LEV = tdims%k_start+1, tdims%k_end
       JDUST_DIV1(LEV)=JDUST_DIV1(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV2(LEV)=JDUST_DIV2(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV3(LEV)=JDUST_DIV3(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV4(LEV)=JDUST_DIV4(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV5(LEV)=JDUST_DIV5(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV6(LEV)=JDUST_DIV6(LEV-1)+THETA_OFF_SIZE
      ENDDO

! Ammonium nitrate scheme
      JNITR_ACC (tdims%k_start)= SI(117,Sect_No,im_index) !Accumulation nitrate MMR
      JNITR_DISS(tdims%k_start)= SI(118,Sect_No,im_index) !Dissolved nitrate MMR 
      DO LEV = tdims%k_start+1, tdims%k_end 
        JNITR_ACC(LEV) = JNITR_ACC(LEV-1)+THETA_OFF_SIZE  
        JNITR_DISS(LEV) = JNITR_DISS(LEV-1)+THETA_OFF_SIZE  
      ENDDO 

! HadCM2 sulphate loading patterns
      JHadCM2_SO4(1)= SI(160,Sect_No,im_index)
      DO LEV=2, NSULPAT  ! nsulpat hard wired to 2 in clmchfcg_scenario_mod
        JHadCM2_SO4(LEV)=JHadCM2_SO4(LEV-1)+THETA_FIELD_SIZE
      ENDDO

! Add for Carbon cycle
      J_CO2FLUX = SI(250,Sect_No,im_index)
      J_CO2_EMITS  = SI(251,Sect_No,im_index)
      JCO2(tdims%k_start)      = SI(252,Sect_No,im_index)
      DO LEV= tdims%k_start+1, tdims%k_end
        JMURK_SOURCE(LEV) = JMURK_SOURCE(LEV-1)+THETA_FIELD_SIZE
        JMURK(LEV) = JMURK(LEV-1)+THETA_OFF_SIZE

! For Sulphur Cycle variables
        JSO2(LEV)=JSO2(LEV-1)+THETA_OFF_SIZE
        JDMS(LEV)=JDMS(LEV-1)+THETA_OFF_SIZE
        JSO4_AITKEN(LEV)=JSO4_AITKEN(LEV-1)+THETA_OFF_SIZE
        JSO4_ACCU(LEV)=JSO4_ACCU(LEV-1)+THETA_OFF_SIZE
        JSO4_DISS(LEV)=JSO4_DISS(LEV-1)+THETA_OFF_SIZE
        JH2O2(LEV)=JH2O2(LEV-1)+THETA_OFF_SIZE
        JSO2_NATEM(LEV)=JSO2_NATEM(LEV-1)+THETA_FIELD_SIZE
        JOH(LEV) = JOH(LEV-1)+THETA_FIELD_SIZE
        JHO2(LEV) = JHO2(LEV-1)+THETA_FIELD_SIZE
        JNH3(LEV)      = JNH3(LEV-1)+THETA_OFF_SIZE
        JSOOT_NEW(LEV) = JSOOT_NEW(LEV-1)+THETA_OFF_SIZE
        JSOOT_AGD(LEV) = JSOOT_AGD(LEV-1)+THETA_OFF_SIZE
        JSOOT_CLD(LEV) = JSOOT_CLD(LEV-1)+THETA_OFF_SIZE
        JBMASS_NEW(LEV) = JBMASS_NEW(LEV-1)+THETA_OFF_SIZE
        JBMASS_AGD(LEV) = JBMASS_AGD(LEV-1)+THETA_OFF_SIZE
        JBMASS_CLD(LEV) = JBMASS_CLD(LEV-1)+THETA_OFF_SIZE
        JOCFF_NEW(LEV) = JOCFF_NEW(LEV-1)+THETA_OFF_SIZE
        JOCFF_AGD(LEV) = JOCFF_AGD(LEV-1)+THETA_OFF_SIZE
        JOCFF_CLD(LEV) = JOCFF_CLD(LEV-1)+THETA_OFF_SIZE
        JH2O2_LIMIT(LEV)=JH2O2_LIMIT(LEV-1)+THETA_FIELD_SIZE
        JO3_CHEM(LEV)=JO3_CHEM(LEV-1)+THETA_FIELD_SIZE

! For Carbon Cycle variables
        JCO2(LEV)=JCO2(LEV-1)+THETA_OFF_SIZE
      END DO
      !CABLE: pointers into d1 (e.g. jTSOIL_TILE) included via argptra.h etc
      call cable_set_atm_pointers( SI, NITEMS, NSECTS, N_INTERNAL_MODEL, &
                                 Sect_No,im_index, & 
                                 jTSOIL_TILE, jSMCL_TILE, jSTHF_TILE,            &
                                 jSNOW_DEPTH3L, jSNOW_MASS3L, jSNOW_TMP3L,       &
                                 jSNOW_RHO3L, jSNOW_RHO1L, jSNOW_AGE, jSNOW_FLG3L )

      JUSER_ANC19 = SI(319,Sect_No,im_index)
      JUSER_ANC20 = SI(320,Sect_No,im_index)
      JUSER_MULT1(1)  = SI(321,Sect_No,im_index)
      JUSER_MULT2(1)  = SI(322,Sect_No,im_index)
      JUSER_MULT18(1) = SI(338,Sect_No,im_index)
      JUSER_MULT19(1) = SI(339,Sect_No,im_index)
      JUSER_MULT20(1) = SI(340,Sect_No,im_index)

! Set for multi-level user ancillaries
      DO LEV=2,MODEL_LEVELS
        JUSER_MULT1(LEV)  = JUSER_MULT1(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT2(LEV)  = JUSER_MULT2(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT18(LEV) = JUSER_MULT18(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT19(LEV) = JUSER_MULT19(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT20(LEV) = JUSER_MULT20(LEV-1)+THETA_FIELD_SIZE
      END DO

! Tiled vegetation and triffid
      JFRAC_TYP     = SI(216,Sect_No,im_index) ! surface type fractions
      JFRAC_CON1    = SI(442,Sect_No,im_index) ! surface type fractions
      JFRAC_CON2    = SI(443,Sect_No,im_index) ! surface type fractions
      JFRAC_CON3    = SI(444,Sect_No,im_index) ! surface type fractions
      JFRAC_CON4    = SI(445,Sect_No,im_index) ! surface type fractions
      JFRAC_CON5    = SI(446,Sect_No,im_index) ! surface type fractions
      JFRAC_CON6    = SI(447,Sect_No,im_index) ! surface type fractions
      JFRAC_CON7    = SI(448,Sect_No,im_index) ! surface type fractions
      JFRAC_CON8    = SI(449,Sect_No,im_index) ! surface type fractions
      JFRAC_CON9    = SI(450,Sect_No,im_index) ! surface type fractions
      JLAI_PFT      = SI(217,Sect_No,im_index) ! leaf area index of PFTs
      JCANHT_PFT    = SI(218,Sect_No,im_index) ! canopy height of PFTs
      JDISTURB      = SI(219,Sect_No,im_index) ! Veg disturbed fraction
      JSOIL_ALB     = SI(220,Sect_No,im_index) ! Snow-free soil albedo
      JOBS_ALB_SW   = SI(243,Sect_No,im_index) ! Observed snow-free SW albedo 
      JOBS_ALB_VIS  = SI(244,Sect_No,im_index) ! Observed snow-free VIS albedo 
      JOBS_ALB_NIR  = SI(245,Sect_No,im_index) ! Observed snow-free NIR albedo 
      JSOIL_CARB    = SI(223,Sect_No,im_index) ! Soil carbon content
      JSOIL_CARB1   = SI(466,Sect_No,im_index) ! Soil carbon content DPM
      JSOIL_CARB2   = SI(467,Sect_No,im_index) ! Soil carbon content RPM
      JSOIL_CARB3   = SI(468,Sect_No,im_index) ! Soil carbon content BIO
      JSOIL_CARB4   = SI(469,Sect_No,im_index) ! Soil carbon content HUM
      JNPP_PFT_ACC  = SI(224,Sect_No,im_index) ! Accumulated NPP on PFTs
      JG_LF_PFT_ACC = SI(225,Sect_No,im_index) ! Accumulated leaf
!                                              ! turnover rate on PFTs
      JG_PHLF_PFT_ACC=SI(226,Sect_No,im_index) ! Accumulat. phenological
!                                              ! leaf turnover rate PFTs
      JRSP_W_PFT_ACC= SI(227,Sect_No,im_index) ! Accum. wood resp PFTs
      JRSP_S_ACC    = SI(228,Sect_No,im_index) ! Accumulated soil resp
      JRSP_S_ACC1 = SI(470,Sect_No,im_index)  ! Soil respiration DPM
      JRSP_S_ACC2 = SI(471,Sect_No,im_index)  ! Soil respiration RPM
      JRSP_S_ACC3 = SI(472,Sect_No,im_index)  ! Soil respiration BIO
      JRSP_S_ACC4 = SI(473,Sect_No,im_index)  ! Soil respiration HUM
      JCAN_WATER_TILE=SI(229,Sect_No,im_index) ! Canopy water content
!                                              ! on tiles
      JCATCH_TILE   = SI(230,Sect_No,im_index) ! Canopy capacity on
!                                              ! tiles
      JRGRAIN_TILE  = SI(231,Sect_No,im_index) ! Snow grain size on
!                                              ! tiles
      JTSTAR_TILE   = SI(233,Sect_No,im_index) ! Tiled surface temp
      JZ0_TILE      = SI(234,Sect_No,im_index) ! Tiled surface roughness
      JZ0H_TILE     = SI(246,Sect_No,im_index) ! Tiled surface thermal 
!                                              ! roughness
! Stash number for snow tile not finalised yet
      JSNODEP_TILE  = SI(240,Sect_No,im_index) ! Tiled snow depth
      JINFIL_TILE   = SI(236,Sect_No,im_index) ! Max tile infilt rate
! Stash codes for DORL, LW_DOWN, SW_TILE not finalised yet
      JDOLR         = SI(239,Sect_No,im_index) ! TOA surface up LW
      JLW_DOWN      = SI(238,Sect_No,im_index) ! Surface down LW
      JSW_TILE      = SI(237,Sect_No,im_index) ! Surface net SW on tiles

! MORUSES urban scheme
      JURBHGT             = SI(494,Sect_No,im_index) ! Building height
      JURBHWR             = SI(495,Sect_No,im_index) ! Height to width
      JURBWRR             = SI(496,Sect_No,im_index) ! Width ratio
      JURBDISP            = SI(497,Sect_No,im_index) ! Displacement height
      JURBZTM             = SI(498,Sect_No,im_index) ! Effective roughness
                                                     ! length for momentum
      JURBALBWL           = SI(499,Sect_No,im_index) ! Wall albedo
      JURBALBRD           = SI(500,Sect_No,im_index) ! Road albedo
      JURBEMISW           = SI(501,Sect_No,im_index) ! Wall emmissivity
      JURBEMISR           = SI(502,Sect_No,im_index) ! Road emmissivity

! River routing fields
      JRIV_SEQUENCE  = SI(151,Sect_No,im_index) ! River sequence
      JRIV_DIRECTION = SI(152,Sect_No,im_index) ! River Direction
      JRIV_STORAGE   = SI(153,Sect_No,im_index) ! River Water Storage
      JTOT_SURFROFF  = SI(155,Sect_No,im_index) ! Acc. surface runoff
      JTOT_SUBROFF   = SI(156,Sect_No,im_index) ! Acc. sub-surf runoff
! Set pointer for inland basin outflow
      JRIV_INLANDATM    = SI(511,Sect_No,im_index)
       !Inland basin outflow



! required for energy correction
      JNET_FLUX=SI(222,Sect_No,im_index)    ! store for energy flux
      JNET_MFLUX=SI(235,Sect_No,im_index)   ! store for moisture flux

!Fields carried forward from previous version.
      JTSTAR_ANOM    = SI(39,Sect_No,im_index)

! JULES version 2 prognostics
      JSNOWDEPTH      = SI(376,Sect_No,im_index) ! Snow depth on ground on tiles (m)
      JRHO_SNOW_GRND  = SI(377,Sect_No,im_index) ! Snowpack bulk density (kg/m3)
      JNSNOW          = SI(380,Sect_No,im_index) ! Number of snow layers on ground on tiles
      JDS             = SI(381,Sect_No,im_index) ! Snow layer thickness (m)
      JSICE           = SI(382,Sect_No,im_index) ! Snow layer ice mass on tiles (Kg/m2)
      JSLIQ           = SI(383,Sect_No,im_index) ! Snow layer liquid mass on tiles (Kg/m2)
      JTSNOWLAYER     = SI(384,Sect_No,im_index) ! Snow layer temperature (K)
      JRHO_SNOW       = SI(385,Sect_No,im_index) ! Snow layer densities (kg/m3)
      JRGRAINL        = SI(386,Sect_No,im_index) ! Snow layer grain size on tiles (microns)

! FLake lake scheme prognostics
      JLAKE_DEPTH  = SI(291,Sect_No,im_index) ! lake depth (m)
      JLAKE_FETCH  = SI(292,Sect_No,im_index) ! typical wind fetch (m)
      JLAKE_T_MEAN = SI(293,Sect_No,im_index) ! lake mean temperature (K)
      JLAKE_T_MXL  = SI(294,Sect_No,im_index) ! lake mixed-layer temperature (K)
      JLAKE_T_ICE  = SI(295,Sect_No,im_index) ! temperature at upper boundary of lake ice (K)
      JLAKE_H_MXL  = SI(296,Sect_No,im_index) ! lake mixed-layer depth (m)
      JLAKE_H_ICE  = SI(297,Sect_No,im_index) ! lake ice thickness (m)
      JLAKE_SHAPE  = SI(298,Sect_No,im_index) ! thermocline shape factor
      JLAKE_G_DT   = SI(299,Sect_No,im_index) ! lake ht.flx / dT (W m-2 K-1)

! Set all pointers referencing Lateral Boundary Conditions

      Sect_No=31  ! LBC section

      JOROG_LBC     = SI(1,Sect_No,im_index)
      JU_LBC        = SI(2,Sect_No,im_index)
      JV_LBC        = SI(3,Sect_No,im_index)
      JW_LBC        = SI(4,Sect_No,im_index)
      JRHO_LBC      = SI(5,Sect_No,im_index)
      JTHETA_LBC    = SI(6,Sect_No,im_index)
      JQ_LBC        = SI(7,Sect_No,im_index)
      JQCL_LBC      = SI(8,Sect_No,im_index)
      JQCF_LBC      = SI(9,Sect_No,im_index)
      JEXNER_LBC    = SI(10,Sect_No,im_index)
      JU_ADV_LBC    = SI(11,Sect_No,im_index)
      JV_ADV_LBC    = SI(12,Sect_No,im_index)
      JW_ADV_LBC    = SI(13,Sect_No,im_index)
      JQCF2_LBC     = SI(14,Sect_No,im_index)
      JQRAIN_LBC    = SI(15,Sect_No,im_index)
      JQGRAUP_LBC   = SI(16,Sect_No,im_index)
      JCF_BULK_LBC  = SI(17,Sect_No,im_index)
      JCF_LIQUID_LBC= SI(18,Sect_No,im_index)
      JCF_FROZEN_LBC= SI(19,Sect_No,im_index)
      JMURK_LBC     = SI(20,Sect_No,im_index)

      JDUST_DIV1_LBC = SI(23,Sect_No,im_index)
      JDUST_DIV2_LBC = SI(24,Sect_No,im_index)
      JDUST_DIV3_LBC = SI(25,Sect_No,im_index)
      JDUST_DIV4_LBC = SI(26,Sect_No,im_index)
      JDUST_DIV5_LBC = SI(27,Sect_No,im_index)
      JDUST_DIV6_LBC = SI(28,Sect_No,im_index)
      JSO2_LBC       = SI(29,Sect_No,im_index)
      JDMS_LBC       = SI(30,Sect_No,im_index)
      JSO4_AITKEN_LBC= SI(31,Sect_No,im_index)
      JSO4_ACCU_LBC  = SI(32,Sect_No,im_index)
      JSO4_DISS_LBC  = SI(33,Sect_No,im_index)
      JNH3_LBC       = SI(35,Sect_No,im_index)
      JSOOT_NEW_LBC  = SI(36,Sect_No,im_index)
      JSOOT_AGD_LBC  = SI(37,Sect_No,im_index)
      JSOOT_CLD_LBC  = SI(38,Sect_No,im_index)
      JBMASS_NEW_LBC = SI(39,Sect_No,im_index)
      JBMASS_AGD_LBC = SI(40,Sect_No,im_index)
      JBMASS_CLD_LBC = SI(41,Sect_No,im_index)
      JOCFF_NEW_LBC  = SI(42,Sect_No,im_index)
      JOCFF_AGD_LBC  = SI(43,Sect_No,im_index)
      JOCFF_CLD_LBC  = SI(44,Sect_No,im_index)
      JNITR_ACC_LBC  = SI(45,Sect_No,im_index)
      JNITR_DISS_LBC = SI(46,Sect_No,im_index)

      JU_LBC_TEND     = SI(257,Sect_No,im_index)
      JV_LBC_TEND     = SI(258,Sect_No,im_index)
      JW_LBC_TEND     = SI(259,Sect_No,im_index)
      JRHO_LBC_TEND   = SI(260,Sect_No,im_index)
      JTHETA_LBC_TEND = SI(261,Sect_No,im_index)
      JQ_LBC_TEND     = SI(262,Sect_No,im_index)
      JQCL_LBC_TEND   = SI(263,Sect_No,im_index)
      JQCF_LBC_TEND   = SI(264,Sect_No,im_index)
      JEXNER_LBC_TEND = SI(265,Sect_No,im_index)
      JU_ADV_LBC_TEND = SI(266,Sect_No,im_index)
      JV_ADV_LBC_TEND = SI(267,Sect_No,im_index)
      JW_ADV_LBC_TEND = SI(268,Sect_No,im_index)
      JQCF2_LBC_TEND   = SI(269,Sect_No,im_index)
      JQRAIN_LBC_TEND  = SI(270,Sect_No,im_index)
      JQGRAUP_LBC_TEND = SI(271,Sect_No,im_index)
      JCF_BULK_LBC_TEND  = SI(272,Sect_No,im_index)
      JCF_LIQUID_LBC_TEND= SI(273,Sect_No,im_index)
      JCF_FROZEN_LBC_TEND= SI(274,Sect_No,im_index)
      JMURK_LBC_TEND   = SI(275,Sect_No,im_index)

      JDUST_DIV1_LBC_TEND = SI(276,Sect_No,im_index)
      JDUST_DIV2_LBC_TEND = SI(277,Sect_No,im_index)
      JDUST_DIV3_LBC_TEND = SI(278,Sect_No,im_index)
      JDUST_DIV4_LBC_TEND = SI(279,Sect_No,im_index)
      JDUST_DIV5_LBC_TEND = SI(280,Sect_No,im_index)
      JDUST_DIV6_LBC_TEND = SI(281,Sect_No,im_index)
      JSO2_LBC_TEND       = SI(282,Sect_No,im_index)
      JDMS_LBC_TEND       = SI(283,Sect_No,im_index)
      JSO4_AITKEN_LBC_TEND= SI(284,Sect_No,im_index)
      JSO4_ACCU_LBC_TEND  = SI(285,Sect_No,im_index)
      JSO4_DISS_LBC_TEND  = SI(286,Sect_No,im_index)
      JNH3_LBC_TEND       = SI(288,Sect_No,im_index)
      JSOOT_NEW_LBC_TEND  = SI(289,Sect_No,im_index)
      JSOOT_AGD_LBC_TEND  = SI(290,Sect_No,im_index)
      JSOOT_CLD_LBC_TEND  = SI(291,Sect_No,im_index)
      JBMASS_NEW_LBC_TEND = SI(292,Sect_No,im_index)
      JBMASS_AGD_LBC_TEND = SI(293,Sect_No,im_index)
      JBMASS_CLD_LBC_TEND = SI(294,Sect_No,im_index)
      JOCFF_NEW_LBC_TEND  = SI(295,Sect_No,im_index)
      JOCFF_AGD_LBC_TEND  = SI(296,Sect_No,im_index)
      JOCFF_CLD_LBC_TEND  = SI(297,Sect_No,im_index)
      JNITR_ACC_LBC_TEND  = SI(298,Sect_No,im_index)
      JNITR_DISS_LBC_TEND = SI(299,Sect_No,im_index)

! Set pointers for Tracer prognostics
! Tracer prognostics are now in section 33, not section 0
      sect_no = 33   ! tracers section
      JVAR=0         ! JVAR+1 is the current tracer to be found
      IF (TR_VARS >  0) THEN
        DO IVAR=A_TRACER_FIRST,A_TRACER_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTRACER(trdims_xstl%k_start,JVAR) = SI(IVAR,sect_no,im_index)
            ! Set up array containing stash item codes for active 
            ! tracer prognostics (1:TR_VARS)
            A_TR_StashItem(JVAR) = IVAR
            DO LEV = trdims_xstl%k_start+1, trdims_xstl % k_end
              JTRACER(LEV,JVAR)=JTRACER(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_TR_INDEX(IVAR-A_TRACER_FIRST+1)=JVAR
          ELSE
            ! If tracer not active, set value to -1
            A_TR_INDEX(IVAR-A_TRACER_FIRST+1) = -1
          END IF
        END DO
      ELSE   
        ! Ensure a sensible address even if no tracers
        JTRACER(trdims_xstl%k_start,1)=1
      ENDIF
      IF(JVAR /= TR_VARS) THEN
        WRITE(6,*) 'STATMPT: TR_VARS and SI are inconsistent'
        WRITE(6,*) 'TR_VARS=',TR_VARS,' .     But, SI implies :',JVAR
        CMESSAGE=  'STATMPT: TR_VARS and SI  inconsistent, see output'
        ICODE=100
        GOTO 9999 ! error return
      END IF

! UKCA tracer prognostics are in section 34.
      sect_no = ukca_sect   ! UKCA tracers section 
      JVAR=0         ! JVAR+1 is the current tracer to be found
      IF (TR_UKCA >  0) THEN
        DO IVAR=A_UKCA_FIRST,A_UKCA_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTR_UKCA(trdims_xstl%k_start,JVAR) = SI(IVAR,sect_no,im_index)
            ! Set up array containing stash item codes for active
            ! tracer prognostics (1:TR_UKCA)
            UKCA_TR_StashItem(JVAR) = IVAR
            DO LEV = trdims_xstl%k_start+1, trdims_xstl%k_end
              JTR_UKCA(LEV,JVAR)=JTR_UKCA(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_UKCA_INDEX(IVAR-A_UKCA_FIRST+1)=JVAR
          ELSE
            ! If tracer not active, set value to -1
            A_UKCA_INDEX(IVAR-A_UKCA_FIRST+1) = -1
          END IF
        END DO
      ELSE
        ! Ensure a sensible address when no tracers
        JTR_UKCA(trdims_xstl%k_start,1)=1
      ENDIF

      IF(JVAR /= TR_UKCA) THEN
        WRITE(6,*) 'STATMPT: TR_UKCA and SI are inconsistent'
        WRITE(6,*) 'TR_UKCA=',TR_UKCA,' .     But, SI implies :',JVAR
        CMESSAGE=  'STATMPT: TR_UKCA and SI  inconsistent, see output'
        ICODE=100
        GOTO 9999 ! error return
      END IF

      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         JO3_UKCA(tdims%k_start)   =SI(ukca_item_sulpc(1),sect_no,im_index)
         JHNO3_UKCA(tdims%k_start) =SI(ukca_item_sulpc(2),sect_no,im_index)
         JH2O2_UKCA(tdims%k_start) =SI(ukca_item_sulpc(3),sect_no,im_index)
         JOH_UKCA(tdims%k_start)   =SI(ukca_item_sulpc(4),sect_no,im_index) 
         JHO2_UKCA(tdims%k_start)  =SI(ukca_item_sulpc(5),sect_no,im_index) 
      END IF
      
! Set pointers for Tracer lateral boundary data

      ! Free tracer lbc data is in section 36
      ! The tracer lbc stash item codes in A_TR_LBC_StashItem
      ! are set up in INBOUNDA
 
      sect_no = 36   ! tracer lbcs section
      jvar = 0
      If (TR_LBC_VARS > 0) Then                                         
        Do ivar= A_TRACER_FIRST,A_TRACER_LAST
          If (SI(ivar,sect_no,im_index) /= 1) Then 
            jvar = jvar + 1
            JTRACER_LBC(jvar)        = SI(ivar,sect_no,im_index)
            JTRACER_LBC_TEND(jvar)   = SI(256+ivar,sect_no,im_index)
            A_TR_LBC_StashItem(jvar) = ivar
          End If
        End Do                                          
      Else   
        ! Ensure a sensible address even if no tracer lbcs   
        JTRACER_LBC(1)=1                                           
        JTRACER_LBC_TEND(1)=1 
      End If ! IF (TR_LBC_VARS > 0)                                     

      ! Set up array (1:TR_VARS) pointing to location in TRACER_LBC 
      ! array for each prognostic tracer (A_tr_active_lbc_index)
      ! This allows tracer prognostics to be active even if
      ! there are no lbcs, and lbc data to be ignored if the 
      ! tracer prognostic is not active. If there is no lbc data
      ! for a particular tracer, the array is set to -1.
       
      DO ivar= 1,TR_VARS
        A_tr_active_lbc_index(ivar) = -1 
        DO jvar= 1,TR_LBC_VARS
          IF (A_TR_LBC_StashItem(jvar) == A_TR_StashItem(ivar)) THEN
            A_tr_active_lbc_index(ivar) = jvar                  
            Write(6,*) 'A_tr_active_lbc_index:',ivar,jvar,              &
                  A_TR_LBC_StashItem(jvar),A_TR_StashItem(ivar),        &
                  A_tr_active_lbc_index(ivar)     
          END IF
        END DO
      END DO


      ! UKCA tracer lbc data is in section 37
      ! The tracer lbc stash item codes in UKCA_TR_LBC_StashItem
      ! are set up in INBOUNDA
 
      sect_no = 37   ! UKCA tracer lbcs section
      jvar = 0
      IF (TR_LBC_UKCA > 0) THEN                                         
        DO ivar= A_UKCA_FIRST,A_UKCA_LAST
          IF (SI(ivar,sect_no,im_index) /= 1) THEN 
            jvar = jvar + 1
            JTR_UKCA_LBC(jvar)        = SI(ivar,sect_no,im_index)
            JTR_UKCA_LBC_TEND(jvar)   = SI(256+ivar,sect_no,im_index)
            UKCA_TR_LBC_StashItem(jvar) = ivar
          END IF
        END DO                                          
      ELSE   
        ! Ensure a sensible address even if no tracer lbcs   
        JTR_UKCA_LBC(1)=1                                           
        JTR_UKCA_LBC_TEND(1)=1 
      END IF ! IF (TR_LBC_UKCA > 0)                                     

     ! Set up array (1:TR_UKCA) pointing to location in TRACER_LBC_UKCA 
     ! array for each prognostic tracer (UKCA_tr_active_lbc_index)
     ! This allows tracer prognostics to be active even if
     ! there are no lbcs, and lbc data to be ignored if the 
     ! tracer prognostic is not active. If there is no lbc data
     ! for a particular tracer, the array is set to -1.
       
      Do ivar= 1,TR_UKCA
        UKCA_tr_active_lbc_index(ivar) = -1 
        Do jvar= 1,TR_LBC_UKCA
          If (UKCA_TR_LBC_StashItem(jvar) == UKCA_TR_StashItem(ivar))   &
          Then
            UKCA_tr_active_lbc_index(ivar) = jvar                  
            Write(6,*) 'UKCA_tr_active_lbc_index:',ivar,jvar,           &
                  UKCA_TR_LBC_StashItem(jvar),UKCA_TR_StashItem(ivar),  &
                  UKCA_tr_active_lbc_index(ivar)     
          End If
        End Do
      End Do

! Set pointers to level dependent constants for atmosphere.
      JETATHETA      =1              ! eta_theta_levels(0:model_levels)
      JETARHO        =JETATHETA+model_levels+1 ! eta_rho_levels
      JRHCRIT        =JETARHO+model_levels+1   ! rhcrit
      JSOIL_THICKNESS=JRHCRIT+model_levels+1   ! soil level depths
! For height definition z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      Jzseak_theta   =JSOIL_THICKNESS+model_levels+1 ! zsea (theta levs)
      JCk_theta      =Jzseak_theta   +model_levels+1 ! C    (theta levs)
      Jzseak_rho     =JCk_theta      +model_levels+1 ! zsea (rho levs)
      JCk_rho        =Jzseak_rho     +model_levels+1 ! C    (rho levs)
      
      IF (A_LEN2_ROWDEPC > 0 .AND. A_LEN2_COLDEPC > 0) THEN
! Set pointers to Row dependent constants for atmosphere.
        JPHI_INPUT_P       =1
        JPHI_INPUT_V       =JPHI_INPUT_P + A_LEN1_ROWDEPC
! Set pointers to Col dependent constants for atmosphere.      
        JLAMBDA_INPUT_P    =1
        JLAMBDA_INPUT_U    =JLAMBDA_INPUT_P + A_LEN1_COLDEPC
      ELSE
        JPHI_INPUT_P       =1
        JPHI_INPUT_V       =1
        JLAMBDA_INPUT_P    =1
        JLAMBDA_INPUT_U    =1
      END IF

! The following if block should only be required for test purposes
! during the transition of vn5.2. This ensures old values for
! setting blev/bulev/brlev on old style dumps (before 'smooth'
! algorithm introduced and new definitions introduced).
      if(A_LEN2_LEVDEPC  <=  4) then ! ie before 5.2 change
         JCk_rho=jetarho
         Jck_theta=jetatheta
         Jzseak_rho=jetarho
         Jzseak_theta=jetatheta
         write(6,*) 'SETCONA_CTL: WARNING. This dump has not been '//   &
         'reconfigured with new 5.2 set of level dependent constants:'
         write(6,*) 'PP headers will revert to 5.0/5.1 style '//        &
         'blev/bulev/brlev definitions'
      endif

 9999 CONTINUE ! ERROR GOTO point.
      IF (lhook) CALL dr_hook('SET_ATM_POINTERS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SET_ATM_POINTERS
