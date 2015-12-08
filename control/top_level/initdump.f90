! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INITDUMP -------------------------------------------
!LL
!LL Purpose:  To read atmosphere dumps, and to calculate
!LL additional constants based on the dump header information.
!LL
!LL Extra constants needed for cloud types calulated within SETDCFLD
!LL
!LL Programming Standard : UM documentation paper no. 3
!LL
!LL                 U.M. Documentation paper no F3
!LL
!LLEND--------------------------------------------------------------
!
!*L Arguments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      SUBROUTINE INITDUMP(                                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
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
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                 sm_ident,ICODE,CMESSAGE)

      USE dynamics_input_mod, ONLY : l_endgame

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE model_file
      USE IO
      USE atm_fields_bounds_mod

      USE UM_ParVars
      USE Control_Max_Sizes
      USE Decomp_DB
      USE init_from_jules_namelists_mod, ONLY: init_from_jules_namelists
      USE nvegparm, ONLY : z0_nvg
      USE dust_parameters_mod, ONLY: z0_soil
      USE nstypes, ONLY : soil,npft
      USE lookup_addresses
      USE nlstgen_mod, ONLY: dump_packim
      USE Submodel_Mod


      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

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
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=1

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end
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

      INTEGER   sm_ident  ! Sub-model indicator
      INTEGER   ICODE     ! Return code
      INTEGER   NFTIN     ! FTN number for read
      INTEGER   NFTSWAP   ! FTN number for swapping radiation incrs
      INTEGER   SEGSTART  ! Pointer to start of radiation incrs
      INTEGER   SEGEND    ! Pointer to end of radiation incrs
      INTEGER   LEN_IO    ! Length of data transferred
      INTEGER   ERROR     ! Error code returned by OPEN

      REAL      A_IO      ! IO completion code

      CHARACTER(LEN=80)   CMESSAGE  ! Error message

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
! --------------------- Comdeck: CHISTORY ----------------------------
!
!  Purpose: COMMON block for history data needed by top level (C0)
!           routines, and passed from run to run.  Mostly set by
!           the User Interface.
!
!           Note that CHISTORY *CALLs ALL individual history comdecks
!
! --------------------------------------------------------------------
!
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations


      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: model_data_time(6)

      ! Indicator of operational run type
      INTEGER :: run_indic_op

      ! Final target date for the run
      INTEGER :: run_resubmit_target(6)

      ! Last field written/read per FT unit
      INTEGER :: ft_lastfield(20:nunits)

      ! Number of automatically-resubmitted job chunks
      ! Used to name output file
      INTEGER :: run_job_counter

! History Common Block for overall model integers variables.

      COMMON /ihisto/                                                 &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

      NAMELIST /nlihisto/                                             &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

  CHARACTER(LEN=8) ::  run_type             ! Type of run
  CHARACTER(LEN=1) ::  ft_active(20:nunits) ! "Y" if file partly written

  LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

  ! History Common Block for overall model character variables.

  COMMON /chisto/                                     &
     run_type,                                        &
     newrun, ft_active

  NAMELIST /nlchisto/                                 &
     run_type,                                        &
     ft_active

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! This file belongs in section: Top Level
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations
      ! History block copy of A_STEP held in file CTIME
      INTEGER :: h_stepim(n_internal_model_max)

      ! No of means activated
      INTEGER :: mean_offsetim(n_internal_model_max)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: offset_dumpsim(n_internal_model_max)

      ! No of mean periods chosen
      INTEGER :: mean_numberim(n_internal_model_max)

      ! Indicators used to correct logical units are used for
      ! atmos partial sum dump I/O
      INTEGER :: run_meanctl_indicim(4,n_internal_model_max)

      ! History Common Block for generic model integer variables.

      COMMON /ihistg/                                         &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         mean_numberim, run_meanctl_indicim

      NAMELIST /nlihistg/                                     &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         run_meanctl_indicim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character variables for
!              managing dump names
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 
!
!   Type declarations
!
! For keeping old restart dump name between creation of new restart dump
! and successful completion of climate means and history file.
CHARACTER(LEN=256) :: save_dumpname_im(n_internal_model_max)
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
! Blank name
CHARACTER(LEN=256) :: blank_file_name
!
! History Common Block for generic model characters variables.
!
COMMON /chistg/save_dumpname_im, checkpoint_dump_im, blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

! CHISTG end
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200

!
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
!
! Description:
!   This comdeck declares an integer variable 'ModelType' whose value
!   determines whether a model run is global, limited area or zonal.
!   The values of ModelType associated with each of the run types are
!   defined by integer parameters which are also declared below.
!    ModelType is set in subroutine SETLOGIC.
!
      ! Value used to represent the global model configuration
      INTEGER,PARAMETER:: GlobalModel=1

      ! Value used to represent the limited area model configuration
      INTEGER,PARAMETER:: LimitedAreaModel=2

      ! Value used to represent the 'periodic in x' model config
      ! INTEGER,PARAMETER:: ZonalModel=2

! Global scalars:
      INTEGER     ModelType  ! Integer switch which is equated to one
                             ! of the above parameters in a model run,
                             ! and so determines the configuration

! COMMON blocks:
      COMMON /RunType/ ModelType

! C_GLOBAL end
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
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
      LOGICAL                                                           &
         L_A_DUMP                                                       &
        ,L_A_PROG_ONLY 
                            ! ) Switches set if only prognostic
                            ! ) fields to be read in.

      INTEGER   i,j,ij    ! Loop counters
      INTEGER                                                           &
         LEN2_LOOKUP                                                    &
        ,LEN_DATA                                                       &
        ,N_PROG_FLDS                                                    &
                            ! No of prognostic fields in dumps
        ,N_PROG_LOOKUP                                                  &
        ,LEN_PROG                                                       &
        ,TOT_LEN_DATA                                                   &
        ,D1_ADDR_SUBMODEL_ID  ! submodel id in D1_ADDR array
      INTEGER :: A_MPP_ADDR(A_LEN2_LOOKUP),                             &
                 A_MPP_LEN(A_LEN2_LOOKUP)
      INTEGER info  ! return code from GC operations

      INTEGER                                                           &
              IFLD                                                      &
                            ! Loop variable
             ,IM_IDENT                                                  &
                            ! internal model identifier
             ,IM_INDEX      ! internal model index for STASH arrays

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER :: filelength

!*---------------------------------------------------------------------
!L   Internal Structure

      IF (lhook) CALL dr_hook('INITDUMP',zhook_in,zhook_handle)

!L 1.0 Read atmosphere dump and initialise atmosphere model.
      IF (sm_ident == atmos_sm) THEN

!L 1.1 Open unit for atmosphere dump, read fixed length header
!L     and set buffer length
        NFTIN = 21

      IF (H_STEPim(a_im) == 0) THEN
! In the first timestep, start from the ASTART dump
        CALL model_file_open(nftin,ft_environ(nftin),                          &
            len_ft_envir(nftin),ioOpenReadOnly,ioNameInEnv)

      ELSE
! Not the first timestep, so a continuation run. Restart from the last dump

! Store the name of the restart dump in save_dumpname_im. Required if job is
! set up to delete restart dump when next dump and history files have
! been safely written.
        filelength = LEN_TRIM(checkpoint_dump_im(a_im))
        save_dumpname_im(a_im) = checkpoint_dump_im(a_im)
        CALL model_file_open(nftin,checkpoint_dump_im(a_im),filelength,        &
            ioOpenReadOnly)
      END IF

! DEPENDS ON: read_flh
        CALL READ_FLH (NFTIN,A_FIXHD,LEN_FIXHD,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN 
          IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
          RETURN
        END IF


        CALL SETPOS (NFTIN,0,ICODE)

!       Test if atmos dump.
        L_A_DUMP = A_FIXHD(5) == 1 .AND. A_FIXHD(2) == atmos_sm

!       Test if only prognostic fields to be read in
        L_A_PROG_ONLY = L_A_DUMP .AND. H_STEPim(a_im) == 0

!       Get no of prognostic fields in atmos dump
        N_PROG_FLDS = A_FIXHD(153)

!       Check N_PROG_FLDS has been set.
        IF (N_PROG_FLDS == IMDI) THEN
          WRITE (6,*) '  '
          WRITE (6,*) ' No of prognostic fields not set in FIXHD(153)'
          WRITE (6,*) ' Run RECONFIGURATION to set FIXHD(153)'
          CMESSAGE = 'INITDUMP: FIXHD(153) not set in atmos dump'
          ICODE = 101
          GO TO 9999  !  Return
        END IF

!       Check N_PROG_FLDS matches with A_PROG_LOOKUP set up by the UI
        N_PROG_LOOKUP = A_PROG_LOOKUP
        IF (N_PROG_FLDS /= N_PROG_LOOKUP) THEN
          WRITE (6,*) ' '
          WRITE (6,*) ' Mismatch in no of prognostic fields.'
          WRITE (6,*) ' No of prog fields in Atmos dump ',N_PROG_FLDS
          WRITE (6,*) ' No of prog fields expected      ',N_PROG_LOOKUP
          WRITE (6,*) ' '
          WRITE (6,*) ' Run RECONFIGURATION to get correct no of',      &
                      ' prognostic fields in atmos dump'
          WRITE (6,*) ' or'
          WRITE (6,*) ' Check/Reset experiment in User Interface'
          WRITE (6,*) ' '
          CMESSAGE = 'INITDUMP: Wrong no of atmos prognostic fields'
          ICODE = 102
          GO TO 9999  !  Return
        END IF

! Initialise D1 to prevent uninitialised data in unused rows of U fields
! *DIR$ CACHE_BYPASS D1
          DO I = 1,LEN_TOT
            D1(I)=0.0
          END DO
!       Determine no of fields to be read in
        IF (L_A_PROG_ONLY) THEN

!         Prognostic fields only
          LEN2_LOOKUP = N_PROG_FLDS
          LEN_PROG = A_PROG_LEN
          TOT_LEN_DATA = A_LEN_DATA

          LEN_DATA = LEN_PROG

          WRITE (6,*) ' '
          WRITE (6,*) ' Read in ',N_PROG_FLDS,' prognostic fields.'


!      INITIALISE DIAGNOSTIC AREA OF D1 TO RMDI
       DO I = LEN_DATA+1, TOT_LEN_DATA
         D1(I)=RMDI
       END DO

        ELSE

!         All fields.
          LEN2_LOOKUP = A_LEN2_LOOKUP
          LEN_DATA    = A_LEN_DATA

          WRITE (6,*) ' '
          WRITE (6,*) ' Read in all ',LEN2_LOOKUP,' fields.'

        END IF

! Ensure that domain decomposition is consistent with submodel

      CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)


!L 1.2 Call READDUMP to read atmosphere dump.

      D1_ADDR_SUBMODEL_ID = SUBMODEL_FOR_SM(atmos_sm)

! DEPENDS ON: um_readdump
      CALL UM_READDUMP(NFTIN, A_FIXHD, LEN_FIXHD,                       &
          A_INTHD, A_LEN_INTHD,                                         &
          A_REALHD, A_LEN_REALHD,                                       &
          A_LEVDEPC, A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,                    &
          A_ROWDEPC, A_LEN1_ROWDEPC, A_LEN2_ROWDEPC,                    &
          A_COLDEPC, A_LEN1_COLDEPC, A_LEN2_COLDEPC,                    &
          A_FLDDEPC, A_LEN1_FLDDEPC, A_LEN2_FLDDEPC,                    &
          A_EXTCNST, A_LEN_EXTCNST,                                     &
          A_DUMPHIST, LEN_DUMPHIST,                                     &
          A_CFI1, A_LEN_CFI1,                                           &
          A_CFI2, A_LEN_CFI2,                                           &
          A_CFI3, A_LEN_CFI3,                                           &
          A_LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                             &
          A_MPP_LOOKUP,MPP_LEN1_LOOKUP,                                 &
               atmos_sm,                                                &
               NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                          &
               D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                        &
               LEN_DATA,D1,                                             &
          .TRUE.)

      CALL MODEL_FILE_CLOSE(NFTIN,FT_ENVIRON(NFTIN),                          &
                      LEN_FT_ENVIR(NFTIN),0,0,ICODE)

        IF (ICODE  >   0) THEN
          IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
          RETURN
        END IF

! Check validity of integer header data and print out information.
! Pass through the global numbers so the validity check works
! glsize(1) is the global ROW_LENGTH
! glsize(2) is the global ROWS
! DEPENDS ON: pr_inhda
        CALL PR_INHDA (A_INTHD, A_LEN_INTHD, glsize(1,fld_type_p),      &
           glsize(2,fld_type_p),MODEL_LEVELS,WET_LEVELS,TR_LEVELS,      &
           ST_LEVELS, SM_LEVELS, BL_LEVELS,                             &
       TR_VARS, ICODE, CMESSAGE)

        IF (ICODE  >   0) THEN
          IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
          RETURN
        END IF

! Check validity of real header data and print out information.
! DEPENDS ON: pr_rehda
        CALL PR_REHDA (A_REALHD, A_LEN_REALHD)

        IF ( l_endgame ) THEN
! DEPENDS ON: check_idealise_4A
        CALL CHECK_IDEALISE_4A                                          &
        ( A_LEN_INTHD, A_LEN_REALHD, A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,    &
          A_INTHD, A_REALHD, A_LEVDEPC )
        ELSE
! DEPENDS ON: check_idealise
        CALL CHECK_IDEALISE                                             &
        ( A_LEN_INTHD, A_LEN_REALHD, A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,    &
          A_INTHD, A_REALHD, A_LEVDEPC )
        END IF

        IF (L_A_PROG_ONLY) THEN

!         Need to pass field address and length info to ADDRESS_CHECK
          DO I=1,LEN2_LOOKUP
            A_MPP_ADDR(I) = A_MPP_LOOKUP(P_NADDR,I)
            A_MPP_LEN(I)  = A_MPP_LOOKUP(P_LBLREC,I)
          END DO

! DEPENDS ON: address_check
          CALL ADDRESS_CHECK (A_LOOKUP,A_MPP_ADDR,                      &
            A_MPP_LEN,LEN1_LOOKUP,LEN2_LOOKUP,                          &
                              SI,NITEMS,NSECTS,LEN_DATA,                &
                              ICODE,CMESSAGE)

          IF (ICODE  >   0) THEN
            IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
            RETURN
          END IF

        END IF

! ------------------------------------------------------------
!       Check that packing codes in lookup table is consistent
!       with packing requested in dump_packim.
! --------------------------------------------------------------

! DEPENDS ON: check_dump_packing
        CALL CHECK_DUMP_PACKING (                                       &
             A_FIXHD, LEN_FIXHD,                                        &
             A_LOOKUP, LEN1_LOOKUP, LEN2_LOOKUP,                        &
             dump_packim(sm_ident), ATMOS_IM )

!       Reset A_FIXHD to correspond to Output Dump
        A_FIXHD(152) = A_LEN2_LOOKUP
        A_FIXHD(160) = A_FIXHD(150) + LEN1_LOOKUP*A_LEN2_LOOKUP
        A_FIXHD(161) = global_A_LEN_DATA

!L 1.3 Call SET_ATM_POINTERS to set integer pointers to data in
!L     atmosphere dump and secondary storage area in D1 array.
! DEPENDS ON: set_atm_pointers
      CALL SET_ATM_POINTERS(                                            &
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

      IF (ICODE  >   0) THEN
        IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
        RETURN
      END IF

! Initialise Jules
      CALL init_from_jules_namelists( land_field, ntiles, sm_levels, &
                                      nice ,nice_use )
!----------------------------------------------------------------------
! Set bare soil roughness for use in 1 tile dust scheme
!----------------------------------------------------------------------
! fix to enable CRUNs for two bin dust with 1 surface tile
! this mirrors code in sparm. Currently this is a simple real number
! in the future this may become an array and thus will need to be 
! updated in line with the sparm code.

z0_soil=z0_nvg(soil-npft)

! SETCONA is now called from INITIAL
! Set ELF flag
        ELF=(A_FIXHD(4) == 3.OR.A_FIXHD(4) == 103)

      END IF

 9999 CONTINUE

      IF (lhook) CALL dr_hook('INITDUMP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INITDUMP
