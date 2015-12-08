! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INBOUNDA
!LL
!LL Purpose : Takes as input,the code defining whether updates of
!LL  boundary data are required.  The physical files required are
!LL  identified, and the headers lookup tables are read into common
!LL  blocks.  Reads the update intervals from the boundary datasets.
!LL  Where the update interval is in months or years, the check will be
!LL  made daily.  
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL Logical components covered : C720
!LL
!LL System task : C7
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LLEND
!
!*L  Arguments
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input

      SUBROUTINE INBOUNDA(                                              &
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
! ARGBND Control data calculated from NAMELIST-
        NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
        FIXHD_BOUNDA,INTHD_BOUNDA,LOOKUP_BOUNDA,LOOKUP_COMP_BOUNDA,     &
        REALHD_BOUNDA,                                                  &
     &   A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA,                             &
     &   A_LEN1_ROWDEPCDA,A_LEN2_ROWDEPCDA,                             &
     &   A_LEN1_COLDEPCDA,A_LEN2_COLDEPCDA)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE atm_fields_bounds_mod
      USE check_iostat_mod

      USE Submodel_Mod
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      USE lookup_addresses
      USE dust_parameters_mod, ONLY:             l_dust,                &
           l_dust_div1,       l_dust_div2,       l_dust_div3,           &
           l_dust_div4,       l_dust_div5,       l_dust_div6
      USE um_input_control_mod,  ONLY:                                  &
           l_dust_div1_lbc,                     l_dust_div2_lbc,        &
                             l_dust_div3_lbc,                           &
           l_dust_div4_lbc,                     l_dust_div5_lbc,        &
                             l_dust_div6_lbc,   l_so2,                  &
           l_so2_lbc,        l_dms,             l_dms_lbc,              &
           l_so4_aitken,     l_so4_aitken_lbc,  l_so4_accu,             &
           l_so4_accu_lbc,   l_so4_diss,        l_so4_diss_lbc,         &
           l_nh3,            l_nh3_lbc,         l_soot_new,             &
           l_soot_new_lbc,   l_soot_agd,        l_soot_agd_lbc,         &
           l_soot_cld,       l_soot_cld_lbc,    l_bmass_new,            &
           l_bmass_new_lbc,  l_bmass_agd,       l_bmass_agd_lbc,        &
           l_bmass_cld,      l_bmass_cld_lbc,   l_ocff_new,             &
           l_ocff_new_lbc,   l_ocff_agd,        l_ocff_agd_lbc,         &
           l_ocff_cld,       l_ocff_cld_lbc,    l_nitr_acc,             &
           l_nitr_acc_lbc,   l_nitr_diss,       l_nitr_diss_lbc,        &
           l_soot,           l_ocff,            lcal360,                &
           model_domain
      USE model_file, ONLY: model_file_open

      USE domain_params, ONLY: mt_bi_cyclic_lam, mt_cyclic_lam
      USE lbc_read_data_mod, ONLY: albc_num, albc_swapstep,             &
                                   albc2_starttime_steps,               &
                                   current_lbc_step, rimweightsa,       &
                                   bouncnst
      USE mphys_inputs_mod, ONLY:                                       &
           l_mcr_qcf2,       l_mcr_qcf2_lbc,    l_mcr_qgraup,           &
           l_mcr_qgraup_lbc, l_mcr_qrain_lbc,   l_mcr_qrain

      USE cloud_inputs_mod, ONLY: l_pc2, l_pc2_lbc
      
      USE murk_inputs_mod,  ONLY: l_murk, l_murk_lbc
      
      USE nlstcall_mod, ONLY : Num_ALBCs

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

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
! TYPBND - needs TYPSIZE included first
      ! Headers from atmosphere boundary data sets
      ! Second index of header arrays = 1 Lateral boundary data
      ! = 2 Lower boundary data
      INTEGER :: FIXHD_BOUNDA(LEN_FIXHD,1)      ! Fixed header
      INTEGER :: INTHD_BOUNDA(A_LEN_INTHD,1)    ! Integer header

        ! Lookups for the first set of LBCs in the LBC file
      INTEGER :: LOOKUP_BOUNDA(LEN1_LOOKUP,RIM_LOOKUPSA)

      ! Varying items from the LOOKUP table for the entire LBC file
      INTEGER :: LOOKUP_COMP_BOUNDA(LEN1_LBC_COMP_LOOKUP,BOUND_LOOKUPSA)

      REAL ::  REALHD_BOUNDA(A_LEN_REALHD,1)   ! Real header

      !  Control data calculated from namelist
      INTEGER :: RIM_STEPSA      ! Set by IN_BOUND
      INTEGER :: NBOUND_LOOKUP(2)
      COMMON/CBND/                                                      &
       RIM_STEPSA
! TYPBND end

      INTEGER :: A_LEN1_LEVDEPCDA   ! IN : copy of A_LEN1_LEVDEPC
      INTEGER :: A_LEN2_LEVDEPCDA   ! IN : copy of A_LEN2_LEVDEPC
      INTEGER :: A_LEN1_ROWDEPCDA   ! IN : copy of A_LEN1_ROWDEPC
      INTEGER :: A_LEN2_ROWDEPCDA   ! IN : copy of A_LEN2_ROWDEPC
      INTEGER :: A_LEN1_COLDEPCDA   ! IN : copy of A_LEN1_COLDEPC
      INTEGER :: A_LEN2_COLDEPCDA   ! IN : copy of A_LEN2_COLDEPC
      
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

!    Local variables

      INTEGER                                                           &
     &        I,                                                        &
     &        J,                                                        &
     &        J1,                                                       &
     &        lbc_num,                                                  &
     &        START_BLOCK,                                              &
     &        NFTIN,                                                    &
     &        NFTENV,                                                   &
                                   ! Index in FT_ENV arrays
     &        im_index,                                                 &
                                   ! Internal model index
     &        ELAPSED_DAYS,                                             &
                                   ! Days since basis time
     &        ELAPSED_SECS,                                             &
                                   ! Secs since basis time
     &        CURRENT_TIME_DAYS,                                        &
                                   ! No. of days to current time
     &        CURRENT_TIME_SECS,                                        &
                                   ! No. of secs-in-day to current time
     &        DAYS_TO_DATA_START,                                       &
                                   ! Days  to start of boundary data
     &        SECS_TO_DATA_START,                                       &
                                   ! Secs  to start of boundary data
     &        DAYS_TO_DATA_END,                                         &
                                   ! Days  to end   of boundary data
     &        SECS_TO_DATA_END,                                         &
                                   ! Secs  to end   of boundary data
     &        STEPS_TO_DATA_START,                                      &
                                   ! Steps to start of boundary data
     &        STEPS_TO_DATA_END,                                        &
                                   ! Steps to end   of boundary data
     &        RIM_STEPSA_OLD,                                           &
                                   ! Data interval for last bndy file
     &        STEPS_TO_BDI_START,                                       &
                                   ! Steps to start/end of
     &        STEPS_TO_BDI_END,                                         &
                                   ! current boundary data interval
     &        BASIS_TO_DATA_START_STEPS,                                &
                                         ! Steps from basis time to
                                         ! start of boundary data
     &        ITEM_BOUNDA(RIM_LOOKUPSA)  ! Boundary updatable item list

      LOGICAL :: THIS_ALBC_FOR_BDI_END ! True if same boundary file to
                                       ! be used for end of current
                                       ! boundary data interval

      REAL                                                                 &
     &      A_LEVDEPC_BO(A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA),               &
     &      A_ROWDEPC_BO(MAX(A_LEN1_ROWDEPCDA,1),MAX(A_LEN2_ROWDEPCDA,1)), &
     &      A_COLDEPC_BO(MAX(A_LEN1_COLDEPCDA,1),MAX(A_LEN2_COLDEPCDA,1))
      
      INTEGER FULL_LOOKUP_BOUNDA(LEN1_LOOKUP,BOUND_LOOKUPSA)

      INTEGER, Parameter :: DUMMY =1
      INTEGER, Parameter :: sect36 = 36000
      INTEGER, Parameter :: sect37 = 37000

      INTEGER             :: ErrorStatus      ! Return code
      CHARACTER (Len=256) :: CMESSAGE         ! Error message
      CHARACTER (Len=*), Parameter :: RoutineName = 'inbounda'

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

      LOGICAL, SAVE :: L_FirstCall = .TRUE.

      INTEGER :: k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('INBOUNDA',zhook_in,zhook_handle)

!L      Internal Structure

      ErrorStatus=0
      CMESSAGE=' '

      IF (L_FirstCall) THEN

!  Rewind file and read namelist for boundary updating constants

        REWIND 5
        READ (UNIT=5, NML=BOUNCNST, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist BOUNCNST")
        

!L  1.0 Initialise variables in COMMON/CBND/, lest undefined for some
!L      choices of boundary updating, eg. -DEF,GLOBAL but not DEF,FLOOR,
!L      as used in section 2 to set BOUNDARY_STEPSim for atmos
        RIM_STEPSA = 0
!  Initialise BNDARY_OFFSETim & BOUNDARY_STEPSim in CTIME comdeck
        DO I=1,INTERNAL_ID_MAX
          BNDARY_OFFSETim(I) = 0
          BOUNDARY_STEPSim(I) = 0
        END DO

      END IF ! (L_FirstCall)

!L  1.1 Update interval for lateral boundaries for atmosphere
!L      Read headers and test whether boundary updating required

      IF (l_fixed_lbcs .OR. l_force_lbc .OR.                          & 
          model_domain == mt_bi_cyclic_lam .OR.                       &
          model_domain == mt_cyclic_lam) THEN
        RIM_STEPSA=0
      ELSE

!L      Open input boundary file and read headers

        NBOUND_LOOKUP(1)=1

        ! Use same unit number for first and second boundary files
        NFTIN = 125

        IF (ALBC_num == 2) THEN
          NFTENV = 126 ! Env Var for 2nd boundary file
        ELSE
          NFTENV = 125 ! Env Var for 1st boundary file
        END IF


        CALL MODEL_FILE_OPEN(NFTIN,FT_ENVIRON(NFTENV),                        &
            LEN_FT_ENVIR(NFTENV),0,0,ErrorStatus)
        IF(ErrorStatus /= 0) THEN
          Write (Cmessage,*)                                            &
     &    'Failure opening boundary file : No ', ALBC_num

          Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF

!       Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
        CALL READ_FLH(NFTIN,FIXHD_BOUNDA(1,1),                          &
     &                       LEN_FIXHD,ErrorStatus,CMESSAGE)
        IF (ErrorStatus >  0) THEN
          WRITE (Cmessage,*)                                            &
     &          'INBOUNDA : Error in READ_FLH for BOUNDA(1,1)'

           Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF

!       Check for negative dimensions
        IF (FIXHD_BOUNDA(101,1) <= 0) FIXHD_BOUNDA(101,1)=1
        IF (FIXHD_BOUNDA(106,1) <= 0) FIXHD_BOUNDA(106,1)=1
        IF (FIXHD_BOUNDA(111,1) <= 0) FIXHD_BOUNDA(111,1)=1
        IF (FIXHD_BOUNDA(112,1) <= 0) FIXHD_BOUNDA(112,1)=1
        IF (FIXHD_BOUNDA(116,1) <= 0) FIXHD_BOUNDA(116,1)=1
        IF (FIXHD_BOUNDA(117,1) <= 0) FIXHD_BOUNDA(117,1)=1
        IF (FIXHD_BOUNDA(121,1) <= 0) FIXHD_BOUNDA(121,1)=1
        IF (FIXHD_BOUNDA(122,1) <= 0) FIXHD_BOUNDA(122,1)=1
        IF (FIXHD_BOUNDA(151,1) <= 0) FIXHD_BOUNDA(151,1)=1
        IF (FIXHD_BOUNDA(152,1) <= 0) FIXHD_BOUNDA(152,1)=1
        IF (FIXHD_BOUNDA(161,1) <= 0) FIXHD_BOUNDA(161,1)=1

!       Check if sufficient space allocated for LOOKUP table
        IF (FIXHD_BOUNDA(152,1) >  BOUND_LOOKUPSA) THEN
      write(6,*)' INBOUNDA; not enough space for LBC lookup headers.'
      write(6,*)'           try increasing value specified in umui'
      write(6,*)'           window atmos_Infile_Options_Headers'
          write(CMESSAGE,*)                                             &
     &        'INBOUNDA: Insufficient space for Lookup Table'
          Errorstatus = 2

          Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF


        CALL SETPOS (NFTIN,0,ErrorStatus)
        IF (ErrorStatus >  0) THEN
          WRITE (6,*) 'INBOUNDA: Problem with SETPOS for BOUNDA(1,1)'
          WRITE (6,*) 'ErrorStatus ',ErrorStatus,' NFTIN ',NFTIN
          Write(Cmessage,*) 'Problem with Setpos for Bounda(1,1)'

          Call Ereport(RoutineName,Errorstatus,Cmessage)
        ENDIF

! DEPENDS ON: readhead
        CALL READHEAD(NFTIN,                                            &
     &                FIXHD_BOUNDA(1,1),LEN_FIXHD,                      &
     &                INTHD_BOUNDA(1,1),FIXHD_BOUNDA(101,1),            &
     &                REALHD_BOUNDA(1,1),FIXHD_BOUNDA(106,1),           &
     &                A_LEVDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(111,1),FIXHD_BOUNDA(112,1),          &
     &                A_ROWDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(116,1),FIXHD_BOUNDA(117,1),          &
     &                A_COLDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(121,1),FIXHD_BOUNDA(122,1),          &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                FULL_LOOKUP_BOUNDA,                               &
     &                FIXHD_BOUNDA(151,1),FIXHD_BOUNDA(152,1),          &
     &                FIXHD_BOUNDA(161,1),                              &
     &                START_BLOCK,ErrorStatus,CMESSAGE)

        IF (ErrorStatus >  0) THEN
          WRITE (6,*) 'INBOUNDA: Problem with READHEAD for BOUNDA(1,1)'
          WRITE (6,*) 'ErrorStatus ',ErrorStatus,' CMESSAGE ',CMESSAGE
          Write(cmessage,*) 'Problem with READHEAD for BOUNDA(1,1)'

          Call Ereport(RoutineName,errorStatus,cmessage)
        ENDIF

! Copy the first set of headers into LOOKUP_BOUNDA
        DO i=1,RIM_LOOKUPSA
          DO j=1,LEN1_LOOKUP
            LOOKUP_BOUNDA(j,i)=FULL_LOOKUP_BOUNDA(j,i)
          ENDDO ! j
        ENDDO ! i

! Copy all the varying items from the LOOKUP into COMP_LOOKUP_BOUNDA
        DO i=1,BOUND_LOOKUPSA

          LOOKUP_COMP_BOUNDA(LBCC_LBYR,i)=FULL_LOOKUP_BOUNDA(LBYR,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBMON,i)=FULL_LOOKUP_BOUNDA(LBMON,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBDAT,i)=FULL_LOOKUP_BOUNDA(LBDAT,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBHR,i)=FULL_LOOKUP_BOUNDA(LBHR,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBMIN,i)=FULL_LOOKUP_BOUNDA(LBMIN,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBSEC,i)=FULL_LOOKUP_BOUNDA(LBSEC,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBEGIN,i)=FULL_LOOKUP_BOUNDA(LBEGIN,i)
          LOOKUP_COMP_BOUNDA(LBCC_NADDR,i)=FULL_LOOKUP_BOUNDA(NADDR,i)

        ENDDO ! i

! Check validity of headers

! Integer headers
        IF (FIXHD_BOUNDA(100,1)  >   0) THEN

          IF (INTHD_BOUNDA(6,1) /= glsize(1,fld_type_p)) THEN
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'ROW_LENGTH from INTHD: ',INTHD_BOUNDA(6,1)
            WRITE(6,*) 'Model ROW_LENGTH: ',glsize(1,fld_type_p)

            ErrorStatus=3
            Write(CMESSAGE,*)'Integer header (ROW_LENGTH) mismatch'

            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          ENDIF

          IF (INTHD_BOUNDA(7,1) /= glsize(2,fld_type_p)) THEN
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'Number of rows from INTHD: ',INTHD_BOUNDA(7,1)
            WRITE(6,*) 'Model number of rows: ',glsize(2,fld_type_p)

            ErrorStatus=4
            Write(CMESSAGE,*)'Integer header (N_ROWS) mismatch'

            Call Ereport(Routinename,Errorstatus,cmessage)
          ENDIF

          If ( INTHD_BOUNDA(17,1) /= A_INTHD(17) ) Then
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'LBC   : Height Generator Method : ',            &
     &                  INTHD_BOUNDA(17,1)
            WRITE(6,*) 'Model : Height Generator Method : ',            &
     &                  A_INTHD(17)

            ErrorStatus=5
            Write (CMESSAGE,*) 'INBOUNDA : Mis-match in height ',       &
     &                         'generator method.'

            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          End If

          If ( INTHD_BOUNDA(24,1) /= A_INTHD(24) ) Then
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'LBC  : First rho level with constant height: ', &
     &                  INTHD_BOUNDA(24,1)
            WRITE(6,*) 'Model: First rho level with constant height: ', &
     &                  A_INTHD(24)

            ErrorStatus=6
            Write (CMESSAGE,*) 'INBOUNDA : Mis-match in height ',       &
     &                         'generator method.'

            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          End if

        ENDIF ! IF (FIXHD_BOUNDA(100,1)  >   0)

! Real constants

        IF (FIXHD_BOUNDA(105,1)  >   0) THEN

          ! Check real headers in LBC file against those in dump to
          ! 32 bit accuracy
          DO j=1,6
            IF (LNER(REALHD_BOUNDA(j,1),A_REALHD(j))) THEN
              WRITE(6,*) 'LBC Real Header mismatch at position ',j,':'
              WRITE(6,*) 'Value from LBC file is ',REALHD_BOUNDA(j,1)
              WRITE(6,*) 'Value from model dump is ',A_REALHD(j)

              ErrorStatus=7
              Write(CMESSAGE,*) 'INBOUNDA : Real header mismatch'

              Call Ereport(RoutineName,ErrorStatus,Cmessage)
            ENDIF
          ENDDO ! j

          IF (LNER(REALHD_BOUNDA(16,1),A_REALHD(16))) THEN
            WRITE(6,*) 'LBC Real Header mismatch at position 16 :'
            WRITE(6,*) 'LBC file : Height at top of model ',            &
     &      REALHD_BOUNDA(16,1)
            WRITE(6,*) 'Model    : Height at top of model ',            &
     &      A_REALHD(16)

            ErrorStatus=8
            Write(CMESSAGE,*) 'INBOUNDA : Model Height mismatch'

            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          ENDIF

        ENDIF ! IF (FIXHD_BOUNDA(105,1)  >   0)

! Level dependent constants

        If (FIXHD_BOUNDA(110,1) > 0) Then

! ---------------------------------
! Check eta values for theta levels
! ---------------------------------

          DO j=1, model_levels+1
            If ( LNER( A_LEVDEPC_BO(j,1),                               &
     &                 A_LEVDEPC(jetatheta+j-1) ) ) Then

              WRITE(6,*) 'LBC Level Dependent Constants Mismatch'
              WRITE(6,*) 'Eta values for theta levels mismatch for ',   &
     &                   'level ',J
              WRITE(6,*) 'Value from LBC file   : ',A_LEVDEPC_BO(j,1)
              WRITE(6,*) 'Value from model dump : ',                    &
     &                    A_LEVDEPC(jetatheta+j-1)
              ErrorStatus=9
              Write(CMESSAGE,*)                                         &
     &              'INBOUNDA : Level dependent constant mismatch'

              Call Ereport(RoutineName,ErrorStatus,Cmessage)

            End If
          End Do ! j

! -------------------------------
! Check eta values for rho levels
! -------------------------------
          DO j=1, model_levels
            If ( LNER( A_LEVDEPC_BO(j,2),                               &
     &                 A_LEVDEPC(jetarho+j-1) ) ) Then

              WRITE(6,*) 'LBC Level Dependent Constants Mismatch'
              WRITE(6,*) 'Eta values for rho levels mismatch for ',     &
     &                   'level ',J
              WRITE(6,*) 'Value from LBC file   : ',A_LEVDEPC_BO(j,2)
              WRITE(6,*) 'Value from model dump : ',                    &
     &                    A_LEVDEPC(jetarho+j-1)
              ErrorStatus=10
              Write(CMESSAGE,*)                                         &
     &            'INBOUNDA : Level dependent constant mismatch'

              Call Ereport(RoutineName,ErrorStatus,Cmessage)

            End If
          End Do ! j

        End If ! If (FIXHD_BOUNDA(110,1) > 0)

! Check that for variable resolution input, the LBC file has the same
! grid as the input dump
        IF( (a_fixhd(115) > 0 .AND. a_fixhd(120) > 0) .OR.              &
            (fixhd_bounda(115,1) > 0 .AND. fixhd_bounda(120,1) > 0) ) THEN
          ! Input dump is variable resolution
          IF (fixhd_bounda(116,1) /= a_fixhd(116) .OR.                  &
              fixhd_bounda(117,1) /= a_fixhd(117) .OR.                  &
              fixhd_bounda(121,1) /= a_fixhd(121) .OR.                  &
              fixhd_bounda(122,1) /= a_fixhd(122) ) THEN
              
            WRITE(cmessage,'(A)') "Variable resolution grid " //        &
              "dimensions in LBC file do not match input dump"
            errorstatus = -10
            CALL ereport(RoutineName, errorstatus, cmessage)
          END IF

          DO j = 1, a_fixhd(117)
            DO i = 1, a_fixhd(116)
              k = (j-1) * a_fixhd(116) + i
              IF (lner(a_rowdepc(k), a_rowdepc_bo(i,j))) THEN
                WRITE(cmessage,'(A)') "Row dependent constants " //     &
              "in LBC file do not match input dump"
                errorstatus = -11
                CALL ereport(RoutineName, errorstatus, cmessage)
              END IF
            END DO
          END DO

          DO j = 1, a_fixhd(122)
            DO i = 1, a_fixhd(121)
              k = (j-1) * a_fixhd(121) + i
              IF (lner(a_coldepc(k), a_coldepc_bo(i,j))) THEN
                WRITE(cmessage,'(A)') "Column dependent constants " //  &
              "in LBC file do not match input dump"
                errorstatus = -12
                CALL ereport(RoutineName, errorstatus, cmessage)
              END IF
            END DO
          END DO
        END IF ! If variable resolution input dump




!L      Set update interval
!       If update interval includes months or years, a 360 day
!       calender assumed.

        ! Save previous value of RIM_STEPSA:
        IF (.NOT.L_FirstCall) RIM_STEPSA_OLD = RIM_STEPSA

        RIM_STEPSA=((FIXHD_BOUNDA(35,1)*8640+FIXHD_BOUNDA(36,1)*720     &
     &   +FIXHD_BOUNDA(37,1)*24+FIXHD_BOUNDA(38,1))*3600                &
     &   +FIXHD_BOUNDA(39,1)*60+FIXHD_BOUNDA(40,1))                     &
     &   *STEPS_PER_PERIODim(a_im)/SECS_PER_PERIODim(a_im)

        ! Check that RIM_STEPSA has not changed:
        IF (.NOT.L_FirstCall) THEN
          IF (RIM_STEPSA /= RIM_STEPSA_OLD) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &      'Boundary updating period (RIM_STEPSA) has changed from ',  &
     &      RIM_STEPSA_OLD, ' to ', RIM_STEPSA, '. Not allowed!'

            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF
        END IF

        ! Initialise Current_LBC_Step

        If (L_FirstCall) Then
          Current_LBC_Step = 1 + STEPim(atmos_im)
          If (PrintStatus >= PrStatus_Normal) Then
            write (6,*) ' INBOUNDA : Timestep ',STEPim(a_im),           &
     &                  ' Current_LBC_Step ',Current_LBC_Step
          End If
        End If

        IF (Num_ALBCs == 2 .AND. L_FirstCall) THEN

          ! Calculate step on which to swap boundary files:
          ALBC_SwapStep = ALBC2_StartTime_steps - RIM_STEPSA

          ! Check that ALBC_SwapStep is non-negative:
          IF (ALBC_SwapStep < 0) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Step on which to swap boundary files (ALBC_SwapStep = ', &
     &        ALBC_SwapStep, ') must be non-negative'

            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

          ! Check that ALBC_SwapStep is a multiple of RIM_STEPSA:
          IF (MOD(ALBC_SwapStep, RIM_STEPSA) /= 0) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Step on which to swap boundary files (ALBC_SwapStep = ', &
     &        ALBC_SwapStep, ') must be a a multiple of the boundary '//&
     &        'updating period (RIM_STEPSA = ', RIM_STEPSA, ')'
            ErrorStatus=101

            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

        END IF

      END IF

!L  2   Set interval for setting any boundary field

      BOUNDARY_STEPSim(a_im) = RIM_STEPSA

!L  3    Check LOOKUP Table

      J1=0
      IF (.NOT. (model_domain == mt_bi_cyclic_lam .OR.                  &
          model_domain == mt_cyclic_lam    .OR.                         &
          l_force_lbc .OR. l_fixed_lbcs) ) THEN
        J1=FIXHD_BOUNDA(152,1)

        IF(FIXHD_BOUNDA(150,1) >  0) THEN

! Set up list of variables expected to be boundary updated.
          ITEM_BOUNDA(1)  = 31001 ! Orography
          ITEM_BOUNDA(2)  = 31002 ! U
          ITEM_BOUNDA(3)  = 31003 ! V
          ITEM_BOUNDA(4)  = 31004 ! W
          ITEM_BOUNDA(5)  = 31005 ! Density
          ITEM_BOUNDA(6)  = 31006 ! Potential temperature
          ITEM_BOUNDA(7)  = 31007 ! Specific humidity
          ITEM_BOUNDA(8)  = 31008 ! QCL
          ITEM_BOUNDA(9)  = 31009 ! QCF
          ITEM_BOUNDA(10) = 31010 ! Exner
          ITEM_BOUNDA(11) = 31011 ! U_Adv
          ITEM_BOUNDA(12) = 31012 ! V_Adv
          ITEM_BOUNDA(13) = 31013 ! W Adv

! Include additional microphysics lbcs if active

          lbc_num = 13
          IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31014
          ENDIF
          IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31015
          ENDIF
          IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31016
          ENDIF

          ! Setup for additional cloud fraction lbcs if active
          If (L_pc2_lbc) Then
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31017 ! cf_bulk
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31018 ! cf_liquid
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31019 ! cf_frozen
          EndIf

          ! Include murk aerosol lbcs if in input lbc file
          If (L_murk_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31020 ! murk aerosol
          EndIf

          ! Include dust lbcs if in input lbc file
          If (L_dust_div1_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31023 ! dust_div1
          EndIf
          ! Include dust lbcs if in input lbc file
          If (L_dust_div2_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31024 ! dust_div2
          EndIf
          ! Include dust lbcs if in input lbc file
          If (L_dust_div3_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31025 ! dust_div3
          EndIf
          ! Include dust lbcs if in input lbc file
          If (L_dust_div4_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31026 ! dust_div4
          EndIf
          ! Include dust lbcs if in input lbc file
          If (L_dust_div5_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31027 ! dust_div5
          EndIf
          ! Include dust lbcs if in input lbc file
          If (L_dust_div6_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31028 ! dust_div6
          EndIf

          ! Include so2 lbcs if in input lbc file
          If (L_so2_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31029 ! so2
          EndIf

          ! Include dms lbcs if in input lbc file
          If (L_dms_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31030 ! dms
          EndIf

          ! Include so4_aitken lbcs if in input lbc file
          If (L_so4_aitken_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31031 ! so4_aitken
          EndIf
          ! Include so4_accu lbcs if in input lbc file
          If (L_so4_accu_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31032 ! so4_accu
          EndIf
          ! Include so4_diss lbcs if in input lbc file
          If (L_so4_diss_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31033 ! so4_diss
          EndIf
  
          ! Include nh3 lbcs if in input lbc file
          If (L_nh3_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31035 ! nh3
          EndIf

          ! Include soot_new lbcs if in input lbc file
          If (L_soot_new_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31036 ! soot_new
          EndIf
          ! Include soot_agd lbcs if in input lbc file
          If (L_soot_agd_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31037 ! soot_agd
          EndIf
          ! Include soot_cld lbcs if in input lbc file
          If (L_soot_cld_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31038 ! soot_cld
          EndIf
 
         ! Include bmass_new lbcs if in input lbc file
          If (L_bmass_new_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31039 ! bmass_new
          EndIf
          ! Include bmass_agd lbcs if in input lbc file
          If (L_bmass_agd_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31040 ! bmass_agd
          EndIf
          ! Include bmass_cld lbcs if in input lbc file
          If (L_bmass_cld_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31041 ! bmass_cld
          EndIf

          ! Include ocff_new lbcs if in input lbc file
          If (L_ocff_new_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31042 ! ocff_new
          EndIf
          ! Include ocff_agd lbcs if in input lbc file
          If (L_ocff_agd_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31043 ! ocff_agd
          EndIf
          ! Include ocff_cld lbcs if in input lbc file
          If (L_ocff_cld_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31044 ! ocff_cld
          EndIf

          ! Include nitr_acc lbcs if in input lbc file
          If (L_nitr_acc_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31045 ! nitr_new
          EndIf
          ! Include nitr_diss lbcs if in input lbc file
          If (L_nitr_diss_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31046 ! nitr_diss
          EndIf

          ! Free tracer LBCs
          If (TR_LBC_VARS > 0) Then                                 
                                                                    
            ! Set ITEM_BOUNDA to STASH item codes read in from the 
            ! file so that the user only needs to specify the total 
            ! number of tracer lbcs

            Do i = 1, TR_LBC_VARS 
              ITEM_BOUNDA(lbc_num+i)=sect36 + A_TR_lbc_stashitem(i)
              Write(6,*) 'INBOUNDA (free):',lbc_num+i,ITEM_BOUNDA(lbc_num+i) 
            End Do                                               
            lbc_num = lbc_num + TR_LBC_VARS
                                                                    
          End If              

!          IF (TR_VARS  >   0) THEN

            ! Find STASH item no. for each tracer in use
!            i=0
!            im_index=internal_model_index(A_IM)
!            DO j = A_TRACER_FIRST,A_TRACER_LAST
!              IF (SI(j,0,im_index) /= 1) THEN  ! tracer is in use
!                i=i+1
!                ITEM_BOUNDA(lbc_num+i)=31100+(j-A_TRACER_FIRST)
!              ENDIF
!            ENDDO ! j

!          ENDIF ! IF (TR_VARS  >   0)

          ! UKCA tracer LBCs
          If (TR_LBC_UKCA > 0) Then                                 
                                                                    
            ! Set ITEM_BOUNDA to STASH item codes read in from the 
            ! file so that the user only needs to specify the total 
            ! number of tracer lbcs
            Do i = 1, TR_LBC_UKCA 
              ITEM_BOUNDA(lbc_num+i)=sect37 + UKCA_TR_lbc_stashitem(i)
              Write(6,*) 'INBOUNDA (ukca):',lbc_num+i,ITEM_BOUNDA(lbc_num+i)   
            End Do                                               
                                                                    
          End If              


! DEPENDS ON: chk_look_bounda
          CALL CHK_LOOK_BOUNDA(                                         &
     &      ITEM_BOUNDA,FULL_LOOKUP_BOUNDA,                             &
! ARGBND Control data calculated from NAMELIST-
        NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
        FIXHD_BOUNDA,INTHD_BOUNDA,LOOKUP_BOUNDA,LOOKUP_COMP_BOUNDA,     &
        REALHD_BOUNDA,                                                  &
     &                         ErrorStatus,CMESSAGE)

          IF (ErrorStatus  /=  0) Then
! Use message returned by chk_look_bounda in call to ereport

            CALL ereport(RoutineName,ErrorStatus,Cmessage)
          End If

!L Find start position in lookup tables

          ! Get days/seconds since basis time:
! DEPENDS ON: time2sec
          CALL TIME2SEC (I_YEAR,             I_MONTH,                   &
     &                   I_DAY,              I_HOUR,                    &
     &                   I_MINUTE,           I_SECOND,                  &
     &                   BASIS_TIME_DAYS,    BASIS_TIME_SECS,           &
     &                   ELAPSED_DAYS,       ELAPSED_SECS,              &
     &                   LCAL360)

          ! Get current model time in days/seconds-in-day:
          CURRENT_TIME_DAYS = BASIS_TIME_DAYS + ELAPSED_DAYS +          &
     &                       (BASIS_TIME_SECS + ELAPSED_SECS)/86400

          CURRENT_TIME_SECS = MOD(BASIS_TIME_SECS+ELAPSED_SECS, 86400)

          ! Get days/seconds to start of boundary data:
! DEPENDS ON: time2sec
          CALL TIME2SEC (FIXHD_BOUNDA(21,1), FIXHD_BOUNDA(22,1),        &
     &                   FIXHD_BOUNDA(23,1), FIXHD_BOUNDA(24,1),        &
     &                   FIXHD_BOUNDA(25,1), FIXHD_BOUNDA(26,1),        &
     &                   CURRENT_TIME_DAYS,  CURRENT_TIME_SECS,         &
     &                   DAYS_TO_DATA_START, SECS_TO_DATA_START,        &
     &                   LCAL360)

          ! Get days/seconds to end of boundary data:
! DEPENDS ON: time2sec
          CALL TIME2SEC (FIXHD_BOUNDA(28,1), FIXHD_BOUNDA(29,1),        &
     &                   FIXHD_BOUNDA(30,1), FIXHD_BOUNDA(31,1),        &
     &                   FIXHD_BOUNDA(32,1), FIXHD_BOUNDA(33,1),        &
     &                   CURRENT_TIME_DAYS,  CURRENT_TIME_SECS,         &
     &                   DAYS_TO_DATA_END,   SECS_TO_DATA_END,          &
     &                   LCAL360)

          ! Get steps to start of boundary data:
! DEPENDS ON: tim2step
          CALL TIM2STEP (DAYS_TO_DATA_START,                            &
     &                   SECS_TO_DATA_START,                            &
     &                   STEPS_PER_PERIODim(a_im),                      &
     &                   SECS_PER_PERIODim(a_im),                       &
     &                   STEPS_TO_DATA_START)

          ! Get steps to end of boundary data:
! DEPENDS ON: tim2step
          CALL TIM2STEP (DAYS_TO_DATA_END,                              &
     &                   SECS_TO_DATA_END,                              &
     &                   STEPS_PER_PERIODim(a_im),                      &
     &                   SECS_PER_PERIODim(a_im),                       &
     &                   STEPS_TO_DATA_END)

          ! Get steps from basis time to start of boundary data:
          BASIS_TO_DATA_START_STEPS = STEPim(a_im) + STEPS_TO_DATA_START

          ! Check that the above is a multiple of RIM_STEPSA only
          ! when two boundary files exist
          IF ((Num_ALBCs == 2) .AND.                                    &
     &        (MOD(BASIS_TO_DATA_START_STEPS, RIM_STEPSA) /= 0)) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Steps from basis time to start of boundary data (',      &
     &        BASIS_TO_DATA_START_STEPS, ') must be a multiple of '//   &
     &        'the boundary updating period (', RIM_STEPSA, ')'

            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

          ! There are two situations we need to cater for here:
          !
          !   1. There is no boundary data in memory.
          !   2. The boundary data valid at the start of the current
          !      boundary data interval is already in memory, but
          !      data for the end of the interval needs to be read in.
          !
          ! The first case applies if and only if we are on the first
          ! call to this routine.
          !
          ! In the first case, we need to make sure that there is
          ! boundary data valid at the start of the current boundary
          ! data interval, and then point to it. We also need to make
          ! sure that data is available for the end of the interval. The
          ! only exception to the latter is the case where the data for
          ! the end of the interval is to come from a different boundary
          ! file, which may be the case if two boundary files are being
          ! used.
          !
          ! In the second case, we need to make sure that there is
          ! boundary data valid at the end of the current boundary data
          ! interval, and then point to it.

          ! Steps to start of current boundary data interval:
          STEPS_TO_BDI_START = MOD(BASIS_TO_DATA_START_STEPS,           &
     &                             RIM_STEPSA)

          ! Steps to end   of current boundary data interval:
          STEPS_TO_BDI_END   = STEPS_TO_BDI_START + RIM_STEPSA

          ! Initialise BNDARY_OFFSET for use in CTIME comdeck
          BNDARY_OFFSETim(a_im) = MOD(-BASIS_TO_DATA_START_STEPS,       &
     &                                 RIM_STEPSA)

          IF (L_FirstCall) THEN ! Case 1

            IF (STEPS_TO_DATA_START <= STEPS_TO_BDI_START .AND.         &
     &          STEPS_TO_DATA_END   >= STEPS_TO_BDI_START) THEN
              NBOUND_LOOKUP(1) = (-STEPS_TO_DATA_START / RIM_STEPSA)    &
     &                         * (RIM_LOOKUPSA-1) + 2
            ELSE IF (STEPS_TO_DATA_START > STEPS_TO_BDI_START) THEN
              ErrorStatus = 101
              WRITE (CMessage,*)                                        &
     &          'Boundary data starts after start of current '//        &
     &          'boundary data interval'

              CALL EReport (RoutineName, ErrorStatus, CMessage)
            ELSE
              ErrorStatus = 102
              WRITE (CMessage,*)                                        &
     &          'Boundary data ends before start of current '//         &
     &          'boundary data interval'

              CALL EReport (RoutineName, ErrorStatus, CMessage)
            END IF

            ! Is the data for the end of the current boundary data
            ! interval to come from the same boundary file?
            THIS_ALBC_FOR_BDI_END =.TRUE.
            IF (NUM_ALBCs == 2) THEN
              IF (STEPim(a_im) <= ALBC_SwapStep) THEN
                THIS_ALBC_FOR_BDI_END =.FALSE.
              END IF
            END IF

            ! If so, check that the data exists:
            IF (THIS_ALBC_FOR_BDI_END) THEN
              IF (STEPS_TO_DATA_START > STEPS_TO_BDI_END) THEN
                ErrorStatus = 101
                WRITE (CMessage,*)                                      &
     &            'Boundary data starts after end of current '//        &
     &            'boundary data interval'

                CALL EReport (RoutineName, ErrorStatus, CMessage)
              ELSE IF (STEPS_TO_DATA_END < STEPS_TO_BDI_END) THEN
                ErrorStatus = 102
                WRITE (CMessage,*)                                      &
     &            'Boundary data ends before end of current '//         &
     &            'boundary data interval'

                CALL EReport (RoutineName, ErrorStatus, CMessage)
              END IF
            END IF

          ELSE ! Not first call (Case 2)

            IF (STEPS_TO_DATA_START <= STEPS_TO_BDI_END .AND.           &
     &          STEPS_TO_DATA_END   >= STEPS_TO_BDI_END) THEN
              NBOUND_LOOKUP(1) = (1 - STEPS_TO_DATA_START / RIM_STEPSA) &
     &                         * (RIM_LOOKUPSA-1) + 2
            ELSE IF (STEPS_TO_DATA_START > STEPS_TO_BDI_END) THEN
              ErrorStatus = 101
              WRITE (CMessage,*)                                        &
     &          'Boundary data starts after end of current '//          &
     &          'boundary data interval'

              CALL EReport (RoutineName, ErrorStatus, CMessage)
            ELSE
              ErrorStatus = 102
              WRITE (CMessage,*)                                        &
     &          'Boundary data ends before end of current '//           &
     &          'boundary data interval'

              CALL EReport (RoutineName, ErrorStatus, CMessage)
            END IF

          END IF ! (L_FirstCall)

        END IF

      END IF ! lateral boundary

      L_FirstCall = .FALSE.

!L  4   End of routine

      IF (lhook) CALL dr_hook('INBOUNDA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INBOUNDA

