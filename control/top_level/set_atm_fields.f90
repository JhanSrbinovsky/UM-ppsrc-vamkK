! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Points atmosphere fields to the appropriate sections of D1
!
! Subroutine Interface: 
SUBROUTINE Set_Atm_Fields ( &
! "argptra.h" contains jpointers
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
! "argsts.h" contains SI (STASH index array); used to check tracers
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      D1, LD1, ID1 )

       USE atm_fields_mod   ! atmosphere fields
       USE field_length_mod ! field_length function 
       USE atm_fields_bounds_mod
       USE ancil_info, only: nsmax

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE rimtypes
USE lbc_mod
USE um_input_control_mod,  ONLY: l_sulpc_online_oxidants
USE rad_input_mod, ONLY: lexpand_ozone, lexpand_tpps_ozone
USE ukca_option_mod, ONLY: l_ukca
USE Submodel_Mod


       USE chsunits_mod, ONLY : nunits

   use cable_data_mod, only : cable_set_atm_fields
IMPLICIT NONE
!
! Description: 
!   Routine to point atmosphere fields to the appropriate sections of D1.
!   After calling this subroutine, the fields can be used directly without
!   referring to D1
!
! Method: 
!   Assuming SET_ATM_POINTERS has been called beforehand, this subroutine
!   points each field to an area of D1 starting at the corresponding
!   "jpointer" (at the first level) and ending at the "jpointer" (at the 
!   last level) plus the size of a single level of that field. 
!
!   Tracers are dealt with differently:   First the number of active tracers
!   is computed so that the correct sections of the corresponding tracer
!   jpointers can be used in pointing the tracer fields to D1.  If no tracers
!   are active then the fields are pointed to a dummy array
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code description: 
!   Language:  Fortran 90.
!   This code is written to UM programming standards UMDP3 vn 8.3.
!

! Subroutine arguments

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
! constants (NUNITS) in "chsunits.h" are needed by "ccontrol.h"
! constants (L_3D_CCA) in "ccontrol.h" are needed to determine field sizes below
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

! constants (A_TRACER_FIRST, etc.) in "ctracera.h" are needed for tracers
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

      REAL,    TARGET, INTENT(IN) :: D1(LEN_TOT)
      LOGICAL, TARGET, INTENT(IN) :: LD1(LEN_TOT)
      INTEGER, TARGET, INTENT(IN) :: ID1(LEN_TOT)

! Local variables

      INTEGER :: nTracer ! loop counter over available tracers
      INTEGER :: nActiveTracers ! number of tracers actually being used

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! End of header

! 1.0 Start of subroutine code; point fields to D1

!    atmospheric primary variables
      IF (lhook) CALL dr_hook('SET_ATM_FIELDS',zhook_in,zhook_handle)

     Exner      => D1(JEXNER_RHO_LEVELS(pdims%k_start) :                       &
                      JEXNER_RHO_LEVELS(pdims%k_start) +                       &
                      field_length(theta_points,single_halo,pdims%k_end+1) -1)

     Exner_Surf => D1(JEXNERSURF : JEXNERSURF +                                &
                      field_length(theta_points,single_halo,1) -1)

     DryRho     => D1(JDRYRHO(pdims%k_start) : JDRYRHO(pdims%k_start) +        &
                      field_length(theta_points,single_halo,                   &
                                   pdims%k_end-pdims%k_start+1) -1)

     Etadot     => D1(JETADOT(wdims%k_start) : JETADOT(wdims%k_start) +        &
                      field_length(theta_points,single_halo,                   &
                                   wdims%k_end-wdims%k_start+1) -1)

     THETAV     => D1(JTHETAV(tdims%k_start) : JTHETAV(tdims%k_start) +        &
                      field_length(theta_points,single_halo,                   &
                                   tdims%k_end-tdims%k_start+1) -1)

     PSI_W_SURF => D1(JPSIWS : JPSIWS + field_length(theta_points,no_halo,1) -1)

     PSI_W_LID  => D1(JPSIWL : JPSIWL + field_length(theta_points,no_halo,1) -1)

     m_v        => D1(JMV(qdims%k_start) : JMV(qdims%k_start) +                &
                      field_length(theta_points,single_halo,                   &
                              qdims%k_end-qdims%k_start+1) -1)

     m_cl       => D1(JMCL(qdims%k_start) : JMCL(qdims%k_start) +              &
                      field_length(theta_points,single_halo,                   &
                                   qdims%k_end-qdims%k_start+1) -1)

     m_cf       => D1(JMCF(qdims%k_start) : JMCF(qdims%k_start) +              &
                      field_length(theta_points,single_halo,                   &
                                   qdims%k_end-qdims%k_start+1) -1)

     m_cf2      => D1(JMCF2(qdims%k_start) : JMCF2(qdims%k_start) +            &
                      field_length(theta_points,single_halo,                   &
                                   qdims%k_end-qdims%k_start+1) -1)

     m_gr       => D1(JMGRAUP(qdims%k_start) : JMGRAUP(qdims%k_start) +        &
                      field_length(theta_points,single_halo,                   &
                              qdims%k_end-qdims%k_start+1) -1)

     m_r        => D1(JMRAIN(qdims%k_start) : JMRAIN(qdims%k_start) +          &
                      field_length(theta_points,single_halo,                   &
                                   qdims%k_end-qdims%k_start+1) -1)
    
     

     U         => D1(JU(udims%k_start) : JU(udims%k_start) + &
      field_length(u_points,single_halo,udims%k_end-udims%k_start+1) -1)

     V         => D1(JV(udims%k_start) : JV(udims%k_start) + &
      field_length(v_points,single_halo,vdims%k_end-vdims%k_start+1) -1)

     THETA     => D1(JTHETA(tdims%k_start) : JTHETA(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)

     Q         => D1(JQ(qdims%k_start) : JQ(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     QCF       => D1(JQCF(qdims%k_start) : JQCF(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     TSTAR     => D1(JTSTAR : JTSTAR+field_length(theta_points,no_halo,1) -1)
     LAND      => D1(JLAND  : JLAND +field_length(theta_points,no_halo,1) -1)
     OROGRAPHY => D1(JOROG  : JOROG +field_length(theta_points,no_halo,1) -1)

     W         => D1(JW(wdims_s%k_start) : JW(wdims_s%k_start) + &
      field_length(theta_points,single_halo,wdims%k_end-wdims%k_start+1) -1)

     RHO       => D1(JRHO(pdims%k_start) : JRHO(pdims%k_start)+ &
      field_length(theta_points,single_halo,pdims%k_end) -1)

     QCL       => D1(JQCL(qdims%k_start) : JQCL(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     QCF2      => D1(JQCF2(qdims%k_start) : JQCF2(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     QRAIN     => D1(JQRAIN(qdims%k_start) : JQRAIN(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     QGRAUP    => D1(JQGRAUP(qdims%k_start) : JQGRAUP(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     EXNER_RHO_LEVELS => D1(JEXNER_RHO_LEVELS(pdims%k_start) : &
                            JEXNER_RHO_LEVELS(pdims%k_start) + &
                      field_length(theta_points,single_halo,pdims%k_end+1) -1)

     E_TRB     => D1(JE_TRB(tdims%k_start) : JE_TRB(tdims%k_start)+ &
      field_length(theta_points,extended_halo,tdims%k_end-tdims%k_start+1) -1)
     TSQ_TRB   => D1(JTSQ_TRB(tdims%k_start) : JTSQ_TRB(tdims%k_start)+ &
      field_length(theta_points,extended_halo,tdims%k_end-tdims%k_start+1) -1)
     QSQ_TRB   => D1(JQSQ_TRB(tdims%k_start) : JQSQ_TRB(tdims%k_start)+ &
      field_length(theta_points,extended_halo,tdims%k_end-tdims%k_start+1) -1)
     COV_TRB   => D1(JCOV_TRB(tdims%k_start) : JCOV_TRB(tdims%k_start)+ &
      field_length(theta_points,extended_halo,tdims%k_end-tdims%k_start+1) -1)
     ZHPAR_SHCU=> D1(JZHPAR_SHCU : JZHPAR_SHCU +                        &
      field_length(theta_points,no_halo,1) -1)

!    Coastal Tiling
     FRAC_LAND  => D1(JFRAC_LAND:JFRAC_LAND  + &
                                         field_length(land_points,no_halo,1) -1)
     TSTAR_LAND => D1(JTSTAR_LAND:JTSTAR_LAND+ &
                                         field_length(theta_points,no_halo,1)-1)
     TSTAR_SEA  => D1(JTSTAR_SEA:JTSTAR_SEA  + &
                                         field_length(theta_points,no_halo,1)-1)
     TSTAR_SICE => D1(JTSTAR_SICE:JTSTAR_SICE+ &
                                         field_length(theta_points,no_halo,1)-1)
     TSTAR_SICE_CAT => D1(JTSTAR_SICE_CAT:JTSTAR_SICE_CAT+ &
                                  field_length(theta_points,no_halo,nice_use)-1)

!    SeaIce and Land albedos
     SICE_ALB => D1(JSICE_ALB : JSICE_ALB + &
                                         field_length(theta_points,no_halo,1)-1)
     LAND_ALB => D1(JLAND_ALB : JLAND_ALB + &
                                         field_length(theta_points,no_halo,1)-1)

!    Large-Scale hydrology
     TI_MEAN   => D1(JTI_MEAN:JTI_MEAN+field_length(land_points,no_halo,1) -1)
     TI_SIG    => D1(JTI_SIG:JTI_SIG  +field_length(land_points,no_halo,1) -1)
     FEXP      => D1(JFEXP:JFEXP      +field_length(land_points,no_halo,1) -1)
     GAMMA_INT => D1(JGAMMA_INT:JGAMMA_INT + &
                                       field_length(land_points,no_halo,1) -1)
     FSFC_SAT  => D1(JFSFC_SAT:JFSFC_SAT   + &
                                       field_length(land_points,no_halo,1) -1)
     F_WETLAND => D1(JF_WETLAND:JF_WETLAND + &
                                       field_length(land_points,no_halo,1) -1)
     WATER_TABLE => D1(JWATER_TABLE:JWATER_TABLE+ &
                                       field_length(land_points,no_halo,1) -1)

     STHZW   => D1(JSTHZW  : JSTHZW  + field_length(land_points,no_halo,1) -1)
     A_FSAT  => D1(JA_FSAT : JA_FSAT + field_length(land_points,no_halo,1) -1)
     C_FSAT  => D1(JC_FSAT : JC_FSAT + field_length(land_points,no_halo,1) -1)
     A_FWET  => D1(JA_FWET : JA_FWET + field_length(land_points,no_halo,1) -1)
     C_FWET  => D1(JC_FWET : JC_FWET + field_length(land_points,no_halo,1) -1)

!    Optional atmospheric primary variables

     ZH        => D1(JZH    : JZH    +field_length(theta_points,no_halo,1) -1)
     ddmfx     => D1(jddmfx : jddmfx +field_length(theta_points,no_halo,1) -1)

     U_ADV     => D1(JU_ADV(udims%k_start) : JU_ADV(udims%k_start)+ &
                          field_length(u_points,extended_halo,udims%k_end) -1)
     V_ADV     => D1(JV_ADV(vdims%k_start) : JV_ADV(vdims%k_start)+ &
                          field_length(v_points,extended_halo,vdims%k_end) -1)

     W_ADV     => D1(JW_ADV(wdims_s%k_start) : JW_ADV(wdims_s%k_start)+ &
                  field_length(theta_points,extended_halo,wdims_s%k_end+1) -1)

     NTML      => D1(JNTML  : JNTML +field_length(theta_points,no_halo,1) -1)
     NBDSC     => D1(JNBDSC : JNBDSC+field_length(theta_points,no_halo,1) -1)
     NTDSC     => D1(JNTDSC : JNTDSC+field_length(theta_points,no_halo,1) -1)
     CUMULUS   => D1(JCUMULUS : JCUMULUS+field_length(theta_points,no_halo,1)-1)

     T1_SD     => D1(JT1_SD : JT1_SD+field_length(theta_points,no_halo,1) -1)
     Q1_SD     => D1(JQ1_SD : JQ1_SD+field_length(theta_points,no_halo,1) -1)

     CF_AREA   => D1(JCF_AREA(1) : JCF_AREA(1)+ &
      field_length(theta_points,no_halo      ,qdims%k_end) -1)

     CF_BULK   => D1(JCF_BULK(qdims%k_start) : JCF_BULK(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     CF_LIQUID => D1(JCF_LIQUID(qdims%k_start) : JCF_LIQUID(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     CF_FROZEN => D1(JCF_FROZEN(qdims%k_start) : JCF_FROZEN(qdims%k_start)+ &
      field_length(theta_points,extended_halo,qdims%k_end-qdims%k_start+1) -1)

     ! size of cca varies according to L_3D_CCA - n_cca_lev set in dervsize
     CCA => D1(JCCA(1):JCCA(1) + &
                              field_length(theta_points,no_halo,n_cca_lev) -1)

     CCB          => D1(JCCB   : JCCB  +field_length(theta_points,no_halo,1) -1)

     CCT          => D1(JCCT   : JCCT  +field_length(theta_points,no_halo,1) -1)

     CCLWP        => D1(JCCLWP : JCCLWP+field_length(theta_points,no_halo,1) -1)

     DEEP_FLAG    => D1(JDEEPFLAG : JDEEPFLAG+ &
                                    field_length(theta_points,no_halo,1) -1)
     Past_precip  => D1(JPASTPRECIP : JPASTPRECIP + &
                                    field_length(theta_points,no_halo,1) -1)
     Past_conv_ht => D1(JPASTCONVHT : JPASTCONVHT + &
                                    field_length(theta_points,no_halo,1) -1)
     CANOPY_WATER => D1(JCANOPY_WATER : JCANOPY_WATER+ &
                                     field_length(land_points,no_halo,1) -1)

     LCBASE       => D1(JLCBASE : JLCBASE + &
                                    field_length(theta_points,no_halo,1) -1)

     CCW_RAD      => D1(JCCW_RAD(1) : JCCW_RAD(1) + &
                              field_length(theta_points,no_halo,qdims%k_end) -1)

!    Secondary Fields in D1

     EXNER_THETA_LEVELS => D1(JEXNER_THETA_LEVELS(tdims%k_start):     &
                              JEXNER_THETA_LEVELS(tdims%k_start)+     &
                      field_length(theta_points,single_halo,          &
                                   tdims%k_end-tdims%k_start+1) -1)

     P => D1(JP(pdims%k_start):JP(pdims%k_start) +                    &
                    field_length(theta_points,single_halo,            &
                                  pdims%k_end+1) -1)

     P_THETA_LEVELS => D1(JP_THETA_LEVELS(tdims%k_start):             &
                          JP_THETA_LEVELS(tdims%k_start)+             &
                      field_length(theta_points,single_halo,          &
                                   tdims%k_end-tdims%k_start+1) -1)

     PSTAR =>   D1(JPSTAR : JPSTAR +field_length(theta_points,no_halo,1) -1)

     SW_INCS => D1(JSW_INCS(0) : JSW_INCS(0) + &
                       field_length(theta_points, no_halo, model_levels+1)-1)

     LW_INCS => D1(JLW_INCS(0) : JLW_INCS(0) + &
                       field_length(theta_points, no_halo, model_levels  )-1)

!    Direct PAR flux for STOCHEM
     DIRPAR => D1(JDIRPAR : JDIRPAR+field_length(theta_points,no_halo,1) -1)

!    Soil Fields
     SMCL => D1(JSMCL(1):JSMCL(1)+ &
                              field_length(land_points,no_halo,sm_levels)-1)
     DEEP_SOIL_TEMP => D1(J_DEEP_SOIL_TEMP(1):J_DEEP_SOIL_TEMP(1)+ &
                              field_length(land_points,no_halo,sm_levels)-1)
     VOL_SMC_WILT   => D1(JVOL_SMC_WILT:JVOL_SMC_WILT+ &
                                      field_length(land_points,no_halo,1)-1)
     VOL_SMC_CRIT   => D1(JVOL_SMC_CRIT:JVOL_SMC_CRIT+ &
                                      field_length(land_points,no_halo,1)-1)
     VOL_SMC_SAT    => D1(JVOL_SMC_SAT:JVOL_SMC_SAT+ &
                                      field_length(land_points,no_halo,1)-1)
     SAT_SOIL_COND  => D1(JSAT_SOIL_COND:JSAT_SOIL_COND+ &
                                      field_length(land_points,no_halo,1)-1)
     THERM_CAP  =>D1(JTHERM_CAP:JTHERM_CAP + &
                                     field_length(land_points,no_halo,1) -1)
     THERM_COND =>D1(JTHERM_COND:JTHERM_COND+ &
                                      field_length(land_points,no_halo,1)-1)
     CLAPP_HORN =>D1(JCLAPP_HORN:JCLAPP_HORN+ &
                                      field_length(land_points,no_halo,1)-1)
     SAT_SOILW_SUCTION => D1(JSAT_SOILW_SUCTION:JSAT_SOILW_SUCTION+ &
                                      field_length(land_points,no_halo,1)-1)
     STHU => D1(JSTHU(1):JSTHU(1)+ &
                              field_length(land_points,no_halo,sm_levels)-1)
     STHF => D1(JSTHF(1):JSTHF(1)+ &
                              field_length(land_points,no_halo,sm_levels)-1)

!    Roughness lenght of sea points
     Z0         => D1(JZ0:JZ0           +field_length(theta_points,no_halo,1)-1)
     GS         => D1(JGS:JGS           +field_length(land_points,no_halo,1)-1)

      call cable_set_atm_fields( D1, LEN_TOT, land_points,no_halo,sm_levels,ntiles, &
                                 jTSOIL_TILE, jSMCL_TILE, jSTHF_TILE,            &
                                 jSNOW_DEPTH3L, jSNOW_MASS3L, jSNOW_TMP3L,       &
                                 jSNOW_RHO3L, jSNOW_RHO1L, jSNOW_AGE, jSNOW_FLG3L )
!    Orography Fields
     OROG_SIL => D1(JOROG_SIL : JOROG_SIL+field_length(land_points,no_halo,1)-1)
     OROG_HO2 => D1(JOROG_HO2 : JOROG_HO2+field_length(land_points,no_halo,1)-1)
     OROG_SD  => D1(JOROG_SD  : JOROG_SD +field_length(land_points,no_halo,1)-1)
     OROG_GRAD_X  => D1(JOROG_GRAD_X : JOROG_GRAD_X+ &
                                          field_length(land_points,no_halo,1)-1)
     OROG_GRAD_Y  => D1(JOROG_GRAD_Y : JOROG_GRAD_Y+ &
                                          field_length(land_points,no_halo,1)-1)
     OROG_UNFILT  => D1(JOROG_UNFILT : JOROG_UNFILT+ &
                                          field_length(land_points,no_halo,1)-1)
     OROG_GRAD_XX => D1(JOROG_GRAD_XX : JOROG_GRAD_XX+ &
                                          field_length(land_points,no_halo,1)-1)
     OROG_GRAD_XY => D1(JOROG_GRAD_XY : JOROG_GRAD_XY+ &
                                          field_length(land_points,no_halo,1)-1)
     OROG_GRAD_YY => D1(JOROG_GRAD_YY : JOROG_GRAD_YY+ &
                                          field_length(land_points,no_halo,1)-1)

!    Sea/Sea Ice Fields
     U_SEA => D1(JU_SEA : JU_SEA    +field_length(u_points,no_halo,1) -1)
     V_SEA => D1(JV_SEA : JV_SEA    +field_length(v_points,no_halo,1) -1)
     ICE_FRACTION  => D1(JICE_FRACTION  : JICE_FRACTION + &
                        field_length(theta_points_sea_only,no_halo,1) -1)
     ICE_THICKNESS => D1(JICE_THICKNESS : JICE_THICKNESS+ &
                        field_length(theta_points_sea_only,no_halo,1) -1)
     TI => D1(JTI : JTI+field_length(theta_points_sea_only,no_halo,1) -1)
     ICE_FRACT_CAT => D1(JICE_FRACT_CAT : JICE_FRACT_CAT+ &
                              field_length(theta_points,no_halo,nice) -1)
     ICE_THICK_CAT => D1(JICE_THICK_CAT : JICE_THICK_CAT+ &
                              field_length(theta_points,no_halo,nice) -1)
     TI_CAT => D1(JTI_CAT : JTI_CAT+field_length(theta_points,no_halo,nice) -1)
     ICE_K_CAT => D1(JICE_K_CAT : JICE_K_CAT + &
                              field_length(theta_points,no_halo,nice) -1)
     U_0_P => D1(JU_0_P : JU_0_P+field_length(theta_points,no_halo,1) -1)
     V_0_P => D1(JV_0_P : JV_0_P+field_length(theta_points,no_halo,1) -1)

!    Snow Fields
     SNODEP => D1(JSNODEP : JSNODEP+field_length(theta_points,no_halo,1) -1)
     SNODEP_SEA => D1(JSNODEP_SEA : JSNODEP_SEA+ &
                           field_length(theta_points_sea_only,no_halo,1) -1)
     SNODEP_SEA_CAT => D1(JSNODEP_SEA_CAT : JSNODEP_SEA_CAT+ &
                                 field_length(theta_points,no_halo,nice) -1)
! SNSOOT may not be used as of vn6.6
     SNSOOT => D1(JSNSOOT : JSNSOOT+field_length(theta_points,no_halo,1) -1)
     CATCH_SNOW => D1(JCATCH_SNOW : JCATCH_SNOW+ &
                                field_length(land_points,no_halo,ntiles) -1)
     SNOW_GRND => D1(JSNOW_GRND : JSNOW_GRND+ &
                                field_length(land_points,no_halo,ntiles) -1)

!    Decoupled screen temperatures
     TScrnDcl_TILE => D1(JTScrnDcl_TILE : JTScrnDcl_TILE+ &
                                     field_length(land_points,no_halo,ntiles) -1)
     TScrnDcl_SSI => D1(JTScrnDcl_SSI : JTScrnDcl_SSI+ &
                                    field_length(theta_points,no_halo,1) -1)
     tStbTrans => D1(JtStbTrans : JtStbTrans+ &
                                    field_length(theta_points,no_halo,1) -1)

!    OZONE (has extra surface level for V-AT-POLES)
     IF (lexpand_ozone) THEN
!      Ozone held as zonal averages, i.e. one value per row
       o3 => d1(jozone(o3dims2%k_start):jozone(o3dims2%k_start)+ &
                                  rows * (o3dims2%k_end-o3dims2%k_start+1) -1)
     ELSE
       o3 => d1(jozone(o3dims2%k_start):jozone(o3dims2%k_start)+ &
        field_length(ozone_points,no_halo,o3dims2%k_end-o3dims2%k_start+1) -1)
     END IF

!    Tropopause-based Ozone 
     IF (tpps_ozone_levels > 0) THEN
       IF (lexpand_tpps_ozone) THEN
         tppsozone => d1(jtppsozone(o3dims2%k_start): &
                       jtppsozone(o3dims2%k_start) + rows*tpps_ozone_levels -1)
       ELSE
         tppsozone => d1(jtppsozone(o3dims2%k_start): &
                                               jtppsozone(o3dims2%k_start)+ &
                       field_length(ozone_points,no_halo,tpps_ozone_levels) -1)
       END IF
     ELSE
       tppsozone => dummy_field
     END IF

!    Ozone tracer field and cariolle parameters
     OZONE_TRACER => &
               D1(JOZONE_TRACER(tdims%k_start):JOZONE_TRACER(tdims%k_start)+ &
       field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     O3_PROD_LOSS => &
               D1(JO3_PROD_LOSS(tdims%k_start):JO3_PROD_LOSS(tdims%k_start)+ &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)
     O3_P_L_VMR   => &
               D1(JO3_P_L_VMR(tdims%k_start)  :JO3_P_L_VMR(tdims%k_start)+   &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)
     O3_VMR       => &
               D1(JO3_VMR (tdims%k_start)     :JO3_VMR(tdims%k_start)+       &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)
     O3_P_L_TEMP  => &
               D1(JO3_P_L_TEMP(tdims%k_start) :JO3_P_L_TEMP(tdims%k_start)+  &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)
     O3_TEMP      => &
               D1(JO3_TEMP(tdims%k_start)     :JO3_TEMP(tdims%k_start)+      &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)    
     O3_P_L_COLO3 => &
               D1(JO3_P_L_COLO3(tdims%k_start):JO3_P_L_COLO3(tdims%k_start)+ &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)
     O3_COLO3     => &
               D1(JO3_COLO3(tdims%k_start)    :JO3_COLO3(tdims%k_start)+     &
                             (ROWS*(tdims%k_end-tdims%k_start+1) ) -1)

!    Sources and Aerosol Ancillaries
     MURK_SOURCE => D1(JMURK_SOURCE(tdims%k_start) : &
                                               JMURK_SOURCE(tdims%k_start)+ &
            field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     SO2_EM      => D1(JSO2_EM : JSO2_EM + &
                                      field_length(theta_points,no_halo,1) -1)
     DMS_EM      => D1(JDMS_EM : JDMS_EM + &
                                      field_length(theta_points,no_halo,1) -1)
     MURK        => D1(JMURK(tdims%k_start) : JMURK(tdims%k_start)+ &
        field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)

!    Sulphur cycle
     SO2         => D1(JSO2(tdims%k_start) : JSO2(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DMS         => D1(JDMS(tdims%k_start) : JDMS(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SO4_AITKEN  => D1(JSO4_AITKEN(tdims%k_start) : &
                                                JSO4_AITKEN(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SO4_ACCU    => D1(JSO4_ACCU(tdims%k_start) : JSO4_ACCU(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SO4_DISS    => D1(JSO4_DISS(tdims%k_start) : JSO4_DISS(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     H2O2        => D1(JH2O2(tdims%k_start) : JH2O2(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     NH3         => D1(JNH3(tdims%k_start) : JNH3(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SOOT_NEW    => D1(JSOOT_NEW(tdims%k_start) : JSOOT_NEW(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SOOT_AGD    => D1(JSOOT_AGD(tdims%k_start) : JSOOT_AGD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SOOT_CLD    => D1(JSOOT_CLD(tdims%k_start) : JSOOT_CLD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     BMASS_NEW   => D1(JBMASS_NEW(tdims%k_start) : &
                                                 JBMASS_NEW(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     BMASS_AGD   => D1(JBMASS_AGD(tdims%k_start) : &
                                                 JBMASS_AGD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     BMASS_CLD   => D1(JBMASS_CLD(tdims%k_start) : &
                                                 JBMASS_CLD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     OCFF_NEW    => D1(JOCFF_NEW(tdims%k_start) : JOCFF_NEW(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     OCFF_AGD    => D1(JOCFF_AGD(tdims%k_start) : JOCFF_AGD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     OCFF_CLD    => D1(JOCFF_CLD(tdims%k_start) : JOCFF_CLD(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     SO2_NATEM   => D1(JSO2_NATEM(tdims%k_start) : &
                                                 JSO2_NATEM(tdims%k_start)+ &
      field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     OH          => D1(JOH(tdims%k_start) : JOH(tdims%k_start)+ &
      field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     HO2         => D1(JHO2(tdims%k_start) : JHO2(tdims%k_start)+ &
      field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     H2O2_LIMIT  => D1(JH2O2_LIMIT(tdims%k_start) :&
                                              JH2O2_LIMIT(tdims%k_start)+ &
      field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     O3_CHEM     => D1(JO3_CHEM(tdims%k_start) : JO3_CHEM(tdims%k_start)+ &
      field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     SO2_HILEM   => D1(JSO2_HILEM : JSO2_HILEM+ &
      field_length(theta_points,no_halo,1) -1)
     NH3_EM      => D1(JNH3_EM : JNH3_EM+ &
      field_length(theta_points,no_halo,1) -1)
     SOOT_EM     => D1(JSOOT_EM : JSOOT_EM+ &
      field_length(theta_points,no_halo,1) -1)
     SOOT_HILEM  => D1(JSOOT_HILEM : JSOOT_HILEM+ &
      field_length(theta_points,no_halo,1) -1)
     BMASS_EM    => D1(JBMASS_EM : JBMASS_EM+ &
      field_length(theta_points,no_halo,1) -1)
     BMASS_HILEM => D1(JBMASS_HILEM : JBMASS_HILEM+ &
      field_length(theta_points,no_halo,1) -1)
     OCFF_EM     => D1(JOCFF_EM : JOCFF_EM+ &
      field_length(theta_points,no_halo,1) -1)
     OCFF_HILEM  => D1(JOCFF_HILEM : JOCFF_HILEM+ &
      field_length(theta_points,no_halo,1) -1)
     DMS_CONC    => D1(JDMS_CONC : JDMS_CONC+ &
      field_length(theta_points,no_halo,1) -1)
     DMS_OFLUX   => D1(JDMS_OFLUX : JDMS_OFLUX+ &
      field_length(theta_points,no_halo,1) -1)

! 
! Ammonium nitrate scheme: 
     NITR_ACC  => D1(JNITR_ACC(tdims%k_start) : JNITR_ACC(tdims%k_start)+ &  
       field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)  
     NITR_DISS => D1(JNITR_DISS(tdims%k_start) : JNITR_DISS(tdims%k_start)+ &  
       field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)  

! 
! Aerosol climatologies
     ARCLBIOG_BG => D1(JARCLBIOG_BG(tdims%k_start) : &
                           JARCLBIOG_BG(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLBIOM_FR => D1(JARCLBIOM_FR(tdims%k_start) : & 
                           JARCLBIOM_FR(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLBIOM_AG => D1(JARCLBIOM_AG(tdims%k_start) : & 
                           JARCLBIOM_AG(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLBIOM_IC => D1(JARCLBIOM_IC(tdims%k_start) : & 
                           JARCLBIOM_IC(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLBLCK_FR => D1(JARCLBLCK_FR(tdims%k_start) : & 
                           JARCLBLCK_FR(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLBLCK_AG => D1(JARCLBLCK_AG(tdims%k_start) : & 
                           JARCLBLCK_AG(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLSSLT_FI => D1(JARCLSSLT_FI(tdims%k_start) : & 
                           JARCLSSLT_FI(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLSSLT_JT => D1(JARCLSSLT_JT(tdims%k_start) : & 
                           JARCLSSLT_JT(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLSULP_AC => D1(JARCLSULP_AC(tdims%k_start) : & 
                           JARCLSULP_AC(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLSULP_AK => D1(JARCLSULP_AK(tdims%k_start) : & 
                           JARCLSULP_AK(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLSULP_DI => D1(JARCLSULP_DI(tdims%k_start) : & 
                           JARCLSULP_DI(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B1 => D1(JARCLDUST_B1(tdims%k_start) : & 
                           JARCLDUST_B1(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B2 => D1(JARCLDUST_B2(tdims%k_start) : & 
                           JARCLDUST_B2(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B3 => D1(JARCLDUST_B3(tdims%k_start) : & 
                           JARCLDUST_B3(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B4 => D1(JARCLDUST_B4(tdims%k_start) : & 
                           JARCLDUST_B4(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B5 => D1(JARCLDUST_B5(tdims%k_start) : & 
                           JARCLDUST_B5(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDUST_B6 => D1(JARCLDUST_B6(tdims%k_start) : & 
                           JARCLDUST_B6(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLOCFF_FR => D1(JARCLOCFF_FR(tdims%k_start) : & 
                           JARCLOCFF_FR(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLOCFF_AG => D1(JARCLOCFF_AG(tdims%k_start) : & 
                           JARCLOCFF_AG(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLOCFF_IC => D1(JARCLOCFF_IC(tdims%k_start) : & 
                           JARCLOCFF_IC(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)

     ARCLDLTA_DL => D1(JARCLDLTA_DL(tdims%k_start) : & 
                           JARCLDLTA_DL(tdims%k_start) + &
        field_length(theta_points,no_halo,tdims%k_end-tdims%k_start+1) -1)
     
!    Mineral Dust Scheme
     SOIL_CLAY  => D1(JSOIL_CLAY : JSOIL_CLAY+ &
      field_length(theta_points,no_halo,1) -1)
     SOIL_SILT  => D1(JSOIL_SILT : JSOIL_SILT+ &
      field_length(theta_points,no_halo,1) -1)
     SOIL_SAND  => D1(JSOIL_SAND : JSOIL_SAND+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL1 => D1(JDUST_MREL1 : JDUST_MREL1+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL2 => D1(JDUST_MREL2 : JDUST_MREL2+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL3 => D1(JDUST_MREL3 : JDUST_MREL3+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL4 => D1(JDUST_MREL4 : JDUST_MREL4+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL5 => D1(JDUST_MREL5 : JDUST_MREL5+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_MREL6 => D1(JDUST_MREL6 : JDUST_MREL6+ &
      field_length(theta_points,no_halo,1) -1)
     DUST_DIV1  => D1(JDUST_DIV1(tdims%k_start) : JDUST_DIV1(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DUST_DIV2  => D1(JDUST_DIV2(tdims%k_start) : JDUST_DIV2(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DUST_DIV3  => D1(JDUST_DIV3(tdims%k_start) : JDUST_DIV3(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DUST_DIV4  => D1(JDUST_DIV4(tdims%k_start) : JDUST_DIV4(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DUST_DIV5  => D1(JDUST_DIV5(tdims%k_start) : JDUST_DIV5(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)
     DUST_DIV6  => D1(JDUST_DIV6(tdims%k_start) : JDUST_DIV6(tdims%k_start)+ &
      field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)

!    Carbon Cycle
     CO2FLUX   => D1(J_CO2FLUX : J_CO2FLUX+ &
         field_length(theta_points,no_halo,1) -1)

     CO2_EMITS => D1(J_CO2_EMITS : J_CO2_EMITS+ &
         field_length(theta_points,no_halo,1) -1)

     CO2       => D1(JCO2(tdims%k_start):JCO2(tdims%k_start) + &
         field_length(theta_points,single_halo,tdims%k_end-tdims%k_start+1) -1)

!    level dependent constants
     zseak_theta => D1(JZSEAK_THETA : JZSEAK_THETA+(model_levels+1) -1)
     Ck_theta    => D1(JCK_THETA    : JCK_THETA   +(model_levels+1) -1)
     zseak_rho   => D1(JZSEAK_RHO   : JZSEAK_RHO  +(model_levels+1) -1)
     Ck_rho      => D1(JCK_RHO      : JCK_RHO     +(model_levels+1) -1)

!    User ancillaries
     USER_ANC1   => D1(JUSER_ANC1  : JUSER_ANC1 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC2   => D1(JUSER_ANC2  : JUSER_ANC2 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC3   => D1(JUSER_ANC3  : JUSER_ANC3 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC4   => D1(JUSER_ANC4  : JUSER_ANC4 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC5   => D1(JUSER_ANC5  : JUSER_ANC5 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC6   => D1(JUSER_ANC6  : JUSER_ANC6 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC7   => D1(JUSER_ANC7  : JUSER_ANC7 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC8   => D1(JUSER_ANC8  : JUSER_ANC8 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC9   => D1(JUSER_ANC9  : JUSER_ANC9 + &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC10  => D1(JUSER_ANC10 : JUSER_ANC10+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC11  => D1(JUSER_ANC11 : JUSER_ANC11+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC12  => D1(JUSER_ANC12 : JUSER_ANC12+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC13  => D1(JUSER_ANC13 : JUSER_ANC13+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC14  => D1(JUSER_ANC14 : JUSER_ANC14+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC15  => D1(JUSER_ANC15 : JUSER_ANC15+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC16  => D1(JUSER_ANC16 : JUSER_ANC16+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC17  => D1(JUSER_ANC17 : JUSER_ANC17+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC18  => D1(JUSER_ANC18 : JUSER_ANC18+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC19  => D1(JUSER_ANC19 : JUSER_ANC19+ &
      field_length(theta_points,no_halo,1) -1)
     USER_ANC20  => D1(JUSER_ANC20 : JUSER_ANC20+ &
      field_length(theta_points,no_halo,1) -1)
     USER_MULT1  => D1(JUSER_MULT1(1)  : JUSER_MULT1(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT2  => D1(JUSER_MULT2(1)  : JUSER_MULT2(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT3  => D1(JUSER_MULT3(1)  : JUSER_MULT3(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT4  => D1(JUSER_MULT4(1)  : JUSER_MULT4(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT5  => D1(JUSER_MULT5(1)  : JUSER_MULT5(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT6  => D1(JUSER_MULT6(1)  : JUSER_MULT6(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT7  => D1(JUSER_MULT7(1)  : JUSER_MULT7(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT8  => D1(JUSER_MULT8(1)  : JUSER_MULT8(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT9  => D1(JUSER_MULT9(1)  : JUSER_MULT9(1) + &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT10 => D1(JUSER_MULT10(1) : JUSER_MULT10(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT11 => D1(JUSER_MULT11(1) : JUSER_MULT11(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT12 => D1(JUSER_MULT12(1) : JUSER_MULT12(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT13 => D1(JUSER_MULT13(1) : JUSER_MULT13(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT14 => D1(JUSER_MULT14(1) : JUSER_MULT14(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT15 => D1(JUSER_MULT15(1) : JUSER_MULT15(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT16 => D1(JUSER_MULT16(1) : JUSER_MULT16(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT17 => D1(JUSER_MULT17(1) : JUSER_MULT17(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT18 => D1(JUSER_MULT18(1) : JUSER_MULT18(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT19 => D1(JUSER_MULT19(1) : JUSER_MULT19(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)
     USER_MULT20 => D1(JUSER_MULT20(1) : JUSER_MULT20(1)+ &
      field_length(theta_points,no_halo,model_levels) -1)

!    Tiled vegetation and triffid
     FRAC_TYP =>D1(JFRAC_TYP:JFRAC_TYP  + & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON1=>D1(JFRAC_CON1:JFRAC_CON1+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON2=>D1(JFRAC_CON2:JFRAC_CON2+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON3=>D1(JFRAC_CON3:JFRAC_CON3+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON4=>D1(JFRAC_CON4:JFRAC_CON4+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON5=>D1(JFRAC_CON5:JFRAC_CON5+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON6=>D1(JFRAC_CON6:JFRAC_CON6+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON7=>D1(JFRAC_CON7:JFRAC_CON7+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON8=>D1(JFRAC_CON8:JFRAC_CON8+ & 
                                       field_length(land_points,no_halo,1)-1)
     FRAC_CON9=>D1(JFRAC_CON9:JFRAC_CON9+ & 
                                       field_length(land_points,no_halo,1)-1)
     LAI_PFT  =>D1(JLAI_PFT:JLAI_PFT    + & 
                                       field_length(land_points,no_halo,1) -1)
     CANHT_PFT =>D1(JCANHT_PFT:JCANHT_PFT+ & 
                                       field_length(land_points,no_halo,1)-1)
     DISTURB_VEG =>  D1(JDISTURB:JDISTURB+ & 
                                       field_length(land_points,no_halo,1)-1)
     SOIL_ALB  =>  D1(JSOIL_ALB:JSOIL_ALB+ & 
                                       field_length(land_points,no_halo,1)-1)
     OBS_ALB_SW  =>  D1(JOBS_ALB_SW:JOBS_ALB_SW+ &  
                                       field_length(land_points,no_halo,1)-1) 
     OBS_ALB_VIS  =>  D1(JOBS_ALB_VIS:JOBS_ALB_VIS+ &  
                                       field_length(land_points,no_halo,1)-1) 
     OBS_ALB_NIR  =>  D1(JOBS_ALB_NIR:JOBS_ALB_NIR+ &  
                                       field_length(land_points,no_halo,1)-1) 
     SOIL_CARB =>D1(JSOIL_CARB:JSOIL_CARB+ & 
                                       field_length(land_points,no_halo,1)-1)
     SOIL_CARB1 => D1(JSOIL_CARB1:JSOIL_CARB1+ &
                                       field_length(land_points,no_halo,1) -1)
     SOIL_CARB2 => D1(JSOIL_CARB2:JSOIL_CARB2+ &
                                       field_length(land_points,no_halo,1) -1)
     SOIL_CARB3 => D1(JSOIL_CARB3:JSOIL_CARB3+ &
                                       field_length(land_points,no_halo,1) -1)
     SOIL_CARB4 => D1(JSOIL_CARB4:JSOIL_CARB4+ &
                                       field_length(land_points,no_halo,1) -1)
     NPP_PFT_ACC    => D1(JNPP_PFT_ACC:JNPP_PFT_ACC+ &
                                       field_length(land_points,no_halo,1) -1)
     G_LF_PFT_ACC   => D1(JG_LF_PFT_ACC:JG_LF_PFT_ACC+ &
                                       field_length(land_points,no_halo,1) -1)
     G_PHLF_PFT_ACC => D1(JG_PHLF_PFT_ACC:JG_PHLF_PFT_ACC+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_W_PFT_ACC  => D1(JRSP_W_PFT_ACC:JRSP_W_PFT_ACC+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_S_ACC  => D1(JRSP_S_ACC:JRSP_S_ACC+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_S_ACC1 => D1(JRSP_S_ACC1:JRSP_S_ACC1+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_S_ACC2 => D1(JRSP_S_ACC2:JRSP_S_ACC2+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_S_ACC3 => D1(JRSP_S_ACC3:JRSP_S_ACC3+ &
                                       field_length(land_points,no_halo,1) -1)
     RSP_S_ACC4 => D1(JRSP_S_ACC4:JRSP_S_ACC4+ &
                                       field_length(land_points,no_halo,1) -1)
     CAN_WATER_TILE => D1(JCAN_WATER_TILE:JCAN_WATER_TILE+ &
                                       field_length(land_points,no_halo,1) -1)
     CATCH_TILE  => D1(JCATCH_TILE:JCATCH_TILE+ &
                                       field_length(land_points,no_halo,1) -1)
     RGRAIN_TILE => D1(JRGRAIN_TILE:JRGRAIN_TILE+ &
                                       field_length(land_points,no_halo,1) -1)
     TSTAR_TILE  => D1(JTSTAR_TILE:JTSTAR_TILE+ &
                                         field_length(land_points,no_halo,1) -1)
     Z0_TILE     => D1(JZ0_TILE:JZ0_TILE+field_length(land_points,no_halo,1) -1)
     Z0H_TILE    => D1(JZ0H_TILE:JZ0H_TILE+field_length(land_points,no_halo,1) &
                      -1)
     SNODEP_TILE => D1(JSNODEP_TILE:JSNODEP_TILE+ &
                                         field_length(land_points,no_halo,1) -1)
     INFIL_TILE  => D1(JINFIL_TILE:JINFIL_TILE+ &
                                         field_length(land_points,no_halo, 1)-1)
     DOLR_FIELD  => D1(JDOLR:JDOLR      +field_length(theta_points,no_halo,1)-1)
     LW_DOWN     => D1(JLW_DOWN:JLW_DOWN+field_length(theta_points,no_halo,1)-1)
     SW_TILE_RTS => D1(JSW_TILE:JSW_TILE+field_length(land_points,no_halo, 1)-1)

! MORUSES - new urban two-tile scheme
      hgt   => d1(jurbhgt  :jurbhgt   +field_length(land_points,no_halo,1) -1)
      ! building height
      hwr   => d1(jurbhwr  :jurbhwr   +field_length(land_points,no_halo,1) -1)
      ! height to width
      wrr   => d1(jurbwrr  :jurbwrr   +field_length(land_points,no_halo,1) -1)
      ! width ratio
      disp  => d1(jurbdisp :jurbdisp  +field_length(land_points,no_halo,1) -1)
      ! displacement height
      ztm   => d1(jurbztm  :jurbztm   +field_length(land_points,no_halo,1) -1)
      !
      albwl => d1(jurbalbwl:jurbalbwl +field_length(land_points,no_halo,1) -1)
      ! wall albedo
      albrd => d1(jurbalbrd:jurbalbrd +field_length(land_points,no_halo,1) -1)
      ! road albedo
      emisw => d1(jurbemisw:jurbemisw +field_length(land_points,no_halo,1) -1)
      ! wall emissivity
      emisr => d1(jurbemisr:jurbemisr +field_length(land_points,no_halo,1) -1)
      ! road emissivity

!    River routing fields
     RIV_SEQUENCE  => D1(JRIV_SEQUENCE : JRIV_SEQUENCE+ &
                                     field_length(river_points,no_halo,1) -1)
     RIV_DIRECTION => D1(JRIV_DIRECTION : JRIV_DIRECTION+ &
                                     field_length(river_points,no_halo,1) -1)
     RIV_STORAGE   => D1(JRIV_STORAGE : JRIV_STORAGE+ &
                                     field_length(river_points,no_halo,1) -1)
     TOT_SURFROFF  => D1(JTOT_SURFROFF : JTOT_SURFROFF+ &
                                     field_length(land_points,no_halo,1) -1)
     TOT_SUBROFF   => D1(JTOT_SUBROFF : JTOT_SUBROFF+ &
                                     field_length(land_points,no_halo,1) -1)
     RIV_INLANDATM => D1(JRIV_INLANDATM : JRIV_INLANDATM+ &
                                     field_length(land_points,no_halo,1)  -1)
     ! these are uninitialised upon entering ATM_STEP
     RIV_IAREA     => dummy_field !D1(1:1+row_length*rows)
     RIV_SLOPE     => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWOBS1  => dummy_field !D1(1:1+row_length*rows)
     RIV_INEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_JNEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_LAND      => dummy_field !D1(1:1+row_length*rows)
     RIV_SUBSTORE  => dummy_field !D1(1:1+row_length*rows)
     RIV_SURFSTORE => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWIN    => dummy_field !D1(1:1+row_length*rows)
     RIV_BFLOWIN   => dummy_field !D1(1:1+row_length*rows)
 
!    Required for water conservation correction due to lake evaporation 
     ACC_LAKE_EVAP => D1(JACC_LAKE_EVAP:JACC_LAKE_EVAP                  & 
                   +field_length(theta_points,no_halo,1) -1)

!    Fields to be retained in dumps for coupled models using OASIS    
     C_SOLAR => D1(JC_SOLAR : JC_SOLAR + &
                    field_length(theta_points,no_halo,1) -1)

     C_BLUE =>  D1(JC_BLUE : JC_BLUE + &
                    field_length(theta_points,no_halo,1) -1)

     C_LONGWAVE => D1(JC_LONGWAVE : JC_LONGWAVE + &
                    field_length(theta_points,no_halo,1) -1)

     C_TAUX => D1(JC_TAUX : JC_TAUX + &
                    field_length(u_points,no_halo,1) -1)

     C_TAUY => D1(JC_TAUY : JC_TAUY + &
                    field_length(v_points,no_halo,1) -1)

     C_W10 => D1(JC_W10 : JC_W10 + &
                    field_length(theta_points,no_halo,1) -1)

     C_SENSIBLE => D1(JC_SENSIBLE : JC_SENSIBLE + &
                    field_length(theta_points,no_halo,1) -1)

     C_SUBLIM =>  D1(JC_SUBLIM : JC_SUBLIM + &
                    field_length(theta_points,no_halo,nice_use) -1)

     C_EVAP =>  D1(JC_EVAP : JC_EVAP + &
                    field_length(theta_points,no_halo,1) -1)

     C_FCONDTOPN => D1(JC_FCONDTOPN : JC_FCONDTOPN + &
                    field_length(theta_points,no_halo,nice) -1)

     C_TOPMELTN => D1(JC_TOPMELTN : JC_TOPMELTN + &
                    field_length(theta_points,no_halo,nice) -1)

     C_LSRAIN =>  D1(JC_LSRAIN : JC_LSRAIN + &
                    field_length(theta_points,no_halo,1) -1)

     C_LSSNOW =>  D1(JC_LSSNOW : JC_LSSNOW + &
                    field_length(theta_points,no_halo,1) -1)

     C_CVRAIN =>  D1(JC_CVRAIN : JC_CVRAIN + &
                    field_length(theta_points,no_halo,1) -1)

     C_CVSNOW =>  D1(JC_CVSNOW : JC_CVSNOW + &
                    field_length(theta_points,no_halo,1) -1)

     C_RIVEROUT => D1(JC_RIVEROUT : JC_RIVEROUT + &
                    field_length(theta_points,no_halo,1) -1)

     C_CALVING => D1(JC_CALVING : JC_CALVING + &
                    field_length(theta_points,no_halo,1) -1)
     
! JULES 2 prognostics
      SNOWDEPTH      =>  D1(JSNOWDEPTH     :JSNOWDEPTH     + &
        field_length(land_points,no_halo,ntiles) -1)
      RHO_SNOW_GRND  =>  D1(JRHO_SNOW_GRND :JRHO_SNOW_GRND + &
        field_length(land_points,no_halo,ntiles) -1)
      NSNOW          =>  D1(JNSNOW         :JNSNOW         + &
        field_length(land_points,no_halo,ntiles) -1)
      DS             =>  D1(JDS            :JDS            + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)
      SICE           =>  D1(JSICE          :JSICE          + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)
      SLIQ           =>  D1(JSLIQ          :JSLIQ          + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)
      TSNOWLAYER     =>  D1(JTSNOWLAYER    :JTSNOWLAYER    + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)
      RHO_SNOW       =>  D1(JRHO_SNOW      :JRHO_SNOW      + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)
      RGRAINL        =>  D1(JRGRAINL       :JRGRAINL       + &
        field_length(land_points,no_halo,ntiles*nsmax) -1)

! FLake lake scheme prognostics
      LAKE_DEPTH     =>  D1(JLAKE_DEPTH     :JLAKE_DEPTH     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_FETCH     =>  D1(JLAKE_FETCH     :JLAKE_FETCH     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_T_MEAN    =>  D1(JLAKE_T_MEAN    :JLAKE_T_MEAN    + &
        field_length(land_points,no_halo,1) -1)
      LAKE_T_MXL     =>  D1(JLAKE_T_MXL     :JLAKE_T_MXL     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_T_ICE     =>  D1(JLAKE_T_ICE     :JLAKE_T_ICE     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_H_MXL     =>  D1(JLAKE_H_MXL     :JLAKE_H_MXL     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_H_ICE     =>  D1(JLAKE_H_ICE     :JLAKE_H_ICE     + &
        field_length(land_points,no_halo,1) -1)
      LAKE_SHAPE     =>  D1(JLAKE_SHAPE     :JLAKE_SHAPE     + &
        field_length(land_points,no_halo,1))
      LAKE_G_DT      =>  D1(JLAKE_G_DT      :JLAKE_G_DT      + &
        field_length(land_points,no_halo,1))

!    Required for energy correction
     NET_FLUX  => D1(JNET_FLUX:JNET_FLUX  + &
                                     field_length(theta_points,no_halo,1) -1)
     NET_MFLUX => D1(JNET_MFLUX:JNET_MFLUX+&
                                     field_length(theta_points,no_halo,1) -1)

!    Fields carried forward from previous version
     TSTAR_ANOM => D1(JTSTAR_ANOM : JTSTAR_ANOM + &
                                     field_length(theta_points,no_halo,1) -1)

!    lateral boundary conditions

     OROG_LBC  => D1(JOROG_LBC : JOROG_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)*1  -1)

     U_LBC     => D1(JU_LBC : JU_LBC + &
       LENRIMA(fld_type_u,halo_type_extended,1)* &
                                       (udims_l%k_end-udims_l%k_start+1) -1)

     V_LBC     => D1(JV_LBC : JV_LBC +  &
       LENRIMA(fld_type_v,halo_type_extended,1)* &
                                       (vdims_l%k_end-vdims_l%k_start+1) -1)

     ! model_levels+1
     W_LBC     => D1(JW_LBC : JW_LBC + &
      LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)* &
                                       (wdims_l%k_end-wdims_l%k_start+1) -1)

     RHO_LBC   => D1(JRHO_LBC : JRHO_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (pdims_l%k_end-pdims_l%k_start+1) -1)

     THETA_LBC => D1(JTHETA_LBC : JTHETA_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (tdims_l%k_end-tdims_l%k_start+1) -1)

     Q_LBC     => D1(JQ_LBC : JQ_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QCL_LBC   => D1(JQCL_LBC : JQCL_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QCF_LBC   => D1(JQCF_LBC : JQCF_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     ! model_levels+1
     EXNER_LBC => D1(JEXNER_LBC : JEXNER_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)*(pdims_l%k_end+1) -1) 

     U_ADV_LBC => D1(JU_ADV_LBC : JU_ADV_LBC + &
       LENRIMA(fld_type_u,halo_type_extended,1)* &
                                        (udims_l%k_end-udims_l%k_start+1) -1)

     V_ADV_LBC => D1(JV_ADV_LBC : JV_ADV_LBC + &
       LENRIMA(fld_type_v,halo_type_extended,1)* &
                                        (vdims_l%k_end-vdims_l%k_start+1) -1)

     ! model_levels+1
     W_ADV_LBC => D1(JW_ADV_LBC : JW_ADV_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (wdims_l%k_end-wdims_l%k_start+1) -1) 

     QCF2_LBC  => D1(JQCF2_LBC : JQCF2_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QRAIN_LBC => D1(JQRAIN_LBC : JQRAIN_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QGRAUP_LBC  => D1(JQGRAUP_LBC : JQGRAUP_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_BULK_LBC => D1(JCF_BULK_LBC : JCF_BULK_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_LIQUID_LBC => D1(JCF_LIQUID_LBC : JCF_LIQUID_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_FROZEN_LBC => D1(JCF_FROZEN_LBC : JCF_FROZEN_LBC + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     MURK_LBC  => D1(JMURK_LBC : JMURK_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV1_LBC => D1(JDUST_DIV1_LBC : JDUST_DIV1_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV2_LBC => D1(JDUST_DIV2_LBC : JDUST_DIV2_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV3_LBC => D1(JDUST_DIV3_LBC : JDUST_DIV3_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV4_LBC => D1(JDUST_DIV4_LBC : JDUST_DIV4_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV5_LBC => D1(JDUST_DIV5_LBC : JDUST_DIV5_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV6_LBC => D1(JDUST_DIV6_LBC : JDUST_DIV6_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO2_LBC  => D1(JSO2_LBC : JSO2_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     DMS_LBC => D1(JDMS_LBC : JDMS_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_AITKEN_LBC => D1(JSO4_AITKEN_LBC : JSO4_AITKEN_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_ACCU_LBC => D1(JSO4_ACCU_LBC : JSO4_ACCU_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_DISS_LBC => D1(JSO4_DISS_LBC : JSO4_DISS_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     NH3_LBC => D1(JNH3_LBC : JNH3_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_NEW_LBC => D1(JSOOT_NEW_LBC : JSOOT_NEW_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_AGD_LBC => D1(JSOOT_AGD_LBC : JSOOT_AGD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_CLD_LBC => D1(JSOOT_CLD_LBC : JSOOT_CLD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_NEW_LBC => D1(JBMASS_NEW_LBC : JBMASS_NEW_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_AGD_LBC => D1(JBMASS_AGD_LBC : JBMASS_AGD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_CLD_LBC => D1(JBMASS_CLD_LBC : JBMASS_CLD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_NEW_LBC => D1(JOCFF_NEW_LBC : JOCFF_NEW_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_AGD_LBC => D1(JOCFF_AGD_LBC : JOCFF_AGD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_CLD_LBC => D1(JOCFF_CLD_LBC : JOCFF_CLD_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     NITR_ACC_LBC => D1(JNITR_ACC_LBC : JNITR_ACC_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     NITR_DISS_LBC => D1(JNITR_DISS_LBC : JNITR_DISS_LBC + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                       (tdims_s%k_end-tdims_s%k_start+1) -1)

     U_LBC_TEND     => D1(JU_LBC_TEND : JU_LBC_TEND + &
       LENRIMA(fld_type_u,halo_type_extended,1)* &
                                       (udims_l%k_end-udims_l%k_start+1) -1)

     V_LBC_TEND     => D1(JV_LBC_TEND : JV_LBC_TEND + &
       LENRIMA(fld_type_v,halo_type_extended,1)* &
                                       (vdims_l%k_end-vdims_l%k_start+1) -1)

     ! model_levels+1
     W_LBC_TEND     => D1(JW_LBC_TEND : JW_LBC_TEND + &
      LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)* &
                                       (wdims_l%k_end-wdims_l%k_start+1) -1)

     RHO_LBC_TEND   => D1(JRHO_LBC_TEND : JRHO_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (pdims_l%k_end-pdims_l%k_start+1) -1)

     THETA_LBC_TEND => D1(JTHETA_LBC_TEND : JTHETA_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (tdims_l%k_end-tdims_l%k_start+1) -1)

     Q_LBC_TEND     => D1(JQ_LBC_TEND : JQ_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QCL_LBC_TEND   => D1(JQCL_LBC_TEND : JQCL_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     QCF_LBC_TEND   => D1(JQCF_LBC_TEND : JQCF_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                       (qdims_l%k_end-qdims_l%k_start+1) -1)

     ! model_levels+1
     EXNER_LBC_TEND => D1(JEXNER_LBC_TEND : JEXNER_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)*(pdims_l%k_end+1) -1)

     U_ADV_LBC_TEND => D1(JU_ADV_LBC_TEND : JU_ADV_LBC_TEND + &
       LENRIMA(fld_type_u,halo_type_extended,1)* &
                                        (udims_l%k_end-udims_l%k_start+1) -1)

     V_ADV_LBC_TEND => D1(JV_ADV_LBC_TEND : JV_ADV_LBC_TEND + &
       LENRIMA(fld_type_v,halo_type_extended,1)* &
                                        (vdims_l%k_end-vdims_l%k_start+1) -1)

     ! model_levels+1
     W_ADV_LBC_TEND => D1(JW_ADV_LBC_TEND : JW_ADV_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (wdims_l%k_end-wdims_l%k_start+1) -1)

     QCF2_LBC_TEND => D1(JQCF2_LBC_TEND : JQCF2_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     QRAIN_LBC_TEND => D1(JQRAIN_LBC_TEND : JQRAIN_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     QGRAUP_LBC_TEND => D1(JQGRAUP_LBC_TEND : JQGRAUP_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_BULK_LBC_TEND => D1(JCF_BULK_LBC_TEND : JCF_BULK_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_LIQUID_LBC_TEND => D1(JCF_LIQUID_LBC_TEND : JCF_LIQUID_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     CF_FROZEN_LBC_TEND => D1(JCF_FROZEN_LBC_TEND : JCF_FROZEN_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_extended,1)* &
                                        (qdims_l%k_end-qdims_l%k_start+1) -1)

     MURK_LBC_TEND => D1(JMURK_LBC_TEND : JMURK_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV1_LBC_TEND => D1(JDUST_DIV1_LBC_TEND : JDUST_DIV1_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV2_LBC_TEND => D1(JDUST_DIV2_LBC_TEND : JDUST_DIV2_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV3_LBC_TEND => D1(JDUST_DIV3_LBC_TEND : JDUST_DIV3_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV4_LBC_TEND => D1(JDUST_DIV4_LBC_TEND : JDUST_DIV4_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV5_LBC_TEND => D1(JDUST_DIV5_LBC_TEND : JDUST_DIV5_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DUST_DIV6_LBC_TEND => D1(JDUST_DIV6_LBC_TEND : JDUST_DIV6_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO2_LBC_TEND  => D1(JSO2_LBC_TEND : JSO2_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     DMS_LBC_TEND => D1(JDMS_LBC_TEND : JDMS_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_AITKEN_LBC_TEND => D1(JSO4_AITKEN_LBC_TEND : JSO4_AITKEN_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_ACCU_LBC_TEND => D1(JSO4_ACCU_LBC_TEND : JSO4_ACCU_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SO4_DISS_LBC_TEND => D1(JSO4_DISS_LBC_TEND : JSO4_DISS_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     NH3_LBC_TEND => D1(JNH3_LBC_TEND : JNH3_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_NEW_LBC_TEND => D1(JSOOT_NEW_LBC_TEND : JSOOT_NEW_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_AGD_LBC_TEND => D1(JSOOT_AGD_LBC_TEND : JSOOT_AGD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     SOOT_CLD_LBC_TEND => D1(JSOOT_CLD_LBC_TEND : JSOOT_CLD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_NEW_LBC_TEND => D1(JBMASS_NEW_LBC_TEND : JBMASS_NEW_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_AGD_LBC_TEND => D1(JBMASS_AGD_LBC_TEND : JBMASS_AGD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     BMASS_CLD_LBC_TEND => D1(JBMASS_CLD_LBC_TEND : JBMASS_CLD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_NEW_LBC_TEND => D1(JOCFF_NEW_LBC_TEND : JOCFF_NEW_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_AGD_LBC_TEND => D1(JOCFF_AGD_LBC_TEND : JOCFF_AGD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     OCFF_CLD_LBC_TEND => D1(JOCFF_CLD_LBC_TEND : JOCFF_CLD_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     NITR_ACC_LBC_TEND => D1(JNITR_ACC_LBC_TEND : JNITR_ACC_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)

     NITR_DISS_LBC_TEND => D1(JNITR_DISS_LBC_TEND : JNITR_DISS_LBC_TEND + &
       LENRIMA(fld_type_p,halo_type_single,1)* &
                                        (tdims_s%k_end-tdims_s%k_start+1) -1)


! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle and ammonium nitrate scheme (these are in Section 33):
      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         OH_UKCA   => D1(JOH_UKCA(tdims%k_start) : JOH_UKCA(tdims%k_start)+ &
                       field_length(theta_points,no_halo,     &
                                            tdims%k_end-tdims%k_start+1) -1)
         H2O2_UKCA => D1(JH2O2_UKCA(tdims%k_start): JH2O2_UKCA(tdims%k_start)+ &
                       field_length(theta_points,single_halo, &
                                            tdims%k_end-tdims_s%k_start+1) -1)
         HO2_UKCA  => D1(JHO2_UKCA(tdims%k_start): JHO2_UKCA(tdims%k_start)+ &
                       field_length(theta_points,no_halo,     &
                                            tdims%k_end-tdims%k_start+1) -1)
         O3_UKCA   => D1(JO3_UKCA(tdims%k_start)   : JO3_UKCA(tdims%k_start)+ &
                       field_length(theta_points,single_halo, &
                                            tdims%k_end-tdims_s%k_start+1) -1)
         HNO3_UKCA => D1(JHNO3_UKCA(tdims%k_start):JHNO3_UKCA(tdims%k_start)+ &
                       field_length(theta_points,single_halo, &
                                            tdims%k_end-tdims_s%k_start+1) -1)
      ELSE 
        OH_UKCA   => dummy_field
        H2O2_UKCA => dummy_field
        HO2_UKCA  => dummy_field
        O3_UKCA   => dummy_field
        HNO3_UKCA => dummy_field
      END IF  

! 1.1 point tracer fields to D1

      ! find out how many tracers are active
      nActiveTracers=0
      DO nTracer=A_TRACER_FIRST,A_TRACER_LAST
        IF (SI(nTracer,33,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer

      IF (nActiveTracers /= 0) THEN

        ! set the pointer to the appropriate section of D1
        TRACER => D1( JTRACER(tdims%k_start,A_TRACER_FIRST) : &
                      JTRACER(tdims%k_start,A_TRACER_FIRST) + &
          field_length(theta_points,single_halo,tr_levels)*nActiveTracers -1)

      ELSE
        ! or set it to something non-null if there are no active tracers
        TRACER => dummy_field
      END IF

      ! do the same for section 34 (UKCA) tracers     
      nActiveTracers=0
      DO nTracer=A_UKCA_FIRST,A_UKCA_LAST
        IF (SI(nTracer,34,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer 

      IF (nActiveTracers /= 0) THEN

        TRACER_UKCA => D1( JTR_UKCA(tdims%k_start,A_UKCA_FIRST) : &
                           JTR_UKCA(tdims%k_start,A_UKCA_FIRST) + &
          field_length(theta_points,single_halo,tr_levels)*nActiveTracers -1)

      ELSE
        TRACER_UKCA => dummy_field
      END IF

     ! find out how many free tracer LBCs are active
      nActiveTracers=0
      DO nTracer=1,TR_LBC_VARS
        IF (SI(nTracer,36,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer

      IF (nActiveTracers /= 0) THEN
        ! set the pointer to the appropriate section of D1
        TRACER_LBC => D1(JTRACER_LBC(1) : &
         JTRACER_LBC(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*(tr_levels-tdims%k_start+1)))
        TRACER_LBC_TEND => D1(JTRACER_LBC_TEND(1) : &
         JTRACER_LBC_TEND(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*(tr_levels-tdims%k_start+1)))
      ELSE
        ! or set it to something non-null if there are no active tracer LBCs
        TRACER_LBC => dummy_field
        TRACER_LBC_TEND => dummy_field
      END IF

      ! find out how many UKCA tracer LBCs are active
      nActiveTracers=0
      DO nTracer=1,TR_LBC_UKCA
        IF (SI(nTracer,37,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer

      IF (nActiveTracers /= 0) THEN
        ! set the pointer to the appropriate section of D1
        TRACER_UKCA_LBC => D1(JTR_UKCA_LBC(1) : &
         JTR_UKCA_LBC(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*(tr_levels-tdims%k_start+1)))
        TRACER_UKCA_LBC_TEND => D1(JTR_UKCA_LBC_TEND(1) : &
         JTR_UKCA_LBC_TEND(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*(tr_levels-tdims%k_start+1)))
      ELSE
        ! set to something non-null if there are no active ukca tracer LBCs
        TRACER_UKCA_LBC => dummy_field
        TRACER_UKCA_LBC_TEND => dummy_field
      END IF
      IF (lhook) CALL dr_hook('SET_ATM_FIELDS',zhook_out,zhook_handle)
      RETURN

END SUBROUTINE Set_Atm_Fields
