! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL Purpose : To provide the interface for PHY_DIAG
!LL
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL                       version no. 1, dated 15/01/90
!LL
!LL Logical components covered : D4
!LL
!LL System task : P0
!LL
!LL Documentation : Unified Model Documentation Paper No P0
!LL
!LLEND---------------------------------------------------------------
!*L Arguments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: stash

      SUBROUTINE ST_DIAG2( INT16,                                       &
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
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
     &                    t_inc_pres,                                   & 
     &                    q_inc_pres,  qcl_inc_pres, qcf_inc_pres,      & 
     &                    cf_inc_pres, cfl_inc_pres, cff_inc_pres,      & 
     &                    t_dini,                                       & 
     &                    q_dini,  qcl_dini, qcf_dini,                  & 
     &                    cf_dini, cfl_dini, cff_dini,                  &
     &                    ICODE,CMESSAGE)

! needed by typ_atm_fields.h:
  USE atm_fields_bounds_mod

Use trignometric_mod,  Only : sec_theta_latitude
      Use ancil_info, only:                                             &
     &                      nsmax

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      USE phy_diag_mod, ONLY: phy_diag
      USE um_input_control_mod,  ONLY: model_domain, npmsl_height,      &
           l_pmsl_sor
      USE Submodel_Mod


  USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

!*L
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
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start typ_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "typptra.h", and requires that "ctracera.h" is
!  also included.  
!
! This file belongs in section: Top Level
!
! Note: including this file requires lbc_mod.
!
! ENDGame

REAL ::                                                                 &
  exner_surf(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end),                            &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows),                                          &
  etadot(wdims_s%i_start:wdims_s%i_end,                                 &
         wdims_s%j_start:wdims_s%j_end,                                 &
         wdims_s%k_start:wdims_s%k_end),                                &
  m_v   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cl  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_r   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_gr  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf2 (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end)
REAL :: thetav(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)
REAL :: dryrho(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,                           &
               pdims_s%k_start:pdims_s%k_end)

REAL ::  Exner      (pdims_s%i_start:pdims_s%i_end,                     &
                 pdims_s%j_start:pdims_s%j_end,                         &
                 pdims_s%k_start:pdims_s%k_end+1)

      ! 1.1: Data variables stored in primary space.
      real :: U    (udims_s%i_start:udims_s%i_end,                       &
                    udims_s%j_start:udims_s%j_end,                       &
                    udims_s%k_start:udims_s%k_end)

      real :: V    (vdims_s%i_start:vdims_s%i_end,                       &
                    vdims_s%j_start:vdims_s%j_end,                       &
                    vdims_s%k_start:vdims_s%k_end)
      real :: E_TRB  (tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: TSQ_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: QSQ_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: COV_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: ZHPAR_SHCU(row_length, rows)

      real :: W    (wdims_s%i_start:wdims_s%i_end,                       &
                    wdims_s%j_start:wdims_s%j_end,                       &
                    wdims_s%k_start:wdims_s%k_end)

      real :: RHO  (pdims_s%i_start:pdims_s%i_end,                       &
                    pdims_s%j_start:pdims_s%j_end,                       &
                    pdims_s%k_start:pdims_s%k_end)

      real :: THETA(tdims_s%i_start:tdims_s%i_end,                       &
                    tdims_s%j_start:tdims_s%j_end,                       &
                    tdims_s%k_start:tdims_s%k_end)

      real :: Q     (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCL   (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCF   (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCF2  (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QRAIN (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QGRAUP(qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

! EXNER_RHO_LEVELS is still 1D
! Exner pressure on rho levels
      real :: EXNER_RHO_LEVELS(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))   

! The next few fields are targets as they are used in multi-variate 
! swap_bounds calls

      real, target ::                                                    &
            U_ADV(udims_l%i_start:udims_l%i_end,                         &
                  udims_l%j_start:udims_l%j_end,                         &
                  udims_l%k_start:udims_l%k_end)

      real, target ::                                                    &
            V_ADV(vdims_l%i_start:vdims_l%i_end,                         &
                  vdims_l%j_start:vdims_l%j_end,                         &
                  vdims_l%k_start:vdims_l%k_end)
      real, target ::                                                    &
            W_ADV(wdims_l%i_start:wdims_l%i_end,                         &
                  wdims_l%j_start:wdims_l%j_end,                         &
                  wdims_l%k_start:wdims_l%k_end)

      ! 1.2: Data variables stored in secondary space.
      ! P is still 1D
      real :: P( (pdims_s%i_end-pdims_s%i_start+1)* &
                 (pdims_s%j_end-pdims_s%j_start+1)* &
                 (pdims_s%k_end+1) )

      ! Pressure on theta levels
      real :: P_THETA_LEVELS(tdims_s%i_start:tdims_s%i_end,              &
                             tdims_s%j_start:tdims_s%j_end,              &
                             tdims_s%k_start:tdims_s%k_end)

      ! Exner pressure on theta levels
      real :: EXNER_THETA_LEVELS(tdims_s%i_start:tdims_s%i_end,          &
                                 tdims_s%j_start:tdims_s%j_end,          &
                                 tdims_s%k_start:tdims_s%k_end)

      ! 1.3: Cloud Fields
      real :: CCW_RAD(qdims%i_start:qdims%i_end,                         &
                      qdims%j_start:qdims%j_end,                         &
                                    qdims%k_end)

      ! CCAs size is dependant on L_3D_CCA
      ! N_CCA_LEV will be set to the correct value (either wet_levels or 1)

      real :: CCA(qdims%i_start:qdims%i_end,                             &
                  qdims%j_start:qdims%j_end,                             &
                                  n_cca_lev)
      real :: CF_AREA  (qdims%i_start:qdims%i_end,                       &
                        qdims%j_start:qdims%j_end,                       &
                                      qdims%k_end)
      real :: CF_BULK  (qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID(qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN(qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)

      ! 1.4: Soil Ancillary fields
      real :: DEEP_SOIL_TEMP(land_field,sm_levels)
      real :: SMCL(land_field,sm_levels)
      real :: STHU(land_field,sm_levels)
      real :: STHF(land_field,sm_levels)

      ! 1.5: Radiation Increments
      real :: SW_INCS(rdims2%i_start:rdims2%i_end,                       &
                      rdims2%j_start:rdims2%j_end,                       &
                      rdims2%k_start:rdims2%k_end)   ! SW radiation increments
      real :: LW_INCS(tdims%i_start:tdims%i_end,                         &
                      tdims%j_start:tdims%j_end,                         &
                      tdims%k_start:tdims%k_end)   ! LW radiation increments

      ! PAR radiation increment
      real :: DIRPAR(row_length,rows)

      ! 1.6: Ozone
      real :: O3(o3dims%i_start:o3dims%i_end)
!  tropopause-based ozone
      real :: TPPSOZONE(o3dims%i_start:o3dims%i_end)
!  Cariolle ozone tracer variables
      real :: OZONE_TRACER(tdims_s%i_start:tdims_s%i_end,                &
                           tdims_s%j_start:tdims_s%j_end,                &
                           tdims_s%k_start:tdims_s%k_end)
      real :: O3_PROD_LOSS(1,rows,model_levels)
      real :: O3_P_L_VMR  (1,rows,model_levels)
      real :: O3_VMR      (1,rows,model_levels)
      real :: O3_P_L_TEMP (1,rows,model_levels)
      real :: O3_TEMP     (1,rows,model_levels)
      real :: O3_P_L_COLO3(1,rows,model_levels)
      real :: O3_COLO3    (1,rows,model_levels)
      
! TRACERS are still 1D
      ! Tracer and aerosol fields
      ! TRACERS are dealt w/ differently
      ! these are the maximum sizes:
      real :: TRACER              (tr_levels*theta_off_size*(a_tracer_last-a_tracer_first))
      real :: TRACER_UKCA         ((tdims%k_end-tdims%k_start+1)*theta_off_size*tr_ukca)
      real :: TRACER_LBC          (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_LBC_TEND     (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_UKCA_LBC     (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_ukca)
      real :: TRACER_UKCA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_ukca)

      real :: MURK_SOURCE(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)

      real :: MURK       (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV1  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV2  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV3  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV4  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV5  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV6  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO2        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DMS        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_AITKEN (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_ACCU   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_DISS   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: H2O2       (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: NH3        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_NEW  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_AGD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_CLD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO2_NATEM  (tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)

      real :: OH        (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: HO2       (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: H2O2_LIMIT(tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: O3_CHEM   (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: CO2       (tdims_s%i_start:tdims_s%i_end,                  &
                         tdims_s%j_start:tdims_s%j_end,                  &
                         tdims_s%k_start:tdims_s%k_end)

! USER_MULT<N> are still 1D
! 1.8: Multi-level user ancillary fields

      real :: USER_MULT1 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT2 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT3 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT4 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT5 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT6 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT7 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT8 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT9 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT10(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT11(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT12(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT13(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT14(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT15(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT16(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT17(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT18(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT19(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT20(oneddims%i_start:oneddims%i_end)


      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      real ::      OROG_LBC(LENRIMA(fld_type_p,halo_type_extended,1))
      real ::         U_LBC(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::         V_LBC(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::         W_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::       RHO_LBC(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end)
      real ::     THETA_LBC(LENRIMA(fld_type_p,halo_type_extended,1),tdims_s%k_start:tdims_s%k_end)
      real ::         Q_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCL_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCF_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::      QCF2_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     QRAIN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::    QGRAUP_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::   CF_BULK_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     EXNER_LBC(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end+1)
      real ::     U_ADV_LBC(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::     V_ADV_LBC(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::     W_ADV_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::      MURK_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV1_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV2_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV3_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV4_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV5_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV6_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       SO2_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       DMS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: SO4_AITKEN_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_ACCU_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_DISS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::        NH3_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   NITR_ACC_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  NITR_DISS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)

      ! Lateral Boundary Condition tendencies
      real ::         U_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::         V_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::         W_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::       RHO_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end)
      real ::     THETA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),tdims_s%k_start:tdims_s%k_end)
      real ::         Q_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCL_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCF_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::      QCF2_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     QRAIN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::    QGRAUP_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::   CF_BULK_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     EXNER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end+1)
      real ::     U_ADV_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::     V_ADV_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::     W_ADV_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::      MURK_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV1_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV2_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV3_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV4_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV5_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV6_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       SO2_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       DMS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: SO4_AITKEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_ACCU_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_DISS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::        NH3_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   NITR_ACC_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  NITR_DISS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)

      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
      real :: TSTAR     (row_length,rows)
      logical :: LAND      (row_length,rows)
      real :: TSTAR_ANOM(row_length,rows)

      ! 2.15: Fields for coastal tiling
      real :: FRAC_LAND (land_field)
      real :: TSTAR_LAND(row_length,rows)
      real :: TSTAR_SEA (row_length,rows)
      ! these fields are targets - see note below relating to TI etc.
      real, target :: TSTAR_SICE(row_length,rows,1)
      real, target :: TSTAR_SICE_CAT(row_length,rows,nice_use)

      ! Set pointers for sea-ice and land albedos
      real :: SICE_ALB(row_length,rows)
      real :: LAND_ALB(row_length,rows)

      ! 2.2: Data variables stored in secondary space.

      real :: PSTAR(row_length,rows)

      ! 2.3: Cloud fields
      real ::      LCBASE(row_length,rows)
      real ::         CCB(row_length,rows)
      real ::         CCT(row_length,rows)
      real ::       CCLWP(row_length,rows)
      real ::   DEEP_FLAG(row_length,rows)
      real :: PAST_PRECIP(row_length,rows)
      real :: PAST_CONV_HT(row_length,rows)

      ! 2.4: Boundary layer fields

      real :: ZH(row_length,rows)
 
!     Convective downdraught mass-flux at cloud base 
      REAL :: ddmfx(row_length, rows) 

      ! Standard deviation of turbulent fluctuations of layer 1 temperature
      real :: T1_SD(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      real :: Q1_SD(row_length,rows)

      ! Decoupled screen-level temperatures
      REAL ::  TScrnDcl_SSI(row_length,rows)
      REAL :: TScrnDcl_TILE(land_field,ntiles)
      REAL ::     tStbTrans(row_length,rows)

      ! Number of model levels in the  turbulently mixed layer
      INTEGER :: NTML(row_length,rows)

      ! Top level for turb mixing in any decoupled Sc layer
      INTEGER :: NTDSC(row_length,rows)

      ! Bottom level for turb mixing in any decoupled Sc layer
      INTEGER :: NBDSC(row_length,rows)

      LOGICAL :: CUMULUS(row_length,rows)

      ! 2.4: Soil Ancillary fields

      real :: SAT_SOILW_SUCTION(land_field)
      real :: THERM_CAP        (land_field)
      real :: THERM_COND       (land_field)
      real :: VOL_SMC_CRIT     (land_field)
      real :: VOL_SMC_WILT     (land_field)
      real :: VOL_SMC_SAT      (land_field)
      real :: SAT_SOIL_COND    (land_field)
      real :: CLAPP_HORN       (land_field)

      ! 2.5: Other surface fields
      real :: CANOPY_WATER(land_field)
      real :: Z0          (row_length,rows)
      real :: GS          (land_field)

      ! 2.6: Orographic Ancillary fields

      real :: OROGRAPHY   (row_length,rows)
      real :: OROG_SD     (land_field)
      real :: OROG_SIL    (land_field)
      real :: OROG_HO2    (land_field)
      real :: OROG_GRAD_X (land_field)
      real :: OROG_GRAD_Y (land_field)
      real :: OROG_UNFILT (land_field)
      real :: OROG_GRAD_XX(land_field)
      real :: OROG_GRAD_XY(land_field)
      real :: OROG_GRAD_YY(land_field)

      ! 2.7: Sea/Sea Ice

      real :: U_SEA(row_length,rows)
      real :: V_SEA(row_length,n_rows)
      real :: U_0_P(row_length,rows)
      real :: V_0_P(row_length,rows)

! these fields are targets b/c there are pointers w/in ATM_STEP
! which point to one or another of them depending on the ice category selection
      real, target :: TI            (row_length,rows,1)
      real, target :: ICE_FRACTION  (row_length,rows,1)
      real, target :: ICE_THICKNESS (row_length,rows,1)
      real, target :: TI_CAT        (row_length,rows,nice) 
      real, target :: ICE_FRACT_CAT (row_length,rows,nice)
      real, target :: ICE_THICK_CAT (row_length,rows,nice)
      real, target :: SNODEP_SEA    (row_length,rows,1)
      real, target :: SNODEP_SEA_CAT(row_length,rows,nice_use)

! Effective conductivity is not a target as field only exists in D1 if nice>1
      real :: ICE_K_CAT     (row_length,rows,nice)

      ! 2.8: Snow

      real :: SNODEP         (row_length,rows)
      real :: CATCH_SNOW     (land_field,ntiles)
      real :: SNOW_GRND      (land_field,ntiles)

      ! SNSOOT may not be used as of vn6.6
      real :: SNSOOT(row_length,rows)  ! Snow soot content

      ! 2.9: aerosol emission fields, including mineral dust parent soil props

      real :: SOIL_CLAY (row_length,rows)
      real :: SOIL_SILT (row_length,rows)
      real :: SOIL_SAND (row_length,rows)
      real :: DUST_MREL1(row_length,rows)
      real :: DUST_MREL2(row_length,rows)
      real :: DUST_MREL3(row_length,rows)
      real :: DUST_MREL4(row_length,rows)
      real :: DUST_MREL5(row_length,rows)
      real :: DUST_MREL6(row_length,rows)

      real :: SO2_EM     (row_length,rows)
      real :: DMS_EM     (row_length,rows)
      real :: SO2_HILEM  (row_length,rows)
      real :: NH3_EM     (row_length,rows)
      real :: SOOT_EM    (row_length,rows)
      real :: SOOT_HILEM (row_length,rows)
      real :: BMASS_EM   (row_length,rows)
      real :: BMASS_HILEM(row_length,rows)
      real :: OCFF_EM    (row_length,rows)
      real :: OCFF_HILEM (row_length,rows)
      real :: DMS_CONC   (row_length,rows)
      real :: DMS_OFLUX  (row_length,rows)

! USER_ANC<N> fields are still 1D
      ! 2.10: User ancillary fields
      real :: USER_ANC1(row_length*rows*1)
      real :: USER_ANC2(row_length*rows*1)
      real :: USER_ANC3(row_length*rows*1)
      real :: USER_ANC4(row_length*rows*1)
      real :: USER_ANC5(row_length*rows*1)
      real :: USER_ANC6(row_length*rows*1)
      real :: USER_ANC7(row_length*rows*1)
      real :: USER_ANC8(row_length*rows*1)
      real :: USER_ANC9(row_length*rows*1)
      real :: USER_ANC10(row_length*rows*1)
      real :: USER_ANC11(row_length*rows*1)
      real :: USER_ANC12(row_length*rows*1)
      real :: USER_ANC13(row_length*rows*1)
      real :: USER_ANC14(row_length*rows*1)
      real :: USER_ANC15(row_length*rows*1)
      real :: USER_ANC16(row_length*rows*1)
      real :: USER_ANC17(row_length*rows*1)
      real :: USER_ANC18(row_length*rows*1)
      real :: USER_ANC19(row_length*rows*1)
      real :: USER_ANC20(row_length*rows*1)

      !   2.11: Store arrays for energy correction calculation
      real :: NET_FLUX (row_length,rows)
      real :: NET_MFLUX(row_length,rows)

      !   2.12: Tiled Vegetation and Triffid fields
      real :: FRAC_TYP (land_field)
      real :: FRAC_CON1(land_field)  ! fraction of broadleaf tree
      real :: FRAC_CON2(land_field)  ! fraction of needleleaf tree
      real :: FRAC_CON3(land_field)  ! fraction of C3 grass
      real :: FRAC_CON4(land_field)  ! fraction of C4 grass
      real :: FRAC_CON5(land_field)  ! fraction of shrub
      real :: FRAC_CON6(land_field)  ! fraction of urban
      real :: FRAC_CON7(land_field)  ! fraction of water
      real :: FRAC_CON8(land_field)  ! fraction of soil
      real :: FRAC_CON9(land_field)  ! fraction of ice
      real :: LAI_PFT  (land_field)
      real :: CANHT_PFT(land_field)
      real :: DISTURB_VEG(land_field)
      REAL :: soil_alb(land_field)   ! albedo of underlying soil
      REAL :: obs_alb_sw(land_field) ! Observed snow-free SW albedo 
      REAL :: obs_alb_vis(land_field)! Observed snow-free VIS albedo 
      REAL :: obs_alb_nir(land_field)! Observed snow-free NIR albedo 
      real, target :: SOIL_CARB(land_field)
      real, target :: SOIL_CARB1(land_field)
      real, target :: SOIL_CARB2(land_field)
      real, target :: SOIL_CARB3(land_field)
      real, target :: SOIL_CARB4(land_field)
      real :: NPP_PFT_ACC(land_field)
      real :: G_LF_PFT_ACC(land_field)
      real :: G_PHLF_PFT_ACC(land_field)
      real :: RSP_W_PFT_ACC(land_field)
      real, target :: RSP_S_ACC(land_field)
      real, target :: RSP_S_ACC1(land_field)
      real, target :: RSP_S_ACC2(land_field)
      real, target :: RSP_S_ACC3(land_field)
      real, target :: RSP_S_ACC4(land_field)
      real :: CAN_WATER_TILE(land_field,ntiles)
      real :: CATCH_TILE(land_field,ntiles)
      real :: INFIL_TILE(land_field,ntiles)
      real :: RGRAIN_TILE(land_field,ntiles)
      real :: SNODEP_TILE(land_field,ntiles)
      real :: TSTAR_TILE(land_field,ntiles) 
      real :: Z0_TILE(land_field,ntiles)  
      real :: z0h_tile(land_field,ntiles)  
      real :: DOLR_FIELD(row_length,rows) 
      real :: LW_DOWN(row_length,rows) 
      real :: SW_TILE_RTS(land_field,ntiles)

! MORUSES - new urban two-tile scheme
      real :: hgt(land_field)        ! Building height
      real :: hwr(land_field)        ! Height to width
      real :: wrr(land_field)        ! Width ratio
      real :: disp(land_field)       ! Displacement height
      real :: ztm(land_field)        ! Effective roughness length for momentum
      real :: albwl(land_field)      ! Wall albedo
      real :: albrd(land_field)      ! Road albedo
      real :: emisw(land_field)      ! Wall emmissivity
      real :: emisr(land_field)      ! Road emmissivity

!! REMOVED SLAB AS PART OF VN7.0
!      !   2.13: Slab Model

!   2.14: Carbon cycle fields
      real :: CO2FLUX(row_length,rows)
      real :: CO2_EMITS(row_length,rows)

!   2.15: Fields carried forward from previous version.
!         May not be required
!      real, pointer :: SURF_RESIST_NIT(:)  ! Surface resistance on
!                                    ! non-ice tiles
!      real, pointer :: ROOT_DEPTH_NIT(:)   ! Root depth on non-ice tiles
!      real, pointer :: Z0V_TYP(:)          ! Vegetative Roughness length on
!                                    ! tiles
!      real, pointer :: ICE_EDGE(:)
!      real, pointer :: OROG_TENDENCY(:)    ! Orographic tendencies
!      real, pointer :: OROG_SD_TENDENCY(:) ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
!      real, pointer :: ETATHETA(:)
!      real, pointer :: ETARHO(:)
!      real, pointer :: RHCRIT(:)
!      real, pointer :: SOIL_THICKNESS(:)

      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      real :: zseak_theta(0:model_levels)
      real :: Ck_theta(0:model_levels)
      real :: zseak_rho(model_levels)
      real :: Ck_rho(model_levels)   

      ! 2.16: Fields for large-scale hydrology scheme.
      real :: TI_MEAN(land_field)
      real :: TI_SIG(land_field)
      real :: FEXP(land_field)
      real :: GAMMA_INT(land_field)
      real :: WATER_TABLE(land_field)
      real :: FSFC_SAT(land_field)
      real :: F_WETLAND(land_field)

      real :: STHZW(land_field)
      real :: A_FSAT(land_field)
      real :: C_FSAT(land_field)
      real :: A_FWET(land_field)
      real :: C_FWET(land_field)

      ! 2.17: Fields for River routing.
      real :: RIV_SEQUENCE (river_row_length,river_rows)
      real :: RIV_DIRECTION(river_row_length,river_rows)
      real :: RIV_STORAGE  (river_row_length,river_rows)
      real :: TOT_SURFROFF (land_field)
      real :: TOT_SUBROFF  (land_field)
      real :: RIV_INLANDATM(land_field)
! Fields for water conservation correction due to lake evaporation 
      real :: ACC_LAKE_EVAP(row_length,rows) 
      ! Fields for grid-to-grid river routing (river routing 2A)
      real :: RIV_IAREA    (row_length,rows)   ! Drainage area
      real :: RIV_SLOPE    (row_length,rows)   ! Grid-cell slope
      real :: RIV_FLOWOBS1 (row_length,rows)   ! Initial values of flow
      real :: RIV_INEXT    (row_length,rows)   ! Flow direction (x)
      real :: RIV_JNEXT    (row_length,rows)   ! Flow direction (y)
      real :: RIV_LAND     (row_length,rows)   ! Land-type (land/river/sea)
      real :: RIV_SUBSTORE (row_length,rows)   ! Subsurface storage
      real :: RIV_SURFSTORE(row_length,rows)   ! Surface storage
      real :: RIV_FLOWIN   (row_length,rows)   ! Surface inflow
      real :: RIV_BFLOWIN  (row_length,rows)   ! Subsurface inflow

! Fields used when coupling using OASIS.
      real :: C_SOLAR   (row_length,rows)      ! CPL solar radn
      real :: C_BLUE    (row_length,rows)      ! CPL blue radn
      real :: C_LONGWAVE(row_length,rows)      ! CPL lw radn
      real :: C_TAUX    (row_length,rows)      ! CPL taux 
      real :: C_TAUY    (row_length,rows)      ! CPL tauy 
      real :: C_W10     (row_length,rows)      ! CPL 10m wind   
      real :: C_SENSIBLE(row_length,rows)      ! CPL sensible ht flx
      real :: C_SUBLIM  (row_length,rows,nice_use)! CPL sublim rate
      real :: C_EVAP    (row_length,rows)      ! CPL Evap rate
! Next field is the conductive flux through the ice that is used to
! force the ice model.  If using multilayers in ice (l_sice_multilayers=T), 
! this is the surface downwards heat flux (surf_ht_flux_sice), otherwise
! this is the conductive heat flux through all the ice (sea_ice_htf) 
! previously known as 'botmelt'.
      real :: C_FCONDTOPN(row_length,rows,nice)! CPL Multi-cat fcondtop
      real :: C_TOPMELTN(row_length,rows,nice) ! CPL Multi-cat tmlt
      real :: C_LSRAIN  (row_length,rows)      ! CPL Lg scl rain rate
      real :: C_LSSNOW  (row_length,rows)      ! CPL Lg scl snow rate
      real :: C_CVRAIN  (row_length,rows)      ! CPL Cnvctv rain rate
      real :: C_CVSNOW  (row_length,rows)      ! CPL Cnvctv snur rate
      real :: C_RIVEROUT(row_length,rows)      ! CPL Riv outflow->ocn
      real :: C_CALVING (row_length,rows)      ! CPL Iceberg Calving->ocn

!   2.18: JULES variables
      REAL :: SNOWDEPTH(    land_field,ntiles)
                             ! Snow depth on ground on tiles (m)
      REAL :: RHO_SNOW_GRND(land_field,ntiles)
                             ! Snowpack bulk density (kg/m3)
      REAL :: NSNOW(        land_field,ntiles)
                             ! Number of snow layers on ground on tiles
                             ! NOTE that this is converted to an integer.
      REAL :: DS(        land_field,ntiles,nsmax)
                             ! Snow layer thickness (m)
      REAL :: SICE(      land_field,ntiles,nsmax)
                             ! Snow layer ice mass on tiles (Kg/m2)
      REAL :: SLIQ(      land_field,ntiles,nsmax)
                             ! Snow layer liquid mass on tiles (Kg/m2)
      REAL :: TSNOWLAYER(land_field,ntiles,nsmax)
                             ! Snow layer temperature (K)
      REAL :: RHO_SNOW(  land_field,ntiles,nsmax)
                             ! Snow layer densities (kg/m3)
      REAL :: RGRAINL(   land_field,ntiles,nsmax)
                             ! Snow layer grain size on tiles (microns)

!         FLake lake scheme variables
      REAL :: lake_depth( land_field)
      REAL :: lake_fetch( land_field)
      REAL :: lake_t_mean(land_field)
      REAL :: lake_t_mxl( land_field)
      REAL :: lake_t_ice( land_field)
      REAL :: lake_h_mxl( land_field)
      REAL :: lake_h_ice( land_field)
      REAL :: lake_shape( land_field)
      REAL :: lake_g_dt(  land_field)

      ! UKCA oxidant fields
      real :: OH_UKCA ( tdims%i_start:tdims%i_end,                       &
                        tdims%j_start:tdims%j_end,                       &
                        tdims%k_start:tdims%k_end)
      real :: HO2_UKCA( tdims%i_start:tdims%i_end,                       &
                        tdims%j_start:tdims%j_end,                       &
                        tdims%k_start:tdims%k_end)
      real :: H2O2_UKCA(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: O3_UKCA  (tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: HNO3_UKCA(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)

      ! Aerosol climatologies
      real :: ARCLBIOG_BG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_IC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBLCK_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBLCK_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSSLT_FI(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSSLT_JT(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_AC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_AK(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_DI(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B1(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B2(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B3(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B4(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B5(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B6(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_IC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDLTA_DL(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      
      ! Ammonium nitrate aerosol
      real :: NITR_ACC (tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: NITR_DISS(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)

! End typ_atm_fields.h
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

      REAL                                                              & 
     &        t_inc_pres(row_length,rows,wet_levels)                    & 
! Temperature increment due to adiabatic temperature forcing (K) 
     &,       q_inc_pres(row_length,rows,wet_levels)                    & 
! Vapour increment due to adiabatic temperature forcing (kg kg-1) 
     &,       qcl_inc_pres(row_length,rows,wet_levels)                  & 
! Liquid increment due to adiabatic temperature forcing (kg kg-1) 
     &,       qcf_inc_pres(row_length,rows,wet_levels)                  & 
! Ice increment due to adiabatic temperature forcing (kg kg-1) 
     &,       cf_inc_pres(row_length,rows,wet_levels)                   & 
! Total cloud fraction increment due to adiabatic forcing (no units) 
     &,       cfl_inc_pres(row_length,rows,wet_levels)                  & 
! Liquid cloud fraction increment due to adiabatic forcing (no units) 
     &,       cff_inc_pres(row_length,rows,wet_levels)                  & 
! Ice cloud fraction increment due to adiabatic forcing (no units) 
     &,       t_dini(row_length,rows,wet_levels)                        & 
! Temperature increment due to initiation (K) 
     &,       q_dini(row_length,rows,wet_levels)                        & 
! Vapour increment due to initiation (kg kg-1) 
     &,       qcl_dini(row_length,rows,wet_levels)                      & 
! Liquid increment due to initiation (kg kg-1) 
     &,       qcf_dini(row_length,rows,wet_levels)                      & 
! Ice increment due to initiation (kg kg-1) 
     &,       cf_dini(row_length,rows,wet_levels)                       & 
! Total cloud fraction increment due to initiation (no units) 
     &,       cfl_dini(row_length,rows,wet_levels)                      & 
! Liquid cloud fraction increment due to initiation(no units) 
     &,       cff_dini(row_length,rows,wet_levels)
! Ice cloud fraction increment due to initiation (no units) 

      INTEGER                                                           &
     &        INT16                                                     &
                                ! STASHWORK size = STASH_MAXLEN(16)
     &       ,ICODE             ! Out return code : 0 Normal exit
!                               !                 :>0 Error exit

      CHARACTER(LEN=80)                                                    &
     &        CMESSAGE          ! Out error message if ICODE > 0

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

! C_ETA_PMSL start
! ETA_PMSL is the ETA value used to determine the model level
! used in PMSL reduction calculations.
      REAL,PARAMETER:: ETA_PMSL=0.795
      INTEGER LEVNO_PMSL_CALC   !  Model level for PMSL reduction calc.
      COMMON /PMSLCALC/ LEVNO_PMSL_CALC
! C_ETA_PMSL end

!L Locally dynamically allocated work area

      REAL                                                              &
     &        STASHWORK(INT16),                                         &
     &        T_P_PRESS(NUM_STASH_LEVELS),                              &
     &        HTS_PRESS(NUM_STASH_LEVELS),                              &
     &        RHice_press(NUM_STASH_LEVELS),                            &
     &        RHwat_press(NUM_STASH_LEVELS),                            &
     &        WBPT_PRESS(NUM_STASH_LEVELS),                             &
     &        TH_ADV_PRESS(NUM_STASH_LEVELS),                           &
     &        PRESS_LEVS(NUM_STASH_LEVELS),                             &
     &        TR_PRESS(TR_VARS+1,NUM_STASH_LEVELS),                     &
     &        T_m(THETA_FIELD_SIZE,MODEL_LEVELS),                       &
                                              ! T on model levels
     &        hts_theta(theta_field_size,model_levels),                 &
                                              ! hts on theta levels
     &        hts_rho(theta_field_size,model_levels),                   &
                                              ! hts on rho levels
     &        qc(theta_field_size,wet_levels),                          &
                                              ! cloud water content
     &        qT(theta_field_size,wet_levels)
                                              ! total specific humidity

! Local variables

      INTEGER                                                           &
     &        I,J,                                                      &
     &        NI,                                                       &
     &        K,                                                        &
     &        ISL,                                                      &
     &        BL,                                                       &
     &        TL,                                                       &
     &        LEVEL,                                                    &
     &        LAST_POINT,                                               &
     &        FIRST_POINT                                               &
     &       ,ITR                                                       &
     &       ,STASH_TR_FIRST                                            &
     &       ,STASH_TR_LAST                                             &
     &       ,im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index                                                  &
                            !  Internal Model Index for Stash arrays
     &      ,item                                                       &
                            ! STASH item
     &      ,sect           ! STASH section
      PARAMETER( sect = 16 ) ! for 'physics' end of timestep diagnostics

      INTEGER                                                           &
     &        T_P_LEVS,                                                 &
     &        HTS_LEVS, H2_P_LEVS,                                      &
     &        RHice_levs,RHwat_levs,                                    &
     &        WBPT_LEVS,                                                &
     &        TH_ADV_P_LEVS                                             &
     &       ,TR_PRESS_LEVS(TR_VARS+1)                                  &
     &       ,TR_THETA_FIELD_SIZE   ! size for DA of tracer fields

      INTEGER                                                           &
     &        H2_IND(NUM_STASH_LEVELS)

      INTEGER                                                           &
     &        PT_TRACER(TR_VARS+1)
      INTEGER                                                           &
     &        PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,PT209,    &
     &  PT210,PT211,PT212,PT213,PT214,PT215,PT216,PT217,PT218,PT219,    &
     &  PT220,PT221,PT222,PT223,PT224,PT225,PT255,PT256

      REAL                                                              &
     &        EW_SPACE,                                                 &
     &        NS_SPACE

      LOGICAL                                                           &
     &        SF_TRACER(TR_VARS+1)

      ! 3d work array for calculating separate +/- PC2 increments
      ! only to be allocated if it is needed.
      REAL, ALLOCATABLE ::  work3d(:,:,:)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!L Initialisation
      IF (lhook) CALL dr_hook('ST_DIAG2',zhook_in,zhook_handle)
        FIRST_POINT = 1
        LAST_POINT  = THETA_FIELD_SIZE

!L Internal Structure:

!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)

! Allocate work3 array if calculating +/- increments for cfl,cff,qcl,qcf.
      IF (sf(140,sect) .OR. sf(141,sect) .OR.                           &
          sf(142,sect) .OR. sf(143,sect) .OR.                           &
          sf(146,sect) .OR. sf(147,sect) .OR.                           &
          sf(148,sect) .OR. sf(149,sect) .OR.                           &
          sf(150,sect) .OR. sf(151,sect) .OR.                           &
          sf(156,sect) .OR. sf(157,sect) ) THEN
        ALLOCATE ( work3d(row_length,rows,wet_levels) )
      END IF

!L----- Calculate additional diagnostic quantities------------------

!L------------ Extract required pressures for T_P ----------------------
      NS_SPACE=A_REALHD(2)
      EW_SPACE=A_REALHD(1)

      ISL=STINDEX(1,203,16,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            T_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,T_P_LEVS
              T_P_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
        T_P_LEVS=1
      END IF

!L------------ Extract required pressures for Heights ------------------

      ISL=STINDEX(1,202,16,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            HTS_LEVS=STASH_LEVELS(1,NI)
            DO K =1,HTS_LEVS
              HTS_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
        HTS_LEVS=1
      END IF


!L------------ Extract required pressures for R H wrt ice --------------

      ISL=STINDEX(1,204,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            RHice_levs=STASH_LEVELS(1,NI)
            DO K =1,RHIce_levs
              RHice_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for RHice'
            ICODE=1
            IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for RHice'
          ICODE=1
          IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
          RETURN
        END IF
      ELSE
        RHice_Levs=1
      END IF

!L------------ Extract required pressures for R H wrt water ------------
      ISL=STINDEX(1,256,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            RHwat_levs=STASH_LEVELS(1,NI)
            DO K =1,RHwat_levs
              RHwat_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            CMESSAGE='ST_DIAG2 : Level not pressure for RHwat'
            ICODE=1
            IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
            RETURN
          END IF
        ELSE
          CMESSAGE='ST_DIAG2 : Level not a levels list for RHwat'
          ICODE=1
          IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
          RETURN
        END IF
      ELSE
        RHwat_levs=1
      END IF
!L------------ Extract required pressures for Wet bulb pot temp---------

      ISL=STINDEX(1,205,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            WBPT_LEVS=STASH_LEVELS(1,NI)
            DO K =1,WBPT_LEVS
              WBPT_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for WBPT'
            ICODE=1
            IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for WBPT'
          ICODE=1
          IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
          RETURN
        END IF
      ELSE
        WBPT_LEVS=1
      END IF
!L------------ Extract required pressures for Thermal advection --------

      ISL=STINDEX(1,219,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            TH_ADV_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,TH_ADV_P_LEVS
              TH_ADV_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for THADV'
            ICODE=1
            IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for THADV'
          ICODE=1
          IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
          RETURN
        END IF
      ELSE
        TH_ADV_P_LEVS=1
      END IF

!L------------ Extract required pressures for Tracers ------------------

      STASH_TR_FIRST=226
      STASH_TR_LAST=254
!     ITR is a count of tracers found to be using this diagnostic
      ITR=0
      IF (TR_VARS >  0) THEN

!     Initialize PT_TRACER and SF_TRACER arrays
        DO I=1,TR_VARS
          PT_TRACER(I)=1
          SF_TRACER(I)=.FALSE.
        END DO

      DO J=STASH_TR_FIRST,STASH_TR_LAST
        ISL=STINDEX(1,J,16,im_index)
        IF(ISL >  0) THEN
          IF(STLIST(10,ISL) <  0) THEN
            IF(STLIST(11,ISL) == 2) THEN
                ITR=ITR+1
              NI=-STLIST(10,ISL)
                TR_PRESS_LEVS(ITR)=STASH_LEVELS(1,NI)
                DO K=1,TR_PRESS_LEVS(ITR)
                  TR_PRESS(ITR,K)=STASH_LEVELS(K+1,NI)/1000.0
              ENDDO
                PT_TRACER(ITR)=SI(J,16,im_index)
                SF_TRACER(ITR)=SF(J,16)
            ELSE
               CMESSAGE='ST_DIAG2 : Level not pressure for Tracers'
              ICODE=1
              IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
              RETURN
            END IF
          ELSE
             CMESSAGE='ST_DIAG2 : Level not a levels list for Tracers'
            ICODE=1
            IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
            RETURN
          END IF
        END IF
      END DO
      END IF

!     Set last (or only) values in tracer pointer arrays
      PT_TRACER(TR_VARS+1)=1
      SF_TRACER(TR_VARS+1)=.FALSE.

!     Set size of TR_P_FIELD_DA depending on whether any tracers
!     are using this diagnostic. (Used for dynamic allocation of
!     tracer arrays in phy_diag).
      IF (ITR >  0) THEN
        TR_THETA_FIELD_SIZE=THETA_FIELD_SIZE
      ELSE
        TR_THETA_FIELD_SIZE=1

      END IF

!L----- check height available for calculation of height**2 -----------


      IF (SF(224,16)) THEN
        IF (.NOT.SF(202,16)) THEN
          CMESSAGE='ST_DIAG2 : ERROR h**2 requires H at same timestep'
          ICODE=1
          GOTO 9999
         ELSE
      ISL=STINDEX(1,224,16,im_index)
          IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            H2_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,H2_P_LEVS
              PRESS_LEVS(K)=STASH_LEVELS(K+1,NI)/1000.0
              DO I=1,HTS_LEVS
                IF(PRESS_LEVS(K) == HTS_PRESS(I)) THEN
                  H2_IND(k)=I
                ENDIF
              ENDDO
            ENDDO
          ELSE
            H2_P_LEVS=1
          END IF
         ENDIF
      ELSE
        H2_P_LEVS=1
      END IF


      PT201=SI(201,16,im_index)
      PT202=SI(202,16,im_index)
      PT203=SI(203,16,im_index)
      PT204=SI(204,16,im_index)
      PT205=SI(205,16,im_index)
      PT206=SI(206,16,im_index)
      PT207=SI(207,16,im_index)
      PT208=SI(208,16,im_index)
      PT209=SI(209,16,im_index)
      PT210=SI(210,16,im_index)
      PT211=SI(211,16,im_index)
      PT212=SI(212,16,im_index)
      PT213=SI(213,16,im_index)
      PT214=SI(214,16,im_index)
      PT215=SI(215,16,im_index)
      PT216=SI(216,16,im_index)
      PT217=SI(217,16,im_index)
      PT218=SI(218,16,im_index)
      PT219=SI(219,16,im_index)
      PT220=SI(220,16,im_index)
      PT221=SI(221,16,im_index)
      PT222=SI(222,16,im_index)
      PT223=SI(223,16,im_index)
      PT224=SI(224,16,im_index)
      PT225=SI(225,16,im_index)
      PT255=SI(255,16,im_index)
      PT256=SI(256,16,im_index)

! Initialise STASHWORK as PHY_DIAG avoids MPP halo calculations
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT16
        STASHWORK(I)=0.
      ENDDO

      CALL Phy_diag(                                                    &
! Primary data: in
     & PSTAR,P,RHO,U,V,W                                                &
     &,THETA,Q,QCL,QCF                                                  &
     &,p_theta_levels                                                   &
     &,exner_rho_levels,exner_theta_levels                              &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
     &,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
     &,Model_domain,npmsl_height,L_pmsl_sor                             &
! Pressure levels for output arrays: in
     &,HTS_PRESS,T_P_PRESS                                              &
     &,RHice_press,RHwat_press,WBPT_PRESS                               &
! Flags to request each diagnostic output field: in
     &,SF(  4,16),SF(201,16),SF(202,16),SF(203,16),SF(204,16)           &
     &,SF(205,16),SF(206,16),SF(207,16),SF(222,16),SF(255,16)           &
     &,SF(256,16)                                                       &
! Diagnostics lengths: in
     &,HTS_LEVS,T_P_LEVS,RHice_levs,RHwat_levs,WBPT_LEVS                &
! Diagnostic arrays: out
     &,T_m                                                              &
     &,hts_theta,STASHWORK(PT202),STASHWORK(PT203)                      &
     &,STASHWORK(PT204),STASHWORK(PT205),qc,qT,STASHWORK(PT222)         &
     &,hts_rho,STASHWORK(PT256)                                         &
     &)



      item = 4      !  (4,16) is T on model levels
      IF( SF(item,sect) ) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),T_m,         &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

! PC2 cloud scheme diagnostics for condensation due to vertical             
! advection                                                                 
                                                                            
      item = 181      !  (181,16)  Temperature change from condensation     
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d             
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),t_inc_pres,  &      
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                                            
      item = 182      !  (182,16)  Vapour change from condensation          
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),q_inc_pres,  &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                                            
      item = 183      !  (183,16)  Liquid change from condensation 
      IF( SF(item,sect) ) THEN 
! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qcl_inc_pres,&    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)
                                                                            
      item = 150      !  (150,16)  Liquid change condensation positive
      IF( SF(item,sect) ) THEN 
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,qcl_inc_pres(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)
                                                                            
      item = 151      !  (151,16)  Liquid change condensation negative
      IF( SF(item,sect) ) THEN 
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,qcl_inc_pres(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)

      item = 184      !  (184,16)  Ice change from condensation             
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                               
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qcf_inc_pres,&    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             

      item = 192      !  (192,16)  Total cloud change from condensation     
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                          
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cf_inc_pres, &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                                            
      item = 193      !  (193,16)  Liquid cloud change from condensation    
      IF( SF(item,sect) ) THEN 
! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cfl_inc_pres,&    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect) 
 
      item = 156      !  (156,16)  Liquid cloud change cond'n positive
      IF( SF(item,sect) ) THEN 
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,cfl_inc_pres(i,j,k))
            END DO
          END DO
        END DO

! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect) 
 
      item = 157      !  (157,16)  Liquid cloud change cond'n positive
      IF( SF(item,sect) ) THEN 
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,cfl_inc_pres(i,j,k))
            END DO
          END DO
        END DO

! DEPENDS ON: copydiag_3d                                                
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)  

      item = 194      !  (194,16)  Ice cloud change from condensation       
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                                     
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cff_inc_pres,&    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                                            
! PC2 cloud scheme diagnostics for condensation due to initialisation 
                                                                            
      item = 161      !  (161,16)  Temperature change from initialisation     
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                              
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),t_dini,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                                            
      item = 162      !  (162,16)  Vapour change from initialisation         
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                          
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),q_dini,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
 
      item = 163      !  (163,16)  Liquid change from initialisation        
      IF( SF(item,sect) ) THEN   
! DEPENDS ON: copydiag_3d                                                   
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qcl_dini,    &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
 
      item = 140      !  (140,16)  Liquid change initialisation positive       
      IF( SF(item,sect) ) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,qcl_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                   
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)  
 
      item = 141      !  (141,16)  Liquid change initialisation negative       
      IF( SF(item,sect) ) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,qcl_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                   
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect) 

      item = 164      !  (164,16)  Ice change from initialisation            
      IF( SF(item,sect) ) THEN  
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qcf_dini,    &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               

      ENDIF   ! SF(  item,sect)      

      item = 142      !  (142,16)  Ice change initialisation positive
      IF( SF(item,sect) ) THEN  
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,qcf_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               

      ENDIF   ! SF(  item,sect)

      item = 143      !  (143,16)  Ice change initialisation negative
      IF( SF(item,sect) ) THEN  
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,qcf_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               

      ENDIF   ! SF(  item,sect)
 
      item = 172      !  (172,16)  Total cloud change from initialisation      
      IF( SF(item,sect) ) THEN                                              
! DEPENDS ON: copydiag_3d                                                 
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cf_dini,     &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect)                                             
                                                         
      item = 173      ! (173,16)  Liquid cloud change from initialisation    
      IF( SF(item,sect) ) THEN  
! DEPENDS ON: copydiag_3d                                                 
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cfl_dini,    &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               
                                                                            
      ENDIF   ! SF(  item,sect) 

      item = 146      !  (146,16)  Liquid cloud change init'n positive
      IF( SF(item,sect) ) THEN  
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,cfl_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               

      ENDIF   ! SF(  item,sect)

      item = 147      !  (147,16)  Liquid cloud change init'n negative
      IF( SF(item,sect) ) THEN  
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,cfl_dini(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d                                                  
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      &    
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &    
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &    
     &        stash_levels,num_stash_levels+1,                          &    
     &        im_ident,sect,item,                                       &    
     &        icode,cmessage)                                               

      ENDIF   ! SF(  item,sect)

      item = 174      !  (174,16)  Ice cloud change from initialisation      
      IF( SF(item,sect) ) THEN 
! DEPENDS ON: copydiag_3d                                                 
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),cff_dini,    & 
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       & 
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       & 
     &        stash_levels,num_stash_levels+1,                          & 
     &        im_ident,sect,item,                                       & 
     &        icode,cmessage) 
 
      ENDIF   ! SF(  item,sect)   

      item = 148     !  (148,16)  Ice cloud change init'n positive     
      IF( SF(item,sect) ) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,cff_dini(i,j,k))
            END DO
          END DO
        END DO

! DEPENDS ON: copydiag_3d                                                 
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      & 
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       & 
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       & 
     &        stash_levels,num_stash_levels+1,                          & 
     &        im_ident,sect,item,                                       & 
     &        icode,cmessage) 
 
      ENDIF   ! SF(  item,sect) 

      item = 149     !  (149,16)  Ice cloud change init'n negative     
      IF( SF(item,sect) ) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,cff_dini(i,j,k))
            END DO
          END DO
        END DO

! DEPENDS ON: copydiag_3d                                                 
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),work3d,      & 
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       & 
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       & 
     &        stash_levels,num_stash_levels+1,                          & 
     &        im_ident,sect,item,                                       & 
     &        icode,cmessage) 
 
      ENDIF   ! SF(  item,sect) 

      item = 201    !  (201,16) is geopotential height on theta levels

      IF( SF(item,sect) ) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),hts_theta,   &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      item = 206    !  (206,16) Cloud water content (qc) on q levels

      IF( SF(item,sect) ) THEN
! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qc,          &
     &        row_length,rows,wet_levels,0,0,0,0, at_extremity,         &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      item = 207    !  (207,16) Total specific humidity (qT) on q levels

      IF( SF(item,sect) ) THEN
! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),qT,          &
     &        row_length,rows,wet_levels,0,0,0,0, at_extremity,         &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      item = 255    !  (255,16) is geopotential height on rho levels

      IF( SF(item,sect) ) THEN
! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),hts_rho,     &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      IF (sf(140,sect) .OR. sf(141,sect) .OR.                           &
          sf(142,sect) .OR. sf(143,sect) .OR.                           &
          sf(146,sect) .OR. sf(147,sect) .OR.                           &
          sf(148,sect) .OR. sf(149,sect) .OR.                           &
          sf(150,sect) .OR. sf(151,sect) .OR.                           &
          sf(156,sect) .OR. sf(157,sect) ) THEN
        DEALLOCATE ( work3d )
      END IF


! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,16,STASHWORK,                                &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &           ICODE,CMESSAGE)

 9999 CONTINUE
      IF (lhook) CALL dr_hook('ST_DIAG2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ST_DIAG2
