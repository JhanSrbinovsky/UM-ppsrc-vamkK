! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE write_sulpc_oxidants_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE write_sulpc_oxidants(                                               &
  oh_conc_in, h2o2_mmr_in, ho2_conc_in, o3_mmr_in, hno3_mmr_in,                &
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
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
  first_atmstep_call                                                           &
  )

!---------------------------------------------------------------------
! Purpose: To write to the D1 array the concentrations of oxidants
!          after depletion in the sulphur cycle.
!          Called by Atm_Step only if L_SULPC_ONLINE_OXIDANTS is set
!          to .TRUE.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v8 programming standards
!
! Documentation:
!
!---------------------------------------------------------------------

!     needed by typ_atm_fields.h:
USE atm_fields_bounds_mod

USE ancil_info, ONLY: nsmax
USE ukca_d1_defs
USE UKCA_CONSTANTS, ONLY: c_oh, c_ho2

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE lbc_mod
USE Submodel_Mod
IMPLICIT NONE

!
! Global variables:
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
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
! C_RMOL start
!
! Description:
! This contains the molar universal gas constant (8.314 J K-1 MOL-1)
!
!
      REAL,PARAMETER:: RMOL = 8.314 ! molar gas const

! C_RMOL end

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

! Variables with intent IN:
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
REAL, INTENT(in) :: o3_mmr_in(row_length,rows,model_levels)
REAL, INTENT(in) :: h2o2_mmr_in(row_length,rows,model_levels)
REAL, INTENT(in) :: oh_conc_in(row_length,rows,model_levels)
REAL, INTENT(in) :: ho2_conc_in(row_length,rows,model_levels)
REAL, INTENT(in) :: hno3_mmr_in(row_length,rows,model_levels)
LOGICAL, INTENT(in) :: first_atmstep_call

! Local arrays:
REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_theta_levels_local
REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_theta_levels_local
REAL, DIMENSION(:,:,:), ALLOCATABLE :: theta_local
REAL, DIMENSION(:,:,:), ALLOCATABLE :: t
REAL, DIMENSION(:,:,:), ALLOCATABLE :: air_density

REAL :: oh_mmr_no_halo(row_length,rows,model_levels)
!                     OH mass mixing ratio (kg[OH]/kg[air])
REAL  :: ho2_mmr_no_halo(row_length,rows,model_levels)
!                    HO2 mass mixing ratio (kg[HO2]/kg[air])
REAL, DIMENSION(:,:,:), ALLOCATABLE :: o3_mmr
!                     O3 mass mixing ratio (kg[O3]/kg[air])
REAL, DIMENSION(:,:,:), ALLOCATABLE :: h2o2_mmr
!                     H2O2 mass mixing ratio (kg[H2O2]/kg[air])
REAL, DIMENSION(:,:,:), ALLOCATABLE :: oh_mmr
!                     OH mass mixing ratio (kg[OH]/kg[air])
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ho2_mmr
!                     HO2 mass mixing ratio (kg[HO2]/kg[air])
REAL, DIMENSION(:,:,:), ALLOCATABLE :: hno3_mmr
!                     HNO3 mass mixing ratio (kg[HNO3]/kg[air])

! Arrays with intent OUT:
! Local scalars:
INTEGER :: m_atm_modl ! Atmosphere model identifier
INTEGER :: section    ! Section number of diagnostic
INTEGER :: item       ! Item number of diagnostic
INTEGER :: halo_typ   ! Halo type of diagnostic
INTEGER :: error_code ! For passing to subroutine EREPORT
INTEGER, SAVE :: x_halo_oh     ! } Haloes
INTEGER, SAVE :: x_halo_o3     ! }
INTEGER, SAVE :: x_halo_h2o2   ! } on UKCA
INTEGER, SAVE :: x_halo_ho2    ! }
INTEGER, SAVE :: x_halo_hno3   ! }
INTEGER, SAVE :: y_halo_oh     ! } oxidants
INTEGER, SAVE :: y_halo_o3     ! }
INTEGER, SAVE :: y_halo_h2o2   ! } in D1
INTEGER, SAVE :: y_halo_ho2    ! }
INTEGER, SAVE :: y_halo_hno3   ! }
INTEGER :: i          ! Loop variable
INTEGER :: j          ! Loop variable
INTEGER :: k          ! Loop variable

CHARACTER(LEN=72) :: cmessage        ! Error message

! Denotes whether or not each UKCA oxidant has been found
! in D1 array
LOGICAL, DIMENSION(5) :: l_ukca_found 

IF (lhook) CALL dr_hook('WRITE_SULPC_OXIDANTS',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Determine halo sizes for UKCA oxidant arrays:

IF (first_atmstep_call) THEN
  l_ukca_found(:) = .FALSE.             ! Initalise to F
  m_atm_modl = submodel_for_sm(a_im)
  DO i=1,no_obj_d1(m_atm_modl)
    section = d1_addr(d1_section,i,m_atm_modl)
    IF (section == ukca_sect) THEN
      item = d1_addr(d1_item,i,m_atm_modl)
      IF (item == ukca_item_sulpc(1)) THEN
        halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
        x_halo_o3 = halosize(1,halo_typ)
        y_halo_o3 = halosize(2,halo_typ)
        l_ukca_found(1) = .TRUE.
      ELSE IF (item == ukca_item_sulpc(2)) THEN
        halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
        x_halo_hno3 = halosize(1,halo_typ)
        y_halo_hno3 = halosize(2,halo_typ)
        l_ukca_found(2) = .TRUE.
      ELSE IF (item == ukca_item_sulpc(3)) THEN
        halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
        x_halo_h2o2 = halosize(1,halo_typ)
        y_halo_h2o2 = halosize(2,halo_typ)
        l_ukca_found(3) = .TRUE.
      ELSE IF (item ==  ukca_item_sulpc(4)) THEN
        halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
        x_halo_oh = halosize(1,halo_typ)
        y_halo_oh = halosize(2,halo_typ)
        l_ukca_found(4) = .TRUE.
      ELSE IF (item ==  ukca_item_sulpc(5)) THEN
        halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
        x_halo_ho2 = halosize(1,halo_typ)
        y_halo_ho2 = halosize(2,halo_typ)
        l_ukca_found(5) = .TRUE.
      END IF ! IF condition on item
    END IF   ! IF condition on section
  END DO

  DO i=1,5
    IF (.NOT. l_ukca_found(i)) THEN
      ! If any UKCA oxidants cannot be found in D1 array, return an error
      ! message and stop.
      cmessage = 'At least one UKCA oxidant not found in D1 '                  &
        &                     //'array.'
      error_code = (ukca_sect * 1000) + ukca_item_sulpc(i)
      CALL ereport('WRITE_SULPC_OXIDANTS',error_code,cmessage)
    END IF
  END DO
END IF ! IF(first_atmstep_call)

! Allocate array dimensions:

ALLOCATE(o3_mmr(1-x_halo_o3:row_length+x_halo_o3,                              &
  1-y_halo_o3:rows+y_halo_o3, 1:model_levels))
ALLOCATE(hno3_mmr(1-x_halo_hno3:row_length+x_halo_hno3,                        &
  1-y_halo_hno3:rows+y_halo_hno3, 1:model_levels))
ALLOCATE(h2o2_mmr(1-x_halo_h2o2:row_length+x_halo_h2o2,                        &
  1-y_halo_h2o2:rows+y_halo_h2o2, 1:model_levels))
ALLOCATE(oh_mmr(1-x_halo_oh:row_length+x_halo_oh,                              &
  1-y_halo_oh:rows+y_halo_oh, 1:model_levels))
ALLOCATE(ho2_mmr(1-x_halo_ho2:row_length+x_halo_ho2,                           &
  1-y_halo_ho2:rows+y_halo_ho2, 1:model_levels))

! Initialise arrays to zero:
o3_mmr(:,:,:) = 0.0
hno3_mmr(:,:,:) = 0.0
h2o2_mmr(:,:,:) = 0.0
oh_mmr(:,:,:) = 0.0
ho2_mmr(:,:,:) = 0.0

! Extract from D1 the arrays necessary for calculation of air density:
! (air density is needed for conversion of molecules/cm3 to MMR)
ALLOCATE(p_theta_levels_local(1-offx:row_length+offx,                          &
  1-offy:rows+offy, model_levels))
p_theta_levels_local = RESHAPE(p_theta_levels,                                 &
  (/ row_length+(2*offx),rows+(2*offy),                                        &
  model_levels /) )
ALLOCATE(theta_local(1-offx:row_length+offx,1-offy:rows+offy,                  &
  model_levels))
theta_local = RESHAPE(theta,                                                   &
  (/row_length+(2*offx), rows+(2*offy), model_levels/))
ALLOCATE(exner_theta_levels_local(1-offx:row_length+offx,                      &
  1-offy:rows+offy, model_levels))
exner_theta_levels_local = RESHAPE(exner_theta_levels,                         &
  (/ row_length+(2*offx), rows+(2*offy), model_levels /) )
ALLOCATE(t(row_length,rows,model_levels))
ALLOCATE(air_density(row_length,rows,model_levels))

! Convert OH and HO2 from molecules/cm3 to mass-mixing ratios, to
! be written to D1 array:
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      t(i,j,k)=theta_local(i,j,k)*exner_theta_levels_local(i,j,k)
      air_density(i,j,k) = p_theta_levels_local(i,j,k) * avogadro              &
        * 1.0e-06 / (rmol * t(i,j,k))
      oh_mmr_no_halo(i,j,k) = oh_conc_in(i,j,k) * c_oh                         &
        / air_density(i,j,k)
      ho2_mmr_no_halo(i,j,k) = ho2_conc_in(i,j,k) * c_ho2                      &
        / air_density(i,j,k)
    END DO
  END DO
END DO

! Deallocate arrays that are no longer required:
DEALLOCATE(p_theta_levels_local)
DEALLOCATE(theta_local)
DEALLOCATE(exner_theta_levels_local)
DEALLOCATE(t)
DEALLOCATE(air_density)

! Arrays to be written to D1:
o3_mmr(1:row_length,1:rows,1:model_levels)                                     &
  = o3_mmr_in(1:row_length,1:rows,1:model_levels)
hno3_mmr(1:row_length,1:rows,1:model_levels)                                   &
  = hno3_mmr_in(1:row_length,1:rows,1:model_levels)
h2o2_mmr(1:row_length,1:rows,1:model_levels)                                   &
  = h2o2_mmr_in(1:row_length,1:rows,1:model_levels)
oh_mmr(1:row_length,1:rows,1:model_levels)                                     &
  = oh_mmr_no_halo(1:row_length,1:rows,1:model_levels)
ho2_mmr(1:row_length,1:rows,1:model_levels)                                    &
  = ho2_mmr_no_halo(1:row_length,1:rows,1:model_levels)

! Write to D1 array:
o3_ukca(1-offx:row_length+offx, 1-offy:rows+offy,                              &
  1:model_levels) =                                                            &
  o3_mmr(1-x_halo_o3:row_length+x_halo_o3,                                     &
  1-y_halo_o3:rows+y_halo_o3, 1:model_levels)

hno3_ukca(1-offx:row_length+offx, 1-offy:rows+offy,                            &
  1:model_levels) =                                                            &
  hno3_mmr(1-x_halo_hno3:row_length+x_halo_hno3,                               &
  1-y_halo_hno3:rows+y_halo_hno3, 1:model_levels)

h2o2_ukca(1-offx:row_length+offx, 1-offy:rows+offy,                            &
  1:model_levels) =                                                            &
  h2o2_mmr(1-x_halo_h2o2:row_length+x_halo_h2o2,                               &
  1-y_halo_h2o2:rows+y_halo_h2o2, 1:model_levels)

oh_ukca(1:row_length, 1:rows, 1:model_levels) =                                &
  oh_mmr(1-x_halo_oh:row_length+x_halo_oh,                                     &
  1-y_halo_oh:rows+y_halo_oh, 1:model_levels)

ho2_ukca(1:row_length, 1:rows, 1:model_levels) =                               &
  ho2_mmr(1-x_halo_ho2:row_length+x_halo_ho2,                                  &
  1-y_halo_ho2:rows+y_halo_ho2, 1:model_levels)

! Deallocate temporary arrays
DEALLOCATE(o3_mmr)
DEALLOCATE(hno3_mmr)
DEALLOCATE(h2o2_mmr)
DEALLOCATE(oh_mmr)
DEALLOCATE(ho2_mmr)

IF (lhook) CALL dr_hook('WRITE_SULPC_OXIDANTS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE write_sulpc_oxidants
END MODULE write_sulpc_oxidants_mod
