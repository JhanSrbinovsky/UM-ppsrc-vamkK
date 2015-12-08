! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ End of Timestep climate diagnostics control routine
!
! Subroutine Interface:
      SUBROUTINE ST_DIAG3( STASHWORK,INT30,                             &
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
     &    energy_corr_now,                                              &
     &    inc_u, inc_v, inc_w, inc_t,                                   &
     &    inc_q, inc_qcl, inc_qcf, inc_rho, sin_v_latitude,             &
          inc_qrain, inc_qgraup, inc_qcf2,                              &
          wet_to_dry_n,                                                 &
     &                    ICODE,CMESSAGE)

      Use level_heights_mod
      Use trignometric_mod,  Only : sin_theta_latitude,                 &
     &     cos_theta_latitude, sin_theta_longitude, cos_theta_longitude
      Use rad_mask_trop_mod, Only : max_trop_level,min_trop_level
      Use ancil_info, only:                                             &
     &                      nsmax

      USE atm_fields_bounds_mod, ONLY:                                  &  
        o3dims,                                                         & 
        oneddims,                                                       & 
        pdims_s,                                                        &  
        qdims,   qdims_s, qdims_l,                                      & 
        rdims2,                                                         & 
        tdims,   tdims_s, tdims_l,                                      &  
        udims_s, udims_l,                                               &  
        vdims_s, vdims_l,                                               &  
        wdims,   wdims_s, wdims_l              

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      USE check_prod_levs_mod, ONLY: check_prod_levs
      USE eot_diag_mod, ONLY: eot_diag
      USE um_input_control_mod,  ONLY: model_domain
      USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
      USE Submodel_Mod

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE


! Description:
! Calculates end of timestep diagnostics (held in STASH section 30)
! of which there are three types:
! 1) Model level fields and products
! 2) Pressure level fields and products ie U V T UU VV TT UT VT
! 3) Single level fields eg mountain torque and column integrations
!    eq atmospheric mass
!
! Method:
!   This routine sets controls the EoT diagnostics, and sets up
!   array sizes etc for a call to EOT_DIAG which calculates the
!   fields. STASHWORK is filled with the required data, and
!   passed up to atm_step2, ready for a call to stash.
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!  Standard diagnostics are assigned a number, which are:
!  1 U wind, 2 V wind, 3 W wind, 4 Temperature, 5 Q, 6 RH, 7 Z, 8 Omega
!
!
!  Further user and technical documentation can be found under:
!  http://www-hc/~hadsu/
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: stash
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Global variables (*CALLed COMDECKs etc...):

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
! C_ETA_PMSL start
! ETA_PMSL is the ETA value used to determine the model level
! used in PMSL reduction calculations.
      REAL,PARAMETER:: ETA_PMSL=0.795
      INTEGER LEVNO_PMSL_CALC   !  Model level for PMSL reduction calc.
      COMMON /PMSLCALC/ LEVNO_PMSL_CALC
! C_ETA_PMSL end

      Integer, parameter ::                                             &
     &  npress_diags=8                                                  &
                                 ! No. of diags on pressure levels
     &, nmodel_diags=7           ! No. of diags on modellevels

!   Scalar  arguments with intent(in):

!   Array  arguments with intent(in):
      real energy_corr_now                                              &
      , inc_u(udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end)                            &
      , inc_v(vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
      , inc_w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,      &
              wdims%k_start:wdims%k_end)                                &
      , inc_t(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
      , inc_q(qdims_l%i_start:qdims_l%i_end,                            &
              qdims_l%j_start:qdims_l%j_end,                            &
              1:qdims_l%k_end)                                          &
      , inc_qcl(qdims_l%i_start:qdims_l%i_end,                          &
                qdims_l%j_start:qdims_l%j_end,                          &
                1:qdims_l%k_end)                                        &
      , inc_qcf(qdims_l%i_start:qdims_l%i_end,                          &
                qdims_l%j_start:qdims_l%j_end,                          &
                1:qdims_l%k_end)                                        &
      , inc_qrain(qdims_l%i_start:qdims_l%i_end,                        &
                  qdims_l%j_start:qdims_l%j_end,                        &
                  1:qdims_l%k_end)                                      &
      , inc_qgraup(qdims_l%i_start:qdims_l%i_end,                       &
                   qdims_l%j_start:qdims_l%j_end,                       &
                   1:qdims_l%k_end)                                     &
      , inc_qcf2(qdims_l%i_start:qdims_l%i_end,                         &
                 qdims_l%j_start:qdims_l%j_end,                         &
                 1:qdims_l%k_end)                                       &
      , inc_rho(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)                          &
      , wet_to_dry_n(qdims_s%i_start:qdims_s%i_end,                     &
                     qdims_s%j_start:qdims_s%j_end,                     &
                     qdims_s%k_start:qdims_s%k_end)      
      

!   Scalar arguments with intent(InOut):
      INTEGER                                                           &
     &        INT30,                                                    &
                                ! Dummy for STASH_MAXLEN(30)
     &        ICODE              ! Out return code : 0 Normal exit
                                !                 : >0 Error exit

      CHARACTER(LEN=256)                                                     &
     &        CMESSAGE          ! Out error message if ICODE > 0

!   Array  arguments with intent(InOut):
      REAL                                                              &
     &  STASHWORK(*)

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     &  U_M_LEVS,                                                       &
     &  V_M_LEVS,                                                       &
     &  W_M_LEVS,                                                       &
     &  T_M_LEVS,                                                       &
     &  Q_M_LEVS,                                                       &
     &  Z_M_LEVS,                                                       &
     &  KE_M_LEVS,                                                      &
     &  OM_M_LEVS,                                                      &
     &  wbig_m_levs,                                                    &
     &  wbig2_m_levs,                                                   &
     &  dry_mass_m_levs,                                                &
     &  Rh_m_levs,                                                      &
     &  U_MM_LEVS,                                                      &
     &  V_MM_LEVS,                                                      &
     &  W_MM_LEVS,                                                      &
     &  T_MM_LEVS,                                                      &
     &  Q_MM_LEVS,                                                      &
     &  Z_MM_LEVS,                                                      &
     &  KE_MM_LEVS,                                                     &
     &  TEOT_M_LEVS,                                                    &
     &  U_INC_LEVS,                                                     &
     &  V_INC_LEVS,                                                     &
     &  W_INC_LEVS,                                                     &
     &  T_INC_LEVS,                                                     &
     &  Q_INC_LEVS,                                                     &
     &  QCL_INC_LEVS,                                                   &
     &  QCF_INC_LEVS,                                                   &
        qrain_inc_levs,                                                 &
        qgraup_inc_levs,                                                &
        qqcf2_inc_levs,                                                 &
     &  RHO_INC_LEVS,                                                   &
     &  U2_INC_LEVS,                                                    &
     &  V2_INC_LEVS,                                                    &
     &  W2_INC_LEVS,                                                    &
     &  T2_INC_LEVS,                                                    &
     &  Q2_INC_LEVS,                                                    &
     &  QCL2_INC_LEVS,                                                  &
     &  QCF2_INC_LEVS,                                                  &
     &  RHO2_INC_LEVS,                                                  &
     &  U_P_LEVS,                                                       &
     &  V_P_LEVS,                                                       &
     &  W_P_LEVS,                                                       &
     &  T_P_LEVS,                                                       &
     &  Q_P_LEVS,                                                       &
     &  RH_P_LEVS,                                                      &
     &  Z_P_LEVS,                                                       &
     &  OM_P_LEVS,                                                      &
     &  HEAVY_P_LEVS,                                                   &
     &  TV_P_LEVS,                                                      &
     &  TVOM_P_LEVS,                                                    &
     &  VSTARBAR_P_LEVS,                                                &
     &  WSTARBAR_P_LEVS,                                                &
     &  FY_P_LEVS,                                                      &
     &  FZ_P_LEVS,                                                      &
     &  DIVF_P_LEVS,                                                    &
     &  HEATFLUX_P_LEVS,                                                &
     &  MOMFLUX_P_LEVS,                                                 &
     &  FIELD1_P_LEVS,                                                  &
     &  FIELD2_P_LEVS,                                                  &
     &  n_levels                                                        &
     &  ,im_ident                                                       &
                                !  Internal Model Identifier
     &  ,im_index               !  Internal Model Index for Stash Arrays

      INTEGER                                                           &
     &  PT001,PT002,PT003,PT004,PT005,PT006,PT007,PT008,                &
     &  PT101,PT102,PT103,PT104,PT105,PT106,PT107,                      &
     &  PT111,PT112,PT113,PT114,PT115,                                  &
     &  PT171,PT172,PT173,PT174,PT175,PT176,PT177,PT178,                &
     &  PT181,PT182,PT183,PT184,PT185,PT186,PT187,PT188,pt189,pt190,    &
        pt191,                                                          &
     &  PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,                &
     &  PT301,PT302,PT303,                                              &
     &  PT310,PT311,PT312,PT313,PT314,PT315,PT316,                      &
     &  PT401,PT402,PT403,PT404,PT405,PT406,PT407,PT408,PT409,PT410,    &
     &  PT411,PT412,PT413,PT414,PT415,PT416,PT417,PT418,PT419,PT420,    &
     &  PT421,PT422,PT423,PT424,PT425,PT426,PT427,PT428,PT429,PT430,    &
     &  PT431,PT432,PT433,PT434,PT435,PT436,PT437,PT438,PT439           &
     & ,PT440,PT441,pt442,PT451,PT452,PT453,PT454

!   Array  arguments with intent(out):
      REAL                                                              &
     &  U_PRESS(NUM_STASH_LEVELS)                                       &
     &  ,V_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,W_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,T_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,Q_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,RH_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,Z_PRESS(NUM_STASH_LEVELS)                                      &
                                   ! Z pressure levels
     &  ,OM_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,HEAVY_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,PRESS_LEVS(NUM_STASH_LEVELS)                                   &
     &  ,VSTARBAR_PRESS(NUM_STASH_LEVELS)                               &
        ,WSTARBAR_PRESS(NUM_STASH_LEVELS)                               &
        ,FY_PRESS(NUM_STASH_LEVELS)                                     &
        ,FZ_PRESS(NUM_STASH_LEVELS)                                     &
        ,DIVF_PRESS(NUM_STASH_LEVELS)                                   &
        ,HEATFLUX_PRESS(NUM_STASH_LEVELS)                               &
        ,MOMFLUX_PRESS(NUM_STASH_LEVELS)

      LOGICAL U_M_LIST(MODEL_LEVELS)                                    &
     &  ,V_M_LIST(MODEL_LEVELS)                                         &
     &  ,W_M_LIST(MODEL_LEVELS)                                         &
     &  ,T_M_LIST(MODEL_LEVELS)                                         &
     &  ,Q_M_LIST(MODEL_LEVELS)                                         &
     &  ,Z_M_LIST(MODEL_LEVELS)                                         &
     &  ,KE_M_LIST(MODEL_LEVELS)                                        &
     &  ,OM_M_LIST(MODEL_LEVELS)                                        &
     &  ,wbig_m_list(MODEL_LEVELS)                                      &
     &  ,wbig2_m_list(MODEL_LEVELS)                                     &
     &  ,RH_m_list(MODEL_LEVELS)                                        &
     &  ,DRY_MASS_M_LIST(MODEL_LEVELS)                                  &
     &  ,U_MM_LIST(MODEL_LEVELS)                                        &
     &  ,V_MM_LIST(MODEL_LEVELS)                                        &
     &  ,W_MM_LIST(MODEL_LEVELS)                                        &
     &  ,T_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Q_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Z_MM_LIST(MODEL_LEVELS)                                        &
     &  ,KE_MM_LIST(MODEL_LEVELS)                                       &
     &  ,TEOT_M_LIST(MODEL_LEVELS)                                      &
     &  ,U_INC_LIST(MODEL_LEVELS)                                       &
     &  ,V_INC_LIST(MODEL_LEVELS)                                       &
     &  ,W_INC_LIST(MODEL_LEVELS)                                       &
     &  ,T_INC_LIST(MODEL_LEVELS)                                       &
     &  ,Q_INC_LIST(MODEL_LEVELS)                                       &
     &  ,QCL_INC_LIST(MODEL_LEVELS)                                     &
     &  ,QCF_INC_LIST(MODEL_LEVELS)                                     &
        ,qrain_inc_list(model_levels)                                   &
        ,qgraup_inc_list(model_levels)                                  &
        ,qqcf2_inc_list(model_levels)                                   &
     &  ,RHO_INC_LIST(MODEL_LEVELS)                                     &
     &  ,U2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,V2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,W2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,T2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,Q2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,QCL2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,QCF2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,RHO2_INC_LIST(MODEL_LEVELS)

      logical prod_m_list(model_levels,nmodel_diags,nmodel_diags)

      INTEGER                                                           &
     &  prod_IND(NUM_STASH_LEVELS*2,npress_diags,npress_diags),         &
     &  prod_p_levs(npress_diags,npress_diags)

! Local parameters:

! Local scalars:
      INTEGER                                                           &
     &  I,                                                              &
     &  NI,                                                             &
     &  K,                                                              &
     &  ISL,                                                            &
     &  BL,                                                             &
     &  TL,                                                             &
     &  LEVEL,                                                          &
     &  code30,nii,nij,j

      CHARACTER(LEN=2) cfield1,cfield2

! Local dynamic arrays:
      REAL                                                              &
     &  FIELD1_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,FIELD2_PRESS(NUM_STASH_LEVELS)

      REAL:: sin_v_latitude (row_length, n_rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

!     Set to atmosphere internal model
      IF (lhook) CALL dr_hook('ST_DIAG3',zhook_in,zhook_handle)
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)



!--------------------Extract Reqd Model Levels for U-------------
      ISL=STINDEX(1,001,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_M_LIST(I)) U_M_LEVS=U_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_M_LIST(I)=.FALSE.
        ENDDO
        U_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V-------------
      ISL=STINDEX(1,002,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_M_LIST(I)) V_M_LEVS=V_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_M_LIST(I)=.FALSE.
        ENDDO
        V_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W-------------
      ISL=STINDEX(1,003,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_M_LIST(I)) W_M_LEVS=W_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_M_LIST(I)=.FALSE.
        ENDDO
        W_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T-------------
      ISL=STINDEX(1,004,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_M_LIST(I)) T_M_LEVS=T_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_M_LIST(I)=.FALSE.
        ENDDO
        T_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q-------------
      ISL=STINDEX(1,005,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_M_LIST(I)) Q_M_LEVS=Q_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_M_LIST(I)=.FALSE.
        ENDDO
        Q_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z-------------
      ISL=STINDEX(1,006,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_M_LIST(I)) Z_M_LEVS=Z_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_M_LIST(I)=.FALSE.
        ENDDO
        Z_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for KE-------------
      ISL=STINDEX(1,007,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_M_LIST(I)) KE_M_LEVS=KE_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_M_LIST(I)=.FALSE.
        ENDDO
        KE_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for Omega-------------
      ISL=STINDEX(1,008,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    OM_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        OM_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (OM_M_LIST(I)) OM_M_LEVS=OM_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          OM_M_LIST(I)=.FALSE.
        ENDDO
        OM_M_LEVS=1
      ENDIF
!--------------Extract Reqd Model Levels for Dry mass weighting----
      ISL=STINDEX(1,115,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &  dry_mass_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        DRY_MASS_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (DRY_MASS_M_LIST(I)) DRY_MASS_M_LEVS=DRY_MASS_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          DRY_MASS_M_LIST(I)=.FALSE.
        ENDDO
        DRY_MASS_M_LEVS=1
      ENDIF

!-Extract Reqd model levels for products and store in prod_m_list

      do i=1,nmodel_diags
        do j=1,nmodel_diags
          code30=i*10+j
          IF (SF(code30,30)) THEN
            ISL=STINDEX(1,code30,30,im_index)
! DEPENDS ON: set_levels_list
            CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),     &
     &        PROD_M_LIST(1,i,j),                                       &
     &        STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
          ELSE
            DO k=1,MODEL_LEVELS
              prod_M_LIST(k,i,j)=.FALSE.
            ENDDO
          ENDIF
        enddo
      enddo


!--------------------Extract Reqd Model Levels for wbig----------
      ISL=STINDEX(1,112,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig_M_LIST(I)) Wbig_M_LEVS=Wbig_M_LEVS+1
        ENDDO
      ELSE
        Wbig_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for wbig2----------
      ISL=STINDEX(1,114,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig2_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig2_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig2_M_LIST(I)) Wbig2_M_LEVS=Wbig2_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Wbig2_M_LIST(I)=.FALSE.
        ENDDO
         Wbig2_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for RH ----------
      ISL=STINDEX(1,113,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RH_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RH_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RH_M_LIST(I)) RH_M_LEVS=RH_M_LEVS+1
        ENDDO
      ELSE
        RH_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U mass weighted
      ISL=STINDEX(1,101,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_MM_LIST(I)) U_MM_LEVS=U_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_MM_LIST(I)=.FALSE.
        ENDDO
        U_MM_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V mass weighted
      ISL=STINDEX(1,102,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_MM_LIST(I)) V_MM_LEVS=V_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_MM_LIST(I)=.FALSE.
        ENDDO
        V_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W mass weighted
      ISL=STINDEX(1,103,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_MM_LIST(I)) W_MM_LEVS=W_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_MM_LIST(I)=.FALSE.
        ENDDO
        W_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T mass weighted
      ISL=STINDEX(1,104,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_MM_LIST(I)) T_MM_LEVS=T_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_MM_LIST(I)=.FALSE.
        ENDDO
        T_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q mass weighted
      ISL=STINDEX(1,105,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_MM_LIST(I)) Q_MM_LEVS=Q_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_MM_LIST(I)=.FALSE.
        ENDDO
        Q_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z mass weighted
      ISL=STINDEX(1,106,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_MM_LIST(I)) Z_MM_LEVS=Z_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_MM_LIST(I)=.FALSE.
        ENDDO
        Z_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for KE mass weighted
      ISL=STINDEX(1,107,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_MM_LIST(I)) KE_MM_LEVS=KE_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_MM_LIST(I)=.FALSE.
        ENDDO
        KE_MM_LEVS=1
      ENDIF


!--------------------Extract Reqd Model Levels for T at EOT
      ISL=STINDEX(1,111,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    TEOT_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        TEOT_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (TEOT_M_LIST(I)) TEOT_M_LEVS=TEOT_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          TEOT_M_LIST(I)=.FALSE.
        ENDDO
        TEOT_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T increment
      ISL=STINDEX(1,181,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_INC_LIST(I)) T_INC_LEVS=T_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_INC_LIST(I)=.FALSE.
        ENDDO
        T_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment
      ISL=STINDEX(1,182,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_INC_LIST(I)) Q_INC_LEVS=Q_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_INC_LIST(I)=.FALSE.
        ENDDO
        Q_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment
      ISL=STINDEX(1,183,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL_INC_LIST(I)) QCL_INC_LEVS=QCL_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCL_INC_LIST(I)=.FALSE.
        ENDDO
        QCL_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment
      ISL=STINDEX(1,184,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF_INC_LIST(I)) QCF_INC_LEVS=QCF_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCF_INC_LIST(I)=.FALSE.
        ENDDO
        QCF_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment
      ISL=STINDEX(1,185,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_INC_LIST(I)) U_INC_LEVS=U_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_INC_LIST(I)=.FALSE.
        ENDDO
        U_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment
      ISL=STINDEX(1,186,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_INC_LIST(I)) V_INC_LEVS=V_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_INC_LIST(I)=.FALSE.
        ENDDO
        V_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W increment
      ISL=STINDEX(1,187,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_INC_LIST(I)) W_INC_LEVS=W_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_INC_LIST(I)=.FALSE.
        ENDDO
        W_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment
      ISL=STINDEX(1,188,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO_INC_LIST(I)) RHO_INC_LEVS=RHO_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          RHO_INC_LIST(I)=.FALSE.
        ENDDO
        RHO_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for qrain increment
      isl=stindex(1,189,30,im_index)
      IF(isl >  0) THEN
! DEPENDS ON: set_levels_list
        CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
               qrain_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
        qrain_inc_levs = 0
        DO i=1,model_levels
          IF (qrain_inc_list(i)) qrain_inc_levs=qrain_inc_levs+1
        END DO
      ELSE
        DO i=1,model_levels
          qrain_inc_list(i)=.FALSE.
        END DO
        qrain_inc_levs = 1
      END IF

!--------------------Extract Reqd Model Levels for qgraup increment
      isl=stindex(1,190,30,im_index)
      IF(isl >  0) THEN
! DEPENDS ON: set_levels_list
        CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
               qgraup_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
        qgraup_inc_levs = 0
        DO i=1,model_levels
          IF (qgraup_inc_list(i)) qgraup_inc_levs=qgraup_inc_levs+1
        END DO
      ELSE
        DO i=1,model_levels
          qgraup_inc_list(i)=.FALSE.
        END DO
        qgraup_inc_levs = 1
      END IF

!--------------------Extract Reqd Model Levels for qcf2 increment
      isl=stindex(1,191,30,im_index)
      IF(isl >  0) THEN
! DEPENDS ON: set_levels_list
        CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
               qqcf2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
        qqcf2_inc_levs = 0
        DO i=1,model_levels
          IF (qqcf2_inc_list(i)) qqcf2_inc_levs=qqcf2_inc_levs+1
        END DO
      ELSE
        DO i=1,model_levels
          qqcf2_inc_list(i)=.FALSE.
        END DO
        qqcf2_inc_levs = 1
      END IF

!--------------------Extract Reqd Model Levels for T increment **2
      ISL=STINDEX(1,171,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T2_INC_LIST(I)) T2_INC_LEVS=T2_INC_LEVS+1
        ENDDO
      ELSE
        T2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment **2
      ISL=STINDEX(1,172,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q2_INC_LIST(I)) Q2_INC_LEVS=Q2_INC_LEVS+1
        ENDDO
      ELSE
        Q2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment **2
      ISL=STINDEX(1,173,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL2_INC_LIST(I)) QCL2_INC_LEVS=QCL2_INC_LEVS+1
        ENDDO
      ELSE
        QCL2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment **2
      ISL=STINDEX(1,174,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF2_INC_LIST(I)) QCF2_INC_LEVS=QCF2_INC_LEVS+1
        ENDDO
      ELSE
        QCF2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment **2
      ISL=STINDEX(1,175,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U2_INC_LIST(I)) U2_INC_LEVS=U2_INC_LEVS+1
        ENDDO
      ELSE
        U2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment **2
      ISL=STINDEX(1,176,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V2_INC_LIST(I)) V2_INC_LEVS=V2_INC_LEVS+1
        ENDDO
      ELSE
        V2_INC_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for W increment **2
      ISL=STINDEX(1,177,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W2_INC_LIST(I)) W2_INC_LEVS=W2_INC_LEVS+1
        ENDDO
      ELSE
        W2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment **2
      ISL=STINDEX(1,178,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO2_INC_LIST(I)) RHO2_INC_LEVS=RHO2_INC_LEVS+1
        ENDDO
      ELSE
        RHO2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Pressures for U_COMP_P-------------

      ISL=STINDEX(1,201,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        U_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,U_P_LEVS
          U_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        U_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for V_COMP_P-------------

      ISL=STINDEX(1,202,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        V_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,V_P_LEVS
          V_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        V_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for W_COMP_P-------------

      ISL=STINDEX(1,203,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        W_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,W_P_LEVS
          W_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        W_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for T on wind grid-------

      ISL=STINDEX(1,204,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        T_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,T_P_LEVS
          T_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        T_P_LEVS=1
      END IF


!--------------------Extract Reqd Pressures for q on wind grid-------

      ISL=STINDEX(1,205,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Q_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Q_P_LEVS
          Q_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Q_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for RH-----

      ISL=STINDEX(1,206,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        RH_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,RH_P_LEVS
          RH_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        RH_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for geopotential height-----

      ISL=STINDEX(1,207,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Z_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Z_P_LEVS
          Z_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Z_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for Omega---

      ISL=STINDEX(1,208,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        OM_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,OM_P_LEVS
          OM_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        OM_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for Vstarbar---

      ISL=STINDEX(1,310,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        VSTARBAR_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,VSTARBAR_P_LEVS
          VSTARBAR_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        VSTARBAR_P_LEVS=1
      END IF


!--------------------Extract Reqd Pressures for Wstarbar---

      ISL=STINDEX(1,311,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        WSTARBAR_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,WSTARBAR_P_LEVS
          WSTARBAR_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        WSTARBAR_P_LEVS=1
      END IF

!----Extract Reqd Pressures for EP Flux (Phi component)---

      ISL=STINDEX(1,312,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        FY_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,WSTARBAR_P_LEVS
          FY_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        FY_P_LEVS=1
      END IF

!----Extract Reqd Pressures for EP Flux (z component)---

      ISL=STINDEX(1,313,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        FZ_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,WSTARBAR_P_LEVS
          FZ_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        FZ_P_LEVS=1
      END IF

!-----Extract Reqd Pressures for Divergence of EP Flux---

      ISL=STINDEX(1,314,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        DIVF_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,WSTARBAR_P_LEVS
          DIVF_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        DIVF_P_LEVS=1
      END IF

!-----Extract Reqd Pressures for Meridional heat flux---

      ISL=STINDEX(1,315,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        HEATFLUX_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,HEATFLUX_P_LEVS
          HEATFLUX_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        HEATFLUX_P_LEVS=1
      END IF


!-----Extract Reqd Pressures for Meridional momentum flux---

      ISL=STINDEX(1,316,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        MOMFLUX_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,MOMFLUX_P_LEVS
          MOMFLUX_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        MOMFLUX_P_LEVS=1
      END IF


!--------------------Extract Reqd Pressures for products---


      do i=1,npress_diags
        do j=1,npress_diags
          if (i  ==  1) then
            cfield1='u'
            field1_press=u_press
            field1_p_levs=u_p_levs
          else if (i  ==  2) then
            cfield1='v'
            field1_press=v_press
            field1_p_levs=v_p_levs
          else if (i  ==  3) then
            cfield1='w'
            field1_press=w_press
            field1_p_levs=w_p_levs
          else if (i  ==  4) then
            cfield1='t'
            field1_press=t_press
            field1_p_levs=t_p_levs
          else if (i  ==  5) then
            cfield1='q'
            field1_press=q_press
            field1_p_levs=q_p_levs
          else if (i  ==  6) then
            cfield1='rh'
            field1_press=rh_press
            field1_p_levs=rh_p_levs
          else if (i  ==  7) then
            cfield1='z'
            field1_press=z_press
            field1_p_levs=z_p_levs
          else if (i  ==  8) then
            cfield1='om'
            field1_press=om_press
            field1_p_levs=om_p_levs
          endif

          if (j  ==  1) then
            cfield2='u'
            field2_press=u_press
            field2_p_levs=u_p_levs
          else if (j  ==  2) then
            cfield2='v'
            field2_press=v_press
            field2_p_levs=v_p_levs
          else if (j  ==  3) then
            cfield2='w'
            field2_press=w_press
            field2_p_levs=w_p_levs
          else if (j  ==  4) then
            cfield2='t'
            field2_press=t_press
            field2_p_levs=t_p_levs
          else if (j  ==  5) then
            cfield2='q'
            field2_press=q_press
            field2_p_levs=q_p_levs
          else if (j  ==  6) then
            cfield2='rh'
            field2_press=rh_press
            field2_p_levs=rh_p_levs
          else if (j  ==  7) then
            cfield2='z'
            field2_press=z_press
            field2_p_levs=z_p_levs
          else if (j  ==  8) then
            cfield2='om'
            field2_press=om_press
            field2_p_levs=om_p_levs
          endif
          code30=200+i*10+j
          IF (SF(code30,30)) THEN
            IF ((.NOT.SF(200+i,30)).OR.(.NOT.SF(200+j,30))) THEN
              CMESSAGE=' '
              CMESSAGE='ST_DIAG3 : '//cfield1//'*'//cfield2//           &
     &          ' error '//cfield1//' and '//cfield2//                  &
     &          ' must be requested'
              ICODE=1
              GOTO 9999
            ELSE
              NI=-STLIST(10,STINDEX(1,code30,30,im_index))
              NII=-STLIST(10,STINDEX(1,200+i,30,im_index))
              NIJ=-STLIST(10,STINDEX(1,200+j,30,im_index))

              CALL CHECK_PROD_LEVS(                                     &
     &          NUM_STASH_LEVELS,STASH_LEVELS,NI,                       &
     &          FIELD1_P_LEVS,FIELD1_PRESS,                             &
     &          FIELD2_P_LEVS,FIELD2_PRESS,                             &
     &          prod_P_LEVS(i,j),prod_IND(1,i,j))
            END IF
          ELSE
            prod_P_LEVS(i,j)=1
          ENDIF
        enddo
      enddo

!     Do checks for other pressure fields

      if (sf(302,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.sf(205,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature'//                     &
     &          ' error T*Q must also be requested'
        ICODE=1
        GOTO 9999
      endif

      if (sf(303,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.                                 &
     &  sf(205,30).and.sf(208,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature*omega'//               &
     &          ' error T*Q and omega must also be requested'
        ICODE=1
        GOTO 9999
      endif



!--------------------Extract Reqd Pressures for Heavyside function---

      ISL=STINDEX(1,301,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        HEAVY_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,HEAVY_P_LEVS
          HEAVY_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        HEAVY_P_LEVS=1
      END IF

!Lset up level information for tv
      if (sf(302,30)) then
        if (sf(204,30).and.sf(205,30)) then
          tv_p_levs=prod_p_levs(4,5)
        else
          tv_p_levs=1
          CMESSAGE='ST_DIAG3 : T and Q must be selected to calculate '//&
     &      'virtual temperature'
          ICODE=1
          GOTO 9999
        endif
      else
        tv_p_levs=1
      endif

!Lset up level information for tv*om
      if (sf(303,30)) then
        if (sf(204,30).and.sf(205,30).and.sf(208,30)) then
          tvom_p_levs=prod_p_levs(4,5)
        else
          tvom_p_levs=1
          CMESSAGE='ST_DIAG3 : T,Q and Omega'//                         &
     &      'must be selected to calculate '//                          &
     &      'virtual temperature'
          ICODE=1
          GOTO 9999
        endif
      else
        tvom_p_levs=1
      endif



!-------------------Set up Pointers for STASHWORK -------------------
      PT001=SI(001,30,im_index)
      PT002=SI(002,30,im_index)
      PT003=SI(003,30,im_index)
      PT004=SI(004,30,im_index)
      PT005=SI(005,30,im_index)
      PT006=SI(006,30,im_index)
      PT007=SI(007,30,im_index)
      PT008=SI(008,30,im_index)
      PT101=SI(101,30,im_index)
      PT102=SI(102,30,im_index)
      PT103=SI(103,30,im_index)
      PT104=SI(104,30,im_index)
      PT105=SI(105,30,im_index)
      PT106=SI(106,30,im_index)
      PT107=SI(107,30,im_index)
      PT111=SI(111,30,im_index)
      PT112=SI(112,30,im_index)
      PT113=SI(113,30,im_index)
      PT114=SI(114,30,im_index)
      PT115=SI(115,30,im_index)
      PT171=SI(171,30,im_index)
      PT172=SI(172,30,im_index)
      PT173=SI(173,30,im_index)
      PT174=SI(174,30,im_index)
      PT175=SI(175,30,im_index)
      PT176=SI(176,30,im_index)
      PT177=SI(177,30,im_index)
      PT178=SI(178,30,im_index)
      PT181=SI(181,30,im_index)
      PT182=SI(182,30,im_index)
      PT183=SI(183,30,im_index)
      PT184=SI(184,30,im_index)
      PT185=SI(185,30,im_index)
      PT186=SI(186,30,im_index)
      PT187=SI(187,30,im_index)
      PT188=SI(188,30,im_index)
      pt189=si(189,30,im_index)
      pt190=si(190,30,im_index)
      pt191=si(191,30,im_index)
      PT201=SI(201,30,im_index)
      PT202=SI(202,30,im_index)
      PT203=SI(203,30,im_index)
      PT204=SI(204,30,im_index)
      PT205=SI(205,30,im_index)
      PT206=SI(206,30,im_index)
      PT207=SI(207,30,im_index)
      PT208=SI(208,30,im_index)
      PT301=SI(301,30,im_index)
      PT302=SI(302,30,im_index)
      PT303=SI(303,30,im_index)
      PT310=SI(310,30,im_index)
      PT311=SI(311,30,im_index)
      PT312=SI(312,30,im_index)
      PT313=SI(313,30,im_index)
      PT314=SI(314,30,im_index)
      PT315=SI(315,30,im_index)
      PT316=SI(316,30,im_index)
      PT401=SI(401,30,im_index)
      PT402=SI(402,30,im_index)
      PT403=SI(403,30,im_index)
      PT404=SI(404,30,im_index)
      PT405=SI(405,30,im_index)
      PT406=SI(406,30,im_index)
      PT407=SI(407,30,im_index)
      PT408=SI(408,30,im_index)
      PT409=SI(409,30,im_index)
      PT410=SI(410,30,im_index)
      PT411=SI(411,30,im_index)
      PT412=SI(412,30,im_index)
      PT413=SI(413,30,im_index)
      PT414=SI(414,30,im_index)
      PT415=SI(415,30,im_index)
      PT416=SI(416,30,im_index)
      PT417=SI(417,30,im_index)
      PT418=SI(418,30,im_index)
      PT419=SI(419,30,im_index)
      PT420=SI(420,30,im_index)
      PT421=SI(421,30,im_index)
      PT422=SI(422,30,im_index)
      PT423=SI(423,30,im_index)
      PT424=SI(424,30,im_index)
      PT425=SI(425,30,im_index)
      PT426=SI(426,30,im_index)
      PT427=SI(427,30,im_index)
      PT428=SI(428,30,im_index)
      PT429=SI(429,30,im_index)
      PT430=SI(430,30,im_index)
      PT431=SI(431,30,im_index)
      PT432=SI(432,30,im_index)
      PT433=SI(433,30,im_index)
      PT434=SI(434,30,im_index)
      PT435=SI(435,30,im_index)
      PT436=SI(436,30,im_index)
      PT437=SI(437,30,im_index)
      PT438=SI(438,30,im_index)
      PT439=SI(439,30,im_index)
      PT440=SI(440,30,im_index)
      PT441=SI(441,30,im_index)
      PT442=SI(442,30,im_index)
      PT451=SI(451,30,im_index)
      PT452=SI(452,30,im_index)
      PT453=SI(453,30,im_index)
      PT454=SI(454,30,im_index)

! Initialise STASHWORK array because EOT_DIAG does not initialise halos
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT30
        STASHWORK(I)=0.
      ENDDO

      CALL Eot_diag(                                                    &
! Primary data: in
     &  PSTAR,P,RHO,U,V,W                                               &
     &  ,THETA,Q,QCL,QCF                                                &
     &  ,p_theta_levels ,exner_rho_levels                               &
     &  ,exner_theta_levels,energy_corr_now                             &
     &  ,inc_u, inc_v, inc_w, inc_t ,inc_q, inc_qcl, inc_qcf, inc_rho   &
        ,inc_qrain,inc_qgraup,inc_qcf2                                  &
        ,wet_to_dry_n                                                   &
! Grid sizes and definition: in
     &  ,rows,n_rows,row_length,model_levels,wet_levels,bl_levels       &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,sin_theta_latitude,cos_theta_latitude                          &
     &  ,sin_v_latitude                                                 &
     &  ,sin_theta_longitude,cos_theta_longitude                        &
     &  ,Model_domain ,delta_lambda ,delta_phi,global_row_length        &
! Grid coordinates: in
     &  ,A_REALHD(rh_deltaEW),A_REALHD(rh_deltaNS)                      &
     &  ,A_REALHD(rh_baselat),A_REALHD(rh_baselong)                     &
     &  ,A_REALHD(rh_rotlat),A_REALHD(rh_rotlong)                       &
! Pressure levels for output arrays: in
     &  ,U_PRESS,V_PRESS,W_PRESS,T_PRESS,Q_PRESS,RH_PRESS               &
     &  ,Z_PRESS,OM_PRESS,HEAVY_PRESS,VSTARBAR_PRESS,WSTARBAR_PRESS     &
     &  ,FY_PRESS,FZ_PRESS,DIVF_PRESS,HEATFLUX_PRESS,MOMFLUX_PRESS      &
! Flags to request each diagnostic output field: in
     &  ,SF(001,30),SF(002,30),SF(003,30),SF(004,30),SF(005,30)         &
     &  ,SF(006,30),SF(007,30),SF(008,30)                               &
     &  ,SF(101,30),SF(102,30),SF(103,30),SF(104,30),SF(105,30)         &
     &  ,SF(106,30),SF(107,30)                                          &
     &  ,SF(111,30),SF(112,30),SF(114,30),SF(113,30),SF(115,30)         &
     &  ,SF(181,30),SF(182,30),SF(183,30),SF(184,30),SF(185,30)         &
     &  ,SF(186,30),SF(187,30),SF(188,30),sf(189,30),sf(190,30)         &
        ,sf(191,30)                                                     &
     &  ,SF(171,30),SF(172,30),SF(173,30),SF(174,30),SF(175,30)         &
     &  ,SF(176,30),SF(177,30),SF(178,30)                               &
     &  ,SF(201,30),SF(202,30),SF(203,30)                               &
     &  ,SF(204,30),SF(205,30),SF(206,30),SF(207,30),SF(208,30)         &
     &  ,SF(301,30),SF(302,30),SF(303,30)                               &
     &  ,SF(310,30),SF(311,30),SF(312,30),SF(313,30),SF(314,30)         &
     &  ,SF(315,30),SF(316,30)                                          &
     &  ,SF(401,30),SF(402,30),SF(403,30),SF(404,30),SF(405,30)         &
     &  ,SF(406,30),SF(407,30),SF(408,30),SF(409,30),SF(410,30)         &
     &  ,SF(411,30),SF(412,30),SF(413,30),SF(414,30),SF(415,30)         &
     &  ,SF(416,30),SF(417,30),SF(418,30),SF(419,30),SF(420,30)         &
     &  ,SF(421,30),SF(422,30),SF(423,30),SF(424,30),SF(425,30)         &
     &  ,SF(426,30),SF(427,30),SF(428,30),SF(429,30),SF(430,30)         &
     &  ,SF(431,30),SF(432,30),SF(433,30),SF(434,30),SF(435,30)         &
     &  ,SF(436,30),SF(437,30),SF(438,30),SF(439,30)                    &
     &  ,SF(440,30),SF(441,30)                                          &
     &  ,SF(451,30),SF(452,30),SF(453,30),SF(454,30)                    &
        ,SF(442,30)                                                     &
! Diagnostics lengths: in
     &  ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
     &  ,ke_m_levs,om_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs    &
     &  ,q_mm_levs,z_mm_levs,ke_mm_levs,teot_m_levs,wbig_m_levs         &
     &  ,wbig2_m_levs,RH_m_levs,dry_mass_m_levs                         &
     &  ,t_inc_levs,q_inc_levs,qcl_inc_levs,qcf_inc_levs                &
        ,qrain_inc_levs,qgraup_inc_levs,qqcf2_inc_levs                  &
     &  ,u_inc_levs,v_inc_levs,w_inc_levs,rho_inc_levs                  &
     &  ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
     &  ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
     &  ,u_p_levs,v_p_levs,w_p_levs, t_p_levs,q_p_levs,rh_p_levs        &
     &  ,z_p_levs,om_p_levs,heavy_p_levs,tv_p_levs,tvom_p_levs          &
     &  ,vstarbar_p_levs,wstarbar_p_levs                                &
     &  , Fy_p_levs, Fz_p_levs, divF_p_levs                             &
     &  , heatflux_p_levs, momflux_p_levs                               &
     &  ,prod_p_levs                                                    &
                     ! pressure levels required for each product
     &  ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
     &  ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
     &  ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
     &  ,om_m_list,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list    &
     &  ,z_mm_list,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list      &
     &  ,RH_m_list,dry_mass_m_list                                      &
     &  ,t_inc_list,q_inc_list,qcl_inc_list,qcf_inc_list                &
        ,qrain_inc_list,qgraup_inc_list,qqcf2_inc_list                  &
     &  ,u_inc_list,v_inc_list,w_inc_list,rho_inc_list                  &
     &  ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
     &  ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
     &  ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
     &  ,prod_ind                                                       &
                  ! index of pressure levels required for product
! Diagnostic arrays: out
     &  ,SI(1,30,im_index) ,Stashwork                                   &
     &  ,Stashwork(pt001),Stashwork(pt002),Stashwork(pt003)             &
     &  ,Stashwork(pt004),Stashwork(pt005),Stashwork(pt006)             &
     &  ,Stashwork(pt007),Stashwork(pt008)                              &
     &  ,Stashwork(pt101),Stashwork(pt102),Stashwork(pt103)             &
     &  ,Stashwork(pt104),Stashwork(pt105),Stashwork(pt106)             &
     &  ,Stashwork(pt107),Stashwork(pt111),Stashwork(pt112)             &
     &  ,Stashwork(pt114),Stashwork(pt113),Stashwork(pt115)             &
     &  ,Stashwork(pt181),Stashwork(pt182),Stashwork(pt183)             &
     &  ,Stashwork(pt184),Stashwork(pt185),Stashwork(pt186)             &
     &  ,Stashwork(pt187),Stashwork(pt188),Stashwork(pt189)             &
        ,Stashwork(pt190),Stashwork(pt191)                              &
     &  ,Stashwork(pt171),Stashwork(pt172),Stashwork(pt173)             &
     &  ,Stashwork(pt174),Stashwork(pt175),Stashwork(pt176)             &
     &  ,Stashwork(pt177),Stashwork(pt178)                              &
     &  ,Stashwork(pt201),Stashwork(pt202),Stashwork(pt203)             &
     &  ,Stashwork(pt204),Stashwork(pt205),Stashwork(pt206)             &
     &  ,Stashwork(pt207),Stashwork(pt208)                              &
     &  ,Stashwork(pt301),Stashwork(pt302),Stashwork(pt303)             &
     &  ,Stashwork(pt310),Stashwork(pt311)                              &
     &  ,Stashwork(pt312),Stashwork(pt313),Stashwork(pt314)             &
     &  ,Stashwork(pt315),Stashwork(pt316)                              &
     &  ,Stashwork(pt401),Stashwork(pt402),Stashwork(pt403)             &
     &  ,Stashwork(pt404),Stashwork(pt405),Stashwork(pt406)             &
     &  ,Stashwork(pt407),Stashwork(pt408),Stashwork(pt409)             &
     &  ,Stashwork(pt410),Stashwork(pt411),Stashwork(pt412)             &
     &  ,Stashwork(pt413),Stashwork(pt414),Stashwork(pt415)             &
     &  ,Stashwork(pt416),Stashwork(pt417),Stashwork(pt418)             &
     &  ,Stashwork(pt419),stashwork(pt420),Stashwork(pt421)             &
     &  ,Stashwork(pt422),stashwork(pt423),Stashwork(pt424)             &
     &  ,Stashwork(pt425),stashwork(pt426),Stashwork(pt427)             &
     &  ,Stashwork(pt428),stashwork(pt429),Stashwork(pt430)             &
     &  ,Stashwork(pt431),stashwork(pt432),Stashwork(pt433)             &
     &  ,Stashwork(pt434),stashwork(pt435),Stashwork(pt436)             &
     &  ,Stashwork(pt437),stashwork(pt438),Stashwork(pt439)             &
     &  ,Stashwork(pt440),Stashwork(pt441)                              &

     &  ,Stashwork(pt451),Stashwork(pt452),Stashwork(pt453)             &
     &  ,Stashwork(pt454),Stashwork(pt442))

      IF(ICODE /= 0) THEN
        IF (lhook) CALL dr_hook('ST_DIAG3',zhook_out,zhook_handle)
        RETURN
      ENDIF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,30,STASHWORK,                                &
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
      IF (lhook) CALL dr_hook('ST_DIAG3',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ST_DIAG3
