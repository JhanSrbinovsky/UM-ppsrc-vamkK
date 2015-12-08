! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!  To gather required fields from D1, then call the UK Chemistry
!  and Aerosols submodel components.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Method:
! 1) Interface with atmosphere model:
!   Prognostic fields are selected from D1 using the item codes
!   defined in the routine UKCA_SetD1Defs, these include tracers
!   from section 34 (UKCA).
!   Diagnostic fields are selected from D1 using a tag of 98, these
!   need to be set in the STASH panel of the UMUI.
! 2) UKCA routines are called depending on the scheme selected:
!    - Woodward dust scheme
!    - Emissions routine
!    - photolysis routine
!    - Chemistry control routine
! 3) Updated tracer arrays are returned to the D1 array
!
! CONTAINED subroutines:  GETD1FLDS, PUTD1FLDS
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_MAIN1(                                           &
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
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
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
         Idummy)

      USE dynamics_input_mod, ONLY: l_endgame
      USE dynamics_grid_mod,  ONLY: l_vatpoles

      USE atm_fields_bounds_mod
      USE bl_option_mod,     ONLY: alpha_cd
      USE conversions_mod,   ONLY: pi
      USE ASAD_MOD,          ONLY: asad_mod_final, method, nnaf,        &
                                   dtime, interval, advt, kcdt,         &
                                   FCH4, FCO2, FH2, FN2, FO2,           &
                                   spt, ntrkx, jpspt
      USE UKCA_D1_DEFS
      USE UKCA_TRACER_VARS
      USE UKCA_CSPECIES
      USE UKCA_CONSTANTS
      USE UKCA_TROPOPAUSE
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE UKCA_CHEM1_DAT
      USE UKCA_CHEM_STD_TROP,   ONLY: nhet_std_trop
      USE UKCA_CHEM_TROPISOP,   ONLY: nhet_tropisop
      USE UKCA_TRACE_GAS_MIXRATIO
      USE ASAD_FLUX_DAT,        ONLY: ASAD_LOAD_DEFAULT_FLUXES
      USE ASAD_CHEM_FLUX_DIAGS, ONLY: ASAD_ALLOCATE_CHEMDIAG ,          &
                                      ASAD_SETSTASH_CHEMDIAG,           &
                                      ASAD_INIT_CHEMDIAG,               &
                                      ASAD_TENDENCY_STE,                &
                                      ASAD_MASS_DIAGNOSTIC,             &
                                      ASAD_OUTPUT_TRACER,               &
                                      ASAD_FLUX_PUT_STASH,              &
                                      calculate_STE,                    &
                                      calculate_tendency,               &
                                      L_asad_use_chem_diags,            &
                                      L_asad_use_STE,                   &
                                      L_asad_use_tendency,              &
                                      L_asad_use_mass_diagnostic,       &
                                      L_asad_use_output_tracer,         &
                                      L_asad_use_trop_mask,             &
                                      asad_tropospheric_mask
      USE ukca_option_mod,      ONLY: l_ukca_set_trace_gases,           &
                                      l_ukca_chem, l_ukca_trop,         &
                                      l_ukca_aerchem, l_ukca_raq,       &
                                      l_ukca_intdd, l_ukca_advh2o,      &
                                      l_ukca_mode, l_ukca_strattrop,    &
                                      l_ukca_strat, l_ukca_stratcfc,    &
                                      l_ukca_trophet, l_ukca_std_trop,  &
                                      l_ukca_het_psc, l_ukca_sa_clim,   &
                                      l_ukca_use_background_aerosol,    &
                                      l_ukca_dust, l_ukca_achem,        &
                                      l_ukca_qch4inter, l_ukca_ageair,  &
                                      l_ukca_prescribech4,              &
                                      l_ukca_useumuivals,               &
                                      l_ukca_arg_act,                   &
                                      l_ukca_h2o_feedback,              &
                                      jpctr, jpspec, jppj, jpdd, jpdw,  &
                                      i_mode_setup, l_ukca,             &
                                      i_ukca_photol
      USE ukca_photo_scheme_mod,ONLY: i_ukca_fastj,                     &
                                      i_ukca_fastjx
      USE level_heights_mod,    ONLY: r_theta_levels, r_rho_levels,     &
                                      eta_theta_levels
      USE trignometric_mod,     ONLY: true_longitude,                   &
                                      sin_theta_longitude,              &
                                      sin_theta_latitude,               &
                                      sin_v_latitude,                   &
                                      cos_v_latitude,                   &
                                      FV_cos_theta_latitude,            &
                                      cos_theta_longitude,              &
                                      tan_theta_latitude
      USE dyn_coriolis_mod,     ONLY: f3_at_u
      USE earth_constants_mod,  ONLY: g, two_omega
      USE dust_parameters_mod,  ONLY: ndiv
      USE nstypes,              ONLY: ntype, npft, soil
      USE dec_spec
      USE spec_sw_lw
      USE rad_input_mod,        ONLY: a_sw_radstep_prog,l_sec_var,l_eqt
      USE rad_com_mod,          ONLY: n2ommr, c12mmr,c11mmr,c113mmr,    &
                                      hcfc22mmr,ch4mmr
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE turb_diff_mod,        ONLY: q_pos_method, qlimit,             &
                    l_qpos_diag_pr, qpos_diag_limit, q_pos_tracer_method
      USE ereport_mod,          ONLY: ereport
      USE PrintStatus_mod
      USE Field_Types
      USE UM_ParParams
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon,      &
                            io3_trop_map, io3_trop_map_masscon
      USE problem_mod, ONLY: standard, monsoon, dynamical_core,         &
                             idealised_problem, standard_namelist
      USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
      USE um_input_control_mod,  ONLY:                                  &
           model_domain,                                                &
           h_sect,          lcal360
           
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
      USE Submodel_Mod
      USE nlstcall_mod, ONLY : ltimer

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
! River routing
! TYPATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of TYPAOCPL.h)
!
      REAL :: XPA(AOCPL_ROW_LENGTH+1)  ! Atmosphere TP longitude coordina
      REAL :: XUA(0:AOCPL_ROW_LENGTH)  ! Atmosphere U longitude coordinat
      REAL :: XVA(AOCPL_ROW_LENGTH+1)  ! Atmosphere V longitude coordinat
      REAL :: YPA(AOCPL_P_ROWS)        ! Atmosphere TP latitude coordinat
      REAL :: YUA(AOCPL_P_ROWS)        ! Atmosphere U latitude coordinate
      REAL :: YVA(0:AOCPL_P_ROWS)      ! Atmosphere V latitude coordinate
! TYPATCPL end
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
! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10
! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11
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
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
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

      INTEGER, INTENT(INOUT)     :: Idummy       !

! ErrorStatus

      INTEGER                    :: errcode=0     ! Error flag (0 = OK)
      CHARACTER(LEN=256)             :: cmessage      ! Error return message

! Local scalars

      INTEGER    :: Nfields           ! fields counter
      INTEGER    :: Nukca_Fields      ! No. fields requested
      INTEGER    :: len               ! local dimension
      INTEGER    :: levs              ! number of levels
      INTEGER    :: StashCode         ! stash code
      INTEGER    :: section           ! stash section
      INTEGER    :: item              ! stash item
      INTEGER    :: addr              ! stash item
      INTEGER    :: field_typ         ! Field type
      INTEGER    :: halo_typ          ! halo type
      INTEGER    :: grid_typ          ! grid type
      INTEGER    :: tag               ! stash tag
      INTEGER    :: ptd1              ! D1 pointer
      INTEGER    :: I,I1,I2           ! loop variables
      INTEGER    :: ii                ! loop variables
      INTEGER    :: J,J1,J2           ! loop variables
      INTEGER    :: K,L,N             ! loop variables
      INTEGER    :: KK                ! loop counter
      INTEGER    :: m_atm_modl        ! sub model
      INTEGER    :: GET_FLD_TYPE      ! UM function
      INTEGER    :: timestep_number   ! no. of atmos timesteps since bas
      INTEGER    :: jradon_em         ! alternative D1 pointer for rn em
      INTEGER    :: nd_o3             ! size of um ozone array
      INTEGER    :: index2            ! Indices for controlling
      INTEGER    :: index3            ! albedo calculation
      INTEGER    :: A_Steps_per_hr    ! Atmos steps per hour
      INTEGER    :: Radsteps_per_day  ! Radiation TS per day
      INTEGER    :: A_SW_RADSTEP_new  ! To take acc of radn 3A/3C 
      INTEGER    :: icnt              ! counter
      INTEGER    :: icode=0           ! local error status 
      INTEGER    :: im_index          ! internal model index 
      INTEGER, SAVE :: n_day          ! test whether new day has been reached
      INTEGER    :: het_dimn          ! 1st dimension of het_rates array
      INTEGER    :: asad_findreaction ! integer function
      INTEGER, SAVE :: wetox_in_aer   ! set for wet oxidation in MODE (=1)
                                      ! or UKCA (=0)
      INTEGER, SAVE :: uph2so4inaer   ! update H2SO4 tracer in MODE (=1)
                                      ! or UKCA (=0) depending on chemistry

      INTEGER, SAVE :: nhet_value     ! number of reactions on aerosol surf

      REAL       :: timestep                 ! atmosphere model TS
      REAL       :: min_surf_albedo=0.1      ! fill-in value
      REAL       :: fx                       ! interpolation fraction
      REAL       :: eq_time                  ! Equation of time
      REAL       :: Sindec                   ! Solar declination
      REAL       :: secondssincemidnight     ! day time
      REAL       :: SCS                      ! solar constant scaling factor
      REAL       :: ssmn_incr                ! sec. since midnight incr. counter
      
      CHARACTER(LEN=10)  :: prods(2)              ! Products
      LOGICAL, SAVE :: first       = .TRUE.  ! true only on 1st call
      LOGICAL, SAVE :: firstchem   = .TRUE.  ! true only on 1st call - chemistry
      LOGICAL       :: L_moses_ii                 ! T if MOSES_II in use

! Local arrays

! MPP-related Arrays
! Table of number of points on a row
      INTEGER, DIMENSION(0:nproc-1) :: g_row_length
! Table number of rows in theta field
      INTEGER, DIMENSION(0:nproc-1) :: g_rows

! Tile index etc
      INTEGER, DIMENSION(land_points,ntype) :: tile_index !
      INTEGER, DIMENSION(ntype)             :: tile_pts   ! No of tile p

! arrays filled from D1, allocated dynamically with correct halo read
!  in from addressing array
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: conv_cloud_base
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: conv_cloud_top
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: kent
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: kent_dsc

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: land_sea_mask

      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: all_tracers
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: trmol_post_atmstep
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: chem_diags
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: mode_diags
      REAL, DIMENSION(:,:),     ALLOCATABLE :: het_rates
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: strat_fluxdiags
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: all_emissions
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: aircraftems   ! aircraft NOx
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SO2_volc_3D   ! volcanic SO2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SO2_biom_3D   ! biomass SO2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: BC_biom_3D    ! biomass BC
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: OC_biom_3D    ! biomass OC
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: tracer1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: chem_diag1
      REAL, DIMENSION(:,:),   ALLOCATABLE :: emission1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: theta
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_rho_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_theta_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: totnodens         ! density in molecs/m^3
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_theta_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_rho_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: pv_on_theta_mlevs
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: q
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: qsvp
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rel_humid_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rho_r2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: qcl
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: qcf
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_ppn_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_liq_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_liq_water
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_ice_content
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_rain3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_snow3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_ppn3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_rain3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_snow3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_ppn3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_cloud_amount
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: area_cloud_fraction
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: so4_aitken
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: so4_accum
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: sulphate_od  ! optical depth
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wet_h2o2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wet_o3
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delh2so4_chem  ! Rate of net chem Pdn H2SO4
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_drydep
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wetdep
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: um_ozone   ! i.e. from D1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: um_ozone3d ! after O3_to_3D
      REAL, DIMENSION(:),     ALLOCATABLE :: um_ozone1d ! ip to O3_to_3D
      REAL, DIMENSION(:,:), ALLOCATABLE   :: pstar
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tstar
      REAL, DIMENSION(:,:), ALLOCATABLE   :: U_scalar_10m
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: dust_flux
! dust emissions from 6 bin scheme
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Soil_Layer_Moisture
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tile_Frac
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tstar_tile
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Snow_tile
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Frac_types
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Rough_length
! surf_albedo is interpolated at every timestep from land_albedo_all
! which is calculated on radiation timesteps using the SW fluxes.
      REAL, DIMENSION(:,:), ALLOCATABLE   :: land_albedo
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: land_albedo_all  ! all r
      REAL, DIMENSION(:,:), ALLOCATABLE   :: surf_albedo  ! interpolated
      REAL, DIMENSION(:,:), ALLOCATABLE   :: net_surf_SW
      REAL, DIMENSION(:,:), ALLOCATABLE   :: tot_surf_SW
      REAL, DIMENSION(:), ALLOCATABLE     :: fland      ! land fraction
      REAL, DIMENSION(:,:), ALLOCATABLE   :: climoz2d !climatological O3
      REAL, DIMENSION(:,:), ALLOCATABLE   :: zbl      ! BL height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: surf_hf
      REAL, DIMENSION(:,:), ALLOCATABLE   :: seaice_frac
      REAL, DIMENSION(:,:), ALLOCATABLE   :: conv_cloud_lwp
      REAL, DIMENSION(:,:), ALLOCATABLE   :: u_s     ! surf frict velocity
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: stcon       ! stomatal cond
      REAL, DIMENSION(:,:), ALLOCATABLE   :: laift_lp    ! LAI
      REAL, DIMENSION(:,:), ALLOCATABLE   :: canhtft_lp  ! canopy height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: z0tile_lp
      REAL, DIMENSION(:,:), ALLOCATABLE   :: canwctile_lp ! canopy WC
      REAL, DIMENSION(:,:), ALLOCATABLE   :: ch4_wetl_emiss
      REAL, DIMENSION(:,:), ALLOCATABLE   :: theta_latitude ! gridbox lat
      REAL, DIMENSION(:,:), ALLOCATABLE   :: v_latitude     ! boundary lat
      REAL, DIMENSION(:,:), ALLOCATABLE   :: sinv_latitude  ! sin(boundary lat)
      REAL, DIMENSION(:,:), ALLOCATABLE   :: tropopause_height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: cos_zenith_angle
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE   :: int_zenith_angle ! integral over sza      
      REAL, DIMENSION(:,:), ALLOCATABLE   :: ml_depth
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rhokh_mix
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: dtrdz_charney_grid
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: we_lim
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: zrzi
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: we_lim_dsc
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t_frac_dsc
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: zrzi_dsc
      REAL, DIMENSION(:,:), ALLOCATABLE   :: zhsc
      REAL, DIMENSION(:,:), ALLOCATABLE   :: rb_dust_div1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: bl_tke
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: vertvel
! solid angle of grid cell saved for mass calculation
      REAL, SAVE, ALLOCATABLE :: mass (:,:,:)
      REAL, SAVE, ALLOCATABLE :: solid_angle(:,:)
      REAL, SAVE, ALLOCATABLE :: z_half(:,:,:)
      REAL, SAVE, ALLOCATABLE :: z_half_alllevs(:,:,:)
      REAL, SAVE, ALLOCATABLE :: volume(:,:,:)   ! gridbox volume
      REAL, SAVE, ALLOCATABLE :: area(:,:,:)     ! gridbox area
      REAL, SAVE, ALLOCATABLE :: delta_r(:,:,:)  ! delta radius

! Local dynamic arrays
      REAL, ALLOCATABLE :: STASHwork34(:)
      REAL, ALLOCATABLE :: STASHwork38(:)
      REAL, ALLOCATABLE :: STASHwork50(:)

! Derived data arrays, allocatable if necessary
      REAL, DIMENSION(row_length,rows)       :: land_fraction  !
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: t_theta_levels ! temp on
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: Thick_bl_levels    ! val
!                                               thickness in metres of t
!                                               distance from boundaries
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: p_layer_boundaries
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: dj     ! photolysis rates
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: so4_sa ! aerosol surface area
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: T_chem
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: Q_chem

      REAL :: delta_latitude
      REAL :: delta_longitude

      INTEGER, SAVE :: fastj_levels    ! number of levels on which to do FAST-J
      REAL          :: z_top_of_model  ! top of model

! Variables for COS initialisation
      REAL :: y
      REAL :: ymin
      REAL :: ymax
      REAL :: pmin
      REAL :: pmax
      REAL :: lp
      REAL :: lpmin
      REAL :: lpmax

      REAL, SAVE :: lambda_aitken, lambda_accum ! parameters for computation
                                                ! of surface area density
      LOGICAL :: do_chemistry

      ! required for ASAD Flux Diagnostics
      INTEGER :: ierr, stashsize
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: fluxdiag_all_tracers

! Mixing ratios for longlived tracers (from WMO, 2002)
! Short-lived stratospheric tracers are forced to 0 at the surface.
      REAL, PARAMETER :: my_n2o_mmr = 4.7842E-7
      REAL, PARAMETER :: f12_mmr    = 2.2145E-9
      REAL, PARAMETER :: f11_mmr    = 1.3294E-9
      REAL, PARAMETER :: mebr_mmr   = 39.1206E-12
      REAL, PARAMETER :: hcl_mmr    = 1.5125E-10
      REAL, PARAMETER :: f113_mmr   = 5.3072E-10
      REAL, PARAMETER :: h22_mmr    = 4.1801E-10
      REAL, PARAMETER :: ccl4_mmr   = 5.3158E-10
      REAL, PARAMETER :: mecl_mmr   = 9.4133E-10
      REAL, PARAMETER :: meccl3_mmr = 2.3041E-10
      REAL, PARAMETER :: h1211_mmr  = 2.83858e-11
      REAL, PARAMETER :: h1301_mmr  = 1.59722e-11
      REAL, PARAMETER :: h2_mmr     = 3.4520E-8
      REAL, PARAMETER :: cos_mmr    = 5.2E-10  
      REAL, PARAMETER :: ch4_mmr    = 9.8530E-7
      REAL, PARAMETER :: ch2br2_mmr = 18.0186e-12  ! 3 pptv

      INTEGER, PARAMETER :: max_zero_age = 10    ! Level below which age of
                                                 ! air set to zero

      INTEGER :: klev1        ! First Level for T/Tr/Q grid
      INTEGER :: tr_levs      ! Levels for Tracer grid
      INTEGER :: th_levs      ! Levels for Theta grid
      INTEGER :: q_levs       ! Levels for Q grid
      INTEGER :: w_levs       ! Levels for W grid

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 
      REAL(KIND=jprb)               :: zhook_handle 

!- End of header

! timestep information
      IF (lhook) CALL dr_hook('UKCA_MAIN1',zhook_in,zhook_handle)

      z_top_of_model = A_realhd(rh_z_top_theta)

      m_atm_modl = SUBMODEL_FOR_SM(a_im)   ! submodel code
      timestep = SECS_PER_STEPim(atmos_im) ! timestep in seconds
      timestep_number = STEPim(atmos_im)   ! no. of steps since basis
      A_Steps_per_hr = 3600*STEPS_PER_PERIODim(a_im)/                  &
                            SECS_PER_PERIODim(a_im)
      radsteps_per_day=A_Steps_per_hr*24/a_sw_radstep_prog

! Initialise levels for Tr/T/Q/W grids
      klev1   = tdims%k_start
      tr_levs = trdims_xstl%k_end - trdims_xstl%k_start + 1
      th_levs = tdims%k_end - tdims%k_start + 1
      q_levs  = qdims%k_end - qdims%k_start + 1
      w_levs  = wdims%k_end - wdims%k_start + 1

! Set mpp arrays from parvars.h information

      DO i= 0,nproc-1
        g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
        g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
      END DO ! i processors

! Allocate diagnostic space for STASH
      ALLOCATE (STASHwork34(STASH_maxlen(34,A_im)))
      ALLOCATE (STASHwork38(STASH_maxlen(38,A_im)))
      ALLOCATE (STASHwork50(STASH_maxlen(50,A_im)))
      STASHwork34=0.0
      STASHwork38=0.0
      STASHwork50=0.0

! Initialise ASAD constants - take UM values if requested 
!  (values calculated in atmos_physics1 called from atm_step)
! code initially from asad_cinit (called from ukca_iniasad), 
! but moved here since may need to be updated every timestep
      IF (L_ukca_set_trace_gases) THEN
        fco2 = um_co2_for_ukca/c_co2
        fh2  = um_h2_for_ukca/c_h2
        fn2  = um_n2_for_ukca/c_n2
        fo2  = um_o2_for_ukca/c_o2
        fch4 = um_ch4_for_ukca/c_ch4
      ELSE ! present day values
        fco2 = 350.0e-6
        fh2  = 5.0e-7
        fn2  = 0.78084
        fo2  = 0.20945
        fch4 = 1.76e-6
      END IF
      
      IF (first) THEN
        IF (PrintStatus >= PrStatus_Oper) THEN
          WRITE(6,*) 'First call to ukca_main'
          WRITE(6,*) 'timestep = ',timestep
          WRITE(6,*) 'timestep number = ',timestep_number
        END IF

        ! Initialise chemical definition arrays - needs to be here in case we are
        ! using ASAD
        CALL UKCA_CHEM1_INIT()


        IF (L_ukca_chem) THEN
           ! need to set n_chem_diags here, since it is used by ASAD_INRATS which is called from UKCA_INIASAD
           n_chem_diags = 22
          dtime = timestep
! DEPENDS ON: ukca_iniasad
          CALL UKCA_INIASAD(theta_field_size)

! Initialise chemical diagnostics for N-R chemistry; set L_asad_use_chem_diags
          IF (.NOT. (L_ukca_trop .OR. L_ukca_aerchem .OR. L_ukca_raq))  &
                THEN
             CALL ASAD_LOAD_DEFAULT_FLUXES() 
             ! check to see if chemical diagnostics are required, and 
             ! set up stash for them.
            CALL ASAD_SETSTASH_CHEMDIAG(row_length,rows,model_levels)
            ! initialise chemical diagnostics
            IF (L_asad_use_chem_diags) &
                 call asad_init_chemdiag(icode)
          END IF

        END IF  ! L_ukca_chem
        IF (PrintStatus >= PrStatus_Oper) THEN
          WRITE(6,*) 'First call to ukca_main'
          WRITE(6,*) 'timestep = ',timestep
          WRITE(6,*) 'timestep number = ',timestep_number
        END IF
!
! Check consistency of dry deposition scheme and number of tiles
        IF ( (ntiles /= ntype) .AND. l_ukca_intdd ) THEN
          WRITE(6,*) 'Cannot use interactive dry dep. You have ',      &
          ntiles,' surface tiles but should have ',ntype
          errcode = 1 
          cmessage='UKCA: Cannot use interactive dry dep'

          CALL EREPORT('UKCA_MAIN',errcode,cmessage)
        END IF

! Define D1 definitions, emissions and diagnostics required
! DEPENDS ON: ukca_setd1defs
        CALL UKCA_SETD1DEFS(row_length,rows,n_rows, model_levels,      &
                 bl_levels,tr_levels,wet_levels,land_points,sm_levels, &
                 ntiles,tr_ukca)

        ! UKCA_CALC_CSPECIES needs to be called after SETD1DEFS
        IF (L_ukca_chem) THEN
          ALLOCATE(c_species(n_chem_tracers+n_aero_tracers))
          ALLOCATE(c_na_species(jpspec-jpctr))
          CALL UKCA_CALC_CSPECIES()

          errcode = 0

          DO i=1,jpctr        ! to allow for case when jpctr < n_chem_tracers...
            IF (c_species(i) < 0.001) THEN
! check that ratio of molec. wts. array has been initialised.
              errcode = errcode +1
              cmessage=' c_species array is zero for '//advt(i)
              CALL EREPORT('UKCA_MAIN',-1*errcode,cmessage)
            END IF
          END DO

          IF (errcode > 0) THEN
            cmessage=' c_species array has zero values'//               &
                         ', check UKCA_CSPECIES routine'
              CALL EREPORT('UKCA_MAIN',errcode,cmessage)
          END IF

          DO i=1,nnaf
            IF (c_na_species(i) < 0.001) THEN
! check that ratio of molec. wts. array has been initialised.
              cmessage=' c_na_species array contains zero value'
              CALL EREPORT('UKCA_MAIN',i,cmessage)
            END IF
          END DO
! L_ukca_advh2o *MUST* be .TRUE. if using H2O as 'TR'
          IF ((n_h2o > 0) .AND. (.NOT. L_ukca_advh2o)) THEN
             cmessage=                                                  &
                  ' H2O is defined as a tracer, but L_UKCA_ADVH2O=F'
             CALL EREPORT('UKCA_MAIN',n_h2o,cmessage)
          END IF
        END IF  ! L_ukca_chem


! Define UKCA-mode setup based on value of i_mode_setup from UMUI
! below are called from modules UKCA_SETUP_INDICES and UKCA_MODE_SETUP
        IF (L_ukca_mode) THEN
          IF (PrintStatus >= PrStatus_Oper) WRITE(6,*)                  &
            'MODE aerosol, i_mode_setup=',i_mode_setup
          CALL UKCA_MODE_IMSCAVCOFF
          IF (PrintStatus >= PrStatus_Oper)                             &
            WRITE(6,*) 'Set up impaction scavenging coeffs,NROW,NCOLL=',&
                                                         NROW,NCOLL
          IF (PrintStatus >= PrStatus_Diag) THEN
            DO I=1,NROW
             WRITE(6,*) 'I,RADDROP(I)=',I,RADDROP(I)
             DO J=1,NCOLL
              WRITE(6,*) 'I,J,COLLEFF4(J,I)=',I,J,COLLEFF4(J,I)
             END DO
            END DO
          END IF
          CALL UKCA_MODE_TRACER_INDICES
          IF(i_mode_setup == 1) THEN
            CALL UKCA_INDICES_SV1
            CALL UKCA_INDICES_SUSS_4MODE
            CALL UKCA_MODE_SUSS_4MODE
          ELSE IF(i_mode_setup == 2) THEN
            CALL UKCA_INDICES_ORGV1_SOto3
            CALL UKCA_INDICES_SUSSBCOC_5MODE
            CALL UKCA_MODE_SUSSBCOC_5MODE
          ELSE IF(i_mode_setup == 3) THEN
            CALL UKCA_INDICES_ORGV1_SOto3
            CALL UKCA_INDICES_SUSSBCOC_4MODE
            CALL UKCA_MODE_SUSSBCOC_4MODE
          ELSE IF(i_mode_setup == 4) THEN
            CALL UKCA_INDICES_ORGV1_SOto6
            CALL UKCA_INDICES_SUSSBCOCSO_5MODE
            CALL UKCA_MODE_SUSSBCOCSO_5MODE
          ELSE IF(i_mode_setup == 5) THEN
            CALL UKCA_INDICES_ORGV1_SOto6
            CALL UKCA_INDICES_SUSSBCOCSO_4MODE
            CALL UKCA_MODE_SUSSBCOCSO_4MODE
          ELSE IF(i_mode_setup == 6) THEN
            CALL UKCA_INDICES_NOCHEM
            CALL UKCA_INDICES_DUonly_2MODE
            CALL UKCA_MODE_DUonly_2MODE
!!        ELSE IF(i_mode_setup == 7) THEN
!!          CALL UKCA_INDICES_NOCHEM
!!          CALL UKCA_INDICES_DUonly_3MODE
!!          CALL UKCA_MODE_DUonly_3MODE
!!        ELSE IF(i_mode_setup == 8) THEN
!!          CALL UKCA_INDICES_ORGV1_SOto3
!!          CALL UKCA_INDICES_SUSSBCOCDU_7MODE
!!          CALL UKCA_MODE_SUSSBCOCDU_7MODE
!!        ELSE IF(i_mode_setup == 9) THEN
!!          CALL UKCA_INDICES_ORGV1_SOto3
!!          CALL UKCA_INDICES_SUSSBCOCDU_4MODE
!!          CALL UKCA_MODE_SUSSBCOCDU_4MODE
          ELSE
            cmessage=' i_mode_setup has unrecognised value'
            write(6,*) cmessage,i_mode_setup
            CALL EREPORT('UKCA_MAIN',i_mode_setup,cmessage)
          END IF       ! i_mode_setup
          WRITE(6,*) 'Set up aerosol etc., NTRAER,NBUDAER=',            &
                                           NTRAER,NBUDAER
          DO I=1,NMODES
           WRITE(6,*) 'I,MODE(I),DDPLIM0(I),DDPLIM1(I),SIGMAG(I)=',     &
            I,MODE(I),DDPLIM0(I),DDPLIM1(I),SIGMAG(I)
           DO J=1,NCP
            WRITE(6,*) 'I,J,COMPONENT(I,J),MFRAC_0(I,J)=',              &
                    I,J,COMPONENT(I,J),MFRAC_0(I,J)
           END DO
          END DO
        END IF    ! L_ukca_mode
      END IF ! first

! Loop through all the objects in D1 to find address and characteristics
!  of the requested items
      IF (first) THEN
        IF (PrintStatus >= PrStatus_Oper) THEN
          WRITE(6,*) 'UKCA Prognostics and Diagnostics from D1:'
          WRITE(6,*) 'section,item,levels,length,address,halo_type,',    &
                 'grid_type,field_type'
        END IF
        DO i=1,no_obj_D1(m_atm_modl)
          section = d1_addr(D1_section,i,m_atm_modl)
          item = d1_addr(D1_item,i,m_atm_modl)
          levs = d1_addr(d1_no_levels,i,m_atm_modl)
          len = d1_addr(d1_length,i,m_atm_modl)
          addr    = d1_addr(d1_address,i,m_atm_modl)
          halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
          grid_typ = d1_addr(d1_grid_type,i,m_atm_modl)
! DEPENDS ON: get_fld_type
          field_typ = get_fld_type(grid_typ)
          DO J=1,Nukca_D1items
            IF (UkcaD1Codes(J)%section == section .AND.                &
                UkcaD1Codes(J)%item == item .AND.                      &
                d1_addr(d1_object_type,i,m_atm_modl) == prognostic     &
                .AND. UkcaD1codes(J)%Prognostic) THEN
              UkcaD1Codes(J)%n_levels=levs
              UkcaD1Codes(J)%address=addr
              UkcaD1Codes(J)%length=len
              UkcaD1Codes(J)%halo_type=halo_typ
              UkcaD1Codes(J)%grid_type=grid_typ
              UkcaD1Codes(J)%field_type=field_typ
              IF (PrintStatus >= PrStatus_Diag)                        &
                WRITE(6,*) 'P',j, section,item,levs,len,addr,          &
                             halo_typ,grid_typ,field_typ
            END IF
          END DO
        END DO

!       Diagnostics loop through all stashlist items, check
!       for tag value, read adress etc into UkcaD1Codes.

        DO L=1,totitems
          tag=STLIST(st_macrotag,L)-(1000*(STLIST(st_macrotag,L)/1000))
          IF (tag == 98) THEN
            ptd1      = STLIST(st_D1pos,L)
            section   = d1_addr(d1_section,  ptd1,m_atm_modl)
            item      = d1_addr(d1_item,     ptd1,m_atm_modl)
            levs      = d1_addr(d1_no_levels,ptd1,m_atm_modl)
            len       = d1_addr(d1_length,   ptd1,m_atm_modl)
            addr      = d1_addr(d1_address,  ptd1,m_atm_modl)
            halo_typ  = d1_addr(d1_halo_type,ptd1,m_atm_modl)
            grid_typ  = d1_addr(d1_grid_type,ptd1,m_atm_modl)
! DEPENDS ON: get_fld_type
            field_typ = get_fld_type(grid_typ)
            stashcode = section*1000 + item
            DO J=1,Nukca_D1items
              IF (UkcaD1Codes(J)%section == section .AND.              &
                  UkcaD1Codes(J)%item    == item    .AND.              &
                  .NOT. UkcaD1Codes(J)%Prognostic) THEN
                UkcaD1Codes(J)%n_levels  = levs
                UkcaD1Codes(J)%address   = addr
                UkcaD1Codes(J)%length    = len
                UkcaD1Codes(J)%halo_type = halo_typ
                UkcaD1Codes(J)%grid_type = grid_typ
                UkcaD1Codes(J)%field_type= field_typ
                IF (PrintStatus >= PrStatus_Diag)                      &
                  WRITE(6,*) 'D',section,item,levs,len,addr,           &
                            halo_typ,grid_typ,field_typ
              END IF
            END DO
          END IF  ! tag
        END DO    ! L

! set constants needed in sulphur chemistry
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
         (L_ukca_stratcfc))  THEN
          lambda_aitken = 3./rho_so4 / chi /                            &
                         rad_ait * EXP(-2.5 * (ALOG(sigma)**2))         &
                         * 0.01 ! revert from m^-1 to cm^2/cm^3
          lambda_accum  = 3./rho_so4 / chi /                            &
                         rad_acc * EXP(-2.5 * (ALOG(sigma)**2))         &
                         * 0.01 ! revert from m^-2 to cm^2/cm^3
      END IF

! set up timestep counting. Interval depends on solver: IMPACT and Rosenbrock
! are run every dynamical timestep. Newton-Raphson is run every 2 / 3
! timesteps for a 30/20 minutes dynamical timestep.

        IF (method == 3) THEN
          interval = INT(kcdt/INT(timestep)) 
        ELSE
          interval = 1
        END IF

      END IF  ! first

! allocate ASAD chemical diagnostics and load the requested diagnostics
      IF (L_asad_use_chem_diags)  THEN
        CALL ASAD_ALLOCATE_CHEMDIAG(row_length,rows)
      END IF


! decide whether to do chemistry
! interval and chemical timestep is set in asad_cinit
      do_chemistry = (MOD(timestep_number, interval) == 0)


! Check if all items selected have been identified in D1
      errcode = 0
      DO i=1,Nukca_D1items
        IF (UkcaD1codes(i)%address == IMDI .AND.                       &
            UkcaD1codes(i)%required) THEN
          errcode  = errcode + 1
          cmessage = 'Item address not found in D1 array:  '
          WRITE (6,'(a37,i5)') cmessage (1:37),                        &
            UkcaD1Codes(i)%section*1000 + UkcaD1Codes(i)%item
        END IF
      END DO
      IF (errcode > 0) THEN
          cmessage = 'Some item addresses not found in D1 array'
          CALL EREPORT('UKCA_MAIN', errcode, cmessage)
      END IF

!     Copy fields from D1 array into named item, allocate arrays

      DO i=1,Nukca_D1items
        IF (UkcaD1Codes(i)%required) THEN
          CALL GETD1FLDS(i)
        END IF
      END DO



!     Calculate any derived variables, SO2 fluxes are required in chemistry_ctl

      IF (L_ukca_chem) THEN
        i1 = LBOUND(theta,dim=1)
        i2 = UBOUND(theta,dim=1)
        j1 = LBOUND(theta,dim=2)
        j2 = UBOUND(theta,dim=2)
        ALLOCATE(t_theta_levels(I1:I2,J1:J2,model_levels))
        ALLOCATE(rel_humid_frac(I1:I2,J1:J2,wet_levels))
        ALLOCATE(qsvp(1:row_length,1:rows,wet_levels))
        ALLOCATE(Thick_bl_levels(1:row_length,1:rows,bl_levels))
        ALLOCATE(tile_frac(land_points,ntiles))              ! see below
        ALLOCATE(p_layer_boundaries(row_length,rows,0:model_levels))
        ALLOCATE(surf_albedo(row_length,rows))
        ALLOCATE(totnodens(row_length,rows,model_levels))
        ALLOCATE(ls_ppn3d(row_length,rows,model_levels))
        ALLOCATE(conv_ppn3d(row_length,rows,model_levels))
        ALLOCATE(so4_sa(row_length, rows, model_levels)) ! sulphate area density
        ALLOCATE(delSO2_wet_h2o2(row_length,rows,model_levels))
        ALLOCATE(delSO2_wet_o3(row_length,rows,model_levels))
        ALLOCATE(delh2so4_chem(row_length,rows,model_levels))
        ALLOCATE(delSO2_drydep(row_length,rows,model_levels))
        ALLOCATE(delSO2_wetdep(row_length,rows,model_levels))
        IF (L_ukca_trophet) THEN 
          het_dimn=theta_field_size*model_levels 
        ELSE 
          het_dimn=1 
        END IF 
        
        IF (L_ukca_std_trop) THEN
           nhet_value = nhet_std_trop
        ELSE ! TROPISOP
           nhet_value = nhet_tropisop
        END IF
        ALLOCATE(het_rates(het_dimn,nhet_value))

! Initialise fluxes to zero as may not be filled by chemistry
        delSO2_wet_h2o2(:,:,:)=0.0
        delSO2_wet_o3(:,:,:)=0.0
        delh2so4_chem(:,:,:)=0.0
        delSO2_drydep(:,:,:)=0.0
        delSO2_wetdep(:,:,:)=0.0
      END IF

      IF (.NOT. ALLOCATED(strat_fluxdiags)) THEN 
         IF (.NOT. L_ukca_stratflux) n_strat_fluxdiags=0 
         IF (n_strat_fluxdiags > 0) THEN
            ALLOCATE(strat_fluxdiags(row_length,rows,model_levels,          & 
                 n_strat_fluxdiags)) 
            strat_fluxdiags=0.0 
         END IF
      END IF 

! Required in call to UKCA_MODE, but may be unallocated
      IF (.NOT. ALLOCATED(mode_diags)) THEN
        IF (.NOT. L_ukca_mode) n_mode_diags=1
        ALLOCATE(mode_diags(row_length,rows,model_levels,               &
                 n_mode_diags))
        mode_diags=0.0
      END IF

! Required in call to UKCA_EMISSION_CTL, but may be unallocated
      IF (.NOT. ALLOCATED(SO2_volc_3D)) THEN
        ALLOCATE(SO2_volc_3D(row_length,rows,model_levels))
      END IF
      IF (.NOT. ALLOCATED(SO2_biom_3D)) THEN
        ALLOCATE(SO2_biom_3D(row_length,rows,model_levels))
      END IF
      IF (.NOT. ALLOCATED(BC_biom_3D)) THEN
        ALLOCATE(BC_biom_3D(row_length,rows,model_levels))
      END IF
      IF (.NOT. ALLOCATED(OC_biom_3D)) THEN
        ALLOCATE(OC_biom_3D(row_length,rows,model_levels))
      END IF


      IF (first .AND. L_ukca_chem) THEN

        ALLOCATE(theta_latitude(row_length, rows))
        ALLOCATE(v_latitude(row_length, 0:rows))
        ALLOCATE(sinv_latitude(row_length, 0:rows))
        ALLOCATE(solid_angle(row_length, rows))
        ALLOCATE(volume(row_length,rows,model_levels))
        ALLOCATE(delta_r(row_length,rows,model_levels))
        ALLOCATE(area(row_length,rows,model_levels))
        ALLOCATE(mass(row_length, rows, model_levels))
        ALLOCATE(p_tropopause(row_length,rows))
        ALLOCATE(tropopause_level(row_length,rows))
        ALLOCATE(theta_trop(row_length,rows))
        ALLOCATE(pv_trop(row_length,rows))
        ALLOCATE(L_troposphere(row_length,rows,model_levels))

        IF ( l_vatpoles ) THEN
! For V-At-Poles, V latitude already bounds the grid-cells on both sides
! For PEs at N-Pole:  N_ROWS = ROWS + 1. For other PEs: N_ROWS = ROWS

        IF ( at_extremity (PNorth) )  THEN
          v_latitude(:,0:rows) =                            &
             ASIN(sin_v_latitude(:,1:n_rows))
        ELSE
          v_latitude(:,0:rows-1) =                           &
             ASIN(sin_v_latitude(:,1:n_rows))
          v_latitude(:,rows) = v_latitude(:,rows-1) + delta_phi
        END IF

        ELSE
        theta_latitude  = ASIN(sin_theta_latitude)
        v_latitude(:,0:rows-1) = theta_latitude - 0.5*delta_phi
        v_latitude(:,rows) = v_latitude(:,rows-1) + delta_phi
        END IF  ! vatpoles 

        WHERE (v_latitude < -0.5*pi) v_latitude = -0.5*pi
        WHERE (v_latitude >  0.5*pi) v_latitude =  0.5*pi
        sinv_latitude = SIN(v_latitude)
                     
        solid_angle = delta_lambda *                                    & 
            (sinv_latitude(:,1:rows) - sinv_latitude(:,0:rows-1)) 

        DO k=2,model_levels-1
          delta_r(:,:,k) = (r_rho_levels(1:row_length,1:rows,k+1) -     &
            r_rho_levels(1:row_length,1:rows,k))
        END DO
        delta_r(:,:,1) = (r_rho_levels(1:row_length,1:rows,2) -         &
          r_theta_levels(1:row_length,1:rows,0))
        delta_r(:,:,model_levels) =                                     &
          (r_theta_levels(1:row_length,1:rows,model_levels) -           &
           r_rho_levels(1:row_length,1:rows,model_levels))

        DO k=1,model_levels
          volume(:,:,k) = r_theta_levels(1:row_length,1:rows,k) *       &
            r_theta_levels(1:row_length,1:rows,k) *                     &
            delta_r(:,:,k) * delta_phi * delta_lambda *                 &
            FV_cos_theta_latitude(1:row_length,1:rows)
          area(1:row_length,1:rows,k) = 2*pi*                           &
            r_theta_levels(1:row_length,1:rows,k)**2 *                  &
            (sinv_latitude(1:row_length,1:rows) -                       &
             sinv_latitude(1:row_length,0:rows-1))/                     &
             global_row_length
        END DO

        IF (ALLOCATED(theta_latitude)) DEALLOCATE(theta_latitude)
        IF (ALLOCATED(v_latitude)) DEALLOCATE(v_latitude)
        IF (ALLOCATED(sinv_latitude)) DEALLOCATE(sinv_latitude)

        ALLOCATE(land_albedo_all(row_length,rows,Radsteps_per_day))
        land_albedo_all=min_surf_albedo

        ALLOCATE(z_half(row_length,rows,bl_levels))
        ALLOCATE(z_half_alllevs(row_length,rows,model_levels))
       
        DO k = 1,model_levels
          DO j = 1, rows
            DO i= 1, row_length
              z_half_alllevs(i,j,k)=                                    &
                           r_rho_levels(i,j,k)-r_theta_levels(i,j,0)
            END DO
          END DO
        END DO
        DO k = 1,bl_levels
          z_half(:,:,k) = z_half_alllevs(:,:,k)
        END DO
      
      END IF     ! first and L_ukca_chem

! Update land_albedo_all
!  index2 and index3 are the slices of the albedo array which 
!  we use for the interpolation in time 
!  index 2 should be 1 for timesteps 1 to a_sw_radstep_prog  
!  index 2 should be 2 for timesteps a_sw_radstep_prog+1 etc 

      index2 = 1+MOD((timestep_number-1)/a_sw_radstep_prog,Radsteps_per_day) 

! index3 is index2+1 unless index2 = Radsteps_per_day when it is 1 
      index3 = 1+MOD(index2,Radsteps_per_day) 

! Note that L_SW_Radiate is false on first timestep which is not what 
!  we want. Code our own test 
      IF (MOD((timestep_number-1),a_sw_radstep_prog) == 0) THEN  ! Update land_albedo_all 
        IF ( .NOT. ALLOCATED(land_albedo) ) THEN
! If coastal tiling is off then we do not have the land_albedo array allocated
          I1=LBOUND(tot_surf_sw,1)
          I2=UBOUND(tot_surf_sw,1)
          J1=LBOUND(tot_surf_sw,2)
          J2=UBOUND(tot_surf_sw,2)
          ALLOCATE(land_albedo(I1:I2,J1:J2))
        END IF
        WHERE (tot_surf_sw > 0.0 .AND. tot_surf_sw < 1.5e3)
          land_albedo=(tot_surf_sw-net_surf_sw)/tot_surf_sw
        ELSEWHERE
          land_albedo=min_surf_albedo
        ENDWHERE
        land_albedo_all(:,:,index2)=land_albedo(:,:)
      END IF

! Interpolate surf_albedo to timestep value

      IF (L_ukca_chem) THEN
!      fx is the interpolation weight for index2 
!      fx=1 for timesteps 1, 1+a_sw_radstep_prog etc 
!      fx=1/12 for timesteps a_sw_radstep_prog-1, 2*a_sw_radstep_prog-1 etc 
 
      fx=1.-REAL(MOD((timestep_number-1),a_sw_radstep_prog))/REAL(a_sw_radstep_prog)
      surf_albedo(:,:)=fx*land_albedo_all(:,:,index2)+                  &
                       (1.0-fx)*land_albedo_all(:,:,index3)
      IF (ALLOCATED(net_surf_sw)) DEALLOCATE(net_surf_sw)
      IF (ALLOCATED(tot_surf_sw)) DEALLOCATE(tot_surf_sw)

! read and interpolate aerosol surface area density
      IF (do_chemistry .AND. L_ukca_het_psc) THEN 
         IF (L_ukca_sa_clim) THEN
! DEPENDS ON: ukca_read_aerosol
            CALL ukca_read_aerosol(i_year, i_month, row_length, rows,   &
                 model_levels,                                          &
                 sin_theta_latitude(1:row_length, 1:rows),              &
                 eta_theta_levels(1:model_levels) * z_top_of_model,     &
                 L_ukca_use_background_aerosol, so4_sa)
! Below uses SO4 tracer from classic scheme for troposphere only
            IF (ALLOCATED(so4_aitken)) THEN
! Calculate surface area density fields from aerosol tracers
               DO k=1,model_levels
                  IF (eta_theta_levels(k) * z_top_of_model < 12000.) THEN 
                     so4_sa(:,:,k) = (                                  &
                      lambda_aitken * so4_aitken(1:row_length,1:rows,k) &
                     +lambda_accum  * so4_accum (1:row_length,1:rows,k))&
                     *rho_r2(1:row_length,1:rows,k)/                    &
                      (r_theta_levels(1:row_length,1:rows,k)**2)
                  END IF
               END DO
            ELSE
               so4_sa = 0.
               IF (first) THEN
                  errcode = -1 
                  cmessage=' Sulphate surface area set to zero'
                  CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
               END IF
            END IF
         ELSE 
           IF (L_ukca_mode ) THEN                             ! Take aerosol from MODE
             IF (PrintStatus >=  PrStatus_Diag) WRITE(6,*)              &
                'UKCA_MAIN: Taking aerosol surface area from MODE Diag.'
             so4_sa(1:row_length,1:rows,1:model_levels) =               &
                     chem_diags(:,:,:,icd_surfarea)
           ELSE
             cmessage=' Sulphate surface area undefined'
             errcode = 1
             CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
           END IF ! l_ukca_mode
         END IF ! L_ukca_sa_clim
      END IF

      t_theta_levels(:,:,1:model_levels)=                              &
           exner_theta_levels(:,:,1:model_levels) *                    &
           theta(:,:,1:model_levels)
      Thick_bl_levels(1:row_length,1:rows,1) = 2.0*                    &
            (r_theta_levels(1:row_length,1:rows,1) -                   &
             r_theta_levels(1:row_length,1:rows,0))
      DO k=2,bl_levels
        Thick_bl_levels(1:row_length,1:rows,k) =                       &
              r_rho_levels(1:row_length,1:rows,k+1) -                  &
              r_rho_levels(1:row_length,1:rows,k)
      END DO
      p_layer_boundaries(1:row_length,1:rows,0)                =       &
              pstar(1:row_length,1:rows)
      p_layer_boundaries(1:row_length,1:rows,model_levels)     =       &
              p_theta_levels(1:row_length,1:rows,model_levels)
      p_layer_boundaries(1:row_length,1:rows,1:model_levels-1) =       &
              p_rho_levels(1:row_length,1:rows,2:model_levels)


! Mass:
! We make the hydrostatic assumption in the diagnostic mass calculation
! here so that mass = -b*(solid_angle/(3*g))*(r_top^3 - r_bottom^3)
! where   b =                                                           &
!    &    (p_layer_boundaries(:,:,k) - p_layer_boundaries(:,:,k-1))/    &
!    &    (r_theta_levels(:,:,k)     - r_theta_levels(:,:,k-1))         &
      
         DO k=1,model_levels
           DO j=1,rows
            DO i=1,row_length
              mass(i,j,k) = (-solid_angle(i,j) / (3.*g)) *              &
            ((p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))/ &
             (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) ) *      &
             (r_theta_levels(i,j,k)**3. - r_theta_levels(i,j,k-1)**3.)
            END DO
          END DO
         END DO

        DO k=1,model_levels
          DO j=1,rows
           DO i=1,row_length
             totnodens(i,j,k) = p_theta_levels(i,j,k)                   &
                         /(zboltz*t_theta_levels(i,j,k))  ! molec/m3
             conv_ppn3d(i,j,k) = conv_rain3d(i,j,k) + conv_snow3d(i,j,k)
             ls_ppn3d(i,j,k)   = ls_rain3d  (i,j,k) + ls_snow3d  (i,j,k)

             conv_ppn3d(i,j,k) = MAX(conv_ppn3d(i,j,k), 0.0)
             ls_ppn3d  (i,j,k) = MAX(ls_ppn3d  (i,j,k), 0.0)
            END DO
         END DO
      END DO
      
      i1=LBOUND(rel_humid_frac,dim=1)
      j1=LBOUND(rel_humid_frac,dim=2)
      DO k = 1, wet_levels
! DEPENDS ON: qsat
          CALL QSAT(rel_humid_frac(i1,j1,k),t_theta_levels(i1,j1,k),    &
                   p_theta_levels(i1,j1,k),SIZE(t_theta_levels(:,:,k)))

          qsvp(:,:,k)=rel_humid_frac(1:row_length,1:rows,k)

! Convert to relative humidity fraction
          rel_humid_frac(1:row_length,1:rows,k) =                       &
                  q(1:row_length,1:rows,k)/                             &
                  rel_humid_frac(1:row_length,1:rows,k)
        END DO
        WHERE (rel_humid_frac < 0.0)
          rel_humid_frac=0.0           ! remove negatives
        ENDWHERE
        WHERE (rel_humid_frac >= 1.0)
          rel_humid_frac=0.999         ! remove values >= 1.0
        ENDWHERE
      END IF    ! L_ukca_chem


! Land fraction - if not already allocated set to 1.
! (If coastal tiling is off all land points are just land)
      IF (.NOT. ALLOCATED(fland) ) THEN
        ALLOCATE(fland(land_points))
        fland(1:land_points)=1.
      END IF

      IF (L_ukca_chem .OR. L_ukca_mode) THEN
! Set up land fraction (required for MODE and DRY DEPN)
      land_fraction(:,:)=0.0
      DO l = 1, land_points
        j = (land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        land_fraction(i,j) = fland(l)
      END DO

! Debug these every TS only  IF (PrintStatus >= PrStatus_Diag) 
        IF (PrintStatus >= PrStatus_Diag) THEN
          DO j=1,n_use_tracers
           WRITE(6,*) 'troptracer: ',mype,j,                           &
             MAXVAL(all_tracers(:,:,:,j)),MINVAL(all_tracers(:,:,:,j))
          END DO
          DO j=1,n_chem_diags
            WRITE(6,*) 'chem_diags: ',mype, j,                         &
             MAXVAL(chem_diags(:,:,:,j)),MINVAL(chem_diags(:,:,:,j))
          END DO
        END IF
! This is a huge code block - consider moving it to it's own subroutine
        IF (first .AND. printstatus >= prstatus_oper) THEN
          WRITE(6,*) ' ==================================='
          WRITE(6,*) '  MAX and MIN of UKCA INPUTS from D1'
          WRITE(6,*) ' ==================================='
          WRITE(6,*) 'Soil_Layer_Moisture: ',mype,                     &
            MAXVAL(soil_layer_moisture),MINVAL(soil_layer_moisture)
          DO j=1,n_chem_emissions
            WRITE(6,*) 'emission: ',mype,j,                            &
             MAXVAL(all_emissions(:,:,j)),MINVAL(all_emissions(:,:,j))
          END DO

          IF (n_strat_fluxdiags > 0) THEN
            DO j=1,n_strat_fluxdiags
              WRITE(6,*) 'strat_fluxdiag: ',mype, j,                     &
                          MAXVAL(strat_fluxdiags(:,:,:,j)),              &
                          MINVAL(strat_fluxdiags(:,:,:,j))
            END DO
          END IF

          IF (L_ukca_dust) THEN
            WRITE(6,*) 'frac_types: ',mype,                              &
                MAXVAL(frac_types),MINVAL(frac_types)
            WRITE(6,*) 'tstar_tile: ',mype,                              &
                MAXVAL(tstar_tile),MINVAL(tstar_tile)
            WRITE(6,*) 'snow_tile: ',mype,MAXVAL(snow_tile),             &
                MINVAL(snow_tile)
          END IF

          WRITE(6,*) 'Rough_length: ',mype,                              &
              MAXVAL(Rough_length),MINVAL(Rough_length)
          WRITE(6,*) 'Thick_bl_levels (level 1): ',mype,                 &
                MAXVAL(Thick_bl_levels(:,:,1)),                          &
                MINVAL(Thick_bl_levels(:,:,1))
          WRITE(6,*) 'conv_cloud_base: ',mype,                           &
                MAXVAL(conv_cloud_base(:,:)),                            &
                MINVAL(conv_cloud_base(:,:))
          WRITE(6,*) 'conv_cloud_top: ',mype,                            &
                MAXVAL(conv_cloud_top(:,:)),                             &
                MINVAL(conv_cloud_top(:,:))
          WRITE(6,*) 'land_sea_mask(1,1): ',land_sea_mask(1,1)
          WRITE(6,*) 'f3_at_u: ',mype,                                   &
                MAXVAL(f3_at_u(:,:)),                                    &
                MINVAL(f3_at_u(:,:))
          WRITE(6,*) 'FV_cos_theta_latitude: ',mype,                     &
                MAXVAL(FV_cos_theta_latitude(:,:)),                      &
                MINVAL(FV_cos_theta_latitude(:,:))

          DO k=1,model_levels
            WRITE(6,*) 'LEVEL: ',k
            WRITE(6,*) 'p_rho_levels:   ',mype,k,                        &
                        MAXVAL(p_rho_levels(1:row_length,1:rows,k)),     &
                        MINVAL(p_rho_levels(1:row_length,1:rows,k))
            WRITE(6,*) 'p_theta_levels: ',mype,k,                        &
                        MAXVAL(p_theta_levels(1:row_length,1:rows,k)),   &
                        MINVAL(p_theta_levels(1:row_length,1:rows,k))
            WRITE(6,*) 't_theta_levels: ',mype,k,                        &
                        MAXVAL(t_theta_levels(1:row_length,1:rows,k)),   &
                        MINVAL(t_theta_levels(1:row_length,1:rows,k))
            WRITE(6,*) 'rho_r2:         ',mype,k,                        &
                        MAXVAL(rho_r2(1:row_length,1:rows,k)),           &
                        MINVAL(rho_r2(1:row_length,1:rows,k))
            WRITE(6,*) 'ls_rain3d:      ',mype,k,                        &
                        MAXVAL(ls_rain3d(1:row_length,1:rows,k)),        &
                        MINVAL(ls_rain3d(1:row_length,1:rows,k))
              WRITE(6,*) 'ls_snow3d:      ',mype,k,                      &
                        MAXVAL(ls_snow3d(1:row_length,1:rows,k)),        &
                        MINVAL(ls_snow3d(1:row_length,1:rows,k))
              WRITE(6,*) 'conv_rain3d:      ',mype,k,                    &
                        MAXVAL(conv_rain3d(1:row_length,1:rows,k)),      &
                        MINVAL(conv_rain3d(1:row_length,1:rows,k))
              WRITE(6,*) 'conv_snow3d:      ',mype,k,                    &
                        MAXVAL(conv_snow3d(1:row_length,1:rows,k)),      &
                        MINVAL(conv_snow3d(1:row_length,1:rows,k))
              WRITE(6,*) 'aircraftems:      ',mype,k,                    &
                          MAXVAL(aircraftems(1:row_length,1:rows,k)),    &
                          MINVAL(aircraftems(1:row_length,1:rows,k))
              WRITE(6,*) 'PV_on_theta_mlevs:      ',mype,k,              &
                       MAXVAL(PV_on_theta_mlevs(1:row_length,1:rows,k)), &
                       MINVAL(PV_on_theta_mlevs(1:row_length,1:rows,k))

            IF (L_ukca_aerchem .OR. L_ukca_achem) THEN
              WRITE(6,*) 'Volcanic SO2 emissions: ',mype,k,              &
                        MAXVAL(SO2_volc_3D(1:row_length,1:rows,k)),      &
                        MINVAL(SO2_volc_3D(1:row_length,1:rows,k))
              WRITE(6,*) 'Biomass BC emissions: ',mype,k,                &
                        MAXVAL(BC_biom_3D(1:row_length,1:rows,k)),       &
                        MINVAL(BC_biom_3D(1:row_length,1:rows,k))
              WRITE(6,*) 'Biomass OC emissions: ',mype,k,                &
                        MAXVAL(OC_biom_3D(1:row_length,1:rows,k)),       &
                        MINVAL(OC_biom_3D(1:row_length,1:rows,k))
            END IF
          END DO   ! k

          DO k=1,wet_levels
            WRITE(6,*) 'WET LEVEL: ',k
            WRITE(6,*) 'q:              ',mype,k,                      &
                      MAXVAL(q(1:row_length,1:rows,k)),                &
                      MINVAL(q(1:row_length,1:rows,k))
            WRITE(6,*) 'qcl:            ',mype,k,                      &
                      MAXVAL(qcl(1:row_length,1:rows,k)),              &
                      MINVAL(qcl(1:row_length,1:rows,k))
            WRITE(6,*) 'qcf:            ',mype,k,                      &
                      MAXVAL(qcf(1:row_length,1:rows,k)),              &
                        MINVAL(qcf(1:row_length,1:rows,k))
            WRITE(6,*) 'rel_humid_frac: ',mype,k,                      &
                      MAXVAL(rel_humid_frac(1:row_length,1:rows,k)),   &
                      MINVAL(rel_humid_frac(1:row_length,1:rows,k))
          END DO
          DO k=1,model_levels+1
            WRITE(6,*) 'EXNER_rho_levels: ',mype,k,                    &
                      MAXVAL(exner_rho_levels(1:row_length,1:rows,k)), &
                      MINVAL(exner_rho_levels(1:row_length,1:rows,k))
          END DO
          DO k=1,bl_levels
            WRITE(6,*) 'BL LEVEL: ',k
            WRITE(6,*) 'rhokh_mix:     ',mype,k,                       &
                    MAXVAL(rhokh_mix(1:row_length,1:rows,k)),          &
                    MINVAL(rhokh_mix(1:row_length,1:rows,k))
            WRITE(6,*) 'dtrdz_charney_grid:     ',mype,k,              &
                MAXVAL(dtrdz_charney_grid(1:row_length,1:rows,k)),     &
                MINVAL(dtrdz_charney_grid(1:row_length,1:rows,k))
            WRITE(6,*) 'z_half:     ',mype,k,                          &
                MAXVAL(z_half(1:row_length,1:rows,k)),                 &
                MINVAL(z_half(1:row_length,1:rows,k))
          END DO
          IF ((i_ukca_photol == i_ukca_fastj) .OR.                      &
              (i_ukca_photol == i_ukca_fastjx)) THEN
            DO k=1,model_levels
              WRITE(6,*) 'so4_aitken: ', mype,k,                       &
                       MAXVAL(so4_aitken(1:row_length,1:rows,k)),      &
                       MINVAL(so4_aitken(1:row_length,1:rows,k))
              WRITE(6,*) 'so4_accum: ', mype,k,                        &
                      MAXVAL(so4_accum(1:row_length,1:rows,k)),        &
                      MINVAL(so4_accum(1:row_length,1:rows,k))
            END DO
            IF ((i_ukca_photol == i_ukca_fastjx) .AND.                 &
              ALLOCATED(sulphate_od)) THEN
              DO k=1,sw_spectrum(1)%n_band
                WRITE(6,*) 'sulphate_od: ', mype,k,                    &
                        MAXVAL(sulphate_od(1:row_length,1:rows,k)),    &
                        MINVAL(sulphate_od(1:row_length,1:rows,k))
              END DO
            END IF
          END IF
          WRITE(6,*) 'kent:     ', mype,MAXVAL(kent),    MINVAL(kent)
          WRITE(6,*) 'kent_dsc: ', mype,MAXVAL(kent_dsc),MINVAL(kent_dsc)
          WRITE(6,*) 'zhsc:     ', mype,MAXVAL(zhsc),    MINVAL(zhsc)
          WRITE(6,*) 'ml_depth: ', mype,MAXVAL(ml_depth),MINVAL(ml_depth)
          DO k=1,3
            WRITE(6,*) 'we_lim: ', mype,k,                             &
                      MAXVAL(we_lim(1:row_length,1:rows,k)),           &
                      MINVAL(we_lim(1:row_length,1:rows,k))
            WRITE(6,*) 't_frac: ', mype,k,                             &
                      MAXVAL(t_frac(1:row_length,1:rows,k)),           &
                      MINVAL(t_frac(1:row_length,1:rows,k))
            WRITE(6,*) 'zrzi: ', mype,k,                               &
                      MAXVAL(zrzi(1:row_length,1:rows,k)),             &
                      MINVAL(zrzi(1:row_length,1:rows,k))
            WRITE(6,*) 'we_lim_dsc: ', mype,k,                         &
                      MAXVAL(we_lim_dsc(1:row_length,1:rows,k)),       &
                      MINVAL(we_lim_dsc(1:row_length,1:rows,k))
            WRITE(6,*) 't_frac_dsc: ', mype,k,                         &
                        MAXVAL(t_frac_dsc(1:row_length,1:rows,k)),     &
                      MINVAL(t_frac_dsc(1:row_length,1:rows,k))
            WRITE(6,*) 'zrzi_dsc: ', mype,k,                           &
                      MAXVAL(zrzi_dsc(1:row_length,1:rows,k)),         &
                      MINVAL(zrzi_dsc(1:row_length,1:rows,k))
          END DO
          WRITE(6,*) 'pstar: ',    mype,MAXVAL(pstar),   MINVAL(pstar)
          WRITE(6,*) 'Tstar: ',    mype,MAXVAL(Tstar),   MINVAL(Tstar)
          WRITE(6,*) 'fland: ',    mype,MAXVAL(fland),   MINVAL(fland)
          WRITE(6,*) 'u_s  : ',    mype,MAXVAL(u_s  ),   MINVAL(u_s  )

        END IF   ! End of IF (first .AND. printstatus >= prstatus_oper) statement

!       Warning statements about possible mismatch between interactive
!       methane emissions and emissions ancillary
        IF (first) THEN
          IF (L_ukca_qch4inter .AND. (.NOT. L_ukca_prescribech4)) THEN
            cmessage = 'CH4 WETLANDS EMS ARE ON - Ancillary SHOULD NOT contain wetland ems'       
            errcode=-1

            CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
          ELSE IF (.NOT. L_ukca_prescribech4) THEN
            cmessage = 'CH4 WETLANDS EMS ARE OFF - Ancillary SHOULD contain wetland ems'
            errcode=-2

            CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
          END IF
        END IF

      END IF   ! End of IF (first) statement

!     Set up tile info
      IF (L_ukca_chem) THEN
! DEPENDS ON: tilepts
        CALL TILEPTS(land_points,frac_types,tile_pts,tile_index)

!     Set tile fractions to 1 if aggregate tiles are used (NTILES=1).
!     Otherwise, set tile fractions to surface type fractions.

        tile_frac=0.0
        IF (ntiles == 1) THEN
          DO l=1,land_points
            tile_frac(l,1) = 1.
          END DO
        ELSE
          DO n=1,ntiles
            DO j=1,tile_pts(n)
              l = tile_index(j,n)
              tile_frac(l,n) = frac_types(l,n)
            END DO
          END DO
        END IF

      END IF    ! L_ukca_chem

!======================================
! AGE OF AIR:
!  n_age is set in ukca_calc_cspecies if chemistry is used (L_ukca_chem=T)
!  AGE OF AIR *WILL* BE USED WHEN CHEMISTRY IS USED
! Increment the Age of air tracer, and set lower levels to zero
! Set the Age of air number if no chemistry present
      IF (first .AND. L_ukca_ageair .AND. .NOT. L_ukca_chem) THEN
         age_loop: DO i=1,Nukca_D1items
            IF ((UkcaD1Codes(I)%item == 150 .AND.                           &
                 UkcaD1Codes(I)%section == UKCA_sect) .AND.                 & 
                 nm_spec(150) == 'AGE OF AIR') THEN
               n_age=i
               EXIT age_loop
            END IF
         END DO age_loop
      END IF
      IF (n_age > 0) THEN
        all_tracers(:,:,1:model_levels,n_age)=                          &
            all_tracers(:,:,1:model_levels,n_age) + timestep
        DO k = 1,max_zero_age
          all_tracers(:,:,k,n_age) = 0.0
        END DO
        ! enforce upper boundary condition
        all_tracers(:,:,model_levels,n_age)=&
             all_tracers(:,:,model_levels-1,n_age)
      END IF
!======================================



      IF (L_ukca_chem) THEN
! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA CHEMISTRY MODEL',5)

! Transform tracers to ensure elemental conservation
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
         (L_ukca_stratcfc))  THEN
! DEPENDS ON: ukca_transform_halogen
        CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                       model_levels, tdims_s%halo_i,tdims_s%halo_j,     &
                       all_tracers(:,:,1:model_levels,:),               &
                       qdims_s%halo_i,qdims_s%halo_j,                   &
                       q(:,:,1:model_levels), .TRUE.,timestep_number)
      END IF


!       Calculate tropopause pressure using a combined
!       theta and PV surface

        CALL UKCA_CALC_TROPOPAUSE(row_length, rows, model_levels,      &
             sin_theta_latitude(1:row_length,1:rows),                  &
             theta(1:row_length,1:rows,1:model_levels),                &
             pv_on_theta_mlevs(1:row_length,1:rows,1:model_levels),    &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
             p_theta_levels(1:row_length,1:rows,1:model_levels),       &
             p_tropopause(1:row_length,1:rows),                        &
             tropopause_level(1:row_length,1:rows))

        IF (L_asad_use_chem_diags .AND. L_asad_use_trop_mask) &
             CALL asad_tropospheric_mask( &
             rows,row_length,model_levels,ierr)



!       Calculate the difference in tracers per timestep in the
!       troposphere (moles/s) due to transport using the
!       trmol_post_atmstep array from this timestep and the
!       trmol_post_chem array from the previous timestep.
!       Do this only for those stratospheric flux diagnostics
!       which are switched on.

        IF (first .AND. n_strat_fluxdiags > 0) THEN

         ALLOCATE(trmol_post_chem(row_length,rows,model_levels,       &
                                    n_chem_tracers))  ! moles
          trmol_post_chem = 0.0

          DO l=1,n_strat_fluxdiags
            DO k=1,model_levels
             DO j=1,rows
                    DO i=1,row_length
                  strat_fluxdiags(i,j,k,l) = 0.0
               END DO
             END DO
           END DO
          END DO

       ELSE IF ((.NOT. first) .AND. n_strat_fluxdiags > 0) THEN

              ALLOCATE(trmol_post_atmstep(row_length,rows,model_levels,     &
                                    n_chem_tracers))  ! moles
         trmol_post_atmstep = 0.0

          DO l=1,n_chem_tracers
            DO k=1,model_levels
              DO j=1,rows
                DO i=1,row_length
                  trmol_post_atmstep(i,j,k,l) = all_tracers(i,j,k,l)    &
                      *totnodens(i,j,k)*volume(i,j,k)/(c_species(l)     &
                      *avogadro)   ! moles
                END DO
              END DO
            END DO
          END DO

          icnt = 0
          DO l=1,n_chem_tracers
            IF (UkcaD1codes(istrat_first+l-1)%item /= imdi) THEN
              icnt = icnt + 1

              DO k=1,model_levels
                DO j=1,rows
                  DO i=1,row_length
                    IF (L_troposphere(i,j,k)) THEN              ! troposphere
                      strat_fluxdiags(i,j,k,icnt) =                     &
                       (trmol_post_atmstep(i,j,k,l)-                    &
                        trmol_post_chem(i,j,k,l))/timestep      ! moles/sec
                    ELSE                                        ! stratosphere
                      strat_fluxdiags(i,j,k,icnt) = 0.0
                    END IF
                  END DO
                END DO
              END DO
            END IF
          END DO     ! n_chem_tracers

        IF (ALLOCATED(trmol_post_atmstep)) DEALLOCATE(trmol_post_atmstep)  ! moles

        END IF      ! End of IF first and n_strat_fluxdiags statement

        IF (.NOT. first) THEN ! DO NOT CALL ON FIRST TIMESTEP
           ! ASAD Flux Diags STE
           IF (L_asad_use_chem_diags .AND. L_asad_use_STE) &
                CALL asad_tendency_STE( &
                row_length,rows,model_levels, &
                n_chem_tracers,all_tracers, &
                totnodens,volume,mass,timestep, &
                calculate_STE,.FALSE.,ierr)
        END IF

        IF (do_chemistry .AND. (L_asad_use_chem_diags &
             .AND. L_asad_use_tendency)) &
             CALL asad_tendency_STE( &
             row_length,rows,model_levels, &
             n_chem_tracers,all_tracers, &
             totnodens,volume,mass,timestep, &
             calculate_tendency,.TRUE.,ierr)



!       Do pre-processing for photolysis scheme

! The following (FAST-J or Fast-jX calculation) only needs to be done when
! chemistry is called:

        IF (do_chemistry) THEN

          ALLOCATE(um_ozone3d(row_length,rows,ozone_levels))
          IF (SIZE(um_ozone,DIM=1) == 1) THEN   ! Use zonal average
            DO i=1,row_length
              um_ozone3d(i,:,:)=um_ozone(1,:,1:ozone_levels)
            END DO
          ELSE                                  ! 3-D ozone specifed
            um_ozone3d(:,:,:)= um_ozone(:,:,1:ozone_levels)
          END IF

!       Convert 2D or 3D ozone field from UM to a 1D field

          nd_o3=SIZE(um_ozone)
          ALLOCATE(um_ozone1D(nd_o3))
          um_ozone1D = RESHAPE(um_ozone,(/nd_o3/))

!       Set the MOSES II flag if 8A boundary layer selected

          IF(h_sect(3)=="08A".OR. h_sect(3)=="08B") THEN
            l_moses_ii = .TRUE.
          ELSE
            l_moses_ii = .FALSE.
          END IF

          ALLOCATE(dj(row_length,rows,model_levels,jppj))
          dj=0.0

          IF ((i_ukca_photol == i_ukca_fastj)) THEN
! DEPENDS ON: timer
            IF (ltimer) CALL TIMER('UKCA FASTJ MODEL  ',5)

! DEPENDS ON: ukca_fastj
            CALL UKCA_FASTJ(                                           &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
            dj(1:row_length,1:rows,1:model_levels,1:jppj),             &
            p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
            pstar(1:row_length,1:rows),                                &
            p_theta_levels(1:row_length,1:rows,1:model_levels),        &
            t_theta_levels(1:row_length,1:rows,1:model_levels),        &
            tstar(1:row_length,1:rows),                                &
            so4_aitken(1:row_length,1:rows,1:model_levels),            &
            so4_accum(1:row_length,1:rows,1:model_levels),             &
            q(1:row_length,1:rows,1:wet_levels),                       &
            qcl(1:row_length,1:rows,1:wet_levels),                     &
            qcf(1:row_length,1:rows,1:wet_levels),                     &
            area_cloud_fraction(1:row_length,1:rows,1:wet_levels),     &
            conv_cloud_lwp(1:row_length,1:rows),                       &
            conv_cloud_top(1:row_length,1:rows),                       &
            conv_cloud_base(1:row_length,1:rows),                      &
            conv_cloud_amount(1:row_length,1:rows,1:wet_levels),       &
            surf_albedo(1:row_length,1:rows),                          &
            nd_o3, um_ozone3d(1:row_length,1:rows,1:ozone_levels),     &
            fland,                                                     &
            a_realhd(rh_z_top_theta),                                  &
            l_moses_ii)

! DEPENDS ON: timer
            IF (ltimer) CALL TIMER('UKCA FASTJ MODEL  ',6)

          END IF    ! i_ukca_fastj

          IF ((i_ukca_photol == i_ukca_fastjx)) THEN
! DEPENDS ON: timer
            IF (ltimer) CALL TIMER('UKCA FASTJX MODEL  ',5)

! DEPENDS ON: ukca_fastjx
            CALL UKCA_FASTJX(                                          &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
            dj(1:row_length,1:rows,1:model_levels,1:jppj),             &
            p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
            pstar(1:row_length,1:rows),                                &
            p_theta_levels(1:row_length,1:rows,1:model_levels),        &
            t_theta_levels(1:row_length,1:rows,1:model_levels),        &
            rho_r2(1:row_length,1:rows,1:model_levels),                &
            z_top_of_model,                                            &
            so4_aitken(1:row_length,1:rows,1:model_levels),            &
            so4_accum(1:row_length,1:rows,1:model_levels),             &
            qcl(1:row_length,1:rows,1:wet_levels),                     &
            qcf(1:row_length,1:rows,1:wet_levels),                     &
            rel_humid_frac(1:row_length,1:rows,:),                     &
            area_cloud_fraction(1:row_length,1:rows,1:wet_levels),     &
            conv_cloud_lwp(1:row_length,1:rows),                       &
            conv_cloud_top(1:row_length,1:rows),                       &
            conv_cloud_base(1:row_length,1:rows),                      &
            conv_cloud_amount(1:row_length,1:rows,1:wet_levels),       &
            surf_albedo(1:row_length,1:rows),                          &
            tropopause_level(1:row_length, 1:rows),                    &
            all_tracers(1:row_length,1:rows,1:model_levels,n_o3),      &
            land_fraction)

! DEPENDS ON: timer
            IF (ltimer) CALL TIMER('UKCA FASTJX MODEL  ',6)
          END IF  ! End of IF (i_ukca_fastjx) statement

!       Clear workspace

         IF (ALLOCATED(um_ozone1D)) DEALLOCATE(um_ozone1D)
         IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)

        END IF   ! do_chemistry

!       Allocate array for wetland ch4 emissions when
!       L_ukca_qch4inter is false. Set emissions to zero.

        IF (.NOT.(L_ukca_qch4inter)) THEN
         ALLOCATE(ch4_wetl_emiss(1:row_length,1:rows))
         ch4_wetl_emiss(1:row_length,1:rows) = 0.0
        END IF

! Emission part. Before calling emissions,
! set up array of lower-boundary mixing ratios for stratospheric species.
! Note that some mixing ratios (CH4, N2O, F11, F12, F113, H22) are defined
! in the RUN_UKCA namelist held in UKCA module ukca_option_mod through the 
! UMUI, others are not. If for some species
! actual emissions (not prescribed lower boundary conditions) are required,
! remove these species from lbc_spec, e.g, for CH4.
! Mass mixing ratios are best estimates, taken from the WMO Scientific
! assessment of ozone depletion (2002), for 2000 conditions.
! H2 corresponds to 0.5 ppmv.
! Select lower boundary MMRs. Either from UMUI or predefined constants
! Note that for bromine compounds (MeBr, H1211, H1301) are increased
! by 24% to account for non-included species and HBr is set to 0.
!
        IF (L_ukca_strat.OR. L_ukca_stratcfc .OR. L_ukca_strattrop) THEN

           IF (L_ukca_useumuivals) THEN
! DEPENDS ON: ukca_scenario_prescribed
              CALL UKCA_SCENARIO_PRESCRIBED(n_boundary_vals, lbc_spec,  &
                   lbc_mmr, .NOT.(L_ukca_stratcfc))
           ELSE
! DEPENDS ON: ukca_scenario_wmoa1
              CALL UKCA_SCENARIO_WMOA1(n_boundary_vals, lbc_spec,       &
                   i_year, i_day_number, lbc_mmr, &
                   .NOT.(L_ukca_stratcfc))
           END IF

           IF ((first .AND. L_ukca_useumuivals) .AND.                   &
                (printstatus >= prstatus_oper)) THEN
              DO ii=1,n_boundary_vals
                 WRITE(6,'(2A15,A3,E12.3)') ' UKCA Lower BC ',          &
                      TRIM(ADJUSTL(lbc_spec(ii))),' = ',lbc_mmr(ii)
              END DO
           END IF

        END IF    ! L_ukca_strat etc

!        IF (do_chemistry) THEN

!       Calculate solar zenith angle

          ALLOCATE(cos_zenith_angle(row_length, rows))
          IF (.NOT. ALLOCATED(int_zenith_angle))  &
             ALLOCATE(int_zenith_angle(row_length, rows))


! DEPENDS ON: solpos
          CALL SOLPOS (PREVIOUS_TIME(7), PREVIOUS_TIME(1),              &
               LCAL360, L_SEC_VAR, L_EqT, Eq_Time, Sindec, SCS)

          secondssincemidnight = REAL(Previous_time(4)*3600             &
                               +      Previous_time(5)*60               &
                               +      Previous_time(6))

! Compute COS(sza) integral only when day changes
          IF (FIRST) n_day = 0
          IF ((n_day /= i_day)) THEN
               n_day = i_day
               int_zenith_angle(:,:) = 0.0
               DO kk = 1, (A_Steps_per_hr*24)
                  ssmn_incr = secondssincemidnight + timestep*(kk-1)
                  CALL UKCA_SOLANG(Sindec, ssmn_incr,                   &
                                   timestep, Eq_Time,                   &
                                   sin_theta_latitude, true_longitude,  &
                                   theta_field_size, cos_zenith_angle )

                  WHERE (cos_zenith_angle > 0.0)
                    int_zenith_angle(:,:) = int_zenith_angle(:,:) +     &
                       cos_zenith_angle(:,:) * timestep
                  ENDWHERE
               END DO
          END IF

! This call needs to be after above loop to prevent errors with value in
!  cos_zenith_angle
! DEPENDS ON: ukca_solang
          CALL UKCA_SOLANG(Sindec, secondssincemidnight,                &
                      timestep, Eq_Time,                                &
                      sin_theta_latitude, true_longitude,               &
                      theta_field_size, cos_zenith_angle )


! DEPENDS ON: ukca_emission_ctl
        CALL UKCA_EMISSION_CTL(                                        &
           n_chem_tracers+n_aero_tracers, n_mode_tracers,              &
           n_use_emissions, n_chem_emissions, timestep,                &
           em_chem_spec(1:n_use_emissions),                            &
           f3_at_u(1:row_length,1:rows),                               &
           r_rho_levels(1:row_length,1:rows,1:model_levels),           &
           r_theta_levels(1:row_length,1:rows,0:model_levels),         &
           sin_theta_latitude(1:row_length,1:rows),                    &
           FV_cos_theta_latitude(1:row_length,1:rows),                 &
           tan_theta_latitude(1:row_length,1:rows),                    &
           cos_zenith_angle(1:row_length,1:rows),                      &
           int_zenith_angle(1:row_length,1:rows),                      &
           true_longitude, delta_lambda, delta_phi,                    &
           i_year, i_month, i_day,                                     &
           tropopause_height(1:row_length,1:rows),                     &
           land_sea_mask(1:row_length,1:rows),                         &
           conv_cloud_base(1:row_length,1:rows),                       &
           conv_cloud_top(1:row_length,1:rows),                        &
           theta(1:row_length,1:rows,1:model_levels),                  &
           q(1:row_length,1:rows,1:wet_levels),                        &
           qcl(1:row_length,1:rows,1:wet_levels),                      &
           qcf(1:row_length,1:rows,1:wet_levels),                      &
           exner_rho_levels(1:row_length,1:rows,1:model_levels+1),     &
           rho_r2(1:row_length,1:rows,1:model_levels),                 &
           p_layer_boundaries(1:row_length,1:rows,0:model_levels),     &
           p_theta_levels(1:row_length,1:rows,1:model_levels),         &
           t_theta_levels(1:row_length,1:rows,1:model_levels),         &
           all_emissions(1:row_length,1:rows,1:n_chem_emissions),      &
           aircraftems(1:row_length,1:rows,1:model_levels),            &
           em_index,                                                   &
           BC_biom_3D(1:row_length,1:rows,1:model_levels),             &
           OC_biom_3D(1:row_length,1:rows,1:model_levels),             &
           SO2_volc_3D(1:row_length,1:rows,1:model_levels),            &
           dust_flux, U_scalar_10m,                                    &
           rough_length,                                               &
           land_fraction,                                              &
           seaice_frac,                                                &
           area,                                                       &
           z_half, alpha_cd(1:bl_levels), ml_depth,                    &
           rhokh_mix,                                                  &
           dtrdz_charney_grid, kent, we_lim(:,:,1:3),                  &
           t_frac(:,:,1:3), zrzi(:,:,1:3), kent_dsc,                   &
           we_lim_dsc(:,:,1:3), t_frac_dsc(:,:,1:3),                   &
           zrzi_dsc(:,:,1:3), zhsc,                                    &
           ch4_wetl_emiss(1:row_length,1:rows),                        &
           all_tracers(1:row_length,1:rows,1:model_levels,             &
                       1:n_chem_tracers+n_aero_tracers),               &
           all_tracers(1:row_length,1:rows,:,n_chem_tracers+           & 
                       n_aero_tracers+1:n_chem_tracers+n_aero_tracers+ & 
                       n_mode_tracers),                                &
           n_boundary_vals, lbc_spec, lbc_mmr,                         &
           mass, totnodens, volume,                                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
           STASH_maxlen(38,A_im),STASHwork38)


        IF (L_asad_use_chem_diags .AND. L_asad_use_mass_diagnostic) &
             CALL asad_mass_diagnostic( &
             row_length, &
             rows, &
             model_levels, &
             mass, &
             ierr)


! Do chemistry calculation here (at chemistry timesteps)
        IF (do_chemistry) THEN

          IF (.NOT. ALLOCATED(T_chem)) & 
               ALLOCATE(T_chem(row_length,rows,model_levels))
          IF (.NOT. ALLOCATED(Q_chem)) & 
               ALLOCATE(Q_chem(row_length,rows,model_levels))
          T_chem(:,:,:) = t_theta_levels(1:row_length,1:rows,:)
          Q_chem(:,:,:) = q(1:row_length,1:rows,:)
 
!======================================
! PASSIVE OZONE - CheST/StratChem ONLY
!  n_passive is set in ukca_calc_cspecies
         ! Do passive ozone
         IF (n_passive > 0) THEN
            ! assumes 360 day calendar
            IF((i_day_number == 91).OR.(i_day_number == 301)) THEN
               IF (i_hour == 0) THEN
                  all_tracers(:,:,1:model_levels,n_passive) =          &
                        all_tracers(:,:,1:model_levels,n_o3)
               END IF
            END IF
         END IF
!======================================

        IF (L_asad_use_chem_diags .AND. L_asad_use_tendency) &
             CALL asad_tendency_STE( &
             row_length,rows,model_levels, &
             n_chem_tracers,all_tracers, &
             totnodens,volume,mass,timestep, &
             calculate_tendency,.TRUE.,ierr)


        IF (interval > 1) THEN
           ! MAY NEED TO CALL A SECOND TIME TO HAVE CORRECT MEANING PERIOD
! DEPENDS ON: ukca_solang
           CALL UKCA_SOLANG(Sindec, secondssincemidnight,              &
                REAL(interval)*timestep, Eq_Time,                      &
                sin_theta_latitude, true_longitude,                    &
                theta_field_size, cos_zenith_angle )
        END IF

        IF ( firstchem .AND. L_ukca_mode) THEN
          IF (L_ukca_aerchem) THEN
! SO2 wet oxidation and H2SO4 updating done in UKCA 
            wetox_in_aer = 0
            uph2so4inaer = 0
          ELSE IF (L_ukca_achem) THEN
! SO2 wet oxidation in UKCA, and H2SO4 updating done in MODE
            wetox_in_aer = 0
            uph2so4inaer = 1
          ELSE                        ! No SO2 oxidation in UKCA chemistry
            wetox_in_aer = 1          ! Sulphur emissions are still needed
            uph2so4inaer = 1
          END IF
          IF (printstatus >= prstatus_oper) THEN
            WRITE(6,*) 'uph2so4inaer is set to: ',uph2so4inaer
            WRITE(6,*) 'wetox_in_aer is set to: ',wetox_in_aer
          END IF
        ELSE
          wetox_in_aer = 0
          uph2so4inaer = 0
        END IF                 ! firstchem etc
 
! DEPENDS ON: ukca_chemistry_ctl
         CALL UKCA_CHEMISTRY_CTL(I_month, I_day_number, I_hour,        &
             I_minute - INT(timestep)/60,                              &
             REAL(interval)*timestep,                                  &
             n_chem_tracers+n_aero_tracers, n_chem_diags,              &
             f3_at_u(1:row_length,1:rows)/two_omega,                   &
             FV_cos_theta_latitude(1:row_length,1:rows),               &
             true_longitude,                                           &
             p_theta_levels(1:row_length,1:rows,1:model_levels),       &
             T_chem(1:row_length,1:rows,1:model_levels),               &
             Q_chem(1:row_length,1:rows,1:model_levels),               &
             qcf(1:row_length,1:rows,1:model_levels),                  &
             qcl(1:row_length, 1:rows,1:model_levels),                 &
             rel_humid_frac(1:row_length,1:rows,1:model_levels),       &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
             r_theta_levels(1:row_length, 1:rows, :),                  &
             z_top_of_model,                                           &
             cos_zenith_angle,                                         &
             all_tracers(1:row_length,1:rows,1:model_levels,           &
                         1:n_chem_tracers+n_aero_tracers),             &
             chem_diags(1:row_length,1:rows,:,:),                      &
             Tstar,                                                    &
             Thick_bl_levels,                                          &
             Rough_length,                                             &
             u_s,                                                      &
             ls_ppn3d, conv_ppn3d,                                     &
             cloud_frac(1:row_length, 1:rows, :),                      &
             dj(:,:,:,:),                                              &
             volume(:,:,:),                                            &
             mass(:,:,:),                                              &
! Extra variables for new dry dep scheme
             land_points, land_index,                                  &
             tile_pts, tile_index, tile_frac,                          &
             zbl, surf_hf, seaice_frac, stcon,                         &
             soil_layer_moisture(:,1), fland,                          &
             laift_lp, canhtft_lp,                                     &
             z0tile_lp, tstar_tile, canwctile_lp,                      &
             pv_on_theta_mlevs,                                        &
             theta(1:row_length,1:rows,1:model_levels),                &
             um_ozone3d,                                               &
             uph2so4inaer,                                             &
             delso2_wet_h2o2,                                          &
             delso2_wet_o3,                                            &
             delh2so4_chem,                                            &
             delso2_drydep,                                            &
             delso2_wetdep,                                            &
             so4_sa                                                    &
          )


          IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)
          IF (ALLOCATED(T_chem)) DEALLOCATE(T_chem)
          IF (ALLOCATED(Q_chem)) DEALLOCATE(Q_chem)

           firstchem = .FALSE.
         END IF ! do_chemistry

        IF (n_strat_fluxdiags > 0) THEN
         DO l=1,n_chem_tracers
           DO k=1,model_levels
             DO j=1,rows
               DO i=1,row_length
                 trmol_post_chem(i,j,k,l) = all_tracers(i,j,k,l)        &
                                     *totnodens(i,j,k)*volume(i,j,k)    &
                                     /(c_species(l)*avogadro)    ! moles
               END DO
             END DO
           END DO
         END DO
        END IF   ! End of IF n_strat_fluxdiags > 0 statement

!  Save all_tracers array for STE calculation in ASAD Flux Diagnostics
        IF (.NOT. ALLOCATED(fluxdiag_all_tracers) .AND. &
             (L_asad_use_chem_diags .AND. L_asad_use_STE)) &
             THEN
           ALLOCATE(fluxdiag_all_tracers(&
                1-offx:row_length+offx,1-offy:rows+offy, &
                1:model_levels,1:n_chem_tracers))
           fluxdiag_all_tracers(:,:,1:model_levels,:) =                  &
                all_tracers(:,:,1:model_levels,1:n_chem_tracers)
        END IF

        IF (do_chemistry .AND. (L_asad_use_chem_diags                   &
             .AND. L_asad_use_tendency))                                &
          CALL ASAD_TENDENCY_STE(row_length,rows,model_levels,          &
                                 n_chem_tracers,all_tracers,            &
                                 totnodens,volume,mass,timestep,        &
                                 calculate_tendency,.FALSE.,ierr)

! Transform halogen/nitrogen/hydrogen species back
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
       (L_ukca_stratcfc)) THEN
         ! first, copy unlumped NO2, BrO and HCl fields into chem_diags(1,2, & 3)
         ! all in VMR, converted to MMR in transform halogen for normal tracers,
         ! so convert again here
         chem_diags(1:row_length,1:rows,1:model_levels,icd_no2) =             &
              all_tracers(1:row_length,1:rows,1:model_levels,n_no2)
         chem_diags(1:row_length,1:rows,1:model_levels,icd_bro) =             &
              all_tracers(1:row_length,1:rows,1:model_levels,n_bro)
         chem_diags(1:row_length,1:rows,1:model_levels,icd_hcl) =             &
              all_tracers(1:row_length,1:rows,1:model_levels,n_hcl)
! DEPENDS ON: ukca_transform_halogen
        CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                        model_levels, tdims_s%halo_i, tdims_s%halo_j,   &
                        all_tracers(:,:,1:model_levels,:),              &
                        qdims_s%halo_i, qdims_s%halo_j,                 &
                        q(:,:,1:model_levels), .FALSE., 0)
        ! now put Cly into chem_diags = (lumped HCl)-(3*CFCl3+2*CF2Cl2)
        ! CONVERT THROUGH TO VMR FOR SUBTRACTION, THEN BACK TO MMR AGAIN
        chem_diags(1:row_length,1:rows,1:model_levels,icd_cly) =              &
              ((all_tracers(1:row_length,1:rows,1:model_levels,n_hcl)   &
              /c_hcl) - (3.0E0*(all_tracers(                            &
              1:row_length,1:rows,1:model_levels,n_cfcl3)/c_cfcl3) +    &
              2.0E0*all_tracers(                                        &
              1:row_length,1:rows,1:model_levels,n_cf2cl2)/c_cf2cl2))   &
              *c_hcl
      END IF

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA CHEMISTRY MODEL',6)
      END IF    ! End of IF (L_ukca_chem) statement

! MODE AEROSOL SCHEME
      IF(L_UKCA_MODE .AND. do_chemistry) THEN

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA AEROSOL MODEL  ',5)

! DEPENDS ON: ukca_aero_ctl
        CALL UKCA_AERO_CTL(i_month, i_day_number, i_hour,               &
             I_minute - INT(timestep)/60,                               &
             REAL(interval)*timestep,                                   &
             model_levels, rows, row_length,                            &
             wet_levels,                                                &
             global_row_length,global_rows,                             &
             n_chem_tracers+n_aero_tracers,                             &
             n_mode_tracers,                                            &
             het_dimn, nhet_std_trop,                                   &
             area,                                                      &
             f3_at_u(1:row_length,1:rows)/two_omega,                    &
             FV_cos_theta_latitude(1:row_length,1:rows),                &
             true_longitude,                                            &
             p_theta_levels(1:row_length,1:rows,1:model_levels),        &
             t_theta_levels(1:row_length,1:rows,:),                     &
             q(1:row_length,1:rows,1:model_levels),                     &
             rel_humid_frac(1:row_length,1:rows,:),                     &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
             all_tracers(1:row_length,1:rows,1:model_levels,            &
                  1:n_chem_tracers+n_aero_tracers),                     &
             all_tracers(1:row_length,1:rows,1:model_levels,            &
                  n_chem_tracers+n_aero_tracers+1:                      &
                  n_chem_tracers+n_aero_tracers+                        &
                  n_mode_tracers),                                      &
             bl_levels,                                                 &
             Tstar(1:row_length,1:rows),                                &
             seaice_frac(1:row_length,1:rows),                          &
             Rough_length(1:row_length,1:rows),                         &
             u_s,                                                       &
             U_scalar_10m,                                              &
             ls_rain3d(1:row_length,1:rows,1:model_levels),             &
             conv_rain3d(1:row_length,1:rows,1:model_levels),           &
             land_fraction,                                             &
             theta_field_size*model_levels,                             &
             delso2_wet_h2o2,                                           &
             delso2_wet_o3,                                             &
             delh2so4_chem,                                             &
             mode_diags,                                                &
             het_rates,                                                 &
             all_emissions(1:row_length,1:rows,                         &
                           1:n_chem_emissions),                         &
             em_index,                                                  &
             SO2_volc_3D(1:row_length,1:rows,1:model_levels),           &
             BC_biom_3D(1:row_length,1:rows,1:model_levels),            &
             OC_biom_3D(1:row_length,1:rows,1:model_levels),            &
             cloud_frac(1:row_length,1:rows,1:model_levels),            &
             cloud_liq_frac(1:row_length,1:rows,1:model_levels),        &
             cloud_liq_water(1:row_length,1:rows,1:model_levels),       &
             offx, offy,                                                &
             r_rho_levels(1:row_length,1:rows,1:model_levels),          &
             r_theta_levels(1:row_length,1:rows,0:model_levels),        &
             z_half, z_half_alllevs, ml_depth, delta_r,                 &
             rhokh_mix,                                                 &
             dtrdz_charney_grid, kent, we_lim(:,:,1:3),                 &
             t_frac(:,:,1:3), zrzi(:,:,1:3), kent_dsc,                  &
             we_lim_dsc(:,:,1:3), t_frac_dsc(:,:,1:3),                  &
             zrzi_dsc(:,:,1:3), zhsc,                                   &
             volume,mass,zbl,                                           &
             uph2so4inaer,                                              &
             wetox_in_aer,                                              &
             chem_diags(:,:,:,icd_cdnc),                                &
             chem_diags(:,:,:,icd_surfarea)                             &             
             )
! Store tropospheric heterogenous rates in chem_diags array
         IF (L_ukca_trophet) THEN 
           IF (n_chem_diags >= 11) THEN
             chem_diags(1:row_length,1:rows,:,10) =                     &
               RESHAPE(het_rates(:,1),(/row_length,rows,model_levels/))
             chem_diags(1:row_length,1:rows,:,11) =                     &
               RESHAPE(het_rates(:,2),(/row_length,rows,model_levels/))
           ELSE
             cmessage='Not enough space for Heterogenous Rates'
             icode = 1
             CALL EREPORT('UKCA_MAIN1',icode,cmessage)
           END IF       ! n_chem_diags >= 11
         END IF        ! L_ukca_trophet

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA AEROSOL MODEL  ',6)

! Call activation scheme if switched on in UMUI/hand edit
        IF(L_ukca_arg_act) THEN

! DEPENDS ON: ukca_activate
          CALL UKCA_ACTIVATE(                                           &
              row_length, rows, model_levels, wet_levels,               &
              bl_levels,                                                &
              theta_field_size,                                         &
              n_mode_tracers,                                           &
              n_mode_diags,  n_chem_diags,                              &
              tr_index,                                                 &
              all_tracers(1:row_length,1:rows,1:model_levels,           &
               n_chem_tracers+n_aero_tracers+1:                         &
               n_chem_tracers+n_aero_tracers+n_mode_tracers),           &
              mode_diags, chem_diags,                                   &
              p_theta_levels(1:row_length,1:rows,1:model_levels),       &
              t_theta_levels(1:row_length,1:rows,1:model_levels),       &
              q(1:row_length,1:rows,1:model_levels),                    &
              qsvp(1:row_length,1:rows,:),                              &
              bl_tke(1:row_length,1:rows,:),                            &
              vertvel(1:row_length,1:rows,0:model_levels-1),            &
              cloud_liq_frac(1:row_length,1:rows,:),                    &
              cloud_liq_water(1:row_length,1:rows,:),                   &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
              STASH_maxlen(38,A_im),STASHwork38)

        END IF   ! L_ukca_arg_act
      END IF   ! L_ukca_mode

! Take a copy of tracer fields if requested for diagnostic purposes
      IF (L_asad_use_chem_diags .AND. L_asad_use_output_tracer)         &
           CALL ASAD_OUTPUT_TRACER(row_length, rows, model_levels,      &
                              n_chem_tracers,                           &
                              all_tracers(:,:,1:model_levels,:), ierr)

      IF ( l_endgame ) THEN
! For ENDGAME Copy level1 to level0 before QPOS/copy to D1
      all_tracers(:,:,0,:) = all_tracers(:,:,1,:)
      q(:,:,0) =  q(:,:,1)
      END IF
      
!     Return fields to D1
      CALL PUTD1FLDS()

! In case water vapour feedback is activated, call Q_POS_CNTL to remove points
! with low water vapour and set polar rows to uniform constants
!! Not required for ENDGAME - done during tracer conservation
      IF ( .NOT. l_endgame )  THEN
       IF (L_ukca_h2o_feedback)                                         &
! DEPENDS ON: q_pos_ctl
        CALL Q_Pos_Ctl(  D1(jq(klev1)), row_length, rows, wet_levels, 1,&
                         bl_levels, global_row_length, global_rows,     &
                         mype, nproc, qdims_s%halo_i, qdims_s%halo_j,   &
                         model_domain,                                  &
                         halo_type_extended, q_pos_method, qlimit,      &
                         l_qpos_diag_pr, qpos_diag_limit,               &
                         'Q call from UKCA_main')

! Remove negatives from all chemical tracers and reset polar rows. Note: This
! call replaces update my_tracer_qpos4.mf77 that would have placed this call
! into atmstep2 subroutine. The call is necessary here to avoid a crash
! if radiative feedback is selected.

! In ENDGAME, Check for -ve values directly in radiation (atmos_phys1)
! DEPENDS ON: q_pos_ctl
        CALL Q_Pos_Ctl(                                                 &
            D1(jtr_ukca(klev1,1)), row_length, rows, tr_ukca*tr_levels, &
                tr_ukca, bl_levels, global_row_length, global_rows,     &
               mype, nproc, tdims_s%halo_i, tdims_s%halo_j,             &
                model_domain,                                           &
                halo_type_single, q_pos_tracer_method,  0.0 ,           &
                l_qpos_diag_pr, qpos_diag_limit,                        &
                'UKCA tracers call from UKCA_main')

      END IF   ! If not ENDGAME 
      
! Add output fields to STASHwork arrays.

! ====================================== 
 
      im_index = internal_model_index(atmos_im) 
 
! MODE diagnostics  (items 38,201 - 38,212 now done in ukca_mode_ems_um)
      icnt=0
      DO l=1,nmax_mode_diags
        IF (UkcaD1codes(imode_first+l-1)%item /= IMDI) THEN
          icnt = icnt + 1
          item = UkcaD1codes(imode_first+l-1)%item
          section = MODE_diag_sect
          IF (sf(item,section) .AND. item > item1_mode_diags+12) THEN
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (stashwork38(si(item,section,im_index)),   &
              mode_diags(:,:,:,icnt),                                   &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,section,im_index)),len_stlist,    &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,section,item,icode,cmessage)

            IF (icode >  0) THEN

              CALL EREPORT('UKCA_MAIN1',section*1000+item,cmessage)
            END IF
          END IF
        END IF
      END DO       ! 1,n_mode_diags

! Stratospheric flux diagnostics 
      icnt = 0 
      DO l=1,nmax_strat_fluxdiags 
        IF (UkcaD1codes(istrat_first+l-1)%item /= IMDI) THEN 
          icnt = icnt + 1 
          item = UkcaD1codes(istrat_first+l-1)%item 
          section = mode_diag_sect 
          IF (sf(item,section)) THEN 
! DEPENDS ON: copydiag_3d 
            Call copydiag_3d (stashwork38(si(item,section,im_index)),   & 
              strat_fluxdiags(:,:,:,icnt),                              & 
              row_length,rows,model_levels,0,0,0,0, at_extremity,       & 
              stlist(1,stindex(1,item,section,im_index)),len_stlist,    & 
              stash_levels,num_stash_levels+1,                          & 
              atmos_im,section,item,icode,cmessage) 

            IF (icode >  0) THEN 
 
              CALL EREPORT('UKCA_MAIN1',section*1000+item,cmessage) 
            END IF 
          ELSE 
            cmessage=' Strat flux item not found by STASH flag array' 
 
            CALL EREPORT('UKCA_MAIN1',section*1000+item,cmessage) 
          END IF 
        END IF 
      END DO       ! 1,nmax_strat_fluxdiags 

! ASAD flux diagnostics: fill STASHwork - needs to be done BEFORE ste calc below
      IF (L_asad_use_chem_diags) THEN
         stashsize = SIZE(STASHwork50)
         CALL asad_flux_put_stash(                                      &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
              stashsize,STASHwork50,                                    &
              row_length,rows,model_levels,do_chemistry)
      END IF


! ASAD Flux Diags Strat-Trop Exchange - done AFTER asad_flux_put_stash
      IF (L_asad_use_chem_diags .AND. L_asad_use_STE) THEN
        CALL asad_tendency_STE(                                         &
                     row_length,rows,model_levels,                      &
                     n_chem_tracers,fluxdiag_all_tracers,               &
                     totnodens,volume,mass,timestep,                    &
                     calculate_STE,.TRUE.,ierr)
      END IF

! DEPENDS ON: stash 
        CALL STASH(a_sm,a_im,34,STASHwork34,                            &
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
          icode,cmessage) 

        IF (icode >  0) THEN 
          cmessage=" Error in STASH"//cmessage 
          errcode=34
          CALL EREPORT('UKCA_MAIN1',errcode,cmessage) 
        END IF 

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,38,STASHwork38,                              &
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
          icode,cmessage)
     
      IF (icode >  0) THEN
        cmessage=" Error in STASH"//cmessage
        errcode=38
        CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
      END IF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,50,STASHwork50,                              &
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
          icode,cmessage)
     
      IF (icode >  0) THEN
        cmessage=" Error in STASH"//cmessage
        errcode=50
        CALL EREPORT('UKCA_MAIN1',errcode,cmessage)
      END IF

!     Deallocate arrays

! These are currently only allocated on first call - see ukca_iniasad
!      CALL ASAD_MOD_FINAL


! POSSIBLE THINGS TO DEALLOCATE
      IF (ALLOCATED(theta_latitude)) DEALLOCATE(theta_latitude)
      IF (ALLOCATED(v_latitude)) DEALLOCATE(v_latitude)
      IF (ALLOCATED(sinv_latitude)) DEALLOCATE(sinv_latitude)
      IF (ALLOCATED(net_surf_sw)) DEALLOCATE(net_surf_sw)
      IF (ALLOCATED(tot_surf_sw)) DEALLOCATE(tot_surf_sw)
      IF (ALLOCATED(trmol_post_atmstep)) DEALLOCATE(trmol_post_atmstep)  ! moles
      IF (ALLOCATED(um_ozone1D)) DEALLOCATE(um_ozone1D)
      IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)
      IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)
      IF (ALLOCATED(T_chem)) DEALLOCATE(T_chem)
      IF (ALLOCATED(Q_chem)) DEALLOCATE(Q_chem)
! from PUTD1FLDS
      IF (ALLOCATED(tracer1)) DEALLOCATE(tracer1)
      IF (ALLOCATED(emission1)) DEALLOCATE(emission1)
      IF (ALLOCATED(chem_diag1)) DEALLOCATE(chem_diag1)
      IF (ALLOCATED(all_tracers)) DEALLOCATE(all_tracers)
      IF (ALLOCATED(fluxdiag_all_tracers)) DEALLOCATE(fluxdiag_all_tracers)
      IF (ALLOCATED(all_emissions)) DEALLOCATE(all_emissions)
      IF (ALLOCATED(chem_diags)) DEALLOCATE(chem_diags)


! ACTUAL THINGS TO DEALLOCATE
      IF (ALLOCATED(totnodens)) DEALLOCATE(totnodens)
      IF (ALLOCATED(um_ozone3d)) DEALLOCATE(um_ozone3d)
      IF (ALLOCATED(conv_cloud_base)) DEALLOCATE(conv_cloud_base)
      IF (ALLOCATED(conv_cloud_top)) DEALLOCATE(conv_cloud_top)
      IF (ALLOCATED(land_sea_mask)) DEALLOCATE(land_sea_mask)
      IF (ALLOCATED(aircraftems)) DEALLOCATE(aircraftems)
      IF (ALLOCATED(theta)) DEALLOCATE(theta)
      IF (ALLOCATED(p_rho_levels)) DEALLOCATE(p_rho_levels)
      IF (ALLOCATED(p_theta_levels)) DEALLOCATE(p_theta_levels)
      IF (ALLOCATED(exner_theta_levels)) DEALLOCATE(exner_theta_levels)
      IF (ALLOCATED(exner_rho_levels)) DEALLOCATE(exner_rho_levels)
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(rho_r2)) DEALLOCATE(rho_r2)
      IF (ALLOCATED(qcl)) DEALLOCATE(qcl)
      IF (ALLOCATED(qcf)) DEALLOCATE(qcf)
      IF (ALLOCATED(ls_ppn_frac)) DEALLOCATE(ls_ppn_frac)
      IF (ALLOCATED(cloud_liq_water)) DEALLOCATE(cloud_liq_water)
      IF (ALLOCATED(cloud_ice_content)) DEALLOCATE(cloud_ice_content)
      IF (ALLOCATED(ls_rain3d)) DEALLOCATE(ls_rain3d)
      IF (ALLOCATED(ls_snow3d)) DEALLOCATE(ls_snow3d)
      IF (ALLOCATED(ls_ppn3d)) DEALLOCATE(ls_ppn3d)
      IF (ALLOCATED(conv_rain3d)) DEALLOCATE(conv_rain3d)
      IF (ALLOCATED(conv_snow3d)) DEALLOCATE(conv_snow3d)
      IF (ALLOCATED(conv_ppn3d)) DEALLOCATE(conv_ppn3d)
      IF (ALLOCATED(conv_cloud_amount)) DEALLOCATE(conv_cloud_amount)
      IF (ALLOCATED(area_cloud_fraction)) DEALLOCATE(area_cloud_fraction)
      IF (ALLOCATED(so4_aitken)) DEALLOCATE(so4_aitken)
      IF (ALLOCATED(so4_accum)) DEALLOCATE(so4_accum)
      IF (ALLOCATED(sulphate_od)) DEALLOCATE(sulphate_od)
      IF (ALLOCATED(dj)) DEALLOCATE(dj)
      IF (ALLOCATED(rel_humid_frac)) DEALLOCATE(rel_humid_frac)
      IF (ALLOCATED(qsvp)) DEALLOCATE(qsvp)
      IF (ALLOCATED(U_scalar_10m)) DEALLOCATE(U_scalar_10m)
      IF (ALLOCATED(pstar)) DEALLOCATE(pstar)
      IF (ALLOCATED(tstar)) DEALLOCATE(tstar)
      IF (ALLOCATED(dust_flux)) DEALLOCATE(dust_flux)
      IF (ALLOCATED(Soil_Layer_Moisture)) DEALLOCATE(Soil_Layer_Moisture)
      IF (ALLOCATED(Tile_Frac)) DEALLOCATE(Tile_Frac)
      IF (ALLOCATED(Tstar_tile)) DEALLOCATE(Tstar_tile)
      IF (ALLOCATED(Snow_tile)) DEALLOCATE(Snow_tile)
      IF (ALLOCATED(Frac_types)) DEALLOCATE(Frac_types)
      IF (ALLOCATED(Rough_length)) DEALLOCATE(Rough_length)
      IF (ALLOCATED(fland)) DEALLOCATE(fland)
      IF (ALLOCATED(zbl)) DEALLOCATE(zbl)
      IF (ALLOCATED(surf_hf)) DEALLOCATE(surf_hf)
      IF (ALLOCATED(seaice_frac)) DEALLOCATE(seaice_frac)
      IF (ALLOCATED(conv_cloud_lwp)) DEALLOCATE(conv_cloud_lwp)
      IF (ALLOCATED(u_s)) DEALLOCATE(u_s)
      IF (ALLOCATED(stcon)) DEALLOCATE(stcon)
      IF (ALLOCATED(laift_lp)) DEALLOCATE(laift_lp)
      IF (ALLOCATED(canhtft_lp)) DEALLOCATE(canhtft_lp)
      IF (ALLOCATED(z0tile_lp)) DEALLOCATE(z0tile_lp)
      IF (ALLOCATED(canwctile_lp)) DEALLOCATE(canwctile_lp)
      IF (ALLOCATED(Thick_bl_levels)) DEALLOCATE(Thick_bl_levels)
      IF (ALLOCATED(STASHwork34)) DEALLOCATE(STASHwork34)
      IF (ALLOCATED(STASHwork38)) DEALLOCATE(STASHwork38)
      IF (ALLOCATED(STASHwork50)) DEALLOCATE(STASHwork50)
      IF (ALLOCATED(p_layer_boundaries)) DEALLOCATE(p_layer_boundaries)
      IF (ALLOCATED(t_theta_levels)) DEALLOCATE(t_theta_levels)
      IF (ALLOCATED(ch4_wetl_emiss)) DEALLOCATE(ch4_wetl_emiss)
      IF (ALLOCATED(delSO2_wet_h2o2)) DEALLOCATE(delSO2_wet_h2o2)
      IF (ALLOCATED(delSO2_wet_o3)) DEALLOCATE(delSO2_wet_o3)
      IF (ALLOCATED(delh2so4_chem)) DEALLOCATE(delh2so4_chem)
      IF (ALLOCATED(delSO2_drydep)) DEALLOCATE(delSO2_drydep)
      IF (ALLOCATED(delSO2_wetdep)) DEALLOCATE(delSO2_wetdep)
      IF (ALLOCATED(SO2_volc_3D)) DEALLOCATE(SO2_volc_3D)
      IF (ALLOCATED(SO2_biom_3D)) DEALLOCATE(SO2_biom_3D)
      IF (ALLOCATED(BC_biom_3D)) DEALLOCATE(BC_biom_3D)
      IF (ALLOCATED(OC_biom_3D)) DEALLOCATE(OC_biom_3D)
      IF (ALLOCATED(so4_sa)) DEALLOCATE(so4_sa)
      IF (ALLOCATED(tropopause_height)) DEALLOCATE(tropopause_height)
      IF (ALLOCATED(strat_fluxdiags)) DEALLOCATE(strat_fluxdiags)
      IF (ALLOCATED(mode_diags)) DEALLOCATE(mode_diags)
      IF (ALLOCATED(cloud_frac)) DEALLOCATE(cloud_frac)
      IF (ALLOCATED(dtrdz_charney_grid)) DEALLOCATE(dtrdz_charney_grid)
      IF (ALLOCATED(kent)) DEALLOCATE(kent)
      IF (ALLOCATED(kent_dsc)) DEALLOCATE(kent_dsc)
      IF (ALLOCATED(land_albedo)) DEALLOCATE(land_albedo)
      IF (ALLOCATED(ml_depth)) DEALLOCATE(ml_depth)
      IF (ALLOCATED(PV_on_theta_mlevs)) DEALLOCATE(PV_on_theta_mlevs)
      IF (ALLOCATED(rhokh_mix)) DEALLOCATE(rhokh_mix)
      IF (ALLOCATED(bl_tke)) DEALLOCATE(bl_tke)
      IF (ALLOCATED(vertvel)) DEALLOCATE(vertvel)
      IF (ALLOCATED(t_frac)) DEALLOCATE(t_frac)
      IF (ALLOCATED(t_frac_dsc)) DEALLOCATE(t_frac_dsc)
      IF (ALLOCATED(um_ozone)) DEALLOCATE(um_ozone)
      IF (ALLOCATED(we_lim)) DEALLOCATE(we_lim)
      IF (ALLOCATED(we_lim_dsc)) DEALLOCATE(we_lim_dsc)
      IF (ALLOCATED(zhsc)) DEALLOCATE(zhsc)
      IF (ALLOCATED(zrzi)) DEALLOCATE(zrzi)
      IF (ALLOCATED(zrzi_dsc)) DEALLOCATE(zrzi_dsc)
      IF (ALLOCATED(cloud_liq_frac)) DEALLOCATE(cloud_liq_frac)

      first=.FALSE.
      IF (lhook) CALL dr_hook('UKCA_MAIN1',zhook_out,zhook_handle)
      RETURN

      CONTAINS
! ######################################################################
      SUBROUTINE GETD1FLDS(N)

        IMPLICIT NONE
! Interface block to allow generic calls
      INTERFACE UKCA_EXTRACT_D1_DATA

! DEPENDS ON: UKCA_EXTRACT_D1_DATA1D
        SUBROUTINE UKCA_EXTRACT_D1_DATA1D(                             &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first,N,X)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
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
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA1D

! DEPENDS ON: UKCA_EXTRACT_D1_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_DATA2D(                             &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first,N,X)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
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
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_INTEGER_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_INTEGER_DATA2D(                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first,N,X)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
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
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          INTEGER, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          INTEGER, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_INTEGER_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_LOGICAL_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_LOGICAL_DATA2D(                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first,N,X)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
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
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          LOGICAL, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          LOGICAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_LOGICAL_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_DATA3D
        SUBROUTINE UKCA_EXTRACT_D1_DATA3D(                             &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first,N,X)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
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
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:,:,:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA3D
      END INTERFACE

      INTERFACE
! DEPENDS ON: UKCA_SET_ARRAY_BOUNDS
        SUBROUTINE UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          USE UKCA_D1_DEFS
          USE Field_Types
          USE UM_ParParams
          USE Control_Max_Sizes
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: N            ! id of array
          INTEGER, INTENT(OUT) :: I1,I2,J1,J2  ! array bounds
        END SUBROUTINE UKCA_SET_ARRAY_BOUNDS
      END INTERFACE


      INTEGER, INTENT(IN)     :: N             ! id of array

      INTEGER, SAVE           :: next_tr       ! count tracers
      INTEGER, SAVE           :: idust         ! count dust divs
      INTEGER, SAVE           :: next_em       ! count emissions
      INTEGER, SAVE           :: next_cd       ! count chem diags
      INTEGER                 :: I1,I2,J1,J2   ! array bounds
      integer :: iq,jq,kq

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)         :: zhook_handle

!     Items in range 1-150 fill all_tracers array in sequence

      IF (lhook) CALL dr_hook('GETD1FLDS',zhook_in,zhook_handle)
      IF (UkcaD1codes(N)%section == UKCA_sect .AND.                    &
          UkcaD1codes(N)%item <= n_all_tracers) THEN
        IF (.NOT. ALLOCATED(tracer1)) THEN
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tracer1(I1:I2,J1:J2,                                &
                         klev1:ukcaD1codes(N)%len_dim3))
          ALLOCATE(all_tracers(I1:I2,J1:J2,                            &
                         klev1:ukcaD1codes(N)%len_dim3,n_use_tracers))
          next_tr = 1
        END IF

        IF (.NOT. ALLOCATED(tr_index)) ALLOCATE(tr_index(n_use_tracers))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,tracer1)
        all_tracers(:,:,:,next_tr) = tracer1(:,:,:)
        tr_index(next_tr)          = UkcaD1codes(N)%item
        next_tr                    = next_tr+1

!     Items in range 150+23 - 150+49 fill all_emissions array in sequence

      ELSE IF ((UkcaD1codes(N)%section == UKCA_ems_sect .AND.          &
          UkcaD1codes(N)%item >= n_emiss_first .AND.                   &
          UkcaD1codes(N)%item <= n_emiss_last ) .OR.                   &
         (UkcaD1codes(N)%section == 0 .AND.                            &
         (UkcaD1codes(N)%item == 58  .OR.                              &
          UkcaD1codes(N)%item == 126 .OR.                              &
          UkcaD1codes(N)%item == 127)) .OR.                            &
         (UkcaD1codes(N)%section == 17 .AND.                           &
          UkcaD1codes(N)%item == 205 )) THEN
        IF (.NOT. ALLOCATED(emission1)) THEN
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(emission1(I1:I2,J1:J2))
          ALLOCATE(all_emissions(I1:I2,J1:J2,n_chem_emissions))
          next_em = 1
        END IF
        IF (.NOT. ALLOCATED(em_index))                                 &
                ALLOCATE(em_index(n_chem_emissions))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,emission1)
        all_emissions(:,:,next_em) = emission1(:,:)
        em_index(next_em)          = UkcaD1codes(N)%item
        next_em                    = next_em+1

! Add 3-D SO2 emissions
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                        &
              UkcaD1codes(N)%item == 121) THEN
        CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        ALLOCATE(SO2_volc_3D(I1:I2,J1:J2,klev1:ukcaD1codes(N)%len_dim3))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,SO2_volc_3D)

! Add 3-D BC biomass emissions
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                        &
              UkcaD1codes(N)%item == 322) THEN
        CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        ALLOCATE(BC_biom_3D(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,BC_biom_3D)

! Add 3-D OC biomass emissions
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                       &
              UkcaD1codes(N)%item == 323) THEN
        CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        ALLOCATE(OC_biom_3D(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,OC_biom_3D)

! Add 3-D SO2 biomass emissions
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                       &
              UkcaD1codes(N)%item == 324) THEN
        CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        ALLOCATE(SO2_biom_3D(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,SO2_biom_3D)

! Add 3-D Aircraft emissions
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                       &
              UkcaD1codes(N)%item == 340) THEN
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(aircraftems(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,aircraftems)

! Add  Dust emissions
      ELSE IF (UkcaD1codes(N)%section == 3 .AND.                       &
              UkcaD1codes(N)%item >= 401 .AND.                         &
              UkcaD1codes(N)%item <= 401+ndiv-1 ) THEN
        IF (.NOT. ALLOCATED(dust_flux)) THEN
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_flux(I1:I2,J1:J2,ndiv))
          idust = 1
        END IF
        CALL UKCA_EXTRACT_D1_DATA(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,dust_flux(:,:,idust))
        idust = idust + 1

!     Prognostics: fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                       &
              UkcaD1codes(N)%prognostic) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(4)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(theta(I1:I2,J1:J2,                                  &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,theta)
        CASE(9)
          ALLOCATE(soil_layer_moisture(ukcaD1codes(N)%len_dim1,        &
                                       ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,soil_layer_moisture)
        CASE(10)
           ! if Q has been allocated then this does not need to be taken again 
           ! and will in fact cause an error if it is.
          IF (.NOT. ALLOCATED(q)) THEN
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
             ALLOCATE(q(I1:I2,J1:J2,                                   &
                  klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,q)

! Fill advected chemical H2O tracer if required
            IF (l_UKCA_advh2o) THEN
                all_tracers(:,:,1:wet_levels,next_tr) =                &
                     q(qdims_s%i_start:qdims_s%i_end,                  &
                       qdims_s%j_start:qdims_s%j_end,1:wet_levels)
              next_tr = next_tr + 1
            END IF
          END IF
        CASE(12)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(qcf(I1:I2,J1:J2,                                    &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,qcf)
        CASE(16)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_lwp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_cloud_lwp)
        CASE(24)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tstar(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,tstar)
        CASE(25)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zbl(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,zbl)
        CASE(26)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(Rough_length(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,Rough_length)
        CASE(30)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(land_sea_mask(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,land_sea_mask)
        CASE(31)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(seaice_frac(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,seaice_frac)
        CASE(60)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(um_ozone(I1:I2,J1:J2,                               &
                             klev1:ukcaD1codes(N)%len_dim3))
          um_ozone = RESHAPE(D1(UkcaD1Codes(N)%address:                &
             UkcaD1Codes(N)%address+UkcaD1Codes(N)%length-1),          &
             (/SIZE(um_ozone,DIM=1),SIZE(um_ozone,DIM=2),              &
             SIZE(um_ozone,DIM=3)/))
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        CASE(103)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(so4_aitken(I1:I2,J1:J2,                             &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,so4_aitken)
        CASE(104)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(so4_accum(I1:I2,J1:J2,                              &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,so4_accum)
       CASE(150)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(vertvel(I1:I2,J1:J2,0:ukcaD1codes(N)%len_dim3-1))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,vertvel)
        CASE(211)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_amount(I1:I2,J1:J2,                      &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_cloud_amount)
        CASE(216)
! DEPENDS ON: ukca_set_array_bounds
          ALLOCATE(frac_types(ukcaD1codes(N)%len_dim1,                 &
                              ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,frac_types)
        CASE(217)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(laift_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,laift_lp)
        CASE(218)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(canhtft_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,canhtft_lp)
        CASE(229)
! DEPENDS ON: ukca_set_array_bounds
           CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
           ALLOCATE(canwctile_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,canwctile_lp)
        CASE(233)
          ALLOCATE(tstar_tile(ukcaD1codes(N)%len_dim1,                 &
                              ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,tstar_tile)
        CASE(234)
! DEPENDS ON: ukca_set_array_bounds
           CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
           ALLOCATE(z0tile_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,z0tile_lp)
        CASE(240)
          ALLOCATE(snow_tile(ukcaD1codes(N)%len_dim1,                  &
                             ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,snow_tile)
        CASE(253)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(rho_r2(I1:I2,J1:J2,                                 &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,rho_r2)
        CASE(254)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(qcl(I1:I2,J1:J2,                                    &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,qcl)
        CASE(255)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(exner_rho_levels(I1:I2,J1:J2,                       &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,exner_rho_levels)
        CASE(265)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(area_cloud_fraction(I1:I2,J1:J2,                    &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,area_cloud_fraction)
        CASE(266)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_frac(I1:I2,J1:J2,                             &
                              klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,cloud_frac)
        CASE(267)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_liq_frac(I1:I2,J1:J2,                         &
                                  klev1:ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,cloud_liq_frac)
        CASE(505)
          ALLOCATE(fland(ukcaD1codes(N)%len_dim1))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,fland)
        CASE(510)
          IF (mod((timestep_number-1),a_sw_radstep_prog) == 0) THEN   ! Only on radiation TS 
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(land_albedo(I1:I2,J1:J2))
            CALL UKCA_EXTRACT_D1_DATA(                                 &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
            first,N,land_albedo)
          END IF

        CASE DEFAULT
          cmessage='Item not found in prognostic case statement'

          CALL EREPORT('GETD1FLDS',ukcaD1codes(N)%item,cmessage)
        END SELECT

! Diagnostics (section 0): fill appropriate array
      ELSE IF (UkcaD1codes(N)%section == 0 .AND.                       &
              .NOT. UkcaD1codes(N)%prognostic) THEN
        SELECT CASE(UkcaD1codes(N)%item)

        CASE(406)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(exner_theta_levels(I1:I2,J1:J2,                     &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,exner_theta_levels)
        CASE(407)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(p_rho_levels(I1:I2,J1:J2,                           &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,p_rho_levels)
        CASE(408)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(p_theta_levels(I1:I2,J1:J2,                         &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,p_theta_levels)
        CASE(409)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(pstar(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,pstar)
        CASE DEFAULT
          cmessage='N not found in diagnostic(0) case statement'

          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (section 1): fill appropriate array
      ELSE IF (UkcaD1codes(N)%section == 1) THEN
        SELECT CASE(UkcaD1codes(N)%item)

        CASE(201)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(net_surf_SW(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,net_surf_SW)
        CASE(235)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tot_surf_SW(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,tot_surf_SW)
        CASE DEFAULT
          cmessage='N not found in diagnostic(1) case statement'

          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (section 2): fill appropriate array
      ELSE IF (UkcaD1codes(N)%section == 2) THEN
        SELECT CASE(UkcaD1codes(N)%item)

        CASE(284)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(sulphate_od(I1:I2,J1:J2,                            &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,sulphate_od)
        CASE DEFAULT
          cmessage='N not found in diagnostic(2) case statement'
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

!     Diagnostics (section 3): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 3) THEN
        SELECT CASE(UkcaD1codes(N)%item)
       CASE(25)
! DEPENDS ON: ukca_set_array_bounds
         CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
         ALLOCATE(ml_depth(I1:I2,J1:J2))
         CALL UKCA_EXTRACT_D1_DATA(                                    &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,ml_depth)
        CASE(60)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(rhokh_mix(I1:I2,J1:J2,                              &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,rhokh_mix)
        CASE(64)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dtrdz_charney_grid(I1:I2,J1:J2,                     &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,dtrdz_charney_grid)
        CASE(65)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(kent(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,kent)
        CASE(66)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(we_lim(I1:I2,J1:J2,                                 &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,we_lim)
        CASE(67)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(t_frac(I1:I2,J1:J2,                                 &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,t_frac)
        CASE(68)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zrzi(I1:I2,J1:J2,                                   &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,zrzi)
        CASE(69)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(kent_dsc(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,kent_dsc)
        CASE(70)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(we_lim_dsc(I1:I2,J1:J2,                             &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,we_lim_dsc)
        CASE(71)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(t_frac_dsc(I1:I2,J1:J2,                             &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,t_frac_dsc)
        CASE(72)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zrzi_dsc(I1:I2,J1:J2,                               &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,zrzi_dsc)
        CASE(73)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zhsc(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,zhsc)
        CASE(230)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(u_scalar_10m(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,u_scalar_10m)
        CASE(217)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(surf_hf(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,surf_hf)
        CASE(462)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(stcon(I1:I2,J1:J2,Ukcad1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,stcon)
        CASE(465)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(u_s(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,u_s)
      CASE(473)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(bl_tke(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,bl_tke)
        CASE DEFAULT
          cmessage='N not found in diagnostic(3) case statement'

          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (section 4): fill appropriate array
      ELSE IF (UkcaD1codes(N)%section == 4) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(205)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_liq_water(I1:I2,J1:J2,                        &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,cloud_liq_water)
        CASE(206)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_ice_content(I1:I2,J1:J2,                      &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,cloud_ice_content)
        CASE(222)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_rain3d(I1:I2,J1:J2,                              &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,ls_rain3d)
        CASE(223)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_snow3d(I1:I2,J1:J2,                              &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,ls_snow3d)
        CASE(227)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_ppn_frac(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,ls_ppn_frac)
        CASE DEFAULT
          cmessage='N not found in diagnostic(4) case statement'
          errcode=1
          WRITE(6,'(A72,I10)') cmessage,n
          CALL EREPORT('GETD1FLDS',errcode,cmessage)
        END SELECT

! Diagnostics (section 5): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 5) THEN
        SELECT CASE(UkcaD1codes(N)%item)
       CASE(227)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_rain3d(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_rain3d)
        CASE(228)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_snow3d(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_snow3d)
        CASE(218)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_base(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_cloud_base)
        CASE(219)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_top(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,conv_cloud_top)
        CASE DEFAULT
          cmessage='N not found in diagnostic(5) case statement'
          errcode=1

          CALL EREPORT('GETD1FLDS',errcode,cmessage)
        END SELECT

! Diagnostics (section 8): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 8) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(242)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ch4_wetl_emiss(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,ch4_wetl_emiss)
        CASE DEFAULT
          cmessage='N not found in diagnostic(8) case statement'
          errcode=1

          CALL EREPORT('GETD1FLDS',errcode,cmessage)
        END SELECT

! Diagnostics (section 15): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 15) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(218)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(PV_on_theta_mlevs(I1:I2,J1:J2,                      &
                   ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,PV_on_theta_mlevs)
        CASE DEFAULT
          cmessage='N not found in diagnostic(8) case statement'
          errcode=1

          CALL EREPORT('GETD1FLDS',errcode,cmessage)
       END SELECT

! Diagnostics (Section 30): fill appropriate array
      ELSE IF (UkcaD1codes(N)%section == 30) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(453)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tropopause_height(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                    &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,tropopause_height)
        CASE DEFAULT
          cmessage='N not found in diagnostic(30) case statement'

          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (Section 34): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == UKCA_sect) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(item1_chem_diags:item2_chem_diags)                ! chemical diagnostics
          IF (.NOT. ALLOCATED(chem_diag1)) THEN
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(chem_diag1(I1:I2,J1:J2,                           &
                     ukcaD1codes(N)%len_dim3))
            ALLOCATE(chem_diags(I1:I2,J1:J2,                           &
                     ukcaD1codes(N)%len_dim3,n_chem_diags))
            next_cd = 1
          END IF
          CALL UKCA_EXTRACT_D1_DATA(                                   &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
           first,N,chem_diag1)
           chem_diags(:,:,:,next_cd) = chem_diag1(:,:,:)
           next_cd                   = next_cd + 1
        CASE DEFAULT
          cmessage='N not found in ukca_section case statement'

          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

      ELSE
        cmessage='N not located in IF statement'

        CALL EREPORT('GETD1FLDS',N,cmessage)
      END IF
      IF (lhook) CALL dr_hook('GETD1FLDS',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE GETD1FLDS
! ######################################################################
       SUBROUTINE PUTD1FLDS

       USE chsunits_mod, ONLY : nunits

       IMPLICIT NONE
! Put tracer and non-tracer species fields back into D1 after UKCA

       REAL, DIMENSION(:), ALLOCATABLE :: data1
       REAL, DIMENSION(:), ALLOCATABLE :: data2
       REAL, DIMENSION(:), ALLOCATABLE :: data3

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle

       IF (lhook) CALL dr_hook('PUTD1FLDS',zhook_in,zhook_handle)

! DEPENDS ON: ukca_set_array_bounds
       CALL UKCA_SET_ARRAY_BOUNDS(1,I1,I2,J1,J2)
       ALLOCATE(data1(ukcaD1codes(1)%length))

       IF (L_ukca_chem) THEN
! DEPENDS ON: ukca_set_array_bounds
       CALL UKCA_SET_ARRAY_BOUNDS(idiag_first,I1,I2,J1,J2)
       IF (.NOT. ALLOCATED(data2)) &
            ALLOCATE(data2(ukcaD1codes(idiag_first)%length))
       END IF

       DO N=1,n_use_tracers
         IF ((ukcaD1codes(N)%item /= 10) .OR.                        &
             (ukcaD1codes(N)%section /= 0)) THEN
            
            ! Check that size of data1 is correct
            IF (ukcaD1codes(N)%length /= SIZE(data1)) THEN
               IF (ALLOCATED(data1)) DEALLOCATE(data1)
               ALLOCATE(data1(ukcaD1codes(N)%length))
            END IF

           tracer1(:,:,:)=all_tracers(:,:,:,N)
           data1=RESHAPE(tracer1,(/SIZE(data1)/))

           D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+           &
              ukcaD1codes(N)%length-1)=data1
         END IF

       END DO   ! N tracer loop

       IF (L_ukca_chem) THEN
       DO N=idiag_first,idiag_last

          ! Check that size of data2 is correct
          IF (ukcaD1codes(N)%length /= SIZE(data2)) THEN
             IF (ALLOCATED(data2)) DEALLOCATE(data2)
             ALLOCATE(data2(ukcaD1codes(N)%length))
          END IF

         chem_diag1(:,:,:)=chem_diags(:,:,:,N-idiag_first+1)
         data2=RESHAPE(chem_diag1,(/SIZE(data2)/))

           D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+           &
              ukcaD1codes(N)%length-1)=data2

         END DO   ! N chemical diagnostics loop
       END IF

       IF (L_ukca_h2o_feedback) THEN
! Copy water vapour back into D1 array.
         N = size(q)
         IF (ALLOCATED(data3)) DEALLOCATE(data3)
         ALLOCATE(data3(N))
         data3=RESHAPE(q,(/N/))

! Update D1
         DO N=n_use_tracers+n_use_emissions+1,                         &
              n_use_tracers+n_use_emissions+n_in_progs
           IF ((ukcaD1codes(N)%item  == 10) .AND.                      &
               (ukcaD1codes(N)%section == 0 ))                         &
             D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+         &
               ukcaD1codes(N)%length-1)=data3
         END DO

         IF (ALLOCATED(data3)) DEALLOCATE(data3)
       END IF

       IF (ALLOCATED(data1)) DEALLOCATE(data1)
       IF (ALLOCATED(data2)) DEALLOCATE(data2)
       IF (ALLOCATED(tracer1)) DEALLOCATE(tracer1)
       IF (ALLOCATED(emission1)) DEALLOCATE(emission1)
       IF (ALLOCATED(chem_diag1)) DEALLOCATE(chem_diag1)
       IF (ALLOCATED(all_tracers)) DEALLOCATE(all_tracers)
       IF (ALLOCATED(all_emissions)) DEALLOCATE(all_emissions)
       IF (ALLOCATED(chem_diags)) DEALLOCATE(chem_diags)

       IF (lhook) CALL dr_hook('PUTD1FLDS',zhook_out,zhook_handle)
       RETURN

       END SUBROUTINE PUTD1FLDS
! ######################################################################

       END SUBROUTINE UKCA_MAIN1
