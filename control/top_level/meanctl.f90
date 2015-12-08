! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!      Subroutine: MEANCTL-----------------------------------------
!
!      Purpose: To accumulate partial sums and create time-meaned data
!
!      Programming standard: UM Doc Paper 3
!
!      External documentation: UMDP C5 - Control of means calculations
!
!      ------------------------------------------------------------
!      Interface and arguments:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

SUBROUTINE meanctl (                                              &
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
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
  ind_im,meanlev,icode,cmessage)
!
      USE filenamelength_mod,       ONLY : &
          filenamelength
      USE yomhook,                  ONLY : lhook, dr_hook
      USE parkind1,                 ONLY : jprb, jpim
      USE io_configuration_mod,     ONLY : io_field_padding
      USE IO
      USE model_file,               ONLY : model_file_close, model_file_open
      USE MPPIO_job_control,        ONLY : jobCntrl     
      USE MPPIO_job_control_common, ONLY : JC_delete
      USE atm_fields_bounds_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE Decomp_DB
      USE um_input_control_mod,     ONLY : lcal360
      USE lookup_addresses
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,    &
          dumpfreqim, meanfreqim, mean_reftimeim, pp_len2_meanim
      USE Submodel_Mod

      USE nlstcall_mod, ONLY : lmean, &
                               lclimrealyr, &
                               expt_id, &
                               job_id, &
                               model_status, &
                               time_convention, &
                               lclimrealyr

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

!   Arguments

!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
! Interface variables initialised through INTFCNSTA
! Namelist read in the interface control routine INTF_CTL.

      INTEGER                                                           &
        intf_row_length                                                 &
                         ! Interface field row length
       ,intf_p_rows                                                     &
                         ! Interface field no of rows
       ,intf_p_levels                                                   &
                         ! Interface field no of levels
       ,intf_q_levels                                                   &
                         ! Interface field no of wet levels
       ,intf_tr_levels                                                  &
                         ! Interface field no of tracer levels
       ,intfwidtha                                                      &
                         ! Width of interface zone (atmosphere)
       ,intf_exthalo_ns                                                 &
                         ! Extended Halo in NS direction
       ,intf_exthalo_ew                                                 &
                         ! Extended Halo in EW direction
       ,a_intf_start_hr                                                 &
                         ! ) Start and End time in
       ,a_intf_freq_hr                                                  &
                         ! ) hours, Frequency in h,m,s for which
       ,a_intf_freq_mn                                                  &
                         ! ) atmosphere interface data
       ,a_intf_freq_sc                                                  &
                         ! ) is to be generated.
       ,a_intf_end_hr                                                   &
                         ! )
       ,intf_pack                                                       &
                         ! Packing Indicator for boundary data
       ,lbc_stream_a                                                    &
                         ! Output streams in UMUI
       ,lbc_unit_no_a                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
       ,lbc_first_r_rho                                                 &
                         ! First rho level at which height is constant
       ,intf_v_int_order(max_n_intf_a)

      REAL                                                              &
        intf_ewspace                                                    &
                         ! E-W grid spacing (degrees)
       ,intf_nsspace                                                    &
                         ! N-S grid spacing (degrees)
       ,intf_firstlat                                                   &
                         ! Latitude of first row (degrees)
       ,intf_firstlong                                                  &
                         ! Longitude of first row (degrees)
       ,intf_polelat                                                    &
                         ! Real latitude of coordinate pole (degrees)
       ,intf_polelong                                                   &
                         ! Real longitude of coordinate pole (degrees)
       ,lbc_z_top_model                                                 &
                         ! Height of top of model
       ,lbc_q_min                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , lambda_intf_p(max_intf_lbcrow_length, max_n_intf_a)             &
      , lambda_intf_u(max_intf_lbcrow_length, max_n_intf_a)             &    
      , phi_intf_p(max_intf_lbcrows, max_n_intf_a)                      &
      , phi_intf_v(max_intf_lbcrows, max_n_intf_a)

      LOGICAL                                                           &
        intf_vert_interp                                                &
                         ! Switch to request vertical interpolation
       ,lnewbnd          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  intf_l_var_lbc(max_n_intf_a)

! Switch to not rotate if input and output grids have same poles.
      LOGICAL intf_avoid_rot(MAX_N_INTF_A)

! Switch to output LBC for Endgame
      LOGICAL intf_l_eg_out(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      CHARACTER(LEN=256) :: intf_vertlevs

! Files for HorzGrid namelist  
      CHARACTER(LEN=256) :: intf_HorzGrid(max_n_intf_a)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
        intf_ewspace(max_n_intf_a)    ,intf_nsspace(max_n_intf_a),      &
        intf_firstlat(max_n_intf_a)   ,intf_firstlong(max_n_intf_a),    &
        intf_polelat(max_n_intf_a)    ,intf_polelong(max_n_intf_a),     &
        intf_row_length(max_n_intf_a) ,intf_p_rows(max_n_intf_a),       &
        intf_p_levels(max_n_intf_a)   ,intf_q_levels(max_n_intf_a),     &
        intf_tr_levels(max_n_intf_a)  ,intfwidtha(max_n_intf_a),        &
        intf_exthalo_ns(max_n_intf_a) ,intf_exthalo_ew(max_n_intf_a),   &
        a_intf_start_hr(max_n_intf_a) ,a_intf_freq_hr(max_n_intf_a),    &
        a_intf_freq_mn(max_n_intf_a)  ,a_intf_freq_sc(max_n_intf_a),    &
        a_intf_end_hr(max_n_intf_a)   ,                                 & 
        lnewbnd(max_n_intf_a)         ,intf_vert_interp(max_n_intf_a),  &
        intf_pack(max_n_intf_a)       ,lbc_stream_a(max_n_intf_a),      &
        lbc_unit_no_a(max_n_intf_a)   ,lbc_first_r_rho(max_n_intf_a),   &
        lbc_z_top_model(max_n_intf_a) ,                                 &
        intf_vertlevs(max_n_intf_a)   ,lbc_q_min,                       &
        intf_l_var_lbc                ,intf_horzgrid,                   &
        lambda_intf_p                 ,lambda_intf_u,                   &
        phi_intf_p                    ,phi_intf_v,                      &
        intf_avoid_rot                ,intf_v_int_order,                &
        intf_l_eg_out
!---------------------------------------------------------------------
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
! TYPINFA start
!
! This file needs file TYPSIZE to be included

      ! Headers for atmosphere interface data sets
      INTEGER::FIXHD_INTFA(LEN_FIXHD,N_INTF_A)        ! Fixed header
      INTEGER::INTHD_INTFA(PP_LEN_INTHD,N_INTF_A)     ! Integer header
      INTEGER::LOOKUP_INTFA(LEN1_LOOKUP,INTF_LOOKUPSA,N_INTF_A)! Lookups

      REAL :: REALHD_INTFA(PP_LEN_REALHD,N_INTF_A)   ! Real header
      REAL :: LEVDEPC_INTFA(MAX_INTF_MODEL_LEVELS,INTF_LEN2_LEVDEPC,    &
     &  N_INTF_A)
      REAL::ROWDEPC_INTFA(MAX_LBCROWS,INTF_LEN2_ROWDEPC,N_INTF_A)
      REAL::COLDEPC_INTFA(MAX_LBCROW_LENGTH,INTF_LEN2_COLDEPC,N_INTF_A)
 
      ! Eta Values for each area
      REAL :: LBC_ETA_THETA (MAX_INTF_MODEL_LEVELS+1,N_INTF_A)
      REAL :: LBC_ETA_RHO   (MAX_INTF_MODEL_LEVELS  ,N_INTF_A)
! TYPINFA end

  INTEGER                                                       &
     ind_im,                                                    &
                                ! IN Internal model indicator
     meanlev,                                                   &
                                ! INOUT  Mean level indicator
     icode             ! OUT return code; successful=0, error> 0


  CHARACTER(LEN=80)                                                  &
     cmessage          ! OUT Error message if ICODE > 0

!      Common blocks

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

!LL  Comdeck: CMEANCTL -------------------------------------------------
!LL
!LL  Purpose: Contains control information local to the climate meaning
!LL           control subroutines INITMEAN, MEANCTL etc.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Logical components covered: C500
!LL
!LL  Project task: C5
!LL
!LLEND -------------------------------------------------------------
      ! Length of largest IO buffer needed (A/O)
      INTEGER :: IBUFLEN(2)
!
      COMMON /CMEANCTL/ IBUFLEN
!
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

!      Local variables and arrays


  INTEGER                                                           &
     head_len                                                       &
                                ! Length of each record in header
     ,head_size              ! Size of header on disk

  PARAMETER(                                                        &
     head_len=3                                                     &
     )

  INTEGER                                                           &
     head_out(head_len,totitems) ! Header info for output ps file

  INTEGER                                                           &
     ifind,                                                     &
                                ! Loop counter
     step_dumps,                                                &
                                ! Timestep (in multiples of restart
                                !           dump frequency)
     nmeans,                                                    &
                                ! No. of means chosen (fixed, unless
                                ! there is a mean offset)
     residu,                                                    &
                                ! Reference for partial sum data
     means_total,                                               &
                                ! Number of mean dumps this timestep
     nmvals(4),                                                 &
                                ! Absolute meaning periods (in
                                ! multiples of restart dump frequency)
     ps_flag(4),                                                &
                                ! Flag for partial sum updating
     ft_ps(4,2),                                                &
                                ! Unit numbers for dumps of partial
                                ! sum data (read/write)
     periodlen,                                                 &
                                ! Length in days of current period N
     periodlendm,                                               &
                                ! Length in dumps of period 1 or
                                ! in days of current period N>1
     dumps_per_day                                              &
                                ! Number of restart dumps per day
     ,orig_decomp                                               &
                                ! Used to check for change in
     ,new_decomp                                                &
                                ! decomposition
     ,tag                                                       &
                                ! Stash tag
     ,ptd1                                                      &
                                ! Pointer to D1_ADDR information
     ,modnum                                                    &
                                ! Pointer to D1_ADDR submodel
     ,ptl                                                       &
                                ! Pointer to LOOKUP information
     ,global_length                                             &
                                ! Length of global field
     ,citems              ! Counter for no of objects to mean

  INTEGER                                                           &
     index_read,                                                &
                                ! Specific index no for reading
     index_write          ! Specific index no for writing
  INTEGER                                                           &
     ft_read,                                                   &
                                ! Unit number for partial sum read
     ft_write,                                                  &
                                ! Unit number for partial sum write
     ft_to_close,                                                 &
                                ! Unit number of ps file which is 
                                ! soon to be closed
     len_psname

  INTEGER internal_model
  INTEGER                                                           &
     reinit_steps                                               &
                                ! dummy input for GET_NAME
     ,      d1_addr_submodel_id  ! Submodel number in D1_ADDR
  LOGICAL                                                           &
     lmeaninc             ! increment MEANS_TOTAL or not

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

  INTEGER                                                           &
     ie                                                              &
                                ! loop counter over items
     , maxsize                                                       &
                                ! maximum dump output length
     ,lmaxsize                                                       &
                                ! Maximum output length per level
     ,totsize                                                        &
                                ! Size required for partial sum dumps
     ,levsize                ! Size of each level
  INTEGER dump_freq
  INTEGER mean_freq1,mean_freq2,mean_freq3,mean_freq4
  INTEGER stashoutputtime,mean1freqts

  INTEGER :: i_dump_period_days ! Dump period in days
  INTEGER :: i_dump_period_secs ! Dump period in seconds
  INTEGER :: i_sec_of_today     ! Number of seconds in day of current timestep

  !      Character data

  CHARACTER(LEN=filenamelength) :: psname_read,psname_write,psname_to_close
  CHARACTER(LEN=1)              :: letter_3   ! dummy input for GET_NAME
  CHARACTER(LEN=filenamelength) :: filename
  ! Name of pipe to server

  INTEGER :: copy_len                  ! length of temp. copy
  REAL, ALLOCATABLE :: d1_save(:)      ! temporary save space

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MEANCTL',zhook_in,zhook_handle)

  cmessage=' '
  letter_3='a'
  reinit_steps=0

  !      Define mode of use for means program

  IF(ind_im == 1)THEN
    WRITE(6,'(A)')'MEANCTL: ***** Called in ATMOSPHERIC mode *****'
  END IF
!
! ----------------------------------------------------------------------
!      Find out which mean datasets need to be created
!      on this timestep (if any) and set MEANS_TOTAL accordingly
! ----------------------------------------------------------------------
!

  icode=0
  means_total=0

!      Initially check validity of call to subroutine

  IF(dumpfreqim(ind_im) == 0)THEN ! Is mean dump production off?
    icode=1
    cmessage='MEANCTL: Invalid call to subroutine'
    WRITE(6,'(A,I3,A,I5)') 'meanctl: dumpfreq(',ind_im,')= ', &
       dumpfreqim(ind_im)
    GO TO 9999
  ELSE IF(MOD(stepim(ind_im),dumpfreqim(ind_im)) /= 0)THEN
    icode=2        ! This is not a dumping timestep
    cmessage='MEANCTL: Incorrect timestep to call subroutine'
    WRITE(6,'(A)') 'MEANCTL: STEP is not a multiple of DUMPFREQ'
    WRITE(6,'(A,I3,A,I5,A,I3,A,I5)')          &
       '         step(',   ind_im,            &
       ')= ',              stepim(ind_im),    &
       ' dumpfreq(',       ind_im,            &
       ')= ',              dumpfreqim(ind_im)
    GO TO 9999
  ELSE
    step_dumps = (stepim(ind_im)/dumpfreqim(ind_im))+             &
       offset_dumpsim(ind_im)
  END IF

!      Pick up number of means chosen from history file (MEAN_NUMBERim)
!      or number determined by the offset from the reference time whilst
!      the staggered start of means production unwinds  (MEAN_OFFSETim)

  DO ifind=1,mean_numberim(ind_im)
    IF(ifind == 1)nmvals(ifind) = meanfreqim(ifind,ind_im)
    IF(ifind >  1)nmvals(ifind) = meanfreqim(ifind,ind_im)*       &
       nmvals(ifind-1)
  END DO

  IF (lclimrealyr) THEN
    dumps_per_day=(24*3600*steps_per_periodim(ind_im))/             &
       (dumpfreqim(ind_im)*secs_per_periodim(ind_im))
  END IF

  IF(mean_offsetim(ind_im) == mean_numberim(ind_im))THEN
    nmeans=mean_numberim(ind_im)
  ELSE  
! there is a non-zero offset so increment mean_offset
! when run is partway into the current mean period.
    IF (lclimrealyr) THEN
! Turn on monthly/seasonal/annual meaning once past start of 
! first full month/season/year.
! initmean would have turned it on already if run started at 
! start of a month/season/year
      i_dump_period_days = dumpfreqim(a_im)/steps_per_periodim(a_im)
      i_dump_period_secs = (86400*dumpfreqim(a_im))/steps_per_periodim(a_im) - &
         86400 * i_dump_period_days
      i_sec_of_today = i_hour * 3600 + i_minute * 60 +i_second
      
      IF (i_sec_of_today == i_dump_period_secs .AND.          &
         (i_day - i_dump_period_days) == 1) THEN
        IF (mean_offsetim(ind_im)  ==  0) THEN
          
          mean_offsetim(ind_im) = 1 ! months
        END IF
      
        IF (mean_offsetim(ind_im)  ==  1) THEN
          IF (MOD((i_month-mean_reftimeim(2,ind_im)),3) == 0) THEN
            mean_offsetim(ind_im) = 2 ! seasons
          END IF
        END IF
        
        IF (mean_offsetim(ind_im)  ==  2) THEN
          IF (i_month == mean_reftimeim(2,ind_im)) THEN
            mean_offsetim(ind_im) = 3 ! years
          END IF  ! of test on ifind, i_day, etc.
        END IF
      END if
    ELSE  ! for 360d year, use array of no. of mean dumps
      DO ifind=mean_offsetim(ind_im)+1,mean_numberim(ind_im)
        IF(step_dumps <  0)THEN  ! mean ref time not reached yet
          residu=1-nmvals(ifind)
        ELSE                     ! mean ref time has been passed
          residu=1
        END IF
        IF(MOD(step_dumps,nmvals(ifind)) == residu)THEN
          mean_offsetim(ind_im)=mean_offsetim(ind_im)+1
        END IF
      END DO
    END IF
    nmeans=mean_offsetim(ind_im)
    WRITE(6,'(A,I3,A,I5)')' mean_offset(',ind_im,')=',mean_offsetim(ind_im)
  END IF  ! end of setting NMEANS
!
!       If no processing is required (because of staggered
!       start in means production) then skip to end of subroutine
!
  IF(nmeans == 0)THEN
    icode=-1
    cmessage='MEANCTL: No accumulation/meaning done this step'
    WRITE(6,'(A)')'MEANCTL: No accm/meaning due to staggered start'
    GO TO 9999
  END IF

!
!       Output message whilst stagger unwinds
!
  DO ifind=1,mean_numberim(ind_im)
    IF(nmeans <  ifind)THEN
      WRITE(6,'(A,I1,A,A)')'meanctl: period_',ifind,' mean not activated',   &
         ' because of staggered start in means production'
    END IF
  END DO

  DO ifind=1,nmeans
    ps_flag(ifind)=0
    lmeaninc=.FALSE. ! Initialise to avoid false positive results

! Find out if the end of any meaning period has been reached. If so,
! increment MEANS_TOTAL
    IF (lclimrealyr) THEN
      IF ((i_day == 1) .AND. (i_hour == 0) .AND.         &
         (i_minute == 0) .AND. (i_second == 0)) THEN
! Climate meaning periods for Gregorian calendar must always
! be aligned with start of the month
        IF ((ifind == 1) .AND.                       &
           (stepim(ind_im) > steps_per_periodim(ind_im))) THEN
          lmeaninc=.TRUE.                               ! monthly
        ELSE IF ((ifind == 2) .AND.                  &
           ! seasonal
           (stepim(ind_im) > steps_per_periodim(ind_im)) .AND.   &
           MOD((i_month-mean_reftimeim(2,ind_im)),3) == 0) THEN
          lmeaninc=.TRUE.
        ELSE IF ((ifind == 3) .AND.                  &
           ! annual
           (stepim(ind_im) > steps_per_periodim(ind_im)) .AND.   &
           MOD((i_month-mean_reftimeim(2,ind_im)),12) == 0) THEN
          lmeaninc=.TRUE.
        END IF
      END IF
    ELSE 
      ! 360-day calendar
      IF (MOD(step_dumps,nmvals(ifind)) == 0) THEN
        lmeaninc=.TRUE.
      END IF
    END IF  !
    IF (lmeaninc) THEN
      means_total=means_total+1
      WRITE(6,'(A,I1,A)') 'meanctl: period_',ifind,    &
         ' mean to be created this timestep'
    END IF  ! end of IF test on lmeaninc

! Find out if run is one period(N-1) into the period(N), in which case
! set PS_FLAG(IFIND)=1 because there will not already be a file for
! that partial sum on the disk, so ACUMPS must get it from D1.

    IF (lclimrealyr) THEN
      IF (ifind  ==  1) THEN
        ! Convert dump periods into seconds
        i_dump_period_secs = (86400*dumpfreqim(a_im))/steps_per_periodim(a_im)
        ! How many seconds are we into the current month
        i_sec_of_today = (i_day - 1) * 86400 + i_hour * 3600 + &
           i_minute * 60 + i_second
        ! If equal, we are one mean dump period into the period 1 mean period
        IF (i_sec_of_today == i_dump_period_secs) THEN
          ps_flag(1)=1                              ! months
        END IF
      ELSE IF (ifind  ==  2 .AND.                             &
         ! First mean dump of a period 2 period is at the start of the month
         (i_day  ==  1) .AND. (i_hour  ==  0) .AND.           &
         (i_minute == 0) .AND. (i_second == 0) .AND.          &
         MOD((i_month-1-mean_reftimeim(2,ind_im)),3) == 0)    &
         THEN
        ps_flag(2)=1                                   ! seasons
      ELSE IF (ifind  ==  3 .AND.                             &
         ! First mean dump of a period 3 period is at the start of the month
         (i_day  ==  1) .AND. (i_hour  ==  0) .AND.           &
         (i_minute == 0) .AND. (i_second == 0) .AND.          &
         MOD(((i_month-3)-mean_reftimeim(2,ind_im)),12) == 0) &
         THEN
        ps_flag(3)=1                                   ! years
      END IF
    ELSE  ! for 360d year, check using array of no. of mean dumps
      IF(step_dumps <  0)THEN
        IF(ifind == 1)residu=1-nmvals(ifind)
        IF(ifind /= 1)residu=nmvals(ifind-1)-nmvals(ifind)
      ELSE
        IF(ifind == 1)residu=1
        IF(ifind /= 1)residu=nmvals(ifind-1)
      END IF
      IF(MOD(step_dumps,nmvals(ifind)) == residu)THEN
        ps_flag(ifind)=1
      END IF
    END IF  ! end of test on lclimrealyr and setting of PS_FLAG


  END DO       ! end of loop over IFIND from 1 to NMEANS

!      Set up unit numbers for partial sum dumps

  IF (ind_im == 1) THEN
    ft_ps(1,1)=23
    ft_ps(1,2)=24
  END IF
  DO ifind=2,4
    ft_ps(ifind,1)=25
    ft_ps(ifind,2)=26
  END DO
  ! Units must still be alternated as 3 may be open at a time
  ft_ps(3,1)=45
  ft_ps(3,2)=46

  ! check STASH items output time against restart dump period
  dump_freq=dumpfreqim(ind_im)     ! Unit: timestep
  mean_freq1=meanfreqim(1,ind_im)
  mean_freq2=meanfreqim(2,ind_im)
  mean_freq3=meanfreqim(3,ind_im)
  mean_freq4=meanfreqim(4,ind_im)
  ! Calculating mean_freq1 in unit of timestep
  mean1freqts=mean_freq1*dump_freq
  ! Create header for partial sum file
  IF (ind_im == 1) THEN
    ! Set pointer to submodel info in D1_ADDR
    modnum=submodel_for_sm(ind_im)
    citems=0
    DO ie=1,totitems
      tag=stlist(st_macrotag,ie)/1000
      ! Get pointer to element in D1_ADDR array
      ptd1=stlist(st_d1pos,ie)
      IF(tag /= 0.AND.                                              &
         stlist(s_modl,ie) == d1_addr(d1_imodl,ptd1,modnum))THEN
        stashoutputtime=stlist(st_freq_code,ie)  ! Unit: timestep
        IF (stashoutputtime >= mean1freqts) THEN
          WRITE(6,'(A)')'ERROR: STASH output time frequency greater'
          WRITE(6,'(A)')'       than climate mean period frequency!'
          WRITE(6,'(A,I5,A,A,I5,A)')                           &
             'stash output frequency=',stashoutputtime,'(ts)'  &
             ,'  period 1 period=',mean1freqts,'(ts)'
          icode=4
          cmessage='MEANCTL: ERROR - STASHoutput > CMEAN period'
          GO TO 9999
        END IF
! Tagged for meaning and submodel information matches.
! Counter for no of variables that will be processed
        citems=citems+1
        ptl=d1_addr(d1_lookup_ptr,ptd1,modnum)
        global_length=stlist(st_dump_level_output_length,ie)
        global_length=global_length+1
        head_out(1,citems)=global_length
        head_out(2,citems)=1
        head_out(3,citems)=(((global_length+io_field_padding)       &
           /io_field_padding)*io_field_padding)
        IF (MOD((a_lookup(lbpack,ptl)),10)  ==  2) THEN
          IF (a_lookup(data_type,ptl)  ==  1) THEN
            ! Data to be packed so reduce data length accordingly
            global_length=(global_length+1)/2
            head_out(1,citems)=global_length
            head_out(2,citems)=2
            head_out(3,citems)=(((global_length+io_field_padding)   &
               /io_field_padding)*io_field_padding)
          END IF
        END IF
      END IF
    END DO
  END IF

! Find the largest output field size to dimension I/O buffer array
! and the total size of storage required to set file size
  lmaxsize=1
  head_size=citems*head_len+io_field_padding+2
  head_size=(head_size/io_field_padding)*io_field_padding
  totsize=head_size
  citems=0

  DO ie=1,totitems
    tag=stlist(st_macrotag,ie)/1000
    ! Get pointer to element in D1_ADDR array
    ptd1=stlist(st_d1pos,ie)
    IF (tag /= 0.AND.                                             &
       stlist(s_modl,ie) == d1_addr(d1_imodl,ptd1,modnum)) THEN
      ! Tagged for meaning and submodel information matches.
      citems=citems+1
      lmaxsize=MAX(lmaxsize,                                      &
         stlist(st_dump_level_output_length,ie))
      totsize=totsize+                                            &
         head_out(3,citems)*d1_addr(d1_no_levels,ptd1,modnum)
    END IF
  END DO
  ! Adding 1 extra space allowing for checksum
  lmaxsize=lmaxsize+1

  lmaxsize=((lmaxsize+io_field_padding)/io_field_padding)             &
     *io_field_padding

  ! Dump-packed fields will be rounded up to io_field_padding. Therefore,
  ! when unpacked, they will be rounded up to 2*io_field_padding. Add
  ! more space to account for this:
  lmaxsize = lmaxsize + io_field_padding

!
!  *********************************************************************
!                   LOGICAL SUB-PROCESS C51
!      Start of default process: updating period_1 partial sum data
!  *********************************************************************
!

  index_read=run_meanctl_indicim(1,ind_im)
  index_write=3-run_meanctl_indicim(1,ind_im)
  ft_read=ft_ps(1,index_read)
  ft_write=ft_ps(1,index_write)

  !      Temporary check on unit numbers

  IF (ps_flag(1) /= 1) THEN
    WRITE(6,'(A,I5)') 'Period_1 data read from unit number ',ft_read
  END IF
  WRITE(6,'(A,I5)') 'Period_1 data written to unit number ',ft_write

!
!                            STEP 1
!      Update or create period_1 partial sum data and write out
!      to period_1 partial sum dump
!

  orig_decomp=current_decomp_type
  new_decomp=orig_decomp

  IF (ind_im==1 .AND. orig_decomp/=decomp_standard_atmos) THEN
    new_decomp = decomp_standard_atmos
  ELSE IF (ind_im==2 .AND. orig_decomp/=decomp_nowrap_ocean) THEN
    new_decomp = decomp_nowrap_ocean
  END IF
  IF (new_decomp  /=  orig_decomp) THEN
    CALL change_decomposition(new_decomp,icode)
    IF (icode /= 0) THEN
      WRITE(6,'(A)')'ERROR : MEANCTL'
      WRITE(6,'(A,I5)')'Failed to change decomposition to ',new_decomp
      cmessage='MEANCTL1: Failed to change decomposition'
      GO TO 9999
    END IF
  END IF
  IF (means_total >  0) THEN
!
!                  STEP 1
!      Copy instantaneous dump to temporary memory space
!
!  Strictly, only climate mean tagged diagnostics need to be saved.
    IF (ind_im == 1) THEN
      copy_len = a_len_data+(model_levels+1)*theta_field_size
      ALLOCATE( d1_save(copy_len) )
      d1_save(1:copy_len) = d1(1:copy_len)
    END IF


  END IF

!--preset the file lengths prior to the open
!--pass the new file length to the I/O routines
  CALL set_dumpfile_length(ft_read , totsize)
  CALL set_dumpfile_length(ft_write, totsize)

!      Open input and output partial sum files (preassigned names)

  CALL model_file_open(ft_read,ft_environ(ft_read),                     &
     len_ft_envir(ft_read),1,0,icode)
  IF (icode /= 0) GO TO 9999

  CALL model_file_open(ft_write,ft_environ(ft_write),                   &
     len_ft_envir(ft_write),1,0,icode)
  IF (icode /= 0) GO TO 9999
  !--zero the file lengths
  CALL set_dumpfile_length(ft_read , 0)
  CALL set_dumpfile_length(ft_write, 0)

  IF (ind_im == 1) THEN
    d1_addr_submodel_id = submodel_for_sm(1)
    ! DEPENDS ON: acumps
    CALL acumps(                                                  &
       no_obj_d1(d1_addr_submodel_id),                            &
       d1_addr(1,1,d1_addr_submodel_id),                          &
       a_len_data,d1,                                             &
       lmaxsize,means_total,                                      &
       ps_flag(1),ft_read,ft_write,lclimrealyr,meanlev,           &
       i_month,i_year,                                            &
       head_out,head_len,head_size,                               &
       stepim(ind_im),citems,a_fixhd(12),                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       icode,cmessage)
  END IF

!      Check return code from ACUMPS

  IF (icode /= 0) THEN
    GO TO 9999
  END IF

!      Close input and output partial sum files


  CALL model_file_close(ft_read,ft_environ(ft_read),                    &
     len_ft_envir(ft_read),0,0,icode)

  CALL model_file_close(ft_write,ft_environ(ft_write),                  &
     len_ft_envir(ft_write),0,0,icode)

!      Update RUN_MEANCTL_INDICim for period_1 data

  run_meanctl_indicim(1,ind_im)=3-run_meanctl_indicim(1,ind_im)


!  *********************************************************************
!      End of default process: updating period_1 partial sum data
!  *********************************************************************

  IF (means_total >  0) THEN

!  *********************************************************************
!                   LOGICAL SUB-PROCESS C52
!      Start of means processing and updating of subsequent
!      partial sum dump
!  *********************************************************************

    DO meanlev=1,means_total

      index_read=run_meanctl_indicim(meanlev,ind_im)
      ft_read=ft_ps(meanlev,index_read)

      !      Temporary check on unit number

      WRITE(6,'(a,I1,a,I5)') 'period_',meanlev,' data read:unit number ',ft_read
!
!                              STEP 2
!      Generate period_N time-meaned data and store in main data block
!
!       If real-period meaning selected, find length of current period
      IF (lclimrealyr) THEN
        ! DEPENDS ON: setperlen
        CALL setperlen (meanlev,i_month,i_year,periodlen)
        IF (meanlev == 1) THEN ! divisor only needs to be in terms
          ! of restart dumps for Period_1
          periodlendm=periodlen*dumps_per_day
        ELSE
          periodlendm=periodlen
        END IF
      END IF
      IF (ind_im == 1) THEN

        !      Open input partial sum file (preassigned or calculated name)

        !--pass the new file length to the I/O routines
        CALL set_dumpfile_length(ft_read , totsize)
        IF (meanlev == 1) THEN

          CALL model_file_open(ft_read,ft_environ(ft_read),         &
             len_ft_envir(ft_read),1,0,icode)
          IF (icode /= 0) GO TO 9999
        ELSE
          ! DEPENDS ON: get_name
          CALL get_name(expt_id,job_id,ind_im,meanlev,index_read,   &
             reinit_steps,'s',letter_3,                  &
             model_status,time_convention,               &
             0.0,psname_read,icode,cmessage,lcal360)
          IF (icode >  0) GO TO 9999
          len_psname=LEN(psname_read)

          CALL model_file_open(ft_read,psname_read,len_psname,1,1,    & 
             icode)
          IF (icode /= 0) GO TO 9999
        END IF
        !--unset the file length in the I/O routines
        CALL set_dumpfile_length(ft_read , 0)

        IF (lclimrealyr) THEN  ! Real-period meaning selected
          ! DEPENDS ON: meanps
          CALL meanps(                                                &
             no_obj_d1(d1_addr_submodel_id),                          &
             d1_addr(1,1,d1_addr_submodel_id),                        &
             a_len_data,d1,                                           &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          periodlendm                                                 &
             )
        ELSE            ! 360d year meaning selected
          ! DEPENDS ON: meanps
          CALL meanps(                                              &
             no_obj_d1(d1_addr_submodel_id),                        &
             d1_addr(1,1,d1_addr_submodel_id),                      &
             a_len_data,d1,                                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          meanfreqim(meanlev,ind_im)                                &
             )

        END IF            ! end of test on lclimrealyr
      END IF            ! end of test on IND_IM == 1

      !      Check return code from MEANPS

      IF (icode /= 0) THEN
        WRITE(6,'(A,I3)') 'MEANCTL: Error at mean period ',meanlev
        GO TO 9999
      END IF

      psname_to_close=psname_read

      CALL setpos(ft_read,0,icode)
      ft_to_close=ft_read
!
!                           STEP 3.1
!      Calculate mean diagnostics and extract PPfields from mean data
!
      IF (ind_im == 1) THEN
        ! Find the largest output field size to dimension I/O buffer array
        maxsize=1

        DO ie=1,totitems
          tag=stlist(st_macrotag,ie)/1000
          IF (MOD(tag/(2**(meanlev-1)),2)  ==  1) THEN
            maxsize=MAX(maxsize,                                    &
               stlist(st_dump_level_output_length,ie))
          END IF
        END DO
        ! Extract mean diagnostics
        ! DEPENDS ON: meandiag
        CALL meandiag (                                             &
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
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
        ind_im,meanlev,pp_len2_meanim(meanlev,ind_im),step_dumps,  &
           nmvals(meanlev),                                       &
           maxsize,                                               &
           icode,cmessage)
      END IF
      IF (icode >  0) GO TO 9999
!
!                           STEP 4
!      Check to see if period_N+1 partial sum data needs to updated
!      or created.
!      If so, proceed and write out to period_N+1 partial sum dump
!

      IF (meanlev /= nmeans) THEN

        index_read=run_meanctl_indicim(meanlev+1,ind_im)
        index_write=3-run_meanctl_indicim(meanlev+1,ind_im)
        ft_read=ft_ps(meanlev+1,index_read)
        ft_write=ft_ps(meanlev+1,index_write)

        !      Temporary check on unit numbers

        IF (ps_flag(meanlev+1) /= 1) THEN
          WRITE(6,'(A,I1,A,I5)')                                &
          'period_',meanlev+1,' data read:unit number ',ft_read
        END IF
        WRITE(6,'(A,I1,A,I5)') 'period_',meanlev+1,             &
           ' data written:unit number ',&
           ft_write

        IF (ind_im == 1) THEN

          !      Open input and output partial sum files (calculated names)

          ! DEPENDS ON: get_name
          CALL get_name(expt_id,job_id,ind_im,meanlev+1,index_read, &
             reinit_steps,'s',letter_3,                  &
             model_status,time_convention,               &
             0.0,psname_read,icode,cmessage,lcal360)
          IF (icode >  0) GO TO 9999
          !--pass the new file length to the I/O routines
          CALL set_dumpfile_length(ft_read , totsize)
          len_psname=LEN(psname_read)

          CALL model_file_open(ft_read,psname_read,len_psname,1,1,icode)
          IF (icode /= 0) GO TO 9999
          !--unset the file length in the I/O routines
          CALL set_dumpfile_length(ft_read , 0)

          ! DEPENDS ON: get_name
          CALL get_name(expt_id,job_id,ind_im,meanlev+1,index_write,&
             reinit_steps,'s',letter_3,                  &
             model_status,time_convention,               &
             0.0,psname_write,icode,cmessage,lcal360)
          IF (icode >  0) GO TO 9999

          !--pass the new file length to the I/O routines
          CALL set_dumpfile_length(ft_write, totsize)
          len_psname=LEN(psname_write)

          CALL model_file_open(ft_write,psname_write,len_psname,1,1,icode)
          IF (icode /= 0) GO TO 9999
          !--unset the file length in the I/O routines
          CALL set_dumpfile_length(ft_write, 0)

          d1_addr_submodel_id = submodel_for_sm(1)
          ! DEPENDS ON: acumps
          CALL acumps(                                                       &
             no_obj_d1(d1_addr_submodel_id),                                 &
             d1_addr(1,1,d1_addr_submodel_id),                               &
             a_len_data,d1,                                                  &
             lmaxsize,means_total,                                           &
             ps_flag(meanlev+1),ft_read,ft_write,                            &
             lclimrealyr,meanlev,i_month,i_year,                             &
             head_out,head_len,head_size,                                    &
             stepim(ind_im),citems,a_fixhd(12),                              &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          icode,cmessage)
        END IF

!      Check return code from ACUMPS

        IF (icode /= 0) THEN
          GO TO 9999
        END IF


!      Update RUN_MEANCTL_INDICim for period_N+1 data

        run_meanctl_indicim(meanlev+1,ind_im)=                      &
           3-run_meanctl_indicim(meanlev+1,ind_im)

      END IF

!      Decide disposition of period_N+1 partial sum dumps
!      NB: for restartability it is NOT safe to delete partial 
!      sum files.

      IF (meanlev /= nmeans) THEN
        len_psname=LEN(psname_read)

        CALL setpos(ft_read,0,icode)

        CALL model_file_close(ft_read,psname_read,len_psname,1,0, &
           icode)

        CALL setpos(ft_write,0,icode)

        len_psname=LEN(psname_write)
        CALL model_file_close(ft_write,psname_write,len_psname,1,0, &
           icode)
      END IF

!      Decide disposition of remaining period_N partial sum dump
!      NB: for restartability it is NOT safe to delete partial sum
!      files.

      IF (meanlev >= 2) THEN
        len_psname=LEN(psname_to_close)

        CALL model_file_close(ft_to_close,psname_to_close,len_psname,1,0, &
           icode)
      ELSE

        CALL model_file_close(ft_to_close,ft_environ(ft_to_close), &
           len_ft_envir(ft_to_close),0,0,icode)
      END IF

    END DO ! Loop over meanlev
!
!                     STEP 6
!      Read back instantaneous dump from temporary memory space
!
    IF (ind_im == 1) THEN

      d1( 1:copy_len ) = d1_save( 1:copy_len )
      DEALLOCATE( d1_save )

    END IF
  END IF
!
!  *********************************************************************
!      End of means processing and updating of subsequent
!      partial sum dumps
!  *********************************************************************
!

  orig_decomp=current_decomp_type
  new_decomp=orig_decomp
  IF (ind_im == 1.AND.orig_decomp /= decomp_standard_atmos)THEN
    new_decomp=decomp_standard_atmos
  ELSE IF (ind_im==2.AND.orig_decomp/=decomp_standard_ocean) THEN
    new_decomp=decomp_standard_ocean
  END IF
  IF (new_decomp  /=  orig_decomp) THEN
    CALL change_decomposition(new_decomp,icode)
    IF (icode /= 0) THEN
      WRITE(6,'(A)')'ERROR : MEANCTL'
      WRITE(6,'(A,I5)')'Failed to change decomposition to ',new_decomp
      cmessage='MEANCTL1: Failed to change decomposition'
    END IF
  END IF

9999 CONTINUE

!      Reset MEANLEV to zero

  meanlev=0
  IF (lhook) CALL dr_hook('MEANCTL',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE meanctl
