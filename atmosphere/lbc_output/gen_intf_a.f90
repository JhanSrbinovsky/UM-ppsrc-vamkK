! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generates Atmosphere Lateral Boundary Conditions
!
! Subroutine Interface:

      SUBROUTINE GEN_INTF_A (                                           &
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
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
     &  JINTF,NFTOUT,                                                   &
     &  d1_array,                                                       &
     &  INTERNAL_MODEL,                                                 &
     &  ndustbin_in,                                                    & 
     &  ndustbin_out                                                    &
     & )

      USE Submodel_Mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE filenamelength_mod, ONLY:                                     &
          filenamelength
      USE IO
      USE atm_fields_bounds_mod
      USE dustbin_conversion_mod, ONLY: convert_dust_six_to_two,        &
                                        convert_dust_two_to_six
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      USE rimtypes
      USE lbc_mod
      USE lookup_addresses
      USE dust_parameters_mod, ONLY:             l_dust,                &
           l_dust_div1,       l_dust_div2,       l_dust_div3,           &
           l_dust_div4,       l_dust_div5,       l_dust_div6
      USE um_input_control_mod,  ONLY:                                  &
                              l_dust_div1_lbc,   l_dust_div1_lbc_out,   &
                              l_dust_div2_lbc,   l_dust_div2_lbc_out,   &
                              l_dust_div3_lbc,   l_dust_div3_lbc_out,   &
                              l_dust_div4_lbc,   l_dust_div4_lbc_out,   &
                              l_dust_div5_lbc,   l_dust_div5_lbc_out,   &
                              l_dust_div6_lbc,   l_dust_div6_lbc_out,   &
           l_so2,             l_so2_lbc,         l_so2_lbc_out,         &
           l_dms,             l_dms_lbc,         l_dms_lbc_out,         &
           l_so4_aitken,      l_so4_aitken_lbc,  l_so4_aitken_lbc_out,  &
           l_so4_accu,        l_so4_accu_lbc,    l_so4_accu_lbc_out,    &
           l_so4_diss,        l_so4_diss_lbc,    l_so4_diss_lbc_out,    &
           l_nh3,             l_nh3_lbc,         l_nh3_lbc_out,         &
           l_soot_new,        l_soot_new_lbc,    l_soot_new_lbc_out,    &
           l_soot_agd,        l_soot_agd_lbc,    l_soot_agd_lbc_out,    &
           l_soot_cld,        l_soot_cld_lbc,    l_soot_cld_lbc_out,    &
           l_bmass_new,       l_bmass_new_lbc,   l_bmass_new_lbc_out,   &
           l_bmass_agd,       l_bmass_agd_lbc,   l_bmass_agd_lbc_out,   &
           l_bmass_cld,       l_bmass_cld_lbc,   l_bmass_cld_lbc_out,   &
           l_ocff_new,        l_ocff_new_lbc,    l_ocff_new_lbc_out,    &
           l_ocff_agd,        l_ocff_agd_lbc,    l_ocff_agd_lbc_out,    &
           l_ocff_cld,        l_ocff_cld_lbc,    l_ocff_cld_lbc_out,    &
           l_nitr_acc,        l_nitr_acc_lbc,    l_nitr_acc_lbc_out,    &
           l_nitr_diss,       l_nitr_diss_lbc,   l_nitr_diss_lbc_out,   &
           l_soot,            l_ocff,                                   &
           lcal360

      USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain 
      USE cloud_inputs_mod, ONLY: l_pc2
      USE murk_inputs_mod,  ONLY: l_murk
      USE nlstcall_mod, ONLY : ft_steps, &
                               ltimer

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

!
! Description:
!   <Say what this routine does>
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

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

      Integer  ::  jintf           !  Index to interface area
      Integer  ::  nftout          !  Unit number for interface area
      Integer  ::  internal_model  !  Internal Model Number

      Real,Target ::  d1_array(*)     !  D1 array with prognostic data

      CHARACTER(LEN=80) CMESSAGE ! Error message if ICODE>0

      Integer :: i,j,var   ! Loop indices
      Integer :: lbc_num
      Integer :: lookup_start
      Integer :: len_io
      Integer :: ntime
      Integer :: len_ppname
      Integer :: im_ident      !   Internal model identifier
      Integer :: len_data
      Integer :: len_data_p
      Integer :: len_data_u
      Integer :: len_data_v
      Integer :: len_data_max

      Real A_IO


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

      INTEGER, INTENT(IN) :: ndustbin_in, ndustbin_out

!*---------------------------------------------------------------------
! Local parameters:                                    

       ! lbc stash section numbers                     
       Integer,           Parameter :: Sect31 = 31000
       Integer,           Parameter :: Sect36 = 36000
       Integer,           Parameter :: Sect37 = 37000

! Local variables:                                                   

       integer ifld, ihalo, iside
       integer intf_halosize(2,Nhalo_max)
       integer lbc_row_len, lbc_nrows
       integer len_data_usect
       integer len_data_vsect
       integer len_data_sect
       Integer var1
       Integer rimwidth

       Integer Sect_Prog(intf_lookupsa)
       Integer Item_Prog(intf_lookupsa)
       Integer IPt_D1(intf_lookupsa)
       Integer lbc_halo_type(intf_lookupsa)
       Integer lbc_fld_type (intf_lookupsa)
       Integer lbc_rim_type (intf_lookupsa)
       Integer lbc_level_type(intf_lookupsa)
       Integer lbc_first_level(intf_lookupsa)
       Integer lbc_last_level(intf_lookupsa)
       Integer lbc_levels   (intf_lookupsa)

       Integer lbc_src_halo_type (intf_lookupsa)
       Integer lbc_src_fld_type  (intf_lookupsa)
       Integer lbc_src_levels    (intf_lookupsa)
       Integer lbc_src_level_type(intf_lookupsa)
       Integer lbc_src_first_level(intf_lookupsa)
       Integer lbc_src_last_level(intf_lookupsa)

       Integer lbc_fld_type_interp
       Integer   ::  ErrorStatus
       Character (Len=*), Parameter :: RoutineName= 'Gen_Intf_A'
       Real, Dimension (:,:),    Allocatable :: LBC_Data
       Real, Dimension (:),      Allocatable :: LBC_Data_uv
       Real, Allocatable :: LBC_Coeff1 (:)
       Real, Allocatable :: LBC_Coeff2 (:)
       Real, Allocatable :: Lambda_p_in(:)
       Real, Allocatable :: Phi_p_in(:)
       Real, Allocatable :: Lambda_u_in(:)
       Real, Allocatable :: Phi_v_in(:)
       Real, Allocatable :: Lambda_in(:)
       Real, Allocatable :: Phi_in(:)

       Logical u_or_v
!      -------------------------
       Integer local_row_length
       Integer local_rows
       Integer local_row_length_v
       Integer local_rows_v
       Integer source_halo_x
       Integer source_halo_y
       Integer source_fld_type
       Integer source_halo_type
       Integer source_levels

       Real    source_delta_lat
       Real    source_delta_long
       Real    source_first_lat
       Real    source_first_long
       Real    source_pole_lat
       Real    source_pole_long

       REAL, ALLOCATABLE :: lambda_source(:)
       REAL, ALLOCATABLE :: phi_source(:)

! -----------------------------------------
       Integer orog_global_row_length
       Integer orog_global_rows
       Integer orog_local_row_length
       Integer orog_local_rows
       Integer orog_halo_x
       Integer orog_halo_y
       Integer orog_fld_type
       Integer orog_halo_type

       Real    orog_first_lat
       Real    orog_first_long
       
       REAL, ALLOCATABLE :: lambda_source_orog(:)
       REAL, ALLOCATABLE :: phi_source_orog(:)
       
! -----------------------------------------
       Integer lbc_halo_x
       Integer lbc_halo_y
       
       INTEGER :: lbc_halo_x_coords
       INTEGER :: lbc_halo_y_coords

       Real    lbc_delta_lat
       Real    lbc_delta_long
       Real    lbc_first_lat
       Real    lbc_first_long
       Real    lbc_pole_lat
       Real    lbc_pole_long

       REAL , POINTER :: source_field   (:)
       REAL , POINTER :: source_field_v (:)
       Integer :: lev, jrow
       Integer :: field_size, field_size_v

! Indexes for bottom left & right points in bi-linear interpolation
       Integer, dimension (:), allocatable :: lbc_index_bl
       Integer, dimension (:), allocatable :: lbc_index_br
! Weights for 4 points in bi-linear interpolation

       Real,    dimension (:), allocatable :: lbc_weights_tr
       Real,    dimension (:), allocatable :: lbc_weights_br
       Real,    dimension (:), allocatable :: lbc_weights_tl
       Real,    dimension (:), allocatable :: lbc_weights_bl

       Integer lbc_size, lbc_size_u, lbc_size_v, lbc_size_interp
       Integer max_levels_per_pe
       Integer gather_pe
       Integer i_uv
       Integer halo_x,halo_y,lbc_rows
       INTEGER :: row,pt,ipt,jvar
       Logical l_lbc_u, l_lbc_v, l_lbc_winds
       Logical :: l_calc_lbc_wts
! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! -----------------------------------------------------------
! Stash Codes for LBCs in Section 32 (and section 31 except tracers) 
!
      Integer, Parameter :: lbc_stashcode_orog    = 1 
      Integer, Parameter :: lbc_stashcode_u       = 2 
      Integer, Parameter :: lbc_stashcode_v       = 3 
      Integer, Parameter :: lbc_stashcode_w       = 4 
      Integer, Parameter :: lbc_stashcode_density = 5 
      Integer, Parameter :: lbc_stashcode_theta   = 6 
      Integer, Parameter :: lbc_stashcode_q       = 7 
      Integer, Parameter :: lbc_stashcode_qcl     = 8 
      Integer, Parameter :: lbc_stashcode_qcf     = 9 
      Integer, Parameter :: lbc_stashcode_exner   = 10 
      Integer, Parameter :: lbc_stashcode_u_adv   = 11 
      Integer, Parameter :: lbc_stashcode_v_adv   = 12 
      Integer, Parameter :: lbc_stashcode_w_adv   = 13 
      Integer, Parameter :: lbc_stashcode_qcf2    = 14 
      Integer, Parameter :: lbc_stashcode_qrain   = 15 
      Integer, Parameter :: lbc_stashcode_qgraup  = 16 
      Integer, Parameter :: lbc_stashcode_cf_bulk = 17 
      Integer, Parameter :: lbc_stashcode_cf_liquid = 18 
      Integer, Parameter :: lbc_stashcode_cf_frozen = 19 
      Integer, Parameter :: lbc_stashcode_murk      = 20 
      Integer, Parameter :: lbc_stashcode_free_tracer = 21 
      Integer, Parameter :: lbc_stashcode_ukca_tracer = 22 
      Integer, Parameter :: lbc_stashcode_dust_div1 = 23
      Integer, Parameter :: lbc_stashcode_dust_div2 = 24
      Integer, Parameter :: lbc_stashcode_dust_div3 = 25
      Integer, Parameter :: lbc_stashcode_dust_div4 = 26
      Integer, Parameter :: lbc_stashcode_dust_div5 = 27
      Integer, Parameter :: lbc_stashcode_dust_div6 = 28
      Integer, Parameter :: lbc_stashcode_so2      = 29
      Integer, Parameter :: lbc_stashcode_dms      = 30
      Integer, Parameter :: lbc_stashcode_so4_aitken = 31
      Integer, Parameter :: lbc_stashcode_so4_accu = 32
      Integer, Parameter :: lbc_stashcode_so4_diss = 33
      Integer, Parameter :: lbc_stashcode_nh3      = 35
      Integer, Parameter :: lbc_stashcode_soot_new = 36
      Integer, Parameter :: lbc_stashcode_soot_agd = 37
      Integer, Parameter :: lbc_stashcode_soot_cld = 38
      Integer, Parameter :: lbc_stashcode_bmass_new = 39
      Integer, Parameter :: lbc_stashcode_bmass_agd = 40
      Integer, Parameter :: lbc_stashcode_bmass_cld = 41
      Integer, Parameter :: lbc_stashcode_ocff_new = 42
      Integer, Parameter :: lbc_stashcode_ocff_agd = 43
      Integer, Parameter :: lbc_stashcode_ocff_cld = 44
      Integer, Parameter :: lbc_stashcode_nitr_acc = 45
      Integer, Parameter :: lbc_stashcode_nitr_diss = 46

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_Sec
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
         LBC_DT_Hour, LBC_DT_Min,  LBC_DT_Sec,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_Sec
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
         LBC_VT_Hour, LBC_VT_Min,  LBC_VT_Sec,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------
       Integer :: lbc_rim_size(Nrima_max)
       logical l_vi
       integer n_segs
       integer max_seg_size
       Integer :: prev_src_fld_type
       Integer :: prev_lbc_fld_type
       Integer :: prev_src_halo_type
       Integer :: prev_lbc_halo_type
       integer :: pt_lbc 
       integer :: LBCrow_len 
       integer :: LBCrows      
       Real    :: dlam_wk
       Real    :: dphi_wk 

       LOGICAL :: source_cyclic   ! T : Source Grid is cyclic
       LOGICAL :: source_rotated  ! T : Source Grid is rotated
       LOGICAL :: same_rotation   ! T : Grids have same rotation
       LOGICAL :: interp_on_p     ! T : If we are interpolating on p
       
       LOGICAL :: L_EG_dump       ! T : ENDGAME dump
       LOGICAL :: L_extra_surface_level ! T : extra surface level
                                        ! required for current LBC field
                                        ! (EG out)
       LOGICAL :: L_extra_surface_level_dump ! T : extra surface level
                                             ! present for current field
                                             ! in dump (EG in)

!      ------------------------
      CHARACTER(LEN=filenamelength) :: string         ! work array
      CHARACTER(LEN=14) PPNAME         ! boundary output filename

!*---------------------------------------------------------------------
!     Stash item and section numbers for interface fields
!     Any change to code generating and testing ITEM_INTFA should also
!     consider the corresponding use of ITEM_BOUNDA in INBOUND1/CHKLKBA1
      INTEGER ITEM_INTFA (INTF_LOOKUPSA), SECT_INTFA (INTF_LOOKUPSA)

      ! Stash codes for corresponding input lbc fields. It is these
      ! codes that are written in to the generated lbc file 
      Integer lbc_stash_codes(INTF_LOOKUPSA) 

      INTEGER :: first_dust_lbc_number
      REAL, ALLOCATABLE, TARGET :: dust_div3_out(:)
      REAL, ALLOCATABLE, TARGET :: dust_div4_out(:)
      REAL, ALLOCATABLE, TARGET :: dust_div5_out(:)
      REAL, ALLOCATABLE, TARGET :: dust_div6_out(:)
      
      ! Positions of points on src grid
      REAL :: src_dl_off(nfld_max)
      REAL :: src_dp_off(nfld_max)
 
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!*---------------------------------------------------------------------
!*
!L Internal structure:

       IF (lhook) CALL dr_hook('GEN_INTF_A',zhook_in,zhook_handle)
       ErrorStatus = 0
       CMESSAGE=' '
       NULLIFY(source_field)
       NULLIFY(source_field_v)

       source_cyclic  = .FALSE.
       source_rotated = .FALSE.
       same_rotation  = .FALSE.
       interp_on_p    = .FALSE.

       source_cyclic  = (A_FIXHD(4) < 3)  !  If Global or NH or SH
       source_rotated = (A_REALHD(5)/= 90.0 .or. A_REALHD(6) /= 0.0)
       
       L_EG_dump      = A_FIXHD(9) == 6
       
       ! Not set staggering information for source grid.
      ! src_dl_off(:) = RMDI
      ! src_dp_off(:) = RMDI

       SELECT CASE(a_fixhd(9))
       CASE(3)
         ! Arakawa C with P along the bottom and along the left.
         src_dl_off(fld_type_p) = 0.0
         src_dp_off(fld_type_p) = 0.0
         src_dl_off(fld_type_u) = 0.5
         src_dp_off(fld_type_u) = 0.0
         src_dl_off(fld_type_v) = 0.0
         src_dp_off(fld_type_v) = 0.5
         
       CASE(6)
         ! Arakawa C with v along the bottom and u along the left.
         src_dl_off(fld_type_p) = 0.0
         src_dp_off(fld_type_p) = 0.0
         src_dl_off(fld_type_u) = -0.5
         src_dp_off(fld_type_u) = 0.0
         src_dl_off(fld_type_v) = 0.0
         src_dp_off(fld_type_v) = -0.5
       END SELECT

! Set Data Time

       LBC_DT_Year  = A_FIXHD(21)
       LBC_DT_Month = A_FIXHD(22)
       LBC_DT_Day   = A_FIXHD(23)
       LBC_DT_Hour  = A_FIXHD(24)
       LBC_DT_Min   = A_FIXHD(25)
       LBC_DT_Sec   = A_FIXHD(26)
       LBC_DT_DayNo = A_FIXHD(27)

!      =======================================================
       Src_Grid (fld_type_p,1) = a_realhd(rh_baselat) +                   &
                                 src_dp_off(fld_type_p)*a_realhd(rh_deltaNS)
       Src_Grid (fld_type_p,2) = a_realhd(rh_baselong) +                  &
                                 src_dl_off(fld_type_p)*a_realhd(rh_deltaEW)
       Src_Grid (fld_type_u,1) = a_realhd(rh_baselat) +                   &
                                 src_dp_off(fld_type_u)*a_realhd(rh_deltaNS)
       Src_Grid (fld_type_u,2) = a_realhd(rh_baselong) +                  &
                                 src_dl_off(fld_type_u)*a_realhd(rh_deltaEW)
       Src_Grid (fld_type_v,1) = a_realhd(rh_baselat) +                   &
                                 src_dp_off(fld_type_v)*a_realhd(rh_deltaNS)
       Src_Grid (fld_type_v,2) = a_realhd(rh_baselong) +                  &
                                 src_dl_off(fld_type_v)*a_realhd(rh_deltaEW)



!      =======================================================
       ! Now set staggering for LBC grid

       LBC_Grid (fld_type_p,1) = Phi_intf_p(1,jintf)
       LBC_Grid (fld_type_p,2) = Lambda_intf_p(1,jintf)
       LBC_Grid (fld_type_p,3) = Intf_Row_Length(jintf)
       LBC_Grid (fld_type_p,4) = Intf_P_Rows(jintf)
 
       LBC_Grid (fld_type_u,1) = Phi_intf_p(1,jintf)
       LBC_Grid (fld_type_u,2) = Lambda_intf_u(1,jintf)

       IF (intf_l_eg_out(jintf)) THEN
         LBC_Grid (fld_type_u,3) = Intf_Row_Length(jintf)
         LBC_Grid (fld_type_u,4) = Intf_P_Rows(jintf)
       ELSE
         LBC_Grid (fld_type_u,3) = Intf_Row_Length(jintf) - 1
         LBC_Grid (fld_type_u,4) = Intf_P_Rows(jintf)
       END IF
       LBC_Grid (fld_type_v,1) = Phi_Intf_v(1,jintf)
       LBC_Grid (fld_type_v,2) = Lambda_Intf_p(1,jintf)

       IF (intf_l_eg_out(jintf)) THEN
         LBC_Grid (fld_type_v,3) = Intf_Row_Length(jintf)
         LBC_Grid (fld_type_v,4) = Intf_P_Rows(jintf) + 1
       ELSE
         LBC_Grid (fld_type_v,3) = Intf_Row_Length(jintf)
         LBC_Grid (fld_type_v,4) = Intf_P_Rows(jintf) - 1
       END IF
     
!      =======================================================

!      Set up intf_halosize for this area
       intf_halosize(1,1)=1
       intf_halosize(2,1)=1

       intf_halosize(1,2)=intf_exthalo_ew(jintf)
       intf_halosize(2,2)=intf_exthalo_ns(jintf)

       intf_halosize(1,3)=0
       intf_halosize(2,3)=0
!      =======================================================      

       LBCrow_len = MAX(LBC_Grid (fld_type_p,3),LBC_Grid (fld_type_u,3))
       LBCrows    = MAX(LBC_Grid (fld_type_p,4),LBC_Grid (fld_type_v,4))
! This is being set to the maximum halo size - later on the same
! variable is used store the variables halo size which could be different.

       IF (intf_l_eg_out(jintf)) THEN
       ! Need extra P points around outside edge of domain for
       ! interpolation back to U/V grids for ENDGAME
       ! The modified halo sizes are only to ensure the coordinate vectors
       ! have these extra points.
         lbc_halo_x = intf_exthalo_ew(jintf)+1
         lbc_halo_y = intf_exthalo_ns(jintf)+1
       ELSE
         lbc_halo_x = intf_exthalo_ew(jintf)
         lbc_halo_y = intf_exthalo_ns(jintf)
       END IF
       
       allocate (Lambda_p_in(1-lbc_halo_x: LBCrow_len + lbc_halo_x)) 
       allocate (Phi_p_in(1-lbc_halo_y: LBCrows + lbc_halo_y))
       allocate (Lambda_u_in(1-lbc_halo_x: LBCrow_len + lbc_halo_x)) 
       allocate (Phi_v_in(1-lbc_halo_y: LBCrows + lbc_halo_y))
       allocate (Lambda_in(1-lbc_halo_x: LBCrow_len + lbc_halo_x)) 
       allocate (Phi_in(1-lbc_halo_y: LBCrows + lbc_halo_y))
    
       ! Row length is always u and p are equal.
       do pt_lbc = 1, LBCrow_len
         Lambda_u_in(pt_lbc) = Lambda_intf_u(pt_lbc,jintf)
         Lambda_p_in(pt_lbc) = Lambda_intf_p(pt_lbc,jintf)
       end do
       ! Rows p and v are always different.
       DO pt_lbc = 1, NINT(lbc_grid (fld_type_p,4))
         phi_p_in(pt_lbc) = phi_intf_p(pt_lbc,jintf)
       END DO
       DO pt_lbc = 1, NINT(lbc_grid (fld_type_v,4))
         phi_v_in(pt_lbc) = phi_intf_v(pt_lbc,jintf)
       END DO 
             
! compute lambda and phi intervals dlam_wk and dphi_wk in the 
! halo region


       IF ( intf_l_var_lbc(jintf) ) THEN
         dlam_wk = Lambda_p_in(LBCrow_len)-Lambda_p_in(LBCrow_len-1) 
         dphi_wk = Phi_p_in(NINT(lbc_grid (fld_type_p,4))) &
                   -Phi_p_in(NINT(lbc_grid (fld_type_p,4))-1)
       Else       !regular grid
         dlam_wk =  Intf_ewspace(jintf)
         dphi_wk =  Intf_nsspace(jintf) 
       End if
       
       ! Lets add on any difference between p and v rows.
       DO pt_lbc = NINT(lbc_grid (fld_type_p,4))+1, LBCrows
         phi_p_in(pt_lbc) = phi_p_in(pt_lbc-1) + dphi_wk
       END DO
       DO pt_lbc = NINT(lbc_grid (fld_type_v,4))+1, LBCrows
         phi_v_in(pt_lbc) = phi_v_in(pt_lbc-1) + dphi_wk
       END DO

! compute lambda and phi in the halo region - how lambda and phi are calculated
! can affect how the interpolation weights are calculated since a small
! difference can make all the difference when finding the nearest point.  This
! is especially true near the greenwich meridian where a mod(lambda,360.) is
! taken to make all angles fall between [0,360.0).
                       
       IF ( intf_l_var_lbc(jintf) ) THEN
         do pt_lbc = 1 - lbc_halo_x, 0
           Lambda_p_in(pt_lbc) = Lambda_p_in(1) + (pt_lbc-1)*dlam_wk
           Lambda_u_in(pt_lbc) = Lambda_u_in(1) + (pt_lbc-1)*dlam_wk
         end do 
           
         do pt_lbc = 1 - lbc_halo_y, 0
           Phi_p_in(pt_lbc) = Phi_p_in(1) + (pt_lbc-1)*dphi_wk
           Phi_v_in(pt_lbc) = Phi_v_in(1) + (pt_lbc-1)*dphi_wk
         end do 
            
         do pt_lbc = LBCrow_len + 1, LBCrow_len + lbc_halo_x  
           Lambda_p_in(pt_lbc) = Lambda_p_in(lbcrow_len) +               &
                                 (pt_lbc-lbcrow_len)*dlam_wk
           Lambda_u_in(pt_lbc) = Lambda_u_in(lbcrow_len) +               &
                                 (pt_lbc-lbcrow_len)*dlam_wk
         end do
    
         do pt_lbc =  LBCrows + 1,  LBCrows + lbc_halo_y  
           Phi_p_in(pt_lbc) = Phi_p_in(lbcrows) + (pt_lbc-lbcrows)*dphi_wk
           Phi_v_in(pt_lbc) = Phi_v_in(lbcrows) + (pt_lbc-lbcrows)*dphi_wk 
         end do
       ELSE
         do pt_lbc = 1 - lbc_halo_x, 0
           Lambda_p_in(pt_lbc) = Lambda_p_in(1) + (pt_lbc-1)*dlam_wk
           Lambda_u_in(pt_lbc) = Lambda_u_in(1) + (pt_lbc-1)*dlam_wk
         end do 
           
         do pt_lbc = 1 - lbc_halo_y, 0
           Phi_p_in(pt_lbc) = Phi_p_in(1) + (pt_lbc-1)*dphi_wk
           Phi_v_in(pt_lbc) = Phi_v_in(1) + (pt_lbc-1)*dphi_wk
         end do 
 
         do pt_lbc = LBCrow_len + 1, LBCrow_len + lbc_halo_x  
           Lambda_p_in(pt_lbc) = Lambda_p_in(1) +               &
                                 (pt_lbc-1)*dlam_wk
           Lambda_u_in(pt_lbc) = Lambda_u_in(1) +               &
                                 (pt_lbc-1)*dlam_wk
         end do
    
         do pt_lbc =  LBCrows + 1,  LBCrows + lbc_halo_y  
           Phi_p_in(pt_lbc) = Phi_p_in(1) + (pt_lbc-1)*dphi_wk
           Phi_v_in(pt_lbc) = Phi_v_in(1) + (pt_lbc-1)*dphi_wk 
         end do

       END IF
! ==============================================================

! DEPENDS ON: lbc_grid_sizes
      Call LBC_Grid_Sizes (jintf, intf_l_eg_out(jintf))

! ==============================================================

       im_ident = internal_model

! Logical to indicate if model grid is rotated

!      Item and section codes of variables to be put in LBC file.
       ITEM_INTFA( 1) = lbc_stashcode_orog
       ITEM_INTFA( 2) = lbc_stashcode_u
       ITEM_INTFA( 3) = lbc_stashcode_v
       ITEM_INTFA( 4) = lbc_stashcode_w
       ITEM_INTFA( 5) = lbc_stashcode_density
       ITEM_INTFA( 6) = lbc_stashcode_theta
       ITEM_INTFA( 7) = lbc_stashcode_q
       ITEM_INTFA( 8) = lbc_stashcode_qcl
       ITEM_INTFA( 9) = lbc_stashcode_qcf
       ITEM_INTFA(10) = lbc_stashcode_exner
       ITEM_INTFA(11) = lbc_stashcode_u_adv
       ITEM_INTFA(12) = lbc_stashcode_v_adv
       ITEM_INTFA(13) = lbc_stashcode_w_adv
       Do i = 1, 13
         SECT_INTFA(i) = 32
       End Do
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qcf2
         SECT_INTFA(lbc_num) = 32
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qrain
         SECT_INTFA(lbc_num) = 32
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qgraup
         SECT_INTFA(lbc_num) = 32
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       IF (L_pc2) THEN
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_bulk
         SECT_INTFA(lbc_num) = 32
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_liquid
         SECT_INTFA(lbc_num) = 32
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_frozen
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for murk aerosol lbcs if active
       IF (L_murk) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_murk
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for dust lbcs if active
       ! Note that the presence of a dust bin in the input file (set by 
       ! l_dust_div1 etc) no longer has any bearing on whether that bin 
       ! generates those LBCs
       If (L_DUST_DIV1_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div1
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV2_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div2
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV3_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div3
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV4_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div4
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV5_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div5
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV6_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dust_div6
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for so2 lbcs if active
       IF (L_so2 .AND. L_SO2_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_so2
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for dms lbcs if active
       IF (L_dms .AND. L_DMS_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_dms
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for so4_aitken lbcs if active
       IF (L_so4_aitken .AND. L_SO4_AITKEN_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_so4_aitken
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for so4_accu lbcs if active
       IF (L_so4_accu .AND. L_SO4_ACCU_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_so4_accu
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for so4_diss lbcs if active
       IF (L_so4_diss .AND. L_SO4_DISS_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_so4_diss
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for nh3 lbcs if active
       IF (L_nh3 .AND. L_NH3_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_nh3
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for soot_new lbcs if active
       IF (L_soot_new .AND. L_SOOT_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_soot_new
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for soot_agd lbcs if active
       IF (L_soot_agd .AND. L_SOOT_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_soot_agd
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for soot_cld lbcs if active
       IF (L_soot_cld .AND. L_SOOT_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_soot_cld
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for bmass_new lbcs if active
       IF (L_bmass_new .AND. L_BMASS_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_bmass_new
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for bmass_agd lbcs if active
       IF (L_bmass_agd .AND. L_BMASS_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_bmass_agd
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for bmass_cld lbcs if active
       IF (L_bmass_cld .AND. L_BMASS_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_bmass_cld
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for ocff_new lbcs if active
       IF (L_ocff_new .AND. L_OCFF_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_ocff_new
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for ocff_agd lbcs if active
       IF (L_ocff_agd .AND. L_OCFF_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_ocff_agd
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for ocff_cld lbcs if active
       IF (L_ocff_cld .AND. L_OCFF_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_ocff_cld
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for nitr_acc lbcs if active
       If (L_nitr_acc .AND. L_NITR_ACC_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_nitr_acc
         SECT_INTFA(lbc_num) = 32
       EndIf
       ! Setup for nitr_diss lbcs if active
       If (L_nitr_diss .AND. L_NITR_DISS_LBC_OUT) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_nitr_diss
         SECT_INTFA(lbc_num) = 32
       EndIf

       ! Setup for free tracers lbcs if active
       IF (tr_vars > 0) THEN
         Do i = 1, tr_vars
           lbc_num = lbc_num + 1
           ITEM_INTFA(lbc_num) = lbc_stashcode_free_tracer
           SECT_INTFA(lbc_num) = 32
         End Do
       End If

       ! Setup for ukca tracers lbcs if active
       IF (tr_ukca > 0) THEN
         Do i = 1, tr_ukca
           lbc_num = lbc_num + 1
           ITEM_INTFA(lbc_num) = lbc_stashcode_ukca_tracer
           SECT_INTFA(lbc_num) = 32
         End Do
       End If

!      Stash codes of variables to be put in LBC file. 
!      31000+lbc_stashcode_* set to section 32xxx and need to change to 
!       section 31 stashcodes. This assumes that each field
!       has the same item number in both section 31 and section 32.

       lbc_stash_codes( 1) = sect31+lbc_stashcode_orog 
       lbc_stash_codes( 2) = sect31+lbc_stashcode_u    
       lbc_stash_codes( 3) = sect31+lbc_stashcode_v    
       lbc_stash_codes( 4) = sect31+lbc_stashcode_w    
       lbc_stash_codes( 5) = sect31+lbc_stashcode_density
       lbc_stash_codes( 6) = sect31+lbc_stashcode_theta
       lbc_stash_codes( 7) = sect31+lbc_stashcode_q    
       lbc_stash_codes( 8) = sect31+lbc_stashcode_qcl  
       lbc_stash_codes( 9) = sect31+lbc_stashcode_qcf  
       lbc_stash_codes(10) = sect31+lbc_stashcode_exner
       lbc_stash_codes(11) = sect31+lbc_stashcode_u_adv
       lbc_stash_codes(12) = sect31+lbc_stashcode_v_adv
       lbc_stash_codes(13) = sect31+lbc_stashcode_w_adv
       ! Setup for additional microphysics lbcs if active    
       lbc_num = 13                                          
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active              
         lbc_num = lbc_num + 1                               
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_qcf2
       END IF                                                 
       IF (L_mcr_qrain) THEN  ! qrain lbcs active            
         lbc_num = lbc_num + 1                               
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_qrain
       END IF                                                 
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active          
         lbc_num = lbc_num + 1                               
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_qgraup
       END IF                                                 
       ! Setup for additional cloud fraction lbcs if active  
       IF (L_pc2) THEN                                       
         lbc_num = lbc_num+1                                 
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_cf_bulk
         lbc_num = lbc_num+1                                 
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_cf_liquid
         lbc_num = lbc_num+1                                 
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_cf_frozen
       END IF                           

       ! Murk aerosol lbcs if active
       IF (L_murk) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_murk
       EndIf

       ! Setup for dust lbcs if active
       If (L_DUST_DIV1_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div1
         first_dust_lbc_number = lbc_num
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV2_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div2
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV3_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div3
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV4_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div4
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV5_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div5
       EndIf
       ! Setup for dust lbcs if active
       If (L_DUST_DIV6_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dust_div6
       EndIf
 
       ! Setup for so2 lbcs if active
       IF (L_so2 .AND. L_SO2_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_so2
       EndIf
       ! Setup for dms lbcs if active
       IF (L_dms .AND. L_DMS_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_dms
       EndIf
       ! Setup for so4_aitken lbcs if active
       IF (L_so4_aitken .AND. L_SO4_AITKEN_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_so4_aitken
       EndIf
       ! Setup for so4_accu lbcs if active
       IF (L_so4_accu .AND. L_SO4_ACCU_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_so4_accu
       EndIf
       ! Setup for so4_diss lbcs if active
       IF (L_so4_diss .AND. L_SO4_DISS_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_so4_diss
       EndIf
       ! Setup for nh3 lbcs if active
       IF (L_nh3 .AND. L_NH3_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_nh3
       EndIf
 
       ! Setup for soot_new lbcs if active
       IF (L_soot_new .AND. L_SOOT_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_soot_new
       EndIf
       ! Setup for soot_agd lbcs if active
       IF (L_soot_agd .AND. L_SOOT_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_soot_agd
       EndIf
       ! Setup for soot_cld lbcs if active
       IF (L_soot_cld .AND. L_SOOT_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_soot_cld
       EndIf
 
       ! Setup for bmass_new lbcs if active
       IF (L_bmass_new .AND. L_BMASS_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_bmass_new
       EndIf
       ! Setup for bmass_agd lbcs if active
       IF (L_bmass_agd .AND. L_BMASS_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_bmass_agd
       EndIf
       ! Setup for bmass_cld lbcs if active
       IF (L_bmass_cld .AND. L_BMASS_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_bmass_cld
       EndIf
 
       ! Setup for ocff_new lbcs if active
       IF (L_ocff_new .AND. L_OCFF_NEW_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_ocff_new
       EndIf
       ! Setup for ocff_agd lbcs if active
       IF (L_ocff_agd .AND. L_OCFF_AGD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_ocff_agd
       EndIf
       ! Setup for ocff_cld lbcs if active
       IF (L_ocff_cld .AND. L_OCFF_CLD_LBC_OUT) THEN
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_ocff_cld
       EndIf

       ! Setup for nitr_acc lbcs if active
       If (L_nitr_acc .AND. L_NITR_ACC_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_nitr_acc
       EndIf
       ! Setup for nitr_diss lbcs if active
       If (L_nitr_diss .AND. L_NITR_DISS_LBC_OUT) Then
         lbc_num = lbc_num + 1
         lbc_stash_codes(lbc_num) = sect31+lbc_stashcode_nitr_diss
       EndIf

       ! Setup for tracer lbcs if active
       ! Tracer lbcs have separate stash sections
       !  sect36 for free tracer input lbcs
       !  sect37 for UKCA tracer input lbcs
       ! This assumes that each tracer lbc stash item code corresponds 
       ! to the equivalent tracer prognostic stash item code.
       
       ! Free tracers
       IF (tr_vars > 0) THEN
         Do i = 1, tr_vars
           lbc_num = lbc_num + 1
           lbc_stash_codes(lbc_num) = sect36 + A_TR_StashItem(i)
         End Do
       End If

       ! UKCA tracers
       IF (tr_ukca > 0) THEN
         Do i = 1, tr_ukca
           lbc_num = lbc_num + 1
           lbc_stash_codes(lbc_num) = sect37 + ukca_tr_stashitem(i)
         End Do
       End If

!      Stash Codes for prognostics corresponding to above variables.
       ITEM_PROG( 1) = 33
       ITEM_PROG( 2) = 2
       ITEM_PROG( 3) = 3
       ITEM_PROG( 4) = 150
       ITEM_PROG( 5) = 253
       ITEM_PROG( 6) = 4
       ITEM_PROG( 7) = 10
       ITEM_PROG( 8) = 254
       ITEM_PROG( 9) = 12
       ITEM_PROG(10) = 255
       ITEM_PROG(11) = 256
       ITEM_PROG(12) = 257
       ITEM_PROG(13) = 258
       ! Set Section code for primary prognostics to 0 for LBC_SRC_SETUP
       Do i = 1, 13
         SECT_PROG(i) = 0
       End Do
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 271
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 272
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 273
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       IF (L_pc2) THEN
         lbc_num = lbc_num+1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 266
         lbc_num = lbc_num+1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 267
         lbc_num = lbc_num+1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 268
       EndIf
       ! Setup for murk aerosol lbcs if active
       IF (L_murk) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 90
       EndIf

       ! Setup for dust lbcs if active
       If (L_dust_div1_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 431
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div2_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 432
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div3_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 433
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div4_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 434
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div5_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 435
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div6_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 436
       EndIf
       ! Setup for so2 lbcs if active
       IF (L_so2 .AND. L_so2_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 101
       EndIf
       ! Setup for dms lbcs if active
       IF (L_dms .AND. L_dms_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 102
       EndIf
       ! Setup for so4_aitken lbcs if active
       IF (L_so4_aitken .AND. L_so4_aitken_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 103
       EndIf
       ! Setup for so4_accu lbcs if active
       IF (L_so4_accu .AND. L_so4_accu_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 104
       EndIf
       ! Setup for so4_diss lbcs if active
       IF (L_so4_diss .AND. L_so4_diss_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 105
       EndIf
       ! Setup for nh3 lbcs if active
       IF (L_nh3 .AND. L_nh3_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 107
       EndIf
       ! Setup for soot_new lbcs if active
       IF (L_soot_new .AND. L_soot_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 108
       EndIf
       ! Setup for soot_agd lbcs if active
       IF (L_soot_agd .AND. L_soot_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 109
       EndIf
       ! Setup for soot_cld lbcs if active
       IF (L_soot_cld .AND. L_soot_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 110
       EndIf
       ! Setup for bmass_new lbcs if active
       IF (L_bmass_new .AND. L_bmass_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 111
       EndIf
       ! Setup for bmass_agd lbcs if active
       IF (L_bmass_agd .AND. L_bmass_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 112
       EndIf
       ! Setup for bmass_cld lbcs if active
       IF (L_bmass_cld .AND. L_bmass_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 113
       EndIf
       ! Setup for ocff_new lbcs if active
       IF (L_ocff_new .AND. L_ocff_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 114
       EndIf
       ! Setup for ocff_agd lbcs if active
       IF (L_ocff_agd .AND. L_ocff_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 115
       EndIf
       ! Setup for ocff_cld lbcs if active
       IF (L_ocff_cld .AND. L_ocff_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 116
       EndIf

       ! Setup for nitr_acc lbcs if active
       If (L_nitr_acc .AND. L_nitr_acc_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 117
       EndIf
       ! Setup for nitr_diss lbcs if active
       If (L_nitr_diss .AND. L_nitr_diss_lbc_out) Then
         lbc_num = lbc_num + 1
         SECT_PROG(lbc_num) = 0
         ITEM_PROG(lbc_num) = 118
       EndIf

       ! Setup for free tracer arrays if active
       IF (tr_vars > 0) THEN 
         Do i = 1, tr_vars
           lbc_num = lbc_num + 1
           SECT_PROG(lbc_num) = 33
           ITEM_PROG(lbc_num) = A_TR_StashItem(i)
         End Do
       EndIf

       ! Setup for ukca tracer arrays if active
       IF (tr_ukca > 0) THEN 
         Do i = 1, tr_ukca
           lbc_num = lbc_num + 1
           SECT_PROG(lbc_num) = 34
           ITEM_PROG(lbc_num) = ukca_tr_stashitem(i)
         End Do
       EndIf

!      Addresses in D1 pointing to above prognostics.
       IPT_D1( 1) = JOROG
       IPT_D1( 2) = JU(udims%k_start)
       IPT_D1( 3) = JV(vdims%k_start)
       IPT_D1( 4) = JW(wdims_s%k_start)
       IPT_D1( 5) = JRHO(pdims%k_start)
       IPT_D1( 6) = JTHETA(tdims%k_start)
       IPT_D1( 7) = JQ(qdims%k_start)
       IPT_D1( 8) = JQCL(qdims%k_start)
       IPT_D1( 9) = JQCF(qdims%k_start)
       IPT_D1(10) = JEXNER_RHO_LEVELS(pdims%k_start)
       IPT_D1(11) = JU_ADV(udims%k_start)
       IPT_D1(12) = JV_ADV(vdims%k_start)
       IPT_D1(13) = JW_ADV(wdims_s%k_start)
       
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQCF2(qdims%k_start)
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQRAIN(qdims%k_start)
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQGRAUP(qdims%k_start)
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       IF (L_pc2) THEN
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_BULK(qdims%k_start)
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_LIQUID(qdims%k_start)
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_FROZEN(qdims%k_start)
       EndIf
       ! Setup for murk aerosol lbcs if active
       IF (L_murk) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JMURK(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div1_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV1(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div2_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV2(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div3_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV3(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div4_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV4(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div5_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV5(tdims%k_start)
       EndIf
       ! Setup for dust lbcs if active
       If (L_dust_div6_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDUST_DIV6(tdims%k_start)
       EndIf
       ! Setup for so2 lbcs if active
       IF (L_so2 .AND. L_so2_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSO2(tdims%k_start)
       EndIf
       ! Setup for dms lbcs if active
       IF (L_dms .AND. L_dms_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JDMS(tdims%k_start)
       EndIf
       ! Setup for so4_aitken lbcs if active
       IF (L_so4_aitken .AND. L_so4_aitken_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSO4_AITKEN(tdims%k_start)
       EndIf
       ! Setup for so4_accu lbcs if active
       IF (L_so4_accu .AND. L_so4_accu_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSO4_ACCU(tdims%k_start)
       EndIf
       ! Setup for so4_diss lbcs if active
       IF (L_so4_diss .AND. L_so4_diss_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSO4_DISS(tdims%k_start)
       EndIf
       ! Setup for nh3 lbcs if active
       IF (L_nh3 .AND. L_nh3_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JNH3(tdims%k_start)
       EndIf
       ! Setup for soot_new lbcs if active
       IF (L_soot_new .AND. L_soot_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSOOT_NEW(tdims%k_start)
       EndIf
       ! Setup for soot_agd lbcs if active
       IF (L_soot_agd .AND. L_soot_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSOOT_AGD(tdims%k_start)
       EndIf
       ! Setup for soot_cld lbcs if active
       IF (L_soot_cld .AND. L_soot_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JSOOT_CLD(tdims%k_start)
       EndIf
       ! Setup for bmass_new lbcs if active
       IF (L_bmass_new .AND. L_bmass_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JBMASS_NEW(tdims%k_start)
       EndIf
       ! Setup for bmass_agd lbcs if active
       IF (L_bmass_agd .AND. L_bmass_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JBMASS_AGD(tdims%k_start)
       EndIf
       ! Setup for bmass_cld lbcs if active
       IF (L_bmass_cld .AND. L_bmass_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JBMASS_CLD(tdims%k_start)
       EndIf
       ! Setup for ocff_new lbcs if active
       IF (L_ocff_new .AND. L_ocff_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JOCFF_NEW(tdims%k_start)
       EndIf
       ! Setup for ocff_agd lbcs if active
       IF (L_ocff_agd .AND. L_ocff_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JOCFF_AGD(tdims%k_start)
       EndIf
       ! Setup for ocff_cld lbcs if active
       IF (L_ocff_cld .AND. L_ocff_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JOCFF_CLD(tdims%k_start)
       EndIf

       ! Setup for nitr_acc lbcs if active
       If (L_nitr_acc .AND. L_nitr_acc_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JNITR_ACC(tdims%k_start)
       EndIf
       ! Setup for nitr_diss lbcs if active
       If (L_nitr_diss .AND. L_nitr_diss_lbc_out) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JNITR_DISS(tdims%k_start)
       EndIf

       ! 
       ! Setup for free tracer arrays if active
       ! This will need checking since we had issues with tracer dimensions in
       ! the LAM branch for Endgame.
       IF (tr_vars > 0) THEN
         Do i = 1, tr_vars
           lbc_num = lbc_num + 1
           IPT_D1(lbc_num) = JTRACER(1,i)
         End Do
       EndIf

       ! Setup for ukca tracer arrays if active
       IF (tr_ukca > 0) THEN
         Do i = 1, tr_ukca
           lbc_num = lbc_num + 1
           IPT_D1(lbc_num) = JTR_UKCA(1,i)
         End Do
       EndIf

       lbc_rim_size(rima_type_norm) = IntfWidthA(jintf)
       lbc_rim_size(rima_type_orog) = IntfWidthA(jintf)

!L 2.0 Update Information in Headers

!L     Open boundary output file if reinitialised during run

      IF (ft_steps(NFTOUT) >  0) THEN
        STRING=MODEL_FT_UNIT(NFTOUT)
        PPNAME=STRING(18:31)
        LEN_PPNAME=LEN(PPNAME)

        CALL FILE_OPEN(NFTOUT,PPNAME,LEN_PPNAME,1,1,ErrorStatus)
        IF (ErrorStatus /= 0) THEN
          Write (CMessage,*) 'Error opening preassigned boundary file'

          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        ENDIF

!      Determine position where to Buffer out data to

       NTIME=FT_LASTFIELD(NFTOUT)+1
      ELSE
       NTIME=FT_LASTFIELD(NFTOUT)+1

      ENDIF

! MakeBC only outputs orography for the first timestep
       if (ntime == 1) THEN
         Var1 = 1  ! Include orography
       else
         Var1 = 2  ! Do not include orography
       endif

!L 2.1 Fixed length header
!     FIXHD_INTFA(152,JINTF) =  INTF_LOOKUPSA * NTIME
      FIXHD_INTFA(152,JINTF) =  INTF_LOOKUPSA +                         &
     &                         (INTF_LOOKUPSA-VAR1+1)*(NTIME-1)

!L 2.2 Integer Constants
      INTHD_INTFA(3,JINTF) = NTIME

!L 2.3 LOOKUP Table
!     Determine position in LOOKUP table
      LOOKUP_START=FIXHD_INTFA(150,JINTF) +                             &
     &             FIXHD_INTFA(151,JINTF) *                             &
     &            (FIXHD_INTFA(152,JINTF) - (INTF_LOOKUPSA-VAR1+1) ) - 1
!     Check that there is enough space for this entry in LOOKUP table
      IF (FIXHD_INTFA(150,JINTF)+                                       &
     &    FIXHD_INTFA(151,JINTF)*FIXHD_INTFA(152,JINTF) >               &
     &   FIXHD_INTFA(160,JINTF)) THEN
        Write (CMessage,*)                                              &
     &  'Insufficient space for headers in boundary dataset'
        ErrorStatus = 1

        Call Ereport ( RoutineName, ErrorStatus, CMessage )
      ENDIF

! =======================================================

! DEPENDS ON: lbc_setup
      Call LBC_SetUp (                                                  &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &     item_intfa,                                                  &
     &     lbc_fld_type,                                                &
     &     lbc_halo_type,                                               &
     &     lbc_rim_type,                                                &
     &     lbc_level_type,                                              &
     &     lbc_first_level,                                             &
     &     lbc_last_level,                                              &
     &     lbc_levels,                                                  &
     &     intf_lookupsa,                                               &
     &     jintf,                                                       &
     &     im_ident                                                     &
     & )

! =============================================================
!      Hardwire level type for now
!      level type for orog lbcs = 5, so ok
!      level type for rest      = 0, but needs to be 
!      1 for rho-levels or 2 for theta-levels 
!      so we know which prognostics are on rho or theta levels
       lbc_level_type (2)  = 1
       lbc_level_type (3)  = 1
       lbc_level_type (4)  = 2
       lbc_level_type (5)  = 1
       lbc_level_type (6)  = 2
       lbc_level_type (7)  = 2
       lbc_level_type (8)  = 2
       lbc_level_type (9)  = 2
       lbc_level_type (10) = 1
       lbc_level_type (11) = 1
       lbc_level_type (12) = 1
       lbc_level_type (13) = 2
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       IF (L_pc2) THEN
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
       EndIf
       ! Setup for murk aerosol lbcs if active
       IF (L_murk) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div1 .AND. L_dust_div1_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div2 .AND. L_dust_div2_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div3 .AND. L_dust_div3_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div4 .AND. L_dust_div4_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div5 .AND. L_dust_div5_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dust lbcs if active
       IF (L_dust_div6 .AND. L_dust_div6_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for so2 lbcs if active
       IF (L_so2 .AND. L_so2_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for dms lbcs if active
       IF (L_dms .AND. L_dms_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for so4_aitken lbcs if active
       IF (L_so4_aitken .AND. L_so4_aitken_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for so4_accu lbcs if active
       IF (L_so4_accu .AND. L_so4_accu_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for so4_diss lbcs if active
       IF (L_so4_diss .AND. L_so4_diss_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for nh3 lbcs if active
       IF (L_nh3 .AND. L_nh3_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for soot_new lbcs if active
       IF (L_soot_new .AND. L_soot_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for soot_agd lbcs if active
       IF (L_soot_agd .AND. L_soot_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for soot_cld lbcs if active
       IF (L_soot_cld .AND. L_soot_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for bmass_new lbcs if active
       IF (L_bmass_new .AND. L_bmass_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for bmass_agd lbcs if active
       IF (L_bmass_agd .AND. L_bmass_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for bmass_cld lbcs if active
       IF (L_bmass_cld .AND. L_bmass_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for ocff_new lbcs if active
       IF (L_ocff_new .AND. L_ocff_new_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for ocff_agd lbcs if active
       IF (L_ocff_agd .AND. L_ocff_agd_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for ocff_cld lbcs if active
       IF (L_ocff_cld .AND. L_ocff_cld_lbc_out) THEN
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf

       ! Setup for nitr_acc lbcs if active
       If (L_nitr_acc .AND. L_nitr_acc_lbc_out) Then
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf
       ! Setup for nitr_diss lbcs if active
       If (L_nitr_diss .AND. L_nitr_diss_lbc_out) Then
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf

       ! Setup for tracer lbcs if active
       IF (tr_vars > 0) THEN
         Do i = 1, tr_vars
           lbc_num = lbc_num + 1
           lbc_level_type (lbc_num) = 2
         End Do
       EndIf

       ! Setup for UKCA tracer lbcs if active
       IF (tr_ukca > 0) THEN
         Do i = 1, tr_ukca
           lbc_num = lbc_num + 1
           lbc_level_type (lbc_num) = 2
         End Do
       EndIf

! ===============================================================

! DEPENDS ON: lbc_src_setup
      CALL LBC_Src_Setup (                                              &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &     sect_prog,                                                   &
     &     item_prog,                                                   &
     &     lbc_src_fld_type,                                            &
     &     lbc_src_halo_type,                                           &
     &     lbc_src_level_type,                                          &
     &     lbc_src_first_level,                                         &
     &     lbc_src_last_level,                                          &
     &     lbc_src_levels,                                              &
     &     intf_lookupsa,                                               &
     &     im_ident                                                     &
     & )

!  ==================================================

! DEPENDS ON: lbc_validity_time
      Call LBC_Validity_Time (LCAL360)

!  ==================================================
! DEPENDS ON: lbc_setup_lookup
      Call LBC_Setup_LookUp (                                           &
     &     lookup_intfa(1,1,jintf)                                      &
     &,    Len1_lookup                                                  &
     &,    intf_lookupsa                                                &
     &,    Item_intfa                                                   &
     &,    lbc_stash_codes                                              &
     &,    jintf                                                        &
     &,    ntime                                                        &
     &,    fixhd_intfa(1,jintf)                                         &
     &,    lbc_rim_size                                                 &
     &,    lbc_rim_type                                                 &
     &,    lbc_halo_type                                                &
     &,    lbc_fld_type                                                 &
     &,    intf_halosize                                                &
     &,    lbc_levels                                                   &
     & )
!  ==================================================


! Convert between dust bin schemes if appropriate 
! (ndustbin_in/out are both zero if dust is switched off)
       IF (ndustbin_in == 6 .AND. ndustbin_out == 2) THEN


! Ascertain length of each dust bin in the d1 array (assumes that 
!   all are the same length as dustbin 1)
         source_fld_type = lbc_src_fld_type (first_dust_lbc_number)
         source_levels    = lbc_src_levels(first_dust_lbc_number)
         source_halo_type = lbc_src_halo_type(first_dust_lbc_number)

         field_size = lasize(1,source_fld_type,source_halo_type) *      &
     &                lasize(2,source_fld_type,source_halo_type)

! j pointers should contain all the pointers into d1 for each level.  Makebc
! however only has setup the first pointer so we can only safely use that one.
! Convert 6- to 2-bin
         CALL convert_dust_six_to_two(field_size*source_levels,           &
         d1_array(jdust_div1(1):jdust_div1(1)+field_size*source_levels),  &
         d1_array(jdust_div2(1):jdust_div2(1)+field_size*source_levels),  &
         d1_array(jdust_div3(1):jdust_div3(1)+field_size*source_levels),  &
         d1_array(jdust_div4(1):jdust_div4(1)+field_size*source_levels),  &
         d1_array(jdust_div5(1):jdust_div5(1)+field_size*source_levels),  &
         d1_array(jdust_div6(1):jdust_div6(1)+field_size*source_levels) )

       ELSE IF (ndustbin_in == 2 .AND. ndustbin_out == 6) THEN

! Ascertain length of each dust bin in the d1 array (assumes that 
!   all are the same length as dustbin 1)
         source_fld_type = lbc_src_fld_type (first_dust_lbc_number)
         source_levels    = lbc_src_levels(first_dust_lbc_number)
         source_halo_type = lbc_src_halo_type(first_dust_lbc_number)

         field_size = lasize(1,source_fld_type,source_halo_type) *      &
     &                lasize(2,source_fld_type,source_halo_type)

         ALLOCATE(dust_div3_out(field_size*source_levels))
         ALLOCATE(dust_div4_out(field_size*source_levels))
         ALLOCATE(dust_div5_out(field_size*source_levels))
         ALLOCATE(dust_div6_out(field_size*source_levels))

! Convert 2- to 6-bin
         CALL convert_dust_two_to_six(field_size*source_levels,           &
         d1_array(jdust_div1(1):jdust_div1(1)+field_size*source_levels),  &
         d1_array(jdust_div2(1):jdust_div2(1)+field_size*source_levels),  &
         dust_div3_out,  &
         dust_div4_out,  &
         dust_div5_out,  &
         dust_div6_out )

       END IF


!L 3.0 Main loop to generate LBCs

       prev_src_fld_type  = -1
       prev_lbc_fld_type  = -1
       prev_src_halo_type = -1
       prev_lbc_halo_type = -1

       lbc_fld_type_interp = 1

       DO j = var1, intf_lookupsa

!      Test if weights need to be calculated

         l_calc_lbc_wts =                                               &
     &     ( lbc_src_fld_type(j)  /= prev_src_fld_type ) .or.           &
     &     ( lbc_fld_type(j)      /= prev_lbc_fld_type ) .or.           &
!    &     ( lbc_src_halo_type(j) /= prev_src_halo_type) .or.
     &     ( lbc_halo_type(j)     /= prev_lbc_halo_type)

         l_lbc_u     = lbc_fld_type(j) == fld_type_u
         l_lbc_v     = lbc_fld_type(j) == fld_type_v
         l_lbc_winds = l_lbc_u .or. l_lbc_v

         len_data      = lookup_intfa(lblrec,j,jintf)
         len_data_sect = lookup_intfa(lbnrec,j,jintf)

         lbc_size = len_data / lbc_levels(j)

         lbc_size_interp =                                              &
     &   lbc_interp_lenrima(lbc_fld_type(j),lbc_halo_type(j))
     
     ! Logic to determine if this field uses an extra surface level
         L_extra_surface_level = .FALSE.
         IF (intf_l_eg_out(jintf) .AND. lbc_first_level(j) == 0) THEN
            L_extra_surface_level =                                     &
     &        (lbc_src_first_level(j) /= lbc_first_level(j))
         END IF
         
         L_extra_surface_level_dump = .FALSE.
         IF (L_EG_dump .AND. lbc_src_first_level(j) == 0) THEN
            L_extra_surface_level_dump =                                &
              (lbc_src_first_level(j) /= lbc_first_level(j))
         END IF

         if (l_lbc_u) then

! allocate actual size for u and v in lbc_data

           allocate ( lbc_data (lbc_size_interp*lbc_levels(j),2) )
           jvar = 1
           lbc_size_u = lbc_size
         ELSE IF (l_lbc_v) then
!          space allocated with u-component
           jvar = 2
           lbc_size_v = lbc_size
         else

! need to cater for rounding to sector boundary

           len_data_max  = max ( len_data, len_data_sect )

           allocate ( lbc_data (len_data_max,1) )
           jvar = 1
         endif

         lbc_data(:,jvar) = 0.0

         source_fld_type = lbc_src_fld_type (j)
 

         global_row_length = glsize(1,source_fld_type)
         global_rows       = glsize(2,source_fld_type)

         local_row_length  = blsize(1,source_fld_type)
         local_rows        = blsize(2,source_fld_type)

         local_row_length_v = blsize(1,fld_type_v)
         local_rows_v       = blsize(2,fld_type_v)

         source_halo_type = lbc_src_halo_type(j)
         source_halo_x    = halosize(1,source_halo_type)
         source_halo_y    = halosize(2,source_halo_type)
         source_levels    = lbc_src_levels(j)

         IF (l_lbc_v) THEN
           allocate ( lbc_coeff1(lbc_size_interp) )
           allocate ( lbc_coeff2(lbc_size_interp) )
         ELSE
           allocate ( lbc_coeff1(1) )
           allocate ( lbc_coeff2(1) )
         END IF

         max_levels_per_pe = ((source_levels-1)/nproc)+1

         gather_pe = 0

         l_vi = (intf_vert_interp(jintf)  .and. j /= 1)

! If we are ENDGAME then in lbc_grid_sizes we set up whether we need extra uv
! points to allow interpolation from p.  Lets copy what we do for ND since we
! need extra P around the outside to interpolate onto u/v.
         if (intf_l_eg_out(jintf)) then
           if (l_lbc_u) then
             i_uv = 1
           else if (l_lbc_v) then
             i_uv = 2
           else
             i_uv = 0
           endif
         else
           if (l_lbc_u) then
             i_uv = 1
           else if (l_lbc_v) then
             i_uv = 2
           else
             i_uv = 0
           endif
         endif

         source_delta_lat  = a_realhd(rh_deltaNS)
         source_delta_long = a_realhd(rh_deltaEW)
         source_first_lat  = Src_Grid(Source_fld_type,1)
         source_first_long = Src_Grid(Source_fld_type,2)
         source_pole_lat   = a_realhd(rh_rotlat)
         source_pole_long  = a_realhd(rh_rotlong)

! Set up source grid  

         ALLOCATE ( lambda_source(global_row_length) )
         ALLOCATE ( phi_source(global_rows) )

! Determine if the input file is variable resolution or not
!   If variable resolution (we need a better way of determining this)...
         IF ( a_fixhd(116) > 0 .AND. a_fixhd(121) > 0 ) THEN
         
  ! Fixed header: 
  !   116 first dimension of row dependent constants
  !   117 second dimension of row dependent constants
  !   121 first dimension of column dependent constants
  !   122 second dimension of column dependent constants    
          
  !   Row dep constants: 1:global_rows = P grid points
  !                     global_rows:2*global_rows = V grid points
  !   Col dep constants: 1:global_row_length = P grid points
  !                     global_row_length:2*global_row_length = U grid points
  

           IF ( source_fld_type == fld_type_p ) THEN
            ! P grid
             WRITE(6,*) 'P grid'
             DO pt = 1, global_row_length
               lambda_source(pt) = a_coldepc(pt)
             END DO
             DO row = 1, global_rows
               phi_source(row) = a_rowdepc(row)
             END DO

           ELSE IF ( source_fld_type == fld_type_u ) THEN
            ! U grid
             WRITE(6,*) 'U grid'
             DO pt = 1, global_row_length
               lambda_source(pt) = a_coldepc(pt+a_fixhd(121))
             END DO
             DO row = 1, global_rows
               phi_source(row) = a_rowdepc(row)
             END DO
           
           ELSE IF ( source_fld_type == fld_type_v ) THEN
            ! V grid
             WRITE(6,*) 'V grid'
             DO pt = 1, global_row_length
               lambda_source(pt) = a_coldepc(pt)
             END DO
             DO row = 1, global_rows
               phi_source(row) = a_rowdepc(row+a_fixhd(116))
             END DO
           
           ELSE
             CMessage = 'Invalid source grid type'
             ErrorStatus = 9   ! Fatal Error
             CALL ereport ( RoutineName, ErrorStatus, CMessage )
                        
           END IF
         
         ELSE     ! Fixed grid

           DO pt = 1, global_row_length
             lambda_source(pt) = source_first_long + source_delta_long * &
                 (pt-1)
           END DO
           DO row = 1, global_rows
             phi_source(row)   = source_first_lat  + source_delta_lat  * &
                 (row-1)
           END DO
           
         END IF
 
         same_rotation  = ( source_pole_lat == Intf_polelat(jintf) .AND. &
                            source_pole_long == Intf_polelong(jintf) )
         IF ( .NOT. same_rotation ) THEN
! If we are not on the same rotation then we rotate source when needed.
           interp_on_p = source_rotated
         END IF

! Always rotate if we havent explicitly asked to avoid it.
         IF ( .NOT. intf_avoid_rot(jintf) ) THEN
! Always rotate source if source is rotated (even if on same rotation)
           interp_on_p   = .TRUE.
           same_rotation = .FALSE.
         END IF

! We use the P field grid if we are using P grid. 
! Otherwise we use the field grid.
         IF (interp_on_p) THEN
           lbc_first_lat  = LBC_Grid(lbc_fld_type_interp,1)
           lbc_first_long = LBC_Grid(lbc_fld_type_interp,2)
         ELSE
           lbc_first_lat  = LBC_Grid(lbc_fld_type(j),1)
           lbc_first_long = LBC_Grid(lbc_fld_type(j),2)
         END IF
         
         ! ND:
         !     The interpolation grid size is always that of
         !     the P grid, although this causes us to use
         !     extra U and V values if there is no rotation
         !     required (as then the winds are not interpolated
         !     to the P grid).
         
         ! EG:  
         !     We work with a P grid that has extra points
         !     around the outside edge of the domain. Rather 
         !     than re-defining the row length and number of 
         !     rows variables, we will use local offset 
         !     variables to access these extra points.
         
           lbc_row_len    = NINT(LBC_Grid(lbc_fld_type_interp,3))
           lbc_rows       = NINT(LBC_Grid(lbc_fld_type_interp,4))

         lbc_delta_lat  = Intf_nsspace(jintf)
         lbc_delta_long = Intf_ewspace(jintf)
         lbc_pole_lat   = Intf_polelat(jintf)
         lbc_pole_long  = Intf_polelong(jintf)
         lbc_halo_x     = Intf_halosize(1,lbc_halo_type(j))
         lbc_halo_y     = Intf_halosize(2,lbc_halo_type(j))

         IF ((l_lbc_u .OR. l_lbc_v) .AND. intf_l_eg_out(jintf)) THEN
           lbc_halo_x_coords = lbc_halo_x + 1
           lbc_halo_y_coords = lbc_halo_y + 1
         ELSE
           lbc_halo_x_coords = lbc_halo_x
           lbc_halo_y_coords = lbc_halo_y
         END IF

         rimwidth = lbc_rim_size(lbc_rim_type(j))

! Orography field
         orog_fld_type          = lbc_src_fld_type(1)
         orog_local_row_length  = blsize(1,orog_fld_type)
         orog_local_rows        = blsize(2,orog_fld_type)
         orog_global_row_length = glsize(1,orog_fld_type)
         orog_global_rows       = glsize(2,orog_fld_type)
         orog_first_lat         = src_grid(orog_fld_type,1)
         orog_first_long        = src_grid(orog_fld_type,2)

         orog_halo_type = lbc_src_halo_type(1)
         orog_halo_x    = halosize(1,orog_halo_type)
         orog_halo_y    = halosize(2,orog_halo_type)
         

         ALLOCATE (lambda_source_orog(orog_global_row_length) )
         ALLOCATE (phi_source_orog(orog_global_rows) )

! Set up orography grid
!        If variable resolution
         IF ( a_fixhd(116) > 0 .AND. a_fixhd(121) > 0 ) THEN

! Orography is done on the P grid so no need to test for u/v
         
           DO pt = 1, orog_global_row_length
             lambda_source_orog(pt) = a_coldepc(pt)
           END DO
           DO row = 1, orog_global_rows
             phi_source_orog(row) = a_rowdepc(row)
           END DO
         
         ELSE ! Fixed grid

           DO pt = 1, orog_global_row_length
            lambda_source_orog(pt) = orog_first_long + source_delta_long * &
                 (pt-1)
           END DO
           DO row = 1, orog_global_rows
             phi_source_orog(row)   = orog_first_lat  + source_delta_lat  * &
                 (row-1)
           END DO
         END IF 
         
! vi
         IF (l_vi) THEN
           max_seg_size = ((lbc_size_interp-1)/nproc) + 1
           n_segs       = ((lbc_size_interp-1)/max_seg_size) + 1
         Else
           max_seg_size = 1
           n_segs       = 1
         End If

         If (l_calc_lbc_wts) Then

           If ( allocated(lbc_index_bl)   ) deallocate (lbc_index_bl)
           If ( allocated(lbc_index_br)   ) deallocate (lbc_index_br)
           If ( allocated(lbc_weights_tr) ) deallocate (lbc_weights_tr)
           If ( allocated(lbc_weights_br) ) deallocate (lbc_weights_br)
           If ( allocated(lbc_weights_bl) ) deallocate (lbc_weights_bl)
           If ( allocated(lbc_weights_tl) ) deallocate (lbc_weights_tl)
           allocate ( lbc_index_bl  (lbc_size_interp) )
           allocate ( lbc_index_br  (lbc_size_interp) )
           allocate ( lbc_weights_tr(lbc_size_interp) )
           allocate ( lbc_weights_br(lbc_size_interp) )
           allocate ( lbc_weights_bl(lbc_size_interp) )
           allocate ( lbc_weights_tl(lbc_size_interp) )

          End If

! Point field to prognostic in d1 array to pass down easier to other routines.
! Only required if source grid is LAM as model winds need to be
! un-rotated before generating the LBCs.
! Note that copy is done for all variables ; makes it easier to
! generalise code
         If (.not. l_lbc_v .and. .not. l_lbc_u) Then
 
           field_size = lasize(1,source_fld_type,source_halo_type) *    &
     &                  lasize(2,source_fld_type,source_halo_type)


           IF (ndustbin_in == 2 .AND. ndustbin_out == 6 .AND. &
               sect_intfa(j) == 32) THEN
             SELECT CASE (item_intfa(j))
             CASE(lbc_stashcode_dust_div3)
               source_field => dust_div3_out
             CASE(lbc_stashcode_dust_div4)
               source_field => dust_div4_out
             CASE(lbc_stashcode_dust_div5)
               source_field => dust_div5_out
             CASE(lbc_stashcode_dust_div6)
               source_field => dust_div6_out
             CASE DEFAULT
               source_field => &
                 d1_array(ipt_d1(j):ipt_d1(j)+field_size*source_levels)
             END SELECT
           ELSE

!          Point prognostic to d1_array for everything not winds
             source_field => &
               d1_array(ipt_d1(j):ipt_d1(j)+field_size*source_levels)
           END IF
           ! all other fields
           DO pt_lbc = 1-lbc_halo_x_coords, LBCrow_len + lbc_halo_x_coords
             lambda_in(pt_lbc) = lambda_p_in(pt_lbc)
           END DO
           DO pt_lbc = 1-lbc_halo_y_coords, LBCrows + lbc_halo_y_coords
             phi_in(pt_lbc)    = phi_p_in(pt_lbc)
           END DO

         Else If (l_lbc_u) Then
           Allocate (source_field((local_row_length+ &
                                   2*source_halo_x)* &
                                  (local_rows+2*source_halo_y)* &
                                  source_levels))
!          Copy prognostic from d1_array into local space - for
!          winds, the u-component is copied here.

! DEPENDS ON: mbc_copy_field
           Call mbc_copy_field(d1_array(ipt_d1(j)), source_levels,      &
                               local_rows+2*source_halo_y,              &
                               local_row_length+2*source_halo_x,        &
                               source_field)
           field_size = lasize(1,source_fld_type,source_halo_type) *    &
     &                  lasize(2,source_fld_type,source_halo_type)

       
!          This pass is for a u-component so copy the v-component
!          as well. Store in source_field_v.

            Allocate (source_field_v((local_row_length_v+ &
                                      2*source_halo_x)* &
                                     (local_rows_v+2*source_halo_y)* &
                                     source_levels))

! DEPENDS ON: mbc_copy_field
           Call mbc_copy_field(d1_array(ipt_d1(j+1)), source_levels,   &
                               local_rows_v+2*source_halo_y,           &
                               local_row_length_v+2*source_halo_x,     &
                               source_field_v)

           field_size_v = lasize(1,fld_type_v,source_halo_type) *       &
     &                    lasize(2,fld_type_v,source_halo_type)
   
           If (source_rotated .AND. .NOT. same_rotation) Then
!          The model grid is a LAM grid. On the pass through for the
!          u-component, unrotate the model winds.

! lbc_unrotate_model_winds requires P-point co-ordinates as input, so sending 
! it the orography which is on the P-grid

! DEPENDS ON: lbc_unrotate_model_winds
           Call LBC_Unrotate_Model_Winds (                              &
     &          source_field                                            &
     &,         source_field_v                                          &
     &,         field_size                                              &
     &,         field_size_v                                            &
     &,         row_length                                              &
     &,         rows                                                    &
     &,         n_rows                                                  &
     &,         source_levels                                           &
     &,         source_halo_type                                        &
     &,         source_pole_lat                                         &
     &,         source_pole_long,                                       &
     &          src_grid(fld_type_p,1),                                 &
     &          src_grid(fld_type_p,2),                                 &
     &          source_delta_lat                                        &
     &,         source_delta_long                                       &
     &,         L_EG_dump )
           End If
! Use phi p grid for U grid by default.
           DO pt_lbc = 1-lbc_halo_y_coords, LBCrows + lbc_halo_y_coords
             phi_in(pt_lbc) = phi_p_in(pt_lbc)
           END DO

           If (interp_on_p) Then
! We use the lambda p grid if we are unrotating winds.
             DO pt_lbc = 1-lbc_halo_x_coords, LBCrow_len + lbc_halo_x_coords
               Lambda_in(pt_lbc) = Lambda_p_in(pt_lbc)
             END DO
           Else 
! use lambda_u grid since we are not unrotating winds and do a direct
! interpolation
             DO pt_lbc = 1-lbc_halo_x_coords, LBCrow_len + lbc_halo_x_coords
               Lambda_in(pt_lbc) = Lambda_u_in(pt_lbc)
             END DO

           End If ! source rotated

         Else If (l_lbc_v) Then

           DO pt_lbc = 1-lbc_halo_x_coords, LBCrow_len + lbc_halo_x_coords
             lambda_in(pt_lbc) = lambda_p_in(pt_lbc)
           END DO

! Use p grid for v grid if unrotated winds
           IF( interp_on_p ) THEN
             DO pt_lbc = 1-lbc_halo_y_coords, LBCrows + lbc_halo_y_coords
               phi_in(pt_lbc) = phi_p_in(pt_lbc)
             END DO
           ELSE ! Direct interpolation so we can use correct grid.
             DO pt_lbc = 1-lbc_halo_y_coords, LBCrows + lbc_halo_y_coords
               phi_in(pt_lbc) = phi_v_in(pt_lbc)
             END DO
           END IF ! rotate source
!          If not rotated grid (or chose to rotate grid) and l_lbc_u.
!          The winds are un-rotated during the pass through for the
!          u-component and the v-component is stored in source_field_v.
!      
!          Before calling make_lbcs for the v-component, copy
!          source_field_v into source_field.
 
           Allocate (source_field((local_row_length+ &
                                   2*source_halo_x)* &
                                  (local_rows+2*source_halo_y)* &
                                  source_levels))
! DEPENDS ON: mbc_copy_field
           Call mbc_copy_field(source_field_v, source_levels,        &
                               local_rows+2*source_halo_y,           &
                               local_row_length+2*source_halo_x,     &
                               source_field)
         End If ! On type of lbc (u,v or anything else)

! DEPENDS ON: make_lbcs
         CALL make_lbcs (                                               &
! Prognostic
              source_field(1)                                           &
      ,       local_row_length                                          &
      ,       local_rows                                                &
      ,       global_row_length                                         &
      ,       global_rows                                               &
      ,       source_halo_x                                             &
      ,       source_halo_y                                             &
      ,       source_levels                                             &
      ,       source_fld_type                                           &
      ,       source_halo_type                                          &
      ,       lambda_source                                             &
      ,       phi_source                                                &
      ,       source_pole_lat                                           &
      ,       source_pole_long                                          &
      ,       source_cyclic                                             &
      ,       source_rotated                                            &
! Orography Field
      ,       d1_array(ipt_d1(1))                                       &
      ,       orog_local_row_length                                     &
      ,       orog_local_rows                                           &
      ,       orog_global_row_length                                    &
      ,       orog_global_rows                                          &
      ,       orog_fld_type                                             &
      ,       orog_halo_type                                            &
      ,       orog_halo_x                                               &
      ,       orog_halo_y                                               &
      ,       lambda_source_orog                                        &
      ,       phi_source_orog                                           &
! LBC field
      ,       lbc_row_len                                               &
      ,       lbc_rows                                                  &
      ,       lbc_levels(j)                                             &
      ,       lbc_delta_lat                                             &
      ,       lbc_delta_long                                            &
      ,       lbc_first_lat                                             &
      ,       lbc_first_long                                            &
      ,       lbc_pole_lat                                              &
      ,       lbc_pole_long                                             &
      ,       lbc_halo_x                                                &
      ,       lbc_halo_y                                                &
      ,       lbc_halo_x_coords                                         &
      ,       lbc_halo_y_coords                                         &
      ,       rimwidth                                                  &
      ,       lbc_data(1,jvar)                                          &
      ,       lbc_coeff1                                                &
      ,       lbc_coeff2                                                &
      ,       lbc_size_interp                                           &
      ,       lbc_src_level_type(j)                                     &
      ,       same_rotation                                             &
      ,       intf_l_eg_out(jintf)                                      &
      ,       L_extra_surface_level_dump                                &
      ,       Lambda_in(1-lbc_halo_x_coords:lbc_row_len+lbc_halo_x_coords)&
      ,       Phi_in(1-lbc_halo_y_coords:lbc_rows+lbc_halo_y_coords)      &
! hi - indexes and weights
      ,       l_calc_lbc_wts                                            &
      ,       lbc_index_bl  , lbc_index_br                              &
      ,       lbc_weights_tr, lbc_weights_br                            &
      ,       lbc_weights_bl, lbc_weights_tl                            &
! vi
      ,       max_seg_size                                              &
      ,       n_segs                                                    &
      ,       max_levels_per_pe                                         &
      ,       gather_pe                                                 &
      ,       l_vi                                                      &
! vi - src
      ,       model_levels                                              &
      ,       a_inthd(ih_height_gen)                                    &
      ,       a_inthd(ih_1_c_rho_level)                                 &
      ,       a_realhd(rh_z_top_theta)                                  &
      ,       a_levdepc(jetatheta)                                      &
      ,       a_levdepc(jetarho)                                        &
      ,       lbc_src_first_level(j)                                    &
      ,       lbc_src_last_level(j)                                     &
! vi - lbc 
      ,       intf_v_int_order(jintf)                                   &
      ,       intf_p_levels(jintf)                                      &
      ,       inthd_intfa(ih_height_gen,jintf)                          &
      ,       lbc_first_r_rho(jintf)                                    &
      ,       lbc_z_top_model(jintf)                                    &
      ,       lbc_eta_theta(1,jintf)                                    &
      ,       lbc_eta_rho(1,jintf)                                      &
      ,       lbc_first_level(j)                                        &
      ,       lbc_last_level(j)                                         &
      ,       i_uv                                                      &
       )
! Now deallocate pointer created for winds (we point to d1 for other fields)
        If( l_lbc_u) Then
          DEALLOCATE(source_field)
        Else If( l_lbc_v) Then
          DEALLOCATE(source_field)
          DEALLOCATE(source_field_v)
          NULLIFY(source_field_v)
        End If

        NULLIFY(source_field)
  
        If (.not. l_lbc_v) Then
          deallocate (lbc_coeff1)
          deallocate (lbc_coeff2)
        Endif

        DEALLOCATE ( phi_source_orog )
        DEALLOCATE ( lambda_source_orog )
        DEALLOCATE ( phi_source )
        DEALLOCATE ( lambda_source )


        if (l_lbc_v) Then
! We might not have had to rotate source grid but we only need to rotate
! winds here if we are not on the same rotation.
          if ( .NOT. same_rotation ) THEN

!     LBCs now available for both u and v

! Input to rotate_winds : u and v are on A grid, ie, u and v on P grid

! DEPENDS ON: lbc_rotate_winds
            call lbc_rotate_winds (                                     &
     &         lbc_data(1,1)                                            &
                                !  u lbcs on p grid
     &,        lbc_data(1,2)                                            &
                                !  v lbcs on p grid
     &,        lbc_size_interp                                          &
                                !  size of u & v lbcs on A grid
     &,        lbc_levels(j)                                            &
                                !  no of lbc levels
     &,        lbc_coeff1                                               &
                                !  \ coefficients to
     &,        lbc_coeff2                                               &
                                !  / rotate winds
     & )

! Output from rotate_winds :
! u and v still on A grid. Winds are now w.r.t the rotated grid.

          END IF ! .NOT. same_rotation
 
          deallocate (lbc_coeff1)
          deallocate (lbc_coeff2)

! Interpolate u from A grid to C u-grid

          len_data_u     = lbc_size_u*lbc_levels(j)
          len_data_usect = lookup_intfa(lbnrec,j-1,jintf)

          len_data_max  = max ( len_data_u, len_data_usect )


          allocate ( lbc_data_uv(len_data_max) )

! DEPENDS ON: lbc_u_a_to_c
          call lbc_u_a_to_c (                                           &
               lbc_data(1,1)                                            &
                                !  u lbcs on p grid
      ,        lbc_data_uv                                              &
                                !  u lbcs on u grid
      ,        lbc_size_interp                                          &
                                !  field size on p grid
      ,        lbc_size_u                                               &
                                !  field size on u grid
      ,        lbc_levels(j)                                            &
                                !  no of levels
      ,        intf_row_length(jintf)                                   &
      ,        intf_p_rows(jintf)                                       &
      ,        lbc_rim_size(lbc_rim_type(j))                            &
      ,        intf_halosize(1,2)                                       &
      ,        intf_halosize(2,2)                                       &
      ,        intf_l_var_lbc(jintf)                                    &
      ,        same_rotation                                            &
      ,        intf_l_eg_out(jintf)                                     &
      ,        Lambda_p_in                                              &
      ,        Lambda_u_in                                              &
       )

! DEPENDS ON: lbc_writflds
        call lbc_writflds (nftout,lbc_data_uv,len_data_max,             &
     &                     lookup_intfa(1,j-1,jintf),                   &
     &                     fixhd_intfa(1,jintf), ltimer )

        DEALLOCATE(lbc_data_uv)

! Interpolate v from A grid to C v-grid
          len_data_v     = lbc_size_v*lbc_levels(j)
          len_data_vsect = lookup_intfa(lbnrec,j,jintf)
          len_data_max   = max ( len_data_v, len_data_vsect )

          ALLOCATE(lbc_data_uv(len_data_max))


! DEPENDS ON: lbc_v_a_to_c
          call lbc_v_a_to_c (                                           &
               lbc_data(1,2)                                            &
                                !  v lbcs on p grid
      ,        lbc_data_uv                                              &
                                !  v lbcs on v grid
      ,        lbc_size_interp                                          &
                                !  field size on p grid
      ,        lbc_size_v                                               &
                                !  field size on v grid
      ,        lbc_levels(j)                                            &
                                !  no of levels
      ,        intf_row_length(jintf)                                   &
      ,        intf_p_rows(jintf)                                       &
      ,        lbc_rim_size(lbc_rim_type(j+1))                          &
      ,        intf_halosize(1,2)                                       &
      ,        intf_halosize(2,2)                                       &
      ,        intf_l_var_lbc(jintf)                                    &
      ,        same_rotation                                            &
      ,        intf_l_eg_out(jintf)                                     &
      ,        Phi_p_in                                                 &
      ,        Phi_v_in                                                 &
       )

! DEPENDS ON: lbc_writflds
          call lbc_writflds (nftout,lbc_data_uv,len_data_max,           &
     &                       lookup_intfa(1,j,jintf),                   &
     &                       fixhd_intfa(1,jintf), ltimer )
          deallocate (lbc_data_uv)
          deallocate (lbc_data)

        END IF   !  IF (l_lbc_v)

        IF ( .NOT. l_lbc_winds ) THEN

! DEPENDS ON: lbc_post_interp_transf
          call lbc_post_interp_transf (                                 &
               lbc_size,                                                &
               lbc_first_level(j),                                      &
               lbc_last_level(j),                                       &
               sect_intfa(j),                                           &
               item_intfa(j),                                           &
               lbc_data,                                                &
               L_extra_surface_level,                                   &
               L_vi)

! DEPENDS ON: lbc_writflds
          call lbc_writflds (nftout,lbc_data,len_data_max,              &
     &                       lookup_intfa(1,j,jintf),                   &
     &                       fixhd_intfa(1,jintf), ltimer )
          deallocate (lbc_data)

        END IF !  .NOT. l_lbc_winds

        prev_src_fld_type  = lbc_src_fld_type(j)
        prev_lbc_fld_type  = lbc_fld_type(j)
        prev_src_halo_type = lbc_src_halo_type(j)
        prev_lbc_halo_type = lbc_halo_type(j)

       END DO !  j = var1, intf_lookupsa

       IF ( allocated(lbc_index_bl)   ) DEALLOCATE (lbc_index_bl)
       IF ( allocated(lbc_index_br)   ) DEALLOCATE (lbc_index_br)
       IF ( allocated(lbc_weights_tr) ) DEALLOCATE (lbc_weights_tr)
       IF ( allocated(lbc_weights_br) ) DEALLOCATE (lbc_weights_br)
       IF ( allocated(lbc_weights_bl) ) DEALLOCATE (lbc_weights_bl)
       IF ( allocated(lbc_weights_tl) ) DEALLOCATE (lbc_weights_tl)
       
       IF ( allocated(dust_div3_out) ) DEALLOCATE(dust_div3_out)
       IF ( allocated(dust_div4_out) ) DEALLOCATE(dust_div4_out)
       IF ( allocated(dust_div5_out) ) DEALLOCATE(dust_div5_out)
       IF ( allocated(dust_div6_out) ) DEALLOCATE(dust_div6_out)

        DEALLOCATE (Lambda_p_in)
        DEALLOCATE (Phi_p_in)
        DEALLOCATE (Lambda_u_in)
        DEALLOCATE (Phi_v_in)
        
!L 4.0 Write out headers/data

!L 4.1 Fixed length header


        CALL SETPOS (NFTOUT,0,ErrorStatus)

        CALL BUFFOUT(NFTOUT,FIXHD_INTFA(1:,JINTF),LEN_FIXHD,LEN_IO,A_IO)

        IF (A_IO /= -1.0 .or. LEN_IO /= LEN_FIXHD) THEN

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',LEN_FIXHD

          Write (CMessage,*) 'Error in BUFFOUT of fixed length header.'
          ErrorStatus = 2

          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF

!L 4.2 Integer constants


        CALL BUFFOUT (NFTOUT,INTHD_INTFA(1:,JINTF),                     &
     &                PP_LEN_INTHD,LEN_IO,A_IO)

        IF (A_IO /= -1.0 .or. LEN_IO /= PP_LEN_INTHD) THEN

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',PP_LEN_INTHD

          Write (CMessage,*) 'Error in BUFFOUT of integer header.'
          ErrorStatus = 3

          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF

!L 4.3 PP headers in LOOKUP table

        CALL SETPOS(NFTOUT,LOOKUP_START,ErrorStatus)

        CALL BUFFOUT(NFTOUT,LOOKUP_INTFA(1:,VAR1:,JINTF),               &
     &               LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1),LEN_IO,A_IO)

        IF(A_IO /= -1.0.OR.                                             &
     &     LEN_IO /= LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1)) THEN

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',                  &
     &    LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1)

          Write (CMessage,*) 'Error in BUFFOUT of lookup headers.'
          ErrorStatus = 4

          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF
       
!L     Close boundary output file if reinitialised during run
      IF (ft_steps(NFTOUT) >  0) THEN
        LEN_PPNAME=LEN(PPNAME)

        CALL FILE_CLOSE(NFTOUT,PPNAME,LEN_PPNAME,1,0,ErrorStatus)
      END IF

!L     Update FT_LASTFIELD
      FT_LASTFIELD(NFTOUT) = FT_LASTFIELD(NFTOUT) + 1

      IF (lhook) CALL dr_hook('GEN_INTF_A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GEN_INTF_A
