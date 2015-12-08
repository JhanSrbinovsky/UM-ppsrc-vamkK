! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:
SUBROUTINE makebc_setd1point( max_progs,len_ppindex,no_tracers,    &
           no_tr_ukca,tr_levels_in, lbcdiag_no_tracerstart,           &
           lbcdiag_no_tracerend, lbcdiag_no_trstart_ukca,          &
           pplengths,pplengths_tracer,pplengths_tr_ukca,           &
           lbcdiag_no, tracers_active, tr_ukca_active,             &
           jorog, ju, jv, jw, jrho, jtheta, jq, jqcl, jqcf,        &
           jexner, ju_adv, jv_adv, jw_adv, jqcf2, jqrain, jqgraup, &
           jmurk, jcf_bulk, jcf_liquid, jcf_frozen,                &
           jdust_div1, jdust_div2, jdust_div3, jdust_div4,         &
           jdust_div5, jdust_div6, jso2, jso4_aitken,              &
           jso4_accu, jso4_diss, jdms, jnh3, jsoot_new,            &
           jsoot_agd, jsoot_cld, jbmass_new, jbmass_agd,           &
           jbmass_cld, jocff_new, jocff_agd, jocff_cld,            &
           jnitr_acc, jnitr_diss, calc_length, d1_pointers,        &
           jtracer,jtr_ukca,halo_size_out)

      USE makebc_constants_mod

      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Decomp_DB
      ! The following is required for cintfa.h
      USE Control_Max_Sizes, ONLY : &
          max_n_intf_a,             &
          max_intf_lbcrow_length,   &
          max_intf_lbcrows,         &
          max_intf_levels
      USE cppxref_mod, ONLY:        &
          ppx_grid_type,            &
          ppx_halo_type,            & 
          ppx_lb_code,              &
          ppx_lt_code,              &
          ppx_lv_code 

      USE Submodel_Mod

      IMPLICIT NONE

!
! Description: Set D1 pointers and jpointers
!
! Method:
! Loop over the maximum number of prognostic fields
! to get the field and grid type for each field and calculate the size
! of this field. Use the size of each field to set up d1_pointers
! Finally use the value in d1_pointers to set up the jpointers
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
!
! Declarations
! Global Variables :

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
!
! Subroutine arguments
!
! Scalar arguments with INTENT(IN)
INTEGER, INTENT(IN) ::   max_progs   ! size of the lbcdiag_no array
INTEGER, INTENT(IN) ::   len_ppindex ! Dimension of pp_index
! position in lbcdiag_no array where tracers begin/end
INTEGER, INTENT(IN) ::   no_tracers
INTEGER, INTENT(IN) ::   no_tr_ukca
INTEGER, INTENT(IN) ::   tr_levels_in
INTEGER, INTENT(IN) ::   lbcdiag_no_tracerstart   
INTEGER, INTENT(IN) ::   lbcdiag_no_tracerend
INTEGER, INTENT(IN) ::   lbcdiag_no_trstart_ukca
!
! Array  arguments with INTENT(IN)
INTEGER, INTENT(IN) :: pplengths(len_ppindex,n_internal_model)   
INTEGER, INTENT(IN) :: pplengths_tracer(no_tracers,1)
INTEGER, INTENT(IN) :: pplengths_tr_ukca(no_tr_ukca,1)
INTEGER, INTENT(IN) :: lbcdiag_no(max_progs)
!
! active tracers
INTEGER, INTENT(IN) :: tracers_active(max_tracers)
INTEGER, INTENT(IN) :: tr_ukca_active(max_tr_ukca)
! Scalar arguments with INTENT(INOUT)
!   Pointers 
INTEGER, INTENT(INOUT) :: jorog      !pointer to orography
INTEGER, INTENT(INOUT) :: ju         !pointer to u wind component
INTEGER, INTENT(INOUT) :: jv         !pointer to v wind component
INTEGER, INTENT(INOUT) :: jw         !pointer to w wind component (vert.)
INTEGER, INTENT(INOUT) :: jrho       !pointer to density
INTEGER, INTENT(INOUT) :: jtheta     !pointer to theta
INTEGER, INTENT(INOUT) :: jq         !pointer to specific humidity
INTEGER, INTENT(INOUT) :: jqcl       !pointer to qcl liquid water
INTEGER, INTENT(INOUT) :: jqcf       !pointer to qcf frozen water
INTEGER, INTENT(INOUT) :: jexner     !pointer to exner pressure
INTEGER, INTENT(INOUT) :: ju_adv     !pointer to u advection
INTEGER, INTENT(INOUT) :: jv_adv     !pointer to v advection
INTEGER, INTENT(INOUT) :: jw_adv     !pointer to w advection
!   pc2 and murk pointers
INTEGER, INTENT(INOUT) :: jqcf2      !pointer to qcf2 - cloud ice (cry)
INTEGER, INTENT(INOUT) :: jqrain     !pointer to qrain - rain
INTEGER, INTENT(INOUT) :: jqgraup    !pointer to qgraup - graupel horiz.
INTEGER, INTENT(INOUT) :: jmurk      !pointer to murk
INTEGER, INTENT(INOUT) :: jcf_bulk   !pointer to cf_bulk
INTEGER, INTENT(INOUT) :: jcf_liquid !pointer to cf_liquid
INTEGER, INTENT(INOUT) :: jcf_frozen !pointer to cf_frozen
!   Aerosol pointers 
INTEGER, INTENT(INOUT) :: jdust_div1      ! dust mmr, division 1
INTEGER, INTENT(INOUT) :: jdust_div2      ! dust mmr, division 2
INTEGER, INTENT(INOUT) :: jdust_div3      ! dust mmr, division 3
INTEGER, INTENT(INOUT) :: jdust_div4      ! dust mmr, division 4
INTEGER, INTENT(INOUT) :: jdust_div5      ! dust mmr, division 5
INTEGER, INTENT(INOUT) :: jdust_div6      ! dust mmr, division 6
INTEGER, INTENT(INOUT) :: jso2            ! sulphur dioxide gas
INTEGER, INTENT(INOUT) :: jso4_aitken     ! Aitken mode sulphate aer
INTEGER, INTENT(INOUT) :: jso4_accu       ! accumulation mode sulpha
INTEGER, INTENT(INOUT) :: jso4_diss       ! dissloved  sulphate aero
INTEGER, INTENT(INOUT) :: jdms            ! dimethyl sulphide gas
INTEGER, INTENT(INOUT) :: jnh3            ! ammonia gas mmr
INTEGER, INTENT(INOUT) :: jsoot_new       ! fresh soot mmr
INTEGER, INTENT(INOUT) :: jsoot_agd       ! aged soot mmr
INTEGER, INTENT(INOUT) :: jsoot_cld       ! soot in cloud mmr
INTEGER, INTENT(INOUT) :: jbmass_new      ! fresh biomass mmr
INTEGER, INTENT(INOUT) :: jbmass_agd      ! aged biomass mmr
INTEGER, INTENT(INOUT) :: jbmass_cld      ! cloud biomass mmr
INTEGER, INTENT(INOUT) :: jocff_new       ! fresh OCFF mmr
INTEGER, INTENT(INOUT) :: jocff_agd       ! aged OCFF mmr
INTEGER, INTENT(INOUT) :: jocff_cld       ! OCFF in cloud mmr
INTEGER, INTENT(INOUT) :: jnitr_acc       ! Accumulation nitrate aerosol
INTEGER, INTENT(INOUT) :: jnitr_diss      ! Dissolved nitrate aerosol
!
!   calc_length is the length required to store all prognostics
INTEGER, INTENT(INOUT) ::  calc_length
! d1_pointers - the position of each diagnostic in the fields file.
!
INTEGER, INTENT(INOUT) :: d1_pointers(max_progs)
! Pointer to free and UKCA tracer(s)
!
INTEGER, INTENT(INOUT) :: jtracer(tr_levels_in,no_tracers+1)  
INTEGER, INTENT(INOUT) :: jtr_ukca(tr_levels_in,no_tr_ukca+1)  

! Array  arguments with INTENT(OUT)
! halo_size_out - halo size of all prognostics
INTEGER, INTENT(OUT) :: halo_size_out(max_progs,3)

! Local variables
!   Local scalars :

INTEGER :: loop, k
INTEGER :: tracer_no, tracer_no_ukca
INTEGER :: cumil_length_d1
INTEGER :: field_size_x
INTEGER :: field_size_y
INTEGER :: field_size_z
INTEGER :: field_size
INTEGER :: iproc=0

INTEGER :: lb_code
INTEGER :: lt_code
INTEGER :: lv_code
INTEGER :: bottom_level
INTEGER :: top_level
INTEGER :: grid_code
INTEGER :: fld_type
INTEGER :: halo_type

INTEGER :: section
INTEGER :: item
!
! Function definitions
! Declare functions exppxi and get_fld_type
INTEGER :: exppxi
INTEGER :: get_fld_type
INTEGER :: errorstatus
CHARACTER(LEN=80) ::  cmessage   !error message
!
! Local parameters 
CHARACTER (LEN=80), PARAMETER :: RoutineName='makebc_setd1point'
!
! End of header
! *********************************************************************
! 
! Initialise variables prior to entering the loop
cumil_length_d1=0

tracer_no=0
tracer_no_ukca=0

! Loop over the maximum number of prognostic fields
! (all section zero fields whether on or not, plus
! the number of free and UKCA tracer fields requested)
!
DO loop=1,max_progs

! Check that this LBC is required before working on it
  IF (lbcdiag_no(loop) /= -1) THEN

!
! find correct section and item number
! Easy for section zero items - section is zero and
! item is stored in lbcdiag_no(loop)
! More complicated for Free tracers and UKCA tracers
!
! Free tracers
    IF (loop >= lbcdiag_no_tracerstart .AND.                     &
         loop <= lbcdiag_no_tracerend) THEN
        ! Zero item if you've just started doing the free tracers
        ! item is the last item number found to be 
        ! present in tracers_active
        IF (loop == lbcdiag_no_tracerstart) item=0
        !
        ! Free tracers are in section 33
        section=sec_freetr
        !
        ! Increment tracer_no
        tracer_no=tracer_no+1
        !
        ! loop from item+1 to max_tracers
        ! until we find which is the next active tracer
        !
        DO k=item+1,max_tracers
          !
          ! If have requested this tracer
          ! set a_tr_stashitem and item to k 
          ! and exit loop 
          ! set item so that next search for an 
          ! active tracer starts from next possible
          ! active tracer
          !
          IF (tracers_active(k)==1) THEN
            
            a_tr_stashitem(tracer_no)=k
            
            item=k
            EXIT
          END IF
        END DO
      ! UKCA tracers
      ELSE IF (loop >=  lbcdiag_no_trstart_ukca) THEN
        ! Zero item if you've just started doing the UKCA tracers
        IF (loop == lbcdiag_no_trstart_ukca) item=0
        !
        ! UKCA tracers are in section 34
        section=sec_ukcatr
        !
        ! Increment tracer_no_ukca
        tracer_no_ukca=tracer_no_ukca+1
        
        ! loop until we find which is the next active tracer
        DO k=item+1,max_tr_ukca
          !
          ! If have requested this tracer
          ! set a_trukca_stashitem and item to k 
          ! and exit loop over tracer_no_ukca
          IF (tr_ukca_active(k)==1) THEN
            !
            ukca_tr_stashitem(tracer_no_ukca)=k
            item=k
            EXIT
          END IF
        END DO
!  
! All other items are in section 0 and item number is in lbcsdiag_no
      ELSE
        section=0
        item=lbcdiag_no(loop)
      END IF
!
! Having got the section and item codes now work out the size of this field
! Need to calculate all dimensions using UM functions
!
      ! Extract bottom level code from STASHmaster record.
      lb_code=exppxi(1,section,item,ppx_lb_code,               &
                        errorstatus,cmessage)
      IF (errorstatus /= 0) THEN
        CALL ereport(routinename,errorstatus,cmessage)
      END IF

      ! Extract top level code from STASHmaster record.
      lt_code=exppxi(1,section,item,ppx_lt_code,               &
                        errorstatus,cmessage)
      IF (errorstatus /= 0) THEN
        CALL ereport(routinename,errorstatus,cmessage)
      END IF
 
      ! Extract level type code from STASHmaster record.
      lv_code=exppxi(1,section,item,ppx_lv_code,               &
                        errorstatus,cmessage)
      IF (errorstatus /= 0) THEN
        CALL ereport(routinename,errorstatus,cmessage)
      END IF

      ! Similar code can be found in lbc_src_setup.
      ! If level type is not 5 (i.e. not single level) then convert bottom and
      ! top level codes to values we can use to compare.
      IF (lv_code /= 5) THEN
! DEPENDS ON: levcod
        CALL levcod(lb_code, bottom_level,errorstatus,cmessage)
! DEPENDS ON: levcod
        CALL levcod(lt_code, top_level,   errorstatus,cmessage)
      ELSE
        ! If single level set the top and bottom levels both to 1.
        bottom_level = 1
        top_level    = 1
      END IF

! Use exppxi to get the field/grid type from the ppx array
! DEPENDS ON: exppxi
      grid_code=exppxi(1,section,item,ppx_grid_type,               &
                        errorstatus,cmessage)
      IF (errorstatus /= 0) THEN

        CALL ereport(routinename,errorstatus,cmessage)
      END IF
! DEPENDS ON: exppxi
      halo_type=exppxi(1,section,item,ppx_halo_type,               &
                        errorstatus,cmessage)
      IF (errorstatus /= 0) THEN

        CALL ereport(routinename,errorstatus,cmessage)
      END IF
!
      halo_size_out(loop,1)= &
          decompDB(decomp_smexe)%halosize(1,halo_type)
      halo_size_out(loop,2)= &
          decompDB(decomp_smexe)%halosize(1,halo_type)
      halo_size_out(loop,3)= &
          decompDB(decomp_smexe)%halosize(1,halo_type)

! use get_fld_type to work out which field type it is
! DEPENDS ON: get_fld_type
      fld_type=get_fld_type(grid_code)

      field_size_x=decompDB(decomp_smexe)%g_blsize(1,fld_type,iproc) &
          +(2*decompDB(decomp_smexe)%halosize(1,halo_type))
      field_size_y=decompDB(decomp_smexe)%g_blsize(2,fld_type,iproc) &
          +(2*decompDB(decomp_smexe)%halosize(2,halo_type))
!
      IF (loop >= lbcdiag_no_tracerstart .AND.                       &
          loop <= lbcdiag_no_tracerend) THEN
            field_size_z=pplengths_tracer(tracer_no,1)
      ELSE IF (loop >= lbcdiag_no_trstart_ukca) THEN
            field_size_z=pplengths_tr_ukca(tracer_no_ukca,1)
      ELSE
        field_size_z=pplengths(item,1)
      END IF

      IF (field_size_z /= top_level - bottom_level + 1) THEN
        WRITE (cmessage,'(A,I4,I4)')                                 &
          "Incorrect number of levels available for ", section, item
        errorstatus = 30
        CALL ereport ( RoutineName ,ErrorStatus, CMessage)
      END IF
!
! Calculate size of the field as product of x,y and z 
! dimensions
     field_size=field_size_x*field_size_y*field_size_z

! d1 pointer for this field is the cumulative
! length of all previous fields plus 1
     d1_pointers(loop)=cumil_length_d1+1
!
! Set tracer jpointers here - other pointers set outside the loop
!
     IF (loop >= lbcdiag_no_tracerstart .AND.                     &
         loop <= lbcdiag_no_tracerend) THEN
          jtracer(1,tracer_no)=d1_pointers(loop)
     END IF
     IF (loop >= lbcdiag_no_trstart_ukca) THEN
       jtr_ukca(1,tracer_no_ukca)=d1_pointers(loop)
     END IF
!
! Increase cumil_length_d1 by the size of this field
     cumil_length_d1=cumil_length_d1 + field_size

! Work out the length required to store each prognostic and add
! up to create a cumulative total
     calc_length=calc_length+field_size

   END IF ! Prognostic present in dump
   
 END DO   ! Loop over the maximum number of prognostic fields
! *********************************************************************
! set jpointers from d1_pointers created previously
! use index from module makebc_constants_mod
!
ju         = d1_pointers(lbcindex_u)
jv         = d1_pointers(lbcindex_v)
jtheta     = d1_pointers(lbcindex_theta)
jq         = d1_pointers(lbcindex_q)
jqcf       = d1_pointers(lbcindex_qcf)
jorog      = d1_pointers(lbcindex_orog)
jmurk      = d1_pointers(lbcindex_murk)
jw         = d1_pointers(lbcindex_w)
jrho       = d1_pointers(lbcindex_rho)
jqcl       = d1_pointers(lbcindex_qcl)
jexner     = d1_pointers(lbcindex_exner)
ju_adv     = d1_pointers(lbcindex_u_adv)
jv_adv     = d1_pointers(lbcindex_v_adv)
jw_adv     = d1_pointers(lbcindex_w_adv)
jcf_bulk   = d1_pointers(lbcindex_cf_bulk)
jcf_liquid = d1_pointers(lbcindex_cf_liquid)
jcf_frozen = d1_pointers(lbcindex_cf_frozen)
jqcf2      = d1_pointers(lbcindex_qcf2)
jqrain     = d1_pointers(lbcindex_qrain)
jqgraup    = d1_pointers(lbcindex_qgraup)
jdust_div1 = d1_pointers(lbcindex_dust_div1)
jdust_div2 = d1_pointers(lbcindex_dust_div2)
jdust_div3 = d1_pointers(lbcindex_dust_div3)
jdust_div4 = d1_pointers(lbcindex_dust_div4)
jdust_div5 = d1_pointers(lbcindex_dust_div5)
jdust_div6 = d1_pointers(lbcindex_dust_div6)
jso2       = d1_pointers(lbcindex_so2)
jso4_aitken= d1_pointers(lbcindex_so4_aitken)
jso4_accu  = d1_pointers(lbcindex_so4_accu)
jso4_diss  = d1_pointers(lbcindex_so4_diss)
jdms       = d1_pointers(lbcindex_dms)
jnh3       = d1_pointers(lbcindex_nh3)
jsoot_new  = d1_pointers(lbcindex_soot_new)
jsoot_agd  = d1_pointers(lbcindex_soot_agd)
jsoot_cld  = d1_pointers(lbcindex_soot_cld)
jbmass_new = d1_pointers(lbcindex_bmass_new)
jbmass_agd = d1_pointers(lbcindex_bmass_agd)
jbmass_cld = d1_pointers(lbcindex_bmass_cld)
jocff_new  = d1_pointers(lbcindex_ocff_new)
jocff_agd  = d1_pointers(lbcindex_ocff_agd)
jocff_cld  = d1_pointers(lbcindex_ocff_cld)
jnitr_acc  = d1_pointers(lbcindex_nitr_acc)
jnitr_diss = d1_pointers(lbcindex_nitr_diss)
            

END SUBROUTINE makebc_setd1point

