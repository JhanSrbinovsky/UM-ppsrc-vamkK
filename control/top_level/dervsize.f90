! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Program: DERVSIZE -------------------------------------------------
!LL
!LL  Purpose: Calculate extra sizes required for dynamic allocation of
!LL           main memory in the model, derived from sizes passed by
!LL           READSIZE into the top level program UM_SHELL. These are
!LL           local sizes for each pe, except where explicitly stated,
!LL           having previously called decomposition routines.
!LL
!LL  Programming standard: UM Doc Paper 3
!LL
!LL  External documentation: On-line UM document C1 - Dynamic allocation
!LL                          of primary fields
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE DERVSIZE(                                              &
                   ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE Decomp_DB
      USE domain_params
      USE lbc_mod
      USE dust_parameters_mod, ONLY:               l_dust,              &
           l_dust_div1,       l_dust_div2,         l_dust_div3,         &
           l_dust_div4,       l_dust_div5,         l_dust_div6
      USE um_input_control_mod, ONLY:                                   &
           model_domain,                                                &
           l_dust_div1_lbc,   l_dust_div1_lbc_out,                      &
           l_dust_div2_lbc,   l_dust_div2_lbc_out,                      &
           l_dust_div3_lbc,   l_dust_div3_lbc_out,                      &
           l_dust_div4_lbc,   l_dust_div4_lbc_out,                      &
           l_dust_div5_lbc,   l_dust_div5_lbc_out,                      &
           l_dust_div6_lbc,   l_dust_div6_lbc_out, l_so2,               &
           l_so2_lbc,         l_so2_lbc_out,       l_dms,               &
           l_dms_lbc,         l_dms_lbc_out,       l_so4_aitken,        &
           l_so4_aitken_lbc,  l_so4_accu,          l_so4_aitken_lbc_out,&
           l_so4_accu_lbc,    l_so4_accu_lbc_out,  l_so4_diss,          &
           l_so4_diss_lbc,    l_so4_diss_lbc_out,  l_nh3,               &
           l_nh3_lbc,         l_nh3_lbc_out,       l_soot_new,          &
           l_soot_new_lbc,    l_soot_new_lbc_out,  l_soot_agd,          &
           l_soot_agd_lbc,    l_soot_agd_lbc_out,  l_soot_cld,          &
           l_soot_cld_lbc,    l_soot_cld_lbc_out,  l_bmass_new,         &
           l_bmass_new_lbc,   l_bmass_new_lbc_out, l_bmass_agd,         &
           l_bmass_agd_lbc,   l_bmass_agd_lbc_out, l_bmass_cld,         &
           l_bmass_cld_lbc,   l_bmass_cld_lbc_out, l_ocff_new,          &
           l_ocff_new_lbc,    l_ocff_new_lbc_out,  l_ocff_agd,          &
           l_ocff_agd_lbc,    l_ocff_agd_lbc_out,  l_ocff_cld,          &
           l_ocff_cld_lbc,    l_ocff_cld_lbc_out,  l_nitr_acc,          &
           l_nitr_acc_lbc,    l_nitr_acc_lbc_out,  l_nitr_diss,         &
           l_nitr_diss_lbc,   l_nitr_diss_lbc_out, l_soot,              &
           l_ocff 

      USE mphys_inputs_mod, ONLY:                                        &
           l_mcr_qcf2,        l_mcr_qrain,         l_mcr_qgraup,        &
           l_mcr_qcf2_lbc,    l_mcr_qrain_lbc,     l_mcr_qgraup_lbc 

      USE cloud_inputs_mod, ONLY: l_pc2, l_pc2_lbc
      USE murk_inputs_mod,  ONLY: l_murk, l_murk_lbc
      USE cv_run_mod,       ONLY: l_3d_cca, l_ccrad

      IMPLICIT NONE
!
!  Argument list and comdecks
!
      INTEGER ICODE             ! OUT - Return code
      CHARACTER(LEN=80) CMESSAGE     ! OUT - Error message

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
! For STASH sizes 
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
!  Local variables
      INTEGER                                                           &
        iside                                                           &
                     ! loop counter for LBC sides
      , ifld                                                            &
                     ! loop counter for field types
      , ihalo                                                           &
                     ! loop counter for halo types
      , iproc                                                           &
                     ! loop counter for processors
      , irim                                                            &
                     ! loop counter for rim types
      , info                                                            &
                     ! return code from GCOM
      , lbc_row_len                                                     &
                     ! length of row of LBC
      , lbc_nrows                                                       &
                     ! number of rows in LBC
      , num_optional_lbcs_in                                            &
                              ! no. of optional lbc fields in input
      , num_optional_lbcs_out ! no. of optional lbc fields output

!*----------------------------------------------------------------------

      INTEGER nohalo_IMT,nohalo_JMT,glob_IMT,glob_JMT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('DERVSIZE',zhook_in,zhook_handle)
      ICODE=0

!     Initialise number of optional in/out lbcs to zero
      num_optional_lbcs_in  = 0
      num_optional_lbcs_out = 0
!L
!L   Atmosphere Boundary Datasets.
!L   2nd dimension of Level Dependent Constants.
      INTF_LEN2_LEVDEPC=4
!L   2nd dimension of Row/Col Dependent Constants.
      INTF_LEN2_ROWDEPC=2
      INTF_LEN2_COLDEPC=2      
!L
! Section only needed in model.

!L   Sizes applicable to all resolutions
      THETA_FIELD_SIZE=ROW_LENGTH*ROWS
! One less V row on the Northern most processors
      IF (at_extremity(PNorth) .and.                                    &
               model_domain  /=  mt_bi_cyclic_lam) THEN
        IF (l_vatpoles) THEN
          N_ROWS=ROWS+1
        ELSE
          N_ROWS=ROWS-1
        END IF ! vatpoles
      ELSE
        N_ROWS=ROWS
      ENDIF
      U_FIELD_SIZE=THETA_FIELD_SIZE
      V_FIELD_SIZE=ROW_LENGTH*N_ROWS
      theta_off_size   = (row_length + 2*offx)   * (rows   + 2*offy)
      theta_halo_size  = (row_length + 2*halo_i) * (rows   + 2*halo_j)
      u_off_size       = (row_length + 2*offx)   * (rows   + 2*offy)
      u_halo_size      = (row_length + 2*halo_i) * (rows   + 2*halo_j)
      v_off_size       = (row_length + 2*offx)   * (n_rows + 2*offy)
      v_halo_size      = (row_length + 2*halo_i) * (n_rows + 2*halo_j)

!     Grid bounds settings
      CALL atm_fields_bounds_init(offx,offy,halo_i,halo_j,              &
                  row_length,rows,n_rows,model_levels,wet_levels,       &
                  tr_levels,bl_levels,ozone_levels)

!     No of levels for Convective Cloud Amount.
      IF (L_3D_CCA .OR. L_CCRAD) THEN
        ! This needs to be the number of wet levels from level 1 - even for EG.
        N_CCA_LEV = qdims%k_end
      ELSE
        N_CCA_LEV = 1
      ENDIF
      IF(PrintStatus >= PrStatus_Normal) THEN
        WRITE(6,*)                                                      &
        'DERVSIZE: Number of levels for convective clouds is ',         &
        N_CCA_LEV
      ENDIF

      ! ----------------------------------------------------------------
      ! Count number of optional lateral boundary fields expected in
      ! input dependent on whether the _lbc logicals are true or false
      ! ----------------------------------------------------------------

      ! Additional microphysics variables (ice crystals, rain, graupel)
      If (L_mcr_qcf2_lbc) Then  ! qcf2 lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If
      If (L_mcr_qrain_lbc) Then  ! qrain lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If
      If (L_mcr_qgraup_lbc) Then  ! qgraup lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If

      ! Cloud fractions for PC2 (3 fields: bulk, liquid and frozen)
      If (L_pc2_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 3
      EndIf

      ! Murk aerosol
      If (L_murk_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! Dust
      If (L_dust_div1_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      If (L_dust_div2_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      If (L_dust_div3_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      If (L_dust_div4_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      If (L_dust_div5_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      If (L_dust_div6_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! SO2
      If (L_so2_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! DMS
      If (L_dms_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! SO4_Aitken
      If (L_so4_aitken_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! SO4_ACCU
      If (L_so4_accu_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! SO4_DISS
      If (L_so4_diss_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! NH3
      If (L_NH3_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! SOOT_NEW
      If (L_soot_new_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! SOOT_AGD
      If (L_soot_agd_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! SOOT_CLD
      If (L_soot_cld_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! BMASS_NEW
      If (L_bmass_new_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! BMASS_AGD
      If (L_bmass_agd_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! BMASS_CLD
      If (L_bmass_cld_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! OCFF_NEW
      If (L_ocff_new_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! OCFF_AGD
      If (L_ocff_agd_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! OCFF_CLD
      If (L_ocff_cld_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! NITR_ACC
      If (L_nitr_acc_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf
      ! NITR_DISS
      If (L_nitr_diss_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! ----------------------------------------------------------------
      ! Count number of optional lateral boundary fields to write out
      ! dependent on whether the prognostics are active or not
      ! ----------------------------------------------------------------

      ! Additional microphysics variables (ice crystals, rain, graupel)
      If (L_mcr_qcf2) Then  ! qcf2 lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If
      If (L_mcr_qrain) Then  ! qrain lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If
      If (L_mcr_qgraup) Then  ! qgraup lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If

      ! Cloud fractions for PC2 (3 fields: bulk, liquid and frozen)
      If (L_pc2) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 3
      EndIf

      ! Murk aerosol
      If (L_murk) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! Dust
      If (L_dust_div1 .AND. L_dust_div1_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      If (L_dust_div2 .AND. L_dust_div2_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      If (L_dust_div3 .AND. L_dust_div3_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      If (L_dust_div4 .AND. L_dust_div4_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      If (L_dust_div5 .AND. L_dust_div5_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      If (L_dust_div6 .AND. L_dust_div6_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! SO2
      If (L_so2 .AND. L_so2_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! DMS
      If (L_dms .AND. L_dms_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! SO4_Aitken
      If (L_so4_aitken .AND. L_so4_aitken_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! SO4_ACCU
      If (L_so4_accu .AND. L_so4_accu_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out+ 1
      EndIf
      ! SO4_DISS
      If (L_so4_diss .AND. L_so4_diss_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! NH3
      If (L_NH3 .AND. L_NH3_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! SOOT_NEW
      If (L_soot_new .AND. L_soot_new_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! SOOT_AGD
      If (L_soot_agd .AND. L_soot_agd_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! SOOT_CLD
      If (L_soot_cld .AND. L_soot_cld_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! BMASS_NEW
      If (L_bmass_new .AND. L_bmass_new_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! BMASS_AGD
      If (L_bmass_agd .AND. L_bmass_agd_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! BMASS_CLD
      If (L_bmass_cld .AND. L_bmass_cld_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! OCFF_NEW
      If (L_ocff_new .AND. L_ocff_new_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! OCFF_AGD
      If (L_ocff_agd .AND. L_ocff_agd_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! OCFF_CLD
      If (L_ocff_cld .AND. L_ocff_cld_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! NITR_ACC
      If (L_nitr_acc .AND. L_nitr_acc_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf
      ! NITR_DISS
      If (L_nitr_diss .AND. L_nitr_diss_lbc_out) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! ----------------------------------------------------------------


      A_LEN1_LEVDEPC=MODEL_LEVELS+1
! We use the global values here
      IF (l_vatpoles) THEN
        A_LEN1_ROWDEPC= global_rows + 1
      ELSE
        A_LEN1_ROWDEPC= global_rows
      END IF
      A_LEN1_COLDEPC= global_row_length
      A_LEN1_FLDDEPC= A_LEN1_ROWDEPC * A_LEN1_COLDEPC

!L Number of atmosphere model interface lookups
      intf_lookupsa = 13 + num_optional_lbcs_out



!
! Calculate sizes of Lateral Boundary Conditions Data:
! LENRIMA(fld_type,halo_type,rim_type)
!    size of a single level of LBC data for a given field type and
!    halo type and rim_type
! LBC_SIZEA(side,fld_type,halo_type,rim_type)
!    size of a single edge of LBC data for a given edge (North,
!    East,South,West) and a given field type and halo type and rim type
! LBC_STARTA(side,fld_type,halo_type,rim_type)
!    start address in LBC array of a given edge for a given
!    field type and halo type and rim_type
!

      IF (.NOT.ASSOCIATED(g_lbc_starta)) &
          ALLOCATE(g_lbc_starta          &
          (4,nfld_max,nhalo_max,nrima_max,0:nproc_max-1))

      DO irim=1,Nrima_max
      DO ifld=1,Nfld_max        ! loop over field types
        DO ihalo=1,NHalo_max    ! loop over halo types

          LENRIMA(ifld,ihalo,irim)=0
          global_LENRIMA(ifld,ihalo,irim)=0

          DO iside=1,4          ! loop over North,East,South,West

            DO iproc=0,nproc-1
              g_LBC_STARTA(iside,ifld,ihalo,irim,iproc)=0
            ENDDO
! First calculate the global_LENRIMA values - these are the sizes
! of the complete LBC as stored on disk before it is decomposed
! over processors.

            IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
! North or South boundaries

              IF ( l_vatpoles ) THEN
              lbc_row_len=glsize(1,ifld)
              ELSE
              IF (ifld  ==  fld_type_u) THEN
                lbc_row_len=glsize(1,ifld)-1
              ELSE
                lbc_row_len=glsize(1,ifld)
              END IF
              END IF  ! vatpoles

              lbc_row_len=lbc_row_len+2*halosize(1,ihalo)

              lbc_nrows=halosize(2,ihalo)+RIMWIDTHA(irim)


            ELSE
! East or West boundaries
              lbc_row_len=halosize(1,ihalo)+RIMWIDTHA(irim)

              lbc_nrows=glsize(2,ifld)-2*RIMWIDTHA(irim)

            ENDIF ! North/South or East/West boundary

            IF (RIMWIDTHA(irim)  >   0) THEN

              global_LBC_SIZEA(iside,ifld,ihalo,irim)=                  &
                lbc_row_len*lbc_nrows

              IF (iside  ==  1) THEN
                global_LBC_STARTA(iside,ifld,ihalo,irim)=1
              ELSE
                global_LBC_STARTA(iside,ifld,ihalo,irim)=               &
                  global_LBC_STARTA(iside-1,ifld,ihalo,irim)+           &
                  global_LBC_SIZEA(iside-1,ifld,ihalo,irim)
              ENDIF

              global_LENRIMA(ifld,ihalo,irim)=                          &
                global_LENRIMA(ifld,ihalo,irim)+                        &
                global_LBC_SIZEA(iside,ifld,ihalo,irim)

            ELSE ! No LBCs if RIMWIDTH is  <=  0)

               global_LBC_SIZEA(iside,ifld,ihalo,irim)=0
               global_LBC_STARTA(iside,ifld,ihalo,irim)=1

            ENDIF ! IF (RIMWIDTHA  >   0)

! Now calculate local LENRIMA values, and the associated arrays
! LBC_SIZEA and LBC_STARTA

            IF (at_extremity(iside) .AND.                               &
                (RIMWIDTHA(irim)  >   0)) THEN
! This processor is at the edge of the grid

              IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
! North or South boundaries
! North/South boundaries can contain the corners of the LBCs

! Length of rows includes the halos.
                IF (.NOT. l_vatpoles) THEN
! New Dynamics only: For U fields, there is one less point as the last 
! point along each row is ignored. This only applies to the last processor 
! along the row
                IF ((ifld  ==  fld_type_u) .AND.                        &
                    (at_extremity(PEast))) THEN
                  lbc_row_len=lasize(1,ifld,ihalo)-1
                ELSE
                  lbc_row_len=lasize(1,ifld,ihalo)
                END IF
                ELSE
                  lbc_row_len=lasize(1,ifld,ihalo)
                END IF   ! vatpoles

! And the number of rows is the size of the halo plus the rimwidth
                lbc_nrows=halosize(2,ihalo)+RIMWIDTHA(irim)

              ELSE
! East or West boundaries

! Length of row is the size of the halo plus the rimwidth
                lbc_row_len=halosize(1,ihalo)+RIMWIDTHA(irim)

! Number of rows is the number of rows on this PE minus any
! rows which are looked after by the North/South boundaries
! (ie. the corners).
                lbc_nrows=lasize(2,ifld,ihalo)
                IF (at_extremity(PNorth))                               &
                  lbc_nrows=lbc_nrows-halosize(2,ihalo)-RIMWIDTHA(irim)
                IF (at_extremity(PSouth))                               &
                  lbc_nrows=lbc_nrows-halosize(2,ihalo)-RIMWIDTHA(irim)

              ENDIF ! North/South or East/West boundary

              LBC_SIZEA(iside,ifld,ihalo,irim)=lbc_row_len*lbc_nrows

            ELSE
! This processor is not at the edge of the grid, or RIMWIDTH is
! zero (indicating no LBCs)
              LBC_SIZEA(iside,ifld,ihalo,irim)=0

            ENDIF

! LBC_STARTA contains the offset in the LBC array for each side
! (North,East,South,West) piece of data

            IF (iside  ==  1) THEN
              LBC_STARTA(iside,ifld,ihalo,irim)=1
            ELSE
              LBC_STARTA(iside,ifld,ihalo,irim)=                        &
                LBC_STARTA(iside-1,ifld,ihalo,irim)+                    &
                LBC_SIZEA(iside-1,ifld,ihalo,irim)
            ENDIF
            g_LBC_STARTA(iside,ifld,ihalo,irim,mype)=                   &
              LBC_STARTA(iside,ifld,ihalo,irim)

           LENRIMA(ifld,ihalo,irim)=LENRIMA(ifld,ihalo,irim)+           &
                                    LBC_SIZEA(iside,ifld,ihalo,irim)

          ENDDO ! iside
        ENDDO ! ihalo
      ENDDO ! ifld
      ENDDO ! irim

! Now do some comms so that g_LBC_STARTA is filled with the
! value of LBC_STARTA on every processor

      CALL GC_IMAX(4*Nfld_max*NHalo_max*Nrima_max*nproc,nproc,info,     &
                   g_LBC_STARTA)

! And set up a few other variables

      ! Includes one-off orography field at start
      ! Will be updated later, when we know the number of tracer LBCs
      RIM_LOOKUPSA = 13 + num_optional_lbcs_in 

!LL LAM DERIVED SIZES END
      IF (lhook) CALL dr_hook('DERVSIZE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DERVSIZE
