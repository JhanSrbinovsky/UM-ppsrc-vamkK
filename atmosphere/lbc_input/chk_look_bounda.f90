! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine CHK_LOOK_BOUNDA
!
! Purpose : Cross checks values in LOOKUP records of boundary data
!           with model run values
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input

      SUBROUTINE CHK_LOOK_BOUNDA(                                       &
     &  ITEM_LIST,FULL_LOOKUP_BOUNDA,                                   &
! ARGBND Control data calculated from NAMELIST-
        NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
        FIXHD_BOUNDA,INTHD_BOUNDA,LOOKUP_BOUNDA,LOOKUP_COMP_BOUNDA,     &
        REALHD_BOUNDA,                                                  &
     &                           ICODE,CMESSAGE)

      USE Submodel_Mod
      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParVars
      USE Control_Max_Sizes
      USE rimtypes
      USE lbc_mod
      USE lookup_addresses

      USE cppxref_mod, ONLY:                                            &
          ppx_grid_type,ppx_halo_type,                                  &
          ppx_lv_code,ppx_lb_code, ppx_lt_code

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

      INTEGER                                                           &
     &  ITEM_LIST(RIM_LOOKUPSA)                                         &
                                   ! IN: STASH codes of expected items
     &, FULL_LOOKUP_BOUNDA(LEN1_LOOKUP,BOUND_LOOKUPSA)                  &
                                   ! IN: Full LOOKUP record for LBCs
     &, ICODE                     ! OUT : Return code

      CHARACTER(LEN=80)                                                    &
     &  CMESSAGE                  ! OUT : Error message

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

! Functions
      INTEGER                                                           &
     &  EXPPXI                                                          &
     &, GET_FLD_TYPE

! Local variables
      INTEGER                                                           &
     &  variable                                                        &
                          ! Loop counter for variable
     &, code                                                            &
                          ! item_code value for variable
     &, model                                                           &
                          ! model number for variable
     &, section                                                         &
                          ! section number for variable
     &, item                                                            &
                          ! section number for variable
     &, grid_type                                                       &
                          ! grid code for variable
     &, fld_type                                                        &
                          ! P,U or V for variable
     &, halo_type                                                       &
                          ! halo type for variable
     &, level_type                                                      &
                          ! what type of level for variable
     &, bottom_level_code                                               &
                          ! bottom level code
     &, bottom_level                                                    &
                          ! bottom level for variable
     &, top_level_code                                                  &
                          ! top level code
     &, top_level                                                       &
                          ! top level for variable
     &, n_levels_expected                                               &
                          ! number of levels expected
     &, n_levels_lbc                                                    &
                          ! number of levels in file
     &, halo_x_expected                                                 &
                          ! expected size of halo in x
     &, halo_x_lbc                                                      &
                          ! actual size of halo in x
     &, halo_y_expected                                                 &
                          ! expected size of halo in y
     &, halo_y_lbc                                                      &
                          ! actual size of halo in y
     &, size_x_expected                                                 &
                          ! expected size of field in x
     &, size_x                                                          &
                          ! actual size of field in x
     &, size_y_expected                                                 &
                          ! expected size of field in y
     &, size_y                                                          &
                          ! actual size of y
     &, rim_type                                                        &
                          ! Type of RIMWIDTH
     &, rimwidth_expected                                               &
                          ! expected rimwidth
     &, rimwidth_lbc                                                    &
                          ! actual rimwidth
     &, size_expected     ! Expected size

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CHK_LOOK_BOUNDA',zhook_in,zhook_handle)

!=====================================================================

! 1.0 Check that the first field in the file is the orography field

      IF (LOOKUP_BOUNDA(ITEM_CODE,1)  /=  31001) THEN ! Orography
        ICODE=1
        CMESSAGE='CHK_LOOK_BOUNDA : No orography in LBC file'
        GOTO 9999
      ENDIF

! 2.0 Check for the expected number of variables for each time.
!     Bear in mind that RIM_LOOKUPSA contains an extra field
!     (orography) which only occurs at the start of the field.
!     So the actual number of fields written out at each LBC output
!     time is actually RIM_LOOKUPSA-1

      IF (FULL_LOOKUP_BOUNDA(ITEM_CODE,2)  /=                           &
     &    FULL_LOOKUP_BOUNDA(ITEM_CODE,2+RIM_LOOKUPSA-1)) THEN
        WRITE(6,*) 'Wrong number of LBC variables found in LBC file'
        WRITE(6,*) 'Expecting record ',2+RIM_LOOKUPSA-1,' to contain ', &
     &             'STASH item code ',FULL_LOOKUP_BOUNDA(ITEM_CODE,2)
        WRITE(6,*) 'But found item code ',                              &
     &             FULL_LOOKUP_BOUNDA(ITEM_CODE,2+RIM_LOOKUPSA-1)
        ICODE=2
        CMESSAGE='CHK_LOOK_BOUNDA : Wrong number of LBC fields'
        GOTO 9999
      ENDIF

! 3.0 Now check the header record for each required variable

      DO variable=1,RIM_LOOKUPSA

        code=ITEM_LIST(variable)  ! item_code value for variable

        IF (LOOKUP_BOUNDA(ITEM_CODE,variable)  /=  code) THEN
          WRITE(6,*) 'Unexpected field in LBC file'
          WRITE(6,*) 'Field ',variable,' was expected to be ',          &
     &               code,' but found ',                                &
     &               LOOKUP_BOUNDA(ITEM_CODE,variable)

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=3
          CMESSAGE='CHK_LOOK_BOUNDA : Unexpected field in LBC file'
          GOTO 9999
        ENDIF

        model=LOOKUP_BOUNDA(MODEL_CODE,variable)
        item=MOD(LOOKUP_BOUNDA(ITEM_CODE,variable),1000)
        section=(LOOKUP_BOUNDA(ITEM_CODE,variable)-item)/1000

! DEPENDS ON: exppxi
        grid_type=EXPPXI(model,section,item,ppx_grid_type,              &
     &                   ICODE, CMESSAGE)
! DEPENDS ON: get_fld_type
        fld_type=GET_FLD_TYPE(grid_type)

! DEPENDS ON: exppxi
        halo_type=EXPPXI(model,section,item,ppx_halo_type,              &
     &                   ICODE, CMESSAGE)

! DEPENDS ON: exppxi
        level_type=EXPPXI(model,section,item,ppx_lv_code,               &
     &                    ICODE, CMESSAGE)
        IF (level_type  ==  5) THEN
          n_levels_expected=1
        ELSE
! DEPENDS ON: exppxi
          bottom_level_code=EXPPXI(model,section,item,ppx_lb_code,      &
     &                          ICODE, CMESSAGE)
! DEPENDS ON: exppxi
          top_level_code=EXPPXI(model,section,item,ppx_lt_code,         &
     &                          ICODE, CMESSAGE)
! DEPENDS ON: levcod
          CALL LEVCOD(bottom_level_code,bottom_level,ICODE,CMESSAGE)
! DEPENDS ON: levcod
          CALL LEVCOD(top_level_code,top_level,ICODE,CMESSAGE)
          n_levels_expected=top_level-bottom_level+1
        ENDIF

        halo_x_expected=halosize(1,halo_type)
        halo_y_expected=halosize(2,halo_type)
        size_x_expected=glsize(1,fld_type)
        size_y_expected=glsize(2,fld_type)

! This is dependent on how the LBC is setup.  At the moment we dont correct.
!        IF (fld_type  ==  fld_type_v .OR. fld_type == fld_type_p) THEN
        IF ( .NOT. l_vatpoles .AND. (fld_type  ==  fld_type_u) ) THEN
          size_x_expected=size_x_expected-1
        ENDIF

        IF (LOOKUP_BOUNDA(ITEM_CODE,variable)  ==  31001) THEN
          ! Orography
          rim_type=rima_type_orog
        ELSE
          rim_type=rima_type_norm
        ENDIF

        rimwidth_expected=RIMWIDTHA(rim_type)

        halo_x_lbc=MOD(LOOKUP_BOUNDA(LBUSER3,variable),100)
        halo_y_lbc=MOD(LOOKUP_BOUNDA(LBUSER3,variable)-halo_x_lbc,      &
     &                 10000)/100
        rimwidth_lbc=MOD((LOOKUP_BOUNDA(LBUSER3,variable)-              &
     &                 halo_x_lbc-halo_y_lbc*100),1000000)/10000

        n_levels_lbc=LOOKUP_BOUNDA(LBHEM,variable)-100

        size_expected=global_LENRIMA(fld_type,halo_type,rim_type)*      &
     &                n_levels_expected

        IF (n_levels_lbc  /=  n_levels_expected) THEN
          WRITE(6,*) 'Wrong number of levels for LBC field ',variable
          WRITE(6,*) 'Expected ',n_levels_expected,' levels but found ',&
     &               n_levels_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=4
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong number of levels'
          GOTO 9999
        ENDIF

        IF ((halo_x_lbc  /=  halo_x_expected) .OR.                      &
     &      (halo_y_lbc  /=  halo_y_expected)) THEN
          WRITE(6,*) 'Incorrect halos for LBC field ',variable
          WRITE(6,*) 'Expected halo_x= ',halo_x_expected,               &
     &               ' and halo_y= ',halo_y_expected
          WRITE(6,*) 'but found halo_x= ',halo_x_lbc,                   &
     &               ' and halo_y= ',halo_y_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=5
          CMESSAGE='CHK_LOOK_BOUNDA : Incorrect halos'
          GOTO 9999
        ENDIF

        IF (rimwidth_lbc  /=  rimwidth_expected) THEN
          WRITE(6,*) 'Wrong RIMWIDTH for LBC field ',variable
          WRITE(6,*) 'Expected RIMWIDTH= ',rimwidth_expected,           &
     &               'but found RIMWIDTH= ',rimwidth_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=6
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong RIMWIDTH'
          GOTO 9999
        ENDIF

        size_x=LOOKUP_BOUNDA(LBNPT,variable)
        size_y=LOOKUP_BOUNDA(LBROW,variable)

        IF ((size_x  /=                                                 &
     &       size_x_expected) .OR.                                      &
     &      (size_y  /=                                                 &
     &       size_y_expected)) THEN
          WRITE(6,*) 'Incorrect dimensions for LBC field ',variable
          WRITE(6,*) 'Expected ROW_LENGTH= ',size_x_expected,' and ',   &
     &               'ROWS= ',size_y_expected
          WRITE(6,*) 'But found ROW_LENGTH= ',                          &
     &               size_x,' and ROWS= ',size_y

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=7
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong dimensions'
          GOTO 9999
        ENDIF

        IF (LOOKUP_BOUNDA(LBLREC,variable)  /=  size_expected) THEN
          WRITE(6,*) 'Wrong size for LBC field ',variable
          WRITE(6,*) 'Expected size was ',size_expected,' but found ',  &
     &               LOOKUP_BOUNDA(LBLREC,variable)

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=8
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong size'
          GOTO 9999
        ENDIF

      ENDDO ! variable

 9999 CONTINUE

      IF (lhook) CALL dr_hook('CHK_LOOK_BOUNDA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE CHK_LOOK_BOUNDA
