! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! Takes a general decomposed field on many processors and gathers it
! to a single processor. Also supports the use case of a deferred gather 
! to an IO server, when the final target of the gather is output to a dump
! file.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

SUBROUTINE general_gather_field (                               &
  local_field, global_field,                                    &
  local_size, global_size,                                      &
  levels,                                                       &
  addr_info,                                                    &
  gather_pe,                                                    &
  asyncOutputHandle,                                            &
  icode, cmessage)

  USE IOS_Constants
  USE IOS_stash_common
  USE IOS_model_geometry, ONLY :                                &
      ios_dump_init
  USE IOS_dump
  USE mask_compression, ONLY :                                  &
      compress_to_mask,                                         &
      expand_from_mask
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE gcom_mod, ONLY: gc_none
  USE UM_ParVars
  USE printStatus_mod, ONLY : &
      printstatus,            &
      prstatus_diag

  USE cppxref_mod   ! 25+ ppx_XXX items required
  IMPLICIT NONE

! Subroutine arguments:

  INTEGER, INTENT(IN)  ::  global_size         ! Size of GLOBAL FIELD
  INTEGER, INTENT(IN)  ::  levels              ! How many levels of data to do
  INTEGER, INTENT(IN)  ::  gather_pe(levels)   ! Which PE to gather each level to
  INTEGER, INTENT(IN)  ::  asyncOutputHandle   ! Indicates this gather should 
                                               ! go to an IOS directly via 
                                               ! an async op 
                                               ! -1  : no IOS op 
                                               ! >=0 : handle of the async op
  INTEGER, INTENT(OUT) ::  local_size          ! Size of LOCAL_FIELD
  INTEGER, INTENT(OUT) ::  icode               ! Return code, 0=OK
  
! Required for dimensioning ADDR_INFO
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

  INTEGER, INTENT(IN)  ::  addr_info(d1_list_len) ! Addressing info about field
  
  REAL, INTENT(IN)     ::  local_field(*)         ! My local part of field
  REAL, INTENT(OUT)    ::  global_field(global_size,*) 
                                                  ! Array to gather field to

  CHARACTER(LEN=80)    ::  cmessage               ! Error message

! Parameters and common blocks
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
!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!

      LOGICAL                                                           &
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)                                  &
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

      COMMON /Atmos_LSM_Common/                                         &
     &  atmos_landmask                                                  &
     &, atmos_landmask_local                                            &
     &, atmos_number_of_landpts

! End of comdeck ATM_LSM
! Start i_stgfld

! Description:
!   This file contains an interface to STASH_GATHER_FIELD and
!   must be included whenever this routines is used so as to
!   get declarations of optional arguments correct.
!
      INTERFACE
        SUBROUTINE STASH_GATHER_FIELD (                                 &
     &    LOCAL_FIELD , GLOBAL_FIELD ,                                  &
     &    LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                              &
     &    GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,   &
     &    GRIDTYPE_CODE ,HALO_TYPE,                                     &
     &    GATHER_PE,                                                    &
     &    DATA_EXTRACTED,                                               &
     &    PACKING, IM_IDENT, PACKING_TYPE,                        &
     &    NUM_OUT,                                                      &
     &    COMP_ACCRCY, loc_RMDI,                                        &
     &    ICODE, CMESSAGE)

        INTEGER, INTENT(IN) ::                                          &
     &    LOCAL_SIZE                                                    &
                          ! IN: size of level of LOCAL_FIELD
     &  , GLOBAL_SIZE                                                   &
                          ! IN: size of level of GLOBAL_FIELD
     &  , LEVELS                                                        &
                          ! IN: number of levels
     &  , GLOBAL_NORTH                                                  &
                          ! IN: specification of subdomain boundaries
     &  , GLOBAL_EAST_IN                                                &
                          ! IN: ""
     &  , GLOBAL_SOUTH                                                  &
                          ! IN: ""
     &  , GLOBAL_WEST                                                   &
                          ! IN: ""
     &  , GRIDTYPE_CODE                                                 &
                          ! IN: indicates the type of grid output
     &  , HALO_TYPE                                                     &
                          ! IN: type of halo on this field
     &  , GATHER_PE       ! IN: the PE to gather the global field to

        INTEGER, INTENT(OUT) ::                                         &
     &    ICODE           ! OUT: return code, 0=OK
!
! Optional Arguments to handle the COEX packing if necessary
!
        LOGICAL, INTENT(IN), OPTIONAL ::   &
     &    PACKING
                          ! IN: Set .true. if packing of the input
                          !     field is to be packed

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    IM_IDENT        ! IN: Internal model identifier

        INTEGER, INTENT(INOUT), OPTIONAL ::                             &
     &    PACKING_TYPE    ! IN/OUT: This flag is zero on input,
                          !         then stash packing is selected,
                          !         and the routine returns the
                          !         packing flag.
                          !
                          !         If the variable is set to 1 on input
                          !         then 32-bit packing for dumpfiles
                          !         is selected

        INTEGER, INTENT(OUT), OPTIONAL ::                               &
     &    NUM_OUT         ! OUT: Number of 32-bit IBM words in the
                          !      Packed field for WDGOS packing

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    COMP_ACCRCY     ! IN: Packing Accuracy in Power of 2

        REAL, INTENT(IN), OPTIONAL ::                                   &
     &    loc_RMDI        ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
        LOGICAL, INTENT(IN) ::                                          &
     &    DATA_EXTRACTED  ! IN: TRUE if the data in LOCAL_FIELD has
                          !     already been extracted, or FALSE if
                          !     the extraction must be done here.

        REAL, INTENT(IN) ::                                             &
     &    LOCAL_FIELD(LOCAL_SIZE,LEVELS)
                          ! IN : local data

        REAL, INTENT(OUT) ::                                            &
     &    GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
                          ! OUT : (PE GATHER_PE only) - full gathered
                          !       field

        CHARACTER(LEN=*), INTENT(OUT) ::                                   &
     &    CMESSAGE        ! OUT: Error message if ICODE .NE. 0

        END SUBROUTINE STASH_GATHER_FIELD
      END INTERFACE
! End i_stgfld

! Local variables
INTEGER   ::                                                      &
  grid_type                                                       &
            ! grid type of field being gathered
, grid_code                                                       &
            ! ppx grid code of field being gathered
, halo_type                                                       &
            ! halo type of field being gathered
, k                                                               &
       ! loop over levels
, my_k                                                            &
       ! value of k for GLOBAL_FIELD on this PE
, my_k_temp                                                       &
       ! Temp. my_k for safe references
, info                                                            &
       ! return code from GCOM routines
, dummy                                                           &
       ! dummy variables - ignored return values
, north,south,east,west                                           &
            ! domain limits for STASH output
, mean_type                                                       &
            ! spatial meaning type on diagnostic
, iproc                                                           &
            ! Loop variable over processors
, i         ! Loop variable for debugging

INTEGER       :: get_fld_type  ! function for finding field type
INTEGER       :: ios_subdomain   
LOGICAL, SAVE :: haveInitialisedIOSLandMask=.FALSE.
REAL   ::                                                         &
  buf_expand(max2dfieldsize)                                      &
, buf_expand_local(max2dfieldsize)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!===================================================================

IF (lhook) CALL dr_hook('GENERAL_GATHER_FIELD',zhook_in,zhook_handle)

grid_code=addr_info(d1_grid_type)
! DEPENDS ON: get_fld_type
grid_type=get_fld_type(grid_code)
halo_type=addr_info(d1_halo_type)

!-------------------------------------------------------------------

! Timeseries data
IF ((addr_info(d1_object_type)  ==  diagnostic) .AND.             &
    ((addr_info(d1_proc_no_code)  ==  st_time_series_code).OR.    &
    (addr_info(d1_proc_no_code)  ==  st_time_series_mean))) THEN

  ! Multi-level fields not supported here
  IF (levels  >   1) THEN
    WRITE(6,*) 'GENERAL_GATHER_FIELD : Cannot have more than ',   &
        '1 level for gathering timeseries data.'
    icode=1
    cmessage='GENERAL_GATHER_FIELD : Multi-level timeseries field'
    GO TO 9999
  END IF

  IF (asyncOutputHandle /= -1) THEN
    WRITE(6,*) 'GENERAL_GATHER_FIELD : Cannot have time series ', &
        'in async dumps using IOS'
    icode=1
    cmessage= &
        'GENERAL_GATHER_FIELD : timeseries field with IOS not allowed'
    GO TO 9999
  END IF

  ! Copy the data to GLOBAL_FIELD from PE 0
  IF (mype  ==  0) THEN
    info=gc_none
    CALL gc_rsend(99,global_size,gather_pe,info,global_field,     &
        local_field)
  END IF

  IF (mype  ==  gather_pe(1)) THEN
    info=gc_none
    CALL gc_rrecv(99,global_size,0,info,global_field,local_field)
  END IF

  local_size=global_size

!-------------------------------------------------------------------

! Surface (land points only) fields

ELSE IF (grid_code  ==  ppx_atm_compressed) THEN
   
  IF (.NOT.haveInitialisedIOSLandMask.AND. &
      (asyncOutputHandle /= -1)) THEN
     
    IF (mype == 0 .AND. printstatus >= prstatus_diag) THEN
       WRITE(6,*)'Initialise landmask on IO Server...'
    END IF
    CALL ios_dump_init(atmos_landmask,&
        glsize(1,grid_type)*glsize(2,grid_type)&
        )
    haveInitialisedIOSLandMask=.TRUE.
  END IF

  my_k=0
  DO k=1,levels
! Unpack the local field out to full (local) field size and
! put this into the array buf_expand_local

    CALL expand_from_mask(buf_expand_local,          &
        local_field(1+(k-1)*local_land_field),    &
        atmos_landmask_local,                     &
        lasize(1,grid_type,halo_type)*            &
        lasize(2,grid_type,halo_type),            &
        dummy)
  
! Now gather in all the processors local fields into the global
! field (array buf_expand) or send to IO server

    IF (asyncOutputHandle /= -1) THEN 
      CALL ios_dump_pack_data ( &
       asyncOutputHandle, &                 ! Communicatiosn handle
       buf_expand_local, &                  ! Data In
       arg_grid_type=grid_code, &           ! field type
       arg_halo_type=halo_type, &           ! halo type 
       arg_preprocess_flag=ios_stash_preprocess, &! preprocess active
       arg_subdomain_flag=IOS_Full_Field, & ! subdomain control
       arg_landmask_compress=.TRUE. &
       )
    ELSE
! DEPENDS ON: gather_field
      CALL gather_field(buf_expand_local,buf_expand,    &
          lasize(1,grid_type,halo_type),                &
          lasize(2,grid_type,halo_type),                &
          glsize(1,grid_type),glsize(2,grid_type),      &
          grid_type,halo_type,                          &
          gather_pe(k),gc_all_proc_group,               &
          icode,cmessage)

      IF (icode  /=  0) THEN
        WRITE(6,*)                                                    &
            'GENERAL_GATHER_FIELD : Error detected in call to ',      &
            'GATHER_FIELD (land point field)'
        WRITE(6,*) 'Return code : ',icode
        WRITE(6,*) 'Error message : ',cmessage
        
        icode=1
        cmessage='GENERAL_GATHER_FIELD : Error gathering field'
        GO TO 9999
      END IF

! And now pack the global field (buf_expand) back to land points
! and put into the array GLOBAL_FIELD.

      IF (mype  ==  gather_pe(k)) THEN
        my_k=my_k+1

        CALL compress_to_mask(buf_expand,global_field(1,my_k), &
            atmos_landmask,                                   &
            glsize(1,grid_type)*glsize(2,grid_type),          &
            dummy)
      END IF
    END IF
  END DO ! k : loop over levels

  local_size = local_land_field

!-------------------------------------------------------------------

! Atmosphere/ocean Lateral boundary fields

ELSE IF                                                           &
  (grid_code  ==  ppx_ocn_rim  .OR.                               &
   grid_code  ==  ppx_atm_rim) THEN

    WRITE(6,*)                                                    &
      'GENERAL_GATHER_FIELD : Error attempting to gather ',       &
      'boundary data. Data not stored in dump from vn5.3'
    icode=3
    cmessage='GENERAL_GATHER_FIELD : Error gathering LBCs'
    GO TO 9999

!-------------------------------------------------------------------

! "Normal" fields

ELSE IF                                                           &
!     atmosphere grids
  ((grid_code  ==  ppx_atm_tall)  .OR.                            &
   (grid_code  ==  ppx_atm_tland) .OR.                            &
   (grid_code  ==  ppx_atm_tsea)  .OR.                            &
   (grid_code  ==  ppx_atm_uall) .OR.                             &
   (grid_code  ==  ppx_atm_cuall) .OR.                            &
   (grid_code  ==  ppx_atm_cvall) .OR.                            &
   (grid_code  ==  ppx_atm_tzonal) .OR.                           &
   (grid_code  ==  ppx_atm_uzonal) .OR.                           &
   (grid_code  ==  ppx_atm_tmerid) .OR.                           &
   (grid_code  ==  ppx_atm_compressed) .OR.                       &
   (grid_code  ==  ppx_atm_ozone) .OR.                            &
   (grid_code  ==  ppx_atm_river) .OR.                            &
!     ocean grids
   (grid_code  ==  ppx_ocn_tall) .OR.                             &
   (grid_code  ==  ppx_ocn_tfield) .OR.                           &
   (grid_code  ==  ppx_ocn_tzonal) .OR.                           &
   (grid_code  ==  ppx_ocn_uzonal) .OR.                           &
   (grid_code  ==  ppx_ocn_tmerid) .OR.                           &
   (grid_code  ==  ppx_ocn_uall) .OR.                             &
   (grid_code  ==  ppx_ocn_cuall) .OR.                            &
   (grid_code  ==  ppx_ocn_ufield) .OR.                           &
   (grid_code  ==  ppx_ocn_umerid) .OR.                           &
   (grid_code  ==  ppx_ocn_cvall) .OR.                            &
!     Full WAM Wave Model Grid - Land Sea Mask only
   (grid_code  ==  ppx_wam_all))                                  &
   THEN

  local_size=addr_info(d1_length)/addr_info(d1_no_levels)
  IOS_Subdomain=IOS_Full_Field

  IF (addr_info(d1_object_type)  ==  diagnostic) THEN
    north=addr_info(d1_north_code)
    south=addr_info(d1_south_code)
    east=addr_info(d1_east_code)
    west=addr_info(d1_west_code)

    mean_type=addr_info(d1_gridpoint_code)/10
    IF (mean_type  ==  2) THEN ! zonal mean
      IOS_Subdomain=IOS_Partial_Field
      east=west
    ELSE IF (mean_type  ==  3) THEN ! meridional mean
      IOS_Subdomain=IOS_Partial_Field
      north=south
    ELSE IF (mean_type  >=  4) THEN ! field/global mean
      IOS_Subdomain=IOS_Partial_Field
      east=west
      north=south
    END IF

  ELSE
    north=glsize(2,grid_type)
    west=1
    east=glsize(1,grid_type)
    south=1
  END IF

  my_k=0
  DO k=1,levels

    IF (mype  ==  gather_pe(k)) THEN
      my_k      = my_k+1
      my_k_temp = my_k
    ELSE
      my_k_temp = 1
    END IF

    IF (asyncOutputHandle /= -1) THEN

      CALL ios_dump_pack_data (            &
          asyncOutputHandle,               & ! Communications handle
          local_field(1+(k-1)*local_size : &
          (k  )*local_size),               & ! Data In
          arg_grid_type=grid_code,         & ! field type
          arg_halo_type=halo_type,         & ! halo type 
          arg_preprocess_flag=             & 
          ios_stash_preprocess,            & ! preprocess active
          arg_subdomain_flag=IOS_Subdomain,& ! subdomain control
          arg_S_boundary=south,            &
          arg_N_boundary=north,            &
          arg_W_boundary=west,             &
          arg_E_boundary=east,             &
          arg_landmask_compress=.FALSE.    &
          )
    ELSE
!       STASH_GATHER_FIELD can distribute whole fields, or subarea
!       fields

! DEPENDS ON: stash_gather_field
      CALL stash_gather_field(                                          &
          local_field(1+(k-1)*local_size),global_field(1,my_k_temp),    &
          local_size,global_size,1,                                     &
          north,east,south,west,                                        &
          grid_code,halo_type,gather_pe(k),.TRUE.,                      &
          icode=icode, cmessage=cmessage)
      
      IF (icode  /=  0) THEN
        WRITE(6,*)                                                      &
            'GENERAL_GATHER_FIELD : Error detected in call to ',        &
            'STASH_GATHER_FIELD'
        WRITE(6,*) 'Return code : ',icode
        WRITE(6,*) 'Error message : ',cmessage

        icode=4
        cmessage='GENERAL_GATHER_FIELD : Error gathering field'
        GO TO 9999
      END IF
    END IF
  END DO ! k : loop over levels

!-------------------------------------------------------------------
! Any other type of field
ELSE

  icode=10
  cmessage='GENERAL_GATHER_FIELD : Field type not recognized'

END IF

 9999 CONTINUE

IF (lhook) CALL dr_hook('GENERAL_GATHER_FIELD',zhook_out,zhook_handle)
RETURN

END SUBROUTINE general_gather_field

