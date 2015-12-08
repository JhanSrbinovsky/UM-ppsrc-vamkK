! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Obtain UKCA-MODE input to UKCA_RADAER from D1.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_get(            &
                           ierr        &
  ,                        cmessage    &
  ,                        first_call, &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                           ukca_radaer &
  )

USE Submodel_Mod
USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE atm_fields_bounds_mod, ONLY: tdims_s
USE ukca_radaer_struct_mod
      
USE UM_ParVars
IMPLICIT NONE

!
! Arguments with intent(in)
!
!
! Error indicator (0 is OK, >0 error)
!
INTEGER ierr

!
! Error message if ierr is larger than 0
!
CHARACTER (len=256) :: cmessage

!
! Logical indicating first call: complete setup has to be done
!
LOGICAL first_call

!
! D1 and related variables
!
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
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer

!
! Local variables
!
INTEGER i, &
        j
INTEGER i_obj     

!
! Atmosphere submodel index
!
INTEGER m_atm_modl

!
! STASH section where UKCA diagnostics are expected to reside.
!  Mass-mixing ratios/numbers and misc diagnostics (densities,
!  volumes) can be in different sections.
!
INTEGER, PARAMETER :: ukca_section_mmr = 34
INTEGER, PARAMETER :: ukca_section_oth = 38

!
! Logical for checking whether all diags have been found
!
LOGICAL l_diag_missing
 
!    
! Local buffer read from D1 and its size.
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: buffer
INTEGER buffer_size

!
! Halo sizes
!
INTEGER halo_x
INTEGER prev_halo_x
INTEGER halo_y
INTEGER prev_halo_y

! Tracer Levels
INTEGER tr_levs

!
! Variables for looking for tagged diagnostics in D1.
!
INTEGER tag, ptd1, section, item, levs, len, addr, halotyp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
   
!
!
!

IF (lhook) CALL dr_hook('UKCA_RADAER_GET', zhook_in, zhook_handle)

ierr = 0

m_atm_modl = SUBMODEL_FOR_SM(a_im)

! Set tracer levels
tr_levs = tdims_s%k_end - tdims_s%k_start + 1
!
! For first call, check that all diagnostics required are available
! at the expected dimensions. Retain their D1 address for later
! calls.
!

IF (first_call) THEN
  
  !
  ! Initialise all d1 addresses in mode and component
  ! information structures to "not found" (-1)
  !
  DO j = 1, ukca_radaer%n_mode
    ukca_radaer%d1_address_dry(j) = -1
    ukca_radaer%d1_address_wet(j) = -1
    ukca_radaer%d1_address_rho(j) = -1
  END DO ! j
  
  DO i = 1, ukca_radaer%n_cpnt
    ukca_radaer%d1_address_mmr(i) = -1
    ukca_radaer%d1_address_cvl(i) = -1
  END DO ! i
  
  !
  ! Go through D1 and retain some information about the
  ! diagnostics we need. They should be tagged with tag 98.
  !
  DO i_obj = 1, totitems
    tag = STLIST(st_macrotag,i_obj) - &
          (1000*(STLIST(st_macrotag,i_obj)/1000))

    !
    ! Only consider items with the expected tag.
    !
    IF (tag == 98) THEN
      
      ptd1    = STLIST(st_D1pos, i_obj)
      section = d1_addr(d1_section,   ptd1, m_atm_modl)
      item    = d1_addr(d1_item,      ptd1, m_atm_modl)
      levs    = d1_addr(d1_no_levels, ptd1, m_atm_modl)
      len     = d1_addr(d1_length,    ptd1, m_atm_modl)
      addr    = d1_addr(d1_address,   ptd1, m_atm_modl)
      halotyp = d1_addr(d1_halo_type, ptd1, m_atm_modl)

    
      IF (section == ukca_section_mmr) THEN
      
        !
        ! We are in the section where the component mass-mixing
        ! ratios and modal numbers should reside.
        !
        
        DO i = 1, ukca_radaer%n_cpnt
      
          IF (item == ukca_radaer%stashcode_mmr(i)) THEN

            ukca_radaer%d1_address_mmr(i) = addr
            ukca_radaer%d1_nlevs_mmr(i)   = levs
            ukca_radaer%d1_length_mmr(i)  = len
            ! Mass-mixing ratios have halos
            ukca_radaer%d1_halo_type_mmr(i) = halotyp
        
          END IF 
      
        END DO ! i (components)
      
        DO j = 1, ukca_radaer%n_mode
      
          IF (item == ukca_radaer%stashcode_nbr(j)) THEN
          
            ukca_radaer%d1_address_nbr(j) = addr
            ukca_radaer%d1_nlevs_nbr(j)   = levs
            ukca_radaer%d1_length_nbr(j)  = len
            ! Modal numbers have halos
            ukca_radaer%d1_halo_type_nbr(j) = halotyp
          
          END IF
      
        END DO ! j (modes)
      
      END IF ! ukca_section_mmr
    
      IF (section == ukca_section_oth) THEN
      
        !
        ! We are in the section where dry and wet diameters,
        ! modal densities and volumes, and component volumes
        ! should reside.
        ! This is also where component volumes are expected.
        !
      
        DO j = 1, ukca_radaer%n_mode
        
          IF (item == ukca_radaer%stashcode_dry(j)) THEN
          
            ukca_radaer%d1_address_dry(j) = addr
            ukca_radaer%d1_nlevs_dry(j)   = levs
            ukca_radaer%d1_length_dry(j)  = len
        
          END IF
        
          IF (ukca_radaer%l_soluble(j)) THEN
        
            IF (item == ukca_radaer%stashcode_wet(j)) THEN
          
              ukca_radaer%d1_address_wet(j) = addr
              ukca_radaer%d1_nlevs_wet(j)   = levs
              ukca_radaer%d1_length_wet(j)  = len
          
            ELSE IF (item == ukca_radaer%stashcode_wtv(j)) THEN
            
              ukca_radaer%d1_address_wtv(j) = addr
              ukca_radaer%d1_nlevs_wtv(j)   = levs
              ukca_radaer%d1_length_wtv(j)  = len
          
            END IF
        
          END IF
        
          IF (item == ukca_radaer%stashcode_rho(j)) THEN

            ukca_radaer%d1_address_rho(j) = addr
            ukca_radaer%d1_nlevs_rho(j)   = levs
            ukca_radaer%d1_length_rho(j)  = len
          
          END IF
        
        END DO ! j (modes)
      
        DO i = 1, ukca_radaer%n_cpnt
      
          IF (item == ukca_radaer%stashcode_cvl(i)) THEN
          
            ukca_radaer%d1_address_cvl(i) = addr
            ukca_radaer%d1_nlevs_cvl(i)   = levs
            ukca_radaer%d1_length_cvl(i)  = len

          END IF
        
        END DO ! i (components)
      
      END IF ! ukca_section_oth
    
    END IF ! tag
    
  END DO ! i_obj (D1 items)
  
  !
  ! Check that all diagnostics have been found.
  !
  l_diag_missing = .false.
  
  DO j = 1, ukca_radaer%n_mode
    
    IF (ukca_radaer%d1_address_dry(j) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                        &
            ukca_radaer%stashcode_dry(j), ' not found in D1.'
      l_diag_missing = .true.
    
    ELSE IF (ukca_radaer%l_soluble(j) .and.                       &
        ukca_radaer%d1_address_wet(j) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                            &
             ukca_radaer%stashcode_wet(j), ' not found in D1.'
      l_diag_missing = .true.
    
    ELSE IF (ukca_radaer%l_soluble(j) .and.                       &
        ukca_radaer%d1_address_wtv(j) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                          &
              ukca_radaer%stashcode_wtv(j),' not found in D1.'
      l_diag_missing = .true.
              
    ELSE IF (ukca_radaer%d1_address_rho(j) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                          &
              ukca_radaer%stashcode_rho(j),' not found in D1.'
      l_diag_missing = .true.
    
    ELSE IF (ukca_radaer%d1_address_nbr(j) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                          &
               ukca_radaer%stashcode_nbr(j),' not found in D1.'
      l_diag_missing = .true.
    
    END IF
    
  END DO ! j
  
  IF (l_diag_missing) THEN
    ierr = 702
    cmessage =                                                    &
      'ukca_radaer_get: Diag(s) needed for UKCA are missing from D1.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
  DO i = 1, ukca_radaer%n_cpnt
    
    IF (ukca_radaer%d1_address_mmr(i) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                          &
           ukca_radaer%stashcode_mmr(i),' not found in D1.'
      l_diag_missing = .true.
      
    ELSE IF (ukca_radaer%d1_address_cvl(i) == -1) THEN
      WRITE(6,'(A,I6,A)') 'Diagnostic ',                          &
           ukca_radaer%stashcode_cvl(i), ' not found in D1.'
      l_diag_missing = .true.
      
    END IF
    
  END DO ! i
  
  IF (l_diag_missing) THEN
    ierr = 702
    cmessage =                                                    &
      'ukca_radaer_get: Diag(s) needed for UKCA are missing from D1.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
END IF ! first_call

!
! At this point, we know where to look for the requested D1 items.
! Retrieve. The retrieval takes place into a temporary buffer that
! is reshaped.
!
buffer_size = row_length * rows * model_levels
ALLOCATE(buffer(row_length, rows, model_levels))

DO j = 1, ukca_radaer%n_mode

  !
  ! Modal dry diameter.
  ! Check number of levels and total length.
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_dry(j) /= model_levels) THEN
    ierr = 705
    WRITE(6,'(2(A,I5))')'Expecting ',model_levels,' levels, got ',&
                ukca_radaer%d1_nlevs_dry(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  IF (ukca_radaer%d1_length_dry(j) /= buffer_size) THEN
    ierr = 706
    WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,              &
             ' elements, got ', ukca_radaer%d1_length_dry(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  buffer = reshape(D1(ukca_radaer%d1_address_dry(j) :             &
                      ukca_radaer%d1_address_dry(j) +             &
                      ukca_radaer%d1_length_dry(j) - 1),          &
                   (/row_length, rows, model_levels/))
  ukca_radaer%dry_diam(:,:,:,j) = buffer(:,:,:)

  IF (ukca_radaer%l_soluble(j)) THEN

    !
    ! Modal wet diameter.
    ! Check number of levels and total length.
    ! If ok, read into buffer and reshape into target array.
    !
    IF (ukca_radaer%d1_nlevs_wet(j) /= model_levels) THEN
      ierr = 705
      WRITE(6,'(2(A,I5))') 'Expecting ', model_levels,            &
            ' levels, got ', ukca_radaer%d1_nlevs_wet(j)
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected number of levels in D1 diag.'
      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
      END IF
      RETURN
    END IF

    IF (ukca_radaer%d1_length_wet(j) /= buffer_size) THEN
      ierr = 706
      WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,            &
               ' elements, got ',ukca_radaer%d1_length_wet(j)
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected total size of D1 diag.'
      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
      END IF
      RETURN
    END IF

    buffer = reshape(D1(ukca_radaer%d1_address_wet(j) :           &
                        ukca_radaer%d1_address_wet(j) +           &
                        ukca_radaer%d1_length_wet(j) - 1),        &
                     (/row_length, rows, model_levels/))
    ukca_radaer%wet_diam(:,:,:,j) = buffer(:,:,:)
    
    !
    ! Volume of water.
    ! Check number of levels and total length.
    ! If ok, read into buffer and reshape into target array.
    !
    IF (ukca_radaer%d1_nlevs_wtv(j) /= model_levels) THEN
      ierr = 705
      WRITE(6,'(2(A,I5))') 'Expecting ', model_levels,            &
            ' levels, got ',ukca_radaer%d1_nlevs_wtv(j)
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected number of levels in D1 diag.'
      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
      END IF
      RETURN
    END IF

    IF (ukca_radaer%d1_length_wtv(j) /= buffer_size) THEN
      ierr = 706 
      WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,            &
               ' elements, got ',ukca_radaer%d1_length_wtv(j)
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected total size of D1 diag.'
      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
      END IF
      RETURN
    END IF
    
    buffer = reshape(D1(ukca_radaer%d1_address_wtv(j) :           &
                        ukca_radaer%d1_address_wtv(j) +           &
                        ukca_radaer%d1_length_wtv(j) - 1),        &
                    (/row_length, rows, model_levels/))
    ukca_radaer%modal_wtv(:,:,:,j) = buffer(:,:,:)
    
  ELSE
    
    !
    ! Insoluble modes: copy the dry diameter into the wet 
    ! diameter, for consistency, and set the volume of water
    ! to zero.
    !
    ukca_radaer%wet_diam(:,:,:,j) = ukca_radaer%dry_diam(:,:,:,j)
    ukca_radaer%modal_wtv(:,:,:,j) = 0.0E+00
    
  END IF ! l_soluble

  !
  ! Modal density.
  ! Check number of levels and total length.
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_rho(j) /= model_levels) THEN
    ierr = 705
    WRITE(6,'(2(A,I5))') 'Expecting ', model_levels,             &
            ' levels, got ',ukca_radaer%d1_nlevs_rho(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  IF (ukca_radaer%d1_length_rho(j) /= buffer_size) THEN
    ierr = 706
    WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,            &
            ' elements, got ', ukca_radaer%d1_length_rho(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  buffer = reshape(D1(ukca_radaer%d1_address_rho(j) :             &
                      ukca_radaer%d1_address_rho(j) +             &
                      ukca_radaer%d1_length_rho(j) - 1),          &
                   (/row_length, rows, model_levels/))
  ukca_radaer%modal_rho(:,:,:,j) = buffer(:,:,:)
  
  !
  ! Initialise modal volume to that of water (the latter has been
  ! set to zero for insoluble modes above).
  ! The other components will be added in the loop below.
  !
  ukca_radaer%modal_vol(:,:,:,j) = ukca_radaer%modal_wtv(:,:,:,j)
  
END DO ! j

!
! Variables depending on components now.
!
DO i = 1, ukca_radaer%n_cpnt
  
  !
  ! Component fractional volumes.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_cvl(i) /= model_levels) THEN
    ierr = 705
    WRITE(6,'(2(A,I5))') 'Expecting ', model_levels,             &
           ' levels, got ', ukca_radaer%d1_nlevs_cvl(i)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  IF (ukca_radaer%d1_length_cvl(i) /= buffer_size) THEN
    ierr = 706
    WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,            &
         ' elements, got ', ukca_radaer%d1_length_cvl(i)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  buffer = reshape(D1(ukca_radaer%d1_address_cvl(i) :             &
                      ukca_radaer%d1_address_cvl(i) +             &
                      ukca_radaer%d1_length_cvl(i) - 1),          &
                   (/row_length, rows, model_levels/))
  ukca_radaer%comp_vol(:,:,:,i) = buffer(:,:,:)        
  
END DO ! i

!
! Update the volume of each mode by adding the volume
! of each component within that mode.
!
DO j = 1, ukca_radaer%n_mode

  DO i = 1, ukca_radaer%n_cpnt_in_mode(j)
  
    ukca_radaer%modal_vol(:,:,:,j) = &
      ukca_radaer%modal_vol(:,:,:,j) + &
      ukca_radaer%comp_vol(:,:,:,ukca_radaer%i_cpnt_index(i, j))
  
  END DO ! i

END DO ! j

!
! Deallocate the buffer. It will be reallocated to account 
! for halos in mass mixing ratios and modal number.
!
DEALLOCATE(buffer)
buffer_size = 0

prev_halo_x = -1
prev_halo_y = -1

DO i = 1, ukca_radaer%n_cpnt

  !
  ! Mass mixing ratio.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_mmr(i) /= tr_levs ) THEN
    ierr = 707
    WRITE(6,'(2(A,I5))') 'Expecting ', tr_levs, ' levels, got ',  &
                ukca_radaer%d1_nlevs_mmr(i)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
  halo_x = halosize(1, ukca_radaer%d1_halo_type_mmr(i))
  halo_y = halosize(2, ukca_radaer%d1_halo_type_mmr(i))
  
  IF (halo_x /= prev_halo_x .and.                                 &
      halo_y /= prev_halo_y) THEN
    !
    ! Buffer needs to be allocated/reallocated.
    !
    IF (buffer_size /= 0) THEN
      DEALLOCATE(buffer)
    END IF
    
    buffer_size = (row_length + 2*halo_x) *                       &
                  (rows + 2*halo_y) * tr_levs
    ALLOCATE(buffer(1-halo_x:row_length+halo_x,                   &
                    1-halo_y:rows+halo_y,                         &
                    1:tr_levs))
    
    prev_halo_x = halo_x
    prev_halo_y = halo_y
    
  END IF

  IF (ukca_radaer%d1_length_mmr(i) /= buffer_size) THEN
    ierr = 708
    WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,              &
            ' elements, got ',ukca_radaer%d1_length_mmr(i)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  buffer = reshape(D1(ukca_radaer%d1_address_mmr(i) :             &
                      ukca_radaer%d1_address_mmr(i) +             &
                      ukca_radaer%d1_length_mmr(i) - 1),          &
                   (/ row_length+(2*halo_x), rows+(2*halo_y),     &
                      tr_levs /))

  ukca_radaer%mix_ratio(:, :, :, i) =                             &
                buffer(1:row_length, 1:rows, 1:model_levels)

END DO ! i

DO j = 1, ukca_radaer%n_mode

  !
  ! Modal number concentrations.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_nbr(j) /= tr_levs) THEN
    ierr = 707
    WRITE(6,'(2(A,I5))') 'Expecting ', tr_levs, ' levels, got ',  &
                ukca_radaer%d1_nlevs_nbr(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF
  
  halo_x = halosize(1, ukca_radaer%d1_halo_type_nbr(j))
  halo_y = halosize(2, ukca_radaer%d1_halo_type_nbr(j))
  
  IF (halo_x /= prev_halo_x .and.                                 &
      halo_y /= prev_halo_y) THEN
    !
    ! Buffer needs to be allocated/reallocated.
    !
    IF (buffer_size /= 0) THEN
      DEALLOCATE(buffer)
    END IF
    
    buffer_size = (row_length + 2*halo_x) *                       &
                  (rows + 2*halo_y) * tr_levs
    ALLOCATE(buffer(1-halo_x:row_length+halo_x,                   &
                    1-halo_y:rows+halo_y,                         &
                    1:tr_levs))
    
    prev_halo_x = halo_x
    prev_halo_y = halo_y
    
  END IF

  IF (ukca_radaer%d1_length_nbr(j) /= buffer_size) THEN
    ierr = 708
    WRITE(6,'(2(A,I10))') 'Expecting ', buffer_size,              &
         ' elements, got ', ukca_radaer%d1_length_nbr(j)
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) THEN
      CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  buffer = reshape(D1(ukca_radaer%d1_address_nbr(j) :             &
                      ukca_radaer%d1_address_nbr(j) +             &
                      ukca_radaer%d1_length_nbr(j) - 1),          &
                   (/ row_length+(2*halo_x), rows+(2*halo_y),     &
                      tr_levs /))

  ukca_radaer%modal_nbr(:, :, :, j) =                             &
                buffer(1:row_length, 1:rows, 1:model_levels)

END DO ! j

DEALLOCATE(buffer)

IF (lhook) CALL dr_hook('UKCA_RADAER_GET', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_get
