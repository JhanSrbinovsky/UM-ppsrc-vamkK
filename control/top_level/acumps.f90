! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL    Subroutine:
!LL    ACUMPS
!LL
!LL    Purpose:
!LL    To accumulate partial sums of climate mean tagged diagnostics
!LL    and create dumps containing them. Also to overwrite the D1
!LL    diagnostice with the partial sum for use by MEANPS. This
!LL    saves MEANPS having to reread the partial sum dump.
!LL
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      SUBROUTINE ACUMPS(                                                &
     &  N_OBJS_D1,D1_ADDR                                               &
     &  ,LEN_DATA,D1                                                    &
     &  ,MAXSIZE,MEANS_TOTAL                                            &
     &  ,FLAG,NFTIN,NFTOUT,LCLIMREALYR,MEANLEV                          &
     &  ,I_MONTH,I_YEAR                                                 &
     &  ,HEAD_OUT,HEAD_LEN,HEAD_SIZE,                                   &
     &  TIMESTEP,CMITEMS,FIXHD12,                                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &  ICODE,CMESSAGE)
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE UM_ParVars
      USE lookup_addresses
      USE Submodel_Mod

      IMPLICIT NONE
!
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
      INTEGER                                                           &
     &  N_OBJS_D1               !IN No objects in D1 array

      INTEGER                                                           &
     &  D1_ADDR(D1_LIST_LEN,N_OBJS_D1) !IN Addressing of D1 array

      INTEGER                                                           &
     &  MAXSIZE,                                                        &
                                ! IN dimension of largest data field
     &  LEN_DATA,                                                       &
                                ! IN Length of model data
     &  FLAG,                                                           &
                                ! IN Flag for reading partial sum dump
     &  NFTIN,                                                          &
                                ! IN Unit no for reading partial sums
     &  NFTOUT,                                                         &
                                ! IN Unit no for writing partial sums
     &  ICODE,                                                          &
                                ! OUT Return code; successful=0
                                !                  error>0
     &  MEANLEV,                                                        &
                                ! IN level of climate meaning
     &  MEANS_TOTAL,                                                    &
                                ! IN Indicates a meaning period
     &  I_MONTH,                                                        &
                                ! IN Current model time (months)
     &  I_YEAR,                                                         &
                                ! IN Current model time (years)
     &  FIXHD12,                                                        &
                                ! IN Version of model
     &  CMITEMS,                                                        &
                                ! IN Number of items being meaned
     &  TIMESTEP                ! IN Submodel timestep
!
      CHARACTER(LEN=80)                                                   &
     &       CMESSAGE             ! OUT Error message if ICODE>0
!
     REAL                                                               &
     &  D1(LEN_DATA)            ! IN/OUT Real equivalence of data block
!
     LOGICAL                                                            &
     &  LCLIMREALYR             ! IN Real-period climate meaning
!
!      Common blocks
!
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
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
      INTEGER                                                           &
     &  HEAD_LEN                                                        &
     &  ,HEAD_SIZE

      INTEGER                                                           &
     &  HEAD_OUT(HEAD_LEN,TOTITEMS)                                     &
                                    ! IN Header contains packing
                                !    info for output ps file
     &  ,HEAD_BUF(HEAD_SIZE)

! Header formatted as follows:
! HEAD_OUT(1,*): No of words per level in field
! HEAD_OUT(2,*): 2 for packed, 1 for unpacked
! HEAD_OUT(3,*): No of words per level on disk

!
!      Local variables
!
      INTEGER                                                           &
     &  I,J,K                                                           &
                                ! Loop indices
     &  ,LEN_IO                                                         &
                                ! Actual IO length
     &  ,CITEMS                                                         &
                                ! Count variable
     &  ,PERIODLEN                                                      &
                                ! Current meaning period in days
     &  ,TAG                                                            &
                                ! Stash tag
     &  ,PTD1                                                           &
                                ! Pointer to D1_ADDR information
     &  ,address                                                        &
                                ! Address in local D1
     &  ,levels                                                         &
                                ! Number of levels per diagnostic
     &  ,length                                                         &
                                ! Length of each level in local D1
     &  ,global_length                                                  &
                                ! Length of global field
     &  ,offset                 ! Indexing offset for WORK array
!
      INTEGER                                                           &
     &  HEADER(2)                                                       &
                                ! Initial header
     &  ,HEAD_IN(HEAD_LEN,CMITEMS) ! Packing info for input ps file
                                ! Will differ from HEAD_OUT if packing
                                ! codes have changed mid-run

      REAL                                                              &
     &  IOSTAT,                                                         &
                                ! IO error code
     &  REALPERIODLEN,                                                  &
                                ! explicitly real equivalent
                                ! of PERIODLEN
     &  AWORK,                                                          &
                                ! Accumulative sum for WORK array
     &  CWORK,                                                          &
                                ! Accumulative sum for WORK array
     &  CKWORK,                                                         &
                                ! Checksum for WORK array
     &  CKWORKO                 ! Packed CKWORK
!
!      Local arrays
!
      REAL                                                              &
     &  D1_DATA(MAXSIZE)        ! Work area for fields
      REAL                                                              &
     &  WORK(MAXSIZE+4)        ! Work area and IO buffer

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Align for well-formed io
!dir$ cache_align work

!
      IF (lhook) CALL dr_hook('ACUMPS',zhook_in,zhook_handle)
      IF (ICODE /= 0) GOTO 9999

! Arrays sent to BUFFIN/BUFFOUT need to be cache aligned for
! well-formed io to work, but the cache_align directive wasn't working
! correctly for WORK. The following adds an offset to the index of
! WORK which resolves the problem.
! This is an offset from 0, so an offset of 1 really means no offset!
      offset=1

!   Set up variables needed for weighting accumulations if real-period
!   climate meaning is selected. Partial sums are normalised elsewhere.

      if (lclimrealyr) then
! DEPENDS ON: setperlen
        call setperlen(meanlev,i_month,i_year,periodlen)
        realperiodlen=real(periodlen)
      endif

! STEP 1: Read in headers of previous partial sum and write out
!         header of new.

      IF (FLAG /= 1) THEN       ! PS data exist on disk
! Read headers for input partial sum file

        CALL BUFFIN(NFTIN,HEAD_BUF,HEAD_SIZE,LEN_IO,IOSTAT)
        IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_SIZE)THEN
          WRITE(6,*)'ACUMPS: Error reading header: IO code ',           &
     &      IOSTAT,' on unit ',NFTIN
          WRITE(6,*)'Words requested ',HEAD_SIZE,                       &
     &      ' Words read ',LEN_IO
          ICODE=1
          CMESSAGE='ACUMPS: BUFFIN error - see output'
          GOTO 9999
        ENDIF
! Transfer header information from buffer to header arrays
        HEADER(1)=HEAD_BUF(1) ! Timestep of creation
        HEADER(2)=HEAD_BUF(2) ! Number of records
        K=3
        DO I=1,CMITEMS
          DO J=1,HEAD_LEN
            HEAD_IN(J,I)=HEAD_BUF(K)
            K=K+1
          ENDDO
        ENDDO

        IF (HEADER(1) >= TIMESTEP.OR.HEADER(2) /= CMITEMS)THEN
          WRITE(6,*)'ACUMPS1: Partial sum file inconsistent'
          WRITE(6,*)'PS file holds ', &
              HEADER(2),' items and written at STEP ',HEADER(1)
          WRITE(6,*)'Expected timestep should be < ',TIMESTEP
          WRITE(6,*)'Expected number of items ',CMITEMS
          CMESSAGE='ACUMPS1: Partial sum file inconsistent. See Output'
          ICODE=2
          GOTO 9999
        ENDIF
      ELSE
! No input sum, so initialise header array
        DO I=1,HEAD_SIZE
          HEAD_BUF(I)=0
        ENDDO
      ENDIF

! Write headers for new partial sum file
! Transfer information to io buffer
      HEAD_BUF(1)=TIMESTEP
      HEAD_BUF(2)=CMITEMS
      K=3
      DO I=1,CMITEMS
        DO J=1,HEAD_LEN
          HEAD_BUF(K)=HEAD_OUT(J,I)
          K=K+1
        ENDDO
      ENDDO

      CALL BUFFOUT(NFTOUT,HEAD_BUF,HEAD_SIZE,LEN_IO,IOSTAT)
      IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_SIZE)THEN
        WRITE(6,*)'ACUMPS: Error writing header: IO code ',             &
     &    IOSTAT,' on unit ',NFTOUT
        WRITE(6,*)'Words requested ',HEAD_SIZE,                         &
     &    ' Words written ',LEN_IO
        ICODE=4
        CMESSAGE='ACUMPS: BUFFOUT error - see output'
        GOTO 9999
      ENDIF

! STEP 2 : Loop over all STASH items. For each tagged item, gather
!          current data to D1_DATA array, read partial sum into WORK
!          array (if there is a partial sum), sum the two and write
!          out to new partial sum file.
!           Also, if this is a meaning period, overwrite the field
!          in D1 with the complete sum, to be picked up by MEANPS.

!     Start of loop over STASH items
      CITEMS=0
      DO K=1,TOTITEMS
        TAG=STLIST(st_macrotag,K)/1000
        PTD1=STLIST(st_d1pos,K)
        IF(TAG/=0)THEN
          IF(STLIST(s_modl,k)==D1_ADDR(d1_imodl,PTD1))THEN
! Object tagged for climate meaning and in relevant internal model
          address=D1_ADDR(d1_address,PTD1)
          levels=D1_ADDR(d1_no_levels,PTD1)
          length=D1_ADDR(d1_length,PTD1)/levels
          global_length=STLIST(st_dump_level_output_length,K)
          CITEMS=CITEMS+1
          DO J=1,levels
! Copy current field from D1 to D1_DATA
! by gathering full field to pe0
! DEPENDS ON: general_gather_field
            CALL GENERAL_GATHER_FIELD(                                  &
     &        D1(address),D1_DATA,length,                               &
     &        global_length,1,                                          &
     &        D1_ADDR(1,PTD1),0,                                        &
     &        -1,ICODE,CMESSAGE)
            IF (ICODE /= 0) GOTO 9999
            DO I=global_length+1,MAXSIZE
              D1_DATA(I)=0.
            ENDDO
! Set initial value for AWORK and CWORK
            AWORK=0.0
            CWORK=0.0
! If partial sum exists on disk, read it in and add to current field
            IF (FLAG /= 1) THEN ! PS data exist on disk
! Read in one level of partial sum field
              IF (HEAD_IN(2,CITEMS)  ==  2) THEN

! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*HEAD_IN(3,CITEMS) 32 bit words using BUFFIN32_f77 (because the 
! array is 64 bit)

! DEPENDS ON : buffin32_f77
                 CALL BUFFIN32_f77(NFTIN,WORK(offset:), &
                      2*HEAD_IN(3,CITEMS),LEN_IO,IOSTAT)

! And then halve LEN_IO to satisfy tests against HEAD_IN(3,CITEMS)
                LEN_IO = LEN_IO/2
              ELSE ! For non-packed data
                 CALL BUFFIN(NFTIN,WORK(offset:),HEAD_IN(3,CITEMS)        &
     &            ,LEN_IO,IOSTAT)
              ENDIF
              IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_IN(3,CITEMS))THEN
                WRITE(6,*)'ACUMPS: Error reading partial sum IO code ', &
     &            IOSTAT,' on unit ',NFTIN
                WRITE(6,*)'Words requested ',HEAD_IN(3,CITEMS),         &
     &            ' Words read ',LEN_IO
                ICODE=6
                CMESSAGE='ACUMPS: BUFFIN error - see output'
                GOTO 9999
              ENDIF
              IF (mype == 0) THEN
! Valid data exists on pe0 only
! Unpack if data on disk was packed
                IF (HEAD_IN(2,CITEMS) == 2)THEN
! DEPENDS ON: expand32b
                  CALL EXPAND32B(GLOBAL_LENGTH+1,WORK(offset),FIXHD12)
! Calculate a checksum
                  DO I=1,GLOBAL_LENGTH
                    AWORK=AWORK+WORK(I+offset-1)
                  END DO
                  CKWORK=AWORK/INT(global_length)
! Pack and umpack checksum to force it losing precision in order to do
! the comparison
! DEPENDS ON: pack21
                  CALL PACK21(1,CKWORK,CKWORKO)
! DEPENDS ON: expand32b
                  CALL EXPAND32B(1,CKWORKO,FIXHD12)
                ELSE
                  DO I=1,GLOBAL_LENGTH
                    AWORK=AWORK+WORK(I+offset-1)
                  END DO
                  CKWORKO=AWORK/INT(global_length)
                ENDIF
                IF(CKWORKO /= WORK(global_length+offset)) THEN
                  WRITE(6,*)'ERROR: checksum failure in climate mean'
                  WRITE(6,*)'Section ',D1_ADDR(d1_section,PTD1), &
                     ' item ',D1_ADDR(d1_item,PTD1)
                  WRITE(6,*) 'This can be due to invalid values in field, ', &
                     'or corruption of partial sum file'
                  WRITE(6,*)'Remove or fix diagnostic, and rerun'
                  ICODE=4
                  CMESSAGE='ACUMPS: Diagnostic error. See output for item no.'
                  GOTO 9999
                ENDIF
! Sum with field in D1 - Scale data if 365 day calendar
                IF (LCLIMREALYR)THEN
                  DO I=1,global_length
                    IF (WORK(I+offset-1) == RMDI)THEN
                      D1_DATA(I)=RMDI
                    ELSE
                      D1_DATA(I)=WORK(I+offset-1)+                      &
     &                  (realperiodlen*D1_DATA(I))
                    ENDIF
                  END DO
                ELSE
! 360 day calendar
                  DO I=1,global_length
                    IF (WORK(I+offset-1) == RMDI)THEN
                      D1_DATA(I)=RMDI
                    ELSE
                      D1_DATA(I)=WORK(I+offset-1)+D1_DATA(I)
                    ENDIF
                  END DO
                ENDIF
              ENDIF
            ELSE
! First data for this period - no partial sum to add
              IF (LCLIMREALYR)THEN
! Scale initial data if 365 day calendar
                IF (mype == 0) THEN
                  DO I=1,global_length
                    IF (D1_DATA(I) /= RMDI)THEN
                      D1_DATA(I)=realperiodlen*D1_DATA(I)
                    ENDIF
                  END DO
                ENDIF
              ENDIF
            ENDIF ! End of adding PS data

!         Write out sum to PS file
            IF (mype == 0) THEN
! Copy data to WORK array, packing if necessary
              IF (HEAD_OUT(2,CITEMS) == 2)THEN
                DO I=HEAD_OUT(1,CITEMS),MAXSIZE
                  WORK(I+offset-1)=0.
                ENDDO
! DEPENDS ON: pack21
                CALL PACK21(GLOBAL_LENGTH+1,D1_DATA,                    &
                            WORK(offset))
! DEPENDS ON: expand32b
                CALL EXPAND32B(GLOBAL_LENGTH+1,WORK(offset),FIXHD12)
                DO I=1,GLOBAL_LENGTH
                  CWORK=CWORK+WORK(I+offset-1)
                END DO
                WORK(global_length+offset)=CWORK/INT(global_length)
! DEPENDS ON: pack21
                CALL PACK21(GLOBAL_LENGTH+1,WORK(offset),               &
                            WORK(offset) )

! If data not packed, calculate checksum straight away
              ELSE
                DO I=1,GLOBAL_LENGTH
                  WORK(I+offset-1)=D1_DATA(I)
                  CWORK=CWORK+WORK(I+offset-1)
                ENDDO
                WORK(global_length+offset)=CWORK/INT(global_length)
              ENDIF
            ENDIF


! Output partial sum to file
            IF (HEAD_OUT(2,CITEMS)  ==  2) THEN

! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*HEAD_OUT(3,CITEMS) 32 bit words using BUFFOUT32_F77 (because
! we have a supplied 64 bit array)

! DEPENDS ON : buffout32_f77
               CALL BUFFOUT32_f77(NFTOUT,WORK(offset:),                  &
                    2*HEAD_OUT(3,CITEMS),LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against HEAD_OUT(3,CITEMS)
               LEN_IO = LEN_IO/2
            ELSE
! For non-packed data

              CALL BUFFOUT(NFTOUT,WORK(offset:),                        &
     &          HEAD_OUT(3,CITEMS),LEN_IO,IOSTAT)
            ENDIF
            IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_OUT(3,CITEMS))THEN
              WRITE(6,*)'ACUMPS: Error writing partial sum. Code ',     &
     &          IOSTAT,' on unit ',NFTOUT
              WRITE(6,*)'Words requested ',HEAD_OUT(3,CITEMS),          &
     &          ' Words written ',LEN_IO
              ICODE=7
              CMESSAGE='ACUMPS: BUFFOUT error - see output'
              GOTO 9999
            ENDIF
            IF (MEANS_TOTAL /= 0)THEN
! Overwrite field in D1 with partial sum for use by MEANPS
              IF (mype == 0)then
! Pack and unpack for bit comparison with old system
                IF (HEAD_OUT(2,CITEMS) == 2)THEN
                  DO I=1,HEAD_OUT(1,CITEMS)
                    D1_DATA(I)=WORK(I+offset-1)
                  ENDDO
! DEPENDS ON: expand32b
                  CALL EXPAND32B(GLOBAL_LENGTH,D1_DATA,FIXHD12)
                ENDIF
              ENDIF
! DEPENDS ON: general_scatter_field
              CALL GENERAL_SCATTER_FIELD(                               &
     &         D1(address),D1_DATA,LENGTH,global_length,1,              &
     &          D1_ADDR(1,PTD1),0,ICODE,CMESSAGE)
              IF (ICODE /= 0) GOTO 9999
            ENDIF
            address=address + length ! Point to next level
          ENDDO                 ! End loop over levels
          ENDIF
        ENDIF                   ! End tagged for meaning
      END DO                    ! End of loop over STASH list

 9999 CONTINUE
      IF (lhook) CALL dr_hook('ACUMPS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ACUMPS
