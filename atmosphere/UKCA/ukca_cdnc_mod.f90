! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  Module to define type ukca_cdnc_struct, the structure
!                used by UKCA_CDNC, with intialisation and read routines
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      MODULE UKCA_CDNC_MOD
     
      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      IMPLICIT NONE
      PUBLIC
      SAVE

! Main structure holding all the variables needed 
! to retrieve Cloud Droplet Number Concentration
! from D1, as calculated by UKCA_ACTIVATE

      TYPE UKCA_CDNC_STRUCT
           
        ! Stash code, D1 address, number of levels and total length 
        ! of the D1 fields corresponding to CDNC
        INTEGER :: stashcode_cdnc
        INTEGER :: d1_address_cdnc
        INTEGER :: d1_nlevs_cdnc
        INTEGER :: d1_length_cdnc
       
        ! Stash code, D1 address, number of levels and total length 
        ! of the D1 fields corresponding to CDNC3
        INTEGER :: stashcode_cdnc3
        INTEGER :: d1_address_cdnc3
        INTEGER :: d1_nlevs_cdnc3
        INTEGER :: d1_length_cdnc3
              
        ! Cloud droplet number concentration (m^-3)
        REAL, ALLOCATABLE, DIMENSION(:, :, :) :: cdnc
        ! <Cloud droplet number concentration^-1/3> 
        REAL, ALLOCATABLE, DIMENSION(:, :, :) :: cdnc3
      
      END TYPE UKCA_CDNC_STRUCT

      ! Dimensions of CDNC arrays - will be either (row_length, rows,  
      ! model_levels) or (1,1,1) depending on whether aerosol indirect 
      ! effects are activated or not  
      INTEGER :: cdnc_dim1 
      INTEGER :: cdnc_dim2 
      INTEGER :: cdnc_dim3

! Structure for UKCA CDNC/radiation/ppn interaction
      TYPE (ukca_cdnc_struct) :: ukca_cdnc

      INTERFACE UKCA_CDNC_INIT
        MODULE PROCEDURE UKCA_CDNC_INIT
      END INTERFACE UKCA_CDNC_INIT

      INTERFACE UKCA_CDNC_GET
        MODULE PROCEDURE UKCA_CDNC_GET
      END INTERFACE UKCA_CDNC_GET

      CONTAINS

! ######################################################################

    SUBROUTINE UKCA_CDNC_INIT(ierr,cmessage)
      
    IMPLICIT NONE

    INTEGER, INTENT(INOUT)                  :: ierr        ! Error indicator (0 is OK, >0 error)
    CHARACTER(LEN=256), INTENT(INOUT)       :: cmessage    ! Error message

! Local variables

! The following prognostics are expected in section 34.
    INTEGER, PARAMETER :: stashc_cdnc  = 162
    INTEGER, PARAMETER :: stashc_cdnc3 = 163

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_INIT',zhook_in,    &
                             zhook_handle)

! Initialise STASH code for CDNC prognostics


    ierr = 0
    ukca_cdnc%stashcode_cdnc = stashc_cdnc
    ukca_cdnc%stashcode_cdnc3 = stashc_cdnc3
    IF (PrintStatus >= PrStatus_Diag) THEN 
      WRITE(6,*) 'ukca_cdnc%stashcode_cdnc: ',ukca_cdnc%stashcode_cdnc
      WRITE(6,*) 'ukca_cdnc%stashcode_cdnc3: ',ukca_cdnc%stashcode_cdnc3
    END IF
      
    IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_INIT',zhook_out,   &
                            zhook_handle)
    RETURN
    END SUBROUTINE UKCA_CDNC_INIT

! ######################################################################

      SUBROUTINE UKCA_CDNC_GET(                                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                 first_call,                            &
                                 ierr,                                  &
                                 cmessage)
      
      USE UM_ParVars
      USE atm_fields_bounds_mod,  ONLY: tdims_s
      USE ukca_d1_defs,  ONLY: ukca_sect
      USE Submodel_Mod
      IMPLICIT NONE

      INTEGER, INTENT(OUT)            :: ierr        ! Error indicator (0 is OK, >0 error)
      CHARACTER(len=256), INTENT(OUT) :: cmessage    ! Error message
      LOGICAL, INTENT(IN)             :: first_call  ! Indicates first call:
!                                                    !  complete setup has to be done
      
! D1 and related variables
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

! Local variables
      INTEGER :: i
      INTEGER :: j
      INTEGER :: i_obj     
           
      ! Atmosphere submodel index
      INTEGER :: m_atm_modl  

      ! Logical for checking whether all prognostics have been found
      LOGICAL :: l_missing
       
      ! Expected size
      INTEGER :: buffer_size

      ! Tracer Levels
      INTEGER :: tr_levs
     
      ! Temporary buffer
      REAL, ALLOCATABLE :: temp_buf(:,:,:)
 
      ! Variables for tagged prognostics in D1.
      INTEGER :: section       ! stash section
      INTEGER :: item          ! stash item
      INTEGER :: levs          ! No of levels
      INTEGER :: len           ! length of field
      INTEGER :: addr          ! address in D1
      INTEGER :: halotyp       ! halo type

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',zhook_in,   &
                               zhook_handle)
      
      ierr = 0
      m_atm_modl = SUBMODEL_FOR_SM(a_im)
      tr_levs  = tdims_s%k_end - tdims_s%k_start + 1
      buffer_size = row_length * rows * tr_levs
   
! For first call, check that all needed prognostics are available
! at the expected dimensions. Retain their D1 address for later calls.

      IF (first_call) THEN
        
        ! Initialise d1 address in information structure to "not found" (-1)
        ukca_cdnc%d1_address_cdnc = -1
        ukca_cdnc%d1_address_cdnc3 = -1

        ! Go through D1 and retain some information about the
        ! prognostics we need.
        DO i=1,no_obj_D1(m_atm_modl)
          section = d1_addr(D1_section,i,m_atm_modl)
          IF (section /= ukca_sect) CYCLE
          item    = d1_addr(D1_item,i,m_atm_modl)
          levs    = d1_addr(d1_no_levels,i,m_atm_modl)
          len     = d1_addr(d1_length,i,m_atm_modl)
          addr    = d1_addr(d1_address,i,m_atm_modl)
          halotyp = d1_addr(d1_halo_type,i,m_atm_modl)

          
          IF (item == ukca_cdnc%stashcode_cdnc) THEN
            IF (halotyp /= halo_type_no_halo) THEN
              ierr = 700
              cmessage = 'ukca_cdnc_get: non-zero halo for CDNC'
              IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',    &
                                      zhook_out, zhook_handle)
              RETURN
            END IF
            ukca_cdnc%d1_address_cdnc = addr
            ukca_cdnc%d1_nlevs_cdnc   = levs
            ukca_cdnc%d1_length_cdnc  = len
            IF (PrintStatus >= PrStatus_Diag) THEN 
              WRITE(6,*) 'UKCA_CDNC: address in D1 ',addr
            END IF
          END IF 
               
          IF (item == ukca_cdnc%stashcode_cdnc3) THEN
            IF (halotyp /= halo_type_no_halo) THEN
              ierr = 701
              cmessage = 'ukca_cdnc_get: non-zero halo for CDNC3'
              IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',    &
                                      zhook_out, zhook_handle)
              RETURN
            END IF
            ukca_cdnc%d1_address_cdnc3 = addr
            ukca_cdnc%d1_nlevs_cdnc3   = levs
            ukca_cdnc%d1_length_cdnc3  = len
            IF (PrintStatus >= PrStatus_Diag) THEN 
              write(6,*) 'UKCA_CDNC3: address in D1 ',addr
            END IF
          END IF 

          IF (ukca_cdnc%d1_address_cdnc > 0 .AND.                       & 
              ukca_cdnc%d1_address_cdnc3 > 0) EXIT

        END DO ! i_obj (D1 items)
        
! Check that all prognostics have been found.
        l_missing = .false.
        
        IF (ukca_cdnc%d1_address_cdnc == -1) THEN
           WRITE(6, *) 'Prognostic ', ukca_cdnc%stashcode_cdnc,         &
                        'not found in D1.'
           l_missing = .true.              
        END IF
        
        IF (ukca_cdnc%d1_address_cdnc3 == -1) THEN
           WRITE(6, *) 'Prognostic ', ukca_cdnc%stashcode_cdnc3,        &
                        'not found in D1.'
           l_missing = .true.
        END IF
        
        IF (l_missing) THEN
           ierr = 702
           cmessage =                                                   &
            'ukca_cdnc_get: Prog(s) needed for UKCA are missing from D1'
           IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',       &
                                           zhook_out, zhook_handle)
           RETURN
        END IF 

! Check number of levels and total size of arrays in D1
        IF (ukca_cdnc%d1_nlevs_cdnc /= tr_levs) THEN
           ierr = 705
           WRITE(6,'(A,I4,A,I4)')'Expecting ', tr_levs, ' levels, got ',&
                        ukca_cdnc%d1_nlevs_cdnc
           cmessage =                                                   &
              'ukca_cdnc_get: Unexpected number of levels in D1 prog.'
           IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',       &
                                   zhook_out, zhook_handle)
           RETURN
        END IF

        IF (ukca_cdnc%d1_length_cdnc /= buffer_size) THEN
          ierr = 706
          WRITE(6, *) 'Expecting ', buffer_size, ' elements, got ',     &
                      ukca_cdnc%d1_length_cdnc
          cmessage =                                                    &
            'ukca_cdnc_get: Unexpected total size of D1 prog.'
          IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',        &
                                   zhook_out, zhook_handle)
          RETURN
        END IF
        
        IF (ukca_cdnc%d1_nlevs_cdnc3 /= tr_levs) THEN
          ierr = 705
           WRITE(6,'(A,I4,A,I4)')'Expecting ', tr_levs, ' levels, got ',&
                       ukca_cdnc%d1_nlevs_cdnc3
          cmessage =                                                    &
             'ukca_cdnc_get: Unexpected number of levels in D1 prog.'
          IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',        &
                                  zhook_out, zhook_handle)
          RETURN
        END IF
    
        IF (ukca_cdnc%d1_length_cdnc3 /= buffer_size) THEN
          ierr = 706
          WRITE(6, *) 'Expecting ', buffer_size, ' elements, got ',     &
                      ukca_cdnc%d1_length_cdnc3
          cmessage =                                                    &
            'ukca_cdnc_get: Unexpected total size of D1 prog.'
          IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',        &
                                  zhook_out, zhook_handle)
          RETURN
        END IF
    
      END IF ! first_call

      ALLOCATE(temp_buf(row_length,rows,tr_levs))
      temp_buf(:,:,:) = 0.0

! Cloud droplet number concentration.
      temp_buf(:,:,:)      =  RESHAPE(D1(ukca_cdnc%d1_address_cdnc :    &
                                      ukca_cdnc%d1_address_cdnc +       &
                                      ukca_cdnc%d1_length_cdnc - 1),    &
                                     (/row_length, rows, tr_levs/))
      ukca_cdnc%cdnc(:,:,:) = temp_buf(:,:,1:model_levels)

! <Cloud droplet number concentration^-1/3>
      temp_buf(:,:,:) = 0.0
      temp_buf(:,:,:)      = RESHAPE(D1(ukca_cdnc%d1_address_cdnc3 :    &
                                       ukca_cdnc%d1_address_cdnc3 +     &
                                       ukca_cdnc%d1_length_cdnc3 - 1),  &
                                     (/row_length, rows, tr_levs/))
      ukca_cdnc%cdnc3(:,:,:) = temp_buf(:,:,1:model_levels)

      DEALLOCATE(temp_buf)

      IF (PrintStatus >= PrStatus_Diag) THEN 
        WRITE(6,*) 'CDNC_GET:'
        WRITE(6,*) 'ukca_cdnc%d1_address_cdnc',ukca_cdnc%d1_address_cdnc
        WRITE(6,*) 'ukca_cdnc%d1_address_cdnc3',                        &
                    ukca_cdnc%d1_address_cdnc3
        IF (allocated(ukca_cdnc%cdnc)) WRITE(6,*)                       &
             'ukca_cdnc%cdnc is allocated'
        WRITE(6,*) 'cdnc: ',size(ukca_cdnc%cdnc)
        WRITE(6,*) 'cdnc: ',maxval(ukca_cdnc%cdnc),                     &
                            minval(ukca_cdnc%cdnc)
        WRITE(6,*) 'cdnc3: ',size(ukca_cdnc%cdnc3)
        WRITE(6,*) 'cdnc3: ',maxval(ukca_cdnc%cdnc3),                   &
                             minval(ukca_cdnc%cdnc3)
      END IF


      IF (lhook) CALL dr_hook('UKCA_CDNC_MOD.UKCA_CDNC_GET',zhook_out,  &
                               zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CDNC_GET

END MODULE UKCA_CDNC_MOD
