! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Description:
!  The code `nudges' model variables towards meteorological analyses
!  in order to reproduce `real weather'.
!  Top level routine for Nudging. Calls routines to read in ECMWF data
!  fields and nudge model variables (T,U,V) towards ECMWF.
!  Diagnostics related to nudging are written into STASHwork to be updated
!  to STASH by ATM_STEP

!  The nudged model was developed as part of the UKCA project.
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

!  Available for non-commerical use at present. Please cite
!
!  Telford, P. J., Braesicke, P., Morgenstern, O., and Pyle, J. A.:
!  Technical Note: Description and assessment of a nudged version of
!  the new dynamics Unified Model, Atmos. Chem. Phys., 8, 1701-1712, 2008.

!  Method:
! 1) Interface with atmosphere model:
!    Prognostic fields are passed in from Atm_Step. The nudged variables
!    are passed back to Atm_Step as arguments.
! 2) Sets up various quantities (eg filenames of analyses )
!    before calling routines to nudge T/theta, U & V
! 3) The diagnostics from the nudging are put into STASH using
!    copy_diag() and stashwork array.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! CONTAINED subroutines:
!     NUDGING_NETCDF_LOADER
!     NUDGING_ALLOC_DATA
!     NUDGING_GETD1_DATA
!     NUDGING_UPDATE_STASH
!     NUDGING_DEALLOC_DATA

!  Code Description:
!    Language:  FORTRAN 90 (formatted)

! ######################################################################

! Subroutine Interface:
SUBROUTINE nudging_main1(                                        &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end

! in time stepping information.
  timestep, i_year, i_month, i_day, i_hour, i_minute,            &
  i_second, timestep_number,                                     &

! in data fields.
  theta, u, v, p_rho_levels, exner_theta_levels, p_theta_levels, &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end

  st_work)

USE nudging_control              ! Include control module for nudging
USE nudging_d1_defs              ! Use own way of accessing D1 array
USE atm_fields_bounds_mod, ONLY:  udims, vdims, tdims, pdims,&
      udims_s, vdims_s, tdims_s, pdims_s

USE trignometric_mod             ! Include latitude and longitude information
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE Field_Types
USE UM_ParVars
USE control_max_sizes
USE Submodel_Mod
IMPLICIT NONE

!***********************************************************

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

!*************************************************************************
! arguments with intent in. ie: input variables.

! Parameters

! Submodel indices
INTEGER  m_atm_modl, im_index

! model parameters
REAL                                                              &
  timestep
! time information for current timestep
INTEGER                                                           &
  i_year, i_month, i_day                                          &
, i_hour, i_minute, i_second                                      &
, timestep_number

! Diagnostics info
REAL  st_work(*)             ! STASH workspace for section

! Data arrays
REAL, INTENT(INOUT) ::                                 &
!
 u( udims_s%i_start : udims_s%i_end,                   &
    udims_s%j_start : udims_s%j_end,                   &
    udims_s%k_start : udims_s%k_end ),                 &
!
 v( vdims_s%i_start : vdims_s%i_end,                   &
   vdims_s%j_start : vdims_s%j_end,                    &
   vdims_s%k_start : vdims_s%k_end ),                  &
!
 theta( tdims_s%i_start : tdims_s%i_end,               &
        tdims_s%j_start : tdims_s%j_end,               &
        1:model_levels )

REAL, INTENT(IN) ::                                    &
!
 p_rho_levels( pdims_s%i_start : pdims_s%i_end,        &
               pdims_s%j_start : pdims_s%j_end,        &
               pdims_s%k_start : pdims_s%k_end + 1 ),  &
!
 p_theta_levels( tdims_s%i_start : tdims_s%i_end,      &
                 tdims_s%j_start : tdims_s%j_end,      &
                 tdims_s%k_start : tdims_s%k_end ),    &
!
 exner_theta_levels( tdims_s%i_start : tdims_s%i_end,  &
                     tdims_s%j_start : tdims_s%j_end,  &
                     tdims_s%k_start : tdims_s%k_end )

!**************************************************************************
! Begin local Header

! Miscellaneous variables
LOGICAL, SAVE  :: l_first = .TRUE.        ! Is this the first time
INTEGER        :: a_steps_per_hr          ! no. of atm. steps per hour
INTEGER        :: i,j,k,ii,jj,kk,jnext     ! Loop variables

!***********************************************************************************
! Data_type information from D1
INTEGER                   :: grid_type_theta,field_type_theta
INTEGER                   :: grid_type_u,field_type_u
INTEGER                   :: grid_type_v,field_type_v
! Derived variables
REAL, ALLOCATABLE  :: p_ugrid_rho_levels(:,:,:)  ! pressure on u grid
REAL, ALLOCATABLE  :: p_vgrid_rho_levels(:,:,:)  ! pressure on v grid

! Tropopause pressure (and derivations) 
REAL, ALLOCATABLE  :: trop_pressure(:,:)       ! tropopause pressure on T grid 
REAL, ALLOCATABLE  :: trop_pressure_ugrid(:,:) ! tropopause pressure on u grid 
REAL, ALLOCATABLE  :: trop_pressure_vgrid(:,:) ! tropopause pressure on v grid 

! Variables associated with the above
INTEGER   :: theta_row_length_min   ! Min theta column (no halo)
INTEGER   :: theta_row_length_max   ! Max theta column 
INTEGER   :: theta_rows_min         ! Min theta row (no halo)
INTEGER   :: theta_rows_max         ! Max theta row 

INTEGER   :: u_row_length_min       ! Min U column 
INTEGER   :: u_row_length_max       ! Max U column 
INTEGER   :: u_rows_min             ! Min U row 
INTEGER   :: u_rows_max             ! Max U row 

INTEGER   :: v_row_length_min       ! Min V column 
INTEGER   :: v_row_length_max       ! Max V column 
INTEGER   :: v_rows_min             ! Min V row 
INTEGER   :: v_rows_max             ! Max V row 

REAL      :: frac                   ! Frac between ECMWF timesteps
INTEGER   :: proc_row_length_min    ! min column in PE
INTEGER   :: proc_trow_length_max   ! max column T in PE
INTEGER   :: proc_urow_length_max   ! max column U in PE
INTEGER   :: proc_vrow_length_max   ! max column V in PE
INTEGER   :: proc_rows_min          ! min row in PE
INTEGER   :: proc_trows_max         ! max row in PE
INTEGER   :: proc_urows_max         ! max row in PE
INTEGER   :: proc_vrows_max         ! max row in PE

!************************************************************
! Diagnostics

REAL, ALLOCATABLE   :: diag_theta_data   (:,:,:)     ! T diag 1: ECMWf T on mod levs
REAL, ALLOCATABLE   :: diag_theta_model  (:,:,:)     ! T diag 2: model T on mod levs
REAL, ALLOCATABLE   :: diag_theta_ntend  (:,:,:)     ! T diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_theta_mtend  (:,:,:)     ! T diag 4: model tend tendency
REAL, ALLOCATABLE   :: diag_theta_relax  (:,:,:)     ! T diag 5: relax par
REAL, ALLOCATABLE   :: diagsq_theta_ntend(:,:,:)     ! T diag 6: nudging tend^2
REAL, ALLOCATABLE   :: diagsq_theta_mtend(:,:,:)     ! T diag 7: model tend^2

REAL, ALLOCATABLE   :: diag_u_data       (:,:,:)     ! u diag 1: ECMWf U on mod levs
REAL, ALLOCATABLE   :: diag_u_model      (:,:,:)     ! u diag 2: model U on mod levs
REAL, ALLOCATABLE   :: diag_u_ntend      (:,:,:)     ! u diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_u_mtend      (:,:,:)     ! u diag 4: model tendency
REAL, ALLOCATABLE   :: diag_u_relax      (:,:,:)     ! u diag 5: relax par
REAL, ALLOCATABLE   :: diagsq_u_ntend    (:,:,:)     ! u diag 6: nudging tend^2
REAL, ALLOCATABLE   :: diagsq_u_mtend    (:,:,:)     ! u diag 7: model tend^2

REAL, ALLOCATABLE   :: diag_v_data       (:,:,:)     ! v diag 1: ECMWf V on mod levs
REAL, ALLOCATABLE   :: diag_v_model      (:,:,:)     ! v diag 2: model V on mod levs
REAL, ALLOCATABLE   :: diag_v_ntend      (:,:,:)     ! v diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_v_mtend      (:,:,:)     ! v diag 4: model tendency
REAL, ALLOCATABLE   :: diag_v_relax      (:,:,:)     ! v diag 5: relax par

!**********************************************************************
! Variables associated with the netcdf files

REAL, SAVE, ALLOCATABLE   :: temp_file_data        (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: temp_file_data2       (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: u_file_data           (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: u_file_data2          (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: v_file_data           (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: v_file_data2          (:,:,:,:)

REAL, SAVE, ALLOCATABLE   :: surf_logp_file_data   (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logp_file_data2  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpu_file_data  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpu_file_data2 (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpv_file_data  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpv_file_data2 (:,:,:,:)

!     Dimensions of Netcdf Grid
INTEGER                   :: file_dims(totdims)   ! Array of all dim lengths
INTEGER                   :: file_dims2(totdims)
INTEGER, SAVE             :: file_row_length      ! no. of (T) pts in a row
INTEGER, SAVE             :: file_rows            ! no of (T) rows
INTEGER, SAVE             :: file_levels          ! no. of model levels
INTEGER, SAVE             :: file_timesteps       ! no. of timesteps
INTEGER, SAVE             :: file_urow_length     ! no. of (U) pts in a row
INTEGER, SAVE             :: file_vrows           ! no of (V) rows
INTEGER, SAVE             :: timestep1            ! Timestep of first data point
INTEGER, SAVE             :: timestep2            ! Timestep of second data poin
CHARACTER(LEN=256)        :: dataname1            ! first file+pathname
CHARACTER(LEN=256)        :: dataname2            ! second file+pathname

INTEGER                   :: icode                ! return code

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_MAIN',zhook_in,zhook_handle)

!*************************************************************************
! End of the Header

! In original code debugging flag was a parameter set in the nudging
! control module. Change this so set based on UMs UMUI based debug flag.
! Set once here to minimise changes though if were to rewrite from scratch
! would base on the PrintStatus method.
IF(l_first) THEN
  SELECT CASE(printstatus)

  CASE(prstatus_min)
    glob_debug = 0
  CASE(prstatus_normal)
    glob_debug = 5
  CASE(prstatus_oper)
    glob_debug = 10
  CASE(prstatus_diag)
    glob_debug = 15
  END SELECT
END IF

! Standard subroutine entry comment
! Standard Entry Comment
IF(glob_debug > 0) THEN
  WRITE(OUT,'(A)')'NUDGING_MAIN: Entering routine '
END IF

!**********************************************************

!     Calculate various miscellanea before we start

!**********************************************************

! Confirm that the namelist Run_nudging has been specified properly & read
IF ( ndg_relax_uvalue == 0. .OR. ndg_relax_vvalue == 0. .OR.           &
    ndg_relax_tvalue == 0. .OR. ndg_lev_bottom   <  1  .OR.            &
    ndg_lev_top < 1 .OR. ndg_on_lev_bottom < 1         .OR.            &
    ndg_on_lev_top < 1 .OR.  ndg_datapath == ' '       .OR.            &
    ndg_lev_top > model_levels .OR. ndg_lev_bottom > model_levels .OR. &
    ndg_strat_fac < 0. .OR. ndg_strat_fac > 1.         .OR.            &
    ndg_hours_perdata < 1 .OR. ndg_analysis_source > 3 )   THEN

  icode = 999

  CALL ereport ('NUDGING_MAIN',icode,                                  &
     'Improper values found in namelist: RUN_Nudging')
END IF

! Standard Entry Comment
IF(glob_debug > 15) THEN
  WRITE(OUT,'(A)') 'NUDGING_MAIN: Beginning preparations '
END IF

!************************************************
! Load the submodel code
m_atm_modl   = submodel_for_sm(a_im)
im_index     = internal_model_index(atmos_im)

! Load the number of dynamic timesteps per hour
a_steps_per_hr  = INT(3600./timestep)

!*******************************************************
! Load where the local (processor) indices are in the global array
! using Datastart and Blsize from ParVars
proc_row_length_min  = datastart(1)
proc_trow_length_max = proc_row_length_min + blsize(1,fld_type_p) -1
proc_urow_length_max = proc_row_length_min + blsize(1,fld_type_u) -1
proc_vrow_length_max = proc_row_length_min + blsize(1,fld_type_v) -1

proc_rows_min  = datastart(2)
proc_trows_max = proc_rows_min + blsize(2,fld_type_p) -1
proc_urows_max = proc_rows_min + blsize(2,fld_type_u) -1
proc_vrows_max = proc_rows_min + blsize(2,fld_type_v) -1

! If in debug mode write out where we are on the grid
IF(glob_debug > 15) THEN
  WRITE(OUT,'(A,8(1X,I4))' )                                          &
   ' NUDGING_MAIN: Minimum and Maximum Indices are',                   &
   proc_row_length_min, proc_trow_length_max, proc_urow_length_max,    &
   proc_vrow_length_max, proc_rows_min, proc_trows_max, proc_urows_max,&
   proc_vrows_max
END IF

! Define location for debug prints --approx centre of local domain
dbg_x = (proc_trow_length_max - proc_row_length_min)/2
dbg_y = (proc_trows_max - proc_rows_min)/2
dbg_z = (ndg_lev_top - ndg_lev_bottom)/2

! Define the size of each grid type we are using (no halo).

u_row_length_min = udims%i_start
u_row_length_max = udims%i_end 
u_rows_min       = udims%j_start
u_rows_max       = udims%j_end 

v_row_length_min = vdims%i_start
v_row_length_max = vdims%i_end 
v_rows_min       = vdims%j_start
v_rows_max       = vdims%j_end 

theta_row_length_min = tdims%i_start
theta_row_length_max = tdims%i_end 
theta_rows_min       = tdims%j_start
theta_rows_max       = tdims%j_end 

! If in debug mode provide information about dimensions
IF(glob_debug > 15) THEN
  WRITE(OUT,'(A,4(1X,I4))')                                            &
    ': NUDGING_MAIN: Local u array bounds = ',                         &
    u_row_length_min,  u_row_length_max,                               &
    u_rows_min,  u_rows_max

   WRITE(OUT,'(A,4(1X,I4))')                                           &
    ': NUDGING_MAIN: Local v array bounds = ',                         &
    v_row_length_min,  v_row_length_max,                               &
    v_rows_min,  v_rows_max

  WRITE(OUT, '(A,4(1X,I4))')                                           &
    ': NUDGING_MAIN: Local theta array bounds = ',                     &
    theta_row_length_min,  theta_row_length_max,                       &
    theta_rows_min,  theta_rows_max

END IF ! debug output

!*******************************************************

! Allocate required variables - diagnostics & local arrays
CALL nudging_alloc_data

! Load required diagnostic & u,v,theta grid information
! Need to access D1 for this
CALL nudging_getd1_data

IF(glob_debug > 15) THEN
  WRITE(OUT,'(A)') 'NUDGING_MAIN: Downloaded Diagnostics'
END IF


!********************************************************************
! Calculate pressure on U and V grids on rho levels

p_ugrid_rho_levels(:,:,:) = 0.
p_vgrid_rho_levels(:,:,:) = 0.


! pressure on u grid
j = pdims%j_start
DO jj = u_rows_min, u_rows_max
  i = pdims%i_start
  DO ii = u_row_length_min, u_row_length_max

! DEPENDS ON: nudging_interpolation_zero_d  
    CALL nudging_interpolation_zero_d(    &   
          trop_pressure(i,j),             & ! starting variable
          trop_pressure(i+1,j),           & ! finishing variable
          0.5,                            & ! always at mid-point
          trop_pressure_ugrid(ii,jj),     & ! interpolated variable
          0,                              & ! Linear interpolation
          no_debug)                         ! Debug flag  

    DO k=1, model_levels
! DEPENDS ON: nudging_interpolation_zero_d
       CALL nudging_interpolation_zero_d(   &
             p_rho_levels(i,j,k),           & ! starting variable
             p_rho_levels(i+1,j,k),         & ! finishing variable
             0.5,                           & ! always at mid-point
             p_ugrid_rho_levels(ii,jj,k),   & ! interpolated variable
             0,                             & ! Linear interpolation
             no_debug)                        ! Debug Flag
    END DO   ! levels
    i = i + 1
  END DO     ! row_length
  j = j + 1
END DO       ! rows

! pressure on v grid
j = pdims%j_start
DO jj = v_rows_min, v_rows_max

! For EG/V-at-poles: last V point outside P-grid
! Copy/ interpolate between last P-row itself
  jnext = MIN( j+1, pdims_s%j_end )
  i = pdims%i_start
  DO ii = v_row_length_min, v_row_length_max

! DEPENDS ON: nudging_interpolation_zero_d
    CALL nudging_interpolation_zero_d(     &  
          trop_pressure(i,j),              &  ! starting variable 
          trop_pressure(i,jnext),          &  ! finishing variable  
          0.5,                             &  ! always at mid-point  
          trop_pressure_vgrid(ii,jj),      &  ! interpolated variable  
          0,                               &  ! Linear interpolation  
          no_debug)                           ! Debug flag

    DO k=1, model_levels
! DEPENDS ON: nudging_interpolation_zero_d
      CALL nudging_interpolation_zero_d(   &
              p_rho_levels(i,j,k),         &  ! starting variable
              p_rho_levels(i,jnext, k),    &  ! finishing variable
              0.5,                         &  ! always at mid-point
              p_vgrid_rho_levels(ii,jj,k), &  ! interpolated variable
              0,                           &  ! Linear interpolation
              no_debug)                       ! Debug flag
    END DO     ! levels
    i = i + 1
  END DO       ! row_length
  j = j + 1
END DO         ! rows

! Drop a line letting us know where we are
IF(glob_debug > 10) THEN
  WRITE(OUT,'(A)')' NUDGING_MAIN: Loaded variables from D1 array'
END IF

!****************************************************************
! Load netcdf files and obtain basic information

! Get fraction between the two ECMWF timesteps of the model timestep
! DEPENDS ON: nudging_getfrac
CALL nudging_getfrac(                 &
 a_steps_per_hr,                      & ! No. of model timesteps per hour
 ndg_hours_perdata,                   & ! No. of hours per data timestep
 i_hour,                              & ! Model hour
 i_minute,                            & ! Model minute
 frac,                                & ! Fraction between data timesteps
 no_debug)                              ! Debug Flag

! If in debug mode then write out the value of frac
IF(glob_debug > 15) THEN
  WRITE(OUT,'(A,F12.4)')                                           &
   'NUDGING_MAIN: Fraction between data steps equals ', frac
END IF

! If we have started a new data timestep then reload data
! NB this assumes that the data timestep is made up of an
! integral no. of atm steps
IF(frac == 0.OR.l_first) CALL nudging_netcdf_loader(glob_debug)

! If in debug mode drop a line letting us know where we are
IF(glob_debug > 10) THEN
  WRITE(OUT,'(A)')'NUDGING_MAIN: Loaded data from netcdf files '
END IF

!**********************************************************
! Call the routine that controls the nudging

! DEPENDS ON: nudging_nudging_control
CALL nudging_nudging_control(             &
  temp_name,                              &
  grid_type_theta,                        & ! Grid type
  field_type_theta,                       & ! Field type
  ana_row_length,                         & ! Global row length
  ana_rows,                               & ! Global rows
  file_levels,                            & ! Netcdf file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_trow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_trows_max,                         & ! Maximum row
  timestep,                               & ! model timestep
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  base_lambda,                            &
  delta_lambda,                           &
  base_phi,                               &
  delta_phi,                              &
  ana_base_lambda,                        &
  ana_delta_lambda,                       &
  ana_base_phi,                           &
  ana_delta_phi,                          &
  temp_file_data       (:,:,:,timestep1), &
  temp_file_data2      (:,:,:,timestep2), &
  surf_logp_file_data  (:,:,1,timestep1), &
  surf_logp_file_data2 (:,:,1,timestep2), &
  frac,                                   & ! Fraction between data
  model_levels,                           & ! model levels
  p_theta_levels (                        & ! model pressure
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max,          &
    1:model_levels ),                     &
  trop_pressure(                          & ! Tropopause pressure
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max),         &
  theta (                                 & ! nudged variable (theta)
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max,          &
    1:model_levels ),                     &
  diag_theta_data,                        & ! theta diag 1
  diag_theta_model,                       & ! theta diag 2
  diag_theta_ntend,                       & ! theta diag 3
  diag_theta_mtend,                       & ! theta diag 4
  diag_theta_relax,                       & ! theta diag 5
  glob_debug                              & ! Debug flag
  )

!********************************************************************
! Zonal Wind
! DEPENDS ON: nudging_nudging_control
CALL nudging_nudging_control(             &
  u_name,                                 &
  grid_type_u,                            & ! Grid type
  field_type_u,                           & ! Field type
  ana_row_length,                         & ! Global row length
  ana_rows,                               & ! Global rows
  file_levels,                            & ! Netcfd file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_urow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_urows_max,                         & ! Maximum row
  timestep,                               & ! model timestep
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  (base_lambda+0.5*delta_lambda),         &
  delta_lambda,                           &
  base_phi,                               &
  delta_phi,                              &
  ana_base_lambda_u,                      &
  ana_delta_lambda,                       &
  ana_base_phi,                           &
  ana_delta_phi,                          &
  u_file_data(:,:,:,timestep1),           &
  u_file_data2(:,:,:,timestep2),          &
  surf_logpu_file_data(:,:,1,timestep1),  &
  surf_logpu_file_data2(:,:,1,timestep2), &
  frac,                                   & ! fraction between data
  model_levels,                           & ! model levels
  p_ugrid_rho_levels                      & ! model pressure
   ( u_row_length_min:u_row_length_max,   &
     u_rows_min:u_rows_max,               &
     1:model_levels ),                    &
  trop_pressure_ugrid(                    & ! Tropopause pressure
     u_row_length_min:u_row_length_max,   &
     u_rows_min:u_rows_max),              &
  u (                                     & ! nudged variable (u)
   u_row_length_min:u_row_length_max,     &
   u_rows_min:u_rows_max,                 &
   1:model_levels ),                      &
  diag_u_data,                            & ! u diag 1
  diag_u_model,                           & ! u diag 2
  diag_u_ntend,                           & ! u diag 3
  diag_u_mtend,                           & ! u diag 4
  diag_u_relax,                           & ! u diag 5
  glob_debug                              & ! Debug flag
  )

!********************************************************************
! Meridional Wind
! DEPENDS ON: nudging_nudging_control
CALL nudging_nudging_control(             &
  v_name,                                 &
  grid_type_v,                            & ! Grid type
  field_type_v,                           & ! Field type
  ana_row_length,                         & ! global row length
  (ana_rows-1),                           & ! global rows
  file_levels,                            & ! Netcdf file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_vrow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_vrows_max,                         & ! Maximum row
  timestep,                               & ! model timestep
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  base_lambda,                            &
  delta_lambda,                           &
  (base_phi+0.5*delta_phi),               &
  delta_phi,                              &
  ana_base_lambda,                        &
  ana_delta_lambda,                       &
  ana_base_phi_v,                         &
  ana_delta_phi,                          &
  v_file_data(:,:,:,timestep1),           &
  v_file_data2(:,:,:,timestep2),          &
  surf_logpv_file_data(:,:,1,timestep1),  &
  surf_logpv_file_data2(:,:,1,timestep2), &
  frac,                                   & ! fraction between data
  model_levels,                           & ! Model levels
  p_vgrid_rho_levels                      & ! Model pressure
   ( v_row_length_min:v_row_length_max,   &
     v_rows_min:v_rows_max,               &
     1:model_levels ),                    &
  trop_pressure_vgrid(                    & ! Tropopause pressure
     v_row_length_min:v_row_length_max,   &
     v_rows_min:v_rows_max),              &
  v (                                     & ! nudged variable (v)
    v_row_length_min:v_row_length_max,    &
    v_rows_min:v_rows_max,                &
    1:model_levels ),                     &
  diag_v_data,                            & ! v diag 1
  diag_v_model,                           & ! v diag 2
  diag_v_ntend,                           & ! v diag 3
  diag_v_mtend,                           & ! v diag 4
  diag_v_relax,                           & ! v diag 5
  glob_debug                              & ! Debug flag
  )

! Drop a line letting us know where we are
IF(glob_debug > 10) THEN
  WRITE(OUT,'(A)')' NUDGING_MAIN: Performed nudging routines'
END IF

!****************************************************
! Wrap up at the end

!******************************************************************
! Synchronize theta for the polar rows
! DEPENDS ON: nudging_sync_poles
CALL nudging_sync_poles(                  &
  row_length,                             &  ! local row length
  rows,                                   &  ! local rows
  model_levels,                           &  ! model levels
  tdims_s%halo_i, tdims_s%halo_j,         &  ! x, y halo
  theta,                                  &  ! Variable synchonising
  sin_theta_latitude(1,1),                &  ! min latitude
  sin_theta_latitude(row_length,rows),    &  ! max latitude
  no_debug)                                  ! Debug flag

! Drop a line letting us know where we are
IF(glob_debug > 15) THEN
  WRITE(OUT,'(A)') 'NUDGING_MAIN: Synchronised at poles'
END IF

! Upload variables to the STASHWORK array
CALL nudging_update_stash

! Deallocate arrays used
CALL nudging_dealloc_data

! No longer the first time so set to false
l_first = .FALSE.

! Standard subroutine exit comment
IF(glob_debug > 0) WRITE(OUT,'(A)') ' Leaving NUDGING_MAIN'

IF (lhook) CALL dr_hook('NUDGING_MAIN',zhook_out,zhook_handle)
!*******************************************************

RETURN

!******************************************************
! Contained subroutines

CONTAINS

!*********************************************************************!
!                                                                     !
!     This routine loads data from netcdf files                       !
!     Separated from the main routine for the sake of tidiness        !
!                                                                     !
!*********************************************************************!
SUBROUTINE nudging_netcdf_loader(debug)

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE Field_Types
  USE UM_ParVars
  IMPLICIT NONE

! Local debug flag
INTEGER, INTENT(IN)  :: debug
REAL(KIND=jprb)      :: zhook_handle

!***************************************
! End of Header

IF (lhook) CALL dr_hook('NUDGING_NETCDF_LOADER',zhook_in,zhook_handle)

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,'(A)') ' NUDGING_NETCDF_LOADER: Entered Routine'
END IF

! Load filenames for file immediately before and after model timestep
! Don't worry about ECMWF naming convention, only necessary to change
! if you get data from another source with another naming convention
! DEPENDS ON: nudging_getfilename
CALL nudging_getfilename(               &
  i_year,                               &  ! Model year
  i_month,                              &  ! Model month
  i_day,                                &  ! Model day
  i_hour,                               &  ! Model hour
  dataname1,                            &  ! Return datafile name 1
  dataname2,                            &  ! Return datafile name 2
  timestep1,                            &  ! Return timestep 1
  timestep2,                            &  ! Return timestep 2
  0,                                    &  ! Select ECMWF naming convention
  glob_debug)                              ! Debug Flag

! If in debug mode drop a line giving the filenames
IF(debug > 15) THEN
  WRITE(OUT,'(A,A)')'NUDGING_NETCDF_LOADER: filename 1: ',dataname1
  WRITE(OUT,'(A,A)')'NUDGING_NETCDF_LOADER: filename 2: ',dataname2
END IF

! Load dimensions of the netcdf file
! DEPENDS on: NUDGING_NETCDF_DIMREADER
CALL nudging_netcdf_dimreader(         &
  dataname1,                           &  ! Filename
  ana_row_length,                      &  ! Expected length of rows
  ana_rows,                            &  ! Expected no. of rows
  totdims,                             &  ! Number of dimensions
  dim_names,                           &  ! Array with name of dimensions
  file_dims,                           &  ! Array containing size of dimensions
  .TRUE.,                              &  ! Check that (x/y) dimensions  match
  no_debug)                               ! DEbug Flag

! Load dimensions of the 2nd netcdf file
! DEPENDS on: NUDGING_NETCDF_DIMREADER
CALL nudging_netcdf_dimreader(         &
  dataname2,                           &  ! Filename
  ana_row_length,                      &  ! Expected length of rows
  ana_rows,                            &  ! Expected length of rows
  totdims,                             &  ! Number of dimensions
  dim_names,                           &  ! Array with name of dimensions
  file_dims2,                          &  ! Array containing size of dimensions
  .TRUE.,                              &  ! Check that (x/y) dimensions  match
  no_debug)                               ! DEbug Flag

! This covers switch from ECMWF-60 to ECMWF-91 levels
! Not perfect solution, but not a big problem
! To solve it would have to interpolate in z then time
! which is inefficient, so live with this small error
IF(file_dims(dimindex_lev) /= file_dims2(dimindex_lev)) THEN

  dataname2 = dataname1

  WRITE(OUT,'(5A,2(I4,A))')'NUDGING_MAIN: WARNING- overwriting ',      &
   dataname2, 'with', dataname1, 'as different number of levels (',    &
  file_dims(dimindex_lev), ' AND ', file_dims2(dimindex_lev), ')'

END IF

! Convert the dimensions to indiviual ints
file_row_length   = file_dims( dimindex_lon)      ! (theta) rows
file_rows         = file_dims( dimindex_lat)      ! (theta) row length
file_levels       = file_dims( dimindex_lev)      ! model levels
file_timesteps    = file_dims( dimindex_tim)      ! timesteps
file_urow_length  = file_dims( dimindex_ulon)     ! (U) rows
file_vrows        = file_dims( dimindex_vlat)     ! (V) row length

! If in debug mode write the dimensions of the netcdf file
IF(debug > 15) THEN
  WRITE(OUT,'(A,7(1X,I4))') 'NUDGING_MAIN: NetCDF dimensions are ',    &
   file_row_length, file_rows, file_levels, file_timesteps,            &
   file_urow_length, file_vrows, a_steps_per_hr
END IF

! Clear data arrays
IF(ALLOCATED(temp_file_data))        DEALLOCATE(temp_file_data)
IF(ALLOCATED(temp_file_data2))       DEALLOCATE(temp_file_data2)
IF(ALLOCATED(u_file_data))           DEALLOCATE(u_file_data)
IF(ALLOCATED(u_file_data2))          DEALLOCATE(u_file_data2)
IF(ALLOCATED(v_file_data))           DEALLOCATE(v_file_data)
IF(ALLOCATED(v_file_data2))          DEALLOCATE(v_file_data2)
IF(ALLOCATED(surf_logp_file_data))   DEALLOCATE(surf_logp_file_data)
IF(ALLOCATED(surf_logp_file_data2))  DEALLOCATE(surf_logp_file_data2)
IF(ALLOCATED(surf_logpu_file_data))  DEALLOCATE(surf_logpu_file_data)
IF(ALLOCATED(surf_logpu_file_data2)) DEALLOCATE(surf_logpu_file_data2)
IF(ALLOCATED(surf_logpv_file_data))  DEALLOCATE(surf_logpv_file_data)
IF(ALLOCATED(surf_logpv_file_data2)) DEALLOCATE(surf_logpv_file_data2)

! Allocate the data files using the dimensions loaded
ALLOCATE(temp_file_data                                                &
  (ana_row_length, ana_rows,     file_levels, file_timesteps))
ALLOCATE(temp_file_data2                                               &
  (ana_row_length, ana_rows,     file_levels, file_timesteps))
ALLOCATE(u_file_data                                                   &
  (ana_row_length, ana_rows,     file_levels, file_timesteps))
ALLOCATE(u_file_data2                                                  &
  (ana_row_length, ana_rows,     file_levels, file_timesteps))
ALLOCATE(v_file_data                                                   &
  (ana_row_length, (ana_rows-1), file_levels, file_timesteps))
ALLOCATE(v_file_data2                                                  &
  (ana_row_length, (ana_rows-1), file_levels, file_timesteps))
ALLOCATE(surf_logp_file_data                                           &
  (ana_row_length, ana_rows,     1,           file_timesteps))
ALLOCATE(surf_logp_file_data2                                          &
  (ana_row_length, ana_rows,     1,           file_timesteps))
ALLOCATE(surf_logpu_file_data                                          &
  (ana_row_length, ana_rows,     1,           file_timesteps))
ALLOCATE(surf_logpu_file_data2                                         &
  (ana_row_length, ana_rows,     1,           file_timesteps))
ALLOCATE(surf_logpv_file_data                                          &
  (ana_row_length, (ana_rows-1), 1,           file_timesteps))
ALLOCATE(surf_logpv_file_data2                                         &
  (ana_row_length, (ana_rows-1), 1,           file_timesteps))

! If in debug mode drop a line letting us know where we are
IF(debug > 15) THEN
  WRITE(OUT,'(A)') ' NUDGING_NETCDF_LOADER: ALLOCATED DATA ARRAYS'
END IF

! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(         &
  dataname1,                             & ! Filename
  temp_name,                             & ! Variable name
  ana_row_length,                        & ! Row length
  ana_rows,                              & ! Number of rows
  file_levels,                           & ! Number of levels
  file_timesteps,                        & ! Number of timesteps
  temp_file_data,                        & ! Variable
  no_debug)                                ! Debug flag

! Obtain second value from the second data (netcdf) file
! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(        &
  dataname2,                            & ! Filename
  temp_name,                            & ! Variable name
  ana_row_length,                       & ! Row length
  ana_rows,                             & ! Number of rows
  file_levels,                          & ! Number of levels
  file_timesteps,                       & ! Number of timesteps
  temp_file_data2,                      & ! Variable
  no_debug)                               ! Debug flag

! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(        &
  dataname1,                            & ! Filename
  u_name,                               & ! Variable name
  ana_row_length,                       & ! Row length
  ana_rows,                             & ! Number of rows
  file_levels,                          & ! Number of levels
  file_timesteps,                       & ! Number of timesteps
  u_file_data,                          & ! Variable
  no_debug)                               ! Debug flag

! Obtain second value from the second data (netcdf) file
! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(       &
  dataname2,                           & ! Filename
  u_name,                              & ! Variable name
  ana_row_length,                      & ! Row length
  ana_rows,                            & ! Number of rows
  file_levels,                         & ! Number of levels
  file_timesteps,                      & ! Number of timesteps
  u_file_data2,                        & ! Variable
  no_debug)                              ! Debug flag

! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(       &
  dataname1,                           & ! Filename
  v_name,                              & ! Variable name
  ana_row_length,                      & ! Row length
  (ana_rows-1),                        & ! Number of rows
  file_levels,                         & ! Number of levels
  file_timesteps,                      & ! Number of timesteps
  v_file_data,                         & ! Variable
  no_debug)                              ! Debug flag

! Obtain second value from the second data (netcdf) file
! DEPENDS ON: nudging_netcdf_vareader_4d
CALL nudging_netcdf_vareader_4d(      &
  dataname2,                          & ! Filename
  v_name,                             & ! Variable name
  ana_row_length,                     & ! Row length
  (ana_rows-1),                       & ! Number of rows
  file_levels,                        & ! Number of levels
  file_timesteps,                     & ! Number of timesteps
  v_file_data2,                       & ! Variable
  no_debug)                             ! Debug flag

! If in debug mode drop a line letting us know where we are
IF(debug > 15) THEN
  WRITE(OUT,'(A)') ' NUDGING_NETCDF_LOADER: LOADED MAIN DATA ARRAYS'
END IF

! If we are using hybrid pressure levels load surface pressure
IF(ndg_analysis_source == 0) THEN

! Choose data convention (in one surf P is 4D another 3D)
! Note in 4D case height is 1, but still needs treating as
! a dimension
  IF(data_source == 0) THEN

! LOad data from first file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname1,                      &  ! Filename
      surfp_name,                     &  ! variable name
      ana_row_length,                 &  ! row length
      ana_rows,                       &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logp_file_data,            &  ! Variable
      no_debug)                          ! Debug flag

! Load data from second file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname2,                      &  ! filename
      surfp_name,                     &  ! variable name
      ana_row_length,                 &  ! row length
      ana_rows,                       &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logp_file_data2,           &  ! Variable
      no_debug)

! LOad data from first file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname1,                      &  ! Filename
      usurfp_name,                    &  ! variable name
      ana_row_length,                 &  ! row length
      ana_rows,                       &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logpu_file_data,           &  ! Variable
      no_debug)                          ! Debug flag

! Load data from second file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname2,                      &  ! filename
      usurfp_name,                    &  ! variable name
      ana_row_length,                 &  ! row length
      ana_rows,                       &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logpu_file_data2,          &  ! Variable
      no_debug)

! LOad data from first file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname1,                      &  ! Filename
      vsurfp_name,                    &  ! variable name
      ana_row_length,                 &  ! row length
      (ana_rows-1),                   &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logpv_file_data,           &  ! Variable
      no_debug)                          ! Debug flag

! Load data from second file
! DEPENDS ON: nudging_netcdf_vareader_4d
    CALL nudging_netcdf_vareader_4d(  &
      dataname2,                      &  ! filename
      vsurfp_name,                    &  ! variable name
      ana_row_length,                 &  ! row length
      (ana_rows-1),                   &  ! rows
      1,                              &  ! model levels
      file_timesteps,                 &  ! Timesteps
      surf_logpv_file_data2,          &  ! Variable
      no_debug)

  ELSE

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname1,                      &  ! Filename
      surfp_name,                     &  ! Variable names
      ana_row_length,                 &  ! Row length
      ana_rows,                       &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logp_file_data(:,:,1,:),   &  ! Variable
      no_debug)                          ! Debug flag

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname2,                      &  ! Filename
      surfp_name,                     &  ! Variable names
      ana_row_length,                 &  ! Row length
      ana_rows,                       &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logp_file_data2(:,:,1,:),  &  ! Variable
      no_debug)                          ! Debug flag

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname1,                      &  ! Filename
      usurfp_name,                    &  ! Variable names
      ana_row_length,                 &  ! Row length
      ana_rows,                       &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logpu_file_data(:,:,1,:),  &  ! Variable
      no_debug)                          ! Debug flag

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname2,                      &  ! Filename
      usurfp_name,                    &  ! Variable names
      ana_row_length,                 &  ! Row length
      ana_rows,                       &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logpu_file_data2(:,:,1,:), &  ! Variable
      no_debug)

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname1,                      &  ! Filename
      vsurfp_name,                    &  ! Variable names
      ana_row_length,                 &  ! Row length
      (ana_rows-1),                   &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logpv_file_data(:,:,1,:),  &  ! Variable
      no_debug)                          ! Debug flag

! Load data from first file
! DEPENDS ON: nudging_netcdf_vareader_3d
    CALL nudging_netcdf_vareader_3d(  &
      dataname2,                      &  ! Filename
      vsurfp_name,                    &  ! Variable names
      ana_row_length,                 &  ! Row length
      (ana_rows-1),                   &  ! rows
      file_timesteps,                 &  ! Timesteps
      surf_logpv_file_data2(:,:,1,:), &  ! Variable
      no_debug)

  END IF ! data source
END IF ! Analysis source

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,'(A)') 'NUDGING_NETCDF_LOADER: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_NETCDF_LOADER',zhook_out,zhook_handle)

END SUBROUTINE nudging_netcdf_loader

!************************************************************!
!                                                            !
!     This routine allocates required diagnotics and         !
!     temporary / local arrays                               !
!                                                            !
!********************************************* **************!
SUBROUTINE nudging_alloc_data

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE Field_Types
  USE UM_ParVars
  IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_ALLOC_DATA',zhook_in,zhook_handle)

IF(glob_debug > 15) THEN
  WRITE(OUT,'(A)') ' NUDGING_ALLOC_DATA: Allocating diag variables'
END IF

! Allocate the theta diag1 array
ALLOCATE(diag_theta_data(                                               &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_data(:,:,:) = 0.0

! Allocate the theta diag2 array
ALLOCATE(diag_theta_model(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_model(:,:,:) = 0.0

! Allocate the theta diag3 array
ALLOCATE(diag_theta_mtend(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_mtend(:,:,:) = 0.0

! Allocate the theta diag4 array
ALLOCATE(diag_theta_ntend(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_ntend(:,:,:) = 0.0

! Allocate the theta diag5 array
ALLOCATE(diag_theta_relax(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_relax(:,:,:) = 0.0

! Allocate the u diag1 array
ALLOCATE(diag_u_data(                                                   &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_data(:,:,:) = 0.0

! Allocate the u diag2 array
ALLOCATE(diag_u_model(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_model(:,:,:) = 0.0

! Allocate the u diag3 array
ALLOCATE(diag_u_mtend(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_mtend(:,:,:) = 0.0

! Allocate the u diag4 array
ALLOCATE(diag_u_ntend(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_ntend(:,:,:) = 0.0

! Allocate the u diag5 array
ALLOCATE(diag_u_relax(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_relax(:,:,:) = 0.0

! Allocate the v diag1 array
ALLOCATE(diag_v_data(                                                   &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_data(:,:,:) = 0.0

! Allocate the v diag2 array
ALLOCATE(diag_v_model(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_model(:,:,:) = 0.0

! Allocate the v diag3 array
ALLOCATE(diag_v_mtend(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_mtend(:,:,:) = 0.0

! Allocate the v diag4 array
ALLOCATE(diag_v_ntend(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_ntend(:,:,:) = 0.0

! Allocate the v diag5 array
ALLOCATE(diag_v_relax(                                                 &
      v_row_length_min:v_row_length_max,                               &
      v_rows_min:v_rows_max,                                           &
      1:model_levels))

diag_v_relax(:,:,:) = 0.0

! Allocate pressure on rho levels on the u and v grids
ALLOCATE(p_ugrid_rho_levels(                                           &
      u_row_length_min : u_row_length_max,                             &
      u_rows_min : u_rows_max,                                         &
      1 : model_levels))

p_ugrid_rho_levels(:,:,:) = 0.0

ALLOCATE(p_vgrid_rho_levels(                                           &
     v_row_length_min : v_row_length_max,                              &
     v_rows_min : v_rows_max,                                          &
     1 : model_levels))

p_vgrid_rho_levels(:,:,:) = 0.0

 ! Allocate the tropopause pressure array 
ALLOCATE(trop_pressure(                                                &  
     theta_row_length_min:theta_row_length_max,                        &
     theta_rows_min:theta_rows_max))  

trop_pressure(:,:) = 0.0  

! Allocate the tropopause pressure array  
ALLOCATE(trop_pressure_ugrid(                                          &  
     u_row_length_min:u_row_length_max,                                & 
     u_rows_min:u_rows_max))  

trop_pressure_ugrid(:,:) = 0.0  

! Allocate the tropopause pressure array  
ALLOCATE(trop_pressure_vgrid(                                          &
     v_row_length_min:v_row_length_max,                                & 
     v_rows_min:v_rows_max))  

trop_pressure_vgrid(:,:) = 0.0 

IF (lhook) CALL dr_hook('NUDGING_ALLOC_DATA',zhook_out,zhook_handle)

END SUBROUTINE nudging_alloc_data

!************************************************************!
!                                                            !
!     This routine downloads specific variables from the     !
!     D1 array into specific arrays                          !
!     Section & Item numbers are currently set in            !
!     NUDGING_D1_DEFS module                                 !
!                                                            !
!********************************************* **************!
SUBROUTINE nudging_getd1_data

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE Field_Types
  USE UM_ParVars
  IMPLICIT NONE

INTEGER           :: d1_addr_st,d1_var_len
                    ! start & length of variable in D1
REAL,ALLOCATABLE  :: buff_temp(:)
INTEGER           :: debug
REAL(KIND=jprb)   :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_GETD1_DATA',zhook_in,zhook_handle)

debug = glob_debug

! Get grid type information for prognostic variables
! Currently required in some nudging subroutines

! Loop over sections & Items
DO i = 1, no_obj_d1(m_atm_modl)

  IF ( d1_addr(d1_section,i,m_atm_modl) == sect_prog ) THEN
! Section for prognostic
    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
! Item no.s
      CASE(u_index)
        grid_type_u = d1_addr(d1_grid_type,i,m_atm_modl)
        field_type_u = fld_type_u  

      CASE(v_index)
        grid_type_v = d1_addr(d1_grid_type,i,m_atm_modl)
        field_type_v = fld_type_v

      CASE(theta_index)
        grid_type_theta = d1_addr(d1_grid_type,i,m_atm_modl)
        field_type_theta = fld_type_p

      CASE DEFAULT

! Do nothing
    END SELECT

  ELSE IF( d1_addr(d1_section,i,m_atm_modl) == sect_tropht ) THEN
!Section for climate diagnostics (for tropopause pressure)
    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
! Item no. for tropopause pressure (set in nudging_d1_defs module)

     CASE( trop_pres_index)
       d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
       d1_var_len = d1_addr(d1_length,i,m_atm_modl)

       IF( d1_var_len /= row_length*rows ) THEN
          WRITE(nmessage,'(3(A,1X,I9))')                        &
           'Mismatch in variable size :STASH vs trop pressure', &
            trop_pres_index,'. Expected',                       &
            row_length*rows,': Actual', d1_var_len

            CALL ereport('NUDGE_GET_D1',trop_pres_index,nmessage)
       END IF

       ALLOCATE(buff_temp(row_length*rows))

       buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

       trop_pressure(:,:) =                                      &
       RESHAPE(buff_temp(:),(/row_length,rows/))

       DEALLOCATE(buff_temp)

     END SELECT

! Extract values for post-nudged variables of previous timestep
! This should ensure bit comparability for CRUNs. Otherwise these would be
! initialised to zero and diag_model_tendency ( pre-nudging - previous_post-nudged)
! will get a different value
  ELSE IF ( d1_addr(d1_section,i,m_atm_modl) == sect_nudge ) THEN
! Section for nudging

    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
! Item no.s (set in nudging_d1_defs module)

! Theta diag 2
      CASE( tdiag_model_index )

        d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
        d1_var_len = d1_addr(d1_length,i,m_atm_modl)

        IF( d1_var_len /= row_length*rows*model_levels ) THEN
          WRITE(nmessage,'(3(A,1X,I9))')                        &
          'Mismatch in variable size :STASH vs nudging diag',   &
           tdiag_model_index,'. Expected',                      &
           row_length*rows*model_levels,': Actual', d1_var_len

         CALL ereport('NUDGE_GET_D1',tdiag_model_index,nmessage)
        END IF

        ALLOCATE(buff_temp(row_length*rows*model_levels))

        buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

        diag_theta_model(:,:,:) =                                &
         RESHAPE(buff_temp(:),(/row_length,rows,model_levels/))

        DEALLOCATE(buff_temp)

! u diag 2
      CASE( udiag_model_index )

        d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
        d1_var_len = d1_addr(d1_length,i,m_atm_modl)

        IF( d1_var_len /= row_length*rows*model_levels ) THEN
          WRITE(nmessage,'(3(A,1X,I9))')                         &
          'Mismatch in variable size :STASH vs nudging diag',    &
           udiag_model_index,'. Expected',                       &
           row_length*rows*model_levels,': Actual', d1_var_len

         CALL ereport('NUDGE_GET_D1',udiag_model_index,nmessage)
        END IF

        ALLOCATE(buff_temp(row_length*rows*model_levels))

        buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

        diag_u_model(:,:,:) =                                &
         RESHAPE(buff_temp(:),(/row_length,rows,model_levels/))

        DEALLOCATE(buff_temp)

! v diag 2
      CASE( vdiag_model_index )

        d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
        d1_var_len = d1_addr(d1_length,i,m_atm_modl)

        IF( d1_var_len /= row_length*n_rows*model_levels ) THEN
          WRITE(nmessage,'(3(A,1X,I9))')                         &
          'Mismatch in variable size :STASH vs nudging diag',    &
           vdiag_model_index,'. Expected',                       &
           row_length*rows*model_levels,': Actual', d1_var_len

         CALL ereport('NUDGE_GET_D1',vdiag_model_index,nmessage)
        END IF

        ALLOCATE(buff_temp(row_length*n_rows*model_levels))

        buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

        diag_v_model(:,:,:) =                                &
         RESHAPE(buff_temp(:),(/row_length,n_rows,model_levels/))

        DEALLOCATE(buff_temp)

    END SELECT ! Items in Nudging Section

  END IF       ! If Prognostics/Clim Diags/ Nudging Section

END DO         ! Loop over D1 items

IF (lhook) CALL dr_hook('NUDGING_GETD1_DATA',zhook_out,zhook_handle)

END SUBROUTINE nudging_getd1_data

!************************************************************!
!                                                            !
!     This routine copies the nudging diagnostics into       !
!     STASH using the copy_diag routine                      !
!                                                            !
!********************************************* **************!
SUBROUTINE nudging_update_stash

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE Field_Types
  USE UM_ParVars
  IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_UPDATE_STASH',zhook_in,zhook_handle)

! Check whether Nudging diagnostics and Tropopause pressure
! have been requested in the STASH
IF ( .NOT. sf(trop_pres_index,sect_tropht)  .OR.                 &
     .NOT. sf(tdiag_data_index,sect_nudge)  .OR.                 &
     .NOT. sf(tdiag_model_index,sect_nudge)  .OR.                 &
     .NOT. sf(tdiag_ntend_index,sect_nudge)  .OR.                 &
     .NOT. sf(tdiag_mtend_index,sect_nudge)  .OR.                 &
     .NOT. sf(tdiag_relax_index,sect_nudge)  .OR.                 &
     .NOT. sf(udiag_data_index,sect_nudge)   .OR.                 &
     .NOT. sf(udiag_model_index,sect_nudge)  .OR.                 &
     .NOT. sf(udiag_ntend_index,sect_nudge)  .OR.                 &
     .NOT. sf(udiag_mtend_index,sect_nudge)  .OR.                 &
     .NOT. sf(vdiag_data_index,sect_nudge)   .OR.                 &
     .NOT. sf(vdiag_model_index,sect_nudge)  .OR.                 &
     .NOT. sf(vdiag_ntend_index,sect_nudge)  .OR.                 &
     .NOT. sf(vdiag_mtend_index,sect_nudge)  .OR.                 &
     .NOT. sf(vdiag_relax_index,sect_nudge)  ) THEN

  icode = sect_nudge

  CALL ereport('NUDGING_MAIN',icode,                              &
    'Required Nudging diags or Tropopause pressure not in STASH requests' )
END IF

icode = 0

! Theta Diag 1
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(tdiag_data_index,sect_nudge,im_index)),        &
      diag_theta_data(:,:,:),                                   &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_data_index,sect_nudge,im_index)),&
      len_stlist, stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_data_index,icode,nmessage)

IF ( icode > 0 )                                                &

  CALL ereport('NUDGE COPY DIAG',tdiag_data_index,nmessage)

! Theta Diag 2
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(tdiag_model_index,sect_nudge,im_index)),       &
      diag_theta_model(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_model_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,tdiag_model_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',tdiag_model_index,nmessage)

! Theta Diag 3
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(tdiag_ntend_index,sect_nudge,im_index)),       &
      diag_theta_ntend(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_ntend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',tdiag_ntend_index,nmessage)

! Theta Diag 4
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(tdiag_mtend_index,sect_nudge,im_index)),       &
      diag_theta_mtend(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,tdiag_mtend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',tdiag_mtend_index,nmessage)

! Theta Diag 5
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(tdiag_relax_index,sect_nudge,im_index)),       &
      diag_theta_relax(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_relax_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_relax_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',tdiag_relax_index,nmessage)

! U Diag 1
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(udiag_data_index,sect_nudge,im_index)),        &
      diag_u_data(:,:,:),                                       &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_data_index,sect_nudge,im_index)),&
      len_stlist,stash_levels,num_stash_levels+1,               &
      atmos_im,sect_nudge,udiag_data_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',udiag_data_index,nmessage)

! U Diag 2
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(udiag_model_index,sect_nudge,im_index)),       &
      diag_u_model(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_model_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_model_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',udiag_model_index,nmessage)

! U Diag 3
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(udiag_ntend_index,sect_nudge,im_index)),       &
      diag_u_ntend(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_ntend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',udiag_ntend_index,nmessage)

! U Diag 4
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(udiag_mtend_index,sect_nudge,im_index)),       &
      diag_u_mtend(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_mtend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',udiag_mtend_index,nmessage)

! U Diag 5
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(udiag_relax_index,sect_nudge,im_index)),       &
      diag_u_relax(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_relax_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_relax_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',udiag_relax_index,nmessage)

! V Diag 1
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(vdiag_data_index,sect_nudge,im_index)),        &
      diag_v_data(:,:,:),                                       &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_data_index,sect_nudge,im_index)),&
      len_stlist, stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,vdiag_data_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',vdiag_data_index,nmessage)

! V Diag 2
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(vdiag_model_index,sect_nudge,im_index)),       &
      diag_v_model(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_model_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,vdiag_model_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',vdiag_model_index,nmessage)

! V Diag 3
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(vdiag_ntend_index,sect_nudge,im_index)),       &
      diag_v_ntend(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_ntend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',vdiag_ntend_index,nmessage)

! V Diag 4
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(vdiag_mtend_index,sect_nudge,im_index)),       &
      diag_v_mtend(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_mtend_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',vdiag_mtend_index,nmessage)

! V Diag 5
! DEPENDS ON: copydiag_3d
CALL copydiag_3d (                                              &
      st_work(si(vdiag_relax_index,sect_nudge,im_index)),       &
      diag_v_relax(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_relax_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_relax_index,icode,nmessage)

IF ( icode > 0 )                                                &

   CALL ereport('NUDGE COPY DIAG',vdiag_relax_index,nmessage)

IF(glob_debug > 15) THEN
  WRITE(OUT,'(A)') ' NUDGING_MAIN: Finished Uploading data'
END IF

IF (lhook) CALL dr_hook('NUDGING_UPDATE_STASH',zhook_out,zhook_handle)

END SUBROUTINE nudging_update_stash

!************************************************************!
!                                                            !
!     This routine deallocates allocated diagnotics and      !
!     temporary / local arrays                               !
!                                                            !
!********************************************* **************!
SUBROUTINE nudging_dealloc_data

  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE Field_Types
  USE UM_ParVars
  IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_DEALLOC_DATA',zhook_in,zhook_handle)

IF(glob_debug > 10) THEN
  WRITE(OUT,'(A)') 'NUDGING_DEALLOC_DATA: Deallocating diag variables'
END IF

! Deallocate derived arrays
IF(ALLOCATED(p_ugrid_rho_levels))  DEALLOCATE(p_ugrid_rho_levels)
IF(ALLOCATED(p_vgrid_rho_levels))  DEALLOCATE(p_vgrid_rho_levels)
IF(ALLOCATED(trop_pressure))       DEALLOCATE(trop_pressure) 
IF(ALLOCATED(trop_pressure_ugrid)) DEALLOCATE(trop_pressure_ugrid) 
IF(ALLOCATED(trop_pressure_vgrid)) DEALLOCATE(trop_pressure_vgrid) 

! Deallocate diagnostic arrays
IF(ALLOCATED(diag_theta_data))     DEALLOCATE(diag_theta_data)
IF(ALLOCATED(diag_theta_model))    DEALLOCATE(diag_theta_model)
IF(ALLOCATED(diag_theta_ntend))    DEALLOCATE(diag_theta_ntend)
IF(ALLOCATED(diag_theta_mtend))    DEALLOCATE(diag_theta_mtend)
IF(ALLOCATED(diag_theta_relax))    DEALLOCATE(diag_theta_relax)

IF(ALLOCATED(diag_u_data))         DEALLOCATE(diag_u_data)
IF(ALLOCATED(diag_u_model))        DEALLOCATE(diag_u_model)
IF(ALLOCATED(diag_u_ntend))        DEALLOCATE(diag_u_ntend)
IF(ALLOCATED(diag_u_mtend))        DEALLOCATE(diag_u_mtend)
IF(ALLOCATED(diag_u_relax))        DEALLOCATE(diag_u_relax)

IF(ALLOCATED(diag_v_data))         DEALLOCATE(diag_v_data)
IF(ALLOCATED(diag_v_model))        DEALLOCATE(diag_v_model)
IF(ALLOCATED(diag_v_ntend))        DEALLOCATE(diag_v_ntend)
IF(ALLOCATED(diag_v_mtend))        DEALLOCATE(diag_v_mtend)
IF(ALLOCATED(diag_v_relax))        DEALLOCATE(diag_v_relax)

IF (lhook) CALL dr_hook('NUDGING_DEALLOC_DATA',zhook_out,zhook_handle)

END SUBROUTINE nudging_dealloc_data

!##########################################################################

END SUBROUTINE nudging_main1

