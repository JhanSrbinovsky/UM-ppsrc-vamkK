! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: UM_SHELL -------------------------------------------------
!
!  Purpose: Control routine for the Atm Model.
!           Acquires size information needed for dynamic allocation of
!           configuration-dependent arrays and calls U_MODEL (the
!           master control routine) to allocate the arrays and perform
!           the top-level control functions and timestepping.
!
!  External documentation: On-line UM document C1 - The top-level
!                          dynamic allocation
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

SUBROUTINE um_shell

  USE dynamics_input_mod, ONLY : &
      l_endgame
  USE UM_Config, ONLY : &
      appInit, &
      exe_UM
  USE mpl, ONLY :                                                              &
      mpl_max_processor_name,                                                  &
      mpl_thread_multiple,                                                     &
      mpl_thread_serialized,                                                   &
      mpl_thread_funneled,                                                     &
      mpl_thread_single





!$ USE omp_lib           ! Note OpenMP sentinel

  USE filenamelength_mod, ONLY :                                               &
      filenamelength, datawnamelength, runidnamelength

  USE atm_fields_bounds_Mod
  USE Halo_Exchange, ONLY :                                                    &
      Halo_Exchange_Init
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  USE io,                 ONLY :                                               &
      ioInit,                                                                  &
      ioShutdown,                                                              &
      file_close,                                                              &
      is_Unit_Open
  USE IOS_Client_Queue
  USE ios
  USE IOS_Constants,      ONLY :                                               &
      IOS_maxServers
  USE IOS_Init,           ONLY :                                               &
      IOS_Setup
  USE IOS_Stash_Common,   ONLY :                                               &
      isUsingAsyncStash,                                                       &
      isUSingAsyncDumps
  USE IOS_Stash, ONLY :                                                        &
      IOS_stash_client_fini
  USE IOS_Model_Geometry, ONLY :                                               &
      IOS_Client_Geometry_Init
  USE MPPIO_job_control,  ONLY :                                               &
      jobCntrlInit,                                                            &
      jobCntrl
  USE MPPIO_job_control_common, ONLY :                                         &
      jc_wakeup
  USE mppio_file_utils, ONLY :                                                 &
      mppio_file_utils_init
  USE model_file,     ONLY :                                                   &
      model_file_managed,                                                      &
      model_file_close
  USE mpp_conf_mod,       ONLY :                                               &
      extended_halo_size_ew,                                                   &
      extended_halo_size_ns
  USE ereport_mod, ONLY : ereport,ereport_finalise
  USE UM_ParVars
  USE Decomp_DB
  USE PrintStatus_mod

  USE gcom_mod, ONLY : gc_alltoall_multi, gc_alltoall_version
  USE decomp_params, ONLY :                                                    &
      decomp_standard_atmos,                                                   &
      decomp_unset
  USE rimtypes
  USE lbc_mod
  USE um_input_control_mod,  ONLY:  model_domain, l_oasis, h_sect,             &
       maxsects
  USE ppxlook_mod, ONLY : ppxrecs
! version_mod items required by cstash.h
  USE version_mod, ONLY :                                                      &
      nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,                    &
      nlevp, npslevp, npslistp, outfile_s, outfile_e
  USE Submodel_Mod, ONLY: internal_model_index, a_im
  USE chsunits_mod, ONLY : nunits
  IMPLICIT NONE
!
!  Local parameters
!
  CHARACTER(LEN=*) RoutineName
  PARAMETER (RoutineName = 'UM_SHELL')
!
!  Local variables
!
  INTEGER icode       ! Work - Internal return code
  INTEGER(KIND=integer32) :: icode_OASIS ! 32-bit OASIS return code 
  INTEGER istatus     ! RETURN STATUS FROM OPEN
  INTEGER loop_pe
  INTEGER loop_pe_start
  INTEGER loop_pe_stop

  CHARACTER(LEN=filenamelength) ::  filename  = "dummy filename"
               ! RETURN FILENAME FROM GET_FILE
  CHARACTER(LEN=512) :: cmessage ! Work - Internal error message
  INTEGER :: atm_nprocx          ! number of procs EW for atmosphere
  INTEGER :: atm_nprocy          ! number of procs NS for atmosphere
  INTEGER :: err                 ! error return from FORT_GET_ENV
  INTEGER :: nproc_um_npes       ! Total number of atmos PEs (i.e. EW*NS)

  CHARACTER(LEN=10) c_thread          ! to get nproc_x and nproc_y from
  CHARACTER(LEN=8) c_nproc            ! to get nproc_x and nproc_y from
  
  ! environment variables.
  
  ! to hold the name of the parallel
  ! executable script
  CHARACTER(LEN=filenamelength)  :: dummy_env = "dummy path"      
 
  ! file to write stdout to
  CHARACTER(LEN=filenamelength)  :: stdout_filename = "dummy stdout"
  
  ! base of filename
  CHARACTER(LEN=filenamelength)  :: stdout_basename = "dummy stdout"
  
  ! File to write STASH requests to
  CHARACTER(LEN=filenamelength)  :: stash_filename  = "dummy stash"
  
  ! value of $DATAW
  CHARACTER(LEN=datawnamelength) :: dataw_char = " "

  ! value of $RUNID
  CHARACTER(LEN=runidnamelength) :: runid_char = " "
!
  INTEGER um_nam_max_seconds
!
  CHARACTER(LEN=8) c_nam_max_seconds
!
!  Configuration-dependent sizes for dynamic arrays
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
!
!  Super array sizes for dynamic allocation in U_MODEL/U_MODEL_4A
!
! TYPSZSP super arrays lengths not dependent on sub-models
      INTEGER :: SPD1_LEN
      INTEGER :: SPSTS_LEN
      INTEGER :: SPBND_LEN
! TYPSZSP end
! TYPSZSPA super arrays lengths (atmosphere)
      INTEGER :: A_SPDUM_LEN
      INTEGER :: A_SPPTR_LEN
      INTEGER :: A_SPCON_LEN
      INTEGER :: A_SPINF_LEN
      INTEGER :: A_SPANC_LEN
      INTEGER :: A_SPBND_LEN
      INTEGER :: A_SPSTS_LEN
! TYPSZSPA end
! TYPSZSPC uper arrays lengths (atmosphere-ocean coupled)
      INTEGER :: AO_SPCPL_LEN
! TYPSZSPC end
!
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
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

!L
!
!  Localized sizes for ocean decomposition:
  INTEGER                                                                      &
      row_length_oce , rows_oce

  INTEGER sect_err, rnl_err, um_rnl_skip

  CHARACTER(LEN=8) c_um_rnl_skip
  CHARACTER(LEN=8) ch_date2   !  Date returned from date_and_time
  CHARACTER(LEN=10) ch_time2  !  Time returned from date_and_time

!   Fortran unit numbers
  INTEGER                      :: nftppxref
  INTEGER                      :: nftstmstu
  data nftppxref/22/,nftstmstu/2/
! Variables needed to close all the files
  INTEGER                      :: i

! Variables for IO Server setup
  LOGICAL                      :: isIOServer
  INTEGER                      :: numIOServers
  INTEGER                      :: atm_comm
  INTEGER                      :: errorCode
  CHARACTER(LEN=32)                 :: c_io_pes
  CHARACTER(LEN=10)                 :: thread_level_setc
  INTEGER                      :: thread_level_set

  INTEGER :: dummy_comm     ! Dummy communicator for OASIS

  CHARACTER(LEN=mpl_max_processor_name) :: env_myhost
  INTEGER                               :: env_myhost_len

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  REAL(KIND=jprb)               :: zhook_rendez_vous
  LOGICAL(KIND=jpim)            :: luser_comm

!      IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
  cmessage = ' '
  numIOServers = 0  ! Initialise IO server variable

! For IBM platform, turn on signal handling
!L----------------------------------------------------------------------
!----------------------------------------------------------------------
! 1.0 Initialise Message Passing Libraries
!

! Get the atmosphere decomposition
  CALL fort_get_env('UM_THREAD_LEVEL',15,c_thread,10,err)
  IF (err  /=  0) THEN
    WRITE(6,'(A,A)') 'Warning: Environment variable UM_THREAD_LEVEL has ',     &
        'not been set.'
    WRITE(6,'(A)') 'Setting thread_level to multiple'
    thread_level_setc = 'MULTIPLE'
  ELSE
    READ(c_thread,'(A10)') thread_level_setc
  END IF

  SELECT CASE (thread_level_setc)
  CASE ('MULTIPLE')
    thread_level_set = mpl_thread_multiple
  CASE ('SERIALIZED')
    thread_level_set = mpl_thread_serialized
  CASE ('FUNNELED')
    thread_level_set = mpl_thread_funneled
  CASE ('SINGLE')
    thread_level_set = mpl_thread_single
  CASE default
    WRITE(6,'(A,A,A)') 'Warning: Thread level ', thread_level_setc,            &
        ' not recognised, setting to MULTIPLE.'
    thread_level_set = mpl_thread_multiple
  END SELECT

  CALL fort_get_env('UM_ATM_NPROCX',13,c_nproc,8,err)
  IF (err  /=  0) THEN
    cmessage = 'Warning: Environment variable UM_ATM_NPROCX has '              &
        //'not been set. Setting nproc_x to 1.'
        ! Can't ereport before gc_init, so writing to stdout manually
    WRITE(6,'(A)') cmessage
    atm_nprocx=1
  ELSE
    READ(c_nproc,'(I4)') atm_nprocx
  END IF
  CALL fort_get_env('UM_ATM_NPROCY',13,c_nproc,8,err)
  IF (err  /=  0) THEN
    cmessage = 'Warning: Environment variable UM_ATM_NPROCY has '              &
        //'not been set. Setting nproc_y to 1.'
        ! Can't ereport before gc_init, so writing to stdout manually
    WRITE(6,'(A)') cmessage
    atm_nprocy=1
  ELSE
    READ(c_nproc,'(I4)') atm_nprocy
  END IF

! Calculate total number of atmos processors required:
  nproc_um_npes = atm_nprocx * atm_nprocy

! The total number of processors required (nproc_max) is determined
!  by a call to gc_init/gc_init_thread:

! Standard UM GCOM initialisation when OASIS is not used

  CALL gc_init_thread(mype,nproc_max, thread_level_set)

      ! Use DrHook over global communicator
  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Get number of I/O PEs from the environment
  CALL fort_get_env('FLUME_IOS_NPROC',15,c_io_pes,5,err)
  IF (err /= 0) THEN
        ! If not specified, try to work out a valid number
    numIOServers=nproc_max-nproc_um_npes
    icode=-10
    WRITE(cmessage,'(A,A,I4)')                                                 &
        'FLUME_IOS_NPROC environment variable not set',                        &
        ', I/O PE count set to ',numIOServers
    CALL ereport(routinename,icode,cmessage)
  ELSE
    READ (c_io_pes,'(I5)') numIOServers
  END IF

  IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
    icode=-10
    WRITE(cmessage,'(A,I4)')                                                   &
        'I/O PE count is outside allowed range: ',numIOServers
    CALL ereport(routinename,icode,cmessage)

        ! try to work out a valid number
    numIOServers=nproc_max-nproc_um_npes
    IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
      icode=10
      WRITE(cmessage,'(A)')                                                    &
          'A valid I/O PE count could not be set'
      CALL ereport(routinename,icode,cmessage)
    END IF
  ELSE
    IF (mype == 0) THEN
      WRITE(6,'(A,I2,A)') 'Enabling ', numIOServers, ' I/O PEs'
    END IF
  END IF

! Check the breakdown of processors requested by environment variables
! (UM_ATM_NPROCX * UM_ATM_NPROCY + FLUME_IOS_NPROC) matches the total
! number obtained by GCOM
  IF (nproc_max /= nproc_um_npes + numIOServers) THEN
    icode = 100
    WRITE(cmessage,'(A,i7,A,i7,A)')                                            &
        'UM started on ', nproc_max,                                           &
        ' PEs but ',                                                           &
        nproc_um_npes+numIOServers,                                            &
        ' asked for. Please adjust decomposition'
    CALL ereport(routinename,icode,cmessage)
  END IF

  IF (nproc_max  <   0) THEN
    WRITE(6,'(A)') 'Parallel initialisation failed'
    GO TO 9999
  ELSE

! Get values of DATAW and RUNID from environment
    CALL fort_get_env('DATAW',5,dataw_char,datawnamelength,err)
    IF (err  /=  0) THEN
      cmessage = 'UMSHELL : Failed to get value of $DATAW'
      GO TO 9999
    END IF
    CALL fort_get_env('RUNID',5,runid_char,runidnamelength,err)
    IF (err  /=  0) THEN
      cmessage = 'UMSHELL : Failed to get value of $RUNID'
      GO TO 9999
    END IF

! Send output to unique filename on every PE

    CALL fort_get_env('UM_STDOUT_FILE',14,stdout_basename,                     &
        filenamelength,err)
    IF (err  /=  0) THEN
! Environment variable UM_STDOUT_FILE has not been set, so we will
! construct a default stdout_basename of $DATAW/pe_output/$RUNID.fort6.pe
      stdout_basename=TRIM(dataw_char)//'/pe_output/'//                        &
          TRIM(runid_char)//'.fort6.pe'
    END IF

    stash_filename=TRIM(dataw_char)//'/'// TRIM(runid_char)//                  &
        '.stash'

! Now add PE number (mype) to stdout_basename to get the complete
! stdout_filename for this PE.

    IF (mype  <   10) THEN
      WRITE(stdout_filename,'(A,I1)') TRIM(stdout_basename),mype
    ELSE IF (mype  <   100) THEN
      WRITE(stdout_filename,'(A,I2)') TRIM(stdout_basename),mype
    ELSE IF (mype  <   1000) THEN
      WRITE(stdout_filename,'(A,I3)') TRIM(stdout_basename),mype
    ELSE IF (mype  <   10000) THEN
      WRITE(stdout_filename,'(A,I4)') TRIM(stdout_basename),mype
    ELSE IF (mype  <   100000) THEN
      WRITE(stdout_filename,'(A,I5)') TRIM(stdout_basename),mype
    ELSE
      WRITE(stdout_filename,'(A,I6)') TRIM(stdout_basename),mype
    END IF

! and close unit 6, then reopen to new filename

    CLOSE(6)
    OPEN(6,FILE=stdout_filename)
! Force a close with a delete action - so if there is an existing
! unit6 output file it will be deleted, and the output from this
! run will go to a fresh file
    CLOSE(6,status='DELETE')
    OPEN(6,FILE=stdout_filename)

    IF (mype == 0) THEN
      CLOSE(200)
      OPEN(200,FILE=stash_filename,status='REPLACE')
    END IF

    CALL appInit(exe_UM)

! Set GCOM to use the alternative version of RALLTOALLE
! throughout the run
    CALL gc_setopt(gc_alltoall_version, gc_alltoall_multi, err)

    WRITE(6,'(/,A,/)') '**************************** PROCESSOR '//             &
        'INFORMATION ****************************'
    WRITE(6,'(I7,A)') nproc_max,' Processors initialised.'
    CALL MPL_Get_processor_name(env_myhost, env_myhost_len, err)
    IF (err /= 0) THEN
      WRITE(6,'(A,I7)') 'I am PE ',mype
    ELSE
      WRITE(6,'(A,I7,A,A)') 'I am PE ',mype,' on ', TRIM(env_myhost)
    END IF
! Only want OpenMP section executing if OpenMP is compiled in,
! so protect by sentinal
!$OMP PARALLEL DEFAULT(NONE)
!$OMP MASTER
!$  WRITE(6,'(A,I2,A)') 'I am running with ',                                  &
!$      omp_get_num_threads(),' thread(s).'
!$  WRITE(6,'(A,I6)') 'OpenMP Specification: ',openmp_version
!$OMP END MASTER
!$OMP END PARALLEL
  END IF

!
! DEPENDS ON: timer
  CALL timer(RoutineName,1)

!----------------------------------------------------------------------
!
!    Open unit 5. All runtime control variables subsequently read in
!    from UNIT 5 by namelist.
!
  CALL get_file(5,filename,filenamelength,icode)
  OPEN(5,file=filename,iostat=istatus)
  IF (istatus /= 0) THEN
    icode=500
    WRITE(6,'(A)') ' ERROR OPENING FILE ON UNIT 5'
    WRITE(6,'(A,A)') ' FILENAME =',filename
    WRITE(6,'(A,I6)') ' IOSTAT =',istatus
    GO TO 9999
  END IF
!L------------------------------------------------------------------
!L 0.1 Get submodel/internal model components of model run.
!L
  icode=0
! DEPENDS ON: um_submodel_init
  CALL UM_Submodel_Init(ICODE)
  IF (icode  /=  0) THEN
    cmessage = 'Error calling UM_Submodel_init'
    GO TO 9999
  END IF


  IF (mype == 0) THEN
    CALL date_and_time(ch_date2, ch_time2)
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '********************************************'//&
                   '***********************************'
    WRITE(6,'(21A)')                                               &
      '**************** Start of UM RUN Job : ',                   &
      ch_time2(1:2),  ':',                                         &
      ch_time2(3:4),  ':',                                         &
      ch_time2(5:6),  ' on ',                                      &
      ch_date2(7:8),  '/',                                         & 
      ch_date2(5:6),  '/',                                         &
      ch_date2(1:4),  ' *****************'
    WRITE(6,'(A,/)') '*****************************************'// &
                   '**************************************'
  END IF


! Initialise I/O subsystem 
  CALL ioInit()

! Start I/O Server if required
      ! On exit the IO slave PEs have finished the whole job so can
      ! go to the end of the routine. FIXME check cpu counts against deco
  isIOServer=IOS_Setup(numIOServers)

      ! After IOS_SETUP the communicator may have changed
  CALL GC_Get_Communicator(atm_comm, errorCode)
  CALL MPL_Comm_Rank(atm_comm, mype, errorCode)

  IF (.NOT.isIOServer) THEN
    WRITE(6,'(A,I7)')'Running Atmospheric code as pe ',mype

!----------------------------------------------------------------------
!L
!L    Open unit 8 for server requests and send wakeup message
!L
    CALL mppio_file_utils_init()
    CALL jobCntrlInit()
! WAKEUP doesn't have a file argument, so provide a dummy string
    CALL jobCntrl(jc_wakeup,' little Susie, wake up!')

!L----------------------------------------------------------------------
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
    CALL InitPrintStatus()

! ----------------------------------------------------------------------
!  Read Control file on standard input.
!
! DEPENDS ON: readcntl
    CALL readcntl ( 5,icode,cmessage )
    IF (icode  >   0) GO TO 9999

!L----------------------------------------------------------------------
!L Allow Override of namelist input in operational environment.
!L
! DEPENDS ON: oper_emergency
    CALL Oper_Emergency

!----------------------------------------------------------------------
!L Call READLSTA to read namelists to control atmosphere integration
!L and diagnostic point print.
    IF ( l_endgame ) THEN
! DEPENDS ON: readlsta_4a
      CALL readlsta_4A()
    ELSE
! DEPENDS ON: readlsta
      CALL readlsta()
    END IF

!----------------------------------------------------------------------
!L 1.1 Get configuration-dependent sizes needed for dynamic allocation.
!L
! DEPENDS ON: readsize
    CALL readsize()

!L   Read history and control files for NRUN; also interim control
!L   file for CRUN, and housekeeping file for operational run.
!L
! DEPENDS ON: um_setup
    CALL um_setup(icode,cmessage)

    IF (icode >  0) GO TO 9999

! Decompose atmosphere data and find new local data size

    CALL decompose(decomp_standard_atmos,                                      &
        global_row_length,global_rows,model_levels,                            &
        river_rows, river_row_length,                                          &
        model_domain,                                                          &
        atm_nprocx, atm_nprocy,                                                &
        extended_halo_size_ew,                                                 &
        extended_halo_size_ns,                                                 &
        rimwidtha, nrima_max,row_length,rows )

! Set up the atmosphere decomposition in PARVARS
    CALL change_decomposition(decomp_standard_atmos,icode)

! Now we have a decomposition initialise the halo swap module
    CALL Halo_Exchange_Init()

    IF (icode  /=  0) GO TO 9999

! Output range of gridpoints for each PE
    IF (PrintStatus >= PrStatus_Diag) THEN
      loop_pe_start=mype
      loop_pe_stop =mype
      IF (mype  ==  0) THEN
        loop_pe_start=0
        loop_pe_stop=atm_nprocx*atm_nprocy-1
      END IF

      WRITE(6,'(A)')''
      WRITE(6,'(A)')'Range of gridpoints for processing elements:'
      WRITE(6,'(A)')''
      WRITE(6,'(A7,A2,A15,A2,A15)')                                            &
          '     PE',' |','  West -   East',' |',' South -  North'
      WRITE(6,'(A7,A2,A15,A2,A15)')                                            &
          '-------','-+','---------------','-+','---------------'
      DO loop_pe = loop_pe_start,loop_pe_stop
        WRITE(6,'(I7,A2,I6,A3,I6,A2,I6,A3,I6)')                                &
            loop_pe,' |',                                                      &
            g_datastart_f(1,1,loop_pe),' - ',                                  &
            g_datastart_f(1,1,loop_pe)+g_blsize(1,1,loop_pe)-1,' |',           &
            g_datastart_f(2,1,loop_pe),' - ',                                  &
            g_datastart_f(2,1,loop_pe)+g_blsize(2,1,loop_pe)-1
      END DO
      WRITE(6,'(A7,A2,A15,A2,A15)')                                            &
          '-------','-+','---------------','-+','---------------'

    END IF

! Call DERVSIZE (the call in READSIZE has been deleted)

    icode=0
! DEPENDS ON: dervsize
    CALL dervsize(icode,cmessage)
    IF (icode  /=  0) GO TO 9999

!     Ensure that domain decomposition is set for Atmosphere
    CALL change_decomposition (decomp_standard_atmos,icode)
    IF (icode /= 0) THEN
      WRITE(6,'(A,A)') ' Error returned in CHANGE_DECOMPOSITION',              &
          ' before DERV_LAND_FIELD.'
      WRITE(6,'(A,I5)') ' Error code ',icode
      WRITE(cmessage,'(A)') 'UM_SHELL : Error in CHANGE_DECOMPOSITION.'
      GO TO 9999   !  Exit
    END IF

!     For MPP jobs, calculate the no of land-points on each PE.
! DEPENDS ON: derv_land_field
    CALL derv_land_field (icode,cmessage)
    IF (icode >  0) THEN
      WRITE(6,'(A)') 'Error returned from DERV_LAND_FIELD.'
      WRITE(6,'(A,I5)') 'Error code ',icode
      GO TO 9999   !  Exit
    END IF
        
! Derive lengths involved with output boundary files - atmos.
! DEPENDS ON: derv_intf_a
    CALL derv_intf_a (max_intf_model_levels,max_lbcrow_length,max_lbcrows,     &
                      n_intf_a)
    
    IF (ICODE >  0) GO TO 9999

    WRITE(6,'(A)') '********************************************'//            &
        '***********************************'
!-----------------------------------------------------------------------
! 1.2 Call STASH_PROC: top level control routine for processing of
!                      STASH requests and STASH addressing.

! Open STASHmaster file(s) and count number of records
!   This number is assigned to ppxRecs and used to dynamically
!   allocate the PPX_ arrays in which stash master records are held
    ppxRecs = 1
    icode   = 0
    IF (internal_model_index(a_im) >  0)                                       &
! DEPENDS ON: hdppxrf
        CALL hdppxrf                                                           &
        (nftppxref,'STASHmaster_A',icode,cmessage)
    IF (icode /= 0) GO TO 9999
! Add number of user stash records
! DEPENDS ON: hdppxrf
    CALL hdppxrf(0,'             ',icode,cmessage)

    IF (icode  <   0) THEN
      IF (mype  ==  0) THEN
        WRITE(0,'(A)') 'WARNING : Problem in STASHmaster file(s)'
        WRITE(0,'(A)') '        ',TRIM(cmessage)
      END IF
    ELSE IF (icode  >   0) THEN
      IF (mype  ==  0) THEN
        WRITE(0,'(A)') 'ERROR : Problem in STASHmaster files(s)'
        WRITE(0,'(A)') '      ',TRIM(cmessage)
      END IF
      GO TO 9999  ! Always abort on fatal error.
    END IF

! DEPENDS ON: stash_proc
    CALL stash_proc(nftppxref,nftstmstu,.FALSE.,                               &
        icode,cmessage  )
    IF (icode >  0) GO TO 9999

! Total number of entries (N_PPXRECS) in STASH-addresses array IN_S has
!  obtained by WSTLST in STASH_PROC. Reset ppxRecs to equal this value.
! This is used to dynamically
!  allocate the ppx look-up arrays PPXI, PPXC in U_MODEL/UM_MODEL_4A.

    ppxRecs = n_ppxrecs

!  Set H_SECT to values found in ATMOS_SR in COMDECK MODEL
!  Needs tidied at a later version to be done in a less messy way.
! DEPENDS ON: set_h_sect
    CALL set_h_sect(h_sect,maxsects)

!L----------------------------------------------------------------------
!L 1.3 Calculate addresses of super arrays passed down for dynamic
!L     allocation.
!L
    icode=0
! DEPENDS ON: um_index
    CALL um_index(                                                             &
! ARGSZSP super arrays lengths not dependent on sub-models
     &  SPD1_LEN,SPSTS_LEN,SPBND_LEN,                                   &
! ARGSZSP end
! ARGSZSPA super arrays lengths (atmosphere)
        A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
        A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
! ARGSZSPC super arrays lengths (atmosphere-ocean coupled)
        AO_SPCPL_LEN                                                   &
! ARGSZSPC end
        ,icode,cmessage)

    IF (icode >  0) GO TO 9999

!L----------------------------------------------------------------------
!L 1.5 Set up geometry of the model on the IO servers if needed
!L
!L * Note there is a global barrier (INCLUDING THE IO SERVERS) here *
    IF (isUsingAsyncStash() .OR. isUsingAsyncDumps()) THEN
          ! We need to tell the submodels about our geometry
      CALL IOS_Client_Geometry_Init(decomp_standard_atmos)
    END IF

!L----------------------------------------------------------------------
!L 2. Call U_MODEL/U_MODEL_4A master routine to allocate the main data arrays
!L    and do the calculations.
!L
    IF ( l_endgame ) THEN
! DEPENDS ON: u_model_4A
      CALL u_model_4A (                                               &
          nftppxref,nftstmstu,                                        &
! ARGSZSP super arrays lengths not dependent on sub-models
     &  SPD1_LEN,SPSTS_LEN,SPBND_LEN,                                   &
! ARGSZSP end
! ARGSZSPA super arrays lengths (atmosphere)
        A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
        A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
! ARGSZSPC super arrays lengths (atmosphere-ocean coupled)
        AO_SPCPL_LEN                                                   &
! ARGSZSPC end
         )  
    ELSE
! DEPENDS ON: u_model
      CALL u_model (                                                  &
          nftppxref,nftstmstu,                                        &
! ARGSZSP super arrays lengths not dependent on sub-models
     &  SPD1_LEN,SPSTS_LEN,SPBND_LEN,                                   &
! ARGSZSP end
! ARGSZSPA super arrays lengths (atmosphere)
        A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
        A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
! ARGSZSPC super arrays lengths (atmosphere-ocean coupled)
        AO_SPCPL_LEN                                                   &
! ARGSZSPC end
         ) 
    END IF




! Make sure all the Files are properly Closed. We call model_file_close not
! file_close to ensure that any open files which are PP and hence have
! cached lookups, have the lookups commited to disk before the close.
    DO i=20, nunits
      IF (model_file_managed(i)) THEN
        WRITE(6,'(A,I3,A)')'Managed unit ',i,' is open at end of run, closing'
        CALL model_file_close(i, "dummy_name")
      ELSE IF (is_unit_open(i)) THEN
        WRITE(6,'(A,I3,A)')'Unmanaged unit ',i,' is open at end of run, closing'
        CALL file_close(i, "dummy_name")
      ELSE
        WRITE(6,'(A,I3,A)')'Unit ',i,' is closed'
      END IF
    END DO
!
! Close IO Server
!
    IF (L_IOS_Active()) THEN
      IF (.NOT. isIOServer) THEN

! Shut down IO services on rank 0
        IF (mype==0)                                                           &
            CALL IOS_Shutdown()

! Shut down async stash on all ranks
            IF (isUsingAsyncStash().OR.isUsingAsyncDumps()) &
                CALL IOS_stash_client_fini()

      END IF
    END IF

  END IF ! (am an io server, atmos and io rejoin here)
9999  CONTINUE
  CALL ereport_finalise( )

  CLOSE(5)

  CALL ioShutdown()

  IF (global_rank == 0) THEN
    CALL date_and_time(ch_date2, ch_time2)
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '********************************************'//&
                   '***********************************'
    WRITE(6,'(21A)')                                               &
      '***************** End of UM RUN Job : ',                    &
      ch_time2(1:2),  ':',                                         &
      ch_time2(3:4),  ':',                                         &
      ch_time2(5:6),  ' on ',                                      &
      ch_date2(7:8),  '/',                                         & 
      ch_date2(5:6),  '/',                                         &
      ch_date2(1:4),  ' ******************'
    WRITE(6,'(A,/)') '*****************************************'// &
                   '**************************************'
  END IF
  
  IF (ICODE /= 0) THEN

    CALL Ereport(RoutineName,ICODE,Cmessage)
  END IF

! Time a barrier to ensure dr_hook sees the time for the IOS and atmos to sync up
  IF (lhook) CALL dr_hook('UM_SHELL:RENDEZ-VOUS',zhook_in,zhook_rendez_vous)
  CALL MPL_BARRIER(global_comm,icode)
  IF (lhook) CALL dr_hook('UM_SHELL:RENDEZ-VOUS',zhook_out,zhook_rendez_vous)

! DEPENDS ON: timer
  CALL timer(RoutineName,2)

  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)


      ! Only bother calling GC_EXIT  when PRISM isnt going to shut
      ! things down for this component.
! Close down parallel process communication
  CALL gc_exit()

!     IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

  WRITE(6,'(A,I7,A)') 'Process ',global_rank,' has exited.'

! Close STASH filename unit
  CLOSE(200)

  RETURN

END SUBROUTINE um_shell
