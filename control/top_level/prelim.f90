! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+Construct preliminary STASH list of user requests

! Subroutine Interface:

SUBROUTINE prelim(nrecs,                 &
  ntimes,nlevels,errorstatus,cmessage)

USE rad_input_mod, ONLY: a_lw_radstep_diag, a_sw_radstep_diag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParParams
USE um_input_control_mod, ONLY : lcal360,          &
     phenol_period, triffid_period 
USE cv_param_mod,       ONLY: a_conv_step
USE river_inputs_mod, ONLY: river_step
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE Submodel_Mod
USE version_mod, ONLY: nprofdp, nlevlstsp, npslevp
USE stextend_mod, ONLY: llistty, npos_ts, list_s, itim_s,       &
                        rlevlst_s, levlst_s, lenplst

USE cppxref_mod, ONLY:                                          &
    ppx_version_mask, ppx_grid_type, ppx_lev_flag,              &
    ppx_timavail_code, ppx_space_code, ppx_opt_code,            &
    ppx_lbvc_code, ppx_ptr_code, ppx_item_number,               &
    ppx_lb_code, ppx_lt_code, ppx_lv_code,                      &
    ppx_pf_code, ppx_pl_code, ppx_pt_code

! version_mod items required by cstash.h and stextend.h
USE version_mod, ONLY:                                          &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,         &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,       &
    nlevp, npslevp, npslistp, outfile_s, outfile_e
USE nlstcall_mod, ONLY : model_basis_time, &
                         lmean

USE chsunits_mod, ONLY : nunits

IMPLICIT NONE

!  Description:
!  Constructs a preliminary STASH list of user requests. Uses interim
!  pointer system, by means of the "extra entry" NELEMP+1 in the LIST_S
!  array. At this stage, the input levels encompass all possible levels.
!  Called by STPROC.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 8.2

!  Global variables:
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

! Subroutine arguments


!   Scalar arguments with intent(out):

INTEGER nrecs
INTEGER ntimes
INTEGER nlevels ! Total no. of sets of levs for diags (inpt+outp)
CHARACTER(LEN=80) cmessage

! ErrorStatus:
INTEGER errorstatus

! Local scalars:
LOGICAL      model_lev
LOGICAL      lmask
LOGICAL      lmean_pp
LOGICAL      loffset
INTEGER      totimp
INTEGER      i
INTEGER      ibot1
INTEGER      idiag
INTEGER      idomlev
INTEGER      idom_l
LOGICAL      ldum
INTEGER      ifirst
INTEGER      ifirst1
INTEGER      ilast
INTEGER      ilast1
INTEGER      im
INTEGER      imd
INTEGER      iplof
INTEGER      modl_l
INTEGER      isec_l
INTEGER      item_l
INTEGER      itim_l
INTEGER      itim
INTEGER      itop1
INTEGER      iuse_l
INTEGER      ix1
INTEGER      ix2
INTEGER      iy1
INTEGER      iy2
INTEGER      jlev
INTEGER      lev_offset
INTEGER      lbvc
INTEGER      imax          ! to find max of times-table
INTEGER      itimlst       ! column of times-table
INTEGER      item_chk
INTEGER      river_step_ts ! river routing period in timesteps
INTEGER      start_days,start_secs
INTEGER      end_days,end_secs
INTEGER      ref_days,ref_secs
INTEGER      period

CHARACTER (len=256)          :: cmessage2
CHARACTER (len=*), PARAMETER :: routinename='PRELIM'

! Function and subroutine calls:
LOGICAL  disct_lev
INTEGER  exppxi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header ------------------------------------------------------

! river_step is in seconds - calculate in timesteps

IF (lhook) CALL dr_hook('PRELIM',zhook_in,zhook_handle)
river_step_ts = river_step*(real(steps_per_periodim(a_im))           &
                  /real(secs_per_periodim(a_im)))

! 0.1  Store output-times tables in array ITIM_S

IF (ntimes == 0) THEN
  DO i = 1, nproftp
    IF (iopt_t(i) == 2 .AND. modl_t(i) >  0) THEN

! Profile has output times list
!  MODL_T(I) labels internal model for times list
      DO itim = 1, itim_t(i)

! DEPENDS ON: totimp
        itim_s(itim,i) = totimp(iser_t(itim,i), unt3_t(i), modl_t(i))
        IF (itim_s(itim,i)  ==  -999) THEN
          errorstatus = 100
          WRITE (cmessage,'(a,a,i3)')                                          &
            'PRELIM:TOTIMP:Error in time period conversion',                   &
            ' output times table no.=', i
          WRITE(6,*) cmessage
          GO TO 9999
        END IF
      END DO
      itim_s(itim_t(i)+1,i) = -1
    ELSE
      itim_s(1,i) = -1
    END IF
  END DO
  ntimes = nproftp
END IF

! 0.2  Store output levels lists in array LEVLST_S

lev_offset = nlevels ! Initialised to 0 before entering this routine

! Loop over domain profiles in STASH basis file
DO i = 1, ndprof
  IF (levb_d(i) == -1) THEN

! There is a levels list for this dom prof
    IF (iopl_d(i) == 1 .OR. iopl_d(i) == 2 .OR.                      &
                            iopl_d(i) == 6    ) THEN

! Levs list contains model levs - list type is integer
      llistty(i+lev_offset) = 'I'
    ELSE

! Not model levs - list type real
      llistty(i+lev_offset) = 'R'
    END IF

! LEVT_D(I) = no. of levs in list 'I'
    levlst_s(1, i+lev_offset) = levt_d(i)

! Levels list 'I' was read into (R)LEVLST_D(J,I), J=1,LEVT_D(I),
!  by RDBASIS.
!  Transfer this levels list to (R)LEVLST_S(J,I+LEV_OFFSET),
!  J=2,LEVT_D(I)+1.

    DO jlev = 1, levt_d(i)
      IF (iopl_d(i) == 1 .OR. iopl_d(i) == 2 .OR.                    &
                              iopl_d(i) == 6    ) THEN

!         Model levels
        levlst_s(jlev+1, i+lev_offset) = levlst_d(jlev, i)
      ELSE IF (iopl_d(i) /= 5) THEN

!         Real levels
        rlevlst_s(jlev+1, i+lev_offset) = rlevlst_d(jlev, i)
      END IF
    END DO

    iplof = i + lev_offset

!   Sort this levels list into correct order (if not already in order)
! DEPENDS ON: levsrt
    CALL levsrt( llistty(   iplof),  levlst_s(1, iplof),             &
                levlst_s(2, iplof), rlevlst_s(2, iplof))
  ELSE

! No levels list, i.e., the output from this diag. is on a
!    contiguous range of model levels
    levlst_s(1, i+lev_offset) = 0
    rlevlst_s(1, i+lev_offset) = 0
  END IF
END DO  !  Domain profiles

nlevels = ndprof+lev_offset  ! NDPROF = no. of sets of input levels

IF (nlevels > nlevlstsp) THEN
  WRITE(6,*)                                               &
    'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
  cmessage=                                                &
    'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
  GO TO 9999
END IF

! Section 1. MAIN LOOP - loop over diag requests in STASH basis file

IF (ndiag > 0) THEN

  DO idiag = 1, ndiag

    modl_l = modl_b(idiag)
    isec_l = isec_b(idiag)
    item_l = item_b(idiag)
    idom_l = idom_b(idiag)
    iuse_l = iuse_b(idiag)
    itim_l = itim_b(idiag)

    item_chk = 0

! DEPENDS ON: exppxi
    item_chk = exppxi(modl_l, isec_l, item_l, ppx_item_number,       &
                      errorstatus, cmessage)
    IF (item_chk /= item_l) THEN
      WRITE(cmessage2,*) 'Diagnostic discarded ',                    &
            ' model ', modl_l, ' section ', isec_l, ' item', item_l, &
            ' No stashmaster record'
      errorstatus = -10

      CALL ereport(routinename, errorstatus, cmessage2)

! Make this item null
      item_b(idiag) = 0
      GO TO 999
    END IF
    IF (itim_l /= 0) THEN       ! If the diag is not a null request

! Section 1.0  Extract data required for STASH processing from PPXI

      IF (nrecs == nrecdp) THEN
        WRITE(cmessage2,*)                                           &
          'TOO MANY STASH LIST ENTRIES, REQUEST DENIED',             &
          ' (M,S,I)',                                                &
          modl_l, isec_l, item_l

        errorstatus = -20

        CALL ereport(routinename, errorstatus, cmessage2)
        GO TO 999
      END IF

! DEPENDS ON: exppxi
      vmsk    = exppxi(modl_l ,isec_l ,item_l,ppx_version_mask ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      ispace  = exppxi(modl_l ,isec_l ,item_l,ppx_space_code   ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      itima   = exppxi(modl_l ,isec_l ,item_l,ppx_timavail_code,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      igp     = exppxi(modl_l ,isec_l ,item_l,ppx_grid_type    ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      ilev    = exppxi(modl_l ,isec_l ,item_l,ppx_lv_code      ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      ibot    = exppxi(modl_l ,isec_l ,item_l,ppx_lb_code      ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      itop    = exppxi(modl_l ,isec_l ,item_l,ppx_lt_code      ,     &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      iflag   = exppxi(modl_l ,isec_l ,item_l,ppx_lev_flag     ,     &
                       errorstatus, cmessage)
      DO i = 1, 6

! DEPENDS ON: exppxi
        iopn(i) = exppxi(modl_l ,isec_l ,item_l, ppx_opt_code+i-1 ,  &
                         errorstatus, cmessage)
      END DO
! DEPENDS ON: exppxi
      ipseudo = exppxi(modl_l ,isec_l ,item_l,ppx_pt_code         ,  &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      ipfirst = exppxi(modl_l ,isec_l ,item_l,ppx_pf_code         ,  &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      iplast  = exppxi(modl_l ,isec_l ,item_l,ppx_pl_code         ,  &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      ptr_prog= exppxi(modl_l ,isec_l ,item_l,ppx_ptr_code        ,  &
                       errorstatus, cmessage)
! DEPENDS ON: exppxi
      lbvc    = exppxi(modl_l ,isec_l ,item_l,ppx_lbvc_code       ,  &
                       errorstatus, cmessage)

! Check availability of diagnostic
! DEPENDS ON: tstmsk
      CALL tstmsk(modl_l, isec_l, lmask, ldum, errorstatus, cmessage)
      IF (.NOT. lmask) THEN
        WRITE(cmessage2,*)                                           &
          'DIAGNOSTIC NOT AVAILABLE TO THIS VERSION ',               &
          'REQUEST DENIED ',                                         &
          '(M,S,I)',                                                 &
          modl_l, isec_l, item_l

        errorstatus = -30

        CALL ereport(routinename, errorstatus, cmessage2)
        GO TO 999
      END IF

      nrecs = nrecs+1

      list_s(st_model_code  , nrecs) = modl_l
      list_s(st_sect_no_code, nrecs) = isec_l
      list_s(st_item_code   , nrecs) = item_l

! Prelim pointer for 'child' records
      list_s(nelemp+1       , nrecs) = nrecs
      list_s(st_lookup_ptr  , nrecs) = -1

! Set input code for STASH requests:
!  =0 Use primary or secondary field:       D1(SI(item,section,model))
!  =1 Use field in diagnostic space: STASHwork(SI(item,section,model))
!  =-j Use diagnostic at D1(LIST_S(st_output_addr,j))
      IF ( (ispace == 2) .OR. (ispace == 4)                                &
        .OR. (ispace == 7) .OR. (ispace == 8) .OR. (ispace == 9)) THEN
        list_s(st_input_code, nrecs) = 0
      ELSE
        list_s(st_input_code, nrecs) = 1
      END IF

      IF ( (itima >= 5) .AND. (itima <= 12) ) THEN
        lmean_pp = .TRUE.
      ELSE
        lmean_pp = .FALSE.
      END IF


! 1.1   Expand the domain profile ---------------------------

!   Averaging and Weighting
      im = imsk_d(idom_l)
      IF ( (igp ==  2) .OR. (igp ==  3)    .OR.                      &
           (igp == 12) .OR. (igp == 13))   THEN

! Diags only available over land/sea
        IF ( (imsk_d(idom_l)  ==  1)         .AND.                   &
             (igp == 3 .OR. igp == 13) )     THEN

! Diag requested over land+sea, only available over sea
          im = 3
        ELSE IF ( (imsk_d(idom_l)  ==  1)     .AND.                  &
                  (igp == 2 .OR. igp == 12) ) THEN

! Diag requested over land+sea, only available over land
          im = 2
        ELSE IF ( (imsk_d(idom_l)  ==  2)    .AND.                    &
                (igp == 3 .OR. igp == 13) ) THEN

! Diag requested over land, only available over sea
          WRITE(6,*) 'PRELIM: CHANGED TO SEA DIAG'
          WRITE(6,*) 'MODEL,SECTION,ITEM ',                          &
                     modl_l, isec_l, item_l
          im = 3
        ELSE IF ( (imsk_d(idom_l)  ==  3)    .AND.                   &
                  (igp == 2 .OR. igp == 12) ) THEN

! Diag requested over sea, only available over land
          WRITE(6,*) 'PRELIM: CHANGED TO LAND DIAG'
          WRITE(6,*) 'MODEL,SECTION,ITEM ',                          &
                     modl_l, isec_l, item_l
          im = 2
        END IF
      END IF

      list_s(st_gridpoint_code, nrecs) = im+10*imn_d(idom_l)
      list_s(st_weight_code   , nrecs) =       iwt_d(idom_l)

!   Horizontal area
!    - convert lat/long spec to row/column numbers if appropriate;
!    - convert lat/long spec to equatorial lat/long if appropriate.
      IF (iopa_d(idom_l) == 1) THEN

! Full domain
! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, -90, 0, 360,                                    &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 2 ) THEN

! N Hemis
! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, 0, 0, 360,                                      &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 3 ) THEN

! S Hemis
! DEPENDS ON: lltorc
        CALL lltorc(igp, 0, -90, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 4 ) THEN

! 90N-30N
! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, 30, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 5 ) THEN

! 30S-90S
! DEPENDS ON: lltorc
        CALL lltorc(igp, -30, -90, 0, 360,                                   &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 6 ) THEN

! 30N-00N
! DEPENDS ON: lltorc
        CALL lltorc(igp, 30, 00, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 7 ) THEN

! 00S-30S
! DEPENDS ON: lltorc
        CALL lltorc(igp, 00, -30, 0, 360,                                    &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 8 ) THEN

! 30N-30S
! DEPENDS ON: lltorc
        CALL lltorc(igp, 30, -30, 0, 360,                                    &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs), list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 9 ) THEN

! Other lat/long spec
! DEPENDS ON: lltorc
        CALL lltorc(igp, inth_d(idom_l), isth_d(idom_l),                     &
              iwst_d(idom_l), iest_d(idom_l),                                &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),    &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == 10) THEN

! Grid point spec
! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, -90, 0, 360, iy1, iy2, ix1, ix2)
        list_s(st_north_code, nrecs) = MIN(inth_d(idom_l), iy2)
        list_s(st_south_code, nrecs) = MIN(isth_d(idom_l), iy2)
        list_s(st_west_code , nrecs) = MIN(iwst_d(idom_l), ix2)
        list_s(st_east_code , nrecs) = MIN(iest_d(idom_l), ix2)
      ELSE
        WRITE(cmessage2,*) 'INVALID DOMAIN AREA OPTION=',            &
                           iopa_d(idom_l),                           &
                           '(M,S,I)',                                &
                           modl_l, isec_l, item_l

        errorstatus = -35

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF

! Input level setting
! DEPENDS ON: disct_lev
      model_lev = disct_lev(ilev, errorstatus, cmessage)
      IF (model_lev) THEN

! Model levels
! Set bottom level
! DEPENDS ON: levcod
        CALL levcod(ibot, ibot1, errorstatus, cmessage)

! Set top level
! DEPENDS ON: levcod
        CALL levcod(itop, itop1, errorstatus, cmessage)

! Contig. range of model levels
        IF (iflag == 0) THEN
          list_s(st_input_bottom, nrecs) = ibot1
          list_s(st_input_top   , nrecs) = itop1

! Non-contig. levels list
        ELSE IF (iflag == 1) THEN
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  1
        END IF
      ELSE

! Non-model levels
        IF (ilev == 3) THEN

!  Pressure levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  2
        ELSE IF (ilev == 4) THEN

!  Height levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  3
        ELSE IF (ilev == 5) THEN

!  Special levels
          list_s(st_input_bottom, nrecs) = 100
          list_s(st_input_top   , nrecs) = lbvc
        ELSE IF (ilev == 7) THEN

!  Theta levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  4
        ELSE IF (ilev == 8) THEN

!  PV levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  5
        ELSE IF (ilev == 9) THEN

!  Cloud threshold levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  6
        END IF
      END IF

! Output level specification
! DEPENDS ON: disct_lev
      model_lev = disct_lev(ilev, errorstatus, cmessage)
      IF (model_lev) THEN

! Model levels
        IF (levb_d(idom_l) >= 0) THEN

! Contiguous range of model levels
          list_s(st_output_bottom, nrecs) = MAX(levb_d(idom_l), ibot1)
          list_s(st_output_top   , nrecs) = MIN(levt_d(idom_l), itop1)
          IF ( (levb_d(idom_l) <  ibot1)   .OR.                              &
               (levt_d(idom_l) >  itop1) ) THEN
            WRITE(cmessage2,*)                                               &
              'DIAGNOSTIC HAS LEVEL RANGE OUT OF BOUNDS; CORRECTED ',        &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -40

            CALL ereport(routinename, errorstatus, cmessage2)
          END IF
          IF ( (  ts_d(idom_l) ==   'Y')      .AND.                          &
               ( (levb_d(idom_l) <  ibot1)    .OR.                           &
                 (levt_d(idom_l) >  itop1) )  ) THEN
            WRITE(cmessage2,*)                                               &
              'TIME SERIES DOMAIN',                                          &
              'HAS INCONSISTENT LEVELS; DIAGNOSTIC IGNORED',                 &
              ' (M,S,I) ',                                                   &
              modl_l, isec_l, item_l

            errorstatus = -50

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs-1
            GO TO 999
          END IF
          IF ( (levt_d(idom_l) <  ibot1)  .OR.                               &
               (levb_d(idom_l) >  itop1)) THEN
            WRITE(cmessage2,*)                                               &
              'DIAGNOSTIC HAS TOP/BOT LEVELS INCONSISTENT; DIAG IGNORED',    &
              ' (M,S,I) ',                                                   &
              modl_l, isec_l, item_l

            errorstatus = -60

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs-1
            GO TO 999
          END IF
        ELSE

! Non-contig. list of model levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 1
        END IF
      ELSE

! Non-model levels
        IF (ilev == 5) THEN

! Special level
          list_s(st_output_bottom, nrecs) = 100
          list_s(st_output_top   , nrecs) = lbvc
        ELSE IF (ilev == 3) THEN

! Pressure levels
          list_s(st_output_bottom,nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   ,nrecs) = 2
        ELSE IF (ilev == 4) THEN

! Height levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 3
        ELSE IF (ilev == 7 ) THEN

! Theta levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 4
        ELSE IF (ilev == 8 ) THEN

! PV levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 5
        ELSE IF (ilev == 9 ) THEN

! Cloud threshold levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 6
        ELSE
          WRITE(cmessage2,*) 'DOMAIN LEVEL OPTION=',IOPL_D(IDOM_L),          &
            'DIAG IGNORED. (M,S,I) ',                                        &
            modl_l, isec_l, item_l

          errorstatus=-70

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
      END IF

! Output pseudo-levels level setting
      IF (ipseudo /= plt_d(idom_l)) THEN
        WRITE(cmessage2,*)                                                   &
          'DIAGNOSTIC HAS ',                                                 &
          'INVALID PSEUDO LEVEL TYPE; IGNORED.',                             &
          ' (M,S,I)',                                                        &
          modl_l, isec_l, item_l

        errorstatus = -80

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF
      list_s(st_pseudo_in,nrecs) = 0  !(This is set in INPUTL)
      IF (ipseudo >  0) THEN

! Pseudo levels list for this diagnostic
        list_s(st_pseudo_out, nrecs) = plpos_d(idom_l)
        lenplst(plpos_d(idom_l))     = pllen_d(idom_l)
        ifirst = pslist_d(1, plpos_d(idom_l))
        ilast  = pslist_d(pllen_d(idom_l), plpos_d(idom_l))

! Check pseudo level limits
! DEPENDS ON: pslims
        CALL pslims(ipfirst, iplast, ifirst1, ilast1)
        IF (ifirst <  ifirst1) THEN
          WRITE(cmessage2,*)                                                 &
            'DIAGNOSTIC HAS ',                                               &
            'FIRST PSEUDO LEVEL TOO LOW; IGNORED.',                          &
            '(M,S,I) ',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -90

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
        IF (ilast >  ilast1) THEN
          WRITE(cmessage2,*)                                                 &
            'DIAGNOSTIC HAS ',                                               &
            'LAST PSEUDO LEVEL TOO HIGH; IGNORED',                           &
            ' (M,S,I)',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -95

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
      ELSE
        list_s(st_pseudo_out,nrecs) = 0
      END IF

! Time-series domain profiles
      IF (ts_d(idom_l) == 'Y') THEN

! Pointer for location of time series
        list_s(st_series_ptr, nrecs) = npos_ts(idom_l)
      ELSE
        list_s(st_series_ptr, nrecs) = 0
      END IF

! 1.2   Expand the useage profile --------------------------

      IF (locn_u(iuse_l) == 5) THEN                    ! PP file

        IF (lmean_pp) THEN
          list_s(st_output_code, nrecs) = -27
          list_s(st_macrotag,    nrecs) =  0
        ELSE
          WRITE(6,*)                                                           &
            'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST HAS ',            &
            'OUTPUT DESTINATION CODE 5 (CLIMATE MEAN PP FILE) ',               &
            'BUT DIAGNOSTIC IS NOT A CLIMATE MEAN; REQUEST IGNORED'

          WRITE(cmessage2,*)                                                   &
            'DIAGNOSTIC IS NOT CLIMATE MEAN, BUT SENT TO MEAN FILE:',          &
            'IGNORED. (M,S,I) ',                                               &
            modl_l, isec_l, item_l

          errorstatus = -100

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF

      ELSE IF (lmean_pp) THEN

        WRITE(6,*)                                                             &
          'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST IS A ',             &
          'CLIMATE MEAN - SHOULD HAVE OUTPUT DESTINATION CODE 5 ',             &
          '(CLIMATE MEAN PP FILE); REQUEST IGNORED'

        WRITE (cmessage2,*)                                                    &
          'DIAGNOSTIC IS CLIMATE MEAN, INCORRECT DESTINATION:',                &
          'IGNORED. (M,S,I) ',                                                 &
          modl_l, isec_l, item_l

        errorstatus = -110

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999

      ELSE IF (locn_u(iuse_l) == 3) THEN                ! PP file

        list_s(st_output_code, nrecs) = -iunt_u(iuse_l)
        list_s(st_macrotag,    nrecs) = 0

      ELSE IF (locn_u(iuse_l) == 1) THEN ! Dump store: set user tag

        list_s(st_output_code, nrecs) = 1
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l)

      ELSE IF (locn_u(iuse_l) == 6) THEN ! Secondary dump store:
                                         !             set user tag
        list_s(st_output_code, nrecs) = 2
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l)

      ELSE IF (locn_u(iuse_l) == 2) THEN ! Climate mean: tag set
                                         !   1000*(time mean tag)

        list_s(st_output_code, nrecs) = 1
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l) * 1000

      ELSE IF (locn_u(iuse_l) == 4) THEN ! Printed output

        list_s(st_output_code, nrecs) = 7
        list_s(st_macrotag,    nrecs) = 0

      ELSE

        WRITE(cmessage2,*)                                                     &
          'INVALID USEAGE OPTION=', LOCN_U(IUSE_L),                            &
          ': DIAGNOSTIC IGNORED. (M,S,I) ',                                    &
          modl_l, isec_l, item_l

        errorstatus = -120

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999

      END IF

! 1.3   Expand the time profile ------------------------------

! Initialise as single time field

!   Set time processing record

      IF (lmean_pp) THEN
        IF (ityp_t(itim_l) /= 1) THEN
          WRITE(cmessage2,*)                                                   &
            'CLIMATE MEANS MUST NOT BE TIME PROCESSED.',                       &
            '(M,S,I) ',                                                        &
            modl_l, isec_l, item_l

          errorstatus = -130

          CALL ereport(routinename, errorstatus, cmessage2)
        END IF
        list_s(st_proc_no_code, nrecs) = 1
      ELSE
        list_s(st_proc_no_code, nrecs) = ityp_t(itim_l)
      END IF

! Initialise offset to 0
      list_s(st_offset_code,nrecs)=0

!   Set period record
      IF (ityp_t(itim_l) == 1 .OR. lmean_pp) THEN        ! No period
        list_s(st_period_code, nrecs) = 0
      ELSE IF ( (intv_t(itim_l) == -1)  .AND.                                  &
                (ityp_t(itim_l) == 2) ) THEN
        list_s(st_period_code, nrecs) = -1
      ELSE
        list_s(st_period_code, nrecs) =                                        &
! DEPENDS ON: totimp
          totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
        IF (list_s(st_period_code,nrecs)  ==  -999) THEN
          errorstatus = 101
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                                &
            'PRELIM:TOTIMP:Error in time period conversion',                   &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          WRITE(6,*) cmessage
          GO TO 9999
        END IF
      END IF

      IF ( lmean_pp .AND. (iopt_t(itim_l) /= 1) ) THEN
        WRITE(cmessage2,*)                                                     &
          'CLIMATE MEANS MUST USE STANDARD FREQUENCY. DIAG IGNORED.',          &
          '(M,S,I) ',                                                          &
          modl_l, isec_l, item_l

        errorstatus = -140

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF

      IF (iopt_t(itim_l) == 1 .OR. iopt_t(itim_l) == 3) THEN

! Regular output times
        list_s(st_freq_code, nrecs)=                                           &
! DEPENDS ON: totimp
          totimp(ifre_t(itim_l), unt3_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 102
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                                &
            'PRELIM:TOTIMP:Error in time period conversion',                   &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          WRITE(6,*) cmessage
          GO TO 9999
        END IF

! Time profile is specified with a start and end date, between which,
! the diagnostic is written out.
        IF (iopt_t(itim_l) == 3) THEN

! Calculate the difference in (integer number of) hours between the
! start date in the time profile and the model basis time.
! DEPENDS ON: time2sec
          CALL time2sec(isdt_t(1, itim_l), isdt_t(2, itim_l),                  &
                        isdt_t(3, itim_l), isdt_t(4, itim_l),                  &
                        isdt_t(5, itim_l), isdt_t(6, itim_l),                  &
                        0, 0, start_days, start_secs, lcal360)

! DEPENDS ON: time2sec
          CALL time2sec(model_basis_time(1), model_basis_time(2),              &
                        model_basis_time(3), model_basis_time(4),              &
                        model_basis_time(5), model_basis_time(6),              &
                        0, 0, ref_days, ref_secs, lcal360)

          period = (start_days-ref_days) * 24 +                                &
                       INT((start_secs-ref_secs) / 3600.0)

! Convert the start date (hours) into timesteps since basis time.
          list_s(st_start_time_code, nrecs) = totimp(period, 'H ', modl_l)

        ELSE

          list_s(st_start_time_code, nrecs) =                                  &
! DEPENDS ON: totimp
            totimp(istr_t(itim_l), unt3_t(itim_l), modl_l)
        END IF
        IF (list_s(st_start_time_code, nrecs)  ==  -999) THEN
          errorstatus = 103
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                                &
            'PRELIM:TOTIMP:Error in time period conversion',                   &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          WRITE(6,*) cmessage
          GO TO 9999
        END IF
        IF (iend_t(itim_l) == -1) THEN
          list_s(st_end_time_code, nrecs) = -1
        ELSE
          IF (iopt_t(itim_l) == 3) THEN

! Calculate the difference in (integer number of) hours between the
! end date in the time profile and the model basis time.
! DEPENDS ON: time2sec
            CALL time2sec(iedt_t(1, itim_l), iedt_t(2, itim_l),                &
                          iedt_t(3, itim_l), iedt_t(4, itim_l),                &
                          iedt_t(5, itim_l), iedt_t(6, itim_l),                &
                          0, 0, end_days, end_secs, lcal360)

! DEPENDS ON: time2sec
            CALL time2sec(model_basis_time(1), model_basis_time(2),            &
                          model_basis_time(3), model_basis_time(4),            &
                          model_basis_time(5), model_basis_time(6),            &
                          0, 0, ref_days, ref_secs, lcal360)

! End date uses ceiling to calculate the next largest integer secs so that
! the end date specified in the UMUI is inclusive.
            period = (end_days-ref_days) * 24 +                                &
                       CEILING((end_secs-ref_secs) / 3600.0)

! Convert the end date (hours) into timesteps since basis time.
            list_s(st_end_time_code, nrecs)= totimp(period, 'H ', modl_l)
          ELSE

            list_s(st_end_time_code, nrecs)=                                   &
! DEPENDS ON: totimp
              totimp(iend_t(itim_l), unt3_t(itim_l), modl_l)
          END IF
        END IF
        IF (list_s(st_end_time_code, nrecs)  ==  -999) THEN
          errorstatus = 104
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                                &
            'PRELIM:TOTIMP:Error in time period conversion',                   &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          WRITE(6,*) cmessage
          GO TO 9999
        END IF

!   Set end time to -1 if output requested to end of run

!   Correct start time for radiation, periodic convection, leaf
!   phenology and vegetation competition and river routing
        IF ( (itima == 2) .AND. (a_lw_radstep_diag /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_lw_radstep_diag )
          list_s(st_start_time_code, nrecs) =                                  &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 3) .AND. (a_sw_radstep_diag /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_sw_radstep_diag )
          list_s(st_start_time_code, nrecs) =                                  &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 13) .AND. (a_conv_step /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_conv_step )
          list_s(st_start_time_code, nrecs) =                                &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 14) .AND. (phenol_period /= 1) ) THEN
          imd = MOD(list_s(st_start_time_code, nrecs), phenol_period)
          list_s(st_start_time_code, nrecs) =                                &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 15) .AND. (triffid_period /= 1) ) THEN
          imd = MOD(list_s(st_start_time_code, nrecs), triffid_period)
          list_s(st_start_time_code, nrecs)=                                 &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 16).AND.(river_step_ts /= 1)) THEN
          imd=mod(list_s(st_start_time_code,nrecs),river_step_ts)
          list_s(st_start_time_code,nrecs)=                            &
            list_s(st_start_time_code,nrecs)-imd
          loffset=.TRUE.
        ELSE
          loffset = .FALSE.
        END IF
      ELSE IF (iopt_t(itim_l) == 2) THEN

! List of specified output times
        list_s(st_freq_code, nrecs) = -itim_l
      ELSE
        WRITE(cmessage2,*)                                                   &
          'INVALID OUTPUT TIMES CODE. DIAG IGNORED.',                        &
          '(M,S,I) ',                                                        &
          modl_l, isec_l, item_l

        errorstatus = -150

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF

      IF (lmean_pp) list_s(st_freq_code, nrecs) = 1

      IF ( (list_s(st_proc_no_code, nrecs) >  1)   .AND.                     &
           (list_s(st_proc_no_code, nrecs) <= 6) ) THEN

! Other than single time field
        IF (nrecs >= nrecdp) THEN
          WRITE(cmessage2,*)                                                 &
            'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                     &
            '(M,S,I) ',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -160

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF

        DO i = 1, nelemp + 1          ! Copy stash list forward
          list_s(i, nrecs+1) = list_s(i, nrecs)
        END DO

        IF (loffset) THEN         ! Rad or conv timesteps,
                                    !       1 already added
          list_s(st_start_time_code, nrecs+1) =                              &
            list_s(st_start_time_code, nrecs+1) - 1
          IF (list_s(st_period_code,nrecs) /= -1) THEN

! Offsets are added to start time
            list_s(st_offset_code, nrecs) =                                  &
! DEPENDS ON: totimp
              totimp(ioff_t(itim_l), unt2_t(itim_l), modl_l)
            IF (list_s(st_offset_code, nrecs)  ==  -999) THEN
              errorstatus = 1
              WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                          &
                'PRELIM:TOTIMP:Error in time period conversion',             &
                ' model=',modl_l, ' section=',isec_l, ' item=',item_l
              WRITE(6,*) cmessage
              GO TO 9999
            END IF
            list_s(st_start_time_code, nrecs) =                              &
              list_s(st_start_time_code, nrecs) -                            &
              list_s(st_period_code, nrecs)     +                            &
              list_s(st_offset_code, nrecs)
          ELSE
            list_s(st_start_time_code,nrecs) = 1
          END IF

        ELSE

          IF (list_s(st_period_code, nrecs) /= -1) THEN

! Offsets are added to start time
            list_s(st_offset_code, nrecs) =                                  &
! DEPENDS ON: totimp
              totimp(ioff_t(itim_l), unt2_t(itim_l), modl_l)
            IF (list_s(st_offset_code, nrecs)  ==  -999) THEN
               errorstatus = 1
              WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                          &
                'PRELIM:TOTIMP:Error in time period conversion',             &
                ' model=',modl_l, ' section=',isec_l, ' item=',item_l
              WRITE(6,*) cmessage
              GO TO 9999
            END IF
            list_s(st_start_time_code, nrecs) =                              &
              list_s(st_start_time_code, nrecs) -                            &
              list_s(st_period_code, nrecs) + 1 +                            &
              list_s(st_offset_code, nrecs)
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        END IF

        IF (list_s(st_start_time_code, nrecs) <  1) THEN
          WRITE(cmessage2,*)                                                 &
            'DIAGNOSTIC START TIME BEFORE PERIOD, SETTING TO 1.',            &
            '(M,S,I) ',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -170

          CALL ereport(routinename, errorstatus, cmessage2)
          list_s(st_start_time_code, nrecs) = 1
        END IF

! Check if offset corresponds to a valid timestep.
! If not, reject the diagnostic.
        SELECT CASE (itima)
        CASE (2)   ! Long-Wave Radiation
          IF (MOD( list_s(st_offset_code, nrecs), a_lw_radstep_diag ) /= 0)  &
          THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH LW_RAD STEP.',                     &
              ' DIAG IGNORED. (M,S,I)',                                      &
              modl_l, isec_l, item_l

            errorstatus = -180

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (3)   ! Short-Wave Radiation
          IF (MOD( list_s(st_offset_code, nrecs), a_sw_radstep_diag ) /= 0)  &
          THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH SW_RAD STEP.',                     &
              ' DIAG IGNORED. (M,S,I)',                                      &
              modl_l, isec_l, item_l

            errorstatus = -190

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        CASE (13)   ! Convection
          IF (MOD( list_s(st_offset_code, nrecs), a_conv_step ) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH CONVECT STEP.',                    &
              ' DIAG IGNORED. (M,S,I)',                                      &
              modl_l, isec_l, item_l

            errorstatus = -200

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (14)   ! Leaf Phenology
          IF (MOD( list_s(st_offset_code, nrecs), phenol_period ) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH PHENOL PERIOD.',                   &
              ' DIAG IGNORED. (M,S,I)',                                      &
              modl_l, isec_l, item_l

            errorstatus = -210

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (15)   ! Triffid
          IF (MOD( list_s(st_offset_code, nrecs), triffid_period ) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH TRIFFID PERIOD.',                  &
              ' DIAG IGNORED. (M,S,I)',                                      &
              modl_l, isec_l, item_l

            errorstatus = -220

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (16)   ! Trip river routing model
          IF (MOD(list_s(st_offset_code,nrecs),river_step_ts)/= 0) THEN
            WRITE(cmessage2,*)                                               &
              'OFFSET DOES NOT AGREE WITH RIVER ROUTING PERIOD.',            &
     &        ' DIAG IGNORED. (M,S,I)',                                      &
     &        MODL_L,ISEC_L,ITEM_L
            errorstatus = -230

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GOTO 999
          END IF
        END SELECT

        list_s(st_proc_no_code ,nrecs+1) = 1

        list_s(st_input_bottom ,nrecs+1) = list_s(st_output_bottom, nrecs)

        list_s(st_input_top    ,nrecs+1) = list_s(st_output_top   , nrecs)

        list_s(st_input_code   ,nrecs+1) = -nrecs
        list_s(st_output_code  ,nrecs)   = 1
        list_s(st_series_ptr   ,nrecs+1) = 0
        list_s(nelemp+1        ,nrecs+1) = nrecs + 1

        list_s(st_freq_code,nrecs)=                                          &
! Frequency
! DEPENDS ON: totimp
          totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 105
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                              &
            'PRELIM:TOTIMP:Error in time period conversion',                 &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          WRITE(6,*) cmessage
          GO TO 9999
        END IF

!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

        IF (itima == 2) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_lw_radstep_diag
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_lw_radstep_diag) /= 0) &
          THEN
            WRITE(cmessage2,*)                                               &
              'INCORRECT SAMPLING FOR LW_RADSTEP. FREQ=',                    &
              list_s(st_freq_code, nrecs), ':IGNORED.',                      &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -225

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 3) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_sw_radstep_diag
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_sw_radstep_diag) /= 0) &
          THEN
            WRITE(cmessage2,*)                                               &
              'INCORRECT SAMPLING FOR SW_RADSTEP. FREQ=',                    &
              list_s(st_freq_code,nrecs), ':IGNORED.',                       &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -230

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 13) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_conv_step
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_conv_step) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'INCORRECT SAMPLING FOR CONV_STEP. FREQ=',                     &
              list_s(st_freq_code, nrecs), ':IGNORED.',                      &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -240

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 14) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = phenol_period
          ELSE IF (MOD(list_s(st_freq_code, nrecs), phenol_period) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'INCORRECT SAMPLING FOR PHENOL_PERIOD. FREQ=',                 &
              list_s(st_freq_code,nrecs), ':IGNORED.',                       &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -250

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 15) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = triffid_period
          ELSE IF                                                            &
            (MOD(list_s(st_freq_code, nrecs), triffid_period) /= 0) THEN
            WRITE(cmessage2,*)                                               &
              'INCORRECT SAMPLING FOR TRIFFID_PERIOD. FREQ=',                &
              list_s(st_freq_code,nrecs), ':IGNORED.',                       &
              '(M,S,I) ',                                                    &
              modl_l, isec_l, item_l

            errorstatus = -260

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF(itima == 16) THEN
          IF (list_s(st_freq_code,nrecs) == 1) THEN
              list_s(st_freq_code,nrecs)=river_step_ts
          ELSE IF                                                            &
     &      (MOD(list_s(st_freq_code,nrecs),river_step_ts) /= 0) THEN
            WRITE(cmessage2,*)                                               &
     &        'INCORRECT SAMPLING FOR RIVER_STEP_TS. FREQ=',                 &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                       &
     &        '(M,S,I) ',                                                    &
     &        modl_l,isec_l,item_l

            errorstatus=-270

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GO TO 999
          END IF
        END IF

! For the NRECS item an end_time_code needs to be set if we
! are dealing with a times table rather than  regular diagn.
! This should be the maximum timestep in the time list. The list
! should be ready sorted (and thus maximum is last member) but
! will run through and find maximum to be on the safe side.

        imax = 0
        itimlst = -list_s(st_freq_code, nrecs+1)

        IF (itimlst  >   0) THEN      ! List *not* regular
          DO i = 1, itim_t(itimlst)
            IF (imax  <   itim_s(i, itimlst)) THEN
              imax = itim_s(i, itimlst)
            END IF
          END DO

          list_s(st_end_time_code, nrecs) = imax

        END IF

!   Period

        IF ((intv_t(itim_l) == -1).AND.(ityp_t(itim_l) == 2)) THEN
          list_s(st_period_code, nrecs) = -1
        ELSE
          list_s(st_period_code, nrecs)=                                     &
! DEPENDS ON: totimp
            totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
          IF (list_s(st_period_code, nrecs)  ==  -999) THEN
            errorstatus = 106
            WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                            &
              'PRELIM:TOTIMP:Error in time period conversion',               &
              ' model=',modl_l, ' section=',isec_l, ' item=',item_l
            WRITE(6,*) cmessage
            GO TO 9999
          END IF
        END IF

!   Add the record - unless the output destination is the dump,
!                      and output at the accumulating period
        IF (    locn_u(iuse_l) >  2                                          &
            .OR.                                                             &
             ( (list_s(st_freq_code  , nrecs+1 )    /=                       &
                list_s(st_period_code, nrecs   ) )                           &
                                                   .AND.                     &
               (list_s(st_start_time_code, nrecs+1) /=                       &
                list_s(st_end_time_code  , nrecs+1) )   )  ) THEN

! No tag for parent
          list_s(st_macrotag, nrecs) = 0
          nrecs = nrecs + 1
        END IF

      ELSE IF (list_s(st_proc_no_code, nrecs) == 8) THEN

! Option of "daily" mean timeseries
        IF (nrecs >= nrecdp) THEN
          WRITE(cmessage2,*)                                                 &
            'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                     &
            '(M,S,I) ',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -270

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs - 1
          GO TO 999
        END IF

! Special case where 2 extra records required
!  Record 1 - time mean only no spatial processing
!  Record 2 - timeseries formed extracting from record 1
!  Record 3 - extract timeseries from dump ie record 2

        DO i = 1, nelemp + 1      ! Copy stash list forward
          list_s(i, nrecs+1) = list_s(i, nrecs)
          list_s(i, nrecs+2) = list_s(i, nrecs)
        END DO

        IF (loffset) THEN         ! Rad or conv timesteps,
                                  !  1 already added
          list_s(st_start_time_code, nrecs+2) =                              &
            list_s(st_start_time_code, nrecs+2) - 1
          IF (list_s(st_period_code, nrecs) /= -1) THEN
            list_s(st_start_time_code, nrecs) =                              &
              list_s(st_start_time_code, nrecs) -                            &
              list_s(st_period_code, nrecs)
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        ELSE

          IF (list_s(st_period_code, nrecs) /= -1) THEN
            list_s(st_start_time_code, nrecs) =                              &
              list_s(st_start_time_code, nrecs) -                            &
              list_s(st_period_code, nrecs) + 1
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        END IF

        IF (list_s(st_start_time_code, nrecs) <  1) THEN
          WRITE(cmessage2,*)                                                 &
            'START TIME BEFORE PERIOD, SETTING TO 1',                        &
            '(M,S,I) ',                                                      &
            modl_l, isec_l, item_l

          errorstatus = -280

          CALL ereport(routinename, errorstatus, cmessage2)
          list_s(st_start_time_code, nrecs) = 1
        END IF

        list_s(st_proc_no_code ,nrecs)   = 3  ! time mean
        list_s(st_proc_no_code ,nrecs+1) = 8  ! timseries special case
        list_s(st_proc_no_code ,nrecs+2) = 1  !  extract

! Reset first record to no area weight or spatial processing
! ie first record just controls time meaning

        list_s(st_gridpoint_code, nrecs) = 1
        list_s(st_weight_code, nrecs) = 0

        list_s(st_input_bottom ,nrecs+1) =                                   &
          list_s(st_output_bottom, nrecs)
        list_s(st_input_bottom ,nrecs+2) =                                   &
          list_s(st_output_bottom, nrecs+1)

        list_s(st_input_top    ,nrecs+1) =                                   &
          list_s(st_output_top   , nrecs)
        list_s(st_input_top    ,nrecs+2) =                                   &
          list_s(st_output_top   , nrecs+1)

        list_s(st_input_code   ,nrecs+1) = -nrecs
        list_s(st_input_code   ,nrecs+2) = -nrecs - 1
        list_s(st_output_code  ,nrecs  ) = 1
        list_s(st_output_code  ,nrecs+1) = 1  ! dump
        list_s(st_series_ptr   ,nrecs+2) = 0
        list_s(st_series_ptr   ,nrecs)   = 0
        list_s(nelemp+1        ,nrecs+1) = nrecs + 1
        list_s(nelemp+1        ,nrecs+2) = nrecs + 2

!  definition 8 implies frequency of time mean over every timestep
        list_s(st_freq_code, nrecs) = 1
        list_s(st_freq_code, nrecs+1) =                                      &
! Frequency
! DEPENDS ON: totimp
            totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 107
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                              &
            'PRELIM:TOTIMP:Error in time period conversion',                 &
            ' model=',modl_l, ' section=',isec_l, ' item=',item_l
          GO TO 9999
        END IF


!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

        IF (itima == 2) THEN
          list_s(st_freq_code, nrecs) = a_lw_radstep_diag
        ELSE IF (itima == 3) THEN
          list_s(st_freq_code, nrecs) = a_sw_radstep_diag
        ELSE IF (itima == 13) THEN
          list_s(st_freq_code, nrecs) = a_conv_step
        ELSE IF (itima == 14) THEN
          list_s(st_freq_code, nrecs) = phenol_period
        ELSE IF (itima == 15) THEN
          list_s(st_freq_code, nrecs) = triffid_period
        ELSE IF(ITIMA == 16) THEN
          list_s(st_freq_code,nrecs)=river_step_ts
        END IF

!   Period
! time mean over sampling period
        list_s(st_period_code, nrecs) =                                      &
! DEPENDS ON: totimp
          totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_period_code, nrecs)  ==  -999) THEN
          errorstatus = 108
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                              &
            'PRELIM:TOTIMP:Error in time period conversion',                 &
            ' model=',modl_l,' section=',isec_l,' item=',item_l
          GO TO 9999
        END IF

! period for timeseries recycle period
        list_s(st_period_code, nrecs+1) =                                    &
! DEPENDS ON: totimp
          totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
        IF (list_s(st_period_code, nrecs+1)  ==  -999) THEN
          errorstatus = 109
          WRITE (cmessage,'(a,a,i2,a,i3,a,i3)')                              &
            'PRELIM:TOTIMP:Error in time period conversion',                 &
            ' model=',modl_l,' section=',isec_l,' item=',item_l
          GO TO 9999
        END IF


! st_start_time for 2 record should be period for first record
! unless offset from start of run. Note value independent of logical
!  OFFSET

        IF (list_s(st_period_code, nrecs) /= -1) THEN
          IF (loffset) THEN
            list_s(st_start_time_code, nrecs+1) =                            &
              list_s(st_start_time_code, nrecs+1) -                          &
              list_s(st_period_code, nrecs+1)     +                          &
              list_s(st_freq_code, nrecs+1) - 1
          ELSE
            list_s(st_start_time_code, nrecs+1) =                            &
              list_s(st_start_time_code, nrecs+1) -                          &
              list_s(st_period_code, nrecs+1)     +                          &
              list_s(st_freq_code, nrecs+1)
          END IF
        ELSE
          list_s(st_start_time_code, nrecs+1) = 1
        END IF


!   Add both record
        list_s(st_macrotag, nrecs)=0
        nrecs = nrecs + 2

      END IF       ! Other than single time field

    END IF         ! Diag request not null - ITIM_L /= 0
999 CONTINUE
  END DO         ! Loop over diagnostic requests

END IF         ! NDIAG >  0

! DEPENDS ON: pslcom
CALL pslcom(nrecs)    ! Compress out unused pseudo levels lists

9999 CONTINUE

IF (lhook) CALL dr_hook('PRELIM',zhook_out,zhook_handle)
RETURN

END SUBROUTINE prelim
