! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INITIAL (ENDGAME VERSION) -------------------------------
!LL
!LL  Purpose: Initialises the model ready for integration/assimilation.
!LL This involves reading the model control files and setting up STASH,
!LL reading the initial or restart dump,
!LL initialising the ancillary, boundary and interface
!LL field control routines and updating the ancillary fields on restart
!LL if time to do so, exchanging coupling fields and swapping dumps (if
!LL a coupled model), and initialising the assimilation package if
!LL necessary.  Subsidiary control routines are called to perform these
!LL functions.
!LL
!LL  Programming standard: UM Doc Paper 3, version 8 (01/06/2007)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE INITIAL_4A(                                            &
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
          ,                                                             &
! ARGSP super arrays not dependent on sub-models
     &  SPD1,SPSTS,SPBND,                                               &
! ARGSP end
! ARGSPA super arrays  (atmosphere)
        A_SPDUM,A_SPPTR,A_SPCON,A_SPINF,A_SPANC,A_SPBND,A_SPSTS,        &
! ARGSPA end
! ARGSPC super array    (atmosphere-ocean coupled)
        AO_SPCPL,                                                       &
! ARGSPC end
           ngrgas,grgas_addr,                                           &
           internal_model,submodel,ngroup,meanlev)

use atm_fields_mod
Use atm_fields_bounds_Mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


USE switches, ONLY : l_aggregate,   & ! JULES tile-aggregation switch
                     l_snow_albedo

USE filenamelength_mod, ONLY:                                     &
          filenamelength


! IAU scheme:
USE IAU_mod, ONLY :    &
    L_IAU,             &
    IAU_FirstCallTS,   &
    IAU_LastCallTS,    &
    L_IAU_DumpTS0State

USE setup_iau_mod, ONLY: setup_iau

! JULES variables to be allocated (if no land points)
USE jules_mod, ONLY :  clapp_levs                                     &
                     , sathh_levs                                     &
                     ,  hcap_levs                                     &
                     ,  hcon_levs                                     &
                     ,satcon_levs                                     &
                     ,smvccl_levs                                     &
                     ,smvcwt_levs                                     &
                     ,smvcst_levs

! Declarations:
USE Submodel_Mod
USE IO
USE um_types
USE eg_alpha_mod
USE eg_alpha_ramp_mod
USE ref_pro_mod        ! needed due to the inclusion of cruntimc
USE eg_parameters_mod  ! needed due to the inclusion of cruntimc
USE Control_Max_Sizes
USE UM_ParVars
USE decomp_DB
USE ereport_mod,      ONLY : ereport     
USE integrity_mod,    ONLY : integrity_test, integrity_test_ghash
USE PrintStatus_mod
! Make sure um_sleep is picked up from portio.
USE io_dependencies

USE dynamics_input_mod, ONLY:    IntRand_seed
USE dynamics_testing_mod, ONLY:  L_Physics,L_perturb_IC_theta
USE um_input_control_mod,  ONLY:                                      &
     model_domain,        l_oasis,         oasis_couple_freq,         &
     l_veg_fracs
USE mphys_inputs_mod, ONLY:                                           &
     l_mcr_qcf2,          l_mcr_qrain,     l_mcr_qgraup,              &
     l_mcr_qcf2_lbc,      l_mcr_qrain_lbc, l_mcr_qgraup_lbc
USE river_inputs_mod, ONLY: river_step, l_rivers
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE lbc_read_data_mod, ONLY: albc_num, albc2_starttime_steps, albc_swapstep
USE ukca_option_mod, ONLY: l_ukca
USE nlstgen_mod,  ONLY: secs_per_periodim, steps_per_periodim

USE nlstcall_mod, ONLY : ncpu, &
                         Num_ALBCs, &
                         ALBC2_StartTime_mins, &
                         linterface, &
                         ltimer, &
                         model_status, &
                         model_assim_mode, &
                         run_assim_mode

USE chsunits_mod, ONLY : nunits
IMPLICIT NONE

! Common blocks:

!CL Super array lengths
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
!CL Super arrays
! TYPSPD1 super array of D1 (NOTE: only single component array)
      REAL :: SPD1(SPD1_LEN)
! TYPSPD1 end
! TYPSPDUA super array of dump arrays (atmosphere)
      REAL :: A_SPDUM(A_SPDUM_LEN)

      ! super array of auxilliary stash arrays (atmos)
      REAL :: A_SPSTS(A_SPSTS_LEN)
! TYPSPDUA end
! TYPSPST super array of STASH arrays
      REAL :: SPSTS(SPSTS_LEN)
! TYPSPST end
! TYPSPPTA super array of pointers (atmosphere)
      REAL :: A_SPPTR(A_SPPTR_LEN)
! TYPSPPTA end
! TYPSPCOA super array of constants arrays (atmosphere)
      REAL :: A_SPCON(A_SPCON_LEN)
! TYPSPCOA end
! TYPSPINA super array of output interface arrays (atmosphere)
      REAL :: A_SPINF(A_SPINF_LEN)
! TYPSPINA end
! TYPSPANA super array of ancillary file arrays (atmosphere)
      REAL :: A_SPANC(A_SPANC_LEN)
! TYPSPANA end
! TYPSPBO  super array of input boundary arrays (not sub-model specific)
      REAL :: SPBND(SPBND_LEN)
! TYPSPBO end
! TYPSPBOA super array of input boundary arrays (atmosphere)
      REAL A_SPBND(A_SPBND_LEN)
! TYPSPBOA end
! TYPSPCPL super array of coupling arrays (atmosphere-ocean)
      REAL :: AO_SPCPL(AO_SPCPL_LEN)
! TYPSPCPL end
!L       Model sizes
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
!L       Addresses of arrays within super arrays
!---------------------Start of SPINDEX--------------------------------
!L
!L --------------- D1    array  -----------------------------------
!L (now includes extra copies of atmos/ocean D1 for MPP)
!L
      INTEGER     IXD1_LEN       ! No. of arrays
      PARAMETER(  IXD1_LEN =     4)
      INTEGER     IXD1           ! Addresses of arrays
      COMMON/    CIXD1/      IXD1(IXD1_LEN)
!L
!L --------------- Dump headers--------------------------
!L
!L ATMOS
      INTEGER   A_IXDUM_LEN       ! No. of arrays
      PARAMETER(A_IXDUM_LEN = 14)
      INTEGER   A_IXDUM           ! Addresses of arrays
      COMMON/  CA_IXDUM/ A_IXDUM(A_IXDUM_LEN)
      INTEGER   A_IXSTS_LEN
      PARAMETER(A_IXSTS_LEN = 12)
      INTEGER   A_IXSTS
      COMMON/  CA_IXSTS/ A_IXSTS(A_IXSTS_LEN)
!L
!L --------------- STASH arrays -----------------------------------
!L
      INTEGER     IXSTS_LEN       ! No. of arrays
      PARAMETER(  IXSTS_LEN = 14)
      INTEGER     IXSTS           ! Addresses of arrays
      COMMON/    CIXSTS/      IXSTS(IXSTS_LEN)
!L
!L --------------- Pointers in D1 array and row/level dependent ---
!L --------------- constants
!L ATMOS
      INTEGER   A_IXPTR_LEN       ! No. of arrays
      PARAMETER(A_IXPTR_LEN = 239)
      INTEGER   A_IXPTR           ! Addresses of arrays
      COMMON/  CA_IXPTR/ A_IXPTR(A_IXPTR_LEN)
!L
!L --------------- Pre-calculated arrays of constants -------------
!L
!L ATMOS
      INTEGER   A_IXCON_LEN       ! No. of arrays
! A_IXCON_LEN is so low b/c the bulk of constants that previously were
! stored in um_index.a are now allocated directly in SETCONA 
      PARAMETER(A_IXCON_LEN = 3)
      INTEGER   A_IXCON           ! Addresses of arrays
      COMMON/  CA_IXCON/ A_IXCON(A_IXCON_LEN)
!L
!L --------------- Headers for output interface datasets (boundary
!L                 conditions out)
!L ATMOS
      INTEGER   A_IXINF_LEN       ! No. of arrays
      PARAMETER(A_IXINF_LEN = 9 )
      INTEGER   A_IXINF           ! Addresses of arrays
      COMMON/  CA_IXINF/ A_IXINF(A_IXINF_LEN)
!L
!L --------------- Headers for ancillary files -------------------
!L
!L ATMOS
      INTEGER   A_IXANC_LEN       ! No. of arrays
      PARAMETER(A_IXANC_LEN = 4)
      INTEGER   A_IXANC           ! Addresses of arrays
      COMMON/  CA_IXANC/ A_IXANC(A_IXANC_LEN)
!L
!L --------------- Headers from input boundary files -------------
!L
!L NOT SUB-MODEL DEPENDENT
      INTEGER     IXBND_LEN       ! No. of arrays
      PARAMETER(  IXBND_LEN = 1)
      INTEGER     IXBND           ! Addresses of arrays
      COMMON/    CIXBND/ IXBND(IXBND_LEN)
!L
!L ATMOS
      INTEGER   A_IXBND_LEN       ! No. of arrays
      PARAMETER(A_IXBND_LEN = 5)
      INTEGER   A_IXBND           ! Addresses of arrays
      COMMON/  CA_IXBND/ A_IXBND(A_IXBND_LEN)
!L
!L --------------- Constant arrays needed for atmosphere-ocean----
!L --------------- coupling
!L
      INTEGER   AO_IXCPL_LEN      ! No. of arrays
      PARAMETER(AO_IXCPL_LEN = 10)
      INTEGER   AO_IXCPL          ! Addresses of arrays
      COMMON/ CAO_IXCPL/ AO_IXCPL(AO_IXCPL_LEN)
!L
!-------------------End of SPINDEX------------------------------------
! --------------------- Comdeck: CHISTORY ----------------------------
!
!  Purpose: COMMON block for history data needed by top level (C0)
!           routines, and passed from run to run.  Mostly set by
!           the User Interface.
!
!           Note that CHISTORY *CALLs ALL individual history comdecks
!
! --------------------------------------------------------------------
!
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations


      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: model_data_time(6)

      ! Indicator of operational run type
      INTEGER :: run_indic_op

      ! Final target date for the run
      INTEGER :: run_resubmit_target(6)

      ! Last field written/read per FT unit
      INTEGER :: ft_lastfield(20:nunits)

      ! Number of automatically-resubmitted job chunks
      ! Used to name output file
      INTEGER :: run_job_counter

! History Common Block for overall model integers variables.

      COMMON /ihisto/                                                 &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

      NAMELIST /nlihisto/                                             &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

  CHARACTER(LEN=8) ::  run_type             ! Type of run
  CHARACTER(LEN=1) ::  ft_active(20:nunits) ! "Y" if file partly written

  LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

  ! History Common Block for overall model character variables.

  COMMON /chisto/                                     &
     run_type,                                        &
     newrun, ft_active

  NAMELIST /nlchisto/                                 &
     run_type,                                        &
     ft_active

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! This file belongs in section: Top Level
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations
      ! History block copy of A_STEP held in file CTIME
      INTEGER :: h_stepim(n_internal_model_max)

      ! No of means activated
      INTEGER :: mean_offsetim(n_internal_model_max)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: offset_dumpsim(n_internal_model_max)

      ! No of mean periods chosen
      INTEGER :: mean_numberim(n_internal_model_max)

      ! Indicators used to correct logical units are used for
      ! atmos partial sum dump I/O
      INTEGER :: run_meanctl_indicim(4,n_internal_model_max)

      ! History Common Block for generic model integer variables.

      COMMON /ihistg/                                         &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         mean_numberim, run_meanctl_indicim

      NAMELIST /nlihistg/                                     &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         run_meanctl_indicim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character variables for
!              managing dump names
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 
!
!   Type declarations
!
! For keeping old restart dump name between creation of new restart dump
! and successful completion of climate means and history file.
CHARACTER(LEN=256) :: save_dumpname_im(n_internal_model_max)
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
! Blank name
CHARACTER(LEN=256) :: blank_file_name
!
! History Common Block for generic model characters variables.
!
COMMON /chistg/save_dumpname_im, checkpoint_dump_im, blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

! CHISTG end
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200

!
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

! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! This file belongs in section: Top Level

!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:

!------------------   End of Physics   ---------------------------------

      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(model_levels_max) ! Timestep of max w
      INTEGER :: time_div_max(model_levels_max) ! Timestep of max div
      INTEGER :: time_div_min(model_levels_max) ! Timestep of min div
      INTEGER :: time_lapse_min(model_levels_max) ! Timestep of min
      INTEGER :: time_max_shear(model_levels_max) !Timestep max shear
      INTEGER :: time_max_wind(model_levels_max) ! Timestep of max wind
      INTEGER :: time_KE_max(model_levels_max) ! Timestep of max KE
      INTEGER :: time_KE_min(model_levels_max) ! Timestep of min KE
      INTEGER :: time_noise_max(model_levels_max) ! Timestep of max

      REAL:: frictional_timescale(model_levels_max) ! For idealised case
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(0:model_levels_max) ! Max w at a level
      REAL :: max_div_run(model_levels_max) ! Max divergence at a level
      REAL :: min_div_run(model_levels_max) ! Min divergence at a level
      REAL :: min_lapse_run(model_levels_max) ! Min dtheta/dz at a level
      REAL :: max_shear_run(model_levels_max) ! Max shear at a level
      REAL :: max_wind_run(model_levels_max) ! Max wind at a level
      REAL :: max_KE_run(model_levels_max)   ! Max KE at a level
      REAL :: min_KE_run(model_levels_max)   ! Min KE at a level
      REAL :: max_noise_run(model_levels_max) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA
!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
        rpemax, rpemin, ipesum, rpesum,                                 &
        max_w_run, min_theta1_run, dtheta1_run,                         &
        max_div_run, min_div_run, min_lapse_run,                        &
        max_shear_run, max_wind_run,                                    &
        max_noise_run, max_KE_run, min_KE_run,                          &
        time_KE_max, time_KE_min,                                       &
        time_w_max, time_div_max, time_div_min, time_lapse_min,         &
        time_max_shear, time_max_wind,                                  &
        time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(model_levels_max)
      REAL :: friction_level(model_levels_max)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                              &
       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
       SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
       base_frictional_timescale, SuHe_sigma_cutoff,                   &
       L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
       SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: vert_grid_ratio
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(model_levels_max)
      REAL :: rho_ref(model_levels_max)
      REAL :: exner_ref(model_levels_max + 1)
      REAL :: q_ref(model_levels_max)
      REAL :: u_ref(model_levels_max)
      REAL :: v_ref(model_levels_max)
      REAL :: z_orog_print(0:model_levels_max)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: qforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: uforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: vforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind

! ENDGAME
      REAL :: T_surface
      REAL :: Eccentricity
      ! Following two variables used only if L_rotate_grid=.true.
      REAL :: grid_NP_lon ! Longitude (degrees) of grid's north pole
      REAL :: grid_NP_lat ! Latitude (degrees) of grid's north pole
      REAL :: AA_jet_u0   ! See QJRMS 133,1605--1614
      REAL :: AA_jet_A    !
      REAL :: theta_pert
      REAL :: ring_height
      REAL :: angular_velocity ! Planet's angular velocity
      REAL :: T0_P, T0_E ! deep atmosphere baroclinic wave surface temperatures
      INTEGER :: Trefer_number
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: AA_jet_m  ! See QJRMS 133,1605--1614
      INTEGER :: AA_jet_n  !
      INTEGER :: chain_number ! Run continuation number
      LOGICAL :: L_rotate_grid    ! .true. for rotating North pole of grid
      LOGICAL :: L_baro_Perturbed ! Used for baroclinic test to specify
                                  ! pert or steady
      LOGICAL :: L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,  &
                 L_HeldSuarez2_drag,                                        &
                 L_baro_inst, L_isothermal, L_exact_profile, L_balanced,    &
                 L_solid_body
      LOGICAL :: L_deep_baro_inst ! deep atmosphere baroclinic wave switch          


      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type
      INTEGER :: b_const, k_const ! deep atmosphere baroclinic wave parameters

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                              &
       h_o, h_o_actual, h_o_per_step,                                  &
       lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
       Witch_power, plat_size_x, plat_size_y,                          &
       height_domain, delta_x, delta_y, big_factor, mag, vert_grid_ratio, &
       first_theta_height, thin_theta_height, p_surface,               &
       theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
       u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
       ujet_lat, ujet_width,                                           &
       t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
       u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
       z_orog_print, grow_steps,                                       &
       surface_type, grid_number, grid_flat,                           &
       tprofile_number, qprofile_number, uvprofile_number,             &
       num_profile_data, num_uvprofile_data,                           &
       tforce_option, qforce_option, uvforce_option,                   &
       num_tforce_levels, num_tforce_times,                            &
       num_qforce_levels, num_qforce_times,                            &
       num_uvforce_levels, num_uvforce_times,                          &
       L_pforce, pforce_option, num_pforce_times,                      &
       first_constant_r_rho_level_new,                                 &
       big_layers, transit_layers, mod_layers,                         &
       zprofile_data, tprofile_data, qprofile_data,                    &
       z_uvprofile_data, uprofile_data, vprofile_data,                 &
       tforce_time_interval, qforce_time_interval,                     &
       uvforce_time_interval, newtonian_timescale,                     &
       z_tforce_data, z_qforce_data, z_uvforce_data,                   &
       tforce_data, qforce_data, uforce_data, vforce_data,             &
       tforce_data_modlev, qforce_data_modlev,                         &
       uforce_data_modlev, vforce_data_modlev,                         &
       pforce_time_interval, p_surface_data,                           &
       L_initialise_data,                                              &
       L_perturb_t, perturb_magnitude_t,                               &
       L_perturb_q, perturb_magnitude_q,                               &
       L_perturb_correlate_tq,                                         &
       L_perturb_correlate_vert,                                       &
       L_perturb_correlate_time,                                       &
       perturb_type, perturb_height,                                   &
       L_constant_dz, L_polar_wind_zero,                               &
       L_wind_balance, L_rotate_winds,                                 &
       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
       L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
       L_pressure_balance, L_vert_Coriolis,                            &
       cool_rate, L_force, L_force_lbc,                                &
       zprofile_orog, idl_interp_option, hf,                           &
       L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
       L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
       idl_bubble_option, idl_bubble_max                               &
      , idl_bubble_height, idl_bubble_width, idl_bubble_depth          &
      , idl_bubble_xoffset,idl_bubble_yoffset                          &
      , L_idl_bubble_saturate,                                         &
       L_rotating, L_fixed_lbcs, L_code_test,                          &
       L_baroclinic, L_cyclone,                                        &
       L_damp, L_geo_for, L_bomex,                                     &
       DMPTIM, HDMP, ZDMP,                                             &
       u_geo, v_geo,                                                   &
!ENDGAME
       T_surface, chain_number,                                        &
       Trefer_number,                                                  &
       tstep_plot_frequency, tstep_plot_start, Eccentricity,           &
       L_rotate_grid, grid_NP_lon, grid_NP_lat,                        &
       AA_jet_m, AA_jet_n, AA_jet_u0, AA_jet_A, L_baro_Perturbed,      &
       L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,       &
       L_HeldSuarez2_drag,                                             &
       L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const, k_const,    &      
       ring_height, theta_pert, L_isothermal,                          &
       L_exact_profile, L_balanced, L_solid_body, angular_velocity
! CRUNTIMC end
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end
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
!
! Description:
!   This comdeck declares an integer variable 'ModelType' whose value
!   determines whether a model run is global, limited area or zonal.
!   The values of ModelType associated with each of the run types are
!   defined by integer parameters which are also declared below.
!    ModelType is set in subroutine SETLOGIC.
!
      ! Value used to represent the global model configuration
      INTEGER,PARAMETER:: GlobalModel=1

      ! Value used to represent the limited area model configuration
      INTEGER,PARAMETER:: LimitedAreaModel=2

      ! Value used to represent the 'periodic in x' model config
      ! INTEGER,PARAMETER:: ZonalModel=2

! Global scalars:
      INTEGER     ModelType  ! Integer switch which is equated to one
                             ! of the above parameters in a model run,
                             ! and so determines the configuration

! COMMON blocks:
      COMMON /RunType/ ModelType

! C_GLOBAL end
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
! --------------------- COMDECK LBC_COUP --------------------------
!    Description:
!       This COMDECK stores the variables connected with the
!       Parallel Running between the Global and Mesoscale models.
!       The Mesoscale has to be held back until there are sufficient
!       Boundary Conditions (BCs) to proceed.
!
      logical l_lbc_coup        ! T : Global/Meso Coupling on
      integer um_lbc_stream     ! Output Stream generating BCs.
      
      ! UM 6.5 -  MODEL_ANALYSIS_HRS changed to REAL - 
      !                requires LBC_FC_HRS to REAL also
      real    lbc_fc_hrs        ! Forecast time w.r.t analysis time
      integer um_lbc_wait       ! Wait time between re-tries if BCs
                                ! not available.
      integer um_lbc_wait_max   ! Maximum wait time.
      CHARACTER(LEN=256) lbc_filename ! Name of file that communicates between
                                 ! Global and Meso.

      COMMON /LBC_COUP/ l_lbc_coup, um_lbc_stream, lbc_fc_hrs,          &
     &                  um_lbc_wait, um_lbc_wait_max, lbc_filename
! ----------------------------------------------------------------------

! Subroutine arguments:

! ENDGAME-only
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
! End of ENDGAME-only

      INTEGER internal_model ! OUT internal model identifier:
!                            !   1:Atmos; 2:Ocean; 3:Slab ; etc
      INTEGER submodel       ! OUT submodel partition (dump) identifier:
!                            !   1:Atmos; 2:Ocean; etc
      INTEGER NGROUP         ! OUT   - No of steps in "group"n
      INTEGER MEANLEV        ! OUT - Mean level indicator

! 3-D fields of species to be passed down to radiation
      INTEGER, INTENT(IN) :: ngrgas
      INTEGER, INTENT(OUT) :: grgas_addr(ngrgas)

! Local variables:

      INTEGER  IMEAN      ! Loop index over mean periods
      INTEGER  I          ! Loop index
      INTEGER  ISM        ! Loop index over submodels
      INTEGER LL        ! Counter
      INTEGER LEN_FILENAME ! Length of FILENAME array
      INTEGER Dummyarg  ! Not used, needed to end arg list
      INTEGER SECS_PER_STEPA ! Atmos timestep length in seconds
      
      CHARACTER(LEN=filenamelength) :: filename
      CHARACTER(LEN=2) ENVALUE
!
      integer len_wait_tot     ! Total wait time for boundary data
      CHARACTER(LEN=8) ch_date2     ! Date from date_and_time
      CHARACTER(LEN=10) ch_time2    ! Time from date_and_time
      integer info             ! Return Code from GCOM routine.
      integer lbc_ntimes       ! No of BC's in communication file
      integer ms_ntimes        ! No of BC's required in mesoscale
      INTEGER um_lbc_wait_usec ! No of microseconds to wait

      LOGICAL                                                           &
         not_enough_lbcs,                                                 &
                          ! True if more LBCs are required
         read_more_values  ! True if we are to read more values in


! Number of dust bins in LBC generation, disabled if not using MakeBC
      INTEGER :: ndustbin_in, ndustbin_out 


! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER(LEN=256) Cmessage    ! Error message
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='INITIAL_4A')
! Mixing ratio code
      Logical, Parameter ::                                             &
         l_mr_iau = .FALSE.  ! Use mixing ratio code (if available)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header ------------------------------------------------------

      IF (lhook) CALL dr_hook('INITIAL_4A',zhook_in,zhook_handle)

!
! this #if is required because parts of the integrity test module
! can only be compiled if the configuration includes certain libraries
! and excludes some modules and include files from the dependency 
! search. Making it a runtime logical would simply risk it being set without
! the prerequisit override being applied.

      ICODE=0
      Cmessage=''

!L
!L----------------------------------------------------------------------
!L
!L 1.2 Set FT units as inactive on first step of the integration
!L     and set last field written/read to zero for each unit
!L
      IF (H_STEPim(a_im) == 0) THEN
        DO I=20,NUNITS
          FT_ACTIVE(I)='N'
          FT_LASTFIELD(I)=0
        ENDDO
      ENDIF
!L
!L 1.3 Option to write RADINCS array to file cache2 on unit 16 removed
!L
      IF (PrintStatus  >=  PrStatus_Normal) THEN
      WRITE(6,*) 'Fast i/o of radincs directly to core memory'
      ENDIF
!L
!L Number of CPUs attached to program (ncpu) is hard-wired to 1.
!L
      ncpu = 1
!L
!L---------------------------------------------------------------------
!L 2. Initialise STASH control arrays from STASH control file.
!L---------------------------------------------------------------------


! Note that NSECTS=NSECTP, NITEMS=NITEMP : set in WSTLST

! DEPENDS ON: initctl
      CALL INITCTL(                                                     &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
         ICODE,CMESSAGE )

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITCTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF


!L
!L----------------------------------------------------------------------
!L 3. Read appropriate submodel partition dump to memory.  If coupled,
!L    page out the D1 part of each dump to its 'swap' file and read the
!L    other dump(s) into memory.  Write temporary history file if dumps
!L    read successfully on timestep 0.
!L
!L
!L 3.1  Loop over submodel partition dumps
!L
      DO ISM=1,N_SUBMODEL_PARTITION

        submodel=SUBMODEL_PARTITION_LIST(ISM)

! DEPENDS ON: initdump
        CALL INITDUMP(                                                  &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
     &               submodel,ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITDUMP'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDDO  ! ISM over submodel partitions



! SET_ATM_FIELDS points the fields in atm_fields_mod to the correct parts of D1
! It should be called after INITDUMP (which calls SET_ATM_POINTERS)
! and before INITDIAG (which calls ST_DIAG<n>).

! DEPENDS ON: set_atm_fields
      CALL Set_Atm_Fields(                                              &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
         SPD1(IXD1(2)), SPD1(IXD1(2)), SPD1(IXD1(2)))

      IF (mype  ==  0) THEN
! DEPENDS ON: eg_check_atm_fields
        CALL eg_check_Atm_Fields(                                       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
           SPD1(IXD1(2)), SPD1(IXD1(2)), SPD1(IXD1(2)))
      END IF
!L
!L      Set RUN indicator in atmosphere dump header
!L
! DEPENDS ON: set_run_indic_op
      CALL SET_RUN_INDIC_OP(                                            &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
         ICODE,CMESSAGE)
      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'Failure in call to SET_RUN_INDIC_OP'
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF
!L
!L 3.3  On NRUN initialise dump LOOKUP headers associated with
!L      diagnostic fields with the bare essentials needed to read and
!L      write dumps - the rest to be filled in by STASH during the run.
!L
      IF (H_STEPim(a_im)  == 0) THEN
! DEPENDS ON: inithdrs
        CALL INITHDRS(                                                  &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
     &                ICODE,CMESSAGE)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITHDRS'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDIF    ! End test for NRUN
!
      if(L_perturb_IC_theta .AND. H_STEPim(a_im) == 0) then
! DEPENDS ON: perturb_theta_ctl
       call perturb_theta_ctl(                                          &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

          row_length, rows, model_levels,          &
          global_row_length, global_rows,          &
          model_domain, at_extremity,              &
          offx, offy, IntRand_Seed,   datastart    &
          )
     endif

!L
!L 3.4 Write out temporary copy of D1 for current submodel
!L

      icode=0
      IF (L_WRIT_INIT) THEN

        DO ISM=1,N_SUBMODEL_PARTITION

          submodel=SUBMODEL_PARTITION_LIST(ISM)
        if (submodel  ==  atmos_sm) then
! DEPENDS ON: dumpctl
           CALL DUMPCTL (                                               &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
     &          atmos_sm,0,.TRUE.,'atminitial',0,                       &
     &          ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to DUMPCTL (Atmos)'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF

        endif
        ENDDO ! ISM

      END IF

!L----------------------------------------------------------------------
!L 6.  Initialise means program control block
!L
      DO ISM=1,N_SUBMODEL_PARTITION

        submodel=SUBMODEL_PARTITION_LIST(ISM)

! DEPENDS ON: initmean
        CALL INITMEAN(                                                  &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
           submodel,ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITMEAN for submodel ',       &
             ISM
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDDO ! ISM over submodel partition dumps
!L----------------------------------------------------------------------
!L 4. Set up other control blocks and physical constants
!L
!L 4.1  Initialise the model time and check that history file data time
!L      matches dump(s); set derived time/step information
!L

! DEPENDS ON: inittime
      CALL INITTIME(                                                    &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
         submodel,ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITTIME'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 4.2  Write up temporary history file after successfully reading
!L      initial dumps on timestep 0 and setting model_data_time if
!L      assimilation run, to allow CRUN from initial dumps.
!L
      IF (mype  ==  0) THEN
        IF (H_STEPim(a_im) == 0) THEN
! DEPENDS ON: temphist
          CALL TEMPHIST(XHIST_UNIT,ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to TEMPHIST'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
        END IF
      ENDIF
!L
!L 4.3  Set up control block for updating ancillary fields
!L

! DEPENDS ON: inancctl
      CALL INANCCTL(                                                    &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Ancillary file arrays(atmosphere)-------------
      A_SPANC(A_IXANC( 1)),A_SPANC(A_IXANC( 2)),A_SPANC(A_IXANC( 3)),   &
      A_SPANC(A_IXANC( 4)),                                             &
     &           ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INANCCTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 4.4  Set up control block for updating boundary fields

      IF (model_domain /= mt_global) THEN

      ALBC_num = 1

      ! Setup for use of two atmos boundary files:
      IF (Num_ALBCs == 2) THEN

        ! Check namelist variable ALBC2_StartTime_mins:
        SECS_PER_STEPA = NINT(SECS_PER_STEPim(a_im))

        IF (MOD(ALBC2_StartTime_mins*60, SECS_PER_STEPA) /= 0) THEN
          ICODE = 1
          WRITE (CMessage,*)                                            &
             'ALBC2_StartTime_mins (', ALBC2_StartTime_mins,             &
             ') does not coincide with the start of a timestep'
          CALL EReport (RoutineName, ICODE, CMessage)
        ELSE
          ! Convert into steps:
          ALBC2_StartTime_steps = ALBC2_StartTime_mins*60/SECS_PER_STEPA
        END IF

        ! If on a continuation run, we may be able to go straight to
        ! the second boundary file:
        IF (STEPim(a_im) >= ALBC2_StartTime_steps) ALBC_num = 2

      END IF

      ! NOTE: If using two boundary files, coupling is only activated
      !       for the second.
      IF (l_lbc_coup .AND.                                              &
         .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

        ms_ntimes=2  ! For now we assume the minimum requirement -
                     ! two sets of boundary conditions are required.
                     ! Once INBOUND has been called we may find that
                     ! we need more than this, in which case we will
                     ! jump back to line 100

      ENDIF
      
      END IF  ! not global

 100  CONTINUE  ! Return here if IN_BOUND has been called and there
                ! are insufficient boundary conditions to proceed.
                ! This condition may arise in CRUNS


      IF (model_domain /= mt_global) THEN

      IF (l_lbc_coup .AND.                                              &
         .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

        CALL DATE_AND_TIME(ch_date2, ch_time2)

        IF (PrintStatus  >=  PrStatus_Normal) THEN
          WRITE(6,*) 'LBC_COUP: ',                                      &
             ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',   &
             ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),          &
             ' Wait to call INBOUND in INITIAL.'
        ENDIF ! (PrintStatus  >=  PrStatus_Normal)

! DEPENDS ON: timer
        Call Timer ('LBC_WAIT',3)

        IF (mype  ==  0) THEN  ! Only do the check on the communication
                               ! file on PE 0
          WRITE(6,*) 'ms_ntimes in INITIAL ',ms_ntimes

          not_enough_lbcs=.TRUE.
          len_wait_tot=0

          DO WHILE (not_enough_lbcs) ! loop until the communication
                                     ! file indicates we have enough
                                     ! lbcs to start running with
            read_more_values=.TRUE.

            CLOSE (190) ! We need to close and reopen the communication
                        ! file to ensure we see its latest state
            OPEN (190,FILE=lbc_filename,ACTION="read",IOSTAT=ICODE)

            IF (ICODE  /=  0) THEN ! Something went wrong with the OPEN
              WRITE(6,*) 'Return code from OPEN of communication file ',&
                 'in INITIAL: ',ICODE
              ICODE=401
              WRITE(CMESSAGE,*) 'INITIAL : OPEN failed for Unit 190'

              read_more_values=.FALSE. ! Jump out of the loops
              not_enough_lbcs=.FALSE.  ! and call Ereport
            ENDIF ! IF (ICODE  /=  0)

            DO WHILE (read_more_values) ! loop while we need to read
                                        ! more values from the
                                        ! communication file

              READ (190,*,IOSTAT=ICODE) lbc_ntimes ! read the next value

              IF (ICODE  /=  0) THEN ! Problem with the READ...

                WRITE(6,*) 'ms : Return code from READ ',ICODE

                IF (len_wait_tot  >=  um_lbc_wait_max) THEN

                  ! Maximum wait time has been exceeded.
                  ! Insufficient Boundary Conditions to proceed.
                  ! Likely cause is delay in job generating the BC's.

                  WRITE(6,*) 'ms: Maximum wait time reached ',          &
                     'after ',um_lbc_wait_max,' seconds.'
                  ICODE=402
                  WRITE(CMESSAGE,*) 'INITIAL : Failed to find ',        &
                     'required value in LBC_FILE.'

                  read_more_values=.FALSE. ! Jump out of the loops
                  not_enough_lbcs=.FALSE.  ! and call Ereport

                ELSE ! We've not exceeded the time limit

                  ! Insufficient BCs to proceed.
                  ! Wait for um_lbc_wait seconds before another
                  ! attempt to read the file to see if more BCs
                  ! have been written.

                  WRITE(6,*) 'ms: Wait for ',um_lbc_wait,               &
                     ' seconds and retry.'

                  ! Calculate microseconds from seconds.
                  um_lbc_wait_usec = 1000000*um_lbc_wait
                  CALL um_sleep(um_lbc_wait_usec)

                  len_wait_tot=len_wait_tot+um_lbc_wait

                  WRITE(6,*) 'ms : Total wait so far ',len_wait_tot,    &
                     ' seconds.'

                  read_more_values=.FALSE.
                  ! No more values in the file as it stands, so this
                  ! forces us to close and reopen the file and read
                  ! again to see if any new values have arrived.

                ENDIF ! IF (len_wait_tot  >=  um_lbc_wait_max)

              ELSE ! the READ(190) was successful - now we must
                   ! interpret the value pulled from the file

                IF (lbc_ntimes  >   1000) THEN
                  ! First value in the file is always >1000

                  read_more_values=.TRUE. ! Go round and read next value

                ELSEIF (lbc_ntimes  <   ms_ntimes) THEN
                  ! Value is not required. Proceed to read next value

                  WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,           &
                     ' read in.'//                              &
                     ' lbc_ntimes >= ',ms_ntimes,               &
                     ' is required. Read next value.'
                  read_more_values=.TRUE.

                ELSEIF (lbc_ntimes  >=  ms_ntimes) THEN
                  ! Required value read in. Sufficient BCs to proceed

                  WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,           &
                     ' read in.'//                              &
                     ' lbc_ntimes >= ',ms_ntimes,               &
                     ' is required. Proceed.'

                  CALL DATE_AND_TIME (ch_date2, ch_time2)

                  IF (PrintStatus  >=  PrStatus_Normal) THEN
                    WRITE(6,*) 'LBC_COUP: ',                            &
                       ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),  &
                       ' on ',                                             &
                       ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),  &
                       ' Proceed to call INBOUND in INITIAL.'
                  ENDIF ! (PrintStatus  >=  PrStatus_Normal)

                  read_more_values=.FALSE. ! Don't read any more values
                  not_enough_lbcs=.FALSE.  ! Don't re-examine the file

                ENDIF ! Test for value of lbc_ntimes

              ENDIF ! IF (ICODE  /=  0)

            ENDDO ! DO WHILE (read_more_values)

          ENDDO ! DO WHILE (not_enough_lbcs)

        ENDIF ! IF (mype  ==  0)

! DEPENDS ON: timer
        Call Timer ('LBC_WAIT',4)

        CALL GC_IBCAST(100,1,0,nproc,info,ICODE) ! PE 0 broadcasts ICODE
                                                 ! to all processors
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'INITIAL - Error detected in LBC coupling.'
          WRITE(6,*) 'ICODE: ',ICODE
          WRITE(6,*) 'CMESSAGE : ',CMESSAGE

          CALL Ereport(RoutineName,ICODE,CMESSAGE)
        ENDIF ! IF (ICODE  /=  0)

        CALL GC_IBCAST(100,1,0,nproc,info,lbc_ntimes) ! PE 0 broadcasts
                                                      ! lbc_ntimes to
                                                      ! all processors
      ENDIF ! (l_lbc_coup .AND.
            !  .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1))

! DEPENDS ON: in_bound
      CALL IN_BOUND(                                                    &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Input boundary arrays ------------------------
        SPBND(  IXBND( 1)),                                             &
!L --------------- Input boundary arrays(atmosphere)-------------
     &A_SPBND(A_IXBND( 1)),A_SPBND(A_IXBND( 2)),A_SPBND(A_IXBND( 3)),   &
     &A_SPBND(A_IXBND( 4)),A_SPBND(A_IXBND( 5)),                        &
         A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                                 &
         A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,                                 &
         A_LEN1_COLDEPC,A_LEN2_COLDEPC,                                 &
                                          ! for dynamic array
     &                   ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to IN_BOUND'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      IF (l_lbc_coup .AND.                                              &
         .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

!       Now that IN_BOUND has been called for the first time
!       double check that there are sufficient BCs to proceed.
!
!       Determine which boundary data is required to proceed
!       the next period.
        if (boundary_stepsim(a_im) >  0) then
          IF (ALBC_num == 2) THEN
            ms_ntimes = 1 + ( (stepim(a_im)-ALBC_SwapStep)              &
               / boundary_stepsim(a_im) )
          ELSE
            ms_ntimes = 2 + (stepim(a_im)/boundary_stepsim(a_im))
          END IF
        endif

        if (lbc_ntimes <  ms_ntimes) then

!         There are insufficient BCs to proceed. Go back and wait
!         for sufficient BCs to proceed.
          WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,                   &
             ' lbc_ntimes >= ',ms_ntimes,' is required. '//     &
             'Insufficient BCs to proceed. Wait.'

          GOTO 100

        endif

      endif  ! (l_lbc_coup .AND.
             ! .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1))

!
!    4.4.1  Update atmosphere boundary fields at step zero
!

      IF (BOUNDARY_STEPSim(a_im)  /=  0) THEN ! If there are BCs

        IF (l_lbc_coup .AND.                                            &
           .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                            ! coupling

          IF (PrintStatus  >=  PrStatus_Normal) THEN
            WRITE(6,*) 'LBC_COUP: ',                                    &
               ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),        &
               ' on ',                                                   &
               ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),        &
               ' Proceed to call UPBOUND in INITIAL.'
          ENDIF ! (PrintStatus  >=  PrStatus_Normal)

        ENDIF

! DEPENDS ON: up_bound
        CALL UP_BOUND(submodel,                                         &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Input boundary arrays ------------------------
        SPBND(  IXBND( 1)),                                             &
!L --------------- Input boundary arrays(atmosphere)-------------
     &A_SPBND(A_IXBND( 1)),A_SPBND(A_IXBND( 2)),A_SPBND(A_IXBND( 3)),   &
     &A_SPBND(A_IXBND( 4)),A_SPBND(A_IXBND( 5)),                        &
     &              ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'INITIAL: Failure in call to atmosphere UP_BOUND'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

        IF (Num_ALBCs == 2) THEN

          ! If the data for the start and end of the current boundary
          ! data interval is to come from different boundary files, we
          ! need to make additional calls to INBOUNDA/UP_BOUND:
          IF (ALBC_num == 1 .AND. STEPim(a_im) >= ALBC_SwapStep) THEN

            IF (PrintStatus >= PrStatus_Normal) THEN
              WRITE (6,*) ''
              WRITE (6,*) 'INITIAL: Swapping to 2nd atmos boundary file'
              WRITE (6,*) ''
            END IF

            ALBC_num = 2

            IF (l_lbc_coup) THEN

              ms_ntimes = 1
              GOTO 100

            ELSE

! DEPENDS ON: inbounda
              CALL INBOUNDA(                                            &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Input boundary arrays ------------------------
        SPBND(  IXBND( 1)),                                             &
!L --------------- Input boundary arrays(atmosphere)-------------
     &A_SPBND(A_IXBND( 1)),A_SPBND(A_IXBND( 2)),A_SPBND(A_IXBND( 3)),   &
     &A_SPBND(A_IXBND( 4)),A_SPBND(A_IXBND( 5)),                        &
     &          A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                          &
     &          A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,                          &
     &          A_LEN1_COLDEPC,A_LEN2_COLDEPC)

! DEPENDS ON: up_bound
              CALL UP_BOUND(submodel,                                   &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Input boundary arrays ------------------------
        SPBND(  IXBND( 1)),                                             &
!L --------------- Input boundary arrays(atmosphere)-------------
     &A_SPBND(A_IXBND( 1)),A_SPBND(A_IXBND( 2)),A_SPBND(A_IXBND( 3)),   &
     &A_SPBND(A_IXBND( 4)),A_SPBND(A_IXBND( 5)),                        &
     &          ICODE,CMESSAGE)


              IF (ICODE /= 0) THEN
                WRITE(6,*)                                              &
                   'INITIAL: Failure in call to atmosphere UP_BOUND'
                CALL EReport (RoutineName, ICODE, CMESSAGE)
              END IF

            END IF ! (l_lbc_coup)

          END IF ! (ALBC_num == 1 .AND. STEPim(a_im) >= ALBC_SwapStep)

        END IF ! (Num_ALBCs == 2)

      ENDIF ! IF (BOUNDARY_STEPSim(a_im)  /=  0)

      END IF  ! .NOT. GLOBAL

!
!    4.4.2  Call SETCONA_4A
!
! VATPOLES settings
!
! Set loop index bounds for arrays defined on p, u, v points
!
!      Select Case(model_domain)
!         Case(mt_global)
!            i_p_start = 1
!            i_p_end   = row_length 
!            j_p_start = 1 
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0 
!            j_v_end   = n_rows - 1
!         Case(mt_cyclic_lam)
!            i_p_start = 1
!            i_p_end   = row_length
!            j_p_start = 1
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0
!            j_v_end   = n_rows - 1
!         Case(mt_bi_cyclic_lam)
!            i_p_start = 1
!            i_p_end   = row_length 
!            j_p_start = 1 
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0 
!            j_v_end   = rows - 1
!         Case(mt_lam)
!            i_p_start = 1
!            i_p_end   = row_length - 1 
!            j_p_start = 1 
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length - 1
!            j_v_start = 0 
!            j_v_end   = n_rows - 1
!         Case Default
!            Print*,'Invalid domain type :',model_domain
!            Stop
!         End Select
!
!      i_u_start = 0
!      i_u_end   = row_length-1 
!      j_u_start = 1 
!      j_u_end   = rows
!      k_p_start = 1
!      k_p_end   = model_levels
!      k_w_start = 0
!      k_w_end   = model_levels
!
! End of VATPOLES settings

! DEPENDS ON: setcona_ctl_4A
      CALL SETCONA_CTL_4A(                                                 &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
        icode,cmessage)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'INITIAL: Failure in call to SETCONA_4A'
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF

!L
!L  4.5  Set up control block for writing interface fields.
!L

! DEPENDS ON: intf_ctl
      CALL  INTF_CTL (                                                  &
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),                                     &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

         ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INTF_CTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 4.6  Initialise physical constants used in main physics
!L      packages - includes radiation lookup tables
!L

      IF (L_Physics) THEN

! DEPENDS ON: initphys
      CALL INITPHYS(ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITPHYS'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDIF ! on L_Physics

      IF (L_emcorr) THEN ! Energy correction enabled
!L--------------------------------------------------------------------
!L 4.7  Initialise total atmospheric energy & energy correction
!L--------------------------------------------------------------------
!L  Only done for a new run and then only if the header values in the
!L  dump are set to missing data. If the arrays were set at the
!L  beginning of all NRUNS (new runs) then bit reproducibility
!L  would be lost for short reruns.
!L  The energy correction applied in any day comes from the
!L  change calculated for the previous day. Initialisation for a NRUN
!L  sets the energy correction to zero for the first day.
!L
!L    A_REALHD(18) - total energy flux into the atmosphere
!L                   (used to accumulate change in energy throughout
!L                   a day from physics). Value at start of run zero
!L    A_REALHD(19) - total mass of the atmosphere (wet + dry)
!L                   Note this is not conserved by the dynamics
!L                   The model from UM 5.0 only conserves dry mass.
!L    A_REALHD(20) - total energy of the atmosphere calculated from
!L                   fields in dump.
!L    A_REALHD(21) - energy correction evaluated from previous day of
!L                   run (ie previous run or needs setting to zero).

      IF (STEPim(a_im) == 0 ) THEN


! DEPENDS ON: init_emcorr
         CALL INIT_EMCORR(                                              &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
            ICODE,CMESSAGE)


        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INIT_EMCORR'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

        END IF      ! end test on T+0
!
      END IF    !    LEMCORR

!------------------------------------------------------------------------
!  Section 4.8 - Initialise external halos
!                Some T+0 diagnostics and IAU calculations require
!                initialised values of external halos. LBCs will not
!                have been applied yet, so a simple fill of the appropriate
!                variables will suffice to ensure these calculations
!                use sensible values at the boundaries
!------------------------------------------------------------------------
      IF (model_domain == mt_lam) THEN
        CALL fill_external_halos(u, row_length, rows, model_levels,       &
                                 offx, offy)
        CALL fill_external_halos(v, row_length, n_rows, model_levels,     &
                                 offx, offy)
        CALL fill_external_halos(theta, row_length, rows, model_levels+1, &
                                 offx, offy)
        CALL fill_external_halos(rho, row_length, rows, model_levels,     &
                                 offx, offy)
      END IF

!------------------------------------------------------------------------
!  Section 4.9 - Initialise external halos
!                In LAM's there is the possibility that the extra moisture
!                fields are used but the LBC files do not hold values for
!                them.
!                This code will initialise the fields to sensible values
!                on the external halos in this case.
!------------------------------------------------------------------------
      IF (model_domain == mt_lam) THEN
        IF (l_mcr_qrain .AND. .NOT. l_mcr_qrain_lbc)                     &
        CALL set_external_halos(qrain, row_length, rows, wet_levels+1,   &
                                halo_i, halo_j, 0.0)
        IF (l_mcr_qgraup .AND. .NOT. l_mcr_qgraup_lbc)                   &
        CALL set_external_halos(qgraup, row_length, rows, wet_levels+1,  &
                                halo_i, halo_j, 0.0)
        IF (l_mcr_qcf2 .AND. .NOT. l_mcr_qcf2_lbc)                       &
        CALL set_external_halos(qcf2 ,row_length, rows, wet_levels+1,    &
                                halo_i, halo_j, 0.0)
      END IF

!L----------------------------------------------------------------------
!L 5. Set timestep group control switches for initial step
!L
!L 5.1 Set timestep control switches for initial step
!L

! DEPENDS ON: settsctl
      CALL SETTSCTL (                                                   &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),                                     &
                   internal_model,.TRUE.,meanlev,icode,cmessage)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to SETTSCTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF


!L
!L 5.2 Initialise PP files at step 0
!L

! DEPENDS ON: ppctl_init
      CALL PPCTL_INIT(                            &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),                                     &
         submodel,ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to PPCTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 5.3  Initialise assimilation package (not if assimilation completed)
!L
      IF ( (ASSIM_STEPSim(a_im)+ASSIM_EXTRASTEPSim(a_im) >  0  .AND.    &
         (model_assim_mode == "Atmosphere")      .AND.                  &
         (run_assim_mode   == "Atmosphere")      .AND.                  &
         STEPim(a_im)  <   ASSIM_FIRSTSTEPim(a_im) +                   &
         ASSIM_STEPSim(a_im) + ASSIM_EXTRASTEPSim(a_im))     &
         ) THEN
! DEPENDS ON: in_acctl
      CALL in_acctl(                                                    &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

     &                  ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to IN_ACCTL'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
      ENDIF

!----------------------------------------------------------------------
! [6.1]: Incremental Analysis Update (IAU).
!----------------------------------------------------------------------

      IF (L_IAU) CALL Setup_IAU

      IF (L_IAU .AND. STEPim(a_im) >= IAU_FirstCallTS .AND. &
                      STEPim(a_im) <= IAU_LastCallTS) THEN

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('IAU',3)
! DEPENDS ON: iau
        CALL iau (                                      &
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
                   l_mr_iau,                            & ! in
                   u,      v,          w,               & ! inout
                   u_adv,  v_adv,      w_adv,           & ! inout
                   theta,  rho,        murk,            & ! inout
                   q,      qCL,        qCF,             & ! inout
                   TStar,  TStar_tile, Deep_soil_temp,  & ! inout
                   smcl,   tstar_anom,                  & ! inout 
                   Pstar,  p,                           & ! inout
                   p_theta_levels,                      & ! inout
                   exner_rho_levels,                    & ! inout
                   exner_theta_levels,                  & ! inout
                   snodep,                              & ! inout
                   cf_area,                             & ! inout
                   cf_bulk,                             & ! inout
                   cf_liquid,                           & ! inout
                   cf_frozen,                           & ! inout
                   dust_div1, dust_div2, dust_div3,     & ! inout
                   dust_div4, dust_div5, dust_div6,     & ! inout
                   ozone_tracer )                         ! inout

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('IAU',4)

        ! If required, write out TS0 dump, including any IAU increments
        !  which may have just been added:
        IF (STEPim(a_im) == 0 .AND. L_IAU_DumpTS0State) THEN

! DEPENDS ON: dumpctl
          CALL DUMPCTL (                                                &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
     &                   atmos_sm, MEANLEV, .FALSE., '           ', 0,  &
     &                   ICODE, CMESSAGE)

          IF (ICODE /= 0) THEN
            WRITE (6,*) 'Failure in call to DUMPCTL (Atmos)'
            CALL Ereport(RoutineName, ICODE, CMESSAGE)
          END IF

        END IF

      END IF

!L----------------------------------------------------------------------
!L
!L 7.1  Get derived diagnostics from start fields (atmosphere)
!L
      IF (STEPim(a_im) == 0) THEN

! DEPENDS ON: initdiag
        CALL INITDIAG(                                                  &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
     & Dummyarg)


          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INITDIAG'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF
!
! 7.2 Code to update the boundaries now moved forwards to 4.4.1
!
!L
!L 7.3 Update ancillary fields in dump if start time corresponds to
!L     an ancillary field update time. Also done at T+0 with values
!L     updated to half a period back from first standard update time
!L     to ensure reproducibility between long runs and new runs
!L     started from dump at any time.
!L

      IF (ANCILLARY_STEPSim(a_im) >  0) THEN
        IF (STEPim(a_im) == 0 .OR.                                      &
           MOD(STEPim(a_im),ANCILLARY_STEPSim(a_im)) == 0)             &
! DEPENDS ON: up_ancil
           CALL UP_ANCIL (                                                &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Ancillary file arrays(atmosphere)-------------
      A_SPANC(A_IXANC( 1)),A_SPANC(A_IXANC( 2)),A_SPANC(A_IXANC( 3)),   &
      A_SPANC(A_IXANC( 4)),                                             &
     &                  submodel,                                       &
     &                  ICODE,CMESSAGE)
      ENDIF

          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to UP_ANCIL'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
!L
!L
!L 7.3.1 Initialize tiled prognostics, gridbox mean vegetation
!L       parameters and TRIFFID accumulation prognostics.
!L
!-----------------------------------------------------------------------
! Set l_aggregate appropriately. REVIEW FOR FLEXIBLE TILES?
!-----------------------------------------------------------------------
      IF(NTILES==1)THEN
        L_AGGREGATE=.TRUE.
      ELSE
        L_AGGREGATE=.FALSE.
      ENDIF

      IF (L_VEG_FRACS .AND. STEPim(a_im) == 0 ) THEN

!  Skip INIT_VEG if LAND_FIELD=0 for this PE.
        IF (LAND_FIELD  >   0) THEN
! DEPENDS ON: init_veg
          CALL INIT_VEG(STEPim(a_im),                                   &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

             ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INIT_VEG'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
        ELSE
          WRITE(6,*)'INITIAL; skip INIT_VEG, LAND_FIELD=0 for this PE'
        END IF

      END IF
!
!L
!L 7.3.3 Ensure that convective cloud cover and liquid water path
!L       are consistent with convective cloud base & top. (Corrects
!L       for occasional problems caused by reconfiguration.)
!L
      IF (STEPim(a_im) == 0) THEN

! DEPENDS ON: init_cnv
         CALL INIT_CNV(                                                 &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

            ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INIT_CNV'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF

      END IF
!
! River routing
       IF(L_RIVERS)THEN
! Initialise the step counter for river routing.
        IF (STEPim(a_im) == 0) THEN

! DEPENDS ON: init_riv
         CALL INIT_RIV(                                                 &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

            ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF

        END IF
      ENDIF                   ! L_RIVERS
!
!L
!L 7.3.4 ! initialization of radiative feedback
!L  Initialise and check address array to feed chemical tracers from
!L  UKCA into radiation scheme.  This needs to be done even for a CRUN
!L  so no need to check for STEPim(a_im) == 0.
!L
      grgas_addr = -1
      IF (L_UKCA) THEN

! DEPENDS ON: init_radukca
        CALL INIT_RADUKCA(                                             &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

           ngrgas,grgas_addr)

      END IF
!!
!L
!L 7.4 Generate interface fields at step zero if required
!L
      IF (linterface .AND. STEPim(a_im) == 0) THEN
        ndustbin_in  = 6 
        ndustbin_out = 6 
! DEPENDS ON: gen_intf
        CALL GEN_INTF (                                                 &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),                                     &
     &          submodel,ndustbin_in,ndustbin_out,ICODE,CMESSAGE) 

          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to GEN_INTF'
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF

! 7.5 Initialise JULES variables
!
! Check if land-surface fields are required.
      IF(land_field > 0) THEN
! DEPENDS ON: jules_init
        CALL jules_init(                                              &
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
                         land_field,ntiles,sm_levels                  &
                        ,nice,nice_use                                &
                        ,frac_typ                                     &
                        ,l_snow_albedo                                &
                        ,rgrain_tile,snodep_tile,deep_soil_temp       &
                        ,snow_grnd,nsnow,rgrainl,rho_snow_grnd        &
                        ,sice,sliq,snowdepth,ds,tsnowlayer  )
       ELSE
! No land field but we should still allocate these.  Better if we had split the
! routine to allocate arrays from initialising them but will do for now.
         ALLOCATE( clapp_levs(1,1))
         ALLOCATE( sathh_levs(1,1))
         ALLOCATE(  hcap_levs(1,1))
         ALLOCATE(  hcon_levs(1,1))
         ALLOCATE(satcon_levs(1,1))
         ALLOCATE(smvccl_levs(1,1))
         ALLOCATE(smvcwt_levs(1,1))
         ALLOCATE(smvcst_levs(1,1))
      END IF ! land-surface fields required

!L----------------------------------------------------------------------
!L 8. If coupled model, initialise addresses of coupling fields,
!L    and if model has restarted at the end of a coupling period
!L    exchange coupling fields and swap data (full ocean model)
!L    or both models are at step 0, exchange coupling fields and
!L    swap data (in sense O-A at step 0).
!L

      IF (L_OASIS) THEN

      ! If Running with OASIS3 or OASIS3-MCT we need to set up 
      ! pointers to relevant coupling fields. 


! DEPENDS ON: oasis_inita2o
      Call oasis_inita2o(                                               &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

         icode,                                            &
         cmessage)
      If (icode/=0) Then
        Write(6,*) 'Failure in call to inita2nemo'
        Call Ereport(RoutineName,icode,cmessage)
      Endif


      ENDIF

! Set up the ATMOS grid lat./long.coords for regridding to river
! routing grid
      IF(L_RIVERS)THEN
! DEPENDS ON: init_a2t_4A
        CALL INIT_A2T_4A(                                               &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(114)
      ! tpps_ozone: A_IXPTR(115)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(58)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),A_SPPTR(A_IXPTR(96)),  &
A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)),A_SPPTR(A_IXPTR(99)),                &
A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)),A_SPPTR(A_IXPTR(102)),             &
A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),A_SPPTR(A_IXPTR(105)),             &
A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),A_SPPTR(A_IXPTR(108)),             &
A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),A_SPPTR(A_IXPTR(111)),             &
A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),A_SPPTR(A_IXPTR(114)),             &
A_SPPTR(A_IXPTR(115)),                                                         &
      ! Biomass aerosol + lbc
A_SPPTR(A_IXPTR(116)),A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)), &
A_SPPTR(A_IXPTR(119)),A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)), &
A_SPPTR(A_IXPTR(122)),A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(125)),A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),&
A_SPPTR(A_IXPTR(128)),A_SPPTR(A_IXPTR(129)),A_SPPTR(A_IXPTR(130)),&
A_SPPTR(A_IXPTR(131)),A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(134)),A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),&
A_SPPTR(A_IXPTR(137)),A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),&
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
      ! Cloud fractions lbcs
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
! Add pointer for direct PAR flux
! Pointers displaced by 2 following removal of STOCHEM O3 and CH4
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
      ! Pointers for UKCA oxidant fields (162-165)
      ! Convective cloud fields (166-167)
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),                      &
      ! Pointers for aerosol climatologies (176-196)
A_SPPTR(A_IXPTR(178)),A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),&
A_SPPTR(A_IXPTR(181)),A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),&
A_SPPTR(A_IXPTR(184)),A_SPPTR(A_IXPTR(185)),A_SPPTR(A_IXPTR(186)),&
A_SPPTR(A_IXPTR(187)),A_SPPTR(A_IXPTR(188)),A_SPPTR(A_IXPTR(189)),&
A_SPPTR(A_IXPTR(190)),A_SPPTR(A_IXPTR(191)),A_SPPTR(A_IXPTR(192)),&
A_SPPTR(A_IXPTR(193)),A_SPPTR(A_IXPTR(194)),A_SPPTR(A_IXPTR(195)),&
A_SPPTR(A_IXPTR(196)),A_SPPTR(A_IXPTR(197)),A_SPPTR(A_IXPTR(198)),&
      !Fossil-fuel organic carbon aerosol + lbc - back 2
A_SPPTR(A_IXPTR(199)),A_SPPTR(A_IXPTR(200)),A_SPPTR(A_IXPTR(201)),&
A_SPPTR(A_IXPTR(202)),A_SPPTR(A_IXPTR(203)),A_SPPTR(A_IXPTR(204)),&
A_SPPTR(A_IXPTR(205)),A_SPPTR(A_IXPTR(206)),A_SPPTR(A_IXPTR(207)),&
      !Ammonium nitrate aerosol - back 2
A_SPPTR(A_IXPTR(208)),A_SPPTR(A_IXPTR(209)),A_SPPTR(A_IXPTR(210)),&
A_SPPTR(A_IXPTR(211)),A_SPPTR(A_IXPTR(212)),A_SPPTR(A_IXPTR(213)),&
A_SPPTR(A_IXPTR(214)),&
      !TKE based turbulence scheme - back 2
A_SPPTR(A_IXPTR(215)),A_SPPTR(A_IXPTR(216)),                      &
A_SPPTR(A_IXPTR(217)),                                            &
      !ENDGame
A_SPPTR(A_IXPTR(218)),A_SPPTR(A_IXPTR(219)),A_SPPTR(A_IXPTR(220)),&
A_SPPTR(A_IXPTR(221)),A_SPPTR(A_IXPTR(222)),A_SPPTR(A_IXPTR(223)),&
A_SPPTR(A_IXPTR(224)),A_SPPTR(A_IXPTR(225)),A_SPPTR(A_IXPTR(226)),&
A_SPPTR(A_IXPTR(227)),A_SPPTR(A_IXPTR(228)),A_SPPTR(A_IXPTR(229)),&
!{CABLE: pointers for introduced tiled prognostics to UM STASH, I/O  
A_SPPTR(A_IXPTR(230)),A_SPPTR(A_IXPTR(231)),A_SPPTR(A_IXPTR(232)),&
A_SPPTR(A_IXPTR(233)),A_SPPTR(A_IXPTR(234)),A_SPPTR(A_IXPTR(235)),&
A_SPPTR(A_IXPTR(236)),A_SPPTR(A_IXPTR(237)),A_SPPTR(A_IXPTR(238)),&
A_SPPTR(A_IXPTR(239)),                                            &
!}CABLE 

! ARTATCPL start
! Description: Super Arrays containing gridline coordinates for
! interpolation and area-averaging between atmosphere and
! river-routing grids (Part of ARTAOCPL.h)
!
!L --------------- (Atmosphere-TRIP) coupling arrays  -----------
!L ---------------Lat., Long. values of Atmosphere --------------
      AO_SPCPL(AO_IXCPL( 5)),AO_SPCPL(AO_IXCPL( 6)),                    &
      AO_SPCPL(AO_IXCPL( 7)),AO_SPCPL(AO_IXCPL( 8)),                    &
      AO_SPCPL(AO_IXCPL( 9)),AO_SPCPL(AO_IXCPL( 10)),                   &
! END ARTATCPL
           ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF
!L----------------------------------------------------------------------
!L 9. Print formatted diagnostics from initial dump
!L

! Printctl: printing of climate global and zonal diagnostics is no
!           longer supported

!L----------------------------------------------------------------------
!L 10. Initialisation complete - return to master routine
!L
! Check that operational model running has finished
! initialisation and write a message to the operator
      IF(mype == 0) THEN
         IF(model_status  ==  'Operational') THEN
! DEPENDS ON: operatormessage
            CALL OperatorMessage(nproc)
         ENDIF
      ENDIF
 999  CONTINUE

      IF (lhook) CALL dr_hook('INITIAL_4A',zhook_out,zhook_handle)
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE INITIAL_4A
