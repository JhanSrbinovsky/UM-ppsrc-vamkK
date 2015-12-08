! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: U_MODEL_4A (ENDGAME VERSION) -------------------------
!LL
!LL  Purpose: High level control program for the Unified Model
!LL           (master routine).  Calls lower level control routines
!LL           according to top level switch settings. Called by
!LL           top level routine UMSHELL which provides dimension sizes
!LL           for dynamic allocation of data arrays.
!LL
!LL  Programming standard: UM Doc Paper 3, version 8.3
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Top Level

      SUBROUTINE U_MODEL_4A(                                            &
     &       NFT,NFTU,                                                  &
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

! D1 replacement module

use atm_fields_mod
use earth_constants_mod
use atm_fields_bounds_mod
use integrity_mod

! OASIS Modules
USE oasis_atm_data_mod
  

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IO, ONLY             : io_timestep
USE model_file, ONLY : storeAllLookups, SynchroniseAll
USE um_types
USE filenamelength_mod, ONLY : filenamelength

! ensure that module variables cannot go out of scope
USE horiz_grid_mod
USE ref_pro_mod
USE departure_pts_mod
USE fields_rhs_mod
USE metric_terms_mod
USE helmholtz_const_matrix_mod
USE coriolis_mod
USE gravity_mod
USE eg_parameters_mod
USE wet_to_dry_n_calc_mod, ONLY : wet_to_dry_n

USE Control_Max_Sizes
USE decomp_DB
USE UM_ParVars
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod

! For managing dump deletion
USE MPPIO_job_control,        ONLY : jobcntrl     
USE MPPIO_job_control_common, ONLY : jc_delete

USE um_input_control_mod,  ONLY:                   &
     l_oasis,                                      &
     oasis_couple_freq
USE river_inputs_mod, ONLY: l_rivers
USE ppxlook_mod
USE acp_namel_mod, ONLY: l_ac
USE ukca_option_mod, ONLY: l_ukca

USE Submodel_Mod
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim


USE nlstcall_mod, ONLY : model_basis_time, &
                         LPP, &
                         ldump, &
                         lmean, &
                         lprint, &
                         linterface, &
                         lexit, &
                         ljobrelease, &
                         lancillary, &
                         lboundary, &
                         ltimer, &
                         ft_select

  USE chsunits_mod, ONLY : nunits

  use cable_data_mod, ONLY : set_endstep_umodel

IMPLICIT NONE

!*L  Interface and arguments: ------------------------------------------
!L       Sizes of super arrays
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
!L
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
!L
!L       Addresses of component arrays within super arrays
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
!L
!*----------------------------------------------------------------------
!
!  Common blocks
!
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
!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
! Interface variables initialised through INTFCNSTA
! Namelist read in the interface control routine INTF_CTL.

      INTEGER                                                           &
        intf_row_length                                                 &
                         ! Interface field row length
       ,intf_p_rows                                                     &
                         ! Interface field no of rows
       ,intf_p_levels                                                   &
                         ! Interface field no of levels
       ,intf_q_levels                                                   &
                         ! Interface field no of wet levels
       ,intf_tr_levels                                                  &
                         ! Interface field no of tracer levels
       ,intfwidtha                                                      &
                         ! Width of interface zone (atmosphere)
       ,intf_exthalo_ns                                                 &
                         ! Extended Halo in NS direction
       ,intf_exthalo_ew                                                 &
                         ! Extended Halo in EW direction
       ,a_intf_start_hr                                                 &
                         ! ) Start and End time in
       ,a_intf_freq_hr                                                  &
                         ! ) hours, Frequency in h,m,s for which
       ,a_intf_freq_mn                                                  &
                         ! ) atmosphere interface data
       ,a_intf_freq_sc                                                  &
                         ! ) is to be generated.
       ,a_intf_end_hr                                                   &
                         ! )
       ,intf_pack                                                       &
                         ! Packing Indicator for boundary data
       ,lbc_stream_a                                                    &
                         ! Output streams in UMUI
       ,lbc_unit_no_a                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
       ,lbc_first_r_rho                                                 &
                         ! First rho level at which height is constant
       ,intf_v_int_order(max_n_intf_a)

      REAL                                                              &
        intf_ewspace                                                    &
                         ! E-W grid spacing (degrees)
       ,intf_nsspace                                                    &
                         ! N-S grid spacing (degrees)
       ,intf_firstlat                                                   &
                         ! Latitude of first row (degrees)
       ,intf_firstlong                                                  &
                         ! Longitude of first row (degrees)
       ,intf_polelat                                                    &
                         ! Real latitude of coordinate pole (degrees)
       ,intf_polelong                                                   &
                         ! Real longitude of coordinate pole (degrees)
       ,lbc_z_top_model                                                 &
                         ! Height of top of model
       ,lbc_q_min                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , lambda_intf_p(max_intf_lbcrow_length, max_n_intf_a)             &
      , lambda_intf_u(max_intf_lbcrow_length, max_n_intf_a)             &    
      , phi_intf_p(max_intf_lbcrows, max_n_intf_a)                      &
      , phi_intf_v(max_intf_lbcrows, max_n_intf_a)

      LOGICAL                                                           &
        intf_vert_interp                                                &
                         ! Switch to request vertical interpolation
       ,lnewbnd          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  intf_l_var_lbc(max_n_intf_a)

! Switch to not rotate if input and output grids have same poles.
      LOGICAL intf_avoid_rot(MAX_N_INTF_A)

! Switch to output LBC for Endgame
      LOGICAL intf_l_eg_out(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      CHARACTER(LEN=256) :: intf_vertlevs

! Files for HorzGrid namelist  
      CHARACTER(LEN=256) :: intf_HorzGrid(max_n_intf_a)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
        intf_ewspace(max_n_intf_a)    ,intf_nsspace(max_n_intf_a),      &
        intf_firstlat(max_n_intf_a)   ,intf_firstlong(max_n_intf_a),    &
        intf_polelat(max_n_intf_a)    ,intf_polelong(max_n_intf_a),     &
        intf_row_length(max_n_intf_a) ,intf_p_rows(max_n_intf_a),       &
        intf_p_levels(max_n_intf_a)   ,intf_q_levels(max_n_intf_a),     &
        intf_tr_levels(max_n_intf_a)  ,intfwidtha(max_n_intf_a),        &
        intf_exthalo_ns(max_n_intf_a) ,intf_exthalo_ew(max_n_intf_a),   &
        a_intf_start_hr(max_n_intf_a) ,a_intf_freq_hr(max_n_intf_a),    &
        a_intf_freq_mn(max_n_intf_a)  ,a_intf_freq_sc(max_n_intf_a),    &
        a_intf_end_hr(max_n_intf_a)   ,                                 & 
        lnewbnd(max_n_intf_a)         ,intf_vert_interp(max_n_intf_a),  &
        intf_pack(max_n_intf_a)       ,lbc_stream_a(max_n_intf_a),      &
        lbc_unit_no_a(max_n_intf_a)   ,lbc_first_r_rho(max_n_intf_a),   &
        lbc_z_top_model(max_n_intf_a) ,                                 &
        intf_vertlevs(max_n_intf_a)   ,lbc_q_min,                       &
        intf_l_var_lbc                ,intf_horzgrid,                   &
        lambda_intf_p                 ,lambda_intf_u,                   &
        phi_intf_p                    ,phi_intf_v,                      &
        intf_avoid_rot                ,intf_v_int_order,                &
        intf_l_eg_out
!---------------------------------------------------------------------
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
!L
! Data kind parameters
!
!L
!L  DYNAMIC ALLOCATION OF SUPER ARRAYS:
!L
!L       Main D1 data array
! TYPSPD1 super array of D1 (NOTE: only single component array)
      REAL :: SPD1(SPD1_LEN)
! TYPSPD1 end
!L
!L       STASH related arrays
! TYPSPST super array of STASH arrays
      REAL :: SPSTS(SPSTS_LEN)
! TYPSPST end
!L
!L       Dump headers and lookups
! TYPSPDUA super array of dump arrays (atmosphere)
      REAL :: A_SPDUM(A_SPDUM_LEN)

      ! super array of auxilliary stash arrays (atmos)
      REAL :: A_SPSTS(A_SPSTS_LEN)
! TYPSPDUA end
!L
!L       Pointers (addresses) of model variables and constants
! TYPSPPTA super array of pointers (atmosphere)
      REAL :: A_SPPTR(A_SPPTR_LEN)
! TYPSPPTA end
!L Maximum sizes of fields limited by User Interface
!L
!L       Model derived constants arrays
! TYPSPCOA super array of constants arrays (atmosphere)
      REAL :: A_SPCON(A_SPCON_LEN)
! TYPSPCOA end
!L
!L       Generation of output interface fields
! TYPSPINA super array of output interface arrays (atmosphere)
      REAL :: A_SPINF(A_SPINF_LEN)
! TYPSPINA end
!L
!L       Updating of model from ancillary files
! TYPSPANA super array of ancillary file arrays (atmosphere)
      REAL :: A_SPANC(A_SPANC_LEN)
! TYPSPANA end
!L
!L       Boundary updating for Limited Area Models
! TYPSPBO  super array of input boundary arrays (not sub-model specific)
      REAL :: SPBND(SPBND_LEN)
! TYPSPBO end
! TYPSPBOA super array of input boundary arrays (atmosphere)
      REAL A_SPBND(A_SPBND_LEN)
! TYPSPBOA end
!L
!L       Coupled model arrays (atmosphere-river routing)
! TYPSPCPL super array of coupling arrays (atmosphere-ocean)
      REAL :: AO_SPCPL(AO_SPCPL_LEN)
! TYPSPCPL end
!L
!  Sizes for allocation of AC scheme arrays
! CSIZEOBS start
! The variables involved in dimensioning observational data arrays
! in the data assimilation section are stored in this comdeck.

      ! For ATMOSPHERE assimilations the values are computed in the
      ! initialisation routine INITAC and then passed to the main
      ! routine AC.

      INTEGER :: A_MAX_NO_OBS    !  No of observations in AC Obs files.
      INTEGER :: A_MAX_OBS_SIZE  !  No of obs values in AC Obs files.

      COMMON /CSIZEOBS/ A_MAX_NO_OBS, A_MAX_OBS_SIZE
! CSIZEOBS end
      INTEGER :: obs_flag_len,obs_len
      INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_FLAG
      REAL,    ALLOCATABLE, DIMENSION(:) :: OBS
!
!  Local variables
!
      INTEGER internal_model    ! Work - Internal model identifier
      INTEGER internal_model_prev!Work - Previous internal model ident
      INTEGER submodel          ! Work - Submodel id for dump partition
      INTEGER submodel_prev     ! Work - Previous submodel dump id
      INTEGER NGROUP            ! Work - Number of steps in "group"
      INTEGER MEANLEV           ! Work - Mean level indicator
      INTEGER IABORT            ! Work - Internal return code
      INTEGER I_STEP            ! Work - Loop counter over timesteps
      INTEGER G_THETA_FIELD_SIZE                                        &
                                   ! Sizes for MPP dynamic allocation
     &       ,G_IMTJMT          ! in A-O coupling routines

! Number of dust bins in LBC generation, disabled if not using MakeBC 
      INTEGER :: ndustbin_in, ndustbin_out  
!
! River routing
      INTEGER G_RIVER_FIELD_SIZE   ! Sizes for MPP dynamic allocation
!
      LOGICAL lexitNOW          ! Work - Immediate exit indicator
      INTEGER NFT           ! Unit no. for standard STASHmaster files
      INTEGER NFTU          ! Do. user STASH files (for GET_FILE)
      INTEGER RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER I,J,K,IDX     ! Loop counters
      INTEGER CO2_DIMA,                                                 &
                                   ! CO2 array dimensions
     &        CO2_DIMO,                                                 &
     &        CO2_DIMO2
      INTEGER DMS_DIMA,                                                 &
                                   ! DMS array dimensions
     &        DMS_DIMO,                                                 &
     &        DMS_DIMO2
      Integer info   ! Return code from GCom routines

      integer len_runid       !  No of chars in RUNID

      CHARACTER(LEN=4) runtype_char  !  Run Type (ie. NRUN, CRUN)

! 3-D fields of species to be passed down to radiation
      INTEGER, PARAMETER :: ngrgas = 8
      INTEGER, SAVE :: grgas_addr(ngrgas)
      
      LOGICAL :: put_step ! True when we're going to do a put
                      ! for OASIS purposes, false if we use dummy data.
      LOGICAL :: get_step ! True if we need to get data from OASIS, false
                      ! if we perform a get but discard results.
      LOGICAL :: cpl_update_step ! True when we need to make sure
                             ! the D1 prognostic data is updated
                             ! in preparation for dump creation or
                             ! coupling actions.


! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER(LEN=256) Cmessage    ! Error message
      CHARACTER(LEN=*) RoutineName

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      PARAMETER (   RoutineName='U_MODEL_4A')

! ENDGAME-only declarations
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1
      !First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150
      !Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
      !Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is -1.
      ! Similarly for A_UKCA_INDEX.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      ! A_TR_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: A_TR_StashItem(A_MAX_TRVARS)

      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)
      ! UKCA_tr_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: UKCA_tr_StashItem(A_MAX_UKCAVARS) 

      ! A_TR_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: A_TR_LBC_StashItem(A_MAX_TRVARS) 
      INTEGER :: A_TR_active_lbc_index(A_MAX_TRVARS) 

      ! UKCA_tr_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: UKCA_tr_LBC_StashItem(A_MAX_UKCAVARS) 
      INTEGER :: UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)

      COMMON/ATRACER/A_TR_INDEX, A_TR_StashItem,                        &
     &               A_TR_LBC_StashItem, A_TR_active_lbc_index,         &
     &               A_UKCA_INDEX, UKCA_tr_StashItem,                   &
     &               UKCA_tr_LBC_StashItem, UKCA_tr_active_lbc_index

! CTRACERA end

!
! ENDGame prognostic variables (not included in the start dump)
!



      REAL, parameter ::   inv_r_squared = 1. / Earth_Radius**2
! End of ENDGAME-only declarations

      IF (lhook) CALL dr_hook('U_MODEL_4A',zhook_in,zhook_handle)


      ALLOCATE(PPXI(ppxRecs,PPXREF_CODELEN))
      ALLOCATE(PPXPTR                                                    &
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS))


      ICODE=0
      CMESSAGE=''

!L----------------------------------------------------------------------
!L 0. Start Timer call for U_MODEL_4A (NB: not conditional on ltimer)
!L
! DEPENDS ON: timer
      IF (ltimer) CALL TIMER('U_MODEL_4A ',5)

! DEPENDS ON: lbc_coup_setup
      CALL lbc_coup_setup(ltimer)

      ICODE=0

!  Routine GETPPX_PART reads those ppxref records which correspond to
!  entries in the stash list into the ppx look-up arrays PPXI, PPXC.
!  It also sets the ppx pointer array PPXPTR. The lengths of PPXI, PPXC
!  have been dynamically allocated to the value of ppxRecs.

! Initialise row number in PPXI, PPXC arrays
      RowNumber = 1

! Initialise lookup and pointer array
      DO I=1,ppxRecs
        DO J=1,PPXREF_CODELEN
          PPXI(I,J)=0
        END DO
        DO J=1,PPXREF_CHARLEN
          PPXC(I,J) = ' '
        END DO
      END DO
      DO I = 1,N_INTERNAL_MODEL
        DO J   = 0,PPXREF_SECTIONS
          DO K = 1,PPXREF_ITEMS
            PPXPTR(I,J,K)=0
          END DO
        END DO
      END DO

! Read in STASHmaster records
      IF (INTERNAL_MODEL_INDEX(A_IM) >  0) THEN
! DEPENDS ON: getppx_part
      CALL GETPPX_PART(NFT,NFTU,'STASHmaster_A',A_IM,RowNumber,         &
     &                        ICODE,CMESSAGE)
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
      END IF

!L----------------------------------------------------------------------
!L 1. General initialisation of control and physical data blocks
!L

      ICODE=0

! DEPENDS ON: timer
      CALL timer('INITIAL ',5)
      CALL init_coriolis()

! DEPENDS ON: initial_4A
      CALL INITIAL_4A(                                                  &
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
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

! DEPENDS ON: timer
      CALL timer('INITIAL ',6)

! DEPENDS ON: icopenoutput
      CALL ICOPENOUTPUT(runtype_char)

!  Allocate AC scheme arrays using sizes from AC_INIT
      IF (L_AC) THEN
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
! 2048 gives enough space in WORK in RDOBS2 when no or very few obs.
        obs_len = A_MAX_OBS_SIZE+2048
        ALLOCATE (OBS(obs_len))
      WRITE(6,*)'U_MODEL_4A - OBS arrays allocated with sizes ',        &
     & A_MAX_NO_OBS,A_MAX_OBS_SIZE+2048
      ELSE
        A_MAX_NO_OBS = 1
        A_MAX_OBS_SIZE = 1
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
        obs_len = A_MAX_OBS_SIZE
        ALLOCATE (OBS(obs_len))
      WRITE(6,*)'U_MODEL_4A - OBS arrays allocated with length 1'
      END IF

!L----------------------------------------------------------------------
!L 2. Check for nothing-to-do
!L

! DEPENDS ON: exitchek
      CALL EXITCHEK( internal_model, lexitNOW)

      IF (LEXITNOW) GO TO 9999

!L----------------------------------------------------------------------
!L 3. Start group of timesteps

      IF (l_oasis) THEN

        put_step=.FALSE.
        get_step=.FALSE.
        cpl_update_step=.FALSE.

      END IF

      WRITE(6,'(A,F10.2,A)') 'Model running with timestep ',           & 
                             secs_per_stepim(a_im),' seconds' 

!L----------------------------------------------------------------------
!L 3. Start group of timesteps
!L
   1  CONTINUE
!L----------------------------------------------------------------------
!L 3.1. Start main timestep loop
!L
!L 3.1.1 Increment model time ..

! DEPENDS ON: incrtime
       CALL INCRTIME (                                                   &
!L --------------- Dump headers (atmosphere)-------------
      A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
      A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
      A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
      A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
      A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
           internal_model,ICODE,CMESSAGE)
       IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)


! Keep tabs on PRISM PUT/GET timesteps.
! At the moment we say a put and get timestep are one and the same.
      IF (l_oasis) THEN

! Is this a genuine exchange timestep.
        put_step = (MOD(stepim(atmos_im),oasis_couple_ts).eq.1)
        get_step = (MOD(stepim(atmos_im),oasis_couple_ts).eq.1)
        cpl_update_step=(MOD(stepim(atmos_im),oasis_couple_ts).eq.0)

! Perform coupling exchanges relating to TS-1

! DEPENDS ON: oasis3_geto2a
        CALL oasis3_geto2a(                                         &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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
           get_step)

! DEPENDS ON: oasis3_puta2o
        CALL oasis3_puta2o(                                         &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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
           put_step)

! Advance date ready for next timestep (if there is one)
! DEPENDS ON: OASIS3_ADVANCE_DATE
          CALL oasis3_advance_date()

      END IF ! l_oasis=true

!L 3.1.2 .. set timestep control switches
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
               internal_model,.FALSE.,meanlev,icode,cmessage)
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

!L 3.1.3 If PPfile initialisation time call PP control routine
!L          for instantaneous data (MEANLEV=0)
      IF (LPP) THEN

! DEPENDS ON: ppctl_reinit
        CALL PPCTL_REINIT(                                              &
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
           internal_model,ICODE,CMESSAGE)
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

! Send an EOT message to IOS
! He may or may not act on this to purge outstanding items
      CALL IO_TIMESTEP()


!L       Integrate atmosphere or ocean by 1 timestep
          IF (internal_model == atmos_im) THEN

! Synchronize before the timestep starts
      CALL GC_GSYNC(nproc,info)

! River routing
      IF (L_RIVERS) THEN
! Get 'global' atmos horizontal domain sizes from database
! in DECOMPDB to set dynamic allocation in ATM_STEP_4A for River routing
! on PE 0                                                   .
        G_THETA_FIELD_SIZE=                                             &
     &  decompDB(decomp_standard_atmos)%glsize(1,fld_type_p) *          &
     &  decompDB(decomp_standard_atmos)%glsize(2,fld_type_p)
        G_RIVER_FIELD_SIZE=                                             &
     &  decompDB(decomp_standard_atmos)%glsize(1,fld_type_r) *          &
     &  decompDB(decomp_standard_atmos)%glsize(2,fld_type_r)
      ELSE
        G_THETA_FIELD_SIZE=1
        G_RIVER_FIELD_SIZE=1
      END IF

!CABLE:
   call set_endstep_umodel( TARGET_END_STEPim(a_im) )

! DEPENDS ON: atm_step_4A
         CALL ATM_STEP_4A (                                             &
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
!L --------------- Input boundary arrays ------------------------
        SPBND(  IXBND( 1)),                                             &
!L --------------- Input boundary arrays(atmosphere)-------------
     &A_SPBND(A_IXBND( 1)),A_SPBND(A_IXBND( 2)),A_SPBND(A_IXBND( 3)),   &
     &A_SPBND(A_IXBND( 4)),A_SPBND(A_IXBND( 5)),                        &
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! ENDGame prognostic variables 
    exner,                                                              &
! River routing
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
     & G_THETA_FIELD_SIZE,                                              &
     & G_RIVER_FIELD_SIZE,                                              &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     & ngrgas,grgas_addr)
         
      IF (L_ukca) THEN
! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA_MAIN1',5)
! DEPENDS ON: ukca_main1
        CALL UKCA_MAIN1(                                                &
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

     &   I)                  ! dummy to terminate call
! DEPENDS ON: timer
        IF (ltimer) CALL TIMER('UKCA_MAIN1',6)
      END IF

!L Generate Atmosphere lateral boundary values

      IF (linterface) THEN

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER ('GEN_INTF',3)

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
     &       submodel, ndustbin_in, ndustbin_out, ICode, CMessage) 

! DEPENDS ON: timer
        IF (ltimer) CALL TIMER ('GEN_INTF',4)

        IF (ICode /= 0) THEN
          CALL Ereport(RoutineName, ICode ,Cmessage)
        END IF

      END IF  ! linterface

         END IF        ! internal_model = atmos_im

         IF (l_oasis) THEN

! Applicable to OASIS3 and OASIS3-MCT

           IF (cpl_update_step) THEN
!---------------------------------------------------------------------
! Ensure atmos coupling data in prognostic areas of D1 are up to date.
! The logic here is that cpl_update_step will be TRUE on the timestep
! BEFORE coupling is due to take place.

! Newly generated coupling data is intercepted after ATM_STEP and
! copied to  D1 prognostics.

! On the next timestep i.e. a coupling timestep, the prognostic
! contents will be sent to the other components in the coupling process.

! If there is no subsequent timestep (i.e. if this is the last model
! timestep then the D1 contents will be written to the dump
! ready for any future restart). There is no "end-of-model"
! coupling exchange.
!---------------------------------------------------------------------

! DEPENDS ON: oasis_updatecpl
             CALL oasis_updatecpl(                              &
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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
              cmessage)

           END IF
         END IF


!L 3.1.4 If dump time, call dump control routine
          IF (ldump) THEN
            IF (ltimer) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL',5)

            END IF
! DEPENDS ON: dumpctl
            CALL DUMPCTL (                                              &
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
     &          submodel,MEANLEV,.false.,'           ',0,               &
     &          ICODE,CMESSAGE)
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

            IF (ltimer) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL ',6)
            END IF
          END IF
!L 3.1.5 If printed output time, call print control routine
          IF (lprint) THEN

          IF (PrintStatus >= PrStatus_Oper) THEN
         WRITE(6,*) RoutineName,':Warning, Printing of climate global ' &
     &             ,'and zonal diagnostics no longer supported'
          END IF  ! PrintStatus test

          END IF
!L 3.1.6.1 Release job to process output created so far, if selected
          IF (ljobrelease) THEN
! Write all pp files' cached lookups to ensure that current data is 
! is in the files output stream
            CALL StoreAllLookups()

! Flush/Sync all pp files' units to ensure that data is committed from the 
! application
            CALL SynchroniseAll()

! DEPENDS ON: jobctl
            CALL JOBCTL(internal_model,ICODE,CMESSAGE)

            IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)
          END IF
!L 3.1.7 If partial sum/mean creation time, call means control routine
!L       (calls mean PPfield and diagnostic print routines internally)
          IF (lmean) THEN

! DEPENDS ON: meanctl
            CALL MEANCTL (                                              &
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

     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- Derived constants (atmosphere)--------
       A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),                                     &
     &                  submodel,MEANLEV,ICODE,CMESSAGE)

            IF (ICODE >  0) THEN
              CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF
          END IF

! 3.1.8 Update history file once dumping and any meaning is complete
  IF (ldump) THEN
! checkpoint_dump_im now successfully written. Update history file to restart
! from this instead of superseded dump stored in save_dumpname. Only delete 
! superseded dump when this has been done.

! DEPENDS ON: set_history_values
    CALL set_history_values
    ! Write history to backup location before overwriting main file
    IF (mype  ==  0) THEN
! DEPENDS ON: temphist
      CALL temphist(thist_unit,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(6,*)routinename,':Failure writing temporary restart file'
        WRITE(6,*)'Check for problems and restart from main file'
        CALL ereport(routinename,icode,cmessage)
      END IF
      CALL temphist(xhist_unit,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(6,*)routinename,':Failure writing main restart file'
        WRITE(6,*)'Check for problems and restart from temporary file'
        WRITE(6,*)'by overwriting main file with temporary file'
        CALL ereport(routinename,icode,cmessage)
      ELSE
        ! Main restart file successfully written, so delete backup
! DEPENDS ON: del_hist
        CALL del_hist(thist_unit)
      END IF

      IF (ft_select(22) == "Y") THEN
        ! History files written so can delete superseded dump.
        IF (save_dumpname_im(a_im) /= blank_file_name) THEN
          CALL jobcntrl(jc_delete,save_dumpname_im(a_im))
        END IF
        save_dumpname_im(a_im)=checkpoint_dump_im(a_im)
      END IF ! ft_select(22) == "Y"
    END IF ! IF (mype == 0)
  END IF ! IF (ldump)

!L 3.1.9 If exit check time, check for immediate exit
  IF (lexit) THEN

! DEPENDS ON: exitchek
    CALL EXITCHEK(internal_model, lexitNOW)

    IF (lexitNOW) THEN
      IF (.NOT.ldump) THEN

        WRITE(6,*)routinename,                                         &
           ': Warning: exiting at a period that is not a dump period'
        WRITE(6,*)'Therefore continuing the run will rerun preceding timesteps'
        WRITE(6,*)'This is inefficient and can cause restart problems'

      END IF
      GO TO 9999
    END IF
  END IF
!L 3.1.10 Update ancillary fields if necessary
          IF (lancillary) THEN

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
     &                   submodel,                                      &
     &                   ICODE,CMESSAGE)

      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
          END IF
!L 3.1.11 Update boundary fields if necessary
          IF (lboundary) THEN
             ! DEPENDS ON: lbc_coup_update
             CALL lbc_coup_update( &
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
                  submodel,ICODE,CMESSAGE)
          END IF
!L
!L      End main timestep loop
!L----------------------------------------------------------------------

      GO TO 1
!
 9999 CONTINUE

!L----------------------------------------------------------------------
!L 4. Exit processing: Output error messages and perform tidy-up
!L

      IF (l_oasis) THEN
! DEPENDS ON: oasis_tidy
        CALL oasis_tidy(                       &
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

        icode,cmessage)
      END IF


! DEPENDS ON: iccloseoutput
      CALL ICCLOSEOUTPUT()

!L 4.1 Exit processing: If abnormal completion, output error message
      IABORT = ICODE
      IF (ICODE /= 0) THEN

        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF
!L 4.2 Exit processing: Perform tidy-up

! DEPENDS ON: exitproc
      CALL EXITPROC(ICODE,CMESSAGE)

!L 4.3 Exit processing: If error in exit processing, output error mess
      IF (ICODE /= 0) THEN

        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      DEALLOCATE(PPXI)
      DEALLOCATE(PPXPTR) 

!L----------------------------------------------------------------------
!L 5. Complete Timer call and return
!L
      ICODE=IABORT
! DEPENDS ON: timer
      IF (ltimer) CALL TIMER('U_MODEL_4A ',6)

      IF (lhook) CALL dr_hook('U_MODEL_4A',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE U_MODEL_4A
