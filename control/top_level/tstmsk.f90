! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Checks version mask and option code in a ppx record
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Subroutine Interface:

SUBROUTINE tstmsk(modl,isec,lmask,ladres,errorstatus,cmessage)

! JULES
USE ancil_info, ONLY: nsmax

USE rimtypes
USE lbc_mod

USE g_wave_input_mod, ONLY:                                       &
    l_gwd,l_use_ussp

USE dust_parameters_mod, ONLY: l_dust, l_dust_diag, l_twobin_dust
USE um_input_control_mod, ONLY:                                            &
    l_use_arclocff,                                                        &
    l_soot   ,l_sulpc_so2   ,l_so2_natem  ,l_use_arcldust    ,             &
    l_nh3_em ,l_dms_em   ,l_sulpc_ozone ,l_sulpc_nh3  ,l_use_arclsulp    , &
    l_biomass,l_use_biogenic            ,              l_use_arclsslt    , &
    l_ocff   ,l_dms_ointer  ,l_sulpc_dms  ,l_use_arclblck    ,             &
    l_use_arclbiom,                                                        &
    l_nitrate,l_so2_surfem ,l_use_arcldlta    ,                            &
    l_oasis, l_oasis_icecalve, l_soot_surem,                               &
              l_triffid  ,               l_soot_hilem ,l_use_seasalt_pm  , &
                          l_bmass_surem ,l_bmass_hilem,                    &
                          l_so2_hilem  ,                                   &
    l_veg_fracs,l_ocff_surem  ,l_ocff_hilem ,l_co2_interactive ,           &
                          l_use_sulpc_indirect_sw,l_use_seasalt_direct   , &
                          l_sulpc_online_oxidants,l_use_seasalt_indirect , &
    l_endgame
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: l_cld_area, l_pc2, l_rhcpt
USE river_inputs_mod, ONLY: l_rivers, l_inland
USE cv_run_mod, ONLY: l_ccrad, l_3d_cca, l_conv_hist, l_param_conv
USE murk_inputs_mod,  ONLY: l_murk, l_murk_source

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ancilcta_namelist_mod, ONLY: l_sstanom
USE umsections_mod, ONLY: atmos_sr
USE bl_option_mod,  ONLY: isrfexcnvgust
USE rad_input_mod, ONLY: l_use_cariolle, l_use_tpps_ozone, l_use_grad_corr,&
     l_orog_unfilt, i_ozone_int

! JULES
USE switches, ONLY: l_flake_model, l_spec_albedo, l_albedo_obs, iscrntdiag,&
                    l_snow_albedo, l_sice_multilayers, l_top, can_model,   &
                    l_ctile, l_aggregate, i_aggregate_opt
USE switches_urban, ONLY : l_urban2t

USE PrintStatus_mod
USE UM_ParParams
USE domain_params
USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon,      &
                      io3_trop_map, io3_trop_map_masscon
! version_mod items required by cstash.h and model.h
USE version_mod, ONLY: nproftp, nprofdp, nprofup, ndiagpm,        &
                       ntimep, NTimSerP, nlevp, npslevp,          &
                       npslistp, outfile_s, outfile_e, nsectp


USE Submodel_Mod
IMPLICIT NONE

! Description:
!   Determines whether a diagnostic is available to a particular
!   (internal model,section) version. Also checks the option code IOPN
!   Called by INACTR, PRELIM.

! Method:

! The decimal value of the version mask, VMSK, was read from the
! ppxref file by GETPPX, each (model, section, item). The version
! number of the (model, section) - NMASK - is obtained from the
! H_VERS array.

! The procedure for checking diagnostic availability is as follows:
! (1) Check whether the relevant internal model is included in the
!     submodel configuration, by examining the INTERNAL_MODEL_LIST
!     array; if not, return.
! (2) Check whether NMASK=0 - this implies that the diag is
!      unavailable to any version. If so, return.
! (3) Check whether the diag is available to the specified version.
!     The 'version mask' binary code in the STASHmaster specifies
!     which versions the diag is available to. Eg., if it is available
!     to vns. 1 and 3, but not vn.2, the version mask is 101. VMSK
!     is the decimal equivalent of the version mask (5 in this example).
!     TSTMSK therefore calculates the quantity

!             IMOD = MOD(VMSK,2**NMASK)/2**(NMASK-1)

!       If IMOD=1, the diag is available; if IMOD=0, it isn't.

!   The option code is checked in accordance with the ppxref option
!   code definitions.

!   Note for code developers adding new tests: LMASK is initialised to
!   .TRUE. at top of deck. A series of tests is made and if any one
!   test fails LMASK is reset to .FALSE. Therefore do not reinitialise
!   LMASK to .TRUE. anywhere otherwise you may inadvertently overwrite
!   a preceding .FALSE. setting or you may mislead future code writers
!   who may do so.
!    For this reason, all unnecessary LMASK=.TRUE. were removed at 5.1

!  Code description:
!  Language: Fortran 90.
!  This code is written to UM programming standards version 8.3.

!  System component covered:
!  System task:               Sub-Models Project

! Global variables:
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
! MODEL Defines model-dependent quantities used by data addressing and
! STASH
!
! submodel_mod must be used before this one
! VERSION_MOD module is required for nsectp, outfile_s and outfile_e
!
      INTEGER, PARAMETER :: AASSETS    = 9
      INTEGER, PARAMETER :: MEAD_TYPES = 4
      INTEGER, PARAMETER :: A_MAX_TRVARS=150 !Max.no.of tracers allowed
      INTEGER, PARAMETER :: A_MAX_UKCAVARS=150 ! Max.no.of UKCA allowed
      INTEGER, PARAMETER :: MAX_AOBS=100

      REAL :: H_A_EWSPACE
      REAL :: H_A_NSSPACE
      REAL :: H_A_FIRSTLAT
      REAL :: H_A_FIRSTLONG
      REAL :: H_A_POLELAT
      REAL :: H_A_POLELONG

      INTEGER :: H_A_GROUP
      INTEGER :: H_OROG_ROUGH
      INTEGER :: A_ASSMGRPS
      INTEGER :: NUM_PVPR

      LOGICAL :: A_RECON
      LOGICAL :: H_OROG_GRAD
      LOGICAL :: ATMODS
      LOGICAL :: CMODS
      LOGICAL :: LMESO

      LOGICAL :: TRACER_A (0:A_MAX_TRVARS)
      LOGICAL :: TR_UKCA_A (0:A_MAX_UKCAVARS)
      LOGICAL :: AASSET   (AASSETS)
      INTEGER :: AASPF    (AASSETS)
      INTEGER :: AASPL    (AASSETS)
      INTEGER :: AOBINC   (MAX_AOBS)
      INTEGER :: AOBGRP   (MAX_AOBS)
      INTEGER :: RUN_TARGET_END( 6)

      COMMON/MODELA/ H_A_EWSPACE,H_A_NSSPACE,H_A_FIRSTLAT,H_A_FIRSTLONG,&
     &  H_A_POLELAT,H_A_POLELONG,A_ASSMGRPS,NUM_PVPR ,A_RECON,H_A_GROUP,&
     &  H_OROG_GRAD,ATMODS,CMODS,LMESO,TRACER_A,TR_UKCA_A,              &
     &  AASSET,AASPF,AASPL

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX)
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/                                                     &
     &  LPRIM,LDUMP,LSECD,LWORK,NHEAD,LEXTRA,LPRIM_O2,LPrimIM,LDumpIM,  &
     &  LSecdIM,NHeadSub,ITEM_MAX_REQ,ITEM_MAX_ALL,NSERBLK_S,NSERREC_S, &
     &  NLEVL_S,NMAXLEV_S,NPSLISTS_S,NMAXPSL_S,LSTUSER,NRECS_S,NTIMES_S,&
     &  NHEAD_FILE,                                                     &
     &  global_LPRIM,global_LPrimIM,global_LDUMP,global_LDumpIM
      CHARACTER(LEN=1)  H_ATMOS
      CHARACTER(LEN=1)  H_FLOOR
      CHARACTER(LEN=1)  H_STRAT
      CHARACTER(LEN=1)  H_GLOBAL(N_INTERNAL_MODEL_MAX         )
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/ H_ATMOS,H_GLOBAL,H_FLOOR,H_STRAT

      COMMON/HVERS/ H_VERS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER


! Variables read in by namelist and used in SETMODL
      INTEGER      OCAAA 
      REAL         EWSPACEA,NSSPACEA
      REAL         FRSTLATA,FRSTLONA

      LOGICAL      ZonAvOzone
      LOGICAL      ZonAvTppsOzone
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      OCALB
      REAL         POLELATA
      REAL         POLELONA
      INTEGER      SWBND
      INTEGER      TCA(A_MAX_TRVARS)
      INTEGER      TCA_LBC(A_MAX_TRVARS)  ! =1 if tracer in lbc file 
      INTEGER      TC_UKCA(A_MAX_UKCAVARS)
      INTEGER      TC_LBC_UKCA(A_MAX_UKCAVARS) ! =1 if tr in lbc file 
      INTEGER      StLevGWdrag
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/                                                  &
     &  RUN_TARGET_END,                                                 &
     &  OCAAA,EWSPACEA,POLELATA,FRSTLATA,LATS,                          &
     &  NSSPACEA,POLELONA,FRSTLONA,LONS,                                &
     &  SWBND,LWBND,                                                    &
     &  ZonAvOzone ,ZonAvTppsOzone, AOBINC,  AOBGRP,                    &
     &  StLevGWdrag, BotVDiffLev,TopVDiffLev,                           &
     &  OCALB,TCA,TCA_LBC,TC_UKCA,TC_LBC_UKCA


      CHARACTER(LEN=1) :: LFLOOR
      CHARACTER(LEN=1) :: OROGR
      CHARACTER(LEN=1) :: SWMCR
      CHARACTER(LEN=1) :: MESO

      COMMON/STSHCHAR/                                                  &
     &     LFLOOR,                                                      &
     &  OROGR,   SWMCR, MESO

      NAMELIST/STSHCOMP/                                                &
        RUN_TARGET_END,                                                 &
        OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA  ,LATS   ,         &
                     NSSPACEA    ,POLELONA ,FRSTLONA  ,LONS   ,         &
        SWBND       ,LWBND                            ,OROGR  ,         &
        ZonAvOzone  ,SWMCR       ,MESO     ,                            &
        OCALB       ,LFLOOR      ,AOBINC   ,TCA,                        &
        TCA_LBC     ,TC_UKCA     ,TC_LBC_UKCA   ,AOBGRP          
  
! MODEL end

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER modl    ! Internal model number
INTEGER isec    ! Section number

!   Scalar arguments with intent(out):
LOGICAL lmask   ! T if diag is available to specified version
LOGICAL ladres  ! T if diag is available for addressing primary
CHARACTER(LEN=80) cmessage

! Local scalars
INTEGER nmask   ! Version number for (model,section)
INTEGER imod    ! Determines whether diag is available
INTEGER twonm   ! Used in calculation of IMOD
INTEGER twonm1  ! Used in calculation of IMOD
INTEGER count
INTEGER i
INTEGER n0,n1,n2,n2n1,n3,n4,n5,n6,n7,n8,n9,n10
INTEGER n11,n12,n13,n14,n15,n16,n17,n18,n19,n20
INTEGER n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
INTEGER n10n9
INTEGER sum_iopn

! Hard-wire var_recon to false if not the reconfiguration
LOGICAL, PARAMETER :: var_recon=.FALSE.
! ErrorStatus
INTEGER errorstatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook('TSTMSK',zhook_in,zhook_handle)
sum_iopn=iopn(1)+iopn(2)+iopn(3)+iopn(4)+iopn(5)+iopn(6)
!--------------------------------------------------------------------
! Check whether internal model is included in submodel configuration
!--------------------------------------------------------------------
count = 0
DO i = 1,n_internal_model_max
  IF (internal_model_list(i) == modl) THEN
    lmask =.TRUE.
    ladres=.TRUE.
    count = count + 1
  END IF
END DO
IF (count == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
END IF

!--------------------------------------------------------------------
! Check whether diagnostic is unavailable to any version
!--------------------------------------------------------------------
nmask=h_vers(modl,isec)
IF (nmask == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
END IF

! Determine whether the diag is available to the specified version
twonm  = 2**nmask
twonm1 = 2**(nmask-1)
imod   = MOD(vmsk,twonm)/twonm1
IF(imod == 1) THEN
  lmask =.TRUE.
ELSE IF(imod == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
ELSE
  WRITE(6,*)'S: TSTMSK INVALID DECODING OF VMSK',vmsk
  WRITE(6,*)'s: ... imod=',imod,'     NMASK=',nmask
END IF

!-------------------------------------------------------------------
! Check option codes
!-------------------------------------------------------------------
lmask=.TRUE.
!-------------------------------------------------------------------
! Var reconfiguration now is not restricted to primary fields.
! Therefore need to perform on all sections.
!-------------------------------------------------------------------
IF ( var_recon ) THEN
  n17 = MOD((iopn(4)/10),10)
  IF (n17 /= 1) THEN
    lmask = .FALSE.
  END IF
END IF

IF((isec == 0) .OR. (isec == 31) .OR.                           &
   (isec == 32) .OR. (isec == 33) .OR.                          &
   (isec == 34) .OR. (isec == 36) .OR. (isec == 37)) THEN
! Atmosphere primary field
  IF(sum_iopn /= 0) THEN
    n2n1=MOD (iopn(1),100)
! Up to A_MAX_TRVARS=150 free tracers now allowed in section 33
! Up to A_MAX_UKCAVARS=150 UKCA tracers now allowed in section 34
! Up to 150 free tracer lbcs in section 36
! Up to 150 UKCA tracer lbcs in section 37
    IF (isec == 33 .OR. isec == 34 .OR. isec == 36              &
   .OR. isec == 37) THEN
      n2n1=MOD (iopn(1),1000)
    END IF
    n3  =MOD((iopn(1)/100),10)
    n4  =MOD((iopn(1)/1000),10)
    n5  =MOD((iopn(1)/10000),10)
    n6  =MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    n8  =MOD((iopn(2)/100),10)
! n10n9 is a 2 digit option code for CLASSIC aerosols
    n10n9 =MOD((iopn(2)/1000),100)
    n11 =MOD( iopn(3),10)
    n12 =MOD((iopn(3)/10),10)
    n13 =MOD((iopn(3)/100),10)
    n14 =MOD((iopn(3)/1000),10)
    n15 =MOD((iopn(3)/10000),10)
    n16 =MOD( iopn(4),10)
    n17 =MOD((iopn(4)/10),10)
    n18 =MOD((iopn(4)/100),10)
    n19 =MOD((iopn(4)/1000),10)
    n20 =MOD((iopn(4)/10000),10)
    n21 =MOD( iopn(5),10)
    n22 =MOD((iopn(5)/10),10)
    n23 =MOD((iopn(5)/100),10)
    n24 =MOD((iopn(5)/1000),10)
    n25 =MOD((iopn(5)/10000),10)
    n26 =MOD( iopn(6),10)
    n27 =MOD((iopn(6)/10),10)
    n28 =MOD((iopn(6)/100),10)
    n29 =MOD((iopn(6)/1000),10)
    n30 =MOD((iopn(6)/10000),10)
    IF ((n2n1 == 99) .AND.                                    &
        (isec /= 33) .AND. (isec /= 34) .AND.                 &
        (isec /= 36) .AND. (isec /= 37)) THEN
      lmask=.FALSE.

    ELSE IF (isec == 33 .OR. isec == 36) THEN
      IF (n2n1 > 0 .AND. n2n1 <= a_max_trvars) THEN
        IF ((isec == 33 .AND. .NOT. tracer_a(n2n1)) .OR.      &
            (isec == 36 .AND. tca_lbc(n2n1) == 0)) THEN
          lmask=.FALSE.
        END IF
      END IF

    ELSE IF (isec == 34 .OR. isec == 37) THEN
      IF (n2n1 > 0 .AND. n2n1 <= a_max_ukcavars) THEN
        IF ((isec == 34 .AND. .NOT. tr_ukca_a(n2n1)) .OR.     &
            (isec == 37 .AND. tc_lbc_ukca(n2n1) == 0)) THEN
          lmask=.FALSE.
        END IF
      END IF

! Make n30 a metacharacter to allow more option codes
    ELSE IF (n30 == 0) THEN

!             n3=1 not to be used because of tracers
!             n3=2 not to be used
      IF((n3 == 5).AND.(.NOT.l_top)) THEN
        lmask=.FALSE.   ! TOPMODEL hydrology
      ELSE IF((n4 == 1).AND.(.NOT.l_use_biogenic)) THEN
        lmask=.FALSE.   ! biogenic aerosol
      ELSE IF((n4 == 2).AND.(.NOT.l_use_arclbiom)) THEN
        lmask=.FALSE.   ! biomass burning aerosol
      ELSE IF((n4 == 3).AND.(.NOT.l_use_arclblck)) THEN
        lmask=.FALSE.   ! black carbon aerosol
      ELSE IF((n4 == 4).AND.(.NOT.l_use_arclsslt)) THEN
        lmask=.FALSE.   ! sea salt aerosol
      ELSE IF((n4 == 5).AND.(.NOT.l_use_arclsulp)) THEN
        lmask=.FALSE.   ! sulphate aerosol
      ELSE IF((n4 == 6).AND.(.NOT.l_use_arcldust)) THEN
        lmask=.FALSE.   ! dust aerosol
      ELSE IF((n4 == 7).AND.(.NOT.l_use_arclocff)) THEN
        lmask=.FALSE.   ! organic carbon fossil fuel
      ELSE IF((n4 == 8).AND.(.NOT.l_use_arcldlta)) THEN
        lmask=.FALSE.   ! delta aerosol
      ELSE IF((n4 > 8)) THEN
        lmask=.FALSE.  ! these are not yet defined
      ELSE IF((n5 == 1).AND.(orogr == 'N')) THEN
        lmask=.FALSE.   ! orographic roughness
      ELSE IF((n5 == 2).AND.(.NOT.h_orog_grad)) THEN
        lmask=.FALSE.   ! orographic gradient
      ELSE IF((n5 == 3).AND.(.NOT.l_use_grad_corr)) THEN
        lmask=.FALSE.   ! gradient correction for SW radiation
      ELSE IF((n5 == 4).AND.(.NOT.l_orog_unfilt)) THEN
        lmask=.FALSE.   ! unfiltered orography for horizon angles
      ELSE IF((n6 == 1) .AND.                                   &
          (rimwidtha(rima_type_norm) == 0 .OR. h_global(a_im) == 'Y')) THEN
        lmask=.FALSE.   ! limited area model boundary condition
      ELSE IF((n6 == 2).AND.(h_floor == 'N')) THEN
        lmask=.FALSE.   ! lower boundary - growing orography (redundant)
      ELSE IF((n7 == 1).AND.(.NOT.l_oasis)) THEN
        lmask=.FALSE.   ! OASIS coupling to ocean model
      ELSE IF((n7 == 2)) THEN
        lmask=.FALSE.   ! other coupling fields currently excluded
      ELSE IF((n7 == 3).AND.(.NOT.(l_oasis.AND.l_oasis_icecalve))) THEN
        lmask=.FALSE.   ! oasis iceberg calving
      ELSE IF((n7 == 7).AND.(.NOT.(l_dms_ointer))) THEN
        lmask=.FALSE.   ! coupling for DMS ocean flux
      ELSE IF((n8 == 1).AND.(.NOT.l_sstanom)) THEN
        lmask=.FALSE.   ! SST anomaly
      ELSE IF((n8 == 2).AND.(iscrntdiag /= 2)) THEN
        lmask=.FALSE.
      ELSE IF((n8 == 3).AND.(isrfexcnvgust == 0)) THEN
        lmask=.FALSE.   ! effect of convective downdraughts on surface exchange
      ELSE IF((n8 == 4).AND.(.NOT.l_murk)) THEN
        lmask=.FALSE.   ! total aerosol (murk) for visibility
      ELSE IF((n8 == 5).AND.(.NOT.l_murk_source)) THEN
        lmask=.FALSE.   ! total aerosol emissions
      ELSE IF((n8 == 6).AND.(.NOT.l_snow_albedo)                &
                       .AND.(nsmax == 0)) THEN
        lmask=.FALSE.   ! snow albedo
      ELSE IF((n8 == 7).AND.(atmos_sr(3) /= '1A')) THEN
        lmask=.FALSE.   ! tke closure
      ELSE IF ( ( n8  ==  8) .AND. (atmos_sr(14) == '0A') ) THEN
        lmask=.FALSE.   ! energy adjustment scheme

!
! n10n9 is used for all CLASSIC aerosol related prognostics
! Sulphur cycle (1-19)
      ELSE IF((n10n9 == 1).AND.(.NOT.l_sulpc_so2)) THEN
        lmask=.FALSE.   ! sulphur dioxide cycle
      ELSE IF((n10n9 == 2).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_so2_surfem)) ) THEN
        lmask=.FALSE.   ! surface SO2 emissions
      ELSE IF((n10n9 == 3).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_so2_hilem)) ) THEN
        lmask=.FALSE.   ! high level SO2 emissions
      ELSE IF((n10n9 == 4).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_so2_natem)) ) THEN
        lmask=.FALSE.   ! natural SO2 emissions
      ELSE IF((n10n9 == 5).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_sulpc_dms)) ) THEN
        lmask=.FALSE.   ! dimethyl sulphide in SO2 cycle
      ELSE IF((n10n9 == 6).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_sulpc_dms)                  &
                        .OR.  (.NOT.l_dms_em)) ) THEN
        lmask=.FALSE.   ! dimethyl sulphide emissions
      ELSE IF((n10n9 == 7).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (l_sulpc_online_oxidants)           &
                        .OR.  (.NOT.l_sulpc_ozone)) ) THEN
        lmask=.FALSE.   ! offline ozone oxidant in SO2 cycle
      ELSE IF((n10n9 == 8).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_sulpc_ozone)                &
                        .OR.  (.NOT.l_sulpc_nh3)) )  THEN
        lmask=.FALSE.   ! ozone and ammonia in SO2 cycle
      ELSE IF((n10n9 == 9).AND.((.NOT.l_sulpc_so2)                &
                        .OR.  (.NOT.l_sulpc_ozone)                &
                        .OR.  (.NOT.l_sulpc_nh3)                  &
                        .OR.  (.NOT.l_nh3_em)) )  THEN
        lmask=.FALSE.   ! ammonia emissions and O3 in SO2 cycle
      ELSE IF((n10n9 == 10).AND.((.NOT.l_sulpc_so2)               &
                        .OR.  (l_sulpc_online_oxidants)) )  THEN
        lmask=.FALSE.   ! offline oxidants (not ozone)
! Soot (21 - 23)
      ELSE IF((n10n9 == 21).AND.(.NOT.l_soot))  THEN
        lmask=.FALSE.   ! soot scheme
      ELSE IF((n10n9 == 22).AND.((.NOT.l_soot)                    &
                        .OR.  (.NOT.l_soot_surem)) )  THEN
        lmask=.FALSE.   ! soot scheme with surface emissions
      ELSE IF((n10n9 == 23).AND.((.NOT.l_soot)                    &
                       .OR.  (.NOT.l_soot_hilem)) )  THEN
        lmask=.FALSE.   ! soot scheme with high level emissions
! Biomass (24 - 26)
      ELSE IF ((n10n9 == 24).AND.(.NOT.l_biomass)) THEN
        lmask=.FALSE.   ! biomass scheme
      ELSE IF ((n10n9 == 25).AND.((.NOT.l_biomass)                 &
                       .OR. (.NOT.l_bmass_surem)) ) THEN
        lmask=.FALSE.   ! biomass scheme with surface emissions
      ELSE IF ((n10n9 == 26).AND.((.NOT.l_biomass)                 &
                         .OR. (.NOT.l_bmass_hilem)) ) THEN
        lmask=.FALSE.   ! biomass scheme with high level emissions
! 2 bin and 6 bin Mineral dust
      ELSE IF ((n10n9 == 27).AND.(.NOT.l_dust)) THEN
        lmask=.FALSE.   ! mineral dust scheme (prognostic) 
! 6 bin (and currently also 2 bin) mineral dust diagnosis
      ELSE IF ((n10n9 == 28).AND.(.NOT.l_dust).AND.(.NOT.l_dust_diag)) THEN  
        lmask=.FALSE.   ! mineral dust scheme (diagnostic lifting only) 
! Nitrate (31)
      ELSE IF ((n10n9 == 31).AND.(.NOT.l_nitrate)) THEN
        lmask=.FALSE.   ! nitrate scheme
! OCFF (35-37)
      ELSE IF ((n10n9 == 35).AND.(.NOT.l_ocff)) THEN
        lmask=.FALSE.   ! organic carbon fossil fuel scheme
      ELSE IF ((n10n9 == 36).AND.((.NOT.l_ocff)                    &
                         .OR. (.NOT.l_ocff_surem)) ) THEN
        lmask=.FALSE.   ! OCFF with surface emissions
      ELSE IF ((n10n9 == 37).AND.((.NOT.l_ocff)                    &
                         .OR. (.NOT.l_ocff_hilem)) ) THEN
        lmask=.FALSE.   ! OCFF with high level emissions
! 6 bin only Mineral dust
      ELSE IF ((n10n9 == 38).AND.(l_twobin_dust                    &
                         .OR. (.NOT.l_dust)) ) THEN
        lmask=.FALSE.   ! mineral dust scheme (prognostic) 
! End of classic aerosol section

      ELSE IF (n11 == 1) THEN
        lmask=.FALSE.   ! basic vegetation scheme, no longer supported
      ELSE IF((n11 == 2).AND.(.NOT.l_veg_fracs)) THEN
        lmask=.FALSE.   ! fractional vegetation scheme
      ELSE IF((n11 == 3).AND.                                   &
              (.NOT.l_veg_fracs.OR..NOT.l_triffid) ) THEN
        lmask=.FALSE.   ! TRIFFID vegetation scheme
      ELSE IF((n11 == 4).AND.((.NOT.(l_veg_fracs.AND.l_albedo_obs)).OR.      &
              l_spec_albedo)) THEN
        lmask=.FALSE.   ! SW albedo_obs
      ELSE IF((n11 == 5).AND.(.NOT.(l_veg_fracs.AND.l_albedo_obs.AND.        &
              l_spec_albedo))) THEN
        lmask=.FALSE.   ! VIS and NIR albedo_obs
      ELSE IF( (n11 == 6) .AND.                                              &
        ( .NOT.(l_aggregate.AND.(i_aggregate_opt==1)) ) ) THEN                
        lmask=.FALSE.   ! No separate prognostic for thermal roughness
      ELSE IF((n12 == 3) .AND. (.NOT.l_mcr_qcf2)) THEN
        lmask=.FALSE.   ! QCF2 is prognostic
      ELSE IF((n12 == 4) .AND. (.NOT.l_mcr_qrain)) THEN
        lmask=.FALSE.   ! QRAIN is prognostic
      ELSE IF((n12 == 5) .AND. (.NOT.l_mcr_qgraup)) THEN
        lmask=.FALSE.   ! QGRAUP is prognostic
      ELSE IF((n12 == 6) .AND. (.NOT.l_pc2)) THEN
        lmask=.FALSE.   ! PC2 clod fraction boundary values out (secn.32)
      ELSE IF ((n13 == 1).AND.(l_3d_cca)) THEN
        lmask=.FALSE.    !  convective cloud amount is 2D
      ELSE IF ((n13 == 2).AND.(.NOT.l_3d_cca)) THEN
        lmask=.FALSE.    !  convective cloud amount is 3D
      ELSE IF ((n13 == 3).AND.(.NOT.l_ccrad)) THEN
        lmask=.FALSE.   ! CCRad scheme
      ELSE IF ((n13 == 5).AND.(.NOT.(l_conv_hist))) THEN
        lmask=.FALSE.   ! fields for convection history
      ELSE IF ((n15 == 3).AND.(.NOT.l_co2_interactive)) THEN
        lmask=.FALSE.   ! carbon cycle scheme
      ELSE IF ((n16 == 2).AND.(nsmax==0)) THEN
        lmask=.FALSE.   ! JULES snow scheme with multiple snow layers
      ELSE IF ((n16 == 3).AND.(.NOT.l_snow_albedo)) THEN
        lmask=.FALSE.   ! snow soot
      ELSE IF ((n16 == 4).AND.(.NOT.l_urban2t)) THEN
        lmask=.FALSE.   ! URBAN-2T schemes including MORUSES
      ELSE IF ((n16 == 5).AND.(.NOT.l_flake_model)) THEN
        lmask=.FALSE.   ! FLake lake scheme
! n17 is left out here, as it is for all sections not just prognostics.
      ELSE IF ((n18 == 1).AND.(.NOT.l_ctile)) THEN
        lmask=.FALSE.   ! coastal tiling scheme
      ELSE IF ((n18 == 3).AND.(.NOT.l_rivers)) THEN
        lmask=.FALSE.   ! river routing scheme
      ELSE IF ((n18 == 4).AND.(.NOT.l_inland)) THEN
        lmask=.FALSE.   ! inland re-routing scheme
      ELSE IF ((n18 == 6).AND.(nice_use > 1)) THEN
        lmask=.FALSE.   ! single sea ice category
      ELSE IF ((n18 == 7).AND.(nice == 1)) THEN
        lmask=.FALSE.   ! sea ice categories
      ELSE IF ((n18 == 8).AND.(nice_use == 1)) THEN
        lmask=.FALSE.   ! sea ice categories used fully 
      ELSE IF ((n18 == 9).AND.(.NOT.l_sice_multilayers)) THEN
        lmask=.FALSE.   ! multilayer sea ice scheme
      ELSE IF ((n19 == 1).AND.(.NOT.l_use_tpps_ozone)) THEN
        lmask=.FALSE.   ! tropopause ozone scheme
      ELSE IF ((n20 == 1).AND.(can_model /= 4)                  &
                         .AND.(nsmax == 0)) THEN
        lmask=.FALSE.   ! snow canopy scheme
      ELSE IF (n25 == 3) THEN
        lmask=.FALSE.   ! Direct PAR prognostic not currently used by any scheme
      ELSE IF ((n28 == 1).AND.(.NOT.l_use_cariolle)) THEN
        lmask=.FALSE.   ! cariolle ozone scheme
      ELSE IF ((n29 == 1).AND.(.NOT.l_endgame)) THEN
        lmask=.FALSE.   ! Disable ENDGame prognostics in ND run
      END IF

    ELSE IF (n30 == 1) THEN
! Can be used for a new set of option codes with n30=1
! when those with n30=0 are fully used up.
    END IF ! n30
  END IF ! SUM_IOPN

END IF ! isec

! Print out details of any tracer boundary conditions.
IF ( (isec == 36) .AND. printstatus >= prstatus_diag) THEN
  WRITE(6,*) 'TSTMSK:',isec,n2n1,tca_lbc(n2n1),lmask
END IF
IF ( (isec == 37)  .AND. printstatus >= prstatus_diag) THEN
  WRITE(6,*) 'TSTMSK:',isec,n2n1,tc_lbc_ukca(n2n1),lmask
END IF
! End of atmos primary block

! Atmos diagnostics
IF (isec == 1) THEN
! Shortwave radiation
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1)    ,10)
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n4 = MOD((iopn(1)/1000),10)
    n5 = MOD((iopn(1)/10000),10)
    n6 = MOD( iopn(2),10)
    IF ((n1 == 1).AND.(h_global(a_im) /= 'Y')) THEN
      lmask=.FALSE.
    ELSE IF ((n4  ==  1).AND.(.NOT. l_use_sulpc_indirect_sw)) THEN
      lmask = .FALSE.
    ELSE IF ((n5  ==  1).AND.                                     &
            (.NOT.(l_use_seasalt_direct                           &
              .OR. l_use_seasalt_indirect))) THEN
      lmask = .FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 2) THEN
! Longwave radiation
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1)   ,10)
    IF ((n1 == 1).AND.                                            &
        (i_ozone_int /= io3_trop_map_masscon)) THEN
      lmask=.FALSE.
    END IF
    n6 = MOD( iopn(2),10)
    IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 3) THEN
! Boundary layer
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1)    ,10)
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n4 = MOD((iopn(1)/1000),10)
    n6 = MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    IF ((n1 == 1).AND.(orogr == 'N')) THEN
      lmask=.FALSE.
    END IF
    IF ((n2 == 1).AND.(.NOT.l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2).AND.(.NOT.l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3).AND.(.NOT.l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4).AND.(.NOT.l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5).AND.(.NOT.l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6).AND.(.NOT.l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7).AND.(.NOT.l_nitrate)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 8).AND.(.NOT.l_dust).AND.(.NOT.l_dust_diag)) THEN
      lmask=.FALSE.
    END IF
    IF ((n3 == 1).AND.(.NOT.l_co2_interactive)) THEN
      lmask=.FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    IF ((n4 == 1) .AND. (nice == 1))THEN
      lmask = .FALSE.
    ELSEIF ((n4 == 2) .AND. (nice_use == 1)) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust)) THEN
      lmask = .FALSE.
    ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT.l_dust)) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF(isec == 4) THEN
! Large-scale precipitation
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n6 = MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    IF ((n2 == 1).AND.(.NOT.l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2).AND.(.NOT.l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3).AND.(.NOT.l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4).AND.(.NOT.l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5).AND.(.NOT.l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6).AND.(.NOT.l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7).AND.(.NOT.l_nitrate)) THEN
      lmask=.FALSE.
    END IF
    IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust)) THEN
      lmask = .FALSE.
    ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT.l_dust)) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 5) THEN
! Convection
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n6 = MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    IF ((n2 == 1).AND.(.NOT.l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ( (n3 == 1).AND.(l_3d_cca) ) THEN
      lmask=.FALSE.
    ELSE IF ( (n3 == 2).AND.(.NOT.l_3d_cca) ) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2).AND.(.NOT.l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3).AND.(.NOT.l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4).AND.(.NOT.l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5).AND.(.NOT.l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6).AND.(.NOT.l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7).AND.(.NOT.l_nitrate)) THEN
      lmask=.FALSE.
    END IF
    IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust)) THEN
      lmask = .FALSE.
    ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT.l_dust)) ) THEN
      lmask = .FALSE.
    END IF
  END IF
  n4 = MOD((iopn(1)/1000),10) 
  IF ( (n4 /= 1).AND.(.NOT. l_param_conv) ) THEN
    lmask=.FALSE.
  ENDIF

ELSE IF (isec == 6) THEN
! Gravity Wave Drag parametrization
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    IF ((n2 == 1).AND.(.NOT.l_gwd)) THEN
      lmask=.FALSE.
    ELSE IF ((n3 == 1).AND.(.NOT.l_use_ussp)) THEN
      lmask=.FALSE.
    END IF
  END IF

ELSE IF(isec == 8) THEN
! Hydrology
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1),10)
    n22 = MOD((iopn(5)/10),10)
    IF ((n1 == 1).AND.(.NOT.l_top)) THEN
      lmask=.FALSE.
    ELSE IF ((n22 == 1).AND.(.NOT.l_rivers)) THEN
          lmask=.FALSE.    !River routing switched off
    ELSE IF ((n22 == 2).AND.(.NOT.l_inland)) THEN
          lmask=.FALSE.    !Inland re-routing switched off
    END IF
  END IF

ELSE IF(isec == 9) THEN
! Cloud parametrization
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n6 = MOD( iopn(2),10)
    IF ((n2 == 1).AND.(.NOT.l_cld_area)) THEN
      lmask=.FALSE.
    ELSE IF ((n3 == 1).AND.(.NOT.l_rhcpt)) THEN
      lmask=.FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF(isec == 12) THEN
! Dynamics Advection
  IF (sum_iopn /= 0) THEN
    n6 = MOD( iopn(2),10)
    IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF(isec == 16) THEN
! Extra physics
  IF (sum_iopn /= 0) THEN
    n2n1 = MOD(iopn(1),100)
    n6 = MOD( iopn(2),10)
    IF ((n2n1 /= 0).AND.(.NOT.tracer_a(n2n1))) THEN
      lmask=.FALSE.
    END IF
    IF ( ( n6  ==  1 ) .AND. ( .NOT. l_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 17) THEN
! CLASSIC Aerosol section
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1),     10)
    n2 = MOD((iopn(1)/10),10)
    ! Sulphur cycle diagnostics
    IF ((n1 == 1).AND.(.NOT.l_sulpc_so2)) THEN
      lmask=.FALSE.
    ! Soot diagnostics
    ELSE IF ((n1 == 2).AND.(.NOT.l_soot)) THEN
      lmask=.FALSE.
    ! Biomass aerosol diagnostics
    ELSE IF ((n1 == 3).AND.(.NOT.l_biomass)) THEN
      lmask=.FALSE.
    ! Dust aerosol diagnostics
    ELSE IF ((n1 == 4).AND.(.NOT.l_dust)) THEN
      lmask=.FALSE.
    ! OCFF aerosol diagnostics
    ELSE IF ((n1 == 5).AND.(.NOT.l_ocff)) THEN
      lmask=.FALSE.
    ! SOA aerosol diagnostics
    ELSE IF ((n1 == 6).AND.(.NOT.l_use_biogenic)) THEN
      lmask=.FALSE.
    ! Sea-salt PM diagnostics
    ELSE IF ((n1 == 7).AND.(.NOT.l_use_seasalt_pm)) THEN
      lmask=.FALSE.
    ! Nitrate aerosol diagnostics
    ELSE IF ((n1 == 8).AND.(.NOT.l_nitrate)) THEN
      lmask=.FALSE.
    ! aerosol PM10/PM2.5 diagnostics - any of the above mean it is available
    ELSE IF ((n1 == 9).AND.(.NOT.                                   &
           (l_sulpc_so2      .OR. l_soot .OR. l_biomass      .OR.   &
            l_dust           .OR. l_ocff .OR. l_use_biogenic .OR.   &
            l_use_seasalt_pm .OR. l_nitrate) ) ) THEN
      lmask=.FALSE.
    END IF
    ! The following option codes are usually used 
    ! in conjunection with n1 ==1
    ! DMS diagnostics 
    IF ((n2 == 1).AND.(.NOT.l_sulpc_dms)) THEN
      lmask=.FALSE.
    ! Diagnostics which depend on oxidation of sulphur by ozone
    ELSE IF ((n2 == 2).AND.(.NOT.l_sulpc_ozone)) THEN
      lmask=.FALSE.
    END IF
  END IF

ELSE IF(isec == 18) THEN
! Data assimilation
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1),10)
    n2 = MOD((iopn(1)/10),10)
  END IF

ELSE IF ((isec >= 1).AND.(isec <= 20)) THEN  ! BUT NOT 1,3,18
  IF (sum_iopn /= 0) THEN
    WRITE(6,*)                                                    &
   'MESSAGE FROM ROUTINE TSTMSK: UNEXPECTED OPTION CODE',         &
    iopn(1)
    WRITE(6,*)'IN ATMOSPHERE SECTION ',isec
  END IF

ELSE IF ((isec >= 21).AND.(isec <= 24)) THEN
! Atmos climate mean diagnostics - redundant

ELSE IF (isec  ==  26) THEN
!  River Routing Diagnostics
  IF (sum_iopn /= 0) THEN
    n22 = MOD((iopn(5)/10),10)
    IF ((n22 == 1).AND.(.NOT.l_rivers)) THEN
      lmask=.FALSE.    !  River routing switched off
    ELSE IF ((n22 == 2).AND.(.NOT.l_inland)) THEN
      lmask=.FALSE.    !Inland re-routing switched off
    END IF
  END IF
END IF
! End of atmos diagnostic block

IF(lmask) THEN
  ladres=.TRUE.
  IF ((ispace == 3).OR.(ispace == 5)) THEN
    lmask=.FALSE.
  END IF
ELSE
  ladres=.FALSE.
END IF

 9999 CONTINUE

IF (lhook) CALL dr_hook('TSTMSK',zhook_out,zhook_handle)
RETURN

END SUBROUTINE tstmsk
