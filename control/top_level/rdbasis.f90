! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Read the basis file information
! Subroutine Interface:

SUBROUTINE rdbasis                                                      &
   (iu,cmessage,ErrorStatus)

  USE ereport_mod, ONLY: ereport
  USE check_iostat_mod
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE filenamelength_mod, ONLY :                                        &
     filenamelength
  USE stextend_mod, ONLY: NPOS_TS, NRECS_TS
  USE um_parparams
  USE lbc_mod
  USE missing_data_mod, ONLY:  rmdi, imdi
  USE umsections_mod
  USE version_mod, ONLY:                                                &
      nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,               &
      nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,             &
      nlevp, npslevp, npslistp, outfile_s, outfile_e,                   &
      tsdp

  USE Submodel_Mod
  IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code description:
! FORTRAN 90
! Written to UM programming standards UMDP3  vn 8.3

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

! Subroutine arguments

! Scalar arguments with intent(in):
  INTEGER      :: iu       ! Unit no. of stash basis file
  INTEGER      :: ie
! Scalar arguments with intent(out):
  CHARACTER(LEN=80) :: cmessage      ! Error return message

! Error status:
  INTEGER :: ErrorStatus ! Error return code
  INTEGER :: Readstatus  ! read  return code

! Local variables:
  LOGICAL :: model_lev  !TRUE for model levels
  LOGICAL :: lswap ! Indicates unfinished sort in the bubble sort block.
  LOGICAL :: lswap_diag  ! Current element requires a swap in the bubble sort.
  INTEGER :: i,j,k,l
  INTEGER :: idum, imod
  INTEGER :: IOSTAT
  INTEGER :: ntsrecs    !Counter for time series records
  INTEGER :: list_req_t1(3) ! Temporary list for request's STASHMaster data.
  INTEGER :: list_req_t2(3) ! Temporary list for request's STASHMaster data.
  INTEGER :: list_n_pro_t1(3) ! Temporary list for request's profile indices.
  INTEGER :: list_n_pro_t2(3) ! Temporary list for request's profile indices.
  CHARACTER(LEN=8) :: list_pro_t1(3) ! Temporary list for request's profiles.
  CHARACTER(LEN=8) :: list_pro_t2(3) ! Temporary list for request's profiles.
  CHARACTER (LEN=filenamelength) :: filename = "dummy file"


! Namelist STREQ: STASH requests
  INTEGER     :: isec, item       ! stash section and item
  CHARACTER(LEN=8) :: dom_name         ! domain profile name
  CHARACTER(LEN=8) :: use_name         ! usage profile name
  CHARACTER(LEN=8) :: tim_name         ! time profile name
  CHARACTER(LEN=20):: package          ! name of package
  
  NAMELIST/STREQ/ isec,item,dom_name,tim_name,use_name,package


! Namelist TIME: Time profiles
  CHARACTER(LEN=2) :: unt1,unt2,unt3
  INTEGER     :: ityp,isam,intv,iopt,itim
  INTEGER     :: istr,iend,ifre,ioff,itimes,iser(ntimep)
  INTEGER     :: isdt(6), iedt(6)
  
  NAMELIST/time/ityp,isam,intv,unt1  ,unt2,unt3,iopt                    &
               ,istr,iend,ifre,ioff,itimes,iser,tim_name                &
               ,isdt,iedt


! Namelist DOMAIN: Domain profiles
  INTEGER     ::  iopl                !Level type code
  INTEGER     ::  ilevs               !Flag for range/selected model levels
  INTEGER     ::  levb,levt           !Bottom/top levels for range
  INTEGER     ::  plt                 !Pseudo level type code
  INTEGER     ::  iopa                !Horizontal domain type code
  INTEGER     ::  inth,isth,iest,iwst !Horiz domain limits (IOPA=9,10)
  INTEGER     ::  imsk                !Grid type code
  INTEGER     ::  imn                 !Spatial meaning code
  INTEGER     ::  iwt                 !Weighting code
  INTEGER     ::  levlst (nlevp)      !Levels lists array: integer
  REAL        ::  rlevlst(nlevp)      !real
  INTEGER     ::  pslist (npslevp)    !pseudo
  CHARACTER(LEN=1) ::  ts                        !Flag for time series domain
  INTEGER     ::  tsnum                     !No. of time ser doms in prof
  INTEGER     ::  tblim (tsdp),ttlim (tsdp) !TS dom limits (top/bottom)
  REAL        ::  tblimr(tsdp),ttlimr(tsdp) !ditto for real levels
  INTEGER     ::  tnlim (tsdp),tslim (tsdp) !TS dom limits (N/S)
  INTEGER     ::  twlim (tsdp),telim (tsdp) !TS dom limits (E/W)

  NAMELIST/domain/iopl ,ilevs ,levb  ,levt  ,plt    ,iopa ,imsk ,       &
                  imn  ,iwt   ,ts    ,levlst,rlevlst,dom_name ,         &
                  inth ,isth  ,iest  ,iwst  ,pslist ,                   &
                  tsnum,tblim ,ttlim ,tnlim ,tslim  ,telim,twlim,       &
                  tblimr,ttlimr

! Namelist USE: Usage profiles
  INTEGER :: locn,iunt
  NAMELIST/USE/locn,iunt,use_name

! Function and subroutine calls:
  LOGICAL  disct_lev

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!- End of Header ------------------------------------------------------

! Initialisation
  IF (lhook) CALL dr_hook('RDBASIS',zhook_in,zhook_handle)
  ndiag   =0
  ntprof  =0
  ndprof  =0
  nuprof  =0
  ntsrecs =0

  DO i = 1,ndiagpm
    modl_b(i)=0
    isec_b(i)=0
    item_b(i)=0
    itim_b(i)=0
    idom_b(i)=0
    iuse_b(i)=0
  END DO
  
  itimes   =0
  DO i=1,nproftp
    timpro(i)='        '
    ityp_t(i)=0
    intv_t(i)=0
    unt1_t(i)='  '
    isam_t(i)=0
    unt2_t(i)='  '
    iopt_t(i)=0
    istr_t(i)=0
    iend_t(i)=0
    DO j=1,6
      isdt_t(j, i)=0
      iedt_t(j, i)=0
    END DO
    ifre_t(i)=0
    ioff_t(i)=0
    unt3_t(i)='  '
    itim_t(i)=0
    modl_t(i)=1          ! hard wired to set model to atmosphere 
    DO j = 1,ntimep
      iser_t(j,i)=0
    END DO
  END DO

  DO i=1,nprofdp
    dompro  (i)='        '
    iopl_d  (i)=0
    levb_d  (i)=0
    levt_d  (i)=0
    plt_d   (i)=0
    iopa_d  (i)=0
    inth_d  (i)=0
    isth_d  (i)=0
    iest_d  (i)=0
    iwst_d  (i)=0
    imsk_d  (i)=0
    imn_d   (i)=0
    iwt_d   (i)=0
    ts_d    (i)=' '
    pllen_d (i)=0
    plpos_d (i)=0
    ilev_d  (i)=0
  END DO
  DO i = 1,nlevp
    DO j = 1,nprofdp
      levlst_d (i,j)=0
      rlevlst_d(i,j)=0.0
    END DO
  END DO
  DO i = 1,npslevp
    DO j = 1,npslistp
      pslist_d(i,j)=0
    END DO
  END DO

  DO i = 1,nprofup
    usepro(i)='        '
    locn_u(i)=0
    iunt_u(i)=0
  END DO

! Get name for stash control file
  CALL GET_FILE(iu,filename,filenamelength,ErrorStatus)

! Open stash control file
  OPEN(iu,FILE=filename,IOSTAT=ErrorStatus)
  IF (ErrorStatus >  0) THEN
     WRITE(6,'(a)')'RDBASIS : Failed in OPEN of Stash Control File'
     GOTO 9999
  ELSE IF (ErrorStatus <  0) THEN
     WRITE(6,'(a)')'RDBASIS : Warning message on OPEN of Stash Control File'
     WRITE(6,'(a,i5)')'IOSTAT= ',ErrorStatus
  END IF

! ---------------------------------------------

  ReadStatus = 0

! read in the domain profiles
  DO WHILE (ReadStatus == 0)

! Initialise
    dom_name ='        '
    iopl =0
    levb =0
    levt =0
    ilevs=0
    levlst (1)= imdi
    rlevlst (1)= rmdi
    DO j = 2,nlevp
      levlst (j)= 0
      rlevlst(j)= rmdi
    END DO
    DO j = 1,npslevp
      pslist(j)=0
    END DO
    DO j = 1,tsdp
      tblim (j)=0
      ttlim (j)=0
      tblimr(j)=0.0
      ttlimr(j)=0.0
      tnlim (j)=0
      tslim (j)=0
      telim (j)=0
      twlim (j)=0
    END DO

    READ (UNIT=iu, NML=domain, IOSTAT=ReadStatus)
! only call check_iostat if ReadStatus is > 0
! we do not want EOF to report a warning.

    IF (Readstatus > 0) THEN
      CALL check_iostat(ReadStatus, "namelist domain" )
      
    ELSE IF (Readstatus == 0) THEN
    
      ndprof = ndprof + 1  
      
      IF (ndprof >  nprofdp) THEN
        WRITE(6,'(a)') 'ERROR IN STASHC:'
        WRITE(6,'(a,i5)') 'NUMBER OF DOMAIN PROFILES EXCEEDS LIMIT OF ',&
              nprofdp
        WRITE(6,'(a)') 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF

! Check for errors in levels lists
! DEPENDS ON: disct_lev
      model_lev=disct_lev(iopl,ErrorStatus,cmessage)
      
      IF (model_lev) THEN
! Model levels
        IF (ilevs == 2) THEN
          IF ( levlst(1) == imdi ) THEN
            WRITE(6,'(a,i5,a)')                                         &
                 'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '        &
                 ,ndprof,' HAS NO ENTRIES'
            cmessage='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
            ErrorStatus=1
            GO TO 9999
          END IF
        END IF
      ELSE IF (iopl /= 5) THEN
        IF (rlevlst(1) == rmdi) THEN
          WRITE(6,'(a,i5,a)')                                           &
               'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '          &
               ,ndprof,' HAS NO ENTRIES'
          cmessage='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
          ErrorStatus=1
          GO TO 9999
        END IF
      END IF
      
! Profile name
        dompro(ndprof)=dom_name
! Store level type code in IOPL_D
        iopl_d(ndprof)=iopl
! DEPENDS ON: disct_lev
        model_lev=disct_lev(iopl,ErrorStatus,cmessage)
        
      IF (model_lev) THEN
! Integer levels
        ilev_d(ndprof)=ilevs
        IF (ilevs == 1) THEN
!  Range of model levels
          levb_d(ndprof)=levb
          levt_d(ndprof)=levt
        END IF
        IF (ilevs == 2) THEN
! List of selected model levels
          levb_d(ndprof)=-1
          DO j=1,nlevp
            levlst_d(j,ndprof) = levlst(j)
            IF (levlst(j) >  0) THEN
! Store no. of levels in LEVT_D(ndprof)
              levt_d(ndprof)=levt_d(ndprof)+1
            END IF
          END DO
        END IF
      ELSE IF (iopl /= 5) THEN
! Real levels
        levb_d(ndprof)=-1
        DO j=1,nlevp
          rlevlst_d(j,ndprof) = rlevlst(j)
          IF (rlevlst(j) /=  rmdi) THEN
! Store no. of levels in LEVT_D(NDPROF)
            levt_d(ndprof)=levt_d(ndprof)+1
          END IF
        END DO
      END IF
! Store pseudo level type code in PLT_D
      plt_d (ndprof)=plt
      IF (plt >  0) THEN
! Domain profile 'NDPROF' has pseudo levels list
! Count total no. of pseudo levels lists
        npslists = npslists + 1
! Store list in column 'NPSLISTS' of PSLIST_D
        DO j=1,npslevp
          pslist_d(j,npslists) = pslist (j)
! PPLEN_D(NDPROF) stores length of ps.lev.list for domain 'ndprof'
          IF (pslist(j) >  0) THEN
            pllen_d (ndprof)        = pllen_d(ndprof) + 1
          END IF
        END DO
! PLPOS(NDPROF) stores the column no. in PSLIST_D for dom. prof. 'NDPROF'
        plpos_d(ndprof) = npslists
      END IF
! Store horizontal domain type in IOPA_D
      iopa_d(ndprof)=iopa
      IF (iopa == 9.OR.iopa == 10) THEN
! Specified area
        inth_d(ndprof)=inth
        isth_d(ndprof)=isth
        iest_d(ndprof)=iest
        iwst_d(ndprof)=iwst
      END IF
      imsk_d(ndprof)=imsk ! Gridpoint option
      imn_d (ndprof)=imn  ! Meaning option
      iwt_d (ndprof)=iwt  ! Weighting option
      ts_d  (ndprof)=ts   ! Time series domain
      IF (ts_d(ndprof)  ==  'Y') THEN
! This domain profile has a block of time series domains
! Store time series data for current profile in _TS arrays
        nseries    = nseries+1        ! Time series block number:
        npos_ts(ndprof) = nseries          !      used as a pointer
        nrecs_ts(nseries) = tsnum     ! No. of records in ts block
        nserrec_s  = nserrec_s+tsnum  ! Cumulative total ts records
        IF (nserrec_s <= ntimserp) THEN
          DO j = 1,tsnum
            IF (j <= tsdp) THEN
              ntsrecs = ntsrecs+1
! DEPENDS ON: disct_lev
              model_lev=disct_lev(iopl,ErrorStatus,cmessage)
              IF (model_lev) THEN
                blim_ts (ntsrecs)= tblim (j)
                tlim_ts (ntsrecs)= ttlim (j)
              ELSE IF (iopl /= 5) THEN
                blimr_ts(ntsrecs)= tblimr(j)
                tlimr_ts(ntsrecs)= ttlimr(j)
              END IF
                nlim_ts(ntsrecs) = tnlim(j)
                slim_ts(ntsrecs) = tslim(j)
                elim_ts(ntsrecs) = telim(j)
                wlim_ts(ntsrecs) = twlim(j)
            ELSE
              WRITE(6,'(a,a,i5,a,i5,a,a)')                              &
                   'MESSAGE FROM ROUTINE RDBASIS: ',                    &
                   'no. of time series in domain profile ',ndprof,      &
                   ' exceeds allowed limit of ',tsdp,' some will be',   &
                   ' ignored'
            END IF
          END DO
        ELSE
          WRITE(6,'(a,a,i5,a)')                                         &
               'TIMSER: total no. of time series requested exceeds ',   &
               'allowed limit of ',ntimserp,'; some will be ignored.'
        END IF
      ELSE
        npos_ts  (ndprof)=-1
      END IF
    END IF
  END DO
  
    nserblk_s = nseries

! Rewind stash control file ahead of next namelist searches.
  REWIND(iu)

  ReadStatus = 0

! Time profile namelists
  DO WHILE (ReadStatus == 0)

    tim_name='        '
    isam=0
    intv=0
    iopt=0
    istr=0
    iend=0
    ifre=0
    ioff=0
    UNT1=''
    UNT2=''
    UNT3=''
    DO j = 1,ntimep
      iser(j)=0
    END DO
    DO j = 1,6
      isdt(j) = 0
      iedt(j) = 0
    END DO

    READ (UNIT=iu, NML=time, IOSTAT=ReadStatus)
! only call check_iostat if ReadStatus is > 0
! we do not want EOF to report a warning.

    IF (Readstatus > 0) THEN
      CALL check_iostat(ReadStatus, "namelist time" )
      
    ELSE IF (readstatus == 0) THEN   
    
      ntprof = ntprof + 1

      IF (ntprof >  nproftp) THEN
        WRITE(6,'(a)') 'ERROR IN STASHC:'
        WRITE(6,'(a,i5)') 'NUMBER OF TIME PROFILES EXCEEDS LIMIT OF ',  &
               nproftp
        WRITE(6,'(a)') 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF

      timpro(ntprof)=tim_name
      ityp_t(ntprof)=ityp
      IF (ityp /= 1) THEN
! Diagnostic output is time-processed
        isam_t(ntprof)=isam  !Sampling frequency
        unt2_t(ntprof)=unt2
        intv_t(ntprof)=intv  !Processing interval
        unt1_t(ntprof)=unt1
      END IF
! Diag. output time option
      iopt_t(ntprof)=iopt
      unt3_t(ntprof)=unt3
      IF (iopt == 1) THEN
! Regular output time interval
        istr_t(ntprof)=istr
        iend_t(ntprof)=iend
        ifre_t(ntprof)=ifre
        ioff_t(ntprof)=ioff
      ELSE IF (iopt == 2) THEN
! Specified list of output times
! Length of times table
        itim_t(ntprof)=itimes
! Times table
        DO j = 1,itimes
          iser_t(j,ntprof)=iser(j)
        END DO
      ELSE IF (iopt == 3) THEN
        DO j = 1,6
          isdt_t(j,ntprof)=isdt(j)
          iedt_t(j,ntprof)=iedt(j)
        END DO
        ifre_t(ntprof)=ifre
      END IF
    END IF
  END DO


! Rewind stash control file
  REWIND(iu)

  ReadStatus = 0
  
!Usage profile namelists
  DO WHILE (ReadStatus == 0)
  
    use_name='        '
    locn=0
    iunt=0
      
    READ (UNIT=iu, NML=USE, IOSTAT=ReadStatus)  
! only call check_iostat if ReadStatus is > 0
! we do not want EOF to report a warning.

    IF (Readstatus > 0) THEN
      CALL check_iostat(ReadStatus, "namelist use" )
      
    ELSE IF ( readstatus == 0) THEN
    
      nuprof = nuprof + 1
    
      IF (nuprof >  nprofup) THEN
        WRITE(6,'(a)') 'ERROR IN STASHC:'
        WRITE(6,'(a,i5)') 'NUMBER OF USAGE PROFILES EXCEEDS LIMIT OF ',&
              nprofup
        WRITE(6,'(a)') 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF
   
      usepro(nuprof)=use_name
      locn_u(nuprof)=locn
      iunt_u(nuprof)=iunt
      
    END IF
  END DO 
  
! Rewind stash control file
  REWIND(iu)


!--------------------------------------------------

  ReadStatus = 0  
  
  imod=1   ! hard-wired to atmosphere as only sub model

! read in the stash requests
  DO WHILE (ReadStatus == 0)

    isec=0
    item=0
    
    READ (UNIT=iu, NML=streq, IOSTAT=ReadStatus)    
! only call check_iostat if ReadStatus is > 0
! we do not want EOF to report a warning.

    IF (Readstatus > 0) THEN
      CALL check_iostat(ReadStatus, "namelist streq" )
      
    ELSE IF (readstatus ==0 ) THEN
    
      ndiag = ndiag + 1
      
      IF (ndiag  >  nrecdp ) THEN
        WRITE(cmessage,'(a,a,i5,a)') 'NUMBER OF DIAGNOSTIC REQUESTS ',  &
                                     'EXCEEDS LIMIT OF ',               &
                                     nrecdp ,' SOME HAVE BEEN IGNORED'
        errorstatus=-15
        CALL ereport("rdbasis",errorstatus,cmessage)    
        EXIT
      END IF

      modl_b(ndiag)=imod
      isec_b(ndiag)=isec
      item_b(ndiag)=item
      
      DO j=1,ntprof
        IF (tim_name == timpro(j)) THEN
          itim_b(ndiag)= j
          EXIT
        END IF
      END DO 
      
! check time profile has been set and is one of the time profiles
! already read in from the STASHC file.
      
      IF (itim_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED time profile request ',   &
                               tim_name
        errorstatus=16
        CALL ereport("rdbasis",errorstatus,cmessage)                         
      END IF 
        
      
      DO j=1,ndprof
        IF (dom_name == dompro(j)) THEN
          idom_b(ndiag)= j
          EXIT
        END IF
      END DO 
      
! check domain profile has been set and is one of the domain profiles
! already read in from the STASHC file.
      
      IF (idom_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED domain profile request ', &
                               dom_name
        errorstatus=17
        CALL ereport("rdbasis",errorstatus,cmessage)                         
      END IF       
        
        
      DO j=1,nuprof 
        IF (use_name == usepro(j)) THEN
          iuse_b(ndiag)= j
          EXIT
        END IF
      END DO  
   
! check usage profile has been set and is one of the usage profiles
! already read in from the STASHC file.
      
      IF (iuse_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED usage profile request ',  &
                               use_name
        errorstatus=18
        CALL ereport("rdbasis",errorstatus,cmessage)                         
      END IF    
   
    END IF
   
   
  END DO

! Rewind stash control file
  REWIND(iu)


  ! Now sort requests by modl_b, isec_b, item_b, and profile names.
  ! Use a bubble sort.

  lswap = .TRUE.
  
  DO WHILE (lswap)
  
    ! Loop over STREQ records and reorder based on their data.
    ! Continue until no more swaps occur (lswap is .FALSE.).
    lswap = .FALSE.

    DO i = ndiag-1, 1, -1
      j = i + 1

      ! Construct temporary lists used for sorting and temporary storage.
      ! list_req_t1, list_req_t2 hold STASHMaster request data
      ! list_pro_t1, list_pro_t2 hold profile name data
      ! list_n_pro_t1, list_n_pro_t2 hold the profile index data

      list_req_t1 = (/ modl_b(i), isec_b(i), item_b(i) /)
      list_req_t2 = (/ modl_b(j), isec_b(j), item_b(j) /)
      list_n_pro_t1 = (/ idom_b(i), itim_b(i), iuse_b(i) /)
      list_n_pro_t2 = (/ idom_b(j), itim_b(j), iuse_b(j) /)
      list_pro_t1 = (/ dompro(idom_b(i)), &
                       timpro(itim_b(i)), &
                       usepro(iuse_b(i)) /)
      list_pro_t2 = (/ dompro(idom_b(j)), &
                       timpro(itim_b(j)), &
                       usepro(iuse_b(j)) /)

      lswap_diag = .FALSE.
      DO k = 1, 6
        IF (k < 4) THEN
          IF (list_req_t1(k) > list_req_t2(k)) THEN
              lswap_diag = .TRUE.
          ELSE IF (list_req_t1(k) < list_req_t2(k)) THEN
              EXIT
          END IF
        ELSE IF (k > 3) THEN
          IF (list_pro_t1(k-3) > list_pro_t2(k-3)) THEN
              lswap_diag = .TRUE.
          ELSE IF (list_pro_t1(k-3) < list_pro_t2(k-3)) THEN
              EXIT
          END IF
        END IF
      END DO ! End loop over sort properties
      
      IF (lswap_diag) THEN
        lswap=.TRUE.
        modl_b(i) = list_req_t2(1)
        modl_b(j) = list_req_t1(1)
        isec_b(i) = list_req_t2(2)
        isec_b(j) = list_req_t1(2)
        item_b(i) = list_req_t2(3)
        item_b(j) = list_req_t1(3)
        idom_b(i) = list_n_pro_t2(1)
        idom_b(j) = list_n_pro_t1(1)
        itim_b(i) = list_n_pro_t2(2)
        itim_b(j) = list_n_pro_t1(2)
        iuse_b(i) = list_n_pro_t2(3)
        iuse_b(j) = list_n_pro_t1(3)
      END IF
    END DO  ! End loop over modl_b, etc elements

  END DO ! End while loop


! Set defaults
  ocaaa = imdi
  lwbnd = imdi
  swbnd = imdi
  
  ewspacea = rmdi
  nsspacea = rmdi
  frstlona = rmdi
  frstlata = rmdi
  polelona = rmdi
  polelata = rmdi

  lats = rmdi
  lons = rmdi

  orogr = ' '
  swmcr = ' '
  lfloor = ' '
  meso  = ' '

  zonavozone = .FALSE.
  ZonAvTppsOzone = .FALSE.

  READ (UNIT=5, NML=stshcomp, IOSTAT=ErrorStatus)
  CALL check_iostat(ErrorStatus, "namelist STSHCOMP")
  ZonAvTppsOzone=ZonAvOzone

! Set the number of tracers being used by counting the non-zero elements
! of the arrays that list them
  TR_LBC_VARS=0
  TR_UKCA=0
  TR_LBC_UKCA=0
  DO I=1,A_MAX_TRVARS
    IF(TCA_LBC(I) /= 0) TR_LBC_VARS = TR_LBC_VARS + 1
  END DO
  DO I=1,A_MAX_UKCAVARS
    IF(TC_UKCA(I) /= 0) TR_UKCA = TR_UKCA +1
    IF(TC_LBC_UKCA(I) /= 0) TR_LBC_UKCA = TR_LBC_UKCA + 1
  END DO

! Set rim parameters to 1 unless this is a limited-area run
  IF(OCAAA == 1)THEN
    RIMWIDTHA=1
    NRIM_TIMESA=1
  ELSE IF(OCAAA == 2 .OR. OCAAA == 3 .OR. OCAAA == 4)THEN
    NRIM_TIMESA=MAX(NRIM_TIMESA,1)
  END IF

! Update RIM_LOOKUPSA with the number of tracer LBCs 
  RIM_LOOKUPSA = RIM_LOOKUPSA + TR_LBC_VARS + TR_LBC_UKCA
! Calculate total number of headers in the LBC file
  BOUND_LOOKUPSA=RIM_LOOKUPSA+(NRIM_TIMESA-1)*(RIM_LOOKUPSA-1)
! Number of atmosphere model interface lookups
  intf_lookupsa = intf_lookupsa + tr_vars + tr_ukca

  CLOSE(UNIT=iu,STATUS='KEEP',IOSTAT=IOSTAT)

9999 CONTINUE
  IF (lhook) CALL dr_hook('RDBASIS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE rdbasis
