! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:

      SUBROUTINE INPUTL(NRECS,                                          &
     &                   NLEVELS,ErrorStatus,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Decomp_DB

      USE ppxlook_mod, ONLY: ppxref_sections, ppxref_items
      USE cppxref_mod, ONLY:                                            &
          ppx_halo_type, ppx_grid_type, ppx_lev_flag,                   &
          ppx_pf_code, ppx_pl_code, ppx_pt_code, ppx_lv_code
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod
      USE stextend_mod, ONLY: LLISTTY, INDX_S, IN_S, LIST_S,  &
                              RLEVLST_S, LEVLST_S, LENPLST

      IMPLICIT NONE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!
! Global variables:
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

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER NRECS

!   Scalar arguments with intent(out):
      INTEGER NLEVELS   ! total no. of sets of stash levels

!   Scalar arguments with intent(out):

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      LOGICAL MODEL_LEV
      CHARACTER(LEN=80) CMESSAGE
      LOGICAL LADD
      LOGICAL LDUPLL
      INTEGER I,IL,ILIN
      INTEGER ISTART,IEND
      INTEGER MODL
      INTEGER ISEC
      INTEGER IITM
      INTEGER IP_IN
      INTEGER IX1,IX2
      INTEGER IY1,IY2
      INTEGER IZ_IN
      INTEGER LEN_IN
      INTEGER LEN_PRIMIN
      INTEGER NDUPLL
      INTEGER NLEVIN
      INTEGER LENO
      INTEGER IPF,IPL
! local versions of the global subdomain limits
      INTEGER local_IX1,local_IX2,local_IY1,local_IY2
      INTEGER                                                           &
     &        orig_decomp                                               &
                               ! MPP decomposition before start
     &       ,decomp_type                                               &
                               ! decomposition type
     &       ,sm_ident         ! submodel identifier

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      INTEGER  EXPPXI

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      EXTERNAL EXPPXI

!- End of Header ----------------------------------------------------

      IF (lhook) CALL dr_hook('INPUTL',zhook_in,zhook_handle)

      orig_decomp = current_decomp_type

      DO MODL=1,N_INTERNAL_MODEL_MAX
!
!    Ensure that domain decomposition is consistent with submodel
!
      sm_ident = SUBMODEL_PARTITION_INDEX(MODL)
      IF(sm_ident == atmos_sm) THEN
         decomp_type = decomp_standard_atmos
      ELSE                            ! No decomposition defined:
         decomp_type = orig_decomp    !  return to original
      ENDIF

      CALL CHANGE_DECOMPOSITION(decomp_type,ErrorStatus)

      IF(ErrorStatus >  0) THEN
         CMESSAGE='INPUTL: ERROR in changing MPP decomposition'
         write(6,*) CMESSAGE
         GOTO 9999
      ENDIF
      DO ISEC=0,PPXREF_SECTIONS
      DO IITM=1,PPXREF_ITEMS
        IF(INDX_S(2,MODL,ISEC,IITM) >= 1) THEN
! At least one stash rec
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL,ISEC,IITM,ppx_grid_type     ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL,ISEC,IITM,ppx_lv_code       ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IFLAG   = EXPPXI(MODL,ISEC,IITM,ppx_lev_flag      ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL,ISEC,IITM,ppx_pt_code       ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPFIRST = EXPPXI(MODL,ISEC,IITM,ppx_pf_code       ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPLAST  = EXPPXI(MODL,ISEC,IITM,ppx_pl_code       ,             &
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        HALO_TYPE  = EXPPXI(MODL,ISEC,IITM,ppx_halo_type,               &
     &                                         ErrorStatus,CMESSAGE)

          ISTART=       INDX_S(1,MODL,ISEC,IITM)   ! Pos of 1st rec
          IEND  =ISTART+INDX_S(2,MODL,ISEC,IITM)-1 ! Pos of last rec
! Diagnostics with input on levels list (IFLAG=1),
!  rather than on all possible levels
          IF((IFLAG  == 1   ).AND.                                      &
                                                    !Input on lev list
     &       (ISTART == IEND).AND.                                      &
                                                    !Only 1 stash rec
     &       (LIST_S(st_output_bottom,ISTART) <  0))                    &
                                                    !Output on lev list
     &        THEN
! Only one stash record for this m,s,i - output levels list is
!  the same as the input levels list
            LIST_S(st_input_bottom ,ISTART)=                            &
     &      LIST_S(st_output_bottom,ISTART)
          ELSE IF (IFLAG == 1.AND.ILEV /= 5) THEN
! Input on levels list & more than one stash request -
!  construct input levels list
            NLEVELS=NLEVELS+1
            IF (NLEVELS >  NLEVLSTSP) THEN
              WRITE(6,*) 'ERROR IN ROUTINE INPUTL:'
              WRITE(6,*) 'TOO MANY STASH LEVELS LISTS REQUESTED ',      &
     &                   'ARRAYS WILL BE OVERWRITTEN'
              WRITE(6,*) 'REDUCE NUMBER OF LEVELS LISTS'
              ErrorStatus=1
              GO TO 9999
            END IF
! Construct input levels list: this is the combined list of all
!  the output levels for all the stash requests for this m,s,i
            NLEVIN=1
! Set levels list type - real or integer
! DEPENDS ON: disct_lev
            MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
            IF (.NOT.MODEL_LEV) THEN
! Non-model levels - real
              LLISTTY(NLEVELS)='R'
! Model levels - integer
            ELSE
              LLISTTY(NLEVELS)='I'
            END IF
! Loop over stash recs for this m,s,i
            DO I=ISTART,IEND
! Pointer for input level list
              LIST_S(st_input_bottom,I)=-NLEVELS
              IF (LIST_S(st_output_bottom,I) <  0) THEN
! There is an output levels list:
!  For each of the levels in the output levels lists for the stash
!   record I, find out whether this level is already present in the
!   input levels list NLEVELS constructed so far.
!   If it is, set LADD=F. Otherwise, LADD=T.
! Loop over output levels and check each one
                DO IL=2,LEVLST_S(1,-LIST_S(st_output_bottom,I))+1
                  LADD=.TRUE.
                  IF(NLEVIN >  1) THEN
                    DO ILIN=2,NLEVIN
                      IF(LIST_S(st_output_top,I) /= 1) THEN
! Non-model levels: real
                        IF(RLEVLST_S(IL,-LIST_S(st_output_bottom,I))    &
     &                   ==                                             &
     &                     RLEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                      ELSE
! Model levels: integer
                        IF( LEVLST_S(IL,-LIST_S(st_output_bottom,I))    &
     &                   ==                                             &
     &                      LEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                      END IF
                    END DO
                  END IF

! If LADD=T, add level 'IL' from stash record 'I' output levels list
!  to input levels list NLEVELS
                  IF (LADD) THEN
                    NLEVIN=NLEVIN+1
                    IF(LIST_S(st_output_top,I) /= 1) THEN
                      RLEVLST_S(NLEVIN,NLEVELS)=                        &
     &                RLEVLST_S(IL,-LIST_S(st_output_bottom,I))
                    ELSE
                      LEVLST_S(NLEVIN,NLEVELS)=                         &
     &                LEVLST_S(IL,-LIST_S(st_output_bottom,I))
                    END IF
                  END IF
                END DO     ! Loop over levels

              ELSE
! Contiguous range of model levels for output, rather than list
!  Compare output levels range for stash record I with input levs
!   range NLEVELS. Any of the output levels not already present
!   in the input range is added to the input list.
                DO IL=LIST_S(st_output_bottom,I),                       &
     &                LIST_S(st_output_top   ,I)
                  LADD=.TRUE.
                  DO ILIN=2,NLEVIN
                    IF(IL == LEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                  END DO
                  IF(LADD) THEN
                    NLEVIN=NLEVIN+1
                    LEVLST_S(NLEVIN,NLEVELS)=IL
                  END IF
                END DO
              END IF   !  Levels list/range
            END DO     !  Loop over stash recs

! Record no. of levels in input list just constructed
            LEVLST_S(1,NLEVELS)=NLEVIN-1

            IF(NLEVIN-1 == 0) THEN
              WRITE(6,*) 'ORDINARY LEVEL'
              WRITE(6,*) 'ISEC=',ISEC
              WRITE(6,*) 'IITM=',IITM
              WRITE(6,*) 'NLEVELS=',NLEVELS
              DO I=ISTART,IEND
                WRITE(6,*) 'I=',I
                WRITE(6,*) 'LIST_S(st_output_bottom)=',                 &
     &                      LIST_S(st_output_bottom,I)
                IF (LIST_S(st_output_bottom,I) <  0) THEN
                  DO IL=2,LEVLST_S(1,-LIST_S(st_output_bottom,I))+1
                    WRITE(6,*) 'IL=',IL
                    WRITE(6,*)                                          &
     &              'LEVLST',LEVLST_S(IL,-LIST_S(st_output_bottom,I))
                  END DO
                ELSE
                  WRITE(6,*)                                            &
     &            'LIST_S(st_output_top=',LIST_S(st_output_top,I)
                END IF
              END DO
            END IF
! Sort levels list
! DEPENDS ON: levsrt
            CALL LEVSRT(LLISTTY(  NLEVELS), LEVLST_S(1,NLEVELS),        &
     &                 LEVLST_S(2,NLEVELS),RLEVLST_S(2,NLEVELS))

! Determine whether this levels list is a duplicate of another list
! DEPENDS ON: duplevl
            CALL DUPLEVL(NLEVELS,LDUPLL,NDUPLL)
            IF (LDUPLL) THEN
! Duplicate list at NDUPLL - reset pointer and reduce NLEVELS by 1
              NLEVELS=NLEVELS-1
              DO I=ISTART,IEND
                LIST_S(st_input_bottom,I)=-NDUPLL
              END DO
            END IF
          END IF   !Levels lists

! Pseudo levels lists
          IF((IFLAG  == 1   ).AND.                                      &
     &      ((ISTART == IEND).OR.(IPSEUDO == 0)) ) THEN
! Either no pseudo levels or only one request:
! Input pseudo levels list equals output list
            LIST_S(st_pseudo_in,ISTART)=LIST_S(st_pseudo_out,ISTART)
          ELSE IF (IFLAG == 1) THEN
! Input pseudo levels list with more than one request
            NPSLISTS=NPSLISTS+1
            IF(NPSLISTS >  NPSLISTP) THEN
              WRITE(6,*) 'ERROR IN ROUTINE INPUTL:'
              WRITE(6,*)                                                &
     &       'TOO MANY STASH PSEUDO LEVELS LISTS REQUESTED ',           &
     &       'ARRAYS WILL BE OVERWRITTEN'
              WRITE(6,*) 'REDUCE NUMBER OF PSEUDO LISTS'
              ErrorStatus=1
              GO TO 9999
            END IF
! Construct input pseudo list: combined list of all output
!  pseudo levels for all stash requests for this m,s,i
            NLEVIN=0
            DO I=ISTART,IEND
              LIST_S(st_pseudo_in,I)=NPSLISTS
              DO IL=1,LENPLST(LIST_S(st_pseudo_out,I))
                LADD=.TRUE.
                IF(NLEVIN >  0) THEN
                  DO ILIN=1,NLEVIN
                    IF( PSLIST_D(IL,LIST_S(st_pseudo_out,I)) ==         &
     &                  PSLIST_D(ILIN,NPSLISTS)) LADD=.FALSE.
                  END DO
                END IF
                IF(LADD) THEN
                  NLEVIN=NLEVIN+1
                  PSLIST_D(NLEVIN,NPSLISTS)=                            &
     &            PSLIST_D(IL,LIST_S(st_pseudo_out,I))
                END IF
              END DO
            END DO
            LENPLST(NPSLISTS)=NLEVIN

            IF(NLEVIN == 0) THEN
              WRITE(6,*) 'PSEUDO LEVEL'
              WRITE(6,*) 'ISEC=',ISEC
              WRITE(6,*) 'IITM=',IITM
              WRITE(6,*) 'NPSLISTS=',NPSLISTS
              DO I=ISTART,IEND
                WRITE(6,*) 'I=',I
                WRITE(6,*) 'LENPLST=',LENPLST(LIST_S(st_pseudo_out,I))
                DO IL=1,LENPLST(LIST_S(st_pseudo_out,I))
                  WRITE(6,*) 'IL=',IL
                  WRITE(6,*)                                            &
     &            'PSLIST_D',PSLIST_D(IL,LIST_S(st_pseudo_out,I))
                END DO
              END DO
            END IF
! Sort input pseudo levels list.  The REAL argument is really just a dummy
! since we are processing integer levels here.
! DEPENDS ON: levsrt
            CALL LEVSRT('I',LENPLST(NPSLISTS),PSLIST_D(1,NPSLISTS),     &
     &                  REAL(PSLIST_D(1:LENPLST(NPSLISTS),NPSLISTS)))
! Find out if duplicate
! DEPENDS ON: duppsll
            CALL DUPPSLL(LDUPLL,NDUPLL)
            IF(LDUPLL) THEN
! Duplicate pseudo list at NDUPLL
              NPSLISTS=NPSLISTS-1
              DO I=ISTART,IEND
                LIST_S(st_pseudo_in,I)=NDUPLL
              END DO
            END IF
          ELSE IF (IFLAG == 0.AND.IPSEUDO /= 0) THEN
! Input pseudo levels list contains all possible pseudo levels for
!  this diagnostic
            NPSLISTS=NPSLISTS+1
            DO I=ISTART,IEND
              LIST_S(st_pseudo_in,I)=NPSLISTS
! Decode first & last pseudo level codes from stash master
! DEPENDS ON: pslevcod
              CALL PSLEVCOD(IPFIRST,IPF,'F',ErrorStatus,CMESSAGE)
! DEPENDS ON: pslevcod
              CALL PSLEVCOD(IPLAST ,IPL,'L',ErrorStatus,CMESSAGE)
! Construct list
              DO NLEVIN = IPF,IPL
                PSLIST_D(NLEVIN,NPSLISTS)=NLEVIN
              END DO
            END DO
            LENPLST(NPSLISTS)=IPL-IPF+1
          END IF   ! Pseudo levels

! Calculate horizontal factor for input length
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,IY1,IY2,IX1,IX2)

! Convert from global to local subdomain limits
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE., .TRUE.,                 &
     &                                  IGP,halo_type,mype,             &
     &                                  IY1,IX2,IY2,IX1,                &
     &                                  local_IY1,local_IX2,            &
     &                                  local_IY2,local_IX1)
        IX1=local_IX1
        IX2=local_IX2
        IY1=local_IY1
        IY2=local_IY2
! All sub-model grids: atmos/ocean/wave now ordered S->N
           LEN_IN=(IX2-IX1+1)*(IY2-IY1+1)

! Calculate vertical levels factor for input length
          IF(ILEV /= 5) THEN
! More than one level
            IF(LIST_S(st_input_bottom,ISTART) <  0) THEN
! Level list
              IZ_IN=LEVLST_S(1,-LIST_S(st_input_bottom,ISTART))
            ELSE
! Range of model levs
              IZ_IN=LIST_S(st_input_top   ,ISTART)-                     &
     &              LIST_S(st_input_bottom,ISTART)+1
            END IF
          ELSE
! Single level input
            IZ_IN=1
          END IF

! Calculate pseudo levels factor for input length
          IF(IPSEUDO /= 0) THEN
            IP_IN=LENPLST(LIST_S(st_pseudo_in,ISTART))
          ELSE
            IP_IN=1
          END IF

! Calculate input length for this diag. and store in LIST_S
! Input_code <  0 means that a diag already processed into D1 is being
!   reprocessed, so input length of child diag equals output length of
!   parent.
! Otherwise, the input len is given by the product of the appropriate
!   x-,y-,z-, and p-dimensions.
          DO I=ISTART,IEND
            IF(LIST_S(st_input_code  ,I) >= 0) THEN
               LIST_S(st_input_length,I)=LEN_IN*IZ_IN*IP_IN
            ELSE
               LIST_S(st_input_length ,I)=                              &
     &         LIST_S(st_output_length,-LIST_S(st_input_code,I))
            END IF
 ! Store model no. in last element of LIST_S - for ADDRES
            LIST_S(NELEMP+1,I)=MODL
          END DO

! Recalculate input length for non-primary (length unchanged for
! most cases) and store in IN_S array.
          IF (ISEC /= 0) THEN
            IF ((IGP /= 31).AND.(IGP /= 32))THEN
! DEPENDS ON: addrln
              CALL ADDRLN(IGP,halo_type,LEN_PRIMIN,local_data)
              IN_S(2,MODL,ISEC,IITM)=LEN_PRIMIN*IZ_IN*IP_IN
            END IF
          END IF

        END IF ! At least one stash record for m,s,i

      END DO   ! Items
      END DO   ! Sections
      END DO   ! Models
!
      CALL CHANGE_DECOMPOSITION(orig_decomp,ErrorStatus)

      IF(ErrorStatus >  0) THEN
         CMESSAGE='INPUTL: ERROR in original MPP decomposition'
         write(6,*) CMESSAGE
         GOTO 9999
      ENDIF

 9999 CONTINUE

      IF (lhook) CALL dr_hook('INPUTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INPUTL
