! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Set the STASH addresses for D1
! Subroutine Interface:
      SUBROUTINE ADDRES(                                                &
                        NRECS,                                          &
                        ErrorStatus,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Decomp_DB
      USE ppxlook_mod, ONLY: ppxref_items, ppxref_sections, ppxptr
      USE cppxref_mod, ONLY:                                            &
          ppx_version_mask, ppx_grid_type, ppx_halo_type, ppx_lev_flag, &
          ppx_lb_code, ppx_lt_code, ppx_lv_code,                        &
          ppx_pf_code, ppx_pl_code, ppx_pt_code,                        &
          ppx_ptr_code, ppx_opt_code, ppx_space_code
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod
      USE stextend_mod, ONLY: INDX_S, IN_S, D1_PADDR, N_OBJ_D1,  &
                 MAX_D1_LEN, LIST_S, ITIM_S, LEVLST_S, LENPLST,  &
                     d1_type, d1_im, d1_extra_info, diag, seco

      IMPLICIT NONE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
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
      CHARACTER(LEN=80) CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      INTEGER TOTIMP
      INTEGER Im_ident  !Internal model identifier (absolute - CSMID)
      INTEGER Im_index  !Internal model index (expt. dependent)
      INTEGER Sm_ident  !Submodel identifier (absolute)
      INTEGER ISEC
      INTEGER ISEC_loop
      INTEGER IITM
      INTEGER RLEVS
      INTEGER RADDRESS
      INTEGER PIrow
      INTEGER I,J
      INTEGER IFIRST
      INTEGER IFREQ
      INTEGER IHOURS
      INTEGER ILAST
      INTEGER IREC
      INTEGER IH,IL,IP,IT
      INTEGER LWORK_S(N_SUBMODEL_PARTITION_MAX)
      INTEGER ICODE ! return from CHANGE_DECOMPOSITION

! Local arrays:
!    Submodel definitions array: stores list of Im_index's
!     for each submodel partition
      INTEGER                                                           &
       SM_def(N_SUBMODEL_PARTITION_MAX,N_INTERNAL_MODEL_MAX)

! Function and subroutine calls:
      INTEGER  EXPPXI

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      EXTERNAL EXPPXI

!- End of Header ----------------------------------------------------


! 1.  Set STASHIN addresses and input lengths for primary fields

!   The address loop for primary fields is performed for each
!   internal model in turn. Hence, each internal model's primary
!   data occupies a contiguous block in D1. The order of these blocks
!   is the same as the order of the internal models given in the
!   array INTERNAL_MODEL_LIST.
!   User-defined prognostics are included in this primary addressing
!   routine, since they are incorporated into the ppxref lookup
!   arrays PPXI, PPXC in routine GETPPX.

!   Initialisation
      IF (lhook) CALL dr_hook('ADDRES',zhook_in,zhook_handle)
      N_OBJ_D1_MAX=0
      DO I = 1,N_SUBMODEL_PARTITION_MAX
        N_OBJ_D1(I)=0
      END DO

      DO I = 1,N_SUBMODEL_PARTITION_MAX
      DO J = 1,N_INTERNAL_MODEL_MAX
        SM_def(I,J) = 0
      END DO
      END DO

!   Obtain submodel definitions and store in SMdef array
      DO Im_index = 1,N_INTERNAL_MODEL
!   Submodel ident.
         Sm_ident =   SUBMODEL_FOR_IM(Im_index)
!   Internal model index
         SM_def(Sm_ident,Im_index) = Im_index
      END DO

!   Primary address loop

!     Loop over submodel partitions
      DO Sm_ident = 1,N_SUBMODEL_PARTITION_MAX

!       Initialise LEXTRA
        LEXTRA(Sm_ident)=0

!     Initialise address for reconfiguration
        RADDRESS = 1

!      Loop over internal models for each SM partition
        DO Im_index = 1,N_INTERNAL_MODEL

!       Test whether current SM contains this IM
          IF (SM_def(Sm_ident,Im_index) >  0) THEN

!        Obtain internal model identifier
            Im_ident   = INTERNAL_MODEL_LIST(Im_index)

! Set the correct decomposition in PARVARS

            ICODE=0

            IF (Im_ident  ==  A_IM) THEN
              IF (current_decomp_type  /=  decomp_standard_atmos) &
              CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

            ELSE  ! unsupported decomposition type
              WRITE(6,'(A,A)') 'ADDRES1 : Error - Only atmosphere ', &
                   'submodel is currently supported for UM code.'
              ErrorStatus=-1
              CMESSAGE='Unsupported submodel for UM code'
              GOTO 9999
            END IF

            IF (ICODE  /=  0) THEN
              WRITE(6,'(A,A)') 'ADDRES1 : Error - Could not set ', &
                'decomposition for selected submodel.'
              ErrorStatus=-2
              CMESSAGE='Unsupported decomposition selected for UM code'
              GOTO 9999
            END IF

!        Initialise primary data lengths for TYPSIZE
            IF (Im_ident == A_IM) A_PROG_LEN=0
            IF (Im_ident == A_IM) A_PROG_LOOKUP=0
            PIrow  = 0
            DO ISEC_loop = 1,7   ! Currently there are seven sections
                                 ! that contain "primary" type fields
            IF (ISEC_loop == 1) ISEC=0  ! section zero primary
            IF (ISEC_loop == 2) ISEC=31 ! LBC Input (treated as primary)
            IF (ISEC_loop == 3) ISEC=32 ! LBC Output(treated as primary)
            IF (ISEC_loop == 4) ISEC=33 ! Tracers (treated as primary)
            IF (ISEC_loop == 5) ISEC=34 ! UKCA tracers (like primary)
            IF (ISEC_loop == 6) ISEC=36 ! Free Tracer lbcs
            IF (ISEC_loop == 7) ISEC=37 ! UKCA Tracer lbcs
!       Loop over section zero items
              DO IITM   = 1,PPXREF_ITEMS
!   Check whether there is a primary field corresponding
!         to this item number
                IF (PPXPTR(Im_index,ISEC,IITM) /= 0) THEN
! DEPENDS ON: exppxi
                  VMSK    = EXPPXI(Im_ident,ISEC,IITM,ppx_version_mask, &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code,   &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  IGP     = EXPPXI(Im_ident,ISEC,IITM,ppx_grid_type,    &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  ILEV    = EXPPXI(Im_ident,ISEC,IITM,ppx_lv_code,      &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  IBOT    = EXPPXI(Im_ident,ISEC,IITM,ppx_lb_code,      &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  ITOP    = EXPPXI(Im_ident,ISEC,IITM,ppx_lt_code,      &
                              ErrorStatus,CMESSAGE)
                  DO I=1,6
! DEPENDS ON: exppxi
                    IOPN(I)=EXPPXI(Im_ident,ISEC,IITM,ppx_opt_code+I-1, &
                              ErrorStatus,CMESSAGE)
                  END DO
! DEPENDS ON: exppxi
                  IFLAG   = EXPPXI(Im_ident,ISEC,IITM,ppx_lev_flag,     &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  IPSEUDO = EXPPXI(Im_ident,ISEC,IITM,ppx_pt_code,      &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  IPFIRST = EXPPXI(Im_ident,ISEC,IITM,ppx_pf_code,      &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  IPLAST  = EXPPXI(Im_ident,ISEC,IITM,ppx_pl_code,      &
                              ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
                  halo_type = EXPPXI(Im_ident,ISEC,IITM,ppx_halo_type,  &
                              ErrorStatus,CMESSAGE)
                  IF((ISPACE == 2).OR.(ISPACE == 3).OR.(ISPACE == 9)    &
                 .OR.(ISPACE == 4).OR.(ISPACE == 5).OR.(ISPACE == 10)   &
                 .OR.(ISPACE == 8)) THEN ! Primary variable
! DEPENDS ON: primary
                    CALL PRIMARY(ISEC,IITM,Im_index,Im_ident,Sm_ident,  &
                         RLEVS,RADDRESS,PIrow,ErrorStatus,CMESSAGE)
                  END IF
                END IF  !  PPXPTR(m,s,i)  /=  0
              END DO    !  Loop over items
            END DO      !  ISEC_loop : Loop over sections
          END IF        !  test whether SM contains IM
        END DO          !  Loop over Im_index
      END DO            !  Loop over SM partitions

! LOOKUP array lengths for TYPSIZE
      A_PROG_LOOKUP = NHEAD(A_IM)
! Primary data lengths for TYPSIZE
      A_PROG_LEN = LPrimIM(A_IM)
      WRITE(6,'(/,A)') '********************************************'// &
                       '***********************************'
      WRITE(6,'(A,I8)') ' ADDRES : A_PROG_LOOKUP = ',A_PROG_LOOKUP
      WRITE(6,'(A,I8)') ' ADDRES : A_PROG_LEN    = ',A_PROG_LEN
      WRITE(6,'(A,2I6,2I9)')' ADDRES : NHEAD,LPrimIM = ',NHEAD,LPrimIM

! 2. Loop through stash list to set output addresses and
!                 header positions for diagnostics
      DO IREC=1,NRECS

! Read internal model number from stash list. Stash list has already
! been ordered by internal model, section, item. Thus, all the atmos
! diagnostic addressing will be done first, followed by the slab
! addressing in the case of a slab model.
        Im_ident = LIST_S(st_model_code,IREC)
! Obtain submodel partition id.
        Sm_ident = SUBMODEL_PARTITION_INDEX(Im_ident)

! Set output address relative to D1
        IF(LIST_S(st_output_code,IREC) == 1) THEN

! Diagnostic output to dump rather than direct output pp file
!   Add the output length for this diag to LDUMP; total length of
!   dump so far = LPRIM + LDUMP; hence obtain the start address for
!   the output from the next diagnostic to be stored in dump.

          LIST_S(st_output_addr,IREC)                                   &
                   = LPRIM(Sm_ident)+LDUMP(Sm_ident)+1
! Information for preliminary D1 addressing array
          N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
          IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
            D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=diag
            D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
            D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IREC
          END IF
          LIST_S(st_dump_output_addr,IREC)=                             &
                   global_LPRIM(Sm_ident)+global_LDUMP(Sm_ident)+1
          LDUMP(Sm_ident)                                               &
                   = LDUMP(Sm_ident)+LIST_S(st_output_length,IREC)
          LDumpIM(Im_ident)                                             &
                   = LDumpIM(Im_ident)+LIST_S(st_output_length,IREC)
          global_LDUMP(Sm_ident)=                                       &
            global_LDUMP(Sm_ident)+LIST_S(st_dump_output_length,IREC)
          global_LDUMPIM(Sm_ident)=                                     &
            global_LDUMPIM(Im_ident)+LIST_S(st_dump_output_length,IREC)

          IF(LIST_S(st_output_bottom,IREC) == 100) THEN
! Special levels
            RLEVS=1
          ELSE IF(LIST_S(st_series_ptr,IREC) /= 0) THEN
! Time series domain
            RLEVS=1
          ELSE IF(LIST_S(st_gridpoint_code,IREC) >= 10                  &
             .AND.LIST_S(st_gridpoint_code,IREC) <  20) THEN
! Vertical ave.
            RLEVS=1
          ELSE IF(LIST_S(st_output_bottom,IREC) <  0) THEN
! Levels list
            RLEVS=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
          ELSE
! Range of model levels
            RLEVS=LIST_S(st_output_top   ,IREC)                         &
                 -LIST_S(st_output_bottom,IREC)+1
          END IF

          IF (LIST_S(st_pseudo_out,IREC) >  0) THEN
! Pseudo levels
            RLEVS=RLEVS*LENPLST(LIST_S(st_pseudo_out,IREC))
          END IF

! Set position of pp lookup header in the dump
          LIST_S(st_lookup_ptr,IREC)=NHeadSub(Sm_ident)+1

! Increment NHEAD (there is one pp header for each level at
!  which a diagnostic is output
          NHEAD   (Im_ident)=NHEAD   (Im_ident)+RLEVS
          NHeadSub(Sm_ident)=NHeadSub(Sm_ident)+RLEVS

        ELSE IF(LIST_S(st_output_code,IREC) == 2) THEN

! Secondary data in D1.
! Compute and store secondary data lengths. Start address for
! secondary data is determined below, after total dump
! diagnostic length has been found.

          LIST_S(st_output_addr,IREC)=LSECD(Sm_ident)+1
          LSECD(Sm_ident)                                               &
         =LSECD(Sm_ident)+LIST_S(st_output_length,IREC)
          LSecdIM(Im_ident)                                             &
         =LSecdIM(Im_ident)+LIST_S(st_output_length,IREC)
! Set pointer for pp header
          LIST_S(st_lookup_ptr,IREC)=-1

        ELSE IF(LIST_S(st_output_code,IREC) <  0) THEN

! Diagnostic output to PP file

! Compute no. of pp headers for this diagnostic
!   = output levels * pseudo output levels * output times

!   No. of levels
          IF(LIST_S(st_output_bottom,IREC) == 100) THEN
! Special levels
            IL=1
          ELSE IF(LIST_S(st_series_ptr,IREC) /= 0) THEN
! Time series dom
            IL=1
          ELSE IF(LIST_S(st_gridpoint_code,IREC) >= 10                  &
            .AND.LIST_S(st_gridpoint_code,IREC) <  20) THEN
! Vertical average
            IL=1
          ELSE IF(LIST_S(st_output_bottom,IREC) <  0) THEN
! Levels list
            IL=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
          ELSE
! Range of mod levs
            IL=LIST_S(st_output_top,IREC)                               &
             -LIST_S(st_output_bottom,IREC)+1
          END IF

!   No. of pseudo levels
          IF (LIST_S(st_pseudo_out,IREC) >  0) THEN
            IP=LENPLST(LIST_S(st_pseudo_out,IREC))
          ELSE
            IP=1
          END IF

!   No. of output times
          IF(LIST_S(st_freq_code,IREC) >  0) THEN
            IFIRST=LIST_S(st_start_time_code,IREC)
            IFREQ =LIST_S(st_freq_code      ,IREC)
            IF(LIST_S(st_end_time_code,IREC) == -1) THEN
! Output to continues to end of run
              IHOURS=1+8760*RUN_TARGET_END(1)                           &
                      + 744*RUN_TARGET_END(2)                           &
                      +  24*RUN_TARGET_END(3)                           &
                      +     RUN_TARGET_END(4)
! DEPENDS ON: totimp
              ILAST=TOTIMP(IHOURS,'H ',Im_ident)
              IF (ILAST  ==  -999) THEN
                errorStatus = 1
                cmessage = 'TOTIMP:UNEXPECTED TIME UNIT '//             &
                    'or IRREGULAR DUMPS FOR DUMP FREQUENCY'
                GO TO 9999
              END IF
            ELSE
! Last output time before end of run
              ILAST=LIST_S(st_end_time_code,IREC)
            END IF

            IT= 1 + (ILAST-IFIRST)/IFREQ
            IF (IT <  0) THEN
              IT=0
              WRITE(6,'(A)')                                            &
            ' Output time error detected in routine ADDRESS:'
              WRITE(6,'(A)')                                            &
            ' Output time starts after specified end of run'
              WRITE(6,'(A,4I6)')                                        &
            ' STASH record no.,MODEL,SECTION,ITEM as follows: ',        &
                              IREC, LIST_S(st_model_code,IREC),         &
                                    LIST_S(st_sect_code ,IREC),         &
                                    LIST_S(st_item_code ,IREC)
              WRITE(6,'(A,I6)') 'OUTPUT CODE: ',                        &
                                    LIST_S(st_output_code,IREC)
            END IF
          ELSE
! Times table in STASH_times array
            IT=1
            DO I=1,NTIMEP
              IF (ITIM_S(I,-LIST_S(st_freq_code,IREC)) == -1) THEN
                IT=I-1
                GOTO 260
              END IF
            END DO
 260        CONTINUE
          END IF
! No. of output "headers" - (levels)*(pseudo-levels)*(output times)
          IH=IL*IP*IT
! Assign output unit no. (nn) to (st_output_addr)
          LIST_S(st_output_addr,IREC)=-LIST_S(st_output_code,IREC)
! Assign no. of output headers to NHEAD_FILE(nn)
          NHEAD_FILE(LIST_S(st_output_addr,IREC))=                      &
          NHEAD_FILE(LIST_S(st_output_addr,IREC)) + IH
        ELSE IF (LIST_S(st_output_code,IREC) == 0) THEN
! Inactive record, not output
          LIST_S(st_output_addr,IREC)=-LIST_S(st_output_code,IREC)
        ELSE
          WRITE(6,'(A)') 'ERROR detected in routine ADDRESS '
          WRITE(6,'(A)') 'ILLEGAL OUTPUT CODE FOR STASH RECORD '
          WRITE(6,'(A,4I6)')                                            &
        ' STASH record no.,MODEL,SECTION,ITEM as follows: ',            &
                         IREC, LIST_S(st_model_code,IREC),              &
                               LIST_S(st_sect_code ,IREC),              &
                               LIST_S(st_item_code ,IREC)
        END IF

      END DO      ! End of loop over records for D1 addressing


!     Correct the addressing of SPACE=9 items from being relative
!     to start of LEXTRA space to being relative to start of dump

!     Loop over submodel partitions
      DO  Sm_ident = 1,N_SUBMODEL_PARTITION_MAX

!       Loop over internal models for each SM partition
        DO Im_index = 1,N_INTERNAL_MODEL

!         Test whether current SM contains this IM
          IF (SM_def(Sm_ident,Im_index) >  0) THEN

!           Obtain internal model identifier
            Im_ident   = INTERNAL_MODEL_LIST(Im_index)

            DO ISEC_loop=1,7

              IF (ISEC_loop  ==  1) ISEC=0  ! section zero primary
              IF (ISEC_loop  ==  2) ISEC=31 ! LBC Input (primary)
              IF (ISEC_loop  ==  3) ISEC=32 ! LBC Output (primary)
              IF (ISEC_loop  ==  4) ISEC=33 ! Tracers (primary)
              IF (ISEC_loop  ==  5) ISEC=34 ! UKCA tracers (primary)
              IF (ISEC_loop  ==  6) ISEC=36 ! Free Tracer lbcs
              IF (ISEC_loop  ==  7) ISEC=37 ! UKCA Tracer lbcs
              DO IITM   = 1,PPXREF_ITEMS
!             Check whether there is a primary field corresponding
                IF (PPXPTR(Im_index,ISEC,IITM) /= 0) THEN
! DEPENDS ON: exppxi
                ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code,     &
                  ErrorStatus,CMESSAGE)
                IF (IN_S(1,Im_ident,ISEC,IITM) /= 0                     &
                                                      ! item is active
                      .AND. ISPACE == 9) THEN
                  IN_S(1,Im_ident,ISEC,IITM)=IN_S(1,Im_ident,ISEC,IITM)+&
     &              LPRIM(Sm_ident)+LDUMP(Sm_ident)
                END IF
                END IF
              END DO
            END DO ! ISEC_loop
          END IF
        END DO
      END DO


! Set secondary data addresses relative to start of D1
      DO IREC=1,NRECS
        Im_ident = LIST_S(st_model_code,IREC)
        Sm_ident = SUBMODEL_PARTITION_INDEX(Im_ident)

        IF (LIST_S(st_output_code,IREC) == 2) THEN
          LIST_S(st_output_addr,IREC)  =LIST_S(st_output_addr,IREC)     &
        + LPRIM(Sm_ident)+LDUMP(Sm_ident)+LEXTRA(Sm_ident)
! Information for preliminary D1 addressing array
          N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
          IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN) THEN
            D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=seco
            D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
            D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IREC
          ENDIF
        END IF
      END DO

! 3.  Set input addresses and work lengths for non-primary
!            fields (i.e., ISPACE=0 or 7)
      DO Im_ident=1,N_INTERNAL_MODEL_MAX
        Sm_ident=  SUBMODEL_PARTITION_INDEX(Im_ident)
        DO ISEC    =0,PPXREF_SECTIONS
! Re-initialise sectional work lengths
          DO I=1,N_SUBMODEL_PARTITION_MAX
            LWORK_S(I)=0
          END DO
          DO IITM  =1,PPXREF_ITEMS
            IF(INDX_S(2,Im_ident,ISEC,IITM) >  0) THEN
! Item in STASH list
! Obtain space code & section zero point-back code
!   from ppxref lookup array
! DEPENDS ON: exppxi
          ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code   ,        &
                                                ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          PTR_PROG= EXPPXI(Im_ident,ISEC,IITM,ppx_ptr_code     ,        &
                                                ErrorStatus,CMESSAGE)
! Compute length of work space required
              IF (ISPACE == 0) THEN
! STASH_WORK address & length
                IN_S(1,Im_ident,ISEC,IITM)=LWORK_S(Sm_ident)+1
                LWORK_S(Sm_ident)=LWORK_S(Sm_ident)                     &
                                  +IN_S(2,Im_ident,ISEC,IITM)
              ELSE IF (ISPACE == 7) THEN
! Point-back to primary space in section 0
                IN_S(1,Im_ident,ISEC,IITM    )                          &
               =IN_S(1,Im_ident,0   ,PTR_PROG)
                IN_S(2,Im_ident,ISEC,IITM    )                          &
               =IN_S(2,Im_ident,0   ,PTR_PROG)
              END IF
            END IF
          END DO   ! Items

! Find max sectional work length for each submodel partition
          DO I=1,N_SUBMODEL_PARTITION_MAX
            LWORK(I)=MAX(LWORK(I),LWORK_S(I))
          END DO

        END DO     ! Sections
        IF(Sm_ident /= 0)THEN
!       Save the maximum value for dimensioning full D1 address array
          N_OBJ_D1_MAX=MAX(N_OBJ_D1_MAX,N_OBJ_D1(Sm_ident))
          WRITE(6,'(I12,A,I6)') N_OBJ_D1(Sm_ident),  &
                                ' D1 items in submodel ',Sm_ident
        END IF
      END DO     ! Models
      IF(N_OBJ_D1_MAX >  MAX_D1_LEN)THEN
        WRITE(6,'(A)')'ADDRES1: No of items in D1 exceeds maximum allowed:'
        WRITE(6,'(A,I12,A,I12)') 'Number allowed ',MAX_D1_LEN,  &
                                 ' Number requested ',N_OBJ_D1_MAX
        WRITE(6,'(A)')'Modify the COMDECK STEXTEND to increase'
        WRITE(6,'(A)')'MAX_D1_LEN parameter as required'
        WRITE(6,'(A)')'Such a change can be safely made'
        CMESSAGE='ADDRES1: No of D1 items exceeds max: See output'
        ErrorStatus=1
      END IF

      WRITE(6,'(A,/)') '********************************************'// &
                       '***********************************'
 9999 CONTINUE
      IF (lhook) CALL dr_hook('ADDRES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ADDRES
