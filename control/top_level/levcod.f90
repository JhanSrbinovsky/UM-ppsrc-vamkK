! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Decode the STASH level code

!******************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to rcf_level_code_mod.F90 and rcf_address_mod.F90
!******************************************************************

! Subroutine Interface:
      SUBROUTINE LEVCOD(ILIN,ILOUT,ErrorStatus,CMESSAGE)

      USE rad_input_mod, ONLY: h_swbands,h_lwbands
      USE dynamics_input_mod, ONLY: l_endgame

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE domain_params
! version_mod items required by model.h
      USE version_mod, ONLY: nsectp, outfile_e, outfile_s
      USE Submodel_Mod

      IMPLICIT NONE
! Description:
!   Sets ILOUT to an appropriate level size according to the value of IL
!   Level sizes are parametrised in comdeck MODEL.
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
      INTEGER ILIN    ! Model level code

!   Scalar arguments with intent(out):
      INTEGER ILOUT   ! An actual level
      CHARACTER(LEN=80) CMESSAGE

! Local scalars:
      INTEGER I
      INTEGER J

! ErrorStatus:
      INTEGER ErrorStatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

      IF (lhook) CALL dr_hook('LEVCOD',zhook_in,zhook_handle)
      IF(ILIN == 1) THEN
! First atmos level
        ILOUT=1
      ELSE IF(ILIN == 2) THEN
! Top atmos level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 3) THEN
! Top wet level
        ILOUT=WET_LEVELS
      ELSE IF(ILIN == 4) THEN
        ILOUT=MODEL_LEVELS-1
      ELSE IF(ILIN == 5) THEN
! First boundary layer level
        ILOUT=MIN(1,BL_LEVELS)
      ELSE IF(ILIN == 6) THEN
! Last boundary layer level
        ILOUT=BL_LEVELS
      ELSE IF(ILIN == 7) THEN
        ILOUT=BL_LEVELS+1
      ELSE IF(ILIN == 8) THEN
! First soil level
        ILOUT=MIN(1,ST_LEVELS)
      ELSE IF(ILIN == 9) THEN
! Last soil level
        ILOUT=ST_LEVELS
      ELSE IF(ILIN == 10) THEN
! First tracer level
        ILOUT=MODEL_LEVELS-TR_LEVELS+1
      ELSE IF(ILIN == 11) THEN
! Last tracer level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 12) THEN
        ILOUT=MODEL_LEVELS+1
      ELSE IF(ILIN == 13) THEN
! First gravity wave drag level
        ILOUT=StLevGWdrag
      ELSE IF(ILIN == 14) THEN
! Last gravity wave drag level
        ILOUT=MODEL_LEVELS
      ELSE IF(ILIN == 15) THEN
        ILOUT=BotVDiffLev
      ELSE IF(ILIN == 16) THEN
        ILOUT=TopVDiffLev-1
      ELSE IF(ILIN == 17) THEN
        ILOUT=TopVDiffLev
      ELSE IF(ILIN == 18) THEN
        ILOUT=BL_LEVELS-1
      ELSE IF(ILIN == 19) THEN
        ILOUT=MODEL_LEVELS+1
      ELSE IF(ILIN == 20) THEN
        ILOUT=MIN(2,ST_LEVELS)
      ELSE IF(ILIN == 21) THEN
        ILOUT=1
      ELSE IF(ILIN == 22) THEN
!       ILOUT=NLEVSO - ocean levels no longer supported
        ILOUT=1
      ELSE IF(ILIN == 23) THEN
        ILOUT=OZONE_LEVELS
      ELSE IF(ILIN == 24) THEN
        ILOUT=MODEL_LEVELS*H_SWBANDS
      ELSE IF(ILIN == 25) THEN
        ILOUT=(MODEL_LEVELS+1)*H_SWBANDS
      ELSE IF(ILIN == 26) THEN
        ILOUT=WET_LEVELS*H_SWBANDS
      ELSE IF(ILIN == 27) THEN
        ILOUT=MODEL_LEVELS*H_LWBANDS
      ELSE IF(ILIN == 28) THEN
        ILOUT=(MODEL_LEVELS+1)*H_LWBANDS
      ELSE IF(ILIN == 29) THEN
        ILOUT=WET_LEVELS*H_LWBANDS
      ELSE IF(ILIN == 30) THEN
        ILOUT=2
      ELSE IF(ILIN == 32) THEN
        ILOUT=H_SWBANDS
      ELSE IF(ILIN == 33) THEN
        ILOUT=H_LWBANDS
      ELSE IF(ILIN == 34) THEN
        ILOUT=SM_LEVELS
      ELSE IF(ILIN == 35) THEN
        ILOUT=CLOUD_LEVELS
      ELSE IF(ILIN == 38) THEN
! Surface level on Charney-Phillips theta grid
        ILOUT=0
      ELSE IF(ILIN == 39) THEN
! Number of ISCCP simulator levels
        ILOUT=7
      ELSE IF(ILIN == 40) THEN
! Fields which are different between ENDGAME and New Dynamics regarding the
! whether it starts at theta level 0 in ENDGAME but 1 in New Dynamics.
        IF (l_endgame) THEN
          ILOUT=0
        ELSE
          ILOUT=1
        END IF
      ELSE
        WRITE(6,'(A35,I6)')'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ILIN
        WRITE(CMESSAGE,'(A35,I6)')'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ILIN
        ErrorStatus=1
      END IF

      IF (lhook) CALL dr_hook('LEVCOD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LEVCOD
