! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Defines submodel and section/version configuration
!
! Subroutine Interface:
      SUBROUTINE SETMODL(ErrorStatus,CMESSAGE)


      USE rad_input_mod, ONLY:                                          &
          A_SW_RadStep,                                                 &
          A_LW_RadStep,                                                 &
          A_SW_RadStep_diag,                                            &
          A_SW_Radstep_prog,                                            &
          A_LW_RadStep_diag,                                            &
          A_LW_Radstep_prog,                                            &
          H_SWBands,                                                    &
          H_LWBands

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE UM_ParVars
      USE Decomp_DB
      USE domain_params
      USE version_mod, ONLY: nsectp 
      USE umsections_mod, ONLY: atmos_sr
      USE um_input_control_mod,  ONLY: phenol_period, triffid_period
! version_mod items required by model.h
      USE version_mod, ONLY: nsectp, outfile_e, outfile_s
      USE nlstgen_mod, ONLY : meanfreqim
      USE Submodel_Mod

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

! Array arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE    ! Error return message

! Error status:
      INTEGER        ErrorStatus ! +ve = fatal error


! Local scalars
      REAL    ASteps !Atmos timesteps per hour
      INTEGER I,J
      INTEGER Im,Is
      INTEGER Obs_type(MAX_AOBS)

! Local parameters
      INTEGER NOASFLDS
      PARAMETER(NOASFLDS=4)

! Function and subroutine calls:
      INTEGER  CTOITST

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

!  Define submodel configuration

      IF (lhook) CALL dr_hook('SETMODL',zhook_in,zhook_handle)
      H_ATMOS = 'N'

      DO I = 1,N_INTERNAL_MODEL_MAX
        IF (INTERNAL_MODEL_LIST(I) == ATMOS_IM) THEN
          H_ATMOS ='Y'
        END IF
      END DO

        DO I = 1,N_INTERNAL_MODEL_MAX
          H_VERS(I,0)= 1
        END DO

      IF (H_ATMOS  ==  'Y') THEN
        Im = ATMOS_IM
        DO Is = 1,NSECTP
! DEPENDS ON: ctoitst
          H_VERS(Im,Is)= CTOITST(ATMOS_SR(Is)(1:1))
        END DO
      END IF
! Submodel independent
      DO I=1,N_INTERNAL_MODEL_MAX
        MEAN_NUMBER(I)=0
        DO J=1,4
          IF (meanfreqim(J,I) >  0) THEN
            MEAN_NUMBER(I)=MEAN_NUMBER(I)+1
          END IF
        END DO
      END DO

      Im = ATMOS_IM
      IF (H_ATMOS  ==  'Y') THEN
! Atmos model included

       IF (OCALB == 1) THEN
          H_FLOOR=lfloor
          H_STRAT='N'
        ELSE
          H_FLOOR='Y'
          H_STRAT='Y'
        END IF

        DO I=1,A_MAX_TRVARS   ! Up to A_MAX_TRVARS=150 free tracers
          IF (TCA(I) == 0) THEN
            TRACER_A(I)=.FALSE.
          ELSE
            TRACER_A(I)=.TRUE.
          END IF
        END DO

        DO I=1,A_MAX_UKCAVARS ! Up to A_MAX_UKCAVARS=150 UKCA tracers
          IF (TC_UKCA(I) == 0) THEN
            TR_UKCA_A(I)=.FALSE.
          ELSE
            TR_UKCA_A(I)=.TRUE.
          END IF
        END DO

! Set switches & pseudo level limits for atmos assimilation diags
        DO I=1,AASSETS
          AASSET  (I)=.FALSE.
          Obs_type(I)=0
        END DO
        DO I=1,MAX_AOBS
!   Obtain first digit of obs type
          Obs_type(I)=AOBINC(I)/100
          IF(Obs_type(I) > 0)THEN
            AASSET(Obs_type(I))=.TRUE.
            IF (I >  1) THEN
              IF (Obs_type(I) /= Obs_type(I-1)) THEN
                AASPF (Obs_type(I))=AOBGRP(I)
                AASPL (Obs_type(I))=AOBGRP(I)
              END IF
            ELSE
              AASPF (Obs_type(I))=AOBGRP(I)
              AASPL (Obs_type(I))=AOBGRP(I)
            END IF
            IF (J >  AASSETS) J=AASSETS
            IF (I >  1      ) THEN
              IF (Obs_type(I) == Obs_type(I-1)) THEN
                IF (AOBGRP(I) <  AASPF(Obs_type(I)))                    &
                    AASPF(Obs_type(I))=AOBGRP(I)
                IF (AOBGRP(I) >  AASPL(Obs_type(I)))                    &
                    AASPL(Obs_type(I))=AOBGRP(I)
              END IF
            END IF
          END IF
        END DO

        H_OROG_GRAD = (H_VERS(Im,6) /= 0)


        IF (OCAAA == 1) THEN
! Atmos global model
          H_GLOBAL(A_IM)='Y'

          H_A_EWSPACE=360.0/                                            &
            decompDB(decomp_standard_atmos)%glsize(1,fld_type_p)
          H_A_NSSPACE=180.0/                                            &
            (decompDB(decomp_standard_atmos)%glsize(2,fld_type_p)-1)

          H_A_FIRSTLAT=-90.0       ! S to N
          H_A_FIRSTLONG=0.0
          H_A_POLELAT=90.0
          H_A_POLELONG=0.0
          LMESO=.FALSE.
        ELSE IF (OCAAA == 2 .OR. OCAAA == 3 .OR. OCAAA == 4 ) THEN
! Atmos LAM
          H_GLOBAL(A_IM)='N'
          LMESO=(MESO == 'Y')
          H_A_EWSPACE=EWSPACEA
          H_A_NSSPACE=NSSPACEA
          H_A_FIRSTLAT=FRSTLATA
          H_A_FIRSTLONG=FRSTLONA
          IF (H_A_FIRSTLONG <  0) H_A_FIRSTLONG=H_A_FIRSTLONG+360.0
          H_A_POLELAT=POLELATA
          H_A_POLELONG=POLELONA
        ELSE IF (OCAAA == 5) THEN
! Atmos single column
          H_GLOBAL(A_IM)='N'
          H_A_EWSPACE=360.
          H_A_NSSPACE=180.
          H_A_FIRSTLAT=LATS
          H_A_FIRSTLONG=LONS
          H_A_POLELAT=90.0
          H_A_POLELONG=0.0
          LMESO=.FALSE.
        ELSE IF (OCAAA /= 0) THEN
          write(6,*)                                                    &
         'Setmodl: UNEXPECTED ATMOSPHERIC AREA CODE OCAAA',OCAAA
        END IF

        LEXTRA(A_SM)=(1+MODEL_LEVELS+2*WET_LEVELS)*ROWS*ROW_LENGTH

      ELSE   ! Atmosphere model not included

        ZonAvOzone    =.FALSE.
        H_OROG_GRAD   =.FALSE.
        LMESO         =.FALSE.
        H_FLOOR       ='N'
        H_STRAT       ='N'
        H_GLOBAL(A_IM)='N'

        LAND_FIELD    =0
        ROWS          =0
        MODEL_LEVELS  =0
        WET_LEVELS    =0
        CLOUD_LEVELS  =0
        TR_LEVELS     =0
        TR_VARS       =0
        DO I=1,A_MAX_TRVARS   ! Up to A_MAX_TRVARS=150 free tracers
          TRACER_A(I) =.FALSE.
        END DO
        DO I=1,A_MAX_UKCAVARS ! Up to A_MAX_UKCAVARS=150 UKCA tracers
          TR_UKCA_A(I) =.FALSE.
        END DO
        StLevGWdrag   =0
        BotVDiffLev   =0
        TopVDifflev   =0
        A_SW_RADSTEP  =0
        A_LW_RADSTEP  =0
        A_SW_RADSTEP_DIAG  = 0
        A_SW_RADSTEP_PROG  = 0
        A_LW_RADSTEP_DIAG  = 0
        A_LW_RADSTEP_PROG  = 0
        H_OROG_ROUGH  =0
!       A_CONV_STEP   =0       ! Now hard-wired parameter
        H_SWBANDS     =0
        H_LWBANDS     =0
        H_A_EWSPACE   =0.
        H_A_NSSPACE   =0.
        H_A_FIRSTLAT  =0.
        H_A_FIRSTLONG =0.
        H_A_POLELAT   =0.
        H_A_POLELONG  =0.
!       PHENOL_PERIOD =0.      ! Now hard-wired parameter
!       TRIFFID_PERIOD=0.      ! Now hard-wired parameter
      END IF


      CLOSE(UNIT=10)
      IF (lhook) CALL dr_hook('SETMODL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SETMODL
