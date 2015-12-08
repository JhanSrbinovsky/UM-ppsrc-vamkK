! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Convert lat/long specification to row/column specification
! Subroutine Interface:

      SUBROUTINE LLTORC(GRID_CODE,                                      &
     &                  NORTH_LAT_IN,SOUTH_LAT_IN,                      &
     &                  WEST_LONG_IN,EAST_LONG_IN,                      &
     &                  START_ROW,END_ROW,START_COL,END_COL             &



     &                 )
      USE Submodel_Mod
      USE Ereport_Mod, ONLY : Ereport
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE UM_ParVars
      USE Decomp_DB
      USE grdtypes_mod, ONLY: gt_unset, gt_atmos, gt_ocean, gt_wave,     &
                              gt_thetamass, gt_velocity, gt_u_c, gt_v_c, &
                              gt_hybrid, gt_river, gt_allpts, gt_land,   &
                              gt_sea, gt_full, gt_zonal, gt_meridional,  &
                              gt_ozone, gt_scalar, gt_compressed, gt_lbc,&
                              gt_nocyclic, gt_optcyclic, gt_cyclic
! version_mod items required by model.h
      USE version_mod, ONLY: nsectp, outfile_e, outfile_s
      IMPLICIT NONE

! Description:
!   Uses the gridpoint code, the lat/long spec of the required area,
!   and the model area sizes to calculate the model row/column
!   numbers for the area.
!   Fix for the ocean model, since N is at the bottom.
!   Called by PRELIM, INPUTL, ADDRES.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
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


      INTEGER                                                           &
     &  GRID_CODE                                                       &
                    ! IN : Grid type code from STASHmaster

     &, NORTH_LAT_IN                                                    &
                       ! IN : Latitude of Northern boundary
     &, SOUTH_LAT_IN                                                    &
                       ! IN : Latitude of Southern boundary
     &, WEST_LONG_IN                                                    &
                       ! IN : Longitude of Western boundary
     &, EAST_LONG_IN                                                    &
                       ! IN : Longitude of Eastern boundary

     &, START_ROW                                                       &
                    ! OUT : Row number of start of area
     &, END_ROW                                                         &
                    ! OUT : Row number of end of area
     &, START_COL                                                       &
                    ! OUT : Column number of start of area
     &, END_COL     ! OUT : Column number of end of area


! Local variables
      INTEGER                                                           &
     &  MODEL_TYPE                                                      &
                    ! model type of grid
     &, CONTENT                                                         &
                    ! content type of grid
     &, COVERAGE                                                        &
                    ! coverage of grid
     &, DOMAIN                                                          &
                    ! domain of grid
     &, CYCLIC                                                          &
                    ! does grid contain cyclic wrap columns
     &, fld_type    ! P, U or V field type

      INTEGER                                                           &
     &  NORTH_LAT                                                       &
                   ! Modifiable copies
     &, SOUTH_LAT                                                       &
                   ! of the input arguments
     &, WEST_LONG                                                       &
                   ! so that the format can
     &, EAST_LONG                                                       &
                   ! be changed if necessary

     &, NROWS                                                           &
                   ! number of rows
     &, NCOLS      ! number of columns

      LOGICAL                                                           &
     &  LAM_MODEL  ! True if this is a LAM configuration

      REAL                                                              &
     &  POLE_LAT                                                        &
                   ! Latitude of rotated pole (LAM only)
     &, POLE_LONG  ! Longitude of rotated pole (LAM only)

      INTEGER                                                           &
     &  ROW_ORDERING                                                    &
                      ! ordering North->South or South->North
     &, North_to_South                                                  &
                       ! indicates North->South ordering
     &, South_to_North ! indicates South->North ordering

      PARAMETER                                                         &
     &  (North_to_South=1,South_to_North=2)

      REAL                                                              &
     &  START_LAT                                                       &
                   ! Starting latitude
     &, START_LONG                                                      &
                   ! Starting longitude
     &, DEL_LAT                                                         &
                   ! Latitude grid spacing
     &, DEL_LONG                                                        &
                   ! Longitude grid spacing
     &, R_START_ROW                                                     &
                    ! first row number
     &, R_END_ROW                                                       &
                    ! last row number
     &, R_START_COL                                                     &
                    ! first column number
     &, R_END_COL   ! last column number

! Function call
      INTEGER GET_FLD_TYPE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------

! Copy input arguments
      IF (lhook) CALL dr_hook('LLTORC',zhook_in,zhook_handle)
      NORTH_LAT=NORTH_LAT_IN
      SOUTH_LAT=SOUTH_LAT_IN
      WEST_LONG=WEST_LONG_IN
      EAST_LONG=EAST_LONG_IN

      IF (WEST_LONG  <   0) WEST_LONG=WEST_LONG+360
      IF (EAST_LONG  <=  0) EAST_LONG=EAST_LONG+360

! Get information about grid type

! DEPENDS ON: gt_decode
      CALL GT_DECODE(GRID_CODE,                                         &
     &               MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)

! Get information of field type for PARVARS variables

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(GRID_CODE)

! Find start latitude and longitude of full grid, grid spacing and
! size of grid

      IF (MODEL_TYPE  ==  gt_atmos) THEN

        START_LAT=H_A_FIRSTLAT
        START_LONG=H_A_FIRSTLONG
        DEL_LAT=H_A_NSSPACE
        DEL_LONG=H_A_EWSPACE
        NROWS=decompDB(decomp_standard_atmos)%glsize(2,fld_type)
        NCOLS=decompDB(decomp_standard_atmos)%glsize(1,fld_type)
        ROW_ORDERING=South_to_North

        IF (H_GLOBAL(A_IM) == 'N') THEN   ! LAM Configuration
          LAM_MODEL=.TRUE.
          POLE_LAT=H_A_POLELAT
          POLE_LONG=H_A_POLELONG
        ELSE
          LAM_MODEL=.FALSE.
        ENDIF

      ENDIF

! Ensure that DEL_LAT is always positive
! (ie. doesn't imply a direction)
      IF (DEL_LAT  <   0) DEL_LAT=-DEL_LAT

! This assumes the mass grid. The start latitude and longitude
! may be offset for velocity/U/V fields

      IF (CONTENT  ==  gt_velocity) THEN  ! B grid U/C
        START_LAT=START_LAT+DEL_LAT/2
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF (CONTENT  ==  gt_U_C) THEN   ! C grid U
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF( CONTENT  ==  gt_V_C) THEN   ! C grid V
        START_LAT=START_LAT-DEL_LAT/2
      ENDIF

! This assumes full domain. Now take account of zonal,meridional
! and scalar fields

      IF (DOMAIN  ==  gt_zonal) THEN
        NCOLS=1
      ELSEIF (DOMAIN  ==  gt_meridional) THEN
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_scalar) THEN
        NCOLS=1
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_ozone) THEN
        IF (ZonAvOzone) NCOLS=1
      ENDIF

! This is the global sizes - this may be all we need

      IF ((NORTH_LAT  ==  90) .AND. (SOUTH_LAT  ==  -90) .AND.          &
     &    (WEST_LONG  ==  0) .AND. (EAST_LONG  ==  360)) THEN

        START_ROW=1
        END_ROW=NROWS
        START_COL=1
        END_COL=NCOLS

      ELSE ! Not a simple case

! If this is a LAM configuration we need to transform the latitudes
! and longitudes relative to the rotated pole

        IF (LAM_MODEL) THEN
! DEPENDS ON: lltoll
          CALL LLTOLL(                                                  &
     &        NORTH_LAT,SOUTH_LAT,EAST_LONG,WEST_LONG,                  &
     &        POLE_LAT,POLE_LONG)
        ENDIF

! Make sure that DEL_LAT has a sign consistent with the
! grid ordering. Further up we have ensured that DEL_LAT
! is positive. Now we change it negative if necessary.
!               ( North->South => -ive
!                 South->North => +ive )

        IF (ROW_ORDERING  ==  North_to_South) DEL_LAT=-DEL_LAT

! Calculate the start and end rows

        IF (ROW_ORDERING  ==  North_to_South) THEN

          R_START_ROW=1.0   + (NORTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (START_LAT-SOUTH_LAT)/DEL_LAT

        ELSE ! South->North

          R_START_ROW=1.0   + (SOUTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (NORTH_LAT-START_LAT)/DEL_LAT

        ENDIF

        START_ROW=MAX(1.0,R_START_ROW)
        END_ROW=MIN(REAL(NROWS),R_END_ROW)

        IF (START_LONG  <   0) START_LONG=START_LONG+360.0

        IF ((REAL(START_LONG) + DEL_LONG*(NCOLS-1))  <=  360.0)         &
     &  THEN  ! If the total model domain doesn't cross
              ! the zero meridian

          R_START_COL=1.0 +   (WEST_LONG-START_LONG)/DEL_LONG
          R_END_COL=  1.999 + (EAST_LONG-START_LONG)/DEL_LONG

        ELSE ! model domain crosses the zero meridian

          IF (REAL(WEST_LONG)  <   START_LONG) THEN
            ! If the start is before the zero meridian

            R_START_COL=1.0  + (WEST_LONG + 360.0 - START_LONG)/        &
     &                         DEL_LONG

          ELSE ! the start lies after the meridian

            R_START_COL=1.0 + (REAL(WEST_LONG)-START_LONG)/             &
     &                        DEL_LONG

          ENDIF

          IF (REAL(EAST_LONG)  <   START_LONG) THEN
            ! If the end is after the zero meridian

            R_END_COL=1.0 + (EAST_LONG + 360.0 - START_LONG)/           &
     &                      DEL_LONG

          ELSE ! the end lies before the zero meridian

            R_END_COL=1.0 + (REAL(EAST_LONG) - START_LONG)/             &
     &                      DEL_LONG

          ENDIF

        ENDIF

        START_COL=MIN(REAL(NCOLS),MAX(1.0,R_START_COL))
        END_COL=MIN(REAL(NCOLS),MAX(1.0,R_END_COL))

      ENDIF

      IF (lhook) CALL dr_hook('LLTORC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LLTORC
