! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module to contain photolysis rate arrays
!   Contains subroutines: ukca_strat_photol_init
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
     MODULE UKCA_DISSOC


     USE yomhook, ONLY: lhook, dr_hook
     USE parkind1, ONLY: jprb, jpim
     USE UM_ParVars
     IMPLICIT NONE

     REAL, ALLOCATABLE, SAVE :: AJHNO3 (:)
     REAL, ALLOCATABLE, SAVE :: AJPNA  (:)
     REAL, ALLOCATABLE, SAVE :: AJH2O2 (:)
     REAL, ALLOCATABLE, SAVE :: AJ2A   (:)
     REAL, ALLOCATABLE, SAVE :: AJ2B   (:)
     REAL, ALLOCATABLE, SAVE :: AJ3    (:)
     REAL, ALLOCATABLE, SAVE :: AJ3A   (:)
     REAL, ALLOCATABLE, SAVE :: AJCNITA(:)
     REAL, ALLOCATABLE, SAVE :: AJCNITB(:)
     REAL, ALLOCATABLE, SAVE :: AJBRNO3(:)
     REAL, ALLOCATABLE, SAVE :: AJBRCL (:)
     REAL, ALLOCATABLE, SAVE :: AJOCLO (:)
     REAL, ALLOCATABLE, SAVE :: AJCL2O2(:)
     REAL, ALLOCATABLE, SAVE :: AJHOCL (:)
     REAL, ALLOCATABLE, SAVE :: AJNO   (:)
     REAL, ALLOCATABLE, SAVE :: AJNO2  (:)
     REAL, ALLOCATABLE, SAVE :: AJN2O5 (:)
     REAL, ALLOCATABLE, SAVE :: AJNO31 (:)
     REAL, ALLOCATABLE, SAVE :: AJNO32 (:)
     REAL, ALLOCATABLE, SAVE :: AJBRO  (:)
     REAL, ALLOCATABLE, SAVE :: AJHCL  (:)
     REAL, ALLOCATABLE, SAVE :: AJN2O  (:)
     REAL, ALLOCATABLE, SAVE :: AJHOBR (:)
     REAL, ALLOCATABLE, SAVE :: AJF11  (:)
     REAL, ALLOCATABLE, SAVE :: AJF12  (:)
     REAL, ALLOCATABLE, SAVE :: AJH2O  (:)
     REAL, ALLOCATABLE, SAVE :: AJCCL4 (:)
     REAL, ALLOCATABLE, SAVE :: AJF113 (:)
     REAL, ALLOCATABLE, SAVE :: AJF22  (:)
     REAL, ALLOCATABLE, SAVE :: AJCH3CL(:)
     REAL, ALLOCATABLE, SAVE :: AJC2OA (:)
     REAL, ALLOCATABLE, SAVE :: AJC2OB (:)
     REAL, ALLOCATABLE, SAVE :: AJMHP  (:)
     REAL, ALLOCATABLE, SAVE :: AJCH3BR(:)
     REAL, ALLOCATABLE, SAVE :: AJMCFM (:)
     REAL, ALLOCATABLE, SAVE :: AJCH4  (:)
     REAL, ALLOCATABLE, SAVE :: AJF12B1(:)
     REAL, ALLOCATABLE, SAVE :: AJF13B1(:)
     REAL, ALLOCATABLE, SAVE :: AJCOF2 (:)
     REAL, ALLOCATABLE, SAVE :: AJCOFCL(:)
     REAL, ALLOCATABLE, SAVE :: AJCO2  (:)
     REAL, ALLOCATABLE, SAVE :: AJCOS  (:)
     REAL, ALLOCATABLE, SAVE :: AJHONO (:)
     REAL, ALLOCATABLE, SAVE :: AJMENA (:)
     REAL, ALLOCATABLE, SAVE :: AJCHBR3(:)
     REAL, ALLOCATABLE, SAVE :: AJDBRM (:)
     REAL, ALLOCATABLE, SAVE :: AJCS2  (:)
     REAL, ALLOCATABLE, SAVE :: AJH2SO4(:)
     REAL, ALLOCATABLE, SAVE :: AJSO3  (:)

     CONTAINS

!  UKCA stratospheric photolysis module
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
!---------------------------------------------------------------------------
! Subroutine UKCA_STRAT_PHOTOL_INIT
!------------------------------------------------------------------------
!
! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.

      SUBROUTINE UKCA_STRAT_PHOTOL_INIT
      USE UM_ParVars
      IMPLICIT NONE

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

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_DISSOC:UKCA_STRAT_PHOTOL_INIT',zhook_in, &
                              zhook_handle)

      IF (.NOT. ALLOCATED(AJHNO3 )) ALLOCATE(AJHNO3 (theta_field_size))
      IF (.NOT. ALLOCATED(AJPNA  )) ALLOCATE(AJPNA  (theta_field_size))
      IF (.NOT. ALLOCATED(AJH2O2 )) ALLOCATE(AJH2O2 (theta_field_size))
      IF (.NOT. ALLOCATED(AJ2A   )) ALLOCATE(AJ2A   (theta_field_size))
      IF (.NOT. ALLOCATED(AJ2B   )) ALLOCATE(AJ2B   (theta_field_size))
      IF (.NOT. ALLOCATED(AJ3    )) ALLOCATE(AJ3    (theta_field_size))
      IF (.NOT. ALLOCATED(AJ3A   )) ALLOCATE(AJ3A   (theta_field_size))
      IF (.NOT. ALLOCATED(AJCNITA)) ALLOCATE(AJCNITA(theta_field_size))
      IF (.NOT. ALLOCATED(AJCNITB)) ALLOCATE(AJCNITB(theta_field_size))
      IF (.NOT. ALLOCATED(AJBRNO3)) ALLOCATE(AJBRNO3(theta_field_size))
      IF (.NOT. ALLOCATED(AJBRCL )) ALLOCATE(AJBRCL (theta_field_size))
      IF (.NOT. ALLOCATED(AJOCLO )) ALLOCATE(AJOCLO (theta_field_size))
      IF (.NOT. ALLOCATED(AJCL2O2)) ALLOCATE(AJCL2O2(theta_field_size))
      IF (.NOT. ALLOCATED(AJHOCL )) ALLOCATE(AJHOCL (theta_field_size))
      IF (.NOT. ALLOCATED(AJNO   )) ALLOCATE(AJNO   (theta_field_size))
      IF (.NOT. ALLOCATED(AJNO2  )) ALLOCATE(AJNO2  (theta_field_size))
      IF (.NOT. ALLOCATED(AJN2O5 )) ALLOCATE(AJN2O5 (theta_field_size))
      IF (.NOT. ALLOCATED(AJNO31 )) ALLOCATE(AJNO31 (theta_field_size))
      IF (.NOT. ALLOCATED(AJNO32 )) ALLOCATE(AJNO32 (theta_field_size))
      IF (.NOT. ALLOCATED(AJBRO  )) ALLOCATE(AJBRO  (theta_field_size))
      IF (.NOT. ALLOCATED(AJHCL  )) ALLOCATE(AJHCL  (theta_field_size))
      IF (.NOT. ALLOCATED(AJN2O  )) ALLOCATE(AJN2O  (theta_field_size))
      IF (.NOT. ALLOCATED(AJHOBR )) ALLOCATE(AJHOBR (theta_field_size))
      IF (.NOT. ALLOCATED(AJF11  )) ALLOCATE(AJF11  (theta_field_size))
      IF (.NOT. ALLOCATED(AJF12  )) ALLOCATE(AJF12  (theta_field_size))
      IF (.NOT. ALLOCATED(AJH2O  )) ALLOCATE(AJH2O  (theta_field_size))
      IF (.NOT. ALLOCATED(AJCCL4 )) ALLOCATE(AJCCL4 (theta_field_size))
      IF (.NOT. ALLOCATED(AJF113 )) ALLOCATE(AJF113 (theta_field_size))
      IF (.NOT. ALLOCATED(AJF22  )) ALLOCATE(AJF22  (theta_field_size))
      IF (.NOT. ALLOCATED(AJCH3CL)) ALLOCATE(AJCH3CL(theta_field_size))
      IF (.NOT. ALLOCATED(AJC2OA )) ALLOCATE(AJC2OA (theta_field_size))
      IF (.NOT. ALLOCATED(AJC2OB )) ALLOCATE(AJC2OB (theta_field_size))
      IF (.NOT. ALLOCATED(AJMHP  )) ALLOCATE(AJMHP  (theta_field_size))
      IF (.NOT. ALLOCATED(AJCH3BR)) ALLOCATE(AJCH3BR(theta_field_size))
      IF (.NOT. ALLOCATED(AJMCFM )) ALLOCATE(AJMCFM (theta_field_size))
      IF (.NOT. ALLOCATED(AJCH4  )) ALLOCATE(AJCH4  (theta_field_size))
      IF (.NOT. ALLOCATED(AJF12B1)) ALLOCATE(AJF12B1(theta_field_size))
      IF (.NOT. ALLOCATED(AJF13B1)) ALLOCATE(AJF13B1(theta_field_size))
      IF (.NOT. ALLOCATED(AJCOF2 )) ALLOCATE(AJCOF2 (theta_field_size))
      IF (.NOT. ALLOCATED(AJCOFCL)) ALLOCATE(AJCOFCL(theta_field_size))
      IF (.NOT. ALLOCATED(AJCO2  )) ALLOCATE(AJCO2  (theta_field_size))
      IF (.NOT. ALLOCATED(AJCOS  )) ALLOCATE(AJCOS  (theta_field_size))
      IF (.NOT. ALLOCATED(AJHONO )) ALLOCATE(AJHONO (theta_field_size))
      IF (.NOT. ALLOCATED(AJMENA )) ALLOCATE(AJMENA (theta_field_size))
      IF (.NOT. ALLOCATED(AJCHBR3)) ALLOCATE(AJCHBR3(theta_field_size))
      IF (.NOT. ALLOCATED(AJDBRM )) ALLOCATE(AJDBRM (theta_field_size))
      IF (.NOT. ALLOCATED(AJCS2  )) ALLOCATE(AJCS2  (theta_field_size))
      IF (.NOT. ALLOCATED(AJH2SO4)) ALLOCATE(AJH2SO4(theta_field_size))
      IF (.NOT. ALLOCATED(AJSO3  )) ALLOCATE(AJSO3  (theta_field_size))
      
      IF (lhook) CALL dr_hook('UKCA_DISSOC:UKCA_STRAT_PHOTOL_INIT',zhook_out, &
                              zhook_handle)
      RETURN
      END SUBROUTINE UKCA_STRAT_PHOTOL_INIT

      END MODULE UKCA_DISSOC
