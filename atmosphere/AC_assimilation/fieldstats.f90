! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate basic statistics for a horizontal field

MODULE fieldstats_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE FieldStats (              &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                        Field_len,   & ! in
                        Field,       & ! in
                        grid_type,   & ! in
                        halo_type,   & ! in
                        Global_max,  & ! out
                        Global_min,  & ! out
                        Global_mean, & ! out
                        Global_RMS )   ! out

! Description:
!
!   Calculate basic statistics for a horizontal field on u, v, theta or land
!   points.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE trignometric_mod, ONLY : &
    cos_theta_latitude,      &
    cos_v_latitude

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE Control_Max_Sizes
USE cppxref_mod, ONLY :       &
    ppx_atm_tall,             &
    ppx_atm_cuall,            &
    ppx_atm_cvall,            &
    ppx_atm_compressed

IMPLICIT NONE

! Common blocks:

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
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Subroutine arguments:

INTEGER, INTENT(IN)  :: Field_len
REAL,    INTENT(IN)  :: Field(Field_len) ! Horizontal field
INTEGER, INTENT(IN)  :: grid_type
INTEGER, INTENT(IN)  :: halo_type
REAL,    INTENT(OUT) :: Global_max
REAL,    INTENT(OUT) :: Global_min
REAL,    INTENT(OUT) :: Global_mean      ! Area-weighted mean
REAL,    INTENT(OUT) :: Global_RMS       ! Area-weighted RMS

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FieldStats'

! Local variables:

INTEGER :: i, j, p
INTEGER :: ICode
INTEGER :: xm, ym
INTEGER :: xh, yh
INTEGER :: Expected_len

LOGICAL :: ThetaRows
LOGICAL :: VPoints
LOGICAL :: LandPoints

REAL :: Value
REAL :: Weight

REAL :: Global_sum
REAL :: Global_sumsq
REAL :: Global_sumwts

REAL :: ReshapedField(row_length,rows+1)
REAL :: WeightedField(row_length,rows+1)
REAL :: FieldWeights (row_length,rows+1)

REAL :: new_sum(3)

CHARACTER(LEN=320) :: CMessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook('FIELDSTATS',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! [1]: Do some checks based on grid and halo type.
!-------------------------------------------------------------------------------

ThetaRows  = .FALSE.
VPoints    = .FALSE.
LandPoints = .FALSE.

IF (grid_type == ppx_atm_tall .OR. &
    grid_type == ppx_atm_cuall) THEN
  ThetaRows = .TRUE.
  xm = row_length
  ym = rows
ELSE IF (grid_type == ppx_atm_cvall) THEN
  VPoints = .TRUE.
  xm = row_length
  ym = n_rows
ELSE IF (grid_type == ppx_atm_compressed) THEN
  LandPoints = .TRUE.
  Expected_len = field_len
ELSE
  ICode = 1
  CMessage = 'Grid type not catered for.'
  CALL EReport ( RoutineName, ICode, CMessage )
END IF

IF (.NOT.LandPoints) THEN

  IF (halo_type == halo_type_single) THEN
    xh = offx
    yh = offy
  ELSE IF (halo_type == halo_type_extended) THEN
    xh = halo_i
    yh = halo_j
  ELSE IF (halo_type == halo_type_no_halo) THEN
    xh = 0
    yh = 0
  ELSE
    ICode = 1
    CMessage = 'Invalid halo type.'
    CALL EReport ( RoutineName, ICode, CMessage )
  END IF

  Expected_len = (xm + 2*xh) * (ym + 2*yh)

END IF

IF (Field_len /= Expected_len) THEN
  ICode = 1
  CMessage(1:80)   = 'Supplied local field length incorrect.'
  CMessage(81:160) = ''
  WRITE (CMessage(161:240),*) '  supplied len:  ', Field_len
  WRITE (CMessage(241:320),*) '  expected len:  ', Expected_len
  CALL EReport ( RoutineName, ICode, CMessage )
END IF

!-------------------------------------------------------------------------------
! [2]: Get local statistics for fields on u, v, or theta points.
!-------------------------------------------------------------------------------

IF (.NOT.LandPoints) THEN

  DO j = 1, rows
    DO i = 1, row_length
      ReshapedField(i,j) = 0.0
      WeightedField(i,j) = 0.0
      FieldWeights (i,j) = 0.0
    END DO
  END DO

  IF (halo_type == halo_type_no_halo) THEN
    DO j = 1, ym
      DO i = 1, xm
        ReshapedField(i,j) = Field(i + (j-1)*xm)
      END DO
    END DO
  ELSE
    DO j = 1, ym
      DO i = 1, xm
        ReshapedField(i,j) = Field(i+xh + (yh+j-1)*(xm+2*xh))
      END DO
    END DO
  END IF

  IF (ThetaRows) THEN
    DO j = 1, ym
      DO i = 1, xm
        FieldWeights(i,j) = cos_theta_latitude(i,j)
      END DO
    END DO
  ELSE
    DO j = 1, ym
      DO i = 1, xm
        FieldWeights(i,j) = cos_v_latitude(i,j)
      END DO
    END DO
  END IF

  global_max = -HUGE(1.0)
  global_min =  HUGE(1.0)
  new_sum(:) = 0.0

  ! Ignore undefined values in statistics
  DO j = 1, ym
    DO i = 1, xm
      IF (ReshapedField(i,j) == rmdi) FieldWeights(i,j) = 0.0
      global_min = MIN(global_min, ReshapedField(i,j) )
      global_max = MAX(global_max, ReshapedField(i,j) )
      WeightedField(i,j) = ReshapedField(i,j) &
                         * FieldWeights (i,j)
      new_sum(1) = new_sum(1) + WeightedField(i,j)
      new_sum(2) = new_sum(2) + ReshapedField(i,j) &
                              * WeightedField(i,j)
      new_sum(3) = new_sum(3) + FieldWeights (i,j)
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [3]: Get local statistics for fields on land points.
!-------------------------------------------------------------------------------

IF (LandPoints) THEN

  global_max = -HUGE(1.0)
  global_min =  HUGE(1.0)
  new_sum(:) = 0.0

  DO p = 1, land_field
    i = 1 + MOD(land_index(p)-1, row_length)
    j = 1 + (land_index(p)-1) / row_length
    Value      = Field(p)
    Weight     = cos_theta_latitude(i,j)
    global_max = MAX(global_max, Value)
    global_min = MIN(global_min, Value)
    new_sum(1) = new_sum(1) + Weight * Value
    new_sum(2) = new_sum(2) + Weight * Value**2
    new_sum(3) = new_sum(3) + Weight
  END DO

END IF

!-------------------------------------------------------------------------------
! [4]: Get global statistics.
!-------------------------------------------------------------------------------

CALL gc_rmax(1, nproc, Icode, Global_max)
CALL gc_rmin(1, nproc, Icode, Global_min)
CALL gc_rsum(3, nproc, Icode, new_sum)

Global_mean =       new_sum(1) / new_sum(3)
Global_RMS  = SQRT( new_sum(2) / new_sum(3) )

IF (lhook) CALL dr_hook('FIELDSTATS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE FieldStats
END MODULE fieldstats_mod
