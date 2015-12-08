! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Read in and check a single IAU increment field

MODULE readiaufield_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ReadIAUField (                   &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                          FixHd,            & ! in
                          Len1Lookup,       & ! in
                          Len2Lookup,       & ! in
                          Lookup,           & ! inout
                          IncNum,           & ! in
                          FieldNum,         & ! in
                          LocFldLen,        & ! in
                          FirstCallThisInc, & ! in
                          Field )             ! out

! Method:
!
!   1. Check field dimensions for compatibility with the model.
!   2. Read in field from IAU increment file.
!   3. If required, reset polar rows to their mean values.
!   4. If required, write basic field stats to standard output.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE dynamics_grid_mod, ONLY: l_vatpoles

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IAU_mod, ONLY :   &
    IAU_unit,         &
    IAU_NumFldCodes,  &
    IAU_FldCodes,     &
    IAU_FldDescs,     &
    L_IAU_IncDiags,   &
    L_IAU_ResetPoles

USE global_2d_sums_mod, ONLY: &
    global_2d_sums

USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE Control_Max_Sizes
USE domain_params
USE lookup_addresses
USE um_input_control_mod,  ONLY: model_domain

USE fieldstats_mod, ONLY: fieldstats
USE cppxref_mod, ONLY: ppx_grid_type, ppx_atm_tall, ppx_halo_type
USE ppxlook_mod, ONLY: ppxi,ppxptr
USE Submodel_Mod
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

! Subroutine arguments:

INTEGER, INTENT(IN)    :: FixHd(Len_FixHd) ! Fixed-length header
INTEGER, INTENT(IN)    :: Len1Lookup       ! First  dimension of lookup header
INTEGER, INTENT(IN)    :: Len2Lookup       ! Second dimension of lookup header

INTEGER, INTENT(INOUT) :: Lookup(Len1Lookup,Len2Lookup) ! Lookup header

INTEGER, INTENT(IN)    :: IncNum           ! Increment file number
INTEGER, INTENT(IN)    :: FieldNum         ! Field number in IAU increment file
INTEGER, INTENT(IN)    :: LocFldLen        ! Local length of field
LOGICAL, INTENT(IN)    :: FirstCallThisInc ! Is this the first call for this
                                           ! increment file?
REAL,    INTENT(OUT)   :: Field(LocFldLen) ! Local part of field

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ReadIAUField'

! Local variables:

INTEGER :: i
INTEGER :: isec
INTEGER :: item
INTEGER :: Code
INTEGER :: Code_orig
INTEGER :: Level
INTEGER :: BufLen
INTEGER :: Inc_rows
INTEGER :: Inc_cols
INTEGER :: Model_rows
INTEGER :: Model_cols
INTEGER :: Model_size
INTEGER :: ICode
INTEGER :: IM_index
INTEGER :: ppx_row
INTEGER :: halo_type_orig
INTEGER :: s_addr
INTEGER :: e_addr
INTEGER :: grid_type_orig

LOGICAL :: Code_changed

REAL :: polar_sum(1)
REAL :: polar_row(row_length, 1)
REAL :: polar_mean
REAL :: Global_max
REAL :: Global_min
REAL :: Global_mean
REAL :: Global_RMS

CHARACTER(LEN=800) :: CMessage
CHARACTER(LEN=7)   :: FldDesc

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook('READIAUFIELD',zhook_in,zhook_handle)

Code   = Lookup(ITEM_CODE, FieldNum)
Level  = Lookup(LBLEV,     FieldNum)
BufLen = Lookup(LBLREC,    FieldNum)

!-------------------------------------------------------------------------------
! [1]: Check field dimensions for compatibility with the model.
!-------------------------------------------------------------------------------

Inc_rows = Lookup(LBROW, FieldNum)
Inc_cols = Lookup(LBNPT, FieldNum)

Model_rows = global_rows
Model_cols = global_row_length

IF (l_vatpoles) THEN
IF (Code == 3) Model_rows = global_rows+1 ! v
ELSE
IF (Code == 3) Model_rows = global_rows-1 ! v
END IF ! vatpoles

Model_size = Model_rows * Model_cols

! Check field dimension compatibility with model
IF (Inc_rows == 0 .AND. Inc_cols == 0) THEN
  ! grid-type on Land/Sea points only

  IF (BufLen /= global_land_field) THEN

    ICode = 1
    CMessage(1:80)   = 'Field dimension mis-match.'
    CMessage(81:160) = ''
    WRITE (CMessage(161:240),*) '  Increment no.:     ', IncNum
    WRITE (CMessage(241:320),*) '  Field no.:         ', FieldNum
    WRITE (CMessage(321:400),*) '  STASH code:        ', Code
    WRITE (CMessage(401:480),*) '  Level:             ', Level
    WRITE (CMessage(481:560),*) '  Data land-points:  ', BufLen
    WRITE (CMessage(561:640),*) '  Model land-points: ', global_land_field
    WRITE (CMessage(641:800),*) ''

    CALL EReport (RoutineName, ICode, CMessage)

  END IF

ELSE

  IF (Inc_rows /= Model_rows .OR.  &
      Inc_cols /= Model_cols) THEN

    ICode = 1
    CMessage(1:80)   = 'Field dimension mis-match.'
    CMessage(81:160) = ''
    WRITE (CMessage(161:240),*) '  Increment no.:  ', IncNum
    WRITE (CMessage(241:320),*) '  Field no.:      ', FieldNum
    WRITE (CMessage(321:400),*) '  STASH code:     ', Code
    WRITE (CMessage(401:480),*) '  Level:          ', Level
    WRITE (CMessage(481:560),*) '  Increment rows: ', Inc_rows
    WRITE (CMessage(561:640),*) '  Increment cols: ', Inc_cols
    WRITE (CMessage(641:720),*) '  Model     rows: ', Model_rows
    WRITE (CMessage(721:800),*) '  Model     cols: ', Model_cols

    CALL EReport (RoutineName, ICode, CMessage)

  END IF
END IF

!-------------------------------------------------------------------------------
! [2]: Read in field from IAU increment file.
!-------------------------------------------------------------------------------

! readflds makes use of the PPXI data for the field being read in. However,
! this will only have been set up for qT if it has been requested as a
! STASH diagnostic. To get around this, we temporarily change its STASH code so
! that it is read in as if it were specific humidity.
IF (Code == 16207 .OR. Code == 18001) THEN
  Code_orig = Code
  Code      = 10
  Lookup(ITEM_CODE, FieldNum) = Code
  Code_changed = .TRUE.
ELSE
  Code_changed = .FALSE.
END IF

! Temporarily change PPXI halo entry so that haloes are not read:
IM_index       = INTERNAL_MODEL_INDEX(A_IM)
isec           = Code / 1000
item           = MOD(Code, 1000)
ppx_row        = PPXPTR(IM_index, isec, item)
halo_type_orig = PPXI(ppx_row, ppx_halo_type)
PPXI(ppx_row, ppx_halo_type) = halo_type_no_halo
! When SMC increments are written as ancillaries (i.e. all theta points)
! we must temporarily adjust the grid-type in the PPXI array
grid_type_orig = PPXI(ppx_row, ppx_grid_type)
IF (code == 9 .AND. buflen == Model_size)                               &
    ppxi(ppx_row, ppx_grid_type) = ppx_atm_tall

! DEPENDS ON: readflds
CALL readflds ( IAU_unit,   & ! in
                1,          & ! in
                FieldNum,   & ! in
                Lookup,     & ! in
                Len1Lookup, & ! in
                Field,      & ! out
                BufLen,     & ! in
                FixHd,      & ! in
                ICode,      & ! out
                CMessage )    ! out

IF (ICode > 0) THEN
  WRITE (6,*) 'ReadIAUField: Error reading IAU field no. ', FieldNum
  CALL EReport (RoutineName, ICode, CMessage)
END IF

! Restore PPXI halo entry:
PPXI(ppx_row, ppx_halo_type) = halo_type_orig

! Restore STASH code:
IF (Code_changed) THEN
  Code = Code_orig
  Lookup(ITEM_CODE, FieldNum) = Code
END IF

!-------------------------------------------------------------------------------
! [3]: If required, reset polar rows to their mean values.
!-------------------------------------------------------------------------------

IF (L_IAU_ResetPoles          .AND. &
    Model_domain == mt_global .AND. &
    PPXI(ppx_row, ppx_grid_type) == ppx_atm_tall) THEN

  IF (at_extremity(PNorth)) THEN

    s_addr = 1 + row_length * (rows-1)
    e_addr = s_addr + row_length - 1

    polar_row(:,1) = Field(s_addr:e_addr)

    CALL global_2d_sums(polar_row, row_length, 1, 0, 0, 1, &
                        polar_sum, gc_proc_row_group)

    polar_mean = polar_sum(1) / REAL(global_row_length)

    Field(s_addr:e_addr) = polar_mean

  END IF ! (at_extremity(PNorth))

  IF (at_extremity(PSouth)) THEN

    s_addr = 1
    e_addr = s_addr + row_length - 1

    polar_row(:,1) = Field(s_addr:e_addr)

    CALL global_2d_sums(polar_row, row_length, 1, 0, 0, 1, &
                        polar_sum, gc_proc_row_group)

    polar_mean = polar_sum(1) / REAL(global_row_length)

    Field(s_addr:e_addr) = polar_mean

  END IF ! (at_extremity(PSouth))

END IF

!-------------------------------------------------------------------------------
! [4]: If required, write basic field stats to standard output.
!-------------------------------------------------------------------------------

IF (L_IAU_IncDiags) THEN

  IF (FirstCallThisInc) THEN
    WRITE (6,*) ''
    WRITE (6,'(A,I2.2,A)') &
                ' ReadIAUField: Summary of fields for IAU increment no.', &
                IncNum, ':'
    WRITE (6,*) ''
    WRITE (6,*) '  Field   Level  Max          Min         ' &
                              //' Mean         RMS'
    WRITE (6,*) '  -----   -----  ---          ---         ' &
                              //' ----         ---'
  END IF

  CALL FieldStats (                               &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                    LocFldLen,                    & ! in
                    Field,                        & ! in
                    PPXI(ppx_row, ppx_grid_type), & ! in
                    halo_type_no_halo,            & ! in
                    Global_max,                   & ! out
                    Global_min,                   & ! out
                    Global_mean,                  & ! out
                    Global_RMS )                    ! out

  ! Get field description:
  DO i = 1, IAU_NumFldCodes
    IF (IAU_FldCodes(i) == Code) FldDesc = IAU_FldDescs(i)
  END DO

  WRITE (6,'(3A,I4,A,4(A,ES12.5))')                                     &
    '   ', FldDesc, ' ', Level, ' ', ' ', Global_max,  ' ', Global_min, &
                                     ' ', Global_mean, ' ', Global_RMS

END IF

! Restore PPXI grid_type entry:
PPXI(ppx_row, ppx_grid_type) = grid_type_orig

IF (lhook) CALL dr_hook('READIAUFIELD',zhook_out,zhook_handle)
RETURN

END SUBROUTINE ReadIAUField
END MODULE readiaufield_mod
