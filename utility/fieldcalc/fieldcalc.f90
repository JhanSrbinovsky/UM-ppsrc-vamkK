! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program for post-processing of UM output

PROGRAM FieldCalc

! Description:
!   This program performs a range of functions required for the post-
!   processing of UM output.  It can rewrite some or all of the input
!   fields, calculate diagnostics using standard level and model level
!   fields and scale wind fields.
!   The functions it performs are dependent on the contents of the user-
!   defined namelist control file.
!
! Method:
!   After setting up all the input and output files, the program loops
!   through each entry in the namelist file.  Each entry contains an
!   "Action" field, which is checked against the list of possible
!   actions using an IF construct.  Possible actions include reading
!   input fields into memory, using generic subroutines (e.g. Sum, Dif)
!   to manipulate fields already in memory, calculating diagnostics from
!   fields which may or may not be in memory and writing fields out to
!   file.
!   Various checks and constructs have been added to reduce IO and
!   memory overheads.  See online documentation for more details.
!
! Input Files:
!   Input files are defined using the following UNIX environment
!   variables:-
!   FIELDCALC_NL - Namelist file used to control the functions performed
!                  be Fieldcalc.  Essential for the running of the
!                  program.
!   UNIT20   - UM Fieldsfile containing standard level fields at the
!              required forecast ranges.  Essential for most
!              diagnostics.  Failure to open causes a fatal error.
!   UNIT21   - UM Fieldsfile containing model level fields at the
!              required forecast ranges.  Essential for most aviation
!              diagnostics.  Failure to open causes a fatal error.
!   UNIT22   - UM Fieldsfile containing standard level fields at
!              forecast ranges 6-12hrs earlier.  Failure to open will
!              result in a non-fatal error, but all accumulation
!              diagnostics will be unavailable.
!   OROGFILE - UM Fieldsfile containing the model orography field.
!              Failure to open will result in a non-fatal error, but
!              some aviation diagnostics will be unavailable.
!
! Output Files:
!   Output files are defined using the following UNIX environment
!   variables:-
!   UNIT30    - UM Fieldsfile.  Should contain only standard level flds.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO
USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  LSMField
USE FldCodes_mod, ONLY:   &
  ST_Urho,   ST_Vrho,     &
  ST_Htheta,              &
  ST_BlkCld, ST_ConCld,   &
  ST_CClBP,  ST_CClTP,    &
  ST_LCClBP, ST_LCClTP,   &
  ST_CPNRT,               &
  ST_Prho,   ST_Ptheta,   &
  ST_GWSU,   ST_GWSV,     &
  ST_Ttheta, ST_Pstar,    &
  ST_Ustd,   ST_Vstd,     &
  ST_Dust,   ST_TScreen,  &
  ST_VisIncPpn,           &
  ST_Dust1, ST_Dust2,     &
  ST_Dust3, ST_Dust4,     &
  ST_Dust5, ST_Dust6,     &
  ST_Orog,   LV_Surface,  &
  ST_Htheta, ST_Exnerrho, &
  ST_TropP,  ST_LSM,      &
  ProcessNo

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

! Maximum number of fields to be written
USE maxfields_mod, ONLY: MaxFldsOut

USE stashmaster_mod, ONLY: init_stm

USE ereport_mod, ONLY : ereport, ereport_finalise

USE filenamelength_mod, ONLY : & 
    filenamelength

USE UM_Config, ONLY : &
    appInit, &
    exe_fieldcalc

IMPLICIT None

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

! Local constants:

CHARACTER(LEN=*), PARAMETER :: ProgName = "FieldCalc"

INTEGER, PARAMETER :: MaxFlds    = 550
INTEGER, PARAMETER :: MaxWrite   = 20
INTEGER, PARAMETER :: NLUnit     = 81

INTEGER, PARAMETER :: StdPosn = 20          ! These constants mark
INTEGER, PARAMETER :: UPosn   = 40          ! certain areas of the
INTEGER, PARAMETER :: VPosn   = 95          ! Fields array for certain
INTEGER, PARAMETER :: PPosn   = 150         ! variables.  This makes it
INTEGER, PARAMETER :: TPosn   = 205         ! easier to reuse model
INTEGER, PARAMETER :: ZPosn   = 260         ! fields.

! Dust locations
INTEGER, PARAMETER :: DSTPosn = 315

! Internal field position for items required for ZTD
! The space required is the number of model levels (including ground level).
! Currently 71 is the maximum.
! This reuses the Dust locations since they should not be needed at the same
! time.
INTEGER, PARAMETER :: ZPosnRho  = 315
INTEGER, PARAMETER :: ExnerPosn = 390
INTEGER, PARAMETER :: QPosnTheta= 465

INTEGER :: UM_NumLevs  ! Number of model levels
INTEGER, PARAMETER :: Global = 0

! List of model levels that Fieldcalc will work for
INTEGER, PARAMETER :: Model_Levels_38 = 38
INTEGER, PARAMETER :: Model_Levels_50 = 50
INTEGER, PARAMETER :: Model_Levels_70 = 70

! Independent of number of model levels
INTEGER, PARAMETER :: MX_Zsea_Upr   = 20000 !   Upr lim if no.levs /= 38

! Variables for model level ranges required to derive diagnostics.

INTEGER :: TP_NumLevs_Gl = IMDI   ! TropHeight Fields:
INTEGER :: TP_ZeroLev_Gl = IMDI   !
INTEGER :: IT_NumLevs_Gl = IMDI   ! Isotherm Fields:
INTEGER :: IT_ZeroLev_Gl = IMDI   !
INTEGER :: CT_NumLevs_Gl = IMDI   ! Contrail Fields:
INTEGER :: CT_ZeroLev_Gl = IMDI   !
INTEGER :: MX_NumLevs_Gl = IMDI   ! MaxWind Fields:
INTEGER :: MX_ZeroLev_Gl = IMDI   !
INTEGER :: CA_NumLevs_Gl = IMDI   ! CAT Fields:
INTEGER :: CA_ZeroLev_Gl = IMDI   !
INTEGER :: MW_NumLevs_Gl = IMDI   ! MtnStress Fields:
INTEGER :: MW_ZeroLev_Gl = IMDI   !

INTEGER :: WT_NumLevs_Gl = IMDI   ! WAFC CAT turb:
INTEGER :: WT_ZeroLev_Gl = IMDI   !
INTEGER :: IC_NumLevs_Gl = IMDI   ! Icing Fields (original alg):
INTEGER :: IC_ZeroLev_Gl = IMDI   !
INTEGER :: LI_NumLevs_Gl = IMDI   ! Icing Fields:
INTEGER :: LI_ZeroLev_Gl = IMDI   !
INTEGER :: ICT_NumLevs_Gl = IMDI  ! In (Layer) cloud turb Fields:
INTEGER :: ICT_ZeroLev_Gl = IMDI  !
INTEGER :: CB_NumLevs_Gl = IMDI   ! CB Fields:
INTEGER :: CB_ZeroLev_Gl = IMDI   !
INTEGER :: DT_NumLevs_Gl = IMDI   ! Dust concs
INTEGER :: DT_ZeroLev_Gl = IMDI   ! 


! Number of standard levels required to derive diagnostics.

INTEGER, PARAMETER :: CA_NumStd     = 6   !  CAT Fields
INTEGER, PARAMETER :: MW_NumStd     = 2   !  MtnStress Fields

INTEGER, PARAMETER :: WT_NumStd     = 10  !  WAFC CAT Turb (Shear)

! Note that if the value of WTMW_NumStd is changed, the values set for
! MW_PRef must also be updated, since it's dimension is derived from
! WTMW_NumStd
INTEGER, PARAMETER :: WTMW_NumStd   = 8   !  WAFC CAT Turb (MW)

INTEGER, PARAMETER :: ICT_NumStd    = 5   !  In (Layer) Cloud Turb
INTEGER, PARAMETER :: CB_NumStd     = 1   !  CB Field


! Additional variables for WAFC_TURB action
! This is not a parameter since we may want to change it within the action.
INTEGER :: MW_n_PRef = WTMW_NumStd / 2 ! Number of pressure levels
                                       ! where MW predictor required

REAL, PARAMETER   :: MW_PRef(WTMW_NumStd / 2) =    &
                     (/ 30000.0, 25000.0, 20000.0, 15000.0 /)
                                                  ! Array of pressure levels where
                                                  !   MW predictor required (Pa)


! Local variables:

INTEGER :: P_Levels                    ! Number of model levels
INTEGER :: GridType                    ! Horizontal Grid Type
INTEGER :: NumFlds                     ! No. fields to perform action on
INTEGER :: NumLevs                     ! Number of Levels
INTEGER :: ZeroLev                     ! Level below lowest required
INTEGER :: i, j, k                     ! Loop counter
INTEGER :: ErrorStatus                 ! Program status monitor
INTEGER :: ReadStatus                  ! Namelist read status
REAL    :: Zsea                        ! Maximum height temporary
INTEGER :: n                           ! Index counter
CHARACTER(LEN=9)  :: Action     = "DUMMY"
CHARACTER(LEN=9)  :: PrevAction = "DUMMY"
CHARACTER(LEN=7)  :: LevType    = "STD"
CHARACTER(LEN=6)  :: PackType   = "NONE"
CHARACTER(LEN=80) :: ErrMessage = ""
CHARACTER(LEN=10) :: RunEnv     = ""
CHARACTER(LEN=filenamelength) :: NLFileName  = ""
CHARACTER(LEN=filenamelength) :: gribfile    = "unset"
CHARACTER(LEN=filenamelength) :: stdin_fname = "unset"
LOGICAL :: PPHdrMod = .FALSE.
TYPE(PP_Field_type)  :: OrogField      ! Orography PP-Field
TYPE(UM_Header_type) :: OrogHdr_in     ! UM Headers:  orography,
TYPE(UM_Header_type) :: StdLevHdr_in   !   standard level input
TYPE(UM_Header_type) :: ModLevHdr_in   !   model level input
TYPE(UM_Header_type) :: StdPrvHdr_in   !   previous std level input
TYPE(UM_Header_type) :: StdLevHdr_out  !   output

INTEGER :: Source  (MaxFlds)           ! Position in Fields array
INTEGER :: Store   (MaxFlds)           ! Position in Fields array
INTEGER :: STCode  (MaxFlds)           ! STASH code     (LBUSER(4))
INTEGER :: MO8Level(MaxFlds)           ! MetO8 Level    (LBLEV)
INTEGER :: FCTime  (MaxFlds)           ! Forecast Range (LBFT)
INTEGER :: LBProc  (MaxFlds)           ! Variable type  (LBPROC)
INTEGER :: USource (MaxFlds)           ! Model level source field
INTEGER :: VSource (MaxFlds)           !   postions in Fields array.
INTEGER :: PSource (MaxFlds)           !   Used in conjuction with UPosn
INTEGER :: TSource (MaxFlds)           !   etc. to minimise IO.
INTEGER :: ZSource (MaxFlds)
INTEGER :: MinsPastHr(MaxFlds)         ! Minutes past whole hour for COPYFLDS
INTEGER :: DSource (MaxFlds)           ! Dust concentration
INTEGER :: MinsDif (MaxFlds)           ! Difference in minutes in headers.

! Specify a range of criteria to match.
INTEGER :: MO8LevelLwr  (MaxFlds)
INTEGER :: MO8LevelUpr  (MaxFlds)
INTEGER :: FCTimeStart  (MaxFlds)
INTEGER :: FCTimeFreq   (MaxFlds)
INTEGER :: FCTimeEnd    (MaxFlds)

! Items required for ZTD
INTEGER :: ZSourceRho  (MaxFlds)
INTEGER :: ExnerSource (MaxFlds)
INTEGER :: QSourceTheta(MaxFlds)

REAL    :: Factor  (MaxFlds)           ! Input for diagnostics
REAL    :: PackAcc (MaxFlds)           ! Packing Accuracy

! Additional local variables for icing actions
INTEGER :: NumLyrs                     ! No. layers to output fields for
INTEGER :: LyrMO8L (MaxFlds)           ! LBLEV values assigned to layers
REAL    :: LyrLwrB (MaxFlds)           ! Layer lower boundaris (BRLEV)
REAL    :: LyrUprB (MaxFlds)           ! Layer upper boundaries (BLEV)

! Additional variables for tropopause height
INTEGER :: TropSTCode

! Additional variables for SST perturbations
INTEGER :: ens_member

INTEGER :: me_gc
INTEGER :: nproc_gc

TYPE(PP_Header_type) :: NewPPHdr(MaxWrite)
TYPE(PP_Field_type)  :: Fields(MaxFlds)! Input and output fields
TYPE(PP_Field_type)  :: TempFields(MaxFlds)! Input and output fields

Character (Len=80)   :: Cmessage = ' '   ! Error Message
Character (Len=8)    :: c_ens_member     ! Char string for Ensemble member
CHARACTER (LEN=*), PARAMETER :: routinename = 'FIELDCALC'
NAMELIST / PPHdrModNL / &      ! Header modification flag namelist
  PPHdrMod
NAMELIST / CommandNL /  &      ! Command namelist - for program control
  Action,               &
  NumFlds,              &
  Source,               &
  Store,                &
  STCode,               &
  MO8Level,             &
  FCTime,               &
  LBProc,               &
  LevType,              &
  PackType,             &
  PackAcc,              &
  MinsPastHr,           &
  MinsDif,              &
  Factor,               &
  NumLyrs,              &
  LyrMO8L,              &
  LyrLwrB,              &
  LyrUprB,              &
  MO8LevelLwr,          &
  MO8LevelUpr,          &
  FCTimeStart,          &
  FCTimeFreq,           &
  FCTimeEnd


NAMELIST / NewPPHdrNL / &      ! PP Header namelist - for writing fields
  NewPPHdr

! End of Header --------------------------------------------------------

! Initialise errorstatus
errorstatus = statusOK

! DEPENDS ON: timer
CALL Timer( ProgName, 1 )

! DEPENDS ON: initprintstatus
CALL InitPrintStatus
CALL gc_init(' ',me_gc,nproc_gc)
CALL appInit(exe_fieldcalc)
CALL ioInit()

!-----------------------------------------------------------------------
! 1 PROGRAM SETUP
!-----------------------------------------------------------------------

DO i=1,MaxFlds
  NULLIFY( Fields(i) % RData )
  NULLIFY( TempFields(i) % RData )
  Fields(i)     % ArrayPos  = i
  TempFields(i) % ArrayPos  = i
  Fields(i)     % LookupPos = 0
  TempFields(i) % LookupPos = 0
END DO
NULLIFY( OrogField % RData )
OrogField % ArrayPos = 9999
NULLIFY( LSMField % RData )
LSMField % ArrayPos = 9999

!--------------------------------------------
! 1.1.1 Initialise GRIB file from Env Var
!--------------------------------------------

Call FORT_GET_ENV('GRIBFILE',8, gribfile,filenamelength, ErrorStatus)

If ( ErrorStatus /=  StatusOk) THEN
  Write(Cmessage,'(A)') 'GRIBFILE has not been set.  Cannot use GRIB actions'
  ErrorStatus = StatusWarning
  Call Ereport( ProgName, ErrorStatus, Cmessage )
  gribfile='unset'
End If

!-------------------------
! 1.2 Open Namelist File
!-------------------------
CALL FORT_GET_ENV( "FIELDCALC_NL", 12, NLFileName, filenamelength, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
  ErrorStatus = StatusFatal    ! Having no namelist is a fatal error

  CALL EReport( ProgName, ErrorStatus, &   ! Cannot run without namelist
                "Cannot read env FIELDCALC_NL" )
END IF
WRITE(6,'(2A)') "Namelist File : ", TRIM(NLFileName)
OPEN( UNIT   = NLUnit,        &
      ACCESS = "SEQUENTIAL",  &
      ACTION = "READ",        &
      FILE   = NLFileName,    &
      FORM   = "FORMATTED",   &
      IOSTAT = ErrorStatus,   &
      STATUS = "OLD" )
IF ( ErrorStatus /= StatusOK ) THEN
  ErrorStatus = StatusFatal

  CALL EReport( ProgName, ErrorStatus, &   ! Cannot run without namelist
               "Cannot open Namelist file" )
END IF

!-----------------------------
! 1.3 Initialise Input Files
!-----------------------------
StdLevHdr_in % UnitNum = 20
StdLevHdr_in % FileNameEnv = "UNIT20"
! Make sure this is null since we check later on.  Possibly not needed if we
! initialise all pointers to NULL.
NULLIFY(StdLevHdr_in % lookup)
ModLevHdr_in % UnitNum = 21
ModLevHdr_in % FileNameEnv = "UNIT21"
StdPrvHdr_in % UnitNum = 22
StdPrvHdr_in % FileNameEnv = "UNIT22"


! We can only process GRIBFILE or Fieldsfile but not both.
Call FORT_GET_ENV(StdLevHdr_in % FileNameEnv,6, stdin_fname,&
                  filenamelength, ErrorStatus)

! We dont have the environment variable so we must be processing GRIB input
! only.
IF (ErrorStatus == StatusOK) THEN
  
! DEPENDS ON: read_umhdr
  CALL Read_UMHdr( StdLevHdr_in, .FALSE., ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning
    CALL EReport( ProgName, ErrorStatus,     &    ! error - no input
                 "Error Reading Standard Level Header" )
  END IF
  
! DEPENDS ON: read_umhdr
  CALL Read_UMHdr( ModLevHdr_in, .FALSE., ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning
    CALL EReport( ProgName, ErrorStatus,     &    ! error - no input
                 "Error Reading Model Level Header" )
  ! Get number of model levels from Integer Constants, Word 8
    P_Levels = IMDI
  
  ! Get grid type from Fixed Header, Word 4
    GridType = IMDI
  ELSE
  ! Get number of model levels from Integer Constants, Word 8
    P_Levels = ModLevHdr_in % IntC(8)
  
  ! Get grid type from Fixed Header, Word 4
    GridType = ModLevHdr_in % FixHd(4)
  END IF
  
! DEPENDS ON: read_umhdr
  CALL Read_UMHdr( StdPrvHdr_in, .FALSE., ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning
    CALL EReport( ProgName, ErrorStatus,     &    ! Non-fatal error
                 "Error Reading Previous Standard Level Header" )
    ErrorStatus = StatusOK
    NULLIFY( StdPrvHdr_in % Lookup )
  END IF
  
  !---------------------------
  ! 1.3.1 Read Orography Field
  !---------------------------
  OrogHdr_in % UnitNum = 23
  OrogHdr_in % FileNameEnv = "OROGFILE"
! DEPENDS ON: read_umhdr
  CALL Read_UMHdr( OrogHdr_in, .FALSE., ErrorStatus )
  IF ( ErrorStatus == StatusOK ) THEN
    NumFlds       =       1
    MO8Level  (1) = LV_Surface
    FCTime    (1) =       0
    LBProc    (1) =      -1
    MinsPastHr(1) =       0
    MinsDif   (1) =      -1
    Store     (1) =       1
 
  ! The argument for MinsPastHr has been set to 0 manually because the
  ! namelist hasn't been read yet.
  
  ! First get orography
    STCode    (1) = ST_Orog
! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTime,                  &
                  LBProc, MinsPastHr, MinsDif, Store, PPHdrMod, OrogHdr_in,    &
                  OrogField, ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      ErrorStatus = StatusWarning
      CALL EReport( ProgName, ErrorStatus,                 &
                    "Error Reading Orography Field - &
                    &some diagnostics will be unavailable" )
      NULLIFY( OrogField % RData )
    END IF
  
  ! Now get land-sea mask
    STCode (1) = ST_LSM
! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTIME,                  &
                  LBProc, MinsPastHr, MinsDif, Store, PPHdrMod, OrogHdr_In,    &
                  LSMField, ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      ErrorStatus = StatusWarning
      CALL EReport( ProgName, ErrorStatus,                 &
                   "Error Reading LSM Field - &
                   &land-packed fields will not be unpacked." )
      NULLIFY( LSMField % RData )
    END IF
  END IF
  ErrorStatus = StatusOk
  CALL File_Close ( OrogHdr_in % UnitNum,                &
                    OrogHdr_in % FileNameEnv,            &
                    LEN_TRIM(OrogHdr_in % FileNameEnv),  &
                    0, 0, ErrorStatus )
  
  
  !---------------------------------------------------------------
  ! 1.4.1 Initialise model level ranges now that P_Levels is known
  !---------------------------------------------------------------
  
  Select Case ( P_Levels )
  
    Case ( Model_Levels_38 )
  
      UM_NumLevs    = 38    ! 38 Model Levels
  
      TP_NumLevs_Gl = 27    ! TropHeight Fields:
      TP_ZeroLev_Gl = 5     !   Model Levels 6-32
      IT_NumLevs_Gl = 32    ! Isotherm Fields:
      IT_ZeroLev_Gl = 0     !   Model Levels 1-32
      CT_NumLevs_Gl = 30    ! Contrail Fields:
      CT_ZeroLev_Gl = 0     !   Model Levels 1-30
      MX_NumLevs_Gl = 32    ! MaxWind Fields:
      MX_ZeroLev_Gl = 0     !   Model Levels 1-32
      CA_NumLevs_Gl = 24    ! CAT Fields:
      CA_ZeroLev_Gl = 5     !   Model Levels 6-29
      MW_NumLevs_Gl = 15    ! MtnStress Fields:
      MW_ZeroLev_Gl = 13    !   Model Levels 14-28
  
      WT_NumLevs_Gl = 32    ! WAFC CAT turb ...
      WT_ZeroLev_Gl = 0     !   Model Levels 1-32
      IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
      IC_ZeroLev_Gl = 0     !   Model Levels 1-24
      LI_NumLevs_Gl = 33    ! Icing Fields:
      LI_ZeroLev_Gl = 0     !   Model Levels 1-33
      ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
      ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
      CB_NumLevs_Gl = 32    ! CB Fields:
      CB_ZeroLev_Gl = 0     !   Model Levels 1-32
  
      DT_NumLevs_Gl = 10    ! Dust concs
      DT_ZeroLev_Gl = 0     !   Model Levels 1-10
  
  
    Case ( Model_Levels_50 )
  
      UM_NumLevs    = 50    ! 50 Model Levels
  
      TP_NumLevs_Gl = 28    ! TropHeight Fields:
      TP_ZeroLev_Gl = 5     !   Model Levels 6-33
      IT_NumLevs_Gl = 33    ! Isotherm Fields:
      IT_ZeroLev_Gl = 0     !   Model Levels 1-33
      CT_NumLevs_Gl = 30    ! Contrail Fields:
      CT_ZeroLev_Gl = 0     !   Model Levels 1-30
      MX_NumLevs_Gl = 33    ! MaxWind Fields:
      MX_ZeroLev_Gl = 0     !   Model Levels 1-33
      CA_NumLevs_Gl = 24    ! CAT Fields:
      CA_ZeroLev_Gl = 5     !   Model Levels 6-29
      MW_NumLevs_Gl = 15    ! MtnStress Fields:
      MW_ZeroLev_Gl = 13    !   Model Levels 14-28
  
      WT_NumLevs_Gl = 33    ! WAFC CAT turb ...
      WT_ZeroLev_Gl = 0     !   Model Levels 1-33
      IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
      IC_ZeroLev_Gl = 0     !   Model Levels 1-24
      LI_NumLevs_Gl = 33    ! Icing Fields:
      LI_ZeroLev_Gl = 0     !   Model Levels 1-33
      ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
      ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
      CB_NumLevs_Gl = 33    ! CB Fields:
      CB_ZeroLev_Gl = 0     !   Model Levels 1-33
  
      DT_NumLevs_Gl = 10    ! Dust concs
      DT_ZeroLev_Gl = 0     !   Model Levels 1-10
  
  
    Case ( Model_Levels_70 )
  
      UM_NumLevs    = 70    ! 70 Model Levels
  
      TP_NumLevs_Gl = 46    ! TropHeight Fields:
      TP_ZeroLev_Gl = 8     !   Model Levels 9-54
      IT_NumLevs_Gl = 54    ! Isotherm Fields:
      IT_ZeroLev_Gl = 0     !   Model Levels 1-54
      CT_NumLevs_Gl = 50    ! Contrail Fields:
      CT_ZeroLev_Gl = 0     !   Model Levels 1-50
      MX_NumLevs_Gl = 54    ! MaxWind Fields:
      MX_ZeroLev_Gl = 0     !   Model Levels 1-54
      CA_NumLevs_Gl = 41    ! CAT Fields:
      CA_ZeroLev_Gl = 8     !   Model Levels 9-49
      MW_NumLevs_Gl = 25    ! MtnStress Fields:
      MW_ZeroLev_Gl = 22    !   Model Levels 23-47
      WT_NumLevs_Gl = 54    ! WAFC CAT turb ...
  
      WT_ZeroLev_Gl = 0     !   Model Levels 1-54
      IC_NumLevs_Gl = 41    ! Icing Fields:
      IC_ZeroLev_Gl = 0     !   Model Levels 1-41
      LI_NumLevs_Gl = 54    ! Icing Fields:
      LI_ZeroLev_Gl = 0     !   Model Levels 1-54
      ICT_NumLevs_Gl= 41    ! In (Layer) cloud turb Fields:
      ICT_ZeroLev_Gl= 0     !   Model Levels 1-41
      CB_NumLevs_Gl = 54    ! CB Fields:
      CB_ZeroLev_Gl = 0     !   Model Levels 1-54
  
      DT_NumLevs_Gl = 16    ! Dust concs
      DT_ZeroLev_Gl = 0     !   Model Levels 1-16
  
    Case (IMDI)
      UM_NumLevs = IMDI
      ErrorStatus = StatusWarning
      CALL EReport( ProgName, ErrorStatus,    &  ! error - no output
                   "Number of levels not catered "//&
                   "for but trying to continue." )
  
    Case Default
      UM_NumLevs    = -1    ! Unknown Model Levels
  
      TP_NumLevs_Gl = 0    ! TropHeight Fields:
      TP_ZeroLev_Gl = 0    !   Model Levels
      IT_NumLevs_Gl = 0    ! Isotherm Fields:
      IT_ZeroLev_Gl = 0    !   Model Levels
      CT_NumLevs_Gl = 0    ! Contrail Fields:
      CT_ZeroLev_Gl = 0    !   Model Levels
      MX_NumLevs_Gl = 0    ! MaxWind Fields:
      MX_ZeroLev_Gl = 0    !   Model Levels
      CA_NumLevs_Gl = 0    ! CAT Fields:
      CA_ZeroLev_Gl = 0    !   Model Levels
      MW_NumLevs_Gl = 0    ! MtnStress Fields:
      MW_ZeroLev_Gl = 0    !   Model Levels
      WT_NumLevs_Gl = 0    ! WAFC CAT turb ...
  
      WT_ZeroLev_Gl = 0    !   Model Levels
      IC_NumLevs_Gl = 0    ! Icing Fields:
      IC_ZeroLev_Gl = 0    !   Model Levels
      LI_NumLevs_Gl = 0    ! Icing Fields:
      LI_ZeroLev_Gl = 0    !   Model Levels
      ICT_NumLevs_Gl= 0    ! In (Layer) cloud turb Fields:
      ICT_ZeroLev_Gl= 0    !   Model Levels
      CB_NumLevs_Gl = 0    ! CB Fields:
      CB_ZeroLev_Gl = 0    !   Model Levels
  
      DT_NumLevs_Gl = 0    ! Dust concs
      DT_ZeroLev_Gl = 0    !   Model Levels
  
      ! Warn about unknown model levels.
      ErrorStatus = StatusWarning
      CALL EReport( ProgName, ErrorStatus,    &
                    "Number of model levels not catered for." )
  
  End Select
  
  Write (6,'(A)')      ' '
  Write (6,'(A,I7,A)') ' FieldCalc has been set up for a ',UM_NumLevs,    &
                       ' level model.'
  Write (6,'(A)')      ' '
ELSE
  Write (6,'(A)') ' '
  Write (6,'(A)') ' FieldCalc has been set up only to process input GRIB files.'
  Write (6,'(A)') ' '
  errorstatus = statusOK
END IF


!-----------------------------
! 1.5.1 Initialise output file
!-----------------------------
StdLevHdr_out % UnitNum = 30
StdLevHdr_out % FileNameEnv = "UNIT30"

! Initialise with existing file only if available.
IF (ASSOCIATED(StdLevHdr_in % Lookup) .AND. gribfile == 'unset') THEN
! DEPENDS ON: new_umhdr
  CALL New_UMHdr( StdLevHdr_in, MaxFldsOut, StdLevHdr_out, ErrorStatus )
END IF

IF ( ErrorStatus /= StatusOK ) THEN
  ErrorStatus = StatusFatal

  CALL EReport( ProgName, ErrorStatus,    &    ! Fatal error - no output
               "Error Setting up Output File." )
END IF

!--------------------------
! 1.6 Get other variables
!--------------------------
READ( UNIT = NLUnit,       &
       NML = PPHdrModNL,   &
    IOSTAT = ReadStatus )
IF ( ReadStatus /= StatusOK ) THEN

  CALL EReport( ProgName, ReadStatus, &
               "Error found in PPHdrModNL" )
END IF
WRITE(6,'(A,L7)') "Header modification = ", PPHdrMod
WRITE(6,'(A,L7)') "Global model = ", (GridType == Global)

!--------------------------
! 1.7 Read in STASHmaster
!--------------------------
CALL init_stm

!-----------------------------------------------------------------------
! 2 MAIN DIAGNOSTIC PROCESSING
!-----------------------------------------------------------------------

DO WHILE ( Action /= "END" )

 !----------------------------
 ! 2.1 Read Command Namelist
 !----------------------------
  NumFlds   = 0
  Source    = Store
  Store(:)  = 0
  LBProc(:) = -1
  FCTime(:) = -1
  MinsPastHr(:) = 00
  MinsDif(:)    = -1
  Factor(:) = 0.0
  LevType   = "STD"
  PackType  = "NONE"
  MO8Level    (:) = -1
  MO8LevelLwr (:) = -1
  MO8LevelUpr (:) = -1
  FCTimeStart (:) = -1
  FCTimeFreq  (:) =  1
  FCTimeEnd   (:) = -1
 
  READ( UNIT = NLUnit,       &
         NML = CommandNL,    &
      IOSTAT = ReadStatus )
  IF      ( ReadStatus == EndofFile ) THEN
    WRITE(6,'(A)') "Unexpectedly reached end of Namelist file - quitting"
    GO TO 9999
  ELSE IF ( ReadStatus /= StatusOK ) THEN
    ErrorStatus = ReadStatus

    CALL EReport( ProgName, ReadStatus, &
                 "Error found in CommandNL" )
  END IF

  WRITE(6,'(A8,A9,A2,I3,A1)') "Action: ", Action, " (", NumFlds, ")"
  
  ! gribfile should only be a file for non-GRIB actions.
  IF (gribfile /= 'unset') THEN
    IF ( Action /= "DEGRIBIFY" .AND. &
         Action /= "GRIBIFY"   .AND. &
         Action /= "END" ) THEN
      Errorstatus = StatusFatal
      CALL ereport(routinename, errorstatus, &
                  'GRIBFILE has been set but using non-GRIB action.')
    END IF
  END IF

  ! Actions which depend on orography
  SELECT CASE(Action)
  CASE ("TROPHGHT","ISOTHERM","CONTRAIL","CATURB","WAFC_TURB",         &
        "CLD_TURB","DUST_CONC","ZTD")
    IF ( .NOT. ASSOCIATED( OrogField % RData ) ) THEN
      ErrorStatus = StatusWarning

      CALL EReport( ProgName, ErrorStatus,                      &
                   "Orography not available - cannot calculate action." )
      ErrorStatus = StatusWarning
    END IF
  END SELECT
  IF ( Action == "PRSYM" ) THEN
    NumFlds = 5
  END IF
  ! Actions which cannot store to the location as an input
  SELECT CASE(Action)
  CASE ("WINDSPEED", "DIVERG", "VORTIC", "THERMADV", "SNOWPROB", "PRSYM")
    DO i = 1, NumFlds
      IF (Store(1) == Source(i)) THEN
        ErrorStatus = StatusWarning
        CALL EReport( ProgName, ErrorStatus,                      &
                   "Cannot store in a source location." )
        ErrorStatus = StatusWarning
      END IF
    END DO
  CASE ("ICAO_HT")
    DO i = 1, NumFlds
      IF (ANY(Store(1:NumFlds) == Source(i))) THEN
        ErrorStatus = StatusWarning
        CALL EReport( ProgName, ErrorStatus,                      &
                   "Cannot store in a source location." )
        ErrorStatus = StatusWarning
      END IF
    END DO
  CASE ("WXCODE")
    DO i = 1, 12
      IF (ANY(Store(1:2) == Source(i))) THEN
        ErrorStatus = StatusWarning
        CALL EReport( ProgName, ErrorStatus,                      &
                   "Cannot store in a source location." )
        ErrorStatus = StatusWarning
      END IF
    END DO
  CASE ("TROPHGHT", "ISOTHERM", "CONTRAIL", "MAXWIND", "TOPBASE", "CATURB", &
        "MWTURB", "WAFC_TURB", "ICING", "CLD_TURB", "CB_ACT", "DUST_CONC",  &
        "ZTD")
    IF (ANY(Store(:) > MaxWrite)) THEN
      ErrorStatus = StatusWarning
      CALL EReport( ProgName, ErrorStatus,                      &
                   "Cannot store in a location greater than 20" )
      ErrorStatus = StatusWarning
    END IF
  END SELECT

  ! Dependencies on levtype
  SELECT CASE(levtype)
  CASE ("STDPREV")
    IF ( .NOT. ASSOCIATED( StdPrvHdr_In % Lookup ) ) THEN
      ErrorStatus = StatusWarning

      CALL EReport( ProgName, ErrorStatus,                      &
                   "STDPREV levtype not available - cannot calculate action." )
      ErrorStatus = StatusWarning
    END IF

  END SELECT
 
! Only make actions available if input file type is.
  IF (Levtype == "STD" .AND. Action /= "END") THEN
    IF (.NOT. ASSOCIATED(StdLevHdr_in % lookup)) THEN
      WRITE(6,'(A)') "Command ignored due to filetype unavailable"
      Errorstatus = StatusWarning
    END IF
  ELSE IF (Levtype == "MOD") THEN
    IF (.NOT. ASSOCIATED(ModLevHdr_in % lookup)) THEN
      WRITE(6,'(A)') "Command ignored due to filetype unavailable"
      Errorstatus = StatusWarning
    END IF
  ELSE IF (Levtype == "STDPREV") THEN
    IF (.NOT. ASSOCIATED(StdPrvHdr_in % lookup)) THEN
      WRITE(6,'(A)') "Command ignored due to filetype unavailable"
      Errorstatus = StatusWarning
    END IF
  ELSE IF (Levtype == "GRIB") THEN
    IF (gribfile == 'unset') THEN
      WRITE(6,'(A)') "Command ignored due to filetype unavailable"
      Errorstatus = StatusWarning
    END IF
  END IF
       
  IF      ( Action == "CLEAR_ERR" ) THEN
    ! 2.1.1 Clear previous errors
    ErrorStatus = StatusOK
    ErrMessage = ""

  ELSE IF ( Action == "END" )       THEN
   ! 2.1.2 End of namelist
    ErrorStatus = StatusOK
    WRITE(6,'(A)') "All actions completed successfully."

  ELSE IF ( ErrorStatus /= StatusOK ) THEN
   ! 2.1.3 Previous Error - skip actions until error cleared
    WRITE(6,'(A)') "Command ignored due to previous error."
    IF ( Action == "WRITEFLDS" ) THEN
      NewPPHdr(1:NumFlds) % LBUser7 = ProcessNo
      ! Read the NewPPHdrNL that will come up after WRITEFLDS
      READ( UNIT = NLUnit,      &
            NML  = NewPPHdrNL )     ! Don't worry about errors
    END IF

 !---------------------------------
 ! 2.2 Actions for performing I/O
 !---------------------------------
  ELSE IF ( Action == "COPYALL" )   THEN

   ! 2.2.1 Rewrite entire contents of input file without modification
   !       Save & unpack specified fields
   !       Scale specified fields before writing
    IF      ( LevType == "STD" ) THEN
! DEPENDS ON: copyall
      CALL CopyAll( NumFlds, MaxFlds, STCode,                              &
                    MO8Level,  MO8LevelLwr,  MO8LevelUpr,                  &
                    FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,            &
                    LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,  &
                    StdLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE IF ( LevType == "MOD" ) THEN
! DEPENDS ON: copyall
      CALL CopyAll( NumFlds, MaxFlds, STCode,                              &
                    MO8Level,  MO8LevelLwr,  MO8LevelUpr,                  &
                    FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,            &
                    LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,  &
                    ModLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE
      ErrorStatus = StatusWarning

      CALL EReport( ProgName, ErrorStatus, &
                    "COPYALL : Unknown level type " // LevType )
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "COPYFLDS" )  THEN

   ! 2.2.2 Rewrite specified fields
   !       Save & unpack those with Store /= 0
    IF      ( LevType == "STD" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode,                              &
                      MO8Level,  MO8LevelLwr,  MO8LevelUpr,                  &
                      FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,            &
                      LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,  &
                      StdLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE IF ( LevType == "MOD" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode,                              &
                      MO8Level,  MO8LevelLwr,  MO8LevelUpr,                  &
                      FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,            &
                      LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,  &
                      ModLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE IF ( LevType == "STDPREV" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode,                              &
                      MO8Level,  MO8LevelLwr,  MO8LevelUpr,                  &
                      FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,            &
                      LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,  &
                      StdPrvHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE
      ErrorStatus = StatusWarning

      CALL EReport( ProgName, ErrorStatus, &
                    "COPYFLDS : Unknown level type " // LevType )
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "GETFLDS" )   THEN

   ! 2.2.3 Read fields
    IF      ( LevType == "STD" )     THEN
! DEPENDS ON: getflds
      CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                     LBProc, MinsPastHr, MinsDif, Store, PPHdrMod,      &
                     StdLevHdr_in, Fields, ErrorStatus )
    ELSE IF ( LevType == "MOD" )     THEN
! DEPENDS ON: getflds
      CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                     LBProc, MinsPastHr, MinsDif, Store, PPHdrMod,      &
                     ModLevHdr_in,Fields, ErrorStatus )
    ELSE IF ( LevType == "STDPREV" ) THEN
      IF ( ASSOCIATED( StdPrvHdr_in % Lookup ) ) THEN
! DEPENDS ON: getflds
        CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                       LBProc, MinsPastHr, MinsDif, Store, PPHdrMod,      &
                       StdPrvHdr_in, Fields, ErrorStatus )
      ELSE
        WRITE(6,'(A)') "GETFLDS : Cannot read Previous Standard Level file"
        ErrorStatus = StatusWarning
      END IF
    ELSE
      WRITE(6,'(A)') "GETFLDS : Unknown level type ", LevType
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "PACK" )      THEN

   ! 2.2.4 Pack fields
! DEPENDS ON: packflds
    CALL PackFlds( NumFlds, MaxFlds, PackType, Source, PackAcc, &
                   Fields, ErrorStatus )

  ELSE IF ( Action == "WRITEFLDS" ) THEN

   ! 2.2.5 Write fields

    NewPPHdr(1:NumFlds) = Fields(Source(1:NumFlds)) % Hdr
! Set LBuser7 to default value (ProcessNo) before reading namelist
    NewPPHdr(1:NumFlds) % LBUser7 = ProcessNo
    READ( UNIT = NLUnit,     &
           NML = NewPPHdrNL, &
        IOSTAT = ReadStatus )
    IF ( ReadStatus /= StatusOK ) THEN
      ErrorStatus = ReadStatus

      CALL EReport( ProgName, ReadStatus, "Error found in NewPPHdrNL" )
    END IF

    Fields(Source(1:NumFlds)) % Hdr = NewPPHdr(1:NumFlds)

! DEPENDS ON: writeflds
    CALL WriteFlds( NumFlds, MaxFlds, MaxFldsOut, PackType, Source, &
                    PackAcc, Fields, StdLevHdr_out, ErrorStatus )

 !---------------------------------
 ! 2.3 General Arithmetic Actions
 !---------------------------------
  ELSE IF ( Action == "SUM" )       THEN

   ! 2.3.1 Sum two or more fields
! DEPENDS ON: sum
    CALL Sum( Fields(Source(1)), Fields(Source(2)),     &
              Fields(Store(1)), ErrorStatus )
    DO i = 3,NumFlds
! DEPENDS ON: sum
      CALL Sum( Fields(Source(i)), Fields(Store(1)),&
                Fields(Store(1)), ErrorStatus )
    END DO

  ELSE IF ( Action == "DIF" )       THEN

   ! 2.3.2 Take difference of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,'(A)') "Must have 2 fields for differencing"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: dif
      CALL Dif( Fields(Source(1)), Fields(Source(2)),   &
                Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "ADD" )       THEN

   ! 2.3.3 Add constant to fields
    DO i = 1,NumFlds
! DEPENDS ON: add
      CALL Add( Factor(1), Fields(Source(i)),   &
                Fields(Store(i)), ErrorStatus )
    END DO
  ELSE IF ( Action == "MAXIMUM" )       THEN
! DEPENDS ON: maximum
    CALL maximum( Fields(Source(1)), Fields(Source(2)),   &
                  Fields(Store(1)), ErrorStatus )
   ! Perform operation on remaining fields
    DO i = 3,NumFlds
! DEPENDS ON: add
      CALL maximum( Fields(Source(i)), Fields(Store(1)),   &
                    Fields(Store(1)), ErrorStatus )
    END DO
  ELSE IF ( Action == "MINIMUM" )       THEN
! DEPENDS ON: minimum
    CALL minimum( Fields(Source(1)), Fields(Source(2)),   &
                  Fields(Store(1)), ErrorStatus )
   ! Perform operation on remaining fields.
    DO i = 3,NumFlds
! DEPENDS ON: minimum
      CALL minimum( Fields(Source(i)), Fields(Store(1)),   &
                    Fields(Store(1)), ErrorStatus )
    END DO

  ELSE IF ( Action == "SCALE" )     THEN

   ! 2.3.4 Multipy fields by a constant
    DO i = 1,NumFlds
! DEPENDS ON: scale
      CALL Scale( Factor(1), Fields(Source(i)),   &
                  Fields(Store (i)), ErrorStatus )
    END DO
  ELSE IF ( Action == "WINDSPEED" ) THEN

   ! 2.3.5 Calculate magnitude of a vector (2 flds) (e.g. - wind speed)
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,'(A)') "Must have a U and a V field for WindSpeed"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: vecmag
      CALL VecMag( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "MULTIPLY" )     THEN
    ! 2.3.6 multiply scalar fields
    IF ( NumFlds < 2 ) THEN
      WRITE(6,'(A)') "Must have at least 2 fields for multiplying"
      ErrorStatus = StatusWarning
    ELSE

! DEPENDS ON: multiply
      CALL Multiply( 1, Fields(Source(1)), Fields(Source(2)),     &
                     Fields(Store(1)), ErrorStatus )
      DO i = 3,NumFlds
! DEPENDS ON: multiply
        CALL Multiply( 1, Fields(Source(i)), Fields(Store(1)),    &
                       Fields(Store(1)), ErrorStatus )
      END DO
    END IF


 !---------------------------------------------------
 ! 2.4 Diagnostic Actions Requiring Standard Levels
 !---------------------------------------------------
  ELSE IF ( Action == "DIVERG" )    THEN

   ! 2.4.1 Calculate the Divergence of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,'(A)') "Must have a U and a V field for Diverg"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: diverg
      CALL Diverg( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "VORTIC" )    THEN

   ! 2.4.2 Calculate the Vorticity of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,'(A)') "Must have a U and a V field for Vortic"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: vortic
      CALL Vortic( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "THERMADV" )  THEN

   ! 2.4.3 Thermal Advection
    IF ( NumFlds /= 3 ) THEN
      WRITE(6,'(A)') "Must have U, V and T fields for Thermal Advection"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: thermadv
      CALL ThermAdv( Fields(Source(1)), Fields(Source(2)),   &
                     Fields(Source(3)),                      &
                     Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "ICAO_HT" )   THEN

   ! 2.4.4 Convert pressure to ICAO height
    DO i = 1,NumFlds
! DEPENDS ON: icaoheight
      CALL ICAOHeight( Fields(Source(i)),    &
                       Fields(Store(i)), ErrorStatus )
    END DO

  ELSE IF ( Action == "SNOWPROB" )  THEN

   ! 2.4.5 Snow Probability
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,'(A)') "Must have two height fields for Snow Probability"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: snowprob
      CALL SnowProb( Fields(Source(1)), Fields(Source(2)),  &
                     Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "PRSYM" )     THEN

   ! 2.4.6 Total Precipitation Code
! DEPENDS ON: prsym
    CALL PrSym( Fields(Source(1)), Fields(Source(2)),     &
                Fields(Source(3)), Fields(Source(4)),     &
                Fields(Source(5)),                        &
                Fields(Store(1)), ErrorStatus )

  ELSE IF ( Action == "WXCODE" )    THEN

   ! 2.4.7 Present Weather Code
! DEPENDS ON: wxcode
    CALL WXCode( Fields(Source( 1)), Fields(Source( 2)),  &
                 Fields(Source( 3)), Fields(Source( 4)),  &
                 Fields(Source( 5)), Fields(Source( 6)),  &
                 Fields(Source( 7)), Fields(Source( 8)),  &
                 Fields(Source( 9)), Fields(Source(10)),  &
                 Fields(Source(11)), Fields(Source(12)),  &
                 Fields(Store(1)), ErrorStatus )
    IF ( NumFlds == 2 ) THEN
! DEPENDS ON: prsym
      CALL PrSym( Fields(Source(1)), Fields(Source(2)),   &
                  Fields(Source(3)), Fields(Source(4)),   &
                  Fields(Source(5)),                      &
                  Fields(Store(2)), ErrorStatus )
    END IF

  ELSE IF ( Action == "SSTPERT" )  THEN

  ! 2.4.8 Generate random SST perturbation
  !       Random seed set using model date and ENS_MEMBER
    IF ( NumFlds /= 13 ) THEN
      WRITE(6,'(A)') "Must have 13 fields for SST perturbations"
      WRITE(6,'(A)') " ** Initial T* field"
      WRITE(6,'(A)') " ** 12 monthly DeltaSST fields"
      ErrorStatus = StatusWarning
    ELSE

      ! 2.4.8.1 Initialise random seed: Read ENS_MEMBER from Env Var
      CALL FORT_GET_ENV('ENS_MEMBER',10, c_ens_member, 8, ErrorStatus)
      IF ( ErrorStatus /=  0) THEN
        WRITE(Cmessage,'(A)') 'ENS_MEMBER has not been set: Using = 1'
        ErrorStatus = StatusWarning
        CALL ereport( ProgName, ErrorStatus, Cmessage )
        ens_member = 1
      ELSE
        READ( c_ens_member, '(I8)') ens_member
      END IF
      IF ( Factor(1) <=  0) THEN
        WRITE(Cmessage,'(A)') 'SST perturbation alpha factor <= 0'
        ErrorStatus = StatusFatal
        CALL ereport( ProgName, ErrorStatus, Cmessage )
      END IF

! DEPENDS ON: sstpert
      CALL sstpert( Factor(1), Fields(Source(1)), ens_member,           &
                    Fields(Source(2):Source(13)),                       &
                    Fields(Store(1)), ErrorStatus )
    END IF


 !------------------------------------------------
 ! 2.5 Diagnostic Actions Requiring Model Levels
 !------------------------------------------------
  ELSE IF ( Action == "TROPHGHT" )  THEN

   ! 2.5.1 Tropopause Temperature, Pressure and Height
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run TROPHGHT on a 38, 50 or 70 level model"
      ErrorStatus = StatusWarning
    END IF

    ! For GRIB2 a new max tropopause is needed and was given a new
    ! STASH code.  The value of the maximum trop. is given in the
    ! factor array and a new STCode supplied.  Defaults to old GRIB1 if
    ! not specified.
    IF (Factor(1) > 0.0) THEN
      TropSTCode = STCode(1)
    ELSE
      ! Old Limits
      Factor(1) = 10100.0
      Factor(2) = 199.0
      Factor(3) = 16180.0
      TropSTCode = ST_TropP
    END IF

    NumLevs = TP_NumLevs_Gl
    ZeroLev = TP_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                                &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )

    ! Do calculations
! DEPENDS ON: tropheight
    CALL TropHeight( NumLevs,                                          &
                     TempFields(PSource(1):PSource(NumLevs)),          &
                     TempFields(TSource(1):TSource(NumLevs)),          &
                     TempFields(ZSource(1):ZSource(NumLevs)),          &
                     TropSTCode, Factor(1:3),                          &
                     Fields(Store(1)), Fields(Store(2)),               &
                     Fields(Store(3)), ErrorStatus )

  ELSE IF ( Action == "ISOTHERM" )  THEN

   ! 2.5.2 Isotherm for temperature given in Factor
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run ISOTHERM on a 38, 50, or 70 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = IT_NumLevs_Gl
    ZeroLev = IT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                                &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )
    ! Read pstar
    STCode  (1) = ST_Pstar
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,                &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: isotherm
    CALL IsoTherm( NumLevs, Factor(1), OrogField, TempFields(StdPosn+1), &
                   TempFields(PSource(1):PSource(NumLevs)),              &
                   TempFields(TSource(1):TSource(NumLevs)),              &
                   TempFields(ZSource(1):ZSource(NumLevs)),              &
                   Fields(Store(1)), Fields(Store(2)),                   &
                   ErrorStatus )

  ELSE IF ( Action == "CONTRAIL" )  THEN

   ! 2.5.3 Upper and Lower Contrail Limit
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run CONTRAIL on a 38, 50, or 70 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = CT_NumLevs_Gl
    ZeroLev = CT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                                &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )
     ! Read pstar
    STCode  (1) = ST_Pstar
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,                &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: contrail
    CALL Contrail( NumLevs, TempFields(StdPosn + 1),                   &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   TempFields(TSource(1):TSource(NumLevs)),            &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   Fields(Store(1)), Fields(Store(2)),                 &
                   ErrorStatus )

  ELSE IF ( Action == "MAXWIND" )   THEN

   ! 2.5.4 Maximum Wind U&V Values and Pressure level on B grid
    IF ( P_Levels == UM_NumLevs ) THEN
      NumLevs = MX_NumLevs_Gl
    ELSE
      ! Find lowest level with Zsea >= MX_Zsea_Upr
      NumLevs = 1
      n = ModLevHdr_in % Len1LevDepC * 6
      DO i = 1, ModLevHdr_in % Len1LevDepC
        Zsea = ModLevHdr_in % LevDepC(i + n)
        IF ( (Zsea < MX_Zsea_Upr) .AND. (Zsea > 0) ) THEN
          NumLevs = NumLevs + 1
        END IF
      END DO
      WRITE(6,'(A)') "MaxWind levels reset to ", NumLevs
    END IF

    ZeroLev = MX_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read U @ 300hPa to use as a template for the B-grid
    STCode   = ST_Ustd
    MO8Level = 300
    Source   = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,                &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in, TempFields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read rho level u-component of wind & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, USource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, VSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, PSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Do calculations

    ! Store1/2/3 : Max Wind U/V/P

! DEPENDS ON: maxwind
    CALL MaxWind( NumLevs,                                             &
                  TempFields(USource(1):USource(NumLevs)),             &
                  TempFields(VSource(1):VSource(NumLevs)),             &
                  TempFields(PSource(1):PSource(NumLevs)),             &
                  Fields(Store(1)),                                    &
                  Fields(Store(2)),                                    &
                  Fields(Store(3)),                                    &
                  ErrorStatus )

  ELSE IF ( Action == "TOPBASE" )   THEN

    ! 2.5.5 Max wind Top and Base
    ! Action MAXWIND must be called first
    IF (PrevAction /= "MAXWIND") THEN
      CALL EReport( ProgName, ErrorStatus, &
                    "TOPBASE requires MAXWIND action immediately before." )
      ErrorStatus = StatusWarning
      GOTO 9990
    END IF


    ! Source1/2/3 : Max Wind U/V/P
    ! Store1/2    : Max Wind Base/Top

! DEPENDS ON: topbase
    CALL TopBase( NumLevs,                                      &
                  TempFields(USource(1):USource(NumLevs)),      &
                  TempFields(VSource(1):VSource(NumLevs)),      &
                  TempFields(PSource(1):PSource(NumLevs)),      &
                  Fields(Source(1)),                            &
                  Fields(Source(2)),                            &
                  Fields(Source(3)),                            &
                  Fields(Store(1)),                             &
                  Fields(Store(2)),                             &
                  ErrorStatus )

  ELSE IF ( Action == "CATURB" )    THEN

   ! 2.5.6 Clear Air Turbulence Predictor
    IF ( (P_Levels /= UM_NumLevs ) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,'(A)') "Can only run CATURB for 38, 50, or 70 level Global field"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = CA_NumLevs_Gl
    ZeroLev = CA_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read standard level fields
    STcode  (1:CA_NumStd) = (/ ST_Ustd, ST_Ustd, ST_Ustd,              &
                               ST_Vstd, ST_Vstd, ST_Vstd /)
    MO8Level(1:CA_NumStd) = (/     300,     250,     200,              &
                                   300,     250,     200 /)
    Source  (1:CA_NumStd) = StdPosn + (/ (i, i=1,CA_NumStd) /)
! DEPENDS ON: getflds
    CALL GetFlds( CA_NumStd, MaxFlds, STCode, MO8Level, FCTime,        &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in, TempFields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)
    ! Read rho level u-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, USource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, VSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, PSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Calculate rho level heights using B(UV)-grid orography
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
    ! Interpolate Orography to B-grid
! DEPENDS ON: ctobgrid
    CALL CtoBgrid( OrogField, TempFields(StdPosn+1) % Hdr,             &
                   TempFields(StdPosn+10), ErrorStatus )

! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, TempFields(StdPosn+10),                   &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )

    ! Do calculations
    DO i = 1,3
! DEPENDS ON: caturb
      CALL CATurb( NumLevs,                                         &
                   TempFields(StdPosn+i), TempFields(StdPosn+3+i),  &
                   TempFields(USource(1):USource(NumLevs)),         &
                   TempFields(VSource(1):VSource(NumLevs)),         &
                   TempFields(PSource(1):PSource(NumLevs)),         &
                   TempFields(ZSource(1):ZSource(NumLevs)),         &
                   Fields(Store(i)), ErrorStatus )
    END DO
! DEPENDS ON: maxcaturb
    CALL MaxCATurb( Fields(Store(1)), Fields(Store(2)),             &
                    Fields(Store(3)),                               &
                    Fields(Store(4)), Fields(Store(5)),             &
                    ErrorStatus )

  ELSE IF ( Action == "MWTURB" )    THEN

   ! 2.5.7 Mountain Wave Turbulence Predictor
    IF ( (P_Levels /= UM_NumLevs ) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,'(A)') "Can only run MWTURB for 38, 50, or 70 level Global field"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = MW_NumLevs_Gl
    ZeroLev = MW_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read standard level fields
    Source  (1:MW_NumStd) = StdPosn + (/3,6/)
    STCode  (1:MW_NumStd) = (/ ST_Ustd, ST_Vstd /)
    MO8Level(1:MW_NumStd) = (/     200,     200 /)
! DEPENDS ON: getflds
    CALL GetFlds( MW_NumStd, MaxFlds, STCode, MO8Level, FCTime,        &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in,TempFields,ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read u-component of gravity wave stress & interpolate to B-grid
    STCode  (1:NumLevs) = ST_GWSU
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, USource, PPHdrMod,            &
                    TempFields(StdPosn+3) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read v-component of gravity wave stress & interpolate to B-grid
    STCode  (1:NumLevs) = ST_GWSV
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, VSource, PPHdrMod,            &
                    TempFields(StdPosn+3) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: mtnstress
    CALL MtnStress( NumLevs,                                           &
                    TempFields(StdPosn+3), TempFields(StdPosn+6),      &
                    TempFields(USource(1):USource(NumLevs)),           &
                    TempFields(Vsource(1):VSource(NumLevs)),           &
                    Fields(Store(1)), ErrorStatus )

  ELSE IF ( Action == "WAFC_TURB" )    THEN

   ! 2.5.8 WAFC CAT turb Predictor (includes shear & mountain wave turb)
    IF ( (P_Levels /= UM_NumLevs) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,'(A)') "Can only run WAFC_TURB on a 38/50/70 level Global field"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = WT_NumLevs_Gl
    ZeroLev = WT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Parameters used in mountain wave calculation
    IF (ALL(Factor(1:3) == 0.0)) THEN
      ! Original values if nothing set in Factor.
      Factor(1) = 0.0625
      Factor(2) = (8.0/0.1875)
      Factor(3) = 1.833333
    END IF

    !---Calculate shear turbulence predictor (TI1)--------------
    ! This will be calculated for 400, 350, 300, 250, 200, 150 hPa & put into
    ! positions Store(i) to Store(i+6).

    ! Read standard level fields.  Specified now in the namelist
    ! Get information from numflds and if not specified set defaults.
    IF (NumFlds == 0) THEN
      NumFlds = WT_NumStd / 2
      MO8Level(1:NumFlds) = (/ 400, 300, 250, 200, 150 /)
    END IF

    STcode  (1:NumFlds)           = (/ ( ST_Ustd, i=1,NumFlds ) /)
    STcode  (NumFlds+1:2*NumFlds) = (/ ( ST_Vstd, i=1,NumFlds ) /)
    ! Copy MO8Level from numflds to end
    MO8Level(NumFlds+1:2*NumFlds) = (/ ( MO8Level(i), i=1,NumFlds ) /)

    ! To make the next bit easier (i.e. to mirror previous code) multiply
    ! NumFlds by 2.
    NumFlds = 2*NumFlds

    Source  (1:NumFlds) = StdPosn + (/ (i, i=1,NumFlds) /)
! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in, TempFields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)


    ! Read rho level u-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, USource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, VSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,        &
                    MinsPastHr, MinsDif, PSource, PPHdrMod,            &
                    TempFields(StdPosn+1) % Hdr,ModLevHdr_in,          &
                    TempFields,ErrorStatus )

    ! Calculate rho level heights using B(UV)-grid orography
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
    ! Interpolate Orography to B-grid
! DEPENDS ON: ctobgrid
    CALL CtoBgrid( OrogField, TempFields(StdPosn+1) % Hdr,             &
                   TempFields(StdPosn+NumFlds+1), ErrorStatus )

! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, TempFields(StdPosn+NumFlds+1),            &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )

    ! Step through each pressure level & calculate TI1 Index
    DO i = 1, NumFlds / 2
! DEPENDS ON:wafc_caturb
      CALL WAFC_CATurb( NumLevs,                                       &
                   TempFields(StdPosn+i),                              &
                   TempFields(StdPosn+(NumFlds/2)+i),                  &
                   TempFields(USource(1):USource(NumLevs)),            &
                   TempFields(VSource(1):VSource(NumLevs)),            &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   Fields(Store(i)),ErrorStatus )
    END DO


    !-------------------------------------------------------------
    ! Fields only available to calculate mountain wave part up to T+24
    ! After that just duplicate the shear turbulence fields

    IF (FCTime(1) <= 24) THEN

    !---Calculate mountain wave turbulence predictor --------------
    ! This will be calculated for 300, 250, 200, 150 hPa & put into
    ! positions Store(6) to Store(9).
    ! The standard level fields u and v are re-used from the shear CAT part 
    ! rather than re-reading them in. The positions for the model level fields
    ! in the shear CAT part are re-used to read in the required model fields
    ! for the MW predictor calculation.
      NumLevs = MW_NumLevs_Gl
      ZeroLev = MW_ZeroLev_Gl
      FCTime  (1:NumLevs) = FCTime(1)

      !Read model level fields (u and v component of gravity wave stress etc)
      MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

      ! Read u-component of gravity wave stress & interpolate to B-grid
      STCode  (1:NumLevs) = ST_GWSU
      USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,      &
                      MinsPastHr, MinsDif, USource, PPHdrMod,          &
                      TempFields(StdPosn+3) % Hdr,ModLevHdr_in,        &
                      TempFields,ErrorStatus )

      ! Read v-component of gravity wave stress & interpolate to B-grid
      STCode  (1:NumLevs) = ST_GWSV
      VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,      &
                      MinsPastHr, MinsDif, VSource, PPHdrMod,          &
                      TempFields(StdPosn+3) % Hdr,ModLevHdr_in,        &       
                      TempFields,ErrorStatus )

      ! Read pressure & interpolate to B-grid
      STCode  (1:NumLevs) = ST_Ptheta
      PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,      &
                      MinsPastHr, MinsDif, PSource, PPHdrMod,          &
                      TempFields(StdPosn+3) % Hdr,ModLevHdr_in,        &
                      TempFields,ErrorStatus)

      ! Loop through and do calculations for each reference pressure
      ! level required defined in array MW_PRef(i)
      ! Store results in the  Fields elements after the CAT fields
      ! (i.e. in the locations after Store((WT_NumStd/2)*2)   )

      ! Loop through and do calculations for each reference pressure
      ! level required defined in array MW_PRef(i)
      ! Store results in the  Fields elements after the CAT fields
      ! (i.e. in the locations after Store((WT_NumStd/2)*2)   )

      ! Need to loop through fields to see which pressure levels we have.
      j = 0
      DO i = 1, NumFlds/2
        SELECT CASE (TempFields(StdPosn+i) % Hdr % MO8Level)
        CASE(300, 250, 200, 150)
          j = j + 1
        END SELECT
      END DO

      ! Set number of pressure levels needed.
      MW_n_Pref = j

      j = 1
      DO i = 1, NumFlds/2
        SELECT CASE (TempFields(StdPosn+i) % Hdr % MO8Level)
        CASE(300, 250, 200, 150)
! DEPENDS ON: mtnstress_pref
          CALL MtnStress_PRef( NumLevs,                                &
                               Factor(1:4),                            &
                               MW_PRef(j),                             &
                               TempFields(StdPosn+i),                  &
                               TempFields(StdPosn+(NumFlds/2)+i),      &
                               TempFields(USource(1):USource(NumLevs)),&
                               TempFields(Vsource(1):VSource(NumLevs)),&
                               TempFields(Psource(1):PSource(NumLevs)),&
                               Fields( Store(j+(NumFlds/2)) ),         &
                               ErrorStatus)
          j = j + 1
        END SELECT
      END DO


      ! Combine CAT and MW fields by taking the maximum of the fields at each gridpoint
      ! for each level where both CAT and MW turbulence are produced (i.e.
      ! 300, 250, 200, 150 hPa (400 hPa no produced for MW turbulence)).

      j = 1
      DO i = 1, NumFlds/2
        SELECT CASE (TempFields(StdPosn+i) % Hdr % MO8Level)
! The number in this case statement should equal the MW_n_Pref.
        CASE(300, 250, 200, 150)
!CAT turbulence field
!MW turbulence field
!combined field
! DEPENDS ON: maximum
          CALL Maximum( Fields( Store(i) ),                        &
                        Fields( Store( (NumFlds/2)+j ) ),          &
                        Fields( Store( (NumFlds/2)+MW_n_PRef+j ) ),&
                        ErrorStatus )
          j = j + 1
        END SELECT
      END DO


      !Duplicate the combined CAT/MW fields so we have fields for each level
      !labelled with mean and maximum CAT stash codes

      !for 400mb,350mb (contains CAT only) (generally above 300mb)

      DO i = 1, NumFlds/2
        SELECT CASE (TempFields(StdPosn+i) % Hdr % MO8Level)
        CASE(301:)
! DEPENDS ON: dup_cat
          CALL Dup_CAT( Fields(Store(i)), &             !input CAT/MW field
               Fields(Store((NumFlds/2)+(2*MW_n_PRef)+i)), & !output fld
               ErrorStatus)
        END SELECT
      END DO

      !Loop through 300mb, 250mb, 200mb, 150mb combined CAT/MW fields
      j = 1
      DO i = 1, NumFlds/2
        SELECT CASE (TempFields(StdPosn+i) % Hdr % MO8Level)
        CASE(300, 250, 200, 150)
!input CAT/MW field
!output duplicate field
! DEPENDS ON: dup_cat
          CALL Dup_CAT( Fields(Store((NumFlds/2)+MW_n_PRef+j)),       &
                        Fields(Store((NumFlds/2)+(2*MW_n_PRef)+i)),   &
                        ErrorStatus )
          j = j + 1
        END SELECT
      END DO

    ELSE

      !Mountain wave fields not available - just duplicate shear CAT fields
      DO i=1, NumFlds/2

! DEPENDS ON: dup_cat
        CALL Dup_CAT( Fields( Store(i) ),            &  !input CAT field
                      Fields( Store((NumFlds/2)+i) ),&  !output duplicate field
                      ErrorStatus )

      END DO

    ENDIF

  ELSE IF ( Action == "OLD_ICING" )  THEN

   ! 2.5.9 WAFC Icing Predictor (Large scale cloud layer Icing algorithm LIBLKCLD)
   !       BuLK CLouD fraction for cloud in a given temperature range

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run LIBLKCLD on a 38, 50, or 70 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = IC_NumLevs_Gl
    ZeroLev = IC_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, USource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read theta level pressure in Pa
    STCode  (1:NumLevs) = ST_Ptheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Find theta level gridboxes in the given temp range
! DEPENDS ON: temperature_mask
    CALL Temperature_Mask( Factor(1:2), NumLevs,                       &
                   TempFields(TSource(1):TSource(NumLevs)),            &
                   TempFields(VSource(1):VSource(NumLevs)),            &
                   ErrorStatus )

    ! Calculate theta level cloud fraction for the temp range
! DEPENDS ON: multiply
    CALL Multiply( NumLevs,                                            &
                   TempFields(USource(1):USource(NumLevs)),            &
                   TempFields(VSource(1):VSource(NumLevs)),            &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   ErrorStatus )

    ! Calculate pressure layer cloud fraction for the temp range
! DEPENDS ON: iconplyr
    CALL ICOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L(1:NumLyrs),                                 &
                   LyrLwrB(1:NumLyrs),                                 &
                   LyrUprB(1:NumLyrs),                                 &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   Fields(Store(1):Store(2*NumLyrs)),                  &
                   ErrorStatus )



  ELSE IF ( Action == "ICING" )  THEN

!-----------------------------------------------------------------------------!
!                 Relative Humidity based icing diagnostic                    !
!-----------------------------------------------------------------------------!
! Method:                                                                     !
! 1. Read in pressure, temperature and specific humidity on model levels.     !
!    Calculate saturated mixing ratio from P & T and then combine this with   !
!    specific humidity to get RH (detailed in Rel_humdy subroutine). Then     !
!    determine RH on pressure levels.                                         !
! 2. Read in cloud fraction on model levels. Create mask where cloud is       !
!    present. Determine cloud fraction on pressure levels.                    !
! 3. Create -20 to 0C temperature mask. Convert on to pressure levels.        !
! 4. Multiply the three pressure level arrays to give RH where cloud is       !
!    present and temperature is between 0 and -20C                            !
!-----------------------------------------------------------------------------!

   ! 2.5.9 WAFC Icing Predictor
   !       Relative humidity for cloud in a given temperature range

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run Icing algorithm on a 38/50/70 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = LI_NumLevs_Gl
    ZeroLev = LI_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)


    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

!----------------------------------------------------------------
!---- Step 1: Read in pressure,temperature, specific humidity.---
!---- Calculate RH and interpolate on to pressure levels      ---
!----------------------------------------------------------------

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )


    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )


    ! Read theta level specific humidity
    STCode  (1:NumLevs) = ST_Htheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, USource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Calculate saturated mixing ratio on model levels
! DEPENDS ON: sat_mratio
    CALL Sat_MRatio( NumLevs,                                          &
                     TempFields(PSource(1):PSource(NumLevs)),          &
                     TempFields(TSource(1):TSource(NumLevs)),          &
                     TempFields(VSource(1):VSource(NumLevs)),          &
                     ErrorStatus )

    ! Calculate relative humidity on model levels
! DEPENDS ON: rel_humdy
    CALL Rel_humdy( NumLevs,                                           &
                    TempFields(USource(1):USource(NumLevs)),           &
                    TempFields(VSource(1):VSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    IF (ErrorStatus == StatusOk) THEN
      DO i=1,NumLevs
        DEALLOCATE ( TempFields(USource(i)) % RData )
        DEALLOCATE ( TempFields(VSource(i)) % RData )
        NULLIFY( TempFields(USource(i)) % RData )
        NULLIFY( TempFields(VSource(i)) % RData )
      END DO
    END IF

    ! Calculate relative humidity on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L, LyrLwrB, LyrUprB,                          &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   TempFields(USource(1):USource(NumLyrs)),            &
                   ErrorStatus )
    ! Use MIN to make sure relative humidity is in range 0->1 otherwise we get
    ! unrealistic ICING potential.
    DO k=1, NumLyrs
      DO j=1, Fields(USource(k)) % Hdr % NumRows
        DO i=1, Fields(USource(k)) % Hdr % NumCols
          TempFields(USource(k)) % Rdata(i,j) =                        &
            MIN(TempFields(USource(k)) % Rdata(i,j),1.0)
        END DO
      END DO
    END DO

! Nullify bit of the array that's about to be reused
    IF (ErrorStatus == StatusOk) THEN
      DO i=1, NumLevs
        DEALLOCATE ( TempFields(ZSource(i)) % RData )
        NULLIFY( TempFields(ZSource(i)) % RData )
      END DO
    END IF


!----------------------------------------------------------------
!---- Step 2: Determine temperature on pressure levels        ---
!---- by calulating mean temperature across atmospheric layer ---
!----------------------------------------------------------------

    ! Determine temperature on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L, LyrLwrB, LyrUprB,                          &
                   TempFields(TSource(1):TSource(NumLevs)),            &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   TempFields(ZSource(1):ZSource(NumLyrs)),            &
                   ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    IF (ErrorStatus == StatusOk) THEN
      DO i=1, NumLevs
        DEALLOCATE( TempFields(TSource(i)) % RData )
        NULLIFY( TempFields(TSource(i)) % RData )
      END DO
    END IF

    ! Calculate which gridboxes have temperatures between -20 and 0C
! DEPENDS ON: temperature_mask
    CALL Temperature_Mask( Factor(1:2), NumLyrs,                       &
                           TempFields(ZSource(1):ZSource(NumLyrs)),    &
                           TempFields(TSource(1):TSource(NumLyrs)),    &
                           ErrorStatus )

!----------------------------------------------------------------
!---- Step 3: Read in cloud fraction on model levels         ----
!---- Determine cloud fraction on pressure levels            ----
!----------------------------------------------------------------

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, ZSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )


    ! Determine cloud fraction on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L, LyrLwrB, LyrUprB,                          &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   TempFields(PSource(1):PSource(NumLevs)),            &
                   TempFields(VSource(1):VSource(NumLyrs)),            &
                   ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    IF (ErrorStatus == StatusOk) THEN
      DO i=1, NumLevs
        DEALLOCATE( TempFields(ZSource(i)) % RData )
        NULLIFY( TempFields(ZSource(i)) % RData )
      END DO
    END IF


    ! Calculate which gridboxes have cloud present
! DEPENDS ON: cloud_mask
    CALL Cloud_Mask( (/ 0.0, 1.1 /), NumLyrs,                          &
                     TempFields(VSource(1):VSource(NumLyrs)),          &
                     TempFields(ZSource(1):ZSource(NumLyrs)),          &
                     ErrorStatus )

! Transfer masked cloud mask back into VSource section, has to be done this way
! due to rdata being a pointer.  If we perform field1 = field2 then the field1
! rdata will point to field2 and will not copy the data.
    IF (ErrorStatus == StatusOk) THEN
      DO i=1, NumLyrs
        TempFields(VSource(i)) % Hdr       = TempFields(ZSource(i)) % Hdr
        TempFields(VSource(i)) % Rdata     = TempFields(ZSource(i)) % Rdata
        TempFields(VSource(i)) % LookupPos = TempFields(ZSource(i)) % LookupPos

        DEALLOCATE ( TempFields(ZSource(i)) % RData )
        NULLIFY( TempFields(ZSource(i)) % RData )
      END DO
    END IF

!-----------------------------------------------------------------
!---- Step 4: Apply temperature and cloud masks to             ---
!----         relative humidity                                ---
!-----------------------------------------------------------------

    ! Apply temperature mask to relative humidity fields
! DEPENDS ON: multiply
    CALL Multiply( NumLyrs,                                            &
                   TempFields(USource(1):USource(NumLyrs)),            &
                   TempFields(TSource(1):TSource(NumLyrs)),            &
                   TempFields(ZSource(1):ZSource(NumLyrs)),            &
                   ErrorStatus )


    ! Apply cloud mask to relative humidity fields masked by temperature
! DEPENDS ON: multiply
    CALL Multiply( NumLyrs,                                            &
                   TempFields(ZSource(1):ZSource(NumLyrs)),            &
                   TempFields(VSource(1):VSource(NumLyrs)),            &
                   Fields(Store(1):Store(NumLyrs)),                    &
                   ErrorStatus )


!---------------------------------------------------------------
!---- Step 5: Duplicate icing predictor field and set       ----
!---- fieldscodes in the field headers.                     ----
!---------------------------------------------------------------

    DO i=1, NumLyrs
! DEPENDS ON: dup_ice
      CALL Dup_ice( Fields(Store(i)),                                  &
                    Fields(Store(NumLyrs+i)),                          &
                    ErrorStatus )
    END DO


!-----------------------------------------------------------------------------
!--- Output:                                                                 !
!---   Fields(Store(1):Store(NumLyrs) = fields marked as mean icing          !
!---   Fields(Store(NumLyrs+1):Store(2*NumLyrs) = fields marked as max icing !
!-----------------------------------------------------------------------------


  ELSE IF ( Action == "CLD_TURB" )  THEN

   ! 2.5.10 In-cloud turbulence algorithm

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run LIBLKCLD on a 38, 50, or 70 level model"
      ErrorStatus = StatusWarning
    END IF

    !Read in model fields
    NumLevs = ICT_NumLevs_Gl
    ZeroLev = ICT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read in temperature fields and calculate
    ! equivalent potential temperature fields

    ! Read theta level pressure in Pa
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )


    !Call subroutine to calculate equivalent potential temperature
    !on model levels (put into VSource space in Fields)

    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: eq_pot_temp
    CALL Eq_pot_temp( NumLevs,                                         &
                      TempFields(TSource(1):TSource(NumLevs)),         &
                      TempFields(PSource(1):PSource(NumLevs)),         &
                      TempFields(VSource(1):VSource(NumLevs)),         &
                      ErrorStatus )

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, USource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )


    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                                &
                    TempFields(PSource(1):PSource(NumLevs)),           &
                    TempFields(ZSource(1):ZSource(NumLevs)),           &
                    ErrorStatus )


   ! Do calculations

    MO8Level(1:ICT_NumStd) = (/ 700, 600, 500, 400, 300/)

    DO i = 1, ICT_NumStd

! DEPENDS ON: cld_turb
      CALL cld_turb(                                                   &
                     ! pressure in hPa (in)
                     MO8Level(i),                                      &
                     ! no. model levels (in)
                     NumLevs,                                          &
                     ! cloud on model level fields (in)
                     TempFields(USource(1):USource(NumLevs)),          &
                     ! equivalnt potential temp (in)
                     TempFields(VSource(1):VSource(NumLevs)),          &
                     ! corresponding pressure (in)
                     TempFields(PSource(1):PSource(NumLevs)),          &
                     ! corresponding height (in)
                     TempFields(ZSource(1):ZSource(NumLevs)),          &
                     ! output incloud turb predictor
                     Fields(Store(i)),                                 &
                     ErrorStatus )

      !Duplicate the field and label as the maximum
! DEPENDS ON: dup_ict
      CALL Dup_ICT( Fields(Store(i)),                       &  !input in-cloud turb field
                    Fields(Store(ICT_NumStd+i)),            &  !output duplicate field
                    ErrorStatus )
    END DO

  ELSE IF ( Action == "CB_ACT" )  THEN

   ! 2.5.11 Cumulonimbus (CB) fields

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,'(A)') "Can only run LIBLKCLD on a 38, 50, or 70 level model"
      ErrorStatus = StatusWarning
    END IF
    IF ( Factor(1) == 0.0 ) THEN
      ! Not set so use default (i.e. use fill_n_dspec)
      ! Following signifies whether to use fill_n_dspec
      ! (-1.0 = false, everything else = true)
      Factor(1) = 1.0
    END IF

! Read in standard level fields

! Fields(10) : Convective Cloud Base Pressure : 5/207
! Fields(11) : Convective Cloud Top  Pressure : 5/208
! Fields(12) : Convective Precipitation Rate  : 5/205
! Fields(13) : Lowest Conv Cloud Top Pressure : 5/223

! Each field is only on one level: level 8888

    STcode  (1:CB_NumStd) = (/ ST_CPNRT /)
    MO8Level(1:CB_NumStd) = (/ 8888 /)
    Source  (1:CB_NumStd) = 9 + (/ (i, i=1,CB_NumStd) /)
    FCTime  (1:CB_NumStd) = FCTime(1)
    ! Set to use mean value.
    LBProc  (1:CB_NumStd) = (/ 128 /)

! DEPENDS ON: getflds
    CALL GetFlds( CB_NumStd, MaxFlds, STCode, MO8Level, FCTime,       &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,      &
                  StdLevHdr_in, TempFields, ErrorStatus )


! Read in model level fields

    NumLevs = CB_NumLevs_Gl
    ZeroLev = CB_ZeroLev_Gl
! Reset LBProc to 0 for all subsequent fields.
    LBProc(1:NumLevs) = 0
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

!   TSource set later, ZSource not used here
!   TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
!   ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

    ! Read theta level bulk cloud fraction (0/266)

    STCode  (1:NumLevs) = ST_BlkCld
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs,MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, MinsDif, USource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read Convective Cloud Amount (5/212)

    STCode  (1:NumLevs) = ST_ConCld
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs,MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, MinsDif, VSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read theta level pressure in Pa (0/408)

    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Read theta level Temperature in K (16/004)

    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, TSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

! Do calculations

! Setup Store array - could be read from Namelist if required.
    IF ( ALL(Store == 0) ) THEN
      Store(1:9) = (/ (i, i=1,9) /)
    END IF

! DEPENDS ON: convact
    CALL ConvAct( NumLevs,                                             &
                  Factor(1),                                           &
                  TempFields(Source(1)),                               &
                  TempFields(USource(1):USource(NumLevs)),             &
                  TempFields(VSource(1):VSource(NumLevs)),             &
                  TempFields(PSource(1):PSource(NumLevs)),             &
                  TempFields(TSource(1):TSource(NumLevs)),             &
                  Fields(Store(1)),                                    &
                  Fields(Store(2)),                                    &
                  Fields(Store(3)),                                    &
                  Fields(Store(4)),                                    &
                  Fields(Store(5)),                                    &
                  Fields(Store(6)),                                    &
                  Fields(Store(7)),                                    &
                  Fields(Store(8)),                                    &
                  Fields(Store(9)),                                    &
                  Errorstatus)

  ELSE IF ( Action == "DUST_CONC" )  THEN

    ! 3.0.1 Convert the dust concentration into g/m3

    IF ( P_Levels /= UM_NumLevs ) THEN
       WRITE(6,'(A)') "Can only run DUST_CONC on a 38/50/70 level model"
       ErrorStatus = StatusWarning
    END IF

    ! for the dust calc have chosen to read in bottom theta level

    NumLevs = DT_NumLevs_Gl
    ZeroLev = DT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read rho level pressure  (0/407)

    STCode  (1:NumLevs) = ST_Prho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, USource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Calculate rho level heights

    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                                &
              TempFields(USource(1):USource(NumLevs)),                 &
              TempFields(ZSource(1):ZSource(NumLevs)),                 &
              ErrorStatus )

    ! Read total dust concentration (17/257)

    STCode  (1:NumLevs) = ST_Dust
    DSource (1:NumLevs) = DSTPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, DSource, PPHdrMod,      &
                  ModLevHdr_in, TempFields, ErrorStatus )

    ! Calculate surface dust concentration

! DEPENDS ON: dust_con
    CALL DUST_CON (Numlevs, OrogField,                                 &
                   TempFields(ZSource(1):ZSource(NumLevs)),            &
                   TempFields(DSource(1):DSource(NumLevs)),            &
                   Fields(Store(1)),Fields(Store(2)),                  &
                   ErrorStatus )

  ELSE IF ( Action == "DUST_VIS" )  THEN
 
    ! Read in StdLev fields
    NumFlds = 2
    FCTime  (1:NumFlds) = FCTime(1)
    ! model visibility including precipitation
    STCode  (1) = ST_VisIncPpn
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
    ! 1.5m temperature
    STCode  (2) = ST_TScreen
    MO8Level(2) = LV_Surface
    Source  (2) = StdPosn + 2

! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in, TempFields, ErrorStatus )

    ! Read in MdlLev fields
    NumFlds = 7
    FCTime  (3:9) = FCTime(1)
    ! surface pressure
    STCode  (3) = ST_Pstar
    MO8Level(3) = LV_Surface
    Source  (3) = StdPosn + 3 
    ! Level 1 dust fields
    STCode  (4) = ST_Dust1 
    STCode  (5) = ST_Dust2 
    STCode  (6) = ST_Dust3 
    STCode  (7) = ST_Dust4 
    STCode  (8) = ST_Dust5 
    STCode  (9) = ST_Dust6 
    MO8Level(4) = 1
    MO8Level(5) = 1
    MO8Level(6) = 1
    MO8Level(7) = 1
    MO8Level(8) = 1
    MO8Level(9) = 1
    Source  (4) = StdPosn + 4
    Source  (5) = StdPosn + 5
    Source  (6) = StdPosn + 6
    Source  (7) = StdPosn + 7
    Source  (8) = StdPosn + 8
    Source  (9) = StdPosn + 9

! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode(3:9), MO8Level(3:9),        &
                  FCTime(3:9),                                         &
                  LBProc, MinsPastHr, MinsDif, Source(3:9), PPHdrMod,  &
                  ModLevHdr_in, TempFields, ErrorStatus )

! DEPENDS ON: dust_vis
    CALL DUST_VIS (TempFields(Source(1)),TempFields(Source(2)),        &
                   TempFields(Source(3)),TempFields(Source(4)),        &
                   TempFields(Source(5)),TempFields(Source(6)),        &
                   TempFields(Source(7)),TempFields(Source(8)),        &
                   TempFields(Source(9)),                              &
                   Fields(Store(1)),Fields(Store(2)),                  &
                   ErrorStatus )
!
! Temporary duplicate for dust_vis for 2 bin scheme.
! This can be removed and the above action simplified when the model 
! is changed to output the extinction directly from STASH.
!
  ELSE IF ( Action == "DUST_VIS2" )  THEN
 
    ! Read in StdLev fields
    NumFlds = 2
    FCTime  (1:NumFlds) = FCTime(1)
    ! model visibility including precipitation
    STCode  (1) = ST_VisIncPpn
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
    ! 1.5m temperature
    STCode  (2) = ST_TScreen
    MO8Level(2) = LV_Surface
    Source  (2) = StdPosn + 2

! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, MinsDif, Source, PPHdrMod,       &
                  StdLevHdr_in, TempFields, ErrorStatus )

    ! Read in MdlLev fields
    NumFlds = 3
    FCTime  (3:5) = FCTime(1)
    ! surface pressure
    STCode  (3) = ST_Pstar
    MO8Level(3) = LV_Surface
    Source  (3) = StdPosn + 3 
    ! Level 1 dust fields
    STCode  (4) = ST_Dust1 
    STCode  (5) = ST_Dust2 
    MO8Level(4) = 1
    MO8Level(5) = 1
    Source  (4) = StdPosn + 4
    Source  (5) = StdPosn + 5

! DEPENDS ON: getflds
    CALL GetFlds( NumFlds, MaxFlds, STCode(3:5), MO8Level(3:5),        &
                  FCTime(3:5),                                         &
                  LBProc, MinsPastHr, MinsDif, Source(3:5), PPHdrMod,  &
                  ModLevHdr_in, TempFields, ErrorStatus )

! DEPENDS ON: dust_vis2
    CALL DUST_VIS2(TempFields(Source(1)),TempFields(Source(2)),        &
                   TempFields(Source(3)),TempFields(Source(4)),        &
                   TempFields(Source(5)),                              &
                   Fields(Store(1)),Fields(Store(2)),                  &
                   ErrorStatus )

! ZTD diagnostic
  ELSE IF ( Action == "ZTD" ) THEN

! algorithm uses all levels in the model
! gives model level 1

    NumLevs = UM_numlevs
    ZeroLev = 0

    DO i = 1, NumFlds
      ! Set Met08 levels - loop increases them from 1 to
      ! Numlevs (=um_numlevs)
      MO8Level(1:NumLevs) = (/ (j+ZeroLev, j=1,NumLevs) /)

      ! Read theta level pressure (0/408)
      STCode  (1:NumLevs) = ST_Ptheta
      PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)

!DEPENDS ON: getflds
      CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level,                &
                    (/ ( FCTime(i),j=1,NumLevs ) /),                   &
                    LBProc, MinsPastHr, MinsDif, PSource, PPHdrMod,    &
                    ModLevHdr_in, TempFields, ErrorStatus )


      ! Read specific humidity on theta (0/10)
      STCode  (1:NumLevs) = ST_Htheta
      QSourceTheta(1:NumLevs) = QPosnTheta + MO8Level(1:NumLevs)

!DEPENDS ON: getflds
      CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level,                &
                    (/ ( FCTime(i),j=1,NumLevs ) /),                   &
                    LBProc, MinsPastHr, MinsDif, QSourceTheta,         &
                    PPHdrMod, ModLevHdr_in, TempFields, ErrorStatus )

      ! Read Rho level Exner (0/255)
      ! Exner goes from rho level 1 to TOP+1
      MO8Level(1:NumLevs+1) = (/ (i+ZeroLev, i=1,NumLevs+1) /)
      STCode  (1:NumLevs+1) = ST_ExnerRho

      ExnerSource (1:NumLevs+1) = ExnerPosn + MO8Level(1:NumLevs+1)

!DEPENDS ON: getflds
      CALL GetFlds( NumLevs+1, MaxFlds, STCode, MO8Level,              &
                    (/ ( FCTime(i),j=1,NumLevs+1 ) /),                 &
                    LBProc, MinsPastHr, MinsDif, ExnerSource,          &
                    PPHdrMod, ModLevHdr_in, TempFields, ErrorStatus )

      ZSourceRho (1:NumLevs+1) = ZPosnRho + MO8Level(1:NumLevs+1)

      ! Call calczflds to calculate rho level heights (add 1 to NumLevs)
!DEPENDS ON: calczflds
      CALL CalcZFlds( NumLevs+1, OrogField,                                  &
                      TempFields(ExnerSource(1):ExnerSource(NumLevs+1)),     &
                      TempFields(ZSourceRho(1):ZSourceRho(NumLevs+1)),       &
                      ErrorStatus )

      ! Call zenithdelay
!DEPENDS ON: zenithdelay
      CALL zenithdelay( NumLevs, StdLevHdr_in % IntC(9),                     &
                        TempFields(PSource(1):PSource(NumLevs)),             &
                        TempFields(ExnerSource(1):ExnerSource(NumLevs+1)),   &
                        TempFields(ZSourceRho(1):ZSourceRho(NumLevs+1)),     &
                        TempFields(QSourceTheta(1):QSourceTheta(NumLevs)),   &
                        Fields(Store(i)), ErrorStatus )
    END DO

  ELSE IF ( Action == "THINFLD" ) THEN

!DEPENDS ON: thinfield
      CALL thinfield(Fields(Source(1)),Fields(Store(1)),Factor(1:3),ErrorStatus)

  ELSE IF ( Action == "CUTOUT" ) THEN

    IF (ALL(Factor(1:2) == 0.0)) THEN
      ! Original values if nothing set in Factor - fixed grid (2,2)
      Factor(1) = 2.0
      Factor(2) = 2.0
    END IF

!DEPENDS ON: cutout  
    CALL cutout(StdLevHdr_in, StdLevHdr_out, Factor(1:6), pphdrmod)
  ELSE IF ( Action == "GRIBIFY" ) THEN

    IF (gribfile == "unset") THEN
      WRITE(6,'(A)') "GRIBFILE env. var. not set.  Cannot continue."
      Errorstatus = StatusWarning
      GOTO 9990
    END IF

    IF ( LevType == "STD" ) THEN
      ! DEPENDS ON: gribify
      CALL gribify(StdLevHdr_in, gribfile, errorstatus)
    ELSE IF ( Levtype == "MOD" ) THEN
      ! DEPENDS ON: gribify
      CALL gribify(ModLevHdr_in, gribfile, errorstatus)
    ELSE
      WRITE(6,'(3A)') "Cannot process ", Levtype, " level type input."
      errorstatus=statuswarning
    END IF
 
  ELSE IF ( Action == "DEGRIBIFY" ) THEN

    IF (gribfile == "unset") THEN
      WRITE(6,'(A)') "GRIBFILE env. var. not set.  Cannot continue."
      Errorstatus = StatusWarning
      GOTO 9990
    END IF

    IF      ( LevType == "GRIB" ) THEN
! DEPENDS ON: degribify
      CALL degribify(StdLevHdr_out, MaxFldsOut, gribfile, packType, errorstatus)
    ELSE
      WRITE(6,'(3A)') "Cannot process ", Levtype, " level type input."
      errorstatus=statuswarning
    END IF

  ELSE IF ( Action == "CLEAR_TMP" ) THEN
    ! Deallocate any Temporary field space created within the actions.
    DO i = 1, MaxFlds
      IF ( ASSOCIATED( TempFields(i) % RData ) ) THEN
        DEALLOCATE( TempFields(i) % RData )
        NULLIFY( TempFields(i) % RData )
      END IF
    END DO

  ELSE

   ! Unrecognised Action call
    ErrorStatus = StatusWarning

    CALL EReport( ProgName, ErrorStatus, &
                  "Command not recognised.  Action = " // Action )
    ErrorStatus = StatusWarning

  END IF
9990 CONTINUE
  ! Check if current / previous error has caused Action to be aborted
  IF (ErrorStatus /= StatusOK) THEN
    ErrorStatus = StatusWarning

    CALL EReport( ProgName, ErrorStatus, &
                  "Error - skipping to next Action.")
    ErrorStatus = StatusWarning
  END IF

  PrevAction = Action
END DO

9999 CONTINUE

!-----------------------------------------------------------------------
! 3 CLOSE INPUT FILES
!-----------------------------------------------------------------------
IF ( ASSOCIATED( OrogField % RData ) ) THEN
  DEALLOCATE( OrogField % RData )
  NULLIFY( OrogField % RData )
END IF
IF ( ASSOCIATED( LSMField % RData ) ) THEN
  DEALLOCATE( LSMField % RData )
  NULLIFY( LSMField % RData )
END IF

DO i = 1, MaxFlds
  IF ( ASSOCIATED( Fields(i) % RData ) ) THEN
    DEALLOCATE( Fields(i) % RData )
    NULLIFY( Fields(i) % RData )
  END IF

  IF ( ASSOCIATED( TempFields(i) % RData ) ) THEN
    DEALLOCATE( TempFields(i) % RData )
    NULLIFY( TempFields(i) % RData )
  END IF
END DO

CLOSE( NLUnit )
IF ( ASSOCIATED( StdLevHdr_in % Lookup ) ) THEN
CALL File_Close ( StdLevHdr_in % UnitNum,                &
                  StdLevHdr_in % FileNameEnv,            &
                  LEN_TRIM(StdLevHdr_in % FileNameEnv),  &
                  0, 0, ErrorStatus )
END IF
IF ( ASSOCIATED( ModLevHdr_in % Lookup ) ) THEN
CALL File_Close ( ModLevHdr_in % UnitNum,                &
                  ModLevHdr_in % FileNameEnv,            &
                  LEN_TRIM(ModLevHdr_in % FileNameEnv),  &
                  0, 0, ErrorStatus )
END IF
IF ( ASSOCIATED( StdPrvHdr_in % Lookup ) ) THEN
  CALL File_Close ( StdPrvHdr_in % UnitNum,              &
                    StdPrvHdr_in % FileNameEnv,          &
                    LEN_TRIM(StdPrvHdr_in % FileNameEnv),&
                    0, 0, ErrorStatus )
END IF
!-----------------------------------------------------------------------
! 4 CLOSE OUTPUT FILE
!-----------------------------------------------------------------------
! Write out updated lookup table
! DEPENDS ON: writelookup
CALL WriteLookup( StdLevHdr_out, ErrorStatus )
WRITE(6,'(I7,A)') StdLevHdr_out % NumFlds, " fields written."

IF ( ASSOCIATED( StdLevHdr_out % Lookup ) ) THEN
  CALL File_Close ( StdLevHdr_out % UnitNum,               &
                    StdLevHdr_out % FileNameEnv,           &
                    LEN_TRIM(StdLevHdr_out % FileNameEnv), &
                    0, 0, ErrorStatus )
END IF
CALL ioShutdown()

CALL ereport_finalise( )

! DEPENDS ON: timer
CALL Timer( ProgName, 2 )

END PROGRAM FieldCalc
