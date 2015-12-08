! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to read requested fields from a UM fieldsfile.

SUBROUTINE ReadFld( UMHdr,           &  ! in
                    LDecode,         &  ! in
                    PPHdrMod,        &  ! in
                    Field,           &  ! inout
                    ErrorStatus )       ! inout

! Description:
!   This subroutine reads a given pp-field from a UM fieldsfile.
!   The location of the field in the lookup should already be set within
!   Field, and memory allocated to store the data.
!
! Method:
!   The type of file is ascertained from the pp-header, and the field
!   read in.  If the field is packed and LDecode is true, the field is
!   decompressed.  The data type is checked for compatibility with the
!   program.  If the input field is of a type requiring header
!   modification and PPHdrMod is true, the appropriate modification is
!   done.
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
  UM_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN) :: UMHdr    ! Fieldsfile header info
LOGICAL, INTENT(IN) :: LDecode               ! Unpack? indicator
LOGICAL, INTENT(IN) :: PPHdrMod              ! Mod of certain headers

TYPE(PP_Field_type), INTENT(INOUT) :: Field  ! Output field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ReadFld"
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
! Parameters to define values used as datatype in lookups.
INTEGER, PARAMETER :: Datatype_real = 1
INTEGER, PARAMETER :: Datatype_int  = 2
INTEGER, PARAMETER :: Datatype_log  = 3

! Local Variables:
INTEGER :: NumCrayWords      ! no of values in an input field
INTEGER :: i                 ! local counter
INTEGER :: DataLen           ! Length of a particular field
INTEGER :: InputAddr         ! Word Address in call SETPOS
INTEGER :: Addr              ! Address of a field in the data store
INTEGER :: PackType          ! Packing type N1 of LBPACK
INTEGER :: PackType_I        ! Packing type N1 of LBPACK in loop
INTEGER :: DataType          ! Input data type (integer, real etc.)
INTEGER :: Len_IO            ! Length of data written
INTEGER :: Year, Month, Date           ! Temp variables:
INTEGER :: Hour, Min, Sec, DayNo       !  used in date / time
INTEGER :: EndTimeDays, EndTimeSecs    !  conversion for
INTEGER :: DataTimeDays, DataTimeSecs  !  time-mean &
INTEGER :: FCRangeSecs                 !  accumulated fields
REAL :: Err_IO               ! Error Code
REAL :: AMDI                 ! Missing data indicator
LOGICAL :: Lcal360 = .FALSE. ! 30-day month indicator for time routines
CHARACTER(LEN=80) :: ErrMessage

INTEGER, ALLOCATABLE :: ITempData(:,:)

! End of header --------------------------------------------------------

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Decode LBPACK code
PackType  = MOD(Field % Hdr % LBPack,10)
! Length on disk is what we want to read in.
NumCrayWords = Field % Hdr % LBLRec
! Data position is where we need to read the file from
InputAddr    = Field % Hdr % DataPos

! For unknown/not well-formed records this might be zero, cannot guarantee it
! will work.
IF ( Field % Hdr % lbnrec == 0 ) THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, "LBNREC is 0 - will try to continue.")
! Reset ErrorStatus since we are continuing.
  ErrorStatus = StatusOK
END IF

! We are using lblrec as a indicator of how much space is used on disk.  Really
! should use lbnrec but will use lblrec for now.
IF ( MOD(Field % Hdr % lbpack,10) == 2) THEN
  ! We have 32 bit packing and this is probably a dump.  LBLREC will contain
  ! the number of points not the amount of memory required to store this data.
  ! Lets convert to number of 64 bit words.
  NumCrayWords = (NumCrayWords + 1)/2
END IF

CALL SETPOS( UMHdr % UnitNum, InputAddr, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( .NOT. LDecode ) THEN
  ! Read data without unpacking or number conversion
  CALL BUFFIN( UMHdr % UnitNum, Field % RData, NumCrayWords, &
               Len_IO, Err_IO )

ELSE
  ! Want real, unpacked data for calculations
  DataType = Field % Hdr % LBUser1
  IF ( DataType == Datatype_real ) THEN

    ! Real field
    CALL BUFFIN( UMHdr % UnitNum, Field % RData, NumCrayWords, &
                                Len_IO, Err_IO )
    ! We cannot use Packtype here since we might be land packed.  Packtype is
    ! more the type of numeric packing used - e.g. land packing is 120)
    IF (Field % Hdr % LBPack > 0 ) THEN        ! Is the field packed?

      AMDI = Field % Hdr % BMDI
      ! Compare with standard RMDI
      IF ( AMDI /= RMDI ) THEN
        WRITE(6,*)" WARNING non-standard MDI in use"
      END IF
! DEPENDS ON: unpackflds
      CALL UnPackFlds( PackType, Field % Hdr % NumCols, &
                       Field % Hdr % NumRows, Field,    &
                       ErrorStatus )
      IF ( ErrorStatus /= StatusOK ) THEN
        ErrorStatus = StatusWarning

        CALL EReport( RoutineName, ErrorStatus, &
                      "ReadFld - could not unpack input field" )
        ErrorStatus = StatusWarning
        GO TO 9999
      END IF
    END IF

  ELSE IF ( DataType == Datatype_int .OR. DataType == Datatype_log ) THEN

    ! Integer field
    ALLOCATE( ITempData(NumCrayWords,1) )
    IF ( PackType > 0 ) THEN
      ErrorStatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, &
                    "Cannot handle packed integer/logical data fields" )
      ErrorStatus = StatusWarning
      GO TO 9999
    END IF
    CALL BUFFIN( UMHdr % UnitNum, ITempData, NumCrayWords, &
                                Len_IO, Err_IO )
    Field % RData = REAL( ITempData )

    DEALLOCATE( ITempData )
    Field % Hdr % LBUser1 = 1   ! Field is now real

  ELSE

    WRITE( ErrMessage, '(A24,I2)')   &
                  "Unrecognised data type: ", Field % Hdr % LBUser1
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999

  END IF
END IF

! Perform optional header modifications
IF ( PPHdrMod ) THEN
  IF ( Field % Hdr % MO8Type == 58 ) THEN
   ! Modify MetO8 code for 1.5m max/min/instantaneous temperature
    IF      ( Field % Hdr % LBPROC == 4096 ) THEN
      WRITE(6,*) "Changing MetO8 code from 58 to 157"
      Field % Hdr % MO8Type = 157  ! Minimum temperature
    ELSE IF ( Field % Hdr % LBPROC == 8192 ) THEN
      WRITE(6,*) "Changing MetO8 code from 58 to 156"
      Field % Hdr % MO8Type = 156  ! Maximum temperature
    END IF
  END IF
  IF ( Field % Hdr % LBTIM /= 11 ) THEN
   ! Modify date/time for time-mean / accumulated fields
    IF ( Field % Hdr % ValidYear > 0 ) THEN
      WRITE(6,*) "Modifying date / time header"
     ! Move validity time to correct position in header
      Field % Hdr % ValidYear  = Field % Hdr % DataYear
      Field % Hdr % ValidMonth = Field % Hdr % DataMonth
      Field % Hdr % ValidDate  = Field % Hdr % DataDate
      Field % Hdr % ValidHour  = Field % Hdr % DataHour
      Field % Hdr % ValidMin   = Field % Hdr % DataMin
      Field % Hdr % ValidSec   = Field % Hdr % DataSec
     ! Calculate data time : (validity time - FCRange)
      Year  = Field % Hdr % DataYear
      Month = Field % Hdr % DataMonth
      Date  = Field % Hdr % DataDate
      Hour  = Field % Hdr % DataHour
      Min   = Field % Hdr % DataMin
      Sec   = Field % Hdr % DataSec
      FCRangeSecs = Field % Hdr % FCRange * 3600
! DEPENDS ON: time2sec
      CALL TIME2SEC( Year, Month, Date, Hour, Min, Sec,  &
                     0, 0, EndTimeDays, EndTimeSecs, Lcal360 )
! DEPENDS ON: time_df
      CALL TIME_DF( EndTimeDays,  EndTimeSecs, 0, FCRangeSecs,  &
                    DataTimeDays, DataTimeSecs )
! DEPENDS ON: sec2time
      CALL SEC2TIME( 0, 0, DataTimeDays, DataTimeSecs,   &
                     Year, Month, Date, Hour, Min, Sec, DayNo, Lcal360 )
     ! Put data time into correct position in header
      Field % Hdr % DataYear  = Year
      Field % Hdr % DataMonth = Month
      Field % Hdr % DataDate  = Date
      Field % Hdr % DataHour  = Hour
      Field % Hdr % DataMin   = Min
      Field % Hdr % DataSec   = Sec
    ELSE
      ErrorStatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus,  &
                    "Missing temporal data - cannot recalculate data time" )
    END IF
  END IF
END IF

9999 CONTINUE

END SUBROUTINE ReadFld

