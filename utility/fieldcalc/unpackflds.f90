! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for packing and unpacking

SUBROUTINE UnPackFlds( PackType,        &  ! in
                       NumCols,         &  ! in
                       NumRows,         &  ! in
                       Field,           &  ! inout
                       ErrorStatus )       ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  LenWord,                &
  PP_Field_type,          &
  PP_Header_type,         &
  LSMField
USE mask_compression, ONLY: expand_from_mask
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: PackType
INTEGER, INTENT(IN) :: NumCols
INTEGER, INTENT(IN) :: NumRows

TYPE(PP_Field_type), INTENT(INOUT) :: Field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "UnPackFlds"
LOGICAL, PARAMETER :: Expand = .FALSE.

! Local Variables:
INTEGER :: FldSize
INTEGER :: idum                     ! Dummy integer
INTEGER :: NumWords                 ! Size of compressed field
CHARACTER(LEN=80) :: ErrMessage     ! Error message returned by COEX

! Local Arrays:
REAL :: WorkArray(NumCols*NumRows,1) ! array used for un_packing

LOGICAL :: PackMask        ! Do we have some sort of mask
LOGICAL :: LandMask        ! Is it a landmask?
INTEGER :: Num_Land_Points ! Number of land points

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

NumWords = Field % Hdr % LBLREC
FldSize  = NumCols * NumRows
PackMask = MOD(Field % Hdr % LBPack/10,10) == 2
LandMask = MOD(Field % Hdr % LBPack/100,10) == 1

IF ( PackType == 1 ) THEN
 ! WGDOS packing
! DEPENDS ON: coex
  CALL COEX( WorkArray,  FldSize,  Field % RData, &
             NumWords,   NumCols,  NumRows,       &
             idum,       idum,     Expand,        &
             Field % Hdr % BMDI,   LenWord,       &
             ErrorStatus, ErrMessage)
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, &
                  "COEX: " // ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
ELSE IF ( PackType == 2 ) THEN
! 32 bit packing (in dumps)
  WorkArray = Field % RData
! DEPENDS ON: expand32b
  CALL expand32b( numwords , WorkArray,  Field % Hdr % LBuser7 )

ELSE IF ( PackType == 3 ) THEN
 !  GRIB packing
  !CALL DEGRIB( Field % RData, WorkArray, FldSize, NumWords,    &
  !             Field % Hdr, Field % Hdr % BMDI, FldSize, LenWord)
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "GRIB packing not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

ELSE IF ( PackType == 4 ) THEN
 ! Run length encoded data
  !CALL RUNLEN_DECODE( WorkArray, FldSize, Field % RData, NumWords,  &
  !                    Field % Hdr % BMDI, ErrorStatus, ErrMessage)
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "RUNLEN packing not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

ELSE
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "UNPACK - PackType not recognised." )
  ErrorStatus = StatusWarning
  GO TO 9999

ENDIF

Field % RData = WorkArray

! We now unpack land mask

IF (PackMask)  THEN
  ! Check we have LSM field available
  IF (.NOT. ASSOCIATED(LSMField % Rdata)) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, &
                  "LSM required to unpack land packed fields." )
    ErrorStatus = StatusWarning
  ELSE
    IF (LandMask) THEN
      CALL expand_from_mask(Field % RData, WorkArray, LSMField % RData /= 0.0, &
                            FldSize, Num_Land_Points)
      Field % Hdr % LBPACK = Field % Hdr % LBPACK - 100 ! Data unpacked land
    ELSE
      ! For sea-mask packing we can just reverse the land packing mask.
      CALL expand_from_mask(Field % RData, WorkArray, LSMField % RData == 0.0, &
                            FldSize, Num_Land_Points)
      Field % Hdr % LBPACK = Field % Hdr % LBPACK - 200 ! Data unpacked sea
    END IF

    Field % Hdr % LBPACK = Field % Hdr % LBPACK - 20 ! Data unpacked mask
  END IF
END IF

Field % Hdr % LBPACK = Field % Hdr % LBPACK - PackType ! Data unpacked
Field % Hdr % LBLREC = FldSize
Field % Hdr % LBUser1 = 1  ! Data type is now real
Field % Hdr % BACC = -99.0


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE UnPackFlds

