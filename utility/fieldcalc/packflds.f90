! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for packing and unpacking


!=======================================================================
SUBROUTINE PackFlds( NumFlds,     &  ! in
                     MaxFlds,     &  ! in
                     PackType,    &  ! in
                     Source,      &  ! in
                     PackAcc,     &  ! in
                     Fields,      &  ! inout
                     ErrorStatus )   ! inout

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
  PP_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
CHARACTER(LEN=*), INTENT(IN) :: PackType
INTEGER, INTENT(IN) :: Source (NumFlds)
REAL,    INTENT(IN) :: PackAcc(NumFlds)

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "PackFlds"
LOGICAL, PARAMETER :: Compress = .TRUE.

! Local variables:
INTEGER :: i, ifld
INTEGER :: PackCode
INTEGER :: NumWords
INTEGER :: Num32BitWds
REAL, POINTER :: CompData(:,:)
CHARACTER(LEN=80) :: ErrMessage

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF( PackType == "WGDOS" ) THEN
  PackCode = 1
  DO i=1, NumFlds
    ifld = Source(i)
    Fields(ifld) % Hdr % BACC = PackAcc(i)
    ALLOCATE( CompData( Fields(ifld) % Hdr % LBLREC, 1 ) )
! DEPENDS ON: coex
    CALL COEX( Fields(ifld) % RData,         &
               Fields(ifld) % Hdr % LBLREC,  &
               CompData,                     &
               Fields(ifld) % Hdr % LBLREC,  &
               Fields(ifld) % Hdr % NumCols, &
               Fields(ifld) % Hdr % NumRows, &
               Num32BitWds,                  &
               INT(PackAcc(i)),              &
               Compress,                     &
               Fields(ifld) % Hdr % BMDI,    &
               LenWord,                      &
               ErrorStatus, ErrMessage )
    IF( ErrorStatus /= StatusOK ) THEN
      ErrorStatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, ErrMessage )
      ErrorStatus = StatusWarning
      GO TO 9999
    END IF
    NumWords = ( Num32BitWds-1 + (LenWord/32) ) * 32/LenWord
    DEALLOCATE( Fields(ifld) % RData )
    ALLOCATE( Fields(ifld) % RData( NumWords, 1 ) )
    Fields(ifld) % Hdr % LBLREC = NumWords
    Fields(ifld) % RData( 1:NumWords, 1 ) = CompData( 1:NumWords, 1 )
    DEALLOCATE( CompData )
    Fields(ifld) % Hdr % LBPACK = Fields(ifld) % Hdr % LBPACK + 1
  END DO
ELSE IF( PackType == "CRAY32" ) THEN
  WRITE(6,*) "CRAY32 Packing not supported - leaving unpacked"
ELSE IF( PackType == "GRIB" ) THEN
  WRITE(6,*) "GRIB Packing not supported - leaving unpacked"
ELSE IF( PackType == "RUNLEN" ) THEN
  WRITE(6,*) "RUNLEN Packing not supported - leaving unpacked"
ELSE
  WRITE(6,*) "PackType not recognised - leaving unpacked"
END IF

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE PackFlds

!=======================================================================
