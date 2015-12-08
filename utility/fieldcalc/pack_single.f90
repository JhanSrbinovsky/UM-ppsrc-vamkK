! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for packing and unpacking

SUBROUTINE Pack_Single( PackAcc,     &  ! in
                        PackType,    &  ! in
                        XpndField,   &  ! in
                        CompField,   &  ! inout
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
REAL, INTENT(IN)    :: PackAcc
CHARACTER(LEN=*), INTENT(IN) :: PackType
TYPE(PP_Field_type), INTENT(IN) :: XpndField

TYPE(PP_Field_type), INTENT(INOUT) :: CompField
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Pack_Single"
LOGICAL, PARAMETER :: Compress = .TRUE.

! Local variables:
INTEGER :: PackCode
INTEGER :: NumWords
INTEGER :: Num32BitWds
CHARACTER(LEN=80) :: ErrMessage

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Can only pack REAL data
IF (XpndField % Hdr % LBUser1 /= 1) THEN
  WRITE(6,'(A,I7)') "Cannot pack non-REAL data in STCode ", &
                    XpndField % Hdr % STCode
  errorstatus = StatusWarning
  GO TO 9999
END IF


CompField % Hdr = XpndField % Hdr
IF ( ASSOCIATED( Compfield % RData ) ) THEN
  DEALLOCATE( CompField % RData )
  NULLIFY( CompField % RData )
END IF


IF( TRIM(PackType) == "WGDOS" ) THEN


  PackCode = 1
  CompField % Hdr % BACC = PackAcc
  ALLOCATE( CompField % RData( CompField % Hdr % LBLREC, 1 ) )
! DEPENDS ON: coex
  CALL COEX( XpndField % RData,            &
             XpndField % Hdr % LBLREC,     &
             CompField % RData,            &
             CompField % Hdr % LBLREC,     &
             XpndField % Hdr % NumCols,    &
             XpndField % Hdr % NumRows,    &
             Num32BitWds,                  &
             INT(PackAcc),                 &
             Compress,                     &
             XpndField % Hdr % BMDI,       &
             LenWord,                      &
             ErrorStatus, ErrMessage )
  IF( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
  NumWords = ( Num32BitWds-1 + (LenWord/32) ) * 32/LenWord
  CompField % Hdr % LBLREC = NumWords
  CompField % Hdr % LBPACK = XpndField % Hdr % LBPACK + 1
ELSE IF( TRIM(PackType) == "CRAY32" ) THEN
  WRITE(6,'(A)') "CRAY32 Packing not supported - leaving unpacked"
  ErrorStatus = StatusWarning
ELSE IF( TRIM(PackType) == "GRIB" ) THEN
  WRITE(6,'(A)') "GRIB Packing not supported - leaving unpacked"
  ErrorStatus = StatusWarning
ELSE IF( TRIM(PackType) == "RUNLEN" ) THEN
  WRITE(6,'(A)') "RUNLEN Packing not supported - leaving unpacked"
  ErrorStatus = StatusWarning
ELSE
  WRITE(6,'(A)') "PackType not recognised - leaving unpacked"
  ErrorStatus = StatusWarning
END IF

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Pack_Single

