! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for copying fields from one UM fieldsfile to another

SUBROUTINE CopyAll( NumFlds,     &  ! in
                    MaxFlds,     &  ! in
                    STCode,      &  ! in
                    MO8Level,    &  ! in
                    MO8LevelLwr, &  ! in
                    MO8LevelUpr, &  ! in
                    FCTime,      &  ! in
                    FCTimeStart, &  ! in
                    FCTimeFreq,  &  ! in
                    FCTimeEnd,   &  ! in
                    LBProc,      &  ! in
                    MinsPastHr,  &  ! in
                    MinsDif,     &  ! in
                    Factor,      &  ! in
                    Store,       &  ! in
                    PPHdrMod,    &  ! in
                    UMHdr_in,    &  ! in
                    UMHdr_out,   &  ! inout
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
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode     (MaxFlds)
INTEGER, INTENT(INOUT) :: MO8Level    (MaxFlds)
INTEGER, INTENT(IN)    :: MO8LevelLwr (MaxFlds)
INTEGER, INTENT(IN)    :: MO8LevelUpr (MaxFlds)
INTEGER, INTENT(INOUT) :: FCTime      (MaxFlds)
INTEGER, INTENT(IN)    :: FCTimeStart (MaxFlds)
INTEGER, INTENT(IN)    :: FCTimeFreq  (MaxFlds)
INTEGER, INTENT(IN) :: FCTimeEnd  (MaxFlds)
INTEGER, INTENT(IN) :: LBProc     (MaxFlds)
INTEGER, INTENT(IN) :: MinsPastHr (MaxFlds)
INTEGER, INTENT(IN) :: MinsDif    (MaxFlds)
INTEGER, INTENT(IN) :: Store      (MaxFlds)
REAL,    INTENT(IN) :: Factor     (MaxFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr_in

TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr_out
TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CopyAll"

! Local Variables:
INTEGER :: NumLookups
CHARACTER(LEN=80) :: ErrMessage

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Allocate space for field
NumLookups = UMHdr_in % NumFlds
WRITE(6,*) "Copying ", NumLookups, " Fields"

! Check there is space in the lookup table
IF ( NumLookups + UMHdr_out % NumFlds > UMHdr_out % Len2Lookup ) THEN
  WRITE(6,'(A45)') "Insufficient space allocated for Lookup Table"
  WRITE(6,'(A60)') "Increase the value of MaxFldsOut in maxfields_mod to correct"
  ErrorStatus = StatusFatal

  CALL EReport( RoutineName, ErrorStatus, &
               "Insufficient space allocated for Lookup Table" )
END IF

! By specifying more fields than max we are signalling to copyflds to just copy
! everything.
CALL CopyFlds(NumFlds+MaxFlds+1, MaxFlds, StCode,                      &
              MO8Level, MO8LevelLwr, MO8LevelUpr,                      &
              FCTime, FCTimeStart, FCTimeFreq, FCTimeEnd,              &
              LBProc, MinsPastHr, MinsDif, Factor, Store, PPHdrMod,    &
              UMHdr_in, UMHdr_out, Fields,ErrorStatus)

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CopyAll

