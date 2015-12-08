! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to set up a new UM header in Fieldcalc

SUBROUTINE New_UMHdr( UMHdr,        &  ! in
                      MaxFlds,      &  ! in
                      NewUMHdr,     &  ! inout
                      ErrorStatus )    ! inout

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

USE IO
USE IO_Mod, ONLY: &
  LenFixHd,       &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN)    :: UMHdr
INTEGER, INTENT(IN) :: MaxFlds

TYPE(UM_Header_type), INTENT(INOUT) :: NewUMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "New_UMHdr"
CHARACTER(LEN=*), PARAMETER :: PPFormat = "p"

! Local Variables:
INTEGER :: WordAddress
CHARACTER(LEN=80) :: ErrMessage
LOGICAL :: l_append
CHARACTER(LEN=2)  :: c_append

! End of header --------------------------------------------------------

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Call FORT_GET_ENV('APPENDFILE',10,c_append,1,errorstatus)

l_append = .FALSE.

IF (errorstatus == 0) THEN
  IF (c_append(1:1) == '1') THEN
    WRITE(6,'(A)') 'INFO: Appending to ouput file.'
    l_append = .TRUE.
  END IF
END IF

IF (l_append) THEN
! DEPENDS ON: read_umhdr
  CALL read_umhdr(newUmhdr, l_append, errorstatus)
  IF (errorstatus /= statusok) THEN
    errorstatus = statusfatal
  END IF
! No goto end since read_umhdr should have setuop everything.
  GOTO 9999
END IF


CALL File_Open ( NewUMHdr % UnitNum,                     &
                 NewUMHdr % FileNameEnv,                 &
                 LEN_TRIM(NewUMHdr % FileNameEnv),       &
                 1, 0, ErrorStatus)

IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, &
                "Failed to open output file" )
  ErrorStatus = StatusFatal
  GO TO 9999
END IF

!----------------------------------------------------
! Use header provided as start point for new header
!----------------------------------------------------
WordAddress = 0   ! Start of file
CALL SETPOS ( NewUMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! Setup new length of lookup.
NewUMHdr % Len2Lookup = MaxFlds

! DEPENDS ON: init_pp
CALL INIT_PP ( NewUMHdr % UnitNum,     &
               PPFormat,               &
               UMHdr % FixHd(151),     &
               NewUMHdr % Len2Lookup,  &
               UMHdr % FixHD,          &
               UMHdr % IntC,           &
               UMHdr % RealC,          &
               UMHdr % LevDepC,        &
               UMHdr % RowDepC,        &
               UMHdr % ColDepC,        &
               UMHdr % LenIntC,        &
               UMHdr % LenRealC,       &
               UMHdr % Len1LevDepC,    &
               UMHdr % Len2LevDepC,    &
               UMHdr % Len1RowDepC,    &
               UMHdr % Len2RowDepC,    &
               UMHdr % Len1ColDepC,    &
               UMHdr % Len2ColDepC,    &
               UMHdr % LenIntC,        &
               UMHdr % LenRealC,       &
               ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-------------------------------------
! Read the header written by INIT_PP
!-------------------------------------
ALLOCATE( NewUMHdr % FixHd(LenFixHd) )

WordAddress = 0   ! Start of file
CALL SETPOS ( NewUMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: read_flh
CALL READ_FLH ( NewUMHdr % UnitNum, & ! in
                NewUMHdr % FixHd,   & ! out
                LenFixHd,           & ! in
                ErrorStatus,        & ! out
                ErrMessage )          ! out
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: setup_umhdr
CALL Setup_UMHdr( NewUMHdr )
NewUMHdr % NumFlds = 0
NewUMHdr % IntC    = UMHdr % IntC
NewUMHdr % RealC   = UMHdr % RealC
NewUMHdr % LevDepC = UMHdr % LevDepC
NewUMHdr % RowDepC = UMHdr % RowDepC
NewUMHdr % ColDepC = UMHdr % ColDepC
ALLOCATE( NewUMHdr % Lookup(NewUMHdr % Len2Lookup) )

9999 CONTINUE

END SUBROUTINE New_UMHdr

