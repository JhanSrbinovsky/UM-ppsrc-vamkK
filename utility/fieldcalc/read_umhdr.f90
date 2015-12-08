! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to read a UM Fieldsfile header

SUBROUTINE Read_UMHdr ( UMHdr,        &  ! inout
                        l_append,     &  ! in
                        ErrorStatus )    ! inout

! Description:
!   This is a subroutine to open, read and store the contents of a UM
!   fieldsfile header.
!
! Method:
!   The environment variable containing the file name is held within the
!   derived type UM_Header_type.  This is used to open the file.  The
!   fixed length header is read, and UMHdr set up accordingly.  Space is
!   allocated for the lookup table, and the rest of the UM header read
!   in.  Memory is conserved by only allocating enough memory for the
!   exact number of lookup entries.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO
USE IO_Mod, ONLY:         &
  LenFixHd,               &
  UM_Header_type,         &
  buffin_um_lookup
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
LOGICAL, INTENT(IN)    :: l_append
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Read_UMHdr"

! Local Variables:
INTEGER :: i
INTEGER :: WordAddress
INTEGER :: Len_IO
REAL :: Err_IO

INTEGER, POINTER :: Lookup(:,:)

CHARACTER(LEN=80) :: ErrMessage

! End of Header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
ENDIF

IF (l_append) THEN
  CALL File_Open ( UMHdr % UnitNum,                        &
                   UMHdr % FileNameEnv,                    &
                   LEN_TRIM(UMHdr % FileNameEnv),          &
                   1, 0, ErrorStatus)
ELSE
  CALL File_Open ( UMHdr % UnitNum,                        &
                   UMHdr % FileNameEnv,                    &
                   LEN_TRIM(UMHdr % FileNameEnv),          &
                   0, 0, ErrorStatus)
END IF

IF ( ErrorStatus /= StatusOK ) THEN
  ! Make this a non-Fatal error and try and continue.
  ErrorStatus = -1*ABS(ErrorStatus)

  CALL EReport( RoutineName, ErrorStatus, &
                "Failure to open input file" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-----------------------------------------------------------------------
! Allocate and read in Fixed_Length_Header
!-----------------------------------------------------------------------
ALLOCATE (UMHdr % FixHd(LenFixHd))

WordAddress = 0   ! Start of file
CALL SETPOS ( UMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! Is this just to check how much space should be allocated?
! DEPENDS ON: read_flh
CALL READ_FLH ( UMHdr % UnitNum, UMHdr % FixHd, LenFixHd,   &
                ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: setup_umhdr
CALL Setup_UMHdr( UMHdr )
! Need to change the data allocation
ALLOCATE ( Lookup (UMHdr % Len1Lookup, UMHdr % Len2Lookup) )

WordAddress = 0
CALL SETPOS ( UMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: readhead
CALL READHEAD ( UMHdr % UnitNum, UMHdr % FixHd, LenFixHd,           &
       UMHdr % IntC,      UMhdr % LenIntC,                          &
       UMhdr % RealC,     UMhdr % LenRealC,                         &
       UMhdr % LevDepC,   UMhdr % Len1LevDepC, UMhdr % Len2LevDepC, &
       UMhdr % RowDepC,   UMhdr % Len1RowDepC, UMhdr % Len2RowDepC, &
       UMhdr % ColDepC,   UMhdr % Len1ColDepC, UMhdr % Len2ColDepC, &
       UMhdr % FldsOfC,   UMhdr % Len1FldsOfC, UMhdr % Len2FldsOfC, &
       UMhdr % ExtraC,    UMhdr % LenExtraC,                        &
       UMhdr % HistFile,  UMhdr % LenHistFile,                      &
       UMhdr % CompFldI1, UMhdr % LenCompFldI1,                     &
       UMhdr % CompFldI2, UMhdr % LenCompFldI2,                     &
       UMhdr % CompFldI3, UMhdr % LenCompFldI3,                     &
       Lookup, UMhdr % Len1Lookup, UMhdr % Len2Lookup,              &
       UMhdr % LenData,                                             &
       UMHdr % StartData,                                           &
       ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-----------------------------------------------------------------------
! Deallocate blank lookup entries
!-----------------------------------------------------------------------
! Count Lookup entries
UMHdr % NumFlds = 0
DO i = 1,UMHdr % Len2Lookup
  IF ( Lookup(1, UMHdr%NumFlds+1) == -99 ) THEN
    IF (.NOT. l_append) THEN
      UMHdr % Len2Lookup = UMHdr % NumFlds
    END IF
    EXIT  ! Exit loop
  ELSE
    UMHdr % NumFlds = UMHdr % NumFlds+1
  END IF
END DO
DEALLOCATE( Lookup )

ALLOCATE( UMHdr % Lookup( UMHdr % Len2Lookup ) )
CALL buffin_um_lookup(UMHdr, Len_IO, Err_IO )
IF ( Err_IO /= -1.0 ) THEN   ! Non-standard error code

  CALL EReport( RoutineName, ErrorStatus, "Failure in BUFFIN" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

9999 CONTINUE

END SUBROUTINE Read_UMHdr

