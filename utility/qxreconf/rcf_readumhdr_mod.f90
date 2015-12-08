! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Read in the header from a UM dump.

Module Rcf_ReadUMhdr_Mod

!  Subroutine Rcf_ReadUMhdr - read a header from the dump
!
! Description:
!   Read in a model header from a UM dump.
!
! Method:
!   Call READ_FLH (UM routine) to read in Fixed-length Header.
!   Call Rcf_Allochdr to allocate array sizes for whole header, and
!   Call READHEAD (UM routine) to read in whole header.
!
!  Based on VAR code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE Rcf_ReadUMhdr (  UMhdr )

Use Rcf_UMhead_Mod, Only : &
    LenFixHd,          &
    UM_Header_type

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_AllocHdr_Mod, Only : &
    Rcf_AllocHdr

USE IO
IMPLICIT NONE

! Subroutine arguments
TYPE (UM_header_type), INTENT (INOUT) ::  UMhdr ! UM header from dump

! Local constants
CHARACTER (LEN=*), PARAMETER :: RoutineName='Rcf_ReadUMhdr'

! Local variables
INTEGER             ::  ReturnCode   ! Return code from UM routines
INTEGER             ::  WordAddress  ! Position on file, used in SETPOS

INTEGER             ::  FFLookup
INTEGER             ::  Len_IO
INTEGER, ALLOCATABLE::  Lookup(:)
REAL                ::  A_IO
CHARACTER (Len=80)  ::  Cmessage

! Function & Subroutine interfaces:
EXTERNAL READ_FLH, READHEAD   ! UM routines (Fortran 77)


!-----------------------------------------------------------------------
!  Section 0.1:  Initialisations
!-----------------------------------------------------------------------

CMessage   = ' '
ReturnCode = 0


!-----------------------------------------------------------------------
!  Section 1.    Allocate & read in Fixed-Length Header from Model
!                dump file
!-----------------------------------------------------------------------

ALLOCATE (UMhdr % FixHd(LenFixHd))

WordAddress = 0

CALL SETPOS (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /= 0) THEN
  Cmessage = 'SETPOS failure'
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF

! DEPENDS ON: read_flh
CALL READ_FLH (                   &
  UMhdr % UnitNum,  &   ! in
  UMhdr % FixHd,    &   ! out
  LenFixHd,         &   ! in
  ReturnCode,       &   ! out
  CMessage )            ! out

IF (ReturnCode /= 0 ) THEN
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF


!-----------------------------------------------------------------------
!  Section 2.    Allocate space for the UM header from the Model
!                dump file
!-----------------------------------------------------------------------

! First need to count Lookups if dealing with a FieldsFile
IF ( Umhdr % FixHd( 5 ) == 3 ) Then        ! FieldsFile

  ! Check if Lookups exist first
  FFLookup = 0
  IF ( Umhdr % FixHd( 150 ) > 0 ) THEN

    Call SETPOS( UMhdr % UnitNum, UMhdr % FixHd( 150 ), ReturnCode )

    IF ( ReturnCode /= 0 ) THEN
      Cmessage = 'SETPOS failure'
      CALL Ereport( RoutineName, ReturnCode, Cmessage )
    END IF

    ALLOCATE( Lookup( UMhdr % FixHd( 151 ) ) )
    Lookup(:) = 0

    ! Count the Lookups until the end ...
    DO

      CALL BUFFIN( UMhdr % UnitNum, Lookup, UMhdr % FixHd( 151 ), &
              Len_IO, a_IO )

      IF ( Lookup(1) == -99 ) THEN
        EXIT
      ELSE
        FFLookup = FFLookup + 1
      ENDIF

      IF ( FFLookup > UMhdr % Fixhd(152) ) THEN
        Cmessage = 'FieldsFile has more Lookups than specified in fixed header'
        ReturnCode = 1
        CALL Ereport( RoutineName, ReturnCode, Cmessage )
      END IF
    END DO

    DEALLOCATE( Lookup )
  END IF

  CALL Rcf_AllocHdr( UMhdr, FFLookup )
ELSE
  CALL Rcf_AllocHdr ( UMhdr )        ! inout
END IF


!-----------------------------------------------------------------------
!  Section 3.    Read in complete UM header from the Model dump file
!-----------------------------------------------------------------------

WordAddress = 0

CALL SETPOS (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /= 0) THEN
  Cmessage = 'SETPOS failure'
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF

! DEPENDS ON: readhead
CALL READHEAD (                                                &
  UMhdr % UnitNum,                                             &
  UMhdr % FixHd,     LenFixHd,                                 &
  UMhdr % IntC,      UMhdr % LenIntC,                          &
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
  UMhdr % Lookup,    UMhdr % Len1Lookup,  UMhdr % Len2Lookup,  &
  UMhdr % LenData,                                             &
  UMhdr % StartData,                                           &
  ReturnCode,        CMessage )

IF (ReturnCode /= 0) THEN
  Call Ereport( RoutineName, ReturnCode, Cmessage )
END IF

WordAddress = 0

CALL SETPOS (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /=  0) THEN
  Cmessage = 'SETPOS failure'
  Call Ereport( RoutineName, ReturnCode, Cmessage )
END IF

END SUBROUTINE  Rcf_ReadUMhdr

END MODULE Rcf_ReadUMhdr_Mod
