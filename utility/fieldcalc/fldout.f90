! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to write fields to a UM fieldsfile


!=======================================================================

SUBROUTINE FldOut( Field,       &  ! inout
                   UMHdr,       &  ! inout
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

USE IO
USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE io_configuration_mod, ONLY: &
    io_field_padding
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(INOUT) :: Field
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "FldOut"

! Local Variables:
INTEGER :: LookupPos
INTEGER :: DataAddress
INTEGER :: DataLen
INTEGER :: Len_IO
REAL :: Err_IO
REAL, Allocatable :: Field_Out (:,:)     !  Field written to disk

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Find and set next point to write data
IF( UMHdr % NumFlds == 0) THEN
  DataAddress = UMHdr % StartData - 1  ! First field
ELSE
  ! Start of prev fld + len of prev fld
  DataAddress = UMHdr % Lookup(UMHdr % NumFlds) % DataPos +  &
                UMHdr % Lookup(UMHdr % NumFlds) % LBNREC
END IF

CALL SETPOS ( UMHdr % UnitNum, DataAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN

  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! Set New Position in lookup - add to end of lookup
UMHdr % NumFlds = UMHdr % NumFlds + 1
LookupPos = UMHdr % NumFlds
DataLen = Field % Hdr % LBLREC

! update lookup for field - assume LBLREC is correct, but LBNREC isn't
! LBLREC is incorrect - it should contain the full grid size
Field % Hdr % DataPos = DataAddress
Field % Hdr % LBUser2 = DataAddress
Field % Hdr % LBNREC  = io_field_padding *              &  ! round up
                      ( (DataLen + io_field_padding-1) / io_field_padding )
UMHdr % Lookup( LookupPos ) = Field % Hdr

! Field % RData contains the exact length of data (DataLen)
! Field_Out is rounded up to the next sector boundary

Allocate (Field_Out ( Field % Hdr % LBNREC, 1 ) )

! Copy data and set padding to next sector boundary to zero

Field_Out (1:DataLen, 1) = Field % RData (1:DataLen, 1)
Field_Out (DataLen+1:Field % Hdr % LBNREC, 1) = 0.0

CALL BUFFOUT( UMHdr % UnitNum, Field_Out,               &
              Field % Hdr % LBNREC, Len_IO, Err_IO )

! If checking the error code returned by BUFFOUT remember Err_IO == -1.0
! means write successful

Deallocate (Field_Out)

9999 CONTINUE

END SUBROUTINE FldOut

