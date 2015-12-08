! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to write lookup table to file

SUBROUTINE WriteLookup( UMHdr,        & ! in
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
  UM_Header_type,         &
  buffout_um_lookup
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN) :: UMHdr

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WriteLookup"

! Local variables:
INTEGER :: Len_IO       ! Length of buffer written by BUFFOUT
INTEGER :: Residue      ! Space for lookup not required

REAL    :: Err_IO       ! BUFFOUT error status

! End of header --------------------------------------------------------

IF ( .NOT. ASSOCIATED( UMHdr % Lookup ) ) THEN
  ErrorStatus = StatusWarning
  CALL EReport( RoutineName, ErrorStatus, "No fields to write." )
  GOTO 9999
END IF


CALL buffout_um_lookup( UMHdr, Len_IO, Err_IO )

IF ( Err_IO /= -1.0 ) THEN      ! BUFFOUT returns -1.0 when OK
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, "Failure in BUFFOUT" )
END IF

! Should flush buffer here?

9999 CONTINUE

END SUBROUTINE WriteLookup

