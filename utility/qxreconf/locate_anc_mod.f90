! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Locate position of record for Ancillary field

Module Locate_Anc_Mod

!  Subroutine Locate_Anc_Field  - locates record for ancillary field
!  Subroutine Locate_Anc_File  -  locates record for ancillary file
!
! Description:
!    Locates ANCILmaster records for ancillary fields and files
!
! Method:
!    Given the Ancillary Reference/File Number, it searches through the
!    ancillary records to locate the record. If it is not found,
!    it will abort.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Locate_Anc_Field ( AncRefNo, Irec )

Use Ancil_mod, Only :     &
    anc_record,           &
    ancrecs

Use Ereport_mod, Only : &
    Ereport

IMPLICIT NONE

! Arguments
Integer :: AncRefNo  ! Ancil Ref No to search for
Integer :: Irec      ! Position of Ancil Ref No in Anc_Record

! Local variables
Integer                      :: i
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Locate_Anc_Field'
Character (Len=80)           :: Cmessage

ErrorStatus = 0
Irec = 0

Do i = 1, ancRecs

  If (AncRefNo == anc_record(i) % ancil_ref_number) then
    Irec = i
    exit
  Endif

Enddo

If (Irec == 0) Then
  Write (Cmessage, '(A, I4)') &
  'Cannot find ancillary record for Ancillary Reference No ',AncRefNo
  ErrorStatus = 1
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Return
End Subroutine Locate_Anc_Field

Subroutine Locate_Anc_File ( AncFileNo, Irec )

Use Ancil_mod, Only :     &
    ancFiles,             &
    anc_file

Use Ereport_mod, Only : &
    Ereport

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Arguments
Integer :: AncFileNo  ! Ancil File No to search for
Integer :: Irec       ! Position of Ancil File No in Anc_Record

! Local variables
Integer                      :: i
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Locate_Anc_Field'
Character (Len=80)           :: Cmessage

ErrorStatus = 0
Irec = 0

Do i = 1, ancFiles

  If (AncFileNo == anc_file(i) % anc_file_number) then
    Irec = i
    exit
  Endif

Enddo

If (Irec == 0) Then
  Write (Cmessage, '(A, I4)') &
  'Cannot find ancillary record for Ancillary File No ',AncFileNo
  ErrorStatus = 1
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Return
End Subroutine Locate_Anc_File

End Module Locate_Anc_Mod
