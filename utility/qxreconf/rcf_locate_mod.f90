! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Find a field in a field list by stash item number

Module Rcf_Locate_mod

!  Subroutine Rcf_Locate - locates a field in a field array
!
! Description:
! Fortran 90 replacement for locate1, using structures rather
! than a whole bunch of arrays. Will return the position of the
! field in the array if it exists otherwise will die gracefully.
! Optional argument zero_ok allows return of 0 in certain circumstances
!
! Method:
!   Loop through the array until the correct field is found
!   in section 0.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   11/10/01   Limit search to section 0. P.Selwood.
!   6.2   10/11/05   Extend to section 34 as well. R Barnes
!   6.1   02/06/04   Extend to section 33 as well. R Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Locate( section, item, fields, nfields, pos, zero_ok_arg )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Ereport_Mod, Only : &
    Ereport

IMPLICIT NONE

! Arguments
Integer, Intent(In)                :: nfields
Integer, Intent(In)                :: section
Integer, Intent(In)                :: item
Integer, Intent(Out)               :: pos
Type (field_type), Intent(In)      :: fields( nfields )
Logical, Intent(In), Optional      :: zero_ok_arg


! Local variables
Integer                            :: i
Integer                            :: ErrorStatus
Logical                            :: zero_ok
Character (Len=*), Parameter       :: RoutineName='Rcf_Locate'
Character (Len=80)                 :: Cmessage

zero_ok = .FALSE.
If ( Present( zero_ok_arg ) ) Then
  zero_ok = zero_ok_arg
End If

pos = 0

Do i = 1, nfields
  If (fields(i) % stashmaster % item    == item .AND. &
      fields(i) % stashmaster % section == section ) Then
    pos = i
    Exit
  End If
End Do

If (pos == 0 .AND. .NOT. zero_ok) Then
  ErrorStatus = 10
  Write(Cmessage, '(A, I3, I5)') 'Unable to Rcf_Locate field with STASHcode ', &
                                  section, item
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Return
End Subroutine Rcf_Locate
End Module Rcf_Locate_Mod
