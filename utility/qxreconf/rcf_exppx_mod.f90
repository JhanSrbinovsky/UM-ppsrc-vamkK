! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Extract a STASHmaster record

Module Rcf_Exppx_Mod

!  Subroutine Rcf_Exppx   - extract a STASHmaster record.
!
! Description:
!   Given a model, section and item this routine returns a pointer
!   to the specified STASHmaster record
!
! Method:
!   The relevant row in the record array (STM_record) is found
!   by the ppxptr index.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Function Rcf_Exppx( Im_ident, section, item, NoFindArg, stash_in_arg )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Ppx_Info_mod, Only : &
    STM_record_type,     &
    STM_record,          &
    STM_record_in,       &
    ppxptr,              &
    ppxptr_in

IMPLICIT NONE
! Function arguments:
Integer, Intent(In) :: Im_ident    ! Internal model identifier
Integer, Intent(In) :: section     ! STASH section no.
Integer, Intent(In) :: item        ! STASH item no.
Logical, Optional, Intent(In) :: NoFindArg ! Warning with no entry found
Logical, Optional, Intent(In) :: stash_in_arg ! using input STASHmaster

! Function value (out)
Type (STM_record_type), Pointer :: Rcf_Exppx


! Local variables
Character (Len=*), Parameter :: RoutineName = 'Rcf_Exppx'
Character (Len=80)           :: Cmessage
Integer                      :: row         ! Row no. in PPXI array
Integer                      :: ErrorStatus !+ve = fatal error
Integer, Pointer             :: ppxptr_ptr(:,:,:) ! Pointer for ppxptr/_in
Logical                      :: NoFind  !local value of NoFindArg
Logical                      :: stash_in  ! local value of stash_in_arg
Type(STM_record_type), Pointer :: STM_record_ptr(:) ! Pointer for STM_record/_in

If ( Present (NoFindArg)) Then
  NoFind = NoFindArg
Else
  NoFind = .FALSE.
End If
If ( Present (stash_in_arg) ) Then
  stash_in = stash_in_arg
Else
  stash_in = .FALSE.
End If

! Set up pointers dependent on input/output stash
If ( stash_in ) Then
  ppxptr_ptr     => ppxptr_in
  STM_record_ptr => STM_record_in
Else
  ppxptr_ptr     => ppxptr
  STM_record_ptr => STM_record
End If

Nullify( Rcf_Exppx )

ErrorStatus = 0
If ( Im_ident <= 0 .OR. section < 0 .OR. item <= 0) Then
  If (Im_ident <= 0) Write(6,*) 'Exppxi: INVALID Im_ident'
  If (section  <  0) Write(6,*) 'Exppxi: INVALID SECTION NO.'
  If (item     <=  0) Write(6,*) 'Exppxi: INVALID ITEM NO.'
  Write(6,*) &
      'Im_ident ',Im_ident,' section ',section,' item ',item
  ErrorStatus=1
  Cmessage='ERROR Exppx: INVALID STASH RECORD ID'

  Call Ereport( RoutineName, ErrorStatus, Cmessage )

Else

  ! Obtain row no. in PPXI array
  row = ppxptr_ptr(Im_ident,section,item)

  ! Obtain required data value
  If ( row > 0 ) Then
    Rcf_Exppx => STM_record_ptr( row )
  Else

    If ( NoFind ) Then
      Write (Cmessage, '(A, I5, A, I3, A, I2)')                        &
            'STASH entry not defined : item ', item,                   &
            ' section ', section,' model ',  Im_ident
      ErrorStatus = -2
    Else
      Write (Cmessage, '(A, I5, A, I3, A, I2, A)')                     &
            'Cant find required STASH item ', item,                    &
            ' section ', section,' model ',  Im_ident, ' in STASHmaster'
      ErrorStatus = 2
    End If

    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
End If

NULLIFY( ppxptr_ptr )
NULLIFY( STM_record_ptr )

Return
End Function Rcf_Exppx

End Module Rcf_Exppx_Mod
