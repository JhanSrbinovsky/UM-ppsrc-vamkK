! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module to allocate and deallocate space for field data

Module Rcf_Alloc_Field_Mod
IMPLICIT NONE
!  Subroutine Rcf_Alloc_Field       - allocates space
!  Subroutine Rcf_Dealloc_Field     - frees up space
!
! Description:
!  Allocates and deallocates space for field data
!
! Method:
!  Space is allocated to data pointers based on the data type from
!  stashmaster
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   19/05/03   Perform correct association check, and
!                    initialize allocated arrays to MDI. T.White
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Alloc_Field( field )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Ereport_mod, Only : &
    Ereport

USE cppxref_mod, ONLY: &
    ppx_type_real,     &
    ppx_type_int,      &
    ppx_type_log

IMPLICIT NONE

! Arguments
Type (field_type), Intent(InOut)   :: field

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local data
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Rcf_Alloc_Field'
Character (Len=80)           :: Cmessage

!---------------------------------------------------------------
! Space allocated is local - however must be of the right `type'
!---------------------------------------------------------------

Select Case (field % stashmaster % data_type)
  Case (ppx_type_real)
    If ( .NOT. Associated( field % data) )  Then
      Allocate( field % data( field % level_size, field % levels ) )
      field % data = RMDI
    End If
    Nullify( field % data_int )
    Nullify( field % data_log )

  Case (ppx_type_int)
    If ( .NOT. Associated( field % data_int) ) Then
      Allocate( field % data_int( field % level_size, field % levels ) )
      field % data_int = IMDI
    End If
    Nullify( field % data )
    Nullify( field % data_log )

  Case (ppx_type_log)
    If ( .NOT. Associated( field % data_log) ) Then
      Allocate( field % data_log( field % level_size, field % levels ) )
      field % data_log = .False.
    End If
    Nullify( field % data )
    Nullify( field % data_int )

  Case Default
    Cmessage = 'Field data space cannot be allocated as datatype is '//&
        'not recognised'
    ErrorStatus = -10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

Return

End Subroutine Rcf_Alloc_Field

!--------------------------------------------------------------------
! Dealloc_field frees up *any* allocated data
!--------------------------------------------------------------------

Subroutine Rcf_Dealloc_Field( field )

Use Rcf_Field_Type_Mod, Only : &
    field_type

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

Type (field_type), Intent(InOut)    :: field

If ( Associated( field % data ) ) Then
  Deallocate( field % data )
  Nullify( field % data )
End If

If ( Associated( field % data_int ) ) Then
  Deallocate( field % data_int )
  Nullify( field % data_int )
End If

If ( Associated( field % data_log ) ) Then
  Deallocate( field % data_log )
  Nullify( field % data_log )
End If

Return

End Subroutine Rcf_Dealloc_Field
End Module Rcf_Alloc_Field_mod
