! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for dynamic allocation of Fortran Unit Number

Module Rcf_FortranIO_Mod

!  Subroutine Rcf_Get_Unit     - obtains a unit number for IO
!  Subroutine Rcf_Free_Unit    - frees up a unit number again
!
! Description:
!   The module contains data and routines to manage the dynamic
!   allocation and freeing of unit numbers for Fortran IO.
!
! Method:
!   An array UnitsInUse holds a list of all unit numbers already
!   allocated through the modules subroutines. This is based on
!   code used by the VAR group.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   11/10/01   Added Max_Filename_Len parameter.  R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Integer, Parameter      :: LB               =  10
Integer, Parameter      :: UB               =  64
Integer, Parameter      :: Max_Filename_Len = 256

Integer, Private        :: i

Logical, Save           :: UnitsInUse( LB : UB ) &
                           = (/(.FALSE. , i = LB, UB)/)

Contains

!-------------------------------------------------------------------
  Subroutine Rcf_Get_Unit( nft )

  Use Ereport_Mod, Only :  &
      Ereport

IMPLICIT NONE

  ! Arguments
  Integer, Intent(Out)         :: nft

  ! Local Variables
  Character (Len=*), Parameter :: RoutineName = 'FortranIO'
  Integer                      :: ErrorStatus
  Character (Len=80)           :: Cmessage

  nft = LB

  Do While ( UnitsInUse( nft ) .And. nft <= UB)
   nft = nft + 1
  End Do

  If (nft == UB) Then
    Cmessage = 'Cannot obtain free Fortran file-unit!'
    ErrorStatus = 1
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  UnitsInUse( nft ) = .TRUE.

  Return
  End Subroutine Rcf_Get_Unit
!--------------------------------------------------------------------

  Subroutine Rcf_Free_Unit( nft )

    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE

  ! Arguments
  Integer, Intent(In)     :: nft

  UnitsInUse(nft) = .FALSE.

  Return
  End Subroutine Rcf_Free_Unit

!--------------------------------------------------------------------

End Module Rcf_FortranIO_Mod
