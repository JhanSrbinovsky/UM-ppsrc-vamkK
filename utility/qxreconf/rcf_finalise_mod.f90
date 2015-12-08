! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ finalise the reconfiguration - ties up loose ends

Module Rcf_Finalise_Mod

!   Subroutine Rcf_Finalise - final tidyup tasks
!
! Description:
! This module performs tidy up tasks including closing files and
! calling Timer for output.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Finalise( hdr_in, hdr_out )

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Free_Unit

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars, Only : &
    mype

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid

Use Rcf_Recon_Mod, Only : &
    Grib

USE PrintStatus_mod, Only : &
    LTimer

USE IO
USE model_file, Only : model_file_close

IMPLICIT NONE

! Arguments
Type (um_header_type), Intent(InOut)  :: hdr_in
Type (um_header_type), Intent(InOut)  :: hdr_out

! Local variables
Integer                               :: err
Integer                               :: ErrorStatus
Character (Len=*), Parameter          :: RoutineName='Rcf_Finalise'
Character (Len=80)                    :: Cmessage

!-------------------------------------------------------------------
! Close open files
!-------------------------------------------------------------------
err = 0
! Output dump

Call Model_File_Close( hdr_out % UnitNum, 'ASTART', 6, 0, 0, err)

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Close Output Dump'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Rcf_Free_Unit( hdr_out % UnitNum )

! Input dump
If (Input_Grid % Rotated .AND. .NOT. Grib ) Then

  Call File_Close( hdr_in % UnitNum, 'RECONTMP', 8, 0, 0, err)
Else

  Call File_Close( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)
End If

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Close Input Dump'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Rcf_Free_Unit( hdr_in % UnitNum )

!-------------------------------------------------------------------
! Produce Timer output
!-------------------------------------------------------------------
If (LTimer) Then
! DEPENDS ON: timer
  Call Timer("Reconfigure",2)
End If

Return
End Subroutine Rcf_Finalise
End Module Rcf_Finalise_Mod
