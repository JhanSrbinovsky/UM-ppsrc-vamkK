! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ read the RECON namelist

Module Rcf_readnl_Recon_Mod

!  Subroutine Rcf_Readnl_Recon - read the RECON namelist
!
! Description:
!   Read the recon namelist and set Output_Grid variables accordingly.
!
! Method:
!   Variables read into the Rcf_Recon_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_readnl_recon( nft )

Use Rcf_Recon_Mod

Use Rcf_Grid_Type_Mod, Only : Output_Grid
USE UM_ParVars, Only   : mype
!USE LBC_Mod

USE nlsizes_namelist_mod, ONLY : &
    model_levels

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, Only : &
    mype

Use Ereport_Mod, Only : &
    Ereport

IMPLICIT NONE

! Arguments
Integer                      :: nft     ! Unit number

! Local variables
Integer                      :: ErrorStatus
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Readnl_Recon'

! Initialisation
Grib             = .FALSE.
grib2ff          = .FALSE.
Var_Recon        = .FALSE.
Uars             = .FALSE.
conv_levels      = 0
w_zero_start     = -1
w_zero_end       = -1
q_min            = 0.0
use_smc_stress   = .FALSE.

! Read namelist
Read( Unit = nft, Nml = RECON )
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write( 6 , Nml = RECON )
End If

! Set variables as required
Output_Grid % conv_levels       = conv_levels

! If either w_zero_start or w_zero_end are -1 then set to
! model_levels (from same namelist) as default value.
If (w_zero_start == -1) Then
  w_zero_start = model_levels
End If
If (w_zero_end == -1) Then
  w_zero_end = model_levels
End If
If (w_zero_start == 0) Then  !  W at surface is always set to zero.
  w_zero_start = 1
End If

! Check validity of w_zero_start and w_zero_end
If ( (w_zero_start > w_zero_end) .or.                             &
     (w_zero_start < 0 .or. w_zero_start > model_levels) .or.     &
     (w_zero_end   < 0 .or. w_zero_end   > model_levels) ) Then
  write (6,*) 'w_zero_start ',w_zero_start,' w_zero_end ',w_zero_end
  ErrorStatus = 10
  write (CMessage, '(A)') 'w_zero_start and/or w_zero_end set incorrectly.'
  Call Ereport ( RoutineName, ErrorStatus, cmessage)
End If

Return
End Subroutine Rcf_readnl_recon
End Module Rcf_readnl_Recon_Mod
