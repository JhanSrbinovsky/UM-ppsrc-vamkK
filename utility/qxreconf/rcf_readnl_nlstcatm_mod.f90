! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read the NLSTCATM namelist

MODULE Rcf_readnl_nlstcatm_Mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Nlstcatm - Read the NLSTCATM namelist
!
! Description:
!   Read the NLSTCATM namelist of Atmosphere model control variables.
!
! Method:
!   Reads variables into the Rcf_CntlAtm_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CONTAINS

Subroutine rcf_readnl_nlstcatm (nft)

Use Rcf_CntlAtm_Mod

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, Only : &
    mype

USE um_input_control_mod

IMPLICIT NONE

integer nft

H_Sect = ' '

READ(Unit = nft, Nml = nlstcatm)
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  write(Unit = 6, Nml = nlstcatm)
End If

! For test purposes set the lsp items which have been moved into the run_precip
! namelist.  This is due to the fact that run_precip is not read into the rcf, 
! but the items are needed in TSTMSK.  Set to value contained in N48 1 std 
! job

Return
END SUBROUTINE  Rcf_readnl_nlstcatm
END MODULE Rcf_readnl_nlstcatm_Mod
