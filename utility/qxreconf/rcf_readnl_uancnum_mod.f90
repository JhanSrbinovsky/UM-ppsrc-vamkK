! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Read the UANCNUM namelist

MODULE Rcf_readnl_Uancnum_Mod

!  Subroutine Rcf_Readnl_uancnum - reads the uancnum namelist
!
! Description:
!    Reads the uancnum namelist for user ancilmaster processing.
!
! Method:
!    Reads the variables for the Ancil_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CONTAINS

SUBROUTINE Rcf_readnl_Uancnum( nft )

USE Ancil_mod, ONLY :    &
    n_uancil,            &
    uancfils

USE PrintStatus_mod, ONLY : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, ONLY : &
    mype

IMPLICIT NONE

INTEGER :: nft               ! Unit number
INTEGER :: i                 ! Looper

! Initialise and read UANCNUM namelist : number of user
!   ancilmaster files and total no. of user ancil records

! Initialisation
n_uancil     = 0
DO i = 1,20
  uancfils(i)='        '
END DO

! For ROSE - disabling this code as the name list isn't suitable for
! ROSE in its present form. Have commented out the code rather than deleting it,
! in case there is a need to reactivate it in the future - although the relation
! between name list and user interface would have to be re-designed.

! Read namelist
! Read( Unit = nft, Nml = UANCNUM )

! If (PrintStatus >= PrStatus_Oper .AND. mype == 0) then
!   Write ( unit = 6 , Nml = UANCNUM )
! Endif

RETURN
END SUBROUTINE Rcf_readnl_Uancnum
END MODULE Rcf_readnl_Uancnum_Mod
