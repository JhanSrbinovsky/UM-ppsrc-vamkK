! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the USTSNUM namelist

Module Rcf_readnl_Ustsnum_Mod

!  Subroutine Rcf_Readnl_Ustsnum - reads the ustsnum namelist
!
! Description:
!   Reads information for user stashmaster processing.
!
! Method:
!   Read the namelist for variables in rcf_ppx_info_mod
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_readnl_Ustsnum( nft )

Use Rcf_Ppx_Info_mod, Only : &
    n_ustash,                &
    n_ustash_in,             &
    nrecs_ustash,            &
    nrecs_ustash_in,         &
    ustsfils,                &
    ustsfils_in,             &
    ustsnum

USE PrintStatus_mod, ONly : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, Only : &
    mype

Implicit None

Integer nft               ! Unit number
Integer i                 ! Looper

! Open stash control file and read USTSNUM namelist: number of user
!   stash files and total no. of user stash records

! Initialisation
Do I = 1,20
  USTSFILS(I)='        '
  USTSFILS_IN(I)='        '
End Do


! Read namelist
Read( Unit = nft, Nml = USTSNUM )
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write( 6, Nml = USTSNUM )
End If

! For now lets copy variables in namelist USTSNUM to USTSNUM_IN ones 
! but we will have to possibly do something clever if we want ENDGAME <-> ND 
! to have different user stash.
n_ustash_in     = n_ustash
nrecs_ustash_in = nrecs_ustash
DO i = 1,20
  ustsfils_in(i)  = ustsfils(i)
END DO

Return
End Subroutine Rcf_readnl_Ustsnum
End Module Rcf_readnl_Ustsnum_Mod
