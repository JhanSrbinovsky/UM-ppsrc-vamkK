! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Read the HEADERS namelist

Module rcf_readnl_headers_mod

!  Subroutine Rcf_Readnl_Headers - Read the HEADERS namelist
!
! Description:
!   Read the HEADERS namelist for header overrides.
!
! Method:
!   Variables read into Rcf_Headers_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine rcf_readnl_headers( nft )

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, Only : &
    mype

Use rcf_headers_mod, Only : &
    Fixhd,                  &
    Inthd,                  &
    Relhd,                  &
    headers

Implicit None

! Arguments
Integer, Intent (In)    :: nft

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

! Set defaults = MDI
FixHd(:) = IMDI
IntHd(:) = IMDI
RelHd(:) = RMDI

! Read Namelist
Read( Unit=nft, Nml = headers )

! Write out namelist for diagnostic - Just during development use
! a low output level
If (PrintStatus >= PrStatus_Oper) Then
  If (mype == 0) Then
    Write (Unit=6, Nml=headers)
  End If
End If

Return
End Subroutine rcf_readnl_headers
End Module rcf_readnl_headers_mod
