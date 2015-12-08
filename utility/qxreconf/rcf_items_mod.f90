! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Arrays for the items namelists

Module Rcf_Items_Mod

! Description:
!    Arrays for the items namelists
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

Integer, Allocatable :: Source_Array(:)
Integer, Allocatable :: Sctn_Array(:)
Integer, Allocatable :: Item_Array(:)
Integer, Allocatable :: area_array(:)
Integer, Allocatable :: upas_array(:)
Integer, Allocatable :: upaa_array(:)
Real,    Allocatable :: uprc_array(:)
Character (len=80), Allocatable :: upaf_array(:)

Integer, Save :: Num_Items

End Module Rcf_Items_Mod
