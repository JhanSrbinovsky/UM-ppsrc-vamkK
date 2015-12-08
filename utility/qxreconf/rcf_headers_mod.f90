! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Define headers namelist

Module rcf_headers_mod

! Description:
!    Define headers namelist
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

Integer, Save :: FixHd(256)
Integer, Save :: IntHd(100)
Real, Save    :: RelHd(100)

Namelist /Headers/ FixHd,IntHd,RelHd
End Module rcf_headers_mod
