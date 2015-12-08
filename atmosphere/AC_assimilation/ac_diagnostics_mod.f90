! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Diagnostics necessary for AC assimilation scheme.


Module ac_diagnostics_Mod

! Description:
! This module is used to pass diagnostics from physics routines
! to the AC assimilation scheme.
!
! Method:
! By allocating and filling the arrays in routines control/top_level/
! diagnostics_adv/_conv/_lscld/_lsrain, and using them in ac_ctl.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! Diagnostic Arrays

Real, Allocatable  ::  LSRR(:)  ! large scale rain rate
Real, Allocatable  ::  LSSR(:)  ! large scale snow rate
Real, Allocatable  ::  CVRR(:)  ! convective rain rate
Real, Allocatable  ::  CVSR(:)  ! convective snow rate
Real, Allocatable  ::  CONVCC(:,:)  ! convective cloud cover
Real, Allocatable  ::  CF_LSC(:,:)  ! bulk cloud fraction after LSC
Real, Allocatable  ::  TINC_CVN(:,:)! Temp incr across conv
Real, Allocatable  ::  TINC_PPN(:,:)! Temp incr across LS precip
Real, Allocatable  ::  QCL_LSC(:,:) ! cloud liq water after LSC
Real, Allocatable  ::  QCL_ADV(:,:) ! cloud liq water after advection

End Module ac_diagnostics_Mod
