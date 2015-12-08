! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for LAM rotation in dynamics diagnostics section.

Module rot_coeff_Mod

! Description:
! This module is used to hold rotation coefficient values set in atmos_init.
!
! Method:
! The rot_coeff1/2 arrays are calculated in routine control/top_level/
! set_rot_coeff called from atmos_init, 
! and used in dynamics_diagnostics routines.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Arguments

! Rotation Coefficient Arrays

Real, Allocatable  ::  rot_coeff1 (:,:) 
Real, Allocatable  ::  rot_coeff2 (:,:)
Integer            ::  LAM_max_cfl(2)  ! Max cfl values for LAM BCs

End Module rot_coeff_Mod
