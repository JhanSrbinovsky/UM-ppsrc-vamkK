! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for East-West coefficients in diffusion scheme.

Module diff_coeff_Mod

! Description:
! This module is used to hold diffusion coefficients set in atmos_init.
!
! Method:
! The diff_coeff_u/v arrays are calculated in routine control/top_level/
! set_diff called from atmos_init, 
! and used in diffusion_and_filtering routines.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! East-West diffusion coefficient Arrays

Real, Allocatable  ::  diff_coeff_u (:,:) 
Real, Allocatable  ::  diff_coeff_v (:,:)

End Module diff_coeff_Mod
