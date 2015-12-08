! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Constants necessary for radiation mask and trop levels in radiation scheme.

MODULE rad_mask_trop_mod

! Description:
! This module is used to hold values set in atmos_init.

! Method:
! The arrays are calculated in routine control/top_level/
! set_radmask called from atmos_init,
! and used in short_ and long_wave_radiation routines.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Radiation masking arrays and tropopause levels

REAL,    ALLOCATABLE  ::  es_space_interp(:,:,:)
LOGICAL, ALLOCATABLE  ::  rad_mask (:,:)
INTEGER  ::  min_trop_level  ! first reference model level about 700h
INTEGER  ::  max_trop_level  ! first reference model level about 50hp

END MODULE rad_mask_trop_mod
