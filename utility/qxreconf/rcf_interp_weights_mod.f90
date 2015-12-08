! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Defines interpolation weight variables

MODULE Rcf_Interp_Weights_Mod

! Description:
! This module defines the interpolation data-types, containing
! weights etc for interpolations.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!-------------------------------------------------------------------
! Horizontal interpolation information
!-------------------------------------------------------------------

! Some parameters used for array sizing for weights
INTEGER, PARAMETER        :: Idim = 4
INTEGER, PARAMETER        :: Icof = 3

! These are the horizontal interpolation methods that could be
! used by reconfiguration.

INTEGER, PARAMETER        :: bilinear          = 1
INTEGER, PARAMETER        :: area_weighted     = 2
INTEGER, PARAMETER        :: nearest_neighbour = 3

!------------------
! Read in via namelist
INTEGER, SAVE            :: h_int_method = imdi  ! Interpolation method
LOGICAL, SAVE            :: l_limit_rotations = .FALSE.
LOGICAL, SAVE            :: smcp_int_nearest_neighbour = .FALSE.
! -----------------

LOGICAL, SAVE            :: h_int_active   ! Interpolation on or off (p grid)
LOGICAL, SAVE            :: h_int_active_u ! Interpolation on or off (u grid)
LOGICAL, SAVE            :: h_int_active_v ! Interpolation on or off (v grid)
LOGICAL, SAVE            :: l_same_rotation = .FALSE.

  ! Weights and indexes etc for bi-linear interp
INTEGER, SAVE, ALLOCATABLE, TARGET   :: bl_index_b_l( :, : )
INTEGER, SAVE, ALLOCATABLE, TARGET   :: bl_index_b_r( :, : )
REAL, SAVE, ALLOCATABLE, TARGET      :: weight_t_r( :, : )
REAL, SAVE, ALLOCATABLE, TARGET      :: weight_b_r( :, : )
REAL, SAVE, ALLOCATABLE, TARGET      :: weight_t_l( :, : )
REAL, SAVE, ALLOCATABLE, TARGET      :: weight_b_l( :, : )

  ! Weights and indexes etc for area-weighted interp
    ! Fixed size
REAL, SAVE                   :: aw_area_box( idim )

    ! Allocatables
INTEGER, SAVE, ALLOCATABLE, TARGET   :: aw_index_targ_lhs( :, : )
INTEGER, SAVE, ALLOCATABLE, TARGET   :: aw_index_targ_top( :, : )
REAL, SAVE,   ALLOCATABLE, TARGET    :: aw_colat_t( :, : )
REAL, SAVE,  ALLOCATABLE, TARGET     :: aw_long_l( :, : )

  ! Rotation coefficients
REAL, SAVE, ALLOCATABLE      :: Coeff1( : )
REAL, SAVE, ALLOCATABLE      :: Coeff2( : )
REAL, SAVE, ALLOCATABLE      :: Coeff3( : )
REAL, SAVE, ALLOCATABLE      :: Coeff4( : )

END MODULE Rcf_Interp_Weights_Mod

