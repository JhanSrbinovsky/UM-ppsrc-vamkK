! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares allocatable arrays for turbulent diffusion 
!
!  Code Owner: See Unified Model Code Owner's HTML page
!  This file belongs in section: Top Level

MODULE turb_diff_ctl_mod

! Subgrid turbulence scheme pointers.

IMPLICIT NONE

INTEGER  :: turb_ustart    !  start u/th row on this processor
INTEGER  :: turb_uend      !  end u/th row on this processor
INTEGER  :: turb_vstart    !  start v row on this processor
INTEGER  :: turb_vend      !  end v row on this processor
INTEGER  :: turb_phi_st    !  start row for phi diff on this processor
INTEGER  :: turb_phi_end   !  end row for phi diff on this processor

! Subgrid turbulence scheme arrays.


REAL, ALLOCATABLE, TARGET :: visc_m(:,:,:)  ! diffusion coefficient for
                                            !  momentum
REAL, ALLOCATABLE, TARGET :: visc_h(:,:,:)  ! diffusion coefficient for 
                                            !  heat and moisture
                       ! level 1 value is dummy for use in diagnostics
REAL, ALLOCATABLE ::  rneutml_sq(:,:,:)     ! mixing length scale (m)
REAL, ALLOCATABLE ::  shear(:,:,:)
REAL, ALLOCATABLE ::  delta_smag(:,:)       ! grid size for Smag (m)
REAL, ALLOCATABLE ::  max_diff(:,:)         ! max diffusion coeff for run  

END MODULE turb_diff_ctl_mod
