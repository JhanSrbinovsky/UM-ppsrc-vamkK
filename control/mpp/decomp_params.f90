! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Magic numbers  defining decompositions

MODULE decomp_params

! Description:
!    This data module contains magic numbers defining decompositions
!    for MPP components
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP
  
  IMPLICIT NONE


  INTEGER, PARAMETER :: max_decomps           = 3

! not set decomp
  INTEGER, PARAMETER :: decomp_unset          = 0

! dead decomps that are still referenced, but not actually used
  INTEGER, PARAMETER :: decomp_standard_ocean = -1
  INTEGER, PARAMETER :: decomp_nowrap_ocean   = -2

! decomps for rcf
  INTEGER, PARAMETER :: decomp_rcf_input      = 1
  INTEGER, PARAMETER :: decomp_rcf_output     = 2

! decomps for atmos
  INTEGER, PARAMETER :: decomp_standard_atmos = 1
  INTEGER, PARAMETER :: decomp_standard_wave  = 2

!decomps for small execs
  INTEGER, PARAMETER :: decomp_smexe          = 3

END MODULE decomp_params
