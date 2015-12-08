! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing
MODULE c_grid2grid_mod

IMPLICIT NONE
! C_GRID2GRID start
! Description: Parameters for LAM river routing
! Author: Vicky Bell, CEH Wallingford

      REAL, PARAMETER :: cland  = 0.2  ! land wave speed (m/s)
      REAL, PARAMETER :: criver = 0.2  ! subsurf river wave speed (m/s)
      REAL, PARAMETER :: cbland  = 0.18! subsurf land wave speed (m/s)
      REAL, PARAMETER :: cbriver = 0.18! subsurf river wave speed (m/s)
      REAL, PARAMETER :: runoff_factor = 0.7
!                                      ! runoff volume factor
      REAL, PARAMETER :: retl = 0.     ! return flow (land squares) (<1)
      REAL, PARAMETER :: retr = 0.15   ! return flow (river squares) (<1)
      REAL, PARAMETER :: slfac = 0.    ! slope factor (not used yet)
      INTEGER, PARAMETER :: a_thresh = 1 ! threshold area
! END C_GRID2GRID
!---------------------------------------------------------------------

END MODULE c_grid2grid_mod
