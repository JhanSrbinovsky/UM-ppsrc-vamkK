! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc.
!
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.

MODULE atmos_max_sizes
  
  IMPLICIT NONE

! Maximum sector size for I/O







  INTEGER, PARAMETER :: row_length_max   = 2048 ! maximum row length
  INTEGER, PARAMETER :: rows_max         = 1600 ! max no of rows


  ! maximum permitted size of a halo
  INTEGER, PARAMETER :: Max_Halo_Size    = 10




  INTEGER, PARAMETER :: model_levels_max = 250 ! max no of total levels


  INTEGER, PARAMETER :: horiz_dim_max=max(row_length_max,rows_max)
  INTEGER, PARAMETER :: wet_levels_max   = 250 ! max no of wet levels
  INTEGER, PARAMETER :: Max2DFieldSize   = row_length_max*rows_max
  INTEGER, PARAMETER :: MaxHaloArea      = horiz_dim_max*Max_Halo_Size

END MODULE atmos_max_sizes
