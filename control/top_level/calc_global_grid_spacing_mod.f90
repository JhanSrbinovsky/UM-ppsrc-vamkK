! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose : Calculate useful grid data (spacing and recipricals
!
!   Language: FORTRAN 90
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No P0
!
!  -------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

MODULE calc_global_grid_spacing_mod
IMPLICIT NONE
CONTAINS
  SUBROUTINE  calc_global_grid_spacing

  USE horiz_grid_mod
  USE um_parvars,    ONLY : halo_i, halo_j
  USE proc_info_mod, ONLY: global_row_length,global_rows

  IMPLICIT NONE 
  INTEGER i

  DO i = 1-halo_i, global_row_length+halo_i-1
    glob_dxi1_p(i)  = glob_xi1_p(i+1) - glob_xi1_p(i)
    glob_rdxi1_p(i) = 1.0/glob_dxi1_p(i)
  END DO

  DO i = -halo_i, global_row_length+halo_i-1
    glob_dxi1_u(i)  = glob_xi1_u(i+1) - glob_xi1_u(i)
    glob_rdxi1_u(i) = 1.0/glob_dxi1_u(i)
  END DO

  DO i = 1-halo_j, global_rows+halo_j-1
    glob_dxi2_p(i)  = glob_xi2_p(i+1) - glob_xi2_p(i)
    glob_rdxi2_p(i) = 1.0/glob_dxi2_p(i)
  END DO

  DO i = -halo_j, global_rows+halo_j-1
    glob_dxi2_v(i)  = glob_xi2_v(i+1) - glob_xi2_v(i)
    glob_rdxi2_v(i) = 1.0/glob_dxi2_v(i)
  END DO

  END SUBROUTINE calc_global_grid_spacing
END MODULE calc_global_grid_spacing_mod
