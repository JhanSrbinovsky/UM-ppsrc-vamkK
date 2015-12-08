! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module holding logical that selects underlying grid.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dynamics_grid_mod

IMPLICIT NONE

LOGICAL :: L_vatpoles = .FALSE.     ! currently only set true when 
                                    ! l_endgame is true
                                    ! not of concern to standalone jules
                                    ! but required for compilation

END MODULE dynamics_grid_mod
