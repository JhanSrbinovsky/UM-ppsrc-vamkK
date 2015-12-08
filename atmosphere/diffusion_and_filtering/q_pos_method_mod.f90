! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Module q_pos_method_mod - parameters for Q_POS methods
!----------------------------------------------------------------------
MODULE q_pos_method_mod

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90

IMPLICIT NONE

INTEGER, PARAMETER :: q_pos_original = 1      ! the old "method 1"
INTEGER, PARAMETER :: q_pos_local    = 2      ! the old "method 2" with local 
                                              ! scavenging
INTEGER, PARAMETER :: q_pos_reset    = 3      ! Reset to minimum
INTEGER, PARAMETER :: q_pos_column   = 4      ! conserve within column
INTEGER, PARAMETER :: q_pos_level    = 5      ! conserve within level
INTEGER, PARAMETER :: q_pos_hybrid   = 6      ! hygrid of column and level

END MODULE q_pos_method_mod
