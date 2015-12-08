! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_vert_damp_mod

! Description: height dependent vertical damping
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

REAL, SAVE, ALLOCATABLE ::  mu_w(:,:,:) ! vertical damping
                                        ! coefficient

END MODULE eg_vert_damp_mod
