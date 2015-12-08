! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_horz_drag_mod

! Description: height dependent horizontal drag
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

REAL, SAVE, ALLOCATABLE ::  cd_u(:,:,:) ! horizontal drag
                                        ! coefficient on u points
REAL, SAVE, ALLOCATABLE ::  cd_v(:,:,:) ! horizontal drag
                                        ! coefficient on v points
REAL, SAVE, ALLOCATABLE ::  r_u_store(:,:,:) ! store the source terms in case 
                                             ! of an implict/explicit mix
REAL, SAVE, ALLOCATABLE ::  r_v_store(:,:,:) 

LOGICAL, SAVE :: l_impl_horz_drag
LOGICAL, SAVE :: l_expl_horz_drag

END MODULE eg_horz_drag_mod
