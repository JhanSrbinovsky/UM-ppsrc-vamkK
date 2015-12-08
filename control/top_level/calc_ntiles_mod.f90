! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE calc_ntiles_mod

IMPLICIT NONE
!
! Description:
!  Jules only: Calculates the value of ntiles based on the input variables
!
! Method:
!
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standard UMDP3 vn8.2.
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Top Level

CONTAINS

SUBROUTINE calc_ntiles(l_aggregate,npft,nnvg,ntiles)

IMPLICIT NONE

LOGICAL, INTENT(IN) ::                                                 &
       l_aggregate

INTEGER, INTENT(IN) ::                                                 &
       npft,                                                           &
       nnvg

INTEGER, INTENT(OUT)::                                                 &
       ntiles

! Calculate ntiles based on whether we-re using an aggregate surface scheme
IF(l_aggregate) THEN
  ntiles=1
ELSE
  ntiles = npft + nnvg
END IF

END SUBROUTINE calc_ntiles

END MODULE calc_ntiles_mod
