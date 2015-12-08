! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  var_grid_mod

      MODULE var_grid_mod
      IMPLICIT NONE

! Description:  Variable resolution grid arrays

! Method:   These are calculated in set_var_grid via setcona. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


      REAL,  ALLOCATABLE  :: glambda_p(:)
      REAL,  ALLOCATABLE  :: lambda_p(:)
      REAL,  ALLOCATABLE  :: phi_p(:,:)
      REAL,  ALLOCATABLE  :: glambda_u(:)
      REAL,  ALLOCATABLE  :: lambda_u(:)
      REAL,  ALLOCATABLE  :: phi_v(:,:)
      REAL,  ALLOCATABLE  :: gdlambda_p(:)
      REAL,  ALLOCATABLE  :: dlambda_p(:)
      REAL,  ALLOCATABLE  :: dphi_p(:,:)
      REAL,  ALLOCATABLE  :: gdlambda_u(:)
      REAL,  ALLOCATABLE  :: dlambda_u(:)
      REAL,  ALLOCATABLE  :: dphi_v(:,:)
      REAL,  ALLOCATABLE  :: grecip_dlamp(:)
      REAL,  ALLOCATABLE  :: recip_dlamp(:)
      REAL,  ALLOCATABLE  :: recip_dphip(:,:)
      REAL,  ALLOCATABLE  :: grecip_dlamu(:)
      REAL,  ALLOCATABLE  :: recip_dlamu(:)
      REAL,  ALLOCATABLE  :: recip_dphiv(:,:)
      REAL,  ALLOCATABLE  :: wt_lambda_p(:)
      REAL,  ALLOCATABLE  :: wt_phi_p(:,:)
      REAL,  ALLOCATABLE  :: wt_lambda_u(:)
      REAL,  ALLOCATABLE  :: wt_phi_v(:,:)
 
      END MODULE var_grid_mod
