! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  var_cubic_mod

      MODULE var_cubic_mod
      IMPLICIT NONE

! Description:  Pre-calculated variable resolution coefficients
!               for cubic Lagrange interpolation
!
! Method:   These are calculated in set_coeff_lagrange which is
!           called from set_var_grid via setcona. 
!           The array values are set in
!           init_cubic
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      REAL,  ALLOCATABLE  :: lambda_p_rm(:)
      REAL,  ALLOCATABLE  :: lambda_p_rp(:)
      REAL,  ALLOCATABLE  :: lambda_u_rm(:)
      REAL,  ALLOCATABLE  :: lambda_u_rp(:)
      REAL,  ALLOCATABLE  :: phi_p_rm(:,:)
      REAL,  ALLOCATABLE  :: phi_p_rp(:,:)
      REAL,  ALLOCATABLE  :: phi_v_rm(:,:)
      REAL,  ALLOCATABLE  :: phi_v_rp(:,:)
      REAL,  ALLOCATABLE  :: recip_lambda_p_m(:)
      REAL,  ALLOCATABLE  :: recip_lambda_p(:)
      REAL,  ALLOCATABLE  :: recip_lambda_p_p(:)
      REAL,  ALLOCATABLE  :: recip_lambda_p_p2(:)
      REAL,  ALLOCATABLE  :: recip_lambda_u_m(:)
      REAL,  ALLOCATABLE  :: recip_lambda_u(:)
      REAL,  ALLOCATABLE  :: recip_lambda_u_p(:)
      REAL,  ALLOCATABLE  :: recip_lambda_u_p2(:)
      REAL,  ALLOCATABLE  :: recip_phi_p_m(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_p(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_p_p(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_p_p2(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_v_m(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_v(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_v_p(:,:)
      REAL,  ALLOCATABLE  :: recip_phi_v_p2(:,:)
 
      END MODULE var_cubic_mod

