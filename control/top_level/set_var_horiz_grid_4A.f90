MODULE set_var_horiz_grid_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
IMPLICIT NONE
CONTAINS
SUBROUTINE set_var_horiz_grid(l_cartesian,                                  &
                              a_len2_coldepc,a_len2_rowdepc,                &   
                              Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in )

USE conversions_mod

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY : ereport
USE PrintStatus_mod

USE integrity_mod
USE horiz_grid_mod

USE atm_fields_bounds_mod
USE proc_info_mod, ONLY: datastart=>l_datastart,model_domain,         &
                         global_row_length,global_rows, me
USE um_parvars, ONLY : halo_i, halo_j
USE domain_params

IMPLICIT NONE
!
! Description: Computes the variable resolution horizontal ENDGame grid.
!  
!
! Method: Two separate cases have been coded: (i) even number of points and
!         (ii) odd number points across the domain. In (ii) the grid is
!         construction ensures symmetry.
!
! Documentation: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
!
LOGICAL, INTENT(IN) :: l_cartesian

INTEGER, INTENT(IN) :: a_len2_coldepc,a_len2_rowdepc

REAL, INTENT(IN) ::                                                &
    Lambda_p_in(global_row_length),                                &
    Lambda_u_in(MIN(global_row_length,                             &
                         global_row_length *a_len2_coldepc +1)),   &
    Phi_p_in(global_rows),                                         &
    Phi_v_in(MIN(global_rows+1, global_rows *a_len2_rowdepc +1))


! Local variables
LOGICAL :: l_ideal_stretch
INTEGER :: i, j, i_mid, j_mid, gi, gj

INTEGER :: Nxi1H, Nxi2H, is

REAL :: glob_xi1(-2*halo_i:2*(global_row_length+halo_i)), &
        glob_xi2(-2*halo_j:2*(global_rows+1+halo_j))

REAL :: Dx, xstretch, ystretch, domain_length, domain_width

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER alloc_stat,ErrorStatus


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('SET_VAR_HORIZ_GRID',zhook_in,zhook_handle)

!
! If any of the delta's are zero then the idealised namelist
! has not been set for the variable resolution grid
!
IF( delta_xi1_L*delta_xi1_H*delta_xi2_L*delta_xi2_H == 0.0) THEN
   l_ideal_stretch = .FALSE.
ELSE
   l_ideal_stretch = .TRUE.
END IF

! domain set up:
!
!  <-    X1   -><-    Y    -><-    Z     -><-    Y    -><-    X2   ->
!  <- DxLow   -><-  DxVar  -><-   DxHigh -><-  DxVar  -><-   DxLow ->
IF ( .NOT. l_ideal_stretch ) THEN
  IF ( PrintStatus >= PrStatus_Normal ) THEN
    WRITE(6,fmt='(A)') 'Reading variable resolution grid from start dump'
  END IF
ELSE
  IF ( PrintStatus >= PrStatus_Normal ) THEN
    WRITE(6,fmt='(A)') 'Setting up varaible resolution grid'
    WRITE(6,fmt='(A,2E16.8)') 'Low/High resolution delta lambda: ', &
                               delta_xi1_L,delta_xi1_H
    WRITE(6,fmt='(A,2E16.8)') 'Low/High resolution delta phi: ',    &
                               delta_xi2_L,delta_xi2_H
  END IF
  IF ( MIN(delta_xi1_H,delta_xi1_L) < 1e-3 ) THEN
     ErrorStatus = MIN(delta_xi1_H,delta_xi1_L)
     CALL Ereport("set_var_horiz_grid_4A", ErrorStatus,                 &
                  "Grid spacing in lambda direction is too small/negative")  
  END IF
  IF ( MIN(delta_xi2_H,delta_xi2_L) < 1e-3 ) THEN
     ErrorStatus = MIN(delta_xi2_H,delta_xi2_L)
     CALL Ereport("set_var_horiz_grid_4A", ErrorStatus,                 &
                  "Grid spacing in phi direction is too small/negative")  
  END IF
! Compute number of high resolution cells
  Nxi1H = (global_row_length-1) - 2*Nxi1L - 2*Nxi1V
  Nxi2H = (global_rows-1)       - 2*Nxi2L - 2*Nxi2V

  IF ( Nxi1H < 0 ) THEN
     ErrorStatus = Nxi1H
     CALL Ereport("set_var_horiz_grid_4A", ErrorStatus,                 &
                "Negative number of high res cells in lambda")
  END IF 
  IF ( Nxi2H < 0 ) THEN
     ErrorStatus = Nxi2H
     CALL Ereport("set_var_horiz_grid_4A", ErrorStatus,                 &
                  "Negative number of high res cells in phi")
  END IF
END IF


IF ( .NOT. l_ideal_stretch ) THEN
! Read in global stretched grid
  DO i=1,global_row_length-1
    glob_xi1_u(i)   = lambda_u_in(i+1)
    glob_xi1_p(i)   = lambda_p_in(i)
  END DO
  delta_xi1 = glob_xi1_p(2) - glob_xi1_p(1)

  DO i = -halo_i,0
    glob_xi1_u(i)   = lambda_u_in(1) + i*delta_xi1
  END DO
  DO i = 1-halo_i,0
    glob_xi1_p(i)   = lambda_p_in(1) + (i-1)*delta_xi1
  END DO
  
  DO i = global_row_length,global_row_length+halo_i
    glob_xi1_u(i) = lambda_u_in(global_row_length) + &
                    (i-global_row_length+1)*delta_xi1
  END DO
  DO i = global_row_length,global_row_length+halo_i
    glob_xi1_p(i) = lambda_p_in(global_row_length) + &
                    (i-global_row_length)*delta_xi1
  END DO

  DO i=1,global_rows
    glob_xi2_v(i) = phi_v_in(i+1)
    glob_xi2_p(i) = phi_p_in(i)
  END DO
  delta_xi2 = glob_xi2_p(2) - glob_xi2_p(1)
  
  DO i = -halo_j,0
    glob_xi2_v(i) = phi_v_in(1) + i*delta_xi2
  END DO
  DO i = 1-halo_j,0
    glob_xi2_p(i) = phi_p_in(1) + (i-1)*delta_xi2
  END DO
  
  DO i = global_rows+1,global_rows+halo_j+1
    glob_xi2_v(i) = phi_v_in(global_rows+1) + (i-global_rows)*delta_xi2
  END DO
  DO i = global_rows+1,global_rows+halo_j
    glob_xi2_p(i)   = phi_p_in(global_rows) + (i-global_rows)*delta_xi2
  END DO
ELSE
! Set up stretched grid based on namelist

! Compute geometric stretching factors
  xstretch = (delta_xi1_H/delta_xi1_L)**(1.0/(REAL(2.0*Nxi1V)+1.0))
  ystretch = (delta_xi2_H/delta_xi2_L)**(1.0/(REAL(2.0*Nxi2V)+1.0))

! Adjust grid if necessary
  IF( model_domain == mt_global .OR.                                    &
      model_domain == mt_cyclic_lam ) THEN

     base_xi1  = 0.0
     delta_xi1 = delta_xi1_H

     IF( model_domain == mt_global ) THEN
        delta_xi2 = delta_xi2_H
        base_xi2 = 0.0
     END IF
  END IF

! xi 1
  glob_xi1(-2*halo_i) = base_xi1 - halo_i*delta_xi1_L 
  DO i=1-2*halo_i,2*Nxi1L
     glob_xi1(i) = glob_xi1(i-1) +  0.5*delta_xi1_L
  END DO
  dx = 0.5*delta_xi1_L
  is = 2*Nxi1L+1
  DO i = is,is + 2*Nxi1V
    dx = dx*xstretch
    glob_xi1(i) = glob_xi1(i-1) + dx
  END DO
  is = is + 2*Nxi1V + 1
  DO i = is,is + 2*Nxi1H
    glob_xi1(i) = glob_xi1(i-1) + 0.5*delta_xi1_H
  END DO
  dx = 0.5*delta_xi1_H
  is = is + 2*Nxi1H + 1
  DO i = is,is + 2*Nxi1V
    dx = dx/xstretch
    glob_xi1(i) = glob_xi1(i-1) + dx
  END DO
  is = is + 2*Nxi1V + 1
  DO i = is, 2*(global_row_length+halo_i)
    glob_xi1(i) = glob_xi1(i-1) +  0.5*delta_xi1_L
  END DO

  glob_xi1_u(-halo_i) = glob_xi1(-2*halo_i)
  j = -2*halo_i + 1
  DO i = 1-halo_i,global_row_length+halo_i
     glob_xi1_p(i) = glob_xi1(j)
     glob_xi1_u(i) = glob_xi1(j+1)
     j = j + 2
  END DO

! xi 2
  glob_xi2(-2*halo_j) = base_xi2 - halo_j*delta_xi2_L 
  DO i=1-2*halo_j,2*Nxi2L
     glob_xi2(i) = glob_xi2(i-1) +  0.5*delta_xi2_L
  END DO
  dx = 0.5*delta_xi2_L
  is = 2*Nxi2L+1
  DO i = is,is + 2*Nxi2V
    dx = dx*ystretch
    glob_xi2(i) = glob_xi2(i-1) + dx
  END DO
  is = is + 2*Nxi2V + 1
  DO i = is,is + 2*Nxi2H
    glob_xi2(i) = glob_xi2(i-1) + 0.5*delta_xi2_H
  END DO
  dx = 0.5*delta_xi2_H
  is = is + 2*Nxi2H + 1
  DO i = is,is + 2*Nxi2V
    dx = dx/ystretch
    glob_xi2(i) = glob_xi2(i-1) + dx
  END DO
  is = is + 2*Nxi2V + 1
  DO i = is, 2*(global_rows+halo_j)
    glob_xi2(i) = glob_xi2(i-1) +  0.5*delta_xi2_L
  END DO

  glob_xi2_v(-halo_j) = glob_xi2(-2*halo_j)
  j = -2*halo_j + 1
  DO i = 1-halo_j,global_rows+halo_j
     glob_xi2_p(i) = glob_xi2(j)
     glob_xi2_v(i) = glob_xi2(j+1)
     j = j + 2
  END DO
  
END IF

domain_length = glob_xi1_u(global_row_length) - glob_xi1_u(0)
domain_width  = glob_xi2_v(global_rows)       - glob_xi2_v(0)


!----------------------------------------------------------------------

! If not a cartesian grid convert to Radians and reset domain to [2*pi,pi]
IF( .NOT. l_cartesian ) THEN

   xi1_u = pi_over_180*xi1_u
   xi1_p = pi_over_180*xi1_p

   xi2_v = pi_over_180*xi2_v
   xi2_p = pi_over_180*xi2_p

   glob_xi1_p = pi_over_180*glob_xi1_p
   glob_xi1_u = pi_over_180*glob_xi1_u
   glob_xi2_p = pi_over_180*glob_xi2_p
   glob_xi2_v = pi_over_180*glob_xi2_v

! Reset grid descriptors (degrees -> radians)

   base_xi1  = pi_over_180*base_xi1
   base_xi2  = pi_over_180*base_xi2
   delta_xi1 = pi_over_180*delta_xi1
   delta_xi2 = pi_over_180*delta_xi2

END IF

IF ( PrintStatus >= PrStatus_Diag ) THEN
  WRITE(6,'(A)') 'Global xi1 arrays'
  WRITE(6,fmt='(I4,2E16.8)') -halo_i,-9999999.99,glob_xi1_u(-halo_i) 
  DO i=1-halo_i, global_row_length+halo_i
    WRITE(6,fmt='(I4,2E16.8)') i,glob_xi1_p(i),glob_xi1_u(i) 
  END DO
  WRITE(6,'(A)') ' '
  WRITE(6,'(A)') 'Global xi2 arrays'
    WRITE(6,fmt='(I4,2E16.8)') -halo_j,-99999999.99,glob_xi2_v(-halo_j) 
  DO j=1-halo_j, global_rows+halo_j
    WRITE(6,fmt='(I4,2E16.8)') j,glob_xi2_p(j),glob_xi2_v(j) 
  END DO
END IF

! Copy to local processor grid
DO i=pdims_l%i_start,pdims_l%i_end
  gi = datastart(1) + i - 1
  xi1_p(i) = glob_xi1_p(gi) 
END DO

DO i=udims_l%i_start,udims_l%i_end
  gi = datastart(1) + i - 1
  xi1_u(i) = glob_xi1_u(gi) 
END DO

DO j=pdims_l%j_start,pdims_l%j_end
  gj = datastart(2) + j - 1
  xi2_p(j) = glob_xi2_p(gj)
END DO

DO j=vdims_l%j_start,vdims_l%j_end
  gj = datastart(2) + j - 1
  xi2_v(j) = glob_xi2_v(gj)
END DO

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                    xi1_p,           SIZE(xi1_p),           'xi1_p',  &
                    xi1_u,           SIZE(xi1_u),           'xi1_u',  &
                    xi2_p,           SIZE(xi2_p),           'xi2_p',  &
                    xi2_v,           SIZE(xi2_v),           'xi2_v'   &
                    )

IF (lhook) CALL dr_hook('SET_VAR_HORIZ_GRID',zhook_out,zhook_handle)

END SUBROUTINE set_var_horiz_grid
END MODULE set_var_horiz_grid_mod
