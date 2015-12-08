! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE set_horiz_interp_consts_mod
IMPLICIT NONE

CONTAINS
SUBROUTINE set_horiz_interp_consts()

USE conversions_mod

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY : ereport
USE PrintStatus_mod

USE integrity_mod
USE horiz_grid_mod
USE interp_grid_const_mod

USE atm_fields_bounds_mod
USE proc_info_mod, ONLY: datastart=>l_datastart,model_domain,         &
                         global_row_length,global_rows, me
USE um_parvars, ONLY : halo_i, halo_j
USE domain_params

IMPLICIT NONE
!
! Description: Derives constants from the variable horizontal ENDGame grid
!              for use in the cubic Lagrange interpolation.
!
! Documentation: ENDGame formulation version 1.01
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
!

! Local variables
INTEGER :: i, j, i_mid, j_mid, gi, gj

INTEGER :: Nxi1H, Nxi2H, is

REAL :: s, t

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('SET_HORIZ_INTERP_CONST',zhook_in,zhook_handle)


ALLOCATE( sig_xi1_p(pdims_l%i_start+1:pdims_l%i_end-2),                   &
          sig_xi1_u(udims_l%i_start+1:udims_l%i_end-2),                   &
          sig_xi2_p(pdims_l%j_start+1:pdims_l%j_end-2),                   &
          sig_xi2_v(vdims_l%j_start+1:vdims_l%j_end-2) )

ALLOCATE( tau_xi1_p(pdims_l%i_start+1:pdims_l%i_end-2),                   &
          tau_xi1_u(udims_l%i_start+1:udims_l%i_end-2),                   &
          tau_xi2_p(pdims_l%j_start+1:pdims_l%j_end-2),                   &
          tau_xi2_v(vdims_l%j_start+1:vdims_l%j_end-2) )

ALLOCATE( q1_p(-1:2,pdims_l%i_start+1:pdims_l%i_end-2),                   &
          q1_u(-1:2,udims_l%i_start+1:udims_l%i_end-2),                   &
          q2_p(-1:2,pdims_l%j_start+1:pdims_l%j_end-2),                   &
          q2_v(-1:2,vdims_l%j_start+1:vdims_l%j_end-2) )


! Calculate ratios and copy to local position

! E-W

DO i = pdims_l%i_start+1, pdims_l%i_end-2
   gi = datastart(1) + i - 1

   sig_xi1_p(i) =-(glob_xi1_p(gi) - glob_xi1_p(gi-1))                     &
                 /(glob_xi1_p(gi+1) - glob_xi1_p(gi))

   tau_xi1_p(i) = (glob_xi1_p(gi+2) - glob_xi1_p(gi+1))                   &
                 /(glob_xi1_p(gi+1) - glob_xi1_p(gi))   + 1.0
END DO

DO i = udims_l%i_start+1, udims_l%i_end-2
   gi = datastart(1) + i - 1

   sig_xi1_u(i) =-(glob_xi1_u(gi) - glob_xi1_u(gi-1))                     &
                 /(glob_xi1_u(gi+1) - glob_xi1_u(gi))

   tau_xi1_u(i) = (glob_xi1_u(gi+2) - glob_xi1_u(gi+1))                   &
                 /(glob_xi1_u(gi+1) - glob_xi1_u(gi))   + 1.0

END DO

!N-S

DO j = pdims_l%j_start+1, pdims_l%j_end-2
   gj = datastart(2) + j - 1

   sig_xi2_p(j) =-(glob_xi2_p(gj) - glob_xi2_p(gj-1))                     &
                 /(glob_xi2_p(gj+1) - glob_xi2_p(gj))

   tau_xi2_p(j) = (glob_xi2_p(gj+2) - glob_xi2_p(gj+1))                   &
                 /(glob_xi2_p(gj+1) - glob_xi2_p(gj))   + 1.0

END DO

DO j = vdims_l%j_start+1, vdims_l%j_end-2
   gj = datastart(2) + j - 1

   sig_xi2_v(j) =-(glob_xi2_v(gj) - glob_xi2_v(gj-1))                     &
                 /(glob_xi2_v(gj+1) - glob_xi2_v(gj))

   tau_xi2_v(j) = (glob_xi2_v(gj+2) - glob_xi2_v(gj+1))                   &
                 /(glob_xi2_v(gj+1) - glob_xi2_v(gj))   + 1.0

END DO

! Now calculate normalization constants

DO i = pdims_l%i_start+1, pdims_l%i_end-2
   s          = sig_xi1_p(i)
   t          = tau_xi1_p(i)
   q1_p(-1,i) = 1.0/(s*(s-1.0)*(s-t))
   q1_p(0,i)  =-1.0/(s*t)
   q1_p(1,i)  = 1.0/((1.0-s)*(1.0-t))
   q1_p(2,i)  = 1.0/((t-s)*t*(t-1.0))
END DO

DO i = udims_l%i_start+1, udims_l%i_end-2
   s          = sig_xi1_u(i)
   t          = tau_xi1_u(i)
   q1_u(-1,i) = 1.0/(s*(s-1.0)*(s-t))
   q1_u(0,i)  =-1.0/(s*t)
   q1_u(1,i)  = 1.0/((1.0-s)*(1.0-t))
   q1_u(2,i)  = 1.0/((t-s)*t*(t-1.0))
END DO

DO j = pdims_l%j_start+1, pdims_l%j_end-2
   s          = sig_xi2_p(j)
   t          = tau_xi2_p(j)
   q2_p(-1,j) = 1.0/(s*(s-1.0)*(s-t))
   q2_p(0,j)  =-1.0/(s*t)
   q2_p(1,j)  = 1.0/((1.0-s)*(1.0-t))
   q2_p(2,j)  = 1.0/((t-s)*t*(t-1.0))
END DO

DO j = vdims_l%j_start+1, vdims_l%j_end-2
   s          = sig_xi2_v(j)
   t          = tau_xi2_v(j)
   q2_v(-1,j) = 1.0/(s*(s-1.0)*(s-t))
   q2_v(0,j)  =-1.0/(s*t)
   q2_v(1,j)  = 1.0/((1.0-s)*(1.0-t))
   q2_v(2,j)  = 1.0/((t-s)*t*(t-1.0))
END DO

IF (lhook) CALL dr_hook('SET_HORIZ_INTERP_CONST',zhook_out,zhook_handle)

END SUBROUTINE set_horiz_interp_consts
END MODULE set_horiz_interp_consts_mod
