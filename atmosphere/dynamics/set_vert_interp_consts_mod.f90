! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE set_vert_interp_consts_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE set_vert_interp_consts(eta_rho_levels,eta_theta_levels,model_levels)

USE interp_grid_const_mod
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

!
! Description: Derives constants from the ENDGame vertical grid
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER, INTENT(IN)    :: model_levels
REAL,    INTENT(IN)    :: eta_rho_levels(model_levels)
REAL,    INTENT(IN)    :: eta_theta_levels(0:model_levels)
REAL                   :: p(-2:3), tmp, eta_rho_levels_e(0:model_levels+1)
INTEGER                :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SET_HORIZ_INTERP_CONST',zhook_in,zhook_handle)

ALLOCATE(q3_p(-1:2,0:model_levels+1), q3_w(-1:2,0:model_levels) )
ALLOCATE(q5_p(-2:3,0:model_levels+1), q5_w(-2:3,0:model_levels) )

! ste to a stupid value for debugging

q3_p = -9999999.99999
q3_w = -9999999.99999
q5_p = -9999999.99999
q5_w = -9999999.99999

! Set cubic values

eta_rho_levels_e(1:model_levels) = eta_rho_levels(:)
eta_rho_levels_e(0) = 0.0
eta_rho_levels_e(model_levels+1) = 1.0

DO k = 1, model_levels-1
   p(-1)      = eta_rho_levels_e(k-1)
   p(0)       = eta_rho_levels_e(k)
   p(1)       = eta_rho_levels_e(k+1)
   p(2)       = eta_rho_levels_e(k+2)

   tmp        = ( p(-1) - p(0) )*( p(-1) - p(1) )*( p(-1) - p(2))
   q3_p(-1,k) = 1.0/tmp

   tmp        = ( p(0) - p(-1) )*( p(0) - p(1) )*( p(0) - p(2))
   q3_p(0,k)  = 1.0/tmp

   tmp        = ( p(1) - p(-1) )*( p(1) - p(0) )*( p(1) - p(2))
   q3_p(1,k)  = 1.0/tmp

   tmp        = ( p(2) - p(-1) )*( p(2) - p(0) )*( p(2) - p(1))
   q3_p(2,k)  = 1.0/tmp
END DO

DO k = 1, model_levels-2
   p(-1)      = eta_theta_levels(k-1)
   p(0)       = eta_theta_levels(k)
   p(1)       = eta_theta_levels(k+1)
   p(2)       = eta_theta_levels(k+2)

   tmp        = ( p(-1) - p(0) )*( p(-1) - p(1) )*( p(-1) - p(2))
   q3_w(-1,k) = 1.0/tmp

   tmp        = ( p(0) - p(-1) )*( p(0) - p(1) )*( p(0) - p(2))
   q3_w(0,k)  = 1.0/tmp

   tmp        = ( p(1) - p(-1) )*( p(1) - p(0) )*( p(1) - p(2))
   q3_w(1,k)  = 1.0/tmp

   tmp        = ( p(2) - p(-1) )*( p(2) - p(0) )*( p(2) - p(1))
   q3_w(2,k)  = 1.0/tmp
END DO


! Set quintic values

DO k = 2, model_levels-2
   p(-2)      = eta_rho_levels_e(k-2)
   p(-1)      = eta_rho_levels_e(k-1)
   p(0)       = eta_rho_levels_e(k)
   p(1)       = eta_rho_levels_e(k+1)
   p(2)       = eta_rho_levels_e(k+2)
   p(3)       = eta_rho_levels_e(k+3)

   tmp        = ( p(-2) - p(-1) )*( p(-2) - p(0) )*( p(-2) - p(1))             &
               *( p(-2) - p(2) )*( p(-2) - p(3) )
   q5_p(-2,k) = 1.0/tmp

   tmp        = ( p(-1) - p(-2) )*( p(-1) - p(0) )*( p(-1) - p(1))             &
               *( p(-1) - p(2) )*( p(-1) - p(3) )
   q5_p(-1,k) = 1.0/tmp

   tmp        = ( p(0) - p(-2) )*( p(0) - p(-1) )*( p(0) - p(1))               &
               *( p(0) - p(2) )*( p(0) - p(3) )
   q5_p(0,k) = 1.0/tmp

   tmp        = ( p(1) - p(-2) )*( p(1) - p(-1) )*( p(1) - p(0))               &
               *( p(1) - p(2) )*( p(1) - p(3) )
   q5_p(1,k) = 1.0/tmp

   tmp        = ( p(2) - p(-2) )*( p(2) - p(-1) )*( p(2) - p(0))               &
               *( p(2) - p(1) )*( p(2) - p(3) )
   q5_p(2,k) = 1.0/tmp

   tmp        = ( p(2) - p(-2) )*( p(3) - p(-1) )*( p(3) - p(0))               &
               *( p(3) - p(1) )*( p(3) - p(2) )
   q5_p(3,k) = 1.0/tmp

END DO

DO k = 2, model_levels-3
   p(-2)      = eta_theta_levels(k-2)
   p(-1)      = eta_theta_levels(k-1)
   p(0)       = eta_theta_levels(k)
   p(1)       = eta_theta_levels(k+1)
   p(2)       = eta_theta_levels(k+2)
   p(3)       = eta_theta_levels(k+3)

   tmp        = ( p(-2) - p(-1) )*( p(-2) - p(0) )*( p(-2) - p(1))             &
               *( p(-2) - p(2) )*( p(-2) - p(3) )
   q5_w(-2,k) = 1.0/tmp

   tmp        = ( p(-1) - p(-2) )*( p(-1) - p(0) )*( p(-1) - p(1))             &
               *( p(-1) - p(2) )*( p(-1) - p(3) )
   q5_w(-1,k) = 1.0/tmp

   tmp        = ( p(0) - p(-2) )*( p(0) - p(-1) )*( p(0) - p(1))               &
               *( p(0) - p(2) )*( p(0) - p(3) )
   q5_w(0,k) = 1.0/tmp

   tmp        = ( p(1) - p(-2) )*( p(1) - p(-1) )*( p(1) - p(0))               &
               *( p(1) - p(2) )*( p(1) - p(3) )
   q5_w(1,k) = 1.0/tmp

   tmp        = ( p(2) - p(-2) )*( p(2) - p(-1) )*( p(2) - p(0))               &
               *( p(2) - p(1) )*( p(2) - p(3) )
   q5_w(2,k) = 1.0/tmp

   tmp        = ( p(2) - p(-2) )*( p(3) - p(-1) )*( p(3) - p(0))               &
               *( p(3) - p(1) )*( p(3) - p(2) )
   q5_w(3,k) = 1.0/tmp

END DO

IF (lhook) CALL dr_hook('SET_HORIZ_INTERP_CONST',zhook_out,zhook_handle)

END SUBROUTINE set_vert_interp_consts
END MODULE set_vert_interp_consts_mod
