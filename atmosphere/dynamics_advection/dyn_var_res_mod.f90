! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for variable resolution advection scheme.

Module dyn_var_res_Mod

! Description:
! This module is used to hold variable resolution values set in atmos_init.
!
! Method:
! The arrays are calculated in routine control/top_level/
! set_var_res called from atmos_init, 
! and used in variable resolution dynamics_advection routines.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Arguments

! VarRes Array co-ordinates in radians
REAL, Allocatable :: glambda_p ( : )
REAL, Allocatable :: glambda_u ( : )
REAL, Allocatable :: lambda_p ( : )
REAL, Allocatable :: lambda_u ( : )
REAL, Allocatable :: phi_p ( : , : )
REAL, Allocatable :: phi_v( : , : )
REAL, Allocatable :: gdlambda_p ( : )
REAL, Allocatable :: gdlambda_u ( : )
REAL, Allocatable :: dphi_p ( : , : )
REAL, Allocatable :: dphi_v( : , : )
REAL, Allocatable :: grecip_dlamp(:)
REAL, Allocatable :: grecip_dlamu(:)
REAL, Allocatable :: recip_dphip( : , : )
REAL, Allocatable :: recip_dphiv( : , : )
REAL, Allocatable :: wt_lambda_p(:)
REAL, Allocatable :: wt_lambda_u(:)
REAL, Allocatable :: wt_phi_p( : , : )
REAL, Allocatable :: wt_phi_v( : , : )
REAL, Allocatable :: lambda_p_rm(:)
REAL, Allocatable :: lambda_p_rp(:)
REAL, Allocatable :: lambda_u_rm(:)
REAL, Allocatable :: lambda_u_rp(:)
REAL, Allocatable :: phi_p_rm( : , : )
REAL, Allocatable :: phi_p_rp( : , : )
REAL, Allocatable :: phi_v_rm( : , : )
REAL, Allocatable :: phi_v_rp( : , : )
REAL, Allocatable :: Recip_lambda_p_m ( : )
REAL, Allocatable :: Recip_lambda_p_0 ( : )
REAL, Allocatable :: Recip_lambda_p_p ( : )
REAL, Allocatable :: Recip_lambda_p_p2( : )
REAL, Allocatable :: Recip_lambda_u_m ( : )
REAL, Allocatable :: Recip_lambda_u_0 ( : )
REAL, Allocatable :: Recip_lambda_u_p ( : )
REAL, Allocatable :: Recip_lambda_u_p2( : )
REAL, Allocatable :: Recip_phi_p_m ( : , : )
REAL, Allocatable :: Recip_phi_p_0 ( : , : )
REAL, Allocatable :: Recip_phi_p_p ( : , : )
REAL, Allocatable :: Recip_phi_p_p2( : , : )
REAL, Allocatable :: Recip_phi_v_m ( : , : )
REAL, Allocatable :: Recip_phi_v_0 ( : , : )
REAL, Allocatable :: Recip_phi_v_p ( : , : )
REAL, Allocatable :: Recip_phi_v_p2( : , : )

End Module dyn_var_res_Mod
