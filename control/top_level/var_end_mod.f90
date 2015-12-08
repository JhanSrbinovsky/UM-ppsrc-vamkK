! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  var_end_mod

MODULE var_end_mod

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

! Description:  End values for variable resolution grid
!
! Method:   These are calculated in set_var_grid via setcona. 
!           
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

REAL  :: lambda_p_end   = rmdi
REAL  :: lambda_u_end   = rmdi
REAL  :: phi_p_end      = rmdi
REAL  :: phi_v_end      = rmdi
REAL  :: dlambda_p_end  = rmdi 
REAL  :: dlambda_u_end  = rmdi
REAL  :: dphi_p_end     = rmdi
REAL  :: dphi_v_end     = rmdi

END MODULE var_end_mod
