MODULE eg_parameters_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame departure points
!  
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

IMPLICIT NONE

REAL :: pole_consts(4)

!
! Flag for making linear solver x-y symmetric
LOGICAL :: L_symm_solve  
!
! Preconditioner variables ! Not used yet!
!
REAL    :: SSOR_rlx
INTEGER :: Nterms, Pre_its

INTEGER :: n_filt_ND = 8  ! number of initial filter sweeps when 
                          ! a ND dump is used
                          
LOGICAL :: interp_dpt_pt = .FALSE.
LOGICAL :: interp_log_rho = .FALSE.
LOGICAL :: total_conv_outer = .FALSE.
LOGICAL :: total_conv_inner = .FALSE.

LOGICAL :: L_slice = .FALSE.     ! SLICE conservative semi-Lagrangian scheme

LOGICAL :: l_rho_av_zz = .TRUE.
LOGICAL :: l_cartesian

END MODULE eg_parameters_mod
