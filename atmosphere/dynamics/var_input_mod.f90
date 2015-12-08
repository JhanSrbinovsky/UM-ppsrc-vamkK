! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input settings as used for variable resolution
!   
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

! Method:
!   Switches are initialised to false and read in from the
!   UMUI. The module may then be used directly where the switches
!   are needed within the dynamics semi_lagrangian code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE var_input_mod

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE

!  variable resolution control
LOGICAL :: L_regular = .FALSE.  ! true if NOT variable resolution

INTEGER :: lam_var = imdi   ! number of variable res. lambda intervals
INTEGER :: phi_var = imdi   ! number of variable res. phi intervals

REAL :: var_ratio  = rmdi   ! grid-stretcing ratio for var grid
REAL :: lam_ratio  = rmdi   ! scaling of original grid to high res grid
REAL :: phi_ratio  = rmdi   ! scaling of original grid to high res grid
REAL :: lam_frac   = rmdi   ! proportion of reg. points in West of domain
REAL :: phi_frac   = rmdi   ! proportion of reg. points in South of domain 

END MODULE  var_input_mod
