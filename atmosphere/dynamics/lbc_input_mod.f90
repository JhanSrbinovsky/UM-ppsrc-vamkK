! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input settings as used for lbcs and trapping
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

MODULE lbc_input_mod

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

LOGICAL :: L_LBC_balance = .TRUE.   ! T: impose vertically balanced Exners
                                    !    and rho in LBCs (set w=0)
                                    ! F: leave Exners alone

! Switch for new (correct) lbcs developed in 2008
LOGICAL :: L_lbc_new     = .TRUE.   ! .true. for new lbc switch
LOGICAL :: L_lbc_old     = .FALSE.  ! false for new lbc treatment
                           ! set in setcon
                           
LOGICAL :: L_transparent = .FALSE.  ! default is .false.
                                    ! .true. for transparent lbcs 
                           ! Does not apply lbcs to horizontal winds 
                           ! in the solver
LOGICAL :: L_blend         ! default is .true. for blending lbcs
                           ! for moisture, tracers and w fields

INTEGER :: n_rims_to_do  = imdi     ! rim size for LAM

END MODULE  lbc_input_mod
