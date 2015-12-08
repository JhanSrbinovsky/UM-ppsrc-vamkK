! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with convection.

MODULE cv_cntl_mod


  USE domain_params
  IMPLICIT NONE
  SAVE

! Description:
!   Module containing logical switches from &CNTLATM which control the
!   convection code. 
!
! Method:
!   Switches are declared initialised with the corresponding values set
!   in CNTLATM. The mod cannot use the same names as in &CNTLATM because
!   the include file cntlatm.h occurs in too many routines to make it
!   practical.

!   This mod should only be used in atm_step.F90 (UM) and scm_main.F90 (SCM)
!   where the CNTLATM variables are copied to the module
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:

  Logical :: lcv_3d_ccw    = .FALSE. ! Requires l_ccrad=.TRUE.
                                     ! .TRUE. : Radiation to use 3d ccw
                                     !          profile passed to it from
                                     !          convection.
                                     ! .FALSE.: Radiation constructs
                                     !          mean CCW profile from cclwp,
                                     !          ccb and cct as in original.

END MODULE cv_cntl_mod
