! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Semi-Lagrangian trajectory information.

      MODULE sl_depart_mod
      IMPLICIT NONE

! Description:
!          Semi-Lagrangian trajectory information for time step          
!          that is shared by more than 1 subtoutine.          
!
! Method:
!         The required parameters are set
!         in the routine which uses them first and are deallocated
!         when no longer needed in the time step.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


      REAL,  ALLOCATABLE :: depart_lambda_w(:,:,:)
      REAL,  ALLOCATABLE :: depart_phi_w(:,:,:)
      REAL,  ALLOCATABLE :: depart_r_w(:,:,:)

      END MODULE sl_depart_mod
