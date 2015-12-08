! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE pr_block4_mod

      IMPLICIT NONE

! Description:      Parameters used to define corners of
!                   diagnostic print domain
!
! Method:
!               The required parameters are initialised at start
!               and are used by subroutines print_block4_real
!               and print_block4_int
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      INTEGER  :: dom_w
      INTEGER  :: dom_eo    ! no offset yet
      INTEGER  :: dom_s
      INTEGER  :: dom_no    ! no offset yet
      INTEGER  :: ix
      INTEGER  :: jy

      END MODULE pr_block4_mod
