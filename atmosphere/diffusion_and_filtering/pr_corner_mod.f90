! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

Module PR_CORNER_MOD

! Description:      Parameters used to define corners of 
!                   diagnostic print domain
!
! Method:
!               The required parameters are initialised at start 
!               and are used by subroutine  print_dom_corn
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
      IMPLICIT NONE

      Integer  :: dom_w
      Integer  :: dom_eo    ! no offset yet
      Integer  :: dom_s
      Integer  :: dom_no    ! no offset yet
      Integer  :: ix
      Integer  :: jy

End Module PR_CORNER_MOD
