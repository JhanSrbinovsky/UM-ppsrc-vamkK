! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module containing variables which are common to the tcs warm routines
!
MODULE tcs_common_warm

  USE tcs_class_scales, only:             &
     scales_conv

  USE tcs_class_interface, only:          &
     interface_input

  IMPLICIT NONE
  !
  ! Description:
  !   This routine holds a common set of variables for the tcs
  ! warm routines
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !
  
  TYPE(scales_conv)     :: scales     !
  TYPE(interface_input) :: cb     ! Cloud base
  TYPE(interface_input) :: cb_p1  ! Level above cloud base
  TYPE(interface_input) :: cb_m1  ! Level below cloud base
  TYPE(interface_input) :: inv    ! Inversion base
  TYPE(interface_input) :: inv_p1 ! Level above inv
  TYPE(interface_input) :: inv_m1 ! Level below inv
  
END MODULE tcs_common_warm
