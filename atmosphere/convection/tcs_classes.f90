! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing a collection of tcs classes and overloading 
! allocation and deallocation subroutines
!
MODULE tcs_classes

  USE tcs_class_interface
  USE tcs_class_similarity
  USE tcs_class_scales
  USE tcs_class_cloud

  IMPLICIT NONE
  !
  ! Description:
  ! This module defines the tcs warm rain "cloud_input" derived type and 
  ! subroutines for allocating and deallocating an instance of this class.
  !
  !
  ! Method:
  !  Creates new interfaces to overloads tcs_allocate and tcs_deallocate
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.

  INTERFACE tcs_allocate
     MODULE PROCEDURE allocate_scales
     MODULE PROCEDURE allocate_similarity
     MODULE PROCEDURE allocate_interface_input
     MODULE PROCEDURE allocate_cloud_input
  END INTERFACE

  INTERFACE tcs_deallocate
     MODULE PROCEDURE deallocate_scales
     MODULE PROCEDURE deallocate_similarity
     MODULE PROCEDURE deallocate_interface_input
     MODULE PROCEDURE deallocate_cloud_input
  END INTERFACE

END MODULE tcs_classes
