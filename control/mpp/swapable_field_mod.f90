! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP
MODULE swapable_field_mod

IMPLICIT NONE
  TYPE swapable_field_pointer_type
    INTEGER :: field_type
    INTEGER :: levels
    INTEGER :: rows
    LOGICAL :: vector
    REAL, POINTER :: field(:, :, :)
    REAL, POINTER :: field_2d(:, :)
  END TYPE swapable_field_pointer_type
END MODULE swapable_field_mod
