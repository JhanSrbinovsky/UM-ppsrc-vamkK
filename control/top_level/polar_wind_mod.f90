! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ MODULE polar_wind_mod

MODULE polar_wind_mod
  IMPLICIT NONE

! Description:   Parameters needed to derive polar u-component
!                from v next to pole.                    
!
! Method:
!         The required pointers are set in the routine init_polar wind
!         which is called from control/top_level/setcona 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

  REAL, ALLOCATABLE ::  mag_vector_np (:)
  REAL, ALLOCATABLE ::  dir_vector_np (:)
  REAL, ALLOCATABLE ::  mag_vector_sp (:)
  REAL, ALLOCATABLE ::  dir_vector_sp (:)
  
END MODULE polar_wind_mod
