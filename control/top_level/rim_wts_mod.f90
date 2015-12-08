! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  rim_wts_mod

      MODULE rim_wts_mod
      IMPLICIT NONE

! Description:  Rim weights for moisture and tracers

! Method:   These are calculated in set_rim_wts via setcona. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      REAL,  ALLOCATABLE :: rim_wts(:)
 
      END MODULE rim_wts_mod
