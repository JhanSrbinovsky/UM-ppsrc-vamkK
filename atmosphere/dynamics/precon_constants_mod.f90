! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! GCR solver names

! Description:
!   Module containing GCR parameter options
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

! Method:
! dynamics code options (magic number names)

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE precon_constants_mod

IMPLICIT NONE

INTEGER, PARAMETER ::      no_precon                = 0 
INTEGER, PARAMETER ::      vert_precon              = 1
INTEGER, PARAMETER ::      vert_plus_xyz_ADI_precon = 2
INTEGER, PARAMETER ::      xyz_ADI_precon           = 3
INTEGER, PARAMETER ::      vert_plus_xz_ADI_precon  = 4
INTEGER, PARAMETER ::      xz_ADI_precon            = 5
INTEGER, PARAMETER ::      Dufort_Frankel_precon    = 6
  
END MODULE precon_constants_mod
