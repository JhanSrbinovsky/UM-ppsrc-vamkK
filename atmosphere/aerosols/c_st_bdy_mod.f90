! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_st_bdy_mod

IMPLICIT NONE
! Start c_st_bdy
!
! Contains resistance factors for dry deposition of 3 modes of soot
!
!
! Rb(soot modes)/Rb(H2O)
       REAL, PARAMETER :: RESB_FreshSoot = 1840.77
       REAL, PARAMETER :: RESB_AgedSoot = 1840.77
       REAL, PARAMETER :: RESB_SootInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_Soot = 0.0

! End c_st_bdy

END MODULE c_st_bdy_mod
