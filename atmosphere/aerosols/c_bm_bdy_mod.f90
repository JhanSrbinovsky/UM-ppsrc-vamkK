! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE c_bm_bdy_mod

IMPLICIT NONE
! Start c_bm_bdy
!
! Contains resistance factors for dry deposition of biomass smoke
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Rb(smoke modes)/Rb(H2O)
       REAL, PARAMETER :: RESB_FreshBmass = 2492.4
       REAL, PARAMETER :: RESB_AgedBmass = 2970.7
       REAL, PARAMETER :: RESB_BmassInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_Bmass = 0.0

! End c_bm_bdy

END MODULE c_bm_bdy_mod
