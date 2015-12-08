! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_ocff_bdy_mod

IMPLICIT NONE
! Start c_ocff_bdy
!
! Contains resistance factors for dry deposition of fossil fuel
! organic carbon aerosols
!
! Rb(smoke modes)/Rb(H2O)
       REAL, PARAMETER :: RESB_Freshocff = 2492.0
       REAL, PARAMETER :: RESB_Agedocff  = 2970.0
       REAL, PARAMETER :: RESB_ocffInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_ocff = 0.0

! End c_ocff_bdy

END MODULE c_ocff_bdy_mod
