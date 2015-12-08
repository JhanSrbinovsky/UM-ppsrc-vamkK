! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE aodtype_mod

IMPLICIT NONE
! Start AODTYPE

! Description:
!   Sets aerosol type numbers for aerosol optical depth diags.
!
! Current Code Owner: Nicolas Bellouin
!
! A given aerosol type may gather several aerosol components
! (see AERCMP3A)
!
      INTEGER, PARAMETER :: IP_TYPE_ALLAOD   = 0
      INTEGER, PARAMETER :: IP_TYPE_SULPHATE = 1
      INTEGER, PARAMETER :: IP_TYPE_DUST     = 2
      INTEGER, PARAMETER :: IP_TYPE_SEASALT  = 3
      INTEGER, PARAMETER :: IP_TYPE_SOOT     = 4
      INTEGER, PARAMETER :: IP_TYPE_BIOMASS  = 5
      INTEGER, PARAMETER :: IP_TYPE_BIOGENIC = 6
      INTEGER, PARAMETER :: IP_TYPE_OCFF     = 7
      INTEGER, PARAMETER :: IP_TYPE_DELTA    = 8
      INTEGER, PARAMETER :: IP_TYPE_NITRATE  = 9
      INTEGER, PARAMETER :: IP_TYPE_TWOBDUST = 10
! End AODTYPE

END MODULE aodtype_mod
