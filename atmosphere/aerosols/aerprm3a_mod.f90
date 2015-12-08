! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE aerprm3a_mod

IMPLICIT NONE
! AERPRM3A aerosol parametrizations for two-stream radiation code.

      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_DRY=1
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_MOIST=2
      INTEGER,PARAMETER:: IP_AEROSOL_UNPARAMETRIZED=3 ! Observational
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_PHF_DRY=4
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_PHF_MOIST=5

! AERPRM3A end

END MODULE aerprm3a_mod
