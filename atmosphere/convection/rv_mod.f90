! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
MODULE rv_mod

IMPLICIT NONE

! The parameter held below in this module is ONLY use by 4A convection scheme 
! to retain bit comparison. 
! Newer versions of convection use the value of rv as set in
! atmos_constants mod which is consistent with the other UM constants.
! RV gas constant for water vapour (J/kg/K)
      REAL,PARAMETER:: RV = 461.1
! RV end

END MODULE rv_mod
