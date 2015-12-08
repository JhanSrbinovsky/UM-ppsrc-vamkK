! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Data module containing vertical interpolation control variables

Module Rcf_V_Int_Ctl_Mod

! Description:
!  Contains variables used for control of vertical interpolation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE missing_data_mod, ONLY : imdi

IMPLICIT NONE

LOGICAL                 :: v_int_active
LOGICAL                 :: v_int_active_soil   ! switch for soil interp
INTEGER                 :: v_int_order = imdi

End Module Rcf_V_Int_Ctl_Mod
