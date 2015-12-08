! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Convection history prognostics - constants

MODULE cv_hist_constants_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing constants used in the calculation of the
!   convective history prognostics.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
!------------------------------------------------------------------------------

! Decay period for deep convection flag

REAL, PARAMETER :: decay_period = 10800. ! 3 hours (s)


!------------------------------------------------------------------------------

END MODULE cv_hist_constants_mod
