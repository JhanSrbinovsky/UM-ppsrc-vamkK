! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing error codes for Fieldcalc

MODULE Err_Mod

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! Error codes
INTEGER, PARAMETER :: StatusOK      =  0
INTEGER, PARAMETER :: StatusWarning = -9
INTEGER, PARAMETER :: StatusFatal   =  9
INTEGER, PARAMETER :: EndofFile     = -1

END MODULE Err_Mod
