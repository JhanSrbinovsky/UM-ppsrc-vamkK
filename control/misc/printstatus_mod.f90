! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Define variables for level of output
! Description:
!   Variables to control timer and level of standard output
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc.

MODULE PrintStatus_mod

  IMPLICIT NONE

  INTEGER, PARAMETER   :: PrStatus_Min    = 1  ! Minimum output
  INTEGER, PARAMETER   :: PrStatus_Normal = 2  ! Short informative output
  INTEGER, PARAMETER   :: PrStatus_Oper   = 3  ! Full informative output
  INTEGER, PARAMETER   :: PrStatus_Diag   = 4  ! Extra Diagnostic output
  
  ! The variable that controls message output - default setting
  ! is to 2 - Full informative output.
  
  INTEGER, SAVE        :: PrintStatus     = PrStatus_Normal


END MODULE PrintStatus_mod

