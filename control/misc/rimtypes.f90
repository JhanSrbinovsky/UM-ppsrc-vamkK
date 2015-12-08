! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP.
!
! PARAMETERs defining the RIM_TYPE characteristics

MODULE rimtypes
  IMPLICIT NONE
  INTEGER, PARAMETER :: Nrima_max = 1
  
  INTEGER, PARAMETER :: rima_type_norm=1  ! Normal field
  INTEGER, PARAMETER :: rima_type_orog=1  ! Orography field
  
END MODULE rimtypes
