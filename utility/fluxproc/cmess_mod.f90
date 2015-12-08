! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small Execs
MODULE cmess_mod

IMPLICIT NONE
!----------------------------------------------------------------------
! comdeck: CMESS
! Purpose: declares and stores values used for warning, error and
!          standard log messages
!----------------------------------------------------------------------
! declarations:
! unit numbers for debugging, error, warning and standard log output
      INTEGER, PARAMETER :: OutUnitDbg = 90  
                            ! output unit for debugging information
      INTEGER, PARAMETER :: UnWarn     = 92
      INTEGER, PARAMETER :: UnErr      = 91
      INTEGER, PARAMETER :: UnStd      = 93

! parameters for start of message
      CHARACTER(LEN=9), PARAMETER :: CWarn= 'Warning: '
      CHARACTER(LEN=7), PARAMETER :: CErr = 'ERROR: '
      CHARACTER(LEN=1), PARAMETER :: CStd = ' '

! local character variable for subroutine name
      CHARACTER(LEN=20) CSub
!----------------------------------------------------------------------

END MODULE cmess_mod
