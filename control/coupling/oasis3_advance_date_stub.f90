
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_ADVANCE_DATE()

  ! Description: This subroutine is a stub routine for
  !              OASIS3_ADVANCE_DATE. Run time logic
  !              should prevent it from actually being called.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Coupling
  !
  !=============================================================


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

  CHARACTER(LEN=52)   :: Message
  CHARACTER(LEN=20)   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
  !   0 = Normal exit
  ! +ve = Fatal Error
  ! -ve = Warning

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
 

  IF (lhook) CALL dr_hook('OASIS3_ADVANCE_DATE',zhook_in,zhook_handle)
  ErrorStat   = 1
  RoutineName = 'oasis3_advance_date_stub'
  Message     = 'OASIS3 Routines unavailable - see output.'

  WRITE (6,'(A)') '**ERROR**: oasis3_advance_date unavailable.'
  WRITE (6,'(A)') 'Check OASIS3 cpp key is set'

  CALL ereport(RoutineName, ErrorStat, Message)

  IF (lhook) CALL dr_hook('OASIS3_ADVANCE_DATE',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE OASIS3_ADVANCE_DATE
