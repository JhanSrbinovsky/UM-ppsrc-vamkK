! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine so there is something to compile for section A33_1A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Free Tracers
      SUBROUTINE DUMMY_TRACER(ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

      INTEGER                       ::  ICODE
      CHARACTER (Len=*)             ::  CMESSAGE
      CHARACTER (Len=*), Parameter  ::  RoutineName='DUMMY_TRACER'

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('DUMMY_TRACER',zhook_in,zhook_handle)
      Write(6,*)'**ERROR**: DUMMY_TRACER should not be callable'
      CMESSAGE = 'Routine DUMMY_TRACER should not be callable'
      ICODE = 1

      CALL EReport(RoutineName,ICODE,CMESSAGE)

      IF (lhook) CALL dr_hook('DUMMY_TRACER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DUMMY_TRACER
