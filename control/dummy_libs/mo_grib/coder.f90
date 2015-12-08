
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Dummy routine to resolve externals

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
SUBROUTINE coder()
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

INTEGER                       ::  icode
CHARACTER (len=80)            ::  cmessage
CHARACTER (len=* ), PARAMETER ::  routinename='CODER'

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('CODER',zhook_in,zhook_handle)

cmessage = 'Routine not available.  Please check options.'
icode = 1
CALL ereport(routinename,icode,cmessage)

IF (lhook) CALL dr_hook('CODER',zhook_out,zhook_handle)
RETURN
END SUBROUTINE coder
