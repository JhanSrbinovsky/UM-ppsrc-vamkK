! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_closeinput_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_closeinput()

! Closes the input RPSEED stream if it is open


 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE domain_params
 USE stochastic_physics_run_mod,  ONLY: l_stphseed_read, l_stphseed_write
 IMPLICIT NONE

! Local Variables
 INTEGER :: icode ! Error return code

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_CLOSEINPUT',zhook_in,zhook_handle)
 icode=0

 CLOSE(UNIT=149,IOSTAT=icode)
 IF(icode > 0) THEN
   WRITE(6,*) "RPSEED IN : Failed to close input file"
 ELSE
   WRITE(6,*) "RPSEED IN : Input file closed"
 END IF

 IF (lhook) CALL dr_hook('STPH_CLOSEINPUT',zhook_out,zhook_handle)
 RETURN
END SUBROUTINE stph_closeinput

END MODULE stph_closeinput_mod
