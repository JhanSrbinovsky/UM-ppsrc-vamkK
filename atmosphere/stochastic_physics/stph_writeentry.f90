! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_writeentry_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_writeentry( seed, seedsize)

! This routine is passed the integer random seed in its
! argument and it then writes this to the output file if
! it is open.


 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE domain_params
 USE stochastic_physics_run_mod, ONLY: l_stphseed_read, l_stphseed_write
 IMPLICIT NONE

! Arguments
 INTEGER, INTENT(IN) :: seedsize       ! The size of the random seed
 INTEGER, INTENT(IN) :: seed(seedsize) ! The random seed

! Local Variables
 INTEGER :: i                          ! loop variable

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY',zhook_in,zhook_handle)

! If the unit is open the write the value to the file stream
 DO i = 1, seedsize
   WRITE(149,*) seed(i)
 END DO

 IF (lhook) CALL dr_hook('STPH_WRITEENTRY',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE stph_writeentry
END MODULE stph_writeentry_mod
