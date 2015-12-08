! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
MODULE stph_readentry_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE stph_readentry( seed, seedsize)

! This routine reads in the random seed which is the passed
! up to the calling routine as a subroutine argument.


 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim
 USE filenamelength_mod, ONLY:                                          &
     filenamelength

 USE PrintStatus_mod
 USE UM_ParVars
 USE domain_params
 USE stochastic_physics_run_mod, ONLY: l_stphseed_read, l_stphseed_write
 IMPLICIT NONE

! Arguments
 INTEGER, INTENT(IN)     :: seedsize       ! The size of the random seed
 INTEGER(8), INTENT(OUT) :: seed(seedsize) ! The random seed

! Local Variables
 INTEGER           :: i               ! loop variable
 INTEGER           :: icode           ! 0 if successful read
 LOGICAL           :: OPENED, NAMED
 INTEGER           :: rlength
 CHARACTER(LEN=filenamelength) :: fname
 CHARACTER(LEN=9)  :: dir_acc, rw_sts, r_sts

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_READENTRY',zhook_in,zhook_handle)

! If the unit is open the write the value to the file stream
 IF  (l_stphseed_read) THEN
   IF (mype == 0) THEN
     DO i=1,seedsize
       INQUIRE(UNIT=149, NAMED=NAMED, OPENED=OPENED, NAME=fname,        &
               READ=r_sts, READWRITE=rw_sts, RECL=rlength,              &
               DIRECT=dir_acc)
       IF (i == 1) THEN
         IF (printstatus  ==  prstatus_diag) THEN
           WRITE(6,*) 'Readwrite=', rw_sts
           WRITE(6,*) 'Rec length=', rlength
           WRITE(6,*) 'read=', r_sts
           WRITE(6,*) 'Direct_acces=', dir_acc
           WRITE(6,*) 'opened=', OPENED
           WRITE(6,*) 'named=', NAMED
           WRITE(6,*) 'fname=',fname
           WRITE(6,*) 'seedsize=',seedsize
         END IF  ! Print
       END IF  ! i=1
       READ(149,'(256I20)',ADVANCE='yes') seed(i)
     END DO    ! seedsize
   END IF     ! mype
 END IF       ! L_STPHSEED_READ
 IF (lhook) CALL dr_hook('STPH_READENTRY',zhook_out,zhook_handle)
 RETURN

END SUBROUTINE stph_readentry
END MODULE stph_readentry_mod
