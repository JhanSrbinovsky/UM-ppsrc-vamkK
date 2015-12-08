! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine ABORT_IO-----------------------------------------------
!
!    Purpose:  Prints out message and stops execution of program.
!              Called if ICODE  /=  0
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Dump I/O

SUBROUTINE abort_io(string,cmessage,icode,nft)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
IMPLICIT NONE

INTEGER                                                           &
 icode                                                            &
         !IN Code returned by UM routines
,nft     !IN Unit no being processed

CHARACTER(LEN=*)                                                     &
 string                                                           &
         !IN Subroutine name and position
,cmessage!IN Message returned by UM routines

INTEGER, PARAMETER :: stdout = 0  ! Standard output unit number.
INTEGER, PARAMETER :: stderr = 6  ! Standard error unit number.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('ABORT_IO',zhook_in,zhook_handle)
WRITE(stdout,*) 'Processor ',mype,' calling ABORT'
WRITE(stderr,*) 'Processor ',mype,' calling ABORT'

WRITE(stdout,'('' Error detected in subroutine '',A)') string
WRITE(stderr,'('' Error detected in subroutine '',A)') string

IF (nft /= 0) THEN
  WRITE(stdout,'('' while doing i/o on unit'',I3)') nft
  WRITE(stderr,'('' while doing i/o on unit'',I3)') nft
END IF
WRITE(stdout,'(A)') TRIM(cmessage)
WRITE(stderr,'(A)') TRIM(cmessage)

WRITE(stdout,'('' ICODE='',I6)') icode

IF (lhook) CALL dr_hook('ABORT_IO',zhook_out,zhook_handle)
CALL gc_abort(mype,nproc,cmessage)

CALL ereport('ABORT_IO', icode, cmessage)

IF (lhook) CALL dr_hook('ABORT_IO',zhook_out,zhook_handle)

END SUBROUTINE abort_io

