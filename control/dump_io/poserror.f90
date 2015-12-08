! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE POSERROR---------------------------------------
!
!    Purpose:
!             Prints out a message when position of a data block as
!             pointed to by fixed length header differs from actual
!             position in model dump.
!
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O

SUBROUTINE poserror(string,start_block,head_pos,head_address)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 start_block                                                      &
              !IN Actual position of data block
,head_pos                                                         &
              !IN Position in FIXHD of pointer
,head_address !IN Position in file pointed to by FIXHD(HEAD_POS)

CHARACTER(LEN=80) string  !IN Description of block

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('POSERROR',zhook_in,zhook_handle)
WRITE(6,'('' ******FATAL ERROR WHEN READING MODEL DUMP******'')')
WRITE(6,'('' Conflict between start position of '',A)')string
WRITE(6,'(A,i3,A,i9)')                                            &
 ' Block and Pointer in Fixed Length Header: fixhd(',             &
 head_pos,') =',head_address  
WRITE(6,'('' Current position in file ='',i9,'' words in'')')     &
start_block
WRITE(6,'('' ***********************************************'')')

IF (lhook) CALL dr_hook('POSERROR',zhook_out,zhook_handle)
RETURN
END SUBROUTINE poserror

