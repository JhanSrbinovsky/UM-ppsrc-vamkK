! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Description:
!  This provides a wrapper allowing 32 bit writes to be made 
!  from a supplied 64 bit buffer. Provided to support legacy code.

SUBROUTINE buffout32_f77(unit,array,numWords,wordsWritten,errorCode)
  USE UM_Types, ONLY : real32
  USE IO
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
  INTEGER, INTENT(IN)  :: numWords     ! No. of words to read
  REAL(KIND=real32)    :: array(numWords)! Array to read into
  REAL,    INTENT(OUT) :: errorCode    ! Return code
  INTEGER, INTENT(OUT) :: wordsWritten ! No. of words read
  INTEGER(KIND=jpim), &
      PARAMETER        :: zhook_in  = 0
  INTEGER(KIND=jpim), &
      PARAMETER        :: zhook_out = 1
  REAL(KIND=jprb)      :: zhook_handle

  IF (lhook) CALL dr_hook('BUFFOUT32_F77',zhook_in,zhook_handle)
  CALL buffo32(unit,array(1:numWords),numWords,wordsWritten,errorCode)
  IF (lhook) CALL dr_hook('BUFFOUT32_F77',zhook_out,zhook_handle)

END SUBROUTINE buffout32_f77
