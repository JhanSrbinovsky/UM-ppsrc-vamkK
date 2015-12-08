! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STACCUM  -------------------------------------------------
!LL
!LL  Purpose: Accumulates fields within STASH (temporal service routine)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D721
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE STACCUM(fieldin,result,size,masking,amdi)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
      INTEGER size              ! IN size of fieldin and result.
      REAL fieldin(size)        ! IN  input field to be processed
      REAL result(size)         ! OUT where accum is done.
      LOGICAL masking           ! IN true if masked (ie. MDI possible)
      REAL amdi                 ! IN missing data indicator
!* ---------------------------------------------------------------------
!
!  Local variables
!
      INTEGER i                 ! loop count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!L----------------------------------------------------------------------
!L 1.1 loop over array size, if either result or fieldin is amdi, set
!L    result to amdi else accumulate it in result.
!L
      IF (lhook) CALL dr_hook('STACCUM',zhook_in,zhook_handle)
      IF (masking) THEN
        DO i=1,size
          IF ((result(i) /= amdi).and.(fieldin(i) /= amdi)) THEN
            result(i)=result(i)+fieldin(i)
          ELSE
            result(i)=amdi
          ENDIF
        ENDDO
      ELSE
!L
!L 1.2 loop over array size, accumulate result without checking for
!L     missing data
!L
        DO i=1,size
          result(i)=result(i)+fieldin(i)
        ENDDO
      ENDIF
!
      IF (lhook) CALL dr_hook('STACCUM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STACCUM
