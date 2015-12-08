! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STMAX ----------------------------------------------------
!LL
!LL  Purpose: Computes the point-by-point maximum in time of a field
!LL           by comparing the field at the current time with the
!LL           maximum so far (STASH TEMPORAL service routine)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D722
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH
      SUBROUTINE STMAX(fieldin,result,size,masking,amdi)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
      INTEGER size             ! IN size of fieldin and result.
      REAL fieldin(size)       ! IN input field
      REAL result(size)        ! OUT output field (maximum)
      LOGICAL masking          ! IN true if masked (ie. MDI possible)
      REAL amdi                ! IN missing data indicator
!*----------------------------------------------------------------------
!
! Local variables
!
      INTEGER i ! loop count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------------
!L 1.1 loop over array size, if either result or fieldin is amdi set
!L     result to amdi, else set result to maximum of fieldin and result
!L
      IF (lhook) CALL dr_hook('STMAX',zhook_in,zhook_handle)
      IF (masking) THEN
        DO i=1,size
          IF ((result(i) /= amdi).and.(fieldin(i) /= amdi)) THEN
            result(i)=max(result(i),fieldin(i))
          ELSE
            result(i)=amdi
          ENDIF
        ENDDO
      ELSE
!L
!L 1.2 loop over array size, set result to maximum of fieldin and result
!L     without checking for missing data
!L
        DO i=1,size
          result(i)=max(result(i),fieldin(i))
        ENDDO
      ENDIF
!
      IF (lhook) CALL dr_hook('STMAX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STMAX
