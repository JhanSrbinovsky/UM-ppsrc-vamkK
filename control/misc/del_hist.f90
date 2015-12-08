! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: DEL_HIST -------------------------------------------
!
!  Purpose: delete a history file -called if problems writing
!           out partial sums.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc


      SUBROUTINE del_hist(unit)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE MPPIO_file_utils, ONLY :                                      &
          file_delete

      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: unit   ! history file unit number
      INTEGER                       :: icode  ! error from get_file
      CHARACTER(LEN=filenamelength) :: filename

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('DEL_HIST',zhook_in,zhook_handle)
      CALL get_file(unit,filename,filenamelength,icode)
      CALL file_delete(TRIM(filename))
      IF (lhook) CALL dr_hook('DEL_HIST',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE DEL_HIST
!
