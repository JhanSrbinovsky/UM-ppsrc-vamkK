! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!       Routines that are used to manage a fortran unit file stream 
!       which records a list of iteration counts.
!
!       Code Owner: See Unified Model Code Owners HTML page
!       This file belongs in section: Top Level

!       Subroutine ICCloseOutput
!       ------------------------
!       Closes the output stream if it is open

      SUBROUTINE ICCloseOutput()


          USE yomhook, ONLY: lhook, dr_hook
          USE parkind1, ONLY: jprb, jpim
          USE domain_params
          USE um_input_control_mod,  ONLY: l_icount
          IMPLICIT NONE

!       Local Variables
        integer :: icode ! Error return code

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('ICCLOSEOUTPUT',zhook_in,zhook_handle)
        if (l_icount) then
          close(unit=152,iostat=icode)
          l_icount = .false.
          if (icode > 0) then
            write(6,*) "ITERCOUNT : Failed to close output file"
          end if
        end if
        IF (lhook) CALL dr_hook('ICCLOSEOUTPUT',zhook_out,zhook_handle)
        RETURN
      END SUBROUTINE ICCloseOutput
