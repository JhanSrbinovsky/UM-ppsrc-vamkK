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

!      Subroutine ICWriteEntry(value)
!      ------------------------------
!      This routine is passed the number of iterations in its
!      argument and it then writes this to the output file if
!      it is open.

      SUBROUTINE ICwriteEntry(value)

        USE yomhook, ONLY: lhook, dr_hook
        USE parkind1, ONLY: jprb, jpim
        USE domain_params
        USE um_input_control_mod,  ONLY: l_icount
        IMPLICIT NONE
!      Arguments
        integer, intent(in) :: value ! The number of iterations

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

!      If the unit is open the write the value to the file stream
        IF (lhook) CALL dr_hook('ICWRITEENTRY',zhook_in,zhook_handle)
        if (l_icount) then
          write(152,*) value
        end if
        IF (lhook) CALL dr_hook('ICWRITEENTRY',zhook_out,zhook_handle)
        RETURN
      END SUBROUTINE ICwriteEntry
