! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -------------------- SUBROUTINE UM_FORT_FLUSH ---------------------
!  
!   Purpose: A wrapper script for the flush intrinsic
!
!   -------------------------------------------------------------------
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Misc

      SUBROUTINE UM_FORT_FLUSH(lunit,icode)

!     Required if using the NAG compiler




      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE um_types
      Implicit None

!     The subroutine's arguments, whatever the compiler
      INTEGER, INTENT(IN)  :: lunit
      INTEGER, INTENT(OUT) :: icode
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      !32-bit arguments needed for some platforms.






      IF (lhook) CALL dr_hook('UM_FORT_FLUSH',zhook_in,zhook_handle)

!     If on NAG or X1, require two 32 bit arguments to flush.
      CALL flush(lunit)
      icode = 0
      IF (lhook) CALL dr_hook('UM_FORT_FLUSH',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UM_FORT_FLUSH
