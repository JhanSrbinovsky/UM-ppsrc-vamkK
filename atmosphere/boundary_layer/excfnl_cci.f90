! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE EXCFNL_CCI---------------------------------------------
!
!  Purpose: Compute Compressed Index
!
!  Code Description:
!    Language: FORTRAN 77 + common extensions.
!    This code is written to UMDP3 v6 programming standards.
!
!  Documentation: UMDP No.24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!--------------------------------------------------------------------
SUBROUTINE excfnl_cci (                                                 &
 c_len, to_do, ind_todo                                                 &
 )

  USE atm_fields_bounds_mod, ONLY: pdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: c_len
  LOGICAL, INTENT(INOUT) :: to_do(pdims%i_end*pdims%j_end)
  INTEGER, INTENT(INOUT) :: ind_todo(pdims%i_end*pdims%j_end)

! local variables
  INTEGER                      :: n,m

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('EXCFNL_CCI',zhook_in,zhook_handle)
  m = 0
!          Compress index for main loops
!CDIR Nodep
  DO n = 1, c_len
    IF(to_do(n))   THEN
      m=m+1
      to_do(m)    = to_do(n)
      ind_todo(m) = ind_todo(n)
    END IF
  END DO
  c_len = m

  IF (lhook) CALL dr_hook('EXCFNL_CCI',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE excfnl_cci
