! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Deallocate the SW extinction diagnostics

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!    FORTRAN 95

!-----------------------------------------------------------------------

SUBROUTINE deallocate_swdiag_ext(j_sw)

  USE sw_diag_mod
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: j_sw              ! call to SW radiation

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! Extinction diagnostics

  IF (lhook) CALL dr_hook('DEALLOCATE_SWDIAG_EXT',zhook_in,zhook_handle)
  IF (ASSOCIATED(sw_diag(j_sw)%cloud_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%cloud_extinction)
  IF (ASSOCIATED(sw_diag(j_sw)%cloud_weight_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%cloud_weight_extinction)
  IF (ASSOCIATED(sw_diag(j_sw)%ls_cloud_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%ls_cloud_extinction)
  IF (ASSOCIATED(sw_diag(j_sw)%ls_cloud_weight_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%ls_cloud_weight_extinction)
  IF (ASSOCIATED(sw_diag(j_sw)%cnv_cloud_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%cnv_cloud_extinction)
  IF (ASSOCIATED(sw_diag(j_sw)%cnv_cloud_weight_extinction)) &
    DEALLOCATE(sw_diag(j_sw)%cnv_cloud_weight_extinction)

  IF (lhook) CALL dr_hook('DEALLOCATE_SWDIAG_EXT',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE deallocate_swdiag_ext
