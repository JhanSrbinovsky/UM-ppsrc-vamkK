! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected AND a full boundary layer treatment is not
!           performed.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer

SUBROUTINE bl_lsp( bl_levels,qcf,q,t )

  USE atm_fields_bounds_mod, ONLY:                                      &
   tdims, qdims

  USE atmos_constants_mod, ONLY: cp

  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::                                                &
    bl_levels             ! IN   Number of boundary layer levels

  REAL, INTENT(INOUT) ::                                                &
    qcf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,            &
        bl_levels),                                                     &
                                   ! INOUT Ice water content
    q(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,              &
      bl_levels),                                                       &
                                   ! INOUT
!                                  IN    Vapour+liquid+ice content
!                                  OUT   Vapour+liquid content
    t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      bl_levels)                   ! INOUT
!                                  IN    Liquid ice temperature
!                                  OUT   Liquid temperature
! Temporary Space
  INTEGER ::                                                            &
          i,                                                            &
                                 ! Counter over points
          j,                                                            &
                                 ! Counter over points
          k                ! Counter over boundary layer levels
  REAL :: newqcf              ! Temporary variable for QCF
  REAL :: lsrcp                 ! IN Latent heat of sublimation / Cp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  PARAMETER( lsrcp=((lc+lf)/cp) )

  IF (lhook) CALL dr_hook('BL_LSP',zhook_in,zhook_handle)
  DO k = 1, bl_levels
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
! Convert Q (vapour+liquid+ice) to (vapour+liquid)
        q(i,j,k)=q(i,j,k)-qcf(i,j,k)
! Check that Q is not negative
        IF (q(i,j,k)  <   0.0) THEN
! Evaporate ice to keep Q positive, but don't let ice go negative
! itself
          newqcf=MAX(qcf(i,j,k)+q(i,j,k),0.0)
          q(i,j,k)=q(i,j,k)+(qcf(i,j,k)-newqcf)
          qcf(i,j,k)=newqcf
        END IF
! Adjust T from T liquid ice to T liquid
        t(i,j,k)=t(i,j,k)+lsrcp*qcf(i,j,k)
      END DO
    END DO
  END DO
! End the subroutine
  IF (lhook) CALL dr_hook('BL_LSP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bl_lsp
