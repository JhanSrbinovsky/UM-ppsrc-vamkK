! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_update_fields---------------------------------------
!
!  Purpose: To integrate the prognostic variables appearing
!           in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_update_fields(bl_levels,coef,dfm,prod,disp_coef,field)

  USE atm_fields_bounds_mod, ONLY: pdims, pdims_l, tdims, tdims_s
  USE timestep_mod, ONLY: timestep
  USE mym_option_mod, ONLY: l_my_extra_level, tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                   ! Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
     coef
                   ! factor for the diffusion coefficients to those for
                   ! momentum

  REAL, INTENT(IN) ::                                                   &
     dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,   &
         bl_levels),                                                    &
                   ! diffusion coefficients for momentum
     prod(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! production term
     disp_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tke_levels)
                   ! coefficients of dissipation term

! Intent INOUT Variables
  REAL, INTENT(INOUT) ::                                                &
     field(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end, &
           bl_levels)
                   ! field to be integrated

! Local variables
  INTEGER ::                                                            &
     i, j, k, k_start
                   ! Loop indexes

  REAL ::                                                               &
     aa(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),&
     bb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),&
     cc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels),&
     qq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,tke_levels)
                  ! coefficients of tri-diagonal equations

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_UPDATE_FIELDS',zhook_in,zhook_handle)

  ! Calculate the coefficients of tri-diagonal eqs. due to diffusion
! DEPENDS ON: mym_diff_matcoef
  CALL mym_diff_matcoef(bl_levels, coef, dfm, aa, bb, cc)

  IF (l_my_extra_level) THEN
    k_start = 1
  ELSE
    k_start = 2
  END IF

  DO k = k_start, tke_levels
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        aa(i, j, k) = - aa(i, j, k) * timestep
        bb(i, j, k) = 1.0 - bb(i, j, k) * timestep                      &
                            + timestep * disp_coef(i, j, k)
        cc(i, j, k) = - cc(i, j, k) * timestep
        qq(i, j, k) = field(i, j, k) + timestep * prod(i, j, k)
      END DO
    END DO
  END DO

! Solve the tri-diagonal equations
! DEPENDS ON: mym_implic
  CALL mym_implic(tke_levels, k_start, tke_levels, aa, bb, cc, qq)

  DO k = k_start, tke_levels
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        field(i, j, k) = qq(i, j, k)
      END DO
    END DO
  END DO
  IF (lhook) CALL dr_hook('MYM_UPDATE_FIELDS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_update_fields

