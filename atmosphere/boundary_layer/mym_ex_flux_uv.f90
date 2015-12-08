! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_ex_flux_uv-----------------------------------------
!
!  Purpose: To calculate momentum fluxes in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_ex_flux_uv(                                              &
        i_start, i_end, j_start, j_end, offx, offy, bl_levels,          &
        rdz_u_v, rhokm_u_v, rhogamuv_uv, u_v, tau_xy_fd_uv,             &
        tau_x_y, tau_grad, tau_count_grad)

  USE bl_option_mod, ONLY: formdrag, explicit_stress
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     i_end,                                                             &
                   ! Local number of points on a row
     j_end,                                                             &
                   ! Local number of rows in a theta field
     i_start,j_start,                                                   &
                   ! Starting index of arrays
     offx,                                                              &
                   ! Local number of rows in a theta field
     offy,                                                              &
                   ! Local number of rows in a theta field
     bl_levels
                   ! Max. no. of "boundary" levels
   
  REAL, INTENT(IN) ::                                                   &
     rdz_u_v (i_start:i_end,j_start:j_end, 2:bl_levels),            &
                   ! Reciprocal of the vertical
                   ! distance from level K-1 to
                   ! level K. (K > 1) on wind levels
     rhokm_u_v (i_start:i_end,j_start:j_end, bl_levels),            &
                   ! Exchange coefficients for
                   ! momentum, on UV-grid with
                   ! first and last j_end ignored.
                   ! for K>=2, between rho level K and K-1.
                   ! i.e. assigned at theta level K-1
     rhogamuv_uv(i_start:i_end, j_start:j_end, 2:bl_levels),        &
                   ! Counter Gradient Term for U or V
                   ! defined on UV-grid
     u_v(i_start-offx:i_end+offx,j_start-offy:j_end+offy,bl_levels),&
                   ! Westerly_Southerly component of wind.
     tau_xy_fd_uv(i_start:i_end,j_start:j_end, bl_levels)
                   ! X/Y-component of form-drag stress
                   !    at a UV point

! Intent INOUT Variables
  REAL, INTENT(INOUT) ::                                                &
     tau_x_y (i_start:i_end,j_start:j_end, bl_levels)
                   ! explicit x_y-component of
                   ! turbulent stress at levels
                   ! k-1/2; eg. TAUX(,1) is surface
                   ! stress. UV-grid, 1st and last j_end
                   ! set to "missing data". (N/sq m)

! Intent OUT Variables
  REAL, INTENT(OUT) ::                                                  &
     tau_grad(i_start:i_end,j_start:j_end,bl_levels),               &
                   ! k*du/dz grad stress (kg/m/s2)
     tau_count_grad(i_start:i_end,j_start:j_end,bl_levels)
                   ! Counter gradient stress (kg/m/s2)

! LOCAL VARIABLES.
  INTEGER ::                                                            &
     i, j, k

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_EX_FLUX_UV',zhook_in,zhook_handle)

  k=1
  DO j = j_start, j_end
    DO i = i_start, i_end
      tau_grad(i,j,k) = 0.0
      tau_count_grad(i,j,k) = 0.0
    END DO
  END DO

  DO k = 2, bl_levels
    DO j = j_start, j_end
      DO i = i_start, i_end

        tau_grad(i,j,k) = rhokm_u_v(i,j,k) *                            &
                       ( u_v(i,j,k) - u_v(i,j,k-1) ) *rdz_u_v(i,j,k)
        tau_count_grad(i,j,k) = rhogamuv_uv(i, j, k)
        tau_x_y(i,j,k) = tau_grad(i,j,k) + tau_count_grad(i,j,k)

      END DO
    END DO
  END DO

! Add explicit orographic stress, noting that the surface stress
! is to be added later

  IF (formdrag  ==  explicit_stress) THEN
    DO k = 2, bl_levels
      DO j = j_start, j_end
        DO i = i_start, i_end
          tau_x_y(i,j,k) = tau_x_y(i,j,k) + tau_xy_fd_uv(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook('MYM_EX_FLUX_UV',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_ex_flux_uv
