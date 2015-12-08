! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate explicit flux of momentum in u or v direction
!
!  Programming standard: UMDP 3
!
!  Documentation: UM Documentation Paper No 24.
!
! SUBROUTINE EX_FLUX_UV
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE ex_flux_uv (                                                 &
  i_start,i_end, j_start,j_end,offx,offy, bl_levels,                    &
  u_v, zh, rdz_u_v, rhokm_u_v, f_ngstress_uv, tau_xy_fd_uv,             &
  tau_x_y, tau_grad,tau_non_grad                                        &
  )

  USE bl_option_mod, ONLY : max_stress_grad, ng_stress,                 &
       BrownGrant97_limited, formdrag, explicit_stress

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

  INTEGER, INTENT(IN) ::                                                &
    j_end,offx,offy,                                                    &
    i_end,i_start,j_start,                                              &
                               ! IN No. of points in latitude row.
    bl_levels
                               ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.

  REAL, INTENT(IN) ::                                                   &
    rdz_u_v (i_start:i_end,j_start:j_end, 2:bl_levels),                 &
!                                IN Reciprocal of the vertical
!                                   distance from level K-1 to
!                                   level K. (K > 1) on wind levels
    rhokm_u_v (i_start:i_end,j_start:j_end, bl_levels),                 &
!                                IN Exchange coefficients for
!                                   momentum, on UV-grid with
!                                   first and last j_end ignored.
!                                   for K>=2 (from KMKH).
    f_ngstress_uv(i_start:i_end,j_start:j_end,2:bl_levels),             &
!                                IN dimensionless function for
                               !    non-gradient wind stress,
                               !    either U or V depending on call
    u_v (i_start-offx:i_end+offx,j_start-offy:j_end+offy,bl_levels),    &
!                                IN Westerly_Southerly component of
!                                   wind.
    tau_xy_fd_uv(i_start:i_end,j_start:j_end, bl_levels),               &
                               ! IN X/Y-component of form-drag stress
                               !    at a UV point
    zh(i_start:i_end,j_start:j_end)    ! IN non-local BL depth

! INOUT variables
  REAL, INTENT(INOUT) ::                                                &
    tau_x_y (i_start:i_end,j_start:j_end, bl_levels)
!                                OUT explicit x_y-component of
!                                    turbulent stress at levels
!                                    k-1/2; eg. TAUX(,1) is surface
!                                    stress. UV-grid, 1st and last j_end
!                                    set to "missing data". (N/sq m)

! ARGUMENTS WITH INTENT OUT. IE: INPUT VARIABLES CHANGED ON OUTPUT.
  REAL, INTENT(OUT) ::                                                  &
    tau_grad(i_start:i_end,j_start:j_end,bl_levels),                    &
!                                OUT k*du/dz grad stress (kg/m/s2)
    tau_non_grad(i_start:i_end,j_start:j_end,bl_levels)
!                                OUT Non-grad stress (kg/m/s2)

  INTEGER ::                                                            &
    i,                                                                  &
    j,                                                                  &
    k,                                                                  &
    error

  REAL ::                                                               &
    tau_surf(i_start:i_end,j_start:j_end),                              &
                               ! Explicit surface stress
    sign_tau,                                                           &
                               ! Sign of surface stress
    bl_stress_grad

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('EX_FLUX_UV',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!  0.  Check that the scalars input to define the grid are consistent.
!-----------------------------------------------------------------------

  error = 0

!-----------------------------------------------------------------------
!  1.  Calculate "explicit" surface fluxes of momentum
!-----------------------------------------------------------------------

      ! Set diagnostics to zero in level 1

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, sign_tau, bl_stress_grad)

  k=1

!$OMP DO SCHEDULE(STATIC)
  DO j = j_start, j_end
    DO i = i_start, i_end
      tau_grad(i,j,k) = 0.0
      tau_non_grad(i,j,k) = 0.0
      tau_surf(i,j) = tau_x_y(i,j,1)
    END DO
  END DO
!$OMP END DO

  IF (ng_stress == BrownGrant97_limited) THEN
        ! Limit the explicitly calculated surface stress used to scale
        ! the non-gradient stress parametrization such that the implied
        ! stress gradient across the BL is less than MAX_STRESS_GRAD.
        ! This has been found by experimentation to be sufficient to
        ! stop large T increments being generated in the dynamics
        ! solver, via large increments to W

!$OMP DO SCHEDULE(STATIC)
    DO j = j_start, j_end
      DO i = i_start, i_end
        sign_tau = SIGN(1.0, tau_x_y(i,j,1) )
        bl_stress_grad = ABS( tau_x_y(i,j,1) )/zh(i,j)
        bl_stress_grad = MIN( max_stress_grad, bl_stress_grad )
        tau_surf(i,j) = sign_tau * zh(i,j) * bl_stress_grad
      END DO
    END DO
!$OMP END DO
  END IF


!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = j_start, j_end
      DO i = i_start, i_end

        tau_grad(i,j,k) = rhokm_u_v(i,j,k) *                            &
                        ( u_v(i,j,k) - u_v(i,j,k-1) ) *rdz_u_v(i,j,k)
        tau_non_grad(i,j,k) = f_ngstress_uv(i,j,k) * tau_surf(i,j)
        tau_x_y(i,j,k) = tau_grad(i,j,k) + tau_non_grad(i,j,k)

! Add explicit orographic stress, noting that the surface stress
! is to be added later
        IF (formdrag  ==  explicit_stress) THEN
          tau_x_y(i,j,k) = tau_x_y(i,j,k) + tau_xy_fd_uv(i,j,k)
        END IF

      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL 

  IF (lhook) CALL dr_hook('EX_FLUX_UV',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ex_flux_uv
