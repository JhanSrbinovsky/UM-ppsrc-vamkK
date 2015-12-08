! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_length----------------------------------------------
!
!  Purpose: To calculate mixing length in the MY model.
!           The square root of TKE, required by this subroutine and
!           elsewhere, is also returned
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_length(                                                  &
      row_length, rows, halo_i, halo_j, bl_levels,                      &
      qke, z_uv, z_tq, dbdz, r_mosurf, fb_surf,                         &
      qkw, el)

  USE atmos_constants_mod, ONLY: vkman
  USE mym_const_mod, ONLY: my_alpha4, one_third, elt_min, my_alpha1,    &
                           my_alpha2, my_alpha3
  USE mym_option_mod, ONLY: tke_levels, my_z_limit_elb
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     row_length,                                                        &
                   ! Local number of points on a row
     rows,                                                              &
                   ! Local number of rows in a theta field
     halo_i,                                                            &
                   ! Size of halo in i direction.
     halo_j,                                                            &
                   ! Size of halo in j direction.
     bl_levels
                   ! Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
     qke(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
                                                          bl_levels),   &
                   ! twice of TKE (denoted to q**2) on theta level K-1
     z_uv(row_length,rows,bl_levels),                                   &
                   ! Z_UV(*,K) is height of rho level k
     z_tq(row_length,rows,bl_levels),                                   &
                   ! Z_TQ(*,K) is height of theta level k.
     dbdz(row_length,rows,2:tke_levels),                                &
                   ! Buoyancy gradient across layer
                   ! interface interpolated to theta levels.
                   ! (:,:,K) represents the value on theta level K-1
     r_mosurf(row_length, rows),                                        &
                   ! reciprocal of Monin-Obukhov Length
     fb_surf(row_length,rows)
                   ! Surface buoyancy flux over
                   ! density (m^2/s^3)

! Intent OUT Variables
  REAL, INTENT(OUT) ::                                                  &
     qkw(row_length, rows, tke_levels),                                 &
                   ! q=sqrt(qke) on theta level K-1
     el(row_length, rows, tke_levels)
                   ! mixing length on theta level K-1

! Local variables

  INTEGER ::                                                            &
     i, j, k
                   ! Loop indexes
  REAL ::                                                               &
     qdz,                                                               &
                   ! q times vertical grid space
     alp32,                                                             &
                   ! combined constants (alpha3 / alpha2)
     rbv,                                                               &
                   ! reciprocal of Brunt-Vaisala frequency
     elb,                                                               &
                   ! mixing length related to buoyancy (L_B)
     els,                                                               &
                   ! mixing length related to surface (L_S)
     zeta
                   ! non-dimensional length (height over MO length)
  REAL ::                                                               &
     elt(row_length, rows),                                             &
                   ! mixing length related to vertical distribution
                   ! of TKE (L_T)
     vsc(row_length, rows)
                   ! work arrays


  REAL, PARAMETER ::                                                    &
     zmax = 1.0,                                                        &
                  ! constant used in calculating els
     cns = 2.7
                  ! constant used in calculating els

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_LENGTH',zhook_in,zhook_handle)

  DO j = 1, rows
    DO i = 1, row_length
      elt(i, j) = 0.0
      vsc(i, j) = 0.0
    END DO
  END DO

  DO k = 1, tke_levels
    DO j = 1, rows
      DO i = 1, row_length
        qkw(i, j, k) = SQRT(MAX(qke(i, j, k), 1.e-20))
      END DO
    END DO
  END DO

  ! vertical integration of qz and q
  ! Here, elt is still vertical integration of qz
  ! and vsc is that of q
  DO k = 2, tke_levels
    DO j = 1, rows
      DO i = 1, row_length
        qdz = qkw(i, j, k) * (z_uv(i, j, k) - z_uv(i, j, k - 1))
        elt(i, j) = elt(i, j) + qdz * z_tq(i, j, k - 1)
        vsc(i, j) = vsc(i, j) + qdz
      END DO
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
      elt(i, j) = MAX(my_alpha1 * elt(i, j) / (vsc(i, j) + 1.e-10),     &
                      elt_min)
      vsc(i, j) = (elt(i, j) * MAX(fb_surf(i, j), 0.0)) ** one_third
    END DO
  END DO

  alp32 = my_alpha3 / my_alpha2
  DO k = 2, tke_levels
    DO j = 1, rows
      DO i = 1, row_length
        IF (dbdz(i, j, k) > 0.0) THEN
          rbv = 1.0 / SQRT(dbdz(i, j, k))
          elb = my_alpha2 * qkw(i, j, k) * rbv                          &
                    * (1.0 + alp32 * SQRT(vsc(i, j) * rbv / elt(i, j)))
        ELSE
          elb = 1.0e10
        END IF

        IF (z_tq(i, j, k - 1) > my_z_limit_elb) THEN
          elb = MIN(elb, z_uv(i, j, k) - z_uv(i, j, k - 1))
        END IF

        zeta = z_tq(i, j, k - 1) * r_mosurf(i, j)
        IF (zeta > 0.0) THEN
          els = vkman * z_tq(i, j, k - 1)                               &
                     / (1.0 + cns * MIN(zeta, zmax))
        ELSE
          els = vkman * z_tq(i, j, k - 1)                               &
                     * MIN((1.0 - my_alpha4 * zeta) ** 0.2, 2.0)
        END IF
        el(i, j, k) = elb / ( elb / elt(i, j) + elb / els + 1.0)
      END DO
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
      el(i, j, 1) = el(i, j, 2)
    END DO
  END DO
  IF (lhook) CALL dr_hook('MYM_LENGTH',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_length

