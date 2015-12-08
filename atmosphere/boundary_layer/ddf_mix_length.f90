! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE ddf_mix_length------------------------------------------
!
!  Purpose: To calculate the mixing length in the first order closure
!           model
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE ddf_mix_length(                                              &
      row_length, rows, halo_i, halo_j, bl_levels,                      &
      z_uv, z_tq, dbdz, r_mosurf, fb_surf, h_pbl, e_trb,                &
      elm, coef_ce, ekw)

  USE mym_option_mod, ONLY: tke_dlen,                                   &
              my_length, ddf_length, non_local_like_length,             &
              l_tke_dlen_blackadar, tke_levels              
  USE atmos_constants_mod, ONLY : vkman        
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent In Variables
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
     z_uv(row_length,rows,bl_levels),                                   &
                    ! Z_UV(*,K) is height of u level k
     z_tq(row_length,rows,bl_levels),                                   &
                    ! IN Z_TQ(*,K) is height of theta level k.
                    ! Cloud ice (kg per kg air)
     dbdz(row_length,rows,tke_levels),                                  &
                    ! Buoyancy gradient across layer
                    ! interface interpolated to theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     r_mosurf(row_length, rows),                                        &
                    ! reciprocal of Monin-Obkhov length
     fb_surf(row_length,rows),                                          &
                    ! Surface flux buoyancy over density (m^2/s^3)
     h_pbl(row_length, rows),                                           &
                    ! height of PBL determined by vertical profile
                    ! of SL
     e_trb(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,             &
                                                    bl_levels)
                    ! TKE defined on theta levels K-1

  REAL, INTENT(OUT) ::                                                  &
     elm(row_length, rows, tke_levels),                                 &
                    ! mixing length
     coef_ce(row_length, rows, tke_levels),                             &
                    ! coefficient appeared in a dissipation term
     ekw(row_length, rows, tke_levels)
                    ! SQRT(e_trb)

  ! Local variables
  INTEGER :: i, j, k
                   ! loop counter

  REAL ::                                                               &
     rbv,                                                               &
                    ! reciprocal of Brunt-Vaisala frequency
     elb,                                                               &
                    ! mixing length driven by buoyancy
     els,                                                               &
                    ! mixing length driven by surface
     delta_z
                    ! vertical grid spacing

  REAL ::                                                               &
     qke(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
                                                    bl_levels),         &
                    ! twice of TKE (denoted to q**2) on theta level K-1
     qkw(row_length, rows, tke_levels)
                    ! q=sqrt(qke) on theta level K-1

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('DDF_MIX_LENGTH',zhook_in,zhook_handle)

  DO k = 2, tke_levels
    DO j = 1, rows
      DO i = 1, row_length
        ekw(i, j, k) = SQRT(MAX(e_trb(i, j, k), 1.0e-20))
      END DO
    END DO
  END DO

  IF (tke_dlen == my_length) THEN
    DO k = 2, tke_levels
      DO j = 1, rows
        DO i = 1, row_length
          qke(i, j, k) = 2.0 * e_trb(i, j, k)
        END DO
      END DO
    END DO
! DEPENDS ON: mym_length
    CALL mym_length(                                                    &
          row_length, rows, halo_i, halo_j, bl_levels,                  &
          qke, z_uv, z_tq, dbdz, r_mosurf, fb_surf,                     &
          qkw, elm)
  ELSE IF (tke_dlen == ddf_length                                       &
     .OR. tke_dlen == non_local_like_length) THEN
    DO k = 2, tke_levels
      DO j = 1, rows
        DO i = 1, row_length
          delta_z = z_uv(i, j, k) - z_uv(i, j, k - 1)
          IF (dbdz(i, j, k) > 0.0) THEN
            rbv = 1.0 / SQRT(dbdz(i, j, k))
            elb = MAX(MIN(0.76 * ekw(i, j, k) * rbv,                    &
                          delta_z), 1.e-10)
          ELSE
            elb = delta_z
          END IF
          elm(i, j, k) = elb
        END DO
      END DO
    END DO

    IF (tke_dlen == non_local_like_length) THEN
      DO k = 2, tke_levels
        DO j = 1, rows
          DO i = 1, row_length
            IF (z_tq(i, j, k - 1) < h_pbl(i, j) ) THEN
              elm(i, j, k) = 0.25 * 1.8 * h_pbl(i, j)                   &
                     * (1.0 - EXP(                                      &
                               -4.0 * z_tq(i, j, k - 1)/h_pbl(i, j))    &
                        - 0.0003 * EXP(                                 &
                              8.0 * z_tq(i, j, k - 1) / h_pbl(i, j)))
            END IF
          END DO
        END DO
      END DO
    END IF  ! if tke_dlen == non_local_like_length

    IF (l_tke_dlen_blackadar) THEN
      DO k = 2, tke_levels
        DO j = 1, rows
          DO i = 1, row_length
            els = vkman * z_tq(i, j, k - 1)
            elm(i, j, k) = els / (1.0 + els / elm(i, j, k))
          END DO
        END DO
      END DO
    END IF
  END IF

  ! for diagnostics
  DO j = 1, rows
    DO i = 1, row_length
      elm(i, j, 1) = elm(i, j, 2)
    END DO
  END DO

  IF (tke_dlen == non_local_like_length) THEN
    DO k = 2, tke_levels
      DO j = 1, rows
        DO i = 1, row_length
          coef_ce(i, j, k) = 0.41
        END DO
      END DO
    END DO
  ELSE
    DO k = 2, tke_levels
      DO j = 1, rows
        DO i = 1, row_length
          coef_ce(i, j, k) = 0.19 + 0.74 * elm(i, j, k)                 &
                                / (z_uv(i, j, k) - z_uv(i, j, k - 1))
        END DO
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook('DDF_MIX_LENGTH',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ddf_mix_length
