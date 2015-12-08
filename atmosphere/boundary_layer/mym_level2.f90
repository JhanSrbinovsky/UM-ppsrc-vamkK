! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_level2----------------------------------------------
!
!  Purpose: To calculate fundamental values such as non-dimensional
!           diffusion coefficients in the MY model with level 2
!           scheme.
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
SUBROUTINE mym_level2(                                                  &
      bl_levels, dbdz, dvdzm, gm, gh, sm, sh)

  USE atm_fields_bounds_mod, ONLY: tdims
  USE mym_const_mod, ONLY: ri1, ri2, ri3, ri4, rfc, rf1, rf2, shc, smc
  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                    ! Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
     dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:tke_levels),                                                &
                    ! Buoyancy gradient across layer
                    ! interface interpolated to theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels)
                    ! Modulus of wind shear at theta levels.
                    ! (:,:,K) repserents the value on theta level K-1

! Intent OUT Variables
  REAL, INTENT(OUT) ::                                                  &
     gm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! square of wind shear on theta level K-1
                   ! (a denominator of gradient Richardson number)
     gh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! - buoyancy gradient on theta level K-1
                   ! (a numerator of gradient Richardson number)
     sm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! Non-dimensional diffusion coefficients for
                   ! momentum from level 2 scheme
                   ! defined on theta level K-1
     sh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels)
                   ! Non-dimensional diffusion coefficients for
                   ! scalars from level 2 scheme
                   ! define on theta level K-1
! Local variables
  INTEGER ::                                                            &
     i, j, k
                   ! Loop indexes

  REAL ::                                                               &
     ri,                                                                &
                   ! gradient Richardson Number
     rf
                   ! flux Richardson Number

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_LEVEL2',zhook_in,zhook_handle)

  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        gm(i, j, k) =   dvdzm(i, j, k) * dvdzm(i, j, k)
        gh(i, j, k) = - dbdz(i, j, k)
        !   Gradient Richardson number
        ri = - gh(i, j, k) / MAX( gm(i, j, k), 1.0e-10 )
        !   Flux Richardson number
        rf = MIN(ri1 * (ri + ri2 - SQRT(ri ** 2 - ri3 * ri + ri4)),     &
                 rfc )
        sh(i, j, k) = shc * (rfc - rf) / (1.0 - rf)
        sm(i, j, k) = smc * (rf1 - rf) / (rf2 - rf) * sh(i, j, k)
      END DO
    END DO
  END DO

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      gm(i, j, 1) = 0.0
      gh(i, j, 1) = 0.0
      sh(i, j, 1) = 0.0
      sm(i, j, 1) = 0.0
    END DO
  END DO

  IF (lhook) CALL dr_hook('MYM_LEVEL2',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_level2

