! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!   The basic two-stream coefficients in the differential
!   equations are calculated on the assumption that scattering
!   can be ignored this routine is therefore only suitable
!   for use in the IR region. These coefficients are then used
!   to determine the transmission and reflection coefficients.
!   Coefficients for determining the solar or infra-red source
!   terms are also calculated.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE two_coeff_fast_lw(n_profile                                  &
    , i_layer_first, i_layer_last                                       &
    , l_ir_source_quad, tau                                             &
    , trans, source_coeff                                               &
    , nd_profile, nd_layer, id_lt, id_lb, nd_source_coeff               &
    )


  USE realtype_rd, ONLY: RealK
  USE vectlib_mod, ONLY: exp_v
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_lt                                                             &
!       Topmost declared layer of optical depths
    , id_lb                                                             &
!       Bottom declared layer of optical depths
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to process
    , i_layer_last
!       Last layer to process
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Use a quadratic source function

! Optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, id_lt: id_lb)
!       Optical depths


! Coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      trans(nd_profile, nd_layer)                                       &
!       Diffuse transmission coefficient
    , source_coeff(nd_profile, nd_layer, nd_source_coeff)
!       Source coefficients in two-stream equations


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  INTEGER :: n_in
!       Number of elements for vector exponential

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r
!       The square root of the above

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('TWO_COEFF_FAST_LW',zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(tau(1, id_lt))
  sq_eps_r=SQRT(eps_r)

  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
      trans(l, i)=-1.66e+00_RealK*tau(l, i)
    END DO
    DO l=n_profile+1, nd_profile
      trans(l, i)=0.0e+00_RealK
    END DO
  END DO
  n_in=nd_profile*(i_layer_last-i_layer_first+1)
  CALL exp_v(n_in                                                       &
    , trans(1, i_layer_first), trans(1, i_layer_first))

  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
      source_coeff(l, i, ip_scf_ir_1d)                                  &
        =(1.0e+00_RealK-trans(l, i)+sq_eps_r)                           &
        /(1.66e+00_RealK*tau(l, i)+sq_eps_r)
    END DO
  END DO

  IF (l_ir_source_quad) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        source_coeff(l, i, ip_scf_ir_2d)                                &
          =-(1.0e+00_RealK+trans(l, i)                                  &
          -2.0e+00_RealK*source_coeff(l, i, ip_scf_ir_1d))              &
          /(1.66e+00_RealK*tau(l, i)+sq_eps_r)
      END DO
    END DO
  END IF


  IF (lhook) CALL dr_hook('TWO_COEFF_FAST_LW',zhook_out,zhook_handle)

END SUBROUTINE two_coeff_fast_lw
