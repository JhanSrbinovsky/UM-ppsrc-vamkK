! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate transmission and reflection coefficients.
!
! Method:
!    Straightforward.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE trans_source_coeff(n_profile                                 &
     , i_layer_first, i_layer_last                                      &
     , isolir, l_ir_source_quad                                         &
     , tau, sum, diff, lambda, sec_0                                    &
     , gamma_up, gamma_down                                             &
     , trans, reflect, trans_0, source_coeff                            &
     , nd_profile                                                       &
     , id_op_lt, id_op_lb, id_trs_lt, id_trs_lb                         &
     , nd_source_coeff                                                  &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE vectlib_mod, ONLY : exp_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , id_op_lt                                                          &
!       Topmost declared layer for optical properties
    , id_op_lb                                                          &
!       Bottom declared layer for optical properties
    , id_trs_lt                                                         &
!       Topmost declared layer for transmission coefficients
    , id_trs_lb                                                         &
!       Bottom declared layer for transmission coefficients
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last
!       Last layer to consider

! Algorithmic control
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Quadratic source in infra-red
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region

! Optical properties of the layer
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, id_op_lt: id_op_lb)                               &
!       Optical depths of layers
    , sum(nd_profile, id_op_lt: id_op_lb)                               &
!       Sum of alpha_1 and alpha_2
    , diff(nd_profile, id_op_lt: id_op_lb)                              &
!       Difference of alpha_1 and alpha_2
    , lambda(nd_profile, id_op_lt: id_op_lb)                            &
!       Lambda
    , sec_0(nd_profile)                                                 &
!       Secant of solar zenith angle
    , gamma_up(nd_profile, id_op_lt: id_op_lb)                          &
!       Basic solar coefficient for upward radiation
    , gamma_down(nd_profile, id_op_lt: id_op_lb)
!       Basic solar coefficient for downward radiation

! Transmission and reflection coefficients and coefficients for
! source terms.
  REAL (RealK), INTENT(OUT) ::                                          &
      trans(nd_profile, id_trs_lt: id_trs_lb)                           &
!       Diffuse transmission coefficient
    , reflect(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Diffuse reflection coefficient
    , trans_0(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Direct transmission coefficient
    , source_coeff(nd_profile, id_trs_lt: id_trs_lb                     &
        , nd_source_coeff)
!       Source coefficients


! Local variables
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      gamma                                                             &
!       Gamma
    , exponential                                                       &
!       Exponential of scaled optical depth
    , gamma2                                                            &
!       Gamma squared
    , exponential2
!       Exponential squared

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      sq_eps_r                                                          &
!       The square root of a real number that is negligible
!       compared to 1.0.
    , tmp_inv                                                           &
!       Temporary work variable
    , tol
!       The tolerance used for switching to the asymptotic form
!       of the quadratic source function term.

  REAL (RealK) ::                                                       &
      temp(nd_profile)
!       Temporay variable required if using exp_v


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('TRANS_SOURCE_COEFF',zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  sq_eps_r=SQRT(EPSILON(sq_eps_r))
  tol=SQRT(sq_eps_r)

! Determine the diffuse transmission and reflection coefficients.

  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
      temp(l)=-lambda(l,i)*tau(l,i)
    END DO
    CALL exp_v(n_profile, temp, temp)
    DO l=1, n_profile
      exponential=temp(l)
      exponential2=exponential*exponential
      gamma=(sum(l, i)-lambda(l, i)) / (sum(l, i)+lambda(l, i))
      gamma2=gamma*gamma
      tmp_inv=1.0_RealK / ( 1.0_RealK - exponential2*gamma2 )
      trans(l, i)=exponential*(1.0_RealK-gamma2)*tmp_inv
      reflect(l, i)=gamma*(1.0_RealK-exponential2)*tmp_inv
    END DO
  END DO



  IF (isolir == ip_solar) THEN

!   Calculate the direct transmission and the source coefficients
!   for the solar beam: in the solar case these are
!   the coefficients which will multiply the direct flux at the
!   top of the layer to give the source terms for the upward
!   diffuse flux and the total downward flux.

    DO i=i_layer_first, i_layer_last
     DO l=1, n_profile
        temp(l) = -tau(l,i)*sec_0(l)
     END DO
     CALL exp_v(n_profile,temp,trans_0(1,i))
     DO l=1, n_profile
        source_coeff(l, i, ip_scf_solar_up)                             &
          =(gamma_up(l, i)-reflect(l, i)                                &
          *(1.0_RealK+gamma_down(l, i)))                                &
          -gamma_up(l, i)*trans(l, i)*trans_0(l, i)
        source_coeff(l, i, ip_scf_solar_down)=trans_0(l, i)             &
          *(1.0_RealK+gamma_down(l, i)                                  &
          -gamma_up(l, i)*reflect(l, i))                                &
          -(1.0_RealK+gamma_down(l, i))*trans(l, i)
      END DO
    END DO


  ELSE IF (isolir == ip_infra_red) THEN

!   In the case of infra-red radiation, the first source
!   coefficient holds the multiplier for the first difference
!   of the Planckian function across the layer, and the second
!   that for the second difference.

    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile

!       A tolerance is added to the numerator and the denomiator
!       to avoid ill-conditioning at small optical depths.

        source_coeff(l, i, ip_scf_ir_1d)=(1.0_RealK-trans(l, i)         &
          +reflect(l, i)+sq_eps_r)                                      &
          /(sq_eps_r+tau(l, i)*sum(l, i))

      END DO
    END DO


    IF (l_ir_source_quad) THEN

!     Quadratic correction to source function.
!     This correction is very ill-conditioned for
!     small optical depths so the asymptotic form is then used.

      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          IF (tau(l, i) > tol) THEN
            source_coeff(l, i, ip_scf_ir_2d)                            &
              =-2.0_RealK                                               &
              *(1.0_RealK-trans(l, i)-reflect(l, i)+sq_eps_r)           &
              /(diff(l, i)*tau(l, i)+sq_eps_r)
          ELSE
            source_coeff(l, i, ip_scf_ir_2d)                            &
              =-2.0_RealK+diff(l, i)*tau(l, i)
          END IF
          source_coeff(l, i, ip_scf_ir_2d)                              &
            =-(1.0_RealK+reflect(l, i)+trans(l, i)                      &
            +source_coeff(l, i, ip_scf_ir_2d))                          &
            /(sum(l, i)*tau(l, i)+sq_eps_r)
        END DO
      END DO

    END IF

  END IF


  IF (lhook) CALL dr_hook('TRANS_SOURCE_COEFF',zhook_out,zhook_handle)

END SUBROUTINE trans_source_coeff
