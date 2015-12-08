! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calculate precipation for warm rain tcs convection
!
MODULE tcs_precip
  

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! Module to calculate precipation for warm rain tcs convection
  !
  ! Method:
  !   <reference to documentation to go here, once available>
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !

CONTAINS

  SUBROUTINE calc_precip( n_xx, nlev, conv_type                      &
     ,                     ql_ad, dz_inv, rho_inv                    &
     ,                     eta_rho,eta_theta                         &
     ,                     rainfall, wqr_inv, precip_product_inv     &
     ,                     wqr, precip_product_th, precip_product_uv )

    !-----------------------------------------------------------------
    !
    ! Description:
    !   Calculate the warm rain precipitation
    !
    !-----------------------------------------------------------------

    USE tcs_parameters_warm,     ONLY:                               &
       ql_t, epsilon_rain, beta_rain, gamma_rain, kappa_rain
    USE tcs_constants,           ONLY:                               &
       g
    USE tcs_common_warm ,        ONLY:                               &
       scales

    IMPLICIT NONE
    !-----------------------------------------------------------------
    ! Subroutine Arguments
    !-----------------------------------------------------------------
    INTEGER, INTENT(in) ::                                            &
       n_xx                                                           &
                                ! No. of congestus convection points
       ,  nlev                                                           
    ! No. model levels

    INTEGER, INTENT(in) ::                                            &
       conv_type(n_xx)                                  
    ! Integer index describing convective type:
    !    1=non-precipitating shallow
    !    2=drizzling shallow
    !    3=warm congestus

    REAL, INTENT(in) ::                                               &
       ql_ad(n_xx)                                                    &
                                ! adiabatic liquid water content (kg/kg)
       , dz_inv(n_xx)                                                 &
                                ! dz across inversion (m)
       , rho_inv(n_xx)                                                &
                                ! density at inversion (kg/m3)
       , eta_rho(:,:)                                                 & 
                                ! eta rho levels
       , eta_theta(:,:)         
    ! eta theta levels

    REAL, INTENT(out) ::                                              &
       rainfall(n_xx)                                                 &
                                ! surface rainfall rate
       , wqr_inv(n_xx)                                                & 
                                ! flux of w'qr' at inversion
       , precip_product_inv(n_xx)                                     &
                                ! Integral of precip_product over inv
       , wqr(n_xx,nlev)                                               & 
                                !  w'qr'
       , precip_product_th(n_xx,nlev)                                 &
                                ! precipitation production rate
       , precip_product_uv(n_xx,nlev) ! precipitation production rate


    !-----------------------------------------------------------------
    ! Variables defined locally
    !-----------------------------------------------------------------
    REAL  ::                                                          &
       FRACTION(n_xx)                                                 &
                      ! fraction of ensemble with lwc above threshold
       , ql_up_inv(n_xx)                                              &
                      ! mean liquid water content of ensemble at base
                      ! of inversion
       , p0(n_xx)                                                     &
                      ! precip productio p0 value
       , fr_inv(n_xx)    
                      ! gravitional flux of precipitation at inversion

    REAL  ::                                                          &
       fr(n_xx,nlev)  ! gravitional flux of precipitation


    INTEGER :: neta   ! number of levels in eta_rho/eta_theta
    INTEGER :: i,k    ! loop counters

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('TCS_PRECIP:CALC_PRECIP',zhook_in,zhook_handle)
    neta=SIZE(eta_rho(1,:))

    !-----------------------------------------------------------------
    ! 1.0 Initialise arrays
    !-----------------------------------------------------------------
    precip_product_th(:,:) = 0.0
    precip_product_uv(:,:) = 0.0
    fr(:,:) = 0.0
    wqr(:,:) = 0.0
    FRACTION(:) = 0.0
    rainfall(:) = 0.0
    p0(:) = 0.0
    fr_inv(:) = 0.0
    wqr_inv(:) = 0.0
    precip_product_inv(:) = 0.0

    !-----------------------------------------------------------------
    ! Calculate liquid water content of buoyant updraughts at the
    ! base of the inversion
    !-----------------------------------------------------------------
    ql_up_inv(:) = (1.5/g)*((scales%wstar_up(:)/scales%mb(:))**0.5)* &
       scales%wstar_up(:)*scales%wstar_up(:)/scales%zcld(:)

    ! Only produce precip if adiabatic liquid water content is above
    ! threshold.  The fraction of the ensemble producing rain will
    ! further depend on whether the mean ql is above the threshold.
    WHERE (ql_ad(:) > ql_t .AND. ql_up_inv(:) > ql_t)
      FRACTION(:) = 1.0 - ql_t*ql_t/(ql_up_inv(:)*ql_ad(:))
    ELSEWHERE  (ql_ad(:) > ql_t .AND. ql_up_inv(:) <= ql_t)
      FRACTION(:) = (1.0-ql_t/ql_ad(:))*(1.0-ql_t/ql_ad(:))/         &
         (1.0-ql_up_inv(:)/ql_ad(:))
    END WHERE

    WHERE (ql_ad(:) > ql_t .AND. conv_type(:) > 1) 
      ! Only produce precip if adiabatic
      ! liquid water content is above
      ! threshold 
      ! and is of precipitating type
      !-------------------------------------------------------------
      ! surface precipitation rate
      !-------------------------------------------------------------
      rainfall(:) = epsilon_rain*scales%mb(:)*ql_up_inv(:)*FRACTION(:)

      p0(:) = rainfall(:)/(beta_rain*(1.-0.5*beta_rain)              &
         *scales%zcld(:)*scales%zcld(:) +                            &
         0.5*beta_rain*scales%zcld(:)*dz_inv(:))

      !-------------------------------------------------------------
      ! fr (= downwards settling flux of rain) at inversion
      !-------------------------------------------------------------
      fr_inv(:) = gamma_rain*rainfall(:)

      !-------------------------------------------------------------
      ! wqr at inversion
      !-------------------------------------------------------------
      wqr_inv(:) = gamma_rain*rainfall(:)                            &
         - 0.5*p0(:)*beta_rain*scales%zcld(:)*dz_inv(:)   

      precip_product_inv(:) = kappa_rain*0.5*p0(:)*beta_rain       &
         *scales%zcld(:)*dz_inv(:)/rho_inv(:)
    END WHERE


    !----------------------------------------------------------------
    ! calculate in cloud profiles of precipitation
    !----------------------------------------------------------------
    !--------------------------
    ! on rho levels
    !--------------------------
    DO k=1,neta
      DO i=1,n_xx
        IF (eta_rho(i,k) <= beta_rain)THEN
          ! Lower portion of cloud layer
          precip_product_uv(i,k) = eta_rho(i,k)*scales%zcld(i)      &
             *p0(i)
          fr(i,k) = rainfall(i)
          wqr(i,k) = fr(i,k) - rainfall(i)                          &
             + precip_product_uv(i,k)*eta_rho(i,k)*scales%zcld(i)/2.

        ELSEIF ( eta_rho(i,k) <= 1.0 + SPACING(eta_rho(i,k)) )THEN
          ! Upper portion of cloud layer
          precip_product_uv(i,k) = p0(i)*beta_rain                  &
             *scales%zcld(i)
          fr(i,k) = (((eta_rho(i,k)-beta_rain)/(1. - beta_rain))    &
             *(gamma_rain - 1.) + 1. )*rainfall(i)
          wqr(i,k) = fr(i,k) - rainfall(i)                          &
             + p0(i)*beta_rain*scales%zcld(i)                       &
             * (eta_rho(i,k) - beta_rain/2.)                        &
             *scales%zcld(i)  
        END IF

        !--------------------------
        ! on theta levels
        !--------------------------

        IF ( eta_theta(i,k) <= beta_rain )THEN
          precip_product_th(i,k) = eta_theta(i,k)*scales%zcld(i)    &
             *p0(i)
        ELSEIF ( eta_theta(i,k) <= 1.0 + SPACING(eta_rho(i,k)) )THEN
          precip_product_th(i,k) = p0(i)*beta_rain                  &
             *scales%zcld(i)
        END IF
      END DO
    END DO
    IF (lhook) CALL dr_hook('TCS_PRECIP:CALC_PRECIP',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_precip

END MODULE tcs_precip
