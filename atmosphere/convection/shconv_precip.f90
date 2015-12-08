! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow precipitation
!

      Subroutine shconv_precip( n_sh, nlev, max_cldlev                  &
     &,                     zcld, mb, wstar_up, ql_ad, dz_inv, rho_inv  &
     &,                     eta_half,eta                                &
     &,                     rainfall, wqr_inv, precip_product_inv       &
     &,                     wqr, precip_product_th, precip_product_uv )

!
! Purpose:
!   Calculate the shallow precipitation
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v6 programming standards
!
      USE earth_constants_mod, ONLY: g

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      Integer, intent(in) ::                                            &
     &   n_sh                                                           &
                         ! No. of shallow convection points
     &,  nlev                                                           &
                         ! No. model levels
     &,  max_cldlev      ! Maximum number of cloud levels

      Real, intent(in) ::                                               &
     &  zcld(n_sh)                                                      &
                               ! depth of convective cloud (m)
     &, mb(n_sh)                                                        &
                               ! cloud base mass flux
     &, wstar_up(n_sh)                                                  &
                               ! convective velocity scale, cloud (m/s)
     &, ql_ad(n_sh)                                                     &
                               ! adiabatic liquid water content (kg/kg)
     &, dz_inv(n_sh)                                                    &
                               ! dz across inversion (m)
     &, rho_inv(n_sh)                                                   &
                               ! density at inversion (kg/m3)
     &, eta_half(n_sh,nlev)                                             &
                               ! eta  uv levels
     &, eta(n_sh,nlev)         ! eta
!
! Arguments with intent INOUT:
!
!          NONE

!
! Arguments with intent OUT:
!

      Real, intent(out) ::                                              &
     &  rainfall(n_sh)                                                  &
                                  ! surface rainfall rate
     &, wqr_inv(n_sh)                                                   & 
                                  ! flux of w'qr' at inversion
     &, precip_product_inv(n_sh)                                        &
                                  ! Integral of precip_product over inv
                                  ! *kappa_rain/density
     &, wqr(n_sh,nlev)                                                  & 
                                  !  w'qr'
     &, precip_product_th(n_sh,nlev)                                    &
                                     ! precipitation production rate
     &, precip_product_uv(n_sh,nlev) ! precipitation production rate


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      Real  ::                                                          &
     & fraction(n_sh)                                                   &
                         ! fraction of ensemble with lwc above threshold
     &,ql_up_inv(n_sh)                                                  &
                        ! mean liquid water content of ensemble at base
                        ! of inversion
     &,p0(n_sh)                                                         &
                        ! precip productio p0 value
     &, fr_inv(n_sh)    ! gravitional flux of precipitation at inversion
      Real  ::                                                          &
     & fr(n_sh,nlev) ! gravitional flux of precipitation

!
! Loop counters
!

      Integer :: i,k


      Real, parameter :: ql_t=0.001       ! threshold for precipitation
      Real, parameter :: epsilon_rain=0.5 ! coefficient for surface
                                          ! precipitation calculation.
      Real, parameter :: beta_rain=0.5    ! coefficient for precip
                                          ! production.
      Real, parameter :: gamma_rain=0.5   ! coefficient for rainfall
                                          ! rate at base of inversion.
      Real, parameter :: kappa_rain=1.0   ! coefficient for rainfall
                                          ! rate at base of inversion.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SHCONV_PRECIP',zhook_in,zhook_handle)
      Do k=1,nlev
        Do i = 1,n_sh
          precip_product_th(i,k) = 0.0
          precip_product_uv(i,k) = 0.0
          wqr(i,k) = 0.0
        End do
      End do

! Calculate liquid water content of buoyant updraughts at the base of
! the inversion
      Do i = 1, n_sh
        ql_up_inv(i) = (1.5/g)*((wstar_up(i)/mb(i))**0.5)*              &
     &                 wstar_up(i)*wstar_up(i)/zcld(i)

      End do

! Does precipitation occur ?

      Do i = 1, n_sh

! case of no precipitation

        If (ql_t >= ql_ad(i)) then
          fraction(i) = 0.0
          rainfall(i) = 0.0
          p0(i) = 0.0
          fr_inv(i) = 0.0
          wqr_inv(i) = 0.0
          precip_product_inv(i) = 0.0
        Else
         If (ql_t > ql_up_inv(i)) then
           fraction(i) = (1.0-ql_t/ql_ad(i))*(1.0-ql_t/ql_ad(i))/       &
     &                     (1.0-ql_up_inv(i)/ql_ad(i))
         Else
           fraction(i) = 1.0 - ql_t*ql_t/(ql_up_inv(i)*ql_ad(i))
         End If

! surface precipitation rate

         rainfall(i) = epsilon_rain*mb(i)*ql_up_inv(i)*fraction(i)

! p0
         p0(i) = rainfall(i)/(beta_rain*(1.-0.5*beta_rain)              &
     &                   *zcld(i)*zcld(i) +                             &
     &                          0.5*beta_rain*zcld(i)*dz_inv(i))
! fr at inversion

         fr_inv(i) = gamma_rain*rainfall(i)
! wqr at inversion
!    wqr = rainfall -0.5 p0 beta zcld dz_inv +(gamma-1)rainfall
! simplified to
         wqr_inv(i) = gamma_rain*rainfall(i)                                       &
     &                        - 0.5*p0(i)*beta_rain*zcld(i)*dz_inv(i)   


         precip_product_inv(i) =kappa_rain*0.5*p0(i)*beta_rain          &
     &                               *zcld(i)*dz_inv(i)/rho_inv(i)

        End If
      End Do

!-----------------------------------------------------------------------
! calculate in cloud profiles of precipitation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Rain water flux wqr in cloud integrate precip_product & dfr/dz
!-----------------------------------------------------------------------

      Do k=1,max_cldlev+1
        Do i = 1,n_sh
          If (eta_half(i,k) <= 1.0) then

            If (eta_half(i,k) < beta_rain) then
              precip_product_uv(i,k) = eta_half(i,k)*zcld(i)*p0(i)
              fr(i,k) = rainfall(i)
              wqr(i,k) = fr(i,k) - rainfall(i)+precip_product_uv(i,k)   &
     &                          *0.5*eta_half(i,k)*zcld(i)
            Else
              precip_product_uv(i,k) = p0(i)*beta_rain*zcld(i)
              fr(i,k) = rainfall(i)*                                          &
     &             (1.0 -gamma_rain*beta_rain+eta_half(i,k)*(gamma_rain-1.0)) &
     &                         /(1.0-beta_rain)
              wqr(i,k) = fr(i,k) - rainfall(i)                             &
     &                                +p0(i)*beta_rain*zcld(i)*zcld(i)     &
     &                                      *(eta_half(i,k)-0.5*beta_rain)         
            End If
          Else  ! above  cloud
            precip_product_uv(i,k) = 0.0
            fr(i,k) = 0.0
            wqr(i,k) = 0.0
          End If
        End do
      End do

! theta levels

      Do k=1,max_cldlev+1
        Do i = 1,n_sh
          If (eta(i,k) <= 1.0) then
            If (eta(i,k) < beta_rain) then
              precip_product_th(i,k) = eta(i,k)*zcld(i)*p0(i)
            Else
              precip_product_th(i,k) = p0(i)*beta_rain*zcld(i)
            End If
          Else  ! above  cloud
            precip_product_th(i,k) = 0.0
          End If
        End do
      End do

!-----------------------------------------------------------------------
!    End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_PRECIP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE shconv_precip
