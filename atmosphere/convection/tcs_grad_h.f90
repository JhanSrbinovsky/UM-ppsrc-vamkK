! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calculate dX/dz at time t+1 where X could be h or thetav 
! or some other field. Required as part of flux calculations eg w'qt' 
! and w'thetal'
!
MODULE tcs_grad_h


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! Module to calculate the gradient of the turbulent fluxes
  ! w'theta' and w'q' on in-cloud levels
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

  SUBROUTINE calc_grad_h(n_xx, max_cldlev, nlev, timestep              &
     ,                      nclev                                      &
     ,                      z_rho, z_theta                             &
     ,                      r2rho, r2rho_theta                         &
     ,                      x_t, k_f                                   &
     ,                      wx_inv                                     &
     ,                      kdXdz )

    IMPLICIT NONE
    !------------------------------------------------------------------
    ! Subroutine Arguments
    !------------------------------------------------------------------
    !
    ! Arguments with intent IN:
    !
    INTEGER, INTENT(in) ::                                            &
       n_xx                                                           &
                         ! No. of congestus convection points
       , max_cldlev                                                   &
                         ! Maximum number of convective cloud levels
       , nlev                                                         &
                         ! Maximum number of convective cloud levels
       , nclev(n_xx)   
                         ! cloud levels for point

    REAL, INTENT(in) ::                                               &
       timestep                 ! timestep for convection
    REAL, INTENT(in) ::                                               &
       x_t(n_xx,nlev)                                                 &
                                ! field X at time t
       , z_rho(n_xx,nlev)                                             &
                                ! height of rho levels above surface
       , z_theta(n_xx,nlev)                                           &
                                ! height of theta levels
       , r2rho(n_xx,nlev)                                             &
                                ! r2rho  on rho levels
       , r2rho_theta(n_xx,nlev)                                       &
                                ! r2rho on theta levels
       , k_f(n_xx,nlev)                                               &
                                ! K on rho levels ?
       , wx_inv(n_xx)           
                                ! wx at inversion
    !
    ! Arguments with intent INOUT:
    !
    !     None
    !
    ! Arguments with intent OUT:
    !
    REAL, INTENT(out) ::                                              &
       kdXdz(n_xx,nlev)  ! gradient of field X at t+1

    !-----------------------------------------------------------------------
    ! Variables defined locally
    !-----------------------------------------------------------------------
    REAL ::                                                           &
       a(n_xx,max_cldlev)                                             &
       , b(n_xx,max_cldlev)                                           &
       , c(n_xx,max_cldlev)                                           &
       , x_tp1(n_xx,max_cldlev)                                       &
                                ! X at T+1
       , h_t(n_xx,max_cldlev)                                         &
       , dz_lower, dz_upper, dz_rho

    !-------------------------
    ! Loop counters
    !-------------------------
    INTEGER :: i,k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !-----------------------------------------------------------------------
    ! 1.0 Initialise arrays
    !-----------------------------------------------------------------------
    ! Initialise kdXdz output array

    IF (lhook) CALL dr_hook('TCS_GRAD_H:CALC_GRAD_H',zhook_in,zhook_handle)

    DO k = 1,nlev
      DO i = 1,n_xx
        kdXdz(i,k) = 0.0
      END DO
    END DO

    !-----------------------------------------------------------------------
    ! 2.0 Implict calculation of X and dX/dz
    !-----------------------------------------------------------------------
    ! Note -assume fluxes w'X' are zero at cloud base and inversion
    ! This assumption implies values of H at levels above and below
    ! cloud base are not required
    !-----------------------------------------------------------------------
    ! Problem w'X' not zero at inversion
    ! At present assuming -KdX/dz part is zero - may not be true
    !
    ! w'X'(inv) = -K(inv)[x_t(above)-x_t(below)]/[z(above)-z(below)]
    !                          +w'X'cb FNG(inv)
    !-----------------------------------------------------------------------

    k=1
    DO i=1,n_xx
      dz_rho =( z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

      dz_upper = z_theta(i,k+1) - z_theta(i,k)

      C(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                              &
         *timestep/(dz_rho*dz_upper)

      A(i,k) = 0.0            ! w'h' flux=0 at cloud base
      B(i,k) = 1. -A(i,k) - C(i,k)
      h_t(i,k) = x_t(i,k)
    END DO

    DO k=2,max_cldlev
      DO i=1,n_xx
        IF ((nclev(i)-1) >= k) THEN
          dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

          dz_upper = z_theta(i,k+1) - z_theta(i,k)
          dz_lower = z_theta(i,k) - z_theta(i,k-1)

          C(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                          &
             *timestep/(dz_rho*dz_upper)
          A(i,k) = -k_f(i,k)*r2rho(i,k)                              &
             *timestep/(dz_rho*dz_lower)
          B(i,k) = 1.0 -A(i,k) - C(i,k)

          h_t(i,k) = x_t(i,k)          !
        ELSE IF (nclev(i) == k) THEN

          dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
          dz_lower = z_theta(i,k) - z_theta(i,k-1)
          C(i,k) = 0.0
          A(i,k) = -k_f(i,k)*r2rho(i,k)                              &
             *timestep/(dz_rho*dz_lower)
          B(i,k) = 1. -A(i,k) - C(i,k)
          h_t(i,k) = x_t(i,k)
        ELSE
          ! elements not required in calculation (zero)
          C(i,k) = 0.0
          A(i,k) = 0.0
          B(i,k) = 0.0
          h_t(i,k) = 0.0

        END IF

      END DO
    END DO

    ! CALCULATE NEW TIMESTEP X (moist static energy) using tridiag
    !
    ! DEPENDS ON: tridiag_all
    CALL TRIDIAG_all(max_cldlev,n_xx,nclev,A,B,C,h_t,x_tp1)


    !
    ! Take account of flux at inversion
    !

    DO i=1,n_xx
      k=nclev(i)
      dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
      x_tp1(i,k)=x_tp1(i,k)                                          &
         - wx_inv(i)*r2rho(i,k+1)*timestep/dz_rho

    END DO

    !
    ! Calculate gradient of X from t+1 values
    ! Only applies to cloud interior.
    !
    !  kdXdz= KdX/dz   ie Kterm
    ! Holding values from cloud base to last in cloud level
    ! (ie not inversion)

    DO i = 1,n_xx
      kdXdz(i,1) = 0.0     ! cloud base k_f=0.0
    END DO


    DO k = 2,max_cldlev
      DO i = 1,n_xx
        IF (k <= nclev(i)) THEN
          dz_lower = z_theta(i,k) - z_theta(i,k-1)
          kdXdz(i,k) =k_f(i,k)*(x_tp1(i,k) - x_tp1(i,k-1))           &
             /dz_lower
        END IF
      END DO
    END DO
    IF (lhook) CALL dr_hook('TCS_GRAD_H:CALC_GRAD_H',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_grad_h

END MODULE tcs_grad_h
