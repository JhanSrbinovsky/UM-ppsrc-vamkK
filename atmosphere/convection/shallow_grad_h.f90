! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  New Shallow convection scheme - Calculation of dX/dz at T+1
!
      SUBROUTINE shallow_grad_h(n_sh, max_cldlev, nlev, timestep        &
     &,                      nclev                                      &
     &,                      z_rho, z_theta                             &
     &,                      r2rho, r2rho_theta                         &
     &,                      x_t, k_f                                   &
     &,                      wx_inv                                     &
     &,                      kdXdz )

!
! Purpose:
!   Calculate dX/dz at time t+1 where X could be h or thetav or some
!   other field. Required as part of flux calculations eg w'qt' and
!   w'thetal'
!
!   Called by SHALLOW_CONV5A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      integer, intent(in) ::                                            &
     &  n_sh                                                            &
                      ! No. of shallow convection points
     &, max_cldlev                                                      &
                      ! Maximum number of convective cloud levels
     &, nlev                                                            &
                      ! Maximum number of convective cloud levels
     &, nclev(n_sh)   ! cloud levels for point

      real, intent(in) ::                                               &
     &   timestep         ! timestep for convection
      real, intent(in) ::                                               &
     &  x_t(n_sh,nlev)                                                  &
                               ! field X at time t
     &, z_rho(n_sh,nlev)                                                &
                               ! height of rho levels above surface
     &, z_theta(n_sh,nlev)                                              &
                               ! height of theta levels
     &, r2rho(n_sh,nlev)                                                &
                               ! r2rho  on rho levels
     &, r2rho_theta(n_sh,nlev)                                          &
                               ! r2rho on theta levels
     &, k_f(n_sh,nlev)                                                  &
                               ! K on rho levels ?
     &, wx_inv(n_sh)           ! wx at inversion
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     & kdXdz(n_sh,nlev)  ! gradient of field X at t+1

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------


      integer :: i,k             ! loop counters

      real ::                                                           &
     &  a(n_sh,max_cldlev)                                              &
     &, b(n_sh,max_cldlev)                                              &
     &, c(n_sh,max_cldlev)                                              &
     &, x_tp1(n_sh,max_cldlev)                                          &
                                   ! X at T+1
     &, h_t(n_sh,max_cldlev)                                            &
     &, dz_lower, dz_upper, dz_rho

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
! External routines called:
!
      EXTERNAL                                                          &
     & TRIDIAG_all

!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------
! Initialise kdXdz output array

       IF (lhook) CALL dr_hook('SHALLOW_GRAD_H',zhook_in,zhook_handle)
       Do k = 1,nlev
         Do i = 1,n_sh
            kdXdz(i,k) = 0.0
         End do
       End do

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
       Do i=1,n_sh
           dz_rho =( z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

           dz_upper = z_theta(i,k+1) - z_theta(i,k)

           C(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                            &
     &                                *timestep/(dz_rho*dz_upper)

           A(i,k) = 0.0            ! w'h' flux=0 at cloud base
           B(i,k) = 1. -A(i,k) - C(i,k)
           h_t(i,k) = x_t(i,k)
       End do

       Do k=2,max_cldlev
         Do i=1,n_sh
           If ((nclev(i)-1) >= k) Then
             dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)

             dz_upper = z_theta(i,k+1) - z_theta(i,k)
             dz_lower = z_theta(i,k) - z_theta(i,k-1)

             C(i,k) = -k_f(i,k+1)*r2rho(i,k+1)                          &
     &                                     *timestep/(dz_rho*dz_upper)
             A(i,k) = -k_f(i,k)*r2rho(i,k)                              &
     &                                    *timestep/(dz_rho*dz_lower)
             B(i,k) = 1.0 -A(i,k) - C(i,k)

             h_t(i,k) = x_t(i,k)          !
           Else if (nclev(i) == k) Then

             dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
             dz_lower = z_theta(i,k) - z_theta(i,k-1)
             C(i,k) = 0.0
             A(i,k) = -k_f(i,k)*r2rho(i,k)                              &
     &                                   *timestep/(dz_rho*dz_lower)
             B(i,k) = 1. -A(i,k) - C(i,k)
             h_t(i,k) = x_t(i,k)
           Else
! elements not required in calculation (zero)
             C(i,k) = 0.0
             A(i,k) = 0.0
             B(i,k) = 0.0
             h_t(i,k) = 0.0

           Endif

         End do
       End do

! CALCULATE NEW TIMESTEP X (moist static energy) using tridiag
!
! DEPENDS ON: tridiag_all
       CALL TRIDIAG_all(max_cldlev,n_sh,nclev,A,B,C,h_t,x_tp1)


!
! Take account of flux at inversion
!

       Do i=1,n_sh
         k=nclev(i)
         dz_rho =(z_rho(i,k+1)-z_rho(i,k))*r2rho_theta(i,k)
         x_tp1(i,k)=x_tp1(i,k)                                          &
     &                    - wx_inv(i)*r2rho(i,k+1)*timestep/dz_rho

       End do

!
! Calculate gradient of X from t+1 values
! Only applies to cloud interior.
!
!
!  kdXdz= KdX/dz   ie Kterm
! Holding values from cloud base to last in cloud level
! (ie not inversion)

       Do i = 1,n_sh
         kdXdz(i,1) = 0.0     ! cloud base k_f=0.0
       End do


       Do k = 2,max_cldlev
         Do i = 1,n_sh
           if (k <= nclev(i)) then
             dz_lower = z_theta(i,k) - z_theta(i,k-1)
             kdXdz(i,k) =k_f(i,k)*(x_tp1(i,k) - x_tp1(i,k-1))           &
     &                    /dz_lower
           End if
         End do
       End do


      IF (lhook) CALL dr_hook('SHALLOW_GRAD_H',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE shallow_grad_h
