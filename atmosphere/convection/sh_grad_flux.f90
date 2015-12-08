! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  New Shallow convection scheme - Calculation of turbulent fluxes
!
      SUBROUTINE sh_grad_flux(n_sh, ntra, nlev, trlev, maxlev, l_tracer &
     &,                      r2rho,r2rho_th,dr_across_th                &
     &,                      wthetav, wthetal, wqt,wtracer              &
     &,                      dwthetav_dz, dwthetal_dz, dwqt_dz          &
     &,                      dwtracer_dz)

!
! Purpose:
!   Calculate gradient of the turbulent fluxes
!   w'theta' and w'q' on in cloud levels
!
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
     &, ntra                                                            &
                  ! No. of tracer
     &, maxlev                                                          &
                  ! Maximum number of levels where gradient is non zero
     &, nlev                                                            &
                  ! Maximum number of convective cloud levels
     &, trlev     ! Maximum number of tracer levels

      logical, intent(in) ::                                            &
     &  l_tracer  ! true - tracers present


      real, intent(in) ::                                               &
     &  r2rho(n_sh,nlev)                                                &
                                 ! radius**2 density on rho lev (kg/m)
     &, r2rho_th(n_sh,nlev)                                             &
                                 ! radius**2 density on theta lev (kg/m)
     &, dr_across_th(n_sh,nlev)  !  thickness on theta levels (m)

!
! fluxes  all held on rho levels
!
      real, intent(in) ::                                               &
     &  wthetav(n_sh,nlev)                                              & 
                                ! w'thetav'
     &, wthetal(n_sh,nlev)                                              & 
                                ! w'theta'
     &, wqt(n_sh,nlev)                                                  & 
                                ! w'qt'
     &, wtracer(n_sh,trlev,ntra)  ! w'tracer' (kgm/kg/s)
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     & dwthetav_dz(n_sh,nlev)                                           &
                               ! dwthetav/dz  on theta levels
     &,dwthetal_dz(n_sh,nlev)                                           &
                               ! dwthetal/dz  on theta levels
     &,dwqt_dz(n_sh,nlev)                                               &
                               ! dwqt/dz      on theta levels
     &,dwtracer_dz(n_sh,trlev,ntra)
                               ! dwtracer/dz   on theta levels (kg/kg/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------


      integer :: i,k,ktra         ! loop counters
      integer :: max_trlev

      real ::                                                           &
     & rdz                                                              &
                 ! 1/(dz*r2*rho)
     &, r2_kp1                                                          &
                   ! r**2 rho at k+1 levels
     &, r2_k       ! r**2 rho at k levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! Model constants:
!


!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------

       IF (lhook) CALL dr_hook('SH_GRAD_FLUX',zhook_in,zhook_handle)

       Do k = 1,nlev
         Do i = 1,n_sh
           dwqt_dz(i,k) = 0.0
           dwthetal_dz(i,k) = 0.0
           dwthetav_dz(i,k) = 0.0
         End do
       End do

!-----------------------------------------------------------------------
! 2.0 Level loop to calculate gradients of fluxes
!-----------------------------------------------------------------------

       Do k = 1,maxlev    ! problem if nlev (no rho theta)
         Do i = 1,n_sh

           rdz  =1./(dr_across_th(i,k)*r2rho_th(i,k))


           r2_kp1 = r2rho(i,k+1)
           r2_k   = r2rho(i,k)

           dwqt_dz(i,k)     = (r2_kp1*wqt(i,k+1)-r2_k*wqt(i,k))*rdz

           dwthetal_dz(i,k) = (r2_kp1*wthetal(i,k+1)                    &
     &                                        -r2_k*wthetal(i,k))*rdz

           dwthetav_dz(i,k)= (r2_kp1*wthetav(i,k+1)                     &
     &                                        -r2_k*wthetav(i,k))*rdz
         End do
       End do

!-----------------------------------------------------------------------
! 3.0 Tracers
!-----------------------------------------------------------------------


       If (l_tracer) then

! Check max levels ?
! May be problems if tracers on less levels < maxlev ?
         max_trlev = maxlev
         If (trlev  <=  maxlev) then
           max_trlev = trlev-1
         End If

         Do ktra = 1,ntra
           Do k = 1,trlev
             Do i = 1,n_sh
               dwtracer_dz(i,k,ktra) = 0.0
             End do
           End do

           Do k = 1,max_trlev  ! problem if nlev (no rho theta)
             Do i = 1,n_sh

               rdz  =1./(dr_across_th(i,k)*r2rho_th(i,k))
               r2_kp1 = r2rho(i,k+1)
               r2_k   = r2rho(i,k)
               dwtracer_dz(i,k,ktra) = (r2_kp1*wtracer(i,k+1,ktra)      &
     &                                  -r2_k*wtracer(i,k,ktra))*rdz

             End do
           End do
         End do
       End If
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SH_GRAD_FLUX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE sh_grad_flux
