! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calculate the gradient of the turbulent fluxes
! w'theta' and w'q' on in-cloud levels
!
MODULE tcs_grad_flux


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

  SUBROUTINE calc_grad_flux(n_xx, ntra, nlev, trlev, maxlev, l_tracer &
     ,                      r2rho,r2rho_th,dr_across_th               &
     ,                      wthetavl, wthetal, wqt,wtracer            &
     ,                      dwthetavl_dz, dwthetal_dz, dwqt_dz        &
     ,                      dwtracer_dz)

    IMPLICIT NONE
    !--------------------------------------------------------------------
    ! Subroutine Arguments
    !--------------------------------------------------------------------
    !
    ! Arguments with intent IN:
    !
    INTEGER, INTENT(in) ::                                             &
       n_xx                                                            &
                                ! No. of congestus convection points
       , ntra                                                          &
                                ! No. of tracer
       , maxlev                                                        &
                                ! Maximum number of levels where 
                                ! gradient is non zero
       , nlev                                                          &
                                ! Maximum number of convective cloud levels
       , trlev     
                                ! Maximum number of tracer levels

    LOGICAL, INTENT(in) ::                                             &
       l_tracer                 ! true - tracers present


    REAL, INTENT(in) ::                                                &
       r2rho(n_xx,nlev)                                                &
                                ! radius**2 density on rho lev (kg/m)
       , r2rho_th(n_xx,nlev)                                           &
                                ! radius**2 density on theta lev (kg/m)
       , dr_across_th(n_xx,nlev)  
                                !  thickness on theta levels (m)

    !
    ! fluxes  all held on rho levels
    !
    REAL, INTENT(in) ::                                                &
       wthetavl(n_xx,nlev)                                             & 
                                ! w'thetavl'
       , wthetal(n_xx,nlev)                                            & 
                                ! w'theta'
       , wqt(n_xx,nlev)                                                & 
                                ! w'qt'
       , wtracer(n_xx,trlev,ntra)  
                                ! w'tracer' (kgm/kg/s)
    !
    ! Arguments with intent INOUT:
    !
    !     None
    !
    ! Arguments with intent OUT:
    !
    REAL, INTENT(out) ::                                                &
       dwthetavl_dz(n_xx,nlev)                                          &
                                ! dwthetavl/dz  on theta levels
       ,dwthetal_dz(n_xx,nlev)                                          &
                                ! dwthetal/dz  on theta levels
       ,dwqt_dz(n_xx,nlev)                                              &
                                ! dwqt/dz      on theta levels
       ,dwtracer_dz(n_xx,trlev,ntra)
                                ! dwtracer/dz   on theta levels (kg/kg/s)

    !-----------------------------------------------------------------------
    ! Variables defined locally
    !-----------------------------------------------------------------------

    INTEGER :: max_trlev

    REAL ::                                                            &
       rdz                                                             &
                                ! 1/(dz*r2*rho)
       , r2_kp1                                                        &
                                ! r**2 rho at k+1 levels
       , r2_k       
                                ! r**2 rho at k levels


    !-------------------------
    ! Loop counters
    !-------------------------
    INTEGER :: i,k,ktra  

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !-----------------------------------------------------------------------
    ! 1.0 Initialise arrays - set all functions to zero
    !-----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('TCS_GRAD_FLUX:CALC_GRAD_FLUX',zhook_in,zhook_handle)

    DO k = 1,nlev
      DO i = 1,n_xx
        dwqt_dz(i,k) = 0.0
        dwthetal_dz(i,k) = 0.0
        dwthetavl_dz(i,k) = 0.0
      END DO
    END DO

    !-----------------------------------------------------------------------
    ! 2.0 Level loop to calculate gradients of fluxes
    !-----------------------------------------------------------------------

    DO k = 1,maxlev    ! problem if nlev (no rho theta)
      DO i = 1,n_xx

        rdz  =1./(dr_across_th(i,k)*r2rho_th(i,k))


        r2_kp1 = r2rho(i,k+1)
        r2_k   = r2rho(i,k)

        dwqt_dz(i,k)     = (r2_kp1*wqt(i,k+1)-r2_k*wqt(i,k))*rdz

        dwthetal_dz(i,k) = (r2_kp1*wthetal(i,k+1)                      &
           -r2_k*wthetal(i,k))*rdz

        dwthetavl_dz(i,k)= (r2_kp1*wthetavl(i,k+1)                     &
           -r2_k*wthetavl(i,k))*rdz
      END DO
    END DO


    !-----------------------------------------------------------------------
    ! 3.0 Tracers
    !-----------------------------------------------------------------------


    IF (l_tracer) THEN

      ! Check max levels ?
      ! May be problems if tracers on less levels < maxlev ?
      max_trlev = maxlev
      IF (trlev  <=  maxlev) THEN
        max_trlev = trlev-1
      END IF

      DO ktra = 1,ntra
        DO k = 1,trlev
          DO i = 1,n_xx
            dwtracer_dz(i,k,ktra) = 0.0
          END DO
        END DO

        DO k = 1,max_trlev  ! problem if nlev (no rho theta)
          DO i = 1,n_xx

            rdz  =1./(dr_across_th(i,k)*r2rho_th(i,k))
            r2_kp1 = r2rho(i,k+1)
            r2_k   = r2rho(i,k)
            dwtracer_dz(i,k,ktra) = (r2_kp1*wtracer(i,k+1,ktra)      &
               -r2_k*wtracer(i,k,ktra))*rdz

          END DO
        END DO
      END DO
    END IF
    IF (lhook) CALL dr_hook('TCS_GRAD_FLUX:CALC_GRAD_FLUX',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_grad_flux

END MODULE tcs_grad_flux
