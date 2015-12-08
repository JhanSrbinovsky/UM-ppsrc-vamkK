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
MODULE tcs_wql

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  !   Module to calculate the liquid water flux
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

  SUBROUTINE calc_wql_warm(n_xx, nlev, max_cldlev, ncld_thlev, wql_cb, &
     wql_inv, cld_in, t_cld, sim, wql_cld, dqsatdt_cld) 

    USE tcs_constants,          ONLY :                 &
        Lc, Rv, repsilon, g
    USE tcs_class_similarity,   ONLY :                 &
        similarity
    USE tcs_class_cloud,        ONLY :                 &
        cloud_input
    USE tcs_common_warm,        ONLY :                 &
        scales

    IMPLICIT NONE
    !----------------------------------------------------------------
    ! Subroutine Arguments
    !----------------------------------------------------------------
    INTEGER, INTENT(in) ::                                            &
       n_xx                                                           &
                          ! No. of congestus convection points
       , nlev                                                         &
                          ! No. of model layers
       , max_cldlev       
                          ! maximum noumber of in cloud levels

    INTEGER, INTENT(in) ::                                            &
       ncld_thlev(n_xx)   ! number of theta cloud levels

    REAL, INTENT(in)    ::                                            &
       wql_cb(n_xx)                                                   &
                          ! flux of wql across cloud base
       , wql_inv(n_xx)
                          ! flux of wql across inversion


    TYPE(cloud_input), INTENT(in) :: cld_in


    REAL, INTENT(in) :: t_cld(n_xx,nlev) ! temperature in cloud

    TYPE(similarity), INTENT(in) :: sim  ! Similarity functions

    REAL, INTENT(inout) ::                                            &
       wql_cld(n_xx,nlev)     ! wql in cloud (includes cloud base)

    REAL, INTENT(out) ::                                              &
       dqsatdt_cld(n_xx,max_cldlev+1)    ! dqsat/dt

    !----------------------------------------------------------------
    ! Variables defined locally
    !----------------------------------------------------------------
    REAL ::                                                           &
       thetav_cld(n_xx,max_cldlev+1)                                  &
       , temp(n_xx)

    !-------------------------
    ! Loop counters
    !-------------------------
    INTEGER :: i,k,itop

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !----------------------------------------------------------------
    ! 1.0  Calculate dqsat/dT  in cloud
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    !   dqsatdt on cloud levels  (expression for mixing ratio)
    !----------------------------------------------------------------

    IF (lhook) CALL dr_hook('TCS_WQL:CALC_WQL_WARM',zhook_in,zhook_handle)

    DO k = 1,max_cldlev+1
      DO i = 1,n_xx
        dqsatdt_cld(i,k) = cld_in%qse(i,k)*lc                             &
           /(Rv*t_cld(i,k)*t_cld(i,k))
        thetav_cld(i,k)=cld_in%theta(i,k)*(1.0+cld_in%q_mix(i,k)/repsilon)&
           /(1.+cld_in%q_mix(i,k))
      END DO
    END DO

    !----------------------------------------------------------------
    ! 2.0 Calculate  in cloud values of wql on uv levels using
    !      similarity expression
    !----------------------------------------------------------------
    ! values for new water flux parametrisation

    k=1
    DO i=1,n_xx
      wql_cld(i,k) = wql_cb(i)
      temp(i) =scales%root_mb_new_o_wsc(i)*(scales%wstar_up(i)**3)/(scales%zcld(i)*g)
    END DO

    ! taking thetav/theta=1.

    DO k=2,max_cldlev+1
      DO i=1,n_xx
        wql_cld(i,k) = 0.5*temp(i)*sim%fql_func_rho(i,k)* &
           cld_in%theta(i,k)/thetav_cld(i,k)
      END DO
    END DO

    ! value at inversion

    DO i=1,n_xx
      itop=ncld_thlev(i)+1
      wql_cld(i,itop) = wql_inv(i)
    END DO
    IF (lhook) CALL dr_hook('TCS_WQL:CALC_WQL_WARM',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_wql_warm
    
END MODULE tcs_wql
