! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate CMT increments to U and V.
!
MODULE tcs_cmt_incr_6a


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! This module calculates the cmt increments to U and V winds for 
  ! the tcs warm rain convection schem
  !
  ! Method:
  !   <reference to documentation to go here, once available>
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.2 programming standards.


CONTAINS

  SUBROUTINE calc_cmt_incr(n_xx, nlevs,ntpar_max                       &
     ,                            r2rho, r2rho_th                      &
     ,                            dr_across_rh                         &
     ,                            UW,VW                                &
                                ! output arguements
     ,                            dubydt,dvbydt)

    IMPLICIT NONE
    !-------------------------------------------------------------------
    ! Subroutine Arguments
    !-------------------------------------------------------------------
    !
    ! Arguments with intent IN:
    !
    INTEGER, INTENT(in) ::                                             &
       n_xx                                                            &
       ! total number of congestus convective points
       , nlevs                                                         &
       ! number of model levels
       , ntpar_max      
       ! maximum cloud top level 
    
    REAL, INTENT(in) ::                                                &
       r2rho(n_xx,nlevs)                                               &
                                ! r2 rho on rho levels
       , r2rho_th(n_xx,nlevs)                                          &
                                ! r2 rho on theta levels
       , dr_across_rh(n_xx,nlevs)                                      &
                                ! dr across rho levels
       , UW(n_xx,nlevs)                                                &
                                ! U-Component of stress profile (M2)
       , VW(n_xx,nlevs)     
    ! V-Component of stress profile (M2)
    !
    ! Arguments with intent OUT:
    !
    REAL, INTENT(out) ::                                               &
       dubydt(n_xx,nlevs+1)                                            &
                                ! Tendency in U (MS-2)
       , dvbydt(n_xx,nlevs+1)      
                                ! Tendency in V (MS-2)

    !-----------------------------------------------------------------------
    ! Variables defined locally
    !-----------------------------------------------------------------------
    !
    REAL ::                                                            &
       rhodz                    ! r2 rho * dz

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !-------------------------
    ! Loop counters
    !-------------------------
    INTEGER :: i,k 

    !
    !-----------------------------------------------------------------------
    ! Input to this routine contains stress profile of uw & vw on levels
    ! from the surface to the top of the inversion.
    !-----------------------------------------------------------------------
    !
    ! Lowest uv level surface stress zero
    !
    IF (lhook) CALL dr_hook('TCS_CMT_INCR_6A:CALC_CMT_INCR',zhook_in,zhook_handle)

    k=1
    DO i=1,n_xx
      rhodz  = r2rho(i,k)*dr_across_rh(i,k)
      dubydt(i,K) = -uw(i,k)*r2rho_th(i,k)/rhodz
      dvbydt(i,K) = -vw(i,k)*r2rho_th(i,k)/rhodz
    END DO

    !
    ! Mixed layer and Cloud layer increments
    !
    DO K=2,ntpar_max+1    
      DO i=1,n_xx

        rhodz  = r2rho(i,k)*dr_across_rh(i,k)
        dubydt(i,K) =                                                 &
           -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,K-1))/rhodz
        dvbydt(i,K) =                                                 &
           -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,K-1))/rhodz
      END DO
    END DO
    IF (lhook) CALL dr_hook('TCS_CMT_INCR_6A:CALC_CMT_INCR',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_cmt_incr

END MODULE tcs_cmt_incr_6a
