! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate the fluxes through cloud base.
!
MODULE tcs_cloudbase

  ! Working note: need to sort out level and jump calculations done
  ! here with those done in the main routine


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! This module defines the tcs warm rain "cloud_input" derived type and 
  ! subroutines for allocating and deallocating an instance of this class.
  !
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

  SUBROUTINE calc_cloudbase( n_xx, ntra, l_tracer,                    &
       wthetav_surf, dtracer_cb,                                      &
       dthetav_cb,dthetal_cb, drt_cb,                                 &
       wtheta_plus_cb, wthetal_cb,wq_plus_cb,                         &
       wqt_cb,wql_cb, wh_cb,                                          &
       wthetav_cb, wthetav_cb_m2, wtracer_cb)
!-----------------------------------------------------------------------
!
! Description:
!   Calculate the fluxes at cloud base
!-----------------------------------------------------------------------
  
    USE tcs_parameters_warm,      ONLY :     &
       wql_cb_a1, jump_cb_a1, wthv_cb_a1, sat_a1
    USE tcs_constants,            ONLY :     &
       lc_o_cp, c_virtual, lc, rv
    USE tcs_common_warm ,         ONLY :     &
       scales, cb, cb_m1
    
    IMPLICIT NONE
    !------------------------------------------------------------------
    ! Subroutine Arguments
    !------------------------------------------------------------------
    !
    ! Arguments with intent IN:
    !

      INTEGER, INTENT(in) ::                                           &
     &   n_xx                                                          &
                         ! No. of congestus convection points
     &,  ntra            ! No. of tracers

      LOGICAL, INTENT(in) ::                                           &
     &   l_tracer        ! true for tracers

      REAL, INTENT(in) ::                                              &
     &  wthetav_surf(n_xx)                                             & 
                               ! w'thetav' at surface
     &, dtracer_cb(n_xx,ntra)                                          &
                               ! change in tracers across cloud base
                               ! (kg/kg)
     &, dthetav_cb(n_xx)                                                
                               ! dthetav across cloud base
      !
      ! Arguments with intent INOUT:
      !
      !          NONE
      
      !
      ! Arguments with intent OUT:
      !
      
      REAL, INTENT(out) ::                                              &
     &  wthetal_cb(n_xx)                                                & 
                                  !  w'thetal' at cloud base
     &, wqt_cb(n_xx)                                                    & 
                                  !  w'qt'  at cloud base
     &, wql_cb(n_xx)                                                    & 
                                  !  w'ql' at cloud base
     &, wh_cb(n_xx)                                                     & 
                                  !  w'h' at cloud base
     &, wq_plus_cb(n_xx)                                                & 
                                  !  w'q' at cloud base
     &, wtheta_plus_cb(n_xx)                                            & 
                                  !  w'theta' at cloud base
     &, wthetav_cb(n_xx)                                                & 
                                  !  w'thetav' at cloud base
     &, wthetav_cb_m2(n_xx)                                             & 
                                  !  w'thetav' at cloud base using
                                  !  method 2
     &, dthetal_cb(n_xx)                                                &
                               ! dthetal across cloud base
     &, drt_cb(n_xx)                                                    &
                               ! drt across cloud base
     &, wtracer_cb(n_xx,ntra)  ! w'tracer' at cloud base

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      REAL                                                              &
     &  T_plus(n_xx)                                                    &
                               ! T at + side of cloud base
     &, T_minus(n_xx)                                                   
                               ! T at - side of cloud base

! temporary variables
      REAL                                                              &
     &  term_a(n_xx)                                                    &
     &, term_c(n_xx)                                                    &
     &, term_d(n_xx)                                                    &
     &, r_rsat_term(n_xx)                                               &
                         ! r-rsat   term
     &, drsatdt_topt(n_xx)      
                         ! drsat/dT at the top of the transition region

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      !-------------------------
      ! Loop counters
      !-------------------------
      INTEGER :: i,ktra  

      !-----------------------------------------------------------------
      ! Calculate Some derived fields
      !-----------------------------------------------------------------
      IF (lhook) CALL dr_hook('TCS_CLOUDBASE:CALC_CLOUDBASE',zhook_in,zhook_handle)

      t_plus(:)  = cb%theta(:)*cb%exner_theta(:)
      t_minus(:) = cb_m1%theta(:)*cb_m1%exner_theta(:)

      
      !-----------------------------------------------------------------
      ! The convective parametrization scheme assumes that the
      ! enviromental air has qcl=0 and qcf=0. i.e thetal=theta and rt=r
      ! (this may not be the case in the model).
      !-----------------------------------------------------------------
      dthetal_cb(:) = cb%theta(:) - cb_m1%theta(:)
      
      drt_cb(:)     = cb%q_mix(:) - cb_m1%q_mix(:)

      !-----------------------------------------------------------------
      ! assume w'ql'+ is value at cloud base
      !-----------------------------------------------------------------

      wql_cb(:) = wql_cb_a1*wthetav_surf(:)                            &
           /(lc_o_cp/cb%exner_theta(:)-(1.+c_virtual)*cb%theta(:))

      !-----------------------------------------------------------------
      ! Note that this is the value of wthetav on the lower
      ! interface of cloudbase and should be negative.  Here the
      ! value should be the same as that for wthetavl since there is
      ! no condensed water.
      !-----------------------------------------------------------------
      wthetav_cb(:) = -jump_cb_a1*scales%mb(:)*dthetav_cb(:) ! = wthv_minus
      
      wthetav_cb_m2(:) = -wthv_cb_a1*wthetav_surf(:)  ! = wthv_minus

      !-----------------------------------------------------------------
      ! expression for dqsat/dT
      !-----------------------------------------------------------------
      drsatdt_topt(:) = Lc *cb%qse(:) /(Rv * t_plus(:)*t_plus(:))

      r_rsat_term(:) = cb%q_mix(:) - cb%qse(:)

      !-----------------------------------------------------------------
      ! Invert saturation condition to get wqt and wthetal
      !-----------------------------------------------------------------
      term_a(:) = wthetav_cb(:)
      
      term_c(:) = wql_cb(:)*(1.+lc_o_cp*drsatdt_topt)/sat_a1
      
      term_d(:) = scales%mb(:)*r_rsat_term
      
      wqt_cb(:) = (term_c(:) - term_d(:)                               &
           + cb%exner_theta(:)*drsatdt_topt(:)*term_a(:) )/            &
           (1.+c_virtual*cb%exner_theta(:)*cb_m1%theta(:)*drsatdt_topt)

      wthetal_cb(:) = term_a(:) - c_virtual*cb_m1%theta(:)*wqt_cb(:)

      wq_plus_cb(:) = wqt_cb(:) - wql_cb(:)

      wtheta_plus_cb(:) = wthetal_cb(:)                                &
                           + lc_o_cp*wql_cb(:)/cb_m1%exner_rho(:)


!      ! Working note: The moist static energy flux through cloudbase
!      ! needs some further investigation from LEM results.
!      ! Working note: not used in this version - fqt_temp was 
!      ! passed through the tcs_temp_module
!      wh_cb(:) = wthetav_cb_m2(:) +                                    &
!           wh_cb_a1*(lc_o_cp-c_virtual*cb_m1%theta(:))*fqt_temp

      ! Note expect wthetal & wtheta to be negative at cloud base and wqt
      ! wq and wql to be positive.
      
      !-----------------------------------------------------------------------
      ! tracers - assumed to behave as other conserved variables
      
      IF (l_tracer) THEN
        DO ktra = 1,ntra
          wtracer_cb(:,ktra) = -jump_cb_a1*scales%mb(:)*dtracer_cb(:,ktra)
        END DO
      END IF
      IF (lhook) CALL dr_hook('TCS_CLOUDBASE:CALC_CLOUDBASE',zhook_out,zhook_handle)
      RETURN
      
    END SUBROUTINE calc_cloudbase
    
  END MODULE tcs_cloudbase
