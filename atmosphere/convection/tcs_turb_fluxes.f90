! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to calclulate the turbulent fluxes w'theta' and w'q' 
! on in-cloud levels
!
MODULE tcs_turb_fluxes


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  !
  ! Description:
  ! Module to calclulate the turbulent fluxes w'theta' and w'q' 
  ! on in-cloud levels
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

  SUBROUTINE calc_turb_fluxes(n_xx,ntra, max_cldlev, nlev              &
     ,                      max_cldtrlev, trlev                        &
     ,                      l_tracer                                   &
     ,                      kterm_thv, kterm_h, kterm_tracer           &
     ,                      sim                                        &
     ,                      wql_cld, mf_cld, cld_in                    &
     ,                      dqsatdt_cld                                &
     ,                      precip_product_cld,wup_h_cld               &
     ,                      scales,wthetav_cb,wh_cb,wtracer_cb         &
     ,                      wtheta, wq, wthetavl, wthetav              &
     ,                      wthetal, wqt, wh, wtracer)

    USE tcs_parameters_warm,    ONLY:                                  &
       icong_options, wthvl_factor
    USE tcs_constants,          ONLY:                                  &
        lc_o_cp, c_virtual, g
    USE tcs_class_similarity,   ONLY:                                  &
        similarity
    USE tcs_class_scales,       ONLY:                                  &
        scales_conv
    USE tcs_class_cloud,        ONLY:                                  &
        cloud_input

    IMPLICIT NONE
    !------------------------------------------------------------------
    ! Subroutine Arguments
    !------------------------------------------------------------------
    !
    ! Arguments with intent IN:
    !
    INTEGER, INTENT(in) ::                                             &
       n_xx                                                            &
                                ! No. of congestus convection points
       , ntra                                                          &
                                ! No. of tracers
       , max_cldlev                                                    &
                                ! Maximum number of convective cloud levels
       , nlev                                                          &
                                ! Maximum number of convective cloud levels
       , max_cldtrlev                                                  &
                                ! Maximum number of convective cloud levels
       , trlev                                                           
    ! Maximum number of tracer levels

    LOGICAL, INTENT(in) ::                                             &
       l_tracer      ! true - tracers present

    REAL, INTENT(in) ::                                                &
       kterm_thv(n_xx,nlev)                                            &
                                ! thetav gradient (kterm)
       ,  kterm_h(n_xx,nlev)
    ! moist static energy(h) gradient
    ! (kterm)

    TYPE(similarity) :: sim   ! similarity functions

    REAL, INTENT(in) :: &
       wql_cld(n_xx,nlev)                                              &
                                ! wql on uv cloud lev (similarity func)
       , mf_cld(n_xx,nlev) 
    ! mass_flux on theta levels in cloud

    TYPE(cloud_input), INTENT(in) :: cld_in
    ! input fields on cloud levels

    REAL, INTENT(in) ::                                                &
       dqsatdt_cld(n_xx,nlev)                                          &
                                !dqsat/dT on theta levels in cloud
       , precip_product_cld(n_xx,nlev)                                 &
                                ! precip production rho levels
       , wup_h_cld (n_xx,nlev)                                         &
                                ! wup on rho levels
       , kterm_tracer(n_xx,nlev,ntra)
    ! tracer gradient (kterm)
    !                                   dimensioned on nlev not trlev

    TYPE(scales_conv), INTENT(in) :: scales

    REAL, INTENT(in) ::                                                &
       wthetav_cb(n_xx)                                                &
                                ! w'thetav' at cloud base (K m/s)
       , wh_cb(n_xx)                                                   &
                                ! w'h' at cloud base (K m/s)
       , wtracer_cb(n_xx,ntra)   ! w'tracer' at cloud base (kg/kg m/s)
    !
    ! Arguments with intent INOUT:
    !
    !     None
    !
    ! Arguments with intent OUT:
    !
    REAL, INTENT(out) ::                                               &
       wtheta(n_xx,nlev)                                               &
                                ! wtheta flux on cloud levels  (K m/s)
       ,wq(n_xx,nlev)                                                  &
                                ! wq flux on cloud levels      (kg/kg m/s)
       ,wthetavl(n_xx,nlev)                                            &
                                ! wthetavl flux on cloud levels  (K m/s)
       ,wthetav(n_xx,nlev)                                             &
                                ! wthetav flux on cloud levels  (K m/s)
       ,wthetal(n_xx,nlev)                                             &
                                ! wthetal flux on cloud levels  (K m/s)
       ,wqt(n_xx,nlev)                                                 &
                                ! wqt flux on cloud levels      (kg/kg m/s)
       ,wh(n_xx,nlev)                                                  &
                                ! wh flux on cloud levels      (K m/s)
       ,wtracer(n_xx,trlev,ntra) ! wtracer flux on cloud levels(kg/kg m/s)

    !-----------------------------------------------------------------------
    ! Variables defined locally
    !-----------------------------------------------------------------------
    REAL ::                                                           &
       wtheta_vl_cb(n_xx)                                             
    ! wtheta_vl flux at cloud base.

    ! temporary stores
    REAL ::                                                           &
       fngterm2,  A_term                                              &
       , div_term

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
    ! Initialise dhdz

    IF (lhook) CALL dr_hook('TCS_TURB_FLUXES:CALC_TURB_FLUXES',zhook_in,zhook_handle)

    DO k = 1,max_cldlev
      DO i = 1,n_xx
        wq(i,k) = 0.0
        wqt(i,k) = 0.0
        wtheta(i,k) = 0.0
        wthetavl(i,k) = 0.0
        wthetav(i,k) = 0.0
        wthetal(i,k) = 0.0
      END DO
    END DO


    !-----------------------------------------------------------------------
    ! 2.0 Level loop to calculate in cloud values of fluxes on uv levels
    !-----------------------------------------------------------------------
    !  Numbering of levels
    !  level 1 is cloud base for uv - in cloud levels
    !
    !-----------------------------------------------------------------------
    DO i = 1,n_xx

      wtheta_vl_cb(i) = wthetav_cb(i)

    END DO


    DO k = 1,max_cldlev
      DO i = 1,n_xx

        ! thetav flux
        fngterm2 = wtheta_vl_cb(i)*sim%fng_func_rho(i,k)


        ! additional term for precip added to original w'thetav'
        wthetavl(i,k) = -kterm_thv(i,k) + fngterm2                          &
           + (lc_o_cp/cld_in%exner_theta(i,k)-c_virtual*cld_in%theta(i,k))  &
           *scales%zcld(i)*(wup_h_cld(i,k)/scales%wstar_up(i))              &
           *precip_product_cld(i,k)/cld_in%rho(i,k)

      END DO
    END DO
    ! Here we've split up the loop so that wthetavl can be renormalized if 
    ! required


    IF  (icong_options == 7) THEN 
      ! Rescale wthetavl using TKE scaling
      DO i = 1,n_xx
        wthetavl(i,1:max_cldlev)=wthetavl(i,1:max_cldlev)                   &
           *max_cldlev/SUM(wthetavl(i,1:max_cldlev))                        &
           *SIGN(1.,wthetavl(i,1:max_cldlev))                               &
           *wthvl_factor*(cld_in%theta(i,1:max_cldlev)/g)                   &
           *(scales%mb(i)/scales%wstar_up(i))**(0.5)                        &
           *scales%wstar_up(i)**3/scales%zcld(i)
      END DO
    END IF

    DO k = 1,max_cldlev
      DO i = 1,n_xx            

        ! calculations depend on interpolation for q and qsat
        ! preferred method
        IF(sim%gql_func_rho(i,k) /= 0.0) THEN
          A_term = (1.+lc_o_cp*dqsatdt_cld(i,k))*wql_cld(i,k)               &
             /sim%gql_func_rho(i,k) - mf_cld(i,k)                           &
             *(cld_in%q_mix(i,k) - cld_in%qse(i,k))
        ELSE  ! not in cloud
          a_term=0.0
        END IF

        ! Should I use mf on th or uv levels?

        div_term=1.+c_virtual*cld_in%exner_theta(i,k)*cld_in%theta(i,k)     &
           *dqsatdt_cld(i,k)

        wthetav(i,k)=wthetavl(i,k)                                          &
           +(lc_o_cp/cld_in%exner_theta(i,k)-c_virtual*cld_in%theta(i,k))   &
           *wql_cld(i,k)


        wthetal(i,k)=(wthetavl(i,k)-c_virtual*cld_in%theta(i,k)*a_term)     &
           /div_term

        wqt(i,k) =( A_term + wthetavl(i,k)*cld_in%exner_theta(i,k)          &
           *dqsatdt_cld(i,k))/div_term


        wqt(i,k)=MAX(0.,wqt(i,k))

!
! Currenly no using moist static energy anywhere, 
! so no need to calulate this
!
!        wh(i,k) = -kterm_h(i,k) + wh_cb(i)*sim%fng_func_rho(i,k)
      END DO
    END DO


    ! ----------------------------------------------------------------------
    ! 3.0 Tracers
    ! ----------------------------------------------------------------------

    IF (l_tracer) THEN

      ! initialise
      DO ktra = 1,ntra
        DO k = 1,max_cldtrlev
          DO i = 1,n_xx
            wtracer(i,k,ktra) = 0.0
          END DO
        END DO
      END DO

      DO ktra = 1,ntra
        DO k = 1,max_cldtrlev
          DO i = 1,n_xx

            ! non-gradient term
            fngterm2 = wtracer_cb(i,ktra)*sim%fng_func_rho(i,k)

            wtracer(i,k,ktra) = -kterm_tracer(i,k,ktra) + fngterm2

          END DO
        END DO
      END DO
    END IF
    IF (lhook) CALL dr_hook('TCS_TURB_FLUXES:CALC_TURB_FLUXES',zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE calc_turb_fluxes

END MODULE tcs_turb_fluxes
