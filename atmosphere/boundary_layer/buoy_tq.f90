! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE BUOY_TQ
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! PURPOSE: To calculate buoyancy parameters on p,T,q-levels
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!
SUBROUTINE buoy_tq (                                                    &
! IN dimensions/logicals
 bl_levels, lq_mix_bl,                                                  &
! IN fields
 p,t,q,qcf,qcl,cf,                                                      &
! OUT fields
 bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                     &
 )

  USE atm_fields_bounds_mod, ONLY: tdims
  USE water_constants_mod, ONLY: lc, lf, tm
  USE earth_constants_mod, ONLY: g
  USE atmos_constants_mod, ONLY: cp, r, repsilon, c_virtual
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

  LOGICAL, INTENT(IN) :: lq_mix_bl  ! IN switch for using mixing ratios
  INTEGER, INTENT(IN) ::                                                &
   bl_levels              ! IN No. of atmospheric levels for which

  REAL, INTENT(IN) ::                                                   &
   p(tdims%i_start:tdims%i_end,                                         &
     tdims%j_start:tdims%j_end,bl_levels),                              &
                                      ! IN Pressure at pressure points.
   t(tdims%i_start:tdims%i_end,                                         &
     tdims%j_start:tdims%j_end,bl_levels),                              &
                                      ! IN Temperature (K). At P points
   q(tdims%i_start:tdims%i_end,                                         &
     tdims%j_start:tdims%j_end,bl_levels),                              &
                                ! IN Sp humidity (kg water per kg air).
   qcl(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                                ! IN Cloud liq water (kg per kg air).
   qcf(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                                ! IN Cloud liq water (kg per kg air).
   cf(tdims%i_start:tdims%i_end,                                        &
      tdims%j_start:tdims%j_end, bl_levels)! IN Cloud fraction (decimal).

! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

  REAL, INTENT(OUT) ::                                                  &
   bq(tdims%i_start:tdims%i_end,                                        &
      tdims%j_start:tdims%j_end,bl_levels),                             &
                                ! OUT A buoyancy parameter for clear air
   bt(tdims%i_start:tdims%i_end,                                        &
      tdims%j_start:tdims%j_end,bl_levels),                             &
                                ! OUT A buoyancy parameter for clear air
   bq_cld(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,bl_levels),                         &
                              ! OUT A buoyancy parameter for cloudy air
   bt_cld(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,bl_levels),                         &
                              ! OUT A buoyancy parameter for cloudy air
   bq_gb(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,bl_levels),                          &
                                ! OUT A grid-box mean buoyancy parameter
   bt_gb(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,bl_levels),                          &
                                ! OUT A grid-box mean buoyancy parameter
   a_qs(tdims%i_start:tdims%i_end,                                      &
         tdims%j_start:tdims%j_end,bl_levels),                          &
                              ! OUT Saturated lapse rate factor
   a_dqsdt(tdims%i_start:tdims%i_end,                                   &
           tdims%j_start:tdims%j_end,bl_levels),                        &
                              ! OUT Saturated lapse rate factor
   dqsdt(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,bl_levels)
                              ! OUT Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.

  REAL ::                                                               &
   qs(tdims%i_start:tdims%i_end,                                        &
      tdims%j_start:tdims%j_end)      ! WORK Saturated mixing ratio.

  INTEGER ::                                                            &
    i,j,                                                                &
    k

  REAL ::                                                               &
    bc

  REAL :: etar,grcp,lcrcp,lfrcp,ls,lsrcp

  PARAMETER (                                                           &
   etar=1.0/(1.0-repsilon),                                             &
                                ! Used in buoyancy parameter BETAC.
   grcp=g/cp,                                                           &
                                ! Used in DZTL, FTL calculations.
   lcrcp=lc/cp,                                                         &
                                ! Latent heat of condensation / CP.
   lfrcp=lf/cp,                                                         &
                                ! Latent heat of fusion / CP.
   ls=lc+lf,                                                            &
                                ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                                ! Latent heat of sublimation / CP.
  )

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BUOY_TQ',zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 1.  Loop round levels.
!-----------------------------------------------------------------------
  DO k = 1, bl_levels
!-----------------------------------------------------------------------
! 1.1 Calculate saturated specific humidity at pressure and
!     temperature of current level.
!-----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs,t(1,1,k),p(1,1,k),tdims%i_end*tdims%j_end,lq_mix_bl)

    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end

!-----------------------------------------------------------------------
! 1.2 Calculate buoyancy parameters BT and BQ, required for the
!     calculation of stability.
!-----------------------------------------------------------------------

        bt(i,j,k) = 1.0/t(i,j,k)
        bq(i,j,k) =                                                     &
          c_virtual/(1.0+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))

        IF (t(i,j,k)  >   tm) THEN
          dqsdt(i,j,k) = (repsilon * lc * qs(i,j))                      &
                       / ( r * t(i,j,k) * t(i,j,k) )
!                      ...  (Clausius-Clapeyron) for T above freezing

          a_qs(i,j,k) = 1.0 / (1.0 + lcrcp*dqsdt(i,j,k))

          a_dqsdt(i,j,k) = a_qs(i,j,k) * dqsdt(i,j,k)

          bc = lcrcp*bt(i,j,k) - etar*bq(i,j,k)

        ELSE
          dqsdt(i,j,k) = (repsilon * ls * qs(i,j))                      &
                       / ( r * t(i,j,k) * t(i,j,k) )
!                      ...  (Clausius-Clapeyron) for T below freezing

          a_qs(i,j,k) = 1.0 / (1.0 + lsrcp*dqsdt(i,j,k))

          a_dqsdt(i,j,k) = a_qs(i,j,k) * dqsdt(i,j,k)

          bc = lsrcp*bt(i,j,k) - etar*bq(i,j,k)

        END IF

!-----------------------------------------------------------------------
! 1.3 Calculate in-cloud buoyancy parameters.
!-----------------------------------------------------------------------

        bt_cld(i,j,k) = bt(i,j,k) - a_dqsdt(i,j,k) * bc
        bq_cld(i,j,k) = bq(i,j,k) + a_qs(i,j,k) * bc

!-----------------------------------------------------------------------
! 1.4 Calculate grid-box mean buoyancy parameters.
!-----------------------------------------------------------------------

        bt_gb(i,j,k) = bt(i,j,k) +                                      &
                       cf(i,j,k)*( bt_cld(i,j,k) - bt(i,j,k) )
        bq_gb(i,j,k) = bq(i,j,k) +                                      &
                       cf(i,j,k)*( bq_cld(i,j,k) - bq(i,j,k) )

      END DO ! p_points,j
    END DO ! p_points,i
  END DO ! bl_levels

  IF (lhook) CALL dr_hook('BUOY_TQ',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE buoy_tq
