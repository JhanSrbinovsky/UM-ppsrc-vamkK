! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
!  Purpose: Scavenge Sulphur Cycle tracers by convective precipitation
!
MODULE scnscv2_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE scnscv2( timestep,                                     &
           rho,                                                   &
           q, qcl, qcf, tracer,                                   &
           ccldbase, ccldtop,                                     &
           rainrate, snowrate,                                    &
           l_scav_below_cloud,                                    &
           k_rain, k_snow,                                        &
           accu_scav_tr                                           &
           )

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   tdims, qdims, tdims_s, pdims_s

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!---------------------------------------------------------------------------
! Description: 
!  Scavenge Sulphur Cycle tracers by convective precipitation
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosol
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

!---------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::                &
  ccldbase(tdims%i_end,tdims%j_end)   & !convective cloud base
 ,ccldtop(tdims%i_end,tdims%j_end)      !convective cloud top

REAL, INTENT(IN) ::                   &
  rho(pdims_s%i_start:pdims_s%i_end,  & !  density *r*r (kg/m)
      pdims_s%j_start:pdims_s%j_end,  & !
      pdims_s%k_start:pdims_s%k_end)  & !
 ,q(qdims%i_start:qdims%i_end,        & ! water vapour(kg/kg)
    qdims%j_start:qdims%j_end,        & !
                1:qdims%k_end)        & !
 ,qcl(qdims%i_start:qdims%i_end,      & ! cloud liquid water (kg/kg)
      qdims%j_start:qdims%j_end,      & !
                  1:qdims%k_end)      & !
 ,qcf(qdims%i_start:qdims%i_end,      & ! cloud ice water (kg/kg)
      qdims%j_start:qdims%j_end,      & !
                  1:qdims%k_end)      & !
 ,rainrate(tdims%i_end,tdims%j_end)   & !conv rain rate at surface kg/m2
 ,snowrate(tdims%i_end,tdims%j_end)     !conv snow rate at surface kg/m2


REAL, INTENT(IN) ::        &
  timestep                 &  !timestep in secs
 ,k_rain                   &  !scavenging rate coeff for rain
 ,k_snow                      !scavenging rate coeff for snow

LOGICAL, INTENT(IN) ::     &
  l_scav_below_cloud          !control for scavenging levels


REAL,INTENT(INOUT) ::                               & ! Tracer (kg/kg)
  tracer(tdims_s%i_start:tdims_s%i_end,             &
         tdims_s%j_start:tdims_s%j_end,             &
         tdims_s%k_start:tdims_s%k_end)


REAL , INTENT(OUT) ::                      &
  accu_scav_tr(tdims%i_end, tdims%j_end)     !column total of scvnged trcr

! Local variables

INTEGER ::                    &
  start_level                 & ! lowest level for scavenging
 ,i, j, k                       ! Loop variables

REAL ::                      &
  termr                      & ! to assist calcn of scav rate
 ,terms                      & !
 ,dm                         & ! mass p.u.area of air in layer
 ,delta_tr                   & ! tracer increment due to scvnging
 ,rho1, rho2                   ! air densities

REAL ::                      &
     totrate                   ! total scav rate

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SCNSCV2',zhook_in,zhook_handle)



! Initialise ACCU_SCAV_TR array to zero before adding accumulations
  DO j=1,tdims%j_end
    DO i=1,tdims%i_end
      accu_scav_tr(i,j)=0.0
    END DO
  END DO

! Calculate total scavenging rate

DO j=1,tdims%j_end
  DO i=1,tdims%i_end

    IF (ccldtop(i,j) >  0) THEN
      
      ! Set up START_LEVEL for scavenging
      IF (l_scav_below_cloud) THEN
        start_level = 1
      ELSE
        start_level = ccldbase(i,j)
      END IF

      IF (rainrate(i,j) <= 0.0) THEN    !check for negative ppn
        termr=0.0
      ELSE
        termr=k_rain*rainrate(i,j)
      END IF

      IF (snowrate(i,j) <= 0.0) THEN
        terms=0.0
      ELSE
        terms=k_snow*snowrate(i,j)
      END IF

      ! Calculate TOTRATE, *3600.0 because K_RAIN and K_SNOW are derived for
      !  ppn rates in mm/hr, but model values are kg/m2/s (cf CON_SCAV)

      totrate=(termr+terms)*3600.0*timestep

      ! Increase TOTRATE to obtain rate in cloudy part of grid box
      ! Assume CCA=0.05

      totrate=totrate / 0.05

      ! Calculate amount of tracer scavenged and add to column total

      DO k = start_level, ccldtop(i,j)

        ! Calculate proportion of tracer mixing ratio scavenged out
        delta_tr=tracer(i,j,k)*(1.0-EXP(-totrate))

        ! Reduce DELTA_TR to allow for non_cloudy part of grid box
        delta_tr = delta_tr * 0.05

        ! Calculate mass of air per unit area in layer for conversion of tracer
        !  mixing ratio increment to mass p.u.a. for STASH.
        ! Avoid calculations at top model level and assume no scavenging there.

        IF (k  <   tdims%k_end) THEN
          ! Remove the r squared factor from rho before interpolation
          rho1=rho(i,j,k)  /( r_rho_levels(i,j,k)   * r_rho_levels(i,j,k) )
          rho2=rho(i,j,k+1)/( r_rho_levels(i,j,k+1) * r_rho_levels(i,j,k+1) )

          ! DM = density (interpolated on to theta levels) * delta r
          dm= rho2 * ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) ) +    &
              rho1 * ( r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k) )

          ! Special case for lowest layer to get correct mass
          IF (k  ==  1) THEN
             dm= dm * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))   &
                    / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
          END IF

          ! Convert DM to DRY density if level is wet
          IF (k  <=  qdims%k_end) THEN
            dm = dm * (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
          END IF

        ELSE
          delta_tr = 0.0
          dm = 0.0
        END IF

        ! Increment column total mass p.u.a. of scavenged tracer

        accu_scav_tr(i,j)=accu_scav_tr(i,j)+delta_tr*dm

        ! Decrement tracer mixing ratio

        tracer(i,j,k)=tracer(i,j,k)-delta_tr

      END DO                     !end k loop

    END IF                       !END CCLDTOP > 0 TEST

  END DO                         !End i loop
END DO                           !End j loop

IF (lhook) CALL dr_hook('SCNSCV2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE scnscv2
END MODULE scnscv2_mod
