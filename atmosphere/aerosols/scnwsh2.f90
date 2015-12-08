! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!  Purpose: Scavenge SO2 by convective precipitation
MODULE scnwsh2_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE scnwsh2(                                               &
           timestep,                                              &
           rho,                                                   &
           q, qcl, qcf, s_mmr,                                    &
           ccldbase, ccldtop,                                     &
           rainrate, accu_scav_tr                                 &
           )


! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                  &
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
 ,rainrate(tdims%i_end,tdims%j_end)     !conv rain rate at surface kg/m2

REAL, INTENT(IN) ::        &
  timestep                    !timestep in secs


REAL,INTENT(INOUT) ::                              &
  s_mmr(tdims_s%i_start:tdims_s%i_end,             & !mmr S in SO2
        tdims_s%j_start:tdims_s%j_end,             &
        tdims_s%k_start:tdims_s%k_end)

REAL , INTENT(OUT) ::                      &
  accu_scav_tr(tdims%i_end, tdims%j_end)     !column total of scvnged trcr


!  Local variables

INTEGER :: i, j, k                !loop counters
INTEGER :: start_level            !lowest level for scavenging

REAL ::                      &
  termr                      & ! to assist calcn of scav rate
 ,terms                      & !
 ,dm                         & ! mass p.u.area of air in layer
 ,delta_tr                   & ! tracer increment due to scvnging
 ,rho1, rho2                 & ! air densities
 ,s_ppbv                     & ! S in SO2 in ppbv
 ,rain_mmph                    ! RAIN rate in mm/hour


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


REAL, PARAMETER ::            &
  mmr_to_ppbv = 0.9033e9      & ! 28.966 E9/32.066
                                ! Factor to convert S mmr to ppbv
                                !  (mol wt air / mol wt S) x 10**9
 ,s_thold = 0.3065            & ! Threshold value determining form of
                                ! scavenging function (ppbv)
 ,l0 = 6.5e-5                 & ! Scavenging parameter for S <= S_THOLD (s-1)
 ,l1 = 2.955e-5                 ! Scavenging parameter for S >  S_THOLD (s-1)

IF (lhook) CALL dr_hook('SCNWSH2',zhook_in,zhook_handle)


! Initialise ACCU_SCAV_TR array to zero before adding accumulations
  DO j=1,tdims%j_end
    DO i=1,tdims%i_end
      accu_scav_tr(i,j)=0.0
    END DO
  END DO

! Set level for scavenging to begin
start_level=1

! Calculate scavenging rate if convective cloud present

DO j=1,tdims%j_end
  DO i=1,tdims%i_end

    IF (ccldtop(i,j)  >   0) THEN

    ! Calculate scavenging rate using function dependent on amount of S
    ! present:  d(ln(S))/dt = - L1 * (R/S_THOLD)**2/3  if S  <=  S_THOLD
    !           d(ln(S))/dt = - L1 * (R/S      )**2/3  if S  >   S_THOLD

    ! Calculation of scavenging rate requires S in PPBV
    ! and rainrate in mm/hr

    ! Convert RAINRATE to RAIN_MMPH and increase rate to obtain rate in
    ! cloudy part of grid box assuming CCA=0.05

      rain_mmph = MAX(0.0, rainrate(i,j) * 3600.0 / 0.05)

      DO k=start_level, ccldtop(i,j)

        ! Convert S_MMR to S_PPBV
        s_ppbv = s_mmr(i,j,k) * mmr_to_ppbv

        IF (s_ppbv <= s_thold) THEN        !Use fn for low S
          termr = l0 * ( rain_mmph )**0.666667
        ELSE                               !Use fn for high S
          termr = l1 * ( rain_mmph/s_ppbv )**0.666667
        END IF

        ! Calculate amount of tracer mixing ratio scavenged out per tstep

        delta_tr = s_mmr(i,j,k)*(1.0-EXP(-termr*timestep))

        ! Reduce DELTA_TR to allow for non_cloudy part of grid box
        ! assuming CCA=0.05

        delta_tr = delta_tr * 0.05

        ! Calculate mass of air per unit area in layer for conversion of tracer
        ! mixing ratio increment to mass p.u.a. for STASH
        ! Avoid calculations at top model level and assume no scavenging there.

        IF (k  <   tdims%k_end) THEN
          ! Remove the r squared factor from rho before interpolation
          rho1= rho(i,j,k)  /( r_rho_levels(i,j,k)   * r_rho_levels(i,j,k) )
          rho2= rho(i,j,k+1)/( r_rho_levels(i,j,k+1) * r_rho_levels(i,j,k+1) )

          ! DM = density (interpolated on to theta levels) * delta r
          dm = rho2 * ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) ) +    &
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

        s_mmr(i,j,k)=s_mmr(i,j,k)-delta_tr

      END DO                    !END k loop

    END IF                      !END conv cloud present condn

  END DO                        !End i loop
END DO                          !End j loop

IF (lhook) CALL dr_hook('SCNWSH2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE scnwsh2
END MODULE scnwsh2_mod
