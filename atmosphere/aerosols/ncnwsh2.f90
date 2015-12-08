! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
!  Purpose: Scavenge NH3 by convective precipitation
!
MODULE ncnwsh2_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ncnwsh2(                                               &
           timestep,                                              &
           rho,                                                   &
           q, qcl, qcf, n_mmr,                                    &
           ccldbase, ccldtop,                                     &
           rainrate, accu_scav_nh3                                &
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

! ---------------------------------------------------------------------
! Description:
!    Scavenge NH3 by convective precipitation
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: aerosol
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ---------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::                                                            &
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


REAL, INTENT(IN) ::                   &
  timestep                              !timestep in secs

REAL,INTENT(INOUT) ::                              &
  n_mmr(tdims_s%i_start:tdims_s%i_end,             & !mmr N in NH3
        tdims_s%j_start:tdims_s%j_end,             &
        tdims_s%k_start:tdims_s%k_end)


REAL , INTENT(OUT) ::                      &
  accu_scav_nh3(tdims%i_end, tdims%j_end)    !column total of scvnged NH3


!  Local variables

INTEGER :: i, j, k                !loop counters
INTEGER :: start_level            !lowest level for scavenging

REAL ::                      &
  termr                      & !scavenging rate s-1
 ,dm                         & !mass p.u.area of air in layer
 ,delta_tr                   & !tracer increment due to scvng
 ,rho1, rho2                 & !air densities
 ,rain_mmph                    ! RAIN rate in mm/hour

REAL, PARAMETER:: l0        = 1.0e-4 !Scavenging parameter (s-1)
REAL, PARAMETER:: cca_const = 0.05   !Assumed conv. cloud amount

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NCNWSH2',zhook_in,zhook_handle)



! Initialise ACCU_SCAV_NH3 array to zero before adding accumulations
  DO j=1,tdims%j_end
    DO i=1,tdims%i_end
      accu_scav_nh3(i,j)=0.0
    END DO
  END DO

! Set level for scavenging to begin
start_level=1

! Calculate scavenging rate if convective cloud present

DO j=1,tdims%j_end
  DO i=1,tdims%i_end

    IF (ccldtop(i,j)  >   0) THEN

      ! Calculate scavenging rate using formula: d(ln(N))/dt = - L0 * R**2/3

      ! Calculation of scavenging rate requires rainrate Rin mm/hr.
      ! Convert RAINRATE to RAIN_MMPH and increase rate to obtain rate in
      ! cloudy part of grid box assuming CCA=0.05

      rain_mmph = MAX(0.0, rainrate(i,j) * 3600.0 / cca_const)

      DO k=start_level, ccldtop(i,j)

        termr = l0 * ( rain_mmph )**0.666667

        ! Calculate amount of tracer mixing ratio scavenged out per tstep

        delta_tr = n_mmr(i,j,k)*(1.0-EXP(-termr*timestep))

        ! Reduce DELTA_TR to allow for non_cloudy part of grid box
        ! assuming CCA=0.05

        delta_tr = delta_tr * cca_const

        ! Calculate mass of air per unit area in layer for conversion of tracer
        ! mixing ratio increment to mass p.u.a. for STASH
        ! Avoid calculations at top model level and assume no scavenging there.

        IF (k  <  tdims%k_end) THEN
          ! Remove the r squared factor from rho before interpolation
          rho1= rho(i,j,k)  /( r_rho_levels(i,j,k)   * r_rho_levels(i,j,k) )
          rho2= rho(i,j,k+1)/( r_rho_levels(i,j,k+1) * r_rho_levels(i,j,k+1) )

          ! DM = density (interpolated on to theta levels) * delta r
          dm = rho2 * ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) ) +   &
               rho1 * ( r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k) )

          ! Special case for lowest layer to get correct mass
          IF (k  ==  1) THEN
            dm= dm * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))     &
                   / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
          END IF

          ! Convert DM to DRY density if level is wet
          IF (k  <=   qdims%k_end) THEN
            dm = dm * (1.0 - q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))
          END IF

        ELSE
          delta_tr = 0.0
          dm = 0.0
        END IF

        ! Increment column total mass p.u.a. of scavenged tracer

        accu_scav_nh3(i,j)=accu_scav_nh3(i,j)+delta_tr*dm

        ! Decrement tracer mixing ratio

        n_mmr(i,j,k)=n_mmr(i,j,k)-delta_tr

      END DO                    !END k loop

    END IF                      !END conv cloud present condn

  END DO                        !End i loop
END DO                          !End j loop

IF (lhook) CALL dr_hook('NCNWSH2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ncnwsh2

END MODULE ncnwsh2_mod
