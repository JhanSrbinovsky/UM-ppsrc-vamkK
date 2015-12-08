! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
MODULE rainout_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE rainout(                                                            &
  row_length, rows,                                                            &
  off_x, off_y, halo_i, halo_j,                                                &
  model_levels, wet_model_levels,                                              &
  rho_r2, q,                                                                   &
  qcf_remain, qcl_remain,                                                      &
  qcf_previous, qcl_previous,                                                  &
  ls_rain3d, ls_snow3d,                                                        &
  timestep,                                                                    &
  aero_incloud,                                                                &
  aero_accum,                                                                  &
  rnout_aero)

!---------------------------------------------------------------------
! Purpose: Makes the rain-out of dissolved aerosols during large-scale
!          precipitations. The rained-out mixing ratio depends on the
!          conversion rate of condensed water to precipitated water.
!          Re-evaporation is taken into account by transfering some
!          of the dissolved aerosol mass to the accumulation mode mass.
!          This is computed as being proportional to the amount of
!          precipitated water that re-evaporates. The routine also
!          accounts for those droplets that shrink without re-evap
!          completely by reevaporating only half of the bulk sulphate
!          concentration (unless rain is completely re-evaporated).
!
!          Called by microphys_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20
!
!---------------------------------------------------------------------

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

! Arguments with intent IN:

INTEGER :: row_length
INTEGER :: rows
INTEGER :: off_x                !EW size of std. halo
INTEGER :: off_y                !NS size of std. halo
INTEGER :: halo_i               !EW extended halo
INTEGER :: halo_j               !NS extended halo
INTEGER :: model_levels
INTEGER :: wet_model_levels

REAL:: rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels) !density*r*r

REAL :: qcl_remain(row_length,rows,wet_model_levels)
REAL :: qcf_remain(row_length,rows,wet_model_levels)
REAL :: qcf_previous(row_length,rows,wet_model_levels)
REAL :: qcl_previous(row_length,rows,wet_model_levels)
REAL :: q(row_length,rows,wet_model_levels)

! the LS_xxxx3D arrays are the large-scale precipitation fluxes.
! The flux is defined at the bottom of the layer. It is the flux that
! falls OUT of the layer on which it is defined. Equivalently, it is the
! flux that falls INTO the layer below the one in which it is defined.
REAL :: ls_rain3d(row_length,rows,wet_model_levels)
REAL :: ls_snow3d(row_length,rows,wet_model_levels)

REAL :: timestep ! timestep in seconds

! Arguments with intent IN/OUT:
!   mass mixing ratio of in-cloud and accumulation/aged aerosol
REAL :: aero_incloud(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
REAL :: aero_accum(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)


! Arguments with intent OUT (diagnostics):
REAL :: rnout_aero(row_length,rows)      !tracer removed kg/m2/ts

!  Local variables
INTEGER ::    i, j, k                    !loop variables

! fraction of droplets that shrinks without re-evaporating
REAL,PARAMETER :: fraction = 0.5

! Dissolved aerosol is completely dissolved
REAL,PARAMETER :: frac_aero_scav = 1.0

REAL :: delta_aero              ! amount of aerosol removed from grid box
REAL :: dm_dry                  ! mass p.u.area of dry air in layer
REAL :: rhodr                   ! density multiplied by layer thick.
REAL :: tot_precip3d            ! total ppn in grid box
REAL :: rho1                    ! densities at adjacent model levels
REAL :: rho2
REAL :: beta                    ! see comments below
REAL :: scav_exp                ! Exponential term in scavanged MMR
REAL :: smallp                  ! small +ve number, negligible compared to 1
REAL :: log_smallp

IF (lhook) CALL dr_hook('RAINOUT',zhook_in,zhook_handle)

! Initialisation
! ... epsilon() is defined as almost negligible, so eps/100 is negligible
smallp = EPSILON(1.0) / 100.0
log_smallp = LOG( smallp) 

DO j = 1, rows
  DO i = 1, row_length
    rnout_aero(i, j) = 0.0
  END DO
END DO

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k,    &
!$OMP& rho1, rho2, rhodr, dm_dry, beta, delta_aero, scav_exp)
DO j = 1, rows
  DO i = 1, row_length

    ! Loop on vertical levels, from the top of the atmosphere
    ! to the surface. We assume no scavenging at the topmost
    ! level. Note that wet_model_levels <= model_levels (= in
    ! most cases)

    DO k = wet_model_levels - 1, 1, -1

      ! compute the mass of air per unit area for conversion
      ! of mass mixing ratio to mass

      ! first, get the density from the array rho * r^2
      rho1 = rho_r2(i,j,k)/ ( r_rho_levels(i,j,k) *                            &
        r_rho_levels(i,j,k) )
      rho2 = rho_r2(i,j,k+1)  / ( r_rho_levels(i,j,k+1)   *                    &
        r_rho_levels(i,j,k+1)   )

      ! RHODR is the density (interpolated onto theta levels) multiplied
      ! by delta_r
      rhodr = rho2 * ( r_theta_levels(i,j,k) -                                 &
        r_rho_levels(i,j,k)   )                                                &
        +                                                                      &
        rho1 * ( r_rho_levels(i,j,k+1) -                                       &
        r_theta_levels(i,j,k) )

      ! special case for the lowest layer
      IF (k  ==  1) THEN
        rhodr = rhodr *                                                        &
          (r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1))                      &
          /(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
      END IF

      ! conversion to dry air density
      dm_dry = rhodr * (1.0 - q(i,j,k) -                                       &
        qcl_previous(i,j,k) - qcf_previous(i,j,k))


      ! compute the conversion rate of cloud water into precipitated water
      beta = ls_rain3d(i,j,k)   + ls_snow3d(i,j,k)   -                         &
        ls_rain3d(i,j,k+1) - ls_snow3d(i,j,k+1)

      IF (qcl_previous(i,j,k)+qcf_previous(i,j,k)  > smallp) THEN
        beta = beta / rhodr / (qcl_previous(i,j,k)+                            &
          qcf_previous(i,j,k))
        beta = MAX( 0.0, beta)

        ! get the scavenged mass mixing ratio
        ! here, delta_aero is <= 0
        ! This exponential can underflow if FRAC_AERO_SCAV is too large
        ! Only calculate if result is significant relative to 1.0
        IF ( frac_aero_scav * beta * timestep < - log_smallp ) THEN
          scav_exp = EXP( -frac_aero_scav * beta * timestep)
        ELSE
          scav_exp = 0.0
        END IF
        delta_aero = aero_incloud(i,j,k) * ( scav_exp - 1.0 )

        aero_incloud(i,j,k) = aero_incloud(i,j,k) + delta_aero

        ! the cumulated rain-out for the current level
        ! (by convention, it is positive)
        rnout_aero(i,j) = rnout_aero(i,j) - delta_aero * dm_dry
      END IF

      ! re-evaporation now.
      ! compute the fraction of precipitated water that re-evaporates

      beta = ls_rain3d(i,j,k)   + ls_snow3d(i,j,k) -                           &
        ls_rain3d(i,j,k+1) - ls_snow3d(i,j,k+1)
      IF (beta  <   0.0) THEN
        beta = beta / ( ls_rain3d(i,j,k+1) + ls_snow3d(i,j,k+1) )
      END IF

      ! in the following, if beta > 0, it is set to 0

      IF ( ls_rain3d(i,j,k)   + ls_snow3d(i,j,k)  ==  0.0) THEN
        ! total re-evaporation
        beta = MIN( MAX(0.0, -beta), 1.0)
      ELSE
        ! non total re-evaporation
        ! fraction is here to take into account those
        ! raindrops that shrinks without evaporating totally
        ! (the value of fraction is somehow arbitrary)
        beta = MIN( MAX(0.0, -beta)*fraction, 1.0)
      END IF

      ! fraction of mass mixing ratio that re-evaporates
      ! into the accumulation/aged mode
      ! DELTA_AERO is >= 0
      delta_aero = beta * rnout_aero(i,j) / dm_dry

      aero_accum(i,j,k) = aero_accum(i,j,k) + delta_aero

      rnout_aero(i,j) = (1.0-beta) * rnout_aero(i,j)

    END DO ! k

  END DO ! i
END DO ! j
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('RAINOUT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE rainout

END MODULE rainout_mod
