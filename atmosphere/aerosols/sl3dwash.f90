! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
MODULE sl3dwash_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE sl3dwash(                                                           &
  row_length, rows,                                                            &
  off_x, off_y, halo_i, halo_j,                                                &
  model_levels, wet_model_levels,                                              &
  timestep,                                                                    &
  rho_r2,                                                                      &
  q, qcl, qcf, s_mmr,                                                          &
  ls_rain3d,                                                                   &
  lscav_tr                                                                     &
  )

!----------------------------------------------------------------------
! Purpose: Scavenging of SO2 by large-scale ppn.
!
!          Called by microphys_ctl if Sulphur Cycle is on.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation:  UMDP 20
!
!----------------------------------------------------------------------
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
INTEGER :: off_x                ! EW size of std. halo
INTEGER :: off_y                ! NS size of std. halo
INTEGER :: halo_i               ! EW extended halo
INTEGER :: halo_j               ! NS extended halo
INTEGER :: model_levels
INTEGER :: wet_model_levels


! density*r*r
REAL :: rho_r2        (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
                           model_levels)
REAL :: ls_rain3d(row_length,rows,wet_model_levels)        ! rain rate (kg/m2/s)
REAL :: q(row_length,rows,wet_model_levels)                ! water vapour(kg/kg)
REAL :: qcl(row_length,rows,wet_model_levels)
REAL :: qcf(row_length,rows,wet_model_levels)

REAL :: timestep                                           ! timestep in secs

! Arguments with intent IN/OUT:
! mmr S in SO2
REAL :: s_mmr(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

! Arguments with intent OUT (diagnostics):
REAL lscav_tr(row_length,rows)                            ! accumulated scavenged tr

! Local constants

! Factor to convert S mmr to ppbv (mol wt air / mol wt S) x 10**9
REAL,PARAMETER :: mmr_to_ppbv = 0.9033e9
! Threshold value determining form of scavenging function
REAL,PARAMETER :: s_thold     = 0.3065                
REAL,PARAMETER :: l0          = 6.5e-5    !Scavenging parameter for S <= S_THOLD
REAL,PARAMETER :: l1          = 2.955e-5  !Scavenging parameter for S >  S_THOLD

! Local variables
INTEGER ::  i, j, k             ! loop counters

REAL :: termr                   !scavenging rate s-1
REAL :: delta_tr                !tracer increment due to scavnging
REAL :: dm                      !mass p.u.area of air in layer
REAL :: s_ppbv                  !S in SO2 in ppbv
REAL :: rain_mmph               !RAIN rate in mm/hour
! Densities at adjacent model levels
REAL :: rho1
REAL :: rho2

LOGICAL :: l_do_scav            ! for controlling scavenging

IF (lhook) CALL dr_hook('SL3DWASH',zhook_in,zhook_handle)

!   Initialise LSCAV_TR to zero before doing rainout
DO j=1,rows
  DO i=1,row_length
    lscav_tr(i,j) = 0.0
  END DO
END DO


DO j=1,rows
  DO i=1,row_length

    l_do_scav=.TRUE.

    ! Loop over wet levels from surface up, and only scavenge at points
    ! from which ppn reaches the ground

    DO k=1,wet_model_levels

      IF (l_do_scav) THEN

        IF (ls_rain3d(i,j,k)  >   0) THEN

          ! Calculation of scavenging rate requires S in PPBV
          ! and rainrate in mm/hr
          ! Convert S_MMR to S_PPBV
          s_ppbv = s_mmr(i,j,k)*mmr_to_ppbv

          ! Convert LS_RAIN3D to RAIN_MMPH
          rain_mmph = ls_rain3d(i,j,k)*3600.0

          ! Calculate scavenging rate using function dependent on amount of S
          ! present:  d(ln(S))/dt = - L1 * (R/S_THOLD)**2/3  if S  <=  S_THOLD
          !           d(ln(S))/dt = - L1 * (R/S      )**2/3  if S  >   S_THOLD

          IF (rain_mmph <= 0.0) THEN       !Check for neg ppn
            termr = 0.0
          ELSE

            IF (s_ppbv <= s_thold) THEN    !Use fn for low S
              termr = l0 * ( rain_mmph )**0.666667
            ELSE                           !Use fn for high S
              termr = l1 * ( rain_mmph/s_ppbv )**0.666667
            END IF

          END IF

          ! Calculate amount of tracer scavenged out per tstep

          delta_tr=s_mmr(i,j,k)*(1.0-EXP(-termr*timestep))

          ! Calculate mass of air per unit area in layer for conversion of tracer
          ! mixing ratio increment to mass p.u.a. for STASH
          ! Avoid calculations at top model level and assume no scavenging there.

          IF (k  <   model_levels) THEN
            ! Remove the r squared factor from rho_r2 before interpolation
            rho1= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *                        &
              r_rho_levels(i,j,k) )
            rho2= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *                    &
              r_rho_levels(i,j,k+1) )

            ! DM = density (interpolated on to theta levels) * delta r
            dm= rho2 * ( r_theta_levels(i,j,k) -                               &
              r_rho_levels(i,j,k) ) +                                          &
              rho1 * ( r_rho_levels(i,j,k+1) -                                 &
              r_theta_levels(i,j,k) )

            ! Special case for lowest layer to get correct mass
            IF (k  ==  1) THEN
              dm = dm *                                                        &
                (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))                    &
                /(r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
            END IF

            ! Convert DM to DRY density if level is wet
            IF (k  <=  wet_model_levels) THEN
              dm = dm * (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
            END IF

          ELSE
            delta_tr = 0.0
            dm = 0.0
          END IF

          ! Increment accumulated scavenged tracer in column, multiplying by DM
          ! to convert mmr to mass per unit area.

          lscav_tr(i,j)=lscav_tr(i,j)+delta_tr*dm

          ! Decrement tracer mixing ratio

          s_mmr(i,j,k)=s_mmr(i,j,k)-delta_tr

        ELSE                    !level with zero ppn found

          l_do_scav=.FALSE.     !stop scavenging in this column

        END IF                  !END LS_RAIN3D > 0 condition

      END IF                    !END L_DO_SCAV condition

    END DO                      !END K loop

  END DO                        !End i loop
END DO                          !End j loop

IF (lhook) CALL dr_hook('SL3DWASH',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sl3dwash

END MODULE sl3dwash_mod
