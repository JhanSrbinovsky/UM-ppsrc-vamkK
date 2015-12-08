! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Perform nucleation scavenging of aged smoke to cloud smoke.
MODULE bmassnuclscav_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE bmassnuclscav(                                                      &
  ! Arguments IN
  rows, row_length, off_x, off_y, halo_i, halo_j,                              &
  model_levels, wet_model_levels, timestep,                                    &
  cloudf, qcl, qcf,                                                            &
  ! Arguments IN
  bmass_agd, bmass_cld,                                                        &
  ! Arguments OUT
  delta_bmassnuclscav                                                          &
  )

! Purpose:
!  To perform nucleation scavenging of aged biomass smoke to form
!   the third mode of smoke, smoke in cloud water.
!
!   Called by Aero_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE c_bm_chm_mod, ONLY: cloudtau, evaptau, nuctau, thold
IMPLICIT NONE

!  includes parameters for rate of smoke nucleation/evaporation
!  includes gas constant, R

! Arguments with intent IN:

INTEGER :: rows                 ! no. of rows
INTEGER :: row_length           ! no. of pts along a row
INTEGER :: off_x                ! size of small halo in i
INTEGER :: off_y                ! size of small halo in j
INTEGER :: halo_i               ! EW halo size
INTEGER :: halo_j               ! NS halo size
INTEGER :: model_levels         ! no. of model levels
INTEGER :: wet_model_levels     ! no. of wet model levels

REAL :: cloudf(row_length,rows,wet_model_levels)        ! Decimal cloud fraction
! Cloud liquid water
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
! Cloud frozen water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
! Atmos model timestep
REAL :: timestep                                        

! Arguments with intent IN:
REAL ::bmass_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)
                                !mmr of aged smoke
! mmr of smoke-in-cloud
REAL ::bmass_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)    

! Arguments with intent OUT:

! cloud smoke increment due to nucleation scavenging
REAL :: delta_bmassnuclscav(row_length,rows,model_levels)

! Local variables:

INTEGER :: i,j,k  ! loop counters

REAL :: delta_nuc(row_length,rows,model_levels)   ! Increment to cloud smoke
REAL :: delta_evap(row_length,rows,model_levels)  ! Increment to aged smoke
REAL :: qctotal(row_length,rows,wet_model_levels) ! Total cloud water
REAL :: clear_frac                                ! Clear fraction of grid box (1.0 - cloudf)

REAL :: evaptime         ! timescale for cloud droplets to evaporate
REAL :: nuctime          ! timescale for particles to enter a cloud and nucleate.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('BMASSNUCLSCAV',zhook_in,zhook_handle)


!-----------------------------------------------------------------------
! 1. Initialise increments to zero
!-----------------------------------------------------------------------
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_nuc(i,j,k)=0.0
      delta_evap(i,j,k)=0.0
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
! 2. Release smoke from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any smoke-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
      IF (qctotal(i,j,k)  <   thold) THEN
        delta_evap(i,j,k) = bmass_cld(i,j,k)
        !                      evaporate all the cloud smoke in this grid box
      ELSE IF (cloudf(i,j,k)  <   0.95) THEN
        evaptime=evaptau + 0.5*cloudtau
        delta_evap(i,j,k) = (1.0 - EXP(-timestep/evaptime)) *                  &
          bmass_cld(i,j,k)
      ELSE
        delta_evap(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO

!  Also evaporate any smoke-in-cloud on non-wet model levels

IF (wet_model_levels  <   model_levels) THEN
  DO k=wet_model_levels+1,model_levels
    DO j=1,rows
      DO i=1,row_length
        delta_evap(i,j,k) = bmass_cld(i,j,k)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged smoke particles by cloud droplets.
!    It is assumed that aged smoke particles can nucleate cloud
!    droplets.
!-----------------------------------------------------------------------


DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length

      clear_frac = 1.0 - cloudf(i,j,k)

      IF ((qctotal(i,j,k)  >=  thold) .AND.                                    &
        (cloudf(i,j,k)  >  0.0)) THEN

        nuctime=nuctau+((clear_frac*cloudtau)/(2.0*cloudf(i,j,k)))
        delta_nuc(i,j,k)=(1.0-EXP(-timestep/nuctime))*bmass_agd(i,j,k)

      END IF ! Test on qctotal and thold
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
! 4. Calculate total increment for output
!-----------------------------------------------------------------------

DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_bmassnuclscav(i,j,k) =                                             &
        delta_nuc(i,j,k) - delta_evap(i,j,k)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('BMASSNUCLSCAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE bmassnuclscav
END MODULE bmassnuclscav_mod
