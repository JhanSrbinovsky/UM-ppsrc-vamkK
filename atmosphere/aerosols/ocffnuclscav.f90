! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Perform nucleation scavenging of aged ocff to cloud ocff
MODULE ocffnuclscav_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ocffnuclscav(                                                       &
  ! Arguments IN
  rows, row_length, off_x, off_y, halo_i, halo_j,                              &
  model_levels, wet_model_levels, timestep,                                    &
  cloudf, qcl, qcf,                                                            &
  ! Arguments INOUT
  ocff_agd, ocff_cld,                                                          &
  ! Arguments OUT
  delta_ocffnuclscav                                                           &
  )

USE atmos_constants_mod, ONLY: r
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE c_ocff_chm_mod, ONLY: cloud_tau, evap_tau, nuc_tau, thold
IMPLICIT NONE
!---------------------------------------------------------------------
! Description:
!   To perform nccleation scavenging of aged ocff to form
!   the third mode of ocff, ocff in cloud water.
!
!   Called by Aero_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20
!
!---------------------------------------------------------------------

! Global variables (#include statements etc):

!  includes parameters for rate of ocff nucleation/evaporation

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Arguments with intent IN:

INTEGER :: rows                         !no. of rows
INTEGER :: row_length                   !no. of pts along a row
INTEGER :: off_x                        !size of small halo in i
INTEGER :: off_y                        !size of small halo in j
INTEGER :: halo_i                       !EW halo size
INTEGER :: halo_j                       !NS halo size
INTEGER :: model_levels                 !no. of model levels
INTEGER :: wet_model_levels             !no. of wet model levels

REAL :: cloudf(row_length,rows,wet_model_levels)               !Decimal cloud fraction
!Cloud liquid water
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
!Cloud frozen water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,wet_model_levels)
REAL :: timestep                                               !Atmos model timestep

! Arguments with intent IN:

! mmr of aged OCFF
REAL :: ocff_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)
! mmr of OCFF-in-cloud
REAL :: ocff_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)    

! Arguments with intent OUT:

! cloud ocff increment due to nucleation scavenging
REAL :: delta_ocffnuclscav(row_length,rows,model_levels)

! Local variables:

INTEGER :: i,j,k  ! loop counters

REAL :: delta_nuc(row_length,rows,model_levels)    !Increment to cloud OCFF
REAL :: delta_evap(row_length,rows,model_levels)   !Increment to aged OCFF
REAL :: qctotal(row_length,rows,wet_model_levels)  !Total cloud water
REAL :: clear_frac              !Clear fraction of grid box (1.0 - cloudf)
REAL :: evaptime                !timescale for cloud droplets to evaporate
REAL :: nuctime                 !timescale for particles to enter a cloud and nucleate.

IF (lhook) CALL dr_hook('OCFFNUCLSCAV',zhook_in,zhook_handle)

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
! 2. Release ocff from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any ocff-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------
DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
      IF (qctotal(i,j,k) .LT. thold) THEN
        delta_evap(i,j,k) = ocff_cld(i,j,k)
        !                      evaporate all the cloud ocff in this grid box
      ELSE IF (cloudf(i,j,k) .LT. 0.95) THEN
        evaptime=evap_tau + 0.5*cloud_tau
        delta_evap(i,j,k) = (1.0 - EXP(-timestep/evaptime)) *                  &
          ocff_cld(i,j,k)
      ELSE
        delta_evap(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO

!  Also evaporate any ocff-in-cloud on non-wet model levels

IF (wet_model_levels .LT. model_levels) THEN
  DO k=wet_model_levels+1,model_levels
    DO j=1,rows
      DO i=1,row_length
        delta_evap(i,j,k) = ocff_cld(i,j,k)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged ocff particles by cloud droplets.
!    It is assumed that the ocff particles can nucleate nucleat
!    droplets.

!-----------------------------------------------------------------------


DO k=1,wet_model_levels
  DO j=1,rows
    DO i=1,row_length

      clear_frac = 1.0 - cloudf(i,j,k)

      IF ((qctotal(i,j,k) .GE. thold) .AND.                                    &
        (cloudf(i,j,k) .GT.0.0)) THEN


        nuctime=nuc_tau+((clear_frac*cloud_tau)/(2.0*cloudf(i,j,k)))
        delta_nuc(i,j,k)=(1.0-EXP(-timestep/nuctime))*ocff_agd(i,j,k)

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
      delta_ocffnuclscav(i,j,k) =                                              &
        delta_nuc(i,j,k) - delta_evap(i,j,k)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('OCFFNUCLSCAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ocffnuclscav
END MODULE ocffnuclscav_mod
