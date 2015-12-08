! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE MASS_CALC
!
! Purpose:
!   To calculate mass of air in a model layer
!
!   Called by mcr_ctl
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20
!----------------------------------------------------------------------
MODULE mass_calc_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE mass_calc(row_length, rows, model_levels, wet_model_levels,         &
  r_rho_levels, r_theta_levels, timestep, rho_r2, q, qcl, qcf, dm ) 

! calculates mass of (dry) air per square metre

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with intent IN:

INTEGER :: row_length
INTEGER :: rows            
INTEGER :: model_levels    
INTEGER :: wet_model_levels
INTEGER :: ccldbase(row_length,rows)      !convective cloud base
INTEGER :: ccldtop(row_length,rows)       !convective cloud top

REAL :: r_rho_levels(row_length,rows,model_levels)
REAL :: r_theta_levels(row_length,rows,0:model_levels)
REAL :: rho_r2(row_length,rows,model_levels)              ! density*r*r
REAL :: q(row_length,rows,wet_model_levels)               ! water vapour(kg/kg)
REAL :: qcl(row_length,rows,wet_model_levels)
REAL :: qcf(row_length,rows,wet_model_levels)

REAL ::  timestep                                         ! timestep in secs

! Arguments with intent OUT :
REAL :: dm(row_length, rows, model_levels)                ! mass of air

! Local variables

INTEGER :: i, j, k                                        ! Loop variables
INTEGER :: toplev                                         ! level loop limit

REAL    :: rho1, rho2                                     ! air densities

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('MASS_CALC',zhook_in,zhook_handle)
DO k = 1,model_levels-1
  DO j = 1,rows
    DO i = 1,row_length
      ! Remove the r squared factor from rho before interpolation
      rho1 = rho_r2(i,j,k) /                                                   &
        ( r_rho_levels(i,j,k) * r_rho_levels(i,j,k) )
      rho2 = rho_r2(i,j,k+1)/                                                  &
        ( r_rho_levels(i,j,k+1) *  r_rho_levels(i,j,k+1) )
      ! DM = density (interpolated on to theta levels) * delta r
      dm(i,j,k) = rho2 * ( r_theta_levels(i,j,k) -                             &
        r_rho_levels(i,j,k) ) +                                                &
        rho1 * ( r_rho_levels(i,j,k+1) -                                       &
        r_theta_levels(i,j,k) )

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !MODEL_LEVELS

! Special case for lowest layer to get correct mass
DO j = 1,rows
  DO i = 1,row_length
    dm(i,j,1) =                                                                &
      dm(i,j,1) * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0)) /              &
      (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
  END DO !ROW_LENGTH
END DO !ROWS

! Convert DM to DRY density if level is wet
toplev=model_levels - 1
IF (wet_model_levels  <   toplev) toplev = wet_model_levels
DO k = 1,toplev
  DO j = 1,rows
    DO i = 1,row_length
      dm(i,j,k) = dm (i,j,k) /                                                 &
        (1.0 + q(i,j,k) + qcl(i,j,k) + qcf(i,j,k))
    END DO !ROW_LENGTH
  END DO !ROWS
END DO !WET_MODEL_LEVELS

! Top level
DO j = 1,rows
  DO i = 1,row_length
    dm(i,j,model_levels) = 0.0
  END DO !ROW_LENGTH
END DO !ROWS



IF (lhook) CALL dr_hook('MASS_CALC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE mass_calc
END MODULE mass_calc_mod
