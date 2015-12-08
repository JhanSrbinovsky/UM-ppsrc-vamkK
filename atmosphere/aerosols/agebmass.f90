! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Convert a proportion of fresh smoke to aged smoke
MODULE agebmass_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE agebmass(                                                           &
  ! Arguments IN
  row_length, rows, off_x, off_y,                                              &
  model_levels, timestep, bmass_new,                                           &
  ! Arguments OUT
  delta_agebmass  )

! Purpose:
!   To convert a proportion of the fresh smoke to aged smoke. This
!   conversion takes place as an exponential decay with an e-folding
!   time of 1.0 days.
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

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with intent IN:

INTEGER ::  row_length         ! no. of pts along a row
INTEGER ::  rows               ! no. of rows
INTEGER ::  off_x              ! size of small halo in i
INTEGER ::  off_y              ! size of small halo in j.
INTEGER ::  model_levels       ! no. of model levels
REAL    ::  timestep           ! timestep

! Arguments with intent IN:

!mmr of fresh smoke
REAL  ::  bmass_new(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

! Arguments with intent OUT:
!smoke increment due to ageing
REAL  :: delta_agebmass(row_length, rows, model_levels) 

! Local variables:

INTEGER            ::  i,j,k  ! loop counters
! conversion rate. Equivalent to 1/(6.0 hours expressed in seconds)
REAL, PARAMETER    ::  rate= 4.6296e-5   
REAL               ::  a   ! Local workspace, equal to 1.0-exp(-rate*timestep)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Cycle through all points in the field, calculating the amount
! of smoke converted on each point on this timestep.

IF (lhook) CALL dr_hook('AGEBMASS',zhook_in,zhook_handle)
a = 1.0 - EXP(-rate*timestep)

DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_agebmass(i,j,k) = a * bmass_new(i,j,k)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('AGEBMASS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE agebmass

END MODULE agebmass_mod
