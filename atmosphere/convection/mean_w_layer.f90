! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the average w over a layer
!
MODULE mean_w_layer_mod

IMPLICIT NONE

CONTAINS

! Subroutine Interface:
SUBROUTINE mean_w_layer(nunstable,row_length,rows,model_levels,               &
               k_start, index_i, index_j,                                     &
               depth, z_full_c, z_half_c, w, dmass_theta,                     &
               w_avg)

! ------------------------------------------------------------------------------
! Description:
!   Calculate the average w over a layer
!  
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: &
  nunstable            & ! Number of parcel ascents
 ,row_length           & ! Local number of points on a row
 ,rows                 & ! Local number of rows in a theta field
 ,model_levels           ! Number of model levels


INTEGER, INTENT(IN) :: &
  k_start(nunstable)   & ! Level above which require average
 ,index_i(nunstable)   & ! column number of unstable points
 ,index_j(nunstable)     ! row number of unstable points


REAL, INTENT(IN)    :: &
  depth                  ! depth of layer (m)

REAL, INTENT(IN)    ::                   &
  z_full_c(nunstable, model_levels)      & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)      & ! Height uv lev    (m)
 ,w(row_length,rows,0:model_levels)      & ! vertical velocity (m/s)
 ,dmass_theta(nunstable,model_levels)      ! r**2rho*dr on theta levels

REAL, INTENT(OUT)   :: &
  w_avg(nunstable)       ! Average value of w in required layer (m/s)


! Local variables

INTEGER :: ii, k,i, j            !Loop counters 

REAL ::                &
  mass(nunstable)        ! mass of required layer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------------

IF (lhook) CALL dr_hook('mean_w_layer',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise arrays
DO ii=1,nunstable
  w_avg(ii) = 0.0
  mass(ii)  = 0.0
END DO

DO k=1,model_levels-1
  DO ii=1,nunstable
    i = index_i(ii)   
    j = index_j(ii)   
    IF (k >= k_start(ii).AND.                                              &
        z_full_c(ii,k) <= (z_half_c(ii,k_start(ii)+1)+depth)) THEN

       mass(ii)  = mass(ii) + dmass_theta(ii,k)
       w_avg(ii) = w_avg(ii) + w(i,j,k)*dmass_theta(ii,k)

    END IF
  END DO
END DO

DO ii=1,nunstable
  IF (mass(ii)  >  0.0 ) THEN
    w_avg(ii) = w_avg(ii)/mass(ii)
  END IF
END DO

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('mean_w_layer',zhook_out,zhook_handle)

END SUBROUTINE mean_w_layer

END MODULE mean_w_layer_mod
