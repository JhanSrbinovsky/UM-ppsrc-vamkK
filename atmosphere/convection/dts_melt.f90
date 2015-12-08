! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! calculates the melting rate
!
! Subroutine Interface:
SUBROUTINE dts_melt(nlev,n_dp,dts_ntpar,ntparmax,                     &
                    rho_theta,dr_across_th,temperature,               &
                    depint,subint,freezint,qcfint,                    &
                    melt,meltint)

! Modules used

USE water_constants_mod, ONLY: tm

USE dts_fitpars_mod, ONLY:                                            &
  tsig,  taumelt

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!  Calculates the melting rate
!  
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::  &
  nlev                  & ! No. of model layers
 ,n_dp                  & ! No. convecting points
 ,dts_ntpar(n_dp)       & ! Top level of initial parcel ascent
 ,ntparmax                ! Max ntpar across all conv points
      
REAL, INTENT(IN) ::         &
  rho_theta(n_dp,nlev)      & ! Density on rho levels (kg/m3)
 ,temperature(n_dp,nlev)    & ! temperature in K
 ,dr_across_th(n_dp,nlev)   & ! Thickness of theta levels (m)
 ,depint(n_dp)              & ! Vertical integral of deposition
 ,subint(n_dp)              & ! Vertical integral of sublimation
 ,freezint(n_dp)            & ! Vertical integral of freezing
 ,qcfint(n_dp)                ! Vertical integral of qcf


REAL, INTENT(OUT) ::   &
  melt(n_dp,nlev)      &  ! melting rate
 ,meltint(n_dp)           ! vertical integral of melting rate

! ------------------------------------------------------------------------------
! Local variables
INTEGER :: &
  i_dp,k      ! Loop counters
        
REAL ::         &
  tval          & ! temperature
 ,rtsig2        & ! 1/(tsig*tsig)
 ,frac          & ! 
 ,factor            

REAL ::              &
  icefallint(n_dp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!==============================================================================

  IF (lhook) CALL dr_hook('DTS_MELT',zhook_in,zhook_handle)
  melt(:,:) = 0.0
  meltint(:) = 0.0       
  icefallint(:) = 0.0
  rtsig2 = 1.0/(tsig*tsig)

  icefallint(:) = depint(:)-subint(:)-qcfint(:)-freezint(:)

  DO k=1,ntparmax 
    DO i_dp=1,n_dp

      IF(k <= dts_ntpar(i_dp) .AND. k < nlev) THEN
        tval = temperature(i_dp,k)
                
        IF(tval > tm) THEN
                    
          melt(i_dp,k)  = exp(-(0.5*(tval-tm)**2)*rtsig2) / taumelt
          meltint(i_dp) = meltint(i_dp) + melt(i_dp,k)*rho_theta(i_dp,k)   &
                                                   *dr_across_th(i_dp,k)

        END IF
      END IF
    END DO
  END DO

     
! renormalise to ensure that only melt as much ice as was produced (and fell)
       
  DO i_dp=1,n_dp
    IF(meltint(i_dp) > 0.0) THEN
      factor = icefallint(i_dp)/meltint(i_dp) 
      DO k=1,ntparmax
        melt(i_dp,k) = melt(i_dp,k)*factor
      END DO
      meltint(i_dp) = icefallint(i_dp)
    END IF

! NOT USED at present
! Check to see if surface temperature is warm enough for all to have melted
!   tval = tm+3.*tsig     !nb this condition is a bit arbitrary...
!
!   IF(temperature(i_dp,1) < tval) THEN
!      frac = (temperature(i_dp,1)-tm)/(tval-tm) 
!nb strictly this should be the ratio of the meltint and the integral over
! all T > 273.14K, but doing something a bit arbitrary for the moment
!
!     DO k=1,ntparmax
!nb     melt(i_dp,k) = melt(i_dp,k)*frac
!     END DO
! nb  meltint(i_dp) = meltint(i_dp)*frac
!    END IF

  END DO       ! loop over points 
  IF (lhook) CALL dr_hook('DTS_MELT',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_melt

