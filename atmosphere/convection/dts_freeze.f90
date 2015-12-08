! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the freezing rate
!

SUBROUTINE dts_freeze(nlev,n_dp,dts_ntpar,ntparmax                      &
                     ,z_theta,zlcl,klcl,temperature                     &
                     ,qrainmax,rho_theta,z_rho,freeze,freezint)


USE water_constants_mod, ONLY: tm
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! calculates the melting rate  - SET TO ZERO AT PRESENT
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) ::  &
  nlev                  & ! No. of model layers
 ,n_dp                  & ! No. convecting points
 ,dts_ntpar(n_dp)       & ! Top level of initial parcel ascent
 ,ntparmax              & ! Max ntpar across all conv points
 ,klcl(n_dp)              ! Model level corresponding to LCL
      
REAL, INTENT(IN)    ::    &
  temperature(n_dp,nlev)  & ! temperature in K
 ,z_theta(n_dp,nlev)      & ! Height of theta levels (m)
 ,zlcl(n_dp)              & ! Height of lifting condesnation level (m)
 ,qrainmax(n_dp)          & ! Scaling param for max rain amount
 ,rho_theta(n_dp,nlev)    & ! Density on theta levels (kg/m3)
 ,z_rho(n_dp,nlev)          ! Height of rho levels (m)  

REAL, INTENT(OUT) ::      &
  freeze(n_dp,nlev)       & ! freezing rate     
 ,freezint(n_dp)            ! Column integral of freezing rate

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k             ! Loop counters
        
REAL ::   &
  tval    &
 ,taumelt &
 ,tsig    &
 ,m2 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!====================================================

  IF (lhook) CALL dr_hook('DTS_FREEZE',zhook_in,zhook_handle)

  freezint(:) = 0.0
  freeze(:,:) = 0.0

  taumelt = 400.0 ! s 
  tsig = 7 ! K 


! THIS APPEARS NOT TO DO ANYTHING ? RAS 9/12/09

  DO k=1,ntparmax 
    DO i_dp=1,n_dp
      IF(z_theta(i_dp,k) >= zlcl(i_dp) .and. k <= dts_ntpar(i_dp)) THEN
        tval = temperature(i_dp,k)
        m2 = exp(-((tval-tm)**2)/2./tsig/tsig) / taumelt

!                 IF(tval < tm) THEN
!NB!_think_ this term should probably be zero because deposition term
! is large
!22/9/08              freeze(i_dp,k) = m2*qrainmax(i_dp)
!                 END IF

       END IF
!              freezint(i_dp) = freezint(i_dp) + freeze(i_dp,k)         &
!                   *rho_theta(i_dp,k)*(z_rho(i_dp,k+1)-z_rho(i_dp,k)) 

    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_FREEZE',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_freeze

