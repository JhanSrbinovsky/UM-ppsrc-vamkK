! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the gradient of thetav and the moist adiabatic thetav
!

SUBROUTINE dts_dthvdz(n_dp,nlev,dts_ntpar,ntparmax                         &
                      ,qse,exner_layer_centres,z_theta,dr_across_rh,thetav &
                      ,storethvp,dthvdz,dthvdz_m,dqsedz)
           
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! subroutine to calculate the gradient of thetav, and the gradient of
! the moist adiabatic thetav. For use in dts_wthv, and
! dts_w_variance
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  n_dp                 & ! No. of deep convection points
 ,nlev                 & ! No. of model layers
 ,dts_ntpar(n_dp)      & ! Top of parcel ascent for deep TCS
 ,ntparmax 

REAL, INTENT(IN)    ::             &
  qse(n_dp,nlev)                   & ! qsat of environment (kg/kg)
 ,exner_layer_centres(n_dp,0:nlev) & ! Exner at theta levels
 ,z_theta(n_dp,nlev)               & ! heights of theta levels (m)
 ,dr_across_rh(n_dp,nlev)          & ! thickness of rho layers (m)
 ,thetav(n_dp,nlev)                  ! thetav (K)

REAL, INTENT(INOUT) ::  &
  storethvp(n_dp,nlev)    ! parcel virtual pot temp
         
REAL, INTENT(OUT) ::    &
  dthvdz_m(n_dp,nlev)   & ! d(thetav)/dz of parcel thetav   (K/m)
 ,dthvdz(n_dp,nlev)     & ! d(thetav)/dz                    (K/m)
 ,dqsedz(n_dp,nlev)       ! dqse/dz on rho levels  (kg/kg/m)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::       &
  i_dp,k             ! loop counters 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


 
!-----------------------------------------------------------------------
! diffuses back to a moist adiabat
 IF (lhook) CALL dr_hook('DTS_DTHVDZ',zhook_in,zhook_handle)
  dthvdz_m(:,:) = 0.0 
  dthvdz(:,:)   = 0.0 
  dqsedz(:,:)   = 0.0
    
        
!  k = 1   ! Special case  Note initialisation of arrays above does this
! will use the flux on rho level 2 instead
!  DO i_dp=1,n_dp
!    dthvdz(i_dp,k)   = 0.0
!    dthvdz_m(i_dp,k) = 0.0
!  END DO

! All other levels

  DO k=2,ntparmax 
    DO i_dp=1,n_dp
                    
      dthvdz(i_dp,k) = (thetav(i_dp,k)-thetav(i_dp,k-1))/dr_across_rh(i_dp,k)

      dqsedz(i_dp,k) = (qse(i_dp,k)   -qse(i_dp,k-1)   )/dr_across_rh(i_dp,k)

!NB! or use value derived from parcel ascent in cape calculation:
      dthvdz_m(i_dp,k) = (storethvp(i_dp,k)-storethvp(i_dp,k-1))/          &
                                                  dr_across_rh(i_dp,k)

! Fudge to ensure that this never drops beneath zero -- should really be done
! using a more accurate parcel ascent...
      IF(dthvdz_m(i_dp,k) < 0.0) THEN
        dthvdz_m(i_dp,k)  = 0.0
        storethvp(i_dp,k) = storethvp(i_dp,k-1)
      END IF

    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_DTHVDZ',zhook_out,zhook_handle)
  RETURN


END SUBROUTINE dts_dthvdz
