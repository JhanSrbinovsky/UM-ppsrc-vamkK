! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ calculates sublimation rate
!

SUBROUTINE dts_sublimation(n_dp,nlev,dts_ntpar,ntparmax,               &
                           dr_across_th,rho_theta,q,temperature,qse,   &
                           deposition,depint,freezint,qcfint,          &
                           sublimation,subint)

USE dts_cntl_mod, ONLY:                                                  &
        dts_gamma_fac

USE water_constants_mod, ONLY: tm
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! Calculates sublimation rate
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  nlev                 & ! No. of model layers
 ,n_dp                 & ! Number of deep points
 ,dts_ntpar(n_dp)      & ! Top level 
 ,ntparmax               ! Maximum top cloud level

REAL, INTENT(IN)    ::          &
  dr_across_th(n_dp,nlev)       & ! thickness of theta layers (m)
 ,rho_theta(n_dp,nlev)          & ! Density on theta levels (kg/m3)
 ,q(n_dp,nlev)                  & ! Model mixing ratio  (kg/kg)
 ,temperature(n_dp,nlev)        & ! Temperature (K)
 ,qse(n_dp,nlev)                & ! qsat   (kg/kg) 
 ,deposition(n_dp,nlev)         & ! Deposition
 ,depint(n_dp)                  & ! Vertical integral of deposition 
 ,freezint(n_dp)                & ! vertical integral of freezing
 ,qcfint(n_dp)                    ! Vertical integral of qcf

REAL, INTENT(OUT) ::            &
  sublimation(n_dp,nlev)        & ! sublimation       
 ,subint(n_dp)                    ! Vertical integral of sublimation

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k             ! Loop counters

REAL ::             &
  gamma_fac         &
 ,tmpint

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! The amplitude of the sublimation term is given by assuming equilibrium 
! for the ice budget, i.e.:
!
!      dqice/dt = dep - sub - fall - d<w'qice'>/dz - PC2 = 0        
!
! and that as a result 
!  int (d<w'qice'>/dz) dz = int( dep - sub - fall - PC2)dz
!
! and since we know that <w'qice'> must be zero at the top of the
! convective layer, we can say that  int(dep - sub - fall - PC2)dz= 0 
! and therefore: int(sub dz) = int(dep-fall-PC2)dz
!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('DTS_SUBLIMATION',zhook_in,zhook_handle)
  sublimation(:,:) = 0.0
  subint(:) = 0.0
        
! SAJS: variable which may end up being a function of height
! SAJS: this is needed in the experimental case where sublimation isn't 
! renormalised 
       
!orig  gamma_fac = 0.125 ! 0.4 ! designed not to fall below 60% of deposition

! Sensitivity tests
  gamma_fac = dts_gamma_fac
              
  DO K=1,ntparmax
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp) .AND. k < nlev) THEN 
        IF(temperature(i_dp,k) < tm) THEN 

        ! Say sublimation looks a bit like deposition profile, but with
        ! enhanced values in low rh areas
          IF(qse(i_dp,k) > 0.0) THEN 
            sublimation(i_dp,k) = deposition(i_dp,k)*                       &
                                    (1.-gamma_fac*q(i_dp,k)/qse(i_dp,k)) 
          ELSE
            write(6,*) 'qisat alarm: ',qse(i_dp,k)
            sublimation(i_dp,k) = deposition(i_dp,k)
          END IF
          IF(sublimation(i_dp,k) < 0.0) THEN
            sublimation(i_dp,k) = 0.0
          END IF
                    
          subint(i_dp) = subint(i_dp) + sublimation(i_dp,k)* &
                                dr_across_th(i_dp,k)*rho_theta(i_dp,k)
        END IF
      END IF

    END DO ! i_dp
  END DO ! k

! Renormalise if sublimation term is greater than creation terms 
  DO i_dp=1,n_dp
    tmpint = depint(i_dp)+freezint(i_dp)-qcfint(i_dp)
    IF(subint(i_dp) > tmpint) THEN
 
      IF(tmpint > 0.0) THEN
        DO K=1,nlev 
          sublimation(i_dp,k) = sublimation(i_dp,k)*tmpint/subint(i_dp)
        END DO
        subint(i_dp) = tmpint
      ELSE
        sublimation(i_dp,:) = 0.0
        subint(i_dp) = 0.0
      END IF
    END IF
  END DO
  IF (lhook) CALL dr_hook('DTS_SUBLIMATION',zhook_out,zhook_handle)
  RETURN


END SUBROUTINE dts_sublimation
