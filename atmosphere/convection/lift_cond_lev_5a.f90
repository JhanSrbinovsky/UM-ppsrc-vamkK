! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Lifting condensation level calculation
!
! Subroutine Interface:
SUBROUTINE lift_cond_lev (npnts, nlev, k_plume,                          &
                          l_mixing_ratio,                                & 
                          pstar, q, T,                                   &
                          p_theta_lev, exner_rho, z_rho,                 &
                          T_lcl, p_lcl, z_lcl, qsat_lcl )

USE atmos_constants_mod, ONLY: kappa, pref, repsilon, recip_kappa

USE cv_diag_param_mod, ONLY:                                             &
    a_bolton, b_bolton, c_bolton, d_bolton

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates the lifting condensation level (LCL) temperature,
!   pressure and height.
!
!  Is designed to work on compressed arrays for just a selected set of points.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,nlev                   ! Number of model levels for calculations

INTEGER, INTENT(IN) :: &
  k_plume(npnts)         ! Starting model level for plume ascent

LOGICAL, INTENT(IN) :: & 
  l_mixing_ratio         ! .TRUE. if input q a mixing ratio otherwise 
                         ! assumes input q a specific humidity.
  
REAL, INTENT(IN) ::       &
  pstar(npnts)            & ! Surface pressure (Pa)
 ,q(npnts,nlev)           & ! water vapour on model levels (kg/kg)
 ,T(npnts,nlev)           & ! Temperature on model levels (K)
 ,p_theta_lev(npnts,nlev) & ! Pressure on theta levels (Pa)
 ,exner_rho(npnts,nlev)   & ! Exner Pressure on rho levels 
 ,z_rho(npnts,nlev)         ! Hieght of rho levels  (m)

REAL, INTENT(OUT) ::      &
  T_lcl(npnts)            & ! Temperature of LCL  (K)
 ,p_lcl(npnts)            & ! Pressure of LCL  (Pa)
 ,z_lcl(npnts)            & ! Height of LCL  (m)
 ,qsat_lcl(npnts)           ! qsaturation at zlcl (i.e. cloud base) kg/kg

!-------------------------------------------------------------------------------
! Local variables

INTEGER ::               & 
  i,k                      ! loop counter

REAL ::                  &
  exner_lcl(npnts)       & ! Exner pressure at LCL
 ,exner_surf(npnts)        ! Exner pressure at surface

REAL ::                  &
  factor                 & ! factor used in interpolation
 ,vap_press                ! vapour pressure

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('LIFT_COND_LEV',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
!     using approximations from Bolton (1980)
!-------------------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon       q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio
!-------------------------------------------------------------------------------

If (l_mixing_ratio) THEN   ! expression for mixing ratio

  DO i=1, npnts

    vap_press = 0.01*q(i,k_plume(i)) * p_theta_lev(i,k_plume(i))            &
                                      / (repsilon+q(i,k_plume(i)) )
    IF (vap_press  >   0.0) THEN
      T_lcl(i) = a_bolton + b_bolton/                                       &
                              (c_bolton*LOG(T(i,k_plume(i)))                &
                                         - LOG(vap_press) - d_bolton )

      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                &
                     ( T_lcl(i) / T(i,k_plume(i)) )**recip_kappa

    ELSE
      p_lcl(i) = pstar(i)
    END IF
 
  END DO     ! i loop

ELSE          ! expression for specific humidity

  DO i=1, npnts  
    vap_press = q(i,k_plume(i)) *                                           &
                       p_theta_lev(i,k_plume(i)) / ( 100.0*repsilon )
    IF (vap_press  >   0.0) THEN
      T_lcl(i) = a_bolton + b_bolton/                                       &
                         (c_bolton*LOG(T(i,k_plume(i)))                     &
                                         - LOG(vap_press) - d_bolton )
      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                &
                        ( T_lcl(i) / T(i,k_plume(i)) )**recip_kappa
    ELSE
      p_lcl(i) = pstar(i)
    END IF
 
  END DO

END IF ! test on l_mixing_ratio  

! work out qsat at LCL
! DEPENDS ON: qsat_mix
Call qsat_mix(qsat_lcl,T_lcl,p_lcl,npnts,l_mixing_ratio)
      

!-----------------------------------------------------------------------
! Accurate calculation of height of LCL using exner_lcl rather than p_lcl
! as UM interpolates exner linearly in height but not pressure.
!-----------------------------------------------------------------------

DO i=1, npnts
  exner_lcl(i)  = (p_lcl(i)/pref)**kappa
  exner_surf(i) = (pstar(i)/pref)**kappa     
END DO
      

k = 1         

  DO i=1, npnts
    IF ( exner_lcl(i) >= exner_surf(i)) THEN 
        z_lcl(i) = 0.0           ! at or below surface
    ELSE IF (exner_lcl(i) < exner_surf(i)                                 &
                     .and. exner_lcl(i) > exner_rho(i,k)) THEN
      factor= (exner_rho(i,k) - exner_lcl(i))/                            &
                        (exner_rho(i,k) - exner_surf(i))
      z_lcl(i) = (1.0-factor)*z_rho(i,k)  
    END IF
  END DO         ! points

DO k=2,nlev

  DO i=1, npnts
    IF (exner_lcl(i) >= exner_rho(i,k)                                    &
                        .and. exner_lcl(i) < exner_rho(i,k-1) ) THEN    
      factor= (exner_rho(i,k) - exner_lcl(i))/                            &
                        (exner_rho(i,k) - exner_rho(i,k-1))
      z_lcl(i) = (1.0-factor)*z_rho(i,k)+factor*z_rho(i,k-1)
    END IF   
  END DO         ! points
  
END DO           ! level loop

! Check z_lcl not less than a minimum value 
! Note z_lcl currently used by diurnal cycle diagnosis and 
! also by Deep turbulence scheme. The diurnal cycle diagnosis will be 
! happy with a value of zero but this will not work for the deep 
! turbulence scheme.
! Set z_lcl lowest model depth - fix may imply some model resolution dependence.

DO i=1, npnts 
  z_lcl(i) =MAX(z_lcl(i),z_rho(i,2))
END DO 


!-------------------------------------------------------------------------------
IF (lhook) CALL dr_hook('LIFT_COND_LEV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lift_cond_lev
