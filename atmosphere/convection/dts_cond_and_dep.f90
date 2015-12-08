! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ calculates the condensation
!

SUBROUTINE dts_cond_and_dep(n_dp,nlev,ntparmax,dts_ntpar              &
                            ,zlcl,zfr,rho,rho_theta,dr_across_rh      &
                            ,dr_across_th,qsat_moist_ad,qse,dqsedz,q  &
                            ,massfl_rho,temperature,z_theta,z_rho     &
                            ,condensation,deposition,condint,depint)
       
USE dts_fitpars_mod, ONLY :                                           &
  tcritplume

USE conversions_mod, ONLY: zerodegc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! calculates the condensation
!
! Weaknesses:
! -----------
! -existence of 'fudge' instead of a better approx for massfl
! -extreme cut between condensation and deposition
!
! Other information: 
! ------------------
! called from MICROPHYSICS PARAMETER section of deep_turb_conv
!
!   Called by GLUE_CONV.
!
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
 ,ntparmax                ! Max ntpar across all conv points
     
REAL, INTENT(IN)    ::                &
  temperature(n_dp,nlev)     & ! temperature in K
 ,z_theta(n_dp,nlev)         & ! height of theta levels(m)
 ,z_rho(n_dp,nlev)         & ! height of rho levels(m)
 ,zlcl(n_dp)         & ! height of lifting condensation level
 ,zfr(n_dp)          &
 ,qse(n_dp,nlev)             &
 ,dqsedz(n_dp,nlev)          &
 ,q(n_dp,nlev)               &
 ,massfl_rho(n_dp,nlev)          & ! on thetat levels
 ,rho(n_dp,nlev)             & ! density on rho levels  (kg/m3)
 ,rho_theta(n_dp,nlev)       &  ! density on theta levels (kg/m3)
 ,dr_across_rh(n_dp,nlev)    &  ! thickness of rho layers  (m)
 ,dr_across_th(n_dp,nlev)    &  ! thickness of theta layers  (m)
 ,qsat_moist_ad(n_dp,nlev)    ! qsat along parcel ascent

REAL, INTENT(OUT) ::         &
  condensation(n_dp,nlev)    & ! condensation rate on theta levs kg/kg/s
 ,deposition(n_dp,nlev)      & ! deposition rate on theta levs kg/kg/s
 ,condint(n_dp)              & ! Integral of condensation
 ,depint(n_dp)                 ! integral of deposition

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k              ! loop counters

REAL  ::           &
  fac_ctod         &
 ,frac             &
 ,cond(n_dp,nlev)  &  ! condensation
 ,dep(n_dp,nlev)      ! deposition

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('DTS_COND_AND_DEP',zhook_in,zhook_handle)

  condensation(:,:) = 0.0
  deposition(:,:) = 0.0
  cond(:,:) = 0.0
  dep(:,:) = 0.0
  condint(:) = 0.0
  depint(:) = 0.0

! Calculates -m dqsat(along a moist adiabat)/dz
! will reset qsat_moist_ad back to tmp at the end of the routine

  DO k=2,ntparmax
    DO i_dp=1,n_dp
      ! only condense above zlcl
            
      IF(k <= dts_ntpar(i_dp) .AND. z_rho(i_dp,k) >= zlcl(i_dp)          &
                                                 .AND. k < nlev) THEN 
! (dividing by density because the mass flux contains it)
! for the level just above zlcl, use the gradient from the layer above
! to determine dqsat/dz, and only apply to the
! fraction of the layer above zlcl (as determined by frac)
! this part is designed to avoid too great a step function just above
! the lcl
        IF(qsat_moist_ad(i_dp,k) > 0.) THEN 
          IF(z_rho(i_dp,k-1) < zlcl(i_dp)) THEN
            frac = (z_rho(i_dp,k)-zlcl(i_dp))/dr_across_th(i_dp,k-1) 
                
! slightly different gradient if crossing lcl

            cond(i_dp,k)=-frac*(massfl_rho(i_dp,k)/rho(i_dp,k))*           &
                         (qsat_moist_ad(i_dp,k+1)*rho_theta(i_dp,k+1)      &
                         -qsat_moist_ad(i_dp,k)*rho_theta(i_dp,k))/        &
                         dr_across_rh(i_dp,k)/rho(i_dp,k) ! on rho levels

! once more than one level above zlcl, simply apply -m dqsat/dz
          ELSE
            cond(i_dp,k)=-(massfl_rho(i_dp,k)/rho(i_dp,k))*                &
                         (qsat_moist_ad(i_dp,k)*rho_theta(i_dp,k)          &
                         -qsat_moist_ad(i_dp,k-1)*rho_theta(i_dp,k-1))/    &
                         dr_across_rh(i_dp,k)/rho(i_dp,k) ! on rho levels

          END IF

          ! safety check -- set to zero if less than zero
          ! otherwise could end up with a spiral of cooling
          IF(cond(i_dp,k) < 0.0) THEN
            cond(i_dp,k) = 0.0
          END IF

        END IF

! linear decrease in condensation down to -10C, with linearly
! increasing deposition in this range
! nb need to ensure consistent with plume liquid water content in
! dts_pc2,dts_cape

        IF(z_theta(i_dp,k) > zfr(i_dp)) THEN
! this should be zero at T = 273.14 K, and 1 at T = 263.14K
          fac_ctod = (temperature(i_dp,k)-ZeroDegC)/tcritplume
          IF(fac_ctod > 0.0 .AND. fac_ctod <= 1.0) THEN

            dep(i_dp,k) = fac_ctod*cond(i_dp,k)
            cond(i_dp,k) = (1.0-fac_ctod)*cond(i_dp,k)
          ELSE
            IF(fac_ctod <= 0.0) THEN ! ie warmer than 273.15 K
              dep(i_dp,k) = 0.0
              !cond remains unchanged
            END IF
            IF(fac_ctod > 1.0) THEN ! ie colder than -10C

              dep(i_dp,k) = cond(i_dp,k)   !   *(1.+rhice) Possible alternative
              cond(i_dp,k) = 0.0
            END IF
          END IF

        END IF 
      END IF ! qsat_moist_ad > 0.
    END DO ! i_dp
  END DO ! k

! now put onto theta levels
  DO k=1,nlev-1
    DO i_dp=1,n_dp
      condensation(i_dp,k) = 0.5*(cond(i_dp,k+1)*rho(i_dp,k+1)              &
                                 +cond(i_dp,k)  *rho(i_dp,k))/rho_theta(i_dp,k)  
      deposition(i_dp,k) = 0.5*(dep(i_dp,k+1)*rho(i_dp,k+1)                 &
                               +dep(i_dp,k)  *rho(i_dp,k)   )/rho_theta(i_dp,k) 

! vertical integrals of condensation and deposition:

      condint(i_dp) = condint(i_dp) + condensation(i_dp,k)                  &
                    *rho_theta(i_dp,k)*dr_across_th(i_dp,k)
                  
      depint(i_dp) = depint(i_dp) + deposition(i_dp,k)                      &
                    *rho_theta(i_dp,k)*dr_across_th(i_dp,k)

    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_COND_AND_DEP',zhook_out,zhook_handle)
  RETURN


END SUBROUTINE dts_cond_and_dep
