! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! calculate the velocity variance from the transport term
!
! Subroutine Interface:
SUBROUTINE dts_w_variance(iconvclass,n_dp,nlev,dts_ntpar,ntparmax,         &
                          z_theta,z_rho,rho,rho_theta,zlcl,zfr,ztop,       &
                          q,qse, diffmax,mb,mfr,w2lcl,wcld,wfr,wall,       &
                          scalefac1,scalefac2,scalefac3,thetav,storethvp,  &
                          ww,wwrho,massfl,massfl_rho)


! Modules used
USE dts_cntl_mod, ONLY:                                                      &
       idts_rhfac, dts_dif_pow, idts_zsc, idts_dampfac

USE dts_fitpars_mod, ONLY:                                                 &
  wamp

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! What this subroutine does:
! --------------------------
! Calculates the velocity variance: 
!1) estimate dw3dz from a scaled profile
!2) integrate up to find w3
!3) work out variance by assuming that w2 propto w3^(2/3)
!4) get massfl from w3 
! Want variance on rho levels for input into mse flux, so put dw3dz
! onto theta levels
!
! Inputs:
! Outputs:
! Weaknesses:
! -----------
! Jump at freezing level; behaviour above the freezing level
! 
! Other information: called from FLUX CALCULATION section
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  nlev                 & ! No. of model layers
 ,n_dp                 & ! No. of deep convection points
 ,dts_ntpar(n_dp)      & ! Top level of initial parcel ascent
 ,ntparmax             & ! Maximum of ntpar across all conv points
 ,iconvclass(n_dp)       ! NOT USED
       
REAL, INTENT(IN) ::     &
  z_theta(n_dp,nlev)    & ! height of theta levels
 ,z_rho(n_dp,nlev)      & ! height of rho levels(m)
 ,zlcl(n_dp)            & ! height of lifting condensation level (m)
 ,zfr(n_dp)             & ! height of freezing level (m) 
 ,ztop(n_dp)            & ! height of cloud top (m)
 ,rho(n_dp,nlev)        & ! density on rho levels
 ,rho_theta(n_dp,nlev)  & ! density on theta levels
 ,mb(n_dp)              & ! mass flux at cloud base
 ,mfr(n_dp)             & ! mass flux at freezing level    
 ,w2lcl(n_dp)           & ! velocity variance at the lcl
 ,wcld(n_dp)            & ! convective velocity scale (m/s)
 ,wfr(n_dp)             & ! conv vel scale above freezing (m/s)
 ,wall(n_dp)            & ! conv vel scale for whole layer (m/s)
 ,scalefac1(n_dp)       &
 ,scalefac2(n_dp)       &
 ,scalefac3(n_dp)       &
 ,diffmax(n_dp)         &
 ,thetav(n_dp,nlev)     & ! thetav (K)
 ,storethvp(n_dp,nlev)  &
 ,q(n_dp,nlev)          & ! water vapour (kg/kg)
 ,qse(n_dp,nlev)          ! saturation (kg/kg)

! Arguments with intent OUT:

REAL, INTENT(OUT) ::   &
  ww(n_dp,nlev)        & ! velocity variance on theta levels
 ,wwrho(n_dp,nlev)     & ! vel var on rho levels
 ,massfl(n_dp,nlev)    & ! mass flux on theta levels (m/s)
 ,massfl_rho(n_dp,nlev)  ! Mass flux on rho levels (m/s)

! Local variables

INTEGER ::        & 
  i_dp, k         & ! loop counters
, iuseoldscaling  &
, idampw2 
        
REAL ::                &
  zsc                  &  ! scaled z
 ,diff(n_dp,nlev)      &  
 ,scaleprof(n_dp,nlev) &
 ,tau                  &
 ,rhfac                &
 ,dampfac

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
        

!==============================================

! Initialise variables to zero
      
  IF (lhook) CALL dr_hook('DTS_W_VARIANCE',zhook_in,zhook_handle)

  ww(:,:) = 0.0
  wwrho(:,:) = 0.0
  massfl(:,:) = 0.0
  massfl_rho(:,:) = 0.0
  diff(:,:) = 0.0
      
  iuseoldscaling = 0 
! 1 for old scaling (different scaling above and below freezing level)
! 0 for new (single scaling for whole profile)

! idampw this is a parameter that allows the velocity variance to be linearly
! damped in a very crude way: multiply it by a linear function
! which is 1 up to the lcl, and decays linearly to zero at
! cloud top. This is not ideal, but may address the warm
! anomaly that is appearing at upper levels in aquaplanet runs
! (20/8/09).

  idampw2 = 1 ! 0 not to do this, 1 to apply dangerous damping term 


!=====================================================================
!3) work out variance

        
  DO k=1,nlev
    DO i_dp=1,n_dp
                
      IF(k <= dts_ntpar(i_dp)) THEN
! shape of w2 profile determined by undilute parcel buoyancy excess
        IF(diffmax(i_dp) > 0.0) THEN 
          diff(i_dp,k) = (storethvp(i_dp,k)-thetav(i_dp,k))/diffmax(i_dp)

          IF(diff(i_dp,k) <= 0.0) THEN
            diff(i_dp,k) = 0.0
          ELSE
!orig        diff(i_dp,k) = diff(i_dp,k)**1.5
! Sensitivity tests pass in dts_dif_pow
            diff(i_dp,k) = diff(i_dp,k)**dts_dif_pow

          END IF
        ELSE
          diff(i_dp,k) = 0.0
        END IF
                 
! scaling profile 
! This scaling profile could be modified
        IF(iuseoldscaling == 1) THEN 

          IF(z_theta(i_dp,k) < zfr(i_dp)) THEN 
            scaleprof(i_dp,k) = scalefac1(i_dp)
            IF(wcld(i_dp) <= 0.0) scaleprof(i_dp,k) = 0.0
          ELSE
              
            IF(z_theta(i_dp,k) < ztop(i_dp)) THEN
              zsc = (ztop(i_dp)-z_theta(i_dp,k))/  &
                    (ztop(i_dp)-zfr(i_dp))
            ELSE
              IF(z_theta(i_dp,k) == ztop(i_dp)) THEN
                zsc = 1.0
              END IF
            END IF
! simple interpolation between the two scale factors -- ending with
            ! scalefac2 at the top 
            scaleprof(i_dp,k) = (scalefac1(i_dp)*zsc +  &
                                         scalefac2(i_dp)*(1.0-zsc))
          END IF

        ELSE ! iuseoldscaling
          scaleprof(i_dp,k) = scalefac3(i_dp)
        END IF ! iuseoldscaling = 1

        IF(z_theta(i_dp,k) > ztop(i_dp)) THEN
          ! an adjustment timescale based on the whole layer
          scaleprof(i_dp,k) = 0.0
        END IF
        tau = 0.0
        IF(wall(i_dp) > 0.0) THEN
          tau = ztop(i_dp)/wall(i_dp)
        END IF
        scaleprof(i_dp,k) = scaleprof(i_dp,k)*tau

      END IF  ! k <= dts_ntpar(i_dp)
    END DO
  END DO

  DO k=1,nlev
    DO i_dp=1,n_dp
      IF(k < dts_ntpar(i_dp)) THEN 
        IF(z_theta(i_dp,k) > zlcl(i_dp)) THEN 

         ! let the boundary layer value decay with height according to zsc:

         ! Sensitivity tests - will slow code
          SELECT CASE (idts_zsc)
            CASE (0)  ! default original code

              zsc = exp(-(z_theta(i_dp,k)-zlcl(i_dp))  /zlcl(i_dp)) 

            CASE (1)  ! 

              zsc = exp(-0.5*(z_theta(i_dp,k)-zlcl(i_dp)) /zlcl(i_dp))  

            CASE (2)  ! 

              zsc = exp(-2.0*(z_theta(i_dp,k)-zlcl(i_dp)) /zlcl(i_dp)) 

          END SELECT

          ww(i_dp,k) = wamp*scaleprof(i_dp,k)*diff(i_dp,k) + zsc*w2lcl(i_dp)

          IF(ww(i_dp,k) < 0.0) THEN
            ww(i_dp,k) = 0.0
          END IF
                    
        ELSE
! below the lcl let it decrease to zero with a tanh function -- this
! purposefully excludes the 'bulge' of w2 in the boundary layer

          ww(i_dp,k) =  w2lcl(i_dp)*                                        &
                      (0.5*(1.0001+tanh(6.*(z_theta(i_dp,k)-0.5*zlcl(i_dp)) &
                                            /(0.5*zlcl(i_dp)) )) )**(2./3.)
          IF(ww(i_dp,k) < 0.0) THEN
            ww(i_dp,k) = 0.0
          END IF
        END IF ! diff < 0.0

      END IF
    END DO
  END DO
       
! apply a damping to w2 to reduce its strength at upper levels -- this
! is somewhat arbitrary, however. Ideally one would put in a
! dissipation term into the calculation of w2.
! Sensitivity
  idampw2=idts_dampfac  

  IF(idampw2 == 1) THEN   ! Current default
    DO k=1,nlev
      DO i_dp=1,n_dp
        IF(z_theta(i_dp,k) > zlcl(i_dp) .AND.            &
           z_theta(i_dp,k) <= ztop(i_dp)) THEN
          dampfac = (ztop(i_dp)-z_theta(i_dp,k))         &
                   /(ztop(i_dp)-zlcl(i_dp))
          ww(i_dp,k) = ww(i_dp,k)*dampfac
        END IF
      END DO
    END DO

  END IF

  IF(idampw2 == 2) THEN   ! 0.5 at top
    DO k=1,nlev
      DO i_dp=1,n_dp
        IF(z_theta(i_dp,k) > zlcl(i_dp) .AND.            &
           z_theta(i_dp,k) <= ztop(i_dp)) THEN
          dampfac = 1.0 -((z_theta(i_dp,k)-zlcl(i_dp))   &
                         /(ztop(i_dp)-zlcl(i_dp)))
          ww(i_dp,k) = ww(i_dp,k)*dampfac
        END IF
      END DO
    END DO

  END IF
  IF(idampw2 == 3) THEN   ! as 1 but raise to power of 1.5
    DO k=1,nlev
      DO i_dp=1,n_dp
        IF(z_theta(i_dp,k) > zlcl(i_dp) .AND.            &
           z_theta(i_dp,k) <= ztop(i_dp)) THEN
          dampfac = ((ztop(i_dp)-z_theta(i_dp,k))        &
                   /(ztop(i_dp)-zlcl(i_dp)) ) **1.5
          ww(i_dp,k) = ww(i_dp,k)*dampfac
        END IF
      END DO
    END DO

  END IF

! ww is on theta levels, so now create a variance on rho levels by 
! interpolating ww
  DO k=1,nlev
    IF (k == 1) THEN
      DO i_dp=1,n_dp

      ! implicitly saying ww(i_dp,0) = 0. 
        wwrho(i_dp,1) = 0.5*(ww(i_dp,1)*rho_theta(i_dp,1))/rho(i_dp,1)
      END DO
    ELSE
      DO i_dp=1,n_dp
        wwrho(i_dp,k) = 0.5*(ww(i_dp,k)*rho_theta(i_dp,k)+ww(i_dp,k-1)   &
                            *rho_theta(i_dp,k-1))/rho(i_dp,k)

      END DO
    END IF
  END DO 

!====================================================================
!4) estimate the mass flux
! first on theta levels:

  DO k=1,nlev
    DO i_dp=1,n_dp

      ! Sensitivity tests    - will slow code
      SELECT CASE(idts_rhfac)
        CASE(0)     ! Default original code    

!NB rhfac is a fudge to reduce the mass flux in low humidity environments
! physically it is designed to represent an area decrease of buoyant cores
          rhfac =q(i_dp,k)/qse(i_dp,k) ! was q/qsat

        CASE(1)     ! fixed value       
          rhfac = 1.0

        CASE(2)     !  (1+rh)/2
          rhfac = (1.0+q(i_dp,k)/qse(i_dp,k))*0.5

        CASE(3)     !  rh^2
          rhfac = (q(i_dp,k)*q(i_dp,k))/(qse(i_dp,k)*qse(i_dp,k))

      END SELECT  

      IF(ww(i_dp,k) > 0.0) THEN
        massfl(i_dp,k) = rhfac*mb(i_dp)*(ww(i_dp,k)/w2lcl(i_dp))**0.5
      END IF   
    END DO
  END DO

! calculating the mass flux on rho levels (for use in dts_qflux and
! dts_cond_and_dep)
  DO k=1,nlev

    IF (k == 1) THEN
      DO i_dp=1,n_dp
        massfl_rho(i_dp,1) = 0.5*(massfl(i_dp,1)*rho_theta(i_dp,1))      &
                               /rho(i_dp,1)
      END DO
    ELSE
      DO i_dp=1,n_dp
        massfl_rho(i_dp,k) = 0.5*(massfl(i_dp,k)*rho_theta(i_dp,k)       &
                                 +massfl(i_dp,k-1)*rho_theta(i_dp,k-1))  &
                               /rho(i_dp,k)
      END DO
    END IF 
  END DO
  IF (lhook) CALL dr_hook('DTS_W_VARIANCE',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_w_variance
