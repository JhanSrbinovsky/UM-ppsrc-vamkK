! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  Calculate parameters used to scale the profiles
!

SUBROUTINE dts_flux_par(n_dp,nlev,ntml,freeze_lev,km40,land_mask,            &
                        q,qse,temperature,dr_across_th,zfr,w_max,zlcl,ztop,  &
                        timestep,wstar,cape_below_fr,                        &
                        cape_above_fr,cape_whole_layer,ql_ad,h_ad,           &
                        mb,mfr,w2lcl,wcld,wfr,wall,cstar,                    &
                 qrainmax,qsatsurf,qvatfr,sigma,scalefac1,scalefac2,scalefac3)

USE dts_cntl_mod, ONLY:                                                        &
      idts_mb, dts_mb

USE dts_fitpars_mod, ONLY:                                                   &
  a1, a2, a3, b1, b2, b3, c1, c2, e1, e2, f1, f2, land1, land2

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! Purpose: calculates a number of parameters that are used to scale the
! flux profiles 
!
! Variants that could be tried:
! Other information: called from SCALING PARAMETER section
!
! Called by: deep_turb_conv
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
   n_dp                & ! number of deep points 
  ,nlev                & ! number of levels 
  ,freeze_lev(n_dp)    & ! Level index for freezing level
  ,ntml(n_dp)          & ! start level for convection
  ,km40(n_dp)            ! level index for -40C

LOGICAL, INTENT(IN) :: & 
  land_mask(n_dp)        ! Land/sea mask .true. if land

REAL, INTENT(IN) :: &
  timestep             ! model convection timestep

REAL, INTENT(IN)    ::        &
  wstar(n_dp)                 & ! sub cloud velocity scale (m/s)
 ,cape_below_fr(n_dp)         & ! CAPE below freezing level
 ,cape_above_fr(n_dp)         & ! CAPE above freezing level
 ,cape_whole_layer(n_dp)      & ! total CAPE
 ,ql_ad(n_dp)                 & !
 ,h_ad(n_dp)                  & !
 ,ztop(n_dp)                  & ! cloud top height (m)
 ,zlcl(n_dp)                  & ! height of lifting condensation level (m)
 ,zfr(n_dp)                   & ! height of freezing level (m)
 ,w_max(n_dp)                 & ! Max w in column
 ,temperature(n_dp,nlev)      & ! temperature in K
 ,dr_across_th(n_dp,nlev)     & ! layer depth of theta layers (m)
 ,qse(n_dp,nlev)              & ! qsat environment (kg/kg)
 ,q(n_dp,nlev)                  ! water vapour (kg/kg)
           
REAL, INTENT(OUT) ::  &
  mb(n_dp)            & ! mass flux at cloud base
 ,mfr(n_dp)           & ! mass flux at freezing level
 ,w2lcl(n_dp)         & ! velocity variance at lcl
 ,wcld(n_dp)          & ! velocity scale for convective layer
 ,wfr(n_dp)           & ! velocity scale above
 ,wall(n_dp)          & ! velocity scale for whole layer freezing level
 ,cstar(n_dp)         & ! scale used for c-e (may change...) NOT USED ?
 ,sigma(n_dp)         & ! the fr cloud area at cloud base
 ,qrainmax(n_dp)      & ! scaling param for max rain amount
 ,qsatsurf(n_dp)      & ! qsat at the surface (kg/kg)
 ,qvatfr(n_dp)        & ! q at freezing level (kg/kg)
 ,scalefac1(n_dp)     & ! scale factor below freezing lev
 ,scalefac2(n_dp)     & ! scale factor above freezing lev
 ,scalefac3(n_dp)       ! scale factor for whole layer

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER ::    &
  i_dp,k       ! loop counters

REAL ::       &      
  smallcape   &
 ,w_fac       &
 ,cfl_limit   &
 ,temp_value  &  
 ,temp_value2  

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DTS_FLUX_PAR',zhook_in,zhook_handle)

  smallcape = 1.0e-6 ! a very small cape value

! mass flux at cloud base
! Sensitivity test options

  SELECT CASE (idts_mb)
    CASE (0)  ! original code     
        
      DO i_dp=1,n_dp
        mb(i_dp) = a1*wstar(i_dp) + a2
! for low values of wstar use the shallow closure mb = 0.04 wstar
! allows mb to go to zero smoothly with wstar
        IF(a3*wstar(i_dp) > mb(i_dp)) THEN
          mb(i_dp) = a3*wstar(i_dp)
        END IF 
      END DO

    CASE (1)  !  dts_mb *wstar
        
      DO i_dp=1,n_dp
        mb(i_dp) = dts_mb*wstar(i_dp)
      END DO

  END SELECT

! velocity variance at cloud base

  temp_value = a2/(a3-a1)
  temp_value2 = (temp_value)**2

  DO i_dp=1,n_dp
    w2lcl(i_dp) = c1*wstar(i_dp)**2 + c2
! at the same point as the mb transitions onto a shallower slope,
! linearly interpolate w2lcl down to zero too
! this happens when a1*wstar + a2 = a3 wstar; ie wstar = a2/(a3-a1)
! use the value that w2lcl had at this value to determine the gradient

    
    IF(wstar(i_dp) < temp_value) THEN
      w2lcl(i_dp) = (c1*temp_value2+c2)*(wstar(i_dp)**2/(temp_value2))
    END IF
  END DO


! wcld -- the velocity scale for the convective layer

  DO i_dp=1,n_dp
    IF(cape_below_fr(i_dp) > smallcape) THEN 
      wcld(i_dp) = (mb(i_dp)*cape_below_fr(i_dp))**(1./3.)
      !nb or cape_whole_layer?
    ELSE
      wcld(i_dp) = 0.0
    END IF
  END DO

! wall -- velocity scale for the whole layer

  DO i_dp=1,n_dp
    IF(cape_whole_layer(i_dp) > smallcape) THEN 
      wall(i_dp) = (mb(i_dp)*cape_whole_layer(i_dp))**(1./3.)
    ELSE 
      wall(i_dp) = 0.0
    END IF
  END DO

! mass flux at the freezing level

  mfr(:) = b1+b2*wall(:)

! allow it to go to zero linearly with wall
  DO i_dp=1,n_dp
    IF(mfr(i_dp) < b3*wall(i_dp)) THEN
      mfr(i_dp) = b3*wall(i_dp)
    END IF 
  END DO

  DO i_dp=1,n_dp
    IF(cape_above_fr(i_dp) > smallcape) THEN 
      wfr(i_dp) = (mfr(i_dp)*cape_above_fr(i_dp))**(1./3.) 
    ELSE
      wfr(i_dp) = 0.0
    END IF
  END DO

!Initialise fields
  scalefac1(:) = 0.0
  scalefac2(:) = 0.0
  scalefac3(:) = 0.0

  DO i_dp = 1,n_dp 
    IF(wcld(i_dp) > 0.0 .and. zfr(i_dp) > zlcl(i_dp)) THEN 
      scalefac1(i_dp) =  ( (mb(i_dp)/wcld(i_dp))**0.5 )*                &
                                  (wcld(i_dp)**3)/zfr(i_dp)
    END IF
    !nb this second scale factor may need changing
    IF(wfr(i_dp) > 0.0 .and. ztop(i_dp) > zfr(i_dp)) THEN 
      scalefac2(i_dp) =  ( (mfr(i_dp)/wfr(i_dp))**0.5 )*                &
                                    (wfr(i_dp)**3)/ztop(i_dp) 
              ! Scaled with: (ztop(i_dp)-zfr(i_dp)), but dangerous when
              ! ztop close to zfr
    END IF
    IF(wall(i_dp) > 0.0 .and. ztop(i_dp) > zlcl(i_dp)) THEN 
      scalefac3(i_dp) =  ( (mb(i_dp)/wall(i_dp))**0.5 )*(wall(i_dp)**3) &
                                                   /ztop(i_dp)
    END IF
 
  END DO

! cstar nb not used, so setting to zero
!! calculate cstar cf Alan's latest paper 11/10/07
!mb(:)*ql_ad(:)/(h_ad(:)-zlcl(:))

  cstar(:) = 0.0 

! sigma - fractional area at cloud base

  sigma(:) = f1*wstar(:)**f2

! qrainmax (may want this to be land/sea dependent)
  DO i_dp=1,n_dp
    IF(land_mask(i_dp)) THEN 
      ! land param -- nb this is probably overly simplistic...
      qrainmax(i_dp) = land1 + land2*sigma(i_dp) 
    ELSE
      ! sea param 
      qrainmax(i_dp) = e1 + e2*sigma(i_dp)
    END IF
    IF(qrainmax(i_dp) < 0.0) THEN
      qrainmax(i_dp) = 0.0
    END IF
  END DO
     
! qsatsurf
  qsatsurf(:) = 0.0

! qvatfr
  qvatfr(:) = 0.0
  IF (lhook) CALL dr_hook('DTS_FLUX_PAR',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_flux_par
