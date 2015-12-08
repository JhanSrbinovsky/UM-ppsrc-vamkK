! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  Uses a very simple gradient formula for the water vapour flux, and
! derive theta from difference between mse and q flux
!

SUBROUTINE dts_qflux(n_dp,nlev,dts_ntpar,ntparmax,z_theta,z_rho,rho   &
                 ,rho_theta,dr_across_rh, dr_across_th                &
                 ,wwrho,massfl,massfl_rho,wstar,theta                 &
                 ,q,qsewat,zlcl,ztop,wq0,w2lcl,timestep,wqv,q_withflux)

USE dts_fitpars_mod, ONLY:                                            &
  h1, j1


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! Uses a very simple gradient formula for the water vapour flux, and
! derive theta from difference between mse and q flux
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  nlev                 & ! No. of model layers
 ,n_dp                 & ! No. of deep convection points
 ,dts_ntpar(n_dp)      & ! Top level of initial parcel ascent
 ,ntparmax  

REAL, INTENT(IN)    ::     &
  z_theta(n_dp,nlev)       & ! height of theta levels (m)
 ,z_rho(n_dp,nlev)         & ! height of rho levels(m)
 ,zlcl(n_dp)               & ! Lifting condensation level (m) 
 ,ztop(n_dp)               & ! top of convection (m) 
 ,theta(n_dp,nlev)         & ! potential temperature (K)
 ,q(n_dp,nlev)             & ! water vapour (kg/kg)
 ,wwrho(n_dp,nlev)         & ! velocity variance (m^2/s^-2)
 ,massfl_rho(n_dp,nlev)    & ! mass flux on rho levels
 ,massfl(n_dp,nlev)        & ! mass flux on theta levels
 ,wq0(n_dp)                & ! surface q flux
 ,w2lcl(n_dp)              & ! w2 at lcl
 ,wstar(n_dp)              & ! B-L convective velocity scale (m/s)
 ,timestep                 & ! timstep  (s)
 ,rho(n_dp,nlev)           & ! Density on rho levels (kg/m3) 
 ,rho_theta(n_dp,nlev)     & ! Density on theta levels (kg/m3)
 ,dr_across_rh(n_dp,nlev)  & ! rho layer thickness (m)
 ,dr_across_th(n_dp,nlev)  & ! theta layer thickness (m)
 ,qsewat(n_dp,nlev)          ! qsaturation wrt water (kg/kg)

REAL, INTENT(OUT) ::     &
  wqv(n_dp,nlev)         & ! water vapour flux -- estimated in this
 ,q_withflux(n_dp,nlev)    ! new incremented qv


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  k,i_dp           & ! loop counters
 ,ntparmaxp1         

INTEGER ::         &
  nlevs_arr(n_dp)

REAL ::                 &
  dqvdz(n_dp,nlev)      & ! vertical gradient of mse (on rho levels)
 ,aa(n_dp,nlev)         & ! h(k-1) coefficients
 ,bb(n_dp,nlev)         & ! h(k) coefficients
 ,cc(n_dp,nlev)         & ! h(k+1) coefficients
 ,rr(n_dp,nlev)         & ! r.h.s. of implicit mse equn
 ,kprof(n_dp,nlev)      & ! h1 ww
 ,dznk (n_dp,nlev)      & ! layer thicknesses 1 
 ,drmqsedz(n_dp,nlev)   &
 ,qse_valrho(n_dp,nlev) & ! qsat on rho levels (kg/kg)
 ,q_rho(n_dp,nlev)        ! q on rho levels (kg/kg)


REAL ::        &
  temp_term    & ! try to speed up code  
 ,dz_term

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
        
! Initialise fields
  IF (lhook) CALL dr_hook('DTS_QFLUX',zhook_in,zhook_handle)
  q_rho(:,:) = 0.0
  q_withflux(:,:) = q(:,:)
  dqvdz(:,:) = 0.0
  wqv(:,:) = 0.0
  qse_valrho(:,:) = 0.0
  drmqsedz(:,:) = 0.0

  nlevs_arr(:) = nlev

! Calculate q gradient
  DO k=2,ntparmax+1
    DO i_dp=1,n_dp
      ! dqvdz is on rho levels
      dqvdz(i_dp,k)= (q(i_dp,k)-q(i_dp,k-1))/dr_across_rh(i_dp,k)

      qse_valrho(i_dp,k) = 0.5*(qsewat(i_dp,k)*rho_theta(i_dp,k)+        &
                                  qsewat(i_dp,k-1)*rho_theta(i_dp,k-1))/ &
                                    rho(i_dp,k)

      ! below the lcl, this term should represent the water vapour excess as
      ! a result of non-local lifting (taken to be the water
      ! vapour value at the lowest theta level) 
            
                 
      q_rho(i_dp,k) = 0.5*(q(i_dp,k)*rho_theta(i_dp,k)+                &
                                q(i_dp,k-1)*rho_theta(i_dp,k-1))/      &
                                 rho(i_dp,k)            

    END DO
  END DO
       
!nb  should perhaps interpolate between a surface value here
  DO i_dp=1,n_dp
    temp_term = rho(i_dp,1)/rho_theta(i_dp,1)
    q_rho(i_dp,1)      = q(i_dp,1)*temp_term
    qse_valrho(i_dp,1) = qsewat(i_dp,1)*temp_term
  END DO

! now put onto theta levels
  ntparmaxp1 = ntparmax+1 
  IF (ntparmaxp1 > nlev) THEN
    ntparmaxp1 = nlev-1 
  END IF 

  k=1 
! Uses Correct heights - don't want bottom layer thick here as using values
! valid at bottom rho level
    DO i_dp=1,n_dp
      drmqsedz(i_dp,k)=(wq0(i_dp)/(w2lcl(i_dp)*qse_valrho(i_dp,1)))*     &
                      (wwrho(i_dp,k+1)*qse_valrho(i_dp,k+1)*rho(i_dp,k+1)&
                       -wwrho(i_dp,k)*qse_valrho(i_dp,k)*rho(i_dp,k))/   & 
                      (z_rho(i_dp,k+1)-z_rho(i_dp,k))

    END DO
  DO k=2,ntparmaxp1 
    DO i_dp=1,n_dp
      drmqsedz(i_dp,k)=(wq0(i_dp)/(w2lcl(i_dp)*qse_valrho(i_dp,1)))*     &
                      (wwrho(i_dp,k+1)*qse_valrho(i_dp,k+1)*rho(i_dp,k+1)&
                       -wwrho(i_dp,k)*qse_valrho(i_dp,k)*rho(i_dp,k))/   & 
                                dr_across_th(i_dp,k)
    END DO
  END DO
                   
! rho factor because dx/dt = -1/rho (d rho wX/dz)
  DO k=1,nlev 
    DO i_dp=1,n_dp
      kprof(i_dp,k) = h1*(1.0-z_rho(i_dp,k)/ztop(i_dp))*wwrho(i_dp,k) &
                               *rho(i_dp,k)     ! on rho levels
    END DO
  END DO
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

! set kprof to zero above convective top
  DO i_dp=1,n_dp
    kprof(i_dp,dts_ntpar(i_dp):nlev) = 0.0
  END DO

!array ends set here
  dznk(:,1) = z_theta(:,1)


! nb this would benefit from someone checking...
! this is intended for above the convecting layer
  aa(:,:) = 0.0
  bb(:,:) = 1.0 
  cc(:,:) = 0.0
  rr(:,:) = 0.0
! now deal with convective layer
  ntparmaxp1 = ntparmax 
  IF (ntparmaxp1 > nlev) THEN
    ntparmaxp1 = nlev-1 
  END IF 

  DO k=2,ntparmaxp1 
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN 

        temp_term = timestep/dr_across_th(i_dp,k)/rho_theta(i_dp,k)

        ! k-1 coef:
        aa(i_dp,k) = -kprof(i_dp,k)*temp_term/dr_across_rh(i_dp,k)
                                   
        ! k+1 coef:
        cc(i_dp,k) = -kprof(i_dp,k+1)*temp_term/dr_across_rh(i_dp,k+1)

        ! k coef:  
        bb(i_dp,k) = 1.0 - aa(i_dp,k) -cc(i_dp,k)
                      
      END IF
    END DO
  END DO

       
! now do k=1 level
  k=1
  DO i_dp=1,n_dp
    temp_term =  timestep/dr_across_th(i_dp,k)/rho_theta(i_dp,1)
    dz_term    = 1.0/dr_across_rh(i_dp,2)

    aa(i_dp,1) = 0.0

    cc(i_dp,1) = - temp_term*kprof(i_dp,2)*dz_term

    bb(i_dp,1) = 1.0 +temp_term*kprof(i_dp,2)*dz_term                
 
  END DO

  rr(:,:) = q(:,:)-j1*timestep*drmqsedz(:,:)/rho_theta(:,:)
        
      
  !DEPENDS ON: tridiag_all
  CALL tridiag_all(nlev,n_dp,nlevs_arr,aa,bb,cc,rr,q_withflux)

       
! Assemble the approximation for the mse flux, wmse
  ntparmaxp1 = ntparmax+1
  IF (ntparmaxp1 > nlev) THEN
    ntparmaxp1 = nlev-1 
  END IF 
 !On theta levels
 ! Now specifying in boundary layer as well
  DO k=2,ntparmax+1
    DO i_dp=1,n_dp

      wqv(i_dp,k) = -kprof(i_dp,k)*(q_withflux(i_dp,k)-q_withflux(i_dp,k-1))  &
                                      /dr_across_rh(i_dp,k)/rho(i_dp,k)       &
                     + j1*(wq0(i_dp)/w2lcl(i_dp)/qse_valrho(i_dp,1))*         &
                             wwrho(i_dp,k)*qse_valrho(i_dp,k)

    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_QFLUX',zhook_out,zhook_handle)
  RETURN
       
        
END SUBROUTINE dts_qflux
