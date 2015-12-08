! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  calculate the mse flux
!

SUBROUTINE dts_mseflux(n_dp,nlev,dts_ntpar,ntparmax,                       &
           z_theta,z_rho,rho,rho_theta,dr_across_rh,dr_across_th,          &
           wwrho,wstar,theta,q,mse,zlcl,klcl,ztop,wq0,wth0,timestep,       &
           wmse,mse_withflux,h1all)

USE atmos_constants_mod, ONLY: cp

USE dts_fitpars_mod, ONLY:                                                 &
  h1, h2 

USE water_constants_mod, ONLY: lc
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! calculate the mse flux
! called from FLUX CALCULATION section of deep_turb_conv
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! What this subroutine does:
! --------------------------
! calculates the moist static energy flux
! see report
!
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  nlev                    & ! No. of model layers
 ,n_dp                    & ! No. of deep convection points
 ,dts_ntpar(n_dp)         & ! Top level of initial parcel ascent
 ,ntparmax                &
 ,klcl(n_dp)              

REAL, INTENT(IN)    ::   &
  z_theta(n_dp,nlev)     & ! height of theta levels (m)
 ,z_rho(n_dp,nlev)       & ! height of rho levels(m)
 ,ztop(n_dp)             & ! top of convection (m)
 ,zlcl(n_dp)             & ! Lifting condensation level (m) 
 ,theta(n_dp,nlev)       & ! potential temperature (K)
 ,q(n_dp,nlev)           & ! water vapour  (kg/kg)
 ,mse(n_dp,nlev)         & ! moist static energy   (K)
 ,wwrho(n_dp,nlev)       & ! velocity variance   (m^2/s^-2)
 ,wth0(n_dp)             & ! surface theta flux
 ,wq0(n_dp)              & ! surface q flux
 ,wstar(n_dp)            & ! B-L convective velocity scale (m/s)
 ,timestep               & ! timstep  (s)
 ,rho(n_dp,nlev)         & ! Density on rho levels (kg/m3) 
 ,rho_theta(n_dp,nlev)   & ! Density on theta levels (kg/m3) 
 ,dr_across_rh(n_dp,nlev)& ! rho layer thickness (m)
 ,dr_across_th(n_dp,nlev)  ! theta layer thickness (m)
        
! Arguments with intent OUT:
REAL, INTENT(OUT)    ::    &
  mse_withflux(n_dp,nlev)  & ! new incremented mse
 ,wmse(n_dp,nlev)          & ! mse flux -- estimated in this routine 
 ,h1all(n_dp,nlev)           !


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  k,i_dp ! loop counters

INTEGER ::         &
  nlevs_arr(n_dp)             ! needed for tridiag   

REAL ::                &
  dhdz(n_dp,nlev)      & ! vertical gradient of mse (on rho levels)
 ,wh0(n_dp)            & ! mse flux at surface
 ,func(n_dp,nlev)      & ! fixed functional profile
 ,aa(n_dp,nlev)        & ! h(k-1) coefficients
 ,bb(n_dp,nlev)        & ! h(k) coefficients
 ,cc(n_dp,nlev)        & ! h(k+1) coefficients
 ,rr(n_dp,nlev)        & ! r.h.s. of implicit mse equn
 ,kprof(n_dp,nlev)     & ! h1 ww
 ,gprof(n_dp,nlev)     & ! h2 ww func
 ,dkdz (n_dp,nlev)     & ! gradient of kprof  (NOT USED)
 ,dgdz(n_dp,nlev)        ! gradient of gprof

REAL ::            &
  temp_term

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
       

!-----------------------------------------------------------------------

! Set fields to zero

  IF (lhook) CALL dr_hook('DTS_MSEFLUX',zhook_in,zhook_handle)

  mse_withflux(:,:) = 0.0
  dhdz(:,:) = 0.0
  wmse(:,:) = 0.0
       
  nlevs_arr(:) = nlev

! Calculate surface mse flux
  DO i_dp=1,n_dp
    wh0(i_dp) = wth0(i_dp) + wq0(i_dp)*lc/cp
  END DO

! Calculate mse gradient and set up functional form for 'func'
  DO k=1,ntparmax+1
    DO i_dp=1,n_dp
      ! dhdz is on rho levels
      IF(k > 1) THEN
        dhdz(i_dp,k)= (mse(i_dp,k)-mse(i_dp,k-1))/dr_across_rh(i_dp,k)
      END IF

      ! Func is on rho levels -- a very simple function, designed to be
      ! smooth across the boundary layer
              
      func(i_dp,k)=(wh0(i_dp)/wstar(i_dp)**2)*(1.0-z_rho(i_dp,k)/ztop(i_dp))
    END DO
  END DO


! The next stage is to solve for:
! dh/dt = -d/dz( w'h')
! This needs doing implicitly (because of the gradient term) via a 
! tridiagonal matrix
! The following profiles are required
!        kprof(:,:) = h1*wwrho(:,:)*rho(:,:)      ! on rho levels
!        gprof(:,:) = h2*wwrho(:,:)*func(:,:)*rho(:,:) ! on rho levels
        
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBbb
!NB brutal fudge to ensure that negative flux at top has only a small
! impact, will also need to do similar for wthv_Smth
! This really seems to have a positive impact on the run though...

  h1all(:,:) = h1
  DO k=1,nlev
    DO i_dp=1,n_dp
      h1all(i_dp,k) = h1*(ztop(i_dp)-z_rho(i_dp,k))/ztop(i_dp)
             
! crude way to ensure height well above lcl nb could be better
!     IF(wmse(i_dp,k) < 0. .AND. z_rho(i_dp,k) > ztop(i_dp)*0.5) THEN
!       h1all(i_dp,k) = h1/10.
!     END IF
    END DO
  END DO

  kprof(:,:) = h1all(:,:)*wwrho(:,:)*rho(:,:)      ! on rho levels
  gprof(:,:) = h2*wwrho(:,:)*func(:,:)*rho(:,:)    ! on rho levels

!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

! set kprof to zero above convective top
  DO i_dp=1,n_dp
    kprof(i_dp,dts_ntpar(i_dp):nlev) = 0.0
    gprof(i_dp,dts_ntpar(i_dp):nlev) = 0.0
  END DO

! Assemble the approximation for the mse flux, wmse on theta levels
  DO k=1,ntparmax+1
    DO i_dp=1,n_dp

! Now specifying in boundary layer as well
!     IF(k >= klcl(i_dp)) THEN 
      wmse(i_dp,k) = wwrho(i_dp,k)*(h2*func(i_dp,k) -                      &
                                              h1all(i_dp,k)*dhdz(i_dp,k))     
!     END IF
    END DO
  END DO



! dmsedt_flux = (uu-mse)/timestep

!array ends set here
!  dkdz(:,1) = 0.
!  dgdz(:,1) = 0.  !set anyway so not needed
   dgdz(:,nlev) = 0.0 



!    DO i_dp=1,n_dp
!
! Find vertical gradients of kprof, gprof, putting onto theta levels 
!           dkdz(i_dp,1:nlev-1)=(kprof(i_dp,2:nlev)-kprof(i_dp,1:nlev-1))/&
!                               (z_rho(i_dp,2:nlev)-z_rho(i_dp,1:nlev-1))
!    END DO

! Find vertical gradients of gprof, putting onto theta levels 
! Note level 1 not strictly correct as using layer thickness but not used later.
  DO k=1,nlev-1
    DO i_dp=1,n_dp

      dgdz(i_dp,k)=(gprof(i_dp,k+1)-gprof(i_dp,k))/dr_across_th(i_dp,k)

    END DO
  END DO

! This is intended for above the convecting layer
  aa(:,:) = 0.0
  bb(:,:) = 1.0 
  cc(:,:) = 0.0
  rr(:,:) = 0.0
! now deal with convective layer
  DO k=2,ntparmax+1
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp) .AND. k < nlev) THEN 
        temp_term = timestep/dr_across_th(i_dp,k)/rho_theta(i_dp,k)

        ! k-1 coef:
        aa(i_dp,k) = -kprof(i_dp,k)*temp_term/dr_across_rh(i_dp,k)
             
        ! k+1 coef:
        cc(i_dp,k) = -kprof(i_dp,k+1)*temp_term/dr_across_rh(i_dp,k+1)
        ! k coef:  
        bb(i_dp,k) = 1.0 - cc(i_dp,k) - aa(i_dp,k)

      END IF
    END DO
  END DO
        
! now do k=1 level

  DO i_dp=1,n_dp
    aa(i_dp,1) = 0.0

    cc(i_dp,1) = - timestep*kprof(i_dp,2)/dr_across_th(i_dp,1)/         &
                                dr_across_rh(i_dp,2)/rho_theta(i_dp,1)
    bb(i_dp,1) = 1.0 - cc(i_dp,1)

  END DO

  DO k=2,nlev
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN
        rr(i_dp,k) = mse(i_dp,k)-dgdz(i_dp,k)*timestep               
      ELSE
        rr(i_dp,k) = mse(i_dp,k)
      END IF
    END DO
  END DO

! Now do k=1 level, which is treated slightly differently because
! rho level 1 isn't strictly a flux level, so need flux to be
! applied from the surface:
  DO i_dp=1,n_dp
    rr(i_dp,1) = mse(i_dp,1)-                                              &
                   gprof(i_dp,2)*timestep/z_rho(i_dp,2)/rho_theta(i_dp,1)
  END DO

  !DEPENDS ON: tridiag_all
  CALL tridiag_all(nlev,n_dp,nlevs_arr,aa,bb,cc,rr,mse_withflux)
  IF (lhook) CALL dr_hook('DTS_MSEFLUX',zhook_out,zhook_handle)
  RETURN


END SUBROUTINE dts_mseflux
