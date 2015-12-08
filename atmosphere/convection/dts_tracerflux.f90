! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE dts_tracerflux_mod


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

CONTAINS
!
!+  calculate the tracer flux
!
! Subroutine Interface:
SUBROUTINE dts_tracerflux(n_dp,nlev,trlev,ntra,dts_ntpar,ntparmax,         &
                          z_theta,z_rho,dr_across_rh,dr_across_th,rho,     &
                          rho_theta,wwrho,wstar,tracer,zlcl,ztop,timestep, &
                          dtrabydt)

! Modules 

USE dts_fitpars_mod, ONLY :                                                &
  h1, h2

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Calculates the tracer flux in the same way that the moist static
! energy flux was calculated
! Other information: 
! ------------------
! Called from flux calculation section of deep_turb_conv
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
  n_dp                 & ! No. of deep convection points
 ,nlev                 & ! No. of model layers
 ,trlev                & ! No. of tracer levels
 ,ntra                 & ! No. of tracers
 ,dts_ntpar(n_dp)      & ! Top level of initial parcel ascent
 ,ntparmax            

REAL, INTENT(IN) ::       &
  z_theta(n_dp,nlev)      & ! height of theta levels (m)
 ,z_rho(n_dp,nlev)        & ! height of rho levels(m)
 ,dr_across_rh(n_dp,nlev) & ! Thickness of rho layers (m)
 ,dr_across_th(n_dp,nlev) & ! Thickness of theta layers (m)
 ,rho(n_dp,nlev)          & ! Density on rho levels (kg/m3)
 ,rho_theta(n_dp,nlev)    & ! Density on theta levels  (kg/m3)
 ,wwrho(n_dp,nlev)        & ! velocity variance  (m^2/s^-2)
 ,wstar(n_dp)             & ! B-L convective velocity scale (m/s)
 ,tracer(n_dp,trlev,ntra) & ! Tracers on model theta levels (kg/kg)
 ,ztop(n_dp)              & ! top of convection (m)
 ,zlcl(n_dp)              & ! Lifting condensation level (m) 
 ,timestep                  ! timstep  (s)
        
REAL, INTENT(OUT) ::          &
  dtrabydt(n_dp,nlev,ntra)     ! Increment to tracer (only nlev)
                               ! due to convection routine (kg/kg/s)

! ------------------------------------------------------------------------------
! Local variables
! ------------------------------------------------------------------------------

INTEGER ::        & 
  k,i_dp,itra     & ! loop counters
 ,ntparmaxp1      & ! ntparmax +1
 ,nlev_arr(n_dp)    ! number of levels
        
REAL ::                    &
  func(n_dp,nlev)          & ! fixed functional profile
 ,aa(n_dp,nlev)            & ! h(k-1) coefficients
 ,bb(n_dp,nlev)            & ! h(k) coefficients
 ,cc(n_dp,nlev)            & ! h(k+1) coefficients
 ,rr(n_dp,nlev)            & ! r.h.s. of implicit tracer equn
 ,uu(n_dp,nlev)            & ! h on timestep t+1
 ,kprof(n_dp,nlev)         & ! h1 ww
 ,gprof(n_dp,nlev)         & ! h2 ww func
 ,dgdz(n_dp,nlev)          & ! gradient of gprof
 ,rdzthkm1(n_dp,nlev)      & ! 1/dzthkm1
 ,h1all(n_dp,nlev)         &
 ,tracer_withflux(n_dp,trlev,ntra) ! new incremented tracer

! NOT USED at present may require in future ?
!REAL ::                    &
!  dtradz(n_dp,trlev,ntra)  & ! vertical gradient of tracer (on rho levels)
! ,wtracer(n_dp,nlev,ntra)   ! tracer flux 



REAL ::                &
  temp_term            & ! intermediate cal
 ,rtimestep            & ! 1/timestep
 ,temp_array(n_dp)       ! precal divison to save CPU

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------------

IF (lhook) CALL dr_hook('DTS_TRACERFLUX_MOD:DTS_TRACERFLUX',zhook_in,zhook_handle)
  IF(trlev /= nlev) THEN
    write(6,*) 'WARNING: trlev != nlev',trlev,nlev
  END IF

! calculate to save CPU as division more costly than multiplication

  rtimestep = 1.0/timestep

! Set fields to zero
  tracer_withflux(:,:,:) = 0.0
  dtrabydt(:,:,:) = 0.0

!  dtradz(:,:,:)   = 0.0     ! not used in final calculations
!  wtracer(:,:,:)  = 0.0     ! not used in fil calculations

  nlev_arr(:)=nlev

  DO k=1,nlev
    DO i_dp=1,n_dp
      rdzthkm1(i_dp,k) = 1.0/dr_across_rh(i_dp,k)
    END DO
  END DO

        
! Calculate tracer gradient and set up functional form for 'func'
! dtradz is on rho levels  - NOT currently used except for wtracer
!  DO itra=1,ntra
!    DO k=2,trlev
!      DO i_dp=1,n_dp
!        dtradz(i_dp,k,itra)=                                                 &
!                  (tracer(i_dp,k,itra)-tracer(i_dp,k-1,itra))*rdzthkm1(i_dp,k)
!
!      END DO
!    END DO
!  END DO

! Func is on rho levels -- a very simple function, designed to be
! smooth across the boundary layer
! NB set to zero for time being
  func(:,:) = 0.0

! Possible alternative 
!  DO k=1,ntparmax+1
!    DO i_dp=1,n_dp
!      func(i_dp,k) = (1.0/wstar(i_dp)**2)*(1.0-z_rho(i_dp,k)/ztop(i_dp))
!    END DO
!  END DO
       

! The next stage is to solve for:
! dtracer/dt = -d/dz( w'tracer')
! This needs doing implicitly (because of the gradient term) via a
! tridiagonal matrix

  h1all(:,:) = 0.0     ! initialise as use whole array later       

  DO k=1,nlev
    DO i_dp=1,n_dp       
      IF( (ztop(i_dp) > z_theta(i_dp,1)) .AND. (k < dts_ntpar(i_dp)) ) THEN
        h1all(i_dp,k) = h1*(ztop(i_dp)-z_rho(i_dp,k))/ztop(i_dp)
      END IF  
    END DO
  END DO

  kprof(:,:) = h1all(:,:)*wwrho(:,:)*rho(:,:)      ! on rho levels

!  gprof(:,:) = h2*wwrho(:,:)*func(:,:)*rho(:,:)    ! on rho levels 
  gprof(:,:) = 0.0     ! as func set to zero


! Assemble the approximation for the tracer flux, wtracer
! Calculated but not used at present - so now commented out to save CPU
!  DO itra=1,ntra
!    DO k=1,ntparmax
!      DO i_dp=1,n_dp
!        !On theta levels
!        ! Now specifying in boundary layer as well
!        wtracer(i_dp,k,itra) = wwrho(i_dp,k)*(h2*func(i_dp,k) -        &
!                                       h1all(i_dp,k)*dtradz(i_dp,k,itra))     
!      END DO
!    END DO
!  END DO



! dtracerdt_flux = (uu-tracer)/timestep
!array ends set here
!  dgdz(:,1) = 0.0
! Commented out dgdz calculated as gprof = 0.0
!  DO k=1,nlev-1 
!    DO i_dp=1,n_dp
!      dgdz(i_dp,k)=(gprof(i_dp,k+1)-gprof(i_dp,k))/dr_across_th(i_dp,k)
!    END DO
!  END DO

  dgdz(:,:) = 0.0   ! as gprof is zero

! This is intended for above the convecting layer
  aa(:,:) = 0.0
  bb(:,:) = 1.0 
  cc(:,:) = 0.0
     
! now deal with convective layer

  ntparmaxp1 = ntparmax+1 
  IF (ntparmaxp1 >= nlev) THEN
    ntparmaxp1 = nlev-1 
  END IF 
 
  DO k=2,ntparmaxp1
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN
        temp_term=timestep/dr_across_th(i_dp,k)/rho_theta(i_dp,k)

        ! k-1 coef:
        aa(i_dp,k) = -kprof(i_dp,k)*temp_term*rdzthkm1(i_dp,k)          

        ! k+1 coef:
        cc(i_dp,k) = -kprof(i_dp,k+1)*temp_term*rdzthkm1(i_dp,k+1) 
 
        ! k coef:  
        bb(i_dp,k) = 1.0 -aa(i_dp,k) -cc(i_dp,k)
      END IF
    END DO
  END DO
        
! Level k=1
  DO i_dp=1,n_dp

    temp_term = timestep/dr_across_th(i_dp,1)/rho_theta(i_dp,1)

    aa(i_dp,1) = 0.0

    cc(i_dp,1) = - kprof(i_dp,2)*temp_term*rdzthkm1(i_dp,2)

    bb(i_dp,1) = 1.0 -cc(i_dp,1)

!    temp_array(i_dp) = gprof(i_dp,2)*temp_term  ! gprof = 0.0

  END DO

! Start loop over tracers
  DO itra=1,ntra
    !  rr(:,:) = 0.0       ! set all values anyway so not required
    uu(:,:) = 0.0 
    DO k=2,nlev
      DO i_dp=1,n_dp
        IF(k <= dts_ntpar(i_dp)) THEN
          rr(i_dp,k) = tracer(i_dp,k,itra)-dgdz(i_dp,k)*timestep
        ELSE
          rr(i_dp,k) = tracer(i_dp,k,itra)
        END IF
      END DO
    END DO
    ! now do k=1 level, which is treated slightly differently because
    ! rho level 1 isn't strictly a flux level, so need flux to be
    ! applied from the surface:

    DO i_dp=1,n_dp
      ! rr(i_dp,1) = tracer(i_dp,1,itra) - temp_array(i_dp)                  
      ! Altered code as gprof currently set to zero
      rr(i_dp,1) = tracer(i_dp,1,itra)                 
    END DO

    ! Using version with no individual vector length info          
    !DEPENDS ON: tridiag_all_n
    CALL tridiag_all_n(nlev,n_dp,aa,bb,cc,rr,uu)
           
    DO k=1,nlev
      DO i_dp=1,n_dp
          tracer_withflux(i_dp,k,itra) = uu(i_dp,k)
          dtrabydt(i_dp,k,itra) = (tracer_withflux(i_dp,k,itra)       &
                                  -tracer(i_dp,k,itra))*rtimestep

      END DO
    END DO

  END DO ! itra
  IF (lhook) CALL dr_hook('DTS_TRACERFLUX_MOD:DTS_TRACERFLUX',zhook_out,zhook_handle)
  RETURN


END SUBROUTINE dts_tracerflux

END MODULE dts_tracerflux_mod
