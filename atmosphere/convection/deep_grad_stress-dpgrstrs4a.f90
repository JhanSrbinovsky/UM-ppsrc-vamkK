! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the gradient component of the stress profile for deep CMT
!
SUBROUTINE deep_grad_stress(np_field,npnts,nconv,nlevs,nlcl,ntop,      &
                            nterm,cu_term,cu_tend,                     &
                            ue,ve,visc,phalf,p,rho,timestep,           &
                            ! Output
                            uw,vw)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Description : 
!   To calculate the gradient component of the stress profile for deep
!   convection. Calculation an be done explicitly or implicitly. 
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,nterm                & ! Number of points terminating
 ,cu_term(nterm)       & ! Indices for terminating points
 ,cu_tend(nterm)         ! Index of points in output array

REAL, INTENT(IN)    :: & 
  ue(nlevs,nconv)      & ! Environment U-wind component (m/s)
 ,ve(nlevs,nconv)      & ! Environment V-wind component (m/s)
 ,visc(nlevs,nconv)    & ! Viscosity  
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (hPa)
 ,p(nlevs,nconv)       & ! Pressure on model levels (hPa)
 ,timestep             & ! Model timestep (s)
 ,rho(nlevs,nconv)       ! Density, model uv levels (kg/m3/s)

REAL, INTENT(OUT) ::   &
  uw(nlevs,nconv)      & ! U-component of viscous stress
 ,vw(nlevs,nconv)        ! V-component of viscous stress

INTEGER       :: &
  i,j,k,m,n      &  ! loop counters
  ,nlev             ! Number of levels

REAL ::           &
  a(nlevs)        & ! Implicit solver variables
 ,b(nlevs)        & !
 ,c(nlevs)        & !
 ,u_t(nlevs)      & ! Current U wind
 ,u_tp1(nlevs)    & ! U Wind at T+1
 ,v_t(nlevs)      & ! Current V wind
 ,v_tp1(nlevs)      ! V wind at T+1

REAL ::                 &
  up(nlevs,nconv)       & ! In cloud U-wind component (m/s) 
 ,vp(nlevs,nconv)       & ! In cloud V-wind component (m/s) 
 ,ue_tp1(nlevs,nconv)   & ! Implicitly updated U wind
 ,ve_tp1(nlevs,nconv)     ! Implicitly updated V wind


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEEP_GRAD_STRESS',zhook_in,zhook_handle)

!------------------------------------------------------------------------

! Use implicit timestepping

DO i=1,nterm
  j=cu_term(i)
  nlev=0

! Calculate components of tridiagonal matirx and construct a vector of
! current timestep wind components

  DO k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    IF(k == nlcl(j)) THEN
      a(nlev)=0.0
      c(nlev)=-visc(k+1,j)*timestep/                               &
              ((p(k+1,j)-p(k,j))*(phalf(k+1,j)-phalf(k,j)))
    ELSE IF(k <= ntop(j)) THEN
      a(nlev)=-visc(k,j)*timestep/                                 &
              ((p(k,j)-p(k-1,j))*(phalf(k+1,j)-phalf(k,j)))
      c(nlev)=-visc(k+1,j)*timestep/                               &
              ((p(k+1,j)-p(k,j))*(phalf(k+1,j)-phalf(k,j)))
    ELSE IF(k == ntop(j)+1) THEN
      a(nlev)=-visc(k,j)*timestep/                                 &
              ((p(k,j)-p(k-1,j))*(phalf(k+1,j)-phalf(k,j)))
      c(nlev)=0.0
    END IF
    b(nlev)=1.0-a(nlev)-c(nlev)
    u_t(nlev)=ue(k,j)
    v_t(nlev)=ve(k,j)
  END DO

! Calculate new timestep wind components using tridiag

! DEPENDS ON: tridiag
  CALL tridiag(a,b,c,u_t,u_tp1,nlev)
! DEPENDS ON: tridiag
  CALL tridiag(a,b,c,v_t,v_tp1,nlev)

! Store updated wind components for later

  nlev=0
  DO k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    ue_tp1(k,j)=u_tp1(nlev)
    ve_tp1(k,j)=v_tp1(nlev)
  END DO
END DO      ! nterm

! Calculate stress profiles

DO i=1,nterm
  m=cu_term(i)
  j=nlcl(m)
  uw(j,m)=0.0
  vw(j,m)=0.0
  DO k=j+1,ntop(m)+1
    uw(k,m)=-visc(k,m)*(ue_tp1(k,m)-ue_tp1(k-1,m))/(p(k-1,m)-p(k,m))
    vw(k,m)=-visc(k,m)*(ve_tp1(k,m)-ve_tp1(k-1,m))/(p(k-1,m)-p(k,m))
  END DO
  uw(ntop(m)+2,m)=0.0
  vw(ntop(m)+2,m)=0.0
END DO

IF (lhook) CALL dr_hook('DEEP_GRAD_STRESS',zhook_out,zhook_handle)

RETURN
END SUBROUTINE deep_grad_stress
