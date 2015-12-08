! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the gradient component of the stress due to shallow convection
!

SUBROUTINE shallow_grad_stress(npnts,n_cumulus,nterm,nlevs,       &
                       cu_ind,nlcl,ntop,mb,wsc,wstr,zcld,plcl,    &
                       ptop,p,phalf,rho,ue,ve,timestep,           &
                      ! Outputs
                       uw,vw)

USE earth_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------
! Description:
!   
! Calculates the gradient component of the stress due to shallow 
! convection.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Total number of points in segment
 ,n_cumulus            & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_ind(nterm)        & ! Indices for terminating points
 ,nlcl(n_cumulus)      & ! Lifting condensation level
 ,ntop(n_cumulus)        ! Top level of convection

REAL, INTENT(IN)    ::     & 
  mb(n_cumulus)            & ! Cloud base mass flux for shallow cu (m/s)
 ,wsc(n_cumulus)           & ! Cloud-layer velocity scale (m/s)
 ,wstr(n_cumulus)          & ! Mixed layer velocity scale (m/s)
 ,zcld(npnts)              & ! cloud layer depth (m)
 ,plcl(n_cumulus)          & ! Pressure at lifting condensation level (Pa)
 ,ptop(n_cumulus)          & ! Pressure at top of cloud layer (Pa)
 ,p(nlevs,n_cumulus)       & ! Pressure on model levels (Pa)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,ue(nlevs,n_cumulus)      & ! U-component of mean wind (m/s)
 ,ve(nlevs,n_cumulus)      & ! V-component of mean wind (m/s)
 ,timestep                   ! Model timestep (s)

REAL, INTENT(OUT) ::       &
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)        ! V-component of stress (m2/s2)

! Local variables

INTEGER ::       &
  i,j,k,m,nlev         ! Loop counters


REAL ::                     &
  w(nlevs,n_cumulus)        & ! Non-dimensional plume vertical velocity
 ,mass(nlevs,n_cumulus)     & ! Mass flux profile (m/s)
 ,visc(nlevs,nterm)         & ! Viscosity profile (m2/s)
 ,a(nlevs)                  & ! Tridiagonal matrix elements
 ,b(nlevs)                  & !
 ,c(nlevs)                  & !
 ,u_t(nlevs)                & ! Current velocity vectors (m/s)
 ,v_t(nlevs)                & !
 ,u_tp1(nlevs)              & ! After timestep velocity vectors
 ,v_tp1(nlevs)              & !
 ,ue_tp1(nlevs,n_cumulus)   & ! After timestep velocity vectors
 ,ve_tp1(nlevs,n_cumulus)   & ! After timestep velocity vectors
 ,w02                       & !
 ,p_depth                   & !
 ,zeta                      & !
 ,entr_sc                   & !
 ,exp_k                     & !
 ,exp_kp1                   & !
 ,dz                        & !
 ,dz12                        !

! Parameters (in future consider putting in a module).

REAL, PARAMETER ::  &
  a_stress=0.3      & !
 ,a_w02=10.24         !


REAL f_w
EXTERNAL f_w

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------
IF (lhook) CALL dr_hook('SHALLOW_GRAD_STRESS',zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate vertical velocity profile in updraughts

DO i=1,nterm
  j=cu_ind(i)
  w02=a_w02*(mb(j)*wstr(j)**2)**0.6667
  p_depth=ptop(j)-plcl(j)
  DO k=nlcl(j),ntop(j)
    zeta=(phalf(k,j)-plcl(j))/p_depth
! DEPENDS ON: f_w
    w(k,i)=SQRT(w02+wsc(j)**2*f_w(zeta))/wsc(j)
  END DO
END DO

! Calculates mass flux profile
! Uses fractional detrainment=1.3  fractional entrainment

DO i=1,nterm
  j=cu_ind(i)
  entr_sc=0.04*wsc(j)/(mb(j)*zcld(j))
  p_depth=ptop(j)-plcl(j)
  mass(nlcl(j),i)=mb(j)
  exp_k=EXP(-(phalf(nlcl(j),j)-plcl(j))/p_depth)
  DO k=nlcl(j),ntop(j)-1
    exp_kp1=EXP(-(phalf(k+1,j)-plcl(j))/p_depth)
    zeta=(1.0-1.3)*entr_sc*p_depth*(exp_kp1-exp_k)/(g*rho(k+1,j))
    mass(k+1,i)=mass(k,i)*EXP(zeta)
    exp_k=exp_kp1
  END DO
END DO

! Calculate the eddy viscosity profile

DO i=1,nterm
  j=cu_ind(i)
  DO k=nlcl(j)+1,ntop(j)+1
    IF(k <  ntop(j)) THEN
      visc(k,i)=a_stress*mass(k,i)*w(k,i)*zcld(j)
    ELSE IF(k == ntop(j)) THEN
      visc(k,i)=0.162*mb(j)*zcld(j)
    ELSE IF(k == ntop(j)+1) THEN
      dz=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      visc(k,i)=0.09*mb(j)*dz
    END IF
  END DO
END DO

! Calculate gradient component of stress

! Use implicit timestepping

DO i=1,nterm
  j=cu_ind(i)
  nlev=0

! Calculate components of tridiagnol matrix and construct vector of
! current timestep wind components

  DO k=nlcl(j),ntop(j)+1
    nlev=nlev+1
    IF(k == nlcl(j)) THEN
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k+1,j)-p(k,j))/(g*(rho(k+1,j)+rho(k,j))/2.0)
      a(nlev)=0.0
      c(nlev)=-visc(k+1,i)*timestep/(dz*dz12)
    ELSE IF(k <= ntop(j)) THEN
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      a(nlev)=-visc(k,i)*timestep/(dz*dz12)
      dz12=-(p(k+1,j)-p(k,j))/(g*(rho(k+1,j)+rho(k,j))/2.0)
      c(nlev)=-visc(k+1,i)*timestep/(dz*dz12)
    ELSE IF(k == ntop(j)+1) THEN
      dz=-(phalf(k+1,j)-phalf(k,j))/(g*rho(k,j))
      dz12=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
      a(nlev)=-visc(k,i)*timestep/(dz*dz12)
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
END DO   ! nterm

! Calculate stress profiles

DO i=1,nterm
  j=cu_ind(i)
  m=nlcl(j)
  uw(m,j)=0.0
  vw(m,j)=0.0
!CDIR NODEP
  DO k=m+1,ntop(j)+1
    dz=-(p(k,j)-p(k-1,j))/(g*(rho(k,j)+rho(k-1,j))/2.0)
    uw(k,j)=-visc(k,i)*(ue_tp1(k,j)-ue_tp1(k-1,j))/dz
    vw(k,j)=-visc(k,i)*(ve_tp1(k,j)-ve_tp1(k-1,j))/dz
  END DO
  uw(ntop(j)+2,j)=0.0
  vw(ntop(j)+2,j)=0.0
END DO

IF (lhook) CALL dr_hook('SHALLOW_GRAD_STRESS',zhook_out,zhook_handle)

RETURN
END SUBROUTINE shallow_grad_stress
