! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the cloud base stress for deep CMT
!

SUBROUTINE deep_ngrad_stress(np_field,npnts,nconv,nterm,nlevs,         &
                             nlcl,ntop,cu_term,cu_comp,cu_tend,        &
                             pstar,uw0,vw0,zlcl,ue,ve,visc,            &
                             mass,p,phalf,rho,timestep,                &
                             ! Input/output
                             uw,vw,                                    &
                             ! Output
                             uw_base,vw_base,uw_dp,vw_dp)

USE cv_stash_flg_mod, ONLY:                                            &
    flg_uw_dp, flg_vw_dp

USE earth_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Description : 
!   To calculate the cloud base stress for deep convection and complete
!   the calculation of the stress profile.
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
 ,nterm                & ! Number of points terminating
 ,nlevs                & ! Number of model levels
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,cu_comp(npnts)       & ! Index array for convecting points
 ,cu_term(nterm)       & ! Indices for terminating points
 ,cu_tend(nterm)         ! Index of points in output array

REAL, INTENT(IN)    :: & 
  pstar(npnts)         & ! surface pressure
 ,uw0(npnts)           & ! Surface shear stress x-component (m2/s2)
 ,vw0(npnts)           & ! Surface shear stress x-component (m2/s2)
 ,zlcl(npnts)          & ! Height of LCL (m)
 ,ue(nlevs,nconv)      & ! Environment U-wind component (m/s)
 ,ve(nlevs,nconv)      & ! Environment V-wind component (m/s)
 ,visc(nlevs,nconv)    & ! Viscosity  
 ,mass(nlevs,nconv)    & ! Updraught mass flux (Pa/s)
 ,p(nlevs,nconv)       & ! Pressure on model levels (hPa)
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (hPa)
 ,rho(nlevs,nconv)     & ! Density, model uv levels (kg/m3/s)
 ,timestep               ! Model timestep (s)

REAL, INTENT(INOUT) ::    &
  uw(nlevs,nconv)         & ! U-component of stress
 ,vw(nlevs,nconv)           ! V-component of stress

REAL, INTENT(OUT) ::      &
  uw_base(nconv)          & ! cloud base U-component of stress
 ,vw_base(nconv)          & ! cloud base V-component of stress
 ,uw_dp(np_field,nlevs)   & ! U-component of stress for stash
 ,vw_dp(np_field,nlevs)     ! V-component of stress for stash


! Local variables

INTEGER ::       &
  i              & ! local array index
 ,k              & ! Level index
 ,j              & ! Indexes points in compressed input arrays
 ,m              & ! 
 ,n                ! Indexes points in uncompressed input arrays

REAL ::            &
  a_0(nterm)       & !
 ,a_u(nterm)       & ! Coefficients neeed for evaluating in cloud wind at cloud
 ,a_v(nterm)       & ! base
 ,omg2_jump(nterm) & ! Cloud base jump in Y component of vorticity
 ,omg1_jump(nterm) & ! Cloud base jump in X component of vorticity
 ,mb(nterm)        & ! Cloud base mass flux (m/s)
 ,dz               & !
 ,dzp1             & !
 ,beta             & !
 ,du               & !
 ,dv               & !
 ,dz1              & !
 ,zeta             & !
 ,a                & !
 ,b                & !
 ,c                & !
 ,t                  !

! In future consider putting parameters in a module
REAL, PARAMETER ::   &
  gamma=1.63         & !
 ,delta=2.0          & !
 ,top_press=15000.0    !

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEEP_NGRAD_STRESS',zhook_in,zhook_handle)

!------------------------------------------------------------------------
! Convert cloud-base mass flux from Pa/s to m/s

DO i=1,nterm
  j=cu_term(i)
  n=cu_comp(i)
  k=nlcl(j)
  mb(i)=mass(k,j)/g
END DO

! Calculate jumps in vorticity components across cloud-base.
! 'Implicit technique' assumes du, dv vary as exp(-t/tau) through timestep
! needed because explicit calculation can lead to instability in du, dv 
! under some circumstances.

DO i=1,nterm
  j=cu_term(i)
  n=cu_comp(i)
  dz=-(p(nlcl(j),j)-phalf(nlcl(j),j))/(g*rho(nlcl(j),j))
  du=(ue(nlcl(j),j)-ue(nlcl(j)-1,j))
  dv=(ve(nlcl(j),j)-ve(nlcl(j)-1,j))
  dz1=-(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
  zeta=-(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/25000.0
  b=(1.0/zlcl(n)-(EXP(-zeta)-1.0)/dz1)
  a=zlcl(n)*mb(i)*b/(delta*dz)
  c=(b*(1.0-gamma/delta)-1.0/zlcl(n))*uw0(n)

  t=-LOG((c*(1.0-EXP(-a*timestep))/a+                              &
          du*(EXP(-a*timestep)-1.0))/                              &
          ((c-a*du)*timestep))/a
  omg2_jump(i)=(c*(1.0-EXP(-a*t))/a+du*EXP(-a*t))/dz
  c=(b*(1.0-gamma/delta)-1/zlcl(n))*vw0(n)
  t=-LOG((c*(1.0-EXP(-a*timestep))/a+                              &
          dv*(EXP(-a*timestep)-1.0))/                              &
          ((c-a*dv)*timestep))/a
  omg1_jump(i)=-(c*(1.0-EXP(-a*t))/a+dv*EXP(-a*t))/dz
END DO

! Calculate cloud base stress components. Note factor of g to convect
! back to Pa/s.

DO i=1,nterm
  j=cu_term(i)
  n=cu_comp(i)
  uw_base(j)=g*(zlcl(n)*(-mb(i)*omg2_jump(i)-                      &
                 gamma*uw0(n)/zlcl(n))/delta+uw0(n))
  vw_base(j)=g*(zlcl(n)*(mb(i)*omg1_jump(i)-                       &
                 gamma*vw0(n)/zlcl(n))/delta+vw0(n))
END DO

! Calculate total stress

! Calculate stress profiles the function beta was again tuned to TOGA-COARE
! CRM simulation

DO i=1,nterm
  m=cu_term(i)
  n=cu_tend(i)
  j=nlcl(m)
  uw(j,m)=uw_base(m)
  vw(j,m)=vw_base(m)

! below cloud base
  DO k=1,nlcl(m)-1
    beta=(phalf(k,m)-pstar(m))/(phalf(j,m)-pstar(m))
    uw(k,m)=uw(k,m)+beta*uw_base(m)
    vw(k,m)=vw(k,m)+beta*vw_base(m)
  END DO

  DO k=j+1,ntop(m)+1
    beta=EXP(((phalf(k,m)-phalf(j,m))/25000.0))
    uw(k,m)=uw(k,m)+beta*uw_base(m)
    vw(k,m)=vw(k,m)+beta*vw_base(m)
  END DO
  uw(ntop(m)+2,m)=0.0
  vw(ntop(m)+2,m)=0.0

! Stash diagnostics (NOTE diagnostics are in m/s for direct comparison
! with shallow convection stresses)

  IF(flg_uw_dp) THEN
    DO k=j,ntop(m)+2
      uw_dp(n,k)=uw(k,m)/g
    END DO
  END IF
  IF(flg_vw_dp) THEN
    DO k=j,ntop(m)+2
      vw_dp(n,k)=vw(k,m)/g
    END DO
  END IF
END DO ! nterm

IF (lhook) CALL dr_hook('DEEP_NGRAD_STRESS',zhook_out,zhook_handle)


RETURN
END SUBROUTINE deep_ngrad_stress
