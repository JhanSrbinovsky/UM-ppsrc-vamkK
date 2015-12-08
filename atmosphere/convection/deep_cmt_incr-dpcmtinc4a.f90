! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculates increments to U and V due to deep convection.

SUBROUTINE deep_cmt_incr(np_field,npnts,nconv,nlevs,nterm,        &
                         nlcl,ntop,cu_term,cu_tend,               &
                         zlcl,phalf,p,rho,                        &
                         uw_base,vw_base,uw,vw,                   &
                         ! Output
                         dubydt,dvbydt)

USE earth_constants_mod, ONLY:g
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Description : 
!  Calculates increments to U and V due to deep convection. 
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
!
!------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,nconv                & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,nlcl(nconv)          & ! Lifting condensation level
 ,ntop(nconv)          & ! Top level of convection
 ,cu_term(nterm)       & ! Indices for terminating points
 ,cu_tend(nterm)         ! Indices of points in output array

REAL, INTENT(IN)    :: & 
  zlcl(npnts)          & ! Height of LCL (m)
 ,phalf(nlevs,nconv)   & ! Pressure on model half levels (Pa)
 ,p(nlevs,nconv)       & ! Pressure on model levels (Pa)
 ,rho(nlevs,nconv)     & ! Density, model uv levels (kg/m3/s)
 ,uw_base(nconv)       & ! cloud base U-component of stress
 ,vw_base(nconv)       & ! cloud base V-component of stress
 ,uw(nlevs,nconv)      & ! U-component of stress (Pam/s2)
 ,vw(nlevs,nconv)        ! V-component of stress (Pam/s2)

REAL, INTENT(OUT) ::      &
  dubydt(np_field,nlevs)  & ! U increment (m/s2)
 ,dvbydt(np_field,nlevs)    ! V increment (m/s2)

! Local variables

INTEGER ::       &
  i,j,k,m,n         ! Loop counters

REAL ::          &
  dudt_bl        & ! U CMT tendency in subcloud-layer (m/s2)
 ,dvdt_bl        & ! V CMT tendency in subcloud-layer (m/s2)
 ,dp               ! Pressure difference between adjacent half levels (Pa)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------
IF (lhook) CALL dr_hook('DEEP_CMT_INCR',zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate U and V wind increments by differentiating stress profile

DO i=1,nterm
  m=cu_term(i)
  n=cu_tend(i)
  j=nlcl(m)

! CMT tendencies in the subcloud layer are constant with height

  dudt_bl=-uw(j,m)/(g*zlcl(n))
  dvdt_bl=-vw(j,m)/(g*zlcl(n))
  DO k=1,j-1
    dubydt(n,k)=dudt_bl
    dvbydt(n,k)=dvdt_bl
  END DO
  DO k=j,ntop(m)+1
    dp=-(phalf(k+1,m)-phalf(k,m))
    dubydt(n,k)=-(uw(k+1,m)-uw(k,m))/dp
    dvbydt(n,k)=-(vw(k+1,m)-vw(k,m))/dp
  END DO
END DO

IF (lhook) CALL dr_hook('DEEP_CMT_INCR',zhook_out,zhook_handle)

RETURN
END SUBROUTINE deep_cmt_incr
