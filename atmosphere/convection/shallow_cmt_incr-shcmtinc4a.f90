! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates increments to U and V due to shallow convection
!
SUBROUTINE shallow_cmt_incr(np_field,npnts,n_cumulus,nlevs,nterm, &
                            cu_ind,cu_full,nlcl,ntop,uw,vw,phalf, &
                            rho,zlcl,                             &
                            ! Output
                            dubydt,dvbydt)
                            
USE earth_constants_mod, ONLY: g                            
                            
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------
! Description:
!   Calculates increments to U and V due to shallow convection.
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,n_cumulus            & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_ind(nterm)        & ! Indices for terminating points
 ,cu_full(nterm)       & ! Indices of points in output array
 ,nlcl(n_cumulus)      & ! Lifting condensation level
 ,ntop(n_cumulus)        ! Top level of convection


REAL, INTENT(IN)    ::     & 
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)      & ! V-component of stress (m2/s2)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,zlcl(npnts)                ! Height of LCL (m)

REAL, INTENT(OUT) ::      &
  dubydt(np_field,nlevs)  & ! U increment (m/s2)
 ,dvbydt(np_field,nlevs)    ! V increment (m/s2)

! Local variables

INTEGER ::       &
  i,j,k,m,n         ! Loop counters


REAL ::          &
  dudt_bl        & ! U CMT tendency in subcloud-layer (m/s2)
 ,dvdt_bl        & ! V CMT tendency in subcloud-layer (m/s2)
 ,rhodz            ! Density times height difference between half level

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------
IF (lhook) CALL dr_hook('SHALLOW_CMT_INCR',zhook_in,zhook_handle)
!------------------------------------------------------------------------

DO i=1,nterm
  j=cu_ind(i)
  m=cu_full(i)
  dudt_bl=-uw(nlcl(j),j)/zlcl(m)
  dvdt_bl=-vw(nlcl(j),j)/zlcl(m)

! Mixed layer increments (assumed constant in mixed layer)

  DO k=1,nlcl(j)-1
    dubydt(m,k)=dudt_bl/rho(k,j)
    dvbydt(m,k)=dvdt_bl/rho(k,j)
  END DO

! Cloud layer increments

  DO k=nlcl(j),ntop(j)
    rhodz=-(phalf(k+1,j)-phalf(k,j))/g
    dubydt(m,k)=-(uw(k+1,j)-uw(k,j))/rhodz
    dvbydt(m,k)=-(vw(k+1,j)-vw(k,j))/rhodz
  END DO
END DO

IF (lhook) CALL dr_hook('SHALLOW_CMT_INCR',zhook_out,zhook_handle)

RETURN
END SUBROUTINE shallow_cmt_incr
