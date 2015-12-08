! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! Calculates in cloud gradient stress for deep CMT
!
SUBROUTINE deep_turb_grad_stress(n_dp, nlev, max_cldlev                &
     ,                       nclev                                     &
     ,                       timestep                                  &
     ,                       zcld                                      &
     ,                       mf_cld, w_up, k_func                      &
     ,                       u,v                                       &
     ,                       r2rho,r2rho_th,dr_across_rh,dr_across_th  &
  ! output arguements
     ,                       uw_cld,vw_cld)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


! ------------------------------------------------------------------------------
! Description:
!   This routine calculates the in-cloud gradient part of the deep CMT stress
!   using turbulence ideas. This version is designed for use with the mass flux
!   convection scheme. Note at this stage this version is different from the
!   shallow turbulence version in its assumptions at cloud top.
!
!  Solves the equation
!      du/dt  = -d(uw)/dz for the gradient part of uw.
!
!  The above equation is solved implicitly
!
!  The gradient part of uw = - mf(wup/wcld)*K(eta)du/dz
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  n_dp                 & ! No. of deep convection points
, nlev                 & ! No. of model layers
, max_cldlev           & ! Maximum number of cloud levels
, nclev(n_dp)            ! number of cloud levels

REAL, INTENT(IN) ::    &
  timestep               ! model timestep (s)

REAL, INTENT(IN) ::       &
  zcld(n_dp)              & ! cloud layer depth (m) for CMT
, mf_cld(n_dp,nlev)       & ! mass flux in (th lev) (m/s)
, w_up(n_dp,nlev)         & ! wup/wcld
, k_func(n_dp,nlev)       & ! k function (dependent on eta)
, u(n_dp,nlev)            & ! U-component of mean wind (m/s)
, v(n_dp,nlev)            & ! V-component of mean wind (m/s)
, r2rho(n_dp,nlev)        & ! r2*rho on rho levels (kg/m)
, r2rho_th(n_dp,nlev)     & ! r2*rho on theta levels (kg/m)
, dr_across_rh(n_dp,nlev) & ! thickness of rho layers (m)
, dr_across_th(n_dp,nlev)   ! thickness of theta layers (m)



! Arguments with intent out:

REAL, INTENT(OUT) ::    &
  uw_cld(n_dp,nlev)     & ! U-component of stress (m2/s2)
, vw_cld(n_dp,nlev)       ! V-component of stress (m2/s2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: i,k    ! counters

REAL ::                       &
  visc(n_dp,max_cldlev+1)     & ! viscosity profile (m2/s)
, a(n_dp,max_cldlev+1)        & ! tridiagonal matrix elements
, b(n_dp,max_cldlev+1)        &
, c(n_dp,max_cldlev+1)        &
, ue_tp1(n_dp,max_cldlev+1)   & ! after timestep velocity vectors
, ve_tp1(n_dp,max_cldlev+1)   &
, dz,dz12

INTEGER ::           &
  max_cldlev1          ! max cloud levels plus 1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DEEP_TURB_GRAD_STRESS',zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! Calculate the eddy viscosity profile

!        K'= K(z/zcld)*mf*zcld*wup/wcld

!   on theta levels - ie uw stress levels
! Note no special conditions for values at cloud top unlike shallow
! code.
!-----------------------------------------------------------------------

max_cldlev1 = max_cldlev+1

!CDIR NOUNROLL
DO k=1,max_cldlev1  ! from level above cloud base
  DO i=1,n_dp

    IF (k <=  nclev(i)) THEN         ! in cloud

      visc(i,k)=k_func(i,k)*mf_cld(i,k)*w_up(i,k)*zcld(i)

    ELSE
      visc(i,k)=0.0
    END IF
  END DO
END DO


! Calculate GRADIENT component of stress. Use implicit timestepping

k=1
  DO i=1,n_dp
    dz   = dr_across_rh(i,k)*r2rho(i,k)
    dz12 = dr_across_th(i,k)
    a(i,k) = 0.0
    c(i,k) = -visc(i,k)*r2rho_th(i,k)*timestep/(dz*dz12)
    b(i,k) = 1.0 - a(i,k) - c(i,k)
  END DO

!CDIR NOUNROLL
DO k=2,max_cldlev1
  DO i=1,n_dp
    dz = dr_across_rh(i,k)*r2rho(i,k)

    IF (k < nclev(i)) THEN

      dz12 = dr_across_th(i,k-1)
      a(i,k) = -visc(i,k-1)*r2rho_th(i,k-1)*timestep/(dz*dz12)

      dz12 = dr_across_th(i,k)
      c(i,k) = -visc(i,k)*r2rho_th(i,k)*timestep/(dz*dz12)

    ELSE IF (k == nclev(i)) THEN

      dz12 = dr_across_th(i,k-1)
      a(i,k) = -visc(i,k-1)*r2rho_th(i,k-1)*timestep/(dz*dz12)
      c(i,k) = 0.0

    ELSE    ! elements not required in calculation (zero)

      c(i,k) = 0.0
      a(i,k) = 0.0

    END IF

    b(i,k) = 1.0 - a(i,k) - c(i,k)

  END DO
END DO


! Calculate NEW timestep wind conponents using tridiagonal matrix solver

! DEPENDS ON: tridiag_all
CALL tridiag_all(max_cldlev1,n_dp,nclev,a,b,c,u,ue_tp1)
! DEPENDS ON: tridiag_all
CALL tridiag_all(max_cldlev1,n_dp,nclev,a,b,c,v,ve_tp1)



! Calculate stress profile -Kdu/dz from latest u/v values


!CDIR NOUNROLL
DO k=1,max_cldlev1
  DO i=1,n_dp

    IF (k < nclev(i)) THEN   ! in cloud

      dz = dr_across_th(i,k)
      uw_cld(i,k)=-visc(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))/dz
      vw_cld(i,k)=-visc(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))/dz

    ELSE        ! not in cloud  set to zero

      uw_cld(i,k) = 0.0
      vw_cld(i,k) = 0.0

    END IF

  END DO
END DO

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DEEP_TURB_GRAD_STRESS',zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE deep_turb_grad_stress
