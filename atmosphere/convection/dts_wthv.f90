! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! applies a gradient term to the thetav profile at the freezing level
!
! Subroutine Interface:
SUBROUTINE dts_wthv(n_dp,nlev,ntparmax,dts_ntpar,                          &
                    z_theta,z_rho,dr_across_rh,dr_across_th,rho,rho_theta, &
                    ztop,zlcl,thetav,dthvdz,dthvdz_m,wwrho,timestep,       & 
                    wthv,dthvdt_flux)

! Modules 
USE dts_fitpars_mod, ONLY :                                                &
  h1, h4

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!  Applies a gradient term to the thetav profile at the freezing level
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
  n_dp                 & ! Number of deep points
 ,nlev                 & ! number of model levels
 ,dts_ntpar(n_dp)      & ! top of cloud
 ,ntparmax               ! maximum cloud top level
           
           
REAL, INTENT(IN) ::       &
  wwrho(n_dp,nlev)        & ! Velocity variance on rho levels (m2/s2?)
 ,z_rho(n_dp,nlev)        & ! height of model rho levels (m)
 ,z_theta(n_dp,nlev)      & ! height of model theta levels (m)
 ,dr_across_rh(n_dp,nlev) & ! thickness of rho layers (m)
 ,dr_across_th(n_dp,nlev) & ! thickness of theta layers (m)
 ,zlcl(n_dp)              & ! Lifting condensation level (m)
 ,ztop(n_dp)              & ! cloud top (m)
 ,thetav(n_dp,nlev)       & ! thetav (K)
 ,timestep                & ! convection timestep (s)
 ,dthvdz_m(n_dp,nlev)     & ! d(thv)/dz   (K/m)
 ,dthvdz(n_dp,nlev)       & ! d(thv)/dz   (K/m)
 ,rho(n_dp,nlev)          & ! density on model rho levels (kg/m3)
 ,rho_theta(n_dp,nlev)      ! density on model theta levels (kg/m3)       

REAL, INTENT(OUT) ::      &
  wthv(n_dp,nlev)         & ! wthv flux
 ,dthvdt_flux(n_dp,nlev)    ! smoothing increment to theta_v 

! Local variables

INTEGER ::        & 
  i_dp,k,ival           ! loop counters

INTEGER ::        & 
  nlevs_arr(n_dp)   ! Number of levels

REAL ::           &
  aa(n_dp,nlev)   &
 ,bb(n_dp,nlev)   &
 ,cc(n_dp,nlev)   &
 ,rr(n_dp,nlev)   &
 ,uu(n_dp,nlev)   &
 ,krho(n_dp,nlev) & ! diffusion coefficient
 ,gg(n_dp,nlev)     !smoothing factor 

REAL ::  &
  temp_term      ! 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!============================================================================
! will next want to relax back to these quantities

IF (lhook) CALL dr_hook('DTS_WTHV',zhook_in,zhook_handle)
  krho(:,:) = 0.0
  wthv(:,:) = 0.0
  dthvdt_flux(:,:) = 0.0    

! Need number of level to solve for in each column of the tridiagonal solver
  nlevs_arr(:) = nlev

!  gg(:,:) = 0.0    set to 1.0 below


! diffusion coefficient
  DO k=1,nlev
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN 
                 
        ! term at high levels
        krho(i_dp,k) = h1*wwrho(i_dp,k)*rho(i_dp,k)*                          &
                                           (1.0-z_rho(i_dp,k)/ztop(i_dp))
      END IF
      ! set to have same K profile as mse flux...

      ! a smoothing factor gg to go in front of dthvdz_m
      gg(i_dp,k) = 1.0
      IF(z_rho(i_dp,k) < zlcl(i_dp)) THEN 
        gg(i_dp,k) = 0.0
      END IF
      IF(z_rho(i_dp,k) < 2*zlcl(i_dp) .AND.                                   &
         z_rho(i_dp,k) > zlcl(i_dp)) THEN 

        gg(i_dp,k) = (z_rho(i_dp,k)-zlcl(i_dp))/zlcl(i_dp)

      END IF
    END DO   ! i_dp
  END DO     ! k
        
! this is for above the convective layer
  aa(:,:) = 0.0
  bb(:,:) = 1.0 
  cc(:,:) = 0.0
  uu(:,:) = 0.0
  rr(:,:) = thetav(:,:) ! updated in next part for convectivelayer

       
  DO k=2,ntparmax+1
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp) .AND. k < nlev) THEN 

        temp_term = timestep/(rho_theta(i_dp,k)*dr_across_th(i_dp,k))

        ! k+1 coefficient
        cc(i_dp,k) = -krho(i_dp,k+1)*temp_term/dr_across_rh(i_dp,k+1)

        ! k-1 coefficient
        aa(i_dp,k) = -krho(i_dp,k)*temp_term/dr_across_rh(i_dp,k)
         
        ! k coefficient
        bb(i_dp,k) = 1.0 -cc(i_dp,k) -aa(i_dp,k)

    
        IF(z_theta(i_dp,k) > zlcl(i_dp)) THEN

          rr(i_dp,k) = thetav(i_dp,k)- temp_term*                              &
                         (krho(i_dp,k+1)*(gg(i_dp,k+1)*dthvdz_m(i_dp,k+1)+h4)  &
                        - krho(i_dp,k)  *(gg(i_dp,k)  *dthvdz_m(i_dp,k)  +h4))

        ELSE
          rr(i_dp,k) = thetav(i_dp,k)
        END IF ! z >< zlcl

      END IF ! k <= ntpar
    END DO ! i_dp
  END DO ! k

! now deal with k=1 case:
  DO i_dp=1,n_dp

    temp_term = timestep*krho(i_dp,2)/(z_rho(i_dp,2)                        &
                       *dr_across_rh(i_dp,2)*rho_theta(i_dp,1))

    aa(i_dp,1) = 0.0
    bb(i_dp,1) = 1.0 + temp_term 

    cc(i_dp,1) = - temp_term
    rr(i_dp,1) = thetav(i_dp,1)
  END DO

  !DEPENDS ON: tridiag_all
  CALL tridiag_all(nlev,n_dp,nlevs_arr,aa,bb,cc,rr,uu)

  DO k=1,ntparmax
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN 
        dthvdt_flux(i_dp,k) = (uu(i_dp,k)-thetav(i_dp,k))/timestep 
      END IF
    END DO
  END DO

! assemble what the contribution to the flux is:
  DO k=1,ntparmax
    DO i_dp=1,n_dp           
      wthv(i_dp,k)=krho(i_dp,k)*(h4-(dthvdz(i_dp,k)                         &
                                  -gg(i_dp,k)*dthvdz_m(i_dp,k)))/rho(i_dp,k)
    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_WTHV',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_wthv
