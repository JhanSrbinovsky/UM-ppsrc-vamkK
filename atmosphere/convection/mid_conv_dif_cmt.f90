! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
SUBROUTINE mid_conv_dif_cmt (npnts, nmid, nlev, nmax_layer,             &
                             index_mid, l_mid_all,                      &
                             timestep,                                  &
                             u, v, r_theta, r_rho,                      &
                             z_theta, z_rho, rho, rho_theta,            &
                             p_layer_boundaries,                        &
                             mass_flux, entrain_up,                     &
                             dubydt, dvbydt, uw, vw )

USE cv_param_mod, ONLY: cpress_term
USE earth_constants_mod, ONLY: g
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Diffusive mid-level CMT scheme.
!   Calculates increments to U and V wind components due to mid-level
!   convective momentum transport.
!   Note several lots of mid-level convection can occur in the same column
!   making it slightly more difficult to process the information.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
! ------------------------------------------------------------------------------

! Arguments:
! arguments with INTENT in

INTEGER, INTENT(IN) ::   &
  npnts                  &  ! No. of points
 ,nmid                   &  ! No. of points with mid level convection
 ,nlev                   &  ! No. of model levels
 ,nmax_layer                ! Maximum number of mid-level layers

INTEGER, INTENT(IN) :: &
  index_mid(nmid)        ! index of mid points in full arrays

LOGICAL, INTENT(IN) :: &
  l_mid_all(npnts)       ! true if column has mid-level convection

REAL, INTENT(IN) ::    &
  timestep               ! timestep (s)

REAL, INTENT(IN) ::                &
  u(npnts,nlev)                    & ! u component of wind   (m/s)
 ,v(npnts,nlev)                    & ! v component of wind   (m/s)
 ,r_theta(npnts,0:nlev)            & ! radius of theta levels (m)
 ,r_rho(npnts,nlev)                & ! radius of rho levels (m)
 ,z_theta(npnts,nlev)              & ! height of theta levels (m)
 ,z_rho(npnts,nlev)                & ! height of rho levels (m)
 ,p_layer_boundaries(npnts,0:nlev) & ! pressure on layer boundaries (Pa)
 ,rho(npnts,nlev)                  & ! density on rho levels (kg/m3)
 ,rho_theta(npnts,nlev)            & ! density on theta levels (kg/m3)
 ,mass_flux(npnts,nlev)            & ! mass flux (pa/s)
 ,entrain_up(npnts,nlev)             ! entrained fraction * mass_flux (Pa/s)

! arguments with INTENT out

REAL, INTENT(OUT) ::    &
  dubydt(npnts,nlev+1)  &   ! increment to u component of wind   (m/s/s)
 ,dvbydt(npnts,nlev+1)  &   ! increment to v component of wind   (m/s/s)
 ,uw(npnts,nlev)        &   ! uw Stress profile (N/m2)
 ,vw(npnts,nlev)            ! vw Stress profile (N/m2)

! Local declarations:

INTEGER ::           &
  i,j,k, kk          & ! loop counters
 ,jmid               & ! loop counter
 ,jbase                ! pointer for cloud base

INTEGER ::                  &
  nc_layers(nmid)           & ! number of mid- level layers
 ,nlev_mid(nmid)            & ! total number of levels with increments
 ,ic_base(nmid,nmax_layer)  & ! locations of cloud bases
 ,ic_top(nmid,nmax_layer)     ! locations of cloud tops


LOGICAL ::           &
  convecting(nmid)       ! Indicates level convecting

REAL ::                &
  cfact                & ! (1. - C)
 ,dz                   & ! layer depth
 ,dz_rho               & ! layer depth
 ,dz_theta             & ! layer depth
 ,recip_epsilon        & ! 1/(entrainment rate)
 ,recip_g              & ! 1/g
 ,length_scale         & ! depth in cloud or 1/entrain_rate
 ,grid_cor_fac           ! Correction factor to fluxes due passing through
                         ! a gridbox area increasing with height

REAL ::                      &
  massf(nmid,nlev)           & ! mass flux for mid columns in kg/m2/s
 ,entrain_rate(nmid,nlev)    & ! entrainment rate (/m)
 ,terma(nmid,nlev)           & ! term A
 ,termb(nmid,nlev)           & ! term b
 ,termc(nmid,nlev)           & ! term c
 ,k_dif(nmid,nlev)           & ! Diffusion coefficient/dz
 ,ue_tp1(nmid,nlev)          & ! U at T+1
 ,ve_tp1(nmid,nlev)          & ! V at T+1
 ,u_com(nmid,nlev)           & ! U compressed to just levels altered
 ,v_com(nmid,nlev)           & ! V compressed to just levels altered
 ,z_cbase(nmid,nmax_layer)     ! cloud base for each layer (m)

REAL ::                &
  dubydt2(npnts,nlev)  &   ! increment to u component of wind   (m/s/s)
 ,dvbydt2(npnts,nlev)      ! increment to v component of wind   (m/s/s)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------------------
IF (lhook) CALL dr_hook('MID_CONV_DIF_CMT',zhook_in,zhook_handle)
!-------------------------------------------------------------------------------
! 1.0 Initialise arrays particularly those output.

!    1-C  term

cfact = 1.0 -cpress_term
!
recip_g = 1.0/g

DO k=1,nlev
  DO i=1,npnts
    uw(i,k) = 0.0
    vw(i,k) = 0.0
    ue_tp1(i,k) = 0.0
    ve_tp1(i,k) = 0.0
  END DO
END DO
DO k=1,nlev+1
  DO i=1,npnts
    dubydt(i,k) = 0.0
    dubydt(i,k) = 0.0
  END DO
END DO

!-------------------------------------------------------------------------------
! Extract mass flux and entrainment rates for mid-level convection.
! Need to convert incoming mass flux from Pa/s to kg/m2/s
! Need to convert incoming entrain_up = entrain_rate * dz * mass_flx
! Note mid-level convection only starts from level 2.

DO k=2,nlev
  DO i=1,nmid
    jmid =index_mid(i)
    massf(i,k) = mass_flux(jmid,k)*recip_g

    ! Check mass flux is greater than and very small number.
    IF (mass_flux(jmid,k) >= 1.0e-10) THEN
      dz = z_rho(jmid,k+1)-z_rho(jmid,k)
      entrain_rate(i,k) = entrain_up(jmid,k)/(mass_flux(jmid,k)*dz)

    ELSE
      entrain_rate(i,k) = 0.0
      massf(i,k) = 0.0
    END IF              ! test on mass flux >~ 0
  END DO                ! Mid points loop
END DO                  ! level loop


!-------------------------------------------------------------------------------
! Work out which levels are convecting and how many mid-level layers
! For each layer work out cloud base and cloud top.
!-------------------------------------------------------------------------------

DO i=1,nmid
  convecting(i) = .FALSE.
  nc_layers(i) = 0
END DO

DO k=2,nlev
  DO i=1,nmid
    jmid =index_mid(i)

    IF (convecting(i)) THEN ! Convection so check for top of convection
                            ! i.e. no mass flux

      IF(massf(i,k) == 0.0) THEN
        ic_top(i,nc_layers(i)) = k-1
        convecting(i) = .FALSE.
      END IF

    ELSE      ! Not convecting check for base i.e. a mass flux

      IF(massf(i,k) > 0.0) THEN
        nc_layers(i) = nc_layers(i) + 1
        ic_base(i,nc_layers(i)) = k
        z_cbase(i,nc_layers(i)) = z_rho(jmid,k)
        convecting(i) = .TRUE.
      END IF

    END IF              ! test on convecting
  END DO                ! Mid points loop
END DO                  ! level loop

!-------------------------------------------------------------------------------
! Evaluate terms for tridiagnol matrix
! Only winds for levels 1 - nlev therefore can only calculate wind stress
! for level 2 (no convection from level 1) to nlev-1.
!-------------------------------------------------------------------------------

DO i=1,nmid
  nlev_mid(i)=nlev
END DO

! fill whole array with 1. along diagnol so U(t+1) = u(t) and extract winds

DO k=1,nlev
  DO i=1,nmid
    jmid =index_mid(i)         ! location in full array
    terma(i,k) = 0.0
    termb(i,k) = 1.0
    termc(i,k) = 0.0
    u_com(i,k) = u(jmid,k)     ! winds at t
    v_com(i,k) = v(jmid,k)
    k_dif(i,k) = 0.0           ! set to zero
  END DO
END DO

DO k=2,nlev-1
  DO i=1,nmid
    jmid =index_mid(i)         ! location in full array
    jbase = 1

    DO j=1,nc_layers(i)
      IF ( k >= ic_base(i,j) .AND. k < ic_top(i,j)+2) THEN
        jbase = j
      END IF
    END DO              ! over mid-level layers

    IF (k == ic_base(i,jbase)) THEN

      recip_epsilon = 1.0/entrain_rate(i,k)
      length_scale = MIN(recip_epsilon,(z_theta(jmid,k)-z_cbase(i,jbase)) )

      dz_rho   = z_rho(jmid,k+1) - z_rho(jmid,k)
      dz_theta = z_theta(jmid,k) - z_theta(jmid,k-1)
      k_dif(i,k) = -1.0*massf(i,k)*cfact*length_scale/dz_rho
      grid_cor_fac = r_theta(jmid,k)*r_theta(jmid,k)/                  &
                                   (r_rho(jmid,k)*r_rho(jmid,k))

!           terma(i,k) = 0.0      ! already set
      termc(i,k) = timestep*k_dif(i,k)*grid_cor_fac                    &
                                            /(rho(jmid,k)*dz_theta)
      termb(i,k) = 1.0 -terma(i,k)-termc(i,k)


    ELSE IF (k > ic_base(i,jbase) .AND. k < ic_top(i,jbase)+1) THEN

      recip_epsilon = 1.0/entrain_rate(i,k)
      length_scale = MIN(recip_epsilon,(z_theta(jmid,k)-z_cbase(i,jbase)) )

      dz_rho   = z_rho(jmid,k+1) - z_rho(jmid,k)
      dz_theta = z_theta(jmid,k) - z_theta(jmid,k-1)

      k_dif(i,k) = -1.0*massf(i,k)*cfact*length_scale/dz_rho

      grid_cor_fac = r_theta(jmid,k)*r_theta(jmid,k)/                  &
                                   (r_rho(jmid,k)*r_rho(jmid,k))

      termc(i,k) = timestep*k_dif(i,k)*grid_cor_fac                    &
                                            /(rho(jmid,k)*dz_theta)

      grid_cor_fac =r_theta(jmid,k-1)*r_theta(jmid,k-1)/               &
                                    (r_rho(jmid,k)*r_rho(jmid,k))
      ! k_dif(i,k-1) already known
      terma(i,k) = timestep*k_dif(i,k-1)*grid_cor_fac                  &
                                            /(rho(jmid,k)*dz_theta)

      termb(i,k) = 1.0 -terma(i,k)-termc(i,k)

    ELSE IF (k == ic_top(i,jbase)+1) THEN

      ! For this level k_dif is zero as massf = 0.0

      dz_theta = z_theta(jmid,k) - z_theta(jmid,k-1)

      grid_cor_fac =r_theta(jmid,k-1)*r_theta(jmid,k-1)/               &
                                    (r_rho(jmid,k)*r_rho(jmid,k))
      terma(i,k) = timestep*k_dif(i,k-1)*grid_cor_fac                  &
                                            /(rho(jmid,k)*dz_theta)

!           termc(i,k) = 0.0         ! already set
      termb(i,k) = 1.0 -terma(i,k)-termc(i,k)

    END IF                  ! test on level number

  END DO
END DO

!-------------------------------------------------------------------------------
! Solve for winds at time T+1 by solving tridiagonal matrix
!-------------------------------------------------------------------------------

! DEPENDS ON: tridiag_all
CALL tridiag_all(nlev,nmid,nlev_mid,terma,termb,termc,u_com,ue_tp1)

! DEPENDS ON: tridiag_all
CALL tridiag_all(nlev,nmid,nlev_mid,terma,termb,termc,v_com,ve_tp1)

!-------------------------------------------------------------------------------
! Calculate wind increments and wind stress profiles for output.
!-------------------------------------------------------------------------------

DO k=2,nlev-1
  DO i=1,nmid
    jmid =index_mid(i)

    dubydt(jmid,k) = (ue_tp1(i,k)-u(jmid,k))/timestep
    dvbydt(jmid,k) = (ve_tp1(i,k)-v(jmid,k))/timestep

    uw(jmid,k) = k_dif(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))
    vw(jmid,k) = k_dif(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))

  END DO              ! mid point loop
END DO                ! level loop

!-------------------------------------------------------------------------------
IF (lhook) CALL dr_hook('MID_CONV_DIF_CMT',zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
END SUBROUTINE mid_conv_dif_cmt
