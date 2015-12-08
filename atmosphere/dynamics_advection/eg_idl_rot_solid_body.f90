! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_rot_solid_body_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_rot_solid_body(                                     &
                      l_shallow, l_rotate_grid,                       &
                      row_length, rows, n_rows, halo_i, halo_j,       &
                      offx, offy, model_levels,                       &
                      xi3_at_rho, xi3_at_u, xi3_at_v,                 &
                      intw_w2rho,  intw_rho2w,                        &
                      csxi1_p, csxi1_u, snxi1_p, snxi1_u,             &
                      csxi2_p, snxi2_p,                               &
                      t_surface, p_surface, two_omega,                &
                      grid_np_lon, grid_np_lat,                       &
                      aa_jet_u0, aa_jet_a, aa_jet_m, aa_jet_n,        &
                      f1_comp, f2_comp, f3_comp,                      &
                      u, v, w, u_adv, v_adv, w_adv,                   &
                      exner_star, theta, rho, exner)

USE eg_idl_1d_profs_mod

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi
USE atmos_constants_mod, ONLY : r,kappa,p_zero
USE earth_constants_mod, ONLY : earth_radius, g 
USE level_heights_mod,   ONLY : r_theta_levels
USE atm_fields_bounds_mod
USE timestep_mod, ONLY : timestep_number

IMPLICIT NONE
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for exact solid body rotation solution of deep atmosphere equations.
!   Isothermal case only. Version based on Williamson shallow water tests.
!  
!
! Method:
!   Simply set fields using analytic formulae. May need to balance
!   pressure field with respect to model discretisation?
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



LOGICAL, INTENT (IN) :: l_shallow, l_rotate_grid

! Array Dimensions
INTEGER, INTENT(IN) :: row_length   ! no. points on a row
INTEGER, INTENT(IN) :: rows         ! no. rows
INTEGER, INTENT(IN) :: n_rows       ! no. rows of v
INTEGER, INTENT(IN) :: halo_i       ! Size of halo x-direction
INTEGER, INTENT(IN) :: halo_j       ! Size of halo y-direction
INTEGER, INTENT(IN) :: offx         ! Small halo x-direction
INTEGER, INTENT(IN) :: offy         ! Small halo y-direction
INTEGER, INTENT(IN) :: model_levels ! number of model levels

! Physical Constants
REAL, INTENT(IN) :: t_surface,                                        &
                    p_surface
REAL, INTENT(IN) :: two_omega

! Idealised options

! Vertical co-ordinate information
REAL, INTENT(IN) ::                                                   &
  xi3_at_rho(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
               model_levels),                                         &
  xi3_at_u(-halo_i:row_length-1+halo_i, 1-halo_j:rows+halo_j,         &
         model_levels),                                               &
  xi3_at_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,        &
           model_levels)

REAL, INTENT(IN) ::                                                   &
  intw_w2rho(model_levels,2),  intw_rho2w(model_levels,2)

REAL, INTENT(IN) ::                                                   &
  csxi1_p(1-halo_i:row_length+halo_i),                                &
  csxi1_u(-halo_i:row_length-1+halo_i),                               &
  csxi2_p(1-halo_j:rows+halo_j),                                      &
  snxi1_p(1-halo_i:row_length+halo_i),                                &
  snxi1_u(-halo_i:row_length-1+halo_i),                               &
  snxi2_p(1-halo_j:rows+halo_j)

! Rotated grid
REAL, INTENT(IN)    :: grid_np_lon,                                   &
                       grid_np_lat

! Staniforth     & White
REAL, INTENT(IN)    :: aa_jet_u0,                                     &
                       aa_jet_a
INTEGER, INTENT(IN) :: aa_jet_m,                                      &
                       aa_jet_n

! Components of Coriolis
REAL, INTENT(INOUT) ::                                                &
  f1_comp(row_length,rows),                                           &
  f2_comp(row_length,rows),                                           &
  f3_comp(row_length,rows)

! Wind components
REAL, INTENT (OUT) ::                                                 &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),           &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),         &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL, INTENT (OUT) ::                                                 &
  u_adv(-halo_i:row_length-1+halo_i, 1-halo_j:rows+halo_j,            &
        model_levels),                                                &
  v_adv(1-halo_i:row_length+halo_i, -halo_j:n_rows+halo_j-1,          &
        model_levels),                                                &
  w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
        0:model_levels)

! Thermodynamic fields
REAL, INTENT (INOUT) ::                                               &
  theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),          &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

REAL, INTENT (INOUT) ::                                               &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy)


! Local Variables

INTEGER :: i, j, k    ! Loop indices

! Required for rotated grid
REAL :: lon_p, lat_p

REAL :: t0, s, u0_const, f_sb, dz
REAL :: sin_lat, cos_lat, sin_lon, cos_lon
REAL :: cos_lat_p, sin_lat_p, cos_lon_p, sin_lon_p
REAL :: exner_surface, theta_bar

REAL :: exner_1d(0:model_levels+1), rho_1d(model_levels)
REAL :: thetav_1d(0:model_levels), grav_1d(0:model_levels)
REAL :: dtheta_dz1(3)
REAL :: u0, a
INTEGER :: m, n
REAL :: t_surf,  p_surf


!- End of header

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_ROT_SOLID_BODY',zhook_in,zhook_handle)


IF (l_rotate_grid) THEN

! Rotated pole case
  lon_p=grid_np_lon*pi/180.0
  lat_p=grid_np_lat*pi/180.0

  sin_lon_p=SIN(lon_p)
  cos_lon_p=COS(lon_p)

  sin_lat_p=SIN(lat_p)
  cos_lat_p=COS(lat_p)

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      sin_lon=snxi1_p(i)*cos_lon_p-csxi1_p(i)*sin_lon_p
      cos_lon=csxi1_p(i)*cos_lon_p+snxi1_p(i)*sin_lon_p

      sin_lat=snxi2_p(j)*sin_lat_p-cos_lon*csxi2_p(j)*cos_lat_p

      cos_lat=SQRT(1.0-sin_lat*sin_lat)

! Components of Coriolis
      IF (l_shallow) THEN
        f1_comp(i,j)=0.0
        f2_comp(i,j)=0.0
      ELSE
        f1_comp(i,j)=two_omega*sin_lon*cos_lat_p
        f2_comp(i,j)=two_omega*(csxi2_p(j)*sin_lat_p +                &
                                cos_lon*snxi2_p(j)*cos_lat_p)
      END IF
      f3_comp(i,j)=two_omega*sin_lat

    END DO
  END DO

END IF

IF (timestep_number /= 1) THEN

  IF (lhook) CALL dr_hook('EG_IDL_ROT_SOLID_BODY',zhook_out,zhook_handle)
  RETURN

END IF






dtheta_dz1 = 0.0
Do k = 0, model_levels
   grav_1d(k)    = g*( earth_radius/r_theta_levels(1,1,k) )**2
Enddo

t_surf = t_surface
p_surf = p_surface

      CALL eg_idl_1d_profs(                                           &
                       exner_1d,                                      &
                       rho_1d,                                        &
                       thetav_1d,                                     &
                       2,                                             &
                       t_surf, p_surf,                                &
                       1.0, 1.0,                                      &
                       intw_rho2w, intw_w2rho,                        &
                       r_theta_levels(1,1,:), xi3_at_rho(1,1,:),      &
                       1.0, earth_radius, dtheta_dz1, grav_1d,        &
                       model_levels)

t0=t_surface

u0=aa_jet_u0
a=aa_jet_a
m=aa_jet_m
n=aa_jet_n

u0_const=aa_jet_u0*(aa_jet_u0+two_omega*earth_radius)/(r*t0)

IF (l_rotate_grid) THEN

! Rotated pole case
  lon_p=grid_np_lon*pi/180.0
  lat_p=grid_np_lat*pi/180.0

  sin_lon_p=SIN(lon_p)
  cos_lon_p=COS(lon_p)

  sin_lat_p=SIN(lat_p)
  cos_lat_p=COS(lat_p)

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      sin_lon=snxi1_p(i)*cos_lon_p-csxi1_p(i)*sin_lon_p
      cos_lon=csxi1_p(i)*cos_lon_p+snxi1_p(i)*sin_lon_p

      sin_lat=snxi2_p(j)*sin_lat_p-cos_lon*csxi2_p(j)*cos_lat_p

      cos_lat=SQRT(1.0-sin_lat*sin_lat)

! Theta at surface
      s=cos_lat
      f_sb=0.5*u0_const*s**2
! Modify p_surface to take account of orography
      exner_star(i,j)=p_surface*EXP(f_sb)*                             &
                      EXP(-(r_theta_levels(i,j,0)-earth_radius)*       &
                      g/(r*t0 ))
      exner_surface=(exner_star(i,j)/p_zero)**kappa
      theta(i,j,0)=t0/exner_surface

! Modify hydrostatic pressure, theta and rho
      DO k=1,model_levels
        IF (l_shallow) THEN
          s = cos_lat
        ELSE
          s=xi3_at_rho(i,j,k)*cos_lat/earth_radius
        END IF
        f_sb=0.5*u0_const*s**2
        dz=xi3_at_rho(i,j,k)-earth_radius
        exner(i,j,k)=exner_1d(k)*EXP(kappa*f_sb)
        theta(i,j,k)=2.0*t0/exner(i,j,k)-theta(i,j,k-1)
        theta_bar  = (intw_w2rho(k,2)*theta(i,j,k-1)                  &
                     +intw_w2rho(k,1)*theta(i,j,k))
        rho(i,j,k) = p_zero*exner(i,j,k)**((1.0-kappa)/kappa) /       &
                     (r*theta_bar)
      END DO
    END DO
  END DO

  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

      sin_lon=snxi1_u(i)*cos_lon_p-csxi1_u(i)*sin_lon_p
      cos_lon=csxi1_u(i)*cos_lon_p+snxi1_u(i)*sin_lon_p

      DO k=1,model_levels
        u(i,j,k)=aa_jet_u0*xi3_at_u(i,j,k)*(csxi2_p(j)*sin_lat_p +    &
                 cos_lon*snxi2_p(j)*cos_lat_p)/earth_radius
        u_adv(i,j,k)=u(i,j,k)
      END DO

    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end

      sin_lon=snxi1_p(i)*cos_lon_p-csxi1_p(i)*sin_lon_p
      cos_lon=csxi1_p(i)*cos_lon_p+snxi1_p(i)*sin_lon_p

      DO k=1,model_levels
        v(i,j,k)=-aa_jet_u0*xi3_at_v(i,j,k)*sin_lon*cos_lat_p /       &
                  earth_radius
        v_adv(i,j,k)=v(i,j,k)
      END DO

    END DO
  END DO

  w(:,:,:)=0.0
  w_adv(:,:,:)=0.0

ELSE

! Modify hydrostatic pressure, theta and rho
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      s=csxi2_p(j)
      f_sb=(u0*s**m)**2*(1.0/(2*m)-2*(a*s)**n/(2*m+n) +               &
                         0.5*(a*s)**(2*n)/(m+n)) +                    &
           two_omega*earth_radius*u0*s**(m+1)*(1.0/(m+1) -            &
           (a*s)**n/(m+n+1))
      f_sb=f_sb/(r*t0)
! Modify p_surface to take account of orography
      exner_star(i,j)=p_surface*EXP(f_sb)*                             &
                      EXP(-(r_theta_levels(i,j,0)-earth_radius)*       &
                      g/(r*t0 ))
      exner_surface=(exner_star(i,j)/p_zero)**kappa
      theta(i,j,0)=t0/exner_surface
    END DO
  END DO

  DO k=1,model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (l_shallow) THEN
          s = csxi2_p(j)
        ELSE
          s=xi3_at_rho(i,j,k)*csxi2_p(j)/earth_radius
        END IF
        f_sb=(u0*s**m)**2*(1.0/(2*m)-2*(a*s)**n/(2*m+n) +             &
                           0.5*(a*s)**(2*n)/(m+n)) +                  &
             two_omega*earth_radius*u0*s**(m+1)*(1.0/(m+1) -          &
             (a*s)**n/(m+n+1))
        f_sb=f_sb/(r*t0)
        dz=xi3_at_rho(i,j,k)-earth_radius
        exner(i,j,k)=exner_1d(k)*EXP(kappa*f_sb)
        theta(i,j,k)=2.0*t0/exner(i,j,k)-theta(i,j,k-1)
        theta_bar=0.5*(theta(i,j,k-1)+theta(i,j,k))
        rho(i,j,k) = p_zero*exner(i,j,k)**((1.0-kappa)/kappa) /       &
                     (r*theta_bar)
      END DO
    END DO
  END DO

! Set wind components
  DO k=1,model_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        IF (l_shallow) THEN
          s=csxi2_p(j)
        ELSE
          s=xi3_at_u(i,j,k)*csxi2_p(j)/earth_radius
        END IF
        u(i,j,k) = aa_jet_u0*s**m*(1.0-(a*s)**n)
        u_adv(i,j,k)=u(i,j,k)
      END DO
    END DO
  END DO

  v(:,:,:)=0.0
  v_adv(:,:,:)=0.0
  w(:,:,:)=0.0
  w_adv(:,:,:)=0.0

END IF
! w(row_length/2,rows/2,1) = 0.001

IF (lhook) CALL dr_hook('EG_IDL_ROT_SOLID_BODY',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_rot_solid_body
END MODULE eg_idl_rot_solid_body_mod
