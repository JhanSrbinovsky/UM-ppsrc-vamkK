! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_deep_baroclinic_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_deep_baroclinic(                                      &
                      row_length, rows, n_rows, halo_i, halo_j,         &
                      offx, offy, model_levels,                         &
                      xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,     &
                      g_theta, u, v, w, u_adv, v_adv, w_adv,            &
                      theta, rho, exner, exner_star, L_baro_perturbed,  &
                      T0_E, T0_P, b_const, k_const, l_shallow,          &
                      l_rotate_grid, grid_np_lon, grid_np_lat,          &
                      f1_comp, f2_comp, f3_comp)
 
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R, kappa, p_zero
USE earth_constants_mod, ONLY: Earth_radius, g, omega
USE eg_idl_deep_baro_mod
USE eg_idl_1d_profs_mod
USE horiz_grid_mod
USE conversions_mod

USE proc_info_mod,       ONLY : n_proc,                                 &
                                global_row_length,                      &
                                at_extremity,                           &
                                gc_proc_row_group,                      &
                                model_domain,                           &
                                me
USE field_types,         ONLY : fld_type_u, fld_type_v
USE eg_v_at_poles_mod
USE UM_ParParams
USE atm_fields_bounds_mod
USE eg_parameters_mod, ONLY : pole_consts
USE eg_swap_bounds_mod

USE timestep_mod, ONLY : timestep_number

IMPLICIT NONE
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for baroclinic wave solution of deep atmosphere equations.
!   Isothermal case only. Version based on Staniforth deep tests.
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
 
 
! Array Dimensions-f_
INTEGER, INTENT(IN) :: row_length   ! no. points on a row
INTEGER, INTENT(IN) :: rows         ! no. rows
INTEGER, INTENT(IN) :: n_rows       ! no. rows of v
INTEGER, INTENT(IN) :: halo_i       ! Size of halo x-direction
INTEGER, INTENT(IN) :: halo_j       ! Size of halo y-direction
INTEGER, INTENT(IN) :: offx         ! Small halo x-direction
INTEGER, INTENT(IN) :: offy         ! Small halo y-direction
INTEGER, INTENT(IN) :: model_levels ! number of model levels
  
! Idealised options
LOGICAL, INTENT(IN) :: L_baro_perturbed, l_shallow

! Deep baroclinic options
REAL,    INTENT(IN) :: T0_E, T0_P 
INTEGER, INTENT(IN) :: b_const, k_const
 
! Vertical co-ordinate information
REAL, INTENT(IN) ::                                                     &
  xi3_at_theta(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
               0:model_levels),                                         &
  xi3_at_rho(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
               model_levels),                                           &
  xi3_at_u(-halo_i:row_length-1+halo_i, 1-halo_j:rows+halo_j,           &
         model_levels),                                                 &
  xi3_at_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,          &
         model_levels)

! Gravitational acceleration
REAL, INTENT(IN) ::                                                     &
  g_theta(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)
 
! Wind components
REAL, INTENT (OUT) ::                                                   &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),             &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),           &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)
 
 
REAL, INTENT (OUT) ::                                                   &
  u_adv(-halo_i:row_length-1+halo_i, 1-halo_j:rows+halo_j,              &
        model_levels),                                                  &
  v_adv(1-halo_i:row_length+halo_i, -halo_j:n_rows+halo_j-1,            &
        model_levels),                                                  &
  w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,               &
        0:model_levels)
 
! Thermodynamic fields
REAL, INTENT (INOUT) ::                                                 &
  theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)
 
REAL, INTENT (INOUT) ::                                                 &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy)
 
! Local parameters needed for call to eg_idl_1d_profs:
!
! Profile number used in eg_idl_1d_profs
INTEGER, PARAMETER :: Tprofile_number = 4
! Cause eg_idl_1d_profs to get height information from xi3 arrays
REAL,    PARAMETER :: z_top_of_model  = 1.0
! Array needed by call to eg_idl_1d_profs but not used
REAL,    PARAMETER :: dtheta_dz1(3) = (/0.0, 0.0, 0.0/)

! Local Variables

! Loop indices
INTEGER :: i, j, k 
           
! First guess for Newton solver
REAL    :: T_1st_guess
REAL    :: p_1st_guess
 
REAL    :: exner_pro(0:model_levels+1)

REAL :: tau2_term, T_temp, r_height, temp, tau1_u, tau2_u, tau2_int_u

REAL :: taper

REAL :: xi2_rot, sn_lonmlonp, cs_lonmlonp, cs_Rot, sn_Rot,            &
        lat_p, lon_p, cos_lat_p, sin_lat_p, two_omega,                &
        u_at_v, v_at_u, XX, YY, xi1_rot

REAL ::                                                               &
  u_rot(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  v_rot(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels)

! Rotated grid
REAL, INTENT(IN)    :: grid_np_lon,                                   &
                       grid_np_lat
LOGICAL, INTENT(IN) :: l_rotate_grid
! Components of Coriolis
REAL, INTENT(INOUT) ::                                                &
  f1_comp(row_length,rows),                                           &
  f2_comp(row_length,rows),                                           &
  f3_comp(row_length,rows)

REAL :: div_max,div_min,vort_max,vort_min

INTEGER :: i_err

INTEGER :: pert_type = 3

! External functions
REAL, EXTERNAL   :: eg_baro_pert
REAL, EXTERNAL   :: eg_deep_baro_pert


! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_IDL_DEEP_BAROCLINIC',zhook_in,zhook_handle)

! Setup constants

 T0 = 0.5*(T0_E + T0_P)
 B = (T0 - T0_P)/(T0*T0_P)
 C = (k_const+2.0)/2.0 * (T0_E - T0_P)/(T0_E*T0_P)
 H = R*T0/g
 k_const_save = k_const
 b_const_save = b_const
 l_shallow_save = l_shallow

! Section to rotate grid
IF (l_rotate_grid) THEN
! Rotated pole case
  IF (me == 0) WRITE(6,fmt='(A,2E16.8)')                              &
               'Rotated grid pole at (lat,lon): ',                    &
                        grid_np_lat,grid_np_lon

  lon_p=grid_np_lon*pi/180.0
  lat_p=grid_np_lat*pi/180.0

  sin_lat_p=SIN(lat_p)
  cos_lat_p=COS(lat_p)

  two_omega = 2.0*omega

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! Components of Coriolis
      IF (l_shallow) THEN
        f1_comp(i,j)=0.0
        f2_comp(i,j)=0.0
      ELSE
        f1_comp(i,j) = -two_omega*cos_lat_p*SIN(xi1_p(i))
        f2_comp(i,j) =  two_omega*(sin_lat_p*COS(xi2_p(j)) -            &
                                   cos_lat_p*SIN(xi2_p(j))*COS(xi1_p(i)))
      END IF
      f3_comp(i,j)=two_omega*(sin_lat_p*SIN(xi2_p(j)) +                 &
                              cos_lat_p*COS(xi2_p(j))*COS(xi1_p(i)))
    END DO
  END DO
END IF
 


IF (timestep_number /= 1) THEN

  IF (lhook) CALL dr_hook('EG_IDL_DEEP_BAROCLINIC',zhook_out,zhook_handle)
  RETURN

END IF



! Put each column into hydrostatic balance with temperature of
! analytic steady-state
T_1st_guess=T0
p_1st_guess=1.0e5
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF (l_rotate_grid) THEN
      xi2_rot = ASIN( SIN(xi2_p(j))*sin_lat_p                           &
                     +COS(xi2_p(j))*COS(xi1_p(i))*cos_lat_p)
    ELSE
      xi2_rot = xi2_p(j)
    END IF

! Setup tau fields
    DO k = 0, model_levels
      IF ( k == 0 ) THEN
        r_height = xi3_at_theta(i,j,0)+Earth_radius
        r_height = Earth_radius
      ELSE
!         r_height = xi3_at_rho(i,j,k)
        r_height = xi3_at_theta(i,j,k)
      END IF
      temp = ((r_height-Earth_radius)/(b_const*H))**2

      tau1(k) = 1.0/T0*EXP(eg_lapse/T0*(r_height-Earth_radius))         &
              + B*(1.0 - 2.0*temp)*EXP(-(temp))

      tau2(k) = C*(1.0-2.0*temp)*EXP(-temp) 

      tau1_int(k) = 1.0/eg_lapse                                        &
                   *(EXP(eg_lapse/T0*(r_height-Earth_radius))-1.0)      &
                  + B*(r_height-Earth_radius)*EXP(-temp)

      tau2_int(k) = C*(r_height-Earth_radius)*EXP(-temp)
    END DO
    CALL eg_idl_1d_profs(exner_pro, rho(i,j,:), theta(i,j,:),           &
                         Tprofile_number, T_1st_guess, p_1st_guess,     &
                         xi1_p(i), xi2_rot,                             &
                         intw_rho2w, intw_w2rho,                        &
                         xi3_at_theta(i,j,:), xi3_at_rho(i,j,:),        &
                         z_top_of_model, Earth_radius, dtheta_dz1,      &
                         g_theta(i,j,:), model_levels)
 
    exner(i,j,:)=exner_pro(1:model_levels+1)
 
    exner_star(i,j) = exner_pro(0)
  END DO
END DO

! ! Balanced pressure field
IF (L_baro_Perturbed) THEN
  DO k = 1, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (l_rotate_grid) THEN
          xi2_rot = ASIN( SIN(xi2_p(j))*sin_lat_p                       &
                         +COS(xi2_p(j))*COS(xi1_p(i)*cos_lat_p))
          XX = - sin_lat_p*COS(xi2_p(j))*SIN(xi1_p(i))
          YY = SIN(xi2_rot)*cos_lat_p - COS(xi2_p(j))*COS(xi1_p(i))
          IF( ABS(YY) < 1.0E-5 ) THEN
            IF ( ABS (XX) < 1.0E-5 ) THEN
              xi1_rot = lon_p  + pi 
            ELSE
              xi1_rot = lon_p  + 2.0*pi 
            END IF
          ELSE
            xi1_rot = lon_p + ATAN2(XX,YY) + pi 
          END IF            
        ELSE
          xi2_rot = xi2_p(j)
          xi1_rot = xi1_u(i) 
        END IF

        r_height = xi3_at_rho(i,j,k)

        temp = ((r_height-Earth_radius)/(b_const*H))**2
  
        tau1_u = 1.0/T0*EXP(eg_lapse/T0*(r_height-Earth_radius))        &
                + B*(1.0 - 2.0*temp)*EXP(-(temp))
        tau2_u = C*(1.0-2.0*temp)*EXP(-temp) 
        tau2_int_u = C*(r_height-Earth_radius)*EXP(-temp)

        IF ( L_shallow ) r_height =  Earth_radius
        ! Unperturbed wind field
        tau2_term = (r_height/Earth_radius*COS(xi2_rot)) ** k_const     &
                  - (k_const/(k_const + 2.0))*                          &
                    (r_height/Earth_radius*COS(xi2_rot))                &
                                                   **(k_const + 2.0)

        T_temp = (Earth_radius/r_height)**2 /(tau1_u - tau2_u*tau2_term) 

        IF ( xi3_at_rho(i,j,k) - Earth_radius <= taper_top ) THEN
          taper = 1.0                                                   &
                - 3.0*((xi3_at_rho(i,j,k)-Earth_radius)/taper_top)**2   &
                + 2.0*((xi3_at_rho(i,j,k)-Earth_radius)/taper_top)**3 
                    
        ! DEPENDS ON: eg_deep_baro_pert
          exner(i,j,k) = exner(i,j,k) * EXP(-taper                      &
                 *eg_deep_baro_pert(xi1_rot,xi2_rot,pert_type,3)*       &
                  2.0*omega*SIN(xi2_rot)/(R*T_temp))
        ! Recompute density
          rho(i,j,k) = p_zero/(R*                                       &
                      (intw_w2rho(k,1)*theta(i,j,k)                     &
                     + intw_w2rho(k,2)*theta(i,j,k-1) ))                &
                     *(exner(i,j,k))**((1.0-kappa)/kappa)  
        END IF
      END DO
    END DO
  END DO
END IF
 
! Set wind components
DO k=1,model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

    IF (l_rotate_grid) THEN
      xi2_rot = ASIN( SIN(xi2_p(j))*sin_lat_p                           &
                     +COS(xi2_p(j))*COS(xi1_u(i))*cos_lat_p)
      XX = - sin_lat_p*COS(xi2_p(j))*SIN(xi1_u(i))
      YY = SIN(xi2_rot)*cos_lat_p - COS(xi2_p(j))*COS(xi1_u(i))
      IF( ABS(YY) < 1.0E-5 ) THEN
        IF ( ABS (XX) < 1.0E-5 ) THEN
          xi1_rot = lon_p  + pi 
        ELSE
          xi1_rot = lon_p  + 2.0*pi 
        END IF
      ELSE
        xi1_rot = lon_p + ATAN2(XX,YY) + pi 
      END IF
    ELSE
      xi2_rot = xi2_p(j)
      xi1_rot = xi1_u(i)   
    END IF

      r_height = xi3_at_u(i,j,k)

      temp = ((r_height-Earth_radius)/(b_const*H))**2
  
      tau1_u = 1.0/T0*EXP(eg_lapse/T0*(r_height-Earth_radius))          &
              + B*(1.0 - 2.0*temp)*EXP(-(temp))
      tau2_u = C*(1.0-2.0*temp)*EXP(-temp) 
      tau2_int_u = C*(r_height-Earth_radius)*EXP(-temp)


      IF ( L_shallow ) r_height =  Earth_radius
      ! Unperturbed wind field
      tau2_term = (r_height/Earth_radius*COS(xi2_rot)) ** k_const       &
                - (k_const/(k_const + 2.0))*                            &
                  (r_height/Earth_radius*COS(xi2_rot)) **(k_const + 2.0)

      T_temp = (Earth_radius/r_height) ** 2 /(tau1_u - tau2_u*tau2_term)

! This is the wind proxy: U
      u(i,j,k) = (g/Earth_radius)*k_const*tau2_int_u                    &
                *((r_height/Earth_radius*COS(xi2_rot))**(k_const-1)     &
                - (r_height/Earth_radius*COS(xi2_rot))**(k_const+1))    &
                 *T_temp

      u(i,j,k) = - omega*r_height*COS(xi2_rot)                          &
                 + SQRT( (omega*r_height*COS(xi2_rot))**2               &
                        + r_height*COS(xi2_rot)*u(i,j,k) )

! Add perturbation
      IF (L_baro_Perturbed) THEN
        IF ( xi3_at_u(i,j,k) - Earth_radius <= taper_top ) THEN
          taper = 1.0-3.0*((xi3_at_u(i,j,k)-Earth_radius)/taper_top)**2 &
                    + 2.0*((xi3_at_u(i,j,k)-Earth_radius)/taper_top)**3  
! DEPENDS ON: eg_deep_baro_pert
          u(i,j,k) = u(i,j,k) + taper                                   &
                    *eg_deep_baro_pert(xi1_rot,xi2_rot,pert_type,1)

! Alternative option
! ! DEPENDS ON: eg_deep_baro_pert
!           psi_p = taper*eg_deep_baro_pert(xi1_u(i),xi2_v(j),pert_type,3)
! ! DEPENDS ON: eg_deep_baro_pert
!           psi_m = taper*eg_deep_baro_pert(xi1_u(i),xi2_v(j-1),pert_type,3)
!          
!           u(i,j,k) = u(i,j,k) - 1.0/Earth_radius                      &
!                            *(psi_p-psi_m)/(xi2_v(j)-xi2_v(j-1))
        END IF
      END IF

      u_adv(i,j,k)=u(i,j,k)
    END DO
  END DO
END DO

v=0.0
IF ( L_baro_Perturbed .AND. pert_type /= 1 ) THEN
DO k=1,model_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
    IF (l_rotate_grid) THEN
      xi2_rot = ASIN( SIN(xi2_v(j))*sin_lat_p                           &
                     +COS(xi2_v(j))*COS(xi1_p(i))*cos_lat_p)
      XX = - sin_lat_p*COS(xi2_v(j))*SIN(xi1_p(i))
      YY = SIN(xi2_rot)*cos_lat_p - COS(xi2_v(j))*COS(xi1_p(i))
      IF( ABS(YY) < 1.0E-5 ) THEN
        IF ( ABS (XX) < 1.0E-5 ) THEN
          xi1_rot = lon_p  + pi 
        ELSE
          xi1_rot = lon_p  + 2.0*pi 
        END IF
      ELSE
        xi1_rot = lon_p + ATAN2(XX,YY) + pi 
      END IF
    ELSE
      xi2_rot = xi2_v(j)
      xi1_rot = xi1_p(i)
    END IF

! Add perturbation
      IF ( xi3_at_v(i,j,k) - Earth_radius <= taper_top ) THEN
        taper = 1.0 - 3.0*((xi3_at_v(i,j,k)-Earth_radius)/taper_top)**2 &
                    + 2.0*((xi3_at_v(i,j,k)-Earth_radius)/taper_top)**3  
! DEPENDS ON: eg_deep_baro_pert
        v(i,j,k) =  taper*eg_deep_baro_pert(xi1_rot,xi2_rot,pert_type,2)

! Alternative option
! ! DEPENDS ON: eg_deep_baro_pert
!           psi_p = taper*eg_deep_baro_pert(xi1_u(i),xi2_v(j),pert_type,3)
! ! DEPENDS ON: eg_deep_baro_pert
!           psi_m = taper*eg_deep_baro_pert(xi1_u(i-1),xi2_v(j),pert_type,3)
!          
!           v(i,j,k) = v(i,j,k) + 1.0/(Earth_radius*COS(xi2_v(j)))         &
!                      *(psi_p-psi_m)/(xi1_u(i)-xi1_u(i-1))
        END IF
    END DO
  END DO
END DO
END IF

! Rotate wind vectors
IF (l_rotate_grid) THEN

  CALL eg_swap_bounds(u,udims_s,fld_type_u,.TRUE.)
  CALL eg_swap_bounds(v,vdims_s,fld_type_v,.TRUE.)

  DO k=1,model_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        xi2_rot = ASIN( SIN(xi2_p(j))*sin_lat_p                         &
                       +COS(xi2_p(j))*COS(xi1_u(i))*cos_lat_p)

        sn_lonmlonp = - COS(xi2_p(j))*SIN(xi1_u(i))/COS(xi2_rot)
        cs_lonmlonp = (SIN(xi2_rot)*cos_lat_p                           &
                     - COS(xi2_p(j))*COS(xi1_u(i)))                     &
                     /(COS(xi2_rot)*sin_lat_p)

        cs_Rot = - COS(xi1_u(i))*cs_lonmlonp                            &
                 - SIN(xi1_u(i))*sn_lonmlonp*sin_lat_p
        sn_Rot = COS(xi1_u(i))*SIN(xi2_rot)*sn_lonmlonp - SIN(xi1_u(i)) &
               *(COS(xi2_rot)*cos_lat_p                                 &
               + SIN(xi2_rot)*sin_lat_p*cs_lonmlonp)

! Add perturbation
        v_at_u = 0.0
        IF (L_baro_Perturbed) THEN
          IF ( xi3_at_u(i,j,k) - Earth_radius <= taper_top ) THEN
            XX = - sin_lat_p*COS(xi2_p(j))*SIN(xi1_u(i))
            YY = SIN(xi2_rot)*cos_lat_p - COS(xi2_p(j))*COS(xi1_u(i))
            IF( ABS(YY) < 1.0E-5 ) THEN
              IF ( ABS (XX) < 1.0E-5 ) THEN
                xi1_rot = lon_p  + pi 
              ELSE
                xi1_rot = lon_p  + 2.0*pi 
              END IF
            ELSE
              xi1_rot = lon_p + ATAN2(XX,YY) + pi 
            END IF

            taper = 1.0 - 3.0*((xi3_at_u(i,j,k)-Earth_radius)           &
                    /taper_top)**2                                      &
                   + 2.0*((xi3_at_u(i,j,k)-Earth_radius)/taper_top)**3  
! DEPENDS ON: eg_deep_baro_pert
            v_at_u = taper*eg_deep_baro_pert(xi1_rot,xi2_rot,pert_type,2)
          END IF
        END IF

        u_rot(i,j,k) = u(i,j,k)*cs_Rot + v_at_u*sn_Rot
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        xi2_rot = ASIN( SIN(xi2_v(j))*sin_lat_p                         &
                       +COS(xi2_v(j))*COS(xi1_p(i))*cos_lat_p)
        sn_lonmlonp = - COS(xi2_v(j))*SIN(xi1_p(i))/COS(xi2_rot)
        cs_lonmlonp = (SIN(xi2_rot)*cos_lat_p                           &
                     - COS(xi2_v(j))*COS(xi1_p(i)))                     &
                     /(COS(xi2_rot)*sin_lat_p)

        cs_Rot = - COS(xi1_p(i))*cs_lonmlonp                            &
                 - SIN(xi1_p(i))*sn_lonmlonp*sin_lat_p
        sn_Rot = COS(xi1_p(i))*SIN(xi2_rot)*sn_lonmlonp - SIN(xi1_p(i)) &
               *(COS(xi2_rot)*cos_lat_p                                 &
              + SIN(xi2_rot)*sin_lat_p*cs_lonmlonp)

! Compute u at point xi1_rot, xi2_rot
        XX = - sin_lat_p*COS(xi2_v(j))*SIN(xi1_p(i))
        YY = SIN(xi2_rot)*cos_lat_p - COS(xi2_v(j))*COS(xi1_p(i))
      IF( ABS(YY) < 1.0E-5 ) THEN
        IF ( ABS (XX) < 1.0E-5 ) THEN
          xi1_rot = lon_p  + pi 
        ELSE
          xi1_rot = lon_p  + 2.0*pi 
        END IF
      ELSE
        xi1_rot = lon_p + ATAN2(XX,YY) + pi 
      END IF

        r_height = xi3_at_v(i,j,k)

        temp = ((r_height-Earth_radius)/(b_const*H))**2
  
        tau1_u = 1.0/T0*EXP(eg_lapse/T0*(r_height-Earth_radius))        &
                + B*(1.0 - 2.0*temp)*EXP(-(temp))
        tau2_u = C*(1.0-2.0*temp)*EXP(-temp) 
        tau2_int_u = C*(r_height-Earth_radius)*EXP(-temp)

        IF ( L_shallow ) r_height =  Earth_radius
      ! Unperturbed wind field
        tau2_term = (r_height/Earth_radius*COS(xi2_rot))**k_const       &
                  - (k_const/(k_const + 2.0))*                          &
                    (r_height/Earth_radius*COS(xi2_rot))**(k_const + 2.0)

        T_temp = (Earth_radius/r_height)**2 /(tau1_u -tau2_u*tau2_term)

! This is the wind proxy: U
        u_at_v = (g/Earth_radius)*k_const*tau2_int_u                    &
                  *((r_height/Earth_radius*COS(xi2_rot))**(k_const-1)   &
                  - (r_height/Earth_radius*COS(xi2_rot))**(k_const+1))  &
                   *T_temp

        u_at_v = - omega*r_height*COS(xi2_rot)                          &
                   + SQRT( (omega*r_height*COS(xi2_rot))**2             &
                          + r_height*COS(xi2_rot)*u_at_v )

! Add perturbation
        IF (L_baro_Perturbed) THEN
          IF ( xi3_at_v(i,j,k) - Earth_radius <= taper_top ) THEN
            taper = 1.0 - 3.0*((xi3_at_v(i,j,k)-Earth_radius)           &
                   /taper_top)**2                                       &
                  + 2.0*((xi3_at_v(i,j,k)-Earth_radius)/taper_top)**3  
! DEPENDS ON: eg_deep_baro_pert
            u_at_v = u_at_v + taper                                     &
                     *eg_deep_baro_pert(xi1_rot,xi2_rot,pert_type,1)
          END IF
        END IF

        v_rot(i,j,k) = - u_at_v*sn_Rot + v(i,j,k)*cs_Rot
      END DO
    END DO
  END DO

  IF( model_domain == mt_global ) THEN
    IF( at_extremity(PSouth) ) THEN

      CALL eg_v_at_poles(u_rot,v_rot, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)
    END IF

    IF( at_extremity(PNorth) ) THEN

      CALL eg_v_at_poles(u_rot,v_rot, -1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)

    END IF

  END IF

  DO k=1,model_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
          u(i,j,k) = u_rot(i,j,k)
          u_adv(i,j,k)=u(i,j,k)
      END DO
    END DO
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        v(i,j,k) = v_rot(i,j,k)
        v_adv(i,j,k)=v(i,j,k)
      END DO
    END DO
  END DO

END IF

! Compute Max/Min initial divergence/vorticity on level 1
k=1
div_max  = -1.0
div_min  = 1.0
vort_max = -1.0
vort_min = 1.0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end  
      temp = ((u(i,j,k) - u(i-1,j,k))/(xi1_u(i)-xi1_u(i-1)) &
            +(COS(xi2_v(j))*v(i,j,k)                        &
            - COS(xi2_v(j-1))*v(i,j-1,k))                   &
                                    /(xi2_v(j)-xi2_v(j-1))) &
            *1.0/(Earth_radius*COS(xi2_p(j)))
      IF ( temp > div_max ) div_max = temp
      IF ( temp < div_min ) div_min = temp
    END DO
  END DO
  DO j = vdims%j_start+1, vdims%j_end-1
    DO i = udims%i_start+1, udims%i_end-1
      temp = ((v(i+1,j,k) - v(i,j,k))/(xi1_p(i+1)-xi1_p(i)) &
            +(COS(xi2_p(j+1))*u(i,j+1,k)                    &
            - COS(xi2_p(j))*u(i,j,k))                       &
                                   /(xi2_p(j)-xi2_p(j-1)))  &
             *1.0/(Earth_radius*COS(xi2_v(j)))
      IF ( temp > vort_max ) vort_max = temp
      IF ( temp < vort_min ) vort_min = temp
    END DO
  END DO

  CALL gc_rmax(1,n_proc,i_err,div_max)
  CALL gc_rmin(1,n_proc,i_err,div_min)
  CALL gc_rmax(1,n_proc,i_err,vort_max)
  CALL gc_rmin(1,n_proc,i_err,vort_min)

  IF (me == 0) THEN
    WRITE(6,fmt='(A,2E16.8)') 'Max/Min initial divergence = ',            &
                          div_max,div_min
    WRITE(6,fmt='(A,2E16.8)') 'Max/Min initial vorticity  = ',            &
                          vort_max,vort_min
  END IF

  v_adv = v
  w     = 0.0
  w_adv = 0.0
 
  IF (lhook) CALL dr_hook('EG_IDL_DEEP_BAROCLINIC',zhook_out,zhook_handle)
 
  END SUBROUTINE eg_idl_deep_baroclinic
  END MODULE eg_idl_deep_baroclinic_mod
