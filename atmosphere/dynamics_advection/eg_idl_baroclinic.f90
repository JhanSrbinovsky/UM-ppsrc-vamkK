! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_baroclinic_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_baroclinic(                                           &
                      row_length, rows, n_rows, halo_i, halo_j,         &
                      offx, offy, model_levels,                         &
                      xi3_at_theta, xi3_at_rho, xi3_at_u,g_theta,       &
                      u, v, w, u_adv, v_adv, w_adv,                     &
                      theta, rho, exner, exner_star, L_baro_perturbed)
 
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R, kappa, p_zero
USE earth_constants_mod, ONLY: Earth_radius, g
USE eg_idl_baro_mod,     ONLY: T0
USE eg_idl_1d_profs_mod
USE horiz_grid_mod
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
LOGICAL, INTENT(IN) :: L_baro_perturbed
 
! Vertical co-ordinate information
REAL, INTENT(IN) ::                                                     &
  xi3_at_theta(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
               0:model_levels),                                         &
  xi3_at_rho(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
               model_levels),                                           &
  xi3_at_u(-halo_i:row_length-1+halo_i, 1-halo_j:rows+halo_j,           &
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
 
REAL        :: eta_p_u(model_levels)

! Local parameters needed for call to eg_idl_1d_profs:
!
! Profile number used in eg_idl_1d_profs
INTEGER, PARAMETER :: Tprofile_number = 3
! Cause eg_idl_1d_profs to get height information from xi3 arrays
REAL,    PARAMETER :: z_top_of_model  = 1.0
! Array needed by call to eg_idl_1d_profs but not used
REAL,    PARAMETER :: dtheta_dz1(3) = (/0.0, 0.0, 0.0/)

! Local Variables

! Loop indices
INTEGER :: i,                                                           &
           j,                                                           &
           k

! First guess for Newton solver
REAL    :: T_1st_guess
REAL    :: p_1st_guess

REAL    :: eta_surface
 
REAL    :: exner_pro(0:model_levels+1)
 
! External functions
REAL, EXTERNAL   :: eg_baro_eta_conv
REAL, EXTERNAL   :: eg_baro_u_p
REAL, EXTERNAL   :: eg_baro_pert

IF ( timestep_number /= 1 ) RETURN

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_IDL_BAROCLINIC',zhook_in,zhook_handle)
 
! Put each column into hydrostatic balance with temperature of
! analytic steady-state

T_1st_guess=T0
p_1st_guess=1.0e5
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    CALL eg_idl_1d_profs(exner_pro, rho(i,j,:), theta(i,j,:),           &
                         Tprofile_number, T_1st_guess, p_1st_guess,     &
                         xi1_p(i), xi2_p(j),                            &
                         intw_rho2w, intw_w2rho,                        &
                         xi3_at_theta(i,j,:), xi3_at_rho(i,j,:),        &
                         z_top_of_model, Earth_radius, dtheta_dz1,      &
                         g_theta(i,j,:), model_levels)
 
    exner(i,j,:)=exner_pro(1:model_levels+1)
 
! Do not assume that eta_surface=1.0 (though it should be)
! DEPENDS ON: eg_baro_eta_conv
    eta_surface = eg_baro_eta_conv(xi1_p(i),xi2_p(j),                   &
                                   xi3_at_theta(i,j,0)-Earth_radius)
    exner_star(i,j) = eta_surface**kappa
  END DO
END DO
 
! Set wind components
DO k=1,model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
! DEPENDS ON: eg_baro_eta_conv
      eta_p_u(k) = eg_baro_eta_conv(xi1_u(i),xi2_p(j),                &
                                    xi3_at_u(i,j,k)-Earth_radius)
! DEPENDS ON: eg_baro_u_p
      u(i,j,k) = eg_baro_u_p(xi1_u(i),xi2_p(j),eta_p_u(k))
! Add perturbation
      IF (L_baro_Perturbed) THEN
! DEPENDS ON: eg_baro_pert
        u(i,j,k) = u(i,j,k) + eg_baro_pert(xi1_u(i),xi2_p(j))
      END IF
      u_adv(i,j,k)=u(i,j,k)
    END DO
  END DO
END DO
 
v=0.0
v_adv=0.0
w=0.0
w_adv=0.0
 
IF (lhook) CALL dr_hook('EG_IDL_BAROCLINIC',zhook_out,zhook_handle)
 
END SUBROUTINE eg_idl_baroclinic
END MODULE eg_idl_baroclinic_mod
