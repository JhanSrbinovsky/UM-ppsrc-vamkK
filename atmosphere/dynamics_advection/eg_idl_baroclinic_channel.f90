! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_baroclinic_channel_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_baroclinic_channel(                                   &
                      row_length, rows, n_rows, halo_i, halo_j,         &
                      offx, offy, model_levels,                         &
                      xi3_at_theta, xi3_at_rho, xi3_at_u,g_theta,       &
                      u, v, w, u_adv, v_adv, w_adv,                     &
                      theta, rho, exner, exner_star, L_baro_perturbed)
 
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R, kappa, p_zero
USE earth_constants_mod, ONLY: Earth_radius, g, omega
USE eg_idl_baro_mod
USE eg_idl_1d_profs_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE proc_info_mod, ONLY: global_row_length,global_rows
USE timestep_mod, ONLY : timestep_number

IMPLICIT NONE
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for baroclinic wave test in a Cartesian channel.
!   Isothermal case only. Version based on Ullrich & Jablonowski test
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

REAL    :: u_pert, x0, y0
 
! External functions
REAL, EXTERNAL   :: eg_baro_eta_conv
REAL, EXTERNAL   :: eg_baro_u_p

IF ( timestep_number /= 1 ) RETURN

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook('EG_IDL_BAROCLINIC_CHANNEL',zhook_in,zhook_handle)

!set up constants
L_channel = .TRUE.
Ly = glob_xi2_v(global_rows) - glob_xi2_v(0)
Lx = glob_xi1_u(global_row_length) - glob_xi1_u(0)
Lp = 600.0*1000.0
f0 = 2.0*omega*SIN(baro_phi0*pi/180.0)
beta0 = 0.0 ! beta plane not implemented in EG
! beta0 = 2.0*omega/Earth_radius*COS(baro_phi0*pi/180.0)
x0 = 2000.0*1000.0
y0 = 2500.0*1000.0

!WRITE(6,*) '****************************************'
!WRITE(6,*) '***** Channel Baroclinic wave test *****'
!WRITE(6,*) 'Ly   = ',Ly
!WRITE(6,*) 'Lx   = ',Lx
!WRITE(6,*) 'phi0 = ',baro_phi0
!WRITE(6,*) 'f0   = ',f0
!WRITE(6,*) 'beta0= ',beta0
!WRITE(6,*) 'x0   = ',y0
!WRITE(6,*) 'y0   = ',x0
!WRITE(6,*) '****************************************'

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
        u_pert = EXP(-(  (xi1_u(i) - x0) **2                          &
                       + (xi2_p(j) - y0) **2)/Lp**2)
        u(i,j,k) = u(i,j,k) + u_pert
      END IF
      u_adv(i,j,k)=u(i,j,k)
    END DO
  END DO
END DO
 
v=0.0
v_adv=0.0
w=0.0
w_adv=0.0
 
IF (lhook) CALL dr_hook('EG_IDL_BAROCLINIC_CHANNEL',zhook_out,zhook_handle)
 
END SUBROUTINE eg_idl_baroclinic_channel
END MODULE eg_idl_baroclinic_channel_mod

