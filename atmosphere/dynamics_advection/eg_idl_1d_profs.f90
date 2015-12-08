! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_1d_profs_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_1d_profs(exner_ref_pro, rho_ref_pro, thetav_ref_pro,  &
                           Tprofile_number, T_surface, p_surface,       &
                           xi1_p, xi2_p, intw_rho2w, intw_w2rho,        &
                           coord_at_theta, coord_at_rho,                &
                           depth_factor, surface_offset,                &
                           dtheta_dz1, g_theta_ref, model_levels)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R, p_zero, kappa, recip_kappa

USE ereport_mod, ONLY : ereport

USE eg_idl_deep_baro_mod, ONLY: tau1, tau2, tau1_int, tau2_int,         &
                                k_const_save, L_shallow_save, T0, H, B, &
                                C, b_const_save, eg_lapse
USE earth_constants_mod, ONLY: Earth_radius, g, omega
IMPLICIT NONE
!
! Description:
!  Calculate 1D vertical profiles for exner, thetav and rho.
!  This is basically a wrapper for eg_idl_Newton.
!  
!
! Method:
!
!  Solve discrete equations for hydrostatic balance and equation of
!  state, for given temperature or potential temperature profile.
!  equations solved using Newton-Raphson iteration.
!
!  Returns potential temperature, exner pressure and density.
!
!  Generates thermodynamic profiles at height levels given by
!
!         z = (coord - surface_offset) * depth_factor 
!
!  for coord = coord_at_theta and coord_at_rho.
!
!  Used in two ways:
!
!  (i)  Generate semi-implicit reference profiles with:
!       coord_at_theta = eta_theta_levels
!       coord_at_rho   = eta_rho_levels
!       depth_factor   = height_domain
!       surface_offset = 0.0
!
!  (ii) Generate profiles on a column of the grid:
!       coord_at_theta = xi3_at_theta
!       coord_at_rho   = xi3_at_rho
!       depth_factor   = 1.0
!       surface_offset = Earth_radius
! 
!  (This should be simpified at some point!)

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


INTEGER,INTENT(IN)    :: model_levels
INTEGER,INTENT(IN)    :: Tprofile_number
REAL,   INTENT(INOUT) :: T_surface         ! Surface temperature (K)
REAL,   INTENT(INOUT) :: p_surface         ! Surface pressure (Pa)
REAL,   INTENT(IN)    :: depth_factor
REAL,   INTENT(IN)    :: surface_offset
REAL,   INTENT(IN)    :: xi1_p
REAL,   INTENT(IN)    :: xi2_p
REAL,   INTENT(IN)    :: intw_rho2w(model_levels,2)
REAL,   INTENT(IN)    :: intw_w2rho(model_levels,2)
REAL,   INTENT(IN)    :: coord_at_theta(0:model_levels)
REAL,   INTENT(IN)    :: coord_at_rho(model_levels)
REAL,   INTENT(IN)    :: g_theta_ref(0:model_levels)
REAL,   INTENT(IN)    :: dtheta_dz1(3)

! Variables for output
REAL,   INTENT(OUT)   :: exner_ref_pro(0:model_levels+1) 
REAL,   INTENT(OUT)   :: thetav_ref_pro(0:model_levels)
REAL,   INTENT(OUT)   :: rho_ref_pro(model_levels)

! Local variables
INTEGER               :: k
INTEGER               :: prof_type
REAL                  :: eta(0:model_levels)
REAL                  :: geopotential_factor(0:model_levels)
REAL                  :: levs(0:model_levels)
REAL                  :: T0_ref(0:model_levels)
REAL                  :: u_term(1:model_levels)
REAL                  :: u_temp
REAL                  :: grad

! For Deep baroclinic wave
REAL                  :: r_height, tau2_term

! External functions
REAL, EXTERNAL        :: eg_baro_eta_conv
REAL, EXTERNAL        :: eg_baro_T_p

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_1D_PROFS',zhook_in,zhook_handle)

!
! Set initial guess
!
exner_ref_pro(0)  = (p_surface / p_zero)**kappa
exner_ref_pro(:)  = exner_ref_pro(0)
thetav_ref_pro(:) = T_surface / exner_ref_pro(0)
rho_ref_pro(:)    = 1.0
T0_ref(:)         = T_surface

! When using eta-levels for generating reference profile, the 
! height_domain factor is carried with the gravitational acceleration.
geopotential_factor(:) = depth_factor*g_theta_ref(:)

! Put level data into levs array

! For atmospheric state these are z levels in metres
! For reference state these are eta levels

! Make levels
DO k = 1, model_levels
  levs(k) = coord_at_rho(k) - surface_offset
END DO
levs(0) = coord_at_theta(0) - surface_offset

u_term = 0.0

! Check if Tprofile_number a valid option and set up T0_ref
! if required.

SELECT CASE(Tprofile_number)

  CASE(1)             ! T0_ref is a theta point

    prof_type = 1     ! theta specified
    T0_ref(0) = T_surface

    DO k = 1, model_levels
      T0_ref(k) = T0_ref(k-1) + dtheta_dz1(1) * depth_factor *          &
                                (levs(k)-levs(k-1))
    END DO


  CASE(2)             ! T0_ref is a rho point

    prof_type = 2     ! T specified
    T0_ref(0) = T_surface
    DO k = 1, model_levels
      T0_ref(k) = T0_ref(k-1) - dtheta_dz1(1) * depth_factor *          &
                                (levs(k)-levs(k-1))
    END DO

  CASE(3) ! Jablonowski & Williamson 2006 baroclinic wave test

    prof_type = 2   ! T specified
    DO k = 0, model_levels
!DEPENDS ON: eg_baro_eta_conv
      eta(k) = eg_baro_eta_conv(xi1_p, xi2_p, depth_factor*levs(k))
!DEPENDS ON: eg_baro_T_p
      T0_ref(k) = eg_baro_T_p(xi1_p, xi2_p, eta(k))
! A reasonable first guess:
      exner_ref_pro(k) = (p_zero*eta(k))**kappa
      IF (k >= 1) THEN
        rho_ref_pro(k) = p_zero*eta(k) / (R*T0_ref(k))
      END IF
      thetav_ref_pro(k) = T0_ref(k)/exner_ref_pro(k)
    END DO
! A reasonable first guess:
    exner_ref_pro(model_levels+1)=sqrt(exner_ref_pro(model_levels)**2 / &
                                       exner_ref_pro(model_levels-1))

! Return actual surface temperature and pressure
    p_surface=p_zero*eta(0)  ! Boundary condition for Newton solver
    T_surface=T0_ref(0)

  CASE(4) ! Staniforth Deep atmosphere baroclinic wave test

    prof_type = 2   ! T specified
    DO k = 0, model_levels
      IF ( L_shallow_save ) THEN     
        r_height = Earth_radius
      ELSE
        IF ( k == 0) THEN
          r_height = surface_offset
        ELSE
          r_height = coord_at_theta(k)
          r_height = coord_at_rho(k)
        END IF
      END IF
      tau2_term = (r_height/Earth_radius*COS(xi2_p)) ** k_const_save       &
                - (k_const_save/(k_const_save + 2.0))*                     &
                  (r_height/Earth_radius*COS(xi2_p)) **(k_const_save + 2.0)

      T0_ref(k) = (Earth_radius/r_height)**2 / (tau1(k) - tau2(k)*tau2_term)

      exner_ref_pro(k) = (EXP(-g/R*(tau1_int(k)-tau2_int(k)*tau2_term)))**kappa

      thetav_ref_pro(k) = T0_ref(k)/exner_ref_pro(k)

      IF (k >= 1) THEN
        rho_ref_pro(k) = p_zero*exner_ref_pro(k)**(recip_kappa)/(R*T0_ref(k))

      IF ( .NOT. L_shallow_save ) THEN
! U coriolis term
! This is the wind proxy: U
          u_temp = (g/Earth_radius)*k_const_save*tau2_int(k)                 &
                    *((r_height/Earth_radius*COS(xi2_p))**(k_const_save-1)   &
                    - (r_height/Earth_radius*COS(xi2_p))**(k_const_save+1))  &
                     *T0_ref(k)

          u_temp  = - omega*r_height*COS(xi2_p)                             &
                     + SQRT( (omega*r_height*COS(xi2_p))**2                 &
                            + r_height*COS(xi2_p)*u_temp )

          u_term(k) = u_temp**2/r_height+2.0*omega*u_temp*COS(xi2_p)
        END IF
      END IF

    END DO

! A better guess    
    k = model_levels
    exner_ref_pro(k+1)=exner_ref_pro(k)                                      &
                           + 2.0*(coord_at_theta(k)-coord_at_rho(k))         &
                      /((R/kappa)*thetav_ref_pro(k))*(u_term(k)-g_theta_ref(k))

! Return actual surface temperature and pressure
    p_surface=(p_zero*exner_ref_pro(0)**(recip_kappa))
    T_surface=T0_ref(0)
    
  CASE(5)
  
    prof_type = 2     ! T specified
    T0_ref(0) = T_surface

    DO k = 1, model_levels      
      grad = -6.0/1000.0
      IF( coord_at_rho(k)-surface_offset > 13000.0 ) grad = -2.0/1000.0
      IF( coord_at_rho(k)-surface_offset > 15000.0 ) grad = 0.0
      IF( coord_at_rho(k)-surface_offset > 16000.0 ) grad = 5.0/8000.0
      IF( coord_at_rho(k)-surface_offset > 25000.0 ) grad = 0.0
      T0_ref(k) = T0_ref(k-1) + grad * depth_factor *          &
                                (levs(k)-levs(k-1))
    END DO

  CASE DEFAULT

    Call ereport("eg_idl_1d_profs", Tprofile_number,                    &
                 "Invalid profile type" )

END SELECT

! DEPENDS ON: eg_Newton
  CALL eg_newton(exner_ref_pro, thetav_ref_pro, rho_ref_pro, T0_ref,      &
                 geopotential_factor, levs, intw_rho2w, intw_w2rho,       &
                 p_surface, model_levels, prof_type, u_term)

IF (lhook) CALL dr_hook('EG_IDL_1D_PROFS',zhook_out,zhook_handle)

END SUBROUTINE eg_idl_1d_profs
END MODULE eg_idl_1d_profs_mod
