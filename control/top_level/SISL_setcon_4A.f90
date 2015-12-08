! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_setcon_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sisl_setcon(                                            &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,      &
         offx, offy,  l_datastart, l_slice,                           &
         l_test_tracer, l_cartesian, l_shallow, l_const_grav,         &
         l_impl_horz_drag,l_expl_horz_drag,eccentricity,              &
         height_domain, alpha_u, alpha_v, alpha_w,alpha_theta,        &
         alpha_rho, tprofile_number, trefer_number,dtheta_dz1,        &
         z_top_of_model,g_theta,g_rho, t_surf_ref_in, p_surf_ref_in,  &
         eg_vert_damp_coeff, eg_vert_damp_profile, eta_s,             &
         thetav, rho, t_surface_in, p_surface_in, exner,              &
         alpha_p, l_eliminate_rho,                                    &
         ih, pole_consts)

USE atmos_constants_mod
USE conversions_mod
USE earth_constants_mod
USE init_vert_damp_mod,    ONLY : init_vert_damp
USE eg_init_horz_drag_mod, ONLY : eg_init_horz_drag
USE parkind1,              ONLY : jpim, jprb       !DrHook
USE yomhook,               ONLY : lhook, dr_hook   !DrHook

USE proc_info_mod,         ONLY : me,                                 &
                                  global_row_length,                  &
                                  global_rows,                        &
                                  model_domain,                       &
                                  gc_proc_row_group,                  &
                                  at_extremity
USE timestep_mod,          ONLY : timestep,                           &
                                  timestep_number

USE level_heights_mod,     ONLY : eta_theta_levels,                   &
                                  eta_rho_levels,                     &
                                  xi3_at_theta=>r_theta_levels,       &
                                  xi3_at_rho=>r_rho_levels,           &
                                  xi3_at_u=>r_at_u,                   &
                                  xi3_at_v=>r_at_v,                   &
                                  xi3_at_u_w=>r_at_u_w,               &
                                  xi3_at_v_w=>r_at_v_w
USE eg_helmholtz_mod
USE eg_idl_1d_profs_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE set_metric_terms_4A_mod
USE metric_terms_mod
USE ereport_mod, ONLY : ereport
USE field_types
USE UM_ParParams
USE PrintStatus_mod
USE atm_fields_bounds_mod
USE eg_swap_bounds_mod
USE set_vert_interp_consts_mod
USE set_horiz_interp_consts_mod

USE lookup_table_mod

USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv,        &
                        tp_bv_isoth, tp_dyn_core, tp_dyn_core_lam,    &
                        tp_namelist, tp_dump
IMPLICIT NONE

REAL    eg_vert_damp_coeff, eta_s
INTEGER eg_vert_damp_profile

!
! Description: computes the following constants for the SISL scheme:
!              metric terms, reference profiles, Helmholtz coefficients.
!  
!
! Method:: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,    &
                       l_datastart(3), halo_i, halo_j, offx,      &
                       offy,  tprofile_number,                    &
                       trefer_number
!
LOGICAL, INTENT(IN) :: l_slice, l_test_tracer,                    &
                       l_cartesian, l_shallow, l_const_grav
!
! SI time weights

REAL, INTENT(IN) :: alpha_u, alpha_v, alpha_w, alpha_theta,       &
                    alpha_rho, ih, alpha_p

LOGICAL, INTENT(IN) :: l_eliminate_rho

REAL, INTENT(IN) :: dtheta_dz1(3)

REAL, INTENT(IN) ::  z_top_of_model

! Surface temperature and pressure for refernce profile
REAL, INTENT(IN) :: t_surf_ref_in, p_surf_ref_in

REAL  t_surf_ref2d(1-offx:row_length+offx,                         &
                   1-offy:rows+offy)                               &
,     p_surf_ref2d(1-offx:row_length+offx,                         &
                   1-offy:rows+offy)

! Surface temperature and pressure for initial data
REAL, INTENT(IN) :: t_surface_in, p_surface_in

REAL  t_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)                                   &
,     p_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)

LOGICAL, INTENT(IN)  :: l_impl_horz_drag
LOGICAL, INTENT(IN)  :: l_expl_horz_drag
REAL,    INTENT(IN)  :: eccentricity
REAL,    INTENT(IN)  :: height_domain

REAL, INTENT(INOUT) ::                                            &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

!Output Arrays from this routine
REAL, INTENT(INOUT) ::                                            &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

REAL, INTENT(OUT) :: pole_consts(4)

! Gravity arrays

REAL, INTENT(OUT) ::                                              &
  g_theta (1-offx:row_length+offx, 1-offy:rows+offy,              &
           0:model_levels)                                        &
, g_rho (1-offx:row_length+offx, 1-offy:rows+offy,                &
           model_levels)

! Local variables

INTEGER :: i,j,k, option_bubble, info

! Reference gravity acceleration

REAL :: g_theta_ref (0:model_levels),                             &
        g_rho_ref (model_levels),                                 &
        r_ref_theta(0:model_levels),                              &
        r_ref_rho(model_levels)

REAL :: dtheta_dz1_ref(3)

! Tolerance for initial and reference surface pressures to be considered
! to differ:
REAL, PARAMETER :: p_tol=0.1

! Tolerance for initial and reference surface temperatures to be
! considered to differ:
REAL, PARAMETER :: t_tol=0.01
REAL, PARAMETER :: delta = 1.0e-8
REAL dx, trig_chck, c, d, e, f, tmp_pole(0:row_length-1,4)

! Deep atmosphere
LOGICAL         :: l_deep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SISL_SETCON',zhook_in,zhook_handle)

! set 3D gravity acceleration arrays

!      If( L_Cartesian ) Then
 IF(l_cartesian .OR. l_shallow .OR. l_const_grav) THEN
   l_deep=.FALSE.
   IF (PrintStatus >= PrStatus_Normal) THEN
     ierr = -1
     CALL Ereport ('eg_SISL_setcon',ierr,                             &
                        ' Constant gravity enforced ')
   END IF
   k = 0
   DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
         g_theta(i,j,k)=g
      END DO
   END DO

   DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
         DO i=pdims%i_start, pdims%i_end
            g_theta(i,j,k)=g
            g_rho(i,j,k)  =g
         END DO
      END DO
   END DO
ELSE
   l_deep=.TRUE.
   k = 0
   DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
         g_theta(i,j,k)=g*(earth_radius/xi3_at_theta(i,j,k))**2
      END DO
   END DO

   DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
         DO i=pdims%i_start, pdims%i_end
            g_theta(i,j,k)=g*(earth_radius/xi3_at_theta(i,j,k))**2
            g_rho(i,j,k)  =g*(earth_radius/xi3_at_rho(i,j,k))**2
         END DO
      END DO
   END DO
END IF

CALL eg_swap_bounds(g_theta,tdims_s,fld_type_p,.FALSE.)

CALL eg_swap_bounds(g_rho,pdims_s,fld_type_p,.FALSE.)


!-------------------------------------------------------------------------
! Compute reference heights and gravitational acceleration profiles
!-------------------------------------------------------------------------
r_ref_theta(0) = 0.0
DO k=1, model_levels
    r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
    r_ref_rho(k)   = eta_rho_levels(k) * z_top_of_model
END DO

! set 1D gravity acceleration reference array independent of orography
g_theta_ref(0) = g
IF( l_deep ) THEN
  DO k=1, model_levels
    g_theta_ref(k) = g  * ( earth_radius /                        &
                          ( earth_radius + r_ref_theta(k) ) )**2
    g_rho_ref(k)   = g  * ( earth_radius /                        &
                          ( earth_radius + r_ref_rho(k) ) )**2
  END DO
ELSE
  DO k=1, model_levels
    g_theta_ref(k) = g
    g_rho_ref(k)   = g
  END DO
END IF

!-------------------------------------------------------------------------
! Compute metric terms
!-------------------------------------------------------------------------

CALL set_metric_terms_4A (l_cartesian, l_shallow, eccentricity)

! Calculate constants used for v_at_poles (only done at the
! begining of the run so just use REPROD version

pole_consts = 0.0
tmp_pole    = 0.0
IF (.NOT. l_cartesian .AND. model_domain == mt_global) THEN
  IF( at_extremity(pnorth) .OR. at_extremity(psouth) ) THEN
    DO i = udims%i_start, udims%i_end
      dx = ( xi1_p(i+1) - xi1_p(i) )/(2.0*pi)
      tmp_pole(i,1) = dx*COS(2.0*xi1_u(i))
      tmp_pole(i,2) = dx*SIN(2.0*xi1_u(i))
      tmp_pole(i,3) = dx*SIN(xi1_u(i))
      tmp_pole(i,4) = dx*COS(xi1_u(i))
      END DO
    CALL gcg_rvecsumr(row_length,row_length,1,4,tmp_pole,          &
                      gc_proc_row_group,info,pole_consts)
   c = pole_consts(1)
   d = pole_consts(2)
   e = pole_consts(3)
   f = pole_consts(4)
   trig_chck = SQRT((c + e**2 - f**2)**2 + (d-2.0*e*f)**2) + e**2 + f**2
   IF( trig_chck > 1.0 - delta ) THEN
      CALL Ereport ('eg_SISL_setcon',1,' Problem at the poles ')
   END IF
  END IF
END IF

!-------------------------------------------------------------------------
! Set horizontal linear interpolation weights
!-------------------------------------------------------------------------

DO i=0, row_length
  intw_p2u(i,1) = (xi1_p(i+1)-xi1_u(i)) / (xi1_p(i+1)-xi1_p(i))
  intw_p2u(i,2) = 1.0 - intw_p2u(i,1)
END DO

DO j=0, n_rows
  intw_p2v(j,1) = (xi2_p(j+1)-xi2_v(j)) / (xi2_p(j+1)-xi2_p(j))
  intw_p2v(j,2) = 1.0 - intw_p2v(j,1)
END DO

DO i=1, row_length
  intw_u2p(i,1) = (xi1_u(i)-xi1_p(i)) / (xi1_u(i)-xi1_u(i-1))
  intw_u2p(i,2) = 1.0 -  intw_u2p(i,1)
END DO

DO j=1, rows
  intw_v2p(j,1) = (xi2_v(j)-xi2_p(j)) / (xi2_v(j)-xi2_v(j-1))
  intw_v2p(j,2) = 1.0 -intw_v2p(j,1)
END DO

!-------------------------------------------------------------------------
! Set vertical linear interpolation weights (eq C.10 ENDGame
! formulation vn 1.01)
!-------------------------------------------------------------------------
DO k=1, model_levels
  intw_w2rho(k,1) = ( eta_rho_levels(k)-eta_theta_levels(k-1) ) / &
                    ( eta_theta_levels(k)-eta_theta_levels(k-1) )
  intw_w2rho(k,2) = 1.0 - intw_w2rho(k,1)
END DO

DO k=1, model_levels-1
  intw_rho2w(k,1) = ( eta_theta_levels(k)-eta_rho_levels(k) ) /   &
                    ( eta_rho_levels(k+1)-eta_rho_levels(k) )
  intw_rho2w(k,2) = 1.0 - intw_rho2w(k,1)
END DO

intw_rho2w(model_levels,1) = 0.5
intw_rho2w(model_levels,2) = 0.5  ! as in New dynamics

!-------------------------------------------------------------------------
! Initialise vertical damping
!-------------------------------------------------------------------------
CALL init_vert_damp(eg_vert_damp_profile, eta_s,                  &
                    eg_vert_damp_coeff,height_domain)

CALL eg_init_horz_drag(l_impl_horz_drag,l_expl_horz_drag)

! Set constants used by the interpolation routines

CALL set_vert_interp_consts(eta_rho_levels, eta_theta_levels, model_levels)
CALL set_horiz_interp_consts()
CALL set_look_up_table()

!-------------------------------------------------------------------------
! Compute Helmholtz coefficients
!-------------------------------------------------------------------------

CALL eg_swap_bounds(h1_p,pdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h2_p,pdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h3_p,pdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(phi_at_p,pdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h1_p_eta,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h2_p_eta,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h3_p_eta,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(phi_at_eta,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(h1_xi1_u,udims_s,fld_type_u,.FALSE.)
CALL eg_swap_bounds(h2_xi1_u,udims_s,fld_type_u,.FALSE.)
CALL eg_swap_bounds(h3_xi1_u,udims_s,fld_type_u,.FALSE.)
CALL eg_swap_bounds(phi_at_u,udims_s,fld_type_u,.FALSE.)
CALL eg_swap_bounds(h1_xi2_v,vdims_s,fld_type_v,.FALSE.)
CALL eg_swap_bounds(h2_xi2_v,vdims_s,fld_type_v,.FALSE.)
CALL eg_swap_bounds(h3_xi2_v,vdims_s,fld_type_v,.FALSE.)
CALL eg_swap_bounds(phi_at_v,vdims_s,fld_type_v,.FALSE.)
CALL eg_swap_bounds(deta_xi3,pdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(deta_xi3_u,udims_s,fld_type_u,.FALSE.)
CALL eg_swap_bounds(deta_xi3_v,vdims_s,fld_type_v,.FALSE.)
CALL eg_swap_bounds(deta_xi3_theta,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(dxi1_xi3,tdims_s,fld_type_p,.FALSE.)
CALL eg_swap_bounds(dxi2_xi3,tdims_s,fld_type_p,.FALSE.)


IF (integrity_test) THEN
  CALL update_hash_m(                                                 &
                    phi_at_p,        SIZE(phi_at_p),        'phi_p',  &
                    phi_at_u,        SIZE(phi_at_u),        'phi_u',  &
                    phi_at_v,        SIZE(phi_at_v),        'phi_v',  &
                    phi_at_eta,      SIZE(phi_at_eta),      'phiet',  &
                    h1_p,            SIZE(h1_p),            'h1_p_',  &
                    h1_xi1_u,        SIZE(h1_xi1_u),        'h1x1u',  &
                    h1_xi2_v,        SIZE(h1_xi2_v),        'h1x2v',  &
                    h1_p_eta,        SIZE(h1_p_eta),        'h1pet',  &
                    h2_p,            SIZE(h2_p),            'h2_p_',  &
                    h2_xi1_u,        SIZE(h2_xi1_u),        'h2x1u',  &
                    h2_xi2_v,        SIZE(h2_xi2_v),        'h2xiv',  &
                    h2_p_eta,        SIZE(h2_p_eta),        'h2pet',  &
                    h3_p,            SIZE(h3_p),            'h3_p_',  &
                    h3_xi1_u,        SIZE(h3_xi1_u),        'h3x1u',  &
                    h3_xi2_v,        SIZE(h3_xi2_v),        'h3x2v',  &
                    h3_p_eta,        SIZE(h3_p_eta),        'h3pet',  &
                    intw_u2p,        SIZE(intw_u2p),        'iwu2p',  &
                    intw_v2p,        SIZE(intw_v2p),        'iwv2u',  &
                    intw_p2u,        SIZE(intw_p2u),        'iwp2u',  &
                    intw_p2v,        SIZE(intw_p2v),        'iwp2v',  &
                    intw_rho2w,      SIZE(intw_rho2w),      'iwr2w',  &
                    intw_w2rho,      SIZE(intw_w2rho),      'iww2r',  &
                    deta_xi3,        SIZE(deta_xi3),        'dexi3',  &
                    deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t',  &
                    deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u',  &
                    deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v',  &
                    dxi1_xi3,        SIZE(dxi1_xi3),        'dx1x3',  &
                    dxi2_xi3,        SIZE(dxi2_xi3),        'dx2x3',  &
                    g_theta,         SIZE(g_theta),         'gthet',  &
                    g_rho,           SIZE(g_rho),           'g_rho')
END IF

IF (lhook) CALL dr_hook('EG_SISL_SETCON',zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_setcon
END MODULE eg_sisl_setcon_mod
