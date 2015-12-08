! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_set_init_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_set_init(                                              &
                      row_length, rows, n_rows, halo_i, halo_j,          &
                      offx, offy, model_levels, qprofile_number,         &
                      tprofile_number,                                   &
                      t_surface, p_surface,                              &
                      l_baro_inst, l_HeldSuarez, l_solid_body,           &
                      l_deep_baro_inst, T0_E, T0_P, b_const, k_const,    &
                      l_baro_perturbed, l_shallow, l_const_grav,         &
                      l_isothermal, l_rotate_grid,                       &
                      L_initialise_data, L_cartesian,                    &
                      idl_bubble_option, idl_max_num_bubbles,            &
                      idl_bubble_max, idl_bubble_width,                  &
                      idl_bubble_height,                                 &
                      idl_bubble_depth, idl_bubble_xoffset,              &
                      idl_bubble_yoffset,                                &
                      datastart, global_row_length, global_rows,         &
                      grid_np_lon, grid_np_lat,                          &
                      aa_jet_u0, aa_jet_a, aa_jet_m, aa_jet_n,           &
                      ring_height, theta_pert,                           &
                      u, v, w, u_adv, v_adv, w_adv,  rho, exner, p_star, &
                      u_np1, v_np1, w_np1, etadot_np1,                   &
                      thetav_np1, rho_np1, exner_np1,                    &
                      m_v_np1, m_cl_np1, m_cf_np1,                       &
                      m_r_np1, m_gr_np1, m_cf2_np1, p_star_np1,          &
                      z_top_of_model, dtheta_dz1,                        &
                      l_RK_dps, l_dry, alpha_w, ih,                      &
                      etadot, psi_w_surf, psi_w_lid,                     &
                      thetav,m_v, m_cl, m_cf,m_r, m_gr,                  &
                      m_cf2)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE level_heights_mod
USE atmos_constants_mod
USE earth_constants_mod
USE proc_info_mod,         ONLY : mype=>me, model_domain
USE timestep_mod,          ONLY : timestep,timestep_number
USE ereport_mod,           ONLY : ereport
USE Field_Types

USE eg_idl_baroclinic_mod
USE eg_idl_baroclinic_channel_mod
USE integrity_mod
USE eg_idl_rot_solid_body_mod
USE eg_idl_deep_baroclinic_mod
USE eg_idl_deep_baro_mod
USE eg_idl_initialise_bubble_mod
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE atm_fields_bounds_mod
USE coriolis_mod
USE gravity_mod
USE eg_swap_bounds_mod

USE qprofile_mod, ONLY: qp_dry, qp_qsat, qp_namelist_rh, qp_namelist,      &
                        qp_dump
USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv, tp_bv_isoth,&
                        tp_dyn_core, tp_dyn_core_lam, tp_namelist, tp_dump
USE eg_set_adv_winds_mod                      

USE departure_pts_mod, ONLY :                                            &
                      depart_xi1_u, depart_xi2_u, depart_xi3_u,          &
                      depart_xi1_v, depart_xi2_v, depart_xi3_v,          &
                      depart_xi1_w, depart_xi2_w, depart_xi3_w,          &
                      depart_xi1_rho,depart_xi2_rho,depart_xi3_rho,      &
                      reset_dpt_pts

USE atm_step_local, ONLY : first_atmstep_call

IMPLICIT NONE
!
! Description:
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
REAL ::                                                                 &
  etadot(wdims_s%i_start:wdims_s%i_end,                                 &
         wdims_s%j_start:wdims_s%j_end,                                 &
         wdims_s%k_start:wdims_s%k_end),                                &
  m_v   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cl  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_r   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_gr  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf2 (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end)
REAL :: thetav(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)

INTEGER row_length, rows, n_rows, model_levels
INTEGER offx, offy, halo_i, halo_j

REAL ::                                                                 &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows)

INTEGER global_row_length, global_rows, datastart(3)
INTEGER aa_jet_m, aa_jet_n,qprofile_number,tprofile_number

REAL     t_surface, p_surface, aa_jet_u0, aa_jet_a,               &
         grid_np_lon, grid_np_lat, z_top_of_model,                &
         dtheta_dz1(3), ring_height, theta_pert

LOGICAL l_baro_inst, l_HeldSuarez, l_solid_body,                  &
         l_baro_perturbed, l_shallow, l_const_grav,               &
         l_rotate_grid, l_isothermal, L_initialise_data,          &
         L_cartesian
         
! Deep baroclinic options
LOGICAL, INTENT(IN) :: l_deep_baro_inst 
REAL,    INTENT(IN) :: T0_E, T0_P 
INTEGER, INTENT(IN) :: b_const, k_const      

! Bubble idealised options
INTEGER, INTENT(IN) :: idl_max_num_bubbles
INTEGER, INTENT(IN) :: idl_bubble_option(idl_max_num_bubbles)

REAL, INTENT(IN)  :: idl_bubble_max(idl_max_num_bubbles)
                    ! Bubble maximum amplitude (K)
REAL, INTENT(IN)  :: idl_bubble_width(idl_max_num_bubbles)
                    ! Essentially a scaling factor (m)
REAL, INTENT(IN)  :: idl_bubble_depth(idl_max_num_bubbles)
                    ! Currently not used by ENDGame
REAL, INTENT(IN)  :: idl_bubble_height(idl_max_num_bubbles)
                    ! Height of bubble centre (m)
REAL, INTENT(IN)  :: idl_bubble_xoffset(idl_max_num_bubbles)
         ! Bubble x-offset (normalised units: 0.5 = domain centre)
REAL, INTENT(IN)  :: idl_bubble_yoffset(idl_max_num_bubbles)
         ! Bubble y-offset (normalised units: 0.5 = domain centre)

! Arrays

REAL                                                              &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  u_adv(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,         &
        model_levels),                                            &
  v_adv(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,       &
       model_levels),                                             &
  w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
        0:model_levels),                                          &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),  &
  p_star(1-offx:row_length+offx, 1-offy:rows+offy),               &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL  alpha_w, ih

REAL                                                                  &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
  thetav_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
  exner_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),  &
  p_star_np1(1-offx:row_length+offx, 1-offy:rows+offy),               &
  rho_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
  m_v_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  m_cl_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_cf_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_r_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  m_gr_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_cf2_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)


! Local
INTEGER i, j, k, ierr
REAL    u_at_w, v_at_w


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

LOGICAL l_RK_dps  ! Runge-Kutta departure point flag
LOGICAL l_dry     ! dry thermodynamics flag
LOGICAL initialized

REAL Ttmp, ptmp
REAL sat_pc, sat_grad, sat_min, sat_z0, sat_z1

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_SET_INIT',zhook_in,zhook_handle)


firstas: IF( first_atmstep_call ) THEN
!----------------------------------------------------------------------
! Initialise p_star: Will be computed by the solver at timestep>1
!----------------------------------------------------------------------
  IF(tprofile_number.ne.tp_dump .AND. timestep_number == 1) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        p_star(i,j) = (p_surface / p_zero)**kappa
      END DO
    END DO
  END IF

  IF( L_initialise_data ) THEN

     initialized = .FALSE.

     IF (l_baro_inst) THEN
       IF( model_domain == mt_global ) THEN 
         CALL eg_idl_baroclinic(                                  &
                        row_length, rows, n_rows, halo_i, halo_j, &
                        offx, offy, model_levels,                 &
                        r_theta_levels, r_rho_levels, r_at_u,     &
                        g_theta, u, v, w, u_adv, v_adv, w_adv,    &
                        thetav, rho, exner, p_star,               & 
                        L_baro_perturbed)
         initialized = .TRUE.             
       ELSE IF ( model_domain == mt_cyclic_lam ) THEN
           CALL eg_idl_baroclinic_channel(                        &
                        row_length, rows, n_rows, halo_i, halo_j, &
                        offx, offy, model_levels,                 &
                        r_theta_levels, r_rho_levels, r_at_u,     &
                        g_theta, u, v, w, u_adv, v_adv, w_adv,    &
                        thetav, rho, exner, p_star,               & 
                        L_baro_perturbed)       
         initialized = .TRUE.
       ELSE
           ierr = 1
           CALL ereport("eg_idl_set_init", ierr,                  &
                    " Incorrect domain for baroclinic test" )
       END IF                   
     ELSE IF (l_solid_body) THEN
       CALL eg_idl_rot_solid_body(                                &
                      l_shallow, l_rotate_grid,                   &
                      row_length, rows, n_rows, halo_i, halo_j,   &
                      offx, offy, model_levels,                   &
                      r_rho_levels, r_at_u, r_at_v,               &
                      intw_w2rho, intw_rho2w,                     &
                      csxi1_p, csxi1_u, snxi1_p, snxi1_u,         &
                      csxi2_p, snxi2_p,                           &
                      t_surface, p_surface, two_omega,            &
                      grid_np_lon, grid_np_lat,                   &
                      aa_jet_u0, aa_jet_a, aa_jet_m, aa_jet_n,    &
                      f1_comp, f2_comp, f3_comp,                  &
                      u, v, w, u_adv, v_adv, w_adv,               &
                      p_star, thetav, rho, exner)
       initialized = .TRUE.
     ELSE IF (l_deep_baro_inst) THEN

       CALL eg_init_idl_deep_baro()

       CALL eg_idl_deep_baroclinic(                               &
                      row_length, rows, n_rows, halo_i, halo_j,   &
                      offx, offy, model_levels,                   &
                      r_theta_levels, r_rho_levels,               &
                      r_at_u, r_at_v, g_theta,                    &
                      u, v, w, u_adv, v_adv, w_adv,               &
                      thetav, rho, exner, p_star,                 & 
                      L_baro_perturbed,                           &
                      T0_E, T0_P, b_const, k_const, l_shallow,    &
                      l_rotate_grid, grid_np_lon, grid_np_lat,    &
                      f1_comp, f2_comp, f3_comp)                  
       initialized = .TRUE.
     END IF

     IF ( idl_bubble_option(1) > 0 ) THEN
        CALL eg_idl_initialise_bubble(                            &
                   row_length, rows, halo_i, halo_j,              &
                   offx, offy, model_levels,                      &
                   delta_xi1, delta_xi2, base_xi1, base_xi2,      &
                   xi1_p, xi2_p, r_rho_levels, r_theta_levels,    &
                   intw_w2rho, intw_rho2w,                        &
                   mype,                                          &
                   datastart, global_row_length, global_rows,     &
                   idl_max_num_bubbles, idl_bubble_width(1),      &
                   idl_bubble_height(1), idl_bubble_depth(1),     &
                   idl_bubble_xoffset(1), idl_bubble_yoffset(1),  &
                   idl_bubble_max(1), thetav, rho, exner,         &
                   p_star,                                        &
                   t_surface, p_surface, l_cartesian,             &
                   idl_bubble_option(1)                           &
                   )
        IF ( idl_bubble_option(1) == 2 .AND. timestep_number == 1) THEN
          u = 20.0
          v = 20.0
        END IF
        initialized = .TRUE.
     END IF

     ! initialisation required for some diagnostics
     exner(:,:,model_levels+1) = 1.

     IF ( .NOT.initialized  ) THEN
       ierr = 1  
       CALL ereport( 'EG_SET_INIT', ierr,                         &
              'initialisation requested but no option choosen' )
     END IF
  END IF

!----------------------------------------------------------------------
! COMMUNICATION/LATERAL BOUNDARIES
!----------------------------------------------------------------------

 CALL eg_swap_bounds(u,udims_s,fld_type_u,.TRUE.)
 CALL eg_swap_bounds(v,vdims_s,fld_type_v,.TRUE.)
 CALL eg_swap_bounds(w,wdims_s,fld_type_p,.FALSE.)

!----------------------------------------------------------------------
! Initialise etadot: Will be computed by the solver at timestep>1
!----------------------------------------------------------------------
  tsone : IF(timestep_number == 1) THEN
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

          u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +   &
                                     intw_u2p(i,2)*u(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +     &
                                     intw_u2p(i,2)*u(i,j,k) )

          v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +   &
                                     intw_v2p(j,2)*v(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +     &
                                     intw_v2p(j,2)*v(i,j,k) )

          etadot(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -             &
                             u_at_w*dxi1_xi3(i,j,k)/               &
                                           h1_p_eta(i,j,k) -       &
                             v_at_w*dxi2_xi3(i,j,k)/               &
                                           h2_p_eta(i,j,k) ) /     &
                                             deta_xi3_theta(i,j,k)
         END DO
      END DO
    END DO

    etadot(:,:,0) = 0.0
    etadot(:,:,model_levels) = 0.0

    CALL eg_swap_bounds(etadot,wdims_s,fld_type_p,.FALSE.)

     
    k = 0
    DO j = pdims%j_start, pdims%j_end
       DO i = pdims%i_start, pdims%i_end

            psi_w_surf(i,j)   = Ih*(                                           &
                                (h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*             &
                                dxi1_xi3(i,j,k)*(intw_u2p(i,1)*u(i-1,j,k+1)    &
                                                 +intw_u2p(i,2)*u(i,j,k+1)) +  &
                                (h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*             &
                                dxi2_xi3(i,j,k)*(intw_v2p(j,1)*v(i,j-1,k+1)    &
                                                 +intw_v2p(j,2)*v(i,j,k+1)))

            psi_w_surf(i,j) = psi_w_surf(i,j)/(alpha_w*timestep)

            psi_w_lid(i,j) = 0.0
       END DO
    END DO

    IF( .NOT. l_shallow ) THEN
       k = model_levels
       DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end

               psi_w_surf(i,j) = psi_w_surf(i,j)                               &
                                 -( (intw_u2p(i,1)*u(i-1,j,1)                  &
                                    +intw_u2p(i,2)*u(i,j,1) )**2               &
                                   +(intw_v2p(j,1)*v(i,j-1,1)                  &
                                    +intw_v2p(j,2)*v(i,j,1) )**2               &
                                  )/r_theta_levels(i,j,0)

               psi_w_lid(i,j) = psi_w_lid(i,j)                                 &
                                   -( (intw_u2p(i,1)*u(i-1,j,k)                &
                                      +intw_u2p(i,2)*u(i,j,k) )**2             &
                                     +(intw_v2p(j,1)*v(i,j-1,k)                &
                                      +intw_v2p(j,2)*v(i,j,k) )**2             &
                                    )/r_theta_levels(i,j,k)

          END DO
       END DO
    END IF

  END IF tsone

! DEPENDS ON: swap_bounds
  CALL swap_bounds(p_star,                                                     &
       row_length,rows,1,offx,offy,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(thetav,tdims_s,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(exner,tdims_s,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(rho,pdims_s,fld_type_p,.FALSE.)
  
END IF firstas



! Departure points need to be initialised at the first run, i.e.
! when Runge Kutta departure points are used (Polar rows only, technically)
! because then the departure points are not in the chain dump and the
! departure point algorithm will not compute polar v departure points.
! This is not an issue with interpolated departure points, as they are
! derived from the w departure points. Since this is only called once at
! startup of the model and does not take long it does not matter.

! initialise u,v,w departure points

CALL reset_dpt_pts()

first_ts : IF (timestep_number == 1) THEN


  CALL eg_set_adv_winds(u,v,etadot,                                  &
                      u_adv,v_adv,w_adv,row_length,rows,n_rows,      &
                      model_levels, halo_i, halo_j, l_shallow)



  IF(qprofile_number /= qp_dump) THEN
    m_v   = 0.0
    m_cl  = 0.0
    m_cf  = 0.0
    m_r   = 0.0
    m_gr  = 0.0
    m_cf2 = 0.0
  END IF

  IF(qprofile_number == qp_qsat .AND. .NOT. l_dry) THEN
    WRITE(6,'(A)') 'Setting idealised moisture profile'
    sat_z0 = 1000.0 + earth_radius
    sat_z1 = 20000.0 + earth_radius
    sat_pc = 0.95
    sat_min= 0.01
    sat_grad = (sat_pc - sat_min)/(sat_z1-sat_z0)
  
    DO k = 0,model_levels 
      DO j=1,rows
        DO i=1,row_length
!   Compute T on theta levels        
          IF ( k == 0 ) THEN
            Ttmp = thetav(i,j,k)*(p_star(i,j)/p_zero)**kappa
            ptmp = p_star(i,j)
          ELSE
! Compute p on theta levels
            ptmp = p_zero*(intw_rho2w(k,1)*exner(i,j,k+1) +  &
                           intw_rho2w(k,2)*exner(i,j,k))     &
                           **(1.0/kappa)
            Ttmp = thetav(i,j,k)*(intw_rho2w(k,1)*exner(i,j,k+1) +  &
                                  intw_rho2w(k,2)*exner(i,j,k))
          END IF
! DEPENDS ON: qsat
          CALL qsat(m_v(i,j,k),Ttmp,ptmp,1)         
          IF ( r_theta_levels(i,j,k) < sat_z0 ) THEN
            m_v(i,j,k) = sat_pc*m_v(i,j,k)
          ELSE IF ( r_theta_levels(i,j,k) < sat_z1 ) THEN            
            m_v(i,j,k) = (sat_pc-sat_grad*(r_theta_levels(i,j,k)-sat_z0)) &
                         *m_v(i,j,k)
          ELSE        
            m_v(i,j,k) = 0.0
          END IF
        END DO
      END DO
    END DO
! DEPENDS ON: swap_bounds
    CALL eg_swap_bounds(m_v,tdims_s,fld_type_p,.FALSE.)
    WRITE(6,'(A)') 'Moisture profile of  first point'
    i=1
    j=1
    DO k=1,model_levels
      WRITE(6,fmt='(I4,2E16.8)') k,r_theta_levels(i,j,k),m_v(i,j,k)
    END DO
    WRITE(6,'(A)') ' '
  END IF
END IF first_ts

CALL eg_swap_bounds(u_adv,udims_l,fld_type_u,.TRUE.)
CALL eg_swap_bounds(v_adv,vdims_l,fld_type_v,.TRUE.)
CALL eg_swap_bounds(w_adv,wdims_l,fld_type_p,.FALSE.)

u_np1       = u
v_np1       = v
w_np1       = w
thetav_np1  = thetav
rho_np1     = rho
etadot_np1  = etadot
m_v_np1     = m_v
m_cl_np1    = m_cl
m_cf_np1    = m_cf
m_r_np1     = m_r
m_gr_np1    = m_gr
m_cf2_np1   = m_cf2
exner_np1(:,:,:) = exner(:,:,:)
p_star_np1   = p_star


IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                    exner,           SIZE(exner),           'pi___',  &
                    exner_np1,       SIZE(exner_np1),       'pinp1',  &
                    p_star_np1,      SIZE(p_star_np1),      'p*np1',  &
                    u,               SIZE(u),               'u____',  &
                    v,               SIZE(v),               'v____',  &
                    u_np1,           SIZE(u_np1),           'u_np1',  &
                    v_np1,           SIZE(v_np1),           'v_np1',  & 
                    psi_w_surf,      SIZE(psi_w_surf),      'psiws',  &
                    rho_np1,         SIZE(rho_np1),         'r_np1',  &
                    m_v_np1,         SIZE(m_v_np1),         'mvnp1',  &
                    m_cl_np1,        SIZE(m_cl_np1),        'mclp1',  &
                    m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',  &
                    m_r_np1,         SIZE(m_r_np1),         'mrnp1',  &
                    m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',  &
                    m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',  &
                    thetav_np1,      SIZE(thetav_np1),      'tvnp1',  &
                    w_np1,           SIZE(w_np1),           'w_np1',  &
                    etadot,          SIZE(etadot),          'ed___',  &
                    etadot_np1,      SIZE(etadot_np1),      'ednp1')



IF (lhook) CALL dr_hook('EG_IDL_SET_INIT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_set_init
END MODULE eg_idl_set_init_mod
