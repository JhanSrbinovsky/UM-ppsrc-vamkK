! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_sl_helmholtz_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_sl_helmholtz(                                              &
                 exner_np1, p_star_np1, u_np1, v_np1, w_np1,                   &
                 etadot_np1, rho_np1, thetav_np1, exner_prime_term,            &
                 m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1,    &
                 Ih,  InnIts, pre_type,l_rel_tol, GCR_Diagnostics,             &
                 R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,L_eliminate_rho,      &
                 R_m_v_d, R_m_cl_d, R_m_cf_d,R_m_r_d, R_m_gr_d, R_m_cf2_d,     &
                 tol, tol_sc_fact,n_rows, row_length, rows, model_levels,      &
                 S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf,S_m_cf2, S_m_r,S_m_gr&
                 ,psi_w_surf, psi_w_lid)

      USE eg_alpha_mod,      ONLY : alpha_u,     alpha_v,    alpha_w,          &
                                    alpha_theta, alpha_rho, alpha_p

      USE um_parvars,        ONLY : offx, offy, halo_i, halo_j, datastart
      USE um_parparams

      USE yomhook,           ONLY : lhook, dr_hook
      USE parkind1,          ONLY : jprb, jpim

      USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,          &
                                    xi3_at_theta=>r_theta_levels,              &
                                    xi3_at_rho=>r_rho_levels, xi3_at_u=>r_at_u,&
                                    xi3_at_v=>r_at_v

      USE timestep_mod,      ONLY : timestep
      USE proc_info_mod,     ONLY : mype=>me,global_row_length,global_rows,    &
                                    model_domain, gc_proc_row_group,           &
                                    gc_proc_col_group,                         &
                                    at_extremity, n_proc, n_procx, n_procy

      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE Field_Types

      USE eg_helm_var_rhs_mod
      USE eg_helm_fixd_rhs_mod
      USE eg_add_pressure_mod
      USE eg_helm_rhs_star_mod
      USE integrity_mod
      USE eg_bicgstab_mod
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE metric_terms_mod

      USE gravity_mod
      USE helmholtz_const_matrix_mod
      USE coriolis_mod

      USE domain_params
      USE eg_parameters_mod, ONLY : total_conv_inner

      IMPLICIT NONE

!
! Description: Solves the nonlinear Helmholtz problem
!
! Method: ENDGame formulation version 3.01
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions


      REAL,    PARAMETER                                 :: sc_err_min = 1.0e-4

      INTEGER,                             INTENT(IN)    :: row_length
      INTEGER,                             INTENT(IN)    :: rows
      INTEGER,                             INTENT(IN)    :: model_levels
      INTEGER,                             INTENT(IN)    :: n_rows
      INTEGER   InnIts, pre_type
      INTEGER, INTENT(IN) :: GCR_Diagnostics
                                 ! Switch controlling diagnostic output.
                                 ! 0 = none
                                 ! 1 = initial and final residuals
                                 ! 2 = all
                                 ! 3 = iteration count processing


      REAL,                                INTENT(IN)    :: Ih
      REAL                                               :: tol, tol_sc_fact


      REAL                                                                     &
      psi_w_surf(row_length,rows),                                             &
      psi_w_lid (row_length,rows)


! Logical flags

       LOGICAL, INTENT(IN) ::  l_rel_tol


! Input arrays from the departure point calculation

       REAL, INTENT(IN) ::                                                     &
        R_u_d(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),          &
        R_w_d(row_length,rows,0:model_levels),                                 &
        R_theta_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
        R_rho_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
        R_m_v_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
        R_m_cl_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
        R_m_cf_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
        R_m_r_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
        R_m_gr_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
        R_m_cf2_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

       REAL, INTENT(INOUT) ::                                                  &
        R_v_d(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels)
 
       LOGICAL, PARAMETER  :: L_accel_convergence = .FALSE.


       LOGICAL, INTENT(IN) :: L_eliminate_rho

! Pressure perturbation

      REAL, INTENT(INOUT)  ::                                                  &
       exner_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

      REAL ::                                                                  &
       exner0(1:row_length,1:rows,model_levels+1)

      REAL, INTENT(INOUT)   ::                                                 &
            p_star_np1(1-offx:row_length+offx,1-offy:rows+offy)

! Fields at next time level

      REAL, INTENT(INOUT)   ::                                                 &
            u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),      &
            v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),    &
            rho_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      REAL, INTENT(INOUT)   ::                                                 &
            m_v_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
           m_cl_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
           m_cf_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
            m_r_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
           m_gr_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
          m_cf2_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
         thetav_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
              w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
         etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Output array for plotting 

      REAL, INTENT(OUT)   ::                                                   &
            exner_prime_term(1-offx:row_length+offx,1-offy:rows+offy,          &
                             model_levels)

! Local arrays used for departure points

      REAL                ::                                                   &
           R_u_a(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),       &
           R_v_a(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),     &
           R_w_a(row_length,rows,0:model_levels),                              &
       R_theta_a(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
         R_rho_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
        R_etadot(row_length,rows,0:model_levels),                              &
           R_p_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels) 

       REAL               ::                                                   &
             RHS(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
              Rn(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

       REAL               ::                                                   &
             R_v_North_fixd(1-offx:row_length+offx,model_levels),              &
             R_v_South_fixd(1-offx:row_length+offx,model_levels),              &
              R_v_North_var(1-offx:row_length+offx,model_levels),              &
              R_v_South_var(1-offx:row_length+offx,model_levels)

      INTEGER                              :: k, iter
      REAL                                 :: del_rho
      REAL                                 :: ex_extrema(2)
      REAL                                 :: init_err, fin_err
      INTEGER                              :: no_its
      INTEGER                              :: istat1, istat2

!     Fast Physics source terms                 
      REAL ::                                                                  &
             S_u( -offx:row_length-1+offx,1-offy:  rows  +offy,model_levels),  &
             S_v(1-offx:row_length  +offx, -offy:n_rows-1+offy,model_levels),  &
             S_w(row_length,rows,0:model_levels),                              &
        S_thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
           S_m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
           S_m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
          S_m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
         S_m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

      REAL exner_residual,exner_residual_tmp

      CHARACTER(len=256)            :: Cmessage
      CHARACTER(len=15)             :: Routine = 'eg_sl_helmholtz'
      INTEGER                       :: ICODE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_SL_HELMHOLTZ',zhook_in,zhook_handle)




      IF (integrity_test) THEN

        CALL check_hash_m(                                            & 
                    R_u_d,           SIZE(R_u_d),           'R_u_d',  & !1
                    R_v_d,           SIZE(R_v_d),           'R_v_d',  &
                    R_w_d,           SIZE(R_w_d),           'R_w_d',  &
                    R_theta_d,       SIZE(R_theta_d),       'R_t_d',  &
                    R_rho_d,         SIZE(R_rho_d),         'R_r_d',  &
                    R_m_v_d,         SIZE(R_m_v_d),         'R_mvd',  &
                    R_m_cl_d,        SIZE(R_m_cl_d),        'Rmcld',  &
                    R_m_cf_d,        SIZE(R_m_cf_d),        'Rmcfd',  &
                    R_m_r_d,         SIZE(R_m_r_d),         'Rmr_d',  &
                    R_m_gr_d,        SIZE(R_m_gr_d),        'Rmgrd',  & !10
                    R_m_cf2_d,       SIZE(R_m_cf2_d),       'Rmc2d',  &
                    exner_np1,       SIZE(exner_np1),       'pinp1',  & !20
                    p_star_np1,      SIZE(p_star_np1),      'p*np1',  &
                    u_np1,           SIZE(u_np1),           'u_np1',  &
                    v_np1,           SIZE(v_np1),           'v_np1',  &
                    rho_np1,         SIZE(rho_np1),         'r_np1',  &
                    m_v_np1,         SIZE(m_v_np1),         'mvnp1',  &
                    m_cl_np1,        SIZE(m_cl_np1),        'mclp1',  &
                    m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',  &
                    m_r_np1,         SIZE(m_r_np1),         'mrnp1',  &
                    m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',  & !30
                    m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',  &
                    thetav_np1,      SIZE(thetav_np1),      'tvnp1',  &
                    w_np1,           SIZE(w_np1),           'w_np1',  &
                    etadot_np1,      SIZE(etadot_np1),      'ednp1',  &
                    S_u,             SIZE(S_u),             's_u__',  &
                    S_v,             SIZE(S_v),             's_v__',  &
                    S_w,             SIZE(S_w),             's_w__',  &
                    S_thetav,        SIZE(S_thetav),        's_the',  &
                    S_m_v,           SIZE(S_m_v),           's_m_v',  &
                    S_m_cl,          SIZE(S_m_cl),          's_mcl',  & !40
                    S_m_cf,          SIZE(S_m_cf),          's_mcf',  &
                    S_m_r,           SIZE(S_m_r),           's_m_r',  &
                    S_m_gr,          SIZE(S_m_gr),          's_mgr',  &
                    S_m_cf2,         SIZE(S_m_cf2),         'smcf2')
        CALL check_hash_m(                                            &
                    exner_ref_pro,   SIZE(exner_ref_pro),   'piref',  &
                    thetav_ref_pro,  SIZE(thetav_ref_pro),  'tvref',  &
                    rho_ref_pro,     SIZE(rho_ref_pro),     'r_ref',  &
                    thetav_ref_eta,  SIZE(thetav_ref_eta),  'tvree',  &
                    rho_ref_eta,     SIZE(rho_ref_eta),     'rreta',  &
                    xi1_p,           SIZE(xi1_p),           'xi1_p',  &
                    xi1_u,           SIZE(xi1_u),           'xi1_u',  &
                    xi2_p,           SIZE(xi2_p),           'xi2_p',  &
                    xi2_v,           SIZE(xi2_v),           'xi2_v',  &
                    phi_at_p,        SIZE(phi_at_p),        'phi_p',  &
                    phi_at_u,        SIZE(phi_at_u),        'phi_u',  &
                    phi_at_v,        SIZE(phi_at_v),        'phi_v',  &
                    phi_at_eta,      SIZE(phi_at_eta),      'phiet',  & !60
                    intw_u2p,        SIZE(intw_u2p),        'iwu2p',  &
                    intw_v2p,        SIZE(intw_v2p),        'iwv2u',  &
                    intw_p2u,        SIZE(intw_p2u),        'iwp2u',  &
                    intw_p2v,        SIZE(intw_p2v),        'iwp2v',  &
                    intw_rho2w,      SIZE(intw_rho2w),      'iwr2w',  &
                    intw_w2rho,      SIZE(intw_w2rho),      'iww2r',  &
                    deta_xi3,        SIZE(deta_xi3),        'dexi3',  &
                    deta_xi3_theta,  SIZE(deta_xi3_theta),  'dex3t',  &
                    deta_xi3_u,      SIZE(deta_xi3_u),      'dex3u',  &
                    deta_xi3_v,      SIZE(deta_xi3_v),      'dex3v',  & !70
                    dxi1_xi3,        SIZE(dxi1_xi3),        'dx1x3',  &
                    dxi2_xi3,        SIZE(dxi2_xi3),        'dx2x3',  &
                    h1_p,            SIZE(h1_p),            'h1_p_',  &
                    h1_xi1_u,        SIZE(h1_xi1_u),        'h1x1u',  &
                    h1_xi2_v,        SIZE(h1_xi2_v),        'h1x2v',  &
                    h1_p_eta,        SIZE(h1_p_eta),        'h1pet',  &
                    h2_p,            SIZE(h2_p),            'h2_p_',  &
                    h2_xi1_u,        SIZE(h2_xi1_u),        'h2x1u',  &
                    h2_xi2_v,        SIZE(h2_xi2_v),        'h2xiv',  &
                    h2_p_eta,        SIZE(h2_p_eta),        'h2pet',  & !80
                    h3_p,            SIZE(h3_p),            'h3_p_',  &
                    h3_xi1_u,        SIZE(h3_xi1_u),        'h3x1u',  &
                    h3_xi2_v,        SIZE(h3_xi2_v),        'h3x2v',  &
                    h3_p_eta,        SIZE(h3_p_eta),        'h3pet',  &
                    f1_star,         SIZE(f1_star),         'f1sta',  &
                    f2_star,         SIZE(f2_star),         'f2sta',  &
                    f3_star,         SIZE(f3_star),         'f3sta',  &
                    f1_comp,         SIZE(f1_comp),         'f1cmp',  &
                    f2_comp,         SIZE(f2_comp),         'f2cmp',  &
                    f3_comp,         SIZE(f3_comp),         'f3cmp',  & !90
                    g_theta,         SIZE(g_theta),         'gthet',  &
                    g_rho,           SIZE(g_rho),           'g_rho',  &
                    Csxi1_p,         SIZE(Csxi1_p),         'cxi1p',  &
                    Csxi1_u,         SIZE(Csxi1_u),         'cxi1u',  &
                    Csxi2_p,         SIZE(Csxi2_p),         'cxi2p',  &
                    Csxi2_v,         SIZE(Csxi2_v),         'cxi2v',  &
                    Snxi1_p,         SIZE(Snxi1_p),         'sxi1p',  &
                    Snxi1_u,         SIZE(Snxi1_u),         'sxi1u',  &
                    Snxi2_p,         SIZE(Snxi2_p),         'sxi2p',  &
                    Snxi2_v,         SIZE(Snxi2_v),         'sxi2v')    !100
      END IF

      ICODE = 0

      IF( pre_type < 0 .or. pre_type > 4 ) THEN
         ICODE = 5
         CMESSAGE = 'Invalid preconditioner type'
         CALL Ereport(Routine,ICODE,CMESSAGE)
      END IF
!    Check the precision of the helmholtz matrix and the tolerance for 
!    the Krylov subspace solver is sufficient.
      IF(real_eg_hlm_kind == real32) THEN
         IF(mype == 0 .AND. PrintStatus >= PrStatus_Diag) THEN
            WRITE(6,'(A50)') 'EG_SL_HELMHOLTZ using single precision matrix'
         END IF
         IF(tol <= 1e-5) THEN
            ! Safety first, we're too close to the bone
            cmessage = 'SP Hlm matrix not sufficient for accuracy of solver'
            icode = 3 
            CALL Ereport(Routine,icode,cmessage)
         END IF
      ELSE
         IF(mype == 0 .AND. PrintStatus >= PrStatus_Diag) THEN
            WRITE(6,'(A50)') 'EG_SL_HELMHOLTZ using double precision matrix'
         END IF
      END IF

!     del_rho = 0.0
!     IF ( .NOT. L_SLICE ) del_rho = 1.0
      del_rho = 1.0

! Adjust R_v_d at North and South boundaries for use in a cyclic LAM

      IF( model_domain == mt_cyclic_lam ) THEN
         IF( at_extremity(Psouth) ) THEN
            R_v_d(:,vdims%j_start,:) = 0.0
         END IF
         IF( at_extremity(Pnorth) ) THEN
            R_v_d(:,vdims%j_end,:)   = 0.0
         END IF
      END IF

! Calculate fixed part of Right-Hand-Side
    rhs = 0.0

      CALL eg_helm_fixd_rhs(Rn,eta_rho_levels,row_length, rows, model_levels,&
                            offx, offy)

      
! Calculate pressure deviation from reference profile

      DO k = 1, model_levels
         exner_prime_term(:,:,k) = ( exner_np1(:,:,k) - exner_ref_pro(:,:,k) )
      END DO

      IF( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
         WRITE(6,*)  "********************************************"
         WRITE(6,*)  "*    Linear solve for Helmholtz problem    *"
         WRITE(6,*)  "*   ====================================   *"
      END IF


      IF( total_conv_inner) InnIts=20

      DO iter = 1, InnIts

! Calculate arrival point quantities

         CALL eg_helm_rhs_star(R_u_a, R_v_a, R_w_a, R_etadot,                  &
                            R_p_a, R_theta_a, R_rho_a,                         &
                            R_m_v_d, R_m_cl_d, R_m_cf_d,                       &
                            R_m_r_d, R_m_gr_d, R_m_cf2_d,                      &
                            u_np1, v_np1, w_np1, exner_np1,                    &
                            p_star_np1, exner_prime_term, rho_np1, thetav_np1, &
                            etadot_np1, m_v_np1, m_cl_np1, m_cf_np1,           &
                            m_r_np1, m_gr_np1, m_cf2_np1, R_w_d, del_rho, Ih,  &
                            timestep, alpha_w,row_length, rows, n_rows,        &
                            model_levels,   &
                            iter, S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf,    &
                            S_m_cf2,S_m_r,S_m_gr,psi_w_surf, psi_w_lid)

! Update RHS of Helmholtz operator to include "star" values

         CALL eg_helm_var_rhs(RHS,Rn,Ih, eta_rho_levels,              &
             exner_prime_term,                                        &
             R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a, R_p_a,          &
             R_etadot,row_length, rows, n_rows, model_levels,         &
             offx, offy)

         IF (integrity_test) CALL update_hash_m(rhs,SIZE(rhs),'rhs__')


! Solve the linear system

         CALL eg_bicgstab(exner_prime_term,RHS,tol,pre_type,          &
                   l_rel_tol, row_length, rows, n_rows, model_levels, &
                   sc_err_min, init_err, fin_err, no_its, .false.,    &
                   ICODE)

         IF( ICODE /= 0 ) THEN
            CMESSAGE='Convergence failure in BiCGstab'
            CALL Ereport(Routine,ICODE,CMESSAGE)
         END IF


         IF( PrintStatus == PrStatus_Diag )                           &
             exner0(:,:,:) =   exner_np1(1:row_length,1:rows,:)

! Calculate time level n+1 quantities

         CALL eg_add_pressure(exner_prime_term, exner_np1,p_star_np1, &
                 R_v_South_var, R_v_North_var,                        &
                 u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,&
                 m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,      &
                 m_cf2_np1,   R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,&
                 R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a,             &
                 R_etadot, R_p_a,                                     &
                 alpha_w, timestep,  l_eliminate_rho,                 &
                 Ih,  offx, offy,n_rows, row_length, rows, model_levels)


         IF( PrintStatus == PrStatus_Diag ) THEN
             exner_residual_tmp = 0.

            DO k=1,model_levels+1

              exner_residual = MAXVAL(ABS((exner0(:,:,k)-             &
                                     exner_np1(1:row_length,1:rows,k))&
                                     /exner0(:,:,k)))

              CALL gc_rmax(1,n_proc,istat2,exner_residual)
 
              exner_residual_tmp = MAX( exner_residual_tmp,           &
                                        exner_residual )

            END DO

            exner_residual = exner_residual_tmp
         END IF

         ex_extrema(1) =  MINVAL(exner_prime_term)
         ex_extrema(2) = -MAXVAL(exner_prime_term)

         CALL gc_rmin(2,n_proc,istat1,ex_extrema)

         IF( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
            WRITE(6,"(A13,I3,A29)") " * Inner iteration ",iter,       &
                                  "                     *"
            WRITE(6,"(A35,I4,A6)")                                    &
                          " * No. Of linear solver iterations ",      &
                                       no_its,"     *"
            WRITE(6,"(A17,E13.6,A16)") " * Initial error ", init_err, &
                                       "               *"
            WRITE(6,"(A17,E13.6,A16)") " *   Final error ", fin_err,  &
                                       "               *"

            WRITE(6,'(A21,E13.6,A11)')  " *   Min exner prime ",      &
                                        ex_extrema(1), "          *"
            WRITE(6,'(A21,E13.6,A11)')  " *   Max exner prime ",      &
                                       -ex_extrema(2), "          *"

         END IF

         IF (integrity_test)                                          &
            CALL update_hash_m(                                       &
                   exner_np1,       SIZE(exner_np1),       'pinp1',   &
                   p_star_np1,      SIZE(p_star_np1),      'p*np1',   &
                   u_np1,           SIZE(u_np1),           'u_np1',   &
                   v_np1,           SIZE(v_np1),           'v_np1',   &
                   rho_np1,         SIZE(rho_np1),         'r_np1',   &
                   m_v_np1,         SIZE(m_v_np1),         'mvnp1',   &
                   m_cl_np1,        SIZE(m_cl_np1),        'mclp1',   &
                   m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',   &
                   m_r_np1,         SIZE(m_r_np1),         'mrnp1',   &
                   m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',   &
                   m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',   &
                   thetav_np1,      SIZE(thetav_np1),      'tvnp1',   &
                   w_np1,           SIZE(w_np1),           'w_np1',   &
                   etadot_np1,      SIZE(etadot_np1),      'ednp1')

      tol = tol/tol_sc_fact
      END DO

      IF( mype == 0 .AND. GCR_Diagnostics > 0) THEN
         WRITE(6,*)  "********************************************"
         WRITE(6,*)  " "
      END IF


      IF (lhook) CALL dr_hook('EG_SL_HELMHOLTZ',zhook_out,zhook_handle)

      END SUBROUTINE eg_SL_Helmholtz
      END MODULE eg_sl_helmholtz_mod
