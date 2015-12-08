! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_sl_helmholtz_inc_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_SL_Helmholtz_inc(                                 &
                 exner_np1, p_star_np1, u_np1, v_np1, w_np1,          &
                 etadot_np1, rho_np1, thetav_np1,                     &
                 exner_prime_term, m_v_np1, m_cl_np1, m_cf_np1,       &
                 m_r_np1,m_gr_np1, m_cf2_np1, Ih, InnIts, pre_type,   &
                 l_rel_tol,GCR_Diagnostics,                           &
                 R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,             &
                 R_m_v_d, R_m_cl_d, R_m_cf_d,                         &
                 R_m_r_d, R_m_gr_d, R_m_cf2_d, tol, tol_sc_fact,      &
                 n_rows, row_length,                                  &
                 rows, model_levels, S_u, S_v, S_w, S_thetav,         &
                 S_m_v, S_m_cl, S_m_cf, S_m_cf2, S_m_r, S_m_gr        &
                 ,psi_w_surf, psi_w_lid)

      USE proc_info_mod,     ONLY : datastart=>l_datastart
      USE um_parvars,        ONLY : offx, offy, halo_i, halo_j
      USE eg_alpha_mod,      ONLY : alpha_u,     alpha_v,    alpha_w, &
                                    alpha_theta, alpha_rho, alpha_p
 
      USE yomhook,           ONLY : lhook, dr_hook
      USE parkind1,          ONLY : jprb, jpim

      USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels, &
                                    xi3_at_theta=>r_theta_levels,     &
                                    xi3_at_rho=>r_rho_levels,         &
                                    xi3_at_u=>r_at_u,                 &
                                    xi3_at_v=>r_at_v

      USE timestep_mod,      ONLY : timestep
      USE proc_info_mod,     ONLY : mype=>me,global_row_length,       &
                                    global_rows,                      &
                                    model_domain, gc_proc_row_group,  &
                                    gc_proc_col_group,                &
                                    at_extremity, n_proc, n_procx,    &
                                    n_procy

      USE ereport_mod,       ONLY : ereport
      USE PrintStatus_mod
      USE Field_Types

      USE eg_add_pressure_inc_mod
      USE eg_helm_rhs_inc_mod
      USE integrity_mod
      USE eg_bicgstab_mod
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE metric_terms_mod

      USE gravity_mod
      USE helmholtz_const_matrix_mod
      USE coriolis_mod  

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

      REAL,    PARAMETER                :: sc_err_min = 1.0e-20

      INTEGER,            INTENT(IN)    :: row_length
      INTEGER,            INTENT(IN)    :: rows
      INTEGER,            INTENT(IN)    :: model_levels
      INTEGER,            INTENT(IN)    :: n_rows
      INTEGER,            INTENT(IN)    :: InnIts, pre_type
      REAL,               INTENT(IN)    :: Ih
      REAL                              :: tol, tol_sc_fact

      REAL                                                                     &
      psi_w_surf(row_length,rows),                                             &
      psi_w_lid (row_length,rows)

      INTEGER, INTENT(IN) :: GCR_Diagnostics
                                 ! Switch controlling diagnostic output.
                                 ! 0 = none
                                 ! 1 = initial and final residuals
                                 ! 2 = all
                                 ! 3 = iteration count processing

! Logical flags

       LOGICAL, INTENT(IN) ::  l_rel_tol


! Input arrays from the departure point calculation

      REAL, INTENT(IN) ::    R_u_d(udims_s%i_start:udims_s%i_end,     &
                                   udims_s%j_start:udims_s%j_end,     &
                                   udims_s%k_start:udims_s%k_end)

      REAL, INTENT(IN) ::      R_w_d(wdims%i_start:wdims%i_end,       &
                                     wdims%j_start:wdims%j_end,       &
                                     wdims%k_start:wdims%k_end)

      REAL, INTENT(IN) ::     R_theta_d(tdims_s%i_start:tdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_rho_d(pdims_s%i_start:tdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)


      REAL, INTENT(IN) ::   R_m_v_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cl_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cf_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::   R_m_r_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_gr_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) :: R_m_cf2_d(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      
      REAL, INTENT(INOUT) :: R_v_d(vdims_s%i_start:vdims_s%i_end,     &
                                   vdims_s%j_start:vdims_s%j_end,     &
                                   vdims_s%k_start:vdims_s%k_end)

! Pressure perturbation

      REAL, INTENT(INOUT) :: exner_np1(pdims_s%i_start:pdims_s%i_end, &
                                       pdims_s%j_start:pdims_s%j_end, &
                                       pdims_s%k_start:pdims_s%k_end+1)


      REAL, INTENT(INOUT) :: p_star_np1(pdims_s%i_start:pdims_s%i_end,&
                                        pdims_s%j_start:pdims_s%j_end)

! Fields at next time level

      REAL, INTENT(INOUT) :: u_np1(udims_s%i_start:udims_s%i_end,     &
                                   udims_s%j_start:udims_s%j_end,     &
                                   udims_s%k_start:udims_s%k_end)

      REAL, INTENT(INOUT) :: v_np1(vdims_s%i_start:vdims_s%i_end,     &
                                   vdims_s%j_start:vdims_s%j_end,     &
                                   vdims_s%k_start:vdims_s%k_end)

      REAL, INTENT(INOUT) :: rho_np1(pdims_s%i_start:wdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,   &
                                     pdims_s%k_start:pdims_s%k_end)

      REAL, INTENT(INOUT) ::  m_v_np1(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: m_cl_np1(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: m_cf_np1(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) ::  m_r_np1(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: m_gr_np1(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: m_cf2_np1(tdims_s%i_start:wdims_s%i_end, &
                                       tdims_s%j_start:tdims_s%j_end, &
                                       tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: thetav_np1(tdims_s%i_start:wdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT) :: w_np1(wdims_s%i_start:wdims_s%i_end,     &
                                   wdims_s%j_start:wdims_s%j_end,     &
                                   wdims_s%k_start:wdims_s%k_end)

      REAL, INTENT(INOUT) :: etadot_np1(wdims_s%i_start:wdims_s%i_end,&
                                        wdims_s%j_start:wdims_s%j_end,&
                                        wdims_s%k_start:wdims_s%k_end)

! Output array for plotting 

      REAL, INTENT(OUT)  ::                                           &
         exner_prime_term(pdims_s%i_start:wdims_s%i_end,              &
                          pdims_s%j_start:pdims_s%j_end,              &
                          pdims_s%k_start:pdims_s%k_end)

! Local arrays used for departure points

      REAL                :: R_u_a(udims_s%i_start:udims_s%i_end,     &
                                   udims_s%j_start:udims_s%j_end,     &
                                   udims_s%k_start:udims_s%k_end)

      REAL                :: R_v_a(vdims_s%i_start:vdims_s%i_end,     &
                                   vdims_s%j_start:vdims_s%j_end,     &
                                   vdims_s%k_start:vdims_s%k_end)

      REAL                ::   R_w_a(wdims%i_start:wdims%i_end,       &
                                     wdims%j_start:wdims%j_end,       &
                                     wdims%k_start:wdims%k_end)

      REAL                ::  R_theta_a(tdims_s%i_start:wdims_s%i_end,&
                                        tdims_s%j_start:tdims_s%j_end,&
                                        tdims_s%k_start:tdims_s%k_end)

      REAL                :: R_rho_a(pdims_s%i_start:wdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,   &
                                     pdims_s%k_start:pdims_s%k_end)


      REAL                ::  R_p_a(pdims_s%i_start:wdims_s%i_end,    &
                                    pdims_s%j_start:pdims_s%j_end,    &
                                    pdims_s%k_start:pdims_s%k_end)

      REAL                :: R_etadot_a(wdims%i_start:wdims%i_end,    &
                                        wdims%j_start:wdims%j_end,    &
                                        wdims%k_start:wdims%k_end)  
      REAL ::                                                         &
                       rho_divu_out(pdims_s%i_start:wdims_s%i_end,    &
                                    pdims_s%j_start:pdims_s%j_end,    &
                                    pdims_s%k_start:pdims_s%k_end)

      REAL ::                                                         &
                                RHS(pdims_s%i_start:wdims_s%i_end,    &
                                    pdims_s%j_start:pdims_s%j_end,    &
                                    pdims_s%k_start:pdims_s%k_end)

      REAL ::                                                         &
             R_v_North_fixd(vdims_s%i_start:vdims_s%i_end,            &
                            vdims_s%k_start:vdims_s%k_end),           &
             R_v_South_fixd(vdims_s%i_start:vdims_s%i_end,            &
                            vdims_s%k_start:vdims_s%k_end),           &
              R_v_North_var(vdims_s%i_start:vdims_s%i_end,            &
                            vdims_s%k_start:vdims_s%k_end),           &
              R_v_South_var(vdims_s%i_start:vdims_s%i_end,            &
                            vdims_s%k_start:vdims_s%k_end)

      INTEGER                              :: i, j, k, iter
      REAL                                 :: del_rho
      REAL                                 :: ex_err(2)
      REAL                                 :: init_err, fin_err
      INTEGER                              :: no_its
      INTEGER                              :: istat1
      
      REAL :: RHS_err, tol_save

!     Fast Physics source terms                 
      REAL                :: S_u(udims_s%i_start:udims_s%i_end,       &
                                 udims_s%j_start:udims_s%j_end,       &
                                 udims_s%k_start:udims_s%k_end)

      REAL                :: S_v(vdims_s%i_start:vdims_s%i_end,       &
                                 vdims_s%j_start:vdims_s%j_end,       &
                                 vdims_s%k_start:vdims_s%k_end)

      REAL                ::   S_w(wdims%i_start:wdims%i_end,         &
                                   wdims%j_start:wdims%j_end,         &
                                   wdims%k_start:wdims%k_end)

      REAL                :: S_thetav(tdims_s%i_start:wdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL              ::    S_m_v(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_cl(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_cf(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::    S_m_r(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::   S_m_gr(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL              ::  S_m_cf2(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)
          

      CHARACTER(len=256)            :: Cmessage
      CHARACTER(len=15)             :: Routine = 'eg_sl_helmholtz'
      INTEGER                       :: ICODE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_SL_HELMHOLTZ',zhook_in,zhook_handle)

      ICODE = 0

      IF ( pre_type < 0 .OR. pre_type > 4 ) THEN
        ICODE = 5
        CMESSAGE = 'Invalid preconditioner type'
        CALL Ereport(Routine,ICODE,CMESSAGE)
      END IF
      
      tol_save = tol

!     del_rho = 0.0
!     IF ( .NOT. L_SLICE ) del_rho = 1.0
      del_rho = 1.0

! Adjust R_v_d at North and South boundaries for use in a cyclic LAM

      IF ( model_domain == mt_cyclic_lam ) THEN
        R_v_South_fixd           = R_v_d(:,vdims%j_start,:)
        R_v_North_fixd           = R_v_d(:,vdims%j_end,:)
        R_v_d(:,vdims%j_start,:) = 0.0
        R_v_d(:,vdims%j_end,:)   = 0.0
      END IF

! Calculate pressure deviation from reference profile

      IF ( model_domain == mt_cyclic_lam ) THEN
        ICODE    = 6
        CMESSAGE = 'LAM not available.'
        CALL Ereport(Routine,ICODE,CMESSAGE)
      END IF

      IF ( mype == 0 .AND. GCR_Diagnostics >0 ) THEN
        WRITE(6,*)  "********************************************"
        WRITE(6,*)  "*   iLinear solve for Helmholtz problem    *"
        WRITE(6,*)  "*   ====================================   *"
      END IF     

inner:DO iter = 1, InnIts
      
        DO k = 1, model_levels
          exner_prime_term(:,:,k) = 0.0
        END DO  

! Calculate arrival point quantities1

        CALL eg_helm_rhs_inc(R_u_a, R_v_a, R_w_a,                     &
                            R_p_a, R_theta_a, R_rho_a,                &
                            R_etadot_a, rho_divu_out,                 &
                            R_m_v_d, R_m_cl_d, R_m_cf_d,              &
                            R_m_r_d, R_m_gr_d, R_m_cf2_d,             &
                            u_np1, v_np1, w_np1,                      &
                            exner_np1, p_star_np1,                    &
                            rho_np1, thetav_np1, etadot_np1,          &
                            m_v_np1, m_cl_np1, m_cf_np1,              &
                            m_r_np1, m_gr_np1, m_cf2_np1,             &
                            R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,  &
                            Ih, row_length, rows, n_rows,             &
                            model_levels,model_domain, iter,          &
                            S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf, &
                            S_m_cf2,S_m_r,S_m_gr,                     &
                            RHS,psi_w_surf, psi_w_lid)

! Adjust R_v_a at North and South boundaries for use in a cyclic LAM

        IF ( model_domain == mt_cyclic_lam ) THEN

          R_v_South_var        = R_v_South_fixd +                     &
                                 R_v_a(:,vdims%j_start,:)
          R_v_North_var        = R_v_North_fixd +                     &
                                 R_v_a(:,vdims%j_end,:)

          R_v_a(:,vdims%j_start,:) = 0.0
          R_v_a(:,vdims%j_end,:)   = 0.0

        END IF

         IF ( GCR_Diagnostics > 0 ) THEN
          IF ( PrintStatus >= PrStatus_Normal) THEN
            RHS_err = -1.0
            DO k = 1, model_levels
              DO j = pdims%j_start, pdims%j_end
                DO i = pdims%i_start, pdims%i_end
                  ! Compute error
                  RHS_err = MAX(RHS_err,ABS(RHS(i,j,k)))
                END DO
              END DO
            END DO
          END IF
         END IF

! Solve the linear system
        CALL eg_bicgstab(exner_prime_term,RHS,tol,pre_type,l_rel_tol, &
               row_length, rows, n_rows, model_levels,sc_err_min,     &
               init_err, fin_err, no_its, .true., ICODE)

        IF ( ICODE /= 0 ) THEN
          CMESSAGE='Convergence failure in BiCGstab'
          CALL Ereport(Routine,ICODE,CMESSAGE)
        END IF

        IF ( GCR_Diagnostics > 0 ) THEN
          IF ( PrintStatus >= PrStatus_Normal) THEN

            CALL gc_rmax(1,n_proc,istat1,RHS_err)

            IF ( mype == 0 ) THEN 
              WRITE(6,'(A21,E13.6,A11)')  " *   Error in RHS    ",    &
                                                        RHS_err,      &
                                                 "          *"
            END IF
          END IF
         
          ex_err(1) =-minval(exner_prime_term)
          ex_err(2) = maxval(exner_prime_term)

          CALL gc_rmax(2,n_proc,istat1,ex_err)
          
          IF( mype == 0 ) THEN
            WRITE(6,"(A13,I3,A29)") " * Inner iteration ",iter,       &
                                   "                     *"
            WRITE(6,"(A35,I4,A6)")                                    &
                               " * No. Of linear solver iterations ", &
                                        no_its,"     *"
            WRITE(6,"(A17,E12.6,A16)") " * Initial error ", init_err, &
                                        "               *"
            WRITE(6,"(A17,E12.6,A16)") " *   Final error ", fin_err,  &
                                        "               *"

            WRITE(6,'(A21,E13.6,A11)')  " *   Min exner prime ",      &
                                        -ex_err(1), "          *"
            WRITE(6,'(A21,E13.6,A11)')  " *   Max exner prime ",      &
                                         ex_err(2), "          *"
             WRITE(6,'(A21,E13.6,A11)')  " *   Error in RHS    ",RHS_err, &
                                                 "          *"              
             WRITE(6,'(A21,E13.6,A11)')  " *   Solver tol      ",tol,     &
                                                 "          *" 
          END IF
        END IF

! Calculate time level n+1 quantities
        CALL eg_add_pressure_inc(                                     &
           exner_prime_term,exner_np1,p_star_np1,                     &
           R_v_South_var, R_v_North_var,                              &
           u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,      &
           m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
           R_w_d,R_u_a, R_v_a, R_w_a, R_theta_a, R_etadot_a, R_p_a,   &
           Ih, n_rows, row_length, rows, model_levels)

        tol = tol/tol_sc_fact
        
      END DO inner

      IF( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
         WRITE(6,fmt='(A)')                                           &
               "********************************************"
         WRITE(6,fmt='(A)')  " "
      END IF      
      tol = tol_save

      IF (lhook) CALL dr_hook('EG_SL_HELMHOLTZ',zhook_out,zhook_handle)

      END SUBROUTINE eg_SL_Helmholtz_inc
      END MODULE eg_sl_helmholtz_inc_mod
