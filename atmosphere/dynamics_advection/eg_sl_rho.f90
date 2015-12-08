! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_rho_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_rho(                                                 &
                row_length, rows, n_rows, model_levels, halo_i,       &
                halo_j, offx, offy, datastart, g_i_pe,                &
                l_inc_solver, depart_scheme, depart_order,            &
                high_order_scheme, monotone_scheme,                   &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_cartesian, l_shallow, l_rk_dps,                     &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono, lam_max_cfl,            &
                alpha_rho,                                            &
                etadot, u_np1, v_np1, w_np1,                          &
                etadot_np1, u, v, w,r_rho, r_rho_d,                   &
                r_p_p, r_p_p_d,                                       &
                error_code, hm_rhox, hm_rhoy )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE proc_info_mod,     ONLY : mype=>me,                               &
                              nproc=>n_proc,                          &
                              nproc_x=>n_procx,                       &
                              nproc_y=>n_procy,                       &
                              global_row_length, global_rows,         &
                              at_extremity,gc_proc_row_group,         &
                              gc_proc_col_group,model_domain

USE timestep_mod,      ONLY : timestep
USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,       &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE integrity_mod
USE eg_interpolation_eta_mod
USE departure_point_eta_mod
USE eg_dep_pnt_cart_eta_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE departure_pts_mod
USE metric_terms_mod
USE ereport_mod, ONLY : ereport
USE Field_Types
USE eg_parameters_mod, ONLY : interp_dpt_pt, interp_log_rho, l_rho_av_zz
USE eg_swap_bounds_mod
USE conversions_mod,   ONLY : pi

IMPLICIT NONE
!
! Description:
!   Find rho-grid departure point and compute timelevel n departure
!   point dependent quantity R_rho.
!  
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.

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


INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,        &
                        first_constant_r_rho_level

! MPP options
INTEGER, INTENT(IN) ::                                                &
  halo_i,                                                             &
                     ! Size of halo in i.
  halo_j,                                                             &
                     ! Size of halo in j.
  offx,                                                               &
                     ! Size of small halo in i
  offy,                                                               &
                     ! Size of small halo in j.
  datastart(3),                                                       &
                     ! First gridpoints held by this processor.
  g_i_pe(1-halo_i:global_row_length+halo_i),                          &
                     ! processor on my processor-row
                     ! holding a given value in i direction
  lam_max_cfl(2)     ! Max CFL for a LAM allowed near the
                     ! boundaries


! Loop index bounds for arrays defined on p, u, v points respectively


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme,                                                    &
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.
  depart_scheme,                                                      &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                       &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                           &
                     ! code choosing high order
                     ! interpolation scheme used in Depart routine
  depart_monotone_scheme,                                             &
                     ! code choosing monotone
                     ! interpolation scheme used in Depart routine
  interp_vertical_search_tol,                                         &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono,                                                             &
                   ! True, if interpolation required to be monotone.
  l_depart_high,                                                      &
                   ! True if high order interpolation scheme to
                   ! be used in Depart scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to
                   ! be used in Depart scheme

LOGICAL, INTENT(IN) :: l_cartesian, l_shallow, l_rk_dps
LOGICAL, INTENT(IN) :: l_inc_solver
REAL,    INTENT(IN) :: alpha_rho

! wind components

REAL, INTENT(IN) ::                                                   &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy, model_levels),      &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy, model_levels),    &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
         model_levels),                                               &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
        model_levels),                                                &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
        0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  r_rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n (rotated) w-departure point quantity

REAL, INTENT(OUT) ::                                                  &
  r_rho_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(IN) ::                                                   &
  r_p_p(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
REAL, INTENT(OUT) ::                                                  &
  r_p_p_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! HM_rho terms
REAL, INTENT(IN) ::                                                   &
  hm_rhox(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),     &
  hm_rhoy(1-offx:row_length+offx,-offy:n_rows-1+offy,                 &
          model_levels)


!  VarRes reciprocals needed by interpolation subroutine
      Real                                                              &
        recip_lambda_p_m(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_p_0(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_p_p(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_p_p2(1-halo_i:row_length+halo_i)                   &
      , recip_lambda_u_m(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_u_0(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_u_p(1-halo_i:row_length+halo_i)                    &
      , recip_lambda_u_p2(1-halo_i:row_length+halo_i)                   &
      , recip_phi_p_m(1-halo_j:rows+halo_j)                             &
      , recip_phi_p_0(1-halo_j:rows+halo_j)                             &
      , recip_phi_p_p(1-halo_j:rows+halo_j)                             &
      , recip_phi_p_p2(1-halo_j:rows+halo_j)                            &
      , recip_phi_v_m(1-halo_j:n_rows+halo_j)                           &
      , recip_phi_v_0(1-halo_j:n_rows+halo_j)                           &
      , recip_phi_v_p(1-halo_j:n_rows+halo_j)                           &
      , recip_phi_v_p2(1-halo_j:n_rows+halo_j)

! Local variables

! row_length, number of rows and indexing offset for
! the dep point type computed

INTEGER :: dep_row_len, dep_rows, off_i, off_j, off_k, offz,          &
           number_of_inputs

INTEGER :: i, j, k, del_rho

! Co-ordinates of departure points for rho.

REAL :: alpha

! Halo-ed copy of R_rho

REAL ::                                                               &
   work(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,              &
        model_levels)

REAL :: d_xi1_term, d_xi2_term, deta_term, tmp
REAL :: rdxi1, rdxi2, rdxi3, a_p

REAL :: rho_av, rho_ref_term, rho_switch

! End of header

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_RHO',zhook_in,zhook_handle)

IF (integrity_test)                                                   &
  CALL check_hash_m(                                                  &
                    u_np1,           SIZE(u_np1),           'u_np1',  &
                    v_np1,           SIZE(v_np1),           'v_np1',  &
                    w_np1,           SIZE(w_np1),           'w_np1',  &
                    u,               SIZE(u),               'u_adv',  &
                    v,               SIZE(v),               'v_adv',  &
                    w,               SIZE(w),               'w_adv',  &
                    intw_w2rho,      SIZE(intw_w2rho),      'iww2r',  &
                    eta_theta_levels,SIZE(eta_theta_levels),'etatl',  &
                    eta_rho_levels,  SIZE(eta_rho_levels),  'etarl',  &
                    R_p_p,           SIZE(R_p_p),           'Rpp__',  &
                    exner_ref_pro,   SIZE(exner_ref_pro),   'piref',  &
                    etadot_np1,      SIZE(etadot_np1),      'ednp1',  &
                    xi3_at_rho,      SIZE(xi3_at_rho),      'xi3_r',  &
                    xi1_p,           SIZE(xi1_p),           'xi1_p',  &
                    xi2_p,           SIZE(xi2_p),           'xi2_p' )


dep_row_len = pdims%i_end - pdims%i_start + 1
dep_rows    = pdims%j_end - pdims%j_start + 1
off_i = 0
off_j = 0
off_k = 0
offz  = 0

number_of_inputs = 1
!      alpha = alpha_rho
alpha = 0.5

rho_switch = 1.0
IF ( l_rho_av_zz ) rho_switch = 0.0

IF(.NOT.interp_dpt_pt) THEN
  IF( model_domain /= mt_global ) THEN
    CALL eg_dep_pnt_cart_eta(                                         &
               row_length, rows, n_rows, model_levels, halo_i,        &
               halo_j, offx, offy, mype, nproc, nproc_x, nproc_y,     &
               global_row_length, global_rows, datastart,             &
               at_extremity, g_i_pe, gc_proc_row_group,               &
               gc_proc_col_group, model_domain,                       &
               fld_type_p, dep_row_len, dep_rows, off_i ,off_j, off_k,&
               offz,depart_scheme, depart_order, l_rk_dps,            &
               depart_high_order_scheme,                              &
               depart_monotone_scheme, first_constant_r_rho_level,    &
               interp_vertical_search_tol, check_bottom_levels,       &
               l_depart_high,                                         &
               l_depart_mono, lam_max_cfl, alpha, timestep,           &
               eta_theta_levels, eta_rho_levels,                      &
               xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,          &
               etadot, u_np1, v_np1,                                  &
               w_np1, etadot_np1, u, v, w,                            &
               depart_xi1_rho, depart_xi2_rho, depart_xi3_rho )
  ELSE
    CALL departure_point_eta(                                         &
               row_length, rows, n_rows, model_levels, g_i_pe,        &
               fld_type_p, dep_row_len, dep_rows,                     &
               off_i ,off_j, off_k,offz,l_rk_dps, depart_order,       &
               depart_high_order_scheme, depart_monotone_scheme,      &
               first_constant_r_rho_level, l_depart_high,             &
               l_depart_mono,l_shallow, lam_max_cfl,                  &
               alpha, eta_theta_levels, eta_rho_levels,               &
               xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,          &
               etadot, u_np1, v_np1, etadot_np1,u,v,w,depart_xi1_rho, &
               depart_xi2_rho, depart_xi3_rho )
  END IF
END IF


IF (integrity_test)                                                   &
  CALL update_hash_m(depart_xi1_rho,  size(depart_xi1_rho),  'dxi1r', &
                     depart_xi2_rho,  size(depart_xi2_rho),  'dxi2r', &
                     depart_xi3_rho,  size(depart_xi3_rho),  'dxi3r')

IF (interp_log_rho) THEN
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work(i,j,k) = LOG(r_rho(i,j,k))
      END DO
    END DO
  END DO
ELSE
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work(i,j,k) = r_rho(i,j,k)
      END DO
    END DO
  END DO
END IF

CALL eg_swap_bounds( work,pdims_l,fld_type_p,.FALSE. )

!IF ( l_slice ) THEN

!  r_rho_d = r_rho

!     Call ereport("eg_sl_rho", 1,                                     &
!                  "slice code removed" )
!ELSE
       CALL eg_interpolation_eta(                                     &
                     eta_rho_levels,fld_type_p,                       &
                     number_of_inputs,                                &
                     row_length, rows, model_levels,                  &
                     rows,                                            &
                     row_length, rows, model_levels,                  &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_rho,depart_xi1_rho,depart_xi2_rho,    &
                     mype, nproc, nproc_x, nproc_y,                   &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
!  L_sl_halo_reprod=.false.
                     offx, offy, offx ,offy,error_code,               &
                     work, r_rho_d)
!END IF ! L_slice

IF (interp_log_rho) THEN
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_rho_d(i,j,k) = EXP(r_rho_d(i,j,k))
      END DO
    END DO
  END DO
END IF


IF( .NOT. l_inc_solver ) THEN
   a_p = timestep*alpha_rho
   del_rho = 1.0
!  IF ( l_slice ) del_rho = 0.0

! Use bottom BC etadot(0)=0 (following eg_SISL_Init)
    k = 1
    rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
    DO j=pdims%j_start, pdims%j_end
      rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
      DO i=pdims%i_start, pdims%i_end

          rdxi1 = 1.0/( xi1_u(i) - xi1_u(i-1) )

          d_xi1_term = ( hm_rhox(i,  j,k)*u_np1(i,  j,k) -           &
                         hm_rhox(i-1,j,k)*u_np1(i-1,j,k) )           &
                       *rdxi1

          d_xi2_term = ( hm_rhoy(i,j,  k)*v_np1(i,j  ,k) -           &
                         hm_rhoy(i,j-1,k)*v_np1(i,j-1,k) )           &
                       *rdxi2

          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                        h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*        &
                        rho_ref_eta(i,j,k)*etadot_np1(i,j,k) )        &
                      *rdxi3

          tmp=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

! Calculate: R_rho_d = R_rho_d - rho_ref - alpha*dt*div(rho_ref*U)

          r_rho_d(i,j,k) = r_rho_d(i,j,k) - rho_ref_pro(i,j,k) +      &
                           a_p*(                                      &
                           d_xi1_term+d_xi2_term+deta_term )/tmp

! Add: -alpha*dt*rho_ref*div(U) if NOT SLICE (del_rho=1.0)

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )           &
                       *rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )           &
                       *rdxi2

          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                        h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*        &
                        etadot_np1(i,j,k) ) * rdxi3

          r_rho_d(i,j,k) = r_rho_d(i,j,k)-del_rho*a_p*                &
                           rho_ref_pro(i,j,k) *                       &
                           (d_xi1_term+d_xi2_term+deta_term)/tmp

      END DO
    END DO

!$OMP PARALLEL DO PRIVATE(i,j, tmp, d_xi1_term,d_xi2_term,deta_term,  &
!$OMP& rho_ref_term,rho_av,k,rdxi3,rdxi2,rdxi1) DEFAULT(NONE)         &
!$OMP& SHARED(model_levels,pdims, eta_theta_levels,xi2_v,xi1_u,       &
!$OMP& hm_rhox,u_np1,hm_rhoy,v_np1,h1_p_eta,h2_p_eta,h3_p_eta,        &
!$OMP& deta_xi3_theta,rho_ref_eta,etadot_np1,h1_p,h2_p,h3_p,          &
!$OMP& r_rho_d,rho_ref_pro,a_p,h2_xi1_u,h3_xi1_u,deta_xi3_u,          &
!$OMP& h1_xi2_v,h3_xi2_v,deta_xi3_v, deta_xi3,                        &
!$OMP& intw_w2rho,intw_rho2w,rho_switch,del_rho) SCHEDULE(STATIC)
    DO k=2, model_levels-1
      rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
      DO j=pdims%j_start, pdims%j_end
        rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
        DO i=pdims%i_start, pdims%i_end
          rdxi1      = 1.0/( xi1_u(i) - xi1_u(i-1) )

          d_xi1_term = ( hm_rhox(i,  j,k)*u_np1(i,  j,k) -            &
                         hm_rhox(i-1,j,k)*u_np1(i-1,j,k) )            &
                       *rdxi1

          d_xi2_term = ( hm_rhoy(i,j,  k)*v_np1(i,j  ,k) -            &
                         hm_rhoy(i,j-1,k)*v_np1(i,j-1,k) )            &
                       *rdxi2

          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                       h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*         &
                       rho_ref_eta(i,j,k)*etadot_np1(i,j,k) -         &
                       h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*           &
                       h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                       rho_ref_eta(i,j,k-1)*etadot_np1(i,j,k-1) )     &
                       *rdxi3

          tmp=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

! Calculate: R_rho_d = R_rho_d - rho_ref - alpha*dt*div(rho_ref*U)

          r_rho_d(i,j,k) = r_rho_d(i,j,k) - rho_ref_pro(i,j,k) +      &
                           a_p*(                                      &
                           d_xi1_term+d_xi2_term+deta_term )/tmp

! Add: -alpha*dt*rho_ref*div(U) if NOT SLICE (del_rho=1.0)

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )           &
                       *rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )           &
                       *rdxi2


          deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*              &
                       h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*         &
                       etadot_np1(i,j,k) -                            &
                       h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*           &
                       h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                       etadot_np1(i,j,k-1) ) * rdxi3
                       
          rho_av = intw_w2rho(k,1)*(intw_rho2w(k,1)*rho_ref_pro(i,j,k+1)      &
                                   +intw_rho2w(k,2)*rho_ref_pro(i,j,k)    )   &
                 + intw_w2rho(k,2)*(intw_rho2w(k-1,1)*rho_ref_pro(i,j,k)      &
                                   +intw_rho2w(k-1,2)*rho_ref_pro(i,j,k-1)) 

          rho_ref_term = rho_switch*rho_ref_pro(i,j,k)                &
                       + (1.0 - rho_switch)*rho_av               

          r_rho_d(i,j,k) = r_rho_d(i,j,k)-del_rho*a_p*                &
                           rho_ref_term *                             &
!                            rho_ref_pro(i,j,k) *                       &
                           (d_xi1_term+d_xi2_term+deta_term)/tmp

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    k = model_levels
    rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
    DO j=pdims%j_start, pdims%j_end
      rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
      DO i=pdims%i_start, pdims%i_end
          rdxi1      = 1.0/( xi1_u(i) - xi1_u(i-1) )

          d_xi1_term = ( hm_rhox(i,  j,k)*u_np1(i,  j,k) -            &
                         hm_rhox(i-1,j,k)*u_np1(i-1,j,k) )            &
                       *rdxi1

          d_xi2_term = ( hm_rhoy(i,j,  k)*v_np1(i,j  ,k) -            &
                         hm_rhoy(i,j-1,k)*v_np1(i,j-1,k) )            &
                       *rdxi2

          deta_term = ( -h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*         &
                       h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                       rho_ref_eta(i,j,k-1)*etadot_np1(i,j,k-1) )     &
                       *rdxi3

          tmp=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

! Calculate: R_rho_d = R_rho_d - rho_ref - alpha*dt*div(rho_ref*U)

          r_rho_d(i,j,k) = r_rho_d(i,j,k) - rho_ref_pro(i,j,k) +      &
                           a_p*(                                      &
                           d_xi1_term+d_xi2_term+deta_term )/tmp

! Add: -alpha*dt*rho_ref*div(U) if NOT SLICE (del_rho=1.0)

          d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*             &
                       deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                       h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                       deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )           &
                       *rdxi1

          d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*             &
                       deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                       h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                       deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )           &
                       *rdxi2

          deta_term = ( -h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*         &
                       h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                       etadot_np1(i,j,k-1) ) *rdxi3

          r_rho_d(i,j,k) = r_rho_d(i,j,k)-del_rho*a_p*                &
                           rho_ref_pro(i,j,k) *                       &
                           (d_xi1_term+d_xi2_term+deta_term)/tmp

      END DO
    END DO  
ENDIF

IF (integrity_test)                                                   &
  CALL update_hash_m(R_rho_d,         SIZE(R_rho_d),         'R_r_d', &
                     R_p_p_d,         SIZE(R_p_p_d),         'Rpp_d', &
                     depart_xi1_rho,  SIZE(depart_xi1_rho),  'dxi1r', &
                     depart_xi2_rho,  SIZE(depart_xi2_rho),  'dxi2r', &
                     depart_xi3_rho,  SIZE(depart_xi3_rho),  'dxi3r')

IF (lhook) CALL dr_hook('EG_SL_RHO',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sl_rho
END MODULE eg_sl_rho_mod
