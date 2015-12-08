! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_thermo_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_thermo(                                              &
                row_length, rows, n_rows, model_levels,               &
                halo_i, halo_j, offx, offy,datastart, g_i_pe,         &
                l_inc_solver,  high_order_scheme,                     &
                monotone_scheme, l_high, l_mono,                      &
                alpha_theta,                                          &
                r_theta,  etadot,r_theta_d, error_code )


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

USE mym_option_mod, ONLY:                                             &
                            bdy_tke, l_adv_turb_field, tke_levels

USE atm_fields_bounds_mod
USE integrity_mod
USE eg_interpolation_eta_mod
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE departure_pts_mod
USE Field_Types
USE eg_swap_bounds_mod
USE dynamics_input_mod, ONLY: l_sl_bc_correction

IMPLICIT NONE
!
! Description:
!  Find timelevel n dependent quantity R_theta_d
!  
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.
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

LOGICAL, INTENT(IN) :: l_inc_solver
INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

! MPP options
INTEGER, INTENT(IN) ::                                                &
  halo_i,                                                             &
                     ! Size of halo in i.
  halo_j,                                                             &
                     ! Size of halo in j.
  offx,                                                               &
                     ! Size of small halo in i
  offy
                     ! Size of small halo in j.

 

INTEGER, INTENT(IN) ::                                                &
  datastart(3),                                                       &
                     ! First gridpoints held by this processor.
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction


! Loop index bounds for arrays defined on p, u, v points respectively


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono
                   ! True, if interpolation required to be monotone.
 
REAL, INTENT(IN) :: alpha_theta  ! time-weight for theta eqn
                                 ! and model timestep


! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  r_theta(1-offx:row_length+offx,1-offy:rows+offy,                    &
          0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_theta_d(1-offx:row_length+offx,1-offy:rows+offy,                  &
            0:model_levels)


! Local variables

INTEGER :: i,j,k, number_of_inputs
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

! tmp arrays

REAL :: rdz


REAL :: work(1-halo_i:row_length+halo_i,                              &
              1-halo_j:rows+halo_j,0:model_levels)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_THERMO',zhook_in,zhook_handle)

IF (integrity_test)                                                   &
  CALL check_hash_m(R_theta,       SIZE(R_theta),       'R_t__')


DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
       work(i,j,k) = r_theta(i,j,k)
    END DO
  END DO
END DO

CALL eg_swap_bounds( work,tdims_l,fld_type_p, .FALSE. )

number_of_inputs = 1

! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

CALL eg_interpolation_eta(                                            &
                     eta_theta_levels,fld_type_w,                     &
                     number_of_inputs,                                &
                     row_length, rows, model_levels+1,                &
                     rows,                                            &
                     row_length, rows, model_levels+1,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_w, depart_xi1_w,                      &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     offx, offy, offx ,offy, error_code,              &
                     work, r_theta_d, k_int_linear_in=k_int_linear)


! Compute thetav_ref_eta

IF( .NOT. l_inc_solver ) THEN
! NOTE: etadot = 0 at k=0 and k= model_levels.
!       R_theta_d does not need updating there.

  DO k=1, model_levels-1
    rdz = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k) )
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_theta_d(i,j,k) = r_theta_d(i,j,k)                           &
                          -thetav_ref_pro(i,j,k)                      &
                          +alpha_theta*timestep*etadot(i,j,k)         &
                          *( thetav_ref_eta(i,j,k+1)                  &
                            -thetav_ref_eta(i,j,k)                    &
                           )*rdz
      END DO
    END DO
  END DO

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      k = 0
      r_theta_d(i,j,k) = r_theta_d(i,j,k) - thetav_ref_pro(i,j,k)
      k = model_levels
      r_theta_d(i,j,k) = r_theta_d(i,j,k) - thetav_ref_pro(i,j,k)
    END DO
  END DO
END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(R_theta_d,       SIZE(R_theta_d),       'R_t_d')

IF (lhook) CALL dr_hook('EG_SL_THERMO',zhook_out,zhook_handle)

END SUBROUTINE eg_sl_thermo
END MODULE eg_sl_thermo_mod
