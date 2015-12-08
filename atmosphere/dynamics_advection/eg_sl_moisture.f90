! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_moisture_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_moisture(                                            &
              moisture_array_size,                                    &
              row_length, rows, n_rows, model_levels, halo_i,         &
              halo_j, offx, offy,datastart,g_i_pe,high_order_scheme,  &
              monotone_scheme,  l_high, l_mono,l_pc2, L_mcr_rain,     &
              l_mcr_cf2,l_mcr_graup,                                  &
              r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,          &
              cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
              exner_star,r_m_v_d, r_m_cl_d, r_m_cf_d,r_m_r_d,         &
              r_m_gr_d, r_m_cf2_d,cf_star, cfl_star, cff_star,        &
              error_code )

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

USE cloud_inputs_mod, ONLY: l_fixbug_pc2_mixph
USE atm_fields_bounds_mod
USE integrity_mod
USE eg_interpolation_eta_mod
USE horiz_grid_mod
USE metric_terms_mod
USE departure_pts_mod
USE ereport_mod, ONLY : ereport
USE Field_Types
USE dynamics_input_mod, ONLY: l_sl_bc_correction

IMPLICIT NONE
!
! Description:
!   Find departure point timelevel n dependent quantity R_theta_d.
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



! Model dimensions

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

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
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction


! Loop index bounds for arrays defined on p, u, v points respectively


INTEGER, INTENT(IN) :: moisture_array_size

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

LOGICAL, INTENT(IN) :: l_pc2, l_mcr_rain, l_mcr_cf2, l_mcr_graup


! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  r_m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  r_m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  r_m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  r_m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  r_m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  r_m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL, INTENT(OUT) ::                                                  &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy, model_levels)

REAL, INTENT (INOUT) ::                                               &
        cf_star (1-offx:row_length+offx,                              &
                 1-offy:rows+offy, 0:model_levels),                   &
        cfl_star(1-offx:row_length+offx,                              &
                 1-offy:rows+offy, 0:model_levels),                   &
        cff_star(1-offx:row_length+offx,                              &
                 1-offy:rows+offy, 0:model_levels)

REAL, INTENT (INOUT) :: cf_bulk  (qdims_l%i_start:qdims_l%i_end,      &
                        qdims_l%j_start:qdims_l%j_end,                &
                        qdims_l%k_start:qdims_l%k_end)
REAL, INTENT (INOUT) :: cf_liquid(qdims_l%i_start:qdims_l%i_end,      &
                        qdims_l%j_start:qdims_l%j_end,                &
                        qdims_l%k_start:qdims_l%k_end)
REAL, INTENT (INOUT) :: cf_frozen(qdims_l%i_start:qdims_l%i_end,      &
                        qdims_l%j_start:qdims_l%j_end,                &
                        qdims_l%k_start:qdims_l%k_end)

REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end)


INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_m_v_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  r_m_cl_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  r_m_cf_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  r_m_r_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  r_m_gr_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  r_m_cf2_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Local variables

INTEGER :: i,j,k, number_of_inputs
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

! tmp & dummy arrays


REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,      &
                0:model_levels, moisture_array_size )

REAL :: super_array_out(1-offx:row_length+offx,1-offy:rows+offy,          &
                0:model_levels, moisture_array_size)


INTEGER :: count


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_MOISTURE',zhook_in,zhook_handle)

DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
       super_array(i,j,k,1)  = r_m_v(i,j,k)
       super_array(i,j,k,2)  = r_m_cl(i,j,k)
       super_array(i,j,k,3)  = r_m_cf(i,j,k)
    END DO
  END DO
END DO

count = 3
IF( L_pc2 ) THEN
  IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,count+1) = cf_liquid(i,j,k) +                &
                                       cfl_star (i,j,k) +                &
                                       cf_frozen(i,j,k) +                &
                                       cff_star (i,j,k) -                &
                                       cf_bulk(i,j,k) -                  &
                                       cf_star (i,j,k)
          super_array(i,j,k,count+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,count+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
          super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
  ELSE
! original method, advect the bulk cloud fraction
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,count+1) = cf_bulk(i,j,k)   + cf_star(i,j,k)
          super_array(i,j,k,count+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,count+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
          super_array(i,j,k,count+4) = exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      super_array(i,j,0,count+1) = super_array(i,j,1,count+1)
      super_array(i,j,0,count+2) = super_array(i,j,1,count+2)
      super_array(i,j,0,count+3) = super_array(i,j,1,count+3)
      super_array(i,j,0,count+4) = exner_theta_levels(i,j,0)
    END DO
  END DO
  count = count + 4
END IF

IF( L_mcr_rain ) THEN
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
         super_array(i,j,k,count+1)   = r_m_r(i,j,k)
      END DO
    END DO
  END DO
count = count + 1
END IF

IF( l_mcr_cf2 ) THEN
  DO k=0, model_levels
     DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
           super_array(i,j,k,count+1) = r_m_cf2(i,j,k)
        END DO
     END DO
  END DO
count = count + 1
END IF

IF( l_mcr_graup ) THEN
  DO k=0, model_levels
     DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
           super_array(i,j,k,count+1) = r_m_gr(i,j,k)
        END DO
     END DO
  END DO
count = count + 1
END IF


! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                       &
             super_array, row_length, rows,                             &
             moisture_array_size*(model_levels+1),                      &
             halo_i, halo_j, fld_type_p, .false. )

number_of_inputs = moisture_array_size

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
                     model_domain,l_high, l_mono,                     &
                     depart_xi3_w, depart_xi1_w,                      &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     offx, offy, offx ,offy, error_code,              &
                     super_array, super_array_out,                    &
                     k_int_linear_in=k_int_linear)

! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                       &
             super_array_out, row_length, rows,                         &
             moisture_array_size*(model_levels+1),                      &
             offx, offy, fld_type_p, .false. )


! Unpack super array
DO k=0, model_levels
  DO j=pdims_s%j_start, pdims_s%j_end
    DO i=pdims_s%i_start, pdims_s%i_end
       r_m_v_d(i,j,k)  = super_array_out(i,j,k,1)
       r_m_cl_d(i,j,k) = super_array_out(i,j,k,2)
       r_m_cf_d(i,j,k) = super_array_out(i,j,k,3)
    END DO
  END DO
END DO

count = 3
IF( L_pc2 ) THEN
  IF (l_fixbug_pc2_mixph) THEN
! new method, advect the mixed phase cloud fraction
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end
        cf_star(i,j,0) = super_array_out(i,j,1,count+2)             &
                        +super_array_out(i,j,1,count+3)             &
                        -super_array_out(i,j,1,count+1)
        cfl_star(i,j,0)= super_array_out(i,j,1,count+2)
        cff_star(i,j,0)= super_array_out(i,j,1,count+3)
      END DO
    END DO
    DO k=1, model_levels
      DO j=pdims_s%j_start, pdims_s%j_end
        DO i=pdims_s%i_start, pdims_s%i_end
          cf_star(i,j,k) = super_array_out(i,j,k,count+2)          &
                          +super_array_out(i,j,k,count+3)          &
                          -super_array_out(i,j,k,count+1)
          cfl_star(i,j,k)   = super_array_out(i,j,k,count+2)
          cff_star(i,j,k)   = super_array_out(i,j,k,count+3)
          exner_star(i,j,k) = super_array_out(i,j,k,count+4)
        END DO
      END DO
    END DO
  ELSE
! original method, advect the bulk cloud fraction
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end
        cf_star(i,j,0) = super_array_out(i,j,1,count+1)
        cfl_star(i,j,0)= super_array_out(i,j,1,count+2)
        cff_star(i,j,0)= super_array_out(i,j,1,count+3)
      END DO
    END DO
    DO k=1, model_levels
      DO j=pdims_s%j_start, pdims_s%j_end
        DO i=pdims_s%i_start, pdims_s%i_end
          cf_star(i,j,k)    = super_array_out(i,j,k,count+1)
          cfl_star(i,j,k)   = super_array_out(i,j,k,count+2)
          cff_star(i,j,k)   = super_array_out(i,j,k,count+3)
          exner_star(i,j,k) = super_array_out(i,j,k,count+4)
        END DO
      END DO
    END DO
  END IF !l_fixbug_pc2_mixph
  count = count + 4
END IF

IF( L_mcr_rain ) THEN
DO k=0, model_levels
  DO j=pdims_s%j_start, pdims_s%j_end
    DO i=pdims_s%i_start, pdims_s%i_end
       r_m_r_d(i,j,k) = super_array_out(i,j,k,count+1)
    END DO
  END DO
END DO
count = count + 1
END IF

IF( l_mcr_cf2 ) THEN
  DO k=0, model_levels
     DO j=pdims_s%j_start, pdims_s%j_end
        DO i=pdims_s%i_start, pdims_s%i_end
           r_m_cf2_d(i,j,k) = super_array_out(i,j,k,count+1)
        END DO
     END DO
  END DO
count = count + 1
END IF

IF( l_mcr_graup ) THEN
  DO k=0, model_levels
     DO j=pdims_s%j_start, pdims_s%j_end
        DO i=pdims_s%i_start, pdims_s%i_end
           r_m_gr_d(i,j,k) = super_array_out(i,j,k,count+1)
        END DO
     END DO
  END DO
END IF

IF (integrity_test) THEN
  CALL update_hash_m(R_m_v_d,         SIZE(R_m_v_d),         'R_mvd', &
                     R_m_cl_d,        SIZE(R_m_cl_d),        'Rmcld', &
                     R_m_cf_d,        SIZE(R_m_cf_d),        'Rmcfd')

  IF( L_mcr_rain )                                                    &
    CALL update_hash_m( R_m_r_d,         SIZE(R_m_r_d),      'Rmr_d') 
  IF( l_mcr_graup )                                                   &
    CALL update_hash_m( R_m_gr_d,        SIZE(R_m_gr_d),     'Rmgrd')
  IF( l_mcr_cf2   )                                                   &
    CALL update_hash_m( R_m_cf2_d,       SIZE(R_m_cf2_d),    'Rmc2d')

  IF( L_pc2 ) THEN
    CALL update_hash_m( cf_star,         SIZE(cf_star),      'cfstr', &
                      cfl_star,         SIZE(cfl_star),      'cflst', &
                      cff_star,         SIZE(cff_star),      'cffst', &
                      exner_star,       SIZE(exner_star),    'exstr')
  END IF
END IF

IF (lhook) CALL dr_hook('EG_SL_MOISTURE',zhook_out,zhook_handle)

END SUBROUTINE eg_sl_moisture
END MODULE eg_sl_moisture_mod
