! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE update_moisture_mod
      IMPLICIT NONE 
       
! Description:         updates moisture to new time level
!        
! Method: 
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Atmosphere Dynamics Solver
!        
! Code description:
! Language: Fortran 95.        
! This code is written to UMDP3 standards.        
      CONTAINS
      SUBROUTINE update_moisture(m_v_np1, m_cl_np1, m_cf_np1, m_r_np1,&
                              m_gr_np1, m_cf2_np1,                    &
                              R_m_v_d, R_m_cl_d, R_m_cf_d,            &
                              R_m_r_d, R_m_gr_d, R_m_cf2_d,           &
                              R_m_v_a, R_m_cl_a, R_m_cf_a,            &
                              R_m_r_a, R_m_gr_a, R_m_cf2_a)



      USE atmos_constants_mod, ONLY : p_zero, kappa, R
      USE yomhook,             ONLY : lhook, dr_hook
      USE parkind1,            ONLY : jprb, jpim
      USE eg_helmholtz_mod
      USE proc_info_mod,       ONLY : global_row_length,              &
                                      at_extremity,                   &
                                      gc_proc_row_group,              &
                                      model_domain,                   &
                                      me
      USE horiz_grid_mod
      USE ref_pro_mod
      USE atm_fields_bounds_mod
      USE UM_ParParams
      USE metric_terms_mod
      USE helmholtz_const_matrix_mod
      USE coriolis_mod
      USE gravity_mod
      USE eg_parameters_mod,   ONLY : pole_consts

      IMPLICIT NONE


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


      REAL, INTENT(IN) ::   R_m_v_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cl_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_cf_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::   R_m_r_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) ::  R_m_gr_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(IN) :: R_m_cf2_a(tdims_s%i_start:tdims_s%i_end,    &
                                    tdims_s%j_start:tdims_s%j_end,    &
                                    tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::    m_v_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::   m_cl_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::   m_cf_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::    m_r_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::   m_gr_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(OUT) ::  m_cf2_np1(tdims_s%i_start:tdims_s%i_end,  &
                                      tdims_s%j_start:tdims_s%j_end,  &
                                      tdims_s%k_start:tdims_s%k_end)

      INTEGER i, j, k
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UPDATE_MOISTURE',zhook_in,zhook_handle)

!$OMP PARALLEL DO PRIVATE(i,j)
      DO k = tdims_s%k_start,tdims_s%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            m_v_np1(i,j,k)   = R_m_v_d(i,j,k)   + R_m_v_a(i,j,k)
            m_cl_np1(i,j,k)  = R_m_cl_d(i,j,k)  + R_m_cl_a(i,j,k)
            m_cf_np1(i,j,k)  = R_m_cf_d(i,j,k)  + R_m_cf_a(i,j,k)
            m_r_np1(i,j,k)   = R_m_r_d(i,j,k)   + R_m_r_a(i,j,k)
            m_gr_np1(i,j,k)  = R_m_gr_d(i,j,k)  + R_m_gr_a(i,j,k)
            m_cf2_np1(i,j,k) = R_m_cf2_d(i,j,k) + R_m_cf2_a(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO


! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_v_np1,tdims%i_end-tdims%i_start+1,           &
                               tdims%j_end-tdims%j_start+1,           &
                               tdims%k_end-tdims%k_start+1,           &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_cl_np1,tdims%i_end-tdims%i_start+1,          &
                                tdims%j_end-tdims%j_start+1,          &
                                tdims%k_end-tdims%k_start+1,          &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_cf_np1,tdims%i_end-tdims%i_start+1,          &
                                tdims%j_end-tdims%j_start+1,          &
                                tdims%k_end-tdims%k_start+1,          &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_r_np1,tdims%i_end-tdims%i_start+1,           &
                               tdims%j_end-tdims%j_start+1,           &
                               tdims%k_end-tdims%k_start+1,           &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_gr_np1,tdims%i_end-tdims%i_start+1,          &
                                tdims%j_end-tdims%j_start+1,          &
                                tdims%k_end-tdims%k_start+1,          &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

! DEPENDS ON: swap_bounds
      CALL swap_bounds(m_cf2_np1,tdims%i_end-tdims%i_start+1,         &
                                 tdims%j_end-tdims%j_start+1,         &
                                 tdims%k_end-tdims%k_start+1,         &
                    tdims_s%halo_i, tdims_s%halo_j,fld_type_p,.FALSE.)

   IF (lhook) CALL dr_hook('UPDATE_MOISTURE',zhook_out,zhook_handle)

   END SUBROUTINE update_moisture

   END MODULE update_moisture_mod
