! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Subroutine Interface:
MODULE eg_set_adv_winds_mod
IMPLICIT NONE

CONTAINS
SUBROUTINE eg_set_adv_winds(u,v,w,u_adv,v_adv,w_adv,                  &
                           row_length,rows,n_rows,model_levels,       &
                           halo_i, halo_j, l_shallow)

USE proc_info_mod,       ONLY : model_domain
USE integrity_mod
USE metric_terms_mod,    ONLY : h1_xi1_u, h2_xi2_v
USE atm_fields_bounds_mod
USE eg_swap_bounds_mod
USE earth_constants_mod, ONLY : earth_radius
USE level_heights_mod,   ONLY : r_at_u,r_at_v  
USE domain_params
USE field_types
USE parkind1,            ONLY : jpim, jprb       !DrHook
USE yomhook,             ONLY : lhook, dr_hook   !DrHook

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_shallow

INTEGER, INTENT(IN) :: row_length,rows,n_rows,model_levels,           &
                       halo_i, halo_j

REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,                  &
                      udims_s%j_start:udims_s%j_end,                  &
                      udims_s%k_start:udims_s%k_end),                 &
                    v(vdims_s%i_start:vdims_s%i_end,                  &
                      vdims_s%j_start:vdims_s%j_end,                  &
                      vdims_s%k_start:vdims_s%k_end),                 &
                    w(wdims_s%i_start:wdims_s%i_end,                  &
                      wdims_s%j_start:wdims_s%j_end,                  &
                      wdims_s%k_start:wdims_s%k_end)
                      
REAL, INTENT(INOUT) :: u_adv(udims_l%i_start:udims_l%i_end,           &
                             udims_l%j_start:udims_l%j_end,           &
                             udims_l%k_start:udims_l%k_end),          &
                       v_adv(vdims_l%i_start:vdims_l%i_end,           &
                             vdims_l%j_start:vdims_l%j_end,           &
                             vdims_l%k_start:vdims_l%k_end),          &
                       w_adv(wdims_l%i_start:wdims_l%i_end,           &
                             wdims_l%j_start:wdims_l%j_end,           &
                             wdims_l%k_start:wdims_l%k_end)

INTEGER :: i,j,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('EG_SET_ADV_WINDS',zhook_in,zhook_handle)

IF( model_domain /= mt_global ) THEN
      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               u_adv(i,j,k) = u(i,j,k)/h1_xi1_u(i,j,k)
            END DO
         END DO

         DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
               v_adv(i,j,k) = v(i,j,k)/h2_xi2_v(i,j,k)
            END DO
         END DO
      END DO
ELSE
   IF( l_shallow ) THEN

      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               u_adv(i,j,k) = u(i,j,k)/earth_radius
            END DO
         END DO
      END DO

      DO k = 1, model_levels
         DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
               v_adv(i,j,k) = v(i,j,k)/earth_radius
            END DO
         END DO
      END DO

   ELSE

      DO k = 1, model_levels
         DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
               u_adv(i,j,k) = u(i,j,k)/r_at_u(i,j,k)
            END DO
         END DO
      END DO

      DO k = 1, model_levels
         DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
               v_adv(i,j,k) = v(i,j,k)/r_at_v(i,j,k)
            END DO
         END DO
      END DO

   END IF
END IF

DO k = 0, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      w_adv(i,j,k) = w(i,j,k)
    END DO
  END DO
END DO

CALL eg_swap_bounds(u_adv,udims_l,fld_type_u,.TRUE.)
CALL eg_swap_bounds(v_adv,vdims_l,fld_type_v,.TRUE.)
CALL eg_swap_bounds(w_adv,wdims_l,fld_type_p,.FALSE.)

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  u_adv,           SIZE(u_adv),           'u_adv',    &
                  v_adv,           SIZE(v_adv),           'v_adv',    &
                  w_adv,           SIZE(w_adv),           'w_adv')

IF (lhook) CALL dr_hook('EG_SET_ADV_WINDS',zhook_out,zhook_handle)

END SUBROUTINE eg_set_adv_winds
END MODULE eg_set_adv_winds_mod

