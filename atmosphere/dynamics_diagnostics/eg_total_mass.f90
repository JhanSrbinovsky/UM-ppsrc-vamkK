! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate sum of global mass
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics Advection

   MODULE eg_total_mass_mod
   IMPLICIT NONE

   REAL, SAVE :: total_rho_init

   CONTAINS

      FUNCTION eg_total_gr(rho)

      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
      USE UM_ParVars,            ONLY: gc_proc_row_group,              &
                                       gc_proc_col_group, nproc
      USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
      USE eg_helmholtz_mod,      ONLY: ec_vol
         
      USE global_2d_sums_mod, ONLY : global_2d_sums 
      USE level_heights_mod,  ONLY : xi3_at_rho=>r_rho_levels
      USE gravity_mod,        ONLY : g_rho
      USE earth_constants_mod,   ONLY: g, earth_radius

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Input
      REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end)
! Output
      REAL    :: eg_total_gr
         
! Local

      INTEGER :: i, j, k
      REAL    :: pe_mass(1)
      REAL    :: row_sum(pdims%j_start:pdims%j_end)
      REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end)
                          

! For GCG subroutines (GCOM)
      INTEGER :: lvl    ! Local Vector Length
      INTEGER :: lsl    ! Local Sum Length (the length of the 
                        ! subsection to be summed for each vector)
      INTEGER :: lso    ! Local Sum Offset (element where the 
                        ! summation starts)
      INTEGER :: nv     ! Number of Vectors
      INTEGER :: istat  ! Status of rsum. 0 is OK

! For GC subroutines (GCOM)
      INTEGER :: gc_len ! Number of elements in message

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_in,zhook_handle)

    
      k_sum = 0.0
      DO k = pdims%k_start,pdims%k_end
         DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end
               k_sum(i,j) = k_sum(i,j) + g_rho(i,j,k)*rho(i,j,k)       &
                           *(xi3_at_rho(i,j,k)-earth_radius)           &
                           *ec_vol(i,j,k)
            END DO
         END DO
      END DO


      CALL global_2d_sums(k_sum,pdims%i_end,pdims%j_end,0,0,1,pe_mass)

      eg_total_gr = pe_mass(1)

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_out,zhook_handle)

      END FUNCTION eg_total_gr


      FUNCTION eg_total_mass(rho)

      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
      USE UM_ParVars,            ONLY: gc_proc_row_group,              &
                                       gc_proc_col_group, nproc
      USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
      USE eg_helmholtz_mod,      ONLY: ec_vol
         
      USE global_2d_sums_mod, ONLY: global_2d_sums 
      USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Input
      REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end)
! Output
      REAL    :: eg_total_mass
         
! Local

      INTEGER :: i, j, k
      REAL    :: pe_mass(1)
      REAL    :: row_sum(pdims%j_start:pdims%j_end)
      REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end)
                          

! For GCG subroutines (GCOM)
      INTEGER :: lvl    ! Local Vector Length
      INTEGER :: lsl    ! Local Sum Length (the length of the 
                        ! subsection to be summed for each vector)
      INTEGER :: lso    ! Local Sum Offset (element where the 
                        ! summation starts)
      INTEGER :: nv     ! Number of Vectors
      INTEGER :: istat  ! Status of rsum. 0 is OK

! For GC subroutines (GCOM)
      INTEGER :: gc_len ! Number of elements in message

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_in,zhook_handle)

    
      k_sum = 0.0
      DO k = pdims%k_start,pdims%k_end
         DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end
               k_sum(i,j) = k_sum(i,j) + rho(i,j,k)*ec_vol(i,j,k)
            END DO
         END DO
      END DO

      CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_mass)

      eg_total_mass = pe_mass(1)

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_out,zhook_handle)


      END FUNCTION eg_total_mass


      FUNCTION eg_total_gr_r(rho)

      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
      USE UM_ParVars,            ONLY: gc_proc_row_group,              &
                                       gc_proc_col_group, nproc
      USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
      USE eg_helmholtz_mod,      ONLY: ec_vol
       
      USE global_2d_sums_mod, ONLY : global_2d_sums 
      USE level_heights_mod,  ONLY : xi3_at_rho=>r_rho_levels
      USE gravity_mod,        ONLY : g_rho
      USE earth_constants_mod,   ONLY: g, earth_radius

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Input
      REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end)
! Output
      REAL    :: eg_total_gr_r
         
! Local

      INTEGER :: i, j, k
      REAL    :: pe_mass(1)
      REAL    :: row_sum(pdims%j_start:pdims%j_end)
      REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end)
                        
! For GCG subroutines (GCOM)
      INTEGER :: lvl    ! Local Vector Length
      INTEGER :: lsl    ! Local Sum Length (the length of the 
                        ! subsection to be summed for each vector)
      INTEGER :: lso    ! Local Sum Offset (element where the 
                        ! summation starts)
      INTEGER :: nv     ! Number of Vectors
      INTEGER :: istat  ! Status of rsum. 0 is OK

! For GC subroutines (GCOM)
      INTEGER :: gc_len ! Number of elements in message

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_in,zhook_handle)

    
      k_sum = 0.0
      DO k = pdims%k_start,pdims%k_end
         DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end
               k_sum(i,j) = k_sum(i,j) + g_rho(i,j,k)*rho(i,j,k)       &
                           *(xi3_at_rho(i,j,k)-earth_radius)           &
                           *xi3_at_rho(i,j,k)*ec_vol(i,j,k)
            END DO
         END DO
      END DO


      CALL global_2d_sums(k_sum,pdims%i_end,pdims%j_end,0,0,1,pe_mass)

      eg_total_gr_r = pe_mass(1)

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_out,zhook_handle)

      END FUNCTION eg_total_gr_r


      FUNCTION eg_total_mass_r(rho)

      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
      USE UM_ParVars,            ONLY: gc_proc_row_group,              &
                                       gc_proc_col_group, nproc
      USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
      USE eg_helmholtz_mod,      ONLY: ec_vol
        
      USE global_2d_sums_mod, ONLY: global_2d_sums 
      USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Input
      REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end)
! Output
      REAL    :: eg_total_mass_r
         
! Local

      INTEGER :: i, j, k
      REAL    :: pe_mass(1)
      REAL    :: row_sum(pdims%j_start:pdims%j_end)
      REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end)
                          

! For GCG subroutines (GCOM)
      INTEGER :: lvl    ! Local Vector Length
      INTEGER :: lsl    ! Local Sum Length (the length of the 
                        ! subsection to be summed for each vector)
      INTEGER :: lso    ! Local Sum Offset (element where the 
                        ! summation starts)
      INTEGER :: nv     ! Number of Vectors
      INTEGER :: istat  ! Status of rsum. 0 is OK

! For GC subroutines (GCOM)
      INTEGER :: gc_len ! Number of elements in message

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_in,zhook_handle)

    
      k_sum = 0.0
      DO k = pdims%k_start,pdims%k_end
         DO j = pdims%j_start,pdims%j_end
            DO i = pdims%i_start,pdims%i_end
               k_sum(i,j) = k_sum(i,j)                                 &
                           + rho(i,j,k)*xi3_at_rho(i,j,k)*ec_vol(i,j,k)
            END DO
         END DO
      END DO


      CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_mass)

      eg_total_mass_r = pe_mass(1)

      IF (lhook) CALL dr_hook('EG_TOTAL_MASS',zhook_out,zhook_handle)

      END FUNCTION eg_total_mass_r

      END MODULE
