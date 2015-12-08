! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE eg_inner_prod_mod
      CONTAINS
      FUNCTION eg_inner_prod(x,y)

      USE global_2d_sums_mod,   ONLY : global_2d_sums
      USE yomhook,              ONLY : lhook, dr_hook
      USE parkind1,             ONLY : jprb, jpim
      USE atm_fields_bounds_mod
      USE proc_info_mod,        ONLY : n_proc,gc_proc_col_group, gc_proc_row_group

      IMPLICIT NONE

!
! Description: Funtion to calculate the inner product of two fields
!              for use by the Linear Solvers
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions


! Input vectors

      REAL, INTENT(IN)   ::                                           &
           x(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end),                          &
           y(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end)

      REAL                          :: eg_Inner_Prod

      INTEGER                       :: i, j, k
      REAL                          :: temp(1)
      REAL                          :: k_sum(pdims%i_end,pdims%j_end)
      INTEGER                       :: rows, row_length, levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER j1,j2,jj,len1,len2

      IF (lhook) CALL dr_hook('EG_INNER_PROD',zhook_in,zhook_handle)

      row_length = pdims%i_end
      rows       = pdims%j_end
      levels     = pdims%k_end

      k_sum = 0.0

      len1 = 4
      len2 = ( rows - 1 )/len1 + 1

! The following code is the simple algorithm that is OpenMPed below and
! is given for reference.
!      DO k = 1, levels
!        DO j = 1, rows
!           DO i = 1, row_length
!               k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)
!            END DO
!         END DO
!      END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(j1,j2,i,j,k)
      DO jj = 1, len1
         j1 = 1 + (jj-1)*len2
         j2 = MIN( j1+len2-1, rows)
         DO k = 1, levels-1, 2
            DO j = j1, j2
               DO i = 1, row_length
                  k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)        &
                                          + x(i,j,k+1)*y(i,j,k+1)
               ENDDO
            ENDDO
         ENDDO
         IF( MOD(levels,2) == 1 ) THEN
            k = levels
            DO j = j1, j2
               DO i = 1, row_length
                  k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

      CALL global_2d_sums(k_sum,pdims%i_end,pdims%j_end , 0, 0, 1, temp)

      eg_Inner_Prod = temp(1)

      IF (lhook) CALL dr_hook('EG_INNER_PROD',zhook_out,zhook_handle)

      END FUNCTION eg_Inner_Prod
      END MODULE
