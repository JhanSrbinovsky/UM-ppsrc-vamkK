! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

   MODULE gmres1_coef_mod

   CONTAINS

   FUNCTION gmres1_coef(s,t,err)

   USE global_2d_sums_mod,   ONLY : global_2d_sums
   USE yomhook,              ONLY : lhook, dr_hook
   USE parkind1,             ONLY : jprb, jpim
   USE atm_fields_bounds_mod

   IMPLICIT NONE

!
! Description: Function to calculate the inner product of two fields
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

      REAL,               INTENT(IN)   ::                             &
           s(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end),                          &
           t(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end,                           &
             pdims_s%k_start:pdims_s%k_end)

      REAL                          :: gmres1_coef
      REAL,             INTENT(OUT) :: err
      REAL                          :: temp(3)

      INTEGER                       :: i, j, k
      INTEGER                       :: rows, row_length, levels
      REAL                          :: k_sum(pdims%i_end,pdims%j_end,3)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER j1,j2,jj,len1,len2

      IF (lhook) CALL dr_hook('GMRES1_COEF',zhook_in,zhook_handle)

      row_length = pdims%i_end
      rows       = pdims%j_end
      levels     = pdims%k_end

      k_sum = 0.0

      len1 = 4
      len2 = ( rows - 1 )/len1 + 1

! The following code is the straightforward implementation of the 
! OpenMP code following (for reference)
!      DO k = 1, levels
!        DO j = 1, rows
!            DO i = 1, row_length
!              k_sum(i,j,1) = k_sum(i,j,1) + t(i,j,k)**2
!              k_sum(i,j,2) = k_sum(i,j,2) + t(i,j,k)*s(i,j,k)
!              k_sum(i,j,3) = k_sum(i,j,3) + s(i,j,k)**2
!            END DO
!         END DO
!      END DO 
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(j1,j2,i,j,k)               &
!$OMP&            SHARED(len1, len2, rows, levels, row_length,        &
!$OMP&            k_sum, t, s)
      DO jj = 1, len1

        j1 = 1 + (jj-1)*len2
        j2 = MIN( j1+len2-1, rows)

        DO k = 1, levels-1, 2
          DO j = j1, j2
            DO i = 1, row_length
              k_sum(i,j,1) = k_sum(i,j,1)+t(i,j,k)**2 + t(i,j,k+1)**2
              k_sum(i,j,2) = k_sum(i,j,2)+t(i,j,k)*s(i,j,k)           &
                                         +t(i,j,k+1)*s(i,j,k+1)
              k_sum(i,j,3) = k_sum(i,j,3)+s(i,j,k)**2 + s(i,j,k+1)**2
            END DO
          END DO
        END DO

        IF( MOD(levels,2) == 1 ) THEN
          k = levels
          DO j = j1, j2
            DO i = 1, row_length
              k_sum(i,j,1) = k_sum(i,j,1) + t(i,j,k)**2
              k_sum(i,j,2) = k_sum(i,j,2) + t(i,j,k)*s(i,j,k)
              k_sum(i,j,3) = k_sum(i,j,3) + s(i,j,k)**2
            END DO
          END DO
        END IF
      END DO
!$OMP END PARALLEL DO

      CALL global_2d_sums(k_sum, row_length, rows, 0, 0, 3, temp)

      gmres1_coef = temp(2)/temp(1)

      err         = temp(3) + gmres1_coef**2 *temp(1)                 &
                            - 2.0*gmres1_coef*temp(2)

      IF (lhook) CALL dr_hook('GMRES1_COEF',zhook_out,zhook_handle)

   END FUNCTION gmres1_coef
END MODULE gmres1_coef_mod
