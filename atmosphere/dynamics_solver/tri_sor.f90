      MODULE tri_sor_mod
      IMPLICIT NONE
      CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE Tri_sor(Px,x,NoIts)


      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE Field_Types
      USE helmholtz_const_matrix_mod
      USE UM_ParVars, ONLY : datastart
      USE eg_swap_bounds_mod

      IMPLICIT NONE

!
! Description: SOR type preconditioner using the 
!              tridiagonal part of the matrix
!              + red-black ordering to improve parallelization
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

      INTEGER,  INTENT(IN)    :: NoIts

! Input vector x and output vector Px (i.e. x multiplied by the
! Pre/Post-conditioner).

      REAL,     INTENT(INOUT) ::                                      &
                                 x(pdims_s%i_start:pdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)

      REAL,     INTENT(OUT)   ::                                      &
                                Px(pdims_s%i_start:pdims_s%i_end,     &
                                   pdims_s%j_start:pdims_s%j_end,     &
                                   pdims_s%k_start:pdims_s%k_end)

!     Local variables

      REAL    :: rlx
      INTEGER :: io, jo, ioff, i, j, k, it, isweep
      INTEGER :: i_start, i_end, j_start, j_end

      REAL    :: tmp(pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%k_start:pdims_s%k_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('TRI_SOR',zhook_in,zhook_handle)

      io = MOD(datastart(1) + datastart(2), 2)
      rlx = 1.5

      DO k = pdims_s%k_start,pdims_s%k_end
        Px(:,:,k) = 0.0
      END DO

      DO it = 1, NoIts

! Build RHS using the red/black ordering

        DO isweep = 0, 1

          jo      = io + isweep + 1

! This needs fixing so that the swap_bounds call can be moved

          ioff    = 0 !1 - isweep
          i_start = pdims%i_start - ioff
          j_start = pdims%j_start - ioff
          i_end   = pdims%i_end   + ioff
          j_end   = pdims%j_end   + ioff

! Now solve tridiagonal matrix

! Forward sweep

!$OMP PARALLEL DO PRIVATE(i,k,ioff, tmp)
          DO j = j_start, j_end

            ioff = MOD(jo+j,2)

            IF( it == 1 .AND. isweep == 0 ) THEN
              k = 1
              DO i = i_start+ioff, i_end, 2
                tmp(i,k) = x(i,j,k)
              END DO

              DO k = 2, pdims_s%k_end-1
                DO i = i_start+ioff, i_end,2
                  tmp(i,k) = x(i,j,k)

                  tmp(i,k) = (tmp(i,k) - Hlm_Ld(i,j,k)*tmp(i,k-1))    &
                             *Hd_k(i,j,k)
                END DO
              END DO

              k = pdims_s%k_end
              DO i = i_start+ioff, i_end,2

                tmp(i,k) = x(i,j,k)
                tmp(i,k) = (tmp(i,k) - Hlm_Ld(i,j,k)*tmp(i,k-1))      &
                             *Hd_k(i,j,k)
              END DO
            ELSE
              k = 1
              DO i = i_start+ioff, i_end, 2
                tmp(i,k) = x(i,j,k) - Px(i,j,k)                       &
                            - Hlm_Lw(i,j,k)*Px(i-1,j,k)               &
                            - Hlm_Le(i,j,k)*Px(i+1,j,k)               &
                            - Hlm_Ls(i,j,k)*Px(i,j-1,k)               &
                            - Hlm_Ln(i,j,k)*Px(i,j+1,k)               &
                              - Hlm_Lu(i,j,k)*Px(i,j,k+1)
              END DO

              DO k = 2,pdims_s%k_end-1
                DO i = i_start+ioff, i_end,2
                  tmp(i,k) = x(i,j,k) -  Px(i,j,k)                    &
                              - Hlm_Lw(i,j,k)*Px(i-1,j,k)             &
                              - Hlm_Le(i,j,k)*Px(i+1,j,k)             &
                              - Hlm_Ls(i,j,k)*Px(i,j-1,k)             &
                              - Hlm_Ln(i,j,k)*Px(i,j+1,k)             &
                              - Hlm_Ld(i,j,k)*Px(i,j,k-1)             &
                              - Hlm_Lu(i,j,k)*Px(i,j,k+1)

                  tmp(i,k) = (tmp(i,k) - Hlm_Ld(i,j,k)*tmp(i,k-1))    &
                             *Hd_k(i,j,k)
                END DO
              END DO

              k = pdims_s%k_end
              DO i = i_start+ioff, i_end,2
                tmp(i,k) = x(i,j,k) - Px(i,j,k)                       &
                         - Hlm_Lw(i,j,k)*Px(i-1,j,k)                  &
                         - Hlm_Le(i,j,k)*Px(i+1,j,k)                  &
                         - Hlm_Ls(i,j,k)*Px(i,j-1,k)                  &
                         - Hlm_Ln(i,j,k)*Px(i,j+1,k)                  &
                         - Hlm_Ld(i,j,k)*Px(i,j,k-1)

                tmp(i,k) = (tmp(i,k) - Hlm_Ld(i,j,k)*tmp(i,k-1))      &
                             *Hd_k(i,j,k)
              END DO
            ENDIF

! Backward sweep

            k = pdims_s%k_end
            DO i = i_start+ioff, i_end, 2
              Px(i,j,k) = Px(i,j,k) + rlx*tmp(i,k)
            END DO

            DO k = pdims_s%k_end-1, 1, -1
              DO i = i_start+ioff, i_end, 2
                tmp(i,k) = tmp(i,k) - Hu_k(i,j,k)*tmp(i,k+1)
                Px(i,j,k)  = Px(i,j,k)  + rlx*tmp(i,k)
              END DO
            END DO

          END DO
!$OMP END PARALLEL DO
          IF( it < NoIts .OR. isweep == 0 ) THEN

              CALL eg_swap_bounds(Px,pdims_s,fld_type_p,.FALSE.)

          ENDIF
        END DO

      END DO

    IF (lhook) CALL dr_hook('TRI_SOR',zhook_out,zhook_handle)

    END SUBROUTINE tri_sor

END MODULE
