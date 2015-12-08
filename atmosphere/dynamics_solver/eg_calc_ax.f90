! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_calc_ax_mod
      IMPLICIT NONE
      CONTAINS

      SUBROUTINE eg_Calc_Ax(Ax,x)


      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE Field_Types
      USE helmholtz_const_matrix_mod
      USE eg_swap_bounds_mod
      IMPLICIT NONE

!
! Description: This subroutine multiplies the input field
!              by the Constant part of the Helmholtz operator
!
! Method: ENDGame formulation version 3.02
!
! Method:
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! x is the input vector and Ax is the result of multiplying
! x by the matrix A

      REAL,     INTENT(INOUT) ::                                               &
           Ax(pdims_s%i_start:pdims_s%i_end,                                   &
              pdims_s%j_start:pdims_s%j_end,                                   &
              pdims_s%k_start:pdims_s%k_end),                                  &
            x(pdims_s%i_start:pdims_s%i_end,                                   &
              pdims_s%j_start:pdims_s%j_end,                                   &
              pdims_s%k_start:pdims_s%k_end)


      INTEGER                       :: i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      
      IF (lhook) CALL dr_hook('EG_CALC_AX',zhook_in,zhook_handle)

      CALL eg_swap_bounds(x,pdims_s,fld_type_p,.false.)

!$OMP PARALLEL DO PRIVATE(i,j)
      DO k = pdims_s%k_start,pdims_s%k_end
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            Ax(i,j,k) =  x(i,j,k) +                                            &
                            Hlm_Lw(i,j,k)*x(i-1,j,k)+Hlm_Le(i,j,k)*x(i+1,j,k)  &
                           +Hlm_Ls(i,j,k)*x(i,j-1,k)+Hlm_Ln(i,j,k)*x(i,j+1,k)
          END DO
        END DO


        IF( k == 1 ) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              Ax(i,j,k) = Ax(i,j,k) + Hlm_Lu(i,j,k)*x(i,j,k+1)
            END DO
          END DO
            
        ELSE IF( k == pdims_s%k_end) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              Ax(i,j,k) = Ax(i,j,k) + Hlm_Ld(i,j,k)*x(i,j,k-1)
            END DO
          END DO
        ELSE
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              Ax(i,j,k) = Ax(i,j,k) +                                         &
                              Hlm_Ld(i,j,k)*x(i,j,k-1)+Hlm_Lu(i,j,k)*x(i,j,k+1)
            END DO
          END DO
        END IF
      END DO
!$OMP END PARALLEL DO

      IF (lhook) CALL dr_hook('EG_CALC_AX',zhook_out,zhook_handle)

      END SUBROUTINE eg_calc_ax
      END MODULE eg_calc_ax_mod
