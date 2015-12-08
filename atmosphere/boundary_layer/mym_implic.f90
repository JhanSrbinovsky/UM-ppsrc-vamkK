! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_implic----------------------------------------------
!
!  Purpose: To solve the tri-diagonal equations for the prognostic
!           variables in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_implic(levels, kst, ken, aa, bb, cc, qq)

  USE atm_fields_bounds_mod, ONLY: pdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::                                                &
     levels,                                                            &
                    ! number of levels of variables to be solved
     kst,                                                               &
                    ! index of start level to be solved
     ken
                    ! index of emd level to be solved

  REAL, INTENT(INOUT) ::                                                &
     aa(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),   &
                    ! coefficients of fields on level K-1
                    ! in the tri-diagonal equation
     bb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),   &
                    ! coefficients on fields level K
                    ! in the tri-diagonal equation
     cc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels),   &
                    ! coefficients on fields level K+1
                    ! in the tri-diagonal equation
     qq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, levels)
                    ! right hand side of the tri-diagonal equation

! Local variables
  INTEGER ::                                                            &
     i, j, k
                    ! Loop indexes

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  ! Solve from top to bottom
  IF (lhook) CALL dr_hook('MYM_IMPLIC',zhook_in,zhook_handle)

  DO k = ken, kst + 1, -1
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        aa(i, j, k - 1) = aa(i, j, k - 1) * bb(i, j, k)
        bb(i, j, k - 1) = bb(i, j, k - 1) * bb(i, j, k)                 &
             - aa(i, j, k) * cc(i, j, k - 1)
        qq(i, j, k - 1) = qq(i, j, k - 1) * bb(i, j, k)                 &
             - qq(i, j, k) * cc(i, j, k - 1)
      END DO
    END DO
  END DO

  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      qq(i, j, kst) = qq(i, j, kst) / bb(i, j, kst)
    END DO
  END DO

  ! Solve from bottom to top
  DO k = kst + 1, ken
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        qq(i, j, k) = (qq(i, j, k) - aa(i, j, k) *                      &
             qq(i, j, k - 1)) / bb(i, j, k)
      END DO
    END DO
  END DO
  IF (lhook) CALL dr_hook('MYM_IMPLIC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_implic

