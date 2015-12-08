! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates rho in the SCM.
!
! SUBROUTINE Calc_Rho

SUBROUTINE calc_rho                                                           &
  ( model_levels, rows, row_length, exner_rho_levels, exner_theta_levels      &
  , theta, p                                                                  &
! ARGUMENTS Out
  , rho )

  USE atmos_constants_mod, ONLY: r
  USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

  IMPLICIT NONE

!
! Description:  To calculate rho for the Single Column Model
!
!
! Method:
!   It follows the code used in the main UM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!
! Code description:
!    Programming standards :
!      Fortran 90, Written to UM coding standards
!      as specified in UMDP 3, vn8.2

! Inputs

  INTEGER ::                                                        &
    model_levels                                                    &
  , rows                                                            &
  , row_length

  REAL ::                                                           &
    exner_rho_levels(row_length, rows, model_levels+1)              &
  , exner_theta_levels(row_length, rows, model_levels)              &
  , p(row_length, rows,model_levels+1)                              &
  , theta(row_length, rows,model_levels)  ! Temperature(K)


! Outputs
  REAL ::                                                           &
    rho(row_length, rows, model_levels)

! Local variables
  INTEGER ::                                                        &
   i,j,k

  REAL ::                                                           &
    temp1

! 3-d work arrays
  REAL ::                                                           &
    weight_upper(row_length, rows, model_levels)                    &
  , weight_lower(row_length, rows, model_levels)

!----------------------------------------------------------------------
!      Calculate r**2 rho from equation of state
!----------------------------------------------------------------------

  k=1
  DO j=1, rows
    DO i=1, row_length
      rho(i,j,k) = r_rho_levels(i,j,k) * r_rho_levels(i,j,k) *      &
                   p(i,j,k) /                                       &
                   (R * theta(i,j,k) * exner_rho_levels(i,j,k))
    END DO
  END DO

! set up vertical interpolation weights between theta levels and rho
! levels
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        weight_upper(i,j,k) = (r_rho_levels(i,j,k)                  &
                               - r_theta_levels(i,j,k-1) )          &
                            / (r_theta_levels(i,j,k)                &
                               - r_theta_levels(i,j,k-1) )
        weight_lower(i,j,k) = 1.0 - weight_upper(i,j,k)
      END DO
    END DO
  END DO

  DO k=2, model_levels
    DO j=1, rows
      DO i=1, row_length
        temp1 = weight_upper(i,j,k) * theta(i,j,k) +                &
                weight_lower(i,j,k) * theta(i,j,k-1)
        rho(i,j,k) = r_rho_levels(i,j,k)*r_rho_levels(i,j,k) *      &
                   p(i,j,k) /(R * temp1 * exner_rho_levels(i,j,k))
      END DO
    END DO
  END DO


  RETURN

END SUBROUTINE calc_rho

