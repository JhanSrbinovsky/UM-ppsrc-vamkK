! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates exner in the SCM.
!
!  SUBROUTINE calc_press

SUBROUTINE calc_press                                                         &
! In data
  ( model_levels, wet_levels, rows, row_length, p, theta, q                   &
  , l_calc_exner, l_calc_rho                                                  &
! InOut
  , rho                                                                       &
! Out data
  , exner_theta_levels, exner_rho_levels, p_theta_levels, rp, rp_theta        &
  , p_star )

USE earth_constants_mod, ONLY : g

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE atmos_constants_mod, ONLY:  kappa, p_zero

  IMPLICIT NONE

!
! Description:  To calculate exner, pressure on theta levels,
!               the reciprocol or pressure and rho (if required)
!               for the Single Column Model
!
!
! Method:
!   It follows the code used in the main UM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code description:
!   FORTRAN 90
!   This code is written to UM programming standards version 8

! Inputs

  INTEGER ::                                                        &
    model_levels                                                    &
  , wet_levels                                                      &
  , rows                                                            &
  , row_length

  REAL ::                                                           &
    p(row_length,rows,model_levels+1)   &
  , q(row_length,rows,wet_levels)       &
  , theta(row_length,rows,model_levels)  ! Temperature(K)

  LOGICAL ::      &
    l_calc_exner  &
                   ! If true then exner is calculated from p,
                   ! If false then p is calculated from exner
  , l_calc_rho     ! If true then rho is calculated

! InOut
  REAL ::                                             &
    rho(row_length,rows,model_levels)

! Outputs
  REAL ::                                             &
    exner_rho_levels(row_length,rows,model_levels+1)  &
  , exner_theta_levels(row_length,rows, model_levels) &
  , p_theta_levels(row_length,rows,model_levels)      &
  , rp_theta(row_length,rows,model_levels)            &! reciprocol pressure
  , rp(row_length,rows,model_levels+1)                &! reciprocol pressure
  , p_star(row_length,rows)                             ! surface pressure

  CHARACTER(LEN=80) ::                                                 &
         Cmessage              ! Error message if ICODE >0
  CHARACTER(LEN=*), PARAMETER :: RoutineName = 'Calc_Press'

  INTEGER ::                                                        &
    icode

! Local variables
  INTEGER ::                                                        &
   i,j,k

  REAL ::                                                           &
    constant

!----------------------------------------------------------------------

! 0.1 Initialise secondary arrays.
! calculate p from exner_rho_levels if l_calc_exner is false
! calculate exner_rho_levels from p if l_calc_exner is true
! [halos required for diagnostic calculations at T+0 in INITDIAG.]
  constant = 1./ kappa
  IF (l_calc_exner) THEN
    DO k=1, model_levels+1
      DO j=1, rows
        DO i=1, row_length
          exner_rho_levels(i,j,k) =                                 &
                   (p(i,j,k)/p_zero)**(1.0/constant)
        END DO
      END DO
    END DO
  ELSE
    DO k=1, model_levels+1
      DO j=1, rows
        DO i=1, row_length
          p(i,j,k) = (exner_rho_levels(i,j,k) ** constant) *        &
                              p_zero
        END DO
      END DO
    END DO
  END IF

! calculate Exner at theta levels
! DEPENDS ON: calc_exner_at_theta
  CALL calc_exner_at_theta                                                    &
    ( r_theta_levels, r_rho_levels, exner_rho_levels, row_length, rows        &
    , model_levels, 0, 0, 0, 0, exner_theta_levels, .FALSE. )

! calculate p at theta_levels
! DEPENDS ON: calc_p_from_exner
  CALL calc_p_from_exner                                                      &
    ( p_theta_levels, row_length, rows, model_levels, 0, 0                    &
    , exner_theta_levels, .FALSE. )

! Calculate rho if required (ie if this is the start of the run).

  IF (l_calc_rho) THEN
! DEPENDS ON: calc_rho
    CALL calc_rho                                                             &
      ( model_levels, rows, row_length, exner_rho_levels, exner_theta_levels  &
      , theta, p, rho )
  END IF

! calculate p_star using rho and p on model levels
! DEPENDS ON: calc_p_star
  CALL calc_p_star                                                            &
    ( r_theta_levels, r_rho_levels, p, rho, row_length, rows, model_levels    &
    , 0, 0, 0, 0, p_star )

! calculate reciprocol pressures

  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        rp(i,j,k) = 1.0/p(i,j,k)
        rp_theta(i,j,k) = 1.0/p_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  DO j=1, rows
    DO i=1, row_length
      rp(i,j,model_levels+1) = 1.0/p(i,j,model_levels+1)
    END DO
  END DO

  RETURN

END SUBROUTINE calc_press

