! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates the heights in the SCM.
!
! Subroutine Interface:

SUBROUTINE calc_levels                                                        &
! Input data
  ( orog, height_gen_method                                                   &
  , bl_levels, model_levels, rows, row_length )

  USE vertnamelist_mod, ONLY:                                                 &
      z_top_of_model,  first_constant_r_rho_level 

  USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels,                  &
                             eta_theta_levels, eta_rho_levels


  USE earth_constants_mod, ONLY: earth_radius

  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

!
! Description:  To set up r_theta_levels and r_rho_levels for the
!               Single Column Model
!
!
! Method:
!   It follows the code used in the main UM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!

! Code description:
!      Fortran 90, Written to UM coding standards
!      as specified in UMDP 3, vn8.2

!     INCLUDED COMDECKS

! Inputs
  INTEGER ::                                                        &
    bl_levels                                                       &
  , model_levels                                                    &
  , rows                                                            &
  , row_length                                                      &
  , height_gen_method

  REAL ::                                                           &
    orog(row_length, rows)

  CHARACTER(LEN=80) ::                                                 &
         Cmessage              ! Error message if ICODE >0

  CHARACTER(LEN=*), PARAMETER :: RoutineName = 'Calc_Levels'

  INTEGER :: ErrorStatus

! Local variables
  REAL ::                                                           &
    r_ref_theta(model_levels)                                       &
  , r_ref_rho(model_levels)

  INTEGER ::                                                        &
    i,j,k

  INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
  INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

!----------------------------------------------------------------------
!     Set up heights
!----------------------------------------------------------------------

! Set reference profile

  DO k=1, model_levels
    r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
    r_ref_rho(k)   = eta_rho_levels(k)   * z_top_of_model
  END DO

! Set bottom level, ie orography
  DO j=1, rows
    DO i=1, row_length
      r_theta_levels(i,j,0) = orog(i,j) + earth_radius
    END DO
  END DO
! For constant levels set r to be a constant on the level
  DO k=first_constant_r_rho_level, model_levels
    DO j=1, rows
      DO i=1, row_length
        r_theta_levels(i,j,k) = earth_radius + r_ref_theta(k)
        r_rho_levels(i,j,k)   = earth_radius + r_ref_rho(k)
      END DO
    END DO
  END DO

  SELECT CASE( height_gen_method)
    CASE( height_gen_original)
! The original version of height generation used in the SI dynamics
!
! For boundary layer levels set depth to be constant.
      DO k=1, bl_levels
        DO j=1, rows
          DO i=1, row_length
            r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                       r_ref_theta(k)
            r_rho_levels(i,j,k)   = r_theta_levels(i,j,0) +         &
                                       r_ref_rho(k)
          END DO
        END DO
      END DO
! For intemediate levels use linear relaxation to constant value.
! set orographic heights.
      DO k=bl_levels+1, first_constant_r_rho_level-1
        DO j=1, rows
          DO i=1, row_length

            r_rho_levels(i,j,k) =                                   &
              ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
                r_theta_levels(i,j,bl_levels) ) *                   &
              ( eta_rho_levels(k) -                                 &
                eta_theta_levels(bl_levels) )/                      &
              (eta_rho_levels(first_constant_r_rho_level) -         &
               eta_theta_levels(bl_levels) )                        &
              +  r_theta_levels(i,j,bl_levels)

            r_theta_levels(i,j,k) =                                 &
              ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
                r_theta_levels(i,j,bl_levels) ) *                   &
              ( eta_theta_levels(k) -                               &
                eta_theta_levels(bl_levels) ) /                     &
              ( eta_rho_levels(first_constant_r_rho_level) -        &
                eta_theta_levels(bl_levels) )                       &
              +  r_theta_levels(i,j,bl_levels)

          END DO
        END DO
      END DO

    CASE( height_gen_smooth )
! A smooth quadratic height generation
      DO k=1, first_constant_r_rho_level-1
        DO j=1, rows
          DO i=1, row_length
          r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +  &
           earth_radius + orog(i,j) * (1.0 - eta_rho_levels(k)        &
                /eta_rho_levels(first_constant_r_rho_level))**2

          r_theta_levels(i,j,k) = eta_theta_levels(k) *               &
               z_top_of_model + Earth_radius + Orog(i,j) *            &
               (1.0 - eta_theta_levels(k) /                           &
                eta_rho_levels(first_constant_r_rho_level))**2
          END DO
        END DO
      END DO

    CASE default
      ErrorStatus = 10
      WRITE (Cmessage,*) 'Unrecognised height generation method - ',&
                         'Dump needs to be reconfigured'

      CALL ereport( RoutineName, ErrorStatus, Cmessage )
  END SELECT

  RETURN

END SUBROUTINE calc_levels

