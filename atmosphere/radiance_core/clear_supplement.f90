! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate clear-sky fluxes.
!
! Method:
!   This subroutine is called after fluxes including clouds have
!   been calculated to find the corresponding clear-sky fluxes.
!   The optical properties of the column are already known.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE clear_supplement(ierr, n_profile, n_layer                    &
    , i_solver_clear                                                    &
    , trans_free, reflect_free, trans_0_free, source_coeff_free         &
    , isolir, flux_inc_direct, flux_inc_down                            &
    , s_down_free, s_up_free                                            &
    , albedo_surface_diff, albedo_surface_dir                           &
    , source_ground                                                     &
    , l_scale_solar, adjust_solar_ke                                    &
    , flux_direct_clear, flux_total_clear                               &
    , nd_profile, nd_layer, nd_source_coeff                             &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_source_coeff
!       Size allocated for layers

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir                                                            &
!       Spectral region
    , i_solver_clear
!       Solver for clear fluxes
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar beam
  REAL (RealK), INTENT(IN) ::                                           &
      trans_free(nd_profile, nd_layer)                                  &
!       Transmission coefficients
    , reflect_free(nd_profile, nd_layer)                                &
!       Reflection coefficients
    , trans_0_free(nd_profile, nd_layer)                                &
!       Direct transmission coefficients
    , source_coeff_free(nd_profile, nd_layer, nd_source_coeff)          &
!       Coefficients in source terms
    , albedo_surface_diff(nd_profile)                                   &
!       Diffuse albedo
    , albedo_surface_dir(nd_profile)                                    &
!       Direct albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , source_ground(nd_profile)                                         &
!       Ground source function
    , adjust_solar_ke(nd_profile, nd_layer)
!       Scaling of solar beam

  REAL (RealK), INTENT(INOUT) ::                                        &
      s_down_free(nd_profile, nd_layer)                                 &
!       Downward source
    , s_up_free(nd_profile, nd_layer)
!       Upward source

  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_clear(nd_profile, 0: nd_layer)                        &
!       Clear direct flux
    , flux_total_clear(nd_profile, 2*nd_layer+2)
!       Clear total fluxes


! Dummy variabales.
  INTEGER                                                               &
      n_equation
!       Number of equations
  REAL (RealK) ::                                                       &
      a5(nd_profile, 5, 2*nd_layer+2)                                   &
!       Pentadiagonal matrix
    , b(nd_profile, 2*nd_layer+2)                                       &
!       Rhs of matrix equation
    , work_1(nd_profile, 2*nd_layer+2)
!       Working array for solver

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=256)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'clear_supplement'


  IF (lhook) CALL dr_hook('CLEAR_SUPPLEMENT',zhook_in,zhook_handle)

! The source functions only need to be recalculated in the visible.
  IF (isolir == ip_solar) THEN
! DEPENDS ON: solar_source
    CALL solar_source(n_profile, n_layer                                &
      , flux_inc_direct                                                 &
      , trans_0_free, source_coeff_free                                 &
      , l_scale_solar, adjust_solar_ke                                  &
      , flux_direct_clear                                               &
      , s_down_free, s_up_free                                          &
      , nd_profile, nd_layer, nd_source_coeff                           &
      )
  END IF


! Select an appropriate solver for the equations of transfer.
  IF (i_solver_clear == ip_solver_pentadiagonal) THEN

!   Calculate the elements of the matrix equations.
! DEPENDS ON: set_matrix_pentadiagonal
    CALL set_matrix_pentadiagonal(n_profile, n_layer                    &
      , trans_free, reflect_free                                        &
      , s_down_free, s_up_free                                          &
      , albedo_surface_diff, albedo_surface_dir                         &
      , flux_direct_clear(1, n_layer), flux_inc_down                    &
      , source_ground                                                   &
      , a5, b                                                           &
      , nd_profile, nd_layer                                            &
      )
    n_equation=2*n_layer+2

! DEPENDS ON: band_solver
    CALL band_solver(n_profile, n_equation                              &
      , 2, 2                                                            &
      , a5, b                                                           &
      , flux_total_clear                                                &
      , work_1                                                          &
      , nd_profile, 5, 2*nd_layer+2                                     &
      )

  ELSE IF (i_solver_clear == ip_solver_homogen_direct) THEN

!   Solve for the fluxes in the column directly.
! DEPENDS ON: solver_homogen_direct
    CALL solver_homogen_direct(n_profile, n_layer                       &
      , trans_free, reflect_free                                        &
      , s_down_free, s_up_free                                          &
      , isolir, albedo_surface_diff, albedo_surface_dir                 &
      , flux_direct_clear(1, n_layer), flux_inc_down                    &
      , source_ground                                                   &
      , flux_total_clear                                                &
      , nd_profile, nd_layer                                            &
      )

  ELSE

     cmessage = '*** Error: The solver specified for clear-sky fluxes ' &
        //'is not valid.'
     ierr=i_err_fatal
     CALL ereport(RoutineName, ierr, cmessage)

  END IF


  IF (lhook) CALL dr_hook('CLEAR_SUPPLEMENT',zhook_out,zhook_handle)

END SUBROUTINE clear_supplement
