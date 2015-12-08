! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_P_star

      Subroutine Calc_P_star(                                           &
     &                  r_theta_levels, r_rho_levels, p, rho,           &
     &                  row_length, rows, model_levels,                 &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  p_star)

! Purpose:
!          Calculates surface pressure.
!
! Method:
!          Is described in ;
!
!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE earth_constants_mod, ONLY: g

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y      ! Size of small halo in j.

! Physical constants

      REAL ::                                                           &
        rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
              model_levels)                                             &
      , p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
            model_levels+1)                                             &
      , r_theta_levels(1-halo_i:row_length+halo_i,                      &
                         1-halo_j:rows+halo_j,0:model_levels)           &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j, model_levels)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  p_star (row_length, rows)

! Local Variables.

      Integer                                                           &
     &  i, j         ! Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate pressure at surface.
! ----------------------------------------------------------------------

! Surface pressure calculated by hydrostatic extrapolation from bottom
! model level.

      IF (lhook) CALL dr_hook('CALC_P_STAR',zhook_in,zhook_handle)
      Do j = 1, rows
        Do i = 1, row_length
          p_star(i,j) = p(i,j,1) + g * rho(i,j,1)                       &
     &                               * (r_rho_levels(i,j,1) -           &
     &                                  r_theta_levels(i,j,0) ) /       &
     &                                 (r_rho_levels(i,j,1) *           &
     &                                  r_rho_levels(i,j,1) )

        End Do
      End Do

! End of routine
      IF (lhook) CALL dr_hook('CALC_P_STAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Calc_P_star
