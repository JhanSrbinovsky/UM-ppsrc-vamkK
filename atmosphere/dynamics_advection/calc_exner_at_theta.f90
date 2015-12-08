! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_Exner_at_theta

      Subroutine Calc_Exner_at_theta(                                   &
     &                      r_theta_levels, r_rho_levels, exner_rho,    &
     &                      row_length, rows, model_levels,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      exner_theta_levels, L_include_halos)

! Purpose:
!          Calculates exner pressure at theta levels.
!
! Method:
!          Linear interpolation in height of exner (at rho levels)
!
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE atm_fields_bounds_mod, ONLY:                                  &
           tdims_s, tdims_l, pdims_l, pdims_s

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
      LOGICAL                                                           &
     &  L_include_halos  ! If .TRUE. then include halo regions
                         ! when performing the calculations


      REAL, INTENT(IN) :: r_rho_levels                                  &
                                 (pdims_l%i_start:pdims_l%i_end,        &
                                  pdims_l%j_start:pdims_l%j_end,        &
                                  pdims_l%k_start:pdims_l%k_end)
                                  
      REAL, INTENT(IN) :: r_theta_levels                                &
                                 (tdims_l%i_start:tdims_l%i_end,        &
                                  tdims_l%j_start:tdims_l%j_end,        &
                                                0:tdims_l%k_end)           
      REAL, INTENT(IN) :: exner_rho                                     &
                                 (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end + 1)


! Arguments with Intent OUT. ie: Output variables.

      REAL ::              exner_theta_levels                           &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &,  i_start,i_end                                                  &
                         ! Loop bounds
     &,  j_start,j_end   ! Loop bounds

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate pressure at desired theta levels.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CALC_EXNER_AT_THETA',zhook_in,zhook_handle)
      IF (L_include_halos) THEN
        i_start=1-off_x
        i_end=row_length+off_x
        j_start=1-Off_y
        j_end=rows+Off_y
      ELSE
        i_start=1
        i_end=row_length
        j_start=1
        j_end=rows
      ENDIF

      Do k = 1, model_levels - 1
        Do j = j_start, j_end
          Do i = i_start, i_end
            exner_theta_levels(i,j,k) = ( exner_rho(i,j,k) *            &
     &                                (r_rho_levels(i,j,k+1) -          &
     &                                 r_theta_levels(i,j,k) ) +        &
     &                                exner_rho(i,j,k+1) *              &
     &                                (r_theta_levels(i,j,k) -          &
     &                                 r_rho_levels(i,j,k) ) ) /        &
     &                                (r_rho_levels(i,j,k+1) -          &
     &                                 r_rho_levels(i,j,k) )
          End Do
        End Do
      End Do

        k = model_levels

!AM  extra pressure level above top theta level is same height above
!AM  as is the pressure below - hence weights are 0.5 and there is
!AM  no need to store the r_rho_level for the extra pressure
      Do j = j_start, j_end
        Do i = i_start, i_end
            exner_theta_levels(i,j,k) = 0.5 *                           &
     &                      ( exner_rho(i,j,k) + exner_rho(i,j,k+1) )
        End Do
      End Do

! End of routine
      IF (lhook) CALL dr_hook('CALC_EXNER_AT_THETA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Calc_Exner_at_theta
