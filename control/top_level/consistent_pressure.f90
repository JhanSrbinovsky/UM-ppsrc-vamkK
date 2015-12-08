! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Given the exner field at rho levels derive all other pressure fields
!
! Called by atm_step. 
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           r_theta_levels, r_rho_levels, rho,                     &
     &           p, pstar, p_theta_levels,exner_theta_levels)          

      USE earth_constants_mod, ONLY: g

      USE atmos_constants_mod, ONLY: kappa, recip_kappa, p_zero

      USE atm_fields_bounds_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE

      Integer, Intent(IN) :: row_length
      Integer, Intent(IN) :: rows
      Integer, Intent(IN) :: model_levels
      Integer, Intent(IN) :: offx        
      Integer, Intent(IN) :: offy        
      Integer, Intent(IN) :: halo_i      
      Integer, Intent(IN) :: halo_j      

      Real, Intent(IN) :: exner_rho_levels                              &
                                 (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end + 1)

      Real, Intent(IN) :: rho    (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end)

      Real, Intent(IN) :: r_rho_levels                                  &
                                 (pdims_l%i_start:pdims_l%i_end,        &
                                  pdims_l%j_start:pdims_l%j_end,        &
                                  pdims_l%k_start:pdims_l%k_end)

      Real, Intent(IN) :: r_theta_levels                                &
                                 (tdims_l%i_start:tdims_l%i_end,        &
                                  tdims_l%j_start:tdims_l%j_end,        &
                                                0:tdims_l%k_end)

      Real, Intent(OUT) :: p     (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end + 1)

      Real, Intent(OUT) :: pstar(row_length,rows)

      Real, Intent(OUT) :: p_theta_levels                               &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

      Real, Intent(OUT) :: exner_theta_levels                           &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 

      Integer :: i,j,k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Calculate pressure from Exner.
      IF (lhook) CALL dr_hook('CONSISTENT_PRESSURE',zhook_in,zhook_handle)
      Do k = pdims_s%k_start, pdims_s%k_end + 1
        Do j = pdims_s%j_start, pdims_s%j_end
!CDIR NODEP
          Do i = pdims_s%i_start, pdims_s%i_end
            p(i,j,k)= (exner_rho_levels(i,j,k)**recip_kappa) * p_zero
          End Do
        End Do
      End Do

! Halos updated 
! DEPENDS ON: swap_bounds
      call Swap_Bounds(P,                                               &
     &                 row_length, rows, model_levels+1,                &
     &                 offx, offy, fld_type_p, .false.)

! DEPENDS ON: calc_p_star
      Call Calc_P_star(                                                 &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   P, RHO,                                        &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   PSTAR )

! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta(                                         &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   EXNER_RHO_LEVELS,                              &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   EXNER_THETA_LEVELS, .FALSE.)

! Calculate pressure from Exner at theta levels.
! DEPENDS ON: calc_p_from_exner
      call Calc_P_from_Exner(                                           &
                         P_THETA_LEVELS,                                &
                         row_length, rows,                              &
                         tdims_s%k_end-tdims_s%k_start+1,               &
                         offx, offy,                                    &
                         EXNER_THETA_LEVELS,.FALSE.)

      IF (lhook) CALL dr_hook('CONSISTENT_PRESSURE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Consistent_pressure
