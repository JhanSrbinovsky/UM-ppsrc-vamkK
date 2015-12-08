! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE CAL_ENG_MASS_CORR--------------------------------------
!  
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              - TO CALCUATE THE NECESSARY CORRECTION TO
!                TEMPERATURE TO CONSERVE TOTAL ENERGY
!  
!    NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!  
!    DOCUMENTATION :
!  
!----------------------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Energy Correction

      SUBROUTINE cal_eng_mass_corr_4A (                                 &
                                    a_energysteps,                      &
                                    sum_energy_fluxes,                  &
                                    tot_mass_init,tot_energy_init,      &
                                    energy_corr,tot_energy_final)

!
      USE atmos_constants_mod,   ONLY: cp, r
      USE earth_constants_mod,   ONLY: earth_radius
      USE global_2d_sums_mod,    ONLY: global_2d_sums

      USE conversions_mod,       ONLY: pi
      USE yomhook,               ONLY: lhook, dr_hook
      USE parkind1,              ONLY: jprb, jpim
      USE UM_ParVars,            ONLY: proc_all_group=>gc_all_proc_group
      USE proc_info_mod,         ONLY: me

      USE timestep_mod,          ONLY: timestep
      USE atm_fields_bounds_mod, ONLY: pdims

      USE horiz_grid_mod,       ONLY : delta_lambda => delta_xi1,       &
                                       delta_phi    => delta_xi2
      IMPLICIT NONE

!

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      INTEGER ::  a_energysteps
                          ! number of timesteps per energy period
      REAL ::                                                           &
        sum_energy_fluxes(pdims%i_end,pdims%j_end) ! total energy flux
!                                                  ! into atmosphere
!
      REAL :: tot_energy_init      ! IN TOTAL ENERGY OF ATMOSPHERE
                                   ! AT START OF energy correction period
!
      REAL :: tot_mass_init        ! IN TOTAL dry MASS OF ATMOSPHERE
                                   ! AT START OF MODEL
!
      REAL :: tot_energy_final     ! IN TOTAL ENERGY OF ATMOSPHERE
                                   ! AT END OF energy corr period
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!
      REAL :: energy_corr          ! INOUT ENERGY CORRECTION FACTOR
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL :: chg_energy          ! CHANGE IN ENERGY
!
      REAL :: error_energy,                                             &
                               ! ERROR IN ENERGY CALCULATION
           factor,                                                      &
                               ! grid resolution factor
           tot_fluxes,                                                  &
                               ! global total energy flux
           flux_corr           ! flux correction
!
      REAL :: tot_fluxes_sum(1)  ! for summation
!----------------------------------------------------------------------
! INTERNAL LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER :: i,j                ! LOOP COUNTER

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!
!
!*---------------------------------------------------------------------
!
!======================================================================
! CALCULATION OF TEMPERATURE CHANGE TO BE ADDED OVER THE NEXT DAY
! TO CORRECT FOR ERROR IN ENERGY BUDGET DURING PRESENT DAY
!======================================================================
!

! Do global sum of heat fluxes added to the atmosphere during period
! This sum of fluxes is accumulated in D1(jnet_flux)
! Note no halo on input array, grid_type 1, halo_type 3

      IF (lhook) CALL dr_hook('CAL_ENG_MASS_CORR_4A',zhook_in,zhook_handle)
      tot_fluxes = 0.0

      CALL global_2d_sums(sum_energy_fluxes, pdims%i_end, pdims%j_end,  &
                          0, 0, 1, tot_fluxes_sum, proc_all_group)

      tot_fluxes = tot_fluxes_sum(1)

! add on energy correction added during last energy correction period

      factor     = delta_lambda*delta_phi
      tot_fluxes = tot_fluxes*factor

      tot_fluxes = tot_fluxes +                                         &
           (cp-r)*energy_corr*tot_mass_init*a_energysteps*timestep
!
!----------------------------------------------------------------------
! CALCULATE ENERGY CHANGE DURING THE energy period
!----------------------------------------------------------------------
!
      chg_energy = tot_energy_final - tot_energy_init
!
!
!----------------------------------------------------------------------
! CALCULATE ERROR = DIFFERENCE BETWEEN CHANGE IN TOTAL ENERGY
! DURING THE period AND THAT EXPECTED DUE TO THE FLUXES OF ENERGY INTO
! THE ATMOSPHERE
!-----------------------------------------------------------------------
!
      error_energy = tot_fluxes - chg_energy
!
!
!-----------------------------------------------------------------------
! CALCULATE TEMPERATURE CHANGE TO BE APPLIED OVER THE
! NEXT period TO CORRECT FOR THE ERROR
! Now use Cv instead of Cp
!-----------------------------------------------------------------------

      energy_corr = error_energy / ((cp-r)*tot_mass_init)
      flux_corr = error_energy/                                         &
      (4.*pi*earth_radius*earth_radius*a_energysteps*timestep)
!
!----------------------------------------------------------------------
! CALCULATE RATE OF TEMPERATURE CHANGE WHICH NEEDS TO BE
! APPLIED OVER NEXT PERIOD TO CORRECT FOR ERROR
!----------------------------------------------------------------------
!
      energy_corr = energy_corr / (a_energysteps*timestep)
!
!
!----------------------------------------------------------------------
! DIAGNOSTICS
!----------------------------------------------------------------------
!


      IF (me  ==  0)                                                    &
        WRITE(6, fmt='(7(1X,A,E13.5,A,/))')                             &
            'FINAL TOTAL ENERGY          = ',TOT_ENERGY_FINAL,   ' J/ ',&
            'INITIAL TOTAL ENERGY        = ',TOT_ENERGY_INIT,    ' J/ ',&
            'CHG IN TOTAL ENERGY O. P.   = ',CHG_ENERGY,         ' J/ ',&
            'FLUXES INTO ATM OVER PERIOD = ',TOT_FLUXES,         ' J/ ',&
            'ERROR IN ENERGY BUDGET      = ',ERROR_ENERGY,       ' J/ ',&
            'TEMP CORRECTION OVER A DAY  = ',ENERGY_CORR*86400.0,' K  ',&
            'TEMPERATURE CORRECTION RATE = ',ENERGY_CORR,        ' K/S',&
            'FLUX CORRECTION (ATM)       = ',flux_corr,          ' W/M2'

!----------------------------------------------------------------------
! Note dry mass correction applied in eng_mass_cal routine
!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CAL_ENG_MASS_CORR_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE CAL_ENG_MASS_CORR_4A
      
