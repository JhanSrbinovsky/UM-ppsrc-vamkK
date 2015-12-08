! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_LBC_Reset
!

      Subroutine IDL_LBC_Reset(                                         &
     &             Lbc_size_g, model_levels                             &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_force_times, force_time_interval                 &
     &,            force_data_modlev, fld)


! Purpose: Reset values to the forcing values to a field
!
! Method:  Interpolates idealised forcing data to the current time
!          Calculates the domain average of the field for each level
!          Calculates the difference between the domain avg profile and
!          the 1D forcing profile.
!          Set the LBCs to this increment
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None


! Variables with Intent (In)

      Integer, Intent(In) ::                                            &
     &  Lbc_size_g                                                      &
                            ! Number of points
     &, model_levels                                                    &
                            ! Number of model levels
     &, num_force_times                                                 &
                            ! Number of times in forcing array
     &, timestep_number                                                 &
                            ! Timestep number in run
     &, max_model_levels                                                &
                            ! Maximum number of model levels
     &, max_num_force_times ! Maximum number of forcing times

      Real, Intent(In) ::                                               &
     &  timestep                                                        &
                               ! Timestep interval (secs)
     &, force_time_interval    ! Interval between forcing data

      ! Forcing data on model levels
      Real, Intent(In) ::                                               &
     & force_data_modlev(max_model_levels, max_num_force_times)

! Variables with Intent (InOut)

      Real, Intent(InOut) ::                                            &
     & fld(lbc_size_g, model_levels)     ! Field increment


! Local variables

      Integer                                                           &
     &  i,k                                                             &
                   ! Loop counters
     &, t_before                                                        &
                   ! Time index in forcing array before current time
     &, t_after    ! Time index in forcing array after current time

      Real                                                              &
     &  fld_avg                                                         &
                              ! Global average of field on a level
     &, current_time                                                    &
                              ! Time since start of run in seconds
     &, weight                                                          &
                              ! Weight for linear time interpolation
     &, force_data_modlev_int ! Time interpolated forcing data

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------
! Setup time interpolation variables
!---------------------------------------------------------------------

      ! Setup time interpolation variables
      ! IF (lhook) CALL dr_hook('IDL_LBC_RESET',zhook_in,zhook_handle)
      current_time = timestep*timestep_number
      t_before     = INT(current_time/force_time_interval) + 1
      t_after      = t_before + 1

      ! Set t_after appropriately if end of forcing dataset
      ! coincides with end of model run
      If (t_after  ==  num_force_times+1) t_after = t_before

      ! Set linear interpolation weights
      weight   = 1.0 - MOD(current_time,force_time_interval)            &
     &                        / force_time_interval

!---------------------------------------------------------------------
! Loop over levels, calculate domain average
!---------------------------------------------------------------------

      Do k = 1, model_levels

        !-------------------------------------------------------------
        ! Interpolate forcing data to current time
        !-------------------------------------------------------------

        force_data_modlev_int = force_data_modlev(k,t_before)*weight    &
     &      + force_data_modlev(k,t_after)*(1-weight)

        !-------------------------------------------------------------
        ! Set term to forcing data
        !-------------------------------------------------------------

        Do i = 1, Lbc_size_g
          fld(i,k) = force_data_modlev_int
        End Do

      End Do  ! loop over levels


! End of routine.
      IF (lhook) CALL dr_hook('IDL_LBC_RESET',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_LBC_Reset
