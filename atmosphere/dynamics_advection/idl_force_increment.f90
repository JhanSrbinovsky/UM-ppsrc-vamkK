! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_FORCE_INCREMENT
!

      Subroutine IDL_Force_Increment(                                   &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_force_times, force_time_interval                 &
     &,            per_sec_factor                                       &
     &,            force_data_modlev                                    &
     &,            fld_incr)


! Purpose: To apply a forcing increment to a field
!
! Method:  Interpolates forcing data increment (rate of change)
!          to the current time and adds to a model field
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

      Integer                                                           &
     &  row_length                                                      &
                            ! Number of points on a row
     &, rows                                                            &
                            ! Number of rows
     &, model_levels                                                    &
                            ! Number of model levels
     &, off_x                                                           &
                            ! Size of halo in i
     &, off_y                                                           &
                            ! Size of halo in j.
     &, num_force_times                                                 &
                            ! Number of times in forcing array
     &, timestep_number                                                 &
                            ! Model timestep number in run
     &, max_model_levels                                                &
                            ! Max number of model levels
     &, max_num_force_times ! Max number of times in forcing data

      Real                                                              &
     &  timestep                                                        &
                            ! Length of timestep in seconds
     &, force_time_interval                                             &
                            ! Forcing data time interval
     &, per_sec_factor      ! Factor converting data in namelist
                            ! to increment per second

      ! Forcing data interpolated to model levels
      Real                                                              &
     &  force_data_modlev(max_model_levels, max_num_force_times)


! Variables with Intent (InOut)

      Real                                                              &
     &  fld_incr(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &         model_levels)     ! field increment


! Local variables

      Integer                                                           &
     &  i,j,k                                                           &
                   ! loop counters
     &, t_before                                                        &
                   ! time index in forcing array before current time
     &, t_after    ! time index in forcing array after current time

      Real                                                              &
     &  current_time                                                    &
                              ! since start of run in seconds
     &, weight                                                          &
                              ! weight for linear time interpolation
     &, force_data_modlev_int ! time interpolated forcing data

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_FORCE_INCREMENT',zhook_in,zhook_handle)

      ! Setup time interpolation variables
      current_time = timestep*timestep_number
      t_before     = INT(current_time/force_time_interval) + 1
      t_after      = t_before + 1

      ! Set t_after appropriately if end of forcing dataset
      ! coincides with end of model run
      If (t_after  ==  num_force_times+1) t_after = t_before

      ! Set linear interpolation weights
      weight   = 1.0 - MOD(current_time,force_time_interval)            &
     &                        / force_time_interval

      Do k = 1, model_levels

        ! Interpolate forcing data to current time
        force_data_modlev_int = force_data_modlev(k,t_before)*weight    &
     &      + force_data_modlev(k,t_after)*(1-weight)

        ! Convert forcing data to increment per timestep
        force_data_modlev_int = force_data_modlev_int                   &
     &                          * per_sec_factor * timestep

        ! Add increment to field
        Do j = 1,rows
          Do i = 1,row_length
            fld_incr(i,j,k) = fld_incr(i,j,k) + force_data_modlev_int
          End Do
        End Do

      End Do  ! loop over levels


! End of routine.
      IF (lhook) CALL dr_hook('IDL_FORCE_INCREMENT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Force_Increment

