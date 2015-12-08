! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Force_Lbc
!

      Subroutine IDL_Force_Lbc(                                         &
     &             row_length, rows, off_x, off_y                       &
     &,            halo_i, halo_j,  LENRIMA                             &
     &,            timestep, timestep_number                            &
     &,            model_levels, wet_levels                             &
     &,            max_model_levels, max_num_force_times                &
     &,            u_lbc, v_lbc, theta_lbc                              &
     &,            q_lbc, u_adv_lbc, v_adv_lbc                          &
     &,            exner_rho_levels_LBC                                 &
     &,            r_theta_levels, r_rho_levels                         &
     &,            eta_theta_levels, eta_rho_levels                     &
     &,            height_domain, theta_surface                         &
     &,            pforce_option                                        &
     &,            tforce_option, qforce_option, uvforce_option         &
     &,            num_pforce_times                                     &
     &,            num_tforce_times, num_qforce_times                   &
     &,            num_uvforce_times                                    &
     &,            pforce_time_interval                                 &
     &,            tforce_time_interval, qforce_time_interval           &
     &,            uvforce_time_interval                                &
     &,            p_surface_data                                       &
     &,            tforce_data_modlev, qforce_data_modlev               &
     &,            uforce_data_modlev, vforce_data_modlev               &
     &,            newtonian_timescale )


      USE earth_constants_mod, ONLY: g, earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE c_idl_force_options_mod, ONLY: no_forcing, force_increment,   &
                                         force_relax, force_reset
      IMPLICIT NONE

! Purpose: To apply idealised LBC forcing.
!
! Method: Resets LBC data to profile data timeseries specified in
!         idealised namelist.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions
!   This code is written to UMDP3 programming standards.

! Parameters required for dimensioning some of the arguments

! Variables with Intent (In)


      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, off_x                                                           &
                         ! Size of small halo in i
     &, off_y                                                           &
                         ! Size of small halo in j.
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_levels                                                      &
     &, num_pforce_times                                                &
                            ! No. of times in P forcing data
     &, num_tforce_times                                                &
                            ! No. of times in T forcing data
     &, num_qforce_times                                                &
                            ! No. of times in Q forcing data
     &, num_uvforce_times                                               &
                            ! No. of times in UV forcing data
     &, pforce_option                                                   &
     &, tforce_option                                                   &
     &, qforce_option                                                   &
     &, uvforce_option                                                  &
     &, timestep_number                                                 &
                             ! Model timestep number in run&,
     &, max_model_levels                                                &
                             ! Max number of model levels
     &, max_num_force_times                                             &
                             ! Max number of times in forcing data
     &, lenrima(Nfld_max,NHalo_max)
                             ! IN : Size of single level of LBC

      Real, Intent(In) ::                                               &
     &  timestep                                                        &
                             ! Length of timestep in seconds
     &, pforce_time_interval                                            &
                             ! Forcing data time interval
     &, tforce_time_interval                                            &
                             ! Forcing data time interval
     &, qforce_time_interval                                            &
     &, uvforce_time_interval                                           &
     &, height_domain                                                   &
     &, theta_surface                                                   &
     &, newtonian_timescale  ! Relaxation timescale

      ! Forcing data interpolated to model levels
      Real, Intent(In) ::                                               &
     &  p_surface_data(max_num_force_times)                             &
     &, tforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, qforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, uforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, vforce_data_modlev(max_model_levels, max_num_force_times)

! Variables with Intent (InOut)

      Real, Intent(InOut) ::                                            &
     &  u_lbc(lenrima(fld_type_u,halo_type_extended),                   &
     &            model_levels)                                         &
     &, v_lbc(lenrima(fld_type_v,halo_type_extended),                   &
     &            model_levels)                                         &
     &, theta_lbc(lenrima(fld_type_p,halo_type_extended),               &
     &            model_levels)                                         &
     &, q_lbc(lenrima(fld_type_p,halo_type_extended),                   &
     &        wet_levels)                                               &
     &, u_adv_lbc(lenrima(fld_type_u,halo_type_extended),               &
     &            model_levels)                                         &
     &, v_adv_lbc(lenrima(fld_type_v,halo_type_extended),               &
     &            model_levels)                                         &
     &, exner_rho_levels_lbc(lenrima(fld_type_p,halo_type_extended),    &
     &            model_levels+1)                                       &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, r_theta_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j, &
     &         0: model_levels)                                         &
     &, r_rho_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,   &
     &         model_levels)


! Local variables
      Integer i, k
      Integer Lbc_sizes_u, Lbc_sizes_v, Lbc_sizes_p

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_FORCE_LBC',zhook_in,zhook_handle)
      Lbc_sizes_p = lenrima(fld_type_p,halo_type_extended)
      Lbc_sizes_u = lenrima(fld_type_u,halo_type_extended)
      Lbc_sizes_v = lenrima(fld_type_v,halo_type_extended)

      !----------------------------------------------
      !              Temperature forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_p, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_tforce_times, tforce_time_interval               &
     &,            tforce_data_modlev, theta_lbc )

      !----------------------------------------------
      !              Humidity forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_p, wet_levels                              &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_qforce_times, qforce_time_interval               &
     &,            qforce_data_modlev, q_lbc )

!      !----------------------------------------------
!      !              Exner Pressure
!      !   (Reset data to specified forcing data)
!      !----------------------------------------------
!      ! Not yet included
!
!        Call IDL_LBC_Pressure(
!     &             row_length, rows, halo_i, halo_j
!     &,            Lbc_sizes_p
!     &,            model_levels, wet_levels
!     &,            theta_surface
!     &,            r_theta_levels, r_rho_levels
!     &,            eta_theta_levels, eta_rho_levels
!     &,            theta_lbc, q_lbc
!     &,            timestep, timestep_number
!     &,            max_num_force_times
!     &,            num_pforce_times, pforce_time_interval
!     &,            p_surface_data
!     &,            exner_rho_levels_LBC )


      !----------------------------------------------
      !           Horizontal Wind Forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_u, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            uforce_data_modlev, u_lbc )

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_v, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            vforce_data_modlev, v_lbc )

      !----------------------------------------------
      !             Horizontal Wind (_adv)
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_u, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            uforce_data_modlev, u_adv_lbc )

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_v, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            vforce_data_modlev, v_adv_lbc )


! End of routine.

      IF (lhook) CALL dr_hook('IDL_FORCE_LBC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Force_Lbc
