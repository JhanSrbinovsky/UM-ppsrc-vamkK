! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Force
!
      Subroutine IDL_Force(                                             &
     &                     row_length                                   &
     &,                    rows, model_levels, timestep                 &
     &,                    delta_lambda, delta_phi, model_domain        &
     &,                    lambda_half_width, phi_half_width            &
     &,                    p_max, p_top, p_bottom, p_arbitrary          &
     &,                    lambda_heat_centre, phi_heat_centre          &
     &,                    max_heat_per_day, newtonian_timescale        &
     &,                    SuHe_newtonian_timescale_ka                  &
     &,                    SuHe_newtonian_timescale_ks                  &
     &,                    SuHe_pole_equ_deltaT, SuHe_static_stab       &
     &,                    SuHe_level_weight, SuHe_sigma_cutoff         &
     &,                    L_SH_williamson, SuHE_relax                  &
     &,                    L_damp, L_geo_for                            &
     &,                    L_bomex                                      &
     &,                    DMPTIM, HDMP, ZDMP                           &
     &,                    u_geo, v_geo                                 &
     &,                    off_x, off_y, l_datastart                    &
     &,                    at_extremity                                 &
     &,                    exner_theta_levels, p, theta_star            &
     &,                    p_theta_levels, p_star                       &
     &,                    theta                                        &
     &,                    cos_theta_latitude,sin_theta_latitude        &
     &,                    cool_rate, theta_surface                     &
     &,                    timestep_number                              &
     &,                    max_model_levels, max_num_force_times        &
     &,                    q, q_star                                    &
     &,                    u, v, R_u, R_v, u_ref, v_ref                 &
     &,                    theta_ref, n_rows                            &
     &,                    q_ref                                        &
     &,                    eta_theta_levels, eta_rho_levels             &
     &,                    height_domain                                &
     &,                    global_row_length, global_rows               &
     &,                    gc_all_proc_group, nproc                     &
     &,                    tforce_option, qforce_option, uvforce_option &
     &,                    num_tforce_times, num_qforce_times           &
     &,                    num_uvforce_times                            &
     &,                    tforce_time_interval, qforce_time_interval   &
     &,                    uvforce_time_interval                        &
     &,                    tforce_data_modlev, qforce_data_modlev       &
     &,                    uforce_data_modlev, vforce_data_modlev       &
     &,                    r_rho_levels, r_theta_levels                 &
     &,                    halo_i, halo_j                               &
     &,                    q_inc_subs, th_inc_subs                      &
     &,                    q_inc_ls, th_inc_ls                          &
     &,                    u_inc_dmp, q_inc_dmp, th_inc_dmp             &
     &,                    v_inc_dmp                                    &
     &,                    f3_at_u, f3_at_v                             &
     &,                    L_physics, problem_number, L_force)


! Purpose: top level deck for idealised forcing
!
! Method:
!          Is described in ;
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
      USE c_idl_force_options_mod, ONLY: no_forcing, force_increment,   &
                                         force_relax, force_reset
      USE problem_mod, ONLY: standard, monsoon, dynamical_core,         &
                             idealised_problem, standard_namelist
      Implicit None


      Integer                                                           &
     &  model_domain                                                    &
     &, problem_number

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in v field
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
                         ! Size of small halo in i
     &, off_y                                                           &
                         ! Size of small halo in j.
     &, l_datastart(3)                                                  &
     &, global_row_length                                               &
                             ! Points per global row
     &, global_rows                                                     &
                             ! No of global (theta) rows
     &, gc_all_proc_group                                               &
                             ! Group identifier for all processors
     &, nproc                                                           &
                             ! Number of processors
     &, timestep_number                                                 &
                             ! Timestep number since start of run
     &, halo_i, halo_j


      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &         model_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &         model_levels)

      Real                                                              &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, p_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
     &         model_levels)                                            &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      model_levels)                                               &
     &, q_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &      model_levels)                                               &
     &, u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

      Real, Intent(In) ::                                               &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, f3_at_u (1-off_x:row_length+off_x,                              &
     &                      1-off_y:rows+off_y)                         &
     &, f3_at_v (1-off_x:row_length+off_x,                              &
     &                      1-off_y:n_rows+off_y)

      Real                                                              &
     & u_g(model_levels)                                                &
     &,v_g(model_levels)           
                        ! Geostrophic winds

      Real                                                              &
     &  u_ref(model_levels)                                             &
                              ! u profile for use in lbcs
     &, v_ref(model_levels)                                             &
                              ! u_adv profile for use in lbcs
     &, theta_ref(model_levels)                                         &
     &, q_ref(model_levels)                                             &
     &, p_star(row_length, rows)                                        &
     &, cos_theta_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y)&
     &, sin_theta_latitude(row_length, rows)                            &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  height_domain

! Forcing variables

      Integer                                                           &
     &  max_model_levels                                                &
     &, max_num_force_times                                             &
     &, tforce_option                                                   &
     &, qforce_option                                                   &
     &, uvforce_option                                                  &
     &, num_tforce_times                                                &
     &, num_qforce_times                                                &
     &, num_uvforce_times                                               &
     &, SuHe_relax

! Loop counters
      Integer                                                           &
     & i, j, k

      Real                                                              &
     &  tforce_time_interval                                            &
     &, qforce_time_interval                                            &
     &, uvforce_time_interval                                           &
     &, tforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, qforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, uforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, vforce_data_modlev(max_model_levels, max_num_force_times)


! Monsoon variables

      Real                                                              &
     &  lambda_half_width                                               &
     &, phi_half_width                                                  &
     &, p_max                                                           &
     &, p_top                                                           &
     &, p_bottom                                                        &
     &, p_arbitrary                                                     &
     &, lambda_heat_centre                                              &
     &, phi_heat_centre                                                 &
     &, max_heat_per_day                                                &
     &, newtonian_timescale

      Real :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      Real :: u_geo, v_geo         ! Geostrophic wind

      Logical                                                           &
     &  L_physics                                                       &
     &, L_force                                                         &
     &, L_bomex                                                         &
     &, L_damp                                                          &
     &, L_geo_for


      Real                                                              &
     &  timestep                                                        &
     &, cool_rate                                                       &
     &, theta_surface                                                   &
     &, const_SST

! Suarez Held variables

      Real                                                              &
     &  SuHe_newtonian_timescale_ka                                     &
     &, SuHe_newtonian_timescale_ks                                     &
     &, SuHe_pole_equ_deltaT                                            &
     &, SuHe_static_stab                                                &
     &, SuHe_level_weight(model_levels)                                 &
     &, SuHe_sigma_cutoff     
     

      Logical                                                           &
     &  L_SH_williamson


      ! Local variables

      Real per_sec_factor  ! factor converting data in namelist
                           ! to increment per second

      Real                                                              &
     & q_inc_subs(row_length, rows, model_levels)                       &
     &,th_inc_subs(row_length, rows, model_levels)                      &
                                  ! subsidence increments
     &,q_inc_ls(row_length, rows, model_levels)                         &
     &,th_inc_ls(row_length, rows, model_levels)                        &
                                  ! large scale increments
     &,u_inc_dmp(row_length, rows, model_levels)                        &
     &,q_inc_dmp(row_length, rows, model_levels)                        &
     &,th_inc_dmp(row_length, rows, model_levels) 
                                  !Damping incs
      Real                                                              &
     & v_inc_dmp(row_length, n_rows, model_levels)

      Real z_uv   ! Height of uv levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_FORCE',zhook_in,zhook_handle)
      If (problem_number  ==  monsoon ) Then

! DEPENDS ON: force_monsoon
        Call Force_monsoon(                                             &
     &                         row_length                               &
     &,                        rows, model_levels, timestep             &
     &,                        delta_lambda, delta_phi, model_domain    &
     &,                        lambda_half_width, phi_half_width        &
     &,                        p_max, p_top, p_bottom, p_arbitrary      &
     &,                        lambda_heat_centre, phi_heat_centre      &
     &,                        max_heat_per_day, newtonian_timescale    &
     &,                        off_x, off_y, l_datastart                &
     &,                        at_extremity                             &
     &,                        exner_theta_levels, p, theta_star )

      Else If (problem_number  ==  dynamical_core) Then

! DEPENDS ON: force_suarez_held
        Call Force_Suarez_Held(                                         &
     &                         row_length, rows, model_levels           &
     &,                        timestep                                 &
     &,                        SuHe_newtonian_timescale_ka              &
     &,                        SuHe_newtonian_timescale_ks              &
     &,                        SuHe_pole_equ_deltaT                     &
     &,                        SuHe_static_stab                         &
     &,                        SuHe_level_weight                        &
     &,                        SuHe_sigma_cutoff                        &
     &,                        L_SH_williamson, SuHe_relax              &
     &,                        cos_theta_latitude                       &
     &,                        sin_theta_latitude, off_x, off_y         &
     &,                        exner_theta_levels                       &
     &,                        p_theta_levels, p_star                   &
     &,                        theta, theta_star )

      Else If (problem_number  ==  idealised_problem) Then

!-----------------------------------------------------------------------
!                       Temperature forcing
!-----------------------------------------------------------------------

        !----------------------------------------------
        ! Increment (Add forcing increment to a field)
        !----------------------------------------------
        If (tforce_option  ==  force_increment) Then

          ! Convert from K/day to K/sec
          per_sec_factor = 1./86400.

! DEPENDS ON: idl_force_increment
          Call IDL_Force_Increment(                                     &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_tforce_times, tforce_time_interval               &
     &,            per_sec_factor                                       &
     &,            tforce_data_modlev, theta_star)


        !----------------------------------------------
        ! Relaxation (Newtonian relaxation to specified forcing data)
        !----------------------------------------------
        Else If (tforce_option  ==  force_relax) Then

! DEPENDS ON: idl_force_relax
          Call IDL_Force_Relax(                                         &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            off_x, off_y, global_row_length, global_rows         &
     &,            nproc, timestep, timestep_number                     &
     &,            max_model_levels, max_num_force_times                &
     &,            newtonian_timescale                                  &
     &,            num_tforce_times, tforce_time_interval               &
     &,            tforce_data_modlev, theta, theta_star)


        !----------------------------------------------
        ! Reset (Reset data to specified forcing data)
        !----------------------------------------------
!        Else If (tforce_option  ==  force_reset) Then
!
!          Call IDL_Force_Reset()

        End If  ! on tforce_option

!-----------------------------------------------------------------------
!                         Humidity Forcing
!-----------------------------------------------------------------------

        !----------------------------------------------
        ! Increment (Add forcing increment to a field)
        !----------------------------------------------
        If (qforce_option  ==  force_increment) Then

          ! No conversion of input data
          per_sec_factor = 1.
! DEPENDS ON: idl_force_increment
          Call IDL_Force_Increment(                                     &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_qforce_times, qforce_time_interval               &
     &,            per_sec_factor                                       &
     &,            qforce_data_modlev, q_star)


        !----------------------------------------------
        ! Relaxation (Newtonian relaxation to specified forcing data)
        !----------------------------------------------
        Else If (qforce_option  ==  force_relax) Then

! DEPENDS ON: idl_force_relax
          Call IDL_Force_Relax(                                         &
     &             row_length, rows, model_levels                       &
     &,            halo_i, halo_j                                       &
     &,            off_x, off_y, global_row_length, global_rows         &
     &,            nproc, timestep, timestep_number                     &
     &,            max_model_levels, max_num_force_times                &
     &,            newtonian_timescale                                  &
     &,            num_qforce_times, qforce_time_interval               &
     &,            qforce_data_modlev, q, q_star)


        !----------------------------------------------
        ! Reset (Reset data to specified forcing data)
        !----------------------------------------------
!        Else If (qforce_option  ==  force_reset) Then
!
!          Call IDL_Force_Reset()

        End If  ! on qforce_option

!-----------------------------------------------------------------------
!                     Horizontal Wind Forcing
!-----------------------------------------------------------------------

        !----------------------------------------------
        ! Increment (Add forcing increment to a field)
        !----------------------------------------------
        If (uvforce_option  ==  force_increment) Then

          ! No conversion of input data
          per_sec_factor = 1.
! DEPENDS ON: idl_force_increment
          Call IDL_Force_Increment(                                     &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            per_sec_factor                                       &
     &,            uforce_data_modlev, R_u)

! DEPENDS ON: idl_force_increment
          Call IDL_Force_Increment(                                     &
     &             row_length, n_rows, model_levels                     &
     &,            off_x, off_y                                         &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            per_sec_factor                                       &
     &,            vforce_data_modlev, R_v)

        !----------------------------------------------
        ! Relaxation (Newtonian relaxation to specified forcing data)
        !----------------------------------------------
        Else If (uvforce_option  ==  force_relax) Then

! DEPENDS ON: idl_force_relax
          Call IDL_Force_Relax(                                         &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            off_x, off_y, global_row_length, global_rows         &
     &,            nproc, timestep, timestep_number                     &
     &,            max_model_levels, max_num_force_times                &
     &,            newtonian_timescale                                  &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            uforce_data_modlev, u, R_u)

! DEPENDS ON: idl_force_relax
          Call IDL_Force_Relax(                                         &
     &             row_length, n_rows, model_levels                     &
     &,            off_x, off_y                                         &
     &,            off_x, off_y, global_row_length, global_rows         &
     &,            nproc, timestep, timestep_number                     &
     &,            max_model_levels, max_num_force_times                &
     &,            newtonian_timescale                                  &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            vforce_data_modlev, v, R_v)

        !----------------------------------------------
        ! Reset (Reset data to specified forcing data)
        !----------------------------------------------
!        Else If (uvforce_option  ==  force_reset) Then
!
!          Call IDL_Force_Reset()

        End If  ! on uvforce_option

      Do k=1, model_levels
        Do j=1, rows
          Do i=1, row_length
            th_inc_subs(i,j,k)=theta_star(i,j,k)
            q_inc_subs(i,j,k)=q_star(i,j,k)
          End Do
        End Do
      End Do

! subsidence and large scale forcing
! DEPENDS ON: idl_force_subs
          Call IDL_Force_subs(                                          &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            halo_i, halo_j                                       &
     &,            timestep                                             &
     &,            r_theta_levels, r_rho_levels                         &
     &,            L_bomex                                              &
     &,            q_star, theta_star                                   &
     &,            q, theta                                             &
     &,            th_inc_subs, q_inc_subs                              &
     &,            th_inc_ls, q_inc_ls )


! Call damping layer code
      If (L_damp) Then

        Do k=1, model_levels

          Do j=1, rows
            Do i=1, row_length
              u_inc_dmp(i,j,k)=R_u(i,j,k)
              th_inc_dmp(i,j,k)=theta_star(i,j,k)
              q_inc_dmp(i,j,k)=q_star(i,j,k)
            End Do
          End Do

          Do j=1, n_rows
            Do i=1, row_length
              v_inc_dmp(i,j,k)=R_v(i,j,k)
            End Do
          End Do

        End Do ! loop over k

          
! DEPENDS ON: idl_damp
        Call IDL_Damp(                                                  &
     &             row_length, rows, n_rows, model_levels               &
     &,            off_x, off_y                                         &
     &,            halo_i, halo_j                                       &
     &,            timestep, height_domain                              &
     &,            eta_theta_levels,  eta_rho_levels                    &
     &,            DMPTIM, HDMP, ZDMP                                   &
     &,            u, v, theta, q                                       &
     &,            u_ref, v_ref, theta_ref, q_ref                       &
     &,            R_u, R_v, theta_star, q_star)

        Do k=1, model_levels
          Do j=1, rows
            Do i=1, row_length
              u_inc_dmp(i,j,k)=R_u(i,j,k)-u_inc_dmp(i,j,k)
              th_inc_dmp(i,j,k)=theta_star(i,j,k)-th_inc_dmp(i,j,k)
              q_inc_dmp(i,j,k)=q_star(i,j,k)-q_inc_dmp(i,j,k)
            End Do
          End Do

          Do j=1, n_rows
            Do i=1, row_length
              v_inc_dmp(i,j,k)=R_v(i,j,k)-v_inc_dmp(i,j,k)
            End Do
          End Do

        End Do ! loop over k
      End If ! L_damp

      If  (L_geo_for) Then
        ! Geostrophic forcing

        Do k=1, model_levels
          u_g(k)=u_geo
          v_g(k)=v_geo

          ! Height of uv levels
          z_uv= r_rho_levels(1,1,k) - r_theta_levels(1,1,0)

          If (L_bomex) Then
            ! BOMEX Geostrophic forcing
            u_g(k)=-10.0 + z_uv * 1.8e-3
            v_g(k)=0.0
          End If ! L_bomex
         
          Do j=1, rows
            Do i=1, row_length
              R_u(i,j,k)=R_u(i,j,k) - f3_at_u(i,j)*v_g(k)*timestep
            End Do
          End Do

          Do j=1, n_rows
            Do i=1, row_length
              R_v(i,j,k)=R_v(i,j,k) + f3_at_v(i,j)*u_g(k)*timestep
            End Do
          End Do

        End Do  ! loop over model_levels
      End If  ! L_geo_for

      End if     !  problem number

! End of routine.
      IF (lhook) CALL dr_hook('IDL_FORCE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Force

