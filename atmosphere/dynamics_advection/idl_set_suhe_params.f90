! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Set_SuHe_params

      Subroutine IDL_Set_SuHe_params                                    &
     &                      (row_length, rows, model_levels             &
     &,                     off_x, off_y, halo_i, halo_j                &
     &,                     me, n_proc, at_extremity                    &
     &,                     l_datastart, all_proc_group                 &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     p, p_theta_levels                           &
     &,                     base_frictional_timescale                   &
     &,                     friction_level                              &
     &,                     SuHe_sigma_cutoff, SuHe_level_weight )

! Purpose:
!          Sets intitial data profiles for dynamical core
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.


      USE earth_constants_mod, ONLY: g, earth_radius

      USE atmos_constants_mod, ONLY: r

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Include physical constants

       Integer                                                          &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j               ! Size of halo in j direction.

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc                                                          &
                   ! Total number of processors
     &, all_proc_group                                                  &
                       ! Group identifier for all processors.
     &, l_datastart(2)       ! First gridpoints held by this processor

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)

! Primary Arrays used in all models
      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    model_levels+1)                                               &
     &, p_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
     &    model_levels+1)

! Suarez Held variables

      Real                                                              &
     &  SuHe_pole_equ_deltaT                                            &
     &, SuHe_static_stab                                                &
     &, SuHe_sigma_cutoff                                               &
     &, base_frictional_timescale                                       &
     &, friction_level(model_levels)                                    &
     &, SuHe_level_weight(model_levels)

! local variables

      Integer                                                           &
     &  i, j, k                                                         &
     &, info

      Real                                                              &
     &  temp1                                                           &
     &, delta_z                                                         &
     &, goverRT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External routines

! ----------------------------------------------------------------------
! Section 1. Initialise parameters
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('IDL_SET_SUHE_PARAMS',zhook_in,zhook_handle)
      temp1 = 280.0
      goverRT = g/(R*temp1)

! set level weights and friction_level - only do on pe0 then broadcast
      If (me  ==  0) Then
        Do k = 1, model_levels
          temp1 = ( exp ( goverRT *                                     &
     &                 (r_theta_levels(1,1,0) - r_theta_levels(1,1,k))) &
     &                 - SuHe_sigma_cutoff)  /                          &
     &                                       (1.0 - SuHe_sigma_cutoff)
          SuHe_level_weight(k) = max (0.0, temp1)
        End Do
        Do k = 1, model_levels
          temp1 = ( exp ( goverRT *                                     &
     &                 (r_theta_levels(1,1,0) - r_rho_levels(1,1,k)))   &
     &                 - SuHe_sigma_cutoff)  /                          &
     &                                       (1.0 - SuHe_sigma_cutoff)
          friction_level(k) = base_frictional_timescale *               &
     &                            max (0.0, temp1)
        End Do

! Print out model levels
        print*,' At grid-point 1,1 '
        print*,'Height fields (minus Earth_radius) and theta profile'
        Do k = model_levels, 1, -1
          temp1 = r_theta_levels(1,1,k) - r_theta_levels(1,1,k-1)
          print*,' Theta-Level ',k,' height '                           &
     &,    r_theta_levels(1,1,k) - Earth_radius, ' delta_z = ',temp1
        End Do

        k = 0
        print*,' Theta-Level ',k,' height '                             &
     &,    r_theta_levels(1,1,k) - Earth_radius

        k = model_levels
        temp1 = r_theta_levels(1,1,k) - r_rho_levels(1,1,k) +           &
     &        r_theta_levels(1,1,k) - Earth_radius

        k = model_levels + 1
        print*,' Level ',k,' Nominal rho-height ', temp1                &
     &,  ' Pressure = ',p(1,1,k)
        Do k = model_levels , 1, -1
          print*,' Rho-level ',k,' height '                             &
     &,    r_rho_levels(1,1,k) - Earth_radius                           &
     &,    ' Pressure = ',p(1,1,k)
        End Do

        Do k = 1, model_levels / 2
          print*,' SuHe_level_weight(',k,') = ',SuHe_level_weight(k)
        End Do
        Do k = 1, model_levels / 2
          print*,'friction_level(',k,') = ',                            &
     &                    friction_level(k)
        End Do
      End If  !(me  ==  0)

      If (n_proc  >   1) Then
        Call gcg_ssync(all_proc_group, info)
        Call gcg_rbcast(201, model_levels, 0, all_proc_group, info,     &
     &                    friction_level)
        Call gcg_rbcast(201, model_levels, 0, all_proc_group, info,     &
     &                    SuHe_level_weight)
      endif

      IF (lhook) CALL dr_hook('IDL_SET_SUHE_PARAMS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Set_SuHe_params
