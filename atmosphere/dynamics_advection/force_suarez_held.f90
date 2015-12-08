! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Force_Suarez_Held

      Subroutine Force_Suarez_Held                                      &
     &                        (row_length, rows, model_levels           &
     &,                        timestep                                 &
     &,                        SuHe_newtonian_timescale_ka              &
     &,                        SuHe_newtonian_timescale_ks              &
     &,                        SuHe_pole_equ_deltaT                     &
     &,                        SuHe_static_stab                         &
     &,                        SuHe_level_weight                        &
     &,                        SuHe_sigma_cutoff                        &
     &,                        L_SH_Williamson, SuHe_relax              &
     &,                        cos_theta_latitude                       &
     &,                        sin_theta_latitude, off_x, off_y         &
     &,                        exner_theta_levels                       &
     &,                        p_theta_levels, p_star                   &
     &,                        theta, theta_star )

! Purpose:
!          Provides temperature relaxation for Suarez Held test problem.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.


      USE earth_constants_mod, ONLY: g
      USE atmos_constants_mod, ONLY: R, recip_kappa
      USE ereport_mod, ONLY: ereport
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Include physical constants

      Logical                                                           &
     &  L_SH_williamson

      Integer                                                           &
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
     &, SuHe_relax   ! Switch to choose form of relaxtion

      Real                                                              &
     &  timestep

      Real                                                              &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &         model_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &         model_levels)                                            &
     &, p_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
     &         model_levels)                                            &
     &, p_star(row_length, rows)                                        &
     &, cos_theta_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y)&
     &, sin_theta_latitude(row_length, rows)

! Suarez Held variables

      Real                                                              &
     &  SuHe_newtonian_timescale_ka                                     &
     &, SuHe_newtonian_timescale_ks                                     &
     &, SuHe_pole_equ_deltaT                                            &
     &, SuHe_static_stab                                                &
     &, SuHe_level_weight(model_levels)                                 &
     &, SuHe_sigma_cutoff

! local variables

      Integer                                                           &
     &  i, j, k

      Real                                                              &
        temp1                                                           &
      , temp2                                                           &
      , theta_eq(row_length, rows,model_levels)                         &
      , newtonian_timescale                                             

! Williamson variables
! 1. parameters. NB: changing these WILL change answers
      Real                                                              &
     &  p_d                                                             &
     &, p_pl                                                            &
     &, lapse_rate_d                                                    &
     &, lapse_rate_i                                                    &
     &, delta_phi_0                                                     &
     &, A_will                                                          &
     &, phi_0                                                           &
     &, T_0

! 2. derived variables
      Real                                                              &
     &  p_i                                                             &
     &, p_eq                                                            &
     &, power_d                                                         &
     &, power_i                                                         &
     &, latitude                                                        &
     &, p_lim

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER :: ErrorStatus
      CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'force_suarez_held'
      CHARACTER(LEN=80)             :: cmessage
      
! No External routines

! ----------------------------------------------------------------------
! Section 1.  Calculate temperature relaxation forcing.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FORCE_SUAREZ_HELD',zhook_in,zhook_handle)

      If (L_SH_williamson) Then
        p_d = 10000.0                 ! 100 hpa
        p_pl = 200.0                  ! limited to 2hpa
        delta_phi_0 = pi / 12.0       ! 15 degrees converted to radians
        A_will = 2.65 / delta_phi_0
        phi_0 = 60.0 * pi / 180.      ! 60 degrees converted to radians
        T_0 = 200.0                   ! temperature at 100 hpa in
                                      !Held-Suarez
        lapse_rate_d = 0.002          ! 2 K /km
        lapse_rate_i = -0.003345      ! 3.345 K /km

        p_eq = p_d
        power_d = R * lapse_rate_d / g
        power_i = R * lapse_rate_i / g

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              p_lim = p_theta_levels(i,j,k)

! restrict Williamson temperature to not change with height above 2hpa
! as values get very unphysical.
              If (p_lim  <   200.0) p_lim =200.0

              If (p_lim  <=  p_d) Then
                temp1 = T_0 * (p_lim / p_d) ** power_d

                latitude = acos(cos_theta_latitude(i,j))
                p_i = p_eq - (p_eq - p_pl) * 0.5 *                      &
     &                       (1+tanh(A_will*(abs(latitude)-phi_0)))

                If (p_lim  <=  p_i) Then
                  temp2 = T_0 *                                         &
     &                    ((p_lim / p_i) ** power_i -                   &
     &                     1.0)
                  temp1 = temp1 + temp2
                End If
! convert temperature to potential temperature
                theta_eq(i,j,k) = temp1 / exner_theta_levels(i,j,k)

              Else
! calculate Theta_eq.
                temp1 = 200. / exner_theta_levels(i,j,k)
                temp2 = 315. - SuHe_pole_equ_deltaT                     &
     &                    * sin_theta_latitude(i,j)                     &
     &                    * sin_theta_latitude(i,j)                     &
     &                    - SuHe_static_stab                            &
     &                    * log(exner_theta_levels(i,j,k))              &
     &                    * recip_kappa                                 &
     &                    * cos_theta_latitude(i,j)                     &
     &                    * cos_theta_latitude(i,j)

                theta_eq(i,j,k) = max(temp1, temp2)
              End If
            End Do
          End Do
        End Do
      Else
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
! calculate Theta_eq.
              temp1 = 200. / exner_theta_levels(i,j,k)
              temp2 = 315. - SuHe_pole_equ_deltaT                       &
     &                    * sin_theta_latitude(i,j)                     &
     &                    * sin_theta_latitude(i,j)                     &
     &                    - SuHe_static_stab                            &
     &                    * log(exner_theta_levels(i,j,k))              &
     &                    * recip_kappa                                 &
     &                    * cos_theta_latitude(i,j)                     &
     &                    * cos_theta_latitude(i,j)

              theta_eq(i,j,k) = max(temp1, temp2)
            End Do
          End Do
        End Do
      End If

      if(SuHe_relax  ==  1)then
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
! calculate relaxation term.
            newtonian_timescale = SuHe_newtonian_timescale_ka           &
     &                          + ( SuHe_newtonian_timescale_ks -       &
     &                              SuHe_newtonian_timescale_ka )       &
     &                          * cos_theta_latitude(i,j) ** 4          &
     &                          * SuHe_level_weight(k)

            theta_star(i,j,k) = theta_star(i,j,k) - timestep *          &
     &                                    newtonian_timescale *         &
     &                           (theta_star(i,j,k) - theta_eq(i,j,k))


          End Do
        End Do
      End Do

      elseif (SuHe_relax  ==  2)then

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
! calculate relaxation term.
! temp1 is SuHe_level_weight(k)
        temp1 = (p_theta_levels(i,j,k)/p_star(i,j) - SuHe_sigma_cutoff) &
     &                                      / (1.0 - SuHe_sigma_cutoff)
              newtonian_timescale = SuHe_newtonian_timescale_ka         &
     &                          + ( SuHe_newtonian_timescale_ks -       &
     &                              SuHe_newtonian_timescale_ka )       &
     &                          * cos_theta_latitude(i,j) ** 4          &
     &                          * max(0.0, temp1)

              theta_star(i,j,k) = theta_star(i,j,k) - timestep *        &
     &                                    newtonian_timescale *         &
     &                           (theta_star(i,j,k) - theta_eq(i,j,k))
            End Do
          End Do
        End Do

      else
        ErrorStatus = 1
        cmessage = "Value of SuHe_relax not supported - put correct "// &
                   "value of 1 or 2 in idealised namelist"
        CALL ereport(RoutineName, ErrorStatus, cmessage)
      Endif   !  SuHe_relax  ==  1
      IF (lhook) CALL dr_hook('FORCE_SUAREZ_HELD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Force_Suarez_Held

