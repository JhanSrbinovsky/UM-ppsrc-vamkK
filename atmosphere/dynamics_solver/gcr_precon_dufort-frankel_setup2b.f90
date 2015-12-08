! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_D_F_setup_2B
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
      SUBROUTINE GCR_precon_D_F_setup_2B(                               &
     &             row_length, rows, n_rows, model_levels, model_domain,&
     &             eta_theta_levels, eta_rho_levels,                    &
     &             FV_sec_theta_latitude, dlambda_p,                    &
     &             recip_dlamp, recip_dlamu, recip_dphip, recip_dphiv,  &
     &             HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,                  &
     &             HM_Czz, HM_Cz, HM_C3, HM_C4,                         &
     &             weight_upper, weight_lower,                          &
     &             offx, offy, halo_i, halo_j, at_extremity,            &
     &             i_start, i_stop, j_start, j_stop, j_begin, j_end,    &
     &             g_row_len, proc_row_group,                           &
     &             L_regular, GCR_ADI_pseudo_timestep,                  &
     &             HM_Cxx1_Cxx2, HM_Cyy1_Cyy2,                          &
     &             minus_diag_timestep, FV_sec_theta_timestep,          &
     &             a0_z_df, a1_z_df, factor_z_df)

! Purpose:
!      Setup for Implicit Dufort-Frankel type algorithm pre-conditioner
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE global_2d_sums_mod, ONLY: global_2d_sums
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent In. ie: Input variables.

      Integer, Intent(In) ::                                            &
     &  model_domain                                                    &
     &, row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, offx                                                            &
     &, offy                                                            &
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, proc_row_group

      Real, Intent(In) ::                                               &
     &  FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

! interpolation weights for moving between theta levels and rho levels
      Real, Intent(In) ::                                               &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

! Coefficients of the elliptic operator
      Real, Intent(In) ::                                               &
     &  HM_Cxx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Czz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_C3 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_C4 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! vertical co-ordinate information
      Real, Intent(In) ::                                               &
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

!  VarRes horizontal co-ordinate spacing.
      Real, Intent(In) ::                                               &
     &  recip_dlamp(1-halo_i : row_length + halo_i)                     &
     &, recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, dlambda_p(1-halo_i : row_length+halo_i)

      Real, Intent(In) ::                                               &
     &  GCR_ADI_pseudo_timestep

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! parallel variables   IN

! loop bounds set in PE_Helmholtz
      Integer, Intent(In) ::                                            &
     &  i_start, i_stop, j_start, j_stop, j_begin, j_end                &
     &, g_row_len

      Logical, Intent(In) ::                                            &
     &  L_regular
                         ! true for regular resolution

! Arguments with Intent Out.
      Real, Intent(Out) ::                                              &
     &  a0_z_df(row_length,rows,model_levels)                           &
     &, a1_z_df(row_length,rows,model_levels)                           &
     &, factor_z_df(row_length,rows,model_levels)

      Real, Intent(Out) ::                                              &
     &  HM_Cxx1_Cxx2                                                    &
     &       (1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_Cyy1_Cyy2                                                    &
     &       (1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, minus_diag_timestep(row_length,rows,model_levels)               &
     &, FV_sec_theta_timestep(row_length,rows)

! Local Variables.

      Integer                                                           &
     &  i,j,k,info

      Real                                                              &
     &  fact1, fact2, fact3

! Local arrays.

      Real                                                              &
     &  a2_z_df(row_length,rows,model_levels)                           &
     &, minus_diag_component(row_length,rows,model_levels)              &
     &, l_s_poles(row_length,model_levels)                              &
     &, sum_s(model_levels)                                             &
     &, l_n_poles(row_length,model_levels)                              &
     &, sum_n(model_levels)

      Real                                                              &
     &  recip_g_row_len

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     No External Routines:

!-----------------------------------------------------------------------
! Section 1.0 Setup
!-----------------------------------------------------------------------

      IF (lhook)                                                       &
       CALL dr_hook('GCR_PRECON_D_F_SETUP_2B',zhook_in,zhook_handle)
      If (L_regular) Then
        Do k = 1, model_levels
          Do j = j_begin - 1, j_end + 1
            Do i = i_start - 1, i_stop + 1
              HM_Cxx1_Cxx2(i,j,k) = HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k)
              HM_Cyy1_Cyy2(i,j,k) = HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i = i_start, i_stop
              minus_diag_component(i,j,k) =                             &
     &           (HM_Cxx1_Cxx2(i  ,j,k) + HM_Cxx1_Cxx2(i-1,j,k))        &
     &         + (HM_Cyy1_Cyy2(i,j  ,k) + HM_Cyy1_Cyy2(i,j-1,k))
            End Do
          End Do
        End Do
      Else
        Do k = 1, model_levels
          Do j = j_begin - 1, j_end + 1
            Do i = i_start - 1, i_stop + 1
              HM_Cxx1_Cxx2(i,j,k) = HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k)     &
     &               * recip_dlamp(i)
              HM_Cyy1_Cyy2(i,j,k) = HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k)     &
     &               * recip_dphip(i,j)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i = i_start, i_stop
              minus_diag_component(i,j,k) =                             &
     &           (HM_Cxx1_Cxx2(i  ,j,k) + HM_Cxx1_Cxx2(i-1,j,k))        &
     &               * recip_dlamu(i-1)                                 &
     &         + (HM_Cyy1_Cyy2(i,j  ,k) + HM_Cyy1_Cyy2(i,j-1,k))        &
     &               * recip_dphiv(i,j-1)
            End Do
          End Do
        End Do
      End If ! L_regular

! ----------------------------------------------------------------------
! Section 1.1 Poles in Global Model
! ----------------------------------------------------------------------

! average the value and add on constant term, note any other polar
! point will do as all values are the same.

      If (model_domain == mt_global) Then

        If ( L_regular ) Then
          recip_g_row_len = 1. / g_row_len
        Else
          recip_g_row_len = 1. / (2. * Pi)
        End If

        If (at_extremity(PSouth)) Then
          j = 1
          If( L_regular ) then
            Do k = 1, model_levels
              Do i = 1, row_length
                l_s_poles(i,k) = HM_Cyy1_Cyy2(i,j,k)
              End Do
            End Do
          Else !  variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                l_s_poles(i,k) = HM_Cyy1_Cyy2(i,j,k)                    &
     &              * dlambda_p(i) * recip_dphiv(i,j+1)
              End Do
            End Do
          End if ! L_regular

          CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_s,                      &
                              proc_row_group)

          j = 1
          Do k = 1, model_levels
            minus_diag_component(1,j,k) = sum_s(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              minus_diag_component(i,j,k) = minus_diag_component(1,j,k)
            End Do
          End Do
        End If !  at_extremity(PSouth)

        If (at_extremity(PNorth)) Then
          j = rows-1
          If( L_regular ) then
            Do k = 1, model_levels
              Do i = 1, row_length
                l_n_poles(i,k) = HM_Cyy1_Cyy2(i,j,k)
              End Do
            End Do
          Else !  variable resolution
            Do k = 1, model_levels
              Do i = 1, row_length
                l_n_poles(i,k) = HM_Cyy1_Cyy2(i,j,k)                    &
     &              * dlambda_p(i) * recip_dphiv(i,j+1)
              End Do
            End Do
          End if ! L_regular

          CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_n,                      &
                              proc_row_group)

          j = rows
          Do k = 1, model_levels
            minus_diag_component(1,j,k) = sum_n(k) * recip_g_row_len
! Copy answer at one point to all others
            Do i =2,row_length
              minus_diag_component(i,j,k) = minus_diag_component(1,j,k)
            End Do
          End Do
        End If  ! at_extremity(PNorth)

      End If ! model_domain == mt_global

      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            minus_diag_component(i,j,k) = minus_diag_component(i,j,k)   &
     &        * FV_sec_theta_latitude(i,j)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
! Section 1.2 Set up for vertical operator
!-----------------------------------------------------------------------
      Do k = 1, model_levels
        If (k == 1) Then
          fact1 = 1.d0/(eta_theta_levels(k) - eta_theta_levels(k-1))
          fact2 = 1.d0/(eta_rho_levels(k+1) - eta_rho_levels(k))
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1_z_df(i,j,k) = fact2 * (fact1 * HM_Czz(i,j,k) +         &
     &                         HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0_z_df(i,j,k) = - a1_z_df(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else If (k == model_levels) Then
          fact1 = 1.d0/(eta_theta_levels(k) - eta_theta_levels(k-1))
          fact3 = 1.d0/(eta_rho_levels(k) - eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a2_z_df(i,j,k) = fact3 * (fact1 * HM_Czz(i,j,k-1) -       &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0_z_df(i,j,k) = - a2_z_df(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        Else
          fact1 = 1.d0/(eta_theta_levels(k) - eta_theta_levels(k-1))
          fact2 = 1.d0/(eta_rho_levels(k+1) - eta_rho_levels(k))
          fact3 = 1.d0/(eta_rho_levels(k) - eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1_z_df(i,j,k) = fact2 * (fact1 * HM_Czz(i,j,k) +         &
     &                         HM_C3(i,j,k) * HM_Cz(i,j,k) *            &
     &                                weight_upper(i,j,k) )
              a2_z_df(i,j,k) = fact3 * (fact1 * HM_Czz(i,j,k-1) -       &
     &                       HM_C3(i,j,k) * HM_Cz(i,j,k-1) *            &
     &                                weight_lower(i,j,k) )
              a0_z_df(i,j,k) =                                          &
     &        - a2_z_df(i,j,k) - a1_z_df(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        End If

      End Do

      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            a0_z_df(i,j,k) = minus_diag_component(i,j,k)                &
     &                               * (1.d0 + GCR_ADI_pseudo_timestep) &
     &                      - a0_z_df(i,j,k) * GCR_ADI_pseudo_timestep
          End Do
        End Do
      End Do
      Do k = 1, model_levels - 1
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            a1_z_df(i,j,k) = -a1_z_df(i,j,k) * GCR_ADI_pseudo_timestep
          End Do
        End Do
      End Do
      Do k = 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            a2_z_df(i,j,k) = -a2_z_df(i,j,k) * GCR_ADI_pseudo_timestep
          End Do
        End Do
      End Do

      Do j = j_start, j_stop
        Do i = i_start, i_stop
          a0_z_df(i,j,1) = 1.d0/a0_z_df(i,j,1)
        End Do
      End Do

      Do k = 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            factor_z_df(i,j,k) = a2_z_df(i,j,k) * a0_z_df(i,j,k-1)
            a0_z_df(i,j,k) = 1.d0                                       &
     &       /(a0_z_df(i,j,k) - factor_z_df(i,j,k)*a1_z_df(i,j,k-1))
          End Do
        End Do
      End Do

      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            minus_diag_timestep(i,j,k) = minus_diag_component(i,j,k)    &
     &            * (1.d0 + GCR_ADI_pseudo_timestep)
          End Do
        End Do
      End Do

      Do j = j_start, j_stop
        Do i = i_start, i_stop
          FV_sec_theta_timestep(i,j) =                                  &
     &      GCR_ADI_pseudo_timestep * FV_sec_theta_latitude(i,j)
        End Do
      End Do

      IF (lhook)                                                       &
      CALL dr_hook('GCR_PRECON_D_F_SETUP_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_D_F_setup_2B
