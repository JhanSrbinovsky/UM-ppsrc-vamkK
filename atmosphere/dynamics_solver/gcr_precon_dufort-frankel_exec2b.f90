! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_D_F_exec_2B
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
      SUBROUTINE GCR_precon_D_F_exec_2B(                                &
     &           r,                                                     &
     &           row_length, rows, n_rows, model_levels, model_domain,  &
     &           recip_dlamu, recip_dphiv, dlambda_p,                   &
     &           offx, offy, halo_i, halo_j, at_extremity, L_regular,   &
     &           i_start, i_stop, j_start, j_stop, j_begin, j_end,      &
     &           g_row_len, proc_row_group,                             &
     &           GCR_ADI_pseudo_timestep, GCR_n_ADI_pseudo_timesteps,   &
     &           HM_Cxx1_Cxx2, HM_Cyy1_Cyy2,                            &
     &           minus_diag_timestep, FV_sec_theta_timestep,            &
     &           a0_z_df, a1_z_df, factor_z_df,                         &
     &           Soln)


! Purpose:
!      Calculates Implicit Dufort-Frankel type algorithm pre-conditioner
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

      USE global_2d_sums_mod, ONLY : global_2d_sums
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
     &  r(row_length,rows,model_levels)

!  VarRes horizontal co-ordinate spacing.
      Real, Intent(In) ::                                               &
     &  recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, dlambda_p(1-halo_i : row_length+halo_i)

      Real, Intent(In) ::                                               &
     &  a0_z_df(row_length,rows,model_levels)                           &
     &, a1_z_df(row_length,rows,model_levels)                           &
     &, factor_z_df(row_length,rows,model_levels)

      Real, Intent(In) ::                                               &
     &  HM_Cxx1_Cxx2                                                    &
     &       (1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_Cyy1_Cyy2                                                    &
     &       (1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, minus_diag_timestep(row_length,rows,model_levels)               &
     &, FV_sec_theta_timestep(row_length,rows)

      Real, Intent(In) ::                                               &
     &  GCR_ADI_pseudo_timestep

      Integer, Intent(In) ::                                            &
     &  GCR_n_ADI_pseudo_timesteps
                                     ! Number of pseudo timesteps
                                     ! to perform.
                                     
! loop bounds set in PE_Helmholtz
      Integer, Intent(In) ::                                            &
     &  i_start, i_stop, j_start, j_stop, j_begin, j_end                &
     &, g_row_len

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular
                         ! true for regular resolution

! Arguments with Intent Out.
      Real, Intent(Out) ::                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                ! SOLUTION.

! Local Variables.

      Integer                                                           &
     &  i, j, k, pseudo_timestep_number, info

! Local arrays.
      Real                                                              &
     &  Lh_of_Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)&
     &, RHS(row_length,rows,model_levels)                               &
     &, l_s_poles(row_length,model_levels)                              &
     &, sum_s(model_levels)                                             &
     &, l_n_poles(row_length,model_levels)                              &
     &, sum_n(model_levels)

      Real                                                              &
     &  recip_g_row_len

      Logical :: L_first_print_parm = .true.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     No External Routines:

!-----------------------------------------------------------------------
! Section 0 Setup
!-----------------------------------------------------------------------
      IF (lhook)                                  &
        CALL dr_hook('GCR_PRECON_D_F_EXEC_2B',zhook_in,zhook_handle)
      If (L_first_print_parm) Then
        Write(6,*)'GCR precon Dufort-Frankel parameters'
        Write(6,*)'GCR_ADI_pseudo_timestep   =', &
     &      GCR_ADI_pseudo_timestep
        Write(6,*)'GCR_n_ADI_pseudo_timesteps=', &
     &      GCR_n_ADI_pseudo_timesteps
        L_first_print_parm = .false.
      End If

!-----------------------------------------------------------------------
! loop start
!-----------------------------------------------------------------------
      Do pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

        If (pseudo_timestep_number == 1) Then

          Do k = 1, model_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                Soln(i,j,k) = 0.0d0
              End Do
            End Do
          End Do

          Do k = 1, model_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
                Lh_of_Soln(i,j,k) = 0.0d0
              End Do
            End Do
          End Do

        Else

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     Soln,row_length,rows,model_levels,           &
     &                     offx, offy, fld_type_p, .false.)

          If (model_domain == mt_lam) Then
            If (at_extremity(PWest)) Then
              Do k = 1, model_levels
                Do j = j_start, j_stop
                  Soln(1,j,k) = Soln(i_start,j,k)
                End Do
              End Do
            End If
            If (at_extremity(PEast)) Then
              Do k = 1, model_levels
                Do j = j_start, j_stop
                  Soln(row_length,j,k) = Soln(i_stop,j,k)
                End Do
              End Do
            End If
            If (at_extremity(PSouth)) Then
              Do k = 1, model_levels
                Do i = 1, row_length
                  Soln(i,1,k) = Soln(i,j_start,k)
                End Do
              End Do
            End If
            If (at_extremity(PNorth)) Then
              Do k = 1, model_levels
                Do i = 1, row_length
                  Soln(i,rows,k) = Soln(i,j_stop,k)
                End Do
              End Do
            End If
          End If ! model_domain

! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(Soln,row_length,rows,                &
     &                             model_levels,offx,offy)

! ----------------------------------------------------------------------
! Section 1.1 Calculate horizontal derivative terms
! ----------------------------------------------------------------------
          If (L_regular) Then
            Do k = 1, model_levels
              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  Lh_of_Soln(i,j,k) =                                   &
     &              (                                                   &
     &                HM_Cxx1_Cxx2(i  ,j,k)*(Soln(i+1,j,k)-Soln(i,j,k)) &
     &              - HM_Cxx1_Cxx2(i-1,j,k)*(Soln(i,j,k)-Soln(i-1,j,k)) &
     &              )                                                   &
     &            + (                                                   &
     &                HM_Cyy1_Cyy2(i,j  ,k)*(Soln(i,j+1,k)-Soln(i,j,k)) &
     &              - HM_Cyy1_Cyy2(i,j-1,k)*(Soln(i,j,k)-Soln(i,j-1,k)) &
     &              )
                End Do
              End Do
            End Do
          Else
            Do k = 1, model_levels
              Do j = j_begin, j_end
                Do i = i_start, i_stop
                  Lh_of_Soln(i,j,k) =                                   &
     &              (                                                   &
     &                HM_Cxx1_Cxx2(i  ,j,k)*(Soln(i+1,j,k)-Soln(i,j,k)) &
     &              - HM_Cxx1_Cxx2(i-1,j,k)*(Soln(i,j,k)-Soln(i-1,j,k)) &
     &              ) * recip_dlamu(i-1)                                &
     &            + (                                                   &
     &                HM_Cyy1_Cyy2(i,j  ,k)*(Soln(i,j+1,k)-Soln(i,j,k)) &
     &              - HM_Cyy1_Cyy2(i,j-1,k)*(Soln(i,j,k)-Soln(i,j-1,k)) &
     &              ) * recip_dphiv(i,j-1)
                End Do
              End Do
            End Do
          End If ! L_regular

! ----------------------------------------------------------------------
! Section 1.2 Poles in Global Model
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
                    l_s_poles(i,k) = (Soln(i,j+1,k) - Soln(i,j,k))      &
     &                  * HM_Cyy1_Cyy2(i,j,k)
                  End Do
                End Do
              Else !  variable resolution
                Do k = 1, model_levels
                  Do i = 1, row_length
                    l_s_poles(i,k) = (Soln(i,j+1,k) - Soln(i,j,k))      &
     &                  * HM_Cyy1_Cyy2(i,j,k)                           &
     &                  * dlambda_p(i) * recip_dphiv(i,j+1)
                  End Do
                End Do
              End if ! L_regular

              CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,       &
                                  model_levels, sum_s,                  &
                                  proc_row_group)

              j = 1
              Do k = 1, model_levels
                Lh_of_Soln(1,j,k) = +sum_s(k) * recip_g_row_len
! Copy answer at one point to all others
                Do i =2,row_length
                  Lh_of_Soln(i,j,k) = Lh_of_Soln(1,j,k)
                End Do
              End Do
            End If !  at_extremity(PSouth)

            If (at_extremity(PNorth)) Then
              j = rows-1
              If( L_regular ) then
                Do k = 1, model_levels
                  Do i = 1, row_length
                    l_n_poles(i,k) = (Soln(i,j+1,k) - Soln(i,j,k))      &
     &                  * HM_Cyy1_Cyy2(i,j,k)
                  End Do
                End Do
              Else !  variable resolution
                Do k = 1, model_levels
                  Do i = 1, row_length
                    l_n_poles(i,k) = (Soln(i,j+1,k) - Soln(i,j,k))      &
     &                  * HM_Cyy1_Cyy2(i,j,k)                           &
     &                  * dlambda_p(i) * recip_dphiv(i,j+1)
                  End Do
                End Do
              End if ! L_regular

              CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,       &
                                  model_levels, sum_n,                  &
                                  proc_row_group)

              j = rows
              Do k = 1, model_levels
                Lh_of_Soln(1,j,k) = -sum_n(k) * recip_g_row_len
! Copy answer at one point to all others
                Do i =2,row_length
                  Lh_of_Soln(i,j,k) = Lh_of_Soln(1,j,k)
                End Do
              End Do
            End If  ! at_extremity(PNorth)

          End If ! model_domain == mt_global

        End If ! pseudo_timestep_number == 1

! ----------------------------------------------------------------------
! Section 1.3   Set Right Hand Side
! ----------------------------------------------------------------------
        Do k = 1, model_levels
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              RHS(i,j,k) = (Lh_of_Soln(i,j,k) - r(i,j,k))               &
     &           * FV_sec_theta_timestep(i,j)                           &
     &           + Soln(i,j,k) * minus_diag_timestep(i,j,k)
            End Do
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2  Solve vertical equation
! ----------------------------------------------------------------------

        Do k = 1, model_levels
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              Soln(i,j,k) = RHS(i,j,k)
            End Do
          End Do
        End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

        Do k= 2, model_levels
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              Soln(i,j,k) = Soln(i,j,k) -                               &
     &                       factor_z_df(i,j,k)*Soln(i,j,k-1)
            End Do
          End Do
        End Do

! Back substitute to get solution.

        Do j = j_start, j_stop
          Do i = i_start, i_stop
            Soln(i,j,model_levels)  = a0_z_df(i,j,model_levels) *       &
     &                                Soln(i,j,model_levels)
          End Do
        End Do

        Do k= model_levels-1, 1, -1
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              Soln(i,j,k)  = a0_z_df(i,j,k) *                           &
     &      ( Soln(i,j,k) - a1_z_df(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do

      End Do ! pseudo_timestep_number

      IF (lhook)                              &
       CALL dr_hook('GCR_PRECON_D_F_EXEC_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_D_F_exec_2B
