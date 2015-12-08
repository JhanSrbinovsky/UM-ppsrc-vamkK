! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_1_setup

      Subroutine GCR_precon_1_setup(                                    &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_C3, HM_C4,                 &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     model_domain,                                &
     &                     weight_upper, weight_lower,                  &
     &                     a0, a1, factor,                              &
     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity,                                &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     halo_i,halo_j                                &
     &                     )

! Purpose:
!          Calculates constant terms used inpre-conditioning operator
!          applied to field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain

      Integer                                                           &
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i                      &
     &                 ,1-halo_j:rows+halo_j,0:model_levels)            &
     &, r_rho_levels (1-halo_i:row_length+halo_i                        &
     &                 ,1-halo_j:rows+halo_j,model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k                                                           &
     &, istat                                                           &
     &, j0                                                              &
     &, j1                                                              &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      Real                                                              &
     &  factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

!  parallel variables
      Integer                                                           &
     &  proc_row_group

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local arrays.

      Real                                                              &
     &  a2(row_length,rows,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!   External routines.

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_PRECON_1_SETUP',zhook_in,zhook_handle)

        if(at_extremity(PSouth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
        j0= 2
      Else
        j0 = 1
      End If
        if(at_extremity(PNorth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
        j1 = rows - 1
      Else
        j1 = rows
      End If

      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) Then
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
      Else If(model_domain  ==  mt_lam) Then
! Solve over interior points only.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
        if(at_extremity(PEast)) i_stop = row_length-1
        if(at_extremity(PWest)) i_start = 2
      Elseif (model_domain  ==  mt_cyclic_LAM) then
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
      End If

      Do k = 1, model_levels
        If ( k  ==  1) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0(i,j,k) = - a1(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else if (k  ==  model_levels) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        Else
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) *                 &
     &                                weight_upper(i,j,k) )
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - a1(i,j,k)  - HM_C4(i,j,k)
            End Do
          End Do
        End If

      End Do

! Include horizontal second derivative terms in matrix.

! Interior points - Common to all model_domains.

      Do k = 1, model_levels
        Do j = j0, j1
          Do i = i_start, i_stop
            a0(i,j,k) = a0(i,j,k) - FV_sec_theta_latitude(i,j) *        &
     &                  (HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) +            &
     &                   HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) +                &
     &                   HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) +                &
     &                   HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) )
          End Do
        End Do
      End Do

! Treatment of North/South boundary.

      If (Model_domain  ==  mt_Global )Then

        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_s_component(i,k)= - HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*              &
     &                                HM_Cyy2(i,rows-1,k)
            End Do
          End Do
        End If

        If(at_extremity(PSouth))then
          CALL global_2d_sums(sum_s_component, row_length, 1, 0, 0,     &
                              model_levels, sum_s, proc_row_group)

          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)            &
     &                   /global_row_length
            Do i = 1,row_length
              a0(i,1,k) = a0(i,1,k) + sum_s(k)
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0,     &
                              model_levels, sum_n, proc_row_group)

          Do k = 1, model_levels
            sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)         &
     &                   /global_row_length
            Do i = 1,row_length
              a0(i,rows,k) = a0(i,rows,k) + sum_n(k)
            End Do
          End Do
        End If

      End If

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

! for global parallel code a0 only valid on interior points , NOT halos
      Do j = j_start, j_stop
        Do i = i_start, i_stop
          a0(i,j,1) = 1./a0(i,j,1)
        End Do
      End Do

      Do k= 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            factor(i,j,k) = a2(i,j,k) * a0(i,j,k-1)
            a0(i,j,k) = 1./(a0(i,j,k) - factor(i,j,k)*a1(i,j,k-1))
          End Do
        End Do
      End Do

!     end of routine GCR_precon_1_setup

      IF (lhook) CALL dr_hook('GCR_PRECON_1_SETUP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_1_setup
