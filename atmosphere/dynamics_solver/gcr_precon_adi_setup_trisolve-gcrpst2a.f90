! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_setup_trisol

      Subroutine GCR_precon_ADI_setup_trisol(j_start,j_end,             &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_C3, HM_C4,                 &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     model_domain, GCR_precon_option,             &
     &                     ADI_pseudo_timestep,                         &
     &                     rescale,                                     &
     &                     weight_upper, weight_lower,                  &
     &                     a0_z, a1_z, factor_z,                        &
     &                     a0_za, a1_za, factor_za,                     &
     &                     a0_y, a1_y, factor_y,                        &

! new terms
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1, bv_soln_1_term2,            &

     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity, me,                            &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     proc_col_group,halo_i,halo_j)

! Purpose:
!          Calculates ADI pre-conditioning matrix coefficients.
!          This code only for global model at the moment.
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


      USE global_2d_sums_mod, ONLY: global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      
      USE precon_constants_mod, ONLY:                                   &
          no_precon,vert_precon,vert_plus_xyz_ADI_precon,xyz_ADI_precon,&
          vert_plus_xz_ADI_precon,xz_ADI_precon,Dufort_Frankel_precon
      
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer,intent(in)                            :: j_start,j_end

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain                                                    &
     &, GCR_precon_option  ! 0 = no preconditioning
                           ! 1 = Vertical block pre-conditioner
                           ! 2 = 1 iteration of 1 followed by 3D ADI
                           ! 3 = 3D ADI only
                           ! 4 = 1 iteration of 1 followed by xz ADI
                                ! 5 = xz ADI only

      Integer                                                           &
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j,me

      Integer                                                           &
     &  proc_row_group                                                  &
     &, proc_col_group

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
     &  ADI_pseudo_timestep

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, rescale (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
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
     &  a0_za(row_length,rows,model_levels)                             &
     &, a1_za(row_length,rows,model_levels)                             &
     &, factor_za(row_length,rows,model_levels)                         &
     &, a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)

      Real                                                              &
     &  a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)

      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &

     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k, istat

      Real                                                              &
     &  factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3                                                        &
     &, recip_ADI_pseudo_timestep

      Integer                                                           &
     &  i_start                                                         &
     &, i_end

! Local arrays.

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      Real                                                              &
     &  a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  a2_za(row_length,rows,model_levels)                             &
     &, a2_y(row_length,rows,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!    External routines.
      External                                                          &
     &  mpp_tri_solve_setup

!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_SETUP_TRISOL',zhook_in,zhook_handle)
      recip_ADI_pseudo_timestep = 1./ ADI_pseudo_timestep

!kk   Bring in j_start and j_end as arguments

!-----------------------------------------------------------------------
!     Section 2. Set-up matrix coefficients in x direction.
!-----------------------------------------------------------------------

      If (model_domain  ==  mt_Global) Then
        i_start = 1
        i_end = row_length
        If (at_extremity(PWest) ) Then
          i_start = 2
        End If
        If (at_extremity(PEast) ) Then
          i_end = row_length-2
        End If

! Form matrix

! Calculate coefficients
        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = 1, row_length
              a_plus_x(i,j,k) = HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) *         &
     &                      FV_sec_theta_latitude(i,j)
              a_minus_x(i,j,k) = HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) *    &
     &                      FV_sec_theta_latitude(i,j)
              a_central_x(i,j,k) = - a_plus_x(i,j,k) - a_minus_x(i,j,k) &
     &                      - recip_ADI_pseudo_timestep                 &
     &                       * rescale(i,j,k)                           &
     &                      - HM_C4(i,j,k)
            End Do
          End Do
        End Do ! end loop over levels

! set any coefficients not covered by general loop
        Do k = 1, model_levels
          Do j = 1, j_start-1
            Do i = 1, row_length
              a_central_x(i,j,k) = 0.0
              a_plus_x(i,j,k) = 0.0
              a_minus_x(i,j,k) = 0.0
            End Do
          End Do
        End Do ! end loop over levels
        Do k = 1, model_levels
          Do j = j_end+1, rows
            Do i = 1, row_length
              a_central_x(i,j,k) = 0.0
              a_plus_x(i,j,k) = 0.0
              a_minus_x(i,j,k) = 0.0
            End Do
          End Do
        End Do ! end loop over levels

! DEPENDS ON: mpp_tri_solve_setup
        Call mpp_tri_solve_setup(                                       &
     &                     row_length, rows, model_levels,              &
     &                     j_start, j_end,                              &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_central_x, a_plus_x, a_minus_x,            &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2)

!          Else
! Limited area not implemented so ignored.
      End If

!-----------------------------------------------------------------------
!     Section 3. Setup matrix coefficients in vertical direction.
!-----------------------------------------------------------------------

! Form matrix

      Do k = 1, model_levels
        If ( k  ==  1) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                       eta_rho_levels(k))
          Do j = 1, rows
            Do i = 1, row_length
              a1_za(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +     &
     &                      HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0_za(i,j,k) = - a1_za(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else if (k  ==  model_levels) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                    eta_rho_levels(k-1))
          Do j = 1, rows
            Do i = 1, row_length
              a2_za(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -   &
     &                      HM_C3(i,j,k) * HM_Cz(i,j,k-1) *             &
     &                                weight_lower(i,j,k) )
              a0_za(i,j,k) = - a2_za(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        Else
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = 1, rows
            Do i = 1, row_length
              a1_za(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +     &
     &                       HM_C3(i,j,k) * HM_Cz(i,j,k) *              &
     &                                weight_upper(i,j,k) )
              a2_za(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -   &
     &                       HM_C3(i,j,k) * HM_Cz(i,j,k-1) *            &
     &                                weight_lower(i,j,k) )
              a0_za(i,j,k) = -a2_za(i,j,k)- a1_za(i,j,k)- HM_C4(i,j,k)
            End Do
          End Do
        End If
      End Do

      If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.          &
     &    GCR_precon_option  ==  vert_plus_xz_ADI_precon)               &
     &   Then
! need to set up coefficients for block vertical solve as well
! Copy a1_za to a1_z
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              a1_z(i,j,k) = a1_za(i,j,k)
            End Do
          End Do
        End Do

! Include horizontal second derivative terms in matrix.
! For model domain = 1 only

        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = 1, row_length
              a0_z(i,j,k) = a0_za(i,j,k)-FV_sec_theta_latitude(i,j) *   &
     &                  (HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) +            &
     &                   HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) +                &
     &                   HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) +                &
     &                   HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) )
            End Do
          End Do
        End Do

! polar boundaries
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
                              model_levels, sum_s,                      &
                              proc_row_group)

          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)            &
     &                   /global_row_length
            Do i = 1,row_length
              a0_z(i,1,k) = a0_za(i,1,k) + sum_s(k)
            End Do
          End Do
        End If

        If(at_extremity(PNorth))then
          CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0,     &
                              model_levels, sum_n,                      &
                              proc_row_group)

          Do k = 1, model_levels
            sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)         &
     &                   /global_row_length
            Do i = 1,row_length
              a0_z(i,rows,k) = a0_za(i,rows,k) + sum_n(k)
            End Do
          End Do
        End If

      End If ! on precon options

! For ADI coefficients add on relaxation pseudo_timestep
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            a0_za(i,j,k) = a0_za(i,j,k)                                 &
     &                    - recip_ADI_pseudo_timestep                   &
     &                    * rescale(i,j,k)
          End Do
        End Do
      End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

      Do j = 1, rows
        Do i = 1, row_length
          a0_za(i,j,1) = 1./a0_za(i,j,1)
        End Do
      End Do
      Do k= 2, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            factor_za(i,j,k) = a2_za(i,j,k) * a0_za(i,j,k-1)
            a0_za(i,j,k) = 1./(a0_za(i,j,k) - factor_za(i,j,k)          &
     &                      * a1_za(i,j,k-1))
          End Do
        End Do
      End Do

      If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.          &
     &    GCR_precon_option  ==  vert_plus_xz_ADI_precon)               &
     &   Then
! Repeat for block vertical pre-conditioner
! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

        Do j = 1, rows
          Do i = 1, row_length
            a0_z(i,j,1) = 1./a0_z(i,j,1)
          End Do
        End Do
        Do k= 2, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              factor_z(i,j,k) = a2_za(i,j,k) * a0_z(i,j,k-1)
              a0_z(i,j,k) = 1./(a0_z(i,j,k) - factor_z(i,j,k)           &
     &                      * a1_z(i,j,k-1))
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!     Section 4. Set-up matrix coefficients in y direction.
!-----------------------------------------------------------------------

      If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.          &
     &    GCR_precon_option  ==  xyz_ADI_precon) Then
! first do all non-polar points
        If (model_domain  ==  mt_Global) Then

          Do k = 1, model_levels
            Do j = j_start, j_end
              Do i = 1, row_length
                a1_y(i,j,k) = HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) *           &
     &                        FV_sec_theta_latitude(i,j)
                a2_y(i,j,k) = HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) *       &
     &                        FV_sec_theta_latitude(i,j)
                a0_y(i,j,k) = - a2_y(i,j,k) - a1_y(i,j,k)               &
     &                      - recip_ADI_pseudo_timestep                 &
     &                       * rescale(i,j,k)                           &
     &                      - HM_C4(i,j,k)
              End Do
            End Do

          End Do ! end loop over levels

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

          If (at_extremity(PSouth)) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                a0_y(i,2,k) = 1./a0_y(i,2,k)
              End Do
            End Do
          End If

!          Do k= 1, model_levels
!            Do j = 3, rows-1
!              Do i = 1, row_length
!                factor_y(i,j,k) = a2_y(i,j,k) * a0_y(i,j-1,k)
!                a0_y(i,j,k) = 1./(a0_y(i,j,k) - factor_y(i,j,k)
!     &                          *a1_y(i,j-1,k))
!              End Do
!            End Do
!          End Do

! DEPENDS ON: solve_forward_y1
           Call solve_forward_y1(a0_y, a1_y, a2_y, factor_y,            &
     &                         row_length, rows, model_levels)

        End If

      End If ! on precon option

!     end of routine GCR_precon_ADI_setup

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_SETUP_TRISOL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_ADI_setup_trisol
