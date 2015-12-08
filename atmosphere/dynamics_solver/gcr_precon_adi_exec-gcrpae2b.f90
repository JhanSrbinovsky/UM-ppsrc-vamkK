! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_exec_2B

      SUBROUTINE GCR_precon_ADI_exec_2B(                                &
     &                     r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,       &
     &                     HM_Czz, HM_Cz,                               &
     &                     HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,              &
     &                     HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
     &                     HM_C2n, HM_C3, HM_C4, HM_C5,                 &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     first_constant_r_rho_level,                  &
     &                     first_constant_r_rho_level_m1,               &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     FV_sec_theta_latitude,                       &
     &                     FV_cos_theta_latitude,                       &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     dlambda_p, dphi_p, dlambda_u, dphi_v,        &
     &                     recip_dlamp, recip_dphip,                    &
     &                     recip_dlamu, recip_dphiv,                    &
     &                     wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,&
     &                     model_domain, GCR_precon_option,             &
     &                     GCR_adi_add_full_soln,                       &
     &                     ADI_pseudo_timestep,                         &
     &                     GCR_n_ADI_pseudo_timesteps,                  &
     &                     weight_upper, weight_lower,                  &
     &                     a0_z, a1_z, factor_z,                        &
     &                     a0_y, a1_y, factor_y,                        &
     &                     a0_x, a1_x, factor_x,                        &
     &                     F_vector_x, G_vector_x,                      &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     n_proc, global_row_length,                   &
     &                     at_extremity, L_regular, me,                 &
     &                     n_procx, n_procy, proc_row_group,            &
     &                     j_begin, j_end,                              &
     &                     i_start, i_stop, j_start, j_stop,            &
     &                     Soln)

! Purpose:
!          Calculates ADI pre-conditioning operator applied to field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE conversions_mod, ONLY: pi
      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      
      USE precon_constants_mod, ONLY:                                   &
          no_precon,vert_precon,vert_plus_xyz_ADI_precon,xyz_ADI_precon,&
          vert_plus_xz_ADI_precon,xz_ADI_precon,Dufort_Frankel_precon
      
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)
     &, model_domain                                                    &
     &, GCR_n_ADI_pseudo_timesteps                                      &
     &, GCR_precon_option  ! 0 = no preconditioning
                           ! 1 = Vertical block pre-conditioner
                           ! 2 = 1 iteration of 1 followed by 3D ADI
                           ! 3 = 3D ADI only
                           ! 4 = 1 iteration of 1 followed by xz ADI
                                ! 5 = xz ADI only

      Logical                                                           &
     &  GCR_adi_add_full_soln ! true then use full equation on RHS
                            ! on second and subsequent ADI timesteps

      Integer                                                           &
     &  offx, offy                                                      &
     &, halo_i, halo_j                                                  &
     &, j_begin, j_end                                                  &
                                     ! loop bounds set in PE_Helmholtz
     &, i_start, i_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop                                                 &
                                     ! loop bounds set in PE_Helmholtz
     &, global_row_length                                               &
     &, me                     !  processor id

      Integer                                                           &
     &  proc_row_group                                                  &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular  ! True for regular resolution

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      Real                                                              &
     &  ADI_pseudo_timestep 

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
           ! vertical co-ordinate information
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cxy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Cyx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C2n (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C5 (1-offx:row_length+offx,1-offy:rows+offy,                 &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-offx:row_length+offx,1-offy:rows+offy,                &
     &         first_constant_r_rho_level_m1)                           &
     &, r(row_length,rows,model_levels)                                 &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)  &
     &, FV_cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
                  !  VarRes horizontal co-ordinate spacing.
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, dlambda_p(1-halo_i : row_length+halo_i)                         &
     &, dlambda_u(1-halo_i : row_length+halo_i)                         &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-offy:rows+offy)            &
     &, dphi_v(1-halo_i:row_length+halo_i, 1-offy:n_rows+offy)          &
     &, recip_dlamp(1-halo_i : row_length + halo_i)                     &
     &, recip_dlamu(1-halo_i : row_length + halo_i)                     &
     &, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j) &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      Real                                                              &
     &solno(1-offx:row_length+offx+1,1-offy:rows+offy,model_levels)

! Local Variables.

      Real                                                              &
     &  term_1                                                          &
     &, term_2                                                          &
     &, factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

      Integer                                                           &
     &  i,j,k                                                           &
     &, pseudo_timestep_number                                          &
     &, i_end                                                           &
     &, istat                                                           &
     &, ibase                                                           &
     &, len

      Real                                                              &
     &  recip_ADI_pseudo_timestep


! Local arrays.

      Real                                                              &
     &  RHS(row_length,rows,model_levels)                               &
     &, a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)                          &
     &, a0_x(row_length+1,rows,model_levels)                            &
     &, a1_x(row_length+1,rows,model_levels)                            &
     &, factor_x(row_length+1,rows,model_levels)                        &
     &, F_vector_x(row_length+1,rows,model_levels)                      &
     &, G_vector_x(2,rows,model_levels)                                 &
     &, soln_prev(row_length,rows,model_levels)                         &
     &, soln_end(rows,model_levels)                                     &
     &, L_of_soln(row_length,rows,model_levels)                         &
     &, work ( 1-offx:row_length+offx, 1-offy:rows+offy )

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External routines.
      External                                                          &
     &   GCG_rbcast

!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_EXEC_2B',zhook_in,zhook_handle)

      ibase = (me/n_procx) * n_procx

      Do pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

! Form matrix
        If (pseudo_timestep_number == 1) Then

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do  !  k = 1, model_levels

        Else  !  pseudo_timestep_number > 1
! swop soln on all procesors
! DEPENDS ON: swap_bounds
          Call swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

        If (GCR_adi_add_full_soln) Then
! DEPENDS ON: gcr_elliptic_operator_2B
          CALL GCR_Elliptic_Operator_2B(                                &
     &                           soln, row_length, rows, n_rows,        &
     &                           model_levels, model_domain,            &
     &                           first_constant_r_rho_level,            &
     &                           first_constant_r_rho_level_m1,         &
     &                           eta_theta_levels, eta_rho_levels,      &
     &                           FV_cos_theta_latitude,                 &
     &                           lambda_p, phi_p,lambda_u, phi_v,       &
     &                           dlambda_p, dphi_p,                     &
     &                           dlambda_u, dphi_v,                     &
     &                           recip_dlamp, recip_dphip,              &
     &                           recip_dlamu, recip_dphiv,              &
     &                           wt_lambda_p, wt_phi_p,                 &
     &                           wt_lambda_u, wt_phi_v,                 &
     &                           HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,    &
     &                           HM_Czz, HM_Cz,                         &
     &                           HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,        &
     &                           HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,    &
     &                           HM_C2n, HM_C3, HM_C4, HM_C5,           &
     &                           weight_upper, weight_lower,            &
     &                           L_of_soln,                             &
     &                           offx, offy, halo_i,halo_j,             &
     &                           at_extremity, proc_row_group,          &
     &                           global_row_length,                     &
     &                           i_start, i_stop, j_start, j_stop,      &
     &                           j_begin, j_end,                        &
     &                           L_regular  )

! add on to RHS and save solution from previous timestep.
            Do k = 1, model_levels
              Do j = 1, rows
                Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)      &
     &                       - L_of_soln(i,j,k)                         &
     &                      * FV_sec_theta_latitude(i,j)
                  soln_prev(i,j,k) = soln(i,j,k)
                End Do
              End Do
            End Do !  k = 1, model_levels

          Else  !  .NOT. GCR_adi_add_full_soln
! add on constant terms and save solution from previous timestep.
            Do k = 1, model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j) +  &
     &                                    HM_C4(i,j,k) * soln(i,j,k)
                  soln_prev(i,j,k) = soln(i,j,k)
                End Do
              End Do
            End Do  !  k = 1, model_levels

! Add on z, y terms
            Do k = 1, model_levels
! Include Vertical derivative term in right-hand-side
              If ( k == 1) Then
                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_2 = 1. / (eta_rho_levels(k+1) -                  &
     &                         eta_rho_levels(k))

                Do j = 1, rows
                  Do i = 1, row_length

! Calculate upper derivative
                    term_1 = ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                      * factor_2

! Calculate lower derivative. This is zero by boundary condition.

                  RHS(i,j,k) = RHS(i,j,k) -                             &
     &                              ( factor_1 * term_1 *               &
     &                               HM_Czz(i,j,k) + HM_C3(i,j,k) *     &
     &                               term_1 * HM_Cz(i,j,k) )
                  End Do
                End Do

              Else if (k == model_levels) Then

                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_3 = 1. / (eta_rho_levels(k) -                    &
     &                         eta_rho_levels(k-1))

                Do j = 1, rows
                  Do i = 1, row_length

! Calculate upper derivative. This is zero by boundary condition.

! Calculate lower derivative.
                    term_2 = ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                          * factor_3

                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                              ( factor_1 *                        &
     &                                term_2 * HM_Czz(i,j,k-1)          &
     &                               - HM_C3(i,j,k) *                   &
     &                                term_2 * HM_Cz(i,j,k-1) *         &
     &                                weight_lower(i,j,k) )

                  End Do
                End Do

              Else  !  1 < k < model_levels

                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_2 = 1. / (eta_rho_levels(k+1) -                  &
     &                         eta_rho_levels(k))
                factor_3 = 1. / (eta_rho_levels(k) -                    &
     &                         eta_rho_levels(k-1))
                Do j = 1, rows
                  Do i = 1, row_length

! Calculate upper derivative
                    term_1 = ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                        * factor_2

! Calculate lower derivative
                    term_2 = ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                        * factor_3

! Calculate operator
                    RHS(i,j,k) = RHS(i,j,k) -                           &
     &                              ( factor_1 *                        &
     &                              ( term_1 * HM_Czz(i,j,k) -          &
     &                                term_2 * HM_Czz(i,j,k-1) )        &
     &                               + HM_C3(i,j,k) *                   &
     &                              ( term_1 * HM_Cz(i,j,k) *           &
     &                                weight_upper(i,j,k) +             &
     &                                term_2 * HM_Cz(i,j,k-1) *         &
     &                                weight_lower(i,j,k) ) )

                  End Do
                End Do

              End If  ! k == 1

              if( L_regular) then
                Do j = j_begin, j_end
                  Do i = 1, row_length
                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                         ( HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k) *      &
     &                           ( soln(i,j,k) - soln(i,j+1,k) ) +      &
     &                         HM_Cyy1(i,j-1,k) * HM_Cyy2(i,j-1,k) *    &
     &                            ( soln(i,j,k)- soln(i,j-1,k) ) ) *    &
     &                          FV_sec_theta_latitude(i,j)
                  End Do
                End Do
              else  ! variable resolution
                Do j = j_begin-1, j_end
                  Do i = 1, row_length
                    work(i,j) = HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k) *       &
     &                         ( soln(i,j+1,k) - soln(i,j,k) ) /        &
     &                                         dphi_p(i,j)
                  End Do
                End Do
               Do j = j_begin, j_end
                  Do i = 1, row_length
                  RHS(i,j,k) = RHS(i,j,k) - (work(i,j) - work(i,j-1)) * &
     &                                    FV_sec_theta_latitude(i,j) /  &
     &                                                   dphi_v(i,j)
                  End Do
                End Do
              endif ! L_regular


            End Do ! k = 1, model_levels

            If(at_extremity(PSouth))then

              if( L_regular) then
                Do k = 1, model_levels
                  Do i = 1,row_length
                    sum_s_component(i,k) = HM_Cyy1(i,1,k) *             &
     &                                     HM_Cyy2(i,1,k) *             &
     &                                   ( soln(i,1,k) - soln(i,2,k) )
                  End Do
                End Do ! k = 1, model_levels
              else  ! variable resolution
                Do k = 1, model_levels
                  Do i = 1,row_length
                    sum_s_component(i,k) = HM_Cyy1(i,1,k) *             &
     &                                     HM_Cyy2(i,1,k) *             &
     &                                       dlambda_p(i) *             &
     &                                 ( soln(i,1,k) - soln(i,2,k) )
                  End Do
                End Do ! k = 1, model_levels
              endif ! L_regular

              CALL global_2d_sums(sum_s_component, row_length, 1, 0, 0, &
                                  model_levels, sum_s,                  &
                                  proc_row_group)

              if( L_regular) then
                Do k = 1, model_levels
                  sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) /    &
     &                                  global_row_length
                  Do i = 1,row_length
                    RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                  End Do
                End Do  !  k = 1, model_levels
              else  ! variable resolution
                Do k = 1, model_levels
                  sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) /    &
     &                        (2.0 * Pi * dphi_v(1,2) * dphi_v(1,2))
                  Do i = 1,row_length
                    RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                  End Do
                End Do  !  k = 1, model_levels
              endif ! L_regular
 
            End If  !  at_extremity(PSouth)

            If(at_extremity(PNorth))then

              if( L_regular) then
                Do k = 1, model_levels
                  Do i = 1,row_length
                    sum_n_component(i,k) = - HM_Cyy1(i,rows-1,k) *      &
     &                                       HM_Cyy2(i,rows-1,k) *      &
     &                              (soln(i,rows-1,k)-soln(i,rows,k))
                  End Do
                EndDo ! k = 1, model_levels
              else  ! variable resolution
                Do k = 1, model_levels
                  Do i = 1,row_length
                    sum_n_component(i,k) =  HM_Cyy1(i,rows-1,k) *       &
     &                                      HM_Cyy2(i,rows-1,k) *       &
     &                                             dlambda_p(i) *       &
     &                           ( soln(i,rows,k) - soln(i,rows-1,k) )
                  End Do
                EndDo ! k = 1, model_levels
              endif ! L_regular

              CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0, &
                                  model_levels, sum_n,                  &
                                  proc_row_group)

              if( L_regular) then
                Do k = 1, model_levels
                  sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) / &
     &                                  global_row_length
                  Do i = 1,row_length
                    RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                  End Do
                End Do  !  k = 1, model_levels
              else  ! variable resolution
                Do k = 1, model_levels
                  sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) / &
     &                  (2.0 * Pi * dphi_v(1,n_rows) * dphi_v(1,n_rows))
                  Do i = 1,row_length
                    RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                  End Do
                End Do  !  k = 1, model_levels
              endif ! L_regular
                
            End If  !  at_extremity(PNorth)

            if( L_regular) then
              Do k = 1, model_levels
                Do j = j_begin, j_end
                  Do i = 1, row_length
                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                           ( HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k) *    &
     &                             ( soln(i,j,k) - soln(i+1,j,k) ) +    &
     &                         HM_Cxx1(i-1,j,k) * HM_Cxx2(i-1,j,k) *    &
     &                            ( soln(i,j,k)- soln(i-1,j,k) ) ) *    &
     &                          FV_sec_theta_latitude(i,j)
                  End Do
                End Do
              End Do  ! k = 1, model_levels
            else  ! variable resolution
              Do k = 1, model_levels
                Do j = j_begin, j_end
                  Do i = 0, row_length
                    work(i,j) = HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k) *       &
     &                             (soln(i,j,k) - soln(i+1,j,k)) /      &
     &                                       dlambda_p(i)
                  End Do
                End Do
                Do j = j_begin, j_end
                  Do i = 1, row_length
                    RHS(i,j,k) = RHS(i,j,k) -                           &
     &                                FV_sec_theta_latitude(i,j) *      &
     &                               ( work(i,j) - work(i-1,j) ) /      &
     &                                        dlambda_u(i-1)
                  End Do
                End Do
              End Do  ! k = 1, model_levels
            endif ! L_regular

          End If  !  GCR_adi_add_full_soln

        End If ! pseudo_timestep_number == 1

!-----------------------------------------------------------------------
!     Section 2. Perform ADI sweep in lambda direction.
!                Only implemented for cyclic domains.
!-----------------------------------------------------------------------

        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i= 1, row_length
              soln(i,j,k)= RHS(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels

! solve cyclic tridiagonal system for solution
! reduce matrix to upper diagonal form.

      IF (n_procx >  1) then

        do k = 1, model_levels
          do j = 1-offy, rows+offy
            do i = 1-offx, row_length+offx
              solno(i,j,k) = soln(i,j,k)
            enddo
          enddo
        enddo  !  k = 1, model_levels

! DEPENDS ON: solve_forward_x2
        Call solve_forward_x2(solno, factor_x,                          &
     &                      row_length, rows, model_levels,             &
     &                       j_begin, j_end, offx, offy)

! back substitute to get solution.

        If (at_extremity(PEast)) Then
          Do k = 1, model_levels
            Do j = j_begin, j_end
              solno(row_length-1,j,k) = solno(row_length-1,j,k) *       &
     &                                 a0_x(row_length-1,j,k)
            End Do
          End Do
        End If  ! at_extremity(PEast)

! DEPENDS ON: solve_backward_x1
        CALL solve_backward_x1(solno, a0_x, a1_x,                       &
     &                       row_length, rows, model_levels,            &
     &                       j_begin, j_end, offx, offy)

        do k = 1, model_levels
          do j = 1-offy, rows+offy
            do i = 1-offx, row_length+offx
              soln(i,j,k) = solno(i,j,k)
            enddo
          enddo
        enddo ! k = 1, model_levels

! solve for last element of field.
! first broadcast soln(1,j,k) to all row processors from
! western most (first on row) processor
        If (at_extremity(PWest)) Then
          Do k = 1, model_levels
            Do j = j_begin, j_end
              soln_end(j,k) = soln(1,j,k)
            End Do
          End Do
        End If  !  at_extremity(PWest)

        If (n_procx  >   1) Then
          len = model_levels * rows
          Call gcg_rbcast(201, len, ibase,                              &
     &                    proc_row_group, istat, soln_end)
        End If

        i_end = row_length

        If (at_extremity(PEast)) Then

          Do k = 1, model_levels
            Do j = j_begin, j_end
              soln(row_length,j,k) = ( soln(row_length,j,k) -           &
     &                                G_vector_x(1,j,k)*soln_end(j,k)-  &
     &                                G_vector_x(2,j,k) *               &
     &                                soln(row_length-1,j,k) )          &
     &                                *a0_x(row_length,j,k)
            End Do
          End Do ! k = 1, model_levels

          i_end = row_length - 1

! solve for all other elements of field.
! first broadcast soln(row_length,j,k) to all row processors from
! eastern most (last on row) processor
          Do k = 1, model_levels
            Do j = j_begin, j_end
              soln_end(j,k) = soln(row_length,j,k)
            End Do
          End Do  !  k = 1, model_levels

        End If  ! at_extremity(PEast)

        If (n_procx > 1) Then
          len = model_levels * rows
          Call gcg_rbcast(201, len, ibase+n_procx-1,                    &
     &                    proc_row_group, istat, soln_end)
        End If !  n_procx > 1

        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i= 1, i_end
              soln(i,j,k) = soln(i,j,k) - F_vector_x(i,j,k) *           &
     &                                    soln_end(j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels

      ELSE
! solve cyclic tridiagonal system for solution
! reduce matrix to upper diagonal form.

            Do i= 2, row_length - 1
              Do k = 1, model_levels
                Do j = j_begin, j_end
             soln(i,j,k) = soln(i,j,k) - factor_x(i,j,k)*soln(i-1,j,k)
                End Do
              End Do
            End Do
! back substitute to get solution.

            Do k = 1, model_levels
              Do j = j_begin, j_end
                soln(row_length-1,j,k) = soln(row_length-1,j,k) *       &
     &                                  a0_x(row_length-1,j,k)
              End Do
            End Do
            Do i= row_length - 2, 1, -1
              Do k = 1, model_levels
                Do j = j_begin, j_end
                  soln(i,j,k) = a0_x(i,j,k) * (soln(i,j,k) -            &
     &                                    a1_x(i,j,k)*soln(i+1,j,k))
                End Do
              End Do
            End Do  !  k = 1, model_levels

! solve for last element of field.

            Do k = 1, model_levels
              Do j = j_begin, j_end
                soln(row_length,j,k) = ( soln(row_length,j,k) -         &
     &                             G_vector_x(1,j,k) * soln(1,j,k) -    &
     &                                G_vector_x(2,j,k) *               &
     &                                soln(row_length-1,j,k) )          &
     &                             *a0_x(row_length,j,k)
              End Do
            End Do  !  k = 1, model_levels

! solve for all other elements of field.

            Do k = 1, model_levels
              Do j = j_begin, j_end
                Do i= 1, row_length - 1
                  soln(i,j,k) = soln(i,j,k) - F_vector_x(i,j,k) *       &
     &                                       soln(row_length,j,k)
                End Do
              End Do
            End Do  !  k = 1, model_levels

        ENDIF   !n_procx

! swop soln on all procesors
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &                   soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

! Modify RHS to contain x term correction

        if( L_regular) then
          Do k = 1, model_levels
            Do j = j_begin, j_end
              Do i = 1, row_length
                RHS(i,j,k) = RHS(i,j,k) + (                             &
     &                       HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i+1,j,k)) +           &
     &                       HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)          &
     &                       *(soln(i,j,k)- soln(i-1,j,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do !  k = 1, model_levels
        else ! variable resolution
          Do k = 1, model_levels
            Do j = j_begin, j_end
              Do i = 0, row_length
                 work(i,j) =  HM_Cxx1(i,j,k) * HM_Cxx2(i,j,k) *         &
     &                         ( soln(i+1,j,k) - soln(i,j,k) ) /        &
     &                                      dlambda_p(i)
              End Do
            End Do
            Do j = j_begin, j_end
              Do i = 1, row_length
                RHS(i,j,k) = RHS(i,j,k) - FV_sec_theta_latitude(i,j) *  &
     &                                    (work(i,j) - work(i-1,j) ) /  &
     &                                            dlambda_u(i)
              End Do
            End Do
          End Do !  k = 1, model_levels
        endif ! L_regular

!-----------------------------------------------------------------------
!     Section 3. Perform ADI sweep in phi direction.
!-----------------------------------------------------------------------

        If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.        &
     &      GCR_precon_option  ==  xyz_ADI_precon) Then

          Do k = 1, model_levels
            Do j = j_begin, j_end
              Do i = 1, row_length
                soln(i,j,k)= RHS(i,j,k)
              End Do
            End Do
          End Do  ! k = 1, model_levels

      If ( n_procy > 1 ) then
! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

! DEPENDS ON: solve_forward_y2
          Call solve_forward_y2(soln, factor_y,                         &
     &                        row_length, rows, model_levels,           &
     &                        offx, offy)

! Back substitute to get solution.

          If(at_extremity(PNorth)) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                Soln(i,rows-1,k) = a0_y(i,rows-1,k) *                   &
     &                             soln(i,rows-1,k)
              End Do
            End Do  !  k = 1, model_levels
          End If  !  at_extremity(PNorth)

! DEPENDS ON: solve_backward_y1
          CALL solve_backward_y1(soln, a0_y, a1_y,                     &
     &                         row_length, rows, model_levels,          &
     &                         offx, offy)

! set soln at poles to zero
          If(at_extremity(PNorth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows,k) = 0.
              End Do
            End Do
          End If  !  at_extremity(PNorth)
          If(at_extremity(PSouth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,1,k) = 0.
              End Do
            End Do
          End If  !  at_extremity(PSouth)

      else ! n_procy = 1

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

            Do k= 1, model_levels
              Do j = 3, rows - 1
                Do i = 1, row_length
              soln(i,j,k) = soln(i,j,k) - factor_y(i,j,k)*soln(i,j-1,k)
                End Do
              End Do
            End Do  !  k= 1, model_levels

! Back substitute to get solution.

            Do k = 1, model_levels
              Do i = 1, row_length
                Soln(i,rows-1,k) = a0_y(i,rows-1,k) *                   &
     &                                   soln(i,rows-1,k)
              End Do
            End Do  !  k = 1, model_levels
            Do k = 1, model_levels
              Do j = rows-2, 2, -1
                Do i = 1, row_length
                  Soln(i,j,k) = a0_y(i,j,k) * ( soln(i,j,k) -           &
     &                               a1_y(i,j,k) * Soln(i,j+1,k) )
                End Do
              End Do
            End Do  !  k = 1, model_levels

          If(at_extremity(PNorth)) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                Soln(i,rows,k) = 0.
              End Do
            End Do
          End If  !  at_extremity(PNorth)

          If(at_extremity(PSouth)) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                Soln(i,1,k) = 0.
              End Do
            End Do
          End If  !  at_extremity(PSouth)

      endif ! n_procy > 1

! swap solution on all procesors
! DEPENDS ON: swap_bounds
          Call swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

! add on y direction correction
          if( L_regular ) then
            Do k = 1, model_levels
              Do j = j_begin, j_end
                Do i = 1, row_length
                  RHS(i,j,k) = RHS(i,j,k) +                             &
     &                         ( HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k) *      &
     &                           ( soln(i,j,k) - soln(i,j+1,k) ) +      &
     &                           HM_Cyy1(i,j-1,k) * HM_Cyy2(i,j-1,k) *  &
     &                              ( soln(i,j,k)- soln(i,j-1,k) ) ) *  &
     &                            FV_sec_theta_latitude(i,j)
                End Do
              End Do
            End Do !  k = 1, model_levels
          else  ! variable resolution
            Do k = 1, model_levels
              Do j = j_begin-1, j_end
                Do i = 1, row_length
                  work(i,j) =  HM_Cyy1(i,j,k) * HM_Cyy2(i,j,k) *        &
     &                         ( soln(i,j+1,k) - soln(i,j,k) ) /        &
     &                                         dphi_p(i,j)
                End Do
              End Do
              Do j = j_begin, j_end
                Do i = 1, row_length
                  RHS(i,j,k) = RHS(i,j,k) - FV_sec_theta_latitude(i,j) *&
     &                                     ( work(i,j) - work(i,j-1) ) /&
     &                                                 dphi_v(i,j)
                End Do
              End Do
            End Do !  k = 1, model_levels
          endif ! L_regular

          If(at_extremity(PSouth))then

            if( L_regular ) then
              Do k = 1, model_levels
                Do i = 1,row_length
                  sum_s_component(i,k) = HM_Cyy1(i,1,k) *               &
     &                                   HM_Cyy2(i,1,k) *               &
     &                                   ( soln(i,1,k) - soln(i,2,k) )
                End Do
              End Do  ! k = 1, model_levels
            else  ! variable resolution
              Do k = 1, model_levels
                Do i = 1,row_length
                  sum_s_component(i,k) = HM_Cyy1(i,1,k) *               &
     &                                   HM_Cyy2(i,1,k) *               &
     &                                     dlambda_p(i) *               &
     &                                   ( soln(i,1,k) - soln(i,2,k) ) /&
     &                                                 dphi_p(i,2)
                End Do
              End Do  ! k = 1, model_levels
            endif ! L_regular

            CALL global_2d_sums(sum_s_component, row_length, 1, 0, 0,   &
                                model_levels, sum_s,                    &
                                proc_row_group)

            if( L_regular ) then
              Do k = 1, model_levels
                sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1) /      &
     &                                global_row_length
                Do i = 1,row_length
                  RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                End Do
              End Do  ! k = 1, model_levels
            else  ! variable resolution
              Do k = 1, model_levels
                sum_s(k) = sum_s(k) / (Pi * dphi_v(1,2) * dphi_v(1,2))
                Do i = 1,row_length
                  RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                End Do
              End Do  ! k = 1, model_levels
            endif ! L_regular

          End If  !  at_extremity(PSouth)

          If(at_extremity(PNorth))then
            if( L_regular ) then
              Do k = 1, model_levels
                Do i = 1,row_length
                  sum_n_component(i,k) = - HM_Cyy1(i,rows-1,k) *        &
     &                                     HM_Cyy2(i,rows-1,k) *        &
     &                               ( soln(i,rows-1,k)-soln(i,rows,k) )
                End Do
              End Do ! k = 1, model_levels
            else  ! variable resolution
              Do k = 1, model_levels
                Do i = 1,row_length
                sum_n_component(i,k) =  HM_Cyy1(i,rows-1,k) *           &
     &                                  HM_Cyy2(i,rows-1,k) *           &
     &                                         dlambda_p(i) *           &
     &                           ( soln(i,rows,k) - soln(i,rows-1,k) ) /&
     &                                           dphi_p(1,rows)
                End Do
              End Do ! k = 1, model_levels
            endif ! L_regular

            CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0,   &
                                model_levels, sum_n,                    &
                                proc_row_group)

            if( L_regular ) then
              Do k = 1, model_levels
                sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows) /   &
     &                                global_row_length
                Do i = 1,row_length
                  RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                End Do
              End Do  ! k = 1, model_levels
            else  ! variable resolution
              Do k = 1, model_levels
                sum_n(k) = sum_n(k) / ( Pi * dphi_v(1,n_rows) *         &
     &                                       dphi_v(1,n_rows) )
                Do i = 1,row_length
                  RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                End Do
              End Do  ! k = 1, model_levels
            endif ! L_regular

          End If  !  at_extremity(PNorth)

        End If ! GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.
               ! GCR_precon_option  ==  xyz_ADI_precon)

!-----------------------------------------------------------------------
!     Section 4. Perform ADI sweep in vertical direction.
!-----------------------------------------------------------------------

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

        Do k = 2, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              RHS(i,j,k) = RHS(i,j,k) - factor_z(i,j,k)*RHS(i,j,k-1)
            End Do
          End Do
        End Do  !  k = 2, model_levels

! Back substitute to get solution.

        Do j = 1, rows
          Do i = 1, row_length
            Soln(i,j,model_levels) = a0_z(i,j,model_levels) *           &
     &                               RHS(i,j,model_levels)
          End Do
        End Do
        Do k = model_levels-1, 1, -1
          Do j = 1, rows
            Do i = 1, row_length
              Soln(i,j,k) = a0_z(i,j,k) * ( RHS(i,j,k) -                &
     &                            a1_z(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do  !  k = model_levels-1, 1, -1

! add on correction to previous solution
        If ( pseudo_timestep_number > 1 ) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                Soln(i,j,k) = Soln(i,j,k) + Soln_prev(i,j,k)
              End Do
            End Do
          End Do  !  k = 1, model_levels
        End If  !  pseudo_timestep_number > 1

      End Do ! pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

!     end of routine GCR_precon_ADI_exec_2B

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_EXEC_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_ADI_exec_2B
