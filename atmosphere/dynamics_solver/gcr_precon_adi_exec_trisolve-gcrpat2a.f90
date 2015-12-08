! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_ADI_exec_trisolve

      Subroutine GCR_precon_ADI_exec_trisol(j_start,j_end,              &
     &                     r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,       &
     &                     HM_Czz, HM_Cz,                               &
     &                     HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,              &
     &                     HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,          &
     &                     HM_C2n, HM_C3, HM_C4, HM_C5,                 &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     first_constant_r_rho_level,                  &
     &                     first_constant_r_rho_level_m1,               &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     FV_sec_theta_latitude,                       &
     &                     FV_cos_theta_latitude,                       &
     &                     model_domain, GCR_precon_option,             &
     &                     GCR_adi_add_full_soln,                       &
     &                     ADI_pseudo_timestep,                         &
     &                     GCR_n_ADI_pseudo_timesteps, rescale,         &
     &                     weight_upper, weight_lower,                  &
     &                     a0_z, a1_z, factor_z,                        &
     &                     a0_y, a1_y, factor_y,                        &

! new
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1, bv_soln_1_term2,            &

     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity, neighbour, me,                 &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     proc_col_group,halo_i,halo_j,                &
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


      USE global_2d_sums_mod, ONLY: global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      
      USE precon_constants_mod, ONLY:                                   &
          no_precon,vert_precon,vert_plus_xyz_ADI_precon,xyz_ADI_precon,&
          vert_plus_xz_ADI_precon,xz_ADI_precon,Dufort_Frankel_precon
      
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      integer, intent(in)                   :: j_start,j_end

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
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j, me                                              &
     &, neighbour(4)

      Integer                                                           &
     &  proc_row_group                                                  &
     &, proc_col_group

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

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
     &  r_theta_levels (1-halo_i:row_length+halo_i                      &
     &                 ,1-halo_j:rows+halo_j,0:model_levels)            &
     &, r_rho_levels (1-halo_i:row_length+halo_i                        &
     &                 ,1-halo_j:rows+halo_j,model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
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
     &, rescale (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)  &
     &, FV_cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
     &  a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)

      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &
     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Local Variables.

      Real                                                              &
     &  term_1                                                          &
     &, term_2                                                          &
     &, factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

      Integer                                                           &
     &  i,j,k                                                           &
     &, pseudo_timestep_number

      Integer                 :: istat

      Integer, Parameter :: jblock=4

      integer :: len1
      integer :: len2
      integer :: jj  
      integer :: j1  
      integer :: j2  
      integer :: my_j_start
      integer :: my_j_stop 

      Real                                                              &
     &  recip_ADI_pseudo_timestep

! Local arrays.

!kk   1. dimension row_length+1 to avoid copy in Mpp_tri_solve_exec
      real,dimension(row_length+1,rows,model_levels)           :: RHS
!kk   To save copying, a secon array with j_start:j_end is required
      real,dimension(row_length+1,j_start:j_end,model_levels)  :: RHS_c

      Real                                                              &
     &  a0_z(row_length,rows,model_levels)                              &
     &, a1_z(row_length,rows,model_levels)                              &
     &, factor_z(row_length,rows,model_levels)                          &
     &, a0_y(row_length,rows,model_levels)                              &
     &, a1_y(row_length,rows,model_levels)                              &
     &, factor_y(row_length,rows,model_levels)                          &
     &, soln_prev(row_length,rows,model_levels)                         &
     &, L_of_soln(row_length,rows,model_levels)

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
!     Section 1. Set initial guess to solution and initialise
!                variables.
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_EXEC_TRISOL',zhook_in,zhook_handle)
      recip_ADI_pseudo_timestep = 1./ ADI_pseudo_timestep

      Do pseudo_timestep_number = 1, GCR_n_ADI_pseudo_timesteps

! Form matrix
        If (pseudo_timestep_number  ==  1) Then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)         &
!$OMP& SHARED(model_levels,j_start,j_end ,row_length,r,RHS_c,           &
!$OMP&        FV_sec_theta_latitude,RHS,rows)
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = 1, row_length
                RHS_c(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
!           If needed, set boundary rows directly in RHS
            Do j = 1,j_start-1
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
            Do j = j_end+1, rows
              Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do
!$OMP END PARALLEL DO
        Else
! swop soln on all procesors
! DEPENDS ON: swap_bounds
        Call swap_bounds(soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

        If (GCR_adi_add_full_soln) Then
! DEPENDS ON: gcr_elliptic_operator
          Call GCR_Elliptic_Operator(                                   &
     &                           soln, row_length,                      &
     &                           rows, model_levels, model_domain,      &
     &                           first_constant_r_rho_level,            &
     &                           first_constant_r_rho_level_m1,         &
     &                           eta_theta_levels, eta_rho_levels,      &
     &                           r_theta_levels, r_rho_levels,          &
     &                           FV_cos_theta_latitude,                 &
     &                           HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,    &
     &                           HM_Czz, HM_Cz,                         &
     &                           HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,        &
     &                           HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,    &
     &                           HM_C2n, HM_C3, HM_C4, HM_C5,           &
     &                           weight_upper, weight_lower,            &
     &                           L_of_soln,                             &
     &                           offx, offy, at_extremity, n_rows,      &
     &                           global_row_length, n_proc,             &
     &                           proc_row_group, halo_i,halo_j          &
     &                           )

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
            End Do

          Else
! add on constant terms and save solution from previous timestep.
            Do k = 1, model_levels
              Do j = 1, rows
                Do i = 1, row_length
                RHS(i,j,k) = r(i,j,k) * FV_sec_theta_latitude(i,j)      &
     &                       + HM_C4(i,j,k) * soln(i,j,k)
                soln_prev(i,j,k) = soln(i,j,k)
                End Do
              End Do
            End Do

! Add on z, y terms
            Do k = 1, model_levels
! Include Vertical derivative term in right-hand-side
              If ( k  ==  1) Then
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
     &                              ( factor_1 *                        &
     &                               term_1 * HM_Czz(i,j,k)             &
     &                              + HM_C3(i,j,k) *                    &
     &                               term_1 * HM_Cz(i,j,k) )
                  End Do
                End Do

              Else if (k  ==  model_levels) Then
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

              Else
                factor_1 = 1. / (eta_theta_levels(k) -                  &
     &                         eta_theta_levels(k-1))
                factor_2 = 1. / (eta_rho_levels(k+1) -                  &
     &                         eta_rho_levels(k))
                factor_3 = 1. / (eta_rho_levels(k) -                    &
     &                         eta_rho_levels(k-1))
                Do j = 1, rows
!dir$ split
!dir$ unroll 4
                  Do i = 1, row_length

                    RHS(i,j,k) = RHS(i,j,k) -                           &
     &                       ( soln(i,j,k+1) - soln(i,j,k) )            &
     &                        * factor_2 *                              &
     &                              ( factor_1 *                        &
     &                                 HM_Czz(i,j,k)                    &
     &                               + HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k) *                   &
     &                                weight_upper(i,j,k) )
                  Enddo
                Enddo
                Do j = 1, rows
!dir$ split
!dir$ unroll 4
                  Do i = 1, row_length
                    RHS(i,j,k) = RHS(i,j,k) +                           &
     &                       ( soln(i,j,k) - soln(i,j,k-1) )            &
     &                        * factor_3 *                              &
     &                              ( factor_1 *                        &
     &                                  HM_Czz(i,j,k-1)                 &
     &                               - HM_C3(i,j,k) *                   &
     &                                 HM_Cz(i,j,k-1) *                 &
     &                                weight_lower(i,j,k)  )

                  End Do
                End Do

              End If

              Do j = j_start, j_end
!dir$ split
!dir$ unroll 4
                Do i = 1, row_length
                  RHS(i,j,k) = RHS(i,j,k) + (                           &
     &                       HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i,j+1,k)) +           &
     &                       HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k)          &
     &                       *(soln(i,j,k)- soln(i,j-1,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
                End Do
              End Do

            End Do ! end loop over levels

            If(at_extremity(PSouth))then
              Do k = 1, model_levels
                Do i = 1,row_length
                sum_s_component(i,k)= HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)     &
     &                                *(soln(i,1,k)-soln(i,2,k))
                End Do
              End Do
            End If
            If(at_extremity(PNorth))then
              Do k = 1, model_levels
                Do i = 1,row_length
                sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*            &
     &                                  HM_Cyy2(i,rows-1,k)             &
     &                          *(soln(i,rows-1,k)-soln(i,rows,k))
                End Do
              End Do
            End If

            If(at_extremity(PSouth))then
              CALL global_2d_sums(sum_s_component, row_length, 1, 0, 0, &
                                  model_levels, sum_s,                  &
                                  proc_row_group)
              Do k = 1, model_levels
                sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)        &
     &                     /global_row_length
                Do i = 1,row_length
                  RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
                End Do
              End Do
            End If

            If(at_extremity(PNorth))then
              CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0, &
                                  model_levels, sum_n,                  &
                                  proc_row_group)
              Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)       &
     &                     /global_row_length
                Do i = 1,row_length
                  RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
                End Do
              End Do
            End If


           Do k = 1, model_levels
             Do j = j_start, j_end
               Do i = 1, row_length
                 RHS_c(i,j,k) = RHS(i,j,k) + (                          &
     &                     HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)                &
     &                     *(soln(i,j,k) - soln(i+1,j,k)) +             &
     &                     HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)            &
     &                     *(soln(i,j,k)- soln(i-1,j,k)))               &
     &                     *FV_sec_theta_latitude(i,j)
               End Do
             End Do
           End Do

          End If

        End If ! on first timestep or not

!-----------------------------------------------------------------------
!     Section 2. Perform ADI sweep in lambda direction.
!                Only implemented for cyclic domains.
!-----------------------------------------------------------------------

! copy rhs into L_of_Soln and use this to pass into subroutine.
! required since RHS passed into subroutine is changed on output and
! this routine requires an unchanged RHS later on.
!         Do k = 1, model_levels
!           Do j = j_start, j_end
!             Do i= 1, row_length
!               L_of_soln(i,j,k)= RHS(i,j,k)
!             End Do
!           End Do
!         End Do

! DEPENDS ON: mpp_tri_solve_exec
          Call mpp_tri_solve_exec(                                      &
     &                     row_length, rows, model_levels,              &
     &                     offx, offy, j_start, j_end,                  &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2,                             &
     &                     L_of_soln, RHS_c, soln)

! swop soln on all procesors
! DEPENDS ON: swap_bounds
        Call swap_bounds(                                               &
     &                   soln,row_length,rows,model_levels,             &
     &                   offx,offy,fld_type_p,.FALSE.)

! Modify RHS to contain x term correction

!$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(NONE) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,j_start,j_end ,row_length,RHS_c,             &
!$OMP& HM_Cxx1,RHS,HM_Cxx2,soln,FV_sec_theta_latitude)
       Do k = 1, model_levels
         Do j = j_start, j_end
           Do i = 1, row_length
             RHS(i,j,k) = RHS_c(i,j,k) + (                              &
     &                       HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k)              &
     &                       *(soln(i,j,k) - soln(i+1,j,k)) +           &
     &                       HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k)          &
     &                       *(soln(i,j,k)- soln(i-1,j,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
           End Do
         End Do
       End Do
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!     Section 3. Perform ADI sweep in phi direction.
!-----------------------------------------------------------------------

      If (GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.          &
     &    GCR_precon_option  ==  xyz_ADI_precon )                       &
     &    Then
          Do k = 1, model_levels
            Do i= 1, row_length
              Do j = j_start, j_end
                soln(i,j,k)= RHS(i,j,k)
              End Do
            End Do
          End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

!          Do k= 1, model_levels
!            Do j = 3, rows - 1
!              Do i = 1, row_length
!                soln(i,j,k)= soln(i,j,k)-
!     &                       factor_y(i,j,k)*soln(i,j-1,k)
!              End Do
!            End Do
!          End Do
! DEPENDS ON: solve_forward_y2
          Call solve_forward_y2(soln, factor_y,                         &
     &                        row_length, rows, model_levels,           &
     &                        offx, offy)

! Back substitute to get solution.

          If(at_extremity(PNorth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows-1,k) = a0_y(i,rows-1,k) *                   &
     &                             soln(i,rows-1,k)
              End Do
            End Do
          End If

!          Do k= 1, model_levels
!            Do j = rows-2, 2, -1
!              Do i = 1, row_length
!                Soln(i,j,k) = a0_y(i,j,k) * ( soln(i,j,k) -
!     &                              a1_y(i,j,k) * Soln(i,j+1,k) )
!              End Do
!            End Do
!          End Do

! DEPENDS ON: solve_backward_y1
          Call solve_backward_y1(soln, a0_y, a1_y,                      &
     &                         row_length, rows, model_levels,          &
     &                         offx, offy)

! set soln at poles to zero
          If(at_extremity(PNorth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,rows,k) = 0.
              End Do
            End Do
          End If
          If(at_extremity(PSouth)) Then
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,1,k) = 0.
              End Do
            End Do
          End If

! swop soln on all procesors
! DEPENDS ON: swap_bounds
          CALL swap_bounds(                                             &
     &                     soln,row_length,rows,model_levels,           &
     &                     offx,offy,fld_type_p,.FALSE.)

! add on y direction correction
          Do k = 1, model_levels
!CDIR unroll=4
            Do j = j_start, j_end
              Do i = 1, row_length
                RHS(i,j,k) = RHS(i,j,k) + (                             &
     &                     HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k)                &
     &                       *(soln(i,j,k) - soln(i,j+1,k)) +           &
     &                       HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k)          &
     &                       *(soln(i,j,k)- soln(i,j-1,k)))             &
     &                       *FV_sec_theta_latitude(i,j)
              End Do
            End Do
          End Do

          If(at_extremity(PSouth))then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_s_component(i,k)= HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)     &
     &                                *(soln(i,1,k)-soln(i,2,k))
              End Do
            End Do
          End If
          If(at_extremity(PNorth))then
            Do k = 1, model_levels
              Do i = 1,row_length
                sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*            &
     &                                  HM_Cyy2(i,rows-1,k)             &
     &                           *(soln(i,rows-1,k)-soln(i,rows,k))
              End Do
            End Do
          End If

          If(at_extremity(PSouth))then
            CALL global_2d_sums(sum_s_component, row_length, 1, 0, 0,   &
                                model_levels, sum_s,                    &
                                proc_row_group)

            Do k = 1, model_levels
              sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)          &
     &                     /global_row_length
              Do i = 1,row_length
                RHS(i,1,k) = RHS(i,1,k) + sum_s(k)
              End Do
            End Do
          End If

          If(at_extremity(PNorth))then
            CALL global_2d_sums(sum_n_component, row_length, 1, 0, 0,   &
                                model_levels, sum_n,                    &
                                proc_row_group)

            Do k = 1, model_levels
              sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)       &
     &                     /global_row_length
              Do i = 1,row_length
                RHS(i,rows,k) = RHS(i,rows,k) + sum_n(k)
              End Do
            End Do
          End If

        End If ! on precon option
!-----------------------------------------------------------------------
!     Section 4. Perform ADI sweep in vertical direction.
!-----------------------------------------------------------------------

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

! This code has been cache blocked by 4.
! Some OpenMP added.
      my_j_start=1
      my_j_stop=rows
!
      len1=jblock
      len2=(my_j_stop-my_j_start+1)/len1+1
!
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j,k,i,jj,j1,j2)&
!$OMP&  SHARED(my_j_start,my_j_stop,len1,len2,model_levels,row_length,  &
!$OMP&   factor_z,RHS,Soln,a0_z,a1_z)
      Do jj=1, len1
        j1=my_j_start+(jj-1)*len2
        j2=min(my_j_start+jj*len2-1, my_j_stop)
        Do k= 2, model_levels
          Do j = j1, j2
            Do i = 1, row_length
              RHS(i,j,k) = RHS(i,j,k) - factor_z(i,j,k)*RHS(i,j,k-1)
            End Do
          End Do
        End Do

! Back substitute to get solution.

        Do j = j1, j2
          Do i = 1, row_length
            Soln(i,j,model_levels) = a0_z(i,j,model_levels) *           &
     &                               RHS(i,j,model_levels)
          End Do
        End Do

        Do k= model_levels-1, 1, -1
          Do j = j1, j2
            Do i = 1, row_length
              Soln(i,j,k) = a0_z(i,j,k) * ( RHS(i,j,k) -                &
     &                            a1_z(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do
      End Do  ! jj
!$OMP END PARALLEL DO

! add on correction to previous solution
        If (pseudo_timestep_number  >   1) Then
          Do j = 1, rows
!CDIR unroll=8
            Do k= 1, model_levels
              Do i = 1, row_length
                Soln(i,j,k) = Soln(i,j,k) + Soln_prev(i,j,k)
              End Do
            End Do
          End Do
        End If

      End Do ! end loop over number of pseudo timesteps

!     end of routine GCR_precon_ADI_exec

      IF (lhook) CALL dr_hook('GCR_PRECON_ADI_EXEC_TRISOL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_ADI_exec_trisol
