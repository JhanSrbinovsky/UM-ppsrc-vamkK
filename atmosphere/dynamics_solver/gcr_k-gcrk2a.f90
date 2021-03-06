! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_k

      Subroutine GCR_k(                                                 &
                       rows, row_length, model_levels,                  &
                       model_domain, timestep_number, CycleNo,          &
                       GCR_Diagnostics, GCR_max_its, GCR_min_its,       &
                       GCR_sum_its, GCR_its_switch, GCR_its_avg_step,   &
                       GCR_max_time, GCR_min_time,                      &
                       GCR_Max_iterations, GCR_Error_Tolerance,         &
                       GCR_Restart_value, GCR_Abs_Tolerance,            &
                       GCR_use_Abs_Tol, GCR_zero_init_guess,            &
                       GCR_use_residual_Tol,                            &
                       GCR_precon_option, GCR_ADI_Pseudo_timestep,      &
                       GCR_n_ADI_pseudo_timesteps,                      &
                       GCR_adi_add_full_soln, L_gcr_fast_x,             &
                       first_constant_r_rho_level,                      &
                       first_constant_r_rho_level_m1,                   &
                       delta_lambda, delta_phi,                         &
                       eta_theta_levels, eta_rho_levels,                &
                       r_theta_levels, r_rho_levels,                    &
                       FV_sec_theta_latitude, cos_v_latitude,           &
                       FV_cos_theta_latitude,                           &
                       HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,              &
                       HM_Czz, HM_Cz, HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,   &
                       HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,              &
                       HM_C2, HM_C3, HM_C4, HM_C5,                      &
                       HM_RHS, rescale,                                 &
                       weight_upper, weight_lower, soln,                &
                       offx, offy, me, n_rows, n_proc,                  &
                       n_procx, n_procy,                                &
                       at_extremity, neighbour, l_datastart,            &
                       proc_row_group, g_row_len, proc_col_group,       &
                       number_dif_horiz_points, halo_i, halo_j,         &
                       global_rows,                                     &
                       g_rows, g_row_length,                            &
                       ldump                                            &
                       )

! Purpose:
!          Solves Helmholtz equation using GCR(k) Method with
!          preconditioning of solution.
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


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE precon_constants_mod, ONLY:                                   &
          no_precon,vert_precon,vert_plus_xyz_ADI_precon,xyz_ADI_precon,&
          vert_plus_xz_ADI_precon,xz_ADI_precon,Dufort_Frankel_precon
      
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
        row_length                                                      &
                         ! number of point on a row.
      , rows                                                            &
                         ! number of rows.
      , model_levels     ! number of model levels.

      Integer                                                           &
        offx,offy,n_rows,n_proc,halo_i,halo_j

      Integer                                                           &
        model_domain                                                    &
                         ! holds integer code for model domain
                         ! 1 = global
                         ! 2 = limited area
                         ! 3 = periodic in x limited area
      , first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
      , first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)
      , timestep_number                                                 &
      , CycleNo

      Integer                                                           &
        GCR_Restart_value                                               &
                           ! After how many iterations do we restart
      , GCR_Max_iterations                                              &
                           ! Maximum number of iterations to do
      , GCR_Diagnostics                                                 &
                           ! Switch controlling diagnostic output.
                           ! 0 = none
                           ! 1 = initial and final residuals
                           ! 2 = all
                           ! 3 = iteration count processing
      , GCR_its_switch                                                  &
                            ! Iterations analysis switch
      , GCR_max_its                                                     &
                            ! Max iterations this period
      , GCR_min_its                                                     &
                            ! Min iterations this period
      , GCR_max_time                                                    &
                            ! Timestep number for max GCR its
      , GCR_min_time                                                    &
                            ! Timestep number for min GCR its
      , GCR_its_avg_step(3)                                             &
                            ! Iterations analysis step now
      , GCR_sum_its                                                     &
                           ! Sum iterations over test period

      , GCR_precon_option                                               &
                           ! 0 = no preconditioning
                           ! 1 = block diagonal vertical pre-conditioner
                           ! 2 = 1 iteration of 1 followed by 3D ADI
                           ! 3 = 3D ADI only
                           ! 4 = 1 iteration of 1 followed by xz ADI
                           ! 5 = xz ADI only
      , GCR_n_ADI_pseudo_timesteps   ! Number of ADI pseudo timesteps
                                     ! to perform.

      Real                                                              &
        GCR_Error_Tolerance                                             &
                            ! Acceptable error tolerance relative to
                            ! initial Error
      , GCR_Abs_Tolerance                                               &
                          ! Acceptable tolerance on absolute error
      , GCR_ADI_pseudo_timestep

      ! True if this is a dumping period. If true saved variables, are
      ! reset so they compare with a run that is restarted from the dump
      Logical, Intent(In) :: ldump

      Logical                                                           &
        GCR_use_abs_Tol                                                 &
                            ! True if absolute tolerance to be used
      , GCR_use_residual_Tol                                            &
                             ! True if residual tolerance to be used
      , GCR_zero_init_guess                                             &
                            ! True if initial guess to solution is
                            ! zero.
      , GCR_adi_add_full_soln                                           &
                              ! true then use full equation on RHS
                            ! on second and subsequent ADI timesteps
      , L_gcr_fast_x        ! true then user faster non reproducible
                            !      code

      Real                                                              &
        HM_Cxx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cxx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cxy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cxy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cyy1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cyy2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cyx1 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Cyx2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)  &
      , HM_Czz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
      , HM_Cz (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
      , HM_C2 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
      , HM_C3 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
      , HM_C4 (1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
      , HM_C5 (1-offx:row_length+offx,1-offy:rows+offy,                 &
               first_constant_r_rho_level_m1)                           &
      , HM_Cxz (1-offx:row_length+offx,1-offy:rows+offy,                &
               first_constant_r_rho_level_m1)                           &
      , HM_Cyz (1-offx:row_length+offx,1-offy:rows+offy,                &
               first_constant_r_rho_level_m1)                           &
      , HM_Cxp (1-offx:row_length+offx,1-offy:rows+offy,                &
               first_constant_r_rho_level_m1)                           &
      , HM_Cyp (1-offx:row_length+offx,1-offy:rows+offy,                &
               first_constant_r_rho_level_m1)

      Real                                                              &
        rescale (1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
        weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
                 model_levels)                                          &
      , weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
                 model_levels)

      Real                                                              &
        HM_RHS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                  !right-hand-side

      Real                                                              &
            ! Trigonometric functions
        FV_sec_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy) &
                              ! finite volume secant array
      , FV_cos_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy) &
                              ! finite volume cosine array
      , cos_v_latitude (1-offx:row_length+offx,1-offy:n_rows+offy)

      Real                                                              &
           ! vertical co-ordinate information
        r_theta_levels (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,&
                        0:model_levels)                                 &
      , r_rho_levels (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,  &
                      model_levels)                                     &
      , eta_theta_levels(0:model_levels)                                &
      , eta_rho_levels(model_levels)

      Real                                                              &
           ! horizontal grid-spacings
        delta_lambda                                                    &
      , delta_phi

! parallel variables   IN
      Integer                                                           &
       l_datastart(3)                                                   &
      ,  n_procx, n_procy,number_dif_horiz_points                       &
      , neighbour(4), proc_row_group,g_row_len,me,proc_col_group        &
      , global_rows

      Integer                                                           &
        g_rows(0:n_proc-1)                                              &
      , g_row_length(0:n_proc-1)


      Logical                                                           &
        at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid



! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
        soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                     !solution, an initial guess
                                           ! must be supplied.

! Arguments with Intent OUT. ie: variables Output only

! Local Variables.

      Integer                                                           &
        i, j, k, kl                                                     &
                     ! Loop indices
      , n_iterations                                                    &
                     ! number of iterations completed
      , n_restart                                                       &
                     ! number of steps completed before restarting
      , i_start                                                         &
      , i_stop                                                          &
      , j_start                                                         &
      , j_stop

! ----------------------- include file: GCR_ITER_DIM -------------------
! Description: Introduce integer parameters to:
!              (1) declare the maximum allowed number of dynamics 
!                  cycles (timestep iterations)
!              (2) declare the number of GCR Iterations analysis steps
!
!
      INTEGER, PARAMETER :: MAX_NUMCYCLES = 10 ! Max number of 
                                               ! dynamics cycles
      INTEGER, PARAMETER :: GCR_ANAL_STEPS = 3 ! Number of GCR Iterations
                                               ! analysis step

! To speed up the code, GCR_iterations is the first iteration for which
! a convergence test is performed.
! GCR_iterations_min is the minimum number of iterations over the model
! run
! N.B. to get bit-comparison between NRUNs and CRUNs we need to reset
! GCR_iterations and GCR_iterations_min every time we write out a dump.

      Integer, Parameter :: GCR_iterations_initial_val = 1
      Integer, Parameter :: GCR_iterations_min_initial_val = 1000

      Integer, Save :: GCR_iterations(max_numcycles) = &
         GCR_iterations_initial_val
                                ! Start with 1 iteration
                                ! up to 10 timestepping iterations are
                                ! accounted which is the allowed
                                ! maximum. 
                         

      Integer, Save :: GCR_iterations_min(max_numcycles) = &
         GCR_iterations_min_initial_val  
                                ! Previous min

       integer nrp1

      Logical                                                           &
        L_converged

      Real                                                              &
        Beta                                                            &
                     ! coefficient of GCR(k)
      , Error_Norm                                                      &
                     ! Norm of the Error
      , Abs_Norm                                                        &
                       ! Norm of the Increment added to the
                             ! current solution.
      , Initial_Error_Norm                                              &
                             ! Norm of the Initial Error
      , avg_its                                                         &
                  ! iteration average over test period
      , avg_its_denom  ! Averaging period in iterations analysis

! Local arrays

      Real                                                              &
        p(1-offx:row_length+offx,1-offy:rows+offy,model_levels,         &
          0:GCR_Restart_value)                                          &
                      ! preconditioned residual
      , r(row_length,rows,model_levels)                                 &
                      ! residual
      , Alpha(0:GCR_Restart_value-1)                                    &
      , L_of_p(row_length,rows,model_levels,                            &
               0:GCR_Restart_value)                                     &
      , a0_y(row_length,rows,model_levels)                              &
      , a1_y(row_length,rows,model_levels)                              &
      , factor_y(row_length,rows,model_levels)                          &
      , a0_za(row_length,rows,model_levels)                             &
      , a1_za(row_length,rows,model_levels)                             &
      , factor_za(row_length,rows,model_levels)                         &
      , a0_z(row_length,rows,model_levels)                              &
      , a1_z(row_length,rows,model_levels)                              &
      , factor_z(row_length,rows,model_levels)                          &
      , a0_x(row_length+1,rows,model_levels)                            &
      , a1_x(row_length+1,rows,model_levels)                            &
      , factor_x(row_length+1,rows,model_levels)                        &
      , F_vector_x(row_length+1,rows,model_levels)                      &
      , G_vector_x(2,rows,model_levels)                                 &
      , HM_C2n (1-offx:row_length+offx,1-offy:rows+offy,                &
               first_constant_r_rho_level_m1)                           &
      , init_error_mean(global_rows)

      Real                                                              &
        bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
      , bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
      , bv_factor_forward(2,2*n_procx,rows,model_levels)                &
      , bv_factor_backward(2*n_procx,rows,model_levels)                 &
      , recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
      , bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
      , recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
        bv_soln_n_term(rows, model_levels)                              &
      , bv_soln_1_term1(rows, model_levels)                             &
      , bv_soln_1_term2(rows, model_levels)

      real,allocatable,dimension(:,:,:)       :: factor_forward
      real,allocatable,dimension(:,:,:)       :: factor_backward
      integer                                 :: ilen
      integer                                 :: j_start_c,j_end_c

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1.   Set Initial values.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_K',zhook_in,zhook_handle)
      If (.not. GCR_use_residual_Tol .and. .not. GCR_use_abs_tol) Then
        GCR_use_residual_Tol = .true.
      End If

      If (model_domain  ==  mt_global .or.                              &
          model_domain  ==  mt_bi_cyclic_LAM) Then
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start=1
        j_stop = rows
      Else If(model_domain  ==  mt_lam)then
! Solve over interior points only.
        if(at_extremity(PEast))then
          i_stop = row_length  - 1
        Else
          i_stop = row_length
        End If
        if(at_extremity(PWest))then
          i_start = 2
        Else
          i_start = 1
        End If
        if(at_extremity(PSouth))then
          j_start = 2
        Else
          j_start = 1
        End If
        if(at_extremity(PNorth))then
          j_stop = rows - 1
        Else
          j_stop = rows
        End If
      Elseif (model_domain  ==  mt_cyclic_LAM) then
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        If (at_extremity(PSouth)) Then
          j_start = 2
        else
          j_start = 1
        End If
        If (at_extremity(PNorth)) Then
          j_stop = rows - 1
        else
          j_stop = rows
        End If
      End If

      If(model_domain == mt_global) then
        j_start_c = 1
        j_end_c   = rows
        If (at_extremity(PSouth) ) Then
          j_start_c = 2
        End If
        If (at_extremity(PNorth) ) Then
          j_end_c   = rows-1
        End If
      Else ! dimension factor_ arrays to 1, if not mt_global
        j_start_c = 1
        j_end_c = 1
        ilen = 1
      End if

!kk   Array can only be allocated after j_startc and j_end_c have
!     been computed. Exact dimension is necessary to allow loop
!     collapsing in Mpp_tri_solve_exec

      ilen = row_length+1
      allocate(factor_forward(ilen,j_start_c:j_end_c,model_levels))
      allocate(factor_backward(ilen,j_start_c:j_end_c,model_levels))

! DEPENDS ON: gcr_two_norm
      Call GCR_Two_Norm(                                                &
                        HM_RHS, row_length, rows,                       &
                        model_levels, model_domain, Error_Norm,         &
                        offx, offy, at_extremity, n_proc,               &
                        proc_col_group, proc_row_group,                 &
                        number_dif_horiz_points, l_datastart            &
                        )

!     Check RHS non-zero
      If (Error_Norm  >   tiny(1.0)) then   ! if RHS is non-zero

! set HM_C2n = HM_C2 / delta r
! this saves repeating the exercise ever time gcr_elliptic_op is
! called.

      Do k = 1, first_constant_r_rho_level - 1
        Do j= 0, rows+1
          Do i= 0, row_length+1
            HM_C2n(i,j,k) = HM_C2(i,j,k)                                &
                          / (r_rho_levels(i,j,k+1) -                    &
                             r_rho_levels(i,j,k))
          End Do
        End Do
      End Do

      If ( .not. GCR_zero_init_guess) Then
! Elliptic operator applied to initial guess at solution.
! store answer in r.

! DEPENDS ON: gcr_elliptic_operator
        Call GCR_Elliptic_Operator(                                     &
                                 soln, row_length,                      &
                                 rows, model_levels, model_domain,      &
                                 first_constant_r_rho_level,            &
                                 first_constant_r_rho_level_m1,         &
                                 eta_theta_levels, eta_rho_levels,      &
                                 r_theta_levels, r_rho_levels,          &
                                 FV_cos_theta_latitude,                 &
                                 HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,    &
                                 HM_Czz, HM_Cz,                         &
                                 HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,        &
                                 HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,    &
                                 HM_C2n, HM_C3, HM_C4, HM_C5,           &
                                 weight_upper, weight_lower, r,         &
                                 offx, offy, at_extremity, n_rows,      &
                                 g_row_len, n_proc,                     &
                                 proc_row_group, halo_i,halo_j          &
                                 )

! Initial Residual is elliptic operator applied to initial guess at
! solution minus right-hand-side

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              r(i,j,k) = r(i,j,k) - HM_RHS(i,j,k) *                     &
                                    FV_cos_theta_latitude(i,j)
            End Do
          End Do
        End Do

      Else

! Initial Residual is minus right-hand-side

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              r(i,j,k) = - HM_RHS(i,j,k) *                              &
                           FV_cos_theta_latitude(i,j)
            End Do
          End Do
        End Do

      End If

! set-up preconditioning array
      If (GCR_precon_option  ==  vert_precon .or.                       &
          model_domain  ==  mt_lam) Then
! DEPENDS ON: gcr_precon_1_setup
        Call GCR_precon_1_setup(                                        &
                        HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,             &
                        HM_Czz, HM_Cz, HM_C3, HM_C4,                    &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        row_length, rows, model_levels,                 &
                        FV_sec_theta_latitude,                          &
                        model_domain,                                   &
                        weight_upper, weight_lower,                     &
                        a0_z, a1_z, factor_z,                           &
                        offx, offy, n_proc, g_row_len,                  &
                        at_extremity,                                   &
                        n_procx, n_procy, proc_row_group,               &
                        halo_i,halo_j                                   &
                          )
      Else If ( (GCR_precon_option  /=  no_precon) .and.                &
        (GCR_precon_option  /=  vert_precon) ) Then
! Note: also sets up block vertical pre-conditioner if
! precon_option eq vert_plus_xyz_ADI_precon or vert_plus_xz_ADI_precon
        If (L_gcr_fast_x) Then

! DEPENDS ON: gcr_precon_adi_setup_trisol
          Call GCR_precon_adi_setup_trisol(j_start_c,j_end_c,           &
                        HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,             &
                        HM_Czz, HM_Cz, HM_C3, HM_C4,                    &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        row_length, rows, model_levels,                 &
                        FV_sec_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_ADI_pseudo_timestep,                        &
                        rescale,                                        &
                        weight_upper, weight_lower,                     &
                        a0_z, a1_z, factor_z,                           &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &

! a_plus_x, a_minus_x replaced by a1_x, a0_x to save memeory
                        a1_x, a0_x,                                     &
                        factor_forward,factor_backward,                 &
                        recip_a_central_x,                              &
                        recip_bv_a_matrix_diag, bv_a_matrix_sup,        &
                        bv_a_matrix_0, bv_a_matrix_np1,                 &
                        bv_factor_forward, bv_factor_backward,          &
                        bv_soln_n_term,                                 &
                        bv_soln_1_term1, bv_soln_1_term2,               &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, me,                               &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group, halo_i, halo_j)
        Else
! DEPENDS ON: gcr_precon_adi_setup
          Call GCR_precon_adi_setup(                                    &
                        HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,             &
                        HM_Czz, HM_Cz, HM_C3, HM_C4,                    &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        row_length, rows, model_levels,                 &
                        FV_sec_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_ADI_pseudo_timestep,                        &
                        rescale,                                        &
                        weight_upper, weight_lower,                     &
                        a0_z, a1_z, factor_z,                           &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &
                        a0_x, a1_x, factor_x,                           &
                        F_vector_x, G_vector_x,                         &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, me,                               &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group, halo_i, halo_j)
        End If
      End If

! Call Pre-conditioner to calculate p from r.

      If (GCR_precon_option  ==  no_precon) Then
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              p(i,j,k,0) = r(i,j,k)
            End Do
          End Do
        End Do

      Else If (GCR_precon_option  ==  vert_precon .or.                  &
               GCR_precon_option  ==  vert_plus_xyz_ADI_precon .or.     &
               GCR_precon_option  ==  vert_plus_xz_ADI_precon           &
               .or. model_domain  ==  mt_lam) Then

! DEPENDS ON: gcr_precon_1_exec
        Call GCR_precon_1_exec(                                         &
                        r, model_domain,                                &
                        row_length, rows, model_levels,                 &
                        FV_sec_theta_latitude,                          &
                        a0_z, a1_z, factor_z, p(1-offx,1-offy,1,0)      &
!!      parallel variables added
                       ,offx,offy, at_extremity                         &
                        )

      Else
        If (L_gcr_fast_x) Then
! DEPENDS ON: gcr_precon_adi_exec_trisol
        Call GCR_precon_adi_exec_trisol(j_start_c,j_end_c,              &
                        r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                        HM_Czz, HM_Cz,                                  &
                        HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,                 &
                        HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,             &
                        HM_C2n, HM_C3, HM_C4, HM_C5,                    &
                        row_length, rows, n_rows, model_levels,         &
                        first_constant_r_rho_level,                     &
                        first_constant_r_rho_level_m1,                  &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        FV_sec_theta_latitude,                          &
                        FV_cos_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_adi_add_full_soln,                          &
                        GCR_ADI_pseudo_timestep,                        &
                        GCR_n_ADI_pseudo_timesteps, rescale,            &
                        weight_upper, weight_lower,                     &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &
                        a1_x, a0_x,                                     &
                        factor_forward,factor_backward,                 &
                        recip_a_central_x,                              &
                        bv_a_matrix_0, bv_a_matrix_np1,                 &
                        recip_bv_a_matrix_diag, bv_a_matrix_sup,        &
                        bv_factor_forward, bv_factor_backward,          &
                        bv_soln_n_term,                                 &
                        bv_soln_1_term1, bv_soln_1_term2,               &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, neighbour, me,                    &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group,halo_i,halo_j,                   &
                        p(1-offx,1-offy,1,0))

        Else
! DEPENDS ON: gcr_precon_adi_exec
          Call GCR_precon_adi_exec(                                     &
                        r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                        HM_Czz, HM_Cz,                                  &
                        HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,                 &
                        HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,             &
                        HM_C2n, HM_C3, HM_C4, HM_C5,                    &
                        row_length, rows, n_rows, model_levels,         &
                        first_constant_r_rho_level,                     &
                        first_constant_r_rho_level_m1,                  &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        FV_sec_theta_latitude,                          &
                        FV_cos_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_adi_add_full_soln,                          &
                        GCR_ADI_pseudo_timestep,                        &
                        GCR_n_ADI_pseudo_timesteps, rescale,            &
                        weight_upper, weight_lower,                     &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &
                        a0_x, a1_x, factor_x,                           &
                        F_vector_x, G_vector_x,                         &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, neighbour, me,                    &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group,halo_i,halo_j,                   &
                        p(1-offx,1-offy,1,0))
        End If

      End If

! DEPENDS ON: swap_bounds
      Call swap_bounds(                                                 &
             p, row_length, rows, model_levels, offx, offy,             &
             fld_type_p, .false.)

! Elliptic operator applied to p

! DEPENDS ON: gcr_elliptic_operator
      Call GCR_Elliptic_Operator(                                       &
                                 p(1-offx,1-offy,1,0), row_length,      &
                                 rows, model_levels, model_domain,      &
                                 first_constant_r_rho_level,            &
                                 first_constant_r_rho_level_m1,         &
                                 eta_theta_levels, eta_rho_levels,      &
                                 r_theta_levels, r_rho_levels,          &
                                 FV_cos_theta_latitude,                 &
                                 HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,    &
                                 HM_Czz, HM_Cz,                         &
                                 HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,        &
                                 HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,    &
                                 HM_C2n, HM_C3, HM_C4, HM_C5,           &
                                 weight_upper, weight_lower,            &
                                 L_of_p(1,1,1,0),                       &
                                 offx, offy, at_extremity, n_rows,      &
                                 g_row_len, n_proc,                     &
                                 proc_row_group, halo_i, halo_j         &
                                 )

! Calculate initial error norm.

      If (GCR_use_residual_Tol ) Then
! DEPENDS ON: gcr_two_norm
        Call GCR_Two_Norm(                                              &
                        r, row_length, rows, model_levels,              &
                        model_domain, Initial_Error_Norm,               &
                        0, 0, at_extremity, n_proc,                     &
                        proc_col_group, proc_row_group,                 &
                        number_dif_horiz_points, l_datastart )
      Else
        Initial_Error_Norm = 0.0
        Error_Norm = 0.0
      End If

      If (GCR_Diagnostics == 2) Then
! DEPENDS ON: gcr_error_print
        Call GCR_error_print(                                           &
                             r, proc_row_group,                         &
                             g_row_len, global_rows,                    &
                             g_rows, g_row_length,                      &
                             n_proc, n_procx, n_procy, me,              &
                             row_length, rows, model_levels,            &
                             model_domain, at_extremity,                &
                             init_error_mean,0)
      End If  !  GCR_Diagnostics == 2

      If (GCR_use_abs_Tol ) Then
! Calculate Abs Norm.
! first calculate term to norm, store in p(,,n_restart+1) as workspace

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              p(i,j,k,1) = abs( r(i,j,k) * FV_sec_theta_latitude(i,j)   &
                                / rescale(i,j,k) )
            End Do
          End Do
        End Do

! DEPENDS ON: gcr_calc_abs_norm
        Call GCR_calc_abs_norm(                                         &
                              p(1-offx,1-offy,1,1), offx, offy,         &
                              n_proc, model_domain,                     &
                              at_extremity,                             &
                              row_length, rows,                         &
                              model_levels, Abs_Norm)

      If (GCR_Diagnostics == 2) Then
! DEPENDS ON: gcr_abs_print
          Call GCR_abs_print(                                           &
                           p(1-offx,1-offy,1,1), offx, offy,            &
                           proc_row_group,                              &
                           g_row_len, global_rows,                      &
                           g_rows, g_row_length,                        &
                           n_proc, n_procx, n_procy, me,                &
                           row_length, rows,                            &
                           model_levels, model_domain,                  &
                           at_extremity)
      End If  !  GCR_Diagnostics == 2

      Else
! Initialise to zero
        Abs_Norm = 0.
      End If

! If diagnostics requested print out certain information

      If (GCR_Diagnostics == 2) Then
        If(me == 0)then
          write (6,*) ' '
          write (6,*) ' =============================================='
          write (6,*) ' Elliptic Solver GCR(',GCR_Restart_value,') '    &
                     ,' requested. '
          write (6,*) ' maximum number of iterations is : ',            &
                        GCR_Max_iterations
          write (6,*) ' maximum number of solution vectors before '
          write (6,*) ' restarting is : ',GCR_Restart_value
          write (6,*) ' Effective number of iteration is thus : ',      &
                        GCR_Max_iterations * GCR_Restart_value
          If (GCR_use_residual_Tol ) Then
            write (6,*) ' Residual Convergence criteria is : ',         &
                      GCR_Error_Tolerance,' x Initial Error Norm. '
          End If
          If (GCR_use_abs_Tol ) Then
            write (6,*) ' Absolute error norm convergence less than ',  &
                        GCR_Abs_Tolerance
          End If
        End If

      End If  !  GCR_Diagnostics == 2

      If ( GCR_Diagnostics ==1 .or. GCR_Diagnostics ==2 ) Then
        If (me == 0) then
          write (6,*) ' =============================================='
          If (GCR_use_residual_Tol ) Then
            write (6,*) ' Initial Error Norm is : ', Initial_Error_Norm
          End If
          If (GCR_use_abs_Tol ) Then
            write (6,*) ' initial Absolute Norm : ',Abs_Norm
          End If
        End If
      End If ! GCR_Diagnostics ==1 .or. GCR_Diagnostics ==2

! ----------------------------------------------------------------------
! Section 2.   Solve via iterative method
! ----------------------------------------------------------------------

      n_iterations = 0
      L_converged  = .false.

      Do        !     while n_iterations < GCR_Max_iterations

        n_iterations = n_iterations + 1   !  update iteration counter

        If (GCR_Diagnostics == 2) Then
         if(me == 0)then
          write (6,*) ' '
          write (6,*) ' Iteration ',n_iterations
          write (6,*) ' ============================================='
          write (6,*) ' '
         endif
        End If  !  GCR_Diagnostics == 2

        n_restart = 0

! Loop over the desired k solution space vectors, k is held in
! GCR_Restart_value
        Do     !    while n_restart < GCR_Restart_value

! Calculate beta coefficient

! DEPENDS ON: gcr_coefficient
          Call GCR_Coefficient(                                         &
                     r, L_of_p(1,1,1,n_restart), row_length,            &
                     rows, model_levels, model_domain, Beta,            &
                     at_extremity, n_proc,                              &
                     proc_col_group, proc_row_group,                    &
                     l_datastart )

! Add on Beta * p to current solution and
! add on Beta * L_of_p to current residual to get residual at next
! iteration level.

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(j_start, j_stop,i_start, i_stop,model_levels,             &
!$OMP&        soln,beta,p,r,L_of_p,n_iterations,GCR_iterations,         &
!$OMP&        CycleNo,GCR_use_abs_Tol,rows,row_length,rescale,          &
!$OMP&        FV_sec_theta_latitude,n_restart)
          Do k = 1, model_levels
            Do j = j_start, j_stop
!dir$ unroll 8
              Do i = i_start, i_stop
                soln(i,j,k) = soln(i,j,k) + Beta * p(i,j,k,n_restart)
                r(i,j,k) = r(i,j,k) + Beta * L_of_p(i,j,k,n_restart)
              End Do
            End Do

!  Test for convergence only if n_iterations > GCR_iterations(CycleNo)  
! Calculate Abs Norm.
! first calculate term to norm, store in p(,,n_restart+1) as workspace

            If ( n_iterations > GCR_iterations(CycleNo)               &
                              .and. GCR_use_abs_Tol) Then
              Do j = 1, rows
                Do i = 1, row_length
                  p(i,j,k,n_restart+1) = abs ( r(i,j,k)                 &
                                       * FV_sec_theta_latitude(i,j)     &
                                       / rescale(i,j,k) )
                End Do
              End Do
            End If
          End Do
!$OMP END PARALLEL DO

          If ( n_iterations > GCR_iterations(CycleNo) ) Then

            If (GCR_use_abs_Tol ) Then

! DEPENDS ON: gcr_calc_abs_norm
              Call GCR_calc_abs_norm(                                   &
                              p(1-offx,1-offy,1,n_restart+1),           &
                              offx, offy, n_proc, model_domain,         &
                              at_extremity,                             &
                              row_length, rows,                         &
                              model_levels, Abs_Norm)

! Check for convergence
              If ( Abs_Norm < GCR_Abs_Tolerance ) L_converged = .true.

            ElseIf (GCR_use_residual_Tol ) Then
! Calculate Error Norm

! DEPENDS ON: gcr_two_norm
              Call GCR_Two_Norm(                                        &
                            r, row_length, rows,                        &
                            model_levels, model_domain, Error_Norm,     &
                            0, 0, at_extremity, n_proc,                 &
                            proc_col_group, proc_row_group,             &
                            number_dif_horiz_points, l_datastart )

! Check for convergence
              If ( Error_Norm <                                         &
                          Initial_Error_Norm * GCR_Error_Tolerance )    &
                 L_converged = .true.

            End If !  GCR_use_abs_Tol

            If (GCR_Diagnostics == 2) Then

! Output value of Error_Norm as diagnostic if required.
! Output Increment Norm if required.
              If(me == 0)then

                If (GCR_use_residual_Tol .and. GCR_use_abs_Tol) Then

            write (6,*) 'Step ',n_restart+1,' Error Norm : ',Error_Norm &
                        ,' Absolute Norm : ',Abs_Norm
                Else If (GCR_use_residual_Tol) Then
            write (6,*) 'Step ',n_restart+1,' Error Norm : ',Error_Norm
                Else If (GCR_use_abs_Tol) Then
            write (6,*) 'Step ',n_restart+1                             &
                        ,' Absolute Norm : ',Abs_Norm
                End If
              End If
            End If  !  GCR_Diagnostics == 2

          EndIf ! n_iterations > GCR_iterations(CycleNo)

          if ( L_converged ) EXIT

! Calculate q by applying pre-conditioner to current value of r.
! store in next available p location.

           If (GCR_precon_option  ==  no_precon) Then
              Do k = 1, model_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    p(i,j,k,n_restart+1) = r(i,j,k)
                  End Do
                End Do
              End Do

            Else If (GCR_precon_option  ==  vert_precon                 &
                     .or. model_domain  ==  mt_lam) Then

! DEPENDS ON: gcr_precon_1_exec
              Call GCR_precon_1_exec(                                   &
                            r, model_domain,                            &
                            row_length, rows, model_levels,             &
                            FV_sec_theta_latitude,                      &
                            a0_z, a1_z, factor_z,                       &
                            p(1-offx,1-offy,1,n_restart+1),             &
                            offx, offy, at_extremity )

            Else
              If (L_gcr_fast_x) Then
! DEPENDS ON: gcr_precon_adi_exec_trisol
                Call GCR_precon_adi_exec_trisol(j_start_c,j_end_c,      &
                        r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                        HM_Czz, HM_Cz,                                  &
                        HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,                 &
                        HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,             &
                        HM_C2n, HM_C3, HM_C4, HM_C5,                    &
                        row_length, rows, n_rows, model_levels,         &
                        first_constant_r_rho_level,                     &
                        first_constant_r_rho_level_m1,                  &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        FV_sec_theta_latitude,                          &
                        FV_cos_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_adi_add_full_soln,                          &
                        GCR_ADI_pseudo_timestep,                        &
                        GCR_n_ADI_pseudo_timesteps, rescale,            &
                        weight_upper, weight_lower,                     &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &
                        a1_x, a0_x,                                     &
                        factor_forward,factor_backward,                 &
                        recip_a_central_x,                              &
                        bv_a_matrix_0, bv_a_matrix_np1,                 &
                        recip_bv_a_matrix_diag, bv_a_matrix_sup,        &
                        bv_factor_forward, bv_factor_backward,          &
                        bv_soln_n_term,                                 &
                        bv_soln_1_term1, bv_soln_1_term2,               &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, neighbour, me,                    &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group,halo_i,halo_j,                   &
                        p(1-offx,1-offy,1,n_restart+1))
              Else
! DEPENDS ON: gcr_precon_adi_exec
                Call GCR_precon_adi_exec(                               &
                        r, HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
                        HM_Czz, HM_Cz,                                  &
                        HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,                 &
                        HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,             &
                        HM_C2n, HM_C3, HM_C4, HM_C5,                    &
                        row_length, rows, n_rows, model_levels,         &
                        first_constant_r_rho_level,                     &
                        first_constant_r_rho_level_m1,                  &
                        eta_theta_levels, eta_rho_levels,               &
                        r_theta_levels, r_rho_levels,                   &
                        FV_sec_theta_latitude,                          &
                        FV_cos_theta_latitude,                          &
                        model_domain, GCR_precon_option,                &
                        GCR_adi_add_full_soln,                          &
                        GCR_ADI_pseudo_timestep,                        &
                        GCR_n_ADI_pseudo_timesteps, rescale,            &
                        weight_upper, weight_lower,                     &
                        a0_za, a1_za, factor_za,                        &
                        a0_y, a1_y, factor_y,                           &
                        a0_x, a1_x, factor_x,                           &
                        F_vector_x, G_vector_x,                         &
                        offx,offy,n_proc, g_row_len,                    &
                        at_extremity, neighbour, me,                    &
                        n_procx,n_procy,proc_row_group,                 &
                        proc_col_group,halo_i,halo_j,                   &
                        p(1-offx,1-offy,1,n_restart+1))
              End If

            End If

! DEPENDS ON: swap_bounds
            Call swap_bounds(                                           &
                               p(1-offx,1-offy,1,n_restart+1),          &
                               row_length,rows,model_levels,            &
                               offx, offy, fld_type_p, .false.)

! Calculate operator applied to this q value, and store in next
! available L_of_p location.

! DEPENDS ON: gcr_elliptic_operator
            Call GCR_Elliptic_Operator(                                 &
                            p(1-offx,1-offy,1,n_restart+1), row_length, &
                                      rows, model_levels, model_domain, &
                                      first_constant_r_rho_level,       &
                                      first_constant_r_rho_level_m1,    &
                                      eta_theta_levels, eta_rho_levels, &
                                      r_theta_levels, r_rho_levels,     &
                                      FV_cos_theta_latitude,            &
                                      HM_Cxx1, HM_Cxx2, HM_Cyy1,        &
                                      HM_Cyy2, HM_Czz, HM_Cz,           &
                                      HM_Cxz, HM_Cyz, HM_Cxp, HM_Cyp,   &
                                      HM_Cxy1, HM_Cxy2,                 &
                                      HM_Cyx1, HM_Cyx2,                 &
                                      HM_C2n, HM_C3, HM_C4, HM_C5,      &
                                      weight_upper, weight_lower,       &
                                      L_of_p(1,1,1,n_restart+1),        &
                                      offx, offy, at_extremity, n_rows, &
                                      g_row_len, n_proc,                &
                                      proc_row_group,                   &
                                      halo_i, halo_j )

! Calculate all Alpha values for this new L_of_p

            Do kl = 0, n_restart

! DEPENDS ON: gcr_coefficient
              Call GCR_Coefficient(                                     &
                               L_of_p(1,1,1,n_restart+1),               &
                               L_of_p(1,1,1,kl), row_length,            &
                               rows, model_levels, model_domain,        &
                               Alpha(kl),                               &
                               at_extremity, n_proc,                    &
                               proc_col_group, proc_row_group,          &
                               l_datastart)

            End Do

! Calculate value of p at n_restart + 1, and
! Calculate value of L_of_p at n_restart + 1

            nrp1=n_restart+1

            If(n_restart == 1) Then

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(row_length,rows,model_levels,nrp1,p,alpha,l_of_p)
              Do k = 1, model_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    p(i,j,k,nrp1) = p(i,j,k,nrp1) +                     &
                                           Alpha(0) * p(i,j,k,0) +      &
                                           Alpha(1) * p(i,j,k,1)
                    L_of_p(i,j,k,nrp1) =L_of_p(i,j,k,nrp1) +            &
                                  Alpha(0) * L_of_p(i,j,k,0) +          &
                                  Alpha(1) * L_of_p(i,j,k,1)
                  End Do
                End Do
              End Do !k=1,model_levels
!$OMP END PARALLEL DO
            Else
!!!! putting the openmp before the kl loop slows the code down 0.20 secs
!!!!  becomes 0.30 secs and the openmp factor goes from 1.88 to 1.16
            Do kl = 0, n_restart

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP&  SHARED(p,alpha,nrp1,row_length,rows,model_levels,L_of_p,kl)
              Do k = 1, model_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    p(i,j,k,nrp1) = p(i,j,k,nrp1) +                     &
                                           Alpha(kl) * p(i,j,k,kl)
                    L_of_p(i,j,k,nrp1) =                                &
                                  L_of_p(i,j,k,nrp1) +                  &
                                  Alpha(kl) * L_of_p(i,j,k,kl)
                  End Do
                End Do
              End Do  !  k = 1, model_levels
!$OMP END PARALLEL DO
            End Do  !  kl = 0, n_restart

            End If   ! n_restart == 1

            n_restart = n_restart + 1   ! as last step update counter
            if( n_restart == GCR_Restart_value ) EXIT

        End Do ! n_restart < GCR_Restart_value

! After end of each Restart section reset values at end index back to
! zero, unless solution has converged or last iteration reached.

        If ( L_converged .or. n_iterations == GCR_Max_iterations) EXIT

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i)         &
!$OMP& SHARED(p,L_of_p,GCR_Restart_value,row_length,rows,model_levels)
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                p(i,j,k,0) = p(i,j,k,GCR_Restart_value)
                L_of_p(i,j,k,0) = L_of_p(i,j,k,GCR_Restart_value)
              End Do
            End Do
          End Do
!$OMP END PARALLEL DO

      End Do    !    while n_iterations < GCR_Max_iterations

!  Set minimum number of iterations for next step
!  ~80% of the minimum number of iterations should be fairly safe
!  ~40% is used for cycles (timestepping iterations) after the first
!       as these tend to vary more

      If ( CycleNo == 1 ) Then 
        GCR_iterations_min(1) = MIN(GCR_iterations_min(1),              &
                                                     n_iterations)
        GCR_iterations(1) = (4*GCR_iterations_min(1))/5
      Else
        GCR_iterations_min(CycleNo) = MIN(GCR_iterations_min(CycleNo),  &
                                                     n_iterations)
        GCR_iterations(CycleNo) = (2*GCR_iterations_min(CycleNo))/5
      End If

! DEPENDS ON: icwriteentry
      call ICwriteEntry(n_iterations)

! Output converged or not message

      If (GCR_Diagnostics > 0) Then
        If (GCR_Diagnostics < 3) Then
          If ( me == 0 ) Then
            If (L_converged) Then
              write (6,*) ' GCR(',GCR_Restart_value,') converged in '   &
                                        ,n_iterations,' iterations. '
          If (GCR_Diagnostics  ==  1 .and. GCR_use_residual_Tol) Then
            write (6,*) ' Final Error Norm : ',Error_Norm
          End If
          If (GCR_Diagnostics  ==  1 .and. GCR_use_abs_Tol ) Then
            write (6,*) ' Final Absolute Norm : ',Abs_Norm
          End If
          write (6,*) ' =============================================='
            Else  !  L_converged is false
          write (6,*) ' GCR(',GCR_Restart_value,') failed to converge', &
                     ' in ',n_iterations,' iterations. '
          If (GCR_Diagnostics  ==  1 .and. GCR_use_residual_Tol) Then
            write (6,*) ' Final Error Norm : ',Error_Norm
          End If
          If (GCR_Diagnostics  ==  1 .and. GCR_use_abs_Tol ) Then
            write (6,*) ' Final Absolute Norm : ',Abs_Norm
          End If
          write (6,*) ' =============================================='
            End If  !  L_converged
          End If  !  me == 0

        Else     ! GCR_Diagnostics = 3
          GCR_sum_its = GCR_sum_its + n_iterations
          if( GCR_max_its < n_iterations)then
              GCR_max_its = n_iterations
              GCR_max_time = timestep_number
          endif !  GCR_max_its < n_iterations
          if( GCR_min_its > n_iterations)then
             GCR_min_its = n_iterations
             GCR_min_time = timestep_number
          endif !  GCR_min_its > n_iterations
          if( n_iterations == GCR_Max_iterations ) then
              write (6,996) n_iterations, timestep_number
              write (6,*)                                               &
                    '  ****   WARNING This run is likely to FAIL *****'
          elseif( timestep_number  ==  GCR_its_avg_step(GCR_its_switch))&
             then
            if( GCR_its_switch ==1 ) then
              avg_its_denom = float( GCR_its_avg_step(GCR_its_switch) )
            else
              avg_its_denom = float( GCR_its_avg_step(GCR_its_switch) - &
                                 GCR_its_avg_step(GCR_its_switch - 1) )
            endif ! GCR_its_switch ==1
            avg_its = float( GCR_sum_its ) / avg_its_denom
            if( GCR_its_switch == 1 ) then
              write (6,997) GCR_its_avg_step(1), avg_its
              write (6,999) GCR_max_its, GCR_max_time,                  &
                         GCR_min_its, GCR_min_time
            else ! GCR_its_switch > 1
              write (6,998) GCR_its_avg_step(GCR_its_switch - 1),       &
                         GCR_its_avg_step(GCR_its_switch), avg_its
              write (6,999) GCR_max_its, GCR_max_time,                  &
                         GCR_min_its, GCR_min_time
            endif ! GCR_its_switch == 1
            GCR_its_switch = GCR_its_switch + 1
            if( GCR_its_switch > 3 ) GCR_its_switch = 3
            GCR_sum_its = 0   !  Re-initialise GCR_sum_its
            GCR_max_its = 0   !  Re-initialise GCR_max_its
            GCR_min_its = 1000   !  Re-initialise GCR_min_its
            GCR_max_time = 0  !  Re-initialise GCR_max_time
            GCR_min_time = 0  !  Re-initialise GCR_min_time
          elseif( mod(timestep_number,                                  &
                      GCR_its_avg_step(GCR_its_switch))== 0) then
            avg_its_denom = float ( GCR_its_avg_step(GCR_its_switch) )
            avg_its = float( GCR_sum_its ) / avg_its_denom
            write (6,998)                                               &
                  timestep_number-GCR_its_avg_step(GCR_its_switch)      &
                 ,timestep_number, avg_its
            write (6,999) GCR_max_its, GCR_max_time,                    &
                       GCR_min_its, GCR_min_time
            GCR_sum_its = 0    !  Re-initialise GCR_sum_its
            GCR_max_its = 0    !  Re-initialise GCR_max_its
            GCR_min_its = 1000    !  Re-initialise GCR_min_its
            GCR_max_time = 0   !  Re-initialise GCR_max_time
            GCR_min_time = 0   !  Re-initialise GCR_min_time
          End If ! timestep_number == GCR_its_avg_step(GCR_its_switch)
        End If   ! GCR_Diagnostics < 3
      End If   !  GCR_Diagnostics > 0

      If (GCR_Diagnostics == 2) Then
! DEPENDS ON: gcr_error_print
        Call GCR_error_print(                                           &
                           r, proc_row_group,                           &
                           g_row_len, global_rows,                      &
                           g_rows, g_row_length,                        &
                           n_proc, n_procx, n_procy, me,                &
                           row_length, rows, model_levels,              &
                           model_domain, at_extremity,                  &
                           init_error_mean,1)

        If (GCR_use_abs_Tol ) Then
! Calculate Abs Norm.
! first calculate term to norm, store in p(,,n_restart+1) as workspace

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                p(i,j,k,1) = abs( r(i,j,k) * FV_sec_theta_latitude(i,j) &
                                / rescale(i,j,k) )
              End Do
            End Do
          End Do

! DEPENDS ON: gcr_abs_print
          Call GCR_abs_print(                                           &
                           p(1-offx,1-offy,1,1), offx, offy,            &
                           proc_row_group,                              &
                           g_row_len, global_rows,                      &
                           g_rows, g_row_length,                        &
                           n_proc, n_procx, n_procy, me,                &
                           row_length, rows,                            &
                           model_levels, model_domain,                  &
                           at_extremity)

        End If ! GCR_use_abs_Tol
      End If ! GCR_Diagnostics == 2

      else                      ! if RHS is zero
        soln(:,:,:) = 0.0       ! set soln to zero
        L_converged = .true.
        If (GCR_Diagnostics  >   0) Then
          If (me == 0) Then
            write (6,*)                                                 &
            ' RHS zero so GCR(',GCR_Restart_value,') not needed '
          End If
        End If
        write (6,*) ' =============================================='
      End If
      deallocate(factor_forward)
      deallocate(factor_backward)


 996  FORMAT('  WARNING iterations = ',I5,' at timestep ',I6)
 997  FORMAT(' Between start and timestep ',I6,                         &
             ' average iterations = ', F10.3)
 998  FORMAT(' Between timestep ', I6,' and ', I6,                      &
             ' average iterations =  ', F10.3)
 999  FORMAT(' Iterations: Max =',I4,' at timestep ', I6,               &
             '. Min =',I4,' at timestep ', I6)

      if (ldump) then
        ! Ensure these values are the same in the next time step whether
        ! or not the run continues or is restarted from the dump that is
        ! about to be written.
        gcr_iterations(cycleNo) = gcr_iterations_initial_val
        GCR_iterations_min(cycleNo) = gcr_iterations_min_initial_val
      end if

      IF (lhook) CALL dr_hook('GCR_K',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE GCR_k
