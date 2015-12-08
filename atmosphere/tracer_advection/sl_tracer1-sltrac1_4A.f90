! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Tracer1
!
MODULE SL_Tracer1_mod
IMPLICIT NONE
CONTAINS

      SUBROUTINE SL_Tracer1_4A(                                         &
                            super_array_size,                           &
                            eta_theta_levels,                           &
                            r_rho_levels, r_theta_levels,               &
                            p_star, p, p_theta_levels,                  &
                            rho, l_tracer1_non_hydro,                   &
                            row_length, rows, n_rows, model_levels,     &
                            halo_lam, halo_phi,                         &
                            FV_cos_theta_latitude,                      &
                            me, n_proc, n_procx, n_procy,               &
                            halo_i, halo_j, l_datastart,                &
                            g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,      &
                            group_2dcomm, at_extremity,                 &
                            global_row_length, global_rows,             &
                            proc_row_group, proc_col_group,             &
                            proc_all_group, off_x, off_y,               &
                            high_order_scheme_moist,                    &
                            monotone_scheme_moist,                      &
                            model_domain, L_high_moist,                 &
                            L_mono_moist,                               &
                            check_bottom_levels,                        &
                            interp_vertical_search_tol,                 &
                            first_constant_r_rho_level,                 &
                            depart_lambda, depart_phi,                  &
                            depart_r,                                   &
                            CO2, L_CO2_interactive,                     &
                            Murk, L_murk_advect,                        &
                            dust_div1,dust_div2,dust_div3,              &
                            dust_div4,dust_div5,dust_div6,              &
                            L_dust,                                     &
                            soot_new, soot_agd, soot_cld, L_soot,       &
                            bmass_new, bmass_agd, bmass_cld,            &
                            l_biomass,                                  &
                            ocff_new, ocff_agd, ocff_cld, l_ocff,       &
                            so2, so4_aitken, so4_accu,                  &
                            so4_diss, nh3, dms,                         &
                            L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,      &
                            nitr_acc, nitr_diss, L_nitrate,             &
                            tracers, tr_levels, tr_vars,                &
                            tr_lbc_vars,tracer_lbcs,                    &
                            A_max_trvars,A_tr_active_lbc_index,         &
                            tracers_ukca, tr_ukca,                      &
                            L_use_cariolle, ozone_tracer,               &
                            tr_lbc_ukca,ukca_tracer_lbcs,               &
                            a_max_ukcavars,                             &
                            ukca_tr_active_lbc_index,                   &
                            rimwidth,rimweights,                        &
                            lenrim,lbc_size,lbc_start,                  &
                            Error_Code)

! Purpose:
!          Performs semi-Lagrangian advection of tracers
!          (based on sl_thermo)
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE dust_parameters_mod, ONLY: l_twobin_dust
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE atm_fields_bounds_mod

      USE eg_interpolation_eta_mod
      USE Field_Types
      USE proc_info_mod,      ONLY: gc_proc_row_group, gc_proc_col_group
      USE horiz_grid_mod,     ONLY: xi1_p, xi2_p, xi1_p, xi2_p
      USE dynamics_input_mod, ONLY: l_sl_bc_correction

      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER                                                           &
        row_length                                                      &
                     ! number of points on a row
      , rows                                                            &
                     ! number of rows.
      , n_rows                                                          &
                     ! number of v-rows.
      , model_levels                                                    &
                     ! Number of model levels.
      , me                                                              &
                     ! My processor number
      , n_proc                                                          &
                     ! Total number of processors
      , n_procx                                                         &
                     ! Number of processors in longitude
      , n_procy                                                         &
                     ! Number of processors in latitude
      , halo_i                                                          &
                     ! Size of large halo in i.
      , halo_j                                                          &
                     ! Size of large halo in j.
      , off_x                                                           &
                     ! Size of small halo in i
      , off_y                                                           &
                     ! Size of small halo in j.
      , l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
      , proc_row_group                                                  &
                       ! Group id for processors on the same row
      , proc_col_group                                                  &
                       ! Group id for processors on the same column
      , proc_all_group                                                  &
                       ! Group id for all processors
      , max_look                                                        &
                             ! max size of look-up arrays for searches
      , global_row_length                                               &
                            ! global number of points on a row
      , global_rows                                                     &
                            ! global number of rows
      , g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                             ! processor on my processor-row
                             ! holding a given value in i direction
      , g_j_pe(1-halo_j:global_rows      +halo_j)                       &
                             ! processor on my processor-column
                             ! holding a given value in j direction
      , rimwidth                                                        &
                       ! Width of boundaries in LBCs     
      , lenrim                                                          &
                       ! Size of single level of LBC data
      , lbc_size(4)                                                     &
                       ! Size of each side of LBC data   
      , lbc_start(4)                                                    
                       ! Start of each side in LBC data 

      REAL   rimweights(rimwidth) 
                       ! Weights to apply to the LBCs    


      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand
      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: size_int_mult ! error check size for comms on demand

      INTEGER                                                           &
        tr_levels                                                       &
      , tr_vars                                                         &
      , tr_ukca              ! No of UKCA tracers

      INTEGER, INTENT(IN) :: tr_lbc_vars    ! No. of active tracer lbcs
      INTEGER, INTENT(IN) :: tr_lbc_ukca
                                       ! No. of active UKCA tracer lbcs
      INTEGER, INTENT(IN) :: A_max_trvars
                                       ! Max number of tracer vars
      INTEGER, INTENT(IN) :: A_max_ukcavars
                                       ! Max number of UKCA tracer vars
      INTEGER, INTENT(IN) ::                                            &
        A_tr_active_lbc_index(A_max_trvars) 
                                       ! Active tracer lbcs index
      INTEGER, INTENT(IN) ::                                            &
        UKCA_tr_active_lbc_index(A_max_ukcavars) 
                                       ! Active UKCA tracer lbcs index

      INTEGER                                                           &
        first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      LOGICAL                                                           &
       l_tracer1_non_hydro

      INTEGER                                                           &
        high_order_scheme_moist                                         &
                                 ! a code saying which high order
                           ! scheme to use for moist variables.
      , monotone_scheme_moist                                           &
                              ! a code saying which monotone
                           ! scheme to use for moist variables.
      , interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
      , check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      INTEGER                                                           &
        model_domain     ! holds integer code for model domain

      LOGICAL                                                           &
        L_high_moist                                                    &
                       ! True, if high order interpolation required
                       !       for moist variables.
      , L_mono_moist
                       ! True, if interpolation required to be monotone
                       !       for moist variables.

      REAL, INTENT(IN) ::                                               &
        tracer_lbcs(lenrim,tr_levels,tr_lbc_vars)      !Tracer lbcs
      REAL, INTENT(IN) ::                                               &
        ukca_tracer_lbcs(lenrim,tr_levels,tr_lbc_ukca) !UKCA Tracer lbcs

      LOGICAL, INTENT(IN) ::                                            &
        L_CO2_interactive                                               &
      , L_murk_advect                                                   &
      , L_DUST                                                          &
      , L_Soot, L_biomass, L_ocff, L_nitrate                            &
      , L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms                           &
      , L_USE_CARIOLLE

      REAL, INTENT(INOUT) :: co2                                      &
                        (tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: murk                                    &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_new                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_agd                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: soot_cld                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so2                                     & 
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_aitken                              &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_accu                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: so4_diss                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: nh3                                     &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dms                                     &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div1                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div2                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div3                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div4                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div5                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: dust_div6                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_new                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_agd                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: bmass_cld                               &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: ocff_new                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: ocff_agd                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: ocff_cld                                &
                         (tdims_s%i_start:tdims_s%i_end,              &
                          tdims_s%j_start:tdims_s%j_end,              &
                          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  :: nitr_acc                                &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT)  :: nitr_diss                               &
                       (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT(INOUT)  ::                                         &
       tracers(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
                tr_levels,tr_vars)                                    &
      , tracers_ukca(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end,tr_ukca)           &
! Add cariolle specific parameters for ozone tracer     
            , ozone_tracer(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end) 

      REAL                                                            &
        depart_lambda (wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end)                     &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
      , depart_phi (wdims%i_start:wdims%i_end,                        &
                    wdims%j_start:wdims%j_end,                        &
                    wdims%k_start:wdims%k_end)                        &
                                                    ! Phi Co-ordinate of
                                                    ! co-ordinate of
                                                    ! departure point.
      , depart_r (wdims%i_start:wdims%i_end,                          &
                  wdims%j_start:wdims%j_end,                          &
                  wdims%k_start:wdims%k_end)     ! Vertical
                                                 ! co-ordinate of
                                                 ! departure point.

      REAL                                                              &
        r_rho_levels (1-halo_i:row_length+halo_i,                       &
                      1-halo_j:rows+halo_j, model_levels)               &
      , r_theta_levels (1-halo_i:row_length+halo_i,                     &
                        1-halo_j:rows+halo_j, 0:model_levels)           &
      , eta_theta_levels(0:model_levels)

      REAL                                                              &
        p_star (row_length, rows)                                       &
      , p (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)  &
      , p_theta_levels (1-off_x:row_length+off_x,                       &
                        1-off_y:rows+off_y, model_levels)               &
      ,  rho   (1-off_x:row_length+off_x,                               &
                 1-off_y:rows+off_y, model_levels)                      &
      ,  drkp1, drk

      REAL                                                              &
                      ! Trig functions.
        FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
                               1-off_y:rows+off_y)

! look-up table halos
       INTEGER                                                          &
        halo_lam                                                        &
      , halo_phi

! Arguments with INTENT OUT. ie: Output variables.

      INTEGER                                                           &
        Error_Code     ! Non-zero on exit if error detected.

      LOGICAL                                                           &
        at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid



! Local Variables.

! scalars

      INTEGER                                                           &
        i, j, k                                                         &
                    ! Loop indices
      , temp                                                            &
      , counter                                                           &
      , tr_start                                                        &
      , lbc_index      ! pointer to tracer for active tracer lbcs

      INTEGER :: k_int_linear ! Linear interpolation is used at departure
                              ! points in this layer and below.
                              ! (Optional argument for subroutine
                              !  eg_interpolation_eta.)

      LOGICAL                                                           &
        L_do_halos                                                      &
                                 ! update the halos?
      , L_do_boundaries          ! update the boundaries?

! arrays

      REAL                                                              &
        work(row_length, rows, model_levels)                            &
      , work4(row_length, rows)

      INTEGER ::  super_array_size, array_size_count

      REAL                                                            &
        super_array(       tdims_l%i_start:tdims_l%i_end,             &
                           tdims_l%j_start:tdims_l%j_end,             &
                           tdims_l%k_start:tdims_l%k_end,             &
                           super_array_size)                          &
      , data_out_super(    tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                           tdims%k_start:tdims%k_end,                 &
                           super_array_size)

      INTEGER                                                           &
        i_out (row_length, rows, model_levels)                          &
      , j_out (row_length, rows, model_levels)

      REAL                                                              &
        weight_lambda (row_length, rows, model_levels)                  &
      , weight_phi    (row_length, rows, model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
!  Section 0.    Initialise array_size_count
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SL_TRACER1_4A',zhook_in,zhook_handle)
      array_size_count=0

! ----------------------------------------------------------------------
! Section 3.1  Set delta_p at data points, ensure positive.
! ----------------------------------------------------------------------
! Store in work

       IF(l_tracer1_non_hydro)THEN

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, drkp1, drk)

        k = 1
!$OMP DO  SCHEDULE(STATIC) 
        DO j = 1, rows
          DO i = 1, row_length
            drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
            drk   = r_theta_levels(i,j,k)   - r_theta_levels(i,j,k-1)
            work(i,j,k) = rho(i,j,k+1)*drkp1 + rho(i,j,k)  *drk
            work(i,j,k) = work(i,j,k)/r_theta_levels(i,j,k)**2
          END DO
        END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) 
        DO k = 2, model_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
              drk   = r_rho_levels(i,j,k)     - r_theta_levels(i,j,k-1)
              work(i,j,k) = rho(i,j,k+1)*drkp1  + rho(i,j,k)*drk
              work(i,j,k) = work(i,j,k)/r_theta_levels(i,j,k)**2
            END DO
          END DO
        END DO
!$OMP END DO

        k = model_levels
!$OMP DO SCHEDULE(STATIC) 
        DO j = 1, rows
          DO i = 1, row_length
              drk   = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              work(i,j,k)   = rho(i,j,k)  *drk
              work(i,j,k) = work(i,j,k)/r_theta_levels(i,j,k)**2
          END DO
        END DO
!$OMP END DO

!$OMP END PARALLEL

       ELSE

        k = 1
        DO j = 1, rows
          DO i = 1, row_length
            work(i,j,k) = p_star(i,j) - p(i,j,2)
          END DO
        END DO
        DO k = 2, model_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              work(i,j,k) = p(i,j,k) - p(i,j,k+1)
            END DO
          END DO
        END DO
        k = model_levels
        DO j = 1, rows
          DO i = 1, row_length
            work(i,j,k) = p(i,j,k) - p_theta_levels(i,j,k)
          END DO
        END DO
       END IF         ! l_tracer1_non_hydro

! ----------------------------------------------------------------------
! Section 3.2.1  carbon cycle.
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

      IF(l_CO2_interactive)THEN

!$OMP MASTER 
        array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count) = CO2(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO 
      END IF  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 3.2.2.1  Soot cycle.
! ----------------------------------------------------------------------
      IF (l_Soot) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count) =  soot_new(i,j,k)
              super_array(i,j,k,array_size_count+1) =  soot_agd(i,j,k)
              super_array(i,j,k,array_size_count+2) =  soot_cld(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO 

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER

!$OMP BARRIER
      END IF  ! l_soot

! ----------------------------------------------------------------------
! Section 3.2.2.2  Biomass aerosol.
! ----------------------------------------------------------------------
      IF (l_Biomass) THEN

!$OMP MASTER 
        array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count) =  bmass_new(i,j,k)
              super_array(i,j,k,array_size_count+1) =  bmass_agd(i,j,k)
              super_array(i,j,k,array_size_count+2) =  bmass_cld(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER

!$OMP BARRIER
      END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 3.2.3  sulphur cycle.
! ----------------------------------------------------------------------
      IF (l_Sulpc_so2) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER
        
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count) =  so4_aitken(i,j,k)
              super_array(i,j,k,array_size_count+1) =  so4_accu(i,j,k)
              super_array(i,j,k,array_size_count+2) =  so4_diss(i,j,k)
              super_array(i,j,k,array_size_count+3) =  so2(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO 

!$OMP MASTER
        array_size_count=array_size_count +3
!$OMP END MASTER

!$OMP BARRIER

        IF(L_sulpc_nh3)THEN

!$OMP MASTER 
          array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_array(i,j,k,array_size_count) =  nh3(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        END IF
        IF(L_sulpc_dms)THEN

!$OMP MASTER
          array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_array(i,j,k,array_size_count) =  dms(i,j,k)
              END DO
            END DO
          END DO
!$OMP END DO
        END IF

      END IF  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 3.2.4  mineral dust.
! ----------------------------------------------------------------------
      IF (L_DUST) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO J = 1, ROWS
            DO I = 1, ROW_LENGTH
              super_array(i,j,k,array_size_count) =  dust_div1(i,j,k)
              super_array(i,j,k,array_size_count+1) =  dust_div2(i,j,k)
              IF(.NOT.l_twobin_dust) THEN
                super_array(i,j,k,array_size_count+2) =  dust_div3(i,j,k)
                super_array(i,j,k,array_size_count+3) =  dust_div4(i,j,k)
                super_array(i,j,k,array_size_count+4) =  dust_div5(i,j,k)
                super_array(i,j,k,array_size_count+5) =  dust_div6(i,j,k)
              END IF
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        IF(l_twobin_dust) THEN
          array_size_count=array_size_count + 1
        ELSE
          array_size_count=array_size_count + 5
        END IF
!$OMP END MASTER

!$OMP BARRIER

      END IF ! L_DUST

!$OMP END PARALLEL

! ----------------------------------------------------------------------
! New addition   Fossil-fuel organic carbon
! ----------------------------------------------------------------------
      IF (L_ocff) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count)   =  ocff_new(i,j,k)
              super_array(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)
              super_array(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +2
      END IF  ! L_ocff

! ----------------------------------------------------------------------
! Section 3.2.4.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (L_USE_CARIOLLE) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO J = 1, ROWS
            DO I = 1, ROW_LENGTH
              super_array(i,j,k,array_size_count)   = OZONE_TRACER(i,j,k)
            END DO
          END DO
        END DO
      END IF ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! Section 3.2.4.2   Ammonium nitrate aerosol
! ----------------------------------------------------------------------
      IF (L_nitrate) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count)   =  nitr_acc(i,j,k)
              super_array(i,j,k,array_size_count+1) =  nitr_diss(i,j,k)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +1
      END IF  ! L_nitrate

! ----------------------------------------------------------------------
! ANY NEW NAMED AEROSOL SPECIES SHOULD BE ADDED HERE ^^^
! and in the same location/order in sl_tracer2 and tr_set_phys 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 3.2.5  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
       IF ( tr_vars > 0 ) THEN
        DO counter=1,tr_vars
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
             super_array(i,j,k,array_size_count) =  tracers(i,j,k,counter)
              END DO
            END DO
          END DO
        END DO

       END IF ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 3.2.6  UKCA tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
       IF ( tr_ukca > 0 ) THEN
        DO counter=1,tr_ukca
          array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                super_array(i,j,k,array_size_count) =                   &
                                 tracers_ukca(i,j,k,counter)
              END DO
            END DO
          END DO
        END DO
       END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 3.2.99  Murk cycle. This must be the last Full level field in
!                             the super_array
! ----------------------------------------------------------------------
      IF (L_Murk_advect) THEN
        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,array_size_count) =  murk(i,j,k)
            END DO
          END DO
        END DO
      END IF  ! L_Murk_advect

! ----------------------------------------------------------------------
! Call SWAPBOUNDS and set external halos for all arrays
! ----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
      CALL Swap_Bounds(                                                 &
             super_array,                                               &
             row_length, rows,                                          &
             (tdims%k_end-tdims%k_start+1)*array_size_count,            &
             halo_i, halo_j, fld_type_p,  .FALSE. )

      IF (model_domain == mt_lam) THEN
! DEPENDS ON: set_external_halos
        CALL set_external_halos(                                        &
              super_array,                                              &
              row_length, rows,                                         &
              (tdims%k_end-tdims%k_start+1)*array_size_count,           &
              halo_i, halo_j, 0.0)
      END IF

!------------------------------------------------------------
! having defined array now do interpolation
!------------------------------------------------------------


      IF(super_array_size >= 1)THEN

          super_array(:,:,0,1:super_array_size)=                      &
          super_array(:,:,1,1:super_array_size)

! Set layers over which linear interpolation is used
          IF (l_sl_bc_correction) THEN
            k_int_linear=2
          ELSE
            k_int_linear=1
          END IF

          CALL eg_interpolation_eta(                                  &
                     eta_theta_levels,fld_type_w,                     &
                     super_array_size,                                &
                     row_length, rows,tdims%k_end-tdims%k_start+1,    &
                     rows,                                            &
                     row_length, rows,tdims%k_end-tdims%k_start+1,    &
                     high_order_scheme_moist, monotone_scheme_moist,  &
                     model_domain, l_high_moist,                      &
                     l_mono_moist,                                    &
                     depart_r,depart_lambda,depart_phi,               &
                     me, n_proc, n_procx, n_procy,                    &
                     halo_i, halo_j,                                  &
                     global_row_length, l_datastart, at_extremity,    &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, off_x ,off_y,error_code,                   &
                     super_array,data_out_super,                      &
                     k_int_linear_in=k_int_linear)

      END IF


! ----------------------------------------------------------------------
! Section 3.2.1  carbon cycle.
! ----------------------------------------------------------------------

      array_size_count=0

      IF(l_CO2_interactive)THEN

        array_size_count=array_size_count +1

        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              co2(i,j,k) = data_out_super(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 3.2.2.1  Soot cycle.
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i ,j ,k)
      IF (l_Soot) THEN

!$OMP MASTER 
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              soot_new(i,j,k) = data_out_super(i,j,k,array_size_count)
              soot_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
              soot_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
!$OMP END DO 

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER
      END IF  ! l_soot

! ----------------------------------------------------------------------
! Section 3.2.2.2  Biomass aerosol.
! ----------------------------------------------------------------------
      IF (l_Biomass) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
             bmass_new(i,j,k) = data_out_super(i,j,k,array_size_count)
             bmass_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
             bmass_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER

      END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 3.2.3  sulphur cycle.
! ----------------------------------------------------------------------
      IF (l_Sulpc_so2) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              so4_aitken(i,j,k) = data_out_super(i,j,k,array_size_count)
              so4_accu(i,j,k) = data_out_super(i,j,k,array_size_count+1)
              so4_diss(i,j,k) = data_out_super(i,j,k,array_size_count+2)
              so2(i,j,k) = data_out_super(i,j,k,array_size_count+3)
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +3
!$OMP END MASTER 
!$OMP BARRIER

        IF(L_sulpc_nh3)THEN

!$OMP MASTER
          array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                nh3(i,j,k) = data_out_super(i,j,k,array_size_count)
              END DO
            END DO
          END DO
!$OMP END DO
        END IF
        IF(L_sulpc_dms)THEN

!$OMP MASTER
          array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
            DO k = tdims%k_start, tdims%k_end
              DO j = 1, rows
                DO i = 1, row_length
                dms(i,j,k) = data_out_super(i,j,k,array_size_count)
                END DO
              END DO
            END DO
!$OMP END DO
        END IF

        END IF  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 3.2.4  mineral dust.
! ----------------------------------------------------------------------
        IF (L_DUST) THEN
!$OMP MASTER 
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
          DO k = tdims%k_start, tdims%k_end
            DO J = 1, ROWS
              DO I = 1, ROW_LENGTH
          dust_div1(i,j,k) = data_out_super(i,j,k,array_size_count)
          dust_div2(i,j,k) = data_out_super(i,j,k,array_size_count+1)
                IF(.NOT.l_twobin_dust) THEN
          dust_div3(i,j,k) = data_out_super(i,j,k,array_size_count+2)
          dust_div4(i,j,k) = data_out_super(i,j,k,array_size_count+3)
          dust_div5(i,j,k) = data_out_super(i,j,k,array_size_count+4)
          dust_div6(i,j,k) = data_out_super(i,j,k,array_size_count+5)
                END IF
              END DO
            END DO
          END DO
!$OMP END DO 
 
!$OMP MASTER
          IF(l_twobin_dust) THEN
            array_size_count=array_size_count + 1
          ELSE
            array_size_count=array_size_count + 5
          END IF
!$OMP END MASTER 
!$OMP BARRIER
        END IF ! L_DUST

! ----------------------------------------------------------------------
! New addition  Fossil-fuel organic carbon aerosol.
! ----------------------------------------------------------------------
      IF (L_ocff) THEN

!$OMP MASTER
        array_size_count=array_size_count +1
!$OMP END MASTER 
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
             ocff_new(i,j,k) = data_out_super(i,j,k,array_size_count)
             ocff_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
             ocff_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
            END DO
          END DO
        END DO
!$OMP END DO

!$OMP MASTER
        array_size_count=array_size_count +2
!$OMP END MASTER 
!$OMP BARRIER

      END IF  ! L_ocff

!$OMP END PARALLEL 
! ----------------------------------------------------------------------
! Section 3.2.4.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
        IF (L_USE_CARIOLLE) THEN

        array_size_count=array_size_count +1
          DO k = tdims%k_start, tdims%k_end
            DO J = 1, ROWS
              DO I = 1, ROW_LENGTH
          OZONE_TRACER(i,j,k) = data_out_super(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END IF ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! New addition  Ammonium nitrate aerosol
! ----------------------------------------------------------------------
      IF (L_nitrate) THEN

        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
             nitr_acc(i,j,k) = data_out_super(i,j,k,array_size_count)
             nitr_diss(i,j,k)= data_out_super(i,j,k,array_size_count+1)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +1

      END IF  ! L_nitrate

! ----------------------------------------------------------------------
! ANY NEW NAMED AEROSOL SPECIES SHOULD BE ADDED HERE ^^^
! and in the same location/order in sl_tracer2 and tr_set_phys
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 3.2.5  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF ( tr_vars > 0 ) THEN

        DO counter=1,tr_vars
        array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS 
            DO J = 1, ROWS
              DO I = 1, ROW_LENGTH
         tracers(i,j,k,counter) = data_out_super(i,j,k,array_size_count)
              END DO
            END DO
          END DO
        END DO
      END IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 3.2.6  UKCA tracers  (tr_levels = model_levels)
! ----------------------------------------------------------------------
      IF ( tr_ukca > 0 ) THEN

        DO counter = 1, tr_ukca
        array_size_count = array_size_count + 1
          DO k = tdims%k_start, tdims%k_end
            DO j = 1, rows
              DO i = 1, row_length
                tracers_ukca(i,j,k,counter) =                             &
                                 data_out_super(i,j,k,array_size_count)
              END DO
            END DO
          END DO  ! Loop over levels
        END DO !  counter = 1, tr_ukca
      END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 3.2.99  Murk cycle. This must be the last Full level field in
!                            the super_array
! ----------------------------------------------------------------------
      IF (L_Murk_advect) THEN

        array_size_count=array_size_count +1
        DO k = tdims%k_start, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              murk(i,j,k) = data_out_super(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      END IF  ! L_Murk_advect

! END of routine.
      IF (lhook) CALL dr_hook('SL_TRACER1_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SL_Tracer1_4A

END MODULE SL_Tracer1_mod
