! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_filter_Ctl

      Subroutine NI_filter_Ctl(                                         &
     &                      theta, u, v, w, Exner, rho,                 &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v, delta_lambda, delta_phi,    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      max_filter_rows, u_sweeps, v_sweeps,        &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      top_filt_start, top_filt_end,               &
     &                      up_diff, max_upd_levels,                    &
     &                      horizontal_level, me,                       &
     &                      global_row_length, global_rows,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, model_domain,      &
     &                      proc_row_group, proc_col_group,             &
     &                      L_polar_filter, L_pofil_new,                &
     &                      L_pfcomb, L_pftheta, L_pfuv,                &
     &                      L_pfw, L_pfexner, L_diff_exner,             &
     &                      L_diff_thermo, L_diff_wind, L_diff_w,       &
     &                      L_pofil_hadgem2, exner_theta_levels,        &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &                      STASHwork)

! Purpose:
!          Filter fields on rows adjacent to pole and
!          set polar v components to wave 1.
!
! Method:
!          Using 1-2-1 Shapiro filter to zap 2-grid waves
!          T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      Use swapable_field_mod, Only :                                    &
          swapable_field_pointer_type
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
! Modules for z-level diffusion
      USE atm_fields_bounds_mod
      USE level_heights_mod,    ONLY : eta_theta_levels, eta_rho_levels 
      USE pofil_zlevel_mod
      USE conversions_mod,      ONLY: pi_over_180
      USE Submodel_Mod

      Implicit None

! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_polar_filter                                                  &
                       ! switch for polar filter
     &, L_pofil_new                                                     &
                       ! switch for new polar filter
     &, L_pfcomb                                                        &
                       ! switch for combined polar filter/diffusion
     &, L_pftheta                                                       &
                       ! switch for polar filter for theta
     &, L_pfw                                                           &
                       ! switch for polar filter for w
     &, L_pfuv                                                          &
                       ! switch for polar filter for horizontal winds
     &, L_pfexner                                                       &
                       ! switch for polar filter for Exner pressure
     &, L_diff_Exner                                                    &
                       ! switch for diffusion of Exner pressure
     &, L_diff_thermo                                                   &
                       ! switch for horiz. diffusion of theta
     &, L_diff_wind                                                     &
                       ! switch for horiz. diffusion of u,v
     &, L_diff_w       ! switch for horiz. diffusion of w
      Logical :: L_pofil_hadgem2  ! use hadgem2 polar filtering settings

      Integer                                                           &
     &  max_filter_rows                                                 &
                           ! max array size for u_begin etc
     &, max_upd_levels                                                  &
                         ! max no. levels for upper diffusion
     &, global_u_filter                                                 &
                         ! number of filter sweeps; 0=diffusion only
     &, global_v_filter                                                 &
                         ! number of filter sweeps; 0=diffusion only
     &, row_length                                                      &
                         ! number of points on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
     &, global_row_length                                               &
     &, first_constant_r_rho_level                                      &
     &, first_constant_r_rho_level_m1                                   &
     &, global_rows                                                     &
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, l_datastart(3)                                                  &
                             ! First gridpoints held by this processor
     &, me                                                              &
                  !  this processor
     &, model_domain

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)
      Real                                                              &
     &  cos_theta_longitude(row_length,rows)                            &
                                             ! cos of longitude
     &, sin_theta_longitude(row_length,rows)                            &
                                             ! sin of longitude
     &, sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, cos_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, cos_v_latitude(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y ) &
     &, sec_v_latitude(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y ) &
     &, up_diff(max_upd_levels)                                         &
                                 !  upper-level diffusion coefficients
!    diffusion coefficient for u/theta rows
     &, diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y )     &
!    diffusion coefficient for u/theta rows
     &, diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y )

      Real                                                              &
           ! limits for polar filtering (in radians)
     &  north_lat_limit                                                 &
     &, south_lat_limit                                                 &
     &, polar_filter_coefficient                                        &
                                       !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit                                           &
                              ! max latitude at which filter can start
     &, diff_coeff_phi                                                  &
                         ! NS diffusion coeff
     &, diff_coeff_thermo                                               &
                            ! NS diffusion coeff thermo
     &, diff_coeff_wind     ! NS diffusion coeff u,v



      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, diff_order_thermo                                               &
                                ! diffusion order
     &, diff_order_wind                                                 &
                              ! diffusion order
     &, proc_col_group                                                  &
     &, polar_filter_n_sweeps                                           &
                                     ! number of sweeps of filter to do
     &, horizontal_level                                                &
                                     ! steep slope control
     &, top_filt_start                                                  &
                             ! start level for upper-level diffusion
     &, top_filt_end         ! end level for upper-level diffusion

      Integer                                                           &
     &  neighbour(4)                                                    &
     &, u_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, v_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, u_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, u_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, v_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, v_end(0:max_filter_rows)    ! row pointers for 1-2-1 filter

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN/OUT.

      Real, Target ::                                                   &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)                                                &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)

      Real ::                                                           &
     &  exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels + 1)                                        &
     &, rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &         model_levels )

!   Array arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields


! Local variables
      Integer                                                           &
     &  i, j, k                                                         &
                          ! loop counters
     &, j_start, j_stop                                                 &
                          ! Loop bounds
     &, gi                                                              &
                       ! pointer relative to origin
     &, active_levels  ! number of levels for upper-level diffusion

      Integer :: i_field  ! counter for swapbounds

      Real                                                              &
     &  pole_term     ! factor at pole calculated once

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)

      Real,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     &  r_theta_at_u                                                    &
     &, r_theta_at_v                                                    &
     &, r_uv_b

      Real :: u_inc (row_length, rows, model_levels) 
      Real :: v_inc (row_length, n_rows, model_levels) 
      Real :: w_inc (row_length, rows, 0:model_levels) 
      Real :: T_inc (row_length, rows, model_levels) 
      Real :: exner_inc (row_length, rows, model_levels+1) 
      Integer :: sect
      Integer :: item
      Integer :: im_index      !  internal model index for STASH arrays
      INTEGER :: Errorstatus = 0  ! initial value for error code
      CHARACTER (LEN=80) :: CMessage !  Error message

      Type(swapable_field_pointer_type) :: fields_to_swap(4)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      REAL ::                                                           &
        diff_coeff_phi_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y), &
                            !    NS diffusion coefficient for u/theta rows
        diff_coeff_phi_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)
                            !    NS diffusion coefficient for v rows

      REAL :: xi1_p(1-halo_i:row_length+halo_i),                        &
              xi1_u(1-halo_i:row_length+halo_i),                        &
              xi2_p(1-halo_j:rows+halo_j),                              &
              xi2_v(1-halo_j:n_rows+1+halo_j)

      REAL :: csphi_theta(1-halo_j:rows+halo_j),                        &
              csphi_v(1-halo_j:n_rows+1+halo_j)

      Real,DIMENSION(:,:,:),ALLOCATABLE :: r_at_v_new  

      Integer :: istat
      Real    :: temp
      Real    :: k_sum(row_length,rows),  k_sum_v(row_length,n_rows) 
      Real    :: row_sum(rows) 

      INTEGER :: gj
      REAL    :: base_lambda, base_phi

! ----------------------------------------------------------------------
! 0.1  Section for setting up diagnostic arrays
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('NI_FILTER_CTL',zhook_in,zhook_handle)
      sect=13
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''

      IF (sf(0,sect)) THEN
        IF(sf(381,sect))t_inc(:,:,:)= theta(1:row_length, 1:rows, :)
        IF(sf(385,sect))u_inc(:,:,:)=u(1:row_length, 1:rows, :)
        IF(sf(386,sect))v_inc(:,:,:)=v(1:row_length, 1:n_rows, :)
        IF(sf(387,sect))w_inc(:,:,:)=w(1:row_length, 1:rows, :)
        IF(sf(388,sect))exner_inc(:,:,:)= exner(1:row_length, 1:rows, :)
      ENDIF

! ----------------------------------------------------------------------
! 1.0  Section for original polar filter
! ----------------------------------------------------------------------

      If (L_polar_filter) Then

! DEPENDS ON: polar_filter
        Call Polar_filter(                                              &
     &                      theta, u, v, w,                             &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      global_row_length, delta_lambda,            &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, proc_row_group)

      End If    !  L_polar_filter

! ----------------------------------------------------------------------
! 2.0  Section for new polar filter
! ----------------------------------------------------------------------

      ALLOCATE ( r_theta_at_u(1-off_x:row_length+off_x,                 &
     &                        1-off_y:rows+off_y, 0:model_levels) )
      ALLOCATE ( r_theta_at_v(1-off_x:row_length+off_x,                 &
     &                        1-off_y:n_rows+off_y, 0:model_levels) )

      Do k = 0, model_levels
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
! end loop bound OK since r_theta_levels is halo_i > off_x
            r_theta_at_u(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
     &                                  r_theta_levels(i  ,j,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

      Do k = 0, model_levels
        Do j = 1-off_y, n_rows+off_y
! end loop bound OK since r_theta_levels is halo_j > off_y
          Do i = 1-off_x, row_length+off_x
            r_theta_at_v(i,j,k) = .5 * ( r_theta_levels(i,j+1,k) +      &
     &                                   r_theta_levels(i,j  ,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

! ----------------------------------------------------------------------
! 2.1  polar filter/diffusion theta/w  polar filter Exner
! ----------------------------------------------------------------------

! Swap_bounds are done inside the polar filter sweeps
      if( L_pftheta .or. L_diff_thermo )then
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                theta, fld_type_p, 0, 0, 0, 1, 0, 0,              &
     &                model_levels, model_levels, model_levels,         &
     &                first_constant_r_rho_level_m1,                    &
     &                rows, n_rows, rows, row_length,                   &
     &                r_theta_levels, r_rho_levels,                     &
     &                r_theta_at_u, r_theta_at_v,                       &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_thermo, .false.,             &
     &                L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_th_field
        Call pofil_th_field(                                            &
     &                  theta, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end, u_sweeps,      &
     &                  global_u_filter, model_levels, horizontal_level,&
     &                  diff_order_thermo, diff_coeff_thermo,           &
     &                  diff_coeff_u, L_diff_thermo )
      endif !  L_pofil_new
      endif !L_pftheta .or. L_diff_thermo

      if( L_pfw .or. L_diff_w )then
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                w(1-off_x,1-off_y, 1), fld_type_p,                &
     &                0, 0, 0, 1, 0, 0,                                 &
     &                model_levels, model_levels, model_levels-1,       &
     &                first_constant_r_rho_level_m1,                    &
     &                rows, n_rows, rows, row_length,                   &
     &                r_theta_levels, r_rho_levels,                     &
     &                r_theta_at_u, r_theta_at_v,                       &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_w, .false.,                  &
     &                L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_th_field
        Call pofil_th_field(                                            &
     &                  w(1-off_x,1-off_y, 1),                          &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end, u_sweeps,      &
     &                  global_u_filter,                                &
     &                  model_levels - 1, horizontal_level,             &
     &                  diff_order_thermo, diff_coeff_thermo,           &
     &                  diff_coeff_u, L_diff_w )
      endif !  L_pofil_new
      endif ! L_pfw .or. L_diff_w

      if( L_pfexner .and. L_pofil_new )then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                Exner, fld_type_p, 0, 0, 1, 0, 1, 1,              &
     &                model_levels+1, model_levels, model_levels+1,     &
     &                first_constant_r_rho_level,                       &
     &                rows, n_rows, rows, row_length,                   &
     &                r_rho_levels, r_theta_levels,                     &
     &                r_at_u, r_at_v,                                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_Exner, .false.,              &
     &                L_pofil_hadgem2 )
      endif !  L_pfexner .and. L_pofil_new

      if( L_pfuv .or. L_diff_wind )then

        ALLOCATE ( r_uv_b(1-off_x:row_length+off_x,                     &
     &                    1-off_y:n_rows+off_y, model_levels) )

!  r needed at centre of grid but no swap bound option for grid-type
!      so fill required halos explicitly
!      (fortunately have large halos to start)
        Do k = 1, model_levels
          Do j = 1-off_y, n_rows+off_y
! end loop bound OK since r_rho_levels is halo_j > off_y
            Do i = 1-off_x, row_length+off_x
              r_uv_b(i,j,k) = .5 * ( r_at_u(i,j,k) + r_at_u(i,j+1,k) )
            End Do
          End Do
        End Do ! k = 1, model_levels

! ----------------------------------------------------------------------
! 2.2  polar filter/diffusion u/v
! ----------------------------------------------------------------------
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                 u, fld_type_u, 1, 0, 1, 0, 1, 1,                 &
     &                 model_levels, model_levels, model_levels,        &
     &                 first_constant_r_rho_level,                      &
     &                 rows, n_rows, rows, row_length,                  &
     &                 r_at_u, r_theta_at_u, r_rho_levels, r_uv_b,      &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 sec_theta_latitude, cos_v_latitude, 0.0,         &
     &                 model_domain, at_extremity, n_procy,             &
     &                 max_filter_rows, global_u_filter,                &
     &                 u_sweeps, u_begin, u_end,                        &
     &                 horizontal_level, diff_coeff_phi,                &
     &                 diff_coeff_u, L_diff_wind, .true.,               &
     &                 L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_u
        Call pofil_u(                                                   &
     &              u, r_at_u, r_theta_levels, r_rho_levels,            &
     &              off_x, off_y, halo_i, halo_j,                       &
     &              sec_theta_latitude, cos_v_latitude,                 &
     &              me, n_procy, delta_lambda, delta_phi,               &
     &              rows, n_rows, row_length, model_levels,             &
     &              max_filter_rows, u_begin, u_end, u_sweeps,          &
     &              global_u_filter,                                    &
     &              horizontal_level, diff_order_wind,                  &
     &              diff_coeff_wind, diff_coeff_u, L_diff_wind )
      endif !  L_pofil_new

      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
      Call pofil_new(                                                   &
     &               v, fld_type_v, 0, 1, 1, 0, 1, 1,                   &
     &               model_levels, model_levels, model_levels,          &
     &               first_constant_r_rho_level,                        &
     &               n_rows, rows, rows, row_length,                    &
     &               r_at_v, r_theta_at_v, r_uv_b, r_rho_levels,        &
     &               off_x, off_y, off_x, off_y,                        &
     &               halo_i, halo_j, off_x, off_y,                      &
     &               off_x, off_y, halo_i, halo_j,                      &
     &               sec_v_latitude, cos_theta_latitude, 0.0,           &
     &               model_domain, at_extremity, n_procy,               &
     &               max_filter_rows, global_v_filter,                  &
     &               v_sweeps, v_begin, v_end,                          &
     &               horizontal_level, diff_coeff_phi,                  &
     &               diff_coeff_v, L_diff_wind, .true.,                 &
     &               L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_v
      Call pofil_v(                                                     &
     &              v, r_at_v, r_theta_levels, r_rho_levels,            &
     &              off_x, off_y, halo_i, halo_j,                       &
     &              sec_v_latitude, cos_theta_latitude,                 &
     &              l_datastart, global_row_length,                     &
     &              delta_lambda, delta_phi,                            &
     &              me, at_extremity, proc_row_group, n_procy,          &
     &              rows, n_rows, row_length, model_levels,             &
     &              max_filter_rows, v_begin, v_end, v_sweeps,          &
     &              global_v_filter,                                    &
     &              horizontal_level, diff_order_wind,                  &
     &              diff_coeff_wind, diff_coeff_v, L_diff_wind )
      endif !  L_pofil_new

      DEALLOCATE ( r_uv_b )
      DEALLOCATE ( r_theta_at_u )
      DEALLOCATE ( r_theta_at_v )

! ----------------------------------------------------------------------
! Section 2.3  Set u wind at poles since v has changed
!              only for global and if upper level diffusion is inactive
! ----------------------------------------------------------------------

      If ( model_domain == mt_global .and.                              &
     &     top_filt_start > model_levels ) Then

! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                           v,                                     &
     &                           sin_theta_longitude,                   &
     &                           cos_theta_longitude, row_length,       &
     &                           n_rows, model_levels,                  &
     &                           mag_vector_np, dir_vector_np,          &
     &                           mag_vector_sp, dir_vector_sp,          &
     &                           off_x, off_y, global_row_length,       &
     &                           proc_row_group, at_extremity)

        If (at_extremity(PSouth) ) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              u(i,1,k) = - mag_vector_sp(k) * sin ( (gi-.5)*            &
     &                     delta_lambda - dir_vector_sp(k))
            End Do
          End Do  ! k = 1, model_levels
        End If   !  at_extremity(PSouth)
        If (at_extremity(PNorth) ) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              u(i,rows,k) = mag_vector_np(k) *                          &
     &                          sin ( (gi-.5)*delta_lambda -            &
     &                              dir_vector_np(k))
            End Do
          End Do ! k = 1, model_levels
        End If  ! at_extremity(PNorth

      endIf ! model_domain = mt_global and top_filt_start > model_levels


! u field has been changed at poles, therefore need swap_bounds call
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   u, row_length, rows, model_levels,             &
     &                   off_x, off_y, fld_type_u, .true.)

      endif ! L_pfuv .or. L_diff_wind

! ----------------------------------------------------------------------
! 3.0  Section for upper-level diffusion
! ----------------------------------------------------------------------

      If ( top_filt_start < model_levels + 1 ) then

        active_levels = top_filt_end - top_filt_start + 1
        pole_term = 8. / ( global_row_length * delta_phi )
        j_start = 1
        j_stop = rows
        If (model_domain == mt_Global) Then
          If (at_extremity(PSouth)) j_start = 2
          If (at_extremity(PNorth)) j_stop = rows - 1
        End If

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 theta, fld_type_p, 0,                            &
     &                 model_levels, active_levels,                     &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 w(1-off_x, 1-off_y, 1), fld_type_p, 0,           &
     &                 model_levels - 1, active_levels - 1,             &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 u, fld_type_u, 0,                                &
     &                 model_levels, active_levels,                     &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

        j_start = 1
        j_stop = n_rows

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 v, fld_type_v, 1,                                &
     &                 model_levels, active_levels,                     &
     &                 n_rows, rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! ----------------------------------------------------------------------
! Section 3.1  Set u wind at poles since v has changed
! ----------------------------------------------------------------------

        If (model_domain == mt_Global) Then

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                             v,                                   &
     &                             sin_theta_longitude,                 &
     &                             cos_theta_longitude, row_length,     &
     &                             n_rows, model_levels,                &
     &                             mag_vector_np, dir_vector_np,        &
     &                             mag_vector_sp, dir_vector_sp,        &
     &                             off_x, off_y, global_row_length,     &
     &                             proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,1,k) = - mag_vector_sp(k) * sin ( (gi-.5)*          &
     &                       delta_lambda - dir_vector_sp(k))
              End Do
            End Do  ! k = 1, model_levels
          End If   !  at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,rows,k) = mag_vector_np(k) *                        &
     &                            sin ( (gi-.5)*delta_lambda -          &
     &                                dir_vector_np(k))
              End Do
            End Do ! k = 1, model_levels
          End If  ! at_extremity(PNorth

        endIf !model_domain == mt_Global

        i_field = 0

        i_field = i_field + 1
        fields_to_swap(i_field) % field      => theta(:,:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = model_levels 
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field      => w(:,:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = model_levels+1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field      => u(:,:,:)
        fields_to_swap(i_field) % field_type = fld_type_u
        fields_to_swap(i_field) % levels     = model_levels
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field      => v(:,:,:)
        fields_to_swap(i_field) % field_type = fld_type_v
        fields_to_swap(i_field) % levels     = model_levels
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .TRUE.

! DEPENDS ON: swap_bounds_mv
        CALL swap_bounds_mv(fields_to_swap, i_field, row_length,        &
                            off_x, off_y)

      endif !  top_filt_start < model_levels + 1

!*******************  increments for STASH   **************************

      IF(sf(0,sect)) THEN

! T increment
        item = 381 
        IF(sf(item,sect)) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_inc(i,j,k) = (theta(i,j,k) - T_inc(i,j,k))            &
     &                                  *exner_theta_levels(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        T_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! u wind increment
        item=385
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_inc(i,j,k) = u(i,j,k) - u_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        u_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! v wind increment
        item=386
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_inc(i,j,k) = v(i,j,k) - v_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        v_inc,                                                    &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! w wind increment
        item=387
        IF( sf(item,sect) ) THEN

          DO k=0,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_inc(i,j,k) = w(i,j,k) - w_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        w_inc,                                                    &
     &        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! exner increment
        item=388
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels+1
            DO j=1,rows
              DO i=1,row_length
                exner_inc(i,j,k) = exner(i,j,k) - exner_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        exner_inc,                                                &
     &        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

      ENDIF ! sf(0,sect) 

      IF (lhook) CALL dr_hook('NI_FILTER_CTL',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE NI_filter_Ctl

