! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! purpose: Interface to Atmospheric Physics GWD Schemes.
!         Scheme 1: Orographic flow blocking and gravity wave scheme.
!         Scheme 2: Non-orographic ultra-simple spectral gw scheme.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: gravity wave drag


MODULE NI_gwd_ctl_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE NI_gwd_ctl (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y                                    &
     &,  global_row_length,n_proc, n_procy, proc_row_group              &
     &, at_extremity, neighbour                                         &

! model dimensions.
     &, row_length, rows, n_rows, land_points                           &
     &, model_levels                                                    &

! Model switches
     &, model_domain                                                    &
! trig arrays
     &, sin_theta_longitude, sin_theta_latitude                         &

! in coordinate information
     &, delta_lambda,delta_phi,true_latitude                            &
     &, exner_theta_levels                                              & 

! in time stepping information.
     &, timestep, timestep_number                                       &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork                                                        &
! SCM diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &
! in data fields.
     &, u, v                                                            &
     &, land_sea_mask, p_layer_boundaries                               &
     &, rho, theta_latest, sd_orog_land, orog_grad_xx_land              &
     &, orog_grad_xy_land, orog_grad_yy_land, land_index                &

! in/out
     &, R_u, R_v, T_inc                                                 &

! error information
     &, Error_code  )

! Definitions of prognostic variable array sizes
      USE atm_fields_bounds_mod, ONLY:                                  &
     &   udims, vdims, wdims, tdims, pdims                              &
     &,  udims_s, vdims_s, wdims_s, tdims_s, pdims_s                     

! Model level heights from centre of Earth
      USE level_heights_mod, ONLY:  &
     &   r_theta_levels             &  ! Radii on theta levels (m)
     &,  r_rho_levels                  ! Radii on rho levels (m)

      USE g_wave_input_mod, ONLY : l_gwd,                &      
                                   L_use_ussp,           &      
                                   l_gwd_40km,           &
                                   l_nonhydro,           &
                                   l_dynbeta,            &
                                   l_taus_scale,         &
                                   l_fix_gwsatn,         &
                                   sat_scheme,           &
                                   i_gwd_vn,             &      
                                   i_gwd_vn_4a,          &
                                   i_gwd_vn_5a,          &
                                   kay_gwave,            &
                                   gwd_frc,              &
                                   gwd_fsat,             &
                                   gsharp,               &
                                   fbcd,                 &
                                   l_smooth,             &
                                   l_gw_heating
                                   
      USE g_wave_4a_mod, ONLY: g_wave_4a
      USE g_wave_5a_mod, ONLY: g_wave_5a
      USE gw_ussp_mod, ONLY : gw_ussp

      USE ereport_mod, ONLY : Ereport

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Submodel_Mod
      IMPLICIT NONE

!     Fixed starting theta-level (e.g. for P_layer_boundaries)
      INTEGER, PARAMETER :: tkfix0start    = 0
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
      INTEGER, PARAMETER :: tkfix1start    = 1
!
! arguments with intent in. ie: input variables.

! Parallel setup variables
      INTEGER                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                   ! number of points on a row
     &, proc_row_group                                                  &
                   ! Group id for processors on the same row
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)   ! Array with the Ids of the four neighbours
                       ! in the horizontal plane

      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, south
                         ! east or west of the processor grid


! Model dimensions
      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, land_points

! Model switches
      INTEGER                                                           &
     &  model_domain

! model parameters
      REAL                                                              &
     &  timestep                                                       
  
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

! Diagnostics info
       REAL                                                             &
     & STASHwork(*) ! STASH workspace for section 6 (GW Drag)

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN)  ::  &
    nSCMDpkgs                ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN)  ::  &
    L_SCMDiags(nSCMDpkgs)    ! Logicals for SCM diagnostics packages

! Data arrays

      INTEGER                                                           &
     &  land_index (land_points)            ! set from land_sea_mask

      REAL ::                                                           &
     & u(udims_s%i_start:udims_s%i_end,                                 &
     &   udims_s%j_start:udims_s%j_end,                                 &
     &   udims_s%k_start:udims_s%k_end)                                 &
            !primary u field (ms**-1)
     &,v(vdims_s%i_start:vdims_s%i_end,                                 &
     &   vdims_s%j_start:vdims_s%j_end,                                 &
     &   vdims_s%k_start:vdims_s%k_end)                                 & 
            !primary v field (ms**-1)
     &,rho(pdims_s%i_start:pdims_s%i_end,                               &  
           pdims_s%j_start:pdims_s%j_end,                               & 
           pdims_s%k_start:pdims_s%k_end)  
            !density *r*r (kg/m)
      REAL                                                              &
     &  theta_latest(tdims%i_start:tdims%i_end,                         &
     &               tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end)

      REAL                                                              &
     &  p_layer_boundaries(tdims%i_start:tdims%i_end,                   &
     &                     tdims%j_start:tdims%j_end,                   &
     &                     tkfix0start:tdims%k_end)                     &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, sd_orog_land (land_points)                                      &
                                   ! orog/qrparm.orog.stdev
     &, orog_grad_xx_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxx
     &, orog_grad_xy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxy
     &, orog_grad_yy_land(land_points) ! orog/qrparm.orog.sigmayy


      REAL ::                                                           &
     & exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                & 
     &                    tdims_s%j_start:tdims_s%j_end,                & 
     &                    tdims_s%k_start:tdims_s%k_end)   
                                       ! Exner on theta level

      LOGICAL                                                           &
     &  land_sea_mask(row_length, rows)

! Co-ordinate arrays
      REAL                                                              &
     ! local vertical co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, true_latitude(row_length, rows)                      

      REAL                                                              &
     &  sin_theta_longitude (row_length, rows)                          &
     &, sin_theta_latitude  (row_length, rows)

! time information for current timestep
      INTEGER                                                           &
     &  timestep_number


! arguments with intent in/out. ie: input variables changed on output.

! arguments with intent out. ie: output variables.
      REAL       ::                                                     &
     &  r_u(udims_s%i_start:udims_s%i_end,                              & 
     &      udims_s%j_start:udims_s%j_end,                              &
     &      udims_s%k_start:udims_s%k_end)                              &
                                           !u wind increment diagnostic
     &, r_v(vdims_s%i_start:vdims_s%i_end,                              & 
     &      vdims_s%j_start:vdims_s%j_end,                              &
     &      vdims_s%k_start:vdims_s%k_end)                              &
                                           !v wind increment diagnostic
     &, T_inc(tdims%i_start:tdims%i_end,                                &
     &        tdims%j_start:tdims%j_end,                                &
     &        1:tdims%k_end)               !Temperature increment       

      INTEGER                                                           &
     &  Error_code


! local variables
      LOGICAL                                                           &
     &  stress_ud_on                                                    &
     &, stress_ud_p_on                                                  &
     &, stress_vd_on                                                    &
     &, stress_ud_satn_on                                               &
     &, stress_vd_satn_on                                               &
     &, stress_ud_wake_on                                               &
     &, stress_vd_wake_on                                               &
     &, du_dt_satn_on                                                   &
     &, du_dt_satn_p_on                                                 &
     &, dv_dt_satn_on                                                   &
     &, du_dt_wake_on                                                   &
     &, dv_dt_wake_on                                                   &
     &, GWSPEC_EFLUX_ON                                                 &
     &, GWSPEC_EFLUX_P_ON                                               &
     &, GWSPEC_SFLUX_ON                                                 &
     &, GWSPEC_WFLUX_ON                                                 &
     &, GWSPEC_WFLUX_P_ON                                               &
     &, GWSPEC_NFLUX_ON                                                 &
     &, GWSPEC_EWACC_ON                                                 &
     &, GWSPEC_EWACC_P_ON                                               &
     &, GWSPEC_NSACC_ON                                                 &
     &, u_s_d_on                                                        &
     &, v_s_d_on                                                        &
     &, nsq_s_d_on                                                      &
     &, fr_d_on                                                         &
     &, bld_d_on                                                        &
     &, bldt_d_on                                                       &
     &, num_lim_d_on                                                    &
     &, num_fac_d_on                                                    &
     &, tausx_d_on                                                      &
     &, tausy_d_on                                                      &
     &, taus_scale_d_on                                                 &
     &, orog_slope_d_on                                                 &
     &, orog_anis_d_on                                                  &
     &, orog_dir_d_on                                                   &
     &, L_u_incr_gwd                                                    &
     &, L_v_incr_gwd                                                    &
      , L_t_incr_gwd

      INTEGER                                                           &
        i,j,k      ! loop counters 

      INTEGER                                                           &
        errorstatus      ! Return code : 0 Normal Exit : >0 Error

      CHARACTER(LEN=256)                                                &
        cmessage         ! Error message if return code >0


! Local data arrays

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      REAL,DIMENSION(:,:,:),ALLOCATABLE::                               &
     &  u_incr_diagnostic                                               &
                              ! u wind increment for STASH
     &, v_incr_diagnostic                                               &
                              ! v wind increment for STASH
     &, t_incr_diagnostic                                               &
                              ! temperature increment for STASH
     &, stress_ud                                                       &
     &, stress_vd                                                       &
     &, stress_ud_satn                                                  &
     &, stress_vd_satn                                                  &
     &, stress_ud_wake                                                  &
     &, stress_vd_wake                                                  &
     &, du_dt_satn                                                      &
     &, dv_dt_satn                                                      &
     &, du_dt_wake                                                      &
     &, dv_dt_wake                                                      &
     &, GWSPEC_EFLUX                                                    &
     &, GWSPEC_SFLUX                                                    &
     &, GWSPEC_WFLUX                                                    &
     &, GWSPEC_NFLUX                                                    &
     &, GWSPEC_EWACC                                                    &
     &, GWSPEC_NSACC

      REAL,DIMENSION(:,:),ALLOCATABLE::                                 &
     &  u_s_d                                                           &
     &, v_s_d                                                           &
     &, nsq_s_d                                                         &
     &, fr_d                                                            &
     &, bld_d                                                           &
     &, bldt_d                                                          &
     &, num_lim_d                                                       &
     &, num_fac_d                                                       &
     &, tausx_d                                                         &
     &, tausy_d                                                         &
     &, taus_scale_d                                                    &
     &, orog_slope_d                                                    &
     &, orog_anis_d                                                     & 
     &, orog_dir_d

! Diagnostic land_point array sizes
      INTEGER                                                           &
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_satn                                            &
     &,points_stress_vd_satn                                            &
     &,points_stress_ud_wake                                            &
     &,points_stress_vd_wake                                            &
     &,points_du_dt_satn                                                &
     &,points_dv_dt_satn                                                &
     &,points_du_dt_wake                                                &
     &,points_dv_dt_wake                                                &
     &,points_u_s_d                                                     &
     &,points_v_s_d                                                     &
     &,points_nsq_s_d                                                   &
     &,points_fr_d                                                      &
     &,points_bld_d                                                     &
     &,points_bldt_d                                                    &
     &,points_num_lim_d                                                 &
     &,points_num_fac_d                                                 &
     &,points_tausx_d                                                   &
     &,points_tausy_d                                                   &
     &,points_taus_scale_d                                              &
     &,points_orog_slope_d                                              &
     &,points_orog_anis_d                                               &
     &,points_orog_dir_d

      INTEGER ::            &
        u_rows,             & ! rows on u grid
        v_rows                ! rows on v grid

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
! Error reporting variables
      CHARACTER (Len=*), PARAMETER ::                                   &
     &  RoutineName='NI_gwd_ctl'

! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches
!-----------------------------------------------------------------------------
! Section 0 - initialisation and setup
!-----------------------------------------------------------------------------   
! Work out number of rows for U and V grids - will depend on whether
! ENDGame or not.

    u_rows = udims%j_end - udims%j_start + 1
    v_rows = vdims%j_end - vdims%j_start + 1
! ----------------------------------------------------------------------
! Section GWD.1 Set stash diagnostic switches
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('NI_GWD_CTL',zhook_in,zhook_handle)

! General case of the atmosphere model, ie : with stash.
      GWSPEC_EFLUX_ON = sf(101,6)
      GWSPEC_EFLUX_P_ON = sf(111,6)
      GWSPEC_SFLUX_ON = sf(102,6)
      GWSPEC_WFLUX_ON = sf(103,6)
      GWSPEC_WFLUX_P_ON = sf(113,6)
      GWSPEC_NFLUX_ON = sf(104,6)
      GWSPEC_EWACC_ON = sf(105,6)
      GWSPEC_EWACC_P_ON = sf(115,6)
      GWSPEC_NSACC_ON = sf(106,6)
      L_t_incr_gwd      = sf(181,6)
      L_u_incr_gwd      = sf(185,6)
      L_v_incr_gwd      = sf(186,6)
      stress_ud_on      = sf(201,6)
      stress_ud_p_on    = sf(241,6)
      stress_vd_on      = sf(202,6)
      du_dt_satn_on     = sf(207,6)
      du_dt_satn_p_on   = sf(247,6)
      dv_dt_satn_on     = sf(208,6)
      u_s_d_on          = sf(214,6)
      v_s_d_on          = sf(215,6)
      nsq_s_d_on        = sf(216,6)
      fr_d_on           = sf(217,6)
      bld_d_on          = sf(218,6)
      bldt_d_on         = sf(222,6)
      stress_ud_satn_on = sf(223,6)
      stress_vd_satn_on = sf(224,6)
      stress_ud_wake_on = sf(227,6)
      stress_vd_wake_on = sf(228,6)
      du_dt_wake_on     = sf(231,6)
      dv_dt_wake_on     = sf(232,6)
      num_lim_d_on      = sf(233,6)
      num_fac_d_on      = sf(234,6)
      tausx_d_on        = sf(235,6)
      tausy_d_on        = sf(236,6)
      taus_scale_d_on   = sf(237,6)
      orog_slope_d_on   = sf(248,6)
      orog_anis_d_on    = sf(249,6)
      orog_dir_d_on     = sf(250,6)


! ----------------------------------------------------------------------
! Section GWD.1
! ----------------------------------------------------------------------

      IF (error_code  ==  0 ) THEN
! Save R_u    before updating
      IF ( l_u_incr_gwd) THEN  ! STASHflag set
        ALLOCATE ( u_incr_diagnostic(udims%i_start:udims%i_end,      & 
                                     udims%j_start:udims%j_end,      & 
                                     udims%k_start:udims%k_end)  )
        DO k=1,model_levels
         DO j=udims%j_start,udims%j_end
          DO i=udims%i_start,udims%i_end
            u_incr_diagnostic(i,j,k) = R_u(i,j,k)
          ENDDO ! i
         ENDDO ! j
        ENDDO ! k
      ELSE
        ALLOCATE ( u_incr_diagnostic(1,1,1) ) 

      END IF                    ! on STASHflag

! Save R_v    before updating
      IF ( l_v_incr_gwd) THEN  ! STASHflag set
        ALLOCATE ( v_incr_diagnostic(vdims%i_start:vdims%i_end,      & 
                                     vdims%j_start:vdims%j_end,      &
                                     vdims%k_start:vdims%k_end)  )
        DO k=1,model_levels
         DO j=vdims%j_start,vdims%j_end
           DO i=vdims%i_start,vdims%i_end
            v_incr_diagnostic(i,j,k) = R_v(i,j,k)
          ENDDO ! i
         ENDDO ! j
        ENDDO ! k
      ELSE
        ALLOCATE ( v_incr_diagnostic(1,1,1) )
      END IF                    ! on STASHflag

! Save T_inc   before updating
      IF ( l_t_incr_gwd) THEN  ! STASHflag set
        ALLOCATE ( t_incr_diagnostic(tdims%i_start:tdims%i_end,      & 
                                     tdims%j_start:tdims%j_end,      &
                                     tkfix1start:tdims%k_end)  )
        DO k=tkfix1start,model_levels
         DO j=tdims%j_start,tdims%j_end
           DO i=tdims%i_start,tdims%i_end
            t_incr_diagnostic(i,j,k) = T_inc(i,j,k)
          END DO ! i
         END DO ! j
        END DO ! k
      ELSE
        ALLOCATE ( t_incr_diagnostic(1,1,1) )
      END IF                    ! on STASHflag

      IF ( stress_ud_on .OR. stress_ud_p_on) THEN  ! STASHflag set
        ALLOCATE ( stress_ud(row_length,u_rows,0:model_levels) )
        points_stress_ud = land_points
      ELSE
        points_stress_ud = 1
        ALLOCATE ( stress_ud(1,1,1) )
      END IF

      IF ( stress_vd_on ) THEN  ! STASHflag set
        ALLOCATE ( stress_vd(row_length,v_rows,0:model_levels) )
        points_stress_vd = land_points
      ELSE
        points_stress_vd = 1
        ALLOCATE ( stress_vd(1,1,1) )
      END IF

      IF ( du_dt_satn_on .OR. du_dt_satn_p_on ) THEN  ! STASHflag set
        ALLOCATE ( du_dt_satn(row_length,u_rows,model_levels) )
        points_du_dt_satn = land_points
      ELSE
        points_du_dt_satn = 1
        ALLOCATE ( du_dt_satn(1,1,1) )
      END IF

      IF ( dv_dt_satn_on ) THEN  ! STASHflag set
        ALLOCATE ( dv_dt_satn(row_length,v_rows,model_levels) )
        points_dv_dt_satn = land_points
      ELSE
        points_dv_dt_satn = 1
        ALLOCATE ( dv_dt_satn(1,1,1) )
      END IF

      IF ( u_s_d_on ) THEN  ! STASHflag set
        ALLOCATE ( u_s_d(row_length,rows) )
        points_u_s_d = land_points
      ELSE
        points_u_s_d = 1
        ALLOCATE ( u_s_d(1,1) )
      END IF

      IF ( v_s_d_on ) THEN  ! STASHflag set
        ALLOCATE ( v_s_d(row_length,rows) )
        points_v_s_d = land_points
      ELSE
        points_v_s_d = 1
        ALLOCATE ( v_s_d(1,1) )
      END IF

      IF ( nsq_s_d_on ) THEN  ! STASHflag set
        ALLOCATE ( nsq_s_d(row_length,rows) )
        points_nsq_s_d = land_points
      ELSE
        points_nsq_s_d = 1
        ALLOCATE ( nsq_s_d(1,1) )
      END IF

      IF ( fr_d_on ) THEN  ! STASHflag set
        ALLOCATE ( fr_d(row_length,rows) )
        points_fr_d = land_points
      ELSE
        points_fr_d = 1
        ALLOCATE ( fr_d(1,1) )
      END IF

      IF ( bld_d_on ) THEN  ! STASHflag set
        ALLOCATE ( bld_d(row_length,rows) )
        points_bld_d = land_points
      ELSE
        points_bld_d = 1
        ALLOCATE ( bld_d(1,1) )
      END IF

      IF ( bldt_d_on ) THEN  ! STASHflag set
        ALLOCATE ( bldt_d(row_length,rows) )
        points_bldt_d = land_points
      ELSE
        points_bldt_d = 1
        ALLOCATE ( bldt_d(1,1) )
      END IF

      IF ( num_lim_d_on ) THEN  ! STASHflag set
        ALLOCATE ( num_lim_d(row_length,rows) )
        points_num_lim_d = land_points
      ELSE
        points_num_lim_d = 1
        ALLOCATE ( num_lim_d(1,1) )
      END IF

      IF ( num_fac_d_on ) THEN  ! STASHflag set
        ALLOCATE ( num_fac_d(row_length,rows) )
        points_num_fac_d = land_points
      ELSE
        points_num_fac_d = 1
        ALLOCATE ( num_fac_d(1,1) )
      END IF

      IF ( stress_ud_satn_on ) THEN  ! STASHflag set
        ALLOCATE ( stress_ud_satn(row_length,u_rows,0:model_levels) )
        points_stress_ud_satn = land_points
      ELSE
        points_stress_ud_satn = 1
        ALLOCATE ( stress_ud_satn(1,1,1) )
      END IF

      IF ( stress_vd_satn_on ) THEN  ! STASHflag set
        ALLOCATE ( stress_vd_satn(row_length,v_rows,0:model_levels) )
        points_stress_vd_satn = land_points
      ELSE
        points_stress_vd_satn = 1
        ALLOCATE ( stress_vd_satn(1,1,1) )
      END IF

      IF ( stress_ud_wake_on ) THEN  ! STASHflag set
        ALLOCATE ( stress_ud_wake(row_length,u_rows,0:model_levels) )
        points_stress_ud_wake = land_points
      ELSE
        points_stress_ud_wake = 1
        ALLOCATE ( stress_ud_wake(1,1,1) )
      END IF

      IF ( stress_vd_wake_on ) THEN  ! STASHflag set
        ALLOCATE ( stress_vd_wake(row_length,v_rows,0:model_levels) )
        points_stress_vd_wake = land_points
      ELSE
        points_stress_vd_wake = 1
        ALLOCATE ( stress_vd_wake(1,1,1) )
      END IF

      IF ( du_dt_wake_on ) THEN  ! STASHflag set
        ALLOCATE ( du_dt_wake(row_length,u_rows,model_levels) )
        points_du_dt_wake = land_points
      ELSE
        points_du_dt_wake = 1
        ALLOCATE ( du_dt_wake(1,1,1) )
      END IF

      IF ( dv_dt_wake_on ) THEN  ! STASHflag set
        ALLOCATE ( dv_dt_wake(row_length,v_rows,model_levels) )
        points_dv_dt_wake = land_points
      ELSE
        points_dv_dt_wake = 1
        ALLOCATE ( dv_dt_wake(1,1,1) )
      END IF

      IF ( tausx_d_on ) THEN  ! STASHflag set
        ALLOCATE ( tausx_d(row_length,u_rows) )
        points_tausx_d = land_points
      ELSE
        points_tausx_d = 1
        ALLOCATE ( tausx_d(1,1) )
      END IF

      IF ( tausy_d_on ) THEN  ! STASHflag set
        ALLOCATE ( tausy_d(row_length,v_rows) )
        points_tausy_d = land_points
      ELSE
        points_tausy_d = 1
        ALLOCATE ( tausy_d(1,1) )
      END IF

      IF ( taus_scale_d_on ) THEN  ! STASHflag set
        ALLOCATE ( taus_scale_d(row_length,rows) )
        points_taus_scale_d = land_points
      ELSE
        points_taus_scale_d = 1
        ALLOCATE ( taus_scale_d(1,1) )
      END IF

      IF ( orog_slope_d_on ) THEN  ! STASHflag set
        ALLOCATE ( orog_slope_d(row_length,rows) )
        points_orog_slope_d = land_points
      ELSE
        points_orog_slope_d = 1
        ALLOCATE ( orog_slope_d(1,1) )
      END IF

      IF ( orog_anis_d_on ) THEN  ! STASHflag set
        ALLOCATE ( orog_anis_d(row_length,rows) )
        points_orog_anis_d = land_points
      ELSE
        points_orog_anis_d = 1
        ALLOCATE ( orog_anis_d(1,1) )
      END IF

      IF ( orog_dir_d_on ) THEN  ! STASHflag set
        ALLOCATE ( orog_dir_d(row_length,rows) )
        points_orog_dir_d = land_points
      ELSE
        points_orog_dir_d = 1
        ALLOCATE ( orog_dir_d(1,1) )
      END IF

      IF ( GWSPEC_EFLUX_ON .OR. GWSPEC_EFLUX_P_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_EFLUX(udims%i_start:udims%i_end,               &
     &              udims%j_start:udims%j_end,tkfix1start:tdims%k_end))
      ELSE
        ALLOCATE (GWSPEC_EFLUX(1,1,1))
      END IF
      IF ( GWSPEC_SFLUX_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_SFLUX(vdims%i_start:vdims%i_end,               &
     &              vdims%j_start:vdims%j_end,tkfix1start:tdims%k_end))
      ELSE
        ALLOCATE (GWSPEC_SFLUX(1,1,1))
      END IF
      IF ( GWSPEC_WFLUX_ON .OR. GWSPEC_WFLUX_P_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_WFLUX(udims%i_start:udims%i_end,               &
     &              udims%j_start:udims%j_end,tkfix1start:tdims%k_end))
      ELSE
        ALLOCATE (GWSPEC_WFLUX(1,1,1))
      END IF
      IF ( GWSPEC_NFLUX_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_NFLUX(vdims%i_start:vdims%i_end,               &
     &              vdims%j_start:vdims%j_end,tkfix1start:tdims%k_end))
      ELSE
        ALLOCATE (GWSPEC_NFLUX(1,1,1))
      END IF
      IF ( GWSPEC_EWACC_ON .OR. GWSPEC_EWACC_P_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_EWACC(udims%i_start:udims%i_end,               &
     &             udims%j_start:udims%j_end,udims%k_start:udims%k_end))
      ELSE
        ALLOCATE (GWSPEC_EWACC(1,1,1))
      END IF
      IF ( GWSPEC_NSACC_ON ) THEN  ! STASHflag set
        ALLOCATE (GWSPEC_NSACC(vdims%i_start:vdims%i_end,               &
     &             vdims%j_start:vdims%j_end,vdims%k_start:vdims%k_end))
      ELSE
        ALLOCATE (GWSPEC_NSACC(1,1,1))
      END IF

!-------------------------------------------------------------------
!Section GWD.2b  Call Orographic Gravity Wave Drag Scheme 
!-------------------------------------------------------------------
      IF (L_gwd) THEN

         SELECT CASE ( i_gwd_vn )
            CASE ( i_gwd_vn_4a )

               CALL G_WAVE_4A(                                                &
               theta_latest, u, v, row_length, rows, n_rows,u_rows,v_rows,    &
               off_x, off_y,                                                  &
               global_row_length,n_proc, n_procy, proc_row_group,at_extremity,&
               model_domain, model_levels,                                    &
               rho, delta_lambda,delta_phi,                                   &
               true_latitude,sd_orog_land,                                    &
               orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
               land_index, land_points, timestep, kay_gwave, gwd_frc,         &
               r_u, r_v, T_inc, l_taus_scale, l_fix_gwsatn, l_gwd_40km,       &
               sat_scheme, gwd_fsat,                                          &
               Gsharp, fbcd, l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,   &
               stress_ud      , stress_ud_on     , stress_ud_p_on       ,     &
                                                   points_stress_ud     ,     &
               stress_vd      , stress_vd_on     , points_stress_vd     ,     &
               stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
               stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
               stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
               stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
               du_dt_satn     , du_dt_satn_on    , du_dt_satn_p_on      ,     &
                                                   points_du_dt_satn    ,     &
               dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
               du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
               dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
               u_s_d          , u_s_d_on         , points_u_s_d         ,     &
               v_s_d          , v_s_d_on         , points_v_s_d         ,     &
               nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
               fr_d           , fr_d_on          , points_fr_d          ,     &
               bld_d          , bld_d_on         , points_bld_d         ,     &
               bldt_d         , bldt_d_on        , points_bldt_d        ,     &
               num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
               num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
               tausx_d        , tausx_d_on       , points_tausx_d       ,     &
               tausy_d        , tausy_d_on       , points_tausy_d       ,     &
               taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
               orog_slope_d   , orog_slope_d_on  , points_orog_slope_d  ,     &
               orog_anis_d    , orog_anis_d_on   , points_orog_anis_d   ,     &
               orog_dir_d     , orog_dir_d_on    , points_orog_dir_d    ,     &
               error_code)
     
            CASE ( i_gwd_vn_5a )

               CALL G_WAVE_5A(                                                &
               theta_latest, u, v, row_length, rows, n_rows,u_rows,v_rows,    &
               off_x, off_y,                                                  &
               global_row_length,n_proc, n_procy, proc_row_group,at_extremity,&
               model_domain, model_levels,                                    &
               rho, delta_lambda,delta_phi,                                   &
               true_latitude,sd_orog_land,                                    &
               orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
               land_index, land_points, timestep, kay_gwave, gwd_frc,         &
               r_u, r_v, T_inc, l_taus_scale, l_fix_gwsatn, l_gwd_40km,       &
               sat_scheme, gwd_fsat,                                          &
               Gsharp, fbcd, l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,   &
               stress_ud      , stress_ud_on     , stress_ud_p_on       ,     &
                                                   points_stress_ud     ,     &
               stress_vd      , stress_vd_on     , points_stress_vd     ,     &
               stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
               stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
               stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
               stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
               du_dt_satn     , du_dt_satn_on    , du_dt_satn_p_on      ,     &
                                                   points_du_dt_satn    ,     &
               dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
               du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
               dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
               u_s_d          , u_s_d_on         , points_u_s_d         ,     &
               v_s_d          , v_s_d_on         , points_v_s_d         ,     &
               nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
               fr_d           , fr_d_on          , points_fr_d          ,     &
               bld_d          , bld_d_on         , points_bld_d         ,     &
               bldt_d         , bldt_d_on        , points_bldt_d        ,     &
               num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
               num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
               tausx_d        , tausx_d_on       , points_tausx_d       ,     &
               tausy_d        , tausy_d_on       , points_tausy_d       ,     &
               taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
               orog_slope_d   , orog_slope_d_on  , points_orog_slope_d  ,     &
               orog_anis_d    , orog_anis_d_on   , points_orog_anis_d   ,     &
               orog_dir_d     , orog_dir_d_on    , points_orog_dir_d    ,     &
               error_code)

            CASE DEFAULT ! i_gwd_vn

               errorstatus = 10
               WRITE (cmessage,'(A,A,I6)') 'Gravity wave drag version ',     &
               'not recognised. i_gwd_vn = ',i_gwd_vn
               CALL Ereport ( RoutineName, errorstatus, cmessage)

         END SELECT ! i_gwd_vn

      END IF ! l_gwd

      IF ( L_use_ussp ) THEN
         CALL GW_USSP(MODEL_LEVELS, MODEL_DOMAIN, ROWS, N_ROWS,               &
                OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                     &
                global_row_length,n_proc,n_procy,proc_row_group,at_extremity, &
                R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,             &
                R_U,R_V, T_inc,                                               &
                SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                      &
                THETA_LATEST, RHO, TIMESTEP, U, V,                            &
                l_gw_heating(3),                                              &
                GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,          &
                GWSPEC_EWACC,GWSPEC_NSACC,                                    &
                GWSPEC_EFLUX_ON,GWSPEC_EFLUX_P_ON,GWSPEC_SFLUX_ON,            &
                GWSPEC_WFLUX_ON,GWSPEC_WFLUX_P_ON,GWSPEC_NFLUX_ON,            &
                GWSPEC_EWACC_ON,GWSPEC_EWACC_P_ON,GWSPEC_NSACC_ON)
      END IF ! l_use_ussp

! ----------------------------------------------------------------------
! Section GWD.3 Call GWD diagnostics
! ----------------------------------------------------------------------
        IF(sf(0,6)) THEN ! diagnostics requested this timestep
! DEPENDS ON: diagnostics_gwd
          CALL diagnostics_gwd(                                         &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, u_rows, v_rows                     &
     &,                      off_x, off_y                               &
     &,                      at_extremity                               &
     &,                      u, v, R_u, R_v, T_inc                      &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      t_incr_diagnostic                          &
     &,                      exner_theta_levels                         &
     &,                      stress_ud     ,  stress_vd                 &
     &,                      stress_ud_satn,  stress_vd_satn            &
     &,                      stress_ud_wake,  stress_vd_wake            &
     &,                      du_dt_satn    ,  dv_dt_satn                &
     &,                      du_dt_wake   ,  dv_dt_wake                 &
     &,                      u_s_d, v_s_d, nsq_s_d                      &
     &,                      num_lim_d, num_fac_d                       &
     &,                      fr_d, bld_d, bldt_d                        &
     &, tausx_d, tausy_d, taus_scale_d                                  &
     &, sd_orog_land, land_sea_mask, land_points                        &
     &, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
     &, orog_slope_d, orog_anis_d, orog_dir_d                           &      
     &, GWSPEC_EFLUX, GWSPEC_SFLUX, GWSPEC_WFLUX                        &
     &, GWSPEC_NFLUX, GWSPEC_EWACC, GWSPEC_NSACC,                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork                                                        &
     & )
       END IF            ! on sf(0,6) 

      IF ( ALLOCATED( u_incr_diagnostic ) ) DEALLOCATE ( u_incr_diagnostic )
      IF ( ALLOCATED( v_incr_diagnostic ) ) DEALLOCATE ( v_incr_diagnostic )
      IF ( ALLOCATED( t_incr_diagnostic ) ) DEALLOCATE ( t_incr_diagnostic )
      IF ( ALLOCATED( stress_ud         ) ) DEALLOCATE ( stress_ud         )
      IF ( ALLOCATED( stress_vd         ) ) DEALLOCATE ( stress_vd         )
      IF ( ALLOCATED( du_dt_satn        ) ) DEALLOCATE ( du_dt_satn        )
      IF ( ALLOCATED( dv_dt_satn        ) ) DEALLOCATE ( dv_dt_satn        )
      IF ( ALLOCATED( u_s_d             ) ) DEALLOCATE ( u_s_d             )
      IF ( ALLOCATED( v_s_d             ) ) DEALLOCATE ( v_s_d             )
      IF ( ALLOCATED( nsq_s_d           ) ) DEALLOCATE ( nsq_s_d           )
      IF ( ALLOCATED( fr_d              ) ) DEALLOCATE ( fr_d              )
      IF ( ALLOCATED( bld_d             ) ) DEALLOCATE ( bld_d             )
      IF ( ALLOCATED( bldt_d            ) ) DEALLOCATE ( bldt_d            )
      IF ( ALLOCATED( num_lim_d         ) ) DEALLOCATE ( num_lim_d         )
      IF ( ALLOCATED( num_fac_d         ) ) DEALLOCATE ( num_fac_d         )
      IF ( ALLOCATED( stress_ud_satn    ) ) DEALLOCATE ( stress_ud_satn    )
      IF ( ALLOCATED( stress_vd_satn    ) ) DEALLOCATE ( stress_vd_satn    )
      IF ( ALLOCATED( stress_ud_wake    ) ) DEALLOCATE ( stress_ud_wake    )
      IF ( ALLOCATED( stress_vd_wake    ) ) DEALLOCATE ( stress_vd_wake    )
      IF ( ALLOCATED( du_dt_wake        ) ) DEALLOCATE ( du_dt_wake        )
      IF ( ALLOCATED( dv_dt_wake        ) ) DEALLOCATE ( dv_dt_wake        )
      IF ( ALLOCATED( tausx_d           ) ) DEALLOCATE ( tausx_d           )
      IF ( ALLOCATED( tausy_d           ) ) DEALLOCATE ( tausy_d           )
      IF ( ALLOCATED( taus_scale_d      ) ) DEALLOCATE ( taus_scale_d      )
      IF ( ALLOCATED( orog_slope_d      ) ) DEALLOCATE ( orog_slope_d      )
      IF ( ALLOCATED( orog_anis_d       ) ) DEALLOCATE ( orog_anis_d       )
      IF ( ALLOCATED( orog_dir_d        ) ) DEALLOCATE ( orog_dir_d        )
      IF ( ALLOCATED( GWSPEC_EFLUX      ) ) DEALLOCATE ( GWSPEC_EFLUX      )
      IF ( ALLOCATED( GWSPEC_SFLUX      ) ) DEALLOCATE ( GWSPEC_SFLUX      )
      IF ( ALLOCATED( GWSPEC_WFLUX      ) ) DEALLOCATE ( GWSPEC_WFLUX      )
      IF ( ALLOCATED( GWSPEC_NFLUX      ) ) DEALLOCATE ( GWSPEC_NFLUX      )
      IF ( ALLOCATED( GWSPEC_EWACC      ) ) DEALLOCATE ( GWSPEC_EWACC      )
      IF ( ALLOCATED( GWSPEC_NSACC      ) ) DEALLOCATE ( GWSPEC_NSACC      )

      END IF ! on error code equal to zero

      IF (lhook) CALL dr_hook('NI_GWD_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_gwd_ctl

END MODULE NI_gwd_ctl_mod
