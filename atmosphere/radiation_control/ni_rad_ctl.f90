! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface between Atmos_physics1 and the radiation code
!
! Purpose:
!   This glue routine has been inserted in order to
!   manage the 3Z version of the radiation code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
      Subroutine NI_RAD_CTL (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc            &
     &, n_procx, n_procy, neighbour, g_rows, g_row_length               &
     &, me, global_cloud_top                                            &

! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_model_levels, bl_levels                       &
     &, Ozone_levels, cloud_levels, N_cca_levels                        &
     &, NTILES, LAND_FIELD, nice_use, DUST_DIM1, DUST_DIM2              &
     &, biogenic_dim1                                                   &
     &, biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2       &
     &, bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3         &
     &, co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2        &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &
     &, ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2                &
     &, ukca_dim1, ukca_dim2, n_ukca_mode, n_ukca_cpnt                  &

! Model switches
     &, model_domain                                                    &
     &, L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                    &
     &, L_CAL360                                                        &
     &, L_emcorr                                                        &
     &, Ltimer, L_ssice_albedo                                          &
     &, L_snow_albedo                                                   &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a,l_cice_alb   &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_co2_interactive                                               &
     &, L_ukca                                                          &
     &, L_USE_ARCL                                                      &
     &, L_DUST, l_sulpc_so2, l_soot, l_biomass, l_ocff, l_nitrate       &
      , L_cosp                                                          &

! model Parameters
     &, min_trop_level, max_trop_level                                  &
     &, Ntot_land, Ntot_sea                                             &

! in coordinate information
     &, rho_r2                                                          &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
     &, trindx                                                          &

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag, radiation_tstep_prog                      &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &
     &,         PREVIOUS_TIME                                           &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork1,                                                      &
     & STASHwork2                                                       &

! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs,L_SCMDiags                                            &

! in data fields.
     &, p_star                                                          &
     &, p_layer_boundaries, p_layer_centres                             &
     &, p, p_theta_levels                                               &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, fland, land0p5                                   &
     &, T_SURF,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE_CAT,AREA_CLOUD_FRACTION  &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, biogenic, so4_aitken, so4_accu, so4_diss                        &
     &, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
     &, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss, aerosol, arcl&
     &, ukca_radaer, sea_salt_film, sea_salt_jet, co2_3D                &
     &, cos_zenith_angle                                                &
     &, can_rad_mod                                                     &
      , frac_control, n_drop_pot                                        &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, snow_depth, snow_depth_sea_cat, ice_fract, ice_fract_cat        &
     &, ice_thick_cat, rgrain, soot                                     &
     &, cca, ccb, cct, cclwp, ccw, lcbase                               &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height                                   &
     &, T_trop_level, T_trop_height, zh                                 &
      , land_index, albsoil, albobs_sw, albobs_vis, albobs_nir, lai     &
      , snow_tile, frac, tstar_tile, z0_tile                            &
     &, dOLR_rts, LW_down, SW_tile_rts                                  &
     &, u_1, v_1                                                        &
     &, land_alb,sice_alb                                               &
     &, ES_SPACE_INTERP, RAD_MASK                                       &

! in/out
     &, T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                      &
     &, qcf2_n, qrain_n, qgraup_n                                       &
     &, T_inc, q_inc, qcl_inc, cf_inc, cfl_inc                          &
     &, sum_eng_fluxes                                                  &

! out.
      , photosynth_act_rad, rad_hr, dOLR, SW_tile                       &

! COSP arguments
      , cosp_gbx                                                        &
! error information
     &, Error_code  )

      USE rad_input_mod
      USE ukca_option_mod, ONLY: l_ukca_radaer
      USE ukca_radaer_struct_mod
      USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
       
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE cosp_types_mod, ONLY: cosp_gridbox
      USE Submodel_Mod
      USE nstypes
      USE switches, ONLY: l_mod_barker_albedo, l_ctile
      IMPLICIT NONE

! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
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
     &, global_rows                                                     &
                           ! number of global rows
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, me         ! My processor number

!     Global topmost cloudy level
      INTEGER, INTENT(IN) :: global_cloud_top

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid
!
! Model dimensions
!
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, Ozone_levels                                                    &
     &, cloud_levels                                                    &
     &, N_cca_levels                                                    &
     &, ntiles                                                          &
     &, land_field                                                      &
     &, nice_use                                                        &
                  ! number of sea ice categories used in radiation
     &, DUST_DIM1                                                       &
                  !dimensions of mineral dust arrays
     &, DUST_DIM2                                                       &
                  !dimensions of biogenic aerosol arrays
     &, biogenic_dim1                                                   &
     &, biogenic_dim2                                                   &
     &, sulp_dim1                                                       &
                                ! dimensions of S Cyc arrays
     &, sulp_dim2                                                       &
     &, soot_dim1                                                       &
                                ! dimensions of soot arrays
     &, soot_dim2                                                       &
     &, bmass_dim1                                                      &
                                ! dimensions of biomass arrays
     &, bmass_dim2                                                      &
     &, ocff_dim1                                                       &
     &, ocff_dim2                                                       &
                                ! dimensions of fossil-fuel arrays
     &, salt_dim1                                                       &
                                ! dimensions of sea-salt arrays
     &, salt_dim2                                                       &
     &, salt_dim3                                                       &
     &, co2_dim2                                                        &
                               ! dimensions of CO2 array
     &, co2_dim_len                                                     &
     &, co2_dim_row                                                     &
                               ! dimensions of aerosol clim for NWP
     &, arcl_dim1                                                       &
     &, arcl_dim2                                                       &
                               ! dimensions of nitrate arrays
     &, nitrate_dim1                                                    &
     &, nitrate_dim2                                                    &
                               ! dimensions of UKCA_RADAER arrays
     &, ukca_dim1                                                       &
     &, ukca_dim2                                                       &
     &, n_ukca_mode                                                     &
     &, n_ukca_cpnt
     
!
! Model switches
!
      Integer                                                           &
     &  model_domain

      Logical                                                           &
     &  L_Rad_Step                                                      &
                             ! true on radiation timestep (3A)
     &, L_Rad_Step_diag                                                 &
                             ! true if fast radiation timestep (3C)
     &, L_Rad_Step_prog                                                 &
                             ! true if slow radiation timestep (3C)
     &, L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_DUST                                                          &
                   ! mineral dust available for use
     &, L_sulpc_so2                                                     &
                   ! Sulphur cycle available for use
     &, L_soot                                                          &
                   ! soot is available for use
     &, L_biomass                                                       &
                   ! Biomass smoke is available for use
     &, L_ocff                                                          &
                   ! OCFF is available for use
     &, L_nitrate                                                       &
                   ! nitrate aerosol is available for use
     &, L_emcorr                                                        &
                    ! true if energy correction scheme is to be used.
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_ssice_albedo                                                  &
                       ! Switch on the effect of snow on sea-ice albedo
     &, L_sice_meltponds                                                &
                         ! true if sea ice meltponds albedo scheme
     &, L_sice_scattering                                               &
                          ! true if sea ice internal scattering scheme
     &, L_sice_hadgem1a                                                 &
                        ! true if HadGEM1 albedo bug corrected
     &, l_cice_alb                                                      &
                       ! true if sea ice CICE albedo scheme
     &, L_snow_albedo                                                   &
                       ! True if spectral albedo scheme selected
     &, l_mcr_qcf2                                                      &
                       ! Use second ice category
     &, l_mcr_qrain                                                     &
                       ! Use prognostic rain
     &, l_mcr_qgraup                                                    &
                       ! Use graupel
     &, L_pc2                                                           &
                       ! True if using the PC2 cloud scheme
     &, L_mixing_ratio                                                  &
                       ! True if using mixing ratios
     &, L_MURK_RAD                                                      &
                           ! True if using radiative effects of 'murk'.
     &, L_co2_interactive                                               &
                           ! Controls the use of 3D CO2 field
     &, L_ukca                                                          
                           ! Switch for UKCA sub-model

!     Switch for COSP
      LOGICAL :: L_cosp

! model parameters
      Real                                                              &
     &  timestep                                                        &
     &, radiation_timestep                                              &
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                                           

      Integer                                                           &
     &  min_trop_level                                                  &
                        ! Lowest permitted level for the tropopause
!                       ! used for radiative purposes.
     &, max_trop_level                                                  &
                        ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.
     &, can_rad_mod     ! which canopy light mod are we using using Q                                 

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
!     Missing number indicators

! Diagnostics info
      Real                                                              &
     & STASHwork1(*)                                                    &
                         ! STASH workspace
     &,STASHwork2(*)     ! STASH workspace


! Data arrays
      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, p_star(row_length, rows)                                        &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y, model_levels)                &
     &, p(1-off_x:row_length+off_x,                                     &
     &    1-off_y:rows+off_y, model_levels)                             &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                     1-off_y:rows+off_y, model_levels)

      Real                                                              &
     &  cca (row_length, rows, N_cca_levels)                            &
     &, cclwp(row_length, rows)                                         &
                                ! condensed water path (KG/M**2)
     &, area_cloud_fraction(row_length, rows, wet_model_levels)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)                                        &
     &, TSTAR_LAND(ROW_LENGTH,ROWS)                                     &
                                    ! IN Land mean sfc temperature (K)
     &, TSTAR_SEA(ROW_LENGTH,ROWS)                                      &
                                    ! IN Open sea sfc temperature (K).
     &, TSTAR_SICE_CAT(ROW_LENGTH,ROWS,NICE_USE) 
                                    ! IN Sea-ice sfc temperature (K).

      REAL ::                                                           &
        ice_fract_cat(row_length, rows, nice_use),                      &
!         Area fraction of sea ice categories
        ice_thick_cat(row_length, rows, nice_use)
!         Effective thickness of each sea ice categories


      logical                                                           &
     &  land_sea_mask(row_length, rows)                                 &
     &, land0p5(row_length, rows)                                       &
!         A mask set to .TRUE. if the fraction of land in the grid-box
!         exceeds 0.5.
     &, RAD_MASK(row_length, rows)
!         A mask which ensures a chequerboard pattern of radiation
!         calculations over the whole domain (not just one PE)


      Real                                                              &
     &  fland(land_field)
!         Fractional amount of land at each land point

      Real                                                              &
     &  Ntot_land                                                       &
                             ! Number of droplets over land / m-3
     &, Ntot_sea             ! Number of droplets over sea / m-3

      Real                                                              &
     &  snow_depth (row_length, rows)                                   &
                                      ! snow/qrclim.snow.(month)
     &, snow_depth_sea_cat (row_length, rows, nice_use) 
                                      !  snow depth on sea ice

      Real                                                              &
     &  ice_fract (row_length, rows)                                    &
                                     ! ice/qrclim.ice.(month)
     &, soot (row_length, rows)

      Real                                                              &
     &  albsoil(land_field)                                             &
      , albobs_sw(land_field)                                           &
      , albobs_vis(land_field)                                          &
      , albobs_nir(land_field)                                          & 
     &, lai(land_field, npft)                                           &
     &, rgrain(land_field, ntiles)                                      &
     &, snow_tile(land_field, ntiles)                                   &
     &, frac(land_field, ntype)                                         &
     &, frac_control(land_field, ntype)                                 &
     &, tstar_tile(land_field, ntiles)                                  &
     &, z0_tile(land_field, ntiles)                                     &
     &, dOLR_rts(row_length, rows)                                      &
                                        ! TOA - surface upward LW
     &, LW_down(row_length, rows)                                       &
                                        ! Surface downward LW
     &, SW_tile_rts(land_field, ntiles)                                 &
                                        ! Surface net SW on land tiles
     &, land_alb(row_length, rows)                                      &
                                        ! Mean land albedo
     &, sice_alb(row_length, rows)                                      &
                                        ! Mean sea-ice albedo
     &, ES_SPACE_INTERP(4, row_length, rows)
!              Coeffs for spatial interpolation of radiation quantities
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Model switch for each species
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Array indices of components
      Integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      ! Mass-mixing ratios
      Real                                                              &
     &  arcl(row_length, rows, model_levels, n_arcl_compnts)

! UKCA_RADAER structure: Interaction between UKCA-MODE aerosols
!                        and radiation
      TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer

      Integer                                                           &
     &  land_index(land_field)

      Real                                                              &
     &  ozone(row_length, rows, ozone_levels)                           &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)                                  &
     &, u_1(salt_dim1, salt_dim2)                                       &
     &, v_1(salt_dim1, salt_dim2)                                       &
     &, SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, dirpar_inc(row_length, rows)                                    &
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, co2_3D(row_length, rows, model_levels)                          &
      , cos_zenith_angle(row_length, rows)                              &
      , n_drop_pot(row_length, rows, model_levels)
! chemical greenhouse gas fields
      INTEGER, INTENT(IN) :: ngrgas
      REAL, INTENT(IN) ::                                               &
     &        grgas_field(row_length,rows,model_levels,ngrgas)

! Co-ordinate arrays
      Real                                                              &
     &  rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          model_levels)                                           &
                               ! Density*radius^2
     &, delta_lambda                                                    &
     &, delta_phi

! arguments with intent in. ie: input variables.


      Real   , intent(in) :: ccw   (row_length, rows, wet_model_levels)
      Integer, intent(in) :: lcbase(row_length, rows)


! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number                                                 &
     &,         PREVIOUS_TIME(7)

      Real                                                              &
     &  recip_timestep

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

!     Level of tropopause
      INTEGER, INTENT(IN) :: trindx(row_length, rows)

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs          ! No of diagnostics packages

      Logical                                                           &
                           ! Logicals for diagnostics packages
     &  L_SCMDiags(nSCMDpkgs)

      Integer                                                           &
     &  Error_code

! Variables with intent (in/out)

      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, qcf2_n(row_length, rows, wet_model_levels)                      &
                                                     ! 2nd ice prog
     &, qrain_n(row_length, rows, wet_model_levels)                     &
                                                     ! Rain prognostic
     &, qgraup_n(row_length, rows, wet_model_levels)                    &
                                                     ! Graupel
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, cff_n(row_length, rows, wet_model_levels)                       &
     &, T_inc(row_length, rows, model_levels)                           &
     &, q_inc(row_length, rows, wet_model_levels)                       &
     &, qcl_inc(row_length, rows, wet_model_levels)                     &
     &, cf_inc(row_length, rows, wet_model_levels)                      &
     &, cfl_inc(row_length, rows, wet_model_levels)                     &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
     &, sum_eng_fluxes(row_length, rows)

      Real, Intent(InOut) ::                                            &
     &  so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &                                               model_levels)      &
     &, so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, biogenic(row_length, rows, model_levels)                        &
     &, aerosol(row_length, rows, model_levels)                         
!
! Variables with intent out
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, rad_hr(row_length, rows, bl_levels, 2)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
     &, dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
     &, SW_tile(land_field, ntiles) ! Surface net SW on land tiles
! COSP input variables variables
      TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx

! Temporary variables used for running radiation in diagnostic mode:
      REAL, ALLOCATABLE :: SW_incs_tmp(:, :, :)
      REAL, ALLOCATABLE :: LW_incs_tmp(:, :, :)
      REAL, ALLOCATABLE :: dirpar_inc_tmp (:, :)
      REAL, ALLOCATABLE :: dOLR_rts_tmp(:, :)
      REAL, ALLOCATABLE :: LW_down_tmp(:, :)
      REAL, ALLOCATABLE :: SW_tile_rts_tmp(:, :)
      REAL, ALLOCATABLE :: land_alb_tmp(:, :)
      REAL, ALLOCATABLE :: sice_alb_tmp(:, :)

      REAL, ALLOCATABLE :: T_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: q_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: qcl_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: qcf_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: qcf2_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: qrain_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: qgraup_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: cf_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: cfl_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: cff_n_tmp(:, :, :)
      REAL, ALLOCATABLE :: T_inc_tmp(:, :, :)
      REAL, ALLOCATABLE :: q_inc_tmp(:, :, :)
      REAL, ALLOCATABLE :: qcl_inc_tmp(:, :, :)
      REAL, ALLOCATABLE :: cf_inc_tmp(:, :, :)
      REAL, ALLOCATABLE :: cfl_inc_tmp(:, :, :)
      REAL, ALLOCATABLE :: sum_eng_fluxes_tmp(:, :)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('NI_RAD_CTL',zhook_in,zhook_handle)
      IF (lrad_diag_mode) THEN

        ALLOCATE(SW_incs_tmp(row_length, rows, 0:model_levels+1)  )
        ALLOCATE(LW_incs_tmp(row_length, rows, 0:model_levels)    )
        ALLOCATE(dirpar_inc_tmp (row_length, rows)                )
        ALLOCATE(dOLR_rts_tmp(row_length, rows)                   )
        ALLOCATE(LW_down_tmp(row_length, rows)                    )
        ALLOCATE(SW_tile_rts_tmp(land_field, ntiles)              )
        ALLOCATE(land_alb_tmp(row_length, rows)                   )
        ALLOCATE(sice_alb_tmp(row_length, rows)                   )

        ALLOCATE(T_n_tmp(row_length, rows, model_levels)          )
        ALLOCATE(q_n_tmp(row_length, rows, wet_model_levels)      )
        ALLOCATE(qcl_n_tmp(row_length, rows, wet_model_levels)    )
        ALLOCATE(qcf_n_tmp(row_length, rows, wet_model_levels)    )
        IF (l_mcr_qcf2) THEN
          ALLOCATE(qcf2_n_tmp(row_length, rows, wet_model_levels) )
        ENDIF
        IF (l_mcr_qrain) THEN
          ALLOCATE(qrain_n_tmp(row_length, rows, wet_model_levels))
        ENDIF
        IF (l_mcr_qgraup) THEN
          ALLOCATE(qgraup_n_tmp(row_length, rows, wet_model_levels))
        ENDIF
        ALLOCATE(cf_n_tmp(row_length, rows, wet_model_levels)     )
        ALLOCATE(cfl_n_tmp(row_length, rows, wet_model_levels)    )
        ALLOCATE(cff_n_tmp(row_length, rows, wet_model_levels)    )
        ALLOCATE(T_inc_tmp(row_length, rows, model_levels)        )
        ALLOCATE(q_inc_tmp(row_length, rows, wet_model_levels)    )
        ALLOCATE(qcl_inc_tmp(row_length, rows, wet_model_levels)  )
        ALLOCATE(cf_inc_tmp(row_length, rows, wet_model_levels)   )
        ALLOCATE(cfl_inc_tmp(row_length, rows, wet_model_levels)  )
        ALLOCATE(sum_eng_fluxes_tmp(row_length, rows)             )

! For diagnostic mode, save a copy of the Intent(in/out) variables.

        SW_incs_tmp = SW_incs
        LW_incs_tmp = LW_incs
        dirpar_inc_tmp = dirpar_inc
        dOLR_rts_tmp = dOLR_rts
        LW_down_tmp = LW_down
        SW_tile_rts_tmp = SW_tile_rts
        land_alb_tmp = land_alb
        sice_alb_tmp = sice_alb

        T_n_tmp = T_n
        q_n_tmp = q_n
        qcl_n_tmp = qcl_n
        qcf_n_tmp = qcf_n
        cf_n_tmp = cf_n
        cfl_n_tmp = cfl_n
        cff_n_tmp = cff_n
        IF (l_mcr_qcf2) THEN
          qcf2_n_tmp = qcf2_n
        ENDIF
        IF (l_mcr_qrain) THEN
          qrain_n_tmp = qrain_n
        ENDIF
        IF (l_mcr_qgraup) THEN
          qgraup_n_tmp = qgraup_n
        ENDIF
        T_inc_tmp = T_inc
        q_inc_tmp = q_inc
        qcl_inc_tmp = qcl_inc
        cf_inc_tmp = cf_inc
        cfl_inc_tmp = cfl_inc
        sum_eng_fluxes_tmp = sum_eng_fluxes

      END IF

! DEPENDS ON: glue_rad
       CALL Glue_Rad (                                                  &
! Parallel variables
     &    halo_i, halo_j, off_x, off_y, global_row_length, global_rows  &
     &,   proc_row_group, proc_col_group, at_extremity, n_proc, n_procx &
     &,   n_procy, neighbour, g_rows, g_row_length, me                  &
     &,   global_cloud_top                                              &

! model dimensions.
     &,   row_length, rows, n_rows                                      &
     &,   model_levels, wet_model_levels, bl_levels                     &
     &,   Ozone_levels, cloud_levels, N_cca_levels                      &
     &,   NTILES, LAND_FIELD, nice_use, DUST_DIM1, DUST_DIM2            &
     &,   biogenic_dim1                                                 &
     &,   biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2     &
     &,   bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3       &
     &,   co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2      &
     &,   n_arcl_species, n_arcl_compnts, i_arcl_compnts                &
     &,   ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2              &
     &,   ukca_dim1, ukca_dim2, n_ukca_mode, n_ukca_cpnt                &
     
! Model switches
     &,   model_domain                                                  &
     &,   L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                  &
     &,   L_Forcing, L_Timestep, L_Radiance                             &
     &,   L_CAL360, L_SEC_VAR, L_EqT                                    &
      ,   L_INHOM_CLOUD, L_DUST, L_USE_DUST, L_USE_BIOGENIC             &
     &,   l_sulpc_so2, L_use_sulpc_direct                               &
     &,   l_soot, L_use_soot_direct, L_use_soot_indirect                &
     &,   l_biomass, L_use_bmass_direct, L_use_bmass_indirect           &
     &,   l_ocff, L_use_ocff_direct, L_use_ocff_indirect                &
     &,   L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW              &
     &,   l_nitrate, L_use_nitrate_direct, L_use_nitrate_indirect       &
     &,   L_emcorr, L_climat_aerosol, L_clim_aero_hgt                   &
     &,   L_HadGEM1_Clim_Aero, Ltimer, L_ssice_albedo                   &
     &,   L_snow_albedo                                                 &
     &,   L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a,l_cice_alb &
     &,   L_use_seasalt_indirect                                        &
     &,   L_use_seasalt_direct                                          &
     &,   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                         &
     &,   L_pc2, l_mixing_ratio                                         &
     &,   L_MURK_RAD                                                    &
     &,   L_rad_deg                                                     &
     &,   L_co2_interactive                                             &
     &,   L_ukca, l_ukca_radaer                                         &
     &,   L_USE_ARCL, L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX &
      ,   L_cosp                                                        &
     &,   CAN_RAD_MOD                                                   &

! model Parameters
     &,   A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size    &
     &,   INHOM_CLOUD_SW, INHOM_CLOUD_LW, DP_CORR_STRAT, DP_CORR_CONV   &
     &,   CO2_MMR, O2MMR, N2OMMR, CH4MMR                                &
     &,   C11MMR, C12MMR                                                &
     &,   C113MMR, C114MMR, HCFC22MMR, HFC125MMR, HFC134AMMR            &
     &,   alpham, alphac, alphab, dtice                                 &
     &,   dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta                    &
     &,   min_trop_level, max_trop_level                                &
     &,   Ntot_land, Ntot_sea, aero_bl_levels                           &

! in coordinate information
     &,   rho_r2                                                        &
     &,   delta_lambda, delta_phi                                       &
     &,   lat_rot_NP, long_rot_NP                                       &
     &,   trindx                                                        &

! in time stepping information.

     &,   timestep, radiation_timestep                                  &
     &,   radiation_tstep_diag, radiation_tstep_prog                    &
     &,   val_year, val_day_number, val_hour, val_minute                &
     &,   val_second, timestep_number                                   &
     &,   PREVIOUS_TIME                                                 &

! diagnostic info
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ARGSTS end
     &   STASHwork1                                                     &
     &,  STASHwork2                                                     &
!
! SCM diagnostics switches (dummy in full um)
     &,   nSCMDpkgs,L_SCMDiags                                          &

! in data fields.
     &,   p_star                                                        &
     &,   p_layer_boundaries, p_layer_centres                           &
     &,   p, p_theta_levels                                             &
     &,   exner_rho_levels, exner_theta_levels                          &
     &,   land_sea_mask, fland,land0p5,l_ctile                          &
      ,   T_SURF,TSTAR_SEA,TSTAR_SICE_CAT,AREA_CLOUD_FRACTION           &
     &,   DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6   &
     &,   biogenic, so4_aitken, so4_accu, so4_diss                      &
     &,   soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld &
     &,   ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss             &
     &,   aerosol, arcl                                                 &
     &,   ukca_radaer, sea_salt_film, sea_salt_jet, co2_3D              &
      ,   frac_control, n_drop_pot                                      &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,   snow_depth, snow_depth_sea_cat, ice_fract, ice_fract_cat      &
     &,   ice_thick_cat, rgrain, soot                                   &
     &,   cca, ccb, cct, cclwp, ccw, lcbase                             &
     &,   ozone, SW_incs, LW_incs,dirpar_inc                            &
     &,   O3_trop_level, O3_trop_height                                 &
     &,   T_trop_level, T_trop_height, zh                               &
      ,   land_index, albsoil, albobs_sw, albobs_vis, albobs_nir, lai   &
      ,   snow_tile, frac, tstar_tile, z0_tile                          &
     &,   dOLR_rts, LW_down, SW_tile_rts                                &
     &,   land_alb,sice_alb                                             &
     &,   ES_SPACE_INTERP, A_SW_radstep, A_LW_radstep                   &
     &,   A_SW_radstep_diag, A_SW_radstep_prog                          &
     &,   A_LW_radstep_diag, A_LW_radstep_prog, RAD_MASK                &

! in/out
     &,   T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                    &
     &,   qcf2_n, qrain_n, qgraup_n                                     &
     &,   T_inc, q_inc, qcl_inc, cf_inc, cfl_inc                        &
     &,   sum_eng_fluxes                                                &
     &,   cos_zenith_angle                                              &

! out.
      ,   photosynth_act_rad, rad_hr, dOLR, SW_tile                     &

! COSP arguments
      ,   cosp_gbx                                                      &
! error information
     &,   Error_code  )


      IF (lrad_diag_mode) THEN

! For diagnostic mode, the Intent(out) variables are set to zero.

        photosynth_act_rad = 0.0
        rad_hr = 0.0
        dOLR = 0.0
        SW_tile = 0.0

! The Intent(in/out) variables are reset to their original values.

        SW_incs = SW_incs_tmp
        LW_incs = LW_incs_tmp
        dirpar_inc = dirpar_inc_tmp
        dOLR_rts = dOLR_rts_tmp
        LW_down = LW_down_tmp
        SW_tile_rts = SW_tile_rts_tmp
        land_alb = land_alb_tmp
        sice_alb = sice_alb_tmp

        T_n = T_n_tmp
        q_n = q_n_tmp
        qcl_n = qcl_n_tmp
        qcf_n = qcf_n_tmp
        cf_n = cf_n_tmp
        cfl_n = cfl_n_tmp
        cff_n = cff_n_tmp
        IF (l_mcr_qcf2) THEN
          qcf2_n = qcf2_n_tmp
        ENDIF
        IF (l_mcr_qrain) THEN
          qrain_n = qrain_n_tmp
        ENDIF
        IF (l_mcr_qgraup) THEN
          qgraup_n = qgraup_n_tmp
        ENDIF
        T_inc = T_inc_tmp
        q_inc = q_inc_tmp
        qcl_inc = qcl_inc_tmp
        cf_inc = cf_inc_tmp
        cfl_inc = cfl_inc_tmp
        sum_eng_fluxes = sum_eng_fluxes_tmp

        DEALLOCATE(SW_incs_tmp)
        DEALLOCATE(LW_incs_tmp)
        DEALLOCATE(dirpar_inc_tmp)
        DEALLOCATE(dOLR_rts_tmp)
        DEALLOCATE(LW_down_tmp)
        DEALLOCATE(SW_tile_rts_tmp)
        DEALLOCATE(land_alb_tmp)
        DEALLOCATE(sice_alb_tmp)

        DEALLOCATE(T_n_tmp)
        DEALLOCATE(q_n_tmp)
        DEALLOCATE(qcl_n_tmp)
        DEALLOCATE(qcf_n_tmp)
        IF (l_mcr_qcf2) THEN
          DEALLOCATE(qcf2_n_tmp)
        ENDIF
        IF (l_mcr_qrain) THEN
          DEALLOCATE(qrain_n_tmp)
        ENDIF
        IF (l_mcr_qgraup) THEN
          DEALLOCATE(qgraup_n_tmp)
        ENDIF
        DEALLOCATE(cf_n_tmp)
        DEALLOCATE(cfl_n_tmp)
        DEALLOCATE(cff_n_tmp)
        DEALLOCATE(T_inc_tmp)
        DEALLOCATE(q_inc_tmp)
        DEALLOCATE(qcl_inc_tmp)
        DEALLOCATE(cf_inc_tmp)
        DEALLOCATE(cfl_inc_tmp)
        DEALLOCATE(sum_eng_fluxes_tmp)

      END IF

      IF (lhook) CALL dr_hook('NI_RAD_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_RAD_CTL
