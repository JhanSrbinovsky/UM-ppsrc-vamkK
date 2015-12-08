! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for radiation.

! Description:
!   Module containing input switches/settings as used by the radiation code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Method:
!   Switches are initialised to false and read in from the
!   UMUI. The module may then be used directly where the switches
!   are needed within the radiation code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE rad_input_mod

USE rad_com_mod

IMPLICIT NONE

! ----------------
! Control options
! ----------------

LOGICAL :: l_radiation = .FALSE.    !  F: Turns off radiation code

LOGICAL :: l_use_dust = .FALSE.     !  Use mineral dust in rad calculations
LOGICAL :: l_use_biogenic = .FALSE. !  Use biogenic aerosol in radiation code
      
! Use SO4 aerosol from sulphur cycle for direct/indirect effect
! in radiation, the latter for both SW and LW.
LOGICAL :: l_use_sulpc_direct = .FALSE.
LOGICAL :: l_use_sulpc_indirect_sw = .FALSE.
LOGICAL :: l_use_sulpc_indirect_lw = .FALSE.
 
! Indirect radiative effect of sea-salt
LOGICAL :: l_use_seasalt_indirect = .FALSE.

! Direct radiative effect of sea-salt
LOGICAL :: l_use_seasalt_direct = .FALSE.
      
LOGICAL :: l_use_soot_direct = .FALSE.   ! direct radiative effects of soot
LOGICAL :: l_use_soot_indirect = .FALSE. ! indirect effects of soot
      
! Use biomass aerosol for direct/indirect effect in radiation.
LOGICAL :: l_use_bmass_direct = .FALSE.
LOGICAL :: l_use_bmass_indirect = .FALSE.
      
! Use fossil-fuel organic carbon aerosol for direct/indirect
! effect in radiation
LOGICAL :: l_use_ocff_direct = .FALSE.
LOGICAL :: l_use_ocff_indirect = .FALSE.
      
! Use ammonium nitrate aerosol for direct/indirect effect in radiation
LOGICAL :: l_use_nitrate_direct = .FALSE.
LOGICAL :: l_use_nitrate_indirect = .FALSE.
      
! Use aerosol climatologies in radiation instead of prognostic variables
! Set on a species by species basis
LOGICAL :: l_use_arclbiom = .FALSE. ! biomass burning aerosol
LOGICAL :: l_use_arclblck = .FALSE. ! black carbon
LOGICAL :: l_use_arclsslt = .FALSE. ! sea salt
LOGICAL :: l_use_arclsulp = .FALSE. ! sulpahtes
LOGICAL :: l_use_arcldust = .FALSE. ! mineral dust
LOGICAL :: l_use_arclocff = .FALSE. ! organic carbon (fossil fuel)
LOGICAL :: l_use_arcldlta = .FALSE. ! delta aerosol

! Use droplet number from n_drop_pot array:
LOGICAL :: l_use_ndrop = .FALSE.

LOGICAL :: l_use_aod = .FALSE.  ! Aerosol optical depth diagnostic was requested
      
LOGICAL :: L_rad_deg = .FALSE.  ! controls the use of spatial degradation 
!                                 of radiation calc.

! Logicals for different radiation packages
LOGICAL :: l_forcing     = .FALSE. ! Calculate radiative forcings 
LOGICAL :: l_radiance    = .FALSE. ! Calculate radiances          
LOGICAL :: l_timestep    = .FALSE. ! Use new timestepping scheme 
LOGICAL :: l_rad_perturb = .FALSE. ! Use the perturbation version of 
                                   ! the radiative time-stepping
                                   ! Needed in glue_rad.
! Control integer for different radiation packages used by check_run_radiation()
! routine to set l_forcing, l_radiance, l_timestep and l_rad_perturb
INTEGER :: i_rad_extra_call = IMDI 

! Changes to open sea albedo for HadGEM1
LOGICAL :: l_use_spec_sea  = .FALSE.     ! Spectr. dep. sea albedos

! Use modulus of fluxes to remove negative effective extinctions
LOGICAL :: l_mod_k_flux = .FALSE.




! Use a solar zenith angle correction based on optical depth
LOGICAL :: l_rad_szacor = .FALSE.
!                       Needed in glue_rad-rad_ctl3c.
!                       Switch for the solar zenith angle correction to surface
!                       fluxes using the change in optical depth.      

! Scale the condensed water content to simulate
! inhomogeneous clouds
LOGICAL :: l_inhom_cloud = .FALSE.

! Orography correction to SW radiation
LOGICAL :: l_use_orog_corr  = .FALSE.  !  Find gradients from mean orog
LOGICAL :: l_use_grad_corr  = .FALSE.  !  Use ancillary X & Y gradients
! Correction for skyview factor in LW and direct SW
LOGICAL :: l_use_skyview    = .FALSE.
LOGICAL :: l_orog_unfilt    = .FALSE.  !  Use unfiltered ancillary orog 
! Control integer used by check_radiaiton to set the orography correction
! and skyview logicals
INTEGER :: i_rad_topography = IMDI

! ----------------------------

INTEGER :: h_swbands    = imdi   ! Number of shortwave radiation bands
INTEGER :: h_lwbands    = imdi   ! Number of longwave radiation bands
INTEGER :: a_sw_radstep = imdi   ! Number of advection steps per SW step
INTEGER :: a_lw_radstep = imdi   ! Number of advection steps per LW step


! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.
      
INTEGER :: a_sw_radstep_diag = imdi
! Number of advection steps per 'fast' SW step (3C)
INTEGER :: a_lw_radstep_diag = imdi
! Number of advection steps per 'fast' LW step (3C)
INTEGER :: a_sw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)
INTEGER :: a_lw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)

INTEGER :: i_ozone_int = imdi ! Option for interpolation of ozone


! The following three switches REMOVED from run_radiation NL (ROSE project)
! They are set in readlsta.F90, dependent on settings of cusack_aero and 
! cusack_aero_hgt, which have been added to the NL 

! True if climatological aerosol is included.
LOGICAL :: L_climat_aerosol = .FALSE.

! True to use real boundary layer heights to specify the boundary
! layer aerosol.
LOGICAL :: L_clim_aero_hgt = .FALSE.

! Flag to use HadGEM1 setting for climatological aerosols
LOGICAL :: L_HadGEM1_Clim_Aero = .FALSE.

! These two switches ADDED to NL (ROSE project)
INTEGER :: cusack_aero     = IMDI
INTEGER :: cusack_aero_hgt = IMDI

LOGICAL :: lrad_ccrad = .FALSE.
!             Allows access to ccrad code and the logicals
!             lrad_ovrlap and lrad_ccw_scav

LOGICAL :: lrad_diag_mode = .FALSE.
!             Needed in ni_rad_ctl.
!             Switch to set radiation calculations to "diagnostic mode"
!             i.e. the radiation code is isolated from the model
!             evolution (zero energy input) but diagnostics will still
!             be available.

LOGICAL :: l_sw_radiate  ! Activate SW radiation this timestep
LOGICAL :: l_lw_radiate  ! Activate LW radiation this timestep
LOGICAL :: l_sw_radiate_diag
! Activate fast SW radiation this timestep (3C)
LOGICAL :: l_lw_radiate_diag
! Activate fast LW radiation this timestep (3C)
LOGICAL :: l_sw_radiate_prog
! Activate slow SW radiation this timestep (3C)
LOGICAL :: l_lw_radiate_prog
! Activate slow LW radiation this timestep (3C)

! Convert zonal mean ozone to field
LOGICAL :: lexpand_ozone

! convert zonal mean tpps ozone to field
LOGICAL :: lexpand_tpps_ozone

! Tropopause-based Ozone Scheme
LOGICAL :: l_use_tpps_ozone = .FALSE.  !  Use TPPS ozone scheme

!-----------------------------------------------------------
! run_radiation namelists
! ----------------------------------------------------------

LOGICAL :: l_sec_var  = .FALSE.     ! true if using time varying astronomy
LOGICAL :: l_EqT      = .FALSE.     ! True if including the equation of time

LOGICAL :: l_rad_ovrlap   = .FALSE.     
!                           Requires l_ccrad=.TRUE.
!                           Allows Convective and LS Cloud to overlap 
!                           for radiative impacts.
!                           (Experimental, defaulted to FALSE,
!                           requires a hand-edit to change)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY).
!             Current convective/large-scale cloud fractions in the
!             radiation scheme are mutally exclusive. This assumes CCA
!             and the large-scale to overlap with CCA taking dominance.
!             I.E. Large-scale cloud fraction must exceed the convective
!             cloud fraction before having any presence.

LOGICAL :: l_rad_ccw_scav = .FALSE.
!                           Requires l_ccrad=.TRUE. .AND. l_rad_ovrlap=.TRUE.
!                           Allows Convective Cloud Water (CCW) to
!                           compensate for LS Cloud water in overlapping
!                           LS/CCA fractions.
!                           (Experimental, defaulted to FALSE, requires a
!                           hand-edit to change)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY)
!             Allowing the CCA to negate large-scale fractions of lower
!             values means that the large-scale cloud water in the
!             overlapping fraction is lost. This switch will scavenge
!             the large-scale cloud water from the overlapping fraction
!             and combine it with the convective cloud water to
!             conpensate.

LOGICAL :: l_rad_use_clim_volc=.FALSE. 
!                           If .TRUE. use climatological volcanic 
!                           eruption code in climatological aerosol
!                           code      


LOGICAL :: l_rad_snow_emis = .FALSE.
!          Switch to adjust the emissivity in radiation for snow cover.
!          This should eventually be moved to the surface scheme, but
!          JULES cannot currently cope with distinct emissivities
!          for snow-covered surfaces, so the switch currently acts
!          only in radiation and logically belongs here for the present.
LOGICAL :: l_t_land_nosnow = .FALSE.
!          Switch for emissivity of snow used in averaging the
!          surface temperature. Setting this switch to .TRUE. is 
!          deprecated and it is included only for historical reasons.
LOGICAL :: l_quad_t_coast = .FALSE.
!          Switch for quadratic averaging of the surface temperature at
!          coastal points. .FALSE. is deprecated.
LOGICAL :: l_t_rad_solid = .FALSE.
!          Switch to use common soid temperature at coastal points with
!          sea-ice. .TRUE. is deprecated.


! ------------------------------------------

! number of components of clouds
INTEGER,PARAMETER:: npd_cloud_component=4

! number of permitted types of clouds.
INTEGER,PARAMETER:: npd_cloud_type=4

! number of permitted representations of clouds.
INTEGER,PARAMETER:: npd_cloud_representation=4

! number of overlap coefficients for clouds
INTEGER,PARAMETER:: npd_overlap_coeff=18

! number of coefficients for two-stream sources
INTEGER,PARAMETER:: npd_source_coeff=2

! number of regions in a layer
INTEGER,PARAMETER:: npd_region=3


INTEGER :: a_sw_segments = imdi  ! No of batches used in shortwave code
INTEGER :: a_sw_seg_size = -99   ! Size of sw batches, -ve disabled 
INTEGER :: a_lw_segments = imdi  ! No of batches used in longwave code
INTEGER :: a_lw_seg_size = -99   ! Size of lw batches, -ve disabled 

INTEGER :: aero_bl_levels = imdi 
!                          Common number of layers taken to be 
!                          occupied by the boundary-layer
!                          aerosol if the boundary layer
!                          depth is not used to determine the 
!                          number separately at each grid-point
!                          In previous versions of the code,
!                          this was taken to be BL_LEVELS 

INTEGER :: clim_rad_volc_eruption_year = imdi  ! Climatological volcano 
!                                                eruption year 
      
INTEGER :: clim_rad_volc_eruption_month = imdi    ! Climatological volcano 
!                                                eruption month 
      
INTEGER :: rad_mcica_sampling = imdi   ! Version of McICA used (was 1)
!             Needed in open_cloud_gen. Selects the version of McICA
!             used to sample the generated cloud:
!                               0 = full sampling
!                               1 = single sampling
!                               2 = optimal sampling

! --------------------------------
                        
REAL    :: rad_mcica_sigma = rmdi  ! Normalised cloud condensate standard
!                                    deviation for the cloud generator.
!                                    Needed in open_cloud_gen.

INTEGER :: i_cloud_representation = imdi
INTEGER :: i_inhom = imdi
INTEGER :: i_overlap = imdi
INTEGER :: i_fsd = imdi
INTEGER :: i_cloud_representation_2 = imdi
INTEGER :: i_inhom_2 = imdi
INTEGER :: i_overlap_2 = imdi
INTEGER :: i_fsd_2 = imdi

REAL    :: co2_mmr = rmdi          ! CO2 concentration (if constant)

REAL    :: sc = rmdi               ! Solar constant (W/m2)

! Scaling factors to simulate inhomogeneous cloud.
REAL    :: inhom_cloud_sw(npd_cloud_component) = rmdi
REAL    :: inhom_cloud_lw(npd_cloud_component) = rmdi


! Decorrelation pressure scale for large scale cloud
REAL    :: dp_corr_strat = rmdi
      
! Decorrelation pressure scale for convective cloud
REAL    :: dp_corr_conv  = rmdi

      
REAL    :: clim_rad_volc_eruption_weight = rmdi
! Eruption weighting factor for idealised volcanic aerosol. 
! 1.0 is an average 20th century tropical explosive eruption.

REAL    :: aeroscl_csk_clim(5) = (/ rmdi, rmdi, rmdi, rmdi, rmdi /) 
! Scalings for aerosols in Cusack's climatology

! Number of radiation prognostic/diagnostic timesteps per day
INTEGER :: i_sw_radstep_perday_prog = IMDI
INTEGER :: i_lw_radstep_perday_prog = IMDI
INTEGER :: i_sw_radstep_perday_diag = IMDI
INTEGER :: i_lw_radstep_perday_diag = IMDI

! Ozone tracer as input to radiation scheme      
LOGICAL :: l_use_cariolle   = .FALSE.
LOGICAL :: l_use_ozoneinrad = .FALSE.

NAMELIST/RUN_Radiation/ &
       cusack_aero, cusack_aero_hgt, aeroscl_csk_clim,                 &
       a_sw_segments,a_sw_seg_size,a_lw_segments,a_lw_seg_size,co2_mmr,&
       sc, l_sec_var,l_EqT,inhom_cloud_sw,inhom_cloud_lw,dp_corr_strat,&
       rad_mcica_sampling, rad_mcica_sigma,                            &
       dp_corr_conv, aero_bl_levels,                                   &
       l_rad_use_clim_volc, clim_rad_volc_eruption_year,               &
       clim_rad_volc_eruption_month, clim_rad_volc_eruption_weight,    &
       alpham, ssalpham, alphac, ssalphac, alphab, dtice, ssdtice,     &
       dt_bare,dalb_bare_wet,pen_rad_frac,sw_beta,                     &
       n2ommr, ch4mmr, c11mmr, c12mmr,                                 &
       o2mmr, c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr,      &
       is_ncol, i_cloud_representation, i_cloud_representation_2,      &
       i_inhom, i_inhom_2, i_overlap, i_overlap_2, i_fsd, i_fsd_2,     &
       l_rad_snow_emis, l_t_land_nosnow, l_quad_t_coast, l_t_rad_solid,&
       h_swbands, h_lwbands, i_rad_extra_call, i_rad_topography,       &
       l_radiation, l_mod_k_flux, l_rad_deg, l_use_spec_sea,           &
       l_rad_szacor, i_sw_radstep_perday_prog,                         &
       i_lw_radstep_perday_prog,  i_sw_radstep_perday_diag,            &
       i_lw_radstep_perday_diag, l_use_cariolle, l_use_ozoneinrad,     &
       i_ozone_int

! Variables used to calculate radstep
INTEGER, PARAMETER   :: secs_per_day = 86400
INTEGER              :: secs_per_timestep 
INTEGER              :: timesteps_per_day

CONTAINS 

SUBROUTINE check_run_radiation()

! Description:
!   Subroutine to apply logic controls and set control variables based on the 
!   options selected in the run_radiation namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod,  ONLY: ereport
USE nlstgen_mod,  ONLY: secs_per_periodim, steps_per_periodim
USE submodel_mod, ONLY: a_im
IMPLICIT NONE

INTEGER                       :: icode         ! used for ereport
CHARACTER (LEN=80)            :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'check_run_bl'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_RADIATION',zhook_in,zhook_handle)

! Set logicals for different radiation packages based on a control 
! integer choice in namelist
IF (i_rad_extra_call == 0) THEN
  l_forcing     = .FALSE.
  l_timestep    = .FALSE.
  l_rad_perturb = .FALSE.
  l_radiance    = .FALSE.
ELSE IF (i_rad_extra_call == 1) THEN
  l_forcing     = .TRUE.
  l_timestep    = .FALSE.
  l_rad_perturb = .FALSE.
  l_radiance    = .FALSE.
ELSE IF (i_rad_extra_call == 2) THEN
  l_forcing     = .FALSE.
  l_timestep    = .TRUE.
  l_rad_perturb = .TRUE.
  l_radiance    = .FALSE.
ELSE IF (i_rad_extra_call == 3) THEN
  l_forcing     = .FALSE.
  l_timestep    = .FALSE.
  l_rad_perturb = .FALSE.
  l_radiance    = .TRUE.
ELSE
  WRITE (cmessage,'(A70)') 'i_rad_extra_call value invalid, default to ' &
                            // '0: single call to radiation'
  icode = -100
  CALL ereport(RoutineName, icode, cmessage)
END IF

! Set logicals for different radiation topography based on a control
! integer choice in namelist
IF (i_rad_topography == 0) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 1) THEN
  l_use_orog_corr   = .TRUE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 2) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .TRUE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 3) THEN
  l_use_orog_corr   = .TRUE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .TRUE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 4) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .TRUE.
  l_use_skyview     = .TRUE.
  l_orog_unfilt     = .TRUE.
ELSE
  WRITE (cmessage,'(A58)') 'i_rad_topography value invalid, default to ' &
                            // '0: flat surface'
  icode = -100
  CALL ereport(RoutineName, icode, cmessage)
END IF


! Determine number of advection timesteps between the radiation 
! diagnostics/prognostics steps
! nlstcgen not read in by reconfiguration or SCM so calculation of radstep 
! variables not possible.  These items are not needed by the reconfiguration 
! and for the SCM they are calculated in scm_main.

secs_per_timestep = secs_per_periodim(a_im) / steps_per_periodim(a_im)
timesteps_per_day = secs_per_day / secs_per_timestep

a_sw_radstep_prog = timesteps_per_day / i_sw_radstep_perday_prog
a_sw_radstep_diag = timesteps_per_day / i_sw_radstep_perday_diag
a_lw_radstep_prog = timesteps_per_day / i_lw_radstep_perday_prog
a_lw_radstep_diag = timesteps_per_day / i_lw_radstep_perday_diag
IF (i_rad_extra_call == 0) THEN 
  a_sw_radstep_diag = a_sw_radstep_prog
  a_lw_radstep_diag = a_lw_radstep_prog
END IF


IF (lhook) CALL dr_hook('CHECK_RUN_RADIATION',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_radiation

END MODULE rad_input_mod



