! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE stochastic_physics_run_mod

!  Global data module for switches/options concerned with stochastic physics
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Stochastic Physics

  ! Description:
  !   Module containing runtime logicals/options used by the SKEB2 
  !   and PR2 schemes.
  !   and the stochastic_physics_run_setup subroutine to control logic of 
  !   selected options.

  ! Method:
  !   All switches/options which are contained in the &RUN_Stochastic
  !   sub-namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate, but mostly
  !   rmdi values are set to force users to specify values through Rose.

  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.

  !   Note, flags for other stochastic physics routines, which were
  !   embedded within other namelists, have now been moved here.
  
  USE missing_data_mod, ONLY: rmdi, imdi

  IMPLICIT NONE
  !======================================================================
  ! Logical switches set from RUN_Stochastic namelist via UI
  !======================================================================

  LOGICAL :: l_skeb2  = .FALSE.                                         &
             ! switch controls use of SKEB2 scheme
 ,           l_skeb2_psicdisp = .FALSE.                                 &
             ! TRUE = incl streamfunction modulation by convection
 ,           l_skeb2_psisdisp = .FALSE.                                 &
             ! TRUE = incl streamfunction modulation by smagorinksy
 ,           l_skeb2_skeb1disp = .FALSE.                                &
             ! TRUE = incl streamfunction modulation by SKEB1-type
!!! Should not need both of these and only if l_skeb2_psicdisp = TRUE
!!! ,           l_skeb2_cdisp_cape = .FALSE.                               &
             ! TRUE = Calc Dissipation from convection using CAPE
!!! ,           l_skeb2_cdisp_mflx = .FALSE.                               &
             ! TRUE = Calc Dissipation from convection using Mass-Flux
 ,           l_skeb2_velpot = .FALSE.                                   &
             ! TRUE = Calc divergent-wind incr from VelPot forcing
 ,           l_rp2 = .FALSE.                                            &
             ! switch controls use of RP2 scheme
 ,           l_skebsmooth_adv = .FALSE.                                 &
             ! TRUE = Perform advanced smoothing of energy diss fields
 ,           l_skebprint = .FALSE.                                      
             ! TRUE = Print global KE backscattered at each timestep

  ! l_skeb2_cdisp_cape/mflx only valid if l_skeb2_psicdisp is true
  ! and only one or other should be used. Replace with integer type
    INTEGER :: type_cape = 4 ! Calc Conv Dissipation using CAPE
    INTEGER :: type_mflx = 5 ! Calc Conv Dissipation using Mass-Flux
    INTEGER :: skeb2_cdisp = imdi ! type of streamfn modulation by conv

  !======================================================================
  ! Integer options set from RUN_Stochastic namelist via UI
  !======================================================================
  !
  !   N1 and N2 define the range of spherical harmonic orders over
  !   which backscatter is applied (remember: N(N+1)/R*R is the
  !   effective wavenumber squared). 40<N<80 corresponds to horizontal
  !   wavelengths in the range 500 km to 1000 km.
  !
  !   CAREFUL: This value has to be changed to the maximum wavenumber
  !        (the N number) solved by the model version
  !
  !   skeb2_toplev defines the top level of the streamfunction modulation
  !   field. this will change with vertical resolution. At 38 levels this
  !   is set=33 (~ 21.8km) to avoid making changes above the stratosphere
  !   skeb2_botlev defines the bottom level of the streamfunction
  !   modulating field (we may want to avoid changes in the boundary
  !   layer, certainly in level one)
  !

 INTEGER :: stphseed = imdi       ! Control variable for random seed options:
                                  ! 0 => No random seed file
                                  ! 1 => Read random seed from file
                                  ! 2 => Write random seed to file

!!! INTEGER :: n1=5                  ! minimum wavenumber for backscatter
!!! INTEGER :: n2=60                 ! maximum wavenumber for backscatter
!!! INTEGER :: skeb2_toplev=33       ! Top level of SKEB2 calculations
!!! INTEGER :: skeb2_botlev=2        ! Bottom level of SKEB2 calculations
!!! INTEGER :: nsmooth = 5           ! Iteration count for spatial smoothing
!!! INTEGER :: rhcrit_ref_level = 3  ! RHCrit reference level for RP2
!!! INTEGER :: ran_max = 1           ! number of independent RP2 variations
!!! INTEGER :: ran_count = 1         ! counter for random number array in RP
 INTEGER :: n1=imdi               ! minimum wavenumber for backscatter
 INTEGER :: n2=imdi               ! maximum wavenumber for backscatter
 INTEGER :: skeb2_toplev=imdi     ! Top level of SKEB2 calculations
 INTEGER :: skeb2_botlev=imdi     ! Bottom level of SKEB2 calculations
 INTEGER :: nsmooth = imdi        ! Iteration count for spatial smoothing
 INTEGER :: rhcrit_ref_level=imdi ! RHCrit reference level for RP2
 INTEGER :: offx_stph = imdi      ! Halo size used in spatial smoothing
 INTEGER :: offy_stph = imdi      ! Halo size used in spatial smoothing
                                  ! Both set to nsmooth in sthp_setup
 INTEGER :: ran_max = imdi        ! number of independent RP2 variations
 INTEGER :: ran_count = imdi      ! counter for random number array in RP

  !======================================================================
  ! Real values set from  RUN_Stochastic namelist via UMUI
  !======================================================================
 REAL :: tau = rmdi ! decorrelation time (~5.5 hrs is typical)
 REAL :: tot_backscat = rmdi ! global-mean rate of energy backscatter 
                             ! in m**2 s**(-3)
 REAL :: br = rmdi ! backscatter ratio (as fraction of diss. energy)
 REAL :: sdispfac = rmdi ! Multiplication factor for numerical dissipation 
                         ! field (this has been determined empirically)
 REAL :: cdispfac = rmdi ! Multiplication factor for convection dissipation 
                         ! field (this has been determined empirically)
 REAL :: kdispfac = rmdi ! Multiplication factor for SKEB1 (KE) dissipation 
                         ! field (this has been determined empirically)
 REAL :: alphac = rmdi ! Updraught proportion of gridbox (0.2% typical)

!!! REAL :: rv=5 ! rv: Vertical correlation range (in num of levels)
              ! will this change with L70 changes (not in UMUI - unused!)
!!! ,       tau = 2.e04                                                    &
            ! tau is the decorrelation time (~5.5 hrs in this case)
!!! ,       tot_backscat = 1.e-4                                           &
            ! global-mean rate of energy backscatter in m**2 s**(-3)
!!! ,       br = 0.2                                                       &
            ! backscatter ratio (fraction of diss. energy backscattered)
!!! ,       sdispfac = 2.0                                                 &
             ! Multiplication factor for numerical dissipation field
             ! (this has been determined empirically)
!!! ,       cdispfac = 1.0                                                 &
             ! Multiplication factor for convection dissipation field
             ! (this has been determined empirically)
!!! ,       kdispfac = 0.5 
             ! Multiplication factor for SKEB1 (KE) dissipation field
             ! (this has been determined empirically)
!!! ,       alphac = 2.0e-3
            ! Updraught proportion of gridbox (0.2%)

  REAL :: entcoef_min = rmdi        ! Minimum entrainment rate coefficient
  REAL :: entcoef_max = rmdi        ! Maximum entrainment rate coefficient
  REAL :: cape_timescale_min = rmdi ! Minimum CAPE closure timescale
  REAL :: cape_timescale_max = rmdi ! Maximum CAPE closure timescale
!!!  REAL :: entcoef_min = 2.75        ! Minimum entrainment rate coefficient
!!!  REAL :: entcoef_max = 4.0         ! Maximum entrainment rate coefficient
!!!  REAL :: cape_timescale_min = 1800.0 ! Minimum CAPE closure timescale
!!!  REAL :: cape_timescale_max = 3600.0 ! Maximum CAPE closure timescale

  REAL :: rhcrit_max          = rmdi
  REAL :: rhcrit_min          = rmdi
! RHCrit minimum and maximum default values
!!!  REAL :: rhcrit_max          = 0.906
!!!  REAL :: rhcrit_min          = 0.874

!------------------------------------

  REAL :: m_ci                = rmdi
  REAL :: m_ci_max            = rmdi
  REAL :: m_ci_min            = rmdi
! Stochastic physics ice fallspeed multiplier m_ci, 
! with minimum and maximum default values:
!!!  REAL :: m_ci                = 1.000
!!!  REAL :: m_ci_max            = 1.400
!!!  REAL :: m_ci_min            = 0.600

!------------------------------------

! ROSE - BOUNDARY LAYER stochastic items moved here from bl_option_mod 
! Maximum and minimum values for the STPH_RP scheme
! Boundary Layer
REAL :: par_mezcla_max= rmdi ! was 0.5
REAL :: par_mezcla    = rmdi ! was 0.15
REAL :: par_mezcla_min=rmdi ! was 0.05
! Max, mean and min value for the neutral mixing length
REAL :: g0_rp_max = rmdi ! was 20.0
REAL :: g0_rp  = rmdi ! was 10.0
REAL :: g0_rp_min = rmdi ! was 5.0
! Max,mean and min values for the flux profile parameter
REAL :: charnock_max= rmdi ! was 0.026
REAL :: charnock_min= rmdi ! was 0.01
! Max and min values for the charnock parameter
REAL :: lambda_min_rp_max= rmdi ! was 100.0
REAL :: lambda_min_rp = rmdi ! was 40.0
REAL :: lambda_min_rp_min= rmdi ! was 20.0
! Max, mean and min values for the minimum mixing length
REAL :: ricrit_rp_max= rmdi ! was 1.0
REAL :: ricrit_rp = rmdi ! was 1.0
REAL :: ricrit_rp_min= rmdi ! was 0.25
! Max, mean and min values for the critical Ri
REAL :: a_ent_1_rp_max = rmdi ! was 0.4
REAL :: a_ent_1_rp  = rmdi ! was 0.23
REAL :: a_ent_1_rp_min = rmdi ! was 0.1
! Max, mean and min values for the entrainment parameter A1
REAL :: g1_rp_max= rmdi ! was 1.5
REAL :: g1_rp = rmdi ! was 0.85
REAL :: g1_rp_min= rmdi ! was 0.5
! Max, mean and min values for the velocity scale parameter

! Gravity Wave Drag physics parameters 
! maximum and minimum default values
REAL :: kay_gwave_max = rmdi
REAL :: kay_gwave_min = rmdi
REAL :: gwd_frc_max   = rmdi
REAL :: gwd_frc_min   = rmdi
! Values kept from gravity_wave_drag/g_wave_input_mod.F90 for reference
!!!REAL :: kay_gwave_max = 4400.000
!!!REAL :: kay_gwave_min = 2200.000
!!!REAL :: gwd_frc_max   = 6.000
!!!REAL :: gwd_frc_min   = 2.000

  !======================================================================
  ! Arrays from Convection required by SKEB2
  !======================================================================
 REAL, ALLOCATABLE ::                                                   &
       skeb2_up_flux(:, :, :)                                           &
            !  updraught mass flux
 ,     skeb2_dwn_flux(:, :, :)                                          &
            ! downdraught mass flux
 ,     skeb2_cape(:, :)
            ! CAPE

  !======================================================================
  ! Define Namelist RUN_Stochastic
  !======================================================================
 NAMELIST/run_stochastic/                                               &
       l_skeb2                                                          &
 ,     l_rp2                                                            &
 ,     n1                                                               &
 ,     n2                                                               &
 ,     br                                                               &
 ,     tot_backscat                                                     &
 ,     tau                                                              &
 ,     alphac                                                           &
 ,     l_skeb2_psicdisp, l_skeb2_psisdisp, l_skeb2_skeb1disp            &
!!! ,     l_skeb2_cdisp_cape, l_skeb2_cdisp_mflx                           &
 ,     sdispfac, cdispfac, kdispfac, skeb2_cdisp, nsmooth               &
 ,     skeb2_toplev, skeb2_botlev, l_skeb2_velpot, ran_max              &
 ,     rhcrit_ref_level, l_skebsmooth_adv, l_skebprint, stphseed        &
 ,     entcoef_min, entcoef_max, cape_timescale_min, cape_timescale_max &
 ,     m_ci, m_ci_max, m_ci_min, rhcrit_max, rhcrit_min                 &
 ,     kay_gwave_max, kay_gwave_min, gwd_frc_max, gwd_frc_min           &
 ,     par_mezcla_max, par_mezcla, par_mezcla_min                       &
 ,     g0_rp_max, g0_rp, g0_rp_min, charnock_max, charnock_min          &
 ,     lambda_min_rp_max, lambda_min_rp, lambda_min_rp_min              &
 ,     ricrit_rp_max, ricrit_rp, ricrit_rp_min                          &
 ,     a_ent_1_rp_max, a_ent_1_rp, a_ent_1_rp_min                       &
 ,     g1_rp_max, g1_rp, g1_rp_min

  !======================================================================
  ! Array for spatial 1-2-1 smoothing
  !======================================================================
 REAL, ALLOCATABLE, SAVE ::                                             &
       mask_pdamp(:, :)                                                 &
            !  Array of pattern damping coefficients for SKEB2
 ,     mask_smooth(:, :)
            !  Array of smoothing coefficients for SKEB2

!======================================================================
! LOGICAL switches not set in namelist
!======================================================================

LOGICAL  ::  l_stphseed_read = .FALSE.                                  
             ! TRUE = Read in previously specified seed
LOGICAL  ::  l_stphseed_write = .FALSE. 
             ! TRUE = WRITE out seed

CONTAINS 

SUBROUTINE check_run_stochastic()

! Description:
!   Subroutine to apply logic controls and set control variables based on the 
!   options selected in the run_stochastic namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_STOCHASTIC',zhook_in,zhook_handle)

! Set random seed logicals based on integer choice in namelist
IF (stphseed == 0) THEN
   l_stphseed_read  = .FALSE.
   l_stphseed_write = .FALSE.  
ELSE IF (stphseed == 1) THEN
   l_stphseed_read  = .TRUE.
   l_stphseed_write = .FALSE.       
ELSE IF (stphseed == 2) THEN
   l_stphseed_read  = .FALSE.
   l_stphseed_write = .TRUE.       
END IF

IF (lhook) CALL dr_hook('CHECK_RUN_STOCHASTIC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_stochastic

END MODULE stochastic_physics_run_mod
