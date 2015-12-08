! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with convection.

MODULE cv_run_mod

  ! Description:
  !   Module containing runtime logicals/options used by the convection code.
  !
  ! Method:
  !   All switches/options which are contained in the &Run_convection 
  !   sub-namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate.
  !   
  !   Any routine wishing to use these options may do so with the 'Use' 
  !   statement.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v7.4 programming standards.
  !
  ! Declarations:

  USE missing_data_mod, ONLY: RMDI, IMDI

  IMPLICIT NONE
  SAVE

  ! To satisfy the requirements of ROSE, the defaults for all these
  ! switches/variables have been set to FALSE/MDI. 
  ! The original defaults are noted in comments below. 
  ! These values are also indicated in the ROSE help text.

  INTEGER :: i_convection_vn = imdi
                        ! Switch to determine version of convection scheme
                        ! 0 => 0A scheme (i.e. diagnostics only)
                        ! 4 => 4A scheme
                        ! 5 => 5A scheme
                        ! 6 => 6A scheme
  INTEGER, PARAMETER :: i_convection_vn_0a = 0
  INTEGER, PARAMETER :: i_convection_vn_4a = 4
  INTEGER, PARAMETER :: i_convection_vn_5a = 5
  INTEGER, PARAMETER :: i_convection_vn_6a = 6

  !===========================================================================
  ! Logical switches set from UI
  !===========================================================================

  ! Comments for TRUE status

  Logical :: l_fix_udfactor = .FALSE. ! Fix application of UD_FACTOR (lock)

  Logical :: l_cloud_deep   = .FALSE. ! Use depth criterion for convective anvils
                                      ! Original default TRUE
  Logical :: l_mom          = .FALSE. ! Use Convective Momentum Transport

  Logical :: l4a_kterm      = .FALSE. ! Use kterm instead for Deep CMT (4A)
                                      ! Original default TRUE
  Logical :: l_eman_dd      = .FALSE. ! Use Emanuel downdraught scheme (4A)
                          
  Logical :: l_rediagnosis  = .FALSE. ! If more than one convection step per model
                                      ! physics timestep then rediagnose cumulus
                                      ! points by recalling conv_diag before each
                                      ! convective sub-step. 
                                      ! Original default TRUE

  Logical :: l_anvil        = .FALSE. ! Apply anvil scheme to 3D convective cloud.
                                      ! No effect unless l_3d_cca = .TRUE.

  Logical :: l_dcpl_cld4pc2 = .FALSE. ! Decouples section 0/5 cloud properties from
                                      ! each other, use with PC2. May change results
                                      ! if used without PC2.

  LOGICAL :: l_murk_conv    = .FALSE. ! Enables convective mixing of MURK as 
                                      ! a tracer

  LOGICAL :: l_safe_conv    = .FALSE. ! Safer convection, remove negative q 
                                      ! before attempting convection, don't 
                                      ! add increments for ascents with
                                      ! negative CAPE. (5A/6A only)

  LOGICAL :: l_ccrad        = .FALSE. ! Main Logical, will include code
                                      ! connected with CCRad.
                                      ! (including bugfixes)

  LOGICAL :: l_3d_cca       = .FALSE. ! Use 3D convective cloud amount

  LOGICAL :: l_conv_hist    = .FALSE. ! True if 3 extra prognostics holding 
                                      ! convective history information.

  LOGICAL :: l_param_conv   = .FALSE. ! Run time switch for convection scheme
  
  !===========================================================================
  ! Logical switches not set from UI
  !===========================================================================
  
  LOGICAL :: l_phase_lim    = .TRUE.  ! Limits phase change of precipitation
                                      ! so Latent Heat does not cause 
                                      ! Temperature to cross Melting
                                      ! Temperature as in Hadam3 physics.

  LOGICAL :: l_pc2_diag_sh  = .FALSE. ! If true uses  diagnostic convective 
                                      ! shallow cloud in PC2 replacing 
                                      ! prognotsic qcl etc

  LOGICAL :: l_cv_conserve_check = .FALSE. 
                                      ! Diagnostic conservation checks. Goes
                                      ! through energy correction code 
                                      ! (5A/6A only.

  !===========================================================================
  ! Integer options set from UI
  !===========================================================================

  ! Convection integer options set from UMUI Convection Scheme

  Integer :: n_conv_calls       = imdi ! Number of calls to convection
                                       ! per physics timestep
                                       ! Original default 1

  Integer :: a_convect_segments = imdi ! No of batches used in convection
                                       ! Original default 1

  Integer :: a_convect_seg_size = imdi ! Size of convection segments. Can be 
                                       ! specified as an alternative to no. of
                                       ! segments. Original default -99

  Integer :: cld_life_opt       = imdi ! Convective cloud decay time
                                       ! Original default cld_life_constant (0)   
  Integer :: rad_cloud_decay_opt= imdi ! Convective cloud decay
                                       ! Original default rad_decay_off (0)
  Integer :: anv_opt            = imdi ! Anvil cloud basis
                                       ! Original default anv_model_levels (2)

  ! Convection Scheme Options (5A)
  ! NOTE: These options were valid at the time of writing. They are used in
  !       development code(5A) and so very likely to change.
  !       Users should consult the Convection Group for current available
  !       options.

  Integer :: iconv_shallow  = imdi  ! Shallow (Original default 0)
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Turbulence scheme (non-precipitating)
                             !   3: Turbulence scheme (precipitating)

  Integer :: iconv_congestus = imdi ! Congestus (Original default 0)
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use
   
  Integer :: iconv_mid       = imdi ! Mid-level (Original default 0)
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use

  Integer :: iconv_deep      = imdi ! Deep (Original default 0)
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use

  Integer :: deep_cmt_opt    = imdi ! Deep CMT tunings (Original default 0)
                             !   0: 4A Scheme (turbulence based) 
                             !   1: Operational 70 level (modified 4A scheme)
                             !   2: Gregory-Kershaw scheme
                             !   3: New turbulence based scheme.
                             !   4: Future use

  Integer :: mid_cmt_opt     = imdi ! Mid CMT scheme to be used (Orig default 0)
                             !   0: Gregory-Kershaw scheme 
                             !   1: Diffusive scheme
                             
  Integer :: icvdiag    = imdi ! Diagnosis calculation options (Orig default 1)
                             !   0: 4A Scheme (Default)
                             !   1: Improved 4A Scheme
                             !   2: Future/Experimental (dilute parcel)

  INTEGER :: cvdiag_inv = imdi ! Inversion test in convective diagnosis? 5A & 6A
                             ! Original default 1
                             ! When doing an undilute parcel ascent
                             ! Does not apply to 4A code.
                             !   0: No inversion test
                             !   1: Original inversion test (Default 5A and 6A)
                             !   2: Future alternative inversion tests
  
  Integer :: tv1_sd_opt = imdi ! Standard dev of level virtual temperature options
                             ! Original default 0
                             !   0: Assume BL depth of 300m (Default)
                             !   1: Use calculated BL depth
                             !   2: As (2) plus stability dependence and coastal mean


  Integer :: adapt    = imdi ! Adaptive detrainment/entrainment options
                             ! Original default 0
                             !   0: 4A Scheme (Default)
                             !   1: Adaptive detrainment: mid and deep convection
                             !   2: Future/Experimental (En/Detrainment)
                             !   3: Adaptive detrainment: deep convection
                             !   4: Adaptive detrainment: shallow, mid and 
                             !      deep convection
                             !   5: Smoothed adaptive detrainment: 
                             !      mid and deep convection
                             !   6: Smoothed adaptive detrainment: 
                             !      shallow, mid and deep convection
                             !   7: Improved smoothed adaptive detrainment: 
                             !      mid and deep convection
                             !   8: Improved smoothed adaptive detrainment: 
                             !      shallow, mid and deep convection

  INTEGER :: ent_opt_dp = imdi  ! Deep entrainment option (Orig default 0)
                             !   0: original Ap/(p*)^2 
                             !   1: n/z dependence where n=ent_fac
                             !   2: As Ap/(p*)^2 but multiplied by extra factor
                             !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))

  INTEGER :: ent_opt_md = imdi  ! mid entrainment option (Orig default 0)
                             !   0: original Ap/(p*)^2 
                             !   1: n/z dependence where n=ent_fac
                             !   2: As Ap/(p*)^2 but multiplied by extra factor
                             !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))

  Integer :: cape_opt   = imdi  ! CAPE closure options (Orig default 0)
                             !   0: RH based timescale
                             !   1: RH based timescale (timestep limited)
                             !   2: Fixed timescale
                             !   3: RH based timescale (timestep limited)
                             !      plus w test
                             !   4: Area scaled
                             !   5: Experimental

  Integer :: cape_bottom = imdi    ! Start level for w_max in column
                                   ! Original default IMDI

  Integer :: cape_top    = imdi    ! End   level for w_max in column
                                   ! Original default IMDI

  Integer :: sh_pert_opt  = imdi   ! Initial perturbation method for 
                                   ! shallow cumulus (Orig default 0)
                                   ! 0 = Original code
                                   ! 1 = Revised  code

  Integer :: limit_pert_opt = imdi ! Limits convective parcel perturbation
                                   ! to physically sensible values.
                                   ! Orig default 0
                                   ! 0 = original code - no limits
                                   ! 1 = apply limits to main ascent only
                                   ! 2 = apply limits to main ascent and in 
                                   !     the convection diagnosis


  Integer :: dd_opt         = imdi ! Downdraught scheme options 
                                   ! Orig default 0  
                                   ! 0 = Original code
                                   ! 1 = Revised  code

  Integer :: termconv       = imdi ! Original default 0
                                   ! 0 for default
                                   ! 1 for modified termination condition

  Integer :: bl_cnv_mix     = imdi ! Options for mixing convective increments in the BL
                                   ! Original default 0
                                   ! 0: original code 
                                   ! 1: only mix the increments from the initial 
                                   !    parcel perturbation

  Integer :: cnv_wat_load_opt=imdi ! Options for including liquid and frozen water loading
                                   ! in the convective updraught buoyancy calculation
                                   ! Original default 0
                                   ! 0: Do not include water loading (default)
                                   ! 1: Include water loading

  Integer :: cca2d_sh_opt = imdi   ! Method to evaluate cca2d (Shallow)
                                   ! Original default 0
  Integer :: cca2d_md_opt = imdi   ! Method to evaluate cca2d (Mid-level)
                                   ! Original default 0
  Integer :: cca2d_dp_opt = imdi   ! Method to evaluate cca2d (Deep)
                                   ! Original default 0

  Integer :: ccw_for_precip_opt   ! Option controlling critical cloud water 
                                  ! for the formation of precipitation
                                  ! Original default 0
                                  ! 0 - original code
                                  ! 1 - no Dcrit option
                                  ! 2 - Manoj's first function 
                                  ! 3 - Manoj's congestus function

  INTEGER :: plume_water_load     ! Option for water loading in undilute parcel
                                  ! ascent
                                  ! Original default 0
                                  ! 0 - no water removal original undilute
                                  ! 1 - remove any water > 1g/kg
                                  ! 2 - remove any water > profile varying with 
                                  !     qsat 

  INTEGER :: dil_plume_water_load ! Option for water loading in dilute parcel
                                  ! ascent
                                  ! Original default 0
                                  ! 0 - no water removal 
                                  ! 1 - remove any water > 1g/kg 
                                  ! 2 - remove any water > profile varying with 
                                  !     qsat 

  !===========================================================================
  ! Real values set from UI
  !===========================================================================

  Real :: cca_sh_knob = rmdi     ! Scales Shallow cloud fraction (CCRad)
                                 ! Original default 1.0
  Real :: cca_md_knob = rmdi     ! Scales Mid     cloud fraction (CCRad)
                                 ! Original default 1.0
  Real :: cca_dp_knob = rmdi     ! Scales Deep    cloud fraction (CCRad)
                                 ! Original default 1.0
  Real :: ccw_sh_knob = rmdi     ! Scales Shallow cloud water (CCRad)
                                 ! Original default 1.0
  Real :: ccw_md_knob = rmdi     ! Scales Mid     cloud water (CCRad)
                                 ! Original default 1.0
  Real :: ccw_dp_knob = rmdi     ! Scales Deep    cloud water (CCRad)
                                 ! Original default 1.0

  Real :: fixed_cld_life = rmdi ! Fixed convective cloud lifetime decay value
                                ! (seconds) (Original default 7200.0)

  Real :: cca_min = rmdi  ! Threshold value of convective cloud fraction
                          ! below which cca has neglible radiative impact and
                          ! is reset to zero (Original default 0.02)

  Real :: r_det   = rmdi  ! Parameter controlling adaptive detrainment -
                          ! Orig default 0.75 (operational)
                          ! HadGEM1a recommended 0.5

  Real :: cape_min     = rmdi ! Scale dependent min cape
                              ! Original default RMDI
  Real :: w_cape_limit = rmdi ! Test w for scale dependent cape timescale
                              ! Original default RMDI
  Real :: mparwtr      = rmdi ! Maximum value of the function that is used to calculate 
                              ! the maximum convective cloud water/ice in a layer (kg/kg)
                              ! Original default 1.0e-3
  Real :: qlmin        = rmdi ! Minimum value of the function that is used to calculate 
                              ! the maximum convective cloud water/ice in a layer (kg/kg)
                              ! Original default 2.0e-4
  Real :: fac_qsat     = rmdi ! Factor used to scale environmental qsat to give the
                              ! the maximum convective cloud water/ice in a layer
                              ! Original default RMDI
  Real :: mid_cnv_pmin = rmdi ! The minimum pressure (max height) at which mid 
                              ! level convection is allowed to initiate (Pa)
                              ! Original default 0.0
  Real :: amdet_fac    = rmdi ! Factor multiplying (1-rh) for adaptive mixing 
                              ! detrainment rate.
                              ! Original default 1.0
  ! Scales with cca2d to determine convective cloud amount with anvil
  Real :: anvil_factor = rmdi ! x cca2d = max cca, capped at 1.0
                              ! Original default RMDI
  Real :: tower_factor = rmdi ! x cca2d = min cca
                              ! Original default RMDI
  Real :: ud_factor    = rmdi ! Updraught factor used in calculation of 
                              ! convective water path Original default RMDI

  Real :: tice         = rmdi ! Phase change temperature in plume
                              ! Original default 273.15
  Real :: qstice       = rmdi ! Qsat at phase change temperature
                              ! (freeze/melting temperature)
                              ! Original default 3.5E-3
! 4A code only
  Real :: ent_fac      = rmdi ! Factor multiplying entrainment rate - deep & mid
                              ! Original default 1.0
! 5A & 6A code only
  Real :: ent_fac_dp      = rmdi ! Factor multiplying entrainment rate - deep
                                 ! Original default 1.0
  Real :: ent_fac_md      = rmdi ! Factor multiplying entrainment rate - mid-level
                                 ! Original default 1.0
  Real :: ent_dp_power    = rmdi ! Power n for (p/p*)^n for entrainment option 3
                                 ! Original default 2.0
  Real :: ent_md_power    = rmdi ! Power n for (p/p*)^n for entrainment option 3    
                                 ! Original default 2.0
  Real :: cape_timescale  = rmdi ! Timescale for CAPE closure.
                                 ! Original default RMDI
  REAL :: cvdiag_sh_wtest = rmdi ! w for air above shallow convection must be 
                                 ! less than this. (Default value is 0.0)
                                 ! Original default 0.0
! 6A code only
  REAL :: eff_dcfl        = rmdi ! Factor the defines the efficiency by which
                                 ! detrained liquid condensate creates liquid 
                                 ! cloud fraction
  REAL :: eff_dcff        = rmdi ! Factor the defines the efficiency by which
                                 ! detrained frozen condensate creates frozen 
                                 ! cloud fraction

!------------------------------------------------------------------------------
! Switches added for Warm TCS scheme
!------------------------------------------------------------------------------
! The following have been removed from the RUN_convection name list 
!  - defaults are hard wired

  INTEGER :: iwtcs_diag1 = 0  ! Options for WTCS diagnosis
  INTEGER :: iwtcs_diag2 = 0  ! Options for WTCS diagnosis

!------------------------------------------------------------------------------  
! Further changes

  REAL :: qmin_conv = 1.0e-8  ! Minimum allow value of q after convection
                              ! Also used for negative q check.  
  

!------------------------------------------------------------------------------
! Define namelist &Run_Convection read in from CNTLATM control file.
! Changes made to this list will affect both the Full UM and the SCM
!------------------------------------------------------------------------------

        Namelist/Run_Convection/                                              &

        i_convection_vn,                                                      &

        ! Logical switches
        l_mom,    l_fix_udfactor, l4a_kterm,                                  &
        l_eman_dd,l_cloud_deep,                                               &
                  l_rediagnosis,  l_dcpl_cld4pc2,                             &
        l_anvil,  l_murk_conv, l_safe_conv, l_cv_conserve_check,              &
        l_ccrad,              l_3d_cca,               l_conv_hist,            &
        l_param_conv,                                                         &
        ! General scheme options/variables
        n_conv_calls,         a_convect_segments,     a_convect_seg_size,     &
        sh_pert_opt,                                                          &
        dd_opt,               deep_cmt_opt,           mid_cmt_opt,            &
        termconv,             adapt,                  r_det,                  &
        tice,                 qstice,                 ent_fac,                &
        ent_fac_dp,           ent_fac_md,                                     &
        ent_opt_dp,           ent_opt_md,                                     &
        ent_dp_power,         ent_md_power,                                   &
        bl_cnv_mix,                                                           &
        mid_cnv_pmin,         amdet_fac,              ccw_for_precip_opt,     &
        cnv_wat_load_opt,     tv1_sd_opt,             limit_pert_opt,         &

        ! Convective diagnosis options
        icvdiag,              plume_water_load,       dil_plume_water_load,   &
        cvdiag_inv,           cvdiag_sh_wtest,                                &

        ! Cape related options/variables
        cape_opt,             cape_bottom,            cape_top,               &
        cape_timescale,       w_cape_limit,           cape_min,               &

        ! Convective cloud options/variables
        cld_life_opt,         rad_cloud_decay_opt,    cca_min,                &
        fixed_cld_life,       ud_factor,              mparwtr,                &
        qlmin,                fac_qsat,               eff_dcfl,               &
        eff_dcff,                                                             &

        ! CCRad options options/variables
        cca2d_sh_opt,         cca_sh_knob,            ccw_sh_knob,            &
        cca2d_md_opt,         cca_md_knob,            ccw_md_knob,            &
        cca2d_dp_opt,         cca_dp_knob,            ccw_dp_knob,            &

        ! Anvil scheme options
        anvil_factor,         tower_factor,           anv_opt,                &
 
        ! Convection type options
        iconv_shallow,        iconv_mid,              iconv_deep,             &
        iconv_congestus    

!------------------------------------------------------------------------------

END MODULE cv_run_mod
