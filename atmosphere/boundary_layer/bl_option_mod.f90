! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!+ Data module for switches/options concerned with the BL scheme.
! Description:
!   Module containing runtime options/data used by the boundary
!   layer scheme and check_run_bl routine for logic checking the
!   selected options.

! Method:
!   Switches and associated data values used by the boundary layer
!   scheme are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'USE'
!   statement.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer

! Code Description:
!   Language: FORTRAN 90

MODULE bl_option_mod

USE  missing_data_mod, ONLY: rmdi, imdi
USE  control_max_sizes, ONLY: max_number_alpha_cds
USE  visbty_constants_mod, ONLY: calc_prob_of_vis
USE  mym_option_mod, ONLY: l_mono_adv_turb, l_conserv_adv_turb,         &
      high_order_scheme_adv_turb, l_my_initialize,                      &
      l_my_ini_zero, l_local_above_tkelvs, my_z_limit_elb, wb_ng_max,   &
      shcu_levels, my_ini_dbdz_min, monotone_scheme_adv_turb,           &
      l_print_max_tke, l_my_lowest_pd_surf_tqc, tke_levels,             &
      tke_cm_fa, my_lowest_pd_surf, tke_cm_mx, bdy_tke, tke_dlen,       &
      my_z_extra_fact, l_my_prod_adj, my_prod_adj_fact,                 &
      l_my_extra_level, l_my_condense, l_shcu_buoy, l_adv_turb_field

! Declarations:

IMPLICIT NONE

!======================================================================= 
! BL control 
!=======================================================================
! Switch for turning off boundary layer code
LOGICAL :: L_bl = .FALSE. !  F: Turns off boundary layer code  

!======================================================================= 
!   Permissible settings for BL options.
!=======================================================================
INTEGER, PARAMETER :: off = 0  ! Switch disabled
INTEGER, PARAMETER :: on  = 1  ! Switch enabled

!     Options for non-gradient stress following
INTEGER, PARAMETER :: BrownGrant97 = 1
INTEGER, PARAMETER :: BrownGrant97_limited = 2
!       Brown and Grant (1997), version 2 including a limit on its size

!     Options for flux gradient formulation
INTEGER, PARAMETER :: Locketal2000   = 0

!       Flux gradients as in Lock et al. (2000)
INTEGER, PARAMETER :: HoltBov1993 = 1

!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)

!     Options for entrainment enhancement in Sc over Cu
INTEGER, PARAMETER :: Buoyrev_feedback = 1

!     Options for form drag
INTEGER, PARAMETER :: No_drag         = 0
INTEGER, PARAMETER :: Effective_z0    = 1
INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory

INTEGER, PARAMETER :: DynDiag_ZL = 1
INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used as a dynamic criterion in the
!       diagnosis of BL types: version 2 includes changes to
!       cope with BL_LEVELS >> 3km
INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
!       As 2 but only applied to points diagnosed with Cumulus 
!       and strictly for sea points (fland<0.01, cf 0.5)
INTEGER, PARAMETER :: DynDiag_Ribased = 4
!       As 3 but also overrides Cumulus diagnosis if 
!          ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
!       Note that here Ri accounts for gradient adjustment by the 
!       non-local scheme.

!     Options for surface exchange
INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

!     Options for convective boundary layers
INTEGER, PARAMETER ::  UM_std     = 0
INTEGER, PARAMETER ::  neut_cbl   = 1
INTEGER, PARAMETER ::  LEM_conven = 2
INTEGER, PARAMETER ::  LEM_std    = 3

!     Options for stable boundary layers
INTEGER, PARAMETER ::  Long_tails           = 0
INTEGER, PARAMETER ::  Sharpest             = 1
INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
INTEGER, PARAMETER ::  Mes_tails            = 3
INTEGER, PARAMETER ::  Louis_tails          = 4
INTEGER, PARAMETER ::  Depth_based          = 5
INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
INTEGER, PARAMETER ::  LEM_stability        = 7
INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8
INTEGER, PARAMETER ::  Equilibrium_SBL     = 9

!     Options for Prandtl number (in local Ri scheme)
INTEGER, PARAMETER ::  Constant_SBL = 0
INTEGER, PARAMETER ::  LockMailhot2004 = 1  
  
!     Options for Keep_Ri_FA
INTEGER, PARAMETER :: except_disc_inv = 2! Only reduce local K to zero at 
                                         ! discontinuous inversions
!     Options for local_fa
INTEGER, PARAMETER :: to_sharp_across_1km = 1
INTEGER, PARAMETER :: ntml_level_corrn    = 2
INTEGER, PARAMETER :: free_trop_layers    = 3

!     Options for Kprof_cu
INTEGER, PARAMETER :: klcl_entr = 1      ! BL entrainment parametrization
  
INTEGER :: Kprof_cu  = off   
                       ! switch for non-local mixing across the LCL in Cu
                       ! Set when smag-1d blending selected
                       ! klcl_entr (=1) => set K-profile at the LCL using 
                       !     the standard CBL entrianment parametrization
!=======================================================================
! Integer values set from RUN_BL
!=======================================================================
INTEGER :: cbl_op = imdi ! was 0
                       ! Options for convective BL stability functions
                       ! UM_std   (=0) => Standard UM
                       ! LEM_conv (=1) => Conventional LEM
                       ! neut_conv(=2) => Keep neutral stability for Ri<0

INTEGER :: Variable_RiC = imdi ! was 0
                       ! Switch to allow different critical Richardson
                       !   numbers in the diagnosis of BL depth:
                       ! 0 => RiC=1 everywhere
                       ! 1 => RiC=0.25 for SHARPEST
                       !          1 otherwise

INTEGER :: idyndiag = imdi ! was 0
                       ! Switch to use zi/L in the diagnosis of
                       !   shear-driven boundary layers:
                       ! OFF (=0) => not used
                       ! DynDiag_ZL (=1) => over sea uses zi/L but
                       !   gives problems when BL_LEVELS is high
                       ! DynDiag_ZL_corrn (=2) => as 1 but copes with
                       !   high BL_LEVELS
                       ! DynDiag_ZL_CuOnly (=3) => as 2 but only applied 
                       !   to points diagnosed with Cumulus and strictly 
                       !   for sea points (fland<0.01, cf 0.5)
                       ! DynDiag_Ribased (=4) => as 3 but also overrides 
                       !   Cumulus diagnosis if 
                       !     ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
                       !   Note that here Ri accounts for gradient 
                       !   adjustment by the non-local scheme.

INTEGER :: ISrfExCnvGust = imdi ! was 0
                       ! Switch to include the effect of convective
                       ! downdraughts on surface exchange
                       ! OFF (=0) => not used: only boundary-layer
                       !   gustiness is considered (original version)
                       ! IP_SrfExWithCnv (=1) the impact of gustiness
                       !   due to boundary layer eddies is reduced
                       !   relative to the above, but eddies driven
                       !   by convective downdraughts are included

INTEGER :: entr_enhance_by_cu = imdi ! was 0
                       ! For 8C scheme ONLY: switch to enhance
                       !   entrainment in decoupled stratocu over cu
                       ! OFF (=0) => not used
                       ! Buoyrev_feedback (=1) => buoyancy reversal
                       !   feedback enhanced by Cu for 0<D<0.1

INTEGER :: relax_sc_over_cu = imdi ! was 0
                       ! For 8C scheme ONLY: switch to relax the 
                       !   requirements for diagnosis of Sc over Cu
                       !   to be for all cumulus, not just shallow

INTEGER :: fric_heating = imdi ! was 0
                       ! Switch to apply heating source from turbulence
                       ! dissipation
                       ! OFF (=0) => not used
                       ! ON  (=1) => used
INTEGER :: subs_couple_fix = imdi ! was 0
                       ! 1 - suppresses subsidence coupling with 
                       !     entrainment unless w well-behaved
INTEGER :: sg_orog_mixing = imdi ! was 0
                       ! SBL mixing dependent on sg orography
                       !   1 - extending length of SHARPEST tail
                       !       following McCabe and Brown (2007)
                       !   2 - include subgrid drainage shear in Ri and 
                       !       diffusion coefficients, 
                       !   3 - as 2 + orographically enhanced lambda_m 
                       !       (note lambda_h enhancement not included 
                       !       in bdy_expl2 in error but now operational 
                       !       in UKV) and smooth decrease to lambda_min
                       !       above
INTEGER :: entr_smooth_dec = off
                       ! Method of doing entrainment at dsc top
                       ! 0 - old method, include surface terms only
                       ! in weakly coupled, no cumulus situations
                       ! 1 - taper method - no hard limit but reduce
                       ! surface terms according to svl difference
INTEGER :: decfix= on  ! switch for correction to decoupling diagnosis
                       ! ROSE: removed from NL
INTEGER :: stopwe_sbl= on 
                       ! switch for spurious entrainment in SBLs
                       ! ROSE: removed from NL
INTEGER :: trweights1= imdi ! Was "on"
                       ! switch to use implicit weights of 1 for tracers (if on)

INTEGER :: flux_grad= imdi ! Was "Locketal2000" (0)
                       ! switch for revised flux-gradient relationships
INTEGER :: ng_stress= BrownGrant97_limited ! Was "off"
                       ! switch for non-gradient stress (ROSE: removed from NL)
INTEGER :: sbl_op= imdi ! was "Long_tails"
                       !  stable boundary layer option
INTEGER :: ishear_bl = on
                       ! switch for shear-dominated b.l.
                       ! ROSE: removed from NL, set in readlsta
INTEGER :: formdrag= imdi  ! Was Explicit_stress
                       ! switch for orographic form drag

INTEGER :: cor_mo_iter= imdi  ! was off (0)
                       ! switch for MO iteration correction
INTEGER :: non_local_bl= on  ! ROSE: removed from NL
                       ! switch on the non-local scheme
INTEGER :: local_fa= imdi  ! was off (0)
                       ! switch for free atmospheric mixing options
                       ! to_sharp_across_1km (=1) => smoothly switch to
                       !       sharpest across 1km AGL
                       ! ntml_level_corrn (=2) => correct level for 
                       !       ntml_local (downwards)
                       ! free_trop_layers (=3) => as "ntml_level_corrn" but 
                       !       also diagnose FA turbulent layer depths
INTEGER :: Prandtl= imdi  ! was "Constant_SBL"
                       ! switch for Prandtl number options
INTEGER :: Buddy_sea= imdi ! was "off"
                       ! switch to use the wind speed from
                       ! adjacent sea points for the
                       ! sea part of coastal grid points
INTEGER :: FD_stab_dep= imdi  ! Was on (1)
                       ! switch to implement stability
                       ! dependence of orographic form drag
INTEGER :: Keep_Ri_FA= on ! ROSE: was off, removed from NL
                       ! switch to keep local mixing in free atmosphere
INTEGER :: nl_bl_levels= imdi ! was "off"  
                       ! number of levels for non_local sheme
INTEGER :: ISeaZ0T= imdi ! was "Fixed_Z0T"
                       ! option for specifying the thermal roughness length
                       !       at sea points

!=======================================================================
! Real values set from RUN_BL
!=======================================================================
REAL :: WeightLouisToLong = rmdi ! was 0.0
                       ! Weighting of the Louis tail towards long tails:
                       ! 0.0 = Louis tails
                       ! 1.0 = Long tails

REAL :: zhloc_depth_fac = rmdi ! was 0.5
                       ! For idyndiag options DynDiag_Ribased and 
                       ! DynDiag_RiGAbased: the fractional height into 
                       ! the cloud layer reached by the local Ri-based 
                       ! BL depth calculation
REAL :: lambda_min_nml = 40.0
                       ! Minimum value of the mixing length (m)
REAL :: lambda_max_nml = 500.0
                       ! Maximum value of the mixing length (m)
REAL :: lambda_fac = 0.15
                       ! fraction of BL depth for mixing length
!=======================================================================
! Real parameters
!=======================================================================
REAL :: max_stress_grad = 0.05
                       ! Maximum implied stress gradient across the
                       ! boundary layer, used to limit the explicit
                       ! stress applied in non-local scalings (m/s2)

REAL :: dec_thres_cloud = rmdi ! was 0.05
                       ! Decoupling threshold for cloudy boundary
                       ! layers (larger makes decoupling less likely)

REAL :: max_cu_depth = 250.0
                       ! Max height above LCL for K profile in cumulus

REAL :: t_drain = 1800.0
                       ! timescale for drainage flow development (s)

REAL :: h_scale = 1500.0  
                       ! horizontal scale for drainage flows (m)
                       ! Currently hard-wired for UKV run over UK

REAL :: pr_max  = 5.0
                       ! Maximum allowed value of the Prandtl number with 
                       ! the LockMailhot2004 option

REAL :: beta_bl = 0.15
REAL :: beta_fa = 1.0
!                      ! Parameters governing the speed of transition from 
!                      ! 1D BL to Smagorinsky.  
!                      ! beta = 0.15 matches Honnert et al well in the BL 
!                      ! but a faster transition seems appropriate in the 
!                      ! free atmosphere above.
REAL :: RiCrit_sharp = 0.25
                       ! Critical Ri when SHARPEST stability functions used

REAL :: a_grad_adj = 3.26
                       ! parameter used in gradient adjustment

REAL :: max_t_grad = 1.0e-3
                       ! parameter used in gradient adjustment
                       
REAL :: Muw_SBL= 1.5   ! ROSE: removed from NL
REAL :: Mwt_SBL= 1.0   ! ROSE: removed from NL
                       ! Powers to use in prescription of
                       !    equilibrium profiles of stress and
                       !    buoyancy flux in Equilib. SBL model

REAL :: Puns = rmdi ! was 0.5 
REAL :: Pstb = rmdi ! was 2.0  
                    ! parameters for uncond stable numerical solver
                    ! Puns : used in an unstable BL column
                    ! Pstb : used in an stable BL column 
                       
REAL :: orog_drag_param= rmdi ! Was 0.3  
                       ! Drag coefficient for orographic form drag   
                       
REAL :: SeaSalinityFactor= rmdi ! was 1.0
                       ! Scaling of qsat allowing for salinity of sea water                       
 
REAL :: charnock = rmdi
INTEGER :: alpha_cd_batches = imdi
INTEGER :: alpha_Cd_items(max_number_alpha_cds) = imdi
REAL :: alpha_Cd_vals (max_number_alpha_cds) = rmdi 
REAL :: alpha_Cd (max_number_alpha_cds) = rmdi

! Maximum and minimum values for the STPH_RP scheme
! Boundary Layer
! ROSE: these items moved to stochastic_physics_run_mod
                                                     
!=======================================================================
!LOGICALs
!=======================================================================
LOGICAL :: l_use_bl_diag_term= .FALSE.  ! ROSE: removed from NL
      
LOGICAL :: L_us_blsol= .FALSE. ! Switch for stable and non-oscillatory
                               !  BL vertical diffusion scheme (was unset)

LOGICAL :: L_SBLco = .TRUE.   ! Switch for coupled gradient
                              ! method in Equilibrium SBL model
                              ! ROSE: removed from NL 

LOGICAL :: l_lambdam2     = .FALSE. ! LambdaM=2*LambdaH (operational setting)

LOGICAL :: l_full_lambdas = .FALSE. ! Lambdas NOT reduced above NTML_LOCAL+1

! logical for whether to skip calculations based on start of timestep
! quantities when using semi-lagrangian cycling with Endgame
! This is hard-wired to true in readlsta_4a so Endgame always uses it,
! otherwise it is false
! N.B. results should bit-compare whether this logical is true or false
! and so changing it will be a good test of whether new code has been
! added correctly
LOGICAL :: l_quick_ap2 = .FALSE.

!=======================================================================
!run_bl namelist
!=======================================================================    

NAMELIST/RUN_BL/ formdrag, orog_drag_param                                &
     ,sbl_op, cbl_op, cor_mo_iter                                         &
     ,alpha_cd_batches, alpha_cd_items, alpha_cd_vals                     &
     ,trweights1, flux_grad , nl_bl_levels, charnock                      &
     ,l_lambdam2, l_full_lambdas, local_fa, prandtl, seasalinityfactor    &
     ,iseaz0t, buddy_sea, l_us_blsol, puns, pstb                          &
     ,fd_stab_dep                                                         &
     ,weightlouistolong, variable_ric, idyndiag                           &
     ,zhloc_depth_fac, entr_smooth_dec                                    &
     ,entr_enhance_by_cu, dec_thres_cloud, isrfexcnvgust                  &
     ,fric_heating, sg_orog_mixing, subs_couple_fix                       &
     ,bdy_tke, tke_dlen                                                   &
     ,tke_cm_mx, tke_cm_fa ,my_lowest_pd_surf                             &
     ,l_my_condense, l_shcu_buoy, l_adv_turb_field                        &
     ,l_my_prod_adj, my_prod_adj_fact, tke_levels                         &
     ,l_my_initialize, l_my_ini_zero, l_local_above_tkelvs                &
     ,l_print_max_tke, my_ini_dbdz_min                                    &
     ,my_z_limit_elb, wb_ng_max, shcu_levels, relax_sc_over_cu            &
     ,calc_prob_of_vis

CONTAINS

SUBROUTINE check_run_bl()

! Description:
!   Subroutine to apply logic checks and set control variables based on the 
!   options selected in the run_bl namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umsections_mod, ONLY: atmos_sr
USE ereport_mod,    ONLY: ereport

IMPLICIT NONE

INTEGER                       :: icode         ! used for ereport
CHARACTER (LEN=80)            :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'check_run_bl'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_BL',zhook_in,zhook_handle)

! Set switch for turning on boundary layer code, default false
IF (atmos_sr(3) == '1A' .OR.    &
    atmos_sr(3) == '9B' .OR.    & 
    atmos_sr(3) == '9C') THEN
  l_bl  = .TRUE.
END IF

IF (.NOT. l_bl) THEN
  WRITE (cmessage,'(A44)') 'Boundary Layer code not used, l_bl = .FALSE.'
  icode = -100
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (lhook) CALL dr_hook('CHECK_RUN_BL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_bl

END MODULE bl_option_mod

