! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!+ Data module for switches for UM inputs, nlstcatm namelist
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

MODULE um_input_control_mod

USE var_input_mod,    ONLY: l_regular
USE missing_data_mod, ONLY: imdi, rmdi
USE rad_input_mod,    ONLY:                                                 &
     l_use_seasalt_indirect,         l_use_biogenic,                        &
     l_use_seasalt_direct,           l_use_dust,                            &
     l_use_sulpc_indirect_sw,        l_use_sulpc_indirect_lw,               &
     l_use_sulpc_direct,             l_use_ocff_direct,                     &
     l_use_ocff_indirect,            l_use_nitrate_direct,                  &
     l_use_nitrate_indirect,         l_use_soot_direct,                     &
     l_use_soot_indirect,            l_use_bmass_direct,                    &
     l_use_bmass_indirect,           l_use_arclbiom,                        &
     l_use_arclblck,                 l_use_arclsslt,                        &
     l_use_arclsulp,                 l_use_arcldust,                        &
     l_use_arclocff,                 l_use_arcldlta,                        &
     l_use_aod

IMPLICIT NONE

! Use the same number of cloud droplets for 1st and 2nd indirect effects
LOGICAL :: l_consistent_cdnc = .FALSE.

! Use sulphate aerosol in microphysics
LOGICAL :: l_use_sulphate_autoconv

! Use sea-salt aerosol in microphysics
LOGICAL :: l_use_seasalt_autoconv

! Use soot aerosol in microphysics
LOGICAL :: l_use_soot_autoconv

! Use biomass aerosol in microphysics
LOGICAL :: l_use_bmass_autoconv

! Use fossil-fuel organic carbon in microphysics
LOGICAL :: l_use_ocff_autoconv

! Use ammonium nitrate aerosol in microphysics
LOGICAL :: l_use_nitrate_autoconv

! Use autoconversion de-biasing scheme in microphysics
LOGICAL :: l_auto_debias
! Use sulphate aerosol no. in S-cycle
LOGICAL :: l_use_sulphate_sulpc

! Use sea-salt aerosol no. in S-cycle
LOGICAL :: l_use_seasalt_sulpc

! Use soot aerosol no. in S-cycle
LOGICAL :: l_use_soot_sulpc

! Use biomass aerosol no. in S-cycle
LOGICAL :: l_use_bmass_sulpc

! Use fossil-organic carbon aerosol no. in S-cycle
LOGICAL :: l_use_ocff_sulpc

! Use ammonium nitrate aerosol no. in S-cycle
LOGICAL :: l_use_nitrate_sulpc

! Use sea-salt aerosol no. in PM diagnostics
LOGICAL :: l_use_seasalt_pm

!---  Tracers ---
! Aerosol

LOGICAL :: l_dust_div1_lbc, l_dust_div1_lbc_out ! Dust active
LOGICAL :: l_dust_div2_lbc, l_dust_div2_lbc_out
LOGICAL :: l_dust_div3_lbc, l_dust_div3_lbc_out
LOGICAL :: l_dust_div4_lbc, l_dust_div4_lbc_out
LOGICAL :: l_dust_div5_lbc, l_dust_div5_lbc_out
LOGICAL :: l_dust_div6_lbc, l_dust_div6_lbc_out
LOGICAL :: l_so2, l_so2_lbc, l_so2_lbc_out           ! so2 active
LOGICAL :: l_dms, l_dms_lbc, l_dms_lbc_out           ! dms active 
LOGICAL :: l_so4_aitken, l_so4_aitken_lbc            ! so4_aitken 
LOGICAL :: l_so4_aitken_lbc_out                      ! active
LOGICAL :: l_so4_accu, l_so4_accu_lbc, l_so4_accu_lbc_out  ! so4_accu active
LOGICAL :: l_so4_diss, l_so4_diss_lbc, l_so4_diss_lbc_out  ! so4_diss active
LOGICAL :: l_nh3, l_nh3_lbc, l_nh3_lbc_out           ! nh3 active
LOGICAL :: l_soot_new, l_soot_new_lbc, l_soot_new_lbc_out        ! soot active
LOGICAL :: l_soot_agd, l_soot_agd_lbc, l_soot_agd_lbc_out
LOGICAL :: l_soot_cld, l_soot_cld_lbc, l_soot_cld_lbc_out
LOGICAL :: l_bmass_new, l_bmass_new_lbc, l_bmass_new_lbc_out   ! biomass active
LOGICAL :: l_bmass_agd, l_bmass_agd_lbc, l_bmass_agd_lbc_out 
LOGICAL :: l_bmass_cld, l_bmass_cld_lbc, l_bmass_cld_lbc_out
LOGICAL :: l_ocff_new, l_ocff_new_lbc,l_ocff_new_lbc_out    ! fossil fuel active
LOGICAL :: l_ocff_agd, l_ocff_agd_lbc, l_ocff_agd_lbc_out
LOGICAL :: l_ocff_cld, l_ocff_cld_lbc, l_ocff_cld_lbc_out
LOGICAL :: l_nitr_acc, l_nitr_acc_lbc, l_nitr_acc_lbc_out
LOGICAL :: l_nitr_diss, l_nitr_diss_lbc, l_nitr_diss_lbc_out

! For Aero_Ctl (Sulphur cycle or Soot)
INTEGER call_chem_freq     !No. times chem called per atm tstep

! Overall logical for CLASSIC aerosol scheme
! Note that this is currently used (and therefore set) 
! only locally in the boundary layer scheme, but would ideally 
! be used throughout the model
LOGICAL :: L_AERO_CLASSIC

! Sulphur cycle
LOGICAL :: l_sulpc_so2   ! S Cycle: SO2 MMR included
LOGICAL :: l_sulpc_dms   ! S Cycle: DMS MMR included
LOGICAL :: l_sulpc_ozone ! S Cycle: Ozone included for oxidation 
!          of DMS and SO2
LOGICAL :: l_sulpc_so2_o3_nonbuffered ! S Cycle: SO2+O3 reaction NOT
! buffered by NH3.
LOGICAL :: l_sulpc_online_oxidants ! Sulphur Cycle : Use online
! oxidants from UKCA
LOGICAL :: l_sulpc_2_way_coupling  ! Sulphur Cycle : Depleted oxidants
! are passed back to UKCA
LOGICAL :: l_so2_surfem  ! SO2 Surface Emissions
LOGICAL :: l_so2_hilem   ! SO2 High Level Emissions
LOGICAL :: l_so2_natem   ! SO2 Natural Emissions
LOGICAL :: l_dms_em      ! DMS Emissions
LOGICAL :: l_dms_em_inter      ! Interactive DMS Emissions
LOGICAL :: l_dms_ointer        ! DMS emissions from ocean model
LOGICAL :: l_dms_liss_merlivat ! Switches to determine which
LOGICAL :: l_dms_wanninkhof    !    scheme to use for interactive
LOGICAL :: l_dms_nightingale   !    sea-air exchange of DMS
LOGICAL :: l_sulpc_nh3   ! S Cycle: NH3 tracer included
LOGICAL :: l_nh3_em      ! S Cycle: NH3 emiss included

! Soot cycle

LOGICAL :: l_soot                ! Soot included
LOGICAL :: l_soot_surem          ! surface Soot emiss included
LOGICAL :: l_soot_hilem          ! elevated Soot emiss included

! Biomass aerosol

LOGICAL :: l_biomass             ! Biomass aerosol included
LOGICAL :: l_bmass_surem         ! Sfc biomass emiss included
LOGICAL :: l_bmass_hilem         ! Elevated bmass emiss included

! Fossil-fuel organic carbon aerosol

LOGICAL :: l_ocff                ! OCFF aerosol included
LOGICAL :: l_ocff_surem          ! Surface OCFF emissions included
LOGICAL :: l_ocff_hilem          ! Elevated OCFF emiss included

! Ammonium nitrate aerosol

LOGICAL :: l_nitrate             ! Ammonium nitrate aerosol included


! Dust deposition for ocean biology
LOGICAL :: l_dust2ocn        ! pass dust dep fields to the ocean

LOGICAL :: l_q10                  ! control T fn for soil resp

! MOSES II and Triffid logicals--------------------

LOGICAL :: l_veg_fracs          ! Switch for vegetation fractions
LOGICAL :: l_triffid            ! Switch for interactive veg model
LOGICAL :: l_phenol             ! Switch for leaf phenology

! Switch for running TRIFFID in equilibrium mode
LOGICAL :: l_trif_eq

! Switch for starting NRUN mid-way through a TRIFFID calling
! period
LOGICAL :: l_nrun_mid_trif
LOGICAL :: l_disturb      ! Switch for anthropogenic disturbance

! Vegetation:

! Update frequency for leaf phenology (days)
INTEGER :: phenol_period

INTEGER :: triffid_period ! Update frequency for TRIFFID (days)

! Hardwire until re-assessment of whether these need to be
! re-introduced as UMUI-set switches.
! RR old switches needed for addressing but should be defunct?
! RR - can be set with parameters if needed in the interim. ->




!===========================================================================
! LOGICAL - Misc items set from nlstcatm namelist
! need to be moved to relevant section
!===========================================================================

! PMSL smoothing
! True:  Use SOR algorithm
! False: Use Jacobi algorithm
LOGICAL :: l_pmsl_sor   = .TRUE. 

! Switch for interpolated winds in lbcs
! True for advecting winds interpolated in boundary zone
LOGICAL :: l_int_uvw_lbc = .FALSE.

! Boundary layer tracer mixing
LOGICAL :: l_bl_tracer_mix = .FALSE.

! Carbon cycle
! Include surface emissions
LOGICAL :: l_co2_emits        = .FALSE.  

! Hydrology:
LOGICAL :: l_hydrology  = .FALSE.  ! F: Turns off hydrology code

! Methane oxidation
LOGICAL :: l_use_methox = .FALSE.

!===========================================================================
! INTEGER - Misc items set from nlstcatm namelist
! need to be moved to relevant section
!===========================================================================

! Type of problem to be solved
INTEGER :: problem_number = IMDI

! Specification of CO2 absorption
! 1 => Simple method with fixed value.
! 2 => Complex method allowing linear and/or exponential variation.
! 3 => From the interactive carbon cycle.
INTEGER :: i_co2_opt = IMDI

! Time at which data assimilation starts (Hours after Basis Time)
INTEGER :: a_assim_start_min = IMDI

! Time at which data assimilation  ends
INTEGER :: a_assim_end_min   = IMDI

!===========================================================================
! REAL - Misc items set from nlstcatm namelist
! need to be moved to relevant section
!===========================================================================

! PMSL diagnostic 
! Orographic height threshold for new pmsl calculation
! from ATMOS_STASH_Misc in UMUI for Calc_NPMSL routine
REAL    :: npmsl_height = RMDI

!===========================================================================
! LOGICAL Top Level items set from nlstcatm namelist
!===========================================================================

! 360-day calendar. SAVE attribute required.
LOGICAL :: lcal360 = .FALSE. 

! Limited area model uses lateral boundary tendencies
LOGICAL :: l_lateral_boundary = .FALSE.

! OASIS coupling switch
LOGICAL :: l_oasis = .FALSE. 

! OASIS includes iceberg calving ancillary data
LOGICAL :: l_oasis_icecalve = .FALSE. 

! Couple through master PE
LOGICAL :: l_couple_master = .FALSE. 

! Use mixing ratio in atmos_physics1
LOGICAL :: l_mr_physics1 = .FALSE.          

! Use mixing ratio in atmos_physics2
LOGICAL :: l_mr_physics2 = .FALSE.

! Select dynamical core
! TRUE  => ENDGame
! FALSE => New Dynamics
LOGICAL :: l_endgame  = .FALSE.

!===========================================================================
! INTEGER Top Level items set from nlstcatm namelist
!===========================================================================

! Domain of atmosphere model   
INTEGER :: model_domain = imdi  
! 1 => Global Model        
! 2 => Limited Area Model (Classic style) 
! 3 => Limited Area Model (Cyclic boundary conditions - EW only)
! 4 => Limited Area Model (Cyclic boundary conditions - EW and NS)
! 5 => Single Column Model
! 6 => Site Specific Forecast Model (SSFM) 

! Coupling frequency in hours
INTEGER :: oasis_couple_freq = imdi

!===========================================================================
! Top Level items not set in namelist
!===========================================================================

! Output iteration counts
LOGICAL :: l_icount = .FALSE.      

! Max. no. of code sections
INTEGER, PARAMETER :: maxsects = 99 

! Array of code section-versions
CHARACTER(LEN=3) h_sect(0:maxsects) 

! Size of super array holding all tracers
INTEGER ::   super_array_size
INTEGER ::   moisture_array_size

!===========================================================================
! Misc items not set in namelist
!===========================================================================

! Carbon cycle
! Interactive 3D CO2 field for use with carbon cycle model
! This is set by the check_nlstcatm routine based on the value of
! i_co2_opt
LOGICAL :: l_co2_interactive  = .FALSE.

LOGICAL, PARAMETER :: ltleads =.FALSE. ! Switch for Leads temperature.
                      ! If FALSE, they are assumed to be TFS
                      ! Else they are prognostic.
                      ! Default setting: leads temperatures are set to TFS

!===========================================================================
! Define the nlstcatm namelist
!===========================================================================

      NAMELIST/nlstcatm/                                                &
        model_domain,                                                   &
        l_oasis, oasis_couple_freq, l_couple_master, l_oasis_icecalve,  &
        problem_number,                                                 &
        a_assim_start_min, a_assim_end_min,                             &
        npmsl_height,l_pmsl_sor,                                        &
        l_bl_tracer_mix, l_int_uvw_lbc,                                 &
        l_dust_div1_lbc,l_dust_div1_lbc_out,                            &
        l_dust_div2_lbc,l_dust_div2_lbc_out,                            &
        l_dust_div3_lbc,l_dust_div3_lbc_out,                            &
        l_dust_div4_lbc,l_dust_div4_lbc_out,                            &
        l_dust_div5_lbc,l_dust_div5_lbc_out,                            &
        l_dust_div6_lbc,l_dust_div6_lbc_out,                            &
        l_so2,l_so2_lbc,l_so2_lbc_out, l_dms,l_dms_lbc,l_dms_lbc_out,   &
        l_so4_aitken, l_so4_aitken_lbc,l_so4_aitken_lbc_out,            &
        l_so4_accu, l_so4_accu_lbc, l_so4_accu_lbc_out,                 &
        l_so4_diss, l_so4_diss_lbc, l_so4_diss_lbc_out,                 &
        l_nh3, l_nh3_lbc, l_nh3_lbc_out,                                &
        l_soot_new, l_soot_new_lbc, l_soot_new_lbc_out,                 &
        l_soot_agd, l_soot_agd_lbc, l_soot_agd_lbc_out,                 &
        l_soot_cld, l_soot_cld_lbc, l_soot_cld_lbc_out,                 &
        l_bmass_new, l_bmass_new_lbc, l_bmass_new_lbc_out,              &
        l_bmass_agd, l_bmass_agd_lbc, l_bmass_agd_lbc_out,              &
        l_bmass_cld, l_bmass_cld_lbc, l_bmass_cld_lbc_out,              &
        l_ocff_new, l_ocff_new_lbc, l_ocff_new_lbc_out,                 &
        l_ocff_agd, l_ocff_agd_lbc, l_ocff_agd_lbc_out,                 &
        l_ocff_cld, l_ocff_cld_lbc, l_ocff_cld_lbc_out,                 &
        l_nitr_acc, l_nitr_acc_lbc, l_nitr_acc_lbc_out,                 &
        l_nitr_diss, l_nitr_diss_lbc, l_nitr_diss_lbc_out,              &
        l_sulpc_so2,l_sulpc_dms,l_sulpc_ozone,                          &
        l_sulpc_so2_o3_nonbuffered, l_so2_surfem, l_so2_hilem,          &
        l_so2_natem,                                                    &
        l_sulpc_online_oxidants, l_sulpc_2_way_coupling,                &
        l_dms_em, l_dms_em_inter,                                       &
        l_dms_ointer,                                                   &
        l_dms_liss_merlivat, l_dms_wanninkhof, l_dms_nightingale,       &
        l_sulpc_nh3, l_nh3_em, l_soot, l_soot_surem, l_soot_hilem,      &
        l_biomass, l_bmass_surem, l_bmass_hilem,                        &
        l_dust2ocn,                                                     &
        i_co2_opt, l_co2_emits,                                         &
        l_q10, l_veg_fracs, l_triffid, l_phenol,                        &
        l_trif_eq, l_nrun_mid_trif, l_disturb,                          &
        phenol_period, triffid_period,                                  &
        l_use_seasalt_indirect, l_use_biogenic,                         &
        l_use_seasalt_direct, l_use_dust, l_use_sulpc_indirect_sw,      &
        l_use_sulpc_indirect_lw, l_use_sulpc_direct,                    &
        l_use_sulphate_autoconv, l_use_seasalt_autoconv, l_auto_debias, &
        l_use_sulphate_sulpc, l_use_seasalt_sulpc,                      &
        l_ocff, l_ocff_surem, l_ocff_hilem, l_use_ocff_autoconv,        &
        l_use_ocff_sulpc, l_use_ocff_direct, l_use_ocff_indirect,       &
        l_nitrate, l_use_nitrate_direct, l_use_nitrate_indirect,        &
        l_use_nitrate_autoconv, l_use_nitrate_sulpc, l_use_seasalt_pm,  &
        l_consistent_cdnc,                                              &
        call_chem_freq, l_use_soot_direct, l_use_soot_indirect,         &
        l_use_soot_autoconv, l_use_soot_sulpc, l_use_bmass_direct,      &
        l_use_bmass_indirect, l_use_bmass_autoconv, l_use_bmass_sulpc,  &
        l_use_arclbiom, l_use_arclblck,  l_use_arclsslt,                &
        l_use_arclsulp, l_use_arcldust,  l_use_arclocff, l_use_arcldlta,&
        l_use_methox, l_mr_physics1,l_mr_physics2,                      &
        l_hydrology,                                                    &
        l_use_aod,                                                      &
        lcal360,                                                        &
        l_endgame, l_lateral_boundary, l_regular

CONTAINS

SUBROUTINE check_nlstcatm()

! Description:
!   Subroutine to apply logic controls and set control variables based on the 
!   options selected in the nlstcatm namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_NLSTCATM',zhook_in,zhook_handle)

! Set logicals based on integer choice in namelist
IF (i_co2_opt == 3) THEN
  l_co2_interactive = .TRUE.
END IF

IF (lhook) CALL dr_hook('CHECK_NLSTCATM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nlstcatm

END MODULE um_input_control_mod
