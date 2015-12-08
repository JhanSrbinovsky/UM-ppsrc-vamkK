! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Perform a 1-timestep integration of the Atmosphere Model
!
! Subroutine Interface:

SUBROUTINE atm_step_4A(                                               &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGBND Control data calculated from NAMELIST-
        NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
        FIXHD_BOUNDA,INTHD_BOUNDA,LOOKUP_BOUNDA,LOOKUP_COMP_BOUNDA,     &
        REALHD_BOUNDA,                                                  &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ENDGame prognostic variables
exner,                                                                &
! River routing
! ARGATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of ARGAOCPL.h)
!
        XPA, XUA, XVA, YPA, YUA, YVA,                                   &
! END ARGATCPL
 g_p_field,                                                           &
 g_r_field,                                                           &
 obs_flag,obs,obs_flag_len,obs_len,                                   &
 ngrgas,grgas_addr)

USE sl_tracer1_mod
USE tr_reset_mod
USE TR_Set_Phys_mod
USE unpack_tracer_mod
USE atm_fields_bounds_mod
USE rad_input_mod
USE switches, ONLY: l_ctile
USE cv_cntl_mod, ONLY: lcv_3d_ccw
USE mphys_constants_mod, ONLY: ntot_land, ntot_sea
USE run_aerosol_mod
USE level_heights_mod
USE trignometric_mod
USE dyn_coriolis_mod
USE dyn_var_res_mod
USE diff_coeff_mod
USE turb_diff_mod
USE metric_terms_mod
USE rad_mask_trop_mod
USE rot_coeff_mod
USE swapable_field_mod, ONLY : swapable_field_pointer_type
USE timestep_mod
USE atm_step_local
! Logicals and stash diagnostic types needed for SKEB2
USE stochastic_physics_run_mod,  ONLY:                                &
    l_skeb2, l_rp2, rhcrit_max, rhcrit_min, m_ci, m_ci_max, m_ci_min, &
    par_mezcla_max, par_mezcla_min, par_mezcla, &
    g0_rp, g0_rp_max, g0_rp_min, charnock_max, charnock_min, &
    ricrit_rp, ricrit_rp_max, ricrit_rp_min,    &
    lambda_min_rp, lambda_min_rp_max, lambda_min_rp_min,              &
    a_ent_1_rp, a_ent_1_rp_max, a_ent_1_rp_min, &
    g1_rp, g1_rp_max, g1_rp_min, gwd_frc_max, gwd_frc_min, &
    kay_gwave_max, kay_gwave_min
    
USE ancil_info, ONLY: nsmax
USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
USE iau_mod, ONLY : l_iau, iau_firstcallts, iau_lastcallts

USE earth_constants_mod, ONLY: g, earth_radius

USE g_wave_input_mod
USE bl_option_mod

USE Submodel_Mod
USE nstypes, ONLY: ntype

USE conversions_mod
USE atmos_constants_mod
USE water_constants_mod
USE ukca_radaer_struct_mod
USE ukca_cdnc_mod, ONLY: ukca_cdnc_struct, ukca_cdnc,             &
                   ukca_cdnc_get, cdnc_dim1, cdnc_dim2, cdnc_dim3
USE ukca_option_mod, ONLY: L_ukca_chem, L_ukca_useumuivals,       &
                     L_ukca_set_trace_gases, L_ukca_strat,        &
                     L_ukca_strattrop, L_ukca_prescribech4,       &
                     L_ukca_trop, L_ukca_tropisop, L_ukca_raq,    &
                     l_ukca, l_ukca_aie1, l_ukca_aie2,            &
                     l_ukca_radaer

!      USE ancil_info, ONLY: nsmax

USE lbc_read_data_mod, ONLY: rimweightsa
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE eg_destroy_vert_damp_mod, ONLY : eg_destroy_vert_damp
USE eg_destroy_horz_drag_mod, ONLY : eg_destroy_horz_drag
USE eg_dry_static_adj_mod, ONLY : eg_dry_static_adj
USE um_types
USE eg_helmholtz_mod
USE proc_info_mod,            ONLY: n_proc,me

USE eg_idl_forcing_mod
USE eg_r_mod
USE eg_r_s_mod
USE eg_f1sp_mod
USE eg_f1sp_inc_mod
USE eg_thetav_theta_mod
USE eg_info_file_mod
USE eg_sisl_resetcon_mod
USE eg_set_helm_lhs_mod
USE eg_sisl_init_mod
USE eg_sisl_setcon_mod
USE eg_q_to_rhox_mod
USE eg_q_to_mix_mod
USE eg_sl_helmholtz_mod
USE eg_sl_helmholtz_inc_mod
USE eg_sl_full_wind_mod
USE eg_sl_rho_mod
USE eg_sl_thermo_mod
USE eg_sl_moisture_mod
USE eg_idl_set_init_mod
USE eg_coriolis_star_mod
USE eg_dx_diags_mod
USE eg_set_helmholtz_mod
USE eg_idl_ni_init_mod
USE eg_mix_to_q_mod
USE eg_v_at_poles_mod
USE integrity_mod
USE eg_NI_filter_Ctl_mod
USE eg_NI_filter_incs_Ctl_mod
USE eg_dxout_mod
USE horiz_grid_mod
USE ref_pro_mod

USE eg_balance_lbc_values_mod

USE eg_alpha_mod
USE eg_alpha_ramp_mod
USE set_metric_terms_4A_mod

USE Control_Max_Sizes
USE UM_ParVars
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE set_star_zero_level_mod

USE departure_pts_mod
USE fields_rhs_mod
USE helmholtz_const_matrix_mod
USE coriolis_mod
USE gravity_mod
USE eg_parameters_mod

USE print_diag_mod

USE rimtypes
USE lbc_mod
USE eg_conserv_moist_mod
USE eg_conserv_tracers_mod
USE eg_swap_bounds_mod

USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon, io3_trop_map,  &
                      io3_trop_map_masscon
USE problem_mod, ONLY: standard, monsoon, dynamical_core, idealised_problem,&
                       standard_namelist
USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv, tp_bv_isoth, &
                        tp_dyn_core, tp_dyn_core_lam, tp_namelist, tp_dump
USE eng_mass_diag_mod, ONLY : eng_mass_diag
USE diagnostics_adv_mod
USE wet_to_dry_n_calc_mod

USE eg_set_adv_winds_mod
! Stochastic Physics
USE stph_rp2_mod,   ONLY: stph_rp2
USE stph_setup_mod, ONLY: stph_setup
USE stph_skeb2_mod, ONLY: stph_skeb2

USE eg_total_mass_mod
USE dynamics_input_mod, ONLY:                                     &
          L_mix_ratio,L_fix_mass,                                 &
          L_eg_dry_static_adj, tol_sc_fact,                       &
          T_surf_ref,p_surf_ref,NumCycles,L_new_tdisc,            &
          L_LBC_balance,L_regular,L_qwaterload,                   &
          GCR_use_residual_Tol,GCR_diagnostics,GCR_precon_option, &
          GCR_tol_res,GCR_tol_abs,L_lbc_old,                      &
          eg_vert_damp_coeff,L_RK_dps,L_inc_solver, Ih,           &
          L_conserv_smooth_lap,L_eliminate_rho,L_init_Fnm1,       &
          eg_vert_damp_profile, eta_s,                            &
          l_expl_horz_drag,l_impl_horz_drag,                      &
          innits,l_accel_convergence,n_rims_to_do,                &
          L_filter_dump
          
USE dynamics_testing_mod, ONLY:                                   &
          L_Physics,L_Run_With_Physics2,                          &
          L_Backwards, L_dry, L_idealised_data  


USE sl_input_mod, ONLY:                                           &
          L_conserve_tracers,halo_lam,halo_phi,look_lam,look_phi, &
          recip_dlam,recip_dphi,L_conserv,L_mono,L_high,          &
          L_Ritchie_high,L_Ritchie_mono,L_2d_sl_geometry,         &
          L_sl_halo_reprod,high_order_scheme,monotone_scheme,     &
          Ritchie_high_order_scheme,Ritchie_monotone_scheme,      &
          Depart_scheme,Depart_order,interp_vertical_search_tol,  &
          Theta_SL,moist_SL,Wind_SL,rho_SL,l_2dcomm

! Aerosols
USE aero_ctl_mod_4A,          ONLY: aero_ctl_4A
USE get_sulpc_oxidants_mod,   ONLY: get_sulpc_oxidants
USE set_arcl_clim_mod,        ONLY: set_arcl_clim
USE set_arcl_dimensions_mod,  ONLY: set_arcl_dimensions
USE write_sulpc_oxidants_mod, ONLY: write_sulpc_oxidants

USE dust_parameters_mod, ONLY:                                       &
     l_dust,            l_dust_div1,      l_dust_div2,               &
     l_dust_div3,       l_dust_div4,      l_dust_div5,               &
     l_dust_div6
USE um_input_control_mod,  ONLY:                                     &
     model_domain,      problem_number,   l_int_uvw_lbc,             &
     l_dust_div1_lbc,   l_dust_div2_lbc,  l_dust_div3_lbc,           &
     l_dust_div4_lbc,   l_dust_div5_lbc,  l_dust_div6_lbc,           &
     l_so2,             l_so2_lbc,                                   &
     l_dms,             l_dms_lbc,        l_so4_aitken,              &
     l_so4_aitken_lbc,  l_so4_accu,       l_so4_accu_lbc,            &
     l_so4_diss,        l_so4_diss_lbc,   l_nh3,                     &
     l_nh3_lbc,         l_soot_new,       l_soot_new_lbc,            &
     l_soot_agd,        l_soot_agd_lbc,   l_soot_cld,                &
     l_soot_cld_lbc,    l_bmass_new,      l_bmass_new_lbc,           &
     l_bmass_agd,       l_bmass_agd_lbc,  l_bmass_cld,               &
     l_bmass_cld_lbc,   l_ocff_new,       l_ocff_new_lbc,            &
     l_ocff_agd,        l_ocff_agd_lbc,   l_ocff_cld,                &
     l_ocff_cld_lbc,    l_nitr_acc,       l_nitr_acc_lbc,            &
     l_nitr_diss,       l_nitr_diss_lbc,  l_sulpc_so2,               &
     l_sulpc_dms,       l_sulpc_ozone,    l_sulpc_so2_o3_nonbuffered,&
     l_so2_surfem,      l_so2_hilem,      l_so2_natem,               &
     l_dms_em,          l_dms_em_inter,   l_sulpc_2_way_coupling,    &
     l_dms_ointer,      l_dms_wanninkhof, l_sulpc_online_oxidants,   &
     l_dms_nightingale, l_sulpc_nh3,      l_dms_liss_merlivat,       &
     l_nh3_em,          l_soot,           l_soot_surem,              &
     l_soot_hilem,      l_biomass,        l_bmass_surem,             &
     l_bmass_hilem,     l_co2_interactive,l_use_seasalt_autoconv,    &
     l_ocff,            l_ocff_surem,     l_use_seasalt_sulpc,       &
     l_ocff_hilem,      l_use_ocff_sulpc, l_use_sulphate_sulpc,      &
     l_nitrate,         l_use_seasalt_pm, l_use_nitrate_sulpc,       &
     call_chem_freq,    l_use_soot_sulpc, l_use_bmass_sulpc,         &
     l_mr_physics1,     l_mr_physics2,                               &
     lcal360,           moisture_array_size, super_array_size
     
USE domain_params

USE mphys_inputs_mod, ONLY:                                          &
     l_mcr_qcf2,        l_mcr_qcf2_lbc,   l_mcr_qgraup,              &
     l_mcr_qgraup_lbc,  l_mcr_qrain_lbc,  l_mcr_qrain

USE cloud_inputs_mod, ONLY:                                          &
     l_cld_area,        l_acf_cusack,     l_pc2,                     &
     l_pc2_reset,       l_pc2_lbc,        l_acf_brooks,              &
     rhcrit

USE river_inputs_mod, ONLY: river_vel, river_mcoef, i_river_vn

USE eng_corr_inputs_mod, ONLY:                                       &
     l_emcorr,          lmass_corr,       lemq_print,                &
     a_energysteps,     lqt_corr,         lenergy

USE nudging_input_mod, ONLY: l_nudging

USE murk_inputs_mod, ONLY: l_murk, l_murk_advect, l_murk_lbc

USE cosp_input_mod, ONLY: l_cosp


USE nlstcall_mod, ONLY : ldump, &
                         lexit, &
                         ltimer


use cable_data_mod, only : cable_atm_step
USE chsunits_mod, ONLY : nunits

USE science_fixes_mod, ONLY : l_use_old_mass_fix

IMPLICIT NONE

!
! Description: Perform a 1-timestep integration of the Atmosphere Model,
!   including assimilation, physics and dynamics processing.
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! Declarations:

! Global variables (*CALLed COMDECKs etc...):
! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=1

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end
! TYPDUMA needs TYPSIZE included first
!L --------------- Dump headers (atmosphere)-------------
      INTEGER :: A_FIXHD(LEN_FIXHD)    ! fixed length header
      INTEGER :: A_INTHD(A_LEN_INTHD)  ! integer header
      INTEGER :: A_CFI1(A_LEN_CFI1+1)  ! compress field index
      INTEGER :: A_CFI2(A_LEN_CFI2+1)  ! compress field index
      INTEGER :: A_CFI3(A_LEN_CFI3+1)  ! compress field index

      REAL::A_REALHD(A_LEN_REALHD)                    ! real header
      REAL::A_LEVDEPC(A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1)! level  dep const
      REAL::A_ROWDEPC(A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1)! row    dep const
      REAL::A_COLDEPC(A_LEN1_COLDEPC*A_LEN2_COLDEPC+1)! column dep const
      REAL::A_FLDDEPC(A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1)! field  dep const
      REAL::A_EXTCNST(A_LEN_EXTCNST+1)                ! extra constants
      REAL::A_DUMPHIST(LEN_DUMPHIST+1)                ! temp hist file

      ! Meaningful parameter names for integer constants header:
! ----------------------- include file: IHEADAPM -----------------------
! Description: Meaningful parameter names to index A_INTHD array in
!              atmosphere dump, ie INTEGER CONSTANTS, and reduce magic
!              numbers in code.
!
      INTEGER,PARAMETER:: ih_a_step          = 1  ! Timestep no.
      INTEGER,PARAMETER:: ih_rowlength       = 6  ! No. of points E-W
      INTEGER,PARAMETER:: ih_rows            = 7  ! No. of points N-S

      ! No. of model levels (0=surface)
      INTEGER,PARAMETER:: ih_model_levels    = 8

      ! No. of model levels with moisture
      INTEGER,PARAMETER:: ih_wet_levels      = 9

      ! No. of deep soil temperature levels
      INTEGER,PARAMETER:: ih_soilT_levels    = 10

      INTEGER,PARAMETER:: ih_cloud_levels    = 11 ! No. of cloud levels
      INTEGER,PARAMETER:: ih_tracer_levels   = 12 ! No. of tracer levels

      ! No. of boundary layer levels
      INTEGER,PARAMETER:: ih_boundary_levels = 13
      INTEGER,PARAMETER:: ih_N_types         = 15 ! No. of field types

       ! Height generation method
      INTEGER,PARAMETER:: ih_height_gen      = 17

      ! First rho level at which height is constant
      INTEGER,PARAMETER:: ih_1_c_rho_level   = 24

      INTEGER,PARAMETER:: ih_land_points     = 25 ! No. of land points
      INTEGER,PARAMETER:: ih_ozone_levels    = 26 ! No. of ozone levels

      ! No. of deep soil moisture levels
      INTEGER,PARAMETER:: ih_soilQ_levels    = 28

      ! Number of convective cloud levels
      INTEGER,PARAMETER:: ih_convect_levels  = 34
      INTEGER,PARAMETER:: ih_rad_step        = 35 ! Radiation timestep
      INTEGER,PARAMETER:: ih_AMIP_flag       = 36 ! Flag for AMIP run
      INTEGER,PARAMETER:: ih_AMIP_year       = 37 ! First AMIP year
      INTEGER,PARAMETER:: ih_AMIP_month      = 38 ! First AMIP month
      INTEGER,PARAMETER:: ih_AMIP_day        = 49 ! First AMIP day
      INTEGER,PARAMETER:: ih_ozone_month     = 40 ! Current ozone month
      INTEGER,PARAMETER:: ih_SH_zonal_quad   = 41 ! L_SH_zonal_quadratics
      INTEGER,PARAMETER:: ih_SH_zonal_begin  = 42 ! SH_zonal_begin
      INTEGER,PARAMETER:: ih_SH_zonal_period = 43 ! SH_zonal_period
      INTEGER,PARAMETER:: ih_SH_level_weight = 44 ! SuHe_level_weight
      INTEGER,PARAMETER:: ih_SH_sigma_cutoff = 45 ! SuHe_sigma_cutoff
      INTEGER,PARAMETER:: ih_friction_time   = 46 ! frictional_timescale

! IHEADAPM end
      ! Meaningful parameter names for real constants header:
! ----------------------- include file: RHEADAPM -----------------------
! Description: Meaningful parameter names to index A_REALHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.

      ! East-West   grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaEW         = 1

      ! North-South grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaNS         = 2

      ! Latitude  of first p point in degrees
      INTEGER,PARAMETER:: rh_baselat         = 3

      ! Longitude of first p point in degrees
      INTEGER,PARAMETER:: rh_baselong        = 4

      ! Latitude  of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlat          = 5

      ! Longitude of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlong         = 6

      ! Height of top theta level (m)
      INTEGER,PARAMETER:: rh_z_top_theta     =16

      ! total moisture of the atmosphere
      INTEGER,PARAMETER:: rh_tot_m_init      =18

      ! total mass of atmosphere
      INTEGER,PARAMETER:: rh_tot_mass_init   =19

      ! total energy of atmosphere
      INTEGER,PARAMETER:: rh_tot_energy_init =20

      ! energy correction = energy drift
      INTEGER,PARAMETER:: rh_energy_corr     =21

! RHEADAPM end
      ! Meaningful parameter names for fixed header:
! ----------------------- include file: FHEADAPM -----------------------
! Description: Meaningful parameter names to index A_FIXHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.
 
      ! Start of Row Dependent Constant
      INTEGER,PARAMETER:: fh_RowDepCStart   = 115

      ! Start of Col Dependent Constant
      INTEGER,PARAMETER:: fh_ColDepCStart   = 120

! FHEADAPM end
      ! PP headers

      INTEGER :: A_LOOKUP(LEN1_LOOKUP,A_LEN2_LOOKUP) ! lookup heads
      INTEGER :: A_MPP_LOOKUP(MPP_LEN1_LOOKUP,A_LEN2_LOOKUP)
      INTEGER :: a_ixsts(len_a_ixsts)     ! stash index array

      REAL    :: a_spsts(len_a_spsts)     ! atmos stash array
! TYPDUMA end
!-----------------Start of TYPCONA------------------------------------

! Constants for routines independent of resolution.
! Requires Use Control_Max_Sizes for MAX_REQ_THPV_LEVS in cconsts.h
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end

! typcona.h originally contained constants for the atmosphere.
! Almost all of these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, RAD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block did not correspond to constants specified in 
! argcona.h, so it remains here, even though argcona.h has been deleted.
! cderv_trig and CDERIVED in cconsts.h should be moved to modules too.

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1
      !First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150
      !Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
      !Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is -1.
      ! Similarly for A_UKCA_INDEX.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      ! A_TR_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: A_TR_StashItem(A_MAX_TRVARS)

      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)
      ! UKCA_tr_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: UKCA_tr_StashItem(A_MAX_UKCAVARS) 

      ! A_TR_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: A_TR_LBC_StashItem(A_MAX_TRVARS) 
      INTEGER :: A_TR_active_lbc_index(A_MAX_TRVARS) 

      ! UKCA_tr_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: UKCA_tr_LBC_StashItem(A_MAX_UKCAVARS) 
      INTEGER :: UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)

      COMMON/ATRACER/A_TR_INDEX, A_TR_StashItem,                        &
     &               A_TR_LBC_StashItem, A_TR_active_lbc_index,         &
     &               A_UKCA_INDEX, UKCA_tr_StashItem,                   &
     &               UKCA_tr_LBC_StashItem, UKCA_tr_active_lbc_index

! CTRACERA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start typ_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "typptra.h", and requires that "ctracera.h" is
!  also included.  
!
! This file belongs in section: Top Level
!
! Note: including this file requires lbc_mod.
!
! ENDGame

REAL ::                                                                 &
  exner_surf(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end),                            &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows),                                          &
  etadot(wdims_s%i_start:wdims_s%i_end,                                 &
         wdims_s%j_start:wdims_s%j_end,                                 &
         wdims_s%k_start:wdims_s%k_end),                                &
  m_v   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cl  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_r   (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_gr  (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end),                                &
  m_cf2 (qdims_s%i_start:qdims_s%i_end,                                 &
         qdims_s%j_start:qdims_s%j_end,                                 &
         qdims_s%k_start:qdims_s%k_end)
REAL :: thetav(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)
REAL :: dryrho(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,                           &
               pdims_s%k_start:pdims_s%k_end)

REAL ::  Exner      (pdims_s%i_start:pdims_s%i_end,                     &
                 pdims_s%j_start:pdims_s%j_end,                         &
                 pdims_s%k_start:pdims_s%k_end+1)

      ! 1.1: Data variables stored in primary space.
      real :: U    (udims_s%i_start:udims_s%i_end,                       &
                    udims_s%j_start:udims_s%j_end,                       &
                    udims_s%k_start:udims_s%k_end)

      real :: V    (vdims_s%i_start:vdims_s%i_end,                       &
                    vdims_s%j_start:vdims_s%j_end,                       &
                    vdims_s%k_start:vdims_s%k_end)
      real :: E_TRB  (tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: TSQ_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: QSQ_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: COV_TRB(tdims_l%i_start:tdims_l%i_end,                     &
                      tdims_l%j_start:tdims_l%j_end,                     &
                      tdims_l%k_start:tdims_l%k_end)
      real :: ZHPAR_SHCU(row_length, rows)

      real :: W    (wdims_s%i_start:wdims_s%i_end,                       &
                    wdims_s%j_start:wdims_s%j_end,                       &
                    wdims_s%k_start:wdims_s%k_end)

      real :: RHO  (pdims_s%i_start:pdims_s%i_end,                       &
                    pdims_s%j_start:pdims_s%j_end,                       &
                    pdims_s%k_start:pdims_s%k_end)

      real :: THETA(tdims_s%i_start:tdims_s%i_end,                       &
                    tdims_s%j_start:tdims_s%j_end,                       &
                    tdims_s%k_start:tdims_s%k_end)

      real :: Q     (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCL   (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCF   (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QCF2  (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QRAIN (qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

      real :: QGRAUP(qdims_l%i_start:qdims_l%i_end,                      &
                     qdims_l%j_start:qdims_l%j_end,                      &
                     qdims_l%k_start:qdims_l%k_end)

! EXNER_RHO_LEVELS is still 1D
! Exner pressure on rho levels
      real :: EXNER_RHO_LEVELS(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))   

! The next few fields are targets as they are used in multi-variate 
! swap_bounds calls

      real, target ::                                                    &
            U_ADV(udims_l%i_start:udims_l%i_end,                         &
                  udims_l%j_start:udims_l%j_end,                         &
                  udims_l%k_start:udims_l%k_end)

      real, target ::                                                    &
            V_ADV(vdims_l%i_start:vdims_l%i_end,                         &
                  vdims_l%j_start:vdims_l%j_end,                         &
                  vdims_l%k_start:vdims_l%k_end)
      real, target ::                                                    &
            W_ADV(wdims_l%i_start:wdims_l%i_end,                         &
                  wdims_l%j_start:wdims_l%j_end,                         &
                  wdims_l%k_start:wdims_l%k_end)

      ! 1.2: Data variables stored in secondary space.
      ! P is still 1D
      real :: P( (pdims_s%i_end-pdims_s%i_start+1)* &
                 (pdims_s%j_end-pdims_s%j_start+1)* &
                 (pdims_s%k_end+1) )

      ! Pressure on theta levels
      real :: P_THETA_LEVELS(tdims_s%i_start:tdims_s%i_end,              &
                             tdims_s%j_start:tdims_s%j_end,              &
                             tdims_s%k_start:tdims_s%k_end)

      ! Exner pressure on theta levels
      real :: EXNER_THETA_LEVELS(tdims_s%i_start:tdims_s%i_end,          &
                                 tdims_s%j_start:tdims_s%j_end,          &
                                 tdims_s%k_start:tdims_s%k_end)

      ! 1.3: Cloud Fields
      real :: CCW_RAD(qdims%i_start:qdims%i_end,                         &
                      qdims%j_start:qdims%j_end,                         &
                                    qdims%k_end)

      ! CCAs size is dependant on L_3D_CCA
      ! N_CCA_LEV will be set to the correct value (either wet_levels or 1)

      real :: CCA(qdims%i_start:qdims%i_end,                             &
                  qdims%j_start:qdims%j_end,                             &
                                  n_cca_lev)
      real :: CF_AREA  (qdims%i_start:qdims%i_end,                       &
                        qdims%j_start:qdims%j_end,                       &
                                      qdims%k_end)
      real :: CF_BULK  (qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID(qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN(qdims_l%i_start:qdims_l%i_end,                   &
                        qdims_l%j_start:qdims_l%j_end,                   &
                        qdims_l%k_start:qdims_l%k_end)

      ! 1.4: Soil Ancillary fields
      real :: DEEP_SOIL_TEMP(land_field,sm_levels)
      real :: SMCL(land_field,sm_levels)
      real :: STHU(land_field,sm_levels)
      real :: STHF(land_field,sm_levels)

      ! 1.5: Radiation Increments
      real :: SW_INCS(rdims2%i_start:rdims2%i_end,                       &
                      rdims2%j_start:rdims2%j_end,                       &
                      rdims2%k_start:rdims2%k_end)   ! SW radiation increments
      real :: LW_INCS(tdims%i_start:tdims%i_end,                         &
                      tdims%j_start:tdims%j_end,                         &
                      tdims%k_start:tdims%k_end)   ! LW radiation increments

      ! PAR radiation increment
      real :: DIRPAR(row_length,rows)

      ! 1.6: Ozone
      real :: O3(o3dims%i_start:o3dims%i_end)
!  tropopause-based ozone
      real :: TPPSOZONE(o3dims%i_start:o3dims%i_end)
!  Cariolle ozone tracer variables
      real :: OZONE_TRACER(tdims_s%i_start:tdims_s%i_end,                &
                           tdims_s%j_start:tdims_s%j_end,                &
                           tdims_s%k_start:tdims_s%k_end)
      real :: O3_PROD_LOSS(1,rows,model_levels)
      real :: O3_P_L_VMR  (1,rows,model_levels)
      real :: O3_VMR      (1,rows,model_levels)
      real :: O3_P_L_TEMP (1,rows,model_levels)
      real :: O3_TEMP     (1,rows,model_levels)
      real :: O3_P_L_COLO3(1,rows,model_levels)
      real :: O3_COLO3    (1,rows,model_levels)
      
! TRACERS are still 1D
      ! Tracer and aerosol fields
      ! TRACERS are dealt w/ differently
      ! these are the maximum sizes:
      real :: TRACER              (tr_levels*theta_off_size*(a_tracer_last-a_tracer_first))
      real :: TRACER_UKCA         ((tdims%k_end-tdims%k_start+1)*theta_off_size*tr_ukca)
      real :: TRACER_LBC          (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_LBC_TEND     (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_UKCA_LBC     (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_ukca)
      real :: TRACER_UKCA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_ukca)

      real :: MURK_SOURCE(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)

      real :: MURK       (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV1  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV2  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV3  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV4  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV5  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DUST_DIV6  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO2        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: DMS        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_AITKEN (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_ACCU   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO4_DISS   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: H2O2       (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: NH3        (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SOOT_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_NEW  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_AGD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: BMASS_CLD  (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_NEW   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_AGD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: OCFF_CLD   (tdims_s%i_start:tdims_s%i_end,                 &
                          tdims_s%j_start:tdims_s%j_end,                 &
                          tdims_s%k_start:tdims_s%k_end)

      real :: SO2_NATEM  (tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)

      real :: OH        (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: HO2       (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: H2O2_LIMIT(tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: O3_CHEM   (tdims%i_start:tdims%i_end,                      &
                         tdims%j_start:tdims%j_end,                      &
                         tdims%k_start:tdims%k_end)

      real :: CO2       (tdims_s%i_start:tdims_s%i_end,                  &
                         tdims_s%j_start:tdims_s%j_end,                  &
                         tdims_s%k_start:tdims_s%k_end)

! USER_MULT<N> are still 1D
! 1.8: Multi-level user ancillary fields

      real :: USER_MULT1 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT2 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT3 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT4 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT5 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT6 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT7 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT8 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT9 (oneddims%i_start:oneddims%i_end)
      real :: USER_MULT10(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT11(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT12(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT13(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT14(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT15(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT16(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT17(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT18(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT19(oneddims%i_start:oneddims%i_end)
      real :: USER_MULT20(oneddims%i_start:oneddims%i_end)


      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      real ::      OROG_LBC(LENRIMA(fld_type_p,halo_type_extended,1))
      real ::         U_LBC(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::         V_LBC(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::         W_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::       RHO_LBC(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end)
      real ::     THETA_LBC(LENRIMA(fld_type_p,halo_type_extended,1),tdims_s%k_start:tdims_s%k_end)
      real ::         Q_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCL_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCF_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::      QCF2_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     QRAIN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::    QGRAUP_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::   CF_BULK_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     EXNER_LBC(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end+1)
      real ::     U_ADV_LBC(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::     V_ADV_LBC(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::     W_ADV_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::      MURK_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV1_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV2_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV3_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV4_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV5_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV6_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       SO2_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       DMS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: SO4_AITKEN_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_ACCU_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_DISS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::        NH3_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_NEW_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_AGD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_CLD_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   NITR_ACC_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  NITR_DISS_LBC(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)

      ! Lateral Boundary Condition tendencies
      real ::         U_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::         V_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::         W_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::       RHO_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end)
      real ::     THETA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),tdims_s%k_start:tdims_s%k_end)
      real ::         Q_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCL_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::       QCF_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::      QCF2_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     QRAIN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::    QGRAUP_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::   CF_BULK_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_LIQUID_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real :: CF_FROZEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),qdims_l%k_start:qdims_l%k_end)
      real ::     EXNER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),pdims_s%k_start:pdims_s%k_end+1)
      real ::     U_ADV_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)
      real ::     V_ADV_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)
      real ::     W_ADV_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)
      real ::      MURK_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV1_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV2_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV3_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV4_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV5_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: DUST_DIV6_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       SO2_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::       DMS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real :: SO4_AITKEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_ACCU_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SO4_DISS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::        NH3_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   SOOT_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  BMASS_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_NEW_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_AGD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   OCFF_CLD_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::   NITR_ACC_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)
      real ::  NITR_DISS_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),tdims_s%k_start:tdims_s%k_end)

      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
      real :: TSTAR     (row_length,rows)
      logical :: LAND      (row_length,rows)
      real :: TSTAR_ANOM(row_length,rows)

      ! 2.15: Fields for coastal tiling
      real :: FRAC_LAND (land_field)
      real :: TSTAR_LAND(row_length,rows)
      real :: TSTAR_SEA (row_length,rows)
      ! these fields are targets - see note below relating to TI etc.
      real, target :: TSTAR_SICE(row_length,rows,1)
      real, target :: TSTAR_SICE_CAT(row_length,rows,nice_use)

      ! Set pointers for sea-ice and land albedos
      real :: SICE_ALB(row_length,rows)
      real :: LAND_ALB(row_length,rows)

      ! 2.2: Data variables stored in secondary space.

      real :: PSTAR(row_length,rows)

      ! 2.3: Cloud fields
      real ::      LCBASE(row_length,rows)
      real ::         CCB(row_length,rows)
      real ::         CCT(row_length,rows)
      real ::       CCLWP(row_length,rows)
      real ::   DEEP_FLAG(row_length,rows)
      real :: PAST_PRECIP(row_length,rows)
      real :: PAST_CONV_HT(row_length,rows)

      ! 2.4: Boundary layer fields

      real :: ZH(row_length,rows)
 
!     Convective downdraught mass-flux at cloud base 
      REAL :: ddmfx(row_length, rows) 

      ! Standard deviation of turbulent fluctuations of layer 1 temperature
      real :: T1_SD(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      real :: Q1_SD(row_length,rows)

      ! Decoupled screen-level temperatures
      REAL ::  TScrnDcl_SSI(row_length,rows)
      REAL :: TScrnDcl_TILE(land_field,ntiles)
      REAL ::     tStbTrans(row_length,rows)

      ! Number of model levels in the  turbulently mixed layer
      INTEGER :: NTML(row_length,rows)

      ! Top level for turb mixing in any decoupled Sc layer
      INTEGER :: NTDSC(row_length,rows)

      ! Bottom level for turb mixing in any decoupled Sc layer
      INTEGER :: NBDSC(row_length,rows)

      LOGICAL :: CUMULUS(row_length,rows)

      ! 2.4: Soil Ancillary fields

      real :: SAT_SOILW_SUCTION(land_field)
      real :: THERM_CAP        (land_field)
      real :: THERM_COND       (land_field)
      real :: VOL_SMC_CRIT     (land_field)
      real :: VOL_SMC_WILT     (land_field)
      real :: VOL_SMC_SAT      (land_field)
      real :: SAT_SOIL_COND    (land_field)
      real :: CLAPP_HORN       (land_field)

      ! 2.5: Other surface fields
      real :: CANOPY_WATER(land_field)
      real :: Z0          (row_length,rows)
      real :: GS          (land_field)

      ! 2.6: Orographic Ancillary fields

      real :: OROGRAPHY   (row_length,rows)
      real :: OROG_SD     (land_field)
      real :: OROG_SIL    (land_field)
      real :: OROG_HO2    (land_field)
      real :: OROG_GRAD_X (land_field)
      real :: OROG_GRAD_Y (land_field)
      real :: OROG_UNFILT (land_field)
      real :: OROG_GRAD_XX(land_field)
      real :: OROG_GRAD_XY(land_field)
      real :: OROG_GRAD_YY(land_field)

      ! 2.7: Sea/Sea Ice

      real :: U_SEA(row_length,rows)
      real :: V_SEA(row_length,n_rows)
      real :: U_0_P(row_length,rows)
      real :: V_0_P(row_length,rows)

! these fields are targets b/c there are pointers w/in ATM_STEP
! which point to one or another of them depending on the ice category selection
      real, target :: TI            (row_length,rows,1)
      real, target :: ICE_FRACTION  (row_length,rows,1)
      real, target :: ICE_THICKNESS (row_length,rows,1)
      real, target :: TI_CAT        (row_length,rows,nice) 
      real, target :: ICE_FRACT_CAT (row_length,rows,nice)
      real, target :: ICE_THICK_CAT (row_length,rows,nice)
      real, target :: SNODEP_SEA    (row_length,rows,1)
      real, target :: SNODEP_SEA_CAT(row_length,rows,nice_use)

! Effective conductivity is not a target as field only exists in D1 if nice>1
      real :: ICE_K_CAT     (row_length,rows,nice)

      ! 2.8: Snow

      real :: SNODEP         (row_length,rows)
      real :: CATCH_SNOW     (land_field,ntiles)
      real :: SNOW_GRND      (land_field,ntiles)

      ! SNSOOT may not be used as of vn6.6
      real :: SNSOOT(row_length,rows)  ! Snow soot content

      ! 2.9: aerosol emission fields, including mineral dust parent soil props

      real :: SOIL_CLAY (row_length,rows)
      real :: SOIL_SILT (row_length,rows)
      real :: SOIL_SAND (row_length,rows)
      real :: DUST_MREL1(row_length,rows)
      real :: DUST_MREL2(row_length,rows)
      real :: DUST_MREL3(row_length,rows)
      real :: DUST_MREL4(row_length,rows)
      real :: DUST_MREL5(row_length,rows)
      real :: DUST_MREL6(row_length,rows)

      real :: SO2_EM     (row_length,rows)
      real :: DMS_EM     (row_length,rows)
      real :: SO2_HILEM  (row_length,rows)
      real :: NH3_EM     (row_length,rows)
      real :: SOOT_EM    (row_length,rows)
      real :: SOOT_HILEM (row_length,rows)
      real :: BMASS_EM   (row_length,rows)
      real :: BMASS_HILEM(row_length,rows)
      real :: OCFF_EM    (row_length,rows)
      real :: OCFF_HILEM (row_length,rows)
      real :: DMS_CONC   (row_length,rows)
      real :: DMS_OFLUX  (row_length,rows)

! USER_ANC<N> fields are still 1D
      ! 2.10: User ancillary fields
      real :: USER_ANC1(row_length*rows*1)
      real :: USER_ANC2(row_length*rows*1)
      real :: USER_ANC3(row_length*rows*1)
      real :: USER_ANC4(row_length*rows*1)
      real :: USER_ANC5(row_length*rows*1)
      real :: USER_ANC6(row_length*rows*1)
      real :: USER_ANC7(row_length*rows*1)
      real :: USER_ANC8(row_length*rows*1)
      real :: USER_ANC9(row_length*rows*1)
      real :: USER_ANC10(row_length*rows*1)
      real :: USER_ANC11(row_length*rows*1)
      real :: USER_ANC12(row_length*rows*1)
      real :: USER_ANC13(row_length*rows*1)
      real :: USER_ANC14(row_length*rows*1)
      real :: USER_ANC15(row_length*rows*1)
      real :: USER_ANC16(row_length*rows*1)
      real :: USER_ANC17(row_length*rows*1)
      real :: USER_ANC18(row_length*rows*1)
      real :: USER_ANC19(row_length*rows*1)
      real :: USER_ANC20(row_length*rows*1)

      !   2.11: Store arrays for energy correction calculation
      real :: NET_FLUX (row_length,rows)
      real :: NET_MFLUX(row_length,rows)

      !   2.12: Tiled Vegetation and Triffid fields
      real :: FRAC_TYP (land_field)
      real :: FRAC_CON1(land_field)  ! fraction of broadleaf tree
      real :: FRAC_CON2(land_field)  ! fraction of needleleaf tree
      real :: FRAC_CON3(land_field)  ! fraction of C3 grass
      real :: FRAC_CON4(land_field)  ! fraction of C4 grass
      real :: FRAC_CON5(land_field)  ! fraction of shrub
      real :: FRAC_CON6(land_field)  ! fraction of urban
      real :: FRAC_CON7(land_field)  ! fraction of water
      real :: FRAC_CON8(land_field)  ! fraction of soil
      real :: FRAC_CON9(land_field)  ! fraction of ice
      real :: LAI_PFT  (land_field)
      real :: CANHT_PFT(land_field)
      real :: DISTURB_VEG(land_field)
      REAL :: soil_alb(land_field)   ! albedo of underlying soil
      REAL :: obs_alb_sw(land_field) ! Observed snow-free SW albedo 
      REAL :: obs_alb_vis(land_field)! Observed snow-free VIS albedo 
      REAL :: obs_alb_nir(land_field)! Observed snow-free NIR albedo 
      real, target :: SOIL_CARB(land_field)
      real, target :: SOIL_CARB1(land_field)
      real, target :: SOIL_CARB2(land_field)
      real, target :: SOIL_CARB3(land_field)
      real, target :: SOIL_CARB4(land_field)
      real :: NPP_PFT_ACC(land_field)
      real :: G_LF_PFT_ACC(land_field)
      real :: G_PHLF_PFT_ACC(land_field)
      real :: RSP_W_PFT_ACC(land_field)
      real, target :: RSP_S_ACC(land_field)
      real, target :: RSP_S_ACC1(land_field)
      real, target :: RSP_S_ACC2(land_field)
      real, target :: RSP_S_ACC3(land_field)
      real, target :: RSP_S_ACC4(land_field)
      real :: CAN_WATER_TILE(land_field,ntiles)
      real :: CATCH_TILE(land_field,ntiles)
      real :: INFIL_TILE(land_field,ntiles)
      real :: RGRAIN_TILE(land_field,ntiles)
      real :: SNODEP_TILE(land_field,ntiles)
      real :: TSTAR_TILE(land_field,ntiles) 
      real :: Z0_TILE(land_field,ntiles)  
      real :: z0h_tile(land_field,ntiles)  
      real :: DOLR_FIELD(row_length,rows) 
      real :: LW_DOWN(row_length,rows) 
      real :: SW_TILE_RTS(land_field,ntiles)

! MORUSES - new urban two-tile scheme
      real :: hgt(land_field)        ! Building height
      real :: hwr(land_field)        ! Height to width
      real :: wrr(land_field)        ! Width ratio
      real :: disp(land_field)       ! Displacement height
      real :: ztm(land_field)        ! Effective roughness length for momentum
      real :: albwl(land_field)      ! Wall albedo
      real :: albrd(land_field)      ! Road albedo
      real :: emisw(land_field)      ! Wall emmissivity
      real :: emisr(land_field)      ! Road emmissivity

!! REMOVED SLAB AS PART OF VN7.0
!      !   2.13: Slab Model

!   2.14: Carbon cycle fields
      real :: CO2FLUX(row_length,rows)
      real :: CO2_EMITS(row_length,rows)

!   2.15: Fields carried forward from previous version.
!         May not be required
!      real, pointer :: SURF_RESIST_NIT(:)  ! Surface resistance on
!                                    ! non-ice tiles
!      real, pointer :: ROOT_DEPTH_NIT(:)   ! Root depth on non-ice tiles
!      real, pointer :: Z0V_TYP(:)          ! Vegetative Roughness length on
!                                    ! tiles
!      real, pointer :: ICE_EDGE(:)
!      real, pointer :: OROG_TENDENCY(:)    ! Orographic tendencies
!      real, pointer :: OROG_SD_TENDENCY(:) ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
!      real, pointer :: ETATHETA(:)
!      real, pointer :: ETARHO(:)
!      real, pointer :: RHCRIT(:)
!      real, pointer :: SOIL_THICKNESS(:)

      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      real :: zseak_theta(0:model_levels)
      real :: Ck_theta(0:model_levels)
      real :: zseak_rho(model_levels)
      real :: Ck_rho(model_levels)   

      ! 2.16: Fields for large-scale hydrology scheme.
      real :: TI_MEAN(land_field)
      real :: TI_SIG(land_field)
      real :: FEXP(land_field)
      real :: GAMMA_INT(land_field)
      real :: WATER_TABLE(land_field)
      real :: FSFC_SAT(land_field)
      real :: F_WETLAND(land_field)

      real :: STHZW(land_field)
      real :: A_FSAT(land_field)
      real :: C_FSAT(land_field)
      real :: A_FWET(land_field)
      real :: C_FWET(land_field)

      ! 2.17: Fields for River routing.
      real :: RIV_SEQUENCE (river_row_length,river_rows)
      real :: RIV_DIRECTION(river_row_length,river_rows)
      real :: RIV_STORAGE  (river_row_length,river_rows)
      real :: TOT_SURFROFF (land_field)
      real :: TOT_SUBROFF  (land_field)
      real :: RIV_INLANDATM(land_field)
! Fields for water conservation correction due to lake evaporation 
      real :: ACC_LAKE_EVAP(row_length,rows) 
      ! Fields for grid-to-grid river routing (river routing 2A)
      real :: RIV_IAREA    (row_length,rows)   ! Drainage area
      real :: RIV_SLOPE    (row_length,rows)   ! Grid-cell slope
      real :: RIV_FLOWOBS1 (row_length,rows)   ! Initial values of flow
      real :: RIV_INEXT    (row_length,rows)   ! Flow direction (x)
      real :: RIV_JNEXT    (row_length,rows)   ! Flow direction (y)
      real :: RIV_LAND     (row_length,rows)   ! Land-type (land/river/sea)
      real :: RIV_SUBSTORE (row_length,rows)   ! Subsurface storage
      real :: RIV_SURFSTORE(row_length,rows)   ! Surface storage
      real :: RIV_FLOWIN   (row_length,rows)   ! Surface inflow
      real :: RIV_BFLOWIN  (row_length,rows)   ! Subsurface inflow

! Fields used when coupling using OASIS.
      real :: C_SOLAR   (row_length,rows)      ! CPL solar radn
      real :: C_BLUE    (row_length,rows)      ! CPL blue radn
      real :: C_LONGWAVE(row_length,rows)      ! CPL lw radn
      real :: C_TAUX    (row_length,rows)      ! CPL taux 
      real :: C_TAUY    (row_length,rows)      ! CPL tauy 
      real :: C_W10     (row_length,rows)      ! CPL 10m wind   
      real :: C_SENSIBLE(row_length,rows)      ! CPL sensible ht flx
      real :: C_SUBLIM  (row_length,rows,nice_use)! CPL sublim rate
      real :: C_EVAP    (row_length,rows)      ! CPL Evap rate
! Next field is the conductive flux through the ice that is used to
! force the ice model.  If using multilayers in ice (l_sice_multilayers=T), 
! this is the surface downwards heat flux (surf_ht_flux_sice), otherwise
! this is the conductive heat flux through all the ice (sea_ice_htf) 
! previously known as 'botmelt'.
      real :: C_FCONDTOPN(row_length,rows,nice)! CPL Multi-cat fcondtop
      real :: C_TOPMELTN(row_length,rows,nice) ! CPL Multi-cat tmlt
      real :: C_LSRAIN  (row_length,rows)      ! CPL Lg scl rain rate
      real :: C_LSSNOW  (row_length,rows)      ! CPL Lg scl snow rate
      real :: C_CVRAIN  (row_length,rows)      ! CPL Cnvctv rain rate
      real :: C_CVSNOW  (row_length,rows)      ! CPL Cnvctv snur rate
      real :: C_RIVEROUT(row_length,rows)      ! CPL Riv outflow->ocn
      real :: C_CALVING (row_length,rows)      ! CPL Iceberg Calving->ocn

!   2.18: JULES variables
      REAL :: SNOWDEPTH(    land_field,ntiles)
                             ! Snow depth on ground on tiles (m)
      REAL :: RHO_SNOW_GRND(land_field,ntiles)
                             ! Snowpack bulk density (kg/m3)
      REAL :: NSNOW(        land_field,ntiles)
                             ! Number of snow layers on ground on tiles
                             ! NOTE that this is converted to an integer.
      REAL :: DS(        land_field,ntiles,nsmax)
                             ! Snow layer thickness (m)
      REAL :: SICE(      land_field,ntiles,nsmax)
                             ! Snow layer ice mass on tiles (Kg/m2)
      REAL :: SLIQ(      land_field,ntiles,nsmax)
                             ! Snow layer liquid mass on tiles (Kg/m2)
      REAL :: TSNOWLAYER(land_field,ntiles,nsmax)
                             ! Snow layer temperature (K)
      REAL :: RHO_SNOW(  land_field,ntiles,nsmax)
                             ! Snow layer densities (kg/m3)
      REAL :: RGRAINL(   land_field,ntiles,nsmax)
                             ! Snow layer grain size on tiles (microns)

!         FLake lake scheme variables
      REAL :: lake_depth( land_field)
      REAL :: lake_fetch( land_field)
      REAL :: lake_t_mean(land_field)
      REAL :: lake_t_mxl( land_field)
      REAL :: lake_t_ice( land_field)
      REAL :: lake_h_mxl( land_field)
      REAL :: lake_h_ice( land_field)
      REAL :: lake_shape( land_field)
      REAL :: lake_g_dt(  land_field)

      ! UKCA oxidant fields
      real :: OH_UKCA ( tdims%i_start:tdims%i_end,                       &
                        tdims%j_start:tdims%j_end,                       &
                        tdims%k_start:tdims%k_end)
      real :: HO2_UKCA( tdims%i_start:tdims%i_end,                       &
                        tdims%j_start:tdims%j_end,                       &
                        tdims%k_start:tdims%k_end)
      real :: H2O2_UKCA(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: O3_UKCA  (tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: HNO3_UKCA(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)

      ! Aerosol climatologies
      real :: ARCLBIOG_BG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBIOM_IC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBLCK_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLBLCK_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSSLT_FI(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSSLT_JT(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_AC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_AK(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLSULP_DI(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B1(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B2(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B3(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B4(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B5(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDUST_B6(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_FR(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_AG(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLOCFF_IC(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      real :: ARCLDLTA_DL(tdims%i_start:tdims%i_end,                     &
                          tdims%j_start:tdims%j_end,                     &
                          tdims%k_start:tdims%k_end)
      
      ! Ammonium nitrate aerosol
      real :: NITR_ACC (tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)
      real :: NITR_DISS(tdims_s%i_start:tdims_s%i_end,                   &
                        tdims_s%j_start:tdims_s%j_end,                   &
                        tdims_s%k_start:tdims_s%k_end)

! End typ_atm_fields.h
! TYPBND - needs TYPSIZE included first
      ! Headers from atmosphere boundary data sets
      ! Second index of header arrays = 1 Lateral boundary data
      ! = 2 Lower boundary data
      INTEGER :: FIXHD_BOUNDA(LEN_FIXHD,1)      ! Fixed header
      INTEGER :: INTHD_BOUNDA(A_LEN_INTHD,1)    ! Integer header

        ! Lookups for the first set of LBCs in the LBC file
      INTEGER :: LOOKUP_BOUNDA(LEN1_LOOKUP,RIM_LOOKUPSA)

      ! Varying items from the LOOKUP table for the entire LBC file
      INTEGER :: LOOKUP_COMP_BOUNDA(LEN1_LBC_COMP_LOOKUP,BOUND_LOOKUPSA)

      REAL ::  REALHD_BOUNDA(A_LEN_REALHD,1)   ! Real header

      !  Control data calculated from namelist
      INTEGER :: RIM_STEPSA      ! Set by IN_BOUND
      INTEGER :: NBOUND_LOOKUP(2)
      COMMON/CBND/                                                      &
       RIM_STEPSA
! TYPBND end
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

! River routing
! TYPATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of TYPAOCPL.h)
!
      REAL :: XPA(AOCPL_ROW_LENGTH+1)  ! Atmosphere TP longitude coordina
      REAL :: XUA(0:AOCPL_ROW_LENGTH)  ! Atmosphere U longitude coordinat
      REAL :: XVA(AOCPL_ROW_LENGTH+1)  ! Atmosphere V longitude coordinat
      REAL :: YPA(AOCPL_P_ROWS)        ! Atmosphere TP latitude coordinat
      REAL :: YUA(AOCPL_P_ROWS)        ! Atmosphere U latitude coordinate
      REAL :: YVA(0:AOCPL_P_ROWS)      ! Atmosphere V latitude coordinate
! TYPATCPL end
INTEGER                                                           &
  g_p_field                                                       &
                        ! IN : global horiz domain for atmos
, g_r_field             ! IN : global horiz domain for rivers
INTEGER :: obs_flag_len,obs_len
INTEGER :: obs_flag(obs_flag_len)
REAL    :: obs(obs_len)

!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"

! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! This file belongs in section: Top Level

!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:

!------------------   End of Physics   ---------------------------------

      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(model_levels_max) ! Timestep of max w
      INTEGER :: time_div_max(model_levels_max) ! Timestep of max div
      INTEGER :: time_div_min(model_levels_max) ! Timestep of min div
      INTEGER :: time_lapse_min(model_levels_max) ! Timestep of min
      INTEGER :: time_max_shear(model_levels_max) !Timestep max shear
      INTEGER :: time_max_wind(model_levels_max) ! Timestep of max wind
      INTEGER :: time_KE_max(model_levels_max) ! Timestep of max KE
      INTEGER :: time_KE_min(model_levels_max) ! Timestep of min KE
      INTEGER :: time_noise_max(model_levels_max) ! Timestep of max

      REAL:: frictional_timescale(model_levels_max) ! For idealised case
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(0:model_levels_max) ! Max w at a level
      REAL :: max_div_run(model_levels_max) ! Max divergence at a level
      REAL :: min_div_run(model_levels_max) ! Min divergence at a level
      REAL :: min_lapse_run(model_levels_max) ! Min dtheta/dz at a level
      REAL :: max_shear_run(model_levels_max) ! Max shear at a level
      REAL :: max_wind_run(model_levels_max) ! Max wind at a level
      REAL :: max_KE_run(model_levels_max)   ! Max KE at a level
      REAL :: min_KE_run(model_levels_max)   ! Min KE at a level
      REAL :: max_noise_run(model_levels_max) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA
!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
        rpemax, rpemin, ipesum, rpesum,                                 &
        max_w_run, min_theta1_run, dtheta1_run,                         &
        max_div_run, min_div_run, min_lapse_run,                        &
        max_shear_run, max_wind_run,                                    &
        max_noise_run, max_KE_run, min_KE_run,                          &
        time_KE_max, time_KE_min,                                       &
        time_w_max, time_div_max, time_div_min, time_lapse_min,         &
        time_max_shear, time_max_wind,                                  &
        time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(model_levels_max)
      REAL :: friction_level(model_levels_max)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                              &
       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
       SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
       base_frictional_timescale, SuHe_sigma_cutoff,                   &
       L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
       SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: vert_grid_ratio
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(model_levels_max)
      REAL :: rho_ref(model_levels_max)
      REAL :: exner_ref(model_levels_max + 1)
      REAL :: q_ref(model_levels_max)
      REAL :: u_ref(model_levels_max)
      REAL :: v_ref(model_levels_max)
      REAL :: z_orog_print(0:model_levels_max)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: qforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: uforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: vforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind

! ENDGAME
      REAL :: T_surface
      REAL :: Eccentricity
      ! Following two variables used only if L_rotate_grid=.true.
      REAL :: grid_NP_lon ! Longitude (degrees) of grid's north pole
      REAL :: grid_NP_lat ! Latitude (degrees) of grid's north pole
      REAL :: AA_jet_u0   ! See QJRMS 133,1605--1614
      REAL :: AA_jet_A    !
      REAL :: theta_pert
      REAL :: ring_height
      REAL :: angular_velocity ! Planet's angular velocity
      REAL :: T0_P, T0_E ! deep atmosphere baroclinic wave surface temperatures
      INTEGER :: Trefer_number
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: AA_jet_m  ! See QJRMS 133,1605--1614
      INTEGER :: AA_jet_n  !
      INTEGER :: chain_number ! Run continuation number
      LOGICAL :: L_rotate_grid    ! .true. for rotating North pole of grid
      LOGICAL :: L_baro_Perturbed ! Used for baroclinic test to specify
                                  ! pert or steady
      LOGICAL :: L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,  &
                 L_HeldSuarez2_drag,                                        &
                 L_baro_inst, L_isothermal, L_exact_profile, L_balanced,    &
                 L_solid_body
      LOGICAL :: L_deep_baro_inst ! deep atmosphere baroclinic wave switch          


      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type
      INTEGER :: b_const, k_const ! deep atmosphere baroclinic wave parameters

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                              &
       h_o, h_o_actual, h_o_per_step,                                  &
       lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
       Witch_power, plat_size_x, plat_size_y,                          &
       height_domain, delta_x, delta_y, big_factor, mag, vert_grid_ratio, &
       first_theta_height, thin_theta_height, p_surface,               &
       theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
       u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
       ujet_lat, ujet_width,                                           &
       t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
       u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
       z_orog_print, grow_steps,                                       &
       surface_type, grid_number, grid_flat,                           &
       tprofile_number, qprofile_number, uvprofile_number,             &
       num_profile_data, num_uvprofile_data,                           &
       tforce_option, qforce_option, uvforce_option,                   &
       num_tforce_levels, num_tforce_times,                            &
       num_qforce_levels, num_qforce_times,                            &
       num_uvforce_levels, num_uvforce_times,                          &
       L_pforce, pforce_option, num_pforce_times,                      &
       first_constant_r_rho_level_new,                                 &
       big_layers, transit_layers, mod_layers,                         &
       zprofile_data, tprofile_data, qprofile_data,                    &
       z_uvprofile_data, uprofile_data, vprofile_data,                 &
       tforce_time_interval, qforce_time_interval,                     &
       uvforce_time_interval, newtonian_timescale,                     &
       z_tforce_data, z_qforce_data, z_uvforce_data,                   &
       tforce_data, qforce_data, uforce_data, vforce_data,             &
       tforce_data_modlev, qforce_data_modlev,                         &
       uforce_data_modlev, vforce_data_modlev,                         &
       pforce_time_interval, p_surface_data,                           &
       L_initialise_data,                                              &
       L_perturb_t, perturb_magnitude_t,                               &
       L_perturb_q, perturb_magnitude_q,                               &
       L_perturb_correlate_tq,                                         &
       L_perturb_correlate_vert,                                       &
       L_perturb_correlate_time,                                       &
       perturb_type, perturb_height,                                   &
       L_constant_dz, L_polar_wind_zero,                               &
       L_wind_balance, L_rotate_winds,                                 &
       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
       L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
       L_pressure_balance, L_vert_Coriolis,                            &
       cool_rate, L_force, L_force_lbc,                                &
       zprofile_orog, idl_interp_option, hf,                           &
       L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
       L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
       idl_bubble_option, idl_bubble_max                               &
      , idl_bubble_height, idl_bubble_width, idl_bubble_depth          &
      , idl_bubble_xoffset,idl_bubble_yoffset                          &
      , L_idl_bubble_saturate,                                         &
       L_rotating, L_fixed_lbcs, L_code_test,                          &
       L_baroclinic, L_cyclone,                                        &
       L_damp, L_geo_for, L_bomex,                                     &
       DMPTIM, HDMP, ZDMP,                                             &
       u_geo, v_geo,                                                   &
!ENDGAME
       T_surface, chain_number,                                        &
       Trefer_number,                                                  &
       tstep_plot_frequency, tstep_plot_start, Eccentricity,           &
       L_rotate_grid, grid_NP_lon, grid_NP_lat,                        &
       AA_jet_m, AA_jet_n, AA_jet_u0, AA_jet_A, L_baro_Perturbed,      &
       L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,       &
       L_HeldSuarez2_drag,                                             &
       L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const, k_const,    &      
       ring_height, theta_pert, L_isothermal,                          &
       L_exact_profile, L_balanced, L_solid_body, angular_velocity
! CRUNTIMC end
! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10
! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end
!
! Description:
!   This comdeck declares an integer variable 'ModelType' whose value
!   determines whether a model run is global, limited area or zonal.
!   The values of ModelType associated with each of the run types are
!   defined by integer parameters which are also declared below.
!    ModelType is set in subroutine SETLOGIC.
!
      ! Value used to represent the global model configuration
      INTEGER,PARAMETER:: GlobalModel=1

      ! Value used to represent the limited area model configuration
      INTEGER,PARAMETER:: LimitedAreaModel=2

      ! Value used to represent the 'periodic in x' model config
      ! INTEGER,PARAMETER:: ZonalModel=2

! Global scalars:
      INTEGER     ModelType  ! Integer switch which is equated to one
                             ! of the above parameters in a model run,
                             ! and so determines the configuration

! COMMON blocks:
      COMMON /RunType/ ModelType

! C_GLOBAL end
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200


!
! ENDGame prognostic variables (not included in the start dump)
!

REAL rho_r_sq_n  (pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end,                        &
                  pdims_s%k_start:pdims_s%k_end)
REAL rho_r_sq_np1(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end,                        &
                  pdims_s%k_start:pdims_s%k_end)


REAL theta_star(tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end),                         &
!       temorarily stored atmos_phys2 initial state for ENDGame
!        source term computation in eg_R_S:
   theta_star_n(tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end)

REAL :: exner_prime_term(pdims_s%i_start:pdims_s%i_end,                 &
                         pdims_s%j_start:pdims_s%j_end,                 &
                         pdims_s%k_start:pdims_s%k_end)

REAL :: exner_0(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end+1)

REAL                                                                    &
      q_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    qcl_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    qcf_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
     cf_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    cfl_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    cff_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
      m_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    mcl_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
    mcf_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
   mcf2_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
  mrain_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end),                            &
 mgraup_star(qdims_s%i_start:qdims_s%i_end,                             &
             qdims_s%j_start:qdims_s%j_end,                             &
             qdims_s%k_start:qdims_s%k_end)

REAL                                                                    &
 exner_star(edims_s%i_start:edims_s%i_end,                              &
            edims_s%j_start:edims_s%j_end,                              &
            edims_s%k_start:edims_s%k_end),                             &
 frac_control(land_field,ntype)   !Forcing for land surface (3C)


REAL                                                                    &
  rho_n (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end)

REAL  biogenic(row_length, rows, model_levels)


! Subroutine arguments:
! 3-D fields of species to be passed down to radiation
INTEGER, INTENT(IN) :: ngrgas
INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! Local parameters:
CHARACTER(LEN=*) routinename
PARAMETER (   routinename='Atm_Step_4A')



INTEGER, SAVE :: dx_frame_cnt
REAL exner_res,exner_res_tmp

REAL           :: total_rho, mass_fix_factor,                           &
                  total_gr, total_gr_r, total_rho_r,                    &
                  mass_fix_A, mass_fix_B


INTEGER                                                                 &
 errorstatus      ! Return code : 0 Normal Exit : >0 Error

CHARACTER(LEN=256)                                                      &
 cmessage         ! Error message if return code >0


LOGICAL, PARAMETER :: l_test_tracer=.FALSE.

LOGICAL :: l_call_from_solver

! Used for checking if at end of integration
LOGICAL :: lexitnow

INTEGER, SAVE :: timestep_number_chainsafe

! Local arrays for using the aerosol climatology for NWP

! Internal model switches
Logical L_USE_ARCL(NPD_ARCL_SPECIES)


! Array index of each component
INTEGER i_arcl_compnts(npd_arcl_compnts)

INTEGER cloud_tol

REAL                                                              &
  mag_vector_np (model_levels)                                    &
, dir_vector_np (model_levels)                                    &
, mag_vector_sp (model_levels)                                    &
, dir_vector_sp (model_levels)                                    &
, lambda_a (row_length) ! delta_lambda term for polar wind

! LAM LBC tendency

REAL                                                               &
  w_lbc_real_tend(lenrima(fld_type_p,halo_type_extended,           &
                    rima_type_norm),0:model_levels)                &
, exner_lbc_real_tend(lenrima(fld_type_p,halo_type_extended,       &
                    rima_type_norm),model_levels+1)

! Physics arrays needed by dynamics
REAL                                                               & 
          rho_km(rkmdims%i_start:rkmdims%i_end,                    &
                 rkmdims%j_start:rkmdims%j_end,                    &
                 rkmdims%k_start:rkmdims%k_end)                    &
,             ch(chdims %i_start:chdims %i_end,                    &
                 chdims %j_start:chdims %j_end,                    &
                 chdims %k_start:chdims %k_end)                    &
, wet_to_dry_np1(tdims_s%i_start:tdims_s%i_end,                    &
                 tdims_s%j_start:tdims_s%j_end,                    &
                 tdims_s%k_start:tdims_s%k_end)

! arrays holding information to be passed between physics
! routines.

REAL                                                               &
  ls_rain    (row_length, rows)                                    &
, ls_snow    (row_length, rows)                                    &
, micro_tends(row_length, rows, 2, bl_levels)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

! Radiation fields 1. SW & common with LW.
REAL                                                               &
  photosynth_act_rad(row_length, rows)                             &
                                       ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
, rad_hr(row_length, rows, 2, bl_levels)                           &
!                                 BL (LW,SW) rad heating rates
, dolr(row_length,rows)                                            &
!       local field "dolr" is distinguished from "dolr_field" 
!       (in atm_fields_mod)
                             ! TOA - surface upward LW
, sw_tile(land_field,ntiles)                                       &
                             ! Surface net SW on land tiles
, cos_zenith_angle(row_length, rows)


! MPP-related Arrays
INTEGER                                                            &
 g_row_length(0:nproc-1)                                           &
                         ! Table of number of points on a row
,g_rows(0:nproc-1)                                                 &
                         ! Table number of rows in theta field
,g_i_pe(1-halo_i:global_row_length+halo_i)                         &
               ! processor on my
!               processor-row holding a given value in i direction
,g_j_pe(1-halo_j:global_rows+halo_j) ! processor on my
!               processor-row holding a given value in j direction

      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: max_comm_size ! error check size for comms on demand

! Useful diagnostic
INTEGER istat
REAL    w_max, u_at_w, v_at_w


REAL stashwork0_dummy(1) ! STASHwork not defined for section 0,
                         !  but required as dummy argument.

! Local arrays for phys1 and phys2 increments for tracers:
REAL :: super_tracer_phys1(tdims_l%i_start:tdims_l%i_end,          &
                           tdims_l%j_start:tdims_l%j_end,          &
                           tdims_l%k_start:tdims_l%k_end,          &
                           super_array_size)

REAL :: super_tracer_phys2(tdims%i_start:tdims%i_end,              &
                           tdims%j_start:tdims%j_end,              &
                           tdims%k_start:tdims%k_end,              &
                           super_array_size)

REAL :: tracer_phys1(trdims_ltl%i_start:trdims_ltl%i_end,          &
                     trdims_ltl%j_start:trdims_ltl%j_end,          &
                     trdims_ltl%k_start:trdims_ltl%k_end, tr_vars)

REAL :: tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,        &
                     trdims_xstl%j_start:trdims_xstl%j_end,        &
                     trdims_xstl%k_start:trdims_xstl%k_end,        &
                     tr_vars)

REAL :: ukca_tracer_phys1(trdims_ltl%i_start:trdims_ltl%i_end,     &
                          trdims_ltl%j_start:trdims_ltl%j_end,     &
                          trdims_ltl%k_start:trdims_ltl%k_end,     &
                          tr_ukca)

REAL :: ukca_tracer_phys2(trdims_xstl%i_start:trdims_xstl%i_end,   &
                          trdims_xstl%j_start:trdims_xstl%j_end,   &
                          trdims_xstl%k_start:trdims_xstl%k_end,   &
                          tr_ukca)

REAL gs1(land_field)

! local variable
REAL                                                               &
  tot_dry_mass_final                                               &
                      ! mass at end of energy correction period
, tot_energy_final                                                 &
                      ! energy at end of energy correction period
, tot_moist_final                                                  &
                      ! moist at end of energy correction period
, energy_corr_now
                      ! instanteous energy correction

! Monsoon variables
REAL                                                               &
  lambda_half_width                                                &
, phi_half_width                                                   &
, p_max                                                            &
, p_top                                                            &
, p_bottom                                                         &
, p_arbitrary                                                      &
, lambda_heat_centre                                               &
, phi_heat_centre                                                  &
, max_heat_per_day                                                 &
, mons_newtonian_timescale

!  Suarez-Held variables now declared in CRUNTIMC

INTEGER    :: minutos ! LOCAL Store value of timestep in min.

!
! Structure for UKCA/radiation interaction.
!
      TYPE (ukca_radaer_struct), SAVE :: ukca_radaer

REAL                                                               &
  flux_e(row_length, rows)                                         &
                           ! Surface latent heat flux (W/m^2)
, flux_h(row_length, rows)                                         &
                           ! Surface sensible heat flux (W/m^2)
, z0m_scm(row_length, rows)                                        &
                           ! SCM specified z0m (m)
, z0h_scm(row_length, rows)! SCM specified z0m (m)
REAL, DIMENSION(:,:,:), ALLOCATABLE ::                             &
  q_inc_subs, th_inc_subs                                          &
                            ! subsidence increments
, q_inc_ls, th_inc_ls                                              &
                            ! large scale increments
, u_inc_dmp, q_inc_dmp, th_inc_dmp                                 &
                            !Damping incs
, v_inc_dmp
! Tolerance for CycleNo >1
REAL                                                               &
  gcr_run_tol_abs                                                  &
, gcr_run_tol_res

TYPE(swapable_field_pointer_type) :: fields_to_swap(7)  ! multivariate
                                                        ! swapbounds

INTEGER exppxi               ! Function to extract ppxref info

! Local parameters for mixing ratio physics
! Mixing ratios for atmos_physics1 and 2 are defined through namelist
LOGICAL, PARAMETER :: l_mr_qtbalcld=.FALSE. ! Use mr's for qt_bal_cld
LOGICAL, PARAMETER :: l_mr_iau     =.FALSE. ! Use mr's for tfilt_ctl
LOGICAL, PARAMETER :: l_mr_pc2     =.FALSE. ! Use mr's for PC2 routines

INTEGER idx

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL max_r_u,min_r_u ! local diagnostics
REAL max_r_v,min_r_v
REAL max_u, min_u
REAL max_v, min_v
REAL max_w, min_w
REAL max_thetav, min_thetav

REAL max_r_w, min_r_w
REAL max_r_theta, min_r_theta
REAL max_r_rho, min_r_rho
REAL max_r_p_p, min_r_p_p
REAL max_q_star, min_q_star

REAL solver_tolerance
INTEGER ierr

REAL mixing_r ! mixing ratio used in conversion from wet to dry density

LOGICAL, PARAMETER :: couple_app=.TRUE.   ! this parameter can be used
                   ! to force a dynamics only run even with the full physics
                   ! suite executing. For testing purposes. If false it will
                   ! zero all source terms coming from atmos_physics1 &
                   ! atmos_physics2 (as far as dynamics is concerned - mostly).
                   ! It could be set to l_physics, however that would exclude
                   ! the eg_idl_forcing contribution which is enabled by
                   ! disabling physics in the UMUI and setting the
                   ! appropriate IDEALISED namelist logical.

REAL    memory(100)

INTEGER n_filt, cycleno2
REAL tmp_filter1,tmp_filter2

! Input variables to be used by COSP (row_length,rows,model_levels)
!  Convective rainfall
REAL,POINTER,SAVE :: cosp_crain_3d(:,:,:)
!  Convective snowfall
REAL,POINTER,SAVE :: cosp_csnow_3d(:,:,:)

INTEGER :: lbc_size_new

TYPE (array_dims) twoddims
TYPE (array_dims) twoddims_no_halo

LOGICAL, SAVE :: subseq_to_dump = .FALSE.

!- End of header


twoddims=qdims_s
twoddims%k_end = twoddims%k_start
twoddims_no_halo=qdims
twoddims_no_halo%k_end = twoddims%k_start


alpha_changed = .FALSE.
reference_profile_changed = .FALSE.
L_int_uvw_lbc = .FALSE.

IF( GCR_USE_RESIDUAL_TOL ) THEN
solver_tolerance = gcr_tol_res*(tol_sc_fact**(numcycles*innits-1))
ELSE
solver_tolerance = gcr_tol_abs*(tol_sc_fact**(numcycles*innits-1))
END IF

IF (printstatus > prstatus_oper) THEN
  IF(.NOT.couple_app) WRITE(0,fmt='(A)') 'coupling disabled!'


! DEPENDS ON: eg_memory_usage
  IF( first_atmstep_call) THEN
    CALL eg_memory_usage(memory(1))
  ELSE
    CALL eg_memory_usage2(memory(1))
  END IF

  CALL gc_rsum(1,nproc,istat,memory(1))

  WRITE(6,fmt='(A,E32.16)') 'memory in use at atm_step_4A top:',memory(1)

END IF

! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('ATM_STEP_4A',zhook_in,zhook_handle)
! DEPENDS ON: timer
IF (Ltimer) CALL timer('Atm_Step_4A (AS)',5)


errorstatus = 0

!     ENDGame NEEDS the np1 variables
l_new_tdisc = .TRUE.

! DEPENDS ON: Atm_Step_Init
CALL atm_step_init (                                                  &
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
     lambda_a, g_row_length, g_rows, g_i_pe, g_j_pe, flux_e, flux_h,  &
     z0m_scm, z0h_scm, errorstatus )

l_2dcomm = .FALSE.  ! currently not employed by EG
If (l_2dcomm) Then
  size_2dcomm   = nproc - 1
  group_2dcomm  = gc_all_proc_group
  max_comm_size = model_levels * row_length * rows
Else
  size_2dcomm   = nproc_x - 1
  group_2dcomm  = gc_proc_row_group
  max_comm_size = model_levels * global_row_length
End If

check_bottom_levels = MIN(check_bottom_levels, model_levels)
interp_vertical_search_tol = MIN(interp_vertical_search_tol,          &
                                 model_levels/2)

! swapbound the exner from the dump!
IF (first_atmstep_call) THEN
   CALL eg_swap_bounds(exner,wdims_s,fld_type_p,.FALSE.)
   CALL eg_swap_bounds(etadot,tdims_s,fld_type_p,.FALSE.)
END IF
IF (integrity_test) THEN
  IF(me==0) WRITE(0,fmt='(A)')  'Start of timestep fields'
  CALL check_hash_m(                                                  &
                   exner,           SIZE(exner),           'pi___',   & 
                   rho,             SIZE(rho),             'rho__',   & 
                   thetav,          SIZE(thetav),          'tv___',   & 
                   u,               SIZE(u),               'u____',   &
                   v,               SIZE(v),               'v____',   &
                   w,               SIZE(w),               'w____',   &
                   m_v,             SIZE(m_v),             'm_v__',   &
                   m_cl,            SIZE(m_cl),            'm_cl_',   &
                   m_cf,            SIZE(m_cf),            'm_cf_',   &
                   m_r,             SIZE(m_r),             'm_r__',   &
                   m_gr,            SIZE(m_gr),            'm_grp',   &
                   m_cf2,           SIZE(m_cf2),           'm_cf2')
END IF

IF (first_atmstep_call) THEN

   ! allocate storage
   CALL init_ref_pro()
   CALL init_depart_pts()
   CALL init_fields_rhs(l_skeb2)
   CALL init_helmholtz_const_matrix()
   CALL init_gravity()

   CALL eg_init_helmholtz(row_length,rows,n_rows,model_levels,offx,offy)

!  rho elimination
   alpha_p = alpha_rho

   CALL eg_info_file(mype,nproc_x,nproc_y,numcycles,innits,            &
            global_row_length,global_rows,model_levels,                &
            ih,height_domain,l_shallow,l_rotating,l_rotate_grid,       &
            gcr_tol_abs,gcr_precon_option,                             &
            l_inc_solver, l_RK_dps,                                    &
            l_accel_convergence,l_eliminate_rho,                       &
            l_init_fnm1,                                               &
            alpha_u,alpha_v,alpha_w,alpha_theta,                       &
            alpha_rho,alpha_p,l_rotate_winds,l_baroclinic,l_solid_body,&
            l_baro_inst,l_baro_perturbed,                              &
            l_deep_baro_inst, T0_E, T0_P, b_const, k_const,            &
            l_const_grav,l_expl_horz_drag,                             &
            l_impl_horz_drag,l_eg_dry_static_adj,l_fix_mass,           &
            L_conserv_smooth_lap,                                      &
            eg_vert_damp_coeff,eg_vert_damp_profile,eta_s,grid_np_lon, &
            grid_np_lat,aa_jet_u0,l_cartesian,surface_type,            &
            grid_number,tprofile_number,trefer_number,t_surface,h_o,   &
            lambda_fraction,phi_fraction,half_width_x,half_width_y)
END IF

! We do not deallocate the np1 fields until the very end of the run. To avoid
! multiple allocation we have to if-test with timestep_number == 1
!
! DEPENDS ON: Atm_Step_alloc_4A
IF (first_atmstep_call)                                                &
  CALL atm_step_alloc_4A(                                              &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
    cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,           &
    frac_control, r_u, r_v, r_w, errorstatus, 'newtdisc')


! ---------------------------------------------------------------------
! Section 0.1  Initialisation for idealised test problems
!              For standard runs go to end section 0.1
! ---------------------------------------------------------------------
h_print=0.0


IF ( l_initialise_data .AND. first_atmstep_call ) THEN


  CALL eg_idl_ni_init(  pi,model_domain, row_length,                  &
  rows, n_rows, model_levels,wet_levels,tr_vars, tr_levels, bl_levels,&
  first_constant_r_rho_level                                          &
! NOTE: f3_at_v, f3_at_u was replaced by _star. Respective idealised
! need fixing.
  ,cos_theta_latitude, sec_theta_latitude, f3_at_u, f3_at_v,timestep, &
  first_atmstep_call, l_regular,delta_xi1,delta_xi2,delta_lambda,     &
  delta_phi,base_xi1,base_xi2,a_realhd(rh_rotlat),                    &
  a_realhd(rh_rotlong), phi_p, phi_v,r_theta_levels, r_rho_levels,    &
  r_at_u, r_at_v, r_at_u_w, r_at_v_w,z_orog_print, eta_theta_levels,  &
  eta_rho_levels                                                      &
! Multi-processor
  ,offx,offy,halo_i,halo_j, mype, nproc, at_extremity, datastart,     &
  gc_all_proc_group, global_row_length, global_rows,g_rows,           &
  g_row_length, nproc_x                                               &
! Primary fields
  ,thetav(:,:,1:),rho,exner_rho_levels,q, qcl, qcf,qcf2,qrain,qgraup,  &
  cf_bulk, cf_liquid,cf_frozen, u, v, w, u_adv, v_adv, w_adv          &
! Grid info for idealised
  ,a_realhd(rh_z_top_theta), height_domain, big_layers,transit_layers &
  , surface_type, p_surface                                           &
! Profile settings
  ,tprofile_number, qprofile_number, uvprofile_number,brunt_vaisala,  &
  t_surface, dtheta_dz1, height_dz1, u_in, v_in, height_u_in,         &
  ujet_lat, ujet_width, u_ramp_start, u_ramp_end, f_plane, r_plane,q1,&
  num_profile_data, zprofile_data, tprofile_data,qprofile_data,       &
  num_uvprofile_data, z_uvprofile_data,uprofile_data,vprofile_data,   &
  model_levels_max, max_num_profile_data, max_num_force_times,        &
  tforce_option, qforce_option, uvforce_option, num_tforce_levels,    &
  num_tforce_times, num_qforce_levels, num_qforce_times,              &
  num_uvforce_levels, num_uvforce_times, z_tforce_data, tforce_data,  &
  z_qforce_data, qforce_data, z_uvforce_data, uforce_data,vforce_data,&
  tforce_data_modlev, qforce_data_modlev,uforce_data_modlev,          &
  vforce_data_modlev                                                  &
! Dynamical core settings
  ,suhe_pole_equ_deltat, suhe_static_stab,base_frictional_timescale,  &
  frictional_timescale,suhe_sigma_cutoff, suhe_level_weight,          &
  l_sh_williamson                                                     &
! Horizontal function parameters
  ,t_horizfn_number, uv_horizfn_number,t_horizfn_data, l_perturb_t,   &
  perturb_magnitude_t,l_perturb_q, perturb_magnitude_q,               &
  l_perturb_correlate_tq,l_perturb_correlate_vert,                    &
  l_perturb_correlate_time,perturb_type, perturb_height               &
! Profiles for fixed lbcs and sponge zones
  ,u_ref, v_ref, theta_ref, exner_ref, rho_ref,zprofile_orog,         &
  idl_interp_option, hf                                               &
 ! Options
  ,l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,l_constant_dz, l_rotating,   &
  l_fixed_lbcs, l_polar_wind_zero,l_wind_balance, l_rotate_winds,     &
  l_pressure_balance, l_physics,l_dry, l_sponge, l_cartesian,         &
  l_code_test, l_cyclone, l_baroclinic ,h_print, timestep_number,     &
  h_o_actual, grow_steps, h_o_per_step,h_o, grid_number, grid_flat,   &
  first_theta_height,thin_theta_height, big_factor, vert_grid_ratio,  &
  lambda_fraction,phi_fraction,half_width_x, half_width_y,            &
  idl_max_num_bubbles, idl_bubble_option, idl_bubble_max,             &
  idl_bubble_height, idl_bubble_xoffset, idl_bubble_yoffset,          &
  idl_bubble_width, idl_bubble_depth, l_idl_bubble_saturate,a_ixsts,  &
  len_a_ixsts, a_spsts,len_a_spsts,nproc_y, gc_proc_row_group ,       &
  idlsurffluxseaoption, idlsurffluxseaparams, l_flux_bc,              &
  flux_h,flux_e, i_hour, i_minute, i_second, problem_number,orography )

ELSE IF( problem_number  /=  standard .AND.                           &
         first_atmstep_call ) THEN
    IF (.NOT. l_physics .AND. mype  ==  0)THEN
     PRINT*,'Data from dump being used without initialising. '
     PRINT*,'If source of dump is a full model run with orography'
     PRINT*,'then there is no guarantee that run will work'
     PRINT*,'since physics is OFF'
    END IF !.not. L_Physics.and. mype  ==  0

END IF   ! L_initialise_data

!  ********************************************************************
!  Rescale rho. Store rho * r**2 in rho_r_sq for use in
!  03to3D, atmos_physics1/2 and stph_skeb2
!  ********************************************************************

DO k=1, model_levels
  rho_r_sq_n(:,:,k) = rho(:,:,k)
END DO

IF (first_atmstep_call) THEN

  CALL eg_swap_bounds(rho_r_sq_n,pdims_s,fld_type_p,.FALSE.)

  CALL eg_sisl_setcon(                                                  &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,        &
         offx, offy,  datastart, l_slice, l_test_tracer,                &
         l_cartesian, l_shallow, l_const_grav,l_impl_horz_drag,         &
         l_expl_horz_drag,eccentricity, height_domain,                  &
         alpha_u, alpha_v, alpha_w, alpha_theta, alpha_rho,             &
         tprofile_number, trefer_number, dtheta_dz1,                    &
         a_realhd(rh_z_top_theta), g_theta,  g_rho,                     &
         t_surf_ref, p_surf_ref,eg_vert_damp_coeff,eg_vert_damp_profile,&
         eta_s, thetav, rho, t_surface, p_surface, exner,               &
         alpha_p, l_eliminate_rho,                                      &
         ih, pole_consts)

END IF

IF ( timestep_number > 1 ) THEN
! density is recomputed below following the change from thetav to thetavd,
! therefore we do not need to rescale/convert it here at the first timestep

  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end

        mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)      +                &
                                     m_r  (i,j,k)      +                &
                                     m_gr (i,j,k)      +                &
                                     m_cl (i,j,k)      +                &
                                     m_cf (i,j,k)      +                &
                                     m_cf2(i,j,k))     +                &
                    intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                &
                                     m_r  (i,j,k-1)    +                &
                                     m_gr (i,j,k-1)    +                &
                                     m_cl (i,j,k-1)    +                &
                                     m_cf (i,j,k-1)    +                &
                                     m_cf2(i,j,k-1)))

!       convert to dry density
        rho(i,j,k) = rho(i,j,k)/(1.+mixing_r)

!       unscale
        rho(i,j,k) = rho(i,j,k)/(r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
                    
      END DO
    END DO
  END DO

  CALL eg_swap_bounds(rho,pdims_s,fld_type_p,.FALSE.)

END IF


! initialise moisture arrays, metric terms helmholtz coeffs etc

firstAScall: IF ( first_atmstep_call ) THEN

  IF (timestep_number == 1) THEN
  ! this is the conversion from ND. Should not be required anymore
  ! if restarting
  CALL eg_q_to_mix                                                      &
                  (tdims_l,tdims_s,                                     &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup,                                 &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   m_v, m_cl, m_cf                                      &
                  ,m_cf2, m_r, m_gr)

  IF (.NOT.l_mcr_qcf2)   m_cf2(:,:,:) = 0.
  IF (.NOT.l_mcr_qrain)  m_r  (:,:,:) = 0.
  IF (.NOT.l_mcr_qgraup) m_gr (:,:,:) = 0.
  END IF


  tpeqtpdump: IF(tprofile_number == tp_dump) THEN

    ts1 : IF (timestep_number == 1) THEN

! Initialise surface pressure from dump (assumes conversion from p->exner)
      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          exner_surf(i,j) = (pstar(i,j)/p_zero)**kappa
      END DO
    END DO
! DEPENDS ON: swap_bounds
      CALL swap_bounds(exner_surf,row_length,rows,1,offx,offy,                &
                       fld_type_p,.FALSE.)


    IF(pdims%j_end-n_filt_ND+1.lt.pdims%j_start.OR.                          &
       pdims%j_start+n_filt_ND-1.gt.pdims%j_end) THEN

        ErrorStatus = ABS(pdims%j_end-n_filt_ND+1)
        WRITE(Cmessage,fmt='(A,I4.4,A,I4.4,A,I4.4)') &
          "Cannot filter. Reduce NS decomposition.", &
          pdims%j_end,' ',pdims%j_start,' ',n_filt_ND
        CALL Ereport("ATM_STEP_4A", ErrorStatus, Cmessage)

    END IF

     IF( (model_domain == mt_global) .AND. l_filter_dump ) THEN
        IF( mype  ==  0 )                                 &
          WRITE(0,fmt='(A,I4)') 'Filtering initial dump data. n_filt=',n_filt_ND

        DO n_filt = 1, n_filt_ND

          IF(at_extremity(psouth)) THEN
            DO k = pdims%k_start,pdims%k_end + 1
              DO j = pdims%j_start,pdims%j_start+n_filt-1

                tmp_filter1 = exner(pdims%i_start-1,j,k)
 
                DO i = pdims%i_start,pdims%i_end
                   tmp_filter2  = exner(i,j,k)
                exner(i,j,k) = 0.25*(tmp_filter1+2.*exner(i,j,k)+exner(i+1,j,k))
                   tmp_filter1  = tmp_filter2
                END DO
              END DO
            END DO
          END IF

          IF(at_extremity(pnorth)) THEN
            DO k = pdims%k_start,pdims%k_end + 1
              DO j = pdims%j_end-n_filt+1,pdims%j_end

                tmp_filter1 = exner(pdims%i_start-1,j,k)
 
                DO i = pdims%i_start,pdims%i_end
                  tmp_filter2  = exner(i,j,k)
                  exner(i,j,k) = 0.25*(tmp_filter1+2.*exner(i,j,k)           &
                                      +exner(i+1,j,k))
                  tmp_filter1  = tmp_filter2
                END DO
              END DO
            END DO
          END IF

          CALL eg_swap_bounds(exner,wdims_s,fld_type_p,.FALSE.)

          IF(at_extremity(psouth)) THEN
            DO j = pdims%j_start,pdims%j_start+n_filt-1
            tmp_filter1 = exner_surf(pdims%i_start-1,j)
              DO i = pdims%i_start,pdims%i_end
              tmp_filter2   = exner_surf(i,j)
              exner_surf(i,j) = 0.25*(tmp_filter1+2.*exner_surf(i,j)+         &
                                                     exner_surf(i+1,j))
                tmp_filter1   = tmp_filter2
              END DO
            END DO
          END IF

          IF(at_extremity(pnorth)) THEN
            DO j = pdims%j_end-n_filt+1,pdims%j_end
            tmp_filter1 = exner_surf(pdims%i_start-1,j)
              DO i = pdims%i_start,pdims%i_end
              tmp_filter2   = exner_surf(i,j)
              exner_surf(i,j) = 0.25*(tmp_filter1+2.*exner_surf(i,j)          &
                                                    +exner_surf(i+1,j))
                tmp_filter1   = tmp_filter2
              END DO
            END DO
          END IF

! DEPENDS ON: swap_bounds
        CALL swap_bounds(exner_surf,row_length,rows,1,offx,offy,              &
                       fld_type_p,.FALSE.)

        END DO
     END IF


! Here we compute
! the virtual dry thetav from theta and the reinitialise the
! density using the equation of state.

  !====================================================================
  ! this is the conversion from ND. Should not be required anymore
  ! if restarting

    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetav(i,j,k) = theta(i,j,k)*(1.0+m_v(i,j,k)*recip_epsilon)
        END DO
      END DO
    END DO

! Recompute the dry density to satisfy equation of state
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho(i,j,k) = p_zero*exner(i,j,k)**((1.0-kappa)/kappa) /       &
                       (R*(intw_w2rho(k,1)*thetav(i,j,k) +              &
                           intw_w2rho(k,2)*thetav(i,j,k-1)))
        END DO
      END DO
    END DO

! Remove any potential instabillity:
    IF(l_eg_dry_static_adj) THEN
      CALL eg_dry_static_adj(thetav,rho,exner,intw_w2rho)
    END IF

    CALL eg_swap_bounds(rho,pdims_s,fld_type_p,.FALSE.)
    CALL eg_swap_bounds(thetav,tdims_s,fld_type_p,.FALSE.)


! rho has been changed by the above calculations, hence need to recompute
! rho_r_sq
DO k=1, model_levels
  DO j=pdims_s%j_start, pdims_s%j_end
    DO i=pdims_s%i_start, pdims_s%i_end

      mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)      +                  &
                                   m_r  (i,j,k)      +                  &
                                   m_gr (i,j,k)      +                  &
                                   m_cl (i,j,k)      +                  &
                                   m_cf (i,j,k)      +                  &
                                   m_cf2(i,j,k))     +                  &
                  intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                  &
                                   m_r  (i,j,k-1)    +                  &
                                   m_gr (i,j,k-1)    +                  &
                                   m_cl (i,j,k-1)    +                  &
                                   m_cf (i,j,k-1)    +                  &
                                   m_cf2(i,j,k-1)))

!     convert to wet rho
      rho_r_sq_n(i,j,k) = rho(i,j,k)*(1.+mixing_r)    
!     multiply by r**2
      rho_r_sq_n(i,j,k) = rho_r_sq_n(i,j,k)*( r_rho_levels(i,j,k)       &
                                           *r_rho_levels(i,j,k))
    END DO
  END DO
END DO
  !====================================================================


! Fix v at the poles
    IF( model_domain == mt_global ) THEN
      IF( at_extremity(psouth) ) THEN

        CALL eg_v_at_poles(u,v, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)
      END IF

      IF( at_extremity(pnorth) ) THEN

        CALL eg_v_at_poles(u,v,-1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)
      END IF
    END IF

   END IF ts1

    IF (integrity_test)                                               &
      CALL update_hash_m(                                             &
                  exner,           SIZE(exner),           'pi___',    & 
                  u,               SIZE(u),               'u____',    &
                  v,               SIZE(v),               'v____',    &
                  w,               SIZE(w),               'w____',    & 
                  u_adv,           SIZE(u_adv),           'u_adv',    &
                  v_adv,           SIZE(v_adv),           'v_adv',    &
                  w_adv,           SIZE(w_adv),           'w_adv')


  END IF tpeqtpdump

  chain_number=0

  CALL eg_idl_set_init(                                                 &
                      row_length, rows, n_rows, halo_i, halo_j,         &
                      offx, offy, model_levels,qprofile_number,         &
                      tprofile_number,                                  &
                      t_surface, p_surface,                             &
                      l_baro_inst, l_HeldSuarez, l_solid_body,          &
                      l_deep_baro_inst, T0_E, T0_P, b_const, k_const,   &
                      l_baro_perturbed, l_shallow, l_const_grav,        &
                      l_isothermal, l_rotate_grid, L_initialise_data,   &
                      L_cartesian,                                      &
                      idl_bubble_option, idl_max_num_bubbles,           &
                      idl_bubble_max, idl_bubble_width,                 &
                      idl_bubble_height,                                &
                      idl_bubble_depth, idl_bubble_xoffset,             &
                      idl_bubble_yoffset,                               &
                      datastart, global_row_length, global_rows,        &
                      grid_np_lon, grid_np_lat,                         &
                      aa_jet_u0, aa_jet_a, aa_jet_m, aa_jet_n,          &
                      ring_height, theta_pert,                          &
                      u, v, w, u_adv, v_adv, w_adv, rho, exner,         &
                      exner_surf,u_np1, v_np1, w_np1, etadot_np1,       &
                      thetav_np1, rho_np1, exner_np1,                   &
                      m_v_np1, m_cl_np1, m_cf_np1,                      &
                      m_r_np1, m_gr_np1, m_cf2_np1, exner_surf_np1,     &
                      a_realhd(rh_z_top_theta), dtheta_dz1,             &
                      l_RK_dps, l_dry, alpha_w, ih,                     &
                      etadot, psi_w_surf, psi_w_lid,                    &
                      thetav,m_v, m_cl, m_cf,m_r, m_gr,                 &
                      m_cf2)

  idealdata : IF (l_idealised_data) THEN
  ! rho has been changed by the above calculations, hence need to recompute
  ! rho_r_sq
    DO k=1, model_levels
      DO j=pdims_s%j_start, pdims_s%j_end
        DO i=pdims_s%i_start, pdims_s%i_end

          mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)      +                  &
                                       m_r  (i,j,k)      +                  &
                                       m_gr (i,j,k)      +                  &
                                       m_cl (i,j,k)      +                  &
                                       m_cf (i,j,k)      +                  &
                                       m_cf2(i,j,k))     +                  &
                      intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                  &
                                       m_r  (i,j,k-1)    +                  &
                                       m_gr (i,j,k-1)    +                  &
                                       m_cl (i,j,k-1)    +                  &
                                       m_cf (i,j,k-1)    +                  &
                                       m_cf2(i,j,k-1)))

          !     convert to wet rho
          rho_r_sq_n(i,j,k) = rho(i,j,k)*(1.+mixing_r)    
          !     multiply by r**2
          rho_r_sq_n(i,j,k) = rho_r_sq_n(i,j,k)*( r_rho_levels(i,j,k)       &
                                           *r_rho_levels(i,j,k))
        END DO
      END DO
    END DO 
  END IF idealdata
END IF firstAScall

IF (integrity_test_ghash) THEN
  IF ( timestep_number == 1 ) THEN
    IF(me==0) OPEN (eg_unit,file='start_hash.dat')
  ELSE
   IF(me==0) OPEN (eg_unit,file='start_hash.dat', POSITION='APPEND')
  END IF
  IF(me==0) THEN
    WRITE(eg_unit,fmt='(A)')   '############################'
    WRITE(eg_unit,fmt='(A,I5)') 'Global checksum for timestep ',timestep_number
  END IF
  CALL  global_hash(u,      'u         ',udims_s,fld_type_u,eg_unit)
  CALL  global_hash(v,      'v         ',vdims_s,fld_type_v,eg_unit)
  CALL  global_hash(w,      'w         ',wdims_s,fld_type_p,eg_unit)
  CALL  global_hash(etadot, 'etadot    ',wdims_s,fld_type_p,eg_unit)
  CALL  global_hash(thetav, 'thetav    ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(rho,    'rho       ',pdims_s,fld_type_p,eg_unit)
  CALL  global_hash(exner,  'exner     ',pdims_s,fld_type_p,eg_unit)
  CALL  global_hash(rho_r_sq_n,  'rho_r_sq_n',pdims_s,fld_type_p,eg_unit)
  CALL  global_hash(exner_surf,'pstar     ',twoddims,fld_type_p,eg_unit)
  CALL  global_hash(psi_w_surf,'psiwsurf  ',twoddims_no_halo,fld_type_p,eg_unit)
  CALL  global_hash(psi_w_lid,'psiwlid   ',twoddims_no_halo,fld_type_p,eg_unit)
  IF (.NOT. l_dry) THEN
    CALL  global_hash(m_v,  'm_v       ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cl, 'm_cl      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cf, 'm_cf      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_r,  'm_r       ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_gr, 'm_gr      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cf2,'m_cf2     ',qdims_s,fld_type_p,eg_unit)
  END IF
  CALL  global_hash(SO4_AITKEN,'so4 aitken',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SO4_ACCU,'acc       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SO4_DISS,'diss      ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(NH3,'nh3       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SOOT_NEW     ,'soot new  ',tdims_s,fld_type_p,eg_unit) 
  CALL  global_hash(SOOT_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SOOT_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_NEW,'bmass     ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_NEW,'ocff new  ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  IF(me==0) CLOSE (eg_unit)
END IF

!======== LBC Update =========================================================
! ---------------------------------------------------------------------
!   Section 0.2  Update lbcs for LAMs
! ---------------------------------------------------------------------
IF ((ErrorStatus == 0) .and. (model_domain == mt_lam)) THEN
! DEPENDS ON: timer
  IF (Ltimer) CALL timer('AS LAM_LBCS',5)
  
  IF ( first_atmstep_call ) THEN 
    lbc_size=LENRIMA(fld_type_u,halo_type_extended, rima_type_norm)
    lbc_size_new = LENRIMA(fld_type_p,halo_type_extended, rima_type_norm)
    
! Reset u_lbc array to same size at p lbc array    
    LENRIMA(fld_type_u,halo_type_extended, rima_type_norm) =                &
                      LENRIMA(fld_type_p,halo_type_extended, rima_type_norm)
  END IF

!--------------------------------------------------------------
! Idealised UM LBC forcing
!  If active, update lateral boundary arrays to contain
!  idealised namelist profile data interpolated in time.
!--------------------------------------------------------------
  IF (L_initialise_data .and. L_force_lbc) THEN

! DEPENDS ON: eg_idl_force_lbc
    CALL eg_IDL_Force_LBC (                                           &
                      row_length, rows, offx, offy,                   &
                      halo_i, halo_j,                                 &
                      LENRIMA(1,1,rima_type_norm),                    &
                      timestep, timestep_number,                      &
                      model_levels, wet_levels,                       &
                      model_levels_max, max_num_force_times,          &
                      u_lbc, V_LBC,                                   &
                      THETA_LBC,Q_LBC,                                &
                      U_ADV_LBC,V_ADV_LBC,                            &
                      EXNER_LBC,                                      &
                      r_theta_levels, r_rho_levels,                   &
                      eta_theta_levels, eta_rho_levels,               &
                      height_domain, theta_surface,                   &
                      pforce_option,                                  &
                      tforce_option, qforce_option, uvforce_option,   &
                      num_pforce_times,                               &
                      num_tforce_times, num_qforce_times,             &
                      num_uvforce_times,                              &
                      pforce_time_interval,                           &
                      tforce_time_interval, qforce_time_interval,     &
                      uvforce_time_interval,                          &
                      p_surface_data,                                 &
                      tforce_data_modlev, qforce_data_modlev,         &
                      uforce_data_modlev, vforce_data_modlev,         &
                      newtonian_timescale )
  END IF ! on (L_initialise_data .and. L_force_lbc)

  IF ( first_atmstep_call ) THEN
    L_apply_lbcs = .TRUE.
  ELSE
    L_apply_lbcs = .FALSE.
  END IF
  
  IF (model_domain /= mt_global) THEN
  IF ( L_apply_lbcs ) THEN
!--------------------------------------------------------------
!           Update primary fields with LAM LBC data
!--------------------------------------------------------------
!DEPENDS ON: update_lam_lbcs
    CALL update_lam_lbcs(                                              &
         r_rho_levels, r_theta_levels,                                 &
         ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
         TR_VARS,TR_LBC_VARS,TR_LEVELS,                                &
         A_MAX_TRVARS,A_tr_active_lbc_index,                           &
         TR_UKCA,TR_LBC_UKCA,                                          &
         A_MAX_UKCAVARS,UKCA_tr_active_lbc_index,                      &
         OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
         L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
         L_murk, L_murk_lbc,                                           &
         L_LBC_balance, L_int_uvw_lbc,                                 &
         L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,   &
         L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,   &
         L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,   &
         L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,&
         L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,       &
         L_nh3, L_nh3_lbc,                                             &
         L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,       &
         L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,     &
         L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,   &
         L_ocff_new, L_ocff_new_lbc,                                   &
         L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,       &
         L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,     &
         RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
         LENRIMA(1,1,rima_type_norm),                                  &
         LBC_SIZEA(1,1,1,rima_type_norm),                              &
         LBC_STARTA(1,1,1,rima_type_norm),                             &
         THETA_LBC,Q_LBC,QCL_LBC,                                      &
         QCF_LBC,QCF2_LBC,QRAIN_LBC,                                   &
         QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                        &
         CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
         u_lbc,V_LBC,W_LBC,                                            &
         U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                                &
         MURK_LBC,                                                     &
         DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                  &
         DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                  &
         SO2_LBC, DMS_LBC, SO4_AITKEN_LBC,                             &
         SO4_ACCU_LBC,SO4_DISS_LBC,NH3_LBC,                            &
         SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,                     &
         BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                  &
         OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                     &
         NITR_ACC_LBC, NITR_DISS_LBC,                                  &
         TRACER_LBC,TRACER_UKCA_LBC,                                   &
         qdims_s, thetav,m_v,m_cl,m_cf,                                &
         m_cf2, m_r, m_gr,                                             &
         CF_BULK,CF_LIQUID,CF_FROZEN,                                  &
         RHO,EXNER,                                                    &
         U,V,W,                                                        &
         U_ADV,V_ADV,W_ADV,                                            &
         MURK,                                                         &
         DUST_DIV1, DUST_DIV2, DUST_DIV3,                              &
         DUST_DIV4, DUST_DIV5, DUST_DIV6,                              &
         SO2, DMS, SO4_AITKEN,SO4_ACCU,                                &
         SO4_DISS, NH3,                                                &
         SOOT_NEW, SOOT_AGD, SOOT_CLD,                                 &
         BMASS_NEW, BMASS_AGD, BMASS_CLD,                              &
         OCFF_NEW, OCFF_AGD, OCFF_CLD,                                 &
         NITR_ACC, NITR_DISS,                                          &
         DELTA_PHI, DELTA_LAMBDA,                                      &
         BASE_PHI, BASE_LAMBDA,                                        &
         DATASTART, lat_rot_NP,                                        &
         GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                               &
         TRACER, TRACER_UKCA )
! Must now re-calculate the pressure-based variables, namely pressure
! on both rho and theta levels, exner on theta levels and pstar so that
! they are all consistent with the new LBC-updated values of exner on
! rho levels.

! *** NOTE ***
! This is picked up by the later eg_thetav_theta call, but might need 
! a call added here as well
  END IF ! L_apply_lbcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMP. FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------!
! Re-initialise etadot: fixes needed for lateral bcs                    !
!-----------------------------------------------------------------------!
  DO k = 1, model_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

          u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +   &
                                     intw_u2p(i,2)*u(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +     &
                                     intw_u2p(i,2)*u(i,j,k) )

          v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +   &
                                     intw_v2p(j,2)*v(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +     &
                                     intw_v2p(j,2)*v(i,j,k) )

          etadot(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -                  &
                             u_at_w*dxi1_xi3(i,j,k)/                    &
                                           h1_p_eta(i,j,k) -            &
                             v_at_w*dxi2_xi3(i,j,k)/                    &
                                           h2_p_eta(i,j,k) ) /          &
                                             deta_xi3_theta(i,j,k)
      END DO
    END DO
  END DO

  etadot(:,:,0) = 0.0
  etadot(:,:,model_levels) = 0.0

  CALL eg_swap_bounds(etadot,wdims_s,fld_type_p,.FALSE.)

!!!!!!!!!!!!!!!!!!!!!!!!!! END TEMP. FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Lateral boundaries have been updated so update _adv fields_rhs_mod
  CALL eg_set_adv_winds(u,v,etadot,                                    &
                        u_adv,v_adv,w_adv,row_length,rows,n_rows,      &
                        model_levels, halo_i, halo_j, l_shallow)
                        
! time level n fields have been updated so update guess for np1 fields                                                
  END IF   ! .NOT. GLOBAL
  
! DEPENDS ON: timer
  IF (Ltimer) CALL timer('AS LAM_LBCS',6)
END IF     !   model_domain  ==  mt_lam
! ---------------------------------------------------------------------
! End Section 0.2 Update lbcs for LAMs
! ---------------------------------------------------------------------      
!======== LBC Update =========================================================

    CALL eg_swap_bounds(rho,pdims_s,fld_type_p,.FALSE.)

    CALL eg_swap_bounds(  u,udims_s,fld_type_u,.TRUE.) 
    CALL eg_swap_bounds(  v,vdims_s,fld_type_v,.TRUE.) 
    CALL eg_swap_bounds(  w,wdims_s,fld_type_p,.FALSE.) 
    CALL eg_swap_bounds(  etadot,wdims_s,fld_type_p,.FALSE.) 

    exner_np1(:,:,:) = exner(:,:,:)
    u_np1(:,:,:)     = u(:,:,:)
    v_np1(:,:,:)     = v(:,:,:)
    w_np1(:,:,:)     = w(:,:,:)
    etadot_np1(:,:,:)     = etadot(:,:,:)
    exner_surf_np1   = exner_surf  

    rho_np1     = rho
    thetav_np1  = thetav
    m_v_np1     = m_v
    m_cl_np1    = m_cl
    m_cf_np1    = m_cf
    m_r_np1     = m_r
    m_gr_np1    = m_gr
    m_cf2_np1   = m_cf2


IF (integrity_test) THEN
  IF( mype  ==  0 ) WRITE (0,fmt='(A,I5)') '******* Field Diagnostics @ '&
                   ,timestep_number
  CALL update_hash_m(exner,           SIZE(exner),           'pi___',  & 
                   rho,             SIZE(rho),             'rho__',    & 
                   thetav,          SIZE(thetav),          'tv___',    & 
                   u,               SIZE(u),               'u____',    &
                   v,               SIZE(v),               'v____',    &
                   w,               SIZE(w),               'w____',    &
                   exner_np1,       SIZE(exner_np1),       'pinp1',    &
                   u_np1,           SIZE(u_np1),           'u_np1',    &
                   v_np1,           SIZE(v_np1),           'v_np1',    & 
                   u_adv,           SIZE(u_adv),           'u_adv',    &
                   v_adv,           SIZE(v_adv),           'v_adv',    &
                   w_adv,           SIZE(w_adv),           'w_adv',    &
                   psi_w_surf,      SIZE(psi_w_surf),      'psiws',    &
                   rho_np1,         SIZE(rho_np1),         'r_np1',    &
                   m_v_np1,         SIZE(m_v_np1),         'mvnp1',    &
                   m_cl_np1,        SIZE(m_cl_np1),        'mclp1',    &
                   m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',    &
                   m_r_np1,         SIZE(m_r_np1),         'mrnp1',    &
                   m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',    &
                   m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',    &
                   thetav_np1,      SIZE(thetav_np1),      'tvnp1',    &
                   w_np1,           SIZE(w_np1),           'w_np1',    &
                   etadot,          SIZE(etadot),          'ed___',    &
                   etadot_np1,      SIZE(etadot_np1),      'ednp1')
  IF( mype  ==  0 ) WRITE(0,fmt='(A)')                                 &
                      '********************************************'
END IF

CALL eg_alpha_ramp(chain_number)

! DEPENDS ON: timer
  IF (Ltimer) CALL timer ('AS Solver',5)


  CALL eg_sisl_resetcon(                                                &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,        &
         offx, offy,  datastart,l_slice, l_test_tracer,                 &
         l_cartesian, l_shallow, l_const_grav, eccentricity,            &
         height_domain, tprofile_number, trefer_number,                 &
         dtheta_dz1, a_realhd(rh_z_top_theta),t_surf_ref, p_surf_ref,   &
         thetav, rho, t_surface, p_surface, exner_surf,exner,           &
         l_inc_solver,                                                  &
         m_v, m_r, m_gr, m_cl, m_cf, m_cf2,                             &
         f1_comp, f2_comp, f3_comp, ih )

! DEPENDS ON: timer

  CALL eg_set_helm_lhs(row_length, rows, n_rows, model_levels, ih)

  IF (Ltimer) CALL timer ('AS Solver',6)


IF ( first_atmstep_call .OR. subseq_to_dump ) total_rho_init=eg_total_mass(rho)
IF ( subseq_to_dump ) THEN
  CALL eg_swap_bounds(u_adv,udims_l,fld_type_u,.TRUE.)
  CALL eg_swap_bounds(v_adv,vdims_l,fld_type_v,.TRUE.)
  CALL eg_swap_bounds(w_adv,wdims_l,fld_type_p,.FALSE.)
  subseq_to_dump = .FALSE.
END IF

IF ( first_atmstep_call .AND.  timestep_number == 1 ) THEN

     w_max = MAXVAL(ABS(w))
     CALL gc_rmax(1,nproc,istat,w_max)

     IF( mype == 0 .AND. printstatus >= prstatus_diag) THEN
       OPEN(UNIT=eg_unit,FILE='wmax.dat')
       IF (L_fix_mass) THEN
         WRITE(UNIT=eg_unit,fmt='(I7,2E25.16)') 0, w_max, 1.0
       ELSE
         WRITE(UNIT=eg_unit,fmt='(I7,2E25.16)') 0, w_max, total_rho
       END IF
       CLOSE(UNIT=eg_unit)
     END IF

END IF

IF ( first_atmstep_call .AND. chain_number == 0 ) THEN
  dx_frame_cnt = 0
END IF

! Save initial data for plotting if requested
IF( tstep_plot_start >= 0 .AND.                                            &
    timestep_number >= tstep_plot_start) THEN

   IF ( MOD(timestep_number,tstep_plot_frequency) == 0                     &
        .OR. dx_frame_cnt == 0 ) THEN

     CALL eg_dx_diags(dx_frame_cnt, mype,                                  &
                      row_length, rows, n_rows, model_levels,              &
                      offx, offy, halo_i, halo_j,                          &
                      u, v, w, rho, thetav, exner, exner_surf,             &
                      l_dry, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)
   END IF
END IF

!-----------------------------------------------------------------------
! Recompute "starred" Coriolis terms
!-----------------------------------------------------------------------

! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Solver',5)

CALL eg_coriolis_star(rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Solver',6)


! ----------------------------------------------------------------------
! Section 1.0  Initialise semi-lagrangian
! ----------------------------------------------------------------------

! Set timelevel n dependent quantities

r_u(:,:,:)       = 0.0
r_v(:,:,:)       = 0.0
r_w(:,:,:)       = 0.0
r_theta(:,:,:)   = 0.0
r_m_v(:,:,:)     = 0.0
r_m_cl(:,:,:)    = 0.0
r_m_cf(:,:,:)    = 0.0
r_m_r(:,:,:)     = 0.0
r_m_gr(:,:,:)    = 0.0
r_m_cf2(:,:,:)   = 0.0
r_p_p(:,:,:)     = 0.0
r_m_gr_d(:,:,:)  = 0.0
r_m_cf2_d(:,:,:) = 0.0
r_p_p_d(:,:,:)   = 0.0

l_call_from_solver = .FALSE.

IF (l_physics .AND. (.NOT.l_mr_physics1)) THEN

  CALL eg_mix_to_q                                                      &
                  (tdims_l,tdims_s,                                     &
                   m_v, m_cl, m_cf,                                     &
                   m_cf2, m_r, m_gr,                                    &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup)

!-----------------------------------------------------------------------
! Allocate and initialise arrays for total timestep increments
! NOTE it is very important this is always done before any process
!      changes the fields at the beginning of the timestep.
! This also needs to work on q (specific humidities, not mixing ratios)
! thus requiring ugly use of IF test.
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)

! Outside cycle_no loop so no test on cycle_no
! DEPENDS ON: Atm_Step_stash
  CALL atm_step_stash(                                                  &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     errorstatus, 5)

ELSE

! Ensure using q before overwritten by MRs.
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)

! DEPENDS ON: Atm_Step_stash
  CALL atm_step_stash(                                                  &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     errorstatus, 5)


  q  (1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_v (:,:,:)
  qcl(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_cl(:,:,:)
  qcf(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_cf(:,:,:)

  IF (l_mcr_qcf2  )                                                     &
     qcf2  (1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_cf2(:,:,0:model_levels)

  IF (l_mcr_qrain )                                                     &
     qrain (1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_r(:,:,0:model_levels)

  IF (l_mcr_qgraup)                                                     &
     qgraup(1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_gr(:,:,0:model_levels)

END IF

! Compute theta from thetav _AND_ exner_surf _AND_ p_theta_levels
! for use in physics routines
CALL eg_thetav_theta                                                    &
                  (thetav, theta, m_v,                                  &
                   p, pstar, p_theta_levels,                            &
                   exner, exner_surf, exner_theta_levels)


! ---------------------------------------------------------------------

CALL wet_to_dry_n_calc(q,qcl,qcf,qcf2,qrain,qgraup)

! ---------------------------------------------------------------------

h_print=0.0

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',6)

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS Stochastic_Phys',5)
!  --------- UM Section 35---- Stochastic Physics Setup -------------
!  Called only at the first time-step with no arguments
IF (first_atmstep_call) THEN
  IF (l_rp2 .OR. l_skeb2) THEN
    CALL stph_setup ( )
  END IF
END IF
!  --------- UM Section 35---- Stochastic Physics Setup END ---------


!  ----- UM Section 35---- Stochastic Physics Random Parameters -----
! Call to the RANDOM PARAMETERS2 (STPH_RP2) subroutine
IF (l_physics .AND. l_rp2) THEN
  minutos=timestep/60
  IF(MOD(i_hour,3) == 0 .AND. i_minute == minutos) THEN
    IF (printstatus  >=   prstatus_normal) THEN
      WRITE(6,fmt='(A)') 'CALLING RANDOM PARAMETERS2'
    END IF
          CALL stph_rp2(model_levels,                                   &
                  rhcrit,rhcrit_max,rhcrit_min,                         &
                  gwd_frc,gwd_frc_max,gwd_frc_min,                      &
                  kay_gwave,kay_gwave_max,kay_gwave_min,                &
                  m_ci,m_ci_max,m_ci_min,                               &
                  charnock)
  ELSE
   IF (printstatus  >=  prstatus_normal) THEN
       WRITE(6,fmt='(A)') 'NOT CALLING RANDOM PARAMETERS2'
       WRITE(6,fmt='(A)') 'This routine is only called every 3hrs'
   END IF
  END IF
END IF
!  ----- UM Section 35---- Stochastic Physics Random Parameters END ---
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS Stochastic_Phys',6)

! ---------------------------------------------------------------
!    diagnostic printing at beginning of timestep 1
!    (which is really timestep 0 values so lets record it as such)
! ---------------------------------------------------------------
       if( L_diag_print .and. timestep_number ==1 ) then
         CALL Print_diag_4A(U, V, THETA, RHO_r_sq_n, W, Q, QCL, QCF,    &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  offx, offy,timestep_number-1,                   &
     &                  print_step, diag_interval,                      &
     &                  rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
     &                  L_print_pe, L_print_w,                          &
     &                  L_print_wmax, L_print_max_wind,                 &
     &                  L_print_div, L_print_lapse, L_print_theta1,     &
     &                  L_print_shear, L_diag_wind, L_diag_noise,       &
     &                  max_w_run, max_wind_run, min_theta1_run,        &
     &                  dtheta1_run, max_div_run, min_div_run,          &
     &                  min_lapse_run, max_shear_run, time_max_shear,   &
     &                  time_div_max, time_div_min, time_lapse_min,     &
     &                  time_w_max, time_max_wind, time_theta1_min,     &
     &                  max_KE_run, min_KE_run, max_noise_run,          &
     &                  time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       endif     !  timestep_number ==1

!=== Polar filter + diffusion section ==================================
IF ( .NOT.first_atmstep_call ) THEN
! ----------------------------------------------------------------------
! Section 0.4  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------


      If ( L_filter ) then
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Filter',5)        

        IF (printstatus >= prstatus_diag ) THEN
          IF ( mype == 0 ) THEN
            WRITE(6,fmt='(A)') 'Calling polar filter routine'
            IF ( L_pofil_new )  WRITE(6,fmt='(A)')'Using pofil_new'
          END IF    
! section 13:
!         IF( SF(0,13) ) THEN    ! Diagnostics required for this section
! 
! ! Allocate diagnostic space for STASH
!           ALLOCATE (STASHwork13(STASH_maxlen(13,A_im)))
!         
!         END IF

        max_u = MAXVAL(u(udims%i_start:udims%i_end,                     &
                         udims%j_start:udims%j_end,                     &
                         udims%k_start:udims%k_end))
        min_u = MINVAL(u(udims%i_start:udims%i_end,                     &
                         udims%j_start:udims%j_end,                     &
                         udims%k_start:udims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_u)
        CALL gc_rmin(1,n_proc,ierr,min_u)

        max_v = MAXVAL(v(vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end,                     &
                         vdims%k_start:vdims%k_end))
        min_v = MINVAL(v(vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end,                     &
                         vdims%k_start:vdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_v)
        CALL gc_rmin(1,n_proc,ierr,min_v)

        max_w = MAXVAL(w(wdims%i_start:wdims%i_end,                     &
                         wdims%j_start:wdims%j_end,                     &
                         wdims%k_start:wdims%k_end))
        min_w = MINVAL(w(wdims%i_start:wdims%i_end,                     &
                         wdims%j_start:wdims%j_end,                     &
                         wdims%k_start:wdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_w)
        CALL gc_rmin(1,n_proc,ierr,min_w)

        max_thetav = MAXVAL(thetav(tdims%i_start:tdims%i_end,           &
                                   tdims%j_start:tdims%j_end,           &
                                   tdims%k_start:tdims%k_end))
        min_thetav = MINVAL(thetav(tdims%i_start:tdims%i_end,           &
                                   tdims%j_start:tdims%j_end,           &
                                   tdims%k_start:tdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_thetav)
        CALL gc_rmin(1,n_proc,ierr,min_thetav)


          IF (mype == 0) THEN
            WRITE(6,fmt='(A)') '==================================================='
            WRITE(6,fmt='(A)') 'Calling eg_NI_filter_Ctl'
            WRITE(6,fmt='(A)') ' Max/Min before polar filter :'
            WRITE(6,fmt='(A,2E25.10)') '    u = ',max_u, min_u
            WRITE(6,fmt='(A,2E25.10)') '    v = ',max_v, min_v
            WRITE(6,fmt='(A,2E25.10)') '    w = ',max_w, min_w
            WRITE(6,fmt='(A,2E25.10)') 'theta = ',max_thetav, min_thetav
          END IF
        END IF  

       Call eg_NI_filter_Ctl(  thetav, u, v, w, etadot, exner,          &
     &                      exner_theta_levels,                         &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v,                             &
     &                      max_121_rows, u_sweeps, v_sweeps,           &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      top_filt_start, top_filt_end,               &
     &                      up_diff, max_updiff_levels,                 &
     &                      horizontal_level,                           &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc_y, at_extremity, model_domain,        &
     &                      L_pftheta, L_pfuv,                          &
     &                      L_pfw, L_pfexner, L_diff_exner,             &
     &                      L_diff_thermo, L_diff_wind, L_diff_w,       &
     &                      L_pofil_hadgem2, Ltimer,                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &                      xi1_u, xi1_p, xi2_p, xi2_v,                 & 
     &                      pole_consts, gc_proc_row_group,             &
     &                      nproc, gc_proc_col_group,                   &
     &                      global_row_length,                          &
     &                      csxi2_v, csxi2_p, delta_lambda, delta_phi,  &
     &                      Earth_Radius,                               &
     &                      STASHwork13)

IF (integrity_test)                                                     &
CALL update_hash_m(                                                     &
                  u,               SIZE(u),               'u____',      &
                  v,               SIZE(v),               'v____',      &
                  w,               SIZE(w),               'w____')

        IF (printstatus >= prstatus_diag) THEN
          max_u = MAXVAL(u(udims%i_start:udims%i_end,                   &
                           udims%j_start:udims%j_end,                   &
                           udims%k_start:udims%k_end))
          min_u = MINVAL(u(udims%i_start:udims%i_end,                   &
                           udims%j_start:udims%j_end,                   &
                           udims%k_start:udims%k_end))
          CALL gc_rmax(1,n_proc,ierr,max_u)
          CALL gc_rmin(1,n_proc,ierr,min_u)

          max_v = MAXVAL(v(vdims%i_start:vdims%i_end,                     &
                           vdims%j_start:vdims%j_end,                     &
                           vdims%k_start:vdims%k_end))
          min_v = MINVAL(v(vdims%i_start:vdims%i_end,                     &
                           vdims%j_start:vdims%j_end,                     &
                           vdims%k_start:vdims%k_end))
          CALL gc_rmax(1,n_proc,ierr,max_v)
          CALL gc_rmin(1,n_proc,ierr,min_v)

          max_w = MAXVAL(w(wdims%i_start:wdims%i_end,                     &
                           wdims%j_start:wdims%j_end,                     &
                           wdims%k_start:wdims%k_end))
          min_w = MINVAL(w(wdims%i_start:wdims%i_end,                     &
                           wdims%j_start:wdims%j_end,                     &
                           wdims%k_start:wdims%k_end))
          CALL gc_rmax(1,n_proc,ierr,max_w)
          CALL gc_rmin(1,n_proc,ierr,min_w)

          max_thetav = MAXVAL(thetav(tdims%i_start:tdims%i_end,           &
                                     tdims%j_start:tdims%j_end,           &
                                     tdims%k_start:tdims%k_end))
          min_thetav = MINVAL(thetav(tdims%i_start:tdims%i_end,           &
                                     tdims%j_start:tdims%j_end,           &
                                     tdims%k_start:tdims%k_end))
          CALL gc_rmax(1,n_proc,ierr,max_thetav)
          CALL gc_rmin(1,n_proc,ierr,min_thetav)


          IF (mype == 0) THEN
            WRITE(6,fmt='(A)') ' '
            WRITE(6,fmt='(A)') '==================================================='
            WRITE(6,fmt='(A)') ' Max/Min after polar filter :'
            WRITE(6,fmt='(A,2E25.10)') '    u = ',max_u, min_u
            WRITE(6,fmt='(A,2E25.10)') '    v = ',max_v, min_v
            WRITE(6,fmt='(A,2E25.10)') '    w = ',max_w, min_w
            WRITE(6,fmt='(A,2E25.10)') 'theta = ',max_thetav, min_thetav
          END IF
        END IF

      if(L_pfexner .and. L_pofil_new)then
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure ( exner,                           &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           r_theta_levels, r_rho_levels, rho,                     &
     &           p, pstar, p_theta_levels,exner_theta_levels)
      endif  !  (L_pfexner .and. l_pofil_new)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Filter',6)        
      End If    !  L_filter

END IF ! timestep_number > 1 -> filter
!=== End Polar filter + diffusion section ===================================


! ----------------------------------------------------------------------
! Section 1.0  Call Atmospheric Physics1
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Atmos_Phys1 (AP1)',5)

! Biogenic aerosol climatology for the climate and NWP models
IF (l_use_biogenic) THEN
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        biogenic(i,j,k) = arclbiog_bg(i,j,k)
      END DO
    END DO
  END DO
END IF

! Aerosol climatologies - Model switches and climatologies are
! gathered into bigger arrays.

! First set the internal model switches according to the
! CNTLATM settings, and determine how many components we need.
CALL set_arcl_dimensions(l_use_arclbiom, l_use_arclblck,                &
                         l_use_arclsslt, l_use_arclsulp,                &
                         l_use_arcldust, l_use_arclocff,                &
                         l_use_arcldlta, n_arcl_species,                &
                         n_arcl_compnts, l_use_arcl )

! If the aerosol climatology for NWP is used, n_arcl_species
! is larger than 0. In that case, allocate the array gathering
! component mass-mixing ratio and take the values from the
! arrays in arg_atm_fields.h.

IF (n_arcl_species > 0) THEN
  ALLOCATE(arcl(row_length, rows, model_levels, n_arcl_compnts))

  CALL set_arcl_clim(row_length, rows, model_levels,                    &
                     n_arcl_compnts,                                    &
                     ! Internal model switches
                     l_use_arcl,                                        &
                     ! Climatologies from ancillary files
                     arclbiom_fr, arclbiom_ag, arclbiom_ic,             &
                     arclblck_fr, arclblck_ag,                          &
                     arclsslt_fi, arclsslt_jt,                          &
                     arclsulp_ac, arclsulp_ak, arclsulp_di,             &
                     arcldust_b1, arcldust_b2, arcldust_b3,             &
                     arcldust_b4, arcldust_b5, arcldust_b6,             &
                     arclocff_fr, arclocff_ag, arclocff_ic,             &
                     arcldlta_dl,                                       &
                     ! Internal climatology array
                     arcl,                                              &
                     ! Component array indices
                     i_arcl_compnts )

ELSE
   ALLOCATE ( arcl(1,1,1,1) )
END IF

!
! UKCA_RADAER: Obtain current UKCA setup and allocate arrays.
!              Most of the work needs only be made for the first
!              Atm_Step_4A() call.
!
IF (ErrorStatus == 0) THEN

  IF (l_ukca_radaer) THEN

    IF (first_atmstep_call) THEN

! DEPENDS ON: ukca_radaer_init
      CALL ukca_radaer_init(                                            &
                          ErrorStatus,                                  &
                          Cmessage,                                     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                          ukca_radaer )

      IF (ErrorStatus /= 0) THEN
        CALL Ereport("ATM_STEP_4A", ErrorStatus, Cmessage)
      END IF

    END IF

!
! Allocate those arrays that will receive UKCA output.
!
    ALLOCATE(ukca_radaer%mix_ratio(row_length, rows,                    &
                       model_levels, ukca_radaer%n_cpnt))
    ALLOCATE(ukca_radaer%comp_vol(row_length, rows,                     &
                       model_levels, ukca_radaer%n_cpnt))

    ALLOCATE(ukca_radaer%dry_diam(row_length, rows,                     &
                             model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%wet_diam(row_length, rows,                     &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_rho(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_wtv(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_vol(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_nbr(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))

!
! Get current UKCA results
!
! DEPENDS ON: ukca_radaer_get
    CALL ukca_radaer_get(                                               &
                       ErrorStatus,                                     &
                       Cmessage,                                        &
                       first_atmstep_call,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                       ukca_radaer  )

    IF (ErrorStatus /= 0) THEN

      DEALLOCATE(ukca_radaer%mix_ratio)
      DEALLOCATE(ukca_radaer%comp_vol)
      DEALLOCATE(ukca_radaer%dry_diam)
      DEALLOCATE(ukca_radaer%wet_diam)
      DEALLOCATE(ukca_radaer%modal_rho)
      DEALLOCATE(ukca_radaer%modal_wtv)
      DEALLOCATE(ukca_radaer%modal_vol)
      DEALLOCATE(ukca_radaer%modal_nbr)

      CALL Ereport("ATM_STEP_4A", ErrorStatus, Cmessage)

    END IF

  ELSE

! l_ukca_radaer is false. Set array dimensions
! to minimum value and allocate minimum space.
    ukca_radaer%n_mode = 1
    ukca_radaer%n_cpnt = 1

    ALLOCATE(ukca_radaer%mix_ratio(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%comp_vol(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%dry_diam(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%wet_diam(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_rho(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_wtv(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_vol(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_nbr(1, 1, 1, 1))

  END IF ! l_ukca_radaer

END IF ! ErrorStatus

! UKCA_CDNC: Allocate array and obtain Cloud Droplet No. Concentration from D1

      IF (ErrorStatus == 0) THEN

        IF (L_ukca .AND. (l_ukca_aie1 .OR. l_ukca_aie2)) THEN

! Allocate arrays to receive UKCA output.
          ALLOCATE(ukca_cdnc%cdnc(row_length, rows, model_levels))
          ALLOCATE(ukca_cdnc%cdnc3(row_length, rows, model_levels))

! Get current UKCA results
          CALL ukca_cdnc_get(                                           &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                             first_atmstep_call,                        &
                             ErrorStatus,                               &
                             Cmessage)

          IF (ErrorStatus /= 0) THEN
            DEALLOCATE(ukca_cdnc%cdnc)
            DEALLOCATE(ukca_cdnc%cdnc3)
            CALL Ereport("ATM_STEP_4A", ErrorStatus, Cmessage)
          END IF

! Set minimum CDNC value (equivalent to 5 cm-3)  
          WHERE(ukca_cdnc%cdnc < 5.0e06) ukca_cdnc%cdnc = 5.0e06 

          cdnc_dim1 = row_length  
          cdnc_dim2 = rows  
          cdnc_dim3 = model_levels

        ELSE
          ! both L_UKCA_AIE1 and L_UKCA_AIE2 are false.
          ALLOCATE(ukca_cdnc%cdnc(1,1,1))
          ukca_cdnc%cdnc(:,:,:) = 0.0
          ALLOCATE(ukca_cdnc%cdnc3(1,1,1))
          ukca_cdnc%cdnc3(:,:,:) = 0.0

          cdnc_dim1 = 1
          cdnc_dim2 = 1
          cdnc_dim3 = 1

        END IF ! L_UKCA_AIE1 or L_UKCA_AIE2

      END IF ! ErrorStatus = 0


IF (l_physics .AND. errorstatus == 0) THEN

! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init(                                              &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, rho_km, ch, exner_lbc_real_tend,      &
     w_lbc_real_tend, errorstatus, 'tropinit')

! DEPENDS ON: o3_to_3d
  CALL o3_to_3d(lexpand_ozone, i_ozone_int,                             &
                rows, row_length, model_levels, ozone_levels,           &
                halo_i, halo_j, offx, offy, at_extremity,               &
                a_realhd(rh_z_top_theta),                               &
                theta(tdims_s%i_start,tdims_s%j_start,1),               &
                exner_theta_levels(tdims_s%i_start,tdims_s%j_start,1),  &
                rho_r_sq_n,                                             &
                exner_rho_levels,                                       &
                nd_o3, o3,                                              &
                min_trop_level, max_trop_level,                         &
                l_o3_trop_level,l_o3_trop_height,                       &
                l_t_trop_level,l_t_trop_height,                         &
                o3_trop_level,o3_trop_height,                           &
                t_trop_level,t_trop_height,                             &
                gc_proc_row_group,                                      &
                global_row_length,                                      &
                ozone3d(o3dims2%i_start,o3dims2%j_start,1),             &
                errorstatus, cmessage)

  ozone3d(:,:,0)=ozone3d(:,:,1)

! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init(                                              &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, rho_km, ch, exner_lbc_real_tend,      &
     w_lbc_real_tend, errorstatus, 'ozoninit')

END IF !  L_Physics

! DEPENDS ON: Atm_Step_phys_init
CALL atm_step_phys_init(                                                &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, rho_km, ch, exner_lbc_real_tend,      &
     w_lbc_real_tend, errorstatus, 'microphy')

IF (l_physics .AND. errorstatus == 0) THEN

! DEPENDS ON: Atm_Step_diag
  CALL atm_step_diag(                                                   &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star,          &
       cf_star, cfl_star, cff_star, 1)



! this was commented out because ENDGame mixing ratios have already 
! been converted/stored in q. Since this is a mix of 
! allocate/copy/convert this caused an error of not allocated rhst. 
! Now the copy/conversion has been removed from atm_step_alloc_4A 
! when ENDGAME is defined and the allocation can go ahead

! only some code for pc2 in atm_step_alloc_4A is needed
! Convert to mixing ratios from specific humidities if needed
! DEPENDS ON: Atm_Step_alloc_4A
  CALL Atm_Step_alloc_4A(                                               &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
            cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,    &
            frac_control, r_u, r_v, r_w, errorstatus,'MixRatio')

! NB if you are changing the argument list to atmos_physics1, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent.

! NOTE: R_u and  R_v have changed but scm_main has been left untouched
!

! for debugging only, because we triggered the NaN detection in a CRUN below

    q_star  (:,:,:) = 0.
    qcl_star(:,:,:) = 0.
    qcf_star(:,:,:) = 0.

    IF (l_mcr_qcf2  )  qcf2_star  (:,:,:) = 0.
    IF (l_mcr_qrain )  qrain_star (:,:,:) = 0.
    IF (l_mcr_qgraup)  qgraup_star(:,:,:) = 0.
    
! end debugging only

    IF (L_cosp) THEN
!     This is done until a prognostic is developed
      IF (timestep_number == 1) THEN
        nullify(cosp_crain_3d,cosp_csnow_3d)
      END IF
      IF (.not. associated(cosp_crain_3d)) THEN
        IF (PrintStatus >= PrStatus_Oper) THEN
          WRITE(6,'(A,I5)')                                          &
            'COSP in ATM_STEP_4A: allocating cosp_crain_3d in tstep ', &
            timestep_number
        END IF
        ALLOCATE(cosp_crain_3d(row_length,rows,model_levels))
        cosp_crain_3d = 0.0
      END IF
      IF (.not. associated(cosp_csnow_3d)) THEN
        IF (PrintStatus >= PrStatus_Oper) THEN
          WRITE(6,'(A,I5)')                                          &
            'COSP in ATM_STEP_4A: allocating cosp_csnow_3d in tstep ', &
            timestep_number
        END IF
        ALLOCATE(cosp_csnow_3d(row_length,rows,model_levels))
        cosp_csnow_3d = 0.0
      END IF
    ELSE
      IF (.not. associated(cosp_crain_3d)) THEN
         ALLOCATE(cosp_crain_3d(1,1,1))
         cosp_crain_3d = 0.0
      END IF
      IF (.not. associated(cosp_csnow_3d)) THEN
         ALLOCATE(cosp_csnow_3d(1,1,1))
         cosp_csnow_3d = 0.0
      END IF
    END IF
 !CABLE: fudge - this is set in atm_step_phys_init BUt this is called after
   dim_cs2 = LAND_FIELD
   CALL cable_atm_step( mype, & 
            .TRUE., &   ! UM_atm_step=
            .TRUE., & ! fudge L_cable switch, 
            TIMESTEP_NUMBER, &
            TIMESTEP, &
            row_length,     &
            rows, &
            LAND_FIELD, ntiles, sm_levels, dim_cs1, dim_cs2,              &
            true_LATITUDE, true_LONGITUDE )!,          &
            !land_index, &
            !clapp_horn, & ! bexp, &
            !therm_cond, & !hcon, &
            !SAT_SOIL_COND, & ! satcon, &
            !SAT_SOILW_SUCTION, & ! sathh,       &
            !VOL_SMC_sat, & ! smvcst, &
            !VOL_SMC_WILT, & ! smvcwt, &
            !VOL_SMC_crit, & ! smvccl, &
            !soil_alb, & ! albsoil, &
            !lw_down, &
            !cos_zenith_angle, &
            !ls_rain, ls_snow, pstar, CO2_MMR,         &
            !sthu, smcl, sthf, GS, &
            !canopy_water, land_alb )


! DEPENDS ON: atmos_physics1
        Call Atmos_Physics1(                                             &
! Parallel variables
     & global_row_length, global_rows, nproc, nproc_x ,nproc_y           &
     &,g_rows, g_row_length                                              &

! model dimensions
     &, row_length, rows, n_rows, land_field, model_levels, wet_levels   &
     &, bl_levels, st_levels, sm_levels, Ozone_levels, cloud_levels      &
     &, land_ice_points, soil_points, n_cca_lev, ntiles, nice_use        &
     &, salt_dim1, salt_dim2, salt_dim3, tr_levels, tr_ukca              &
      , cdnc_dim1, cdnc_dim2, cdnc_dim3                                  &
     &, co2_dim_len, co2_dim_row, co2_dim_lev                            &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                   &

! model switches
     &,     L_regular, L_lbc_old                                         &
     &,     L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                 &
     &,     LCAL360, Ltimer                                              &
     &,     L_mr_physics1                                                &
     &,     L_ukca_chem, L_ukca_useumuivals, L_ukca_set_trace_gases      &
     &,     L_ukca_strat, L_ukca_strattrop                               &
     &,     L_ukca_prescribech4, L_USE_ARCL                              &

! model Parameters
     &,     RHcrit                                                       &
     &,min_trop_level, max_trop_level                                    &
     
! Pass position of greenhouse gases in tracer_ukca array, for chemical coupling
     &,     ngrgas ,grgas_addr, Ntot_land, Ntot_sea                      &
! parameter for stochastic physics random parameters2
     &,     M_CI                                                         &
! Vertical coordinate levels.
     &, delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                 &
! Time stepping information
     &, I_year, I_day_number, I_hour, I_minute                           &
     &, I_second, PREVIOUS_TIME,                                         &
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &      STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14      &
! Additional variables for SCM diagnostics
     &,     nSCMDpkgs, L_SCMDiags                                        &
! Data Fields.
     &,     THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, RHO_r_sq_n, U, V    &
     &,     P, PSTAR, EXNER, EXNER_THETA_LEVELS, LAND                    &
     &,     P_THETA_LEVELS, FRAC_LAND,frac_control                       &
      ,     ukca_cdnc%cdnc                                               &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,land_index, RGRAIN_TILE,SNSOOT,NTML, CUMULUS                      &
     &,ice_fraction,p_ice_fract_rad, p_ice_thick_rad                     &
     &,CCA, CCB, CCT, CCLWP, CCW_RAD, LCBASE                             &
     &,TSTAR,TSTAR_LAND,TSTAR_SEA,p_tstar_sice                           &
     &,SICE_ALB, LAND_ALB, SNODEP,p_snodep_sice                          &
     &,ozone3D(o3dims2%i_start,o3dims2%j_start,1)                        &
     &, SW_INCS, LW_INCS, DIRPAR                                         &
     &,O3_trop_level, O3_trop_height, T_trop_level, T_trop_height        &
     &,ZH, OROG_SD, OROG_GRAD_XX, OROG_GRAD_XY                           &
     &,OROG_GRAD_YY, CF_AREA, CF_BULK,CF_LIQUID,CF_FROZEN,MURK_SOURCE    &
      ,arcl, soil_alb, obs_alb_sw, obs_alb_vis, obs_alb_nir, lai_pft     & 
      ,snodep_tile, frac_typ,tstar_tile,z0_tile,dolr_field,lw_down       & 
      ,sw_tile_rts,es_space_interp, rad_mask, cos_zenith_angle           & 
! Variables for COSP
      ,cosp_crain_3d,cosp_csnow_3d                                       &
! IN JULES 2 prognostics
      ,snowdepth,lake_h_ice                                              &
! IN/OUT
     &,theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star     &
     &,qgraup_star, cf_star, cfl_star, cff_star, R_u, R_v                &
     &,a_realhd(rh_energy_corr), NET_FLUX, NET_MFLUX                     &
     &,MURK, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5       &
     &,DUST_DIV6, SO2, SO4_AITKEN,SO4_ACCU, SO4_DISS, NH3, SOOT_NEW      &
     &,SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, OCFF_NEW,    &
     & OCFF_AGD, OCFF_CLD, NITR_ACC, NITR_DISS, CO2, TRACER_UKCA         &
     &,biogenic, A_INTHD(23), ukca_radaer                                &
! OUT Fields
     &,     ls_rain, ls_snow, micro_tends, rho                           &
      ,     photosynth_act_rad, rad_hr, dolr, SW_tile                    &
! error information
     &,     ErrorStatus  )


  ! temporarily store r_u etc in the _p2 arrays. 
  ! the _p2 arrays are not used at this point in the code, only
  ! subsequently to eg_f1sp, at which point the stored value is no
  ! longer required. 
  r_u_p2=r_u
  r_v_p2=r_v
  r_w_p2=r_w

! set zero level star quantities for Endgame
! currently set equal to level 1
  CALL set_star_zero_level(                                              &
             theta_star,                                                 &
             q_star,                                                     &
             qcl_star,                                                   &
             qcf_star,                                                   &
             cf_star,                                                    &
             cfl_star,                                                   &
             cff_star,                                                   &
             qcf2_star,                                                  &
             qrain_star,                                                 &
             qgraup_star,                                                &
             L_mcr_qgraup,                                               &
             L_mcr_qrain,                                                &
             L_mcr_qcf2)


  IF (l_mr_physics1) THEN

    m_star  (:,:,:) = q_star  (:,:,:)
    mcl_star(:,:,:) = qcl_star(:,:,:)
    mcf_star(:,:,:) = qcf_star(:,:,:)

    IF (l_mcr_qcf2  )    mcf2_star  (:,:,:) = qcf2_star  (:,:,:)
    IF (l_mcr_qrain )    mrain_star (:,:,:) = qrain_star (:,:,:)
    IF (l_mcr_qgraup)    mgraup_star(:,:,:) = qgraup_star(:,:,:)

! Commented out because ENDGame does not use qstore (yet).
! Reinstate q's by converting from ENDGames mixing ratios.
! Copy contents of q_store (sp hums before atmos_physics1) back into
! the d1(q) variables etc. Include the halo points.
! store physics changes for use in Sl_moist_conserve
! DEPENDS ON: Atm_Step_alloc_4A
!          CALL Atm_Step_alloc_4A( &
!             cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
!             frac_control, r_u, r_v, r_w, errorstatus, 'cpqstore' )

! only needed if we aren't using mr in phys2
  IF (.NOT. l_mr_physics2) THEN
    CALL eg_mix_to_q                                                    &
                  (tdims_l,tdims_s,                                     &
                   m_v, m_cl, m_cf,                                     &
                   m_cf2, m_r, m_gr,                                    &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup)
  END IF


! NOTE: not using store flavours because they were not used to start with.

    ELSE

! temporarily convert _star fields from increments to full fields

    q_star(:,:,:)    = q(1-offx:row_length+offx,1-offy:rows+offy,:)     &
                     + q_star(:,:,:)
    qcl_star(:,:,:)  = qcl(1-offx:row_length+offx,1-offy:rows+offy,:)   &
                     + qcl_star(:,:,:)
    qcf_star(:,:,:)  = qcf(1-offx:row_length+offx,1-offy:rows+offy,:)   &
                     + qcf_star(:,:,:)
    IF (l_mcr_qcf2) THEN 
      qcf2_star(:,:,:)=qcf2(1-offx:row_length+offx,1-offy:rows+offy,:)  &
                      + qcf2_star(:,:,:)
    END IF
    IF (l_mcr_qrain) THEN
      qrain_star(:,:,:)=qrain(1-offx:row_length+offx,1-offy:rows+offy,:)&
                       + qrain_star(:,:,:) 
    END IF
    IF (l_mcr_qgraup) THEN
      qgraup_star(:,:,:)=qgraup(1-offx:row_length+offx,1-offy:rows+offy,:)&
                        + qgraup_star(:,:,:)
    END IF

    CALL eg_q_to_mix                                                    &
                  (tdims_s,tdims_s,                                     &
                   q_star, qcl_star, qcf_star,                          &
                   qcf2_star, qrain_star, qgraup_star,                  &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   m_star, mcl_star, mcf_star                           &
                  ,mcf2_star, mrain_star, mgraup_star)

! now convert back to increments
    q_star(:,:,:)    = q_star(:,:,:)                                    &
                     - q(1-offx:row_length+offx,1-offy:rows+offy,:)
    m_star = m_star - m_v
    qcl_star(:,:,:)  = qcl_star(:,:,:)                                  &
                     - qcl(1-offx:row_length+offx,1-offy:rows+offy,:)
    mcl_star = mcl_star - m_cl
    qcf_star(:,:,:)  = qcf_star(:,:,:)                                  &
                     - qcf(1-offx:row_length+offx,1-offy:rows+offy,:)
    mcf_star = mcf_star - m_cf
    IF (l_mcr_qcf2) THEN
      qcf2_star(:,:,:) = qcf2_star(:,:,:)                               &
                       - qcf2(1-offx:row_length+offx,1-offy:rows+offy,:)
      mcf2_star = mcf2_star - m_cf2
    END IF
    IF (l_mcr_qrain) THEN
      qrain_star(:,:,:)= qrain_star(:,:,:)                              &
                      - qrain(1-offx:row_length+offx,1-offy:rows+offy,:)
      mrain_star = mrain_star - m_r
    END IF
    IF (l_mcr_qgraup) THEN
      qgraup_star(:,:,:)= qgraup_star(:,:,:)                            &
                     - qgraup(1-offx:row_length+offx,1-offy:rows+offy,:)
      mgraup_star = mgraup_star - m_gr
    END IF

  END IF

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)
! DEPENDS ON: Atm_Step_stash
  CALL atm_step_stash( &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       errorstatus, 1)
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',6)

 ELSE  ! L_physics =.false.

! initialise arrays that hold physics increments to zero
! DEPENDS ON: Atm_Step_phys_init
   CALL atm_step_phys_init( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,  &
     cf_star, cfl_star, cff_star, rho_km, ch, exner_lbc_real_tend, &
     w_lbc_real_tend, errorstatus, 'zeroincs')

! initial NULL change atmos_physics1 creates tendencies. So if 
! nothing happens all tendencies are zero.
  theta_star(:,:,:)  = 0.
  q_star    (:,:,:)  = 0.
  qcl_star  (:,:,:)  = 0.
  qcf_star  (:,:,:)  = 0.
  IF (l_mcr_qcf2)  qcf2_star   (:,:,:) = 0.
  IF (l_mcr_qrain)  qrain_star (:,:,:) = 0.
  IF (l_mcr_qgraup) qgraup_star(:,:,:) = 0.
! NOT running any physics code. The stub below initialises starred
! state variables for computation of the source terms.
! This is a "unnecessary" memory copy _provided_ that DR's
! understanding of the physics increments is correct.
!
  CALL eg_idl_forcing(                                                  &
! IN Data Fields.
       u, v, theta, exner_theta_levels, exner,                          &
! IN/OUT
       theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star,   &
       qgraup_star , r_u, r_v,                                          &
       l_expl_horz_drag, L_HeldSuarez, L_HeldSuarez1_drag,              &
       L_HeldSuarez2_drag ,                                             &
! error information
       errorstatus  )

  m_star  (:,:,:) = q_star  (:,:,:)
  mcl_star(:,:,:) = qcl_star(:,:,:)
  mcf_star(:,:,:) = qcf_star(:,:,:)

  IF (l_mcr_qcf2  )    mcf2_star  (:,:,:) = qcf2_star  (:,:,:)
  IF (l_mcr_qrain )    mrain_star (:,:,:) = qrain_star (:,:,:)
  IF (l_mcr_qgraup)    mgraup_star(:,:,:) = qgraup_star(:,:,:)

END IF ! L_Physics

! Deallocate the arrays of UKCA_RADAER
  DEALLOCATE(ukca_radaer%mix_ratio)
  DEALLOCATE(ukca_radaer%comp_vol)
  DEALLOCATE(ukca_radaer%dry_diam)
  DEALLOCATE(ukca_radaer%wet_diam)
  DEALLOCATE(ukca_radaer%modal_rho)
  DEALLOCATE(ukca_radaer%modal_wtv)
  DEALLOCATE(ukca_radaer%modal_vol)
  DEALLOCATE(ukca_radaer%modal_nbr)

! Deallocate arrays for UKCA_CDNC
  DEALLOCATE(ukca_cdnc%cdnc)
  DEALLOCATE(ukca_cdnc%cdnc3)


! Deallocate the array of the aerosol climatology for NWP
DEALLOCATE(arcl)

IF (printstatus >=  prstatus_diag) THEN

  max_r_u     = MAXVAL(R_u(udims%i_start:udims%i_end,                   &
                           udims%j_start:udims%j_end,                   &
                           udims%k_start:udims%k_end))
  min_r_u     = MINVAL(R_u(udims%i_start:udims%i_end,                   &
                           udims%j_start:udims%j_end,                   &
                           udims%k_start:udims%k_end))
  max_r_v     = MAXVAL(R_v(vdims%i_start:vdims%i_end,                   &
                           vdims%j_start:vdims%j_end,                   &
                           vdims%k_start:vdims%k_end))
  min_r_v     = MINVAL(R_v(vdims%i_start:vdims%i_end,                   &
                           vdims%j_start:vdims%j_end,                   &
                           vdims%k_start:vdims%k_end))

  CALL gc_rmax(1,n_proc,ierr,max_r_u)
  CALL gc_rmax(1,n_proc,ierr,max_r_v)
  CALL gc_rmin(1,n_proc,ierr,min_r_u)
  CALL gc_rmin(1,n_proc,ierr,min_r_v)

  IF(couple_app) THEN
    IF(mype==0) THEN
      WRITE(6,fmt='(A)') '============================================   &
                          &========================================'

      WRITE(6,fmt='(A)') 'Slow physics source terms from atmos_physics1:'
      WRITE(6,fmt='(A,2E32.16)') 'r_u      :',min_r_u,max_r_u
      WRITE(6,fmt='(A,2E32.16)') 'r_v      :',min_r_v,max_r_v
    ENDIF
  ELSE
    r_u = 0.
    r_v = 0.
  END IF
END IF


IF (INTEGRITY_TEST)                                                   &
  CALL update_hash_m(r_u,   SIZE(r_u),         'R_u__',               &
                     r_v,   SIZE(r_v),         'R_v__')

! compute the slow physics source terms for theta and moisture
CALL eg_r( theta_star, m_star,mcl_star,mcf_star,mgraup_star, &
           mrain_star,mcf2_star, m_v, theta, couple_app )

CALL eg_sisl_init(                                                    &
         row_length,rows,n_rows,model_levels, l_inc_solver,           &
         l_call_from_solver,.FALSE., ih,g_theta, u, v, w,             &
         thetav, rho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2,exner,        &
         exner_surf, r_u, r_v, r_w, r_theta, r_rho, r_m_v, r_m_cl,    &
         r_m_cf, r_m_r, r_m_gr, r_m_cf2, etadot,  psi_w_surf,psi_w_lid)

! DEPENDS ON: timer
IF (Ltimer) CALL TIMER('AS Atmos_Phys1 (AP1)',6)

IF (l_tracer .AND. errorstatus == 0) THEN

! store physics changes
  CALL tr_set_phys_4A(super_array_size, super_tracer_phys1,           &
                   l_co2_interactive, co2,                            &
                   l_murk_advect, murk,                               &
                   l_soot, soot_new, soot_agd, soot_cld,              &
                   l_sulpc_so2, so2, so4_aitken,                      &
                   so4_accu, so4_diss, l_sulpc_nh3, nh3,              &
                   l_sulpc_dms, dms,                                  &
                   l_dust, dust_div1, dust_div2,  dust_div3,          &
                           dust_div4, dust_div5,  dust_div6,          &
                   l_biomass, bmass_new, bmass_agd, bmass_cld,        &
                   l_ocff, ocff_new, ocff_agd, ocff_cld,              &
                   l_nitrate, nitr_acc, nitr_diss,                    &
                   l_use_cariolle, ozone_tracer,                      &
                   tracer_phys1, tracer,                              &
                   ukca_tracer_phys1, tracer_ukca,                    &
                   row_length, rows,                                  &
                   model_levels, tr_levels, tr_vars, tr_ukca,         &
                   offx, offy, model_domain,                          &
                   .TRUE., tdims_l )

END IF  ! L_tracer and ErrorStatus == 0


! store physics changes for use in Sl_moist_conserve
! only some pc2 code in atm_step_alloc_4A is needed
! DEPENDS ON: Atm_Step_alloc_4A
CALL atm_step_alloc_4A(                                                 &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
         cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,       &
         frac_control, r_u, r_v, r_w, errorstatus, 'alloc_mr' )


! DEPENDS ON: Atm_Step_phys_reset
CALL atm_step_phys_reset(                                               &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
    q_star, qcl_star, qcf_star, r_v, r_u, dolr, theta_star, gs1,        &
    cf_star, cfl_star, cff_star, 'cyclrset')


!=== Polar filter + diffusion of increments section ==================================
! ----------------------------------------------------------------------
! Section 0.4  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------

! Call timer for diffusion code

      If ( L_polar_filter_incs .or. L_filter_incs ) then
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Filter',5)
      IF (printstatus >= prstatus_diag) THEN
        max_r_theta     = MAXVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                         tdims%j_start:tdims%j_end,     &
                                         tdims%k_start:tdims%k_end))
        min_r_theta     = MINVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                         tdims%j_start:tdims%j_end,     &
                                         tdims%k_start:tdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_theta)
        CALL gc_rmin(1,n_proc,ierr,min_r_theta)

        max_r_u     = MAXVAL(r_u(udims%i_start:udims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        min_r_u     = MINVAL(r_u(udims%i_start:tdims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_u)
        CALL gc_rmin(1,n_proc,ierr,min_r_u)

        max_r_v     = MAXVAL(r_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        min_r_v     = MINVAL(r_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_v)
        CALL gc_rmin(1,n_proc,ierr,min_r_v)

        max_r_w     = MAXVAL(r_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        min_r_w     = MINVAL(r_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_w)
        CALL gc_rmin(1,n_proc,ierr,min_r_w)

        max_r_rho     = MAXVAL(r_rho(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        min_r_rho     = MINVAL(r_rho(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_rho)
        CALL gc_rmin(1,n_proc,ierr,min_r_rho)

        max_r_p_p     = MAXVAL(r_p_p(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        min_r_p_p     = MINVAL(r_p_p(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_p_p)
        CALL gc_rmin(1,n_proc,ierr,min_r_p_p)

        IF (mype == 0) THEN
          WRITE(6,fmt='(A)') '==================================================='
          WRITE(6,fmt='(A)') 'Calling eg_NI_filter_incs_Ctl'
          WRITE(6,fmt='(A)') ' Max/Min before polar filter :'
          WRITE(6,fmt='(A,2E25.10)') '    R_u = ',MAX_r_u,MIN_r_u
          WRITE(6,fmt='(A,2E25.10)') '    R_v = ',MAX_r_v,MIN_r_v
          WRITE(6,fmt='(A,2E25.10)') '    R_w = ',MAX_r_w,MIN_r_w
          WRITE(6,fmt='(A,2E25.10)') 'R_theta = ',MAX_r_theta,MIN_r_theta
          WRITE(6,fmt='(A,2E25.10)') '  R_rho = ',MAX_r_rho,MIN_r_rho
          WRITE(6,fmt='(A,2E25.10)') '  R_p_p = ',MAX_r_p_p,MIN_r_p_p
        END IF  
      END IF

      CALL eg_ni_filter_incs_ctl(                                       &
                            r_theta, r_u, r_v, r_w, r_rho, r_p_p,       &
                            s_thetav, s_u, s_v, s_w,                    &
                            thetav, u, v, w, rho, exner,                &
                            row_length, rows, n_rows, model_levels,     &
                            r_theta_levels, r_rho_levels,               &
                            r_at_u, r_at_v,                             &
                            max_121_rows, u_sweeps, v_sweeps,           &
                            global_u_filter, global_v_filter,           &
                            u_begin, u_end, v_begin, v_end,             &
                            diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
                            first_constant_r_rho_level,                 &
                            first_constant_r_rho_level_m1,              &
                            horizontal_level,                           &
                            offx, offy, halo_i, halo_j,                 &
                            nproc_y, at_extremity, model_domain,        &
                            L_polar_filter_incs, L_diff_incs,           &
                            L_pofil_hadgem2, Ltimer,                    &
                            l_eliminate_rho, 0,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                            xi1_u, xi1_p, xi2_p,                        & 
                            pole_consts, gc_proc_row_group,             &
                            global_row_length,                          &
                            csxi2_v, csxi2_p)

      IF (printstatus >= prstatus_diag) THEN
        max_r_theta     = MAXVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                         tdims%j_start:tdims%j_end,     &
                                         tdims%k_start:tdims%k_end))
        min_r_theta     = MINVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                         tdims%j_start:tdims%j_end,     &
                                         tdims%k_start:tdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_theta)
        CALL gc_rmin(1,n_proc,ierr,min_r_theta)

        max_r_u     = MAXVAL(r_u(udims%i_start:udims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        min_r_u     = MINVAL(r_u(udims%i_start:tdims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_u)
        CALL gc_rmin(1,n_proc,ierr,min_r_u)

        max_r_v     = MAXVAL(r_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        min_r_v     = MINVAL(r_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_v)
        CALL gc_rmin(1,n_proc,ierr,min_r_v)

        max_r_w     = MAXVAL(r_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        min_r_w     = MINVAL(r_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_w)
        CALL gc_rmin(1,n_proc,ierr,min_r_w)

        max_r_rho     = MAXVAL(r_rho(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        min_r_rho     = MINVAL(r_rho(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_rho)
        CALL gc_rmin(1,n_proc,ierr,min_r_rho)

        max_r_p_p     = MAXVAL(r_p_p(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        min_r_p_p     = MINVAL(r_p_p(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end,         &
                                     pdims%k_start:pdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_p_p)
        CALL gc_rmin(1,n_proc,ierr,min_r_p_p)

        IF (mype == 0) THEN
          WRITE(6,fmt='(A)') ' '
          WRITE(6,fmt='(A)') ' Max/Min after polar incs filter :'
          WRITE(6,fmt='(A,2E25.10)') '    R_u = ',MAX_r_u,MIN_r_u
          WRITE(6,fmt='(A,2E25.10)') '    R_v = ',MAX_r_v,MIN_r_v
          WRITE(6,fmt='(A,2E25.10)') '    R_w = ',MAX_r_w,MIN_r_w
          WRITE(6,fmt='(A,2E25.10)') 'R_theta = ',MAX_r_theta,MIN_r_theta
          WRITE(6,fmt='(A,2E25.10)') '  R_rho = ',MAX_r_rho,MIN_r_rho
          WRITE(6,fmt='(A,2E25.10)') '  R_p_p = ',MAX_r_p_p,MIN_r_p_p
        END IF      
      END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(r_theta, SIZE(r_theta),     'r_the',             &
                     r_u,     SIZE(r_u),         'R_u__',             &
                     r_v,     SIZE(r_v),         'R_v__',             &
                     R_w,     SIZE(R_w),         'R_w__',             &
                     R_rho,   SIZE(R_rho),       'R_r__',             &
                     R_p_p,   SIZE(R_p_p),       'Rpp__')


! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Filter',6)
      End If    ! L_polar_filter_incs .or. L_filter_incs

!=== End Polar filter + diffusion of increments section ===================

IF ( sf(188,30) ) THEN
  ALLOCATE ( inc_rho(pdims_s%i_start:pdims_s%i_end,            &
                     pdims_s%j_start:pdims_s%j_end,            &
                     pdims_s%k_start:pdims_s%k_end) )
ELSE
  ALLOCATE ( inc_rho( 1,1,1) )
END IF


! ----------------------------------------------------------------------
! Section 0.5 Calculation of coefficients in turbulence scheme
! ----------------------------------------------------------------------

! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
    r_v, r_u, theta_star, q_star, qcl_star, qcf_star,  &
    cf_star, cfl_star, cff_star, rho_km, ch, exner_lbc_real_tend, &
    w_lbc_real_tend, errorstatus, 'turb_cof')

IF (L_subfilter_horiz .OR. L_subfilter_vert) THEN  

! Calculate lambda^2 and S in eg_turb_smagorinsky
! DEPENDS ON: eg_turb_smagorinsky
    CALL eg_turb_smagorinsky(                                           &
                             u, v, w, z0,                               &
                             row_length, rows, n_rows, model_levels,    &
                             r_theta_levels, r_rho_levels,              &
                             xi1_p, xi1_u, xi2_p, xi2_v,                &
                             csxi2_p, csxi2_v )

END IF  !   L_subfilter_horiz or L_subfilter_vert

! ----------------------------------------------------------------------
! Section 1.0 Start outer loop
! ----------------------------------------------------------------------

IF (total_conv_outer) numcycles = 999

exner_res = 1. ! initialise to a "large" value
cycleno=0

outerloop : DO WHILE  (cycleno < numcycles)

cycleno=cycleno+1

IF(total_conv_outer .AND. exner_res < 2.d-3) cycleno = numcycles

!-----------------------------------------------------------------------

! Restore phys1 variables to be used as predictors
! DEPENDS ON: Atm_Step_phys_reset
  CALL atm_step_phys_reset(                                             &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
    q_star, qcl_star, qcf_star, r_v, r_u, dolr, theta_star, gs1,        &
    cf_star, cfl_star, cff_star, 'phy1rest')

!
! here we use the stored value of r_u etc (in the _p2 arrays) which
! have not seen the addition of dynamics increments
!
! DEPENDS ON: Atm_Step_diag
         CALL Atm_Step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
              r_v_p2, r_u_p2, r_w_p2, q_star, qcl_star, qcf_star,       &
              theta_star, cf_star, cfl_star, cff_star, 2)

! ----------------------------------------------------------------------
! Calculate time level n advective momentum quantities.
! ----------------------------------------------------------------------
IF (integrity_test)                                                   &
  CALL  update_hash_m( R_u,            SIZE(R_u),             'R_u__',&
                       R_v,            SIZE(R_v),             'R_v__')

! update LBC's for n+1 variables
      IF (model_domain /= mt_global) THEN
!
! Currently forced to match settings for the standard EG SEUKV job.
! This will change at a later date when enough LAM tests have been performed
! and the best settings have been determined!
!
      IF( .FALSE. .AND. cycleno > 1) THEN
!DEPENDS ON: init_lbc_dynamics
         CALL init_lbc_dynamics(u_np1,v_np1,w_np1, exner_np1,                  &
                                rho_np1, thetav_np1, etadot_np1,               &
                                m_v_np1, m_cl_np1, m_cf_np1,                   &
                                m_r_np1, m_cf2_np1, m_gr_np1,                  &
                                u_lbc, v_lbc, w_lbc, rho_lbc, exner_lbc,       &
                                theta_lbc, q_lbc, qcl_lbc, qcf_lbc,            &
                                qrain_lbc, qcf2_lbc, qgraup_lbc,               &
                                u_lbc_tend, v_lbc_tend, w_lbc_tend,            &
                                rho_lbc_tend, exner_lbc_tend,                  &
                                theta_lbc_tend, q_lbc_tend, qcl_lbc_tend,      &
                                qcf_lbc_tend, qrain_lbc_tend,                  &
                                qcf2_lbc_tend, qgraup_lbc_tend,                &
                                increment_factor, rim_stepsa, timestep_number, &
                                row_length, rows, n_rows, model_levels,        &
                                offx, offy, halo_i, halo_j,                    &
                                l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,         &
                                lenrima(1,1,rima_type_norm),                   &
                                rimwidtha(rima_type_norm), rimweightsa,        &
                                lbc_sizea(1,1,1,rima_type_norm),               &
                                lbc_startA(1,1,1,rima_type_norm),              &
                                at_extremity,                                  &
                                L_do_boundaries, L_do_halos                    &
                               )
      END IF
      END IF

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS S-L Advect (AA)',5)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Full_Wind',5)

  CALL eg_sl_full_wind(                                                 &
               row_length, rows, n_rows, model_levels, halo_i,          &
               halo_j, offx, offy, g_i_pe,depart_scheme,                &
               depart_order, high_order_scheme(wind_sl),                &
               monotone_scheme(wind_sl),                                &
               ritchie_high_order_scheme, ritchie_monotone_scheme,      &
               first_constant_r_rho_level,                              &
               interp_vertical_search_tol, check_bottom_levels,         &
               l_cartesian, l_shallow,                                  &
               l_rk_dps,                                                &
               l_high(wind_sl),  l_mono(wind_sl),  l_ritchie_high,      &
               l_ritchie_mono, lam_max_cfl,                             &
               etadot, u_np1, v_np1, w_np1, etadot_np1,                 &
               u_adv, v_adv, w_adv, r_u, r_v, r_w,                      &
               r_u_d, r_v_d, r_w_d,                                     &
               errorstatus )

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Full_Wind',6)

L_tracer_1_if :  IF ( L_tracer .AND. CycleNo == NumCycles ) THEN  !it1

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AA SL_Tracer',5)

!
! reset tracers to super_tracer_phys1 which contains tracers at level
! n and physics1 increments only. super_tracer_phys1 is the quantity
! to be advected. Furthermore SL_tracer1_4A could be modified to take in
! super_tracer_phys1 directly instead of this multiple packing and
! unpacking of tracers as super_tracer_phys1 remains unchanged during
! the time-step, irrespective of the NumCycles. 
!
      CALL unpack_tracer( super_array_size,                             &
                       super_tracer_phys1, tdims_l,                     &
                       co2, l_co2_interactive,                          &
                       murk, l_murk_advect,                             &
                       soot_new, soot_agd,                              &
                       soot_cld, l_soot,                                &
                       bmass_new, bmass_agd,                            &
                       bmass_cld, l_biomass,                            &
                       ocff_new, ocff_agd, ocff_cld, l_ocff,            &
                       dust_div1,dust_div2,                             &
                       dust_div3,dust_div4,                             &
                       dust_div5,dust_div6, l_dust,                     &
                       so2, so4_aitken,                                 &
                       so4_accu,                                        &
                       so4_diss, nh3, dms,                              &
                       l_sulpc_so2, l_sulpc_nh3, l_sulpc_dms,           &
                       nitr_acc, nitr_diss, l_nitrate,                  &
                       tracer, tr_levels, tr_vars,                      &
                       tracer_ukca, tr_ukca,                            &
                       l_use_cariolle, ozone_tracer, errorstatus)


      CALL sl_tracer1_4A(super_array_size, eta_theta_levels,            &
                       r_rho_levels, r_theta_levels,                    &
                       pstar, p, p_theta_levels,                        &
                       rho_r_sq_n, l_tracer1_non_hydro,                 &
                       row_length, rows, n_rows, model_levels,          &
                       halo_lam, halo_phi,                              &
                       fv_cos_theta_latitude,                           &
                       mype, nproc, nproc_x, nproc_y,                   &
                       halo_i, halo_j, datastart,                       &
                       g_i_pe, g_j_pe, l_2dcomm, size_2dcomm,           &
                       group_2dcomm, at_extremity,                      &
                       global_row_length, global_rows,                  &
                       gc_proc_row_group, gc_proc_col_group,            &
                       gc_all_proc_group, offx, offy,                   &
                       high_order_scheme(moist_sl),                     &
                       monotone_scheme(moist_sl),                       &
                       model_domain, l_high(moist_sl),                  &
                       l_mono(moist_sl),                                &
                       check_bottom_levels,                             &
                       interp_vertical_search_tol,                      &
                       first_constant_r_rho_level,                      &
                       depart_xi1_w, depart_xi2_w,  depart_xi3_w,       &
!
!  d1(jxxx) holds time level n value plus physics1 increment
                       co2, l_co2_interactive,                           &
                       murk, l_murk_advect,                              &
                       dust_div1,dust_div2,                              &
                       dust_div3,dust_div4,                              &
                       dust_div5,dust_div6, l_dust,                      &
                       soot_new, soot_agd,                               &
                       soot_cld, l_soot,                                 &
                       bmass_new, bmass_agd,                             &
                       bmass_cld, l_biomass,                             &
                       ocff_new, ocff_agd, ocff_cld, l_ocff,             &
                       so2, so4_aitken,                                  &
                       so4_accu,                                         &
                       so4_diss, nh3, dms,                               &
                       l_sulpc_so2, l_sulpc_nh3, l_sulpc_dms,            &
                       nitr_acc, nitr_diss, l_nitrate,                   &
                       tracer, tr_levels, tr_vars,                       &
                       tr_lbc_vars,tracer_lbc,                           &
                       a_max_trvars,a_tr_active_lbc_index,               &
                       tracer_ukca, tr_ukca,                             &
                       l_use_cariolle, ozone_tracer,                     &
                       tr_lbc_ukca,tracer_ukca_lbc,                      &
                       a_max_ukcavars,ukca_tr_active_lbc_index,          &
                       rimwidtha(rima_type_norm),rimweightsa,            &
                       lenrima(fld_type_p,halo_type_extended,            &
                               rima_type_norm),                          &
                       lbc_sizea(1,fld_type_p,halo_type_extended,        &
                                 rima_type_norm),                        &
                       lbc_starta(1,fld_type_p,halo_type_extended,       &
                                  rima_type_norm),                       &
                       errorstatus)
     
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AA SL_Tracer',6)

    END IF  L_tracer_1_if



! DEPENDS ON: timer
  IF (Ltimer) CALL TIMER('AA SL_Thermo',5)

  CALL eg_sl_thermo(                                                    &
                row_length, rows, n_rows, model_levels,                 &
                halo_i, halo_j, offx, offy,datastart, g_i_pe,           &
                l_inc_solver,high_order_scheme(theta_sl),               &
                monotone_scheme(theta_sl),l_high(theta_sl),             &
                l_mono(theta_sl), alpha_theta,                          &
                r_theta,etadot_np1,                                     &
                r_theta_d, errorstatus )

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Thermo',6)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Rho',5)

  CALL eg_sl_rho(                                                       &
                row_length, rows, n_rows, model_levels, halo_i,         &
                halo_j, offx, offy, datastart, g_i_pe,                  &
                l_inc_solver, depart_scheme, depart_order,              &
                high_order_scheme(rho_sl), monotone_scheme(rho_sl),     &
                ritchie_high_order_scheme, ritchie_monotone_scheme,     &
                first_constant_r_rho_level,interp_vertical_search_tol,  &
                check_bottom_levels, l_cartesian, l_shallow,            &
                l_rk_dps,                                               &
                l_high(rho_sl),l_mono(rho_sl), l_ritchie_high,          &
                l_ritchie_mono, lam_max_cfl, alpha_rho,                 &
                etadot, u_np1, v_np1, w_np1,etadot_np1, u_adv, v_adv,   &
                w_adv, r_rho, r_rho_d, r_p_p, r_p_p_d, errorstatus,     &
                hm_rhox, hm_rhoy )

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Rho',6)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AA SL_Moisture',5)

  IF (.NOT. l_dry) THEN
    CALL eg_sl_moisture( moisture_array_size,                           &
                row_length, rows, n_rows, model_levels, halo_i,         &
                halo_j, offx, offy, datastart, g_i_pe,                  &
                high_order_scheme(moist_sl),monotone_scheme(moist_sl),  &
                l_high(moist_sl),l_mono(moist_sl),l_pc2, L_mcr_qrain,   &
                l_mcr_qcf2, l_mcr_qgraup,                               &
                r_m_v, r_m_cl, r_m_cf,r_m_r, r_m_gr, r_m_cf2,           &
                cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
                exner_star,r_m_v_d, r_m_cl_d, r_m_cf_d,                 &
                r_m_r_d, r_m_gr_d, r_m_cf2_d,                           &
                cf_star, cfl_star, cff_star,                            &
                errorstatus)
  END IF

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('AA SL_Moisture',6)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS S-L Advect (AA)',6)

! ----------------------------------------------------------------------
! Section 3.0  Call Atmospheric Physics2
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',5)


IF (.NOT. l_run_with_physics2) THEN
  l_physics_store=l_physics
  l_physics=.FALSE.
END IF

IF (l_tracer .AND. cycleno == numcycles ) THEN
!
!  super_tracer_phys2 = tracers ( before physic2 )
!   
      CALL tr_set_phys_4A(super_array_size, super_tracer_phys2,   &
                   l_co2_interactive, co2, l_murk_advect, murk,   &
                   l_soot, soot_new, soot_agd,soot_cld,           &
                   l_sulpc_so2, so2,so4_aitken,so4_accu,so4_diss, &
                   l_sulpc_nh3, nh3,l_sulpc_dms, dms,             &
                   l_dust, dust_div1,dust_div2,dust_div3,         &
                           dust_div4,dust_div5,dust_div6,         &
                   l_biomass, bmass_new, bmass_agd, bmass_cld,    &
                   l_ocff, ocff_new, ocff_agd, ocff_cld,          &
                   l_nitrate, nitr_acc, nitr_diss,                &
                   l_use_cariolle, ozone_tracer,                  &
                   tracer_phys1, tracer,                          &
                   ukca_tracer_phys1, tracer_ukca,                &
                   row_length, rows,                              &
                   model_levels, tr_levels, tr_vars, tr_ukca,     &
                   offx, offy, model_domain,                      &
                   .FALSE., tdims )

END IF  ! L_tracer and CycleNo == NumCycles


!        F1SP computes the predictors for u,v,w aka R_u, R_v, theta_star
!        and moisture fields

IF( l_inc_solver) THEN
  CALL eg_f1sp_inc(theta_star ,m_star ,mcl_star,mcf_star              &
                    ,mcf2_star,mrain_star,mgraup_star                 &
                    ,u_np1,v_np1,w_np1,u,v,w, exner_np1               &
                    ,exner_surf_np1,rho_np1, thetav_np1, etadot_np1   &
                    ,m_v_np1, m_cl_np1, m_cf_np1, m_r_np1             &
                    ,m_gr_np1, m_cf2_np1,   ih, l_slice               &
                    ,row_length,rows,n_rows,model_levels              &
                    ,L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,             &
                     psi_w_surf,psi_w_lid)

ELSE
  CALL eg_f1sp(theta_star ,m_star ,mcl_star,mcf_star                  &
                     ,mcf2_star,mrain_star,mgraup_star                &
                     ,u_np1,v_np1,w_np1,u,v,w, exner_np1              &
                     ,exner_surf_np1,rho_np1, thetav_np1, etadot_np1  &
                     ,m_v_np1, m_cl_np1, m_cf_np1, m_r_np1            &
                     ,m_gr_np1, m_cf2_np1,   ih,l_slice               &
                     ,row_length,rows,n_rows,model_levels             &
                     ,L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,            &
                      psi_w_surf,psi_w_lid)
END IF

IF (l_mr_physics2) THEN
  q_star  (:,:,:) = m_star  (:,:,:)
  qcl_star(:,:,:) = mcl_star(:,:,:)
  qcf_star(:,:,:) = mcf_star(:,:,:)

  IF (l_mcr_qcf2  )    qcf2_star  (:,:,:) = mcf2_star  (:,:,:)
  IF (l_mcr_qrain )    qrain_star (:,:,:) = mrain_star (:,:,:)
  IF (l_mcr_qgraup)    qgraup_star(:,:,:) = mgraup_star(:,:,:)

ELSE
  CALL eg_mix_to_q                                                    &
                  (tdims_s,tdims_s,                                   &
                   m_star, mcl_star, mcf_star,                        &
                   mcf2_star, mrain_star, mgraup_star,                &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   q_star, qcl_star, qcf_star,                        &
                   qcf2_star,qrain_star,qgraup_star                   &
                   )
END IF

! ---------------------------------------------------------------
! Section 2.3 Diagnostics at end of advection
! ----------------------------------------------------------------------
! Apply diagnostics at final cycle only
IF ( cycleno == numcycles ) THEN

! Section 12: 'dynamics advection' based quantities
  IF( SF(0,12) .AND. ErrorStatus == 0 ) THEN

! Allocate diagnostic space for STASH
    ALLOCATE (stashwork12(stash_maxlen(12,a_im)))

    CALL diagnostics_adv(                                             &
                  row_length, rows, n_rows, model_levels, wet_levels, &
! primary wind fields:
                  u, v,                                               &
                  theta, q, qcl, qcf, qrain, qgraup, qcf2,            &
                  cf_bulk, cf_liquid, cf_frozen,                      &
! wind field increments after advection:
                  r_u_p2, r_v_p2, r_w_p2,                             &
! wind field increments before advection (on stashflag):
                  u_incr_diagnostic, v_incr_diagnostic,               &
                  t_incr_diagnostic, q_incr_diagnostic,               &
                  qcl_incr_diagnostic, qcf_incr_diagnostic,           &
                  qrain_incr_diagnostic, qgraup_incr_diagnostic,      &
                  qcf2_incr_diagnostic,                               &
                  cf_incr_diagnostic, cfl_incr_diagnostic,            &
                  cff_incr_diagnostic,                                &
                  theta_star, q_star, qcl_star, qcf_star,             &
                  qrain_star, qgraup_star, qcf2_star,                 &
                  cf_star, cfl_star, cff_star,                        &
                  exner_theta_levels,                                 &
! departure points for w
                  depart_xi1_rho, depart_xi2_rho, depart_xi3_rho,     &
                  r_theta_levels,                                     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                  stashwork12)

! DEPENDS ON: Atm_Step_diag
    CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
             r_v_p2, r_u_p2, r_w_p2, q_star, qcl_star, qcf_star,      &
             theta_star, cf_star, cfl_star, cff_star, 3)
        
! DEPENDS ON: stash
    CALL stash(a_sm,a_im,12,stashwork12,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          errorstatus,cmessage)

    DEALLOCATE (stashwork12)

  END IF !   SF(0,12)

END IF ! CycleNo == NumCycles


!        Save star fields to obtain increments after call to
!        Atmos_Physics2

!        we do not need to save the moisture, because we have those in
!        the ENDGame mixing ration fields anyway!

!        This section can very likely be replace by some of the logic in
!        atm_step_local!

  r_u_p2_n = r_u_p2
  r_v_p2_n = r_v_p2
  r_w_p2_n = r_w_p2

  theta_star_n = theta_star

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',6)

lphysics :  IF (l_physics) THEN
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',5)


! DEPENDS ON: Atm_Step_diag
          CALL Atm_Step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
             r_v_p2, r_u_p2, r_w_p2,                  &
             q_star, qcl_star, qcf_star,  theta_star, &
             cf_star, cfl_star, cff_star, 4)


! NB if you are changing the argument list to atmos_physics2, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent. Note there are two calls to atmos_physics2 in scm_main.

! DEPENDS ON: atmos_physics2
          Call atmos_physics2(                                           &
! Parallel variables
        global_row_length, global_rows, nproc, nproc_x, nproc_y          &
      , g_rows, g_row_length, NumCycles, CycleNo                         &
! field dimensions etc.
      ,     row_length, rows, n_rows, land_field, model_levels, nice     &
      ,     nice_use                                                     &
      ,     wet_levels, bl_levels, st_levels, sm_levels, cloud_levels    &
      ,     land_ice_points, soil_points, n_cca_lev, ntiles, tr_levels   &
      ,     first_constant_r_rho_level, DIM_CS1, DIM_CS2                 &
! IN Substepping information and switches
      ,     L_regular, l_mr_physics2                                     &
      ,     L_dry, L_lbc_old, LCAL360, Ltimer                            &
      ,     L_CLD_AREA                                                   &
! Model Parameters
      ,     RHcrit, CO2_MMR, tr_vars, tr_ukca                            &
! Vertical coordinate levels.
      ,     rho                                                          &
      ,     delta_lambda, delta_phi                                      &
      ,     gdlambda_p(lambda_start), dphi_p, wt_lambda_p, wt_lambda_u   & 
      ,     wt_phi_p, wt_phi_v, lat_rot_NP, long_rot_NP, f3_at_u         &
! Time stepping information
      , I_year, I_day_number, I_hour, I_minute, I_second,                &
! River routing
        AOCPL_ROW_LENGTH,AOCPL_P_ROWS,XPA,XUA,XVA,YPA,YUA,YVA,           &
        G_P_FIELD,G_R_FIELD,A_INTHD(16),lasize(1,fld_type_r,             &
        halo_type_no_halo),lasize(2,fld_type_r,halo_type_no_halo),       &
        glsize(1,fld_type_r),glsize(2,fld_type_r),RIVER_VEL,             &
        RIVER_MCOEF,I_RIVER_VN,RIV_DIRECTION,RIV_SEQUENCE,RIV_STORAGE,   &
!  Add inland basin outflow to arguments
        RIV_INLANDATM,                                                   &
!  Add lake evaporation: 
        ACC_LAKE_EVAP,                                                   &
! Grid-to-grid river routing
        RIV_IAREA,RIV_SLOPE, RIV_FLOWOBS1,RIV_INEXT, RIV_JNEXT,RIV_LAND, &
        RIV_SUBSTORE, RIV_SURFSTORE,RIV_FLOWIN,RIV_BFLOWIN,              &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
            STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,     &
            STASHwork26,                                                 &
! SCM Diagnostics (dummy values in full UM)
        nSCMDpkgs, L_SCMDiags,                                           &
! Data Fields.
       THETA, Q, QCL, QCF, qrain, qgraup, qcf2,                          &
       RHO_r_sq_n, U, V, W, W_ADV, P ,PSTAR,                             &
       EXNER,EXNER_THETA_LEVELS, LAND, P_THETA_LEVELS                    &

! ancillary fields and fields needed to be kept from timestep to
! timestep
      ,land_index, land_ice_index, soil_index, CANOPY_WATER              &
      ,SNODEP, THERM_COND, THERM_CAP, VOL_SMC_CRIT                       &
      ,VOL_SMC_WILT, VOL_SMC_SAT, STHF, STHU                             &
      ,OROG_SIL,OROG_HO2,OROG_SD,ICE_THICKNESS,ICE_FRACTION              &
      ,U_SEA,V_SEA,U_0_P,V_0_P,CCA,CCB                                   &
      ,CCT, CCLWP, CCW_RAD, LCBASE, DEEP_SOIL_TEMP, p_ti, TI, ICE_K_CAT  &
      ,TSTAR,Z0, p_ice_fract, p_ice_thick                                &
      ,SAT_SOIL_COND,SAT_SOILW_SUCTION,CLAPP_HORN                        &
      ,SMCL, T1_SD, Q1_SD, ZH, ddmfx, CF_AREA, CF_BULK, CF_LIQUID        &
      ,CF_FROZEN, ls_rain, ls_snow, micro_tends                          &
      ,photosynth_act_rad, rad_hr, SOIL_CLAY                             &
      ,SOIL_SILT,SOIL_SAND, DUST_MREL1,DUST_MREL2                        &
      ,DUST_MREL3,DUST_MREL4,DUST_MREL5,DUST_MREL6                       &
      ,SO2_HILEM, SO2_EM, NH3_EM, DMS_EM,SOOT_HILEM, SOOT_EM, OCFF_HILEM &
      ,OCFF_EM, CO2_EMITS, CO2FLUX ,deep_flag, past_precip, past_conv_ht &

! IN/OUT
      , theta_star,q_star,qcl_star,qcf_star, qrain_star, qgraup_star     &
      , qcf2_star, cf_star,cfl_star,cff_star                             &
      , r_u_p2, r_v_p2, r_w_p2                                           &
      , net_flux, net_mflux, murk,tracer, tracer_ukca                    &
      , DUST_DIV1, DUST_DIV2, DUST_DIV3 ,DUST_DIV4, DUST_DIV5, DUST_DIV6 &
      , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new          &
      , soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld, ocff_new    & 
      , ocff_agd, ocff_cld, nitr_acc, nitr_diss, co2                     &
! IN/OUT River routing
      , TOT_SURFROFF, TOT_SUBROFF                                        &
! OUT Fields
      , rho_km, cH, NTML, CUMULUS, NBDSC, NTDSC                          &
      , rhcpt, rhc_row_length, rhc_rows                                  &

! Additional variables for MOSES II
      , FRAC_TYP, DISTURB_VEG, CANHT_PFT, LAI_PFT                        &
      , CAN_WATER_TILE, CATCH_TILE, CATCH_SNOW                           &
      , SNOW_GRND, SNODEP_TILE, Z0_TILE, Z0H_TILE, TSTAR_TILE            &
      , INFIL_TILE, RGRAIN_TILE, CS, GS                                  &
      , co2_dim_row, co2_dim_len, A_INTHD(23)                            &
      , STEPim(atmos_im), G_LF_PFT_ACC, G_PHLF_PFT_ACC, NPP_PFT_ACC      &
      , RSP_W_PFT_ACC, RSA                                               &
      , land_pts_trif, npft_trif, dolr, LW_down, SW_tile                 &
      , FRAC_LAND,TSTAR_LAND,TSTAR_SEA,p_tstar_sice,TSTAR_SICE           &
      , SOIL_ALB,cos_zenith_angle                                        &
! INOUT variables for TKE based turbulence schemes 
      , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                     &
! Additional variables required for large-scale hydrology:
      , FEXP,GAMMA_INT,TI_MEAN,TI_SIG, FSFC_SAT,F_WETLAND                &
      , WATER_TABLE,STHZW, A_FSAT,C_FSAT,A_FWET,C_FWET                   &
! IN/OUT JULES 2 prognostics
      , SNOWDEPTH, RHO_SNOW_GRND, NSNOW                                  &
      , DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL                    &
! IN/OUT FLake lake scheme prognostics
      , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                  &
      , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                  &
      , lake_g_dt                                                        &
! Cariolle ozone 
      , OZONE_TRACER                                                     &

! Additional screen-level variables
      , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                           &

! Variables required for COSP (out)
      ,cosp_crain_3d, cosp_csnow_3d                                      &

! SCM and idealised UM surface forcing parameters
      ,L_flux_bc,flux_e,flux_h,L_spec_z0,z0m_scm,z0h_scm,                &

! error information
                            ErrorStatus )
     
                                    
IF (l_tracer .AND. cycleno == numcycles ) THEN
!
! set super_tracer_phys2 such that it contains only the atmos_physics2 contribution.
!     super_tracer_phys2 = tracers - super_tracer_phys2(before atmos_physics2)
!  
  CALL tr_reset_4A (super_array_size, super_tracer_phys2,    &
       l_co2_interactive, co2,  l_murk_advect, murk,         &
       l_soot, soot_new, soot_agd, soot_cld,                 &
       l_sulpc_so2, so2, so4_aitken, so4_accu, so4_diss,     &
       l_sulpc_nh3, nh3, l_sulpc_dms, dms, l_dust,           &
       dust_div1, dust_div2, dust_div3, dust_div4,           &
       dust_div5, dust_div6,                                 &
       l_biomass, bmass_new, bmass_agd, bmass_cld,           &
       l_ocff, ocff_new, ocff_agd, ocff_cld,                 &
       l_nitrate, nitr_acc, nitr_diss,                       &
       l_use_cariolle, ozone_tracer, tracer_phys2, tracer,   &
       ukca_tracer_phys2, tracer_ukca, row_length, rows,     &
       model_levels, tr_levels, tr_vars, tr_ukca, offx, offy ) 
END IF  ! L_tracer and CycleNo == NumCycles
                                                              

! set zero level star quantities for Endgame
! currently set equal to level 1
    CALL set_star_zero_level(                                           &
             theta_star,                                                &
             q_star,                                                    &
             qcl_star,                                                  &
             qcf_star,                                                  &
             cf_star,                                                   &
             cfl_star,                                                  &
             cff_star,                                                  &
             qcf2_star,                                                 &
             qrain_star,                                                &
             qgraup_star,                                               &
             L_mcr_qgraup,                                              &
             L_mcr_qrain,                                               &
             L_mcr_qcf2)


! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',6)

! ---------------------------------------------------------------

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)
! DEPENDS ON: Atm_Step_stash
    CALL atm_step_stash( &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         errorstatus, 2)
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',6)

    IF (.NOT. l_run_with_physics2 ) THEN
      l_physics=l_physics_store
    END IF

!  --------- UM Section 35---- Stochastic Physics SKEB2 -------------

! Call SKEB2 after the physics have been completed

 IF (l_skeb2) THEN
! SKEB2 only called once in the outer loop
   IF ( cycleno == 1 ) THEN

     IF (sf(0,35)) THEN
       ! ALLOCATE diagnostic space for STASH
       ALLOCATE (stashwork35(stash_maxlen(35,a_im)))
     ELSE
       ! Code to ALLOCATE unity arrays when not used
       ALLOCATE(stashwork35(1))
       stashwork35(1) = 0.0
     END IF

! DEPENDS ON: timer
     IF (ltimer) CALL timer ('SKEB2',3)

     CALL stph_skeb2(                                                   &
! in
           row_length, rows, n_rows, model_levels,                      &
           delta_phi, delta_lambda,                                     &
           rho_r_sq_n, u, v,                                            &
! in/out
           r_u_skeb, r_v_skeb,                                          &
! STASH array sf
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
           stashwork35, first_atmstep_call)

! DEPENDS ON: timer
     IF (ltimer) CALL timer ('SKEB2',4)

! Send diagnostics to STASH
     IF (sf(0,35) .AND. errorstatus == 0) THEN

! DEPENDS ON: stash
       CALL stash( a_sm, a_im, 35, stashwork35,                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                 errorstatus, cmessage)

     END IF !(sf(0,35))

     DEALLOCATE (stashwork35)

   END IF  !(cycleno)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Stochastic_Phys',6)
 END IF  !(l_skeb2)

ELSE

  CALL eg_idl_forcing2(                                               &
  u, v,  exner_theta_levels, exner,                                   &
  r_u_p2, r_v_p2, L_HeldSuarez1_drag,L_HeldSuarez2_drag, errorstatus )


END IF lphysics

! DEPENDS ON: timer
IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',5)

! Compute fast physics source terms for ENDGame
CALL eg_r_s(                                                            &
         theta_star,q_star,qcl_star,qcf_star,qcf2_star,                 &
         qrain_star, qgraup_star,                                       &
         theta_star_n,m_star,mcl_star,mcf_star,mcf2_star,               &
         mrain_star,mgraup_star, couple_app, l_skeb2)


! DEPENDS ON: timer
IF (Ltimer) CALL TIMER('AS Atmos_Phys2 (AP2)',6)

! DEPENDS ON: timer
IF (Ltimer) CALL TIMER('AS Diffusion',5)

IF (errorstatus  ==  0 ) THEN  !it2

! -------- DIFFUSION SECTION ------------------

  IF ( L_diff_active ) THEN !  it3

! DEPENDS ON: Atm_Step_diag
    CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star, &
         cf_star, cfl_star, cff_star, 5)

! DEPENDS ON: eg_diff_ctl
    CALL eg_diff_ctl(                                                   &
                     L_backwards,                                       &
                     timestep, pos_timestep, neg_timestep,              &
                     thetav_np1,m_v_np1,u_np1,v_np1,w_np1,rho_r_sq_n,   &
                     exner_rho_levels, exner_theta_levels,              &
                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,      &
                     eta_theta_levels, eta_rho_levels,                  &
                     offx, offy, halo_i, halo_j,                        &
                     row_length, rows, n_rows, model_levels,            &
                     first_constant_r_rho_level,                        &
                     L_tardiff_q, w_conv_limit, tardiffq_factor,        &
                     tardiffq_test, tardiffq_start, tardiffq_end,       &
                     s_thetav, s_m_v, s_u, s_v, s_w )

! ----------------------------------------------------------------------
! Section 13.1 Diagnostics at from diffusion and divergence damping
! ----------------------------------------------------------------------
! Apply diagnostics only at last cycle
  IF ( cycleno == numcycles ) THEN

! section 13:
    IF( sf(0,13) .AND. errorstatus == 0) THEN

! Allocate diagnostic space for STASH
      IF(.NOT. l_filter) ALLOCATE (stashwork13(stash_maxlen(13,a_im)))

! DEPENDS ON: diagnostics_dif
      CALL diagnostics_dif(                                        &
         row_length, rows, n_rows, model_levels, bl_levels,        &
! primary  fields:
         theta, q,                                                 &
! wind field increments after  dif :
         r_u, r_v,                                                 &
! wind field increments before diffusion (on stashflag):
         u_incr_diagnostic, v_incr_diagnostic,                     &
 t_incr_diagnostic,q_incr_diagnostic,                              &
         w_local_mask,                                             &
         theta_star, q_star,                                       &
         exner_theta_levels,                                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         stashwork13)

! DEPENDS ON: Atm_Step_diag
      CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         r_v, r_u, r_w, q_star, qcl_star, qcf_star, theta_star, &
         cf_star, cfl_star, cff_star, 6)

! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)
! DEPENDS ON: Atm_Step_stash
      CALL atm_step_stash( &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         errorstatus, 3)
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',6)

    END IF !   SF(0,13)

    END IF ! CycleNo == NumCycles

  END IF      ! L_diff_active  ei3

END IF       ! ErrorStatus  ==  0  ei2

       IF ( cycleno == numcycles ) THEN

! DEPENDS ON: Atm_Step_alloc_4A
        CALL Atm_Step_alloc_4A(                                         &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
             cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,   &
             frac_control, r_u, r_v, r_w, errorstatus, 'dealsmag' ) 

       END IF ! cycleno == numcycles

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Diffusion',6)

!---------- END DIFFUSION SECTION -------------------------------

!=== Polar filter + diffusion of fast physics section ==================================
! ----------------------------------------------------------------------
! Section 0.4  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------
      If ( L_polar_filter_incs .or. L_filter_incs ) then  !it4
! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Filter',5)
      IF (printstatus >= prstatus_diag) THEN
        max_r_theta     = MAXVAL(s_thetav(tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end,    &
                                          tdims%k_start:tdims%k_end))
        min_r_theta     = MINVAL(s_thetav(tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end,    &
                                          tdims%k_start:tdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_theta)
        CALL gc_rmin(1,n_proc,ierr,min_r_theta)

        max_r_u     = MAXVAL(s_u(udims%i_start:udims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        min_r_u     = MINVAL(s_u(udims%i_start:tdims%i_end,             &
                                 udims%j_start:udims%j_end,             &
                                 udims%k_start:udims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_u)
        CALL gc_rmin(1,n_proc,ierr,min_r_u)

        max_r_v     = MAXVAL(s_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        min_r_v     = MINVAL(s_v(vdims%i_start:vdims%i_end,             &
                                 vdims%j_start:vdims%j_end,             &
                                 vdims%k_start:vdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_v)
        CALL gc_rmin(1,n_proc,ierr,min_r_v)

        max_r_w     = MAXVAL(s_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        min_r_w     = MINVAL(s_w(wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end))
        CALL gc_rmax(1,n_proc,ierr,max_r_w)
        CALL gc_rmin(1,n_proc,ierr,min_r_w)

        IF (mype == 0) THEN
          WRITE(6,fmt='(A)') '==================================================='
          WRITE(6,fmt='(A)') 'Calling eg_NI_filter_incs_Ctl'
          WRITE(6,fmt='(A)') ' Max/Min before polar filter :'
          WRITE(6,fmt='(A,2E25.10)') '    S_u = ',MAX_r_u,MIN_r_u
          WRITE(6,fmt='(A,2E25.10)') '    S_v = ',MAX_r_v,MIN_r_v
          WRITE(6,fmt='(A,2E25.10)') '    S_w = ',MAX_r_w,MIN_r_w
          WRITE(6,fmt='(A,2E25.10)') 'S_theta = ',MAX_r_theta,MIN_r_theta
        END IF
      END IF

      CALL eg_NI_filter_incs_Ctl(                                       &
                            r_theta, r_u, r_v, r_w, r_rho, r_p_p,       &
                            s_thetav, s_u, s_v, s_w,                    &
                            thetav, u, v, w, rho, exner,                &
                            row_length, rows, n_rows, model_levels,     &
                            r_theta_levels, r_rho_levels,               &
                            r_at_u, r_at_v,                             &
                            max_121_rows, u_sweeps, v_sweeps,           &
                            global_u_filter, global_v_filter,           &
                            u_begin, u_end, v_begin, v_end,             &
                            diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
                            first_constant_r_rho_level,                 &
                            first_constant_r_rho_level_m1,              &
                            horizontal_level,                           &
                            offx, offy, halo_i, halo_j,                 &
                            nproc_y, at_extremity, model_domain,        &
                            L_polar_filter_incs, L_diff_incs,           &
                            L_pofil_hadgem2, Ltimer,                    &
                            l_eliminate_rho, cycleno,                   &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                            xi1_u, xi1_p, xi2_p,                        &
                            pole_consts, gc_proc_row_group,             &
                            global_row_length,                          &
                            csxi2_v, csxi2_p)

      IF (integrity_test) THEN
        CALL eg_swap_bounds(s_u,udims_s,fld_type_u,.TRUE.) 
        CALL eg_swap_bounds(s_v,vdims_s,fld_type_v,.TRUE.) 
        CALL eg_swap_bounds(s_w,wdims_s,fld_type_p,.FALSE.) 
        CALL eg_swap_bounds(s_thetav,tdims_s,fld_type_p,.FALSE.) 


        CALL update_hash_m(    S_u,     SIZE(S_u),             's_u__', &
                               S_v,     SIZE(S_v),             's_v__', &
                               S_w,     SIZE(S_w),             's_w__', &
                               S_thetav,SIZE(S_thetav),        's_the')
      END IF

      IF (printstatus >= prstatus_diag) THEN
          max_r_theta     = MAXVAL(s_thetav)
          min_r_theta     = MINVAL(s_thetav)
          CALL gc_rmax(1,n_proc,ierr,max_r_theta)
          CALL gc_rmin(1,n_proc,ierr,min_r_theta)
          max_r_u     = MAXVAL(s_u)
          min_r_u     = MINVAL(s_u)
          CALL gc_rmax(1,n_proc,ierr,max_r_u)
          CALL gc_rmin(1,n_proc,ierr,min_r_u)
          max_r_v     = MAXVAL(s_v)
          min_r_v     = MINVAL(s_v)
          CALL gc_rmax(1,n_proc,ierr,max_r_v)
          CALL gc_rmin(1,n_proc,ierr,min_r_v)
          max_r_w     = MAXVAL(s_w)
          min_r_w     = MINVAL(s_w)
          CALL gc_rmax(1,n_proc,ierr,max_r_w)
          CALL gc_rmin(1,n_proc,ierr,min_r_w)

          IF (mype == 0) THEN
            WRITE(6,fmt='(A)') ' '
            WRITE(6,fmt='(A)') ' Max/Min after polar incs filter :'
            WRITE(6,fmt='(A,2E25.10)') '    S_u = ',MAX_r_u,MIN_r_u
            WRITE(6,fmt='(A,2E25.10)') '    S_v = ',MAX_r_v,MIN_r_v
            WRITE(6,fmt='(A,2E25.10)') '    S_w = ',MAX_r_w,MIN_r_w
            WRITE(6,fmt='(A,2E25.10)') 'S_theta = ',MAX_r_theta,MIN_r_theta
            WRITE(6,fmt='(A)') '==================================================='
          END IF
      END IF

! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Filter',6)
      END IF    ! L_polar_filter_incs .or. L_filter_incs  ei4
!=== End Polar filter + diffusion of increments section ===================


! Obtain increments from Atmos_Physics2 and diffusion
! DEPENDS ON: Atm_Step_alloc_4A
  CALL atm_step_alloc_4A( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
       cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
       frac_control, r_u, r_v, r_w, errorstatus, 'phy2diff' )
       
!------------------------------------------------------------------------
! Set hydrostatic balance in LBC regions
!------------------------------------------------------------------------
      IF (model_domain ==  mt_lam) THEN
! DEPENDS ON: timer
        IF (Ltimer) CALL timer('AS LAM_LBCS',5)       

        IF (L_fixed_lbcs) THEN
          L_update_lbcs = .FALSE.
        ELSE
          L_update_lbcs = .TRUE.
        END IF !  L_fixed_lbcs

        L_do_halos=.FALSE.
        L_do_boundaries=.TRUE.

        IF (RIM_STEPSA  ==  0) THEN
          increment_factor = 0.0
          L_balance = .FALSE.
        ELSE IF ( MOD(Timestep_Number-1,RIM_STEPSA) == 0) THEN
          L_balance = .TRUE.
          increment_factor = 1.0 / RIM_STEPSA
        ELSE
          L_balance = .FALSE.
          increment_factor = 1.0 /                                     &
     &               (RIM_STEPSA - MOD(Timestep_Number-1,RIM_STEPSA))
        END IF

        lbc_size=LENRIMA(fld_type_p,halo_type_extended, rima_type_norm)

        IF ( L_update_lbcs ) THEN
! Obtain Exner tendency to pass to the solver to apply lbc

! Apply additional balance to EXNER_LBC_TEND, RHO_LBC_TEND 
! and W_LBC_TEND only on timestep after up_bound has been called.

          IF (L_LBC_balance .and. L_balance) THEN
! Potential error in which element of RIMWIDTHA to pass in
          WRITE(6,*) 'RIMWIDTHA(1),RIMWIDTHA((rima_type_norm)) ',      &
                      RIMWIDTHA(1),RIMWIDTHA((rima_type_norm)),rima_type_norm
          WRITE(6,*) 'RIMWIDTHA ',RIMWIDTHA

            CALL eg_balance_lbc_values(                                &
           EXNER_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,               &
           Q_LBC_TEND, W_LBC_TEND, W_ADV_LBC_TEND,                     &
           u_lbc_tend, V_LBC_TEND,                                     &
           QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,                  &
           QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                            &
           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                      &
           R_RHO_LEVELS, R_THETA_LEVELS,                               &
           ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I,HALO_J,  &
           LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),      &
           LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),      &
           LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),      &
           LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm), &
           LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm), &
           LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm), &
           RIMWIDTHA(rima_type_norm), N_RIMS_TO_DO,RIMWEIGHTSA,        &
           AT_EXTREMITY,                                               &
           DELTA_PHI, DELTA_LAMBDA,                                    &
           BASE_PHI, BASE_LAMBDA,                                      &
           DATASTART,                                                  &  
           LAT_ROT_NP, GLOBAL_ROW_LENGTH, GLOBAL_ROWS)
          END IF  !  L_LBC_balance and. L_balance

          DO k = 1, MODEL_LEVELS
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k) = increment_factor *              &
     &                         ( EXNER_LBC_TEND(i,k) - EXNER_LBC(i,k) )
            END DO
          END DO

        END IF ! L_update_lbcs
        
! DEPENDS ON: timer
        IF (Ltimer) CALL timer('AS LAM_LBCS',6)

      END IF  ! model_domain
!------------------------------------------------------------------------
! End - Set hydrostatic balance in LBC regions
!------------------------------------------------------------------------

! ----------------------------------------------------------------------
! Form and solve Helmholtz equation and update variables.
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Solver',5)

  IF (printstatus > prstatus_normal) exner_0 = exner_np1

  IF( l_inc_solver ) THEN  
     CALL eg_sl_helmholtz_inc(                                          &
         exner_np1, exner_surf_np1, u_np1, v_np1, w_np1,                &
         etadot_np1, rho_np1, thetav_np1, exner_prime_term,             &
         m_v_np1, m_cl_np1, m_cf_np1,m_r_np1, m_gr_np1, m_cf2_np1,      &
         ih, innits, gcr_precon_option, gcr_use_residual_tol,           &
         GCR_Diagnostics, r_u_d, r_v_d, r_w_d,r_theta_d, r_rho_d,       &
         r_m_v_d, r_m_cl_d, r_m_cf_d,r_m_r_d, r_m_gr_d, r_m_cf2_d,      &
         solver_tolerance, tol_sc_fact,n_rows, row_length, rows,        &
         model_levels,s_u,s_v,s_w,s_thetav,s_m_v,s_m_cl,s_m_cf,s_m_cf2, &
         s_m_r,s_m_gr,psi_w_surf, psi_w_lid)
  ELSE
     CALL eg_sl_helmholtz (                                             &
         exner_np1, exner_surf_np1, u_np1, v_np1, w_np1,                &
         etadot_np1, rho_np1, thetav_np1, exner_prime_term,             &
         m_v_np1, m_cl_np1, m_cf_np1,m_r_np1, m_gr_np1, m_cf2_np1,      &
         ih,  innits, gcr_precon_option, gcr_use_residual_tol,          &
         GCR_Diagnostics, r_u_d, r_v_d, r_w_d,r_theta_d, r_rho_d,       &
         l_eliminate_rho, r_m_v_d, r_m_cl_d, r_m_cf_d,                  &
         r_m_r_d, r_m_gr_d, r_m_cf2_d,solver_tolerance, tol_sc_fact,    &
         n_rows,row_length, rows, model_levels,                         &
         s_u,s_v,s_w,s_thetav,s_m_v,s_m_cl,s_m_cf,s_m_cf2,              &
         s_m_r,s_m_gr,psi_w_surf, psi_w_lid)
  END IF



IF ( printstatus > prstatus_diag ) THEN

exner_res_tmp = 0.0

DO k=pdims_s%k_start,pdims_s%k_end+1
  exner_res = maxval(abs((exner_np1(:,:,k)-exner_0(:,:,k))/exner_0(:,:,k)))
  CALL gc_rmax(1,nproc,istat,exner_res)

  IF (mype == 0 .AND. printstatus > prstatus_normal)                    &
        WRITE(0,fmt='(A,E32.16,2I5)') 'exner residual:',exner_res,cycleno,k

  exner_res_tmp = MAX(exner_res_tmp,exner_res)

END DO

IF (mype == 0 .AND. printstatus > prstatus_normal)                    &
        WRITE(0,fmt='(A,E32.16,I5)') 'exner residual - final:',        &
        exner_res_tmp,cycleno

exner_res=exner_res_tmp

END IF


! DEPENDS ON: timer
      IF (Ltimer) CALL timer ('AS Solver',6)

        If (L_tracer .and. CycleNo == NumCycles) then  !it5
! DEPENDS ON: timer
          If (Ltimer) CALL timer('AS S-L Advect (AA)',5)

! DEPENDS ON: timer
         IF (Ltimer) CALL timer('AA SL_Tracer',5)

!      Cariolle scheme is called to calculate the tracer ozone. All the tracers
!      are calculated at the end of the timestep. The ozone tracer 
!      calculated here will be used in the radiation scheme on the next timestep.

!      Insert if statement here to check that tracer that goes into this call 
!      ozone. Ozone may not be the only tracer.
            IF (PrintStatus  >=  PrStatus_Normal .AND.                   &
                first_atmstep_call .AND. mype == 0 ) THEN

               WRITE(6,fmt='(A,L1)') 'Atm_Step_4A: L_USE_CARIOLLE = ',    &
                                    L_USE_CARIOLLE
            END IF

            If (L_USE_CARIOLLE) then
               
               If (PrintStatus  >=  PrStatus_Normal .AND.                  &
                first_atmstep_call .AND. mype == 0) then

                    WRITE(6,fmt='(A,L1)') 'Atm_Step_4A: Calling Cariolle_o3_psc'
               End If

! DEPENDS ON: cariolle_o3_psc
               CALL cariolle_o3_psc (OZONE_TRACER,                       &
                      O3_PROD_LOSS,   O3_P_L_VMR,                        &
                      O3_VMR,         O3_P_L_TEMP,                       &
                      O3_TEMP,        O3_P_L_COLO3,                      &
                      O3_COLO3,                                          &
                      THETA,                                             &
                      P_THETA_LEVELS,                                    &
                      offx,offy,theta_off_size,                          &
                      rows,row_length,                                   &
                      timestep,                                          &
                      EXNER_THETA_LEVELS,                                &
                      model_levels)      
!    
! Halos updated
               CALL eg_swap_bounds(OZONE_TRACER,                            &
                        pdims_s, fld_type_p, .FALSE.)
            
            else
               If (PrintStatus  >=  PrStatus_Normal .AND.                &
                   first_atmstep_call) then

                   WRITE(6,*) 'Atm_Step_4A: Cariolle scheme not called'
               End If
            End if

! DEPENDS ON: timer
         IF (Ltimer) CALL timer('AA SL_Tracer',6)

! DEPENDS ON: timer
         IF (Ltimer) CALL timer('AS S-L Advect (AA)',6)

        END IF  ! L_tracer ei5


! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Diffusion',5)

!=== Targetted diffusion of m_v_np1 ======================================
      IF ( l_tardiff_q ) THEN
        IF (printstatus >= prstatus_diag) THEN
          max_q_star     = MAXVAL(m_v_np1)
          min_q_star     = MINVAL(m_v_np1)
          CALL gc_rmax(1,n_proc,ierr,max_q_star)
          CALL gc_rmin(1,n_proc,ierr,min_q_star)
          IF ( mype == 0 ) THEN
            WRITE(6,fmt='(A)')                                           &
                     '==================================================='
            WRITE(6,fmt='(A)') 'Calling tardiff for m_v_np1'
            WRITE(6,fmt='(A)') ' Max/Min before polar filter :'
            WRITE(6,fmt='(A,2E25.10)') 'm_v_np1 = ',MAX_q_star,MIN_q_star
          END IF
        END IF  

! DEPENDS ON: eg_tardiff_q_w
         CALL eg_tardiff_q_w(                                            &
                           m_v_np1, w_np1,                               &
                           r_theta_levels, r_rho_levels,                 &
                           offx, offy, halo_i, halo_j,                   &
                           rows, n_rows, row_length,                     &
                           model_levels, model_domain,                   &
                           w_conv_limit, tardiffq_factor, tardiffq_test, &
                           tardiffq_end, .FALSE.,                        &
                           csxi2_p,csxi2_v, nproc_y, L_pofil_hadgem2,    &
                           u_sweeps, u_begin, u_end, max_121_rows,       &
                           horizontal_level,                             &
                           first_constant_r_rho_level_m1)

         CALL eg_swap_bounds(  m_v_np1, tdims_s, fld_type_p, .FALSE.)

         IF (integrity_test)                                             &
           CALL update_hash(m_v_np1,     SIZE(m_v_np1),        'mvnp1')

         IF (printstatus >= prstatus_diag) THEN         
           max_q_star     = MAXVAL(m_v_np1)
           min_q_star     = MINVAL(m_v_np1)
           CALL gc_rmax(1,n_proc,ierr,max_q_star)
           CALL gc_rmin(1,n_proc,ierr,min_q_star)
           IF ( mype == 0 ) THEN
             WRITE(6,fmt='(A)') ' '
             WRITE(6,fmt='(A)') ' Max/Min after polar incs filter :'
             WRITE(6,fmt='(A,2E25.10)') 'm_v_np1 = ',MAX_q_star,MIN_q_star
             WRITE(6,fmt='(A)') '===================================================' 
           END IF
         END IF
      END IF    ! l_tardiff_q
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS Diffusion',6)
!=== End Targetted diffusion of m_v_np1 =====================================

! If a real data run reset lowest levels to reflect ND behaviour
  IF( .NOT. L_idealised_data ) THEN
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        thetav_np1(i,j,0) = thetav_np1(i,j,1)
        m_v_np1(i,j,0)    = m_v_np1(i,j,1)
        m_cl_np1(i,j,0)   = m_cl_np1(i,j,1)
        m_cf_np1(i,j,0)   = m_cf_np1(i,j,1)
        m_r_np1(i,j,0)    = m_r_np1(i,j,1)
        m_gr_np1(i,j,0)   = m_gr_np1(i,j,1)
        m_cf2_np1(i,j,0)  = m_cf2_np1(i,j,1)
      END DO
    END DO
    IF (integrity_test)                                               &
      CALL update_hash_m(                                             &
                  m_v_np1,         SIZE(m_v_np1),         'mvnp1',    &
                  m_cl_np1,        SIZE(m_cl_np1),        'mclp1',    &
                  m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',    &
                  m_r_np1,         SIZE(m_r_np1),         'mrnp1',    &
                  m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',    &
                  m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',    &
                  thetav_np1,      SIZE(thetav_np1),      'tvnp1')
  END IF


! remove any potential instabillity:
  IF(l_eg_dry_static_adj) CALL eg_dry_static_adj( thetav_np1,      &
                                      rho_np1,exner_np1,intw_w2rho)

END DO  outerloop ! end iterations for trajectory calc (Outer loop)

!
! fix dry density to maintain mass conservation of dry air
!
IF (L_fix_mass .AND. model_domain == mt_global) THEN

  total_rho   = eg_total_mass(rho_np1)

  IF( l_use_old_mass_fix ) THEN
    mass_fix_factor=total_rho_init/total_rho
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          rho_np1(i,j,k)=rho_np1(i,j,k)*mass_fix_factor
        END DO
      END DO
    END DO

  ELSE

    total_rho_r = eg_total_mass_r(rho_np1) 
    total_gr    = eg_total_gr  (rho_np1) 
    total_gr_r  = eg_total_gr_r(rho_np1) 

    mass_fix_A  = (total_rho_init*total_gr_r-total_rho_r*total_gr)       &
                 /(total_rho*total_gr_r-total_rho_r*total_gr)

    mass_fix_B  = ((total_rho  -total_rho_init)*total_gr)                &
                 /(total_rho*total_gr_r-total_rho_r*total_gr)


    DO k = pdims%k_start, pdims%k_end
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end

           mass_fix_factor   = mass_fix_A + mass_fix_B*r_rho_levels(i,j,k)
           rho_np1(i,j,k)    = rho_np1(i,j,k)*mass_fix_factor

           mass_fix_factor   = mass_fix_A + mass_fix_B*r_theta_levels(i,j,k)

           thetav_np1(i,j,k) = thetav_np1(i,j,k)/mass_fix_factor

        END DO 
      END DO 
    END DO 

    k = 0 
    DO j = pdims_s%j_start, pdims_s%j_end 
      DO i = pdims_s%i_start, pdims_s%i_end 
        thetav_np1(i,j,0) = thetav_np1(i,j,1)
      END DO 
    END DO 

    mass_fix_factor=total_rho_init/total_rho
   
  END IF

  IF (PrintStatus > PrStatus_Oper .AND. mype ==0) THEN
     WRITE(6,fmt='(A,E16.12)') 'Mass Loss corrected for ', mass_fix_factor
  END IF

END IF


! compute the new wet rho * r**2 

DO k=1, model_levels
  DO j=pdims%j_start-offy, pdims%j_end+offy
    DO i=pdims%i_start-offx, pdims%i_end+offx

        mixing_r=(intw_w2rho(k,1)*(m_v_np1  (i,j,k)      +              &
                                   m_r_np1  (i,j,k)      +              &
                                   m_gr_np1 (i,j,k)      +              &
                                   m_cl_np1 (i,j,k)      +              &
                                   m_cf_np1 (i,j,k)      +              &
                                   m_cf2_np1(i,j,k))     +              &
                  intw_w2rho(k,2)*(m_v_np1  (i,j,k-1)    +              &
                                   m_r_np1  (i,j,k-1)    +              &
                                   m_gr_np1 (i,j,k-1)    +              &
                                   m_cl_np1 (i,j,k-1)    +              &
                                   m_cf_np1 (i,j,k-1)    +              &
                                   m_cf2_np1(i,j,k-1)))

!     convert to wet rho
      rho_r_sq_np1(i,j,k) = rho_np1(i,j,k)*(1.+mixing_r)    

!     multiply by r**2
      rho_r_sq_np1(i,j,k) = rho_r_sq_np1(i,j,k)*( r_rho_levels(i,j,k)   &
                                                 *r_rho_levels(i,j,k))
    END DO
  END DO
END DO

CALL eg_swap_bounds(rho_r_sq_np1,pdims_s,fld_type_p,.FALSE.)

IF ( sf(188,30) ) THEN
  inc_rho(:,:,:) = rho_r_sq_np1(:,:,:) - rho_r_sq_n
END IF


! 
! fix moisture mixing ratios to satisfy mass conservation if 
! L_conserv(moist_SL) = true; if not then the
! eg_correct_moisture routine will simply clip the fields that 
! are below some minimum values [ m >= m_min ]
!   

  IF (Ltimer) CALL timer('EG_CORRECT_QS',5) 
                                                  
  CALL eg_correct_moisture(rho, rho_np1,                          &
       r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,             &
       m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
       s_m_v, s_m_cl, s_m_cf, s_m_r, s_m_gr, s_m_cf2,             & 
       qlimit, L_mcr_qrain, L_mcr_qgraup, L_mcr_qcf2,             &
       L_conserv(moist_SL),  mype, L_conserv_smooth_lap           ) 
    
  IF (Ltimer) CALL timer('EG_CORRECT_QS',6)
     
  
IF ( L_tracer ) THEN 
! 
!  fix tracers to satisfy mass conservation if 
!  L_conserve_tracers = true; if not then the
!  eg_correct_tracers routine will simply clip the fields that 
!  are below some minimum values (t >= t_min (generally zero) )
!   
  
  IF (Ltimer) CALL timer('EG_CORRECT_TRACERS',5) 
                
  CALL eg_correct_tracers(                                        &
                          mype, super_array_size,                 &
                          super_tracer_phys1, super_tracer_phys2, &
                          rho, rho_np1,                           &
                          CO2, L_CO2_interactive,                 &
                          murk, L_murk_advect,                    &
                          soot_new, soot_agd, soot_cld, L_soot,   &
                          bmass_new, bmass_agd, bmass_cld,        &
                          L_biomass,                              &
                          ocff_new, ocff_agd, ocff_cld, l_ocff,   &
                          DUST_DIV1,DUST_DIV2,DUST_DIV3,          &
                          DUST_DIV4,DUST_DIV5,DUST_DIV6,          &
                          L_DUST,                                 &
                          so2, so4_aitken, so4_accu,              &
                          so4_diss, nh3, dms,                     &
                          L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,  &
                          nitr_acc, nitr_diss, L_nitrate,         &
                          L_USE_CARIOLLE, OZONE_TRACER,           &
                          tr_ukca, tracer_ukca,                   &
                          L_conserve_tracers,                     &
                          L_conserv_smooth_lap                    )
                          
  IF (Ltimer) CALL timer('EG_CORRECT_TRACERS',6)
                                                                
END IF  ! L_tracer
   



IF ( errorstatus == 0 ) THEN

! DEPENDS ON: Atm_Step_alloc_4A
 CALL atm_step_alloc_4A( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
      cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
      frac_control, r_u, r_v, r_w, errorstatus,'lbc_updt')

END IF        !  ErrorStatus  ==  0

thetav(:,:,0) = thetav_np1(:,:,0)
   w(:,:,0)   = w_np1(:,:,0)
DO k = 1, model_levels
   u(:,:,k)      = u_np1(:,:,k)
   v(:,:,k)      = v_np1(:,:,k)
   w(:,:,k)      = w_np1(:,:,k)
   rho(:,:,k)    = rho_np1(:,:,k)
   exner(:,:,k)  = exner_np1(:,:,k)
   etadot(:,:,k) = etadot_np1(:,:,k)
   thetav(:,:,k) = thetav_np1(:,:,k)
END DO
k = model_levels+1
exner(:,:,k)  = exner_np1(:,:,k)

DO  k = 0, model_levels
   m_v(:,:,k)    = m_v_np1(:,:,k)
   m_cl(:,:,k)   = m_cl_np1(:,:,k)
   m_cf(:,:,k)   = m_cf_np1(:,:,k)
   m_cf2(:,:,k)  = m_cf2_np1(:,:,k)
   m_r(:,:,k)    = m_r_np1(:,:,k)
   m_gr(:,:,k)   = m_gr_np1(:,:,k)
END DO
exner_surf(:,:)   = exner_surf_np1(:,:)

IF (L_pc2) THEN

  cf_bulk(1:row_length,1:rows,:)=cf_star(1:row_length,1:rows,:)
  cf_liquid(1:row_length,1:rows,:)=cfl_star(1:row_length,1:rows,:)
  cf_frozen(1:row_length,1:rows,:)=cff_star(1:row_length,1:rows,:)

END IF
! ---------------------------------------------------------------------
!   Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------
IF ( model_domain == mt_lam ) THEN
! DEPENDS ON: TIMER
  IF (ltimer) CALL timer('AS LAM_LBCS',5)

  IF (L_Fixed_lbcs) THEN
    L_update_lbcs = .FALSE.
    L_apply_lbcs = .TRUE.
  ELSE IF ( RIM_STEPSA == 0 ) THEN
    L_update_lbcs = .FALSE.
  ELSE
    L_update_lbcs = .TRUE.
  END IF !  L_Fixed_lbcs

  IF ( L_update_lbcs ) THEN
        
! DEPENDS ON: BOUNDVAL
    CALL BOUNDVAL(LENRIMA(1,1,rima_type_norm),                 &
                  L_mcr_qcf2_lbc, L_mcr_qrain_lbc,             &
                  L_mcr_qgraup_lbc, L_pc2_lbc,                 &
                  L_murk_lbc, L_int_uvw_lbc,                   &
                  L_dust_div1_lbc,L_dust_div2_lbc,             &
                  L_dust_div3_lbc,L_dust_div4_lbc,             &
                  L_dust_div5_lbc,L_dust_div6_lbc,             &
                  L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,        &
                  L_so4_accu_lbc,L_so4_diss_lbc,               &
                  L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,     &
                  L_soot_cld_lbc,L_bmass_new_lbc,              &
                  L_bmass_agd_lbc,L_bmass_cld_lbc,             &
                  L_ocff_new_lbc,                              &
                  L_ocff_agd_lbc,L_ocff_cld_lbc,               &
                  L_nitr_acc_lbc, L_nitr_diss_lbc,             &
                  U_LBC, U_LBC_TEND,                           &
                  V_LBC, V_LBC_TEND,                           &
                  W_LBC, W_LBC_TEND,                           &
                  RHO_LBC, RHO_LBC_TEND,                       &
                  THETA_LBC, THETA_LBC_TEND,                   &
                  Q_LBC, Q_LBC_TEND,                           &
                  QCL_LBC, QCL_LBC_TEND,                       &
                  QCF_LBC, QCF_LBC_TEND,                       &
                  QCF2_LBC, QCF2_LBC_TEND,                     &
                  QRAIN_LBC, QRAIN_LBC_TEND,                   &
                  QGRAUP_LBC, QGRAUP_LBC_TEND,                 &
                  CF_BULK_LBC, CF_BULK_LBC_TEND,               &
                  CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,           &
                  CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,           &
                  EXNER_LBC, EXNER_LBC_TEND,                   &
                  U_ADV_LBC, U_ADV_LBC_TEND,                   &
                  V_ADV_LBC, V_ADV_LBC_TEND,                   &
                  W_ADV_LBC, W_ADV_LBC_TEND,                   &
                  MURK_LBC, MURK_LBC_TEND,                     &
                  DUST_DIV1_LBC, DUST_DIV1_LBC_TEND,           &
                  DUST_DIV2_LBC, DUST_DIV2_LBC_TEND,           &
                  DUST_DIV3_LBC, DUST_DIV3_LBC_TEND,           &
                  DUST_DIV4_LBC, DUST_DIV4_LBC_TEND,           &
                  DUST_DIV5_LBC, DUST_DIV5_LBC_TEND,           &
                  DUST_DIV6_LBC, DUST_DIV6_LBC_TEND,           &
                  SO2_LBC, SO2_LBC_TEND,                       &
                  DMS_LBC, DMS_LBC_TEND,                       &
                  SO4_AITKEN_LBC, SO4_AITKEN_LBC_TEND,         &
                  SO4_ACCU_LBC, SO4_ACCU_LBC_TEND,             &
                  SO4_DISS_LBC, SO4_DISS_LBC_TEND,             &
                  NH3_LBC, NH3_LBC_TEND,                       &
                  SOOT_NEW_LBC, SOOT_NEW_LBC_TEND,             &
                  SOOT_AGD_LBC, SOOT_AGD_LBC_TEND,             &
                  SOOT_CLD_LBC, SOOT_CLD_LBC_TEND,             &
                  BMASS_NEW_LBC, BMASS_NEW_LBC_TEND,           &
                  BMASS_AGD_LBC, BMASS_AGD_LBC_TEND,           &
                  BMASS_CLD_LBC, BMASS_CLD_LBC_TEND,           &
                  OCFF_NEW_LBC, OCFF_NEW_LBC_TEND,             &
                  OCFF_AGD_LBC, OCFF_AGD_LBC_TEND,             &
                  OCFF_CLD_LBC, OCFF_CLD_LBC_TEND,             &
                  NITR_ACC_LBC, NITR_ACC_LBC_TEND,             &
                  NITR_DISS_LBC, NITR_DISS_LBC_TEND,           &
                  TRACER_LBC, TRACER_LBC_TEND,                 &
                  TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,       &
                  1, 0, ErrorStatus, CMESSAGE)

  END IF ! L_update_lbcs

!--------------------------------------------------------------
!           Update primary fields with LAM LBC data
!--------------------------------------------------------------
  IF (model_domain /= mt_global) THEN
!DEPENDS ON: update_lam_lbcs
  CALL update_lam_lbcs(                                                &
         r_rho_levels, r_theta_levels,                                 &
         ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
         TR_VARS,TR_LBC_VARS,TR_LEVELS,                                &
         A_MAX_TRVARS,A_tr_active_lbc_index,                           &
         TR_UKCA,TR_LBC_UKCA,                                          &
         A_MAX_UKCAVARS,UKCA_tr_active_lbc_index,                      &     
         OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
         L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
         L_murk, L_murk_lbc,                                           &
         L_LBC_balance, L_int_uvw_lbc,                                 &
         L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,   &
         L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,   &
         L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,   &
         L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,&
         L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,       &
         L_nh3, L_nh3_lbc,                                             &
         L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,       &
         L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,     &
         L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,   &
         L_ocff_new, L_ocff_new_lbc,                                   &
         L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,       &
         L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,     &
         RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
         LENRIMA(1,1,rima_type_norm),                                  &
         LBC_SIZEA(1,1,1,rima_type_norm),                              &
         LBC_STARTA(1,1,1,rima_type_norm),                             &
         THETA_LBC, Q_LBC, QCL_LBC,                                    &
         QCF_LBC, QCF2_LBC, QRAIN_LBC,                                 &
         QGRAUP_LBC, CF_BULK_LBC, CF_LIQUID_LBC,                       &
         CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
         U_LBC, V_LBC, W_LBC,                                          &
         U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,                              &
         MURK_LBC,                                                     &
         DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                  &
         DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                  &
         SO2_LBC, DMS_LBC, SO4_AITKEN_LBC,                             &
         SO4_ACCU_LBC,SO4_DISS_LBC,NH3_LBC,                            &
         SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,                     &
         BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                  &
         OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                     &
         NITR_ACC_LBC, NITR_DISS_LBC,                                  &
         TRACER_LBC, TRACER_UKCA_LBC,                                  &
         qdims_s, thetav, m_v, m_cl, m_cf,                             &
         m_cf2, m_r, m_gr,                                             &
         CF_BULK, CF_LIQUID, CF_FROZEN,                                &
         RHO, EXNER,                                                   &
         U, V, W, U_ADV, V_ADV, W_ADV, MURK,                           &
         DUST_DIV1, DUST_DIV2, DUST_DIV3,                              &
         DUST_DIV4, DUST_DIV5, DUST_DIV6,                              &
         SO2, DMS, SO4_AITKEN,SO4_ACCU,                                &
         SO4_DISS, NH3,                                                &
         SOOT_NEW, SOOT_AGD, SOOT_CLD,                                 &
         BMASS_NEW, BMASS_AGD, BMASS_CLD,                              &
         OCFF_NEW, OCFF_AGD, OCFF_CLD,                                 &
         NITR_ACC, NITR_DISS,                                          &     
         DELTA_PHI, DELTA_LAMBDA,                                      &
         BASE_PHI, BASE_LAMBDA,                                        &
         DATASTART, lat_rot_NP,                                        &
         GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                               &
         TRACER, TRACER_UKCA  ) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMP. FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------!
! Re-initialise etadot: fixes needed for lateral bcs                    !
!-----------------------------------------------------------------------!
  DO k = 1, model_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

          u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +   &
                                     intw_u2p(i,2)*u(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +     &
                                     intw_u2p(i,2)*u(i,j,k) )

          v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +   &
                                     intw_v2p(j,2)*v(i,j,k+1) ) +   &
                   intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +     &
                                     intw_v2p(j,2)*v(i,j,k) )

          etadot(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -                  &
                             u_at_w*dxi1_xi3(i,j,k)/                    &
                                           h1_p_eta(i,j,k) -            &
                             v_at_w*dxi2_xi3(i,j,k)/                    &
                                           h2_p_eta(i,j,k) ) /          &
                                             deta_xi3_theta(i,j,k)
      END DO
    END DO
  END DO

  etadot(:,:,0) = 0.0
  etadot(:,:,model_levels) = 0.0

  END IF  ! .NOT. GLOBAL
! DEPENDS ON: TIMER
  IF (ltimer) CALL timer('AS LAM_LBCS',6)
END IF     !   model_domain  ==  mt_lam
! ---------------------------------------------------------------------
!   End - Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------

CALL eg_set_adv_winds(u,v,etadot,                                      &
                      u_adv,v_adv,w_adv,row_length,rows,n_rows,        &
                      model_levels, halo_i, halo_j, l_shallow)

idx = 1
DO k = pdims_s%k_start,pdims_s%k_end + 1
  DO j = pdims_s%j_start,pdims_s%j_end
    DO i = pdims_s%i_start,pdims_s%i_end

      exner_rho_levels(idx) = exner(i,j,k)
      idx = idx + 1

    END DO
  END DO
END DO

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  etadot,          SIZE(etadot),          'ed___',    &
                  exner,           SIZE(exner),           'pi___',    & 
                  u,               SIZE(u),               'u____',    &
                  v,               SIZE(v),               'v____',    &
                  w,               SIZE(w),               'w____',    & 
                  u_adv,           SIZE(u_adv),           'u_adv',    &
                  v_adv,           SIZE(v_adv),           'v_adv',    &
                  w_adv,           SIZE(w_adv),           'w_adv')


IF (printstatus >= prstatus_diag) THEN
  ! Obtain domain maximum wind and total mass diagnostics
  w_max = MAXVAL(ABS(w))
  CALL gc_rmax(1,nproc,istat,w_max)
END IF

IF (L_fix_mass) THEN
  IF ( mype == 0 .AND. printstatus >= prstatus_diag ) THEN
    OPEN(UNIT=eg_unit,FILE='wmax.dat',POSITION='append')
    WRITE(UNIT=eg_unit,fmt='(I7,2E25.16)') timestep_number,   &
                                           w_max, mass_fix_factor
    CLOSE(UNIT=eg_unit)
  END IF
ELSE

  total_rho = eg_total_mass(rho)
  IF ( mype == 0 .AND. printstatus >= prstatus_diag) THEN
    OPEN(UNIT=eg_unit,FILE='wmax.dat',POSITION='append')
    WRITE(UNIT=eg_unit,fmt='(I7,2E25.16)') timestep_number,   &
                                           w_max, total_rho
    CLOSE(UNIT=eg_unit)
  END IF
END IF

! DEPENDS ON: exitchek
CALL exitchek(atmos_im, lexitnow)

IF (l_mr_pc2 .OR. l_mr_qtbalcld) THEN

  q  (1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_v (:,:,:)
  qcl(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_cl(:,:,:)
  qcf(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)  =       &
                m_cf(:,:,:)

  IF (l_mcr_qcf2  )                                                     &
     qcf2  (1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_cf2(:,:,0:model_levels)

  IF (l_mcr_qrain )                                                     &
     qrain (1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_r(:,:,0:model_levels)

  IF (l_mcr_qgraup)                                                     &
     qgraup(1-offx:row_length+offx, 1-offy:rows+offy,:) =               &
     m_gr(:,:,0:model_levels)

ELSE

  CALL eg_mix_to_q                                                      &
                  (tdims_l,tdims_s,                                     &
                   m_v, m_cl, m_cf,                                     &
                   m_cf2, m_r, m_gr,                                    &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup)

END IF

CALL eg_thetav_theta                                                    &
                  (thetav, theta, m_v,                                  &
                   p, pstar, p_theta_levels,                            &
                   exner, exner_surf, exner_theta_levels)


! ------------------------------------------------------------------
! Section 17  Aerosol Modelling - includes Sulphur Cycle, soot, biomass
! and fossil-fuels organic carbon (OCFF) 
! ------------------------------------------------------------------
!
! DEPENDS ON: timer
      IF (ltimer) CALL timer('AS Aerosol Modelling',5)
!
      aeroif: IF (l_sulpc_so2.OR.l_soot.OR.l_biomass.OR.                &
                  l_ocff.OR.l_nitrate.OR.l_dust) THEN
!
! Allocate diagnostic space for STASH
        ALLOCATE (stashwork17(stash_maxlen(17,a_im)))
!
! Don't call Swap_bounds for fields used in Aero_Ctl_4A
!
        ALLOCATE(o3_mmr  (tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tdims%k_start:tdims%k_end))
        ALLOCATE(hno3_mmr(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tdims%k_start:tdims%k_end))
        ALLOCATE(h2o2_mmr(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tdims%k_start:tdims%k_end))
        ALLOCATE(oh_conc (tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tdims%k_start:tdims%k_end))
        ALLOCATE(ho2_conc(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tdims%k_start:tdims%k_end)) 
!

        CALL get_sulpc_oxidants(                                        &
       l_sulpc_online_oxidants, l_ukca, l_ukca_trop, l_ukca_tropisop,   &
       l_ukca_strattrop, l_ukca_raq, first_atmstep_call,                &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
       o3_mmr(:,:,1:), hno3_mmr(:,:,1:), h2o2_mmr(:,:,1:),              &
       oh_conc(:,:,1:), ho2_conc(:,:,1:))

        CALL aero_ctl_4A(                                                &
! Parallel variables
        halo_i, halo_j, offx, offy, global_row_length, global_rows ,     &
        gc_proc_row_group, gc_proc_col_group,                            &
        at_extremity, nproc, nproc_x, nproc_y  ,                         &
        neighbour, g_rows, g_row_length,  mype ,                         &
! model dimensions
        row_length, rows, n_rows, land_field    ,                        &
        model_levels, wet_levels, bl_levels, n_cca_lev  ,                &
        theta_field_size           ,                                     &
        salt_dim1, salt_dim2, salt_dim3 ,                                &
        aero_dim1, aero_dim2, aero_dim3,                                 &
! model switches
        model_domain, lcal360, l_sec_var, l_eqt, ltimer ,                &
! model parameters
        ntot_land, ntot_sea ,                                            &
! co-ordinate information
        delta_lambda, delta_phi ,                                        &
        lat_rot_np, long_rot_np,                                         &
! time stepping information
        timestep      ,                                                  &
        i_year, i_day_number, i_hour, i_minute ,                         &
        i_second, timestep_number ,                                      &
        previous_time ,                                                  &
        call_chem_freq  ,                                                &
! trig arrays
        sin_theta_longitude, cos_theta_longitude  ,                      &
        fv_cos_theta_latitude ,                                          &
! grid-dependent arrays
            f3_at_u, true_longitude, true_latitude,                      &
!
! data fields in
        u, v, tstar, tstar_sea,                                          &
        theta(:,:,1:), q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:),              &
        rho_r_sq_np1, land, frac_land,                                   &
        p_theta_levels(:,:,1:),                                          &
        exner_rho_levels,                                                &
        exner_theta_levels(:,:,1:),                                      &
        ice_fraction, snodep  ,                                          &
        cf_bulk(:,:,1:),                                                 &
        oh_conc(:,:,1:), h2o2_mmr(:,:,1:), ho2_conc(:,:,1:),             &
        o3_mmr(:,:,1:), hno3_mmr(:,:,1:),                                &
        so2_em, so2_hilem, so2_natem(:,:,1:),                            &
        dms_em, dms_conc, nh3_em ,                                       &
        dms_oflux ,                                                      &
        soot_em, soot_hilem, bmass_em, bmass_hilem, ocff_em, ocff_hilem, &
        land_index,                                                      &
! logicals in
        l_sulpc_so2, l_sulpc_dms, l_sulpc_ozone,                         &
        l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3,                         &
        l_sulpc_online_oxidants,                                         &
        l_sulpc_2_way_coupling,                                          &
        l_use_sulphate_sulpc, l_use_seasalt_sulpc, l_soot,               &
        l_use_soot_sulpc, l_biomass, l_use_bmass_sulpc,                  &
        l_ocff, l_use_ocff_sulpc,                                        &
        l_nitrate, l_use_nitrate_sulpc,                                  &
        l_so2_surfem, l_so2_hilem, l_so2_natem, l_dms_em,                &
        l_dms_em_inter, l_dms_liss_merlivat,                             &
        l_dms_wanninkhof, l_dms_nightingale,                             &
        l_dms_ointer,                                                    &
        l_nh3_em, l_ctile,                                               &
        l_soot_surem, l_soot_hilem, l_bmass_surem, l_bmass_hilem,        &
        l_ocff_surem, l_ocff_hilem, l_use_biogenic,                      &
        l_use_seasalt_direct, l_use_seasalt_indirect,                    &
        l_use_seasalt_autoconv, l_use_seasalt_pm, l_dust,                &
!
! data fields in/out
        so2(:,:,1:), dms(:,:,1:),                                        &
        so4_aitken(:,:,1:), so4_accu(:,:,1:), so4_diss(:,:,1:),          &
        h2o2(:,:,1:), nh3(:,:,1:),                                       &
        soot_new(:,:,1:), soot_agd(:,:,1:), soot_cld(:,:,1:),            &
        bmass_new(:,:,1:), bmass_agd(:,:,1:), bmass_cld(:,:,1:),         &
        ocff_new(:,:,1:), ocff_agd(:,:,1:), ocff_cld(:,:,1:),            &
        biogenic(:,:,1:),                                                &
        nitr_acc(:,:,1:), nitr_diss(:,:,1:),                             &
!
! data fields in
        dust_div1(:,:,1:), dust_div2(:,:,1:), dust_div3(:,:,1:),         &
        dust_div4(:,:,1:), dust_div5(:,:,1:), dust_div6(:,:,1:),         &
!
! data fields out
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
        stashwork17,                                                     &
! error info
        errorstatus )

       IF(l_sulpc_online_oxidants .AND. l_sulpc_2_way_coupling) THEN 
! 
         CALL write_sulpc_oxidants(                                     &
               oh_conc(:,:,1:), h2o2_mmr(:,:,1:), ho2_conc(:,:,1:),     &
               o3_mmr(:,:,1:), hno3_mmr(:,:,1:),                        &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
        first_atmstep_call  ) 
! 
       END IF 
!
! Don't call Swap_bounds for updated fields
!
! Diagnostics STASHed for Aerosol section 17
!
! DEPENDS ON: stash
       CALL stash(a_sm,a_im,17,stashwork17,                              &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMGA is a subset of ARGDUMA, needed for generic interfacing into
! the STASH routine. See TYPDUMA for description of individual
! components
      ! Dump components and lengths
     &  A_FIXHD, A_INTHD,A_LEN_INTHD, A_REALHD,A_LEN_REALHD,            &
     &  A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                        &
     &  A_LOOKUP,A_LEN2_LOOKUP,                                         &
      ! STASH superarray
     &  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts,                       &
! ARGDUMGA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          errorstatus,cmessage)
!
       DEALLOCATE (stashwork17)
!
        IF(ALLOCATED(o3_mmr  )) DEALLOCATE(o3_mmr)
        IF(ALLOCATED(hno3_mmr)) DEALLOCATE(hno3_mmr)
        IF(ALLOCATED(h2o2_mmr)) DEALLOCATE(h2o2_mmr)
        IF(ALLOCATED(oh_conc )) DEALLOCATE(oh_conc)
        IF(ALLOCATED(ho2_conc)) DEALLOCATE(ho2_conc)

     END IF aeroif        ! END L_SULPC_SO2.OR.L_SOOT TEST
!
! DEPENDS ON: timer
     IF (ltimer) CALL TIMER('AS Aerosol Modelling',6)


IF (l_physics) THEN

! Are we using the PC2 cloud scheme?

  IF (l_pc2) THEN

! ----------------------------------------------------------------------
! PC2: Calculate condensation due to changes in temperature resulting
!    from adiabatic changes in pressure (mainly from vertical advection)
!    Also call checking routine from within this subroutine.
! ----------------------------------------------------------------------

    CALL atm_step_alloc_4A( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
        cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
        frac_control, r_u, r_v, r_w, errorstatus, 'allocPC2')

! DEPENDS ON: pc2_pressure_forcing
    CALL pc2_pressure_forcing (                                         &
              p,pstar,p_theta_levels(pdims_s%i_start,pdims_s%j_start,1),&
              rhc_row_length, rhc_rows,                                 &
              timestep, rhcpt, theta, cf_bulk,                          &
              cf_liquid, cf_frozen,                                     &
              q, qcl,qcf,                                               &
              exner_star,                                               &
              EXNER_THETA_LEVELS(pdims_s%i_start,pdims_s%j_start,1),    &
              ccb, cumulus, rhts, tlts, qtts, ptts,cf_area(1,1,1),      &
              t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres,       &
              cf_inc_pres, cfl_inc_pres,                                &
              cff_inc_pres, t_dini, q_dini, qcl_dini, qcf_dini,         &
              cf_dini, cfl_dini, cff_dini, l_mr_pc2, l_acf_cusack,      &
              l_cld_area)

    DEALLOCATE(rhts)
    DEALLOCATE(tlts)
    DEALLOCATE(qtts)
    DEALLOCATE(ptts)

  END IF ! L_pc2

! ----------------------------------------------------------------------
! Add ability to get increments from qt_bal_cld call and output in section 15
! ----------------------------------------------------------------------

  CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star, &
      cf_star, cfl_star, cff_star, 8)

  IF ( (.NOT. l_pc2) .OR. l_pc2_reset ) THEN
! ----------------------------------------------------------------------
! Call cloud scheme to make cloud consistent with moisture fields
! ----------------------------------------------------------------------

! DEPENDS ON: qt_bal_cld
    CALL qt_bal_cld( pstar,                                       &
        p_theta_levels(pdims_s%i_start,pdims_s%j_start,1),p,      &
        theta,exner_theta_levels,                                 &
        q,qcl,qcf,qcf2,                                           &
        rhcpt, rhc_row_length, rhc_rows, bl_levels,               &
        delta_lambda, delta_phi,                                  &
        fv_cos_theta_latitude,                                    &
        l_cld_area, l_acf_cusack, l_acf_brooks,                   &
        l_mcr_qcf2,                                               &
        l_mr_qtbalcld,                                            &
        ntml, cumulus, cf_area, cf_bulk, cf_liquid, cf_frozen,    &
        mype)

    ! calculate changes to T , q etc
    CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star, &
      cf_star, cfl_star, cff_star, 9)

    END IF  ! L_pc2 and L_pc2_reset

    IF (l_run_with_physics2) DEALLOCATE ( rhcpt )

    ! set zero level quantities for Endgame
    ! currently set equal to level 1
    theta    (:,:,0) = theta    (:,:,1)
    q        (:,:,0) = q        (:,:,1)
    qcl      (:,:,0) = qcl      (:,:,1)
    qcf      (:,:,0) = qcf      (:,:,1)
    cf_bulk  (:,:,0) = cf_bulk  (:,:,1)
    cf_liquid(:,:,0) = cf_liquid(:,:,1)
    cf_frozen(:,:,0) = cf_frozen(:,:,1)

    If (l_mcr_qcf2) qcf2(:,:,0) = qcf2(:,:,1)

  END IF   ! L_physics




! ----------------------------------------------------------------------
! Section 9.0a Calculate Mass and energy of atmosphere if required
!              using end of timestep values
!             Only done every energy correction period.
! ----------------------------------------------------------------------
! DEPENDS ON: timer
          IF(ltimer) CALL TIMER('AS Energy mass   ',5)
      IF (L_EMCORR) THEN

! Set energy correction for use in section 30
        energy_corr_now=A_REALHD(rh_energy_corr)

        IF (LENERGY) THEN
! zero fields to be calculated

          tot_energy_final = 0.0
          tot_dry_mass_final = 0.0
          tot_moist_final = 0.0

          CALL eng_mass_diag (                                           &
                            wet_levels,                                  &
                            delta_xi1,delta_xi2,                         &
                            theta , u, v,                                &
                            w, rho_r_sq_np1 , q,                         &
                            qcl, qcf,                                    &
                            wet_to_dry_n,                                &
                            exner_theta_levels,                          &
                            net_mflux,                                   &
                            a_realhd(rh_tot_mass_init),                  &
                            a_realhd(rh_tot_m_init),                     &
! logical to indicate mass and moist correction required
                            Lmass_corr,Lqt_corr,Lemq_print,              &
! energy correction timestep info
                            a_energysteps,timestep,                      &
! IN/OUT  results from calculations
                            tot_energy_final,tot_dry_mass_final,         &
                            tot_moist_final)

! DEPENDS ON: cal_eng_mass_corr_4A
          CALL cal_eng_mass_corr_4A (                                    &
                            a_energysteps,                               &
                            net_flux,                                    &
                            a_realhd(rh_tot_mass_init),                  &
                            a_realhd(rh_tot_energy_init),                &
                            a_realhd(rh_energy_corr),                    &
                            tot_energy_final )

! Swap initial energy and final energy.

          a_realhd(rh_tot_energy_init) = tot_energy_final

! Swap initial moisture and final moisture.

          a_realhd(rh_tot_m_init) = tot_moist_final


        END IF   ! LENERGY

      ELSE
! Set energy correction for use in section 30
        energy_corr_now=0.0
      END IF     ! L_EMCORR
! DEPENDS ON: timer

      IF(ltimer) CALL TIMER('AS Energy mass   ',6)

! DEPENDS ON: timer 
      IF (ltimer) CALL TIMER('AS IAU',5) 
!---------------------------------------------------------------------- 
! Incremental Analysis Update (IAU). 
!---------------------------------------------------------------------- 
 
      IF (L_IAU .AND. STEPim(a_im) >= IAU_FirstCallTS .AND. & 
                      STEPim(a_im) <= IAU_LastCallTS) THEN 


! DEPENDS ON: iau
        CALL iau (                                      & 
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                   l_mr_iau,                            & ! in 
                   u,      v,            w,             & ! inout 
                   u_adv,  v_adv,        w_adv,         & ! inout 
                   theta,  rho_r_sq_np1, murk,          & ! inout 
                   q,      qCL,          qCF,           & ! inout 
                   TStar,  TStar_tile,   Deep_soil_temp,& ! inout 
                   smcl,   tstar_anom,                  & ! inout 
                   Pstar,  p,                           & ! inout 
                   p_theta_levels,                      & ! inout 
                   exner_rho_levels,                    & ! inout 
                   exner_theta_levels,                  & ! inout 
                   snodep,                              & ! inout 
                   cf_area,                             & ! inout 
                   cf_bulk,                             & ! inout 
                   cf_liquid,                           & ! inout 
                   cf_frozen,                           & ! inout 
                   dust_div1, dust_div2, dust_div3,     & ! inout
                   dust_div4, dust_div5, dust_div6,     & ! inout
                   ozone_tracer )                         ! inout 


        ! IAU may have changed wet density. The dry density has to see this
        ! for conversions at the end to be accurate

        CALL eg_q_to_mix                                                &
                  (tdims_l,tdims_s,                                     &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup,                                 &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   m_v_np1, m_cl_np1, m_cf_np1                          &
                  ,m_cf2_np1, m_r_np1, m_gr_np1)

        DO k=1, model_levels
          DO j=pdims%j_start-offy, pdims%j_end+offy
            DO i=pdims%i_start-offx, pdims%i_end+offx

              mixing_r=(intw_w2rho(k,1)*(m_v_np1  (i,j,k)      +        &
                                         m_r_np1  (i,j,k)      +        &
                                         m_gr_np1 (i,j,k)      +        &
                                         m_cl_np1 (i,j,k)      +        &
                                         m_cf_np1 (i,j,k)      +        &
                                         m_cf2_np1(i,j,k))     +        &
                        intw_w2rho(k,2)*(m_v_np1  (i,j,k-1)    +        &
                                         m_r_np1  (i,j,k-1)    +        &
                                         m_gr_np1 (i,j,k-1)    +        &
                                         m_cl_np1 (i,j,k-1)    +        &
                                         m_cf_np1 (i,j,k-1)    +        &
                                         m_cf2_np1(i,j,k-1)))

              !     convert to dry rho
              rho(i,j,k) = rho_r_sq_np1(i,j,k)/(1.+mixing_r)    
 
              !     dividing by r**2
              rho(i,j,k) = rho(i,j,k)/( r_rho_levels(i,j,k)             &
                                       *r_rho_levels(i,j,k))
            END DO
          END DO
        END DO

        CALL eg_swap_bounds(rho,pdims_s,fld_type_p,.FALSE.)
 
      END IF 
 
! DEPENDS ON: timer 
      IF (LTIMER) CALL TIMER('AS IAU',6) 

! Are we using the PC2 cloud scheme to determine area cloud fraction?

  IF (l_pc2 .AND. .NOT. l_pc2_reset) THEN
    IF (.NOT. l_cld_area) THEN

! ----------------------------------------------------------------------
! PC2: Set area cloud fraction to the bulk cloud fraction. Use the
!    D1 arrays directly
! ----------------------------------------------------------------------

    DO k = 1, wet_levels
      DO j = j_start, j_stop
        DO i = i_start, i_stop
          cf_area(i,j,k) = cf_bulk(i,j,k)
        END DO
      END DO
    END DO

    ELSE IF (l_cld_area) THEN

      IF (l_acf_brooks) THEN
        ALLOCATE ( cf_bulk_nohalo(row_length,rows,wet_levels) )
        ALLOCATE ( cf_liquid_nohalo(row_length,rows,wet_levels) )
        ALLOCATE ( cf_frozen_nohalo(row_length,rows,wet_levels) )

! Place bulk, liquid and frozen cloud fractions in halo-free arrays
! Use indexing over the full row and row_length (including any LAM
! boundary rim) since the call to ls_acf_brooks uses this indexing.
        DO k = 1, wet_levels
          DO j = 1, rows
            DO i = 1, row_length
              cf_bulk_nohalo(i,j,k)   = cf_bulk(i,j,k)
              cf_liquid_nohalo(i,j,k) = cf_liquid(i,j,k)
              cf_frozen_nohalo(i,j,k) = cf_frozen(i,j,k)
            END DO
          END DO
        END DO

! DEPENDS ON: ls_acf_brooks
              CALL ls_acf_brooks (                                      &
                 delta_lambda, delta_phi,                               &
                 fv_cos_theta_latitude,                                 &
                 cf_bulk_nohalo, cf_liquid_nohalo,                      &
                 cf_frozen_nohalo, cumulus, cf_area )

        DEALLOCATE ( cf_bulk_nohalo )
        DEALLOCATE ( cf_liquid_nohalo )
        DEALLOCATE ( cf_frozen_nohalo )

      END IF ! L_ACF_Brooks

    END IF ! L_cld_area

  END IF  ! L_pc2 .and. .not. L_pc2_reset
!
!
! convert back to EG:
  IF (l_mr_pc2 .OR. l_mr_qtbalcld) THEN

  m_v (:,:,:)  =                                                        &
     q  (1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)
  m_cl(:,:,:)  =                                                        &
     qcl(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)
  m_cf(:,:,:)  =                                                        &
     qcf(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)

  IF (l_mcr_qcf2  )                                                     &
  m_cf2(:,:,:)  =                                                       &
     qcf2(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)

  IF (l_mcr_qrain )                                                     &
  m_r (:,:,:)  =                                                        &
     qrain(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)

  IF (l_mcr_qgraup)                                                     &
  m_gr(:,:,:)  =                                                        &
     qgraup(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels)

    CALL eg_mix_to_q                                                    &
                  (tdims_l,tdims_s,                                     &
                   m_v, m_cl, m_cf,                                     &
                   m_cf2, m_r, m_gr,                                    &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup)



  ELSE


    CALL eg_q_to_mix                                                    &
                  (tdims_l,tdims_s,                                     &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup,                                 &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   m_v, m_cl, m_cf                                      &
                  ,m_cf2, m_r, m_gr)

  END IF

! Nudging with analysis data
IF ( L_nudging ) THEN 
 
  IF (Ltimer) CALL timer('AS Nudging',5)
  IF ( model_domain /= mt_global ) THEN 
    ErrorStatus = 39 
    cmessage = 'Nudging not yet tested for Limited Area Model' 

    CALL ereport('ATM_STEP',ErrorStatus,cmessage) 
  END IF 

! Allocate stashwork array
  CALL Atm_Step_diag(                                             &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star,    &
       cf_star, cfl_star, cff_star, 39)

! DEPENDS ON: nudging_main1         
  CALL nudging_main1 (                                             & 
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
 
! time stepping information. 
       timestep, i_year, i_month, i_day, i_hour, i_minute,         & 
       i_second, timestep_number,                                  &  
 
! in data fields. 
       theta(:,:,1), u, v, p, exner_theta_levels, p_theta_levels,  & 
! out updated STASH 
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
        STASHwork39) 

    theta(:,:,0) = theta(:,:,1)      ! copy lev 1 to lev 0

! Swap bounds for modified fields
    CALL eg_swap_bounds(theta,tdims_s,fld_type_p,.FALSE.)
    CALL eg_swap_bounds(u,udims_s,fld_type_u,.TRUE.)
    CALL eg_swap_bounds(v,vdims_s,fld_type_v,.TRUE.)

! Copy diagnostics into main D1 aaray 
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',5)
! DEPENDS ON: Atm_Step_stash
          CALL Atm_Step_stash(                                     &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
           errorstatus,39)
! DEPENDS ON: timer
      IF (Ltimer) CALL timer('AS STASH',6)
          
      IF (Ltimer) CALL timer('AS Nudging',6)

   END IF                ! L_nudging 

! need to recompute thetav for EG

    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetav(i,j,k) = theta(i,j,k)*(1.0+m_v(i,j,k)*recip_epsilon)
        END DO
      END DO
    END DO

    CALL eg_swap_bounds(thetav,tdims_s,fld_type_p,.FALSE.)

! simulate dumping 
IF (ldump) THEN
  subseq_to_dump = .TRUE.
  CALL reset_dpt_pts()
END IF

! Scale back for compatabillity with top level code
DO k=1, model_levels
  DO j=pdims%j_start-offy, pdims%j_end+offy
    DO i=pdims%i_start-offx, pdims%i_end+offx

      mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)      +                  &
                                   m_r  (i,j,k)      +                  &
                                   m_gr (i,j,k)      +                  &
                                   m_cl (i,j,k)      +                  &
                                   m_cf (i,j,k)      +                  &
                                   m_cf2(i,j,k))     +                  &
                  intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                  &
                                   m_r  (i,j,k-1)    +                  &
                                   m_gr (i,j,k-1)    +                  &
                                   m_cl (i,j,k-1)    +                  &
                                   m_cf (i,j,k-1)    +                  &
                                   m_cf2(i,j,k-1)))

!     convert to wet rho
      rho(i,j,k) = rho(i,j,k)*(1.+mixing_r)    

!     multiply by r**2
      rho(i,j,k) = rho(i,j,k)*(r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
    END DO
  END DO
END DO

IF (integrity_test) THEN
  IF (me == 0) WRITE(0,fmt='(A)')                                       &
                     'End of  timestep (before chain dump) fields'
  CALL update_hash_m(                                                   &
                   exner,           SIZE(exner),           'pi___',     & 
                   rho,             SIZE(rho),             'rho__',     & 
                   thetav,          SIZE(thetav),          'tv___',     & 
                   u,               SIZE(u),               'u____',     &
                   v,               SIZE(v),               'v____',     &
                   w,               SIZE(w),               'w____',     &
                   exner_np1,       SIZE(exner_np1),       'pinp1',     &
                   u_np1,           SIZE(u_np1),           'u_np1',     &
                   v_np1,           SIZE(v_np1),           'v_np1',     &   
                   u_adv,           SIZE(u_adv),           'u_adv',     &
                   v_adv,           SIZE(v_adv),           'v_adv',     &
                   w_adv,           SIZE(w_adv),           'w_adv',     &
                   psi_w_surf,      SIZE(psi_w_surf),      'psiws',     &
                   rho_np1,         SIZE(rho_np1),         'r_np1',     &
                   m_v,             SIZE(m_v),             'm_v__',     &
                   m_cl,            SIZE(m_cl),            'm_cl_',     &
                   m_cf,            SIZE(m_cf),            'm_cf_',     &
                   m_r,             SIZE(m_r),             'm_r__',     &
                   m_gr,            SIZE(m_gr),            'm_grp',     &
                   m_cf2,           SIZE(m_cf2),           'm_cf2',     &
                   m_v_np1,         SIZE(m_v_np1),         'mvnp1',     &
                   m_cl_np1,        SIZE(m_cl_np1),        'mclp1',     &
                   m_cf_np1,        SIZE(m_cf_np1),        'mcfp1',     &
                   m_r_np1,         SIZE(m_r_np1),         'mrnp1',     &
                   m_gr_np1,        SIZE(m_gr_np1),        'mgrp1',     &
                   m_cf2_np1,       SIZE(m_cf2_np1),       'mcf21',     &
                   thetav_np1,      SIZE(thetav_np1),      'tvnp1',     &
                   w_np1,           SIZE(w_np1),           'w_np1',     &
                   etadot,          SIZE(etadot),          'ed___',     &
                   etadot_np1,      SIZE(etadot_np1),      'ednp1')
END IF



! ----------------------------------------------------------------------
! Section 9.0 Diagnostics at end of timestep
! ----------------------------------------------------------------------
! Calculation of total increments
! NOTE - this must be after all processes which alter model prognostics
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      IF(ltimer) CALL TIMER('AS End TStep Diags',5)

  CALL atm_step_diag( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
      r_v, r_u, r_w, q_star, qcl_star, qcf_star,  theta_star, &
      cf_star, cfl_star, cff_star, 10)

! section 15: 'dynamics' based quantities
IF(      sf(0,15) .AND. errorstatus == 0) THEN

! DEPENDS ON: st_diag1
  CALL st_diag1(stash_maxlen(15,a_im),                             &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
 t_incr_diagnostic,q_incr_diagnostic,qcl_incr_diagnostic,          &
 errorstatus,cmessage)

END IF ! SF(0,15)

! Cleanup arrays holding increments

IF (l_physics .AND. errorstatus == 0) THEN
  DEALLOCATE(t_incr_diagnostic)
  DEALLOCATE(q_incr_diagnostic)
  DEALLOCATE(qcl_incr_diagnostic)
END IF

! section 16: 'physics' based quantities
IF(      sf(0,16) .AND. errorstatus == 0) THEN

! DEPENDS ON: st_diag2
  CALL st_diag2(stash_maxlen(16,a_im),                             &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
 t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres,              &
 cf_inc_pres, cfl_inc_pres, cff_inc_pres,                         &
 t_dini, q_dini, qcl_dini, qcf_dini,                              &
 cf_dini, cfl_dini, cff_dini,                                     &
 errorstatus,cmessage)

END IF !   SF(0,16)

IF (l_pc2 .AND. errorstatus == 0 .AND. l_physics) THEN

  CALL atm_step_alloc_4A( &
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
       cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
       frac_control, r_u, r_v, r_w, errorstatus, 'dllocPC2')

END IF  ! L_pc2 and ErrorStatus == 0


! section 30: climate diagnostics
      IF(      SF(0,30) .AND. ErrorStatus == 0) THEN
! size of diagnostic space
        ALLOCATE (STASHwork30(STASH_maxlen(30,A_im)))
! DEPENDS ON: st_diag3
        CALL St_diag3(STASHwork30,STASH_maxlen(30,A_im),                 &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
          energy_corr_now,                                               &
          inc_u, inc_v, inc_w, inc_t,                                    &
          inc_q, inc_qcl, inc_qcf,                                       &
          inc_rho,sin_v_latitude,                                        &
          inc_qrain, inc_qgraup, inc_qcf2,                               &
          wet_to_dry_n,                                                  &
          ErrorStatus,CMessage)

        DEALLOCATE (STASHwork30) ! Clear space
        DEALLOCATE (inc_q)
        DEALLOCATE (inc_qcl)
        DEALLOCATE (inc_qcf)
        DEALLOCATE (inc_qrain)
        DEALLOCATE (inc_qgraup)
        DEALLOCATE (inc_qcf2)

      END IF ! SF(0,30)

IF( errorstatus == 0) THEN
  DEALLOCATE (inc_t)
  DEALLOCATE (inc_u)
  DEALLOCATE (inc_v)
  DEALLOCATE (inc_w)
  DEALLOCATE (inc_rho)
END IF

! Check error condition
IF(errorstatus >  0) THEN
   CALL ereport(routinename,errorstatus,cmessage)
END IF
! DEPENDS ON: timer
      IF(ltimer) CALL TIMER('AS End TStep Diags',6)

! DEPENDS ON: timer
      IF(ltimer) CALL TIMER('AS STASH',5)
! DEPENDS ON: Atm_Step_stash
CALL atm_step_stash( &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Start arg_atm_fields.h
! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!

! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                &
! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV,                                &
! 1.2: Data variables stored in secondary space.
       P,                                                                    &
! Pressure on theta levels
       P_THETA_LEVELS,                                                       &
! Exner pressure on theta levels
       EXNER_THETA_LEVELS,                                                   &
! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,                          &
! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF,                                     &
! 1.5: Radiation Increments
       SW_INCS, LW_INCS,                                                     &
! PAR radiation increment
       DIRPAR,                                                               &
! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP,          &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3,                                        &
!  tropopause-based ozone
       TPPSOZONE,                                                            &
! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK,                               &
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6,     &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3,                 &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD,        &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, CO2,                         &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5,           &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10,          &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15,      &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20,      &
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC,               &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,        &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC,            &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC,                       &
       DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                     &
       DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                     &
       SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,    &
       NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,               &
       BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                     &
       OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                        &
       NITR_ACC_LBC, NITR_DISS_LBC,                                     &
       TRACER_LBC,                                                      &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,&
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND,           &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND,                                 &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND,        &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND,  &
       MURK_LBC_TEND,                                                   &
       DUST_DIV1_LBC_TEND, DUST_DIV2_LBC_TEND, DUST_DIV3_LBC_TEND,      &
       DUST_DIV4_LBC_TEND, DUST_DIV5_LBC_TEND, DUST_DIV6_LBC_TEND,      &
       SO2_LBC_TEND, DMS_LBC_TEND, SO4_AITKEN_LBC_TEND,                 &
       SO4_ACCU_LBC_TEND, SO4_DISS_LBC_TEND, NH3_LBC_TEND,              &
       SOOT_NEW_LBC_TEND, SOOT_AGD_LBC_TEND, SOOT_CLD_LBC_TEND,         &
       BMASS_NEW_LBC_TEND, BMASS_AGD_LBC_TEND, BMASS_CLD_LBC_TEND,      &
       OCFF_NEW_LBC_TEND, OCFF_AGD_LBC_TEND, OCFF_CLD_LBC_TEND,         &
       NITR_ACC_LBC_TEND, NITR_DISS_LBC_TEND,                           &
       TRACER_LBC_TEND,                                                 &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, TSTAR_SICE_CAT,         &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB,                                                   &
! 2.2: Data variables stored in secondary space.
       PSTAR,                                                                &
! 2.3: Convective cloud fields
       CCB, CCT, CCLWP, deep_flag, past_precip, past_conv_ht,                &
! 2.4: Boundary layer fields
       ZH,                                                                   &
! Standard deviation of turbulent fluctuations of layer 1
       T1_SD,                                                                &
! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD,                                                                &
! Radiative screen-level temperatures
       TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans, &
! Number of model levels in the  turbulently mixed layer
       NTML,                                                                 &
! Top level for turb mixing in any decoupled Sc layer
       NTDSC,                                                                &
! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS,                                                       &
! Convective downdraught mass-flux at cloud base
       ddmfx,                                                                & 
! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT,               &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN,                 &
! 2.5: Other surface fields
       CANOPY_WATER, Z0, GS,                                                 &
! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2,                               &
       OROG_GRAD_X, OROG_GRAD_Y, OROG_UNFILT,                                &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY,                             &
! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS,              &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, ICE_K_CAT,                  &
! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT,         &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND,                                      &
       DUST_MREL1, DUST_MREL2, DUST_MREL3,                                   &
       DUST_MREL4, DUST_MREL5, DUST_MREL6,                                   &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM,               &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX,                           &
! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5,                &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10,               &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15,           &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20,           &
!   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX,                                                  &
!   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5,      &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9,                           &
       LAI_PFT, CANHT_PFT, DISTURB_VEG,                                      &
       SOIL_ALB, OBS_ALB_SW, OBS_ALB_VIS, OBS_ALB_NIR, SOIL_CARB,            &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4,                       &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC,                            &
       RSP_W_PFT_ACC, RSP_S_ACC,                                             &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4,                       &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE,                  &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, Z0H_TILE,                           &
       DOLR_FIELD,                                                           &
       LW_DOWN, SW_TILE_RTS,                                                 &
! Fields for MORUSES - new two-tile urban scheme
       HGT, HWR, WRR, DISP, ZTM, ALBWL, ALBRD, EMISW, EMISR,                 &
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS,                                                   &
!   2.15: Fields carried forward from previous version.
!         May not be required
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho,                             &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP,                                                &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND,                          &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET,                                &
!   Field required: water conservation correction due to lake evaporation 
       ACC_LAKE_EVAP,                                                   & 
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE,                             &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM,                             &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT,             &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN,       &
       C_SOLAR,C_BLUE,C_LONGWAVE,C_TAUX,C_TAUY,C_W10,                        &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_FCONDTOPN,C_TOPMELTN,C_LSRAIN,           &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_CALVING,                      &
!   2.18: JULES
       SNOWDEPTH, RHO_SNOW_GRND,                                             &
       NSNOW,                                                                &
       DS, SICE, SLIQ, TSNOWLAYER, RHO_SNOW, RGRAINL,                        &
! FLake lake scheme
       lake_depth, lake_fetch, lake_t_mean, lake_t_mxl,                      &
       lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape,                      &
       lake_g_dt,                                                            &
! UKCA oxidant fields
       OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA,                                &
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR,      &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK,      &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4,      &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC,      &
       ARCLDLTA_DL,                                                          &
! Convective Cloud Fields
       LCBASE, CCW_RAD,                                                      &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM,                    &
! UKCA tracer LBCs
       TRACER_UKCA_LBC, TRACER_UKCA_LBC_TEND,                                &
! Ammonium nitrate aerosol
       HNO3_UKCA, NITR_ACC, NITR_DISS,                                       &
! TKE based turbulent scheme
       E_TRB, TSQ_TRB,                                                       &
       QSQ_TRB, COV_TRB, ZHPAR_SHCU,                                         &
! ENDGame
       DryRho,EtaDot,ThetaV,psi_w_surf,psi_w_lid,m_v,m_cl,m_cf,m_cf2,m_r,    &
       m_gr,exner_surf,                                                      &
! End arg_atm_fields.h
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         errorstatus, 4)
! DEPENDS ON: timer
      IF(ltimer) CALL TIMER('AS STASH',6)


IF( l_diag_print ) THEN
         CALL print_diag_4A(                                            &
                        u, v, theta, rho_r_sq_np1, w, q, qcl, qcf,      &
                        rows, n_rows, row_length, model_levels,         &
                        offx, offy, timestep_number,                    &
                        print_step, diag_interval,                      &
                        rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
                        L_print_pe, L_print_w,                          &
                        L_print_wmax, L_print_max_wind,                 &
                        L_print_div, L_print_lapse, L_print_theta1,     &
                        L_print_shear, L_diag_wind, L_diag_noise,       &
                        max_w_run, max_wind_run, min_theta1_run,        &
                        dtheta1_run, max_div_run, min_div_run,          &
                        min_lapse_run, max_shear_run, time_max_shear,   &
                        time_div_max, time_div_min, time_lapse_min,     &
                        time_w_max, time_max_wind, time_theta1_min,     &
                        max_KE_run, min_KE_run, max_noise_run,          &
                        time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         IF ( l_flush6 ) CALL um_fort_flush(6,info)
END IF     !  L_diag_print


IF (integrity_test_ghash) THEN
   IF(me==0) OPEN (eg_unit,file='start_hash.dat', POSITION='APPEND')
  IF(me==0) THEN
    WRITE(eg_unit,fmt='(A)')   '############################'
    WRITE(eg_unit,fmt='(A,I5)') 'Global checksum for end of timestep ',&
                                timestep_number
  END IF
  CALL  global_hash(u,      'u         ',udims_s,fld_type_u,eg_unit)
  CALL  global_hash(v,      'v         ',vdims_s,fld_type_v,eg_unit)
  CALL  global_hash(w,      'w         ',wdims_s,fld_type_p,eg_unit)
  CALL  global_hash(etadot, 'etadot    ',wdims_s,fld_type_p,eg_unit)
  CALL  global_hash(thetav, 'thetav    ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(rho,    'tho       ',pdims_s,fld_type_p,eg_unit)
  CALL  global_hash(exner,  'exner     ',pdims_s,fld_type_p,eg_unit)
  CALL  global_hash(rho_r_sq_n,  'rho_r_sq_n',pdims_s,fld_type_p,eg_unit)
  !CALL  global_hash(exner_surf,'m_v       ',qdims_s,fld_type_p,eg_unit)
  !CALL  global_hash(psi_w_surf,'m_v       ',qdims_s,fld_type_p,eg_unit)
  IF (.NOT. l_dry) THEN
    CALL  global_hash(m_v,  'm_v       ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cl, 'm_cl      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cf, 'm_cf      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_r,  'm_r       ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_gr, 'm_gr      ',qdims_s,fld_type_p,eg_unit)
    CALL  global_hash(m_cf2,'m_cf2     ',qdims_s,fld_type_p,eg_unit)
  END IF
  CALL  global_hash(SO4_AITKEN,'so4 aitken',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SO4_ACCU,'acc       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SO4_DISS,'diss      ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(NH3,'nh3       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SOOT_NEW     ,'soot new  ',tdims_s,fld_type_p,eg_unit) 
  CALL  global_hash(SOOT_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(SOOT_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_NEW,'bmass     ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(BMASS_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_NEW,'ocff new  ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_AGD,'agd       ',tdims_s,fld_type_p,eg_unit)
  CALL  global_hash(OCFF_CLD,'cld       ',tdims_s,fld_type_p,eg_unit)
  IF(me==0) CLOSE (eg_unit)
END IF


! ----------------------------------------------------------------------
! Diagnostics at end of timestep
! ----------------------------------------------------------------------

! Logical first_atmstep_call is true on the first call to ATM_STEP_4A
! and is set to false at the end of ATM_STEP_4A (uses the SAVE command)
first_atmstep_call = .FALSE.

IF (Lexitnow) THEN
  CALL eg_destroy_vert_damp()
  CALL eg_destroy_horz_drag()
  CALL eg_destroy_helmholtz()
  CALL destroy_wet_to_dry_n()
  IF (l_skeb2) CALL destroy_r_skeb()
END IF

! DEPENDS ON: timer
IF (Ltimer) CALL timer('Atm_Step_4A (AS)',6)
IF (lhook) CALL dr_hook('ATM_STEP_4A',zhook_out,zhook_handle)

END SUBROUTINE atm_step_4A
