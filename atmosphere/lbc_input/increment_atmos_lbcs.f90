! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine INCREMENT_ATMOS_LBCS
!
! Purpose : Increments atmosphere LBCs
!
! ---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input

      SUBROUTINE INCREMENT_ATMOS_LBCS(                                  &
     &  STEPS_TO_NEXT_UPDATE,                                           &
     &  LENRIM,                                                         &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  TR_LBC_VARS,TR_LEVELS,TR_LBC_UKCA,                              &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
     &  L_murk_lbc, L_int_uvw_lbc,                                      &
     &  L_dust_div1_lbc,L_dust_div2_lbc,                                &
     &  L_dust_div3_lbc,L_dust_div4_lbc,                                &
     &  L_dust_div5_lbc,L_dust_div6_lbc,                                &
     &  L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                           &
     &  L_so4_accu_lbc,L_so4_diss_lbc,                                  &
     &  L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                        &
     &  L_soot_cld_lbc,L_bmass_new_lbc,                                 &
     &  L_bmass_agd_lbc,L_bmass_cld_lbc,                                &
     &  L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,                   &
     &  L_nitr_acc_lbc,L_nitr_diss_lbc,                                 &
     &  U_LBC, U_LBC_TEND,                                              &
     &  V_LBC, V_LBC_TEND,                                              &
     &  W_LBC, W_LBC_TEND,                                              &
     &  RHO_LBC, RHO_LBC_TEND,                                          &
     &  THETA_LBC, THETA_LBC_TEND,                                      &
     &  Q_LBC, Q_LBC_TEND,                                              &
     &  QCL_LBC, QCL_LBC_TEND,                                          &
     &  QCF_LBC, QCF_LBC_TEND,                                          &
     &  QCF2_LBC, QCF2_LBC_TEND,                                        &
     &  QRAIN_LBC, QRAIN_LBC_TEND,                                      &
     &  QGRAUP_LBC, QGRAUP_LBC_TEND,                                    &
     &  CF_BULK_LBC, CF_BULK_LBC_TEND,                                  &
     &  CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,                              &
     &  CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,                              &
     &  EXNER_LBC, EXNER_LBC_TEND,                                      &
     &  U_ADV_LBC, U_ADV_LBC_TEND,                                      &
     &  V_ADV_LBC, V_ADV_LBC_TEND,                                      &
     &  W_ADV_LBC, W_ADV_LBC_TEND,                                      &
     &  MURK_LBC, MURK_LBC_TEND,                                        &
     &  DUST_DIV1_LBC,DUST_DIV1_LBC_TEND,                               &
     &  DUST_DIV2_LBC,DUST_DIV2_LBC_TEND,                               &
     &  DUST_DIV3_LBC,DUST_DIV3_LBC_TEND,                               &
     &  DUST_DIV4_LBC,DUST_DIV4_LBC_TEND,                               &
     &  DUST_DIV5_LBC,DUST_DIV5_LBC_TEND,                               &
     &  DUST_DIV6_LBC,DUST_DIV6_LBC_TEND,                               &
     &  SO2_LBC,SO2_LBC_TEND,                                           &
     &  DMS_LBC,DMS_LBC_TEND,                                           &
     &  SO4_AITKEN_LBC,SO4_AITKEN_LBC_TEND,                             &
     &  SO4_ACCU_LBC,SO4_ACCU_LBC_TEND,                                 &
     &  SO4_DISS_LBC,SO4_DISS_LBC_TEND,                                 &
     &  NH3_LBC,NH3_LBC_TEND,                                           &
     &  SOOT_NEW_LBC,SOOT_NEW_LBC_TEND,                                 &
     &  SOOT_AGD_LBC,SOOT_AGD_LBC_TEND,                                 &
     &  SOOT_CLD_LBC,SOOT_CLD_LBC_TEND,                                 &
     &  BMASS_NEW_LBC,BMASS_NEW_LBC_TEND,                               &
     &  BMASS_AGD_LBC,BMASS_AGD_LBC_TEND,                               &
     &  BMASS_CLD_LBC,BMASS_CLD_LBC_TEND,                               &
     &  OCFF_NEW_LBC,OCFF_NEW_LBC_TEND,                                 &
     &  OCFF_AGD_LBC,OCFF_AGD_LBC_TEND,                                 &
     &  OCFF_CLD_LBC,OCFF_CLD_LBC_TEND,                                 &
     &  NITR_ACC_LBC,NITR_ACC_LBC_TEND,                                 &
     &  NITR_DISS_LBC,NITR_DISS_LBC_TEND,                               &
     &  TRACER_LBCS, TRACER_LBCS_TEND,                                  &
     &  TRACER_UKCA_LBCS, TRACER_UKCA_LBCS_TEND,                        &
     &  RIM_STEPSA,ICODE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                       pdims_s, tdims_s, qdims_l,       &
                                       trdims_s
      USE ereport_mod, ONLY : ereport
      USE UM_ParParams
      IMPLICIT NONE

! Parameters required for argument declarations

! Arguments:

      INTEGER                                                           &
     &  STEPS_TO_NEXT_UPDATE                                            &
                                    ! IN : Number of timesteps before
                                    !      next LBC update
     &, LENRIM(Nfld_max,NHalo_Max)                                      &
                                    ! IN : Size of a level of LBC
     &, MODEL_LEVELS                                                    &
                                    ! IN : Number of model levels
     &, WET_LEVELS                                                      &
                                    ! IN : Number of wet model levels
     &, TR_LBC_VARS                                                     &
                                    ! IN : Number of tracer LBCs
     &, TR_LBC_UKCA                                                     &
                                    ! IN : Number of UKCA tracer LBCs
     &, TR_LEVELS                                                       &
                                    ! IN : Number of tracer levels
     &, RIM_STEPSA                  ! IN : Number of timesteps in LBC
                                    !      period

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2_lbc                                                  &
                         ! true if using second cloud ice lbcs
     &, L_mcr_qrain_lbc                                                 &
                         ! true if using rain lbcs
     &, L_mcr_qgraup_lbc                                                &
                         ! true if using graupel lbcs
     &, L_pc2_lbc                                                       &
                         ! true if using cloud fraction lbcs
     &, L_int_uvw_lbc                                                   &
                          ! true if using interpolated advecting winds 
                          !  in lateral boundaries
     &, L_murk_lbc                                                      &
                          ! true if using murk aerosol lbcs
     &, L_dust_div1_lbc                                                 &
     &, L_dust_div2_lbc                                                 &
     &, L_dust_div3_lbc                                                 &
     &, L_dust_div4_lbc                                                 &
     &, L_dust_div5_lbc                                                 &
     &, L_dust_div6_lbc                                                 &
                         ! true if DUST in lbcs
     &, L_so2_lbc                                                       &
                         ! true if SO2 in lbcs
     &, L_dms_lbc                                                       &
                         ! true if DMS in lbcs
     &, L_so4_aitken_lbc                                                &
                         ! true if SO4_AITKEN in lbcs
     &, L_so4_accu_lbc                                                  &
                         ! true if SO4_ACCU in lbcs           
     &, L_so4_diss_lbc                                                  &
                         ! true if SO4_DISS in lbcs            
     &, L_nh3_lbc                                                       &
                         ! true if NH3 in lbcs       
     &, L_soot_new_lbc                                                  &
     &, L_soot_agd_lbc                                                  &
     &, L_soot_cld_lbc                                                  &
                         ! true if soot in lbcs       
     &, L_bmass_new_lbc                                                 &
     &, L_bmass_agd_lbc                                                 &
     &, L_bmass_cld_lbc                                                 &
                         ! true if biomass in lbcs
     &, L_ocff_new_lbc                                                  &
     &, L_ocff_agd_lbc                                                  &
     &, L_ocff_cld_lbc                                                  &
                         ! true if fossil fuel aerosol in lbcs
     &, L_nitr_acc_lbc                                                  &
     &, L_nitr_diss_lbc   ! true if nitrate aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_TENDs are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

      REAL, INTENT (INOUT) ::                                           &
        u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
              udims_s%k_start:udims_s%k_end)                            &
                                    ! IN/OUT : U LBC
      , u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
                          udims_s%k_start:udims_s%k_end)                &
                                    ! IN : U LBC tendency
      , v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
              vdims_s%k_start:vdims_s%k_end)                            &
                                    ! IN/OUT : V LBC
      , v_lbc_tend(lenrim(fld_type_v,halo_type_extended),               &
                   vdims_s%k_start:vdims_s%k_end)                       &
                                    ! IN : V LBC tendency
      , w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
              wdims_s%k_start:wdims_s%k_end)                            &
                                    ! IN/OUT : V LBC
      , w_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
                   wdims_s%k_start:wdims_s%k_end)                       &
                                    ! IN : V LBC tendency
      , rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                pdims_s%k_start:pdims_s%k_end)                          &
                                    ! IN/OUT : rho LBC
      , rho_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     pdims_s%k_start:pdims_s%k_end)                     &
                                    ! IN : rho LBC tendency
      , theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  tdims_s%k_start:tdims_s%k_end)                        &
                                    ! IN/OUT : theta LBC
      , theta_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! IN : theta LBC tendency
      , q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
              qdims_l%k_start:qdims_l%k_end)                            &
                                    ! IN/OUT : Q LBC
      , q_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
                   qdims_l%k_start:qdims_l%k_end)                       &
                                    ! IN : Q LBC tendency
      , qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! IN/OUT : QCL LBC
      , qcl_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     qdims_l%k_start:qdims_l%k_end)                     &
                                    ! IN : QCL LBC tendency
      , qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! IN/OUT : QCL LBC
      , qcf_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     qdims_l%k_start:qdims_l%k_end)                     &
                                    ! IN : QCL LBC tendency
      , qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
                 qdims_l%k_start:qdims_l%k_end)                         &
                                    ! IN/OUT : QCF2 LBC
      , qcf2_lbc_tend(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN : QCF2 LBC tendency
      , qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  qdims_l%k_start:qdims_l%k_end)                        &
                                    ! IN/OUT : QRAIN LBC
      , qrain_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       qdims_l%k_start:qdims_l%k_end)                   &
                                    ! IN : QRAIN LBC tendency
      , qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
                   qdims_l%k_start:qdims_l%k_end)                       &
                                    ! IN/OUT : QGRAUP LBC
      , qgraup_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                        qdims_l%k_start:qdims_l%k_end)                  &
                                    ! IN : QGRAUP LBC tendency
      , cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
                    qdims_l%k_start:qdims_l%k_end)                      &
                                    ! IN/OUT : CF_BULK LBC
      , cf_bulk_lbc_tend(lenrim(fld_type_p,halo_type_extended),         &
                         qdims_l%k_start:qdims_l%k_end)                 &
                                    ! IN : CF_BULK LBC tendency
      , cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN/OUT : CF_LIQUID LBC
      , cf_liquid_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                                  qdims_l%k_start:qdims_l%k_end)        &
                                    ! IN : CF_LIQUID LBC tendency
      , cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN/OUT : CF_FROZEN LBC
      , cf_frozen_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                           qdims_l%k_start:qdims_l%k_end)               &
                                    ! IN : CF_FROZEN LBC tendency
      , exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  pdims_s%k_start:pdims_s%k_end+1)                      &
                                    ! IN/OUT : Exner LBC
      , exner_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       pdims_s%k_start:pdims_s%k_end+1)                 &
                                         ! IN : Exner LBC tendency
      , u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
                  udims_s%k_start:udims_s%k_end)                        &
                                    ! IN/OUT : u_adv LBC
      , u_adv_lbc_tend(lenrim(fld_type_u,halo_type_extended),           &
                       udims_s%k_start:udims_s%k_end)                   &
                                    ! IN : u_adv LBC tendency
      , v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
                  vdims_s%k_start:vdims_s%k_end)                        &
                                    ! IN/OUT : v_adv LBC
      , v_adv_lbc_tend(lenrim(fld_type_v,halo_type_extended),           &
                       vdims_s%k_start:vdims_s%k_end)                   &
                                    ! IN : v_adv LBC tendency
      , w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
                         wdims_s%k_start:wdims_s%k_end)                 &
                                    ! IN/OUT : W LBC
      , w_adv_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       wdims_s%k_start:wdims_s%k_end)                   &
                                       ! IN : W LBC tendency
      , murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
                 tdims_s%k_start:tdims_s%k_end)                         &
                                    ! IN/OUT : MURK LBC
      , murk_lbc_tend(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end) ! IN : MURK LBC tendency

      REAL, INTENT (INOUT) ::                                           &    
        dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div1 LBC
      , dust_div1_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div1 LBC tendency
      , dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div2 LBC
      , dust_div2_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div2 LBC tendency
      , dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div3 LBC
      , dust_div3_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div3 LBC tendency
      , dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div4 LBC
      , dust_div4_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div4 LBC tendency
      , dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div5 LBC
      , dust_div5_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div5 LBC tendency
      , dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div6 LBC
      , dust_div6_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div6 LBC tendency
      , SO2_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : SO2 LBC
      , SO2_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : SO2 LBC tendency
      , dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : DMS LBC
      , dms_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : DMS LBC tendency
      , SO4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! IN/OUT : SO4_AITKEN LBC
      , SO4_aitken_lbc_tend(lenrim(fld_type_p,halo_type_single),        &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! IN : SO4_AITKEN LBC tendency
      , SO4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : SO4_accU LBC
      , SO4_accu_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : SO4_accU LBC tendency
      , SO4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : SO4_DISS LBC
      , SO4_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : SO4_DISS LBC tendency
      , NH3_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : NH3 LBC
      , NH3_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : NH3 LBC tendency
      , soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_NEW LBC
      , soot_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_NEW LBC tendency
      , soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_agd LBC
      , soot_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_agd LBC tendency
      , soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_CLD LBC
      , soot_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_CLD LBC tendency
      , bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_NEW LBC
      , bmass_new_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_NEW LBC tendency
      , bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_agd LBC
      , bmass_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_agd LBC tendency
      , bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_CLD LBC
      , bmass_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_CLD LBC tendency
      , ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_NEW LBC
      , ocff_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_NEW LBC tendency
      , ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_agd LBC
      , ocff_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_agd LBC tendency
      , ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_CLD LBC
      , ocff_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_CLD LBC tendency
      , nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : nitr_acc LBC
      , nitr_acc_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : nitr_acc LBC tendency
      , nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : nitr_DISS LBC
      , nitr_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : nitr_DISS LBC tendency
      , tracer_lbcs(lenrim(fld_type_p,halo_type_extended),              &
                          trdims_s%k_start:trdims_s%k_end,tr_lbc_vars)  &
                                           ! IN/OUT : Tracer LBCs
      , tracer_lbcs_tend(lenrim(fld_type_p,halo_type_extended),         &
                        trdims_s%k_start:trdims_s%k_end,tr_lbc_vars)    &
                                           ! IN : Tracer LBCs tendency  
      , tracer_ukca_lbcs(lenrim(fld_type_p,halo_type_extended),         &
                        trdims_s%k_start:trdims_s%k_end,tr_lbc_ukca)    &
                                           ! IN/OUT : UKCA Tracer LBCs
      , tracer_ukca_lbcs_tend(lenrim(fld_type_p,halo_type_extended),    &
                             trdims_s%k_start:trdims_s%k_end,tr_lbc_ukca)
                                           ! IN : UKCA Tracer LBCs tend


      INTEGER                                                           &
     &  ICODE                           ! Error code

! Local variables

      REAL                                                              &
     &  increment_factor   ! Factor to increment LBC by

      INTEGER                                                           &
     &  tracer                                                          &
                  ! loop counter over tracer variables
     &, k                                                               &
                  ! loop counter over levels
     &, i         ! loop counter over horizontal dimension

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ---------------------------------------------------------

      IF (lhook) CALL dr_hook('INCREMENT_ATMOS_LBCS',zhook_in,zhook_handle)

! Calculate increment_factor

      IF (STEPS_TO_NEXT_UPDATE  >   RIM_STEPSA) THEN
        WRITE(6,*) 'STEPS_TO_NEXT_UPDATE must be < ',RIM_STEPSA
        WRITE(6,*) 'Value supplied was ',STEPS_TO_NEXT_UPDATE
        ICODE=1

        Call Ereport("INCREMENT_ATMOS_LBCS ", ICODE,                    &
     &               "STEPS_TO_NEXT_UPDATE>BCSTEPS" )
      END IF

      increment_factor=1.0/STEPS_TO_NEXT_UPDATE

      DO k=udims_s%k_start,udims_s%k_end
        DO i=1,LENRIM(fld_type_u,halo_type_extended)
          U_LBC(i,k)=U_LBC(i,k)+increment_factor*                       &
     &                          (U_LBC_TEND(i,k)-U_LBC(i,k))
        END DO ! i
      END DO !k=udims_s%k_start,udims_s%k_end

      DO k=vdims_s%k_start,vdims_s%k_end
        DO i=1,LENRIM(fld_type_v,halo_type_extended)
          V_LBC(i,k)=V_LBC(i,k)+increment_factor*                       &
     &                          (V_LBC_TEND(i,k)-V_LBC(i,k))
        END DO ! i
      END DO !k=vdims_s%k_start,vdims_s%k_end

      DO k=pdims_s%k_start,pdims_s%k_end
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          RHO_LBC(i,k)=RHO_LBC(i,k)+increment_factor*                   &
                                    (RHO_LBC_TEND(i,k)-RHO_LBC(i,k))
        END DO ! i
      END DO ! k=pdims_s%k_start,pdims_s%k_end

      DO k=tdims_s%k_start,tdims_s%k_end
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          THETA_LBC(i,k)=THETA_LBC(i,k)+increment_factor*               &
     &                                  (THETA_LBC_TEND(i,k)-           &
     &                                   THETA_LBC(i,k))
        END DO ! i
      END DO ! k=tdims_s%k_start,tdims_s%k_end


      DO k=qdims_l%k_start,qdims_l%k_end
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          Q_LBC(i,k)=Q_LBC(i,k)+increment_factor*                       &
     &                          (Q_LBC_TEND(i,k)-Q_LBC(i,k))
          QCL_LBC(i,k)=QCL_LBC(i,k)+increment_factor*                   &
     &                              (QCL_LBC_TEND(i,k)-QCL_LBC(i,k))
          QCF_LBC(i,k)=QCF_LBC(i,k)+increment_factor*                   &
     &                              (QCF_LBC_TEND(i,k)-QCF_LBC(i,k))
        END DO ! i
      END DO ! k

      IF (L_mcr_qcf2_lbc) THEN
        DO k=qdims_l%k_start,qdims_l%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QCF2_LBC(i,k) = QCF2_LBC(i,k) + increment_factor *          &
     &                       (QCF2_LBC_TEND(i,k)-QCF2_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_mcr_qrain_lbc) THEN
        DO k=qdims_l%k_start,qdims_l%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QRAIN_LBC(i,k) = QRAIN_LBC(i,k) + increment_factor *        &
     &                       (QRAIN_LBC_TEND(i,k)-QRAIN_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_mcr_qgraup_lbc) THEN
        DO k=qdims_l%k_start,qdims_l%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QGRAUP_LBC(i,k) = QGRAUP_LBC(i,k) + increment_factor *      &
     &                         (QGRAUP_LBC_TEND(i,k)-QGRAUP_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_pc2_lbc) THEN  ! Cloud fractions are in lbcs
        DO k=qdims_l%k_start,qdims_l%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            CF_BULK_LBC(i,k) = CF_BULK_LBC(i,k) + increment_factor *    &
     &                     (CF_BULK_LBC_TEND(i,k)-CF_BULK_LBC(i,k))
            CF_LIQUID_LBC(i,k) = CF_LIQUID_LBC(i,k)+increment_factor*   &
     &                     (CF_LIQUID_LBC_TEND(i,k)-CF_LIQUID_LBC(i,k))
            CF_FROZEN_LBC(i,k) = CF_FROZEN_LBC(i,k)+increment_factor*   &
     &                     (CF_FROZEN_LBC_TEND(i,k)-CF_FROZEN_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF ! L_pc2_lbc

      DO k=wdims_s%k_start,wdims_s%k_end
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          W_LBC(i,k)=W_LBC(i,k)+increment_factor*                       &
     &                          (W_LBC_TEND(i,k)-W_LBC(i,k))
        END DO ! i
      END DO ! k

      DO k=pdims_s%k_start,pdims_s%k_end+1
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
           EXNER_LBC(i,k)=EXNER_LBC(i,k)+increment_factor*              &
     &                                  (EXNER_LBC_TEND(i,k)-           &
     &                                   EXNER_LBC(i,k))
        END DO ! i
      END DO ! k

      IF ( .not. L_int_uvw_lbc) THEN

        DO k = udims_s%k_start, udims_s%k_end
          DO i = 1, LENRIM(fld_type_u,halo_type_extended)
            U_ADV_LBC(i,k)=U_ADV_LBC(i,k)+increment_factor *            &
     &                                  (U_ADV_LBC_TEND(i,k) -          &
     &                                   U_ADV_LBC(i,k))
          END DO ! i
        END DO !k = udims_s%k_start, udims_s%k_end

        DO k = vdims_s%k_start, vdims_s%k_end
          DO i = 1, LENRIM(fld_type_v,halo_type_extended)
            V_ADV_LBC(i,k)=V_ADV_LBC(i,k)+increment_factor *            &
     &                                  (V_ADV_LBC_TEND(i,k) -          &
     &                                   V_ADV_LBC(i,k))
          END DO ! i
        END DO !  k = vdims_s%k_start, vdims_s%k_end
 
        DO k = wdims_s%k_start,wdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            W_ADV_LBC(i,k)=W_ADV_LBC(i,k)+increment_factor *            &
     &                                  (W_ADV_LBC_TEND(i,k) -          &
     &                                   W_ADV_LBC(i,k))
          END DO ! i
        END DO ! k = wdims_s%k_start,wdims_s%k_end

      END IF !  .NOT. L_int_uvw_lbc

      IF (L_murk_lbc) THEN
      ! murk aerosol turned on and murk lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            MURK_LBC(i,k) = MURK_LBC(i,k) + increment_factor *          &
     &                       (MURK_LBC_TEND(i,k)-MURK_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_dust_div1_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV1_LBC(i,k) = DUST_DIV1_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV1_LBC_TEND(i,k)-DUST_DIV1_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_dust_div2_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV2_LBC(i,k) = DUST_DIV2_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV2_LBC_TEND(i,k)-DUST_DIV2_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

     IF (L_dust_div3_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV3_LBC(i,k) = DUST_DIV3_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV3_LBC_TEND(i,k)-DUST_DIV3_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_dust_div4_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV4_LBC(i,k) = DUST_DIV4_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV4_LBC_TEND(i,k)-DUST_DIV4_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

     IF (L_dust_div5_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV5_LBC(i,k) = DUST_DIV5_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV5_LBC_TEND(i,k)-DUST_DIV5_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_dust_div6_lbc) THEN
      ! dust lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           DUST_DIV6_LBC(i,k) = DUST_DIV6_LBC(i,k) + increment_factor * &
     &                 (DUST_DIV6_LBC_TEND(i,k)-DUST_DIV6_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_so2_lbc) THEN
      ! so2 lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SO2_LBC(i,k) = SO2_LBC(i,k) + increment_factor *            &
     &                             (SO2_LBC_TEND(i,k)-SO2_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_dms_lbc) THEN
      ! dms lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            DMS_LBC(i,k) = DMS_LBC(i,k) + increment_factor *            &
     &                             (DMS_LBC_TEND(i,k)-DMS_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_so4_aitken_lbc) THEN
      ! so4_aitken lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
         SO4_AITKEN_LBC(i,k) = SO4_AITKEN_LBC(i,k) + increment_factor * &
     &               (SO4_AITKEN_LBC_TEND(i,k)-SO4_AITKEN_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_so4_accu_lbc) THEN
      ! so4_accu lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SO4_ACCU_LBC(i,k) = SO4_ACCU_LBC(i,k) + increment_factor *  &
     &               (SO4_ACCU_LBC_TEND(i,k)-SO4_ACCU_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_so4_diss_lbc) THEN
      ! so4_diss lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SO4_DISS_LBC(i,k) = SO4_DISS_LBC(i,k) + increment_factor *  &
     &               (SO4_DISS_LBC_TEND(i,k)-SO4_DISS_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_nh3_lbc) THEN
      ! nh3 lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            NH3_LBC(i,k) = NH3_LBC(i,k) + increment_factor *            &
     &               (NH3_LBC_TEND(i,k)-NH3_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_soot_new_lbc) THEN
      ! soot_new lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SOOT_NEW_LBC(i,k) = SOOT_NEW_LBC(i,k) + increment_factor *  &
     &               (SOOT_NEW_LBC_TEND(i,k)-SOOT_NEW_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_soot_agd_lbc) THEN
      ! soot_agd lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SOOT_AGD_LBC(i,k) = SOOT_AGD_LBC(i,k) + increment_factor *  &
     &               (SOOT_AGD_LBC_TEND(i,k)-SOOT_AGD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_soot_cld_lbc) THEN
      ! soot_cld lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            SOOT_CLD_LBC(i,k) = SOOT_CLD_LBC(i,k) + increment_factor *  &
     &               (SOOT_CLD_LBC_TEND(i,k)-SOOT_CLD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_bmass_new_lbc) THEN
      ! bmass_new lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           BMASS_NEW_LBC(i,k) = BMASS_NEW_LBC(i,k) + increment_factor * &
     &               (BMASS_NEW_LBC_TEND(i,k)-BMASS_NEW_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_bmass_agd_lbc) THEN
      ! bmass_agd lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           BMASS_AGD_LBC(i,k) = BMASS_AGD_LBC(i,k) + increment_factor * &
     &               (BMASS_AGD_LBC_TEND(i,k)-BMASS_AGD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_bmass_cld_lbc) THEN
      ! bmass_cld lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
           BMASS_CLD_LBC(i,k) = BMASS_CLD_LBC(i,k) + increment_factor * &
     &               (BMASS_CLD_LBC_TEND(i,k)-BMASS_CLD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_ocff_new_lbc) THEN
      ! ocff_new lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            OCFF_NEW_LBC(i,k) = OCFF_NEW_LBC(i,k) + increment_factor *  &
     &               (OCFF_NEW_LBC_TEND(i,k)-OCFF_NEW_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_ocff_agd_lbc) THEN
      ! ocff_agd lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            OCFF_AGD_LBC(i,k) = OCFF_AGD_LBC(i,k) + increment_factor *  &
     &               (OCFF_AGD_LBC_TEND(i,k)-OCFF_AGD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_ocff_cld_lbc) THEN
      ! ocff_cld lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            OCFF_CLD_LBC(i,k) = OCFF_CLD_LBC(i,k) + increment_factor *  &
     &               (OCFF_CLD_LBC_TEND(i,k)-OCFF_CLD_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_nitr_acc_lbc) THEN
      ! nitr_new lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            NITR_ACC_LBC(i,k) = NITR_ACC_LBC(i,k) + increment_factor *  &
     &               (NITR_ACC_LBC_TEND(i,k)-NITR_ACC_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF

      IF (L_nitr_diss_lbc) THEN
      ! nitr_diss lbcs active
        DO k=tdims_s%k_start,tdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            NITR_DISS_LBC(i,k) = NITR_DISS_LBC(i,k) + increment_factor *  &
     &               (NITR_DISS_LBC_TEND(i,k)-NITR_DISS_LBC(i,k))
          END DO ! i
        END DO ! k
      END IF
      
      DO tracer=1,TR_LBC_VARS
        DO k=trdims_s%k_start,trdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            TRACER_LBCS(i,k,tracer)=TRACER_LBCS(i,k,tracer)+            &
     &        increment_factor*                                         &
     &        (TRACER_LBCS_TEND(i,k,tracer)-TRACER_LBCS(i,k,tracer))
          END DO
        END DO
      END DO

      DO tracer=1,TR_LBC_UKCA
        DO k=trdims_s%k_start,trdims_s%k_end
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            TRACER_UKCA_LBCS(i,k,tracer)=TRACER_UKCA_LBCS(i,k,tracer)+  &
     &        increment_factor*                                         &
     &  (TRACER_UKCA_LBCS_TEND(i,k,tracer)-TRACER_UKCA_LBCS(i,k,tracer))
          END DO
        END DO 
      END DO


 9999 CONTINUE
      IF (lhook) CALL dr_hook('INCREMENT_ATMOS_LBCS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INCREMENT_ATMOS_LBCS
