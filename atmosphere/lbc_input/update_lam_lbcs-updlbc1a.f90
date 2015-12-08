! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Update the Lateral Boundary Conditions (LBCs) of LAM fields
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input

      SUBROUTINE UPDATE_LAM_LBCs(                                       &
        r_rho_levels, r_theta_levels,                                   &
        ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,                 &
        TR_VARS,TR_LBC_VARS,TR_LEVELS,                                  &
        A_MAX_TRVARS,A_tr_active_lbc_index,                             &
        TR_UKCA,TR_LBC_UKCA,                                            &
        A_MAX_UKCAVARS,UKCA_tr_active_lbc_index,                        &
        OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                           &
        L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                   &
        L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc,   &
        L_murk, L_murk_lbc,                                             &
        L_LBC_balance, L_int_uvw_lbc,                                   &
        L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,     &
        L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,     &
        L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,     &
        L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,  &
        L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,         &
        L_nh3, L_nh3_lbc,                                               &
        L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,         &
        L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,       &
        L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,     &
        L_ocff_new, L_ocff_new_lbc,                                     &
        L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,         &
        L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,       &
        RIMWIDTH,RIMWEIGHTS,LENRIM,LBC_SIZE,LBC_START,                  &
        THETA_LBC, Q_LBC, QCL_LBC, QCF_LBC,                             &
        QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
        CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
        RHO_LBC, EXNER_LBC,                                             &
        U_LBC, V_LBC, W_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,           &
        MURK_LBC,                                                       &
        DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                    &
        DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                    &
        SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,   &
        NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,              &
        BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                    &
        OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                       &
        NITR_ACC_LBC, NITR_DISS_LBC,                                    &
        TRACER_LBCS, UKCA_TRACER_LBCS,                                  &
        QTYPE, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                 &
        CF_BULK, CF_LIQUID, CF_FROZEN,                                  &
        RHO, EXNER, U, V, W, U_ADV, V_ADV, W_ADV, MURK,                 &
        DUST_DIV1, DUST_DIV2, DUST_DIV3,                                &
        DUST_DIV4, DUST_DIV5, DUST_DIV6,                                &
        SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,                       &
        NH3, SOOT_NEW, SOOT_AGD, SOOT_CLD,                              &
        BMASS_NEW, BMASS_AGD, BMASS_CLD,                                &
        OCFF_NEW, OCFF_AGD, OCFF_CLD,                                   &
        NITR_ACC, NITR_DISS,                                            &
        DELTA_PHI, DELTA_LAMBDA,                                        &
        BASE_PHI, BASE_LAMBDA,                                          &
        DATASTART, lat_rot_NP,                                          &
        GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                 &
        FREE_TRACERS, UKCA_TRACERS                                      &
           )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE atm_fields_bounds_mod

      USE PrintStatus_mod
      USE UM_ParParams

      USE eg_balance_lbc_values_mod
      USE dynamics_input_mod, ONLY: l_endgame

      IMPLICIT NONE
! Updates all the input fields with their LBCs
!

! Parameters required for dimensioning some of the arguments
! Arguments

      INTEGER                                                           &
        ROW_LENGTH                                                      &
                          ! IN : Length of a model row
      , ROWS                                                            &
                          ! IN : Number of rows for theta,u fields
      , N_ROWS                                                          &
                          ! IN : Number of rows for v fields
      , MODEL_LEVELS                                                    &
                          ! IN : Number of model levels
      , WET_LEVELS                                                      &
                          ! IN : Number of moist model levels
      , OFFX                                                            &
                          ! IN : Size of "single" halo (EW direction)
      , OFFY                                                            &
                          ! IN : Size of "single" halo (NS direction)
      , HALO_I                                                          &
                          ! IN : Size of extended halo (EW direction)
      , HALO_J                                                          &
                          ! IN : Size of extended halo (NS direction)
      , RIMWIDTH                                                        &
                          ! IN : Size of boundary region
      , LENRIM(Nfld_max,NHalo_max)                                      &
                          ! IN : Size of single level of LBC
      , LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
                          ! IN : Size of a side of a LBC
      , LBC_START(4,Nfld_max,NHalo_max)                                 &
                          ! IN : Start of a side in a LBC
      , DATASTART(3)                                                    &
                          ! IN : position of first data point for PE
      , GLOBAL_ROW_LENGTH                                               &
                          ! IN : no. of columns in LAM
      , GLOBAL_ROWS                                                     &
                          ! IN : no. or rows in LAM
      , TR_VARS                                                         &
                          ! Number of tracer prognostics
      , TR_UKCA                                                         &
                          ! Number of UKCA tracer prognostics
      , TR_LBC_VARS                                                     &
                          ! Number of tracer lbcs
      , TR_LBC_UKCA                                                     &
                          ! Number of UKCA tracer lbcs
      , TR_LEVELS                                                       &
                          ! Number of tracer levels
      , A_MAX_TRVARS                                                    &
                          ! Max number of tracer prognostics
      , A_MAX_UKCAVARS                                                  &
                          ! Max number of UKCA tracer prognostics
      , A_tr_active_lbc_index(A_MAX_TRVARS)                             &
                          ! >0 if lbc active
      , UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)
                          ! >0 if ukca lbc active

      Logical, Intent (In) ::                                           &
        AT_EXTREMITY(4)                                                 &
                      ! IN : At an edge?
      , L_mcr_qcf2                                                      &
                      ! true if prognostic 2nd cloud ice active
      , L_mcr_qrain                                                     &
                      ! true if prognostic rain active
      , L_mcr_qgraup                                                    &
                      ! true if prognostic graupel active
      , L_mcr_qcf2_lbc                                                  &
                          ! true if prognostic 2nd cloud ice in lbcs
      , L_mcr_qrain_lbc                                                 &
                          ! true if prognostic rain in lbcs
      , L_mcr_qgraup_lbc                                                &
                          ! true if prognostic graupel in lbcs
      , L_pc2                                                           &
                          ! true if prognostic cloud fracs active
      , L_pc2_lbc                                                       &
                          ! true if prognostic cloud fracs in lbcs
      , L_murk                                                          &
                          ! true if murk aerosol active
      , L_murk_lbc                                                      &
                          ! true if murk aerosol lbcs active
      , L_int_uvw_lbc                                                   &
                          ! true if using interpolated advecting winds 
                          !  in lateral boundaries
      , L_LBC_balance                                                   & 
                          ! true if imposing balance in vertical 
                          !      momentum equation in LBC regions
      , L_dust_div1                                                     &
      , L_dust_div2                                                     &
      , L_dust_div3                                                     &
      , L_dust_div4                                                     &
      , L_dust_div5                                                     &
      , L_dust_div6                                                     &
                          ! true if DUST active
      , L_dust_div1_lbc                                                 &
      , L_dust_div2_lbc                                                 &
      , L_dust_div3_lbc                                                 &
      , L_dust_div4_lbc                                                 &
      , L_dust_div5_lbc                                                 &
      , L_dust_div6_lbc                                                 &
                          ! true if DUST lbcs active
      , L_so2                                                           &
                          ! true if SO2 active
      , L_so2_lbc                                                       &
                          ! true if SO2 lbcs active
      , L_dms                                                           &
                          ! true if DMS active
      , L_dms_lbc                                                       &
                          ! true if DMS lbcs aactive
      , L_so4_aitken                                                    &
                          ! true if SO4_AITKEN active
      , L_so4_aitken_lbc                                                &
                          ! true if SO4_AITKEN lbcs active
      , L_so4_accu                                                      &
                          ! true if SO4_ACCU active           
      , L_so4_accu_lbc                                                  &
                          ! true if SO4_ACCU lbcs active          
      , L_so4_diss                                                      &
                          ! true if SO4_DISS active          
      , L_so4_diss_lbc                                                  &
                          ! true if SO4_DISS lbcs active          
      , L_nh3                                                           &
                          ! true if NH3 active       
      , L_nh3_lbc                                                       &
                          ! true if NH3 lbcs active      
      , L_soot_new                                                      &
      , L_soot_agd                                                      &
      , L_soot_cld                                                      &
                          ! true if soot active  
      , L_soot_new_lbc                                                  &
      , L_soot_agd_lbc                                                  &
      , L_soot_cld_lbc                                                  &
                          ! true if soot lbcs active      
      , L_bmass_new                                                     &
      , L_bmass_agd                                                     &
      , L_bmass_cld                                                     &
                          ! true if biomass active
      , L_bmass_new_lbc                                                 &
      , L_bmass_agd_lbc                                                 &
      , L_bmass_cld_lbc                                                 &
                          ! true if biomass lbcs active
      , L_ocff_new                                                      &
      , L_ocff_agd                                                      &
      , L_ocff_cld                                                      &
                          ! true if fossil fuel aerosol active
      , L_ocff_new_lbc                                                  &
      , L_ocff_agd_lbc                                                  &
      , L_ocff_cld_lbc                                                  &
                          ! true if fossil fuel aerosol lbcs active
      , L_nitr_acc                                                      &
      , L_nitr_diss                                                     &
                          ! true if nitrate aerosol active
      , L_nitr_acc_lbc                                                  &
      , L_nitr_diss_lbc   ! true if nitrate aerosol lbcs active


      Real, Intent(In) ::                                               &
           ! vertical co-ordinate information
        r_theta_levels(1-halo_i:row_length+halo_i,                      &
                       1-halo_j:rows+halo_j,0:model_levels)             &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j, model_levels)

      Real, Intent (In) ::                                              &
        DELTA_PHI                                                       &
                    ! grid spacing (latitude)
      , DELTA_LAMBDA                                                    &
                    ! grid spacing (longitude)
      , BASE_PHI                                                        &
                    ! first latitude
      , BASE_LAMBDA                                                     &
                    ! first longitude
      , lat_rot_NP                                                      &
                    ! lat of rotated pole
      , RIMWEIGHTS(RIMWIDTH)  ! weight to apply to LBC

      REAL, INTENT (INOUT) ::                                           &
        THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
              qdims_l%k_start:qdims_l%k_end)                            &
      , QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
      , QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
      , QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
                qdims_l%k_start:qdims_l%k_end)                          &
      , QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                qdims_l%k_start:qdims_l%k_end)                          &
      , QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
                qdims_l%k_start:qdims_l%k_end)                          &
      , CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
                qdims_l%k_start:qdims_l%k_end)                          &
      , CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
                qdims_l%k_start:qdims_l%k_end)                          &
      , CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
                qdims_l%k_start:qdims_l%k_end)                          &
      , RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                pdims_s%k_start:pdims_s%k_end)                          &
      , EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  pdims_s%k_start:pdims_s%k_end+1)                      &
      , U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
              udims_s%k_start:udims_s%k_end)                            &
      , V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
              vdims_s%k_start:vdims_s%k_end)                            &
      , W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
              wdims_s%k_start:wdims_s%k_end)                            &
      , U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
                  udims_s%k_start:udims_s%k_end)                        &
      , V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
                  vdims_s%k_start:vdims_s%k_end)                        &
      , W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  wdims_s%k_start:wdims_s%k_end)                        &
      , MURK_LBC(LENRIM(fld_type_p,halo_type_single),                   &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV1_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV2_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV3_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV4_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV5_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DUST_DIV6_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SO2_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , DMS_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SO4_AITKEN_LBC(LENRIM(fld_type_p,halo_type_single),             &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SO4_ACCU_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SO4_DISS_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , NH3_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SOOT_NEW_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SOOT_AGD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , SOOT_CLD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , BMASS_NEW_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , BMASS_AGD_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , BMASS_CLD_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , OCFF_NEW_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , OCFF_AGD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , OCFF_CLD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , NITR_ACC_LBC(LENRIM(fld_type_p,halo_type_single),               &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , NITR_DISS_LBC(LENRIM(fld_type_p,halo_type_single),              &
                  tdims_s%k_start:tdims_s%k_end)                        &
      , TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),              &
                   trdims_ltl%k_start:trdims_ltl%k_end,TR_LBC_VARS)     &
      , UKCA_TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),         &
                   trdims_ltl%k_start:trdims_ltl%k_end,TR_LBC_UKCA)

! Type for moisture arrays (different between Endgame and ND).
      TYPE(array_dims), INTENT(IN) :: qtype

      REAL, INTENT (INOUT) ::                                           &
        THETA(tdims_s%i_start:tdims_s%i_end,                            &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end)                            &
      , Q(qtype%i_start:qtype%i_end,                                    &
          qtype%j_start:qtype%j_end,                                    &
          qtype%k_start:qtype%k_end)                                    &
      , QCL(qtype%i_start:qtype%i_end,                                  &
            qtype%j_start:qtype%j_end,                                  &
            qtype%k_start:qtype%k_end)                                  &
      , QCF(qtype%i_start:qtype%i_end,                                  &
            qtype%j_start:qtype%j_end,                                  &
            qtype%k_start:qtype%k_end)                                  &
      , QCF2(qtype%i_start:qtype%i_end,                                 &
             qtype%j_start:qtype%j_end,                                 &
             qtype%k_start:qtype%k_end)                                 &
      , QRAIN(qtype%i_start:qtype%i_end,                                &
              qtype%j_start:qtype%j_end,                                &
              qtype%k_start:qtype%k_end)                                &
      , QGRAUP(qtype%i_start:qtype%i_end,                               &
               qtype%j_start:qtype%j_end,                               &
               qtype%k_start:qtype%k_end)                               &
      , CF_BULK(qdims_l%i_start:qdims_l%i_end,                          &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      , CF_LIQUID(qdims_l%i_start:qdims_l%i_end,                        &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      , CF_FROZEN(qdims_l%i_start:qdims_l%i_end,                        &
              qdims_l%j_start:qdims_l%j_end,                            &
              qdims_l%k_start:qdims_l%k_end)                            &
      , RHO(pdims_s%i_start:pdims_s%i_end,                              &
          pdims_s%j_start:pdims_s%j_end,                                &
          pdims_s%k_start:pdims_s%k_end)                                &
      , EXNER(pdims_s%i_start:pdims_s%i_end,                            &
          pdims_s%j_start:pdims_s%j_end,                                &
          pdims_s%k_start:pdims_s%k_end+1)                              &
      , U(udims_s%i_start:udims_s%i_end,                                &
              udims_s%j_start:udims_s%j_end,                            &
              udims_s%k_start:udims_s%k_end)                            &
      , V(vdims_s%i_start:vdims_s%i_end,                                &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
      , W(wdims_s%i_start:wdims_s%i_end,                                &
          wdims_s%j_start:wdims_s%j_end,                                &
          wdims_s%k_start:wdims_s%k_end)                                &
      , U_ADV(udims_l%i_start:udims_l%i_end,                            &
              udims_l%j_start:udims_l%j_end,                            &
              udims_l%k_start:udims_l%k_end)                            &
      , V_ADV(vdims_l%i_start:vdims_l%i_end,                            &
              vdims_l%j_start:vdims_l%j_end,                            &
              vdims_l%k_start:vdims_l%k_end)                            &
      , W_ADV(wdims_l%i_start:wdims_l%i_end,                            &
          wdims_l%j_start:wdims_l%j_end,                                &
          wdims_l%k_start:wdims_l%k_end)                                &
      , MURK(tdims_s%i_start:tdims_s%i_end,                             &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)

      REAL, INTENT (INOUT) ::                                           &
        DUST_DIV1(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DUST_DIV2(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DUST_DIV3(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DUST_DIV4(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DUST_DIV5(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DUST_DIV6(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SO2(tdims_s%i_start:tdims_s%i_end,                              &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , DMS(tdims_s%i_start:tdims_s%i_end,                              &
          tdims_s%j_start:tdims_s%j_end,                                &
         tdims_s%k_start:tdims_s%k_end)                                 &
      , SO4_AITKEN(tdims_s%i_start:tdims_s%i_end,                       &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SO4_ACCU(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SO4_DISS(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , NH3(tdims_s%i_start:tdims_s%i_end,                              &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SOOT_NEW(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SOOT_AGD(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , SOOT_CLD(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , BMASS_NEW(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , BMASS_AGD(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , BMASS_CLD(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , OCFF_NEW(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , OCFF_AGD(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , OCFF_CLD(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , NITR_ACC(tdims_s%i_start:tdims_s%i_end,                         &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , NITR_DISS(tdims_s%i_start:tdims_s%i_end,                        &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
      , FREE_TRACERS(trdims_s%i_start:trdims_s%i_end,                   &
                     trdims_s%j_start:trdims_s%j_end,                   &
                     trdims_s%k_start:trdims_s%k_end, TR_VARS)          &
      , UKCA_TRACERS(trdims_s%i_start:trdims_s%i_end,                   &
                     trdims_s%j_start:trdims_s%j_end,                   &
                     trdims_s%k_start:trdims_s%k_end,TR_UKCA)

! Local variables

      LOGICAL                                                           &
        L_Do_Boundaries                                                 &
      , L_Do_Halos
      INTEGER                                                           &
        N_RIMS_TO_DO                                                    &
      , lbc_index                                                       &
      , tracer

! This pointer is used to distinguish between EG and ND
! moisture field sizes (avoids changing the argument lists)


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

       IF (lhook) CALL dr_hook('UPDATE_LAM_LBCS',zhook_in,zhook_handle)

       L_Do_Boundaries = .TRUE.
       L_Do_Halos = .TRUE.
       N_RIMS_TO_DO = RIMWIDTH

       IF(L_LBC_balance)THEN

! Reset LBC array values for Exner, rho and w 

       IF(l_endgame) THEN
         CALL EG_BALANCE_LBC_VALUES(                                    &
          EXNER_LBC, RHO_LBC, THETA_LBC,Q_LBC, W_LBC, W_ADV_LBC,        &
          U_LBC, V_LBC,                                                 &
          QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,            &
          L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                        &
          R_RHO_LEVELS, R_THETA_LEVELS,                                 &
          ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I, HALO_J,   &
          LENRIM(fld_type_p,halo_type_extended),                        &
          LENRIM(fld_type_u,halo_type_extended),                        &
          LENRIM(fld_type_v,halo_type_extended),                        &
          LBC_START(1,fld_type_p,halo_type_extended),                   &
          LBC_START(1,fld_type_u,halo_type_extended),                   &
          LBC_START(1,fld_type_v,halo_type_extended),                   &
          RIMWIDTH, N_RIMS_TO_DO, RIMWEIGHTS, AT_EXTREMITY,             &
          DELTA_PHI, DELTA_LAMBDA,                                      &
          BASE_PHI, BASE_LAMBDA,                                        &
          DATASTART, LAT_ROT_NP,                                        &
          GLOBAL_ROW_LENGTH, GLOBAL_ROWS)

       ELSE
! DEPENDS ON: balance_lbc_values
         CALL BALANCE_LBC_VALUES(                                       &
          EXNER_LBC, RHO_LBC, THETA_LBC,Q_LBC, W_LBC, W_ADV_LBC,        &
          U_LBC, V_LBC,                                                 &
          R_RHO_LEVELS, R_THETA_LEVELS,                                 &
          ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I, HALO_J,   &
          LENRIM(fld_type_p,halo_type_extended),                        &
          LENRIM(fld_type_u,halo_type_extended),                        &
          LENRIM(fld_type_v,halo_type_extended),                        &
          LBC_START(1,fld_type_p,halo_type_extended),                   &
          LBC_START(1,fld_type_u,halo_type_extended),                   &
          LBC_START(1,fld_type_v,halo_type_extended),                   &
          RIMWIDTH, N_RIMS_TO_DO, RIMWEIGHTS, AT_EXTREMITY,             &
          DELTA_PHI, DELTA_LAMBDA,                                      &
          BASE_PHI, BASE_LAMBDA,                                        &
          DATASTART, LAT_ROT_NP,                                        &
          GLOBAL_ROW_LENGTH, GLOBAL_ROWS, L_int_uvw_lbc                 &
          )
       END IF

       END IF
       
! Theta
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        tdims  %i_end - tdims  %i_start + 1,                            &
        tdims  %j_end - tdims  %j_start + 1,                            &
        tdims_s%i_end - tdims  %i_end      ,                            &
        tdims_s%j_end - tdims  %j_end      ,                            &
        tdims_s%k_end - tdims_s%k_start + 1,                            &
        fld_type_p,THETA,                                               &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,THETA_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Q
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        tdims  %i_end - tdims  %i_start + 1,                            &
        tdims  %j_end - tdims  %j_start + 1,                            &
        qtype%i_end - tdims  %i_end,                                    &
        qtype%j_end - tdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,Q,                                                   &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,Q_LBC,                                            &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! QCL
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        tdims  %i_end - tdims  %i_start + 1,                            &
        tdims  %j_end - tdims  %j_start + 1,                            &
        qtype%i_end - tdims  %i_end,                                    &
        qtype%j_end - tdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,QCL,                                                 &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,QCL_LBC,                                          &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! QCF
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        tdims  %i_end - tdims  %i_start + 1,                            &
        tdims  %j_end - tdims  %j_start + 1,                            &
        qtype%i_end - tdims  %i_end,                                    &
        qtype%j_end - tdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,QCF,                                                 &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,QCF_LBC,                                          &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! QCF2 - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qcf2 .AND. L_mcr_qcf2_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qtype%i_end - qdims  %i_end,                                    &
        qtype%j_end - qdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,QCF2,                                                &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,QCF2_LBC,                                         &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)
      END IF

! QRAIN - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qrain .AND. L_mcr_qrain_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qtype%i_end - qdims  %i_end,                                    &
        qtype%j_end - qdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,QRAIN,                                               &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,QRAIN_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)
      END IF

! QGRAUP - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qgraup .AND. L_mcr_qgraup_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qtype%i_end - qdims  %i_end,                                    &
        qtype%j_end - qdims  %j_end,                                    &
        qtype%k_end - qtype%k_start + 1,                                &
        fld_type_p,QGRAUP,                                              &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,QGRAUP_LBC,                                       &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)
      END IF

! Cloud fractions - Update if prognostics are active and data in lbcs
      IF (L_pc2 .AND. L_pc2_lbc) THEN

! CF_BULK
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qdims_l%i_end - qdims  %i_end      ,                            &
        qdims_l%j_end - qdims  %j_end      ,                            &
        qdims_l%k_end - qdims_l%k_start + 1,                            &
        fld_type_p,CF_BULK,                                             &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,CF_BULK_LBC,                                      &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! CF_LIQUID
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qdims_l%i_end - qdims  %i_end      ,                            &
        qdims_l%j_end - qdims  %j_end      ,                            &
        qdims_l%k_end - qdims_l%k_start + 1,                            &
        fld_type_p,CF_LIQUID,                                           &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,CF_LIQUID_LBC,                                    &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! CF_FROZEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        qdims  %i_end - qdims  %i_start + 1,                            &
        qdims  %j_end - qdims  %j_start + 1,                            &
        qdims_l%i_end - qdims  %i_end      ,                            &
        qdims_l%j_end - qdims  %j_end      ,                            &
        qdims_l%k_end - qdims_l%k_start + 1,                            &
        fld_type_p,CF_FROZEN,                                           &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,CF_FROZEN_LBC,                                    &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

      END IF ! L_pc2 .AND. L_pc2_lbc

! RHO
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        pdims  %i_end - pdims  %i_start + 1,                            &
        pdims  %j_end - pdims  %j_start + 1,                            &
        pdims_s%i_end - pdims  %i_end      ,                            &
        pdims_s%j_end - pdims  %j_end      ,                            &
        pdims_s%k_end - pdims_s%k_start + 1,                            &
        fld_type_p,RHO,                                                 &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,RHO_LBC,                                          &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! EXNER
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        pdims  %i_end - pdims  %i_start + 1,                            &
        pdims  %j_end - pdims  %j_start + 1,                            &
        pdims_s%i_end - pdims  %i_end      ,                            &
        pdims_s%j_end - pdims  %j_end      ,                            &
        pdims_s%k_end - pdims_s%k_start + 2,                            &
        fld_type_p,EXNER,                                               &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        halo_i,halo_j,EXNER_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! U
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        udims  %i_end - udims  %i_start + 1,                            &
        udims  %j_end - udims  %j_start + 1,                            &
        udims_s%i_end - udims  %i_end      ,                            &
        udims_s%j_end - udims  %j_end      ,                            &
        udims_s%k_end - udims_s%k_start + 1,                            &
        fld_type_u,U,                                                   &
        LENRIM(fld_type_u,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
        LBC_START(1,fld_type_u,halo_type_extended),                     &
        HALO_I, HALO_J,U_LBC,                                           &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! V
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        vdims  %i_end - vdims  %i_start + 1,                            &
        vdims  %j_end - vdims  %j_start + 1,                            &
        vdims_s%i_end - vdims  %i_end      ,                            &
        vdims_s%j_end - vdims  %j_end      ,                            &
        vdims_s%k_end - vdims_s%k_start + 1,                            &
        fld_type_v,V,                                                   &
        LENRIM(fld_type_v,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
        LBC_START(1,fld_type_v,halo_type_extended),                     &
        HALO_I, HALO_J,V_LBC,                                           &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! W
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        wdims  %i_end - wdims  %i_start + 1,                            &
        wdims  %j_end - wdims  %j_start + 1,                            &
        wdims_s%i_end - wdims  %i_end      ,                            &
        wdims_s%j_end - wdims  %j_end      ,                            &
        wdims_s%k_end - wdims_s%k_start + 1,                            &
        fld_type_p,W,                                                   &
        LENRIM(fld_type_p,halo_type_extended),                          &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I, HALO_J,W_LBC,                                           &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

       IF ( .NOT. L_int_uvw_lbc ) THEN

! U_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        udims  %i_end - udims  %i_start + 1,                            &
        udims  %j_end - udims  %j_start + 1,                            &
        udims_l%i_end - udims  %i_end      ,                            &
        udims_l%j_end - udims  %j_end      ,                            &
        udims_l%k_end - udims_l%k_start + 1,                            &
        fld_type_u,                                                     &
        U_ADV,LENRIM(fld_type_u,halo_type_extended),                    &
        LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
        LBC_START(1,fld_type_u,halo_type_extended),                     &
        HALO_I,HALO_J,U_ADV_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! V_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        vdims  %i_end - vdims  %i_start + 1,                            &
        vdims  %j_end - vdims  %j_start + 1,                            &
        vdims_l%i_end - vdims  %i_end      ,                            &
        vdims_l%j_end - vdims  %j_end      ,                            &
        vdims_l%k_end - vdims_l%k_start + 1,                            &
        fld_type_v,                                                     &
        V_ADV,LENRIM(fld_type_v,halo_type_extended),                    &
        LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
        LBC_START(1,fld_type_v,halo_type_extended),                     &
        HALO_I,HALO_J,V_ADV_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

! W_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
        wdims  %i_end - wdims  %i_start + 1,                            &
        wdims  %j_end - wdims  %j_start + 1,                            &
        wdims_l%i_end - wdims  %i_end      ,                            &
        wdims_l%j_end - wdims  %j_end      ,                            &
        wdims_l%k_end - wdims_l%k_start + 1,                            &
        fld_type_p,                                                     &
        W_ADV,LENRIM(fld_type_p,halo_type_extended),                    &
        LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
        LBC_START(1,fld_type_p,halo_type_extended),                     &
        HALO_I,HALO_J,W_ADV_LBC,                                        &
        RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
        L_Do_Boundaries,L_Do_Halos)

       END IF ! .NOT. L_int_uvw_lbc

! MURK Aerosol
      ! If murk aerosol is in use and murk lbcs are active,
      IF (L_murk .AND. L_murk_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,MURK,                                               &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,MURK_LBC,                                            &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV1
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div1 .AND. L_dust_div1_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV1,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV1_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV2
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div2 .AND. L_dust_div2_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV2,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV2_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV3
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div3 .AND. L_dust_div3_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV3,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV3_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV4
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div4 .AND. L_dust_div4_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV4,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV4_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV5
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div5 .AND. L_dust_div5_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV5,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV5_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DUST_DIV6
      ! If dust is in use and dust lbcs are active,
      IF (L_dust_div6 .AND. L_dust_div6_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DUST_DIV6,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DUST_DIV6_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SO2
      ! If so2 is in use and so2 lbcs are active,
      IF (L_so2 .AND. L_so2_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SO2,                                                &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SO2_LBC,                                             &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! DMS
      ! If dms is in use and dms lbcs are active,
      IF (L_dms .AND. L_dms_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,DMS,                                                &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,DMS_LBC,                                             &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SO4_AITKEN
      ! If so4_aitken is in use and so4_aitken lbcs are active,
      IF (L_so4_aitken .AND. L_so4_aitken_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SO4_AITKEN,                                         &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SO4_AITKEN_LBC,                                      &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SO4_ACCU
      ! If so4_accu is in use and so4_accu lbcs are active,
      IF (L_so4_accu .AND. L_so4_accu_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SO4_ACCU,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SO4_ACCU_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SO4_DISS
      ! If so4_diss is in use and so4_diss lbcs are active,
      IF (L_so4_diss .AND. L_so4_diss_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SO4_DISS,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SO4_DISS_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! NH3
      ! If nh3 is in use and nh3 lbcs are active,
      IF (L_nh3 .AND. L_nh3_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,NH3,                                                &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,NH3_LBC,                                             &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SOOT_NEW
      ! If soot_new is in use and soot_new lbcs are active,
      IF (L_soot_new .AND. L_soot_new_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SOOT_NEW,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SOOT_NEW_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SOOT_AGD
      ! If soot_agd is in use and soot_agd lbcs are active,
      IF (L_soot_agd .AND. L_soot_agd_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SOOT_AGD,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SOOT_AGD_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! SOOT_CLD
      ! If soot_cld is in use and soot_cld lbcs are active,
      IF (L_soot_cld .AND. L_soot_cld_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,SOOT_CLD,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,SOOT_CLD_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! BMASS_NEW
      ! If bmass_new is in use and bmass_new lbcs are active,
      IF (L_bmass_new .AND. L_bmass_new_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,BMASS_NEW,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,BMASS_NEW_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! BMASS_AGD
      ! If bmass_agd is in use and bmass_agd lbcs are active,
      IF (L_bmass_agd .AND. L_bmass_agd_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,BMASS_AGD,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,BMASS_AGD_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! BMASS_CLD
      ! If bmass_cld is in use and bmass_cld lbcs are active,
      IF (L_bmass_cld .AND. L_bmass_cld_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,BMASS_CLD,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,BMASS_CLD_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! OCFF_NEW
      ! If ocff_new is in use and ocff_new lbcs are active,
      IF (L_ocff_new .AND. L_ocff_new_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,OCFF_NEW,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,OCFF_NEW_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! OCFF_AGD
      ! If ocff_agd is in use and ocff_agd lbcs are active,
      IF (L_ocff_agd .AND. L_ocff_agd_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,OCFF_AGD,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,OCFF_AGD_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! OCFF_CLD
      ! If ocff_cld is in use and ocff_cld lbcs are active,
      IF (L_ocff_cld .AND. L_ocff_cld_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,OCFF_CLD,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,OCFF_CLD_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! NITR_ACC
      ! If nitr_acc is in use and nitr_acc lbcs are active,
      IF (L_nitr_acc .AND. L_nitr_acc_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,NITR_ACC,                                           &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,NITR_ACC_LBC,                                        &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! NITR_DISS
      ! If nitr_diss is in use and nitr_diss lbcs are active,
      IF (L_nitr_diss .AND. L_nitr_diss_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
         tdims  %i_end - tdims  %i_start + 1,                           &
         tdims  %j_end - tdims  %j_start + 1,                           &
         tdims_s%i_end - tdims  %i_end      ,                           &
         tdims_s%j_end - tdims  %j_end      ,                           &
         tdims_s%k_end - tdims_s%k_start + 1,                           &
         fld_type_p,NITR_DISS,                                          &
         LENRIM(fld_type_p,halo_type_single),                           &
         LBC_SIZE(1,fld_type_p,halo_type_single),                       &
         LBC_START(1,fld_type_p,halo_type_single),                      &
         OFFX,OFFY,NITR_DISS_LBC,                                       &
         RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
         L_Do_Boundaries,L_Do_Halos)
      END IF

! TRACER LBCs
! If tracer lbcs are present, update relevant tracer field
      IF (TR_VARS > 0 .AND. TR_LBC_VARS > 0) THEN 
                                                 
        ! Loop over active tracer prognostics
        DO tracer = 1,TR_VARS                                           

          lbc_index = A_tr_active_lbc_index(tracer)

          ! If tracer lbc data is present, update the lbcs
          IF (lbc_index /= -1) THEN
            
            IF (printstatus >= prstatus_diag) WRITE(6,*)                &
              'UPDATE_LAM_LBCS:',tracer,lbc_index
            
            ! Update tracer field boundaries
! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
             trdims_xstl%i_end - trdims_xstl%i_start + 1,               &
             trdims_xstl%j_end - trdims_xstl%j_start + 1,               &
             trdims_s%i_end    - trdims_xstl%i_end      ,               &
             trdims_s%j_end    - trdims_xstl%j_end      ,               &
             trdims_s%k_end    - trdims_s%k_start + 1   ,               &
             fld_type_p,                                                &
             FREE_TRACERS(trdims_s%i_start:trdims_s%i_end,              &
                          trdims_s%j_start:trdims_s%j_end,              &
                          trdims_s%k_start:trdims_s%k_end,              &
                          tracer),                                      &
             LENRIM(fld_type_p,halo_type_extended),                     &
             LBC_SIZE(1,fld_type_p,halo_type_extended),                 &
             LBC_START(1,fld_type_p,halo_type_extended),                &
             HALO_I,HALO_J,TRACER_LBCS(1,1,lbc_index),                  &
             RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,             &
             L_Do_Boundaries,L_Do_Halos)

          END IF

        END DO                                                          

      END IF  ! on (TR_VARS > 0 .AND. TR_LBC_VARS > 0)

! UKCA TRACER LBCs
! If ukca tracer lbcs are present, update relevant tracer field
      IF (TR_UKCA > 0 .AND. TR_LBC_UKCA > 0) THEN 
                                                 
        ! Loop over active tracer prognostics
        DO tracer = 1,TR_UKCA                                           

          lbc_index = UKCA_tr_active_lbc_index(tracer)

          ! If tracer lbc data is present, update the lbcs
          IF (lbc_index /= -1) THEN
            
            IF (printstatus >= prstatus_diag) WRITE(6,*)                &
              'UPDATE_LAM_LBCS (UKCA):',tracer,lbc_index
            
            ! Update tracer field boundaries
! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
             trdims_xstl%i_end - trdims_xstl%i_start + 1,               &
             trdims_xstl%j_end - trdims_xstl%j_start + 1,               &
             trdims_s%i_end    - trdims_xstl%i_end      ,               &
             trdims_s%j_end    - trdims_xstl%j_end      ,               &
             trdims_s%k_end    - trdims_s%k_start + 1   ,               &
             fld_type_p,                                                &
             UKCA_TRACERS(trdims_s%i_start:trdims_s%i_end,              &
                          trdims_s%j_start:trdims_s%j_end,              &
                          trdims_s%k_start:trdims_s%k_end,              &
                          tracer),                                      &
             LENRIM(fld_type_p,halo_type_extended),                     &
             LBC_SIZE(1,fld_type_p,halo_type_extended),                 &
             LBC_START(1,fld_type_p,halo_type_extended),                &
             HALO_I,HALO_J,UKCA_TRACER_LBCS(1,1,lbc_index),             &
             RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,             &
             L_Do_Boundaries,L_Do_Halos)

          END IF

        END DO                                                          

      END IF  ! on (TR_UKCA > 0 .AND. TR_LBC_UKCA > 0)

       IF (lhook) CALL dr_hook('UPDATE_LAM_LBCS',zhook_out,zhook_handle)
       RETURN
       END SUBROUTINE UPDATE_LAM_LBCs
