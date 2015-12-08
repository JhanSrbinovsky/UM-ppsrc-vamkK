! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine READ_ATMOS_LBCS
!
! Purpose : Reads in all the Atmos LBCs
!
! ---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input

      SUBROUTINE READ_ATMOS_LBCS(                                       &
     &  LENRIM, global_LENRIM,                                          &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  TR_LBC_VARS,TR_LEVELS,TR_LBC_UKCA,                              &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, l_pc2,                   &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
     &  L_murk, L_murk_lbc,                                             &
     &  L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,     &
     &  L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,     &
     &  L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,     &
     &  L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,  &
     &  L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,         &
     &  L_nh3, L_nh3_lbc,                                               &
     &  L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,         &
     &  L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,       &
     &  L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,     &
     &  L_ocff_new, L_ocff_new_lbc,                                     &
     &  L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,         &
     &  L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,       &
     &  LEN1_LOOKUP,LEN_FIXHD,NFTIN, N_LOOKUPS,LOOKUP, FIXHD,           &
     &  U_LBC,V_LBC,W_LBC,RHO_LBC,THETA_LBC,Q_LBC,QCL_LBC,QCF_LBC,      &
     &  QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
     &  CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
     &  EXNER_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,                     &
     &  MURK_LBC,                                                       &
     &  DUST_DIV1_LBC, DUST_DIV2_LBC, DUST_DIV3_LBC,                    &
     &  DUST_DIV4_LBC, DUST_DIV5_LBC, DUST_DIV6_LBC,                    &
     &  SO2_LBC, DMS_LBC, SO4_AITKEN_LBC, SO4_ACCU_LBC, SO4_DISS_LBC,   &
     &  NH3_LBC, SOOT_NEW_LBC, SOOT_AGD_LBC, SOOT_CLD_LBC,              &
     &  BMASS_NEW_LBC, BMASS_AGD_LBC, BMASS_CLD_LBC,                    &
     &  OCFF_NEW_LBC, OCFF_AGD_LBC, OCFF_CLD_LBC,                       &
     &  NITR_ACC_LBC, NITR_DISS_LBC,                                    &
     &  TRACER_LBCS, UKCA_TRACER_LBCS,                                  &
     &  ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                       pdims_s, tdims_s, qdims_l,       &
                                       trdims_s
      USE UM_ParParams
      USE lookup_addresses
      USE dynamics_input_mod, ONLY: l_endgame

      IMPLICIT NONE

! Parameters required for argument declarations

! Arguments:

      INTEGER                                                           &
     &  LENRIM(Nfld_max,NHalo_Max)                                      &
                                    ! IN : Size of a level of LBC
     &, global_LENRIM(Nfld_max,NHalo_Max)                               &
                                    ! IN : Full (disk) size of a
                                    !      level of LBC
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
     &, LEN1_LOOKUP                                                     &
                                    ! IN : Size of LOOKUP header
     &, LEN_FIXHD                                                       &
                                    ! IN : Size of fixed length header
     &, NFTIN                                                           &
                                    ! IN : Unit for LBC file
     &, N_LOOKUPS                                                       &
                                    ! IN  : Number of lookup headers
     &, LOOKUP(LEN1_LOOKUP,N_LOOKUPS)                                   &
                                    ! IN : LOOKUP headers
     &, FIXHD(LEN_FIXHD)            ! IN : Fixed header

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup                                                    &
                          ! true if graupel active
     &, L_mcr_qcf2_lbc                                                  &
                          ! true if second cloud ice lbcs active
     &, L_mcr_qrain_lbc                                                 &
                          ! true if rain lbcs active
     &, L_mcr_qgraup_lbc                                                &
                          ! true if graupel lbcs active
     &, L_pc2                                                           &
                          ! true if prognostic cloud fractions used
     &, L_pc2_lbc                                                       &
                          ! true if cloud fractions in lbcs
     &, L_murk                                                          &
                          ! true if prognostic murk aerosol active
     &, L_murk_lbc                                                      &
                          ! true if murk aerosol in lbcs
     &, L_dust_div1                                                     &
     &, L_dust_div2                                                     &
     &, L_dust_div3                                                     &
     &, L_dust_div4                                                     &
     &, L_dust_div5                                                     &
     &, L_dust_div6                                                     &
                          ! true if DUST active
     &, L_dust_div1_lbc                                                 &
     &, L_dust_div2_lbc                                                 &
     &, L_dust_div3_lbc                                                 &
     &, L_dust_div4_lbc                                                 &
     &, L_dust_div5_lbc                                                 &
     &, L_dust_div6_lbc                                                 &
                          ! true if DUST lbcs active
     &, L_so2                                                           &
                          ! true if SO2 active
     &, L_so2_lbc                                                       &
                          ! true if SO2 lbcs active
     &, L_dms                                                           &
                          ! true if DMS active
     &, L_dms_lbc                                                       &
                          ! true if DMS lbcs aactive
     &, L_so4_aitken                                                    &
                          ! true if SO4_AITKEN active
     &, L_so4_aitken_lbc                                                &
                          ! true if SO4_AITKEN lbcs active
     &, L_so4_accu                                                      &
                          ! true if SO4_ACCU active           
     &, L_so4_accu_lbc                                                  &
                          ! true if SO4_ACCU lbcs active          
     &, L_so4_diss                                                      &
                          ! true if SO4_DISS active          
     &, L_so4_diss_lbc                                                  &
                          ! true if SO4_DISS lbcs active          
     &, L_nh3                                                           &
                          ! true if NH3 active       
     &, L_nh3_lbc                                                       &
                          ! true if NH3 lbcs active      
     &, L_soot_new                                                      &
     &, L_soot_agd                                                      &
     &, L_soot_cld                                                      &
                          ! true if soot active  
     &, L_soot_new_lbc                                                  &
     &, L_soot_agd_lbc                                                  &
     &, L_soot_cld_lbc                                                  &
                          ! true if soot lbcs active      
     &, L_bmass_new                                                     &
     &, L_bmass_agd                                                     &
     &, L_bmass_cld                                                     &
                          ! true if biomass active
     &, L_bmass_new_lbc                                                 &
     &, L_bmass_agd_lbc                                                 &
     &, L_bmass_cld_lbc                                                 &
                          ! true if biomass lbcs active
     &, L_ocff_new                                                      &
     &, L_ocff_agd                                                      &
     &, L_ocff_cld                                                      &
                          ! true if fossil fuel aerosol active
     &, L_ocff_new_lbc                                                  &
     &, L_ocff_agd_lbc                                                  &
     &, L_ocff_cld_lbc                                                  &
                          ! true if fossil fuel aerosol lbcs active
     &, L_nitr_acc                                                      &
     &, L_nitr_diss                                                     &
                          ! true if nitrate aerosol active
     &, L_nitr_acc_lbc                                                  &
     &, L_nitr_diss_lbc    ! true if nitrate aerosol lbcs active

      REAL                                                              &
        U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
              udims_s%k_start:udims_s%k_end)                            &
                                    ! OUT : U LBC
      , V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
              vdims_s%k_start:vdims_s%k_end)                            &
                                    ! OUT : V LBC
      , W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
              wdims_s%k_start:wdims_s%k_end)                            &
                                    ! OUT : V LBC
      , RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                pdims_s%k_start:pdims_s%k_end)                          &
                                    ! OUT : Rho LBC
      , THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  tdims_s%k_start:tdims_s%k_end)                        &
                                    ! OUT : Theta LBC
      , Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
              qdims_l%k_start:qdims_l%k_end)                            &
                                    ! OUT : Q LBC
      , QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! OUT : QCL LBC
      , QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! OUT : QCL LBC
      , QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
                 qdims_l%k_start:qdims_l%k_end)                         &
                                    ! OUT : QCF2 LBC
      , QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  qdims_l%k_start:qdims_l%k_end)                        &
                                    ! OUT : QRAIN LBC
      , QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
                   qdims_l%k_start:qdims_l%k_end)                       &
                                    ! OUT : QGRAUP LBC
      , CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
                    qdims_l%k_start:qdims_l%k_end)                      &
                                    ! OUT : CF_BULK LBC
      , CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! OUT : CF_LIQUID LBC
      , CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! OUT : CF_FROZEN LBC
      , EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                  pdims_s%k_start:pdims_s%k_end+1)                      &
                                    ! OUT : Exner LBC
      , U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
                  udims_s%k_start:udims_s%k_end)                        &
                                    ! OUT : U_ADV LBC
      , V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
                  vdims_s%k_start:vdims_s%k_end)                        &
                                    ! OUT : V_ADV LBC
      , W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
                         wdims_s%k_start:wdims_s%k_end)                 &
                                    ! OUT : W LBC
      , MURK_LBC(LENRIM(fld_type_p,halo_type_single),                   &
                 tdims_s%k_start:tdims_s%k_end)                         &
                                       ! OUT : MURK LBC
      , DUST_DIV1_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV1 LBC
      , DUST_DIV2_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV2 LBC
      , DUST_DIV3_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV3 LBC
      , DUST_DIV4_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV4 LBC
      , DUST_DIV5_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV5 LBC
      , DUST_DIV6_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : DUST_DIV6 LBC
      , SO2_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! OUT : SO2 LBC
      , DMS_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! OUT : DMS LBC
      , SO4_AITKEN_LBC(LENRIM(fld_type_p,halo_type_single),             &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! OUT : SO4_AITKEN LBC
      , SO4_ACCU_LBC(LENRIM(fld_type_p,halo_type_single),               &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! OUT : SO4_ACCU LBC
      , SO4_DISS_LBC(LENRIM(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! OUT : SO4_DISS LBC
      , NH3_LBC(LENRIM(fld_type_p,halo_type_single),                    &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! OUT : NH3 LBC
      , SOOT_NEW_LBC(LENRIM(fld_type_p,halo_type_single),               &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! OUT : SOOT_NEW LBC
      , SOOT_AGD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! OUT : SOOT_AGD LBC
      , SOOT_CLD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! OUT : SOOT_CLD LBC
      , BMASS_NEW_LBC(LENRIM(fld_type_p,halo_type_single),              &
                             tdims_s%k_start:tdims_s%k_end)             &
                                    ! OUT : BMASS_NEW LBC
      , BMASS_AGD_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : BMASS_AGD LBC
      , BMASS_CLD_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : BMASS_CLD LBC
      , OCFF_NEW_LBC(LENRIM(fld_type_p,halo_type_single),               &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! OUT : OCFF_NEW LBC
      , OCFF_AGD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! OUT : OCFF_AGD LBC
      , OCFF_CLD_LBC(LENRIM(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! OUT : OCFF_CLD LBC
      , NITR_ACC_LBC(LENRIM(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! OUT : NITR_ACC LBC
      , NITR_DISS_LBC(LENRIM(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! OUT : NITR_DISS LBC
      , TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),              &
                    trdims_s%k_start:trdims_s%k_end,tr_lbc_vars)        &
                                    ! OUT : Tracer LBCs
      , UKCA_TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),         &
                        trdims_s%k_start:trdims_s%k_end,tr_lbc_ukca)
                                    ! OUT : UKCA Tracer LBCs

      INTEGER                                                           &
     &  ICODE                           ! Error code

      CHARACTER(LEN=80)                                                    &
     &  CMESSAGE                        ! Error message

! Local variables

      INTEGER                                                           &
     &  tracer                                                          &
                 ! loop counter over tracer variables
     &, lbc_num  ! number of active lbcs (used for microphysics lbcs)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ---------------------------------------------------------

      IF (lhook) CALL dr_hook('READ_ATMOS_LBCS',zhook_in,zhook_handle)

! Check the field about to be read is U_LBC - if not, something
! has gone wrong.
      IF ( LOOKUP(ITEM_CODE,1)  /=  31002) THEN
        WRITE(6,*) 'Expected to find U_LBC (31002) in LBC file '
        WRITE(6,*) 'But found ',LOOKUP(ITEM_CODE,1),' instead.'
        ICODE=1
        CMESSAGE='READ_ATMOS_LBCS : Found wrong field in LBC file'
        GOTO 9999
      ENDIF

! Now read in all the fields

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,1),LEN1_LOOKUP,U_LBC,            &
                    global_LENRIM(fld_type_u,halo_type_extended)*       &
                    (udims_s%k_end - udims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading U_LBC'
        GOTO 9999
      ENDIF


! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,2),LEN1_LOOKUP,V_LBC,            &
                    global_LENRIM(fld_type_v,halo_type_extended)*       &
                    (vdims_s%k_end - vdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading V_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,3),LEN1_LOOKUP,W_LBC,            &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (wdims_s%k_end - wdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading W_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,4),LEN1_LOOKUP,RHO_LBC,          &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (pdims_s%k_end - pdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading RHO_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,5),LEN1_LOOKUP,THETA_LBC,        &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (tdims_s%k_end - tdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading THETA_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,6),LEN1_LOOKUP,Q_LBC,            &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (qdims_l%k_end - qdims_l%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading Q_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,7),LEN1_LOOKUP,QCL_LBC,          &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (qdims_l%k_end - qdims_l%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCL_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,8),LEN1_LOOKUP,QCF_LBC,          &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (qdims_l%k_end - qdims_l%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCF_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,9),LEN1_LOOKUP,EXNER_LBC,        &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (pdims_s%k_end - pdims_s%k_start+2),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading EXNER_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,10),LEN1_LOOKUP,U_ADV_LBC,       &
                    global_LENRIM(fld_type_u,halo_type_extended)*       &
                    (udims_s%k_end - udims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading U_ADV_LBC'
        GOTO 9999
      ENDIF


! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,11),LEN1_LOOKUP,V_ADV_LBC,       &
                    global_LENRIM(fld_type_v,halo_type_extended)*       &
                    (vdims_s%k_end - vdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading V_ADV_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,12),LEN1_LOOKUP,W_ADV_LBC,       &
                    global_LENRIM(fld_type_p,halo_type_extended)*       &
                    (wdims_s%k_end - wdims_s%k_start+1),                &
                    FIXHD,                                              &
                    ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading W_ADV_LBC'
        GOTO 9999
      ENDIF

      ! set number of active lbcs for lookup table
      lbc_num = 12

      ! ----------------------------------------------------------------
      ! Read ice crystal (qcf2) lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qcf2_lbc) THEN

        ! qcf2 is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,QCF2_LBC,                             &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCF2_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_mcr_qcf2) THEN

        ! qcf2 prognostic is active but qcf2 lbcs are not in input file
        ! set lbc array to zero
        QCF2_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qcf2, L_mcr_qcf2_lbc

      ! ----------------------------------------------------------------
      ! Read rain lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qrain_lbc) THEN

        ! qrain is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,QRAIN_LBC,                            &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QRAIN_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_mcr_qrain) THEN

        ! qrain is active but qrain lbcs are not in input file
        ! set lbc array to zero
        QRAIN_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qrain, L_mcr_qrain_lbc

      ! ----------------------------------------------------------------
      ! Read graupel lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qgraup_lbc) THEN

        ! qgraup is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,QGRAUP_LBC,                           &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QGRAUP_LBC'
          GOTO 9999
        ENDIF

      ELSE  IF (L_mcr_qgraup) THEN

        ! qgraup is active but qgraup lbcs are not in input file
        ! set lbc array to zero
        QGRAUP_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qgraup, L_mcr_qgraup_lbc

      ! ----------------------------------------------------------------
      ! Read cloud fraction lbcs present in lbc file
      ! ----------------------------------------------------------------

      ! Set cloud fraction lbcs
      IF (L_pc2_lbc) THEN

        ! Cloud fractions variables are in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,CF_BULK_LBC,                          &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_BULK_LBC'
          GOTO 9999
        ENDIF

       lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,CF_LIQUID_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_LIQUID_LBC'
          GOTO 9999
        ENDIF

        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,CF_FROZEN_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (qdims_l%k_end - qdims_l%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_FROZEN_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_pc2) THEN

        ! Cloud fraction prognostics are active but cloud fraction lbcs
        ! are not in input file so set the lbc arrays to zero.
        CF_BULK_LBC(:,:)   = 0.0
        CF_LIQUID_LBC(:,:) = 0.0
        CF_FROZEN_LBC(:,:) = 0.0

      ENDIF  ! on L_pc2, L_pc2_lbc

      ! ----------------------------------------------------------------
      ! Read murk aerosol lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_murk_lbc) THEN   ! murk lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,MURK_LBC,                             &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading MURK_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_murk) THEN

        ! murk prognostic is active but murk lbcs are not in input file
        ! set lbc array to zero
        MURK_LBC(:,:) = 0.0

      ENDIF ! on L_murk_lbc, L_murk

      ! ----------------------------------------------------------------
      ! Read dust lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_dust_div1_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV1_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV1_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV1) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV1_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div1_lbc, L_dust_div1

      !-----------------------------------------------------------------

      IF (L_dust_div2_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV2_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV2_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV2) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV2_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div2_lbc, L_dust_div2

      !-----------------------------------------------------------------

      IF (L_dust_div3_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV3_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV3_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV3) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV3_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div3_lbc, L_dust_div3

      !-----------------------------------------------------------------

      IF (L_dust_div4_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV4_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV4_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV4) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV4_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div4_lbc, L_dust_div4

      !-----------------------------------------------------------------

      IF (L_dust_div5_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV5_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV5_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV5) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV5_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div5_lbc, L_dust_div5

      !-----------------------------------------------------------------

      IF (L_dust_div6_lbc) THEN   ! dust lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DUST_DIV6_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV6_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DUST_DIV6) THEN

        ! dust prognostic is active but dust lbcs are not in input file
        ! set lbc array to zero
        DUST_DIV6_LBC(:,:) = 0.0

      ENDIF ! on L_dust_div6_lbc, L_dust_div6

      ! ----------------------------------------------------------------
      ! Read SO2 lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_so2_lbc) THEN   ! so2 lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SO2_LBC,                              &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SO2_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SO2) THEN

        ! so2 prognostic is active but so2 lbcs are not in input file
        ! set lbc array to zero
        SO2_LBC(:,:) = 0.0

      ENDIF ! on L_so2_lbc, L_so2

      ! ----------------------------------------------------------------
      ! Read DMS lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_dms_lbc) THEN   ! dms lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,DMS_LBC,                              &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading DMS_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_DMS) THEN

        ! dms prognostic is active but dms lbcs are not in input file
        ! set lbc array to zero
        DMS_LBC(:,:) = 0.0

      ENDIF ! on L_dms_lbc, L_dms

      ! ----------------------------------------------------------------
      ! Read SO4_AITKEN lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_so4_aitken_lbc) THEN   ! so4_aitken lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SO4_AITKEN_LBC,                       &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SO4_AITKEN_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SO4_AITKEN) THEN

        ! so4_aitken prognostic is active but so4_aitken lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SO4_AITKEN_LBC(:,:) = 0.0

      ENDIF ! on L_so4_aitken_lbc, L_so4_aitken

      ! ----------------------------------------------------------------
      ! Read SO4_ACCU lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_so4_accu_lbc) THEN   ! so4_accu lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SO4_ACCU_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SO4_ACCU_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SO4_ACCU) THEN

        ! so4_accu prognostic is active but so4_accu lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SO4_ACCU_LBC(:,:) = 0.0

      ENDIF ! on L_so4_accu_lbc, L_so4_accu

      ! ----------------------------------------------------------------
      ! Read SO4_DISS lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_so4_diss_lbc) THEN   ! so4_diss lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SO4_DISS_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SO4_DISS_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SO4_DISS) THEN

        ! so4_diss prognostic is active but so4_diss lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SO4_DISS_LBC(:,:) = 0.0

      ENDIF ! on L_so4_diss_lbc, L_so4_diss

      ! ----------------------------------------------------------------
      ! Read NH3 lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_nh3_lbc) THEN   ! nh3 lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,NH3_LBC,                              &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading NH3_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_NH3) THEN

        ! nh3 prognostic is active but nh3 lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        NH3_LBC(:,:) = 0.0

      ENDIF ! on L_nh3_lbc, L_nh3

      ! ----------------------------------------------------------------
      ! Read soot lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_soot_new_lbc) THEN   ! soot_new lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SOOT_NEW_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SOOT_NEW_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SOOT_NEW) THEN

        ! soot_new prognostic is active but soot_new lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SOOT_NEW_LBC(:,:) = 0.0

      ENDIF ! on L_soot_new, L_soot_new

      !-----------------------------------------------------------------

      IF (L_soot_agd_lbc) THEN   ! soot_agd lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SOOT_AGD_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SOOT_AGD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SOOT_AGD) THEN

        ! soot_agd prognostic is active but soot_agd lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SOOT_AGD_LBC(:,:) = 0.0

      ENDIF ! on L_soot_agd, L_soot_agd

      !-----------------------------------------------------------------
 
      IF (L_soot_cld_lbc) THEN   ! soot_cld lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,SOOT_CLD_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),                &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading SOOT_CLD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_SOOT_CLD) THEN

        ! soot_cld prognostic is active but soot_cld lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        SOOT_CLD_LBC(:,:) = 0.0

      ENDIF ! on L_soot_cld, L_soot_cld

      ! ----------------------------------------------------------------
      ! Read biomass lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_bmass_new_lbc) THEN   ! bmass_new lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,BMASS_NEW_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading BMASS_NEW_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_BMASS_NEW) THEN

        ! bmass_new prognostic is active but bmass_new lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        BMASS_NEW_LBC(:,:) = 0.0

      ENDIF ! on L_bmass_new, L_bmass_new

      !-----------------------------------------------------------------

      IF (L_bmass_agd_lbc) THEN   ! bmass_agd lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,BMASS_AGD_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading BMASS_AGD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_BMASS_AGD) THEN

        ! bmass_agd prognostic is active but bmass_agd lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        BMASS_AGD_LBC(:,:) = 0.0

      ENDIF ! on L_bmass_agd, L_bmass_agd

      !-----------------------------------------------------------------
 
      IF (L_bmass_cld_lbc) THEN   ! bmass_cld lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,BMASS_CLD_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading BMASS_CLD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_BMASS_CLD) THEN

        ! bmass_cld prognostic is active but bmass_cld lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        BMASS_CLD_LBC(:,:) = 0.0

      ENDIF ! on L_bmass_cld, L_bmass_cld

      ! ----------------------------------------------------------------
      ! Read fossil fuel aerosol lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_ocff_new_lbc) THEN   ! ocff_new lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,OCFF_NEW_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading OCFF_NEW_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_OCFF_NEW) THEN

        ! ocff_new prognostic is active but ocff_new lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        OCFF_NEW_LBC(:,:) = 0.0

      ENDIF ! on L_ocff_new, L_ocff_new

      !-----------------------------------------------------------------

      IF (L_ocff_agd_lbc) THEN   ! ocff_agd lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,OCFF_AGD_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading OCFF_AGD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_OCFF_AGD) THEN

        ! ocff_agd prognostic is active but ocff_agd lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        OCFF_AGD_LBC(:,:) = 0.0

      ENDIF ! on L_ocff_agd, L_ocff_agd

      !-----------------------------------------------------------------
 
      IF (L_ocff_cld_lbc) THEN   ! ocff_cld lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,OCFF_CLD_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading OCFF_CLD_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_OCFF_CLD) THEN

        ! ocff_cld prognostic is active but ocff_cld lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        OCFF_CLD_LBC(:,:) = 0.0

      ENDIF ! on L_ocff_cld, L_ocff_cld

      ! ----------------------------------------------------------------
      ! Read nitrate aerosol lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_nitr_acc_lbc) THEN   ! nitr_acc lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,NITR_ACC_LBC,                         &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading NITR_ACC_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_NITR_ACC) THEN

        ! nitr_acc prognostic is active but nitr_acc lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        NITR_ACC_LBC(:,:) = 0.0

      ENDIF ! on L_nitr_acc_lbc, L_nitr_acc

      !-----------------------------------------------------------------

      IF (L_nitr_diss_lbc) THEN   ! nitr_diss lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
                      LEN1_LOOKUP,NITR_DISS_LBC,                        &
                      global_LENRIM(fld_type_p,halo_type_single)*       &
                      (tdims_s%k_end - tdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading NITR_DISS_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_NITR_DISS) THEN

        ! nitr_diss prognostic is active but nitr_diss lbcs 
        ! are not in input file.
        ! Set lbc array to zero
        NITR_DISS_LBC(:,:) = 0.0

      ENDIF ! on L_nitr_diss_lbc, L_nitr_diss

      ! ----------------------------------------------------------------
      ! Read tracer lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      DO tracer=1,TR_LBC_VARS

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num+tracer),               &
                      LEN1_LOOKUP, TRACER_LBCS(1,1,tracer),             &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (trdims_s%k_end-trdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading Tracer LBC ',   &
     &               tracer
          GOTO 9999
        ENDIF

      ENDDO  ! tracer

      !Increment number of lbcs
      lbc_num = lbc_num + TR_LBC_VARS

      ! ----------------------------------------------------------------
      ! Read UKCA tracer lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      DO tracer=1,TR_LBC_UKCA

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num+tracer),               &
                      LEN1_LOOKUP, UKCA_TRACER_LBCS(1,1,tracer),        &
                      global_LENRIM(fld_type_p,halo_type_extended)*     &
                      (trdims_s%k_end-trdims_s%k_start+1),              &
                      FIXHD,                                            &
                      ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading UKCA            &
     &                                  Tracer LBC ',                   &
     &               tracer
          GOTO 9999
        ENDIF

      ENDDO  ! tracer


 9999 CONTINUE

      IF (l_endgame) THEN
! DEPENDS ON: convert_lbcs
        CALL convert_lbcs(theta_lbc, q_lbc, qcl_lbc, qcf_lbc,             &
                           qcf2_lbc, qrain_lbc, qgraup_lbc,               &
                           l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,         &
                           LENRIM(fld_type_p,halo_type_extended) )
      END IF

      IF (lhook) CALL dr_hook('READ_ATMOS_LBCS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READ_ATMOS_LBCS
