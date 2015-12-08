! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Fix the Lateral Boundary Conditions (LBCs) of LAM fields

      SUBROUTINE IDL_FIX_LAM_LBCs(                                      &
     &  ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,                 &
     &  HALO_I,HALO_J,AT_EXTREMITY,                                     &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, l_pc2,                   &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
     &  RIMWIDTH,RIMWEIGHTS,LENRIM,LBC_SIZE,LBC_START,                  &
     &  THETA_LBC, Q_LBC, QCL_LBC, QCF_LBC,                             &
     &  QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
     &  CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
     &  RHO_LBC, EXNER_LBC,                                             &
     &  U_LBC, V_LBC, W_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,           &
     &  theta_eh, rho_eh, exner_rho_levels_eh,                          &
     &  Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                               &
     &  CF_BULK, CF_LIQUID, CF_FROZEN, U_ADV, V_ADV, W_ADV              &
     &  )

! Purpose:
! Fix the LBCs from the input fields
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Input
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

! Parameters required for dimensioning some of the arguments

! Arguments

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                          ! IN : Length of a model row
     &, ROWS                                                            &
                          ! IN : Number of rows for theta,u fields
     &, N_ROWS                                                          &
                          ! IN : Number of rows for v fields
     &, MODEL_LEVELS                                                    &
                          ! IN : Number of model levels
     &, WET_LEVELS                                                      &
                    ! IN : Number of moist model levels
     &, HALO_I                                                          &
                          ! IN : Size of extended halo (EW direction)
     &, HALO_J                                                          &
                          ! IN : Size of extended halo (NS direction)
     &, RIMWIDTH                                                        &
                          ! IN : Size of boundary region
     &, LENRIM(Nfld_max,NHalo_max)                                      &
                          ! IN : Size of single level of LBC
     &, LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
                          ! IN : Size of a side of a LBC
     &, LBC_START(4,Nfld_max,NHalo_max)
                          ! IN : Start of a side in a LBC

      Logical, Intent (In) ::                                           &
     &  AT_EXTREMITY(4)                                                 &
                          ! IN : At an edge?
     &, L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup                                                    &
                          ! true if graupel active
     &, L_pc2                                                           &
                          ! true if prognostic cloud fraction active
     &, L_mcr_qcf2_lbc                                                  &
                          ! true if second cloud ice lbcs active
     &, L_mcr_qrain_lbc                                                 &
                          ! true if rain lbcs active
     &, L_mcr_qgraup_lbc                                                &
                          ! true if graupel lbcs active
     &, L_pc2_lbc         ! true if prognostic cloud frac. lbcs active

      REAL                                                              &
     &  RIMWEIGHTS(RIMWIDTH)  ! IN : weight to apply to LBC

      REAL                                                              &
             ! LBCs : Intent OUT
     &  THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        WET_LEVELS)                                               &
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
     &          WET_LEVELS)                                             &
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &          WET_LEVELS)                                             &
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &          WET_LEVELS)                                             &
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
     &          WET_LEVELS)                                             &
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          MODEL_LEVELS)                                           &
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS+1)                                       &
     &, U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        0:MODEL_LEVELS)                                           &
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            0:MODEL_LEVELS)

      ! Work arrays with extended halos (_eh)
      ! Needed so that external halo values in LAMS can be set correctly
      ! for the lateral boundary arrays.
      REAL, Intent (InOut) ::                                           &
     &  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           model_levels)                                          &
     &, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &           model_levels)                                          &
     &, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                 &
     &           1-halo_j:rows+halo_j, model_levels+1)

      Real, Intent (InOut) ::                                           &
     &  Q(1-HALO_I:ROW_LENGTH+HALO_I,                                   &
     &    1-HALO_J:ROWS+HALO_J,                                         &
     &    WET_LEVELS)                                                   &
     &, QCL(1-HALO_I:ROW_LENGTH+HALO_I,                                 &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QCF(1-HALO_I:ROW_LENGTH+HALO_I,                                 &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QCF2(1-HALO_I:ROW_LENGTH+HALO_I,                                &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QRAIN(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QGRAUP(1-HALO_I:ROW_LENGTH+HALO_I,                              &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_BULK(1-HALO_I:ROW_LENGTH+HALO_I,                             &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_LIQUID(1-HALO_I:ROW_LENGTH+HALO_I,                           &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_FROZEN(1-HALO_I:ROW_LENGTH+HALO_I,                           &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, U_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:ROWS+HALO_J,                                     &
     &        MODEL_LEVELS)                                             &
     &, V_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:N_ROWS+HALO_J,                                   &
     &        MODEL_LEVELS)                                             &
     &, W_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:ROWS+HALO_J,                                     &
     &        0:MODEL_LEVELS)

! Local variables

      LOGICAL                                                           &
     &  L_Do_Boundaries                                                 &
     &, L_Do_Halos
      INTEGER                                                           &
     &  N_RIMS_TO_DO                                                    &
     &, i, j, k    ! loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_FIX_LAM_LBCS',zhook_in,zhook_handle)

      L_Do_Boundaries=.TRUE.
      L_Do_Halos=.TRUE.
      N_RIMS_TO_DO=RIMWIDTH

! THETA
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_p,theta_eh, &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,THETA_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! Q
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,Q,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,Q_LBC,                                            &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! RHO
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_p,rho_eh,   &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,RHO_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! EXNER
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,        &
     &  exner_rho_levels_eh,                                            &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,EXNER_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! U
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_u,U_ADV,    &
     &  LENRIM(fld_type_u,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I, HALO_J,U_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! V
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,N_ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_v,V_ADV,  &
     &  LENRIM(fld_type_v,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I, HALO_J,V_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! U_ADV
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_u,          &
     &  U_ADV,LENRIM(fld_type_u,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I,HALO_J,U_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! V_ADV
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,N_ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_v,        &
     &  V_ADV,LENRIM(fld_type_v,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I,HALO_J,V_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! W
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,W_ADV,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I, HALO_J,W_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! W_ADV
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,        &
     &  W_ADV,LENRIM(fld_type_p,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,W_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
! QCL
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCL,        &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCL_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! QCF
! DEPENDS ON: idl_fix_lateral_boundaries
      CALL IDL_FIX_LATERAL_BOUNDARIES(                                  &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF,        &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCF_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! QCF2
      If (L_mcr_qcf2 .AND. L_mcr_qcf2_lbc) Then  ! qcf2 active
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF2,       &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCF2_LBC,                                         &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

! QRAIN
      If (L_mcr_qrain .AND. L_mcr_qrain_lbc) Then  ! qrain active
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QRAIN,      &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QRAIN_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

! QGRAUP
      If (L_mcr_qgraup .AND. L_mcr_qgraup_lbc) Then  ! qgraup active
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QGRAUP,     &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QGRAUP_LBC,                                       &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

      If (L_pc2 .AND. L_pc2_lbc) Then  ! prog. cloud fraction active

! CF_BULK
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_BULK,    &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_BULK_LBC,                                      &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! CF_LIQUID
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_LIQUID,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_LIQUID_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! CF_FROZEN
! DEPENDS ON: idl_fix_lateral_boundaries
       CALL IDL_FIX_LATERAL_BOUNDARIES(                                 &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_FROZEN,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_FROZEN_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

      End If  ! L_pc2

      IF (lhook) CALL dr_hook('IDL_FIX_LAM_LBCS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_FIX_LAM_LBCs
