! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Update the Lateral Boundaries of each end-of-iteration temp field
!  (e.h u_np1, etc) when the iterative semi-implicit semi-Lagrangian 
!  scheme is used and CycleNo >1. 
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: LBC Input

      SUBROUTINE UPDATE_LBC_ITERSL(                                     &
     &  ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,                 &
     &  OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,L_new_tdisc,               &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc,   &
     &  RIMWIDTH,RIMWEIGHTS,LENRIM,LBC_SIZE,LBC_START,                  &
     &  THETA_LBC, Q_LBC, QCL_LBC, QCF_LBC,                             &
     &  QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
     &  CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
     &  RHO_LBC, EXNER_LBC,                                             &
     &  U_LBC, V_LBC, W_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,           &
     &  THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                        &
     &  CF_BULK, CF_LIQUID, CF_FROZEN,                                  &
     &  RHO, EXNER, U, V, W, U_ADV, V_ADV, W_ADV                        &
     &  )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE
! Updates all the input fields with their LBCs
!
! Parameters required for dimensioning some of the arguments

! Arguments

      Integer, Intent (In) ::                                           &
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
     &, OFFX                                                            &
                          ! IN : Size of "single" halo (EW direction)
     &, OFFY                                                            &
                          ! IN : Size of "single" halo (NS direction)
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
     &, L_mcr_qcf2_lbc                                                  &
                          ! true if prognostic 2nd cloud ice in lbcs
     &, L_mcr_qrain_lbc                                                 &
                          ! true if prognostic rain in lbcs
     &, L_mcr_qgraup_lbc                                                &
                          ! true if prognostic graupel in lbcs
     &, L_pc2_lbc                                                       &
                          ! true if prognostic cloud fracs in lbcs
     &, L_new_tdisc


      Real, Intent(In) :: RIMWEIGHTS(RIMWIDTH)  
                          ! IN : weight to apply to LBC

      Real, Intent(In)::                                                &
     &  THETA_LBC(LENRIM(fld_type_p,halo_type_extended),MODEL_LEVELS)   &
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)         &
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)       &
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)       &
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)      &
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)     &
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)    &
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS)   &
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS) &
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),WET_LEVELS) &
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),MODEL_LEVELS)     &
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),MODEL_LEVELS+1) &
     &, U_LBC(LENRIM(fld_type_u,halo_type_extended),MODEL_LEVELS)       &
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),MODEL_LEVELS)       &
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),0:MODEL_LEVELS)     &
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),MODEL_LEVELS)   &
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),MODEL_LEVELS)   &
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),0:MODEL_LEVELS) 

      Real, Intent (InOut) ::                                           &
     &  THETA(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,MODEL_LEVELS)     &
     &, Q(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)           &
     &, QCL(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)         &
     &, QCF(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)         &
     &, QCF2(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)        &
     &, QRAIN(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)       &
     &, QGRAUP(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,WET_LEVELS)      &
     &, CF_BULK(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS+HALO_J,        &
     &          WET_LEVELS)                                             &
     &, CF_LIQUID(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS+HALO_J,      &
     &            WET_LEVELS)                                           &
     &, CF_FROZEN(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS+HALO_J,      &
     &            WET_LEVELS)                                           &
     &, RHO(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,MODEL_LEVELS)       &
     &, EXNER(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,MODEL_LEVELS+1)   &
     &, U(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,MODEL_LEVELS)         &
     &, V(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:N_ROWS+OFFY,MODEL_LEVELS)       &
     &, W(1-OFFX:ROW_LENGTH+OFFX,1-OFFY:ROWS+OFFY,0:MODEL_LEVELS)       &
     &, U_ADV(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS+HALO_J,          &
     &        MODEL_LEVELS)                                             &
     &, V_ADV(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:N_ROWS+HALO_J,        &
     &        MODEL_LEVELS)                                             &
     &, W_ADV(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS+HALO_J,          &
     &        0:MODEL_LEVELS)

! Local variables

      LOGICAL L_Do_Boundaries, L_Do_Halos
      INTEGER N_RIMS_TO_DO

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

       IF (lhook) CALL dr_hook('UPDATE_LBC_ITERSL',zhook_in,zhook_handle)

       L_Do_Boundaries=.TRUE.
       L_Do_Halos=.TRUE.
       N_RIMS_TO_DO=RIMWIDTH

       IF ( L_NEW_TDISC ) THEN

! RHO
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_p,RHO,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,RHO_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                     

! U
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_u,U,            &
     &  LENRIM(fld_type_u,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I, HALO_J,U_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                    

! V
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,N_ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_v,V,          &
     &  LENRIM(fld_type_v,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I, HALO_J,V_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                      
! W
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS+1,fld_type_p,W,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I, HALO_J,W_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

      END IF ! L_new_tdisc

! Cloud fractions
      IF (L_pc2_lbc) THEN  ! prognostic cloud fracs in lbcs

! CF_BULK
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_BULK,    &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_BULK_LBC,                                      &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! CF_LIQUID
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_LIQUID,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_LIQUID_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                     

! CF_FROZEN
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_FROZEN,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_FROZEN_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

      End If ! L_pc2_lbc

! U_ADV
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_u,          &
     &  U_ADV,LENRIM(fld_type_u,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I,HALO_J,U_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                     

! V_ADV
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,N_ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_v,        &
     &  V_ADV,LENRIM(fld_type_v,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I,HALO_J,V_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                     

! W_ADV
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,        &
     &  W_ADV,LENRIM(fld_type_p,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,W_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)                                     

       IF (lhook) CALL dr_hook('UPDATE_LBC_ITERSL',zhook_out,zhook_handle)
       RETURN
       END SUBROUTINE UPDATE_LBC_ITERSL
