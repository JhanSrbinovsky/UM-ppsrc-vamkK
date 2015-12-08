! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_COPY_LBC_TEND

      SUBROUTINE IDL_COPY_LBC_TEND(                                     &
     &  LENRIM,                                                         &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, l_pc2,                   &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
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
     &  W_ADV_LBC, W_ADV_LBC_TEND                                       &
     &  )

!
! Purpose : Copies LBC array to LBC_TEND array
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

! Parameters required for argument declarations

! Arguments:

      INTEGER                                                           &
     &  LENRIM(Nfld_max,NHalo_Max)                                      &
                                    ! IN : Size of a level of LBC
     &, MODEL_LEVELS                                                    &
                                    ! IN : Number of model levels
     &, WET_LEVELS                  ! IN : Number of wet model levels

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2                                                      &
                         ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                         ! true if rain active
     &, L_mcr_qgraup                                                    &
                         ! true if graupel active
     &, L_pc2    ! true if prognostic cloud fractions are active

       Logical, Intent(In) ::                                           &
     &  L_mcr_qcf2_lbc                                                  &
                          ! true if second cloud ice lbcs active
     &, L_mcr_qrain_lbc                                                 &
                          ! true if rain lbcs active
     &, L_mcr_qgraup_lbc                                                &
                          ! true if graupel lbcs active
     &, L_pc2_lbc         ! true if prognostic cloud frac. lbcs active

      Real, Intent (InOut) ::                                           &
     &  U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! IN/OUT : U LBC
     &, U_LBC_TEND(LENRIM(fld_type_u,halo_type_extended),               &
     &             MODEL_LEVELS)                                        &
                                    ! IN : U LBC tendency
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! IN/OUT : V LBC
     &, V_LBC_TEND(LENRIM(fld_type_v,halo_type_extended),               &
     &             MODEL_LEVELS)                                        &
                                    ! IN : V LBC tendency
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        0:MODEL_LEVELS)                                           &
                                    ! IN/OUT : V LBC
     &, W_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),               &
     &             0:MODEL_LEVELS)                                      &
                                    ! IN : V LBC tendency
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          MODEL_LEVELS)                                           &
                                    ! IN/OUT : Rho LBC
     &, RHO_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               MODEL_LEVELS)                                      &
                                    ! IN : Rho LBC tendency
     &, THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : Theta LBC
     &, THETA_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : Theta LBC tendency
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        WET_LEVELS)                                               &
                                    ! IN/OUT : Q LBC
     &, Q_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),               &
     &             WET_LEVELS)                                          &
                                    ! IN : Q LBC tendency
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCL LBC
     &, QCL_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               WET_LEVELS)                                        &
                                    ! IN : QCL LBC tendency
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCL LBC
     &, QCF_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               WET_LEVELS)                                        &
                                    ! IN : QCL LBC tendency
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCF2 LBC
     &, QCF2_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),            &
     &               WET_LEVELS)                                        &
                                    ! IN : QCF2 LBC tendency
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QRAIN LBC
     &, QRAIN_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &               WET_LEVELS)                                        &
                                    ! IN : QRAIN LBC tendency
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QGRAUP LBC
     &, QGRAUP_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),          &
     &               WET_LEVELS)                                        &
                                    ! IN : QGRAUP LBC tendency
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_BULK LBC
     &, CF_BULK_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),         &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_BULK LBC tendency
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_LIQUID LBC
     &, CF_LIQUID_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_LIQUID LBC tendency
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_FROZEN LBC
     &, CF_FROZEN_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_FROZEN LBC tendency
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS+1)                                       &
                                    ! IN/OUT : Exner LBC
     &, EXNER_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 MODEL_LEVELS+1)                                  &
                                         ! IN : Exner LBC tendency
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : U_ADV LBC
     &, U_ADV_LBC_TEND(LENRIM(fld_type_u,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : U_ADV LBC tendency
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : V_ADV LBC
     &, V_ADV_LBC_TEND(LENRIM(fld_type_v,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : V_ADV LBC tendency
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            0:MODEL_LEVELS)                                       &
                                    ! IN/OUT : W LBC
     &, W_ADV_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 0:MODEL_LEVELS) ! IN : W LBC tendency

! Local variables

      INTEGER                                                           &
     &  k                                                               &
                  ! loop counter over levels
     &, i         ! loop counter over horizontal dimension

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ---------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_COPY_LBC_TEND',zhook_in,zhook_handle)

      DO k=1,MODEL_LEVELS
        DO i=1,LENRIM(fld_type_u,halo_type_extended)
          U_LBC_TEND(i,k)=U_LBC(i,k)
          U_ADV_LBC_TEND(i,k)=U_ADV_LBC(i,k)
        ENDDO ! i

        DO i=1,LENRIM(fld_type_v,halo_type_extended)
          V_LBC_TEND(i,k)=V_LBC(i,k)
          V_ADV_LBC_TEND(i,k)=V_ADV_LBC(i,k)
        ENDDO ! i

        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          RHO_LBC_TEND(i,k)=RHO_LBC(i,k)
          THETA_LBC_TEND(i,k)=THETA_LBC(i,k)
        ENDDO ! i
      ENDDO ! k

      DO k=1,WET_LEVELS
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          Q_LBC_TEND(i,k)=Q_LBC(i,k)
          QCL_LBC_TEND(i,k)=QCL_LBC(i,k)
          QCF_LBC_TEND(i,k)=QCF_LBC(i,k)
        ENDDO ! i
      ENDDO ! k

      If (L_mcr_qcf2 .AND. L_mcr_qcf2_lbc) Then  ! qcf2 active
        Do k=1,WET_LEVELS
          Do i=1,LENRIM(fld_type_p,halo_type_extended)
            QCF2_LBC_TEND(i,k) = QCF2_LBC(i,k)
          End Do ! i
        End Do ! k
      End If

      If (L_mcr_qrain .AND. L_mcr_qrain_lbc) Then  ! qrain active
        Do k=1,WET_LEVELS
          Do i=1,LENRIM(fld_type_p,halo_type_extended)
            QRAIN_LBC_TEND(i,k) = QRAIN_LBC(i,k)
          End Do ! i
        End Do ! k
      End If

      If (L_mcr_qgraup .AND. L_mcr_qgraup_lbc) Then  ! qgraup active
        Do k=1,WET_LEVELS
          Do i=1,LENRIM(fld_type_p,halo_type_extended)
            QGRAUP_LBC_TEND(i,k) = QGRAUP_LBC(i,k)
          End Do ! i
        End Do ! k
      End If

      If (L_pc2 .AND. L_pc2_lbc) Then  ! prog. cloud fraction active
        Do k=1,WET_LEVELS
          Do i=1,LENRIM(fld_type_p,halo_type_extended)
            CF_BULK_LBC_TEND(i,k) = CF_BULK_LBC(i,k)
            CF_LIQUID_LBC_TEND(i,k) = CF_LIQUID_LBC(i,k)
            CF_FROZEN_LBC_TEND(i,k) = CF_FROZEN_LBC(i,k)
          End Do ! i
        End Do ! k
      End If ! L_pc2

      DO k=0,MODEL_LEVELS
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          W_LBC_TEND(i,k)=W_LBC(i,k)
          W_ADV_LBC_TEND(i,k)=W_ADV_LBC(i,k)
        ENDDO ! i
      ENDDO ! k

      DO k=1,MODEL_LEVELS+1
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
           EXNER_LBC_TEND(i,k)=EXNER_LBC(i,k)
        ENDDO ! i
      ENDDO ! k

      IF (lhook) CALL dr_hook('IDL_COPY_LBC_TEND',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_COPY_LBC_TEND
