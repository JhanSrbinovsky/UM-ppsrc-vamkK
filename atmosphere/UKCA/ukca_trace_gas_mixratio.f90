! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to hold species mass mixing ratios from atmos_physics1 (etc) 
!  which will be used to prescribe the surface concentrations
!  in UKCA when the logical L_ukca_prescribech4 is set to true, and when
!  lower boundary conditions are required.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_TRACE_GAS_MIXRATIO

      IMPLICIT NONE
      SAVE
      PUBLIC
    
      REAL :: um_ch4_for_ukca
      REAL :: um_co2_for_ukca
      REAL :: um_n2o_for_ukca
      REAL :: um_o2_for_ukca
      REAL :: um_cfc11_for_ukca
      REAL :: um_cfc12_for_ukca
      REAL :: um_cfc113_for_ukca
      REAL :: um_hcfc22_for_ukca
      REAL :: um_hfc125_for_ukca
      REAL :: um_hfc134a_for_ukca
      REAL :: um_mebr_for_ukca
      REAL :: um_mecl_for_ukca
      REAL :: um_ch2br2_for_ukca
      REAL :: um_h2_for_ukca
      REAL :: um_n2_for_ukca
      REAL :: um_cfc114_for_ukca
      REAL :: um_cfc115_for_ukca
      REAL :: um_ccl4_for_ukca
      REAL :: um_meccl3_for_ukca
      REAL :: um_hcfc141b_for_ukca
      REAL :: um_hcfc142b_for_ukca
      REAL :: um_h1211_for_ukca
      REAL :: um_h1202_for_ukca
      REAL :: um_h1301_for_ukca
      REAL :: um_h2402_for_ukca
      REAL :: um_cos_for_ukca

! These to be initialised from UMUI values
      REAL :: MeBrMMR
      REAL :: MeClMMR
      REAL :: CH2Br2MMR
      REAL :: H2MMR
      REAL :: N2MMR
      REAL :: CFC114MMR
      REAL :: CFC115MMR
      REAL :: CCl4MMR
      REAL :: MeCCl3MMR
      REAL :: HCFC141bMMR
      REAL :: HCFC142bMMR
      REAL :: H1211MMR
      REAL :: H1202MMR
      REAL :: H1301MMR
      REAL :: H2402MMR
      REAL :: COSMMR
        
      INTERFACE UKCA_SET_TRACE_GAS_MIXRATIO
        MODULE PROCEDURE UKCA_SET_TRACE_GAS_MIXRATIO
      END INTERFACE UKCA_SET_TRACE_GAS_MIXRATIO

      CONTAINS

! ######################################################################
! Description:
!  Subroutine to copy species mass mixing ratios from atmos_physics1 
!  and UKCA UMUI settings to module for use later in UKCA as lower
!  boundary conditions or as species concentrations (CO2, O2, H2 or N2)
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Called from Atmos_Physics1
!
! Method:
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_SET_TRACE_GAS_MIXRATIO(                           &
           ch4_mix_ratio, co2_mmr, n2o_mix_ratio,                       &
           o2mmr, cfc11_mix_ratio, cfc12_mix_ratio,                     &
           c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr)

      USE ukca_option_mod,      ONLY: ukca_MeBrMMR, ukca_MeClMMR,       &
                                      ukca_CH2Br2MMR, ukca_H2MMR,       &
                                      ukca_N2MMR, ukca_CFC115MMR,       &
                                      ukca_CCl4MMR, ukca_MeCCl3MMR,     &
                                      ukca_HCFC141bMMR, ukca_H1211MMR,  &
                                      ukca_HCFC142bMMR, ukca_H1202MMR,  &
                                      ukca_H1301MMR, ukca_H2402MMR,     &
                                      ukca_COSMMR, L_ukca_strat,        &
                                      L_ukca_strattrop, L_ukca_stratcfc,&
                                      L_ukca_achem
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE ereport_mod,          ONLY: ereport
      USE Control_Max_Sizes
      IMPLICIT NONE

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!  Mixing ratios from UMUI LW radiation settings  - these may be time varying
      REAL, INTENT(IN) :: ch4_mix_ratio
      REAL, INTENT(IN) :: co2_mmr
      REAL, INTENT(IN) :: n2o_mix_ratio
      REAL, INTENT(IN) :: o2mmr
      REAL, INTENT(IN) :: cfc11_mix_ratio
      REAL, INTENT(IN) :: cfc12_mix_ratio
      REAL, INTENT(IN) :: c113mmr
      REAL, INTENT(IN) :: c114mmr
      REAL, INTENT(IN) :: hcfc22mmr
      REAL, INTENT(IN) :: hfc125mmr
      REAL, INTENT(IN) :: hfc134ammr

      LOGICAL, SAVE :: L_first=.TRUE.
      CHARACTER(LEN=80) :: cmessage
      INTEGER :: ierr

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_SET_TRACE_GAS_MIXRATIO',zhook_in,  &
                               zhook_handle)

! Pick up UMUI values for UKCA gases,  
!  currently only need to do on first timestep
!  any missing data will be flagged below
      IF (L_first) THEN
         MeBrMMR     = ukca_MeBrMMR
         MeClMMR     = ukca_MeClMMR
         CH2Br2MMR   = ukca_CH2Br2MMR
         H2MMR       = ukca_H2MMR
         N2MMR       = ukca_N2MMR
         CFC115MMR   = ukca_CFC115MMR
         CCl4MMR     = ukca_CCl4MMR
         MeCCl3MMR   = ukca_MeCCl3MMR
         HCFC141bMMR = ukca_HCFC141bMMR
         HCFC142bMMR = ukca_HCFC142bMMR
         H1211MMR    = ukca_H1211MMR
         H1202MMR    = ukca_H1202MMR
         H1301MMR    = ukca_H1301MMR
         H2402MMR    = ukca_H2402MMR
         COSMMR      = ukca_COSMMR
      END IF

      IF (ch4_mix_ratio < 0.0) THEN
        IF (L_first) THEN
          cmessage='Missing value for trace_gas_mix_ratio CH4'
          WRITE(6,*) 'Set value for CH4 in UMUI panel -> ',             &
                     'setting to default pre-industrial level'
          ierr = -1
          CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
        END IF
        um_ch4_for_ukca = 4.360e-07
      ELSE
        um_ch4_for_ukca = ch4_mix_ratio
      END IF

      IF (co2_mmr < 0.0) THEN
        IF (L_first) THEN
          cmessage='Missing value for trace_gas_mix_ratio CO2'
          WRITE(6,*) 'Set value for CO2 in UMUI panel -> ',             &
                     'setting to default pre-industrial level'
          ierr = -2
          CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
        END IF
        um_co2_for_ukca = 4.31980e-04
      ELSE
        um_co2_for_ukca = co2_mmr
      END IF

      IF (H2MMR < 0.0) THEN
        IF (L_first) THEN
          cmessage='Missing value for trace_gas_mix_ratio H2'
          WRITE(6,*) 'Set value for H2 in UMUI hand-edit -> ',          &
                     'setting to default pre-industrial level'
          ierr = -14
          CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
        END IF
        um_h2_for_ukca = 3.452e-8 ! as in cinit
      ELSE
        um_h2_for_ukca = H2MMR
      END IF

! These species are only set for stratospheric chemistry schemes
      IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc) THEN

        IF (n2o_mix_ratio < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio N2O'
            WRITE(6,*) 'Set value for N2O in UMUI panel -> ',           &
                       'setting to default pre-industrial level'
            ierr = -3
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_n2o_for_ukca = 4.180e-07
        ELSE
          um_n2o_for_ukca = n2o_mix_ratio
        END IF

        IF (o2mmr < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio O2'
            WRITE(6,*) 'Set value for O2 in UMUI panel -> ',            &
                       'setting to default pre-industrial level'
            ierr = -4
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_o2_for_ukca = 0.23144 ! from cinit
        ELSE
          um_o2_for_ukca = o2mmr
        END IF


        IF (cfc11_mix_ratio  < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CFC11'
            WRITE(6,*) 'Set value for CFC11 in UMUI panel -> ',         &
                       'setting to default pre-industrial level'
            ierr = -5
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cfc11_for_ukca = 0.0
        ELSE
          um_cfc11_for_ukca = cfc11_mix_ratio
        END IF

        IF (cfc12_mix_ratio < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CFC12'
            WRITE(6,*) 'Set value for CFC12 in UMUI panel -> ',           &
                       'setting to pre-industrial level'
            ierr = -6
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cfc12_for_ukca = 0.0
        ELSE
          um_cfc12_for_ukca = cfc12_mix_ratio
        END IF

        IF (c113mmr < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CFC113'
            WRITE(6,*) 'Set value for CFC113 in UMUI panel -> ',        &
                       'setting to default pre-industrial level'
            ierr = -7
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cfc113_for_ukca = 0.0
        ELSE
          um_cfc113_for_ukca = c113mmr
        END IF

        IF (hcfc22mmr < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio HCFC22'
            WRITE(6,*) 'Set value for HCFC22 in UMUI panel -> ',        &
                       'setting to default pre-industrial level'
            ierr = -8
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_hcfc22_for_ukca = 0.0
        ELSE
          um_hcfc22_for_ukca = hcfc22mmr
        END IF

        IF (hfc125mmr < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio HFC125'
            WRITE(6,*) 'Set value for HFC125 in UMUI panel -> ',        &
                       'setting to default pre-industrial level'
            ierr = -9
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_hfc125_for_ukca = 0.0
        ELSE
          um_hfc125_for_ukca = hfc125mmr
        END IF

        IF (hfc134ammr < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio HFC134A'
            WRITE(6,*) 'Set value for HFC134A in UMUI panel -> ',       &
                       'setting to default pre-industrial level'
            ierr = -10
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_hfc134a_for_ukca = 0.0
        ELSE
          um_hfc134a_for_ukca = hfc134ammr
        END IF

        IF (MeBrMMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio MeBr'
            WRITE(6,*) 'Set value for MeBr in UMUI -> ',                &
                       'setting to default pre-industrial level'
            ierr = -11
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_mebr_for_ukca = 1.64025e-11 ! 5pptv
        ELSE
          um_mebr_for_ukca = MeBrMMR
        END IF

        IF (MeClMMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio MeCl'
            WRITE(6,*) 'Set value for MeCl in UMUI -> ',                &
                       'setting to default pre-industrial level'
            ierr = -12
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_mecl_for_ukca = 8.36736e-10 ! 480pptv
        ELSE
          um_mecl_for_ukca = MeClMMR
        END IF

        IF (CH2Br2MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CH2Br2'
            WRITE(6,*) 'Set value for CH2Br2 in UMUI -> ',              &
                       'setting to default pre-industrial level'
            ierr = -13
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_ch2br2_for_ukca = 18.0186e-12 ! 3pptv
        ELSE
          um_ch2br2_for_ukca = CH2Br2MMR
        END IF

        IF (N2MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio N2'
            WRITE(6,*) 'Set value for N2 in UMUI -> ',                  &
                       'setting to default pre-industrial level'
            ierr = -15
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_n2_for_ukca = 0.754682 ! as in cinit
        ELSE
          um_n2_for_ukca = N2MMR
        END IF

        IF (C114MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CFC114'
            WRITE(6,*) 'Set value for CFC-114 in UMUI -> ',             &
                       'setting to default pre-industrial level'
            ierr = -16
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cfc114_for_ukca = 0.0
        ELSE
          um_cfc114_for_ukca = C114MMR
        END IF

        IF (CFC115MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CFC115'
            WRITE(6,*) 'Set value for CFC-115 in UMUI -> ',             &
                       'setting to default pre-industrial level'
            ierr = -17
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cfc115_for_ukca = 0.0
        ELSE
          um_cfc115_for_ukca = CFC115MMR
        END IF

        IF (CCl4MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio CCl4'
            WRITE(6,*) 'Set value for CCl4 in UMUI -> ',                &
                       'setting to default pre-industrial level'
            ierr = -18
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_ccl4_for_ukca = 0.0
        ELSE
          um_ccl4_for_ukca = CCl4MMR
        END IF

        IF (MeCCl3MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio MeCCl3'
            WRITE(6,*) 'Set value for MeCCl3 in UMUI -> ',              &
                       'setting to default pre-industrial level'
            ierr = -18
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_meccl3_for_ukca = 0.0
        ELSE
          um_meccl3_for_ukca = MeCCl3MMR
        END IF

        IF (HCFC141bMMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio HCFC141b'
            WRITE(6,*) 'Set value for HCFC141b in UMUI -> ',            &
                       'setting to default pre-industrial level'
            ierr = -19
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_hcfc141b_for_ukca = 0.0
        ELSE
          um_hcfc141b_for_ukca = HCFC141bMMR
        END IF

        IF (HCFC142bMMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio HCFC142b'
            WRITE(6,*) 'Set value for HCFC142b in UMUI -> ',            &
                       'setting to default pre-industrial level'
            ierr = -20
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_hcfc142b_for_ukca = 0.0
        ELSE
          um_hcfc142b_for_ukca = HCFC142bMMR
        END IF

        IF (H1211MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio H1211'
            WRITE(6,*) 'Set value for H1211 in UMUI hand-edit -> ',     &
                       'setting to default pre-industrial level'
            ierr = -21
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_h1211_for_ukca = 0.0
        ELSE
          um_h1211_for_ukca = H1211MMR
        END IF

        IF (H1202MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio H1202'
            WRITE(6,*) 'Set value for H1202 in UMUI -> ',       &
                       'setting to default pre-industrial level'
            ierr = -22
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_h1202_for_ukca = 0.0
        ELSE
          um_h1202_for_ukca = H1202MMR
        END IF

        IF (H1301MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio H1301'
            WRITE(6,*) 'Set value for H1301 in UMUI -> ',               &
                       'setting to default pre-industrial level'
            ierr = -23
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_h1301_for_ukca = 0.0
        ELSE
          um_h1301_for_ukca = H1301MMR
        END IF

        IF (H2402MMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio H2402'
            WRITE(6,*) 'Set value for H2402 in UMUI -> ',               &
                       'setting to default pre-industrial level'
            ierr = -24
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_h2402_for_ukca = 0.0
        ELSE
          um_h2402_for_ukca = H2402MMR
        END IF

        IF (COSMMR < 0.0) THEN
          IF (L_first) THEN
            cmessage='Missing value for trace_gas_mix_ratio COS  '
            WRITE(6,'(A,A)') 'Set value for COS in UMUI  -> ',                &
                       'setting to default value of 520 pptm'
            ierr = -25
            CALL EREPORT('ukca_set_trace_gas_mixratio',ierr,cmessage)
          END IF
          um_cos_for_ukca = 520.0e-12 
           ! see Montzka et al, JGR, 109, D22302, (2004)
        ELSE
          um_cos_for_ukca = COSMMR
        END IF

      END IF    ! L_ukca_strat etc.


      L_first = .FALSE.

      
      IF (lhook) CALL dr_hook('UKCA_SET_TRACE_GAS_MIXRATIO',zhook_out,  &
                               zhook_handle)
      RETURN
      END SUBROUTINE UKCA_SET_TRACE_GAS_MIXRATIO

      END MODULE ukca_trace_gas_mixratio
