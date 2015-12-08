! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
!    Parameterization based on
!    Method 1) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!           2) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_CALCNUCRATE(NBOX,DTZ,T,S,RH,AIRD,H2SO4,           &
       DELH2SO4_NUCL,SEC_ORG,JRATE,BLN_ON,IBLN,I_NUC_METHOD,            &
       HEIGHT,HTPBLG,S_COND_S)
!---------------------------------------------------------------
!
! Purpose
! -------
! Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
! Parameterization based on
! Method 1) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!        2) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
!       
! Parameters
! ----------
! None
!
! Inputs:
! ------
! NBOX    : Number of grid boxes
! DTZ     : Timestep for nucl/cond competition (s)
! T       : Mid-level temperature (K)
! S       : Specific humidity (kg/kg)
! RH      : Relative humidity (dimensionless 0-1)
! AIRD    : Number density of air (per cc)
! H2SO4   : H2SO4   number density (cm-3)
! SEC_ORG : SEC_ORG number density (cm-3)
! BLN_ON  : Switch for whether BLN is on (1) or off (0)
! IBLN    : Switch for whether BLN is ACT(1), KIN(2), PNAS(3), 
!           EUCAARI-kinetc(4), EUCAARI-org1(5) or EUCAARI-org2(6)
! I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamakki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
!  PNAS(IBLN=3),eucaari-kinetic(IBLN=4),eucaari-org1(IBLN=5),
!  eucaari-org2(IBLN=6)
!  note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate as in PNAS.
! HEIGHT  : Mid-level height of gridbox
! HTPBLG  : Height of boundary-layer in gridbox vertical-column
! S_COND_S :  Condensation sink
!
! Outputs:
! -------
! H2SO4      : H2SO4 number density (cm-3)  (INOUT)
! DELH2SO4_NUCL : Change in H2SO4 due to nucleation (molecules/DTZ)
! JRATE      : H2SO4 depletion rate by nucleation (/cc/s)
!
! Local Variables
! ---------------
! CORT      : Air temperature corrected to lie in range (K)
! CORRH     : Relative humidity corr. to lie in range (0-1)
! JPAN      : Pandis H2SO4/H2O nucleation rate (cm^-3 s^-1)
! JKUL      : Kulmala H2SO4/H2O nucleation rate (cm^-3 s^-1)
! JBLN      : BL nucleation rate (cm^-3 s^-1)
! H2SO4OLD2 : Old number density of H2SO4(g) (cm^-3)
! ACONS     : Constant in calculation of JKUL and JPAN
! BCONS     : Constant in calculation of JKUL and JPAN
! DELTA     : Parameter in calculation of JKUL
! A2        : Parameter in calculation of JKUL
! B2        : Parameter in calculation of JKUL
! C2        : Parameter in calculation of JKUL
! D2        : Parameter in calculation of JKUL
! E2        : Parameter in calculation of JKUL
! NCRIT     : H2SO4 no. dens. giving JKUL=1 cm^-3s^-1 (cm^-3)
! NSULF     : ln(H2SO4/NCRIT)
! PP        : partial pressure of H2SO4 (Pa)
! SVP       : saturation vapour pressure of H2SO4 (Pa)
! RAC       : relative acidity = PP/SVP
! NWP       : water vapour concentration (cm-3)
! XAL       : H2SO4 mole fraction
! THETA     : =chi in Kulmala et al (98) eq 20 = ln(JKUL)
! L1,L2     : Logicals for BLN methods and switches
! i_bhn_method : method for binary nucleation (Kulmala or Vehkamakki)
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ  : Boltzman Constant (kg m2 s-2 K-1 molec-1)
! NMOL    : Number of molecules per particle at nucleation
! CONC_EPS: Threshold for concentration (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NCP     : Number of possible aerosol components
! CP_SU   : Index of component in which sulfate is stored
!
! References
! ----------
! * Pandis et al (1994) "The relationship between
!   DMS flux and CCN concentration in remote marine
!   regions", JGR, vol. 99, no. D8, pp. 16,945-16,957.
!
! * Kulmala et al (1998) "Parameterizations for
!   sulfuric acid/water nucleation rates",
!   JGR, vol. 103, no. D7, pp. 8301-8307.
!
! * Vehkamaki et al (2002) "An improved parameterization for 
!   sulfuric acid-water nucleation rates for tropospheric
!   and stratospheric conditions!, JGR, vol 107, No D22
!   4622, doi:10.1029/2002JD002184
!
!
!
!----------------------------------------------------------------------
      USE UKCA_CONSTANTS,     ONLY: zboltz, nmol, conc_eps
      USE UKCA_MODE_SETUP,    ONLY: ncp, cp_su
      USE ukca_binapara_mod,  ONLY: ukca_binapara
      USE yomhook,            ONLY: lhook, dr_hook
      USE parkind1,           ONLY: jprb, jpim
      USE ereport_mod,        ONLY: ereport
      IMPLICIT NONE

!     Arguments
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: BLN_ON
      INTEGER, INTENT(IN) :: IBLN
      INTEGER, INTENT(IN) :: I_NUC_METHOD
      REAL, INTENT(IN)    :: SEC_ORG(NBOX)
      REAL, INTENT(IN)    :: DTZ
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: S(NBOX)
      REAL, INTENT(IN)    :: RH(NBOX)
      REAL, INTENT(IN)    :: AIRD(NBOX)
      REAL, INTENT(IN)    :: HEIGHT(NBOX)
      REAL, INTENT(IN)    :: HTPBLG(NBOX)
      REAL, INTENT(IN)    :: S_COND_S(NBOX)
      REAL, INTENT(INOUT) :: H2SO4(NBOX)
      REAL, INTENT(OUT)   :: JRATE(NBOX)
      REAL, INTENT(OUT)   :: DELH2SO4_NUCL(NBOX)

!     Local variables
      CHARACTER(LEN=72) :: cmessage     ! error message
      INTEGER, PARAMETER :: i_bhn_method_kulmala   = 1
      INTEGER, PARAMETER :: i_bhn_method_vekhamaki = 2
      INTEGER :: i_bhn_method
      INTEGER :: errcode                ! Variable passed to ereport
      INTEGER :: JL
      REAL :: CORT(NBOX)
      REAL :: CORRH(NBOX)
      REAL :: JPAN
      REAL :: H2SO4OLD2
      REAL :: ACONS
      REAL :: BCONS
      REAL :: DELTA
      REAL :: A2
      REAL :: B2
      REAL :: C2
      REAL :: D2
      REAL :: E2
      REAL :: NCRIT
      REAL :: NSULF
      REAL :: PP
      REAL :: SVP
      REAL :: NWP
      REAL :: XAL
      REAL :: RAC
      REAL :: THETA
      REAL :: JKUL
      REAL :: JBLN
      REAL :: JAPP_BLN
      REAL :: AFAC_ACT
      REAL :: AFAC_KIN
      REAL :: AFAC_PNA
      REAL :: AFAC_KIN2
      REAL :: AFAC_ORG
      REAL :: AFAC_ORG2
      REAL :: AFAC_KIN3
      REAL :: DPBLN
      LOGICAL :: L1,L2

      REAL :: JVEH(NBOX)
      REAL :: RC(NBOX)
      REAL :: JAPP(NBOX)
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_CALCNUCRATE',zhook_in,zhook_handle)

! .. Set nucleation method
!    1=Pandis(94)
!    2=Kulmala(98) with option to revert to activation (IBLN=1) 
!      or kinetic (IBLN=2) in BL, above revert to BLN if BLN_ON=1 
!      -- otherwise just use Kulmala in BL
!    3=use same rate at all levels either activation(IBLN=1), 
!     kinetic(IBLN=2),PNAS(IBLN=3),eucaari-kinetic(IBLN=4),eucaari-org1(IBLN=5),
!     eucaari-org2(IBLN=6) note that if I_NUC_METHOD=3 & IBLN=3 then also add 
!     on Kulmala rate as in PNAS.
!
!     i_bhn_method = i_bhn_method_kulmala   ! set for Kulmala    
      i_bhn_method = i_bhn_method_vekhamaki ! set for Vehkamakki
!
      AFAC_ACT=5.0E-7
      AFAC_KIN=4.0E-13
      AFAC_PNA=5.0E-13
      AFAC_KIN2=2.6E-13 
      AFAC_ORG=4.0E-13 
      AFAC_KIN3=8.2E-15 
      AFAC_ORG2=7.0E-14 
! 
! .. below defines the size used for the BL nucleation rate 
      IF((IBLN >= 1).AND.(IBLN <= 3)) DPBLN=1.5 ! Metzger  (2010,pnas) 
      IF((IBLN >= 4).AND.(IBLN <= 6)) DPBLN=2.0 ! Paasonen (2010,acpd)

! .. local copies of rh and t with bounds
      corrh(:) = rh(:)
      cort(:)  = t(:)
      WHERE (corrh < 0.1) corrh = 0.1
      WHERE (corrh > 0.9) corrh = 0.9
      WHERE (cort < 233.0) cort = 233.0
      WHERE (cort > 298.0) cort = 298.0

! Initialise JRATE and DELH2SO4_NUCL to zero
      JRATE(:)=0.0
      DELH2SO4_NUCL(:)=0.0

      IF (I_NUC_METHOD == 1) THEN
       DO JL=1,NBOX
!       Calculate H2SO4-H2O binary nucleation rate using Pandis(94)
        ACONS=10.0**(7.0 - (64.24 + 4.7*CORRH(JL)))
        BCONS=6.13+1.95*CORRH(JL)
        JPAN=ACONS*(H2SO4(JL)**BCONS)
        IF (JPAN > 1.0E-5) THEN
!        Calculate new H2S04 concentration
         H2SO4OLD2=H2SO4(JL)
         H2SO4(JL)=(H2SO4(JL)**(1.0 - BCONS) +                          &
               (BCONS - 1.0)*NMOL*ACONS*DTZ )                           &
                ** ((1.0/(1.0 - BCONS) ) )
         IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL) = H2SO4OLD2
         IF (H2SO4(JL) < 0.0) H2SO4(JL) = 0.0
!        Calculate change in sulfuric acid
         DELH2SO4_NUCL(JL)=DELH2SO4_NUCL(JL)+H2SO4OLD2-H2SO4(JL)
        END IF
       END DO
      END IF

      IF((I_NUC_METHOD == 2).OR.(I_NUC_METHOD == 3)) THEN

        IF(i_bhn_method == i_bhn_method_vekhamaki) THEN ! use Vehkamakki for BHN
         CALL UKCA_BINAPARA(NBOX,T,RH,H2SO4,JVEH,RC)
!        Vehkamakki et al (2002)
!        It returns nucleation rate JVEH at critical-cluster-radius RC
       END IF

       L2=((I_NUC_METHOD == 3).AND.(IBLN == 3))
! .. L2 is logical for having chosen BLN throughout the column but
! ..                   also have chosen PNAS approach which also includes
! ..                   Kulmala nucleation rate at all levels
       DO JL=1,NBOX

        L1=((I_NUC_METHOD == 2) .AND. (HEIGHT(JL) > HTPBLG(JL) .OR.     &
             BLN_ON == 0))
! .. L1 is logical for having chosen Kulmala  and either BLN is off
! ..                   (i.e. Kulmala over whole column) or 
! ..                   BLN is on and box is above the height of the BL
        IF (L1.OR.L2) THEN ! include Kulmala in JRATE and DELH2SO4_NUCL
!
         IF(i_bhn_method == i_bhn_method_kulmala) THEN 
!       Calculate H2SO4-H2O binary nucleation rate using Kulmala(98)
          IF (H2SO4(JL) > CONC_EPS) THEN
!          Calculate parametrised variables
           DELTA=1.0+(CORT(JL)-273.15)/273.15
           A2=25.1289-(4890.8/CORT(JL))-2.2479*DELTA*CORRH(JL)
           B2=(7643.4/CORT(JL))-1.9712*DELTA/CORRH(JL)
           C2=-1743.3/CORT(JL)
!          Calculate H2SO4 needed to get J=1cm-3 s-1
           NCRIT=EXP (-14.5125 + 0.1335*CORT(JL) -                      &
            10.5462*CORRH(JL) + 1958.4*CORRH(JL)/CORT(JL))
           NSULF=LOG(H2SO4(JL)/NCRIT)
!          Calculate saturation vapour pressure of H2SO4 (Pa)
           SVP=EXP( 27.78492066 - 10156.0/T(JL))
!          Calculate water vapour concentration (cm-3)
           NWP=AIRD(JL)*1.609*S(JL)
!          Calculate partial pressure of H2SO4 (Pa)
           PP=H2SO4(JL)*1.E6*ZBOLTZ*T(JL)
!          Calculate relative acidity
           RAC=PP/SVP
           IF (NWP > 0.0) THEN
            D2=1.2233 - (0.0154*RAC/(RAC+CORRH(JL)))-                   &
                   0.0415*LOG(NWP) + 0.0016*CORT(JL)
           ELSE
            D2=1.2233 - (0.0154*RAC/(RAC+CORRH(JL)))                    &
                   + 0.0016*CORT(JL)
           END IF
           E2 = 0.0102
           XAL=D2+E2*LOG(H2SO4(JL))
! .. NSULF is ln(H2SO4/NCRIT)
           THETA=A2*NSULF+B2*XAL+C2
! .. JKUL is nucleation rate in particles/cc/s
           JKUL=EXP(THETA)
!           write(6,*) 'JKUL (must be > 1.0e-3 to nucleate)=',JKUL
!          Calculate new H2SO4
           IF (JKUL > 1.E-3) THEN
!            write(6,*) 'JKUL>1.0E-3,here1'
            ACONS=(EXP(B2*D2+C2))*(1.0/NCRIT)**A2
            BCONS=A2 + B2*E2
!           Calculate new H2S04 concentrations
            H2SO4OLD2 = H2SO4(JL)
            H2SO4(JL)=( H2SO4(JL)**(1.0-BCONS) +                        &
                    (BCONS - 1.0)*NMOL*ACONS*DTZ)                       &
                    **((1.0 / (1.0 - BCONS)))
!           Calculate change in sulfuric acid
            IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL)=H2SO4OLD2
            IF (H2SO4(JL) < 0.0) H2SO4(JL)=0.0
            DELH2SO4_NUCL(JL)=DELH2SO4_NUCL(JL)+H2SO4OLD2-H2SO4(JL)
!           Calculate depletion rate of H2SO4 (molecules/cc/s)
            JRATE(JL)=JRATE(JL)+NMOL*ACONS*H2SO4OLD2**BCONS
           END IF ! JKUL > 1.0E-3
          END IF ! H2SO4 > CONC_EPS
         END IF ! i_bhn_method=i_bhn_method_kulmala (for Kulmala)
!
         IF(i_bhn_method == i_bhn_method_vekhamaki) THEN 
!         Calculate H2SO4-H2O binary nucleation rate using Vehkamakki

          IF ( (H2SO4(JL)>CONC_EPS).AND.(JVEH(JL)>0.0)                  &
             .AND.(S_COND_S(JL)>0.0).AND.(RC(JL)>0.0) ) THEN
!
           JAPP(JL)=JVEH(JL)*EXP(0.23*(1.0/3.0 - 1.0/(2*RC(JL)))*       &
                    S_COND_S(JL)/(1.E6*H2SO4(JL)*1.E-13))
!
!          Calculate new H2SO4 and particle number density
           IF (JAPP(JL)>1.0e-3) THEN
            JRATE(JL)=JRATE(JL)+JAPP(JL)
!
!           Calculate new H2S04 concentrations
            H2SO4OLD2 = H2SO4(JL)
            H2SO4(JL)= H2SO4(JL)-JRATE(JL)*NMOL*DTZ
!
!  Calculate change in sulfuric acid
            IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL)=H2SO4OLD2
            IF (H2SO4(JL) < 0.0) H2SO4(JL)=0.0
!
            DELH2SO4_NUCL(JL)=DELH2SO4_NUCL(JL)+H2SO4OLD2-H2SO4(JL)
!
           END IF ! if JAPP > 1.0e-3
          END IF ! if H2SO4 > 0 & JVEH > 0 & S_COND_S > 0 & RC > 0
!
         END IF ! if i_bhn_method=i_bhn_method_vekhamaki (for Vehkamakki)

        END IF ! if L1 or L2
!
        IF (.NOT.L1) THEN ! include BLN in JRATE and DELH2SO4_NUCL
! .. .NOT.L1 means either I_NUC_METHOD=3 or I_NUC_METHOD=2 and BLN_ON=1 
!     and in BL
         IF (H2SO4(JL) > CONC_EPS) THEN
          IF(IBLN == 1) JBLN=AFAC_ACT*H2SO4(JL)             ! activtn
          IF(IBLN == 2) JBLN=AFAC_KIN*H2SO4(JL)*H2SO4  (JL) ! kinetic
          IF(IBLN == 3) JBLN=AFAC_PNA*H2SO4(JL)*SEC_ORG(JL) ! as PNAS
          IF(IBLN == 4) JBLN=AFAC_KIN2*H2SO4(JL)*H2SO4  (JL) ! eucaari 1 
          IF(IBLN == 5) JBLN=AFAC_ORG *H2SO4(JL)*SEC_ORG(JL) ! eucaari 2 
          IF(IBLN == 6) JBLN=AFAC_KIN3*H2SO4(JL)*H2SO4  (JL)            &
                            +AFAC_ORG2*H2SO4(JL)*SEC_ORG(JL) ! eucaari 3 
! .. JBLN refers to nucleation rate at DPBLN diameter (nm) 
          JAPP_BLN=JBLN*EXP(0.23*(1.0/3.0 - 1.0/DPBLN)*                 &
                    S_COND_S(JL)/(1.E6*H2SO4(JL)*1E-13))
          IF(JAPP_BLN > 1.0e-3) THEN
           JRATE(JL)=JRATE(JL)+JAPP_BLN
           H2SO4OLD2 = H2SO4(JL)
           H2SO4(JL)=H2SO4(JL)-JAPP_BLN*NMOL*DTZ
!          Calculate change in sulfuric acid
           IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL)=H2SO4OLD2
           IF (H2SO4(JL) < 0.0) H2SO4(JL)=0.0
           DELH2SO4_NUCL(JL)=DELH2SO4_NUCL(JL)+H2SO4OLD2-H2SO4(JL)
          END IF ! JAPP>1.0e-3
         END IF ! H2SO4>CONC_EPS
        END IF ! METHOD=2 and in BL with BLN_ON=1 or METHOD=3

       END DO    ! Loop over NBOX
      END IF     ! METHOD 2 OR 3

!! note I_NUC_METHOD=1 not supported (cause model to stop here)
      IF ( ( I_NUC_METHOD < 2 .OR. I_NUC_METHOD > 3) .OR.               &
           ( IBLN < 1 .OR. IBLN > 6 ) ) THEN
       cmessage='I_NUC_METHOD (2 - 3) or IBLN (1 - 6) '// &
                 'out of valid range'
       errcode = 1
       CALL EREPORT('UKCA_CALCNUCRATE',errcode,cmessage)
      END IF

      IF (lhook) CALL dr_hook('UKCA_CALCNUCRATE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALCNUCRATE
