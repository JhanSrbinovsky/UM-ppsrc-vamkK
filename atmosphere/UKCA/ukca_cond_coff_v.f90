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
!    Subroutine which calculates the condensation coefficients
!    for condensable cpt vapour condensing on a particle with
!    radius RP. Includes a switch to either use Fuchs (1964)
!    or Fuchs-Sutugin (1971).
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
      SUBROUTINE UKCA_COND_COFF_V(NV,MASK,RP,TSQRT,AIRDM3,RHOA,         &
                   MMCG,SE,DMOL,IFUCHS,CC,SINKARR,PMID,T,DIFVOL,IDCMFP)
!----------------------------------------------------------------
!
! Purpose
! -------
! Subroutine which calculates the condensation coefficients
! for condensable cpt vapour condensing on a particle with
! radius RP. Includes a switch to either use Fuchs (1964)
! or Fuchs-Sutugin (1971).
!
! Inputs
! ------
! NV     : total number of values
! MASK   : Mask where to calculate values
! RP     : Radius of aerosol particle
! TSQRT  : Square-root of mid-level air temperature (K)
! AIRDM3 : Number density of air (m3)
! RHOA   : Air density (kg/m3)
! MMCG   : Molar mass of condensing gas (kg per mole)
! SE     : Sticking efficiency [taken as 0.3 from Raes et al 1992]
! DMOL   : Molecular diameter of condensable (m)
! IFUCHS : Switch for Fuchs (1964) or Fuchs-Sutugin (1971)
! PMID   : Centre level pressure (Pa)
! T      : Centre level temperature (K)
! DIFVOL : Diffusion volume for H2SO4 or SEC_ORG
! IDCMFP : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
!
! Outputs
! -------
! CC     : Condensation coeff for cpt onto ptcl (m^3s^-1)
!
! Local Variables
! ---------------
! VEL_CP   : Thermal velocity of condensable gas (ms-1)
! MFP_CP   : Mean free path of condensable gas (m)
! DCOFF_CP : Diffusion coefficient of condensable gas in air (m2s-1)
! KN       : Knudsen number
! FKN      : Correction factor for molecular effects
! AKN      : Corr fac for limitations in interfacial mass transport
! DENOM    : Denominator in Fuchs (1964) condensation coeff expression
! ZZ       : Ratio of molar masses of condensable gas and air
! TERM1,TERM2,.. : Terms in evaluation of condensation coefficient
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
! RA     : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RR     : Universal gas constant = 8.314 K mol^-1 K^-1
! PPI    : 3.1415927...
! AVC    : Avogadros constant (mol-1)
! MM_DA  : Molar mass of air (kg mol-1)
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,   ONLY: zboltz, ra, rr, ppi, avc, mm_da
      USE yomhook,          ONLY: lhook, dr_hook
      USE parkind1,         ONLY: jprb, jpim
      IMPLICIT NONE

! Subroutine interface:
      INTEGER, INTENT(IN) :: NV
      INTEGER, INTENT(IN) :: IFUCHS
      INTEGER, INTENT(IN) :: IDCMFP
      LOGICAL, INTENT(IN) :: MASK(NV)
      REAL, INTENT(IN)    :: RP(NV)
      REAL, INTENT(IN)    :: TSQRT(NV)
      REAL, INTENT(IN)    :: AIRDM3(NV)
      REAL, INTENT(IN)    :: RHOA(NV)
      REAL, INTENT(IN)    :: MMCG
      REAL, INTENT(IN)    :: SE
      REAL, INTENT(IN)    :: DMOL
      REAL, INTENT(IN)    :: PMID(NV)
      REAL, INTENT(IN)    :: T(NV)
      REAL, INTENT(IN)    :: DIFVOL
      REAL, INTENT(OUT)   :: CC(NV)
      REAL, INTENT(OUT)   :: SINKARR(NV)
!
! .. Local variables
      REAL    :: VEL_CP(NV)
      REAL    :: MFP_CP(NV)
      REAL    :: DCOFF_CP(NV)
      REAL    :: KN(NV)
      REAL    :: FKN(NV)
      REAL    :: AKN(NV)
      REAL    :: DENOM(NV)
      REAL    :: ZZ
      REAL    :: TERM1
      REAL    :: TERM2
      REAL    :: TERM3
      REAL    :: TERM4
      REAL    :: TERM5
      REAL    :: TERM6
      REAL    :: TERM7
      REAL    :: TERM8
      REAL    :: DAIR

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_COND_COFF_V',zhook_in,zhook_handle)

      TERM1=SQRT(8.0*RR/(PPI*MMCG))
! .. used in calcn of thermal velocity of condensable gas

      ZZ=MMCG/MM_DA
      TERM2=1.0/(PPI*SQRT(1.0+ZZ)*DMOL*DMOL)
! .. used in calcn of MFP of condensable gas (S&P, pg 457, eq 8.11)

      TERM3=(3.0/(8.0*AVC*DMOL*DMOL))
      TERM4=SQRT((RA*MM_DA*MM_DA/(2.0*PPI))*((MMCG+MM_DA)/MMCG))
      TERM5=TERM3*TERM4
! .. used in calcn of diffusion coefficient of condensable gas
!
      TERM6=4.0e6*PPI
! .. used in calculation of condensation coefficient
!
      TERM7=sqrt((1.0/(MM_DA*1000.0))+(1.0/(MMCG*1000.0)))
      DAIR=19.7 ! diffusion volume of air molecule (Fuller et al, Reid et al)
      TERM8=(DAIR**(1.0/3.0)+DIFVOL**(1.0/3.0))**2
! .. used in new calculation of diffusion coefficient

      CC(:)=0.0
      SINKARR(:)=0.0

! .. Calc. thermal velocity of condensable gas
      WHERE(MASK(:)) VEL_CP(:)=TERM1*TSQRT(:)

! .. Calculate diffusion coefficient of condensable gas
      IF(IDCMFP == 1) WHERE(MASK(:)) DCOFF_CP(:)=TERM5*TSQRT(:)/RHOA(:)
      IF(IDCMFP == 2) WHERE(MASK(:)) DCOFF_CP(:)=                       &
                (1.0E-7*(T(:)**1.75)*TERM7)/(PMID(:)/101325.0*TERM8)

!! n.b. IDCMFP=2 method is as below as coded in GLOMAP-bin:
!!        DCOFF_CP(:)=(1.0E-7*(T(:)**1.75)*sqrt((1.0/(MM_DA*1000.0))+     &
!!          (1.0/(MMCG*1000.0))))/                                        &
!!          (PMID(:)/101325.0*(DAIR**(1.0/3.0)+DIFVOL**(1.0/3.0))**2)

! .. Calc. mean free path of condensable gas
      IF(IDCMFP == 1) WHERE(MASK(:)) MFP_CP(:)=TERM2/AIRDM3(:) ! Gb v1
      IF(IDCMFP == 2) WHERE(MASK(:)) MFP_CP(:)=                         &
                                  3.0*DCOFF_CP(:)/VEL_CP(:) ! as Gb v1_1

      IF(IFUCHS == 1) THEN
!      If IFUCHS=1 use basic Fuchs (1964) expression
       WHERE(MASK(:))
        DENOM(:)=                                                       &
       (4.0*DCOFF_CP(:)/(SE*VEL_CP(:)*RP(:)))+(RP(:)/(RP(:)+MFP_CP(:)))
!       Calc. condensation coefficient
        CC(:)=TERM6*DCOFF_CP(:)*RP(:)/DENOM(:)
        SINKARR(:)=RP(:)*1.0E6/DENOM(:)
       ENDWHERE
      END IF

      IF(IFUCHS == 2) THEN
!      If IFUCHS=2 use basic Fuchs-Sutugin (1971) expression
       WHERE(MASK(:))
!       Calculate Knudsen number of condensable gas w.r.t. particle
        KN(:)=MFP_CP(:)/RP(:)
!       Calc. corr. factor for molecular effects
        FKN(:)=(1.0+KN(:))/(1.0+1.71*KN(:)+1.33*KN(:)*KN(:))
!       Calc. corr. factor for limitations in interfacial mass transport
        AKN(:)=1.0/(1.0+1.33*KN(:)*FKN(:)*(1.0/SE-1.0))
!       Calc. condensation coefficient
        CC(:)=TERM6*DCOFF_CP(:)*RP(:)*FKN(:)*AKN(:)
        SINKARR(:)=RP(:)*1.0E6*FKN(:)*AKN(:)
       ENDWHERE
      END IF

      IF (lhook) CALL dr_hook('UKCA_COND_COFF_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_COND_COFF_V
