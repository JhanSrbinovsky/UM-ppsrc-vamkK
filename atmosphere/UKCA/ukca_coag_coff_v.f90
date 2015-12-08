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
!    Calculate coagulation coefficient for modes coagulating
!    with radii RI, RJ, volumes VI, VJ.
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
      SUBROUTINE UKCA_COAG_COFF_V(NV,MASK,RI,RJ,VI,VJ,RHOI,RHOJ,        &
           MFPA,DVISC,T,PMID,KIJ,COAG_ON,ICOAG)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculate coagulation coefficient for modes coagulating
! with radii RI, RJ, volumes VI, VJ.
!
! 4 possible methods of calculation chosen via switch ICOAG.
!
! 1) Method as Jacobson, "Fundamentals of Atmospheric Modelling", p. 446
!    for Brownian coagulation between ptls in the transition regime
!
!      K_ij=4*PI*(r_i+r_j)*(D_i+D_j)/(TERM1+TERM2)
!
!    where r_i,r_j are (wet) particle radii
!          D_i,D_j are particle diffusion coefficients
!
!          D_i=ACunn_i*kT/(6*PI*r_i*DVISC)
!
!          ACunn_i is the Cunningham slip correction
!                = 1.0+Kn_i*(1.257+0.4*EXP{-1.1/Kn_i})
!
!          Kn_i is the Knudsen number = MFPA/r_i
!
!          MFPA=2*DVISC/(RHOA*vel_i)
!
!          TERM1=(r_i+r_j)/(r_i+r_j+sqrt{delta_i^2+delta_j^2}
!          TERM2=4*(D_i+D_j)/( sqrt{vel_i^2+vel_j^2}*(r_i+r_j) )
!
!          delta_i,delta_j are mean distances from sphere centre after
!                          travelling distance of particle mfp (mfpp).
!          delta_i=( (2*r_i+mfpp_i)^3 - (4*r_i^2 + mfpp_i^2)^{3/2} )
!
!          mfpp_i= mean free path of particle = 2*D_i/PI/vel_i
!
!          vel_i,vel_j are thermal speed of particles
!                vel_i = sqrt{8kT/PI/mass_i}
!
!          where k is Boltzmann's constant
!                mass_i is the mass of a particle = rho_i*vol_i
!
! 2) Method as in HAM/M7, see Stier et al (2005), Vignati et al (2004)
!
!      K_ij=16*pi*RMID*DMID/(4*DMID/VMID/RMID + RMID/(RMID+DELTAMID))
!
!    where DMID,VMID,DELTAMID are D_i, vel_i and delta_i
!          evaluated for a particle with mean radius RMID=(r_i+r_j)/2
!
!          n.b. I have been informed since then that although this
!          is what is described in the Stier et al (2005) and
!          Vignati et al (2004), that in the latest HAM module,
!          DMID,VMID and DELTAMID are now evaluated more completely as:
!
!          DMID={D(r_i)+D(r_j)}/2.0
!          VMID=sqrt{vel_i^2+vel_j^2}
!          DELTAMID=sqrt{delta_i^2+delta_j^2}/sqrt{2.0}
!
! 3) Method as in UM sulfur scheme as in UMDP 20
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                                      MFP*ACunn*((1/r_i)^2+(1/r_j)^2)}
!
!    where k=Boltzmann's constant
!          DVISC=dynamic viscosity of air
!          T=air temperature
!          MFP=MFPA=mean free path of air=6.6e-8*p0*T/(p*T0)
!                             [p0=1.01325e5 Pa,T0=293.15K]
!          ACunn=Cunningham slip correction assumed = 1.59 for all modes
!
! 4) Method as in 3) but with MFPP_i,MFPPj and ACunn_i,ACunn_j and MFPA
!                             calculated explicitly as in 1)
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                         (MFPP_i*ACunn_i/r_i^2)+(MFPP_j*ACunn_j/r_j^2)}
!
! Vector version !!!!!!!!!!!!!!!!!!!!!
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NV     : Number of values
! MASK   : Mask where to calculate values
! RI     : Geometric mean radius for mode I
! RJ     : Geometric mean radius for mode J
! VI     : Volume for ptcl with r=RI
! VJ     : Volume for ptcl with r=RJ
! T      : Air temperature (K)
! PMID   : Air pressure (Pa)
! RHOI   : Density of aerosol particle for mode I (kgm^-3)
! RHOJ   : Density of aerosol particle for mode J (kgm^-3)
! MFPA   : Mean free path of air (m)
! DVISC  : Dynamic viscosity of air (kg m^-1 s^-1)
! COAG_ON: Switch for turning coagulation off
! ICOAG  : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!
! Outputs
! -------
! KIJ    : Coagulation coefficient for I-J (cm^3 s^-1)
!
! Local Variables
! ---------------
! RMID   : Arithmetic mean of I and J mode mean radii (m)
! VMID   : Volume corresponding to radius RMID
! VEL    : Thermal velocity of ptcl of size RMID (m/s)
! DCOEF  : Diffusion coefficient of ptcl of size RMID (m^2 s^-1)
! MFPP   : Mean free path of aerosol ptcl of size RMID (m)
! VELI   : Thermal velocity of ptcl in mode I (m/s)
! DCOEFI : Diffusion coefficient of ptcl in mode I (m^2 s^-1)
! MFPPI  : Mean free path of aerosol ptcl in mode I (m)
! DELI   : Coeff. in calc. of coag kernel for mode I
! VELJ   : Thermal velocity of ptcl in mode J (m/s)
! DCOEFJ : Diffusion coefficient of ptcl in mode J (m^2 s^-1)
! MFPPJ  : Mean free path of aerosol ptcl in mode J (m)
! DELJ   : Coeff. in calc. of coag kernel for mode J
! KNI    : Knudsen number of particle in mode I
! CCI    : Cunningham correction factor for ptcl in mode I
! KNJ    : Knudsen number of particle in mode J
! CCJ    : Cunningham correction factor for ptcl in mode J
! TERMV1 : Vector term in calculation of KIJ
! TERMV2 : Vector term in calculation of KIJ
! TERMV3 : Vector term in calculation of KIJ
! TERMV4 : Vector term in calculation of KIJ
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! PPI    : 3.1415927
! ZBOLTZ : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,     ONLY: ppi, zboltz
      USE parkind1,           ONLY: jprb, jpim
      USE yomhook,            ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface:
      INTEGER, INTENT(IN) :: NV
      INTEGER, INTENT(IN) :: COAG_ON
      INTEGER, INTENT(IN) :: ICOAG
      REAL, INTENT(IN)    :: RI(NV)
      REAL, INTENT(IN)    :: RJ(NV)
      REAL, INTENT(IN)    :: VI(NV)
      REAL, INTENT(IN)    :: VJ(NV)
      REAL, INTENT(IN)    :: RHOI(NV)
      REAL, INTENT(IN)    :: RHOJ(NV)
      REAL, INTENT(IN)    :: MFPA(NV)
      REAL, INTENT(IN)    :: DVISC(NV)
      REAL, INTENT(IN)    :: T(NV)
      REAL, INTENT(IN)    :: PMID(NV)
      LOGICAL, INTENT(IN) :: MASK(NV)
      REAL, INTENT(OUT)   :: KIJ(NV)
!
! Local variables:
      REAL :: DPI(NV)
      REAL :: DPJ(NV)
      REAL :: VELI(NV)
      REAL :: DCOEFI(NV)
      REAL :: MFPPI(NV)
      REAL :: DELI(NV)
      REAL :: KNI(NV)
      REAL :: CCI(NV)
      REAL :: VELJ(NV)
      REAL :: DCOEFJ(NV)
      REAL :: MFPPJ(NV)
      REAL :: DELJ(NV)
      REAL :: KNJ(NV)
      REAL :: CCJ(NV)
      REAL :: VEL(NV)
      REAL :: DCOEF(NV)
      REAL :: MFPP(NV)
      REAL :: KN(NV)
      REAL :: CC(NV)
      REAL :: RTOT(NV)
      REAL :: DTOT(NV)
      REAL :: RMID(NV)
      REAL :: VMID(NV)
      REAL :: RHOMID(NV)
      REAL :: TERMV1(NV)
      REAL :: TERMV2(NV)
      REAL :: TERMV3(NV)
      REAL :: TERMV4(NV)
      REAL :: MFPA_CALC(NV)
      REAL :: TERM1
      REAL :: TERM2
      REAL :: TERM3
      REAL :: TERM4
      REAL :: TERM5

      REAL, PARAMETER :: UKCA_MFP_REF=6.6E-8
      REAL, PARAMETER :: UKCA_PREF_MFP=1.01325E5
      REAL, PARAMETER :: UKCA_TREF_MFP=293.15E0
      REAL :: UKCA_REF

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('UKCA_COAG_COFF_V',zhook_in,zhook_handle)
      KIJ(:) = 0.0
!
      IF(COAG_ON == 0) THEN
        IF (lhook) CALL dr_hook('UKCA_COAG_COFF_V',zhook_out,zhook_handle)
        RETURN
      END IF
!
      TERM1=8.0*ZBOLTZ/PPI
      TERM2=ZBOLTZ/(6.0*PPI)
      TERM3=8.0/PPI
      TERM4=4.0e6*PPI
      TERM5=16.0e6*PPI
!
!  .. For ICOAG=1, calculate KIJ following the full method
!  .. as GLOMAP does for each bin (see header above).
      IF(ICOAG == 1) THEN
       WHERE(MASK(:))

        DPI(:)=2.0*RI(:)
        VELI(:)=SQRT(TERM1*T(:)/(RHOI(:)*VI(:)))
        KNI(:)=MFPA(:)/RI(:)
        CCI(:)=1.0+KNI(:)*(1.257+0.4*EXP(-1.1/KNI(:)))
! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
        DCOEFI(:)=TERM2*CCI(:)*T(:)/(RI(:)*DVISC(:))
        MFPPI(:)=TERM3*DCOEFI(:)/VELI(:)
        DELI(:)=((DPI(:)+MFPPI(:))**3-                                  &
                 SQRT((DPI(:)*DPI(:)+MFPPI(:)*MFPPI(:))**3))/           &
                 (3.0*DPI(:)*MFPPI(:))-DPI(:)
!
        DPJ(:)=2.0*RJ(:)
        VELJ(:)=SQRT(TERM1*T(:)/(RHOJ(:)*VJ(:)))
        KNJ(:)=MFPA(:)/RJ(:)
        CCJ(:)=1.0+KNJ(:)*(1.257+0.4*EXP(-1.1/KNJ(:)))
! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
        DCOEFJ(:)=TERM2*CCJ(:)*T(:)/(RJ(:)*DVISC(:))
        MFPPJ(:)=TERM3*DCOEFJ(:)/VELJ(:)
        DELJ(:)=((DPJ(:)+MFPPJ(:))**3-                                  &
                 SQRT((DPJ(:)*DPJ(:)+MFPPJ(:)*MFPPJ(:))**3))/           &
                 (3.0*DPJ(:)*MFPPJ(:))-DPJ(:)
!
        RTOT(:)=RI(:)+RJ(:)
        DTOT(:)=DCOEFI(:)+DCOEFJ(:)
        TERMV1(:)=                                                      &
         RTOT(:)/(RTOT(:)+SQRT(DELI(:)*DELI(:)+DELJ(:)*DELJ(:)))
        TERMV2(:)=                                                      &
         4.0*DTOT(:)/SQRT(VELI(:)*VELI(:)+VELJ(:)*VELJ(:))/RTOT(:)
        KIJ(:)=TERM4*RTOT(:)*DTOT(:)/(TERMV1(:)+TERMV2(:))
       ENDWHERE
      END IF
!
!  .. For ICOAG=2, calculates KIJ following the M7 method as
!  .. described in Stier et al (2005) in ACP (see header above)
      IF(ICOAG == 2) THEN
       WHERE(MASK(:))
        RMID(:)=0.5*(RI(:)+RJ(:))
        VMID(:)=(PPI/0.75)*RMID(:)**3
        RHOMID(:)=0.5*(RHOI(:)+RHOJ(:))
        VEL(:)=SQRT(TERM1*T(:)/(RHOMID(:)*VMID(:)))
        KN(:)=MFPA(:)/RMID(:)
        CC(:)=1.0+KN(:)*(1.257+0.4*EXP(-1.1/KN(:)))
! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
        DCOEF(:)=TERM2*CC(:)*T(:)/(RMID(:)*DVISC(:))
        TERMV1(:)=4.0*DCOEF(:)/(VEL(:)*RMID(:))
        MFPP(:)=2.0*DCOEF(:)/(PPI*VEL(:))
        TERMV2(:)=RMID(:)/(RMID(:)+MFPP(:))
        KIJ(:)=TERM5*RMID(:)*DCOEF(:)/(TERMV1(:)+TERMV2(:))
       ENDWHERE
      END IF
!
!  .. For ICOAG=3, KIJ is calculated as in the original UM sulphate
!  .. aerosol scheme where the Cunningham slip correction was assumed
!  .. the same for each mode and the mean free path length for
!  .. particles in each mode assumed to that of air
      IF(ICOAG == 3) THEN
       UKCA_REF=UKCA_MFP_REF*UKCA_PREF_MFP/UKCA_TREF_MFP
       WHERE(MASK(:))
        MFPA_CALC(:)=UKCA_REF*T(:)/PMID(:)
        TERMV3(:)=(2.0e6*ZBOLTZ*T(:)/3.0/DVISC(:))
        TERMV4(:)=1.591*UKCA_MFP_REF*(1.0/RI(:)/RI(:)+1.0/RJ(:)/RJ(:))
        KIJ(:)=TERMV3(:)*(RI(:)+RJ(:))*(1.0/RI(:)+1.0/RJ(:)+TERMV4(:))
       ENDWHERE
      END IF
!
!  .. For ICOAG=4, KIJ is calculated as in the original UM sulphate
!  .. aerosol scheme but values of CCI,CCJ,MFPPI,MFPPJ are computed
!  .. rather than just setting MFPPI=MFPPJ=MFPA and CCI=CCJ=1.591
      IF(ICOAG == 4) THEN
       WHERE(MASK(:))
        KNI(:)=MFPA(:)/RI(:)
        CCI(:)=1.0+KNI(:)*(1.257+0.4*EXP(-1.1/KNI(:)))
! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
! .. UM originally used 1.249, 0.42 and -0.87 as in
! .. Jacobson page 445 [following Kasten (1968) reanalysis of Millikan]
        KNJ(:)=MFPA(:)/RJ(:)
        CCJ(:)=1.0+KNJ(:)*(1.257+0.4*EXP(-1.1/KNJ(:)))
! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
! .. UM originally used 1.249, 0.42 and -0.87 as in
! .. Jacobson page 445 [following Kasten (1968) reanalysis of Millikan]
        TERMV3(:)=(2.0e6*ZBOLTZ*T(:)/3.0/DVISC(:))
        TERMV4(:)=                                                      &
         MFPPI(:)*CCI(:)/RI(:)/RI(:)+MFPPJ(:)*CCJ(:)/RJ(:)/RJ(:)
        KIJ(:)=TERMV3(:)*(RI(:)+RJ(:))*(1.0/RI(:)+1.0/RJ(:)+TERMV4(:))
       ENDWHERE
      END IF

      IF (lhook) CALL dr_hook('UKCA_COAG_COFF_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_COAG_COFF_V
