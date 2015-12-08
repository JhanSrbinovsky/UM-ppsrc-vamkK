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
!    Calculates water content of each mode given component
!    concentrations (in air) using ZSR and binary molalities
!    evaluated using water activity data from Jacobson,
!    "Fundamentals of Atmospheric Modelling", page 610 Table B.
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
      SUBROUTINE UKCA_WATER_CONTENT_V(NV,MASK,CL,RH,IONS,WC)

!-----------------------------------------------------------------------
!
! Purpose:
! -------
! Calculates water content of each mode given component
! concentrations (in air) using ZSR and binary molalities
! evaluated using water activity data from Jacobson,
! "Fundamentals of Atmospheric Modelling", page 610 Table B.
!
! Inputs
! ------
! NV     : Total number of gridboxes in domain
! IONS   : Logical indicating presence of each ion
! CL     : Concentration of each ion (moles/cc air)
! RH     : Relative humidity (fraction)
! MASK   : Logical array for where in domain to do calculation.
!
! Outputs
! -------
! WC     : Water content for aerosol (moles/cm3 of air)
!
! Local Variables
! ---------------
! IC     : Loop variable for cations
! IA     : Loop variable for anions
! AW     : Water activity (local copy of RH fraction)
! CLI    : Internal copy of CL
! CLP    : Ion pair concentrations (moles/cc air)
! MB     : Ion pair solution molalities (moles/cc water)
! N      : Ion stoiciometries for each electrolyte
! Z      : Charge for each ion
! Y      : Coefficients in expressions for binary molalities from
!        : Jacobson page 610 (Table B.10) for each electrolyte
! RH_MIN : Lowest rh for which expression is valid
! MOLAL_MAX : Highest molality for which expression is valid.
!
!----------------------------------------------------------------------
      USE UKCA_MODE_SETUP, ONLY: ncation, nanion
      USE parkind1,        ONLY: jprb, jpim
      USE yomhook,         ONLY: lhook, dr_hook
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: NV                       ! No of points
      LOGICAL, INTENT(IN) :: MASK(NV)                 ! Domain mask
      LOGICAL, INTENT(IN) :: IONS(NV,-NANION:NCATION) ! Ion presence switches
      REAL, INTENT(IN)    :: RH(NV)                   ! Relative humidity
      REAL, INTENT(IN)    :: CL(NV,-NANION:NCATION)   ! Ion conc. (mol/cm3 of air)
      REAL, INTENT(OUT)   :: WC(NV)                   ! water content (mol/cm3 of air)

! Local variables
      INTEGER :: I                   ! Loop counter
      INTEGER :: J                   ! Loop counter
      INTEGER :: IC                  ! Cation loop variable
      INTEGER :: IA                  ! Anion  loop variable
      INTEGER :: M                   ! Counter
      INTEGER :: IDX(NV)             ! Index

!WATER ACTIVITY (%RH EXPRESSED AS A FRACTION)
      REAL    :: AW(NV)
      REAL    :: DUM(NV)
!INTERNAL COPY OF CL
      REAL    :: CLI(NV,-NANION:NCATION)
!ION PAIR CONCENTRATIONS
      REAL    :: CLP(NV,NCATION,-NANION:-1)
!ION PAIR BINARY SOLUTION MOLALITIES AT R
      REAL    :: MB(NV,NCATION,-NANION:-1)
!ION STOICIOMETRIES FOR EACH ELECTROLYTE
      REAL    :: N(-NANION:NCATION)
!ION CHARGES
      REAL    :: Z(-NANION:NCATION)
      REAL    :: Y(3,-4:-1,0:7)
      REAL    :: RH_MIN(3,-4:-1)
      REAL    :: MOLAL_MAX(3,-4:-1)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      DATA (Z(I),I=-NANION,NCATION)/1.0,1.0,2.0,1.0,0.0,1.0,1.0,1.0/
!                                    Cl NO3 SO4 HSO4 H2O  H NH4  Na
!     H+ HSO4- (1,-1)
      DATA (Y(1,-1,J),J=0,7)/                                           &
            3.0391387536E1,-1.8995058929E2, 9.7428231047E2,             &
           -3.1680155761E3, 6.1400925314E3,-6.9116348199E3,             &
            4.1631475226E3,-1.0383424491E3/
      DATA RH_MIN(1,-1),MOLAL_MAX(1,-1)/0.0E0,30.4E0/

!     2H+ SO42- (1,-2)
      DATA (Y(1,-2,J),J=0,7)/                                           &
            3.0391387536E1,-1.8995058929E2, 9.7428231047E2,             &
           -3.1680155761E3, 6.1400925314E3,-6.9116348199E3,             &
            4.1631475226E3,-1.0383424491E3/
      DATA RH_MIN(1,-2),MOLAL_MAX(1,-2)/0.0E0,30.4E0/

!     H+ NO3- (1,-3)
      DATA (Y(1,-3,J),J=0,7)/                                           &
            2.306844303E1,-3.563608869E1,-6.210577919E1,                &
            5.510176187E2,-1.460055286E3, 1.894467542E3,                &
           -1.220611402E2, 3.098597737E2/
      DATA RH_MIN(1,-3),MOLAL_MAX(1,-3)/0.0E0,22.6E0/

!     H+ Cl- (1,-4)
      DATA (Y(1,-4,J),J=0,7)/                                           &
            1.874637647E1,-2.052465972E1,-9.485082073E1,                &
            5.362930715E2,-1.223331346E3, 1.427089861E3,                &
           -8.344219112E2, 1.90992437E2/
      DATA RH_MIN(1,-4),MOLAL_MAX(1,-4)/0.0E0,18.5E0/

!     Na+ HSO4- (2,-1)
      DATA (Y(2,-1,J),J=0,7)/                                           &
            1.8457001681E2,-1.6147765817E3, 8.444076586E3,              &
           -2.6813441936E4, 5.0821277356E4,-5.5964847603E4,             &
            3.2945298603E4,-8.002609678E3/
      DATA RH_MIN(2,-1),MOLAL_MAX(2,-1)/1.9E0,158.0E0/

!     2Na+ SO42- (2,-2)
      DATA (Y(2,-2,J),J=0,7)/                                           &
            5.5983158E2,-2.56942664E3,4.47450201E3,                     &
           -3.45021842E3, 9.8527913E2, 0.0E0,                           &
            0.0E0,        0.0E0/
      DATA RH_MIN(2,-2),MOLAL_MAX(2,-2)/58.0E0,13.1E0/

!     Na+ NO3- (2,-3)
      DATA (Y(2,-3,J),J=0,7)/                                           &
            3.10221762E2,-1.82975944E3, 5.13445395E3,                   &
           -8.01200018E3, 7.07630664E3,-3.33365806E3,                   &
            6.5442029E2,  0.0E0/
      DATA RH_MIN(2,-3),MOLAL_MAX(2,-3)/30.0E0,56.8E0/

!     Na+ Cl- (2,-4)
      DATA (Y(2,-4,J),J=0,7)/                                           &
            5.875248E1,-1.8781997E2, 2.7211377E2,                       &
           -1.8458287E2, 4.153689E1, 0.0E0,                             &
            0.0E0,       0.0E0/
      DATA RH_MIN(2,-4),MOLAL_MAX(2,-4)/47.0E0,13.5E0/

!     NH4+ HSO4- (3,-1)
      DATA (Y(3,-1,J),J=0,7)/                                           &
            2.9997156464E2,-2.8936374637E3, 1.4959985537E4,             &
           -4.5185935292E4, 8.110895603E4, -8.4994863218E4,             &
            4.7928255412E4,-1.1223105556E4/
      DATA RH_MIN(3,-1),MOLAL_MAX(3,-1)/6.5E0,165.0E0/

!     2NH4+ SO42- (3,-2)
      DATA (Y(3,-2,J),J=0,7)/                                           &
            1.1065495E2,-3.6759197E2, 5.0462934E2,                      &
           -3.1543839E2, 6.770824E1,  0.0E0,                            &
            0.0E0,       0.0E0/
      DATA RH_MIN(3,-2),MOLAL_MAX(3,-2)/37.0E0,29.0E0/

!     NH4+ NO3- (3,-3)
      DATA (Y(3,-3,J),J=0,7)/                                           &
            3.983916445E3, 1.153123266E4,-2.13956707E5,                 &
            7.926990533E5,-1.407853405E6, 1.351250086E6,                &
           -6.770046795E5, 1.393507324E5/
      DATA RH_MIN(3,-3),MOLAL_MAX(3,-3)/62.0E0,28.0E0/

!     NH4+ Cl- (3,-4)
      DATA (Y(3,-4,J),J=0,7)/                                           &
           -7.110541604E3, 7.217772665E4,-3.071054075E5,                &
            7.144764216E5,-9.840230371E5, 8.03407288E5,                 &
           -3.603924022E5, 6.856992393E4/
      DATA RH_MIN(3,-4),MOLAL_MAX(3,-4)/47.0E0,23.2E0/

      IF (lhook) CALL dr_hook('UKCA_WATER_CONTENT_V',zhook_in,zhook_handle)

      M = 0
      DO I=1,NV
         IF(MASK(I)) THEN
            M = M+1
            IDX(M) = I
         END IF
      END DO

! Copy fractional relative humidity to water activity (local)
      AW(:M)=RH(IDX(:M))

! Write cl to internal variable to be adjusted in this subroutine only
      DO I=-NANION,NCATION
         CLI(:M,I)=CL(IDX(:M),I)
      END DO
! Calculate mole concentrations of hypothetical ion pairs
      DO IC=1,NCATION
         DO IA=-NANION,-1
! ..Calculate stoichiometries for each ion pair
            N(IC)=Z(IA)
            N(IA)=Z(IC)
            IF(Z(IC) == Z(IA).AND.Z(IA) /= 1.)THEN
               N(IC)=N(IC)/Z(IC)
               N(IA)=N(IA)/Z(IA)
            ENDIF
            WHERE(IONS(IDX(:M),IA).AND.IONS(IDX(:M),IC))
! .. Calculate minimum ion pair concentration and subtract from
! .. Ion concentration
               CLP(:M,IC,IA)=MIN(CLI(:M,IC)/N(IC),CLI(:M,IA)/N(IA))
               CLI(:M,IC)=CLI(:M,IC)-N(IC)*CLP(:M,IC,IA)
               CLI(:M,IA)=CLI(:M,IA)-N(IA)*CLP(:M,IC,IA)
            ENDWHERE
         END DO
      END DO
! Calculate binary electrolyte molalities at given aw
      DO IC=1,NCATION
         DO IA=-NANION,-1
            WHERE(AW(:M) < RH_MIN(IC,IA)/1.0e2)
               AW(:M)=RH_MIN(IC,IA)/1.0e2
            ENDWHERE
            MB(:M,IC,IA)=0.
            WHERE(IONS(IDX(:M),IA).AND.IONS(IDX(:M),IC))
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,0)*AW(:M)**0
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,1)*AW(:M)**1
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,2)*AW(:M)**2
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,3)*AW(:M)**3
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,4)*AW(:M)**4
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,5)*AW(:M)**5
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,6)*AW(:M)**6
               MB(:M,IC,IA)=MB(:M,IC,IA)+Y(IC,IA,7)*AW(:M)**7
               MB(:M,IC,IA)=MIN(MB(:M,IC,IA),MOLAL_MAX(IC,IA))
            ENDWHERE
         END DO
      END DO

! Calculate water content (mol/cm3 air)
      DUM(:M)=0.
      DO IC=1,NCATION
         DO IA=-NANION,-1
            WHERE(IONS(IDX(:M),IA).AND.IONS(IDX(:M),IC))
               DUM(:M)=DUM(:M)+CLP(:M,IC,IA)/MB(:M,IC,IA)
            ENDWHERE
         END DO
      END DO
      WC(IDX(:M))=(1./18.E-3)*DUM(:M)

      IF (lhook) CALL dr_hook('UKCA_WATER_CONTENT_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_WATER_CONTENT_V
