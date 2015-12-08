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
!    Subroutine to calculate impaction scavenging of aerosols
!    by falling raindrops.
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
      SUBROUTINE UKCA_IMPC_SCAV(NBOX,ND,MD,                             &
       CRAIN,DRAIN,WETDP,DTC,BUD_AER_MAS)
!----------------------------------------------------------------------
!
!     Subroutine to calculate impaction scavenging of aerosols
!     by falling raindrops
!
!     Uses empirical expression for the terminal velocity of raindrops
!     of Easter & Hales (1984)
!
!     Uses look-up table of aerosol-raindrop collision efficiencies
!     provided by Yan Yin, University of Aberystwyth. These were
!     originally generated for a specific NROW=19-bin raindrop radius
!     grid colliding with the NCOLL=20-bin aerosol size grid currently
!     used in GLOMAP.
!
!     Currently the routine is only set for a NRBINS=7-bin raindrop
!     size grid consisting of the 4th,6th,8th,10th,12th,14th,16th bin
!     centres of the original 19-bin raindrop size grid. These
!     correspond to raindrop radii of 4.0, 10.08, 25.4, 64.0, 161.3,
!     406.4 and 1024.0 microns and a geometric scaling factor of 2.52.
!
!     The raindrop radii and corresponding aerosol collision
!     efficiencies have already been read in from file and are
!     passed into IMPC_SCAV explicitly as the arrays RADDROP(NROW)
!     and COLLEFF4(NCOLL,NROW)
!
!     Raindrop size distribution assumed is Marshall-Palmer distribution
!     as modified by Sekhon & Srivastava (1971) to take into account
!     rainfall intensity.
!
!     MP have dN/dDp = n_0 * exp (-Psi*Dp) with Psi=4.1*p0^-0.21 mm^-1
!
!       where p0 is precipitation rate (mm/hour)
!             Dp is particle diameter (mm)
!             n_0 is a constant = 8000 drops/m3
!
!     SS modified this to set n_0 to also vary with precipitation rate:
!
!             n_0=7000*p0^0.37 m^-3 mm^-1
!         and Psi= 3.8*p0^-0.14 mm^-1
!
!     Scheme developed using the BCS relation of Slinn 1983, scavenging
!     coeffs supplied by Yan Yin, from Flossmann et al (1985)
!
!     This scheme was developed by Kirsty Pringle as part of PhD thesis.
!
!     References
!     ----------
!     * Flossmann, A. I.; Hall, W. D. & Pruppacher, H. R. (1985)
!       "A theoretical study of the wet removal of atmospheric pollutants:
!       Part 1: The redistribution of aerosol particles captured through
!       nucleation and impaction scavenging by growing cloud drops."
!       J. Atmos. Sci., vol. 42, pp. 582--606.
!
!     * Andronache, C. (2003),
!       "Estimated variability of below-cloud aerosol removal by rainfall
!       for observed aerosol size distribution"
!       Atmos. Chem. Phys., vol. 3, pp. 131--143.
!
!     * Gong, S.-L.; Barry, L. A. & Blanchet J.-P. (1997)
!       "Modelling sea-salt aerosols in the atmosphere,
!        1: Model development!"
!        J. Geophys. Res., vol. 102, no. D3, pp 3,805-3,818
!
!     * Easter, R.C. & Hales, J. M. (1983)
!       "Precipitation scavenging, dry deposition and resuspension,"
!       Chapter Interpretation of the OSCAR data for reactive gas
!       scavenging, pp. 649-662.
!
!     * Sekhon & Srivastava (1971)
!       "Doppler observations of drop size distributions in a thunderstorm"
!       J. Atmos. Sci., vol. 28, pp. 983--984.
!
!     Parameters
!     ----------
!
!     Inputs
!     ------
!     NBOX        : Number of grid boxes
!     ND          : Aerosol ptcl number density for size bin (cm^-3)
!     MD          : Avg cpt mass of aerosol ptcl in size bin (particle^-1)
!     CRAIN       : Convective rain rate array (kgm^-2s^-1)
!     DRAIN       : Dynamic rain rate array (kgm^-2s^-1)
!     WETDP       : Wet diameter corresponding to DRYDP (m)
!     DTC         : Time step of process (s)
!
!     Outputs
!     -------
!     ND          : new aerosol number conc (cm^-3)
!     BUD_AER_MAS : Updated aerosol budgets
!
!     Local variables
!     ---------------
!     TOTRAIN     : Total combined rain rate array (conv + dyn) (mm/hr)
!     RNDPDIAM_cm : Diameter of raindrop (cm)
!     RNDPDIAM_mm : Diameter of raindrop (mm)
!     RNDPDIAM_m  : Diameter of raindrop (m)
!     VELDR_cms   : Terminal velocity of raindrop (cm/s)
!     VELDR_ms    : Terminal velocity of raindrop (m/s)
!     SCAV        : Holds scavenging coefficient for each of the 7 rain bins
!     SCAVCOEFF_COUNT : Holds sum of calculated scavenging coeffs
!                        over all rain bins for each aerosol bin
!     SCAVCOEFF   : Holds final calculated total scavenging coeff for each
!                   aerosol bin summed over all 7 rain bins
!     COUNTCOLL   : Holds index of column (aerosol size) in aerosol-raindrop
!                   collision l-u table for calculated aerosol bin wet radius
!     NRAINMAX    : dN/dD_p [=n0 in Seinfeld & Pandis pg. 832] (m^-3 mm-1)
!     NDRAIN      : dN/d(log D_p) * delta(log(D_p)) [D_p=rndrp diam] (m^-3)
!                   n.b. delta(log(D_p))=ln(2.52)=0.924, where 2.52 is the
!                   geometric scaling factor for the raindrop size grid
!     INTERZZ     : Holds (pi/4)*RNDPDIAM_m^2*VELDR_ms*NDRAIN
!     INTERC      : Holds Psi*D_p in Marshall-Palmer distribution
!     INTERB      : Holds term in determination of column in LUT.
!     R1          : mid-point of first aerosol particle size bin (microns)
!     LNR1        : natural logarithm of R1
!     FACTOR      : geometric scaling factor for aerosol particle size grid
!     FC          : Fraction of grid box over which convective precip occurs.
!     NRBINS      : No. of raindrop size bins used in GLOMAP raindrop spectrum
!     DELN1       : Change in aerosol no. conc. by impaction scav by CRAIN
!     DELN2       : Change in aerosol no. conc. by impaction scav by DRAIN
!     DELN        : Change in aerosol no. conc. by impaction scav (total)
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI         : 3.1415927
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES      : Number of modes set
!     NCP         : Number of components set
!     MODE        : Logical variable defining which modes are set.
!     COMPONENT   : Logical variable defining which components are set.
!     COLLEFF4    : Collision eff for aerosol-raindrop (LU table)
!     NCOLL       : No. of columns in LUT (=20) corresp to aerosol grid
!     NROW        : No. of rows in LUT (=19) corresp to raindrop grid
!     RADDROP     : Raindrop radii at rdrop bin mid-pt (microns)
!     NUM_EPS     : Value of NEWN below which don't recalculate MD
!                                                   or carry out process
!     CP_SU       : Index of component in which SO4    cpt is stored
!     CP_BC       : Index of component in which BC     cpt is stored
!     CP_OC       : Index of component in which 1st OC cpt is stored
!     CP_CL       : Index of component in which NaCl   cpt is stored
!     CP_DU       : Index of component in which dust   cpt is stored
!     CP_SO       : Index of component in which 2nd OC cpt is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,     ONLY: PPI
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE yomhook,            ONLY: lhook, dr_hook
      USE parkind1,           ONLY: jprb, jpim
      IMPLICIT NONE

! .. Subroutine interface
      INTEGER, INTENT(IN) :: NBOX
      REAL, INTENT(IN)    :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: DTC
      REAL, INTENT(IN)    :: CRAIN(NBOX)
      REAL, INTENT(IN)    :: DRAIN(NBOX)
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)

! .. Local variables
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: JL
      INTEGER :: JVR
      INTEGER :: IROW
      INTEGER :: ICOLL
      INTEGER :: COUNTCOLL(NBOX,NMODES)
      INTEGER, PARAMETER :: NRBINS=7
      REAL    :: LNR1
      REAL    :: R1
      REAL    :: FACTOR
      REAL    :: INTERZZ(NBOX,NRBINS)
      REAL    :: SCAV(NBOX,NMODES,NRBINS)
      REAL    :: SCAVCOEFF(NBOX,NMODES)
      REAL    :: SCAVCOEFF_COUNT(NBOX,NMODES,NRBINS)
      REAL    :: INTERC(NBOX,NRBINS)
      REAL    :: INTERB
      REAL    :: VELDR_ms(NRBINS)
      REAL    :: VELDR_cms(NRBINS)
      REAL    :: RNDPDIAM_cm(NRBINS)
      REAL    :: RNDPDIAM_mm(NRBINS)
      REAL    :: RNDPDIAM_m(NRBINS)
      REAL    :: NDRAIN(NBOX,NRBINS)
      REAL    :: NRAINMAX(NBOX)
      REAL    :: TOTRAIN(NBOX)
      REAL    :: DELN
      REAL    :: DELN1
      REAL    :: DELN2
      REAL, PARAMETER :: FC=0.3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! .. Here convert rain rate from kg/m2/s to mm/hr
! .. Combine the two rain rates (con and dyn)
      IF (lhook) CALL dr_hook('UKCA_IMPC_SCAV',zhook_in,zhook_handle)
      DO JL=1,NBOX
       TOTRAIN(JL)=(CRAIN(JL)*3600.)+(DRAIN(JL)*3600.)
      ENDDO

! .. In LUT there are 19 Rows but we only have NRBINS=7 raindrop bins,
! .. so use every 2nd row starting from row 4
! .. Only consider raindrops radius >4 and <1024 um
! .. NB/ HAVE TO CHANGE THIS IF NRBINS IS ALTERED

      JVR=0
      DO IROW=4,16,2 ! pick out 7 sizes corresponding to 7 raindrop bins
       JVR=JVR+1
       RNDPDIAM_cm(JVR)=RADDROP(IROW)*2.0/10000.0
       RNDPDIAM_mm(JVR)=RNDPDIAM_cm(JVR)*10.0
       RNDPDIAM_m(JVR)=RNDPDIAM_mm(JVR)/1000.0
! .. RADDROP contains raindrop bin radii in microns
! .. RNDPDIAM_cm contains raindrop bin diameters in cm
! .. RNDPDIAM_mm & RNDPDIAM_m contain rndrp bin diams in mm & m respect.
      END DO

! .. Initialise Arrays
      SCAVCOEFF_COUNT(:,:,:)=0.0
      SCAVCOEFF(:,:)=0.0
      COUNTCOLL(:,:)=0

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        DO JL=1,NBOX
         IF(TOTRAIN(JL) > 0.0)THEN

! .. From aerosol wet radius find which column of LU table is needed
          IF(0.5*WETDP(JL,IMODE) < 0.001E-6)THEN
           COUNTCOLL(JL,IMODE)=1
          ELSEIF(0.5*WETDP(JL,IMODE) > 9.0E-6)THEN
           COUNTCOLL(JL,IMODE)=10
          ELSE
           FACTOR=0.457
           R1=0.001
           LNR1=-6.9
           INTERB=((0.5*WETDP(JL,IMODE)*1.0E6)/R1)
! .. 1.0E6 in above to convert wet radius to um
           COUNTCOLL(JL,IMODE)=INT((LOG(INTERB)/FACTOR)+1.)
          END IF
         END IF
        END DO
       END IF
      END DO

! .. Loop to calculate the number of rain drops at each rdrop bin-centre
! .. following Marshall Palmer distribution
      DO JVR=1,NRBINS
!
! .. Empirical relationship from Easter and Hales (1984) to calulate
! .. terminal velocity of raindrop in cm/s
        IF(RNDPDIAM_cm(JVR) <= 0.10)THEN
          VELDR_cms(JVR)=4055.0*RNDPDIAM_cm(JVR)
        ELSEIF(RNDPDIAM_cm(JVR) > 0.10)THEN
          VELDR_cms(JVR)=13000.0*(RNDPDIAM_m(JVR)**0.5)
        END IF
! .. Convert raindrop velocity from cm/s to m/s
        VELDR_ms(JVR)=VELDR_cms(JVR)/100.0
      END DO

      DO JL=1,NBOX
       IF(TOTRAIN(JL) > 0.0)THEN
        DO JVR=1,NRBINS

! .. Marshall Palmer but with sophistication of NRAINMAX calculated
! .. according to Sekhon & Srivastava (1971) (Seinfeld & Pandis pg 832)

         NRAINMAX(JL)=7000.0*TOTRAIN(JL)**0.37
! NRAINMAX is n_0 in m^-3 per mm

         INTERC(JL,JVR)=3.8*TOTRAIN(JL)**(-0.14)*RNDPDIAM_mm(JVR)
! INTERC is Psi*D_p in MP distribution


         NDRAIN(JL,JVR)=NRAINMAX(JL)*RNDPDIAM_mm(JVR)                   &
                       *0.924*EXP(-INTERC(JL,JVR))
! .. Includes various conversions as Nd given in terms of mm-1
! .. => *diameter(mm)*width of bin
! .. Geometric scaling factor = 2.52 (change if NRBINS is changed)
! .. (ln(2.52)=0.924)
!
! .. Leave NDRAIN in m-3 for use in SCAV calculation below
! .. Convert Raindrop size from diameter in mm to diameter in m

         INTERZZ(JL,JVR)=                                               &
            (PPI/4.0)*((RNDPDIAM_m(JVR)*RNDPDIAM_m(JVR))                &
            *VELDR_ms(JVR)*NDRAIN(JL,JVR))

        END DO ! over NRBINS
       END IF ! TOTRAIN>0
      END DO ! over NBOX

! .. Calculate scavenging coefficients
      DO JL=1,NBOX
       IF(TOTRAIN(JL) > 0.0)THEN
        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN
          IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
! .. Loop over raindrop bins
           DO JVR=1,NRBINS
            IROW=(JVR*2)+2 ! hard-wired for 7 raindrop size bins
! .. COUNTCOLL contains the index of the relevant column of the LU table
! .. for aerosol-rndrp coll. efficiency for each aerosol bin wet radius
            ICOLL=COUNTCOLL(JL,IMODE)
            IF(COUNTCOLL(JL,IMODE) == 0) ICOLL=1
            SCAV(JL,IMODE,JVR)=COLLEFF4(ICOLL,IROW)*INTERZZ(JL,JVR)

! .. For each aerosol bin, sum calculated scavenging coeff over all rain bins
            IF(JVR == 1) THEN
             SCAVCOEFF_COUNT(JL,IMODE,JVR)=SCAV(JL,IMODE,JVR)
            END IF
            IF(JVR > 1) THEN
             SCAVCOEFF_COUNT(JL,IMODE,JVR)=                             &
         SCAVCOEFF_COUNT(JL,IMODE,JVR-1)+SCAV(JL,IMODE,JVR)
            END IF
            IF(JVR == NRBINS) THEN
             SCAVCOEFF(JL,IMODE)=SCAVCOEFF_COUNT(JL,IMODE,JVR)
            END IF

           END DO ! LOOP OVER RAINDROP BINS
          END IF ! IF ND>ND_EPS
         END IF ! IF MODE PRESENT
        END DO ! LOOP OVER MODES
       END IF ! IF TOTRAIN>0
      END DO ! LOOP OVER BOXES

! .. Below calculates the rate of removal
!
! .. Calculate removal (from each bin) following 1st order rate loss
! .. Apply BCS only to the portion of the gridbox where rain is occuring
! .. If dynamic rain then cloud cover = 1.0 => apply BCS to all aerosols
! .. If convective rain then cloud cover = FC = 0.3

      DO JL=1,NBOX
       IF(TOTRAIN(JL) > 0.0)THEN
        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN
          IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
           IF(CRAIN(JL) > 0.0) THEN
            DELN1=FC*ND(JL,IMODE)*(1.0-EXP(-SCAVCOEFF(JL,IMODE)*DTC))
            ND(JL,IMODE)=ND(JL,IMODE)-DELN1
            IF(ND(JL,IMODE) < 0.0) ND(JL,IMODE)=0.0
! .. Convective rain occurs only over fraction FC of box
           ELSE
            DELN1=0.0
           END IF
           IF(DRAIN(JL) > 0.0) THEN
            DELN2=ND(JL,IMODE)*(1.0-EXP(-SCAVCOEFF(JL,IMODE)*DTC))
            ND(JL,IMODE)=ND(JL,IMODE)-DELN2
            IF(ND(JL,IMODE) < 0.0) ND(JL,IMODE)=0.0
! .. Dynamic rain occurs over whole of box
           ELSE
            DELN2=0.0
           END IF

           DELN=DELN1+DELN2

! .. Store cpt impaction scavenging mass fluxes for budget calculations
           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             IF(ICP == CP_SU) THEN
              IF((IMODE == 1).AND.(NMASIMSCSUNUCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSUNUCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSUNUCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 2).AND.(NMASIMSCSUAITSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSUAITSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSUAITSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 3).AND.(NMASIMSCSUACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSUACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSUACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCSUCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSUCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSUCORSOL)+DELN*MD(JL,IMODE,ICP)
             END IF
             IF(ICP == CP_BC) THEN
              IF((IMODE == 2).AND.(NMASIMSCBCAITSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCBCAITSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCBCAITSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 3).AND.(NMASIMSCBCACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCBCACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCBCACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCBCCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCBCCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCBCCORSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 5).AND.(NMASIMSCBCAITINS > 0))               &
        BUD_AER_MAS(JL,NMASIMSCBCAITINS)=                               &
        BUD_AER_MAS(JL,NMASIMSCBCAITINS)+DELN*MD(JL,IMODE,ICP)
             END IF
             IF(ICP == CP_OC) THEN
              IF((IMODE == 1).AND.(NMASIMSCOCNUCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCOCNUCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCOCNUCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 2).AND.(NMASIMSCOCAITSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCOCAITSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCOCAITSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 3).AND.(NMASIMSCOCACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCOCACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCOCACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCOCCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCOCCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCOCCORSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 5).AND.(NMASIMSCOCAITINS > 0))               &
        BUD_AER_MAS(JL,NMASIMSCOCAITINS)=                               &
        BUD_AER_MAS(JL,NMASIMSCOCAITINS)+DELN*MD(JL,IMODE,ICP)
             END IF
             IF(ICP == CP_CL) THEN
              IF((IMODE == 3).AND.(NMASIMSCSSACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSSACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSSACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCSSCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSSCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSSCORSOL)+DELN*MD(JL,IMODE,ICP)
             END IF
             IF(ICP == CP_SO) THEN
              IF((IMODE == 1).AND.(NMASIMSCSONUCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSONUCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSONUCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 2).AND.(NMASIMSCSOAITSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSOAITSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSOAITSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 3).AND.(NMASIMSCSOACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSOACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSOACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCSOCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCSOCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCSOCORSOL)+DELN*MD(JL,IMODE,ICP)
             END IF
             IF(ICP == CP_DU) THEN
              IF((IMODE == 3).AND.(NMASIMSCDUACCSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCDUACCSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCDUACCSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 4).AND.(NMASIMSCDUCORSOL > 0))               &
        BUD_AER_MAS(JL,NMASIMSCDUCORSOL)=                               &
        BUD_AER_MAS(JL,NMASIMSCDUCORSOL)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 6).AND.(NMASIMSCDUACCINS > 0))               &
        BUD_AER_MAS(JL,NMASIMSCDUACCINS)=                               &
        BUD_AER_MAS(JL,NMASIMSCDUACCINS)+DELN*MD(JL,IMODE,ICP)
              IF((IMODE == 7).AND.(NMASIMSCDUCORINS > 0))               &
        BUD_AER_MAS(JL,NMASIMSCDUCORINS)=                               &
        BUD_AER_MAS(JL,NMASIMSCDUCORINS)+DELN*MD(JL,IMODE,ICP)
             END IF
            END IF ! if component present
           END DO ! loop over components
          END IF ! if ND>NUM_EPS
         END IF ! if distribution defined
        END DO ! loop over distributions
       END IF ! if rain is present in box
      END DO ! loop over boxes

      IF (lhook) CALL dr_hook('UKCA_IMPC_SCAV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_IMPC_SCAV
