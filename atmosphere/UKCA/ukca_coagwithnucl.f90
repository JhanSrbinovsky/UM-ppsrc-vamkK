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
!    Calculate the change in number and mass concentration of
!    insoluble modes due to a combination of coagulation and nucleation.
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
      SUBROUTINE UKCA_COAGWITHNUCL(NBOX,ND,MD,MDT,DELGC_NUCL,DTZ,JRATE, &
        AGETERM2,INTRAOFF,INTEROFF,VERBOSE,BUD_AER_MAS,KII_ARR,KIJ_ARR)

!---------------------------------------------------------------
!
! Purpose
! -------
! Calculate the change in number and mass concentration of
! insoluble modes due to a combination of coagulation and nucleation.
!
! Newly nucleated particles assumed to be 100 H2SO4 molecules
! per particle which equates to 5nm diameter.
!
! Coagulation equations are of the form dN/dt = a*N^2 + b*N + c
!
! Can rearrange this to be in the form of an indefinite integral
!
! int_{N_0}^{N} dx/X = int_{t_0}^{t_0+deltat} dt = deltat
!
! where X= A*x^2 + B*x + C
!
! Then solve this analytically (see header to ukca_solvecoagnucl_v.F90)
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NBOX    : Number of grid boxes
! ND      : Initial no. concentration of aerosol mode (ptcls/cc)
! MD      : Initial avg cpt mass conc of aerosol mode (molecules/ptcl)
! MDT     : Initial avg total mass conc of aerosol mode (molecules/ptcl)
! DELGC_NUCL : Change in H2SO4 (g) due to nucleation (molecules/cc)
! DTZ     : nucl/coag/cond time step (s)
! JRATE   : H2SO4-H2O binary nucleation rate (molecules/cc/s)
! INTRAOFF: Switch to turn off intra-modal coagulation
! INTEROFF: Switch to turn off intra-modal coagulation
! VERBOSE : Switch for whether to do various test print statements
! KII_ARR : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
! KIJ_ARR : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
!
! Outputs
! -------
! ND      : Updated no. concentration of aerosol mode (ptcls/cc)
! MD      : Updated avg cpt mass conc of aerosol mode (molecules/ptcl)
! MDT     : Updated avg total mass conc of aerosol mode (molecules/ptcl)
! AGETERM2: Rate of accomodation of material to each insoluble mode
!           as a result of coagulation with smaller soluble modes
!           (in molecules cpt /cm3/DTZ)
!           Used for calculation of ageing rate in UKCA_AGEING.
! BUD_AER_MAS : Aerosol mass budgets
!
! Local variables
! ---------------
! RPI       : Geometric mean radius for mode IMODE (m)
! RPJ       : Geometric mean radius for mode JMODE (m)
! VPI       : Volume of particle of radius RPI (m^3)
! VPJ       : Volume of particle of radius RPJ (m^3)
! NDOLD     : Initial aerosol ptcl number conc (ptcls/cc)
! NDOLD_V   : Initial aerosol ptcl number conc (ptcls/cc) [1D version]
! MDOLD     : Initial aerosol cpt mass/ptcl (molecules/ptcl)
! MDCPOLD   : Initial aerosol cpt mass conc. (molecules/cm3)
! MDCPNEW   : New aerosol cpt mass conc. (molecules/cm3)
! MTRAN     : Stores transfer of cpt masses between modes by coag (/cm3)
! MTRANFMI  : Stores transfer of cpt masses from mode (/cm3)
! MTRANTOI  : Stores transfer of cpt masses to   mode (/cm3)
! MTRANNET  : Stores net transfer of cpt masses for mode (/cm3)
! DELN      : Change in ND (from NDOLD) due to coag & nucl (ptcls/cc)
! A         : Constant in nucl/coag equation
! B         : Constant in nucl/coag equation
! C         : Constant in nucl/coag equation
! BTERM     : Inter-modal coag rate term for single i-j
! KII       : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
! KIJ       : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
! XXX       : Term in exponential (1.0-EXP(-XXX(:)))
! XXX_EPS   : Tolerance for XXX below which don't evaluate exponential
! MASK1,MASK2,.. : Logicals to define domain regions for where loops
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! NMOL      : Number of molecules per particle at nucleation
! CONC_EPS  : Threshold for GC to calc nucl+coag (molecules per cc)
! DN_EPS    : Value of DELN below which do not carry out process
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! COAG_MODE : Defines which mode an IMODE-JMODE coagulation goes into
! MMID      : Mid-point masses for initial radius grid
! MFRAC_0   : Initial mass fraction to set when no particles.
! MODESOL   : Defines whether the mode is soluble or not (1 or 0)
! NUM_EPS   : Value of NEWN below which do not carry out process
! CP_SU     : Index of component in which sulfate is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which OC is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_DU     : Index of component in which dust is stored
! CP_SO     : Index of component in which condensible organic is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,   ONLY: nmol, conc_eps, dn_eps
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE yomhook,          ONLY: lhook, dr_hook
      USE parkind1,         ONLY: jprb, jpim
      IMPLICIT NONE

! Subroutine interface:
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: INTRAOFF
      INTEGER, INTENT(IN) :: INTEROFF
      INTEGER, INTENT(IN) :: VERBOSE
      REAL, INTENT(IN)    :: DELGC_NUCL(NBOX,NCHEMG)
      REAL, INTENT(IN)    :: DTZ
      REAL, INTENT(IN)    :: JRATE(NBOX)
      REAL, INTENT(IN)    :: KII_ARR(NBOX,NMODES)
      REAL, INTENT(IN)    :: KIJ_ARR(NBOX,NMODES,NMODES)
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: AGETERM2(NBOX,4,3,NCP)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)

!  .. Local variables
      INTEGER :: IMODE
      INTEGER :: JMODE
      INTEGER :: ICP
      INTEGER :: JCP
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK1A(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK3(NBOX)
      LOGICAL :: MASK4(NBOX)
      REAL    :: NDOLD(NBOX,NMODES)
      REAL    :: NDOLD_V(NBOX)
      REAL    :: DELN(NBOX)
      REAL    :: MDOLD(NBOX,NMODES,NCP)
      REAL    :: MTRAN(NBOX,NMODES,NMODES,NCP)
      REAL    :: MTRANFMI(NBOX,NMODES,NCP)
      REAL    :: MTRANTOI(NBOX,NMODES,NCP)
      REAL    :: MTRANNET(NBOX)
      REAL    :: MDCPOLD(NBOX)
      REAL    :: MDCPNEW(NBOX)
      REAL    :: A(NBOX)
      REAL    :: B(NBOX)
      REAL    :: C(NBOX)
      REAL    :: BTERM(NBOX)
      REAL    :: KII(NBOX)
      REAL    :: KIJ(NBOX)
      REAL    :: XXX(NBOX)
      REAL    :: DELH2SO4_NUCL(NBOX)
      REAL, PARAMETER :: XXX_EPS=1.0e-3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_COAGWITHNUCL',zhook_in,zhook_handle)

! Copy H2SO4 values from DELGC_NUCL to local variable
      IF(MH2SO4 > 0) THEN
       DELH2SO4_NUCL(:)=DELGC_NUCL(:,MH2SO4)
      ELSE
       DELH2SO4_NUCL(:)=0.0
      END IF

      BTERM(:) = 0.

      DO IMODE=1,NMODES
       NDOLD(:,IMODE)=ND(:,IMODE)
       DO ICP=1,NCP
        MDOLD(:,IMODE,ICP)=MD(:,IMODE,ICP)
        DO JMODE=1,NMODES
         MTRAN(:,IMODE,JMODE,ICP)=0.0
        END DO
       END DO
      END DO

      DO IMODE=1,4
       IF(MODE(IMODE)) THEN

        KII(:)=KII_ARR(:,IMODE) ! copy in pre-calculated KII

        A(:)=0.0
        B(:)=0.0
        C(:)=0.0

        ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
        MASK1(:) =(NDOLD(:,IMODE) > NUM_EPS(IMODE))

        IF(INTRAOFF /= 1) THEN
         WHERE(MASK1(:)) A(:)=-0.5*KII(:)
        ENDIF

        DO JMODE=(IMODE+1),4 ! inter-coag with larger soluble modes
         IF(MODE(JMODE)) THEN

          KIJ(:)=KIJ_ARR(:,IMODE,JMODE) ! copy in pre-calculated KIJ

          ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
          ! and NDOLD(:,JMODE) > NUM_EPS(IMODE)
          MASK2(:) =(NDOLD(:,IMODE) > NUM_EPS(IMODE)).AND.              &
                    (NDOLD(:,JMODE) > NUM_EPS(JMODE))

          IF(INTEROFF /= 1) THEN
           WHERE(MASK2(:))
            BTERM(:)=-KIJ(:)*NDOLD(:,JMODE)
            B(:)=B(:)+BTERM(:)
           ENDWHERE
          END IF
          XXX(:)=-BTERM(:)*DTZ
          MASK4(:)=ABS(XXX(:)) > XXX_EPS
! .. above only evaluates exponential where it is "worth it"
! .. (in cases where XXX is larger than specified tolerance XXX_EPS)
          DO ICP=1,NCP
           IF(COMPONENT(JMODE,ICP)) THEN
!
            WHERE(MASK4(:).AND.MASK2(:))                                &
              MTRAN(:,IMODE,JMODE,ICP)=                                 &
              MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*(1.0-EXP(-XXX(:)))
!
            WHERE((.NOT.MASK4(:)).AND.MASK2(:))                         &
              MTRAN(:,IMODE,JMODE,ICP)=                                 &
              MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*XXX(:)
!
! .. MTRAN(:,IMODE,JMODE,ICP) transfers from IMODE to JMODE
!
           END IF
          END DO ! loop over components
         END IF ! if (MODE(JMODE))
        END DO ! end loop of JMODE over soluble modes

        DO JMODE=(IMODE+4),7 ! inter-coag with larger insoluble modes
         IF(MODE(JMODE)) THEN

          KIJ(:)=KIJ_ARR(:,IMODE,JMODE) ! copy in pre-calculated KIJ

          ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
          ! and NDOLD(:,JMODE) > NUM_EPS
          MASK2(:) =(NDOLD(:,IMODE) > NUM_EPS(IMODE)).AND.              &
                    (NDOLD(:,JMODE) > NUM_EPS(JMODE))

          IF(INTEROFF /= 1) THEN
           WHERE(MASK2(:))
            BTERM(:)=-KIJ(:)*NDOLD(:,JMODE)
            B(:)=B(:)+BTERM(:)
           ENDWHERE
          END IF
! .. calculate MTRAN for sol-ins inter-modal coag to store in AGETERM2
! .. & carry out transfer of mass at end of the subroutine by MTRANNET.
! .. Also transfer number (include sol-ins term in B summation)
          XXX(:)=-BTERM(:)*DTZ
          MASK4(:)=ABS(XXX(:)) > XXX_EPS
! .. above only evaluates exponential where it is "worth it"
! .. (in cases where XXX is larger than specified tolerance XXX_EPS).
          DO ICP=1,NCP
           WHERE(MASK4(:).AND.MASK2(:))                                 &
       MTRAN(:,IMODE,JMODE,ICP)=                                        &
       MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*(1.0-EXP(-XXX(:)))
           WHERE((.NOT.MASK4(:)).AND.MASK2(:))                          &
       MTRAN(:,IMODE,JMODE,ICP)=                                        &
       MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*XXX(:)
          END DO
         END IF
        END DO ! end loop of JMODE over insoluble modes

        ! reset mass if insignificant # of particles
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(.NOT.MASK1(:))                                          &
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
         END IF
        END DO
        WHERE(.NOT.MASK1(:)) MDT(:,IMODE)=MMID(IMODE)

        IF(IMODE == 1) THEN

         WHERE(DELH2SO4_NUCL(:) > CONC_EPS)
          C(:)=DELH2SO4_NUCL(:)/DTZ/NMOL
         ENDWHERE ! if nucl mode & some particles are to be nucleated

        END IF

        MASK1A(:)=MASK1(:).OR.                                          &
                  ((IMODE == 1).AND.(DELH2SO4_NUCL(:) > CONC_EPS))

        NDOLD_V(:)=NDOLD(:,IMODE)
! DEPENDS ON: ukca_solvecoagnucl_v
        CALL UKCA_SOLVECOAGNUCL_V(NBOX,MASK1A,A,B,C,NDOLD_V,DTZ,DELN)

        WHERE(ABS(DELN(:)) > DN_EPS)
         ND(:,IMODE)=NDOLD(:,IMODE)+DELN(:)
        ENDWHERE

       END IF ! if mode is defined
      END DO ! end loop of IMODE over soluble modes
            ! (have updated ND here but not MD,MDT)
            ! (have stored transfer of mass from I->J in MTRAN(I,J)
            ! (updating of MD,MDT done at end of subroutine)
!
! .. Below is the section which has insoluble modes coagulating
! .. either with themselves (intra-modal) or inter-modally
! .. with larger soluble modes (inter-modal with ins neglected).
! .. Also the AGETERM2 is stored for the mass that coagulates
! .. with each insoluble mode from the smaller soluble modes.
!
! Solve for change in no. & mass due to intra and inter-modal coag here.
! Transfer of soluble mass from each soluble mode to each insoluble mode
! by inter-modal coag is stored here in AGETERM2 for use in UKCA_AGEING.
! Inter-modal coagulation with insoluble modes considered inefficient
! and neglected (as in M7).
!
! .. The section below stores the inter-modal coagulation between
! .. soluble modes and larger insoluble modes for calculation of ageing.
! .. Note that mass will be transferred to corresponding soluble mode
! .. (as prescribed by COAG_MODE) at end of coagulation subroutine
! .. because all soluble material is transferred in timestep by ageing
!
      DO IMODE=5,NMODES
       IF(MODE(IMODE)) THEN

        AGETERM2(:,:,IMODE-4,:)=0.0

        DO JMODE=1,4
         IF(MODE(JMODE)) THEN
          DO JCP=1,NCP
           IF(COMPONENT(JMODE,JCP)) THEN
            AGETERM2(:,JMODE,IMODE-4,JCP)=MTRAN(:,JMODE,IMODE,JCP)
           END IF
          END DO
         END IF
        END DO

! .. The section below calculates change in no. & mass due to
! .. intra-modal coag of insoluble modes and inter-modal coag with
! .. larger soluble modes (this does not contribute to AGETERM2 because
! .. it is automatically passed over when it coagulates with the larger
! .. soluble particles). Note that the number is updated here whereas
! .. the mass is updated at the end of the subroutine (stored in MTRAN).

        KII(:)=KII_ARR(:,IMODE) ! copy in pre-calculated KII

        A(:)=0.0
        B(:)=0.0
        C(:)=0.0

        ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
        MASK1(:) =(NDOLD(:,IMODE) > NUM_EPS(IMODE))

        IF(INTRAOFF /= 1) THEN
         WHERE(MASK1(:)) A(:)=-0.5*KII(:)
        ENDIF

        DO JMODE=(IMODE-2),4 ! ins inter-coag with larger soluble modes
         IF(MODE(JMODE)) THEN

          KIJ(:)=KIJ_ARR(:,IMODE,JMODE) ! copy in pre-calculated KIJ

          ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
          ! and NDOLD(:,JMODE) > NUM_EPS
          MASK2(:) =(NDOLD(:,IMODE) > NUM_EPS(IMODE)).AND.              &
                    (NDOLD(:,JMODE) > NUM_EPS(JMODE))

          IF(INTEROFF /= 1) THEN
           WHERE(MASK2(:))
            BTERM(:)=-KIJ(:)*NDOLD(:,JMODE)
            B(:)=B(:)+BTERM(:)
           ENDWHERE
          END IF
! .. calculate MTRAN for ins-sol inter-modal coag for carrying
! .. out transfer of mass at end of the subroutine by MTRANNET.
! .. Also transfer number (include ins-sol term in B summation)
          XXX(:)=-BTERM(:)*DTZ
          MASK4(:)=ABS(XXX(:)) > XXX_EPS
! .. above only evaluates exponential where it is "worth it"
! .. (in cases where XXX is larger than specified tolerance XXX_EPS)
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            WHERE(MASK4(:).AND.MASK2(:))                                &
             MTRAN(:,IMODE,JMODE,ICP)=                                  &
              MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*(1.0-EXP(-XXX(:)))
            WHERE((.NOT.MASK4(:)).AND.MASK2(:))                         &
             MTRAN(:,IMODE,JMODE,ICP)=                                  &
              MDOLD(:,IMODE,ICP)*NDOLD(:,IMODE)*XXX(:)
           END IF
          END DO
         END IF
        END DO ! end loop of JMODE over larger soluble modes

        NDOLD_V(:)=NDOLD(:,IMODE)
! DEPENDS ON: ukca_solvecoagnucl_v
        CALL UKCA_SOLVECOAGNUCL_V(NBOX,MASK1,A,B,C,NDOLD_V,DTZ,DELN)

        WHERE(MASK1(:) .AND. ABS(DELN(:)) > DN_EPS)
         ND(:,IMODE)=NDOLD(:,IMODE)+DELN(:)
        ENDWHERE

        ! reset mass if insignificant # of particles
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(.NOT.MASK1(:))                                          &
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
         END IF
        END DO
        WHERE(.NOT.MASK1(:)) MDT(:,IMODE)=MMID(IMODE)

       END IF ! end of IF(MODE(IMODE))
      END DO ! end loop of IMODE over insoluble modes
!
      DO ICP=1,NCP
       DO IMODE=1,NMODES
        MTRANFMI(:,IMODE,ICP)=0.0
        MTRANTOI(:,IMODE,ICP)=0.0
       END DO
      END DO
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          DO JMODE=1,NMODES
           IF(MODE(JMODE)) THEN
            MTRANFMI(:,IMODE,ICP)=MTRANFMI(:,IMODE,ICP)+                &
                 MTRAN(:,IMODE,JMODE,ICP)
            MTRANTOI(:,COAG_MODE(IMODE,JMODE),ICP)=                     &
                 MTRANTOI(:,COAG_MODE(IMODE,JMODE),ICP)+                &
                 MTRAN(:,IMODE,JMODE,ICP)
           END IF
          END DO
         END IF
        END DO
       END IF
      END DO

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN

        ! Calculations are only done where ND(:,IMODE) > NUM_EPS
        MASK1(:) =(ND(:,IMODE) > NUM_EPS(IMODE))

        MDT(:,IMODE)=0.0
        MTRANNET(:)=0.0
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          IF(IMODE == 1) THEN
           IF(ICP == 1) THEN ! add nucleated mass to sulfate cpt
            MTRANNET(:)=MTRANTOI(:,IMODE,ICP)-                          &
                        MTRANFMI(:,IMODE,ICP)+DELH2SO4_NUCL(:)
           ELSE
            MTRANNET(:)=MTRANTOI(:,IMODE,ICP)-MTRANFMI(:,IMODE,ICP)
           END IF
          ELSE
           MTRANNET(:)=MTRANTOI(:,IMODE,ICP)-MTRANFMI(:,IMODE,ICP)
          END IF
          MDCPOLD(:)=NDOLD(:,IMODE)*MDOLD(:,IMODE,ICP)
          MDCPNEW(:)=MDCPOLD(:)+MTRANNET(:)
! where MDCPNEW<0, set ND to zero (MDT,MD reset at end routine)
          WHERE(MASK1(:) .AND. (MDCPNEW(:) < 0.0))
           ND(:,IMODE)=0.0
           MASK1(:)=.FALSE. ! set false so not used for other icp values
! n.b. MD and MDT will be re-set at end of routine
          ENDWHERE
          MASK3(:)=MASK1(:).AND.(MDCPNEW(:) >= 0.0)
          WHERE(MASK3(:))
           MD(:,IMODE,ICP)=MDCPNEW(:)/ND(:,IMODE)
           MDT(:,IMODE)=MDT(:,IMODE)+MD(:,IMODE,ICP)
          ENDWHERE
          DO JMODE=1,NMODES
           IF(MODE(JMODE)) THEN
            IF((IMODE == 1).AND.(JMODE == 2)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR12 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR12)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR12)+MTRAN(:,IMODE,JMODE,CP_SU)
             ENDIF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR12 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR12)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR12)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR12 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR12)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR12)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=1,2 (from mode 1 to mode 2)
            IF((IMODE == 1).AND.(JMODE == 3)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR13 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR13)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR13)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR13 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR13)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR13)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR13 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR13)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR13)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=1,3 (from mode 1 to mode 3)
            IF((IMODE == 1).AND.(JMODE == 4)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR14 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR14)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR14)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR14 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR14)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR14)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR14 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR14)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR14)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=1,4 (from mode 1 to mode 4)
            IF((IMODE == 1).AND.(JMODE == 5)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR15 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR15)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR15)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR15 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR15)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR15)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR15 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR15)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR15)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=1,5 (from mode 1 to mode 5)
            IF((IMODE == 1).AND.(JMODE == 6)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR16 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR16)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR16)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR16 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR16)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR16)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR16 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR16)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR16)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=1,6 (from mode 1 to mode 6)
            IF((IMODE == 1).AND.(JMODE == 7)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR17 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR17)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR17)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR17 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR17)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR17)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR17 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR17)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR17)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            ENDIF ! IF IMODE,JMODE=1,7 (from mode 1 to mode 7)
            IF((IMODE == 2).AND.(JMODE == 3)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR23 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR23)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR23)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_BC).AND.(NMASCOAGBCINTR23 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGBCINTR23)=                                 &
       BUD_AER_MAS(:,NMASCOAGBCINTR23)+MTRAN(:,IMODE,JMODE,CP_BC)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR23 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR23)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR23)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR23 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR23)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR23)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=2,3 (from mode 2 to mode 3)
            IF((IMODE == 2).AND.(JMODE == 4)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR24 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR24)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR24)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_BC).AND.(NMASCOAGBCINTR24 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGBCINTR24)=                                 &
       BUD_AER_MAS(:,NMASCOAGBCINTR24)+MTRAN(:,IMODE,JMODE,CP_BC)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR24 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR24)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR24)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR24 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR24)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR24)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
            END IF ! IF IMODE,JMODE=2,4 (from mode 2 to mode 4)
            IF((IMODE == 3).AND.(JMODE == 4)) THEN
             IF((ICP == CP_SU).AND.(NMASCOAGSUINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSUINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGSUINTR34)+MTRAN(:,IMODE,JMODE,CP_SU)
             END IF
             IF((ICP == CP_BC).AND.(NMASCOAGBCINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGBCINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGBCINTR34)+MTRAN(:,IMODE,JMODE,CP_BC)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR34)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
             IF((ICP == CP_CL).AND.(NMASCOAGSSINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSSINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGSSINTR34)+MTRAN(:,IMODE,JMODE,CP_CL)
             END IF
             IF((ICP == CP_SO).AND.(NMASCOAGSOINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGSOINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGSOINTR34)+MTRAN(:,IMODE,JMODE,CP_SO)
             END IF
             IF((ICP == CP_DU).AND.(NMASCOAGDUINTR34 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGDUINTR34)=                                 &
       BUD_AER_MAS(:,NMASCOAGDUINTR34)+MTRAN(:,IMODE,JMODE,CP_DU)
             END IF
            END IF ! IF IMODE,JMODE=3,4 (from mode 3 to mode 4)
            IF((IMODE == 5).AND.(JMODE == 3)) THEN
             IF((ICP == CP_BC).AND.(NMASCOAGBCINTR53 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGBCINTR53)=                                 &
       BUD_AER_MAS(:,NMASCOAGBCINTR53)+MTRAN(:,IMODE,JMODE,CP_BC)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR53 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR53)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR53)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
            END IF ! IF IMODE,JMODE=5,3 (from mode 5 to mode 3)
            IF((IMODE == 5).AND.(JMODE == 4)) THEN
             IF((ICP == CP_BC).AND.(NMASCOAGBCINTR54 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGBCINTR54)=                                 &
       BUD_AER_MAS(:,NMASCOAGBCINTR54)+MTRAN(:,IMODE,JMODE,CP_BC)
             END IF
             IF((ICP == CP_OC).AND.(NMASCOAGOCINTR54 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGOCINTR54)=                                 &
       BUD_AER_MAS(:,NMASCOAGOCINTR54)+MTRAN(:,IMODE,JMODE,CP_OC)
             END IF
            END IF ! IF IMODE,JMODE=5,4 (from mode 5 to mode 4)
            IF((IMODE == 6).AND.(JMODE == 4)) THEN
             IF((ICP == CP_DU).AND.(NMASCOAGDUINTR64 > 0)) THEN
              WHERE(MASK3(:))                                           &
       BUD_AER_MAS(:,NMASCOAGDUINTR64)=                                 &
       BUD_AER_MAS(:,NMASCOAGDUINTR64)+MTRAN(:,IMODE,JMODE,CP_DU)
             END IF
            END IF ! IF IMODE,JMODE=6,4 (from mode 6 to mode 4)
           END IF ! IF MODE(JMODE)
          END DO ! LOOP JMODE (those that may have been transferred
                !             to mode JMODE from mode IMODE)
         END IF ! IF COMPONENT(IMODE,ICP)
        END DO ! LOOP OVER COMPONENTS
!
        ! reset mass if insignificant # of particles
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(.NOT.MASK1(:))                                          &
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
         END IF
        END DO
        WHERE(.NOT.MASK1(:)) MDT(:,IMODE)=MMID(IMODE)
!
       END IF ! IF(MODE(IMODE)
      END DO ! LOOP IMODE

      IF (lhook) CALL dr_hook('UKCA_COAGWITHNUCL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_COAGWITHNUCL
