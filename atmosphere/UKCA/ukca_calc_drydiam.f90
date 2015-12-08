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
!   Calculates geometric mean dry diameter for multi-component
!   aerosol population which is lognormally distributed with
!   number concentration ND in mode, component mass concentration
!   MD in mode, component density RHOCOMP and component molecular mass MM
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
      SUBROUTINE UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)

!-------------------------------------------------------------
!
!     Calculates geometric mean dry diameter for multi-component
!     aerosol population which is lognormally distributed with
!     number concentration ND in mode,
!     cpt mass concentration MD in mode,
!     cpt density RHOCOMP and cpt molecular mass MM
!
!     Calculate dry volume per particle using composition info as:
!
!     DVOL = sum_cp { MD(ICP)*MM(ICP)/(AVC*RHOCOMP(ICP)) }
!
!     Where AVC is Avogadro's constant. Then, from Jacobsen,
!     "Fundamentals of Atmospheric Modeling", pp. 412, have
!
!     dry volume conc. = ND*(PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     i.e. DVOL  = (PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     where Dp is the number mean dry diameter,
!     and sigma_g is the geometric standard deviation.
!
!     Then calcaulte Dp as:
!
!     Dp=CBRT( DVOL*(6/PI)/EXP{9.2 log^2(sigma_g)} )
!
!     Inputs
!     ------
!     NBOX      : Number of grid boxes
!     ND        : Aerosol ptcl no. concentration (ptcls per cc)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!     VERBOSE   : Switch for level of verbosity
!
!     Outputs
!     -------
!     DRYDP     : Median particle dry diameter for each mode (m)
!     DVOL      : Median particle dry volume for each mode (m^3)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!
!     Local Variables
!     ---------------
!     None
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI       : 3.1415927
!     AVC       : Avogadro's constant (per mole)
!     MMSUL     : Molar mass of a pure H2SO4 aerosol (kg per mole)
!     RHOSUL    : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES    : Number of aerosol modes
!     NCP       : Number of aerosol components
!     MM        : Molecular mass of each component
!     RHOCOMP   : Densities (dry) of each component (kg/m^3)
!     MODE      : Which modes are being carried
!     COMPONENT : Which components are in each of modes
!     DDPLIM0   : Lower limits for dry diameter for each mode (m)
!     DDPLIM1   : Upper limits for dry diameter for each mode (m)
!     MFRAC_0   : Initial mass fraction to set when no particles.
!     X         : EXP((9/2)*LOG^2(SIGMA_G))
!     NUM_EPS   : Value of NEWN below which do not carry out process
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE yomhook,     ONLY: lhook, dr_hook
      USE parkind1,    ONLY: jprb, jpim
      USE ereport_mod, ONLY: ereport
      IMPLICIT NONE

! Interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: VERBOSE
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(OUT)   :: DRYDP(NBOX,NMODES)
      REAL, INTENT(OUT)   :: DVOL(NBOX,NMODES)

! Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      REAL    :: DDPCUB(NBOX)
      REAL    :: SIXOVRPIX(NMODES)
      LOGICAL :: MASK(NBOX)
      REAL    :: RATIO1(NCP)
      REAL    :: RATIO2(NMODES)
      REAL    :: DP_THRESH1
      REAL    :: DP
      CHARACTER(LEN=72) :: cmessage
      INTEGER           :: errcode

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
! Below is over NCP
      IF (lhook) CALL dr_hook('UKCA_CALC_DRYDIAM',zhook_in,zhook_handle)
      RATIO1(:)=MM(:)/AVC/RHOCOMP(:)
!
! Below is over NMODES
      SIXOVRPIX(:)=6.0/(PPI*X(:))
      RATIO2(:)=MMSUL*MMID(:)/(AVC*RHOSUL)
!
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        MASK(:) =(ND(:,IMODE) > NUM_EPS(IMODE))
        WHERE(MASK(:))
         DVOL(:,IMODE)=0.0
        ELSEWHERE
         DVOL(:,IMODE)=MMID(IMODE)*MMSUL/(AVC*RHOSUL)
        ENDWHERE
!     calculate particle dry volume using composition info
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(MASK(:))                                                &
           DVOL(:,IMODE)=DVOL(:,IMODE)+RATIO1(ICP)*MD(:,IMODE,ICP)
         END IF
        END DO
!     DVOL calculates particle dry volume assuming pure H2SO4
       ELSE
        DVOL(:,IMODE)=RATIO2(IMODE)
       END IF
      END DO

      DO IMODE=1,NMODES
       DDPCUB(:)=SIXOVRPIX(IMODE)*DVOL(:,IMODE)
       DRYDP(:,IMODE)=DDPCUB(:)**(1.0/3.0)
      END DO

! .. also check whether mean diameter too low for mode
!      DO IMODE=1,NMODES
      DO IMODE=1,3 ! only do check for modes 1, 2 and 3
!      DO IMODE=1,1 ! only do check for mode 1
       IF(MODE(IMODE)) THEN
        DO JL=1,NBOX
         DP    =DRYDP(JL,IMODE)
         DP_THRESH1=DDPLIM0(IMODE)*0.1
         IF(DP < DP_THRESH1) THEN
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            MD(JL,IMODE,ICP)=MLO(IMODE)*MFRAC_0(IMODE,ICP)
           ENDIF
          ENDDO
          MDT(JL,IMODE)=MLO(IMODE)
          DVOL(JL,IMODE)=MLO(IMODE)*MMSUL/(AVC*RHOSUL)
          DRYDP(JL,IMODE)=(SIXOVRPIX(IMODE)*DVOL(JL,IMODE))**(1.0/3.0)
         END IF ! if DP<DPTHRESH
        END DO
       END IF
      END DO

      DO IMODE=1,NMODES
       IF((minval(dvol (:,IMODE)) <= 0.0) .OR.                          &
          (minval(drydp(:,IMODE)) <= 0.0)) THEN
        cmessage = ' dvol or drydp <= 0'
        WRITE(6,*) 'In calcdrydiam: drydp (min,max,sum) imode=',imode
        WRITE(6,*) MINVAL(drydp(:,IMODE)),MAXVAL(drydp(:,IMODE)),       &
                   SUM(drydp(:,IMODE))
        WRITE(6,*) 'Location of min: ',MINLOC(drydp(:,IMODE))
        WRITE(6,*) 'Location of max: ',MAXLOC(drydp(:,IMODE))
        WRITE(6,*) 'In calcdrydiam: dvol (min,max,sum) imode=',imode
        WRITE(6,*) MINVAL(dvol(:,IMODE)),MAXVAL(dvol(:,IMODE)),         &
                   SUM(dvol(:,IMODE))
        WRITE(6,*) 'Location of min: ',MINLOC(drydp(:,IMODE))
        WRITE(6,*) 'Location of max: ',MAXLOC(drydp(:,IMODE))
        errcode=1

        CALL EREPORT('ukca_calc_drydiam',errcode,cmessage)
       END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_CALC_DRYDIAM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_DRYDIAM
