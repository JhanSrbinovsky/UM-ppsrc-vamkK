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
!    Calculates in-cloud aerosol wet deposition (nucleation-scavenging).
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
      SUBROUTINE UKCA_RAINOUT(NBOX,ND,MD,MDT,FCONV_CONV,DELTA_Z,RHOA,   &
          CRAIN,DRAIN,CRAIN_UP,DRAIN_UP,CLWC,CLF,T,DTC,BUD_AER_MAS,     &
          INUCSCAV,DRYDP,WETDP)
!---------------------------------------------------------------------------
!
!     Calculates in-cloud aerosl wet deposition (nucleation-scavenging)
!
!     Includes large- (dynamic) & small- scale (convective) precip.
!
!     Include conversion rates FCONV_CONV and FCONV_DYN which represent
!     fraction of condensate this is converted to rain in 6 hours.
!
!     Currently takes FCONV_CONV as input with FCONV_DYN set constant
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Initial no. concentration of aer mode (ptcls/cm3)
!     MD         : Initial avg cpt mass conc of aer mode (molcules/ptcl)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     FCONV_CONV : Fraction of box condensate converted to rain in
!                                                   6 hours (convective)
!     CRAIN      : Rain rate for conv. precip. in box (kgm^-2s^-1)
!     DRAIN      : Rain rate for dynamic precip. in box (kgm^-2s^-1)
!     CRAIN_UP   : Rain rate for conv. precip. in box above (kgm^-2s^-1)
!     DRAIN_UP   : Rain rate for dyn. precip. in box above (kgm^-2s^-1)
!     T          : Air temperature (K)
!     DTC        : Chemistry timestep (s)
!     INUCSCAV   : Switch for scheme for removal by nucl scav
!     RHOA       : Air density (kg/m3)
!     DELTA_Z    : Difference between rho_levels, i.e. layer thickness (m)
!     CLF        : Liquid cloud fraction
!     CLWC       : Cloud liquid water content (kg/kg)
!     DRYDP      : Geometric mean dry diameter for each mode (m)
!     WETDP      : Geometric mean wet diameter for each mode (m)
!
!     Outputs
!     -------
!     ND           : Updated no. concentration in each mode (ptcls/cm3)
!     BUD_AER_MASS : Aerosol mass budgets
!
!     Local Variables
!     ---------------
!     SWITCH       : Rain created in box? (0=none,1=conv,2=dyn,3=both=3)
!     RSCAV        : Scavenging parameters for each mode
!     FCONV_DYN    : Fraction of box condensate converted to rain in
!                                                   6 hours (dynamic)
!     TAU_CONV_DYN : e-folding timescale for conversion of
!                                        condensate to rain (dynamic)
!     TAU_CONV_CONV: e-folding timescale for conversion of
!                                        condensate to rain (convective)
!     FBOX_CONV    : Gridbox fraction over which convective rain occurs
!     TICE         : Temperature below which insoluble aerosol can act
!                    as ice nucleii and hence be removed
!     DELN         : Change in number conc. due to nucleation-scavenging
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES       : Number of possible aerosol modes
!     NCP          : Number of possible aerosol components
!     MODE         : Defines which modes are set
!     COMPONENT    : Defines which components are set in each mode
!     MODESOL      : Defines whether mode is soluble of not (integer)
!     NUM_EPS      : Value of NEWN below which don't recalculate MD or
!                                                    carry out process
!     CP_SU        : Index of component in which H2SO4  cpt is stored
!     CP_BC        : Index of component in which BC     cpt is stored
!     CP_OC        : Index of component in which 1st OC cpt is stored
!     CP_CL        : Index of component in which NaCl   cpt is stored
!     CP_DU        : Index of component in which dust   cpt is stored
!     CP_SO        : Index of component in which 2nd OC cpt is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! .. Subroutine interface
      INTEGER, INTENT(IN) :: NBOX                    ! Number of grid boxes
      INTEGER, INTENT(IN) :: INUCSCAV                ! Switch for scheme for
!                                                      removal by nucl scav
      REAL, INTENT(IN)    :: FCONV_CONV(NBOX)
      REAL, INTENT(IN)    :: DELTA_Z(NBOX)
      REAL, INTENT(IN)    :: RHOA(NBOX)
      REAL, INTENT(IN)    :: CRAIN(NBOX)
      REAL, INTENT(IN)    :: DRAIN(NBOX)
      REAL, INTENT(IN)    :: CRAIN_UP(NBOX)
      REAL, INTENT(IN)    :: DRAIN_UP(NBOX)
      REAL, INTENT(IN)    :: CLWC(NBOX)
      REAL, INTENT(IN)    :: CLF(NBOX)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: DTC
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(INOUT) :: ND (NBOX,NMODES)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD (NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)

! .. Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: SWITCH
      REAL    :: RSCAV(NMODES)
      REAL    :: RSCAVM(NMODES)
      REAL    :: DELN
      REAL    :: DELN0
      REAL    :: NDNEW
      REAL    :: TAU_CONV_CONV
      REAL    :: TAU_CONV_DYN
      REAL    :: BETA
      REAL    :: CRAINDIFF
      REAL    :: DRAINDIFF
      REAL, PARAMETER :: FCONV_DYN=0.9999
      REAL, PARAMETER :: FBOX_CONV=0.3
      REAL, PARAMETER :: TICE=258.0

      REAL    :: LNRATN
      REAL    :: ERFNUM
      REAL    :: FRAC_N
      REAL    :: DERF
      REAL    :: LOG2SG
      REAL    :: LNRATM
      REAL    :: ERFMAS
      REAL    :: FRAC_M
      REAL    :: DM(NCP)
      REAL    :: SCAV1
      REAL    :: SCAV2
      REAL    :: SCAV1M
      REAL    :: SCAV2M

      REAL    :: DP
      REAL    :: DP2
      REAL    :: MAXFRACN
      REAL    :: NSCAVACT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_RAINOUT',zhook_in,zhook_handle)

!      print*,'UKCA_RAINOUT: min,max of CLWC=',minval(CLWC),maxval(CLWC)
!
!!      NSCAVACT=150.0e-9 ! dry radius
      NSCAVACT=103.0e-9 ! dry radius
!!      NSCAVACT=70.0e-9 ! dry radius
!!      NSCAVACT=40.0e-9 ! dry radius
!!      NSCAVACT=20.0e-9 ! dry radius
!
      MAXFRACN=0.90
!
! ..  Calculate timescale for conversion to dynamic rain assuming
! ..  a factor FCONV_DYN is converted to rain in 6 hours
      TAU_CONV_DYN=(-6.0*3600.0)/(LOG(1.0-FCONV_DYN))

!! now follow approach implemented at start of AEROS to better mimic
!! GLOMAP-bin nucleation-scavenging approach with size threshold.

      DO JL=1,NBOX
!
       IF(INUCSCAV == 1) THEN


        IMODE=2

!!        DP=DRYDP(JL,IMODE)
        DP=WETDP(JL,IMODE)
!!        print*,'Aitken-soluble 0.5*DP,NSCAVACT=',(0.5*DP),NSCAVACT


         IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
          LNRATN=LOG(NSCAVACT*2.0/DP)
          ERFNUM=LNRATN/SQRT(2.0)/LOG(SIGMAG(IMODE))
          FRAC_N=0.5*(1.0+DERF(dble(ERFNUM))) ! fraction remaining in mode

          IF(FRAC_N > MAXFRACN) THEN ! if less than 1% of number is larger than size then don't remove
!           print*,'Aitsol:  <1% larger than NSCAVACT. Dont do NUCSCAV'
           SCAV1 =0.0
           SCAV1M=0.0
          ELSE
!           print*,'Aitsol: >=1% larger than NSCAVACT so do NUCSCAV'
           SCAV1 =1.0-FRAC_N

           LOG2SG=LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE))
           DP2=EXP(LOG(DP)+3.0*LOG2SG) ! volume median diameter
!!           print*,'Aitken-soluble 0.5*DP2,NSCAVACT=',(0.5*DP2),NSCAVACT

! .. 2nd threshold determines fraction of number/mass
           LNRATM=LOG(NSCAVACT*2.0/DP2)
           ERFMAS=LNRATM/SQRT(2.0)/LOG(SIGMAG(IMODE))
           FRAC_M=0.5*(1.0+DERF(dble(ERFMAS)))
           SCAV1M=1.0-FRAC_M
          END IF
!!          print*,'Aitken-soluble nmbr: RSCAV =',SCAV1
!!          print*,'Aitken-soluble mass: RSCAVM=',SCAV1M
         ELSE
           SCAV1 =0.0
           SCAV1M=0.0
         END IF  ! IF(ND(JL,IMODE) > NUM_EPS(IMODE))

        IMODE=3

!!        DP=DRYDP(JL,IMODE)
        DP=WETDP(JL,IMODE)
!!        print*,'accumn-soluble 0.5*DP,NSCAVACT=',(0.5*DP),NSCAVACT


         IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
          LNRATN=LOG(NSCAVACT*2.0/DP)
          ERFNUM=LNRATN/SQRT(2.0)/LOG(SIGMAG(IMODE))
          FRAC_N=0.5*(1.0+DERF(dble(ERFNUM))) ! fraction remaining in mode
          IF(FRAC_N > MAXFRACN) THEN ! if less than 10% of number is larger than size then don't remove
!!           print*,'acc-sol:  <1% larger than NSCAVACT. Dont do NUCSCAV'
           SCAV2 =0.0
           SCAV2M=0.0
          ELSE
!!           print*,'acc-sol: >=1% larger than NSCAVACT so do NUCSCAV'
           SCAV2 =1.0-FRAC_N

           LOG2SG=LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE))
           DP2=EXP(LOG(DP)+3.0*LOG2SG) ! volume median diameter
!!           print*,'accumn-soluble 0.5*DP2,NSCAVACT=',(0.5*DP2),NSCAVACT

! .. 2nd threshold determines fraction of number/mass
           LNRATM=LOG(NSCAVACT*2.0/DP2)
           ERFMAS=LNRATM/SQRT(2.0)/LOG(SIGMAG(IMODE))
           FRAC_M=0.5*(1.0+DERF(dble(ERFMAS)))
! .. limit DELM to be max 99.9%
           IF(FRAC_M < 0.001) FRAC_M=0.001

           SCAV2M=1.0-FRAC_M
          END IF
!!          print*,'accumn-soluble nmbr: RSCAV =',SCAV2
!!          print*,'accumn-soluble mass: RSCAV =',SCAV2M

         ELSE
           SCAV2 =0.0
           SCAV2M=0.0
         END IF  ! IF(ND(JL,IMODE) > NUM_EPS(IMODE))

        RSCAV =(/0.00,SCAV1 ,SCAV2 ,1.00,1.00,1.00,1.00/)
        RSCAVM=(/0.00,SCAV1M,SCAV2M,1.00,1.00,1.00,1.00/)
!!        print*,'RSCAV(2),RSCAV(3)=',RSCAV(2),RSCAV(3)
!!        print*,'RSCAVM(2),RSCAVM(3)=',RSCAVM(2),RSCAVM(3)
!!        pause

       END IF ! INUCSCAV=1

       IF(INUCSCAV == 2) THEN ! set number and mass same as M7
        RSCAV =(/0.10,0.25,0.85,0.99,0.20,0.40,0.40/)
        RSCAVM=(/0.10,0.25,0.85,0.99,0.20,0.40,0.40/)
       END IF

       SWITCH=0
       CRAINDIFF=CRAIN(JL)-CRAIN_UP(JL)
       DRAINDIFF=DRAIN(JL)-DRAIN_UP(JL)
       IF(CRAINDIFF > 1.0e-10) SWITCH=1
       IF(DRAINDIFF > 1.0e-10) THEN
        IF(SWITCH == 0) SWITCH=2
        IF(SWITCH == 1) SWITCH=3
       END IF

! .. Switch > 0 only if rain is FORMED in that level
! .. Convective rain = 1, Dynamic rain = 2
! .. No rain occurs at the top level of the atm
!
       IF(SWITCH > 0)THEN

        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN
          IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN

! .. Only apply to soluble modes or insoluble when T<TICE
           IF((MODESOL(IMODE) == 1).OR.(T(JL) < TICE)) THEN

!----------------------------------------------------------------------
!
! This section does removal by small-scale (conv) precipitation
! .. Convective rain : ND -> ND*(1-FCONV_CONV) over 6 hours
!                      only apply over fraction FBOX_CONV
!
            IF(SWITCH == 1) THEN ! if convective rain
             IF(FCONV_CONV(JL) < 1.0) THEN
              TAU_CONV_CONV=(-6.0*3600.0)/(LOG(1.0-FCONV_CONV(JL)))
              DELN0=FBOX_CONV*ND(JL,IMODE)*(1.0-EXP(-DTC/TAU_CONV_CONV))
             ELSE
              DELN0=FBOX_CONV*ND(JL,IMODE)
             END IF
            END IF ! if small-scale (convective) rain

!-----------------------------------------------------------------------
!
! This section does removal by large-scale (dyn.) precipitation
!
! .. Dynamic rain    : ND -> ND*(1-FCONV_DYN) over 6 hours
!                      apply over all of box
!
            IF((SWITCH == 2).OR.(SWITCH == 3)) THEN ! if dynamic rain
             IF ( CLWC(JL) < 1.0e-10 ) THEN
! Old method
              DELN0=ND(JL,IMODE)*(1.0-EXP(-DTC/TAU_CONV_DYN ))
             ELSE
              BETA = DRAINDIFF/DELTA_Z(JL)/RHOA(JL)/CLWC(JL)
!! here, CLWC needs to match units of DRAIN/CRAIN (kg/kg since rain in kg/m2/s)
!! if CLWC is in g/kg then divide by (CLWC(JL)*1.0e-3) rather than CLWC(JL)
!! --- POSSIBLE UNITS MIS_MATCH HERE!!!!! 
!! -- if currently have wrong units, BETA would be a factor 1000 too large
!! -- i.e. wet removal would be a factor 1000 too intense.
!!              DELN0=BETA*CLF(JL)*ND(JL,IMODE)*DTC
              DELN0=ND(JL,IMODE)*CLF(JL)*(1.0-EXP(-DTC*BETA))
!! updated DELN as CLF*(1.0-EXP(-DTC*BETA)) rather than BETA*ND*DTC*CLF
             END IF 

            END IF ! if large-scale (dynamic) rain

!-----------------------------------------------------------------------
!
! multiply DELN0 by scavenging coefficient for delta for number (DELN)
! set DM for each aerosol component based on RSCAVM value
! Note, if INUCSCAV=1, then  RSCAV & RSCAVM are calculated based on size

            DELN=DELN0*RSCAV(IMODE)
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
! .. calculate cpt mass conc to transfer to next mode (use old no/mass)
              DM(ICP)=DELN0*MD(JL,IMODE,ICP)*RSCAVM(IMODE)
             END IF
            END DO

!----------------------------------------------------------------------

! .. calculate updated number concentration due to nucleation-scavening
            NDNEW=ND(JL,IMODE)-DELN

            IF(NDNEW > NUM_EPS(IMODE)) THEN

! .. first remove mass to be transferred from mode IMODE
             MDT(JL,IMODE)=0.0
             DO ICP=1,NCP
              IF(COMPONENT(IMODE,ICP)) THEN
               MD(JL,IMODE,ICP)=                                         &
                   (ND(JL,IMODE)*MD(JL,IMODE,ICP)-DM(ICP))/NDNEW
               MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
              ELSE
               MD(JL,IMODE,ICP)=0.0
              END IF ! COMPONENT(IMODE,ICP)
             END DO

! .. update number concentration following nucleation-scavening
             ND(JL,IMODE)=NDNEW

!-----------------------------------------------------------------------
!
! .. This section stores removal of each cpt mass for budget calculations

             DO ICP=1,NCP
              IF(COMPONENT(IMODE,ICP)) THEN
               IF(ICP == CP_SU) THEN
                IF((IMODE == 1).AND.(NMASNUSCSUNUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUNUCSOL)+DM(ICP)                 
                IF((IMODE == 2).AND.(NMASNUSCSUAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUAITSOL)+DM(ICP)                 
                IF((IMODE == 3).AND.(NMASNUSCSUACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCSUCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUCORSOL)+DM(ICP)                 
               END IF
               IF(ICP == CP_BC) THEN
                IF((IMODE == 2).AND.(NMASNUSCBCAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCAITSOL)+DM(ICP)                 
                IF((IMODE == 3).AND.(NMASNUSCBCACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCBCCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCCORSOL)+DM(ICP)                 
                IF((IMODE == 5).AND.(NMASNUSCBCAITINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCAITINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCAITINS)+DM(ICP)                 
               END IF
               IF(ICP == CP_OC) THEN
                IF((IMODE == 1).AND.(NMASNUSCOCNUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCNUCSOL)+DM(ICP)                 
                IF((IMODE == 2).AND.(NMASNUSCOCAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCAITSOL)+DM(ICP)                 
                IF((IMODE == 3).AND.(NMASNUSCOCACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCOCCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCCORSOL)+DM(ICP)                 
                IF((IMODE == 5).AND.(NMASNUSCOCAITINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCAITINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCAITINS)+DM(ICP)                 
               END IF
               IF(ICP == CP_CL) THEN
                IF((IMODE == 3).AND.(NMASNUSCSSACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSSACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSSACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCSSCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSSCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSSCORSOL)+DM(ICP)                 
               END IF
               IF(ICP == CP_SO) THEN
                IF((IMODE == 1).AND.(NMASNUSCSONUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSONUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSONUCSOL)+DM(ICP)                 
                IF((IMODE == 2).AND.(NMASNUSCSOAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOAITSOL)+DM(ICP)                 
                IF((IMODE == 3).AND.(NMASNUSCSOACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCSOCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOCORSOL)+DM(ICP)                 
               END IF
               IF(ICP == CP_DU) THEN
                IF((IMODE == 3).AND.(NMASNUSCDUACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCDUACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCDUACCSOL)+DM(ICP)                 
                IF((IMODE == 4).AND.(NMASNUSCDUCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCDUCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCDUCORSOL)+DM(ICP)                 
                IF((IMODE == 6).AND.(NMASNUSCDUACCINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCDUACCINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCDUACCINS)+DM(ICP)                 
                IF((IMODE == 7).AND.(NMASNUSCDUCORINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCDUCORINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCDUCORINS)+DM(ICP)                 
               END IF

              END IF ! if component(imode,icp)
             END DO ! loop over components

            END IF ! if ND>NUM_EPS(IMODE)

!----------------------------------------------------------------------

           END IF ! if mode is soluble or T<TICE

          END IF ! IF ND>NUM_EPS
         END IF ! IF MODE is switched on
        END DO ! Loop over modes

       END IF ! IF RAIN IS PRODUCED IN THIS LEVEL

      END DO ! Loop over gridboxes

      IF (lhook) CALL dr_hook('UKCA_RAINOUT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_RAINOUT
