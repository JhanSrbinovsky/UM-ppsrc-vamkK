! *****************************COPYRIGHT*******************************
! (c) [University of Leeds] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! *****************************COPYRIGHT*******************************
!
!  Description:
!     UKCA-MODE aerosol code: interface routine called from
!     UKCA_AERO_CTL to compute first order reaction coefficients
!     for loss to aerosol for use in the chemistry solver.
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

      SUBROUTINE UKCA_TROP_HETCHEM(nbox, nhet, temp, rh, aird, pvol,    &
                                   wetdp, sarea, het_rates)

      USE UKCA_CONSTANTS,       ONLY: avc, rmol, m_ho2, m_n2o5,         &
                                      m_air
      USE UKCA_MODE_SETUP,      ONLY: nmodes, ncp, mode, component,     &
                                      cp_su, cp_bc, cp_oc, cp_cl,       &
                                      cp_du, cp_so
      USE conversions_mod,      ONLY: pi
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE ereport_mod,          ONLY: ereport
      USE Control_Max_Sizes
      USE PrintStatus_mod
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

! Arguments
      INTEGER, INTENT(IN) :: nbox                 ! No of points
      INTEGER, INTENT(IN) :: nhet                 ! No het. chem. rates

      REAL, INTENT(IN)    :: temp(nbox)           ! temperature [K]
      REAL, INTENT(IN)    :: rh(nbox)             ! Relative humidity [fraction]
      REAL, INTENT(IN)    :: aird(nbox)           ! Total number density [molecules/cm^3]
      REAL, INTENT(IN)    :: pvol(nbox,nmodes,ncp) 
! .. Partial volume of soluble components as fraction of
! .. solution volume according to fraction of total mass (including water)
      REAL, INTENT(IN)    :: wetdp(nbox,nmodes)   ! Mean wet radius for the mode [m]
      REAL, INTENT(IN)    :: sarea(nbox,nmodes)   ! Surface area concentn. [cm^2/cm^3]
      REAL, INTENT(OUT)   :: het_rates(nbox,nhet) ! rate coefficients

! Local variables

      INTEGER :: i, j, k, l                   ! Loop variables
      INTEGER :: imode                        ! mode loop counter
      INTEGER :: icp                          ! component loop counter

      REAL    :: molec_weight(nhet)           ! Molecular weights of species
      REAL    :: sticking_cf(nbox,nhet,ncp)   ! Gamma as calculated by CAL_* [unitless]
      REAL    :: total_vol(nbox,nmodes)       ! Total volume of each mode
      REAL    :: frac_vol(nbox,nmodes,ncp)    ! pvol/total_vol for each mode [unitless]
      REAL    :: mode_gamma(nbox,nhet,nmodes) ! Weighted gammas for each mode
      REAL    :: mode_rate(nbox,nhet,nmodes)  ! First order loss rate coefficient
      REAL    :: total_sa(nbox)               ! Total aerosol surface area

      INTEGER           :: errcode            ! Error code
      CHARACTER(LEN=72) :: cmessage           ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM',zhook_in,zhook_handle)

      IF (PrintStatus >= PrStatus_Diag) THEN 
        WRITE(6,*) 'UKCA_TROP_HETCHEM Routine:'
        WRITE(6,*) 'nbox: ',nbox
        WRITE(6,*) 'temp: ',maxval(temp),minval(temp),sum(temp)/        &
                    size(temp)
        WRITE(6,*) 'rh  : ',maxval(rh),minval(rh),sum(rh)/size(rh)
        WRITE(6,*) 'aird: ',maxval(aird),minval(aird),sum(aird)/        &
                    size(aird)
        DO imode=1,nmodes
          IF (mode(imode)) THEN
            WRITE(6,*) 'sarea: ',imode,maxval(sarea(:,imode)),          &
                        minval(sarea(:,imode)),                         &
                        sum(sarea(:,imode))/size(sarea(:,1))
            WRITE(6,*) 'wetdp: ',imode,maxval(wetdp(:,imode)),          &
                        minval(wetdp(:,imode)),                         &
                        sum(wetdp(:,imode))/size(wetdp(:,1))
          END IF
        END DO

        DO imode=1,nmodes
          IF (mode(imode)) THEN
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
               WRITE(6,*) 'pvol: ',imode,icp,maxval(pvol(:,imode,icp)), &
                           minval(pvol(:,imode,icp)),                   &
                           sum(pvol(:,imode,icp))/size(pvol(:,1,1))
              END IF
            END DO
          END IF
        END DO
      END IF

      molec_weight(:) = 9e5
      molec_weight(1) = m_n2o5
      molec_weight(2) = m_ho2
      IF (ANY(molec_weight > 1e5)) THEN
        cmessage = ' Molecular weights not defined for all cases'
        errcode = 1
        CALL EREPORT('UKCA_TROP_HETCHEM',errcode,cmessage)
      END IF

! ========================================================
!  1. Calculate gamma on each aerosol type for 1,nhet
! ========================================================

! Get <gamma> for N2O5 or HO2 reaction, which is a function of aerosol type,
! temperature, and RH.

      IF (PrintStatus >= PrStatus_Diag) WRITE(6,'(A)') 'Sticking_cf:'
      sticking_cf(:,:,:) = 0.0e0
      DO i=1,nhet
        IF ( i == 1 ) THEN
          DO icp=1,ncp      ! No of components
            sticking_cf(:,i,icp) = CAL_N2O5(icp, temp(:), rh(:))
            IF (PrintStatus >= PrStatus_Diag)                           &
            WRITE(6,*) i,icp,maxval(sticking_cf(:,i,icp)),              &
                             minval(sticking_cf(:,i,icp)),              &
                   sum(sticking_cf(:,i,icp))/size(sticking_cf(:,i,icp))
          END DO
        ELSE IF ( i == 2 ) THEN 
          DO icp=1,ncp
            sticking_cf(:,i,icp) = CAL_HO2(icp, temp(:), rh(:))
            IF (PrintStatus >= PrStatus_Diag)                           &
            WRITE(6,*) i,icp,maxval(sticking_cf(:,i,icp))
          END DO
        ELSE
          cmessage = ' Error in calculating sticking_cf' 
          errcode = 1
          CALL EREPORT('UKCA_TROP_HETCHEM',errcode,cmessage)
        END IF
      END DO   ! end of loop over nhet

! So now we have an array sticking_cf(nbox,nhet,aerotypes)
!  with gamma for n2o5 and ho2 for each aerosol type.

! =================================================================
!  2. Calculate mean gamma for each modes based on partial volumes
! =================================================================

      ! Calc fractional volumes

      total_vol(:,:) = 0.0e0
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          DO icp=1,ncp
            IF (component(imode,icp)) THEN
               total_vol(:,imode) = total_vol(:,imode) +                &
                                   pvol(:,imode,icp)
            END IF
          END DO  ! end of loop over components
        END IF
      END DO  ! end of loop over modes

      IF (PrintStatus >= PrStatus_Diag) WRITE(6,'(A)') 'frac_vol:' 
      frac_vol(:,:,:) = 0.0e0
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          DO icp=1,ncp
            IF (component(imode,icp)) THEN
              frac_vol(:,imode,icp) = pvol(:,imode,icp) /               &
                                           total_vol(:,imode)
              IF (PrintStatus >= PrStatus_Diag)                         &
              WRITE(6,*) imode,icp,maxval(frac_vol(:,imode,icp)),       &
                                   minval(frac_vol(:,imode,icp)),       &
                 sum(frac_vol(:,imode,icp))/size(frac_vol(:,imode,icp))
            END IF
          END DO  ! end of loop over components
        END IF
      END DO  ! end of loop over modes

! Calculate mode_gamma by weighting according to frac_vols,
!   e.g. gammaSO4*frac_volSO4, for each component, then sum up over each mode.

      mode_gamma(:,:,:) = 0.0
      DO i=1,nhet
        DO imode=1,nmodes
          IF (mode(imode)) THEN
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                mode_gamma(:,i,imode) = mode_gamma(:,i,imode) +         &
                    (frac_vol(:,imode,icp) * sticking_cf(:,i,icp))
              END IF
            END DO ! end loop over ncp
          END IF
        END DO ! end loop over nmodes
      END DO  ! end loop over nhet

! So now we have mode_gamma(nbox,nhet,nmodes)

! ===========================================================
! 3. Calculate first order rate coefficients for each mode
! ===========================================================

      IF (PrintStatus >= PrStatus_Diag) WRITE(6,'(A)') 'Mode_rate:'
      mode_rate(:,:,:) = 0.0e0
      DO i=1,nhet
        DO imode=1,nmodes
          IF (mode(imode)) THEN
            mode_rate(:,i,imode) = CAL_LOSS1K(sarea(:,imode),           &
                                   wetdp(:,imode), aird(:),             &
                                   mode_gamma(:,i,imode), temp(:),      &
                                   molec_weight(i))
            IF (PrintStatus >= PrStatus_Diag)                           &
              WRITE(6,*) imode,maxval(mode_rate(:,i,imode)),            &
                                 minval(mode_rate(:,i,imode)),          &
                  sum(mode_rate(:,i,imode))/size(mode_rate(:,i,imode))
          ENDIF
        END DO
      END DO

! =====================================================
! 4. Calc average rate for the box by adding mode rates
! =====================================================


      total_sa(:) = 0.0e0
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          total_sa(:) = total_sa(:) + sarea(:,imode)
        END IF
      END DO

      het_rates(:,:)=0.0e0
      DO i=1,nhet
        DO imode=1,nmodes
          IF (mode(imode)) THEN
            WHERE (total_sa > 1e-30)
              het_rates(:,i) = het_rates(:,i) + mode_rate(:,i,imode)
            ENDWHERE
          END IF
        END DO
      END DO


      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(6,*) 'rc_n2o5: ',maxval(het_rates(:,1)),                  &
                             minval(het_rates(:,1)),                    &
                             sum(het_rates(:,1))/size(het_rates(:,1))
        WRITE(6,*) 'rc_ho2:  ',maxval(het_rates(:,2)),                  &
                             minval(het_rates(:,2)),                    &
                             sum(het_rates(:,2))/size(het_rates(:,2))
      END IF

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM',zhook_out,zhook_handle)

      RETURN

! ######################################################################
!                       INTERNAL SUBROUTINES 
      CONTAINS
! ===============================================================
!                N 2 O 5   F U N C T I O N
! ===============================================================

      !==============================================================
      ! Function N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis.
      !
      ! NOTES:
      !==============================================================

      FUNCTION CAL_N2O5(aerotype, temp, rh)

      !==============================================================
      ! Internal function CAL_N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis.
      !==============================================================
      IMPLICIT NONE
      
      ! Arguments
      INTEGER, INTENT(IN) :: aerotype     ! # denoting aerosol cpt
      REAL,    INTENT(IN) :: temp(nbox)   ! Temperature [K]
      REAL,    INTENT(IN) :: rh(nbox)     ! Relative Humidity [fraction]

      ! Function return value
      REAL :: cal_n2o5(nbox)

      ! Local variables
      REAL                :: ttemp(nbox)   ! dummy array for temp
      REAL                :: rh_p(nbox)    ! RH as a percent
      REAL                :: maxtemp       ! maximum

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      !==============================================================
      ! N2O5 begins here!
      !==============================================================

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_N2O5',zhook_in,  &
                              zhook_handle)

      ! Convert RH to % (max = 100%)
! Values already limited to <= 1.0 - see rel_humid_frac in UKCA_MAIN1
      rh_p(:)  = rh(:) * 100.0

      ! Special handling for various aerosols
      SELECT CASE ( aerotype )

         !----------------
         ! Dust
         !----------------
         CASE ( cp_du )        
          
            !======================================================
            ! Based on unpublished Crowley work, was 0.01E0
            !------------------------------------------------------
            ! Now gamma=0.03-0.15 (Mogili et al 2006) so use
            ! a value of 0.1.
            ! Also 0.2 on saharan test dust (Karagulian et al 2006)
            ! though this is now thought to be too high (method)
            ! (hlm 2008)
            !======================================================
            cal_n2o5(:) = 0.1

         !----------------
         ! Sulfate
         !----------------
         CASE ( cp_su )      

            !========================================================
            ! RH dependence from Kane et al., Heterogeneous uptake of 
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470 
            !========================================================
            cal_n2o5(:) = 2.79e-4 + rh_p(:)*(  1.30e-4 +                &
                                 rh_p(:)*( -3.43e-6 +                   &
                                 rh_p(:)*(  7.52e-8 ) ) )

            !========================================================
            ! Temperature dependence factor (Cox et al, Cambridge UK) 
            ! is of the form:
            !
            !          10^( LOG10( G294 ) - 0.04 * ( ttemp - 294 ) )
            ! fact = ------------------------------------------------
            !                     10^( LOG10( G294 ) )
            !
            ! Where G294 = 1e-2 and ttemp is MAX( temp, 282 ).
            ! 
            ! For computational speed, replace LOG10( 1e-2 ) with -2
            ! and replace 10^( LOG10( G294 ) ) with G294 
            !========================================================
!            ttemp(:) = MAX( temp(:), 282E0 ) - 294.0  [MAX isn't array valued]
            maxtemp = 282.0 - 294.0
            WHERE (temp < 282e0)
             ttemp = maxtemp
            ELSEWHERE
             ttemp = temp - 294.0
            ENDWHERE

            ! Apply temperature dependence
            cal_n2o5(:) = cal_n2o5(:)*10.0**(-2.0-4e-2*ttemp(:))/1e-2
                                
         !----------------
         ! Black Carbon
         !----------------
         CASE ( cp_bc )

            !======================================================
            ! From IUPAC
            !======================================================
            cal_n2o5(:) = 0.005e0

         !----------------
         ! Organic Carbon
         !----------------
         CASE ( cp_oc , cp_so )  

            !========================================================
            ! Based on Thornton, Braban and Abbatt, 2003
            ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
            ! of relative humidity, particle phase and particle size
            !========================================================
            WHERE ( rh_p >= 57E0 )
               cal_n2o5(:) = 0.03e0
            ELSEWHERE
               cal_n2o5(:) = rh_p(:) * 5.2e-4
            ENDWHERE

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------           
         CASE ( cp_cl )          

            !=======================================================
            ! Based on IUPAC recommendation
            !-------------------------------------------------------
            ! May'07 lechlm Now ELSE is changed from G=0.005E0 to 
            ! gamma = rh_p * 0.0005E0 with a study by Thornton and
            ! Abbatt 2005, JOURNAL.
            !=======================================================
            WHERE ( rh_p >= 62 )
               cal_n2o5(:) = 0.03e0
            ELSEWHERE
               cal_n2o5(:) = rh_p(:) * 0.0005e0
            ENDWHERE

         !----------------         
         ! Default
         !----------------
         CASE DEFAULT
            cmessage=' Unknown aerosol surface for N2O5 hydrolysis'
            WRITE (6,*) cmessage
            WRITE (6,*) 'AEROSOL TYPE =',aerotype
            CALL EREPORT('UKCA_TROP_HETCHEM:CAL_N2O5',aerotype,cmessage)

      END SELECT   
         
      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_N2O5',zhook_out,   &
                              zhook_handle)

      RETURN
      END FUNCTION CAL_N2O5

! ======================================================================
!                  H O 2   F U N C T I O N
! ======================================================================

      FUNCTION CAL_HO2(aerotype, temp, rh)

      !================================================================
      ! Function HO2 computes the gamma sticking factor
      ! for HO2 self reaction.
      ! 
      ! NOTES:
      ! (1 ) Dust values for RH<50% are needed.
      ! (2 ) Data for organics is from a conference poster that was not
      !        published. Replace with data as it becomes available.
      ! (3 ) Temperature dependence on salt seems fairly good (N.B. this 
      !        is only for SOLID NaCl, not aq). Fits available data 
      !        (though this data was not published to our knowledge).
      !        This temperature dependence fits well for sulphate data,
      !        so has been used here due to lack of data. 
      !        Update when available.
      ! (4 ) Data for various relative humidities available now (04/08)
      !        from Taketani et al 2008 on sulfate aerosol.
      ! (5 ) It should be noted that it may be possible for gamma to 
      !        reach very high values if TMIs are present. However,
      !        it is unclear if there are any observations of this.
      !        See notes on dust in the code further down.
      ! (hlm 3/3/09)
      !================================================================

      IMPLICIT NONE
      
      ! Arguments
      INTEGER, INTENT(IN) :: aerotype         ! # denoting aerosol cpt
      REAL,    INTENT(IN) :: temp(nbox)     ! Temperature [K]
      REAL,    INTENT(IN) :: rh(nbox)       ! Relative Humidity [fraction]

      ! Function return value
      REAL :: cal_ho2(nbox)

      ! Local variables
      REAL                :: rh_p(nbox)     ! dummy array for RH
      REAL                :: factt(nbox)    ! for calc temp dependence
      REAL                :: facth(nbox)    ! for calc rh dependence

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!==============================================================
! HO2 begins here!
!==============================================================

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_HO2',zhook_in,  &
                               zhook_handle)

      ! Convert rh to % (max = 100%)  !! RH is already <= 1.0
      rh_p(:)  = rh(:) * 100.0

      ! Special handling for various aerosols
      SELECT CASE ( aerotype )

         !----------------
         ! Dust
         !----------------
         CASE ( cp_du )

            !=========================================================
            ! Hanel '76 (reported in Dentener etal '96). If rh>50%,
            ! then dust has enough water to allow gamma=0.1.
            !-------------------------------------------------------
            ! N.B. this should be updated as more information becomes
            ! available. Note that if TMIs (such as Cu and Fe) are
            ! present in high enough concentrations, then large uptake
            ! may be seen (gamma~0.8). See Thornton et al, 2008, JGR.
            !=========================================================
            WHERE ( rh_p >= 50 )
               cal_ho2(:) = 0.1e0
            ELSEWHERE
               !update when data available
               cal_ho2(:) = 0.05e0
            ENDWHERE

         !----------------
         ! Sulfate
         !----------------
         CASE ( cp_su )
                                
            !===========================================================
            ! Applying the temperature dependence for NaCl, it seems to 
            ! fit the data very well, so use this.
            ! 
            ! Data from Cooper & Abbatt '96 (Strat conditions) and
            ! Thornton & Abbatt '05.
            !
            ! Now data from Taketani et al 2008.
            ! Measurements on (NH4)2SO4 at 296K, tentatively recommend 
            ! gamma = 0.15, though there is dependence with rh.
            !===========================================================

            ! rh_p factor doesn't work below 35%rh
            WHERE ( rh_p <= 35 )
               cal_ho2(:) = 0.01E0
            ELSEWHERE
              ! Calculate temperature factor. Fix rel to 0.11 
               factt(:) = (5.66E-5 * EXP( 1560 / temp(:) ) ) / 0.11
              ! Calc rh_p factor from Taketani et al 2008
               facth(:) = (-2.86 * EXP(-0.078*rh_p(:)) + 0.192 ) / 0.11
              ! Calc gamma
               cal_ho2(:) = 0.11 * factt(:) * facth(:)
            ENDWHERE

            !Test with recommendations by Thornton and Abbatt 2008
            ! cal_ho2(:) = 0.15

         !----------------
         ! Black Carbon
         !----------------
         CASE ( cp_bc )

            !===========================================================
            ! Saathoff et al 2001. gamma=<0.01 on soot.
            !----------------------------------------------------------
            ! Ivanov (ref as dust) has measurements on soot, but can't
            ! find available anywhere.
            !===========================================================
            cal_ho2(:) = 0.01e0

         !----------------
         ! Organic Carbon
         !----------------
         CASE ( cp_oc , cp_so )

            !===========================================================
            ! Ivanov et al. Conference poster. Range of organic surfaces
            ! gamma=1d-4 to 5d-2 @ room temp.
            !===========================================================
            ! mid-value of available data
            cal_ho2(:) = 2.5e-2

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( cp_cl )
          
            !===========================================================
            ! Taketani etal ?. For NaCl, gamma<0.01 @ 20%rh solid
            ! (at room temp)                   0.05 @ 45%rh aq
            !---------------------------------------------------------
            ! Remorov etal 2002. Study on solid NaCl
            ! gamma=(5.66 pm 3.62)d-5 exp((1560 pm 140)/temp)
            !===========================================================

            ! Temperature dependence (for solid NaCl)
            WHERE ( rh_p >= 62 )
               cal_ho2(:) = 0.05e0
            ELSEWHERE
               cal_ho2(:) = 5.66e-5 * EXP( 1560 / temp(:) )
            ENDWHERE

         !----------------         
         ! Default
         !----------------
         CASE DEFAULT
            cmessage='Unknown aerosol surface for HO2 self reaction'
            WRITE (6,*) cmessage,' aerotype=',aerotype
            WRITE (6,*) 'AEROSOL TYPE =',aerotype
            CALL EREPORT('CAL_HO2.UKCA_TROP_HETCHEM',aerotype,cmessage)

      END SELECT   
         
      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_HO2',zhook_out,  &
                               zhook_handle)
      RETURN
      END FUNCTION CAL_HO2

! =========================================================================
!                F I R S T   O R D E R   L O S S 
! =========================================================================

      FUNCTION CAL_LOSS1K( sarea, wetdp, aird, mode_gamma,              &
                           temp, molec_wt_s)

! =================================================================
!  Subroutine CAL_LOSS1K calculates the 1st-order loss rate of
!  species on wet aerosol surface. (Adapted from lwh, jyl, gmg, djj,
!  7/1/94; bmy, 4/1/03, 4/7/04)
! =================================================================
! The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
! is computed as:
!
!      K [1/s] = sarea / [ wetdp/dfkg + 4./(gamma * xmms) ]        
! 
! where XMMS = Mean molecular speed [cm/s] = sqrt(8R*T/pi*M) for Maxwell 
!              [Formulae gives [m/s] if R in J/(mol.K) and M in kg/mol]
!       DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! NOTES:
! (1 ) Return a default value if RADIUS is zero (i.e. is smaller
!       than a very small number) (bmy, 4/7/04)
!********************************************************************

      IMPLICIT NONE

! Arguments
      REAL, INTENT(IN) :: sarea(nbox)         ! surface area of wet
                                              !  aerosol/vol of air [cm2/cm3]
      REAL, INTENT(IN) :: wetdp(nbox)         ! radius of wet aerosol [cm],
                                              ! order of 0.01-10 um
      REAL, INTENT(IN) :: aird(nbox)          ! density of air [#/cm3]
      REAL, INTENT(IN) :: mode_gamma(nbox)    ! Uptake coefficient [unitless] ~ 0.1
      REAL, INTENT(IN) :: temp(nbox)          ! temperature [K]
      REAL, INTENT(IN) :: molec_wt_s          ! molecular weight of species [g/mole]

! Input from module
!      rmol                ! gas constant [J/(mol.K)]
!      avc                 ! avogadros number
!      m_air               ! Molec Wt. of air kg/mol

! Function return value
      REAL :: cal_loss1k(nbox)                ! Reaction coefficient on aerosol

! Local variables
      REAL     :: dfkg(nbox)                  ! Gas phase diffusion coefficient
      REAL     :: rho_a(nbox)                 ! Air density in [kg/m^3]
      REAL     :: mw_s                        ! molec_weight of species [kg/mol]
      REAL, PARAMETER :: cf=1e-2              ! conversion factor 1/(m/s) to 1/(cm/s)
      REAL, PARAMETER :: dq=4.5e-10           ! Diameter of air molecule [m]

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      !=============================================================
      ! CAL_LOSS1K begins here!
      !=============================================================

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_LOSS1K',zhook_in,  &
                               zhook_handle)

      mw_s = molec_wt_s/1e3
      rho_a(:) = aird(:)*m_air*1e6/avc
      IF (PrintStatus >= PrStatus_Diag)                                 &
        WRITE(6,*) 'RHO_A: ',maxval(rho_a),minval(rho_a)
      dfkg(:)=0.0  ! for debug only
      WHERE (sarea < 0.0 .OR. wetdp < 1e-30 .OR. mode_gamma < 1e-30)

         ! Use default value if AREA is negative
         ! or if RADIUS is either zero or negative
         cal_loss1k = 1.0e-3

      ELSEWHERE

! dfkg = Gas phase diffusion coeff [cm2/s] (order of 0.1)
! see  Bauer et al., 2004, JGR 109, D02304 (eqn 7), with factor of 1e4 => [cm^2/s]
         dfkg(:) = (3e4/(8.0 * avc * dq**2 * rho_a(:))) *               &
              SQRT((rmol*temp(:)*m_air/(2.0*pi))*(mw_s + m_air)/mw_s)
         cal_loss1k(:) = sarea(:)/( wetdp(:)/dfkg(:) +                  &
                 (cf/mode_gamma(:))*SQRT(2.0*pi*mw_s/(rmol*temp(:))))

      ENDWHERE
      IF (PrintStatus >= PrStatus_Diag)                                 &
        WRITE(6,*) 'dfkg: ',maxval(dfkg),minval(dfkg)

      IF (lhook) CALL dr_hook('UKCA_TROP_HETCHEM:CAL_LOSS1K',zhook_out, &
                              zhook_handle)

      RETURN
      END FUNCTION CAL_LOSS1K

      END SUBROUTINE UKCA_TROP_HETCHEM
