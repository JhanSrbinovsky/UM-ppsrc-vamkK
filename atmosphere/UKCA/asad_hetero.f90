! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Dummy heterogeneous chemistry routine.
!
!     The purpose of this routine is to set and return the heterogeneous
!     reaction rates. If the user has heterogeneous chemistry turned on
!     then this subroutine will be called. The user must supply their
!     own version of this routine to compute the heterogeneous rates.
!
!     Note that this subroutine is called repeatedly. It should not
!     therefore be used to do any I/O unless absolutely necessary. The
!     routine inihet is provided to initialise the heterogeneous chemist
!     by reading in files etc.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE ASAD_HETERO(n_points, cld_f, cld_l, rc_het)

      USE ASAD_MOD,        ONLY: t, p, tnd, rk, ih_o3, ih_h2o2, ih_so2, &
                                 ih_hno3, ihso3_h2o2, iho2_h, in2o5_h,  &
                                 iso3_o3, ihso3_o3, iso2_oh, ih2o2_oh,  &
                                 ihno3_oh, spb, sph, nbrkx, nhrkx,      &
                                 jpspb, jpsph, jpeq
      USE ukca_option_mod, ONLY: L_ukca_achem, L_ukca_trophet,          &
                                 L_ukca_tropisop, L_ukca_mode,          &
                                 jpbk, jphk, jpdw
      USE UKCA_CONSTANTS,  ONLY: avc, rhow, m_air, H_plus
      USE parkind1,        ONLY: jprb, jpim
      USE yomhook,         ONLY: lhook, dr_hook
      USE ereport_mod,     ONLY : ereport
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n_points         ! No of spatial points

      REAL, INTENT(IN) :: cld_f(n_points)     ! Cloud fraction
      REAL, INTENT(IN) :: cld_l(n_points)     ! Cloud liquid water (kg/kg)
      REAL, INTENT(IN) :: rc_het(n_points,2)  ! Heterog. Chem. Rates (tropospheric)

!     Local variables

      REAL, PARAMETER    :: qcl_min = 1.0E-12 ! do calcs when qcl > qcl_min
      REAL               :: vr(n_points)      ! volume ratio
      REAL               :: fdiss(n_points, jpdw, jpeq+1)
                                              ! fractional dissociation array
                                              ! final index: 1) dissolved
                                              !              2) 1st dissociation
                                              !              3) 2nd dissociation

      INTEGER            :: asad_findreaction ! integer function
      INTEGER            :: icode             ! Error code
      CHARACTER (LEN=72) :: cmessage          ! Error message
      CHARACTER(LEN=10)       :: prods(2)          ! Products
      LOGICAL, SAVE      :: first = .TRUE.    ! Identifies first call
      LOGICAL            :: todo(n_points)    ! T where cloud frac above threshold

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('ASAD_HETERO',zhook_in,zhook_handle)

!       1. Identify equations and calculate fractional dissociation
!          --------------------------------------------------------

      IF (first) THEN

        IF (L_ukca_achem) THEN
! Check that the indicies of the aqueous arrays are identified
          IF (ih_o3 == 0 .OR. ih_h2o2 == 0 .OR. ih_so2 == 0 .OR.        &
              ih_hno3 == 0 ) THEN
            cmessage=' Indicies for Aqueous chemistry uninitialised'//  &
                     ' - O3, H2O2, SO2, and HNO3 must be made '//       &
                     ' soluble species in CHCH_DEFS array'
            write(6,*) 'ih_o3:   ',ih_o3
            write(6,*) 'ih_h2o2: ',ih_h2o2
            write(6,*) 'ih_hno3: ',ih_hno3
            write(6,*) 'ih_so2:  ',ih_so2
            icode=1
            CALL EREPORT('ASAD_HETERO',icode,cmessage)
          END IF

          IF (L_ukca_trophet .AND. .NOT. L_ukca_mode) THEN
            cmessage=' Tropospheric heterogeneous chemistry is flagged'// &
                   ' but MODE aerosol scheme is not in use'
            icode=1
            CALL EREPORT('ASAD_HETERO',icode,cmessage)
          END IF

          ihso3_h2o2 = 0
          iso3_o3    = 0
          ihso3_o3   = 0
          ih2o2_oh   = 0
          ihno3_oh   = 0
          in2o5_h    = 0
          iho2_h     = 0

          prods = (/'NULL0     ','          '/)
! DEPENDS ON: asad_findreaction
          ihso3_h2o2 = asad_findreaction( 'SO2       ', 'H2O2      ',   &
                                 prods, 2, sph, nhrkx, jphk+1, jpsph )
          prods = (/'NULL1     ','          '/)   ! Identifies HSO3- + O3(aq)
! DEPENDS ON: asad_findreaction
          ihso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',     &
                                   prods, 2, sph, nhrkx, jphk+1, jpsph )
          prods = (/'NULL2     ','          '/)   ! Identifies SO3-- + O3(aq)
! DEPENDS ON: asad_findreaction
          iso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',      &
                                   prods, 2, sph, nhrkx, jphk+1, jpsph )

          prods = (/'H2O       ','HO2       '/)
! DEPENDS ON: asad_findreaction
          ih2o2_oh = asad_findreaction( 'H2O2      ', 'OH        ',     &
                                 prods, 2, spb, nbrkx, jpbk+1, jpspb )
          prods = (/'H2O       ','NO3       '/)
! DEPENDS ON: asad_findreaction
          ihno3_oh = asad_findreaction( 'HONO2     ', 'OH        ',     &
                                 prods, 2, spb, nbrkx, jpbk+1, jpspb )

          IF (ihso3_h2o2 == 0 .OR. iso3_o3 == 0 .OR. ih2o2_oh == 0 .OR. &
              ihso3_o3 == 0 .OR. ihno3_oh == 0 ) THEN
            write(6,*) 'ihso3_h2o2: ',ihso3_h2o2
            write(6,*) 'ihso3_o3: ',ihso3_o3
            write(6,*) 'iso3_o3: ',iso3_o3
            write(6,*) 'ih2o2_oh: ',ih2o2_oh
            write(6,*) 'ihno3_oh: ',ihno3_oh
            cmessage=' Heterogenous chemistry called, but eqns'//       &
                      ' not found - see output'
            icode = 1
            CALL EREPORT('ASAD_HETERO',icode,cmessage)
          END IF
        END IF   ! L_ukca_achem

! Search for tropospheric heterogeneous reactions
        prods = (/'HONO2     ','          '/)
! DEPENDS ON: asad_findreaction
        in2o5_h = asad_findreaction( 'N2O5      ', '          ',        &
                                 prods, 2, sph, nhrkx, jphk+1, jpsph )
        prods = (/'H2O2      ','          '/)
! DEPENDS ON: asad_findreaction
        iho2_h = asad_findreaction( 'HO2       ', '          ',         &
                                 prods, 2, sph, nhrkx, jphk+1, jpsph )

        IF (L_ukca_trophet .AND. (iho2_h == 0 .OR. in2o5_h == 0)) THEN
          write(6,*) 'in2o5_h: ',in2o5_h
          write(6,*) 'iho2_h: ',iho2_h
          cmessage=' Tropospheric heterogenous chemistry is flagged,'// &
                   ' but equations not found - see output'
          icode = 1
          CALL EREPORT('ASAD_HETERO',icode,cmessage)
        END IF   ! iho3_h=0 etc

        first = .FALSE.

      END IF      ! first

      IF (L_ukca_achem .AND. ANY(cld_l > qcl_min))&
        THEN
! DEPENDS ON: ukca_fdiss
        CALL UKCA_FDISS(n_points, qcl_min, t, p, cld_l, fdiss)
        todo(:) = cld_l(:) > qcl_min
      ELSE
        fdiss(:,:,:) = 0.0
        todo(:) = .FALSE.
      END IF

!       2. Calculate heterogeneous rates and reduce rates due to aqueous fraction
!          ----------------------------------------------------------------------

      IF (L_ukca_achem) THEN
        WHERE (todo(:))

! Convert clw in kg/kg to volume ratio
          vr(:) = cld_l(:)*tnd(:)*m_air*(1e6/avc)/rhow

! HSO3- + H2O2(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
        rk(:,ihso3_h2o2) = 2.1295E+14*EXP(-4430.0/t(:))*                &
           (H_plus/(1.0 + 13.0*H_plus))*cld_f(:)*fdiss(:,ih_so2,2)*     &
           fdiss(:,ih_h2o2,1)*1000.0/(avc*vr(:))

! HSO3- + O3(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
        rk(:,ihso3_o3) = 4.0113E+13*EXP(-5530.0/t(:))*                  &
                         cld_f(:)*fdiss(:,ih_so2,2)*fdiss(:,ih_o3,1)*   &
                         1000.0/(avc*vr(:))

! SO3-- + O3(aq) => SO4-- [Kreidenweis et al. (2003), optimised]
        rk(:,iso3_o3) = 7.43E+16*EXP(-5280.0/t(:))*cld_f(:)*            &
                           fdiss(:,ih_so2,3)*fdiss(:,ih_o3,1)*          &
                           1000.0/(avc*vr(:))

! H2O2 + OH ; HNO3 + OH : reduce to take account of dissolved fraction
        rk(:,ih2o2_oh) = rk(:,ih2o2_oh)*                                &
              (1.0 - (fdiss(:,ih_h2o2,1)+fdiss(:,ih_h2o2,2))*cld_f(:))
        rk(:,ihno3_oh) = rk(:,ihno3_oh)*                                &
              (1.0 - (fdiss(:,ih_hno3,1)+fdiss(:,ih_hno3,2))*cld_f(:))
        ELSEWHERE
          rk(:,ihso3_h2o2) = 0.0
          rk(:,ihso3_o3) = 0.0
          rk(:,iso3_o3) = 0.0
        ENDWHERE

      END IF      ! L_ukca_achem

      IF (L_ukca_trophet) THEN
! N2O5 => HNO3 (heterogenous)
        rk(:,in2o5_h) = rc_het(:,1)

! HO2 + HO2 => H2O2 (heterogenous)
        rk(:,iho2_h) = rc_het(:,2)
      ELSE
        IF (in2o5_h > 0) rk(:,in2o5_h) = 0.0
        IF (iho2_h > 0)  rk(:,iho2_h) = 0.0
      END IF

      IF (lhook) CALL dr_hook('ASAD_HETERO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ASAD_HETERO
