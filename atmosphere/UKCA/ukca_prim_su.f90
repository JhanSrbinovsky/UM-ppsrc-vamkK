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
!    Calculates primary particulate sulfate emissions.
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
      SUBROUTINE UKCA_PRIM_SU(row_length, rows, model_levels,           &
        so2_high_level, parfrac, emanso2, emvolconso2, emvolexpso2,     &
        embiomso2, iso2ems, aer_mas_primsu, aer_num_primsu)

!--------------------------------------------------------
!
!     Calculates primary particulate sulfate emissions
!
!     Assumes particulate emissions emitted as lognormal modes
!     with geometric number mean diameters of GMDIAM and geometric
!     standard deviations of GSD.
!
!     Inputs
!     -------
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     MD         : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!     PARFRAC    : Fraction of anthrop. S mass emitted as ptcls
!     EMANSO2    : Anthrop. SO2 ems rates (N_EM_TYPE) (kgSO2/box/s)
!     EMVOLCONSO2: Natural SO2 ems rates (from cont. volc. src)
!                                                       (in kgSO2/box/s)
!     EMVOLEXPSO2(NBOX) Natural SO2 ems rates (from expl. volc. src)
!                                                       (in kgSO2/box/s)
!     EMBIOMSO2(NBOX) Biomass burning SO2 ems rates (kgSO2/box/s)
!     ISO2EMS    : Switch for size distribution of primary sulfate ems
!
!     Outputs
!     -------
!     AER_MAS_PRIMSU:
!     AER_NUM_PRIMSU:
!
!     Local variables
!     ---------------
!     N_EM_TYPES          : Number of emission types for ISO2EMS setting
!     N_EM_MODES          : Number of emission modes for ISO2EMS setting
!     EM_TYPE_MODE_GMDIAM : Geom. mean diam. of ems type & mode (nm)
!     EM_TYPE_MODE_GSD    : Geom. st. dev. of ems type & mode
!     EM_TYPE_MODE_FRAC   : Fraction of ems type which goes to this mode
!     EM_TYPE_MODE_MODE   : Index of mode to emit ems type & mode into
!     N_EM_BIOMSO2_MODES  : Number of emission modes for biomass SO4
!     EM_BIOMSO2_MODE_MODE: Index of mode to emit biomass  SO4 into
!     EM_BIOMSO2_MODE_GMDIAM: Geom. mean diam. of each biomass SO4 mode
!     EM_BIOMSO2_MODE_GSD : Geom. std. dev. of each biomass SO4 mode
!     n_em_volso2_modes   : Number of emission modes for volcanic SO4
!     EM_VOLSO2_MODE_GMDIAM: Geom. mean diam. of each volcanic SO4 mode
!     EM_VOLSO2_MODE_GSD  : Geom. std. dev. of each volcanic SO4 mode
!     EM_VOLSO2_MODE_MODE : Index of mode to emit volcanic SO4 into
!     EM_TYPE_MASS        : Mass emitted in ems type (kgH2SO4/box/s)
!     MODE_EM        : Mode into which emit into
!     LGSD           : Natural log of geometric standard deviation
!     DELN           : New and change in number concentration in
!                       (equivalent kg/m2/s as dry air)
!     PARMASS        : Mass emitted to box (type) (kgH2SO4/box/s)
!     MODEMASS       : Mass emitted to box (type & mode) (kgH2SO4/box/s)
!     MODEVOL        : Vol. conc. emitted (type & mode) (nm3/box/s)
!     FACTOR         : Converts from number to kg units
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI        : 3.1415927
!     AVC        : Avogadro's constant (molecules per mole)
!     DN_EPS     : Value of DELN below which do not carry out process
!     EMS_EPS    : Ems flux value below which don't emit (kg/gridbox/s)
!     MM_DA      : Molar mass of dry air (kg/mol)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES     : Number of modes set
!     MM         : Molar masses of condensable components (kg/mol)
!     RHOCOMP    : Density of each of the aerosol components (kg/m3)
!                                                  or carry out process
!     CP_SU      : index of cpt which sulfate ptcl mass is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     MSOTWO     : Index of MOLWT, WTRATC and S0G for SO2
!     MM_GAS     : Array of molar masses for gas phase species (kg/mol)
!     Various indices for budget terms in AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,       ONLY: ppi, avc, dn_eps, ems_eps, mm_da
      USE UKCA_MODE_SETUP,      ONLY: nmodes, rhocomp, mm, cp_su
      USE UKCA_SETUP_INDICES,   ONLY: msotwo, mm_gas
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE ereport_mod,          ONLY: ereport
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: row_length                  ! Number of columns
      INTEGER, INTENT(IN) :: rows                        ! Number of rows
      INTEGER, INTENT(IN) :: model_levels                ! Number of model levs
      INTEGER, INTENT(IN) :: so2_high_level              ! Level for injection
                                                         !  of stack SO2 emiss
      INTEGER, INTENT(IN) :: iso2ems                     ! Switch to define size
                                                         !  distbn.
      REAL,    INTENT(IN) :: emanso2(row_length,rows,6)       
!  Anthrop. SO2 ems (N_EM_TYPE) (kgSO2/box/s) (1 = low, 2 = high)
      REAL, INTENT(IN)    :: emvolconso2(row_length,rows,model_levels)
!  Natural SO2 ems rates (from cont. volc. src) in (kgSO2/box/s)
      REAL, INTENT(IN)    :: emvolexpso2(row_length,rows,model_levels)
!  Natural SO2 ems rates (from expl. volc. src) in (kgSO2/box/s)
      REAL, INTENT(IN)    :: embiomso2(row_length,rows,model_levels)
!  Biomass burning SO2 ems rates (kgSO2/m^2/s)
      REAL, INTENT(INOUT) :: aer_mas_primsu(row_length,rows,        &
                                                model_levels,nmodes)
!  Aerosol mass change (molecules/cc/tstep)
      REAL, INTENT(INOUT) :: aer_num_primsu(row_length,rows,        &
                                                model_levels,nmodes)
!  Aerosol mass change (ptcls/cc/tstep)
      REAL, INTENT(IN)    :: parfrac
! Fraction of SO2 emissions to particulate

! Local variables
      INTEGER :: errcode                     ! Error code for ereport
      INTEGER :: imode                       ! Loop counter for modes
      INTEGER :: ia                          ! Loop counter for emission types
      INTEGER :: ib                          ! Loop counter for emission modes
      INTEGER :: ilev                        ! Injection level for emissions 
      INTEGER :: em_biomso2_mode_mode(2)     ! Index of mode to emit biomass
                                             !  SO4 into
      INTEGER :: em_volso2_mode_mode(2)      ! Index of mode to emit volcanic
                                             !  SO4 into
      INTEGER :: em_type_mode_mode(3,2)      ! Index of mode to emit ems type
                                             !  and mode into
      INTEGER :: n_em_types                  ! No of emission types
      INTEGER :: n_em_modes                  ! No of emission modes
      INTEGER :: n_em_biomso2_modes          ! No of emission modes for
                                             !  biomass SO4
      INTEGER :: n_em_volso2_modes           ! No of emission modes for
                                             !  volcanic SO4
      INTEGER :: mode_em                     ! Mode into which emit into
      REAL    :: parmass(row_length,rows)    ! Mass emitted to box (type)
                                             !  (kgH2SO4/box/s)
      REAL    :: modemass(row_length,rows)   ! Mass emitted to box (type
                                             !  and mode) (kgH2SO4/box/s)
      REAL    :: modevol(row_length,rows)    ! Vol. conc. emitted (type
                                             !  and mode) (nm3/box/s)
      REAL    :: deln(row_length,rows)       ! Change in number concentration
                                             ! in equivalent kg/m2/s as dry air
      REAL    :: parmass3d(row_length,rows,model_levels) 
                                             ! Mass emitted to box (type)
                                             !  (kgH2SO4/box/s)
      REAL    :: modemass3d(row_length,rows,model_levels)   
                                             ! Mass emitted to box (type
                                             !  and mode) (kgH2SO4/box/s)
      REAL    :: modevol3d(row_length,rows,model_levels)     
                                             ! Vol. conc. emitted (type
                                             !  and mode) (nm3/box/s)
      REAL    :: deln3d(row_length,rows,model_levels)   ! Change in no. conc
                                             ! in equivalent kg/m2/s as dry air
      REAL    :: lgsd                        ! Natural log of geometric
                                             !  standard deviation
      REAL    :: factor                      ! conversion factor (kg/molecule)
      REAL    :: em_type_mode_gmdiam(3,2)    ! Geom. mean diam. of ems type
                                             !  and mode (nm)
      REAL    :: em_type_mode_gsd(3,2)       ! Geom. st. dev. of ems type
                                             !  and mode
      REAL    :: em_type_mode_frac(3,2)      ! Fraction of ems type which goes
                                             !  to this mode
      REAL    :: em_type_mass(row_length,rows,3) ! Mass emitted in ems type
                                             !     (kgH2SO4/box/s)
      REAL    :: em_biomso2_mode_gmdiam(2)   ! Geom. mean diam. of each biomass
                                             !  SO4 mode
      REAL    :: em_biomso2_mode_gsd(2)      ! Geom. std. dev. of each biomass
                                             !  SO4 mode
      REAL    :: em_biomso2_mode_frac(2)     ! Fraction of ems type which goes
                                             !  to this mode
      REAL    :: em_volso2_mode_gmdiam(2)    ! Geom. mean diam. of each
                                             !  volcanic SO4 mode
      REAL    :: em_volso2_mode_gsd(2)       ! Geom. std. dev. of each
                                             !  volcanic SO4 mode
      REAL    :: em_volso2_mode_frac(2)      ! Fraction of ems type which goes
                                             !  to this mode

      CHARACTER(LEN=72) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_PRIM_SU',zhook_in,zhook_handle)

      IF (msotwo <= 0) THEN
        cmessage = ' Error in UKCA_PRIM_SU, MSOTWO<=0'
        errcode = 1
        WRITE(6,'(A72,A7,I10)') cmessage,' MSOTWO=',msotwo
        CALL EREPORT('UKCA_PRIM_SU',errcode,cmessage)
      END IF

      factor = mm_da/avc
! Converts from  molecules or particle per m2/s to kg/m2/s 

      IF (iso2ems == 1) THEN
! Follow BS95 as in GLOMAP version from Spracklen et al (2005)
!
! .. type 1 is GEIA SO2 ems from surface  src, size:BS95, alt: 0-100m
! .. type 2 is GEIA SO2 ems from elevated src, size:BS95, alt: so2_high_level
! ..  (from umui)

        n_em_types = 2
        n_em_modes = 2
        em_type_mode_gmdiam(1,1) = 10.0
        em_type_mode_gmdiam(1,2) = 70.0
        em_type_mode_gmdiam(2,1) = 10.0
        em_type_mode_gmdiam(2,2) = 70.0
        em_type_mode_gsd(1,1) = 1.6
        em_type_mode_gsd(1,2) = 2.0
        em_type_mode_gsd(2,1) = 1.6
        em_type_mode_gsd(2,2) = 2.0
        em_type_mode_frac(1,1) = 0.15
        em_type_mode_frac(1,2) = 0.85
        em_type_mode_frac(2,1) = 0.15
        em_type_mode_frac(2,2) = 0.85
        em_type_mode_mode(1,1) = 2
        em_type_mode_mode(1,2) = 2
        em_type_mode_mode(2,1) = 2
        em_type_mode_mode(2,2) = 2
        em_type_mass(:,:,1) = emanso2(:,:,1)*parfrac*mm(cp_su)/   &
                              mm_gas(msotwo)
        em_type_mass(:,:,2) = emanso2(:,:,2)*parfrac*mm(cp_su)/   &
                              mm_gas(msotwo)
! .. biomass burning not done for ISO2EMS=1 (EMBIOMSO2=0)
        n_em_volso2_modes = 2
        em_volso2_mode_gmdiam(1) = 30.0
        em_volso2_mode_gmdiam(2) = 80.0
        em_volso2_mode_gsd(1) = 1.8
        em_volso2_mode_gsd(2) = 1.8
        em_volso2_mode_frac(1) = 0.5
        em_volso2_mode_frac(2) = 0.5
        em_volso2_mode_mode(1) = 2
        em_volso2_mode_mode(2) = 2
      ELSE IF (iso2ems == 2) THEN
! Follow AEROCOM recommendations for ACB
!
! .. type 1: roads, off-road & domestic, size:TRAFFIC, alt: lowest layer
! .. type 2: industrial, power-plants, size:INDUSTR, alt:100-300m
! .. type 3: shipping, size:INDUSTR, alt:lowest layer

        n_em_types = 3
        n_em_modes = 1
        em_type_mode_gmdiam(1,1) = 30.0
        em_type_mode_gmdiam(2,1) = 1000.0
        em_type_mode_gmdiam(3,1) = 1000.0
        em_type_mode_gsd(1,1) = 1.8
        em_type_mode_gsd(2,1) = 2.0
        em_type_mode_gsd(3,1) = 2.0
        em_type_mode_frac(1,1) = 1.0
        em_type_mode_frac(2,1) = 1.0
        em_type_mode_frac(3,1) = 1.0
        em_type_mode_mode(1,1) = 2
        em_type_mode_mode(2,1) = 3
        em_type_mode_mode(3,1) = 3
        em_type_mass(:,:,1) = (emanso2(:,:,3)+emanso2(:,:,5)+           &
                               emanso2(:,:,6))*parfrac*mm(cp_su)/       &
                               mm_gas(msotwo)
        em_type_mass(:,:,2) = (emanso2(:,:,1)+emanso2(:,:,2))*          &
                              parfrac*mm(cp_su)/mm_gas(msotwo)
        em_type_mass(:,:,3) = emanso2(:,:,4)*parfrac*mm(cp_su)/         &
                              mm_gas(msotwo)
        n_em_biomso2_modes = 1
        em_biomso2_mode_gmdiam(1) = 80.0
        em_biomso2_mode_gsd(1) = 1.8
        em_biomso2_mode_frac(1) = 1.0
        em_biomso2_mode_mode(1) = 2
        n_em_volso2_modes = 2
        em_volso2_mode_gmdiam(1) = 30.0
        em_volso2_mode_gmdiam(2) = 80.0
        em_volso2_mode_gsd(1) = 1.8
        em_volso2_mode_gsd(2) = 1.8
        em_volso2_mode_frac(1) = 0.5
        em_volso2_mode_frac(2) = 0.5
        em_volso2_mode_mode(1) = 2
        em_volso2_mode_mode(2) = 2
      ELSE IF (iso2ems == 3 .OR. iso2ems == 5 .OR. iso2ems == 6 ) THEN
! Follow Stier05 modifications to AEROCOM ACB
!
! .. type 1 is roads, off-road & domestic, size: modifiedTRAFFIC,
!                                          alt: lowest layer
! .. type 2 is industrial, power-plants  , size: modifiedINDUSTR,
!                                          alt: 100-300m, so2_high_level here
! .. type 3 is shipping, size:modifiedINDUSTR, alt: lowest layer
! For the UM, high level SO2 emissions go into type 2, low level into type 3.

        n_em_types = 3
        n_em_modes = 2
        em_type_mode_gmdiam(1,1) = 60.0
        em_type_mode_gmdiam(1,2) = 150.0
        em_type_mode_gmdiam(2,1) = 150.0
        em_type_mode_gmdiam(2,2) = 1500.0
        em_type_mode_gmdiam(3,1) = 150.0
        em_type_mode_gmdiam(3,2) = 1500.0
        em_type_mode_gsd(1,1) = 1.59
        em_type_mode_gsd(1,2) = 1.59
        em_type_mode_gsd(2,1) = 1.59
        em_type_mode_gsd(2,2) = 2.0
        em_type_mode_gsd(3,1) = 1.59
        em_type_mode_gsd(3,2) = 2.0
        em_type_mode_frac(1,1) = 0.5
        em_type_mode_frac(1,2) = 0.5
        em_type_mode_frac(2,1) = 0.5
        em_type_mode_frac(2,2) = 0.5
        em_type_mode_frac(3,1) = 0.5
        em_type_mode_frac(3,2) = 0.5
        em_type_mode_mode(1,1) = 2
        em_type_mode_mode(1,2) = 3
        em_type_mode_mode(2,1) = 3
        em_type_mode_mode(2,2) = 4
        em_type_mode_mode(3,1) = 3
        em_type_mode_mode(3,2) = 4
        em_type_mass(:,:,1) = (emanso2(:,:,3)+emanso2(:,:,5)+           &
                               emanso2(:,:,6))* parfrac*mm(cp_su)/      &
                               mm_gas(msotwo)
        em_type_mass(:,:,2) = (emanso2(:,:,1)+emanso2(:,:,2))*          &
                               parfrac*mm(cp_su)/mm_gas(msotwo)
        em_type_mass(:,:,3) = emanso2(:,:,4)*parfrac*mm(cp_su)/         &
                              mm_gas(msotwo)
        n_em_biomso2_modes = 2
        em_biomso2_mode_gmdiam(1) = 60.0
        em_biomso2_mode_gmdiam(2) = 150.0
        em_biomso2_mode_gsd(1) = 1.59
        em_biomso2_mode_gsd(2) = 1.59
        em_biomso2_mode_frac(1) = 0.5
        em_biomso2_mode_frac(2) = 0.5
        em_biomso2_mode_mode(1) = 2
        em_biomso2_mode_mode(2) = 3
        n_em_volso2_modes = 2
        em_volso2_mode_gmdiam(1) = 60.0
        em_volso2_mode_gmdiam(2) = 150.0
        em_volso2_mode_gsd(1) = 1.59
        em_volso2_mode_gsd(2) = 1.59
        em_volso2_mode_frac(1) = 0.5
        em_volso2_mode_frac(2) = 0.5
        em_volso2_mode_mode(1) = 2
        em_volso2_mode_mode(2) = 3
      ELSE IF(iso2ems == 4) THEN
! GEIA Emissions with AEROCOM recommendations for ACB
!
! .. type 1 is low level emissions, size:TRAFFIC, alt: lowest layer
! .. type 2 is high level emissions, size:INDUSTR, alt: 100-300m

        n_em_types = 2
        n_em_modes = 2
        em_type_mode_gmdiam(1,1) = 60.0
        em_type_mode_gmdiam(1,2) = 150.0
        em_type_mode_gmdiam(2,1) = 150.0
        em_type_mode_gmdiam(2,2) = 1500.0
        em_type_mode_gsd(1,1) = 1.59
        em_type_mode_gsd(1,2) = 1.59
        em_type_mode_gsd(2,1) = 1.59
        em_type_mode_gsd(2,2) = 2.0
        em_type_mode_frac(1,1) = 0.5
        em_type_mode_frac(1,2) = 0.5
        em_type_mode_frac(2,1) = 0.5
        em_type_mode_frac(2,2) = 0.5
        em_type_mode_mode(1,1) = 2
        em_type_mode_mode(1,2) = 3
        em_type_mode_mode(2,1) = 3
        em_type_mode_mode(2,2) = 4
        em_type_mass(:,:,1) = emanso2(:,:,1)*parfrac*mm(cp_su)/         &
                              mm_gas(msotwo)
        em_type_mass(:,:,2) = emanso2(:,:,2)*parfrac*mm(cp_su)/         &
                              mm_gas(msotwo)
        n_em_biomso2_modes = 1
        em_biomso2_mode_gmdiam(1) = 80.0
        em_biomso2_mode_gsd(1) = 1.8
        em_biomso2_mode_frac(1) = 1.0
        em_biomso2_mode_mode(1) = 2
        n_em_volso2_modes = 2
        em_volso2_mode_gmdiam(1) = 30.0
        em_volso2_mode_gmdiam(2) = 80.0
        em_volso2_mode_gsd(1) = 1.8
        em_volso2_mode_gsd(2) = 1.8
        em_volso2_mode_frac(1) = 0.5
        em_volso2_mode_frac(2) = 0.5
        em_volso2_mode_mode(1) = 2
        em_volso2_mode_mode(2) = 2
      ELSE IF (iso2ems  == 7) THEN
! Follow Stier05 modifications to AEROCOM ACB

! .. type 1 includes all types for AEROCOM hindcast

        n_em_types = 2 ! anthro_noship and anthro_ship
        n_em_modes = 2
        em_type_mode_gmdiam(1,1) = 150.0
        em_type_mode_gmdiam(1,2) = 1500.0
        em_type_mode_gmdiam(2,1) = 150.0
        em_type_mode_gmdiam(2,2) = 1500.0
        em_type_mode_gsd(1,1) = 1.59
        em_type_mode_gsd(1,2) = 2.0
        em_type_mode_gsd(2,1) = 1.59
        em_type_mode_gsd(2,2) = 2.0
        em_type_mode_frac(1,1) = 0.5
        em_type_mode_frac(1,2) = 0.5
        em_type_mode_frac(2,1) = 0.5
        em_type_mode_frac(2,2) = 0.5
        em_type_mode_mode(1,1) = 3
        em_type_mode_mode(1,2) = 4
        em_type_mode_mode(2,1) = 3
        em_type_mode_mode(2,2) = 4
        em_type_mass(:,:,1) = (emanso2(:,:,1))*parfrac*mm(cp_su)/       &
                               mm_gas(msotwo)
        em_type_mass(:,:,2) = (emanso2(:,:,2))*parfrac*mm(cp_su)/       &
                               mm_gas(msotwo)
        n_em_biomso2_modes = 2
        em_biomso2_mode_gmdiam(1) = 60.0
        em_biomso2_mode_gmdiam(2) = 150.0
        em_biomso2_mode_gsd(1) = 1.59
        em_biomso2_mode_gsd(2) = 1.59
        em_biomso2_mode_frac(1) = 0.5
        em_biomso2_mode_frac(2) = 0.5
        em_biomso2_mode_mode(1) = 2
        em_biomso2_mode_mode(2) = 3
        n_em_volso2_modes = 2
        em_volso2_mode_gmdiam(1) = 60.0
        em_volso2_mode_gmdiam(2) = 150.0
        em_volso2_mode_gsd(1) = 1.59
        em_volso2_mode_gsd(2) = 1.59
        em_volso2_mode_frac(1) = 0.5
        em_volso2_mode_frac(2) = 0.5
        em_volso2_mode_mode(1) = 2
        em_volso2_mode_mode(2) = 3
      ELSE
        cmessage = ' iso2ems is not in range 1-7 '
        errcode = 1
        WRITE(6,'(A72,I10)') cmessage,iso2ems
        CALL EREPORT('UKCA_PRIM_SU',errcode,cmessage)
      END IF

!----------------------------------------------------------------------
! ..  This section does primary SO4 from SO2 from anthropogenic src

      DO ia = 1,n_em_types                      ! loop through emission types
        parmass(:,:) = em_type_mass(:,:,ia)     ! kgH2SO4/m2/s

! .. Set injection level ...
        IF (ia == 3) THEN
          ilev = so2_high_level                 ! Stack emissions
        ELSE
          ilev = 1
        END IF

        DO ib=1,n_em_modes                      ! loop through emission modes

! .. Calulate natural logs of standard deviation
          lgsd = LOG(em_type_mode_gsd(ia,ib))
! .. Store which mode to emit primary SO4 into
          mode_em = em_type_mode_mode(ia,ib)

! .. Particulate mass emissions (kg_H2SO4 per m2 per s) in mode
          modemass(:,:) = parmass(:,:)*em_type_mode_frac(ia,ib)
! .. Calculate total particle volume (nm3 per m2 per s)
          modevol(:,:) = 1E27*modemass(:,:)/rhocomp(cp_su)
! .. Calculate total particle number (particles per m2 per s) * factor
! ..  factor converts from ptcls/m2/s to equivalent kg/m2/s as dry air
          deln(:,:) = modevol(:,:)*factor/((ppi/6.0)*                   &
            (em_type_mode_gmdiam(ia,ib)**3)*EXP(4.5*lgsd*lgsd))

! .. Sum up each emitted mass into mode via new budget term
          aer_mas_primsu(:,:,ilev,mode_em) =                            &
            aer_mas_primsu(:,:,ilev,mode_em) +                          &
            modemass(:,:)*factor*avc/mm(cp_su)               ! kg/m2/s

! sum up each emitted number into mode via new budget term
          aer_num_primsu(:,:,ilev,mode_em) =                            &
            aer_num_primsu(:,:,ilev,mode_em) + deln(:,:)     ! (equiv)-kg/m2/s

        END DO   ! ib = 1,n_em_modes
      END DO   ! ia = 1,n_em_types

!----------------------------------------------------------------------
! ..  This section does primary SO4 from SO2 from biomass burning src

      parmass3d(:,:,:) = embiomso2(:,:,:)*parfrac*                      &
                         mm(cp_su)/mm_gas(msotwo)            ! in kgH2SO4/m2/s

      DO ib=1,n_em_biomso2_modes          ! loop over No of emission modes

! .. Store which mode to emit primary SO4 into
        mode_em = em_biomso2_mode_mode(ib)

! .. Calculate natural logs of standard deviation
        lgsd = LOG(em_biomso2_mode_gsd(ib))

! .. Particulate mass emissions (kg_H2SO4 per m2 per s) in mode
        modemass3d = parmass3d*em_biomso2_mode_frac(ib)
! .. Calculate total particle volume (nm3 per m2 per s)
        modevol3d = 1e27*modemass3d/rhocomp(cp_su)
! .. Calculate total particle number (per m2 per s) * factor
!     factor converts from number/m2/s to equiv-mass/m2/s
        deln3d = modevol3d*factor/((ppi/6.0)*                           &
           (em_biomso2_mode_gmdiam(ib)**3)*EXP(4.5*lgsd*lgsd))


! sum up each emitted mass   into mode via new budget term
        aer_mas_primsu(:,:,:,mode_em) =                                 &
          aer_mas_primsu(:,:,:,mode_em) +                               &
          modemass3d(:,:,:)*factor*avc/mm(cp_su)             ! kg/m2/s

! sum up each emitted number into mode via new budget term
        aer_num_primsu(:,:,:,mode_em) =                                 &
          aer_num_primsu(:,:,:,mode_em) + deln3d(:,:,:)      ! equiv-mass/m2/s

      END DO    ! loop over ib

!----------------------------------------------------------------------
! ..  This section does primary SO4 from SO2 from volcanic (cont.) src

      parmass3d(:,:,:) = emvolconso2(:,:,:)*parfrac*                    &
                         mm(cp_su)/mm_gas(msotwo)            ! in kgH2SO4/m2/s

      DO ib=1,n_em_volso2_modes

! .. Calculate natural logs of standard deviation
        lgsd = LOG(em_volso2_mode_gsd(ib))

! .. Store which mode to emit primary SO4 into
        mode_em = em_volso2_mode_mode(ib)

! .. Particulate mass emissions (kg_H2SO4 per m2 per s) in mode
        modemass3d = parmass3d*em_volso2_mode_frac(ib)
! .. Calculate total particle volume (nm3 per m2 per s)
        modevol3d = 1E27*modemass3d/rhocomp(cp_su)
! .. Calculate total particle number (per m2 per s) * factor
! ..  factor converts from number/m2/s to equiv-mass/m2/s
        deln3d = modevol3d*factor/((ppi/6.0)*                           &
                   (em_volso2_mode_gmdiam(ib)**3)*EXP(4.5*lgsd*lgsd))

! sum up each emitted mass   into mode via new budget term
        aer_mas_primsu(:,:,:,mode_em) =                                 &
          aer_mas_primsu(:,:,:,mode_em) +                               &
          modemass3d(:,:,:)*factor*avc/mm(cp_su)            ! kg/m2/s

! sum up each emitted number into mode via new budget term
        aer_num_primsu(:,:,:,mode_em) =                                 &
          aer_num_primsu(:,:,:,mode_em) + deln3d(:,:,:)     ! equiv-kg/m2/s

      END DO      ! ib = 1,n_em_volso2_modes

!----------------------------------------------------------------------
! ..  This section does primary SO4 from SO2 from volcanic (expl.) src

      parmass3d(:,:,:) = emvolexpso2(:,:,:)*parfrac*                    &
                       mm(cp_su)/mm_gas(msotwo)             ! in kgH2SO4/m2/s

      DO ib = 1,n_em_volso2_modes

! .. Calculate natural logs of standard deviation
        lgsd = LOG(em_volso2_mode_gsd(ib))

! .. Store which mode to emit primary SO4 into
        mode_em = em_volso2_mode_mode(ib)

! .. Particulate mass emissions (kg_H2SO4 per m2 per s) in mode
        modemass3d = parmass3d*em_volso2_mode_frac(ib)
! .. Calculate total particle volume (nm3 per m2 per s)
        modevol3d = 1e27*modemass3d/rhocomp(cp_su)
! .. Calculate total particle number (per m2 per s) * factor
! ..  factor converts from number/m2/s to equiv-mass/m2/s
        deln3d = modevol3d*factor/((ppi/6.0)*                           &
                   (em_volso2_mode_gmdiam(ib)**3)*EXP(4.5*lgsd*lgsd))


! .. Sum up each emitted mass   into mode via new budget term
        aer_mas_primsu(:,:,:,mode_em) =                                 &
          aer_mas_primsu(:,:,:,mode_em) +                               &
          modemass3d(:,:,:)*factor*avc/mm(cp_su)              ! kg/m2/s

! .. Sum up each emitted number into mode via new budget term
        aer_num_primsu(:,:,:,mode_em) =                                 &
          aer_num_primsu(:,:,:,mode_em) + deln3d(:,:,:)       ! equiv-kg/m2/s

      END DO        !  ib = 1,n_em_volso2_modes

      IF (lhook) CALL dr_hook('UKCA_PRIM_SU',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_PRIM_SU
