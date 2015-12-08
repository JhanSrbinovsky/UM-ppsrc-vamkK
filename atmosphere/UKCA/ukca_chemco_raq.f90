! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Calculate the reaction rate co-efficients for use in the
!  Backward Euler solver with RAQ chemistry (based on STOCHEM ).
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Calculates rate coefficients for RAQ chemistry
!   Units : cm, molecule, s
!   Inputs  : tc, m, h2o, o2
!   Outputs : rc
!   This routine was adapted from the modset stochem_chemistry.mf77
!   and some reaction rate coefficients were corrected following
!   Atkinson (2004, 2006), Atm. Chem. Phys
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_CHEMCO_RAQ (nr, n_pnts, tc, m, h2o, o2, rc)
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nr
      INTEGER, INTENT(IN) :: n_pnts

      REAL, INTENT(IN)    :: tc(n_pnts)      ! Temperature (K)
      REAL, INTENT(IN)    :: m(n_pnts)       ! Density (molecules cm-3)
      REAL, INTENT(IN)    :: h2o(n_pnts)     ! Water concentration (molec cm-3)
      REAL, INTENT(IN)    :: o2(n_pnts)      ! Oxygen concn (molec cm-3)
      REAL, INTENT(OUT)   :: rc(n_pnts,nr)   ! Rate constants

! Local variables

      INTEGER, PARAMETER :: npdep = 10       ! No. press-dep reactions
      REAL,    PARAMETER :: o2frac= 0.20948  ! o2 molar fraction in air
      REAL,    PARAMETER :: n2frac= 0.78084  ! n2 molar fraction in air
!
      INTEGER :: i, j, k                     ! Loop counts
      INTEGER :: r                           ! Reaction type
      INTEGER :: nklo, nkhi, nres ! Reaction nos of klo, khi and location
!                                   of result for pressure-dependent reacs
!
      REAL :: tfc = 1.0 / 300.0  ! Used in k's above
      REAL :: ar                 ! Arrhenius pre-exponential factor
      REAL :: n                  ! Temperature power (K)
      REAL :: ea                 ! Activation energy (as -Ea/R, units: K)
      REAL :: fac1               ! Temporary store
      REAL :: fac2               ! Temporary store
      REAL :: fac3               ! Temporary store
      REAL :: brn                ! Temporary store
      REAL :: f                  ! Appropriate value of Fc factor
      REAL :: fc                 ! Fc factor for R29
      REAL :: fac1v(n_pnts)      ! Temporary store vector
!      
      CHARACTER(LEN=72) :: cmessage  ! Error return message
!
! Set up rate coefficient data type
! Rate constant k = A.(T/300)**n.exp(E/T)
! E is -Ea/R, where Ea is the activation energy and R the ideal gas constant
! Hence Ea has units of K.
! r is the reaction type:
!  0 - no reaction
!  1 - k = A  (temperature-independent)
!  2 - k = A.(T/300)**n
!  3 - k = A.exp(E/T)
!  4 - k = A.(T/300)**n.exp(E/T)
!
      TYPE RR
        REAL    :: a    ! Pre-exponential factor
        REAL    :: n    ! Power for temperature
        REAL    :: e    ! Activation Energy (actually -Ea/R)
        INTEGER :: r    ! Reaction type (see above)
      ENDTYPE RR
!
! Set up type to hold information needed to calculate pressure-dependent
! rate constants
      TYPE PD
        INTEGER :: klo  ! Reaction nr holding low-pressure rate constant
        INTEGER :: khi  ! Reaction nr holding high-pressure rate constant
        INTEGER :: res  ! Rate constant to store final result in
        REAL    :: fc   ! Fc factor
      ENDTYPE PD
!
! Holds A, n, Ea for each rate constant
      TYPE(RR), DIMENSION(:), ALLOCATABLE, SAVE :: rk
      TYPE(RR) :: zero_rate = RR(0.0, 0.0, 0.0, 0)
!
! Pressure-dependent reactions
      TYPE(PD), DIMENSION(npdep) :: pdep = (/                           &
!          klo khi res  Fc
        PD(  4,  5,  5, 0.85),                                          &
! O + NO + M = NO2 + M
        PD( 20, 46, 20, 0.35),                                          &
! NO2 + NO3 + M = N2O5 + M
        PD( 21, 47, 21, 0.4),                                           &
! NO2 + OH + M = HNO3 + M
        PD( 22, 48, 22, 0.6),                                           &
! NO2 + HO2 + M = HO2NO2 + M
        PD( 23, 25, 23, 0.6),                                           &
! HO2NO2 + M = HO2 + NO2 + M
        PD( 29, 49, 29, 0.35),                                          &
! N2O5 + M = NO2 + NO3 + M
        PD( 76, 77, 77, 0.3),                                           &
! CH3COO2 + NO2 + M = CH3COO2NO2 (PAN) + M
        PD( 78, 82, 78, 0.3),                                           &
! CH3COO2NO2 + M = CH3COO2 + NO2 + M
        PD(108,109,109, 0.48),                                          &
! OH + C2H4 + M (+O2) = HOC2H4O2 + M
        PD(125,130,125, 0.5)                                            &
! OH + C3H6 + M (+O2) = HOC3H6O2 + M
       /)
!
      LOGICAL, SAVE :: first = .TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_CHEMCO_RAQ',zhook_in,zhook_handle)

!
      IF (first) THEN     ! Set up array rk on first call

        ALLOCATE(rk(nr))  ! Holds A, n, Ea for each rate constant
!
! Set up rk array 19 lines at a time. 19 continuation lines are Fortran-90
! standard limit (according to warning messages when more are specified)
!
!              A       n     -Ea/R  rtype   Num   Reaction Details
      rk(1:19) = (/                                                     &
        RR( 6.0e-34, -2.6,     0.0,  2),                                &
       !   1 : O + O2 + M = O3 + M 
        zero_rate,                                                      & 
        zero_rate,                                                      &
        RR( 1.0e-31, -1.6,     0.0,  2),                                &
       !   4 : O + NO + M = NO2 + M  klo
        RR( 3.0e-11,  0.3,     0.0,  2),                                &
       !   5 : O + NO + M = NO2 + M  khi 
        RR( 3.2e-11,  0.0,    70.0,  3),                                &
       !   6 : O(1D) + M = O(3P) + M
        RR( 1.8e-11,  0.0,   110.0,  3),                                &
       !   7 : O(1D) + M 2nd expression
        RR( 2.2e-10,  0.0,     0.0,  1),                                &
       !   8 : O(1D) + H2O = 2 OH
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR( 1.4e-12,  0.0, -1310.0,  3),                                &
       !  11 : NO + O3 = NO2 + O2
        RR( 1.4e-13,  0.0, -2470.0,  3),                                &
       !  12 : NO2 + O3 = NO3 + O2
        RR( 1.7e-12,  0.0,  -940.0,  3),                                &
       !  13 : OH + O3 = HO2 + O2
        RR( 2.03e-16, 4.57, 693.0,   4),                                &
       !  14 : HO2 + O3 = OH + O2 + O2
        RR( 1.8e-11,  0.0,   110.0,  3),                                &
       !  15 : NO + NO3 = NO2 + NO2
        RR( 5.5e-12,  0.0,   188.0,  3),                                &
       !  16 : NO2 + O = NO + O2
        RR( 3.6e-12,  0.0,   270.0,  3),                                &
       !  17 : NO + HO2 = OH + NO2
        zero_rate,                                                      &
        RR( 4.5e-14,  0.0, -1260.0,  3)                                 &
        /)  !  19 : NO2 + NO3 = NO + NO2 + O2

      rk(20:38) = (/                                                    &
        RR( 3.6e-30, -4.1,     0.0,  2),                                &
       !  20 : NO2 + NO3 + M = N2O5 + M  klo see R46
        RR( 3.3e-30, -3.0,     0.0,  2),                                &
       !  21 : NO2 + OH + M = HNO3 + M  klo see R47
        RR( 1.8e-31, -3.2,     0.0,  2),                                &
       !  22 : NO2 + HO2 + M = HO2NO2 + M  klo see R48
        RR( 4.1E-5,  0.0, -10650.0,  3),                                &
       !  23 : HO2NO2 + M = HO2 + NO2 + M   klo see R25
        RR( 3.2e-13, 0.0,    690.0,  3),                                &
       !  24 : OH + HO2NO2 = H2O + NO2 + O2
        RR( 4.8E15,  0.0, -11170.0,  3),                                &
       !  25 : HO2NO2 + M = HO2 + NO2 + M   khi see R23
        zero_rate,                                                      &
        RR( 8.5e-13,  0.0, -2450.0,  3),                                &
       !  27 : NO3 + NO3 = NO2 + NO2 + O2
        RR( 3.5e-12,  0.0,  -925.0,  3),                                &
       !  28 : OH + NH3  = NH2 + H2O (reaction not used)
        RR( 1.3e-03,-3.5,-11000.0,  4),                                 &
       !  29 : N2O5 + M = NO2 + NO3 + M  klo
        RR( 4.8e-11,  0.0,   250.0,  3),                                &
       !  30 : HO2 + OH = H2O + O2
        RR( 2.9e-12,  0.0,  -160.0,  3),                                &
       !  31 : OH + H2O2 = H2O + HO2
        RR( 4.2e-12,  0.0,     0.0,  1),                                &
       !  32 : NO3 + HO2 = HNO3 + O2
        RR( 5.5e-12,  0.0, -2000.0,  3),                                &
       !  33 : OH + H2 (+O2) = HO2 + H2O
        RR( 3.5e-12,  0.0,     0.0,  1),                                &
       !  34 : NO3 + HO2 = OH + NO2 + O2
        RR( 2.4e-14,  0.0,   460.0,  3),                                &
       !  35 : OH + HNO3 = NO3 + H2O  See R50,R51
        RR( 2.2e-13,  0.0,   600.0,  3),                                &
       !  36 : HO2 + HO2 = H2O2 + O2
        RR( 1.9e-33,  0.0,   980.0,  3),                                &
       !  37 : HO2 + HO2 (+M) = H2O2 + O2
        RR( 1.4e-21,  0.0,  2200.0,  3)                                 &
       !  38 : HO2 + HO2 (+H2O) = H2O2 + O2
        /)

      rk(39:57) = (/                                                    &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR( 2.66e-12, 0.0,   200.0,  3),                                &
        !  44 : OH + CH3OOH = CH3O2 + H2O
        RR( 1.14e-12, 0.0,   200.0,  3),                                &
        !  45 : OH + CH3OOH = OH + HCHO + H2O
        RR( 1.9e-12,  0.2,     0.0,  2),                                &
        !  46 : NO2 + NO3 + M = N2O5 + M  khi for R20
        RR( 4.1e-11,  0.0,     0.0,  1),                                &
        !  47 : NO2 + OH + M = HNO3 + M  khi see R21
        RR( 4.7e-12,  0.0,     0.0,  1),                                &
        !  48 : NO2 + HO2 + M = HO2NO2 + M  khi see R22
        RR( 9.7e+14,  0.1,-11080.0,  4),                                &
        !  49 : N2O5 + M = NO2 + NO3 + M  khi see R29
        RR( 2.7e-17,  0.0,  2199.0,  3),                                &
        !  50 : OH + HNO3 = NO3 + H2O 2nd term for R35
        RR( 6.5e-34,  0.0,  1335.0,  3),                                &
        !  51 : OH + HNO3 = NO3 + H2O 3rd term for R35
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate /)                                                                      

      rk(58:76) = (/                                                    &
        zero_rate,                                                      &
        RR(1.85E-12,  0.0, -1690.0,  3),                                &
        !  59 : OH + CH4 = CH3O2 + H2O
        RR(2.3e-12,   0.0,   360.0,  3),                                &
        !  60 : NO + CH3O2 +O2 = HCHO + HO2 + NO2
        RR(7.4e-13,   0.0,  -520.0,  3),                                &
        !  61 : CH3O2 + CH3O2 = 2 HCHO + 2 HO2
        RR(1.03e-13,  0.0,   365.0,  3),                                &
        !  62 : Total rate const: CH3O2 + CH3O2 = Products
        RR(7.3e-12,   0.0,  -620.0,  3),                                &
        !  63 : CH3OH + OH = HO2 + HCHO + H2O
        zero_rate,                                                      &
        RR( 4.1e-13,  0.0,   750.0,  3),                                &
        !  65 : CH3O2 + HO2 = CH3OOH + O2
        RR( 5.4e-12,  0.0,   135.0,  3),                                &
        !  66 : OH + HCHO = HO2 + CO + H2O
        RR( 5.8e-16,  0.0,     0.0,  1),                                &
        !  67 : NO3 + HCHO = HO2 + CO + HNO3
        zero_rate,                                                      &
        RR(3.54e-33,  0.0,     0.0,  1),                                &
        !  69 : OH + CO Press dependent term
        RR( 1.5e-13,  0.0,     0.0,  1),                                &
        !  70 : OH + CO = HO2 + CO2
        RR(6.9e-12,   0.0, -1000.0,  3),                                &
        !  71 : OH + C2H6 = C2H5O2
        RR( 2.6e-12,  0.0,     0.0,  1),                                &
        !  72 : C2H5O2 + NO = CH3CHO + HO2 + NO2
        RR( 2.0e-13,  0.0,     0.0,  1),                                &
        !  73 : C2H5O2 + CH3O2 = CH3CHO + 2 HO2 + HCHO + O2
        RR( 4.4e5,    0.0, -3910.0,  3),                                &
        !  74 : branching ratio for R80
        !       CH3O2 + CH3COO2 = 2 HCHO
        RR( 4.4e-12,  0.0,   365.0,  3),                                &
        !  75 : OH + CH3CHO = CH3COO2 + H2O
        RR( 2.7e-28, -7.1,     0.0,  2)                                 &
        !  76 : CH3COO2 + NO2 + M = CH3COO2NO2 + M  klo
        /)

      rk(77:95) = (/                                                    &
        RR( 1.2e-11, -0.9,     0.0,  2),                                &
        !  77 : CH3COO2 + NO2 + M = CH3COO2NO2 + M  khi
        RR( 4.9e-03,  0.0,-12100.0,  3),                                &
        !  78 : CH3COO2NO2 + M = CH3COO2 + NO2 + M  klo
        RR( 2.0e-11,  0.0,     0.0,  1),                                &
        !  79 : CH3COO2 + NO = CH3O2 + CO2 + NO2
        RR( 1.1e-11,  0.0,     0.0,  1),                                &
        !  80 : CH3O2 + CH3COO2 = HCHO + HO2 + CH3O2 + CO2 + O2
        RR(7.9e-13,   2.0,   300.0,  4),                                &
        !  81 : OH + n-C4H10 (+O2) = s-C4H9O2 + H2O
        RR( 5.4e+16,  0.0,-13830.0,  3),                                &
        !  82 : CH3COO2NO2 + M = CH3COO2 + NO2 + M  khi
        RR( 2.54e-12, 0.0,   360.0,  3),                                &
        !  83 : s-C4H9O2 + NO = CH3COC2H5 + HO2 + NO2
        RR( 2.5e-13,  0.0,     0.0,  1),                                &
        !  84 : s-C4H9O2 + CH3O2 = CH3COC2H5+2 HO2 +HCHO+O2
        zero_rate,                                                      &
        RR(1.3e-12,   0.0,   -25.0,  3),                                &
        !  86 : OH + CH3COC2H5 (+O2) = MEKO2 + H2O
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR( 8.8e-12,  0.0, -1320.0,  3),                                &
        !  89 : CH3COCH3 + OH = CH3COCH2O2 + H2O  See R94
        RR( 6.4e-14,  0.0,     0.0,  1),                                &
        !  90 : C2H5O2 + C2H5O2 = C2H5O + C2H5O + O2
        RR( 2.9e-12,  0.0,   500.0,  3),                                &
        !  91 : CH3COO2 + CH3COO2 = 2 CH3O2 + 2 CO2 + O2
        RR( 7.6e-12,  0.0,  -585.0,  3),                                &
        !  92 : OH + C3H8 = i-C3H7O2 + H2O
        RR( 2.7e-12 , 0.0,   360.0,  3),                                &
        !  93 : i-C3H7O2 + NO = NO2 + HO2 + CH3COCH3
        RR( 1.7e-14,  0.0,   423.0,  3),                                &
        !  94 : CH3COCH3 + OH = CH3COCH2O2 + H2O  2nd term for R89
        RR(2.45e-12,  0.0,   360.0,  3)                                 &
        !  95 : CH3COCH2O2 + NO = NO2 + CH3COO2 + HCHO
        /)

      rk(96:114) = (/                                                   &
        RR( 3.8e-12,  0.0,     0.0,  1),                                &
        !  96 : CH3COCH2O2 + CH3O2 = HO2 + 2 HCHO + CH3COO2
        RR( 4.0e-14,  0.0,     0.0,  1),                                &
        !  97 : i-C3H7O2 + CH3O2 = HCHO + 2 HO2 + CH3COCH3
        RR( 9.5e-13,  0.0,  -650.0,  3),                                &
        !  98 : OH + PAN = NO3 + HCHO
        RR( 3.8e-13,  0.0,   900.0,  3),                                &
        !  99 : HO2 + C2H5O2 = C2H5OOH + O2
        RR( 8.0e-12,  0.0,     0.0,  1),                                &
        ! 100 : OH + C2H5OOH = CH3CHO + OH
        RR(1.51e-13,  0.0,  1300.0,  3),                                &
        ! 101 : HO2 + C3H7O2 = i-C3H7OOH + O2
        RR(1.66e-11,  0.0,     0.0,  1),                                &
        ! 102 : OH + i-C3H7OOH = CH3COCH3 + OH
        RR(1.82e-13,  0.0,  1300.0,  3),                                &
        ! 103 : HO2 + s-C4H9O2 = s-C4H9OOH + O2
        RR(2.15e-11,  0.0,     0.0,  1),                                &
        ! 104 : OH + s-C4H9OOH = CH3COC2H5 + OH
        RR(2.54e-12,  0.0,   360.0,  3),                                &
        ! 105 : NO + CH3COCH(O2)CH3 = CH3COO2 + CH3CHO + NO2
        RR( 8.8e-13,  0.0,     0.0,  1),                                &
        ! 106 : CH3O2 + CH3COCH(O2)CH3 =
        !                    HCHO + HO2 + CH3COO2 + CH3CHO + O2
        zero_rate,                                                      &
        RR( 8.6e-29, -3.1,     0.0,  2),                                &
        ! 108 : OH + C2H4 (+O2) = HOC2H4O2  klo
        RR( 9.0e-12,-0.85,     0.0,  2),                                &
        ! 109 : OH + C2H4 (+O2) = HOC2H4O2  khi
        RR( 9.0e-12,  0.0,     0.0,  1),                                &
        ! 110 : HOC2H4O2 + NO = 2 HCHO + HO2 + NO2
        RR( 2.0e-12,  0.0,     0.0,  1),                                &
        ! 111 : CH3O2 + HOC2H4O2 = 3 HCHO + 2 HO2 + O2
        RR( 1.2e-14,  0.0, -2630.0,  3),                                &
        ! 112 : O3 + C2H4 = HCHO + 0.47 CH2O2 + 0.31 CO
        !                 + 0.22 CO2 + 0.31 H2O + 0.13 H2 + 0.20 HO2  
        zero_rate,                                                      &
        zero_rate /)                                                                              

      rk(115:133) = (/                                                  &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR(2.75E-15,  0.0, -1878.0,  3),                                &
         ! 123 : O3 + C3H6 = HCHO + 0.30 CH4 + 0.40 CO
         !    + 0.60 CO2 + 0.28 OH + 0.12 CH3OH + 0.30 HO2 + 0.58 CH3O2    
        RR(2.75E-15,  0.0, -1878.0,  3),                                &
         ! 124 : O3 + C3H6 = CH3CHO + 0.24 H2 + 0.58 CO
         !          + 0.42 CO2 + 0.58 H2O +0.18 HO2     
        RR( 8.0e-27, -3.5,     0.0,  2),                                &
         ! 125 : OH + C3H6 (+O2) = HOC3H6O2  klo see R130
        RR(2.54E-12,  0.0,   360.0,  3),                                &
         ! 126 : HOC3H6O2 + NO = HCHO + HO2 + CH3CHO + NO2
        RR( 6.0e-13,  0.0,     0.0,  1),                                &
         ! 127 : CH3O2 + HOC3H6O2 = 2 HCHO + 2 HO2 + CH3CHO + O2
        RR(7.86e-15,  0.0, -1913.0,  3),                                &
         ! 128 : O3 + C5H8 = MVK + 0.78 CO + 0.22 CH2OO + 0.27 HO2
         !                       + 0.27 OH
        RR(7.56e-16,  0.0, -1521.0,  3),                                &
         ! 129 : O3 + MVK = MGLYOX + 0.76 CO + 0.24 CH2OO + 0.36 HO2
         !                         + 0.36 OH
        RR( 3.0e-11,  -1.0,    0.0,  2),                                &
         ! 130 : OH + C3H6 (+O2) = HOC3H6O2  khi see R125
        RR( 5.2e-13,  0.0,   980.0,  3),                                &
         ! 131 : CH3CO3 + HO2 = 0.3 O3 + 0.8CH3O2 + 0.2CH3COO2
        RR( 5.0e-13,  0.0,     0.0,  1),                                &
         ! 132 : HOIPO2 + CH3O2 = 2HO2 + HCHO + MVK
        RR( 2.0e-12,  0.0,     0.0,  1)                                 &
         ! 133 : HOMVKO2+ CH3O2 = 2HO2 + HCHO + MGLYOX
        /)

      rk(134:152) = (/                                                  &
        RR(2.45e-13,  0.0,  1250.0,  3),                                &
         ! 134 : HOIPO2 + HO2 = ISOOH + O2
        RR( 4.2e-11,  0.0,     0.0,  1),                                &
         ! 135 : ISOOH  + OH = MVK + HCHO + OH
        RR(2.23e-13,  0.0,  1250.0,  3),                                &
         ! 136 : HOMVKO2 + HO2 = MVKOOH + O2
        RR(5.77e-11,  0.0,     0.0,  1),                                &
         ! 137 : MVKOOH + OH = MGLYOX + HCHO + OH
        RR(1.72e-11,  0.0,     0.0,  1),                                &
         ! 138 : MGLYOX + OH = CH3COO2 + CO
        RR(1.14e-11,  0.0,     0.0,  1),                                &
         ! 139 : GLYOX + OH = HO2 + 2 CO
        zero_rate,                                                      &
        RR(2.54e-11,  0.0,   410.0,  3),                                &
         ! 141 : OH + C5H8 = HOIPO2 + H2O
        RR(2.08e-12,  0.0,   180.0,  3),                                &
         ! 142 : HOC5H8O2 + NO = MVK + HO2 + HCHO + NO2
        RR(4.13e-12,  0.0,   452.0,  3),                                &
         ! 143 : OH + MVK = HOMVKO2 + H2O
        RR(2.5e-12,   0.0,   360.0,  3),                                &
         ! 144 : HOMVKO2 + NO = CH3COCHO + CH2O + HO2 + NO2
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR( 5.7e-12,  0.0, -4426.0,  3),                                &
         ! 151 : NO3 + C2H6 = C2H5O2 + HNO3
        RR(2.8e-12,   0.0, -3280.0,  3)                                 &
         ! 152 : NO3 + n-C4H10 = s-C4H9O2 + HNO3
        /)

      rk(153:172) = (/                                                  &
        RR(3.3e-12,   2.0, -2880.0,  4),                                &
         ! 153 : NO3 + C2H4 = C2H4NO3 = CH2(NO3)CHO + HO2
        RR(4.95e-12,  0.0,     0.0,  1),                                &
         ! 154 : CH2(NO3)CHO + OH = HCHO + NO2 + CO2
        RR(4.59e-13,  0.0, -1156.0,  3),                                &
         ! 155 : NO3 + C3H6 = C3H6NO3 = CH3CH(NO3)CHO + HO2
        RR(5.25e-12,  0.0,     0.0,  1),                                &
         ! 156 : CH3CH(NO3)CHO + OH = CH3CHO + NO2 + CO2
        zero_rate,                                                      &
        RR( 1.4E-12,  0.0, -1860.0,  3),                                &
         ! 158 : NO3 + CH3CHO = CH3COO2 + HNO3
        RR(3.03e-12,  0.0,  -446.0,  3),                                &
         ! 159 : NO3 + C5H8 = (NO3)C4H6CHO + HO2
        RR(4.16e-11,  0.0,     0.0,  1),                                &
         ! 160 : (NO3)C4H6CHO + OH = MVK + NO2
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        RR(1.36e-11,  0.0,     0.0,  1),                                &
         ! 170 : OH + oXYL = HO2 + 0.8 MEMALD + 0.8 MGLYOX
        RR( 1.0e-11,  0.0,     0.0,  1),                                &
         ! 171 : OXYL1 + NO2 = ORGNIT
        RR( 5.6e-11,  0.0,     0.0,  1)                                 &
         ! 172 : OH + MEMALD = MEMALD1
        /)

      rk(173:192) = (/                                                  &
        RR( 2.54e-12, 0.0,   360.0,  3),                                &
         ! 173 : MEMALD1 + NO = HO2 + CH3COCHO + CHOCHO + NO2
        RR(1.18E-12,  0.0,   338.0,  3),                                &
         ! 174 : OH + TOLUENE = MEMALD + GLYOX + HO2
        RR( 3.6e-13,  0.0,     0.0,  1),                                &
         ! 175 : OH + TOLUENE = TOLP1
        RR(1.36e-11,  0.0,     0.0,  1),                                &
         ! 176 : OH + oXYL    = OXYL1
        RR( 2.7e-12,  0.0,     0.0,  1),                                &
         ! 177 : OH + ORGNIT  = MEMALD + GLYOX + NO2
        RR( 7.0e-14,  0.0,     0.0,  1),                                &
         ! 178 : NO3 + ORGNIT = MEMALD + GLYOX + 2NO2
        zero_rate,                                                      &
        RR( 2.5e-13,  0.0,  1300.0,  3),                                &
         ! 180 : OXYL1 + HO2 = MGLYOX + MEMALD
        RR( 1.0e-13,  0.0,     0.0,  1),                                &
         ! 181 : MEMALD1 + CH3O2 = 2HO2 + HCHO + MGLYOX + GLYOX     
        RR( 1.0e-11,  0.0,     0.0,  1),                                &
         ! 182 : HO2 + TOLP1 = MEMALD + GLYOX + OH
        RR( 1.0e-11,  0.0,     0.0,  1),                                &
         ! 183 : NO2 + TOLP1 = ORGNIT
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate,                                                      &
        zero_rate /)
!
!
! Check that reaction types are correct
        DO i = 1, nr
          r = rk(i)%r
          ar = rk(i)%a
          n = rk(i)%n
          ea = rk(i)%e
          IF ((r == 1) .AND. (n /= 0.0 .OR. ea /= 0.0)) THEN
            cmessage = 'Faulty rate constant'

            CALL EREPORT('UKCA_CHEMCO', i, cmessage)
          END IF
          IF ((r == 2) .AND. (n == 0.0 .OR. ea /= 0.0)) THEN
            cmessage = 'Faulty rate constant'

            CALL EREPORT('UKCA_CHEMCO', i, cmessage)
          END IF
          IF ((r == 3) .AND. (n /= 0.0 .OR. ea == 0.0)) THEN
            cmessage = 'Faulty rate constant'

            CALL EREPORT('UKCA_CHEMCO', i, cmessage)
          END IF
          IF ((r == 4) .AND. (n == 0.0 .OR. ea == 0.0)) THEN
            cmessage = 'Faulty rate constant'

            CALL EREPORT('UKCA_CHEMCO',i,cmessage)
          END IF
        END DO
!
        first = .FALSE.
!
      END IF  ! If first = .TRUE.
!
      rc = 0.0
!
      DO i = 1, nr
        r  = rk(i)%r
        ar = rk(i)%a
        n  = rk(i)%n
        ea = rk(i)%e
        IF (r == 1) THEN
          rc(:,i) = ar
        ELSE IF (r == 2) THEN
          rc(:,i) = ar * ((tc(:)*tfc) ** n)
        ELSE IF (r == 3) THEN
          rc(:,i) = ar * EXP(ea / tc(:))
        ELSE IF (r == 4) THEN
          rc(:,i) = ar * ((tc(:)*tfc) ** n) * EXP(ea / tc(:))
        END IF
      END DO
!
! Pressure-dependent reactions
!
      DO k = 1, npdep
        nklo = pdep(k)%klo                  ! Reaction no. for klo
        nkhi = pdep(k)%khi                  ! Reaction no. for khi
        nres = pdep(k)%res                  ! Reaction no. for result
        f    = pdep(k)%fc
        IF (f > 0.01) THEN
          fc = f
          DO j = 1, n_pnts
            rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
            brn        = 0.75 - 1.27 * LOG10(fc)
            fac1       = rc(j,nklo) / rc(j,nkhi)    ! khi
            fac2       = rc(j,nklo) / (1.0 + fac1)
            fac3       = 1.0 + ((LOG10(fac1) / brn) ** 2)
            rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
          END DO
        ELSE ! In stochem_chemistry.mf77 f was set to zero for R29
             ! to flag this reaction and recalculate f later.
             ! Following Atkinson (2004), f changed to 0.35, but 
             ! we leave the code below in case of any new change.    
          DO j = 1, n_pnts
            fc         = EXP(-tc(j)/250.0) + EXP(-1050.0/tc(j)) ! For R29 
            rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
            brn        = 0.75 - 1.27 * LOG10(fc)
            fac1       = rc(j,nklo) / rc(j,nkhi)    ! khi
            fac2       = rc(j,nklo) / (1.0 + fac1)
            fac3       = 1.0 + ((LOG10(fac1) / brn) ** 2)
            rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
          END DO
        END IF
      END DO
!
! Do remaining calculations for certain reactions. This section
! needs to be edited by hand
!
! Reaction   7 : O(1D) + M
      rc(:, 7) = (o2frac*rc(:,6) + n2frac*rc(:,7)) * m(:)
!
! Reaction  35 : OH + HNO3
      fac1v(:) = rc(:,51) * m(:)
      rc(:,35) = rc(:,35) + fac1v(:) / (1.0 + fac1v(:) / rc(:,50))
!
! Reaction  36 : HO2 + HO2 (+ M, H2O)
      rc(:,36) = (rc(:,36) + rc(:,37) * m(:)) *                       &
           (1.0 + rc(:,38) * h2o(:))
!
! Reactions  61 and 62 : CH3O2 + CH3O2
! (use both branches in UKCA_DERIV)
      rc(:,62) = rc(:,62) - rc(:,61)
!
! Reaction  70 : OH + CO
      rc(:,70) = rc(:,70) + (rc(:,69) * m(:))
!
! Reactions  74 and 80 : CH3O2 + CH3COO2
! (consider both of them in UKCA_DERIV)
      fac1v(:) = rc(:,74) / (1.0 + rc(:,74))
      rc(:,74) = rc(:,80) * (1.0-fac1v(:))  ! -> 2 HCHO + O2
      rc(:,80) = rc(:,80) * fac1v(:)        ! -> HCHO + HO2 + CH3O2 + CO2 + O2
!
! Reaction 89 and 94: CH3COCH3 + OH = CH3COCH2O2 + H2O
      rc(:,94) = rc(:,89) + rc(:,94)
!
! Reaction   1 : O + O2 (+ M) = O3 + M
      rc(:, 1) = rc(:, 1) * m(:) * o2(:)
!

      IF (lhook) CALL dr_hook('UKCA_CHEMCO_RAQ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CHEMCO_RAQ
