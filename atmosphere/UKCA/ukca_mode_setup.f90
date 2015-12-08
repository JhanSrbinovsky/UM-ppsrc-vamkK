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
!    Module to store MODE setup arrays after call to UKCA_MODE_SETUP_ALL.
!    Contains public subroutines:
!      UKCA_MODE_TRACER_INDICES
!      UKCA_MODE_IMSCAVCOFF
!      UKCA_MODE_ALLCP_4MODE
!      UKCA_MODE_SUSS_4MODE
!      UKCA_MODE_SUSSBCOC_4MODE
!      UKCA_MODE_SUSSBCOC_5MODE
!      UKCA_MODE_SUSSBCOCSO_5MODE
!      UKCA_MODE_SUSSBCOCSO_4MODE
!      UKCA_MODE_DUonly_2MODE
!      UKCA_MODE_DUonly_3MODE (needs to be added at some point)
!      UKCA_MODE_SUSSBCOCDU_7MODE
!      UKCA_MODE_SUSSBCOCDU_4MODE
!    which define modes and components for different components/modes setup.
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
      MODULE UKCA_MODE_SETUP
! ---------------------------------------------------------------------|
!+ Module to contain modes and components
!
! Description:
! To allow use throughout UM, module stores MODE setup arrays after
! call to UKCA_MODE_SETUP_ALL.
!
! Note: currently code is hard-coded so that ordering of modes must
! 1) nucln, 2)   soluble Aitken, 3)   soluble accum, 4)   soluble coarse
!           5) insoluble Aitken, 6) insoluble accum, 7) insoluble coarse
!
! ---------------------------------------------------------------------|

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
      SAVE

      INTEGER, PARAMETER :: nmodes=7  ! No of modes
      INTEGER, PARAMETER :: ncp=6     ! No of components
      INTEGER, PARAMETER :: ncation=3 ! No possible cation species
      INTEGER, PARAMETER :: nanion =4 ! No possible anion species
!
      INTEGER, PARAMETER :: CP_SU=1  ! Index to store SO4    cpt
      INTEGER, PARAMETER :: CP_BC=2  ! Index to store BC     cpt
      INTEGER, PARAMETER :: CP_OC=3  ! Index to store 1st OC cpt
      INTEGER, PARAMETER :: CP_CL=4  ! Index to store NaCl   cpt
      INTEGER, PARAMETER :: CP_DU=5  ! Index to store dust   cpt
      INTEGER, PARAMETER :: CP_SO=6  ! Index to store 2nd OC cpt

! Mode switches (1=on, 0=0ff)
      INTEGER, DIMENSION(nmodes) :: mode_choice
! Component switches (1=on, 0=off)
      INTEGER, DIMENSION(ncp)    :: component_choice
! Components that are soluble
      INTEGER, DIMENSION(ncp)    :: soluble_choice
! Components allowed in each mode (must be consistent with coag_mode)
      INTEGER, DIMENSION(nmodes,ncp) :: component_mode
! Modes resulting when two modes coagulate
      INTEGER, DIMENSION(nmodes,nmodes) :: coag_mode
! Specify which modes are soluble
      INTEGER, DIMENSION(nmodes) :: modesol
!
! Variables for impaction scavenging (as in Pringle, 2006 PhD thesis)
      INTEGER, PARAMETER  :: NCOLL=20 ! # of columns in LUT (aer. bins)
      INTEGER, PARAMETER  :: NROW =19 ! # of rows in LUT (raindrop bins)

! Tracer indices 
      INTEGER :: ii_nd(nmodes)        ! indices in mode_tracers of NCONC
      INTEGER :: ii_md(nmodes,ncp)    ! indices in mode_tracers of CPTMMR

! Molar masses of components (kg mol-1)
      REAL, DIMENSION(ncp)    :: mm
! Mass density of components (kg m^-3)
      REAL, DIMENSION(ncp)    :: rhocomp
! Number of dissociating ions in soluble components
      REAL, DIMENSION(ncp)    :: no_ions
! Lower size limits of geometric mean radius for each mode
      REAL, DIMENSION(nmodes) :: fracbcem
! Fraction of bc ems to go into each mode
      REAL, DIMENSION(nmodes) :: fracocem
! Fraction of oc ems to go into each mode
      REAL, DIMENSION(nmodes) :: ddplim0
! Upper size limits of geometric mean radius for each mode
      REAL, DIMENSION(nmodes) :: ddplim1
! Mid-point of size mode (m)
      REAL, DIMENSION(nmodes) :: DDPMID
! Mid-point masses for initial radius grid
      REAL, DIMENSION(nmodes) :: MMID
! Lo-interf masses for initial radius grid
      REAL, DIMENSION(nmodes) :: MLO
! Hi-interf masses for initial radius grid
      REAL, DIMENSION(nmodes) :: MHI
! Fixed geometric standard deviation for each mode
      REAL, DIMENSION(nmodes) :: sigmag
! EXP((9/2)*LOG^2(SIGMA_G))
      REAL, DIMENSION(nmodes) :: x
! Threshold for number in mode to carry out calculations
      REAL, DIMENSION(nmodes) :: NUM_EPS
! Initial fractions of mass in each mode among components
      REAL, DIMENSION(nmodes,ncp) :: mfrac_0

      REAL, DIMENSION(nrow)       :: raddrop  ! raindrop bins
      REAL, DIMENSION(ncoll,nrow) :: colleff4 ! collision efficiency

! Mode names
      CHARACTER(len=7),DIMENSION(nmodes) :: mode_names
! Component names
      CHARACTER(len=7),DIMENSION(ncp)    :: component_names

! Modes (T/F)
      LOGICAL, DIMENSION(nmodes)     :: mode
! Components set in each mode (T/F)
      LOGICAL, DIMENSION(nmodes,ncp) :: component
! Components which are soluble (T/F)
      LOGICAL, DIMENSION(ncp)        :: soluble

      CONTAINS

      SUBROUTINE UKCA_MODE_TRACER_INDICES

! To set up the tracer indices for number and mass

      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook( &
              'UKCA_MODE_SETUP:UKCA_MODE_TRACER_INDICES',               &
                               zhook_in,zhook_handle)

! .. Template which matches the
! .. 31 prognostic tracers as indexed in mode_tracers to
! .. the ND and MD values (as set using variables ITRA)

      II_ND                =(/ 1, 3, 7,13,19,22,24/)

      II_MD(1:NMODES,CP_SU)=(/ 2, 4, 8,14, 0, 0, 0/)
      II_MD(1:NMODES,CP_BC)=(/ 0, 5, 9,15,20, 0, 0/)
      II_MD(1:NMODES,CP_OC)=(/26, 6,10,16,21, 0, 0/)
      II_MD(1:NMODES,CP_CL)=(/ 0,27,11,17, 0, 0, 0/)
      II_MD(1:NMODES,CP_DU)=(/ 0, 0,12,18, 0,23,25/)
      II_MD(1:NMODES,CP_SO)=(/28,29,30,31, 0, 0, 0/)

      IF (lhook) CALL dr_hook( &
              'UKCA_MODE_SETUP:UKCA_MODE_TRACER_INDICES',               &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_TRACER_INDICES

      SUBROUTINE UKCA_MODE_IMSCAVCOFF

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_IMSCAVCOFF',   &
                               zhook_in,zhook_handle)
      RADDROP =(/   1.0, 1.587,  2.52,   4.0,  6.35, 10.08,             &
                   16.0,  25.4, 40.32,  64.0, 101.6, 161.3,             &
                  256.0, 406.4, 645.1,1024.0,1625.0,2580.0,4096.0/)

      COLLEFF4( 1,1:NROW)=(/ 0.522E+05,0.139E+05,0.328E+04,0.775E+03,   &
                             0.183E+03,0.432E+02,0.102E+02,0.291E+01,   &
                             0.108E+01,0.439E+00,0.201E+00,0.110E+00,   &
                             0.633E-01,0.366E-01,0.815E-02,0.168E-02,   &
                             0.394E-03,0.591E-04,0.132E-04/)

      COLLEFF4( 2,1:NROW)=(/ 0.126E+05,0.373E+04,0.985E+03,0.260E+03,   &
                             0.687E+02,0.182E+02,0.480E+01,0.150E+01,   &
                             0.536E+00,0.224E+00,0.107E+00,0.608E-01,   &
                             0.364E-01,0.236E-01,0.731E-02,0.167E-02,   &
                             0.395E-03,0.592E-04,0.133E-04/)

      COLLEFF4( 3,1:NROW)=(/ 0.445E+04,0.139E+04,0.390E+03,0.110E+03,   &
                             0.308E+02,0.864E+01,0.243E+01,0.783E+00,   &
                             0.270E+00,0.117E+00,0.564E-01,0.325E-01,   &
                             0.203E-01,0.136E-01,0.595E-02,0.166E-02,   &
                             0.396E-03,0.594E-04,0.133E-04/)

      COLLEFF4( 4,1:NROW)=(/ 0.259E+04,0.810E+03,0.227E+03,0.639E+02,   &
                             0.179E+02,0.503E+01,0.141E+01,0.440E+00,   &
                             0.146E+00,0.662E-01,0.309E-01,0.177E-01,   &
                             0.115E-01,0.728E-02,0.446E-02,0.164E-02,   &
                             0.399E-03,0.597E-04,0.134E-04/)

      COLLEFF4( 5,1:NROW)=(/ 0.196E+04,0.602E+03,0.166E+03,0.457E+02,   &
                             0.126E+02,0.346E+01,0.953E+00,0.280E+00,   &
                             0.927E-01,0.410E-01,0.183E-01,0.990E-02,   &
                             0.655E-02,0.444E-02,0.389E-02,0.164E-02,   &
                             0.402E-03,0.604E-04,0.135E-04/)

      COLLEFF4( 6,1:NROW)=(/ 0.192E+04,0.569E+03,0.151E+03,0.401E+02,   &
                             0.106E+02,0.282E+01,0.749E+00,0.206E+00,   &
                             0.689E-01,0.280E-01,0.119E-01,0.580E-02,   &
                             0.384E-02,0.304E-02,0.393E-02,0.165E-02,   &
                             0.406E-03,0.612E-04,0.138E-04/)

      COLLEFF4( 7,1:NROW)=(/ 0.208E+04,0.604E+03,0.156E+03,0.405E+02,   &
                             0.105E+02,0.271E+01,0.703E+00,0.185E+00,   &
                             0.572E-01,0.217E-01,0.903E-02,0.462E-02,   &
                             0.280E-02,0.208E-02,0.398E-02,0.168E-02,   &
                             0.415E-03,0.628E-04,0.142E-04/)

      COLLEFF4( 8,1:NROW)=(/ 0.233E+04,0.636E+03,0.162E+03,0.410E+02,   &
                             0.104E+02,0.264E+01,0.671E+00,0.174E+00,   &
                             0.513E-01,0.184E-01,0.747E-02,0.384E-02,   &
                             0.222E-02,0.161E-02,0.408E-02,0.173E-02,   &
                             0.430E-03,0.657E-04,0.149E-04/)

      COLLEFF4( 9,1:NROW)=(/ 0.235E+04,0.659E+03,0.165E+03,0.412E+02,   &
                             0.103E+02,0.257E+01,0.643E+00,0.168E+00,   &
                             0.490E-01,0.168E-01,0.661E-02,0.326E-02,   &
                             0.188E-02,0.140E-02,0.422E-02,0.180E-02,   &
                             0.452E-03,0.698E-04,0.160E-04/)

      COLLEFF4(10,1:NROW)=(/ 0.165E+04,0.457E+03,0.112E+03,0.277E+02,   &
                             0.680E+01,0.167E+01,0.412E+00,0.106E+00,   &
                             0.304E-01,0.999E-02,0.386E-02,0.186E-02,   &
                             0.124E-02,0.140E-02,0.447E-02,0.193E-02,   &
                             0.491E-03,0.771E-04,0.179E-04/)

      COLLEFF4(11,1:NROW)=(/ 0.899E+03,0.246E+03,0.597E+02,0.145E+02,   &
                             0.352E+01,0.856E+00,0.208E+00,0.524E-01,   &
                             0.145E-01,0.466E-02,0.179E-02,0.860E-03,   &
                             0.719E-03,0.165E-02,0.486E-02,0.213E-02,   &
                             0.554E-03,0.891E-04,0.211E-04/)

      COLLEFF4(12,1:NROW)=(/ 0.117E+04,0.326E+03,0.807E+02,0.200E+02,   &
                             0.496E+01,0.123E+01,0.305E+00,0.777E-01,   &
                             0.219E-01,0.720E-02,0.281E-02,0.137E-02,   &
                             0.941E-03,0.330E-02,0.563E-02,0.255E-02,   &
                             0.686E-03,0.116E-03,0.283E-04/)

      COLLEFF4(13,1:NROW)=(/ 0.130E+04,0.371E+03,0.938E+02,0.237E+02,   &
                             0.601E+01,0.152E+01,0.385E+00,0.979E-01,   &
                             0.276E-01,0.926E-02,0.364E-02,0.173E-02,   &
                             0.101E-02,0.406E-02,0.694E-02,0.327E-02,   &
                             0.930E-03,0.167E-03,0.429E-04/)

      COLLEFF4(14,1:NROW)=(/ 0.118E+04,0.333E+03,0.842E+02,0.213E+02,   &
                             0.537E+01,0.136E+01,0.342E+00,0.876E-01,   &
                             0.250E-01,0.841E-02,0.330E-02,0.153E-02,   &
                             0.801E-03,0.260E-02,0.973E-02,0.490E-02,   &
                             0.152E-02,0.303E-03,0.842E-04/)

      COLLEFF4(15,1:NROW)=(/ 0.774E+03,0.223E+03,0.572E+02,0.147E+02,   &
                             0.378E+01,0.970E+00,0.249E+00,0.636E-01,   &
                             0.180E-01,0.606E-02,0.238E-02,0.110E-02,   &
                             0.658E-03,0.107E-02,0.167E-01,0.940E-02,   &
                             0.335E-02,0.791E-03,0.249E-03/)

      COLLEFF4(16,1:NROW)=(/ 0.372E+01,0.177E+01,0.781E+00,0.345E+00,   &
                             0.153E+00,0.675E-01,0.299E-01,0.130E-01,   &
                             0.624E-02,0.346E-02,0.179E-02,0.142E-02,   &
                             0.164E-02,0.401E-02,0.413E-01,0.277E-01,   &
                             0.124E-01,0.389E-02,0.152E-02/)

      COLLEFF4(17,1:NROW)=(/ 0.234E-18,0.108E-16,0.705E-15,0.462E-13,   &
                             0.302E-11,0.198E-09,0.129E-07,0.844E-06,   &
                             0.568E-04,0.386E-02,0.286E-01,0.448E-01,   &
                             0.569E-01,0.859E-01,0.183E+00,0.165E+00,   &
                             0.108E+00,0.536E-01,0.295E-01/)

      COLLEFF4(18,1:NROW)=(/ 0.902E-37,0.794E-33,0.160E-28,0.324E-24,   &
                             0.655E-20,0.132E-15,0.267E-11,0.540E-07,   &
                             0.816E-03,0.482E-01,0.245E+00,0.372E+00,   &
                             0.436E+00,0.473E+00,0.493E+00,0.497E+00,   &
                             0.414E+00,0.299E+00,0.225E+00/)

      COLLEFF4(19,1:NROW)=(/ 0.275E-30,0.136E-26,0.146E-22,0.156E-18,   &
                             0.167E-14,0.178E-10,0.191E-06,0.202E-02,   &
                             0.203E+00,0.427E+00,0.586E+00,0.669E+00,   &
                             0.708E+00,0.730E+00,0.746E+00,0.738E+00,   &
                             0.679E+00,0.588E+00,0.520E+00/)

      COLLEFF4(20,1:NROW)=(/ 0.136E-33,0.238E-29,0.102E-24,0.436E-20,   &
                             0.186E-15,0.797E-11,0.341E-06,0.143E-01,   &
                             0.722E+00,0.805E+00,0.869E+00,0.902E+00,   &
                             0.915E+00,0.932E+00,0.104E+01,0.927E+00,   &
                             0.904E+00,0.871E+00,0.842E+00/)

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_IMSCAVCOFF',   &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_IMSCAVCOFF
      SUBROUTINE UKCA_MODE_ALLCP_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components with all components
!+ switched on but only 4 modes used.
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Mode names
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_ALLCP_4MODE',  &
                               zhook_in,zhook_handle)
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,1,0/)
! ***n.b. in above have kept all cpts on (not SO) for UM test***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7 but sgacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust   so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
      fracocem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_ALLCP_4MODE',  &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_ALLCP_4MODE
      SUBROUTINE UKCA_MODE_SUSS_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate and sea-salt only in 4 modes.
!+ Uses 10 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSS_4MODE',   &
                              zhook_in,zhook_handle)
! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,0,0,1,0,0/)
! *** n.b. only have h2so4 and nacl cpts on for this setup ***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)
!      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-8,1.0e-7,1.0e-6/)
!      ddplim1=(/1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigmacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
      fracocem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSS_4MODE',   &
                              zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSS_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCDU_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ SO4, sea-salt, bc, oc (secondary & primary combined) & du in 4 modes.
!+ Uses 19 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_4MODE',  &
                              zhook_in,zhook_handle)
! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,1,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
      fracocem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_4MODE',  &
                              zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOCDU_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCDU_7MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ SO4, sea-salt, bc, oc (secondary & primary combined) & du in 7 modes.
!+ Uses 26 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_7MODE', &
                               zhook_in,zhook_handle)

! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,1,1,1/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,1,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8,1.0e-8/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
      fracocem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCDU_7MODE', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOCDU_7MODE
      SUBROUTINE UKCA_MODE_SUSSBCOC_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc & oc (secondary & primary combined) in 4 modes.
!+ Uses 17 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Mode names
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_4MODE', &
                               zhook_in,zhook_handle)
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,0,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
      fracocem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_4MODE', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOC_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOC_5MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc & oc (secondary & primary combined) in 5 modes.
!+ Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Mode names
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_5MODE', &
                               zhook_in,zhook_handle)
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,1,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,0,0/)
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,1,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! Assume other components have same mass density as H2SO4
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
      fracocem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOC_5MODE', &
                              zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOC_5MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCSO_4MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!+ Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_4MODE', &
                              zhook_in,zhook_handle)
! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,0,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,0,1/) ! ***all cpts on except dust***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
      fracocem=(/0.0,1.0,0.0,0.0,0.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into   soluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_4MODE', &
                              zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOCSO_4MODE
      SUBROUTINE UKCA_MODE_SUSSBCOCSO_5MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!+ Uses 23 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_5MODE', &
                               zhook_in,zhook_handle)
! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/1,1,1,1,1,0,0/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/1,1,1,1,0,1/) ! ***all cpts on except dust***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
      fracocem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO
      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_SUSSBCOCSO_5MODE', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_SUSSBCOCSO_5MODE
      SUBROUTINE UKCA_MODE_DUonly_2MODE
! ---------------------------------------------------------------------|
!+ Subroutine to define modes and components for version with
!+ only du cpt in 2 (insoluble) modes.
!+ Uses  4 aerosol tracers
! ---------------------------------------------------------------------|
      USE UKCA_CONSTANTS,   ONLY: avc, rhosul, mmsul, ppi
      IMPLICIT NONE

      INTEGER :: imode
      INTEGER :: icp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_DUONLY_2MODE', &
                               zhook_in,zhook_handle)

! Mode names
      mode_names=(/'nuc_sol','ait_sol','acc_sol','cor_sol',             &
                             'ait_ins','acc_ins','cor_ins'/)
! Mode switches (1=on, 0=0ff)
      mode_choice=(/0,0,0,0,0,1,1/)
! Specify which modes are soluble
      modesol=(/1,1,1,1,0,0,0/)
! Component names
      component_names=                                                  &
       (/'h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org'/)
! Component switches (1=on, 0=off)
      component_choice=(/0,0,0,0,1,0/) ! ***all cpts on except dust***
! Components that are soluble
      soluble_choice=(/1,0,0,1,0,0/)
! Components allowed in each mode (must be consistent with coag_mode)
      component_mode(1,1:ncp)=(/1,0,0,0,0,1/)       !allowed in nuc_sol
      component_mode(2,1:ncp)=(/1,1,1,0,0,1/)       !allowed in ait_sol
      component_mode(3,1:ncp)=(/1,1,1,1,1,1/)       !allowed in acc_sol
      component_mode(4,1:ncp)=(/1,1,1,1,1,1/)       !allowed in cor_sol
      component_mode(5,1:ncp)=(/0,1,1,0,0,0/)       !allowed in ait_ins
      component_mode(6,1:ncp)=(/0,0,0,0,1,0/)       !allowed in acc_ins
      component_mode(7,1:ncp)=(/0,0,0,0,1,0/)       !allowed in cor_ins

! Specify size limits of geometric mean radius for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
      ddplim0=(/1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6/)
      ddplim1=(/1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5/)

! Specify fixed geometric standard deviation for each mode
      sigmag=(/1.59,1.59,1.40,2.0,1.59,1.59,2.0/) ! to match M7, but sigacc=1.4
!      sigmag=(/1.59,1.59,1.59,2.0,1.59,1.59,2.0/) ! to match M7
!      sigmag=(/1.59,1.45,1.4,2.0,1.45,1.4,2.0/) ! to match UM scheme

      DO imode=1,nmodes
        x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
      NUM_EPS=(/1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5/)
!!      NUM_EPS=(/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-2,1.0e-3,1.0e-4/)
!!      NUM_EPS=(/1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3,1.0e-3/)
!
      DO imode=1,nmodes
       ddpmid(imode)=                                                   &
          EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
       mmid(imode)=                                                     &
          (PPI/6.0)*(ddpmid(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mlo (imode)=                                                     &
          (PPI/6.0)*(ddplim0(imode)**3)*(rhosul*avc/mmsul)*X(imode)
       mhi (imode)=                                                     &
          (PPI/6.0)*(ddplim1(imode)**3)*(rhosul*avc/mmsul)*X(imode)
      END DO

! Initial fractions of mass in each mode among components
      mfrac_0(1,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !nucln. soluble
      mfrac_0(2,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !Aitken soluble
      mfrac_0(3,1:ncp)=(/1.0,0.0,0.0,0.0,0.0,0.0/) !accum. soluble
      mfrac_0(4,1:ncp)=(/0.0,0.0,0.0,1.0,0.0,0.0/) !coarse soluble
      mfrac_0(5,1:ncp)=(/0.0,0.5,0.5,0.0,0.0,0.0/) !Aitken insoluble
      mfrac_0(6,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !accum. insoluble
      mfrac_0(7,1:ncp)=(/0.0,0.0,0.0,0.0,1.0,0.0/) !coarse insoluble

! Modes resulting when two modes coagulate
      coag_mode(1,1:nmodes)=(/1,2,3,4,2,3,4/)
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
      coag_mode(2,1:nmodes)=(/2,2,3,4,2,3,4/)
! 2 coagulating with:    1,2,3,4,5,6,7
      coag_mode(3,1:nmodes)=(/3,3,3,4,3,3,4/)
! 3 coagulating with:    1,2,3,4,5,6,7
      coag_mode(4,1:nmodes)=(/4,4,4,4,4,4,4/)
! 4 coagulating with:    1,2,3,4,5,6,7
      coag_mode(5,1:nmodes)=(/2,2,3,4,5,6,7/)
! 5 coagulating with:    1,2,3,4,5,6,7
      coag_mode(6,1:nmodes)=(/3,3,3,4,6,6,7/)
! 6 coagulating with:    1,2,3,4,5,6,7
      coag_mode(7,1:nmodes)=(/4,4,4,4,7,7,7/)
! 7 coagulating with:    1,2,3,4,5,6,7

! Molar masses of components (kg mol-1)
      mm=(/0.098,0.012,0.0168,0.05844,0.100,0.0168/)
!          h2so4  bc     oc    nacl   dust  sec_org
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 C-H ratio)
! assume secondary organic species has mm=mm_oc
! Mass density of components (kg m^-3)
      rhocomp=(/1769.0,1500.0,1500.0,1600.0,2650.0,1500.0/)
! number of dissociating ions in soluble components
      no_ions=(/3.0,0.0,0.0,2.0,0.0,0.0/)
!
      fracbcem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
      fracocem=(/0.0,0.0,0.0,0.0,1.0,0.0,0.0/)
! fractions of primary BC/OC emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
      mode=(mode_choice > 0)
      component=.FALSE.
      soluble=.FALSE.
      DO imode=1,nmodes
        DO icp=1,ncp
          IF(((component_mode(imode,icp) == 1).AND.                     &
                (component_choice(icp) == 1)).AND.                      &
                (mode_choice(imode) == 1)) THEN
             component(imode,icp)=.true.
          END IF
          IF(soluble_choice(icp) == 1) soluble(icp)=.true.
        END DO
      END DO

      IF (lhook) CALL dr_hook('UKCA_MODE_SETUP:UKCA_MODE_DUONLY_2MODE', &
                               zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_MODE_DUonly_2MODE

      END MODULE UKCA_MODE_SETUP
