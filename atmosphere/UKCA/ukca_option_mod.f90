! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold all UKCA variables in RUN_UKCA
!          namelist
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_option_mod

USE missing_data_mod,      ONLY: rmdi, imdi
USE ukca_photo_scheme_mod, ONLY: i_ukca_nophot

IMPLICIT NONE

! Declarations for UKCA sub-model
! -----------------------------------------------------------------------------
! Namelist items

LOGICAL :: l_ukca           =.FALSE. ! True when UKCA is switched on
LOGICAL :: l_ukca_aie1      =.FALSE. ! True when 1st aerosol ind effect required
LOGICAL :: l_ukca_aie2      =.FALSE. ! True when 2nd aerosol ind effect required

! Main chemistry namelist inputs:
INTEGER :: i_ukca_chem         = 0      ! chemistry scheme to use
LOGICAL :: l_ukca_chem_aero    =.FALSE. ! add aerosol precursors to chemistry
LOGICAL :: l_ukca_trophet      =.FALSE. ! T for tropospheric heterogeneous chem
LOGICAL :: l_ukca_mode         =.FALSE. ! True for UKCA-MODE aerosol scheme
LOGICAL :: l_ukca_qch4inter    =.FALSE. ! True for interact wetland CH4 ems
LOGICAL :: l_ukca_useumuivals  =.FALSE. ! True when using UMUI CFC values
LOGICAL :: l_ukca_het_psc      =.FALSE. ! True for Het/PSC chemistry
LOGICAL :: l_ukca_sa_clim      =.FALSE. ! True to use SPARC surface area density
LOGICAL :: l_ukca_h2o_feedback =.FALSE. ! True for H2O feedback from chem
LOGICAL :: l_ukca_rado3        =.FALSE. ! T when using UKCA O3 in radiation
LOGICAL :: l_ukca_radch4       =.FALSE. ! T when using UKCA CH4 in radiation
LOGICAL :: l_ukca_radn2o       =.FALSE. ! T when using UKCA N2O in radiation
LOGICAL :: l_ukca_radf11       =.FALSE. ! T when using UKCA CFC-11 in radn
LOGICAL :: l_ukca_radf12       =.FALSE. ! T when using UKCA CFC-12 in radn
LOGICAL :: l_ukca_radf113      =.FALSE. ! T when using UKCA CFC-113 in radn
LOGICAL :: l_ukca_radf22       =.FALSE. ! T when using UKCA HCFC-22 in radn
LOGICAL :: l_ukca_radaer       =.FALSE. ! Radiative effects of UKCA aerosols
LOGICAL :: l_ukca_intdd        =.FALSE. ! T when using interact dry deposition
LOGICAL :: l_ukca_prescribech4 =.FALSE. ! T when prescribing surface ch4
LOGICAL :: l_ukca_set_trace_gases =.FALSE. ! T to use UM values for fCO2 etc
LOGICAL :: l_ukca_use_background_aerosol =.FALSE. ! use bg aerosol climatology

INTEGER :: dts0 = 300                   ! Default Backward Euler timestep
INTEGER :: nit  = 8                     ! Number of iterations of BE Solver

INTEGER :: i_ukca_photol = i_ukca_nophot ! Photolysis scheme to use

! Directory pathname for 2d photolysis rates
CHARACTER (LEN=256) :: phot2d_dir  = 'phot2d dir is unset'
! Dir pathname for 2d upper boundary data
CHARACTER (LEN=256) :: strat2d_dir = 'strat2d dir is unset'

INTEGER :: fastjx_numwl = imdi        ! No. of wavelengths to use (8, 12, 18)
INTEGER :: fastjx_mode  = imdi        ! 1 = use just 2D above prescutoff, 
                                      ! 2 = merge, 3 = just fastjx) 
REAL :: fastjx_prescutoff             ! Press for 2D stratospheric photolysis
CHARACTER(LEN=256) :: jvspec_file ='jvspec file is unset' ! FastJ spectral file
CHARACTER(LEN=256) :: jvscat_file ='jvscat file is unset' ! FastJ scatter file
CHARACTER(LEN=256) :: jvspec_dir  ='jvspec dir is unset'  ! Dir for jvspec file

! Dir for stratospheric aerosol file
CHARACTER(LEN=256) :: dir_strat_aer  = 'dir_strat_aer is unset'   
! File for stratospheric aerosol file
CHARACTER(LEN=256) :: file_strat_aer = 'file_strat_aer is unset'  

! UKCA_MODE control features:
LOGICAL :: l_ukca_primsu    =.FALSE. ! T for primary sulphate aerosol emissions
LOGICAL :: l_ukca_primss    =.FALSE. ! T for primary sea-salt aerosol emissions
LOGICAL :: l_ukca_primbcoc  =.FALSE. ! T for primary BC and OC aerosol emissions
LOGICAL :: l_ukca_primdu    =.FALSE. ! T for primary dust aerosol emissions
LOGICAL :: l_ukca_use_2dtop =.FALSE. ! T for using 2-D top boundary files
LOGICAL :: l_bcoc_ff        =.FALSE. ! T for primary fossil fuel BC/OC emissions
LOGICAL :: l_bcoc_bf        =.FALSE. ! T for primary biofuel BC/OC emissions
LOGICAL :: l_bcoc_bm        =.FALSE. ! T for primary biomass BC/OC emissions
LOGICAL :: l_mode_bhn_on    =.TRUE.  ! T for binary sulphate nucleation
LOGICAL :: l_mode_bln_on    =.TRUE.  ! T for BL sulphate nucleation
LOGICAL :: l_ukca_arg_act   =.FALSE. ! T when using AR&G aerosol activation
LOGICAL :: l_ukca_sfix      =.FALSE. ! T for diagnosing UKCA CCN at 
                                     ! fixed supersaturation

INTEGER :: i_mode_setup     = imdi        ! Defines MODE aerosol scheme
INTEGER :: i_mode_nzts      = imdi        ! No. of substeps for nucleation/
                                          ! sedimentation
INTEGER :: i_mode_bln_param_method = 1    ! 1=activ; 2=kinetc; 3=PNAS/Metzer
                                          ! 4=EUCAARI-kinetc; 5=EUCAARI-org1
                                          ! 6=EUCAARI-org2
                                          ! maps to IBLN in GLOMAP
! Not included in namelist at present:
INTEGER :: i_mode_nucscav   = 1    ! Choice of nucl. scavenging co-effs
INTEGER :: i_mode_ss_scheme = 1    ! Defines sea-salt emission scheme.


REAL :: mode_parfrac         = rmdi ! Fraction of SO2 emissions as aerosol(%)

REAL :: ukca_MeBrMMR         = rmdi ! UKCA trace gas mixing value
REAL :: ukca_MeClMMR         = rmdi ! UKCA trace gas mixing value
REAL :: ukca_CH2Br2MMR       = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H2MMR           = rmdi ! UKCA trace gas mixing value
REAL :: ukca_N2MMR           = rmdi ! UKCA trace gas mixing value
REAL :: ukca_CFC115MMR       = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_CCl4MMR         = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_MeCCl3MMR       = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_HCFC141bMMR     = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_HCFC142bMMR     = rmdi ! UKCA trace gas mixing value         
REAL :: ukca_H1211MMR        = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_H1202MMR        = rmdi ! UKCA trace gas mixing value   
REAL :: ukca_H1301MMR        = rmdi ! UKCA trace gas mixing value     
REAL :: ukca_H2402MMR        = rmdi ! UKCA trace gas mixing value   
REAL :: ukca_COSMMR          = rmdi ! UKCA trace gas mixing value

! Define the RUN_UKCA namelist

NAMELIST/RUN_UKCA/ l_ukca, l_ukca_aie1, l_ukca_aie2,              &
         i_ukca_chem, l_ukca_chem_aero,                           & 
         i_ukca_photol,                                           &
         l_ukca_mode,                                             &
         l_ukca_qch4inter,                                        &
         l_ukca_useumuivals,                                      &
         l_ukca_het_psc, l_ukca_sa_clim,                          &
         l_ukca_h2o_feedback,                                     &
         l_ukca_rado3, l_ukca_radch4, l_ukca_radn2o,              &
         l_ukca_radf11, l_ukca_radf12, l_ukca_radf113,            &
         l_ukca_radf22, l_ukca_radaer,                            &
         l_ukca_intdd, l_ukca_trophet, l_ukca_prescribech4,       &
         l_ukca_set_trace_gases, l_ukca_use_background_aerosol,   &
         l_ukca_primsu, l_ukca_primss,                            &
         l_ukca_primbcoc, l_ukca_primdu, l_ukca_use_2dtop,        &
         l_bcoc_ff, l_bcoc_bf, l_bcoc_bm, l_mode_bhn_on,          &
         l_mode_bln_on, l_ukca_arg_act,                           &
         l_ukca_sfix, i_mode_setup, i_mode_nzts,                  &
         i_mode_bln_param_method, mode_parfrac, dts0, nit,        &
         jvspec_dir, jvspec_file, jvscat_file, phot2d_dir,        &
         strat2d_dir, fastjx_numwl, fastjx_mode,                  &
         fastjx_prescutoff, dir_strat_aer, file_strat_aer,        &
         ukca_MeBrmmr, ukca_MeClmmr, ukca_CH2Br2mmr, ukca_H2mmr,  &
         ukca_N2mmr, ukca_CFC115mmr, ukca_CCl4mmr,                &
         ukca_MeCCl3mmr, ukca_HCFC141bmmr, ukca_HCFC142bmmr,      &
         ukca_H1211mmr, ukca_H1202mmr, ukca_H1301mmr,             &
         ukca_H2402mmr, ukca_COSmmr

! -----------------------------------------------------------------------------
! These are set in ukca_setup_chem_mod after the namelist is read

LOGICAL :: l_ukca_chem      =.FALSE. ! True when UKCA chemistry is on
LOGICAL :: l_ukca_ageair    =.FALSE. ! True for Age of Air
LOGICAL :: l_ukca_trop      =.FALSE. ! True for tropospheric chemistry (B-E) 
LOGICAL :: l_ukca_raq       =.FALSE. ! True for regional air quality chem (B-E)
LOGICAL :: l_ukca_tropisop  =.FALSE. ! True for trop chemistry + isoprene
LOGICAL :: l_ukca_strat     =.FALSE. ! True for strat+reduced trop chemistry
LOGICAL :: l_ukca_strattrop =.FALSE. ! True for std strat+trop chemistry
LOGICAL :: l_ukca_achem     =.FALSE. ! add aerosol chemistry to scheme (NR)
LOGICAL :: l_ukca_aerchem   =.FALSE. ! True for trop+aerosol chemistry (B-E)
LOGICAL :: l_ukca_advh2o    =.FALSE. ! True when H2O treated as tracer by ASAD

! These schemes are not yet included but logicals used in code so needed here
LOGICAL :: l_ukca_std_trop  =.FALSE. ! True for Standard Trop chemistry (N-R)
LOGICAL :: l_ukca_exttc     =.FALSE. ! True for Extended Trop. chemistry (B-E)
LOGICAL :: l_ukca_stratcfc  =.FALSE. ! True for extended strat chemistry
LOGICAL :: l_ukca_dust      =.FALSE. ! True for UKCA dust scheme

! File  and directory for reference sulphur aerosol file 
! Not currently used as l_use_stratclim in ukca_fastjx is FALSE
CHARACTER(LEN=256) :: dir_reff_sulp  = 'dir_reff_sulp is unset'   
CHARACTER(LEN=256) :: file_reff_sulp = 'file_reff_sulp is unset'  

! Tracers and chemistry integers:

INTEGER :: ukca_int_method  = 0      ! Defines chemical integration method
INTEGER :: jpctr  = 0                ! No. of chemical tracers
INTEGER :: jpspec = 0                ! No. of chemical species
INTEGER :: jpbk   = 0                ! No. of bimolecular reactions
INTEGER :: jptk   = 0                ! No. of termolecular reactions
INTEGER :: jppj   = 0                ! No. of photolytic reactions
INTEGER :: jphk   = 0                ! No. of heterogeneous reactions
INTEGER :: jpnr   = 0                ! jpbk + jptk + jppj + jphk
INTEGER :: jpdd   = 0                ! No. of dry deposited species
INTEGER :: jpdw   = 0                ! No. of wet deposited species

CONTAINS

SUBROUTINE check_run_ukca()

! Description:
!   Subroutine to apply logic checks based on the 
!   options selected in the run_ukca namelist.
!   Note that some chemistry scheme specific checks
!   are done in ukca_setup_chem

USE ukca_photo_scheme_mod, ONLY:  i_ukca_nophot, i_ukca_phot2d, &
                            i_ukca_fastj,  i_ukca_fastjx  
USE ereport_mod,     ONLY: ereport
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER   :: routinename = 'check_run_ukca'
CHARACTER (LEN=70)             :: cmessage   ! Error message
INTEGER                        :: errcode    ! Variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)

! Check photolysis switches
IF ( i_ukca_photol /= i_ukca_nophot .AND. &
     i_ukca_photol /= i_ukca_phot2d .AND. &
     i_ukca_photol /= i_ukca_fastj  .AND. &
     i_ukca_photol /= i_ukca_fastjx ) THEN
   cmessage='Unknown photolysis scheme.'
   errcode=ABS(i_ukca_photol)
   CALL ereport(routinename,errcode,cmessage)
END IF

! GLOMAP-mode related switches
IF ( l_ukca_mode .AND. .NOT. l_ukca ) THEN
   cmessage='Cannot use GLOMAP-mode aerosols without UKCA'
   errcode=1
   CALL ereport(routinename,errcode,cmessage)
END IF 

! Direct aerosol effects
IF ( l_ukca_radaer .AND. .NOT. l_ukca_mode ) THEN
   cmessage='Cannot use RADAER without GLOMAP-mode aerosols'
   errcode=2
   CALL ereport(routinename,errcode,cmessage)
END IF 

! Indirect aerosol effects
IF ( (l_ukca_aie1 .OR. l_ukca_aie2) .AND. .NOT. l_ukca_mode ) THEN
   cmessage='Cannot use AIE without GLOMAP-mode aerosols'
   errcode=2
   CALL ereport(routinename,errcode,cmessage)
END IF 

IF (lhook) CALL dr_hook(routinename,zhook_out,zhook_handle)

END SUBROUTINE check_run_ukca
END MODULE ukca_option_mod
