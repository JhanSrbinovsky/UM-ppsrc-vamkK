! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing standard tropospheric
!    chemistry for the ASAD system. Will configure tropospheric chemistry
!    to include aerosol chemistry and isoprene chemistry as required.
!    The subroutine INIT_STD_TROP constructs the final arrays from the
!    components using logicals defined in the UMUI.
!
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

MODULE UKCA_CHEM_STD_TROP

USE ASAD_MOD,           ONLY: jddepc, jddept
USE UKCA_CHEM_DEFS_MOD, ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t
USE ukca_option_mod,    ONLY: jpspec, jpbk, jptk, jphk, jppj, jpdd, jpdw 
USE Control_Max_Sizes
IMPLICIT NONE
PRIVATE

! Standard tropospheric chemistry with optional: aerosol and isoprene chemistry
! =============================================================================

TYPE(CHCH_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: chch_defs_std_trop
TYPE(RATB_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratb_defs_std_trop
TYPE(RATT_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratt_defs_std_trop
TYPE(RATJ_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratj_defs_std_trop
TYPE(RATH_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: rath_defs_std_trop
REAL,     DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: depvel_defs_std_trop
REAL,     DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: henry_defs_std_trop

! No of heterogenous processes
INTEGER, PARAMETER, PUBLIC :: nhet_std_trop = 2        ! On aerosol surfaces
INTEGER, PARAMETER, PUBLIC :: naq_std_trop  = 2        ! Aqueous phase

! No of dry deposited species
INTEGER, PARAMETER, PUBLIC :: ndry_std_trop = 22       ! Std. Tropospheric

! No of wet deposited species
INTEGER, PARAMETER, PUBLIC :: nwet_std_trop = 15       ! Std. Tropospheric

TYPE(CHCH_T), DIMENSION( 46), PUBLIC :: chch_defs_std_trop1=(/          &
chch_t(  1,'O3        ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t(  1,'O(1D)     ',  1,'SS        ','          ',  0,  0,  0),     &
chch_t(  2,'NO        ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t(  3,'NO3       ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t(  4,'NO2       ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t(  5,'N2O5      ',  2,'TR        ','          ',  1,  1,  0),     &
chch_t(  6,'HO2NO2    ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t(  7,'HONO2     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t(  7,'H2        ',  1,'CT        ','          ',  0,  0,  0),     &
chch_t(  8,'H2O2      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t(  9,'CH4       ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t( 10,'CO        ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t( 10,'CO2       ',  1,'CT        ','          ',  0,  0,  0),     &
chch_t( 11,'HCHO      ',  1,'TR        ','          ',  1,  1,  1),     &
chch_t( 11,'H2O       ',  1,'CF        ','          ',  0,  0,  0),     &
chch_t( 12,'MeOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 13,'HONO      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 13,'O2        ',  1,'CT        ','          ',  0,  0,  0),     &
chch_t( 13,'N2        ',  1,'CT        ','          ',  0,  0,  0),     &
chch_t( 14,'C2H6      ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 15,'EtOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 16,'MeCHO     ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t( 17,'PAN       ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 18,'C3H8      ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 19,'n-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 20,'i-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 21,'EtCHO     ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 22,'Me2CO     ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 23,'MeCOCH2OOH',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 24,'PPAN      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 25,'MeONO2    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 26,'O3S       ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 26,'O(3P)     ',  1,'SS        ','          ',  0,  0,  0),     &
chch_t( 27,'OH        ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 28,'HO2       ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 29,'MeOO      ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 30,'EtOO      ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 31,'MeCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 32,'n-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 33,'i-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 34,'EtCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 35,'MeCOCH2OO ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 36,'O(3P)S    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 37,'O(1D)S    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 38,'OHS       ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 39,'HO2S      ',  1,'TR        ','          ',  0,  1,  0)      &
 /)


TYPE(RATB_T), DIMENSION(  30),PARAMETER :: ratb_defs_a=(/                &
ratb_t('HO2       ','NO        ','OH        ','NO2       ','          ', &
'          ',3.60E-12,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','NO3       ','OH        ','NO2       ','          ', &
'          ',4.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','O3        ','OH        ','O2        ','          ', &
'          ',2.03E-16,  4.57,   -693.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','HO2       ','H2O2      ','          ','          ', &
'          ',2.20E-13,  0.00,   -600.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeOO      ','MeOOH     ','          ','          ', &
'          ',3.80E-13,  0.00,   -780.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeOO      ','HCHO      ','          ','          ', &
'          ',3.80E-13,  0.00,   -780.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','EtOO      ','EtOOH     ','          ','          ', &
'          ',3.80E-13,  0.00,   -900.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeCO3     ','MeCO3H    ','          ','          ', &
'          ',2.08E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ', &
'          ',1.04E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeCO3     ','OH        ','MeOO      ','          ', &
'          ',2.08E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','n-PrOO    ','n-PrOOH   ','          ','          ', &
'          ',1.51E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','i-PrOO    ','i-PrOOH   ','          ','          ', &
'          ',1.51E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ', &
'          ',3.05E-13,  0.00,  -1040.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ', &
'          ',1.25E-13,  0.00,  -1040.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ', &
'          ',1.36E-13,  0.00,  -1250.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ', &
'          ',2.95E-12,  0.00,   -285.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','NO        ','MeONO2    ','          ','          ', &
'          ',2.95E-15,  0.00,   -285.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ', &
'          ',1.30E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','MeOO      ','MeOH      ','HCHO      ','          ', &
'          ',1.03E-13,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','MeOO      ','HO2       ','HO2       ','HCHO      ', &
'HCHO      ',1.03E-13,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','MeCO3     ','HO2       ','HCHO      ','MeOO      ', &
'          ',1.80E-12,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeOO      ','MeCO3     ','MeCO2H    ','HCHO      ','          ', &
'          ',2.00E-13,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ', &
'          ',2.60E-12,  0.00,   -380.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ', &
'          ',2.30E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('EtOO      ','MeCO3     ','MeCHO     ','HO2       ','MeOO      ', &
'          ',4.40E-13,  0.00,  -1070.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ', &
'          ',7.50E-12,  0.00,   -290.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeCO3     ','NO3       ','MeOO      ','CO2       ','NO2       ', &
'          ',4.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('n-PrOO    ','NO        ','EtCHO     ','HO2       ','NO2       ', &
'          ',2.90E-12,  0.00,   -350.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ', &
'          ',2.50E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('i-PrOO    ','NO        ','Me2CO     ','HO2       ','NO2       ', &
'          ',2.70E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000)     &
 /)

TYPE(RATB_T), DIMENSION(  30),PARAMETER :: ratb_defs_b=(/                &
ratb_t('i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ', &
'          ',2.50E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ', &
'          ',6.70E-12,  0.00,   -340.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('EtCO3     ','NO3       ','EtOO      ','CO2       ','NO2       ', &
'          ',4.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ', &
'          ',2.80E-12,  0.00,   -300.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ', &
'          ',2.50E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO        ','NO3       ','NO2       ','NO2       ','          ', &
'          ',1.80E-11,  0.00,   -110.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO        ','O3        ','NO2       ','          ','          ', &
'          ',1.40E-12,  0.00,   1310.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO2       ','O3        ','NO3       ','          ','          ', &
'          ',1.40E-13,  0.00,   2470.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ', &
'          ',2.00E-12,  0.00,   2440.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ', &  ! 40
'          ',1.40E-12,  0.00,   1860.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ', &
'          ',3.46E-12,  0.00,   1862.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('NO3       ','Me2CO     ','HONO2     ','MeCOCH2OO ','          ', &
'          ',3.00E-17,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('N2O5      ','H2O       ','HONO2     ','HONO2     ','          ', &
'          ',2.50E-22,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(3P)     ','O3        ','O2        ','O2        ','          ', &
'          ',8.00E-12,  0.00,   2060.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','CH4       ','OH        ','MeOO      ','          ', &
'          ',1.05E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','CH4       ','HCHO      ','H2        ','          ', &
'          ',7.50E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','CH4       ','HCHO      ','HO2       ','HO2       ', &
'          ',3.45E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','H2O       ','OH        ','OH        ','          ', &
'          ',2.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','N2        ','O(3P)     ','N2        ','          ', &
'          ',2.10E-11,  0.00,   -115.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)     ','O2        ','O(3P)     ','O2        ','          ', &  ! 50
'          ',3.20E-11,  0.00,    -67.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','CH4       ','H2O       ','MeOO      ','          ', &
'          ',1.85E-12,  0.00,   1690.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','C2H6      ','H2O       ','EtOO      ','          ', &
'          ',6.90E-12,  0.00,   1000.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','C3H8      ','n-PrOO    ','H2O       ','          ', &
'          ',7.60E-12,  0.00,    585.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','C3H8      ','i-PrOO    ','H2O       ','          ', &
'          ',7.60E-12,  0.00,    585.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','CO        ','HO2       ','          ','          ', &
'          ',1.44E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','EtCHO     ','H2O       ','EtCO3     ','          ', &
'          ',5.10E-12,  0.00,   -405.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','EtOOH     ','H2O       ','MeCHO     ','OH        ', &
'          ',8.01E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','EtOOH     ','H2O       ','EtOO      ','          ', &
'          ',1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','H2        ','H2O       ','HO2       ','          ', &
'          ',7.70E-12,  0.00,   2100.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','H2O2      ','H2O       ','HO2       ','          ', &   ! 60
'          ',2.90E-12,  0.00,    160.00, 0.000, 0.000, 0.000, 0.000)     &
 /)

TYPE(RATB_T), DIMENSION(  23),PARAMETER :: ratb_defs_c=(/                &
ratb_t('OH        ','HCHO      ','H2O       ','HO2       ','CO        ', &
'          ',5.40E-12,  0.00,   -135.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','HO2       ','H2O       ','          ','          ', &
'          ',4.80E-11,  0.00,   -250.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','HO2NO2    ','H2O       ','NO2       ','          ', &
'          ',1.90E-12,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','HONO2     ','H2O       ','NO3       ','          ', &
'          ',1.50E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','HONO      ','H2O       ','NO2       ','          ', &
'          ',2.50E-12,  0.00,   -260.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ', &
'          ',1.02E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeOOH     ','H2O       ','MeOO      ','          ', &
'          ',1.89E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ', &
'          ',4.00E-13,  0.00,    845.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ', &
'          ',8.80E-12,  0.00,   1320.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ', &   ! 70
'          ',1.70E-14,  0.00,   -420.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeCOCH2OOH','H2O       ','MeCOCH2OO ','          ', &
'          ',1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeCOCH2OOH','OH        ','MGLY      ','          ', &
'          ',8.39E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','MeCHO     ','H2O       ','MeCO3     ','          ', &
'          ',4.40E-12,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','NO3       ','HO2       ','NO2       ','          ', &
'          ',2.00E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','O3        ','HO2       ','O2        ','          ', &
'          ',1.70E-12,  0.00,    940.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','OH        ','H2O       ','O(3P)     ','          ', &
'          ',6.31E-14,  2.60,   -945.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','PAN       ','HCHO      ','NO2       ','H2O       ', &
'          ',3.00E-14,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','PPAN      ','MeCHO     ','NO2       ','H2O       ', &
'          ',1.27E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','n-PrOOH   ','n-PrOO    ','H2O       ','          ', &
'          ',1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','n-PrOOH   ','EtCHO     ','H2O       ','OH        ', &   ! 80
'          ',1.10E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','i-PrOOH   ','i-PrOO    ','H2O       ','          ', &
'          ',1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ', &
'          ',1.66E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(3P)     ','NO2       ','NO        ','O2        ','          ', &
'          ',5.50E-12,  0.00,   -188.00, 0.000, 0.000, 0.000, 0.000)     &
 /)


TYPE(RATB_T), DIMENSION(  83    ),PARAMETER :: ratb_defs_std_trop1=      &
  (/ ratb_defs_a, ratb_defs_b, ratb_defs_c /)


TYPE(RATH_T), DIMENSION(   0) :: rath_defs_std_trop1


TYPE(RATJ_T), DIMENSION(  25) :: ratj_defs_std_trop1=(/                 &
ratj_t('EtOOH     ','PHOTON    ','MeCHO     ','HO2       ','OH        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,       &
ratj_t('H2O2      ','PHOTON    ','OH        ','OH        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jh2o2     ') ,       &
ratj_t('HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhchoa    ') ,       &
ratj_t('HCHO      ','PHOTON    ','H2        ','CO        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhchob    ') ,       &
ratj_t('HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jpna      ') ,       &
ratj_t('HONO2     ','PHOTON    ','OH        ','NO2       ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhono2    ') ,       &
ratj_t('MeCHO     ','PHOTON    ','MeOO      ','HO2       ','CO        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jaceta    ') ,       &
ratj_t('MeCHO     ','PHOTON    ','CH4       ','CO        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jacetb    ') ,       &
ratj_t('N2O5      ','PHOTON    ','NO3       ','NO2       ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jn2o5     ') ,       &
ratj_t('NO2       ','PHOTON    ','NO        ','O(3P)     ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno2      ') ,       &
ratj_t('NO3       ','PHOTON    ','NO        ','O2        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno3a     ') ,       &
ratj_t('NO3       ','PHOTON    ','NO2       ','O(3P)     ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno3b     ') ,       &
ratj_t('O2        ','PHOTON    ','O(3P)     ','O(3P)     ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo2       ') ,       &
ratj_t('O3        ','PHOTON    ','O2        ','O(1D)     ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3a      ') ,       &
ratj_t('O3        ','PHOTON    ','O2        ','O(3P)     ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3b      ') ,       &
ratj_t('PAN       ','PHOTON    ','MeCO3     ','NO2       ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jpan      ') ,       &
ratj_t('PPAN      ','PHOTON    ','EtCO3     ','NO2       ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jpan      ') ,       &
ratj_t('HONO      ','PHOTON    ','OH        ','NO        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhono     ') ,       &
ratj_t('EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jetcho    ') ,       &
ratj_t('Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jaceto    ') ,       &
ratj_t('MeONO2    ','PHOTON    ','HO2       ','HCHO      ','NO2       ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmena     ') ,       &
ratj_t('MeOOH     ','PHOTON    ','HO2       ','HCHO      ','OH        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,       &
ratj_t('n-PrOOH   ','PHOTON    ','HO2       ','EtCHO     ','OH        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,       &
ratj_t('i-PrOOH   ','PHOTON    ','HO2       ','Me2CO     ','OH        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,       &
ratj_t('MeCOCH2OOH','PHOTON    ','HCHO      ','MeCO3     ','OH        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ')         &
  /)
!ratj_t('O3S       ','PHOTON    ','O2        ','O(1D)     ','          ',&
!  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3a      ') ,       &
!ratj_t('O3S       ','PHOTON    ','O2        ','O(3P)     ','          ',&
!  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3b      ')         &
!  /)


TYPE(RATT_T), DIMENSION(  13) :: ratt_defs_std_trop1=(/                 &
ratt_t('HO2       ','HO2       ','H2O2      ','O2        ',             &
  0.0, 1.90E-33, 0.00, -980.0, 0.00E+00,  0.00,  0.0, 0.000, 0.000),    &
ratt_t('HO2       ','NO2       ','HO2NO2    ','m         ',             &
  0.6, 1.80E-31, -3.20,   0.0, 4.70E-12,  0.00,  0.0, 0.000, 0.000) ,   &
ratt_t('HO2NO2    ','m         ','HO2       ','NO2       ',             &
  0.6, 4.10E-05, 0.00, 10650.0, 4.80E+15, 0.00, 11170.0, 0.000, 0.000), &
ratt_t('MeCO3     ','NO2       ','PAN       ','m         ',             &
  0.3, 2.70E-28, -7.10,   0.0, 1.20E-11, -0.90,  0.0, 0.000, 0.000) ,   &
ratt_t('PAN       ','m         ','MeCO3     ','NO2       ',             &
  0.3, 4.90E-03, 0.00, 12100.0, 5.40E+16, 0.00, 13830.0, 0.000, 0.000), &
ratt_t('N2O5      ','m         ','NO2       ','NO3       ',             &
  0.3, 1.30E-03, -3.50, 11000.0, 9.70E+14, 0.10, 11080.0, 0.000, 0.000),&
ratt_t('NO2       ','NO3       ','N2O5      ','m         ',             &
  0.3, 3.60E-30, -4.10,     0.0, 1.90E-12,  0.20, 0.0, 0.000, 0.000) ,  &
ratt_t('O(3P)     ','O2        ','O3        ','m         ',             &
  0.0, 5.70E-34, -2.60,     0.0, 0.00E+00,  0.00, 0.0, 0.000, 0.000) ,  &
ratt_t('OH        ','NO        ','HONO      ','m         ',             &
1420.0,7.40E-31, -2.40,  0.0, 3.30E-11, -0.30, 0.0, 0.000, 0.000) ,     &
ratt_t('OH        ','NO2       ','HONO2     ','m         ',             &
  0.4, 3.30E-30, -3.00,     0.0, 4.10E-11,  0.00, 0.0, 0.000, 0.000) ,  &
ratt_t('OH        ','OH        ','H2O2      ','m         ',             &
  0.5, 6.90E-31, -0.80,     0.0, 2.60E-11,  0.00, 0.0, 0.000, 0.000) ,  &
ratt_t('EtCO3     ','NO2       ','PPAN      ','m         ',             &
  0.3, 2.70E-28, -7.10,     0.0, 1.20E-11, -0.90, 0.0, 0.000, 0.000) ,  &
ratt_t('PPAN      ','m         ','EtCO3     ','NO2       ',             &
  0.4, 1.70E-03,  0.00, 11280.0, 8.30E+16, 0.00, 13940.0, 0.000, 0.000) &
  /)



! Dry deposition array:
! Dry deposition rates are written in a table format, with the columns corresponding
! to times and the rows corresponding to surface type:
!
!           SUMMER                       WINTER
!       day  night  24hr    day  night  24hr
!         x     x     x       x     x     x     Water
!         x     x     x       x     x     x     Forest
!         x     x     x       x     x     x     Grass/Shrub
!         x     x     x       x     x     x     Desert
!         x     x     x       x     x     x     Snow/Ice
!
! '24 hr' corresponds to the 24 hour average deposition velocity, currently not used.
! The units of x are cm s-1

! These have to be in the same order as the species which dry deposit in the
! chch array.

REAL, DIMENSION(6,5,ndry_std_trop) :: depvel_defs_std_trop1=RESHAPE((/   &
        0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  & ! O3
        0.85,  0.30,  0.65,  0.65,  0.25,  0.45,  & 
        0.65,  0.25,  0.45,  0.65,  0.25,  0.45,  & 
        0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  & 
        0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & ! NO
        0.14,  0.01,  0.07,  0.01,  0.01,  0.01,  &
        0.10,  0.01,  0.06,  0.01,  0.01,  0.01,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! NO3
        0.83,  0.04,  0.44,  0.06,  0.04,  0.05,  &
        0.63,  0.06,  0.35,  0.07,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! NO2
        0.83,  0.04,  0.44,  0.06,  0.04,  0.05,  &
        0.63,  0.06,  0.35,  0.07,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! N2O5
        4.00,  3.00,  3.50,  2.00,  1.00,  1.50,  &
        2.50,  1.50,  2.00,  1.00,  1.00,  1.00,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! HO2NO2
        3.20,  1.40,  2.30,  2.00,  1.00,  1.50,  & 
        1.80,  0.80,  1.30,  1.00,  1.00,  1.00,  & 
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! HONO2
        3.20,  1.40,  2.30,  2.00,  1.00,  1.50,  & 
        1.80,  0.80,  1.30,  1.00,  1.00,  1.00,  & 
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! H2O2
        1.25,  0.16,  0.71,  0.28,  0.12,  0.20,  &
        1.25,  0.53,  0.89,  0.83,  0.78,  0.81,  &
        0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  &
        0.32,  0.32,  0.32,  0.32,  0.32,  0.32,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & ! CH4
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & ! CO
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        1.00,  0.60,  0.80,  1.00,  0.60,  0.80,  & ! HCHO
        1.60,  0.40,  1.00,  0.01,  0.01,  0.01,  & 
        0.71,  0.20,  0.45,  0.03,  0.03,  0.03,  & 
        0.40,  0.40,  0.40,  0.00,  0.00,  0.00,  & 
        0.30,  0.30,  0.30,  0.30,  0.30,  0.30,  & 
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! MeOOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  &
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & ! HONO
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! EtOOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  &
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! MeCHO
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & ! PAN
        0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & 
        0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & 
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & 
        0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  &
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! n-PrOOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  &
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! i-PrOOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  &
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! EtCHO
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! MeCOCH2OOH
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & ! PPAN
        0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & 
        0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & 
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & 
        0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  &
        0.07,  0.07,  0.07,  0.07,  0.07,  0.07,  & ! O3S
        1.00,  0.11,  0.56,  0.26,  0.11,  0.19,  &
        1.00,  0.37,  0.69,  0.59,  0.46,  0.53,  &
        0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  &
        0.07,  0.07,  0.07,  0.07,  0.07,  0.07   &
        /),(/  6,  5, ndry_std_trop/) )


! The following formula is used to calculate the effective Henry's Law coef,
! which takes the affects of dissociation and complex formation on a species'
! solubility (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 and 2 above give the data for this gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase the above term is multiplied by
! another factor of 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 and 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
! The data in columns 5 and 6 give the data for a second dissociation,
! e.g for SO2, HSO3^{-}, and SO3^{2-}
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]

REAL, DIMENSION(6,nwet_std_trop) :: henry_defs_std_trop1=RESHAPE((/     &
 0.2000E+01, 0.2000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! NO3
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! N2O5
 0.1300E+05, 0.6900E+04, 0.1000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HO2NO2
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HONO2
 0.8300E+05, 0.7400E+04, 0.2400E-11,-0.3730E+04, 0.0000E+00, 0.0000E+00,&  ! H2O2
 0.3300E+04, 0.6500E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HCHO
 0.3100E+03, 0.5000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeOOH
 0.5000E+02, 0.4900E+04, 0.5600E-03,-0.1260E+04, 0.0000E+00, 0.0000E+00,&  ! HONO
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! EtOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! n-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! i-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeCOCH2OOH
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HO2
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeOO
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00 &  ! HO2S
  /),(/  6, nwet_std_trop/) )


INTERFACE ukca_init_std_trop
  MODULE PROCEDURE ukca_init_std_trop
END INTERFACE ukca_init_std_trop
PUBLIC ukca_init_std_trop

CONTAINS

! ######################################################################

SUBROUTINE UKCA_INIT_STD_TROP

! To initialise chemistry files from specified species and rate arrays
! using the defined logicals. Composes chemistry definition arrays from
! components using logical definitions from RUN_UKCA namelist held in 
! UKCA module ukca_option_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

CHARACTER (LEN=72) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_CHEM_STD_TROP:UKCA_INIT_STD_TROP',zhook_in,   &
                        zhook_handle)

ALLOCATE(chch_defs_std_trop(jpspec))
ALLOCATE(ratb_defs_std_trop(jpbk))
ALLOCATE(ratj_defs_std_trop(jppj))
ALLOCATE(ratt_defs_std_trop(jptk))
ALLOCATE(rath_defs_std_trop(jphk))
ALLOCATE(depvel_defs_std_trop(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_std_trop(6,jpdw))


! Standard tropospheric chemistry
chch_defs_std_trop = chch_defs_std_trop1
ratb_defs_std_trop = (/ ratb_defs_a, ratb_defs_b, ratb_defs_c /)
ratj_defs_std_trop = ratj_defs_std_trop1
ratt_defs_std_trop = ratt_defs_std_trop1
rath_defs_std_trop = rath_defs_std_trop1
depvel_defs_std_trop = depvel_defs_std_trop1
henry_defs_std_trop = henry_defs_std_trop1


IF (lhook) CALL dr_hook('UKCA_CHEM_STD_TROP:UKCA_INIT_STD_TROP',zhook_out,  &
                        zhook_handle)
RETURN
END SUBROUTINE UKCA_INIT_STD_TROP

END MODULE
