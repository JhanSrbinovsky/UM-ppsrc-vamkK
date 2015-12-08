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

MODULE UKCA_CHEM_TROPISOP

USE ASAD_MOD,           ONLY: jddepc, jddept
USE UKCA_CHEM_DEFS_MOD, ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t
USE ukca_option_mod,    ONLY: l_ukca_tropisop, l_ukca_achem, jpspec,  &
                              jpbk, jptk, jphk, jppj, jpdd, jpdw
USE Control_Max_Sizes
IMPLICIT NONE
PRIVATE

! Standard tropospheric chemistry with isoprene and optional: aerosol chemistry
! =============================================================================

TYPE(CHCH_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: chch_defs_tropisop
TYPE(RATB_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratb_defs_tropisop
TYPE(RATT_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratt_defs_tropisop
TYPE(RATJ_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: ratj_defs_tropisop
TYPE(RATH_T), DIMENSION(:), ALLOCATABLE, PUBLIC :: rath_defs_tropisop
REAL,     DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: depvel_defs_tropisop
REAL,     DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: henry_defs_tropisop

! No of heterogenous processes
INTEGER, PARAMETER, PUBLIC :: nhet_tropisop = 2        ! On aerosol surfaces
INTEGER, PARAMETER, PUBLIC :: naq_tropisop  = 3        ! Aqueous phase

! No of dry deposited species
INTEGER, PARAMETER, PUBLIC :: ndry_isop     = 33       ! + C5H8 chemistry
INTEGER, PARAMETER, PUBLIC :: ndry_aer      = 40       ! + Aerosol chemistry

! No of wet deposited species
INTEGER, PARAMETER, PUBLIC :: nwet_isop     = 24       ! + C5H8 chemistry
INTEGER, PARAMETER, PUBLIC :: nwet_aer      = 28       ! + Aerosol chemistry


TYPE(CHCH_T), DIMENSION(56), PUBLIC :: chch_defs_isop=(/                &
chch_t(  1,'O3        ',  1,'TR        ','          ',  1,  1,  0),     &
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
chch_t( 25,'O(3P)     ',  1,'SS        ','          ',  0,  0,  0),     &
chch_t( 26,'C5H8      ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 27,'ISOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 28,'ISON      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 29,'MACR      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 30,'MACROOH   ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 31,'MPAN      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 32,'HACET     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 33,'MGLY      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 34,'NALD      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 35,'HCOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 36,'MeCO3H    ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 37,'MeCO2H    ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 38,'ISO2      ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 39,'MACRO2    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 40,'OH        ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 41,'HO2       ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 42,'MeOO      ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 43,'EtOO      ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 44,'MeCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 45,'n-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 46,'i-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 47,'EtCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 48,'MeCOCH2OO ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 49,'MeOH      ',  1,'TR        ','          ',  1,  1,  0)      &
 /)

TYPE(CHCH_T), DIMENSION( 67), PUBLIC :: chch_defs_aer=(/                &
chch_t(  1,'O3        ',  1,'TR        ','          ',  1,  1,  0),     &
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
chch_t( 25,'O(3P)     ',  1,'SS        ','          ',  0,  0,  0),     &
chch_t( 26,'C5H8      ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 27,'ISOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 28,'ISON      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 29,'MACR      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 30,'MACROOH   ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 31,'MPAN      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 32,'HACET     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 33,'MGLY      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 34,'NALD      ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 35,'HCOOH     ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 36,'MeCO3H    ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 37,'MeCO2H    ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 38,'ISO2      ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 39,'MACRO2    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 40,'DMS       ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 41,'SO2       ',  1,'TR        ','          ',  1,  1,  1),     &
chch_t( 42,'H2SO4     ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 43,'MSA       ',  1,'TR        ','          ',  1,  0,  0),     &
chch_t( 44,'DMSO      ',  1,'TR        ','          ',  1,  1,  0),     &
chch_t( 45,'NH3       ',  1,'TR        ','          ',  1,  1,  1),     &
chch_t( 46,'CS2       ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 47,'COS       ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 48,'H2S       ',  1,'TR        ','          ',  0,  0,  1),     &
chch_t( 49,'OH        ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 50,'HO2       ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 51,'MeOO      ',  1,'TR        ','          ',  0,  1,  0),     &
chch_t( 52,'EtOO      ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 53,'MeCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 54,'n-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 55,'i-PrOO    ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 56,'EtCO3     ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 57,'MeCOCH2OO ',  1,'TR        ','          ',  0,  0,  0),     &
chch_t( 58,'MeOH      ',  1,'TR        ','          ',  1,  1,  1),     &
chch_t( 59,'Monoterp  ',  1,'TR        ','          ',  1,  0,  1),     &
chch_t( 60,'Sec_Org   ',  1,'TR        ','          ',  1,  1,  0)      &
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


TYPE(RATB_T), DIMENSION(34),PARAMETER :: ratb_defs_isop=(/               &
ratb_t('OH        ','C5H8      ','ISO2      ','          ','          ', & ! IUPAC09
'          ',2.70E-11,  0.00,   -390.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('O3        ','C5H8      ','MACR      ','HCHO      ','MACRO2    ', &
'MeCO3     ',3.33E-15,  0.00,   1995.00, 1.950, 1.740, 0.300, 0.300),    & ! Split 1
ratb_t('O3        ','C5H8      ','MeOO      ','HCOOH     ','CO        ', & ! IUPAC09
'H2O2      ',3.33E-15,  0.00,   1995.00, 0.240, 0.840, 0.420, 0.270),    & ! Split 2
ratb_t('O3        ','C5H8      ','HO2       ','OH        ','          ', & ! IUPAC09
'          ',3.33E-15,  0.00,   1995.00, 0.750, 0.750, 0.000, 0.000),    & ! Split 3
ratb_t('NO3       ','C5H8      ','ISON      ','          ','          ', & ! IUPAC09
'          ',3.15E-12,  0.00,    450.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('NO        ','ISO2      ','NO2       ','MACR      ','HCHO      ', & ! Split 1
'HO2       ',2.43E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('NO        ','ISO2      ','ISON      ','          ','          ', & ! Split 2
'          ',1.12E-13,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('HO2       ','ISO2      ','ISOOH     ','          ','          ', &
'          ',2.05E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('ISO2      ','ISO2      ','MACR      ','MACR      ','HCHO      ', &
'HO2       ',2.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','ISOOH     ','MACR      ','OH        ','          ', &
'          ',1.00E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','ISON      ','HACET     ','NALD      ','          ', &
'          ',1.30E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MACR      ','MACRO2    ','          ','          ', &
'          ',1.30E-12,  0.00,   -610.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MACR      ','MACRO2    ','          ','          ', &
'          ',4.00E-12,  0.00,   -380.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ', & ! Split 1
'CO        ',2.13E-16,  0.00,   1520.00, 1.800, 0.900, 0.640, 0.440),    &
ratb_t('O3        ','MACR      ','OH        ','MeCO3     ','          ', & ! Split 2
'          ',2.13E-16,  0.00,   1520.00, 0.380, 0.200, 0.000, 0.000),    &
ratb_t('O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ', & ! Split 1
'CO        ',3.50E-16,  0.00,   2100.00, 1.800, 0.900, 0.640, 0.440),    &
ratb_t('O3        ','MACR      ','OH        ','MeCO3     ','          ', & ! Split 2
'          ',3.50E-16,  0.00,   2100.00, 0.380, 0.200, 0.000, 0.000),    &
ratb_t('NO        ','MACRO2    ','NO2       ','MeCO3     ','HACET     ', & ! Split 1
'CO        ',1.27E-12,  0.00,   -360.00, 2.000, 0.500, 0.500, 0.500),    &
ratb_t('NO        ','MACRO2    ','MGLY      ','HCHO      ','HO2       ', & ! Split 2
'          ',1.27E-12,  0.00,   -360.00, 1.000, 1.500, 1.500, 0.000),    &
ratb_t('HO2       ','MACRO2    ','MACROOH   ','          ','          ', &
'          ',1.82E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('MACRO2    ','MeOO      ','MGLY      ','HACET     ','MeCO3     ', & ! Split 1
'HCHO      ',1.00E-12,  0.00,      0.00, 1.000, 0.750, 0.250, 2.750),    &
ratb_t('MACRO2    ','MeOO      ','HO2       ','CO        ','          ', & ! Split 2
'          ',1.00E-12,  0.00,      0.00, 1.170, 0.250, 0.000, 0.000),    &
ratb_t('MACRO2    ','MACRO2    ','HACET     ','MGLY      ','HCHO      ', & ! Split 1
'CO        ',1.00E-12,  0.00,      0.00, 2.000, 2.000, 1.000, 1.000),    &
ratb_t('MACRO2    ','MACRO2    ','HO2       ','          ','          ', & ! Split 2
'          ',1.00E-12,  0.00,      0.00, 2.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MACROOH   ','MACRO2    ','          ','          ', &
'          ',3.00E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MPAN      ','HACET     ','NO2       ','          ', &
'          ',2.90E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','HACET     ','MGLY      ','HO2       ','          ', &
'          ',3.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MGLY      ','MeCO3     ','CO        ','          ', &
'          ',1.50E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ', &
'          ',3.46E-12,  0.00,   1860.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','NALD      ','HCHO      ','CO        ','NO2       ', &
'          ',4.40E-12,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MeCO3H    ','MeCO3     ','          ','          ', &
'          ',3.70E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MeCO2H    ','MeOO      ','          ','          ', &
'          ',4.00E-13,  0.00,   -200.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','HCOOH     ','HO2       ','          ','          ', &
'          ',4.50E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('OH        ','MeOH      ','HO2       ','HCHO      ','          ', & !  117
'          ',  2.85E-12,  0.00,    345.00, 0.000, 0.000, 0.000, 0.000) &   !  117  IUPAC2007
 /)

TYPE(RATB_T), DIMENSION(  10),PARAMETER :: ratb_defs_aer=(/              &
ratb_t('DMS       ','OH        ','SO2       ','MeOO      ','HCHO      ', &
'          ',9.60E-12,  0.00,    240.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ', &
'          ',3.04E-12,  0.00,   -350.00, 0.600, 0.400, 0.000, 0.000),    &
ratb_t('DMS       ','NO3       ','SO2       ','HONO2     ','MeOO      ', &
'HCHO      ',1.90E-13,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('DMSO      ','OH        ','SO2       ','MSA       ','          ', &
'          ',5.80E-11,  0.00,      0.00, 0.600, 0.400, 0.000, 0.000),    &
ratb_t('CS2       ','OH        ','SO2       ','COS       ','          ', &
'          ',8.80E-16,  0.00,  -2300.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('H2S       ','OH        ','SO2       ','          ','          ', &
'          ',6.00E-12,  0.00,     75.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('COS       ','OH        ','SO2       ','          ','          ', &
'          ',1.10E-13,  0.00,   1200.00, 0.000, 0.000, 0.000, 0.000),    &
ratb_t('Monoterp  ','OH        ','Sec_Org   ','          ','          ', &
'          ',1.20E-11,  0.00,   -444.00, 0.130, 0.000, 0.000, 0.000),    &
ratb_t('Monoterp  ','O3        ','Sec_Org   ','          ','          ', &
'          ',1.01E-15,  0.00,   732.00, 0.130, 0.000, 0.000, 0.000),     &
ratb_t('Monoterp  ','NO3       ','Sec_Org   ','          ','          ', &
'          ',1.19E-12,  0.00,   -925.00, 0.130, 0.000, 0.000, 0.000)     &
 /)

TYPE(RATH_T), DIMENSION(   0) :: rath_defs_std_trop1

TYPE(RATH_T), DIMENSION(   0) :: rath_defs_isop

! The rate coefficients for these reactions are
! defined in the asad_hetero routine. The 'NULL' products serve to identify the
! reactions - these reactions have no gas-phase products.
TYPE(RATH_T), DIMENSION(   5) :: rath_defs_aer=(/                        &
rath_t('SO2       ','H2O2      ','NULL0     ','          ','          ', & !HSO3+H2O2(aq)
'          ', 0.000, 0.000, 0.000, 0.000),                               &
rath_t('SO2       ','O3        ','NULL1     ','          ','          ', & !HSO3+O3(aq)
'          ', 0.000, 0.000, 0.000, 0.000),                               &
rath_t('SO2       ','O3        ','NULL2     ','          ','          ', & !SO3+O3(aq)
'          ', 0.000, 0.000, 0.000, 0.000),                               &
rath_t('N2O5      ','          ','HONO2     ','          ','          ', & ! Heterogenous
'          ', 2.000, 0.000, 0.000, 0.000),                               &
rath_t('HO2       ','          ','H2O2      ','          ','          ', & ! Heterogenous
'          ', 0.500, 0.000, 0.000, 0.000)                                &
/)


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

TYPE(RATJ_T), DIMENSION(10) :: ratj_defs_isop=(/                        &
ratj_t('ISOOH     ','PHOTON    ','OH        ','MACR      ','HCHO      ',&
     'HO2       ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmhp      '),             &
ratj_t('ISON      ','PHOTON    ','NO2       ','MACR      ','HCHO      ',&
     'HO2       ', 0.0, 0.0, 0.0, 0.0, 100.0,'jiprn     '),             &
ratj_t('MACR      ','PHOTON    ','MeCO3     ','HCHO      ','CO        ',&
     'HO2       ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmacr     '),             &
ratj_t('MPAN      ','PHOTON    ','MACRO2    ','NO2       ','          ',&
     '          ', 0.0, 0.0, 0.0, 0.0, 100.0,'jpan      '),             &
ratj_t('MACROOH   ','PHOTON    ','OH        ','HO2       ','OH        ',&
     'HO2       ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmacro    '),             &
ratj_t('MACROOH   ','PHOTON    ','HACET     ','CO        ','MGLY      ',&
     'HCHO      ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmacro    '),             &
ratj_t('HACET     ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',&
     '          ', 0.0, 0.0, 0.0, 0.0, 100.0,'jhacet    '),             &
ratj_t('MGLY      ','PHOTON    ','MeCO3     ','CO        ','HO2       ',&
     '          ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmkal     '),             &
ratj_t('NALD      ','PHOTON    ','HCHO      ','CO        ','NO2       ',&
     'HO2       ', 0.0, 0.0, 0.0, 0.0, 100.0,'jaceta    '),             &
ratj_t('MeCO3H    ','PHOTON    ','MeOO      ','OH        ','          ',&
     '          ', 0.0, 0.0, 0.0, 0.0, 100.0,'jmeco3h   ')              &
 /)

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


TYPE(RATT_T), DIMENSION(2) :: ratt_defs_isop=(/                         &
ratt_t('MACRO2    ','NO2       ','MPAN      ','m         ',             &
  0.3, 2.70E-28, -7.10,     0.0, 1.20E-11, -0.90, 0.0, 0.000, 0.000),   &
ratt_t('MPAN      ','m         ','MACRO2    ','NO2       ',             &
  0.3, 4.90E-03,  0.00, 12100.0, 5.40E+16, 0.00, 13830.0, 0.000, 0.000) &
  /)

TYPE(RATT_T), DIMENSION(1) :: ratt_defs_aer=(/                          &
ratt_t('SO2       ','OH        ','HO2       ','H2SO4     ',             &
  0.6, 3.00E-31, -3.30,     0.0, 1.50E-12, 0.00, 0.0, 0.000, 0.000)     &
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


REAL, DIMENSION(6,5,ndry_isop) :: depvel_defs_isop=RESHAPE((/           &
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
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! ISOOH
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! ISON
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! MACR
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! MACROOH
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & ! MPAN
        0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & 
        0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & 
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & 
        0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! HACET
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! MGLY
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! NALD
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! HCOOH
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! MeCO3H
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! MeCO2H
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! MeOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & 
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & 
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & 
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01   & 
        /),(/  6,  5, ndry_isop/) )

REAL, DIMENSION(6,5,ndry_aer) :: depvel_defs_aer=RESHAPE((/             &
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
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! ISOOH
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! ISON
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! MACR
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! MACROOH
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & ! MPAN
        0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & 
        0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & 
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & 
        0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! HACET
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! MGLY
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & ! NALD
        0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  &
        0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! HCOOH
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & ! MeCO3H
        0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  &
        0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  &
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  &
        0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  &
        0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & ! MeCO2H
        1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  &
        0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  &
        0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  &
        0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  &
        0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & ! SO2
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & ! H2SO4
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & ! MSA
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  &
        1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! DMSO
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  &
        0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & ! NH3
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & ! MeOH
        0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & 
        0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & 
        0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & 
        0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & 
        0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & ! monoterp
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & ! SEC_ORG
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  &
        0.10,  0.10,  0.10,  0.10,  0.10,  0.10   &
        /),(/  6,  5,  ndry_aer/) )


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


REAL, DIMENSION(6,nwet_aer) :: henry_defs_aer=RESHAPE((/                &
 0.1130E-01, 0.2300E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! O3
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
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! ISOOH
 0.3000E+04, 0.7400E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! ISON
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MACROOH
 0.1400E+03, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HACET
 0.3500E+04, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MGLY
 0.6900E+04, 0.5600E+04, 0.1800E-03,-0.1510E+04, 0.0000E+00, 0.0000E+00,&  ! HCOOH
 0.7500E+03, 0.5300E+04, 0.6300E-08, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeCO3H
 0.4700E+04, 0.6000E+04, 0.1800E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeCO2H
 0.1230E+01, 0.3020E+04, 0.1230E-01, 0.2010E+04, 0.6000E-07, 0.1120E+04,&  ! SO2
 0.5000E+05, 0.6425E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! DMSO
 0.1000E+07, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! NH3 (H*)
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HO2
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeOO
 0.2300E+03, 0.4900E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeOH 
 0.1000E+06, 0.1200E+02, 0.0000E-00, 0.0000E+00, 0.0000E+00, 0.0000E+00 &  ! SEC_ORG
  /),(/  6,  nwet_aer/) )

REAL, DIMENSION(6,nwet_isop) :: henry_defs_isop=RESHAPE((/              &
 0.1130E-01, 0.2300E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! O3
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
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! ISOOH
 0.3000E+04, 0.7400E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! ISON
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MACROOH
 0.1400E+03, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HACET
 0.3500E+04, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MGLY
 0.6900E+04, 0.5600E+04, 0.1800E-03,-0.1510E+04, 0.0000E+00, 0.0000E+00,&  ! HCOOH
 0.7500E+03, 0.5300E+04, 0.6300E-08, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeCO3H
 0.4700E+04, 0.6000E+04, 0.1800E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeCO2H
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! HO2
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  ! MeOO
 0.2300E+03, 0.4900E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00 &  ! MeOH           
  /),(/  6, nwet_isop/) )

INTERFACE ukca_init_tropisop
  MODULE PROCEDURE ukca_init_tropisop
END INTERFACE ukca_init_tropisop
PUBLIC ukca_init_tropisop

CONTAINS

! ######################################################################

SUBROUTINE UKCA_INIT_TROPISOP

! To initialise chemistry files from specified species and rate arrays
! using the defined logicals. Composes chemistry definition arrays from
! components using logical definitions from RUN_UKCA namelist held in UKCA 
! module ukca_option_mod

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

CHARACTER (LEN=72) :: cmessage
INTEGER            :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_INIT_TROPISOP',zhook_in,zhook_handle)

ALLOCATE(chch_defs_tropisop(jpspec))
ALLOCATE(ratb_defs_tropisop(jpbk))
ALLOCATE(ratj_defs_tropisop(jppj))
ALLOCATE(ratt_defs_tropisop(jptk))
ALLOCATE(rath_defs_tropisop(jphk))
ALLOCATE(depvel_defs_tropisop(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_tropisop(6,jpdw))

 IF (L_ukca_tropisop) THEN
   IF (L_ukca_achem) THEN
! Standard tropospheric chemistry plus isoprene and aerosol chemistry 
     chch_defs_tropisop = (/chch_defs_aer/)
     ratb_defs_tropisop = (/ratb_defs_a, ratb_defs_b, ratb_defs_c,      &
                            ratb_defs_isop, ratb_defs_aer/)
     ratj_defs_tropisop = (/ratj_defs_std_trop1, ratj_defs_isop/)
     ratt_defs_tropisop = (/ratt_defs_std_trop1, ratt_defs_isop,        &
                            ratt_defs_aer/)
     rath_defs_tropisop = (/rath_defs_std_trop1, rath_defs_isop,        &
                            rath_defs_aer/)
     depvel_defs_tropisop = depvel_defs_aer
     henry_defs_tropisop  = henry_defs_aer
   ELSE
! Standard tropospheric chemistry plus isoprene chemistry
     chch_defs_tropisop = (/ chch_defs_isop /)
     ratb_defs_tropisop = (/ ratb_defs_a, ratb_defs_b, ratb_defs_c,     &
                             ratb_defs_isop /)
     ratj_defs_tropisop = (/ ratj_defs_std_trop1, ratj_defs_isop /)
     ratt_defs_tropisop = (/ ratt_defs_std_trop1, ratt_defs_isop /)
     rath_defs_tropisop = (/ rath_defs_std_trop1, rath_defs_isop /)
     depvel_defs_tropisop = depvel_defs_isop
     henry_defs_tropisop = henry_defs_isop
   ENDIF
 ELSE
   cmessage = ' No relevant definitions to compose chemistry'
   errcode = 1
   CALL EREPORT('UKCA_INIT_TROPISOP',errcode,cmessage)
 ENDIF

IF (lhook) CALL dr_hook('UKCA_INIT_TROPISOP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INIT_TROPISOP

END MODULE UKCA_CHEM_TROPISOP
