! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing chemistry
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
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE UKCA_CHEM1_DAT
USE ukca_option_mod,    ONLY:  L_ukca_aerchem, L_ukca_raq, L_ukca_trop, &
                               L_ukca_tropisop, L_ukca_std_trop,        &
                               L_ukca_strattrop, L_ukca_stratcfc,       &
                               L_ukca_strat, jpspec, jpbk, jptk, jphk,  &
                               jppj, jpdd, jpdw
USE UKCA_CHEM_DEFS_MOD, ONLY:  chch_t, ratb_t, rath_t, ratj_t, ratt_t,  &
                               chch_defs, ratb_defs, rath_defs,         &
                               ratj_defs, ratt_defs, depvel_defs,       &
                               henry_defs
! B-E options:
USE UKCA_CHEM_AER,      ONLY:  chch_defs_aer, ratb_defs_aer,            &
                               ratj_defs_aer, rath_defs_aer,            &
                               ratt_defs_aer, depvel_defs_aer,          &
                               henry_defs_aer

USE UKCA_CHEM_RAQ,      ONLY:  chch_defs_raq, ratb_defs_raq,            &
                               ratj_defs_raq, rath_defs_raq,            &
                               ratt_defs_raq, depvel_defs_raq,          &
                               henry_defs_raq
! N-R options:
USE UKCA_CHEM_STD_TROP, ONLY:  chch_defs_std_trop, ratb_defs_std_trop,  &
                               ratj_defs_std_trop, rath_defs_std_trop,  &
                               ratt_defs_std_trop, depvel_defs_std_trop,&
                               henry_defs_std_trop, ukca_init_std_trop

USE UKCA_CHEM_TROPISOP, ONLY:  chch_defs_tropisop, ratb_defs_tropisop,  &
                               ratj_defs_tropisop, rath_defs_tropisop,  &
                               ratt_defs_tropisop, depvel_defs_tropisop,&
                               henry_defs_tropisop, ukca_init_tropisop

USE UKCA_CHEM_STRAT,    ONLY:  chch_defs_strat, ratb_defs_strat,        &
                               ratj_defs_strat, rath_defs_strat,        &
                               ratt_defs_strat, depvel_defs_strat,      &
                               henry_defs_strat, ukca_init_strat

USE UKCA_CHEM_STRATTROP,ONLY:  chch_defs_strattrop, ratb_defs_strattrop,&
                               ratj_defs_strattrop, rath_defs_strattrop,&
                               ratt_defs_strattrop,                     &
                               depvel_defs_strattrop,                   &
                               henry_defs_strattrop, ukca_init_strattrop
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
PRIVATE
SAVE

PUBLIC ukca_chem1_init

! 1) Standard tropospheric chemistry
! ========================================================

TYPE(CHCH_T), DIMENSION( 46), PUBLIC :: chch_defs_trop=(/              &
chch_t(  1,'O(3P)     ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t(  2,'O(1D)     ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t(  3,'O3        ',  1,'TR        ','          ',  1,  0,  0),    &
chch_t(  4,'NO        ',  1,'TR        ','          ',  1,  0,  1),    &
chch_t(  5,'NO3       ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t(  6,'NO2       ',  1,'TR        ','          ',  1,  0,  1),    &
chch_t(  7,'N2O5      ',  2,'TR        ','          ',  1,  1,  0),    &
chch_t(  8,'HO2NO2    ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t(  9,'HONO2     ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 10,'OH        ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 11,'HO2       ',  1,'SS        ','          ',  0,  1,  0),    &
chch_t( 12,'H2        ',  1,'CT        ','          ',  0,  0,  0),    &
chch_t( 13,'H2O2      ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 14,'CH4       ',  1,'TR        ','          ',  1,  0,  1),    &
chch_t( 15,'CO        ',  1,'TR        ','          ',  1,  0,  1),    &
chch_t( 16,'CO2       ',  1,'CT        ','          ',  0,  0,  0),    &
chch_t( 17,'HCHO      ',  1,'TR        ','          ',  1,  1,  1),    &
chch_t( 18,'MeOO      ',  1,'SS        ','          ',  0,  1,  0),    &
chch_t( 19,'H2O       ',  1,'CF        ','          ',  0,  0,  0),    &
chch_t( 20,'MeOOH     ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 21,'HONO      ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 22,'O2        ',  1,'CT        ','          ',  0,  0,  0),    &
chch_t( 23,'N2        ',  1,'CT        ','          ',  0,  0,  0),    &
chch_t( 24,'C2H6      ',  1,'TR        ','          ',  0,  0,  1),    &
chch_t( 25,'EtOO      ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 26,'EtOOH     ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 27,'MeCHO     ',  1,'TR        ','          ',  1,  0,  1),    &
chch_t( 28,'MeCO3     ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 29,'PAN       ',  1,'TR        ','          ',  1,  0,  0),    &
chch_t( 30,'C3H8      ',  1,'TR        ','          ',  0,  0,  1),    &
chch_t( 31,'n-PrOO    ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 32,'i-PrOO    ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 33,'n-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 34,'i-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 35,'EtCHO     ',  1,'TR        ','          ',  1,  0,  0),    &
chch_t( 36,'EtCO3     ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 37,'Me2CO     ',  1,'TR        ','          ',  0,  0,  1),    &
chch_t( 38,'MeCOCH2OO ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 39,'MeCOCH2OOH',  1,'TR        ','          ',  1,  1,  0),    &
chch_t( 40,'PPAN      ',  1,'TR        ','          ',  1,  0,  0),    &
chch_t( 41,'MeONO2    ',  1,'TR        ','          ',  0,  0,  0),    &
chch_t( 42,'O(3P)S    ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 43,'O(1D)S    ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 44,'O3S       ',  1,'TR        ','          ',  1,  0,  0),    &
chch_t( 45,'OHS       ',  1,'SS        ','          ',  0,  0,  0),    &
chch_t( 46,'HO2S      ',  1,'SS        ','          ',  0,  1,  0)     &
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
ratb_t('NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ', &
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
ratb_t('O(1D)     ','O2        ','O(3P)     ','O2        ','          ', &
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
ratb_t('OH        ','H2O2      ','H2O       ','HO2       ','          ', &
'          ',2.90E-12,  0.00,    160.00, 0.000, 0.000, 0.000, 0.000)     &
 /)

TYPE(RATB_T), DIMENSION(  28),PARAMETER :: ratb_defs_c=(/                &
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
ratb_t('OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ', &
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
ratb_t('OH        ','n-PrOOH   ','EtCHO     ','H2O       ','OH        ', &
'          ',1.10E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','i-PrOOH   ','i-PrOO    ','H2O       ','          ', &
'          ',1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ', &
'          ',1.66E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(3P)     ','NO2       ','NO        ','O2        ','          ', &
'          ',5.50E-12,  0.00,   -188.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('HO2S      ','O3S       ','HO2       ','O2        ','          ', &
'          ',2.03E-16,  4.57,   -693.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('OHS       ','O3S       ','OH        ','O2        ','          ', &
'          ',1.70E-12,  0.00,    940.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)S    ','H2O       ','H2O       ','          ','          ', &
'          ',2.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)S    ','N2        ','O(3P)S    ','N2        ','          ', &
'          ',2.10E-11,  0.00,   -115.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('O(1D)S    ','O2        ','O(3P)S    ','O2        ','          ', &
'          ',3.20E-11,  0.00,    -67.00, 0.000, 0.000, 0.000, 0.000)     &
 /)

TYPE(RATB_T), DIMENSION(  88    ),PARAMETER :: ratb_defs_trop=         &
  (/ ratb_defs_a, ratb_defs_b, ratb_defs_c /)


TYPE(RATH_T), DIMENSION(   0) :: rath_defs_trop


TYPE(RATJ_T), DIMENSION(  20) :: ratj_defs_trop=(/                      &
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
ratj_t('HONO      ','PHOTON    ','OH        ','NO        ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhono     ') ,       &
ratj_t('EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jetcho    ') ,       &
ratj_t('Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jacetone  ') ,       &
ratj_t('MeONO2    ','PHOTON    ','HO2       ','HCHO      ','NO2       ',&
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmena     ')         &
  /)


TYPE(RATT_T), DIMENSION(  14) :: ratt_defs_trop=(/                     &
ratt_t('HO2       ','HO2       ','H2O2      ','O2        ',            &
  0.0, 1.90E-33,  0.00,  -980.0, 0.00E+00,  0.00, 0.0, 0.0, 0.0) ,     &
ratt_t('HO2       ','NO2       ','HO2NO2    ','m         ',            &
  0.6, 1.80E-31, -3.20,     0.0, 4.70E-12,  0.00, 0.0, 0.0, 0.0) ,     &
ratt_t('HO2NO2    ','m         ','HO2       ','NO2       ',            &
  0.6, 4.10E-05,  0.00, 10650.0, 4.80E+15,  0.00, 11170.0, 0.0, 0.0),  &
ratt_t('MeCO3     ','NO2       ','PAN       ','m         ',            &
  0.3, 2.70E-28, -7.10,     0.0, 1.20E-11, -0.90, 0.0, 0.0, 0.0) ,     &
ratt_t('PAN       ','m         ','MeCO3     ','NO2       ',            &
  0.3, 4.90E-03,  0.00, 12100.0, 5.40E+16,  0.00, 13830.0, 0.0, 0.0),  &
ratt_t('N2O5      ','m         ','NO2       ','NO3       ',            &
  0.3, 1.30E-03, -3.50, 11000.0, 9.70E+14,  0.10, 11080.0, 0.0, 0.0),  &
ratt_t('NO2       ','NO3       ','N2O5      ','m         ',            &
  0.3, 3.60E-30, -4.10,     0.0, 1.90E-12,  0.20, 0.0, 0.0, 0.0) ,     &
ratt_t('O(3P)     ','O2        ','O3        ','m         ',            &
  0.0, 5.70E-34, -2.60,     0.0, 0.00E+00,  0.00, 0.0, 0.0, 0.0) ,     &
ratt_t('OH        ','NO        ','HONO      ','m         ',            &
  1420.0, 7.40E-31, -2.40,  0.0, 3.30E-11, -0.30, 0.0, 0.0, 0.0) ,     &
ratt_t('OH        ','NO2       ','HONO2     ','m         ',            &
  0.4, 3.30E-30, -3.00,     0.0, 4.10E-11,  0.00, 0.0, 0.0, 0.0) ,     &
ratt_t('OH        ','OH        ','H2O2      ','m         ',            &
  0.5, 6.90E-31, -0.80,     0.0, 2.60E-11,  0.00, 0.0, 0.0, 0.0) ,     &
ratt_t('EtCO3     ','NO2       ','PPAN      ','m         ',            &
  0.3, 2.70E-28, -7.10,     0.0, 1.20E-11, -0.90, 0.0, 0.0, 0.0) ,     &
ratt_t('PPAN      ','m         ','EtCO3     ','NO2       ',            &
  0.4, 1.70E-03,  0.00, 11280.0, 8.30E+16,  0.00, 13940.0, 0.0, 0.0),  &
ratt_t('O(3P)S    ','O2        ','O3S       ','m         ',            &
  0.0, 5.70E-34, -2.60,     0.0, 0.00E+00,  0.00, 0.0, 0.0, 0.0)       &
  /)


REAL, DIMENSION(  6,  5, 22) :: depvel_defs_trop=RESHAPE((/            &
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,                             &
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                             &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                             &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                             &
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.14,  0.01,  0.07,  0.01,  0.01,  0.01,                             &
  0.10,  0.01,  0.06,  0.01,  0.01,  0.01,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,                             &
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,                             &
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,                             &
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,                             &
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                             &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                             &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                             &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,                             &
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,                             &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                             &
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                             &
  1.00,  0.01,  0.50,  0.01,  0.01,  0.01,                             &
  0.71,  0.03,  0.37,  0.03,  0.03,  0.03,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,                             &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                             &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                             &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                             &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                             &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,                             &
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,                             &
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,                             &
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00,                             &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                             &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                             &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                             &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                             &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,                             &
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,                             &
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                             &
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,                             &
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,                             &
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,                             &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                             &
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,                             &
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,                             &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                             &
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00,                             &
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07,                             &
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                             &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                             &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                             &
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07                              &
  /),(/  6,  5, 22/) )

REAL, DIMENSION(  6, 15) :: henry_defs_trop=RESHAPE((/                  &
 0.2000E+01, 0.2000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! NO3
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! N2O5
 0.1300E+05, 0.6900E+04, 0.1000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! HO2NO2
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! HONO2
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! HO2
 0.8300E+05, 0.7400E+04, 0.2400E-11,-0.3730E+04, 0.0000E+00, 0.0000E+00,&   ! H2O2
 0.3300E+04, 0.6500E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! HCHO
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! MeOO
 0.3100E+03, 0.5000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! MeOOH
 0.5000E+02, 0.4900E+04, 0.5600E-03,-0.1260E+04, 0.0000E+00, 0.0000E+00,&   ! HONO
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! EtOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! n-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! i-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   ! MeCOCH2OOH
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00 &   ! HO2S
  /),(/  6, 15/) )



CONTAINS

SUBROUTINE UKCA_CHEM1_INIT

! To select the correct chemical definitions according
!  to the UKCA logicals set in the UMUI

USE ASAD_MOD,    ONLY: jddept, jddepc
USE ereport_mod, ONLY : ereport
USE Control_Max_Sizes
IMPLICIT NONE

INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: i
INTEGER :: ia
INTEGER :: ierr
INTEGER :: n_ss
INTEGER :: isizes(7)
CHARACTER(LEN=72) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_CHEM1_DAT:UKCA_CHEM1_INIT',zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(chch_defs)) ALLOCATE(chch_defs(jpspec))
IF (.NOT. ALLOCATED(ratb_defs)) ALLOCATE(ratb_defs(jpbk))
IF (.NOT. ALLOCATED(rath_defs)) ALLOCATE(rath_defs(jphk))
IF (.NOT. ALLOCATED(ratj_defs)) ALLOCATE(ratj_defs(jppj))
IF (.NOT. ALLOCATED(ratt_defs)) ALLOCATE(ratt_defs(jptk))
IF (.NOT. ALLOCATED(depvel_defs)) ALLOCATE(depvel_defs(jddept,jddepc,jpdd))
IF (.NOT. ALLOCATED(henry_defs)) ALLOCATE(henry_defs(6,jpdw))

IF (L_ukca_trop) THEN                  ! B-E option
  isizes(1) = SIZE(chch_defs_trop)
  isizes(2) = SIZE(ratb_defs_trop)
  isizes(3) = SIZE(rath_defs_trop)
  isizes(4) = SIZE(ratj_defs_trop)
  isizes(5) = SIZE(ratt_defs_trop)
  isizes(6) = SIZE(depvel_defs_trop,3)
  isizes(7) = SIZE(henry_defs_trop,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_trop
  ratb_defs=ratb_defs_trop
  rath_defs=rath_defs_trop
  ratj_defs=ratj_defs_trop
  ratt_defs=ratt_defs_trop
  depvel_defs=depvel_defs_trop
  henry_defs=henry_defs_trop
ELSEIF (L_UKCA_TROPISOP) THEN
  CALL UKCA_INIT_TROPISOP()
  isizes(1) = SIZE(chch_defs_tropisop)
  isizes(2) = SIZE(ratb_defs_tropisop)
  isizes(3) = SIZE(rath_defs_tropisop)
  isizes(4) = SIZE(ratj_defs_tropisop)
  isizes(5) = SIZE(ratt_defs_tropisop)
  isizes(6) = SIZE(depvel_defs_tropisop,3)
  isizes(7) = SIZE(henry_defs_tropisop,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_tropisop
  ratb_defs=ratb_defs_tropisop
  rath_defs=rath_defs_tropisop
  ratj_defs=ratj_defs_tropisop
  ratt_defs=ratt_defs_tropisop
  depvel_defs=depvel_defs_tropisop
  henry_defs=henry_defs_tropisop
ELSEIF (L_ukca_raq) THEN               ! B-E option
  isizes(1) = SIZE(chch_defs_raq)
  isizes(2) = SIZE(ratb_defs_raq)
  isizes(3) = SIZE(rath_defs_raq)
  isizes(4) = SIZE(ratj_defs_raq)
  isizes(5) = SIZE(ratt_defs_raq)
  isizes(6) = SIZE(depvel_defs_raq,3)
  isizes(7) = SIZE(henry_defs_raq,2)
  CALL CHECK_DIMS
  chch_defs   = chch_defs_raq
  ratb_defs   = ratb_defs_raq
  rath_defs   = rath_defs_raq
  ratj_defs   = ratj_defs_raq
  ratt_defs   = ratt_defs_raq
  depvel_defs = depvel_defs_raq
  henry_defs  = henry_defs_raq
ELSEIF (L_ukca_aerchem) THEN           ! B-E option
  isizes(1) = SIZE(chch_defs_aer)
  isizes(2) = SIZE(ratb_defs_aer)
  isizes(3) = SIZE(rath_defs_aer)
  isizes(4) = SIZE(ratj_defs_aer)
  isizes(5) = SIZE(ratt_defs_aer)
  isizes(6) = SIZE(depvel_defs_aer,3)
  isizes(7) = SIZE(henry_defs_aer,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_aer
  ratb_defs=ratb_defs_aer
  rath_defs=rath_defs_aer
  ratj_defs=ratj_defs_aer
  ratt_defs=ratt_defs_aer
  depvel_defs=depvel_defs_aer
  henry_defs=henry_defs_aer
ELSEIF (L_ukca_std_trop) THEN
  CALL UKCA_INIT_STD_TROP
  isizes(1) = SIZE(chch_defs_std_trop)
  isizes(2) = SIZE(ratb_defs_std_trop)
  isizes(3) = SIZE(rath_defs_std_trop)
  isizes(4) = SIZE(ratj_defs_std_trop)
  isizes(5) = SIZE(ratt_defs_std_trop)
  isizes(6) = SIZE(depvel_defs_std_trop,3)
  isizes(7) = SIZE(henry_defs_std_trop,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_std_trop
  ratb_defs=ratb_defs_std_trop
  rath_defs=rath_defs_std_trop
  ratj_defs=ratj_defs_std_trop
  ratt_defs=ratt_defs_std_trop
  depvel_defs=depvel_defs_std_trop
  henry_defs=henry_defs_std_trop
ELSEIF (L_ukca_strat) THEN
  CALL UKCA_INIT_STRAT()
  isizes(1) = SIZE(chch_defs_strat)
  isizes(2) = SIZE(ratb_defs_strat)
  isizes(3) = SIZE(rath_defs_strat)
  isizes(4) = SIZE(ratj_defs_strat)
  isizes(5) = SIZE(ratt_defs_strat)
  isizes(6) = SIZE(depvel_defs_strat,3)
  isizes(7) = SIZE(henry_defs_strat,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_strat
  ratb_defs=ratb_defs_strat
  rath_defs=rath_defs_strat
  ratj_defs=ratj_defs_strat
  ratt_defs=ratt_defs_strat
  depvel_defs=depvel_defs_strat
  henry_defs=henry_defs_strat
ELSEIF (L_UKCA_STRATCFC) THEN
  cmessage='Chemistry definition arrays not included for STRATCFC'
  errcode=1
  CALL EREPORT('UKCA_CHEM1_INIT',errcode,cmessage)
ELSEIF (L_UKCA_STRATTROP) THEN
  CALL UKCA_INIT_STRATTROP()
  isizes(1) = SIZE(chch_defs_strattrop)
  isizes(2) = SIZE(ratb_defs_strattrop)
  isizes(3) = SIZE(rath_defs_strattrop)
  isizes(4) = SIZE(ratj_defs_strattrop)
  isizes(5) = SIZE(ratt_defs_strattrop)
  isizes(6) = SIZE(depvel_defs_strattrop,3)
  isizes(7) = SIZE(henry_defs_strattrop,2)
  CALL CHECK_DIMS
  chch_defs=chch_defs_strattrop
  ratb_defs=ratb_defs_strattrop
  rath_defs=rath_defs_strattrop
  ratj_defs=ratj_defs_strattrop
  ratt_defs=ratt_defs_strattrop
  depvel_defs=depvel_defs_strattrop
  henry_defs=henry_defs_strattrop
ELSE
  cmessage = ' Unknown chemistry, unable to set chch_defs'
  errcode = 1
  CALL EREPORT('UKCA_CHEM1_INIT',errcode,cmessage)
ENDIF

IF (lhook) CALL dr_hook('UKCA_CHEM1_DAT:UKCA_CHEM1_INIT',zhook_out,zhook_handle)
RETURN

CONTAINS

SUBROUTINE CHECK_DIMS
IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_CHEM1_DAT:CHECK_DIMS',zhook_in,zhook_handle)

  IF (isizes(1) /= jpspec) THEN
    cmessage= ' chch_defs: size not equal to jpspec'
    write(6,*) 'SIZES: ',isizes(1),' AND ',jpspec,' ARE NOT EQUAL!'
    ierr=91
  ELSE IF (isizes(2) /= jpbk) THEN
    cmessage= ' ratb_defs: size not equal to jpbk'
    write(6,*) 'SIZES: ',isizes(2),' AND ',jpbk,' ARE NOT EQUAL!'
    ierr=92
  ELSE IF (isizes(3) /= jphk) THEN
    cmessage= ' rath_defs: size not equal to jphk'
    write(6,*) 'SIZES: ',isizes(3),' AND ',jphk,' ARE NOT EQUAL!'
    ierr=93
  ELSE IF (isizes(4) /= jppj) THEN
    cmessage= ' ratj_defs: size not equal to jppj'
    write(6,*) 'SIZES: ',isizes(4),' AND ',jppj,' ARE NOT EQUAL!'
    ierr=94
  ELSE IF (isizes(5) /= jptk) THEN
    cmessage= ' ratt_defs: size not equal to jptk'
    write(6,*) 'SIZES: ',isizes(5),' AND ',jptk,' ARE NOT EQUAL!'
    ierr=95
  ELSE IF (isizes(6) /= jpdd) THEN
    cmessage= ' depvel_defs: size of dimn=3 not equal to jpdd'
    write(6,*) 'SIZES: ',isizes(6),' AND ',jpdd,' ARE NOT EQUAL!'
    ierr=96
  ELSE IF (isizes(7) /= jpdw) THEN
    cmessage= ' henry_defs: size of dimn=2 not equal to jpdw'
    write(6,*) 'SIZES: ',isizes(7),' AND ',jpdw,' ARE NOT EQUAL!'
    ierr=97
  ELSE
    ierr=0
  ENDIF
  IF (ierr > 1) THEN
    CALL EREPORT('UKCA_CHEM1_INIT.check_dims',ierr,cmessage)
  ENDIF

 IF (lhook) CALL dr_hook('UKCA_CHEM1_DAT:CHECK_DIMS',zhook_out,zhook_handle)
 RETURN
 END SUBROUTINE CHECK_DIMS

END SUBROUTINE ukca_chem1_init

END MODULE
