! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
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
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE UKCA_CHEM_STRATTROP

USE UKCA_CHEM_DEFS_MOD, ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t,  &
                              chch_defs, ratb_defs, rath_defs,         &
                              ratj_defs, ratt_defs, depvel_defs,       &
                              henry_defs
USE ukca_option_mod,    ONLY: l_ukca_strattrop, l_ukca_het_psc,        &
                              l_ukca_trophet, l_ukca_achem, jpspec,    &
                              jpbk, jptk, jphk, jppj, jpdd, jpdw
USE ASAD_MOD,           ONLY: jddepc, jddept
USE Control_Max_Sizes
IMPLICIT NONE
PRIVATE

! Standard troposphere-stratosphere chemistry with optional aerosol chemistry
! ===========================================================================

TYPE(CHCH_T), ALLOCATABLE, PUBLIC :: chch_defs_strattrop(:)
TYPE(RATB_T), ALLOCATABLE, PUBLIC :: ratb_defs_strattrop(:)
TYPE(RATT_T), ALLOCATABLE, PUBLIC :: ratt_defs_strattrop(:)
TYPE(RATJ_T), ALLOCATABLE, PUBLIC :: ratj_defs_strattrop(:)
TYPE(RATH_T), ALLOCATABLE, PUBLIC :: rath_defs_strattrop(:)
REAL,         ALLOCATABLE, PUBLIC :: depvel_defs_strattrop(:,:,:)
REAL,         ALLOCATABLE, PUBLIC :: henry_defs_strattrop(:,:)

! No of heterogenous processes
INTEGER, PARAMETER, PUBLIC :: nhet_strattrop = 5         ! Stratos chemistry
INTEGER, PARAMETER, PUBLIC :: nhet_st_aer    = 3         ! Aqueous phase
INTEGER, PARAMETER, PUBLIC :: nhet_st_tpht   = 2         ! trophet rxns

! No of dry deposited species
INTEGER, PARAMETER, PUBLIC :: ndry_strattrop = 36        ! Stratos chemistry
INTEGER, PARAMETER, PUBLIC :: ndry_st_aer    = 5         ! Aerosol chemistry

! No of wet deposited species
INTEGER, PARAMETER, PUBLIC :: nwet_strattrop = 29        ! Stratos chemistry
INTEGER, PARAMETER, PUBLIC :: nwet_st_aer    = 34        ! Aerosol chemistry

 
! ATA NLA CheST Chemistry v1.2
TYPE(CHCH_T), PUBLIC :: chch_defs_strattrop_chem(1:75)=(/            &
!   1
chch_t(  1,'O(3P)     ',  1,'TR        ','Ox        ',  0,  0,  0),  & 
!   2
chch_t(  2,'O(1D)     ',  1,'SS        ','Ox        ',  0,  0,  0),  & 
!   3 DD: 1,
chch_t(  3,'O3        ',  1,'TR        ','Ox        ',  1,  0,  0),  & 
!   4
chch_t(  4,'N         ',  1,'TR        ','NOx       ',  0,  0,  0),  & 
!   5 DD: 2,
chch_t(  5,'NO        ',  1,'TR        ','NOx       ',  1,  0,  0),  & 
!   6 DD: 3,WD: 1,  
chch_t(  6,'NO3       ',  1,'TR        ','NOx       ',  1,  1,  0),  & 
!   7 DD: 4,      EM: 1     
chch_t(  7,'NO2       ',  1,'TR        ','NOx       ',  1,  0,  1),  & 
!   8 DD: 5,WD: 2,  
chch_t(  8,'N2O5      ',  1,'TR        ','          ',  1,  1,  0),  & 
!   9 DD: 6,WD: 3,          
chch_t(  9,'HO2NO2    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  10 DD: 7,WD: 4,           
chch_t( 10,'HONO2     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  11 DD: 8,WD: 5,          
chch_t( 11,'H2O2      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  12             EM: 2      
chch_t( 12,'CH4       ',  1,'TR        ','          ',  0,  0,  2),  & 
!  13 DD: 9,      EM:  3
chch_t( 13,'CO        ',  1,'TR        ','          ',  1,  0,  3),  & 
!  14 DD:10,WD: 6,EM:  4
chch_t( 14,'HCHO      ',  1,'TR        ','          ',  1,  1,  4),  & 
!  15       WD: 7,  
chch_t( 15,'MeOO      ',  1,'TR        ','          ',  0,  1,  0),  & 
!  16 DD:11,WD: 8,     
chch_t( 16,'MeOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  17     
chch_t( 17,'H         ',  1,'TR        ','HOx       ',  0,  0,  0),  & 
!  18
chch_t( 18,'H2O       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  19
chch_t( 19,'OH        ',  1,'TR        ','HOx       ',  0,  0,  0),  & 
!  20       WD: 9,
chch_t( 20,'HO2       ',  1,'TR        ','HOx       ',  0,  1,  0),  & 
!  21     
chch_t( 21,'Cl        ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  22
chch_t( 22,'Cl2O2     ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  23
chch_t( 23,'ClO       ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  24
chch_t( 24,'OClO      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  25
chch_t( 25,'Br        ',  1,'TR        ','Brx       ',  0,  0,  0),  & 
!  26
chch_t( 26,'BrO       ',  1,'TR        ','Brx       ',  0,  0,  0),  & 
!  27
chch_t( 27,'BrCl      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  28       WD:10, 
chch_t( 28,'BrONO2    ',  1,'TR        ','          ',  0,  1,  0),  & 
!  29       
chch_t( 29,'N2O       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  30 DD:12,WD:11,
chch_t( 30,'HCl       ',  1,'TR        ','          ',  1,  1,  0),  & 
!  31 DD:13,WD:12,
chch_t( 31,'HOCl      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  32 DD:14,WD:13,
chch_t( 32,'HBr       ',  1,'TR        ','          ',  1,  1,  0),  & 
!  33 DD:15,WD:14,
chch_t( 33,'HOBr      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  34       WD:15, 
chch_t( 34,'ClONO2    ',  1,'TR        ','          ',  0,  1,  0),  & 
!  35      
chch_t( 35,'CFCl3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  36
chch_t( 36,'CF2Cl2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  37
chch_t( 37,'MeBr      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  38 DD:16,WD:16,
chch_t( 38,'HONO      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  39             EM: 5
chch_t( 39,'C2H6      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  40
chch_t( 40,'EtOO      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  41 DD:17,WD: 17,
chch_t( 41,'EtOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  42 DD:18,       EM: 6
chch_t( 42,'MeCHO     ',  1,'TR        ','          ',  1,  0,  1),  & 
!  43
chch_t( 43,'MeCO3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  44 DD:19, 
chch_t( 44,'PAN       ',  1,'TR        ','          ',  1,  0,  0),  & 
!  45              EM: 7
chch_t( 45,'C3H8      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  46
chch_t( 46,'n-PrOO    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  47
chch_t( 47,'i-PrOO    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  48 DD:20,WD:18,
chch_t( 48,'n-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  49 DD:21,WD:19,
chch_t( 49,'i-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  50 DD:22,
chch_t( 50,'EtCHO     ',  1,'TR        ','          ',  1,  0,  0),  & 
!  51
chch_t( 51,'EtCO3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  52             EM: 8 
chch_t( 52,'Me2CO     ',  1,'TR        ','          ',  0,  0,  1),  & 
!  53
chch_t( 53,'MeCOCH2OO ',  1,'TR        ','          ',  0,  0,  0),  & 
!  54 DD:23,WD:20,
chch_t( 54,'MeCOCH2OOH',  1,'TR        ','          ',  1,  1,  0),  & 
!  55 DD:24,
chch_t( 55,'PPAN      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  56
chch_t( 56,'MeONO2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  57             EM: 9 
chch_t( 57,'C5H8      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  58
chch_t( 58,'ISO2      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  59 DD:25,WD:21,
chch_t( 59,'ISOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  60 DD:26,WD:22,
chch_t( 60,'ISON      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  61 DD:27,
chch_t( 61,'MACR      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  62
chch_t( 62,'MACRO2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  63 DD:28,WD:23,
chch_t( 63,'MACROOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  64 DD:29,
chch_t( 64,'MPAN      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  65 DD:30,WD:24,
chch_t( 65,'HACET     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  66 DD:31,WD:25,
chch_t( 66,'MGLY      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  67 DD:32,
chch_t( 67,'NALD      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  68 DD:33,WD:26,
chch_t( 68,'HCOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  69 DD:34,WD:27,
chch_t( 69,'MeCO3H    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  70 DD:35,WD:28,
chch_t( 70,'MeCO2H    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  71
chch_t( 71,'H2        ',  1,'TR        ','          ',  0,  0,  0),  & 
!  72 DD:36,WD:29,
chch_t( 72,'MeOH      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  73
chch_t( 73,'CO2       ',  1,'CT        ','          ',  0,  0,  0),  & 
!  74
chch_t( 74,'O2        ',  1,'CT        ','          ',  0,  0,  0),  & 
!  75
chch_t( 75,'N2        ',  1,'CT        ','          ',  0,  0,  0)   & 
  /)

TYPE(CHCH_T), PUBLIC :: chch_defs_strattrop_aer(1:87)=(/             &
!   1
chch_t(  1,'O(3P)     ',  1,'TR        ','Ox        ',  0,  0,  0),  & 
!   2
chch_t(  2,'O(1D)     ',  1,'SS        ','Ox        ',  0,  0,  0),  & 
!   3 DD: 1,WD: 1,
chch_t(  3,'O3        ',  1,'TR        ','Ox        ',  1,  1,  0),  & 
!   4
chch_t(  4,'N         ',  1,'TR        ','NOx       ',  0,  0,  0),  & 
!   5 DD: 2,
chch_t(  5,'NO        ',  1,'TR        ','NOx       ',  1,  0,  0),  & 
!   6 DD: 3,WD: 2,
chch_t(  6,'NO3       ',  1,'TR        ','NOx       ',  1,  1,  0),  & 
!   7 DD: 4,      EM: 1
chch_t(  7,'NO2       ',  1,'TR        ','NOx       ',  1,  0,  1),  & 
!   8 DD: 5,WD: 3,
chch_t(  8,'N2O5      ',  1,'TR        ','          ',  1,  1,  0),  & 
!   9 DD: 6,WD: 4,
chch_t(  9,'HO2NO2    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  10 DD: 7,WD: 5,
chch_t( 10,'HONO2     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  11 DD: 8,WD: 6,
chch_t( 11,'H2O2      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  12 EM: 2
chch_t( 12,'CH4       ',  1,'TR        ','          ',  0,  0,  2),  & 
!  13 DD: 9,      EM:  3
chch_t( 13,'CO        ',  1,'TR        ','          ',  1,  0,  3),  & 
!  14 DD:10,WD: 7,EM:  4
chch_t( 14,'HCHO      ',  1,'TR        ','          ',  1,  1,  4),  & 
!  15       WD: 8,
chch_t( 15,'MeOO      ',  1,'TR        ','          ',  0,  1,  0),  & 
!  16 DD:11,WD: 9,
chch_t( 16,'MeOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  17
chch_t( 17,'H         ',  1,'TR        ','HOx       ',  0,  0,  0),  & 
!  18
chch_t( 18,'H2O       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  19
chch_t( 19,'OH        ',  1,'TR        ','HOx       ',  0,  0,  0),  & 
!  20       WD:10,
chch_t( 20,'HO2       ',  1,'TR        ','HOx       ',  0,  1,  0),  & 
!  21
chch_t( 21,'Cl        ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  22
chch_t( 22,'Cl2O2     ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  23
chch_t( 23,'ClO       ',  1,'TR        ','Clx       ',  0,  0,  0),  & 
!  24
chch_t( 24,'OClO      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  25
chch_t( 25,'Br        ',  1,'TR        ','Brx       ',  0,  0,  0),  & 
!  26
chch_t( 26,'BrO       ',  1,'TR        ','Brx       ',  0,  0,  0),  & 
!  27
chch_t( 27,'BrCl      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  28       WD:11,
chch_t( 28,'BrONO2    ',  1,'TR        ','          ',  0,  1,  0),  & 
!  29
chch_t( 29,'N2O       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  30 DD:12,WD:12,
chch_t( 30,'HCl       ',  1,'TR        ','          ',  1,  1,  0),  & 
!  31 DD:13,WD:13,
chch_t( 31,'HOCl      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  32 DD:14,WD:14,
chch_t( 32,'HBr       ',  1,'TR        ','          ',  1,  1,  0),  & 
!  33 DD:15,WD:15,
chch_t( 33,'HOBr      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  34       WD:16,
chch_t( 34,'ClONO2    ',  1,'TR        ','          ',  0,  1,  0),  & 
!  35
chch_t( 35,'CFCl3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  36
chch_t( 36,'CF2Cl2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  37
chch_t( 37,'MeBr      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  38 DD:16,WD:17,
chch_t( 38,'HONO      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  39             EM: 5
chch_t( 39,'C2H6      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  40
chch_t( 40,'EtOO      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  41 DD:17,WD:18,
chch_t( 41,'EtOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  42 DD:18,      EM: 6
chch_t( 42,'MeCHO     ',  1,'TR        ','          ',  1,  0,  1),  & 
!  43
chch_t( 43,'MeCO3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  44 DD:19,
chch_t( 44,'PAN       ',  1,'TR        ','          ',  1,  0,  0),  & 
!  45             EM: 7
chch_t( 45,'C3H8      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  46
chch_t( 46,'n-PrOO    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  47
chch_t( 47,'i-PrOO    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  48 DD:20,WD:19,
chch_t( 48,'n-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  49 DD:21,WD:20,
chch_t( 49,'i-PrOOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  50 DD:22,
chch_t( 50,'EtCHO     ',  1,'TR        ','          ',  1,  0,  0),  & 
!  51
chch_t( 51,'EtCO3     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  52             EM: 8
chch_t( 52,'Me2CO     ',  1,'TR        ','          ',  0,  0,  1),  & 
!  53
chch_t( 53,'MeCOCH2OO ',  1,'TR        ','          ',  0,  0,  0),  & 
!  54 DD:23,WD:21,
chch_t( 54,'MeCOCH2OOH',  1,'TR        ','          ',  1,  1,  0),  & 
!  55 DD:24,
chch_t( 55,'PPAN      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  56
chch_t( 56,'MeONO2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  57             EM: 9 
chch_t( 57,'C5H8      ',  1,'TR        ','          ',  0,  0,  1),  & 
!  58
chch_t( 58,'ISO2      ',  1,'TR        ','          ',  0,  0,  0),  & 
!  59 DD:25,WD:22,
chch_t( 59,'ISOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  60 DD:26,WD:23,
chch_t( 60,'ISON      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  61 DD:27,
chch_t( 61,'MACR      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  62
chch_t( 62,'MACRO2    ',  1,'TR        ','          ',  0,  0,  0),  & 
!  63 DD:28,WD:24,
chch_t( 63,'MACROOH   ',  1,'TR        ','          ',  1,  1,  0),  & 
!  64 DD:29,      
chch_t( 64,'MPAN      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  65 DD:30,WD:25,
chch_t( 65,'HACET     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  66 DD:31,WD:26,
chch_t( 66,'MGLY      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  67 DD:32,
chch_t( 67,'NALD      ',  1,'TR        ','          ',  1,  0,  0),  & 
!  68 DD:33,WD:27,
chch_t( 68,'HCOOH     ',  1,'TR        ','          ',  1,  1,  0),  & 
!  69 DD:34,WD:28,
chch_t( 69,'MeCO3H    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  70 DD:35,WD:29,
chch_t( 70,'MeCO2H    ',  1,'TR        ','          ',  1,  1,  0),  & 
!  71
chch_t( 71,'H2        ',  1,'TR        ','          ',  0,  0,  0),  & 
!  72 DD:36,WD:30,
chch_t( 72,'MeOH      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  73
chch_t( 73,'CO2       ',  1,'CT        ','          ',  0,  0,  0),  & 
!  74
chch_t( 74,'O2        ',  1,'CT        ','          ',  0,  0,  0),  & 
!  75
chch_t( 75,'N2        ',  1,'CT        ','          ',  0,  0,  0),  & 
!  76
chch_t( 76,'DMS       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  77 DD:37,WD:31,EM:10
chch_t( 77,'SO2       ',  1,'TR        ','          ',  1,  1,  1),  & 
!  78
chch_t( 78,'H2SO4     ',  1,'TR        ','          ',  0,  0,  0),  & 
!  79
chch_t( 79,'MSA       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  80 DD:38,WD:32
chch_t( 80,'DMSO      ',  1,'TR        ','          ',  1,  1,  0),  & 
!  81 DD:39,WD:33,EM:11
chch_t( 81,'NH3       ',  1,'TR        ','          ',  1,  1,  1),  & 
!  82
chch_t( 82,'CS2       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  83
chch_t( 83,'COS       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  84
chch_t( 84,'H2S       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  85
chch_t( 85,'SO3       ',  1,'TR        ','          ',  0,  0,  0),  & 
!  86 DD:40,      EM: 12
chch_t( 86,'Monoterp  ',  1,'TR        ','          ',  1,  0,  1),  & 
!  87 DD:41,WD:34
chch_t( 87,'Sec_Org   ',  1,'TR        ','          ',  1,  1,  0)   & 
  /)

TYPE(RATB_T) :: ratb_defs_strattrop_chem(198)

! reactions found in either Trop or Strat but not both 
TYPE(RATB_T), PARAMETER :: ratb_defs_strattrop_aer(1:15)=(/              &
ratb_t('CS2       ','O(3P)     ','COS       ','SO2       ','CO        ', &
'          ',3.20E-11,  0.00,    650.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('CS2       ','OH        ','COS       ','SO2       ','          ', &
'          ',1.25E-16,  0.00,  -4550.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('DMS       ','OH        ','SO2       ','          ','          ', &
'          ',1.20E-11,  0.00,    260.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('DMS       ','OH        ','MSA       ','SO2       ','          ', &
'          ',3.04E-12,  0.00,   -350.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('DMS       ','NO3       ','SO2       ','          ','          ', &
'          ',1.90E-13,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('DMS       ','O(3P)     ','SO2       ','          ','          ', &
'          ',1.30E-11,  0.00,   -410.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('H2S       ','O(3P)     ','OH        ','SO2       ','          ', &
'          ',9.20E-12,  0.00,   1800.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('H2S       ','OH        ','SO2       ','H2O       ','          ', &
'          ',6.00E-12,  0.00,     75.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('COS       ','O(3P)     ','CO        ','SO2       ','          ', &
'          ',2.10E-11,  0.00,   2200.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('COS       ','OH        ','CO2       ','SO2       ','          ', &
'          ',1.10E-13,  0.00,   1200.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('SO2       ','O3        ','SO3       ','          ','          ', &
'          ',3.00E-12,  0.00,   7000.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('SO3       ','H2O       ','H2SO4     ','H2O       ','          ', &
'          ',8.50E-41,  0.00,  -6540.00, 0.000, 0.000, 0.000, 0.000) ,   &
ratb_t('Monoterp  ','OH        ','Sec_Org   ','          ','          ', &
'          ',1.20E-11,  0.00,   -444.00, 0.130, 0.000, 0.000, 0.000) ,   &
ratb_t('Monoterp  ','O3        ','Sec_Org   ','          ','          ', &
'          ',1.01E-15,  0.00,    732.00, 0.130, 0.000, 0.000, 0.000) ,   &
ratb_t('Monoterp  ','NO3       ','Sec_Org   ','          ','          ', &
'          ',1.19E-12,  0.00,   -925.00, 0.130, 0.000, 0.000, 0.000)     &
  /)

TYPE(RATH_T), ALLOCATABLE :: rath_defs_strattrop_chem(:)

TYPE(RATH_T) :: rath_defs_strattrop_psc(1:nhet_strattrop)=(/             &
rath_t('ClONO2    ','H2O       ','HOCl      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('ClONO2    ','HCl       ','Cl        ','Cl        ','HONO2     ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('HOCl      ','HCl       ','Cl        ','Cl        ','H2O       ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('N2O5      ','H2O       ','HONO2     ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000), &
rath_t('N2O5      ','HCl       ','Cl        ','NO2       ','HONO2     ', &
'          ', 0.000, 0.000, 0.000, 0.000)  &
  /)


! Aerosol chemistry: there are no gas phase products, the 'NULLx' products
!  identify the reactions in asad_hetero
TYPE(RATH_T) :: rath_defs_strattrop_aer(1:nhet_st_aer)=(/                &
!HSO3+H2O2(aq)
rath_t('SO2       ','H2O2      ','NULL0     ','          ','          ', & 
'          ', 0.000, 0.000, 0.000, 0.000),                               &
!HSO3+O3(aq)
rath_t('SO2       ','O3        ','NULL1     ','          ','          ', & 
'          ', 0.000, 0.000, 0.000, 0.000),                               &
!HSO3+O3(aq)
rath_t('SO2       ','O3        ','NULL2     ','          ','          ', & 
'          ', 0.000, 0.000, 0.000, 0.000)                                &
  /)

! used if no Aerosol Chemistry included
TYPE(RATH_T) :: rath_defs_strattrop_null(0)

! Tropospheric heterogenous reactions
TYPE(RATH_T) :: rath_defs_strattrop_trophet(1:nhet_st_tpht)=(/           &
! Heterogenous
rath_t('N2O5      ','          ','HONO2     ','          ','          ', & 
'          ', 2.000, 0.000, 0.000, 0.000),                               &
! Heterogenous
rath_t('HO2       ','          ','H2O2      ','          ','          ', & 
'          ', 0.500, 0.000, 0.000, 0.000)                                &
  /)


 
TYPE(RATJ_T), PARAMETER :: ratj_defs_strattrop01(1:28)=(/                &
ratj_t('EtOOH     ','PHOTON    ','MeCHO     ','HO2       ','OH        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('H2O2      ','PHOTON    ','OH        ','OH        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jh2o2     ') ,  &
ratj_t('HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhchoa    ') ,  &
ratj_t('HCHO      ','PHOTON    ','H2        ','CO        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhchob    ') ,  &
ratj_t('HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jpna67    ') ,  &
ratj_t('HONO2     ','PHOTON    ','OH        ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhono2    ') ,  &
ratj_t('MeCHO     ','PHOTON    ','MeOO      ','HO2       ','CO        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jaceta    ') ,  &
ratj_t('MeCHO     ','PHOTON    ','CH4       ','CO        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jacetb    ') ,  &
ratj_t('MeOOH     ','PHOTON    ','HO2       ','HCHO      ','OH        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('N2O5      ','PHOTON    ','NO3       ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jn2o5     ') ,  &
ratj_t('NO2       ','PHOTON    ','NO        ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jno2      ') ,  &
ratj_t('NO3       ','PHOTON    ','NO        ','O2        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jno3a     ') ,  &
ratj_t('NO3       ','PHOTON    ','NO2       ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jno3b     ') ,  &
ratj_t('O2        ','PHOTON    ','O(3P)     ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jo2       ') ,  &
ratj_t('O3        ','PHOTON    ','O2        ','O(1D)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jo3a      ') ,  &
ratj_t('O3        ','PHOTON    ','O2        ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jo3b      ') ,  &
ratj_t('PAN       ','PHOTON    ','MeCO3     ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jpan      ') ,  &
ratj_t('HONO      ','PHOTON    ','OH        ','NO        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhono     ') ,  &
ratj_t('EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jetcho    ') ,  &
ratj_t('Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jaceto    ') ,  &
ratj_t('n-PrOOH   ','PHOTON    ','EtCHO     ','HO2       ','OH        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('i-PrOOH   ','PHOTON    ','Me2CO     ','HO2       ','OH        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('MeCOCH2OOH','PHOTON    ','MeCO3     ','HCHO      ','OH        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('PPAN      ','PHOTON    ','EtCO3     ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jpan      ') ,  &
ratj_t('MeONO2    ','PHOTON    ','HO2       ','HCHO      ','NO2       ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmena     ') ,  &
ratj_t('ISOOH     ','PHOTON    ','OH        ','MACR      ','HCHO      ', &
     'HO2       ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmhp      ') ,  &
ratj_t('ISON      ','PHOTON    ','NO2       ','MACR      ','HCHO      ', &
     'HO2       ',    0.0,   0.0,   0.0,   0.0, 100.000,'jiprn     ') ,  &
ratj_t('MACR      ','PHOTON    ','MeCO3     ','HCHO      ','CO        ', &
     'HO2       ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmacr     ')    &
  /) 

TYPE(RATJ_T), PARAMETER :: ratj_defs_strattrop02(1:28)=(/                &
ratj_t('MPAN      ','PHOTON    ','MACRO2    ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jpan      ') ,  &
ratj_t('MACROOH   ','PHOTON    ','OH        ','HO2       ','OH        ', &
     'HO2       ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmacro    ') ,  &
ratj_t('MACROOH   ','PHOTON    ','HACET     ','CO        ','MGLY      ', &
     'HCHO      ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmacro    ') ,  &
ratj_t('HACET     ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhacet    ') ,  &
ratj_t('MGLY      ','PHOTON    ','MeCO3     ','CO        ','HO2       ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmkal     ') ,  &
ratj_t('NALD      ','PHOTON    ','HCHO      ','CO        ','NO2       ', &
     'HO2       ',    0.0,   0.0,   0.0,   0.0, 100.000,'jaceta    ') ,  &
ratj_t('MeCO3H    ','PHOTON    ','MeOO      ','OH        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmeco3h   ') ,  &
ratj_t('BrCl      ','PHOTON    ','Br        ','Cl        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jbrcl     ') ,  &
ratj_t('BrO       ','PHOTON    ','Br        ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jbro      ') ,  &
ratj_t('BrONO2    ','PHOTON    ','Br        ','NO3       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jbrna     ') ,  &
ratj_t('BrONO2    ','PHOTON    ','BrO       ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jbrnb     ') ,  &
ratj_t('O2        ','PHOTON    ','O(3P)     ','O(1D)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jo2b      ') ,  &
ratj_t('OClO      ','PHOTON    ','O(3P)     ','ClO       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'joclo     ') ,  &
ratj_t('NO        ','PHOTON    ','N         ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jno       ') ,  &
ratj_t('HOBr      ','PHOTON    ','OH        ','Br        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhobr     ') ,  &
ratj_t('N2O       ','PHOTON    ','N2        ','O(1D)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jn2o      ') ,  &
ratj_t('H2O       ','PHOTON    ','OH        ','H         ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jh2o      ') ,  &
ratj_t('ClONO2    ','PHOTON    ','Cl        ','NO3       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jclna     ') ,  &
ratj_t('ClONO2    ','PHOTON    ','ClO       ','NO2       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jclnb     ') ,  &
ratj_t('HCl       ','PHOTON    ','H         ','Cl        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhcl      ') ,  &
ratj_t('HOCl      ','PHOTON    ','OH        ','Cl        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jhocl     ') ,  &
ratj_t('Cl2O2     ','PHOTON    ','Cl        ','Cl        ','O2        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jcl2o2    ') ,  &
ratj_t('CFCl3     ','PHOTON    ','Cl        ','Cl        ','Cl        ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jcfcl3    ') ,  &
ratj_t('CF2Cl2    ','PHOTON    ','Cl        ','Cl        ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jcfc2     ') ,  &
ratj_t('MeBr      ','PHOTON    ','Br        ','H         ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jmebr     ') ,  &
ratj_t('CH4       ','PHOTON    ','MeOO      ','H         ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jch4      ') ,  &
ratj_t('CO2       ','PHOTON    ','CO        ','O(3P)     ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jco2      ') ,  &
ratj_t('HO2NO2    ','PHOTON    ','OH        ','NO3       ','          ', &
     '          ',    0.0,   0.0,   0.0,   0.0, 100.000,'jpna33    ')    &
  /)
  
TYPE(RATJ_T) :: ratj_defs_strattrop_chem(1:56)=(/                        &
          ratj_defs_strattrop01,ratj_defs_strattrop02 /)
 

TYPE(RATJ_T) :: ratj_defs_strattrop_aer(1:4) = (/                        &
ratj_t('CS2      ','PHOTON    ','COS        ','SO2       ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jcs2      '),      &
ratj_t('COS      ','PHOTON    ','CO         ','SO2       ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jcos      '),      &
ratj_t('H2SO4    ','PHOTON    ','SO3        ','OH        ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jh2so4    '),      &
ratj_t('SO3      ','PHOTON    ','SO2        ','O(3P)     ','          ', &
     '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jso3      ')       &
  /)
   
! Rates taken from  file "UKCA_Reaction_Rates - Termolecular Reactions.csv" 
!on 11:41 2012/05/01 

TYPE(RATT_T) :: ratt_defs_strattrop_chem(1:24)=(/                       &
! T001 JPL 2011  
ratt_t('O(3P)     ','O2        ','O3        ','m         ',     0.0,    &  
  6.00E-34, -2.50,     0.00,  0.00E+00,  0.00,     0.00, 0.000, 0.000), & 
! T002 JPL 2011 
ratt_t('O(3P)     ','NO        ','NO2       ','m         ',     0.6,    & 
  9.00E-32, -1.50,     0.00,  3.00E-11,  0.00,     0.00, 0.000, 0.000), &   
! T003 JPL 2011 
ratt_t('O(3P)     ','NO2       ','NO3       ','m         ',     0.6,    & 
  2.50E-31, -1.80,     0.00,  2.20E-11, -0.70,     0.00, 0.000, 0.000), &    
! T004 JPL 2011 
ratt_t('O(1D)     ','N2        ','N2O       ','m         ',     0.0,    & 
  2.80E-36, -0.90,     0.00,  0.00E+00,  0.00,     0.00, 0.000, 0.000), &  
! T005 JPL 2011 
ratt_t('BrO       ','NO2       ','BrONO2    ','m         ',     0.6,    & 
  5.20E-31, -3.20,     0.00,  6.90E-12,  0.00,     0.00, 0.000, 0.000), &   
! T006 JPL 2011
ratt_t('ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,    & 
  1.60E-32, -4.50,     0.00,  3.00E-12, -2.00,     0.00, 0.000, 0.000), &   
! T007 IUPAC 2005 
ratt_t('Cl2O2     ','m         ','ClO       ','ClO       ',     0.5,    & 
  3.70E-07,  0.00,  7690.00,  1.80E+14,  0.00,  7690.00, 0.000, 0.000), &   
! T008 JPL 2011 
ratt_t('ClO       ','NO2       ','ClONO2    ','m         ',     0.6,    & 
  1.80E-31, -3.40,     0.00,  1.50E-11,  0.00,     0.00, 0.000, 0.000), &  
! T009 JPL 2011 
ratt_t('H         ','O2        ','HO2       ','m         ',     0.6,    & 
  4.40E-32, -1.30,     0.00,  7.50E-11,  0.00,     0.00, 0.000, 0.000), &   
! T010 JPL 2011 see also asad_trimol.F90 
ratt_t('HO2       ','HO2       ','H2O2      ','O2        ',     0.0,    & 
  2.10E-33,  0.00,  -920.00,  0.00E+00,  0.00,     0.00, 0.000, 0.000), & 
! T011 JPL 2011 
ratt_t('HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,    & 
  2.00E-31, -3.40,     0.00,  2.90E-12,  0.00,     0.00, 0.000, 0.000), &  
! T012 IUPAC 2001
ratt_t('HO2NO2    ','m         ','HO2       ','NO2       ',     0.5,    & 
  4.10E-05,  0.00, 10650.00,  4.80E+15,  0.00, 11170.00, 0.000, 0.000), &  
! T013 JPL 2011  
ratt_t('OH        ','NO        ','HONO      ','m         ',     0.6,    & 
  7.00E-31, -2.60,     0.00,  3.60E-11, -0.10,     0.00, 0.000, 0.000), &   
! T014 JPL 2011 
ratt_t('OH        ','NO2       ','HONO2     ','m         ',     0.6,    & 
  1.80E-30, -3.00,     0.00,  2.80E-11,  0.00,     0.00, 0.000, 0.000), &    
! T015 JPL 2011
ratt_t('OH        ','OH        ','H2O2      ','m         ',     0.6,    & 
  6.90E-31, -1.00,     0.00,  2.60E-11,  0.00,     0.00, 0.000, 0.000), &   
! T016 MCMv3.2 
ratt_t('MeCO3     ','NO2       ','PAN       ','m         ',     0.3,    & 
  2.70E-28, -7.10,     0.00,  1.20E-11, -0.90,     0.00, 0.000, 0.000), &   
! T017 MCMv3.2 
ratt_t('PAN       ','m         ','MeCO3     ','NO2       ',     0.3,    & 
  4.90E-03,  0.00, 12100.00,  5.40E+16,  0.00, 13830.00, 0.000, 0.000), &   
! T018 MCMv3.2
ratt_t('EtCO3     ','NO2       ','PPAN      ','m         ',     0.3,    &  
  2.70E-28, -7.10,     0.00,  1.20E-11, -0.90,     0.00, 0.000, 0.000), &  
! T019 MCMv3.2 
ratt_t('PPAN      ','m         ','EtCO3     ','NO2       ',     0.3,    & 
  4.90E-03,  0.00, 12100.00,  5.40E+16,  0.00, 13830.00, 0.000, 0.000), &   
! T020 MVMv3.2
ratt_t('MACRO2    ','NO2       ','MPAN      ','m         ',     0.3,    &  
  2.70E-28, -7.10,     0.00,  1.20E-11, -0.90,     0.00, 0.000, 0.000), &   
! T021 MCMv3.2
ratt_t('MPAN      ','m         ','MACRO2    ','NO2       ',     0.3,    & 
  4.90E-03,  0.00, 12100.00,  5.40E+16,  0.00, 13830.00, 0.000, 0.000), &   
! T022 IUPAC 2002 
ratt_t('NO2       ','NO3       ','N2O5      ','m         ',     0.3,    &  
  3.60E-30, -4.10,     0.00,  1.90E-12,  0.20,     0.00, 0.000, 0.000), &    
! T023 IUPAC 2002
ratt_t('N2O5      ','m         ','NO2       ','NO3       ',     0.3,    & 
  1.30E-03, -3.50, 11000.00,  9.70E+14,  0.10, 11080.00, 0.000, 0.000), &    
! T024 IUPAC 2001 
ratt_t('NO        ','NO        ','NO2       ','NO2       ',     0.0,    & 
  3.30E-39,  0.00,  -530.00,  0.00E+00,  0.00,     0.00, 0.000, 0.000)  &   
  /) 
 
!---------------------------------------------------------------------- 
! NOTES: CheST Termolecular Reactions 
!---------------------------------------------------------------------- 
! T001 O(3P)+O2 -> O3 m JPL 2011   
! T001 IUPAC 2002 recommend k = 5.623E-34*(T/300)^-2.6 (based on 
! T001 weighted mean of kO2+kN2) 
!---------------------------------------------------------------------- 
! T002 O(3P)+NO -> NO2 m JPL 2011   
! T002 IUPAC 2002 recommend k0 = 1.0E-3*(T/300)-1.6*[N2] 
!---------------------------------------------------------------------- 
! T003 O(3P)+NO2 -> NO3 m JPL 2011   
! T003 IUPAC 2002 recommend k0 = 1.3E-31*(T/300)^-1.5 kinf = 
! T003 2.3E-11*(T/300)^0.24 Fc = 0.6 
!---------------------------------------------------------------------- 
! T004 O(1D)+N2 -> N2O m JPL 2011   
! T004 IUPAC 2002 k0=2.8E-36*[N2] 
!---------------------------------------------------------------------- 
! T005 BrO+NO2 -> BrONO2 m JPL 2011   
! T005 IUPAC Fc = 0.55 k0 = 4.2E-31*exp(T/300)^-2.4*[N2] kinf=2.7E-11 
!---------------------------------------------------------------------- 
! T006 ClO+ClO -> Cl2O2 m JPL 2011   
! T006 IUPAC recommend k0 = 2.0E-32*(T/300)^-4*[N2] kinf = 1.0E-11. 
! T006 Fc=0.45 In general this is a problematic rate 
!---------------------------------------------------------------------- 
! T007 Cl2O2+m -> ClO ClO IUPAC 2005   
! T007 No JPL data 
!---------------------------------------------------------------------- 
! T008 ClO+NO2 -> ClONO2 m JPL 2011   
! T008 IUPAC recommend k0 = 1.6E-31*(T/300)^-3.4*[N2] kinf = 7.0E-11 
! T008 Fc=0.4 
!---------------------------------------------------------------------- 
! T009 H+O2 -> HO2 m JPL 2011   
! T009 IUPAC 2009 recommend k0 = 4.3E-32*(T/300)^-1.2*[N2] kinf = 
! T009 9.6E-11 and FC(ent) = 0.5 
!---------------------------------------------------------------------- 
! T010 HO2+HO2 -> H2O2 O2 JPL 2011 see also asad_trimol.F90  
! T010 IUPAC (2001) k = 1.9E-33*exp(980/T)*[N2]. Note that this reaction 
! T010 is special and in the presence of H2O the rate constant needs to 
! T010 be adjusted by: {1 + 1.4E-21*[H2O]*exp(2200/T)} (same H2O 
! T010 expression for JPL). JPL also include the formation of HO2-H2O as 
! T010 a separate species. 
!---------------------------------------------------------------------- 
! T011 HO2+NO2 -> HO2NO2 m JPL 2011   
! T011 IUPAC 2009 k0 = 1.4E-31*(T/300)^-3.1*[N2] kinf = 4.0E-12 Fc = 0.4 
!---------------------------------------------------------------------- 
! T012 HO2NO2+m -> HO2 NO2 IUPAC 2001   
! T012 No JPL data.. NOTE should multiply the expression by [N2] NOT [M] 
!---------------------------------------------------------------------- 
! T013 OH+NO -> HONO m JPL 2011   
! T013 IUPAC 2002 recommend k0 = 7.40E-31*(T/300)^-2.4 kinf = 
! T013 3.3E-11*(T/300)^-0.3 Fc = 0.81 
!---------------------------------------------------------------------- 
! T014 OH+NO2 -> HONO2 m JPL 2011   
! T014 IUPCA 2009 recommend k0 = 3.3E-30*(T/300)^-3.0 kinf = 6.0E-11 Fc 
! T014 = 0.4 
!---------------------------------------------------------------------- 
! T015 OH+OH -> H2O2 m JPL 2011   
! T015 IUPAC 2009 recommend k0 = 6.9E-31*(T/300)^-0.8 kinf = 
! T015 3.9E-11*(T/300)^-0.47 
!---------------------------------------------------------------------- 
! T016 MeCO3+NO2 -> PAN m MCMv3.2   
! T016 Based on IUPAC 2003. JPL recommend k0 = 9.7E-29*(T/300)^-5.6 kinf 
! T016 = 9.3E-12*(T/300)^-1.5 Fc = 0.6 
!---------------------------------------------------------------------- 
! T017 PAN+m -> MeCO3 NO2 MCMv3.2   
! T017 IUPAC 2003 
!---------------------------------------------------------------------- 
! T018 EtCO3+NO2 -> PPAN m MCMv3.2   
! T018 JPL recommend k0 = 9.0E-28*(T/300)^-8.9 kinf = 
! T018 7.7E-12*(T/300)^-0.2 Fc = 0.6. IUPAC 2003 k= 
! T018 1.70E-3*exp(-11280/T) kinf = 8.3E16*exp(-13940/T) Fc = 0.36 
!---------------------------------------------------------------------- 
! T019 PPAN+m -> EtCO3 NO2 MCMv3.2   
! T019 IUPAC 2003 k0 = 1.70E-03*exp(-11280/T) kinf = 
! T019 8.30E+16*exp(-13940/T) Fc = 0.36 
!---------------------------------------------------------------------- 
! T020 MACRO2+NO2 -> MPAN m MVMv3.2   
! T020 - 
!---------------------------------------------------------------------- 
! T021 MPAN+m -> MACRO2 NO2 MCMv3.2   
! T021 IUPAC recommend k = 1.6E16*exp(-13500/T). No JPL data 
!---------------------------------------------------------------------- 
! T022 NO2+NO3 -> N2O5 m IUPAC 2002   
! T022 JPL recommend k0 = 2.0E-30*(T/300)^-4.4 kinf = 
! T022 1.4E-12*(T/300)^-0.7 Fc = 0.6 
!---------------------------------------------------------------------- 
! T023 N2O5+m -> NO2 NO3 IUPAC 2002   
! T023 No JPL data. NOTE there is code in asad_trimol.F90 which looks 
! T023 for this reaction and modifies the rate by an extra factor. I can 
! T023 NOT find why it does this in the literature. 
!---------------------------------------------------------------------- 
! T024 NO+NO -> NO2 NO2 IUPAC 2001   
! T024 No JPL data.. 
!----------------------------------------------------------------------   
   
TYPE(RATT_T) :: ratt_defs_strattrop_aer(1)=(/                           &
ratt_t('SO2       ','OH        ','SO3       ','HO2       ', 0.6,        &
    3.00E-31, -3.30,     0.0, 1.50E-12,  0.00,       0.0, 0.000, 0.000) &
  /)
 
 
REAL :: depvel_defs_strattrop_chem(1:6,1:5,1:ndry_strattrop)


REAL :: depvel_defs_strattrop_aer(1:6,1:5,1:ndry_st_aer)=RESHAPE((/     &
! SO2 (added as in ukca_chem_aer.F90)
  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & !     37.1 
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     37.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     37.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     37.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     37.5
! DMSO    
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !     38.1 
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     38.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     38.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     38.4 
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     38.5
! NH3     
  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & !     39.1 
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     39.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     39.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  & !     39.4 
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     39.5
! monoterp
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     40.1 
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     40.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     40.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     40.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     40.5
! SEC_ORG 
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     41.1 
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     41.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     41.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     41.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10   & !     41.5
  /),(/  6,  5, ndry_st_aer/) )



! The following formula is used to calculate the effective Henry's Law coef,
! which takes the affects of dissociation and complex formation on a species'
! solubility (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 & 2 above give the data for this gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase the above term is multiplied 
! by another factor of 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 & 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
! The data in columns 5 and 6 give the data for a second dissociation,
! e.g for SO2, HSO3^{-}, and SO3^{2-}
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]


REAL :: henry_defs_strattrop_chem(1:6,1:nwet_strattrop)=RESHAPE((/      &
!    1  NO3
 0.2000E+01, 0.2000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    2  N2O5
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!    3  HO2NO2
 0.1300E+05, 0.6900E+04, 0.1000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!    4  HONO2
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!    5  H2O2
 0.8300E+05, 0.7400E+04, 0.2400E-11,-0.3730E+04, 0.0000E+00, 0.0000E+00,&   
!    6  HCHO
 0.3300E+04, 0.6500E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    7  MeOO
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    8  MeOOH
 0.3100E+03, 0.5000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    9  HO2
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   10  BrONO2
 0.2100E+06, 0.8700E+04, 0.1570E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   11  HCl
 0.1900E+02, 0.6000E+03, 0.1000E+05, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   12  HOCl
 0.9300E+03, 0.5900E+04, 0.3200E-07, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   13  HBr
 0.1300E+01, 0.1020E+05, 0.1000E+10, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   14  HOBr
 0.6100E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   15  ClONO2
 0.2100E+06, 0.8700E+04, 0.1570E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   16  HONO
 0.5000E+02, 0.4900E+04, 0.5600E-03,-0.1260E+04, 0.0000E+00, 0.0000E+00,&   
!   17  EtOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   18  n-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   19  i-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   20  MeCOCH2OOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   21  ISOOH
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   22  ISON
 0.3000E+04, 0.7400E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   23  MACROOH
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   24  HACET
 0.1400E+03, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   25  MGLY
 0.3500E+04, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   26  HCOOH
 0.6900E+04, 0.5600E+04, 0.1800E-03,-0.1510E+04, 0.0000E+00, 0.0000E+00,&   
!   27  MeCO3H
 0.7500E+03, 0.5300E+04, 0.6300E-08, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   28  MeCO2H
 0.4700E+04, 0.6000E+04, 0.1800E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   29  MeOH
 0.2300E+03, 0.4900E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00 &   
  /),(/  6, nwet_strattrop/) )


REAL :: henry_defs_strattrop_aer(1:6,1:nwet_st_aer)=RESHAPE((/          &
!    1  O3
 0.1130E-01, 0.2300E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    2  NO3
 0.2000E+01, 0.2000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    3  N2O5
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    4  HO2NO2
 0.1300E+05, 0.6900E+04, 0.1000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!    5  HONO2
 0.2100E+06, 0.8700E+04, 0.2000E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    6  H2O2
 0.8300E+05, 0.7400E+04, 0.2400E-11,-0.3730E+04, 0.0000E+00, 0.0000E+00,&   
!    7  HCHO
 0.3300E+04, 0.6500E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    8  MeOO
 0.2000E+04, 0.6600E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!    9  MeOOH
 0.3100E+03, 0.5000E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   10  HO2
 0.4000E+04, 0.5900E+04, 0.2000E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   11  BrONO2
 0.2100E+06, 0.8700E+04, 0.1570E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   12  HCl
 0.1900E+02, 0.6000E+03, 0.1000E+05, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   13  HOCl
 0.9300E+03, 0.5900E+04, 0.3200E-07, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   14  HBr
 0.1300E+01, 0.1020E+05, 0.1000E+10, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   15  HOBr
 0.6100E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   16  ClONO2
 0.2100E+06, 0.8700E+04, 0.1570E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   17  HONO
 0.5000E+02, 0.4900E+04, 0.5600E-03,-0.1260E+04, 0.0000E+00, 0.0000E+00,&  
!   18  EtOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   19  n-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   20  i-PrOOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   21  MeCOCH2OOH
 0.3400E+03, 0.5700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   22  ISOOH
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   23  ISON
 0.3000E+04, 0.7400E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&  
!   24  MACROOH
 0.1700E+07, 0.9700E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   25  HACET
 0.1400E+03, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   26  MGLY
 0.3500E+04, 0.7200E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   27  HCOOH
 0.6900E+04, 0.5600E+04, 0.1800E-03,-0.1510E+04, 0.0000E+00, 0.0000E+00,&  
!   28  MeCO3H
 0.7500E+03, 0.5300E+04, 0.6300E-08, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   29  MeCO2H
 0.4700E+04, 0.6000E+04, 0.1800E-04, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   30  MeOH
 0.2300E+03, 0.4900E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   31  SO2
 0.1230E+01, 0.3020E+04, 0.1230E-01, 0.2010E+04, 0.6000E-07, 0.1120E+04,&   
!   32  DMSO
 0.5000E+05, 0.6425E+04, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   33  NH3 (H*)
 0.1000E+07, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,&   
!   34  SEC_ORG
 0.1000E+06, 0.1200E+02, 0.0000E-00, 0.0000E+00, 0.0000E+00, 0.0000E+00 &   
  /),(/  6,  nwet_st_aer/) )
  
INTERFACE ukca_init_strattrop
  MODULE PROCEDURE ukca_init_strattrop
END INTERFACE ukca_init_strattrop
PUBLIC ukca_init_strattrop

CONTAINS

SUBROUTINE UKCA_INIT_STRATTROP()

! To initialise chemistry files from specified species and rate arrays
! using the defined logicals. Composes chemistry definition arrays from
! components using logical definitions from RUN_UKCA namelist held in 
! UKCA module ukca_option_mod

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

  INTEGER            :: ierr
  CHARACTER (LEN=72) :: cmessage

TYPE(RATB_T) :: ratb_defs_strattrop01(45)
TYPE(RATB_T) :: ratb_defs_strattrop02(45)
TYPE(RATB_T) :: ratb_defs_strattrop03(45)
TYPE(RATB_T) :: ratb_defs_strattrop04(45)
TYPE(RATB_T) :: ratb_defs_strattrop05(18)

REAL :: depvel_defs_strattrop01(360)
REAL :: depvel_defs_strattrop02(360)
REAL :: depvel_defs_strattrop03(360)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_INIT_STRATTROP',zhook_in,zhook_handle)


! allocate final arrays
ALLOCATE(chch_defs_strattrop(jpspec))
ALLOCATE(ratb_defs_strattrop(jpbk))
ALLOCATE(ratj_defs_strattrop(jppj))
ALLOCATE(ratt_defs_strattrop(jptk))
ALLOCATE(rath_defs_strattrop(jphk))
ALLOCATE(depvel_defs_strattrop(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_strattrop(6,jpdw))


! have to define some arrays here due to problems with some compilers

! Rates taken from  file "UKCA_Reaction_Rates - Bimolecular Reactions.csv" 
! on 12:08 2012/05/01 

ratb_defs_strattrop01 = (/   & 
! B001 JPL2011
ratb_t('Br        ','Cl2O2     ','BrCl      ','Cl        ','O2        ',& 
'          ',  5.90E-12,  0.00,    170.00, 0.000, 0.000, 0.000, 0.000), &    
! B002 JPL2011
ratb_t('Br        ','HCHO      ','HBr       ','CO        ','HO2       ',& 
'          ',  1.70E-11,  0.00,    800.00, 0.000, 0.000, 0.000, 0.000), &    
! B003 JPL2011
ratb_t('Br        ','HO2       ','HBr       ','O2        ','          ',& 
'          ',  4.80E-12,  0.00,    310.00, 0.000, 0.000, 0.000, 0.000), &    
! B004 JPL2011
ratb_t('Br        ','O3        ','BrO       ','O2        ','          ',& 
'          ',  1.60E-11,  0.00,    780.00, 0.000, 0.000, 0.000, 0.000), &    
! B005 JPL2011
ratb_t('Br        ','OClO      ','BrO       ','ClO       ','          ',& 
'          ',  2.60E-11,  0.00,   1300.00, 0.000, 0.000, 0.000, 0.000), &    
! B006 JPL2011
ratb_t('BrO       ','BrO       ','Br        ','Br        ','O2        ',& 
'          ',  2.40E-12,  0.00,    -40.00, 0.000, 0.000, 0.000, 0.000), &    
! B007 JPL2011
ratb_t('BrO       ','ClO       ','Br        ','Cl        ','O2        ',& 
'          ',  2.30E-12,  0.00,   -260.00, 0.000, 0.000, 0.000, 0.000), &   
! B008 JPL2011
ratb_t('BrO       ','ClO       ','Br        ','OClO      ','          ',& 
'          ',  9.50E-13,  0.00,   -550.00, 0.000, 0.000, 0.000, 0.000), &    
! B009 JPL2011
ratb_t('BrO       ','ClO       ','BrCl      ','O2        ','          ',& 
'          ',  4.10E-13,  0.00,   -290.00, 0.000, 0.000, 0.000, 0.000), &    
! B010 JPL2011
ratb_t('BrO       ','HO2       ','HBr       ','O3        ','          ',& 
'          ',  0.00E+00,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B011 JPL2011 
ratb_t('BrO       ','HO2       ','HOBr      ','O2        ','          ',& 
'          ',  4.50E-12,  0.00,   -460.00, 0.000, 0.000, 0.000, 0.000), &    
! B012 JPL2011
ratb_t('BrO       ','NO        ','Br        ','NO2       ','          ',& 
'          ',  8.80E-12,  0.00,   -260.00, 0.000, 0.000, 0.000, 0.000), &    
! B013 JPL2011 
ratb_t('BrO       ','OH        ','Br        ','HO2       ','          ',& 
'          ',  1.70E-11,  0.00,   -250.00, 0.000, 0.000, 0.000, 0.000), &  
! B014 JPL2011
ratb_t('CF2Cl2    ','O(1D)     ','Cl        ','ClO       ','          ',& 
'          ',  1.40E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B015 JPL2011 
ratb_t('CFCl3     ','O(1D)     ','Cl        ','Cl        ','ClO       ',& 
'          ',  2.30E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B016 JPL2011 
ratb_t('Cl        ','CH4       ','HCl       ','MeOO      ','          ',& 
'          ',  7.30E-12,  0.00,   1280.00, 0.000, 0.000, 0.000, 0.000), &    
! B017 IUPAC2006 
ratb_t('Cl        ','Cl2O2     ','Cl        ','Cl        ','Cl        ',& 
'          ',  7.60E-11,  0.00,    -65.00, 0.000, 0.000, 0.000, 0.000), &    
! B018 JPL2011
ratb_t('Cl        ','ClONO2    ','Cl        ','Cl        ','NO3       ',& 
'          ',  6.50E-12,  0.00,   -135.00, 0.000, 0.000, 0.000, 0.000), &    
! B019 JPL2011
ratb_t('Cl        ','H2        ','HCl       ','H         ','          ',& 
'          ',  3.05E-11,  0.00,   2270.00, 0.000, 0.000, 0.000, 0.000), &   
! B020 JPL2011 
ratb_t('Cl        ','H2O2      ','HCl       ','HO2       ','          ',& 
'          ',  1.10E-11,  0.00,    980.00, 0.000, 0.000, 0.000, 0.000), &    
! B021 JPL2011 
ratb_t('Cl        ','HCHO      ','HCl       ','CO        ','HO2       ',& 
'          ',  8.10E-11,  0.00,     30.00, 0.000, 0.000, 0.000, 0.000), &   
! B022 JPL2011 
ratb_t('Cl        ','HO2       ','ClO       ','OH        ','          ',& 
'          ',  3.65E-11,  0.00,    375.00, 0.000, 0.000, 0.000, 0.000), &   
! B023 JPL2011
ratb_t('Cl        ','HO2       ','HCl       ','O2        ','          ',& 
'          ',  1.40E-11,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000), &   
! B024 JPL2011
ratb_t('Cl        ','HOCl      ','Cl        ','Cl        ','OH        ',& 
'          ',  3.40E-12,  0.00,    130.00, 0.000, 0.000, 0.000, 0.000), &   
! B025 JPL2011
ratb_t('Cl        ','MeOOH     ','HCl       ','MeOO      ','          ',& 
'          ',  5.70E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B026 JPL2011
ratb_t('Cl        ','NO3       ','ClO       ','NO2       ','          ',& 
'          ',  2.40E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B027 JPL2011
ratb_t('Cl        ','O3        ','ClO       ','O2        ','          ',& 
'          ',  2.30E-11,  0.00,    200.00, 0.000, 0.000, 0.000, 0.000), &    
! B028 JPL2011
ratb_t('Cl        ','OClO      ','ClO       ','ClO       ','          ',& 
'          ',  3.40E-11,  0.00,   -160.00, 0.000, 0.000, 0.000, 0.000), &    
! B029 JPL2011 
ratb_t('ClO       ','ClO       ','Cl        ','Cl        ','O2        ',& 
'          ',  1.00E-12,  0.00,   1590.00, 0.000, 0.000, 0.000, 0.000), &    
! B030 JPL2011
ratb_t('ClO       ','ClO       ','Cl        ','Cl        ','O2        ',& 
'          ',  3.00E-11,  0.00,   2450.00, 0.000, 0.000, 0.000, 0.000), &   
! B031 JPL2011
ratb_t('ClO       ','ClO       ','Cl        ','OClO      ','          ',& 
'          ',  3.50E-13,  0.00,   1370.00, 0.000, 0.000, 0.000, 0.000), &    
! B032 JPL2011 
ratb_t('ClO       ','HO2       ','HCl       ','O3        ','          ',& 
'          ',  0.00E+00,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B033 JPL2011
ratb_t('ClO       ','HO2       ','HOCl      ','O2        ','          ',& 
'          ',  2.60E-12,  0.00,   -290.00, 0.000, 0.000, 0.000, 0.000), &    
! B034 JPL2011
ratb_t('ClO       ','MeOO      ','Cl        ','HCHO      ','HO2       ',& 
'          ',  3.30E-12,  0.00,    115.00, 0.000, 0.000, 0.000, 0.000), &   
! B035 JPL2011
ratb_t('ClO       ','NO        ','Cl        ','NO2       ','          ',& 
'          ',  6.40E-12,  0.00,   -290.00, 0.000, 0.000, 0.000, 0.000), &    
! B036 JPL2011
ratb_t('ClO       ','NO3       ','Cl        ','O2        ','NO2       ',& 
'          ',  4.60E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B037 JPL2011 
ratb_t('ClO       ','NO3       ','OClO      ','NO2       ','          ',& 
'          ',  0.00E+00,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B038 IUPAC2002
ratb_t('EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ',& 
'          ',  6.70E-12,  0.00,   -340.00, 0.000, 0.000, 0.000, 0.000), &    
! B039 MCMv3.2
ratb_t('EtCO3     ','NO3       ','EtOO      ','CO2       ','NO2       ',& 
'          ',  4.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B040 IUPAC2002
ratb_t('EtOO      ','MeCO3     ','MeCHO     ','HO2       ','MeOO      ',& 
'          ',  4.40E-13,  0.00,  -1070.00, 0.000, 0.000, 0.000, 0.000), &    
! B041 IUPAC2005
ratb_t('EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',&  
'          ',  2.55E-12,  0.00,   -380.00, 0.000, 0.000, 0.000, 0.000), &   
! B042 IUPAC2008 
ratb_t('EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ',& 
'          ',  2.30E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B043 JPL2011
ratb_t('H         ','HO2       ','H2        ','O2        ','          ',& 
'          ',  6.90E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B044 JPL2011
ratb_t('H         ','HO2       ','O(3P)     ','H2O       ','          ',& 
'          ',  1.62E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B045 JPL2011
ratb_t('H         ','HO2       ','OH        ','OH        ','          ',& 
'          ',  7.20E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000)  &  
  /) 
  
ratb_defs_strattrop02 = (/   & 
! B046 JPL2011
ratb_t('H         ','NO2       ','OH        ','NO        ','          ',& 
'          ',  4.00E-10,  0.00,    340.00, 0.000, 0.000, 0.000, 0.000), &   
! B047 JPL2011
ratb_t('H         ','O3        ','OH        ','O2        ','          ',& 
'          ',  1.40E-10,  0.00,    470.00, 0.000, 0.000, 0.000, 0.000), &   
! B048 MCMv3.2* 
ratb_t('HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',& 
'          ',  4.40E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000), &   
! B049 MCMv3.2
ratb_t('HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',& 
'          ',  7.80E-14,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000), &    
! B050 IUPAC2011
ratb_t('HO2       ','EtOO      ','EtOOH     ','          ','          ',& 
'          ',  6.40E-13,  0.00,   -710.00, 0.000, 0.000, 0.000, 0.000), &    
! B051 JPL2011 see also asad_bimol  
ratb_t('HO2       ','HO2       ','H2O2      ','          ','          ',& 
'          ',  3.00E-13,  0.00,   -460.00, 0.000, 0.000, 0.000, 0.000), & 
! B052 Poschl00
ratb_t('HO2       ','ISO2      ','ISOOH     ','          ','          ',& 
'          ',  2.05E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000), &   
! B053 Poschl00
ratb_t('HO2       ','MACRO2    ','MACROOH   ','          ','          ',& 
'          ',  1.82E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000), &  
! B054 IUPAC2009
ratb_t('HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',&  
'          ',  7.80E-14,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000), &   
! B055 IUPAC2009
ratb_t('HO2       ','MeCO3     ','MeCO3H    ','          ','          ',& 
'          ',  2.13E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000), &    
! B056 IUPAC2009
ratb_t('HO2       ','MeCO3     ','OH        ','MeOO      ','          ',& 
'          ',  2.29E-13,  0.00,   -980.00, 0.000, 0.000, 0.000, 0.000), &   
! B057 IUPAC2009
ratb_t('HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',& 
'          ',  9.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B058 IUPAC2009 see also asad_bimol
ratb_t('HO2       ','MeOO      ','HCHO      ','          ','          ',& 
'          ',  3.80E-13,  0.00,   -780.00, 0.000, 0.000, 0.000, 0.000), & 
! B059 IUPAC2009 see also asad_bimol  
ratb_t('HO2       ','MeOO      ','MeOOH     ','          ','          ',&  
'          ',  3.80E-13,  0.00,   -780.00, 0.000, 0.000, 0.000, 0.000), &  
! B060 JPL2011
ratb_t('HO2       ','NO        ','OH        ','NO2       ','          ',& 
'          ',  3.30E-12,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000), &    
! B061 JPL2011
ratb_t('HO2       ','NO3       ','OH        ','NO2       ','          ',& 
'          ',  3.50E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B062 IUPAC2001
ratb_t('HO2       ','O3        ','OH        ','O2        ','          ',& 
'          ',  2.03E-16,  4.57,   -693.00, 0.000, 0.000, 0.000, 0.000), &   
! B063 MCMv3.2
ratb_t('HO2       ','i-PrOO    ','i-PrOOH   ','          ','          ',& 
'          ',  1.51E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000), &   
! B064 MCMv3.2
ratb_t('HO2       ','n-PrOO    ','n-PrOOH   ','          ','          ',& 
'          ',  1.51E-13,  0.00,  -1300.00, 0.000, 0.000, 0.000, 0.000), &    
! B065 Poschl00
ratb_t('ISO2      ','ISO2      ','MACR      ','MACR      ','HCHO      ',& 
'HO2       ',  2.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B066 Poschl00
ratb_t('MACRO2    ','MACRO2    ','HACET     ','MGLY      ','HCHO      ',& 
'CO        ',  1.00E-12,  0.00,      0.00, 2.000, 2.000, 1.000, 1.000), &   
! B067 Poschl00
ratb_t('MACRO2    ','MACRO2    ','HO2       ','          ','          ',& 
'          ',  1.00E-12,  0.00,      0.00, 2.000, 0.000, 0.000, 0.000), &   
! B068 JPL2011
ratb_t('MeBr      ','Cl        ','Br        ','HCl       ','          ',& 
'          ',  1.40E-11,  0.00,   1030.00, 0.000, 0.000, 0.000, 0.000), &    
! B069 JPL2011
ratb_t('MeBr      ','O(1D)     ','Br        ','OH        ','          ',& 
'          ',  1.80E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B070 JPL2011
ratb_t('MeBr      ','OH        ','Br        ','H2O       ','          ',& 
'          ',  2.35E-12,  0.00,   1300.00, 0.000, 0.000, 0.000, 0.000), &   
! B071 IUPAC2002
ratb_t('MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',& 
'          ',  7.50E-12,  0.00,   -290.00, 0.000, 0.000, 0.000, 0.000), &  
! B072 IUPAC2008
ratb_t('MeCO3     ','NO3       ','MeOO      ','CO2       ','NO2       ',& 
'          ',  4.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B073 MCMv3.2
ratb_t('MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',& 
'          ',  2.70E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &   
! B074 MCMv3.2 
ratb_t('MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ',& 
'          ',  2.30E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B075 IUPAC2002
ratb_t('MeOO      ','MeCO3     ','HO2       ','HCHO      ','MeOO      ',& 
'          ',  1.80E-12,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000), &   
! B076 IUPAC2002
ratb_t('MeOO      ','MeCO3     ','MeCO2H    ','HCHO      ','          ',& 
'          ',  2.00E-13,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000), &    
! B077 IUPAC2002 see also asad_bimol
ratb_t('MeOO      ','MeOO      ','HO2       ','HO2       ','HCHO      ',&   
'HCHO      ',  1.03E-13,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000), & 
! B078 IUPAC2002 see also asad_bimol
ratb_t('MeOO      ','MeOO      ','MeOH      ','HCHO      ','          ',& 
'          ',  1.03E-13,  0.00,   -365.00, 0.000, 0.000, 0.000, 0.000), & 
! B079 IUPAC2005
ratb_t('MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',& 
'          ',  2.30E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), & 
! B080 IUPAC2005
ratb_t('MeOO      ','NO        ','MeONO2    ','          ','          ',& 
'          ',  2.30E-15,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &    
! B081 MCMv3.2
ratb_t('MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',&  
'          ',  1.20E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B082 JPL2011
ratb_t('N         ','NO        ','N2        ','O(3P)     ','          ',& 
'          ',  2.10E-11,  0.00,   -100.00, 0.000, 0.000, 0.000, 0.000), &   
! B083 JPL2011
ratb_t('N         ','NO2       ','N2O       ','O(3P)     ','          ',& 
'          ',  5.80E-12,  0.00,   -220.00, 0.000, 0.000, 0.000, 0.000), &  
! B084 JPL2011 
ratb_t('N         ','O2        ','NO        ','O(3P)     ','          ',& 
'          ',  1.50E-11,  0.00,   3600.00, 0.000, 0.000, 0.000, 0.000), &   
! B085 ????
ratb_t('N2O5      ','H2O       ','HONO2     ','HONO2     ','          ',& 
'          ',  2.50E-22,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B086 Poschl00 
ratb_t('NO        ','ISO2      ','ISON      ','          ','          ',& 
'          ',  1.12E-13,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &    
! B087 Poschl00
ratb_t('NO        ','ISO2      ','NO2       ','MACR      ','HCHO      ',& 
'HO2       ',  2.43E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &    
! B088 Poschl00
ratb_t('NO        ','MACRO2    ','MGLY      ','HCHO      ','HO2       ',& 
'          ',  1.27E-12,  0.00,   -360.00, 1.000, 1.500, 1.500, 0.000), &    
! B089 Poschl00 
ratb_t('NO        ','MACRO2    ','NO2       ','MeCO3     ','HACET     ',&  
'CO        ',  1.27E-12,  0.00,   -360.00, 2.000, 0.500, 0.500, 0.500), &   
! B090 JPL2011
ratb_t('NO        ','NO3       ','NO2       ','NO2       ','          ',& 
'          ',  1.50E-11,  0.00,   -170.00, 0.000, 0.000, 0.000, 0.000)  & 
  /) 
  
ratb_defs_strattrop03 = (/   & 
! B091 JPL2011
ratb_t('NO        ','O3        ','NO2       ','          ','          ',& 
'          ',  3.00E-12,  0.00,   1500.00, 0.000, 0.000, 0.000, 0.000), &    
! B092 JPL2011
ratb_t('NO2       ','NO3       ','NO        ','NO2       ','O2        ',& 
'          ',  4.50E-14,  0.00,   1260.00, 0.000, 0.000, 0.000, 0.000), &   
! B093 JPL2011
ratb_t('NO2       ','O3        ','NO3       ','          ','          ',& 
'          ',  1.20E-13,  0.00,   2450.00, 0.000, 0.000, 0.000, 0.000), &    
! B094 JPL2011
ratb_t('NO3       ','Br        ','BrO       ','NO2       ','          ',& 
'          ',  1.60E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B095 IUPAC2007
ratb_t('NO3       ','C5H8      ','ISON      ','          ','          ',& 
'          ',  3.15E-12,  0.00,    450.00, 0.000, 0.000, 0.000, 0.000), &   
! B096 IUPAC2007
ratb_t('NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',&  
'          ',  6.30E-15,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B097 IUPAC2007
ratb_t('NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',& 
'          ',  2.00E-12,  0.00,   2440.00, 0.000, 0.000, 0.000, 0.000), &    
! B098 MCMv3.2
ratb_t('NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',& 
'          ',  3.36E-12,  0.00,   1860.00, 0.000, 0.000, 0.000, 0.000), &    
! B099 IUPAC2007
ratb_t('NO3       ','Me2CO     ','HONO2     ','MeCOCH2OO ','          ',& 
'          ',  3.00E-17,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B100 IUPAC2007
ratb_t('NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ',& 
'          ',  1.40E-12,  0.00,   1860.00, 0.000, 0.000, 0.000, 0.000), &   
! B101 JPL2011
ratb_t('O(1D)     ','CH4       ','HCHO      ','H2        ','          ',& 
'          ',  9.00E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B102 JPL2011
ratb_t('O(1D)     ','CH4       ','HCHO      ','HO2       ','HO2       ',& 
'          ',  3.45E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B103 JPL2011
ratb_t('O(1D)     ','CH4       ','OH        ','MeOO      ','          ',& 
'          ',  1.31E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B104 JPL2011
ratb_t('O(1D)     ','CO2       ','O(3P)     ','CO2       ','          ',& 
'          ',  7.50E-11,  0.00,   -115.00, 0.000, 0.000, 0.000, 0.000), &    
! B105 IUPAC2008
ratb_t('O(1D)     ','H2        ','OH        ','H         ','          ',& 
'          ',  1.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B106 JPL2011
ratb_t('O(1D)     ','H2O       ','OH        ','OH        ','          ',& 
'          ',  1.63E-10,  0.00,    -60.00, 0.000, 0.000, 0.000, 0.000), &    
! B107 JPL2011
ratb_t('O(1D)     ','HBr       ','HBr       ','O(3P)     ','          ',& 
'          ',  3.00E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B108 JPL2011
ratb_t('O(1D)     ','HBr       ','OH        ','Br        ','          ',& 
'          ',  1.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B109 JPL2011
ratb_t('O(1D)     ','HCl       ','H         ','ClO       ','          ',& 
'          ',  3.60E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B110 JPL2011
ratb_t('O(1D)     ','HCl       ','O(3P)     ','HCl       ','          ',& 
'          ',  1.35E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B111 JPL2011
ratb_t('O(1D)     ','HCl       ','OH        ','Cl        ','          ',& 
'          ',  1.01E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B112 JPL2011
ratb_t('O(1D)     ','N2        ','O(3P)     ','N2        ','          ',& 
'          ',  2.15E-11,  0.00,   -110.00, 0.000, 0.000, 0.000, 0.000), &    
! B113 JPL2011
ratb_t('O(1D)     ','N2O       ','N2        ','O2        ','          ',& 
'          ',  4.60E-11,  0.00,    -20.00, 0.000, 0.000, 0.000, 0.000), &    
! B114 JPL2011
ratb_t('O(1D)     ','N2O       ','NO        ','NO        ','          ',& 
'          ',  7.30E-11,  0.00,    -20.00, 0.000, 0.000, 0.000, 0.000), &    
! B115 JPL2011
ratb_t('O(1D)     ','O2        ','O(3P)     ','O2        ','          ',& 
'          ',  3.30E-11,  0.00,    -55.00, 0.000, 0.000, 0.000, 0.000), &    
! B116 JPL2011
ratb_t('O(1D)     ','O3        ','O2        ','O(3P)     ','O(3P)     ',& 
'          ',  1.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B117 JPL2011
ratb_t('O(1D)     ','O3        ','O2        ','O2        ','          ',& 
'          ',  1.20E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B118 JPL2011
ratb_t('O(3P)     ','BrO       ','O2        ','Br        ','          ',& 
'          ',  1.90E-11,  0.00,   -230.00, 0.000, 0.000, 0.000, 0.000), &    
! B119 JPL2011
ratb_t('O(3P)     ','ClO       ','Cl        ','O2        ','          ',& 
'          ',  2.80E-11,  0.00,    -85.00, 0.000, 0.000, 0.000, 0.000), &    
! B120 JPL2011
ratb_t('O(3P)     ','ClONO2    ','ClO       ','NO3       ','          ',& 
'          ',  3.60E-12,  0.00,    840.00, 0.000, 0.000, 0.000, 0.000), &    
! B121 ????
ratb_t('O(3P)     ','H2        ','OH        ','H         ','          ',& 
'          ',  9.00E-18,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B122 JPL2011
ratb_t('O(3P)     ','H2O2      ','OH        ','HO2       ','          ',& 
'          ',  1.40E-12,  0.00,   2000.00, 0.000, 0.000, 0.000, 0.000), &   
! B123 JPL2011
ratb_t('O(3P)     ','HBr       ','OH        ','Br        ','          ',& 
'          ',  5.80E-12,  0.00,   1500.00, 0.000, 0.000, 0.000, 0.000), &    
! B124 JPL2011
ratb_t('O(3P)     ','HCHO      ','OH        ','CO        ','HO2       ',& 
'          ',  3.40E-11,  0.00,   1600.00, 0.000, 0.000, 0.000, 0.000), &    
! B125 JPL2011
ratb_t('O(3P)     ','HCl       ','OH        ','Cl        ','          ',& 
'          ',  1.00E-11,  0.00,   3300.00, 0.000, 0.000, 0.000, 0.000), &    
! B126 IUPAC2001
ratb_t('O(3P)     ','HO2       ','OH        ','O2        ','          ',& 
'          ',  2.70E-11,  0.00,   -224.00, 0.000, 0.000, 0.000, 0.000), &   
! B127 JPL2011
ratb_t('O(3P)     ','HOCl      ','OH        ','ClO       ','          ',& 
'          ',  1.70E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B128 JPL2011
ratb_t('O(3P)     ','NO2       ','NO        ','O2        ','          ',& 
'          ',  5.10E-12,  0.00,   -210.00, 0.000, 0.000, 0.000, 0.000), &    
! B129 IUPAC2009
ratb_t('O(3P)     ','NO3       ','O2        ','NO2       ','          ',& 
'          ',  1.70E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B130 JPL2011
ratb_t('O(3P)     ','O3        ','O2        ','O2        ','          ',& 
'          ',  8.00E-12,  0.00,   2060.00, 0.000, 0.000, 0.000, 0.000), &    
! B131 JPL2011
ratb_t('O(3P)     ','OClO      ','O2        ','ClO       ','          ',& 
'          ',  2.40E-12,  0.00,    960.00, 0.000, 0.000, 0.000, 0.000), &   
! B132 JPL2011
ratb_t('O(3P)     ','OH        ','O2        ','H         ','          ',& 
'          ',  1.80E-11,  0.00,   -180.00, 0.000, 0.000, 0.000, 0.000), &   
! B133 IUPAC2007* 
ratb_t('O3        ','C5H8      ','HO2       ','OH        ','          ',& 
'          ',  3.33E-15,  0.00,   1995.00, 0.750, 0.750, 0.000, 0.000), &   
! B134 IUPAC2007* 
ratb_t('O3        ','C5H8      ','MACR      ','HCHO      ','MACRO2    ',& 
'MeCO3     ',  3.33E-15,  0.00,   1995.00, 1.950, 1.740, 0.300, 0.300), &   
! B135 IUPAC2007* 
ratb_t('O3        ','C5H8      ','MeOO      ','HCOOH     ','CO        ',& 
'H2O2      ',  3.33E-15,  0.00,   1995.00, 0.240, 0.840, 0.420, 0.270)  &   
  /) 
  
ratb_defs_strattrop04 = (/   & 
! B136 IUPAC2007*
ratb_t('O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',& 
'CO        ',  2.13E-16,  0.00,   1520.00, 1.800, 0.900, 0.640, 0.440), &   
! B137 IUPAC2007*
ratb_t('O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',& 
'CO        ',  3.50E-16,  0.00,   2100.00, 1.800, 0.900, 0.640, 0.440), &   
! B138 IUPAC2007*
ratb_t('O3        ','MACR      ','OH        ','MeCO3     ','          ',& 
'          ',  2.13E-16,  0.00,   1520.00, 0.380, 0.200, 0.000, 0.000), &   
! B139 IUPAC2007*
ratb_t('O3        ','MACR      ','OH        ','MeCO3     ','          ',& 
'          ',  3.50E-16,  0.00,   2100.00, 0.380, 0.200, 0.000, 0.000), &   
! B140 JPL2011
ratb_t('OClO      ','NO        ','NO2       ','ClO       ','          ',& 
'          ',  2.50E-12,  0.00,    600.00, 0.000, 0.000, 0.000, 0.000), &    
! B141 IUPAC2007
ratb_t('OH        ','C2H6      ','H2O       ','EtOO      ','          ',& 
'          ',  6.90E-12,  0.00,   1000.00, 0.000, 0.000, 0.000, 0.000), &    
! B142 IUPAC2007 see also asad_bimol
ratb_t('OH        ','C3H8      ','i-PrOO    ','H2O       ','          ',& 
'          ',  7.60E-12,  0.00,    585.00, 0.000, 0.000, 0.000, 0.000), &   
! B143 IUPAC2007 see also asad_bimol
ratb_t('OH        ','C3H8      ','n-PrOO    ','H2O       ','          ',& 
'          ',  7.60E-12,  0.00,    585.00, 0.000, 0.000, 0.000, 0.000), & 
! B144 IUPAC2009 
ratb_t('OH        ','C5H8      ','ISO2      ','          ','          ',&  
'          ',  2.70E-11,  0.00,   -390.00, 0.000, 0.000, 0.000, 0.000), &   
! B145 JPL2011
ratb_t('OH        ','CH4       ','H2O       ','MeOO      ','          ',& 
'          ',  2.45E-12,  0.00,   1775.00, 0.000, 0.000, 0.000, 0.000), &    
! B146 IUPAC2005 see also asad_bimol
ratb_t('OH        ','CO        ','HO2       ','          ','          ',& 
'          ',  1.44E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), & 
! B147 JPL2011
ratb_t('OH        ','ClO       ','HCl       ','O2        ','          ',& 
'          ',  6.00E-13,  0.00,   -230.00, 0.000, 0.000, 0.000, 0.000), &   
! B148 JPL2011 
ratb_t('OH        ','ClO       ','HO2       ','Cl        ','          ',& 
'          ',  7.40E-12,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000), &    
! B149 JPL2011 
ratb_t('OH        ','ClONO2    ','HOCl      ','NO3       ','          ',& 
'          ',  1.20E-12,  0.00,    330.00, 0.000, 0.000, 0.000, 0.000), &   
! B150 IUPAC2007
ratb_t('OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',& 
'          ',  4.90E-12,  0.00,   -405.00, 0.000, 0.000, 0.000, 0.000), &   
! B151 MCMv3.2
ratb_t('OH        ','EtOOH     ','H2O       ','EtOO      ','          ',& 
'          ',  1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &   
! B152 MCMv3.2
ratb_t('OH        ','EtOOH     ','H2O       ','MeCHO     ','OH        ',& 
'          ',  8.01E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B153 JPL2011
ratb_t('OH        ','H2        ','H2O       ','HO2       ','          ',& 
'          ',  2.80E-12,  0.00,   1800.00, 0.000, 0.000, 0.000, 0.000), &   
! B154 IUPAC2001
ratb_t('OH        ','H2O2      ','H2O       ','HO2       ','          ',& 
'          ',  2.90E-12,  0.00,    160.00, 0.000, 0.000, 0.000, 0.000), &  
! B155 IUPAC2007
ratb_t('OH        ','HACET     ','MGLY      ','HO2       ','          ',& 
'          ',  1.60E-12,  0.00,   -305.00, 0.000, 0.000, 0.000, 0.000), &    
! B156 JPL2011 
ratb_t('OH        ','HBr       ','H2O       ','Br        ','          ',& 
'          ',  5.50E-12,  0.00,   -200.00, 0.000, 0.000, 0.000, 0.000), &    
! B157 IUPAC2007
ratb_t('OH        ','HCHO      ','H2O       ','HO2       ','CO        ',& 
'          ',  5.40E-12,  0.00,   -135.00, 0.000, 0.000, 0.000, 0.000), &   
! B158 IUPAC2007
ratb_t('OH        ','HCOOH     ','HO2       ','          ','          ',& 
'          ',  4.50E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B159 JPL2011 
ratb_t('OH        ','HCl       ','H2O       ','Cl        ','          ',& 
'          ',  1.80E-12,  0.00,    250.00, 0.000, 0.000, 0.000, 0.000), &   
! B160 JPL2011 
ratb_t('OH        ','HO2       ','H2O       ','          ','          ',& 
'          ',  4.80E-11,  0.00,   -250.00, 0.000, 0.000, 0.000, 0.000), &   
! B161 IUPAC2007
ratb_t('OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',& 
'          ',  3.20E-13,  0.00,   -690.00, 0.000, 0.000, 0.000, 0.000), &    
! B162 JPL2011
ratb_t('OH        ','HOCl      ','ClO       ','H2O       ','          ',& 
'          ',  3.00E-12,  0.00,    500.00, 0.000, 0.000, 0.000, 0.000), &   
! B163 IUPAC2004
ratb_t('OH        ','HONO      ','H2O       ','NO2       ','          ',& 
'          ',  2.50E-12,  0.00,   -260.00, 0.000, 0.000, 0.000, 0.000), &   
! B164 IUPAC2004 see also asad_bimol 
ratb_t('OH        ','HONO2     ','H2O       ','NO3       ','          ',&  
'          ',  2.40E-14,  0.00,   -460.00, 0.000, 0.000, 0.000, 0.000), & 
! B165 Poschl00 
ratb_t('OH        ','ISON      ','HACET     ','NALD      ','          ',& 
'          ',  1.30E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B166 Poschl00
ratb_t('OH        ','ISOOH     ','MACR      ','OH        ','          ',&  
'          ',  1.00E-10,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B167 IUPAC2007
ratb_t('OH        ','MACR      ','MACRO2    ','          ','          ',& 
'          ',  1.30E-12,  0.00,   -610.00, 0.000, 0.000, 0.000, 0.000), &   
! B168 IUPAC2007
ratb_t('OH        ','MACR      ','MACRO2    ','          ','          ',& 
'          ',  4.00E-12,  0.00,   -380.00, 0.000, 0.000, 0.000, 0.000), &    
! B169 MCMv3.2
ratb_t('OH        ','MACROOH   ','MACRO2    ','          ','          ',& 
'          ',  3.77E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B170 IUPAC2008
ratb_t('OH        ','MGLY      ','MeCO3     ','CO        ','          ',& 
'          ',  1.90E-11,  0.00,   -575.00, 0.000, 0.000, 0.000, 0.000), &    
! B171 IUPAC2006
ratb_t('OH        ','MPAN      ','HACET     ','NO2       ','          ',& 
'          ',  2.90E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B172 IUPAC2007
ratb_t('OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',& 
'          ',  1.70E-14,  0.00,   -423.00, 0.000, 0.000, 0.000, 0.000), &   
! B173 IUPAC2007
ratb_t('OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',& 
'          ',  8.80E-12,  0.00,   1320.00, 0.000, 0.000, 0.000, 0.000), &   
! B174 IUPAC2009
ratb_t('OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',&  
'          ',  4.70E-12,  0.00,   -345.00, 0.000, 0.000, 0.000, 0.000), &   
! B175 MCMv3.2 
ratb_t('OH        ','MeCO2H    ','MeOO      ','          ','          ',&  
'          ',  8.00E-13,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B176 MCMv3.2
ratb_t('OH        ','MeCO3H    ','MeCO3     ','          ','          ',&  
'          ',  3.70E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B177 MCMv3.2
ratb_t('OH        ','MeCOCH2OOH','H2O       ','MeCOCH2OO ','          ',& 
'          ',  1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &   
! B178 MCMv3.2
ratb_t('OH        ','MeCOCH2OOH','OH        ','MGLY      ','          ',&  
'          ',  8.39E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B179 IUPAC2006
ratb_t('OH        ','MeOH      ','HO2       ','HCHO      ','          ',& 
'          ',  2.85E-12,  0.00,    345.00, 0.000, 0.000, 0.000, 0.000), &   
! B180 IUPAC2006
ratb_t('OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',& 
'          ',  4.00E-13,  0.00,    845.00, 0.000, 0.000, 0.000, 0.000)  &   
  /) 
  
ratb_defs_strattrop05 = (/   & 
! B060a added Alex 
ratb_t('HO2       ','NO        ','HONO2     ','          ','          ',& 
'          ',  3.60E-12,  0.00,   -270.00, 0.000, 0.000, 0.000, 0.000) ,&   
! B181 IUPAC2007 
ratb_t('OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',& 
'          ',  2.12E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &   
! B182 IUPAC2007
ratb_t('OH        ','MeOOH     ','H2O       ','MeOO      ','          ',& 
'          ',  1.89E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &   
! B183 IUPAC2009
ratb_t('OH        ','NALD      ','HCHO      ','CO        ','NO2       ',& 
'          ',  4.70E-12,  0.00,   -345.00, 0.000, 0.000, 0.000, 0.000), &   
! B184 JPL2011
ratb_t('OH        ','NO3       ','HO2       ','NO2       ','          ',& 
'          ',  2.20E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B185 JPL2011
ratb_t('OH        ','O3        ','HO2       ','O2        ','          ',& 
'          ',  1.70E-12,  0.00,    940.00, 0.000, 0.000, 0.000, 0.000), &    
! B186 JPL2011
ratb_t('OH        ','OClO      ','HOCl      ','O2        ','          ',& 
'          ',  1.40E-12,  0.00,   -600.00, 0.000, 0.000, 0.000, 0.000), &   
! B187 IUPAC2001
ratb_t('OH        ','OH        ','H2O       ','O(3P)     ','          ',& 
'          ',  6.31E-14,  2.60,   -945.00, 0.000, 0.000, 0.000, 0.000), &   
! B188 MCMv3.2
ratb_t('OH        ','PAN       ','HCHO      ','NO2       ','H2O       ',& 
'          ',  3.00E-14,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B189 MCMv3.2
ratb_t('OH        ','PPAN      ','MeCHO     ','NO2       ','H2O       ',& 
'          ',  1.27E-12,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B190 MCMv3.2 
ratb_t('OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ',& 
'          ',  1.66E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &    
! B191 MCMv3.2 
ratb_t('OH        ','i-PrOOH   ','i-PrOO    ','H2O       ','          ',& 
'          ',  1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &   
! B192 MCMv3.2
ratb_t('OH        ','n-PrOOH   ','EtCHO     ','H2O       ','OH        ',& 
'          ',  1.10E-11,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000), &   
! B193 MCMv3.2
ratb_t('OH        ','n-PrOOH   ','n-PrOO    ','H2O       ','          ',& 
'          ',  1.90E-12,  0.00,   -190.00, 0.000, 0.000, 0.000, 0.000), &    
! B194 IUPAC2005
ratb_t('i-PrOO    ','NO        ','Me2CO     ','HO2       ','NO2       ',& 
'          ',  2.70E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &   
! B195 MCMv3.2
ratb_t('i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',& 
'          ',  2.70E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000), &    
! B196 IUPAC2005
ratb_t('n-PrOO    ','NO        ','EtCHO     ','HO2       ','NO2       ',& 
'          ',  2.90E-12,  0.00,   -350.00, 0.000, 0.000, 0.000, 0.000), &    
! B197 MCM3.2
ratb_t('n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',& 
'          ',  2.70E-12,  0.00,   -360.00, 0.000, 0.000, 0.000, 0.000)  &  
  /) 
  
!---------------------------------------------------------------------- 
! NOTES: CheST Bimolecular Reactions 
!---------------------------------------------------------------------- 
! B001 Br+Cl2O2 -> BrCl Cl O2 JPL2011   
! B001 k(298)=3.3E-12. Both JPL2011 and IUPAC2006 agree. Reaction forms 
! B001 ClOO rather than Cl and O2. 
!---------------------------------------------------------------------- 
! B002 Br+HCHO -> HBr CO HO2 JPL2011   
! B002 k(298)=1.1E-12. Unchanged since JPL2009. No reference by IUPAC. 
! B002 JPL cite two studies of which this value is derived from by a 
! B002 weighted mean. 
!---------------------------------------------------------------------- 
! B003 Br+HO2 -> HBr O2  JPL2011   
! B003 k(298)=1.7E-12. Unchanged since JPL2006. IUPAC2003 recommend 
! B003 k=7.7E-12*exp(-450/temp). 
!---------------------------------------------------------------------- 
! B004 Br+O3 -> BrO O2  JPL2011   
! B004 k(298)=1.2E-12. Unchanged since JPL2009. IUPAC2003 recommend k=1 
! B004 .7E-11exp(-800/temp) 
!---------------------------------------------------------------------- 
! B005 Br+OClO -> BrO ClO  JPL2011   
! B005 k(298)=3.4E-13. Originally from JPL1990 
!---------------------------------------------------------------------- 
! B006 BrO+BrO -> Br Br O2 JPL2011   
! B006 Originally from JPL1997 
!---------------------------------------------------------------------- 
! B007 BrO+ClO -> Br Cl O2 JPL2011   
! B007 k(298)=5.5E-12. Unchanged since JPL2009. IUPAC2004 recommend 
! B007 k=2.9E-12*exp(220/temp) with k(298)=6.1E-12. Both JPL and IUPAC 
! B007 recommended products are Br + ClOO. 
!---------------------------------------------------------------------- 
! B008 BrO+ClO -> Br OClO  JPL2011   
! B008 k(298)=6.0E-12. Unchanged since JPL2009. IUPAC2004 recommend 
! B008 k=1.60E-13*exp(430/temp) with k(298)=6.8E-12 (at 220 K k = 
! B008 1.13E-12 c.f 1.16E-12). This channel accounts for around 60% of 
! B008 the total for this reaction. 
!---------------------------------------------------------------------- 
! B009 BrO+ClO -> BrCl O2  JPL2011   
! B009 k(298)=1.1E-12. Unchanged since JPL2009. IUPAC2004 recommend 
! B009 k=5.8E-13*exp(170/temp) with k(298)=1.0E-12. Branching ratios put 
! B009 this channel at ~8%. 
!---------------------------------------------------------------------- 
! B010 BrO+HO2 -> HBr O3  JPL2011   
! B010 Remove this reaction since both JPL2006 and IUPAC2003 suggest 
! B010 that the branching ratio for this channel is 0.0 
!---------------------------------------------------------------------- 
! B011 BrO+HO2 -> HOBr O2  JPL2011   
! B011 k(298)=2.1E-11. Unchanged since JPL2006. IUPAC2003 recommend k=4 
! B011 .5E-12exp(500/temp) with k(298)=2.4E-11. Both JPL and IUPAC 
! B011 suggest that essentially this is the only channel for this 
! B011 reaction. 
!---------------------------------------------------------------------- 
! B012 BrO+NO -> Br NO2  JPL2011   
! B012 k(298)=2.1E-11. Originally from JPL1982 
!---------------------------------------------------------------------- 
! B013 BrO+OH -> Br HO2  JPL2011   
! B013 k(298)=3.9E-11. Unchanged since JPL2006. IUPAC2003 recommend 
! B013 k=1.8E-11*exp(250/temp) with k(298)=4.1E-11. The temperature 
! B013 dependence comes from the study by Bedjanian et al. 
!---------------------------------------------------------------------- 
! B014 CF2Cl2+O(1D) -> Cl ClO  JPL2011   
! B014 k(298)=1.4E-10. Unchanged since JPL1992. Both JPL1992 and 
! B014 IUPAC2011 agree. The ClO channel is ~87%. 
!---------------------------------------------------------------------- 
! B015 CFCl3+O(1D) -> Cl Cl ClO JPL2011   
! B015 k(298)=2.3E-10. Unchanged since JPL1992. Both JPL1992 and 
! B015 IUPAC2011 agree. The ClO channel is ~80%. 
!---------------------------------------------------------------------- 
! B016 Cl+CH4 -> HCl MeOO  JPL2011   
! B016 k(298)=1.0E-13. IUPAC recommend k=6.6E-12*exp(-1240/temp). 
!---------------------------------------------------------------------- 
! B017 Cl+Cl2O2 -> Cl Cl Cl IUPAC2006   
! B017 k(298)=1.0E-10. Probable error in JPL2011 gives beta=9.4E-11! 
! B017 IUPAC2006 recommend k = 7.6E-11*exp(65/temp) based on study by 
! B017 Ingham et al. 
!---------------------------------------------------------------------- 
! B018 Cl+ClONO2 -> Cl Cl NO3 JPL2011   
! B018 IUPAC similar (6.2 and -145). Actual reaction forms Cl2 rather 
! B018 than 2Cl. 
!---------------------------------------------------------------------- 
! B019 Cl+H2 -> HCl H  JPL2011   
! B019 IUPAC similar (3.9 and 2310) 
!---------------------------------------------------------------------- 
! B020 Cl+H2O2 -> HCl HO2  JPL2011   
! B020 IUPAC identical. 
!---------------------------------------------------------------------- 
! B021 Cl+HCHO -> HCl CO HO2 JPL2011   
! B021 IUPAC similar (8.1 and 34). 
!---------------------------------------------------------------------- 
! B022 Cl+HO2 -> ClO OH  JPL2011   
! B022 IUPAC gives total rate + k2/k. Hard to compare 
!---------------------------------------------------------------------- 
! B023 Cl+HO2 -> HCl O2  JPL2011   
! B023 IUPAC gives total rate + k2/k. Hard to compare 
!---------------------------------------------------------------------- 
! B024 Cl+HOCl -> Cl Cl OH JPL2011   
! B024 No IUPAC rate 
!---------------------------------------------------------------------- 
! B025 Cl+MeOOH -> HCl MeOO  JPL2011   
! B025 IUPAC rate similar (5.9E-11) 
!---------------------------------------------------------------------- 
! B026 Cl+NO3 -> ClO NO2  JPL2011   
! B026 IUPAC rates identical 
!---------------------------------------------------------------------- 
! B027 Cl+O3 -> ClO O2  JPL2011   
! B027 IUPAC rates similar (2.8 and 250) 
!---------------------------------------------------------------------- 
! B028 Cl+OClO -> ClO ClO  JPL2011   
! B028 IUPAC rates similar (3.2 and -170) 
!---------------------------------------------------------------------- 
! B029 ClO+ClO -> Cl Cl O2 JPL2011   
! B029 IUPAC rates identical 
!---------------------------------------------------------------------- 
! B030 ClO+ClO -> Cl Cl O2 JPL2011   
! B030 IUPAC rates identical 
!---------------------------------------------------------------------- 
! B031 ClO+ClO -> Cl OClO  JPL2011   
! B031 IUPAC rates identical 
!---------------------------------------------------------------------- 
! B032 ClO+HO2 -> HCl O3  JPL2011   
! B032 No rate given in IUPAC or JPL. Reaction should probably be 
! B032 removed. 
!---------------------------------------------------------------------- 
! B033 ClO+HO2 -> HOCl O2  JPL2011   
! B033 IUPAC rates similar (2.2 and -340) 
!---------------------------------------------------------------------- 
! B034 ClO+MeOO -> Cl HCHO HO2 JPL2011   
! B034 IUPAC provides a more detailed breakdown of products. Would 
! B034 require a modification of asad_bimol to treat properly. Plus 
! B034 retain policy of using JPL for ClOx reactions. 
!---------------------------------------------------------------------- 
! B035 ClO+NO -> Cl NO2  JPL2011   
! B035 IUPAC rates very similar (6.2 and -295) 
!---------------------------------------------------------------------- 
! B036 ClO+NO3 -> Cl O2 NO2 JPL2011   
! B036 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6 
! B036 going to other channel 
!---------------------------------------------------------------------- 
! B037 ClO+NO3 -> OClO NO2  JPL2011   
! B037 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6 
! B037 going to other channel. If use JPL approach reaction could be 
! B037 removed. 
!---------------------------------------------------------------------- 
! B038 EtCO3+NO -> EtOO CO2 NO2 IUPAC2002   
! B038 No JPL rate 
!---------------------------------------------------------------------- 
! B039 EtCO3+NO3 -> EtOO CO2 NO2 MCMv3.2   
! B039 No IUPAC/JPL rate 
!---------------------------------------------------------------------- 
! B040 EtOO+MeCO3 -> MeCHO HO2 MeOO IUPAC2002   
! B040 No JPL rate 
!---------------------------------------------------------------------- 
! B041 EtOO+NO -> MeCHO HO2 NO2 IUPAC2005   
! B041 No JPL rate 
!---------------------------------------------------------------------- 
! B042 EtOO+NO3 -> MeCHO HO2 NO2 IUPAC2008   
! B042 No JPL rate 
!---------------------------------------------------------------------- 
! B043 H+HO2 -> H2 O2  JPL2011   
! B043 No IUPAC rate 
!---------------------------------------------------------------------- 
! B044 H+HO2 -> O(3P) H2O  JPL2011   
! B044 No IUPAC rate 
!---------------------------------------------------------------------- 
! B045 H+HO2 -> OH OH  JPL2011   
! B045 No IUPAC rate 
!---------------------------------------------------------------------- 
! B046 H+NO2 -> OH NO  JPL2011   
! B046 No IUPAC rate 
!---------------------------------------------------------------------- 
! B047 H+O3 -> OH O2  JPL2011   
! B047 No IUPAC rate 
!---------------------------------------------------------------------- 
! B048 HO2+EtCO3 -> O2 EtCO3H  MCMv3.2*   
! B048 No IUPAC/JPL rate. No exact correspondence in MCM. Take total 
! B048 rate and subtract O3 channel. 
!---------------------------------------------------------------------- 
! B049 HO2+EtCO3 -> O3 EtCO2H  MVMv3.2   
! B049 No IUPAC/JPL rate 
!---------------------------------------------------------------------- 
! B050 HO2+EtOO -> EtOOH   IUPAC2011   
! B050 - 
!---------------------------------------------------------------------- 
! B051 HO2+HO2 -> H2O2   JPL2011 see also asad_bimol  
! B051 IUPAC slightly different (2.2 & -600) . Rate modified in the 
! B051 presence of H2O in asad bimol. 
!---------------------------------------------------------------------- 
! B052 HO2+ISO2 -> ISOOH   Poschl00   
! B052 - 
!---------------------------------------------------------------------- 
! B053 HO2+MACRO2 -> MACROOH   Poschl00   
! B053 - 
!---------------------------------------------------------------------- 
! B054 HO2+MeCO3 -> MeCO2H O3  IUPAC2009   
! B054 - 
!---------------------------------------------------------------------- 
! B055 HO2+MeCO3 -> MeCO3H   IUPAC2009   
! B055 - 
!---------------------------------------------------------------------- 
! B056 HO2+MeCO3 -> OH MeOO  IUPAC2009   
! B056 - 
!---------------------------------------------------------------------- 
! B057 HO2+MeCOCH2OO -> MeCOCH2OOH   IUPAC2009   
! B057 Use measured value rather than the expression from the MCM (which 
! B057 does include T dependence) 
!---------------------------------------------------------------------- 
! B058 HO2+MeOO -> HCHO H2O  IUPAC2009 see also asad_bimol  
! B058 Total rate given. Split between channels controlled by asad_bimol 
!---------------------------------------------------------------------- 
! B059 HO2+MeOO -> MeOOH O2  IUPAC2009 see also asad_bimol  
! B059 Total rate given. Split between channels controlled by asad_bimol 
!---------------------------------------------------------------------- 
! B060 HO2+NO -> OH NO2  JPL2011   
! B060 IUPAC values similar (3.45 and 270) 
!---------------------------------------------------------------------- 
! B061 HO2+NO3 -> OH NO2 O2 JPL2011   
! B061 IUPAC value slightly larger (4 cf 3.5) 
!---------------------------------------------------------------------- 
! B062 HO2+O3 -> OH O2  IUPAC2001   
! B062 Use more sophisticated expression of IUPAC. JPL have simple term 
! B062 1*10^-14 and 490. 
!---------------------------------------------------------------------- 
! B063 HO2+i-PrOO -> i-PrOOH   MCMv3.2   
! B063 No rate in IUPAC or JPL 
!---------------------------------------------------------------------- 
! B064 HO2+n-PrOO -> n-PrOOH   MCMv3.2   
! B064 No rate in IUPAC or JPL 
!---------------------------------------------------------------------- 
! B065 ISO2+ISO2 -> MACR MACR HCHOHO2 Poschl00   
! B065 - 
!---------------------------------------------------------------------- 
! B066 MACRO2+MACRO2 -> HACET MGLY HCHOCO Poschl00   
! B066 - 
!---------------------------------------------------------------------- 
! B067 MACRO2+MACRO2 -> HO2   Poschl00   
! B067 - 
!---------------------------------------------------------------------- 
! B068 MeBr+Cl -> Br HCl  JPL2011   
! B068 - 
!---------------------------------------------------------------------- 
! B069 MeBr+O(1D) -> Br OH  JPL2011   
! B069 - 
!---------------------------------------------------------------------- 
! B070 MeBr+OH -> Br H2O  JPL2011   
! B070 - 
!---------------------------------------------------------------------- 
! B071 MeCO3+NO -> MeOO CO2 NO2 IUPAC2002   
! B071 - 
!---------------------------------------------------------------------- 
! B072 MeCO3+NO3 -> MeOO CO2 NO2 IUPAC2008   
! B072 - 
!---------------------------------------------------------------------- 
! B073 MeCOCH2OO+NO -> MeCO3 HCHO NO2 MCMv3.2   
! B073 No rate in IUPAC or JPL 
!---------------------------------------------------------------------- 
! B074 MeCOCH2OO+NO3 -> MeCO3 HCHO NO2 MCMv3.2   
! B074 No rate in IUPAC or JPL 
!---------------------------------------------------------------------- 
! B075 MeOO+MeCO3 -> HO2 HCHO MeOO IUPAC2002   
! B075 JPL only gives total rates 
!---------------------------------------------------------------------- 
! B076 MeOO+MeCO3 -> MeCO2H HCHO  IUPAC2002   
! B076 JPL only gives total rates 
!---------------------------------------------------------------------- 
! B077 MeOO+MeOO -> HO2 HO2 HCHOHCHO IUPAC2002 see also asad_bimol  
! B077 IUPAC gives total k and k2 as f(T) resulting in complicated 
! B077 expression evaluated in asad_bimol. JPL only gives a total rate 
! B077 (9.5 and -390). 
!---------------------------------------------------------------------- 
! B078 MeOO+MeOO -> MeOH HCHO  IUPAC2002 see also asad_bimol  
! B078 IUPAC gives total k and k2 as f(T) resulting in complicated 
! B078 expression evaluated in asad_bimol. JPL only gives a total rate 
! B078 (9.5 and -390). 
!---------------------------------------------------------------------- 
! B079 MeOO+NO -> HO2 HCHO NO2 IUPAC2005   
! B079 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO 
! B079 and HO2 and 0.1% goes to MeONO2. Where does this partition come 
! B079 from? 
!---------------------------------------------------------------------- 
! B080 MeOO+NO -> MeONO2   IUPAC2005   
! B080 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO 
! B080 and HO2 and 0.1% goes to MeONO2. Where does this partition come 
! B080 from? 
!---------------------------------------------------------------------- 
! B081 MeOO+NO3 -> HO2 HCHO NO2 MCMv3.2   
! B081 No JPL value. 
!---------------------------------------------------------------------- 
! B082 N+NO -> N2 O(3P)  JPL2011   
! B082 No IUPAC value. 
!---------------------------------------------------------------------- 
! B083 N+NO2 -> N2O O(3P)  JPL2011   
! B083 No IUPAC value. 
!---------------------------------------------------------------------- 
! B084 N+O2 -> NO O(3P)  JPL2011   
! B084 IUPAC  values slightly different (5.1 and -198). Different 
! B084 choices of results and fitting. 
!---------------------------------------------------------------------- 
! B085 N2O5+H2O -> HONO2 HONO2  ????   
! B085 Where is this from? Both JPL and IUPAC only give upper limits 
! B085 (which are lower than this value) 
!---------------------------------------------------------------------- 
! B086 NO+ISO2 -> ISON   Poschl00   
! B086 - 
!---------------------------------------------------------------------- 
! B087 NO+ISO2 -> NO2 MACR HCHOHO2 Poschl00   
! B087 - 
!---------------------------------------------------------------------- 
! B088 NO+MACRO2 -> MGLY HCHO HO2 Poschl00   
! B088 Note to accomodate all the products split into two reactions and 
! B088 half the rates. 
!---------------------------------------------------------------------- 
! B089 NO+MACRO2 -> NO2 MeCO3 HACETCO Poschl00   
! B089 Note to accomodate all the products split into two reactions and 
! B089 half the rates. 
!---------------------------------------------------------------------- 
! B090 NO+NO3 -> NO2 NO2  JPL2011   
! B090 Discrepancy with IUPAC (1.8E-11 & beta =170) related to 
! B090 differences in averaging. 
!---------------------------------------------------------------------- 
! B091 NO+O3 -> NO2   JPL2011   
! B091 Considerable discrepancy with IUPAC (1.4E-12 & beta=1310) 
! B091 resulting in differences ~10% between the two. Use JPL as include 
! B091 the extra (room temperature only though) of Stedman & Niki and 
! B091 Bemand et al. Believe difference arises from fitting & averaging 
! B091 vs averaging & fitting. 
!---------------------------------------------------------------------- 
! B092 NO2+NO3 -> NO NO2 O2 JPL2011   
! B092 This is a termolecular reaction. Why is it here? IUPAC don't 
! B092 recognise the reaction and JPL note its existence is not 
! B092 definitively proven. However they do use studies that infer the 
! B092 rates from thermal decomposition of N2O5. 
!---------------------------------------------------------------------- 
! B093 NO2+O3 -> NO3   JPL2011   
! B093 Slight disagreement with IUPAC rates (1.4E-11 & -2470). Use JPL 
! B093 by default. 
!---------------------------------------------------------------------- 
! B094 NO3+Br -> BrO NO2  JPL2011   
! B094 Both IUPAC & JPL recommend approach of Mellouki et al (1989) 
!---------------------------------------------------------------------- 
! B095 NO3+C5H8 -> ISON   IUPAC2007   
! B095 - 
!---------------------------------------------------------------------- 
! B096 NO3+EtCHO -> HONO2 EtCO3  IUPAC2007   
! B096 No JPL rates. No T dependence given buy IUPAC. 
!---------------------------------------------------------------------- 
! B097 NO3+HCHO -> HONO2 HO2 CO IUPAC2007   
! B097 No direct measurements of T dependence. Infer from T dependence 
! B097 of MeCHO + NO3.  JPL values at 298K same as IUPAC. 
!---------------------------------------------------------------------- 
! B098 NO3+MGLY -> MeCO3 CO HONO2 MCMv3.2   
! B098 No rate in JPL/IUPAC. 2.4*IUPAC MeCHO + NO3. 
!---------------------------------------------------------------------- 
! B099 NO3+Me2CO -> HONO2 MeCOCH2OO  IUPAC2007   
! B099 Number is actually an upper limit. Use it as no other number 
! B099 available (JPL don't provide one either) 
!---------------------------------------------------------------------- 
! B100 NO3+MeCHO -> HONO2 MeCO3  IUPAC2007   
! B100 - 
!---------------------------------------------------------------------- 
! B101 O(1D)+CH4 -> HCHO H2  JPL2011   
! B101 IUPAC values slightly lower (7.50E-12) 
!---------------------------------------------------------------------- 
! B102 O(1D)+CH4 -> HCHO HO2 HO2 JPL2011   
! B102 IUPAC values same 
!---------------------------------------------------------------------- 
! B103 O(1D)+CH4 -> OH MeOO  JPL2011   
! B103 IUPAC values slightly lower (1.05E-10) 
!---------------------------------------------------------------------- 
! B104 O(1D)+CO2 -> O(3P) CO2  JPL2011   
! B104 - 
!---------------------------------------------------------------------- 
! B105 O(1D)+H2 -> OH H  IUPAC2008   
! B105 JPL values identical 
!---------------------------------------------------------------------- 
! B106 O(1D)+H2O -> OH OH  JPL2011   
! B106 JPL value dominated by the work of Dunlea & Ravishankra. Slightly 
! B106 lower than the IUPAC values (5% at room temperature) 
!---------------------------------------------------------------------- 
! B107 O(1D)+HBr -> HBr O(3P)  JPL2011   
! B107 Use Wine al for BR to O3P 
!---------------------------------------------------------------------- 
! B108 O(1D)+HBr -> OH Br  JPL2011   
! B108 Use Wine et al O3P BR and assume rest of HBR gos to OH + Br 
!---------------------------------------------------------------------- 
! B109 O(1D)+HCl -> H ClO  JPL2011   
! B109 Use Wine et al for branching ratios. 
!---------------------------------------------------------------------- 
! B110 O(1D)+HCl -> O(3P) HCl  JPL2011   
! B110 Use Wine et al for branching ratios. 
!---------------------------------------------------------------------- 
! B111 O(1D)+HCl -> OH Cl  JPL2011   
! B111 Use Wine et al for branching ratios. 
!---------------------------------------------------------------------- 
! B112 O(1D)+N2 -> O(3P) N2  JPL2011   
! B112 IUPAC identical 
!---------------------------------------------------------------------- 
! B113 O(1D)+N2O -> N2 O2  JPL2011   
! B113 Use JPL values (slightly higher than IUPAC 4.3e-11) as slightly 
! B113 newer and include T dependence 
!---------------------------------------------------------------------- 
! B114 O(1D)+N2O -> NO NO  JPL2011   
! B114 Use JPL values (slightly lower than IUPAC 7.6e-11) as slightly 
! B114 newer and include T dependence 
!---------------------------------------------------------------------- 
! B115 O(1D)+O2 -> O(3P) O2  JPL2011   
! B115 IUPAC values slightly different (3.2 and -67) 
!---------------------------------------------------------------------- 
! B116 O(1D)+O3 -> O2 O(3P) O(3P) JPL2011   
! B116 IUPAC values identical 
!---------------------------------------------------------------------- 
! B117 O(1D)+O3 -> O2 O2  JPL2011   
! B117 IUPAC values identical 
!---------------------------------------------------------------------- 
! B118 O(3P)+BrO -> O2 Br  JPL2011   
! B118 - 
!---------------------------------------------------------------------- 
! B119 O(3P)+ClO -> Cl O2  JPL2011   
! B119 Updated JPL 10-6 
!---------------------------------------------------------------------- 
! B120 O(3P)+ClONO2 -> ClO NO3  JPL2011   
! B120 Updated JPL 10-6 
!---------------------------------------------------------------------- 
! B121 O(3P)+H2 -> OH H  ????   
! B121 Not defined in IUPAC or JPL. Where does this come from? 
!---------------------------------------------------------------------- 
! B122 O(3P)+H2O2 -> OH HO2  JPL2011   
! B122 IUPAC values identical 
!---------------------------------------------------------------------- 
! B123 O(3P)+HBr -> OH Br  JPL2011   
! B123 - 
!---------------------------------------------------------------------- 
! B124 O(3P)+HCHO -> OH CO HO2 JPL2011   
! B124 Value not given in IUPAC (are we sure?) 
!---------------------------------------------------------------------- 
! B125 O(3P)+HCl -> OH Cl  JPL2011   
! B125 - 
!---------------------------------------------------------------------- 
! B126 O(3P)+HO2 -> OH O2  IUPAC2001   
! B126 JPL are the same to one significant figure. 
!---------------------------------------------------------------------- 
! B127 O(3P)+HOCl -> OH ClO  JPL2011   
! B127 - 
!---------------------------------------------------------------------- 
! B128 O(3P)+NO2 -> NO O2  JPL2011   
! B128 Note T dependence in IUPAC is slightly different (198 cf 210) 
!---------------------------------------------------------------------- 
! B129 O(3P)+NO3 -> O2 NO2  IUPAC2009   
! B129 Rate smaller in JPL2011 (1.1e-11) who prefer the results of 
! B129 Graham & Johnston (1978) to Canosa-Mas et al. (1989). Though 
! B129 IUPAC note that the results agree within uncertainties. 
!---------------------------------------------------------------------- 
! B130 O(3P)+O3 -> O2 O2  JPL2011   
! B130 IUPAC identical 
!---------------------------------------------------------------------- 
! B131 O(3P)+OClO -> O2 ClO  JPL2011   
! B131 - 
!---------------------------------------------------------------------- 
! B132 O(3P)+OH -> O2 H  JPL2011   
! B132 IUPAC slightly different (2.4 & -110). JPL include slightly more 
! B132 results and newer results. However main difference probably 
! B132 derives from fitting. 
!---------------------------------------------------------------------- 
! B133 O3+C5H8 -> HO2 OH  IUPAC2007*   
! B133 Total rates from IUPAC. Arbitrarily split into three reactions. 
! B133 Fractional products given by Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B134 O3+C5H8 -> MACR HCHO MACRO2MeCO3 IUPAC2007*   
! B134 Total rates from IUPAC. Arbitrarily split into three reactions. 
! B134 Fractional products given by Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B135 O3+C5H8 -> MeOO HCOOH COH2O2 IUPAC2007*   
! B135 Total rates from IUPAC. Arbitrarily split into three reactions. 
! B135 Fractional products given by Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B136 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*   
! B136 Complicated expression. Rate is average of IUPAC MACR + O3 and 
! B136 MVK + O3. Split into further two reactions to allow all products 
! B136 to be included. Hence rates are IUPAC/4. Fractional Products from 
! B136 Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B137 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*   
! B137 Complicated expression. Rate is average of IUPAC MACR + O3 and 
! B137 MVK + O3. Split into further two reactions to allow all products 
! B137 to be included. Hence rates are IUPAC/4. Fractional Products from 
! B137 Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B138 O3+MACR -> OH MeCO3  IUPAC2007*   
! B138 Complicated expression. Rate is average of IUPAC MACR + O3 and 
! B138 MVK + O3. Split into further two reactions to allow all products 
! B138 to be included. Hence rates are IUPAC/4. Fractional Products from 
! B138 Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B139 O3+MACR -> OH MeCO3  IUPAC2007*   
! B139 Complicated expression. Rate is average of IUPAC MACR + O3 and 
! B139 MVK + O3. Split into further two reactions to allow all products 
! B139 to be included. Hence rates are IUPAC/4. Fractional Products from 
! B139 Poschl et al (2000). 
!---------------------------------------------------------------------- 
! B140 OClO+NO -> NO2 ClO  JPL2011   
! B140 Introduce T dependence of Bemand et al. New measurements by Li et 
! B140 al in good agreement with these 
!---------------------------------------------------------------------- 
! B141 OH+C2H6 -> H2O EtOO  IUPAC2007   
! B141 JPL slightly different (7.22 & 1020). Difference arises from 
! B141 fitting. Use IUPAC by default for organics. 
!---------------------------------------------------------------------- 
! B142 OH+C3H8 -> i-PrOO H2O  IUPAC2007 see also asad_bimol  
! B142 Total k rate (given) from IUPAC2007. Split between the two 
! B142 channels is temperature dependent using measurements of Droege 
! B142 and Tully (1986) and covered by asad_bimol. Note JPL values 
! B142 slightly higher 
!---------------------------------------------------------------------- 
! B143 OH+C3H8 -> n-PrOO H2O  IUPAC2007 see also asad_bimol  
! B143 Total k rate (given) from IUPAC2007. Split between the two 
! B143 channels is temperature dependent using measurements of Droege 
! B143 and Tully (1986) and covered by asad_bimol. Note JPL values 
! B143 slightly higher 
!---------------------------------------------------------------------- 
! B144 OH+C5H8 -> ISO2   IUPAC2009   
! B144 No JPL rate. 
!---------------------------------------------------------------------- 
! B145 OH+CH4 -> H2O MeOO  JPL2011   
! B145 Small difference with IUPAC results (1.85 and 1690) 
!---------------------------------------------------------------------- 
! B146 OH+CO -> HO2   IUPAC2005 see also asad_bimol  
! B146 Base rate modified by density dependence term in asad_bimol. JPL 
! B146 treat as a termolecular reaction. 
!---------------------------------------------------------------------- 
! B147 OH+ClO -> HCl O2  JPL2011   
! B147 - 
!---------------------------------------------------------------------- 
! B148 OH+ClO -> HO2 Cl  JPL2011   
! B148 - 
!---------------------------------------------------------------------- 
! B149 OH+ClONO2 -> HOCl NO3  JPL2011   
! B149 - 
!---------------------------------------------------------------------- 
! B150 OH+EtCHO -> H2O EtCO3  IUPAC2007   
! B150 Latest value includes measurements of Le Crane et al (2005)> 
!---------------------------------------------------------------------- 
! B151 OH+EtOOH -> H2O EtOO  MCMv3.2   
! B151 IUPAC provides a limit 
!---------------------------------------------------------------------- 
! B152 OH+EtOOH -> H2O MeCHO OH MCMv3.2   
! B152 IUPAC provides a limit 
!---------------------------------------------------------------------- 
! B153 OH+H2 -> H2O HO2  JPL2011   
! B153 IUPAC different parameterisation (7.7 and 2100). Same at 298 
!---------------------------------------------------------------------- 
! B154 OH+H2O2 -> H2O HO2  IUPAC2001   
! B154 - 
!---------------------------------------------------------------------- 
! B155 OH+HACET -> MGLY HO2  IUPAC2007   
! B155 Includes temperature dependence of Dillon et al (2006) 
!---------------------------------------------------------------------- 
! B156 OH+HBr -> H2O Br  JPL2011   
! B156 - 
!---------------------------------------------------------------------- 
! B157 OH+HCHO -> H2O HO2 CO IUPAC2007   
! B157 - 
!---------------------------------------------------------------------- 
! B158 OH+HCOOH -> HO2   IUPAC2007   
! B158 - 
!---------------------------------------------------------------------- 
! B159 OH+HCl -> H2O Cl  JPL2011   
! B159 - 
!---------------------------------------------------------------------- 
! B160 OH+HO2 -> H2O   JPL2011   
! B160 IUPAC identical 
!---------------------------------------------------------------------- 
! B161 OH+HO2NO2 -> H2O NO2  IUPAC2007   
! B161 The preferred values are based on the recent and extensive 
! B161 absolute rate study of JimÃÂ©nez et al. (2004) i 
!---------------------------------------------------------------------- 
! B162 OH+HOCl -> ClO H2O  JPL2011   
! B162 - 
!---------------------------------------------------------------------- 
! B163 OH+HONO -> H2O NO2  IUPAC2004   
! B163 IUPAC includes the more modern results of Burkholder. JPL values 
! B163 slightly different (1.8 and 390) 
!---------------------------------------------------------------------- 
! B164 OH+HONO2 -> H2O NO3  IUPAC2004 see also asad_bimol  
! B164 Include rate with no density dependence. Density dependence 
! B164 calculated using asad_bimol 
!---------------------------------------------------------------------- 
! B165 OH+ISON -> HACET NALD  Poschl00   
! B165 Use original MIM rate of Poschl et al (2000) 
!---------------------------------------------------------------------- 
! B166 OH+ISOOH -> MACR OH  Poschl00   
! B166 Use original MIM rate of Poschl et al (2000) 
!---------------------------------------------------------------------- 
! B167 OH+MACR -> MACRO2   IUPAC2007   
! B167 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC MVK 
! B167 rate. 
!---------------------------------------------------------------------- 
! B168 OH+MACR -> MACRO2   IUPAC2007   
! B168 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC 
! B168 MACR rate. 
!---------------------------------------------------------------------- 
! B169 OH+MACROOH -> MACRO2   MCMv3.2   
! B169  No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B170 OH+MGLY -> MeCO3 CO  IUPAC2008   
! B170 - 
!---------------------------------------------------------------------- 
! B171 OH+MPAN -> HACET NO2  IUPAC2006   
! B171 - 
!---------------------------------------------------------------------- 
! B172 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007   
! B172 Two part reaction due to IUPAC recommendation of the form k1+k2 
!---------------------------------------------------------------------- 
! B173 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007   
! B173  Two part reaction due to IUPAC recommendation of the form k1+k2 
!---------------------------------------------------------------------- 
! B174 OH+MeCHO -> H2O MeCO3  IUPAC2009   
! B174 Use total cross section. 
!---------------------------------------------------------------------- 
! B175 OH+MeCO2H -> MeOO   MCMv3.2   
! B175 No data in IUPAC or JPL. Rate k is based on MCM v3.2. Can't find 
! B175 the data used for previous versions of the code. 
!---------------------------------------------------------------------- 
! B176 OH+MeCO3H -> MeCO3   MCMv3.2   
! B176 No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B177 OH+MeCOCH2OOH -> H2O MeCOCH2OO  MCMv3.2   
! B177  No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B178 OH+MeCOCH2OOH -> OH MGLY  MCMv3.2   
! B178  No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B179 OH+MeOH -> HO2 HCHO  IUPAC2006   
! B179 JPL rates same. 
!---------------------------------------------------------------------- 
! B180 OH+MeONO2 -> HCHO NO2 H2O IUPAC2006   
! B180 JPL rates slightly different 8.0E-13 and 1000.  Large 
! B180 uncertainties on measurement. 
!---------------------------------------------------------------------- 
! B181 OH+MeOOH -> H2O HCHO OH IUPAC2007   
! B181 Use total rate from IUPAC * 0.4 (IUPAC BR to HCHO). Note JPL rate 
! B181 is considerably lower. There are discrepancies between measured 
! B181 values of this cross section which JPL and IUPAC reconcile 
! B181 differently. We employ the IUPAC values. 
!---------------------------------------------------------------------- 
! B182 OH+MeOOH -> H2O MeOO  IUPAC2007   
! B182 Use total rate from IUPAC * 0.6 (IUPAC BR to MeOO). Note JPL rate 
! B182 is considerably lower. 
!---------------------------------------------------------------------- 
! B183 OH+NALD -> HCHO CO NO2 IUPAC2009   
! B183 As in Poschl et al use MeCHO + OH rates (but updated). 
!---------------------------------------------------------------------- 
! B184 OH+NO3 -> HO2 NO2  JPL2011   
! B184 Note small discrepancy with IUPAC values (JPL rates 10% higher). 
! B184 Same measurements used 
!---------------------------------------------------------------------- 
! B185 OH+O3 -> HO2 O2  JPL2011   
! B185 Same values in IUPAC. 
!---------------------------------------------------------------------- 
! B186 OH+OClO -> HOCl O2  JPL2011   
! B186 Updated based on the recommended value reported by Gierczak et 
! B186 al. 
!---------------------------------------------------------------------- 
! B187 OH+OH -> H2O O(3P)  IUPAC2001   
! B187 No T dependence given by JPL 
!---------------------------------------------------------------------- 
! B188 OH+PAN -> HCHO NO2 H2O MCMv3.2   
! B188 No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B189 OH+PPAN -> MeCHO NO2 H2O MCMv3.2   
! B189 No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B190 OH+i-PrOOH -> Me2CO OH  MCMv3.2   
! B190 No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B191 OH+i-PrOOH -> i-PrOO H2O  MCMv3.2   
! B191 No data in IUPAC or JPL. Rate k is based on MCM v3.2 
!---------------------------------------------------------------------- 
! B192 OH+n-PrOOH -> EtCHO H2O OH MCMv3.2   
! B192 Rate k is based on MCM v3.2. IUPAC only provides total cross 
! B192 sections. 
!---------------------------------------------------------------------- 
! B193 OH+n-PrOOH -> n-PrOO H2O  MCMv3.2   
! B193 Rate k is based on MCM v3.2. IUPAC only provides total cross 
! B193 sections. 
!---------------------------------------------------------------------- 
! B194 i-PrOO+NO -> Me2CO HO2 NO2 IUPAC2005   
! B194 - 
!---------------------------------------------------------------------- 
! B195 i-PrOO+NO3 -> Me2CO HO2 NO2 MCMv3.2   
! B195 No data in IUPAC or JPL. Take rates from MCM v3.2. 
!---------------------------------------------------------------------- 
! B196 n-PrOO+NO -> EtCHO HO2 NO2 IUPAC2005   
! B196 The recommendation accepts the Arrhenius expression of Eberhard 
! B196 and Howard (1996).  Neglect n-propyl nitrate formation. 
!---------------------------------------------------------------------- 
! B197 n-PrOO+NO3 -> EtCHO HO2 NO2 MCM3.2   
! B197 No data in IUPAC or JPL. Take rates from MCM v3.2. 
!---------------------------------------------------------------------- 
 
ratb_defs_strattrop_chem=(/ &
  ratb_defs_strattrop01,ratb_defs_strattrop02,ratb_defs_strattrop03, &
  ratb_defs_strattrop04,ratb_defs_strattrop05 /)
  

depvel_defs_strattrop01=(/                  &
!  1  O3 (Ganzeveld & Lelieveld (1995) note 1 (modified to same as Guang)           
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  & !      1.1
  0.85,  0.30,  0.65,  0.65,  0.25,  0.45,  & !      1.2
  0.65,  0.25,  0.45,  0.65,  0.25,  0.45,  & !      1.3
  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  & !      1.4
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  & !      1.5
!  2  NO (inferred from NO2 - see Giannakopoulos (1998))                            
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !      2.1
  0.14,  0.01,  0.07,  0.01,  0.01,  0.01,  & !      2.2
  0.10,  0.01,  0.06,  0.01,  0.01,  0.01,  & !      2.3
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !      2.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !      2.5
!  3  NO3 (as NO2)                                                                  
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !      3.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,  & !      3.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,  & !      3.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !      3.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !      3.5
!  4  NO2                                                                           
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !      4.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,  & !      4.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,  & !      4.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !      4.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !      4.5
!  5  N2O5 (as HNO3)                                                                
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      5.1
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,  & !      5.2
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,  & !      5.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      5.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !      5.5
!  6  HO2NO2 (as HNO3)                                                              
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      6.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,  & !      6.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,  & !      6.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      6.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !      6.5
!  7  HNO3 (Zhang et al., 2003)                                                     
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      7.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,  & !      7.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,  & !      7.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      7.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !      7.5
!  8  H2O2                                                                          
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !      8.1
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,  & !      8.2
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,  & !      8.3
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  & !      8.4
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32,  & !      8.5
!  9  CO (see Giannokopoulos (1998))                                                
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !      9.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !      9.2
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !      9.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !      9.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !      9.5
! 10  HCHO (Zhang et al., 2003)                                                     
  1.00,  0.60,  0.80,  1.00,  0.60,  0.80,  & !     10.1
  1.60,  0.40,  1.00,  0.01,  0.01,  0.01,  & !     10.2
  0.71,  0.20,  0.45,  0.03,  0.03,  0.03,  & !     10.3
  0.40,  0.40,  0.40,  0.00,  0.00,  0.00,  & !     10.4
  0.30,  0.30,  0.30,  0.30,  0.30,  0.30,  & !     10.5
! 11  MeOOH                                                                         
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !     11.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & !     11.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & !     11.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     11.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     11.5
! 12  HCl (same as HBr)                                                             
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     12.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !     12.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !     12.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     12.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20   & !     12.5
  /)
  
depvel_defs_strattrop02=(/                  &
! 13  HOCl (same as HOBr)                                                           
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     13.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     13.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     13.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     13.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     13.5
! 14  HBr (note 3)                                                                  
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     14.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !     14.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !     14.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     14.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     14.5
! 15  HOBr (note 3)                                                                 
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     15.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     15.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & !     15.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & !     15.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     15.5
! 16  HONO                                                                          
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     16.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     16.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     16.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     16.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     16.5
! 17  EtOOH (as MeOOH)                                                              
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !     17.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & !     17.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & !     17.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     17.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     17.5
! 18  MeCHO                                                                         
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !     18.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  & !     18.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  & !     18.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     18.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     18.5
! 19  PAN (Zhang et al., 2003)                                                      
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     19.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & !     19.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & !     19.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     19.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  & !     19.5
! 20  n-PrOOH (as MeOOH)                                                            
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !     20.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & !     20.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & !     20.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     20.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     20.5
! 21  i-PrOOH (as MeOOH)                                                            
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !     21.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & !     21.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & !     21.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     21.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     21.5
! 22  EtCHO (as MeCHO)                                                              
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !     22.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  & !     22.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  & !     22.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     22.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     22.5
! 23  MeCOCH2OOH (as MeCO3H - where does the desert data come from!?)               
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !     23.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  & !     23.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  & !     23.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     23.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  & !     23.5
! 24  PPAN (as PAN)                                                                 
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     24.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & !     24.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & !     24.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     24.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00   & !     24.5
  /)
  
depvel_defs_strattrop03=(/                  &
! 25  ISOOH (as MeCO3H) (MIM)                                                       
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !     25.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  & !     25.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  & !     25.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     25.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  & !     25.5
! 26  ISON (as MeCO2H) (MIM) - note 2                                               
  0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & !     26.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  & !     26.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  & !     26.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  & !     26.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  & !     26.5
! 27  MACR (as MeCHO) (MIM)                                                         
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !     27.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  & !     27.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  & !     27.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     27.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     27.5
! 28  MACRO2H (as MeCO3H) (MIM)                                                     
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !     28.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  & !     28.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  & !     28.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     28.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  & !     28.5
! 29  MPAN (as PAN) (MIM)                                                           
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !     29.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,  & !     29.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,  & !     29.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  & !     29.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00,  & !     29.5
! 30  HACET (as MeCO3H) (MIM)                                                       
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !     30.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  & !     30.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  & !     30.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     30.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  & !     30.5
! 31  MGLY (as MeCHO) (MIM)                                                         
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !     31.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  & !     31.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  & !     31.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     31.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     31.5
! 32  NALD (as MeCHO) (MIM)                                                         
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !     32.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,  & !     32.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,  & !     32.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     32.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !     32.5
! 33  HCOOH (from Sander & Crutzen (1996) - note 3) (MIM)                           
  0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & !     33.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  & !     33.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  & !     33.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  & !     33.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  & !     33.5
! 34  MeCO3H (from Giannakopoulos (1998)) (MIM)                                     
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !     34.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,  & !     34.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,  & !     34.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     34.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,  & !     34.5
! 35  MeCO2H (as HCOOH) (MIM)                                                       
  0.50,  0.50,  0.50,  0.25,  0.25,  0.25,  & !     35.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,  & !     35.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,  & !     35.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,  & !     35.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03,  & !     35.5
! 36  MeOH (as MeOOH)                                                               
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !     36.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,  & !     36.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,  & !     36.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  & !     36.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01   & !     36.5
  /)
  

  depvel_defs_strattrop_chem=RESHAPE((/                                 &
                 depvel_defs_strattrop01,depvel_defs_strattrop02,       &
                 depvel_defs_strattrop03/),(/  6,  5, ndry_strattrop/))

  IF (L_ukca_het_psc) THEN
     ALLOCATE(rath_defs_strattrop_chem(nhet_strattrop))
     rath_defs_strattrop_chem = (/rath_defs_strattrop_psc/)
  ELSE
     ALLOCATE(rath_defs_strattrop_chem(0))
     rath_defs_strattrop_chem = (/rath_defs_strattrop_null/)
  END IF
     
 ierr = 0
 IF (L_ukca_strattrop) THEN
   IF (L_ukca_achem) THEN
!  Standard stratospheric chemistry with aerosol chemistry 
     IF (SIZE(chch_defs_strattrop_aer) /= jpspec) ierr = 1
     IF (SIZE(ratb_defs_strattrop_chem) + SIZE(ratb_defs_strattrop_aer) &
          /= jpbk) ierr = 2
     IF (SIZE(ratj_defs_strattrop_chem) + SIZE(ratj_defs_strattrop_aer) &
          /= jppj) ierr = 3
     IF (SIZE(ratt_defs_strattrop_chem) + SIZE(ratt_defs_strattrop_aer) &
          /= jptk) ierr = 4
     IF (L_ukca_het_psc) THEN
       IF (L_ukca_trophet) THEN
         IF (SIZE(rath_defs_strattrop_chem) +                           &
              SIZE(rath_defs_strattrop_aer) +                           &
              SIZE(rath_defs_strattrop_trophet)                         &
              /= jphk) ierr = 5
       ELSE
         IF (SIZE(rath_defs_strattrop_chem) +                           &
              SIZE(rath_defs_strattrop_aer)                             &
              /= jphk) ierr = 5
       END IF
     ELSE
       IF (L_ukca_trophet) THEN
         IF (SIZE(rath_defs_strattrop_aer) +                            &
              SIZE(rath_defs_strattrop_trophet)                         &
              /= jphk) ierr = 5
       ELSE
         IF (SIZE(rath_defs_strattrop_aer)                              &
              /= jphk) ierr = 5
       END IF
     END IF         
     IF (ndry_strattrop + ndry_st_aer /= jpdd) ierr = 6
     IF (nwet_st_aer /= jpdw) ierr = 7
     IF (ierr > 0) THEN
     cmessage = ' Size of chemical definition array incorrect'
     CALL EREPORT('UKCA_INIT_STRATTROP',ierr,cmessage)
   END IF
   chch_defs_strattrop = (/chch_defs_strattrop_aer/)
   ratb_defs_strattrop =                                                &
        (/ratb_defs_strattrop_chem, ratb_defs_strattrop_aer/)
   ratj_defs_strattrop =                                                &
        (/ratj_defs_strattrop_chem, ratj_defs_strattrop_aer/)
   ratt_defs_strattrop =                                                &
        (/ratt_defs_strattrop_chem, ratt_defs_strattrop_aer/)
   IF (L_ukca_het_psc) THEN
     IF (L_ukca_trophet) THEN
       rath_defs_strattrop =                                            &
              (/rath_defs_strattrop_chem,                               &
                rath_defs_strattrop_aer,                                &
                rath_defs_strattrop_trophet/)
     ELSE
       rath_defs_strattrop =                                            &
              (/rath_defs_strattrop_chem,                               &
                rath_defs_strattrop_aer/)
     END IF
   ELSE
     IF (L_ukca_trophet) THEN
       rath_defs_strattrop = (/rath_defs_strattrop_aer,               &
                               rath_defs_strattrop_trophet/)
     ELSE
       rath_defs_strattrop = (/rath_defs_strattrop_aer/)
     END IF
   END IF
   ! dry dep
   depvel_defs_strattrop(:,:,1:ndry_strattrop) =                        &
        depvel_defs_strattrop_chem(:,:,1:ndry_strattrop)
   depvel_defs_strattrop(:,:,                                           &
        ndry_strattrop+1:ndry_strattrop+ndry_st_aer) =                  &
                 depvel_defs_strattrop_aer(:,:,1:ndry_st_aer)
   ! wet dep
   henry_defs_strattrop(:,:) = henry_defs_strattrop_aer(:,1:nwet_st_aer)
   ELSE
! Standard stratospheric chemistry
   IF (SIZE(chch_defs_strattrop_chem) /= jpspec) ierr = 1
   IF (SIZE(ratb_defs_strattrop_chem) /= jpbk)   ierr = 2
   IF (SIZE(ratj_defs_strattrop_chem) /= jppj)   ierr = 3
   IF (SIZE(ratt_defs_strattrop_chem) /= jptk)   ierr = 4
   IF (L_ukca_het_psc) THEN
      IF (SIZE(rath_defs_strattrop_chem) /= jphk)   ierr = 5
   END IF
   IF (ndry_strattrop /= jpdd) ierr = 6
   IF (nwet_strattrop /= jpdw) ierr = 7
   IF (ierr > 0) THEN
     cmessage = ' Size of chemical definition array incorrect'
     CALL EREPORT('UKCA_INIT_STRATTROP',ierr,cmessage)
   END IF
   chch_defs_strattrop   = (/ chch_defs_strattrop_chem /)
   ratb_defs_strattrop   = (/ ratb_defs_strattrop_chem /)
   ratj_defs_strattrop   = (/ ratj_defs_strattrop_chem /)
   ratt_defs_strattrop   = (/ ratt_defs_strattrop_chem /)
   rath_defs_strattrop   = (/ rath_defs_strattrop_chem /)
   depvel_defs_strattrop =  depvel_defs_strattrop_chem
   henry_defs_strattrop  =   henry_defs_strattrop_chem
   ENDIF
 ELSE
   cmessage = 'No relevant definitions to compose chemistry'
   ierr = 1
   CALL EREPORT('UKCA_INIT_STRATTROP',ierr,cmessage)
 ENDIF

IF (lhook) CALL dr_hook('UKCA_INIT_STRATTROP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INIT_STRATTROP

END MODULE UKCA_CHEM_STRATTROP
