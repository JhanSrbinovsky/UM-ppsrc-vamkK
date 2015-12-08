! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ data module defining RECON namelist

MODULE Rcf_Recon_Mod

! Description:
!   Defines variables for the RECON namelist.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


USE Submodel_Mod, ONLY :      &
    Submodel_Ident

USE missing_data_mod, ONLY:   & 
    IMDI

IMPLICIT NONE

INTEGER       :: DUMP_PACK
INTEGER       :: TPPS_OZONE_LEVELS
INTEGER       :: CONV_LEVELS
INTEGER       :: W_ZERO_START
INTEGER       :: W_ZERO_END

INTEGER       :: LEN_FIXHD     = IMDI
INTEGER       :: LEN_DUMPHIST  = IMDI

REAL          :: Q_MIN
REAL          :: PERTURBATION

LOGICAL       :: GRIB
LOGICAL       :: GRIB2FF
LOGICAL       :: UARS
LOGICAL       :: RESET
LOGICAL       :: TRANS
LOGICAL       :: LSPIRAL_S
LOGICAL       :: Var_Recon
LOGICAL       :: USE_SMC_STRESS
LOGICAL       :: polar_check
LOGICAL       :: l_canopy_snow_throughfall = .FALSE.
LOGICAL       :: l_interp_input_only = .FALSE.

! Namelist RECON.

NAMELIST /RECON/                                                  &
 SUBMODEL_IDENT, DUMP_PACK, CONV_LEVELS,                          &
 LEN_FIXHD, LEN_DUMPHIST, GRIB, UARS, RESET,                      &
 PERTURBATION, TRANS,  LSPIRAL_S,                                 &
 Var_Recon, Q_Min, W_Zero_Start, W_Zero_End,                      &
 Use_SMC_Stress, polar_check, grib2ff, l_canopy_snow_throughfall, &
 l_interp_input_only


END MODULE Rcf_Recon_Mod
