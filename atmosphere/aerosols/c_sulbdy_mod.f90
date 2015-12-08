! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_sulbdy_mod

IMPLICIT NONE
! Parameters for dry deposition of Sulphur Cycle tracers
! Rb(SO2)/Rb(H2O)
REAL,PARAMETER:: RESB_SO2=1.53    ! CUBRT(MolWt(SO2)/MolWt(H2O))

! Rb(NH3)/Rb(H2O)
REAL,PARAMETER:: RESB_NH3=0.981

! Rb(SO4 modes)/Rb(H2O)
REAL,PARAMETER:: RESB_SO4_AIT=94.9
REAL,PARAMETER:: RESB_SO4_ACC=2530.0
REAL,PARAMETER:: RESB_SO4_DIS=0.0

! Rs(SO2)/Rs(H2O)
REAL,PARAMETER:: RESS_SO2=1.89    ! SQRT(MolWt(SO2)/MolWt(H2O))

! Rs(NH3)/Rs(H2O)
REAL,PARAMETER:: RESS_NH3=0.972

! Rb(SO4 modes)/Rb(H2O)
REAL,PARAMETER:: RESS_SO4_AIT=0.0    ! Valid for small particles
REAL,PARAMETER:: RESS_SO4_ACC=0.0
REAL,PARAMETER:: RESS_SO4_DIS=0.0

! low limit for canopy conductance
REAL,PARAMETER:: COND_LIM=1.0E-3

! res to dry dep SO2 over snow
REAL,PARAMETER:: R_SNOW=1.0E3     ! s/m, from Pasdro et al 1993

! parameter for snow fraction calcn
REAL,PARAMETER:: ASNOW=0.2         ! m2/kg, from radiation doc.

! C_SULBDY end

END MODULE c_sulbdy_mod
