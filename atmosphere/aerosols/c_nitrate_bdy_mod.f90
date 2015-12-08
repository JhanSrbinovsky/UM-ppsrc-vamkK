! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE c_nitrate_bdy_mod

IMPLICIT NONE
!----------------COMDECK C_NITRATE_BDY------------------------------------
! Parameters for dry deposition of ammonium nitrate tracers.
! Use same vales as for sulphate for the moment.
!
      ! Rb(NH4NO3 modes)/Rb(H2O)
      REAL,PARAMETER:: RESB_NITR_ACC=2530.0
      REAL,PARAMETER:: RESB_NITR_DISS=0.0
!
      ! Rs(NH4NO3 modes)/Rs(H2O)
      REAL,PARAMETER:: RESS_NITR_ACC=0.0
      REAL,PARAMETER:: RESS_NITR_DISS=0.0
!
! C_NITRATE_BDY end

END MODULE c_nitrate_bdy_mod
