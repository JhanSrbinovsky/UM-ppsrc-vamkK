! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Defines namelist containing actually used UM sections, these
! values will be passed to atmos_sr and indep_sr arrays
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards. 


MODULE umsections_mod

   USE version_mod, ONLY: nsectp

   IMPLICIT NONE

!  UM code sections
   CHARACTER (LEN = 2), SAVE  ::  atmos_sr(0:nsectp)
   CHARACTER (LEN = 2), SAVE  ::  indep_sr(0:nsectp)

!  ROSE Atmos sections
   CHARACTER(LEN=2) :: Section_A01 = '  ' ! Radiation SW
   CHARACTER(LEN=2) :: Section_A02 = '  ' ! Radiation LW
   CHARACTER(LEN=2) :: Section_A03 = '  ' ! Boundary Layer
   CHARACTER(LEN=2) :: Section_A04 = '  ' ! LS Precipitation
   CHARACTER(LEN=2) :: Section_A05 = '  ' ! Convection
   CHARACTER(LEN=2) :: Section_A06 = '  ' ! Gravity Wave Drag

   CHARACTER(LEN=2) :: Section_A08 = '  ' ! Hydrology
   CHARACTER(LEN=2) :: Section_A09 = '  ' ! LS Cloud
   CHARACTER(LEN=2) :: Section_A10 = '  ' ! Dynamical Solver
   CHARACTER(LEN=2) :: Section_A11 = '  ' ! Atm Tracer Advection
   CHARACTER(LEN=2) :: Section_A12 = '  ' ! Primary Field Advection
   CHARACTER(LEN=2) :: Section_A13 = '  ' ! Diffusion, Filtering, Moisture
   CHARACTER(LEN=2) :: Section_A14 = '  ' ! Energy Correction
   CHARACTER(LEN=2) :: Section_A15 = '  ' ! Dyn Diag
   CHARACTER(LEN=2) :: Section_A16 = '  ' ! Phy Diag
   CHARACTER(LEN=2) :: Section_A17 = '  ' ! Aerosols
   CHARACTER(LEN=2) :: Section_A18 = '  ' ! Data Assimilation
   CHARACTER(LEN=2) :: Section_A19 = '  ' ! Vegetation Distribution

   CHARACTER(LEN=2) :: Section_A26 = '  ' ! River Routing

   CHARACTER(LEN=2) :: Section_A30 = '  ' ! Climatological Diag
   CHARACTER(LEN=2) :: Section_A31 = '  ' ! Other Ancil LBC
   CHARACTER(LEN=2) :: Section_A32 = '  ' ! Other Ancil LBC
   CHARACTER(LEN=2) :: Section_A33 = '  ' ! Atmos Config Tracer
   CHARACTER(LEN=2) :: Section_A34 = '  ' ! UKCA Prognostics
   CHARACTER(LEN=2) :: Section_A35 = '  ' ! Stochastic
   CHARACTER(LEN=2) :: Section_A36 = '  ' ! LBC Tracers
   CHARACTER(LEN=2) :: Section_A37 = '  ' ! UKCA LBC
   CHARACTER(LEN=2) :: Section_A38 = '  ' ! UKCA 2
   CHARACTER(LEN=2) :: Section_A39 = '  ' ! Nudging
   CHARACTER(LEN=2) :: Section_A50 = '  ' ! UKCA Diagnostics

   CHARACTER(LEN=2) :: Section_A71 = '  ' ! Ozone

!Independent sections
   CHARACTER(LEN=2) :: Section_C70 = '  ' ! cntl lev routine for internal models
   CHARACTER(LEN=2) :: Section_C72 = '  ' ! cntl lev routine for coupled models
   CHARACTER(LEN=2) :: Section_C80 = '  ' ! dump I/O 
   CHARACTER(LEN=2) :: Section_C82 = '  ' ! ancil files sesction
   CHARACTER(LEN=2) :: Section_C84 = '  ' ! stash 
   CHARACTER(LEN=2) :: Section_C92 = '  ' ! limited area and ancil service 
   CHARACTER(LEN=2) :: Section_C94 = '  ' ! further service routine
   CHARACTER(LEN=2) :: Section_C95 = '  ' ! portable I/O
   CHARACTER(LEN=2) :: Section_C96 = '  ' ! MPP service
   CHARACTER(LEN=2) :: Section_C97 = '  ' ! timer code
   CHARACTER(LEN=2) :: Section_C98 = '  ' ! openMP

   NAMELIST/UMSECTIONS/                                                 &
      Section_A01, Section_A02, Section_A03, Section_A04, Section_A05,  &
      Section_A06, Section_A08, Section_A09, Section_A10, Section_A11,  &
      Section_A12, Section_A13, Section_A14, Section_A15, Section_A16,  &
      Section_A17, Section_A18, Section_A19, Section_A26,               &
      Section_A30, Section_A31, Section_A32, Section_A33, Section_A34,  &
      Section_A35, Section_A36, Section_A37, Section_A38, Section_A39,  &
      Section_A50, Section_A71,                                         &
      Section_C70, Section_C72, Section_C80, Section_C82, Section_C84,  &
      Section_C92, Section_C94, Section_C95, Section_C96, Section_C97,  &
      Section_C98

   CONTAINS
   SUBROUTINE assign_umsections 
      IMPLICIT NONE

! Local counter
      INTEGER :: i

! Initialise model config. arrays 
      DO i = 0, nsectp
         atmos_sr(i)='  '
         indep_sr(i)='  '
      END DO

      atmos_sr(1) = Section_A01
      atmos_sr(2) = Section_A02
      atmos_sr(3) = Section_A03
      atmos_sr(4) = Section_A04
      atmos_sr(5) = Section_A05
      atmos_sr(6) = Section_A06  
      atmos_sr(8) = Section_A08
      atmos_sr(9) = Section_A09
      atmos_sr(10) = Section_A10
      atmos_sr(11) = Section_A11
      atmos_sr(12) = Section_A12
      atmos_sr(13) = Section_A13
      atmos_sr(14) = Section_A14
      atmos_sr(15) = Section_A15
      atmos_sr(16) = Section_A16
      atmos_sr(17) = Section_A17
      atmos_sr(18) = Section_A18
      atmos_sr(19) = Section_A19
      atmos_sr(26) = Section_A26
      atmos_sr(30) = Section_A30
      atmos_sr(31) = Section_A31
      atmos_sr(32) = Section_A32
      atmos_sr(33) = Section_A33
      atmos_sr(34) = Section_A34
      atmos_sr(35) = Section_A35
      atmos_sr(36) = Section_A36
      atmos_sr(37) = Section_A37
      atmos_sr(38) = Section_A38
      atmos_sr(39) = Section_A39
      atmos_sr(50) = Section_A50
      atmos_sr(71) = Section_A71

      indep_sr(70)= Section_C70
      indep_sr(72)= Section_C72
      indep_sr(80)= Section_C82
      indep_sr(82)= Section_C82
      indep_sr(84)= Section_C84
      indep_sr(92)= Section_C92
      indep_sr(94)= Section_C94
      indep_sr(95)= Section_C95
      indep_sr(96)= Section_C96
      indep_sr(97)= Section_C97
      indep_sr(98)= Section_C98

   END SUBROUTINE assign_umsections
END MODULE umsections_mod
