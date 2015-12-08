! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the STSHCOMP namelist

MODULE Rcf_readnl_stshcomp_Mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_stshcomp - reads the STSHCOMP namelist.
!
! Description:
!   Read the STSHCOMP namelist for controlling section choices etc.
!
! Method:
!   Variables read into Rcf_Model_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

!+  Read in the STSHCOMP namelists.

SUBROUTINE rcf_readnl_stshcomp (nft)

USE Rcf_Model_Mod, ONLY : &
    A_Max_TrVars,                                     &
    A_Max_UKCAVars,                                   &
    tca,                                              &
    tca_lbc,                                          &
    tc_lbc_ukca,                                      &
    ZonAvTppsOzone,                                   &
    STSHCOMP,                                         &
    meso,                                             &
    ZonAvOzone,                                       &
    ewspacea,     nsspacea,                           &
    frstlona,     frstlata,                           &
    polelona,     polelata,                           &
    delta_lon,    delta_lat
    

USE PrintStatus_mod, ONLY : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, ONLY : &
    mype

USE missing_data_mod, ONLY : &
    RMDI

USE rcf_grid_type_mod, ONLY:  &
    output_grid
USE nlsizes_namelist_mod, ONLY: &
    model_levels, tr_levels, sm_levels, st_levels, tr_vars

USE umsections_mod  ! All sections in use

IMPLICIT NONE

INTEGER nft
INTEGER i
INTEGER tr_lbc_vars
INTEGER tr_lbc_ukca

! Set defaults
MESO   = ' '
ewspacea =  RMDI
nsspacea =  RMDI
frstlona =   0.0
frstlata = -90.0
polelona =   0.0
polelata =  90.0

READ(UNIT = nft, NML = UMSECTIONS)
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(UNIT = 6, NML = UMSECTIONS)
END IF

! Assign UM sections in ROSE format to atmos_sr and indep_sr
 CALL assign_umsections

READ(UNIT = nft, NML = STSHCOMP)
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(UNIT = 6, NML = STSHCOMP)
END IF

! Give some values more user-friendly names
delta_lon = ewspacea
delta_lat = nsspacea

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
 ZonAvTppsOzone=ZonAvOzone

! Set the number of soil moisture levels to zero unless the hydrology scheme
! (section 8) in use is '7A' - Moses-II or '8A' - Jules, then set it to be the
! same as the number of soil temperature levels
sm_levels=0
IF (Atmos_SR(8) == '7A' .OR. Atmos_SR(8) == '8A')THEN
  sm_levels=st_levels
END IF
output_grid % sm_levels = sm_levels

! If atmospheric tracers (section 33) are being used (i.e. not set to '0A'),
! then count the number of tracers actually being used
tr_vars=0
IF (Atmos_SR(33) /= '0A')THEN
  DO i=1,A_Max_TrVars
    IF(tca(i) /= 0)tr_vars = tr_vars + 1
  END DO
END IF

! Find the number of tracers being used by counting the non-zero elements
! of the arrays that list them
tr_lbc_vars=0
tr_lbc_ukca=0
DO i=1,A_Max_TrVars
  IF(tca_lbc(i) /= 0) tr_lbc_vars = tr_lbc_vars + 1
END DO
DO i=1,A_Max_UKCAVars
  IF(tc_lbc_ukca(i) /= 0) tr_lbc_ukca = tr_lbc_ukca + 1
END DO

! If no tracers have been selected and if aerosols (section 17) aren't being
! used and if UKCA (section 34) isn't being used then reset the number of
! tracer levels to 1 - otherwise, use model_levels
IF(tr_vars==0 .AND. tr_lbc_vars==0 .AND. tr_lbc_ukca==0              &
    .AND. atmos_sr(17)=="0A" .AND. atmos_sr(34)/="1A")THEN
  tr_levels=1
ELSE
  tr_levels=model_levels
END IF
output_grid % tr_levels         = tr_levels

RETURN

END SUBROUTINE  Rcf_readnl_stshcomp

END MODULE Rcf_readnl_stshcomp_Mod
