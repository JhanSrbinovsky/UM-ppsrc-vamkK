! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE nudging_input_mod

! Description:
!       Input/namelist control for Nudging scheme.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi,imdi

IMPLICIT NONE

!-------------------- Nudging control parameters ----------------------
!
      REAL      :: ndg_relax_uvalue = rmdi
      REAL      :: ndg_relax_vvalue = rmdi
      REAL      :: ndg_relax_tvalue = rmdi
                 ! Relaxation parameters for u, v, T

      INTEGER   :: ndg_lev_bottom = imdi
      INTEGER   :: ndg_lev_top    = imdi
                 ! Levels that we nudge upon from lowest to highest

      INTEGER   :: ndg_on_lev_bottom = imdi
      INTEGER   :: ndg_on_lev_top    = imdi
                 ! Number of levels to go from none to full strength nudging.

      CHARACTER(LEN=256)   :: ndg_datapath = ' '
                            ! Path to analyses files
      INTEGER   :: ndg_hours_perdata = imdi
                 ! Numbers of hours per analyses data timestep
      INTEGER   :: ndg_analysis_source = imdi
                                  ! Source of Analyses data being used
                                  !  0 = ECMWF on hybrid levels,
                                  !  1 = ECMWF on pressure levels,
                                  !  2 = UM analyses
                                  !  3 = JRA Analyses
      REAL      :: ndg_strat_fac = rmdi
                 ! Factor to weaken nudging in stratosphere
      LOGICAL   :: l_nudging        
                 ! True when Nudging is switched on

      NAMELIST /Run_Nudging/                                            &
        ndg_relax_uvalue,ndg_relax_vvalue,ndg_relax_tvalue,             &
        ndg_lev_bottom,ndg_lev_top,ndg_on_lev_bottom,ndg_on_lev_top,    &
        ndg_datapath,ndg_hours_perdata,ndg_analysis_source,             &
        ndg_strat_fac, l_nudging


END MODULE nudging_input_mod
