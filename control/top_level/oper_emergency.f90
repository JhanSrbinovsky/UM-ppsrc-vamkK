! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
SUBROUTINE oper_emergency

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE Submodel_Mod
      USE nlstcall_mod, ONLY : model_analysis_mins, &
                               run_target_end, &
                               run_assim_mode

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE
!
!    Method:  Reset UMUI determined namelist items to alternate
!             predetermined values dependent on value of Environment
!             variables known to the operational suite

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! System component covered: Control

! Declarations:

! Global variables (*CALLed COMDECKs etc...):
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"

! Print status information

! Local scalars:
CHARACTER(LEN=*) routinename
PARAMETER (   routinename='OperatorEmergency')
INTEGER icode            ! return code from FORT_GET_ENV
CHARACTER(LEN=80)  onlyto3    ! value of EV MES SHORT RUN
CHARACTER(LEN=80)  onlyto9    ! value of EV MOGREPS SHORT RUN
CHARACTER(LEN=80)  onlyto12   ! value of EV MOGREPS SHORT RUN
CHARACTER(LEN=80)  onlyto15   ! value of EV MOGREPS SHORT RUN
CHARACTER(LEN=80)  onlyto18   ! value of EV MOGREPS SHORT RUN
CHARACTER(LEN=80)  onlyto21   ! value of EV MOGREPS SHORT RUN
CHARACTER(LEN=80)  onlyto24   ! value of EV FOAM SHORT RUN
CHARACTER(LEN=80)  onlyto42   ! value of EV EXTENDED UK4
CHARACTER(LEN=80)  onlyto48   ! value of EV GL 2DAY RUN
CHARACTER(LEN=80)  onlyto54   ! value of EV MES ALTERNATE RUN
CHARACTER(LEN=80)  onlyto72   ! value of EV GL 3DAY RUN
CHARACTER(LEN=80)  fastrun    ! value of EV FASTRUN
CHARACTER(LEN=80)  shortstep  ! value of EV SHORTSTEP

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

! Function & Subroutine calls:
EXTERNAL fort_get_env

IF (lhook) CALL dr_hook('OPER_EMERGENCY',zhook_in,zhook_handle)

! Initialise all strings to false.
onlyto3   = "false"
onlyto9   = "false"
onlyto12  = "false"
onlyto15  = "false"
onlyto18  = "false"
onlyto21  = "false"
onlyto24  = "false"
onlyto42  = "false"
onlyto48  = "false"
onlyto54  = "false"
onlyto72  = "false"
fastrun   = "false"
shortstep = "false"

! FASTRUN
CALL fort_get_env('FASTRUN',7,FASTRUN,80,ICODE)
IF (fastrun == 'true'.AND.ICODE == 0)THEN
  run_assim_mode='FastRun   '
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "FASTRUN=true, setting run_assim_mode=",run_assim_mode
  END IF  ! PrintStatus test
END IF

! SHORTSTEP
CALL fort_get_env('SHORTSTEP',9,SHORTSTEP,80,ICODE)
IF (shortstep == 'true'.AND.ICODE == 0)THEN
  IF (printstatus >= prstatus_oper) THEN
    WRITE(6,*) routinename,':Warning, Wind limit re-setting '                  &
      ,'for SHORTSTEP no longer supported'
  END IF  ! PrintStatus test
END IF

IF (model_domain == mt_global) THEN
! UM6.5 : model_analysis_hrs replaced by model_analysis_mins
! ONLYTO72 for GL ATMOS plus 6 hour assm
CALL fort_get_env('ONLYTO72',8,ONLYTO72,80,ICODE)
IF (onlyto72 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                    ! reset all time
  run_target_end(3)=3                    ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0 +3 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*72)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO72=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

! ONLYTO48 for GL ATMOS
CALL fort_get_env('ONLYTO48',8,ONLYTO48,80,ICODE)
IF (onlyto48 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                  ! reset all time
  run_target_end(3)=2                  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*48)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO48=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

! ONLYTO21 for MOGREPS GL ATMOS
CALL fort_get_env('ONLYTO21',8,ONLYTO21,80,ICODE)
IF (onlyto21 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                  ! reset all time
  run_target_end(3)=0                  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+21 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*21)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO21=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

! ONLYTO18 for MOGREPS GL ATMOS
CALL fort_get_env('ONLYTO18',8,ONLYTO18,80,ICODE)
IF (onlyto18 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                  ! reset all time
  run_target_end(3)=0                  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+18 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*18)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO18=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

ELSE 
! ONLYTO3 for MES ATMOS plus 3 hours assm
CALL fort_get_env('ONLYTO3',7,ONLYTO3,80,ICODE)
IF (onlyto3 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0  ! reset all time
  run_target_end(3)=0  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+4 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*4)-                        &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO3=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

! ONLYTO42 for EXTENDED UK4 ATMOS
CALL fort_get_env('ONLYTO42',8,ONLYTO42,80,ICODE)
IF (onlyto42 == 'true' .AND. ICODE == 0) THEN
  run_target_end(:)=0  ! reset all time
  run_target_end(3)=1  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+18 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*42)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO42=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

! ONLYTO54 for ALTERNATE MES ATMOS
CALL fort_get_env('ONLYTO54',8,ONLYTO54,80,ICODE)
IF (onlyto54 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0  ! reset all time
  run_target_end(3)=2  ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+6 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*54)-                       &
    60.0*run_target_end(4)          ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO54=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF
END IF  ! if GLOBAL
!  ONLYTO9 for ATMOS
CALL fort_get_env('ONLYTO9',7,ONLYTO9,80,ICODE)
IF (onlyto9 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                    ! reset all time
  run_target_end(3)=0                    ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+9 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*9)-                        &
    60.0*run_target_end(4)           ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO9=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF
!  ONLYTO12 for ATMOS
CALL fort_get_env('ONLYTO12',8,ONLYTO12,80,ICODE)
IF (onlyto12 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                    ! reset all time
  run_target_end(3)=0                    ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+12 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*12)-                       &
    60.0*run_target_end(4)            ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO12=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF
!  ONLYTO15 for ATMOS
CALL fort_get_env('ONLYTO15',8,ONLYTO15,80,ICODE)
IF (onlyto15 == 'true'.AND.ICODE == 0)THEN
  run_target_end(:)=0                    ! reset all time
  run_target_end(3)=0                    ! days
  run_target_end(4)=REAL(model_analysis_mins)/60.0+15 ! hours
  run_target_end(5)=REAL(model_analysis_mins)+(60.0*15)-                       &
    60.0*run_target_end(4)            ! minutes
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0)WRITE(6,*)                                                   &
      "ONLYTO15=true, setting run_target_end=",run_target_end
  END IF  ! PrintStatus test
END IF

IF (lhook) CALL dr_hook('OPER_EMERGENCY',zhook_out,zhook_handle)
RETURN
END SUBROUTINE oper_emergency
