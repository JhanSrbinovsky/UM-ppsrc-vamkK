! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: Provide top level interface to model/application initialisation
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

MODULE UM_Config

  USE io_dependencies

  USE application_description, ONLY :         &
      exe_UM, exe_RCF, exe_scm,               &
      exe_pumf, exe_set_ancillary_time,       &
      exe_flux_transform, exe_convieee,       &
      exe_combine, exe_merge, exe_cumf,       &
      exe_fluxproc, exe_hreset, exe_fieldcos, &
      exe_fieldop, exe_fieldcalc, exe_convpp, &
      exe_setup, exe_fldmod, exe_hprint,      &
      exe_pptoanc, exe_pickup, exe_makebc,    &
      exe_frames, exe_vomext,                 &
      setApplicationDesc,                     &
      getExeType,                             &
      isParallel,                             &
      isSmallExec,                            &
      reportApplication

  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE appInit(exe)
    USE model_file, ONLY : model_file_init
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: exe
    
    CALL setApplicationDesc(exe)
    CALL model_file_init()
  END SUBROUTINE appInit
END MODULE UM_Config
