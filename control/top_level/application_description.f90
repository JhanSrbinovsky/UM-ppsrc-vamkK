! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: Application_Description -----------------------------------
!
!  Purpose: Provides runtime registry and enquiry functions allowing
!           UM code to determine the executable target (program main),
!           build type, sizes, and parallel nature. This facilitates
!           reduction in preprocessor based compilation in favour of
!           software switching
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

MODULE Application_Description

  USE ereport_mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: exe_unknown           = 1
  INTEGER, PARAMETER :: exe_UM                = 2
  INTEGER, PARAMETER :: exe_RCF               = 3
  INTEGER, PARAMETER :: exe_scm               = 4
  INTEGER, PARAMETER :: exe_pumf              = 5
  INTEGER, PARAMETER :: exe_set_ancillary_time= 6
  INTEGER, PARAMETER :: exe_flux_transform    = 7
  INTEGER, PARAMETER :: exe_convieee          = 8
  INTEGER, PARAMETER :: exe_combine           = 9
  INTEGER, PARAMETER :: exe_merge             = 10
  INTEGER, PARAMETER :: exe_cumf              = 11
  INTEGER, PARAMETER :: exe_fluxproc          = 12
  INTEGER, PARAMETER :: exe_hreset            = 13
  INTEGER, PARAMETER :: exe_fieldcos          = 14
  INTEGER, PARAMETER :: exe_fieldop           = 15
  INTEGER, PARAMETER :: exe_fieldcalc         = 16
  INTEGER, PARAMETER :: exe_convpp            = 17
  INTEGER, PARAMETER :: exe_setup             = 18
  INTEGER, PARAMETER :: exe_fldmod            = 19
  INTEGER, PARAMETER :: exe_hprint            = 20
  INTEGER, PARAMETER :: exe_pptoanc           = 21
  INTEGER, PARAMETER :: exe_pickup            = 22
  INTEGER, PARAMETER :: exe_makebc            = 23
  INTEGER, PARAMETER :: exe_frames            = 24
  INTEGER, PARAMETER :: exe_vomext            = 25
  INTEGER, PARAMETER :: num_exe_types         = 25

  CHARACTER (LEN=*),PARAMETER::&
      UM_name=&
      "Application Unknown     "//&
      "Unified Model           "//&
      "Reconfiguration         "//&
      "Single Column Model     "//&
      "Print UM File           "//&
      "Set Ancilary Time       "//&
      "Flux Transform          "//&
      "Convert IEEE            "//&
      "Combine                 "//&
      "Merge                   "//&
      "Compare UM Files        "//&
      "Flux Process            "//&
      "HReset                  "//&
      "Field Cos               "//&
      "Field Op                "//&
      "Field Calc              "//&
      "Convert PP              "//&
      "Setup                   "//&
      "Field Mod               "//&
      "HPrint                  "//&
      "PP to Ancilary          "//&
      "Pickup                  "//&
      "Make Boundary Conditions"//&
      "Frames                  "//&
      "Vom Extract             "
  !    123456789012345678901234 <- 24

  ! This is the length of each string subsection above
  INTEGER, PARAMETER  :: UM_Name_Len = 24

  INTEGER :: exe_type         =exe_unknown
  LOGICAL :: exe_is_parallel  =.FALSE.
  LOGICAL :: exe_addr_64      =.FALSE.
  LOGICAL :: exe_data_64      =.FALSE.
  
CONTAINS
  
  SUBROUTINE        setApplicationDesc(appType) ! Call once from main.
    USE UM_Types
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: appType
    INTEGER             :: rtype=4
    INTEGER             :: itype=4
    INTEGER             :: ptype=4
    INTEGER             :: errorCode
    CHARACTER(LEN=132)  :: message

    errorCode=1

    IF (appType>0 .AND. appType <= num_exe_types ) THEN
      exe_type=appType  
    ELSE 
      WRITE(message,'(A,I3)') &
          'Tried to set application to unknown value:',appType
      CALL ereport('applcation_description:setApplicationDescription', &
          errorCode,message)
    END IF
    
    rtype= umFortranRealSize()
    itype= umFortranIntegerSize()    
    ptype= umFortranPointerSize()

    IF (rtype==8 .AND. itype==8) THEN
      exe_data_64=.TRUE.
    ELSE IF (rtype==4 .AND. itype==4) THEN
      exe_data_64=.FALSE.
    ELSE
      WRITE(message,'(A)') &
          'Application has different lengths for reals and integers'
      CALL ereport('applcation_description:setApplicationDescription', &
          errorCode,message)
    END IF

    IF (ptype==8) THEN
      exe_addr_64=.TRUE.
    ELSE IF (ptype==4) THEN
      exe_data_64=.FALSE.
    ELSE
      WRITE(message,'(A)') &
          'Application has an unknown pointer size'
      CALL ereport('applcation_description:setApplicationDescription', &
          errorCode,message)
    END IF

    CALL reportApplication()

  END SUBROUTINE setApplicationDesc

  ! Provide the application type
  FUNCTION getExeType() RESULT(r)      
    IMPLICIT NONE
    INTEGER :: r
    r=exe_type
  END FUNCTION getExeType

  ! Is this an MPI program (ie I can call MPL, rather than nproc>1)
  FUNCTION isParallel() RESULT(r)
    IMPLICIT NONE
    LOGICAL :: r
    INTEGER :: istat
    INTEGER :: val   
    INTEGER :: errorCode
    CHARACTER(LEN=132)  :: message
    INTEGER, PARAMETER :: GC_PARALLEL=1
    INTEGER, PARAMETER :: GC_SERIAL=2
    INTEGER, PARAMETER :: GC_IS_PARALLEL=3
    CALL GC_GetOpt(GC_IS_PARALLEL,val,istat)
    IF (val==GC_SERIAL) THEN
      r=.FALSE.
    ELSE IF (val==GC_PARALLEL) THEN
      r=.TRUE.
    ELSE
      errorCode=2
      WRITE(message,'(A,I3)') &
          'GCOM returned an unknown value for the GC_IS_PARALLEL check:', &
          val
      CALL ereport('applcation_description:isParallel', &
          errorCode,message)
    END IF
  END FUNCTION isParallel

  ! Convenience routines, to taste.
  FUNCTION isSmallExec() RESULT(r)
    IMPLICIT NONE 
    LOGICAL :: r
    r=.FALSE.
    IF (exe_type==1) THEN
      r=.TRUE.
    ELSE IF (exe_type>2) THEN
      r=.TRUE.
    END IF

  END FUNCTION isSmallExec

  SUBROUTINE reportApplication()
    IMPLICIT NONE
    INTEGER           :: i
    CHARACTER(LEN=2)  :: ds
    CHARACTER(LEN=2)  :: ps
    WRITE(6,'(A)')' ******************************************'
    WRITE(6,'(A,I2,A,A)')    &
        ' App ID:',exe_type, &
        ', Name:',TRIM(UM_NameStr(exe_type))
    WRITE(6,'(A)')' ------------------------------------------'

    WRITE(ds,'(A)')'32'
    WRITE(ps,'(A)')'32'
    IF (exe_data_64) WRITE(ds,'(A)')'64'
    IF (exe_addr_64) WRITE(ps,'(A)')'64'
    WRITE(6,'(A,A,A,A,A)')'  - Data size is ',ds, &
        ' bit. Program is ',ps,' bit.'
    IF (isParallel()) THEN
      WRITE(6,'(A,I2)')'  - Program is parallel.'
    ELSE
      WRITE(6,'(A,I2)')'  - Program is serial.'
    END IF
    WRITE(6,'(A)')' ******************************************'
    
  END SUBROUTINE reportApplication

  FUNCTION UM_NameStr(i) RESULT(r)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: i
    CHARACTER(LEN=UM_Name_Len)    :: r    
    WRITE(r,'(A)')UM_name((exe_type-1)*UM_Name_Len+1:i*UM_Name_len)
  END FUNCTION UM_NameStr

END MODULE Application_Description
