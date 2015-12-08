! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program: FIELDCOS 
!
! Purpose:
! To read a model dump or direct access fieldsfile and convert it to
! a sequential PP file ready for transfer to a different platform.
!
!  A general note on fieldcos -- When doing a bit compare on
!  the output of fieldcos half words may disagree. This is caused
!  by the extra half word after an odd number of words in a field
!  and is nothing to worry about. (Simon Tett 13/5/92)
!
!
! Programming standard: UM Doc Paper 3, version 1 
!
! Logical components covered: C41
!
! Project task: C4
!
! External documentation: UM documentation paper Y8
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      PROGRAM FIELDCOS
      USE IO
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
          
      USE check_iostat_mod
          
      USE ereport_mod, ONLY : ereport, ereport_finalise
      USE UM_Config, ONLY : &
          appInit, &
          exe_fieldcos
      USE pp_header_manips, ONLY : set_oper
      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

!     arguments for called routines
      CHARACTER                                                         &
     &     CMESSAGE*80                                                  &
                            ! Error message from lower routines
     &    ,FORMAT_OUT*6     ! IBM/IEEE for output format
                            ! GRIB - pure binary grib stash codes
                            ! GRIB1 - pure binary grib - standard codes
                            ! GRIB2 - pure binary grib - Other table 2
      CHARACTER(LEN=filenamelength) :: infile
                            ! Pathname of input file.
      CHARACTER(LEN=*) RoutineName
      PARAMETER (RoutineName='FIELDCOS')
      LOGICAL                                                           &
     &     UNPACK                                                       &
                                    ! indicates whether to unpack
     &    ,OPER                                                         &
                                    ! indicates whether operational
     &    ,MASS                     ! indicates whether its for MASS

      NAMELIST /PACK/ UNPACK,FORMAT_OUT
      NAMELIST /TYPE/ OPER, MASS
      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                                    ! First dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                    ! Size of the LOOKUP on the file
     &    ,PPUNIT                                                       &
                                    ! unit no of required fieldsfile
     &    ,COS_PPUNIT                                                   &
                                    ! unit no of COS output file
     &    ,IEXTRA(10)                                                   &
                                    ! spare for future use
     &    ,ICODE                                                        &
                                    ! return code
     &    ,DATA_ADD                                                     &
                                    ! The word address of the data.
     &    ,IWA                                                          &
                                    ! Word address in call SETPOS
     &    ,LEN_IO                                                       &
                                    ! Length of IO done
     &    ,LEN_FIXHD                ! Length of fixed length header
      PARAMETER(LEN_FIXHD=256)
      INTEGER                                                           &
     &     PP_FIXHD(LEN_FIXHD)      !  Fixed length header
      REAL                                                              &
     &     A_IO                     ! status returned by BUFFIN
      PARAMETER(LEN1_LOOKUP=64)
      DATA UNPACK/.FALSE./
      DATA FORMAT_OUT/'IBM   '/
      DATA OPER/.FALSE./
      DATA MASS/.FALSE./
! defines NUNITS
! CGRIBTAB start
!
! Description:
!  Holds a lookup table for converting between Unified Model
! stash codes and some other grib table of codes
!
! Declarations:
      INTEGER, PARAMETER:: MAX_SECT_GRBTAB=16
      INTEGER, PARAMETER:: MAX_ITEM_GRBTAB=300

      INTEGER :: GRIB_TABLE(0:MAX_SECT_GRBTAB,MAX_ITEM_GRBTAB)

      COMMON /GRIBTAB/ GRIB_TABLE

! CGRIBTAB end
      LOGICAL FLAG                  ! =T/F file exists/not
      COMMON /FLAG_IO/FLAG(NUNITS)  ! needed for BUFFIN check
!*---------------------------------------------------------------------
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,IX                                                           &
                                  ! used as a dummy variable in UNIT
     &    ,ERR                                                          &
     &    ,DIAG_UNIT
     
      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error     
      INTEGER :: me_gc
      INTEGER :: nproc_gc
      LOGICAL LCAL360   ! 360 day calendar switch
!  Initialise LCAL360
      DATA LCAL360 /.FALSE./
      CHARACTER(LEN=filenamelength) :: diagfile

      CALL gc_init(' ',me_gc,nproc_gc)
      CALL appInit(exe_fieldcos)
      CALL ioInit()

!=====================================================================
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!=====================================================================
      READ (UNIT=5, NML=PACK, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist PACK")
      READ (UNIT=5, NML=TYPE, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist TYPE")

      CALL set_oper(oper)

      IF (FORMAT_OUT /= 'IEEE' .AND. MASS) THEN
        
        ICODE    = 1
        CMESSAGE = 'Can only use IEEE format with MASS.'
        CALL EREPORT(RoutineName,ICODE,CMESSAGE)
      END IF

      WRITE(6,*)'  UNPACK  ',UNPACK
      WRITE(6,*)'  FORMAT  ',FORMAT_OUT
      WRITE(6,*)'  OPER    ',OPER
      WRITE(6,*)'  MASS    ',MASS
      DO I=1,10
        IEXTRA(I)=0
      ENDDO

      DIAG_UNIT = 7
      CALL GET_FILE(DIAG_UNIT,DIAGFILE,filenamelength,ICODE)
      OPEN(UNIT=DIAG_UNIT,FILE=DIAGFILE)

      PPUNIT=10
      COS_PPUNIT=11
! -------------------------------------------------------------------
! If FORMAT_OUT is GRIB1 or GRIB2 initialise grib field code
! conversion table
      IF (FORMAT_OUT == 'GRIB1') THEN
! DEPENDS ON: grib_table_init1
        CALL GRIB_TABLE_INIT1
      ELSE IF (FORMAT_OUT == 'GRIB2') THEN
! DEPENDS ON: grib_table_init2
        CALL GRIB_TABLE_INIT2
      ENDIF
!L-------------Read in the FIXED length header------------------------
      CALL GET_FILE(PPUNIT,INFILE,filenamelength,ICODE)
      CALL FILE_OPEN(PPUNIT,INFILE,filenamelength,0,1,ERR)
      FLAG(PPUNIT)=.TRUE.          ! needed for BUFFIN check
      CALL BUFFIN(PPUNIT,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,     &
     &                  LEN_FIXHD)
          CMESSAGE='FIELDCOS : I/O error reading FIXED LENGTH HEADER'
          ICODE=2
          WRITE(6,*)' I/O error reading FIXED LENGTH HEADER'
! DEPENDS ON: abort_io
          CALL ABORT_IO('FIELDCOS',CMESSAGE,ICODE,PPUNIT)
      ENDIF
      DATA_ADD=PP_FIXHD(160)-1 ! Start address for the data.
      IWA= PP_FIXHD(150)-1     ! Start address for the lookup table.
      PP_LEN2_LOOKUP=PP_FIXHD(152)
      WRITE(6,*)' PP_LEN2_LOOKUP  ',PP_LEN2_LOOKUP
      WRITE(6,*)' dump type=',pp_fixhd(5),                              &
     &       ' 3=fieldsfile,1=dump,2=time mean dump,4=ancil,5=bound'
! DEPENDS ON: read_write_fieldcos
      CALL READ_WRITE_fieldcos(PP_LEN2_LOOKUP,LEN1_LOOKUP,DATA_ADD,     &
     &                PP_FIXHD,                                         &
     &                IWA,UNPACK,FORMAT_OUT,PPUNIT,COS_PPUNIT,          &
     &                IEXTRA,OPER,MASS,ICODE,CMESSAGE,LCAL360)
      IF(ICODE /= 0) THEN

        CALL EREPORT(RoutineName,ICODE,CMESSAGE)
      ENDIF

      CALL ioShutdown()

      CALL ereport_finalise( )

      END PROGRAM FIELDCOS

