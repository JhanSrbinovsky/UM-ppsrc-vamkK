! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: FLDMOD --------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      PROGRAM FLDMOD
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE IO
      USE ereport_mod, ONLY: ereport, ereport_finalise
      USE UM_Config, ONLY : &
          appInit, &
          exe_fldmod
      IMPLICIT NONE 
      INTEGER :: me_gc
      INTEGER :: nproc_gc
      
      EXTERNAL DIMENS1
      CHARACTER CMESSAGE*80
      CHARACTER(LEN=filenamelength) :: diagfile
      CHARACTER(LEN=filenamelength) :: infile
      CHARACTER(LEN=8) c_nproc            ! to get nproc_x and nproc_y from
!                                    ! environment variables.
!                         up to an EVEN no for conversion to IBM format
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
                                  !    Length of fixed length header
     &    ,LEN_INTHD                                                    &
     &    ,LEN_REALHD                                                   &
     &    ,LEN1_LEVDPC                                                  &
     &    ,LEN2_LEVDPC                                                  &
     &    ,LEN1_ROWDPC                                                  &
     &    ,LEN2_ROWDPC                                                  &     
     &    ,LEN1_COLDPC                                                  &
     &    ,LEN2_COLDPC                                                  &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,PPUNIT1                                                      &
                                  !OUT unit no of required fieldsfile 1
     &    ,PPUNIT2                                                      &
                                  !OUT unit no of required fieldsfile 2
     &    ,DIAG_UNIT                                                    &
                                  !unit number for diagnostics
     &    ,ICODE                                                        &
                                  !IN  return code
     &    ,ERR
!    LOCAL VARIABLES
      PARAMETER(LEN_FIXHD=256)
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
                                  !IN  Fixed length header
     &    ,IWA                                                          &
                                  !
     &    ,IX                                                           &
                                  !
     &    ,LEN_IO                 !
      REAL                                                              &
     &     A_IO                   !
!

      CHARACTER(LEN=*),PARAMETER :: RoutineName = 'fldmod'

      CALL gc_init(' ',me_gc,nproc_gc)
      CALL appInit(exe_fldmod)
      CALL ioInit()

!    OPEN DIAGNOSTIC FILE
      DIAG_UNIT = 7
      CALL GET_FILE(DIAG_UNIT,DIAGFILE,filenamelength,ICODE)
      OPEN(UNIT=DIAG_UNIT,FILE=DIAGFILE)

!*****************************************************************
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
      icode = 0
      cmessage= '                                         '
!
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
!*****************************************************************
      PPUNIT1=10
      PPUNIT2=11
!*****************************************************************
!     Buffer in the Fixed Length Header and obtain lengths
!*****************************************************************
      CALL GET_FILE(PPUNIT1,INFILE,filenamelength,ICODE)
      CALL FILE_OPEN(PPUNIT1,INFILE,filenamelength,0,1,ICODE)
      CALL BUFFIN(PPUNIT1,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,     &
     &                  LEN_FIXHD)
          CMESSAGE='FFREAD : I/O error reading FIXED LENGTH HEADER'
          ICODE=2
          WRITE(6,*)'umthin1 - I/O error reading FIXED LENGTH HEADER'
          CALL ereport(RoutineName,icode,cmessage)
        ENDIF
      LEN_INTHD=PP_FIXHD(101)
      LEN_REALHD=PP_FIXHD(106)
      LEN1_LEVDPC=PP_FIXHD(111)
      LEN2_LEVDPC=PP_FIXHD(112)
      LEN1_ROWDPC=PP_FIXHD(116)
      LEN2_ROWDPC=PP_FIXHD(117)
      LEN1_COLDPC=PP_FIXHD(121)
      LEN2_COLDPC=PP_FIXHD(122)
      LEN1_LOOKUP=PP_FIXHD(151)
      LEN2_LOOKUP=PP_FIXHD(152)
! DEPENDS ON: dimens1
      CALL DIMENS1(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,        &
     &   LEN1_ROWDPC,LEN2_ROWDPC,LEN1_COLDPC,LEN2_COLDPC,               &
     &   LEN1_LOOKUP,LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,PPUNIT1,PPUNIT2,    &
     &   ICODE,CMESSAGE)
      IF(ICODE /= 0) THEN
        WRITE(7,100) ICODE
        WRITE(7,110) CMESSAGE
      ENDIF

      CALL ioShutdown()

      CALL ereport_finalise( )

 100  FORMAT(' ICODE EQUAL TO ',I2)
 110  FORMAT(A80)

      CALL ioShutdown()

      CALL ereport_finalise( )
      END PROGRAM FLDMOD
