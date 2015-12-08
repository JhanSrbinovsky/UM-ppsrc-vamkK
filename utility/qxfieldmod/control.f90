! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      SUBROUTINE CONTROL(PPUNIT1,PPUNIT2,LEN1_LOOKUP,LEN2_LOOKUP,       &
     &                   LOOKUP,PP_INTHD,LEN_INTHD,                     &
     &                   PP_FIXHD,LEN_FIXHD,ICODE,CMESSAGE,NENT)
     
      USE check_iostat_mod
     
      IMPLICIT NONE
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
     &    ,LEN_INTHD                                                    &
     &    ,LEN_LOOKUP                                                   &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,LOOKNEW(LEN1_LOOKUP,LEN2_LOOKUP)                             &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,PP_INTHD(LEN_INTHD)                                          &
     &    ,LEN_IO                                                       &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2                                                      &
     &    ,NENT
      INTEGER                                                           &
     &     ROW_LENGTH                                                   &
     &    ,P_ROWS                                                       &
     &    ,P_FIELD                                                      &
     &    ,LENBUF                                                       &
     &    ,I                                                            &
     &    ,J
      REAL                                                              &
     &     A_IO

      INTEGER                                                           &
     &    STIME_MOD                                                     &
     &,   ETIME_MOD                                                     &
     &,   NFIELDS_MOD                                                   &
     &,   MTYPE_MOD(500)                                                &
     &,   MLEVS_MOD(500)                                                &
     &,   STIME_SEL                                                     &
     &,   ETIME_SEL                                                     &
     &,   NFIELDS_SEL                                                   &
     &,   MTYPE_SEL(500)                                                &
     &,   MLEVS_SEL(500)                                                &
     &,   STIME_REJ                                                     &
     &,   ETIME_REJ                                                     &
     &,   NFIELDS_REJ                                                   &
     &,   MTYPE_REJ(500)                                                &
     &,   MLEVS_REJ(500)                                                &
     &,   PPUNIT_OROG                                                   &
     &,   STIME_THI                                                     &
     &,   ETIME_THI                                                     &
     &,   NFIELDS_THI                                                   &
     &,   MTYPE_THI(500)                                                &
     &,   MLEVS_THI(500)                                                &
     &,   IXXSTEP_THI(500)                                              &
     &,   IYYSTEP_THI(500)

      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error

      CHARACTER                                                         &
     &    OUTPUT_PACK_TYPE*6

      REAL                                                              &
     &    AMULT(500)                                                    &
     &,   WIND_10M_OROG                                                 &
                                 !  LEVEL ABOVE WHICH 10M WIND FIXED
     &,   WIND_10M_SCALE         !  SCALE APPLIED TO LEVEL 1 WINDS

      LOGICAL                                                           &
     &    MODIFY                                                        &
     &,   REJECT                                                        &
     &,   SELECT                                                        &
     &,   WIND_10M                                                      &
     &,   THIN

      NAMELIST /MODS/                                                   &
     &  MODIFY,STIME_MOD,ETIME_MOD,NFIELDS_MOD,                         &
     &                                      MTYPE_MOD,MLEVS_MOD,AMULT,  &
     &  SELECT,STIME_SEL,ETIME_SEL,NFIELDS_SEL,MTYPE_SEL,MLEVS_SEL,     &
     &  REJECT,STIME_REJ,ETIME_REJ,NFIELDS_REJ,MTYPE_REJ,MLEVS_REJ,     &
     &  WIND_10M,WIND_10M_SCALE,WIND_10M_OROG,PPUNIT_OROG,              &
     &  THIN,STIME_THI,ETIME_THI,NFIELDS_THI,MTYPE_THI,MLEVS_THI,       &
     &                                        IXXSTEP_THI,IYYSTEP_THI,  &
     &  OUTPUT_PACK_TYPE

!-----------------------------------------------------------------------
      CHARACTER(LEN=*) cmessage
      EXTERNAL FIELDS
!
!L---------------------------------------------------------------
!L     init namelist
!L---------------------------------------------------------------
      MODIFY   = .FALSE.
      REJECT   = .FALSE.
      SELECT   = .FALSE.
      WIND_10M = .FALSE.
      THIN=.FALSE.
      STIME_MOD = -99
      ETIME_MOD = -99
      NFIELDS_MOD=0
      STIME_SEL = -99
      ETIME_SEL = -99
      NFIELDS_SEL=0
      STIME_REJ = -99
      ETIME_REJ = -99
      NFIELDS_REJ=0
      STIME_THI = -99
      ETIME_THI = -99
      NFIELDS_THI=0
      DO I=1,500
        MTYPE_MOD(I)=0
        MLEVS_MOD(I)=0
        AMULT(I)=1.0
        MTYPE_SEL(I)=0
        MLEVS_SEL(I)=0
        MTYPE_REJ(I)=0
        MLEVS_REJ(I)=0
        MTYPE_THI(I)=0
        MLEVS_THI(I)=0
        IXXSTEP_THI(I)=2
        IYYSTEP_THI(I)=2
      ENDDO
      WIND_10M_OROG  = -9999.
      WIND_10M_SCALE = .7
      PPUNIT_OROG    = 12
      OUTPUT_PACK_TYPE='WGDOS '
!
!L---------------------------------------------------------------
!L     read namelist
!L---------------------------------------------------------------
      READ (UNIT=5, NML=MODS, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist MODS")
      WRITE(7,MODS)

!L---------------------------------------------------------------
!L     Set up constants
!L---------------------------------------------------------------
      ROW_LENGTH=PP_INTHD(6)
      P_ROWS=PP_INTHD(7)
      P_FIELD=ROW_LENGTH*P_ROWS
      LENBUF=P_FIELD + 512

! DEPENDS ON: fields
      CALL FIELDS(PP_FIXHD,LEN_FIXHD,LENBUF,P_FIELD,                    &
     &             LOOKUP,LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,NENT,          &
     &             STIME_MOD,ETIME_MOD,NFIELDS_MOD,                     &
     &                                       MTYPE_MOD,MLEVS_MOD,AMULT, &
     &             STIME_SEL,ETIME_SEL,NFIELDS_SEL,MTYPE_SEL,MLEVS_SEL, &
     &             STIME_REJ,ETIME_REJ,NFIELDS_REJ,MTYPE_REJ,MLEVS_REJ, &
     &             STIME_THI,ETIME_THI,NFIELDS_THI,MTYPE_THI,MLEVS_THI, &
     &                                         IXXSTEP_THI,IYYSTEP_THI, &
     &             MODIFY,SELECT,REJECT,THIN,OUTPUT_PACK_TYPE,          &
     &             WIND_10M,WIND_10M_OROG,WIND_10M_SCALE,PPUNIT_OROG,   &
     &             PPUNIT1,PPUNIT2,ICODE,CMESSAGE)
 9999 CONTINUE
      RETURN
      END SUBROUTINE CONTROL

