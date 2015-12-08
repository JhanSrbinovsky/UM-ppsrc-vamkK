! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: DIMENS1--------------------------------------------
!
! Purpose: To read a   direct access PP file  and convert it to a
! sequential file read to be passed across to the IBM
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      SUBROUTINE DIMENS1(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,  &
     &   LEN1_ROWDPC,LEN2_ROWDPC, LEN1_COLDPC,LEN2_COLDPC,              &
     &   LEN1_LOOKUP,LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,PPUNIT1,PPUNIT2,    &
     &   ICODE,CMESSAGE)
      IMPLICIT NONE
      EXTERNAL READPP
      CHARACTER(LEN=*) cmessage
      INTEGER                                                           &
     &     LEN_INTHD                                                    &
     &    ,LEN_FIXHD                                                    &
     &    ,LEN_REALHD                                                   &
     &    ,LEN1_LEVDPC                                                  &
     &    ,LEN2_LEVDPC                                                  &
     &    ,LEN1_ROWDPC                                                  &
     &    ,LEN2_ROWDPC                                                  &
     &    ,LEN1_COLDPC                                                  &
     &    ,LEN2_COLDPC                                                  &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!
! DEPENDS ON: readpp
      CALL READPP(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,         &
     &     LEN1_ROWDPC,LEN2_ROWDPC,LEN1_COLDPC,LEN2_COLDPC,LEN1_LOOKUP, &
     &     LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,LOOKUP,LOOKUP,PPUNIT1,        &
     &     PPUNIT2,ICODE,CMESSAGE)
 9999 CONTINUE
      IF(ICODE /= 0) RETURN
      RETURN
      END SUBROUTINE DIMENS1
