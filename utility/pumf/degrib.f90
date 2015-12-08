! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

      SUBROUTINE DEGRIB(FIELD,WORK_ARRAY,IDIM,NUM_CRAY_WORDS,           &
     &                   ILABEL,AMDI,NUM_UNPACK_VALUES,LEN_FULL_WORD)
!  
!    Routine: DEGRIB----------------------------------------------
!  
!    Purpose: Routine to unpack GRIB record from field and return in
!             WORK_ARRAY
!  
!    Reviewer:  Date of review:
!  
!    Tested under compiler: cft77
!    Tested under OS version: UNICOS 7
!  
!   Programming standard: UM Doc Paper 3, version
!  
!   Logucal component number:
!  
!   Project task:
!  
!  
!   Documentation:
!     UMDP
!  
!   -------------------------------------------------------------
      IMPLICIT NONE
      INTEGER                                                           &
     &     IDIM                                                         &
     &    ,NUM_CRAY_WORDS                                               &
     &    ,ILABEL(64)                                                   &
     &    ,NUM_UNPACK_VALUES                                            &
     &    ,LEN_FULL_WORD
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
     &    ,WORK_ARRAY(IDIM)                                             &
     &    ,AMDI

      INTEGER                                                           &
     &     LEN_VERT                                                     &
     &    ,NUM_VERT                                                     &
     &    ,LEN_BITMAP                                                   &
     &    ,NUM_BITMAP                                                   &
     &    ,LEN_Q                                                        &
     &    ,NUM_Q                                                        &
     &    ,WIDTH                                                        &
     &    ,LEN_B0                                                       &
     &    ,LEN_B1                                                       &
     &    ,LEN_B2                                                       &
     &    ,LEN_B3                                                       &
     &    ,LEN_B4                                                       &
     &    ,LEN_BR                                                       &
     &    ,LEN_WRKI                                                     &
     &    ,LEN_WRKI2                                                    &
     &    ,LEN_WRKR
      PARAMETER (                                                       &
     &     LEN_VERT=4                                                   &
     &    ,LEN_Q=4                                                      &
     &    ,LEN_B0=4                                                     &
     &    ,LEN_B1=30                                                    &
     &    ,LEN_B2=20                                                    &
     &    ,LEN_B3=2                                                     &
     &    ,LEN_B4=2                                                     &
     &    ,LEN_BR=20                                                    &
     &    ,LEN_WRKI=500                                                 &
     &    ,LEN_WRKI2=1000                                               &
     &    ,LEN_WRKR=500)
      INTEGER                                                           &
     &     QUASI(LEN_Q)                                                 &
     &    ,BITMAP(IDIM)                                                 &
     &    ,BLOCK0(LEN_B0)                                               &
     &    ,BLOCK1(LEN_B1)                                               &
     &    ,BLOCK2(LEN_B2)                                               &
     &    ,BLOCK3(LEN_B3)                                               &
     &    ,BLOCK4(LEN_B4)                                               &
     &    ,POSN(4)                                                      &
     &    ,WORD                                                         &
     &    ,OFF                                                          &
     &    ,ERROR                                                        &
     &    ,WORKINT(LEN_WRKI)                                            &
     &    ,WORKINT2(LEN_WRKI2)
      REAL                                                              &
     &     VERT_COORDS(LEN_VERT)                                        &
     &    ,BLOCKR(LEN_BR)                                               &
     &    ,WORKRE(LEN_WRKR)

      LEN_BITMAP=IDIM
      ERROR=6
      POSN(1)=0
      OFF=0

! DEPENDS ON: decode

      CALL DECODE(WORK_ARRAY,WORK_ARRAY2,IDIM,NUM_UNPACK_VALUES,        &
     &            VERT_COORDS,LEN_VERT,NUM_VERT,                        &
     &            BITMAP,LEN_BITMAP,NUM_BITMAP,                         &
     &            QUASI,LEN_Q,NUM_Q,                                    &
     &            WIDTH,LEN_FULL_WORD,                                  &
     &            BLOCK0,BLOCK1,BLOCK2,BLOCK3,BLOCK4,BLOCKR,            &
     &            FIELD,IDIM,POSN,WORD,OFF,ERROR,                       &
     &            WORKINT,WORKINT2,WORKRE)

      RETURN
      END SUBROUTINE DEGRIB
