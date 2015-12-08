! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface and arguments: ------------------------------------------
!  Routine: LOGICAL_TO_REAL ------------------------------------------
!
!  Purpose: To convert logical data within FIELD to real data.
!  the data in FIELD.
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  External documentation:
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE LOGICAL_TO_REAL_ffread1a(IDIM,LOGICAL_FIELD,FIELD,     &
     &                           NVALS,ILABEL,ICODE,CMESSAGE)
      USE lookup_addresses
      IMPLICIT NONE
      INTEGER                                                           &
     &     IDIM                                                         &
                                !IN  The full unpacked size of a field
     &    ,ILABEL(45)                                                   &
                                !OUT holds integer part of LOOKUP
     &    ,ICODE                !OUT Non zero for any error
      REAL                                                              &
     &     FIELD(IDIM)          !OUT On Input contains Real data.
      LOGICAL                                                           &
     &     LOGICAL_FIELD(IDIM)  !INOUT On Input contains logical data.
!                               ! contains the un-packed data.
      CHARACTER(LEN=*) cmessage    !OUT Will contain any error mesages.
!*
!     LOCAL  VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! Loop counter
     &    ,NVALS                  ! IN no of values in an input field
!
!
      DO  I=1,NVALS
        IF(LOGICAL_FIELD(I))THEN
          FIELD(I)=1.0
        ELSE
          FIELD(I)=0.0
        ENDIF
      ENDDO
      ILABEL(DATA_TYPE)=1     ! The data type must now be real
      ICODE=0
      RETURN
      END SUBROUTINE LOGICAL_TO_REAL_ffread1a
