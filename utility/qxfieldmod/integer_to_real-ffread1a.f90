! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: INTEGER_TO_REAL ------------------------------------------
!
! Purpose: To convert logical data within FIELD to real data.
! the data in FIELD.
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! Logical components covered: ...
!
! Project task: ...
!
! External documentation:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE INTEGER_TO_REAL_ffread1a(IDIM,INTEGER_FIELD,FIELD,     &
     &                           NVALS,ILABEL,ICODE,CMESSAGE)
      USE lookup_addresses
      IMPLICIT NONE
      INTEGER                                                           &
     &     IDIM                                                         &
                                !IN  The full unpacked size of a field
     &    ,ILABEL(45)                                                   &
                                !OUT holds integer part of LOOKUP
     &    ,ICODE                                                        &
                                !OUT Non zero for any error
     &    ,INTEGER_FIELD(IDIM)                                          &
                                !IN  On input contains integer data.
     &    ,NVALS                !IN no of values in an input field
      
      REAL                                                              &
     &     FIELD(IDIM)          !OUT On Input contains Real data.
!                               ! contains the un-packed data.
      CHARACTER(LEN=*) cmessage    !OUT Will contain any error mesages.
!*
!     LOCAL  VARIABLES
      INTEGER                                                           &
     &     I                    ! Loop counter
!
!
      DO  I=1,NVALS
        FIELD(I)=INTEGER_FIELD(I)
      ENDDO
      ILABEL(DATA_TYPE)=1     ! The data type must now be real
      ICODE=0
      RETURN
      END SUBROUTINE INTEGER_TO_REAL_ffread1a
