! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Routine: INTEGER_TO_REAL 
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Small execs
      SUBROUTINE INTEGER_TO_REAL_fieldcos(IDIM,INTEGER_FIELD,FIELD,     &
     &                           NVALS,ILABEL,ICODE,CMESSAGE)
      USE lookup_addresses

      IMPLICIT NONE
!     arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)         !OUT error mesages.
      INTEGER                                                           &
     &     IDIM                                                         &
                                !IN full unpacked size of a field
     &    ,NVALS                                                        &
                                !IN no of values in an input field
     &    ,INTEGER_FIELD(IDIM)                                          &
                                !IN contains integer data.
     &    ,ILABEL(44)                                                   &
                                !OUT integer part of LOOKUP
     &    ,ICODE                !OUT error code
      REAL                                                              &
     &     FIELD(IDIM)          !OUT contains Real data.
!     Local variables
      INTEGER                                                           &
     &     I                    ! loop counter
!
!

      DO  I=1,NVALS
        FIELD(I)=INTEGER_FIELD(I)
      ENDDO
      ILABEL(DATA_TYPE)=1       ! The data type must now be real
      ICODE=0
      RETURN
      END SUBROUTINE INTEGER_TO_REAL_fieldcos
