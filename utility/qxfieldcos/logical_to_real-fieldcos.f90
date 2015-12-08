! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Subroutine: LOGICAL_TO_REAL ------------------------------------------
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs
      SUBROUTINE LOGICAL_TO_REAL_fieldcos(IDIM,LOGICAL_FIELD,FIELD,     &
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
     &    ,ILABEL(44)                                                   &
                                !OUT integer part of LOOKUP
     &    ,ICODE                !OUT error code
      REAL                                                              &
     &     FIELD(IDIM)          !OUT contains Real data.
      LOGICAL                                                           &
     &     LOGICAL_FIELD(IDIM)  !IN contains logical data.
!                               ! contains the un-packed data.
!     Local variables
      INTEGER                                                           &
     &     I                    ! loop counter
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
      END SUBROUTINE LOGICAL_TO_REAL_fieldcos
