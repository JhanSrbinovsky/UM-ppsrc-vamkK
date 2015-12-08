! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:



!LL  Routine: CHECK_EXTRA ----------------------------------------------
!LL
!LL  Purpose: To check that code is correct for vector
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
!*L  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE CHECK_EXTRA(CODE,DATA_VALUES,ICODE,CMESSAGE)
      IMPLICIT NONE
! Description: Define valid extra data vector types for
!              use in pp fields.
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.4   11/04/02   Original code.  R. Hill
!===============================================================
      INTEGER,PARAMETER :: NO_EXTRA_VECTORS=14 ! Number of valid extra
                                               ! data vector types
      INTEGER :: EXTRA_VECTOR(NO_EXTRA_VECTORS)

      ! Vector type 10 is not currently supported!
      DATA EXTRA_VECTOR/1,2,3,4,5,6,7,8,9,11,12,13,14,15/
      LOGICAL VALID_TYPE ! Flag to indicate valid vector type
      INTEGER I ! LOOP counter
!     arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)         !OUT error message
      INTEGER                                                           &
     &     CODE                                                         &
                                !IN Code to be checked
     &    ,DATA_VALUES                                                  &
                                !IN Number of data values in vector
     &    ,ICODE                !OUT error code
!     Local variables
      INTEGER                                                           &
     &     TYPE

      DATA_VALUES=CODE/1000
      TYPE=CODE-DATA_VALUES*1000
      VALID_TYPE = .FALSE.
      DO I = 1, NO_EXTRA_VECTORS
         IF (TYPE == EXTRA_VECTOR(I)) VALID_TYPE = .TRUE.
      ENDDO

      IF (.NOT.VALID_TYPE) THEN
         ICODE = 1
         CMESSAGE='CHECK_DATA: Unsupported extra data vector code'
         WRITE(6,*) "Unsupported Extra Data code", TYPE
      ENDIF

      RETURN
      END SUBROUTINE CHECK_EXTRA



