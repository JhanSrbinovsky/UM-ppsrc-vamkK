! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Derived datatypes used by IOS

MODULE IOS_types
  USE UM_types
  USE IOS_Constants
  IMPLICIT NONE
!-------------------------------------------------------
! IO Server metadata type definition
!-------------------------------------------------------

    TYPE IOS_metadata_type
       SEQUENCE
   ! Sequence and order of definition are important in comms.
   ! Action
       INTEGER                       :: action
   ! Common to all
       INTEGER                       :: unit
   ! Open/Close
       INTEGER                       :: name_length
       INTEGER                       :: intent
       INTEGER                       :: delete
   ! Read/Write
       INTEGER                       :: data_size
   ! Setpos
       INTEGER                       :: address
       INTEGER                       :: client
       INTEGER                       :: handle
       INTEGER                       :: Originating_Slot

   ! Multi-purpose String
       CHARACTER (LEN=IOS_string_max):: string
    END TYPE IOS_metadata_type

    TYPE IOS_async_object
      SEQUENCE
      INTEGER                           :: seq_id
      INTEGER                           :: request(3)
      INTEGER                           :: state
      INTEGER                           :: pe
      TYPE(IOS_metadata_type)           :: md
      INTEGER(KIND=IOS_TUKind), POINTER :: payload(:)
    END TYPE IOS_async_object

    TYPE IOS_STATUS
      SEQUENCE
      INTEGER                       :: unit
      INTEGER                       :: Owner   ! The rank that opened the file
      INTEGER                       :: position! The byte position in the file
      INTEGER                       :: extent  ! The byte length of the file
      LOGICAL                       :: isOpen
      LOGICAL                       :: isReadOnly
      CHARACTER (LEN=IOS_string_max):: filename
    END TYPE IOS_STATUS    
    INTEGER, PARAMETER              :: ios_state_size=6

  END MODULE IOS_types
