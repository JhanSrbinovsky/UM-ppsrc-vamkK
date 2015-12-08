! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface and arguments: ------------------------------------------
!  Routine: UN_PACK  -------------------------------------------------
!
!  Purpose: To unpack data from the input array FIELD and return
!  the data in FIELD.
!
!  Tested under compiler:   cft77
!  Tested under OS version: UNICOS 5.1
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
      SUBROUTINE UN_PACK_ffread1a(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,  &
     &                   ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)
        USE ereport_mod, ONLY : ereport
      USE lookup_addresses
      IMPLICIT NONE
      INTEGER                                                           &
     &     PACK_TYPE                                                    &
                                !IN  The type of packing used
     &    ,IDIM                                                         &
                                !IN  The full unpacked size of a field
     &    ,ILABEL(45)                                                   &
                                !OUT holds integer part of LOOKUP
     &    ,ICODE                                                        &
                                !OUT Non zero for any error
     &    ,PP_FIXHD(*)          !IN  PPfile fixed length header
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                !INOUT On Input contains data.On output
     &    ,AMDI                 !IN  Missing data indicator.
!                               ! contains the un-packed data.
      CHARACTER(LEN=*) cmessage    !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
      REAL                                                              &
     &     WORK_ARRAY(IDIM)       !WORK array used for un_packing
      INTEGER                                                           &
     &     LEN_FULL_WORD                                                &
                                  ! The length of a FULL_WORD
     &    ,IXX                                                          &
                                  ! Returned X dimension from COEX
     &    ,IYY                                                          &
                                  ! Returned Y dimension from COEX
     &    ,IDUM                                                         &
                                  ! Dummy variable
     &    ,NUM_CRAY_WORDS                                               &
                                  ! IN no of values in an input field
     &    ,I                                                            &
                                  ! Loop counter
     &    ,NUM_UNPACK_VALUES      ! Number of numbers originally packed
!
!
      DATA LEN_FULL_WORD/64/
!
      IF(PACK_TYPE == 1) THEN     ! WGDOS packing
! DEPENDS ON: coex
        CALL COEX(WORK_ARRAY,IDIM,FIELD,NUM_CRAY_WORDS,IXX,IYY,         &
     &  IDUM,IDUM,.FALSE.,AMDI,LEN_FULL_WORD,ICODE,CMESSAGE)
        NUM_UNPACK_VALUES=IXX*IYY
      ELSEIF(PACK_TYPE == 3) THEN !  GRIB packing





        WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'

        CALL EREPORT('UN_PACK', 1000,                                   &
     &   'GRIB unpacking only support on NEC SX6/8')


      ELSEIF(PACK_TYPE == 4) THEN ! Run length encoded data
        NUM_UNPACK_VALUES = ILABEL(LBNPT) * ILABEL(LBROW)
! DEPENDS ON: runlen_decode
        CALL RUNLEN_DECODE(WORK_ARRAY,IDIM,FIELD,NUM_CRAY_WORDS,        &
     &                     AMDI,ICODE,CMESSAGE)
      ELSE
        ICODE=6
        CMESSAGE=' UNPACK - packing type not yet supported'
      ENDIF
      DO I=1,NUM_UNPACK_VALUES
      FIELD(I)=WORK_ARRAY(I)
      END DO
      ILABEL(DATA_TYPE)=1  ! The data type must now be real
      ILABEL(LBPACK)=ILABEL(LBPACK)-PACK_TYPE ! data no longer packed

      IF(ICODE  /=  0) CALL EREPORT("UN_PACK", ICODE, CMESSAGE)
      RETURN
      END SUBROUTINE UN_PACK_ffread1a
