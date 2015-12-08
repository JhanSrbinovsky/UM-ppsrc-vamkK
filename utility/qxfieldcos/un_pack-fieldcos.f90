! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Suboutine: UN_PACK  
!LL
!LL  Purpose: To unpack data from the input array FIELD and return
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
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Small execs
      SUBROUTINE UN_PACK_fieldcos(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,  &
     &          ILABEL,LEN_ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)
        USE ereport_mod, ONLY : ereport
        USE lookup_addresses

        IMPLICIT NONE
!     arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)         !OUT error mesages.
      INTEGER                                                           &
     &     PACK_TYPE                                                    &
                                !INOUT Type of packing used
     &    ,IDIM                                                         &
                                !IN    full unpacked size of a field
     &    ,PP_FIXHD(*)                                                  &
                                !IN    PPfile fixed length header
     &    ,NUM_CRAY_WORDS                                               &
                                !IN    length of input field
     &    ,LEN_ILABEL                                                   &
                                !IN    length of ilabel array
     &    ,ILABEL(LEN_ILABEL)                                           &
                                !INOUT holds integer part of LOOKUP
     &    ,ICODE                !OUT   Non zero for any error
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                !INOUT Input contains packed data.
!                               !      Output contains un-packed data.
     &    ,AMDI                 !IN    Missing data indicator.
!*---------------------------------------------------------------------
!     arguments for called routines
      INTEGER                                                           &
     &     LEN_FULL_WORD                                                &
                                  ! The length of a FULL_WORD
     &    ,IXX                                                          &
                                  ! Returned X dimension from COEX
     &    ,IYY                                                          &
                                  ! Returned Y dimension from COEX
     &    ,IDUM                   ! Dummy variable
      REAL                                                              &
     &     WORK_ARRAY(IDIM)       !WORK array used for un_packing

!     LOCAL  VARIABLES
      INTEGER                                                           &
     &     NUM_UNPACK_VALUES                                            &
                                  ! Number of numbers originally packed
     &    ,I                      ! loop counter
!
!
      DATA LEN_FULL_WORD/64/
!
      IF(PACK_TYPE == 1) THEN     ! WGDOS packing
! DEPENDS ON: coex
        CALL COEX(WORK_ARRAY,IDIM,FIELD,NUM_CRAY_WORDS,IXX,IYY,         &
     &  IDUM,IDUM,.FALSE.,AMDI,LEN_FULL_WORD,ICODE,CMESSAGE)
        NUM_UNPACK_VALUES=IXX*IYY
        ILABEL(LBLREC)=ILABEL(LBROW)*ILABEL(LBNPT)+ILABEL(LBEXT)
      ELSEIF(PACK_TYPE == 2) THEN !  32 Bit CRAY packing
        NUM_CRAY_WORDS=NUM_CRAY_WORDS*2
! DEPENDS ON: expand21
        CALL EXPAND21(NUM_CRAY_WORDS,FIELD,WORK_ARRAY)
        NUM_UNPACK_VALUES=NUM_CRAY_WORDS
      ELSEIF(PACK_TYPE == 3) THEN !  GRIB packing





        WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'

        CALL EREPORT('UN_PACK', 1000,                                   &
     &   'GRIB unpacking only support on NEC SX6/8')


      ELSEIF(PACK_TYPE == 4) THEN !  Run Length encoding
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
      ILABEL(DATA_TYPE)=1  ! data must now be real
      ILABEL(LBPACK)=ILABEL(LBPACK)-PACK_TYPE ! data no longer packed
      PACK_TYPE=0          ! data now not packed
      RETURN
      END SUBROUTINE un_pack_fieldcos
