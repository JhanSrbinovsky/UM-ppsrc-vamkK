! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
! Purpose: To repack data from the input array FIELD and return

      SUBROUTINE RE_PACK(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,           &
     &                   ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
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
     &    ,RLABEL(19)           !    holds real part of LOOKUP
      CHARACTER(LEN=*) cmessage    !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
      REAL                                                              &
     &     WORK_ARRAY(IDIM)                                             &
                                  ! WORK array used for packing
     &    ,AMDI                   ! Missing data indicator.
      INTEGER                                                           &
     &     LEN_FULL_WORD                                                &
                                  ! The length of a FULL_WORD
     &    ,IXX                                                          &
                                  ! X dimension for COEX
     &    ,IYY                                                          &
                                  ! Y dimension for COEX
     &    ,ISC                                                          &
                                  ! Accuracy required for COEX
     &    ,IDUM                                                         &
                                  ! Dummy variable
     &    ,NUM_CRAY_WORDS                                               &
                                  ! IN no of values in an input field
     &    ,NUM_UNPACK_VALUES                                            &
                                  ! Number of numbers originally packed
     &    ,I                                                            &
                                  ! Loop counter
     &    ,PACK_CODE                                                    &
                                  ! Packing actually used
     &    ,GRIB_PACKING           ! OUT - profile for packing
!
!
      DATA LEN_FULL_WORD/64/
!
      AMDI=RLABEL(18)
! Packing actually used (might be different then packing wanted)
      PACK_CODE=PACK_TYPE

      IF(PACK_TYPE == 1) THEN     ! WGDOS packing
        IXX=ILABEL(LBNPT)
        IYY=ILABEL(LBROW)
        ISC=NINT(RLABEL(6))
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,WORK_ARRAY,IDIM,IXX,IYY,                   &
     &  NUM_CRAY_WORDS,ISC,.TRUE.,AMDI,LEN_FULL_WORD,                   &
     &  ICODE,CMESSAGE)
      ELSEIF(PACK_TYPE == 2) THEN !  32 Bit CRAY packing

      ELSEIF(PACK_TYPE == 4) THEN ! Run length encoding
        IXX=ILABEL(LBROW)
        IYY=ILABEL(LBNPT)
! DEPENDS ON: runlen_encode
        CALL RUNLEN_ENCODE(FIELD,IXX*IYY,WORK_ARRAY,IXX*IYY,            &
     &                     NUM_CRAY_WORDS,AMDI,ICODE,CMESSAGE)
        ! Size of run length encoded data is greater than unpacked
        ! field therefore leave field unpacked.
          if (NUM_CRAY_WORDS  >=  IXX*IYY) then
            PACK_CODE = 0
            DO I=1,IXX*IYY
             WORK_ARRAY(I) = FIELD(I)
            END DO
            NUM_CRAY_WORDS = IXX*IYY
          endif

      ELSE
        ICODE=6
        CMESSAGE=' UNPACK - packing type not yet supported'
      ENDIF
      DO I=1,NUM_cray_words
        FIELD(I)=WORK_ARRAY(I)
      END DO
      ILABEL(DATA_TYPE)=1  ! The data type must now be real
      ILABEL(LBPACK)=ILABEL(LBPACK)+PACK_CODE ! data now packed
      RETURN
      END SUBROUTINE RE_PACK
