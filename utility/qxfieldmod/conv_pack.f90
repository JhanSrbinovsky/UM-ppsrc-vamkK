! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs

      SUBROUTINE CONV_PACK(ILABEL,RLABEL,PACK_CODE,                     &
     &                     INPUT_PACK_TYPE,OUTPUT_PACK_TYPE,            &
     &                     FIELD,IDIM,LEN_FIELD,                        &
     &                     PP_FIXHD,ICODE,CMESSAGE)
      USE lookup_addresses
      IMPLICIT NONE

      INTEGER                                                           &
     &     ILABEL(50)                                                   &
     &    ,PACK_CODE                                                    &
     &    ,IDIM                                                         &
     &    ,PP_FIXHD(*)                                                  &
     &    ,LEN_FIELD                                                    &
     &    ,ICODE
      REAL                                                              &
     &     RLABEL(19)                                                   &
     &    ,FIELD(IDIM)
      CHARACTER                                                         &
     &     INPUT_PACK_TYPE*6                                            &
     &    ,OUTPUT_PACK_TYPE*6                                           &
     &    ,CMESSAGE*(*)
      REAL                                                              &
     &     AMDI

      AMDI=RLABEL(18)

! DEPENDS ON: un_pack_ffread1a
      CALL UN_PACK_ffread1a(PACK_CODE,IDIM,FIELD,LEN_FIELD,             &
     &             ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)

      IF(ICODE /= 0) THEN
        WRITE(7,*) ICODE
        WRITE(7,*) CMESSAGE
      ENDIF

      LEN_FIELD = ILABEL(LBROW) * ILABEL(LBNPT)
      WRITE(6,*) INPUT_PACK_TYPE,' NOW UNPACKED'

      IF(OUTPUT_PACK_TYPE == 'NONE  ') THEN
        pack_code=0   ! no repacking needed
      ELSEIF(OUTPUT_PACK_TYPE == 'WGDOS ') THEN
        pack_code=1   ! repack using coex
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ELSEIF(OUTPUT_PACK_TYPE == 'CRAY32') THEN
        pack_code=3   ! repack using cray 32
        WRITE(6,*) 'packing not supported'

      ELSEIF(OUTPUT_PACK_TYPE == 'GRIB  ') THEN
        pack_code=3   ! repack using grib
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ELSEIF(OUTPUT_PACK_TYPE == 'RUNLEN') THEN
        pack_code=4   ! repack using run length encoding
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ENDIF

      IF(ICODE /= 0) THEN
        WRITE(7,*) ICODE
        WRITE(7,*) CMESSAGE
      ENDIF
      ILABEL(LBLREC) = LEN_FIELD

      WRITE(6,*) 'NOW PACKED INTO ',OUTPUT_PACK_TYPE
!     WRITE(6,*)'len_field=',len_field

      RETURN
      END SUBROUTINE CONV_PACK
