! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE WIND_10M_FIX(PDATA,RDATA,PDATA_LEN,                    &
     &                        FCT,ITYPE,LEVEL,IPROJ,PPUNIT1,            &
     &                        WIND_10M_SCALE,WIND_10M_OROG,             &
     &                        MODEL_OROG,ILABEL_OROG,RLABEL_OROG,       &
     &                        IDIM,PACK_CODE,AMDI)
!
!    subroutine to unpack a 10m winds and replace if posible by
!    the level 1 wind scaled using wind_10m_scale
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
      IMPLICIT NONE
      INTEGER IDIM,PDATA_LEN,PACK_CODE
      INTEGER ICODE
      CHARACTER CMESSAGE*80
      REAL RDATA(IDIM),FIELD(IDIM),FIELD1(IDIM),AMDI
      REAL MODEL_OROG(IDIM),RLABEL_OROG(19),RLABEL(19)
      REAL WIND_10M_OROG,WIND_10M_SCALE
      INTEGER PDATA(IDIM),NROW,NCOL,ISC,LWORD
      INTEGER ILABEL_OROG(45),ILABEL(45)
      INTEGER FCT,ITYPE,ITYPE1,LEVEL,LEVEL1,IPROJ,PPUNIT1
      INTEGER IEXTRA(10)
      INTEGER I
      INTEGER IXX,IYY
      LOGICAL OPACK
      DATA LWORD/64/

      DO I=1,10
        IEXTRA(I)=0
      ENDDO

      write(6,*) ' read level1 winds'
      ITYPE1 = 6
      IF(ITYPE == 75) ITYPE1 = 5
      LEVEL1 = 1
! DEPENDS ON: ffread
      CALL FFREAD(IPROJ,FCT,ITYPE1,LEVEL1,PPUNIT1,FIELD1,IDIM,          &
     &                ILABEL,RLABEL,IEXTRA,ICODE,CMESSAGE)

      write(6,*) 'icode=',icode
      IF(ICODE == 0) THEN
        IF(PACK_CODE == 1) THEN
          OPACK=.FALSE.
          WRITE(6,*)'call coex'
! DEPENDS ON: coex
          CALL COEX(FIELD,IDIM,PDATA,IDIM,NROW,NCOL,PDATA_LEN,          &
     &              ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

          WRITE(6,*)'loop field'
          DO I=1,NCOL*NROW
            IF(FIELD(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN
           WRITE(6,*)i,model_orog(i),field(i),field1(i),field1(i)*.8
                FIELD(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO

          OPACK=.TRUE.
! DEPENDS ON: coex
          CALL COEX(FIELD,IDIM,PDATA,IDIM,NROW,NCOL,PDATA_LEN,          &
     &              ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)
        ELSEIF(PACK_CODE == 4) THEN
          OPACK=.FALSE.

! DEPENDS ON: runlen_decode
          CALL RUNLEN_DECODE(FIELD,IDIM,PDATA,IDIM,                     &
     &                       AMDI,ICODE,CMESSAGE )

          DO I=1,PDATA_LEN
            IF(FIELD(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN

                FIELD(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO
          IXX=ILABEL(18)
          IYY=ILABEL(19)
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(FIELD,IXX*IYY,PDATA,IXX*IYY,               &
     &                       PDATA_LEN,AMDI,ICODE,CMESSAGE)
        ELSEIF(PACK_CODE == 0) THEN
          DO I=1,PDATA_LEN
            IF(RDATA(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN
                RDATA(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO
        ELSE
          WRITE(6,*)pack_code,' not yet coded'
        ENDIF
      ELSE
        WRITE(7,*) ICODE
        WRITE(7,*) '10M WIND FAILED WITH THE FOLLOWING REASON:'
        WRITE(7,*) CMESSAGE
      ENDIF

      RETURN
      END SUBROUTINE WIND_10M_FIX
