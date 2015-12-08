! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE SCALE_FIELD(PDATA,RDATA,NPOINTS,SCALE_FACTOR,          &
     &                       LREC,PDATA_LEN,PACK_CODE,AMDI,             &
     &                      PLEN, ICODE, CMESSAGE)
!    subroutine to unpack a field, multiply by a scale factor,
!    then repack data.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      IMPLICIT NONE
      INTEGER NPOINTS,PDATA_LEN,PACK_CODE
      REAL FIELD(NPOINTS),RDATA(NPOINTS),SCALE_FACTOR,AMDI
      INTEGER PDATA(NPOINTS),NROW,NCOL,ISC,LWORD
      INTEGER LREC, PLEN
      INTEGER ICODE
      INTEGER I
      CHARACTER(LEN=80) CMESSAGE
      LOGICAL OPACK
      DATA LWORD/64/
      ICODE=0

!     Initialise FIELD variable
      DO I=1,NPOINTS
        FIELD(I) = 0.0
      ENDDO
      IF(PACK_CODE == 1) THEN
        OPACK=.FALSE.
! DEPENDS ON: coex
        CALL COEX(FIELD,NPOINTS,PDATA,NPOINTS,NROW,NCOL,PDATA_LEN,      &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

        DO I=1,NCOL*NROW
          IF(FIELD(I) /= AMDI) THEN
            FIELD(I) = FIELD(I) * SCALE_FACTOR
          ENDIF
        ENDDO

        OPACK=.TRUE.
! DEPENDS ON: coex
        CALL COEX(FIELD,NPOINTS,PDATA,NPOINTS,NROW,NCOL,PDATA_LEN,      &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

        PLEN = (PDATA_LEN + 1) /2

      ELSEIF(PACK_CODE == 4) THEN

! DEPENDS ON: runlen_decode
        CALL RUNLEN_DECODE(FIELD,NPOINTS,PDATA,LREC,                    &
     &                     AMDI,ICODE,CMESSAGE )
        DO I=1,NPOINTS
          IF(FIELD(I) /= AMDI) THEN
            FIELD(I) = FIELD(I) * SCALE_FACTOR
          ENDIF
        ENDDO

! DEPENDS ON: runlen_encode
        CALL RUNLEN_ENCODE(FIELD,NPOINTS,PDATA,NPOINTS,                 &
     &                     PLEN,AMDI,ICODE,CMESSAGE)
      ELSEIF(PACK_CODE == 0) THEN
        PLEN = LREC
        DO I=1,PDATA_LEN
          IF(RDATA(I) /= AMDI) THEN
            RDATA(I) = RDATA(I) * SCALE_FACTOR
          ENDIF
        ENDDO
      ELSE
        WRITE(6,*)pack_code,' not yet coded'
      ENDIF

      RETURN
      END SUBROUTINE SCALE_FIELD
