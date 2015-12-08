! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE THIN_FIELD(PDATA,RDATA,PDATA_LEN,IXX,IYY,              &
     &                      IXXSTEP,IYYSTEP,IDIM,PACK_CODE,AMDI,        &
     &                      LREC, ICODE, CMESSAGE)
!
!    Subroutine to unpack a field, thin, then repack data.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
!
      IMPLICIT NONE
      INTEGER IDIM,PDATA_LEN,PACK_CODE,IXXSTEP,IYYSTEP
      INTEGER LREC,ICODE
      INTEGER PDATA(IDIM),IXX,IYY,ISC,LWORD
      INTEGER i,j,k,kk,ix1,iy1
      integer countx,county
      REAL RDATA(IDIM),FIELD(IDIM),AMDI
      LOGICAL OPACK
      CHARACTER(LEN=80) CMESSAGE
      DATA LWORD/64/

      ICODE=0
      IF(PACK_CODE == 1) THEN
        OPACK=.FALSE.
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,PDATA,IDIM,IXX,IYY,PDATA_LEN,              &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

! If IXX and IYY are not decreased by 1 then GRDSET ( a PP routine)
! will fail and give the message 'BAD GRID DEFINITION'.
! Unfortunately the same failure occurs if IXX and IYY are decreased
! when a step size of 1 is specified so IXX and IYY will only be
! decreased for step sizes > 1 (in case anyone uses a step size of 1
! instead of using SELECT in the namelist)

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            FIELD(K) = FIELD(KK)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP

        OPACK=.TRUE.
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,PDATA,IDIM,IXX,IYY,PDATA_LEN,              &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)
        IF(ICODE /= 0) THEN
          WRITE(7,*) ICODE
          WRITE(7,*) 'THIN_FIELD FAILED WITH THE FOLLOWING REASON:'
          WRITE(7,*) CMESSAGE
        ENDIF

      ELSE IF(PACK_CODE == 4) THEN

! DEPENDS ON: runlen_decode
        CALL RUNLEN_DECODE(FIELD,IXX*IYY,PDATA,LREC,                    &
     &                     AMDI,ICODE,CMESSAGE )

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            FIELD(K) = FIELD(KK)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP

! DEPENDS ON: runlen_encode
        CALL RUNLEN_ENCODE(FIELD,IXX*IYY,PDATA,IXX*IYY,                 &
     &                     PDATA_LEN,AMDI,ICODE,CMESSAGE)
      ELSE IF(PACK_CODE == 0) THEN

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            RDATA(K) = RDATA(kk)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP
        PDATA_LEN = IXX * IYY

      ELSE
        WRITE(6,*)pack_code,' not yet coded'
      END IF

      RETURN
      END SUBROUTINE THIN_FIELD
