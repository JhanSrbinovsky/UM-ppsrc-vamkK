! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   SUBROUTINE CMPS -----------------
!
!   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!
!   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 3,
!
!  Logical component number: S72
!
!   DOCUMENTATION:  Description of WGDOS in UMDP F3
!   -------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: STASH
      SUBROUTINE CMPS(IX,FIELD,ICOMP,NOP,APREC,IBIT,BASE,RMDI,          &
     &                ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

!     Subroutine arguments
      INTEGER IX,NOP,ICOMP(IX-1),IBIT
      INTEGER ICODE
      REAL FIELD(IX),APREC,BASE,RMDI
      CHARACTER CMESSAGE*80

!     Local variables
      LOGICAL OBTMIS,OBTZER
      INTEGER IMAP(IX),IMIS(IX),IZERO(IX),ITEMP(IX)
      INTEGER IGATH1(IX),IGATH2(IX)
      INTEGER I,IBASE,IMAX,I2,JBIT,JJ,JWORD,NMIS,NPOINT,NZERO
      INTEGER IDIV, jword_end, jbit_end, II, iword

      Real :: Atemp(Ix)
      Real :: Bprec

!     FUNCTIONS
      INTEGER ishift

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!         IX   =  LENGTH OF FIELD
!      FIELD   =  FIELD FOR COMPRESSING
!      ICOMP   =  RETURNED COMPRESSED DATA
!        NOP   =  NUMBER OF WORDS OF COMP FILLED BY THIS CALL
!      APREC   =  PRECISION
!       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
!       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
!       RMDI   =  MISSING DATA INDICATOR
!
!     INITIALISE VARIABLES/TEMP ARRAYS
!
      IF (lhook) CALL dr_hook('CMPS',zhook_in,zhook_handle)
      OBTMIS = .FALSE.
      OBTZER = .FALSE.
      BPREC  = 1./APREC
      DO I = 1,IX
          IMAP(I)  = 1
          IZERO(I) = 0
      END DO
      DO I = 1,IX-1
          ICOMP(I) = 0
      END DO
!
!     SCAN FIELD FOR MISSING DATA AND STORE RESULT IN IMIS,
!     SCAN FIELD FOR ZERO VALUES AND STORE RESULT IN IZERO,
!     SCAN FIELD FOR MINIMUM VALUE (IGNORE MISSING DATA INDICATOR)
!
      BASE  = HUGE(0.0)
      NMIS  = 0
      NZERO = 0
      JJ    = 0

      DO I = 1,IX
       IF (FIELD(I) /= RMDI) THEN
         IMIS(I)=0
       ELSE
         IMIS(I)=1
       ENDIF
      END DO

!
! GET NO. OF NON-RMDI POINTS + COMPRESS INDEX TO REMOVE THEM
!
       JJ  =0
       DO I=1,IX
         IF (FIELD(I) /= RMDI) THEN
           JJ        =JJ+1
           IGATH1(JJ)=I
           IF(FIELD(I) <  BASE) BASE=FIELD(I)
         END IF
! SET BASE VALUE
       END DO
!
      NMIS=IX-JJ
!
      IF(JJ /= 0)THEN
! REMOVE MISSING DATA
!DIR$ IVDEP
       DO I =1,JJ
! Pack data to required multple of precision.
       ATEMP(I)=NINT(FIELD(IGATH1(I))*BPREC)
       END DO
!
! GET NO. OF NON-ZERO (NON-RMDI) POINTS + COMPRESS INDEX TO REMOVE THEM
!
       NPOINT=0
       DO I  =1,JJ
         IF (ATEMP(I) /= 0.0) THEN
           IZERO(IGATH1(I))=0
           NPOINT        =NPOINT+1
           IGATH2(NPOINT)=I
         ELSE
           IZERO(IGATH1(I))=1
         END IF
       END DO
!
       NZERO=JJ-NPOINT
!
      ELSE
! If we have not checked base (if row contains only RMDI or ZERO so JJ == 0)
! then set BASE to -1.  IMAX is set later to -1 as well.
       BASE = -1.0
      ENDIF
 
!
!     CHECK IF A BITMAP FOR MISSING DATA IS NEEDED,
!
      IF (NMIS >  0) THEN
          OBTMIS = .TRUE.
      ELSE
          OBTMIS = .FALSE.
      END IF
!
!     ROUND BASE TO PRECISION REQUIRED
!
      IF (BASE /= 0.0) THEN
          BASE  = BASE*BPREC
          IBASE = NINT(BASE)
          BASE  = IBASE*APREC
      ELSE
          IBASE=0
      END IF
!
!     FIND DIFFERENCE FROM BASE AND SCALE
!     FIND MAXIMUM DIFFERENCE
!
      IMAX = -1
      DO I = 1,JJ
          ITEMP(I) = NINT(ATEMP(I)) - IBASE
          IF(ITEMP(I) <  0) ITEMP(I) = 0
          IF (IMAX <  ITEMP(I)) IMAX = ITEMP(I)
      END DO
!
!     FIND NUMBER OF BITS REQUIRED TO STORE MAX DIFFERENCE
!
      IBIT  = 0
      ! Enable a maximum range of 2^5 bits
      IF (IMAX  >   2.147483E+09) THEN
         ICODE = 2
         CMESSAGE='COEX: Unable to WGDOS pack to this accuracy'
         GOTO 9999
      ELSE
      IF (IMAX >  0) THEN
          I2    = 1
          DO WHILE(IMAX >= I2)
              IBIT  = IBIT + 1
              I2    = I2*2
          ENDDO
      ENDIF
      ENDIF
!
!     SET START POSITION IN OUTPUT ARRAY
!
      JWORD = 1
      JBIT  = 63
!
!     IF BIT-MAPPING FOR MISSING DATA THEN INSERT BITMAP
!
      IF (OBTMIS) THEN
          DO I = 1,IX
              IF (IMIS(I) == 1) THEN
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  IMAP(I)      = 0
              END IF
              IF (JBIT == 0) THEN
                  JBIT  = 63
                  JWORD = JWORD + 1
              ELSE
                  JBIT  = JBIT - 1
              END IF
          END DO
      END IF
!
!     IF WORTHWHILE USE BIT MAP AND COMPRESS OUT ZEROS.
!
      IF (IBIT >  0) THEN
          IF (NZERO >  IX/IBIT) THEN
              OBTZER = .TRUE.
              DO I = 1,IX
                  IF (IZERO(I) == 1) THEN
                      ICOMP(JWORD) = IBCLR(ICOMP(JWORD),JBIT)
                      IMAP(I)      = 0
                  ELSE
                      ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  END IF
                  IF (JBIT == 0) THEN
                      JBIT  = 63
                      JWORD = JWORD + 1
                  ELSE
                      JBIT  = JBIT - 1
                  END IF
              END DO
          ELSE
              OBTZER = .FALSE.
          END IF
!
!         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
!         AND SET POINTERS TO START OF NEXT 32 BIT BOUNDARY.
!
          IF (OBTZER .OR. OBTMIS) THEN
              DO I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
              END DO
              IF (JBIT /= 63) THEN
                  IF (JBIT >= 31) THEN
                      JBIT = 31
!
! We have set bits in the bottom half of the 64-bit word - clear
! them again
!
! DEPENDS ON: ishift
                      icomp(jword)=ishift(ishift(icomp(jword), -32), 32)
                  ELSE
                      JWORD = JWORD + 1
                      JBIT = 63
                  ENDIF
              ELSE
!
! We have set all the bits in the 64-bit word - clear them again
!
                icomp(jword)=0
                JBIT = 63
              ENDIF
          ELSE
!
! We have set all the bits in the 64-bit word - clear them again
!
              icomp(jword)=0
              JBIT = 63
          END IF
!
!         IF BIT MAPPING ZEROS - COMPRESS OUT UNWANTED ZERO VALUES
!        (OTHERWISE PACK ALL NON RMDI DATA (JJ) )
!
          IF (OBTZER) THEN
!DIR$ IVDEP
           DO I= 1,NPOINT
           ITEMP(I)=ITEMP(IGATH2(I))
           END DO
          ELSE
           NPOINT=JJ
          END IF
!
!         MOVE INTO OUTPUT ARRAY USING MINIMUM NUMBER OF BITS REQUIRED
!
          jword_end = jword + NPOINT*IBIT/64
          jbit_end  = jbit  - MOD(NPOINT*IBIT,64)
          IF (MOD(IBIT,32) == 0) THEN
              IDIV = 32
          ELSE IF (MOD(IBIT,16) == 0) THEN
              IDIV = 16
          ELSE IF (MOD(IBIT,8) == 0) THEN
              IDIV = 8
          ELSE IF (MOD(IBIT,4) == 0) THEN
              IDIV = 4
          ELSE IF (MOD(IBIT,2) == 0) THEN
              IDIV = 2
          ELSE
              IDIV = 1
          ENDIF
          DO II = 0, 63/IDIV
              iword = jword
!CDIR NODEP
              DO I = II+1,NPOINT,64/IDIV
!
! Adjust the position of the bits in the first word, and then
! insert the bits into the word
!
! DEPENDS ON: ishift
                  icomp(iword)=ior(ishift(itemp(i), (jbit+1)-ibit),     &
     &                             icomp(iword))
!
! Position the bits correctly for the second word, and 'or' them in
!
                  icomp(iword+1)=                                       &
! DEPENDS ON: ishift
     &             ior(ishift(itemp(i), min(64, 64-(ibit-(jbit+1)))),   &
     &                 icomp(iword+1))
                  iword = iword + IBIT/IDIV
              ENDDO
!
              JBIT = JBIT - IBIT
              IF (JBIT <  0) THEN
                  JWORD = JWORD + 1
                  JBIT  = JBIT + 64
              END IF
          ENDDO
          jword = jword_end
          jbit  = jbit_end
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          END IF

      ELSEIF(IBIT == 0) THEN
!
!         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
!
          IF (OBTMIS) THEN
              DO I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
              END DO
          END IF
      END IF
!
!     CALCULATE LENGTH IN 32 BIT WORDS
!
      NOP = JWORD*64 - JBIT - 1
      NOP = (NOP+31)/32
!
!     SET FLAGS TO INDICATE BIT MAPS
!
      IF (OBTZER) IBIT = IBIT + 128
      IF (OBTMIS) IBIT = IBIT + 32


 9999 CONTINUE
      IF (lhook) CALL dr_hook('CMPS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE CMPS
